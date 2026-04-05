#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex150_lyapunov_pairwise_vs_coulomb3b.py
=======================================
P1 (cd.): ``lambda_max`` (Benettin, wektor styczny) dla tego samego IC
(Burrau), porownanie

  - ``pairwise`` (tylko V2),
  - ``coulomb_3b`` (V2 + dokladne I Coulomb w V3, ``eom_tgp``).

Jacobian ``acc`` wzgledem ``x`` — różnice centralne (kosztowne).

``--quick``: krotszy czas (~30–90 s zaleznie od CPU).

Uruchomienie: ``python ex150_lyapunov_pairwise_vs_coulomb3b.py [--quick]``
"""

from __future__ import annotations

import argparse
import os
import sys
import time

import numpy as np

_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

from nbody.dynamics_backends import build_tgp_integration_pair
from nbody.lyapunov import largest_lyapunov_exponent_benettin, pythagorean_three_body_burrau


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true")
    args = parser.parse_args()

    beta = gamma = 0.08
    soft = 1e-6
    pos0, vel0, C = pythagorean_three_body_burrau()

    acc_pair, _ = build_tgp_integration_pair(
        "pairwise",
        beta=beta,
        gamma=gamma,
        softening=soft,
        include_3body=True,
    )
    acc_3b, _ = build_tgp_integration_pair(
        "coulomb_3b",
        beta=beta,
        gamma=gamma,
        softening=soft,
        include_3body=True,
    )

    if args.quick:
        t_final = 5.0
        dt = 0.012
        renorm = 15
        jac_eps = 1e-4
    else:
        t_final = 12.0
        dt = 0.008
        renorm = 18
        jac_eps = 6e-5

    print(
        "ex150: lambda_max pairwise vs coulomb_3b (Burrau), beta=gamma="
        f"{beta} soft={soft}"
    )
    print(f"  t_final={t_final} dt={dt} renorm_every={renorm}")

    t0 = time.perf_counter()
    lam_p, sp = largest_lyapunov_exponent_benettin(
        pos0,
        vel0,
        C,
        acc_pair,
        t_final=t_final,
        dt=dt,
        renorm_every=renorm,
        jac_eps=jac_eps,
        rng=np.random.default_rng(1501),
    )
    t1 = time.perf_counter()
    lam_3, s3 = largest_lyapunov_exponent_benettin(
        pos0,
        vel0,
        C,
        acc_3b,
        t_final=t_final,
        dt=dt,
        renorm_every=renorm,
        jac_eps=jac_eps,
        rng=np.random.default_rng(1502),
    )
    t2 = time.perf_counter()

    print(f"  pairwise     lambda={lam_p:.6f}  steps={sp}  cpu={t1 - t0:.2f}s")
    print(f"  coulomb_3b   lambda={lam_3:.6f}  steps={s3}  cpu={t2 - t1:.2f}s")
    print(f"  delta(lambda) = {lam_3 - lam_p:.6f}")

    ok = (
        np.isfinite(lam_p)
        and np.isfinite(lam_3)
        and lam_p > 0.02
        and abs(lam_3 - lam_p) > 1e-4
    )
    if not ok:
        raise SystemExit(1)
    print("  PASS (finite, chaotic-scale lambda, 3b != pairwise)")


if __name__ == "__main__":
    main()
