#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex159_lyapunov_coulomb3b_rk4_vs_leapfrog.py
===========================================
P1 (cd.): ``lambda_max`` dla backendu **coulomb_3b** (``eom_tgp``, V2 + I Coulomb
w V3), Burrau N=3 — Jacobian ``a`` wzgledem ``x`` przez **FD**, porownanie
RK4+tangent vs leapfrog+tangent.

Kosztowne (wiele wywolan ``acc_fn`` na krok); ``--quick`` utrzymuje krotki horyzont.

Uruchomienie: ``python ex159_lyapunov_coulomb3b_rk4_vs_leapfrog.py [--quick]``
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
from nbody.lyapunov import (
    largest_lyapunov_exponent_benettin,
    largest_lyapunov_exponent_benettin_leapfrog,
    pythagorean_three_body_burrau,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true")
    args = parser.parse_args()

    beta = gamma = 0.08
    soft = 1e-6
    pos0, vel0, C = pythagorean_three_body_burrau()

    acc_fn, _ = build_tgp_integration_pair(
        "coulomb_3b",
        beta=beta,
        gamma=gamma,
        softening=soft,
        include_3body=True,
    )

    if args.quick:
        # Krotki czas daje duzy rozrzut estymatora; 3.5/0.04 ~ zgodnosc |Δλ| ~0.2
        t_final = 3.5
        dt = 0.04
        renorm = 8
        jac_eps = 1e-4
        delta_tol = 0.45
    else:
        t_final = 6.0
        dt = 0.028
        renorm = 10
        jac_eps = 8e-5
        delta_tol = 0.22

    print(
        "ex159: lambda_max coulomb_3b (Burrau), FD Jacobian, "
        "RK4 tangent vs leapfrog tangent"
    )
    print(f"  beta=gamma={beta} soft={soft} t_final={t_final} dt={dt}")

    t0 = time.perf_counter()
    lam_rk, sr = largest_lyapunov_exponent_benettin(
        pos0,
        vel0,
        C,
        acc_fn,
        t_final=t_final,
        dt=dt,
        renorm_every=renorm,
        jac_eps=jac_eps,
        position_jacobian_fn=None,
        rng=np.random.default_rng(1591),
    )
    t1 = time.perf_counter()
    lam_lf, sl = largest_lyapunov_exponent_benettin_leapfrog(
        pos0,
        vel0,
        C,
        acc_fn,
        t_final=t_final,
        dt=dt,
        renorm_every=renorm,
        jac_eps=jac_eps,
        position_jacobian_fn=None,
        rng=np.random.default_rng(1592),
    )
    t2 = time.perf_counter()

    print(f"  steps_rk={sr} cpu_rk={t1 - t0:.2f}s")
    print(f"  steps_lf={sl} cpu_lf={t2 - t1:.2f}s")
    print(f"  lambda_max RK4+tangent      = {lam_rk:.6f}")
    print(f"  lambda_max leapfrog+tangent  = {lam_lf:.6f}")
    print(f"  |delta| = {abs(lam_rk - lam_lf):.6f}")

    ok = (
        np.isfinite(lam_rk)
        and np.isfinite(lam_lf)
        and lam_rk > 0.02
        and lam_lf > 0.02
        and abs(lam_rk - lam_lf) < delta_tol
    )
    if not ok:
        raise SystemExit(1)
    print(f"  PASS (finite, chaotic-scale, |delta| < {delta_tol})")


if __name__ == "__main__":
    main()
