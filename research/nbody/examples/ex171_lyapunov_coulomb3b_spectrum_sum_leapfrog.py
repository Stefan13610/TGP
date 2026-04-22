#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex171_lyapunov_coulomb3b_spectrum_sum_leapfrog.py
=================================================
P1 (cd.): pelne spektrum Lyapunova **k = 6N** (N=3 => 18) dla backendu
``coulomb_3b`` (``V_2`` + zamkniety irreducible ``I`` Coulomb) + **leapfrog+tangent** +
Jacobian z **roznic skonczonych** (``position_jacobian_fn=None``).

Regresja: **|sum(lambda)|** blisko zera (por. ``ex153`` Newton, ``ex170``
``yukawa_feynman`` + analityczny $J$, ``ex172`` ``pairwise``). Uklad jest Hamiltonowski; symplektyczny
LF + przyblizony $J$ przez FD zwykle utrzymuje bardzo mala sume na skonczonym
czasie.

Uruchomienie:
``python ex171_lyapunov_coulomb3b_spectrum_sum_leapfrog.py [--quick]``
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np

_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

from nbody.dynamics_backends import build_tgp_integration_pair
from nbody.lyapunov import lyapunov_spectrum_benettin_leapfrog, pythagorean_three_body_burrau


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true")
    args = parser.parse_args()

    beta = gamma = 0.08
    soft = 1e-6
    seed = 1710
    pos0, vel0, C = pythagorean_three_body_burrau()
    n = len(C)
    k = 6 * n

    if args.quick:
        # Short horizon + coarse dt → spectrum sum drifts O(1e-3);
        # on different platforms (Ubuntu vs Windows) FP rounding
        # changes chaotic trajectory → accept larger tolerance
        t_final = 1.5
        dt = 0.06
        renorm = 5
        jac_eps = 1e-5
        tol_sum = 2e-3
    else:
        t_final = 4.0
        dt = 0.04
        renorm = 10
        jac_eps = 1e-5
        tol_sum = 1e-5

    acc_c, _ = build_tgp_integration_pair(
        "coulomb_3b",
        beta=beta,
        gamma=gamma,
        softening=soft,
        include_3body=True,
    )

    spec, steps = lyapunov_spectrum_benettin_leapfrog(
        pos0,
        vel0,
        C,
        acc_c,
        n_exponents=k,
        t_final=t_final,
        dt=dt,
        renorm_every=renorm,
        jac_eps=jac_eps,
        position_jacobian_fn=None,
        rng=np.random.default_rng(seed),
    )

    ssum = float(np.sum(spec))
    print(
        "ex171: Lyapunov spectrum sum (coulomb_3b, leapfrog+tangent, FD J)"
    )
    print(
        f"  t_final={t_final} dt={dt} renorm_every={renorm} steps={steps} "
        f"k={k} seed={seed}"
    )
    print(f"  sum(lambda) = {ssum:.6e}  (tol |sum| < {tol_sum})")
    print(f"  lambda[0]={spec[0]:.6f}  lambda[-1]={spec[-1]:.6f}")

    ok = (
        np.all(np.isfinite(spec))
        and abs(ssum) < tol_sum
        and float(spec[0]) > 0.02
    )
    if not ok:
        raise SystemExit(1)
    print("  PASS")


if __name__ == "__main__":
    main()
