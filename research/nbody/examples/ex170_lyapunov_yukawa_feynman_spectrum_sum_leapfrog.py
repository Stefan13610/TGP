#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex170_lyapunov_yukawa_feynman_spectrum_sum_leapfrog.py
======================================================
P1 (cd.): pelne spektrum Lyapunova **k = 6N** (N=3 => 18) dla
``yukawa_feynman`` + **leapfrog+tangent** + **analityczny**
``acceleration_jacobian_yukawa_feynman_analytic``.

Dla ukladu Hamiltona (symplektyczny tangent + jawny symplektyczny Jacobian)
srednia ``sum(lambda)`` po calej bazie powinna byc blisko **0** (cf. ``ex153``
dla Newtona). Porownaj ``ex171`` (``coulomb_3b`` + FD ``J``), ``ex172`` (``pairwise`` + analityczny ``J``) - te same ``dt``/``t``.

Uruchomienie:
``python ex170_lyapunov_yukawa_feynman_spectrum_sum_leapfrog.py [--quick] [--n-quad N]``
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
from nbody.lyapunov import (
    acceleration_jacobian_yukawa_feynman_analytic,
    lyapunov_spectrum_benettin_leapfrog,
    pythagorean_three_body_burrau,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true")
    parser.add_argument(
        "--n-quad",
        type=int,
        default=0,
        help="n_quad Feynman (0 = domyslne wg trybu)",
    )
    args = parser.parse_args()

    beta = gamma = 0.08
    soft = 1e-6
    seed = 1700
    pos0, vel0, C = pythagorean_three_body_burrau()
    n = len(C)
    k = 6 * n

    if args.quick:
        t_final = 1.5
        dt = 0.06
        renorm = 5
        jac_eps = 1e-5
        n_quad = 12
        tol_sum = 1e-4
    else:
        t_final = 4.0
        dt = 0.04
        renorm = 10
        jac_eps = 1e-5
        n_quad = 14
        tol_sum = 1e-5
    if args.n_quad > 0:
        n_quad = int(args.n_quad)

    acc_y, _ = build_tgp_integration_pair(
        "yukawa_feynman",
        beta=beta,
        gamma=gamma,
        softening=soft,
        include_3body=True,
        n_quad_feynman=n_quad,
    )

    def jac_fn(p: np.ndarray, c: np.ndarray, nq: int = n_quad) -> np.ndarray:
        return acceleration_jacobian_yukawa_feynman_analytic(
            p,
            c,
            beta=beta,
            gamma=gamma,
            softening=soft,
            n_quad_feynman=nq,
        )

    spec, steps = lyapunov_spectrum_benettin_leapfrog(
        pos0,
        vel0,
        C,
        acc_y,
        n_exponents=k,
        t_final=t_final,
        dt=dt,
        renorm_every=renorm,
        jac_eps=jac_eps,
        position_jacobian_fn=jac_fn,
        rng=np.random.default_rng(seed),
    )

    ssum = float(np.sum(spec))
    print(
        "ex170: Lyapunov spectrum sum (yukawa_feynman, leapfrog+tangent, "
        f"analytic J, n_quad={n_quad})"
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
