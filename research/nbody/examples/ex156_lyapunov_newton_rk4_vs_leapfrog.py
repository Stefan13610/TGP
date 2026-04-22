#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex156_lyapunov_newton_rk4_vs_leapfrog.py
========================================
P1 (cd.): porownanie ``lambda_max`` (Benettin, jeden wektor styczny) dla
**Newtona** + Burrau: orbita + styczna przez **RK4** vs **leapfrog**.

Te same ``t_final``, ``dt``, ``renorm``, jawny Jacobian Newtona; rozne
ziarna RNG dla wektora stycznego (obie trajektorie bazowe sa identyczne
co do algorytmu — roznia sie dyskretna mapa styczna).

Uruchomienie: ``python ex156_lyapunov_newton_rk4_vs_leapfrog.py [--quick]``
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np

_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

from nbody.dynamics_v2 import forces_newton
from nbody.lyapunov import (
    acceleration_jacobian_newton_softened,
    largest_lyapunov_exponent_benettin,
    largest_lyapunov_exponent_benettin_leapfrog,
    pythagorean_three_body_burrau,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true")
    args = parser.parse_args()

    soft = 1e-6
    G = 1.0
    pos0, vel0, C = pythagorean_three_body_burrau()

    def acc_newton(p: np.ndarray, c: np.ndarray) -> np.ndarray:
        return forces_newton(p, c, G=G, softening=soft) / c[:, None]

    def jac_newton(p: np.ndarray, c: np.ndarray) -> np.ndarray:
        return acceleration_jacobian_newton_softened(
            p, c, G=G, softening=soft
        )

    if args.quick:
        t_final = 10.0
        dt = 0.012
        renorm = 15
        jac_eps = 9e-5
        delta_tol = 0.55
    else:
        t_final = 18.0
        dt = 0.008
        renorm = 18
        jac_eps = 7e-5
        delta_tol = 0.45

    lam_rk, sr = largest_lyapunov_exponent_benettin(
        pos0,
        vel0,
        C,
        acc_newton,
        t_final=t_final,
        dt=dt,
        renorm_every=renorm,
        jac_eps=jac_eps,
        position_jacobian_fn=jac_newton,
        rng=np.random.default_rng(1561),
    )
    lam_lf, sl = largest_lyapunov_exponent_benettin_leapfrog(
        pos0,
        vel0,
        C,
        acc_newton,
        t_final=t_final,
        dt=dt,
        renorm_every=renorm,
        jac_eps=jac_eps,
        position_jacobian_fn=jac_newton,
        rng=np.random.default_rng(1562),
    )

    print("ex156: lambda_max Newton (Burrau), RK4 tangent vs leapfrog tangent")
    print(f"  t_final={t_final} dt={dt} steps_rk={sr} steps_lf={sl}")
    print(f"  lambda_max RK4+tangent     = {lam_rk:.6f}")
    print(f"  lambda_max leapfrog+tangent = {lam_lf:.6f}")
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
