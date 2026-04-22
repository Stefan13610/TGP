#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex149_lyapunov_spectrum_newton.py
=================================
P1 (cd.): pierwsze k wykładników Lapunowa (Benettin + QR) dla **Newtona**
w problemie pitagorejskim — bez TGP (referencja numeryczna).

Uwaga: przy skróconym ``t_final`` (``--quick``) kolejność ``lambda[1],...``
może jeszcze nie odpowiadać asymptotycznemu uporządkowaniu — rośnie
``t_final`` w trybie pełnym.

Domyślnie **leapfrog+tangent** (jak ``ex153``); ``--rk4`` — spektrum z RK4.

Uruchomienie: ``python ex149_lyapunov_spectrum_newton.py [--quick] [--rk4]``
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
    lyapunov_spectrum_benettin,
    lyapunov_spectrum_benettin_leapfrog,
    pythagorean_three_body_burrau,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true")
    parser.add_argument(
        "--rk4",
        action="store_true",
        help="RK4+tangent zamiast leapfrog+tangent",
    )
    args = parser.parse_args()

    soft = 1e-6
    G = 1.0
    pos0, vel0, C = pythagorean_three_body_burrau()

    def acc_newton(p: np.ndarray, c: np.ndarray) -> np.ndarray:
        return forces_newton(p, c, G=G, softening=soft) / c[:, None]

    k = 6
    if args.quick:
        t_final = 14.0
        dt = 0.008
        renorm = 20
        jac_eps = 8e-5
    else:
        t_final = 48.0
        dt = 0.004
        renorm = 22
        jac_eps = 5e-5

    def jac_newton(p: np.ndarray, c: np.ndarray) -> np.ndarray:
        return acceleration_jacobian_newton_softened(
            p, c, G=G, softening=soft
        )

    rng = np.random.default_rng(149)
    use_rk4 = bool(args.rk4)
    if use_rk4:
        spec, steps = lyapunov_spectrum_benettin(
            pos0,
            vel0,
            C,
            acc_newton,
            n_exponents=k,
            t_final=t_final,
            dt=dt,
            renorm_every=renorm,
            jac_eps=jac_eps,
            position_jacobian_fn=jac_newton,
            rng=rng,
        )
        integ = "RK4+tangent"
    else:
        spec, steps = lyapunov_spectrum_benettin_leapfrog(
            pos0,
            vel0,
            C,
            acc_newton,
            n_exponents=k,
            t_final=t_final,
            dt=dt,
            renorm_every=renorm,
            jac_eps=jac_eps,
            position_jacobian_fn=jac_newton,
            rng=rng,
        )
        integ = "leapfrog+tangent"

    print("ex149: Lyapunov spectrum (Newton, Burrau IC), first k=6 exponents")
    print(
        f"  t_final={t_final} dt={dt} steps={steps} "
        f"integrator={integ} jacobian=Newton analytical"
    )
    for j in range(k):
        print(f"  lambda[{j}] = {spec[j]: .6f}")
    ok = np.all(np.isfinite(spec)) and float(spec[0]) > 0.05
    if not ok:
        raise SystemExit(1)
    print("  PASS (finite, lambda[0] > 0.05)")


if __name__ == "__main__":
    main()
