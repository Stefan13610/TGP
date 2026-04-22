#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex153_lyapunov_spectrum_sum_newton.py
=====================================
P1 (cd.): pelne spektrum (k = 6N = 18 dla N=3) dla **Newtona** + Burrau IC.

Dla autonomicznego potencjalowego ukladu Hamiltona suma wszystkich wykładników
Lapunowa (srednia czasowa) powinna dążyć do **0** (zachowanie objetosci /
slad macierzy monodromii symplektycznej).

Test regresyjny: ``|sum(lambda)| < tol``. Domyślnie **leapfrog + styczny**
(symplektyczna baza) — suma pierwszych ``6N`` wykładników bywa blisko **0**;
``--rk4`` włącza RK4 na rozszerzonym ODE (często ``|Σλ|`` większe przy skończonym ``T``).

Uruchomienie:
``python ex153_lyapunov_spectrum_sum_newton.py [--quick] [--rk4] [--out PATH]``
"""

from __future__ import annotations

import argparse
import csv
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
        help="RK4 na rozszerzonym stanie zamiast leapfrog+styczny (porownanie)",
    )
    parser.add_argument(
        "--out",
        default="",
        help="CSV: kolumny j, lambda (oraz wiersz sum)",
    )
    args = parser.parse_args()

    soft = 1e-6
    G = 1.0
    pos0, vel0, C = pythagorean_three_body_burrau()
    n = len(C)
    k = 6 * n

    def acc_newton(p: np.ndarray, c: np.ndarray) -> np.ndarray:
        return forces_newton(p, c, G=G, softening=soft) / c[:, None]

    if args.quick:
        t_final = 10.0
        dt = 0.012
        renorm = 15
        jac_eps = 9e-5
        tol_lf = 0.06
        tol_rk = 0.35
    else:
        t_final = 22.0
        dt = 0.008
        renorm = 20
        jac_eps = 6e-5
        tol_lf = 0.05
        tol_rk = 0.28

    def jac_newton(p: np.ndarray, c: np.ndarray) -> np.ndarray:
        return acceleration_jacobian_newton_softened(
            p, c, G=G, softening=soft
        )

    rng = np.random.default_rng(153)
    use_rk4 = bool(args.rk4)
    tol = tol_rk if use_rk4 else tol_lf

    if not use_rk4:
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
    else:
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

    ssum = float(np.sum(spec))
    print("ex153: full Lyapunov spectrum (Newton), N=3, k=6N=18")
    print(
        f"  t_final={t_final} dt={dt} steps={steps} "
        f"integrator={integ} jacobian=Newton analytical (soft={soft})"
    )
    print(f"  sum(lambda) = {ssum:.6f}  (expect ~0 Hamiltonian)")
    print(f"  lambda[0]={spec[0]:.6f}  lambda[-1]={spec[-1]:.6f}")

    if args.out:
        with open(args.out, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["j", "lambda"])
            for j in range(k):
                w.writerow([str(j), f"{float(spec[j]):.16e}"])
            w.writerow(["sum", f"{ssum:.16e}"])
        print(f"  wrote {args.out}")

    ok = np.all(np.isfinite(spec)) and abs(ssum) < tol
    if not ok:
        raise SystemExit(1)
    print(f"  PASS (|sum| < {tol}, {'RK4' if use_rk4 else 'leapfrog'})")


if __name__ == "__main__":
    main()
