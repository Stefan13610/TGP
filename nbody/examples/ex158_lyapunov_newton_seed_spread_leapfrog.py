#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex158_lyapunov_newton_seed_spread_leapfrog.py
=============================================
P1 (cd.): rozrzut estymatora ``lambda_max`` (Benettin + leapfrog + styczna)
względem **ziarna RNG** inicjalizacji wektora stycznego — ta sama dynamika
bazowa (Newton, Burrau), jawny Jacobian.

Przy skonczonym ``T`` najwiekszy wykladnik jest wyznaczany z losowego kierunku
stycznego; rozrzut miedzy seedami pokazuje szum estymatora (nie fizyczna
niepewnosc samego ``lambda_1``).

Uruchomienie: ``python ex158_lyapunov_newton_seed_spread_leapfrog.py [--quick] [--out PATH]``
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
    largest_lyapunov_exponent_benettin_leapfrog,
    pythagorean_three_body_burrau,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true")
    parser.add_argument(
        "--out",
        default="",
        help="opcjonalny CSV: seed, lambda_max",
    )
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
        t_final = 7.0
        dt = 0.014
        renorm = 14
        jac_eps = 1e-4
        seeds = list(range(15800, 15812))
        std_tol = 0.22
    else:
        t_final = 12.0
        dt = 0.01
        renorm = 16
        jac_eps = 8e-5
        seeds = list(range(15800, 15820))
        std_tol = 0.20

    lams: list[float] = []
    print(
        "ex158: lambda_max spread vs tangent RNG seed "
        "(Newton, leapfrog, analytical J, Burrau)"
    )
    print(f"  t_final={t_final} dt={dt} n_seeds={len(seeds)}")
    for sd in seeds:
        lam, _st = largest_lyapunov_exponent_benettin_leapfrog(
            pos0,
            vel0,
            C,
            acc_newton,
            t_final=t_final,
            dt=dt,
            renorm_every=renorm,
            jac_eps=jac_eps,
            position_jacobian_fn=jac_newton,
            rng=np.random.default_rng(int(sd)),
        )
        lams.append(float(lam))
        print(f"  seed={sd}  lambda_max={lam:.6f}")

    arr = np.array(lams, dtype=float)
    mean = float(np.mean(arr))
    std = float(np.std(arr, ddof=0))
    print(f"  mean={mean:.6f}  std={std:.6f}  min={arr.min():.6f}  max={arr.max():.6f}")

    if args.out:
        with open(args.out, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["seed", "lambda_max"])
            for sd, lam in zip(seeds, lams):
                w.writerow([str(sd), f"{lam:.16e}"])
        print(f"  wrote {args.out}")

    ok = np.all(np.isfinite(arr)) and mean > 0.02 and std < std_tol
    if not ok:
        raise SystemExit(1)
    print(f"  PASS (finite, mean>0.02, std<{std_tol})")


if __name__ == "__main__":
    main()
