#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex166_lyapunov_yukawa_feynman_convergence_grid_csv.py
=====================================================
P1 (cd.): siatka ``(dt, n_quad)`` dla ``yukawa_feynman`` + **analityczny** ``J`` —
wszystkie przebiegi z **tym samym** ``t_final`` i **tym samym** seedem stycznym,
zeby porownac wplyw dyskretyzacji i kwadratury Feynmana na ``lambda_max``.

Wynik: CSV (domyslnie ``examples/_outputs/ex166_yukawa_convergence_grid.csv``).

``--quick``: mniejsza siatka + test: wzgledny rozrzut ``(max-min)/mean`` ponizej
progu. Tryb pelny: wieksza siatka; PASS gdy wszystkie ``lambda`` skonczone i
``> 0.02`` (bez progu rozrzutu - przy dluzszym ``t_final`` roznice moga byc duze).

Uruchomienie: ``python ex166_lyapunov_yukawa_feynman_convergence_grid_csv.py [--quick] [--out PATH]``
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

from nbody.dynamics_backends import build_tgp_integration_pair
from nbody.lyapunov import (
    acceleration_jacobian_yukawa_feynman_analytic,
    largest_lyapunov_exponent_benettin_leapfrog,
    pythagorean_three_body_burrau,
)


def _run_cell(
    *,
    t_final: float,
    dt: float,
    n_quad: int,
    renorm: int,
    seed: int,
) -> tuple[float, int]:
    beta = gamma = 0.08
    soft = 1e-6
    pos0, vel0, C = pythagorean_three_body_burrau()

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

    return largest_lyapunov_exponent_benettin_leapfrog(
        pos0,
        vel0,
        C,
        acc_y,
        t_final=t_final,
        dt=dt,
        renorm_every=renorm,
        jac_eps=1e-5,
        position_jacobian_fn=jac_fn,
        rng=np.random.default_rng(seed),
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true")
    parser.add_argument(
        "--out",
        default="",
        help="sciezka CSV (domyslnie examples/_outputs/ex166_yukawa_convergence_grid.csv)",
    )
    args = parser.parse_args()

    seed = 1660
    if args.quick:
        t_final = 3.5
        renorm = 10
        dts = [0.06, 0.045, 0.03]
        n_quads = [12, 18]
        spread_tol = 0.32
    else:
        t_final = 5.0
        renorm = 12
        dts = [0.05, 0.042, 0.034, 0.028]
        n_quads = [14, 20, 26]
        spread_tol = None

    out_dir = os.path.join(os.path.dirname(__file__), "_outputs")
    os.makedirs(out_dir, exist_ok=True)
    out_path = args.out or os.path.join(out_dir, "ex166_yukawa_convergence_grid.csv")

    rows: list[tuple[float, float, int, int, int, float, int]] = []
    lams: list[float] = []

    print("ex166: grid dt x n_quad (yukawa_feynman, analytic J, same seed)")
    print(f"  t_final={t_final} renorm_every={renorm} seed={seed}")
    print("  dt       n_quad   lambda_max   steps")

    for dt in dts:
        for nq in n_quads:
            lam, steps = _run_cell(
                t_final=t_final,
                dt=float(dt),
                n_quad=int(nq),
                renorm=renorm,
                seed=seed,
            )
            rows.append((t_final, float(dt), int(nq), renorm, seed, float(lam), int(steps)))
            lams.append(float(lam))
            print(f"  {float(dt):.4f}   {int(nq):3d}     {lam: .6f}   {steps}")

    with open(out_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "t_final",
                "dt",
                "n_quad",
                "renorm_every",
                "seed",
                "lambda_max",
                "steps",
            ]
        )
        for row in rows:
            w.writerow(
                [
                    f"{row[0]:.16e}",
                    f"{row[1]:.16e}",
                    row[2],
                    row[3],
                    row[4],
                    f"{row[5]:.16e}",
                    row[6],
                ]
            )

    print(f"  wrote {out_path}")

    ok = all(np.isfinite(l) and l > 0.02 for l in lams)
    if spread_tol is not None and lams:
        mean_l = float(np.mean(lams))
        spread = (max(lams) - min(lams)) / mean_l if mean_l > 0 else float("inf")
        print(f"  relative spread (max-min)/mean = {spread:.4f}  (tol {spread_tol})")
        ok = ok and spread < spread_tol

    if not ok:
        raise SystemExit(1)
    print("  PASS")


if __name__ == "__main__":
    main()
