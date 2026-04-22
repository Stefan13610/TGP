#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex155_lyapunov_grid_beta_softening_csv.py
=========================================
P1 (cd.): siatka ``(beta, softening)`` przy ``gamma=beta`` — TGP **pairwise**,
``lambda_max`` (Burrau IC). Wynik: CSV w ``_outputs/``.

Uzyteczne jako tabelaryczny przeglad wrazliwosci na parametry 2-cialowe
(przy ustalonym IC chaosu).

Uruchomienie: ``python ex155_lyapunov_grid_beta_softening_csv.py [--quick] [--out PATH]``
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

from nbody.dynamics_v2 import analytical_forces_tgp_pairwise
from nbody.lyapunov import largest_lyapunov_exponent_benettin, pythagorean_three_body_burrau


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true")
    parser.add_argument(
        "--out",
        default="",
        help="sciezka CSV (domyslnie examples/_outputs/ex155_beta_softening_grid.csv)",
    )
    args = parser.parse_args()

    pos0, vel0, C = pythagorean_three_body_burrau()

    if args.quick:
        betas = np.array([0.04, 0.08, 0.12], dtype=float)
        soft_list = np.array([3e-7, 1e-6, 3e-6], dtype=float)
        t_final = 6.0
        dt = 0.012
        renorm = 12
        jac_eps = 1e-4
    else:
        betas = np.linspace(0.03, 0.14, 5)
        soft_list = np.logspace(-7.2, -4.5, 5)
        t_final = 10.0
        dt = 0.01
        renorm = 16
        jac_eps = 8e-5

    out_dir = os.path.join(os.path.dirname(__file__), "_outputs")
    os.makedirs(out_dir, exist_ok=True)
    out_path = args.out or os.path.join(out_dir, "ex155_beta_softening_grid.csv")

    print("ex155: lambda_max grid (beta x softening), pairwise TGP, Burrau")
    print(f"  t_final={t_final} dt={dt}  |betas|={len(betas)} |soft|={len(soft_list)}")

    seed = 15500
    rows: list[tuple[float, float, float]] = []
    for bi, beta in enumerate(betas):
        for si, soft in enumerate(soft_list):
            b = float(beta)
            g = b
            s = float(soft)

            def acc_fn(
                p: np.ndarray, c: np.ndarray, bb=b, gg=g, sf=s
            ) -> np.ndarray:
                return analytical_forces_tgp_pairwise(p, c, bb, gg, sf) / c[:, None]

            lam, _st = largest_lyapunov_exponent_benettin(
                pos0,
                vel0,
                C,
                acc_fn,
                t_final=t_final,
                dt=dt,
                renorm_every=renorm,
                jac_eps=jac_eps,
                rng=np.random.default_rng(seed + bi * 31 + si),
            )
            rows.append((b, s, float(lam)))
            print(f"  beta={b:.4f} soft={s:.2e}  lambda={lam:.6f}")

    with open(out_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["beta", "softening", "lambda_max"])
        for b, s, lam in rows:
            w.writerow([f"{b:.16e}", f"{s:.16e}", f"{lam:.16e}"])

    last = rows[-1][2]
    ok = np.isfinite(last) and last > 0.02
    if not ok:
        raise SystemExit(1)
    print(f"ex155: wrote {out_path}  PASS")


if __name__ == "__main__":
    main()
