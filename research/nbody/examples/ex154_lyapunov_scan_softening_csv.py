#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex154_lyapunov_scan_softening_csv.py
====================================
P1 (cd.): skan ``softening`` (regulator numeryczny w V2) dla TGP **pairwise** —
``lambda_max``, IC Burrau. Wynik: CSV w ``examples/_outputs/``.

Uruchomienie: ``python ex154_lyapunov_scan_softening_csv.py [--quick]``

Opcja ``--out PATH`` nadpisuje domyślną ścieżkę pliku wyjściowego.
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
        help="sciezka CSV (domyslnie examples/_outputs/ex154_softening_scan.csv)",
    )
    args = parser.parse_args()

    beta = gamma = 0.08
    pos0, vel0, C = pythagorean_three_body_burrau()

    if args.quick:
        soft_list = np.array([1e-7, 3e-6, 1e-5], dtype=float)
        t_final = 8.0
        dt = 0.01
        renorm = 16
        jac_eps = 9e-5
    else:
        soft_list = np.logspace(-7.5, -4.0, 8)
        t_final = 12.0
        dt = 0.008
        renorm = 20
        jac_eps = 7e-5

    out_dir = os.path.join(os.path.dirname(__file__), "_outputs")
    os.makedirs(out_dir, exist_ok=True)
    out_path = args.out or os.path.join(out_dir, "ex154_softening_scan.csv")

    rows: list[tuple[float, float]] = []
    for i, soft in enumerate(soft_list):
        eps = float(soft)

        def acc_fn(
            p: np.ndarray, c: np.ndarray, bb=beta, gg=gamma, s=eps
        ) -> np.ndarray:
            return analytical_forces_tgp_pairwise(p, c, bb, gg, s) / c[:, None]

        lam, _st = largest_lyapunov_exponent_benettin(
            pos0,
            vel0,
            C,
            acc_fn,
            t_final=t_final,
            dt=dt,
            renorm_every=renorm,
            jac_eps=jac_eps,
            rng=np.random.default_rng(15400 + i),
        )
        rows.append((eps, lam))
        print(f"  softening={eps:.3e}  lambda_max={lam:.6f}")

    with open(out_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["softening", "lambda_max"])
        for eps, lam in rows:
            w.writerow([f"{eps:.16e}", f"{lam:.16e}"])

    last = rows[-1][1]
    ok = np.isfinite(last) and last > 0.02
    if not ok:
        raise SystemExit(1)
    print(f"ex154: wrote {out_path}  PASS")


if __name__ == "__main__":
    main()
