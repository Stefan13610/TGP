#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex151_lyapunov_scan_beta_pairwise.py
====================================
P1 (cd.): skan ``beta`` (przy ``gamma=beta``) dla **TGP pairwise** —
``lambda_max`` w problemie pitagorejskim (te same IC co ex148).

Pokazuje wrazliwosc chaosu na sile czlonu ~1/r^2 w V2 (przy ustalonym soft).

Uruchomienie: ``python ex151_lyapunov_scan_beta_pairwise.py [--quick] [--out PATH]``
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
        help="opcjonalny CSV: beta,lambda_max (domyslnie brak pliku)",
    )
    args = parser.parse_args()

    soft = 1e-6
    pos0, vel0, C = pythagorean_three_body_burrau()

    if args.quick:
        betas = np.array([0.03, 0.06, 0.10], dtype=float)
        t_final = 9.0
        dt = 0.01
        renorm = 18
        jac_eps = 9e-5
    else:
        betas = np.linspace(0.02, 0.14, 7)
        t_final = 14.0
        dt = 0.008
        renorm = 20
        jac_eps = 7e-5

    print("ex151: lambda_max vs beta (pairwise TGP, Burrau IC)")
    print(f"  softening={soft}  t_final={t_final} dt={dt}")
    print("  beta     lambda_max")
    last_lam = float("nan")
    rows: list[tuple[float, float]] = []
    for i, b in enumerate(betas):
        beta = float(b)
        gamma = beta

        def acc_fn(p: np.ndarray, c: np.ndarray, b=beta, g=gamma) -> np.ndarray:
            return analytical_forces_tgp_pairwise(p, c, b, g, soft) / c[:, None]

        lam, _st = largest_lyapunov_exponent_benettin(
            pos0,
            vel0,
            C,
            acc_fn,
            t_final=t_final,
            dt=dt,
            renorm_every=renorm,
            jac_eps=jac_eps,
            rng=np.random.default_rng(15100 + i),
        )
        print(f"  {beta:6.4f}   {lam: .6f}")
        rows.append((beta, float(lam)))
        last_lam = lam

    if args.out:
        with open(args.out, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["beta", "lambda_max"])
            for b, lam in rows:
                w.writerow([f"{b:.16e}", f"{lam:.16e}"])
        print(f"  wrote {args.out}")

    ok = np.isfinite(last_lam) and last_lam > 0.02
    if not ok:
        raise SystemExit(1)
    print("  PASS (last lambda finite and > 0.02)")


if __name__ == "__main__":
    main()
