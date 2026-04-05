#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex161_lyapunov_pairwise_beta_d_rep_scan.py
==========================================
P1 (cd.): skan **β** (przy ``gamma=beta``) dla **TGP pairwise** z **jawnym
Jacobianem** ``acceleration_jacobian_tgp_pairwise_softened`` (leapfrog+tangent).

Dopisuje kolumnę **d_rep**: wewnętrzny pierwiastek siły 2-ciałowej dla
``C_ref = min(C_i)`` (masy Burrau 3,4,5 → ``C_ref=3``), z ``pairwise.force_zeros_2body``.
Wymaga ``β > (9/2) C_ref`` w jednostkach modułu, stąd **większe β** niż w ex151.

Uruchomienie: ``python ex161_lyapunov_pairwise_beta_d_rep_scan.py [--quick] [--out PATH]``
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
from nbody.lyapunov import (
    acceleration_jacobian_tgp_pairwise_softened,
    largest_lyapunov_exponent_benettin_leapfrog,
    pythagorean_three_body_burrau,
)
from nbody.pairwise import force_zeros_2body


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true")
    parser.add_argument(
        "--out",
        default="",
        help="opcjonalny CSV: beta,d_rep,lambda_max",
    )
    args = parser.parse_args()

    soft = 1e-6
    pos0, vel0, C = pythagorean_three_body_burrau()
    c_ref = float(np.min(C))

    if args.quick:
        betas = np.array([14.0, 17.0, 20.0], dtype=float)
        t_final = 4.0
        dt = 0.02
        renorm = 10
    else:
        betas = np.linspace(13.5, 22.0, 6)
        t_final = 6.0
        dt = 0.015
        renorm = 12

    print("ex161: lambda_max vs beta, pairwise TGP + analytic J (Burrau IC)")
    print(f"  softening={soft}  C_ref=min(C)={c_ref}  t_final={t_final} dt={dt}")
    print("  beta     d_rep(inner)   lambda_max")
    last_lam = float("nan")
    rows: list[tuple[float, float, float]] = []
    for i, beta in enumerate(betas):
        beta = float(beta)
        gamma = beta
        roots = force_zeros_2body(c_ref, beta, gamma)
        if roots is None:
            d_rep = float("nan")
        else:
            d_rep = float(roots[0])

        def acc_fn(p: np.ndarray, c: np.ndarray, b=beta, g=gamma) -> np.ndarray:
            return analytical_forces_tgp_pairwise(p, c, b, g, soft) / c[:, None]

        def jac_fn(
            p: np.ndarray, c: np.ndarray, b=beta, g=gamma
        ) -> np.ndarray:
            return acceleration_jacobian_tgp_pairwise_softened(
                p, c, beta=b, gamma=g, softening=soft
            )

        lam, _st = largest_lyapunov_exponent_benettin_leapfrog(
            pos0,
            vel0,
            C,
            acc_fn,
            t_final=t_final,
            dt=dt,
            renorm_every=renorm,
            jac_eps=1e-5,
            position_jacobian_fn=jac_fn,
            rng=np.random.default_rng(16100 + i),
        )
        dr_s = f"{d_rep:10.4f}" if np.isfinite(d_rep) else "       nan"
        print(f"  {beta:6.2f}   {dr_s}   {lam: .6f}")
        rows.append((beta, d_rep, float(lam)))
        last_lam = lam

    if args.out:
        with open(args.out, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["beta", "d_rep_inner", "lambda_max"])
            for b, dr, lam in rows:
                w.writerow(
                    [
                        f"{b:.16e}",
                        "" if not np.isfinite(dr) else f"{dr:.16e}",
                        f"{lam:.16e}",
                    ]
                )
        print(f"  wrote {args.out}")

    ok = np.isfinite(last_lam) and last_lam > 0.02
    if not ok:
        raise SystemExit(1)
    print("  PASS (last lambda finite and > 0.02)")


if __name__ == "__main__":
    main()
