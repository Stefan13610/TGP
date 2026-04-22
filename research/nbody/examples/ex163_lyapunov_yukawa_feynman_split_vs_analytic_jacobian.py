#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex163_lyapunov_yukawa_feynman_split_vs_analytic_jacobian.py
===========================================================
P1 (cd.): **yukawa_feynman** + leapfrog+tangent — porównanie ``lambda_max`` przy
Jacobianie **split** (FD na ``V_3``) vs **w pełni analitycznym**
(``acceleration_jacobian_yukawa_feynman_analytic``).

Te same IC co ``ex160`` / ``ex162``.

Uruchomienie: ``python ex163_lyapunov_yukawa_feynman_split_vs_analytic_jacobian.py [--quick] [--n-quad N]``
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
    acceleration_jacobian_yukawa_feynman_split,
    largest_lyapunov_exponent_benettin_leapfrog,
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
    pos0, vel0, C = pythagorean_three_body_burrau()

    rel_tol = 0.09
    if args.quick:
        t_final = 0.95
        dt = 0.06
        renorm = 4
        jac_eps = 1.2e-4
        n_quad = 10
    else:
        t_final = 1.35
        dt = 0.05
        renorm = 5
        jac_eps = 1e-4
        n_quad = 14
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

    rng = np.random.default_rng(1630)

    def jac_split(p: np.ndarray, c: np.ndarray) -> np.ndarray:
        return acceleration_jacobian_yukawa_feynman_split(
            p,
            c,
            beta=beta,
            gamma=gamma,
            softening=soft,
            n_quad_feynman=n_quad,
            jac_eps=jac_eps,
        )

    def jac_an(p: np.ndarray, c: np.ndarray) -> np.ndarray:
        return acceleration_jacobian_yukawa_feynman_analytic(
            p,
            c,
            beta=beta,
            gamma=gamma,
            softening=soft,
            n_quad_feynman=n_quad,
        )

    lam_sp, s1 = largest_lyapunov_exponent_benettin_leapfrog(
        pos0,
        vel0,
        C,
        acc_y,
        t_final=t_final,
        dt=dt,
        renorm_every=renorm,
        jac_eps=jac_eps,
        position_jacobian_fn=jac_split,
        rng=rng,
    )
    lam_an, s2 = largest_lyapunov_exponent_benettin_leapfrog(
        pos0,
        vel0,
        C,
        acc_y,
        t_final=t_final,
        dt=dt,
        renorm_every=renorm,
        jac_eps=jac_eps,
        position_jacobian_fn=jac_an,
        rng=rng,
    )

    dlam = abs(lam_sp - lam_an)
    lam_scale = max(abs(lam_sp), abs(lam_an), 1e-9)
    rel = dlam / lam_scale
    print(
        "ex163: yukawa_feynman lambda_max — split J (FD V3) vs analytic J "
        f"(n_quad={n_quad})"
    )
    print(f"  t_final={t_final} dt={dt} steps_sp={s1} steps_an={s2}")
    print(f"  lambda_max (split)   = {lam_sp:.6f}")
    print(f"  lambda_max (analytic)= {lam_an:.6f}")
    print(f"  |delta| = {dlam:.6f}  rel={rel:.4f}  (rel_tol {rel_tol})")

    ok = (
        np.isfinite(lam_sp)
        and np.isfinite(lam_an)
        and min(lam_sp, lam_an) > 0.02
        and rel < rel_tol
    )
    if not ok:
        raise SystemExit(1)
    print("  PASS (finite, chaotic-scale, rel < rel_tol)")


if __name__ == "__main__":
    main()
