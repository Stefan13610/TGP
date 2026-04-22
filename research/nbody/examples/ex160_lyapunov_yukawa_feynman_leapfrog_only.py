#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex160_lyapunov_yukawa_feynman_leapfrog_only.py
==============================================
P1 (cd.): pojedynczy estymator ``lambda_max`` dla **yukawa_feynman** z
**leapfrog+tangent** i Jacobianem FD (N=3, Burrau).

Uzupełnia ``ex152`` (tam RK4 na rozszerzonym ODE + porownanie z coulomb_3b).

Uruchomienie: ``python ex160_lyapunov_yukawa_feynman_leapfrog_only.py [--quick] [--n-quad N]``
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
from nbody.lyapunov import largest_lyapunov_exponent_benettin_leapfrog, pythagorean_three_body_burrau


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

    lam, steps = largest_lyapunov_exponent_benettin_leapfrog(
        pos0,
        vel0,
        C,
        acc_y,
        t_final=t_final,
        dt=dt,
        renorm_every=renorm,
        jac_eps=jac_eps,
        position_jacobian_fn=None,
        rng=np.random.default_rng(1600),
    )

    print(
        "ex160: lambda_max yukawa_feynman, leapfrog+tangent, FD Jacobian "
        f"(n_quad={n_quad})"
    )
    print(f"  t_final={t_final} dt={dt} steps={steps}")
    print(f"  lambda_max = {lam:.6f}")

    ok = np.isfinite(lam) and lam > 0.02
    if not ok:
        raise SystemExit(1)
    print("  PASS (finite, lambda > 0.02)")


if __name__ == "__main__":
    main()
