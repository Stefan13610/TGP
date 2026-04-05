#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex157_lyapunov_tgp_pairwise_rk4_vs_leapfrog.py
==============================================
P1 (cd.): ``lambda_max`` dla **TGP pairwise** (Burrau), Jacobian ``a`` wzgledem
``x`` przez **roznicowe centralne** — porownanie **RK4+tangent** vs
**leapfrog+tangent**.

Ten sam ``acc_fn``, rozne mapy styczne (jak ``ex156`` dla Newtona).

Uruchomienie: ``python ex157_lyapunov_tgp_pairwise_rk4_vs_leapfrog.py [--quick]``
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np

_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

from nbody.dynamics_v2 import analytical_forces_tgp_pairwise
from nbody.lyapunov import (
    largest_lyapunov_exponent_benettin,
    largest_lyapunov_exponent_benettin_leapfrog,
    pythagorean_three_body_burrau,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true")
    args = parser.parse_args()

    beta = gamma = 0.08
    soft = 1e-6
    pos0, vel0, C = pythagorean_three_body_burrau()

    def acc_tgp(p: np.ndarray, c: np.ndarray) -> np.ndarray:
        return analytical_forces_tgp_pairwise(p, c, beta, gamma, soft) / c[:, None]

    if args.quick:
        t_final = 8.0
        dt = 0.012
        renorm = 14
        jac_eps = 1e-4
        delta_tol = 0.75
    else:
        t_final = 14.0
        dt = 0.008
        renorm = 18
        jac_eps = 8e-5
        delta_tol = 0.65

    lam_rk, sr = largest_lyapunov_exponent_benettin(
        pos0,
        vel0,
        C,
        acc_tgp,
        t_final=t_final,
        dt=dt,
        renorm_every=renorm,
        jac_eps=jac_eps,
        position_jacobian_fn=None,
        rng=np.random.default_rng(1571),
    )
    lam_lf, sl = largest_lyapunov_exponent_benettin_leapfrog(
        pos0,
        vel0,
        C,
        acc_tgp,
        t_final=t_final,
        dt=dt,
        renorm_every=renorm,
        jac_eps=jac_eps,
        position_jacobian_fn=None,
        rng=np.random.default_rng(1572),
    )

    print(
        "ex157: lambda_max TGP pairwise (Burrau), FD Jacobian, "
        "RK4 tangent vs leapfrog tangent"
    )
    print(f"  beta=gamma={beta} soft={soft} t_final={t_final} dt={dt}")
    print(f"  steps_rk={sr} steps_lf={sl}")
    print(f"  lambda_max RK4+tangent      = {lam_rk:.6f}")
    print(f"  lambda_max leapfrog+tangent  = {lam_lf:.6f}")
    print(f"  |delta| = {abs(lam_rk - lam_lf):.6f}")

    ok = (
        np.isfinite(lam_rk)
        and np.isfinite(lam_lf)
        and lam_rk > 0.02
        and lam_lf > 0.02
        and abs(lam_rk - lam_lf) < delta_tol
    )
    if not ok:
        raise SystemExit(1)
    print(f"  PASS (finite, chaotic-scale, |delta| < {delta_tol})")


if __name__ == "__main__":
    main()
