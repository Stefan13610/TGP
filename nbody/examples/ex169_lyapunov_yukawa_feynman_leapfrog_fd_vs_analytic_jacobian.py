#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex169_lyapunov_yukawa_feynman_leapfrog_fd_vs_analytic_jacobian.py
=================================================================
P1 (cd.): **ta sama** trajektoria leapfroga, **ten sam** seed styczny —
``lambda_max`` z Benettina przy Jacobianie z **pelnych roznic skonczonych**
vs **``acceleration_jacobian_yukawa_feynman_analytic``**.

Uzupelnia ``ex162``--``ex163`` (tam RK4 albo split/FD na V3): tutaj oba biegi
to ``largest_lyapunov_exponent_benettin_leapfrog`` + ``yukawa_feynman``.
Pelne spektrum ``pairwise``: ``ex173``.

Uruchomienie:
``python ex169_lyapunov_yukawa_feynman_leapfrog_fd_vs_analytic_jacobian.py [--quick]``
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
    largest_lyapunov_exponent_benettin_leapfrog,
    pythagorean_three_body_burrau,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true")
    args = parser.parse_args()

    beta = gamma = 0.08
    soft = 1e-6
    seed = 1691
    pos0, vel0, C = pythagorean_three_body_burrau()

    if args.quick:
        t_final = 0.95
        dt = 0.06
        renorm = 4
        jac_eps = 1.2e-4
        n_quad = 10
        rel_tol = 0.008
    else:
        t_final = 1.15
        dt = 0.052
        renorm = 5
        jac_eps = 9e-5
        n_quad = 12
        rel_tol = 0.012

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

    lam_fd, s_fd = largest_lyapunov_exponent_benettin_leapfrog(
        pos0,
        vel0,
        C,
        acc_y,
        t_final=t_final,
        dt=dt,
        renorm_every=renorm,
        jac_eps=jac_eps,
        position_jacobian_fn=None,
        rng=np.random.default_rng(seed),
    )
    lam_an, s_an = largest_lyapunov_exponent_benettin_leapfrog(
        pos0,
        vel0,
        C,
        acc_y,
        t_final=t_final,
        dt=dt,
        renorm_every=renorm,
        jac_eps=jac_eps,
        position_jacobian_fn=jac_fn,
        rng=np.random.default_rng(seed),
    )

    mxl = max(abs(lam_fd), abs(lam_an), 1e-12)
    rel = abs(lam_fd - lam_an) / mxl

    print(
        "ex169: yukawa_feynman leapfrog+tangent - full FD J vs analytic J "
        f"(same seed={seed})"
    )
    print(
        f"  t_final={t_final} dt={dt} renorm={renorm} jac_eps={jac_eps} "
        f"n_quad={n_quad}"
    )
    print(f"  lambda_max (FD)     = {lam_fd:.8f}  steps={s_fd}")
    print(f"  lambda_max (analytic)= {lam_an:.8f}  steps={s_an}")
    print(f"  rel|delta|/max|lambda| = {rel:.6g}  (tol {rel_tol})")

    ok = (
        np.isfinite(lam_fd)
        and np.isfinite(lam_an)
        and lam_fd > 0.02
        and lam_an > 0.02
        and rel < rel_tol
        and s_fd == s_an
    )
    if not ok:
        raise SystemExit(1)
    print("  PASS")


if __name__ == "__main__":
    main()
