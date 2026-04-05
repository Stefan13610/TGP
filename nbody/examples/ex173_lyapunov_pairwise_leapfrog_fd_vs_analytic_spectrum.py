#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex173_lyapunov_pairwise_leapfrog_fd_vs_analytic_spectrum.py
===========================================================
P1 (cd.): **pairwise** (``V_2``) + leapfrog+tangent + **pelne** spektrum
``k=6N`` — porownanie Jacobianu z **pelnych roznic skonczonych** vs
``acceleration_jacobian_tgp_pairwise_softened``. Ten sam seed, ta sama trajektoria.

Uzupelnia ``ex169`` (Yukawa, pojedynczy ``lambda_max``) oraz ``ex172`` (sama suma
spektrum przy analitycznym ``J``).

Uruchomienie:
``python ex173_lyapunov_pairwise_leapfrog_fd_vs_analytic_spectrum.py [--quick]``
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
    acceleration_jacobian_tgp_pairwise_softened,
    lyapunov_spectrum_benettin_leapfrog,
    pythagorean_three_body_burrau,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true")
    args = parser.parse_args()

    beta = gamma = 0.08
    soft = 1e-6
    seed = 1731
    pos0, vel0, C = pythagorean_three_body_burrau()
    n = len(C)
    k = 6 * n

    if args.quick:
        t_final = 0.95
        dt = 0.06
        renorm = 4
        jac_eps = 1.2e-4
        rel_tol = 5e-5
    else:
        t_final = 1.15
        dt = 0.052
        renorm = 5
        jac_eps = 9e-5
        rel_tol = 2e-4

    acc_p, _ = build_tgp_integration_pair(
        "pairwise",
        beta=beta,
        gamma=gamma,
        softening=soft,
    )

    def jac_fn(p: np.ndarray, c: np.ndarray) -> np.ndarray:
        return acceleration_jacobian_tgp_pairwise_softened(
            p, c, beta=beta, gamma=gamma, softening=soft
        )

    spec_fd, s_fd = lyapunov_spectrum_benettin_leapfrog(
        pos0,
        vel0,
        C,
        acc_p,
        n_exponents=k,
        t_final=t_final,
        dt=dt,
        renorm_every=renorm,
        jac_eps=jac_eps,
        position_jacobian_fn=None,
        rng=np.random.default_rng(seed),
    )
    spec_an, s_an = lyapunov_spectrum_benettin_leapfrog(
        pos0,
        vel0,
        C,
        acc_p,
        n_exponents=k,
        t_final=t_final,
        dt=dt,
        renorm_every=renorm,
        jac_eps=jac_eps,
        position_jacobian_fn=jac_fn,
        rng=np.random.default_rng(seed),
    )

    mxa = float(np.max(np.abs(spec_an)))
    rel_max = float(np.max(np.abs(spec_fd - spec_an))) / max(mxa, 1e-12)
    sum_fd = float(np.sum(spec_fd))
    sum_an = float(np.sum(spec_an))

    print(
        "ex173: pairwise k=6N spectrum - FD J vs analytic J "
        f"(same seed={seed})"
    )
    print(
        f"  t_final={t_final} dt={dt} renorm={renorm} jac_eps={jac_eps} k={k}"
    )
    print(f"  max|lambda_fd - lambda_an|/max|lambda_an| = {rel_max:.6e}  (tol {rel_tol})")
    print(f"  sum(lambda) FD={sum_fd:.6e}  analytic={sum_an:.6e}")
    print(f"  lambda[0] FD={spec_fd[0]:.8f}  analytic={spec_an[0]:.8f}  steps={s_fd}/{s_an}")

    ok = (
        np.all(np.isfinite(spec_fd))
        and np.all(np.isfinite(spec_an))
        and rel_max < rel_tol
        and float(spec_fd[0]) > 0.02
        and float(spec_an[0]) > 0.02
        and s_fd == s_an
    )
    if not ok:
        raise SystemExit(1)
    print("  PASS")


if __name__ == "__main__":
    main()
