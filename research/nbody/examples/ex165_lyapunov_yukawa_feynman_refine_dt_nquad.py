#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex165_lyapunov_yukawa_feynman_refine_dt_nquad.py
================================================
P1 (cd.): **zbieżność** estymatu ``lambda_max`` dla ``yukawa_feynman`` przy
**analitycznym** ``J`` - dwa przebiegi z **tym samym** seedem stycznym
(``rng(1650)``), tym samym ``t_final``, lecz różnym ``dt`` i ``n_quad``
(dokładniejszy = mniejszy ``dt``, wyższe ``n_quad``).

Nie dowodzi asymptotyki ``dt->0`` (krotki horyzont); wskazuje, ze po
dopieszczeniu dyskretyzacji ``lambda_max`` nie zmienia sie skokowo o rzad
wielkosci - regresja przeciw przypadkowym regresjom zbyt chropowatej siatki.

Uruchomienie: ``python ex165_lyapunov_yukawa_feynman_refine_dt_nquad.py [--quick]``
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


def _run_lam(
    *,
    t_final: float,
    dt: float,
    n_quad: int,
    renorm: int,
    seed: int,
) -> float:
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

    lam, _st = largest_lyapunov_exponent_benettin_leapfrog(
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
    return float(lam)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true")
    args = parser.parse_args()

    seed = 1650
    if args.quick:
        t_final = 4.0
        renorm = 10
        dt_coarse, nq_coarse = 0.05, 14
        dt_fine, nq_fine = 0.032, 20
        rel_tol = 0.13
    else:
        t_final = 6.0
        renorm = 12
        dt_coarse, nq_coarse = 0.04, 20
        dt_fine, nq_fine = 0.036, 22
        rel_tol = 0.06

    lam_c = _run_lam(
        t_final=t_final,
        dt=dt_coarse,
        n_quad=nq_coarse,
        renorm=renorm,
        seed=seed,
    )
    lam_f = _run_lam(
        t_final=t_final,
        dt=dt_fine,
        n_quad=nq_fine,
        renorm=renorm,
        seed=seed,
    )

    scale = max(abs(lam_c), abs(lam_f), 1e-9)
    rel = abs(lam_c - lam_f) / scale

    print("ex165: yukawa_feynman lambda_max - coarse vs fine (analytic J, same seed)")
    print(f"  t_final={t_final} renorm_every={renorm} seed={seed}")
    print(
        f"  coarse: dt={dt_coarse} n_quad={nq_coarse}  lambda={lam_c:.6f}"
    )
    print(f"  fine:   dt={dt_fine} n_quad={nq_fine}  lambda={lam_f:.6f}")
    print(f"  rel|delta|/max|lambda| = {rel:.4f}  (tol {rel_tol})")

    ok = (
        np.isfinite(lam_c)
        and np.isfinite(lam_f)
        and min(lam_c, lam_f) > 0.02
        and rel < rel_tol
    )
    if not ok:
        raise SystemExit(1)
    print("  PASS (finite, lambda>0.02, relative refinement OK)")


if __name__ == "__main__":
    main()
