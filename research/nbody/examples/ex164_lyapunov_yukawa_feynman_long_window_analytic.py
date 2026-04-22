#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex164_lyapunov_yukawa_feynman_long_window_analytic.py
=====================================================
P1 (cd.): **yukawa_feynman** + leapfrog+tangent z **pełnym analitycznym**
``acceleration_jacobian_yukawa_feynman_analytic`` — **dłuższe okno** czasu i
wyższe ``n_quad`` niż ``ex160`` (krótki bieg z FD).

Cel: stabilniejszy estymator ``lambda_max`` przy rozsądnym koszcie (brak FD na
``J``). Tryb pełny (bez ``--quick``) jeszcze wydłuża integrację — do ręcznych
uruchomień „publication-ready'', nie do najszybszej regresji.

Po Benettinie (leapfrog+styczna): **diagnostyka energii** przez ``rk45_integrate``
na ``[0, t_final]`` dla ``yukawa_feynman`` — tylko informacja (orbita Benettina
to leapfrog, nie RK45); ``solve_ivp`` moze zakonczyc sie ``success=False`` przy
Burrau mimo malego ``max|dE/E0|`` na ``t_eval``.

Uruchomienie: ``python ex164_lyapunov_yukawa_feynman_long_window_analytic.py [--quick] [--n-quad N]``
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
from nbody.dynamics_v2 import rk45_integrate
from nbody.lyapunov import (
    acceleration_jacobian_yukawa_feynman_analytic,
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

    if args.quick:
        t_final = 4.25
        dt = 0.04
        renorm = 10
        n_quad = 16
    else:
        t_final = 8.0
        dt = 0.03
        renorm = 12
        n_quad = 22
    if args.n_quad > 0:
        n_quad = int(args.n_quad)

    acc_y, pot_y = build_tgp_integration_pair(
        "yukawa_feynman",
        beta=beta,
        gamma=gamma,
        softening=soft,
        include_3body=True,
        n_quad_feynman=n_quad,
    )

    def jac_fn(p: np.ndarray, c: np.ndarray) -> np.ndarray:
        return acceleration_jacobian_yukawa_feynman_analytic(
            p,
            c,
            beta=beta,
            gamma=gamma,
            softening=soft,
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
        jac_eps=1e-5,
        position_jacobian_fn=jac_fn,
        rng=np.random.default_rng(1640),
    )

    print(
        "ex164: lambda_max yukawa_feynman, long window, leapfrog+tangent, "
        f"analytic J (n_quad={n_quad})"
    )
    print(f"  t_final={t_final} dt={dt} renorm_every={renorm} steps={steps}")
    print(f"  lambda_max = {lam:.6f}")

    r_e = rk45_integrate(
        pos0,
        vel0,
        C,
        acc_y,
        pot_y,
        t_span=(0.0, float(t_final)),
        n_output=500,
        rtol=1e-9,
        atol=1e-11,
        quiet=True,
    )
    mx_e = float(np.max(np.abs(r_e["energy_error"])))
    print(
        f"  rk45 energy diag yukawa_feynman: success={r_e['success']}  "
        f"max|(E-E0)/E0|={mx_e:.3e}  (reference Hamiltonian flow, not leapfrog orbit)"
    )

    ok = np.isfinite(lam) and lam > 0.02
    if not ok:
        raise SystemExit(1)
    print("  PASS (finite, lambda > 0.02)")


if __name__ == "__main__":
    main()
