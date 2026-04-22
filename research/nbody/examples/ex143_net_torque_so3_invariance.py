#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex143_net_torque_so3_invariance.py
==================================
Potencjał $V(\\{\\mathbf{x}_i\\})$ zależy tylko od $d_{ij}=|\\mathbf{x}_i-\\mathbf{x}_j|$,
więc jest **invarianty względem globalnego $SO(3)$** (obrót wszystkich współrzędnych).
Dla izolowanego układu:

    d/dt sum_i C_i x_i x v_i = sum_i x_i x F_i  ->  0

Sprawdzamy numerycznie ``||sum_i x_i x (C_i a_i)||`` względem skali
``max_i |x_i| * max_i |F_i|`` dla wszystkich ``TGP_INTEGRATION_BACKENDS``.

Zob. ``dynamics_backends.total_torque_about_origin_from_accelerations`` i
[tgp_nbody_lagrangian_eom.tex](../tgp_nbody_lagrangian_eom.tex).
"""

from __future__ import annotations

import os
import sys

import numpy as np

_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

from nbody.dynamics_backends import (
    TGP_INTEGRATION_BACKENDS,
    build_tgp_integration_pair,
    total_torque_about_origin_from_accelerations,
)


def main() -> None:
    rng = np.random.default_rng(143)
    beta = gamma = 1.0
    eps = 1e-6
    n_quad = 20
    n_cases = 12
    worst = 0.0
    worst_meta = ""

    print("ex143: net torque |sum_i x_i x (C_i a_i)| wzgledem skali r*F")
    for backend in TGP_INTEGRATION_BACKENDS:
        acc_fn, _ = build_tgp_integration_pair(
            backend,
            beta=beta,
            gamma=gamma,
            softening=eps,
            include_3body=True,
            n_quad_feynman=n_quad,
        )
        local_worst_ratio = 0.0
        for _ in range(n_cases):
            n = int(rng.integers(2, 7))
            pos = rng.normal(scale=1.05, size=(n, 3))
            C = np.abs(rng.normal(size=n)) + 0.15
            a = acc_fn(pos, C)
            tau = total_torque_about_origin_from_accelerations(pos, a, C)
            tnorm = float(np.linalg.norm(tau))
            F = C[:, None] * a
            fmag = float(np.max(np.linalg.norm(F, axis=1)))
            rmax = float(np.max(np.linalg.norm(pos, axis=1)))
            scale = rmax * fmag
            bound = 1e-8 * max(1.0, scale) + 1e-12
            ratio = tnorm / bound
            local_worst_ratio = max(local_worst_ratio, ratio)
            if tnorm > worst:
                worst = tnorm
                worst_meta = f"{backend} n={n} |tau|={tnorm:.3e}"
        ok_b = local_worst_ratio <= 1.0
        print(
            f"  {backend:22s} max(|tau|/tol_bound) = {local_worst_ratio:.3g}  {'PASS' if ok_b else 'FAIL'}"
        )
        if not ok_b:
            raise SystemExit(1)

    print(f"  global worst |tau| = {worst:.3e}  ({worst_meta})  PASS")


if __name__ == "__main__":
    main()
