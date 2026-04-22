#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex147_com_acceleration_zero.py
==============================
Z ``tgp_nbody_lagrangian_eom.tex``: przy wewnetrznym V od d_ij,

    a_cm = (sum_i C_i a_i) / M_tot = (sum_i F_i) / M_tot = 0.

Sprawdzamy ``||center_of_mass_acceleration(a)||`` dla backendow
``coulomb_3b`` i ``yukawa_feynman`` (losowe konfiguracje).
"""

from __future__ import annotations

import os
import sys

import numpy as np

_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

from nbody.dynamics_backends import (
    build_tgp_integration_pair,
    center_of_mass_acceleration,
)


def main() -> None:
    rng = np.random.default_rng(147)
    beta = gamma = 1.0
    eps = 1e-6
    n_quad = 20
    worst = 0.0

    print("ex147: |a_cm| = |sum C a| / M_tot  (coulomb_3b, yukawa_feynman)")
    for backend in ("coulomb_3b", "yukawa_feynman"):
        acc_fn, _ = build_tgp_integration_pair(
            backend,
            beta=beta,
            gamma=gamma,
            softening=eps,
            include_3body=True,
            n_quad_feynman=n_quad,
        )
        local = 0.0
        for _ in range(16):
            n = int(rng.integers(2, 7))
            pos = rng.normal(scale=1.05, size=(n, 3))
            C = np.abs(rng.normal(size=n)) + 0.14
            a = acc_fn(pos, C)
            acm = center_of_mass_acceleration(a, C)
            local = max(local, float(np.linalg.norm(acm)))
        worst = max(worst, local)
        ok_b = local < 1e-8
        print(f"  {backend:16s} max|a_cm|={local:.3e}  {'PASS' if ok_b else 'FAIL'}")
        if not ok_b:
            raise SystemExit(1)
    print(f"  global worst = {worst:.3e}  PASS")


if __name__ == "__main__":
    main()
