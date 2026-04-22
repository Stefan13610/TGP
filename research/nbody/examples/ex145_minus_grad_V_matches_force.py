#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex145_minus_grad_V_matches_force.py
===================================
Walidacja spojnosci: dla potencjalu ``pot_fn`` z backendu ``coulomb_3b``,
sily z centralnych roznic ``-grad V`` (``forces_from_potential_central_diff``)
powinny zgadzac sie z ``C_i * a_i`` z ``acc_fn``.

To jest **warstwa numeryczna** (koszt O(N) ewaluacji V), ale domyka petle
``V -> F`` dla zamknietego backendu Coulomba na I.

PASS: max |F_fd - F_code| / max |F_code| < 5e-5.
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
    forces_from_potential_central_diff,
)


def main() -> None:
    rng = np.random.default_rng(145)
    beta = gamma = 1.0
    soft = 1e-6
    fd_eps = 8e-6

    acc_fn, pot_fn = build_tgp_integration_pair(
        "coulomb_3b",
        beta=beta,
        gamma=gamma,
        softening=soft,
        include_3body=True,
    )

    worst = 0.0
    for trial in range(6):
        n = int(rng.integers(3, 6))
        pos = rng.normal(scale=1.0, size=(n, 3))
        C = np.abs(rng.normal(size=n)) + 0.2
        F_code = C[:, None] * acc_fn(pos, C)
        F_fd = forces_from_potential_central_diff(pot_fn, pos, C, eps=fd_eps)
        denom = max(1e-30, float(np.max(np.linalg.norm(F_code, axis=1))))
        rel = float(np.max(np.linalg.norm(F_fd - F_code, axis=1)) / denom)
        worst = max(worst, rel)

    ok = worst < 5e-5
    print("ex145: -grad V (central diff) vs C*a  (coulomb_3b)")
    print(f"  max relative |dF|/|F| = {worst:.3e}  {'PASS' if ok else 'FAIL'}")
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
