#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex146_yukawa_feynman_grad_matches_force.py
==========================================
Jak ``ex145``, ale backend ``yukawa_feynman`` (dokladne I_Y + V_2 ze softeningiem).

Ograniczenie do N=3 (jeden triplet) i umiarkowane ``n_quad_feynman`` ze wzgledu na koszt
centralnych roznic (6 ewaluacji V na skladowa i na cialo).

PASS: max |F_fd - F| / max |F| < 1e-4.
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
    rng = np.random.default_rng(146)
    beta = gamma = 1.0
    soft = 1e-6
    fd_eps = 3e-5
    n_quad = 18

    acc_fn, pot_fn = build_tgp_integration_pair(
        "yukawa_feynman",
        beta=beta,
        gamma=gamma,
        softening=soft,
        include_3body=True,
        n_quad_feynman=n_quad,
    )

    worst = 0.0
    for _ in range(5):
        pos = rng.normal(scale=1.0, size=(3, 3))
        C = np.array([0.25, 0.28, 0.22], dtype=float)
        F_code = C[:, None] * acc_fn(pos, C)
        F_fd = forces_from_potential_central_diff(pot_fn, pos, C, eps=fd_eps)
        denom = max(1e-30, float(np.max(np.linalg.norm(F_code, axis=1))))
        rel = float(np.max(np.linalg.norm(F_fd - F_code, axis=1)) / denom)
        worst = max(worst, rel)

    ok = worst < 1e-4
    print("ex146: -grad V (central diff) vs C*a  (yukawa_feynman, N=3)")
    print(f"  n_quad={n_quad}  max relative |dF|/|F| = {worst:.3e}  {'PASS' if ok else 'FAIL'}")
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
