#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex141_net_force_translational_invariance.py
============================================
Weryfikacja numeryczna niezmiennika translacyjnego:

    sum_i F_i = sum_i C_i * a_i  ~  0

dla wszystkich backendow z ``TGP_INTEGRATION_BACKENDS`` (potencjał zależy tylko
od względnych odległości).  Używa ``dynamics_backends.total_force_from_accelerations``.

Zob. [EOM_PROGRAM_NBODY.md](../EOM_PROGRAM_NBODY.md) (punkt o redukcji / stałych ruchu).
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
    total_force_from_accelerations,
)


def main() -> None:
    rng_master = np.random.default_rng(2026)
    beta = gamma = 1.0
    eps = 1e-6
    n_quad = 20
    n_cases = 10
    worst = 0.0
    worst_meta = ""

    print("ex141: net force |sum_i C_i a_i| (wszystkie backendy)")
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
            n = int(rng_master.integers(2, 7))
            pos = rng_master.normal(scale=1.1, size=(n, 3))
            C = np.abs(rng_master.normal(size=n)) + 0.15
            a = acc_fn(pos, C)
            net = total_force_from_accelerations(a, C)
            fn = float(np.linalg.norm(net))
            fscale = float(np.max(np.linalg.norm(C[:, None] * a, axis=1)))
            bound = 1e-9 * max(1.0, fscale) + 1e-12
            ratio = fn / bound
            local_worst_ratio = max(local_worst_ratio, ratio)
            if fn > worst:
                worst = fn
                worst_meta = f"{backend} n={n} |net|={fn:.3e}"
        ok_b = local_worst_ratio <= 1.0
        print(
            f"  {backend:22s} max(|net|/tol_bound) = {local_worst_ratio:.3g}  {'PASS' if ok_b else 'FAIL'}"
        )
        if not ok_b:
            raise SystemExit(1)

    print(f"  global worst |net| = {worst:.3e}  ({worst_meta})  PASS")


if __name__ == "__main__":
    main()
