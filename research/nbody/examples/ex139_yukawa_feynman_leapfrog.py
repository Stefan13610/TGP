#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex139_yukawa_feynman_leapfrog.py
================================
Krótki odcinek dynamiki z backendem **``yukawa_feynman``**:
$V_2$ (z softeningiem) + dokładne $I_Y$ z kwadratury Feynmana 2D
(`three_body_force_exact`).

Kontynuacja programu z [EOM_PROGRAM_NBODY.md](../EOM_PROGRAM_NBODY.md) i
[../dynamics_backends.py](../dynamics_backends.py).

Uwagi:
  - Koszt ~ $O(N^3 \\cdot n_{\\rm quad}^2)$ na ewaluację przyspieszenia; tu $N=3$
    (jeden triplet), `n_quad_feynman` obniżone do 22 dla rozsądnego czasu.
  - Sektor 3B używa twardych odległości euklidesowych (bez softeningu), sektor 2B —
    ze `softening` jak w `dynamics_v2.potential_tgp` (patrz dokumentacja backendu).
  - Krótki horyzont + mały krok: leapfrog ma być stabilny przy umiarkowanym $m_{\\rm sp}d$.

Test: względna zmiana energii całkowitej < 5e-5 na zapisanych punktach czasu.
"""

from __future__ import annotations

import os
import sys
import time

import numpy as np

_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

from nbody.dynamics_backends import build_tgp_integration_pair
from nbody import dynamics_v2


def main() -> None:
    beta = gamma = 1.0
    eps = 1e-5
    n_quad = 22

    a = 3.0
    pos0 = np.array(
        [
            [0.0, 0.0, 0.0],
            [a, 0.02, 0.0],
            [0.5 * a, 0.866 * a, 0.01],
        ],
        dtype=float,
    )
    C = np.array([0.25, 0.25, 0.25], dtype=float)
    vel0 = np.zeros_like(pos0)

    acc_fn, pot_fn = build_tgp_integration_pair(
        "yukawa_feynman",
        beta=beta,
        gamma=gamma,
        softening=eps,
        include_3body=True,
        n_quad_feynman=n_quad,
    )

    t0 = time.perf_counter()
    out = dynamics_v2.leapfrog_integrate(
        pos0,
        vel0,
        C,
        acc_fn,
        pot_fn,
        t_span=(0.0, 0.12),
        dt=8e-5,
        save_every=20,
    )
    elapsed = time.perf_counter() - t0

    E = out["energy"]
    dE = float(np.max(E) - np.min(E))
    rel = dE / (abs(np.mean(E)) + 1e-30)
    ok = rel < 5e-5

    print("ex139: leapfrog + build_tgp_integration_pair('yukawa_feynman')")
    print(f"  n_quad={n_quad}  wall_time={elapsed:.2f}s")
    print(f"  mean(E)={np.mean(E):.6g}  max|dE|/mean|E|={rel:.3g}  {'PASS' if ok else 'FAIL'}")
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
