#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex138_eom_tgp_coulomb_leapfrog.py
=================================
Demonstracja **dwóch warstw** z [EOM_PROGRAM_NBODY.md](../EOM_PROGRAM_NBODY.md):

  1. Analityka + wybór backendu: ``dynamics_backends.build_tgp_integration_pair('coulomb_3b')``.
  2. Czas: ``dynamics_v2.leapfrog_integrate`` z zwróconą parą ``(acc_fn, pot_fn)``.

Nie jest to nowy wynik fizyczny — tylko kanoniczny szkielet pod dalsze backendy
(Yukawa-exact → ``three_body_force_exact``).

Test: krótki horyzont czasu (układ bez pędu początkowego szybko się zwija — sztywny
potencjał; długi leapfrog z stałym krokiem traci energię przy zbliżeniach).
"""

from __future__ import annotations

import os
import sys

import numpy as np

# Katalog roboczy = TGP_v1/ — import pakietu ``nbody`` (wymagane przez łańcuch relative importów).
_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

from nbody.dynamics_backends import build_tgp_integration_pair
from nbody import dynamics_v2


def main() -> None:
    beta = gamma = 1.0
    eps = 1e-5
    # Trójkąt równoboczny w płaszczyźnie z = 0, lekko zaburzony
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
        "coulomb_3b",
        beta=beta,
        gamma=gamma,
        softening=eps,
        include_3body=True,
    )

    # Krótki horyzont: pełny zwój wymaga adaptacyjnego kroku lub RK — tu tylko spójność EOM+potencjału.
    out = dynamics_v2.leapfrog_integrate(
        pos0,
        vel0,
        C,
        acc_fn,
        pot_fn,
        t_span=(0.0, 0.25),
        dt=1e-4,
        save_every=25,
    )
    E = out["energy"]
    dE = float(np.max(E) - np.min(E))
    rel = dE / (abs(np.mean(E)) + 1e-30)
    ok = rel < 1e-6
    print("ex138: leapfrog + build_tgp_integration_pair('coulomb_3b')")
    print(f"  mean(E)={np.mean(E):.6g}  max|dE|/mean|E|={rel:.3g}  {'PASS' if ok else 'FAIL'}")
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
