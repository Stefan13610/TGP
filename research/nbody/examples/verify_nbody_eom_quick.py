#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify_nbody_eom_quick.py
=========================
Lekka regresja skryptow EOM / ``dynamics_backends`` (bez ciezkich ``ex139``/``ex140``).

Uruchom z dowolnego katalogu:

    python nbody/examples/verify_nbody_eom_quick.py

(wymaga ``TGP_v1/`` jako cwd lub poprawnej sciezki do tego pliku).

Zakres **tylko** ``nbody/``; ewentualne rozjazdy z glownym tomem LaTeX scalasz osobno.
"""

from __future__ import annotations

import os
import subprocess
import sys

_SCRIPTS = (
    "ex138_eom_tgp_coulomb_leapfrog.py",
    "ex141_net_force_translational_invariance.py",
    "ex142_yukawa_overlap_quadrature_convergence.py",
    "ex143_net_torque_so3_invariance.py",
    "ex144_conserved_P_L_leapfrog.py",
    "ex145_minus_grad_V_matches_force.py",
    "ex146_yukawa_feynman_grad_matches_force.py",
    "ex147_com_acceleration_zero.py",
)


def main() -> None:
    here = os.path.dirname(os.path.abspath(__file__))
    repo = os.path.abspath(os.path.join(here, "..", ".."))
    py = sys.executable
    failed = []
    for name in _SCRIPTS:
        path = os.path.join(here, name)
        print(f"==> {name}", flush=True)
        r = subprocess.run([py, path], cwd=repo, env={**os.environ, "PYTHONUNBUFFERED": "1"})
        if r.returncode != 0:
            failed.append(name)
    if failed:
        print("FAILED:", ", ".join(failed))
        raise SystemExit(1)
    print("verify_nbody_eom_quick: all PASS")


if __name__ == "__main__":
    main()
