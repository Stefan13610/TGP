#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex144_conserved_P_L_leapfrog.py
===============================
Przy izolowanym układzie i potencjale tylko od $d_{ij}$ oczekujemy zachowania:

  P = sum_i C_i v_i,   L = sum_i C_i x_i x v_i

(w układzie bezwładnym, w którym mierzymy wspolrzedne).  Start z obrotem
sztywnym wokol srodka mas (w plaszczyznie xy), wiec P(0)=0.

Integracja: ``coulomb_3b`` + leapfrog (jak ``ex138``), potem sprawdzenie dryfu
|P| i |L| na zapisanych krokach.

PASS: wzgledny dryf max |P| i max |L| ponizej progow.
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
    total_angular_momentum,
    total_linear_momentum,
)
from nbody import dynamics_v2


def main() -> None:
    beta = gamma = 1.0
    eps = 1e-5
    acc_fn, pot_fn = build_tgp_integration_pair(
        "coulomb_3b",
        beta=beta,
        gamma=gamma,
        softening=eps,
        include_3body=True,
    )

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
    cm = np.sum(C[:, None] * pos0, axis=0) / np.sum(C)
    omega = 0.1
    vel0 = np.zeros_like(pos0)
    for i in range(3):
        r = pos0[i] - cm
        vel0[i] = omega * np.array([-r[1], r[0], 0.0])

    out = dynamics_v2.leapfrog_integrate(
        pos0,
        vel0,
        C,
        acc_fn,
        pot_fn,
        t_span=(0.0, 0.2),
        dt=8e-5,
        save_every=20,
    )
    P0 = total_linear_momentum(out["velocities"][0], C)
    L0 = total_angular_momentum(out["positions"][0], out["velocities"][0], C)
    pnorm0 = float(np.linalg.norm(P0))
    lnorm0 = float(np.linalg.norm(L0))

    p_drift = 0.0
    l_drift = 0.0
    for k in range(len(out["t"])):
        P = total_linear_momentum(out["velocities"][k], C)
        L = total_angular_momentum(out["positions"][k], out["velocities"][k], C)
        p_drift = max(p_drift, float(np.linalg.norm(P - P0)))
        l_drift = max(l_drift, float(np.linalg.norm(L - L0)))

    tol_p = 1e-9 * max(1.0, pnorm0) + 1e-12
    tol_l = 1e-7 * max(1.0, lnorm0) + 1e-14
    ok_p = p_drift < max(tol_p, 1e-8)
    ok_l = l_drift < max(tol_l, 1e-10)

    print("ex144: leapfrog coulomb_3b - dryf P i L")
    print(f"  |P0|={pnorm0:.3e}  max|dP|={p_drift:.3e}  {'PASS' if ok_p else 'FAIL'}")
    print(f"  |L0|={lnorm0:.3e}  max|dL|={l_drift:.3e}  {'PASS' if ok_l else 'FAIL'}")
    if not (ok_p and ok_l):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
