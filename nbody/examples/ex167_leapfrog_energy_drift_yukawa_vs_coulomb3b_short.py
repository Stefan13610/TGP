#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex167_leapfrog_energy_drift_yukawa_vs_coulomb3b_short.py
========================================================
P1 (cd.): **kontrola energii** przy ``leapfrog_integrate`` dla IC Burrau (v=0):
``yukawa_feynman`` vs ``coulomb_3b`` (oba z sektorem V3).

UWAGA: dla tej konfiguracji przy **dluzszym** ``t`` energia moze **eksplodowac**
(bliskie przejscia, miekki regulator). Regresja celowo uzywa **bardzo krotkiego**
okna czasu; nie jest to dowod stabilnosci energetycznej na horyzoncie fizycznym.

Na Burrau z dopasowana energia (cztery galezie jak ``ex178``): ``ex190_leapfrog_energy_drift_matched_h_family_short.py`` (leapfrog), ``ex191_rk45_energy_diag_matched_h_family_short.py`` (RK45), ``ex192_matched_h_leapfrog_vs_rk45_energy_table.py`` (tabela); ``ex193``/``ex194`` (CSV LF/RK + odczyt).

Uruchomienie: ``python ex167_leapfrog_energy_drift_yukawa_vs_coulomb3b_short.py [--quick]``
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
from nbody.dynamics_v2 import leapfrog_integrate
from nbody.lyapunov import pythagorean_three_body_burrau


def _max_rel_energy_drift(
    *,
    backend: str,
    t_final: float,
    dt: float,
    save_every: int,
    n_quad: int,
    beta: float,
    gamma: float,
    softening: float,
) -> float:
    pos0, vel0, C = pythagorean_three_body_burrau()
    acc, pot = build_tgp_integration_pair(
        backend,
        beta=beta,
        gamma=gamma,
        softening=softening,
        include_3body=True,
        n_quad_feynman=n_quad,
    )
    r = leapfrog_integrate(
        pos0,
        vel0,
        C,
        acc,
        pot,
        t_span=(0.0, t_final),
        dt=dt,
        save_every=save_every,
    )
    return float(np.max(np.abs(r["energy_error"])))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true")
    args = parser.parse_args()

    beta = gamma = 0.08
    softening = 1e-6

    if args.quick:
        t_final = 0.25
        dt = 0.02
        save_every = 2
        n_quad = 16
        tol = 0.012
    else:
        t_final = 0.28
        dt = 0.018
        save_every = 3
        n_quad = 18
        tol = 0.02

    backends = ("yukawa_feynman", "coulomb_3b")
    drifts = {}
    for b in backends:
        drifts[b] = _max_rel_energy_drift(
            backend=b,
            t_final=t_final,
            dt=dt,
            save_every=save_every,
            n_quad=n_quad,
            beta=beta,
            gamma=gamma,
            softening=softening,
        )

    print(
        "ex167: leapfrog max |(E-E0)/E0| (Burrau v=0), "
        f"t_final={t_final}, dt={dt}, n_quad={n_quad}"
    )
    for b in backends:
        ok = drifts[b] < tol
        flag = "PASS" if ok else "FAIL"
        print(f"  {b}: max|dE/E0|={drifts[b]:.6g}  (tol {tol})  {flag}")
        if not ok:
            raise SystemExit(1)
    print("ex167: all PASS")


if __name__ == "__main__":
    main()
