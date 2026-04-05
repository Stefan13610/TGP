#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex148_lyapunov_tgp_vs_newton_pairwise.py
========================================
P1 z ``PLAN_ROZWOJU_NBODY.md``: porownanie najwiekszego wykladnika Lapunowa
dla klasycznego problemu Pitagorejskiego (3 ciala, m propr. 3,4,5) miedzy

  - Newton: ``forces_newton`` / C,
  - TGP: wylacznie sektor parami ``analytical_forces_tgp_pairwise`` / C
    (bez V3 — izolowany wplyw modyfikacji 2-cialowej).

Opcjonalnie: skalowanie predkosci w TGP tak, by energia mechaniczna
(przy tych samych pozycjach) rownala sie energii Newtona (to samo x, inne v).

Pelny backend ``yukawa_feynman`` z dopasowaniem ``H`` (losowe ``v`` w galezi Yukawy):
``ex174_lyapunov_newton_vs_yukawa_feynman_matched_energy.py``.
Ta sama para Newton vs pairwise (tylko ``V_2``), leapfrog + jawny ``J``:
``ex177_lyapunov_newton_vs_pairwise_matched_energy.py``.

Uruchomienie: ``python ex148_lyapunov_tgp_vs_newton_pairwise.py [--quick]``
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np

_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

from nbody.dynamics_v2 import (
    analytical_forces_tgp_pairwise,
    forces_newton,
    potential_newton,
    potential_tgp,
)
from nbody.lyapunov import (
    largest_lyapunov_exponent_benettin,
    pythagorean_three_body_burrau,
    random_velocities_for_excess_energy,
    scale_velocities_match_energy,
    total_mechanical_energy,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--quick",
        action="store_true",
        help="krotszy czas i wiekszy dt (regresja ~ sekundy)",
    )
    args = parser.parse_args()

    soft = 1e-6
    beta = 0.08
    gamma = beta
    G = 1.0

    pos0, vel0, C = pythagorean_three_body_burrau()

    def acc_newton(p: np.ndarray, c: np.ndarray) -> np.ndarray:
        return forces_newton(p, c, G=G, softening=soft) / c[:, None]

    def acc_tgp(p: np.ndarray, c: np.ndarray) -> np.ndarray:
        return analytical_forces_tgp_pairwise(p, c, beta, gamma, soft) / c[:, None]

    def pot_n(p: np.ndarray, c: np.ndarray) -> float:
        return potential_newton(p, c, G=G, softening=soft)

    def pot_t(p: np.ndarray, c: np.ndarray) -> float:
        return potential_tgp(p, c, beta, gamma, soft)

    E_new = total_mechanical_energy(pos0, vel0, C, pot_n)
    rng_match = np.random.default_rng(151)
    vel_tgp_matched = scale_velocities_match_energy(
        pos0, vel0, C, pot_t, target_energy=E_new
    )
    if float(np.sum(vel_tgp_matched**2)) < 1e-28:
        vel_tgp_matched = random_velocities_for_excess_energy(
            pos0, C, pot_t, target_energy=E_new, rng=rng_match
        )
    E_tgp_m = total_mechanical_energy(pos0, vel_tgp_matched, C, pot_t)
    if abs(E_tgp_m - E_new) > 1e-6 * max(1.0, abs(E_new)):
        raise RuntimeError("energy match failed for TGP branch")

    if args.quick:
        t_final = 18.0
        dt = 0.008
        renorm = 18
        jac_eps = 8e-5
    else:
        t_final = 72.0
        dt = 0.004
        renorm = 22
        jac_eps = 5e-5

    rng = np.random.default_rng(148)

    lam_n, sn = largest_lyapunov_exponent_benettin(
        pos0,
        vel0,
        C,
        acc_newton,
        t_final=t_final,
        dt=dt,
        renorm_every=renorm,
        jac_eps=jac_eps,
        rng=rng,
    )

    rng2 = np.random.default_rng(149)
    lam_t_samev, st1 = largest_lyapunov_exponent_benettin(
        pos0,
        vel0,
        C,
        acc_tgp,
        t_final=t_final,
        dt=dt,
        renorm_every=renorm,
        jac_eps=jac_eps,
        rng=rng2,
    )

    rng3 = np.random.default_rng(150)
    lam_t_match, st2 = largest_lyapunov_exponent_benettin(
        pos0,
        vel_tgp_matched,
        C,
        acc_tgp,
        t_final=t_final,
        dt=dt,
        renorm_every=renorm,
        jac_eps=jac_eps,
        rng=rng3,
    )

    print("ex148: largest Lyapunov exponent (Benettin + RK4 tangent), Pythagorean 3-body")
    print(f"  softening={soft}  beta=gamma={beta}  t_final={t_final}  dt={dt}")
    print(f"  Newton E0={E_new:.6g}")
    print(f"  TGP same v: E0={total_mechanical_energy(pos0, vel0, C, pot_t):.6g}")
    print(f"  TGP matched E: E={E_tgp_m:.6g} (target {E_new:.6g})")
    print(f"  steps: newton={sn}  tgp_sv={st1}  tgp_m={st2}")
    print(f"  lambda_max Newton:              {lam_n:.6g}")
    print(f"  lambda_max TGP (same v0=0):     {lam_t_samev:.6g}")
    print(f"  lambda_max TGP (E matched):     {lam_t_match:.6g}")

    ok = (
        np.isfinite(lam_n)
        and np.isfinite(lam_t_samev)
        and np.isfinite(lam_t_match)
        and lam_n > 0.01
    )
    if not ok:
        raise SystemExit(1)
    print("  PASS (finite lambda, Newton chaotic lambda > 0.01)")


if __name__ == "__main__":
    main()
