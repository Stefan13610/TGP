#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex180_lyapunov_burrau_newton_pairwise_three_hamiltonians.py
===========================================================
P1 (cd. / P1.C): analog ``ex175`` (Yukawa) i ``ex179`` (Coulomb), lecz dla **TGP
pairwise** (tylko ``V_2``). Na tych samych pozycjach Burrau:

1. **Newton** z ``v=0`` — ``H_N``.
2. **Pairwise** z **tym samym** ``v=0`` — ``H_P`` w ``potential_tgp`` (gleboka
   studnia wzgledem ``H_N``; jak w ``ex148`` wiersz „same v'').
3. **Pairwise** z ``v`` takim, ze ``H_P \approx H_N`` (``scale_velocities_match_energy``,
   ewentualnie ``random_velocities_for_excess_energy`` z seedem ``1770`` jak w
   ``ex177`` / ``ex178``).

Wszystkie biegi: **leapfrog+tangent** + **analityczny** ``J`` (Newton:
``acceleration_jacobian_newton_softened``; pairwise:
``acceleration_jacobian_tgp_pairwise_softened``).

Dla (1) i (3) te same RNG styczne co w ``ex177`` (``1771``, ``1772``). Srodkowy
wiersz: styczny ``1802``.

Zob. ``ex148`` (RK4, dluzszy ``t``), ``ex177`` (para Newton vs matched), ``ex178``.

Uruchomienie:
``python ex180_lyapunov_burrau_newton_pairwise_three_hamiltonians.py [--quick]``
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
    acceleration_jacobian_newton_softened,
    acceleration_jacobian_tgp_pairwise_softened,
    largest_lyapunov_exponent_benettin_leapfrog,
    pythagorean_three_body_burrau,
    random_velocities_for_excess_energy,
    scale_velocities_match_energy,
    total_mechanical_energy,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true")
    args = parser.parse_args()

    soft = 1e-6
    G = 1.0
    beta = gamma = 0.08
    pos0, vel0, C = pythagorean_three_body_burrau()

    if args.quick:
        t_final = 2.0
        dt = 0.05
        renorm = 7
        jac_eps = 1e-5
    else:
        t_final = 4.0
        dt = 0.04
        renorm = 10
        jac_eps = 1e-5

    def pot_newton(p: np.ndarray, c: np.ndarray) -> float:
        return potential_newton(p, c, G=G, softening=soft)

    def pot_pair(p: np.ndarray, c: np.ndarray) -> float:
        return potential_tgp(p, c, beta, gamma, soft)

    H_n = total_mechanical_energy(pos0, vel0, C, pot_newton)
    H_p_same_v = total_mechanical_energy(pos0, vel0, C, pot_pair)

    vel_p_matched = scale_velocities_match_energy(
        pos0, vel0, C, pot_pair, target_energy=H_n
    )
    if float(np.sum(vel_p_matched**2)) < 1e-28:
        vel_p_matched = random_velocities_for_excess_energy(
            pos0, C, pot_pair, target_energy=H_n, rng=np.random.default_rng(1770)
        )
    H_p_match = total_mechanical_energy(pos0, vel_p_matched, C, pot_pair)
    tol_e = 1e-5 * max(1.0, abs(H_n))
    if abs(H_p_match - H_n) > tol_e:
        print(f"  FAIL: |H_P_match-H_n|={abs(H_p_match - H_n):.6e} > {tol_e:.6e}")
        raise SystemExit(1)

    def acc_newton(p: np.ndarray, c: np.ndarray) -> np.ndarray:
        return forces_newton(p, c, G=G, softening=soft) / c[:, None]

    def acc_pair(p: np.ndarray, c: np.ndarray) -> np.ndarray:
        return analytical_forces_tgp_pairwise(p, c, beta, gamma, soft) / c[:, None]

    def jac_newton(p: np.ndarray, c: np.ndarray) -> np.ndarray:
        return acceleration_jacobian_newton_softened(p, c, G=G, softening=soft)

    def jac_pair(p: np.ndarray, c: np.ndarray) -> np.ndarray:
        return acceleration_jacobian_tgp_pairwise_softened(
            p, c, beta=beta, gamma=gamma, softening=soft
        )

    lam_n, sn = largest_lyapunov_exponent_benettin_leapfrog(
        pos0,
        vel0,
        C,
        acc_newton,
        t_final=t_final,
        dt=dt,
        renorm_every=renorm,
        jac_eps=jac_eps,
        position_jacobian_fn=jac_newton,
        rng=np.random.default_rng(1771),
    )
    lam_p0, s0 = largest_lyapunov_exponent_benettin_leapfrog(
        pos0,
        vel0,
        C,
        acc_pair,
        t_final=t_final,
        dt=dt,
        renorm_every=renorm,
        jac_eps=jac_eps,
        position_jacobian_fn=jac_pair,
        rng=np.random.default_rng(1802),
    )
    lam_pm, sm = largest_lyapunov_exponent_benettin_leapfrog(
        pos0,
        vel_p_matched,
        C,
        acc_pair,
        t_final=t_final,
        dt=dt,
        renorm_every=renorm,
        jac_eps=jac_eps,
        position_jacobian_fn=jac_pair,
        rng=np.random.default_rng(1772),
    )

    print(
        "ex180: Burrau x - Newton + pairwise V2 (three Hamiltonian contexts, "
        "leapfrog+analytic J)\n"
        f"  soft={soft} beta=gamma={beta} t_final={t_final} dt={dt} renorm={renorm}"
    )
    print(f"  H_N (Newton, v=0)          = {H_n:.8g}")
    print(f"  H_P (pairwise, same v=0)   = {H_p_same_v:.8g}")
    print(f"  H_P (matched to H_N)       = {H_p_match:.8g}")
    print(f"  lambda_max Newton v=0            = {lam_n:.8g}  steps={sn}")
    print(f"  lambda_max pairwise v=0          = {lam_p0:.8g}  steps={s0}")
    print(f"  lambda_max pairwise H_P=H_N    = {lam_pm:.8g}  steps={sm}")

    ok = (
        np.isfinite(lam_n)
        and np.isfinite(lam_p0)
        and np.isfinite(lam_pm)
        and lam_n > 0.01
        and lam_p0 > 0.02
        and lam_pm > 0.02
        and sn == s0 == sm
    )
    if not ok:
        raise SystemExit(1)
    print("  PASS (finite, chaotic-scale lambdas, identical step counts)")


if __name__ == "__main__":
    main()
