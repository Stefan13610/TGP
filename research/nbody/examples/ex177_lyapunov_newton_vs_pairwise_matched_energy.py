#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex177_lyapunov_newton_vs_pairwise_matched_energy.py
===================================================
P1 (cd. / P1.C): **Newton** vs **TGP pairwise** (tylko ``V_2``, bez ``V_3``) na Burrau ``x``;
w galezi TGP predkosci tak dobrane, by **energia w ``potential_tgp``** rownala sie ``H_N``
przy Newtonie i ``v=0``. Ten sam schemat energetyczny co ``ex148`` (tam RK4+tangent i dluzszy
``t``); tu **leapfrog+tangent** i **analityczny** ``J`` dla obu galezi (Newton:
``acceleration_jacobian_newton_softened``; pairwise:
``acceleration_jacobian_tgp_pairwise_softened``).

Parametry czasu jak ``ex174`` / ``ex176`` (``--quick`` / pelny), zeby porownania ``lambda_max``
miedzy skryptami ``matched_energy`` byly na tej samej siatce.

Pelne ``V_3`` z dopasowaniem ``H``: ``ex174`` (Yukawa), ``ex176`` (``coulomb_3b``).
Trzy konteksty hamiltonowskie pairwise (Newton / pairwise ``v=0`` / ``H_P=H_N``): ``ex180``.

Uruchomienie:
``python ex177_lyapunov_newton_vs_pairwise_matched_energy.py [--quick]``
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

    def pot_tgp(p: np.ndarray, c: np.ndarray) -> float:
        return potential_tgp(p, c, beta, gamma, soft)

    H_n = total_mechanical_energy(pos0, vel0, C, pot_newton)

    vel_t_matched = scale_velocities_match_energy(
        pos0, vel0, C, pot_tgp, target_energy=H_n
    )
    if float(np.sum(vel_t_matched**2)) < 1e-28:
        rng_m = np.random.default_rng(1770)
        vel_t_matched = random_velocities_for_excess_energy(
            pos0, C, pot_tgp, target_energy=H_n, rng=rng_m
        )

    H_t = total_mechanical_energy(pos0, vel_t_matched, C, pot_tgp)
    tol_e = 1e-5 * max(1.0, abs(H_n))
    if abs(H_t - H_n) > tol_e:
        print(f"  FAIL: energy match |H_t-H_n|={abs(H_t - H_n):.6e} > {tol_e:.6e}")
        raise SystemExit(1)

    def acc_newton(p: np.ndarray, c: np.ndarray) -> np.ndarray:
        return forces_newton(p, c, G=G, softening=soft) / c[:, None]

    def acc_tgp(p: np.ndarray, c: np.ndarray) -> np.ndarray:
        return analytical_forces_tgp_pairwise(p, c, beta, gamma, soft) / c[:, None]

    def jac_newton(p: np.ndarray, c: np.ndarray) -> np.ndarray:
        return acceleration_jacobian_newton_softened(p, c, G=G, softening=soft)

    def jac_tgp(p: np.ndarray, c: np.ndarray) -> np.ndarray:
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
    lam_t, st = largest_lyapunov_exponent_benettin_leapfrog(
        pos0,
        vel_t_matched,
        C,
        acc_tgp,
        t_final=t_final,
        dt=dt,
        renorm_every=renorm,
        jac_eps=jac_eps,
        position_jacobian_fn=jac_tgp,
        rng=np.random.default_rng(1772),
    )

    print(
        "ex177: lambda_max Newton vs TGP pairwise (Burrau x, H matched in V2 pot, "
        "leapfrog+analytic J)"
    )
    print(
        f"  soft={soft} beta=gamma={beta} t_final={t_final} dt={dt} renorm={renorm} "
        f"jac_eps={jac_eps}"
    )
    print(f"  H_newton(v=0)={H_n:.8g}  H_pairwise(matched)={H_t:.8g}")
    print(f"  lambda_max Newton:   {lam_n:.8g}  steps={sn}")
    print(f"  lambda_max pairwise: {lam_t:.8g}  steps={st}")

    ok = (
        np.isfinite(lam_n)
        and np.isfinite(lam_t)
        and lam_n > 0.01
        and lam_t > 0.02
        and sn == st
    )
    if not ok:
        raise SystemExit(1)
    print("  PASS (finite, Newton > 0.01, pairwise > 0.02, same step count)")


if __name__ == "__main__":
    main()
