#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex175_lyapunov_burrau_newton_yukawa_three_hamiltonians.py
=========================================================
P1 (cd. / P1.C): na **tych samych pozycjach** Burrau trzy rozne kontekty
hamiltonowskie:

1. **Newton** z ``v=0`` — energia ``H_N`` w potencjale Newtona.
2. **``yukawa_feynman``** z **tym samym** ``v=0`` — energia ``H_Y`` w potencjale Yukawy
   (zwykle bardzo rozna od ``H_N``; to **nie** jest „ten sam stan'' fazowy w obu teorii).
3. **``yukawa_feynman``** z ``v`` z ``random_velocities_for_excess_energy`` tak, by
   ``H_Y \approx H_N`` (jak w ``ex174``).

Wszystkie trzy biegi: **leapfrog+tangent** + **analityczny** ``J`` (Newton:
``acceleration_jacobian_newton_softened``; Yukawa:
``acceleration_jacobian_yukawa_feynman_analytic``).

Cel: regresja + **czytelna diagnostyka** dla porownan Newton vs TGP na wspolnym ``x``:
bez dopasowania energii w galezi Yukawy ``lambda_max`` nie odpowiada temu samemu
skalowaniu co Newton (patrz ``ex174`` vs wiersz 2 ponizej).

Dla (1) i (3) uzyte sa te same RNG styczne co w ``ex174`` (``1741``, ``1742``), wiec
wartosci ``lambda_max`` w tych dwoch wierszach zgadzaja sie z ``ex174`` przy tych
samych ``dt`` / ``t_final`` / ``renorm``.

Newton vs ``coulomb_3b`` przy dopasowanym ``H`` (inny backend, FD ``J``): ``ex176``.
Trzy konteksty hamiltonowskie dla **Coulomba** (jak ten skrypt dla Yukawy): ``ex179``.
Dla **pairwise** ``V_2``: ``ex180_lyapunov_burrau_newton_pairwise_three_hamiltonians.py``.

Uruchomienie:
``python ex175_lyapunov_burrau_newton_yukawa_three_hamiltonians.py [--quick] [--n-quad N]``
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
from nbody.dynamics_v2 import forces_newton, potential_newton
from nbody.lyapunov import (
    acceleration_jacobian_newton_softened,
    acceleration_jacobian_yukawa_feynman_analytic,
    largest_lyapunov_exponent_benettin_leapfrog,
    pythagorean_three_body_burrau,
    random_velocities_for_excess_energy,
    total_mechanical_energy,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--n-quad", type=int, default=0, help="n_quad Feynman (0 = domyslnie)")
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
        n_quad = 12
    else:
        t_final = 4.0
        dt = 0.04
        renorm = 10
        jac_eps = 1e-5
        n_quad = 14
    if args.n_quad > 0:
        n_quad = int(args.n_quad)

    def pot_newton(p: np.ndarray, c: np.ndarray) -> float:
        return potential_newton(p, c, G=G, softening=soft)

    H_n = total_mechanical_energy(pos0, vel0, C, pot_newton)

    acc_y, pot_y = build_tgp_integration_pair(
        "yukawa_feynman",
        beta=beta,
        gamma=gamma,
        softening=soft,
        include_3body=True,
        n_quad_feynman=n_quad,
    )

    H_y_same_v = total_mechanical_energy(pos0, vel0, C, pot_y)

    rng_v = np.random.default_rng(1740)
    vel_y_matched = random_velocities_for_excess_energy(
        pos0, C, pot_y, target_energy=H_n, rng=rng_v
    )
    H_y_match = total_mechanical_energy(pos0, vel_y_matched, C, pot_y)
    tol_e = 1e-5 * max(1.0, abs(H_n))
    if abs(H_y_match - H_n) > tol_e:
        print(f"  FAIL: |H_y_match-H_n|={abs(H_y_match - H_n):.6e} > {tol_e:.6e}")
        raise SystemExit(1)

    def acc_newton(p: np.ndarray, c: np.ndarray) -> np.ndarray:
        return forces_newton(p, c, G=G, softening=soft) / c[:, None]

    def jac_newton(p: np.ndarray, c: np.ndarray) -> np.ndarray:
        return acceleration_jacobian_newton_softened(p, c, G=G, softening=soft)

    def jac_yukawa(p: np.ndarray, c: np.ndarray) -> np.ndarray:
        return acceleration_jacobian_yukawa_feynman_analytic(
            p,
            c,
            beta=beta,
            gamma=gamma,
            softening=soft,
            n_quad_feynman=n_quad,
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
        rng=np.random.default_rng(1741),
    )
    lam_y0, s0 = largest_lyapunov_exponent_benettin_leapfrog(
        pos0,
        vel0,
        C,
        acc_y,
        t_final=t_final,
        dt=dt,
        renorm_every=renorm,
        jac_eps=jac_eps,
        position_jacobian_fn=jac_yukawa,
        rng=np.random.default_rng(1752),
    )
    lam_ym, sm = largest_lyapunov_exponent_benettin_leapfrog(
        pos0,
        vel_y_matched,
        C,
        acc_y,
        t_final=t_final,
        dt=dt,
        renorm_every=renorm,
        jac_eps=jac_eps,
        position_jacobian_fn=jac_yukawa,
        rng=np.random.default_rng(1742),
    )

    print(
        "ex175: Burrau x - three Hamiltonian contexts (leapfrog+analytic J)\n"
        f"  soft={soft} beta=gamma={beta} n_quad={n_quad} "
        f"t_final={t_final} dt={dt} renorm={renorm}"
    )
    print(f"  H_N (Newton, v=0)     = {H_n:.8g}")
    print(f"  H_Y (Yukawa, same v=0)= {H_y_same_v:.8g}")
    print(f"  H_Y (matched to H_N)  = {H_y_match:.8g}")
    print(f"  lambda_max Newton v=0           = {lam_n:.8g}  steps={sn}")
    print(f"  lambda_max Yukawa v=0           = {lam_y0:.8g}  steps={s0}")
    print(f"  lambda_max Yukawa H_Y=H_N     = {lam_ym:.8g}  steps={sm}")

    ok = (
        np.isfinite(lam_n)
        and np.isfinite(lam_y0)
        and np.isfinite(lam_ym)
        and lam_n > 0.01
        and lam_y0 > 0.02
        and lam_ym > 0.02
        and sn == s0 == sm
    )
    if not ok:
        raise SystemExit(1)
    print("  PASS (finite, chaotic-scale lambdas, identical step counts)")


if __name__ == "__main__":
    main()
