#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex179_lyapunov_burrau_newton_coulomb_three_hamiltonians.py
==========================================================
P1 (cd. / P1.C): analog ``ex175`` dla backendu **``coulomb_3b``** zamiast
``yukawa_feynman``. Na tych samych pozycjach Burrau trzy kontekty:

1. **Newton** z ``v=0`` — ``H_N``.
2. **``coulomb_3b``** z **tym samym** ``v=0`` — ``H_C`` w potencjale Coulomba
   (gleboka studnia wzgledem ``H_N``; inna warstwa energetyczna).
3. **``coulomb_3b``** z ``v`` dopasowanym do ``H_N`` (seed prędkości ``1760`` jak w
   ``ex176`` / ``ex178``).

Wszystkie biegi: **leapfrog+tangent**. Newton: **analityczny** ``J``;
``coulomb_3b``: **FD** ``J`` (``jac_eps=1e-4``), jak ``ex159`` / ``ex176``.

Dla (1) i (3) te same RNG styczne co w ``ex176`` (``1761``, ``1762``) — zgodność
``lambda_max`` z wierszami Newton i ``coulomb_3b`` (matched) w ``ex176`` przy tej
samej siatce. Środkowy wiersz: styczny ``1792``.

Porownanie z rodzina matched-H: ``ex176``, ``ex178``.

Uruchomienie:
``python ex179_lyapunov_burrau_newton_coulomb_three_hamiltonians.py [--quick]``
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
    largest_lyapunov_exponent_benettin_leapfrog,
    pythagorean_three_body_burrau,
    random_velocities_for_excess_energy,
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
        jac_eps_newton = 1e-5
        jac_eps_fd = 1e-4
    else:
        t_final = 4.0
        dt = 0.04
        renorm = 10
        jac_eps_newton = 1e-5
        jac_eps_fd = 1e-4

    def pot_newton(p: np.ndarray, c: np.ndarray) -> float:
        return potential_newton(p, c, G=G, softening=soft)

    H_n = total_mechanical_energy(pos0, vel0, C, pot_newton)

    acc_c, pot_c = build_tgp_integration_pair(
        "coulomb_3b",
        beta=beta,
        gamma=gamma,
        softening=soft,
        include_3body=True,
    )

    H_c_same_v = total_mechanical_energy(pos0, vel0, C, pot_c)

    rng_v = np.random.default_rng(1760)
    vel_c_matched = random_velocities_for_excess_energy(
        pos0, C, pot_c, target_energy=H_n, rng=rng_v
    )
    H_c_match = total_mechanical_energy(pos0, vel_c_matched, C, pot_c)
    tol_e = 1e-5 * max(1.0, abs(H_n))
    if abs(H_c_match - H_n) > tol_e:
        print(f"  FAIL: |H_c_match-H_n|={abs(H_c_match - H_n):.6e} > {tol_e:.6e}")
        raise SystemExit(1)

    def acc_newton(p: np.ndarray, c: np.ndarray) -> np.ndarray:
        return forces_newton(p, c, G=G, softening=soft) / c[:, None]

    def jac_newton(p: np.ndarray, c: np.ndarray) -> np.ndarray:
        return acceleration_jacobian_newton_softened(p, c, G=G, softening=soft)

    lam_n, sn = largest_lyapunov_exponent_benettin_leapfrog(
        pos0,
        vel0,
        C,
        acc_newton,
        t_final=t_final,
        dt=dt,
        renorm_every=renorm,
        jac_eps=jac_eps_newton,
        position_jacobian_fn=jac_newton,
        rng=np.random.default_rng(1761),
    )
    lam_c0, s0 = largest_lyapunov_exponent_benettin_leapfrog(
        pos0,
        vel0,
        C,
        acc_c,
        t_final=t_final,
        dt=dt,
        renorm_every=renorm,
        jac_eps=jac_eps_fd,
        position_jacobian_fn=None,
        rng=np.random.default_rng(1792),
    )
    lam_cm, sm = largest_lyapunov_exponent_benettin_leapfrog(
        pos0,
        vel_c_matched,
        C,
        acc_c,
        t_final=t_final,
        dt=dt,
        renorm_every=renorm,
        jac_eps=jac_eps_fd,
        position_jacobian_fn=None,
        rng=np.random.default_rng(1762),
    )

    print(
        "ex179: Burrau x - Newton + coulomb_3b (three Hamiltonian contexts, "
        "leapfrog; Coulomb FD J)\n"
        f"  soft={soft} beta=gamma={beta} t_final={t_final} dt={dt} renorm={renorm}"
    )
    print(f"  H_N (Newton, v=0)        = {H_n:.8g}")
    print(f"  H_C (coulomb, same v=0)  = {H_c_same_v:.8g}")
    print(f"  H_C (matched to H_N)     = {H_c_match:.8g}")
    print(f"  lambda_max Newton v=0           = {lam_n:.8g}  steps={sn}")
    print(f"  lambda_max coulomb_3b v=0       = {lam_c0:.8g}  steps={s0}")
    print(f"  lambda_max coulomb_3b H_C=H_N   = {lam_cm:.8g}  steps={sm}")

    ok = (
        np.isfinite(lam_n)
        and np.isfinite(lam_c0)
        and np.isfinite(lam_cm)
        and lam_n > 0.01
        and lam_c0 > 0.02
        and lam_cm > 0.02
        and sn == s0 == sm
    )
    if not ok:
        raise SystemExit(1)
    print("  PASS (finite, chaotic-scale lambdas, identical step counts)")


if __name__ == "__main__":
    main()
