#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex176_lyapunov_newton_vs_coulomb3b_matched_energy.py
=====================================================
P1 (cd. / P1.C): **Newton** vs **``coulomb_3b``** (``V_2`` + irreducible ``I`` Coulomb
w ``V_3``) na tych samych pozycjach Burrau; w galezi TGP predkosci losowane tak, by
**calkowita energia w potencjale ``coulomb_3b``** rownala sie ``H_N`` przy Newtonie
i ``v=0`` (ten sam schemat co ``ex174``, inny backend niz ``yukawa_feynman``).

- Newton: leapfrog+tangent + **analityczny** ``J`` (``acceleration_jacobian_newton_softened``);
  seed styczny ``1761`` (inny niz ``ex174``, wiec ``lambda_max`` Newtona rozni sie od wiersza Newtona w ``ex174`` przy tym samym ``dt``).
- ``coulomb_3b``: leapfrog+tangent + **FD** ``J`` (``position_jacobian_fn=None``; jak
  ``ex159``, ``ex171``). Dla FD uzyty jest nieco wiekszy ``jac_eps`` niz dla Newtona.

Regresja porownawcza (P1.C), nie teza o tlumieniu chaosu. Zob. tez ``ex174`` (Yukawa Feynman),
``ex175`` (trzy hamiltoniany Yukawy na tym samym ``x``), ``ex177`` (pairwise ``V_2`` + dopasowanie ``H``), ``ex179`` (trzy hamiltoniany Coulomba jak ``ex175``), ``ex180`` (trzy hamiltoniany pairwise), ``ex181`` (siedem wierszy), ``ex182`` (CSV), ``ex183`` (odczyt CSV), ``ex184`` (skan ``t_final``), ``ex185`` (odczyt skanu CSV), ``ex186`` (CSV tabeli ``ex178``), ``ex187`` (odczyt CSV ``ex186``), ``ex188`` (skan ``t_final``), ``ex189`` (odczyt CSV ``ex188``), ``ex190`` (leapfrog ``|dE/E0|``, krotki ``t``), ``ex191`` (RK45), ``ex192`` (tabela LF vs RK), ``ex193`` (CSV), ``ex194`` (odczyt CSV).

Uruchomienie:
``python ex176_lyapunov_newton_vs_coulomb3b_matched_energy.py [--quick]``
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

    rng_v = np.random.default_rng(1760)
    vel_c = random_velocities_for_excess_energy(
        pos0, C, pot_c, target_energy=H_n, rng=rng_v
    )
    H_c = total_mechanical_energy(pos0, vel_c, C, pot_c)
    tol_e = 1e-5 * max(1.0, abs(H_n))
    if abs(H_c - H_n) > tol_e:
        print(f"  FAIL: energy match |H_c-H_n|={abs(H_c - H_n):.6e} > {tol_e:.6e}")
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
    lam_c, sc = largest_lyapunov_exponent_benettin_leapfrog(
        pos0,
        vel_c,
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
        "ex176: lambda_max Newton vs coulomb_3b (Burrau x, H matched in Coulomb pot, "
        "leapfrog; Newton analytic J, coulomb FD J)"
    )
    print(
        f"  soft={soft} beta=gamma={beta} t_final={t_final} dt={dt} renorm={renorm} "
        f"jac_eps Newton={jac_eps_newton} FD={jac_eps_fd}"
    )
    print(f"  H_newton(v=0)={H_n:.8g}  H_coulomb(matched)={H_c:.8g}")
    print(f"  lambda_max Newton:     {lam_n:.8g}  steps={sn}")
    print(f"  lambda_max coulomb_3b: {lam_c:.8g}  steps={sc}")

    ok = (
        np.isfinite(lam_n)
        and np.isfinite(lam_c)
        and lam_n > 0.01
        and lam_c > 0.02
        and sn == sc
    )
    if not ok:
        raise SystemExit(1)
    print("  PASS (finite, Newton > 0.01, coulomb_3b > 0.02, same step count)")


if __name__ == "__main__":
    main()
