#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex174_lyapunov_newton_vs_yukawa_feynman_matched_energy.py
=========================================================
P1 (cd. / P1.C): **Newton** vs **``yukawa_feynman``** na tych samych pozycjach
startowych Burrau; dla Yukawy predkosci losowane tak, by **calkowita energia
w potencjale Yukawy** rownala sie energii Newtona przy ``v=0`` (jak dopasowanie
E w ``ex148``, lecz z pelnym ``V_3`` Feynmana).

Obie galezie: **leapfrog+tangent** + **analityczny** ``J`` (Newton:
``acceleration_jacobian_newton_softened``; Yukawa:
``acceleration_jacobian_yukawa_feynman_analytic``).

To nie jest dowod „tlumienia chaosu” — tylko **porownawcza** regresja na wspolnym
skalowaniu energetycznym przy roznych silach; interpretacja w ``PLAN_ROZWOJU_NBODY.md`` (P1.C).

Trzy konteksty hamiltonowskie na tym samym ``x`` (Newton ``H_N``, Yukawa ``v=0``, Yukawa z ``H_Y=H_N``):
``ex175_lyapunov_burrau_newton_yukawa_three_hamiltonians.py``.
Analog dla backendu ``coulomb_3b`` (FD ``J``): ``ex176_lyapunov_newton_vs_coulomb3b_matched_energy.py``.
Tylko ``V_2`` (pairwise), oba ``J`` analityczne: ``ex177_lyapunov_newton_vs_pairwise_matched_energy.py``.
Zbiorcza tabela (Newton + wszystkie trzy backendy TGP przy ``H=H_N``): ``ex178_lyapunov_matched_energy_family_table.py``; CSV: ``ex186_lyapunov_matched_h_family_csv.py``; odczyt CSV: ``ex187_summarize_ex186_matched_h_csv.py``; skan ``t_final``: ``ex188_lyapunov_matched_h_family_t_final_scan_csv.py``; odczyt skanu: ``ex189_summarize_ex188_matched_h_t_scan_csv.py``. Krotka diagnostyka leapfrog ``|dE/E0|`` (cztery galezie): ``ex190_leapfrog_energy_drift_matched_h_family_short.py``. RK45: ``ex191_rk45_energy_diag_matched_h_family_short.py``. Tabela LF/RK: ``ex192_matched_h_leapfrog_vs_rk45_energy_table.py``; CSV: ``ex193_matched_h_lf_vs_rk45_energy_csv.py``; odczyt: ``ex194_summarize_ex193_matched_h_lf_vs_rk_csv.py``.
Przeglad ``v=0`` vs ``H=H_N`` dla wszystkich backendow: ``ex181_lyapunov_matched_energy_seven_row_overview.py``; CSV: ``ex182_lyapunov_seven_row_overview_csv.py``; odczyt CSV: ``ex183_summarize_ex182_seven_row_csv.py``; skan ``t_final``: ``ex184_lyapunov_seven_row_t_final_scan_csv.py``; odczyt skanu: ``ex185_summarize_ex184_t_scan_csv.py``.
Trzy hamiltoniany Coulomba na tym samym ``x`` (jak ``ex175`` dla Yukawy): ``ex179_lyapunov_burrau_newton_coulomb_three_hamiltonians.py``.
Trzy hamiltoniany pairwise ``V_2``: ``ex180_lyapunov_burrau_newton_pairwise_three_hamiltonians.py``.

Uruchomienie:
``python ex174_lyapunov_newton_vs_yukawa_feynman_matched_energy.py [--quick] [--n-quad N]``
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
    parser.add_argument(
        "--n-quad",
        type=int,
        default=0,
        help="n_quad Feynman (0 = domyslne wg trybu)",
    )
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

    E_newton = total_mechanical_energy(pos0, vel0, C, pot_newton)

    acc_y, pot_y = build_tgp_integration_pair(
        "yukawa_feynman",
        beta=beta,
        gamma=gamma,
        softening=soft,
        include_3body=True,
        n_quad_feynman=n_quad,
    )

    rng_v = np.random.default_rng(1740)
    vel_y = random_velocities_for_excess_energy(
        pos0, C, pot_y, target_energy=E_newton, rng=rng_v
    )
    E_y = total_mechanical_energy(pos0, vel_y, C, pot_y)
    tol_e = 1e-5 * max(1.0, abs(E_newton))
    if abs(E_y - E_newton) > tol_e:
        print(f"  FAIL: energy match |E_y-E_n|={abs(E_y - E_newton):.6e} > {tol_e:.6e}")
        raise SystemExit(1)

    def acc_newton(p: np.ndarray, c: np.ndarray) -> np.ndarray:
        return forces_newton(p, c, G=G, softening=soft) / c[:, None]

    def jac_newton(p: np.ndarray, c: np.ndarray) -> np.ndarray:
        return acceleration_jacobian_newton_softened(p, c, G=G, softening=soft)

    def jac_yukawa(p: np.ndarray, c: np.ndarray, nq: int = n_quad) -> np.ndarray:
        return acceleration_jacobian_yukawa_feynman_analytic(
            p,
            c,
            beta=beta,
            gamma=gamma,
            softening=soft,
            n_quad_feynman=nq,
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
    lam_y, sy = largest_lyapunov_exponent_benettin_leapfrog(
        pos0,
        vel_y,
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
        "ex174: lambda_max Newton vs yukawa_feynman (Burrau x, E matched in Yukawa pot, "
        "leapfrog+analytic J)"
    )
    print(
        f"  soft={soft} beta=gamma={beta} n_quad={n_quad} "
        f"t_final={t_final} dt={dt} renorm={renorm}"
    )
    print(f"  E_newton(v=0)={E_newton:.8g}  E_yukawa(matched)={E_y:.8g}")
    print(f"  lambda_max Newton:          {lam_n:.8g}  steps={sn}")
    print(f"  lambda_max yukawa_feynman:  {lam_y:.8g}  steps={sy}")

    ok = (
        np.isfinite(lam_n)
        and np.isfinite(lam_y)
        and lam_n > 0.01
        and lam_y > 0.02
        and sn == sy
    )
    if not ok:
        raise SystemExit(1)
    print("  PASS (finite, Newton > 0.01, Yukawa > 0.02, same step count)")


if __name__ == "__main__":
    main()
