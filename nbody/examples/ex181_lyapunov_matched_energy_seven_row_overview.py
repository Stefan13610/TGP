#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex181_lyapunov_matched_energy_seven_row_overview.py
===================================================
P1 (cd. / P1.C): **siedem** wierszy ``lambda_max`` na Burrau ``x`` w jednej tabeli:
**Newton** oraz po **dwa** wiersze na kazdy z trzech backendow TGP (**``v=0``** vs
**``H`` dopasowane do ``H_N``**).

Pokazuje w jednym stdout, **dlaczego** porownanie przy tym samym ``v=0`` (glebokie
``H`` w TGP) rozni sie od porownania przy **wspolnym** ``H`` (``ex174``--``ex178``,
``ex175`` / ``ex179`` / ``ex180``).

- Yukawa / pairwise: dopasowanie jak ``ex174`` / ``ex177`` (seed prędkości ``1740`` /
  ``scale`` + ``1770``).
- Coulomb: seed ``1760`` (jak ``ex176`` / ``ex178`` / ``ex179``).

Siatka czasu jak ``ex178``. Styczne **1811**--**1817** (unikalne dla tego skryptu).

Zob. ``ex178`` (tylko wiersze matched ``H``; CSV: ``ex186``; odczyt: ``ex187``; skan ``t_final``: ``ex188``; odczyt: ``ex189``; leapfrog ``|dE/E0|``: ``ex190``; RK45: ``ex191``; tabela: ``ex192``; CSV LF/RK: ``ex193``; odczyt: ``ex194``), ``ex175``, ``ex179``, ``ex180``.
Eksport CSV (te same liczby): ``ex182_lyapunov_seven_row_overview_csv.py``;
odczyt: ``ex183_summarize_ex182_seven_row_csv.py``.
Skan ``t_final`` (CSV): ``ex184_lyapunov_seven_row_t_final_scan_csv.py``; odczyt:
``ex185_summarize_ex184_t_scan_csv.py``.

Uruchomienie:
``python ex181_lyapunov_matched_energy_seven_row_overview.py [--quick] [--n-quad N]``
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import Any

import numpy as np

_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

from nbody.dynamics_backends import build_tgp_integration_pair
from nbody.dynamics_v2 import (
    analytical_forces_tgp_pairwise,
    forces_newton,
    potential_newton,
    potential_tgp,
)
from nbody.lyapunov import (
    acceleration_jacobian_newton_softened,
    acceleration_jacobian_tgp_pairwise_softened,
    acceleration_jacobian_yukawa_feynman_analytic,
    largest_lyapunov_exponent_benettin_leapfrog,
    pythagorean_three_body_burrau,
    random_velocities_for_excess_energy,
    scale_velocities_match_energy,
    total_mechanical_energy,
)


def compute_seven_row_lyapunov_overview(
    *,
    quick: bool,
    n_quad_override: int = 0,
    t_final_override: float | None = None,
) -> tuple[list[tuple[str, float, float, int]], dict[str, Any]]:
    """
    Ta sama tabela co stdout ``ex181``: siedem wierszy (label, H_model, lambda_max, steps).
    ``meta`` ma parametry siatki (m.in. do CSV w ``ex182``).

    ``t_final_override``: nadpisanie horyzontu (``dt`` / ``renorm`` bez zmian — jak w trybie
    ``quick`` / pelnym); uzywane w ``ex184``.
    """
    soft = 1e-6
    G = 1.0
    beta = gamma = 0.08
    pos0, vel0, C = pythagorean_three_body_burrau()

    if quick:
        t_final = 2.0
        dt = 0.05
        renorm = 7
        jac_eps_newton = 1e-5
        jac_eps_analytic = 1e-5
        jac_eps_fd = 1e-4
        n_quad = 12
    else:
        t_final = 4.0
        dt = 0.04
        renorm = 10
        jac_eps_newton = 1e-5
        jac_eps_analytic = 1e-5
        jac_eps_fd = 1e-4
        n_quad = 14
    if n_quad_override > 0:
        n_quad = int(n_quad_override)
    if t_final_override is not None:
        t_final = float(t_final_override)
        if not np.isfinite(t_final) or t_final <= 0:
            raise ValueError(f"invalid t_final_override={t_final_override!r}")

    def pot_newton(p: np.ndarray, c: np.ndarray) -> float:
        return potential_newton(p, c, G=G, softening=soft)

    def pot_pair(p: np.ndarray, c: np.ndarray) -> float:
        return potential_tgp(p, c, beta, gamma, soft)

    H_n = total_mechanical_energy(pos0, vel0, C, pot_newton)
    tol_e = 1e-5 * max(1.0, abs(H_n))

    acc_y, pot_y = build_tgp_integration_pair(
        "yukawa_feynman",
        beta=beta,
        gamma=gamma,
        softening=soft,
        include_3body=True,
        n_quad_feynman=n_quad,
    )
    H_y_v0 = total_mechanical_energy(pos0, vel0, C, pot_y)
    vel_y_m = random_velocities_for_excess_energy(
        pos0, C, pot_y, target_energy=H_n, rng=np.random.default_rng(1740)
    )
    H_y_m = total_mechanical_energy(pos0, vel_y_m, C, pot_y)

    acc_c, pot_c = build_tgp_integration_pair(
        "coulomb_3b",
        beta=beta,
        gamma=gamma,
        softening=soft,
        include_3body=True,
    )
    H_c_v0 = total_mechanical_energy(pos0, vel0, C, pot_c)
    vel_c_m = random_velocities_for_excess_energy(
        pos0, C, pot_c, target_energy=H_n, rng=np.random.default_rng(1760)
    )
    H_c_m = total_mechanical_energy(pos0, vel_c_m, C, pot_c)

    H_p_v0 = total_mechanical_energy(pos0, vel0, C, pot_pair)
    vel_p_m = scale_velocities_match_energy(pos0, vel0, C, pot_pair, target_energy=H_n)
    if float(np.sum(vel_p_m**2)) < 1e-28:
        vel_p_m = random_velocities_for_excess_energy(
            pos0, C, pot_pair, target_energy=H_n, rng=np.random.default_rng(1770)
        )
    H_p_m = total_mechanical_energy(pos0, vel_p_m, C, pot_pair)

    for label, h in (
        ("yukawa matched", H_y_m),
        ("coulomb matched", H_c_m),
        ("pairwise matched", H_p_m),
    ):
        if abs(h - H_n) > tol_e:
            raise ValueError(
                f"|H-{label}-H_n|={abs(h - H_n):.6e} > tol_e={tol_e:.6e}"
            )

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

    def jac_yukawa(p: np.ndarray, c: np.ndarray) -> np.ndarray:
        return acceleration_jacobian_yukawa_feynman_analytic(
            p,
            c,
            beta=beta,
            gamma=gamma,
            softening=soft,
            n_quad_feynman=n_quad,
        )

    rows: list[tuple[str, float, float, int]] = []

    def run_row(
        label: str,
        h_model: float,
        pos: np.ndarray,
        vel: np.ndarray,
        acc,
        jac_eps: float,
        jac_fn,
        rng_seed: int,
    ) -> None:
        lam, st = largest_lyapunov_exponent_benettin_leapfrog(
            pos,
            vel,
            C,
            acc,
            t_final=t_final,
            dt=dt,
            renorm_every=renorm,
            jac_eps=jac_eps,
            position_jacobian_fn=jac_fn,
            rng=np.random.default_rng(rng_seed),
        )
        rows.append((label, h_model, float(lam), int(st)))

    run_row("Newton v=0", H_n, pos0, vel0, acc_newton, jac_eps_newton, jac_newton, 1811)
    run_row("Yukawa v=0", H_y_v0, pos0, vel0, acc_y, jac_eps_analytic, jac_yukawa, 1812)
    run_row("Yukawa H=H_N", H_y_m, pos0, vel_y_m, acc_y, jac_eps_analytic, jac_yukawa, 1813)
    run_row("Coulomb v=0", H_c_v0, pos0, vel0, acc_c, jac_eps_fd, None, 1814)
    run_row("Coulomb H=H_N", H_c_m, pos0, vel_c_m, acc_c, jac_eps_fd, None, 1815)
    run_row("Pairwise v=0", H_p_v0, pos0, vel0, acc_pair, jac_eps_analytic, jac_pair, 1816)
    run_row("Pairwise H=H_N", H_p_m, pos0, vel_p_m, acc_pair, jac_eps_analytic, jac_pair, 1817)

    steps = [r[3] for r in rows]
    if len(set(steps)) != 1:
        raise ValueError(f"step counts differ: {steps}")

    meta: dict[str, Any] = {
        "soft": soft,
        "beta": beta,
        "gamma": gamma,
        "n_quad": n_quad,
        "t_final": t_final,
        "dt": dt,
        "renorm": renorm,
        "H_n": H_n,
        "jac_eps_newton": jac_eps_newton,
        "jac_eps_analytic": jac_eps_analytic,
        "jac_eps_fd": jac_eps_fd,
        "steps": int(steps[0]),
    }
    return rows, meta


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--n-quad", type=int, default=0, help="n_quad Yukawa (0 = domyslnie)")
    args = parser.parse_args()

    try:
        rows, meta = compute_seven_row_lyapunov_overview(
            quick=bool(args.quick), n_quad_override=int(args.n_quad)
        )
    except ValueError as e:
        print(f"  FAIL: {e}")
        raise SystemExit(1)

    soft = float(meta["soft"])
    beta = float(meta["beta"])
    gamma = float(meta["gamma"])
    n_quad = int(meta["n_quad"])
    t_final = float(meta["t_final"])
    dt = float(meta["dt"])
    renorm = int(meta["renorm"])
    st0 = int(meta["steps"])

    print(
        "ex181: seven-row Lyapunov overview (Newton + 3x TGP x {v=0, H=H_N})\n"
        f"  soft={soft} beta=gamma={beta} n_quad(Yukawa)={n_quad} "
        f"t_final={t_final} dt={dt} renorm={renorm}  steps={st0}"
    )
    print("")
    wlab, wH, wlam = 22, 14, 14
    print(f"  {'branch':<{wlab}} {'H_model':>{wH}} {'lambda_max':>{wlam}}")
    print(f"  {'-' * wlab} {'-' * wH} {'-' * wlam}")
    for lab, hm, lm, _ in rows:
        print(f"  {lab:<{wlab}} {hm:>{wH}.8g} {lm:>{wlam}.8g}")

    lams = [r[2] for r in rows]
    ok = (
        all(np.isfinite(lams))
        and lams[0] > 0.01
        and all(x > 0.02 for x in lams[1:])
    )
    if not ok:
        raise SystemExit(1)
    print("")
    print("  PASS (finite, Newton > 0.01, TGP rows > 0.02, uniform steps)")


if __name__ == "__main__":
    main()
