#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex178_lyapunov_matched_energy_family_table.py
=============================================
P1 (cd. / P1.C): **jedna** tabela ``lambda_max`` na Burrau ``x`` przy **wspolnym**
``H_N`` (Newton, ``v=0``): jeden bieg Newtona + trzy galezie TGP z **dopasowanym**
``H`` w swoim potencjale (te same IC predkosci co w ``ex174`` / ``ex176`` / ``ex177``).

| Branch | Potencjal / sila | ``J`` |
|--------|------------------|-------|
| Newton | Newton | analityczny |
| yukawa_feynman | pelny ``V_2+V_3`` Feynman | analityczny |
| coulomb_3b | ``V_2`` + ``I`` Coulomb | FD (``jac_eps`` wiekszy) |
| pairwise | tylko ``V_2`` | analityczny |

Siatka czasu jak ``ex174`` ``--quick`` / tryb pelny. Styczne: **wspolny** Newton
``1781``; TGP ``1782``--``1784`` (inne niz pojedyncze ``ex174``/``ex176``/``ex177``,
wiec liczby ``lambda_max`` moga sie roznic od tych skryptow mimo tych samych ``v``).

Szczegoly pojedynczych par: ``ex174``, ``ex176``, ``ex177``. Trojki kontekstow: ``ex175``, ``ex179``, ``ex180``. Pelny przeglad ``v=0`` vs ``H=H_N``: ``ex181``; CSV: ``ex182``; podsumowanie: ``ex183``; skan ``t_final``: ``ex184``; odczyt skanu: ``ex185``. Ta tabela (CSV): ``ex186_lyapunov_matched_h_family_csv.py``; odczyt: ``ex187_summarize_ex186_matched_h_csv.py``. Skan ``t_final`` (cztery wiersze): ``ex188_lyapunov_matched_h_family_t_final_scan_csv.py``; odczyt: ``ex189_summarize_ex188_matched_h_t_scan_csv.py``. Krotka diagnostyka leapfrog ``|dE/E0|`` (ta sama rodzina IC): ``ex190_leapfrog_energy_drift_matched_h_family_short.py`` (API ``matched_h_family_leapfrog_branches``). RK45 (DOP853) na krotkim ``t``: ``ex191_rk45_energy_diag_matched_h_family_short.py``. Tabela LF vs RK (wspolny ``t_final``): ``ex192_matched_h_leapfrog_vs_rk45_energy_table.py``; CSV: ``ex193_matched_h_lf_vs_rk45_energy_csv.py``; odczyt: ``ex194_summarize_ex193_matched_h_lf_vs_rk_csv.py``.

Uruchomienie:
``python ex178_lyapunov_matched_energy_family_table.py [--quick] [--n-quad N]``
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import Any, Callable

import numpy as np

LeapfrogAcc = Callable[[np.ndarray, np.ndarray], np.ndarray]
LeapfrogPot = Callable[[np.ndarray, np.ndarray], float]

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


def matched_h_family_leapfrog_branches(
    *, quick: bool, n_quad_override: int = 0
) -> tuple[list[tuple[str, np.ndarray, np.ndarray, LeapfrogAcc, LeapfrogPot]], dict[str, Any]]:
    """
    Cztery galezie jak w ``compute_matched_h_family_table``: dla kazdej
    ``(label, pos0, vel, acc_fn, pot_fn)`` zgodne z ``dynamics_v2.leapfrog_integrate``.

    ``meta``: m.in. ``soft``, ``beta``, ``gamma``, ``n_quad``, ``H_n``, ``C``.
    Uzywane w ``ex190`` (krotka diagnostyka energetyczna leapfroga).
    """
    soft = 1e-6
    G = 1.0
    beta = gamma = 0.08
    pos0, vel0, C = pythagorean_three_body_burrau()

    if quick:
        n_quad = 12
    else:
        n_quad = 14
    if n_quad_override > 0:
        n_quad = int(n_quad_override)

    def pot_newton(p: np.ndarray, c: np.ndarray) -> float:
        return potential_newton(p, c, G=G, softening=soft)

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
    vel_y = random_velocities_for_excess_energy(
        pos0, C, pot_y, target_energy=H_n, rng=np.random.default_rng(1740)
    )
    H_y = total_mechanical_energy(pos0, vel_y, C, pot_y)

    acc_c, pot_c = build_tgp_integration_pair(
        "coulomb_3b",
        beta=beta,
        gamma=gamma,
        softening=soft,
        include_3body=True,
    )
    vel_c = random_velocities_for_excess_energy(
        pos0, C, pot_c, target_energy=H_n, rng=np.random.default_rng(1760)
    )
    H_c = total_mechanical_energy(pos0, vel_c, C, pot_c)

    def pot_tgp_pair(p: np.ndarray, c: np.ndarray) -> float:
        return potential_tgp(p, c, beta, gamma, soft)

    vel_p = scale_velocities_match_energy(pos0, vel0, C, pot_tgp_pair, target_energy=H_n)
    if float(np.sum(vel_p**2)) < 1e-28:
        vel_p = random_velocities_for_excess_energy(
            pos0, C, pot_tgp_pair, target_energy=H_n, rng=np.random.default_rng(1770)
        )
    H_p = total_mechanical_energy(pos0, vel_p, C, pot_tgp_pair)

    for label, h in (("yukawa", H_y), ("coulomb", H_c), ("pairwise", H_p)):
        if abs(h - H_n) > tol_e:
            raise ValueError(
                f"|H_{label}-H_n|={abs(h - H_n):.6e} > tol_e={tol_e:.6e}"
            )

    def acc_newton(p: np.ndarray, c: np.ndarray) -> np.ndarray:
        return forces_newton(p, c, G=G, softening=soft) / c[:, None]

    def acc_pair(p: np.ndarray, c: np.ndarray) -> np.ndarray:
        return analytical_forces_tgp_pairwise(p, c, beta, gamma, soft) / c[:, None]

    branches: list[tuple[str, np.ndarray, np.ndarray, LeapfrogAcc, LeapfrogPot]] = [
        ("Newton", pos0, vel0, acc_newton, pot_newton),
        ("yukawa_feynman", pos0, vel_y, acc_y, pot_y),
        ("coulomb_3b", pos0, vel_c, acc_c, pot_c),
        ("pairwise V2", pos0, vel_p, acc_pair, pot_tgp_pair),
    ]
    meta: dict[str, Any] = {
        "soft": soft,
        "beta": beta,
        "gamma": gamma,
        "n_quad": n_quad,
        "H_n": float(H_n),
        "G": G,
        "C": C,
    }
    return branches, meta


def compute_matched_h_family_table(
    *,
    quick: bool,
    n_quad_override: int = 0,
    t_final_override: float | None = None,
) -> tuple[list[tuple[str, float, float, int]], dict[str, Any]]:
    """
    Cztery wiersze jak stdout ``ex178``: (branch, H_model, lambda_max, steps).
    ``H_model`` to energia w **tym** potencjale dla stanu uzytego w danym biegu.

    ``t_final_override``: nadpisanie horyzontu (``dt`` / ``renorm`` bez zmian); uzywane w ``ex188``.
    """
    branches, m = matched_h_family_leapfrog_branches(
        quick=quick, n_quad_override=n_quad_override
    )
    soft = float(m["soft"])
    beta = float(m["beta"])
    gamma = float(m["gamma"])
    G = float(m["G"])
    n_quad = int(m["n_quad"])
    H_n = float(m["H_n"])
    C = m["C"]

    if quick:
        t_final = 2.0
        dt = 0.05
        renorm = 7
        jac_eps_newton = 1e-5
        jac_eps_analytic = 1e-5
        jac_eps_fd = 1e-4
    else:
        t_final = 4.0
        dt = 0.04
        renorm = 10
        jac_eps_newton = 1e-5
        jac_eps_analytic = 1e-5
        jac_eps_fd = 1e-4
    if t_final_override is not None:
        t_final = float(t_final_override)
        if not np.isfinite(t_final) or t_final <= 0:
            raise ValueError(f"invalid t_final_override={t_final_override!r}")

    (_, pos0, vel0, acc_newton, _) = branches[0]
    (_, _, vel_y, acc_y, _) = branches[1]
    (_, _, vel_c, acc_c, _) = branches[2]
    (_, _, vel_p, acc_pair, _) = branches[3]

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
        jac_eps=jac_eps_newton,
        position_jacobian_fn=jac_newton,
        rng=np.random.default_rng(1781),
    )
    lam_y, sy = largest_lyapunov_exponent_benettin_leapfrog(
        pos0,
        vel_y,
        C,
        acc_y,
        t_final=t_final,
        dt=dt,
        renorm_every=renorm,
        jac_eps=jac_eps_analytic,
        position_jacobian_fn=jac_yukawa,
        rng=np.random.default_rng(1782),
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
        rng=np.random.default_rng(1783),
    )
    lam_p, sp = largest_lyapunov_exponent_benettin_leapfrog(
        pos0,
        vel_p,
        C,
        acc_pair,
        t_final=t_final,
        dt=dt,
        renorm_every=renorm,
        jac_eps=jac_eps_analytic,
        position_jacobian_fn=jac_pair,
        rng=np.random.default_rng(1784),
    )

    steps = (sn, sy, sc, sp)
    if len(set(steps)) != 1:
        raise ValueError(f"step counts differ: {steps}")

    pot_y_ref = branches[1][4]
    pot_c_ref = branches[2][4]
    pot_p_ref = branches[3][4]
    H_y = total_mechanical_energy(pos0, vel_y, C, pot_y_ref)
    H_c = total_mechanical_energy(pos0, vel_c, C, pot_c_ref)
    H_p = total_mechanical_energy(pos0, vel_p, C, pot_p_ref)

    rows: list[tuple[str, float, float, int]] = [
        ("Newton", float(H_n), float(lam_n), int(sn)),
        ("yukawa_feynman", float(H_y), float(lam_y), int(sy)),
        ("coulomb_3b", float(H_c), float(lam_c), int(sc)),
        ("pairwise V2", float(H_p), float(lam_p), int(sp)),
    ]
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
        "steps": int(sn),
    }
    return rows, meta


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--n-quad", type=int, default=0, help="n_quad Yukawa (0 = domyslnie)")
    args = parser.parse_args()

    try:
        rows, meta = compute_matched_h_family_table(
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
    H_n = float(meta["H_n"])

    lam_n = rows[0][2]
    lam_y = rows[1][2]
    lam_c = rows[2][2]
    lam_p = rows[3][2]
    sn = rows[0][3]

    print(
        "ex178: matched-H family (one Newton + three TGP), leapfrog Benettin\n"
        f"  soft={soft} beta=gamma={beta} n_quad(Yukawa)={n_quad} "
        f"t_final={t_final} dt={dt} renorm={renorm}\n"
        f"  H_N={H_n:.8g} (all TGP branches matched within tol)"
    )
    print("")
    print(f"  {'branch':<16} {'lambda_max':>14} {'steps':>6}")
    print(f"  {'-'*16} {'-'*14} {'-'*6}")
    sy = rows[1][3]
    sc = rows[2][3]
    sp = rows[3][3]
    print(f"  {'Newton':<16} {lam_n:14.8g} {sn:6d}")
    print(f"  {'yukawa_feynman':<16} {lam_y:14.8g} {sy:6d}")
    print(f"  {'coulomb_3b':<16} {lam_c:14.8g} {sc:6d}")
    print(f"  {'pairwise V2':<16} {lam_p:14.8g} {sp:6d}")

    ok = (
        np.isfinite(lam_n)
        and np.isfinite(lam_y)
        and np.isfinite(lam_c)
        and np.isfinite(lam_p)
        and lam_n > 0.01
        and lam_y > 0.02
        and lam_c > 0.02
        and lam_p > 0.02
        and sn == sy == sc == sp
    )
    if not ok:
        raise SystemExit(1)
    print("")
    print("  PASS (finite, thresholds, identical step counts)")


if __name__ == "__main__":
    main()
