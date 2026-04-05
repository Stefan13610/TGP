#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex192_matched_h_leapfrog_vs_rk45_energy_table.py
==============================================
P1.C (cd.): **jedna tabela** stdout: ``max|(E-E0)/E0|`` dla **tej samej** rodziny IC
co ``ex178`` — kolumna leapfrog (jak ``ex190``) vs kolumna RK45 DOP853 (jak ``ex191``)
przy **wspolnym** ``t_final`` w obu kolumnach (quick: ``0.1``; pelny: ``0.11``), zeby
porownanie bylo na tym samym horyzoncie czasu.

Parametry leapfrog: jak ``ex190`` (``dt``, ``save_every``). RK45: ``rtol``/``atol`` jak
``ex191``; ``n_output`` dopasowany do ``t_final``. Progi PASS: suma warunkow z ``ex190``
oraz ``ex191`` (Newton osobno dla RK45).

API: ``compute_matched_h_lf_vs_rk_energy_rows`` — uzywane w ``ex193`` (CSV).

Szczegoly interpretacji: ``tgp_lyapunov_benettin.tex``; pojedyncze kolumny: ``ex190``, ``ex191``.
Eksport CSV: ``ex193_matched_h_lf_vs_rk45_energy_csv.py``; odczyt: ``ex194_summarize_ex193_matched_h_lf_vs_rk_csv.py``.

Uruchomienie: ``python ex192_matched_h_leapfrog_vs_rk45_energy_table.py [--quick] [--n-quad N]``
"""

from __future__ import annotations

import argparse
import importlib.util
import os
import sys
from typing import Any

import numpy as np

_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

from nbody.dynamics_v2 import leapfrog_integrate, rk45_integrate

_HERE = os.path.dirname(os.path.abspath(__file__))


def _load_ex178():
    path = os.path.join(_HERE, "ex178_lyapunov_matched_energy_family_table.py")
    spec = importlib.util.spec_from_file_location("_ex178_matched_h", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load ex178 from {path}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _lf_max_drift(
    pos0: np.ndarray,
    vel0: np.ndarray,
    C: np.ndarray,
    acc,
    pot,
    *,
    t_final: float,
    dt: float,
    save_every: int,
) -> float:
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


def _rk45_max_drift(
    pos0: np.ndarray,
    vel0: np.ndarray,
    C: np.ndarray,
    acc,
    pot,
    *,
    t_final: float,
    n_output: int,
    rtol: float,
    atol: float,
) -> tuple[float, bool, str]:
    r = rk45_integrate(
        pos0,
        vel0,
        C,
        acc,
        pot,
        t_span=(0.0, float(t_final)),
        n_output=int(n_output),
        rtol=rtol,
        atol=atol,
        quiet=True,
    )
    ee = r["energy_error"]
    if not np.all(np.isfinite(ee)):
        return float("nan"), bool(r["success"]), str(r.get("message", ""))
    return float(np.max(np.abs(ee))), bool(r["success"]), str(r.get("message", ""))


def compute_matched_h_lf_vs_rk_energy_rows(
    *, quick: bool, n_quad_override: int = 0
) -> tuple[list[tuple[str, float, float, bool, str]], dict[str, Any]]:
    """
    Cztery wiersze: ``(branch, lf_max|dE/E0|, rk_max|dE/E0|, rk_success, rk_message)``.
    ``meta`` zawiera parametry siatki i progi (m.in. do CSV w ``ex193``).
    """
    ex178 = _load_ex178()
    branches, m0 = ex178.matched_h_family_leapfrog_branches(
        quick=quick, n_quad_override=int(n_quad_override)
    )

    n_quad = int(m0["n_quad"])
    C = m0["C"]
    rtol = 1e-9
    atol = 1e-11

    if quick:
        t_final = 0.1
        dt = 0.015
        save_every = 2
        lf_tol = 0.017
        n_output = 80
        newton_rk_tol = 0.001
    else:
        t_final = 0.11
        dt = 0.014
        save_every = 2
        lf_tol = 0.012
        n_output = 88
        newton_rk_tol = 0.0015

    rows: list[tuple[str, float, float, bool, str]] = []
    for lab, pos0, vel, acc, pot in branches:
        lf_m = _lf_max_drift(
            pos0, vel, C, acc, pot, t_final=t_final, dt=dt, save_every=save_every
        )
        rk_m, rk_succ, rk_msg = _rk45_max_drift(
            pos0,
            vel,
            C,
            acc,
            pot,
            t_final=t_final,
            n_output=n_output,
            rtol=rtol,
            atol=atol,
        )
        rows.append((lab, lf_m, rk_m, rk_succ, rk_msg))

    meta: dict[str, Any] = {
        "soft": float(m0["soft"]),
        "beta": float(m0["beta"]),
        "gamma": float(m0["gamma"]),
        "n_quad": n_quad,
        "H_N_ref": float(m0["H_n"]),
        "t_final": float(t_final),
        "dt": float(dt),
        "lf_save_every": int(save_every),
        "rk_n_output": int(n_output),
        "rtol": float(rtol),
        "atol": float(atol),
        "lf_max_tol": float(lf_tol),
        "newton_rk_max_tol": float(newton_rk_tol),
    }
    return rows, meta


def matched_h_lf_vs_rk_energy_passes(
    rows: list[tuple[str, float, float, bool, str]], meta: dict[str, Any]
) -> bool:
    lf_tol = float(meta["lf_max_tol"])
    newton_rk_tol = float(meta["newton_rk_max_tol"])
    for lab, lf_m, rk_m, rk_succ, _ in rows:
        if not np.isfinite(lf_m) or lf_m >= lf_tol:
            return False
        if not np.isfinite(rk_m) or not rk_succ:
            return False
        if lab == "Newton" and rk_m > newton_rk_tol:
            return False
    return True


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--n-quad", type=int, default=0, help="n_quad Yukawa (0 = domyslnie)")
    args = parser.parse_args()

    try:
        rows, meta = compute_matched_h_lf_vs_rk_energy_rows(
            quick=bool(args.quick), n_quad_override=int(args.n_quad)
        )
    except ValueError as e:
        print(f"ex192: FAIL: {e}")
        raise SystemExit(1)

    t_final = meta["t_final"]
    dt = meta["dt"]
    save_every = meta["lf_save_every"]
    n_output = meta["rk_n_output"]
    rtol = meta["rtol"]
    atol = meta["atol"]
    n_quad = meta["n_quad"]
    lf_tol = meta["lf_max_tol"]
    newton_rk_tol = meta["newton_rk_max_tol"]

    print(
        "ex192: matched-H family — leapfrog vs RK45 max |(E-E0)/E0| (common t_final)\n"
        f"  t_final={t_final}  leapfrog: dt={dt} save_every={save_every}  "
        f"RK45: n_output={n_output} rtol={rtol:g} atol={atol:g}  n_quad(Yukawa)={n_quad}"
    )

    wbr = 22
    wlf = 14
    wrk = 14
    wok = 6
    print("")
    print(
        f"  {'branch':<{wbr}} {'LF max|dE|':>{wlf}} {'RK max|dE|':>{wrk}} "
        f"{'RK ok':>{wok}}"
    )
    print(f"  {'-' * wbr} {'-' * wlf} {'-' * wrk} {'-' * wok}")

    all_ok = True
    for lab, lf_m, rk_m, rk_succ, rk_msg in rows:
        row_ok = True
        if not np.isfinite(lf_m) or lf_m >= lf_tol:
            row_ok = False
        if not np.isfinite(rk_m) or not rk_succ:
            row_ok = False
        if lab == "Newton" and np.isfinite(rk_m) and rk_m > newton_rk_tol:
            row_ok = False

        flag = "yes" if row_ok else "no"
        if not row_ok:
            all_ok = False
        print(f"  {lab:<{wbr}} {lf_m:>{wlf}.4e} {rk_m:>{wrk}.4e} {flag:>{wok}}")
        if not rk_succ and rk_msg:
            print(f"    (RK45 message: {rk_msg})")

    print("")
    if not all_ok:
        print(
            "ex192: FAIL (see ex190 thresholds for LF; ex191 for RK45 success "
            f"and Newton < {newton_rk_tol})"
        )
        raise SystemExit(1)
    print(
        f"ex192: PASS (LF max|dE/E0| < {lf_tol}; RK45 success; "
        f"Newton RK < {newton_rk_tol})"
    )


if __name__ == "__main__":
    main()
