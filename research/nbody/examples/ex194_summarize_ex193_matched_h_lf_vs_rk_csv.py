#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex194_summarize_ex193_matched_h_lf_vs_rk_csv.py
===============================================
P1.C (cd.): **odczyt** CSV z ``ex193`` (cztery wiersze: leapfrog vs RK45
``max|(E-E0)/E0|`` przy wspolnym ``t_final``, jak ``ex192``). Sprawdza spojnosc
metadanych miedzy wierszami oraz te same progi co ``ex192`` (``lf_max_tol``,
``rk_success``, Newton ``rk_max`` wzgledem ``newton_rk_max_tol``).

Domyslne wejscie: ``examples/_outputs/ex193_matched_h_lf_vs_rk_energy.csv``
(najpierw ``ex193`` lub cala regresja ``verify_nbody_lyapunov_quick``).

Analog roli do ``ex187`` wzgledem ``ex186``.

Uruchomienie:
``python ex194_summarize_ex193_matched_h_lf_vs_rk_csv.py [--quick] [--csv PATH]``
"""

from __future__ import annotations

import argparse
import csv
import os

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_DEFAULT_CSV = os.path.join(_HERE, "_outputs", "ex193_matched_h_lf_vs_rk_energy.csv")

_BRANCH_ORDER = (
    "Newton",
    "yukawa_feynman",
    "coulomb_3b",
    "pairwise V2",
)

_NEED_COLS = frozenset(
    {
        "soft",
        "beta",
        "gamma",
        "n_quad",
        "t_final",
        "dt",
        "lf_save_every",
        "rk_n_output",
        "rtol",
        "atol",
        "lf_max_tol",
        "newton_rk_max_tol",
        "H_N_ref",
        "branch",
        "lf_max_abs_dE_E0",
        "rk_max_abs_dE_E0",
        "rk_success",
    }
)


def _load_rows(path: str) -> list[dict[str, float | int | str]]:
    with open(path, newline="", encoding="utf-8") as f:
        r = csv.DictReader(f)
        if r.fieldnames is None or not _NEED_COLS.issubset(set(r.fieldnames)):
            raise ValueError(f"CSV {path}: expected columns {sorted(_NEED_COLS)}")
        out: list[dict[str, float | int | str]] = []
        for row in r:
            out.append(
                {
                    "soft": float(row["soft"]),
                    "beta": float(row["beta"]),
                    "gamma": float(row["gamma"]),
                    "n_quad": int(row["n_quad"]),
                    "t_final": float(row["t_final"]),
                    "dt": float(row["dt"]),
                    "lf_save_every": int(row["lf_save_every"]),
                    "rk_n_output": int(row["rk_n_output"]),
                    "rtol": float(row["rtol"]),
                    "atol": float(row["atol"]),
                    "lf_max_tol": float(row["lf_max_tol"]),
                    "newton_rk_max_tol": float(row["newton_rk_max_tol"]),
                    "H_N_ref": float(row["H_N_ref"]),
                    "branch": str(row["branch"]),
                    "lf_max_abs_dE_E0": float(row["lf_max_abs_dE_E0"]),
                    "rk_max_abs_dE_E0": float(row["rk_max_abs_dE_E0"]),
                    "rk_success": int(row["rk_success"]),
                }
            )
    return out


def _meta_consistent(rows: list[dict[str, float | int | str]]) -> bool:
    r0 = rows[0]
    keys = (
        "soft",
        "beta",
        "gamma",
        "n_quad",
        "t_final",
        "dt",
        "lf_save_every",
        "rk_n_output",
        "rtol",
        "atol",
        "lf_max_tol",
        "newton_rk_max_tol",
        "H_N_ref",
    )
    for r in rows[1:]:
        for k in keys:
            if k in ("n_quad", "lf_save_every", "rk_n_output"):
                if r[k] != r0[k]:
                    return False
            else:
                if not np.isclose(float(r[k]), float(r0[k]), rtol=0.0, atol=1e-12):
                    return False
    return True


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--csv", default="", help=f"default: {_DEFAULT_CSV}")
    args = parser.parse_args()

    csv_path = args.csv or _DEFAULT_CSV

    print("ex194: read ex193 matched-H LF vs RK energy CSV + validation")
    if not os.path.isfile(csv_path):
        print(f"  FAIL: missing {csv_path}")
        print("  Run: python nbody/examples/ex193_matched_h_lf_vs_rk45_energy_csv.py --quick")
        raise SystemExit(1)

    try:
        rows = _load_rows(csv_path)
    except ValueError as e:
        print(f"  FAIL: {e}")
        raise SystemExit(1)

    if len(rows) != 4:
        print(f"  FAIL: expected 4 data rows, got {len(rows)}")
        raise SystemExit(1)

    by_branch = {str(r["branch"]): r for r in rows}
    if set(by_branch) != set(_BRANCH_ORDER):
        print("  FAIL: branch set mismatch")
        print(f"    expected: {list(_BRANCH_ORDER)}")
        print(f"    got:      {sorted(by_branch)}")
        raise SystemExit(1)

    rows_ord = [by_branch[b] for b in _BRANCH_ORDER]

    if not _meta_consistent(rows_ord):
        print("  FAIL: inconsistent grid metadata across rows")
        raise SystemExit(1)

    lf_tol = float(rows_ord[0]["lf_max_tol"])
    nrk_tol = float(rows_ord[0]["newton_rk_max_tol"])

    for r in rows_ord:
        lf_m = float(r["lf_max_abs_dE_E0"])
        rk_m = float(r["rk_max_abs_dE_E0"])
        ok_rk = int(r["rk_success"]) == 1
        br = str(r["branch"])
        if not np.isfinite(lf_m) or lf_m >= lf_tol:
            print(f"  FAIL: LF drift for {br}: {lf_m} (tol {lf_tol})")
            raise SystemExit(1)
        if not np.isfinite(rk_m) or not ok_rk:
            print(f"  FAIL: RK45 for {br}: max_dE={rk_m} success={ok_rk}")
            raise SystemExit(1)
        if br == "Newton" and rk_m > nrk_tol:
            print(f"  FAIL: Newton RK max|dE|={rk_m} > {nrk_tol}")
            raise SystemExit(1)

    r0 = rows_ord[0]
    print(f"  csv={csv_path}")
    print(
        f"  soft={r0['soft']:.6g} beta=gamma={r0['beta']:.6g} n_quad={r0['n_quad']} "
        f"t_final={r0['t_final']} dt={r0['dt']} lf_save_every={r0['lf_save_every']} "
        f"rk_n_output={r0['rk_n_output']}"
    )
    print("")
    wbr, wlf, wrk = 18, 14, 14
    print(f"  {'branch':<{wbr}} {'LF max|dE|':>{wlf}} {'RK max|dE|':>{wrk}}")
    print(f"  {'-' * wbr} {'-' * wlf} {'-' * wrk}")
    for r in rows_ord:
        print(
            f"  {str(r['branch']):<{wbr}} "
            f"{float(r['lf_max_abs_dE_E0']):>{wlf}.4e} "
            f"{float(r['rk_max_abs_dE_E0']):>{wrk}.4e}"
        )

    print("")
    print("  Interpretation: same four matched-H branches as ex178/ex186; energy")
    print("  drift columns mirror ex192 (common t_final for LF and RK45).")

    print("")
    print("ex194: PASS")


if __name__ == "__main__":
    main()
