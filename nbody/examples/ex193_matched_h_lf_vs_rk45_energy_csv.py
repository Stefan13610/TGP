#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex193_matched_h_lf_vs_rk45_energy_csv.py
========================================
P1.C (cd.): **cztery** wiersze jak ``ex192`` (leapfrog vs RK45 ``max|(E-E0)/E0|`` przy
wspolnym ``t_final``), zapis do CSV w ``examples/_outputs/`` — te same obliczenia co
``compute_matched_h_lf_vs_rk_energy_rows`` w ``ex192_matched_h_leapfrog_vs_rk45_energy_table.py``.

Odczyt i walidacja: ``ex194_summarize_ex193_matched_h_lf_vs_rk_csv.py``.

Uruchomienie:
``python ex193_matched_h_lf_vs_rk45_energy_csv.py [--quick] [--n-quad N] [--out PATH]``
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import os
import sys

_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

_HERE = os.path.dirname(os.path.abspath(__file__))


def _load_ex192():
    path = os.path.join(_HERE, "ex192_matched_h_leapfrog_vs_rk45_energy_table.py")
    spec = importlib.util.spec_from_file_location("_ex192_lf_rk", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load ex192 from {path}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--n-quad", type=int, default=0, help="n_quad Yukawa (0 = domyslnie)")
    parser.add_argument(
        "--out",
        default="",
        help="sciezka CSV (domyslnie examples/_outputs/ex193_matched_h_lf_vs_rk_energy.csv)",
    )
    args = parser.parse_args()

    ex192 = _load_ex192()
    try:
        rows, meta = ex192.compute_matched_h_lf_vs_rk_energy_rows(
            quick=bool(args.quick), n_quad_override=int(args.n_quad)
        )
    except ValueError as e:
        print(f"ex193: FAIL: {e}")
        raise SystemExit(1)

    out_dir = os.path.join(_HERE, "_outputs")
    os.makedirs(out_dir, exist_ok=True)
    out_path = args.out or os.path.join(out_dir, "ex193_matched_h_lf_vs_rk_energy.csv")

    soft = float(meta["soft"])
    beta = float(meta["beta"])
    gamma = float(meta["gamma"])
    n_quad = int(meta["n_quad"])
    t_final = float(meta["t_final"])
    dt = float(meta["dt"])
    lf_se = int(meta["lf_save_every"])
    rk_no = int(meta["rk_n_output"])
    rtol = float(meta["rtol"])
    atol = float(meta["atol"])
    lf_tol = float(meta["lf_max_tol"])
    nrk_tol = float(meta["newton_rk_max_tol"])
    Hn = float(meta["H_N_ref"])

    with open(out_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(
            [
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
            ]
        )
        for lab, lf_m, rk_m, rk_succ, _ in rows:
            w.writerow(
                [
                    f"{soft:.16e}",
                    f"{beta:.16e}",
                    f"{gamma:.16e}",
                    n_quad,
                    f"{t_final:.16e}",
                    f"{dt:.16e}",
                    lf_se,
                    rk_no,
                    f"{rtol:.16e}",
                    f"{atol:.16e}",
                    f"{lf_tol:.16e}",
                    f"{nrk_tol:.16e}",
                    f"{Hn:.16e}",
                    lab,
                    f"{lf_m:.16e}",
                    f"{rk_m:.16e}",
                    1 if rk_succ else 0,
                ]
            )

    print(
        "ex193: CSV export (ex192 matched-H LF vs RK energy, four rows)\n"
        f"  t_final={t_final} dt={dt} lf_save_every={lf_se} rk_n_output={rk_no} "
        f"n_quad(Yukawa)={n_quad}"
    )
    print(f"  wrote {out_path}")

    if not ex192.matched_h_lf_vs_rk_energy_passes(rows, meta):
        print("  FAIL: thresholds (see ex192/ex190/ex191)")
        raise SystemExit(1)
    print("  PASS (same criteria as ex192)")


if __name__ == "__main__":
    main()
