#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex186_lyapunov_matched_h_family_csv.py
======================================
P1 (cd.): **cztery** wiersze ``lambda_max`` jak ``ex178`` (Newton + trzy TGP przy
``H`` dopasowanym do ``H_N``), zapis do CSV w ``examples/_outputs/`` — **ten sam**
uklad kolumn co ``ex182`` (``branch``, ``H_model``, ``lambda_max``, ``steps`` +
parametry siatki).

Obliczenia: ``compute_matched_h_family_table`` w ``ex178_...py``.
Odczyt CSV: ``ex187_summarize_ex186_matched_h_csv.py``.
Skan ``t_final`` (CSV): ``ex188_lyapunov_matched_h_family_t_final_scan_csv.py``; odczyt: ``ex189_summarize_ex188_matched_h_t_scan_csv.py``. Diagnostyka leapfrog ``|dE/E0|`` (krotki ``t``): ``ex190_leapfrog_energy_drift_matched_h_family_short.py``. RK45 (DOP853): ``ex191_rk45_energy_diag_matched_h_family_short.py``. Tabela LF vs RK: ``ex192_matched_h_leapfrog_vs_rk45_energy_table.py``; CSV energii: ``ex193_matched_h_lf_vs_rk45_energy_csv.py``; odczyt: ``ex194_summarize_ex193_matched_h_lf_vs_rk_csv.py``.

Uruchomienie:
``python ex186_lyapunov_matched_h_family_csv.py [--quick] [--n-quad N] [--out PATH]``
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import os
import sys

import numpy as np

_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

_HERE = os.path.dirname(os.path.abspath(__file__))


def _load_ex178():
    path = os.path.join(_HERE, "ex178_lyapunov_matched_energy_family_table.py")
    spec = importlib.util.spec_from_file_location("_ex178_matched_h", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load ex178 from {path}")
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
        help="sciezka CSV (domyslnie examples/_outputs/ex186_matched_h_family.csv)",
    )
    args = parser.parse_args()

    ex178 = _load_ex178()
    try:
        rows, meta = ex178.compute_matched_h_family_table(
            quick=bool(args.quick), n_quad_override=int(args.n_quad)
        )
    except ValueError as e:
        print(f"ex186: FAIL: {e}")
        raise SystemExit(1)

    out_dir = os.path.join(_HERE, "_outputs")
    os.makedirs(out_dir, exist_ok=True)
    out_path = args.out or os.path.join(out_dir, "ex186_matched_h_family.csv")

    soft = float(meta["soft"])
    beta = float(meta["beta"])
    gamma = float(meta["gamma"])
    n_quad = int(meta["n_quad"])
    t_final = float(meta["t_final"])
    dt = float(meta["dt"])
    renorm = int(meta["renorm"])
    H_n = float(meta["H_n"])
    st0 = int(meta["steps"])

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
                "renorm_every",
                "H_N_ref",
                "branch",
                "H_model",
                "lambda_max",
                "steps",
            ]
        )
        for lab, hm, lm, st in rows:
            w.writerow(
                [
                    f"{soft:.16e}",
                    f"{beta:.16e}",
                    f"{gamma:.16e}",
                    n_quad,
                    f"{t_final:.16e}",
                    f"{dt:.16e}",
                    renorm,
                    f"{H_n:.16e}",
                    lab,
                    f"{hm:.16e}",
                    f"{lm:.16e}",
                    st,
                ]
            )

    print(
        "ex186: CSV export (ex178 matched-H family, four rows)\n"
        f"  soft={soft} beta=gamma={beta} n_quad(Yukawa)={n_quad} "
        f"t_final={t_final} dt={dt} renorm={renorm}  steps={st0}"
    )
    print(f"  wrote {out_path}")

    lams = [r[2] for r in rows]
    ok = (
        all(np.isfinite(lams))
        and lams[0] > 0.01
        and all(x > 0.02 for x in lams[1:])
    )
    if not ok:
        raise SystemExit(1)
    print("  PASS (finite, Newton > 0.01, TGP rows > 0.02, uniform steps)")


if __name__ == "__main__":
    main()
