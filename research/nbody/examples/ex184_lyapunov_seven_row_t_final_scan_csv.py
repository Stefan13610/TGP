#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex184_lyapunov_seven_row_t_final_scan_csv.py
===========================================
P1.C (cd.): kilka wartosci ``t_final`` przy **tej samej** siatce ``dt`` / ``renorm_every``
co ``ex181`` (tryb ``--quick`` lub pelny). Dla kazdego ``t_final`` — siedem wierszy
jak ``ex181``/``ex182``; wynik: jeden CSV (7 x len(t_list) wierszy danych).

Obliczenia: ``compute_seven_row_lyapunov_overview(..., t_final_override=...)`` w ``ex181``.
Odczyt CSV: ``ex185_summarize_ex184_t_scan_csv.py``.
Analog czterowierszowy (matched ``H`` jak ``ex178``): ``ex188`` / ``ex189``.

Uruchomienie:
``python ex184_lyapunov_seven_row_t_final_scan_csv.py [--quick] [--n-quad N] [--out PATH]``
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

_T_FINAL_QUICK = (1.25, 2.0)
_T_FINAL_FULL = (2.0, 3.25, 4.0)


def _load_ex181():
    path = os.path.join(_HERE, "ex181_lyapunov_matched_energy_seven_row_overview.py")
    spec = importlib.util.spec_from_file_location("_ex181_lyap_seven_row", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load ex181 from {path}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _pass_slice(rows: list[tuple[str, float, float, int]]) -> bool:
    lams = [r[2] for r in rows]
    if not all(np.isfinite(lams)):
        return False
    if lams[0] <= 0.01:
        return False
    if any(x <= 0.02 for x in lams[1:]):
        return False
    steps = [r[3] for r in rows]
    return len(set(steps)) == 1


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--n-quad", type=int, default=0, help="n_quad Yukawa (0 = domyslnie)")
    parser.add_argument(
        "--out",
        default="",
        help="sciezka CSV (domyslnie examples/_outputs/ex184_seven_row_t_final_scan.csv)",
    )
    args = parser.parse_args()

    t_list = _T_FINAL_QUICK if args.quick else _T_FINAL_FULL
    ex181 = _load_ex181()

    out_dir = os.path.join(_HERE, "_outputs")
    os.makedirs(out_dir, exist_ok=True)
    out_path = args.out or os.path.join(out_dir, "ex184_seven_row_t_final_scan.csv")

    all_csv_rows: list[tuple] = []
    print(
        "ex184: t_final scan -> CSV (seven Lyapunov rows per horizon)\n"
        f"  t_final values ({'quick' if args.quick else 'full'}): {t_list}"
    )

    for tf in t_list:
        try:
            rows, meta = ex181.compute_seven_row_lyapunov_overview(
                quick=bool(args.quick),
                n_quad_override=int(args.n_quad),
                t_final_override=float(tf),
            )
        except ValueError as e:
            print(f"ex184: FAIL at t_final={tf}: {e}")
            raise SystemExit(1)

        if not _pass_slice(rows):
            print(f"ex184: FAIL thresholds or step mismatch at t_final={tf}")
            raise SystemExit(1)

        soft = float(meta["soft"])
        beta = float(meta["beta"])
        gamma = float(meta["gamma"])
        n_quad = int(meta["n_quad"])
        t_final = float(meta["t_final"])
        dt = float(meta["dt"])
        renorm = int(meta["renorm"])
        H_n = float(meta["H_n"])
        st0 = int(meta["steps"])
        lams = [float(r[2]) for r in rows]
        print(
            f"  t_final={t_final:g}  steps={st0}  "
            f"lambda_max in [{min(lams):.5g}, {max(lams):.5g}]"
        )

        for lab, hm, lm, st in rows:
            all_csv_rows.append(
                (
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
                )
            )

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
        for row in all_csv_rows:
            w.writerow(row)

    print(f"  wrote {out_path}  ({len(all_csv_rows)} data rows)")
    print("  PASS (each slice: finite, Newton > 0.01, TGP > 0.02, uniform steps)")


if __name__ == "__main__":
    main()
