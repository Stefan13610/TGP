#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex183_summarize_ex182_seven_row_csv.py
======================================
P1 (cd.): **odczyt i podsumowanie** CSV z ``ex182`` (siedem wierszy ``lambda_max``).
Sprawdza spojnosc metadanych miedzy wierszami, progi jak ``ex181``/``ex182``, drukuje
roznice ``v=0`` vs ``H=H_N`` w obrebie kazdego backendu TGP. Opcjonalnie PNG (slupki;
wymaga ``matplotlib``).

Domyslne wejscie: ``examples/_outputs/ex182_seven_row_overview.csv``
(najpierw ``ex182`` lub cala regresja ``verify_nbody_lyapunov_quick``).

Analog roli do ``ex168`` wzgledem ``ex166``. Wielopunktowy skan ``t_final`` (inny plik CSV):
``ex184_lyapunov_seven_row_t_final_scan_csv.py``; podsumowanie: ``ex185_summarize_ex184_t_scan_csv.py``.
Cztery wiersze matched ``H`` (``ex186``): ``ex187_summarize_ex186_matched_h_csv.py``; skan ``t_final``: ``ex188_lyapunov_matched_h_family_t_final_scan_csv.py``; odczyt: ``ex189_summarize_ex188_matched_h_t_scan_csv.py``; leapfrog ``|dE/E0|``: ``ex190``; RK45: ``ex191``; tabela: ``ex192``; CSV LF/RK: ``ex193``; odczyt: ``ex194``.

``--quick``: stdout + walidacja (bez matplotlib).

Uruchomienie:
``python ex183_summarize_ex182_seven_row_csv.py [--quick] [--csv PATH] [--png PATH]``
"""

from __future__ import annotations

import argparse
import csv
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_DEFAULT_CSV = os.path.join(_HERE, "_outputs", "ex182_seven_row_overview.csv")
_DEFAULT_PNG = os.path.join(_HERE, "_outputs", "ex183_ex182_lambda_bar.png")

_BRANCH_ORDER = (
    "Newton v=0",
    "Yukawa v=0",
    "Yukawa H=H_N",
    "Coulomb v=0",
    "Coulomb H=H_N",
    "Pairwise v=0",
    "Pairwise H=H_N",
)

_NEED_COLS = frozenset(
    {
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
                    "renorm_every": int(row["renorm_every"]),
                    "H_N_ref": float(row["H_N_ref"]),
                    "branch": str(row["branch"]),
                    "H_model": float(row["H_model"]),
                    "lambda_max": float(row["lambda_max"]),
                    "steps": int(row["steps"]),
                }
            )
    return out


def _meta_consistent(rows: list[dict[str, float | int | str]]) -> bool:
    r0 = rows[0]
    for r in rows[1:]:
        if r["n_quad"] != r0["n_quad"] or r["renorm_every"] != r0["renorm_every"]:
            return False
        if r["steps"] != r0["steps"]:
            return False
        if not np.allclose(
            [r["soft"], r["beta"], r["gamma"], r["t_final"], r["dt"], r["H_N_ref"]],
            [r0["soft"], r0["beta"], r0["gamma"], r0["t_final"], r0["dt"], r0["H_N_ref"]],
            rtol=0.0,
            atol=1e-12,
        ):
            return False
    return True


def _try_bar_png(rows_ordered: list[dict[str, float | int | str]], png_path: str) -> bool:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as e:
        print(f"  (matplotlib not available — PNG skipped: {e})")
        return False

    labels = [str(r["branch"]) for r in rows_ordered]
    y = [float(r["lambda_max"]) for r in rows_ordered]
    fig, ax = plt.subplots(figsize=(9.0, 4.2), constrained_layout=True)
    x = np.arange(len(labels))
    ax.bar(x, y, color="steelblue", alpha=0.85)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=35, ha="right")
    ax.set_ylabel("lambda_max")
    ax.set_title("ex182 overview (Burrau x, seven branches)")
    ax.grid(True, axis="y", alpha=0.3)
    os.makedirs(os.path.dirname(os.path.abspath(png_path)) or ".", exist_ok=True)
    fig.savefig(png_path, dpi=140)
    plt.close(fig)
    print(f"  wrote {png_path}")
    return True


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--csv", default="", help=f"default: {_DEFAULT_CSV}")
    parser.add_argument("--png", default="", help=f"default: {_DEFAULT_PNG}")
    args = parser.parse_args()

    csv_path = args.csv or _DEFAULT_CSV
    png_path = args.png or _DEFAULT_PNG

    print("ex183: read ex182 seven-row CSV + summary (+ optional PNG)")
    if not os.path.isfile(csv_path):
        print(f"  FAIL: missing {csv_path}")
        print("  Run: python nbody/examples/ex182_lyapunov_seven_row_overview_csv.py --quick")
        raise SystemExit(1)

    try:
        rows = _load_rows(csv_path)
    except ValueError as e:
        print(f"  FAIL: {e}")
        raise SystemExit(1)

    if len(rows) != 7:
        print(f"  FAIL: expected 7 data rows, got {len(rows)}")
        raise SystemExit(1)

    by_branch = {str(r["branch"]): r for r in rows}
    if set(by_branch) != set(_BRANCH_ORDER):
        print(f"  FAIL: branch set mismatch")
        print(f"    expected: {list(_BRANCH_ORDER)}")
        print(f"    got:      {sorted(by_branch)}")
        raise SystemExit(1)

    rows_ord = [by_branch[b] for b in _BRANCH_ORDER]

    if not _meta_consistent(rows_ord):
        print("  FAIL: inconsistent soft/beta/gamma/t_final/dt/renorm/steps/H_N_ref across rows")
        raise SystemExit(1)

    lams = [float(r["lambda_max"]) for r in rows_ord]
    if not all(np.isfinite(lams)):
        print("  FAIL: non-finite lambda_max")
        raise SystemExit(1)
    if lams[0] <= 0.01:
        print(f"  FAIL: Newton lambda_max={lams[0]:.6g} <= 0.01")
        raise SystemExit(1)
    if any(x <= 0.02 for x in lams[1:]):
        bad = [(b, x) for b, x in zip(_BRANCH_ORDER[1:], lams[1:]) if x <= 0.02]
        print(f"  FAIL: TGP row(s) lambda_max <= 0.02: {bad}")
        raise SystemExit(1)

    Hn = float(rows_ord[0]["H_N_ref"])
    tol_h = 1e-5 * max(1.0, abs(Hn))
    for r in rows_ord:
        br = str(r["branch"])
        hm = float(r["H_model"])
        if br.endswith("H=H_N") or br.startswith("Newton"):
            if abs(hm - Hn) > tol_h:
                print(f"  FAIL: |H_model-H_N_ref| too large for {br}: {abs(hm - Hn):.6e} > {tol_h:.6e}")
                raise SystemExit(1)

    r0 = rows_ord[0]
    print(f"  csv={csv_path}")
    print(
        f"  soft={r0['soft']:.6g} beta=gamma={r0['beta']:.6g} n_quad={r0['n_quad']} "
        f"t_final={r0['t_final']} dt={r0['dt']} renorm={r0['renorm_every']} steps={r0['steps']}"
    )
    print(f"  H_N_ref={Hn:.8g}")
    print("")
    wlab, wlam = 20, 14
    print(f"  {'branch':<{wlab}} {'lambda_max':>{wlam}}")
    print(f"  {'-' * wlab} {'-' * wlam}")
    for r in rows_ord:
        print(f"  {str(r['branch']):<{wlab}} {float(r['lambda_max']):>{wlam}.8g}")

    print("")
    print("  Delta lambda_max (H=H_N minus v=0), same TGP backend:")
    pairs = (
        ("Yukawa", "Yukawa v=0", "Yukawa H=H_N"),
        ("Coulomb", "Coulomb v=0", "Coulomb H=H_N"),
        ("Pairwise", "Pairwise v=0", "Pairwise H=H_N"),
    )
    for name, a, b in pairs:
        la = float(by_branch[a]["lambda_max"])
        lb = float(by_branch[b]["lambda_max"])
        print(f"    {name}:  {lb - la:+.8g}  (v=0 {la:.6g} -> matched {lb:.6g})")

    print("")
    print("  Interpretation:")
    print("  - Rows share the same Burrau positions; v=0 vs H=H_N are different phase-space")
    print("    contexts in the TGP potentials (see ex175/ex179/ex180, ex181).")
    print("  - lambda_max here is a finite-time Benettin estimate on the leapfrog mesh,")
    print("    not an asymptotic Lyapunov exponent.")

    if not args.quick:
        _try_bar_png(rows_ord, png_path)

    print("")
    print("ex183: PASS")


if __name__ == "__main__":
    main()
