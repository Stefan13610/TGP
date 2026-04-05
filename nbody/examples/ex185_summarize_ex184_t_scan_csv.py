#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex185_summarize_ex184_t_scan_csv.py
===================================
P1.C (cd.): **odczyt i walidacja** CSV z ``ex184`` (wiele ``t_final``, kazdy po siedem
wierszy jak ``ex181``). Sprawdza spojnosc metadanych miedzy **wszystkimi** wierszami,
progi jak ``ex183``, drukuje zwięzłą tabele vs ``t_final``. Opcjonalnie PNG: ``lambda_max``
vs ``t_final`` dla siedmiu galezi (``matplotlib``).

Domyslne wejscie: ``examples/_outputs/ex184_seven_row_t_final_scan.csv``
(najpierw ``ex184`` lub ``verify_nbody_lyapunov_quick``).

Analog: ``ex183`` wzgledem ``ex182``. Czterowierszowy skan matched ``H``: ``ex189`` wzgledem ``ex188``.

``--quick``: stdout + walidacja (bez matplotlib).

Uruchomienie:
``python ex185_summarize_ex184_t_scan_csv.py [--quick] [--csv PATH] [--png PATH]``
"""

from __future__ import annotations

import argparse
import csv
import os
from collections import defaultdict

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_DEFAULT_CSV = os.path.join(_HERE, "_outputs", "ex184_seven_row_t_final_scan.csv")
_DEFAULT_PNG = os.path.join(_HERE, "_outputs", "ex185_ex184_lambda_vs_t.png")

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


def _meta_consistent_slice(rows: list[dict[str, float | int | str]]) -> bool:
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


def _global_invariants(rows: list[dict[str, float | int | str]]) -> bool:
    r0 = rows[0]
    keys = ("soft", "beta", "gamma", "n_quad", "dt", "renorm_every", "H_N_ref")
    for r in rows[1:]:
        if r["n_quad"] != r0["n_quad"] or r["renorm_every"] != r0["renorm_every"]:
            return False
        if not np.allclose(
            [r["soft"], r["beta"], r["gamma"], r["dt"], r["H_N_ref"]],
            [r0["soft"], r0["beta"], r0["gamma"], r0["dt"], r0["H_N_ref"]],
            rtol=0.0,
            atol=1e-12,
        ):
            return False
    return True


def _validate_ordered_slice(rows_ord: list[dict[str, float | int | str]]) -> str | None:
    if not _meta_consistent_slice(rows_ord):
        return "inconsistent meta within one t_final slice"
    lams = [float(r["lambda_max"]) for r in rows_ord]
    if not all(np.isfinite(lams)):
        return "non-finite lambda_max"
    if lams[0] <= 0.01:
        return f"Newton lambda_max={lams[0]:.6g} <= 0.01"
    if any(x <= 0.02 for x in lams[1:]):
        return f"TGP lambda_max <= 0.02 in slice t_final={rows_ord[0]['t_final']}"
    Hn = float(rows_ord[0]["H_N_ref"])
    tol_h = 1e-5 * max(1.0, abs(Hn))
    for r in rows_ord:
        br = str(r["branch"])
        hm = float(r["H_model"])
        if br.endswith("H=H_N") or br.startswith("Newton"):
            if abs(hm - Hn) > tol_h:
                return f"|H_model-H_N_ref| for {br}"
    return None


def _try_lines_png(
    t_sorted: list[float],
    by_t_branch: dict[float, dict[str, float]],
    png_path: str,
) -> bool:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as e:
        print(f"  (matplotlib not available — PNG skipped: {e})")
        return False

    fig, ax = plt.subplots(figsize=(8.5, 4.8), constrained_layout=True)
    for br in _BRANCH_ORDER:
        ys = [by_t_branch[tf][br] for tf in t_sorted]
        ax.plot(t_sorted, ys, "o-", linewidth=1.2, markersize=4, label=br)
    ax.set_xlabel("t_final")
    ax.set_ylabel("lambda_max")
    ax.set_title("ex184 scan: seven branches vs horizon")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=7, loc="best")
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

    print("ex185: read ex184 t-scan CSV + summary (+ optional PNG)")
    if not os.path.isfile(csv_path):
        print(f"  FAIL: missing {csv_path}")
        print("  Run: python nbody/examples/ex184_lyapunov_seven_row_t_final_scan_csv.py --quick")
        raise SystemExit(1)

    try:
        rows = _load_rows(csv_path)
    except ValueError as e:
        print(f"  FAIL: {e}")
        raise SystemExit(1)

    if len(rows) < 14 or len(rows) % 7 != 0:
        print(f"  FAIL: expected at least 14 rows and multiple of 7, got {len(rows)}")
        raise SystemExit(1)

    if not _global_invariants(rows):
        print("  FAIL: soft/beta/gamma/n_quad/dt/renorm/H_N_ref differ across CSV rows")
        raise SystemExit(1)

    by_t: dict[float, list[dict[str, float | int | str]]] = defaultdict(list)
    for r in rows:
        by_t[float(r["t_final"])].append(r)

    if len(by_t) < 2:
        print(f"  FAIL: expected at least 2 distinct t_final, got {len(by_t)}")
        raise SystemExit(1)

    t_sorted = sorted(by_t)
    slices: list[tuple[float, list[dict[str, float | int | str]]]] = []
    for tf in t_sorted:
        grp = by_t[tf]
        if len(grp) != 7:
            print(f"  FAIL: t_final={tf:g} has {len(grp)} rows, expected 7")
            raise SystemExit(1)
        by_branch = {str(r["branch"]): r for r in grp}
        if set(by_branch) != set(_BRANCH_ORDER):
            print(f"  FAIL: branch set mismatch at t_final={tf:g}")
            raise SystemExit(1)
        rows_ord = [by_branch[b] for b in _BRANCH_ORDER]
        err = _validate_ordered_slice(rows_ord)
        if err:
            print(f"  FAIL at t_final={tf:g}: {err}")
            raise SystemExit(1)
        slices.append((tf, rows_ord))

    r0 = slices[0][1][0]
    print(f"  csv={csv_path}")
    print(
        f"  soft={r0['soft']:.6g} beta=gamma={r0['beta']:.6g} n_quad={r0['n_quad']} "
        f"dt={r0['dt']} renorm={r0['renorm_every']}  distinct t_final={len(t_sorted)}"
    )
    print(f"  H_N_ref={float(r0['H_N_ref']):.8g}")
    print("")
    wt, ws, wn, wr = 10, 6, 12, 18
    print(
        f"  {'t_final':>{wt}} {'steps':>{ws}} {'Newton':>{wn}} "
        f"{'TGP lambda range':>{wr}}"
    )
    print(f"  {'-' * wt} {'-' * ws} {'-' * wn} {'-' * wr}")

    by_t_branch: dict[float, dict[str, float]] = {}
    for tf, rows_ord in slices:
        lams = [float(r["lambda_max"]) for r in rows_ord]
        st = int(rows_ord[0]["steps"])
        newt = lams[0]
        tgp = lams[1:]
        by_t_branch[tf] = {b: float(rows_ord[i]["lambda_max"]) for i, b in enumerate(_BRANCH_ORDER)}
        print(
            f"  {tf:>{wt}.5g} {st:>{ws}d} {newt:>{wn}.5g} "
            f"[{min(tgp):.5g}, {max(tgp):.5g}]"
        )

    tf_last, last_ord = slices[-1]
    by_b = {str(r["branch"]): r for r in last_ord}
    print("")
    print(f"  Delta lambda_max at t_final={tf_last:g} (H=H_N minus v=0):")
    for name, a, b in (
        ("Yukawa", "Yukawa v=0", "Yukawa H=H_N"),
        ("Coulomb", "Coulomb v=0", "Coulomb H=H_N"),
        ("Pairwise", "Pairwise v=0", "Pairwise H=H_N"),
    ):
        la = float(by_b[a]["lambda_max"])
        lb = float(by_b[b]["lambda_max"])
        print(f"    {name}:  {lb - la:+.8g}")

    print("")
    print("  Interpretation:")
    print("  - lambda_max vs t_final is a finite-time Benettin estimate; slopes need not")
    print("    stabilize on this Burrau + leapfrog grid (see ex166/ex168, PLAN P1.C).")
    print("  - Comparing curves across branches separates Newton from TGP contexts.")

    if not args.quick:
        _try_lines_png(t_sorted, by_t_branch, png_path)

    print("")
    print("ex185: PASS")


if __name__ == "__main__":
    main()
