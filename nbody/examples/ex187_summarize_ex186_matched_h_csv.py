#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex187_summarize_ex186_matched_h_csv.py
======================================
P1 (cd.): **odczyt i podsumowanie** CSV z ``ex186`` (cztery wiersze ``lambda_max``:
Newton + trzy TGP przy ``H`` dopasowanym do ``H_N``, jak ``ex178``). Te same kolumny
co ``ex182``/``ex186``. Sprawdza spojnosc metadanych, progi jak ``ex178``/``ex186``,
wszystkie ``H_model`` wzgledem ``H_N_ref`` (matched family). Opcjonalnie PNG (slupki).

Domyslne wejscie: ``examples/_outputs/ex186_matched_h_family.csv``
(najpierw ``ex186`` lub cala regresja ``verify_nbody_lyapunov_quick``).

Analog roli do ``ex183`` wzgledem ``ex182``. Skan ``t_final``: ``ex188``; odczyt skanu: ``ex189``. Leapfrog ``|dE/E0|``: ``ex190``; RK45: ``ex191``; tabela: ``ex192``; CSV LF/RK: ``ex193``; odczyt: ``ex194``.

``--quick``: stdout + walidacja (bez matplotlib).

Uruchomienie:
``python ex187_summarize_ex186_matched_h_csv.py [--quick] [--csv PATH] [--png PATH]``
"""

from __future__ import annotations

import argparse
import csv
import os

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_DEFAULT_CSV = os.path.join(_HERE, "_outputs", "ex186_matched_h_family.csv")
_DEFAULT_PNG = os.path.join(_HERE, "_outputs", "ex187_ex186_lambda_bar.png")

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
    fig, ax = plt.subplots(figsize=(8.2, 4.0), constrained_layout=True)
    x = np.arange(len(labels))
    ax.bar(x, y, color="darkseagreen", alpha=0.88)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=28, ha="right")
    ax.set_ylabel("lambda_max")
    ax.set_title("ex186 matched-H family (ex178, four branches)")
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

    print("ex187: read ex186 matched-H CSV + summary (+ optional PNG)")
    if not os.path.isfile(csv_path):
        print(f"  FAIL: missing {csv_path}")
        print("  Run: python nbody/examples/ex186_lyapunov_matched_h_family_csv.py --quick")
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
        print("  FAIL: inconsistent soft/beta/gamma/t_final/dt/renorm/steps/H_N_ref")
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
        hm = float(r["H_model"])
        if abs(hm - Hn) > tol_h:
            print(
                f"  FAIL: |H_model-H_N_ref| for {r['branch']}: "
                f"{abs(hm - Hn):.6e} > {tol_h:.6e}"
            )
            raise SystemExit(1)

    r0 = rows_ord[0]
    print(f"  csv={csv_path}")
    print(
        f"  soft={r0['soft']:.6g} beta=gamma={r0['beta']:.6g} n_quad={r0['n_quad']} "
        f"t_final={r0['t_final']} dt={r0['dt']} renorm={r0['renorm_every']} steps={r0['steps']}"
    )
    print(f"  H_N_ref={Hn:.8g} (all rows matched within tol)")
    print("")
    wlab, wlam = 18, 14
    print(f"  {'branch':<{wlab}} {'lambda_max':>{wlam}}")
    print(f"  {'-' * wlab} {'-' * wlam}")
    for r in rows_ord:
        print(f"  {str(r['branch']):<{wlab}} {float(r['lambda_max']):>{wlam}.8g}")

    print("")
    print("  Interpretation:")
    print("  - All four trajectories share the same target Hamiltonian value H_N in")
    print("    their respective potentials (see ex178, ex174/ex176/ex177).")
    print("  - For v=0 vs H=H_N splits per backend, see ex181/ex182 and ex183.")
    print("  - lambda_max is a finite-time Benettin estimate on the leapfrog mesh.")

    if not args.quick:
        _try_bar_png(rows_ord, png_path)

    print("")
    print("ex187: PASS")


if __name__ == "__main__":
    main()
