#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex120_ex115_scan_overlay.py
===========================
Wczytuje wynik **ex115** (`ex115_window_scan.csv`) i rysuje warstwę podsumowującą:

  - punkty (r_L, r_R): PASS vs FAIL (`pass_all`)
  - gwiazda **referencja** [20, 35]
  - gwiazda **najlepsze okno** (min `max_rmse_over_A` wśród PASS)

Jeśli CSV nie istnieje, uruchamia `ex115_tail_window_scan.py` (subprocess).

Wyjścia:
  _outputs/ex120_ex115_overlay.png
  _outputs/ex120_ex115_overlay_meta.csv   (best_r_L, best_r_R, n_pass, ref_pass)

Testy:
  E1: PNG + meta CSV
  E2: wiersz ref (20,35) ma pass_all
  E3: punkt „best” ma pass_all=True

Uruchomienie: python ex120_ex115_scan_overlay.py
"""

from __future__ import annotations

import csv
import io
import subprocess
import sys
import warnings
from pathlib import Path

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
else:
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

import numpy as np

warnings.filterwarnings("ignore")

HERE = Path(__file__).resolve().parent
OUT_DIR = HERE / "_outputs"
CSV115 = OUT_DIR / "ex115_window_scan.csv"
EX115_SCRIPT = HERE / "ex115_tail_window_scan.py"

REF_R_L = 20.0
REF_R_R = 35.0

RESULTS: list[tuple[str, str, str]] = []


def check(cond: bool, label: str, detail: str = "") -> bool:
    status = "PASS" if cond else "FAIL"
    RESULTS.append((label, status, detail))
    icon = "[PASS]" if cond else "[FAIL]"
    line = f"  {icon} {label}"
    if detail:
        line += f"\n         => {detail}"
    print(line)
    return cond


def ensure_ex115_csv() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    if CSV115.is_file() and CSV115.stat().st_size > 50:
        return
    print("  Brak ex115_window_scan.csv — uruchamiam ex115...")
    subprocess.run(
        [sys.executable, str(EX115_SCRIPT)],
        cwd=str(HERE),
        check=True,
    )


def load_rows() -> list[dict]:
    rows = []
    with open(CSV115, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            row["r_L"] = float(row["r_L"])
            row["r_R"] = float(row["r_R"])
            row["width"] = float(row["width"])
            row["skip"] = row["skip"].strip().lower() == "true"
            row["pass_all"] = row["pass_all"].strip().lower() == "true"
            if not row["skip"]:
                row["max_rmse_over_A"] = float(row["max_rmse_over_A"])
                row["max_mean_abs_g1"] = float(row["max_mean_abs_g1"])
                row["max_dphi_vs_ref_rad"] = float(row["max_dphi_vs_ref_rad"])
            rows.append(row)
    return rows


def main():
    print("=" * 70)
    print("EX120: overlay skanu ex115 (PASS + ref + best RMSE/A)")
    print("=" * 70)

    ensure_ex115_csv()
    rows = load_rows()

    valid = [r for r in rows if not r["skip"]]
    adm = [r for r in valid if r["pass_all"]]
    if not adm:
        print("  BŁĄD: brak okien PASS")
        sys.exit(1)

    best = min(adm, key=lambda x: x["max_rmse_over_A"])
    ref_row = next(
        (
            r
            for r in valid
            if abs(r["r_L"] - REF_R_L) < 1e-6 and abs(r["r_R"] - REF_R_R) < 1e-6
        ),
        None,
    )

    print(f"  Okna PASS: {len(adm)} / {len(valid)}")
    print(
        f"  Best RMSE/A: r_L={best['r_L']}, r_R={best['r_R']}, "
        f"max_rmse/A={best['max_rmse_over_A']:.4f}"
    )
    print()

    meta_path = OUT_DIR / "ex120_ex115_overlay_meta.csv"
    with open(meta_path, "w", encoding="utf-8") as f:
        f.write("key,value\n")
        f.write(f"best_r_L,{best['r_L']}\n")
        f.write(f"best_r_R,{best['r_R']}\n")
        f.write(f"best_max_rmse_over_A,{best['max_rmse_over_A']}\n")
        f.write(f"n_pass,{len(adm)}\n")
        f.write(f"ref_pass,{ref_row['pass_all'] if ref_row else 'missing'}\n")

    png_path = OUT_DIR / "ex120_ex115_overlay.png"
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 2, figsize=(11.0, 5.0))

        rL_fail = [r["r_L"] for r in valid if not r["pass_all"]]
        rR_fail = [r["r_R"] for r in valid if not r["pass_all"]]
        rL_ok = [r["r_L"] for r in adm]
        rR_ok = [r["r_R"] for r in adm]

        ax = axes[0]
        ax.scatter(rR_fail, rL_fail, c="lightcoral", s=35, alpha=0.7, label="FAIL", edgecolors="none")
        ax.scatter(rR_ok, rL_ok, c="seagreen", s=45, alpha=0.85, label="PASS", edgecolors="darkgreen", linewidths=0.3)
        ax.scatter([REF_R_R], [REF_R_L], c="cyan", s=200, marker="*", edgecolors="black", linewidths=0.6, zorder=6, label="ref [20,35]")
        ax.scatter(
            [best["r_R"]],
            [best["r_L"]],
            c="gold",
            s=220,
            marker="*",
            edgecolors="black",
            linewidths=0.7,
            zorder=7,
            label="best min(RMSE/A)",
        )
        ax.set_xlabel("r_R")
        ax.set_ylabel("r_L")
        ax.set_title("ex120: klasyfikacja okien (ex115)")
        ax.legend(loc="upper left", fontsize=8)
        ax.grid(True, alpha=0.3)

        ax = axes[1]
        rL_u = sorted(set(r["r_L"] for r in rows))
        rR_u = sorted(set(r["r_R"] for r in rows))
        X, Y = np.meshgrid(rR_u, rL_u)
        Z = np.full((len(rL_u), len(rR_u)), np.nan)
        Zrmse = np.full((len(rL_u), len(rR_u)), np.nan)
        for r in rows:
            if r["skip"]:
                continue
            i = rL_u.index(r["r_L"])
            j = rR_u.index(r["r_R"])
            Z[i, j] = 1.0 if r["pass_all"] else 0.0
            Zrmse[i, j] = r["max_rmse_over_A"]

        pcm = ax.pcolormesh(X, Y, Zrmse, shading="auto", cmap="viridis")
        plt.colorbar(pcm, ax=ax, label="max RMSE/A (triada)")
        ax.scatter([REF_R_R], [REF_R_L], c="white", s=160, marker="*", edgecolors="black", linewidths=0.6, zorder=6)
        ax.scatter([best["r_R"]], [best["r_L"]], c="gold", s=180, marker="*", edgecolors="black", linewidths=0.7, zorder=7)
        ax.set_xlabel("r_R")
        ax.set_ylabel("r_L")
        ax.set_title("Heatmap RMSE/A + ref / best")

        fig.suptitle("ex120 — integracja wizualna ze skanem ex115", fontsize=10)
        plt.tight_layout()
        fig.savefig(png_path, dpi=140, bbox_inches="tight")
        plt.close(fig)
        print(f"  Zapisano: {png_path}")
        print(f"  Zapisano: {meta_path}")
    except Exception as ex:
        print(f"  BŁĄD matplotlib: {ex}")
        sys.exit(1)

    print()

    check(meta_path.is_file(), "E1: meta CSV", str(meta_path))
    check(png_path.is_file(), "E1b: PNG overlay", str(png_path))
    check(ref_row is not None and ref_row["pass_all"], "E2: ref [20,35] PASS", str(ref_row))
    check(best["pass_all"], "E3: best ∈ PASS", f"r_L={best['r_L']} r_R={best['r_R']}")

    n_fail = sum(1 for _a, s, _b in RESULTS if s == "FAIL")
    print("=" * 70)
    print(f"EX120: {len(RESULTS) - n_fail}/{len(RESULTS)} testów PASS")
    print("=" * 70)
    sys.exit(1 if n_fail else 0)


if __name__ == "__main__":
    main()
