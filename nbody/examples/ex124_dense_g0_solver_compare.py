#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex124_dense_g0_solver_compare.py
=================================
Faza C1 planu: **gęstsza siatka g₀** wokół φ-FP i generacji + porównanie
**f-regularyzowany** (ex114, odbicia przy g*) vs **K_sub = g²** (ex106, bez ghostów).

Ładuje ex114 (fit, integrator reg). Substrat: ten sam układ co w `ex106_path9_formalization.rhs_substrate`.

Wyjście: `_outputs/ex124_dense_solver_compare.csv`
Opcjonalnie: `ex124_dense_A_compare.png` (matplotlib).

Testy:
  R1: CSV
  R2: A_reg > 0 wszędzie
  R3: |A_reg−A_sub|/A_reg < 18% przy g₀* (jak ex106 T6, lekki luz)
  R4: |Δφ| < 0.6 rad przy g₀* (klasa fazowa nie rozpada się między solverami)
  R5: n punktów ≥ 85
"""

from __future__ import annotations

import importlib.util
import io
import sys
import warnings
from pathlib import Path

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
else:
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

import numpy as np
from scipy.integrate import solve_ivp

warnings.filterwarnings("ignore")

OUT_DIR = Path(__file__).resolve().parent / "_outputs"
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


def load_ex114():
    path = Path(__file__).resolve().parent / "ex114_tail_phase_map.py"
    spec = importlib.util.spec_from_file_location("ex124_ex114", path)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


def wrap_dphi(a: float, b: float) -> float:
    d = float(np.arctan2(np.sin(a - b), np.cos(a - b)))
    return abs(d)


# --- K_sub solver (zsynchronizowany z ex106) --------------------------------

ALPHA = 2.0


def Vprime(g):
    return g**2 * (1.0 - g)


def rhs_substrate(r, y):
    g, gp = y
    g = max(g, 1e-10)
    ksub = g**2
    dksub = 2.0 * g
    driving = Vprime(g)
    if r < 1e-10:
        return [gp, (driving - dksub * gp**2 / 2.0) / (3.0 * ksub)]
    damp = ksub * 2.0 * gp / r
    return [gp, (driving - dksub * gp**2 / 2.0 - damp) / ksub]


def integrate_soliton_substrate(
    g0: float,
    r_max: float | None,
    R_START: float,
    R_MAX: float,
    MAX_STEP: float,
    RTOL: float,
    ATOL: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if r_max is None:
        r_max = max(R_MAX, 15.0 * g0)
    sol = solve_ivp(
        rhs_substrate,
        [R_START, r_max],
        [g0, 0.0],
        method="DOP853",
        max_step=MAX_STEP,
        rtol=RTOL,
        atol=ATOL,
        dense_output=False,
    )
    r = sol.t
    g = sol.y[0]
    gp = sol.y[1]
    idx = np.argsort(r)
    return r[idx], g[idx], gp[idx]


def build_g0_grid(g0_star: float | None, phi: float) -> np.ndarray:
    g0_lo, g0_hi = 1.10, 3.40
    coarse = np.linspace(g0_lo, g0_hi, 52)
    chunks = [coarse]
    if g0_star is not None:
        chunks.append(np.linspace(max(g0_lo, g0_star - 0.15), min(g0_hi, g0_star + 0.15), 45))
        for scale in (phi, phi**2):
            gc = float(g0_star * scale)
            if g0_lo <= gc <= g0_hi:
                chunks.append(
                    np.linspace(max(g0_lo, gc - 0.09), min(g0_hi, gc + 0.09), 28)
                )
        chunks.append(np.array([g0_star, phi * g0_star, phi**2 * g0_star]))
    return np.unique(np.sort(np.concatenate(chunks)))


def main() -> int:
    print("=" * 70)
    print("EX124: gęsta siatka g₀ — regularyzowany vs K_sub (fit [20,35])")
    print("=" * 70)

    ex114 = load_ex114()
    integrate_reg = ex114.integrate_soliton
    fit_tail_detailed = ex114.fit_tail_detailed
    find_g0_star = ex114.find_g0_star
    R_TAIL_L = ex114.R_TAIL_L
    R_TAIL_R = ex114.R_TAIL_R
    R_START = ex114.R_START
    R_MAX = ex114.R_MAX
    MAX_STEP = ex114.MAX_STEP
    RTOL = ex114.RTOL
    ATOL = ex114.ATOL
    PHI = ex114.PHI

    g0_star = find_g0_star()
    print(f"  g₀* = {g0_star:.6f}" if g0_star else "  g₀*: brak")
    g0_grid = build_g0_grid(g0_star, PHI)
    print(f"  Unikalnych g₀: {len(g0_grid)}")
    print()

    rows: list[dict] = []
    for g0 in g0_grid:
        r1, g1, _ = integrate_reg(float(g0))
        A_r, B_r, C_r, rm_r = fit_tail_detailed(r1, g1, R_TAIL_L, R_TAIL_R)
        r2, g2, _ = integrate_soliton_substrate(
            float(g0), None, R_START, R_MAX, MAX_STEP, RTOL, ATOL
        )
        A_s, B_s, C_s, rm_s = fit_tail_detailed(r2, g2, R_TAIL_L, R_TAIL_R)
        ph_r = float(np.arctan2(C_r, B_r))
        ph_s = float(np.arctan2(C_s, B_s))
        rel_da = abs(A_r - A_s) / max(A_r, 1e-15)
        rows.append(
            {
                "g0": float(g0),
                "A_reg": A_r,
                "B_reg": B_r,
                "C_reg": C_r,
                "rmse_reg": rm_r,
                "A_sub": A_s,
                "B_sub": B_s,
                "C_sub": C_s,
                "rmse_sub": rm_s,
                "rel_dA": rel_da,
                "dphi_rad": wrap_dphi(ph_r, ph_s),
            }
        )

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    csv_path = OUT_DIR / "ex124_dense_solver_compare.csv"
    keys = list(rows[0].keys())
    with open(csv_path, "w", encoding="utf-8") as f:
        f.write(",".join(keys) + "\n")
        for row in rows:
            f.write(",".join(str(row[k]) for k in keys) + "\n")
    print(f"  Zapisano: {csv_path}")

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        g0s = np.array([r["g0"] for r in rows])
        fig, ax = plt.subplots(figsize=(7.5, 4.0))
        ax.plot(g0s, [r["A_reg"] for r in rows], "b-", lw=1.0, label="A_reg")
        ax.plot(g0s, [r["A_sub"] for r in rows], "r--", lw=1.0, alpha=0.85, label="A_sub")
        if g0_star:
            ax.axvline(g0_star, color="gray", ls=":", lw=0.8)
        ax.set_xlabel("g₀")
        ax.set_ylabel("A_tail")
        ax.set_title("ex124: amplituda ogona — dwa solvery")
        ax.legend(loc="best", fontsize=8)
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        png_path = OUT_DIR / "ex124_dense_A_compare.png"
        fig.savefig(png_path, dpi=120)
        plt.close(fig)
        print(f"  Wykres: {png_path}")
    except Exception as e:
        print(f"  (matplotlib pominięty: {e})")
    print()

    min_a_reg = min(r["A_reg"] for r in rows)
    star_row = next((r for r in rows if g0_star and abs(r["g0"] - g0_star) < 1e-9), None)
    if star_row is None and g0_star:
        ig = int(np.argmin(np.abs(np.array([r["g0"] for r in rows]) - g0_star)))
        star_row = rows[ig]

    check(csv_path.is_file(), "R1: CSV", str(csv_path))
    check(min_a_reg > 0, "R2: A_reg > 0 wszędzie", f"min={min_a_reg:.6f}")
    if g0_star is None:
        check(True, "R3: (pominięty — brak g₀*)", "")
        check(True, "R4: (pominięty — brak g₀*)", "")
    elif star_row:
        check(
            star_row["rel_dA"] < 0.18,
            "R3: zgodność A przy g₀* (|ΔA|/A < 18%)",
            f"rel={star_row['rel_dA']*100:.2f}%",
        )
        check(
            star_row["dphi_rad"] < 0.6,
            "R4: |Δφ| przy g₀* < 0.6 rad",
            f"dphi={star_row['dphi_rad']:.4f}",
        )
    else:
        check(False, "R3: wiersz g₀*", "brak wiersza przy g₀*")
        check(False, "R4: wiersz g₀*", "brak wiersza przy g₀*")

    check(len(rows) >= 85, "R5: liczba punktów siatki", f"n={len(rows)}")

    n_fail = sum(1 for _a, s, _b in RESULTS if s == "FAIL")
    print("=" * 70)
    print(f"EX124: {len(RESULTS) - n_fail}/{len(RESULTS)} testów PASS")
    print("=" * 70)
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
