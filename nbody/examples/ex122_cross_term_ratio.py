#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex122_cross_term_ratio.py
==========================
Na triadzie (φ-FP) i wybranych odcinkach r liczy ilorazy **członu krzyżowego**
vs **V'(g)** z prawej strony ODE (ex114):

  R_CV = | (alpha/g) * g'^2 | / | V'(g) |,

oraz pomocniczo R_Ch = |cross| / |g-1| (gdy |g-1| > 1e-6).

Używa g'(r) z **integratora** (solve_ivp), nie różnicowania post-factum.

Wyjścia:
  _outputs/ex122_cross_term_ratios.csv

Testy:
  C1: CSV
  C2: mediana R_CV na [25,32] < 0.5 dla wszystkich trzech generacji
      (luźny próg — „cross nie dominuje nad V'” w środku ogona)
  C3: percentyl 99 R_CV na [20,35] < 45 (max jest zawyżone przy |V'|→0)

Teoria: ex122_tail_cross_term_sketch.md
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

warnings.filterwarnings("ignore")

OUT_DIR = Path(__file__).resolve().parent / "_outputs"

BANDS = [
    ("win_20_35", 20.0, 35.0),
    ("win_25_32", 25.0, 32.0),
]

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
    spec = importlib.util.spec_from_file_location("ex122_ex114", path)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


def ratios_on_band(
    r: np.ndarray,
    g: np.ndarray,
    gp: np.ndarray,
    alpha: float,
    r_lo: float,
    r_hi: float,
) -> tuple[float, float, float, float, int]:
    """Zwraca median_R_CV, max_R_CV, p99_R_CV, median_R_Ch, n_pts."""
    mask = (r >= r_lo) & (r <= r_hi)
    Vp = g[mask] ** 2 * (1.0 - g[mask])
    cross = (alpha / np.maximum(g[mask], 1e-10)) * gp[mask] ** 2
    h = np.abs(g[mask] - 1.0)

    R_CV = np.abs(cross) / np.maximum(np.abs(Vp), 1e-14)
    R_Ch = np.abs(cross) / np.maximum(h, 1e-6)

    finite = np.isfinite(R_CV) & (np.abs(Vp) > 1e-12)
    if np.sum(finite) < 5:
        return float("nan"), float("nan"), float("nan"), float("nan"), 0

    R_CV = R_CV[finite]
    R_Ch = R_Ch[finite]
    return (
        float(np.median(R_CV)),
        float(np.max(R_CV)),
        float(np.percentile(R_CV, 99.0)),
        float(np.median(R_Ch)),
        int(np.sum(finite)),
    )


def main():
    print("=" * 70)
    print("EX122: |cross| / |V'|  oraz  |cross| / |g-1|  (triada)")
    print("=" * 70)

    ex114 = load_ex114()
    ALPHA = ex114.ALPHA
    find_g0_star = ex114.find_g0_star
    integrate_soliton = ex114.integrate_soliton
    PHI = ex114.PHI

    g0_star = find_g0_star()
    if g0_star is None:
        print("  BŁĄD: brak g₀*")
        sys.exit(1)

    triad = [
        ("e", float(g0_star)),
        ("mu", float(PHI * g0_star)),
        ("tau", float(PHI**2 * g0_star)),
    ]

    rows = []
    med_narrow: list[float] = []
    p99_wide: list[float] = []

    for name, g0 in triad:
        r, g, gp = integrate_soliton(g0)
        print(f"  {name}: g₀={g0:.5f}")
        for band_name, r_lo, r_hi in BANDS:
            med_cv, max_cv, p99_cv, med_ch, npt = ratios_on_band(
                r, g, gp, ALPHA, r_lo, r_hi
            )
            rows.append(
                {
                    "gen": name,
                    "band": band_name,
                    "median_R_cross_over_Vp": med_cv,
                    "max_R_cross_over_Vp": max_cv,
                    "p99_R_cross_over_Vp": p99_cv,
                    "median_R_cross_over_h": med_ch,
                    "n_points": npt,
                }
            )
            print(
                f"      [{band_name}]  med|cross/V'|={med_cv:.4f}  max={max_cv:.4f}  "
                f"p99={p99_cv:.4f}  n={npt}"
            )
            if band_name == "win_25_32" and np.isfinite(med_cv):
                med_narrow.append(med_cv)
            if band_name == "win_20_35" and np.isfinite(p99_cv):
                p99_wide.append(p99_cv)
        print()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    csv_path = OUT_DIR / "ex122_cross_term_ratios.csv"
    keys = list(rows[0].keys())
    with open(csv_path, "w", encoding="utf-8") as f:
        f.write(",".join(keys) + "\n")
        for row in rows:
            f.write(",".join(str(row[k]) for k in keys) + "\n")
    print(f"  Zapisano: {csv_path}")
    print()

    max_med_narrow = max(med_narrow) if med_narrow else float("nan")
    max_p99_wide = max(p99_wide) if p99_wide else float("nan")

    check(csv_path.is_file(), "C1: CSV", str(csv_path))
    check(
        len(med_narrow) == 3 and max_med_narrow < 0.5,
        "C2: mediana |cross/V'| < 0.5 na [25,32] (triada)",
        f"max mediana={max_med_narrow:.4f}",
    )
    check(
        len(p99_wide) == 3 and max_p99_wide < 45.0,
        "C3: p99 |cross/V'| < 45 na [20,35] (triada)",
        f"max p99={max_p99_wide:.4f}",
    )

    n_fail = sum(1 for _a, s, _b in RESULTS if s == "FAIL")
    print("=" * 70)
    print(f"EX122: {len(RESULTS) - n_fail}/{len(RESULTS)} testów PASS")
    print("=" * 70)
    sys.exit(1 if n_fail else 0)


if __name__ == "__main__":
    main()
