#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex118_rmse_vs_linear_residual.py
================================
Ta sama triada (φ-FP) i pasma co ex117: porównanie

  zeta = RMS(L[h]) / RMS(h)     (L = h'' + 2h'/r + h)
  rho  = RMSE / A               (fit (g-1)r ~ B cos r + C sin r)

Oczekiwanie: większe odchylenie od równania liniowego → zwykle większa reszta
ansatzu cos/sin (dodatnia korelacja na tym zbiorze 3×2 punktów).

Wyjścia:
  _outputs/ex118_rmse_vs_zeta.csv
  (w konsoli: współczynnik Pearsona rho(zeta, rho))

Testy:
  K1: CSV
  K2: wszystkie rho = RMSE/A < 0.12 w [25,32] (triada)
  K3: Pearson(zeta, rho) > 0.25 (słaba dodatnia korelacja na 6 punktach)

Uruchomienie: python ex118_rmse_vs_linear_residual.py
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


def load_module(name: str, fname: str):
    path = Path(__file__).resolve().parent / fname
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


def main():
    print("=" * 70)
    print("EX118: RMSE/A vs zeta (reszta operatora liniowego)")
    print("=" * 70)

    ex114 = load_module("ex118_ex114", "ex114_tail_phase_map.py")
    ex117 = load_module("ex118_ex117", "ex117_linear_operator_residual.py")
    zeta_ratio = ex117.zeta_ratio
    integrate_soliton = ex114.integrate_soliton
    fit_tail_detailed = ex114.fit_tail_detailed

    g0_star = ex114.find_g0_star()
    if g0_star is None:
        print("  BŁĄD: brak g₀*")
        sys.exit(1)

    PHI = ex114.PHI
    triad = [
        ("e", float(g0_star)),
        ("mu", float(PHI * g0_star)),
        ("tau", float(PHI**2 * g0_star)),
    ]

    rows = []
    zetas: list[float] = []
    rhos: list[float] = []

    for name, g0 in triad:
        r, g, _ = integrate_soliton(g0)
        print(f"  {name}: g₀={g0:.5f}")
        for band_name, r_lo, r_hi in BANDS:
            A, Bf, Cf, rmse = fit_tail_detailed(r, g, r_lo, r_hi)
            rho_fit = rmse / max(A, 1e-15)
            z, rms_h = zeta_ratio(r, g, r_lo, r_hi)
            rows.append(
                {
                    "gen": name,
                    "band": band_name,
                    "r_lo": r_lo,
                    "r_hi": r_hi,
                    "A_tail": A,
                    "rmse_tail": rmse,
                    "rmse_over_A": rho_fit,
                    "zeta": z,
                    "rms_h": rms_h,
                }
            )
            if np.isfinite(z) and np.isfinite(rho_fit):
                zetas.append(z)
                rhos.append(rho_fit)
            print(
                f"      [{band_name}]  zeta={z:.4f}  RMSE/A={rho_fit:.4f}  "
                f"(A={A:.4f})"
            )
        print()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    csv_path = OUT_DIR / "ex118_rmse_vs_zeta.csv"
    keys = list(rows[0].keys())
    with open(csv_path, "w", encoding="utf-8") as f:
        f.write(",".join(keys) + "\n")
        for row in rows:
            f.write(",".join(str(row[k]) for k in keys) + "\n")
    print(f"  Zapisano: {csv_path}")

    zv = np.array(zetas, dtype=float)
    rv = np.array(rhos, dtype=float)
    if len(zv) >= 3 and np.std(zv) > 1e-12 and np.std(rv) > 1e-12:
        pearson = float(np.corrcoef(zv, rv)[0, 1])
    else:
        pearson = float("nan")

    print(f"  Pearson(zeta, RMSE/A) na n={len(zv)} punktach: {pearson:.4f}")
    print()

    rho_narrow = [
        rw["rmse_over_A"]
        for rw in rows
        if rw["band"] == "win_25_32" and np.isfinite(rw["rmse_over_A"])
    ]
    max_rho_narrow = max(rho_narrow) if rho_narrow else float("nan")

    check(csv_path.is_file(), "K1: CSV", str(csv_path))
    check(
        np.isfinite(max_rho_narrow) and max_rho_narrow < 0.12,
        "K2: max(RMSE/A) triady na [25,32] < 0.12",
        f"max={max_rho_narrow:.4f}",
    )
    check(
        np.isfinite(pearson) and pearson > 0.25,
        "K3: Pearson(zeta, RMSE/A) > 0.25",
        f"r={pearson:.4f}",
    )

    n_fail = sum(1 for _a, s, _b in RESULTS if s == "FAIL")
    print("=" * 70)
    print(f"EX118: {len(RESULTS) - n_fail}/{len(RESULTS)} testów PASS")
    print("=" * 70)
    sys.exit(1 if n_fail else 0)


if __name__ == "__main__":
    main()
