#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex116_tail_rL_from_epsilon.py
==============================
Dynamiczny lewy brzeg okna tail: **najmniejsze r_L** (jak najdalej w lewo),
przy którym dla **wszystkich trzech** profili triady jednocześnie

    średnia |g(r)-1| na [r_L, r_R] < ε

Przy ustalonym r_R (domyślnie 35, jak w ex114). To realizuje intencję
ex114_tail_okna_podsumowanie.md §4.2 — okno zdefiniowane przez dynamikę
(ε), a nie ręcznie.

Wyjścia:
  _outputs/ex116_rL_vs_epsilon.csv

Testy:
  M1: monotoniczność — przy mniejszym ε potrzebne większe r_L (ostrzejszy próg)
  M2: dla ε=0.12 istnieje rozwiązanie z szerokością ≥ 8
  M3: r_L(0.12) ≤ 20 + tolerancja (okno [20,35] jest **co najmniej** tak ostre)

Uruchomienie: python ex116_tail_rL_from_epsilon.py  (z examples/)
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

R_R_FIXED = 35.0
MIN_WIDTH = 8.0
MIN_POINTS = 10
R_L_SEARCH_LO = 4.0
R_L_SEARCH_HI = R_R_FIXED - MIN_WIDTH  # 27

EPS_LIST = [0.12, 0.10, 0.08, 0.06, 0.05, 0.04]

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
    spec = importlib.util.spec_from_file_location("ex116_ex114", path)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


def mean_abs_g_minus_one(r: np.ndarray, g: np.ndarray, r_L: float, r_R: float) -> float:
    mask = (r >= r_L) & (r <= r_R)
    if np.sum(mask) < MIN_POINTS:
        return float("inf")
    return float(np.mean(np.abs(g[mask] - 1.0)))


def triad_all_below_epsilon(
    profiles: list[tuple[np.ndarray, np.ndarray]],
    r_L: float,
    r_R: float,
    eps: float,
) -> bool:
    for r, g in profiles:
        m = mean_abs_g_minus_one(r, g, r_L, r_R)
        if m >= eps or not np.isfinite(m):
            return False
    return True


def leftmost_r_L(profiles: list[tuple[np.ndarray, np.ndarray]], r_R: float, eps: float) -> float | None:
    """
    Najmniejsze r_L w [R_L_SEARCH_LO, R_R - MIN_WIDTH], dla którego triada
    spełnia średnią < ε na [r_L, r_R].
    """
    lo, hi = float(R_L_SEARCH_LO), float(r_R - MIN_WIDTH)

    def ok(r_L: float) -> bool:
        return triad_all_below_epsilon(profiles, r_L, r_R, eps)

    if not ok(hi):
        return None
    if ok(lo):
        return lo

    for _ in range(50):
        mid = 0.5 * (lo + hi)
        if ok(mid):
            hi = mid
        else:
            lo = mid
        if hi - lo < 1e-4:
            break
    return hi


def main():
    print("=" * 70)
    print("EX116: r_L z progu ε (śr.|g-1| na [r_L, r_R], triada)")
    print("=" * 70)
    print(f"  r_R = {R_R_FIXED},  min szerokość = {MIN_WIDTH},  ε ∈ {EPS_LIST}")
    print()

    ex114 = load_ex114()
    g0_star = ex114.find_g0_star()
    if g0_star is None:
        print("  BŁĄD: brak g₀*")
        sys.exit(1)

    PHI = ex114.PHI
    integrate_soliton = ex114.integrate_soliton
    triad_g0 = [float(g0_star), float(PHI * g0_star), float(PHI**2 * g0_star)]
    names = ["e", "mu", "tau"]

    profiles: list[tuple[np.ndarray, np.ndarray]] = []
    for name, g0 in zip(names, triad_g0):
        r, g, _ = integrate_soliton(g0)
        profiles.append((r, g))
        print(f"    profil {name}: g₀={g0:.5f}")

    rows = []
    for eps in EPS_LIST:
        r_L_star = leftmost_r_L(profiles, R_R_FIXED, eps)
        if r_L_star is None:
            means = [float("nan")] * 3
            rows.append(
                {
                    "epsilon": eps,
                    "r_L_star": float("nan"),
                    "width": float("nan"),
                    "mean_e": float("nan"),
                    "mean_mu": float("nan"),
                    "mean_tau": float("nan"),
                }
            )
            print(f"  ε={eps:.3f}: brak r_L w [{R_L_SEARCH_LO}, {R_R_FIXED - MIN_WIDTH}]")
            continue

        means = [
            mean_abs_g_minus_one(r, g, r_L_star, R_R_FIXED) for r, g in profiles
        ]
        rows.append(
            {
                "epsilon": eps,
                "r_L_star": r_L_star,
                "width": R_R_FIXED - r_L_star,
                "mean_e": means[0],
                "mean_mu": means[1],
                "mean_tau": means[2],
            }
        )
        print(
            f"  ε={eps:.3f}:  r_L*={r_L_star:.4f}  width={R_R_FIXED - r_L_star:.2f}  "
            f"max(śr.|g-1|)={max(means):.5f}"
        )

    print()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    csv_path = OUT_DIR / "ex116_rL_vs_epsilon.csv"
    keys = list(rows[0].keys())
    with open(csv_path, "w", encoding="utf-8") as f:
        f.write(",".join(keys) + "\n")
        for row in rows:
            f.write(",".join(str(row[k]) for k in keys) + "\n")
    print(f"  Zapisano: {csv_path}")
    print()

    # --- monotoniczność: ε maleje => r_L* nie maleje (ostrzejszy próg nie wcześniej) ---
    r_Ls = [rw["r_L_star"] for rw in rows if np.isfinite(rw["r_L_star"])]
    epss = [rw["epsilon"] for rw in rows if np.isfinite(rw["r_L_star"])]
    mono_ok = True
    for i in range(len(r_Ls) - 1):
        # epss są malejące [0.12, 0.10, ...] — przy mniejszym ε oczekujemy r_L* >= poprzedni
        if r_Ls[i + 1] + 1e-6 < r_Ls[i]:
            mono_ok = False
            break

    row_012 = next((rw for rw in rows if abs(rw["epsilon"] - 0.12) < 1e-9), None)
    r_L_012 = row_012["r_L_star"] if row_012 else float("nan")

    check(csv_path.is_file(), "M1: CSV zapisany", str(csv_path))
    check(
        mono_ok,
        "M2: r_L* rośnie (słabo) gdy ε maleje wzdłuż listy",
        f"eps={epss}, r_L={['%.3f'%x for x in r_Ls]}",
    )
    check(
        row_012 is not None and np.isfinite(r_L_012) and (R_R_FIXED - r_L_012) >= MIN_WIDTH,
        "M3: ε=0.12 daje szerokie okno",
        f"r_L*={r_L_012}, width={R_R_FIXED - r_L_012 if np.isfinite(r_L_012) else float('nan')}",
    )
    # ręczne [20,35] jest „dopuszczalne” jeśli dynamiczne r_L*(0.12) ≤ 20 + slack
    check(
        np.isfinite(r_L_012) and r_L_012 <= 20.5,
        "M4: r_L*(0.12) ≤ 20.5 — protokół [20,35] nie zaczyna się wcześniej niż dynamika",
        f"r_L*={r_L_012:.4f} (ref r_L=20)",
    )

    n_fail = sum(1 for _a, s, _b in RESULTS if s == "FAIL")
    print("=" * 70)
    print(f"EX116: {len(RESULTS) - n_fail}/{len(RESULTS)} testów PASS")
    print("=" * 70)
    sys.exit(1 if n_fail else 0)


if __name__ == "__main__":
    main()
