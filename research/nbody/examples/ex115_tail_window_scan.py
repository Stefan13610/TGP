#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex115_tail_window_scan.py
=========================
Skan prostokątów (r_L, r_R) dla fitu ogona na **triadzie** (g₀*, φ·g₀*, φ²·g₀*).

Dla każdego okna (na tych samych profilach co ex114 — jedna integracja na g₀):
  - max RMSE/A wśród trzech generacji
  - max średnie |g-1| w oknie
  - max |Δφ| względem fazy z **okna referencyjnego** (ex114: [20, 35])

Kolumna `pass_all`: wszystkie trzy warunki jak w ex114_tail_okna_podsumowanie.md:
  RMSE/A < 10%,  śr.|g-1| < 0.12,  |Δφ| vs ref < 0.35 rad

Wyjścia:
  _outputs/ex115_window_scan.csv
  _outputs/ex115_window_scan_admissible.png  (maska + opis)

Uruchomienie z katalogu examples/:
  python ex115_tail_window_scan.py
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

# Progi zsynchronizowane z ex114 / ex114_tail_okna_podsumowanie.md
REF_R_L = 20.0
REF_R_R = 35.0
MAX_RMSE_OVER_A = 0.10
MAX_MEAN_ABS_G1 = 0.12
MAX_DPHI_VS_REF = 0.35
MIN_WINDOW_WIDTH = 8.0
MIN_MASK_POINTS = 10

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
    spec = importlib.util.spec_from_file_location("ex114_tw", path)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


def main():
    print("=" * 70)
    print("EX115: SKAN OKIEN TAIL (r_L, r_R) — triada vs ref [20, 35]")
    print("=" * 70)

    ex114 = load_ex114()
    integrate_soliton = ex114.integrate_soliton
    fit_tail_detailed = ex114.fit_tail_detailed
    find_g0_star = ex114.find_g0_star
    wrap_angle_diff = ex114.wrap_angle_diff
    PHI = ex114.PHI

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    g0_star = find_g0_star()
    if g0_star is None:
        print("  BŁĄD: brak g₀* — przerwanie")
        sys.exit(1)

    triad = [
        ("e", float(g0_star)),
        ("mu", float(PHI * g0_star)),
        ("tau", float(PHI**2 * g0_star)),
    ]

    profiles: list[tuple[str, float, np.ndarray, np.ndarray]] = []
    phi_ref: list[float] = []
    print(f"  g₀* = {g0_star:.6f}")
    for name, g0 in triad:
        r, g, _ = integrate_soliton(g0)
        A, B, C, _rmse = fit_tail_detailed(r, g, REF_R_L, REF_R_R)
        phi_ref.append(float(np.arctan2(C, B)))
        profiles.append((name, g0, r, g))
        print(f"    {name}: g₀={g0:.5f}  φ_ref={np.degrees(phi_ref[-1]):.2f}°")
    print()

    r_L_vals = np.unique(
        np.concatenate([np.arange(15.0, 27.0, 1.5), np.array([REF_R_L])])
    )
    r_R_vals = np.unique(
        np.concatenate([np.arange(28.0, 40.0, 1.5), np.array([REF_R_R])])
    )

    rows = []
    for r_L in r_L_vals:
        for r_R in r_R_vals:
            if r_R - r_L < MIN_WINDOW_WIDTH:
                continue

            max_ra = 0.0
            max_mg = 0.0
            max_dphi = 0.0
            bad = False

            for i, (name, _g0, r, g) in enumerate(profiles):
                mask = (r >= r_L) & (r <= r_R)
                if np.sum(mask) < MIN_MASK_POINTS:
                    bad = True
                    break
                mean_g1 = float(np.mean(np.abs(g[mask] - 1.0)))
                max_mg = max(max_mg, mean_g1)

                A, B, C, rmse = fit_tail_detailed(r, g, float(r_L), float(r_R))
                if A < 1e-14 or not np.isfinite(rmse):
                    bad = True
                    break
                max_ra = max(max_ra, rmse / A)
                phi = float(np.arctan2(C, B))
                max_dphi = max(
                    max_dphi,
                    abs(wrap_angle_diff(phi - phi_ref[i])),
                )

            if bad:
                rows.append(
                    {
                        "r_L": float(r_L),
                        "r_R": float(r_R),
                        "width": float(r_R - r_L),
                        "max_rmse_over_A": float("nan"),
                        "max_mean_abs_g1": float("nan"),
                        "max_dphi_vs_ref_rad": float("nan"),
                        "pass_all": False,
                        "skip": True,
                    }
                )
                continue

            passes = (
                max_ra < MAX_RMSE_OVER_A
                and max_mg < MAX_MEAN_ABS_G1
                and max_dphi < MAX_DPHI_VS_REF
            )
            rows.append(
                {
                    "r_L": float(r_L),
                    "r_R": float(r_R),
                    "width": float(r_R - r_L),
                    "max_rmse_over_A": max_ra,
                    "max_mean_abs_g1": max_mg,
                    "max_dphi_vs_ref_rad": max_dphi,
                    "pass_all": bool(passes),
                    "skip": False,
                }
            )

    csv_path = OUT_DIR / "ex115_window_scan.csv"
    keys = list(rows[0].keys())
    with open(csv_path, "w", encoding="utf-8") as f:
        f.write(",".join(keys) + "\n")
        for row in rows:
            f.write(",".join(str(row[k]) for k in keys) + "\n")
    print(f"  Zapisano: {csv_path}")

    valid = [rw for rw in rows if not rw["skip"]]
    adm = [rw for rw in valid if rw["pass_all"]]
    print(f"  Okna poprawne (wystarczająco punktów): {len(valid)}")
    print(f"  Okna PASS (wszystkie progi): {len(adm)}")
    if adm:
        best = min(adm, key=lambda x: x["max_rmse_over_A"])
        print(
            f"  Najlepsze po max(RMSE/A): r_L={best['r_L']}, r_R={best['r_R']}, "
            f"RMSE/A={best['max_rmse_over_A']:.4f}"
        )
    print()

    # --- wykres: siatka r_L × r_R ---
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        rL_u = sorted(set(rw["r_L"] for rw in rows))
        rR_u = sorted(set(rw["r_R"] for rw in rows))
        mat = np.full((len(rL_u), len(rR_u)), np.nan)
        for rw in rows:
            if rw["skip"]:
                continue
            i = rL_u.index(rw["r_L"])
            j = rR_u.index(rw["r_R"])
            mat[i, j] = 1.0 if rw["pass_all"] else 0.0

        fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.8))

        ax = axes[0]
        X, Y = np.meshgrid(rR_u, rL_u)
        pcm = ax.pcolormesh(X, Y, mat, shading="auto", cmap="RdYlGn", vmin=0, vmax=1)
        ax.axhline(REF_R_L, color="cyan", ls="--", lw=1.0, alpha=0.9)
        ax.axvline(REF_R_R, color="cyan", ls="--", lw=1.0, alpha=0.9)
        ax.plot([REF_R_R], [REF_R_L], "c*", ms=14, label="ref (20,35)")
        ax.set_xlabel("r_R")
        ax.set_ylabel("r_L")
        ax.set_title("Obszar dopuszczalny (1=zielone = pass_all)")
        ax.legend(loc="upper right", fontsize=8)
        plt.colorbar(pcm, ax=ax, label="pass", ticks=[0, 1])

        ax = axes[1]
        mat_rmse = np.full((len(rL_u), len(rR_u)), np.nan)
        for rw in rows:
            if rw["skip"]:
                continue
            i = rL_u.index(rw["r_L"])
            j = rR_u.index(rw["r_R"])
            mat_rmse[i, j] = rw["max_rmse_over_A"]
        pcm2 = ax.pcolormesh(X, Y, mat_rmse, shading="auto", cmap="viridis")
        ax.axhline(REF_R_L, color="white", ls="--", lw=0.8)
        ax.axvline(REF_R_R, color="white", ls="--", lw=0.8)
        ax.set_xlabel("r_R")
        ax.set_ylabel("r_L")
        ax.set_title("max(RMSE/A) na triadzie")
        plt.colorbar(pcm2, ax=ax)

        fig.suptitle(
            "ex115: skan okien tail (progi jak ex114: RMSE/A<10%, |g-1|<0.12, Δφ<0.35 vs [20,35])",
            fontsize=9,
        )
        plt.tight_layout()
        png_path = OUT_DIR / "ex115_window_scan_admissible.png"
        fig.savefig(png_path, dpi=140, bbox_inches="tight")
        plt.close(fig)
        print(f"  Wykres: {png_path}")
    except Exception as ex:
        print(f"  (matplotlib: {ex})")
    print()

    # --- testy ---
    ref_row = next(
        (
            rw
            for rw in rows
            if not rw["skip"]
            and abs(rw["r_L"] - REF_R_L) < 1e-6
            and abs(rw["r_R"] - REF_R_R) < 1e-6
        ),
        None,
    )
    check(csv_path.is_file() and csv_path.stat().st_size > 100, "S1: CSV skanu", str(csv_path))
    check(
        ref_row is not None and ref_row["pass_all"],
        "S2: okno referencyjne [20,35] ∈ pass_all",
        "brak wiersza" if ref_row is None else f"pass={ref_row['pass_all']}",
    )
    check(len(adm) >= 4, "S3: co najmniej 4 dopuszczalne okna", f"n={len(adm)}")

    n_fail = sum(1 for _a, s, _b in RESULTS if s == "FAIL")
    print("=" * 70)
    print(f"EX115: {len(RESULTS) - n_fail}/{len(RESULTS)} testów PASS")
    print("=" * 70)
    sys.exit(1 if n_fail else 0)


if __name__ == "__main__":
    main()
