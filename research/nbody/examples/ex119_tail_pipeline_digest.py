#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex119_tail_pipeline_digest.py
==============================
Jedno uruchomienie spinające **ex114–ex118** (ogon / linearyzacja / RMSE):

  - triada φ-FP, pasma [20,35] i [25,32]
  - zeta (ex117), RMSE/A (ex114 fit)
  - r_L*(ε) dla ε ∈ {0.12, 0.08} (ex116)
  - Pearson(zeta, RMSE/A) na 6 punktach (ex118)

Wyjścia:
  _outputs/ex119_pipeline_digest.csv   (metric, value, unit, note)
  _outputs/ex119_pipeline_digest.png     (scatter + słupki zeta)

Testy:
  D1: CSV + PNG
  D2: Pearson > 0.9 (spójność z ex118 na tym samym protokole)
  D3: max RMSE/A na [25,32] < 0.12

Uruchomienie: python ex119_tail_pipeline_digest.py
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


def load_mod(tag: str, fname: str):
    path = Path(__file__).resolve().parent / fname
    spec = importlib.util.spec_from_file_location(tag, path)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


def main():
    print("=" * 70)
    print("EX119: DIGEST — ogon tail (ex114–ex118)")
    print("=" * 70)

    ex114 = load_mod("ex119_ex114", "ex114_tail_phase_map.py")
    ex116 = load_mod("ex119_ex116", "ex116_tail_rL_from_epsilon.py")
    ex117 = load_mod("ex119_ex117", "ex117_linear_operator_residual.py")

    integrate_soliton = ex114.integrate_soliton
    fit_tail_detailed = ex114.fit_tail_detailed
    find_g0_star = ex114.find_g0_star
    zeta_ratio = ex117.zeta_ratio
    leftmost_r_L = ex116.leftmost_r_L
    R_R_FIXED = ex116.R_R_FIXED

    g0_star = find_g0_star()
    if g0_star is None:
        print("  BŁĄD: brak g₀*")
        sys.exit(1)

    PHI = ex114.PHI
    triad = [
        ("e", float(g0_star)),
        ("mu", float(PHI * g0_star)),
        ("tau", float(PHI**2 * g0_star)),
    ]

    profiles: list[tuple[np.ndarray, np.ndarray]] = []
    for _name, g0 in triad:
        r, g, _ = integrate_soliton(g0)
        profiles.append((r, g))

    digest: list[dict] = []

    def add(metric: str, value: float | str, unit: str = "", note: str = ""):
        digest.append(
            {
                "metric": metric,
                "value": value if isinstance(value, str) else f"{value:.8g}",
                "unit": unit,
                "note": note,
            }
        )

    add("g0_star_fp", g0_star, "1", "φ-FP ex114")
    add("g0_mu", PHI * g0_star, "1", "φ·g0*")
    add("g0_tau", PHI**2 * g0_star, "1", "φ²·g0*")

    r_L_012 = leftmost_r_L(profiles, R_R_FIXED, 0.12)
    r_L_008 = leftmost_r_L(profiles, R_R_FIXED, 0.08)
    add("r_L_star_eps_0.12", r_L_012 if r_L_012 is not None else float("nan"), "1", "ex116")
    add("r_L_star_eps_0.08", r_L_008 if r_L_008 is not None else float("nan"), "1", "ex116")

    zetas: list[float] = []
    rhos: list[float] = []
    labels_scatter: list[str] = []
    zeta_narrow_by_gen: dict[str, float] = {}

    for (name, g0), (r, g) in zip(triad, profiles):
        for band_name, r_lo, r_hi in BANDS:
            A, _B, _C, rmse = fit_tail_detailed(r, g, r_lo, r_hi)
            rho = rmse / max(A, 1e-15)
            z, _rms_h = zeta_ratio(r, g, r_lo, r_hi)
            if np.isfinite(z) and np.isfinite(rho):
                zetas.append(z)
                rhos.append(rho)
                labels_scatter.append(f"{name}\n{band_name}")
            if band_name == "win_25_32" and np.isfinite(z):
                zeta_narrow_by_gen[name] = z
            add(f"zeta_{name}_{band_name}", z, "1", "L[h]/RMS(h)")
            add(f"rmse_over_A_{name}_{band_name}", rho, "1", "fit cos/sin")

    zv = np.array(zetas, dtype=float)
    rv = np.array(rhos, dtype=float)
    if len(zv) >= 3 and np.std(zv) > 1e-14 and np.std(rv) > 1e-14:
        pearson = float(np.corrcoef(zv, rv)[0, 1])
    else:
        pearson = float("nan")
    add("pearson_zeta_rmse_over_A_n6", pearson, "1", "ex118-style")

    max_rho_narrow = max(
        float(digest[i]["value"])
        for i, row in enumerate(digest)
        if row["metric"].endswith("win_25_32") and row["metric"].startswith("rmse_over_A")
    )

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    csv_path = OUT_DIR / "ex119_pipeline_digest.csv"
    with open(csv_path, "w", encoding="utf-8") as f:
        f.write("metric,value,unit,note\n")
        for row in digest:
            f.write(
                f"{row['metric']},{row['value']},{row['unit']},{row['note']}\n"
            )

    print(f"  Zapisano: {csv_path}")
    print(f"  Pearson(zeta, RMSE/A) n=6: {pearson:.4f}")
    print()

    png_path = OUT_DIR / "ex119_pipeline_digest.png"
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.6))

        ax = axes[0]
        ax.scatter(zv, rv, c="tab:blue", s=80, edgecolors="black", linewidths=0.5)
        for i, lab in enumerate(labels_scatter):
            ax.annotate(lab.replace("\n", " "), (zv[i], rv[i]), fontsize=6, xytext=(4, 4), textcoords="offset points")
        if len(zv) >= 2 and np.std(zv) > 1e-14:
            m, b = np.polyfit(zv, rv, 1)
            xs = np.linspace(float(np.min(zv)), float(np.max(zv)), 50)
            ax.plot(xs, m * xs + b, "r--", lw=1.0, alpha=0.8, label=f"lin. fit (r={pearson:.3f})")
            ax.legend(fontsize=8)
        ax.set_xlabel("zeta = RMS(L[h])/RMS(h)")
        ax.set_ylabel("RMSE / A_tail")
        ax.set_title("ex119: zgodność liniówka ↔ ansatz")
        ax.grid(True, alpha=0.3)

        ax = axes[1]
        gens = ["e", "mu", "tau"]
        zs = [zeta_narrow_by_gen.get(g, float("nan")) for g in gens]
        ax.bar(gens, zs, color=["tab:blue", "tab:orange", "tab:green"], edgecolor="black")
        ax.set_ylabel("zeta")
        ax.set_title("Pas [25,32]: reszta L[h] względem h")
        ax.grid(True, axis="y", alpha=0.3)

        fig.suptitle("ex119 — digest łańcucha tail (triada φ-FP)", fontsize=10)
        plt.tight_layout()
        fig.savefig(png_path, dpi=140, bbox_inches="tight")
        plt.close(fig)
        print(f"  Wykres: {png_path}")
    except Exception as ex:
        print(f"  (matplotlib: {ex})")
    print()

    check(csv_path.is_file() and csv_path.stat().st_size > 80, "D1: CSV digest", str(csv_path))
    check(png_path.is_file(), "D1b: PNG digest", str(png_path))
    check(np.isfinite(pearson) and pearson > 0.9, "D2: Pearson > 0.9", f"r={pearson:.4f}")
    check(np.isfinite(max_rho_narrow) and max_rho_narrow < 0.12, "D3: max RMSE/A [25,32]", f"max={max_rho_narrow:.4f}")

    n_fail = sum(1 for _a, s, _b in RESULTS if s == "FAIL")
    print("=" * 70)
    print(f"EX119: {len(RESULTS) - n_fail}/{len(RESULTS)} testów PASS")
    print("=" * 70)
    sys.exit(1 if n_fail else 0)


if __name__ == "__main__":
    main()
