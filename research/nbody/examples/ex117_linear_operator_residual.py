#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex117_linear_operator_residual.py
===================================
Weryfikacja numeryczna linearyzacji z ex117_tail_linearization.md:

    L[h] := h'' + (2/r) h' + h,   h = g - 1.

Profil g(r) z tego samego ODE co ex114. Pochodne na jednolitej siatce r
(spline kubiczny scipy), żeby ograniczyć szum różnic skończonych.

Metryki na [r_a, r_b]:
  zeta = sqrt( mean( L[h]^2 ) ) / sqrt( mean( h^2 ) )   (pomijamy h≈0)

Wyjście: _outputs/ex117_linear_residual.csv

Testy:
  Q1: CSV
  Q2: na [25, 32] dla triady zeta < 2.5 (luźny próg — sprawdza sens linearyzacji)
  Q3: zeta(τ) <= zeta(e) + 1.0 — τ trudniejszy, ale nie dowolnie gorszy na tym odcinku

Uruchomienie: python ex117_linear_operator_residual.py
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

# Dwa pasma: szersze i węższe „głęboki tail”
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
    spec = importlib.util.spec_from_file_location("ex117_ex114", path)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


def zeta_ratio(
    r_raw: np.ndarray,
    g_raw: np.ndarray,
    r_lo: float,
    r_hi: float,
    n_uni: int = 400,
) -> tuple[float, float]:
    """
    Zwraca (zeta, rms_h) dla maski |h|>1e-6; jeśli brak punktów, (nan, nan).
    """
    m = (r_raw >= r_lo - 0.5) & (r_raw <= r_hi + 0.5)
    r_sub = r_raw[m]
    g_sub = g_raw[m]
    if len(r_sub) < 30:
        return float("nan"), float("nan")

    # Jednolita siatka w zakresie rzeczywistych danych wewnątrz [r_lo, r_hi]
    r0 = max(float(np.min(r_sub)), r_lo)
    r1 = min(float(np.max(r_sub)), r_hi)
    if r1 <= r0 + 0.5:
        return float("nan"), float("nan")

    ru = np.linspace(r0, r1, n_uni)
    try:
        from scipy.interpolate import interp1d

        g_ip = interp1d(
            r_sub,
            g_sub,
            kind="cubic",
            bounds_error=False,
            fill_value=np.nan,
        )
        gu = g_ip(ru)
    except Exception:
        return float("nan"), float("nan")

    valid = np.isfinite(gu)
    ru, gu = ru[valid], gu[valid]
    if len(ru) < 40:
        return float("nan"), float("nan")

    h = gu - 1.0
    dr = ru[1] - ru[0]
    hp = np.gradient(h, dr)
    hpp = np.gradient(hp, dr)
    L = hpp + (2.0 / ru) * hp + h

    mask = (ru >= r_lo) & (ru <= r_hi) & (np.abs(h) > 1e-6)
    if np.sum(mask) < 15:
        return float("nan"), float("nan")

    Lm = L[mask]
    hm = h[mask]
    num = float(np.sqrt(np.mean(Lm**2)))
    den = float(np.sqrt(np.mean(hm**2)))
    if den < 1e-12:
        return float("nan"), float("nan")
    return num / den, den


def main():
    print("=" * 70)
    print("EX117: reszta operatora liniowego L[h]=h''+2h'/r+h  (h=g-1)")
    print("=" * 70)

    ex114 = load_ex114()
    g0_star = ex114.find_g0_star()
    if g0_star is None:
        print("  BŁAD: brak g₀*")
        sys.exit(1)

    PHI = ex114.PHI
    integrate_soliton = ex114.integrate_soliton
    triad = [
        ("e", float(g0_star)),
        ("mu", float(PHI * g0_star)),
        ("tau", float(PHI**2 * g0_star)),
    ]

    rows = []
    zeta_tau_25_32 = float("nan")
    zeta_e_25_32 = float("nan")

    for name, g0 in triad:
        r, g, _ = integrate_soliton(g0)
        print(f"  {name}: g₀={g0:.5f}")
        for band_name, r_lo, r_hi in BANDS:
            z, rh = zeta_ratio(r, g, r_lo, r_hi)
            rows.append(
                {
                    "gen": name,
                    "band": band_name,
                    "r_lo": r_lo,
                    "r_hi": r_hi,
                    "zeta": z,
                    "rms_h": rh,
                }
            )
            print(f"      [{band_name}]  zeta = {z:.4f}  (rms|h|={rh:.5f})")
            if band_name == "win_25_32":
                if name == "tau":
                    zeta_tau_25_32 = z
                if name == "e":
                    zeta_e_25_32 = z
        print()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    csv_path = OUT_DIR / "ex117_linear_residual.csv"
    keys = list(rows[0].keys())
    with open(csv_path, "w", encoding="utf-8") as f:
        f.write(",".join(keys) + "\n")
        for row in rows:
            f.write(",".join(str(row[k]) for k in keys) + "\n")
    print(f"  Zapisano: {csv_path}")
    print()

    # Testy na węższym pasie (głębszy asymptotyczny tail)
    zetas_narrow = [
        rw["zeta"] for rw in rows if rw["band"] == "win_25_32" and np.isfinite(rw["zeta"])
    ]
    max_zeta_narrow = max(zetas_narrow) if zetas_narrow else float("nan")

    check(csv_path.is_file(), "Q1: CSV", str(csv_path))
    check(
        np.isfinite(max_zeta_narrow) and max_zeta_narrow < 2.5,
        "Q2: max zeta (triada, [25,32]) < 2.5",
        f"max={max_zeta_narrow:.4f}",
    )
    check(
        np.isfinite(zeta_tau_25_32)
        and np.isfinite(zeta_e_25_32)
        and zeta_tau_25_32 <= zeta_e_25_32 + 1.0,
        "Q3: zeta(tau) <= zeta(e)+1 na [25,32]",
        f"z_e={zeta_e_25_32:.4f} z_tau={zeta_tau_25_32:.4f}",
    )

    n_fail = sum(1 for _a, s, _b in RESULTS if s == "FAIL")
    print("=" * 70)
    print(f"EX117: {len(RESULTS) - n_fail}/{len(RESULTS)} testów PASS")
    print("=" * 70)
    sys.exit(1 if n_fail else 0)


if __name__ == "__main__":
    main()
