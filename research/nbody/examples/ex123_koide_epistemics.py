#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex123_koide_epistemics.py
==========================
Jawne porównanie dwóch domknięć przy **r₂₁ ze Ścieżki 9** (A_tail, φ-FP):

  r_31^(phi2) = (A(phi^2 g0*) / A(g0*))^4
  r_31^(Koide) = rozwiązanie Q_K(r_21_path, r_31) = 3/2

oraz sanity: Koide(R21_PDG) ≈ R31_PDG.

Wykorzystuje ex114 (integracja + fit tail [20,35]).

Wyjścia:
  _outputs/ex123_koide_epistemics.csv

Testy:
  E1: CSV
  E2: |Koide(R21_PDG) - R31_PDG| / R31_PDG < 0.002 (algebra + PDG)
  E3: |r_21_path - R21_PDG| / R21_PDG < 0.02 (φ-FP)
  E4: Q_K(r_21_path, r_31_Koide) ≈ 3/2
  E5: r_31^(phi2) i r_31^(Koide) skończone, > r_21_path

Teoria: ex123_koide_epistemics.md
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

R21_PDG = 206.768
R31_PDG = 3477.48

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
    spec = importlib.util.spec_from_file_location("ex123_ex114", path)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


def koide_r31_from_r21(r21: float) -> float | None:
    """Jak ex113: Q_K = 3/2, wybór fizyczny r_31 > r_21."""
    a = 1.0 + np.sqrt(r21)
    b = 1.0 + r21
    disc = 6.0 * a**2 - 3.0 * b
    if disc < 0:
        return None
    x_plus = 2.0 * a + np.sqrt(disc)
    x_minus = 2.0 * a - np.sqrt(disc)
    r31_plus = x_plus**2
    r31_minus = x_minus**2 if x_minus > 0 else None
    candidates = [r for r in [r31_plus, r31_minus] if r is not None and r > r21]
    if not candidates:
        return None
    return float(min(candidates))


def koide_Q(r21: float, r31: float) -> float:
    return (1.0 + np.sqrt(r21) + np.sqrt(r31)) ** 2 / (1.0 + r21 + r31)


def main():
    print("=" * 70)
    print("EX123: Epistemika — φ²·A_tail vs Koide (to samo r₂₁ ze Ścieżki 9)")
    print("=" * 70)

    ex114 = load_ex114()
    find_g0_star = ex114.find_g0_star
    integrate_soliton = ex114.integrate_soliton
    fit_tail_detailed = ex114.fit_tail_detailed
    PHI = ex114.PHI
    R_TAIL_L = ex114.R_TAIL_L
    R_TAIL_R = ex114.R_TAIL_R

    g0_star = find_g0_star()
    if g0_star is None:
        print("  BŁĄD: brak g₀*")
        sys.exit(1)

    def A_tail(g0: float) -> float:
        r, g, _ = integrate_soliton(g0)
        A, _B, _C, _rmse = fit_tail_detailed(r, g, R_TAIL_L, R_TAIL_R)
        return float(A)

    A_e = A_tail(float(g0_star))
    A_mu = A_tail(float(PHI * g0_star))
    A_tau_phi2 = A_tail(float(PHI**2 * g0_star))

    r21_path = (A_mu / A_e) ** 4
    r31_phi2 = (A_tau_phi2 / A_e) ** 4
    r31_koide = koide_r31_from_r21(r21_path)
    r31_koide_pdg = koide_r31_from_r21(R21_PDG)

    if r31_koide is None:
        print("  BŁĄD: brak rozwiązania Koide dla r_21_path")
        sys.exit(1)

    q_koide_closure = koide_Q(r21_path, r31_koide)
    rel_diff_phi2_k = abs(r31_phi2 - r31_koide) / r31_koide * 100.0

    print()
    print("  --- Wejście ze Ścieżki 9 (φ-FP, okno tail ex114) ---")
    print(f"  g₀*     = {g0_star:.6f}")
    print(f"  A(e)    = {A_e:.6f}")
    print(f"  A(μ)    = {A_mu:.6f}")
    print(f"  A(τ|φ²) = {A_tau_phi2:.6f}   [g₀ = φ²·g₀*]")
    print(f"  r₂₁_path = (A_μ/A_e)⁴     = {r21_path:.4f}   (PDG {R21_PDG})")
    print(f"  r₃₁^(φ²) = (A_τ/A_e)⁴     = {r31_phi2:.4f}   (PDG {R31_PDG})")
    print(f"  r₃₁^(K)  = Koide(r₂₁_path)= {r31_koide:.4f}")
    print(f"  |r₃₁^(φ²)−r₃₁^(K)|/r₃₁^(K) = {rel_diff_phi2_k:.3f}%")
    print()
    print("  --- Sanity: Koide na masach PDG ---")
    print(f"  Koide(R21_PDG) = {r31_koide_pdg:.4f}  (PDG R31 {R31_PDG})")
    print(f"  Q_K(R21_PDG, R31_PDG) = {koide_Q(R21_PDG, R31_PDG):.8f}")
    print()

    rows = [
        ("R21_PDG", R21_PDG, "PDG"),
        ("R31_PDG", R31_PDG, "PDG"),
        ("g0_star_fp", g0_star, "ex114"),
        ("A_tail_e", A_e, "fit"),
        ("A_tail_mu", A_mu, "fit"),
        ("A_tail_tau_phi2", A_tau_phi2, "g0=phi^2 g0*"),
        ("r21_path9", r21_path, "(A_mu/A_e)^4"),
        ("r31_phi2_rule", r31_phi2, "(A_tau/A_e)^4"),
        ("r31_koide_from_r21_path", r31_koide, "Q=3/2"),
        ("rel_diff_phi2_vs_koide_percent", rel_diff_phi2_k, "%"),
        ("Q_K_path_koide_closure", q_koide_closure, "should be 1.5"),
        ("r31_koide_from_R21_PDG", r31_koide_pdg, "sanity"),
    ]

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    csv_path = OUT_DIR / "ex123_koide_epistemics.csv"
    with open(csv_path, "w", encoding="utf-8") as f:
        f.write("metric,value,note\n")
        for name, val, note in rows:
            f.write(f"{name},{val},{note}\n")
    print(f"  Zapisano: {csv_path}")
    print()

    check(csv_path.is_file(), "E1: CSV", str(csv_path))
    check(
        r31_koide_pdg is not None
        and abs(r31_koide_pdg - R31_PDG) / R31_PDG < 0.002,
        "E2: Koide(R21_PDG) ≈ R31_PDG",
        f"delta={abs(r31_koide_pdg - R31_PDG) / R31_PDG * 100:.4f}%",
    )
    check(
        abs(r21_path - R21_PDG) / R21_PDG < 0.02,
        "E3: r_21_path blisko PDG (<2%)",
        f"delta={abs(r21_path - R21_PDG) / R21_PDG * 100:.4f}%",
    )
    check(abs(q_koide_closure - 1.5) < 1e-5, "E4: Q_K(path, r_31^K) = 3/2", f"Q={q_koide_closure:.8f}")
    check(
        r31_phi2 > r21_path and r31_koide > r21_path,
        "E5: r_31^(phi2), r_31^(K) > r_21_path",
        f"r21={r21_path:.2f}",
    )

    n_fail = sum(1 for _a, s, _b in RESULTS if s == "FAIL")
    print("=" * 70)
    print(f"EX123: {len(RESULTS) - n_fail}/{len(RESULTS)} testów PASS")
    print("=" * 70)
    sys.exit(1 if n_fail else 0)


if __name__ == "__main__":
    main()
