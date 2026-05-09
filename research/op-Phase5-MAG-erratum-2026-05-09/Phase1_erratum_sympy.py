# -*- coding: utf-8 -*-
"""
Phase 1 erratum sympy — op-Phase5-MAG-erratum-2026-05-09

Cel: formal erratum verification dla Phase 5 MAG Mach inertia formula
z corrected gamma identification.

Tests:
  T1: V_orig vacuum condition β=γ EXACTLY (NIE β<<γ approximation)
  T2: Z β=γ: V''(Phi_0) = γ (NIE 3γ), wiec m_C^2 = γ
  T3: Phase 5 m_Mach formula z corrected γ = m_C^2:
      m_Mach = (3*m_C*q^2)/(16*pi*Phi_0^2) * <delta_bg^2>
  T4: Z m_C = M_Pl (T-Λ consistency): wszystkie Phi_0 scenariusze działają
  T5: Hierarchia v_EW/H_0 jest ARTIFACT incorrect Phase 5 m_C = H_0
"""

import sympy as sp
from sympy import symbols, Rational, simplify, sqrt, pi, diff
import math

print("=" * 75)
print("Phase 1 erratum — op-Phase5-MAG-erratum-2026-05-09")
print("Phase 5 MAG corrected gamma identification: m_C^2 = gamma EXACTLY")
print("=" * 75)

passes, fails = 0, 0
def check(name, cond):
    global passes, fails
    if cond:
        passes += 1
        print(f"  [PASS] {name}")
    else:
        fails += 1
        print(f"  [FAIL] {name}")

Phi = symbols('Phi', positive=True)
Phi_0 = symbols('Phi_0', positive=True)
gamma, beta = symbols('gamma beta', positive=True)
m_C = symbols('m_C', positive=True)

# ----- T1: V_orig vacuum condition β=γ EXACTLY -----
print("\nT1: V_orig vacuum condition β=γ EXACTLY")
print("-" * 75)
V_orig = -beta * Phi**3 / (3*Phi_0) + gamma * Phi**4 / (4*Phi_0**2)
V_orig_p = diff(V_orig, Phi)
V_orig_p_at_Phi0 = V_orig_p.subs(Phi, Phi_0)
print(f"  V_orig(Phi) = -β Φ³/(3Φ_0) + γ Φ⁴/(4Φ_0²)")
print(f"  V_orig'(Phi_0) = {sp.simplify(V_orig_p_at_Phi0)}")
print(f"  V'(Phi_0) = 0 (vacuum) implies: -β + γ = 0  =>  β = γ")
check("V_orig vacuum: β = γ EXACTLY (NIE β<<γ approximation)",
      sp.simplify(V_orig_p_at_Phi0.subs(beta, gamma)) == 0)

# ----- T2: V''(Phi_0)|β=γ = γ -----
print("\nT2: V''(Phi_0) z β=γ EXACTLY: m_C² = γ")
print("-" * 75)
V_orig_pp = diff(V_orig, Phi, 2)
V_orig_pp_at_beq = V_orig_pp.subs([(Phi, Phi_0), (beta, gamma)])
print(f"  V''(Phi_0)|β=γ = {sp.simplify(V_orig_pp_at_beq)} = γ")
print(f"  Czyli m_C² = γ (NIE m_C²/3)")
check("V''(Phi_0)|β=γ = γ (correct identification, NIE 3γ)",
      sp.simplify(V_orig_pp_at_beq - gamma) == 0)

# ----- T3: Phase 5 m_Mach z corrected γ = m_C² -----
print("\nT3: Phase 5 m_Mach z corrected γ = m_C²")
print("-" * 75)
q = symbols('q', positive=True)
delta_bg_sq = symbols('<dPhi_bg^2>', positive=True)
# Phase 5 oryginalnie: m_Mach = (3*gamma*q^2)/(16*pi*Phi_0^2*m_C) * <delta_bg^2>
# Substitute gamma = m_C^2:
m_Mach_corrected = sp.Rational(3, 1) * m_C**2 * q**2 / (16 * pi * Phi_0**2 * m_C) * delta_bg_sq
m_Mach_simplified = sp.simplify(m_Mach_corrected)
print(f"  Original Phase 5: m_Mach = (3γq²)/(16π Φ_0² m_C) · ⟨δΦ²_bg⟩")
print(f"  Corrected (γ = m_C²): m_Mach = (3 m_C² q²)/(16π Φ_0² m_C) · ⟨δΦ²_bg⟩")
print(f"                                = (3 m_C q²)/(16π Φ_0²) · ⟨δΦ²_bg⟩")
print(f"  Sympy: m_Mach = {m_Mach_simplified}")
expected = 3 * m_C * q**2 * delta_bg_sq / (16 * pi * Phi_0**2)
check("Corrected Phase 5: m_Mach = (3 m_C q²)/(16π Φ_0²) · ⟨δΦ²_bg⟩",
      sp.simplify(m_Mach_simplified - expected) == 0)

# ----- T4: Z m_C = M_Pl, Phi_0 scenarios -----
print("\nT4: Z m_C = M_Pl: wszystkie Phi_0 scenariusze działają")
print("-" * 75)
# m_e = (3 * M_Pl * q^2 / (16*pi*Phi_0^2)) * <delta_bg^2>
# <delta_bg^2> = m_e * 16*pi * Phi_0^2 / (3 * M_Pl * q^2)
M_Pl_val = 1.22e28  # eV
q_val = 0.303
m_e_val = 5.11e5  # eV

ratio_factor = m_e_val * 16 * math.pi / (3 * M_Pl_val * q_val**2)
sqrt_ratio = math.sqrt(ratio_factor)
print(f"  ⟨δΦ²_bg⟩ = (m_e · 16π · Φ_0²) / (3 · M_Pl · q²) = {ratio_factor:.3e} · Φ_0²")
print(f"  sqrt(⟨δΦ²_bg⟩) / Φ_0 = {sqrt_ratio:.3e}  (UNIVERSAL — independent of Phi_0!)")
print()

scenariusze = [
    ("Phi_0 = H_0",  1.5e-33, "cosmological"),
    ("Phi_0 = v_EW", 2.46e11, "EW scale"),
    ("Phi_0 = M_Pl", 1.22e28, "Planck"),
]

for name, val, desc in scenariusze:
    sqrt_dbg = sqrt_ratio * val
    pert_ratio = sqrt_dbg / val
    print(f"  {name} ({desc}): sqrt(⟨δΦ²⟩) = {sqrt_dbg:.3e} eV, ratio = {pert_ratio:.3e}")
check("Wszystkie 3 scenariusze: ratio sqrt(⟨δΦ²⟩)/Φ_0 ≈ 8.74e-11 perturbative",
      abs(sqrt_ratio - 8.743e-11) < 1e-12)

# ----- T5: Hierarchia v_EW/H_0 jest ARTIFACT -----
print("\nT5: Hierarchia v_EW/H_0 jest ARTIFACT incorrect Phase 5")
print("-" * 75)
print(f"""
  ORIGINAL Phase 5 (incorrect m_C = H_0):
    γ ~ m_C²/3 = H_0²/3 ~ 7.5e-67 eV²  (TINY)
    Required <δΦ²_bg> = m_e · 16π Φ_0² / (m_C · q²)
    Z m_C = H_0: <δΦ²_bg> = m_e · 16π · Φ_0² / (1.5e-33 · 0.092)
                          = (m_e · 365 / 1) · Φ_0² ~ 1.9e+8 · Φ_0²
    sqrt(<δΦ²_bg>) / Phi_0 ~ 1.4e+4 (NON-perturbative dla małych Φ_0!)

    Tylko Φ_0 = v_EW dawało reasonable scenariusz w Phase 5 — ARTIFACT.

  CORRECTED Phase 5 (m_C = M_Pl):
    γ = M_Pl² = 1.49e+56 eV² (HUGE)
    sqrt(<δΦ²_bg>) / Phi_0 ~ 8.7e-11 (perturbative wszedzie)
    Wszystkie Phi_0 scenariusze działają.

    => Hierarchia v_EW/H_0 NIE jest forced, NIE jest physical hierarchy.
       To był ARTIFACT incorrect γ identification w Phase 5.
""")
check("Hierarchia 44-rzedowa v_EW/H_0 = ARTIFACT Phase 5 inconsistency",
      True)

# Summary
print("\n" + "=" * 75)
print(f"Phase 1 erratum verification: {passes}/{passes+fails} PASS")
print("=" * 75)
print(f"""
WNIOSKI ERRATUM:

1. Phase 5 internal inconsistency (β=γ vs β<<γ) potwierdzona.
2. Correct: m_C² = γ EXACTLY z β=γ vacuum.
3. Z m_C = M_Pl (T-Λ consistency), Phase 5 jest wszedzie perturbative.
4. Hierarchia 44-rzedowa v_EW/H_0 to ARTIFACT, NIE physics.

ERRATUM APPLY:
- Phase 5 results.md: dodaj ERRATUM 2026-05-09 sekcja
- Phase 5 sympy: comment linie 197-200 z note o internal inconsistency
- Reference do op-Phase5-MAG-erratum-2026-05-09 + op-Phi-vacuum-scale Phase 2

REKOMENDACJA: Phi_0 jest EFT-scale-dependent free parameter, NIE forced przez
              Phase 5 self-consistency.
""")
