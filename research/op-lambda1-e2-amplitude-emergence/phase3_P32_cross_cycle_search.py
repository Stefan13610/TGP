#!/usr/bin/env python3
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
"""
phase3_P32_cross_cycle_search.py
==================================

PURPOSE
-------
λ.1 Phase 3 P3.2: Cross-cycle search dla e² traces w innych TGP sektorach.

Pytanie: czy e² (= Euler² = 7.389) pojawia się w innych sektorach poza
R3 charged lepton mass formula?

SECTORS DO SPRAWDZENIA
----------------------
1. **Neutrina (ζ.1)**: K_ν=1/2 Majorana — może e² lub inny wykładnik?
2. **Newton constant (χ.1)**: G_N derivation — czy zawiera e²?
3. **α-em (dodatekO_u1)**: phase sector — confirm że NIE ma e_Euler
4. **K-taxonomia Brannen**: r=√2, K_lep=2/3 — algebraic structures

PASS CRITERION
--------------
"Identify ≥1 dodatkowy sektor z e² hint LUB strukturalna izolacja
R3-amplitude (e² specific tylko dla charged K=2/3)"

Autor: λ.1 Phase 3 P3.2
Data: 2026-05-02
"""

import math
import numpy as np

E = math.e
E_SQ = E**2
PI = math.pi
PHI = (1 + math.sqrt(5)) / 2

print("=" * 78)
print("  λ.1 P3.2 — Cross-cycle search dla e² traces")
print("=" * 78)
print()


# ----------------------------------------------------------------
# SECTION 1: NEUTRINOS (ζ.1) — K=1/2 Majorana
# ----------------------------------------------------------------
print("=" * 78)
print("  SEKCJA 1: Neutrinos (ζ.1, K_ν=1/2 Majorana)")
print("=" * 78)
print()

# Z ζ.1 results (sek00 + audyt)
M_NU_1 = 0.80   # meV
M_NU_2 = 8.65
M_NU_3 = 50.11

print(f"  Neutrino masses (ζ.1):")
print(f"    m_1 = {M_NU_1} meV, m_2 = {M_NU_2} meV, m_3 = {M_NU_3} meV")
print()

# Mass ratios
R21 = M_NU_2 / M_NU_1
R31 = M_NU_3 / M_NU_1
R32 = M_NU_3 / M_NU_2

print(f"  Mass ratios:")
print(f"    m_2/m_1 = {R21:.4f}")
print(f"    m_3/m_1 = {R31:.4f}")
print(f"    m_3/m_2 = {R32:.4f}")
print()

# Test: czy ratios match formuła c·(g₀_2/g₀_1)^n_ν dla różnych n_ν
# Z Phase 1 L1.4: dla phi-drabinki (g₀_2 = phi·g₀_1), m_2/m_1 = (A_2/A_1)²·phi^n
# Brak znanego g₀_ν1 — search

print(f"  Test: czy istnieje n_ν które fittuje neutrino ratios z phi-drabinką?")
print(f"  (A_ν2/A_ν1)² · phi^n_ν = m_ν2/m_ν1 = {R21:.4f}")
print()
print(f"  Inne kandydaci dla n_ν (jeśli phi-drabinka):")

phi_log = math.log(PHI)
log_R21 = math.log(R21)
log_R31 = math.log(R31)

# Bez znanego g₀_ν1, można tylko skanować ratio (A²/A¹)² ·phi^n = R21
# Załóż że (A_ν2/A_ν1)² ≈ jakaś konstanta — skan różnych wartości
print(f"  {'(A2/A1)^2':>10} | {'n_ν z m_2/m_1':>15} | {'n_ν z m_3/m_1':>15} | {'consistent?':>12}")
print("  " + "-" * 70)

for A_ratio_sq in [0.5, 1.0, 2.0, 5.0, 10.0]:
    n_from_21 = (log_R21 - math.log(A_ratio_sq)) / phi_log
    # For m_3/m_1, need (A_3/A_1)^2 — assume same scaling (A_3/A_1 = phi^2 · A_ratio_sq)
    # Tak naprawdę nie wiemy A_ν2/A_ν1 vs A_ν3/A_ν1 - to dwa niezależne parametry
    # Sprawdz dla "uniform scaling" gdzie A_νi ~ phi^{i-1} · const
    A_31_sq = A_ratio_sq * PHI**2  # naive uniform
    n_from_31 = (log_R31 - math.log(A_31_sq)) / (2*phi_log)
    consistent = abs(n_from_21 - n_from_31) < 0.5
    print(f"  {A_ratio_sq:10.4f} | {n_from_21:15.4f} | {n_from_31:15.4f} | "
          f"{'YES' if consistent else 'NO':>12}")

print()
print(f"  e² = {E_SQ:.4f}")
print(f"  e²/2 = {E_SQ/2:.4f} (charged lepton n)")
print(f"  e²/4 = {E_SQ/4:.4f}")
print(f"  e = {E:.4f}")
print()

# Sprawdz: dla A_ratio_sq = 1 (i.e. A_ν2/A_ν1 = 1, equal amplitudes),
# n_from_21 = log(R21)/log(phi) = log(10.81)/0.481 = 4.95
# Hmm 4.95 ≈ 5 = N_lepton_states? Lub 5e/φ = 5·1.68 = 8.4? Nope.
# Nie wygląda jak e²/2.

n_eq_amp = log_R21 / phi_log
print(f"  Dla A_ν2/A_ν1 = 1 (równe amplitudy): n_ν = log(R21)/log(phi) = {n_eq_amp:.4f}")
print(f"  Compare: e² = {E_SQ:.4f}, e = {E:.4f}, π = {PI:.4f}, 5 = 5")
print()
print(f"  Najbliższy match: n_ν ≈ 5 (diff {(n_eq_amp - 5)/5*100:+.1f}%)")
print(f"  N_lepton states (3 charged + 3 neutrinos = 6, nope), but 5 = ?")
print()

# Sprawdz czy K_ν=1/2 daje natural exponent
# Charged: K_lep = 2/3, exponent n_lep(α=2) = e²/2
# Neutrino: K_ν = 1/2, exponent n_ν(α=2) = ?
# Hipoteza: n(K) = e²·(1-K) ?
# K=2/3: n = e²/3? = 2.46 (nope, charged daje 3.69)
# K=1/2: n = e²/2 = 3.69 (charged value!)
# Nope, no symmetric formula immediately

# Inna hipoteza: n(K) = e²·(1-α/4)·g(K) gdzie g(K) zależy od K
# Charged g(2/3) = 1, więc n(2,2/3) = e²/2
# Co byłoby g(1/2) by dać n_ν matching neutrina?

print(f"  Hipoteza λ.1 cross-test: n_ν = e²/2 · g(K_ν) gdzie g(K=2/3)=1")
print(f"  Numerical n_ν z neutrino masses ≈ {n_eq_amp:.2f}")
print(f"  Implied g(K_ν=1/2) = n_ν / (e²/2) = {n_eq_amp / (E_SQ/2):.4f}")
print()
print(f"  g_neutrino ≈ 1.34 — does it match anything?")
print(f"  Compare: K_lep/K_ν = (2/3)/(1/2) = 4/3 = {4/3:.4f}")
print(f"  Match: 1.34 vs 1.333... DIFF +0.5%!")
print()
print(f"  HIPOTEZA: g(K) = 2/(3K)? Dla K=2/3: g = 2/(3·2/3) = 1 ✓")
print(f"           Dla K=1/2: g = 2/(3·1/2) = 4/3 ≈ 1.333 (consistent z 1.34?)")
print()
print(f"  Albo g(K) = K_lep/K = (2/3)/K? Identyczna formuła.")
print()
print(f"  Test przy g(K) = 2/(3K):")
print(f"    n_ν = e²/2 · 4/3 = (2/3)·e² = {(2/3)*E_SQ:.4f}")
print(f"    PASTA z numerical n_ν = {n_eq_amp:.4f} — match {(((2/3)*E_SQ)/n_eq_amp - 1)*100:+.2f}%")


# ----------------------------------------------------------------
# SECTION 2: Newton constant (χ.1)
# ----------------------------------------------------------------
print()
print("=" * 78)
print("  SEKCJA 2: Newton constant (χ.1)")
print("=" * 78)
print()

print("""
  Z χ.1 program (op-chi1-newton-constant-derivation):
    G_N = g* / (M_TGP² · ξ_grav)

  Gdzie:
    g* = 0.71 (AS NGFP fixed point)
    M_TGP ≈ 0.2845 · M_Pl ≈ 3.473·10¹⁸ GeV
    ξ_grav = TGP-strukturalny dim-less factor (do wyznaczenia)

  Pytanie: czy ξ_grav zawiera e²?
""")

# Z literatury TGP nie znamy ξ_grav explicit, więc spróbujmy
# wyprowadzić numerycznie z obserwacji
G_N_PDG = 6.6743e-11  # m³/(kg·s²)
M_PL_PDG = 1.220890e19  # GeV (reduced Planck mass = 2.435e18)

# Z relacji G_N · M_Pl² ≈ ℏc (units):
# G_N · M_Pl² = 1 (in natural units where ℏ=c=1, M_Pl is reduced)
# So G_N · M_Pl² · c²/ℏ → dimensions match by definition

# χ.1 hipoteza: G_N = g*/(M_TGP² · ξ_grav)
# z M_TGP = 0.2845 M_Pl:
# G_N = g*/((0.2845)² · M_Pl² · ξ_grav)
# G_N · M_Pl² = g*/(0.2845² · ξ_grav)
# 1 = 0.71 / (0.0809 · ξ_grav)
# ξ_grav = 0.71 / 0.0809 = 8.776

xi_grav_inferred = 0.71 / (0.2845**2)
print(f"  Inferred ξ_grav (jeśli G_N·M_Pl² = 1): {xi_grav_inferred:.4f}")
print()
print(f"  Compare with e²:")
print(f"    e² = {E_SQ:.4f}")
print(f"    ξ_grav / e² = {xi_grav_inferred / E_SQ:.4f}")
print(f"    ξ_grav vs e²: diff {(xi_grav_inferred/E_SQ - 1)*100:+.2f}%")
print()

if abs(xi_grav_inferred/E_SQ - 1) < 0.05:
    print(f"  >>> POSSIBLE MATCH: ξ_grav ≈ e² (within 5%)")
else:
    print(f"  Match z e² nie jest exact ({xi_grav_inferred:.2f} vs {E_SQ:.2f})")

# Sprawdz alt
print()
print(f"  Inne kandydaci dla ξ_grav = {xi_grav_inferred:.4f}:")
print(f"    e² · 1.188 = {E_SQ * 1.188:.4f} (lambda*-related?)")
print(f"    π² = {PI**2:.4f} (geom)")
print(f"    e² + 4/3 = {E_SQ + 4/3:.4f}")
print(f"    8.776 (numerical) ≈ 168/19.14 ≈ ? non-trivial")


# ----------------------------------------------------------------
# SECTION 3: α-em (phase sector confirm)
# ----------------------------------------------------------------
print()
print("=" * 78)
print("  SEKCJA 3: α-em (phase sector — NEGATIVE confirm)")
print("=" * 78)
print()

print("""
  Z dodatekO_u1_formalizacja.tex (Phase 6 result):
    α_em^(sub) = Φ₀/(8π·J_phase) ≈ 1/94.09 (Planck)

  Składniki:
    - Φ₀ = TGP fundamental scale (sek00:76)
    - J_phase = phase coupling constant
    - 8π = 2·(4π) z full sphere integration

  e_Euler NIE pojawia się w α_em derivation.

  Confirm: phase sector używa **Φ₀, J, π** (NIE e_Euler).
  λ.1 hipoteza ma być **specific dla amplitude sector**.

  Status: NEGATIVE confirm — e² izolowana w amplitude.
""")


# ----------------------------------------------------------------
# SECTION 4: K-taxonomia Brannen
# ----------------------------------------------------------------
print("=" * 78)
print("  SEKCJA 4: K-taxonomia Brannen (r=√2, K_lep=2/3, K_ν=1/2)")
print("=" * 78)
print()

print("""
  Z Brannen lattice + R3 closure:
    r = √2 (substrate ratio)
    K_lep = 2/3 (charged leptons Koide)
    K_ν = 1/2 (neutrinos Majorana)
    Φ_0 = 24.783 (Brannen lock)
    K_geo · φ⁴ kinetic prefactor

  Sprawdzamy czy e² wyłania się z Brannen K-taxonomy:
""")

# r = √2 — irrational ale algebraic
# Sprawdz różne kombinacje
candidates = [
    ("r² = 2", 2),
    ("K_lep = 2/3", 2/3),
    ("K_ν = 1/2", 1/2),
    ("K_lep + K_ν = 7/6", 7/6),
    ("1/K_ν - 1/K_lep = 2 - 3/2 = 1/2", 1/2),
    ("3·K_lep = 2", 2),
    ("e² · r = 10.45", E_SQ * math.sqrt(2)),
    ("e² · K_lep = e²·2/3 ≈ 4.93", E_SQ * 2/3),
    ("e² · K_ν = e²/2 = 3.69", E_SQ / 2),  # charged n!
    ("Φ_0/K_lep = 24.783/(2/3) = 37.17", 24.783 / (2/3)),
]

print(f"  {'Expression':<35} | {'Value':>10}")
print("  " + "-" * 50)
for name, val in candidates:
    print(f"  {name:<35} | {val:10.4f}")

print()
print(f"  KEY OBSERVATION:")
print(f"    e² · K_ν = e² · 1/2 = e²/2 = 3.6945  ← Charged lepton n(α=2)!")
print(f"    e² · K_lep = e² · 2/3 = 4.926  ← potential neutrino n?")
print()
print(f"  Hipoteza: n(K) = e² · K (alternative do n = e²/2)?")
print(f"    K_lep = 2/3 → n = 2e²/3 = {2*E_SQ/3:.4f} (charged value 3.69 NIE matche)")
print()
print(f"  Lub: n = e² · (1 - K) ?")
print(f"    K_lep = 2/3 → n = e²/3 = {E_SQ/3:.4f} (NIE matche 3.69)")
print(f"    K_ν = 1/2 → n = e²/2 = {E_SQ/2:.4f} (charged value!)")
print()

# Sprawdz: czy może K-charged i K-neutrino są SWAPPED?
# Może wykładnik e²/2 jest dla obu, ale "g₀ scaling" jest inny dla K=2/3 vs K=1/2
# Ta hipoteza została zfalsyfikowana w P1 L1.4

print()
print("  Wniosek z Brannen analysis:")
print("  - e² · K_ν = e²/2 = charged n(α=2) — INTERESUJĄCA zbieżność")
print("  - Może implikować że n(K) = e² · (1 - K) lub similar")
print("  - Charged: K=2/3 → n = e²·(1-2/3) = e²/3 = 2.46 (NOPE, charged daje 3.69)")
print("  - To NIE jest clean Brannen-derivation dla e²")


# ----------------------------------------------------------------
# SECTION 5: Summary cross-cycle
# ----------------------------------------------------------------
print()
print("=" * 78)
print("  SEKCJA 5: Summary — co znaleziono w cross-cycle search?")
print("=" * 78)
print()

print(f"""
  CO ZNALAZŁEM:

  1. NEUTRINA (ζ.1):
     - n_ν z naive fit (równe amplitudy A_ν2/A_ν1 = 1) ≈ 4.95
     - Hipoteza n(K) = e²/2 · 2/(3K) daje n_ν = (2/3)e² = {(2/3)*E_SQ:.4f}
     - Match z numerical: ratio {(2/3)*E_SQ/n_eq_amp:.4f}
     - Status: SUGESTYWNE, ale brak explicit g₀_ν kalibracji
     - PARTIAL MATCH

  2. NEWTON (χ.1):
     - Inferred ξ_grav = 8.78 (z G_N relations)
     - Compare e² = 7.39: diff {(xi_grav_inferred/E_SQ - 1)*100:+.2f}%
     - Status: NIE matche e² directly, ale w bliskości
     - NEGATIVE (no clean match)

  3. α-em (phase sector):
     - 8π/(Φ₀·J_phase·...) — brak e_Euler
     - Status: NEGATIVE confirm — e² izolowane w amplitude
     - NEGATIVE (jako cross-test, ale POSITIVE jako boundary)

  4. K-taxonomia (Brannen):
     - e² · K_ν = e²/2 = charged n(α=2) — interesująca zbieżność
     - Brak clean derivation
     - NEGATIVE (numerologiczne hint)

  PASS criterion P3.2:
    "Identify ≥1 dodatkowy sektor z e² hint LUB strukturalna izolacja
     R3-amplitude (e² specific tylko dla charged K=2/3)"

  Wynik:
    - **Strukturalna izolacja**: phase sector NIE ma e_Euler ✓
      (z dodatekO_u1, confirms boundary)
    - **Hint w neutrinach**: n_ν ≈ (2/3)e² jest sugesytywne ale weak
    - **Newton**: nie matche e² directly (8.78 vs 7.39)

  Status: **PARTIAL PASS** (0.5)
    - Boundary R3-amplitude vs phase sector confirmed (POSITIVE)
    - No clear additional sector z e² (NEGATIVE for expansion)
    - Sugestywne hint w neutrinach, ale wymaga osobnego cyklu

  IMPLIKACJA:
  - λ.1 hipoteza może być zawężona do: "e² jest specific dla R3 charged
    K=2/3 amplitude sector"
  - Cross-cycle linkage do innych sektorów (neutrina, Newton, α) nie
    jest fundamental — to są **separate** mechanisms
""")
