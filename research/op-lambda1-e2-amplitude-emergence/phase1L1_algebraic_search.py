#!/usr/bin/env python3
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
"""
phase1L1_algebraic_search.py
==============================

PURPOSE
-------
λ.1 Phase 1 sub-task L1.1: Algebraic search dla 10/3 w TGP integers.

CONTEXT
-------
Z why_n3/Phase 7: Φ_eff (cosmological 24.66) ≈ (10/3) · e² do diff 0.12%.
Pytanie: czy 10/3 ma natural TGP-algebraic origin?

TGP-INTEGERS (zebrane z core/sek*.tex):
- 3 (N_gen, d_spatial)
- 4 (Hobart-Derrick balance, d_spacetime)
- 7, 8 (z P(φ) = β/7·φ⁷ - γ/8·φ⁸)
- 12 (z V(1) = γ/12)
- 14, 56 (z 3/14 = 12/56 screening)
- 36 (z Φ_eff = 36·Ω_Λ)
- 168 (= |GL(3,F₂)|, Φ_0 bare = 168·Ω_Λ)
- 184 = 168 + 16 (?, niejasne)
- κ = 3/(4·Φ_eff), κ⁻¹ ≈ 32.88

PHYSICAL INTERPRETATIONS testowane dla 10:
- T_4 = sum(1..4) = 10 (triangular number)
- 4·5/2 = 10 (symmetric 2-tensor in 4D = liczba g_μν components)
- Poincare group dim = 6 (Lorentz) + 4 (translations) = 10
- SU(5) generators = 24, ale fundamentalna rep = 5 (nie 10)
- 10 = 8 + 2 = (gluons number) + 2

PHYSICAL INTERPRETATIONS testowane dla 3:
- N_gen (number of fermion generations)
- d_spatial (3 dim spatial)
- 3 from Φ_eff numerator (3/14 = P(1)/V(1))

PASS CRITERION: 10/3 = TGP-rational z explicit fizyczną motywacją,
match do Φ_eff/e² z diff <0.5%.

Autor: λ.1 Phase 1 L1.1
Data: 2026-05-01
"""

import math
from itertools import combinations, product

# Stałe
E = math.e
E_SQ = E**2
PI = math.pi
PHI_GOLDEN = (1 + math.sqrt(5)) / 2

# Φ_eff wartości w TGP (różne anchors)
PHI_EFF_COSMO = 24.66      # = 36·Ω_Λ (z sek00:77)
PHI_EFF_BRANNEN = 24.783    # Brannen lock 2026-05-01 (sek09:1077)
PHI_EFF_25 = 25.0           # Preferred PPN point (sek08a:672)

# Target
TARGET = 10/3               # = 3.333... user's empirical hit

print("=" * 78)
print("  λ.1 L1.1 — Algebraic search dla 10/3 w TGP integers")
print("=" * 78)
print()
print(f"Target: 10/3 = {TARGET:.6f}")
print()
print("Origin: Φ_eff (24.66) ≈ (10/3)·e² z diff 0.12%")
print("Hypoteza: jeśli 10/3 jest TGP-algebraic, e² jest fundamentally derived.")
print()


# ----------------------------------------------------------------
# SECTION 1: TGP integer enumerations
# ----------------------------------------------------------------
print("=" * 78)
print("  SEKCJA 1: TGP integer rationals candidating dla 10/3")
print("=" * 78)
print()

TGP_INTEGERS = {
    "3": 3, "4": 4, "5": 5, "6": 6, "7": 7, "8": 8,
    "10": 10, "12": 12, "14": 14, "15": 15, "16": 16,
    "20": 20, "24": 24, "28": 28, "36": 36, "42": 42,
    "56": 56, "84": 84, "168": 168, "184": 184,
}

# Test rationals with denominators 1..10
candidates = []
for n_name, n in TGP_INTEGERS.items():
    for d in range(1, 11):
        ratio = n / d
        diff = abs(ratio - TARGET) / TARGET * 100
        if diff < 1.0:  # within 1%
            candidates.append((n_name, d, ratio, diff))

candidates.sort(key=lambda x: x[3])

print(f"  {'n':>5} | {'d':>3} | {'n/d':>9} | {'diff% from 10/3':>15}")
print("  " + "-" * 50)
for n_name, d, ratio, diff in candidates[:15]:
    print(f"  {n_name:>5} | {d:>3} | {ratio:9.5f} | {diff:15.4f}%")

print()


# ----------------------------------------------------------------
# SECTION 2: TGP physical interpretations dla 10
# ----------------------------------------------------------------
print("=" * 78)
print("  SEKCJA 2: Physical interpretations dla 10")
print("=" * 78)
print()

interpretations_10 = [
    ("T_4 = sum(1..4) (triangular number)", "topological count"),
    ("4·5/2 = 10 (symmetric 2-tensor 4D)", "g_μν independent components"),
    ("Poincare group dim = 6 + 4", "spacetime symmetry"),
    ("8 + 2 (gluons + EM photons?)", "gauge boson count"),
    ("2·5 = 10 (di-quarks?)", "particle physics"),
    ("3 + 7 (N_gen + (P-V)·56)", "TGP-internal sum"),
    ("4 + 6 (Hobart-Derrick + Lorentz dim)", "TGP×relativity"),
]

for interp, context in interpretations_10:
    print(f"  - {interp}: {context}")

print()
print("  KEY CANDIDATE: 4·5/2 = 10 = N_gravity_DOF (spin-2 tensor 4D)")
print("  Physical motivation:")
print("    - g_μν symmetric tensor w 4D ma 4·(4+1)/2 = 10 niezależnych components")
print("    - To jest dokładnie 'gravity degrees of freedom' przed gauge fixing")
print("    - W TGP grawitacja jest emergent z amplitude sector → naturalna relacja")
print()


# ----------------------------------------------------------------
# SECTION 3: TGP physical interpretations dla 3
# ----------------------------------------------------------------
print("=" * 78)
print("  SEKCJA 3: Physical interpretations dla 3")
print("=" * 78)
print()

interpretations_3 = [
    ("N_gen (number of fermion generations)", "fermion taxonomy"),
    ("d_spatial = 3", "spatial dimensions"),
    ("3 = P(1) numerator factor (γ/7 - γ/8 = γ/56, 12/56 = 3/14)", "TGP algebraic"),
    ("3 = κ·4·Φ_eff (z κ = 3/(4Φ_eff))", "TGP gravity coupling"),
    ("3 colors SU(3)_C", "gauge"),
    ("3 = N_gen = 3 dla R3 + bariera (why_n3)", "R3 result"),
]

for interp, context in interpretations_3:
    print(f"  - {interp}: {context}")

print()
print("  KEY CANDIDATE: 3 = N_gen (R3 result z why_n3)")
print("  Physical motivation:")
print("    - 3 generacje fermionów to fundamental TGP prediction (R3 closure)")
print("    - W amplitude sector emergent fermions = 3 generations")
print("    - 10/3 = (gravity DOF) / (fermion generations)")
print()


# ----------------------------------------------------------------
# SECTION 4: Hipoteza H_10/3 — physical interpretation
# ----------------------------------------------------------------
print("=" * 78)
print("  SEKCJA 4: Hipoteza H_10/3")
print("=" * 78)
print()
print(f"""
  HIPOTEZA H_10/3:
    10/3 = N_gravity_DOF / N_gen
         = (4·5/2) / 3
         = (symmetric 2-tensor 4D components) / (fermion generations)

  Fizyczna interpretacja:
    Φ_eff jest kumulatywnym polem (z bare→IR screening, sek00:77).
    Każdy fermion (3 generacje) "wnosi wkład" z grawitacyjnego coupling
    w 4D (10 components g_μν). Stosunek:
      Φ_eff per generation = (gravity DOF) / N_gen × (jakaś fundamental scale)
      Φ_eff = (10/3) · e²

    Gdzie e² byłaby "fundamental scale per gravity DOF per fermion generation".

  STATUS: SUGESTYWNE — physical motivation jest plausible, ale wymaga
  konkretnej derywacji (nie tylko counting argument).

  Co byłoby potrzebne żeby TO ZAMKNĄĆ:
    1. Pochodzić Φ_eff z explicit substrate calculation
    2. Pokazać że result faktoryzuje się jako (gravity DOF) × (e²) / N_gen
    3. Wyjaśnić dlaczego e² jest "fundamental scale" per DOF
""")


# ----------------------------------------------------------------
# SECTION 5: Numerical verification
# ----------------------------------------------------------------
print("=" * 78)
print("  SEKCJA 5: Numerical verification trzech wariantów Φ_eff")
print("=" * 78)
print()

predicted = TARGET * E_SQ
print(f"  (10/3) · e² = {predicted:.6f}")
print()

for name, phi in [
    ("Φ_eff cosmological (36·Ω_Λ)", PHI_EFF_COSMO),
    ("Φ_eff Brannen lock (sek09:1077)", PHI_EFF_BRANNEN),
    ("Φ_eff = 25 (PPN preferred)", PHI_EFF_25),
]:
    diff = (predicted/phi - 1) * 100
    flag = " <<< CLOSE" if abs(diff) < 0.5 else ""
    print(f"  {name:<40} = {phi}, diff {diff:+.4f}%{flag}")

print()


# ----------------------------------------------------------------
# SECTION 6: Alternative 10/3 sources w TGP
# ----------------------------------------------------------------
print("=" * 78)
print("  SEKCJA 6: Inne potencjalne źródła 10/3 w TGP")
print("=" * 78)
print()

# Test combinations
print(f"  Testy kombinacji TGP integers:")
print(f"  {'expression':<35} | {'value':>10} | {'10/3 match':>11}")
print("  " + "-" * 65)

test_exprs = [
    ("(7+3)/3 = N(P)+N_gen / N_gen", (7+3)/3),
    ("(8+2)/3", (8+2)/3),
    ("4·5/(3·2) = 10/6 = 5/3", (4*5)/(3*2)),  # = 10/3? No, 10/6 = 5/3
    ("(4+6)/3", (4+6)/3),
    ("Tetragonal: 4/(0.4·3)", 4/(0.4*3)),  # =3.333
    ("2·5/3", (2*5)/3),
    ("(P+V)/(P-V) for P=7,V=4: 11/3", (7+4)/3),
    ("Polynomial degree ratio: 8/2.4", 8/2.4),
    ("ζ(2)·something...", PI**2/6),  # = 1.645, far
    ("κ_TGP·Φ_eff·X", 0.030 * 24.66 * E_SQ/4),
]

for name, val in test_exprs:
    diff = abs(val - TARGET) / TARGET * 100
    flag = "  ✓ MATCH" if diff < 0.5 else ""
    print(f"  {name:<35} | {val:10.5f} | {diff:10.3f}%{flag}")

print()


# ----------------------------------------------------------------
# SECTION 7: HONEST conclusion L1.1
# ----------------------------------------------------------------
print("=" * 78)
print("  SEKCJA 7: PASS / FAIL judgment dla L1.1")
print("=" * 78)
print()

print(f"""
  Co znaleziono:

  1. 10/3 = T_4/N_gen ALBO 10/3 = (gravity DOF in 4D) / N_gen
     - 10 = 4·5/2 = liczba components symmetric 2-tensor w 4D
     - 3 = N_gen (z R3 result + N_gen=3 bariera)
     - Kombinacja **fizycznie motywowana** w TGP-context

  2. Match numerczny:
     - (10/3)·e² = {predicted:.4f}
     - Φ_eff_cosmological = 24.66, diff -0.12% ✓
     - Φ_eff_Brannen = 24.783, diff -0.62% (gorzej)

  3. ALE: argument "10 = gravity DOF / N_gen" jest **counting argument**,
     nie **derivation**. Bez konkretnego mechanizmu który łączy:
     (gravity DOF) × (fundamental scale) / (generations) = Φ_eff
     to jest **post-hoc rationalization**, nie mass derivation.

  PASS criterion (L1.1):
    "Znalezienie kombinacji TGP-integers z explicit fizyczną motywacją
     która daje 10/3 z diff <0.5% Φ_eff/e²"

  Spełnione:
    ✓ Kombinacja: 4·5/(2·3) = 10/3 ✓ (algebraic identity)
    ✓ Fizyczna motywacja: gravity DOF / N_gen ✓ (TGP-natural)
    ✓ Match: 0.12% dla Φ_eff cosmological ✓ (<0.5%)

  ALE caveat:
    ⚠️ Argument jest **counting**, nie **derivation**
    ⚠️ Wymagane jest pokazanie że Φ_eff *faktycznie* faktoryzuje się
       jako (10/3)·e² z fizycznego mechanizmu, nie tylko że LICZBOWO matche

  STATUS: **PARTIAL PASS**
    L1.1 spełnia formalne PASS criterion (counting + diff <0.5%)
    ALE motivation jest weak — wymaga L1.3 (partition function) lub
    L1.2 (β-function) żeby ZAMKNĄĆ derywację z fundamentu.

  Recommendation: count L1.1 jako 0.5 PASS (suggestive ale niedomknięte).
""")
