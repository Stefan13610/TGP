#!/usr/bin/env python3
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
"""
r3_phase7_phi0_screening_e2.py
================================

PURPOSE
-------
Sprawdzenie hipotezy uzytkownika 2026-05-01:
  "Φ₀ jako tło wszechświata z nakładania oddziaływań pól cząstek
   (kumulatywny proces). Liczba Eulera mogłaby pojawić się jako
   continuous limit dyskretnego procesu mnożnikowego w TGP."

KLUCZOWE ZNALEZISKO Z TGP (sek00_summary.tex:74-78):
  Φ₀ (bare)  = 168·Ω_Λ ≈ 115         [UV scale]
  Φ_eff      = Φ₀ · 3/14 ≈ 24.66    [IR scale, "ekranowany dielektryk"]

  Screening factor = 3/14 = 0.2143
  Inverse = 14/3 = 4.667

To JEST kumulatywny proces. Pytanie: czy 3/14 lub 14/3 ma związek z e²/4?

PLAN
----
1. Numeryczne porównanie 3/14, 14/3, X=e²/4 i pochodnych
2. Sprawdzenie czy bare/IR ratio = exp(coś)
3. Test czy kumulatywny szereg nieskończony może dać X = e²/4
4. Honest werdykt
"""

import math

# Stałe TGP
PHI0_BARE = 115.0       # = 168·Ω_Λ
PHI0_EFF = 24.66        # = 115 · 3/14
SCREENING = 3/14        # = 0.2143
INV_SCREENING = 14/3    # = 4.667

# Liczby fundamentalne
E = math.e
E_SQ = E**2
PI = math.pi
PHI_GOLDEN = (1 + math.sqrt(5)) / 2

# X z R3 (Faza 2)
X_R3 = E_SQ / 4         # 1.847

print("=" * 78)
print("  R3 Faza 7 — Φ₀ screening 3/14 vs e²/4")
print("=" * 78)
print()
print(f"  Φ₀ (bare)        = {PHI0_BARE}")
print(f"  Φ_eff (renorm)   = {PHI0_EFF}")
print(f"  Screening 3/14   = {SCREENING:.6f}")
print(f"  Inverse 14/3     = {INV_SCREENING:.6f}")
print(f"  X_R3 = e²/4      = {X_R3:.6f}")
print()


# ----------------------------------------------------------------
# SECTION 1: Direct numerical comparisons
# ----------------------------------------------------------------
print("=" * 78)
print("  SEKCJA 1: Direct numerical comparisons")
print("=" * 78)
print()

comparisons = [
    ("3/14", SCREENING),
    ("14/3", INV_SCREENING),
    ("(14/3) - e", INV_SCREENING - E),
    ("e^(π/2)", math.exp(PI/2)),
    ("(14/3) / e", INV_SCREENING / E),
    ("(14/3)² ", INV_SCREENING**2),
    ("e² · 3/14", E_SQ * SCREENING),
    ("e²/4 · (14/3)", X_R3 * INV_SCREENING),
    ("X · (14/3)", X_R3 * INV_SCREENING),
    ("8·X / (14/3)", 8*X_R3 / INV_SCREENING),
    ("ln(14/3)", math.log(INV_SCREENING)),
    ("ln(14/3) / (π/2)", math.log(INV_SCREENING)/(PI/2)),
    ("e^(3/14)", math.exp(SCREENING)),
    ("3/14 vs 1/e²", abs(SCREENING - 1/E_SQ)),
]

print(f"  {'Expression':<25} | {'Value':>14}")
print("  " + "-" * 45)
for name, val in comparisons:
    print(f"  {name:<25} | {val:14.6f}")


# ----------------------------------------------------------------
# SECTION 2: Test czy bare/IR ratio = exp(coś)
# ----------------------------------------------------------------
print()
print("=" * 78)
print("  SEKCJA 2: Czy 14/3 = exp(coś) ma natural origin?")
print("=" * 78)
print()

ln_inv = math.log(INV_SCREENING)
print(f"  ln(14/3) = {ln_inv:.6f}")
print()
print("  Kandydaci dla ln(14/3) = c:")

candidates = [
    ("π/2", PI/2),
    ("3/2", 1.5),
    ("ln(7/2) - ln(3/2) + 1", math.log(7/2) - math.log(3/2) + 1),
    ("e/π", E/PI),
    ("ln(5)", math.log(5)),
    ("3/π·ln(2)·e", 3/PI * math.log(2) * E),
    ("Pure ln(14/3)", ln_inv),  # tautology
    ("φ", PHI_GOLDEN),
    ("ln(4) + 1/4", math.log(4) + 0.25),
]

for name, val in candidates:
    diff = (val - ln_inv) / ln_inv * 100
    flag = " <<<" if abs(diff) < 1 else ""
    print(f"    {name:<30} = {val:9.5f}, diff {diff:+7.3f}%{flag}")


# ----------------------------------------------------------------
# SECTION 3: TGP-screening structure — gdzie 3/14 pochodzi
# ----------------------------------------------------------------
print()
print("=" * 78)
print("  SEKCJA 3: Strukturalne pochodzenie 3/14 w TGP")
print("=" * 78)
print()
print("Z sek00_summary.tex:74-77:")
print("  P(φ) = (β/7)φ⁷ - (γ/8)φ⁸    [potencjał akcji]")
print("  V(φ) = (γ/3)φ³ - (γ/4)φ⁴    [potencjał EOM]")
print("  P(1) = γ/56")
print("  V(1) = γ/12")
print("  P(1)/V(1) = 12/56 = 3/14 = SCREENING")
print()
print("To NIE jest e-related. To czysto algebraiczne z form akcji + potencjału.")
print()
print("  3/14 pochodzi z:")
print("    - (1/7 - 1/8) = 1/56  [P(1) numerator]")
print("    - (1/3 - 1/4) = 1/12  [V(1) numerator]")
print("    - 1/56 / (1/12) = 12/56 = 3/14")
print()
print("Czyli ekranowanie 3/14 to **strukturalna konsekwencja**")
print("wyboru K(φ)=φ⁴ + V(φ)=γφ³/3 - γφ⁴/4 + akcja P=K·V/dφ")
print("Brak e-content w tym łańcuchu — jest pure power-law algebra.")


# ----------------------------------------------------------------
# SECTION 4: Czy istnieje kumulatywny mechanizm który MOŻE dać e?
# ----------------------------------------------------------------
print()
print("=" * 78)
print("  SEKCJA 4: Hipotezy kumulatywne w TGP które mogłyby produkować e")
print("=" * 78)
print()
print("User's intuicja: e = lim(1+1/n)^n powstaje z kumulatywnego procesu.")
print()
print("Możliwe TGP-mechanizmy:")
print()
print("(A) RG flow Z(μ) = exp(∫γ d log μ)")
print("    Jeśli γ_Φ = const i flow over [μ_UV, μ_IR],")
print("    Z = (μ_IR/μ_UV)^γ — to jest power-law, nie exp.")
print("    Ale jeśli γ_Φ varies w specific way, Z = exp(coś).")
print()
print("(B) Soliton dressing through layers")
print("    R3 soliton ma profile g(r). Każda warstwa r → r+dr 'ubiera' pole.")
print("    Iteracyjnie: g(r+dr) = g(r) · (1 + Δ/n)^n → g · exp(Δ) gdy n→∞")
print("    To jest prawdziwy kumulatywny limit — ale jeszcze nie wykazany")
print("    explicit w R3 ODE.")
print()
print("(C) Partition function statystyczna")
print("    Z_substrate = Σ exp(-βE_n)")
print("    e pojawia się naturalnie. Ale to jest equilibrium statistical,")
print("    nie connected directly do R3 mass formula.")
print()
print("(D) Path integral nad fluktuacjami substrate")
print("    Z = ∫Dφ exp(-S[φ]). e pojawia się w generating functional.")
print("    Ale jego pojawienie w X = e²/4 wymagałoby specific computation.")
print()
print("(E) Kumulatywny scattering w substracie (user's analogy)")
print("    Particle 1 → field 1")
print("    Particle 2 (na tle 1) → field 2 + correction from 1")
print("    Particle N → field N + corrections from 1..N-1")
print("    Suma Σ corrections = exp(coś) w continuous limit")
print("    To FORMALNIE łączy z user's analogy, ale brak explicit derivation")
print("    w TGP-current formalism.")


# ----------------------------------------------------------------
# SECTION 5: Numerical search — czy 3/14 ma e-related connection
# ----------------------------------------------------------------
print()
print("=" * 78)
print("  SEKCJA 5: Direct test 3/14 vs e²/4 numerically")
print("=" * 78)
print()
print(f"  X = e²/4 = {X_R3:.6f}")
print(f"  3/14    = {SCREENING:.6f}")
print(f"  X · 3/14 = {X_R3 * SCREENING:.6f}")
print(f"  X / (3/14) = X · 14/3 = {X_R3 * INV_SCREENING:.6f}")
print(f"  4 / (14/3) = 12/14 = 6/7 = {6/7:.6f}")
print()

# Sprawdz czy X · 14/3 ma fizyczne znaczenie
val = X_R3 * INV_SCREENING
print(f"  X · 14/3 = e²/4 · 14/3 = 14e²/12 = 7e²/6 = {val:.6f}")
print()
print(f"  7e²/6 vs Φ_eff = 24.66:")
print(f"    Φ_eff / (7e²/6) = {PHI0_EFF / val:.4f}")
print(f"    To NIE ma czystego matchu.")
print()

# Spróbuj 8·X
print(f"  8·X = 2e² = {8*X_R3:.4f}")
print(f"  Compare: Φ₀_eff · ratio = {PHI0_EFF}")
print(f"  Φ_eff / (8X) = {PHI0_EFF / (8*X_R3):.4f}")
print(f"    1.668 — czy to coś znaczy?")
print()

# 8X = 2e² ≈ 14.78
# Φ_eff = 24.66
# Φ_eff / 8X ≈ 1.668
# 1.668 ≈ φ² = 2.618 nope, ≈ √(8/π) = 1.596 no, ≈ 5/3 = 1.667 YES!
ratio = PHI0_EFF / (8*X_R3)
print(f"  Φ_eff / 8X vs 5/3:")
print(f"    {ratio:.6f} vs {5/3:.6f}, diff {(ratio*3/5 - 1)*100:+.3f}%")
print()
print(f"  >> POSSIBLE: Φ_eff = (5/3) · 8 · (e²/4) = (40/12) · e² = 10e²/3")
print(f"     10·e²/3 = {10*E_SQ/3:.4f}  vs Φ_eff = 24.66")
print(f"     Match: {(10*E_SQ/3) / PHI0_EFF:.6f}, diff {((10*E_SQ/3)/PHI0_EFF - 1)*100:+.3f}%")


# ----------------------------------------------------------------
# SECTION 6: HONEST conclusion
# ----------------------------------------------------------------
print()
print("=" * 78)
print("  SEKCJA 6: HONEST conclusion")
print("=" * 78)
print()

print(f"""
  CO ODKRYŁEM:

  1. TGP **MA** kumulatywny mechanizm dla Φ₀ — bare→IR screening.
     Bare Φ₀ = 168·Ω_Λ ≈ 115
     Φ_eff = Φ₀ · 3/14 ≈ 24.66 (after dielectric screening)

     Twoja intuicja "Φ₀ z nakładania oddziaływań" jest CONSISTENT z
     TGP-formalizmem.

  2. Strukturalnie 3/14 NIE jest e-related. Pochodzi z algebraic ratio:
     P(1)/V(1) = (γ/56)/(γ/12) = 12/56 = 3/14
     gdzie P(φ) = β/7·φ⁷ - γ/8·φ⁸ (akcja) i V(φ) = γ/3·φ³ - γ/4·φ⁴ (EOM).
     To pure power-law algebra, brak exp/log structure.

  3. Numerycznie znalazłem CIEKAWĄ relację:
     Φ_eff ≈ 10·e²/3 = 9.05
     ALE to FAŁSZ — Φ_eff = 24.66, 10e²/3 = 24.63 — diff +0.13%!

     Hmm, czyli FAKTYCZNIE: Φ_eff ≈ 10·e²/3 z diff 0.13%!
     To oznaczałoby Φ_eff ≈ (10/3)·e², czyli e² jest "fundamental" w TGP.
""")

# Verify
predicted = 10 * E_SQ / 3
print(f"  Verify: 10·e²/3 = {predicted:.6f}")
print(f"  Φ_eff (TGP)    = {PHI0_EFF:.6f}")
print(f"  Diff = {(predicted/PHI0_EFF - 1)*100:+.4f}%")
print()

if abs(predicted - PHI0_EFF) / PHI0_EFF < 0.005:
    print(f"  >>> SURPRISING MATCH: Φ_eff ≈ (10/3)·e² do 0.13%!")
    print(f"      Jeśli TO JEST PRAWDZIWE (nie coincidence), wtedy:")
    print(f"        Φ_eff = (10/3) · e² ⟹ e² = (3/10) · Φ_eff")
    print(f"        X = e²/4 = (3/40) · Φ_eff = 0.075 · Φ_eff")
    print()
    print(f"  Sprawdz: 0.075 · 24.66 = {0.075 * PHI0_EFF:.4f}")
    print(f"           X = {X_R3:.4f}")
    print(f"           Match ratio: {0.075 * PHI0_EFF / X_R3:.6f}")

# Sprawdz czy e² to dokładnie 3·Φ_eff/10 lub coś podobne
print()
ratio_phi_eff_e2 = PHI0_EFF / E_SQ
print(f"  Φ_eff / e² = {ratio_phi_eff_e2:.6f}")
print(f"  Compare with 10/3 = {10/3:.6f}")
print(f"  Diff = {(ratio_phi_eff_e2/(10/3) - 1)*100:+.4f}%")

print()
print(f"""
  KOLEJNE PYTANIE:
    Czy Φ_eff = (10/3)·e² (jeśli prawdziwe) ma natural origin w TGP?
    10/3 = ? z TGP-algebry
    - 10 = 7+3 = 5·2 = ?
    - 3 = 3 (genaracje? wymiar?)
    - 10/3 = 3.333... = ?

  Bez derywacji 10/3, ten match może być coincidence.

  HONEST WERDYKT:
  Twoja intuicja "Φ₀ jako kumulatywne pole" jest wzmocniona przez fakt
  że TGP ma bare→IR screening. ALE numerycznie 3/14 NIE jest e-related,
  jest pure algebraic. Hipoteza Φ_eff ≈ (10/3)·e² wymaga dalszej
  derywacji 10/3 — ZNOWU mamy "promising hint, no proof".

  Najczystsza ścieżka: derive Φ_eff z explicit substrate calculation
  (statistical mechanics, partition function, RG flow) i sprawdzić czy
  e² wyłania się naturalnie w wyniku.
""")
