#!/usr/bin/env python3
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
"""
r3_phase6_alpha_em_connection.py
==================================

PURPOSE
-------
Faza 6 (REDO) — sprawdzenie hipotezy uzytkownika 2026-05-01:

  "X = e²/4 przypomina stałą struktury subtelnej w naturalnych jednostkach
   gdzie c=1, ε₀=1, ℏ=1, i 4π → 4 (przez TGP-modyfikację radial measure)."

CONTEXT
-------
Stała struktury subtelnej ma KLUCZOWO RÓŻNE zapisy w róznych systemach
jednostek:

  SI:                α = e²/(4πε₀ℏc) ≈ 1/137.036
  Heaviside-Lorentz: α = e²/(4π)    [ℏ=c=ε₀=1]
  Gaussian CGS:      α = e²/(ℏc) = e²  [ℏ=c=1, no rationalization]
  Atomic units:      α = 1/c

W każdym z tych zapisów "e" to ŁADUNEK elementarny (e_charge), nie liczba
Eulera (e_Euler ≈ 2.71828).

R3 EMPIRICAL
------------
Z Faza 2: X = e²_Euler / 4 ≈ 7.389/4 ≈ 1.847 (numerical fit < 0.1%)

Pytanie uzytkownika:
  Czy "e²/4" w R3 może być INTERPRETACJĄ stałej struktury subtelnej w
  TGP-naturalnych jednostkach gdzie:
  - 4π zostaje 4 (przez TGP-modyfikację — np. radial measure, Z₂ symmetry)
  - "e" to e_charge w TGP-units (NIE e_Euler) — może być różny od standardu
  - varying c, ℏ, ε₀ (TGP-emergent) tłumaczyłyby odchyły od kanonicznej α_SI

PLAN
----
1. Numeryczne porównanie: X_R3 vs α_em w różnych systemach jednostek
2. Algebraic search: w jakiej kalibracji α = e²/4 (zamiast e²/(4π))?
3. TGP-context: czy 4π → 4 ma natural origin (np. radial measure, sphere
   normalization)?
4. Implikacje dla R3 mass formula i emergent Dirac
5. HONEST report: czy hipoteza jest substantialna, czy numerologia?

UWAGA SEMANTYCZNA
-----------------
"e" jest **dwuznaczne**:
  - e_Euler ≈ 2.71828 (matematyczna stała exponential)
  - e_charge w SI ≈ 1.602·10⁻¹⁹ C (ładunek elektronu)
  - e_charge naturalne (Heaviside-Lorentz) ≈ 0.30282 (√(4π·α_SI))
  - e_charge naturalne (Gaussian) ≈ 0.08542 (√α_SI)

Te dwa są **NIEZWIĄZANE** standardowo. Pytanie czy w TGP jest jakiś
"unifying" zapis gdzie się łączą.

Autor: Faza 6 redo (po insight uzytkownika 2026-05-01).
"""

import numpy as np
import math

# Stałe fundamentalne (SI lub naturalne)
E_EULER = math.e                     # liczba Eulera 2.71828
E_EULER_SQ = E_EULER ** 2            # 7.389
PI = math.pi
ALPHA_SI = 7.2973525693e-3           # PDG 2024
INV_ALPHA = 1.0 / ALPHA_SI           # 137.036
PHI_GOLDEN = (1 + math.sqrt(5)) / 2

# X z Faza 2
X_R3 = E_EULER_SQ / 4                # 1.84726


print("=" * 78)
print("  Faza 6 REDO — X = e²/4 jako TGP-modified α-em?")
print("=" * 78)
print()
print(f"Empirical R3:  X = e_Euler² / 4 = {E_EULER_SQ:.4f} / 4 = {X_R3:.6f}")
print(f"Standard α:    α_SI = {ALPHA_SI:.10f} = 1/{INV_ALPHA:.4f}")
print()
print("Pytanie: czy istnieje jednostkowy system gdzie X = α-em?")
print()


# ----------------------------------------------------------------
# SECTION 1: Algebraic representations of fine structure constant
# ----------------------------------------------------------------
print("=" * 78)
print("  SEKCJA 1: Algebraiczne zapisy α w różnych systemach jednostek")
print("=" * 78)
print()

print("  Standard zapisy (gdzie e_c = ładunek elementarny):")
print()
print(f"    α_SI = e_c² / (4π·ε₀·ℏ·c)         = {ALPHA_SI:.6e}")
print()

# Heaviside-Lorentz natural (ℏ=c=ε₀=1)
e_c_HL = math.sqrt(4 * PI * ALPHA_SI)
print(f"    α_HL = e_c²/(4π)  [HL natural]    = {ALPHA_SI:.6e}")
print(f"           e_c (HL)   = √(4π·α_SI)    = {e_c_HL:.6f}")
print()

# Gaussian CGS (ℏ=c=1, no 4π rationalization)
e_c_G = math.sqrt(ALPHA_SI)
print(f"    α_G  = e_c²       [Gaussian]      = {ALPHA_SI:.6e}")
print(f"           e_c (G)    = √α_SI         = {e_c_G:.6f}")
print()

# Atomic units (ℏ=m_e=e_c=1)
print(f"    α_atomic = 1/c                     = {ALPHA_SI:.6e}")
print(f"           c (atomic) = 1/α            = {INV_ALPHA:.4f}")
print()


# ----------------------------------------------------------------
# SECTION 2: Hipoteza A — "e_TGP = e_Euler"
# ----------------------------------------------------------------
print("=" * 78)
print("  SEKCJA 2: Hipoteza A — czy ładunek w TGP units = e_Euler?")
print("=" * 78)
print()
print("  Jeśli 'e' w X = e²/4 to ładunek **w jakiejś TGP-natural calibracji**,")
print("  i równa się liczbie Eulera, sprawdźmy czy istnieje taki system.")
print()

print(f"  Required: e_TGP = e_Euler = {E_EULER:.6f}")
print(f"  Compare with HL natural: e_c (HL) = {e_c_HL:.6f}")
print(f"  Ratio: e_TGP / e_HL = {E_EULER / e_c_HL:.4f}")
print()
print(f"  W HL system, e_c ≈ 0.303. W TGP-hipotetycznym: e_TGP ≈ 2.718.")
print(f"  Ratio ≈ 9 — to JEST ekstremalnie inny ładunek.")
print()
print(f"  STRUKTURALNIE: TGP-natural units mają BARDZO różną normalizację")
print(f"  ładunku. Albo 'e' nie znaczy ładunek (jest to inny obiekt), albo")
print(f"  TGP-natural calibration jest poza standardem.")


# ----------------------------------------------------------------
# SECTION 3: Hipoteza B — "X to α-modified przez 4π → 4"
# ----------------------------------------------------------------
print()
print("=" * 78)
print("  SEKCJA 3: Hipoteza B — X = α z modyfikacją 4π → 4")
print("=" * 78)
print()

print("  Strukturalna obserwacja:")
print(f"    α_HL = e²/(4π) [standardowo, e to ładunek]")
print(f"    X_R3 = e²/4    [empirycznie, e to liczba Eulera]")
print()
print("  Jedyna formalna różnica: 4π → 4 (czynnik 1/π).")
print()
print("  Czy istnieje TGP-sense system gdzie 4π naturalnie staje się 4?")
print()

# Check rationalization factors
print("  Sources 4π w field theory:")
print("    - Surface area S² = 4π·r²")
print("    - Coulomb law w SI ma 1/(4π·ε₀)")
print("    - Heaviside rationalization: ε₀·c² = 1/μ₀, e²/(4πε₀) (coupling)")
print()
print("  Sources '4' (no π):")
print("    - 4 = 2D (spatial) + 1 + 1 (topological, e.g., dim group)")
print("    - 4 = Hobart-Derrick balance point dla TGP (Faza 2!)")
print("    - 4 = liczba wymiarów spacetime (3+1)")
print()
print("  W TGP, α = 4 może być Hobart-Derrick balance point — gdzie")
print("  X = e²/4 reprezentuje 'natural mass-coupling unit'.")


# ----------------------------------------------------------------
# SECTION 4: Numerical comparison X_R3 vs α_SI w różnych conversions
# ----------------------------------------------------------------
print()
print("=" * 78)
print("  SEKCJA 4: Numeryczna analiza X_R3 vs α_SI")
print("=" * 78)
print()

print(f"  X_R3 = {X_R3:.6f}")
print(f"  α_SI = {ALPHA_SI:.6e}")
print(f"  Ratio X / α = {X_R3 / ALPHA_SI:.6f}")
print()

# Test różne conversions
candidates = [
    ("π · X / α", PI * X_R3 / ALPHA_SI),
    ("X / α", X_R3 / ALPHA_SI),
    ("X · α", X_R3 * ALPHA_SI),
    ("4π · X", 4*PI * X_R3),
    ("X · 4π·α (HL)", X_R3 * 4*PI * ALPHA_SI),
    ("(4π·α) / X", 4*PI*ALPHA_SI / X_R3),
    ("e_Euler² · α", E_EULER_SQ * ALPHA_SI),
    ("4·X/(e_Euler²)", 4*X_R3/E_EULER_SQ),  # = 1 by construction
    ("X · π", X_R3 * PI),
    ("X / π", X_R3 / PI),
    ("log(X) / log(α)", math.log(X_R3) / math.log(ALPHA_SI)),
    ("1/(X · α)", 1.0 / (X_R3 * ALPHA_SI)),
]

print(f"  {'Ratio / Conversion':<22} | {'Value':>15}")
print("  " + "-" * 45)
for name, val in candidates:
    print(f"  {name:<22} | {val:15.6e}")


# ----------------------------------------------------------------
# SECTION 5: TGP-modified α with varying c, ℏ, ε₀
# ----------------------------------------------------------------
print()
print("=" * 78)
print("  SEKCJA 5: TGP-modified α z varying c_loc(ψ), ε₀_loc(ψ), ℏ_loc(ψ)")
print("=" * 78)
print()

print("  W TGP (Faza 1+3):")
print("    c_loc(ψ) = c · √A(ψ),  gdzie A(ψ) = (4-3ψ)/ψ")
print("    ψ_e = 0.950, ψ_μ = 1.155, ψ_τ = 1.288")
print()

# Compute c_loc / c dla generacji
def Apsi(psi):
    return (4 - 3*psi) / psi

psi_e = 0.95019
psi_mu = 1.15512
psi_tau = 1.28797

print(f"  {'gen':>5} | {'ψ':>9} | {'A(ψ)':>9} | {'c_loc/c':>9} | {'1/A':>9} | {'A/A(ψ=1)':>10}")
print("  " + "-" * 60)
for gen, psi in [('vac', 1.0), ('e', psi_e), ('μ', psi_mu), ('τ', psi_tau)]:
    A = Apsi(psi)
    print(f"  {gen:>5} | {psi:9.5f} | {A:9.5f} | {math.sqrt(A):9.5f} | "
          f"{1/A:9.5f} | {A/Apsi(1.0):10.5f}")

print()
print("  Hipoteza: α_TGP(ψ) = e_c²/(4π · ε₀_loc · ℏ_loc · c_loc)")
print("  Jeśli ε₀_loc i ℏ_loc varies with ψ tak by COMPENSATE c_loc·A,")
print("  α może pozostać stałe (locally), ale CONNECTION 4π → 4 wymaga")
print("  jakiejś topologicznej redukcji.")


# ----------------------------------------------------------------
# SECTION 6: Algebraic identity X·α = e²·α/4
# ----------------------------------------------------------------
print()
print("=" * 78)
print("  SEKCJA 6: Algebraic identity test")
print("=" * 78)
print()

# X = e²/4, więc X · α_HL = e_Euler²/4 · e_charge²/(4π) = e_Euler²·e_charge²/(16π)
prod = X_R3 * ALPHA_SI
print(f"  X · α_SI = {prod:.6e}")
print()

# X / α_HL gdy α_HL = e_charge²/(4π)
# = (e_Euler²/4) / (e_charge²/(4π)) = π · e_Euler² / e_charge²

# Sprawdz: w TGP-natural gdzie e_charge_TGP = e_Euler:
# α_TGP_HL = e_Euler²/(4π) = X · π / 4 · ... no

# Direct: jeśli α_TGP = X = e²/4 i jednocześnie α_TGP = e_c_TGP²/(4π)
# Wtedy: e_c_TGP² = 4π · X = 4π · 1.847 = 23.21
# e_c_TGP = √23.21 = 4.818
e_c_TGP_HL = math.sqrt(4 * PI * X_R3)
print(f"  Jeśli X = α_TGP w HL natural (α = e_c²/(4π)):")
print(f"    e_c_TGP = √(4π·X) = {e_c_TGP_HL:.4f}")
print(f"    vs e_c (HL) standard = {e_c_HL:.4f}")
print(f"    Ratio e_c_TGP / e_c_HL = {e_c_TGP_HL/e_c_HL:.4f}")
print()

# To by oznaczało że ładunek w TGP units jest 16x większy niż HL
print(f"  Ładunek w TGP-units byłby ~{e_c_TGP_HL/e_c_HL:.1f}x większy niż HL —")
print(f"  niefizyczne dla standardowego ładunku elektronu.")
print()

# Inna interpretacja: X to NIE α direct, ale RELATED przez specyficzne ratio
print("  ALTERNATIVA: X · π = α_HL · π² · ... — szukamy ladnego ratio")
print()

# Sprawdz czy X · α ma jakiś sens
print(f"  X · α_SI = {X_R3 * ALPHA_SI:.6e} = {X_R3 * ALPHA_SI:.10f}")
print(f"  X · α_SI · 137 = {X_R3 * ALPHA_SI * INV_ALPHA:.4f}")
print(f"  → To = X dokładnie! (bo α · α^(-1) = 1)")
print()

# Czy 1/X · 4π = α^{-1}?
test_val = 4 * PI / X_R3
print(f"  4π / X = {test_val:.6f}  vs α^{{-1}} = {INV_ALPHA:.4f}")
print(f"  Ratio = {test_val/INV_ALPHA:.6f}  → NIE matche (137 vs 6.8)")


# ----------------------------------------------------------------
# SECTION 7: Cleaner ratio search
# ----------------------------------------------------------------
print()
print("=" * 78)
print("  SEKCJA 7: Search for any clean integer/rational ratio X : α")
print("=" * 78)
print()

# Look for X = α · K dla małych K
K_test = X_R3 / ALPHA_SI
print(f"  X / α_SI = {K_test:.4f} = ?")
print()

# Check if this is close to clean rational
candidates_K = [
    ("253", 253),
    ("256 = 2^8", 256),
    ("2^8 · π/π", 256),
    ("e^4 · (1/e²)/2", math.e**4 / 2),
    ("π² · e^{0.5}", PI**2 * math.exp(0.5)),
    ("254", 254),
    ("250 · (1+1/100)", 250 * 1.01),
    ("e^(11/2)/e", math.exp(11/2 - 1)),
]

print(f"  K = X/α candidates:")
for name, val in candidates_K:
    diff_pct = (val/K_test - 1) * 100
    flag = " <<<" if abs(diff_pct) < 1 else ""
    print(f"    {name:<25} = {val:10.4f}  diff {diff_pct:+8.3f}%{flag}")
print()
print("  Brak czystego clean integer match.")


# ----------------------------------------------------------------
# SECTION 8: HONEST conclusion
# ----------------------------------------------------------------
print()
print("=" * 78)
print("  SEKCJA 8: HONEST conclusion")
print("=" * 78)
print()

print(f"""
  CO ODKRYŁEM:

  1. STRUKTURALNE PODOBIEŃSTWO: X = e²/4 ma identyczną formę jak α_HL = e²/(4π)
     z jedną formalną różnicą: czynnik 1/π.

  2. NUMERYCZNIE NIE MATCHUJĄ: X_R3 = 1.847, α_SI = 0.0073 (4 rzędy wielkości
     różnicy). To NIE jest po prostu "α w innych jednostkach" — wymagałoby
     ekstremalnej redefinicji ładunku (e_TGP ≈ 9·e_standard) co jest
     niefizyczne.

  3. SEMANTYKA "e": w X = e²/4, "e" to liczba Eulera (z empirical fit).
     W α_HL = e²/(4π), "e" to ładunek elektromagnetyczny. Te dwa są
     **NIEZWIĄZANE** standardowo.

  4. JEDNAK STRUCTURAL HINT: α często pojawia się w kontekście "fundamental
     coupling jednostkowe". Pojawienie się PODOBNEJ struktury w R3 (ładunek/(4π)
     vs liczba/4) jest **sugesytywne**, ale **nie evidence**.

  CO TO MOŻE OZNACZAĆ:

  Hipoteza H1 (NUMEROLOGY): X = e²/4 ma podobną formę do α_HL ale to
  COINCIDENCE form, brak fundamental connection.

  Hipoteza H2 (DEEP ANALOGY): R3 mass formula używa "kanonicznego coupling
  unit" które przypomina α — np. obie pochodzą z wave-function renormalization
  lub partition function evaluation. Wtedy 4π → 4 może być TGP-specyficzna
  redukcja (np. przez Z₂ symmetry).

  Hipoteza H3 (TGP-α emerguje): w TGP, fundamental coupling NIE jest α_SI
  (które wymaga zewnętrznego pomiaru). Może istnieje TGP-natural α' = e²/4
  gdzie 'e' ma fizyczne znaczenie inne niż ładunek elektronu — związane z
  exp(action) lub topological winding.

  STATUS: Hipoteza H2/H3 jest INTERESUJĄCA ale BRAK derywacji.
          Hipoteza H1 (numerologia) NIE WYKLUCZONA.

  CO BYŁOBY POTRZEBNE ŻEBY ZAMKNĄĆ:
  - Pochodzić e²/4 z konkretnego TGP loop integral lub partition function
  - Pokazać że "e" w R3 ma fundamental meaning (nie tylko empirical fit)
  - Połączyć z α_em (R3 powinno predict α_em z innego sektora?)

  HONEST WERDYKT: User's intuicja o α-link jest STRUKTURALNIE SUGESTYWNA
  ale wymaga osobnego cyklu (np. cykl który łączy R3 z op-alpha-fine-structure
  cyclem TGP). Bez tego, to pozostaje **suggestive analogy**, nie derywacja.
""")
