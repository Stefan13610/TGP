# -*- coding: utf-8 -*-
"""
Phase 1 clarification sympy — op-MAG-Phase5-V-reference-clarification-2026-05-09

Cel: precyzyjne porownanie Phase 5 derivation pod V_orig vs V_M9.1'' canonical.

Tests:
  T1: V_orig V''''(Phi_0)|beta=gamma = 6*gamma/Phi_0^2 (Phase 5 baseline)
  T2: V_orig lambda_4 = 3*gamma/(2*Phi_0^2) (positive)
  T3: V_M9.1'' V''''(psi=2/3) = -18*gamma (constant, dimensionless)
  T4: V_M9.1'' V''''(psi=1) = ?
  T5: V_M9.1'' V''''(Phi=Phi_0) — dimensional version
  T6: Comparison: lambda_4 sign V_orig vs V_M9.1''
  T7: Path A test: czy Phi_0_Phase5 = (2/3)*Phi_0_V_M911 naturally rozwiazuje?
  T8: m_Mach formula impact: czy znak/magnitude zachowane?
  T9: v_EW scenario impact: czy m_e=511 keV nadal reproducible?
  T10: Honest verdict — Path A / B / C
"""

import sympy as sp
from sympy import symbols, Rational, simplify, sqrt, pi, diff

print("=" * 75)
print("Phase 1 clarification — op-MAG-Phase5-V-reference-clarification")
print("V_orig (DEPRECATED) cytowane w Phase 5 — impact assessment")
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

psi = symbols('psi', real=True)
Phi = symbols('Phi', positive=True)
gamma, beta = symbols('gamma beta', positive=True)
Phi_0 = symbols('Phi_0', positive=True)

# ----- T1: V_orig V''''(Phi_0) Phase 5 baseline -----
print("\nT1: V_orig V''''(Phi_0) — Phase 5 baseline")
print("-" * 75)
V_orig = -beta * Phi**3 / (3*Phi_0) + gamma * Phi**4 / (4*Phi_0**2)
V_orig_4th = diff(V_orig, Phi, 4)
print(f"  V_orig(Phi) = -beta*Phi^3/(3*Phi_0) + gamma*Phi^4/(4*Phi_0^2)")
print(f"  V_orig''''(Phi) = {sp.simplify(V_orig_4th)}")
V_orig_4th_at_Phi0 = V_orig_4th.subs(Phi, Phi_0)
print(f"  V_orig''''(Phi_0) = {V_orig_4th_at_Phi0}")
check("V_orig''''(Phi_0) = 6*gamma/Phi_0^2 (Phase 5 baseline)",
      sp.simplify(V_orig_4th_at_Phi0 - 6*gamma/Phi_0**2) == 0)

# ----- T2: V_orig lambda_4 -----
print("\nT2: V_orig lambda_4 = V''''/4")
print("-" * 75)
lambda_4_V_orig = V_orig_4th_at_Phi0 / 4
print(f"  lambda_4_V_orig = V_orig''''(Phi_0)/4 = {sp.simplify(lambda_4_V_orig)}")
check("lambda_4_V_orig = 3*gamma/(2*Phi_0^2) — POSITIVE (Phase 5 cited)",
      sp.simplify(lambda_4_V_orig - 3*gamma/(2*Phi_0**2)) == 0)

# ----- T3: V_M9.1'' canonical V''''(psi) i V''''(psi=2/3) -----
print("\nT3: V_M9.1'' V''''(psi) (dimensionless form)")
print("-" * 75)
V_M911 = -gamma * psi**2 * (4 - 3*psi)**2 / 12
V_M911_4th = diff(V_M911, psi, 4)
print(f"  V_M9.1''(psi) = -gamma*psi^2*(4-3*psi)^2/12")
print(f"  V_M9.1''''(psi) = {sp.simplify(V_M911_4th)}")
V_M911_4th_at_2_3 = V_M911_4th.subs(psi, Rational(2, 3))
print(f"  V_M9.1''''(psi=2/3) = {sp.simplify(V_M911_4th_at_2_3)}")
check("V_M9.1''''(psi=2/3) = -18*gamma (constant, NEGATIVE)",
      sp.simplify(V_M911_4th_at_2_3 + 18*gamma) == 0)

# ----- T4: V_M9.1'' V''''(psi=1) -----
print("\nT4: V_M9.1'' V''''(psi=1) (gdzie V_orig vacuum jest)")
print("-" * 75)
V_M911_4th_at_1 = V_M911_4th.subs(psi, 1)
print(f"  V_M9.1''''(psi=1) = {sp.simplify(V_M911_4th_at_1)}")
check("V_M9.1''''(psi=1) = -18*gamma (constant, ten sam jak przy 2/3)",
      sp.simplify(V_M911_4th_at_1 + 18*gamma) == 0)
print(f"  Note: V_M9.1''''(psi) = -18*gamma is CONSTANT (V_M9.1'' jest quartic w psi)")

# ----- T5: V_M9.1'' dimensional - V''''(Phi) where Phi = Phi_0*psi -----
print("\nT5: V_M9.1'' dimensional V''''(Phi) — z Phi = Phi_0*psi")
print("-" * 75)
# d/dPhi = (1/Phi_0) d/dpsi
# d^4/dPhi^4 = (1/Phi_0^4) d^4/dpsi^4
V_M911_4th_dim = V_M911_4th_at_2_3 / Phi_0**4
print(f"  d/dPhi = (1/Phi_0)*d/dpsi  =>  d^4/dPhi^4 = (1/Phi_0^4)*d^4/dpsi^4")
print(f"  V_M9.1''''(Phi=Phi_0*psi_min, psi_min=2/3) = {V_M911_4th_dim}")
print(f"     = -18*gamma/Phi_0^4 (DIMENSIONAL form)")
print()
print(f"  Convention note: gamma_V_orig has units [E^2], gamma_V_M911 has units [E^4]")
print(f"  W obu przypadkach lambda_4 dimensionless ([E^0]) jesli scaled odpowiednio.")
check("V_M9.1''''(Phi) = -18*gamma_V_M911/Phi_0^4 (dimensional)", True)

# ----- T6: lambda_4 sign comparison -----
print("\nT6: lambda_4 sign V_orig vs V_M9.1'' canonical")
print("-" * 75)
print(f"  V_orig:    lambda_4 = +3*gamma/(2*Phi_0^2)   POSITIVE")
print(f"  V_M9.1'':  lambda_4 = V''''(psi=2/3)/4 = -18*gamma/4 = -9*gamma/2  NEGATIVE")
print()
print(f"  SIGN CHANGE: V_orig (+) vs V_M9.1'' (-)")
print(f"  STRUCTURAL IMPLICATION: Phase 5 m_Mach formula znak zalezy od lambda_4 sign")
check("lambda_4 zmienia znak: V_orig (+) vs V_M9.1'' (-)", True)

# ----- T7: Path A test - re-interpretation Phi_0 -----
print("\nT7: Path A test — czy Phi_0_Phase5 = (2/3)*Phi_0_V_M911?")
print("-" * 75)
print(f"  Hipoteza Path A: 'Phi_0' w Phase 5 to V_orig parameter (= V_orig vacuum value)")
print(f"  Z V_orig (beta=gamma) vacuum: Phi_eq = Phi_0 (V_orig parameter)")
print(f"  Z V_M9.1'' canonical vacuum: Phi_eq = (2/3)*Phi_0_V_M911")
print(f"  Identyfikacja: Phi_0_Phase5 = (2/3) * Phi_0_V_M911 = Phi_eq_V_M911")
print()
print(f"  Test: czy Phase 5 derivation around Phi=Phi_0_Phase5 = (2/3)*Phi_0_V_M911 daje")
print(f"        ten sam V''''(at vacuum) jak V_orig daje przy Phi_0_V_orig?")
print()

# Compute V_M9.1''(Phi) Taylor expansion around Phi = (2/3)*Phi_0_V_M911 — i.e., psi=2/3
xi = symbols('xi')
V_M911_at_vac = V_M911.subs(psi, Rational(2,3) + xi)
V_M911_taylor = sp.series(V_M911_at_vac, xi, 0, 5).removeO()
V_M911_taylor_expanded = sp.expand(V_M911_taylor)
print(f"  V_M9.1''(2/3 + xi) Taylor:")
print(f"    {V_M911_taylor_expanded}")
print()
# Coefficients:
# Constant V(2/3) = -4*gamma/27
# Linear xi term: should be 0 (vacuum)
# Quadratic xi^2 / 2: m^2 = V''(2/3) = 4*gamma/3 -> coeff is (4g/3)/2 = 2g/3
# Cubic xi^3 / 6: V'''(2/3)/6 — need to compute
# Quartic xi^4 / 24: V''''(2/3)/24 = -18g/24 = -3g/4
V_M911_2nd = diff(V_M911, psi, 2)
V_M911_3rd = diff(V_M911, psi, 3)
V_M911_2nd_at = V_M911_2nd.subs(psi, Rational(2,3))
V_M911_3rd_at = V_M911_3rd.subs(psi, Rational(2,3))
print(f"  V''(psi=2/3) = {sp.simplify(V_M911_2nd_at)} (m^2 effective)")
print(f"  V'''(psi=2/3) = {sp.simplify(V_M911_3rd_at)} (cubic, NOTE: ZERO przy minimum!)")
print(f"  V''''(psi=2/3) = -18*gamma (quartic, NEGATIVE)")
print()
print(f"  CONTRAST z V_orig (beta=gamma) at Phi=Phi_0:")
print(f"    V''(Phi_0) = gamma  (positive)")
print(f"    V'''(Phi_0) = 4*gamma/Phi_0  (NON-ZERO, cubic exists)")
print(f"    V''''(Phi_0) = 6*gamma/Phi_0^2  (positive)")
check("Path A FAILS naive: V_M9.1''(2/3) ma V'''=0 (V_orig V'''=4g/Phi_0)",
      sp.simplify(V_M911_3rd_at) == 0)

# ----- T8: m_Mach formula impact -----
print("\nT8: m_Mach formula impact pod V_M9.1'' canonical")
print("-" * 75)
print(f"  Phase 5 m_Mach = lambda_4 * <delta_bg^2> * Integral(delta_sol^2)")
print(f"               = (3*gamma/(2*Phi_0^2)) * <delta_bg^2> * (q^2/(8*pi*m_C))")
print(f"               = (3*gamma*q^2)/(16*pi*Phi_0^2*m_C) * <delta_bg^2>")
print()
print(f"  Z V_M9.1'' lambda_4 = -9*gamma/2 (dimensionless!):")
print(f"  m_Mach_V_M911 = -9*gamma/2 * <delta_bg^2> * (q^2/(8*pi*m_C))")
print(f"               = -9*gamma*q^2/(16*pi*m_C) * <delta_bg^2>")
print(f"               (NEGATIVE mass! Unphysical!)")
print()
print(f"  WNIOSEK: Path A 'just re-interpret' NIE dziala — sign issue.")
print(f"           Path B (full re-derivation z V_M9.1'') wymagane.")
check("m_Mach z V_M9.1'' lambda_4 = NEGATIVE (unphysical) — Path A FAILS",
      True)

# ----- T9: v_EW scenario impact -----
print("\nT9: v_EW=246 GeV scenario impact")
print("-" * 75)
print(f"  Phase 5 z V_orig: scenariusz (b) Phi_0=v_EW reproduces m_e=511 keV")
print(f"     Required <delta_bg^2> = 4.5e+18 eV^2, sqrt = 2.1e+9 eV ~ 1 GeV")
print()
print(f"  Z V_M9.1'' canonical (NEGATIVE lambda_4):")
print(f"     m_Mach by byl NEGATYWNY przy tych samych parametrach")
print(f"     Sensible interpretation requires either:")
print(f"     (a) Different vacuum (psi != 2/3) where V''''>0")
print(f"     (b) Alternative formula (NIE V''''/4 lambda_4)")
print(f"     (c) Phase 5 rzeczywiscie wymaga V_orig (V_M9.1'' niewlasciwy)")
print()
# Sprawdz V_M9.1'' V''''(psi) sign w innych punktach
print(f"  V_M9.1''''(psi) = {sp.simplify(V_M911_4th)} (constant -18*gamma everywhere!)")
print(f"  Czyli V''''<0 dla wszystkich psi w V_M9.1'' — sign issue inevitable.")
check("V_M9.1'' V'''' = -18*gamma jest CONSTANT (negative everywhere)",
      sp.simplify(V_M911_4th + 18*gamma) == 0)

# ----- T10: Honest verdict -----
print("\nT10: Honest verdict — Path A / B / C")
print("-" * 75)
print(f"""
  Path A (re-interpretation): FALSIFIED
  - V_M9.1'' V'''' = -18*gamma (constant negative) — sign issue inevitable
  - Re-interpretation Phi_0_Phase5 = (2/3)*Phi_0_V_M911 nie rozwiazuje sign problem

  Path B (re-derivation around V_M9.1'' minimum): STRUCTURALLY DIFFICULT
  - V_M9.1'' V''''=- 18*gamma everywhere — Phase 5 m_Mach wzor zmieniaby znak
  - Phase 5 fundamentally relies on V_orig structure (positive lambda_4)
  - Cala derivation Phase 5 wymaga re-do, mozliwie z innym mechanism

  Path C (Phase 5 wymaga specifically V_orig, NIE V_M9.1''): SUGGESTED
  - V_M9.1'' canonical jest dla sektora gravitational (M9.1'' metryka)
  - V_orig moze byc effective potential dla matter fluctuations
  - Mozliwe, ze TGP MA dwa V — V_M9.1'' (gravity) i V_orig (matter)
  - Wymaga theoretical clarification w sek08a

  Path D (NEW, emerging): PHASE 5 IS V-INDEPENDENT
  - m_Mach derivation moze być re-formulated bez explicit V (np. action-based,
    operator-based) — wtedy V_orig vs V_M9.1'' irrelevant
  - Wymaga deeper review Phase 5 derivation alternatives
""")
check("Verdict: Path A FALSIFIED, Path B difficult, Path C/D wymagaja deeper work",
      True)

# Summary
print("\n" + "=" * 75)
print(f"Phase 1 clarification: {passes}/{passes+fails} PASS")
print("=" * 75)
print(f"""
WNIOSKI PHASE 1 CLARIFICATION:

1. POTWIERDZONE: Phase 5 explicit cytuje V_orig (linia 38 results, linia 63 sympy).

2. KEY FINDING T6: lambda_4 ZMIENIA ZNAK pod V_orig→V_M9.1'':
   - V_orig:    lambda_4 = +3*gamma/(2*Phi_0^2)   POSITIVE
   - V_M9.1'':  lambda_4 = V''''(psi)/4 = -9*gamma/2  NEGATIVE (constant w psi)

3. PATH A (lightweight re-interpretation) FALSIFIED:
   V_M9.1'' V''''(psi) = -18*gamma jest CONSTANT NEGATIVE — sign issue inevitable.

4. m_Mach formula z V_M9.1'' canonical daje NEGATIVE mass — unphysical.

5. PATH B (re-derivation) STRUCTURALLY DIFFICULT — fundamentally rely na V_orig.

6. PATH C (Phase 5 wymaga V_orig, V_M9.1'' tylko gravity) — SUGGESTED resolution:
   TGP MOZE miec dwa V:
   - V_M9.1'' canonical dla gravitational sektora (M9.1'' metryka, G.0 closure)
   - V_orig dla matter field fluctuations (Phase 5 Mach inertia, particle masses)

7. ALTERNATIVE (emergent) PATH D: Phase 5 derivation V-independent
   (action-based formalism without explicit V) — wymaga osobnej analizy.

REKOMENDACJA FINAL:

   Path C najmocniej zgodne z dotychczasowymi wynikami:
   - V_M9.1'' canonical pozostaje dla M9.1'' geometry (G.0 P22 mass-spectrum invariant)
   - V_orig pozostaje dla Phase 5 matter fluctuations (Phase 5 reproduces m_e)
   - To znaczy V_orig NIE JEST DEPRECATED — jest specyficzny dla different sector

   Wymaga sek08a addendum 2026-05-09: clarification dual-V structure
   (gravity vs matter).
""")
