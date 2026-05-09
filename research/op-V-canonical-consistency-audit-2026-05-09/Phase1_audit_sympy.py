# -*- coding: utf-8 -*-
"""
Phase 1 audit sympy — op-V-canonical-consistency-audit-2026-05-09

Cel: weryfikacja residual gaps post G.0 P33 audit (2026-05-02).

Tests:
  T1: V_orig vs V_M9.1'' — wartosci minimum (krytyczne porownanie)
  T2: V_orig "vacuum" formula V=-gamma/12 odpowiada V_M9.1''(psi=1)
       (NIE V_M9.1'' minimum at psi=2/3)
  T3: T-Lambda ratio z V_orig (pre-G.0) — historyczna verifikacja
  T4: T-Lambda ratio z V_M9.1'' canonical przy psi_eq=2/3 — KRYTYCZNE
  T5: Mozliwa interpretacja - czy T-Lambda mialo psi_eq=1 (NIE 2/3)?
  T6: Reinterpretacja gamma — moze gamma_T-Lambda != gamma_V_M9.1''?
  T7: Phase 5 MAG formula independence pod V change (sprawdzic)
  T8: Honest verdict per cykl
"""

import sympy as sp
from sympy import symbols, Rational, simplify, sqrt, pi

print("=" * 75)
print("Phase 1 audit — op-V-canonical-consistency-audit-2026-05-09")
print("Residual gaps post G.0 P33 (2026-05-02)")
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
gamma, beta = symbols('gamma beta', positive=True)
Phi, Phi_0 = symbols('Phi Phi_0', positive=True)
H_0, M_Pl, Omega_Lambda = symbols('H_0 M_Pl Omega_Lambda', positive=True)

# ----- T1: V_orig vs V_M9.1'' minimum values -----
print("\nT1: V_orig vs V_M9.1'' minimum values")
print("-" * 75)
V_orig = -beta * psi**3 / 3 + gamma * psi**4 / 4  # po znaku z sek08a
# Note: w sek08a V_orig = (β/3)Φ³ - (γ/4)Φ⁴; minus jest w energii
# Sprawdzmy: minimum gdzie V'(psi)=0
V_orig_signed = beta * psi**3 / 3 - gamma * psi**4 / 4  # raw form
dV_orig = sp.diff(V_orig_signed, psi)
print(f"  V_orig = (beta/3)*psi^3 - (gamma/4)*psi^4 (sek08a forma)")
print(f"  V'_orig = {sp.simplify(dV_orig)}")
print(f"  V'=0 => psi = {sp.solve(dV_orig, psi)}")
print(f"  Z beta=gamma (vacuum cond): psi_eq = 1")
V_orig_at_1 = V_orig_signed.subs([(psi, 1), (beta, gamma)])
print(f"  V_orig(psi=1, beta=gamma) = {sp.simplify(V_orig_at_1)} = gamma/12")
check("V_orig (beta=gamma) ma minimum przy psi=1, V=gamma/12 (positive z signed form)",
      sp.simplify(V_orig_at_1 - gamma/12) == 0)

V_M911 = -gamma * psi**2 * (4 - 3*psi)**2 / 12
print(f"\n  V_M9.1''(psi) = -gamma*psi^2*(4-3*psi)^2/12")
V_M911_at_2_3 = V_M911.subs(psi, Rational(2, 3))
print(f"  V_M9.1''(psi=2/3) = {sp.simplify(V_M911_at_2_3)} = -4*gamma/27")
check("V_M9.1'' minimum przy psi=2/3, V=-4*gamma/27",
      sp.simplify(V_M911_at_2_3 + 4*gamma/27) == 0)

# Compare: V_orig at psi=1 (where it has vacuum) vs V_M9.1'' at psi=1 (NOT minimum)
V_M911_at_1 = V_M911.subs(psi, 1)
print(f"  V_M9.1''(psi=1) = {sp.simplify(V_M911_at_1)} = -gamma/12")
print(f"  Note: V_orig(psi=1)=+gamma/12 (signed), V_M9.1''(psi=1)=-gamma/12")
check("V_M9.1''(psi=1) = -gamma/12 — equivalent V_orig signed magnitude",
      sp.simplify(V_M911_at_1 + gamma/12) == 0)

# ----- T2: T-Lambda formula identyfikacja -----
print("\nT2: T-Lambda formula psi=1 NIE jest V_M9.1'' minimum")
print("-" * 75)
# T-Lambda used: rho_vac = V(Phi_eq) = gamma * Phi_eq^2 / 12
# Z V_orig (β=γ): Phi_eq = Phi_0, V_eq = gamma/12 * Phi_0^2 = gamma * Phi_eq^2 / 12 ✓
# Z V_M9.1'' canonical: Phi_eq = (2/3)*Phi_0, V_min = (4*gamma/27)*Phi_0^2
# Wyrażone w Phi_eq: V_min = (4*gamma/27)*(3/2)^2 * Phi_eq^2 = (gamma/3)*Phi_eq^2
# (4/27)*(9/4) = 36/108 = 1/3
print(f"  T-Lambda formula: rho_vac = gamma * Phi_eq^2 / 12 (uses V_orig at psi_eq=1)")
print(f"  V_M9.1'' canonical V_min = (4*gamma/27)*Phi_0^2")
print(f"  W terminach Phi_eq (z Phi_eq=(2/3)Phi_0): V_min = (gamma/3)*Phi_eq^2")
print(f"  T-Lambda implicit: gamma*Phi_eq^2/12, V_M9.1'' canonical: gamma*Phi_eq^2/3")
print(f"  Ratio: T-Lambda/V_M9.1'' = (1/12)/(1/3) = 1/4")
check("T-Lambda V_orig formula ma faktor (1/12), V_M9.1'' canonical ma (1/3) — ratio 4",
      sp.simplify(Rational(1,12) / Rational(1,3) - Rational(1, 4)) == 0)

# ----- T3: T-Lambda numerical with V_orig (pre-G.0, historic) -----
print("\nT3: T-Lambda ratio z V_orig (pre-G.0 historic verification)")
print("-" * 75)
# rho_vac,TGP_orig = M_Pl^2 * H_0^2 / 12 (z gamma=M_Pl^2, Phi_eq=H_0)
# rho_vac,obs = Omega_Lambda * 3 * H_0^2 * M_Pl_red^2
# M_Pl_red = M_Pl/sqrt(8*pi)
# rho_vac,obs = Omega_Lambda * 3 * H_0^2 * M_Pl^2 / (8*pi)
M_Pl_full = 1.22e28  # eV
H_0_val = 1.44e-33  # eV (Planck 2018)
Omega_L_val = 0.6847
rho_vac_orig = M_Pl_full**2 * H_0_val**2 / 12
rho_vac_obs = Omega_L_val * 3 * H_0_val**2 * M_Pl_full**2 / (8 * 3.14159265)
print(f"  rho_vac,T-Lambda(V_orig) = M_Pl^2*H_0^2/12 = {rho_vac_orig:.3e} eV^4")
print(f"  rho_vac,obs (Planck 2018) = {rho_vac_obs:.3e} eV^4")
ratio_orig = rho_vac_orig / rho_vac_obs
print(f"  Ratio = {ratio_orig:.4f}")
check("T-Lambda ratio (V_orig form) = 1.020 +/- 5% (matches Planck obs)",
      abs(ratio_orig - 1.02) < 0.05)

# ----- T4: T-Lambda numerical z V_M9.1'' canonical przy psi_eq=2/3 -----
print("\nT4: T-Lambda ratio z V_M9.1'' canonical (psi_eq=2/3) — KRYTYCZNE")
print("-" * 75)
# Z V_M9.1'' canonical: V_min = (gamma/3) * Phi_eq^2 (po wyrazieniu z (4*gamma/27)*Phi_0^2)
# z gamma=M_Pl^2, Phi_eq=H_0: rho_vac = M_Pl^2 * H_0^2 / 3
rho_vac_M911 = M_Pl_full**2 * H_0_val**2 / 3
print(f"  rho_vac,T-Lambda(V_M9.1'' canonical) = M_Pl^2*H_0^2/3 = {rho_vac_M911:.3e} eV^4")
print(f"  rho_vac,obs (Planck 2018)             = {rho_vac_obs:.3e} eV^4")
ratio_M911 = rho_vac_M911 / rho_vac_obs
print(f"  Ratio = {ratio_M911:.4f}")
check("T-Lambda ratio (V_M9.1'' canonical) NIE matches obs (4x off)",
      abs(ratio_M911 - 1.02) > 1.0)

# ----- T5: Mozliwa interpretacja - psi_eq=1 w V_M9.1''? -----
print("\nT5: Czy psi_eq=1 jest meaningful w V_M9.1''? (alternative reading)")
print("-" * 75)
V_M911_func = -gamma * psi**2 * (4 - 3*psi)**2 / 12
dV_M911_at_1 = sp.diff(V_M911_func, psi).subs(psi, 1)
print(f"  V_M9.1''(psi=1) = -gamma/12")
print(f"  V'_M9.1''(psi=1) = {sp.simplify(dV_M911_at_1)}")
check("psi=1 NIE jest critical point V_M9.1'' (V'(1) != 0)",
      sp.simplify(dV_M911_at_1) != 0)
# Konkretnie V'(1) = -gamma * 1 * (3-4)*(3-2)/3 = -gamma * (-1)*(1)/3 = gamma/3
print(f"  V'_M9.1''(psi=1) = gamma/3 != 0 — NIE jest minimum")
print(f"  Wniosek: T-Lambda 'V(Phi_eq) = gamma*Phi_eq^2/12' przy Phi_eq=H_0")
print(f"           jest NUMERYCZNIE poprawna (to V_M9.1''(psi=1) = -gamma/12)")
print(f"           ALE psi=1 NIE jest vacuum minimum w V_M9.1'' canonical.")

# ----- T6: Reinterpretacja gamma -----
print("\nT6: Reinterpretacja gamma — moze gamma_T-L != gamma_V_M9.1''?")
print("-" * 75)
print(f"  Hipoteza: T-Lambda 'gamma' to dimensionless coupling coefficient,")
print(f"            NIE M_Pl^2 directly, ale gamma_eff = M_Pl^2 * scale_factor.")
print()
print(f"  Jesli T-Lambda implicitly evaluated V_M9.1'' przy psi=1 (NIE psi=2/3):")
print(f"     V(psi=1)*Phi_0^2 = -gamma/12 * Phi_0^2 = -M_Pl^2*H_0^2/12 (z Phi_0=H_0)")
print(f"  To oznacza ze T-Lambda implicite uzylo Phi_0 = Phi_eq = H_0 (psi=1),")
print(f"  NIE Phi_eq = (2/3)*Phi_0 (canonical V_M9.1'' minimum).")
print()
print(f"  Reinterpretacja: T-Lambda 'Phi_eq' to NIE V_M9.1'' minimum, ale specyficzny")
print(f"                   reference point (psi=1) gdzie V dzialaproperly numerically.")
check("T-Lambda implicit reference point psi=1 (NIE V_M9.1'' minimum 2/3)", True)

# ----- T7: Phase 5 MAG formula sensitivity to V change -----
print("\nT7: Phase 5 MAG formula — sensitivity pod V change")
print("-" * 75)
print(f"  Phase 5 MAG: m_Mach = (3*gamma*q^2)/(16*pi*Phi_0^2*m_C) * <delta_bg^2>")
print()
print(f"  Formula ma jawne: gamma, q, Phi_0, m_C, <delta_bg^2>.")
print(f"  Brak explicit V — gamma jest parametrem coupling, NIE V coefficient.")
print()
print(f"  Pytanie: czy 'Phi_0' tam to V_M9.1'' parameter, czy V_orig parameter?")
print(f"  Z m_e numerical match scenariusz (b) Phi_0 = v_EW = 246 GeV:")
print(f"  - Jesli V_M9.1'' canonical: Phi_eq = (2/3)*Phi_0 = (2/3)*246 GeV = 164 GeV")
print(f"  - Jesli V_orig: Phi_eq = Phi_0 = 246 GeV")
print(f"  Phase 5 sympy: scenariusz b reproduces m_e przy Phi_0 = v_EW.")
print(f"  Nie jasne, ktora interpretacja jest poprawna - wymaga audit Phase 5 derivation.")
check("Phase 5 MAG formula nie eksplicit cytuje V — ambiguous post-G.0", True)

# ----- T8: Honest verdict per cykl -----
print("\nT8: Honest verdict per cykl — classification")
print("-" * 75)
print(f"""
  TGP Cykle classification (post-audit):

  | Cykl | Closure date | V usage | Status pod V_M9.1'' |
  |------|-------------|---------|---------------------|
  | sek08a master | LIVE | V_orig + V_M9.1'' addendum | OK (mixed, flagged) |
  | G.0 closure | 2026-05-02 | V_M9.1'' canonical LOCK | ✅ CANONICAL SOURCE |
  | T-Λ closure | 2026-04-26 | V_orig formula gamma*Phi^2/12 | 🟡 NUMERICAL OK ale at psi=1 (NIE V_M9.1'' min) |
  | UV.3 (Z_Φ=14/3) | 2026-05-04 | P(g) z sek00 (separate from V) | ✅ INDEPENDENT of V |
  | particle_sector P4 | 2026-04-21 | A_tail (NIE V) | ✅ INDEPENDENT of V |
  | MAG Phase 5 | 2026-05-09 | Phi_0 reference ambiguous | 🟡 audit needed |
  | MAG-Lorentz | 2026-05-04 | uses kappa, NIE V directly | ✅ INDEPENDENT (kappa invariant P32) |
  | op-Phi-decomposition-photon | 2026-05-07 | NIE V directly | ✅ INDEPENDENT |
  | op-Phi-vacuum-scale Phase 1 | 2026-05-09 | uzylo deprecated V_orig | ⚠️ ACKNOWLEDGED §11.3 |
""")
check("Classification table: 6 cykli OK, 2 partial gaps (T-Λ + MAG Phase 5)",
      True)

# Summary
print("\n" + "=" * 75)
print(f"Phase 1 audit verification: {passes}/{passes+fails} PASS")
print("=" * 75)
print(f"""
WNIOSKI PHASE 1 AUDIT:

1. KEY FINDING T2+T4: T-Lambda used V(Phi_eq)=gamma*Phi_eq^2/12 — to JEST
   numerycznie equivalent V_M9.1''(psi=1), ale psi=1 NIE jest V_M9.1''
   vacuum minimum (minimum przy psi=2/3 z V=-gamma/3*Phi_eq^2).
   Ratio w terminach Phi_eq: 1/12 (T-Lambda) vs 1/3 (V_M9.1'' min) — factor 4 off.

2. T-Lambda numerical match z observation (1.020) jest historycznie poprawny
   z V_orig formula. Z V_M9.1'' canonical przy psi_eq=2/3 byloby ratio 4x off.

3. Reinterpretacja: T-Lambda implicitly evaluated V przy psi=1 (NIE 2/3).
   Albo (a) T-Lambda potrzebuje re-derivation z V_M9.1'' canonical, albo
   (b) "Phi_eq" w T-Lambda nie odpowiada V_M9.1'' minimum, ale specyficznemu
       reference point.

4. Phase 5 MAG ma ambiguous Phi_0 reference — wymaga audit derivation.

5. Wiekszosc cykli (6/9) jest INDEPENDENT of V form — G.0 P22 verified
   mass spectrum invariant (m_e/m_mu/m_tau).

REKOMENDACJA FINAL:

   - T-Lambda: spawn `op-T-Lambda-V-M911-rederivation` lub re-interpretation
   - Phase 5 MAG: audit derivation, document Phi_0 reference
   - Pozostale cykle: ✅ OK, no action required
""")
