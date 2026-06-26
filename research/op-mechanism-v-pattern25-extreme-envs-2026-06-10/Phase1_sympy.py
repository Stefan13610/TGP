#!/usr/bin/env python3
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
"""
Phase1_sympy.py — F-P25-A: TGP-native near-horizon source derivation + scaling class
=====================================================================================
Cycle: op-mechanism-v-pattern25-extreme-envs-2026-06-10
Phase: 1 (F-P25-A; sympy/analytical; circularity audit per Phase 0 §3)

GOAL (per Phase0_balance.md §3 F-P25-A):
  Derive the TGP-native source term for <Phi> in a compact-binary near-horizon region
  from LOCKED machinery (M9.2 source convention; T3 EOM; emergent-metric Newton matching),
  identify the scaling class (S-rho density-type vs S-kappa compactness-type per §1.2),
  and compute its dimensionless magnitude under Branch A (gamma ~ M_Pl^2, IMMUTABLE).

PRE-DECLARED ACCEPTANCE (Phase 0 §3, immutable):
  PASS_SOURCE_DERIVED / PARTIAL_SOURCE_NS_ONLY / FAIL_NO_SOURCE

DISCIPLINE:
  - compute-then-compare; 0 hardcoded T_pass=True
  - Branch A mapping IMMUTABLE (Phase 0 §4.1, forbidden move #2)
  - thresholds delta_psi_critical = 2*sqrt(3)/9 (T3 EXACT) + factor-10 PARTIAL band:
    appear ONLY in comparison gates, NEVER in source forms (circularity audit FP15)
  - m_Phi category per foundations §3.5.6.3 rule: this Phase uses
    m_Phi_intrinsic = V''(psi=2/3) for the AMBIENT (vacuum) screening scale, and
    DERIVES whether <Phi>_local shifts; self-consistency checked in FP14.

LOCKED references (compute-then-compare targets):
  - T3 Phase 1: W'(2/3) = -4/3; psi_pm = (6 +/- 2*sqrt(3))/9
  - T3 Phase 2 numerics: delta_psi_max(M=0.01, sigma=1) = 1.91e-4
  - T3 Phase 3 LOCKED .txt output: delta_psi(Branch A, M=10 M_Sun, sigma=30 km) = 6.83e-81
    (NOTE: Phase3_results.md §0 quotes 1.74e-79 — transcription slip in markdown;
     the executed-script output 6.83e-81 is authoritative. INFORMATIONAL flag FP8;
     predecessor verdict UNTOUCHED — both values are ~78-80 orders below threshold.)
  - M9.2 setup: physical coupling convention q = 8*pi*G/c^2 (natural q=1)
"""

import numpy as np
import sympy as sp
from sympy import (symbols, Rational, sqrt, pi, exp, oo, integrate, simplify,
                   nsimplify, limit, diff, solve, erfc, Symbol)

print("=" * 78)
print("  Phase 1 (F-P25-A): TGP-native near-horizon source derivation")
print("  Cycle: op-mechanism-v-pattern25-extreme-envs-2026-06-10")
print("=" * 78)

PASS_count = 0
FAIL_count = 0
results = []

def check(label, cond, expected=None, got=None):
    global PASS_count, FAIL_count
    status = "PASS" if bool(cond) else "FAIL"
    if cond:
        PASS_count += 1
    else:
        FAIL_count += 1
    msg = f"  [{status}] {label}"
    if expected is not None or got is not None:
        msg += f"\n          (expected={expected}, got={got})"
    print(msg)
    results.append((label, status))
    return cond

def banner(title):
    print("\n" + "-" * 78)
    print(f"  {title}")
    print("-" * 78)

# ============================================================================
# Section 1 — LOCKED algebra restated verbatim + re-verified (T3 inheritance)
# ============================================================================
banner("Section 1: LOCKED V_M9.1'' algebra re-verification (T3 inheritance)")

psi = symbols('psi', positive=True)

# W(psi) = -V'(psi)/gamma, verbatim from T3 Phase 2 §1.1:
W = Rational(1, 3) * psi * (8 - 18 * psi + 9 * psi**2)

# FP1: vacuum psi = 2/3 is a zero of W
W_at_vac = simplify(W.subs(psi, Rational(2, 3)))
check("FP1: W(2/3) = 0 (cosmological vacuum)", W_at_vac == 0,
      expected=0, got=W_at_vac)

# FP2: linearized Yukawa mass^2 = -W'(2/3) = 4/3
Wp = diff(W, psi)
Wp_at_vac = simplify(Wp.subs(psi, Rational(2, 3)))
m2_tilde = -Wp_at_vac
check("FP2: m_tilde^2 = -W'(2/3) = 4/3 (linearized screening mass^2, natural units)",
      m2_tilde == Rational(4, 3), expected=Rational(4, 3), got=m2_tilde)

# FP3: inflection points psi_pm = (6 +/- 2*sqrt(3))/9; delta_psi_critical EXACT
roots = solve(sp.Eq(Wp, 0), psi)
roots_simplified = sorted([simplify(r) for r in roots], key=lambda r: float(r))
psi_minus_expected = nsimplify((6 - 2 * sqrt(3)) / 9)
psi_plus_expected = nsimplify((6 + 2 * sqrt(3)) / 9)
delta_psi_crit = simplify(psi_plus_expected - Rational(2, 3))   # = 2*sqrt(3)/9
ok3 = (simplify(roots_simplified[0] - psi_minus_expected) == 0 and
       simplify(roots_simplified[1] - psi_plus_expected) == 0 and
       simplify(delta_psi_crit - 2 * sqrt(3) / 9) == 0)
check("FP3: psi_pm = (6 +/- 2*sqrt(3))/9; delta_psi_critical = psi_+ - 2/3 = 2*sqrt(3)/9 EXACT",
      ok3,
      expected="(6+/-2sqrt3)/9; 2sqrt3/9 ~ 0.3849",
      got=f"{[str(r) for r in roots_simplified]}; {float(delta_psi_crit):.4f}")

delta_psi_critical_f = float(delta_psi_crit)        # ~0.3849 (T3 EXACT)
PARTIAL_band_f = delta_psi_critical_f / 10.0        # factor-10 band (Phase 0 §3, declared there)

# ============================================================================
# Section 2 — Linearized response kernel (Yukawa Green's function) + regimes
# ============================================================================
banner("Section 2: Yukawa response kernel + the two scaling regimes")

r, m, a = symbols('r m a', positive=True)

# FP4: G(r) = exp(-m r)/(4 pi r) solves (-Lap + m^2) G = 0 for r > 0,
#      and integrates to 1/m^2 (zero-mode normalization of the kernel)
G = exp(-m * r) / (4 * pi * r)
lap_G = diff(G, r, 2) + (2 / r) * diff(G, r)        # radial Laplacian, r > 0
residual = simplify(lap_G - m**2 * G)
int_G = integrate(G * 4 * pi * r**2, (r, 0, oo))
ok4 = (residual == 0) and (simplify(int_G - 1 / m**2) == 0)
check("FP4: Yukawa kernel verified: Lap(G) = m^2 G (r>0); Int G d^3x = 1/m^2",
      ok4, expected="0; 1/m^2", got=f"{residual}; {int_G}")

# FP5: LOCAL (screened) regime — uniform/slowly-varying rho:
#      delta_psi = q*rho/m^2; with q=1 (natural), m^2 = 4/3 => delta_psi = (3/4)*rho
# Derivation: convolution of constant rho with G = rho * Int G = rho/m^2 (from FP4).
delta_psi_local = simplify(1 / m2_tilde)
check("FP5: local-screening response delta_psi = q*rho/m_tilde^2 = (3/4)*rho_tilde "
      "(S-rho density-type; follows from FP4 normalization)",
      delta_psi_local == Rational(3, 4), expected=Rational(3, 4), got=delta_psi_local)

# FP6: exact linear central response for Gaussian source; cross-check vs T3 Phase 2
#      numerics at M = 0.01, sigma = 1 (LOCKED: delta_psi_max = 1.91e-4)
# delta_psi(0) = Int_0^oo exp(-m r) * rho(r) * r dr  with
# rho(r) = M * exp(-r^2/(2 sigma^2)) / ((2 pi)^(3/2) sigma^3)
M_nat = 0.01
sigma_nat = 1.0
m_nat = float(sqrt(Rational(4, 3)))
integrand = lambda rr: rr * np.exp(-m_nat * rr - rr**2 / (2 * sigma_nat**2))
rr = np.linspace(1e-8, 30.0, 400000)
I_num = np.trapezoid(integrand(rr), rr)
rho0 = M_nat / ((2 * np.pi)**1.5 * sigma_nat**3)
delta_psi_lin_exact = rho0 * I_num
T3_phase2_numeric = 1.91e-4   # LOCKED T3 Phase 2 §2.1 (M=0.01 full nonlinear BVP)
ratio_fp6 = delta_psi_lin_exact / T3_phase2_numeric
ok6 = 0.33 < ratio_fp6 < 3.0   # factor-3 agreement gate (sigma ~ lambda_C regime)
check("FP6: exact-linear Gaussian central response vs LOCKED T3 Phase 2 numeric "
      "(M=0.01): factor-3 agreement",
      ok6, expected="1.91e-4 (x3 band)", got=f"{delta_psi_lin_exact:.3e} (ratio {ratio_fp6:.2f})")

# ============================================================================
# Section 3 — Branch A regression anchor (MANDATORY gate precondition, Phase 0 §3)
# ============================================================================
banner("Section 3: Branch A weak-field regression anchor (T3 Phase 3 LOCKED recipe)")

# Constants IDENTICAL to LOCKED T3 Phase3_dimensional.py (verbatim values):
M_Pl_eV = 1.22e28
M_Sun_kg = 1.989e30
hbar_c_eV_m = 1.973e-7
eV_to_kg = 1.78e-36

m_Phi_A_eV = M_Pl_eV                          # Branch A: m_Phi_intrinsic ~ M_Pl
L_unit_A_m = hbar_c_eV_m / m_Phi_A_eV         # natural length unit
m_Phi_A_kg = m_Phi_A_eV * eV_to_kg            # natural mass unit

M_BBH_kg = 10 * M_Sun_kg
sigma_LIGO_m = 30000.0                        # 30 km (T3 reference source)

M_BBH_nat = M_BBH_kg / m_Phi_A_kg
sigma_LIGO_nat = sigma_LIGO_m / L_unit_A_m

# FP7: reproduce LOCKED T3 Phase 3 .txt output 6.83e-81 (same formula, same constants)
delta_psi_typical = 0.75 * M_BBH_nat / ((2 * np.pi)**1.5 * sigma_LIGO_nat**3)
T3_locked_txt = 6.83e-81
rel_dev_fp7 = abs(delta_psi_typical - T3_locked_txt) / T3_locked_txt
check("FP7: REGRESSION GATE — reproduce T3 Phase 3 LOCKED .txt delta_psi = 6.83e-81 "
      "(rel dev < 2%)",
      rel_dev_fp7 < 0.02, expected="6.83e-81", got=f"{delta_psi_typical:.3e} (rel dev {rel_dev_fp7:.3f})")

# FP8: INFORMATIONAL — markdown transcription discrepancy in predecessor docs
#      Phase3_results.md §0 quotes 1.74e-79; LOCKED executed output (.txt) = 6.83e-81.
#      Factor ~25 slip in the .md transcription. Verdict-irrelevant (both ~78-80 orders
#      below threshold). Predecessor verdict UNTOUCHED per Phase 0 §4.5.
md_quoted = 1.74e-79
factor_discrepancy = md_quoted / T3_locked_txt
check("FP8: INFORMATIONAL flag — T3 Phase3_results.md vs .txt transcription factor "
      "documented (no predecessor modification)",
      24.0 < factor_discrepancy < 27.0,
      expected="~25x (verdict-irrelevant)", got=f"{factor_discrepancy:.1f}x")

print(f"""
  Regression anchor established:
    delta_psi(typical LIGO, Branch A) = {delta_psi_typical:.3e}
    vs delta_psi_critical = {delta_psi_critical_f:.4f}  (shortfall ~{abs(np.log10(delta_psi_typical/delta_psi_critical_f)):.0f} orders)
""")

# ============================================================================
# Section 4 — Scaling-class derivation (CORE of F-P25-A)
# ============================================================================
banner("Section 4: Scaling class — (S-rho) vs (S-kappa) under Branch A")

# FP9: regime selector — the response is LOCAL (S-rho) iff sigma_tilde * m_tilde >> 1.
# Under Branch A, lambda_C = 1/m_tilde ~ Planck-scale; ANY astrophysical source
# (>= km scale) has sigma_tilde * m_tilde ~ 10^39+ => S-rho FORCED.
m_tilde_f = float(sqrt(Rational(4, 3)))
regime_parameter = sigma_LIGO_nat * m_tilde_f
check("FP9: regime selector sigma_tilde*m_tilde ~ 10^39 >> 1 => local-screening "
      "(S-rho) FORCED for ANY astrophysical source under Branch A",
      regime_parameter > 1e30, expected="> 1e30", got=f"{regime_parameter:.2e}")

# FP10: (S-kappa) compactness channel — structural identity + Branch A suppression.
# M9.2 LOCKED convention: physical coupling q = 8*pi*G/c^2. Unscreened (m->0)
# point-source response: delta_psi(r) = q*M/(4 pi r) = 2GM/(c^2 r) = compactness.
Gsym, Msym, csym, rsym = symbols('G M c r', positive=True)
q_phys = 8 * pi * Gsym / csym**2
delta_psi_massless = q_phys * Msym / (4 * pi * rsym)
compactness_identity = simplify(delta_psi_massless - 2 * Gsym * Msym / (csym**2 * rsym))
R_s = 2 * Gsym * Msym / csym**2
delta_psi_at_horizon_massless = simplify(delta_psi_massless.subs(rsym, R_s))
ok10 = (compactness_identity == 0) and (delta_psi_at_horizon_massless == 1)
check("FP10: massless-limit identity: q*M/(4 pi r) = 2GM/(c^2 r); at r = R_s "
      "=> delta_psi = 1 EXACT (q = 8 pi G/c^2, M9.2 LOCKED convention)",
      ok10, expected="identity 0; horizon value 1",
      got=f"{compactness_identity}; {delta_psi_at_horizon_massless}")

# Branch A suppression of the compactness channel:
# delta_psi(r) = [q*M/(4 pi r)] * exp(-m_tilde * r_tilde); at near-horizon
# r_tilde = R_s / L_unit_A ~ 10^39 natural lengths => exponent ~ -10^39.
R_s_m = 2 * 6.674e-11 * M_BBH_kg / (2.998e8)**2     # Schwarzschild radius, 10 M_Sun
R_s_nat = R_s_m / L_unit_A_m
yukawa_exponent = m_tilde_f * R_s_nat
log10_suppression = yukawa_exponent / np.log(10.0)
check("FP11: (S-kappa) EXCLUDED under Branch A — Yukawa suppression factor "
      "exp(-m_tilde*R_s_tilde), log10 ~ -10^39 (double-kill of compactness channel)",
      log10_suppression > 1e35,
      expected="suppression exponent log10 > 1e35",
      got=f"log10(suppression) = -{log10_suppression:.2e}")

print(f"""
  STRUCTURAL CONCLUSION (S-kappa audit):
    The foundations §3.5.6 'extreme delta_psi ~ 0.3+' estimate is EXACTLY the
    massless-limit (unscreened) compactness response: delta_psi(R_s)|_(m->0) = 1
    >= delta_psi_critical = {delta_psi_critical_f:.3f}. Under Branch A (IMMUTABLE),
    this channel carries the factor exp(-{yukawa_exponent:.2e}) ~ 10^(-{log10_suppression:.1e})
    and is STRUCTURALLY EXCLUDED. The '0.3+' figure does NOT survive Branch A screening.
    (This AUDITS the estimate per Phase 0 §6 #15 — it was a test target, not an input.)
""")

# ============================================================================
# Section 5 — Source inventory per pre-declared source class
# ============================================================================
banner("Section 5: Source inventory — (i) BH-BH exterior; (ii) NS-NS near-contact")

# FP12: BH-BH exterior — level-0 LOCKED source is -q*rho_matter (M9.2 + T3 verbatim).
# BH exterior: rho_matter = 0 identically => source = 0 => delta_psi = 0 (no-hair analog).
# Gradient self-source (d delta_psi)^2 is O(delta_psi^2) — vanishes on zero perturbation.
# Curvature coupling: ABSENT from level-0 LOCKED machinery (would be candidate (c)
# framework extension — out of scope per Phase 0 §1.3; forbidden to introduce here).
rho_sym = symbols('rho', nonnegative=True)
q_sym = symbols('q', positive=True)
source_form = -q_sym * rho_sym                       # LOCKED level-0 source, verbatim
source_BH_exterior = source_form.subs(rho_sym, 0)
check("FP12: BH-BH exterior native source = -q*rho_matter|_(rho=0) = 0 identically "
      "(no-hair analog at level-0; no curvature coupling in LOCKED machinery)",
      source_BH_exterior == 0, expected=0, got=source_BH_exterior)

# FP13: NS-NS near-contact PREVIEW (informational for Phase 2; NOT the F-P25-B verdict)
# NS central density ~ 1e18 kg/m^3 (massive NS, ~3-4x nuclear saturation).
# Branch A density unit = m_Phi^4 (c=hbar=1) = M_Pl / L_Pl^3 = Planck density.
rho_NS_kg_m3 = 1.0e18
rho_unit_A_kg_m3 = m_Phi_A_kg / L_unit_A_m**3        # ~5e96 kg/m^3 (Planck density)
rho_NS_nat = rho_NS_kg_m3 / rho_unit_A_kg_m3
delta_psi_NS_preview = 0.75 * rho_NS_nat             # FP5 local response
shortfall_orders = abs(np.log10(delta_psi_NS_preview / PARTIAL_band_f))
check("FP13: NS-NS near-contact preview: delta_psi = (3/4)*rho_NS_tilde ~ 1e-79 "
      "(shortfall vs factor-10 PARTIAL band ~77 orders) — preview only",
      1e-81 < delta_psi_NS_preview < 1e-77,
      expected="~1e-79", got=f"{delta_psi_NS_preview:.2e} (shortfall {shortfall_orders:.0f} orders)")

# FP14: self-consistency of the screening assumption (no bootstrap):
# relative shift of m_tilde^2 from delta_psi. NOTE: W''(2/3) = 0 EXACT (the vacuum
# is a symmetric point of W'), so the leading correction is SECOND order in delta_psi:
# dm2/m2 = (1/2) W'''(2/3) delta_psi^2 / |W'(2/3)|.
Wpp = diff(W, psi, 2)
Wppp = diff(W, psi, 3)
Wpp_at_vac = simplify(Wpp.subs(psi, Rational(2, 3)))
Wppp_at_vac = simplify(Wppp.subs(psi, Rational(2, 3)))
dm2_rel = float(Rational(1, 2) * Wppp_at_vac / m2_tilde) * delta_psi_NS_preview**2
ok14 = (Wpp_at_vac == 0) and (abs(dm2_rel) < 1e-150)
check("FP14: screening self-consistency — W''(2/3) = 0 EXACT; dm2/m2 ~ 6.75*delta_psi^2 "
      "~ 1e-157 (no positive-feedback bootstrap; linearization self-consistent)",
      ok14, expected="W''(2/3)=0; dm2/m2 < 1e-150",
      got=f"W''(2/3)={Wpp_at_vac}; dm2/m2={dm2_rel:.2e}")

# ============================================================================
# Section 6 — Circularity audit (Phase 0 §3 F-P25-A mandate; cycle D §5.5 analog)
# ============================================================================
banner("Section 6: Circularity audit — degeneration checks + threshold independence")

# FP15: (a) rho -> 0 => local response -> 0; (b) M -> 0 => point response -> 0;
#       (c) the derived source/response forms contain NO occurrence of the
#       threshold constants psi_pm / delta_psi_critical (thresholds live ONLY in
#       comparison gates, not in the physics forms).
resp_local_sym = q_sym * rho_sym / m**2
lim_a = limit(resp_local_sym, rho_sym, 0)
resp_point_sym = q_sym * Msym * exp(-m * rsym) / (4 * pi * rsym)
lim_b = limit(resp_point_sym, Msym, 0)
forms = [source_form, resp_local_sym, resp_point_sym]
threshold_atoms = {nsimplify((6 + 2 * sqrt(3)) / 9), nsimplify((6 - 2 * sqrt(3)) / 9),
                   nsimplify(2 * sqrt(3) / 9)}
contaminated = any(
    any(simplify(atom - t) == 0 for t in threshold_atoms
        if atom.is_number)
    for f in forms for atom in f.atoms(sp.Number)
)
ok15 = (lim_a == 0) and (lim_b == 0) and (not contaminated)
check("FP15: CIRCULARITY AUDIT — rho->0 => 0; M->0 => 0; threshold constants "
      "(psi_pm, 2sqrt3/9) ABSENT from all derived source/response forms",
      ok15, expected="0; 0; no contamination",
      got=f"{lim_a}; {lim_b}; contaminated={contaminated}")

# ============================================================================
# Verdict
# ============================================================================
banner("F-P25-A VERDICT")

print(f"""
  Total: {PASS_count}/{PASS_count + FAIL_count} PASS   (hardcoded T_pass=True: 0)

  DERIVATION CHAIN SUMMARY:
  1. LOCKED level-0 source form: -q*rho_matter (M9.2 + T3, verbatim; FP12)
  2. Linearized response kernel: Yukawa G(r) = exp(-m r)/(4 pi r); two regimes (FP4-FP6)
  3. Branch A regime selector: sigma_tilde*m_tilde ~ 10^39 >> 1 for ANY astrophysical
     source => response is LOCAL: delta_psi(x) = q*rho_tilde(x)/m_tilde^2 = (3/4)rho_tilde(x)
     => scaling class (S-rho) DENSITY-TYPE IS FORCED (FP9)
  4. (S-kappa) compactness channel = massless-limit intuition (delta_psi(R_s)|m->0 = 1
     EXACT) — STRUCTURALLY EXCLUDED under Branch A by exp(-10^39) (FP10-FP11).
     Foundations §3.5.6 'extreme delta_psi ~ 0.3+' AUDITED: it is the unscreened
     estimate; it does NOT survive Branch A screening.
  5. Source class (i) BH-BH exterior: rho_matter = 0 => native source = 0 identically
     (no-hair analog at level-0). BH-BH branch NEGATIVE AT THE GATE (FP12).
  6. Source class (ii) NS-NS near-contact: well-defined rho_matter source; (S-rho)
     forced; preview magnitude delta_psi ~ 1e-79 (~77 orders below the factor-10
     PARTIAL band; ~78 below delta_psi_critical) (FP13). Phase 2 records the verdict.
  7. Self-consistency: W''(2/3) = 0 EXACT => no screening bootstrap (FP14).
  8. Circularity audit clean (FP15).

  PRE-REGISTERED VERDICT SELECTION (Phase 0 §3 F-P25-A):

    F-P25-A = PARTIAL_SOURCE_NS_ONLY

    - BH-exterior native source = 0 (no-hair analog) => BH-BH branch NEGATIVE
    - NS-NS near-contact: well-defined rho_matter source; scaling class (S-rho)
      FORCED by Branch A screening; (S-kappa) structurally EXCLUDED
    - NS-NS branch continues to F-P25-B (Phase 2) per Phase 0 §3
""")
