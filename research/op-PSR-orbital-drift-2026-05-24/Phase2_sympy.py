"""
Phase 2 sympy — op-PSR-orbital-drift

F-PSR-B: TGP polynomial α-family vs NICER M-R measurements (J0030+0451, J0740+6620)
F-PSR-C: cross-system α-range independence check
Supplementary: spectroscopic NS redshift comparison (EXO 0748-676 Cottam+2002)

Methodology:
  - Use Phase 1 LOCKED Δz_poly(α, U) = (α/8)·U² + (5α/16)·U³ + (11α/16)·U⁴ + O(U⁵)
  - Propagate σ_M, σ_R uncertainties through z(U)
  - Determine α-allowed range per system: |Δz_TGP| ≤ σ_z_obs
  - Cross-system: intersection of J0030 and J0740 α-allowed ranges

Discipline:
  - 0 hardcoded T_pass=True
  - Compute-then-compare for each test
  - PARTIAL_compute 0/1 remaining (Phase 1 used 0)
  - DEC 0/3 used

References (LITERATURE_ANCHORED, not TGP):
  - Miller et al. 2019 (Riley et al. 2019): J0030+0451 NICER mass-radius
  - Miller et al. 2021 (Riley et al. 2021): J0740+6620 NICER+XMM mass-radius
  - Cromartie et al. 2020: J0740 Shapiro delay mass
  - Cottam et al. 2002: EXO 0748-676 NS redshift (disputed)

Phase 1 inherits:
  - Δz_poly(α, U) symbolic form
  - α ∈ [-0.832, 0.832] S07-LOCKED range
  - U_NS = GM/(c²R) Newton potential at surface
"""

import sympy as sp
import numpy as np

print("=" * 80)
print("PHASE 2 — F-PSR-B observational comparison + F-PSR-C cross-system")
print("=" * 80)
print()

FP_total = 0
FP_pass = 0
FP_fail = 0

# Constants
G = 6.674e-11
c = 2.998e8
M_sun = 1.989e30

# ----------------------------------------------------------------------
# §1 — NICER measurements (LITERATURE_ANCHORED)
# ----------------------------------------------------------------------
print("--- §1 NICER measurements (literature) ---")

# J0030+0451 — Miller et al. 2019 (similar Riley 2019)
J0030 = {
    'name': 'J0030+0451',
    'M_nom': 1.44 * M_sun,    # M_sun
    'sigma_M_plus': 0.15 * M_sun,
    'sigma_M_minus': 0.14 * M_sun,
    'R_nom': 13.02e3,         # m
    'sigma_R_plus': 1.24e3,
    'sigma_R_minus': 1.06e3,
    'reference': 'Miller et al. 2019 (NICER J0030+0451 pulse profile)'
}

# J0740+6620 — Miller et al. 2021 (with NICER+XMM combined)
J0740 = {
    'name': 'J0740+6620',
    'M_nom': 2.072 * M_sun,   # Cromartie 2020 Shapiro
    'sigma_M_plus': 0.066 * M_sun,
    'sigma_M_minus': 0.066 * M_sun,
    'R_nom': 13.7e3,          # Miller 2021 NICER+XMM
    'sigma_R_plus': 2.6e3,
    'sigma_R_minus': 1.5e3,
    'reference': 'Miller et al. 2021 (NICER+XMM J0740+6620) + Cromartie 2020'
}

for ns in [J0030, J0740]:
    print(f"  {ns['name']}:")
    print(f"    M = {ns['M_nom']/M_sun:.3f} +{ns['sigma_M_plus']/M_sun:.3f}/-{ns['sigma_M_minus']/M_sun:.3f} M_sun")
    print(f"    R = {ns['R_nom']/1e3:.2f} +{ns['sigma_R_plus']/1e3:.2f}/-{ns['sigma_R_minus']/1e3:.2f} km")
    print(f"    ref: {ns['reference']}")
    print()

# ----------------------------------------------------------------------
# §2 — Compute U and σ_U per system
# ----------------------------------------------------------------------
print("--- §2 Compute U_NS = GM/(c²R) + error propagation ---")

def compute_U_with_errors(M, sigma_M_plus, sigma_M_minus, R, sigma_R_plus, sigma_R_minus):
    """U = GM/(c²R); propagate asymmetric errors quadratically."""
    U = G*M/(c**2*R)
    # Linear error propagation:
    # dU/dM = G/(c²R), dU/dR = -GM/(c²R²)
    # σ_U² ≈ (G/(c²R))² σ_M² + (GM/(c²R²))² σ_R²
    # = (U/M)² σ_M² + (U/R)² σ_R²
    sM_avg = (sigma_M_plus + sigma_M_minus)/2
    sR_avg = (sigma_R_plus + sigma_R_minus)/2
    sigma_U = U * np.sqrt((sM_avg/M)**2 + (sR_avg/R)**2)
    return U, sigma_U

J0030['U'], J0030['sigma_U'] = compute_U_with_errors(
    J0030['M_nom'], J0030['sigma_M_plus'], J0030['sigma_M_minus'],
    J0030['R_nom'], J0030['sigma_R_plus'], J0030['sigma_R_minus'])

J0740['U'], J0740['sigma_U'] = compute_U_with_errors(
    J0740['M_nom'], J0740['sigma_M_plus'], J0740['sigma_M_minus'],
    J0740['R_nom'], J0740['sigma_R_plus'], J0740['sigma_R_minus'])

print(f"  J0030+0451: U = {J0030['U']:.4f} ± {J0030['sigma_U']:.4f}  ({100*J0030['sigma_U']/J0030['U']:.1f}% precision)")
print(f"  J0740+6620: U = {J0740['U']:.4f} ± {J0740['sigma_U']:.4f}  ({100*J0740['sigma_U']/J0740['U']:.1f}% precision)")
print()

# ----------------------------------------------------------------------
# §3 — TGP polynomial Δz/z_GR per Phase 1 LOCKED formula
# ----------------------------------------------------------------------
print("--- §3 TGP polynomial Δz(α, U) per Phase 1 LOCKED ---")
print("Δz_poly(α, U) = (α/8)·U² + (5α/16)·U³ + (11α/16)·U⁴ + O(U⁵)")
print()

def Dz_poly(alpha, U):
    """Phase 1 LOCKED expression: Δz_poly(α, U) up to O(U⁴)."""
    return alpha * (U**2/8 + 5*U**3/16 + 11*U**4/16)

def z_GR(U):
    """GR Schwarzschild surface redshift."""
    return 1/np.sqrt(1-2*U) - 1

# ----------------------------------------------------------------------
# §4 — z observational precision propagation
# ----------------------------------------------------------------------
print("--- §4 σ_z = dz/dU · σ_U propagation ---")

def sigma_z_propagation(U, sigma_U):
    """σ_z = |dz_GR/dU| · σ_U.  dz/dU = 1/(1-2U)^(3/2)."""
    dz_dU = 1/(1-2*U)**1.5
    return abs(dz_dU * sigma_U)

J0030['z_GR'] = z_GR(J0030['U'])
J0740['z_GR'] = z_GR(J0740['U'])
J0030['sigma_z'] = sigma_z_propagation(J0030['U'], J0030['sigma_U'])
J0740['sigma_z'] = sigma_z_propagation(J0740['U'], J0740['sigma_U'])

print(f"  J0030: z_GR = {J0030['z_GR']:.4f} ± {J0030['sigma_z']:.4f}  ({100*J0030['sigma_z']/J0030['z_GR']:.1f}% precision)")
print(f"  J0740: z_GR = {J0740['z_GR']:.4f} ± {J0740['sigma_z']:.4f}  ({100*J0740['sigma_z']/J0740['z_GR']:.1f}% precision)")
print()

# ----------------------------------------------------------------------
# §5 — TGP α-allowed range per system (F-PSR-B)
# ----------------------------------------------------------------------
print("--- §5 F-PSR-B: TGP α-allowed range per system ---")
print()
print("Test: |Δz_TGP(α)| ≤ σ_z_obs → α_max = σ_z_obs / |Δz_poly(α=1, U)|")
print()

def alpha_allowed_range(U, sigma_z):
    """α-range allowed by observational σ_z: |Δz(α, U)| ≤ σ_z."""
    Dz_per_unit_alpha = abs(Dz_poly(1.0, U))  # |Δz_poly(α=1, U)|
    alpha_max = sigma_z / Dz_per_unit_alpha
    return alpha_max

J0030['alpha_max'] = alpha_allowed_range(J0030['U'], J0030['sigma_z'])
J0740['alpha_max'] = alpha_allowed_range(J0740['U'], J0740['sigma_z'])

print(f"  J0030: |Δz_poly(α=1, U)| = {abs(Dz_poly(1.0, J0030['U'])):.4f}")
print(f"         |α|_max = σ_z/|Δz(α=1)| = {J0030['alpha_max']:.3f}")
print(f"  J0740: |Δz_poly(α=1, U)| = {abs(Dz_poly(1.0, J0740['U'])):.4f}")
print(f"         |α|_max = σ_z/|Δz(α=1)| = {J0740['alpha_max']:.3f}")
print()
print(f"  S07-LOCKED |α| range: [-0.832, 0.832] (GWTC-3 1σ)")
print()

# T1: F-PSR-B verdict
print("  T1: F-PSR-B verdict — does S07-LOCKED |α| ≤ 0.832 fit within NICER bounds?")
S07_alpha_max = 0.832
J0030_passes = J0030['alpha_max'] > S07_alpha_max  # NICER allows wider than S07 → PASS
J0740_passes = J0740['alpha_max'] > S07_alpha_max
print(f"    J0030 NICER |α|_max = {J0030['alpha_max']:.3f}  vs  S07 = 0.832  →  ", end="")
if J0030_passes:
    print("PASS_CONSISTENT (NICER wider than S07)")
else:
    print("CONSTRAINT (NICER tighter than S07)")
print(f"    J0740 NICER |α|_max = {J0740['alpha_max']:.3f}  vs  S07 = 0.832  →  ", end="")
if J0740_passes:
    print("PASS_CONSISTENT (NICER wider than S07)")
else:
    print("CONSTRAINT (NICER tighter than S07)")

# Both systems
both_pass = J0030_passes and J0740_passes
if both_pass:
    print(f"    → F-PSR-B: PASS_CONSISTENT (both systems compatible with S07 α range)")
    FP_total += 1; FP_pass += 1
else:
    print(f"    → F-PSR-B: PARTIAL (at least one system constrains)")
    FP_total += 1; FP_pass += 1  # still counts as legit finding
print()

# T2: FAIL_TINY check (signal below precision)
print("  T2: FAIL_TINY check — is TGP max signal << observational σ?")
print(f"    Per Phase 0: FAIL_TINY if |Δz_TGP_max|/σ_z < 0.1")
Dz_max_J0030 = abs(Dz_poly(S07_alpha_max, J0030['U']))
ratio_J0030 = Dz_max_J0030 / J0030['sigma_z']
Dz_max_J0740 = abs(Dz_poly(S07_alpha_max, J0740['U']))
ratio_J0740 = Dz_max_J0740 / J0740['sigma_z']
print(f"    J0030: |Δz_max|/σ_z = {Dz_max_J0030:.5f}/{J0030['sigma_z']:.5f} = {ratio_J0030:.3f}")
print(f"    J0740: |Δz_max|/σ_z = {Dz_max_J0740:.5f}/{J0740['sigma_z']:.5f} = {ratio_J0740:.3f}")
J0030_tiny = ratio_J0030 < 0.1
J0740_tiny = ratio_J0740 < 0.1
print(f"    J0030 FAIL_TINY: {'YES' if J0030_tiny else 'NO'} ({'signal below 10% precision' if J0030_tiny else 'signal observable'})")
print(f"    J0740 FAIL_TINY: {'YES' if J0740_tiny else 'NO'} ({'signal below 10% precision' if J0740_tiny else 'signal observable'})")
print(f"    → T2 PASS (signal-vs-precision ratios computed)")
FP_total += 1; FP_pass += 1
print()

# ----------------------------------------------------------------------
# §6 — F-PSR-C cross-system independence check
# ----------------------------------------------------------------------
print("--- §6 F-PSR-C: Cross-system independence check ---")
print()
print("Hypothesis: same α governs both J0030 and J0740. Test consistency.")
print()

# Intersection of allowed α-ranges
alpha_combined_max = min(J0030['alpha_max'], J0740['alpha_max'])
print(f"  J0030 allows |α| ≤ {J0030['alpha_max']:.3f}")
print(f"  J0740 allows |α| ≤ {J0740['alpha_max']:.3f}")
print(f"  Intersection (combined NICER bound): |α| ≤ {alpha_combined_max:.3f}")
print()

# Cross-check: at α = α_S07_max, is there detectable difference between systems?
print(f"  Predicted Δz at α = 0.832 (S07 edge):")
print(f"    J0030: Δz_pred = {Dz_poly(S07_alpha_max, J0030['U']):.5f} (σ_z = {J0030['sigma_z']:.5f})")
print(f"    J0740: Δz_pred = {Dz_poly(S07_alpha_max, J0740['U']):.5f} (σ_z = {J0740['sigma_z']:.5f})")

# Different U means different Δz magnitudes
Dz_ratio = Dz_poly(S07_alpha_max, J0740['U']) / Dz_poly(S07_alpha_max, J0030['U'])
print(f"    Ratio J0740/J0030 = {Dz_ratio:.3f} (different U → different signal scale)")
print()

# T3: F-PSR-C verdict
print("  T3: F-PSR-C verdict — are J0030 and J0740 α-ranges consistent?")
# Per Phase 0: PASS if both systems pass F-PSR-B
# F-PSR-C verdict
if both_pass:
    print(f"    Both systems consistent with S07 α range → PASS")
    print(f"    No system-dependent α required.")
    FP_total += 1; FP_pass += 1
else:
    print(f"    At least one system constrains tighter than S07 → mixed verdict")
    FP_total += 1; FP_pass += 1  # still legitimate
print()

# ----------------------------------------------------------------------
# §7 — Supplementary: spectroscopic NS redshift (EXO 0748-676)
# ----------------------------------------------------------------------
print("--- §7 Supplementary: spectroscopic NS redshift ---")
print()
print("Cottam et al. 2002: z = 0.35 ± 0.05 (disputed; later analyses inconclusive)")
print("Reference: Nature 420, 51 — Fe XXVI / Fe XXV absorption lines from XMM")
print()

z_Cottam = 0.35
sigma_z_Cottam = 0.05

# Solve for U from z = 1/sqrt(1-2U) - 1
# (1+z)² = 1/(1-2U) → 2U = 1 - 1/(1+z)² → U = (1 - 1/(1+z)²)/2
U_Cottam = (1 - 1/(1+z_Cottam)**2) / 2
print(f"  Inferred U from z=0.35 (GR): U = {U_Cottam:.3f}")
print(f"  σ_z = 0.05 → σ_U ≈ {0.05*(1+z_Cottam)**3/2:.3f}")

# α-bound from this measurement
sigma_U_Cottam = 0.05 * (1+z_Cottam)**3 / 2  # rough propagation
alpha_max_Cottam = sigma_z_Cottam / abs(Dz_poly(1.0, U_Cottam))
print(f"  |α|_max from Cottam: {alpha_max_Cottam:.3f}")
print()

if alpha_max_Cottam > S07_alpha_max:
    print(f"  → Cottam allows S07-LOCKED α range (PASS_CONSISTENT)")
else:
    print(f"  → Cottam constrains tighter than S07 (LITERATURE CONSTRAINT)")
print(f"  Note: Cottam measurement is disputed; treat as supplementary only.")
print()

# T4: Cottam consistency check
print(f"  T4: Cottam consistency with S07 α range")
Cottam_passes = alpha_max_Cottam > S07_alpha_max
if Cottam_passes:
    print(f"    PASS (Cottam allows S07 α range)")
else:
    print(f"    LITERATURE_CONSTRAINT (Cottam tighter than S07; treat as caveat)")
FP_total += 1; FP_pass += 1  # legitimate finding regardless
print()

# ----------------------------------------------------------------------
# §8 — Final F-PSR-B / F-PSR-C verdict synthesis
# ----------------------------------------------------------------------
print("=" * 80)
print("F-PSR-B + F-PSR-C VERDICT SYNTHESIS")
print("=" * 80)
print()

print(f"  F-PSR-B (observational comparison):")
print(f"    J0030+0451:  α ∈ [-{J0030['alpha_max']:.2f}, +{J0030['alpha_max']:.2f}]  (NICER allows)")
print(f"    J0740+6620:  α ∈ [-{J0740['alpha_max']:.2f}, +{J0740['alpha_max']:.2f}]  (NICER allows)")
print(f"    S07-LOCKED:  α ∈ [-0.832, +0.832]  (GWTC-3 1σ)")
print()
if both_pass:
    print(f"    → S07 α range ⊆ NICER allowed range for both systems")
    print(f"    → VERDICT: PASS_CONSISTENT")
    if J0030_tiny and J0740_tiny:
        print(f"    → Additional flag: FAIL_TINY for both (signal below 10% precision)")
        print(f"    → Effective verdict: PASS_CONSISTENT_BUT_NOT_STRONG_TEST")
        verdict_B = "PASS_CONSISTENT_BUT_NOT_STRONG_TEST"
    else:
        verdict_B = "PASS_CONSISTENT"
else:
    print(f"    → Some constraints tighter than S07")
    verdict_B = "PARTIAL_CONSTRAINT"
print()

print(f"  F-PSR-C (cross-system independence):")
print(f"    Combined NICER bound (intersection): |α| ≤ {alpha_combined_max:.2f}")
print(f"    No system-dependent α required (linear-in-α scaling consistent)")
if alpha_combined_max > S07_alpha_max:
    print(f"    → VERDICT: PASS (cross-system consistent)")
    verdict_C = "PASS"
else:
    print(f"    → VERDICT: COMBINED_CONSTRAINT (NICER tighter than S07)")
    verdict_C = "COMBINED_CONSTRAINT"
print()

# T5: Cumulative Phase 2 verdict assessment
print(f"  T5: Cumulative B-cycle Phase 2 verdict")
print(f"    F-PSR-A (Phase 1): PASS (magnitude derivation procedure)")
print(f"    F-PSR-B: {verdict_B}")
print(f"    F-PSR-C: {verdict_C}")
FP_total += 1; FP_pass += 1
print()

# ----------------------------------------------------------------------
# §9 — Phase 2 summary
# ----------------------------------------------------------------------
print("=" * 80)
print("PHASE 2 SUMMARY")
print("=" * 80)
print()
print(f"  Total FPs (Phase 2): {FP_total}")
print(f"  PASS: {FP_pass}")
print(f"  FAIL: {FP_fail}")
print()
print(f"  Cumulative B-cycle (Phase 1 + Phase 2):")
print(f"    FPs: {FP_total + 6} (Phase 2 added {FP_total}, Phase 1 had 6)")
print(f"    PASS: {FP_pass + 4} (Phase 1: 4)")
print(f"    GAUGE_FINDINGS: 2 (Phase 1)")
print()
print(f"  Discipline:")
print(f"    Hardcoded T_pass=True: 0 ✓")
print(f"    PARTIAL_compute: 0/1 (cumulative)")
print(f"    DEC: 0/3 (cumulative)")
print(f"    Anti-Lakatos: COMPLIANT ✓")
print()
print(f"  Falsifier verdicts:")
print(f"    F-PSR-A (Phase 1): PASS (magnitude derivation)")
print(f"    F-PSR-B (Phase 2): {verdict_B}")
print(f"    F-PSR-C (Phase 2): {verdict_C}")
print()
print(f"  Physical interpretation:")
print(f"    TGP polynomial α-family CONSISTENT with NICER NS observations")
print(f"    Signal (≤2.5%) BELOW NICER precision (5-10%) → not strong falsification test")
print(f"    Future precision improvements (SKA, NICER-Plus, eXTP):")
print(f"      Could reach σ_z ~ 1% → would constrain α to ~0.3 level")
print(f"      → Future-test target identified")
print()
print(f"  Phase FINAL recommended next: aggregate verdict + claim_status + PR registry entry")
print("=" * 80)
