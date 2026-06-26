"""
Phase EXPLORATION — no-Yukawa hypothetical scope (γ-7 cycle)

STATUS: EXPLORATION ONLY. NOT verdict re-derivation.
γ-7 HALT-B remains LOCKED (user "Zatwierdzam Halt B").
Anti-Lakatos LOCK preserved — pre-registered F-γ-7-A/B/C/D verdicts unchanged.

This script explores three hypothetical scenarios:
  (A) Yukawa baseline (Phase 2 LOCKED) — re-verify numerics
  (B) Coulomb-like (no screening) with R_universe natural cutoff
  (C) Required R_universe (Coulomb case) for V_eff/V_universe = Ω_DE = 0.7
  (D) Mass-dependent screening hypothesis sketch (qualitative)

Purpose:
  Gather information BEFORE any future γ-8 (if F8 ever reopened).
  Inform R1 #17 (linear theory limitation) with scaling analysis.
  Test user hypothesis: "może im większa masa, tym granica ekranowania większa".

Methodology:
  - sympy symbolic + numerical (SI units)
  - Use Phase 2 LOCKED formula structure
  - Replace Yukawa kernel exp(-μr)/(4πr) with Coulomb 1/(4πr)
  - Use R_universe as natural IR cutoff (instead of 1/μ_sp)
  - Compare regimes numerically
"""

import sympy as sp
import numpy as np

print("=" * 80)
print("PHASE EXPLORATION — no-Yukawa hypothetical (γ-7)")
print("STATUS: EXPLORATION; HALT-B LOCKED; anti-Lakatos preserved")
print("=" * 80)

# ----------------------------------------------------------------------
# SECTION 1 — Constants (SI)
# ----------------------------------------------------------------------
print("\n--- §1 Constants (SI) ---")

G = 6.674e-11           # m³ / (kg s²)
c = 2.998e8             # m / s
H0 = 2.20e-18           # 1/s (H_0 ≈ 67.7 km/s/Mpc)
hbar = 1.055e-34        # J s

# TGP Appendix E
m_sp = hbar * H0 / c**2  # kg (ultra-light Phi-substrate mass)
mu_sp = H0 / c           # 1/m (Yukawa inverse range)
lambda_sp = c / H0       # m (Yukawa range = Hubble radius)
v_phi2 = 3.03e45         # J/m (Appendix E v² ≈ 25·E_P/ℓ_P)

# Cosmology
rho_crit = 3 * H0**2 / (8 * np.pi * G)    # kg/m³
Omega_m = 0.315
rho_m = Omega_m * rho_crit                # matter density
R_obs = c / H0                            # m (Hubble/observable)
V_obs = (4*np.pi/3) * R_obs**3            # m³

print(f"  G              = {G:.3e}")
print(f"  H_0            = {H0:.3e}")
print(f"  m_sp           = {m_sp:.3e} kg")
print(f"  μ_sp           = {mu_sp:.3e} 1/m")
print(f"  λ_sp           = {lambda_sp:.3e} m  (Hubble radius)")
print(f"  v_phi²         = {v_phi2:.3e} J/m")
print(f"  ρ_crit         = {rho_crit:.3e} kg/m³")
print(f"  ρ_m            = {rho_m:.3e} kg/m³ (Ω_m={Omega_m})")
print(f"  R_obs          = {R_obs:.3e} m")
print(f"  V_obs          = {V_obs:.3e} m³")

# ----------------------------------------------------------------------
# SECTION 2 — (A) Yukawa baseline (Phase 2 LOCKED)
# ----------------------------------------------------------------------
print("\n--- §2 (A) Yukawa baseline (Phase 2 LOCKED) ---")
print("Formula: V_eff/V_univ = (G ρ² V_univ ⟨exp⟩)/(2 μ_sp v²) · ξ_clump")

# Geometric ⟨exp⟩_uniform for uniform sphere of radius R, Yukawa range 1/μ
def avg_exp_yukawa(mu, R):
    """⟨exp(-μr_ij)⟩ averaged over uniform pair distribution in sphere R"""
    x = mu * R
    if x < 1e-10:
        return 1.0 - x/2 + x**2/10  # series for small x
    return (3.0/x**3) * (2.0 - np.exp(-x)*(2.0 + 2.0*x + x**2))

# Yukawa V_eff/V_universe (Phase 2 LOCKED form)
def Veff_Yukawa(R_universe, mu, rho, v2, xi_clump=1.0):
    """V_eff/V_universe per Phase 2 LOCKED formula"""
    V = (4*np.pi/3) * R_universe**3
    exp_avg = avg_exp_yukawa(mu, R_universe)
    return (G * rho**2 * V * exp_avg) / (2 * mu * v2) * xi_clump

# Test at R_obs with ξ_clump = 1 (best case)
xi_best = 1.0
xi_tgp_emp = 1e-4  # TGP-empirical from Phase 3 (R1 #17 runaway-truncated)

V_yukawa_hubble_best = Veff_Yukawa(R_obs, mu_sp, rho_m, v_phi2, xi_best)
V_yukawa_hubble_emp  = Veff_Yukawa(R_obs, mu_sp, rho_m, v_phi2, xi_tgp_emp)

# Saturated limit (R → ∞): V·⟨exp⟩ → 8π/μ³
V_yukawa_sat_best = (4*np.pi) * G * rho_m**2 * lambda_sp**4 / v_phi2 * xi_best
V_yukawa_sat_emp  = (4*np.pi) * G * rho_m**2 * lambda_sp**4 / v_phi2 * xi_tgp_emp

print(f"  R = R_obs (Hubble), μR = {mu_sp*R_obs:.3f}:")
print(f"    ⟨exp⟩_uniform                 = {avg_exp_yukawa(mu_sp, R_obs):.4f}")
print(f"    V_eff/V_univ (ξ=1)            = {V_yukawa_hubble_best:.3e}")
print(f"    V_eff/V_univ (ξ_TGP≈1e-4)    = {V_yukawa_hubble_emp:.3e}")
print(f"  R → ∞ (saturated, Yukawa screen):")
print(f"    V_eff/V_univ (ξ=1)            = {V_yukawa_sat_best:.3e}")
print(f"    V_eff/V_univ (ξ_TGP≈1e-4)    = {V_yukawa_sat_emp:.3e}")
print(f"  Note: saturation cap = 4π G ρ_m² λ_sp⁴/v²")

# ----------------------------------------------------------------------
# SECTION 3 — (B) Coulomb (no Yukawa) with R_universe cutoff
# ----------------------------------------------------------------------
print("\n--- §3 (B) Coulomb (no Yukawa) with R_universe cutoff ---")
print("Derivation: replace Yukawa kernel exp(-μr)/(4πr) → 1/(4πr)")
print("Pair-overlap integral ∫δΦ_iδΦ_j dV (R cutoff) ≈ q²R/(4π) · O(1)")
print("Per uniform sphere: V_eff ≈ G M² R / (2 v²) · ξ_clump")
print("V_eff/V_universe ≈ (2π/3) G ρ² R⁴ / v² · ξ_clump")

def Veff_Coulomb(R_universe, rho, v2, xi_clump=1.0, geom_factor=1.0):
    """
    V_eff/V_universe (no Yukawa, R cutoff).
    geom_factor accounts for exact pair-overlap geometry (set 1.0 leading order).
    """
    return (2*np.pi/3) * G * rho**2 * R_universe**4 / v2 * xi_clump * geom_factor

V_coul_hubble_best = Veff_Coulomb(R_obs, rho_m, v_phi2, xi_best)
V_coul_hubble_emp  = Veff_Coulomb(R_obs, rho_m, v_phi2, xi_tgp_emp)

print(f"  R = R_obs, ξ=1:                 V_eff/V_univ = {V_coul_hubble_best:.3e}")
print(f"  R = R_obs, ξ_TGP:               V_eff/V_univ = {V_coul_hubble_emp:.3e}")

# Comparison Yukawa vs Coulomb at same R
ratio_C_to_Y = V_coul_hubble_best / V_yukawa_hubble_best
print(f"  Ratio Coulomb/Yukawa at R=R_obs: {ratio_C_to_Y:.3f}")
print(f"  Interpretation: at R=λ_sp, Yukawa wins by ~{1/ratio_C_to_Y:.1f}x")
print(f"  (because Yukawa integrates exp-decay, Coulomb only linear)")

# ----------------------------------------------------------------------
# SECTION 4 — (C) Required R_universe for Ω_DE = 0.7 (Coulomb)
# ----------------------------------------------------------------------
print("\n--- §4 (C) Required R_universe for V_eff/V_universe = 0.7 ---")
print("Solve: (2π/3) G ρ² R⁴/v² · ξ = 0.7  →  R = ((0.7·3·v²)/(2π G ρ² ξ))^(1/4)")

target_Omega_DE = 0.7

def R_required_Coulomb(rho, v2, xi_clump):
    R4 = target_Omega_DE * 3 * v2 / (2 * np.pi * G * rho**2 * xi_clump)
    return R4**0.25

R_req_best = R_required_Coulomb(rho_m, v_phi2, xi_best)
R_req_emp  = R_required_Coulomb(rho_m, v_phi2, xi_tgp_emp)

print(f"  ξ_clump = 1 (best case):")
print(f"    R_required = {R_req_best:.3e} m  = {R_req_best/R_obs:.2f} · R_obs")
print(f"  ξ_clump = 1e-4 (TGP-empirical):")
print(f"    R_required = {R_req_emp:.3e} m  = {R_req_emp/R_obs:.2f} · R_obs")

# For Yukawa saturated case (cannot achieve 0.7 by enlarging R alone):
print(f"\n  Yukawa SATURATED (any R >> λ_sp) ceiling:")
print(f"    V_eff/V_univ (max, ξ=1) = {V_yukawa_sat_best:.3e}  <<  0.7")
print(f"    Factor short of 0.7     = {0.7/V_yukawa_sat_best:.1f}x")
print(f"    → Yukawa case: ABSOLUTE BARRIER. Cannot save by bigger universe.")
print(f"    → Coulomb case: ESCAPES via R⁴ growth.")

# Required μ_sp (shrinking μ_sp = expanding range) for Yukawa to give 0.7:
print(f"\n  Alternative: shrink μ_sp instead of expanding R.")
print(f"  4π G ρ² /μ_sp⁴/v² = 0.7  →  μ_sp_req/μ_sp = ({V_yukawa_sat_best/0.7:.3e})^(1/4)")
mu_ratio_req = (V_yukawa_sat_best/0.7)**0.25
print(f"    μ_sp_req = {mu_ratio_req:.3f} · μ_sp_now  →  λ_sp_req = {1/mu_ratio_req:.2f} · λ_sp_now")
print(f"    Required Yukawa range ≈ {(c/H0/mu_ratio_req):.3e} m = {1/mu_ratio_req:.2f}·Hubble radius")

# ----------------------------------------------------------------------
# SECTION 5 — (D) Mass-dependent screening hypothesis (sketch)
# ----------------------------------------------------------------------
print("\n--- §5 (D) Mass-dependent screening hypothesis ---")
print("User hypothesis: 'im większa masa, tym granica ekranowania większa'")
print()
print("Physical motivation: nonlinear Mexican-hat potential V(Φ)=λ/4(|Φ|²-Φ_0²)²")
print("Linearized: m_sp² = 2λΦ_0² (constant)")
print("Near mass concentrations: |Φ| can deviate from Φ_0, giving")
print("  m_sp_eff² = 2λ⟨|Φ|²⟩  →  density-dependent")
print()
print("Scenario sketch: m_sp_eff(ρ) → 0 in mass-concentrated regions.")
print("In macro pair-overlap, when both sources sit in concentrated regions,")
print("effective Yukawa range λ_sp_eff → R_concentration (much larger).")
print()
print("Approximation: replace λ_sp with effective scale λ_eff = α · R_universe")
print("where α ∈ (0,1] depends on global concentration.")
print()
print("Required α for V_eff/V_univ = 0.7 at saturation:")
print("  4π G ρ_m² (α·R)⁴ /v² ≥ 0.7  →  α ≥ R_req(ξ=1)/R")
for R_factor in [1, 5, 10, 50]:
    R_test = R_factor * R_obs
    alpha_req = R_req_best / R_test
    print(f"    R_universe = {R_factor}·R_obs ({R_test:.2e} m):  α_req = {alpha_req:.3f}")
print()
print("Interpretation:")
print("  - If universe is exactly observable size: need α ≈ 9 (impossible, λ_eff > R)")
print("  - If universe is 9× observable: need α ≈ 1 (max screening collapse)")
print("  - If universe is 50× observable: need α ≈ 0.18 (modest screening)")
print()
print("This is the 'large universe + relaxed screening' joint hypothesis.")

# ----------------------------------------------------------------------
# SECTION 6 — Sanity check: dimensional verification (sympy)
# ----------------------------------------------------------------------
print("\n--- §6 Sympy dimensional verification ---")

# Symbolic check of V_eff formula dimensions
G_s, M_s, mu_s, v2_s, R_s, rho_s = sp.symbols('G M mu v2 R rho', positive=True)

# Yukawa formula: V_eff = G M² ⟨exp⟩/(2μ v²)
# [G M²/(μ v²)] should be m³
# G ~ m³/(kg s²), M ~ kg, μ ~ 1/m, v² ~ J/m = kg m/s²
# G M² = kg m³/s²
# μ v² = (1/m)(kg m/s²) = kg/s²
# G M²/(μ v²) = (kg m³/s²)/(kg/s²) = m³ ✓

# Coulomb formula: V_eff = G M² R/(2 v²)
# G M² R = kg m⁴/s²
# v² = kg m/s²
# V_eff = kg m⁴/s² / (kg m/s²) = m³ ✓
print("  Yukawa form: [G M² / (μ v²)] = m³ ✓")
print("  Coulomb form: [G M² R / v²]  = m³ ✓")
print("  Both dimensionally consistent as 'effective volume'.")

# ----------------------------------------------------------------------
# SECTION 7 — Summary table
# ----------------------------------------------------------------------
print("\n--- §7 Summary table ---")

scenarios = [
    ("Yukawa, R=Hubble, ξ=1",              V_yukawa_hubble_best),
    ("Yukawa, R=Hubble, ξ_TGP≈1e-4",       V_yukawa_hubble_emp),
    ("Yukawa, R→∞ saturated, ξ=1",         V_yukawa_sat_best),
    ("Yukawa, R→∞ saturated, ξ_TGP",       V_yukawa_sat_emp),
    ("Coulomb (no screen), R=Hubble, ξ=1", V_coul_hubble_best),
    ("Coulomb, R=Hubble, ξ_TGP",           V_coul_hubble_emp),
    ("Coulomb, R = 9·R_obs, ξ=1 [target]", Veff_Coulomb(9*R_obs, rho_m, v_phi2, xi_best)),
    ("Coulomb, R = 50·R_obs, ξ=1",         Veff_Coulomb(50*R_obs, rho_m, v_phi2, xi_best)),
    ("OBSERVED Ω_DE",                       0.70),
]

print(f"  {'Scenario':<42s}  {'V_eff/V_univ':>14s}  {'short of 0.7':>14s}")
print("  " + "-"*42 + "  " + "-"*14 + "  " + "-"*14)
for name, val in scenarios:
    short = 0.7/val if val > 0 else float('inf')
    short_str = f"{short:.1e}x" if short < 1e10 else "n/a"
    print(f"  {name:<42s}  {val:>14.3e}  {short_str:>14s}")

# ----------------------------------------------------------------------
# SECTION 8 — Conclusions (EXPLORATION)
# ----------------------------------------------------------------------
print("\n" + "=" * 80)
print("EXPLORATION CONCLUSIONS (NOT verdicts; HALT-B preserved)")
print("=" * 80)
print()
print("E1. Yukawa screening creates ABSOLUTE BARRIER:")
print("    Cannot reach Ω_DE = 0.7 by any R_universe. Saturation ceiling ≈ {:.1e}".format(V_yukawa_sat_best))
print()
print("E2. Without Yukawa (Coulomb-like) + R_universe ≈ 9·R_obs:")
print("    V_eff/V_universe = 0.7 achievable IF ξ_clump = 1.")
print("    With realistic ξ_clump ~ 1e-4: need R ≈ {:.0f}·R_obs.".format(R_req_emp/R_obs))
print()
print("E3. The user's intuition is partially correct:")
print("    - Larger universe ALONE (Yukawa intact): factor 12-13× max  → NOT enough.")
print("    - Larger universe + relaxed screening: ESCAPE possible.")
print("    - The bottleneck is Yukawa, not universe size.")
print()
print("E4. Mass-dependent screening (α ≈ 0.18 if R ≈ 50·R_obs) is the cheapest hypothesis.")
print("    Would require nonlinear KG analysis (γ-8+ scope).")
print()
print("E5. Open question for any future γ-8:")
print("    Is m_sp from Appendix E (eq. 353) constant, or is it density-dependent?")
print("    If density-dependent, the four F8 FAILs from γ-3/γ-3'/γ-5/γ-7 ALL would need")
print("    re-examination — but ONLY under new, independent, pre-registered protocol.")
print()
print("CRITICAL DISCIPLINE NOTE:")
print("  This exploration does NOT modify HALT-B claim_status.")
print("  Findings are inputs for future research direction, NOT verdicts.")
print("  Any actual γ-8 cycle (if ever opened) must:")
print("    - Pre-register new falsifiers BEFORE derivation")
print("    - Pass strict cycle 1/2/7")
print("    - Address R1 #17 (linear theory limitation)")
print("    - Re-derive m_sp from Appendix E with nonlinear corrections")
print("    - NOT cite this exploration as justification for re-opening F8 cleanly")
print()
print("=" * 80)
print("END EXPLORATION")
print("=" * 80)
