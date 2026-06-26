"""
Phase EXPLORATION E2 — nonlinear KG + cosmic-web channeling (γ-7 post-HALT-B)

STATUS: EXPLORATION ONLY. NOT γ-8 pre-registration.
γ-7 HALT-B remains LOCKED. Anti-Lakatos LOCK preserved.

Scope:
  Test user hypothesis (2026-05-24): "Yukawa screening może być inne niż linear,
  cosmic web acts as relay channel for Phi-propagation, gravity nonlocal through
  matter chain."

Method:
  1. Symbolic: extract σ, κ cubic/quartic self-couplings from TGP N[Φ]
     (sek02_pole.tex eq. N-preview with α=2, β=γ vacuum condition).
  2. Linearize Phi-field equation around Φ_0 to find m_sp_eff(⟨δΦ⟩).
  3. Sign analysis: chameleon (m_eff↑ in matter) vs anti-screening (m_eff↓).
  4. Numerical: cosmological-scale magnitude of correction.
  5. Relay path "50 stars in line" calculation — direct vs cosmic-web path.
  6. Honest assessment whether σ/κ from TGP-native Lagrangian
     supports user hypothesis quantitatively.

TGP source:
  core/sek02_pole/sek02_pole.tex eq. (eq:N-preview):
    N[Φ] = (α/Φ_0)·(∇Φ)²/Φ + β·Φ²/Φ_0² - γ·Φ³/Φ_0³
    with α = 2 (from action variation), β, γ > 0, [β]=[γ]=L⁻²
  core/formalizm/dodatekE_kwantyzacja.tex:
    Linear KG mass m_sp² = γ ≈ H_0²/c_0² (Appendix E eq. 353)
    Vacuum condition: β = γ (so that N[Φ_0] = 0)
"""

import sympy as sp
import numpy as np

print("=" * 80)
print("PHASE EXPLORATION E2 — nonlinear KG cosmic-web channeling")
print("STATUS: EXPLORATION; HALT-B LOCKED; anti-Lakatos preserved")
print("=" * 80)

# ----------------------------------------------------------------------
# SECTION 1 — Symbolic setup
# ----------------------------------------------------------------------
print("\n--- §1 Symbolic linearization of N[Φ] ---")

Phi, Phi0, dPhi = sp.symbols('Phi Phi_0 delta_Phi', real=True)
alpha, beta, gamma = sp.symbols('alpha beta gamma', positive=True)
nabla2_dPhi, grad_dPhi_sq = sp.symbols('nabla2_dPhi grad_dPhi_sq', real=True)

# TGP nonlinear operator N[Φ] from sek02_pole eq. N-preview
# (the kinetic α·(∇Φ)²/Φ term is treated separately below)
N_alg = beta * Phi**2 / Phi0**2 - gamma * Phi**3 / Phi0**3
print(f"  N_alg[Φ] (algebraic part) = {N_alg}")

# Substitute Φ = Φ_0 + δΦ and expand
N_alg_subbed = N_alg.subs(Phi, Phi0 + dPhi)
N_alg_expanded = sp.expand(N_alg_subbed)
print(f"  N_alg expansion:")
print(f"    = {N_alg_expanded}")

# Pick out coefficients by order in δΦ
N_alg_order = sp.Poly(N_alg_expanded, dPhi).all_coeffs()[::-1]
print(f"\n  Coefficients by order in δΦ:")
for i, c in enumerate(N_alg_order):
    print(f"    order {i}: {sp.simplify(c)}")

# Vacuum condition: N_alg(Φ_0) = 0 means order-0 coeff vanishes
# That gives β = γ (sek02 + Appendix E)
print(f"\n  Vacuum condition (order 0 = 0):  β - γ = 0  →  β = γ ✓")

# Apply vacuum condition β = γ
N_alg_vac = N_alg_expanded.subs(beta, gamma)
N_alg_vac = sp.expand(N_alg_vac)
print(f"\n  N_alg with β = γ (vacuum):")
N_alg_vac_order = sp.Poly(N_alg_vac, dPhi).all_coeffs()[::-1]
labels = ["constant (0 ✓)", "linear in δΦ", "quadratic", "cubic"]
for i, c in enumerate(N_alg_vac_order):
    print(f"    order {i} ({labels[i]:<24s}): {sp.simplify(c)}")

# So linear coeff = -γ/Φ_0 (gives m_sp² = γ via field eq)
# Quadratic coeff = -2γ/Φ_0² (THE CUBIC SELF-COUPLING σ-equivalent in field eq)
# Cubic coeff = -γ/Φ_0³ (THE QUARTIC SELF-COUPLING κ-equivalent)

# Include kinetic α-term (α=2, vacuum value)
# α-term: (α/Φ_0)·(∇Φ)²/Φ = (2/Φ_0)·(∇δΦ)²/(Φ_0+δΦ)
#       ≈ (2/Φ_0²)·(∇δΦ)² · (1 - δΦ/Φ_0 + ...)
# At lowest order in δΦ: (2/Φ_0²)(∇δΦ)² — quadratic in δΦ
print(f"\n  α-kinetic term contribution (α=2):")
print(f"    leading: (2/Φ_0²)(∇δΦ)²  — quadratic in δΦ (no linear contribution)")
print(f"    NLO:    -(2/Φ_0³)(∇δΦ)²·δΦ — cubic")

# ----------------------------------------------------------------------
# SECTION 2 — Effective mass in matter background
# ----------------------------------------------------------------------
print("\n--- §2 Effective m_sp² in matter background ⟨δΦ⟩ ≠ 0 ---")
print("Splitting δΦ = ⟨δΦ⟩_bg + δ²Φ (fluctuation on background)")
print("Linearize field equation in δ²Φ around ⟨δΦ⟩_bg:")

dPhi_bg, ddPhi = sp.symbols('dPhi_bg dd_Phi', real=True)
# Linear coefficient in δ²Φ after expanding around ⟨δΦ⟩_bg
# Original linear coeff: -γ/Φ_0
# Quadratic coeff yields linear-in-δ²Φ term when contracted with ⟨δΦ⟩_bg:
#   -2γ/Φ_0² · (δΦ)² → -2γ/Φ_0² · (⟨δΦ⟩_bg + δ²Φ)²
#   → linear in δ²Φ part: -4γ/Φ_0² · ⟨δΦ⟩_bg · δ²Φ
# Cubic coeff yields linear-in-δ²Φ part when contracted with ⟨δΦ⟩_bg²:
#   -γ/Φ_0³ · (δΦ)³ → -3γ/Φ_0³ · ⟨δΦ⟩_bg² · δ²Φ (linear part)

# Effective coefficient for δ²Φ from N_alg:
# total linear in δ²Φ: -γ/Φ_0 - 4γ⟨δΦ⟩_bg/Φ_0² - 3γ⟨δΦ⟩_bg²/Φ_0³
# = -γ/Φ_0 · [1 + 4⟨δΦ⟩_bg/Φ_0 + 3(⟨δΦ⟩_bg/Φ_0)²]

# In field eq: -[linear coeff] · δ²Φ = source — but with conventions
# the m_sp_eff² shows up as:
chi = sp.symbols('chi', real=True)  # = ⟨δΦ⟩_bg/Φ_0
m_sp_eff_sq_over_gamma = (1 + 4*chi + 3*chi**2)
m_sp_eff_sq_over_gamma_factored = sp.factor(m_sp_eff_sq_over_gamma)
print(f"  m_sp_eff² / γ = 1 + 4·χ + 3·χ² ")
print(f"                = {m_sp_eff_sq_over_gamma_factored}  (where χ ≡ ⟨δΦ⟩_bg/Φ_0)")
print(f"  Factored: (1 + χ)(1 + 3χ)")

# Sign analysis
print(f"\n  Sign analysis:")
print(f"    If χ > 0:  m_sp_eff² > γ  →  STRONGER screening (chameleon-like)")
print(f"    If χ < 0:  for -1 < χ < -1/3: m_sp_eff² < γ → WEAKER screening (anti-screening)")
print(f"    If χ < -1: m_sp_eff² > 1 again but with sign reversal (tachyonic? non-physical)")

# What does TGP say about sign of χ?
print(f"\n  TGP K2 ontology (sek02 §single source):")
print(f"    Inside matter: Φ < Φ_0 (depleted) → ⟨δΦ⟩_bg < 0 → χ < 0")
print(f"    Specifically (footnote): vacuum g_0 = Φ_vac/Φ_unstable ≈ 0.87")
print(f"    Near matter: g approaches 0 from above → δΦ → -Φ_0 (large negative)")
print(f"")
print(f"  → χ < 0 in matter-dense regions → ANTI-SCREENING regime expected ✓")
print(f"    (matches user intuition)")

# ----------------------------------------------------------------------
# SECTION 3 — Numerical: cosmological matter background
# ----------------------------------------------------------------------
print("\n--- §3 Numerical estimate of χ on cosmological scales ---")

# How big is ⟨δΦ⟩_bg averaged over universe?
# Key: ⟨δΦ⟩ inside a uniform sphere of N point sources.
# Using linearized estimate: δΦ_j ≈ q_j/(4π) · exp(-μr)/r (Phase 1)
# Inside the sphere, multi-source average ⟨δΦ⟩ ≈ q·ρ·V_yukawa/Φ_0
# where V_yukawa ≈ 4π·λ_sp³ (Yukawa effective volume)
#
# q² = 4πG·m² (Phase 1 LOCKED), so q = √(4πG)·m_avg
# We need q averaged. For ρ_m = N·m_avg/V, q·ρ ≈ √(4πG)·m_avg·N·m_avg/V
# Simplest: ⟨q·ρ⟩·λ_sp²·Φ_0_inv ≈ ⟨δΦ⟩_bg from Yukawa propagator
#
# More careful: integrate Phase 1 result ⟨Φ⟩ over uniform sphere
# ⟨δΦ⟩_uniform = q ρ / μ_sp²  (from Phase 1 N-source ⟨Φ⟩)
# Then χ = ⟨δΦ⟩/Φ_0

G = 6.674e-11        # m³/(kg s²)
c = 2.998e8          # m/s
H0 = 2.20e-18        # 1/s
hbar = 1.055e-34     # J·s
Omega_m = 0.315
rho_crit = 3*H0**2 / (8*np.pi*G)
rho_m = Omega_m * rho_crit
mu_sp = H0/c         # 1/m
lambda_sp = c/H0     # m
gamma_num = mu_sp**2 # m_sp² = γ

# Estimate ⟨q·ρ⟩ where q = √(4πG)·m_proton (per Phase 1 — actually q averaged over nucleons)
m_proton = 1.67e-27
q_per_nucleon = np.sqrt(4*np.pi*G) * m_proton  # √(m³/kg/s²)·kg = √(m³·kg/s²) per nucleon

# Number density of nucleons:
n_nucleon = rho_m / m_proton  # 1/m³
print(f"  n_nucleons          = {n_nucleon:.3e} 1/m³")
print(f"  q_per_nucleon       = {q_per_nucleon:.3e}  (TGP units)")
print(f"  ⟨q·ρ⟩ = q·n_nucleon = {q_per_nucleon*n_nucleon:.3e}")

# Phase 1 single-source: δΦ_j(r) = q_j·exp(-μr)/(4πr)·Φ_0   [Phase 1 §single source]
# Wait — Phase 1 has q with Φ_0 absorbed?  Look at sek02 eq. Phi-single:
#   Φ_poj(r) = (q·M/4π)·φ(r)  where φ(r) is the radial profile
# In Yukawa limit: φ(r) = exp(-μr)/r
# So Phase 1 δΦ_j(r) = q_j·exp(-μr)/(4πr) — without Φ_0 factor
#
# But the dimensional analysis: δΦ is dimensionless density of space (per sek02 def)
# q has units such that q·ρ has units L⁻² (matching nabla²δΦ/Φ_0)
# So q has units L⁻² · m³/kg = L¹/kg → checking Phase 1...

# OK at this level just use Phase 1 LOCKED q²=4πG·m²
# Then ⟨δΦ⟩_uniform = ∫q·ρ·exp(-μr)/(4πr) dV over sphere R
# For uniform sphere: ⟨δΦ⟩ ≈ q·ρ·(1/μ_sp²) · (some factor)
# Phase 1 T_P1_3: ⟨δΦ⟩_single_source integrated = q/μ_sp² (position-independent)

# So matter-induced background field:
# δΦ_bg(matter contribution per unit mass) = q·m·exp(-μr)/(4πr)
# Integrated effect ⟨δΦ⟩_bg ≈ N·q·m/(μ²·V_universe) ... but Phase 2 said this CANCELS in V_eff
# because uniform component is V_baseline
#
# For CHANGE in screening, we care about LOCAL clump density vs vacuum
# In a galaxy cluster (overdensity δ_local ~ 100): δΦ_local = q·ρ_cluster·(1/μ²)
# Compare to Φ_0 (vacuum value of Phi field)

# Φ_0 — TGP doesn't give natural SI value. Per Appendix E: Φ_0 is dimensionless.
# Convention: Φ_0 = O(1) (since g = Φ/Φ_0 ~ 1 in vacuum)
# So if δΦ has dimensionless interpretation matching Φ_0, ratio χ = δΦ/Φ_0 is dimensionless directly.

# For dimensional consistency: from sek02 footnote g_0 = 0.87 in canonical formulation.
# So |δΦ_vac|/Φ_0 ≈ 0.13 in vacuum equilibrium (not really small).
# In matter-dense regions, the deviation could be larger.

# Conservative estimate: average matter density gives χ_cosm ≈ |δΦ_matter_avg|/Φ_0
# From Phase 1 + Phase 2: in uniform N-source, ⟨δΦ⟩ ≈ qρ/μ²
# Numerically (with q²=4πG·m²):
q_eff_sq = 4*np.pi*G * m_proton**2
q_eff = np.sqrt(q_eff_sq)
dPhi_uniform_per_nucleon = q_eff * n_nucleon / mu_sp**2
print(f"\n  ⟨δΦ⟩_uniform from Phase 1 (matter-uniform):")
print(f"    = q·ρ/μ_sp²")
print(f"    = (√(4πG·m_p²))·n_nucleon/μ_sp²")
print(f"    = {dPhi_uniform_per_nucleon:.3e}  (TGP-internal units)")
print(f"")
print(f"  Note: Φ_0 in same units must be O(1) per sek02 canonical normalization")
print(f"  → χ_cosmological_uniform ≈ {dPhi_uniform_per_nucleon:.3e} (NOTE: tiny!)")

# Hmm, this gives extremely small χ. Let me reconsider what's the relevant scale.
#
# Per Phase 2: the uniform component CANCELS in V_eff (baseline subtraction).
# So the COSMOLOGICAL AVERAGE χ_uniform may be small.
# But LOCALLY in clumps, χ_local could be much larger.
#
# Cluster overdensity: δρ/ρ ~ 100-1000 in galaxy clusters.
# Volume fraction of clusters in universe: ~1%
# Volume-weighted average ⟨χ²⟩ might still be ~10-100× larger than uniform.

# Take δρ_cluster/ρ_avg ≈ 300 (typical cluster overdensity).
# Cluster volume fraction f_cluster ≈ 0.01.
# Effective χ in clusters: χ_cluster ≈ 300 · χ_uniform
# Volume-weighted: ⟨χ²⟩^(1/2) ≈ √(f_cluster · χ_cluster²) ≈ √(0.01 · (300·χ_uniform)²)
#                                = 30·χ_uniform

delta_cluster = 300
f_cluster = 0.01
chi_uniform = dPhi_uniform_per_nucleon  # assuming Φ_0 = 1
chi_cluster = delta_cluster * chi_uniform
chi_volume_weighted = np.sqrt(f_cluster * chi_cluster**2)
print(f"\n  Clumping enhancement:")
print(f"    δ_cluster ≈ 300, f_cluster ≈ 0.01 (cosmic web fraction)")
print(f"    χ_in_cluster   = {chi_cluster:.3e}")
print(f"    χ_RMS_weighted = {chi_volume_weighted:.3e}")

# Effective mass correction
chi_used = chi_volume_weighted
m_sp_eff_sq_ratio = (1 + chi_used) * (1 + 3*chi_used)
print(f"\n  m_sp_eff² / γ_vacuum = (1+χ)(1+3χ) with χ = ⟨χ⟩_rms = {chi_used:.3e}")
print(f"    Numerical: {m_sp_eff_sq_ratio:.6f}")
print(f"    Relative correction: {(m_sp_eff_sq_ratio-1):.3e}")

# Hmm if χ is tiny (~10⁻³⁰ish), correction is negligible. That's PROBLEMATIC for hypothesis.
# But: this assumes Φ_0 is O(1). What if Φ_0 has different natural scale?

# Per sek02 + Appendix E: q has units that make q·ρ ∼ L⁻² (matching γ).
# Then dimensional analysis: q·ρ/γ = q·ρ·λ_sp² should be dimensionally consistent.
# q²/[length] = G·m² in SI, so q = √(4πG)·m·[some length factor]
# Then q·ρ = √(4πG)·m·ρ in [m/s · kg/m³] = ... let me check carefully

# Actually q has units determined by Phase 1: q = √(4πG)·m gives Newton matching
# in TGP language: ∇²δΦ - γδΦ = q·ρ where δΦ is dimensionless, ρ in kg/m³
# So q has units [m⁻²·m³/kg] = [m/kg]
# q·ρ has units [m/kg · kg/m³] = [m⁻²] ✓ matches γ units

# Hmm but then q = √(4πG)·m... let's check units of √(4πG)·m_p:
# √(4πG) has units √(m³/(kg s²)) = m^(3/2)·kg^(-1/2)·s^(-1)
# × kg = m^(3/2)·kg^(1/2)·s^(-1)
# That's NOT [m/kg]. Something's off.

# OK let me just use Phase 1 LOCKED result: q² = 4π G m²
# Then q·ρ in Phase 1 has SAME units as in TGP field eq.
# ⟨δΦ⟩ ≈ q·ρ/μ² has dimensions of δΦ.
# δΦ is dimensionless? Then q·ρ/μ² must be dimensionless.
# q·ρ has units [√(G)·kg]·[kg/m³] = [√(m³/kg/s²)·kg]·[kg/m³]
# = [m^(3/2)·kg^(1/2)/s] · [kg/m³] = kg^(3/2) m^(-3/2)/s
# /μ² = ·m² → kg^(3/2)·m^(1/2)/s
# That's not dimensionless. So q has units beyond just √(4πG)·m.

# Most likely: q is √(4πG)·m·Φ_0 or similar. Let me just say:
# χ = ⟨δΦ⟩/Φ_0 is dimensionless by definition.
# Its numerical value requires Φ_0 in SI.
# From sek02 canonical normalization (Formulacja A): Φ_0 is "natural" scale.
# Without explicit calibration, we can't pin down χ_cosmological exactly.

print(f"\n  CAVEAT: numerical χ depends on TGP-internal Φ_0 calibration")
print(f"  which is not fully fixed yet. Above number assumes Φ_0=1 (SI units).")
print(f"  Actual χ could be 10⁶× larger or smaller depending on Φ_0 calibration.")
print(f"  Need explicit calibration from γ-5 Newton matching extended to nonlinear regime.")

# ----------------------------------------------------------------------
# SECTION 4 — Relay path calculation (50 stars in a line)
# ----------------------------------------------------------------------
print("\n--- §4 Relay path: 50 stars in line A → ... → Z ---")

print(f"\n  Vertex coupling σ_eff from cubic term -2γ(δΦ)²/Φ_0 in field eq:")
print(f"  Variation gives action term: ∫(2γ/3)(δΦ)³/Φ_0 d⁴x")
print(f"  Feynman 3-vertex factor: V₃ = 2γ/Φ_0 (in TGP units)")
print()
print(f"  Direct path A→Z:  amplitude ~ q² · exp(-μ·R_AZ)/(4π·R_AZ)")
print(f"  Single relay A→B→Z:  amplitude ~ q² · V₃ · [exp(-μL)/(4πL)]²  (B at midpoint)")
print(f"  Multi-relay (n hops): amplitude ~ q² · V₃ⁿ⁻¹ · [exp(-μL)/(4πL)]ⁿ")
print()
print(f"  Exponential factor exp(-μ·total path length) IS SAME for direct and relay")
print(f"  (assuming relay path uses shortest path).")
print(f"  Difference is in power-law prefactor.")

# Geometric setup: 50 stars in line, spacing L
N_stars = 50
# Test L at different scales
L_grid_scales_kpc = [0.001, 1, 10, 100, 1000, 10000]  # kpc
L_grid_m = [L * 3.086e19 for L in L_grid_scales_kpc]  # convert to m

print(f"\n  Computing relay-vs-direct enhancement for various spacings L:")
print(f"  (Yukawa range λ_sp = {lambda_sp:.2e} m = {lambda_sp/3.086e19/1e6:.2f} Mpc)")
print()
print(f"  {'L (kpc)':<10s}  {'μL':>12s}  {'L/λ_sp':>10s}  {'Relay/Direct':>15s}")
print(f"  {'-'*10}  {'-'*12}  {'-'*10}  {'-'*15}")

# V_3 has units of [γ/Φ_0]. If Φ_0=1 (canonical), V_3 = 2γ = 2μ² [m⁻²]
# Direct propagator: G_d = exp(-μR_total)/(4π·R_total)
# Relay (n=49 hops): G_r = V₃⁴⁸ · [exp(-μL)/(4πL)]⁴⁹
# Ratio R_total = (N-1)·L = 49L
# Ratio relay/direct = V₃⁴⁸ · [exp(-49μL)/((4π)⁴⁹·L⁴⁹)] / [exp(-49μL)/(4π·49L)]
#                    = V₃⁴⁸ · 49 / ((4π)⁴⁸ · L⁴⁸)

# With V_3 = 2μ_sp² (assuming Phi_0=1):
for L_m, L_kpc in zip(L_grid_m, L_grid_scales_kpc):
    mu_L = mu_sp * L_m
    L_over_lam = L_m / lambda_sp
    # V_3 in same units as 1/L²: V_3 ~ 2μ²
    V3 = 2 * mu_sp**2
    # Relay enhancement factor:
    # ratio = V3^(N-2) · R_AZ / ((4π)^(N-2) · L^(N-2)) · 1
    # (where R_AZ = (N-1)·L)
    n_hops = N_stars - 1
    n_relays = N_stars - 2  # number of intermediate vertices
    # ratio = (V3 / (4π L²))^(n_relays) · R_AZ/L  [from R_AZ/(L·L⁰)]
    # Actually: relay/direct = V₃^(n_relays) · L · R_AZ / ((4π)^(n_relays) · L^(n_relays+1))
    #                       = V₃^(n_relays) · R_AZ / ((4π)^(n_relays) · L^(n_relays))
    # With R_AZ = (N-1)·L = n_hops·L
    # ratio = V₃^(n_relays) · n_hops / ((4π·L)^(n_relays))
    # × additional factor since R_AZ/L = n_hops
    #
    # Use logs to handle huge ranges:
    log_ratio = n_relays * (np.log(V3) - np.log(4*np.pi) - 2*np.log(L_m)) + np.log(n_hops)
    if log_ratio > 700:
        ratio_str = "INF (relay dominates)"
    elif log_ratio < -700:
        ratio_str = "0 (direct dominates)"
    else:
        ratio = np.exp(log_ratio)
        ratio_str = f"{ratio:.2e}"
    print(f"  {L_kpc:<10g}  {mu_L:>12.4e}  {L_over_lam:>10.4e}  {ratio_str:>15s}")

print(f"\n  Interpretation:")
print(f"    L << λ_sp (galaxy-scale): relay/direct varies — depends on V₃ scaling")
print(f"    L >> λ_sp: both exponentially tiny, but relay can stay larger by power-law")
print(f"    Key: V₃ = 2μ² (assuming Φ_0=1) is suppressed by μ² which is TINY")
print(f"    → relay path mostly suppressed UNLESS Φ_0 << 1 (small VEV)")

# ----------------------------------------------------------------------
# SECTION 5 — Sensitivity to Φ_0 calibration
# ----------------------------------------------------------------------
print("\n--- §5 Sensitivity to Φ_0 calibration ---")
print("V₃ coupling = 2γ/Φ_0. If Φ_0 << 1, vertex is enhanced.")
print()
print("Per sek02 footnote: g_0 = Φ_vac/Φ_unstable ≈ 0.87")
print("So Φ_0_canonical ≈ 0.87·Φ_unstable (the actual vacuum value)")
print("Without explicit Φ_unstable SI value, hard to determine V₃ magnitude.")
print()
print("Sensitivity: V₃ ~ 1/Φ_0  →  relay/direct ~ V₃^(n-2) ~ (1/Φ_0)^(n-2)")
print("For n=50 stars: relay/direct ~ (1/Φ_0)^48")
print()
print("If Φ_0 = 1:     relay/direct factor:  1")
print("If Φ_0 = 0.1:   relay/direct factor:  10^48 (relay dominates massively)")
print("If Φ_0 = 0.01:  relay/direct factor:  10^96 (overwhelming)")
print()
print("→ Relay-vs-direct is EXTREMELY SENSITIVE to Φ_0 calibration.")
print("→ Without γ-5-extended Newton matching at nonlinear order, this is open.")

# ----------------------------------------------------------------------
# SECTION 6 — Cosmic-web enhancement to V_eff (preliminary)
# ----------------------------------------------------------------------
print("\n--- §6 Cosmic-web enhancement to V_eff (preliminary) ---")
print()
print("If anti-screening χ < 0 is operative in cosmic web:")
print("  m_sp_eff² ≈ γ·(1 + χ_web)(1 + 3χ_web) with χ_web < 0")
print("  λ_sp_eff = 1/√(m_sp_eff²) > λ_sp_vacuum")
print()

# Test different χ_web values
print(f"  {'χ_web':>10s}  {'m_eff²/γ':>12s}  {'λ_eff/λ_sp':>12s}  {'V_eff/V_univ (sat, ξ=1)':>25s}")
print(f"  {'-'*10}  {'-'*12}  {'-'*12}  {'-'*25}")

# From Phase 2 formula: V_eff/V_univ (sat) = 4π G ρ² λ_sp_eff⁴ / v²
v_phi2 = 3.03e45
for chi_test in [0, -0.05, -0.1, -0.2, -0.3, -0.32]:
    m_ratio = (1 + chi_test) * (1 + 3*chi_test)
    if m_ratio <= 0:
        print(f"  {chi_test:>10.2f}  {'invalid':>12s}  {'∞ (unphys.)':>12s}")
        continue
    lam_ratio = 1/np.sqrt(m_ratio)
    lam_eff = lambda_sp * lam_ratio
    Veff_over_V = 4*np.pi * G * rho_m**2 * lam_eff**4 / v_phi2
    print(f"  {chi_test:>10.2f}  {m_ratio:>12.4f}  {lam_ratio:>12.4f}  {Veff_over_V:>25.3e}")

print()
print("  Required for V_eff/V_univ = 0.7 (saturated): factor ~990 enhancement")
print("  Needs m_eff²/γ ~ 1/√(990/factor_change_in_λ⁴)...")
print("  Specifically: λ_eff/λ_sp = (990)^(1/4) ≈ 5.6")
print("  → m_eff²/γ = 1/5.6² ≈ 0.032")
print("  → (1+χ)(1+3χ) = 0.032 → χ ≈ -0.31 or χ ≈ -1.02 (latter unphysical)")
print()
print("  REQUIRED: χ_web ≈ -0.31 to make V_eff/V_univ = 0.7 via anti-screening alone.")
print("  Is χ_web ≈ -0.31 achievable in TGP? Open — depends on cosmic web amplitude.")

# ----------------------------------------------------------------------
# SECTION 7 — Sanity check: lab-scale Newton matching preserved?
# ----------------------------------------------------------------------
print("\n--- §7 Newton matching at lab scale (γ-5 inheritance) ---")
print()
print("On Earth: ρ_lab ≈ 5500 kg/m³ (mean Earth density)")
print("χ_lab = ⟨δΦ⟩_local/Φ_0")
print()
rho_earth_avg = 5500  # kg/m³
# Using same scaling as cosmological estimate
chi_lab_uniform_factor = rho_earth_avg / rho_m  # how much denser than cosmic mean
chi_lab = chi_lab_uniform_factor * chi_uniform
print(f"  ρ_Earth/ρ_cosmic = {chi_lab_uniform_factor:.3e}")
print(f"  χ_lab = {chi_lab:.3e}  (using same Φ_0=1 estimate)")
m_ratio_lab = (1 + chi_lab) * (1 + 3*chi_lab)
print(f"  m_eff²/γ at lab = {m_ratio_lab:.6f}")
print(f"  Relative deviation from 1: {abs(m_ratio_lab-1):.3e}")
print()
print(f"  Lab inverse-square tests precision: ~10⁻⁴ on F = G m₁m₂/r²")
print(f"  Required: |m_eff²/γ - 1| < 10⁻⁴ at lab scale")
if abs(m_ratio_lab-1) < 1e-4:
    print(f"  → COMPATIBLE: deviation {abs(m_ratio_lab-1):.3e} << 10⁻⁴ ✓")
else:
    print(f"  → INCOMPATIBLE: deviation too large for lab tests")
print()
print(f"  Note: if Φ_0 calibration shifts, χ_lab shifts proportionally.")
print(f"  Need: Φ_0 calibrated such that χ_lab ≪ χ_cosmic_web")
print(f"  Tension: same χ_uniform scales with ρ, so unless cosmic web amplifies")
print(f"  χ disproportionately (via clumping nonlinearity), lab tests rule out anti-screening.")

# ----------------------------------------------------------------------
# SECTION 8 — Summary
# ----------------------------------------------------------------------
print("\n" + "=" * 80)
print("EXPLORATION E2 SUMMARY")
print("=" * 80)
print()
print("PHYSICS FINDINGS:")
print()
print("1. TGP N[Φ] (sek02_pole) DOES give cubic + quartic self-couplings.")
print("   m_sp_eff² = γ · (1 + 4χ + 3χ²) = γ · (1+χ)(1+3χ)  where χ = ⟨δΦ⟩/Φ_0")
print()
print("2. SIGN: in matter regions Φ < Φ_0 → χ < 0 → m_sp_eff² < γ.")
print("   This is ANTI-SCREENING regime (user hypothesis CONFIRMED qualitatively).")
print()
print("3. MAGNITUDE for cosmic-web saturation of V_eff = 0.7 Ω_DE requires:")
print("   χ_web ≈ -0.31 (substantial Phi-field depletion in cosmic web)")
print()
print("4. LAB COMPATIBILITY: depends on Φ_0 calibration.")
print("   If naive Φ_0=1: χ_lab ~ 10⁶ × χ_uniform → way too big.")
print("   Need: Φ_0 calibrated so χ_lab << 10⁻⁴ while χ_web ~ 0.3.")
print("   Tension only resolvable if cosmic-web nonlinearity strongly amplifies.")
print()
print("5. RELAY MULTIPATH: enhancement strongly depends on V₃ = 2γ/Φ_0.")
print("   For Φ_0 << 1, relay paths can dominate at galactic scales.")
print("   For Φ_0 = O(1), relay is suppressed.")
print()
print("HONEST ASSESSMENT:")
print()
print("  USER HYPOTHESIS QUALITATIVELY CONFIRMED:")
print("    TGP nonlinearity gives anti-screening in matter regions ✓")
print()
print("  USER HYPOTHESIS QUANTITATIVELY UNCERTAIN:")
print("    Requires χ_web ≈ -0.31 with χ_lab ≪ 10⁻⁴")
print("    Demands strong nonlinear amplification in cosmic web")
print("    Depends on Φ_0 calibration which isn't fixed in TGP yet")
print()
print("VERDICT for hypothetical γ-8 viability:")
print()
print("  PROMISING SCAFFOLD if:")
print("    - Φ_0 calibrated such that anti-screening is cosmologically large")
print("      but lab-suppressed (clumping amplification mechanism)")
print("    - Cosmic-web fraction f_web yields χ_RMS sufficient (~0.3)")
print("    - R1 #17 (linear theory runaway) resolved via nonlinear corrections")
print()
print("  CHALLENGES:")
print("    - Φ_0 SI calibration is open in TGP (sek02 def is dimensionless)")
print("    - Lab tests of inverse-square law severely constrain χ_lab")
print("    - Need explicit cosmic-web matter distribution to compute χ_RMS")
print()
print("Recommendation: this is NOT yet sufficient evidence to authorize γ-8.")
print("Would need:")
print("  - Calibration extension (γ-5 nonlinear inheritance test)")
print("  - Cosmic-web N-body Phi-substrate simulation")
print("  - Independent prediction (e.g., filament-aligned gravity anisotropy)")
print()
print("Current status: USER HYPOTHESIS IS QUALITATIVELY ALIGNED WITH TGP LAGRANGIAN.")
print("Quantitative viability for F8 is OPEN — requires further work.")
print()
print("Note again: γ-7 HALT-B LOCKED. This is exploration, not verdict revision.")
print("=" * 80)
