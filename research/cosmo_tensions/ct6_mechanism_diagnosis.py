#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ct6_mechanism_diagnosis.py: Diagnosis of tachyonic backreaction failure.

ct5 RESULT: B_ψ/H₀² ~ 10⁻⁹ — NINE orders of magnitude too small.

ROOT CAUSE: Scale mismatch.
- Tachyonic instability operates at λ_tach = 2π/μ ~ 5632 Mpc (super-Hubble)
- Structure formation occurs at 1-100 Mpc (sub-Hubble)
- At sub-Hubble scales: ∇²ψ >> γψ → tachyonic term negligible
- At super-Hubble scales: no significant matter perturbations to amplify

THIS SCRIPT:
1. Quantifies the scale mismatch
2. Asks: what γ WOULD be needed for sub-Hubble tachyonic instability?
3. Checks if that γ conflicts with observed Λ
4. Explores ALTERNATIVE backreaction channels within TGP
5. Computes Jensen's inequality correction from spatial averaging
6. Considers nonlinear soliton profile contributions
7. Final assessment: which mechanisms survive?
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from scipy.integrate import solve_ivp
import math

print("="*78)
print("  DIAGNOSIS: WHY TACHYONIC BACKREACTION FAILS + ALTERNATIVES")
print("="*78)

# ==========================================================================
# PARAMETERS
# ==========================================================================
H0 = 67.4      # km/s/Mpc (Planck)
c_over_H0 = 4451  # Mpc (c₀/H₀, Hubble radius)
Om_m = 0.315
Om_L = 0.685

# TGP parameters (from Λ_eff = γ/12 = Ω_Λ·3H₀²/(c₀²))
gamma_nominal = 36 * Om_L  # γc₀²/H₀² ≈ 24.66
m_tach_nominal = math.sqrt(gamma_nominal)  # ≈ 4.97 H₀/c₀
lambda_tach = 2 * math.pi / m_tach_nominal * c_over_H0  # Mpc

print(f"\n  NOMINAL TGP PARAMETERS (from Λ requirement):")
print(f"  γc₀²/H₀² = {gamma_nominal:.2f}")
print(f"  m_tach/H₀ = {m_tach_nominal:.3f}")
print(f"  λ_tach = {lambda_tach:.0f} Mpc")
print(f"  c₀/H₀ = {c_over_H0} Mpc")
print(f"  λ_tach/R_Hubble = {lambda_tach/c_over_H0:.2f}")

# ==========================================================================
# PART 1: SCALE MISMATCH QUANTIFICATION
# ==========================================================================
print(f"\n{'='*78}")
print(f"  PART 1: Scale mismatch analysis")
print(f"{'='*78}")

# At what scales does the tachyonic term become important?
# k² c₀²/a² = γc₀² → k = μ = √γ/c₀ = m_tach · H₀/c₀
# In physical wavelength: λ = 2π/k = 2π c₀/(m_tach·H₀) = 2π·R_H/m_tach

# Structure formation scales:
scales_Mpc = [1, 8, 50, 100, 500, 1000, 3000]
print(f"\n  Scale     k·c₀/H₀    k²c₀²/(γc₀²)   Tachyonic?")
print(f"  {'—'*60}")
for R in scales_Mpc:
    k_c0_H0 = c_over_H0 / R  # k in H₀/c₀ units
    ratio = k_c0_H0**2 / gamma_nominal  # k²/μ² — if >> 1, gradient dominates
    tach = "YES" if ratio < 1 else f"NO (gradient {ratio:.0f}× stronger)"
    print(f"  {R:6d} Mpc  {k_c0_H0:8.1f}      {ratio:12.1f}      {tach}")

# ==========================================================================
# PART 2: WHAT γ WOULD MAKE IT WORK?
# ==========================================================================
print(f"\n{'='*78}")
print(f"  PART 2: Required γ for sub-Hubble tachyonic instability")
print(f"{'='*78}")

# For tachyonic instability at scale R: need γ > (c₀/(R·H₀))²
target_scales = [8, 50, 100, 200, 500]
print(f"\n  λ_tach    γ_required    m_tach/H₀    Λ_predicted/Λ_obs")
print(f"  {'—'*60}")
for R in target_scales:
    gamma_req = (c_over_H0 / R)**2
    m_req = math.sqrt(gamma_req)
    Lambda_ratio = gamma_req / gamma_nominal  # Λ ∝ γ
    print(f"  {R:5d} Mpc   {gamma_req:10.0f}      {m_req:8.1f}       {Lambda_ratio:.0f}×")

print(f"""
  CONCLUSION: For tachyonic instability at S₈ scale (8 Mpc):
    γ required = {(c_over_H0/8)**2:.0f} (vs nominal {gamma_nominal:.1f})
    That's {(c_over_H0/8)**2/gamma_nominal:.0f}× too large → Λ would be {(c_over_H0/8)**2/gamma_nominal:.0f}× observed.

  Even for 500 Mpc (barely sub-Hubble):
    γ required = {(c_over_H0/500)**2:.0f} → Λ would be {(c_over_H0/500)**2/gamma_nominal:.0f}× observed.

  There is NO value of γ that simultaneously gives:
    - Correct Λ (requires γ ≈ {gamma_nominal:.0f})
    - Sub-Hubble tachyonic instability (requires γ >> {gamma_nominal:.0f})

  This rules out perturbative tachyonic backreaction as H₀ tension mechanism.""")

# ==========================================================================
# PART 3: JENSEN'S INEQUALITY CORRECTIONS
# ==========================================================================
print(f"\n{'='*78}")
print(f"  PART 3: Jensen's inequality from nonlinear metric averaging")
print(f"{'='*78}")

# TGP metric: ds² = -(c₀²/ψ)dt² + ψ·δᵢⱼdxⁱdxʲ
# If ψ varies spatially: <1/ψ> ≠ 1/<ψ> and <ψ> ≠ ψ₀
#
# For ψ = 1 + δψ with <δψ> = 0:
# <1/ψ> ≈ 1 + <δψ²> - <δψ³> + ... ≈ 1 + σ²_ψ
# <√ψ> ≈ 1 - σ²_ψ/8
# <ψ^(3/2)> ≈ 1 + (3/8)σ²_ψ
#
# Volume element in TGP: √det(g_spatial) = ψ^(3/2)
# <ψ^(3/2)> > (<ψ>)^(3/2) by Jensen (ψ^(3/2) is convex)

# What is σ²_ψ from weak field?
# δψ = 2Φ/c₀² where Φ is Newtonian potential
# σ²_Φ/c⁴ ≈ variance of Newtonian potential

# Gravitational potential variance from power spectrum
# σ²_Φ ≈ ∫ P_Φ(k) k² dk/(2π²)
# At late times: Φ ~ 10⁻⁵ on all scales (Φ doesn't grow in ΛCDM)
sigma_Phi_c2 = 1e-5  # Φ/c² typical value (dimensionless)
sigma2_psi_wf = (2 * sigma_Phi_c2)**2  # from δψ = 2Φ/c₀²
print(f"\n  Weak field: δψ = 2Φ/c₀²")
print(f"  Typical Φ/c₀² ≈ {sigma_Phi_c2:.0e}")
print(f"  σ²_ψ ≈ (2×{sigma_Phi_c2:.0e})² = {sigma2_psi_wf:.1e}")

# Jensen corrections to Friedmann:
# H² = 8πG/3 · <ρ> + Λ/3
# ρ includes contributions from ψ variance through metric
#
# Effective Λ correction from <1/ψ> ≠ 1:
# In lapse: c²_eff = c₀²·<1/ψ> ≈ c₀²(1 + σ²_ψ)
# Fractional change: Δc²/c² ≈ σ²_ψ
delta_c2 = sigma2_psi_wf
print(f"\n  Jensen corrections:")
print(f"  Δ(c²)/c² from <1/ψ> ≈ σ²_ψ = {delta_c2:.1e}")
print(f"  ΔΛ/Λ ≈ σ²_ψ = {delta_c2:.1e}")
print(f"  ΔH₀/H₀ ≈ σ²_ψ × Ω_Λ/2 = {delta_c2 * Om_L / 2:.1e}")
print(f"  → NEGLIGIBLE (need ΔH₀/H₀ ≈ 0.04)")

# ==========================================================================
# PART 4: CLUSTER-SCALE NONLINEAR CONTRIBUTIONS
# ==========================================================================
print(f"\n{'='*78}")
print(f"  PART 4: Nonlinear regime — clusters, halos, voids")
print(f"{'='*78}")

# In collapsed structures (galaxies, clusters), δ >> 1
# The potential Φ/c² is much larger than 10⁻⁵:
# Galaxy: Φ/c² ~ v²/c² ~ (200 km/s / 3×10⁵ km/s)² ~ 4×10⁻⁷ (but it's the TOTAL potential)
# Actually Φ/c² at virial radius:
# Galaxy (M~10¹² M☉, R~200 kpc): Φ = GM/R → Φ/c² ~ 10⁻⁶
# Cluster (M~10¹⁵ M☉, R~2 Mpc): Φ/c² ~ 10⁻⁵
# These are LOCAL values, not cosmic average

# Volume-weighted average:
# Fraction of volume in halos ~ 10⁻³ (halos occupy tiny volume fraction)
# Voids: 80% of volume, δ ~ -0.8, Φ/c² ~ few × 10⁻⁵
# Filaments: 15% of volume, δ ~ 1-10, Φ/c² ~ 10⁻⁵

f_halo = 1e-3   # volume fraction in halos
Phi_halo = 1e-5  # Φ/c² in halos
f_void = 0.80
Phi_void = 2e-5  # deeper void potential
f_filament = 0.15
Phi_filament = 5e-6

# Volume-weighted σ²_ψ in nonlinear regime
sigma2_psi_nl = (f_halo * (2*Phi_halo)**2 +
                  f_void * (2*Phi_void)**2 +
                  f_filament * (2*Phi_filament)**2)
print(f"\n  Volume-weighted σ²_ψ from structures:")
print(f"  Halos ({f_halo:.0e} vol):   Φ/c² ~ {Phi_halo:.0e}, contrib = {f_halo*(2*Phi_halo)**2:.1e}")
print(f"  Voids ({f_void:.2f} vol):   Φ/c² ~ {Phi_void:.0e}, contrib = {f_void*(2*Phi_void)**2:.1e}")
print(f"  Filaments ({f_filament:.2f} vol): Φ/c² ~ {Phi_filament:.0e}, contrib = {f_filament*(2*Phi_filament)**2:.1e}")
print(f"  Total σ²_ψ = {sigma2_psi_nl:.2e}")
print(f"  → Still 10⁻⁹ level — NOT enough")

# ==========================================================================
# PART 5: TEMPORAL DERIVATIVES — STRUCTURE FORMATION RATE
# ==========================================================================
print(f"\n{'='*78}")
print(f"  PART 5: ψ̇ from cosmological evolution of structures")
print(f"{'='*78}")

# Even though δψ ~ 10⁻⁵ is small, the TIME DERIVATIVE matters for backreaction.
# B_ψ = 3<ψ̇²/ψ>
#
# ψ̇ sources:
# 1. From growth of structure: δψ̇ ~ d(2Φ/c²)/dt
# 2. In ΛCDM: dΦ/dt ≈ -H·Φ·(1 - f) where f = d ln D/d ln a ≈ Ω_m^0.55
#    At z=0: f ≈ 0.315^0.55 ≈ 0.525, so (1-f) ≈ 0.475
#    Φ̇/Φ ≈ -0.475·H₀ → Φ decays in Λ-dominated era
# 3. Infall: ψ̇_infall ~ H·f·δψ on scales entering nonlinear regime
#
# Actually, the Newtonian potential is nearly constant in matter era
# and starts decaying in Λ era. The rate of change:
# dΦ/dt ~ -H₀·Φ·0.5 (rough)

Phi_over_c2 = 1e-5
dPhi_dt_H0 = 0.5 * H0 * Phi_over_c2  # H₀ units, dimensional analysis

# ψ̇ = 2Φ̇/c₀² → ψ̇ ~ 2 × 0.5 × H₀ × 10⁻⁵ = 10⁻⁵ H₀
dpsi_dt = 2 * 0.5 * Phi_over_c2  # in H₀ units
B_temporal = 3 * dpsi_dt**2  # B_ψ/H₀²
print(f"  δψ = 2Φ/c₀² ≈ {2*Phi_over_c2:.0e}")
print(f"  δψ̇ = 2Φ̇/c₀² ≈ 2 × 0.5 × H₀ × {Phi_over_c2:.0e} = {dpsi_dt:.0e} H₀")
print(f"  B_ψ/H₀² = 3(δψ̇)² = {B_temporal:.1e}")
print(f"  Required: 0.174")
print(f"  → 10⁻¹⁰ level — NOT enough")

# ==========================================================================
# PART 6: WHAT AMPLITUDE δψ IS NEEDED?
# ==========================================================================
print(f"\n{'='*78}")
print(f"  PART 6: Required δψ amplitude for H₀ tension")
print(f"{'='*78}")

# B_ψ/H₀² = 0.174 needed
# B_ψ = 3<ψ̇²/ψ> ≈ 3<ψ̇²>
# If ψ oscillates with frequency ω and amplitude A:
#   <ψ̇²> ≈ ω²·A²/2
# For tachyonic mode: ω = m_tach·H₀ ≈ 5H₀
#   B_ψ/H₀² ≈ 3 × m² × A²/2 = 3 × 24.7 × A²/2

B_required = 0.174
A_required = math.sqrt(2 * B_required / (3 * gamma_nominal))
print(f"  If ψ oscillates at tachyonic frequency m = {m_tach_nominal:.1f}·H₀:")
print(f"  B_ψ/H₀² = 3·m²·A²/2 = 0.174")
print(f"  A = √(2·B/(3·m²)) = {A_required:.4f}")
print(f"  → Need δψ ~ {A_required:.2e} (i.e., ψ deviates from 1 by {A_required*100:.2f}%)")

# If ψ̇ is set by Hubble timescale instead:
# <ψ̇²> ≈ H₀²·<δψ²>
# B_ψ/H₀² = 3·<δψ²>
A_hubble = math.sqrt(B_required / 3)
print(f"\n  If ψ̇ ~ H₀·δψ (Hubble timescale):")
print(f"  B_ψ/H₀² = 3·<δψ²> = 0.174")
print(f"  δψ_rms = {A_hubble:.4f}")
print(f"  → Need δψ ~ {A_hubble:.1e} (~{A_hubble*100:.0f}% deviation)")

# Compare with what we get from gravity:
print(f"\n  Available from gravity: δψ ~ 2Φ/c₀² ~ 2×10⁻⁵")
print(f"  GAP (tachyonic freq): {A_required/(2e-5):.0f}×")
print(f"  GAP (Hubble freq):    {A_hubble/(2e-5):.0f}×")

# ==========================================================================
# PART 7: ALTERNATIVE — SOLITON POPULATION COSMOLOGY
# ==========================================================================
print(f"\n{'='*78}")
print(f"  PART 7: Soliton population as effective dark energy component")
print(f"{'='*78}")

# In TGP, solitons are fundamental objects (particles).
# A gas of solitons has its own equation of state.
# If solitons have mass m_sol and the "gas" is relativistic or has pressure:
#
# For non-relativistic soliton gas: w = 0 (like CDM)
# For relativistic solitons: w = 1/3 (like radiation)
# For soliton condensate: w could be anything
#
# The KEY question: does a soliton gas contribute to effective Λ?
#
# Soliton has excess volume: δV_sol = ∫(ψ^(3/2) - 1) d³x
# If each soliton contributes δV > 0, then a population of solitons
# modifies the effective spatial curvature → effective Λ

# From soliton analysis (δ_crit = 1.206188):
# For soliton with g₀ = 2 (so δ = g₀ - 1 = 1):
# ψ = g₀ at center, ψ → 1 at infinity
# Excess metric volume: ∫(ψ^(3/2) - 1) 4πr² dr ≈ volume integral of tail

# This is a DIFFERENT mechanism: not backreaction of perturbations,
# but TOPOLOGICAL contribution of solitons to cosmological expansion.

# Simple estimate:
# Each soliton with g₀ ≈ 2 has "core" radius r_c ~ 1/m (in soliton units)
# Volume per soliton: V_sol ~ (4π/3) r_c³
# Excess spatial volume fraction: δV/V ~ n_sol × V_sol × (g₀^(3/2) - 1)
# where n_sol = number density of solitons

# If solitons ARE dark matter particles:
# ρ_DM = n_sol × m_sol → n_sol = ρ_DM / m_sol
# ρ_DM = Ω_DM × ρ_crit = 0.266 × 3H₀²/(8πG)

# The question reduces to: what is the soliton mass in physical units?
# This depends on the TGP length scale ℓ_TGP = c₀/√γ = c₀²/(m_tach·H₀·c₀)

# In Planck units: if ℓ_TGP = ℓ_Planck, then m_sol ~ m_Planck
# But ℓ_TGP = c₀/(√γ · H₀/c₀ · c₀) ... needs care

# Let's be concrete about what we know:
# In TGP units: soliton is a profile ψ(r) with r in units of some length ℓ
# The soliton equation: g'' + g'²/g + 2g'/r + g = 1
# has NO free length scale → the length scale must come from the potential γ

# In physical units: the "natural" length is ℓ = c₀/√(γ·c₀²/H₀²) × (c₀/H₀)
# Wait. Let me be more careful.

# The TGP equation in physical units:
# ψ̈ + 3Hψ̇ + 3ψ̇²/ψ = c₀²[∇²ψ/a² + γ(ψ-1)]
# For static soliton (ψ̈ = 0, Hψ̇ = 0):
# c₀²[∇²ψ + γ(ψ-1)] = 0
# ∇²ψ = -γ(ψ-1)
# With r in physical units: if we define r' = r·√γ, we get
# ∇'²ψ = -(ψ-1)
# So soliton size in physical units: R_sol = 1/√γ (in c₀/H₀ units)
# R_sol_physical = c₀/(H₀·√γ) = c₀/(H₀·m_tach) = 4451/4.97 = 895 Mpc

R_sol_physical = c_over_H0 / m_tach_nominal
print(f"  Soliton radius (from γ): R_sol = c₀/(H₀·m_tach) = {R_sol_physical:.0f} Mpc")
print(f"  → COSMOLOGICAL size solitons!")
print(f"  → NOT particle-like dark matter candidates at this γ")

# But the soliton equation from the brannen_sqrt2 analysis uses
# a DIFFERENT length scale — the microscopic TGP lattice scale.
# The equation g'' + g'²/g + 2g'/r + g = 1 has r in LATTICE units.
# The γ in the cosmological equation is the COSMOLOGICAL coupling,
# which may be completely different from the microscopic soliton equation.

print(f"""
  KEY INSIGHT: There are TWO scales in TGP:
  1. MICROSCOPIC: soliton equation g'' + g'²/g + 2g'/r + g = 1
     Length scale: ℓ_micro (lattice spacing or Planck length)
     Soliton radius: r_sol ~ few × ℓ_micro
     These solitons are PARTICLES (dark matter candidates)

  2. COSMOLOGICAL: ψ̈ + 3Hψ̇ + 3ψ̇²/ψ = c₀²[∇²ψ + γ(ψ-1)]
     Length scale: ℓ_cosmo = c₀/(H₀·√γ) ~ 900 Mpc
     This is the scale where substrate HOMOGENEITY breaks down

  The confusion was treating the cosmological γ as if it controlled
  microscopic soliton physics. These are SEPARATE regimes.

  The cosmological equation's γ gives Λ.
  The microscopic equation has its OWN coupling → particle masses.

  They are connected but NOT the same equation with the same γ.""")

# ==========================================================================
# PART 8: THE REAL BACKREACTION — BUCHERT FRAMEWORK
# ==========================================================================
print(f"\n{'='*78}")
print(f"  PART 8: Buchert-type backreaction in TGP")
print(f"{'='*78}")

# In GR, the Buchert framework shows that spatial averaging of
# an inhomogeneous universe gives effective Friedmann equations
# with backreaction terms Q_D and <R>_D.
#
# Standard GR backreaction is ~10⁻⁵ (Newtonian limit).
# The question: does TGP's nonlinear structure ENHANCE this?

# In TGP, the spatial metric is ψ·δᵢⱼ.
# The scalar curvature of a slice: R = 0 for flat spatial metric (δᵢⱼ)
# BUT ψ·δᵢⱼ has curvature from ψ gradients!
# R^(3) = -2(∇²ψ)/(ψ²) + (∇ψ)²/(2ψ³)  [for conformally flat metric]
# Actually for gᵢⱼ = ψ δᵢⱼ (3D conformal metric):
# R^(3) = -2∇²(ln ψ)/ψ + ... (complicated)

# Let me use the conformal factor properly:
# gᵢⱼ = e^{2α} δᵢⱼ where e^{2α} = ψ → α = ln(ψ)/2
# R^(3) = -4(∇²α + (∇α)²)/e^{2α}
# For ψ ≈ 1 + δψ, α ≈ δψ/2:
# R^(3) ≈ -4(∇²(δψ/2)) = -2∇²(δψ)

# <R^(3)> = -2<∇²(δψ)> = 0 (by periodicity/boundary conditions)
# But <R²> ≠ 0 → <R² - <R>²> contributes to backreaction

# Buchert backreaction for TGP:
# Q_D = 2/3(<θ²> - <θ>²) - 2<σ²>
# where θ = 3H_local, σ = shear

# In TGP: H_local depends on ψ through lapse and expansion:
# H_local = ȧ/a + ψ̇/(2ψ)  [additional expansion from substrate dilation]
# <θ> = 3<H_local> = 3(H + <ψ̇/(2ψ)>) ≈ 3H (if <ψ̇> = 0)
# <θ²> - <θ>² = 9·var(H_local) = 9·<(ψ̇/(2ψ))²>
#              ≈ 9/4 · <ψ̇²> (for ψ ≈ 1)

# This gives: Q_D = 2/3 × 9/4 × <ψ̇²> = 3/2 × <ψ̇²>
# Compare with our B_ψ = 3<ψ̇²/ψ> ≈ 3<ψ̇²>: Q_D = B_ψ/2

dpsi_dot_rms = 1e-5  # H₀ units (from δψ̇ ~ H·δψ ~ H·10⁻⁵)
Q_D = 1.5 * dpsi_dot_rms**2
print(f"  Buchert backreaction in TGP:")
print(f"  Q_D = (3/2)·<ψ̇²> ≈ (3/2)·(H₀·δψ)² = {Q_D:.1e} H₀²")
print(f"  Effective: ΔH²/H² ≈ Q_D/(3H₀²) = {Q_D/3:.1e}")
print(f"  → Same order as before: 10⁻¹⁰, NEGLIGIBLE")

# ==========================================================================
# PART 9: DOES γ HAVE TO BE 36·Ω_Λ?
# ==========================================================================
print(f"\n{'='*78}")
print(f"  PART 9: Re-examining the γ–Λ relationship")
print(f"{'='*78}")

# The claim Λ_eff = γ/12 needs careful re-derivation.
# Let's check what VALUE of γ gives the observed dark energy.
#
# The TGP action (from literature context):
# S_TGP = ∫ √(-g) [R/(16πG) + L_ψ] d⁴x
# L_ψ = ψ⁴/(2c₀²) (∂ψ)² + V(ψ)    [K(ψ) = ψ⁴ kinetic coupling]
# V(ψ) = γ(ψ-1)²/2     [tachyonic potential]
#
# At ψ = 1 (vacuum): L_ψ = V(1) = 0 → NO vacuum energy from potential
# → Λ is NOT from V(ψ=1) but from the ZERO-POINT energy of fluctuations!
#
# Alternatively: the soliton equation g'' + g'²/g + 2g'/r + g = 1
# has g = 1 as trivial solution. The "+g = 1" term means:
# -g + 1 = -(g-1) → restoring toward g=1 (BUT combined with nonlinear kinetic
# term g'²/g, the stability is subtle)
#
# For the cosmological equation, the term γ(ψ-1) in:
# ψ̈ + 3Hψ̇ + 3ψ̇²/ψ = c₀²γ(ψ-1)
# RHS = 0 at ψ = 1 → no cosmological constant from this term alone.

print(f"""
  Re-examination of Λ_eff origin in TGP:

  The potential V(ψ) = γ(ψ-1)²/2 gives V(1) = 0.
  → NO classical vacuum energy at ψ = 1.

  The cosmological constant in TGP must come from:
  A) Zero-point fluctuations of ψ field → quantum calculation needed
  B) The constant term in g + ... = 1 → the "1" represents spatial flatness
  C) A separate cosmological constant NOT related to γ
  D) Boundary conditions in the soliton lattice

  If Λ and γ are INDEPENDENT parameters, then:
  - γ could be much larger than 36·Ω_Λ
  - Tachyonic instability could operate at shorter scales
  - The backreaction mechanism MIGHT work after all

  BUT: the relationship Λ_eff = γ/12 was derived from the TGP action.
  Let's check if there's an alternative derivation...""")

# The soliton equation: g'' + (1/g)(g')² + (2/r)g' + g = 1
# In "cosmological" form (FRW, spatial part):
# The "+g" term becomes "+γψ" in the cosmological equation
# The "= 1" becomes "= γ·1"
# So γ enters as the coupling of ψ to spatial curvature
# And the vacuum ψ = 1 is set by the "= 1" (or "= γ/γ = 1")

# Actually, the original equation has coefficient 1 for both g and 1:
# g'' + g'²/g + 2g'/r = 1 - g = -(g - 1)
# So the "mass²" of fluctuations around g=1 is -1 (in soliton units)
# The physical γ is this mass² converted to physical units:
# γ = (1/ℓ_micro²) in microscopic units, or
# γ = c₀²/ℓ_micro² in c₀,H₀ units

# The KEY question: is ℓ_micro the SAME as ℓ_cosmo?
# If TGP has ONE length scale: ℓ_micro = ℓ_cosmo → γ determined by Λ → super-Hubble
# If TGP has SEPARATE scales: γ_micro ≠ γ_cosmo → different physics

print(f"""
  CRITICAL QUESTION: Does TGP have one or two scales?

  ONE scale: ℓ = c₀/√γ determines BOTH soliton size AND Λ
    → γ fixed by Λ_obs → γ ≈ 25 (H₀/c₀)² → ℓ ~ 900 Mpc
    → Solitons are 900 Mpc objects (NOT particles!)
    → Tachyonic instability is super-Hubble
    → Backreaction mechanism FAILS

  TWO scales: ℓ_micro (soliton/particle) and ℓ_cosmo (expansion)
    → γ_micro >> γ_cosmo
    → Solitons are microscopic (particles) ✓
    → But tachyonic instability uses γ_cosmo (still super-Hubble)
    → Backreaction mechanism still fails, but for different reason
    → HOWEVER: the two γ's might couple → effective γ_eff on intermediate scales

  This is a fundamental architectural question for TGP.""")

# ==========================================================================
# PART 10: QUANTITATIVE SUMMARY — ALL MECHANISMS
# ==========================================================================
print(f"\n{'='*78}")
print(f"  PART 10: QUANTITATIVE SUMMARY — ALL MECHANISMS")
print(f"{'='*78}")

# Collect all computed B_ψ/H₀² values
mechanisms = [
    ("Standard GR backreaction (Buchert)", 1e-10, "From <δ²> ~ 10⁻⁵"),
    ("TGP weak-field perturbative", 1e-9, "From σ²_ψ ~ (2Φ/c₀²)²"),
    ("Tachyonic transfer function (ct5)", 1.04e-9, "Full k-integral with tachyonic pole"),
    ("Temporal derivative (dψ/dt)", 3e-10, "From Φ̇ ~ H₀·Φ"),
    ("Jensen inequality (metric)", 4e-10, "From <1/ψ> - 1/<ψ>"),
    ("Buchert Q_D in TGP", 1.5e-10, "From var(θ_local)"),
    ("Tachyonic saturation estimate (ct3)", 0.03, "Optimistic: m_tach×exp(-damping)"),
    ("Required for H₀ tension", 0.174, "67.4 → 73.0"),
]

print(f"\n  {'Mechanism':<45s}  B_ψ/H₀²      Status")
print(f"  {'—'*75}")
for name, val, note in mechanisms:
    if val > 0.01:
        status = "✓ SUFFICIENT" if val >= 0.17 else "~ MARGINAL"
    elif val > 1e-5:
        status = "✗ too small"
    else:
        status = "✗✗ NEGLIGIBLE"
    print(f"  {name:<45s}  {val:.2e}    {status}")
    print(f"  {'':45s}  ({note})")

print(f"""

  ═══════════════════════════════════════════════════════════════════════
  FINAL ASSESSMENT
  ═══════════════════════════════════════════════════════════════════════

  All RIGOROUS calculations give B_ψ/H₀² ~ 10⁻¹⁰ to 10⁻⁹.
  The ct3 "optimistic" estimate of 0.03 was based on a HOMOGENEOUS
  analysis that ignored the spatial gradient suppression.

  When spatial gradients are properly included (ct5), the tachyonic
  mode only enhances super-Hubble scales where there's no significant
  perturbation source → effect is negligible.

  STATUS: The perturbative backreaction mechanism from 3ψ̇²/ψ
  DOES NOT produce sufficient effect to explain H₀ tension.

  REMAINING AVENUES:

  1. NON-PERTURBATIVE: Soliton population changes effective expansion
     (requires knowing soliton mass in physical units)

  2. RUNNING Λ: If γ runs with scale/energy → Λ_eff(z) ≠ const
     (requires RG analysis of TGP)

  3. SEPARATE SCALES: If γ_micro ≠ γ_cosmo, microscopic soliton
     physics might affect macroscopic expansion through collective effects

  4. MODIFIED DISPERSION: TGP lattice structure → modified photon
     propagation → direct H₀ measurement affected (not actual expansion)

  5. COUPLED DARK SECTOR: If soliton mass depends on ψ_background,
     and ψ_background evolves → effective dark matter mass changes
     → modified growth history → S₈ and w(z) affected

  6. RE-EXAMINE TGP ACTION: Perhaps the K=ψ⁴ kinetic coupling
     was too specific. Other K(ψ) choices might give stronger backreaction.

  The HONEST conclusion: the simple "tachyonic backreaction resolves
  H₀ tension" narrative does not survive quantitative scrutiny.
  TGP's cosmological implications need a DIFFERENT mechanism.
  ═══════════════════════════════════════════════════════════════════════
""")
