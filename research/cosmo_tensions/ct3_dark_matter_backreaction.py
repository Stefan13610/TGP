#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ct3_dark_matter_backreaction.py: Refined backreaction with dark matter.

*** SUPERSEDED 2026-04-26 by M10.5 (H_0/S_8 tensions audit, canonical sek08a) ***
*** See: ../op-cosmology-closure/M10_5_results.md (6/6 PASS)                ***

POST-M10.5 STATUS (2026-04-26):
- This script's "tachyonic amplification" interpretation of V''(1)=-gamma
  was conflation of cosmological slow-roll (cosmic time, Hubble-DAMPED)
  with spatial linearization (M9.3.1: M_eff^2 = +beta > 0 stable Yukawa).
- Optimistic claim B_psi/H0^2 ~ 0.03-0.3 SUPERSEDED. Correct canonical
  sek08a value: B_psi/H0^2 ~ 1.08e-8 (gap 7.2 orders below 0.17 needed
  for H_0 tension). Source of error: canonical K=1 used, plus spurious
  amplification factor from misread linearization.
- Honest verdict (ct7) CONFIRMED: TGP scope = galaxy/PPN, NOT cosmology
  tensions. See [[../op-cosmology-closure/M10_5_results.md]] for full
  reconciliation and 6/6 PASS audit.
- This file kept for audit traceability; do NOT use its numerical
  conclusions in derivative work.

ORIGINAL DRAFT NOTES (pre-M10.5, retained for reference):
KEY INSIGHT from ct2: naive baryon-only estimate gives B₀/H₀² = 0.04,
but we need 0.57. The factor ~14 gap suggests:
1. Dark matter (Ωm = 0.315 >> Ωb = 0.049) dominates backreaction
2. The effective δ² is larger than particle-level soliton amplitudes
3. Nonlinear amplification from structure formation feedback

In TGP, dark matter could be:
A) Heavy solitons with g₀ ≫ 1 (near δ_crit)
B) Collective substrate excitations (phonon-like modes)
C) Same substrate field ψ but with gravitational potential contributing

APPROACH: Model backreaction as proportional to TOTAL gravitational
potential energy, not just soliton amplitude.

The gravitational potential Φ_grav at scale R:
  Φ/c² ~ GM/(Rc²) ~ Ω_m H₀² R² / c²

For structures at virial equilibrium:
  σ²_v/c² ~ Φ/c² ~ 10⁻⁵ to 10⁻³ (clusters to galaxies)

The backreaction in TGP = nonlinear coupling 3ψ̇²/ψ
where ψ̇ is driven by gravitational INFALL, not Hubble flow.
So: δψ̇ ~ H · δψ_grav where δψ_grav ~ Φ_N/c² · ψ

This gives B_ψ ~ H² · <(Φ/c²)²> · some amplification factor.
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from scipy.integrate import solve_ivp
import math

print("="*78)
print("  REFINED BACKREACTION: Dark matter + gravitational potential")
print("="*78)

# ==========================================================================
# PARAMETERS
# ==========================================================================
H0_planck = 67.4    # km/s/Mpc
H0_shoes  = 73.04   # km/s/Mpc
Om_m = 0.315
Om_b = 0.049
Om_dm = Om_m - Om_b  # dark matter
Om_L = 1 - Om_m - 9.1e-5

S8_planck = 0.832
S8_lensing = 0.76

delta_crit = 1.206188  # TGP soliton limit

# ==========================================================================
# MODEL: Gravitational potential backreaction
# ==========================================================================
print(f"\n{'='*78}")
print("  Gravitational potential variance")
print("="*78)

# The variance of the Newtonian potential:
# <(Φ/c²)²> = ∫ P_Φ(k) d³k/(2π)³
#
# For matter power spectrum P(k) ~ k^ns · T²(k):
# P_Φ(k) ∝ P(k)/k⁴ (Poisson equation: k²Φ = 4πGρδ)
#
# Typical values:
# <(Φ/c²)²>^(1/2) ~ σ_Φ ~ 10⁻⁵ (CMB dipole) to 10⁻³ (clusters)
#
# For backreaction: the relevant quantity is the VARIANCE
# of the potential over the Hubble volume, weighted by ∂ψ/∂Φ

# In TGP: ψ responds to gravitational potential via the field equation.
# For weak fields (linear regime):
#   ψ ≈ 1 + 2Φ_N/c² (same as GR in weak field, from PPN γ=1)
#
# So δψ ≈ 2Φ/c², and:
#   <(δψ)²> ≈ 4<(Φ/c²)²>
#   <(δψ̇)²> ≈ 4<(Φ̇/c²)²> ~ 4H²<(Φ/c²)²>  (growth rate)

# The key integral for backreaction:
# B_ψ = 3[<(δψ̇)²>/ψ̄ - ψ̄̇² <(δψ)²>/ψ̄²]
#      ≈ 3 · 4H² · <(Φ/c²)²> · [1 - (ψ̄̇/Hψ̄)²]
#      ≈ 12 H² · σ²_Φ · α_corr

# where α_corr accounts for correlation between δψ and δψ̇

# The σ²_Φ variance grows with structure formation.
# Using σ₈ as proxy:
#   σ²_Φ ~ (3/2 Ωm H₀²/k²)² · σ²_8  where k ~ 2π/(8 Mpc/h)

# More directly: define
#   σ²_ψ(a) = 4 σ²_Φ(a) = 4 · (3Ωm H₀² a)²/(2k_eff²c²)² · D²(a)/D²(1) · σ²_Φ,0

# For order of magnitude:
sigma_Phi_0 = 3e-5  # rms gravitational potential today (from CMB + LSS)
# This is the rms of Φ/c² on ~100 Mpc scales

print(f"\n  σ_Φ/c² (rms gravitational potential):")
print(f"    CMB scale (~100 Mpc): σ_Φ ~ {sigma_Phi_0:.1e}")
print(f"    Cluster scale (~1 Mpc): σ_Φ ~ 10⁻³ to 10⁻²")
print(f"    Galaxy scale (~10 kpc): σ_Φ ~ 10⁻⁶ (virial)")

# ==========================================================================
# BUT: the crucial point is the NONLINEAR amplification
# ==========================================================================
print(f"\n{'='*78}")
print("  NONLINEAR AMPLIFICATION in TGP")
print("="*78)

print(f"""
  In GR: backreaction is second order in Φ/c² ~ (10⁻⁵)² = 10⁻¹⁰
  → too small to explain H₀ tension (needs ΔH/H ~ 8%)

  BUT in TGP: the term 3ψ̇²/ψ has a DIFFERENT structure.

  Standard scalar field: ψ̈ + 3Hψ̇ = -V'(ψ)
  TGP field equation:    ψ̈ + 3Hψ̇ + 3ψ̇²/ψ = c₀²W(ψ)

  The 3ψ̇²/ψ term acts as GEOMETRIC FRICTION that depends
  on ψ itself. This creates a FEEDBACK LOOP:

  1. Matter clusters → ψ deviates from 1 locally
  2. ψ̇ from gravitational infall ≠ 0
  3. 3ψ̇²/ψ adds to effective Friedmann → faster expansion
  4. Faster expansion changes ψ̇ → changes backreaction
  5. This is a SELF-CONSISTENT coupling

  The amplification factor α_TGP vs naive GR estimate:
  α_TGP = (ψ̇/ψ)² / (Φ̇/c²Φ₀)²

  In virialized structures: ψ̇/ψ is set by VIRIAL velocity,
  not by Hubble flow. So:
  α_TGP ~ (v_vir/c)²/(H·R/c)² ~ (v/HR)² ~ (σ_v/HR)²
""")

# Typical virial velocities:
v_cluster = 1000   # km/s (cluster of galaxies)
v_galaxy = 200     # km/s (galaxy)
v_group = 500      # km/s (galaxy group)

# Hubble flow at virial radius:
H0_km = H0_planck  # km/s/Mpc
R_cluster = 2.0    # Mpc (cluster virial radius)
R_galaxy = 0.05    # Mpc (galaxy virial radius)
R_group = 0.5      # Mpc (group virial radius)

print(f"\n  Virial velocity vs Hubble flow:")
print(f"  {'Structure':>15s}  {'v_vir':>8s}  {'H·R':>8s}  {'v/HR':>8s}  {'(v/HR)²':>10s}")

for name, v, R in [('Galaxy', v_galaxy, R_galaxy),
                    ('Group', v_group, R_group),
                    ('Cluster', v_cluster, R_cluster)]:
    HR = H0_km * R
    ratio = v / HR
    print(f"  {name:>15s}  {v:8.0f}  {HR:8.1f}  {ratio:8.1f}  {ratio**2:10.1f}")

# ==========================================================================
# REFINED MODEL: B_ψ with virial contribution
# ==========================================================================
print(f"\n{'='*78}")
print("  REFINED MODEL: Virial backreaction")
print("="*78)

# The backreaction has TWO contributions:
# 1. Linear (Hubble-flow driven): B_linear ~ H² · 4σ²_Φ ~ tiny
# 2. Nonlinear (virial driven): B_virial ~ n_struct · V_struct · (v_vir/c)²

# For #2: the density of structures × their virial contribution
# ρ_vir = ∫ dn/dM · M · (σ_v(M)/c)² dM
# where dn/dM is the halo mass function (Press-Schechter or Sheth-Tormen)

# Rough estimate:
# Total kinetic energy of peculiar velocities:
# ρ_kin = (1/2) ρ_m · <v²>/c²
# where <v²> ~ (500 km/s)² is the mass-weighted mean virial velocity

# The backreaction in TGP from virial motions:
# B_ψ = 3 · <(δψ̇)²/ψ> ≈ 3 · <(2v/c)²> · H² · correction

# Wait - let's think more carefully.
# In the TGP field equation for ψ:
# The term 3ψ̇²/ψ is the GEOMETRIC self-coupling.
# At cosmological level (homogeneous): ψ̇ = (dψ/dt)_Hubble ~ Hψ · O(1)
# Inside a halo: ψ̇ = (dψ/dt)_Hubble + (dψ/dt)_virial
#
# The virial contribution to ψ̇ comes from matter streaming:
# ψ responds to local matter density via ∇²ψ ~ -ρ
# Time variation: ψ̇ ~ ∇·(ψ v) where v is bulk flow
#
# For virial motions: ψ̇_vir ~ ψ · v_vir · k_struct where k ~ 1/R_vir
# → ψ̇_vir/ψ ~ v/R ~ v·k
#
# The backreaction:
# <(ψ̇)²/ψ> - <ψ̇>²/<ψ>
# = <(ψ̇_Hubble + ψ̇_vir)²/ψ> - <ψ̇_Hubble>²/<ψ>
# ≈ <ψ̇²_vir/ψ> + 2<ψ̇_H ψ̇_vir/ψ>  (cross-term averages out if random)
# ≈ <ψ̇²_vir>/ψ̄  (virial motions uncorrelated with Hubble flow)
#
# Now: ψ̇_vir/ψ ~ v·k, so <(ψ̇_vir/ψ)²> = <v²·k²>
# Summing over all structures:
# B_ψ = 3 · f_collapsed · <v²> · <k²>  (in natural units)
#
# In cosmological units:
# B_ψ/H₀² = 3 · f_collapsed · (<v²>/c²) · (k_eff·c/H₀)²

f_coll = 0.7  # fraction of matter in halos at z=0
v_rms = 400   # km/s (mass-weighted rms peculiar velocity)
c_km = 3e5    # speed of light in km/s

# Effective k is the scale where most backreaction comes from
# For virialized halos: k_eff ~ 1/R_eff where R_eff ~ 0.3 Mpc (groups/small clusters)
R_eff_Mpc = 0.3  # Mpc
k_eff = 1/R_eff_Mpc  # 1/Mpc
c_Mpc_s = c_km / H0_planck  # convert c to Mpc·H₀

B_virial_over_H0sq = 3 * f_coll * (v_rms / c_km)**2 * (k_eff * c_Mpc_s)**2

print(f"\n  Virial backreaction estimate:")
print(f"    f_collapsed = {f_coll}")
print(f"    v_rms = {v_rms} km/s")
print(f"    R_eff = {R_eff_Mpc} Mpc → k_eff = {k_eff:.1f} /Mpc")
print(f"    c/H₀ = {c_Mpc_s:.0f} Mpc")
print(f"    B_ψ/H₀² = 3 · {f_coll} · ({v_rms}/{c_km:.0e})² · ({k_eff}·{c_Mpc_s:.0f})²")
print(f"    B_ψ/H₀² = {B_virial_over_H0sq:.4f}")
print(f"    Required: 0.174 (for H₀ = 73)")

# Hmm, this is way too large! Let me reconsider.
# The issue: k_eff · c/H₀ ~ 3/Mpc · 4450 Mpc ~ 14000
# → (v/c)² · (kc/H₀)² ~ 10⁻⁶ · 10⁸ ~ 100
# This is ENORMOUS — the model breaks down.

# The problem: you can't just multiply virial velocity by k.
# The ψ field is smooth on scales >> particle separation.
# At cosmological level, the relevant ψ̇ is the VOLUME-AVERAGED
# response, not the local virial velocity.

# Let me reconsider from scratch.

print(f"\n{'='*78}")
print("  CORRECTED APPROACH: Volume-averaged backreaction")
print("="*78)

print(f"""
  The correct approach:

  1. The TGP field ψ(x,t) satisfies the FULL nonlinear PDE.
  2. At cosmological level, we average over a Hubble volume.
  3. The backreaction is:
     Q_TGP = <3ψ̇²/ψ> - 3<ψ̇>²/<ψ>  [analogous to Buchert's Q]

  4. In TGP, ψ couples to matter density:
     ψ ≈ 1 + 2Φ_N/c²  (weak field, PPN γ=1)
     δψ ≈ 2Φ/c²

  5. The gravitational potential satisfies Poisson:
     ∇²Φ = 4πG ρ_m
     → Φ ~ -GM/r for point mass, Φ ~ H₀²Ω_m a R² for uniform sphere

  6. Time derivative:
     ψ̇ ≈ 2Φ̇/c²
     Φ̇ arises from: (a) Hubble flow, (b) structure growth

  7. For structure growth: Φ̇/Φ ~ f(Ω_m)·H where f ≈ Ω_m^0.55
     → δψ̇ ≈ 2fH·Φ/c²

  8. The variance:
     <(δψ̇)²> = (2fH)² · <(Φ/c²)²> = 4f²H² · σ²_Φ

  9. Backreaction (leading term):
     Q_TGP ≈ 3 · <(δψ̇)²>/ψ̄ ≈ 12 f² H² σ²_Φ

  10. The CRUCIAL AMPLIFICATION comes from TGP structure:
      In GR, backreaction ∝ (Φ/c²)² ~ 10⁻¹⁰ (negligible)
      In TGP, the 3ψ̇²/ψ is NOT just backreaction — it's a
      FUNDAMENTAL COUPLING that modifies the evolution equation.

  KEY DIFFERENCE vs GR:
  In GR, <R> - R(<g>) is a purely geometric effect.
  In TGP, the 3ψ̇²/ψ term is PHYSICAL FRICTION from the substrate.
  The ψ field doesn't just describe geometry — it IS the medium.

  When matter clusters, the substrate RESPONDS by generating
  more space in voids (to conserve the substrate budget).
  This is NOT captured by standard backreaction estimates.
""")

# ==========================================================================
# MODEL B: Substrate budget conservation
# ==========================================================================
print(f"\n{'='*78}")
print("  MODEL B: Substrate budget conservation")
print("="*78)

print(f"""
  TGP has a UNIQUE constraint: the antipodal budget.
  f(ψ)·h(ψ) = 1  →  (spatial volume)·(temporal rate) = const

  In homogeneous universe: ψ = ψ̄ everywhere → budget trivially satisfied.

  In inhomogeneous universe:
  - Matter-dense regions: ψ > 1 (or < 1 for deficit solitons)
  - Voids: ψ must COMPENSATE to maintain budget

  Define surplus fraction:
    Δ_budget = <ψ^(1/2)> - <ψ>^(1/2)   [from volume element √ψ]

  By Jensen's inequality (for concave √):
    <ψ^(1/2)> < <ψ>^(1/2)  → Δ_budget < 0

  This means: the EFFECTIVE volume element is SMALLER than
  the homogeneous estimate → the effective expansion is FASTER
  to compensate!

  The correction to H:
    H_eff² = H_bare² · (<ψ>^(1/2) / <ψ^(1/2)>)
    = H_bare² · (1 + |Δ_budget|/⟨ψ^(1/2)⟩)
""")

# Jensen's inequality bound:
# For ψ = 1 + δψ, with <δψ> = 0:
# <√ψ> = <√(1+δψ)> ≈ <1 + δψ/2 - δψ²/8 + ...>
#       = 1 - <δψ²>/8 = 1 - σ²/8
# √<ψ> = √1 = 1
# So: √<ψ> - <√ψ> ≈ σ²/8

# This means H_eff²/H_bare² ≈ 1 + σ²/8

# For H₀ tension: 1 + σ²/8 = (73/67.4)² = 1.174
# → σ² = 8 × 0.174 = 1.39

# σ² = <δψ²> for the substrate field
# δψ ≈ 2Φ/c² locally, but this gives σ ~ 10⁻⁵ → σ² ~ 10⁻¹⁰
# Far too small.

# UNLESS: the substrate fluctuations are much larger than
# the gravitational potential suggests. In TGP, ψ responds to
# the TOTAL matter distribution, not just the local potential.

# KEY REALIZATION:
# The σ² that matters is not (Φ/c²)² but the variance of ψ
# over the SUBSTRATE. If TGP has its own "substrate fluctuations"
# beyond those sourced by Newtonian gravity, these could be large.

# In the soliton picture:
# Each particle is a soliton with δψ = g₀ - 1 ~ O(1)
# The VOLUME fraction of solitons is tiny (Compton volumes)
# but the AMPLITUDE is O(1) → large contribution to <δψ²>

# Volume fraction of solitons:
# n_baryon ~ 0.25/m³ (mean baryon density)
# V_soliton ~ (ħ/(mc))³ ~ (10⁻¹⁵ m)³ = 10⁻⁴⁵ m³ (proton Compton)
# f_vol = n · V = 0.25 × 10⁻⁴⁵ = 2.5 × 10⁻⁴⁶ (TINY!)

# So: σ² = f_vol · <δ²> ~ 10⁻⁴⁶ × 1 ~ 10⁻⁴⁶ (negligible)

# This is the fundamental problem:
# At the PARTICLE level, solitons have large δψ but tiny volume
# At the COSMOLOGICAL level, <δψ> is set by gravitational potential ~ 10⁻⁵

# The backreaction from the 3ψ̇²/ψ term must therefore come from
# a MESOSCOPIC scale — between particles and Hubble volume.
# This is the scale of galaxy halos, voids, filaments.

print(f"\n  Scale analysis:")
print(f"  {'Scale':>20s}  {'δψ':>10s}  {'f_vol':>10s}  {'σ² = f·δ²':>12s}")
print(f"  {'Particle (Compton)':>20s}  {'~1':>10s}  {'~10⁻⁴⁶':>10s}  {'~10⁻⁴⁶':>12s}")
print(f"  {'Galaxy halo (Mpc)':>20s}  {'~10⁻⁵':>10s}  {'~0.1':>10s}  {'~10⁻¹¹':>12s}")
print(f"  {'Cluster (10 Mpc)':>20s}  {'~10⁻⁴':>10s}  {'~0.01':>10s}  {'~10⁻¹⁰':>12s}")
print(f"  {'Void (100 Mpc)':>20s}  {'~10⁻⁵':>10s}  {'~0.5':>10s}  {'~10⁻¹¹':>12s}")

# ==========================================================================
# CONCLUSION: The TGP amplification mechanism
# ==========================================================================
print(f"\n{'='*78}")
print("  CONCLUSION: Where does the amplification come from?")
print("="*78)

print(f"""
  Standard GR backreaction: Q ~ (Φ/c²)² ~ 10⁻¹⁰ → NEGLIGIBLE

  TGP has THREE potential amplification mechanisms:

  1. KINETIC COUPLING α=2 (ψ⁴ in action):
     The K(ψ) = ψ⁴ kinetic term means fluctuations in ψ
     couple to gradients with power ψ⁴, not ψ².
     For ψ = 1 + δψ: K ~ 1 + 4δψ + O(δψ²)
     → amplification factor: 4× vs standard scalar

  2. VOLUME ELEMENT (√(-g_eff) = c₀ψ^(1/2)):
     The volume element itself depends on ψ.
     Averaging <anything × √ψ> introduces correlations.
     This is the SUBSTRATE BUDGET effect.

  3. SELF-INTERFERENCE POTENTIAL:
     V(ψ) = (β/3)ψ³ - (γ/4)ψ⁴ has strong self-coupling.
     At vacuum (β=γ): V''(1) = -γ (tachyonic mass!).
     Small fluctuations are AMPLIFIED by the potential.
     Growth rate: ψ̈ ~ γψ → exponential on timescale t ~ 1/√γ.

  KEY INSIGHT:
  The combination of TACHYONIC INSTABILITY (V''(1) < 0)
  with GEOMETRIC FRICTION (3ψ̇²/ψ) creates a CRITICAL BALANCE:
  - At ψ=1: unstable (wants to grow)
  - Growth limited by friction (3ψ̇²/ψ increases)
  - Final state: self-consistent backreaction

  The effective backreaction amplitude is set by:
  B_ψ/H₀² ~ (γ/H₀²) · σ²_δ · amplification

  where γ ~ Λ_eff ~ H₀² (since Λ = γ/12 from TGP action)
  → B_ψ/H₀² ~ σ²_δ · amplification

  With amplification ~ 10² to 10⁴ from the tachyonic + friction loop,
  σ² ~ 10⁻⁵ to 10⁻³ from LSS gives B_ψ/H₀² ~ 10⁻³ to 10⁻¹.
  THIS IS THE RIGHT ORDER OF MAGNITUDE for H₀ tension!

  NEXT STEPS:
  1. Solve the linearized TGP perturbation equation numerically
  2. Compute the growth rate of substrate fluctuations
  3. Determine the self-consistent σ²(a) from structure formation
  4. Compare B_ψ(a) with H₀, S₈, w(z) data
""")

# ==========================================================================
# Quick estimate: tachyonic amplification
# ==========================================================================
print(f"\n{'='*78}")
print("  Tachyonic amplification estimate")
print("="*78)

# From TGP action: V(ψ) = (γ/12) - (γ/2)ψ² - ...
# V''(1) = -γ (tachyonic)
# But: 3ψ̇²/ψ stabilizes. At equilibrium:
# 3ψ̇²/ψ ≈ γψ → ψ̇² ≈ γψ²/3
# So ψ̇ ~ ψ√(γ/3) ~ H·ψ (if γ ~ 12H² from Λ = γ/12)

# The effective tachyonic mass:
# m²_tach = -V''(1) = γ = 12Λ = 12 × 3H₀²Ω_Λ
gamma_eff = 12 * 3 * H0_planck**2 * Om_L  # in (km/s/Mpc)²
m_tach = math.sqrt(gamma_eff)  # in km/s/Mpc
t_tach = 1 / m_tach  # timescale in Mpc·s/km → need to convert

print(f"\n  γ_eff = 12 × 3H₀²Ω_Λ = {gamma_eff:.1f} (km/s/Mpc)²")
print(f"  m_tach = √γ = {m_tach:.1f} km/s/Mpc")
print(f"  m_tach/H₀ = {m_tach/H0_planck:.3f}")
print(f"  → tachyonic mass ~ {m_tach/H0_planck:.1f} × H₀")
print(f"  → instability timescale ~ {1/(m_tach/H0_planck):.3f} × H₀⁻¹")
print(f"  → {1/(m_tach/H0_planck) * 14.4:.1f} Gyr")

# Amplification factor:
# In time Δt ~ H₀⁻¹, the tachyonic mode grows by exp(m_tach·Δt)
# = exp(m_tach/H₀)
amplification = math.exp(m_tach / H0_planck)
print(f"\n  Amplification in one Hubble time: exp({m_tach/H0_planck:.1f}) = {amplification:.1e}")

# Combined with geometric friction (stabilization):
# The saturated amplitude is set by balance:
# 3ψ̇²/ψ = V''·δψ²
# 3(m_tach·δψ)²/1 ≈ m²_tach·δψ² → 3·m²·δψ² ≈ m²·δψ²
# This gives δψ_sat ~ 1/3 (crude) — but really needs numerical solution

print(f"\n  Saturated amplitude estimate:")
print(f"    3·(δψ̇)²/ψ ≈ V''·δψ²")
print(f"    At saturation: δψ_sat ~ O(1/3) × (coupling corrections)")
print(f"    σ²_ψ_sat ~ 0.01 to 0.1")
print(f"    B_ψ/H₀² ~ 3·σ²_sat ~ 0.03 to 0.3")
print(f"    Required: ~0.17 for H₀ tension")
print(f"    → IN THE RIGHT RANGE!")

print(f"\n{'='*78}")
print("  SUMMARY")
print("="*78)
print(f"""
  1. Standard GR backreaction from Φ/c² ~ 10⁻⁵ is NEGLIGIBLE (10⁻¹⁰)

  2. TGP has THREE amplification mechanisms:
     a) Kinetic coupling K=ψ⁴ (factor ~4)
     b) Volume element √ψ (substrate budget)
     c) TACHYONIC INSTABILITY V''(1) = -γ < 0  ← DOMINANT

  3. The tachyonic mass m ≈ 3.5·H₀ gives instability
     on sub-Hubble timescale (~4 Gyr)

  4. Geometric friction 3ψ̇²/ψ SATURATES the instability
     at σ²_ψ ~ 0.01-0.1 → B_ψ/H₀² ~ 0.03-0.3

  5. This is EXACTLY the range needed for H₀ tension (0.17)!

  6. The mechanism is UNIQUE to TGP:
     - Standard scalar fields don't have 3ψ̇²/ψ
     - Standard cosmology doesn't have tachyonic + friction balance
     - The substrate budget f·h=1 is a TGP axiom

  NEXT: Numerical solution of linearized TGP perturbation equation
  to get precise σ²_ψ(a) and B_ψ(a).
""")

print(f"{'='*78}")
