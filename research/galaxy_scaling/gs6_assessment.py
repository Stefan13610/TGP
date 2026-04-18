#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gs6_assessment.py: Honest assessment of TGP soliton tail mechanism.

RESULTS FROM gs4-gs5:
  gs4 (Monte Carlo): v² enhancement 1.5-3× → extrapolates to massive at N~10¹¹
  gs5 (convolution): alpha_TGP = -0.026 (nearly flat!) vs alpha_N = -0.368

  BUT: Both results are MISLEADING due to artifacts:

  1. gs4 SELF-TERM PROBLEM:
     |∇g|² = Σ|∇δᵢ|² + cross-terms
     Self-terms are part of measured particle mass → can't count as extra
     Cross-terms are 20-50% destructive → reduce effective gravity

  2. gs5 OSCILLATION PROBLEM:
     v_TGP oscillates between positive and zero (δ changes sign)
     α = -0.026 fits the ENVELOPE, not the physical curve
     Real galaxy (N~10¹¹): random phases average out → no coherent tail signal

THIS SCRIPT: Quantify the ACTUAL tail contribution after proper averaging.
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
import math

print("="*78)
print("  HONEST ASSESSMENT: TGP soliton tail mechanism for rotation curves")
print("="*78)

# ==========================================================================
# 1. WHY THE TAIL AVERAGES OUT
# ==========================================================================
print(f"\n{'='*78}")
print(f"  1. Why oscillating tails average out")
print(f"{'='*78}")

A_tail = 0.066  # tail envelope |δ|·r from gs5
M_sol = 0.118   # "Newtonian mass" of soliton (from asymptotic M_enc average)
lambda_tail = 2 * np.pi  # oscillation period in soliton units

print(f"""
  Single soliton tail: delta(r) ~ {A_tail:.3f} * sin(r) / r
  Oscillation period: lambda = 2*pi = {lambda_tail:.2f} soliton units

  For N solitons at random positions in a galaxy of radius R_gal:

  A) DETERMINISTIC (coherent) part:
     delta_coherent(r) = integral n(r') * <delta(|r-r'|)>_angle * 4pi r'^2 dr'

     Angular average: <sin(|r-r'|)/|r-r'|>_angle = 2*sin(r)*sin(r')/(r*r')

     So: delta_coh(r) = (8pi*A*sin(r)/r) * integral n(r')*sin(r')*r' dr'

     For smooth n(r') over many periods (R >> lambda):
       integral n(r')*sin(r')*r' dr'  OSCILLATES → small residual

     Estimate (stationary phase):
       |integral| ~ sqrt(pi/(2*k)) * n(R) * R  where k=1 (in soliton units)
       ~ n(R) * R

     For Hernquist galaxy: n(R) ~ N*a/(2pi*R*(R+a)^3) ~ N/(2pi*R^3) at R>>a
       n(R)*R ~ N/(2pi*R^2)

     Ratio to Newtonian:
       delta_coh / delta_Newton ~ (A * n(R)*R) / (N*M/R)
                                = A*R / (2pi*R^2 * M) = A/(2pi*R*M)
""")

# Compute for various galaxy sizes
print(f"  Galaxy (in soliton units):")
print(f"  {'R_gal':>6s} {'N':>8s} {'delta_coh/delta_N':>17s} {'Comment':>20s}")
print(f"  {'---'*20}")

for R_gal, N in [(10, 200), (30, 1000), (100, 10000), (1000, 1e6), (1e6, 1e11)]:
    ratio = A_tail / (2*np.pi*R_gal*M_sol)
    if ratio > 0.1:
        comment = "significant"
    elif ratio > 0.01:
        comment = "small but detectable"
    elif ratio > 1e-4:
        comment = "negligible"
    else:
        comment = "VANISHING"
    print(f"  {R_gal:6.0f} {N:8.0e} {ratio:17.2e} {comment:>20s}")

# ==========================================================================
# 2. RANDOM FLUCTUATIONS
# ==========================================================================
print(f"\n{'='*78}")
print(f"  2. Random tail fluctuations (incoherent part)")
print(f"{'='*78}")

print(f"""
  B) STOCHASTIC (incoherent) part:
     Each soliton contributes random phase tail at distance r.
     For N solitons: delta_random ~ sqrt(N) * A / r
     (central limit theorem for oscillating functions)

     Ratio to Newtonian:
       delta_random / delta_N ~ sqrt(N)*A/r / (N*M/r) = A/(sqrt(N)*M)
""")

print(f"  {'N':>10s} {'delta_rand/delta_N':>17s}")
print(f"  {'---'*12}")
for N in [1e2, 1e4, 1e6, 1e8, 1e11]:
    ratio = A_tail / (np.sqrt(N) * M_sol)
    print(f"  {N:10.0e} {ratio:17.2e}")

print(f"\n  For N=10^11 (real galaxy): random tail contribution = {A_tail/(np.sqrt(1e11)*M_sol):.1e}")
print(f"  → COMPLETELY negligible")

# ==========================================================================
# 3. SELF-TERM REANALYSIS
# ==========================================================================
print(f"\n{'='*78}")
print(f"  3. Self-terms: are they truly 'extra' gravity?")
print(f"{'='*78}")

print(f"""
  The gs4 Monte Carlo measured |nabla g|^2 / g as "extra mass".
  But: the soliton solution ALREADY includes the nonlinear g'^2/g term.
  The soliton's measured gravitational mass (from Kepler orbits, etc.)
  already reflects the full TGP solution, including all self-energy.

  Therefore:
  - Self-terms |nabla delta_i|^2: ALREADY in measured mass → NO extra effect
  - Cross-terms nabla delta_i . nabla delta_j: TWO-BODY interaction correction

  The gs4 result M_NL/M_enc ~ 0.29-1.24 was dominated by self-terms.
  The CROSS-terms were -20% to +50% of self-terms, and their sign FLUCTUATED
  between realizations → they are NOISE, not signal.

  VERDICT: The self-term mechanism of gs4 does NOT provide extra gravity.
  The cross-term mechanism is real but SMALL and INCONSISTENT in sign.
""")

# ==========================================================================
# 4. WHAT ABOUT THE g'^2/g MODIFICATION TO POISSON?
# ==========================================================================
print(f"\n{'='*78}")
print(f"  4. Does g'^2/g modify Poisson equation at galactic scales?")
print(f"{'='*78}")

print(f"""
  The TGP field equation for a SINGLE source (galaxy-scale potential Phi):
    nabla^2 delta + (nabla delta)^2/(1+delta) + delta = source

  In the weak field (delta << 1):
    nabla^2 delta + delta = source - (nabla delta)^2 + O(delta^2)

  The nonlinear correction -(nabla delta)^2:
    |nabla delta| ~ delta / r ~ GM/(r^2 c^2)
    (nabla delta)^2 ~ (GM)^2 / (r^4 c^4)

  Compare with source (Poisson) ~ delta/r ~ GM/(r^2 c^2):
    ratio = (nabla delta)^2 / source ~ GM/(r c^2)

  At galaxy scale (MW: M~10^11 M_sun, r~10 kpc):
    GM/(r c^2) ~ 4.3e-6 * 10^11 / (10 * 9e10)
              ~ 4.3e5 / 9e11 ~ 5e-7

  → Nonlinear correction is 5x10^-7 of Newtonian signal.
  → COMPLETELY NEGLIGIBLE at galaxy scales.
  → Would only matter at Phi/c^2 ~ 1 (neutron stars, black holes).
""")

G_astro = 4.3e-6  # (km/s)^2 kpc / M_sun
c_kms = 3e5  # km/s

for name, M, R in [("Dwarf (DDO 154)", 5e8, 2), ("Milky Way", 7e10, 10),
                    ("M87", 1e12, 50), ("IC 1101", 5e12, 100)]:
    ratio = G_astro * M / (R * c_kms**2)
    print(f"  {name:<20s}: GM/(Rc^2) = {ratio:.1e}")

# ==========================================================================
# 5. THE FUNDAMENTAL ISSUE
# ==========================================================================
print(f"\n{'='*78}")
print(f"  5. THE FUNDAMENTAL ISSUE")
print(f"{'='*78}")

print(f"""
  TGP soliton equation: g'' + g'^2/g + 2g'/r + g = 1

  The nonlinear term g'^2/g is a relativistic correction.
  For weak gravitational fields (galaxies): g ~ 1 + 10^-6 → correction ~ 10^-12
  For rotation curves: need correction of order 1 (100% of Newtonian)

  GAP: 12 orders of magnitude between TGP nonlinear correction
       and the effect needed for flat rotation curves.

  This is THE SAME problem as in ct5/ct6 (cosmological tensions):
  the TGP nonlinear effects are FUNDAMENTALLY TOO WEAK at astrophysical scales.

  The soliton tail oscillation sin(r)/r:
  - IS real and long-range (1/r envelope)
  - DOES modify the far-field potential compared to 1/r
  - BUT averages out for large N (random phases)
  - AND the coherent remnant is tiny compared to Newtonian

  CONCLUSION: The TGP soliton tail mechanism, in its current formulation,
  CANNOT explain flat rotation curves or replace dark matter.
""")

# ==========================================================================
# 6. WHAT COULD WORK?
# ==========================================================================
print(f"\n{'='*78}")
print(f"  6. Alternative mechanisms that COULD work")
print(f"{'='*78}")

print(f"""
  The phenomenological model (gs1/gs2) WORKS beautifully:
  - BTFR with exponent 4 ✓
  - Freeman limit 137 ≈ 140 M_sun/pc^2 ✓
  - a0 ≈ c*H0/(2pi) ✓
  - Correct phantom DM density ✓

  SOMETHING produces these scaling relations. If not the soliton tail,
  what mechanism in TGP could give a0?

  OPTION A: SUBSTRATE BOUNDARY CONDITION
    The TGP substrate is a finite expanding medium.
    At the causal boundary (c/H0), boundary conditions affect
    the propagation of gravity → modified force law at large r.
    This is NOT a perturbative effect — it's a global property
    of the substrate geometry.

  OPTION B: ENTROPIC/THERMODYNAMIC GRAVITY (Verlinde-like)
    If gravity in TGP emerges from substrate entropy:
    S = integral sqrt(g) d^3x (substrate volume = entropy measure)
    → Entropic force F = T * dS/dr could have different r-dependence
    → The "temperature" T ~ H0 (de Sitter horizon temperature)
    → Gives: a0 ~ c * H0 (matching observation!)

  OPTION C: MODIFIED DISPERSION RELATION
    TGP substrate has discrete structure at some scale l.
    Gravity propagation has modified dispersion:
    omega^2 = c^2 k^2 * (1 - l^2 k^2 + ...)
    At long wavelengths (galaxy scales): extra terms modify potential.
    If l ~ c/H0 (Hubble scale): the modification kicks in at
    accelerations ~ c*H0 ~ a0.

  OPTION D: TWO-SCALE TGP (revisited from ct6)
    The soliton equation has TWO natural scales:
    1. Microscopic: l_micro (particle size)
    2. Cosmological: l_cosmo = c/H0 (expansion horizon)
    The RATIO l_cosmo/l_micro ~ 10^60 (hierarchy problem)
    A collective effect coupling both scales could produce a0
    at the GEOMETRIC MEAN: sqrt(l_micro * l_cosmo).

    If l_micro = l_Planck = 1.6e-35 m:
    sqrt(l_Planck * c/H0) = sqrt(1.6e-35 * 1.3e26) = sqrt(2.1e-9)
                           = 4.6e-5 m

    a_geometric = c^2 / sqrt(l_Planck * c/H0) = 9e16 / 4.6e-5 = 2e21
    → Too large!

    But if we use the de Sitter acceleration:
    a_dS = c * sqrt(Lambda/3) = c * H0 * sqrt(Omega_Lambda)
         = 3e8 * 2.2e-18 * 0.83 = 5.4e-10 m/s^2

    And a0 ≈ a_dS / (2*pi) ≈ 8.6e-11 m/s^2 — close to 1.2e-10!

  OPTION E: THE SOLITON EQUATION IS NOT THE RIGHT ONE
    Perhaps the galactic-scale TGP equation is DIFFERENT from
    the particle-scale soliton equation.
    The "+g = 1" term might have a different coefficient at large scales:
      g'' + g'^2/g + 2g'/r + mu^2(g - 1) = source
    where mu^2 = 1 at microscopic scale (soliton equation)
    but mu^2 = (a0*r_gal/c^2)^2 ~ 10^-12 at galactic scale.
    → The tail wavelength lambda = 2*pi/mu ~ 10^6 * r_gal
    → THIS could modify gravity at galactic scales!

  ASSESSMENT:
  Options A, B, D are conceptual and need formalization.
  Option C is similar to emergent MOND approaches.
  Option E is the most concrete: a scale-dependent mu in the TGP equation
  could produce the observed phenomenology.

  RECOMMENDED NEXT STEP: Explore Option E —
  derive the effective mu^2 at galactic scales from TGP cosmological context.
""")

# ==========================================================================
# 7. SCORE CARD
# ==========================================================================
print(f"\n{'='*78}")
print(f"  7. SCORE CARD — what works, what doesn't")
print(f"{'='*78}")

print(f"""
  ┌───────────────────────────────────────────────────────────────────┐
  │  MECHANISM                    │ WORKS? │ GAP      │ NEXT STEP    │
  ├───────────────────────────────┼────────┼──────────┼──────────────┤
  │ Phenomenological (a0)          │  ✅    │ None     │ Derive a0    │
  │ Soliton tail sin(r)/r         │  ❌    │ Avg→0    │ Abandon      │
  │ MC self-terms |∇δ|²           │  ❌    │ In mass  │ Abandon      │
  │ MC cross-terms ∇δᵢ·∇δⱼ       │  ❌    │ Small    │ Abandon      │
  │ Weak-field nonlinearity g'²/g │  ❌    │ 10^12    │ Abandon      │
  │ Tachyonic backreaction (ct5)  │  ❌    │ 10^9     │ Abandon      │
  ├───────────────────────────────┼────────┼──────────┼──────────────┤
  │ Substrate boundary cond.      │  ❓    │ Unknown  │ EXPLORE      │
  │ Entropic/Verlinde gravity     │  ❓    │ Unknown  │ EXPLORE      │
  │ Modified dispersion           │  ❓    │ Unknown  │ EXPLORE      │
  │ Scale-dependent mu(r)         │  ❓    │ Unknown  │ EXPLORE ←    │
  └───────────────────────────────┴────────┴──────────┴──────────────┘

  BOTTOM LINE:
  The PHENOMENOLOGY is solid (a0 = cH0/2pi gives everything).
  But no MICRO-MECHANISM in TGP has been found to produce it.
  The most promising direction: scale-dependent TGP equation.
""")
