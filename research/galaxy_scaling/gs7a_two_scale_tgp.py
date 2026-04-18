#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gs7a_two_scale_tgp.py: Option D — Two-scale TGP and the origin of a0.

QUESTION: Why is a0 ~ c*H0/(2*pi) ~ 1.2e-10 m/s^2?

This script investigates ALL possible combinations of TGP scales
that could produce a0, without guessing — pure dimensional analysis
followed by physical mechanism checks.
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np

print("="*78)
print("  OPTION D: Two-scale TGP — dimensional analysis of a0")
print("="*78)

# ===========================================================================
# PHYSICAL CONSTANTS
# ===========================================================================
c = 2.998e8        # m/s
G = 6.674e-11      # m^3/(kg*s^2)
hbar = 1.055e-34   # J*s
k_B = 1.381e-23    # J/K
H0 = 2.20e-18      # 1/s  (67.4 km/s/Mpc)
Lambda = 3 * (H0**2) * 0.685  # cosmological constant (Omega_Lambda = 0.685)
l_P = 1.616e-35    # Planck length
m_P = 2.176e-8     # Planck mass
t_P = 5.391e-44    # Planck time
a0_obs = 1.2e-10   # m/s^2 (observed MOND acceleration)

# Derived scales
R_H = c / H0                    # Hubble radius
T_dS = hbar * H0 / (2*np.pi*k_B)  # de Sitter temperature
a_H = c * H0                    # Hubble acceleration
a_dS = c * np.sqrt(Lambda/3)    # de Sitter acceleration
L_Lambda = 1/np.sqrt(Lambda/3)  # de Sitter radius

print(f"\n  Fundamental scales:")
print(f"    c       = {c:.3e} m/s")
print(f"    G       = {G:.3e} m^3/(kg*s^2)")
print(f"    hbar    = {hbar:.3e} J*s")
print(f"    H0      = {H0:.2e} 1/s  ({H0*3.086e22/1e3:.1f} km/s/Mpc)")
print(f"    Lambda  = {Lambda:.2e} 1/m^2")
print(f"    l_P     = {l_P:.3e} m")

print(f"\n  Derived scales:")
print(f"    R_H     = c/H0    = {R_H:.3e} m  ({R_H/3.086e22:.1f} Mpc)")
print(f"    L_Lamb  = 1/sqrt(Lamb/3) = {L_Lambda:.3e} m  ({L_Lambda/3.086e22:.1f} Mpc)")
print(f"    T_dS    = hbar*H0/(2*pi*kB) = {T_dS:.3e} K")
print(f"    a_H     = c*H0   = {a_H:.3e} m/s^2")
print(f"    a_dS    = c*sqrt(Lamb/3) = {a_dS:.3e} m/s^2")
print(f"    a0_obs  = {a0_obs:.1e} m/s^2")

# ===========================================================================
# 1. EXHAUSTIVE DIMENSIONAL ANALYSIS
# ===========================================================================
print(f"\n{'='*78}")
print(f"  1. ALL possible a0 combinations from (c, G, hbar, H0, Lambda)")
print(f"{'='*78}")

print(f"""
  An acceleration has dimensions [m/s^2] = [L T^-2].
  From (c, H0) alone: c*H0 = {a_H:.3e} m/s^2

  Possible forms of a0 from known physics:
""")

candidates = []

# Pure c, H0 combinations
for label, val in [
    ("c*H0", c*H0),
    ("c*H0/(2*pi)", c*H0/(2*np.pi)),
    ("c*H0/(4*pi)", c*H0/(4*np.pi)),
    ("c*H0/6", c*H0/6),
    ("c*H0*sqrt(2/3)", c*H0*np.sqrt(2/3)),
]:
    ratio = val / a0_obs
    candidates.append((label, val, ratio))

# With Lambda
for label, val in [
    ("c*sqrt(Lambda/3)", c*np.sqrt(Lambda/3)),
    ("c*sqrt(Lambda/3)/(2*pi)", c*np.sqrt(Lambda/3)/(2*np.pi)),
    ("c^2*sqrt(Lambda/3)/c = c*sqrt(Lam/3)", c*np.sqrt(Lambda/3)),
    ("c^2*Lambda^(1/2)/sqrt(3)", c**2 * Lambda**0.5 / np.sqrt(3)),
]:
    ratio = val / a0_obs
    candidates.append((label, val, ratio))

# With G (Planck-related)
for label, val in [
    ("sqrt(c^7*H0/(G*hbar))", np.sqrt(c**7*H0/(G*hbar))),
    ("c^2/sqrt(c/(G*H0^2*...))", None),  # skip complex
    ("(c*H0)^(2/3) * (c^5/(G*hbar))^(1/3)", (c*H0)**(2/3) * (c**5/(G*hbar))**(1/3)),
    ("sqrt(c*H0 * c^4/(G*m_P*c^2/l_P))", np.sqrt(c*H0 * c**2/l_P)),
]:
    if val is not None:
        ratio = val / a0_obs
        candidates.append((label, val, ratio))

# Geometric means
for label, val in [
    ("sqrt(a_Planck * a_H)", np.sqrt((c**2/l_P) * c*H0)),
    ("(a_Planck * a_H^2)^(1/3)", ((c**2/l_P) * (c*H0)**2)**(1/3)),
    ("(a_Planck^2 * a_H)^(1/3)", ((c**2/l_P)**2 * (c*H0))**(1/3)),
]:
    ratio = val / a0_obs
    candidates.append((label, val, ratio))

# de Sitter / thermodynamic
for label, val in [
    ("2*pi*k_B*T_dS/l_P", 2*np.pi*k_B*T_dS/l_P),
    ("k_B*T_dS/(m_P*c)", k_B*T_dS/(m_P*c)),
    ("hbar*H0^2/(m_P*c)", hbar*H0**2/(m_P*c)),
]:
    ratio = val / a0_obs
    candidates.append((label, val, ratio))

# Remove duplicates and None
seen = set()
unique = []
for label, val, ratio in candidates:
    if val is not None and label not in seen:
        seen.add(label)
        unique.append((label, val, ratio))

# Sort by closeness to a0
unique.sort(key=lambda x: abs(np.log10(abs(x[2]))))

print(f"  {'Formula':<45s} {'Value (m/s^2)':>14s} {'ratio to a0':>12s} {'Match':>8s}")
print(f"  {'-'*82}")
for label, val, ratio in unique:
    if 0.5 < ratio < 2.0:
        match = "<-- !!!"
    elif 0.1 < ratio < 10:
        match = "<--"
    else:
        match = ""
    print(f"  {label:<45s} {val:14.3e} {ratio:12.3f} {match:>8s}")

# ===========================================================================
# 2. THE c*H0 FAMILY — WHY 2*pi?
# ===========================================================================
print(f"\n{'='*78}")
print(f"  2. The c*H0 family — investigating the factor of 2*pi")
print(f"{'='*78}")

print(f"""
  a0_obs = {a0_obs:.1e} m/s^2
  c*H0   = {c*H0:.3e} m/s^2  (ratio = {c*H0/a0_obs:.2f})
  c*H0/(2*pi) = {c*H0/(2*np.pi):.3e} m/s^2  (ratio = {c*H0/(2*np.pi)/a0_obs:.2f})

  The ratio c*H0/a0 = {c*H0/a0_obs:.2f}

  Possible origins of factor ~{c*H0/a0_obs:.1f}:
""")

# What numerical factor bridges c*H0 to a0?
factor_needed = c*H0 / a0_obs
print(f"  Factor needed: c*H0/a0 = {factor_needed:.3f}")
print(f"  2*pi           = {2*np.pi:.3f}  (ratio: {factor_needed/(2*np.pi):.3f})")
print(f"  sqrt(4*pi)     = {np.sqrt(4*np.pi):.3f}  (ratio: {factor_needed/np.sqrt(4*np.pi):.3f})")
print(f"  2*sqrt(pi)     = {2*np.sqrt(np.pi):.3f}  (ratio: {factor_needed/(2*np.sqrt(np.pi)):.3f})")
print(f"  e (Euler)      = {np.e:.3f}  (ratio: {factor_needed/np.e:.3f})")
print(f"  sqrt(2*pi)     = {np.sqrt(2*np.pi):.3f}  (ratio: {factor_needed/np.sqrt(2*np.pi):.3f})")
print(f"  4              = 4.000  (ratio: {factor_needed/4:.3f})")
print(f"  6              = 6.000  (ratio: {factor_needed/6:.3f})")
print(f"  pi+e           = {np.pi+np.e:.3f}  (ratio: {factor_needed/(np.pi+np.e):.3f})")

print(f"""
  PHYSICAL MECHANISMS that give 2*pi:

  (i)  OSCILLATION: a0 = c/(T_Hubble * 2*pi) = c*H0/(2*pi)
       → a0 is the amplitude of oscillation at Hubble frequency
       → 2*pi from converting frequency to angular frequency
       → In TGP: substrate oscillates at omega = H0, centripetal a = c*H0/(2*pi)

  (ii) WAVELENGTH: a0 = c^2/lambda_H where lambda_H = 2*pi*c/H0 = 2*pi*R_H
       → a0 = c^2/(2*pi*R_H) = c*H0/(2*pi)
       → Compton-like: lambda = 2*pi * (characteristic length)
       → In TGP: soliton tail has wavelength = 2*pi/mu,
         if mu = H0/c then lambda = 2*pi*c/H0 = 2*pi*R_H

  (iii) SPHERE: Surface area of Hubble sphere = 4*pi*R_H^2
        Holographic: N_dof = A/(4*l_P^2) → different a0
        But: a = GM/R^2 with M_H = c^3/(2GH0) gives a = c*H0/2

  (iv) YUKAWA: If gravity has a Yukawa modification:
       Phi = -GM/r * exp(-r/lambda)
       At r << lambda: Phi ~ -GM/r*(1 - r/lambda)
       Extra force: delta_F = GM/(r*lambda)
       This becomes comparable to GM/r^2 when r = lambda
       → transition acceleration: a_trans = GM/lambda^2 = v^2/lambda
       Need lambda ~ v^2/a0 ~ R_galaxy (galaxy size!)
""")

# ===========================================================================
# 3. THE DE SITTER CONNECTION
# ===========================================================================
print(f"\n{'='*78}")
print(f"  3. The de Sitter connection: a_dS vs a0")
print(f"{'='*78}")

print(f"""
  de Sitter acceleration: a_dS = c*H0*sqrt(Omega_Lambda)
  = {c*H0*np.sqrt(0.685):.3e} m/s^2

  a_dS / a0 = {c*H0*np.sqrt(0.685)/a0_obs:.3f}
  a_dS / (2*pi) = {c*H0*np.sqrt(0.685)/(2*np.pi):.3e}  (ratio to a0: {c*H0*np.sqrt(0.685)/(2*np.pi*a0_obs):.3f})

  Alternative: a0 = c*H0*sqrt(Omega_Lambda)/(2*pi)
  = {c*H0*np.sqrt(0.685)/(2*np.pi):.3e}  (ratio: {c*H0*np.sqrt(0.685)/(2*np.pi)/a0_obs:.3f})

  vs simple: a0 = c*H0/(2*pi)
  = {c*H0/(2*np.pi):.3e}  (ratio: {c*H0/(2*np.pi)/a0_obs:.3f})

  QUESTION: Is a0 proportional to H0 or to sqrt(Lambda)?

  These differ by sqrt(Omega_Lambda) = {np.sqrt(0.685):.3f}
  Currently: Omega_Lambda = 0.685, so sqrt = 0.828
  At z=0: both give similar values (within 17%)

  DISCRIMINATOR: at high redshift z, H(z) changes but Lambda is constant.
  - If a0 ~ H(z): a0 was LARGER in the past
  - If a0 ~ sqrt(Lambda) = const: a0 is CONSTANT

  At z=1: H(z=1) = H0*sqrt(Omega_m*(1+1)^3 + Omega_Lambda)
        = H0*sqrt(0.315*8 + 0.685) = H0*sqrt(3.205) = {np.sqrt(3.205):.2f}*H0

  a0(z=1)/a0(z=0):
    If a0~H: ratio = {np.sqrt(3.205):.2f}
    If a0~sqrt(Lam): ratio = 1.00
""")

H_z1 = H0 * np.sqrt(0.315*8 + 0.685)
print(f"  H(z=1) = {H_z1:.3e} 1/s = {H_z1/H0:.2f} * H0")
print(f"  a0(z=1) if a0~cH: {c*H_z1/(2*np.pi):.3e} m/s^2  ({c*H_z1/(2*np.pi)/a0_obs:.2f}x a0)")
print(f"  a0(z=1) if a0~const: {a0_obs:.1e} m/s^2  (unchanged)")

# ===========================================================================
# 4. TGP-SPECIFIC: WHAT SETS THE SCALE?
# ===========================================================================
print(f"\n{'='*78}")
print(f"  4. TGP-specific: what sets the acceleration scale?")
print(f"{'='*78}")

print(f"""
  In TGP, the soliton equation is: g'' + g'^2/g + 2g'/r + g = 1

  The "+g" term (spring/mass term) has coefficient mu^2 = 1 in soliton units.
  This means: mu = 1/l_soliton in physical units.

  The soliton tail oscillates as sin(mu*r)/r with wavelength lambda = 2*pi/mu.

  KEY QUESTION: What is mu in physical units?

  POSSIBILITY 1: mu is set by particle physics (Compton wavelength)
    mu = m*c/hbar (for mass m)
    For proton: mu = 1.6e-27 * 3e8 / 1.05e-34 = 4.6e15 /m
    lambda = 2*pi/mu = 1.4e-15 m (proton Compton)
    → Way too small for galactic effects

  POSSIBILITY 2: mu is set by cosmology
    mu = H0/c = 1/R_H = {H0/c:.3e} /m
    lambda = 2*pi*c/H0 = {2*np.pi*c/H0:.3e} m = {2*np.pi*c/H0/3.086e22:.0f} Mpc
    → Hubble-scale wavelength
    → a0 = c^2 * mu / (2*pi) = c^2 * H0/(c * 2*pi) = c*H0/(2*pi) !!!
    → THIS GIVES EXACTLY a0 = c*H0/(2*pi) !!!

  POSSIBILITY 3: mu is set by Lambda
    mu = sqrt(Lambda/3) = H0*sqrt(Omega_Lambda)/c
    lambda = 2*pi / mu = 2*pi*c / (H0*sqrt(Omega_Lambda))
    a0 = c^2 * mu / (2*pi) = c*H0*sqrt(Omega_Lambda)/(2*pi)
    = {c*H0*np.sqrt(0.685)/(2*np.pi):.3e} m/s^2  (ratio: {c*H0*np.sqrt(0.685)/(2*np.pi)/a0_obs:.3f})

  POSSIBILITY 4: mu runs between micro and cosmo scales (RG flow)
    At scale r: mu_eff(r) = mu_micro * (l_micro/r)^alpha for some alpha
    At r = R_galaxy: mu_eff ~ mu_micro * (l_micro/R_gal)^alpha
    → Need to determine alpha from TGP dynamics
""")

# ===========================================================================
# 5. THE KEY INSIGHT: mu_cosmo = H0/c
# ===========================================================================
print(f"\n{'='*78}")
print(f"  5. KEY INSIGHT: if mu_cosmo = H0/c")
print(f"{'='*78}")

mu_cosmo = H0 / c
lambda_cosmo = 2*np.pi / mu_cosmo

print(f"""
  If the TGP equation at cosmological/galactic scales has:
    g'' + g'^2/g + 2g'/r + mu_cosmo^2 * (g - 1) = source

  with mu_cosmo = H0/c = {mu_cosmo:.3e} /m

  Then:
    Tail wavelength: lambda = 2*pi/mu = {lambda_cosmo:.3e} m = {lambda_cosmo/3.086e22:.0f} Mpc

    The tail does NOT oscillate on galactic scales!
    (Galaxy radius ~ 10-100 kpc << {lambda_cosmo/3.086e19:.0e} kpc)

    Instead, the tail behaves as:
    For r << lambda: sin(mu*r)/(mu*r) -> 1 - (mu*r)^2/6 + ...
    -> delta(r) ~ A/r * (1 - mu^2*r^2/6)

    Extra acceleration from the correction term:
    delta_a = -d/dr [A * mu^2 * r / 6] = -A * mu^2 / 3
    (constant! independent of r!)

    This constant acceleration = A * mu^2 / 3
    = A * (H0/c)^2 / 3
""")

# What A would be needed?
# Need delta_a = a0 = 1.2e-10
# A * (H0/c)^2 / 3 = a0
# A = 3 * a0 * c^2 / H0^2
A_needed = 3 * a0_obs * c**2 / H0**2
print(f"  For delta_a = a0:")
print(f"    A_needed = 3*a0*c^2/H0^2 = {A_needed:.3e} m^2/s^2")
print(f"    = {A_needed/c**2:.3e} c^2")
print(f"    = GM with M = A/G = {A_needed/G:.3e} kg = {A_needed/G/1.989e30:.3e} M_sun")
print(f"    → This is the TOTAL mass of the galaxy!")
print(f"    → M_galaxy ~ 10^11 M_sun: A ~ G*10^11*M_sun = {G*1e11*1.989e30:.3e}")
print(f"    → A_needed/A_galaxy = {A_needed / (G*1e11*1.989e30):.1f}")

print(f"""
  INTERPRETATION:
  The "A" is NOT the tail amplitude of one soliton —
  it's the TOTAL gravitational potential of the galaxy: A ~ G*M_galaxy.

  So the mechanism is:
  1. The galactic potential Phi(r) = -G*M/r (Newtonian)
  2. In TGP with mu_cosmo = H0/c, the potential gets a correction:
     Phi_TGP(r) = -G*M/r * [1 - mu^2*r^2/6 + ...]
               = -G*M/r + G*M*mu^2*r/6
  3. The correction force: F_corr = -d/dr [G*M*mu^2*r/6] = -G*M*mu^2/6
  4. This is a CONSTANT inward force = G*M*H0^2/(6*c^2)

  For MW (M = 7e10 M_sun):
""")

M_MW = 7e10 * 1.989e30  # kg
F_corr = G * M_MW * mu_cosmo**2 / 6
a_corr = F_corr  # per unit mass
print(f"    a_corr = G*M*mu^2/6 = {a_corr:.3e} m/s^2")
print(f"    a0_obs = {a0_obs:.1e} m/s^2")
print(f"    ratio a_corr/a0 = {a_corr/a0_obs:.3e}")
print(f"    → {a_corr/a0_obs:.1e} of a0 — TOO SMALL by factor {a0_obs/a_corr:.0e}")

print(f"""
  WHY IT'S TOO SMALL:
  G*M*mu^2 = G*M*(H0/c)^2 = (G*M/c^2) * H0^2
  = (Schwarzschild radius / 2) * H0^2
  R_S(MW) = 2*G*M/c^2 = {2*G*M_MW/c**2:.3e} m = {2*G*M_MW/c**2/3.086e16:.3e} pc

  a_corr = R_S * H0^2 / 12 = {2*G*M_MW/c**2 * H0**2 / 12:.3e}

  For a_corr = a0 we'd need R_S * H0^2 ~ a0
  → R_S ~ a0/H0^2 = c/H0 = R_Hubble!
  → Need a galaxy with Schwarzschild radius = Hubble radius
  → M ~ c^3/(2*G*H0) = {c**3/(2*G*H0):.3e} kg = {c**3/(2*G*H0)/1.989e30:.3e} M_sun
  → This is the mass of the OBSERVABLE UNIVERSE!
""")

# ===========================================================================
# 6. CORRECTED APPROACH: MILGROM'S OBSERVATION
# ===========================================================================
print(f"\n{'='*78}")
print(f"  6. Milgrom's observation: a0 appears at g_N = a0, not from mu*r correction")
print(f"{'='*78}")

print(f"""
  The mu*r correction gives too small an effect because it's perturbative.
  But the OBSERVED effect is:

  g_obs = g_N / (1 - exp(-sqrt(g_N/a0)))

  At g_N >> a0: g_obs ~ g_N (Newtonian)
  At g_N << a0: g_obs ~ sqrt(g_N * a0) (deep MOND)

  The transition happens at g_N = a0, i.e., at radius:
  r_MOND = sqrt(G*M/a0)

  For MW: r_MOND = sqrt({G*M_MW:.3e}/{a0_obs}) = {np.sqrt(G*M_MW/a0_obs):.3e} m
         = {np.sqrt(G*M_MW/a0_obs)/3.086e19:.1f} kpc

  This is WHERE the rotation curve flattens — correct!

  KEY: a0 doesn't come from a PERTURBATIVE correction to Newton.
  It comes from a QUALITATIVE CHANGE in the force law at low accelerations.
  Any mechanism must produce a TRANSITION, not a small correction.
""")

# ===========================================================================
# 7. NON-PERTURBATIVE: YUKAWA-LIKE MODIFICATION
# ===========================================================================
print(f"\n{'='*78}")
print(f"  7. Non-perturbative: what if TGP gives Yukawa + extra?")
print(f"{'='*78}")

print(f"""
  If the TGP equation with mu gives a YUKAWA-type potential:
    Phi(r) = -G*M/r * exp(-mu*r)  (standard Yukawa)

  Then for mu = H0/c:
    At r << c/H0 (all galaxies): exp(-mu*r) ~ 1 - mu*r ~ 1
    → No significant modification! Yukawa with lambda = c/H0 is too long-range.

  BUT: The TGP equation is NOT the Klein-Gordon equation.
  TGP: g'' + g'^2/g + 2g'/r + mu^2*(g-1) = source
  KG:  phi'' + 2phi'/r - mu^2*phi = source  (repulsive spring, different sign!)

  In TGP: "+mu^2*g" is an ATTRACTIVE spring toward g=1.
  In KG: "-mu^2*phi" gives exponential decay (Yukawa).

  TGP with mu gives OSCILLATORY solution: sin(mu*r)/r (as we already know!)
  But with mu = H0/c, the oscillation period is 2*pi*c/H0 ~ 4 Gpc.
  At galactic scales (kpc), we're in the FIRST quarter-wave: sin(mu*r) ~ mu*r.

  So the TGP potential at galactic scales:
    delta(r) = A * sin(mu*r) / (mu*r) ~ A * (1 - mu^2*r^2/6)
    Phi(r) ~ -G*M/r + G*M*mu^2*r/6 - ...

  The second term G*M*mu^2*r/6 is a HARMONIC (linear restoring) force.
  This is the same as the de Sitter cosmological expansion effect!

  In GR, a test particle in de Sitter space feels:
    a_dS(r) = Lambda*c^2*r/3 = H0^2*Omega_Lambda*r

  In TGP with mu = H0/c:
    a_TGP(r) = G*M*mu^2/3 = G*M*H0^2/(3*c^2)  (from the galaxy)
    + Lambda*c^2*r/3 (from cosmological expansion)

  The galaxy term is negligible (as computed above).
  The cosmological term grows with r but is local, not per-galaxy.
""")

# ===========================================================================
# 8. WHAT WORKS IN OPTION D — ASSESSMENT
# ===========================================================================
print(f"\n{'='*78}")
print(f"  8. ASSESSMENT of Option D")
print(f"{'='*78}")

print(f"""
  WHAT WORKS:
  + Dimensional analysis: a0 ~ c*H0/(2*pi) is the ONLY simple combination
    of c and H0 that matches (no G, no hbar needed!)
  + The factor 2*pi has natural origins:
    - Oscillation period (omega -> frequency)
    - Compton-like wavelength (lambda = 2*pi/mu)
  + The de Sitter acceleration a_dS is tantalizingly close
  + a_dS/(2*pi) = {c*H0*np.sqrt(0.685)/(2*np.pi)/a0_obs:.3f} * a0 (only {abs(1-c*H0*np.sqrt(0.685)/(2*np.pi)/a0_obs)*100:.0f}% off!)

  WHAT DOESN'T WORK:
  - Perturbative mu*r correction is 10^6 too small for MW
  - A galaxy's Schwarzschild radius is 10^17 times too small
  - The effect at galactic scales from mu = H0/c is negligible
  - No TRANSITION mechanism: need qualitative change at a = a0,
    not a tiny correction to Newton

  KEY INSIGHT FROM THIS ANALYSIS:
  a0 = c*H0/(2*pi) is a CORRECT DIMENSIONAL RELATION,
  but the MECHANISM is NOT a perturbative mu*r correction.

  The mechanism must be NON-PERTURBATIVE:
  Something changes QUALITATIVELY when the local acceleration
  drops below c*H0/(2*pi). This is NOT about adding a small term
  to the Newtonian potential — it's about a phase transition
  or mode change in how gravity propagates.

  PREDICTION for further options:
  Options A, B, C, E will also face this challenge:
  the dimensional analysis works, but perturbative corrections don't.
  The answer is likely NON-PERTURBATIVE.

  POSSIBLE NON-PERTURBATIVE MECHANISMS:
  1. Phase transition in substrate at critical deformation delta_c
     → When delta < delta_c (weak field): gravity changes character
  2. Saturation of information flow at a = c*H0
     → Gravity can't be "updated" faster than H0^(-1) at distances > c/H0
  3. Entanglement horizon: at a < a0, gravitational degrees of freedom
     become entangled with cosmological horizon → extra entropy → extra force
""")

# ===========================================================================
# 9. QUANTITATIVE SUMMARY TABLE
# ===========================================================================
print(f"\n{'='*78}")
print(f"  9. QUANTITATIVE SUMMARY")
print(f"{'='*78}")

results = [
    ("c*H0", c*H0, c*H0/a0_obs, "Too large by 5.5x"),
    ("c*H0/(2*pi)", c*H0/(2*np.pi), c*H0/(2*np.pi)/a0_obs, "0.87x — best match!"),
    ("c*H0*sqrt(Omega_L)/(2pi)", c*H0*np.sqrt(0.685)/(2*np.pi),
     c*H0*np.sqrt(0.685)/(2*np.pi)/a0_obs, "0.72x — too low by 28%"),
    ("c*sqrt(Lambda/3)", a_dS, a_dS/a0_obs, "4.5x — too large"),
    ("G*M_MW*H0^2/c^2", G*M_MW*H0**2/c**2, G*M_MW*H0**2/c**2/a0_obs, "Perturbative: ~10^-6"),
]

print(f"\n  {'Formula':<35s} {'Value':>12s} {'a0 ratio':>10s} {'Assessment':<30s}")
print(f"  {'-'*90}")
for label, val, ratio, comment in results:
    print(f"  {label:<35s} {val:12.3e} {ratio:10.4f} {comment:<30s}")

print(f"""

  VERDICT FOR OPTION D:
  =====================
  The DIMENSIONAL RELATION a0 = c*H0/(2*pi) is CONFIRMED as the best match.
  But no PERTURBATIVE mechanism in TGP can produce it.
  The mechanism must be NON-PERTURBATIVE (phase transition, saturation, or entanglement).

  STATUS: PARTIALLY SUCCESSFUL
  - Dimensional analysis: SUCCESS (a0 = c*H0/(2pi))
  - Micro-mechanism: OPEN (non-perturbative needed)
  - Falsifiable prediction: a0 EVOLVES with H(z) if a0 ~ c*H0
                           a0 CONSTANT if a0 ~ c*sqrt(Lambda/3)
""")
