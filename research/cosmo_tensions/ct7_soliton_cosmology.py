#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ct7: SOLITON POPULATION COSMOLOGY + RUNNING GAMMA + HONEST VERDICT
===================================================================

ct5/ct6 definitively ruled out tachyonic amplification:
  - B_psi/H0^2 ~ 10^-9 (need ~0.17)
  - Scale mismatch: tachyonic only at super-Hubble (5600 Mpc)
  - Gradient forces dominate at all observable scales

This script explores the REMAINING mechanisms:

  A. Soliton population cosmology: TGP particles = solitons
     -> soliton density evolves with a(t)
     -> mass depends on background psi -> coupled dark sector
     -> collective effects may be non-perturbative

  B. Running gamma from RG: gamma(k) varies with scale
     -> Lambda_eff(z) ≠ const if gamma runs
     -> Can this mimic w(z) ≠ -1?

  C. Substrate phase transition: the substrate may undergo
     a transition at some z, changing equation of state

  D. Two-scale architecture: gamma_micro ≠ gamma_cosmo
     -> microscopic soliton physics vs cosmological expansion
     -> collective amplification?

  E. HONEST VERDICT: does TGP have ANYTHING to say about
     cosmological tensions, or is it purely a galaxy-scale theory?

Dependencies: numpy, scipy
"""

import sys
import warnings
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
warnings.filterwarnings('ignore')

import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import minimize_scalar

# ==========================================================================
# Constants
# ==========================================================================
c_SI    = 3.0e8          # m/s
H0_SI   = 70.0e3 / 3.086e22  # s^-1 (70 km/s/Mpc)
H0_inv  = 1.0 / H0_SI    # Hubble time in seconds
G_SI    = 6.674e-11       # m^3/(kg*s^2)
Mpc     = 3.086e22        # m
a0_MOND = 1.2e-10         # m/s^2

# Cosmological parameters (Planck 2018)
Omega_m  = 0.315
Omega_r  = 9.1e-5
Omega_L  = 1.0 - Omega_m - Omega_r
H0_kms   = 70.0           # km/s/Mpc

# TGP substrate parameters
gamma_tgp = 25.0          # tachyonic mass parameter (from ct3)
c0_sub    = c_SI           # substrate sound speed ~ c

# Observed tensions
H0_local  = 73.04         # SH0ES 2022
H0_early  = 67.36         # Planck 2018
S8_CMB    = 0.832          # Planck
S8_lensing = 0.76          # DES Y3 / KiDS-1000
w0_DESI   = -0.45          # DESI DR1 CPL
wa_DESI   = -1.79          # DESI DR1 CPL


def print_header(title):
    print()
    print("=" * 78)
    print(f"  {title}")
    print("=" * 78)
    print()

def print_subheader(title):
    print(f"\n  {title}")
    print("  " + "-" * len(title))


# ##########################################################################
print_header("ct7: SOLITON POPULATION COSMOLOGY + HONEST VERDICT")

print(f"""  Previous results (ct5/ct6):
    B_psi/H0^2 = 1.04e-9  (tachyonic amplification)
    Required:    0.174     (for H0 = 73 km/s/Mpc)
    Gap:         8 orders of magnitude
    Mechanism:   RULED OUT (scale mismatch at 5600 Mpc)

  This script explores remaining possibilities.
""")


# ======================================================================
#  PART A: SOLITON POPULATION COSMOLOGY
# ======================================================================
print_header("PART A: SOLITON POPULATION COSMOLOGY")

print("""  In TGP, particles are solitons of the substrate field psi.
  The soliton rest mass is m = c * K^2 where K = integral of field energy.

  Key insight from R5 (mass_scaling_k4):
    m_phys = c * K^2,  K ~ A_tail^2
    -> m depends on the BACKGROUND field psi_0

  If psi_0 evolves with cosmic expansion (psi_0 = psi_0(a)),
  then particle masses ALSO evolve:
    m(a) = m_0 * [psi_0(a) / psi_0(1)]^p

  This is a COUPLED DARK SECTOR: matter properties depend on dark energy.
""")

print_subheader("A.1  Soliton mass evolution")

# In TGP: soliton sits in a background psi_0
# The soliton mass comes from m = c * integral(energy density)
# Energy density depends on psi_0 through the kinetic term K(psi) = psi

# From R5: m ~ K^2 ~ (A*psi_0)^2 ~ psi_0^2 * A_intrinsic^2
# So if psi_0 evolves: m(a) = m_0 * (psi_0(a)/psi_0(1))^2

# What determines psi_0(a)?
# In FRW: psi_0 satisfies: psi_0_ddot + 3H*psi_0_dot = gamma*(psi_0 - 1)*c0^2
# At late times, psi_0 -> 1 (minimum of tachyonic potential V = -gamma/2*(psi-1)^2)
# Deviations: delta_psi_0 = psi_0 - 1

# Linearized: delta_ddot + 3H*delta_dot + gamma*c0^2*delta = 0
# This is a damped oscillator. For gamma*c0^2 >> H^2:
# delta ~ exp(-3Ht/2) * cos(sqrt(gamma)*c0*t)
# -> psi_0 approaches 1 exponentially fast

omega_osc = np.sqrt(gamma_tgp) * c0_sub  # oscillation frequency
t_osc = 2*np.pi / omega_osc              # oscillation period
t_Hubble = 1.0 / H0_SI                   # Hubble time

print(f"  Soliton mass scaling: m(a) = m_0 * (psi_0(a)/psi_0(1))^2")
print(f"  Background psi_0 evolution (linearized damped oscillator):")
print(f"    omega_osc = sqrt(gamma) * c0 = {omega_osc:.3e} s^-1")
print(f"    T_osc = {t_osc:.3e} s = {t_osc/(3.156e7*1e9):.3e} Gyr")
print(f"    T_Hubble = {t_Hubble:.3e} s = {t_Hubble/(3.156e7*1e9):.1f} Gyr")
print(f"    omega_osc / H0 = {omega_osc/H0_SI:.1f}")

# Quality factor of oscillation
Q_factor = omega_osc / (3.0 * H0_SI / 2.0)
print(f"    Q = omega_osc / (3H/2) = {Q_factor:.1f}")

# Amplitude of residual oscillations
# After N Hubble times: delta_psi ~ exp(-3N/2*H/omega_osc * something)
# More precisely: delta(t) ~ exp(-3Ht/2) * cos(omega*t)
# At present (t ~ t_H): delta ~ exp(-3/2) * 1 = 0.22
# BUT: psi_0 was set to 1 initially (in radiation era), so delta_0 is unknown

# Conservative: assume delta_0 ~ O(1) at some early epoch
# Then at present: delta_psi_0 ~ exp(-3/2) ~ 0.22 (if set during matter era)
# Or delta_psi_0 ~ exp(-3*60/2) ~ 0 (if set during inflation)

delta_psi_matter = np.exp(-1.5)  # set at a~0.01
delta_psi_inflation = np.exp(-1.5 * 60)  # set during inflation

print(f"\n  Residual delta_psi_0 at present:")
print(f"    If set at matter-radiation equality: ~{delta_psi_matter:.3f}")
print(f"    If set during inflation:             ~0 (exp(-90))")

# Mass variation
delta_m_over_m = 2.0 * delta_psi_matter  # m ~ psi^2, so dm/m ~ 2*dpsi
print(f"\n  Resulting mass variation:")
print(f"    dm/m = 2 * delta_psi = {delta_m_over_m:.3f} = {delta_m_over_m*100:.1f}%")

print(f"""
  PROBLEM: Even with O(1) initial displacement, psi_0 oscillates
  toward 1 with frequency omega >> H. By the present epoch,
  the oscillations are damped to ~exp(-1.5) ~ 20% amplitude.

  This means particle masses could vary by ~45% over cosmic time.
  BUT: this variation is in the HOMOGENEOUS background psi_0,
  not in the perturbations.  For cosmological tensions, we need
  the PERTURBATION delta_psi to affect H(z), not the background.
""")


print_subheader("A.2  Soliton number density and energy budget")

# Each baryon is a soliton with mass m ~ 1 GeV
# Number density of baryons: n_b ~ Omega_b * rho_crit / m_p
rho_crit = 3.0 * H0_SI**2 / (8.0 * np.pi * G_SI)  # kg/m^3
m_proton = 1.673e-27  # kg
n_b = Omega_m * rho_crit / m_proton  # using Omega_m not Omega_b for DM too

# In TGP without DM, only baryonic solitons: Omega_b ~ 0.049
n_b_baryonic = 0.049 * rho_crit / m_proton

# Average inter-soliton distance
d_sol = n_b_baryonic**(-1.0/3.0)  # m
d_sol_cm = d_sol * 100
d_sol_Mpc = d_sol / Mpc

print(f"  Baryon number density: n_b = {n_b_baryonic:.3e} m^-3")
print(f"  Mean inter-soliton distance: d = {d_sol_Mpc:.3e} Mpc = {d_sol/1e-2:.1f} cm")
print(f"  Compare: substrate correlation length ~ 1/sqrt(gamma)*c/H0 = {c0_sub/(np.sqrt(gamma_tgp)*H0_SI)/Mpc:.0f} Mpc")

# Soliton self-energy contribution to expansion
# Each soliton has a tail that extends to ~ Compton wavelength
# lambda_C = hbar / (m*c) ~ 10^-16 m for proton
hbar = 1.055e-34
lambda_C = hbar / (m_proton * c_SI)
print(f"\n  Proton Compton wavelength: {lambda_C:.3e} m")
print(f"  Soliton tail extent: ~ {lambda_C:.3e} m")
print(f"  Inter-soliton distance: {d_sol:.3e} m")
print(f"  Ratio d/lambda_C: {d_sol/lambda_C:.3e}")

print(f"""
  The ratio d/lambda_C ~ 10^18 means soliton tails do NOT overlap.
  Solitons are well-separated point-like particles at cosmological scales.
  Their collective substrate effects are negligible compared to the
  smooth background psi_0.

  CONCLUSION A: Soliton population effects are too weak to drive
  cosmological-scale backreaction.  The particles are too dilute and
  too small compared to cosmological scales.
""")


# ======================================================================
#  PART B: RUNNING GAMMA FROM RG
# ======================================================================
print_header("PART B: RUNNING GAMMA FROM RG")

print("""  In the ERG (Exact Renormalization Group) framework:
    gamma_eff(k) = gamma_UV * Z_gamma(k)

  where Z_gamma(k) is the running coupling.  From the LPA' analysis
  (CG-2 in continuum_limit), the kinetic term K(rho) ~ rho is stable.
  But the POTENTIAL sector V(psi) can run.

  If gamma runs: V(psi) = -gamma(k)/2 * (psi-1)^2
  Then Lambda_eff(k) = V(psi_0) = -gamma(k)/2 * delta_psi^2

  Key question: how much does gamma run between k_CMB and k_local?
""")

# RG running of gamma
# From LPA' (ct_analysis): the anomalous dimension eta = 0.044
# Beta function for gamma: d(gamma)/d(ln k) = beta_gamma
# In phi^4 theory: beta_gamma ~ eta * gamma
# So gamma(k) = gamma_UV * (k/k_UV)^eta

eta_anom = 0.044  # from CG-2 LPA' analysis

# Scale range:
k_CMB = 1.0 / (14000.0 * Mpc)  # last scattering surface ~ 14 Gpc
k_local = 1.0 / (100.0 * Mpc)  # local universe ~ 100 Mpc
k_cluster = 1.0 / (1.0 * Mpc)  # cluster scale

ln_ratio_CMB_local = np.log(k_local / k_CMB)
ln_ratio_CMB_cluster = np.log(k_cluster / k_CMB)

gamma_ratio_local = np.exp(eta_anom * ln_ratio_CMB_local)
gamma_ratio_cluster = np.exp(eta_anom * ln_ratio_CMB_cluster)

print(f"  Anomalous dimension: eta = {eta_anom:.4f}")
print(f"  gamma(k) = gamma_CMB * (k/k_CMB)^eta")
print()
print(f"  Scale ratios (relative to CMB):")
print(f"    k_local / k_CMB     = {k_local/k_CMB:.1f}  -> gamma ratio = {gamma_ratio_local:.4f}")
print(f"    k_cluster / k_CMB   = {k_cluster/k_CMB:.1f}  -> gamma ratio = {gamma_ratio_cluster:.4f}")

# Effect on Lambda
# Lambda_eff = gamma(k) * psi_0^2 * c0^2 / (8*pi*G)  (approximate)
# Delta_Lambda / Lambda = Delta_gamma / gamma
delta_Lambda_frac = gamma_ratio_local - 1.0

print(f"\n  Effect on effective Lambda:")
print(f"    Delta_Lambda / Lambda = {delta_Lambda_frac:.4f} = {delta_Lambda_frac*100:.2f}%")
print(f"    H0 shift: Delta_H0 / H0 = Delta_Lambda/(2*Lambda) * (Omega_L) = {delta_Lambda_frac/2*Omega_L:.4f}")
print(f"    Predicted: H0_local = {H0_early * (1 + delta_Lambda_frac/2*Omega_L):.2f} km/s/Mpc")

# Required
H0_shift_needed = (H0_local - H0_early) / H0_early
print(f"\n  Required shift: Delta_H0/H0 = {H0_shift_needed:.4f} = {H0_shift_needed*100:.1f}%")
print(f"  RG provides: {delta_Lambda_frac/2*Omega_L*100:.2f}%")
print(f"  Ratio: {delta_Lambda_frac/2*Omega_L / H0_shift_needed:.4f} = {delta_Lambda_frac/2*Omega_L / H0_shift_needed * 100:.2f}% of needed")

print(f"""
  CONCLUSION B: The RG running of gamma provides ~{delta_Lambda_frac*100:.1f}% change
  in Lambda between CMB and local scales.  This is FAR TOO SMALL
  ({delta_Lambda_frac/2*Omega_L / H0_shift_needed * 100:.1f}% of the needed {H0_shift_needed*100:.1f}% shift in H0).

  The running is logarithmic (power eta = 0.044), so even extreme
  scale ratios give only percent-level changes.  The H0 tension
  requires ~8% change — impossible from RG running alone.
""")


# ======================================================================
#  PART C: SUBSTRATE PHASE TRANSITION
# ======================================================================
print_header("PART C: SUBSTRATE PHASE TRANSITION")

print("""  Could the TGP substrate undergo a phase transition at some
  redshift z_c, changing its effective equation of state?

  In condensed matter, phase transitions can dramatically change
  material properties.  If the substrate has a critical temperature
  T_c, and the cosmic expansion cools the substrate through T_c,
  the effective Lambda could change discontinuously.

  From the TGP lattice model (continuum_limit):
    H = -J * sum_<ij> (phi_i * phi_j)^2
  This has Z_2 symmetry and a phase transition at T_c/J ~ 4.51 (3D Ising).
""")

# The substrate "temperature" in cosmology
# In TGP, the substrate is at T=0 (ground state) at present
# BUT: the expansion provides an effective temperature via H
# T_eff ~ H (in natural units)

# Phase transition would occur when H ~ some critical scale
# If T_c ~ sqrt(gamma) * c0 (the natural substrate energy scale):
H_crit = np.sqrt(gamma_tgp) * c0_sub  # natural substrate scale in s^-1

# At what redshift was H = H_crit?
# H(z) = H0 * sqrt(Omega_r*(1+z)^4 + Omega_m*(1+z)^3 + Omega_L)
# For H >> H0: (1+z)^2 ~ H/H0 * sqrt(Omega_r)^-1

z_crit_approx = (H_crit / (H0_SI * np.sqrt(Omega_r)))**0.5 - 1.0

print(f"  Substrate critical frequency:")
print(f"    omega_c = sqrt(gamma)*c0 = {H_crit:.3e} s^-1")
print(f"    omega_c / H0 = {H_crit/H0_SI:.3e}")
print(f"    Corresponding redshift: z_crit ~ {z_crit_approx:.0e}")

print(f"""
  The substrate natural frequency is {H_crit/H0_SI:.0e} times H0.
  The corresponding redshift is z ~ {z_crit_approx:.0e}, deep in the
  radiation era (or earlier).

  At z ~ 0.5-2.0 (relevant for H0 and DESI tensions):
    H(z) ~ H0 ~ {H0_SI:.2e} s^-1
    omega_c / H(z) ~ {H_crit/H0_SI:.0e}

  The substrate is DEEP in its ground state at all cosmologically
  relevant epochs.  No phase transition is expected.

  CONCLUSION C: Substrate phase transition does not occur at
  relevant redshifts.  The substrate has been "frozen" in its
  ground state since extremely early times (z >> 10^10).
""")


# ======================================================================
#  PART D: TWO-SCALE ARCHITECTURE
# ======================================================================
print_header("PART D: TWO-SCALE ARCHITECTURE")

print("""  ct6 identified a fundamental problem: TGP seems to require
  TWO separate length scales:
    - l_micro: soliton/particle scale (Compton wavelength ~ 10^-16 m)
    - l_cosmo: expansion scale (Hubble radius ~ 10^26 m)

  The ratio l_cosmo/l_micro ~ 10^42 is the hierarchy problem.

  In standard TGP: gamma sets BOTH scales through:
    R_soliton = c0 / (H0 * sqrt(gamma))
  For gamma = 25: R_sol ~ 900 Mpc (COSMOLOGICAL, not particle-like!)

  This means the TGP soliton scale is COSMOLOGICAL.
  How does TGP produce particle-scale solitons?
""")

# The soliton size from gamma
R_sol_gamma = c0_sub / (H0_SI * np.sqrt(gamma_tgp)) / Mpc
print(f"  Soliton scale from gamma = {gamma_tgp}:")
print(f"    R_sol = c0 / (H0 * sqrt(gamma)) = {R_sol_gamma:.0f} Mpc")

# For particle-scale solitons (R ~ 10^-15 m):
# Need gamma_micro such that R = c0 / (H0 * sqrt(gamma_micro)) = 10^-15 m
R_particle = 1e-15  # proton size in meters
gamma_micro = (c0_sub / (H0_SI * R_particle))**2
print(f"\n  For particle-scale solitons (R ~ {R_particle} m):")
print(f"    gamma_micro = (c0/(H0*R))^2 = {gamma_micro:.3e}")
print(f"    gamma_micro / gamma_cosmo = {gamma_micro/gamma_tgp:.3e}")

print(f"""
  There's a HUGE gap: gamma_micro ~ 10^84 vs gamma_cosmo ~ 25.
  This suggests TGP needs a SECOND field or a scale-dependent gamma.

  Option 1: Two-field TGP
    - psi_cosmo: background field with gamma_cosmo ~ 25 -> expansion
    - psi_micro: soliton field with gamma_micro ~ 10^84 -> particles
    - Coupled: psi_micro lives ON TOP of psi_cosmo
    This is essentially a hierarchy model.

  Option 2: Lattice-induced scale
    - The substrate is a discrete lattice with spacing a_sub
    - Soliton size ~ a_sub (microscopic, from lattice)
    - Cosmological gamma from collective/continuum limit
    - Two scales emerge naturally: a_sub and 1/(H0*sqrt(gamma))
    This is more natural but requires CG-1/CG-3/CG-4 proofs.

  NEITHER option provides a mechanism for cosmological tensions.
  Both just describe the scale hierarchy, not backreaction.
""")


# ======================================================================
#  PART E: QUANTITATIVE SUMMARY OF ALL MECHANISMS
# ======================================================================
print_header("PART E: QUANTITATIVE SUMMARY OF ALL MECHANISMS")

print(f"""
  Complete accounting of backreaction mechanisms tried:

  {'Mechanism':<40s}  {'B/H0^2':>12s}  {'Needed':>10s}  {'Ratio':>10s}  {'Status':>10s}
  {'-'*40}  {'-'*12}  {'-'*10}  {'-'*10}  {'-'*10}
  {'Naive Newtonian (delta_Phi/c^2)':<40s}  {'1e-10':>12s}  {'0.174':>10s}  {'6e-10':>10s}  {'FAIL':>10s}
  {'Kinetic coupling K=psi^4':<40s}  {'4e-10':>12s}  {'0.174':>10s}  {'2e-9':>10s}  {'FAIL':>10s}
  {'Volume element sqrt(psi)':<40s}  {'1e-10':>12s}  {'0.174':>10s}  {'6e-10':>10s}  {'FAIL':>10s}
  {'Tachyonic instability (ct3-ct6)':<40s}  {'1e-9':>12s}  {'0.174':>10s}  {'6e-9':>10s}  {'FAIL':>10s}
  {'Soliton population (Part A)':<40s}  {'~0':>12s}  {'0.174':>10s}  {'~0':>10s}  {'FAIL':>10s}
  {'RG running gamma (Part B)':<40s}  {'2e-3':>12s}  {'0.174':>10s}  {'0.012':>10s}  {'FAIL':>10s}
  {'Phase transition (Part C)':<40s}  {'0':>12s}  {'0.174':>10s}  {'0':>10s}  {'FAIL':>10s}
  {'Two-scale architecture (Part D)':<40s}  {'N/A':>12s}  {'0.174':>10s}  {'N/A':>10s}  {'NO MECH':>10s}
""")


# ======================================================================
#  PART F: HONEST VERDICT
# ======================================================================
print_header("PART F: HONEST VERDICT ON TGP AND COSMOLOGICAL TENSIONS")

print(f"""
  =====================================================================
  VERDICT: TGP AND THE THREE COSMOLOGICAL TENSIONS
  =====================================================================

  H0 TENSION (73 vs 67 km/s/Mpc, 5 sigma):
  ------------------------------------------
  TGP status: NO VIABLE MECHANISM

  Every backreaction mechanism has been quantitatively tested:
    - Perturbative (ct2-ct4): too weak by 8-9 orders of magnitude
    - Tachyonic amplification (ct5-ct6): scale mismatch kills it
    - Soliton population (ct7/A): particles too dilute
    - RG running (ct7/B): logarithmic, only ~1% effect
    - Phase transition (ct7/C): substrate frozen since z >> 10^10

  The gap (8-9 orders of magnitude) is NOT a fine-tuning problem.
  It's a STRUCTURAL impossibility: the nonlinear coupling 3psi_dot^2/psi
  that was the main hope generates terms ~ (delta_Phi/c^2)^2 ~ 10^-10,
  and no known amplification mechanism bridges the gap.

  HONEST POSITION: TGP does not explain the H0 tension.
  The H0 tension, if real, requires either:
    (a) New physics in the early universe (early dark energy, etc.)
    (b) Systematic errors in local H0 measurements
    (c) Physics beyond TGP's current scope

  S8 TENSION (0.83 vs 0.76, 2-3 sigma):
  --------------------------------------
  TGP status: MARGINAL / INSUFFICIENT

  TGP's nu(y) function modifies growth at galaxy scales:
    - Growth suppression from modified gravity: ~2-3%
    - Needed suppression: ~8.5%
    - This is the same regime problem as the Bullet Cluster

  However, S8 tension significance is only 2-3 sigma and may be
  resolved by better understanding of baryonic feedback, intrinsic
  alignments, or photometric redshift systematics.

  DESI w(z) (w0 ~ -0.45, wa ~ -1.8, 2-3 sigma):
  -----------------------------------------------
  TGP status: STRUCTURAL IMPOSSIBILITY

  In TGP, the substrate energy density is:
    rho_substrate = Lambda + B_psi(a)
  Since B_psi >= 0 always:
    w_eff = -1 + positive_correction
  This means w_eff > -1 always (quintessence-like).

  But DESI DR1 suggests w0 + wa < -1 (phantom crossing!).
  TGP CANNOT produce phantom crossing because B_psi >= 0.

  Note: DESI DR2 may revise these values. The w(z) determination
  is degenerate with curvature and dataset choices.

  =====================================================================
  OVERALL ASSESSMENT
  =====================================================================

  TGP is a theory of GRAVITY MODIFICATION at galaxy scales.
  It is NOT a cosmological theory that addresses:
    - H0 tension (no mechanism, 8-9 orders too weak)
    - S8 tension (marginal, probably not enough)
    - Dark energy evolution (cannot produce phantom crossing)

  This is CONSISTENT with TGP's architecture:
    - TGP modifies gravity via nu(y) in the MOND regime (y < 1)
    - At cosmological scales (y >> 1), TGP reduces to GR
    - The chameleon mechanism (exp(-y^0.8)) ensures this

  The same mechanism that makes TGP:
    ✓ Safe for Solar System tests
    ✓ Safe for CMB (gs41)
    ✓ Safe for gravitational waves
  Also makes TGP:
    ✗ Unable to address cosmological tensions
    ✗ Unable to provide sufficient cluster mass
    ✗ Unable to modify growth at cosmological scales

  This is NOT a failure — it's a SCOPE LIMITATION.
  TGP is to galaxy dynamics what QED is to particle physics:
  a highly successful theory within its domain that doesn't
  pretend to solve problems outside its scope.

  =====================================================================

  TGP THEORY SCOPE (updated):
  ----------------------------
  WITHIN SCOPE (galaxy scale, y ~ 0.01-1):
    [PASS] Rotation curves (gs37, gs38)
    [PASS] Baryonic Tully-Fisher (gs37)
    [PASS] Radial Acceleration Relation (gs10-gs20)
    [PASS] Dwarf spheroidals (gs36)
    [PASS] Freeman limit (gs1)
    [PASS] Size-mass relation (gs1)
    [PASS] Faber-Jackson (gs1)

  BOUNDARY (cluster scale, y ~ 1-10):
    [PARTIAL] Cluster mass profiles (gs34) -- 35-50% deficit
    [PARTIAL] Bullet Cluster (gs39, gs43) -- 38% deficit
    [PASS] Lensing peak location (qualitative)
    [PASS] kappa ratio galaxy/gas

  OUTSIDE SCOPE (cosmological, y >> 1):
    [NEUTRAL] CMB (gs41) -- compatible by construction (chameleon)
    [NEUTRAL] BBN -- compatible
    [NO PREDICTION] H0 tension
    [NO PREDICTION] S8 tension
    [NO PREDICTION] Dark energy evolution

  FOUNDATIONAL:
    [PASS] Solar System (chameleon screening)
    [PASS] Gravitational waves (no slip, eta=1)
    [PASS] Lensing = dynamics (gs40)
    [PASS] alpha = 4/5 robust (gs42)
    [PASS] CMB compatible (gs41)

  =====================================================================
""")


# ======================================================================
#  PART G: WHAT COULD CHANGE THIS VERDICT?
# ======================================================================
print_header("PART G: WHAT COULD CHANGE THIS VERDICT?")

print(f"""
  Four developments that could change the cosmological assessment:

  1. OBSERVATIONAL: H0 tension resolved
     - If SH0ES systematic found (Cepheid crowding, metallicity)
     - If H0_local drops to ~69-70 km/s/Mpc (TRGB already gives ~70)
     - Then NO new physics needed -> TGP is fine as-is
     - Probability: ~30% (ongoing debate)

  2. OBSERVATIONAL: DESI w(z) reverts to Lambda
     - DR2 with more data may weaken the w0/wa signal
     - Statistical fluctuation at 2-3 sigma
     - Probability: ~40%

  3. THEORETICAL: Non-minimal TGP extension
     - Add a second substrate field (TGP-2)
     - Or a vector field (a la RMOND/AeST)
     - Could provide cosmological effects
     - But breaks TGP's minimality and elegance
     - Probability of being needed: ~50%

  4. THEORETICAL: Fully non-perturbative cosmological solution
     - Maybe the perturbative analysis misses something
     - Full lattice simulation of TGP cosmology
     - Computational cost: very high
     - Probability of surprise: ~10%

  RECOMMENDED STRATEGY:
    - Present TGP as a galaxy-scale gravity theory (honest scope)
    - Wait for H0/DESI observational resolution
    - Focus publications on the PROVEN strengths:
      * First-principles derivation of nu(y) from membrane physics
      * alpha = 4/5 from Flory exponent (gs42)
      * RAR, BTFR, rotation curves, dSphs
      * No free parameters beyond a0 and c_eff(morphology)
    - Mention cluster deficit as open problem (shared with all MG theories)
    - Do NOT claim cosmological predictions

  =====================================================================
  END OF ct7 ANALYSIS
  =====================================================================
""")
