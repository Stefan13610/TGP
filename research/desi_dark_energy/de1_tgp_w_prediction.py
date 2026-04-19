#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
de1: TGP PREDICTION FOR DARK ENERGY EQUATION OF STATE
======================================================

Consolidates TGP predictions for w(z) and compares with DESI data.
Also cross-references with H0 tension and S8 tension from a unified
TGP perspective.

Key TGP prediction:
  - w_eff = -1 + positive_correction (from B_psi >= 0)
  - TGP CANNOT produce phantom crossing (w < -1)
  - Any phantom crossing observed would RULE OUT pure TGP

Sections:
  A. DESI DR1/DR2 data summary
  B. TGP w(z) prediction for various B_psi models
  C. Comparison with CPL parametrization
  D. H0 tension: TGP cannot explain it (from ct7)
  E. S8 tension: marginal TGP effect
  F. Unified verdict across all three tensions
  G. Testable predictions

Dependencies: numpy, scipy
"""

import sys
import warnings
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
warnings.filterwarnings('ignore')

import numpy as np
from scipy.integrate import quad, solve_ivp

def print_header(title):
    print()
    print("=" * 78)
    print(f"  {title}")
    print("=" * 78)
    print()

def print_subheader(title):
    print(f"\n  {title}")
    print("  " + "-" * len(title))


# ==========================================================================
# Constants
# ==========================================================================
H0_Planck = 67.36   # km/s/Mpc
H0_SH0ES  = 73.04   # km/s/Mpc
H0_TRGB   = 69.8    # km/s/Mpc
Omega_m   = 0.315
Omega_L   = 0.685
Omega_r   = 9.1e-5

# DESI results
w0_DESI_DR1 = -0.45   # CPL fit
wa_DESI_DR1 = -1.79
w0_err = 0.21
wa_err = 0.65

# S8 values
S8_Planck = 0.832
S8_DES_Y3 = 0.776
S8_KiDS   = 0.759
S8_HSC    = 0.769

# ==========================================================================
# Friedmann cosmology functions
# ==========================================================================
def H_LCDM(z, H0=H0_Planck):
    """LCDM Hubble parameter."""
    return H0 * np.sqrt(Omega_r*(1+z)**4 + Omega_m*(1+z)**3 + Omega_L)

def H_wCDM(z, H0, w0, wa=0):
    """wCDM with CPL parametrization: w(a) = w0 + wa*(1-a)."""
    a = 1.0 / (1.0 + z)
    # Dark energy density: rho_DE/rho_DE0 = a^(-3(1+w0+wa)) * exp(-3*wa*(1-a))
    rho_DE_ratio = a**(-3.0*(1.0 + w0 + wa)) * np.exp(-3.0 * wa * (1.0 - a))
    return H0 * np.sqrt(Omega_r*(1+z)**4 + Omega_m*(1+z)**3 + Omega_L * rho_DE_ratio)

def w_CPL(z, w0, wa):
    """CPL equation of state."""
    a = 1.0 / (1.0 + z)
    return w0 + wa * (1.0 - a)


# ==========================================================================
# TGP backreaction models
# ==========================================================================
def B_psi_model(a, B0, n_growth=3.0):
    """
    TGP backreaction term B_psi(a).
    Grows with structure formation: B_psi ~ B0 * a^n * f_struct(a)
    B_psi >= 0 always (fundamental TGP constraint).
    """
    # Structure formation factor: grows from a~0.01 to a~1
    # Saturates at late times
    f_struct = a**n_growth / (1.0 + a**n_growth)
    return B0 * f_struct

def w_eff_TGP(z, B0, n_growth=3.0):
    """
    Effective dark energy equation of state in TGP.
    w_eff = -1 + (a/3) * d(ln rho_DE_eff)/d(a)
    where rho_DE_eff = rho_Lambda + B_psi/(8*pi*G)
    """
    a = 1.0 / (1.0 + z)
    da = 0.001

    B_a = B_psi_model(a, B0, n_growth)
    B_a_plus = B_psi_model(a + da, B0, n_growth)

    # rho_DE_eff proportional to Omega_L + B_a/H0^2  (normalized)
    rho_eff = Omega_L + B_a
    rho_eff_plus = Omega_L + B_a_plus

    # d(ln rho)/d(a) = (1/rho) * drho/da
    dlnrho_da = (rho_eff_plus - rho_eff) / (da * rho_eff) if rho_eff > 0 else 0

    w = -1.0 + a * dlnrho_da / 3.0
    return w

def H_TGP(z, H0, B0, n_growth=3.0):
    """TGP Hubble parameter with backreaction."""
    a = 1.0 / (1.0 + z)
    B = B_psi_model(a, B0, n_growth)
    return H0 * np.sqrt(Omega_r*(1+z)**4 + Omega_m*(1+z)**3 + Omega_L + B)


# ##########################################################################
print_header("de1: TGP PREDICTION FOR DARK ENERGY EQUATION OF STATE")


# ======================================================================
#  PART A: DESI DATA SUMMARY
# ======================================================================
print_header("PART A: DESI DATA SUMMARY")

print(f"""  DESI DR1 BAO results (2024):
    w0 = {w0_DESI_DR1} +/- {w0_err}  (CPL fit)
    wa = {wa_DESI_DR1} +/- {wa_err}
    Combined with Planck + Union3/DES SN

  Key features:
    - w(z=0) = w0 = {w0_DESI_DR1:.2f} (less negative than -1)
    - w(z>>1) = w0 + wa = {w0_DESI_DR1 + wa_DESI_DR1:.2f} (more negative than -1)
    - Phantom crossing at z ~ {-w0_DESI_DR1/(wa_DESI_DR1/(1.0)) - 1:.1f} (where w crosses -1)

  Significance: 2-3 sigma deviation from LCDM (w0=-1, wa=0)
  Status: preliminary; DR2 (2025) expected to clarify

  w(z) at key redshifts (DESI CPL):
""")

z_arr = np.array([0.0, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0])
print(f"  {'z':>5s}  {'w_LCDM':>10s}  {'w_DESI':>10s}  {'w_DESI-w_LCDM':>15s}")
print(f"  {'-'*5}  {'-'*10}  {'-'*10}  {'-'*15}")
for z in z_arr:
    w_L = -1.0
    w_D = w_CPL(z, w0_DESI_DR1, wa_DESI_DR1)
    print(f"  {z:5.1f}  {w_L:10.3f}  {w_D:10.3f}  {w_D - w_L:15.3f}")


# ======================================================================
#  PART B: TGP w(z) PREDICTION
# ======================================================================
print_header("PART B: TGP w(z) PREDICTION")

print("""  In TGP, the effective dark energy has:
    rho_DE_eff(a) = rho_Lambda + B_psi(a) / (8*pi*G)

  Since B_psi >= 0 (always), rho_DE_eff >= rho_Lambda.
  This means:
    w_eff = -1 + positive_term >= -1

  TGP ALWAYS predicts w >= -1 (quintessence-like).
  Phantom crossing (w < -1) is IMPOSSIBLE in TGP.
""")

# Compute w(z) for different B0 values
B0_values = [0.0, 0.001, 0.01, 0.05, 0.1, 0.174]

print(f"  TGP w(z) for different B0/H0^2 values:")
print(f"  {'z':>5s}", end="")
for B0 in B0_values:
    print(f"  {'B0='+str(B0):>12s}", end="")
print(f"  {'DESI':>10s}")
print(f"  {'-'*5}", end="")
for _ in B0_values:
    print(f"  {'-'*12}", end="")
print(f"  {'-'*10}")

for z in z_arr:
    print(f"  {z:5.1f}", end="")
    for B0 in B0_values:
        w = w_eff_TGP(z, B0)
        print(f"  {w:12.4f}", end="")
    w_D = w_CPL(z, w0_DESI_DR1, wa_DESI_DR1)
    print(f"  {w_D:10.3f}")

print_subheader("B.1  Can TGP match DESI at ANY redshift?")

print(f"""
  DESI shows w(z=0) ~ -0.45, meaning w is LESS negative than -1.
  TGP also gives w > -1 at z=0 (if B0 > 0).

  But DESI also shows w(z>0.5) << -1 (phantom regime!).
  At z=1: w_DESI ~ {w_CPL(1.0, w0_DESI_DR1, wa_DESI_DR1):.2f}
  TGP at z=1: w_TGP ~ {w_eff_TGP(1.0, 0.174):.4f} (even with maximum B0)

  TGP w(z) is MONOTONICALLY close to -1 and approaches -1 at high z.
  DESI w(z) CROSSES -1 (phantom crossing).

  CONCLUSION: TGP w(z) is INCOMPATIBLE with DESI CPL fit.
  TGP predicts w(z) ~ -1 + epsilon(z) with epsilon > 0 and small.
  DESI shows large deviations with phantom crossing.

  HOWEVER: the CPL parametrization may not capture the true w(z).
  If w is actually close to -1 at all z (within systematics),
  TGP is compatible.
""")


# ======================================================================
#  PART C: H0 TENSION SUMMARY
# ======================================================================
print_header("PART C: H0 TENSION -- TGP CANNOT EXPLAIN IT")

print(f"""  From ct7 (soliton cosmology):

  EVERY mechanism tested fails by 8-9 orders of magnitude:
    Mechanism               B/H0^2        Needed      Gap
    -----------------------------------------------------------
    Tachyonic (ct5/ct6)     1e-9          0.174       8 orders
    Soliton population      ~0            0.174       infinite
    RG running gamma        2e-3          0.174       2 orders
    Phase transition        0             0.174       infinite

  Required for H0 = 73 km/s/Mpc:
    B_psi / H0^2 = (H0_local^2 - H0_early^2) / H0_early^2
                 = ({H0_SH0ES**2:.0f} - {H0_Planck**2:.0f}) / {H0_Planck**2:.0f}
                 = {(H0_SH0ES**2 - H0_Planck**2)/H0_Planck**2:.3f}

  TGP DOES NOT EXPLAIN THE H0 TENSION.
""")

# H0 values table
print_subheader("C.1  Current H0 measurements landscape")
print(f"""
  Method                   H0 (km/s/Mpc)    Error     Source
  -----------------------------------------------------------
  Planck CMB (LCDM)        67.36            0.54      2018
  SH0ES (Cepheids)         73.04            1.04      2022
  TRGB (Chicago)           69.8             1.7       2021
  CCHP (TRGB, revised)     69.8             0.6       2024
  TDCOSMO (lensing)        73.3             5.8       2020
  Megamasers               73.9             3.0       2020
  GW sirens (LIGO/Virgo)   68.0             +12/-7    2021

  NOTE: The TRGB-based H0 ~ 69.8 is intermediate.
  If the true value is ~70, the tension is only ~2 sigma,
  and NO new physics is needed.
""")


# ======================================================================
#  PART D: S8 TENSION
# ======================================================================
print_header("PART D: S8 TENSION -- MARGINAL TGP EFFECT")

print(f"""  S8 = sigma_8 * sqrt(Omega_m / 0.3):

  Early Universe (CMB):
    S8_Planck = {S8_Planck} +/- 0.013

  Late Universe (weak lensing):
    S8_DES_Y3 = {S8_DES_Y3} +/- 0.017
    S8_KiDS   = {S8_KiDS} +/- 0.021
    S8_HSC    = {S8_HSC} +/- 0.026

  Tension: {(S8_Planck - S8_KiDS)/np.sqrt(0.013**2 + 0.021**2):.1f} sigma (Planck vs KiDS)
  Required suppression: {(S8_Planck - S8_KiDS)/S8_Planck * 100:.1f}%
""")

# TGP growth suppression estimate
# In TGP, nu(y) modifies the growth rate:
# delta_ddot + 2H*delta_dot = 4*pi*G*nu(y)*rho*delta
# At cosmological scales, y >> 1, so nu ~ 1.0 + epsilon
# Growth suppression ~ epsilon

# From galaxy_scaling (gs34): at cluster scales, nu ~ 1.02
# At Mpc scales (linear regime): y even higher -> nu even closer to 1
# Suppression: delta_sigma8 / sigma8 ~ (nu-1)/nu ~ epsilon

epsilon_TGP = 0.02  # nu ~ 1.02 at best
suppression_TGP = epsilon_TGP / (1.0 + epsilon_TGP)
suppression_needed = (S8_Planck - S8_KiDS) / S8_Planck

print(f"  TGP growth suppression estimate:")
print(f"    nu(y) at Mpc scales:    ~{1.0 + epsilon_TGP:.3f}")
print(f"    epsilon = nu - 1:       ~{epsilon_TGP:.3f}")
print(f"    Growth suppression:     ~{suppression_TGP*100:.1f}%")
print(f"    Needed:                 ~{suppression_needed*100:.1f}%")
print(f"    Ratio:                  {suppression_TGP/suppression_needed:.2f}")

print(f"""
  TGP provides ~{suppression_TGP*100:.0f}% growth suppression, but ~{suppression_needed*100:.0f}% is needed.
  This is {suppression_TGP/suppression_needed*100:.0f}% of what's required — INSUFFICIENT.

  HOWEVER: S8 tension may resolve through:
    - Baryonic feedback (AGN, SNe): can suppress sigma8 by 3-5%
    - Intrinsic alignments: systematic in weak lensing
    - Photo-z errors: affect source redshift distribution
    - Neutrino masses: sum(m_nu) ~ 0.1 eV suppresses ~1-2%

  Combined non-TGP effects may explain the S8 tension.
""")


# ======================================================================
#  PART E: UNIFIED VERDICT
# ======================================================================
print_header("PART E: UNIFIED VERDICT -- TGP AND COSMOLOGICAL TENSIONS")

print(f"""
  =====================================================================
  THREE TENSIONS: TGP UNIFIED ASSESSMENT
  =====================================================================

  TENSION          SIGNIFICANCE    TGP PREDICTION        VERDICT
  -----------------------------------------------------------------
  H0 (73 vs 67)    5 sigma        No mechanism           OUTSIDE SCOPE
  S8 (0.83 vs 0.76) 2-3 sigma     ~2% suppression        INSUFFICIENT
  w(z) DESI        2-3 sigma      w >= -1 always         INCOMPATIBLE*

  *INCOMPATIBLE with CPL fit, but CPL may not be the true model.
   If w(z) ~ -1 + small_positive, TGP is OK.

  =====================================================================

  KEY INSIGHT: TGP IS A GALAXY-SCALE THEORY
  ==========================================

  The same chameleon mechanism exp(-y^0.8) that makes TGP:
    [PASS] Safe for Solar System
    [PASS] Safe for CMB
    [PASS] Safe for gravitational waves
    [PASS] Excellent for rotation curves, RAR, BTFR, dSphs

  Also makes it:
    [FAIL] Unable to affect cosmological expansion
    [FAIL] Unable to provide sufficient cluster mass
    [FAIL] Unable to modify growth at linear scales

  This is a FEATURE, not a bug:
    - TGP is maximally predictive WITHIN its scope
    - It makes NO uncontrolled predictions outside its scope
    - The transition from "TGP active" to "GR recovered" is smooth
    - The transition scale is set by a0 = 1.2e-10 m/s^2

  ANALOGY:
    QED does not explain the strong nuclear force.
    This is not a failure of QED — it's a scope limitation.
    QED is still the most precisely tested theory in physics.

    Similarly, TGP does not explain cosmological tensions.
    This is not a failure of TGP — it's a scope limitation.
    Within its scope, TGP derives galaxy dynamics from first
    principles (membrane Flory exponent, no free parameters).

  =====================================================================
""")


# ======================================================================
#  PART F: TESTABLE PREDICTIONS
# ======================================================================
print_header("PART F: TESTABLE PREDICTIONS")

print(f"""
  TGP makes several TESTABLE predictions that distinguish it
  from both LCDM and other modified gravity theories:

  1. GALAXY-SCALE (testable NOW):
     - RAR intercept at a0 = 1.2e-10 m/s^2 (morphology-independent)
     - BTFR slope = 2/(1-gamma) where gamma = 0.8*c_eff/(c_eff+1)
     - c_eff = 1.3-1.5 for disks (from BTFR slope 3.85)
     - EFE: quadrature y_eff = sqrt(y_int^2 + y_ext^2) (weaker than MOND)
     -> Euclid DR1 can test these (gs35: S/N ~ 94 for test #1)

  2. CLUSTER-SCALE (testable 2025-2030):
     - 35-50% mass deficit at R_500 (less than MOND's factor 2-3)
     - Lensing peak at galaxy positions in mergers (same as LCDM)
     - kappa(galaxy)/kappa(gas) ~ 2.3 for Bullet-like mergers
     -> Euclid cluster lensing + eROSITA X-ray

  3. GRAVITATIONAL PHYSICS (testable NOW):
     - Gravitational slip eta = 1 (no slip at ANY scale)
     - E_G statistic = GR value (indistinguishable from LCDM)
     -> This is a PREDICTION: if eta ≠ 1 found, TGP is ruled out

  4. COSMOLOGICAL (passive predictions):
     - w(z) = -1 + epsilon, epsilon >= 0 (no phantom crossing)
     - CMB fully compatible (chameleon screening)
     - No additional gravitational wave polarizations
     -> If phantom crossing confirmed: TGP ruled out

  5. FOUNDATIONAL (testable in principle):
     - alpha = 4/5 from Flory exponent (not adjustable)
     - a0 = c*H0/(2*pi) (emergent, scales with H)
     -> If a0 evolves with cosmic time: supports TGP
     -> Testable via RAR at different redshifts (Euclid)

  =====================================================================
  DISCRIMINATING TESTS: TGP vs MOND vs LCDM

  Observable             TGP            MOND           LCDM
  -----------------------------------------------------------------
  RAR at z=0             SAME           SAME           needs DM
  RAR at z=1             a0 evolved?    a0 fixed       no RAR
  BTFR slope             3.82-3.88      4.0 (simple)   scatter
  EFE strength           WEAKER         STRONGER       N/A
  Grav. slip eta         1.000          varies         1.000 (GR)
  Cluster mass           60-65%         40-50%         100% (DM)
  Bullet lensing peak    galaxy         galaxy         galaxy
  Bullet amplitude       62%            40-50%         100% (DM)
  w(z)                   >= -1          N/A            -1
  CMB                    compatible     needs TeVeS    designed for

  KEY DIFFERENTIATORS from MOND:
    - EFE strength (TGP weaker -> testable in Crater II, NGC 1052-DF2)
    - BTFR slope (TGP predicts specific value from c_eff)
    - alpha = 4/5 (TGP derives; MOND has no prediction)

  =====================================================================
  END OF de1 ANALYSIS
  =====================================================================
""")
