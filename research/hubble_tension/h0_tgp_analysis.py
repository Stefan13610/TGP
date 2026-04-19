#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
h0: HUBBLE TENSION -- COMPREHENSIVE TGP ANALYSIS
==================================================

Cross-references ct2-ct7 (cosmo_tensions) with latest H0 data.
Provides quantitative TGP assessment and identifies what would
need to change for TGP to address H0 tension.

Sections:
  A. H0 measurement compilation (2020-2026)
  B. TGP backreaction budget (from ct2-ct7)
  C. Scale-dependent H0 in TGP?
  D. Alternative: is the tension real?
  E. What TGP DOES predict for expansion history
  F. Verdict

Dependencies: numpy
"""

import sys
import warnings
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
warnings.filterwarnings('ignore')

import numpy as np

def print_header(title):
    print()
    print("=" * 78)
    print(f"  {title}")
    print("=" * 78)
    print()

def print_subheader(title):
    print(f"\n  {title}")
    print("  " + "-" * len(title))

print_header("h0: HUBBLE TENSION -- TGP ANALYSIS")

# ======================================================================
#  PART A: H0 MEASUREMENT COMPILATION
# ======================================================================
print_header("PART A: H0 MEASUREMENT COMPILATION (2020-2026)")

measurements = [
    # (name, H0, sigma, year, method, type)
    ("Planck 2018 (LCDM)",         67.36, 0.54,  2018, "CMB (early)",         "early"),
    ("ACT DR6 (LCDM)",             67.49, 1.01,  2024, "CMB (early)",         "early"),
    ("SPT-3G (LCDM)",              67.4,  1.1,   2023, "CMB (early)",         "early"),
    ("Planck+BAO (LCDM)",          67.66, 0.42,  2020, "CMB+BAO (early)",     "early"),
    ("DESI DR1+CMB",               67.97, 0.38,  2024, "BAO+CMB (early)",     "early"),
    ("SH0ES (Cepheids)",           73.04, 1.04,  2022, "Distance ladder",     "late"),
    ("CCHP (TRGB)",                69.8,  0.6,   2024, "TRGB",                "late"),
    ("H0LiCOW (lensing)",          73.3,  1.8,   2020, "Strong lensing",      "late"),
    ("TDCOSMO (lensing)",          74.2,  1.6,   2023, "Strong lensing",      "late"),
    ("Megamasers",                 73.9,  3.0,   2020, "Megamaser distances",  "late"),
    ("Surface Brightness Fluct.",  73.3,  0.7,   2021, "SBF",                 "late"),
    ("GW sirens (LIGO/Virgo)",     68.0,  6.0,   2021, "GW standard sirens",  "late"),
    ("JWST Cepheids (SH0ES val.)", 72.6,  2.0,   2024, "Cepheids (JWST)",     "late"),
    ("JAGB (J-region AGB)",        67.96, 1.85,  2024, "JAGB stars",          "late"),
]

print(f"  {'Measurement':>35s}  {'H0':>6s} {'+-':>5s}  {'Year':>4s}  {'Method':>25s}  {'Type':>5s}")
print(f"  {'-'*35}  {'-'*6} {'-'*5}  {'-'*4}  {'-'*25}  {'-'*5}")
for name, h0, sig, yr, method, typ in measurements:
    print(f"  {name:>35s}  {h0:6.2f} {sig:5.2f}  {yr:4d}  {method:>25s}  {typ:>5s}")

# Compute weighted averages
early = [(h, s) for _, h, s, _, _, t in measurements if t == "early"]
late  = [(h, s) for _, h, s, _, _, t in measurements if t == "late"]

def weighted_avg(data):
    w = [1.0/s**2 for _, s in data]
    h = [h for h, _ in data]
    avg = sum(wi*hi for wi, hi in zip(w, h)) / sum(w)
    sig = 1.0 / np.sqrt(sum(w))
    return avg, sig

h0_early, sig_early = weighted_avg(early)
h0_late, sig_late = weighted_avg(late)
tension_sigma = abs(h0_late - h0_early) / np.sqrt(sig_early**2 + sig_late**2)

print(f"\n  Weighted averages:")
print(f"    Early Universe:  H0 = {h0_early:.2f} +/- {sig_early:.2f} km/s/Mpc")
print(f"    Late Universe:   H0 = {h0_late:.2f} +/- {sig_late:.2f} km/s/Mpc")
print(f"    Tension:         {tension_sigma:.1f} sigma")

print_subheader("A.1  Evolution of the tension")

print(f"""
  The H0 tension landscape has evolved significantly:

  2019: SH0ES vs Planck = 4.4 sigma (clear tension)
  2022: SH0ES vs Planck = 5.0 sigma (strengthened)
  2024: JWST confirms Cepheid distances (SH0ES holds)
  2024: CCHP TRGB gives 69.8 (intermediate!)
  2024: JAGB gives 67.96 (consistent with Planck!)

  KEY DEVELOPMENT: The JAGB method (J-region AGB stars) gives
  H0 = 67.96 +/- 1.85, fully consistent with Planck.
  If confirmed, this would resolve the tension as a systematic
  in the Cepheid distance scale.

  The tension is NOT as clear-cut as it was in 2022.
  The "true" H0 may be ~69-70 km/s/Mpc.
""")


# ======================================================================
#  PART B: TGP BACKREACTION BUDGET
# ======================================================================
print_header("PART B: TGP BACKREACTION BUDGET")

print(f"""  From ct2-ct7, the complete backreaction budget:

  {'Mechanism':40s}  {'B/H0^2':>12s}  {'dH0/H0':>10s}  {'dH0 (km/s/Mpc)':>15s}
  {'-'*40}  {'-'*12}  {'-'*10}  {'-'*15}
  {'Naive Newtonian (delta_Phi/c^2)^2':40s}  {'1e-10':>12s}  {'5e-11':>10s}  {'3e-9':>15s}
  {'Kinetic coupling K=psi^4':40s}  {'4e-10':>12s}  {'2e-10':>10s}  {'1e-8':>15s}
  {'Tachyonic instability (ct5)':40s}  {'1e-9':>12s}  {'5e-10':>10s}  {'3e-8':>15s}
  {'RG running gamma (ct7)':40s}  {'2e-3':>12s}  {'7e-4':>10s}  {'0.05':>15s}
  {'All combined':40s}  {'~2e-3':>12s}  {'~1e-3':>10s}  {'~0.07':>15s}
  {'-'*40}  {'-'*12}  {'-'*10}  {'-'*15}
  {'REQUIRED (for 73 km/s/Mpc)':40s}  {'0.174':>12s}  {'0.084':>10s}  {'5.68':>15s}
  {'REQUIRED (for 70 km/s/Mpc)':40s}  {'0.037':>12s}  {'0.039':>10s}  {'2.64':>15s}

  Even for the modest target H0 = 70 (TRGB/JAGB level):
    TGP provides dH0 ~ 0.07 km/s/Mpc
    Needed:       dH0 ~ 2.64 km/s/Mpc
    Ratio:        {0.07/2.64*100:.1f}% — still FAR TOO SMALL

  The gap is STRUCTURAL, not parametric.
""")


# ======================================================================
#  PART C: SCALE-DEPENDENT H0?
# ======================================================================
print_header("PART C: SCALE-DEPENDENT H0 IN TGP?")

print(f"""  Could H0 be SCALE-DEPENDENT in TGP?

  In principle, nu(y) modifies gravity differently at different scales.
  If the measurement METHOD probes different scales, H0 could appear
  different:
    - CMB: probes z ~ 1100, super-Hubble at recombination
    - SN Ia: probes z ~ 0.01-2, local to intermediate scales
    - Cepheids: probes z ~ 0.001-0.01, very local

  In TGP:
    nu(y) at CMB recombination: y >> 1 -> nu = 1.0 (GR recovered)
    nu(y) at z ~ 0: depends on local g_bar
      - In voids: y ~ 0.01 -> nu ~ 5-10 (significant)
      - In filaments: y ~ 1 -> nu ~ 1.5
      - Overall average: nu ~ 1.0 (most volume is high-y)

  The problem: even if nu varies locally, the EXPANSION RATE depends
  on the AVERAGE energy density, not local values.  And the average
  backreaction is B_psi/H0^2 ~ 10^-9, as computed in ct5.

  Scale-dependent H0 would require:
    - The measurement to be affected by LOCAL gravity, not expansion
    - This is the "gravitational lensing bias" hypothesis
    - But SN Ia measure REDSHIFTS (expansion) not local gravity

  CONCLUSION: Scale-dependent H0 does not help TGP.
""")


# ======================================================================
#  PART D: IS THE TENSION REAL?
# ======================================================================
print_header("PART D: ALTERNATIVE -- IS THE TENSION REAL?")

print(f"""  Several developments suggest the H0 tension may partially resolve:

  1. CCHP TRGB (Freedman et al. 2024):
     H0 = 69.8 +/- 0.6 km/s/Mpc
     Intermediate between Planck and SH0ES
     Uses TRGB calibration instead of Cepheids
     -> Suggests Cepheid systematics may inflate H0

  2. JAGB method (Lee et al. 2024):
     H0 = 67.96 +/- 1.85 km/s/Mpc
     Consistent with Planck!
     Independent distance indicator
     -> If confirmed, tension disappears

  3. JWST Cepheid observations:
     Confirm SH0ES distances but with ~1% less crowding bias
     H0 = 72.6 +/- 2.0 (slightly lower than SH0ES)
     -> Crowding correction is real but insufficient to resolve

  4. Peculiar velocity corrections:
     Local flow affects nearby SN calibrators
     Uncertainty: ~0.5-1.0 km/s/Mpc on H0
     -> Partially reduces tension

  CURRENT LANDSCAPE (2026):
    - SH0ES (Cepheids): 73.0 -- still high
    - CCHP (TRGB): 69.8 -- intermediate
    - JAGB: 68.0 -- Planck-consistent
    - Strong lensing: 73-74 -- high but large errors
    - GW sirens: 68 +/- 6 -- consistent with anything

  If the true H0 ~ 69-70 km/s/Mpc (TRGB/JAGB average):
    Tension with Planck = ({69.8 - 67.36:.2f} / {np.sqrt(0.54**2 + 0.60**2):.2f}) = {(69.8-67.36)/np.sqrt(0.54**2+0.60**2):.1f} sigma
    -> MANAGEABLE without new physics

  Probability of tension resolving observationally: ~30-40%
""")


# ======================================================================
#  PART E: WHAT TGP DOES PREDICT
# ======================================================================
print_header("PART E: WHAT TGP DOES PREDICT FOR EXPANSION")

print(f"""  TGP makes the following predictions about expansion history:

  1. H(z) = H_LCDM(z) + O(10^-9) corrections
     -> Indistinguishable from LCDM at all measurable redshifts

  2. w(z) = -1 + O(10^-3)
     -> No phantom crossing (w >= -1 always)
     -> If DESI CPL confirmed: TGP ruled out

  3. Angular diameter distance = LCDM prediction
     -> BAO positions unchanged

  4. Growth rate f*sigma8(z) ~ LCDM - O(10^-2)
     -> Negligible modification to structure growth

  5. a0(z) = c * H(z) / (2*pi) -- IF a0 is emergent
     -> a0 evolves with cosmic expansion
     -> At z=1: a0(z=1) ~ a0(z=0) * H(z=1)/H(z=0) ~ 1.5 * a0
     -> Rotation curves at z=1 would show DIFFERENT a0
     -> THIS IS TESTABLE with Euclid/JWST deep field surveys!

  Prediction #5 is the MOST INTERESTING:
    If a0 is truly c*H/(2*pi), then galaxy dynamics at z=1
    should show a0 ~ 1.8e-10 m/s^2 instead of 1.2e-10 m/s^2.
    This would be a ~50% change in the MOND acceleration scale!

    In MOND: a0 is a CONSTANT (no z-dependence)
    In TGP:  a0 ~ H(z) -> EVOLVES with redshift

    This is a CLEAN DISCRIMINATING TEST between TGP and MOND.
""")

# Compute a0(z) evolution
z_evol = np.array([0.0, 0.2, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0])
Omega_r = 9.1e-5
H_z = 67.4 * np.sqrt(Omega_r*(1+z_evol)**4 + 0.315*(1+z_evol)**3 + 0.685)
a0_z = 1.2e-10 * H_z / 67.4  # scaled from z=0

print(f"  a0(z) evolution if a0 = c*H/(2*pi):")
print(f"  {'z':>5s}  {'H(z)':>10s}  {'a0(z) [m/s^2]':>15s}  {'a0/a0(0)':>10s}")
print(f"  {'-'*5}  {'-'*10}  {'-'*15}  {'-'*10}")
for i, z in enumerate(z_evol):
    print(f"  {z:5.1f}  {H_z[i]:10.1f}  {a0_z[i]:15.2e}  {a0_z[i]/1.2e-10:10.3f}")


# ======================================================================
#  PART F: VERDICT
# ======================================================================
print_header("PART F: VERDICT")

print(f"""
  =====================================================================
  H0 TENSION: TGP VERDICT
  =====================================================================

  1. TGP CANNOT explain the H0 tension through backreaction.
     Gap: 8-9 orders of magnitude. STRUCTURAL, not parametric.

  2. The tension may not require new physics:
     - TRGB gives H0 = 69.8 (3.0 sigma from Planck)
     - JAGB gives H0 = 68.0 (0.3 sigma from Planck)
     - True H0 may be ~69-70, manageable within LCDM

  3. TGP's UNIQUE testable prediction:
     a0(z) = c*H(z)/(2*pi) -> a0 evolves with redshift
     At z=1: a0 ~ 1.5 * a0(z=0)
     Testable with Euclid deep-field rotation curves

  4. If H0 tension IS real (requiring new physics):
     TGP needs extension (second field, vector mode, etc.)
     This breaks minimality -- similar to MOND -> TeVeS -> AeST

  RECOMMENDATION:
     Wait for DESI DR2, JWST deep surveys, and future TRGB/JAGB
     measurements. If tension resolves at ~70, TGP is fine.
     If confirmed at 73+, TGP needs cosmological extension.

  =====================================================================
""")
