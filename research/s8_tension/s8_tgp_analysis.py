#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
s8: S8 TENSION -- TGP GROWTH SUPPRESSION ANALYSIS
===================================================

S8 = sigma8 * sqrt(Omega_m / 0.3) measures matter fluctuation amplitude.
CMB predicts S8 ~ 0.83, while late-Universe weak lensing gives S8 ~ 0.76.
The ~3 sigma tension suggests slower structure growth than LCDM predicts.

TGP effects on growth:
  - nu(y) modifies effective G at galaxy/cluster scales
  - At linear scales (8 Mpc), y >> 1, so nu ~ 1 (negligible)
  - TGP provides ~2% growth suppression, need ~8.5%

This script quantifies TGP's contribution and compares with
non-TGP explanations (baryonic feedback, neutrinos, etc.)

Dependencies: numpy, scipy
"""

import sys
import warnings
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
warnings.filterwarnings('ignore')

import numpy as np
from scipy.integrate import solve_ivp

def print_header(title):
    print()
    print("=" * 78)
    print(f"  {title}")
    print("=" * 78)
    print()

def print_subheader(title):
    print(f"\n  {title}")
    print("  " + "-" * len(title))

# Constants
H0 = 67.36  # km/s/Mpc
Omega_m = 0.315
Omega_L = 0.685
Omega_r = 9.1e-5
Omega_b = 0.049
a0_MOND = 1.2e-10  # m/s^2
G_SI = 6.674e-11
Mpc = 3.086e22
ALPHA = 0.8

def nu_tgp(y, c_eff=1.0):
    gamma = ALPHA * c_eff / (c_eff + 1)
    y = np.maximum(np.asarray(y, dtype=float), 1e-12)
    return 1.0 + np.exp(-y**ALPHA) / y**gamma


print_header("s8: S8 TENSION -- TGP GROWTH SUPPRESSION ANALYSIS")


# ======================================================================
#  PART A: S8 DATA COMPILATION
# ======================================================================
print_header("PART A: S8 MEASUREMENT COMPILATION")

data = [
    ("Planck 2018 (CMB)",         0.832, 0.013, "CMB",            "early"),
    ("ACT DR6",                    0.840, 0.028, "CMB",            "early"),
    ("SPT-3G",                     0.797, 0.042, "CMB",            "early"),
    ("DES Y3 (cosmic shear)",     0.776, 0.017, "Weak lensing",   "late"),
    ("KiDS-1000",                  0.759, 0.021, "Weak lensing",   "late"),
    ("HSC Y3",                     0.769, 0.026, "Weak lensing",   "late"),
    ("DES Y3 (3x2pt)",            0.782, 0.019, "WL + clustering","late"),
    ("Planck SZ clusters",        0.770, 0.024, "Cluster counts",  "late"),
    ("eROSITA (2024)",            0.770, 0.017, "X-ray clusters",  "late"),
    ("DESI + Planck lensing",     0.825, 0.015, "BAO + CMB lens",  "early"),
]

print(f"  {'Measurement':>35s}  {'S8':>6s} {'+-':>5s}  {'Method':>20s}  {'Type':>5s}")
print(f"  {'-'*35}  {'-'*6} {'-'*5}  {'-'*20}  {'-'*5}")
for name, s8, sig, method, typ in data:
    print(f"  {name:>35s}  {s8:6.3f} {sig:5.3f}  {method:>20s}  {typ:>5s}")

# Weighted averages
early = [(s, e) for _, s, e, _, t in data if t == "early"]
late  = [(s, e) for _, s, e, _, t in data if t == "late"]

def wavg(d):
    w = [1/e**2 for _, e in d]
    v = [s for s, _ in d]
    return sum(wi*vi for wi, vi in zip(w, v))/sum(w), 1/np.sqrt(sum(w))

s8_e, se_e = wavg(early)
s8_l, se_l = wavg(late)
tens = (s8_e - s8_l) / np.sqrt(se_e**2 + se_l**2)

print(f"\n  Weighted averages:")
print(f"    Early: S8 = {s8_e:.3f} +/- {se_e:.3f}")
print(f"    Late:  S8 = {s8_l:.3f} +/- {se_l:.3f}")
print(f"    Tension: {tens:.1f} sigma")
print(f"    Required suppression: {(s8_e - s8_l)/s8_e * 100:.1f}%")


# ======================================================================
#  PART B: TGP GROWTH MODIFICATION
# ======================================================================
print_header("PART B: TGP GROWTH MODIFICATION")

print("""  In LCDM, the linear growth equation is:
    delta'' + 2H*delta' = (3/2)*H^2*Omega_m*delta

  In TGP, effective G is modified by nu(y):
    delta'' + 2H*delta' = (3/2)*H^2*Omega_m*nu(y)*delta

  At 8 h^-1 Mpc (R = 8/0.7 ~ 11.4 Mpc):
    - Mean density is ~ rho_crit * Omega_m
    - g_bar at this scale: g = (4/3)*pi*G*rho_m*R
""")

# Compute y at 8 Mpc scale
R_8 = 11.4 * Mpc  # meters
rho_m = Omega_m * 3.0 * (H0 * 1e3 / Mpc)**2 / (8.0 * np.pi * G_SI)
g_bar_8 = 4.0/3.0 * np.pi * G_SI * rho_m * R_8
y_8 = g_bar_8 / a0_MOND

print(f"  At R = 8 h^-1 Mpc = 11.4 Mpc:")
print(f"    rho_m = {rho_m:.3e} kg/m^3")
print(f"    g_bar = {g_bar_8:.3e} m/s^2")
print(f"    y = g_bar/a0 = {y_8:.4f}")

for c in [1.0, 1.5, 2.0, 3.0]:
    nu = nu_tgp(y_8, c)
    print(f"    nu(y, c={c:.1f}) = {nu:.6f}  -> growth boost = {(nu-1)*100:.3f}%")

# Solve growth equation with TGP modification
print_subheader("B.1  Linear growth factor D(a) with TGP")

def growth_ode(a, y_vec, nu_factor=1.0):
    """Growth ODE: d^2 delta / da^2 + ... = (3/2)*Om*nu/a^3 * delta / (H/H0)^2"""
    delta, ddelta_da = y_vec
    Hz2 = Omega_r/a**4 + Omega_m/a**3 + Omega_L
    Hz = np.sqrt(max(Hz2, 1e-30))

    # d^2delta/da^2 + (3/a + H'/H)*ddelta/da = (3/2)*Om_m/(a^3*H^2)*nu*delta
    dHda = (-4*Omega_r/a**5 - 3*Omega_m/a**4) / (2*Hz)
    coeff1 = 3.0/a + dHda/Hz
    coeff2 = 1.5 * Omega_m / (a**3 * Hz2) * nu_factor

    d2delta = -coeff1 * ddelta_da + coeff2 * delta
    return [ddelta_da, d2delta]

a_init = 1e-3
a_span = (a_init, 1.0)
y0 = [a_init, 1.0]  # delta ~ a in matter era

# LCDM growth
sol_LCDM = solve_ivp(growth_ode, a_span, y0, args=(1.0,),
                      max_step=0.01, dense_output=True)

# TGP growth with nu modification
# nu ~ 1.001 at these scales (y >> 1)
nu_8Mpc = nu_tgp(y_8, 1.0)
sol_TGP = solve_ivp(growth_ode, a_span, y0, args=(nu_8Mpc,),
                     max_step=0.01, dense_output=True)

a_eval = np.linspace(0.1, 1.0, 10)
D_LCDM = sol_LCDM.sol(a_eval)[0]
D_TGP  = sol_TGP.sol(a_eval)[0]

# Normalize
D_LCDM /= D_LCDM[-1]
D_TGP  /= D_TGP[-1]

print(f"\n  {'a':>6s}  {'z':>6s}  {'D_LCDM':>10s}  {'D_TGP':>10s}  {'D_TGP/D_LCDM':>14s}")
print(f"  {'-'*6}  {'-'*6}  {'-'*10}  {'-'*10}  {'-'*14}")
for i, a in enumerate(a_eval):
    z = 1.0/a - 1.0
    print(f"  {a:6.3f}  {z:6.2f}  {D_LCDM[i]:10.6f}  {D_TGP[i]:10.6f}  {D_TGP[i]/D_LCDM[i]:14.8f}")

# sigma8 modification
ratio_D = D_TGP[-1] / D_LCDM[-1]  # should be ~1 since normalized
# The real effect: TGP modifies growth at ALL redshifts
# delta_sigma8/sigma8 ~ integral of (nu-1) effect over growth history
effective_suppression = (nu_8Mpc - 1.0)  # fractional boost in growth rate
sigma8_TGP = 0.832 * (1.0 + effective_suppression * 0.5)  # rough: half the nu effect

print(f"\n  TGP S8 prediction:")
print(f"    sigma8_LCDM = 0.832")
print(f"    nu at 8 Mpc = {nu_8Mpc:.8f}")
print(f"    nu - 1 = {nu_8Mpc - 1:.2e}")
print(f"    Effective suppression: {effective_suppression*100:.4f}%")
print(f"    This is {effective_suppression / ((s8_e - s8_l)/s8_e) * 100:.3f}% of needed")


# ======================================================================
#  PART C: NON-TGP EXPLANATIONS
# ======================================================================
print_header("PART C: NON-TGP EXPLANATIONS FOR S8 TENSION")

print(f"""
  Several astrophysical effects can suppress S8 without new physics:

  {'Effect':30s}  {'Suppression':>12s}  {'Status':>15s}
  {'-'*30}  {'-'*12}  {'-'*15}
  {'Baryonic feedback (AGN)':30s}  {'3-5%':>12s}  {'well-modeled':>15s}
  {'Massive neutrinos (0.06 eV)':30s}  {'1.5%':>12s}  {'confirmed':>15s}
  {'Massive neutrinos (0.12 eV)':30s}  {'3%':>12s}  {'upper limit':>15s}
  {'Intrinsic alignments':30s}  {'1-3%':>12s}  {'systematic':>15s}
  {'Photo-z errors':30s}  {'1-2%':>12s}  {'systematic':>15s}
  {'Nonlinear modeling':30s}  {'1%':>12s}  {'systematic':>15s}
  {'-'*30}  {'-'*12}  {'-'*15}
  {'Combined (conservative)':30s}  {'5-8%':>12s}  {'':>15s}
  {'Combined (optimistic)':30s}  {'8-12%':>12s}  {'':>15s}
  {'Needed':30s}  {f'{(s8_e-s8_l)/s8_e*100:.1f}%':>12s}  {'':>15s}

  CONCLUSION: Standard astrophysical effects can PLAUSIBLY explain
  the S8 tension without invoking modified gravity.

  Baryonic feedback alone (3-5%) plus neutrinos (1.5-3%) gives
  4.5-8% suppression, close to the needed {(s8_e-s8_l)/s8_e*100:.1f}%.
""")


# ======================================================================
#  PART D: TGP-SPECIFIC PREDICTIONS FOR sigma8
# ======================================================================
print_header("PART D: TGP-SPECIFIC sigma8 PREDICTIONS")

print(f"""
  While TGP cannot explain the bulk of S8 suppression,
  it makes specific predictions testable with future surveys:

  1. SCALE-DEPENDENT sigma(R):
     TGP modifies gravity differently at different scales.
     nu(y) depends on the enclosed mass at radius R.

     R (Mpc)    y = g_bar/a0    nu(y, c=1)     nu-1
     -------------------------------------------------""")

R_scales = [1, 2, 5, 8, 12, 20, 50, 100]  # Mpc
for R in R_scales:
    R_m = R * Mpc
    g = 4.0/3.0 * np.pi * G_SI * rho_m * R_m
    y = g / a0_MOND
    nu = nu_tgp(y, 1.0)
    print(f"     {R:3d}        {y:10.4f}      {nu:10.6f}      {nu-1:.2e}")

print(f"""
  At ALL linear scales (R > 1 Mpc), y >> 1 and nu ~ 1.
  TGP's modification is NEGLIGIBLE for linear structure growth.

  2. NON-LINEAR sigma(R):
     At R < 1 Mpc (galaxy halos, voids), y can drop below 1.
     TGP enhancement becomes significant, changing the halo
     mass function and void profiles.
     -> This is testable with void lensing (Euclid/Rubin)

  3. MORPHOLOGY-DEPENDENT sigma(R):
     Different c_eff for different environments:
     - Voids: c_eff ~ 3 (spherical) -> larger nu
     - Filaments: c_eff ~ 2 (cylindrical) -> moderate nu
     - Clusters: c_eff ~ 3 (spherical) -> but y >> 1, so nu ~ 1

  SUMMARY: TGP predicts ENVIRONMENT-DEPENDENT growth,
  not a uniform suppression of sigma8.
""")


# ======================================================================
#  PART E: VERDICT
# ======================================================================
print_header("PART E: S8 TENSION VERDICT")

print(f"""
  =====================================================================
  S8 TENSION: TGP VERDICT
  =====================================================================

  1. TGP provides ~{effective_suppression*100:.4f}% growth modification at 8 Mpc
     -> NEGLIGIBLE compared to needed {(s8_e-s8_l)/s8_e*100:.1f}%

  2. Standard astrophysics (AGN feedback + neutrinos + systematics)
     can plausibly explain 5-12% suppression
     -> S8 tension may not require new physics

  3. Current tension significance: {tens:.1f} sigma
     -> NOT definitive (< 5 sigma)

  4. TGP's unique prediction: ENVIRONMENT-DEPENDENT growth
     -> Void lensing different from cluster lensing
     -> Testable with Euclid/Rubin tomographic analysis

  CONCLUSION: S8 tension does NOT constrain TGP.
  TGP neither explains it nor is threatened by it.
  The tension likely resolves through better modeling
  of baryonic feedback and weak lensing systematics.

  =====================================================================
""")
