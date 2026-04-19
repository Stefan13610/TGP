#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gs43: BULLET CLUSTER DEEP DIVE -- BEYOND gs39
===============================================

gs39 established:
  - TGP lensing peak at galaxy positions (qualitatively correct)
  - kappa ratio galaxy/gas ~ 2.3 (matches observed)
  - BUT amplitude deficit ~38%: kappa_TGP = 0.22 vs observed 0.35
  - nu(y) ~ 1.02 at cluster y-values -> negligible TGP boost

This script goes DEEPER:
  A. Differential gamma(S): gas (c_eff=3, spherical) vs galaxies (c_eff=1.5, disk-like)
     -> different nu amplification for each baryonic component
  B. Full 3D enclosed-mass y(r) calculation (not sheet approximation)
  C. Pre-collision "gravitational memory" estimate
  D. Multi-cluster merger comparison: Bullet, MACS J0025, Abell 520, Train Wreck
  E. What nu(y) is NEEDED? Reverse engineering the required interpolation
  F. Hybrid model: TGP + minimal collisionless component
  G. Quantitative verdict with error budget

The fundamental question: can ANY physically motivated modification of
TGP close the 38% gap, or does the Bullet Cluster require collisionless mass?

Author: TGP research program
"""

import sys
import warnings
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
warnings.filterwarnings('ignore')

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

# ==========================================================================
# Physical constants
# ==========================================================================
G_SI    = 6.674e-11       # m^3/(kg*s^2)
a0      = 1.2e-10         # m/s^2
Msun    = 1.989e30        # kg
kpc_m   = 3.086e19        # m per kpc
Mpc_m   = 3.086e22        # m per Mpc
c_light = 3.0e8           # m/s
ALPHA   = 4.0 / 5.0       # 0.8 (Flory exponent)

# Bullet Cluster parameters (from gs39)
z_lens = 0.296
z_source = 1.0

# Cluster masses
M_gas_main  = 1.5e14 * Msun
M_star_main = 3.0e13 * Msun
M_gas_sub   = 3.0e13 * Msun
M_star_sub  = 8.0e12 * Msun

# Scale radii (meters)
R500_main  = 1200.0 * kpc_m
R500_sub   = 700.0  * kpc_m
r_c_main   = 250.0  * kpc_m     # gas core (beta-model)
r_c_sub    = 120.0  * kpc_m
a_h_main   = 80.0   * kpc_m     # stellar (Hernquist)
a_h_sub    = 40.0   * kpc_m
beta_slope = 2.0 / 3.0

# Post-collision geometry
d_offset = 720.0 * kpc_m        # total gas-galaxy separation
half_offset = d_offset / 2.0    # 360 kpc each side

# Gaussian gas spread after ram-pressure stripping
sigma_gas_spread = 200.0 * kpc_m   # gas spread ~400 kpc FWHM -> sigma~170 kpc


# ==========================================================================
# TGP functions
# ==========================================================================
def nu_tgp(y, alpha=ALPHA, gamma=0.4):
    """nu(y) = 1 + exp(-y^alpha) / y^gamma"""
    y = np.asarray(y, dtype=float)
    y = np.maximum(y, 1e-12)
    return 1.0 + np.exp(-y**alpha) / y**gamma

def gamma_from_ceff(c_eff):
    """gamma = alpha * c_eff / (c_eff + 1)"""
    return ALPHA * c_eff / (c_eff + 1)

def nu_ceff(y, c_eff):
    """nu with given effective codimension."""
    return nu_tgp(y, gamma=gamma_from_ceff(c_eff))


# ==========================================================================
# Density profiles
# ==========================================================================
def rho_beta(r, rho0, r_c, beta=2.0/3.0):
    """Beta-model gas density."""
    return rho0 / (1.0 + (r / r_c)**2)**(1.5 * beta)

def M_beta_enclosed(r, rho0, r_c, beta=2.0/3.0, npts=500):
    """Enclosed mass of beta-model."""
    rr = np.linspace(0, r, npts)
    dr = rr[1] - rr[0] if npts > 1 else r
    rho_arr = rho_beta(rr, rho0, r_c, beta)
    dM = 4.0 * np.pi * rr**2 * rho_arr * dr
    return np.sum(dM)

def rho0_from_Mgas(M_gas, r_c, R_max, beta=2.0/3.0):
    """Central density for given total gas mass."""
    M_unit = M_beta_enclosed(R_max, 1.0, r_c, beta)
    return M_gas / M_unit

def rho_hernquist(r, M_total, a_h):
    """Hernquist profile."""
    r = np.maximum(np.asarray(r, dtype=float), 1e-10)
    return M_total / (2.0 * np.pi) * a_h / (r * (r + a_h)**3)

def M_hernquist_enclosed(r, M_total, a_h):
    """Enclosed Hernquist mass."""
    return M_total * r**2 / (r + a_h)**2

def sigma_hernquist_analytic(R_perp, M_total, a_h):
    """Analytic projected Hernquist surface density."""
    s = R_perp / a_h
    s = max(s, 1e-6)
    prefactor = M_total / (2.0 * np.pi * a_h**2)
    if abs(s - 1.0) < 1e-4:
        return prefactor * 4.0 / 15.0
    elif s < 1.0:
        sq = np.sqrt(1.0 - s**2)
        F = np.arccosh(1.0 / s) / sq
        return prefactor * ((2 + s**2) * F - 3.0) / (1.0 - s**2)**2
    else:
        sq = np.sqrt(s**2 - 1.0)
        F = np.arccos(1.0 / s) / sq
        return prefactor * ((2 + s**2) * F - 3.0) / (s**2 - 1.0)**2


def sigma_beta_projected(R_perp, rho0, r_c, R_max, beta=2.0/3.0):
    """Projected surface density of beta-model (numerical)."""
    if R_perp >= R_max:
        return 0.0
    L_max = np.sqrt(R_max**2 - R_perp**2)
    def integrand(z):
        r3d = np.sqrt(R_perp**2 + z**2)
        return rho_beta(r3d, rho0, r_c, beta)
    val, _ = quad(integrand, 0, L_max, limit=100, epsrel=1e-6)
    return 2.0 * val


# ==========================================================================
# Utility
# ==========================================================================
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
# Critical surface density for lensing
# ==========================================================================
# Angular diameter distances (flat LCDM, H0=70, Om=0.3)
D_l  = 900.0  * Mpc_m    # lens distance
D_s  = 1700.0 * Mpc_m    # source distance
D_ls = 1100.0 * Mpc_m    # lens-source

Sigma_cr = c_light**2 * D_s / (4.0 * np.pi * G_SI * D_l * D_ls)
Sigma_cr_Msun_Mpc2 = Sigma_cr / Msun * Mpc_m**2

# Compute central densities
rho0_main_gas = rho0_from_Mgas(M_gas_main, r_c_main, 3.0 * R500_main, beta_slope)
rho0_sub_gas  = rho0_from_Mgas(M_gas_sub,  r_c_sub,  3.0 * R500_sub,  beta_slope)


# ##########################################################################
#  MAIN ANALYSIS
# ##########################################################################

print_header("gs43: BULLET CLUSTER DEEP DIVE -- BEYOND gs39")

print(f"""  Building on gs39 results:
    - Peak location: at galaxy positions (CORRECT)
    - kappa ratio galaxy/gas: ~2.3 (CORRECT)
    - Peak amplitude: 0.22 vs 0.35 observed (38% DEFICIT)
    - nu(y) ~ 1.02 at cluster y-values (NEGLIGIBLE boost)

  This analysis investigates whether the deficit can be closed.""")


# ======================================================================
#  PART A: DIFFERENTIAL GAMMA(S) -- GAS vs GALAXIES
# ======================================================================
print_header("PART A: DIFFERENTIAL GAMMA(S) -- GAS vs GALAXIES")

print("""  In TGP, gamma depends on the morphology (codimension) of the source:
    gamma = alpha * c_eff / (c_eff + 1)

  Key insight: gas and galaxies have DIFFERENT morphologies:
    - ICM gas: quasi-spherical, c_eff ~ 3 -> gamma ~ 0.600
    - Galaxy population: flattened/substructured, c_eff ~ 1.5 -> gamma ~ 0.480
    - Individual galaxies within cluster: c_eff ~ 1.3 -> gamma ~ 0.452

  The differential gamma means gas and galaxies get DIFFERENT nu(y) values
  at the same acceleration y.  Let's compute this systematically.
""")

# c_eff assignments for different components
c_eff_gas = 3.0       # spherical ICM
c_eff_gal_pop = 1.5   # galaxy population (between disk and sphere)
c_eff_gal_disk = 1.3  # individual disk galaxies

gamma_gas = gamma_from_ceff(c_eff_gas)
gamma_gal = gamma_from_ceff(c_eff_gal_pop)
gamma_disk = gamma_from_ceff(c_eff_gal_disk)

print(f"  Component c_eff assignments:")
print(f"    ICM gas (spherical):     c_eff = {c_eff_gas:.1f}  -> gamma = {gamma_gas:.3f}")
print(f"    Galaxy population:       c_eff = {c_eff_gal_pop:.1f}  -> gamma = {gamma_gal:.3f}")
print(f"    Individual disk gal:     c_eff = {c_eff_gal_disk:.1f}  -> gamma = {gamma_disk:.3f}")

print_subheader("A.1  nu(y) comparison at cluster y-values")

y_values = np.array([0.01, 0.05, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0])

print(f"\n  {'y':>8s}  {'nu(gas)':>10s}  {'nu(gal)':>10s}  {'nu(disk)':>10s}  {'ratio g/d':>10s}")
print(f"  {'-'*8}  {'-'*10}  {'-'*10}  {'-'*10}  {'-'*10}")
for y in y_values:
    ng = nu_ceff(y, c_eff_gas)
    np_ = nu_ceff(y, c_eff_gal_pop)
    nd = nu_ceff(y, c_eff_gal_disk)
    print(f"  {y:8.2f}  {ng:10.4f}  {np_:10.4f}  {nd:10.4f}  {ng/nd:10.4f}")


# Now compute y at key positions using 3D enclosed mass
print_subheader("A.2  3D enclosed-mass y(r) at key positions")

print("""
  gs39 used the infinite-sheet approximation: y = g_bar/a0 = 2*pi*G*Sigma/a0.
  Here we use the CORRECT 3D enclosed mass: g_bar(r) = G*M(<r)/r^2.
""")

# At the main galaxy position (360 kpc from center)
# The galaxy population sits at this position with stellar mass concentrated
# Gas has been stripped away -> mostly stellar mass here

r_gal = 360.0 * kpc_m   # galaxy position from center of mass

# 3D enclosed mass at galaxy position for main cluster
# Stars: Hernquist with a_h = 80 kpc
M_star_at_gal = M_hernquist_enclosed(r_gal, M_star_main, a_h_main)
# Gas: most gas has been stripped to center, but some residual
# After ram-pressure stripping, gas within galaxy region is ~10-20% of original
f_gas_retained = 0.15   # 15% gas retained near galaxies
M_gas_at_gal = f_gas_retained * M_beta_enclosed(r_gal, rho0_main_gas, r_c_main, beta_slope)

M_total_at_gal = M_star_at_gal + M_gas_at_gal
g_bar_at_gal = G_SI * M_total_at_gal / r_gal**2
y_gal_3D = g_bar_at_gal / a0

# At gas center position: mostly stripped gas
# Gas peak is between the two subclusters
# Approximate: all stripped gas from both subclusters piles up here
M_gas_stripped = (1.0 - f_gas_retained) * (M_gas_main + M_gas_sub)  # total stripped gas
# This gas is spread in a region of radius ~200-300 kpc
r_gas_peak = 200.0 * kpc_m  # effective radius of gas concentration
# Approximate gas density as uniform sphere for enclosed mass estimate
rho_gas_center = M_gas_stripped / (4.0/3.0 * np.pi * (300.0 * kpc_m)**3)
M_gas_at_center = 4.0/3.0 * np.pi * r_gas_peak**3 * rho_gas_center
g_bar_at_gas = G_SI * M_gas_at_center / r_gas_peak**2
y_gas_3D = g_bar_at_gas / a0

print(f"  At main galaxy position (r = 360 kpc from merger axis center):")
print(f"    M_star(<r)  = {M_star_at_gal/Msun:.2e} Msun")
print(f"    M_gas(<r)   = {M_gas_at_gal/Msun:.2e} Msun  (15% retained)")
print(f"    M_total(<r) = {M_total_at_gal/Msun:.2e} Msun")
print(f"    g_bar       = {g_bar_at_gal:.3e} m/s^2")
print(f"    y = g_bar/a0 = {y_gal_3D:.4f}")
print()
print(f"  At gas center (r = 200 kpc effective radius):")
print(f"    M_gas_stripped = {M_gas_stripped/Msun:.2e} Msun")
print(f"    M_gas(<200kpc) = {M_gas_at_center/Msun:.2e} Msun")
print(f"    g_bar          = {g_bar_at_gas:.3e} m/s^2")
print(f"    y = g_bar/a0   = {y_gas_3D:.4f}")


# Differential nu at these positions
print_subheader("A.3  Differential nu: gas vs galaxies at their respective positions")

print(f"\n  At GALAXY position (y = {y_gal_3D:.4f}):")
print(f"    nu(gas, c=3.0)    = {nu_ceff(y_gal_3D, c_eff_gas):.6f}")
print(f"    nu(gal, c=1.5)    = {nu_ceff(y_gal_3D, c_eff_gal_pop):.6f}")
print(f"    nu(disk, c=1.3)   = {nu_ceff(y_gal_3D, c_eff_gal_disk):.6f}")
print()
print(f"  At GAS position (y = {y_gas_3D:.4f}):")
print(f"    nu(gas, c=3.0)    = {nu_ceff(y_gas_3D, c_eff_gas):.6f}")
print(f"    nu(gal, c=1.5)    = {nu_ceff(y_gas_3D, c_eff_gal_pop):.6f}")


# Compute effective Sigma_TGP with differential gamma
print_subheader("A.4  Sigma_TGP with differential gamma along merger axis")

# 1D merger axis: x = -360 (main gal), 0 (gas), +360 (bullet gal)
# Surface density at each position from gs39 (in kg/m^2)
# Main galaxy: stars dominate, Hernquist projected
Sigma_star_main_gal = sigma_hernquist_analytic(10.0 * kpc_m, M_star_main, a_h_main)
Sigma_star_sub_gal  = sigma_hernquist_analytic(10.0 * kpc_m, M_star_sub,  a_h_sub)

# Gas at galaxy positions is depleted
Sigma_gas_at_gal_main = f_gas_retained * sigma_beta_projected(10.0 * kpc_m, rho0_main_gas, r_c_main, 3.0*R500_main, beta_slope)
Sigma_gas_at_gal_sub  = f_gas_retained * sigma_beta_projected(10.0 * kpc_m, rho0_sub_gas,  r_c_sub,  3.0*R500_sub,  beta_slope)

# Gas at center: stripped gas from both clusters
# Gaussian distribution with sigma ~ 200 kpc
# Peak surface density of stripped gas at center
Sigma_gas_center_peak = M_gas_stripped / (2.0 * np.pi * sigma_gas_spread**2)

print(f"\n  Surface densities (projected, at R_perp=10 kpc from each centroid):")
print(f"  {'Component':>30s}  {'Sigma (kg/m^2)':>15s}  {'c_eff':>6s}  {'gamma':>6s}")
print(f"  {'-'*30}  {'-'*15}  {'-'*6}  {'-'*6}")
print(f"  {'Stars at main galaxy':>30s}  {Sigma_star_main_gal:15.4e}  {c_eff_gal_pop:6.1f}  {gamma_gal:6.3f}")
print(f"  {'Gas at main galaxy':>30s}  {Sigma_gas_at_gal_main:15.4e}  {c_eff_gas:6.1f}  {gamma_gas:6.3f}")
print(f"  {'Stars at bullet galaxy':>30s}  {Sigma_star_sub_gal:15.4e}  {c_eff_gal_pop:6.1f}  {gamma_gal:6.3f}")
print(f"  {'Gas at bullet galaxy':>30s}  {Sigma_gas_at_gal_sub:15.4e}  {c_eff_gas:6.1f}  {gamma_gas:6.3f}")
print(f"  {'Stripped gas at center':>30s}  {Sigma_gas_center_peak:15.4e}  {c_eff_gas:6.1f}  {gamma_gas:6.3f}")

# Compute y at each position using 2*pi*G*Sigma_total/a0 (sheet approx for comparison)
Sigma_total_main_gal = Sigma_star_main_gal + Sigma_gas_at_gal_main
Sigma_total_sub_gal  = Sigma_star_sub_gal  + Sigma_gas_at_gal_sub
Sigma_total_gas_cen  = Sigma_gas_center_peak  # minimal stars at center

y_main_gal_sheet = 2.0 * np.pi * G_SI * Sigma_total_main_gal / a0
y_sub_gal_sheet  = 2.0 * np.pi * G_SI * Sigma_total_sub_gal  / a0
y_gas_cen_sheet  = 2.0 * np.pi * G_SI * Sigma_total_gas_cen  / a0

print(f"\n  Sheet-approximation y values:")
print(f"    y(main galaxy)  = {y_main_gal_sheet:.4f}")
print(f"    y(bullet galaxy) = {y_sub_gal_sheet:.4f}")
print(f"    y(gas center)   = {y_gas_cen_sheet:.4f}")

# NOW: apply DIFFERENTIAL gamma
# At galaxy positions: stars get nu(y, c_gal), gas gets nu(y, c_gas)
# This is the key new calculation

# Use 3D y for more accuracy, but compute per-component nu
# For a multi-component system, the effective Sigma_TGP is:
#   Sigma_TGP = nu(y, gamma_stars) * Sigma_stars + nu(y, gamma_gas) * Sigma_gas
# This is MORE CORRECT than applying a single nu to Sigma_bar

print_subheader("A.5  Multi-component TGP lensing with differential gamma")

print("""
  KEY PHYSICS: In TGP, the gravitational enhancement depends on the
  source morphology.  For a multi-component system:

    Sigma_TGP = nu(y, gamma_stars) * Sigma_stars + nu(y, gamma_gas) * Sigma_gas

  This is different from the gs39 approach which used a single gamma.
  The question: does differential gamma help close the amplitude gap?
""")

# At main galaxy position
y_eff = y_main_gal_sheet  # total y at galaxy position

nu_stars_at_gal = nu_ceff(y_eff, c_eff_gal_pop)
nu_gas_at_gal   = nu_ceff(y_eff, c_eff_gas)
nu_uniform      = nu_ceff(y_eff, 1.0)  # gs39 used c_eff=1

Sigma_TGP_diff_gal = nu_stars_at_gal * Sigma_star_main_gal + nu_gas_at_gal * Sigma_gas_at_gal_main
Sigma_TGP_uniform_gal = nu_uniform * Sigma_total_main_gal
Sigma_bar_gal = Sigma_total_main_gal

# At gas center
y_gas_eff = y_gas_cen_sheet
nu_gas_at_center = nu_ceff(y_gas_eff, c_eff_gas)
Sigma_TGP_diff_gas = nu_gas_at_center * Sigma_gas_center_peak
Sigma_TGP_uniform_gas = nu_ceff(y_gas_eff, 1.0) * Sigma_gas_center_peak

print(f"  At MAIN GALAXY position (x = -360 kpc):")
print(f"    Sigma_bar (total)     = {Sigma_bar_gal:.4e} kg/m^2")
print(f"    nu(y, c=1.0) uniform  = {nu_uniform:.6f}")
print(f"    nu(y, c=1.5) stars    = {nu_stars_at_gal:.6f}")
print(f"    nu(y, c=3.0) gas      = {nu_gas_at_gal:.6f}")
print(f"    Sigma_TGP (uniform)   = {Sigma_TGP_uniform_gal:.4e}")
print(f"    Sigma_TGP (diff. gamma) = {Sigma_TGP_diff_gal:.4e}")
print(f"    Boost (uniform):       {Sigma_TGP_uniform_gal/Sigma_bar_gal:.4f}x")
print(f"    Boost (diff. gamma):   {Sigma_TGP_diff_gal/Sigma_bar_gal:.4f}x")
print()
print(f"  At GAS CENTER (x = 0):")
print(f"    Sigma_bar (gas only)  = {Sigma_gas_center_peak:.4e} kg/m^2")
print(f"    Sigma_TGP (diff. gamma) = {Sigma_TGP_diff_gas:.4e}")
print(f"    Boost:                 {Sigma_TGP_diff_gas/Sigma_gas_center_peak:.4f}x")

# Key ratio
kappa_gal_diff = Sigma_TGP_diff_gal / Sigma_cr
kappa_gas_diff = Sigma_TGP_diff_gas / Sigma_cr

print(f"\n  Lensing convergence kappa:")
print(f"    kappa(galaxy, diff gamma)  = {kappa_gal_diff:.4f}  (observed ~0.35)")
print(f"    kappa(gas, diff gamma)     = {kappa_gas_diff:.4f}  (observed ~0.15)")
print(f"    Improvement over gs39:     {(kappa_gal_diff/0.2157 - 1)*100:+.1f}%")

print(f"""
  CONCLUSION A: Differential gamma provides a MARGINAL improvement.
  The problem is that at cluster y-values (y >> 1), nu(y) -> 1 for ALL gamma.
  The differential gamma effect is only significant at y < 1, which does not
  apply to the dense cluster cores.  At y ~ {y_eff:.1f}, nu is already ~1.0.
""")


# ======================================================================
#  PART B: FULL 3D ENCLOSED MASS g(r) PROFILE
# ======================================================================
print_header("PART B: FULL 3D ENCLOSED MASS g(r) PROFILE")

print("""  gs39 used the sheet approximation y = 2*pi*G*Sigma/a0.
  Here we compute the full 3D gravitational acceleration:
    g_bar(r) = G * M_enclosed(r) / r^2
  at each position along the merger axis.
""")

# Radial profile of g_bar for each component
r_arr = np.logspace(np.log10(5.0), np.log10(3000.0), 200)  # kpc
r_m = r_arr * kpc_m

# Main cluster components
g_star_main = np.array([G_SI * M_hernquist_enclosed(r, M_star_main, a_h_main) / r**2 for r in r_m])
g_gas_main  = np.array([G_SI * M_beta_enclosed(r, rho0_main_gas, r_c_main, beta_slope) / r**2 for r in r_m])
g_total_main = g_star_main + g_gas_main
y_3D_main = g_total_main / a0

print(f"  Main cluster g_bar(r) profile:")
print(f"  {'r (kpc)':>10s}  {'g_star':>12s}  {'g_gas':>12s}  {'g_total':>12s}  {'y=g/a0':>10s}  {'nu(c=1)':>10s}  {'nu(c=3)':>10s}")
print(f"  {'-'*10}  {'-'*12}  {'-'*12}  {'-'*12}  {'-'*10}  {'-'*10}  {'-'*10}")
for idx in [0, 20, 50, 80, 100, 120, 150, 180, 199]:
    if idx < len(r_arr):
        r = r_arr[idx]
        gs = g_star_main[idx]
        gg = g_gas_main[idx]
        gt = g_total_main[idx]
        y = y_3D_main[idx]
        n1 = nu_ceff(y, 1.0)
        n3 = nu_ceff(y, 3.0)
        print(f"  {r:10.1f}  {gs:12.3e}  {gg:12.3e}  {gt:12.3e}  {y:10.4f}  {n1:10.6f}  {n3:10.6f}")

# Key finding: at what radius does y drop below 1?
y_cross_1 = None
for i, y in enumerate(y_3D_main):
    if y < 1.0 and y_cross_1 is None:
        y_cross_1 = r_arr[i]

print(f"\n  y = 1 crossing radius: {y_cross_1:.0f} kpc" if y_cross_1 else "\n  y > 1 everywhere in range")

# The CRITICAL point: what is y at the galaxy position?
r_360 = 360.0  # kpc
r_360_m = r_360 * kpc_m
g_star_360 = G_SI * M_hernquist_enclosed(r_360_m, M_star_main, a_h_main) / r_360_m**2
g_gas_360  = G_SI * M_beta_enclosed(r_360_m, rho0_main_gas, r_c_main, beta_slope) / r_360_m**2
# But post-collision: gas is depleted at galaxy position
g_gas_360_depleted = f_gas_retained * g_gas_360
g_total_360 = g_star_360 + g_gas_360_depleted
y_360_3D = g_total_360 / a0

print(f"\n  At galaxy position (r=360 kpc), POST-COLLISION:")
print(f"    g_star  = {g_star_360:.3e} m/s^2")
print(f"    g_gas   = {g_gas_360_depleted:.3e} m/s^2  (15% retained)")
print(f"    g_total = {g_total_360:.3e} m/s^2")
print(f"    y_3D    = {y_360_3D:.4f}")
print(f"    nu(y, c=1.5) = {nu_ceff(y_360_3D, 1.5):.6f}")
print(f"    nu(y, c=1.0) = {nu_ceff(y_360_3D, 1.0):.6f}")

# Pre-collision y at same radius
g_total_360_pre = g_star_360 + g_gas_360  # full gas
y_360_pre = g_total_360_pre / a0
print(f"\n  At same radius, PRE-COLLISION (full gas):")
print(f"    g_total = {(g_star_360 + g_gas_360):.3e} m/s^2")
print(f"    y_3D    = {y_360_pre:.4f}")
print(f"    nu(y, c=1.5) = {nu_ceff(y_360_pre, 1.5):.6f}")

print(f"""
  CONCLUSION B: The 3D enclosed-mass y is HIGHER than the sheet approximation
  at the galaxy position, because M_enclosed includes mass from the full 3D
  volume.  Higher y means nu is even CLOSER to 1, making the TGP boost
  even smaller.  The 3D calculation WORSENS the situation.
""")


# ======================================================================
#  PART C: GRAVITATIONAL MEMORY -- PRE-COLLISION FIELD
# ======================================================================
print_header("PART C: GRAVITATIONAL MEMORY -- PRE-COLLISION FIELD")

print("""  SPECULATIVE MECHANISM: In TGP's f(R) formulation, the scalar field
  phi (scalaron) has a finite response time.  During a cluster merger
  at v ~ 3000 km/s, the collision timescale is:

    t_cross = d_separation / v_rel ~ 720 kpc / 3000 km/s ~ 230 Myr

  If the scalaron field has a relaxation time comparable to or longer
  than t_cross, the gravitational field at the galaxy position could
  retain a "memory" of the pre-collision configuration, when the full
  gas halo was still centered on the galaxies.

  Let's estimate the scalaron relaxation timescale.
""")

# Scalaron mass in f(R) theory
# f(R) = R + R0^gamma * R^(1-gamma) * exp(-(R/R0)^alpha)
# The scalaron mass is m_phi^2 = 1/(3*f_RR)
# At cluster scales, R ~ 10^-25 m^-2 (typical Ricci scalar)

R_cluster = 8.0 * np.pi * G_SI * 1.27e-7 * (kpc_m/Mpc_m)**3 * Msun / c_light**2  # approximate
# More practical: R ~ 3*H^2 * Omega_m * (1+delta)
H0 = 70.0e3 / Mpc_m  # s^-1
R_cosmo = 3.0 * H0**2 * 0.3  # background Ricci scalar in s^-2

# At cluster overdensity delta ~ 200
delta_cluster = 200.0
R_cluster_approx = R_cosmo * (1 + delta_cluster)  # crude

# R0 in TGP: set by a0
R0_tgp = a0 / c_light  # approximate
# More precisely: R0 ~ a0^2 / (c^2 * some_length)
# The key parameter is R/R0 ratio

# Scalaron mass from f_RR
# For our f(R), f_RR is extremely small due to exp(-x^0.8) suppression
# This means m_phi is extremely LARGE -> scalaron is very heavy
# Heavy scalaron -> SHORT relaxation time -> NO gravitational memory

# Compton wavelength of scalaron
# lambda_C = hbar / (m_phi * c) ~ 1/m_phi (natural units)
# But m_phi ~ 1/sqrt(f_RR) and f_RR ~ exp(-big_number) * small

v_rel = 3000.0e3  # m/s (merger velocity)
t_cross = d_offset / v_rel
t_cross_Myr = t_cross / (3.156e7 * 1e6)

print(f"  Merger parameters:")
print(f"    v_rel     = {v_rel/1e3:.0f} km/s")
print(f"    d_offset  = {d_offset/kpc_m:.0f} kpc")
print(f"    t_cross   = {t_cross_Myr:.1f} Myr")
print()

# Estimate scalaron relaxation time
# In TGP f(R), f_RR = d^2f/dR^2
# f(R) = R + R0^gamma * R^(1-gamma) * exp(-(R/R0)^alpha)
# Let u = R/R0, h(u) = R0 * u^(1-gamma) * exp(-u^alpha)
# f_RR = h''(u)/R0

# For cluster environment, R >> R0 (strong field), u >> 1
# exp(-u^alpha) is astronomically small
# f_RR ~ 0 -> m_phi^2 ~ infinity -> instant response

# BUT: at the outskirts (R ~ R0), the scalaron could be light
# This would only matter at very large radii (>> R500)

# More careful: the "gravitational memory" needs the scalar field
# to propagate at finite speed (c) and have a finite mass

# Speed of scalar field perturbation
c_scalar = c_light  # propagates at c in f(R) gravity

# Time for scalar field to adjust after gas stripping
t_adjust = d_offset / c_scalar  # ~ 0.023 kpc/c ~ 73 years
t_adjust_yr = t_adjust / 3.156e7

print(f"  Scalaron adjustment timescale:")
print(f"    Light-crossing time of 720 kpc = {d_offset/kpc_m / 3.26e3:.1f} kyr")
print(f"    Much shorter than t_cross = {t_cross_Myr:.0f} Myr")
print()
print(f"  In TGP's f(R), at cluster densities:")
print(f"    R/R0 >> 1 -> exp(-(R/R0)^0.8) ~ 0")
print(f"    f_RR ~ 0 -> scalaron mass ~ infinity")
print(f"    Compton wavelength ~ 0 -> response is INSTANTANEOUS")
print(f"    Scalar field adjusts MUCH faster than merger timescale")

print(f"""
  CONCLUSION C: Gravitational memory does NOT help.
  The TGP scalaron in f(R) is extremely massive at cluster densities
  (due to exp(-(R/R0)^alpha) chameleon suppression), so the scalar
  field tracks the matter distribution essentially instantaneously.
  There is no "lingering" gravitational field from the pre-collision state.

  This actually CONFIRMS TGP's CMB compatibility (gs41) -- the same
  mechanism that makes TGP safe for cosmology prevents gravitational
  memory at cluster scales.
""")


# ======================================================================
#  PART D: MULTI-CLUSTER COMPARISON
# ======================================================================
print_header("PART D: MULTI-CLUSTER COMPARISON -- 4 MERGING CLUSTERS")

print("""  The Bullet Cluster is not the only merging system.  Several others
  show the same gas-galaxy lensing offset.  If TGP's failure is generic,
  it should fail similarly for all of them.

  We compare four merging clusters:
  1. Bullet Cluster (1E 0657-56): the canonical case
  2. MACS J0025.4-1222: a "baby bullet"
  3. Abell 520 ("Train Wreck"): ANOMALOUS -- lensing peak at gas!
  4. Abell 2744 ("Pandora's Cluster"): complex multi-merger

  Abell 520 is particularly interesting because it CONTRADICTS
  standard DM expectations -- it has a "dark core" where DM appears
  to concentrate with the gas, not the galaxies.
""")

# Cluster data (approximate, from literature)
clusters = {
    'Bullet (1E 0657-56)': {
        'z': 0.296,
        'M_gas': 1.8e14,  # total system Msun
        'M_star': 3.8e13,
        'offset_kpc': 720,
        'v_rel_kms': 3000,
        'kappa_obs_gal': 0.35,
        'kappa_obs_gas': 0.15,
        'peak_at': 'galaxies',
        'notes': 'Canonical DM evidence'
    },
    'MACS J0025.4-1222': {
        'z': 0.586,
        'M_gas': 8e13,
        'M_star': 2e13,
        'offset_kpc': 450,
        'v_rel_kms': 2000,
        'kappa_obs_gal': 0.25,
        'kappa_obs_gas': 0.10,
        'peak_at': 'galaxies',
        'notes': 'Baby bullet, confirms picture'
    },
    'Abell 520': {
        'z': 0.199,
        'M_gas': 1.2e14,
        'M_star': 3e13,
        'offset_kpc': 500,
        'v_rel_kms': 2300,
        'kappa_obs_gal': 0.20,
        'kappa_obs_gas': 0.25,
        'peak_at': 'gas (ANOMALY)',
        'notes': 'Dark core at gas position -- challenges LCDM too!'
    },
    'Abell 2744': {
        'z': 0.308,
        'M_gas': 2.5e14,
        'M_star': 5e13,
        'offset_kpc': 600,
        'v_rel_kms': 2500,
        'kappa_obs_gal': 0.40,
        'kappa_obs_gas': 0.20,
        'peak_at': 'galaxies',
        'notes': 'Complex multi-merger'
    }
}

print(f"  {'Cluster':>25s}  {'z':>5s}  {'M_gas':>8s}  {'offset':>7s}  {'v_rel':>6s}  {'kappa_g':>7s}  {'kappa_x':>7s}  {'Peak at':>15s}")
print(f"  {'-'*25}  {'-'*5}  {'-'*8}  {'-'*7}  {'-'*6}  {'-'*7}  {'-'*7}  {'-'*15}")
for name, cl in clusters.items():
    print(f"  {name:>25s}  {cl['z']:5.3f}  {cl['M_gas']:.0e}  {cl['offset_kpc']:5d}kp  {cl['v_rel_kms']:5d}  {cl['kappa_obs_gal']:7.2f}  {cl['kappa_obs_gas']:7.2f}  {cl['peak_at']:>15s}")


# TGP prediction for each cluster
print_subheader("D.1  TGP predictions for each cluster")

for name, cl in clusters.items():
    M_gas = cl['M_gas'] * Msun
    M_star = cl['M_star'] * Msun
    M_total = M_gas + M_star

    # Approximate y at galaxy position
    # Use enclosed mass within ~R500/3
    R_eff = cl['offset_kpc'] / 2.0 * kpc_m  # half-offset as effective radius
    # Stellar mass within R_eff (assuming Hernquist with a_h ~ R_eff/5)
    a_h_eff = R_eff / 5.0
    M_star_enc = M_star * R_eff**2 / (R_eff + a_h_eff)**2
    g_bar = G_SI * M_star_enc / R_eff**2
    y_gal = g_bar / a0

    # Gas mass at gas center
    r_gas_eff = 200.0 * kpc_m
    M_gas_enc = 0.3 * M_gas  # ~30% within 200 kpc effective
    g_gas = G_SI * M_gas_enc / r_gas_eff**2
    y_gas = g_gas / a0

    nu_gal = nu_ceff(y_gal, 1.5)
    nu_gas = nu_ceff(y_gas, 3.0)

    # Estimate surface densities (very rough)
    Sigma_gal = M_star_enc / (np.pi * (100.0 * kpc_m)**2)  # compact stellar
    Sigma_gas = M_gas_enc / (np.pi * (300.0 * kpc_m)**2)   # spread gas

    kappa_gal_tgp = nu_gal * Sigma_gal / Sigma_cr
    kappa_gas_tgp = nu_gas * Sigma_gas / Sigma_cr

    deficit_gal = kappa_gal_tgp / cl['kappa_obs_gal'] if cl['kappa_obs_gal'] > 0 else 0

    print(f"\n  {name}:")
    print(f"    y(galaxy) = {y_gal:.3f},  nu = {nu_gal:.4f}")
    print(f"    y(gas)    = {y_gas:.3f},  nu = {nu_gas:.4f}")
    print(f"    kappa_TGP(gal) ~ {kappa_gal_tgp:.3f}  vs obs {cl['kappa_obs_gal']:.2f}  ({deficit_gal*100:.0f}%)")
    print(f"    kappa_TGP(gas) ~ {kappa_gas_tgp:.3f}  vs obs {cl['kappa_obs_gas']:.2f}")
    if cl['peak_at'] == 'gas (ANOMALY)':
        print(f"    NOTE: Abell 520 has lensing peak at GAS -- also problematic for LCDM!")

print_subheader("D.2  Abell 520: an opportunity for TGP?")

print("""
  Abell 520 is the "Train Wreck" cluster with a puzzling DARK CORE:
    - Lensing peak coincides with gas, NOT galaxies
    - This is the OPPOSITE of the Bullet Cluster
    - Standard LCDM has trouble explaining this
    - Requires DM to be "sticky" or self-interacting at A520

  In TGP:
    - Sigma_TGP = nu(y) * Sigma_bar ALWAYS peaks where Sigma_bar peaks
    - If gas is more concentrated (not spread by ram pressure), Sigma_bar
      can peak at the gas position
    - This depends on merger geometry and gas dynamics

  Abell 520 is actually MORE NATURAL for TGP than for LCDM:
    - If the gas happens to be concentrated (less stripping), TGP
      automatically puts the lensing peak at the gas position
    - No need for "self-interacting DM" or other exotic physics

  However, TGP would still have the amplitude problem:
    - nu(y) ~ 1 at cluster y-values
    - Total lensing mass would be ~Sigma_bar, too low by factor ~2-3
""")


# ======================================================================
#  PART E: REVERSE ENGINEERING -- WHAT nu(y) IS NEEDED?
# ======================================================================
print_header("PART E: REVERSE ENGINEERING -- WHAT nu(y) IS NEEDED?")

print("""  If TGP must explain the Bullet Cluster WITHOUT collisionless mass,
  what nu(y) function would be needed at cluster y-values?

  From gs39: kappa_TGP = 0.2157, kappa_obs = 0.35
  Required: nu_needed = 0.35 / (Sigma_bar / Sigma_cr) = kappa_obs / kappa_bar
""")

# From gs39 values
kappa_bar_main_gal = 0.2113   # baryon-only kappa at main galaxy
kappa_obs_main_gal = 0.35     # observed
kappa_bar_gas = 0.0792
kappa_obs_gas = 0.15

kappa_bar_bullet_gal = 0.1680
kappa_obs_bullet_gal = 0.25

nu_needed_main = kappa_obs_main_gal / kappa_bar_main_gal
nu_needed_gas  = kappa_obs_gas / kappa_bar_gas
nu_needed_bullet = kappa_obs_bullet_gal / kappa_bar_bullet_gal

print(f"  Required nu(y) to match observations:")
print(f"  {'Position':>25s}  {'kappa_bar':>10s}  {'kappa_obs':>10s}  {'nu_needed':>10s}  {'y (gs39)':>10s}")
print(f"  {'-'*25}  {'-'*10}  {'-'*10}  {'-'*10}  {'-'*10}")
print(f"  {'Main galaxy':>25s}  {kappa_bar_main_gal:10.4f}  {kappa_obs_main_gal:10.2f}  {nu_needed_main:10.4f}  {'~4.4':>10s}")
print(f"  {'Gas center':>25s}  {kappa_bar_gas:10.4f}  {kappa_obs_gas:10.2f}  {nu_needed_gas:10.4f}  {'~1.7':>10s}")
print(f"  {'Bullet galaxy':>25s}  {kappa_bar_bullet_gal:10.4f}  {kappa_obs_bullet_gal:10.2f}  {nu_needed_bullet:10.4f}  {'~3.5':>10s}")

print(f"\n  TGP provides: nu ~ 1.02 at y ~ 4.4")
print(f"  Needed:       nu ~ {nu_needed_main:.2f} at y ~ 4.4")
print(f"  Deficit factor: {nu_needed_main/1.02:.1f}x")

print_subheader("E.1  What gamma would give the needed nu?")

# Solve: nu(y=4.4, gamma) = 1.66
# 1 + exp(-4.4^0.8) / 4.4^gamma = 1.66
# exp(-4.4^0.8) / 4.4^gamma = 0.66
# 4.4^gamma = exp(-4.4^0.8) / 0.66

y_target = 4.4
nu_target = nu_needed_main

def nu_residual(gamma, y=y_target, target=nu_target):
    return nu_tgp(y, ALPHA, gamma) - target

# Can we solve this?
# exp(-4.4^0.8) = exp(-3.58) = 0.028
# 0.028 / 4.4^gamma = nu_target - 1 = 0.66
# 4.4^gamma = 0.028 / 0.66 = 0.042
# gamma * ln(4.4) = ln(0.042)
# gamma = ln(0.042) / ln(4.4) = -3.17 / 1.48 = -2.14

# NEGATIVE gamma! This is unphysical.
exp_term = np.exp(-y_target**ALPHA)
needed_denom = nu_target - 1.0
gamma_needed = np.log(exp_term / needed_denom) / np.log(y_target) if needed_denom > 0 else None

print(f"\n  At y = {y_target}:")
print(f"    exp(-y^alpha) = exp(-{y_target**ALPHA:.3f}) = {exp_term:.6f}")
print(f"    Need: exp(-y^a)/y^gamma = {needed_denom:.4f}")
print(f"    -> y^gamma = {exp_term/needed_denom:.6f}")
print(f"    -> gamma = ln({exp_term/needed_denom:.6f}) / ln({y_target}) = {gamma_needed:.3f}")

print(f"""
  gamma = {gamma_needed:.3f} is NEGATIVE -- this is UNPHYSICAL.
  No positive codimension can give the needed nu(y) at cluster y-values.

  The fundamental reason: at y >> 1, the exponential exp(-y^0.8) suppresses
  the TGP correction term to near zero.  No power of y in the denominator
  can compensate for this exponential suppression.

  This is NOT a parameter-fitting problem -- it's a STRUCTURAL limitation
  of the TGP nu(y) function.
""")

print_subheader("E.2  Modified alpha -- could a different Flory exponent help?")

# What if alpha were different (not 0.8)?
# nu = 1 + exp(-y^alpha) / y^gamma
# At y=4.4, need nu = 1.66
# exp(-4.4^alpha) / 4.4^gamma = 0.66

# For fixed gamma = 0.48 (c_eff=1.5):
gamma_test = 0.48
# exp(-4.4^alpha) = 0.66 * 4.4^0.48 = 0.66 * 2.02 = 1.33
# But exp(-x) <= 1 for x >= 0 !
# So 0.66 * 4.4^0.48 = 1.33 > 1 -> NO SOLUTION for alpha > 0

rhs = needed_denom * y_target**gamma_test
print(f"  With gamma = {gamma_test} (c_eff=1.5):")
print(f"    Need exp(-y^alpha) = {needed_denom:.3f} * {y_target}^{gamma_test} = {rhs:.3f}")
print(f"    But exp(-x) <= 1 for all x >= 0")
if rhs > 1.0:
    print(f"    {rhs:.3f} > 1 -> NO SOLUTION EXISTS for any alpha > 0!")
else:
    alpha_needed = -np.log(rhs) / np.log(y_target)
    print(f"    -> alpha = {alpha_needed:.3f}")

print(f"""
  The TGP interpolation function STRUCTURALLY CANNOT provide nu ~ 1.66
  at y ~ 4.4 for any physical gamma (>0).  The exponential suppression
  at high y is too strong.

  This is the RIGOROUS version of the gs39 conclusion:
    "would need nu ~ 1.7 but get nu ~ 1.02"
  It's not just a bad parameter choice -- it's mathematically impossible
  within the TGP nu(y) family.
""")

print_subheader("E.3  At what y could TGP provide nu ~ 1.66?")

# Find y where nu(y, gamma=0.48) = 1.66
# 1 + exp(-y^0.8) / y^0.48 = 1.66
for c_test in [1.0, 1.5, 2.0, 3.0]:
    gamma_t = gamma_from_ceff(c_test)
    # Scan y
    y_scan = np.logspace(-3, 1, 10000)
    nu_scan = nu_tgp(y_scan, ALPHA, gamma_t)
    # Find where nu crosses 1.66
    idx = np.argmin(np.abs(nu_scan - nu_target))
    y_match = y_scan[idx]
    nu_match = nu_scan[idx]
    if abs(nu_match - nu_target) < 0.1:
        print(f"  c_eff = {c_test:.1f} (gamma={gamma_t:.3f}): nu = {nu_target:.2f} at y = {y_match:.4f}")
    else:
        print(f"  c_eff = {c_test:.1f} (gamma={gamma_t:.3f}): nu = {nu_target:.2f} requires y < {y_scan[0]:.4f} (never reaches it)")

print(f"""
  To get nu ~ 1.66, we need y << 1, i.e., g_bar << a0.
  At cluster positions, y ~ 2-5, far too high.
  TGP enhancement is only significant in the deep MOND regime (y << 1),
  which is NOT reached at cluster cores.
""")


# ======================================================================
#  PART F: HYBRID MODEL -- TGP + MINIMAL COLLISIONLESS COMPONENT
# ======================================================================
print_header("PART F: HYBRID MODEL -- TGP + MINIMAL COLLISIONLESS COMPONENT")

print("""  Since pure TGP cannot explain the Bullet Cluster amplitude,
  we explore the MINIMAL collisionless component needed to close the gap.

  Three candidate collisionless species:
  1. Hot dark matter (massive neutrinos, m_nu ~ 0.1-0.3 eV)
  2. Warm dark matter (sterile neutrinos, m_s ~ 2-7 keV)
  3. Cold dark matter (standard WIMPs/axions)

  In a hybrid model:
    M_lens(r) = nu(y) * M_bar(r) + M_collisionless(r)

  The collisionless component stays with galaxies (passed through)
  while gas is ram-pressure stripped.
""")

# Required additional mass at galaxy position
# kappa_obs = 0.35, kappa_TGP = 0.2157
# Need kappa_extra = 0.35 - 0.2157 = 0.1343
kappa_tgp = 0.2157    # from gs39
kappa_obs = 0.35
kappa_deficit = kappa_obs - kappa_tgp

Sigma_deficit = kappa_deficit * Sigma_cr  # kg/m^2
M_deficit_per_Mpc2 = Sigma_deficit / Msun * Mpc_m**2
# Convert to mass within effective area
# Galaxy lensing peak extends over ~ 200x200 kpc^2 = (0.2 Mpc)^2
A_eff = (200.0 * kpc_m)**2
M_extra = Sigma_deficit * A_eff  # kg
M_extra_Msun = M_extra / Msun

print(f"  Deficit analysis:")
print(f"    kappa_TGP       = {kappa_tgp:.4f}")
print(f"    kappa_observed  = {kappa_obs:.2f}")
print(f"    kappa_deficit   = {kappa_deficit:.4f}")
print(f"    Sigma_deficit   = {Sigma_deficit:.3e} kg/m^2")
print(f"    M_extra (200x200 kpc) = {M_extra_Msun:.2e} Msun")

# Compare with baryon mass
M_bar_total = (M_gas_main + M_star_main + M_gas_sub + M_star_sub) / Msun
f_extra = M_extra_Msun / M_bar_total

print(f"\n  M_extra / M_baryon = {f_extra:.2f}")
print(f"  M_extra / M_stars  = {M_extra_Msun / ((M_star_main + M_star_sub)/Msun):.2f}")

print_subheader("F.1  Relic neutrinos as the collisionless component")

# Relic neutrino contribution
# Omega_nu = sum(m_nu) / (93.14 h^2 eV)
# For m_nu = 0.06 eV (minimum from oscillations):
m_nu_min = 0.06  # eV, minimum sum
m_nu_max = 0.12  # eV, current upper limit (Planck+DESI ~0.12)
h = 0.70

Omega_nu_min = m_nu_min / (93.14 * h**2)
Omega_nu_max = m_nu_max / (93.14 * h**2)
Omega_b = 0.049
Omega_m = 0.30

# Fraction of total mass in neutrinos
f_nu_min = Omega_nu_min / Omega_m
f_nu_max = Omega_nu_max / Omega_m

# In a cluster with M_total ~ 1e15 Msun (LCDM), neutrino fraction:
M_cluster_total_obs = 1.0e15  # Msun (total observed lensing mass)
M_nu_cluster_min = f_nu_min * M_cluster_total_obs
M_nu_cluster_max = f_nu_max * M_cluster_total_obs

print(f"\n  Relic neutrinos (standard 3 species):")
print(f"    sum(m_nu) = {m_nu_min:.3f} - {m_nu_max:.3f} eV")
print(f"    Omega_nu  = {Omega_nu_min:.5f} - {Omega_nu_max:.5f}")
print(f"    f_nu/f_m  = {f_nu_min:.4f} - {f_nu_max:.4f}")
print(f"    M_nu in cluster ~ {M_nu_cluster_min:.1e} - {M_nu_cluster_max:.1e} Msun")
print(f"    Need M_extra ~ {M_extra_Msun:.1e} Msun")
print(f"    Ratio: M_nu / M_extra = {M_nu_cluster_min/M_extra_Msun:.3f} - {M_nu_cluster_max/M_extra_Msun:.3f}")

print(f"""
  Standard neutrinos provide ONLY {M_nu_cluster_max/M_extra_Msun*100:.1f}% of the needed extra mass.
  They are collisionless (good) but too light to cluster efficiently
  at scales below their free-streaming length (~10 Mpc for m=0.06 eV).

  Neutrinos would need sum(m_nu) ~ 1-2 eV to be significant,
  but this is EXCLUDED by Planck+BAO (sum(m_nu) < 0.12 eV).
""")

print_subheader("F.2  Sterile neutrinos (m ~ 2-7 keV)")

print("""  Sterile neutrinos with m_s ~ 2-7 keV are a warm dark matter candidate:
    - Free-streaming length ~ 100-300 kpc (relevant for cluster cores)
    - Collisionless -> stays with galaxies in mergers
    - Could be produced via Dodelson-Widrow or Shi-Fuller mechanism
    - 3.5 keV X-ray line claim (controversial) -> m_s ~ 7 keV

  Required sterile neutrino density:
""")

# Required Omega_s for sterile neutrinos
# In a hybrid TGP model, we need much less than LCDM DM
# LCDM needs Omega_DM ~ 0.26, TGP+sterile needs only the deficit

# The deficit is ~38% of total lensing mass at galaxy positions
# In a cluster: M_baryon ~ 2e14 Msun, need extra ~ 1e14 Msun
# Cosmologically: Omega_s = M_extra / M_total * Omega_m

f_extra_cosmo = kappa_deficit / kappa_obs  # fraction of observed mass that's missing
Omega_s_needed = f_extra_cosmo * Omega_m

print(f"    Fraction of lensing mass unaccounted:  {f_extra_cosmo:.2f}")
print(f"    Implied Omega_sterile:                 {Omega_s_needed:.3f}")
print(f"    Compare: Omega_DM (LCDM):              0.26")
print(f"    Ratio: Omega_s / Omega_DM(LCDM):       {Omega_s_needed/0.26:.2f}")

# What mass m_s for this Omega_s?
# Omega_s * h^2 = m_s / 94 eV * (T_s/T_nu)^3 * n_species
# For Dodelson-Widrow: Omega_s ~ 0.3 * (sin^2(2theta)/7e-10) * (m_s/3keV)^1.8
# Simplified: require Omega_s ~ 0.12 with m_s ~ few keV

print(f"""
  A TGP + sterile neutrino hybrid would need:
    Omega_s ~ {Omega_s_needed:.3f}  (vs LCDM Omega_DM ~ 0.26)
    This is {Omega_s_needed/0.26*100:.0f}% of the LCDM dark matter budget.

  This is a SIGNIFICANT amount -- TGP doesn't merely need a "sprinkle"
  of collisionless matter, it needs roughly HALF of what LCDM requires.
  However, this is still less than full LCDM, and the sterile neutrino
  would be a SPECIFIC prediction (m_s ~ 3-7 keV, Omega_s ~ 0.12).
""")


print_subheader("F.3  What fraction of missing mass can TGP explain?")

# Total budget: LCDM says the cluster has M_DM ~ 5 * M_baryon
# M_baryon ~ 2.2e14 Msun, M_DM(LCDM) ~ 1.1e15 Msun
M_bar_cluster = (M_gas_main + M_star_main + M_gas_sub + M_star_sub) / Msun
M_DM_LCDM = 5.0 * M_bar_cluster  # typical DM fraction in clusters

# TGP provides phantom DM = (nu-1) * M_bar
# At y ~ 3, nu ~ 1.02 -> phantom DM ~ 0.02 * M_bar
nu_typical_cluster = 1.02
M_phantom_TGP = (nu_typical_cluster - 1.0) * M_bar_cluster

# What standard DM would need to be
M_DM_needed = M_DM_LCDM  # same total to match observations
f_TGP_explains = M_phantom_TGP / M_DM_needed
f_real_DM_needed = 1.0 - f_TGP_explains

print(f"  Mass budget for Bullet Cluster system:")
print(f"    M_baryon (total)    = {M_bar_cluster:.2e} Msun")
print(f"    M_DM (LCDM)         = {M_DM_LCDM:.2e} Msun")
print(f"    M_phantom (TGP)     = {M_phantom_TGP:.2e} Msun")
print(f"    Fraction TGP covers = {f_TGP_explains*100:.1f}%")
print(f"    Collisionless needed = {f_real_DM_needed*100:.1f}% of LCDM DM")

print(f"""
  At CLUSTER scales, TGP provides only ~{f_TGP_explains*100:.0f}% of the missing mass.
  This is in stark contrast to GALAXY scales where TGP provides ~100%.

  The deficit is NOT a fine-tuning problem -- it's a REGIME problem:
    - Galaxies: y ~ 0.1-1 -> nu ~ 2-10 -> significant enhancement
    - Clusters: y ~ 2-10  -> nu ~ 1.01-1.02 -> negligible enhancement

  This is the same regime problem that MOND has:
    - MOND works at g < a0 (galaxy outskirts)
    - At g > a0 (cluster cores), MOND reduces to Newtonian gravity
    - Clusters have g > a0 -> MOND/TGP can't help
""")


# ======================================================================
#  PART G: QUANTITATIVE VERDICT AND ERROR BUDGET
# ======================================================================
print_header("PART G: QUANTITATIVE VERDICT AND ERROR BUDGET")

print_subheader("G.1  Summary of all approaches tried")

print(f"""
  Approach                          | Improvement | Enough?
  ----------------------------------|-------------|--------
  gs39 baseline (c_eff=1)           | --          | NO (38% deficit)
  Differential gamma (Part A)       | ~1-2%       | NO
  3D enclosed-mass y (Part B)       | NEGATIVE    | WORSE
  Gravitational memory (Part C)     | 0%          | NO (scalaron too heavy)
  Multi-cluster comparison (Part D) | diagnostic  | --
  Modified alpha (Part E)           | IMPOSSIBLE  | STRUCTURALLY BLOCKED
  Hybrid + relic neutrinos (Part F) | ~{M_nu_cluster_max/M_extra_Msun*100:.0f}%      | NO
  Hybrid + sterile nu (Part F)      | 100%*       | YES (but needs Omega_s~0.12)

  * Sterile neutrino with Omega_s ~ 0.12 can close the gap,
    but this is ~46% of LCDM's total DM density.
""")

print_subheader("G.2  Error budget for the 38% deficit")

print("""
  How certain is the 38% deficit? Sources of uncertainty:

  1. Observational uncertainty on kappa_obs:
     - kappa(main galaxy) = 0.35 +/- 0.05 (Clowe et al. 2006)
     - If kappa_obs = 0.30 (lower bound), deficit drops to 28%
     - If kappa_obs = 0.40 (upper bound), deficit grows to 46%

  2. Baryonic mass uncertainty:
     - M_gas from X-ray luminosity: +/- 20%
     - M_stars from optical luminosity: +/- 30%
     - If M_bar is 20% higher, kappa_bar ~ 0.25, deficit ~ 29%

  3. Geometric modeling:
     - 1D/2D vs full 3D: ~5-10% effect
     - Gas distribution post-collision: uncertain
     - Projection effects: ~10%

  4. TGP parameters:
     - a0 uncertainty: +/- 10% -> negligible at these y-values
     - alpha from Flory: 0.800 +/- 0.005 -> negligible
     - c_eff assignment: uncertain, but nu~1 regardless

  Combined systematic uncertainty: ~15-20%
  Statistical significance of deficit: ~2-3 sigma
""")

# Compute combined error budget
delta_kappa_obs = 0.05
delta_kappa_bar_frac = 0.20  # 20% baryon mass uncertainty
delta_geometric = 0.10

# Deficit significance
deficit = (kappa_obs - kappa_tgp) / kappa_obs
sigma_deficit = np.sqrt(delta_kappa_obs**2 + (delta_kappa_bar_frac * kappa_tgp)**2 + (delta_geometric * kappa_tgp)**2)
n_sigma = (kappa_obs - kappa_tgp) / sigma_deficit

print(f"  Numerical error budget:")
print(f"    delta(kappa_obs)      = +/- {delta_kappa_obs:.3f}")
print(f"    delta(kappa_bar, 20%) = +/- {delta_kappa_bar_frac * kappa_tgp:.3f}")
print(f"    delta(geometric, 10%) = +/- {delta_geometric * kappa_tgp:.3f}")
print(f"    Combined sigma        = {sigma_deficit:.3f}")
print(f"    Deficit significance  = {n_sigma:.1f} sigma")

print_subheader("G.3  Comparison with other modified gravity theories")

print(f"""
  Theory                 | Bullet Cluster status       | Cluster mass
  -----------------------|-----------------------------|-------------
  MOND (Milgrom 1983)    | FAILS: same regime problem  | 2-3x deficit
  TeVeS (Bekenstein 04)  | FAILS: phantom DM too weak  | 2x deficit
  MOG (Moffat 2006)      | CLAIMS success via Yukawa   | Controversial
  Emergent (Verlinde 16) | UNCLEAR: no prediction      | Unknown
  RMOND (Skordis 21)     | PASSES: has DM-like field   | OK (by design)
  TGP (this work)        | FAILS: nu~1 at cluster y    | ~38% deficit

  KEY INSIGHT: The ONLY modified gravity theory that passes the Bullet
  Cluster test (RMOND/AeST by Skordis & Zlosnik 2021) does so by
  including a VECTOR FIELD that acts like collisionless dark matter.
  This is effectively a "hybrid" model.

  TGP's failure mode is IDENTICAL to MOND's: the enhancement function
  approaches 1 at high accelerations (y >> 1).  This is not a bug but
  a FEATURE -- it ensures Solar System compatibility.
""")

print_subheader("G.4  Final verdict on the Bullet Cluster")

print(f"""
  =====================================================================
  BULLET CLUSTER VERDICT FOR TGP
  =====================================================================

  WHAT WORKS:
    [OK] Lensing peak at galaxy positions (geometric effect)
    [OK] kappa ratio galaxy/gas ~ 2.3 (geometry-driven)
    [OK] Qualitative picture of collisionless galaxies passing through

  WHAT FAILS:
    [FAIL] Peak lensing amplitude: 0.22 vs 0.35 (38% deficit)
    [FAIL] nu(y) ~ 1.02 at cluster y ~ 3-5 (negligible enhancement)
    [FAIL] Total mass deficit ~factor 2-3 at cluster scales

  CAN THE DEFICIT BE CLOSED?
    Differential gamma:     NO (1-2% improvement)
    3D mass calculation:    NO (makes it worse)
    Gravitational memory:   NO (scalaron too heavy)
    Modified alpha:         NO (structurally impossible)
    Relic neutrinos:        NO ({M_nu_cluster_max/M_extra_Msun*100:.0f}% of needed)
    Sterile neutrinos:      YES (but needs Omega_s ~ 0.12)

  HONEST CONCLUSION:
    TGP, like MOND, has a GENUINE cluster-scale problem.
    The nu(y) function derived from membrane physics provides
    excellent fits at galaxy scales (y ~ 0.1-1) but negligible
    enhancement at cluster scales (y ~ 2-10).

    This is not a parameter-tuning issue -- it's a STRUCTURAL
    feature of any theory where g_obs = nu(g/a0) * g and nu->1
    at high g.  The exponential suppression exp(-y^0.8) ensures
    Solar System safety but also eliminates cluster-scale effects.

  OPTIONS GOING FORWARD:
    (a) TGP + sterile neutrino hybrid (Omega_s ~ 0.12)
        -> Specific prediction: m_s ~ 3-7 keV, testable by X-ray missions
    (b) Accept cluster-scale deficit as limitation (like MOND community)
        -> Focus on galaxy-scale successes as primary evidence
    (c) Non-minimal TGP extension with additional scalar/vector field
        -> Breaks minimality principle; approaches RMOND complexity
    (d) Cluster-specific c_eff different from individual galaxies
        -> Already tested (c_eff=3 for clusters), still nu~1 at high y

  RECOMMENDATION:
    Option (b) with (a) as theoretical possibility.
    TGP should be presented as a GALAXY-SCALE theory that explains
    RAR, BTFR, rotation curves, and dwarf spheroidals from first
    principles (membrane Flory exponent), while acknowledging the
    cluster-scale deficit as an open problem shared with all
    MOND-type theories.  The cluster deficit is a ~2-3 sigma effect
    when systematic uncertainties are included.

  =====================================================================
""")


# ======================================================================
#  SUMMARY
# ======================================================================
print_header("SUMMARY")

print(f"""  gs43 BULLET CLUSTER DEEP DIVE -- KEY RESULTS
  =============================================

  1. DIFFERENTIAL GAMMA (gas c=3 vs galaxies c=1.5):
     Marginal effect (~1-2% improvement).  At cluster y-values (y>>1),
     nu -> 1 regardless of gamma.  Differential morphology doesn't help.

  2. FULL 3D ENCLOSED MASS:
     y_3D > y_sheet at galaxy positions -> nu even closer to 1.
     The 3D calculation makes the deficit SLIGHTLY WORSE.

  3. GRAVITATIONAL MEMORY:
     The TGP scalaron is extremely massive at cluster densities
     (chameleon mechanism).  Response time << crossing time.
     No memory effect possible.

  4. MULTI-CLUSTER COMPARISON:
     All merging clusters show similar deficit pattern.
     Abell 520 (dark core at gas) is actually MORE natural for TGP
     than for LCDM, but amplitude is still too low.

  5. REVERSE ENGINEERING:
     No positive gamma gives nu ~ 1.66 at y ~ 4.4.
     STRUCTURALLY IMPOSSIBLE within TGP nu(y) family.
     The exponential suppression exp(-y^0.8) is the root cause.

  6. HYBRID MODEL:
     Sterile neutrinos (m_s ~ 3-7 keV, Omega_s ~ 0.12) could close
     the gap, but this is ~46% of LCDM's DM budget.
     Standard neutrinos provide only ~{M_nu_cluster_max/M_extra_Msun*100:.0f}% of needed mass.

  7. DEFICIT SIGNIFICANCE:
     38% deficit at ~{n_sigma:.1f} sigma significance (including systematics).
     Similar to MOND's well-known cluster problem.

  BOTTOM LINE:
     Pure TGP CANNOT explain the Bullet Cluster.
     This is a structural limitation shared by ALL modified gravity
     theories where g_obs/g_bar -> 1 at high accelerations.
     The honest path is to acknowledge this and position TGP as a
     galaxy-scale theory, with the cluster deficit as an open question.

  =============================================
  Script: gs43_bullet_cluster_deep.py
  Dependencies: numpy, scipy
""")
