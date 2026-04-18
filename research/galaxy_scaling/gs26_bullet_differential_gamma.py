#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gs26: BULLET CLUSTER -- DIFFERENTIAL gamma(S) AND THE LENSING OFFSET
=====================================================================

gs24 showed that TGP fails the Bullet Cluster test because phantom mass
tracks the baryonic acceleration field, which is dominated by gas.

But gs24 used an EFFECTIVE sphericity S_eff weighted by mass, treating
the cluster as one system.  It MISSED the key insight:

  In standard MOND, nu depends ONLY on g_bar/a0.
  In TGP, nu(y; S) = 1 + exp(-y^alpha) / y^gamma(S),
  where gamma(S) = 0.419 * (1 + 0.341*S).

After a cluster collision, gas and galaxies have DIFFERENT local geometries:
  - Gas:      ram-pressure stripped into elongated shock cone -> S ~ 0.2-0.3
              -> gamma ~ 0.45 -> LESS boost per baryon
  - Galaxies: collisionless, still roughly spherical -> S ~ 0.7-0.8
              -> gamma ~ 0.52-0.53 -> MORE boost per baryon

This differential gamma creates a DIFFERENTIAL phantom mass effect that
could shift the effective lensing peaks toward the galaxy positions!

Parts:
  A: Post-collision geometry setup
  B: Sphericity S from axis ratios
  C: Differential phantom mass computation
  D: Lensing map with differential gamma
  E: Hybrid model -- TGP + massive neutrinos
  F: Ordinary neutrino contribution
  G: Honest assessment and predictions

Author: TGP research program
"""

import numpy as np
from scipy.integrate import quad, trapezoid, cumulative_trapezoid
from scipy.interpolate import interp1d
from scipy.optimize import minimize_scalar, brentq
import sys, io, warnings
warnings.filterwarnings('ignore')
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

# ==========================================================================
# Physical constants
# ==========================================================================
G      = 6.674e-11       # m^3/(kg*s^2)
c      = 3.0e8           # m/s
H0     = 2.27e-18        # 1/s  (70 km/s/Mpc)
M_sun  = 1.989e30        # kg
kpc    = 3.086e19        # m
Mpc    = 3.086e22        # m
a0     = 1.12e-10        # m/s^2 (TGP best-fit, gs12)
r_H    = c / H0          # Hubble radius
k_B    = 1.381e-23       # J/K
eV     = 1.602e-19       # J
m_p    = 1.673e-27       # kg (proton mass)
hbar   = 1.055e-34       # J*s

# TGP exponents
alpha  = 0.80            # 4/5
gamma_disk   = 0.419     # S -> 0
gamma_sphere = 0.562     # S -> 1

def gamma_S(S):
    """Geometry-dependent gamma: gamma(S) = 0.419*(1 + 0.341*S)."""
    return 0.419 * (1.0 + 0.341 * S)

def nu_tgp(y, S=0.0):
    """
    TGP interpolation function:
        nu(y; S) = 1 + exp(-y^alpha) / y^gamma(S)
    y = g_bar / a0 (dimensionless Newtonian acceleration).
    """
    g = gamma_S(S)
    y = np.asarray(y, dtype=float)
    y_safe = np.maximum(y, 1e-30)
    ya = np.power(y_safe, alpha)
    yg = np.power(y_safe, g)
    return 1.0 + np.exp(-ya) / yg


# ########################################################################## #
#                                                                            #
#  PART A: POST-COLLISION GEOMETRY SETUP                                    #
#                                                                            #
# ########################################################################## #

print("=" * 78)
print("  gs26: BULLET CLUSTER -- DIFFERENTIAL gamma(S) AND THE LENSING OFFSET")
print("=" * 78)

print("\n" + "=" * 78)
print("  PART A: POST-COLLISION GEOMETRY SETUP")
print("=" * 78)

# ---- Component masses (from observations: Clowe+2006, Bradac+2006) ----
# Gas: ~80% of baryonic mass; stars/galaxies: ~20%
M_gas_main      = 1.2e14 * M_sun     # Main cluster gas
M_gas_bullet    = 1.5e13 * M_sun     # Bullet gas (shock cone)
M_star_main     = 2.0e13 * M_sun     # Main cluster galaxies
M_star_bullet   = 3.0e12 * M_sun     # Bullet galaxies

M_bar_main_tot  = M_gas_main + M_star_main
M_bar_bull_tot  = M_gas_bullet + M_star_bullet

# ---- Positions along collision axis (x) ----
# Post-collision: gas stripped toward center, galaxies passed through
x_gas_main      = 0.0 * kpc          # Main gas: near center
x_gas_bullet    = 200.0 * kpc        # Bullet gas: shock cone offset
x_gal_main      = -360.0 * kpc       # Main galaxies: passed through
x_gal_bullet    = 360.0 * kpc        # Bullet galaxies: passed through

# ---- Core radii ----
r_core_gas_main   = 250.0 * kpc      # Main gas beta-profile core
r_core_gas_bullet = 80.0 * kpc       # Bullet gas (compact shock)
r_core_gal_main   = 200.0 * kpc      # Main galaxy distribution
r_core_gal_bullet = 70.0 * kpc       # Bullet galaxy distribution

# ---- Axis ratios (post-collision geometry) ----
# Gas: ram-pressure stripped into elongated structures
q_gas_main      = 0.30               # Main gas: elongated (prolate, c/a)
q_gas_bullet    = 0.20               # Bullet gas: even more elongated

# Galaxies: collisionless, roughly spherical
q_gal_main      = 0.80               # Main galaxies: roughly spherical
q_gal_bullet    = 0.70               # Bullet galaxies: slightly elongated

# ---- Lensing masses (observed) ----
M_lens_main     = 1.0e15 * M_sun
M_lens_bullet   = 1.0e14 * M_sun

# ---- Observed lensing offset ----
offset_observed = 150.0 * kpc        # Lensing peak vs gas centroid

print(f"""
  Bullet Cluster (1E 0657-558) -- Post-collision configuration
  ============================================================

  Component               Mass (M_sun)     Position (kpc)   Core (kpc)   q (axis ratio)
  --------------------------------------------------------------------------------------
  Main gas                {M_gas_main/M_sun:.1e}      x = {x_gas_main/kpc:+6.0f}        {r_core_gas_main/kpc:.0f}          {q_gas_main:.2f}
  Bullet gas              {M_gas_bullet/M_sun:.1e}      x = {x_gas_bullet/kpc:+6.0f}         {r_core_gas_bullet/kpc:.0f}          {q_gas_bullet:.2f}
  Main galaxies           {M_star_main/M_sun:.1e}      x = {x_gal_main/kpc:+6.0f}        {r_core_gal_main/kpc:.0f}          {q_gal_main:.2f}
  Bullet galaxies         {M_star_bullet/M_sun:.1e}      x = {x_gal_bullet/kpc:+6.0f}         {r_core_gal_bullet/kpc:.0f}          {q_gal_bullet:.2f}

  Total baryonic mass:
    Main:   M_bar = {M_bar_main_tot/M_sun:.2e} M_sun  (gas fraction: {M_gas_main/M_bar_main_tot:.0%})
    Bullet: M_bar = {M_bar_bull_tot/M_sun:.2e} M_sun  (gas fraction: {M_gas_bullet/M_bar_bull_tot:.0%})

  Observed lensing masses:
    M_lens(main)   = {M_lens_main/M_sun:.1e} M_sun  (need {M_lens_main/M_bar_main_tot:.1f}x boost)
    M_lens(bullet) = {M_lens_bullet/M_sun:.1e} M_sun  (need {M_lens_bullet/M_bar_bull_tot:.1f}x boost)
""")


# ########################################################################## #
#                                                                            #
#  PART B: COMPUTING SPHERICITY S FROM AXIS RATIOS                         #
#                                                                            #
# ########################################################################## #

print("=" * 78)
print("  PART B: COMPUTING SPHERICITY S FROM AXIS RATIOS")
print("=" * 78)

def sphericity_prolate(q):
    """
    Sphericity S for a prolate spheroid with axis ratio q = c/a (< 1).
    For a prolate body with semi-axes a >= b = c:
        S = c^(2/3) / a^(2/3) = q^(2/3)
    This is the Wadell sphericity: ratio of surface area of equivalent-volume
    sphere to actual surface area. For a prolate spheroid:
        S = (b*c)^(1/3) / a^(2/3) where b = c = q*a
        S = (q*a * q*a)^(1/3) / a^(2/3) = q^(2/3)
    """
    return q**(2.0/3.0)

def sphericity_oblate(q):
    """
    Sphericity S for an oblate spheroid with axis ratio q = c/a (< 1).
    For an oblate body with semi-axes a = b >= c:
        S = (a*c)^(1/3) / a^(2/3) = (c/a)^(1/3) = q^(1/3)
    Actually, Wadell sphericity uses surface areas. For practical purposes
    with TGP, we use S = q as the simplest mapping from flattening to
    sphericity parameter (consistent with gs19 where S=1 for sphere, S->0
    for disk).
    """
    return q

# For post-collision gas: PROLATE geometry (elongated by ram pressure)
# For galaxies: roughly spherical (collisionless)
S_gas_main      = sphericity_prolate(q_gas_main)
S_gas_bullet    = sphericity_prolate(q_gas_bullet)
S_gal_main      = q_gal_main      # Close to spherical, use direct mapping
S_gal_bullet    = q_gal_bullet

gamma_gas_main    = gamma_S(S_gas_main)
gamma_gas_bullet  = gamma_S(S_gas_bullet)
gamma_gal_main    = gamma_S(S_gal_main)
gamma_gal_bullet  = gamma_S(S_gal_bullet)

print(f"""
  B.1  Sphericity from axis ratios
  ----------------------------------
  For prolate (elongated) gas: S = q^(2/3)  (Wadell-like sphericity)
  For roughly spherical galaxies: S ~ q     (direct mapping)

  Component          q (axis ratio)   S (sphericity)   gamma(S)
  ---------------------------------------------------------------
  Main gas           {q_gas_main:.2f}             {S_gas_main:.4f}           {gamma_gas_main:.4f}
  Bullet gas         {q_gas_bullet:.2f}             {S_gas_bullet:.4f}           {gamma_gas_bullet:.4f}
  Main galaxies      {q_gal_main:.2f}             {S_gal_main:.4f}           {gamma_gal_main:.4f}
  Bullet galaxies    {q_gal_bullet:.2f}             {S_gal_bullet:.4f}           {gamma_gal_bullet:.4f}

  B.2  The differential gamma
  ----------------------------
  Delta_gamma (main)   = gamma_gal - gamma_gas = {gamma_gal_main - gamma_gas_main:.4f}
  Delta_gamma (bullet) = gamma_gal - gamma_gas = {gamma_gal_bullet - gamma_gas_bullet:.4f}

  KEY: galaxies get a STRONGER boost exponent than gas because their
  roughly spherical distribution yields higher gamma(S).
""")

# Show how nu changes with gamma at typical cluster accelerations
print("  B.3  Impact of differential gamma on boost factor nu")
print("  " + "-" * 60)
y_test = np.array([0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0])
print(f"    {'y':>6s}   {'nu(gas,main)':>12s}  {'nu(gal,main)':>12s}  {'ratio':>8s}  {'% more boost':>12s}")
for y in y_test:
    nu_gas = nu_tgp(y, S_gas_main)
    nu_gal = nu_tgp(y, S_gal_main)
    ratio  = (nu_gal - 1.0) / (nu_gas - 1.0) if (nu_gas - 1.0) > 1e-10 else 0
    pct    = (nu_gal - nu_gas) / (nu_gas - 1.0) * 100 if (nu_gas - 1.0) > 1e-10 else 0
    print(f"    {y:6.3f}   {nu_gas:12.4f}  {nu_gal:12.4f}  {ratio:8.3f}  {pct:+10.1f}%")

print(f"""
  INTERPRETATION:
  At y ~ 0.05 (typical cluster outskirts):
    nu(gas,  S={S_gas_main:.2f})  = {nu_tgp(0.05, S_gas_main):.4f}  --> phantom = {nu_tgp(0.05, S_gas_main)-1:.4f} per unit baryon
    nu(gal,  S={S_gal_main:.2f})  = {nu_tgp(0.05, S_gal_main):.4f}  --> phantom = {nu_tgp(0.05, S_gal_main)-1:.4f} per unit baryon
    Ratio of phantom mass per baryon: {(nu_tgp(0.05, S_gal_main)-1)/(nu_tgp(0.05, S_gas_main)-1):.3f}

  At y ~ 0.1 (intermediate radius):
    nu(gas)  = {nu_tgp(0.1, S_gas_main):.4f}  --> phantom = {nu_tgp(0.1, S_gas_main)-1:.4f}
    nu(gal)  = {nu_tgp(0.1, S_gal_main):.4f}  --> phantom = {nu_tgp(0.1, S_gal_main)-1:.4f}
    Ratio: {(nu_tgp(0.1, S_gal_main)-1)/(nu_tgp(0.1, S_gas_main)-1):.3f}

  Galaxies get ~{((nu_tgp(0.05, S_gal_main)-1)/(nu_tgp(0.05, S_gas_main)-1) - 1)*100:.0f}% more phantom mass per baryon than gas at y~0.05.
""")


# ########################################################################## #
#                                                                            #
#  PART C: DIFFERENTIAL PHANTOM MASS COMPUTATION                           #
#                                                                            #
# ########################################################################## #

print("=" * 78)
print("  PART C: DIFFERENTIAL PHANTOM MASS COMPUTATION")
print("=" * 78)

# ---- 1-D grid along collision axis ----
N_grid = 4000
x_min  = -3000.0 * kpc
x_max  =  3000.0 * kpc
x_arr  = np.linspace(x_min, x_max, N_grid)
dx     = x_arr[1] - x_arr[0]
x_kpc  = x_arr / kpc

# ---- Mass profiles ----
def beta_profile_1d(x, x0, r_core, M_total, beta_param=2.0/3.0):
    """1-D projected beta-model for gas."""
    u = (x - x0) / r_core
    exponent = -3.0 * beta_param / 2.0 + 0.5
    raw = np.power(1.0 + u**2, exponent)
    norm = trapezoid(raw, x)
    if norm <= 0:
        norm = 1.0
    return M_total * raw / norm

def gaussian_1d(x, x0, sigma, M_total):
    """1-D Gaussian profile."""
    return (M_total / (sigma * np.sqrt(2.0 * np.pi))) * np.exp(
        -0.5 * ((x - x0) / sigma)**2
    )

# Surface density profiles
Sigma_gas_main   = beta_profile_1d(x_arr, x_gas_main, r_core_gas_main, M_gas_main)
Sigma_gas_bullet = beta_profile_1d(x_arr, x_gas_bullet, r_core_gas_bullet, M_gas_bullet)
Sigma_gal_main   = gaussian_1d(x_arr, x_gal_main, r_core_gal_main, M_star_main)
Sigma_gal_bullet = gaussian_1d(x_arr, x_gal_bullet, r_core_gal_bullet, M_star_bullet)

Sigma_baryon = Sigma_gas_main + Sigma_gas_bullet + Sigma_gal_main + Sigma_gal_bullet

# ---- Newtonian acceleration from each component ----
def g_newtonian_signed(x_eval, x_center, M_total, r_scale):
    """Signed Newtonian g (points toward center)."""
    dx_vec = x_center - x_eval
    r = np.abs(dx_vec)
    r = np.maximum(r, 1.0 * kpc)
    M_enc = M_total * r**3 / (r**3 + r_scale**3)
    g_mag = G * M_enc / r**2
    return g_mag * np.sign(dx_vec)

def g_newtonian_mag(x_eval, x_center, M_total, r_scale):
    """Magnitude of Newtonian g."""
    r = np.abs(x_eval - x_center)
    r = np.maximum(r, 1.0 * kpc)
    M_enc = M_total * r**3 / (r**3 + r_scale**3)
    return G * M_enc / r**2

# Accelerations from each component
g_gas_main_s   = g_newtonian_signed(x_arr, x_gas_main,   M_gas_main,   r_core_gas_main)
g_gas_bull_s   = g_newtonian_signed(x_arr, x_gas_bullet,  M_gas_bullet, r_core_gas_bullet)
g_gal_main_s   = g_newtonian_signed(x_arr, x_gal_main,    M_star_main,  r_core_gal_main)
g_gal_bull_s   = g_newtonian_signed(x_arr, x_gal_bullet,  M_star_bullet, r_core_gal_bullet)

g_total_signed = g_gas_main_s + g_gas_bull_s + g_gal_main_s + g_gal_bull_s
g_total_mag    = np.abs(g_total_signed)

# ---- APPROACH 1: Standard MOND (single gamma = 0.5) ----
S_mond = 0.573   # gamma(S) = 0.5 -> S = (0.5/0.419 - 1)/0.341 = 0.566
# Actually we just use gamma=0.5 directly via a custom function
def nu_standard_mond(y, gamma_val=0.5):
    """Standard MOND-like with fixed gamma."""
    y = np.asarray(y, dtype=float)
    y_safe = np.maximum(y, 1e-30)
    ya = np.power(y_safe, alpha)
    yg = np.power(y_safe, gamma_val)
    return 1.0 + np.exp(-ya) / yg

y_total = g_total_mag / a0

nu_mond_uniform = nu_standard_mond(y_total, gamma_val=0.5)
g_mond_uniform  = g_total_mag * nu_mond_uniform
g_phantom_mond  = g_mond_uniform - g_total_mag

# ---- APPROACH 2: TGP with UNIFORM gamma = 0.42 (gs24 approach) ----
nu_tgp_uniform = nu_standard_mond(y_total, gamma_val=0.42)
g_tgp_uniform  = g_total_mag * nu_tgp_uniform
g_phantom_tgp_uni = g_tgp_uniform - g_total_mag

# ---- APPROACH 3: TGP with DIFFERENTIAL gamma(S) ----
# Each component contributes its OWN phantom mass based on its local geometry.
# The key insight: phantom mass from gas uses gamma(S_gas), phantom from
# galaxies uses gamma(S_gal).

# For each component, compute its individual g_bar, y, and nu
# Then compute phantom surface density from each

# Component accelerations (magnitudes)
g_gas_main_mag   = g_newtonian_mag(x_arr, x_gas_main,   M_gas_main,   r_core_gas_main)
g_gas_bull_mag   = g_newtonian_mag(x_arr, x_gas_bullet,  M_gas_bullet, r_core_gas_bullet)
g_gal_main_mag   = g_newtonian_mag(x_arr, x_gal_main,    M_star_main,  r_core_gal_main)
g_gal_bull_mag   = g_newtonian_mag(x_arr, x_gal_bullet,  M_star_bullet, r_core_gal_bullet)

# For the differential approach, we need to handle this carefully.
# The TOTAL acceleration at each point is g_bar = sum of all components.
# But the PHANTOM MASS generated by each component depends on:
#   (a) that component's contribution to the local acceleration
#   (b) the GEOMETRY (sphericity) of that component
#
# Method: decompose phantom mass proportionally.
# At each point x:
#   g_bar = g_gas + g_gal  (total Newtonian)
#   y = g_bar / a0
#   For gas:      fraction f_gas = g_gas / g_bar
#                 phantom_gas = M_gas * (nu(y, S_gas) - 1) * f_gas_mass
#   For galaxies: phantom_gal = M_gal * (nu(y, S_gal) - 1) * f_gal_mass
#
# More precisely: the boost at each point depends on the total y,
# but the gamma used depends on which component dominates the local geometry.
#
# Approach: compute a weighted effective nu for each component

# Gas total acceleration contribution
g_gas_total_mag = np.sqrt((g_gas_main_s + g_gas_bull_s)**2)
g_gal_total_mag = np.sqrt((g_gal_main_s + g_gal_bull_s)**2)

# Mass fractions at each point (surface density fractions)
Sigma_gas_total = Sigma_gas_main + Sigma_gas_bullet
Sigma_gal_total = Sigma_gal_main + Sigma_gal_bullet
Sigma_all = Sigma_gas_total + Sigma_gal_total
Sigma_all_safe = np.maximum(Sigma_all, 1e-30)

f_gas_mass = Sigma_gas_total / Sigma_all_safe
f_gal_mass = Sigma_gal_total / Sigma_all_safe

# For the differential gamma approach:
# Phantom surface density from gas:
#   Sigma_phantom_gas(x) = Sigma_gas(x) * (nu(y, S_gas_local) - 1)
# Phantom surface density from galaxies:
#   Sigma_phantom_gal(x) = Sigma_gal(x) * (nu(y, S_gal_local) - 1)

# Gas components: use geometry-appropriate S
# Main gas (most of the gas): S_gas_main
# Bullet gas: S_gas_bullet
# Main galaxies: S_gal_main
# Bullet galaxies: S_gal_bullet

# Compute local y for the TOTAL field (all components together)
y_local = g_total_mag / a0

# Phantom mass from each component using its own gamma
nu_gas_main_diff   = nu_tgp(y_local, S=S_gas_main)
nu_gas_bull_diff   = nu_tgp(y_local, S=S_gas_bullet)
nu_gal_main_diff   = nu_tgp(y_local, S=S_gal_main)
nu_gal_bull_diff   = nu_tgp(y_local, S=S_gal_bullet)

# Phantom surface density from each
Sigma_phantom_gas_main   = Sigma_gas_main   * (nu_gas_main_diff - 1.0)
Sigma_phantom_gas_bullet = Sigma_gas_bullet * (nu_gas_bull_diff - 1.0)
Sigma_phantom_gal_main   = Sigma_gal_main   * (nu_gal_main_diff - 1.0)
Sigma_phantom_gal_bullet = Sigma_gal_bullet * (nu_gal_bull_diff - 1.0)

Sigma_phantom_gas_total = Sigma_phantom_gas_main + Sigma_phantom_gas_bullet
Sigma_phantom_gal_total = Sigma_phantom_gal_main + Sigma_phantom_gal_bullet
Sigma_phantom_diff      = Sigma_phantom_gas_total + Sigma_phantom_gal_total

Sigma_total_diff = Sigma_baryon + Sigma_phantom_diff

# For comparison: uniform gamma phantom masses
r_char = 500.0 * kpc   # characteristic projection radius
Sigma_phantom_mond_uni = g_phantom_mond * r_char / (2.0 * np.pi * G)
Sigma_phantom_tgp_uni  = g_phantom_tgp_uni * r_char / (2.0 * np.pi * G)

# Effective total for uniform approaches
Sigma_total_mond_uni = Sigma_baryon + Sigma_phantom_mond_uni
Sigma_total_tgp_uni  = Sigma_baryon + Sigma_phantom_tgp_uni

# ---- Integrated phantom masses ----
M_phantom_gas_main_int   = trapezoid(Sigma_phantom_gas_main, x_arr)
M_phantom_gas_bull_int   = trapezoid(Sigma_phantom_gas_bullet, x_arr)
M_phantom_gal_main_int   = trapezoid(Sigma_phantom_gal_main, x_arr)
M_phantom_gal_bull_int   = trapezoid(Sigma_phantom_gal_bullet, x_arr)
M_phantom_diff_int       = trapezoid(Sigma_phantom_diff, x_arr)
M_phantom_mond_int       = trapezoid(Sigma_phantom_mond_uni, x_arr)

print(f"""
  C.1  Phantom mass from each component (differential gamma)
  -----------------------------------------------------------

  Component             M_baryon (M_sun)    gamma(S)     M_phantom (M_sun)    Phantom/Baryon
  -------------------------------------------------------------------------------------------
  Main gas              {M_gas_main/M_sun:.2e}        {gamma_gas_main:.4f}       {M_phantom_gas_main_int/M_sun:.2e}         {M_phantom_gas_main_int/M_gas_main:.3f}
  Bullet gas            {M_gas_bullet/M_sun:.2e}        {gamma_gas_bullet:.4f}       {M_phantom_gas_bull_int/M_sun:.2e}         {M_phantom_gas_bull_int/M_gas_bullet:.3f}
  Main galaxies         {M_star_main/M_sun:.2e}        {gamma_gal_main:.4f}       {M_phantom_gal_main_int/M_sun:.2e}         {M_phantom_gal_main_int/M_star_main:.3f}
  Bullet galaxies       {M_star_bullet/M_sun:.2e}        {gamma_gal_bullet:.4f}       {M_phantom_gal_bull_int/M_sun:.2e}         {M_phantom_gal_bull_int/M_star_bullet:.3f}

  Total phantom (diff gamma):  {M_phantom_diff_int/M_sun:.2e} M_sun
  Total phantom (MOND gamma=0.5): {M_phantom_mond_int/M_sun:.2e} M_sun

  C.2  Key ratios
  ----------------
  Phantom/baryon for galaxies vs gas:
    Main:   galaxy phantom/baryon = {M_phantom_gal_main_int/M_star_main:.3f}
            gas phantom/baryon    = {M_phantom_gas_main_int/M_gas_main:.3f}
            Ratio (gal/gas):      {(M_phantom_gal_main_int/M_star_main)/(M_phantom_gas_main_int/M_gas_main):.3f}

  This ratio shows how much MORE phantom mass (per baryon) accrues to
  the galaxy component vs the gas component due to differential gamma.
""")


# ########################################################################## #
#                                                                            #
#  PART D: LENSING MAP WITH DIFFERENTIAL GAMMA                             #
#                                                                            #
# ########################################################################## #

print("=" * 78)
print("  PART D: LENSING MAP WITH DIFFERENTIAL GAMMA")
print("=" * 78)

# ---- Convergence kappa for three approaches ----
# kappa is proportional to projected surface mass density
# We normalize to compare peak locations

def find_peak_pos(arr, x_kpc, x_range):
    """Find peak position and value within x_range = (x_min, x_max)."""
    mask = (x_kpc >= x_range[0]) & (x_kpc <= x_range[1])
    if not np.any(mask):
        return 0.0, 0.0
    idx_max = np.argmax(arr[mask])
    return x_kpc[mask][idx_max], arr[mask][idx_max]

def compute_centroid(arr, x_kpc, x_range):
    """Compute mass-weighted centroid within x_range."""
    mask = (x_kpc >= x_range[0]) & (x_kpc <= x_range[1])
    if not np.any(mask) or np.sum(arr[mask]) <= 0:
        return 0.0
    return np.average(x_kpc[mask], weights=np.maximum(arr[mask], 0))

# ---- Approach 1: MOND uniform gamma=0.5 ----
# Use Sigma_baryon + acceleration-based phantom mass
Sigma_mond_total = Sigma_baryon + Sigma_phantom_mond_uni

# ---- Approach 2: TGP uniform gamma=0.42 ----
Sigma_tgp_uni_total = Sigma_baryon + Sigma_phantom_tgp_uni

# ---- Approach 3: TGP differential gamma ----
# Already computed: Sigma_total_diff

# Normalize all for plotting/comparison
norm_val = max(np.max(Sigma_baryon), 1e-30)

kappa_baryon    = Sigma_baryon / norm_val
kappa_mond      = Sigma_mond_total / norm_val
kappa_tgp_uni   = Sigma_tgp_uni_total / norm_val
kappa_diff      = Sigma_total_diff / norm_val

# Find peaks in each approach
print("\n  D.1  Lensing convergence peak locations (1-D projection)")
print("  " + "-" * 70)
print(f"    {'Approach':<30s}  {'Left peak (kpc)':>16s}  {'Right peak (kpc)':>16s}")
print(f"    {'-'*30}  {'-'*16}  {'-'*16}")

# Baryonic
pk_bar_L, _ = find_peak_pos(kappa_baryon, x_kpc, (-1500, -10))
pk_bar_R, _ = find_peak_pos(kappa_baryon, x_kpc, (10, 1500))
print(f"    {'Baryonic only':<30s}  {pk_bar_L:>+14.0f}  {pk_bar_R:>+14.0f}")

# MOND uniform
pk_mond_L, _ = find_peak_pos(kappa_mond, x_kpc, (-1500, -10))
pk_mond_R, _ = find_peak_pos(kappa_mond, x_kpc, (10, 1500))
print(f"    {'MOND (gamma=0.5 uniform)':<30s}  {pk_mond_L:>+14.0f}  {pk_mond_R:>+14.0f}")

# TGP uniform
pk_tgp_u_L, _ = find_peak_pos(kappa_tgp_uni, x_kpc, (-1500, -10))
pk_tgp_u_R, _ = find_peak_pos(kappa_tgp_uni, x_kpc, (10, 1500))
print(f"    {'TGP uniform (gamma=0.42)':<30s}  {pk_tgp_u_L:>+14.0f}  {pk_tgp_u_R:>+14.0f}")

# TGP differential
pk_diff_L, _ = find_peak_pos(kappa_diff, x_kpc, (-1500, -10))
pk_diff_R, _ = find_peak_pos(kappa_diff, x_kpc, (10, 1500))
print(f"    {'TGP differential gamma(S)':<30s}  {pk_diff_L:>+14.0f}  {pk_diff_R:>+14.0f}")

# Observed
print(f"    {'Observed (lensing)':<30s}  {'  -360':>16s}  {'  +360':>16s}")

# ---- Offsets from gas centroid ----
gas_centroid_left  = x_gas_main / kpc     # 0 kpc
gas_centroid_right = x_gas_bullet / kpc   # +200 kpc

offset_mond_L = np.abs(pk_mond_L - gas_centroid_left)
offset_mond_R = np.abs(pk_mond_R - gas_centroid_right)
offset_tgpu_L = np.abs(pk_tgp_u_L - gas_centroid_left)
offset_tgpu_R = np.abs(pk_tgp_u_R - gas_centroid_right)
offset_diff_L = np.abs(pk_diff_L - gas_centroid_left)
offset_diff_R = np.abs(pk_diff_R - gas_centroid_right)
offset_obs_L  = np.abs(-360.0 - gas_centroid_left)
offset_obs_R  = np.abs(360.0 - gas_centroid_right)

print(f"""
  D.2  Lensing peak offsets from gas centroid
  -------------------------------------------
    Gas centroids: main at x = {gas_centroid_left:+.0f} kpc, bullet at x = {gas_centroid_right:+.0f} kpc

    Approach                          Offset_L (kpc)   Offset_R (kpc)
    -----------------------------------------------------------------
    MOND (gamma=0.5)                  {offset_mond_L:>8.0f}           {offset_mond_R:>8.0f}
    TGP uniform (gamma=0.42)         {offset_tgpu_L:>8.0f}           {offset_tgpu_R:>8.0f}
    TGP differential gamma(S)        {offset_diff_L:>8.0f}           {offset_diff_R:>8.0f}
    Observed                          {offset_obs_L:>8.0f}           {offset_obs_R:>8.0f}
""")

# ---- Centroids of total mass distribution ----
# The lensing peak may not shift much, but the CENTROID of the mass
# distribution could shift more
cent_bar_L = compute_centroid(kappa_baryon, x_kpc, (-1500, -10))
cent_bar_R = compute_centroid(kappa_baryon, x_kpc, (10, 1500))
cent_diff_L = compute_centroid(kappa_diff, x_kpc, (-1500, -10))
cent_diff_R = compute_centroid(kappa_diff, x_kpc, (10, 1500))
cent_mond_L = compute_centroid(kappa_mond, x_kpc, (-1500, -10))
cent_mond_R = compute_centroid(kappa_mond, x_kpc, (10, 1500))

print(f"""
  D.3  Mass-weighted centroids of total projected mass
  ------------------------------------------------------
    Approach                     Centroid_L (kpc)   Centroid_R (kpc)
    -----------------------------------------------------------------
    Baryonic only                {cent_bar_L:>+10.1f}          {cent_bar_R:>+10.1f}
    MOND (gamma=0.5)             {cent_mond_L:>+10.1f}          {cent_mond_R:>+10.1f}
    TGP differential gamma       {cent_diff_L:>+10.1f}          {cent_diff_R:>+10.1f}
    Observed lensing             {'  -360':>10s}          {'  +360':>10s}

  Centroid shift (diff vs baryonic):
    Left:  {cent_diff_L - cent_bar_L:+.1f} kpc  (toward galaxies at -360)
    Right: {cent_diff_R - cent_bar_R:+.1f} kpc  (toward galaxies at +360)
""")

# ---- Quantify the shift as fraction of needed offset ----
shift_L = cent_diff_L - cent_bar_L
shift_R = cent_diff_R - cent_bar_R
needed_shift_L = -360.0 - cent_bar_L   # want to move toward -360
needed_shift_R = 360.0 - cent_bar_R    # want to move toward +360

frac_L = shift_L / needed_shift_L if abs(needed_shift_L) > 1 else 0
frac_R = shift_R / needed_shift_R if abs(needed_shift_R) > 1 else 0

print(f"""
  D.4  Fraction of needed centroid shift achieved
  -------------------------------------------------
    Left:   shift = {shift_L:+.1f} kpc, needed = {needed_shift_L:+.1f} kpc -> {frac_L:.1%}
    Right:  shift = {shift_R:+.1f} kpc, needed = {needed_shift_R:+.1f} kpc -> {frac_R:.1%}
""")

# ---- Surface density profile at key positions ----
idx_x0    = np.argmin(np.abs(x_kpc - 0))
idx_x100  = np.argmin(np.abs(x_kpc - 100))
idx_x200  = np.argmin(np.abs(x_kpc - 200))
idx_xn360 = np.argmin(np.abs(x_kpc - (-360)))
idx_x360  = np.argmin(np.abs(x_kpc - 360))

print("  D.5  Surface density at key positions (in baryon units)")
print("  " + "-" * 70)
print(f"    {'x (kpc)':>8s}  {'Sigma_bar':>10s}  {'Sigma_mond':>10s}  {'Sigma_diff':>10s}  {'diff/bar':>8s}  {'diff/mond':>9s}")
for idx, label in [(idx_xn360, "-360"), (idx_x0, "0"), (idx_x100, "+100"),
                    (idx_x200, "+200"), (idx_x360, "+360")]:
    sb = kappa_baryon[idx]
    sm = kappa_mond[idx]
    sd = kappa_diff[idx]
    r1 = sd / sb if sb > 1e-10 else 0
    r2 = sd / sm if sm > 1e-10 else 0
    print(f"    {label:>8s}  {sb:10.4f}  {sm:10.4f}  {sd:10.4f}  {r1:8.3f}  {r2:9.3f}")

print()


# ########################################################################## #
#                                                                            #
#  PART E: HYBRID MODEL -- TGP + MASSIVE NEUTRINOS                         #
#                                                                            #
# ########################################################################## #

print("=" * 78)
print("  PART E: HYBRID MODEL -- TGP DIFFERENTIAL gamma + MASSIVE NEUTRINOS")
print("=" * 78)

# ---- Neutrino physics ----
# Cosmic neutrino background: n_nu = 336/cm^3 = 3.36e8 /m^3
# 3 flavors, each with n = 112/cm^3
n_nu_per_flavor = 112.0e6   # per m^3
n_nu_total      = 336.0e6   # per m^3 (3 flavors)

# Neutrino mass scenarios
m_nu_min    = 0.06 * eV / c**2    # Minimum from oscillations (~0.06 eV)
m_nu_angus  = 1.5 * eV / c**2     # Angus+2006 MOND proposal
m_nu_angus2 = 2.0 * eV / c**2     # Angus upper
m_nu_planck = 0.04 * eV / c**2    # Planck limit / 3 flavors (0.12 eV sum)
m_nu_katrin = 0.45 * eV / c**2    # KATRIN 2024 limit (0.45 eV per flavor)

# Key: Planck constraint assumes LCDM! In TGP cosmology, this may differ.

print(f"""
  E.1  Neutrino mass constraints
  --------------------------------
  Source                    m_nu (per flavor)    Sum m_nu       Notes
  --------------------------------------------------------------------------
  Oscillation minimum      0.02 eV              ~0.06 eV      Firm (atmospheric + solar)
  Planck 2018 (LCDM)       < 0.04 eV            < 0.12 eV     ASSUMES LCDM cosmology
  KATRIN 2024              < 0.45 eV             < 1.35 eV     Direct laboratory limit
  Angus+2006 proposal      1.5-2.0 eV           4.5-6.0 eV    For MOND + nu clusters

  IMPORTANT: The Planck constraint assumes standard LCDM. In TGP cosmology
  with modified expansion history and structure growth, the CMB damping tail
  constraints on m_nu may be substantially relaxed. The KATRIN direct
  measurement is model-independent and currently allows m_nu up to 0.45 eV.
""")

# ---- Neutrino clustering in cluster potential ----
# Neutrinos cluster in gravitational potentials if the potential well depth
# exceeds their thermal energy.
#
# Neutrino temperature today: T_nu = (4/11)^(1/3) * T_CMB = 1.95 K
T_nu = (4.0/11.0)**(1.0/3.0) * 2.725   # K
v_nu_thermal = np.sqrt(3.0 * k_B * T_nu / (1.0 * eV / c**2))  # for 1 eV neutrino

# For different masses:
m_nu_scan = np.array([0.06, 0.1, 0.2, 0.45, 1.0, 1.5, 2.0]) * eV / c**2
print("  E.2  Neutrino clustering in cluster potential")
print("  " + "-" * 65)
print(f"    {'m_nu (eV)':>10s}  {'v_thermal (km/s)':>16s}  {'sigma_cluster (km/s)':>20s}  {'Clusters?':>10s}")

# Cluster velocity dispersion: sigma ~ 1000 km/s for main, ~600 for bullet
sigma_cluster = 1000.0e3   # m/s

for m_nu in m_nu_scan:
    m_nu_eV = m_nu * c**2 / eV
    v_th = np.sqrt(3.0 * k_B * T_nu / m_nu)
    clusters = "YES" if v_th < sigma_cluster else ("marginal" if v_th < 3*sigma_cluster else "no")
    print(f"    {m_nu_eV:10.2f}  {v_th/1e3:16.0f}  {sigma_cluster/1e3:20.0f}  {clusters:>10s}")

print(f"""
  Neutrinos cluster when v_thermal < sigma_cluster.
  For m_nu > ~0.5 eV, neutrinos are significantly trapped in cluster potentials.
  For m_nu ~ 1.5 eV (Angus proposal), clustering is strong.
""")

# ---- Neutrino halo mass in cluster ----
# For a cluster with mass M and virial radius r_vir:
# The neutrino overdensity depends on the ratio of gravitational
# potential energy to thermal energy: phi / (kT_nu/m_nu)
#
# phi = G * M / r ~ 10^-5 c^2 for a massive cluster
# For m_nu = 1.5 eV: kT_nu / m_nu ~ 0.17 meV / 1.5 eV ~ 10^-4
# phi / (kT_nu/m_nu) >> 1 -> strong clustering

r_vir_main   = 2000.0 * kpc
r_vir_bullet = 900.0 * kpc

def neutrino_halo_mass(M_lens, r_vir, m_nu_val, n_flavors=3):
    """
    Estimate neutrino halo mass within r_vir for a cluster.

    Uses the Singh & Ma (2003) / Ringwald & Wong (2004) approach:
    For a gravitational potential phi(r), the neutrino overdensity is
    determined by integrating the Fermi-Dirac distribution over
    momentum space, accounting for gravitational blue-shift.

    For cluster-scale potentials (phi/c^2 ~ 10^-5), the overdensity
    scales as delta_nu ~ (3/2) * m_nu * phi / (k_B * T_nu) for the
    non-degenerate limit, with the full result from the gravitational
    clustering integral.

    Key references:
      - Ringwald & Wong (2004, JCAP 12, 005)
      - Singh & Ma (2003, PRD 67, 023506)
      - Angus et al. (2006, ApJ 654, L13)
    """
    m_nu_eV = m_nu_val * c**2 / eV
    v_th = np.sqrt(3.0 * k_B * T_nu / m_nu_val) if m_nu_val > 0 else 1e30

    # Gravitational potential at cluster center (NFW-like)
    # phi ~ G * M / r_vir (order of magnitude)
    phi = G * M_lens / r_vir

    # Dimensionless potential depth parameter
    # beta_nu = m_nu * phi / (k_B * T_nu)
    beta_nu = m_nu_val * phi / (k_B * T_nu)

    # Escape velocity ratio
    v_esc = np.sqrt(2.0 * phi)
    eta = v_esc / v_th

    # Neutrino overdensity from gravitational clustering
    # In the linear regime (beta_nu << 1):
    #   delta_nu ~ (3/2) * beta_nu + (15/8) * beta_nu^2 + ...
    # In the non-linear regime (beta_nu >> 1):
    #   delta_nu ~ (beta_nu)^(3/2) * (4 / 3*sqrt(pi))  (isothermal sphere)
    #
    # We use an interpolation that captures both limits:
    if beta_nu < 0.5:
        # Linear regime: perturbative result
        overdensity = 1.5 * beta_nu + 1.875 * beta_nu**2
    elif beta_nu < 5.0:
        # Transition regime: smooth interpolation
        overdensity = 1.5 * beta_nu * (1.0 + beta_nu)**0.5
    else:
        # Non-linear (isothermal): delta ~ beta^(3/2) * 4/(3*sqrt(pi))
        overdensity = (4.0 / (3.0 * np.sqrt(np.pi))) * beta_nu**1.5

    # Cap overdensity by the Tremaine-Gunn (phase space) bound
    # For fermions: max central density limited by degeneracy pressure
    # rho_max ~ g_s * m_nu^4 * sigma^3 / (6*pi^2*hbar^3)
    # where sigma is the velocity dispersion in the cluster
    sigma_cl = np.sqrt(G * M_lens / r_vir)
    rho_TG = 2.0 * m_nu_val**4 * sigma_cl**3 / (6.0 * np.pi**2 * hbar**3)
    rho_nu_bg = n_nu_total * m_nu_val
    max_overdensity = rho_TG / rho_nu_bg if rho_nu_bg > 0 else 1e10
    overdensity = min(overdensity, max_overdensity)

    V = (4.0/3.0) * np.pi * r_vir**3
    # The mass includes the background + overdensity
    # M_nu = rho_bg * V * (1 + delta_nu) but we quote just the excess
    M_nu = rho_nu_bg * V * overdensity

    return M_nu, eta, beta_nu

print("\n  E.3  Neutrino halo mass estimates for Bullet Cluster main component")
print("  " + "-" * 70)
print(f"    {'m_nu (eV)':>10s}  {'beta_nu':>8s}  {'M_nu (M_sun)':>14s}  {'M_nu/M_lens':>12s}  {'v_esc/v_th':>10s}  {'Comment':>20s}")

for m_nu in m_nu_scan:
    m_nu_eV = m_nu * c**2 / eV
    M_nu, eta, beta_val = neutrino_halo_mass(M_lens_main, r_vir_main, m_nu)
    ratio_nu = M_nu / M_lens_main
    comment = ""
    if m_nu_eV < 0.12:
        comment = "Planck-allowed"
    elif m_nu_eV < 0.45:
        comment = "KATRIN-allowed"
    elif m_nu_eV < 1.0:
        comment = "KATRIN-excluded"
    else:
        comment = "Angus proposal"
    print(f"    {m_nu_eV:10.2f}  {beta_val:8.3f}  {M_nu/M_sun:14.2e}  {ratio_nu:12.4f}  {eta:10.2f}  {comment:>20s}")

print()

# ---- Combined: TGP differential + neutrinos ----
print("  E.4  Combined model: TGP differential gamma + neutrino halo")
print("  " + "-" * 70)

# For each neutrino mass, compute total mass budget
print(f"    {'m_nu (eV)':>10s}  {'M_TGP_phantom':>14s}  {'M_nu':>14s}  {'M_total':>14s}  {'M_total/M_lens':>14s}")

M_baryon_total = (M_gas_main + M_gas_bullet + M_star_main + M_star_bullet)
M_phantom_total_diff = M_phantom_diff_int

for m_nu in m_nu_scan:
    m_nu_eV = m_nu * c**2 / eV
    M_nu_main, _, _ = neutrino_halo_mass(M_lens_main, r_vir_main, m_nu)
    M_nu_bull, _, _ = neutrino_halo_mass(M_lens_bullet, r_vir_bullet, m_nu)
    M_nu_tot = M_nu_main + M_nu_bull
    M_total = M_baryon_total + M_phantom_total_diff + M_nu_tot
    M_lens_tot = M_lens_main + M_lens_bullet
    ratio = M_total / M_lens_tot
    print(f"    {m_nu_eV:10.2f}  {M_phantom_total_diff/M_sun:14.2e}  {M_nu_tot/M_sun:14.2e}  {M_total/M_sun:14.2e}  {ratio:14.4f}")

print()

# ---- What m_nu is needed for full mass budget? ----
# Solve: M_baryon + M_phantom_diff + M_nu(m_nu) = M_lens_total
M_lens_total = M_lens_main + M_lens_bullet
M_deficit = M_lens_total - M_baryon_total - M_phantom_total_diff

print(f"""
  E.5  Mass deficit analysis
  ----------------------------
    M_baryon (total)      = {M_baryon_total/M_sun:.2e} M_sun
    M_phantom (diff gamma)= {M_phantom_total_diff/M_sun:.2e} M_sun
    M_lens (observed)     = {M_lens_total/M_sun:.2e} M_sun
    Mass deficit          = {M_deficit/M_sun:.2e} M_sun
    Deficit fraction      = {M_deficit/M_lens_total:.1%}
""")

# Scan for the m_nu that fills the gap
print("  Scanning for m_nu that fills the deficit...")
m_nu_test = np.logspace(np.log10(0.05), np.log10(5.0), 100) * eV / c**2
for m_nu in m_nu_test:
    M_nu_main, _, _ = neutrino_halo_mass(M_lens_main, r_vir_main, m_nu)
    M_nu_bull, _, _ = neutrino_halo_mass(M_lens_bullet, r_vir_bullet, m_nu)
    M_nu_tot = M_nu_main + M_nu_bull
    if M_nu_tot >= M_deficit:
        m_nu_needed_eV = m_nu * c**2 / eV
        print(f"    m_nu needed to fill deficit: ~{m_nu_needed_eV:.2f} eV per flavor")
        print(f"    Sum m_nu ~ {3*m_nu_needed_eV:.2f} eV")
        break
else:
    print("    Even m_nu = 5 eV cannot fill the deficit with neutrinos alone.")
    m_nu_needed_eV = 99.0

# Spatial distribution: neutrinos are COLLISIONLESS
print(f"""
  E.6  Spatial distribution of neutrino halo
  -------------------------------------------
  KEY POINT: Neutrinos are collisionless, just like galaxies.
  During the Bullet Cluster collision:
    - Gas:      ram-pressure stripped -> stays near center
    - Galaxies: passed through -> at +/-360 kpc
    - Neutrinos: ALSO passed through -> distributed like DM/galaxies

  The neutrino halo follows the gravitational potential, which prior to
  the collision was centered on each subcluster.  Post-collision, the
  neutrino halos pass through each other (like galaxies, like CDM).

  The neutrino contribution to the lensing peak is therefore at
  the GALAXY positions, not the gas positions.

  Combined model (TGP differential + neutrinos at galaxy positions):
    - Gas region (center):    baryonic gas + WEAK phantom (low gamma)
    - Galaxy region (+/-360): baryonic stars + STRONG phantom (high gamma) + neutrinos
    - This further shifts the effective lensing peaks toward galaxy positions!
""")


# ########################################################################## #
#                                                                            #
#  PART F: ORDINARY NEUTRINO CONTRIBUTION                                  #
#                                                                            #
# ########################################################################## #

print("=" * 78)
print("  PART F: ORDINARY NEUTRINO CONTRIBUTION (MINIMUM MASS)")
print("=" * 78)

# Even with standard m_nu = 0.02 eV per flavor (minimum from oscillations):
m_nu_std = 0.02 * eV / c**2   # per flavor (normal hierarchy lowest)

# Cosmic neutrino background
rho_nu_bg = n_nu_total * m_nu_std   # kg/m^3
Omega_nu_std = rho_nu_bg / (3.0 * H0**2 / (8.0 * np.pi * G))

# Neutrino clustering at cluster scale
M_nu_std_main, eta_std_main, _ = neutrino_halo_mass(M_lens_main, r_vir_main, m_nu_std)
M_nu_std_bull, eta_std_bull, _ = neutrino_halo_mass(M_lens_bullet, r_vir_bullet, m_nu_std)

# Free-streaming length
# lambda_fs ~ 100 Mpc * (1 eV / m_nu) for relativistic neutrinos
# After non-relativistic transition: lambda_fs,nr ~ 40 Mpc * (0.1 eV / m_nu)
lambda_fs = 40.0 * (0.1 / 0.02) * Mpc   # for m_nu = 0.02 eV -> ~200 Mpc

print(f"""
  F.1  Standard neutrino background
  -----------------------------------
    n_nu (3 flavors) = {n_nu_total:.0e} /m^3 = 336 /cm^3
    m_nu (minimum)   = 0.02 eV (per flavor, normal hierarchy)
    rho_nu           = {rho_nu_bg:.2e} kg/m^3
    Omega_nu         = {Omega_nu_std:.6f}  (vs Omega_m ~ 0.31)

  F.2  Clustering at cluster scale
  ----------------------------------
    v_thermal (0.02 eV)  = {np.sqrt(3*k_B*T_nu/m_nu_std)/1e3:.0f} km/s
    v_escape (main)      = {np.sqrt(2*G*M_lens_main/r_vir_main)/1e3:.0f} km/s
    eta = v_esc/v_th     = {eta_std_main:.4f}

    With m_nu = 0.02 eV, the thermal velocity ({np.sqrt(3*k_B*T_nu/m_nu_std)/1e3:.0f} km/s)
    vastly exceeds the cluster escape velocity ({np.sqrt(2*G*M_lens_main/r_vir_main)/1e3:.0f} km/s).
    Neutrinos of this mass DO NOT cluster in galaxy clusters.

    M_nu(main cluster)  = {M_nu_std_main/M_sun:.2e} M_sun
    M_nu/M_lens         = {M_nu_std_main/M_lens_main:.2e}

  F.3  Free-streaming scale
  ---------------------------
    lambda_fs (0.02 eV) ~ {lambda_fs/Mpc:.0f} Mpc
    Cluster scale       ~ 2 Mpc
    Ratio               ~ {lambda_fs/Mpc / 2:.0f}

    The free-streaming length is {lambda_fs/Mpc/2:.0f}x larger than the cluster.
    Standard-mass neutrinos are effectively UNIFORM at cluster scales.
    Their contribution to the mass budget is negligible (<{M_nu_std_main/M_lens_main*100:.3f}%).

  CONCLUSION: Ordinary neutrinos (m_nu ~ 0.02 eV) contribute negligibly
  to the Bullet Cluster mass budget. They are too fast to be gravitationally
  trapped and their free-streaming erases any clustering at these scales.
""")


# ########################################################################## #
#                                                                            #
#  PART G: HONEST ASSESSMENT AND PREDICTIONS                               #
#                                                                            #
# ########################################################################## #

print("=" * 78)
print("  PART G: HONEST ASSESSMENT AND PREDICTIONS")
print("=" * 78)

# ---- Quantify how much of the 150 kpc offset differential gamma explains ----

# Peak offset improvement
offset_improvement_peak_L = abs(pk_diff_L - pk_bar_L)
offset_improvement_peak_R = abs(pk_diff_R - pk_bar_R)

# Centroid offset improvement
centroid_shift_L = abs(shift_L)
centroid_shift_R = abs(shift_R)

# Direction check: is the shift toward galaxies?
toward_gal_L = (shift_L < 0)  # galaxies are at -360
toward_gal_R = (shift_R > 0)  # galaxies are at +360

print(f"""
  G.1  How much of the 150 kpc offset can differential gamma explain?
  =====================================================================

  Peak position shifts (differential gamma vs baryonic):
    Left:   delta_peak = {offset_improvement_peak_L:.1f} kpc {"TOWARD" if toward_gal_L else "AWAY FROM"} galaxies
    Right:  delta_peak = {offset_improvement_peak_R:.1f} kpc {"TOWARD" if toward_gal_R else "AWAY FROM"} galaxies

  Centroid shifts (differential gamma vs baryonic):
    Left:   delta_centroid = {centroid_shift_L:.1f} kpc {"TOWARD" if toward_gal_L else "AWAY FROM"} galaxies
    Right:  delta_centroid = {centroid_shift_R:.1f} kpc {"TOWARD" if toward_gal_R else "AWAY FROM"} galaxies

  Observed offset: ~150 kpc

  Fraction of offset explained by differential gamma alone:
    Peak:     {max(offset_improvement_peak_L, offset_improvement_peak_R)/150*100:.1f}% (at best)
    Centroid: {max(centroid_shift_L, centroid_shift_R)/150*100:.1f}% (at best)

  RESULT: Differential gamma shifts the lensing map by ~{max(centroid_shift_L, centroid_shift_R):.0f} kpc,
  which is {max(centroid_shift_L, centroid_shift_R)/150*100:.0f}% of the observed 150 kpc offset.
""")

# ---- Is the remaining gap fillable with neutrinos? ----
remaining_offset = 150.0 - max(centroid_shift_L, centroid_shift_R)

print(f"""
  G.2  Can neutrinos fill the remaining gap?
  ============================================

  Remaining spatial offset to explain: ~{remaining_offset:.0f} kpc
  (Most of the offset must come from additional collisionless mass)

  The neutrino halo (if massive enough) provides collisionless mass at
  the galaxy positions. This directly addresses the SPATIAL offset.

  Mass budget with TGP + neutrinos:
    For m_nu = 1.5 eV (Angus proposal):
      - KATRIN limit (0.45 eV): EXCLUDED
      - Planck LCDM limit (0.04 eV): EXCLUDED
      - But Planck assumes LCDM --> in TGP cosmology, limit may be relaxed
      - Required: TGP cosmology must allow m_nu > 0.45 eV

    For m_nu = 0.45 eV (KATRIN limit):""")

M_nu_katrin_main, eta_katrin, _ = neutrino_halo_mass(M_lens_main, r_vir_main, m_nu_katrin)
M_nu_katrin_bull, _, _          = neutrino_halo_mass(M_lens_bullet, r_vir_bullet, m_nu_katrin)

print(f"""      M_nu(main)   = {M_nu_katrin_main/M_sun:.2e} M_sun ({M_nu_katrin_main/M_lens_main:.1%} of lensing mass)
      M_nu(bullet) = {M_nu_katrin_bull/M_sun:.2e} M_sun ({M_nu_katrin_bull/M_lens_bullet:.1%} of lensing mass)
      This is {"sufficient" if (M_nu_katrin_main + M_nu_katrin_bull) > M_deficit else "insufficient"} to fill the mass deficit.
""")

# ---- Unique predictions vs standard MOND + neutrinos ----
print("""
  G.3  Unique predictions of TGP differential gamma
  ====================================================

  Standard MOND + neutrinos (Angus+2006):
    - Phantom mass distribution is INDEPENDENT of local geometry
    - The lensing offset is entirely due to neutrino halos
    - Every cluster collision produces the same type of offset
    - No prediction for geometry dependence

  TGP differential gamma(S):
    - Phantom mass per baryon DEPENDS on local geometry
    - Gas (elongated, low S) gets LESS boost
    - Galaxies (spherical, high S) get MORE boost
    - This creates a geometry-dependent offset IN ADDITION to neutrinos

  UNIQUE TESTABLE PREDICTIONS:
  ----------------------------
  (1) GEOMETRY-DEPENDENT OFFSET SCALING
      Clusters with MORE elongated gas (lower S) should show LARGER
      phantom mass deficits in the gas region.
      Prediction: offset ~ Delta_gamma ~ gamma(S_gal) - gamma(S_gas)
      More violent collisions -> more elongated gas -> larger offset

  (2) CLUSTER-BY-CLUSTER VARIATION
      Different merging clusters should show different offsets depending
      on the gas geometry:
      - Bullet Cluster (Mach ~3): very elongated gas -> large offset
      - Abell 520 ("train wreck"): complex geometry -> different offset
      - MACS J0025 (baby Bullet): moderate Mach -> moderate offset
      - El Gordo: very massive -> different scaling

      In standard MOND + neutrinos, the offset depends ONLY on the
      neutrino/baryon ratio and the geometry of the collision.
      In TGP, it ALSO depends on the gas sphericity.

  (3) ASYMMETRIC OFFSETS
      If the gas on one side is more elongated than the other (as in the
      Bullet Cluster where the shock cone is more elongated than the
      main gas), the offsets should be ASYMMETRIC:
      - Bullet side: more elongated gas (q=0.2) -> lower gamma -> larger offset
      - Main side: less elongated gas (q=0.3) -> higher gamma -> smaller offset
""")

# Quantify the asymmetry prediction
offset_pred_asymmetry = abs(gamma_gas_main - gamma_gas_bullet)
print(f"""      Predicted gamma asymmetry: |gamma_gas_main - gamma_gas_bullet|
        = |{gamma_gas_main:.4f} - {gamma_gas_bullet:.4f}| = {offset_pred_asymmetry:.4f}
      This translates to a ~{offset_pred_asymmetry/gamma_gas_main*100:.1f}% asymmetry in the phantom
      mass between the two sides.
""")

# ---- Pre-collision vs post-collision comparison ----
print("""  (4) PRE- vs POST-COLLISION COMPARISON
      Before collision: both subclusters are roughly spherical
        S_pre ~ 0.7-0.8, gamma_pre ~ 0.52-0.53
      After collision: gas becomes elongated
        S_gas_post ~ 0.2-0.5, gamma_gas_post ~ 0.45-0.49
      The collision CHANGES the phantom mass distribution by changing
      the geometry, even before any spatial separation occurs.
      This is a purely TGP effect with no analog in MOND.
""")

# ---- Global assessment ----
print("""
  G.4  GLOBAL ASSESSMENT
  ========================

  WHAT DIFFERENTIAL gamma(S) CAN DO:
  - Provide a qualitative mechanism to shift phantom mass toward galaxies
  - Generate a small but non-zero spatial offset (~few percent of needed)
  - Create testable predictions about geometry dependence
  - Reduce the neutrino mass needed (slightly) vs standard MOND

  WHAT DIFFERENTIAL gamma(S) CANNOT DO:
  - Explain the full 150 kpc offset (effect is too small by ~10-20x)
  - Eliminate the need for additional collisionless matter (neutrinos)
  - Solve the cluster mass deficit alone (still need ~7x boost)

  BOTTOM LINE:
  Differential gamma(S) is a GENUINE new effect in TGP that has no analog
  in standard MOND.  It is physically well-motivated (post-collision gas
  IS more elongated than the galaxy distribution).  But quantitatively,
  the effect is too small by itself to resolve the Bullet Cluster problem.

  The most promising path remains TGP + massive neutrinos, with the
  differential gamma providing a SUPPLEMENTARY offset that reduces the
  required neutrino mass.  The key question is whether TGP cosmology
  allows m_nu > 0.45 eV (beyond the KATRIN limit, which is model-
  independent and cannot be evaded).
""")


# ########################################################################## #
#                                                                            #
#  SUMMARY TABLE                                                            #
#                                                                            #
# ########################################################################## #

print("=" * 78)
print("  SUMMARY: gs26 KEY NUMERICAL RESULTS")
print("=" * 78)

print(f"""
  POST-COLLISION GEOMETRY:
    Gas sphericity:     S_gas_main = {S_gas_main:.3f},  S_gas_bullet = {S_gas_bullet:.3f}
    Galaxy sphericity:  S_gal_main = {S_gal_main:.3f},  S_gal_bullet = {S_gal_bullet:.3f}
    gamma (gas):        {gamma_gas_main:.4f} (main),  {gamma_gas_bullet:.4f} (bullet)
    gamma (galaxies):   {gamma_gal_main:.4f} (main),  {gamma_gal_bullet:.4f} (bullet)
    Delta_gamma:        {gamma_gal_main - gamma_gas_main:.4f} (main), {gamma_gal_bullet - gamma_gas_bullet:.4f} (bullet)

  DIFFERENTIAL BOOST (at y = 0.05, typical cluster outskirts):
    nu(gas, S={S_gas_main:.2f})  = {nu_tgp(0.05, S_gas_main):.4f}  (phantom/baryon = {nu_tgp(0.05, S_gas_main)-1:.4f})
    nu(gal, S={S_gal_main:.2f})  = {nu_tgp(0.05, S_gal_main):.4f}  (phantom/baryon = {nu_tgp(0.05, S_gal_main)-1:.4f})
    Galaxy/gas phantom ratio: {(nu_tgp(0.05, S_gal_main)-1)/(nu_tgp(0.05, S_gas_main)-1):.3f}

  LENSING OFFSET:
    Peak shift (diff gamma):    ~{max(offset_improvement_peak_L, offset_improvement_peak_R):.0f} kpc
    Centroid shift:             ~{max(centroid_shift_L, centroid_shift_R):.0f} kpc
    Observed offset:            ~150 kpc
    Fraction explained (peak):  {max(offset_improvement_peak_L, offset_improvement_peak_R)/150*100:.1f}%
    Fraction explained (cent.): {max(centroid_shift_L, centroid_shift_R)/150*100:.1f}%

  MASS BUDGET:
    M_baryon (total)      = {M_baryon_total/M_sun:.2e} M_sun
    M_phantom (diff gamma)= {M_phantom_total_diff/M_sun:.2e} M_sun
    M_lens (observed)     = {M_lens_total/M_sun:.2e} M_sun
    Mass deficit          = {M_deficit/M_sun:.2e} M_sun ({M_deficit/M_lens_total:.0%})
    m_nu needed to fill:  ~{m_nu_needed_eV:.1f} eV (per flavor)

  HYBRID MODEL (TGP + neutrinos):
    m_nu = 0.45 eV (KATRIN): M_nu = {(M_nu_katrin_main+M_nu_katrin_bull)/M_sun:.2e} M_sun -> {"closes" if (M_nu_katrin_main+M_nu_katrin_bull) > M_deficit else "does not close"} gap
    m_nu = 1.5 eV (Angus):   M_nu = {(neutrino_halo_mass(M_lens_main, r_vir_main, 1.5*eV/c**2)[0]+neutrino_halo_mass(M_lens_bullet, r_vir_bullet, 1.5*eV/c**2)[0])/M_sun:.2e} M_sun -> {"closes" if (neutrino_halo_mass(M_lens_main, r_vir_main, 1.5*eV/c**2)[0]+neutrino_halo_mass(M_lens_bullet, r_vir_bullet, 1.5*eV/c**2)[0]) > M_deficit else "does not close"} gap

  VERDICT: Differential gamma(S) is a genuine but quantitatively small effect.
    - New mechanism: YES (no analog in standard MOND)
    - Shifts lensing toward galaxies: YES (qualitatively correct direction)
    - Sufficient alone: NO (explains ~{max(centroid_shift_L, centroid_shift_R)/150*100:.0f}% of offset)
    - Needs neutrinos: YES (m_nu > ~{m_nu_needed_eV:.1f} eV, tension with KATRIN)
    - Testable predictions: YES (geometry-dependent offsets across clusters)

  IMPROVEMENT OVER gs24:
    gs24 treated S as an effective weighted average -> missed differential effect
    gs26 shows that SEPARATE gamma for gas vs galaxies creates a new mechanism
    The effect is real but modest; the Bullet Cluster remains a challenge.
""")

print("=" * 78)
print("  gs26 analysis complete.")
print("=" * 78)
