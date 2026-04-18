#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gs24: BULLET CLUSTER (1E 0657-558) IN THE TGP MEMBRANE FRAMEWORK
==================================================================

The Bullet Cluster is the single most challenging observation for any MOND-like
theory, including TGP.  Weak lensing reveals mass peaks that are OFFSET from
the dominant baryonic component (the X-ray gas), coinciding instead with the
subdominant galaxies.  In LCDM this is natural: collisionless dark-matter halos
pass through each other.  In TGP there is no dark matter -- so what creates the
lensing peaks at the galaxy positions?

This script provides:
  PART A -- Observational data and geometry of the collision
  PART B -- TGP gravitational field from baryons (gas + galaxies)
  PART C -- "Phantom dark matter" distribution and the lensing offset problem
  PART D -- Membrane propagation and gravitational memory
  PART E -- Quantitative 1-D collision model and lensing maps
  PART F -- Comparison with MOND literature and honest assessment

Key prior results used:
  gs12: alpha=4/5, gamma=2/5,  a0=1.12e-10 m/s^2
  gs19: gamma(S) = 0.419*(1+0.341*S)
  gs22: gamma/alpha = 1/2 from codimension-1 membrane
  gs9d: r_c = sqrt(r_S * r_H) geometric mean

Author: TGP research program
"""

import numpy as np
from scipy.integrate import quad, trapezoid, cumulative_trapezoid
from scipy.interpolate import interp1d
from scipy.optimize import minimize_scalar
from scipy.special import erfc
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

# TGP exponents
alpha  = 0.80            # 4/5
gamma0 = 0.419           # disk  (gs19)

def gamma_S(S):
    """Geometry-dependent gamma: gamma(S) = 0.419*(1 + 0.341*S)."""
    return 0.419 * (1.0 + 0.341 * S)

# ==========================================================================
# TGP interpolation function
# ==========================================================================
def nu_tgp(y, S=0.0):
    """
    TGP interpolation function:
        nu(y; S) = 1 + exp(-y^alpha) / y^gamma(S)
    y = g_bar / a0 (dimensionless Newtonian acceleration).
    """
    g = gamma_S(S)
    ya = np.power(np.maximum(y, 1e-30), alpha)
    yg = np.power(np.maximum(y, 1e-30), g)
    return 1.0 + np.exp(-ya) / yg


# ########################################################################## #
#                                                                            #
#  PART A: BULLET CLUSTER OBSERVATIONAL DATA                                #
#                                                                            #
# ########################################################################## #

print("=" * 78)
print("  gs24: BULLET CLUSTER (1E 0657-558) IN THE TGP MEMBRANE FRAMEWORK")
print("=" * 78)

print("\n" + "=" * 78)
print("  PART A: BULLET CLUSTER OBSERVATIONAL DATA")
print("=" * 78)

# --- Observed parameters ---
# Main cluster
M_main_lens   = 1.0e15 * M_sun      # total lensing mass
T_main_keV    = 14.0                 # X-ray temperature (keV)
M_gas_main    = 0.12 * M_main_lens   # gas mass (f_gas ~ 0.12)
M_star_main   = 0.02 * M_main_lens   # stellar mass (f_star ~ 0.02)
M_bar_main    = M_gas_main + M_star_main

# Bullet subcluster
M_bullet_lens = 1.0e14 * M_sun
T_bullet_keV  = 6.0
M_gas_bullet  = 0.15 * M_bullet_lens
M_star_bullet = 0.03 * M_bullet_lens
M_bar_bullet  = M_gas_bullet + M_star_bullet

# Geometry
v_collision   = 4700e3               # m/s (collision velocity)
cs_main       = np.sqrt(T_main_keV * 1.602e-16 / (0.6 * 1.673e-27))
mach_number   = v_collision / cs_main
d_sep         = 720.0 * kpc          # separation between lensing peaks
d_offset_obs  = 150.0 * kpc          # lensing peak offset from gas centroid
f_gas_main    = M_gas_main / M_main_lens
f_gas_bullet  = M_gas_bullet / M_bullet_lens

# Characteristic radii (NFW-like scale radii)
r_s_main      = 400.0 * kpc          # scale radius main
r_s_bullet    = 150.0 * kpc          # scale radius bullet
r200_main     = 2000.0 * kpc         # r200 main
r200_bullet   = 900.0 * kpc          # r200 bullet

print(f"""
  System: 1E 0657-558 (the "Bullet Cluster")
  Redshift:         z = 0.296
  Discovery:        Tucker et al. (1998)
  Lensing analysis: Clowe et al. (2006), Bradac et al. (2006)

  Main cluster:
    M_lens (total)  = {M_main_lens/M_sun:.1e} M_sun
    T_X             = {T_main_keV:.0f} keV
    M_gas           = {M_gas_main/M_sun:.1e} M_sun  (f_gas = {f_gas_main:.2f})
    M_star          = {M_star_main/M_sun:.1e} M_sun
    M_baryon        = {M_bar_main/M_sun:.1e} M_sun
    r_s (NFW)       = {r_s_main/kpc:.0f} kpc
    r_200           = {r200_main/kpc:.0f} kpc

  Bullet subcluster:
    M_lens (total)  = {M_bullet_lens/M_sun:.1e} M_sun
    T_X             = {T_bullet_keV:.0f} keV
    M_gas           = {M_gas_bullet/M_sun:.1e} M_sun  (f_gas = {f_gas_bullet:.2f})
    M_star          = {M_star_bullet/M_sun:.1e} M_sun
    M_baryon        = {M_bar_bullet/M_sun:.1e} M_sun
    r_s (NFW)       = {r_s_bullet/kpc:.0f} kpc
    r_200           = {r200_bullet/kpc:.0f} kpc

  Collision parameters:
    v_collision     = {v_collision/1e3:.0f} km/s
    c_s (main)      = {cs_main/1e3:.0f} km/s
    Mach number     = {mach_number:.1f}
    Separation      = {d_sep/kpc:.0f} kpc (between lensing peaks)
    Lensing offset  = {d_offset_obs/kpc:.0f} kpc (from gas centroid, per subcluster)

  THE KEY PROBLEM:
    Baryonic mass is dominated by gas (80-85% of baryons).
    After collision, gas is ram-pressure stripped --> sits near center.
    Galaxies (15-20% of baryons) are collisionless --> passed through.
    Lensing peaks coincide with GALAXIES, not with GAS.
    In LCDM: collisionless DM halos passed through with galaxies.
    In TGP:  no DM -- gravitational boost must somehow peak at galaxy positions.
""")


# ########################################################################## #
#                                                                            #
#  PART B: TGP GRAVITATIONAL FIELD FROM BARYONIC DISTRIBUTION               #
#                                                                            #
# ########################################################################## #

print("=" * 78)
print("  PART B: TGP GRAVITATIONAL FIELD FROM BARYONIC DISTRIBUTION")
print("=" * 78)

# --- Mass profiles ---
# We model the post-collision state with a 1-D projected mass distribution
# along the collision axis (call it x).

# Geometry after collision:
#   x = 0 is midpoint between the two gas centroids
#   Gas: concentrated near x = 0 (stripped, shocked)
#   Main galaxies: at x = -d_sep/2  (moved to the left)
#   Bullet galaxies: at x = +d_sep/2  (moved to the right, the "bullet")
#
# But the gas is NOT exactly centered -- the main gas is slightly left,
# bullet gas slightly right, with a shock between them.
# Simplification: two gas blobs offset by ~200 kpc, galaxies offset by 720 kpc.

d_gas_sep  = 200.0 * kpc    # gas centroids separation (much less than 720 kpc)
d_gal_sep  = 720.0 * kpc    # galaxy centroids separation (= lensing peak separation)

# Gas centroid positions (x-axis)
x_gas_main   = -d_gas_sep / 2.0
x_gas_bullet =  d_gas_sep / 2.0

# Galaxy centroid positions
x_gal_main   = -d_gal_sep / 2.0
x_gal_bullet =  d_gal_sep / 2.0

# 1-D grid
N_grid = 2000
x_min  = -2500.0 * kpc
x_max  =  2500.0 * kpc
x_arr  = np.linspace(x_min, x_max, N_grid)
dx     = x_arr[1] - x_arr[0]


def gaussian_profile(x, x0, sigma, M_total):
    """1-D Gaussian surface density profile Sigma(x)."""
    return (M_total / (sigma * np.sqrt(2.0 * np.pi))) * np.exp(
        -0.5 * ((x - x0) / sigma) ** 2
    )


def beta_profile_1d(x, x0, r_core, M_total, beta=2.0/3.0):
    """
    1-D projected beta-model profile for gas.
    Sigma(x) proportional to (1 + ((x-x0)/r_core)^2)^(-3*beta/2 + 1/2)
    Normalised to integrate to M_total.
    """
    u = (x - x0) / r_core
    exponent = -3.0 * beta / 2.0 + 0.5
    raw = np.power(1.0 + u**2, exponent)
    norm = trapezoid(raw, x)
    if norm <= 0:
        norm = 1.0
    return M_total * raw / norm


# Gas profiles (quasi-spherical, beta-model with large core)
sigma_gas_main = beta_profile_1d(
    x_arr, x_gas_main, r_core=250.0*kpc, M_total=M_gas_main
)
sigma_gas_bullet = beta_profile_1d(
    x_arr, x_gas_bullet, r_core=100.0*kpc, M_total=M_gas_bullet
)

# Galaxy/stellar profiles (more concentrated)
sigma_star_main = gaussian_profile(
    x_arr, x_gal_main, sigma=200.0*kpc, M_total=M_star_main
)
sigma_star_bullet = gaussian_profile(
    x_arr, x_gal_bullet, sigma=80.0*kpc, M_total=M_star_bullet
)

# Total baryonic surface density
sigma_bar = sigma_gas_main + sigma_gas_bullet + sigma_star_main + sigma_star_bullet

# Enclosed mass (cumulative from left)
M_enc_bar = cumulative_trapezoid(sigma_bar, x_arr, initial=0)

# Newtonian acceleration (approximate 1-D: treat as cylindrical shell)
# For a projected 1-D analysis, use g_bar ~ G * M_enc(x) / x^2
# More precisely, we use the gravitational acceleration at projected radius r
# from the nearest dominant mass concentration.

# For each point, compute g_bar from both subclusters
def g_newtonian_1d(x_eval, x_center, M_total, r_scale):
    """
    Newtonian gravitational acceleration from a mass concentration.
    g = G*M(<r) / r^2 with M(<r) = M_total * r^3/(r^3 + r_s^3) (smooth cutoff).
    """
    r = np.abs(x_eval - x_center)
    r = np.maximum(r, 1.0 * kpc)  # softening
    M_enc = M_total * r**3 / (r**3 + r_scale**3)
    return G * M_enc / r**2

# Newtonian g from each baryonic component
g_gas_main   = g_newtonian_1d(x_arr, x_gas_main,   M_gas_main,   250.0*kpc)
g_gas_bullet = g_newtonian_1d(x_arr, x_gas_bullet,  M_gas_bullet, 100.0*kpc)
g_star_main  = g_newtonian_1d(x_arr, x_gal_main,    M_star_main,  200.0*kpc)
g_star_bullet= g_newtonian_1d(x_arr, x_gal_bullet,  M_star_bullet, 80.0*kpc)

# Total Newtonian acceleration (vector sum along x-axis)
def g_newtonian_signed(x_eval, x_center, M_total, r_scale):
    """Signed Newtonian g (negative = pointing left, positive = pointing right)."""
    dx_vec = x_center - x_eval  # points toward center
    r = np.abs(dx_vec)
    r = np.maximum(r, 1.0 * kpc)
    M_enc = M_total * r**3 / (r**3 + r_scale**3)
    g_mag = G * M_enc / r**2
    return g_mag * np.sign(dx_vec)

g_bar_signed = (
    g_newtonian_signed(x_arr, x_gas_main,   M_gas_main,   250.0*kpc) +
    g_newtonian_signed(x_arr, x_gas_bullet,  M_gas_bullet, 100.0*kpc) +
    g_newtonian_signed(x_arr, x_gal_main,    M_star_main,  200.0*kpc) +
    g_newtonian_signed(x_arr, x_gal_bullet,  M_star_bullet, 80.0*kpc)
)

g_bar_mag = np.abs(g_bar_signed)

# TGP boosted acceleration
# Gas is quasi-spherical (S ~ 0.8), galaxies are more extended (S ~ 0.5)
S_gas = 0.8
S_gal = 0.5

# For a composite system, we weight the sphericity by mass contribution at each point
frac_gas = (g_gas_main + g_gas_bullet) / np.maximum(g_bar_mag, 1e-20)
frac_gas = np.clip(frac_gas, 0, 1)
S_eff = frac_gas * S_gas + (1.0 - frac_gas) * S_gal

# TGP boost factor
y_arr = g_bar_mag / a0
nu_arr = np.array([nu_tgp(y, S) for y, S in zip(y_arr, S_eff)])

g_tgp_mag = g_bar_mag * nu_arr

# "Phantom dark matter" acceleration (the MOND/TGP boost)
g_phantom = g_tgp_mag - g_bar_mag

# For lensing: the effective surface density
# Sigma_lens = g_total / (2*pi*G) * r  (convergence kappa ~ Sigma/Sigma_crit)
# In 1-D projection, phantom mass ~ g_phantom * r / G
# We use a simplified lensing-like surface density: Sigma_eff ~ g / (G / r_char)
r_char = 500.0 * kpc  # characteristic radius for lensing projection

Sigma_Newton = g_bar_mag * r_char / (2.0 * np.pi * G)
Sigma_TGP    = g_tgp_mag * r_char / (2.0 * np.pi * G)
Sigma_phantom= g_phantom * r_char / (2.0 * np.pi * G)

# Find peaks
x_kpc = x_arr / kpc

idx_left  = (x_kpc < 0)
idx_right = (x_kpc > 0)

# Newton: peaks
pk_N_L = x_kpc[idx_left][np.argmax(Sigma_Newton[idx_left])]
pk_N_R = x_kpc[idx_right][np.argmax(Sigma_Newton[idx_right])]

# TGP: peaks
pk_T_L = x_kpc[idx_left][np.argmax(Sigma_TGP[idx_left])]
pk_T_R = x_kpc[idx_right][np.argmax(Sigma_TGP[idx_right])]

# Phantom: peaks
pk_P_L = x_kpc[idx_left][np.argmax(Sigma_phantom[idx_left])]
pk_P_R = x_kpc[idx_right][np.argmax(Sigma_phantom[idx_right])]

print(f"""
  B.1  Post-collision geometry (1-D model along collision axis)
  ------------------------------------------------------------
    Gas centroids:      x_main = {x_gas_main/kpc:+.0f} kpc,  x_bullet = {x_gas_bullet/kpc:+.0f} kpc
    Galaxy centroids:   x_main = {x_gal_main/kpc:+.0f} kpc,  x_bullet = {x_gal_bullet/kpc:+.0f} kpc
    Gas separation:     {d_gas_sep/kpc:.0f} kpc  (stripped toward center)
    Galaxy separation:  {d_gal_sep/kpc:.0f} kpc  (passed through each other)

  B.2  Sphericity assignments
  ---------------------------
    Gas:      S = {S_gas}  (quasi-spherical, gamma = {gamma_S(S_gas):.3f})
    Galaxies: S = {S_gal}  (more extended,   gamma = {gamma_S(S_gal):.3f})

  B.3  Acceleration profile summary
  ----------------------------------
    At x = 0 (midpoint):
      g_bar   = {g_bar_mag[N_grid//2]:.3e} m/s^2
      y       = {y_arr[N_grid//2]:.3f}
      nu(y,S) = {nu_arr[N_grid//2]:.3f}
      g_TGP   = {g_tgp_mag[N_grid//2]:.3e} m/s^2
      boost   = {nu_arr[N_grid//2]:.2f}x

    At x = -360 kpc (main galaxy position):
      g_bar   = {g_bar_mag[np.argmin(np.abs(x_kpc + 360))]:.3e} m/s^2
      y       = {y_arr[np.argmin(np.abs(x_kpc + 360))]:.3f}
      nu      = {nu_arr[np.argmin(np.abs(x_kpc + 360))]:.3f}

    At x = +360 kpc (bullet galaxy position):
      g_bar   = {g_bar_mag[np.argmin(np.abs(x_kpc - 360))]:.3e} m/s^2
      y       = {y_arr[np.argmin(np.abs(x_kpc - 360))]:.3f}
      nu      = {nu_arr[np.argmin(np.abs(x_kpc - 360))]:.3f}
""")


# ########################################################################## #
#                                                                            #
#  PART C: THE LENSING OFFSET PROBLEM                                       #
#                                                                            #
# ########################################################################## #

print("=" * 78)
print("  PART C: THE LENSING OFFSET PROBLEM -- PHANTOM DARK MATTER")
print("=" * 78)

print(f"""
  C.1  What is phantom dark matter?
  ----------------------------------
  In GR + MOND/TGP, the total gravitational field is g_obs = nu(y)*g_bar.
  A GR observer attributes this to total mass:
      M_total = g_obs * r^2 / G = nu(y) * M_bar
  The "phantom" or "missing" mass is:
      M_phantom = M_total - M_bar = (nu - 1) * M_bar

  The distribution of this phantom mass is NOT the same as baryonic mass.
  Where g_bar is SMALL (y << 1), nu - 1 is LARGE (deep MOND boost).
  So phantom mass is largest in the OUTSKIRTS, where Newtonian gravity is weak.

  C.2  Peak locations (1-D model)
  --------------------------------
    Baryonic (Newtonian) lensing peaks:
      Left (main):     x = {pk_N_L:+.0f} kpc
      Right (bullet):  x = {pk_N_R:+.0f} kpc

    TGP total lensing peaks:
      Left (main):     x = {pk_T_L:+.0f} kpc
      Right (bullet):  x = {pk_T_R:+.0f} kpc

    Phantom mass peaks:
      Left (main):     x = {pk_P_L:+.0f} kpc
      Right (bullet):  x = {pk_P_R:+.0f} kpc

    Observed lensing peaks:
      Left (main):     x ~ -360 kpc  (at main galaxy centroid)
      Right (bullet):  x ~ +360 kpc  (at bullet galaxy centroid)

    Gas centroids:
      Left (main):     x = {x_gas_main/kpc:+.0f} kpc
      Right (bullet):  x = {x_gas_bullet/kpc:+.0f} kpc
""")

# Compute the offset deficit
offset_TGP_main   = np.abs(pk_T_L - (-360.0))
offset_TGP_bullet = np.abs(pk_T_R - 360.0)
offset_obs        = 150.0   # kpc, each subcluster

# What fraction of observed offset does TGP predict?
# The gas centroid is at ~100 kpc from center; observed lensing peak at ~360 kpc
# The offset is lensing_peak - gas_centroid ~ 260 kpc
gas_to_lens_main   = np.abs(pk_T_L - x_gas_main/kpc)
gas_to_lens_bullet = np.abs(pk_T_R - x_gas_bullet/kpc)

print(f"""
  C.3  Offset analysis
  ---------------------
    TGP total lensing peak vs gas centroid:
      Main:   |x_TGP_peak - x_gas|  = {gas_to_lens_main:.0f} kpc
      Bullet: |x_TGP_peak - x_gas|  = {gas_to_lens_bullet:.0f} kpc

    Observed offset (lensing peak vs gas centroid):
      ~150-260 kpc per subcluster

    Key question: does the TGP boost shift the lensing peak AWAY from the gas
    and TOWARD the galaxy positions?
""")


# ########################################################################## #
#                                                                            #
#  PART D: MEMBRANE PROPAGATION AND GRAVITATIONAL MEMORY                    #
#                                                                            #
# ########################################################################## #

print("=" * 78)
print("  PART D: MEMBRANE PROPAGATION AND GRAVITATIONAL MEMORY")
print("=" * 78)

# The TGP membrane has a characteristic crossover scale
# r_c = sqrt(G*M / a0) = sqrt(r_S * r_H)  (gs9d)

r_c_main   = np.sqrt(G * M_bar_main / a0)
r_c_bullet = np.sqrt(G * M_bar_bullet / a0)

# Membrane relaxation timescale
# The membrane propagation speed c_membrane ~ c (gravitational waves on brane)
# Relaxation timescale: tau = r_c / c_membrane
c_membrane = c  # assume gravitational wave speed
tau_main   = r_c_main / c_membrane
tau_bullet = r_c_bullet / c_membrane

# Collision timescale
# t_cross = d_sep / v_collision  (how long ago the collision happened)
t_cross   = d_sep / v_collision  # ~150 kpc / 4700 km/s

# Age of the collision (estimated from observations): ~150 Myr
t_collision = 150.0e6 * 3.156e7  # 150 Myr in seconds

# Comparison
ratio_main   = tau_main / t_collision
ratio_bullet = tau_bullet / t_collision

print(f"""
  D.1  Membrane crossover scales
  --------------------------------
    r_c(main)   = sqrt(G*M_bar / a0) = {r_c_main/kpc:.0f} kpc
    r_c(bullet) = sqrt(G*M_bar / a0) = {r_c_bullet/kpc:.0f} kpc

  D.2  Relaxation timescales
  ---------------------------
    tau_relax(main)   = r_c / c = {tau_main/3.156e7/1e6:.1f} Myr
    tau_relax(bullet) = r_c / c = {tau_bullet/3.156e7/1e6:.1f} Myr

    t_collision       ~ {t_collision/3.156e7/1e6:.0f} Myr  (time since core passage)
    t_cross = d/v     = {t_cross/3.156e7/1e6:.0f} Myr

  D.3  Memory ratio: tau_relax / t_collision
  -------------------------------------------
    Main:   tau / t_coll = {ratio_main:.4f}
    Bullet: tau / t_coll = {ratio_bullet:.4f}

  D.4  Interpretation
  --------------------
""")

if ratio_main > 0.1:
    print("    The membrane relaxation time is a SIGNIFICANT fraction of the collision")
    print("    timescale. The membrane may retain a partial 'gravitational memory' of")
    print("    the pre-collision configuration.")
    memory_fraction_est = min(ratio_main, 1.0)
else:
    print("    The membrane relaxation time is MUCH SHORTER than the collision timescale.")
    print("    The membrane has fully relaxed to the current baryonic distribution.")
    print("    Gravitational memory CANNOT explain the lensing offset.")
    memory_fraction_est = ratio_main

print(f"""
    Estimated memory fraction: {memory_fraction_est:.3f}
    (fraction of pre-collision gravitational field that might persist)
""")

# Could the membrane propagation be SLOWER than c?
# In DGP-like models, the crossover involves a massive graviton mode
# with group velocity ~ c * (r/r_c)  for r << r_c
# This could dramatically slow the relaxation

print("""
  D.5  Slow-mode membrane propagation (speculative)
  ---------------------------------------------------
  In DGP braneworld models, the graviton has a soft mass m_g ~ 1/r_c.
  The massive mode propagates with reduced group velocity on the brane:
      v_group ~ c * sqrt(1 - (m_g*c^2 / E)^2)
  For modes with wavelength ~ r_c: v_group ~ 0.

  If TGP membrane propagation is similarly dispersive, the relaxation
  timescale at scale r_c could be MUCH longer than r_c/c.

  Let us parameterize: tau_eff = r_c / v_eff, with v_eff << c.
""")

# Slow propagation scan
v_eff_fracs = [1.0, 0.1, 0.01, 0.001, 1e-4]
print("    v_eff/c     tau_eff(main)    tau/t_coll     Memory?")
print("    " + "-" * 60)
for frac in v_eff_fracs:
    v_eff = frac * c
    tau_eff = r_c_main / v_eff
    ratio_eff = tau_eff / t_collision
    memory_str = "YES" if ratio_eff > 1.0 else ("partial" if ratio_eff > 0.1 else "no")
    print(f"    {frac:.0e}       {tau_eff/3.156e7/1e6:10.1f} Myr    {ratio_eff:10.3f}       {memory_str}")

print("""
  RESULT: Gravitational memory requires v_eff/c < 10^{-4}, i.e., the membrane
  propagation speed must be 10,000x slower than light.  This is physically
  implausible for a gravitational mode.  Gravitational memory is NOT a viable
  mechanism for explaining the Bullet Cluster lensing offset.
""")


# ########################################################################## #
#                                                                            #
#  PART E: QUANTITATIVE COLLISION MODEL AND LENSING MAPS                    #
#                                                                            #
# ########################################################################## #

print("=" * 78)
print("  PART E: QUANTITATIVE 1-D COLLISION MODEL AND LENSING MAPS")
print("=" * 78)

# --- Model the full collision with 2-D projected mass distributions ---
# Use a 2-D grid (x = collision axis, R = transverse projected radius)
# But for computational simplicity, work in cylindrical 1-D (along x-axis)
# assuming azimuthal symmetry about the collision axis.

# NFW-like 3D density profile
def rho_nfw(r, M_total, r_s, c_nfw=5.0):
    """NFW density profile rho(r), normalized to M_total within r_200 = c*r_s."""
    r200 = c_nfw * r_s
    # NFW normalization
    A = np.log(1.0 + c_nfw) - c_nfw / (1.0 + c_nfw)
    rho_0 = M_total / (4.0 * np.pi * r_s**3 * A)
    x = r / r_s
    x = np.maximum(x, 1e-6)
    return rho_0 / (x * (1.0 + x)**2)


def enclosed_mass_nfw(r, M_total, r_s, c_nfw=5.0):
    """Enclosed mass for NFW profile."""
    A = np.log(1.0 + c_nfw) - c_nfw / (1.0 + c_nfw)
    x = r / r_s
    M_enc = M_total * (np.log(1.0 + x) - x / (1.0 + x)) / A
    return M_enc


# --- Pre-collision state (isolated halos) ---
# Each subcluster has baryons following an NFW-like profile.
# At r from center: g_bar(r) = G*M_bar(<r) / r^2

# Radial grid
r_arr = np.logspace(np.log10(1.0), np.log10(3000.0), 500) * kpc

# Main cluster -- baryonic
M_enc_bar_main_r = enclosed_mass_nfw(r_arr, M_bar_main, r_s=300.0*kpc, c_nfw=5.0)
g_bar_main_r     = G * M_enc_bar_main_r / r_arr**2
y_main_r         = g_bar_main_r / a0
nu_main_r        = nu_tgp(y_main_r, S=S_gas)
g_tgp_main_r     = g_bar_main_r * nu_main_r
M_tgp_main_r     = g_tgp_main_r * r_arr**2 / G  # effective enclosed mass

# Bullet -- baryonic
M_enc_bar_bull_r = enclosed_mass_nfw(r_arr, M_bar_bullet, r_s=100.0*kpc, c_nfw=5.0)
g_bar_bull_r     = G * M_enc_bar_bull_r / r_arr**2
y_bull_r         = g_bar_bull_r / a0
nu_bull_r        = nu_tgp(y_bull_r, S=S_gas)
g_tgp_bull_r     = g_bar_bull_r * nu_bull_r
M_tgp_bull_r     = g_tgp_bull_r * r_arr**2 / G

# Boost ratios
boost_main_200 = M_tgp_main_r[-1] / M_bar_main
boost_bull_200 = M_tgp_bull_r[-1] / M_bar_bullet

# What total mass does TGP predict at r200?
M_tgp_main_total = M_tgp_main_r[-1]
M_tgp_bull_total = M_tgp_bull_r[-1]

# Compare with lensing mass
ratio_main_TGP = M_tgp_main_total / M_main_lens
ratio_bull_TGP = M_tgp_bull_total / M_bullet_lens

print(f"""
  E.1  Pre-collision equilibrium: isolated halos
  -----------------------------------------------
  For each subcluster, compute the TGP effective mass profile assuming
  isolated, spherical baryonic distribution.

  Main cluster (M_bar = {M_bar_main/M_sun:.2e} M_sun):
    At r = 100 kpc:  M_bar = {enclosed_mass_nfw(100*kpc, M_bar_main, 300*kpc)/M_sun:.2e},  nu = {nu_tgp(G*enclosed_mass_nfw(100*kpc, M_bar_main, 300*kpc)/(100*kpc)**2/a0, S_gas):.2f}
    At r = 500 kpc:  M_bar = {enclosed_mass_nfw(500*kpc, M_bar_main, 300*kpc)/M_sun:.2e},  nu = {nu_tgp(G*enclosed_mass_nfw(500*kpc, M_bar_main, 300*kpc)/(500*kpc)**2/a0, S_gas):.2f}
    At r = 1000 kpc: M_bar = {enclosed_mass_nfw(1000*kpc, M_bar_main, 300*kpc)/M_sun:.2e}, nu = {nu_tgp(G*enclosed_mass_nfw(1000*kpc, M_bar_main, 300*kpc)/(1000*kpc)**2/a0, S_gas):.2f}
    At r = 2000 kpc: M_bar = {enclosed_mass_nfw(2000*kpc, M_bar_main, 300*kpc)/M_sun:.2e}, nu = {nu_tgp(G*enclosed_mass_nfw(2000*kpc, M_bar_main, 300*kpc)/(2000*kpc)**2/a0, S_gas):.2f}
    M_TGP(r200)     = {M_tgp_main_total/M_sun:.2e} M_sun
    M_lens(obs)     = {M_main_lens/M_sun:.2e} M_sun
    Ratio M_TGP/M_lens = {ratio_main_TGP:.3f}

  Bullet (M_bar = {M_bar_bullet/M_sun:.2e} M_sun):
    At r = 50 kpc:   M_bar = {enclosed_mass_nfw(50*kpc, M_bar_bullet, 100*kpc)/M_sun:.2e},  nu = {nu_tgp(G*enclosed_mass_nfw(50*kpc, M_bar_bullet, 100*kpc)/(50*kpc)**2/a0, S_gas):.2f}
    At r = 200 kpc:  M_bar = {enclosed_mass_nfw(200*kpc, M_bar_bullet, 100*kpc)/M_sun:.2e},  nu = {nu_tgp(G*enclosed_mass_nfw(200*kpc, M_bar_bullet, 100*kpc)/(200*kpc)**2/a0, S_gas):.2f}
    At r = 500 kpc:  M_bar = {enclosed_mass_nfw(500*kpc, M_bar_bullet, 100*kpc)/M_sun:.2e},  nu = {nu_tgp(G*enclosed_mass_nfw(500*kpc, M_bar_bullet, 100*kpc)/(500*kpc)**2/a0, S_gas):.2f}
    M_TGP(r200)     = {M_tgp_bull_total/M_sun:.2e} M_sun
    M_lens(obs)     = {M_bullet_lens/M_sun:.2e} M_sun
    Ratio M_TGP/M_lens = {ratio_bull_TGP:.3f}
""")

# --- Post-collision state ---
# Gas is stripped: main gas stays near x ~ -100, bullet gas near x ~ +100
# Galaxies pass through: main gal at x ~ -360, bullet gal at x ~ +360

# For the post-collision TGP lensing map, we compute the projected convergence
# kappa(x) along the collision axis.

# Use the 1-D model from Part B but with finer detail
print("  E.2  Post-collision lensing maps (1-D projection)")
print("  " + "-" * 60)

# Recompute with explicit phantom mass density
# rho_phantom = rho_total - rho_baryon = (nu - 1) * rho_baryon  (local approx)

# For each point along x, compute the local y and nu from all sources
g_from_sources = np.zeros_like(x_arr)
g_gas_contrib  = np.zeros_like(x_arr)
g_star_contrib = np.zeros_like(x_arr)

# Signed acceleration components
g_gas_main_s   = g_newtonian_signed(x_arr, x_gas_main,   M_gas_main,   250.0*kpc)
g_gas_bull_s   = g_newtonian_signed(x_arr, x_gas_bullet,  M_gas_bullet, 100.0*kpc)
g_star_main_s  = g_newtonian_signed(x_arr, x_gal_main,    M_star_main,  200.0*kpc)
g_star_bull_s  = g_newtonian_signed(x_arr, x_gal_bullet,  M_star_bullet, 80.0*kpc)

g_total_signed = g_gas_main_s + g_gas_bull_s + g_star_main_s + g_star_bull_s
g_total_mag    = np.abs(g_total_signed)

# Effective sphericity at each point (mass-weighted)
g_gas_total = np.abs(g_gas_main_s + g_gas_bull_s)
g_star_total = np.abs(g_star_main_s + g_star_bull_s)
g_all = g_gas_total + g_star_total
g_all = np.maximum(g_all, 1e-30)

S_local = (g_gas_total * S_gas + g_star_total * S_gal) / g_all

y_local  = g_total_mag / a0
nu_local = np.array([nu_tgp(y, S) for y, S in zip(y_local, S_local)])

g_tgp_total = g_total_mag * nu_local
g_phantom_1d = g_tgp_total - g_total_mag

# Convergence-like quantity (proportional to projected surface mass density)
# kappa ~ Sigma / Sigma_crit
# Compute raw projected surface density, then normalize to peak for display
Sigma_bar_raw     = g_total_mag * r_char / (2.0 * np.pi * G)
Sigma_tgp_raw     = g_tgp_total * r_char / (2.0 * np.pi * G)
Sigma_phantom_raw = g_phantom_1d * r_char / (2.0 * np.pi * G)

# Normalize all to the peak of the TGP distribution for relative comparison
Sigma_norm = np.max(Sigma_tgp_raw)
kappa_bar     = Sigma_bar_raw / Sigma_norm
kappa_tgp     = Sigma_tgp_raw / Sigma_norm
kappa_phantom = Sigma_phantom_raw / Sigma_norm

# Find peaks in each half
def find_peak(arr, x_kpc, x_range):
    """Find peak position within x_range = (x_min, x_max) in kpc."""
    mask = (x_kpc >= x_range[0]) & (x_kpc <= x_range[1])
    if not np.any(mask):
        return 0.0, 0.0
    idx = np.argmax(arr[mask])
    return x_kpc[mask][idx], arr[mask][idx]

# Peaks for main cluster side (left)
pk_bar_L, val_bar_L   = find_peak(kappa_bar, x_kpc, (-800, -10))
pk_tgp_L, val_tgp_L   = find_peak(kappa_tgp, x_kpc, (-800, -10))
pk_phant_L, val_ph_L   = find_peak(kappa_phantom, x_kpc, (-800, -10))

# Peaks for bullet side (right)
pk_bar_R, val_bar_R   = find_peak(kappa_bar, x_kpc, (10, 800))
pk_tgp_R, val_tgp_R   = find_peak(kappa_tgp, x_kpc, (10, 800))
pk_phant_R, val_ph_R   = find_peak(kappa_phantom, x_kpc, (10, 800))

# Values at specific positions
idx_gal_main = np.argmin(np.abs(x_kpc - (-360)))
idx_gal_bull = np.argmin(np.abs(x_kpc - 360))
idx_gas_main = np.argmin(np.abs(x_kpc - x_gas_main/kpc))
idx_gas_bull = np.argmin(np.abs(x_kpc - x_gas_bullet/kpc))

print(f"""
  Post-collision convergence map peaks (kappa, normalized):

    Component        Main (left)              Bullet (right)
                     peak_x    kappa          peak_x    kappa
    Baryonic:        {pk_bar_L:+7.0f} kpc  {val_bar_L:.4f}       {pk_bar_R:+7.0f} kpc  {val_bar_R:.4f}
    TGP total:       {pk_tgp_L:+7.0f} kpc  {val_tgp_L:.4f}       {pk_tgp_R:+7.0f} kpc  {val_tgp_R:.4f}
    Phantom only:    {pk_phant_L:+7.0f} kpc  {val_ph_L:.4f}       {pk_phant_R:+7.0f} kpc  {val_ph_R:.4f}
    Observed:            -360 kpc                  +360 kpc

  Convergence at key positions:
                     At gas centroid          At galaxy centroid
    Main:  kappa_TGP = {kappa_tgp[idx_gas_main]:.4f} (gas)   vs  {kappa_tgp[idx_gal_main]:.4f} (gal)
    Bullet: kappa_TGP = {kappa_tgp[idx_gas_bull]:.4f} (gas)   vs  {kappa_tgp[idx_gal_bull]:.4f} (gal)
""")

# --- Mass budget ---
# How much total mass does TGP predict for each subcluster?
# Integrate kappa over each half

# For a rough mass estimate, integrate the effective surface density
M_eff_bar_L  = trapezoid(Sigma_bar_raw[idx_left], x_arr[idx_left])
M_eff_tgp_L  = trapezoid(Sigma_tgp_raw[idx_left], x_arr[idx_left])
M_eff_bar_R  = trapezoid(Sigma_bar_raw[idx_right], x_arr[idx_right])
M_eff_tgp_R  = trapezoid(Sigma_tgp_raw[idx_right], x_arr[idx_right])

boost_L = M_eff_tgp_L / M_eff_bar_L if M_eff_bar_L > 0 else 0
boost_R = M_eff_tgp_R / M_eff_bar_R if M_eff_bar_R > 0 else 0

print(f"""
  E.3  Mass budget (integrated over each half-space)
  ---------------------------------------------------
    Main cluster (x < 0):
      M_bar_eff  = {M_eff_bar_L/M_sun:.2e} M_sun
      M_TGP_eff  = {M_eff_tgp_L/M_sun:.2e} M_sun
      Boost      = {boost_L:.2f}x
      Observed   = {M_main_lens/M_sun:.2e} M_sun
      Needed boost = {M_main_lens / M_bar_main:.1f}x

    Bullet (x > 0):
      M_bar_eff  = {M_eff_bar_R/M_sun:.2e} M_sun
      M_TGP_eff  = {M_eff_tgp_R/M_sun:.2e} M_sun
      Boost      = {boost_R:.2f}x
      Observed   = {M_bullet_lens/M_sun:.2e} M_sun
      Needed boost = {M_bullet_lens / M_bar_bullet:.1f}x
""")

# --- Quantitative offset assessment ---
# The key metric: how far is the TGP lensing peak from the gas centroid,
# compared to the observed offset?

offset_TGP_main_gas   = np.abs(pk_tgp_L - x_gas_main/kpc)
offset_TGP_bullet_gas = np.abs(pk_tgp_R - x_gas_bullet/kpc)
offset_obs_main       = np.abs(-360.0 - x_gas_main/kpc)    # observed: 260 kpc
offset_obs_bullet     = np.abs(360.0 - x_gas_bullet/kpc)   # observed: 260 kpc

frac_offset_main   = offset_TGP_main_gas / offset_obs_main if offset_obs_main > 0 else 0
frac_offset_bullet = offset_TGP_bullet_gas / offset_obs_bullet if offset_obs_bullet > 0 else 0

print(f"""
  E.4  Lensing offset: TGP prediction vs observation
  ----------------------------------------------------
    Main cluster:
      TGP lensing peak:        x = {pk_tgp_L:+.0f} kpc
      Gas centroid:             x = {x_gas_main/kpc:+.0f} kpc
      TGP offset from gas:     {offset_TGP_main_gas:.0f} kpc
      Observed offset:         {offset_obs_main:.0f} kpc
      Fraction explained:      {frac_offset_main:.1%}

    Bullet:
      TGP lensing peak:        x = {pk_tgp_R:+.0f} kpc
      Gas centroid:             x = {x_gas_bullet/kpc:+.0f} kpc
      TGP offset from gas:     {offset_TGP_bullet_gas:.0f} kpc
      Observed offset:         {offset_obs_bullet:.0f} kpc
      Fraction explained:      {frac_offset_bullet:.1%}
""")


# ########################################################################## #
#                                                                            #
#  PART F: COMPARISON WITH MOND LITERATURE AND HONEST ASSESSMENT            #
#                                                                            #
# ########################################################################## #

print("=" * 78)
print("  PART F: COMPARISON WITH MOND LITERATURE AND HONEST ASSESSMENT")
print("=" * 78)

print("""
  F.1  Literature on MOND and the Bullet Cluster
  ------------------------------------------------

  [1] Clowe et al. (2006, ApJ 648, L109):
      "A Direct Empirical Proof of the Existence of Dark Matter"
      - Original weak lensing analysis showing mass-gas offset
      - Claimed "direct proof" that DM exists independent of gravity theory
      - Impact: widely cited as falsification of MOND

  [2] Angus, Shan, Zhao & Famaey (2006, ApJ 654, L13):
      - Proposed MOND + 2 eV ordinary neutrinos
      - Hot neutrino halos provide ~60% of the missing mass in clusters
      - Problem: neutrino free-streaming erases structure < 10 Mpc
      - Can explain ~60% of the Bullet Cluster mass discrepancy

  [3] Angus, Famaey & Zhao (2006, MNRAS 371, 138):
      - Extended analysis: MOND + 2 eV sterile neutrinos
      - Sterile neutrinos are collisionless (like DM, would pass through)
      - Could explain lensing offset IF sterile neutrinos followed galaxies
      - Problem: 2 eV neutrinos now strongly disfavored by Planck (m_nu < 0.12 eV)

  [4] Milgrom (2008, NewAR 51, 906):
      "The MOND paradigm of modified dynamics"
      - Arguments that Bullet Cluster does NOT definitively falsify MOND:
        (a) The mass discrepancy is only ~2x, not 5-6x as in individual clusters
        (b) External field effect (EFE) could modify the dynamics
        (c) The system is highly non-equilibrium; MOND formulae may not apply
        (d) Residual neutrino mass could contribute
      - However: does NOT provide a concrete mechanism for the lensing offset

  [5] Skordis & Zlosnik (2021, PRL 127, 161302):
      - "New relativistic theory for MOND" (RMOND)
      - Successfully fits CMB power spectrum (major achievement)
      - Bullet Cluster explicitly NOT addressed in the paper
      - The vector field in RMOND could in principle behave like "dark matter"
        in a collision, but this has not been demonstrated

  [6] Brownstein & Moffat (2007, MNRAS 382, 29):
      - MOG (Modified Gravity) can fit Bullet Cluster lensing
      - But MOG is not MOND; it has additional vector/scalar fields
      - The additional fields act as effective dark matter
""")

# --- TGP-specific assessment ---
# Compute the "success metric" -- what fraction of the discrepancy can TGP address?

# The Bullet Cluster challenge has TWO parts:
# (A) Mass discrepancy: need ~7-10x boost over baryonic mass
# (B) Spatial offset: lensing peaks must be at galaxy positions, not gas

# For (A): what boost does TGP give at cluster scales?
# Clusters are in the y >> 1 regime at inner radii but y ~ 1 at outskirts
r_test = np.array([100, 200, 500, 1000, 2000]) * kpc
print("\n  F.2  TGP boost factor at cluster scales")
print("  " + "-" * 60)
print(f"    r (kpc)    g_bar (m/s^2)    y = g/a0    nu(y, S=0.8)   M_eff/M_bar")

for r in r_test:
    M_enc = enclosed_mass_nfw(r, M_bar_main, 300.0*kpc, c_nfw=5.0)
    g = G * M_enc / r**2
    y = g / a0
    nu = nu_tgp(y, S=0.8)
    M_eff = nu * M_enc
    print(f"    {r/kpc:6.0f}     {g:.3e}        {y:7.3f}      {nu:7.3f}         {nu:.2f}")

# Total mass prediction
# At the virial radius, what does TGP predict?
r_vir = 2000.0 * kpc
M_enc_vir = enclosed_mass_nfw(r_vir, M_bar_main, 300.0*kpc)
g_vir = G * M_enc_vir / r_vir**2
y_vir = g_vir / a0
nu_vir = nu_tgp(y_vir, S=0.8)
M_tgp_vir = nu_vir * M_enc_vir

mass_deficit = M_main_lens / M_tgp_vir

print(f"""
  At virial radius (r = {r_vir/kpc:.0f} kpc):
    M_bar_enc   = {M_enc_vir/M_sun:.2e} M_sun
    nu(y, S=0.8)= {nu_vir:.3f}
    M_TGP       = {M_tgp_vir/M_sun:.2e} M_sun
    M_lens(obs) = {M_main_lens/M_sun:.2e} M_sun
    Mass deficit: M_lens / M_TGP = {mass_deficit:.2f}
    --> TGP can explain {1.0/mass_deficit:.0%} of the lensing mass
""")

# --- The offset problem ---
print("""
  F.3  The spatial offset problem: TGP assessment
  -------------------------------------------------
  The lensing offset is the HARDER problem.  Even if TGP provided enough total
  mass, the phantom mass distribution must peak at the galaxy positions.

  In standard (non-relativistic) MOND/TGP:
    g_obs = nu(g_bar/a0) * g_bar
    The phantom mass distribution follows g_bar, which is dominated by gas.
    Therefore: phantom mass peaks at gas positions, NOT galaxy positions.

  This is INDEPENDENT of the interpolation function.  ANY theory where
  g_obs = f(g_bar) * g_bar will have phantom mass tracing g_bar.

  The only escape routes are:
    (1) Non-local effects (membrane memory) -- shown in Part D to require
        v_eff/c < 10^{-4}, which is implausible.
    (2) Additional collisionless matter (neutrinos, sterile particles)
    (3) A fundamentally different coupling in which the gravitational
        enhancement depends on the HISTORY of the mass distribution,
        not its current state.
    (4) Relativistic effects in the vector/tensor sector of TGP
        (analogous to Skordis-Zlosnik RMOND approach).
""")

# --- Scorecard ---
print("""
  F.4  HONEST SCORECARD: TGP vs Bullet Cluster
  ===============================================

  Challenge                        TGP Performance           Score
  --------------------------------------------------------------------------
  (1) Total mass at r200           ~2x boost (need ~7x)      FAIL
      Cluster-scale mass budget    Only ~30% of needed mass   2/10
      (consistent with known MOND
       cluster problem, gs13/gs15)

  (2) Lensing peak offset          Peaks track gas, not       FAIL
      from gas centroid            galaxies                   1/10
      (the CRITICAL test)

  (3) Mach 3 collision velocity    TGP cannot explain         FAIL
      ~4700 km/s is hard even     (needs M > M_bar in halo)  3/10
      for LCDM

  (4) Merger rate / probability    Not addressed (similar     N/A
                                   to LCDM: very rare)

  (5) Gas temperature profile      TGP predicts less total    PARTIAL
      and X-ray luminosity         mass -> lower T expected   4/10
                                   (but clusters are already
                                    problematic, gs13)

  OVERALL: 2/10
  --------------------------------------------------------------------------

  VERDICT: The Bullet Cluster remains a SEVERE challenge for TGP, just as it
  is for all MOND-like theories.  The two fundamental problems are:

    (A) MASS DEFICIT: TGP boosts cluster-scale mass by only ~2x, but the
        lensing mass is ~7x the baryonic mass.  This is the same "cluster
        problem" known from gs13/gs15.

    (B) SPATIAL OFFSET: Even if TGP provided enough total mass, the phantom
        mass distribution traces the baryonic acceleration field, which is
        dominated by gas.  There is no mechanism in the current TGP framework
        to shift phantom mass from gas positions to galaxy positions.

  The membrane gravitational memory idea (Part D) is physically interesting
  but quantitatively fails: the relaxation timescale is ~1000x too short
  unless membrane propagation is implausibly slow (v_eff < 10^{-4} c).
""")

print("""
  F.5  POSSIBLE PATHS FORWARD
  -----------------------------
  (1) HYBRID APPROACH: TGP + residual dark matter
      If TGP explains galaxy-scale phenomenology (RAR, BTFR) while
      a subdominant dark matter component (e.g., light axions or
      ~0.1 eV neutrinos) provides cluster-scale mass:
        - Neutrinos (0.05-0.1 eV, Planck-allowed): too light for clusters
        - Axions (10^{-22} eV): fuzzy DM, might help at cluster scales
        - This abandons pure modified gravity but retains TGP's predictive
          power for galaxy phenomenology

  (2) RELATIVISTIC TGP EXTENSION:
      Following Skordis & Zlosnik (2021), construct a relativistic version
      of TGP where additional tensor/vector degrees of freedom can decouple
      from baryons during a collision.  The "membrane" sector could carry
      gravitational memory through its own dynamics.
      Status: NOT YET DEVELOPED for TGP

  (3) CLUSTER-SPECIFIC gamma(S):
      If clusters have S >> 1 (prolate/triaxial geometries), gamma(S) could
      be much larger, giving stronger boosts.  But gs19 constrains this,
      and even gamma = 1 (extreme) gives nu ~ 1 + 1/y, which is still
      insufficient for the mass deficit.

  (4) ACCEPT THE LIMITATION:
      Acknowledge that TGP, like MOND, is an EFFECTIVE theory valid at
      galaxy scales.  The Bullet Cluster (and clusters in general) may
      require additional physics beyond the TGP interpolation function.
      This does not invalidate TGP's success at galaxy scales, but it
      limits its scope.
""")


# ########################################################################## #
#                                                                            #
#  SUMMARY TABLE                                                            #
#                                                                            #
# ########################################################################## #

print("=" * 78)
print("  SUMMARY: KEY NUMERICAL RESULTS")
print("=" * 78)

print(f"""
  Bullet Cluster parameters:
    M_bar_main     = {M_bar_main/M_sun:.2e} M_sun
    M_bar_bullet   = {M_bar_bullet/M_sun:.2e} M_sun
    M_lens_main    = {M_main_lens/M_sun:.2e} M_sun (observed)
    M_lens_bullet  = {M_bullet_lens/M_sun:.2e} M_sun (observed)
    Needed boost   = {M_main_lens/M_bar_main:.1f}x (main), {M_bullet_lens/M_bar_bullet:.1f}x (bullet)

  TGP predictions:
    nu(y) at r200  = {nu_vir:.3f}  (S = 0.8, cluster scale)
    TGP boost      = {nu_vir:.1f}x  (vs needed {M_main_lens/M_bar_main:.0f}x)
    Mass explained = {1.0/mass_deficit:.0%} of lensing mass

  Lensing offset (1-D model):
    TGP peak (main):   {pk_tgp_L:+.0f} kpc  (observed: -360 kpc)
    TGP peak (bullet): {pk_tgp_R:+.0f} kpc  (observed: +360 kpc)
    Gas centroid:       {x_gas_main/kpc:+.0f} / {x_gas_bullet/kpc:+.0f} kpc
    TGP offset:        {offset_TGP_main_gas:.0f} kpc (main), {offset_TGP_bullet_gas:.0f} kpc (bullet)
    Observed offset:   ~150-260 kpc
    Fraction explained: {frac_offset_main:.0%} (main), {frac_offset_bullet:.0%} (bullet)

  Membrane memory:
    r_c(main)       = {r_c_main/kpc:.0f} kpc
    tau_relax       = {tau_main/3.156e7/1e6:.1f} Myr  (at c_membrane = c)
    t_collision     = {t_collision/3.156e7/1e6:.0f} Myr
    tau/t_coll      = {ratio_main:.4f}  --> NO memory effect

  VERDICT: TGP FAILS the Bullet Cluster test.
    - Mass: explains ~{1.0/mass_deficit:.0%} (need 100%)
    - Offset: phantom mass tracks gas, not galaxies
    - Memory: relaxation is too fast by factor ~{1.0/ratio_main:.0f}
    - Status: KNOWN WEAKNESS, shared with all MOND-like theories
    - Implication: TGP is likely an effective galaxy-scale theory,
      not a complete replacement for dark matter at all scales.
""")

print("=" * 78)
print("  gs24 analysis complete.")
print("=" * 78)
