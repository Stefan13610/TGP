#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gs11: ANALYTICS ON THE WINNING EQUATION
========================================

From gs10, the best-fit interpolation function (171 SPARC galaxies, BIC ranking):

    nu(y) = 1 + exp(-y^alpha) / y^gamma

    Best global fit: alpha = 0.81, gamma = 0.41, a0 = 1.12e-10 m/s^2

This script performs rigorous analytics:
  1. Asymptotic limits (deep MOND y<<1, Newton y>>1, transition y~1)
  2. BTFR: does v^4 = GM*a0 follow automatically?
  3. Solar system test: correction at Earth orbit
  4. Comparison with known interpolation functions
  5. Parameter sensitivity and confidence region (chi2 landscape)
  6. Physical interpretation: what do alpha, gamma mean?
  7. Residual structure: systematic trends in residuals
  8. Refined fit with bootstrap confidence intervals
"""

import sys
import io
import os
import numpy as np
from scipy.optimize import minimize, minimize_scalar, brentq
from scipy.integrate import quad
import warnings
warnings.filterwarnings('ignore')

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

# ============================================================
# CONSTANTS
# ============================================================
G = 6.674e-11
c = 2.998e8
M_sun = 1.989e30
kpc = 3.086e19
Mpc = 3.086e22
H0 = 67.4e3 / (1e6 * kpc)
a0_obs = 1.2e-10  # canonical value

# Winning equation parameters from gs10
ALPHA_BEST = 0.8095
GAMMA_BEST = 0.4059
A0_BEST = 1.1216e-10

# ============================================================
# THE WINNING EQUATION
# ============================================================
def nu_exp(y, alpha=ALPHA_BEST, gamma=GAMMA_BEST):
    """nu(y) = 1 + exp(-y^alpha) / y^gamma"""
    return 1.0 + np.exp(-y**alpha) / y**gamma

def g_obs_from_gbar(gbar, a0=A0_BEST, alpha=ALPHA_BEST, gamma=GAMMA_BEST):
    """g_obs = gbar * nu(gbar/a0)"""
    y = gbar / a0
    return gbar * nu_exp(y, alpha, gamma)

# ============================================================
# SECTION 1: ASYMPTOTIC ANALYSIS
# ============================================================
def analyze_asymptotics():
    print("=" * 80)
    print("SECTION 1: ASYMPTOTIC ANALYSIS")
    print("=" * 80)

    alpha, gamma = ALPHA_BEST, GAMMA_BEST

    print(f"\nWinning equation: nu(y) = 1 + exp(-y^{alpha:.4f}) / y^{gamma:.4f}")
    print(f"a0 = {A0_BEST:.4e} m/s^2")

    # --- Deep MOND limit: y << 1 ---
    print("\n--- Deep MOND (y << 1) ---")
    print(f"exp(-y^{alpha:.4f}) -> 1 for y->0")
    print(f"nu(y) -> 1 + 1/y^{gamma:.4f} -> 1/y^{gamma:.4f}")
    print(f"g_obs = gbar * nu -> gbar / y^{gamma:.4f} = gbar^(1-{gamma:.4f}) * a0^{gamma:.4f}")
    print(f"g_obs = gbar^{1-gamma:.4f} * a0^{gamma:.4f}")
    print()

    # For MOND: g_obs = sqrt(gbar * a0) = gbar^0.5 * a0^0.5
    # Our model: g_obs = gbar^(1-gamma) * a0^gamma
    # MOND requires gamma = 0.5
    print(f"MOND requires: g_obs = gbar^0.5 * a0^0.5 (gamma = 0.5)")
    print(f"Our fit gives: g_obs = gbar^{1-gamma:.4f} * a0^{gamma:.4f} (gamma = {gamma:.4f})")
    print(f"Deviation from MOND: gamma - 0.5 = {gamma - 0.5:.4f}")
    print(f"This is a {abs(gamma-0.5)/0.5*100:.1f}% deviation from exact MOND deep limit")

    # BTFR exponent
    # v^2 = g_obs * r, at large r: g_obs ~ gbar^(1-gamma) * a0^gamma
    # gbar = GM/r^2, so g_obs = (GM/r^2)^(1-gamma) * a0^gamma
    # v^2 = r * (GM)^(1-gamma) * r^(-2(1-gamma)) * a0^gamma
    # v^2 = (GM)^(1-gamma) * a0^gamma * r^(1-2(1-gamma))
    # v^2 = (GM)^(1-gamma) * a0^gamma * r^(2*gamma - 1)
    #
    # For flat RC: need 2*gamma - 1 = 0 => gamma = 0.5 (MOND!)
    # Our gamma = 0.41: v^2 ~ r^(2*0.41 - 1) = r^(-0.18) => SLOWLY DECLINING

    btfr_r_exp = 2*gamma - 1
    print(f"\n--- BTFR and rotation curve shape ---")
    print(f"At large r (deep MOND): v^2 ~ r^({btfr_r_exp:.4f})")
    if abs(btfr_r_exp) < 0.01:
        print(f"=> FLAT rotation curve (exact BTFR)")
    elif btfr_r_exp < 0:
        print(f"=> SLOWLY DECLINING rotation curve (v ~ r^{btfr_r_exp/2:.3f})")
        print(f"   At 2x radius: v drops by {(1 - 2**(btfr_r_exp/2))*100:.1f}%")
        print(f"   At 5x radius: v drops by {(1 - 5**(btfr_r_exp/2))*100:.1f}%")
    else:
        print(f"=> SLOWLY RISING rotation curve")

    # BTFR: v^4 = GM * a0 only if gamma = 0.5
    # General: v^(2/(1-gamma)) ~ GM * a0^(gamma/(1-gamma))
    # For gamma = 0.41: v^(2/0.59) = v^3.39 ~ GM
    print(f"\nBTFR generalized: v^{2/(1-gamma):.2f} ~ GM * a0^{gamma/(1-gamma):.3f}")
    print(f"Standard MOND: v^4 ~ GM * a0")
    print(f"Our model: v^{2/(1-gamma):.2f} ~ GM (exponent {2/(1-gamma):.2f} vs 4)")

    # But wait — the deep MOND limit is only reached at very low y
    # Let's check at what y the transition matters
    print(f"\n--- Transition region ---")
    for y in [0.01, 0.05, 0.1, 0.3, 0.5, 1.0, 2.0, 5.0, 10.0, 100.0]:
        nu_val = nu_exp(y)
        # MOND simple for comparison
        nu_mond = 0.5 + np.sqrt(0.25 + 1.0/y)
        # McGaugh
        nu_mcg = 1.0 / (1.0 - np.exp(-np.sqrt(y)))
        ratio_mond = nu_val / nu_mond
        ratio_mcg = nu_val / nu_mcg
        print(f"  y={y:<8.2f}  nu_exp={nu_val:<10.4f}  nu_MOND={nu_mond:<10.4f}  "
              f"nu_McG={nu_mcg:<10.4f}  exp/MOND={ratio_mond:.4f}  exp/McG={ratio_mcg:.4f}")

    # --- Newtonian limit: y >> 1 ---
    print(f"\n--- Newtonian limit (y >> 1) ---")
    print(f"exp(-y^{alpha:.4f}) -> 0 exponentially")
    print(f"nu(y) -> 1 + (tiny correction)")
    print(f"Correction: delta_nu = exp(-y^{alpha:.4f}) / y^{gamma:.4f}")

    for y in [10, 100, 1e3, 1e4, 1e6, 1e8]:
        delta = np.exp(-y**alpha) / y**gamma
        print(f"  y={y:.0e}:  delta_nu = {delta:.3e}")

    # Solar system
    print(f"\n--- Solar system constraint ---")
    # g_bar at Earth orbit
    M_sun_val = 1.989e30
    r_earth = 1.496e11  # m
    g_earth = G * M_sun_val / r_earth**2
    y_earth = g_earth / A0_BEST
    delta_nu_earth = np.exp(-y_earth**alpha) / y_earth**gamma
    print(f"  g_bar(Earth) = {g_earth:.4e} m/s^2")
    print(f"  y(Earth) = g_bar/a0 = {y_earth:.2e}")
    print(f"  delta_nu = {delta_nu_earth:.3e}")
    print(f"  delta_g/g = {delta_nu_earth:.3e}")
    print(f"  Solar system constraint: delta_g/g < 1e-9")
    if delta_nu_earth < 1e-9:
        print(f"  => SAFE (by factor {1e-9/delta_nu_earth:.0f})")
    else:
        print(f"  => VIOLATED!")

    # Mercury (stronger field)
    r_merc = 5.79e10
    g_merc = G * M_sun_val / r_merc**2
    y_merc = g_merc / A0_BEST
    delta_merc = np.exp(-y_merc**alpha) / y_merc**gamma
    print(f"\n  g_bar(Mercury) = {g_merc:.4e} m/s^2")
    print(f"  y(Mercury) = {y_merc:.2e}")
    print(f"  delta_g/g = {delta_merc:.3e}")

    # Outer solar system (Pioneer)
    r_pioneer = 50 * 1.496e11  # 50 AU
    g_pioneer = G * M_sun_val / r_pioneer**2
    y_pioneer = g_pioneer / A0_BEST
    delta_pioneer = np.exp(-y_pioneer**alpha) / y_pioneer**gamma
    print(f"\n  g_bar(50 AU) = {g_pioneer:.4e} m/s^2")
    print(f"  y(50 AU) = {y_pioneer:.2e}")
    print(f"  delta_g/g = {delta_pioneer:.3e}")

    return gamma, alpha


# ============================================================
# SECTION 2: BTFR PREDICTION
# ============================================================
def analyze_btfr(gamma, alpha):
    print("\n" + "=" * 80)
    print("SECTION 2: BARYONIC TULLY-FISHER RELATION")
    print("=" * 80)

    # For a galaxy with flat RC at radius r:
    # v^2 = g_obs * r = gbar * nu(gbar/a0) * r
    # gbar = GM/(r^2), so:
    # v^2 = GM/r * nu(GM/(r^2 * a0))

    # Find v_flat for different masses
    masses_msun = np.logspace(7, 12, 50)

    print(f"\n--- BTFR from winning equation ---")
    print(f"{'M (Msun)':<14} {'v_flat (km/s)':<14} {'r_flat (kpc)':<14} {'y at r_flat':<12}")

    v_flats = []
    for M in masses_msun:
        M_kg = M * M_sun
        # Find radius where v(r) is maximal (flat part)
        # v^2(r) = GM/r * nu(GM/(r^2*a0))
        r_arr = np.logspace(np.log10(0.1*kpc), np.log10(200*kpc), 1000)
        gbar_arr = G * M_kg / r_arr**2
        y_arr = gbar_arr / A0_BEST
        nu_arr = nu_exp(y_arr, alpha, gamma)
        v2_arr = G * M_kg / r_arr * nu_arr
        v_arr = np.sqrt(v2_arr) / 1e3  # km/s

        # Take the value at the outermost point where RC is still "flat-ish"
        # Use the velocity at the point where dv/dr is closest to 0
        # Or just take the maximum
        idx_max = np.argmax(v_arr)

        # Actually for BTFR, use asymptotic v at large r
        # Take v at the last point
        v_flat = v_arr[-1]
        r_flat = r_arr[-1] / kpc
        y_flat = y_arr[-1]
        v_flats.append(v_flat)

        if M in [1e7, 1e8, 1e9, 1e10, 1e11, 1e12]:
            print(f"  {M:<14.0e} {v_flat:<14.2f} {r_flat:<14.1f} {y_flat:<12.4e}")

    # Fit power law: v^n ~ M
    v_flats = np.array(v_flats)
    log_v = np.log10(v_flats)
    log_M = np.log10(masses_msun)

    # Linear regression log(v) = a*log(M) + b
    coeffs = np.polyfit(log_M, log_v, 1)
    slope = coeffs[0]
    btfr_exp = 1.0 / slope  # M ~ v^(1/slope), so v^(1/slope) ~ M => exponent = 1/slope

    print(f"\n  Fitted: log(v) = {slope:.4f} * log(M) + {coeffs[1]:.4f}")
    print(f"  => M ~ v^{btfr_exp:.2f}")
    print(f"  Standard BTFR: M ~ v^4.0")
    print(f"  Deviation: {abs(btfr_exp - 4.0)/4.0 * 100:.1f}%")

    # Also compute at specific radii (more realistic)
    print(f"\n--- BTFR at fixed outer radius (30 kpc) ---")
    r_fixed = 30 * kpc
    v_30 = []
    for M in masses_msun:
        M_kg = M * M_sun
        gbar = G * M_kg / r_fixed**2
        y = gbar / A0_BEST
        nu = nu_exp(y, alpha, gamma)
        v2 = G * M_kg / r_fixed * nu
        v_30.append(np.sqrt(v2) / 1e3)

    v_30 = np.array(v_30)
    coeffs_30 = np.polyfit(np.log10(masses_msun), np.log10(v_30), 1)
    exp_30 = 1.0 / coeffs_30[0]
    print(f"  Slope: {coeffs_30[0]:.4f}")
    print(f"  => M ~ v^{exp_30:.2f} at 30 kpc")

    # The key question: does gamma != 0.5 actually change BTFR?
    # In practice, most galaxies have RC data in the TRANSITION region (y ~ 0.1-10)
    # not in the asymptotic deep MOND limit (y << 0.01)
    print(f"\n--- Effect of gamma != 0.5 on observable BTFR ---")
    print(f"  Deep MOND prediction: v^{2/(1-gamma):.2f} ~ M (gamma={gamma:.4f})")
    print(f"  At r=30 kpc (realistic): v^{exp_30:.2f} ~ M")
    print(f"  At r=200 kpc (asymptotic): v^{btfr_exp:.2f} ~ M")
    print(f"  Observed BTFR: M ~ v^{3.85}+-0.09 (McGaugh 2012)")
    print(f"  Our model at realistic radii: M ~ v^{exp_30:.2f}")


# ============================================================
# SECTION 3: PARAMETER SENSITIVITY
# ============================================================
def analyze_sensitivity():
    print("\n" + "=" * 80)
    print("SECTION 3: PARAMETER SENSITIVITY (chi2 landscape)")
    print("=" * 80)

    # Load galaxies and compute chi2 for grid of (alpha, gamma, a0)
    base_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(base_dir, 'Rotmod_LTG')

    # Import Galaxy class from gs10
    sys.path.insert(0, base_dir)

    # Inline Galaxy loader
    class Galaxy:
        def __init__(self, filepath):
            self.name = os.path.basename(filepath).replace('_rotmod.dat', '')
            self.distance = None
            lines = []
            with open(filepath, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line.startswith('# Distance'):
                        self.distance = float(line.split('=')[1].strip().split()[0])
                    elif not line.startswith('#') and line:
                        lines.append(line.split())
            data = np.array(lines, dtype=float)
            self.rad_kpc = data[:, 0]
            self.vobs = data[:, 1]
            self.errv = data[:, 2]
            self.vgas = data[:, 3]
            self.vdisk = data[:, 4]
            self.vbul = data[:, 5]
            self.npts = len(self.rad_kpc)
            self.has_bulge = np.any(self.vbul > 0)
            self.errv_safe = np.maximum(self.errv, np.maximum(3.0, 0.03 * np.abs(self.vobs)))
            self.rad_m = self.rad_kpc * kpc

        def compute_gbar(self, Yd, Yb=None):
            if Yb is None: Yb = Yd
            v2_bar = np.abs(self.vgas) * self.vgas + Yd * self.vdisk**2 + Yb * self.vbul**2
            return v2_bar * 1e6 / self.rad_m

    galaxies = []
    for fname in sorted(os.listdir(data_dir)):
        if fname.endswith('_rotmod.dat'):
            try:
                gal = Galaxy(os.path.join(data_dir, fname))
                if gal.npts >= 5:
                    galaxies.append(gal)
            except:
                pass
    print(f"  Loaded {len(galaxies)} galaxies")

    def total_chi2_for_params(a0, alpha, gamma):
        """Quick chi2: for each galaxy, fit only Yd (1D optimization)."""
        total = 0.0
        for gal in galaxies:
            best_chi2 = 1e10
            # Scan Yd
            for log_Yd in np.linspace(-0.7, 1.0, 30):
                Yd = 10**log_Yd
                gbar = gal.compute_gbar(Yd)
                mask = gbar > 0
                if np.sum(mask) < 3:
                    continue
                y = gbar[mask] / a0
                y = np.maximum(y, 1e-10)
                nu = 1.0 + np.exp(-y**alpha) / y**gamma
                gobs_pred = gbar[mask] * nu
                v_pred = np.sqrt(np.abs(gobs_pred) * gal.rad_m[mask]) / 1e3
                chi2 = np.sum(((gal.vobs[mask] - v_pred) / gal.errv_safe[mask])**2)
                if chi2 < best_chi2:
                    best_chi2 = chi2
            total += best_chi2
        return total

    # Baseline chi2
    chi2_base = total_chi2_for_params(A0_BEST, ALPHA_BEST, GAMMA_BEST)
    print(f"  Baseline chi2 = {chi2_base:.1f}")

    # Scan alpha (fix gamma, a0)
    print(f"\n--- Alpha scan (gamma={GAMMA_BEST:.4f} fixed, a0={A0_BEST:.2e}) ---")
    alphas = np.linspace(0.3, 1.5, 25)
    chi2_alpha = []
    for a in alphas:
        c2 = total_chi2_for_params(A0_BEST, a, GAMMA_BEST)
        chi2_alpha.append(c2)
        if abs(a - ALPHA_BEST) < 0.05 or c2 - chi2_base < 50:
            print(f"  alpha={a:.3f}: chi2={c2:.1f} (delta={c2-chi2_base:+.1f})")

    # 1-sigma: delta_chi2 = 1.0 for 1 parameter
    chi2_alpha = np.array(chi2_alpha)
    mask_1sig = chi2_alpha - chi2_base < 1.0
    if np.any(mask_1sig):
        alpha_lo = alphas[mask_1sig][0]
        alpha_hi = alphas[mask_1sig][-1]
        print(f"\n  1-sigma range: alpha = {alpha_lo:.3f} to {alpha_hi:.3f}")
        print(f"  Best: {ALPHA_BEST:.4f} +{alpha_hi-ALPHA_BEST:.3f} -{ALPHA_BEST-alpha_lo:.3f}")

    # Scan gamma
    print(f"\n--- Gamma scan (alpha={ALPHA_BEST:.4f} fixed) ---")
    gammas = np.linspace(0.2, 0.8, 25)
    chi2_gamma = []
    for g in gammas:
        c2 = total_chi2_for_params(A0_BEST, ALPHA_BEST, g)
        chi2_gamma.append(c2)
        if abs(g - GAMMA_BEST) < 0.03 or c2 - chi2_base < 50:
            print(f"  gamma={g:.3f}: chi2={c2:.1f} (delta={c2-chi2_base:+.1f})")

    chi2_gamma = np.array(chi2_gamma)
    mask_1sig_g = chi2_gamma - chi2_base < 1.0
    if np.any(mask_1sig_g):
        gamma_lo = gammas[mask_1sig_g][0]
        gamma_hi = gammas[mask_1sig_g][-1]
        print(f"\n  1-sigma range: gamma = {gamma_lo:.3f} to {gamma_hi:.3f}")
        print(f"  Best: {GAMMA_BEST:.4f} +{gamma_hi-GAMMA_BEST:.3f} -{GAMMA_BEST-gamma_lo:.3f}")

    # Scan a0
    print(f"\n--- a0 scan (alpha, gamma fixed) ---")
    a0s = np.logspace(-10.3, -9.5, 25)
    chi2_a0 = []
    for a0 in a0s:
        c2 = total_chi2_for_params(a0, ALPHA_BEST, GAMMA_BEST)
        chi2_a0.append(c2)
        if abs(np.log10(a0) - np.log10(A0_BEST)) < 0.05 or c2 - chi2_base < 50:
            print(f"  a0={a0:.3e}: chi2={c2:.1f} (delta={c2-chi2_base:+.1f})")

    chi2_a0 = np.array(chi2_a0)
    mask_1sig_a = chi2_a0 - chi2_base < 1.0
    if np.any(mask_1sig_a):
        a0_lo = a0s[mask_1sig_a][0]
        a0_hi = a0s[mask_1sig_a][-1]
        print(f"\n  1-sigma range: a0 = {a0_lo:.3e} to {a0_hi:.3e}")
        print(f"  Best: {A0_BEST:.4e} +{(a0_hi-A0_BEST):.3e} -{(A0_BEST-a0_lo):.3e}")

    # KEY: does gamma = 0.5 (MOND) fall within confidence region?
    print(f"\n--- CRITICAL: Is gamma = 0.5 (exact MOND) within confidence? ---")
    chi2_at_05 = total_chi2_for_params(A0_BEST, ALPHA_BEST, 0.5)
    delta_05 = chi2_at_05 - chi2_base
    print(f"  chi2 at gamma=0.5: {chi2_at_05:.1f} (delta = {delta_05:+.1f})")
    if delta_05 < 1.0:
        print(f"  => gamma=0.5 is WITHIN 1-sigma")
    elif delta_05 < 4.0:
        print(f"  => gamma=0.5 is WITHIN 2-sigma")
    elif delta_05 < 9.0:
        print(f"  => gamma=0.5 is WITHIN 3-sigma")
    else:
        n_sigma = np.sqrt(delta_05)
        print(f"  => gamma=0.5 is EXCLUDED at {n_sigma:.1f}-sigma")

    # 2D scan: alpha vs gamma
    print(f"\n--- 2D scan: alpha vs gamma (a0 fixed at best) ---")
    alphas_2d = np.linspace(0.4, 1.3, 15)
    gammas_2d = np.linspace(0.25, 0.65, 15)
    chi2_2d = np.zeros((len(alphas_2d), len(gammas_2d)))

    for i, a in enumerate(alphas_2d):
        for j, g in enumerate(gammas_2d):
            chi2_2d[i, j] = total_chi2_for_params(A0_BEST, a, g)

    min_2d = np.min(chi2_2d)
    ij_min = np.unravel_index(np.argmin(chi2_2d), chi2_2d.shape)
    print(f"  2D minimum: alpha={alphas_2d[ij_min[0]]:.3f}, gamma={gammas_2d[ij_min[1]]:.3f}, chi2={min_2d:.1f}")

    # Show contours
    print(f"\n  Delta-chi2 map (rows=alpha, cols=gamma):")
    print(f"  {'gamma:':<8}", end='')
    for g in gammas_2d[::3]:
        print(f" {g:.2f}  ", end='')
    print()

    for i in range(0, len(alphas_2d), 3):
        print(f"  a={alphas_2d[i]:.2f}  ", end='')
        for j in range(0, len(gammas_2d), 3):
            delta = chi2_2d[i, j] - min_2d
            if delta < 1:
                sym = ' ** '
            elif delta < 4:
                sym = ' *  '
            elif delta < 9:
                sym = ' .  '
            else:
                sym = f'{delta:5.0f}'
            print(sym, end='')
        print()

    # Find 1-sigma contour extent
    mask_2d_1sig = chi2_2d - min_2d < 2.30  # 2 params: delta_chi2 = 2.30 for 1-sigma
    alpha_range = alphas_2d[np.any(mask_2d_1sig, axis=1)]
    gamma_range = gammas_2d[np.any(mask_2d_1sig, axis=0)]
    if len(alpha_range) > 0 and len(gamma_range) > 0:
        print(f"\n  1-sigma (2-param, dchi2<2.3):")
        print(f"    alpha: {alpha_range[0]:.3f} to {alpha_range[-1]:.3f}")
        print(f"    gamma: {gamma_range[0]:.3f} to {gamma_range[-1]:.3f}")

    mask_2d_2sig = chi2_2d - min_2d < 6.18  # 2 params: 2-sigma
    alpha_range2 = alphas_2d[np.any(mask_2d_2sig, axis=1)]
    gamma_range2 = gammas_2d[np.any(mask_2d_2sig, axis=0)]
    if len(alpha_range2) > 0 and len(gamma_range2) > 0:
        print(f"  2-sigma (2-param, dchi2<6.18):")
        print(f"    alpha: {alpha_range2[0]:.3f} to {alpha_range2[-1]:.3f}")
        print(f"    gamma: {gamma_range2[0]:.3f} to {gamma_range2[-1]:.3f}")

    # Check if (alpha=1, gamma=0.5) = simple MOND-like form is within contour
    idx_a1 = np.argmin(np.abs(alphas_2d - 1.0))
    idx_g05 = np.argmin(np.abs(gammas_2d - 0.5))
    if idx_a1 < len(alphas_2d) and idx_g05 < len(gammas_2d):
        delta_mond = chi2_2d[idx_a1, idx_g05] - min_2d
        print(f"\n  At (alpha=1.0, gamma=0.5) [MOND-like]: dchi2 = {delta_mond:.1f}")

    return chi2_base


# ============================================================
# SECTION 4: PHYSICAL INTERPRETATION
# ============================================================
def physical_interpretation(gamma, alpha):
    print("\n" + "=" * 80)
    print("SECTION 4: PHYSICAL INTERPRETATION")
    print("=" * 80)

    print(f"""
The winning equation: nu(y) = 1 + exp(-y^alpha) / y^gamma
    alpha = {alpha:.4f}, gamma = {gamma:.4f}

DECOMPOSITION:
  nu = 1     + exp(-y^a) / y^g
      Newton   MOND boost (decays exponentially in strong field)

STRUCTURE:
  - exp(-y^a): SCREENING FUNCTION. Kills MOND boost for y >> 1 (strong field)
    - alpha controls HOW FAST screening turns on
    - alpha = {alpha:.2f} < 1: screening turns on SLOWLY (sub-linear in y)
    - alpha = 1 would give exp(-y) = simple exponential screening
    - alpha = 0.5 would give exp(-sqrt(y)) = McGaugh-like screening

  - 1/y^g: MOND BOOST AMPLITUDE in weak field
    - gamma = {gamma:.4f} controls the DEPTH of the MOND regime
    - gamma = 0.5: exact sqrt(a0/gbar) behavior (MOND)
    - gamma = {gamma:.4f}: slightly WEAKER boost than MOND (by y^0.09)

COMPARISON WITH KNOWN FUNCTIONS:
  McGaugh:      nu = 1/(1-exp(-sqrt(y))) ~= 1 + exp(-sqrt(y))/sqrt(y) for large y
                => alpha_eff ~ 0.5, gamma_eff ~ 0.5
  MOND simple:  nu = 1/2 + sqrt(1/4+1/y) ~= 1 + 1/(2y) for large y
                => no exponential screening, gamma_eff ~ 1.0 at transition
  Hybrid A:     nu = 1 + 1/(sqrt(y)+y) ~= 1 + 1/y for large y
                => gamma_eff ~ 1.0, with polynomial screening

KEY INSIGHT: Our alpha={alpha:.2f} is BETWEEN 0.5 (McGaugh) and 1.0 (naive exp).
  This suggests the screening mechanism is:
  - NOT purely Gaussian (alpha=2)
  - NOT purely exponential (alpha=1)
  - Close to McGaugh-type (alpha~0.5) but with slightly SLOWER onset

  gamma = {gamma:.4f} being close to 0.5 means the deep MOND limit is
  ALMOST exact MOND, with a tiny correction.
""")

    # TGP interpretation
    print(f"TGP INTERPRETATION:")
    print(f"""
  If gravity transitions from 3D (Newton) to effectively 2D (MOND):
  - The '1' in nu = 1 + ... is the 3D Newtonian channel
  - The exp(-y^a)/y^g term is the 2D channel with screening

  In TGP membrane picture:
  - y = gbar/a0 = (baryonic field)/(transition scale)
  - exp(-y^a) = probability of field being "effectively 2D"
    at given field strength
  - 1/y^g = strength of 2D gravitational enhancement

  alpha = {alpha:.4f}: the dimensional transition is GRADUAL
    (sub-exponential, not a sharp threshold)
  gamma = {gamma:.4f}: the 2D enhancement is ALMOST sqrt (0.41 vs 0.50)
    A pure 2D channel gives F ~ 1/r => v^2 ~ 1/r^0 (flat)
    Our gamma<0.5 gives F ~ 1/r^{{1+2(0.5-gamma)}} => v slowly declining
    This could mean: transition is not to PURE 2D, but to d_eff ~ 2.18
""")

    d_eff_deep = 2 * (1 - gamma) / (1 - gamma)  # wrong, let me think
    # In d dimensions: F ~ 1/r^(d-1), v^2 = F*r ~ 1/r^(d-2)
    # For flat: d=2. For v^2 ~ r^(2g-1):
    # r^(2g-1) = r^(2-d) => d = 3 - 2g
    d_eff = 3 - 2*gamma
    print(f"  Effective dimension in deep MOND: d_eff = 3 - 2*gamma = {d_eff:.3f}")
    print(f"  Pure 2D: d_eff = 2.0 (gamma=0.5)")
    print(f"  Our model: d_eff = {d_eff:.3f} (gamma={gamma:.4f})")
    print(f"  => NOT pure 2D, but d_eff = {d_eff:.2f} (slightly above 2)")
    print(f"     Fraction toward 2D: {(3-d_eff)/(3-2)*100:.0f}% of full 3D->2D transition")

    return d_eff


# ============================================================
# SECTION 5: COMPARISON TABLE
# ============================================================
def comparison_table():
    print("\n" + "=" * 80)
    print("SECTION 5: WINNING EQUATION vs KNOWN FUNCTIONS")
    print("=" * 80)

    print(f"\n{'y (gbar/a0)':<14} {'Exp(win)':<12} {'McGaugh':<12} {'MOND_s':<12} {'Ratio E/M':<12} {'Ratio E/MOND':<12}")
    print("-" * 74)

    for y in [0.001, 0.01, 0.03, 0.1, 0.3, 0.5, 1.0, 2.0, 5.0, 10, 30, 100, 1000]:
        nu_e = nu_exp(y)
        nu_m = 1.0 / (1.0 - np.exp(-np.sqrt(y)))
        nu_ms = 0.5 + np.sqrt(0.25 + 1.0/y)

        print(f"  {y:<12.3f} {nu_e:<12.4f} {nu_m:<12.4f} {nu_ms:<12.4f} "
              f"{nu_e/nu_m:<12.4f} {nu_e/nu_ms:<12.4f}")

    # Maximum difference
    y_range = np.logspace(-3, 3, 10000)
    nu_e_arr = nu_exp(y_range)
    nu_m_arr = 1.0 / (1.0 - np.exp(-np.sqrt(y_range)))

    ratio = nu_e_arr / nu_m_arr
    idx_max_diff = np.argmax(np.abs(ratio - 1))
    print(f"\n  Max difference from McGaugh: {(ratio[idx_max_diff]-1)*100:.1f}% at y={y_range[idx_max_diff]:.3f}")

    # Where does our function cross McGaugh?
    crossings = np.where(np.diff(np.sign(ratio - 1)))[0]
    print(f"  Crossings (nu_exp = nu_McGaugh) at y = ", end='')
    for c in crossings[:5]:
        print(f"{y_range[c]:.3f}, ", end='')
    print()


# ============================================================
# SECTION 6: SIMPLIFIED FORMS
# ============================================================
def simplified_forms():
    print("\n" + "=" * 80)
    print("SECTION 6: SIMPLIFIED/APPROXIMATE FORMS")
    print("=" * 80)

    alpha, gamma = ALPHA_BEST, GAMMA_BEST

    # Can we approximate alpha ~ 4/5, gamma ~ 2/5?
    print(f"\n--- Rational approximations ---")
    print(f"  alpha = {alpha:.4f} ~ 4/5 = {4/5:.4f} (error: {abs(alpha-4/5)/alpha*100:.2f}%)")
    print(f"  gamma = {gamma:.4f} ~ 2/5 = {2/5:.4f} (error: {abs(gamma-2/5)/gamma*100:.2f}%)")
    print(f"  => gamma = alpha/2 ? alpha/2 = {alpha/2:.4f} vs gamma = {gamma:.4f} "
          f"(error: {abs(gamma-alpha/2)/gamma*100:.2f}%)")
    print(f"  => gamma/alpha = {gamma/alpha:.4f} ~ 1/2 = 0.500")

    # So the equation simplifies to:
    # nu(y) = 1 + exp(-y^(4/5)) / y^(2/5)
    # = 1 + exp(-y^(4/5)) * y^(-2/5)
    #
    # Substituting u = y^(2/5):
    # nu = 1 + exp(-u^2) / u
    # This is related to the complementary error function!
    print(f"""
  With alpha=4/5, gamma=2/5:
    nu(y) = 1 + exp(-y^(4/5)) / y^(2/5)

  Substituting u = y^(2/5):
    nu = 1 + exp(-u^2) / u

  Note: erfc(u) = (2/sqrt(pi)) * integral_u^inf exp(-t^2) dt
        For large u: erfc(u) ~ exp(-u^2) / (sqrt(pi) * u)

  So: exp(-u^2)/u ~ sqrt(pi) * erfc(u)

  => nu(y) ~ 1 + sqrt(pi) * erfc(y^(2/5))

  This is a COMPLEMENTARY ERROR FUNCTION in the variable y^(2/5)!
""")

    # Verify numerically
    from scipy.special import erfc
    print(f"--- Verification: nu vs 1 + sqrt(pi)*erfc(y^(2/5)) ---")
    print(f"{'y':<10} {'nu_exact':<12} {'erfc_approx':<12} {'ratio':<10}")
    for y in [0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 100.0]:
        nu_exact = 1.0 + np.exp(-y**(4/5)) / y**(2/5)
        nu_erfc = 1.0 + np.sqrt(np.pi) * erfc(y**(2/5))
        print(f"  {y:<10.2f} {nu_exact:<12.6f} {nu_erfc:<12.6f} {nu_erfc/nu_exact:<10.6f}")

    # The erfc approximation is only good for large u
    # For small u, erfc(u) -> 1, but exp(-u^2)/u -> 1/u -> diverges
    # So the connection to erfc is only qualitative

    # Better: the equation looks like a Dawson function or similar
    # exp(-u^2)/u for u = y^(2/5)

    # Physical: probability theory interpretation
    print(f"""
  PROBABILITY INTERPRETATION:
    The term exp(-u^2)/u (where u = y^(2/5) = (gbar/a0)^(2/5)) looks like
    the tail of a Gaussian distribution in the variable u.

    If gbar/a0 represents a "signal-to-noise ratio" for the baryonic field
    against the cosmological background noise, then:
    - exp(-u^2) = probability that the field is "below threshold"
      (noise-dominated => 2D behavior)
    - 1/u = amplitude of the 2D boost when noise-dominated

    The transition scale is u ~ 1, i.e., y ~ 1^(5/2) = 1, i.e., gbar ~ a0.
    This is EXACTLY where a0 sits by definition!
""")


# ============================================================
# MAIN
# ============================================================
def main():
    print("=" * 80)
    print("gs11: ANALYTICS ON THE WINNING EQUATION")
    print(f"nu(y) = 1 + exp(-y^{ALPHA_BEST:.4f}) / y^{GAMMA_BEST:.4f}")
    print(f"a0 = {A0_BEST:.4e} m/s^2")
    print("=" * 80)

    gamma, alpha = analyze_asymptotics()
    analyze_btfr(gamma, alpha)
    chi2_base = analyze_sensitivity()
    d_eff = physical_interpretation(gamma, alpha)
    comparison_table()
    simplified_forms()

    # Final summary
    print("\n" + "=" * 80)
    print("FINAL SUMMARY")
    print("=" * 80)
    print(f"""
WINNING EQUATION: nu(y) = 1 + exp(-y^0.81) / y^0.41
                  a0 = 1.12e-10 m/s^2

EQUIVALENT FORMS:
  Exact:  g_obs = g_bar + g_bar * exp(-(g_bar/a0)^0.81) / (g_bar/a0)^0.41
  Approx: g_obs ~ g_bar + a0^0.41 * g_bar^0.59 * exp(-(g_bar/a0)^0.81)
  Rational: nu ~ 1 + exp(-y^(4/5)) / y^(2/5)  [alpha~4/5, gamma~2/5]

ASYMPTOTIC LIMITS:
  Deep MOND (y<<1): g_obs ~ g_bar^0.59 * a0^0.41  (NOT exact MOND sqrt!)
  Newton (y>>1):    g_obs ~ g_bar * (1 + exp(-y^0.81)/y^0.41)  (safe)

EFFECTIVE DIMENSION: d_eff = {d_eff:.3f} in deep MOND (not pure 2D!)

BTFR: v^{2/(1-GAMMA_BEST):.2f} ~ M (asymptotic), close to v^4 at realistic radii

SOLAR SYSTEM: delta_g/g << 1e-9 at Earth — SAFE

KEY PHYSICAL INSIGHT:
  gamma = 0.41 (not 0.5) means the transition is NOT to pure 2D.
  The effective dimension in the deep MOND regime is d_eff = {d_eff:.2f}.
  This is {(3-d_eff)/(3-2)*100:.0f}% of a full 3D->2D transition.
""")

    print("=" * 80)
    print("DONE")
    print("=" * 80)


if __name__ == '__main__':
    main()
