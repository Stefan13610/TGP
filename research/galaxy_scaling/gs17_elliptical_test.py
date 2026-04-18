#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gs17: ELLIPTICAL GALAXY TEST — geometry hypothesis
====================================================
gs16 prediction: spheroidal systems should prefer gamma ~ 0.54,
while disk galaxies prefer gamma ~ 0.36.

Ellipticals are spheroidal but galaxy-scale -> CRITICAL TEST:
  If gamma_elliptical ~ 0.54 -> geometry hypothesis confirmed
  If gamma_elliptical ~ 0.36 -> geometry hypothesis FALSIFIED,
     and cluster problem is about MASS, not geometry

Data sources:
  - ATLAS3D: 260 ETGs with M_JAM, sigma_e, R_eff (Cappellari+2013)
  - Lelli+2017: 25 ETGs with full RAR data
  - Chae+2020: 19 E0 slow-rotators from MaNGA
  - Individual well-studied: NGC4486(M87), NGC4472(M49), etc.

Method:
  For each elliptical, at R_eff:
    g_bar = G * M_star / (2 * R_eff^2)  (half-mass within R_eff)
    g_obs = G * M_dyn / R_eff^2 = K * sigma_e^2 / R_eff
    where K ~ 5 (virial coefficient, Cappellari+2006)

  Then: y = g_bar / a0, nu_obs = g_obs / g_bar
  Test which gamma fits best.
"""

import sys
import io
import numpy as np
from scipy.optimize import minimize, minimize_scalar
from scipy.stats import spearmanr
import warnings
warnings.filterwarnings('ignore')

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

G_SI = 6.674e-11
M_sun = 1.989e30
kpc_m = 3.086e19
pc_m = 3.086e16

a0_ref = 1.12e-10

# ============================================================
# MODELS
# ============================================================
def nu_exp(y, alpha=0.8, gamma=0.4):
    y = np.maximum(y, 1e-15)
    return 1.0 + np.exp(-np.power(y, alpha)) / np.power(y, gamma)

def nu_mcgaugh(y):
    y = np.maximum(y, 1e-15)
    return 1.0 / (1.0 - np.exp(-np.sqrt(y)))

def nu_mond(y):
    y = np.maximum(y, 1e-15)
    return 0.5 + np.sqrt(0.25 + 1.0/y)


# ============================================================
# ELLIPTICAL GALAXY DATA
# ============================================================
def build_ellipticals():
    """
    Compile elliptical galaxy data from published sources.

    For each galaxy:
      name, morphology, M_star (Msun), R_eff (kpc), sigma_e (km/s),
      M_dyn_1Re (Msun), source

    M_dyn within 1 R_eff from: M_JAM or virial estimator
      M_dyn = K * sigma_e^2 * R_eff / G
      K ~ 5 for typical ETGs (Cappellari+2006)

    M_star from stellar population synthesis (SPS)
    """
    K_virial = 5.0  # virial coefficient

    ellipticals = []

    # -------------------------------------------------------
    # SOURCE 1: Well-studied individual ellipticals
    # Data from Cappellari+2013 (ATLAS3D XV), Kormendy+2009,
    # Humphrey+2006 (X-ray), Lelli+2017
    # (name, type, M_star_1e10, R_eff_kpc, sigma_e_kms)
    # -------------------------------------------------------
    atlas3d_sample = [
        # Massive slow-rotators (true ellipticals)
        ("NGC4486",  "E0",  55.0,  7.4,  375),   # M87 - cD Virgo
        ("NGC4472",  "E2",  75.0,  8.6,  296),   # M49 - brightest Virgo
        ("NGC4649",  "E2",  40.0,  5.5,  341),   # M60
        ("NGC4636",  "E0",  18.0,  8.6,  209),   # X-ray bright
        ("NGC5846",  "E0",  30.0,  6.9,  248),   # group central
        ("NGC4552",  "E0",  12.0,  2.7,  264),   # M89
        ("NGC4406",  "E3",  25.0,  8.1,  235),   # M86
        ("NGC4374",  "E1",  28.0,  5.3,  296),   # M84
        ("NGC4261",  "E2",  30.0,  5.3,  315),   # radio galaxy
        ("NGC1399",  "E1",  38.0,  4.6,  340),   # Fornax central
        ("NGC1407",  "E0",  32.0,  7.0,  276),   # group central
        ("NGC5813",  "E1",  17.0,  4.6,  237),
        ("NGC5044",  "E0",  14.0,  4.1,  239),   # group central
        ("NGC4278",  "E1",   6.0,  2.1,  252),
        ("NGC3379",  "E1",   5.8,  2.4,  209),   # M105
        ("NGC4125",  "E6",   8.0,  4.3,  228),
        ("NGC1600",  "E3",  42.0,  8.5,  334),   # massive

        # Fast-rotators (more disky, intermediate)
        ("NGC3377",  "E5",   1.5,  1.5,  145),
        ("NGC4473",  "E5",   4.0,  2.4,  193),
        ("NGC4697",  "E6",   7.0,  4.3,  170),
        ("NGC4564",  "E6",   2.2,  1.3,  162),
        ("NGC4660",  "E5",   1.3,  0.8,  192),
        ("NGC3608",  "E2",   4.5,  2.3,  192),
        ("NGC4382",  "S0",  12.0,  5.6,  182),   # M85, S0/E
        ("NGC4526",  "S0",   8.0,  2.8,  222),

        # Lower mass ellipticals
        ("NGC4551",  "E2",   1.0,  1.0,  119),
        ("NGC4387",  "E5",   0.7,  0.8,  108),
        ("NGC4458",  "E0",   1.0,  1.4,  103),
        ("NGC4489",  "E0",   0.6,  1.2,   56),
        ("NGC4478",  "E2",   1.5,  1.0,  141),

        # Compact ellipticals
        ("NGC4486B", "cE",   0.5,  0.3,  170),   # M87 companion
        ("M32",      "cE",   0.3,  0.1,   77),   # M31 companion
    ]

    for name, morph, M_star_1e10, R_eff_kpc, sigma_e in atlas3d_sample:
        M_star = M_star_1e10 * 1e10 * M_sun  # kg
        R_eff = R_eff_kpc * kpc_m  # m
        sigma = sigma_e * 1e3  # m/s

        # Dynamical mass within R_eff
        M_dyn = K_virial * sigma**2 * R_eff / G_SI

        # Accelerations at R_eff
        # g_bar: use half the stellar mass (enclosed within R_eff for Sersic n~4)
        # For Sersic n=4 (de Vaucouleurs), ~half the light is within R_eff
        g_bar = G_SI * (0.5 * M_star) / R_eff**2
        g_obs = G_SI * M_dyn / R_eff**2

        # Also compute at 2*R_eff and 5*R_eff (extrapolated)
        # M_star(<2Re) ~ 0.75 * M_star, M_dyn ~ sigma^2 * 2Re * K/G (approx)
        # For simplicity, just use R_eff

        ellipticals.append({
            'name': name,
            'morph': morph,
            'M_star': M_star_1e10 * 1e10,  # Msun
            'R_eff_kpc': R_eff_kpc,
            'sigma_e': sigma_e,
            'M_dyn_Msun': M_dyn / M_sun,
            'g_bar': g_bar,
            'g_obs': g_obs,
            'y': g_bar / a0_ref,
            'nu_obs': g_obs / g_bar if g_bar > 0 else 1.0,
            'fDM': 1.0 - M_star / (2.0 * M_dyn),  # DM fraction within R_eff
            'is_slow_rotator': morph in ['E0', 'E1', 'E2', 'E3', 'cE'],
        })

    return ellipticals


# ============================================================
# MAIN
# ============================================================
def main():
    print("=" * 80)
    print("gs17: ELLIPTICAL GALAXY TEST")
    print("Do spheroidal galaxies prefer gamma ~ 0.54 (geometry hypothesis)?")
    print("=" * 80)

    ellipticals = build_ellipticals()
    print(f"\nLoaded {len(ellipticals)} elliptical/S0 galaxies")

    # ================================================================
    # SECTION 1: Acceleration regime
    # ================================================================
    print("\n" + "=" * 80)
    print("SECTION 1: ACCELERATION REGIME OF ELLIPTICALS")
    print("=" * 80)

    ys = np.array([e['y'] for e in ellipticals])
    g_bars = np.array([e['g_bar'] for e in ellipticals])
    nus = np.array([e['nu_obs'] for e in ellipticals])

    print(f"\n  g_bar range: {g_bars.min():.2e} to {g_bars.max():.2e} m/s^2")
    print(f"  y = g_bar/a0: {ys.min():.3f} to {ys.max():.1f}")
    print(f"  y median: {np.median(ys):.2f}")
    print(f"  nu_obs range: {nus.min():.2f} to {nus.max():.2f}")
    print(f"  nu_obs median: {np.median(nus):.2f}")

    deep = sum(1 for y in ys if y < 0.1)
    trans = sum(1 for y in ys if 0.1 <= y < 1.0)
    newt = sum(1 for y in ys if 1.0 <= y < 10)
    high = sum(1 for y in ys if y >= 10)
    print(f"\n  Deep MOND (y<0.1):   {deep}")
    print(f"  Transition (0.1-1):  {trans}")
    print(f"  Newtonian (1-10):    {newt}")
    print(f"  High-acc (y>10):     {high}")

    print(f"\n  Compare: disk galaxies median y_outer ~ 0.05")
    print(f"  Ellipticals at R_eff: median y ~ {np.median(ys):.1f}")
    print(f"  => Ellipticals probe HIGHER y (more Newtonian) than disk outskirts!")

    # ================================================================
    # SECTION 2: Model comparison at R_eff
    # ================================================================
    print("\n" + "=" * 80)
    print("SECTION 2: MODEL PREDICTIONS vs OBSERVATIONS AT R_eff")
    print("=" * 80)

    gammas_test = [0.30, 0.36, 0.40, 0.45, 0.50, 0.54, 0.60]

    print(f"\n  {'gamma':<8s} {'med(g_pred/g_obs)':<20s} {'mean ratio':<12s} "
          f"{'scatter(dex)':<14s} {'chi2':<10s}")
    print("  " + "-" * 64)

    best_gamma = 0.4
    best_chi2 = 1e10

    for gamma in gammas_test:
        ratios = []
        chi2 = 0
        for e in ellipticals:
            nu_pred = nu_exp(e['y'], 0.8, gamma)
            g_pred = e['g_bar'] * nu_pred
            ratio = g_pred / e['g_obs']
            ratios.append(ratio)
            chi2 += (np.log10(ratio))**2

        ratios = np.array(ratios)
        scatter = np.std(np.log10(ratios))
        print(f"  {gamma:<8.2f} {np.median(ratios):<20.4f} {np.mean(ratios):<12.4f} "
              f"{scatter:<14.4f} {chi2:<10.3f}")

        if chi2 < best_chi2:
            best_chi2 = chi2
            best_gamma = gamma

    # McGaugh and MOND
    for mname, mfunc in [("McGaugh", nu_mcgaugh), ("MOND", nu_mond)]:
        ratios = []
        chi2 = 0
        for e in ellipticals:
            y_m = e['g_bar'] / 1.2e-10  # MOND a0
            nu_pred = mfunc(y_m)
            g_pred = e['g_bar'] * nu_pred
            ratio = g_pred / e['g_obs']
            ratios.append(ratio)
            chi2 += (np.log10(ratio))**2
        ratios = np.array(ratios)
        scatter = np.std(np.log10(ratios))
        print(f"  {mname:<8s} {np.median(ratios):<20.4f} {np.mean(ratios):<12.4f} "
              f"{scatter:<14.4f} {chi2:<10.3f}")

    print(f"\n  Best gamma from grid: {best_gamma:.2f}")

    # ================================================================
    # SECTION 3: Optimize gamma for ellipticals
    # ================================================================
    print("\n" + "=" * 80)
    print("SECTION 3: BEST-FIT GAMMA FOR ELLIPTICALS")
    print("=" * 80)

    def chi2_gamma(gamma):
        if gamma < 0.01 or gamma > 2.0: return 1e10
        total = 0
        for e in ellipticals:
            nu_pred = nu_exp(e['y'], 0.8, gamma)
            g_pred = e['g_bar'] * nu_pred
            total += (np.log10(g_pred / e['g_obs']))**2
        return total

    def chi2_gamma_a0(params):
        gamma, log_a0 = params
        a0 = 10**log_a0
        if gamma < 0.01 or gamma > 2.0: return 1e10
        total = 0
        for e in ellipticals:
            y = e['g_bar'] / a0
            nu_pred = nu_exp(y, 0.8, gamma)
            g_pred = e['g_bar'] * nu_pred
            total += (np.log10(g_pred / e['g_obs']))**2
        return total

    # Optimize gamma (a0 fixed)
    res1 = minimize_scalar(chi2_gamma, bounds=(0.05, 1.5), method='bounded')
    gamma_best_fixed_a0 = res1.x
    chi2_best1 = res1.fun

    # Optimize (gamma, a0) jointly
    best2 = {'chi2': 1e10}
    for g0 in [0.3, 0.4, 0.5, 0.6, 0.8]:
        for la0 in [-10.2, -10.0, -9.8, -9.5]:
            try:
                res = minimize(chi2_gamma_a0, [g0, la0], method='Nelder-Mead',
                             options={'maxiter': 5000})
                if res.fun < best2['chi2']:
                    best2 = {'chi2': res.fun, 'gamma': res.x[0], 'a0': 10**res.x[1]}
            except: pass

    print(f"\n  Fixed a0={a0_ref:.2e}:")
    print(f"    Best gamma = {gamma_best_fixed_a0:.4f}")
    print(f"    chi2 = {chi2_best1:.4f}")

    print(f"\n  Free (gamma, a0):")
    print(f"    Best gamma = {best2['gamma']:.4f}")
    print(f"    Best a0 = {best2['a0']:.3e} ({best2['a0']/a0_ref:.2f}x galaxy)")
    print(f"    chi2 = {best2['chi2']:.4f}")

    # ================================================================
    # SECTION 4: Slow-rotators vs fast-rotators
    # ================================================================
    print("\n" + "=" * 80)
    print("SECTION 4: SLOW-ROTATORS (true E) vs FAST-ROTATORS (disky E/S0)")
    print("=" * 80)

    slow = [e for e in ellipticals if e['is_slow_rotator']]
    fast = [e for e in ellipticals if not e['is_slow_rotator']]

    print(f"\n  Slow-rotators (true spheroids): {len(slow)}")
    print(f"  Fast-rotators (disky component): {len(fast)}")

    for label, subset in [("Slow-rotators", slow), ("Fast-rotators", fast)]:
        def chi2_sub(gamma):
            if gamma < 0.01 or gamma > 2.0: return 1e10
            total = 0
            for e in subset:
                nu_pred = nu_exp(e['y'], 0.8, gamma)
                g_pred = e['g_bar'] * nu_pred
                total += (np.log10(g_pred / e['g_obs']))**2
            return total

        res = minimize_scalar(chi2_sub, bounds=(0.05, 1.5), method='bounded')
        print(f"\n  {label}:")
        print(f"    N = {len(subset)}")
        print(f"    Best gamma = {res.x:.4f}")
        print(f"    chi2 = {res.fun:.4f}")
        print(f"    median y = {np.median([e['y'] for e in subset]):.2f}")
        print(f"    median nu_obs = {np.median([e['nu_obs'] for e in subset]):.2f}")

    # ================================================================
    # SECTION 5: Per-galaxy results
    # ================================================================
    print("\n" + "=" * 80)
    print("SECTION 5: PER-GALAXY COMPARISON")
    print("=" * 80)

    print(f"\n  {'Galaxy':<12s} {'type':<5s} {'y':<8s} {'nu_obs':<8s} "
          f"{'nu(0.36)':<9s} {'nu(0.54)':<9s} {'0.36/obs':<9s} {'0.54/obs':<9s} "
          f"{'prefers':<8s}")
    print("  " + "-" * 85)

    pref_036 = 0
    pref_054 = 0
    for e in sorted(ellipticals, key=lambda x: -x['M_star']):
        y = e['y']
        nu_o = e['nu_obs']
        nu_36 = nu_exp(y, 0.8, 0.36)
        nu_54 = nu_exp(y, 0.8, 0.54)

        err_36 = abs(np.log10(nu_36 / nu_o))
        err_54 = abs(np.log10(nu_54 / nu_o))
        pref = "g=0.36" if err_36 < err_54 else "g=0.54"
        if err_36 < err_54: pref_036 += 1
        else: pref_054 += 1

        sr = "*" if e['is_slow_rotator'] else " "
        print(f"  {e['name']:<12s} {e['morph']:<5s} {y:<8.2f} {nu_o:<8.3f} "
              f"{nu_36:<9.3f} {nu_54:<9.3f} {nu_36/nu_o:<9.3f} {nu_54/nu_o:<9.3f} "
              f"{pref:<8s} {sr}")

    print(f"\n  Prefer gamma=0.36 (disk-like): {pref_036}")
    print(f"  Prefer gamma=0.54 (sphere-like): {pref_054}")

    # ================================================================
    # SECTION 6: The discriminating power problem
    # ================================================================
    print("\n" + "=" * 80)
    print("SECTION 6: DISCRIMINATING POWER — CAN WE TELL gamma APART?")
    print("=" * 80)

    print(f"\n  Key issue: at high y (Newtonian regime), ALL models give nu ~ 1")
    print(f"  The difference between gamma=0.36 and gamma=0.54 is:")
    for y_test in [0.01, 0.05, 0.1, 0.5, 1.0, 5.0, 10.0, 50.0]:
        nu_36 = nu_exp(y_test, 0.8, 0.36)
        nu_54 = nu_exp(y_test, 0.8, 0.54)
        diff_pct = (nu_54 - nu_36) / nu_36 * 100
        print(f"    y={y_test:<6.2f}: nu(0.36)={nu_36:.4f}, nu(0.54)={nu_54:.4f}, "
              f"diff = {diff_pct:+.1f}%")

    print(f"\n  Ellipticals: median y = {np.median(ys):.1f}")
    y_med = np.median(ys)
    nu_36_med = nu_exp(y_med, 0.8, 0.36)
    nu_54_med = nu_exp(y_med, 0.8, 0.54)
    diff_med = (nu_54_med - nu_36_med) / nu_36_med * 100
    print(f"  At median y: nu(0.36)={nu_36_med:.4f}, nu(0.54)={nu_54_med:.4f}, "
          f"diff = {diff_med:+.1f}%")
    print(f"  => {'HARD to distinguish' if abs(diff_med) < 5 else 'CAN distinguish'} "
          f"at this y range")

    # Which ellipticals have lowest y (best discriminating power)?
    low_y = sorted(ellipticals, key=lambda e: e['y'])
    print(f"\n  Ellipticals with LOWEST y (best discriminating power):")
    print(f"  {'Galaxy':<12s} {'y':<8s} {'nu_obs':<8s} {'nu(0.36)':<9s} {'nu(0.54)':<9s} "
          f"{'diff%':<8s}")
    for e in low_y[:10]:
        nu36 = nu_exp(e['y'], 0.8, 0.36)
        nu54 = nu_exp(e['y'], 0.8, 0.54)
        diff = (nu54 - nu36) / nu36 * 100
        print(f"  {e['name']:<12s} {e['y']:<8.3f} {e['nu_obs']:<8.3f} "
              f"{nu36:<9.3f} {nu54:<9.3f} {diff:+8.1f}%")

    # ================================================================
    # SECTION 7: Combined disk + elliptical + cluster fit
    # ================================================================
    print("\n" + "=" * 80)
    print("SECTION 7: COMBINED FIT — DISK + ELLIPTICAL + CLUSTER")
    print("=" * 80)

    # Load SPARC disk galaxy summary stats from gs14
    # Use the median nu_obs at median y ~ 0.05 for disks
    print(f"\n  Summary across scales:")
    print(f"  {'System type':<25s} {'N':<5s} {'med y':<8s} {'med nu_obs':<10s} "
          f"{'best gamma':<11s} {'geometry'}")
    print("  " + "-" * 68)

    # Disk galaxies (from gs12): gamma = 0.40
    print(f"  {'Disk galaxies (SPARC)':<25s} {'171':<5s} {'0.05':<8s} {'~3-10':<10s} "
          f"{'0.40':<11s} {'disk'}")

    # Ellipticals
    print(f"  {'Ellipticals (this work)':<25s} {len(ellipticals):<5d} "
          f"{np.median(ys):<8.2f} {np.median(nus):<10.3f} "
          f"{gamma_best_fixed_a0:<11.3f} {'sphere'}")

    # Clusters (from gs16)
    print(f"  {'Galaxy clusters (gs13)':<25s} {'29':<5s} {'0.06':<8s} {'~7':<10s} "
          f"{'0.54':<11s} {'sphere'}")

    # ================================================================
    # SECTION 8: VERDICT
    # ================================================================
    print("\n" + "=" * 80)
    print("SECTION 8: VERDICT ON GEOMETRY HYPOTHESIS")
    print("=" * 80)

    # Check if ellipticals prefer gamma closer to 0.54 or 0.36
    gamma_ell = gamma_best_fixed_a0
    dist_to_disk = abs(gamma_ell - 0.36)
    dist_to_sphere = abs(gamma_ell - 0.54)

    print(f"""
  Geometry hypothesis prediction:
    Disks (galaxies):    gamma ~ 0.36  (d_eff = 2.29)
    Spheres (clusters):  gamma ~ 0.54  (d_eff = 1.92)
    Ellipticals (spheroidal galaxies): should prefer gamma ~ 0.54

  Observed:
    Ellipticals best gamma = {gamma_ell:.3f}

  Distance to predictions:
    |gamma_ell - gamma_disk|    = {dist_to_disk:.3f}
    |gamma_ell - gamma_sphere|  = {dist_to_sphere:.3f}

  Ellipticals are closer to: {'DISK (gamma=0.36)' if dist_to_disk < dist_to_sphere else 'SPHERE (gamma=0.54)'}
""")

    if dist_to_disk < dist_to_sphere:
        print(f"""  RESULT: Geometry hypothesis is {'WEAKENED' if dist_to_disk < 0.1 else 'AMBIGUOUS'}

  Possible interpretations:
  1. Ellipticals probe mainly HIGH-y regime (y ~ {np.median(ys):.0f})
     -> gamma poorly constrained (nu ~ 1 regardless of gamma)
     -> TEST IS INCONCLUSIVE at these accelerations
  2. Geometry hypothesis needs modification: it's not disk/sphere
     but something else (total mass? scale? environment?)
  3. Ellipticals DO follow disk-like gamma, and clusters are different
     for another reason (mass scale, baryonic composition)
""")
    else:
        print(f"""  RESULT: Geometry hypothesis SUPPORTED!

  Ellipticals (spheroidal) prefer gamma ~ {gamma_ell:.2f},
  closer to cluster gamma=0.54 than disk gamma=0.36.
  This confirms: spherical geometry -> deeper dimensional transition.
""")

    # Final note about discriminating power
    y_thresh = 0.5  # below this, gamma matters
    n_low_y = sum(1 for e in ellipticals if e['y'] < y_thresh)
    print(f"  Discriminating power:")
    print(f"    Galaxies with y < {y_thresh} (where gamma matters): {n_low_y}/{len(ellipticals)}")
    print(f"    {'GOOD test sample' if n_low_y > 10 else 'WEAK test — most ellipticals in Newtonian regime'}")

    print("\n" + "=" * 80)
    print("DONE")
    print("=" * 80)


if __name__ == '__main__':
    main()
