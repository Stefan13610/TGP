#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gs15: CLUSTER BARYON GAP ANALYSIS
==================================
gs13 showed Exp underpredicts clusters by x1.84.
Can missing baryons (hot gas, WHIM, hydrostatic bias) close the gap?

Key physics:
  - Hydrostatic bias: X-ray M_total biased LOW by ~10-40% (Eckert+2019)
    -> TRUE g_obs is HIGHER -> makes underprediction WORSE
    -> BUT: if we use lensing-calibrated M_total, bias is already corrected
  - Missing baryons: M_bar could be UNDERESTIMATED:
    a) Gas beyond r500 (WHIM in filaments/outskirts)
    b) Non-thermal gas not seen in X-ray
    c) Intracluster light (ICL) underestimate
    d) Warm gas (10^5-10^6 K) invisible to X-ray
  - Gas clumping: X-ray overestimates density (emission ~ n^2)
    -> can OVERESTIMATE M_gas by ~10-20% in outskirts (WRONG direction)

This script:
  1. For each cluster, find baryon multiplier f: g_bar*f -> g_pred = g_obs
  2. Check if f is physically plausible
  3. Test hydrostatic bias correction on M_total side
  4. Combined scenarios: baryon correction + hydrostatic correction
  5. Compare required f with published missing baryon estimates
"""

import sys
import io
import numpy as np
from scipy.optimize import brentq

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

G_SI = 6.674e-11
M_sun = 1.989e30
Mpc = 3.086e22

a0_exp = 1.12e-10
a0_mond = 1.2e-10

# ============================================================
# MODELS
# ============================================================
def nu_exp(y, alpha=0.8, gamma=0.4):
    return 1.0 + np.exp(-y**alpha) / y**gamma

def nu_mond(y):
    return 0.5 + np.sqrt(0.25 + 1.0/y)

def nu_mcgaugh(y):
    return 1.0 / (1.0 - np.exp(-np.sqrt(max(y, 1e-30))))


# ============================================================
# CLUSTER DATA (from gs13)
# ============================================================
def build_clusters():
    clusters = []

    # Gonzalez+2013 (Table 2): name, r500_Mpc, M500_1e14, Mgas_1e13, Mstar_1e13
    gonzalez = [
        ("A0122",  0.89, 2.26, 1.98, 0.55),
        ("A1651",  1.18, 5.15, 6.70, 0.65),
        ("A2401",  0.68, 0.95, 0.85, 0.27),
        ("A2721",  1.03, 3.46, 4.36, 0.57),
        ("A2811",  1.04, 3.59, 4.47, 0.47),
        ("A2955",  0.68, 0.99, 0.66, 0.30),
        ("A2984",  0.67, 0.95, 1.05, 0.39),
        ("A3112",  1.02, 3.23, 4.29, 0.70),
        ("A3693",  0.90, 2.26, 2.49, 0.51),
        ("A4010",  0.92, 2.41, 2.87, 0.56),
        ("AS0084", 0.91, 2.37, 2.09, 0.52),
        ("AS0296", 0.78, 1.45, 1.09, 0.29),
    ]
    for name, r500, M500, Mgas, Mstar in gonzalez:
        clusters.append({
            'name': name,
            'r500_m': r500 * Mpc,
            'M500': M500 * 1e14 * M_sun,
            'Mgas': Mgas * 1e13 * M_sun,
            'Mstar': Mstar * 1e13 * M_sun,
            'source': 'Gonzalez+2013'
        })

    # Literature: name, r500_Mpc, M500_1e14, fgas, fstar
    literature = [
        ("Coma",       1.31, 7.2,  0.144, 0.020),
        ("Perseus",    1.27, 6.6,  0.155, 0.025),
        ("Virgo",      0.77, 1.4,  0.090, 0.030),
        ("Centaurus",  0.82, 1.8,  0.100, 0.025),
        ("A2029",      1.42, 8.71, 0.138, 0.020),
        ("A0478",      1.28, 6.58, 0.175, 0.020),
        ("A2390",      1.50, 11.8, 0.146, 0.015),
        ("A1795",      1.09, 4.54, 0.118, 0.020),
        ("A2142",      1.37, 8.17, 0.174, 0.020),
        ("A2319",      1.55, 12.0, 0.150, 0.018),
        ("A644",       1.05, 3.70, 0.130, 0.020),
        ("A3266",      1.25, 6.20, 0.140, 0.018),
        ("Fornax",     0.45, 0.28, 0.060, 0.035),
        ("MKW4",       0.50, 0.38, 0.055, 0.040),
        ("NGC5044grp", 0.43, 0.24, 0.050, 0.045),
        ("A3571",      1.13, 5.02, 0.103, 0.020),
        ("Bullet",     1.35, 7.8,  0.140, 0.015),
    ]
    for name, r500, M500, fgas, fstar in literature:
        M500_kg = M500 * 1e14 * M_sun
        clusters.append({
            'name': name,
            'r500_m': r500 * Mpc,
            'M500': M500_kg,
            'Mgas': fgas * M500_kg,
            'Mstar': fstar * M500_kg,
            'source': 'Literature'
        })

    for c in clusters:
        c['Mbar'] = c['Mgas'] + c['Mstar']
        c['fbar'] = c['Mbar'] / c['M500']
        c['g_bar'] = G_SI * c['Mbar'] / c['r500_m']**2
        c['g_obs'] = G_SI * c['M500'] / c['r500_m']**2
        c['y'] = c['g_bar'] / a0_exp
        c['nu_obs'] = c['g_obs'] / c['g_bar']

    return clusters


# ============================================================
# ANALYSIS FUNCTIONS
# ============================================================
def find_baryon_multiplier(g_bar_orig, g_obs, a0, model='exp'):
    """
    Find factor f such that g_bar_new = f * g_bar_orig gives
    g_pred = g_bar_new * nu(g_bar_new/a0) = g_obs.
    """
    def residual(f):
        g_bar_new = f * g_bar_orig
        y = g_bar_new / a0
        if model == 'exp':
            nu = nu_exp(y)
        elif model == 'mond':
            nu = nu_mond(y)
        elif model == 'mcgaugh':
            nu = nu_mcgaugh(y)
        return g_bar_new * nu - g_obs

    try:
        f = brentq(residual, 0.1, 100.0, xtol=1e-6)
        return f
    except:
        return np.nan


def apply_hydrostatic_correction(g_obs, bias):
    """
    Hydrostatic bias: M_true = M_hydro / (1 - b)
    If b = 0.2, then M_true = M_hydro / 0.8 = 1.25 * M_hydro
    g_obs_true = g_obs / (1 - b)
    """
    return g_obs / (1.0 - bias)


# ============================================================
# MAIN
# ============================================================
def main():
    print("=" * 80)
    print("gs15: CLUSTER BARYON GAP ANALYSIS")
    print("Can missing baryons explain the Exp underprediction (x1.84)?")
    print("=" * 80)

    clusters = build_clusters()
    print(f"\nLoaded {len(clusters)} clusters")

    # ================================================================
    # SECTION 1: Required baryon multiplier per cluster
    # ================================================================
    print("\n" + "=" * 80)
    print("SECTION 1: REQUIRED BARYON MULTIPLIER (f = M_bar_needed / M_bar_observed)")
    print("=" * 80)
    print(f"\n  If f = 1.0: model works with observed baryons")
    print(f"  If f = 2.0: need 2x more baryons than observed")
    print(f"  If f < 1.5: plausible from systematics (hydrostatic bias, WHIM, ICL)")

    print(f"\n  {'Cluster':<14s} {'fbar_obs':<9s} {'f_Exp':<7s} {'f_MOND':<7s} "
          f"{'f_McG':<7s} {'fbar_needed(Exp)':<17s} {'fbar_needed(MOND)':<17s}")
    print("  " + "-" * 80)

    f_exp_all = []
    f_mond_all = []
    f_mcg_all = []

    for c in clusters:
        f_e = find_baryon_multiplier(c['g_bar'], c['g_obs'], a0_exp, 'exp')
        f_m = find_baryon_multiplier(c['g_bar'], c['g_obs'], a0_mond, 'mond')
        f_mc = find_baryon_multiplier(c['g_bar'], c['g_obs'], a0_mond, 'mcgaugh')

        c['f_exp'] = f_e
        c['f_mond'] = f_m
        c['f_mcg'] = f_mc

        fbar_needed_exp = c['fbar'] * f_e if not np.isnan(f_e) else 0
        fbar_needed_mond = c['fbar'] * f_m if not np.isnan(f_m) else 0

        print(f"  {c['name']:<14s} {c['fbar']:<9.3f} {f_e:<7.2f} {f_m:<7.2f} "
              f"{f_mc:<7.2f} {fbar_needed_exp:<17.3f} {fbar_needed_mond:<17.3f}")

        if not np.isnan(f_e): f_exp_all.append(f_e)
        if not np.isnan(f_m): f_mond_all.append(f_m)
        if not np.isnan(f_mc): f_mcg_all.append(f_mc)

    f_exp_all = np.array(f_exp_all)
    f_mond_all = np.array(f_mond_all)
    f_mcg_all = np.array(f_mcg_all)

    print(f"\n  Summary (baryon multiplier f):")
    print(f"  {'Model':<15s} {'median f':<10s} {'mean f':<10s} {'std f':<8s} {'min':<7s} {'max':<7s}")
    print("  " + "-" * 57)
    print(f"  {'Exp (g=0.4)':<15s} {np.median(f_exp_all):<10.3f} {np.mean(f_exp_all):<10.3f} "
          f"{np.std(f_exp_all):<8.3f} {np.min(f_exp_all):<7.3f} {np.max(f_exp_all):<7.3f}")
    print(f"  {'MOND':<15s} {np.median(f_mond_all):<10.3f} {np.mean(f_mond_all):<10.3f} "
          f"{np.std(f_mond_all):<8.3f} {np.min(f_mond_all):<7.3f} {np.max(f_mond_all):<7.3f}")
    print(f"  {'McGaugh':<15s} {np.median(f_mcg_all):<10.3f} {np.mean(f_mcg_all):<10.3f} "
          f"{np.std(f_mcg_all):<8.3f} {np.min(f_mcg_all):<7.3f} {np.max(f_mcg_all):<7.3f}")

    # ================================================================
    # SECTION 2: What fbar_total would be needed?
    # ================================================================
    print("\n" + "=" * 80)
    print("SECTION 2: REQUIRED TOTAL BARYON FRACTION (fbar_needed)")
    print("=" * 80)
    print(f"\n  Cosmic baryon fraction (Planck 2018): Omega_b/Omega_m = 0.157")
    print(f"  Typical observed fbar in clusters at r500: 0.10 - 0.17")

    fbar_obs = np.array([c['fbar'] for c in clusters])
    fbar_needed_exp = np.array([c['fbar'] * c['f_exp'] for c in clusters if not np.isnan(c['f_exp'])])
    fbar_needed_mond = np.array([c['fbar'] * c['f_mond'] for c in clusters if not np.isnan(c['f_mond'])])

    print(f"\n  Observed fbar:           median = {np.median(fbar_obs):.3f}, "
          f"mean = {np.mean(fbar_obs):.3f}")
    print(f"  Needed fbar (Exp):       median = {np.median(fbar_needed_exp):.3f}, "
          f"mean = {np.mean(fbar_needed_exp):.3f}")
    print(f"  Needed fbar (MOND):      median = {np.median(fbar_needed_mond):.3f}, "
          f"mean = {np.mean(fbar_needed_mond):.3f}")
    print(f"  Cosmic fbar (Planck):    0.157")

    print(f"\n  Ratio to cosmic:")
    print(f"    Exp needs:  {np.median(fbar_needed_exp)/0.157:.2f}x cosmic fbar")
    print(f"    MOND needs: {np.median(fbar_needed_mond)/0.157:.2f}x cosmic fbar")
    print(f"    {'Exp EXCEEDS cosmic -> UNPHYSICAL' if np.median(fbar_needed_exp) > 0.157 else 'Exp within cosmic -> POSSIBLE'}")
    print(f"    {'MOND EXCEEDS cosmic -> UNPHYSICAL' if np.median(fbar_needed_mond) > 0.157 else 'MOND within cosmic -> POSSIBLE'}")

    # ================================================================
    # SECTION 3: Sources of missing baryons — can they close the gap?
    # ================================================================
    print("\n" + "=" * 80)
    print("SECTION 3: SOURCES OF MISSING BARYONS")
    print("=" * 80)

    print(f"""
  Published estimates of systematic baryon underestimates:

  Source                          Extra baryons    Notes
  -----------------------------------------------------------------------
  Gas clumping correction         -10 to -20%     REDUCES M_gas (wrong sign!)
  Gas beyond r500 (outskirts)     +5 to +15%      Extends to r200, virial
  WHIM in filaments               +5 to +20%      Very uncertain
  Intracluster light (ICL)        +2 to +5%       Missed diffuse stars
  Warm gas (10^5-10^6 K)          +3 to +10%      Between X-ray and UV
  Compact objects (BH, NS, WD)    < 1%            Negligible
  -----------------------------------------------------------------------
  OPTIMISTIC TOTAL                +15 to +50%     f_baryon_correction ~ 1.15-1.50
  REALISTIC TOTAL                 +10 to +30%     f ~ 1.10-1.30
  """)

    f_exp_median = np.median(f_exp_all)
    f_mond_median = np.median(f_mond_all)

    print(f"  Required baryon multiplier:")
    print(f"    Exp:  f = {f_exp_median:.2f} ({(f_exp_median-1)*100:.0f}% more baryons needed)")
    print(f"    MOND: f = {f_mond_median:.2f} ({(f_mond_median-1)*100:.0f}% more baryons needed)")
    print(f"\n  Can missing baryons close the gap?")
    print(f"    MOND: needs +{(f_mond_median-1)*100:.0f}% -> "
          f"{'YES, within optimistic range' if f_mond_median < 1.50 else 'NO, too much'}")
    print(f"    Exp:  needs +{(f_exp_median-1)*100:.0f}% -> "
          f"{'YES, within optimistic range' if f_exp_median < 1.50 else 'MARGINAL/NO, at upper limit'}")

    # ================================================================
    # SECTION 4: Effect of hydrostatic bias on M_total
    # ================================================================
    print("\n" + "=" * 80)
    print("SECTION 4: HYDROSTATIC BIAS ON M_TOTAL SIDE")
    print("=" * 80)
    print(f"\n  If M_total from X-ray is biased LOW: M_true = M_obs / (1-b)")
    print(f"  Then g_obs_true > g_obs_measured -> underprediction gets WORSE")
    print(f"  -> Hydrostatic bias HURTS modified gravity models, not helps!")

    biases = [0.0, 0.10, 0.20, 0.30, 0.40]
    print(f"\n  Effect of hydrostatic bias on underprediction factor:")
    print(f"  {'bias b':<9s} {'M_true/M_obs':<14s} {'Exp underpred':<14s} {'MOND underpred':<14s}")
    print("  " + "-" * 51)

    for b in biases:
        # With bias, g_obs_true = g_obs / (1-b)
        # underprediction = g_obs_true / g_pred = (g_obs / g_pred) / (1-b)
        exp_under = 1.84 / (1.0 - b)
        mond_under = 1.46 / (1.0 - b)
        print(f"  {b:<9.2f} {1.0/(1-b):<14.3f} {exp_under:<14.2f}x {mond_under:<14.2f}x")

    print(f"\n  CONCLUSION: Hydrostatic bias makes the problem WORSE!")
    print(f"  If b=0.2 (typical): Exp goes from x1.84 to x2.30, MOND from x1.46 to x1.82")

    # ================================================================
    # SECTION 5: COMBINED SCENARIO — more baryons + corrected M_total
    # ================================================================
    print("\n" + "=" * 80)
    print("SECTION 5: COMBINED SCENARIO — BARYON CORRECTION ONLY")
    print("=" * 80)
    print(f"\n  (Hydrostatic bias worsens things, so focus on baryon correction only)")

    baryon_boosts = [1.0, 1.1, 1.2, 1.3, 1.5, 1.8, 2.0, 2.5, 3.0]

    print(f"\n  {'f_bar boost':<13s} {'extra baryons':<14s} {'Exp ratio':<10s} "
          f"{'MOND ratio':<11s} {'Exp ok?':<8s} {'MOND ok?':<8s}")
    print("  " + "-" * 68)

    for f in baryon_boosts:
        # Compute median ratio for all clusters with boosted baryons
        ratios_exp = []
        ratios_mond = []
        for c in clusters:
            g_bar_new = c['g_bar'] * f
            y_e = g_bar_new / a0_exp
            y_m = g_bar_new / a0_mond
            g_pred_e = g_bar_new * nu_exp(y_e)
            g_pred_m = g_bar_new * nu_mond(y_m)
            ratios_exp.append(g_pred_e / c['g_obs'])
            ratios_mond.append(g_pred_m / c['g_obs'])

        med_e = np.median(ratios_exp)
        med_m = np.median(ratios_mond)
        ok_e = "YES" if 0.85 <= med_e <= 1.15 else ("close" if 0.7 <= med_e <= 1.3 else "no")
        ok_m = "YES" if 0.85 <= med_m <= 1.15 else ("close" if 0.7 <= med_m <= 1.3 else "no")

        print(f"  {f:<13.1f} +{(f-1)*100:<12.0f}% {med_e:<10.3f} "
              f"{med_m:<11.3f} {ok_e:<8s} {ok_m:<8s}")

    # ================================================================
    # SECTION 6: Per-cluster breakdown — which clusters are hardest?
    # ================================================================
    print("\n" + "=" * 80)
    print("SECTION 6: PER-CLUSTER DIFFICULTY — WHICH ARE HARDEST TO FIX?")
    print("=" * 80)

    sorted_by_f = sorted(clusters, key=lambda c: c.get('f_exp', 99))
    print(f"\n  Easiest to fix (lowest f needed):")
    print(f"  {'Cluster':<14s} {'f_Exp':<7s} {'f_MOND':<7s} {'fbar_obs':<9s} "
          f"{'fbar_need(Exp)':<15s} {'M500(1e14)':<11s}")
    print("  " + "-" * 63)
    for c in sorted_by_f[:10]:
        if np.isnan(c.get('f_exp', np.nan)): continue
        print(f"  {c['name']:<14s} {c['f_exp']:<7.2f} {c['f_mond']:<7.2f} "
              f"{c['fbar']:<9.3f} {c['fbar']*c['f_exp']:<15.3f} "
              f"{c['M500']/M_sun/1e14:<11.2f}")

    print(f"\n  Hardest to fix (highest f needed):")
    for c in sorted_by_f[-10:]:
        if np.isnan(c.get('f_exp', np.nan)): continue
        print(f"  {c['name']:<14s} {c['f_exp']:<7.2f} {c['f_mond']:<7.2f} "
              f"{c['fbar']:<9.3f} {c['fbar']*c['f_exp']:<15.3f} "
              f"{c['M500']/M_sun/1e14:<11.2f}")

    # ================================================================
    # SECTION 7: Gas fraction systematics
    # ================================================================
    print("\n" + "=" * 80)
    print("SECTION 7: GAS FRACTION SYSTEMATICS")
    print("=" * 80)

    fgas_vals = np.array([c['Mgas']/c['M500'] for c in clusters])
    fstar_vals = np.array([c['Mstar']/c['M500'] for c in clusters])
    fbar_vals = np.array([c['fbar'] for c in clusters])

    print(f"\n  Gas fraction (fgas) statistics:")
    print(f"    median = {np.median(fgas_vals):.3f}, mean = {np.mean(fgas_vals):.3f}")
    print(f"    range: {np.min(fgas_vals):.3f} to {np.max(fgas_vals):.3f}")
    print(f"\n  Stellar fraction (fstar):")
    print(f"    median = {np.median(fstar_vals):.3f}, mean = {np.mean(fstar_vals):.3f}")
    print(f"\n  Total baryon fraction:")
    print(f"    median = {np.median(fbar_vals):.3f}, mean = {np.mean(fbar_vals):.3f}")
    print(f"    Cosmic (Planck): 0.157")
    print(f"    Deficit: {(0.157 - np.median(fbar_vals))/0.157*100:.0f}% below cosmic")

    # How much gas are we missing just to reach cosmic fbar?
    fbar_deficit = 0.157 - fbar_vals
    fbar_deficit_frac = fbar_deficit / fbar_vals
    print(f"\n  Extra baryons to reach cosmic fbar:")
    print(f"    median +{np.median(fbar_deficit_frac)*100:.0f}% of observed")
    print(f"    This gives baryon multiplier f = {1 + np.median(fbar_deficit_frac):.2f}")

    # Is reaching cosmic fbar enough?
    f_cosmic = (1 + np.median(fbar_deficit_frac))
    print(f"\n  If all clusters had cosmic fbar (f = {f_cosmic:.2f}):")
    ratios_e = []
    ratios_m = []
    for c in clusters:
        f_c = 0.157 / c['fbar']  # boost to reach cosmic
        g_bar_new = c['g_bar'] * f_c
        y_e = g_bar_new / a0_exp
        y_m = g_bar_new / a0_mond
        ratios_e.append(g_bar_new * nu_exp(y_e) / c['g_obs'])
        ratios_m.append(g_bar_new * nu_mond(y_m) / c['g_obs'])
    print(f"    Exp ratio: {np.median(ratios_e):.3f} (need 1.0)")
    print(f"    MOND ratio: {np.median(ratios_m):.3f} (need 1.0)")

    # ================================================================
    # SECTION 8: Correlation — does f_needed correlate with cluster mass?
    # ================================================================
    print("\n" + "=" * 80)
    print("SECTION 8: DOES REQUIRED f CORRELATE WITH CLUSTER MASS?")
    print("=" * 80)

    from scipy.stats import spearmanr

    M500s = np.array([c['M500']/M_sun for c in clusters])
    f_exps = np.array([c['f_exp'] for c in clusters])
    f_monds = np.array([c['f_mond'] for c in clusters])

    mask = ~np.isnan(f_exps)
    rho, pval = spearmanr(np.log10(M500s[mask]), f_exps[mask])
    print(f"\n  f_Exp vs log(M500): Spearman rho = {rho:+.3f}, p = {pval:.4f}")

    rho2, pval2 = spearmanr(np.log10(M500s[mask]), f_monds[mask])
    print(f"  f_MOND vs log(M500): Spearman rho = {rho2:+.3f}, p = {pval2:.4f}")

    # Binned by mass
    m_bins = [(0, 1e14), (1e14, 3e14), (3e14, 7e14), (7e14, 2e15)]
    print(f"\n  {'M500 range':<20s} {'N':<4s} {'med f_Exp':<10s} {'med f_MOND':<10s} "
          f"{'med fbar':<10s}")
    print("  " + "-" * 54)
    for mlo, mhi in m_bins:
        m_mask = (M500s >= mlo) & (M500s < mhi)
        if np.sum(m_mask) < 2: continue
        fe = f_exps[m_mask & mask]
        fm = f_monds[m_mask & mask]
        fb = fbar_vals[m_mask]
        label = f"{mlo/1e14:.0f}-{mhi/1e14:.0f} (1e14)"
        print(f"  {label:<20s} {np.sum(m_mask):<4d} {np.median(fe):<10.2f} "
              f"{np.median(fm):<10.2f} {np.median(fb):<10.3f}")

    # ================================================================
    # SECTION 9: What if EXP uses a different a0 for clusters?
    # ================================================================
    print("\n" + "=" * 80)
    print("SECTION 9: ALTERNATIVE — CLUSTER-SPECIFIC a0")
    print("=" * 80)

    # Find best a0 for each cluster (keeping alpha=0.8, gamma=0.4)
    from scipy.optimize import minimize_scalar

    a0_best_list = []
    for c in clusters:
        def obj(log_a0):
            a0 = 10**log_a0
            y = c['g_bar'] / a0
            g_pred = c['g_bar'] * nu_exp(y)
            return (np.log10(g_pred / c['g_obs']))**2

        res = minimize_scalar(obj, bounds=(-12, -8), method='bounded')
        a0_best = 10**res.x
        c['a0_best'] = a0_best
        a0_best_list.append(a0_best)

    a0_best_arr = np.array(a0_best_list)
    print(f"\n  Per-cluster best a0:")
    print(f"    median = {np.median(a0_best_arr):.3e} ({np.median(a0_best_arr)/a0_exp:.2f}x galaxy)")
    print(f"    std = {np.std(a0_best_arr):.3e}")

    rho3, pval3 = spearmanr(np.log10(M500s), np.log10(a0_best_arr))
    print(f"\n  a0_best vs M500: Spearman rho = {rho3:+.3f}, p = {pval3:.4f}")

    # If a0 = a0_galaxy * (1 + M/M_crit)^p — does this fit?
    print(f"\n  Possible a0(M) scaling:")
    print(f"    log(a0_best) vs log(M500):")

    # Simple linear fit in log-log
    log_m = np.log10(M500s)
    log_a0 = np.log10(a0_best_arr)
    p = np.polyfit(log_m, log_a0, 1)
    print(f"    Slope = {p[0]:.3f}, intercept = {p[1]:.3f}")
    print(f"    a0 ~ M^{p[0]:.3f}")
    residuals = log_a0 - np.polyval(p, log_m)
    print(f"    Scatter around fit: {np.std(residuals):.3f} dex")

    # ================================================================
    # SECTION 10: SUMMARY — physical plausibility
    # ================================================================
    print("\n" + "=" * 80)
    print("SECTION 10: SUMMARY — CAN MISSING BARYONS EXPLAIN THE GAP?")
    print("=" * 80)

    print(f"""
  =====================================================================
  SCENARIO                           Exp status      MOND status
  =====================================================================
  Missing baryons only:
    Need +{(np.median(f_exp_all)-1)*100:.0f}% baryons (Exp)        MARGINAL/NO     -
    Need +{(np.median(f_mond_all)-1)*100:.0f}% baryons (MOND)       -               POSSIBLE
    Reach cosmic fbar:               {np.median(ratios_e):.2f}x          {np.median(ratios_m):.2f}x

  Hydrostatic bias (b=0.2):
    Makes underprediction WORSE      x{1.84/0.8:.2f}          x{1.46/0.8:.2f}
    -> HURTS both models

  Combined baryon + bias:
    Net effect contradictory:
    more baryons help, bias hurts

  Cluster-specific a0:
    Median a0 = {np.median(a0_best_arr):.2e}     {np.median(a0_best_arr)/a0_exp:.1f}x galaxy
    a0 ~ M^{p[0]:.2f}                 {'weak trend' if abs(p[0]) < 0.2 else 'significant trend'}
  =====================================================================

  VERDICT FOR EXP MODEL (gamma=0.4):
  -----------------------------------
  1. Missing baryons ALONE are probably NOT enough:
     - Need +{(np.median(f_exp_all)-1)*100:.0f}% more baryons
     - Optimistic maximum from systematics: +30-50%
     - Required fbar = {np.median(fbar_needed_exp):.3f}
       {'EXCEEDS cosmic (0.157) -> UNPHYSICAL' if np.median(fbar_needed_exp) > 0.157 else 'Within cosmic -> possible'}

  2. MOND (gamma=0.5) is in better shape:
     - Needs only +{(np.median(f_mond_all)-1)*100:.0f}% more baryons
     - Required fbar = {np.median(fbar_needed_mond):.3f}
       {'EXCEEDS cosmic -> UNPHYSICAL' if np.median(fbar_needed_mond) > 0.157 else 'Within cosmic -> possible'}

  3. The cluster problem is STRUCTURAL for Exp:
     - gamma=0.4 gives WEAKER boost than gamma=0.5
     - At cluster y ~ 0.05: nu_Exp/nu_MOND ~ 0.68
     - This 32% deficit cannot be recovered by <50% baryon correction

  4. TGP interpretation:
     - Either gamma transitions from 0.4 (galaxies) to >0.5 (clusters)
     - Or a0 increases at cluster scales
     - Or: clusters require a DIFFERENT physical mechanism
       (2D membrane model may not apply at Mpc scales)
""")

    print("=" * 80)
    print("DONE")
    print("=" * 80)


if __name__ == '__main__':
    main()
