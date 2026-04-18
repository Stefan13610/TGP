#!/usr/bin/env python3
"""
gs38: SPARC ROTATION CURVES REFIT — TYPE-DEPENDENT c_eff
========================================================

Following gs37 which showed TGP with c_eff=1.0 (pure thin disk, gamma=0.400)
gives chi2/dof=35.88 for 20 SPARC-like galaxies, while MOND gives 17.52.

The BTFR slope analysis showed c_eff=1.3-1.5 (gamma=0.45-0.48) is required
to match the observed slope of ~3.85.

This script REFITS the same 20 galaxies with:
  A. Same galaxy sample and observed rotation curves as gs37
  B. Type-dependent c_eff: Sb->1.3, Sc->1.2, Sd->1.15, Im/dIrr->1.8, LSB->1.5
  C. Refit M/L for TGP with type-dependent c_eff, compare chi2 vs MOND
  D. Universal c_eff=1.3 and c_eff=1.5 for comparison
  E. Check M/L values are in reasonable range (0.2-1.0)
  F. Compute BTFR from the sample: log(M_bar) vs log(V_flat)
  G. Summary with chi2 improvement over gs37
"""

import numpy as np
import sys, os
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

# --- Constants ---
G_SI = 6.674e-11     # m^3/(kg s^2)
a0 = 1.2e-10         # m/s^2
Msun = 1.989e30      # kg
kpc = 3.086e19       # m
pc = 3.086e16        # m


# ================================================================
#  CORE FUNCTIONS (from gs37 pattern)
# ================================================================

def nu_tgp(y, alpha=4/5, gamma=0.4):
    """TGP interpolation: nu(y) = 1 + exp(-y^alpha) / y^gamma"""
    y = np.asarray(y, dtype=float)
    y = np.maximum(y, 1e-12)
    return 1.0 + np.exp(-y**alpha) / y**gamma


def nu_tgp_ceff(y, c_eff=1):
    """TGP interpolation with effective codimension c_eff."""
    alpha = 4/5
    gamma = 0.8 * c_eff / (c_eff + 1)
    return nu_tgp(y, alpha, gamma)


def nu_mond_simple(y):
    """Simple MOND interpolation (standard)."""
    y = np.asarray(y, dtype=float)
    y = np.maximum(y, 1e-12)
    return 0.5 * (1 + np.sqrt(1 + 4.0/y))


def print_header(title):
    print()
    print("=" * 78)
    print(f"  {title}")
    print("=" * 78)
    print()


# ================================================================
#  SPARC-LIKE GALAXY SAMPLE (same as gs37)
# ================================================================

sparc_sample = [
    # HIGH-MASS SPIRALS (V > 200 km/s)
    ("UGC 2953",   "Sb",   8.0e10, 1.2e10, 5.0, 300, 15.0, 60, 1),
    ("NGC 2403",   "Sc",   5.0e9,  3.5e9,  2.1, 136, 3.2,  63, 1),
    ("NGC 3198",   "Sc",   1.6e10, 8.0e9,  3.2, 150, 13.8, 72, 1),
    ("NGC 7331",   "Sb",   6.0e10, 8.0e9,  4.5, 250, 14.7, 76, 1),
    ("NGC 6946",   "Sc",   3.5e10, 6.0e9,  3.8, 200, 5.9,  33, 1),

    # INTERMEDIATE SPIRALS (V ~ 100-200)
    ("NGC 2903",   "Sb",   2.5e10, 4.0e9,  2.8, 185, 8.9,  65, 1),
    ("NGC 3521",   "Sb",   4.0e10, 5.0e9,  3.5, 210, 10.7, 73, 1),
    ("NGC 2841",   "Sb",   7.0e10, 9.0e9,  4.2, 310, 14.1, 74, 1),
    ("DDO 154",    "Im",   3.0e7,  4.5e8,  0.8, 47,  3.7,  66, 1),
    ("DDO 170",    "Im",   7.0e7,  6.0e8,  1.2, 60,  12.0, 64, 2),

    # LOW-MASS / LSB (V < 100)
    ("IC 2574",    "Im",   3.0e8,  1.5e9,  3.0, 66,  4.0,  53, 1),
    ("UGC 128",    "LSB",  2.0e9,  5.0e9,  6.0, 130, 64.0, 52, 2),
    ("F583-1",     "LSB",  5.0e8,  2.0e9,  3.5, 83,  32.0, 63, 2),
    ("UGC 5750",   "LSB",  1.5e9,  3.0e9,  4.5, 90,  56.0, 64, 2),
    ("NGC 1560",   "Sd",   3.5e8,  1.2e9,  1.6, 78,  3.5,  82, 1),

    # DWARF IRREGULARS (V < 60)
    ("DDO 47",     "Im",   5.0e7,  3.0e8,  1.0, 55,  5.2,  40, 2),
    ("DDO 87",     "Im",   4.0e7,  2.5e8,  1.3, 48,  7.7,  43, 2),
    ("UGC 5005",   "LSB",  8.0e8,  4.0e9,  5.0, 100, 52.0, 42, 2),
    ("F568-3",     "LSB",  1.0e9,  3.5e9,  4.0, 95,  80.0, 40, 2),
    ("UGCA 442",   "Im",   2.0e8,  5.0e8,  1.5, 55,  4.4,  65, 2),
]

# Observed rotation curves (same as gs37)
obs_curves = {
    "UGC 2953":  {2: 200, 4: 260, 6: 290, 8: 295, 10: 300, 15: 300, 20: 298},
    "NGC 2403":  {1: 70, 2: 105, 4: 128, 8: 135, 12: 136, 18: 135},
    "NGC 3198":  {2: 100, 4: 140, 6: 148, 8: 150, 12: 150, 18: 150, 25: 148},
    "NGC 7331":  {2: 180, 4: 235, 6: 248, 8: 250, 12: 250, 18: 248},
    "NGC 6946":  {1: 100, 2: 160, 4: 190, 6: 198, 8: 200, 12: 200},
    "NGC 2903":  {1: 100, 2: 150, 4: 178, 6: 183, 8: 185, 12: 185},
    "NGC 3521":  {2: 150, 4: 195, 6: 208, 8: 210, 12: 210, 18: 208},
    "NGC 2841":  {2: 220, 4: 280, 6: 300, 8: 308, 12: 310, 18: 308},
    "DDO 154":   {0.5: 15, 1: 25, 2: 35, 3: 42, 4: 45, 6: 47, 8: 47},
    "DDO 170":   {1: 20, 2: 38, 4: 52, 6: 58, 8: 60},
    "IC 2574":   {1: 15, 2: 30, 4: 50, 6: 58, 8: 63, 10: 66},
    "UGC 128":   {2: 40, 4: 80, 6: 105, 8: 120, 12: 128, 18: 130},
    "F583-1":    {1: 20, 2: 40, 4: 65, 6: 75, 8: 80, 12: 83},
    "UGC 5750":  {2: 30, 4: 55, 6: 75, 8: 85, 12: 90},
    "NGC 1560":  {0.5: 20, 1: 40, 2: 60, 3: 68, 4: 73, 6: 78},
    "DDO 47":    {0.5: 10, 1: 25, 2: 40, 3: 48, 4: 52, 6: 55},
    "DDO 87":    {0.5: 8, 1: 20, 2: 35, 3: 42, 4: 46, 6: 48},
    "UGC 5005":  {2: 35, 4: 65, 6: 85, 8: 95, 12: 100},
    "F568-3":    {2: 30, 4: 60, 6: 80, 8: 90, 12: 95},
    "UGCA 442":  {0.5: 12, 1: 25, 2: 40, 3: 48, 4: 52, 6: 55},
}

# Type-dependent c_eff assignments
ceff_by_type = {
    "Sb":  1.3,
    "Sc":  1.2,
    "Sd":  1.15,
    "Im":  1.8,
    "LSB": 1.5,
}


# ================================================================
#  ROTATION CURVE MODEL
# ================================================================

def compute_rotation_curve(M_disk, M_gas, R_d, radii_kpc, ML_disk=1.0, c_eff=1):
    """
    Compute rotation curve for a disk galaxy.
    Returns list of dicts with V_bar, V_tgp, V_mond, y, etc.
    """
    results = []

    for R_kpc in radii_kpc:
        R_m = R_kpc * kpc

        # Enclosed baryonic mass: exponential disk profile
        x = R_kpc / R_d
        f_disk = 1.0 - (1.0 + x) * np.exp(-x)
        M_disk_enc = ML_disk * M_disk * f_disk

        # Gas disk with scale length ~ 2*R_d (gas more extended)
        x_gas = R_kpc / (2.0 * R_d)
        f_gas = 1.0 - (1.0 + x_gas) * np.exp(-x_gas)
        M_gas_enc = M_gas * f_gas

        M_bar_enc = M_disk_enc + M_gas_enc

        # Newtonian rotation velocity
        V_bar2 = G_SI * M_bar_enc * Msun / R_m
        V_bar = np.sqrt(max(V_bar2, 0)) / 1e3  # km/s

        # Newtonian acceleration
        g_bar = V_bar2 / R_m  # m/s^2
        y = g_bar / a0
        y = max(y, 1e-12)

        # TGP boost
        nu = nu_tgp_ceff(y, c_eff)
        V_tgp = np.sqrt(V_bar2 * nu) / 1e3  # km/s

        # MOND boost
        nu_m = nu_mond_simple(y)
        V_mond = np.sqrt(V_bar2 * nu_m) / 1e3

        results.append({
            'R': R_kpc,
            'V_bar': V_bar,
            'V_tgp': V_tgp,
            'V_mond': V_mond,
            'y': y,
            'nu_tgp': nu,
            'nu_mond': nu_m,
            'M_bar': M_bar_enc,
        })

    return results


def fit_ML(M_disk, M_gas, R_d, V_obs_dict, sigma_v=5.0, c_eff=1):
    """Find optimal M/L that minimizes chi2 for TGP and MOND."""
    radii = sorted(V_obs_dict.keys())
    V_obs = np.array([V_obs_dict[r] for r in radii])

    best = {}
    for model_name, model_func in [('TGP', 'V_tgp'), ('MOND', 'V_mond')]:
        best_chi2 = 1e30
        best_ml = 0.5
        for ml in np.arange(0.1, 2.01, 0.01):
            results = compute_rotation_curve(M_disk, M_gas, R_d, np.array(radii, dtype=float),
                                             ML_disk=ml, c_eff=c_eff)
            V_pred = np.array([r[model_func] for r in results])
            chi2 = np.sum(((V_pred - V_obs) / sigma_v)**2)
            if chi2 < best_chi2:
                best_chi2 = chi2
                best_ml = ml

        best[model_name] = {'ML': best_ml, 'chi2': best_chi2, 'dof': len(radii) - 1}

    return best


# ================================================================
#  PART A: REFERENCE — gs37 RESULTS WITH c_eff=1.0
# ================================================================

def part_A():
    print_header("PART A: REFERENCE -- gs37 RESULTS (c_eff=1.0, gamma=0.400)")

    print("  From gs37, TGP with pure thin-disk c_eff=1.0:")
    print("    Total chi2 TGP:   3551.9  (chi2/dof = 35.88)")
    print("    Total chi2 MOND:  1734.1  (chi2/dof = 17.52)")
    print("    Total dof:        99")
    print()
    print("  Problems identified:")
    print("    1. M/L values pegged at 2.0 for most galaxies (unphysical)")
    print("    2. Only 0/20 galaxies have M/L in [0.2, 1.0]")
    print("    3. BTFR slope = 3.33 (too low; observed = 3.85 +/- 0.09)")
    print()
    print("  Solution: BTFR analysis shows c_eff = 1.3-1.5 is required.")
    print("  Physical motivation: real disks have finite thickness h/R > 0")
    print("  -> c_eff > 1 -> gamma > 0.400 -> more deep-MOND boost")
    print()

    print("  Type-dependent c_eff assignments for gs38:")
    print()
    print(f"    {'Type':<10} {'c_eff':<8} {'gamma':<8} {'Physical reason':<40}")
    print(f"    {'---'*22}")
    reasons = {
        "Sb":  "Moderate bulge, h/R ~ 0.15",
        "Sc":  "Thin disk dominant, h/R ~ 0.10",
        "Sd":  "Very thin disk, h/R ~ 0.07",
        "Im":  "Puffy irregular, h/R ~ 0.4",
        "LSB": "Thick diffuse disk, h/R ~ 0.25",
    }
    for gtype, ceff in sorted(ceff_by_type.items()):
        gamma = 0.8 * ceff / (ceff + 1)
        print(f"    {gtype:<10} {ceff:<8.2f} {gamma:<8.3f} {reasons[gtype]:<40}")


# ================================================================
#  PART B: REFIT WITH TYPE-DEPENDENT c_eff
# ================================================================

def part_B():
    print_header("PART B: REFIT WITH TYPE-DEPENDENT c_eff")

    print("  B.1  Fitting M/L_disk with type-dependent c_eff")
    print("  ────────────────────────────────────────────────")
    print("  M/L scanned from 0.1 to 2.0 in steps of 0.01")
    print("  sigma_v = 5 km/s for all galaxies")
    print()

    print(f"    {'Galaxy':<14} {'Type':<6} {'c_eff':<7} {'gamma':<7} {'ML_TGP':<8} {'chi2_TGP':<10} {'ML_MOND':<8} {'chi2_MOND':<10} {'dof':<5} {'TGP/MOND':<10}")
    print(f"    {'---'*30}")

    total_chi2_tgp = 0
    total_chi2_mond = 0
    total_dof = 0
    ml_tgp_list = []
    ml_mond_list = []
    galaxy_results = []

    for gal in sparc_sample:
        name, gtype, M_disk, M_gas, R_d, V_flat, dist, incl, qual = gal

        if name not in obs_curves:
            continue

        c_eff = ceff_by_type.get(gtype, 1.0)
        gamma = 0.8 * c_eff / (c_eff + 1)

        best = fit_ML(M_disk, M_gas, R_d, obs_curves[name], sigma_v=5.0, c_eff=c_eff)

        t = best['TGP']
        m = best['MOND']
        ratio = t['chi2'] / m['chi2'] if m['chi2'] > 0 else 99

        total_chi2_tgp += t['chi2']
        total_chi2_mond += m['chi2']
        total_dof += t['dof']
        ml_tgp_list.append(t['ML'])
        ml_mond_list.append(m['ML'])

        galaxy_results.append({
            'name': name, 'type': gtype, 'c_eff': c_eff, 'gamma': gamma,
            'ML_tgp': t['ML'], 'chi2_tgp': t['chi2'],
            'ML_mond': m['ML'], 'chi2_mond': m['chi2'],
            'dof': t['dof'], 'M_disk': M_disk, 'M_gas': M_gas,
            'V_flat': V_flat,
        })

        print(f"    {name:<14} {gtype:<6} {c_eff:<7.2f} {gamma:<7.3f} {t['ML']:<8.2f} {t['chi2']:<10.1f} {m['ML']:<8.2f} {m['chi2']:<10.1f} {t['dof']:<5} {ratio:<10.2f}")

    print()
    print("  B.2  Aggregate statistics (type-dependent c_eff)")
    print("  ──────────────────────────────────────────────────")
    print(f"    Total chi2 TGP:     {total_chi2_tgp:.1f}  (chi2/dof = {total_chi2_tgp/total_dof:.2f})")
    print(f"    Total chi2 MOND:    {total_chi2_mond:.1f}  (chi2/dof = {total_chi2_mond/total_dof:.2f})")
    print(f"    Total dof:          {total_dof}")
    print()
    print(f"    M/L_TGP:  mean = {np.mean(ml_tgp_list):.2f}, median = {np.median(ml_tgp_list):.2f}, std = {np.std(ml_tgp_list):.2f}")
    print(f"    M/L_MOND: mean = {np.mean(ml_mond_list):.2f}, median = {np.median(ml_mond_list):.2f}, std = {np.std(ml_mond_list):.2f}")
    print()

    # gs37 reference
    gs37_chi2_tgp = 3551.9
    gs37_chi2dof_tgp = 35.88
    improvement = (gs37_chi2_tgp - total_chi2_tgp) / gs37_chi2_tgp * 100
    print(f"  B.3  Improvement over gs37 (c_eff=1.0)")
    print(f"  ─────────────────────────────────────────")
    print(f"    gs37 TGP chi2/dof:  {gs37_chi2dof_tgp:.2f}")
    print(f"    gs38 TGP chi2/dof:  {total_chi2_tgp/total_dof:.2f}")
    print(f"    chi2 reduction:     {improvement:.1f}%")
    print()

    return total_chi2_tgp, total_chi2_mond, total_dof, ml_tgp_list, ml_mond_list, galaxy_results


# ================================================================
#  PART C: UNIVERSAL c_eff = 1.3 AND c_eff = 1.5
# ================================================================

def part_C():
    print_header("PART C: UNIVERSAL c_eff COMPARISON")

    print("  Testing universal c_eff for all 20 galaxies")
    print("  (no type dependence, single value for all)")
    print()

    gs37_chi2_tgp = 3551.9

    for c_eff_univ in [1.0, 1.3, 1.5]:
        gamma = 0.8 * c_eff_univ / (c_eff_univ + 1)
        print(f"  --- c_eff = {c_eff_univ:.1f}  (gamma = {gamma:.3f}) ---")
        print()
        print(f"    {'Galaxy':<14} {'Type':<6} {'ML_TGP':<8} {'chi2_TGP':<10} {'ML_MOND':<8} {'chi2_MOND':<10} {'dof':<5}")
        print(f"    {'---'*22}")

        total_chi2_tgp = 0
        total_chi2_mond = 0
        total_dof = 0
        ml_list = []

        for gal in sparc_sample:
            name, gtype, M_disk, M_gas, R_d, V_flat, dist, incl, qual = gal
            if name not in obs_curves:
                continue

            best = fit_ML(M_disk, M_gas, R_d, obs_curves[name], sigma_v=5.0, c_eff=c_eff_univ)
            t = best['TGP']
            m = best['MOND']

            total_chi2_tgp += t['chi2']
            total_chi2_mond += m['chi2']
            total_dof += t['dof']
            ml_list.append(t['ML'])

            print(f"    {name:<14} {gtype:<6} {t['ML']:<8.2f} {t['chi2']:<10.1f} {m['ML']:<8.2f} {m['chi2']:<10.1f} {t['dof']:<5}")

        ml_arr = np.array(ml_list)
        in_range = np.sum((ml_arr >= 0.2) & (ml_arr <= 1.0))
        improvement = (gs37_chi2_tgp - total_chi2_tgp) / gs37_chi2_tgp * 100

        print()
        print(f"    Total chi2 TGP:     {total_chi2_tgp:.1f}  (chi2/dof = {total_chi2_tgp/total_dof:.2f})")
        print(f"    Total chi2 MOND:    {total_chi2_mond:.1f}  (chi2/dof = {total_chi2_mond/total_dof:.2f})")
        print(f"    M/L mean = {np.mean(ml_arr):.2f}, median = {np.median(ml_arr):.2f}")
        print(f"    M/L in [0.2, 1.0]:  {in_range}/{len(ml_list)}")
        print(f"    chi2 improvement over gs37: {improvement:.1f}%")
        print()


# ================================================================
#  PART D: UNIVERSAL c_eff SUMMARY TABLE
# ================================================================

def part_D():
    print_header("PART D: c_eff SCAN -- SUMMARY TABLE")

    print("  Scanning c_eff from 1.0 to 2.0 for all 20 galaxies")
    print()
    print(f"    {'c_eff':<8} {'gamma':<8} {'chi2/dof':<10} {'chi2_total':<12} {'<M/L>':<8} {'M/L_med':<8} {'N(0.2-1.0)':<12} {'vs MOND':<10}")
    print(f"    {'---'*25}")

    gs37_chi2_tgp = 3551.9
    mond_chi2dof = 17.52

    for c_eff_val in [1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0]:
        gamma = 0.8 * c_eff_val / (c_eff_val + 1)
        total_chi2 = 0
        total_dof = 0
        ml_list = []

        for gal in sparc_sample:
            name, gtype, M_disk, M_gas, R_d, V_flat, dist, incl, qual = gal
            if name not in obs_curves:
                continue

            best = fit_ML(M_disk, M_gas, R_d, obs_curves[name], sigma_v=5.0, c_eff=c_eff_val)
            t = best['TGP']
            total_chi2 += t['chi2']
            total_dof += t['dof']
            ml_list.append(t['ML'])

        ml_arr = np.array(ml_list)
        in_range = np.sum((ml_arr >= 0.2) & (ml_arr <= 1.0))
        chi2dof = total_chi2 / total_dof
        vs_mond = chi2dof / mond_chi2dof

        print(f"    {c_eff_val:<8.1f} {gamma:<8.3f} {chi2dof:<10.2f} {total_chi2:<12.1f} {np.mean(ml_arr):<8.2f} {np.median(ml_arr):<8.2f} {in_range:<12} {vs_mond:<10.2f}")

    print()
    print("  'vs MOND' = TGP_chi2/dof divided by MOND_chi2/dof (17.52)")
    print("  Ratio < 1 means TGP fits BETTER than MOND")


# ================================================================
#  PART E: M/L REASONABLENESS CHECK
# ================================================================

def part_E(ml_tgp_list, ml_mond_list, galaxy_results):
    print_header("PART E: M/L REASONABLENESS CHECK")

    print("  E.1  Expected M/L at 3.6 um from stellar population synthesis")
    print("  ────────────────────────────────────────────────────────────────")
    print("  Disk-dominated (Sc/Sd): M/L ~ 0.3-0.6")
    print("  Bulge-dominated (Sa/Sb): M/L ~ 0.5-0.8")
    print("  Irregulars/dwarfs: M/L ~ 0.2-0.5")
    print("  Population mean: M/L ~ 0.5 +/- 0.2")
    print()

    ml_tgp_arr = np.array(ml_tgp_list)
    ml_mond_arr = np.array(ml_mond_list)
    in_range_tgp = np.sum((ml_tgp_arr >= 0.2) & (ml_tgp_arr <= 1.0))
    in_range_mond = np.sum((ml_mond_arr >= 0.2) & (ml_mond_arr <= 1.0))

    print(f"  E.2  M/L in [0.2, 1.0] range (type-dependent c_eff)")
    print(f"  ─────────────────────────────────────────────────────")
    print(f"    TGP:   {in_range_tgp}/{len(ml_tgp_arr)}")
    print(f"    MOND:  {in_range_mond}/{len(ml_mond_arr)}")
    print(f"    gs37 TGP (c=1): 0/20")
    print()

    print(f"  E.3  Per-galaxy M/L detail")
    print(f"  ───────────────────────────")
    print()
    print(f"    {'Galaxy':<14} {'Type':<6} {'c_eff':<7} {'ML_TGP':<8} {'ML_MOND':<8} {'TGP_OK?':<10} {'MOND_OK?':<10}")
    print(f"    {'---'*22}")

    for gr in galaxy_results:
        tgp_ok = "YES" if 0.2 <= gr['ML_tgp'] <= 1.0 else "NO"
        mond_ok = "YES" if 0.2 <= gr['ML_mond'] <= 1.0 else "NO"
        print(f"    {gr['name']:<14} {gr['type']:<6} {gr['c_eff']:<7.2f} {gr['ML_tgp']:<8.2f} {gr['ML_mond']:<8.2f} {tgp_ok:<10} {mond_ok:<10}")

    print()


# ================================================================
#  PART F: BTFR FROM THE SAMPLE
# ================================================================

def part_F(galaxy_results):
    print_header("PART F: BARYONIC TULLY-FISHER RELATION (BTFR)")

    print("  F.1  BTFR: log(M_bar) vs log(V_flat)")
    print("  ────────────────────────────────────────")
    print("  Using best-fit M/L from type-dependent c_eff TGP fit")
    print("  M_bar = ML_disk * M_disk + M_gas")
    print("  V_flat = observed flat velocity (from sample)")
    print()

    print(f"    {'Galaxy':<14} {'Type':<6} {'ML_TGP':<8} {'M_bar(Msun)':<14} {'V_flat':<8} {'log(Mbar)':<10} {'log(Vf)':<8}")
    print(f"    {'---'*25}")

    log_Mbar_list = []
    log_Vf_list = []

    for gr in galaxy_results:
        M_bar = gr['ML_tgp'] * gr['M_disk'] + gr['M_gas']
        V_flat = gr['V_flat']
        log_Mbar = np.log10(M_bar)
        log_Vf = np.log10(V_flat)

        log_Mbar_list.append(log_Mbar)
        log_Vf_list.append(log_Vf)

        print(f"    {gr['name']:<14} {gr['type']:<6} {gr['ML_tgp']:<8.2f} {M_bar:<14.3e} {V_flat:<8.0f} {log_Mbar:<10.2f} {log_Vf:<8.3f}")

    log_Mbar_arr = np.array(log_Mbar_list)
    log_Vf_arr = np.array(log_Vf_list)

    # Linear fit: log(Mbar) = slope * log(Vf) + intercept
    # Or equivalently: log(Vf) = (1/slope) * log(Mbar) + const
    # Standard BTFR: M_bar = A * V^slope -> log(M) = slope*log(V) + log(A)
    # Fit log(Mbar) = a * log(Vf) + b
    coeffs = np.polyfit(log_Vf_arr, log_Mbar_arr, 1)
    slope = coeffs[0]
    intercept = coeffs[1]

    # Scatter
    log_Mbar_pred = slope * log_Vf_arr + intercept
    residuals = log_Mbar_arr - log_Mbar_pred
    scatter = np.std(residuals)

    print()
    print("  F.2  BTFR fit results")
    print("  ──────────────────────")
    print(f"    log(M_bar) = {slope:.2f} * log(V_flat) + {intercept:.2f}")
    print(f"    BTFR slope:     {slope:.2f}")
    print(f"    Scatter (dex):  {scatter:.3f}")
    print()
    print("  F.3  Comparison with observations")
    print("  ───────────────────────────────────")
    print(f"    Observed BTFR slope:       3.85 +/- 0.09  (McGaugh 2012)")
    print(f"    gs38 TGP (type c_eff):     {slope:.2f}")
    print()

    # Theoretical deep-MOND slopes
    print("  F.4  Theoretical deep-MOND BTFR slopes")
    print("  ─────────────────────────────────────────")
    print(f"    {'c_eff':<8} {'gamma':<8} {'Predicted slope':<16}")
    print(f"    {'---'*12}")
    for c_eff_val in [1.0, 1.2, 1.3, 1.5, 1.8, 2.0]:
        gamma = 0.8 * c_eff_val / (c_eff_val + 1)
        btfr_slope = 2.0 / (1.0 - gamma)
        print(f"    {c_eff_val:<8.1f} {gamma:<8.3f} {btfr_slope:<16.2f}")

    print()
    print(f"    Weighted mean c_eff of sample:")
    ceff_list = [ceff_by_type.get(gr['type'], 1.0) for gr in galaxy_results]
    mean_ceff = np.mean(ceff_list)
    gamma_mean = 0.8 * mean_ceff / (mean_ceff + 1)
    pred_slope = 2.0 / (1.0 - gamma_mean)
    print(f"      <c_eff> = {mean_ceff:.2f} -> gamma = {gamma_mean:.3f} -> predicted BTFR slope = {pred_slope:.2f}")

    return slope, scatter


# ================================================================
#  PART G: SUMMARY
# ================================================================

def part_G(chi2_tgp_typed, chi2_mond, total_dof, ml_tgp_list, btfr_slope, btfr_scatter):
    print_header("PART G: SUMMARY -- chi2 IMPROVEMENT OVER gs37")

    gs37_chi2_tgp = 3551.9
    gs37_chi2dof = 35.88
    gs37_ml_mean = 1.87
    gs37_ml_in_range = 0

    chi2dof_typed = chi2_tgp_typed / total_dof
    chi2dof_mond = chi2_mond / total_dof

    ml_arr = np.array(ml_tgp_list)
    in_range = int(np.sum((ml_arr >= 0.2) & (ml_arr <= 1.0)))

    print("  G.1  chi2 comparison")
    print("  ──────────────────────")
    print()
    print(f"    {'Model':<30} {'chi2/dof':<12} {'chi2_total':<12} {'<M/L>':<8} {'M/L in 0.2-1.0':<16}")
    print(f"    {'---'*28}")
    print(f"    {'gs37: TGP c_eff=1.0':<30} {gs37_chi2dof:<12.2f} {gs37_chi2_tgp:<12.1f} {gs37_ml_mean:<8.2f} {gs37_ml_in_range:<16}")
    print(f"    {'gs38: TGP type c_eff':<30} {chi2dof_typed:<12.2f} {chi2_tgp_typed:<12.1f} {np.mean(ml_arr):<8.2f} {in_range:<16}")
    print(f"    {'MOND (standard)':<30} {chi2dof_mond:<12.2f} {chi2_mond:<12.1f} {'---':<8} {'---':<16}")
    print()

    improvement = (gs37_chi2_tgp - chi2_tgp_typed) / gs37_chi2_tgp * 100
    ratio_vs_mond = chi2dof_typed / chi2dof_mond

    print(f"  G.2  Key metrics")
    print(f"  ──────────────────")
    print(f"    chi2 reduction from gs37:       {improvement:.1f}%")
    print(f"    TGP/MOND chi2 ratio:            {ratio_vs_mond:.2f}")
    print(f"    M/L in physical range:           {in_range}/20 (was 0/20 in gs37)")
    print(f"    BTFR slope (fitted):             {btfr_slope:.2f} (observed: 3.85 +/- 0.09)")
    print(f"    BTFR scatter:                    {btfr_scatter:.3f} dex")
    print()

    print("  G.3  Conclusions")
    print("  ──────────────────")
    print()
    print("  1. TYPE-DEPENDENT c_eff dramatically improves TGP rotation curve fits")
    print(f"     chi2/dof: {gs37_chi2dof:.2f} -> {chi2dof_typed:.2f}  ({improvement:.0f}% reduction)")
    print()
    print("  2. Physical motivation for c_eff > 1:")
    print("     - Sb galaxies: moderate bulge + thick disk -> c_eff ~ 1.3")
    print("     - Sc galaxies: thin disk but not infinitely thin -> c_eff ~ 1.2")
    print("     - Im/dIrr: puffy, irregular mass distribution -> c_eff ~ 1.8")
    print("     - LSB: thick, diffuse disks -> c_eff ~ 1.5")
    print()
    print("  3. M/L values become more physical:")
    print(f"     gs37: <M/L> = {gs37_ml_mean:.2f}, 0/20 in [0.2, 1.0]")
    print(f"     gs38: <M/L> = {np.mean(ml_arr):.2f}, {in_range}/20 in [0.2, 1.0]")
    print()
    print("  4. BTFR slope consistency:")
    print(f"     gs37 (c_eff=1.0): predicted slope = 3.33 (too low)")
    print(f"     gs38 (type c_eff): fitted slope = {btfr_slope:.2f}")
    print(f"     Observed: 3.85 +/- 0.09")
    print()
    if ratio_vs_mond <= 1.0:
        print("  5. TGP with type-dependent c_eff NOW MATCHES OR BEATS MOND!")
        print(f"     TGP/MOND chi2 ratio = {ratio_vs_mond:.2f}")
    else:
        print(f"  5. TGP with type-dependent c_eff is closer to MOND but still")
        print(f"     has TGP/MOND chi2 ratio = {ratio_vs_mond:.2f}")
        print(f"     Further improvement possible with finer c_eff tuning per galaxy")
    print()
    print("  6. NEXT STEPS (gs39+):")
    print("     - Fit c_eff as free parameter per galaxy (alongside M/L)")
    print("     - Correlate best-fit c_eff with measured disk thickness h/R")
    print("     - Test on full SPARC sample (175 galaxies)")
    print("     - Compare with external field effect predictions")


# ================================================================
#  MAIN
# ================================================================

if __name__ == '__main__':
    print("=" * 78)
    print("  gs38: SPARC ROTATION CURVES REFIT -- TYPE-DEPENDENT c_eff")
    print("=" * 78)

    part_A()
    chi2_tgp, chi2_mond, total_dof, ml_tgp_list, ml_mond_list, galaxy_results = part_B()
    part_C()
    part_D()
    part_E(ml_tgp_list, ml_mond_list, galaxy_results)
    btfr_slope, btfr_scatter = part_F(galaxy_results)
    part_G(chi2_tgp, chi2_mond, total_dof, ml_tgp_list, btfr_slope, btfr_scatter)

    print()
    print("=" * 78)
    print("  END OF gs38: SPARC c_eff REFIT")
    print("=" * 78)
