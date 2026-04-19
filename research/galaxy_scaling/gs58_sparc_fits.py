#!/usr/bin/env python3
"""
gs58_sparc_fits.py - SPARC rotation curve fits for TGP interpolation functions
Compares TGP standard, MOND simple, and Form 4 nu-functions on 20 benchmark galaxies.
"""

import numpy as np
import os
import sys

# ==============================================================================
# Constants
# ==============================================================================
a0 = 1.2e-10          # m/s^2, MOND acceleration scale
kpc_to_m = 3.0857e19  # meters per kpc
km_to_m = 1.0e3       # meters per km

ML_disk_fixed = 0.5   # M/L for disk at 3.6 micron
ML_bul_fixed = 0.7    # M/L for bulge at 3.6 micron

# ==============================================================================
# Interpolation functions nu(y) where y = g_bar / a0
# ==============================================================================

def nu_tgp(y):
    """TGP standard: nu(y) = 1 + exp(-y^0.8) / y^0.5714"""
    y = np.maximum(y, 1e-30)
    return 1.0 + np.exp(-y**0.8) / y**0.5714

def nu_mond(y):
    """MOND simple: nu(y) = 0.5 + sqrt(0.25 + 1/y)"""
    y = np.maximum(y, 1e-30)
    return 0.5 + np.sqrt(0.25 + 1.0/y)

def nu_form4(y):
    """Form 4 best-fit: nu(y) = 1/(1 - exp(-y^0.497))"""
    y = np.maximum(y, 1e-30)
    val = 1.0 - np.exp(-y**0.497)
    val = np.maximum(val, 1e-30)
    return 1.0 / val

# ==============================================================================
# Galaxy list: name -> filename mapping
# ==============================================================================

GALAXIES = {
    # Dwarf / low-mass (V_flat < 80 km/s)
    'DDO154':  'DDO154_rotmod.dat',
    'DDO168':  'DDO168_rotmod.dat',
    'NGC3741': 'NGC3741_rotmod.dat',
    'IC2574':  'IC2574_rotmod.dat',
    'NGC2366': 'NGC2366_rotmod.dat',
    # Intermediate (80-150 km/s)
    'NGC2403': 'NGC2403_rotmod.dat',
    'NGC3198': 'NGC3198_rotmod.dat',
    'NGC6503': 'NGC6503_rotmod.dat',
    'NGC2976': 'NGC2976_rotmod.dat',
    'NGC5585': 'NGC5585_rotmod.dat',
    # High-mass spirals (150-250 km/s)
    'NGC2903': 'NGC2903_rotmod.dat',
    'NGC3521': 'NGC3521_rotmod.dat',
    'NGC5055': 'NGC5055_rotmod.dat',
    'NGC6946': 'NGC6946_rotmod.dat',
    'NGC7331': 'NGC7331_rotmod.dat',
    # Massive / HSB (>250 km/s)
    'NGC2841': 'NGC2841_rotmod.dat',
    'NGC2998': 'NGC2998_rotmod.dat',
    'NGC5371': 'NGC5371_rotmod.dat',
    'NGC0801': 'NGC0801_rotmod.dat',
    'NGC7814': 'NGC7814_rotmod.dat',
}

MASS_CATEGORY = {
    'DDO154': 'dwarf', 'DDO168': 'dwarf', 'NGC3741': 'dwarf',
    'IC2574': 'dwarf', 'NGC2366': 'dwarf',
    'NGC2403': 'intermediate', 'NGC3198': 'intermediate', 'NGC6503': 'intermediate',
    'NGC2976': 'intermediate', 'NGC5585': 'intermediate',
    'NGC2903': 'high-mass', 'NGC3521': 'high-mass', 'NGC5055': 'high-mass',
    'NGC6946': 'high-mass', 'NGC7331': 'high-mass',
    'NGC2841': 'massive', 'NGC2998': 'massive', 'NGC5371': 'massive',
    'NGC0801': 'massive', 'NGC7814': 'massive',
}

SHOWCASE = ['NGC2403', 'DDO154', 'NGC2841', 'NGC7331']

# ==============================================================================
# Data loading
# ==============================================================================

def load_galaxy(filepath):
    """Load rotation curve data from SPARC file.
    Returns: R_kpc, Vobs, errV, Vgas, Vdisk, Vbul (all in original units)
    """
    R, Vobs, errV, Vgas, Vdisk, Vbul = [], [], [], [], [], []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or len(line) == 0:
                continue
            parts = line.split()
            if len(parts) < 6:
                continue
            R.append(float(parts[0]))
            Vobs.append(float(parts[1]))
            errV.append(float(parts[2]))
            Vgas.append(float(parts[3]))
            Vdisk.append(float(parts[4]))
            Vbul.append(float(parts[5]))
    return (np.array(R), np.array(Vobs), np.array(errV),
            np.array(Vgas), np.array(Vdisk), np.array(Vbul))


def compute_vbar_squared(Vgas, Vdisk, Vbul, ml_disk, ml_bul):
    """Compute V_bar^2 in (km/s)^2."""
    return (Vdisk**2 * ml_disk) + Vgas**2 + (Vbul**2 * ml_bul)


def compute_gbar(Vbar2_kms2, R_kpc):
    """Compute g_bar in m/s^2 from V_bar^2 (km/s)^2 and R (kpc)."""
    Vbar2_ms2 = Vbar2_kms2 * km_to_m**2   # (m/s)^2
    R_m = R_kpc * kpc_to_m                 # meters
    return Vbar2_ms2 / R_m                 # m/s^2


def compute_vpred(nu_func, gbar, R_kpc, Vbar2_kms2):
    """Compute predicted rotation velocity in km/s.
    V_pred^2 = nu(g_bar/a0) * V_bar^2  (in km/s)^2
    """
    y = gbar / a0
    nu_val = nu_func(y)
    vpred2 = nu_val * Vbar2_kms2   # (km/s)^2
    # protect against negative
    vpred2 = np.maximum(vpred2, 0.0)
    return np.sqrt(vpred2)


def compute_chi2(Vobs, Vpred, errV):
    """Compute chi2 = sum((Vobs - Vpred)^2 / errV^2)."""
    # Protect against zero errors
    err = np.maximum(errV, 0.1)
    return np.sum((Vobs - Vpred)**2 / err**2)


# ==============================================================================
# Main analysis
# ==============================================================================

def main():
    data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Rotmod_LTG')

    nu_functions = {
        'TGP':   nu_tgp,
        'MOND':  nu_mond,
        'Form4': nu_form4,
    }
    model_names = ['TGP', 'MOND', 'Form4']

    print("=" * 100)
    print("GS58: SPARC Rotation Curve Fits — TGP vs MOND vs Form 4")
    print("=" * 100)
    print()
    print(f"a0 = {a0:.2e} m/s^2")
    print(f"Fixed M/L_disk = {ML_disk_fixed},  M/L_bul = {ML_bul_fixed}")
    print(f"Number of galaxies: {len(GALAXIES)}")
    print()

    # =========================================================================
    # PART A: Load all galaxies
    # =========================================================================
    print("=" * 100)
    print("PART A: Benchmark Galaxy Sample")
    print("=" * 100)
    print()

    galaxy_data = {}
    for gname, fname in GALAXIES.items():
        fpath = os.path.join(data_dir, fname)
        if not os.path.exists(fpath):
            print(f"  WARNING: {fpath} not found, skipping {gname}")
            continue
        R, Vobs, errV, Vgas, Vdisk, Vbul = load_galaxy(fpath)
        galaxy_data[gname] = {
            'R': R, 'Vobs': Vobs, 'errV': errV,
            'Vgas': Vgas, 'Vdisk': Vdisk, 'Vbul': Vbul,
            'N': len(R),
        }
        print(f"  {gname:10s}  N={len(R):3d} points   R=[{R[0]:.2f}, {R[-1]:.2f}] kpc   "
              f"Vobs=[{Vobs.min():.1f}, {Vobs.max():.1f}] km/s   category={MASS_CATEGORY[gname]}")

    print(f"\n  Loaded {len(galaxy_data)} galaxies successfully.\n")

    # =========================================================================
    # PART B: Fixed M/L fits
    # =========================================================================
    print("=" * 100)
    print("PART B: Fixed M/L Fits (M/L_disk=0.5, M/L_bul=0.7)")
    print("=" * 100)
    print()

    results_fixed = {}  # gname -> {model: chi2_red}

    for gname in GALAXIES:
        if gname not in galaxy_data:
            continue
        d = galaxy_data[gname]
        R, Vobs, errV = d['R'], d['Vobs'], d['errV']
        Vgas, Vdisk, Vbul = d['Vgas'], d['Vdisk'], d['Vbul']
        N = d['N']

        Vbar2 = compute_vbar_squared(Vgas, Vdisk, Vbul, ML_disk_fixed, ML_bul_fixed)
        gbar = compute_gbar(Vbar2, R)

        results_fixed[gname] = {}
        for mname in model_names:
            Vpred = compute_vpred(nu_functions[mname], gbar, R, Vbar2)
            chi2 = compute_chi2(Vobs, Vpred, errV)
            chi2_red = chi2 / N
            results_fixed[gname][mname] = chi2_red

    # Print table
    header = f"{'Galaxy':12s} {'N':>4s}  {'cat':>14s}  {'chi2r_TGP':>10s}  {'chi2r_MOND':>10s}  {'chi2r_Form4':>11s}  {'Winner':>8s}"
    print(header)
    print("-" * len(header))

    for gname in GALAXIES:
        if gname not in results_fixed:
            continue
        d = galaxy_data[gname]
        r = results_fixed[gname]
        best = min(model_names, key=lambda m: r[m])
        print(f"{gname:12s} {d['N']:4d}  {MASS_CATEGORY[gname]:>14s}  "
              f"{r['TGP']:10.3f}  {r['MOND']:10.3f}  {r['Form4']:11.3f}  {best:>8s}")

    print()

    # =========================================================================
    # PART C: Free M/L fits (scan M/L_disk)
    # =========================================================================
    print("=" * 100)
    print("PART C: Free M/L_disk Fits (scan 0.1 to 1.5, step 0.01)")
    print("=" * 100)
    print()

    ml_scan = np.arange(0.1, 1.501, 0.01)

    results_free = {}   # gname -> {model: (best_ml, chi2_red)}

    for gname in GALAXIES:
        if gname not in galaxy_data:
            continue
        d = galaxy_data[gname]
        R, Vobs, errV = d['R'], d['Vobs'], d['errV']
        Vgas, Vdisk, Vbul = d['Vgas'], d['Vdisk'], d['Vbul']
        N = d['N']

        results_free[gname] = {}
        for mname in model_names:
            best_chi2 = np.inf
            best_ml = ML_disk_fixed
            for ml in ml_scan:
                Vbar2 = compute_vbar_squared(Vgas, Vdisk, Vbul, ml, ML_bul_fixed)
                gbar = compute_gbar(Vbar2, R)
                Vpred = compute_vpred(nu_functions[mname], gbar, R, Vbar2)
                chi2 = compute_chi2(Vobs, Vpred, errV)
                if chi2 < best_chi2:
                    best_chi2 = chi2
                    best_ml = ml
            chi2_red = best_chi2 / N
            results_free[gname][mname] = (best_ml, chi2_red)

    # Print table
    header = (f"{'Galaxy':12s} {'N':>4s}  "
              f"{'ML_TGP':>7s} {'chi2r':>7s}  "
              f"{'ML_MOND':>7s} {'chi2r':>7s}  "
              f"{'ML_F4':>7s} {'chi2r':>7s}  "
              f"{'Winner':>8s}")
    print(header)
    print("-" * len(header))

    for gname in GALAXIES:
        if gname not in results_free:
            continue
        d = galaxy_data[gname]
        r = results_free[gname]
        best = min(model_names, key=lambda m: r[m][1])
        print(f"{gname:12s} {d['N']:4d}  "
              f"{r['TGP'][0]:7.2f} {r['TGP'][1]:7.3f}  "
              f"{r['MOND'][0]:7.2f} {r['MOND'][1]:7.3f}  "
              f"{r['Form4'][0]:7.2f} {r['Form4'][1]:7.3f}  "
              f"{best:>8s}")

    print()

    # =========================================================================
    # PART D: Summary statistics
    # =========================================================================
    print("=" * 100)
    print("PART D: Summary Statistics")
    print("=" * 100)
    print()

    # --- Fixed M/L stats ---
    print("--- Fixed M/L (disk=0.5, bul=0.7) ---")
    for mname in model_names:
        vals = [results_fixed[g][mname] for g in results_fixed]
        print(f"  {mname:6s}:  mean chi2_red = {np.mean(vals):8.3f}   "
              f"median = {np.median(vals):8.3f}   std = {np.std(vals):8.3f}")

    wins_fixed = {m: 0 for m in model_names}
    for g in results_fixed:
        best = min(model_names, key=lambda m: results_fixed[g][m])
        wins_fixed[best] += 1
    print(f"\n  Wins (lowest chi2_red):  TGP={wins_fixed['TGP']}  MOND={wins_fixed['MOND']}  Form4={wins_fixed['Form4']}")
    print()

    # --- Free M/L stats ---
    print("--- Free M/L_disk ---")
    for mname in model_names:
        chi2_vals = [results_free[g][mname][1] for g in results_free]
        ml_vals = [results_free[g][mname][0] for g in results_free]
        print(f"  {mname:6s}:  mean chi2_red = {np.mean(chi2_vals):8.3f}   "
              f"median = {np.median(chi2_vals):8.3f}   std = {np.std(chi2_vals):8.3f}")
        print(f"          M/L_disk: mean = {np.mean(ml_vals):5.2f}   "
              f"std = {np.std(ml_vals):5.2f}   range = [{np.min(ml_vals):.2f}, {np.max(ml_vals):.2f}]")

    wins_free = {m: 0 for m in model_names}
    for g in results_free:
        best = min(model_names, key=lambda m: results_free[g][m][1])
        wins_free[best] += 1
    print(f"\n  Wins (lowest chi2_red):  TGP={wins_free['TGP']}  MOND={wins_free['MOND']}  Form4={wins_free['Form4']}")
    print()

    # M/L scatter by category
    print("--- Best-fit M/L_disk by mass category (Free M/L) ---")
    for cat in ['dwarf', 'intermediate', 'high-mass', 'massive']:
        gals = [g for g in results_free if MASS_CATEGORY[g] == cat]
        if not gals:
            continue
        print(f"\n  {cat}:")
        for mname in model_names:
            mls = [results_free[g][mname][0] for g in gals]
            print(f"    {mname:6s}: {', '.join(f'{v:.2f}' for v in mls)}  "
                  f"  mean={np.mean(mls):.2f}  std={np.std(mls):.2f}")
    print()

    # =========================================================================
    # PART E: Detailed comparison for 4 showcase galaxies
    # =========================================================================
    print("=" * 100)
    print("PART E: Detailed Rotation Curve Comparison (4 Showcase Galaxies)")
    print("=" * 100)

    for gname in SHOWCASE:
        if gname not in galaxy_data:
            continue
        d = galaxy_data[gname]
        R, Vobs, errV = d['R'], d['Vobs'], d['errV']
        Vgas, Vdisk, Vbul = d['Vgas'], d['Vdisk'], d['Vbul']
        N = d['N']

        Vbar2 = compute_vbar_squared(Vgas, Vdisk, Vbul, ML_disk_fixed, ML_bul_fixed)
        gbar = compute_gbar(Vbar2, R)

        Vpred = {}
        for mname in model_names:
            Vpred[mname] = compute_vpred(nu_functions[mname], gbar, R, Vbar2)

        print(f"\n  Galaxy: {gname}  (N={N}, category={MASS_CATEGORY[gname]})")
        print(f"  Fixed M/L_disk={ML_disk_fixed}, M/L_bul={ML_bul_fixed}")
        print(f"  chi2_red:  TGP={results_fixed[gname]['TGP']:.3f}  "
              f"MOND={results_fixed[gname]['MOND']:.3f}  "
              f"Form4={results_fixed[gname]['Form4']:.3f}")
        print()
        hdr = f"  {'R(kpc)':>8s}  {'Vobs':>7s}  {'errV':>5s}  {'V_TGP':>7s}  {'V_MOND':>7s}  {'V_Form4':>7s}  {'Vbar':>7s}"
        print(hdr)
        print("  " + "-" * (len(hdr) - 2))
        for i in range(N):
            Vbar = np.sqrt(max(Vbar2[i], 0))
            print(f"  {R[i]:8.2f}  {Vobs[i]:7.1f}  {errV[i]:5.1f}  "
                  f"{Vpred['TGP'][i]:7.1f}  {Vpred['MOND'][i]:7.1f}  "
                  f"{Vpred['Form4'][i]:7.1f}  {Vbar:7.1f}")
        print()

    # =========================================================================
    # PART F: Implications
    # =========================================================================
    print("=" * 100)
    print("PART F: Implications and Analysis")
    print("=" * 100)
    print()

    # F1: Does TGP systematically over/under-predict at certain mass scales?
    print("--- F1: TGP systematic trends by mass scale ---")
    for cat in ['dwarf', 'intermediate', 'high-mass', 'massive']:
        gals = [g for g in results_fixed if MASS_CATEGORY[g] == cat]
        if not gals:
            continue
        tgp_vals = [results_fixed[g]['TGP'] for g in gals]
        mond_vals = [results_fixed[g]['MOND'] for g in gals]
        form4_vals = [results_fixed[g]['Form4'] for g in gals]
        print(f"\n  {cat:14s}:  TGP mean chi2r = {np.mean(tgp_vals):7.3f}  "
              f"MOND = {np.mean(mond_vals):7.3f}  Form4 = {np.mean(form4_vals):7.3f}")

        # Check sign of residuals for TGP
        for g in gals:
            d = galaxy_data[g]
            Vbar2 = compute_vbar_squared(d['Vgas'], d['Vdisk'], d['Vbul'],
                                         ML_disk_fixed, ML_bul_fixed)
            gbar = compute_gbar(Vbar2, d['R'])
            Vpred = compute_vpred(nu_tgp, gbar, d['R'], Vbar2)
            resid = d['Vobs'] - Vpred
            mean_resid = np.mean(resid)
            frac_over = np.sum(resid < 0) / len(resid) * 100
            print(f"    {g:10s}: mean residual = {mean_resid:+6.1f} km/s   "
                  f"({frac_over:.0f}% over-predicted)")

    print()

    # F2: Form 4 vs TGP standard comparison
    print("--- F2: Is Form 4 (delta=0.497) consistently better than TGP (alpha=0.8)? ---")
    f4_better_fixed = sum(1 for g in results_fixed if results_fixed[g]['Form4'] < results_fixed[g]['TGP'])
    f4_better_free = sum(1 for g in results_free if results_free[g]['Form4'][1] < results_free[g]['TGP'][1])
    ntot = len(results_fixed)
    print(f"  Fixed M/L: Form4 better in {f4_better_fixed}/{ntot} galaxies")
    print(f"  Free  M/L: Form4 better in {f4_better_free}/{ntot} galaxies")

    ratio_fixed = []
    for g in results_fixed:
        ratio_fixed.append(results_fixed[g]['Form4'] / max(results_fixed[g]['TGP'], 1e-10))
    print(f"  Mean chi2r ratio (Form4/TGP, fixed M/L): {np.mean(ratio_fixed):.3f}")
    print(f"  Median chi2r ratio: {np.median(ratio_fixed):.3f}")
    print()

    # F3: 1.8-sigma delta tension
    print("--- F3: The 1.8-sigma delta tension (0.497 vs 0.8) ---")
    print("  Comparing fixed-M/L chi2_red differences (TGP - Form4):")
    for cat in ['dwarf', 'intermediate', 'high-mass', 'massive']:
        gals = [g for g in results_fixed if MASS_CATEGORY[g] == cat]
        diffs = [results_fixed[g]['TGP'] - results_fixed[g]['Form4'] for g in gals]
        print(f"    {cat:14s}:  mean(TGP - Form4) = {np.mean(diffs):+7.3f}   "
              f"range = [{np.min(diffs):+.3f}, {np.max(diffs):+.3f}]")
    all_diffs = [results_fixed[g]['TGP'] - results_fixed[g]['Form4'] for g in results_fixed]
    print(f"    {'ALL':14s}:  mean(TGP - Form4) = {np.mean(all_diffs):+7.3f}")
    print(f"  Positive values = Form4 is better (lower chi2r)")
    print()

    # F4: M/L values
    print("--- F4: Best-fit M/L_disk by model (suspicious if <0.3 or >0.8) ---")
    for mname in model_names:
        mls = [results_free[g][mname][0] for g in results_free]
        suspicious = [g for g in results_free
                      if results_free[g][mname][0] < 0.3 or results_free[g][mname][0] > 0.8]
        print(f"\n  {mname:6s}: mean M/L = {np.mean(mls):.2f}  std = {np.std(mls):.2f}")
        if suspicious:
            print(f"    Suspicious M/L galaxies:")
            for g in suspicious:
                ml = results_free[g][mname][0]
                flag = "LOW" if ml < 0.3 else "HIGH"
                print(f"      {g:12s}: M/L = {ml:.2f}  ({flag})")
        else:
            print(f"    No suspicious M/L values.")

    print()
    print("=" * 100)
    print("END OF GS58 ANALYSIS")
    print("=" * 100)


if __name__ == '__main__':
    main()
