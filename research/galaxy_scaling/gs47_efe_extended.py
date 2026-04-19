#!/usr/bin/env python3
"""
gs47_efe_extended.py
====================
Extended statistical analysis of TGP vs MOND External Field Effect.

Extends gs45_efe_discriminator.py from 4 systems to 8, adding:
  - Andromeda XIX (And XIX)
  - Tucana III
  - Segue 1
  - Antlia 2

Performs:
  A. Extended system catalog with literature values
  B. Systematic TGP vs MOND comparison for every system
  C. Global chi-squared statistical analysis
  D. Sensitivity analysis: d(sigma)/d(y_ext) and discriminating power
  E. Future predictions: required precision for 3-sigma discrimination

Physics:
  TGP  EFE: y_eff = sqrt(y_int^2 + y_ext^2)   (quadrature)
  MOND EFE: y_eff = y_int + y_ext              (linear)

  Interpolating function: nu(y) = 1 + exp(-y^(4/5)) / y^gamma
  gamma = 0.8 * c_eff / (c_eff + 1),  c_eff = 2.5 => gamma = 0.5714

  Velocity dispersion (Wolf estimator):
    sigma = sqrt(G * M_bar * nu(y_eff) / (4 * r_half))

Author : TGP collaboration
Date   : 2026-04-19
"""

import sys
sys.stdout.reconfigure(encoding='utf-8', errors='replace')

import numpy as np

# ============================================================
# Physical constants (SI)
# ============================================================
G_SI   = 6.674e-11        # m^3 kg^-1 s^-2
a0     = 1.2e-10           # m/s^2  (MOND acceleration scale)
Msun   = 1.989e30          # kg
kpc    = 3.086e19           # m
pc     = 3.086e16           # m

# Handy conversion
G_pc   = 4.302e-3          # pc (km/s)^2 / Msun  (not used directly, for reference)

# ============================================================
# TGP parameters
# ============================================================
C_EFF  = 2.5
GAMMA  = 0.8 * C_EFF / (C_EFF + 1.0)   # 0.5714...


# ============================================================
# Core physics functions
# ============================================================
def nu(y, gamma=GAMMA):
    """Interpolating function: nu(y) = 1 + exp(-y^(4/5)) / y^gamma."""
    return 1.0 + np.exp(-y**(4.0/5.0)) / y**gamma


def y_eff_mond(y_int, y_ext):
    """MOND: linear addition => strong EFE suppression."""
    return y_int + y_ext


def y_eff_tgp(y_int, y_ext):
    """TGP: quadrature addition => weaker EFE suppression."""
    return np.sqrt(y_int**2 + y_ext**2)


def compute_y_ext(M_host_Msun, D_kpc):
    """External field: y_ext = G*M_host / (D^2 * a0)."""
    D_m  = D_kpc * kpc
    M_kg = M_host_Msun * Msun
    g_ext = G_SI * M_kg / D_m**2
    return g_ext / a0


def compute_y_int(M_bar_Msun, r_half_pc):
    """Internal field: y_int = G*M_bar / (r_half^2 * a0)."""
    M_kg    = M_bar_Msun * Msun
    r_m     = r_half_pc * pc
    g_int   = G_SI * M_kg / r_m**2
    return g_int / a0


def sigma_predict(M_bar_Msun, r_half_pc, y_ext, model='tgp'):
    """Predict velocity dispersion using Wolf mass estimator.

    Wolf estimator: M_1/2 = 4 * sigma^2 * r_half / G
    With M_dyn = nu(y_eff) * M_bar as the dynamical mass:
      sigma = sqrt(G * M_bar * nu(y_eff) / (4 * r_half))

    Returns (sigma_km_s, y_int, y_eff_value, nu_value).
    """
    M_kg   = M_bar_Msun * Msun
    r_m    = r_half_pc * pc
    y_int  = compute_y_int(M_bar_Msun, r_half_pc)

    if model == 'mond':
        ye = y_eff_mond(y_int, y_ext)
    elif model == 'tgp':
        ye = y_eff_tgp(y_int, y_ext)
    elif model == 'newtonian':
        ye = 1e10  # huge y => nu ~ 1
    else:
        raise ValueError(f"Unknown model: {model}")

    nu_val = nu(ye)
    # Wolf estimator: sigma^2 = G * M_dyn / (4 * r_half)
    sigma_ms = np.sqrt(G_SI * nu_val * M_kg / (4.0 * r_m))
    return sigma_ms / 1e3, y_int, ye, nu_val


# ============================================================
# A. Extended System Catalog
# ============================================================
# Host masses
MW_MASS     = 1.5e12   # Msun
M31_MASS    = 1.5e12
NGC1052_MASS = 3e11

systems = [
    dict(
        name        = "Crater II",
        M_bar_Msun  = 1.6e5,
        r_half_pc   = 1066.0,
        sigma_obs   = 2.7,
        sigma_err   = 0.3,
        host_mass   = MW_MASS,
        D_kpc       = 117.0,
        is_upper    = False,
        notes       = "Ultra-diffuse MW dwarf; y_ext ~ y_int => maximal EFE sensitivity",
    ),
    dict(
        name        = "NGC 1052-DF2",
        M_bar_Msun  = 2e8,
        r_half_pc   = 2200.0,
        sigma_obs   = 8.5,
        sigma_err   = 2.3,
        host_mass   = NGC1052_MASS,
        D_kpc       = 80.0,
        is_upper    = False,
        notes       = "Galaxy 'lacking dark matter'; low sigma if confirmed",
    ),
    dict(
        name        = "NGC 1052-DF4",
        M_bar_Msun  = 1.5e8,
        r_half_pc   = 1600.0,
        sigma_obs   = 4.2,
        sigma_err   = 2.2,
        host_mass   = NGC1052_MASS,
        D_kpc       = 80.0,
        is_upper    = False,
        notes       = "Similar to DF2; very low sigma if confirmed",
    ),
    dict(
        name        = "Fornax",
        M_bar_Msun  = 4.3e7,
        r_half_pc   = 710.0,
        sigma_obs   = 11.7,
        sigma_err   = 0.9,
        host_mass   = MW_MASS,
        D_kpc       = 147.0,
        is_upper    = False,
        notes       = "Classical MW dSph; well-measured kinematics",
    ),
    dict(
        name        = "Andromeda XIX",
        M_bar_Msun  = 5e5,
        r_half_pc   = 1700.0,
        sigma_obs   = 4.7,
        sigma_err   = 1.5,
        host_mass   = M31_MASS,
        D_kpc       = 120.0,
        is_upper    = False,
        notes       = "Ultra-diffuse M31 satellite; large r_half => low y_int",
    ),
    dict(
        name        = "Tucana III",
        M_bar_Msun  = 800.0,
        r_half_pc   = 44.0,
        sigma_obs   = 1.5,
        sigma_err   = 0.5,       # treat upper limit as measurement with generous error
        host_mass   = MW_MASS,
        D_kpc       = 25.0,
        is_upper    = True,       # sigma is an upper limit
        notes       = "Ultra-faint MW stream; sigma < 1.5 km/s (upper limit)",
    ),
    dict(
        name        = "Segue 1",
        M_bar_Msun  = 1000.0,
        r_half_pc   = 29.0,
        sigma_obs   = 3.9,
        sigma_err   = 0.8,
        host_mass   = MW_MASS,
        D_kpc       = 23.0,
        is_upper    = False,
        notes       = "Ultra-faint MW dwarf; strong EFE due to proximity to MW",
    ),
    dict(
        name        = "Antlia 2",
        M_bar_Msun  = 1e6,
        r_half_pc   = 2900.0,
        sigma_obs   = 5.7,
        sigma_err   = 1.1,
        host_mass   = MW_MASS,
        D_kpc       = 132.0,
        is_upper    = False,
        notes       = "Giant ultra-diffuse MW satellite; extremely low surface brightness",
    ),
]


# ============================================================
# Analysis functions
# ============================================================
def analyse_system(s):
    """Compute Newtonian, MOND, TGP predictions for one system."""
    y_ext = compute_y_ext(s['host_mass'], s['D_kpc'])
    y_int = compute_y_int(s['M_bar_Msun'], s['r_half_pc'])

    results = {}
    for model in ('newtonian', 'mond', 'tgp'):
        sigma, yi, ye, nv = sigma_predict(
            s['M_bar_Msun'], s['r_half_pc'], y_ext, model=model
        )
        results[model] = dict(sigma=sigma, y_int=yi, y_eff=ye, nu=nv)

    return results, y_ext, y_int


def d_sigma_d_yext(M_bar_Msun, r_half_pc, y_ext, model, dy=1e-4):
    """Numerical derivative d(sigma)/d(y_ext) for sensitivity analysis."""
    # Use central difference
    h = max(dy * y_ext, 1e-6)
    s_plus, _, _, _  = sigma_predict(M_bar_Msun, r_half_pc, y_ext + h, model)
    s_minus, _, _, _ = sigma_predict(M_bar_Msun, r_half_pc, y_ext - h, model)
    return (s_plus - s_minus) / (2.0 * h)


# ============================================================
# Main
# ============================================================
def main():
    sep = "=" * 80
    dash = "-" * 80

    print(sep)
    print("  gs47_efe_extended.py")
    print("  Extended Statistical Analysis: TGP vs MOND External Field Effect")
    print("  8 dwarf galaxy systems -- chi-squared comparison")
    print(sep)

    # ===========================================================
    # PART A: Extended System Catalog
    # ===========================================================
    print(f"\n{'='*80}")
    print("  PART A: Extended System Catalog")
    print(f"{'='*80}\n")

    print(f"  TGP parameters: c_eff = {C_EFF:.1f}, gamma = {GAMMA:.4f}")
    print(f"  nu(y) = 1 + exp(-y^(4/5)) / y^gamma")
    print(f"  Wolf estimator: sigma = sqrt(G * M_bar * nu / (4 * r_half))")
    print(f"  TGP:  y_eff = sqrt(y_int^2 + y_ext^2)")
    print(f"  MOND: y_eff = y_int + y_ext\n")

    print(f"  Host masses:")
    print(f"    MW:       {MW_MASS:.2e} Msun")
    print(f"    M31:      {M31_MASS:.2e} Msun")
    print(f"    NGC 1052: {NGC1052_MASS:.2e} Msun\n")

    print(f"  {'#':<3s} {'Name':<18s} {'M_bar[Msun]':>12s} {'r_half[pc]':>10s} "
          f"{'sigma_obs':>10s} {'sigma_err':>9s} {'D_host[kpc]':>11s} {'Upper?':>6s}")
    print(f"  {'-'*3} {'-'*18} {'-'*12} {'-'*10} {'-'*10} {'-'*9} {'-'*11} {'-'*6}")

    for i, s in enumerate(systems, 1):
        ul = "Yes" if s['is_upper'] else "No"
        print(f"  {i:<3d} {s['name']:<18s} {s['M_bar_Msun']:>12.2e} {s['r_half_pc']:>10.0f} "
              f"{s['sigma_obs']:>10.1f} {s['sigma_err']:>9.1f} {s['D_kpc']:>11.0f} {ul:>6s}")

    # ===========================================================
    # PART B: Systematic TGP vs MOND comparison
    # ===========================================================
    print(f"\n{'='*80}")
    print("  PART B: Systematic TGP vs MOND Comparison")
    print(f"{'='*80}")

    all_results = []

    for s in systems:
        results, y_ext, y_int = analyse_system(s)
        sig_N    = results['newtonian']['sigma']
        sig_M    = results['mond']['sigma']
        sig_T    = results['tgp']['sigma']
        sig_obs  = s['sigma_obs']
        sig_err  = s['sigma_err']

        delta_T  = abs(sig_T - sig_obs)
        delta_M  = abs(sig_M - sig_obs)
        nsig_T   = delta_T / sig_err
        nsig_M   = delta_M / sig_err

        if delta_T < delta_M:
            closer = "TGP"
        elif delta_M < delta_T:
            closer = "MOND"
        else:
            closer = "TIE"

        row = dict(
            name     = s['name'],
            M_bar    = s['M_bar_Msun'],
            r_half   = s['r_half_pc'],
            y_int    = y_int,
            y_ext    = y_ext,
            ratio    = y_ext / y_int if y_int > 0 else np.inf,
            sig_N    = sig_N,
            sig_M    = sig_M,
            sig_T    = sig_T,
            sig_obs  = sig_obs,
            sig_err  = sig_err,
            delta_T  = delta_T,
            delta_M  = delta_M,
            nsig_T   = nsig_T,
            nsig_M   = nsig_M,
            closer   = closer,
            split    = abs(sig_T - sig_M),
            is_upper = s['is_upper'],
            notes    = s['notes'],
            nu_tgp   = results['tgp']['nu'],
            nu_mond  = results['mond']['nu'],
        )
        all_results.append(row)

        print(f"\n  {dash}")
        print(f"  System: {s['name']}")
        print(f"  {dash}")
        print(f"    M_bar       = {s['M_bar_Msun']:.2e} Msun")
        print(f"    r_half      = {s['r_half_pc']:.0f} pc")
        print(f"    sigma_obs   = {sig_obs:.1f} +/- {sig_err:.1f} km/s"
              + (" (UPPER LIMIT)" if s['is_upper'] else ""))
        print(f"    Host dist   = {s['D_kpc']:.0f} kpc,  Host mass = {s['host_mass']:.2e} Msun")
        print(f"    y_int       = {y_int:.6f}")
        print(f"    y_ext       = {y_ext:.6f}")
        print(f"    y_ext/y_int = {row['ratio']:.3f}")
        print(f"    Notes: {s['notes']}")
        print()

        for model, label in [('newtonian', 'NEWTON'), ('mond', 'MOND'), ('tgp', 'TGP')]:
            r = results[model]
            print(f"      {label:<8s}  y_eff = {r['y_eff']:12.6f}   "
                  f"nu = {r['nu']:10.4f}   sigma = {r['sigma']:7.2f} km/s")

        print()
        print(f"    |TGP  - obs| = {delta_T:.2f} km/s  ({nsig_T:.1f} sigma)")
        print(f"    |MOND - obs| = {delta_M:.2f} km/s  ({nsig_M:.1f} sigma)")
        print(f"    TGP-MOND split = {row['split']:.2f} km/s")
        print(f"    Closer model: {closer}")

    # Summary table
    print(f"\n  {sep}")
    print(f"  PART B SUMMARY TABLE")
    print(f"  {sep}\n")

    hdr = (f"  {'System':<18s} {'sig_N':>6s} {'sig_M':>6s} {'sig_T':>6s} "
           f"{'sig_obs':>7s} {'err':>5s} {'split':>6s} "
           f"{'nsig_T':>6s} {'nsig_M':>6s} {'closer':>7s}")
    print(hdr)
    print(f"  {'-'*len(hdr)}")

    for r in all_results:
        print(f"  {r['name']:<18s} {r['sig_N']:6.2f} {r['sig_M']:6.2f} {r['sig_T']:6.2f} "
              f"{r['sig_obs']:7.1f} {r['sig_err']:5.1f} {r['split']:6.2f} "
              f"{r['nsig_T']:6.1f} {r['nsig_M']:6.1f} {r['closer']:>7s}")

    print(f"\n  All sigma values in km/s.")

    # ===========================================================
    # PART C: Statistical Analysis (Chi-squared)
    # ===========================================================
    print(f"\n{'='*80}")
    print("  PART C: Chi-Squared Statistical Analysis")
    print(f"{'='*80}\n")

    # Compute chi-squared for each model
    chi2_tgp  = 0.0
    chi2_mond = 0.0
    chi2_newt = 0.0
    n_used    = 0

    print(f"  {'System':<18s} {'chi2_N':>8s} {'chi2_M':>8s} {'chi2_T':>8s} {'note':>12s}")
    print(f"  {'-'*18} {'-'*8} {'-'*8} {'-'*8} {'-'*12}")

    for r in all_results:
        c2_n = ((r['sig_N'] - r['sig_obs']) / r['sig_err'])**2
        c2_m = ((r['sig_M'] - r['sig_obs']) / r['sig_err'])**2
        c2_t = ((r['sig_T'] - r['sig_obs']) / r['sig_err'])**2

        note = ""
        if r['is_upper']:
            # For upper limits: only penalize if prediction exceeds observed
            if r['sig_N'] <= r['sig_obs']:
                c2_n = 0.0
            if r['sig_M'] <= r['sig_obs']:
                c2_m = 0.0
            if r['sig_T'] <= r['sig_obs']:
                c2_t = 0.0
            note = "(upper lim)"

        chi2_newt += c2_n
        chi2_mond += c2_m
        chi2_tgp  += c2_t
        n_used    += 1

        print(f"  {r['name']:<18s} {c2_n:8.2f} {c2_m:8.2f} {c2_t:8.2f} {note:>12s}")

    print(f"  {'-'*60}")
    print(f"  {'TOTAL':<18s} {chi2_newt:8.2f} {chi2_mond:8.2f} {chi2_tgp:8.2f}")
    print(f"  {'REDUCED (N='+ str(n_used) +')':<18s} "
          f"{chi2_newt/n_used:8.2f} {chi2_mond/n_used:8.2f} {chi2_tgp/n_used:8.2f}")

    dchi2 = chi2_tgp - chi2_mond
    print(f"\n  Delta_chi2 = chi2_TGP - chi2_MOND = {dchi2:.2f}")

    if dchi2 < 0:
        winner = "TGP"
        print(f"  => TGP is statistically preferred (lower total chi2)")
    elif dchi2 > 0:
        winner = "MOND"
        print(f"  => MOND is statistically preferred (lower total chi2)")
    else:
        winner = "TIE"
        print(f"  => Models are statistically tied")

    print(f"\n  Interpretation:")
    print(f"    N_systems = {n_used}")
    print(f"    chi2_Newton = {chi2_newt:.2f}  (reduced = {chi2_newt/n_used:.2f})")
    print(f"    chi2_MOND   = {chi2_mond:.2f}  (reduced = {chi2_mond/n_used:.2f})")
    print(f"    chi2_TGP    = {chi2_tgp:.2f}  (reduced = {chi2_tgp/n_used:.2f})")
    print(f"    Overall winner: {winner} (Delta_chi2 = {dchi2:+.2f})")

    # Per-system winner tally
    tgp_wins  = sum(1 for r in all_results if r['closer'] == 'TGP')
    mond_wins = sum(1 for r in all_results if r['closer'] == 'MOND')
    ties      = sum(1 for r in all_results if r['closer'] == 'TIE')
    print(f"\n  Per-system closer model tally:")
    print(f"    TGP:  {tgp_wins}")
    print(f"    MOND: {mond_wins}")
    print(f"    TIE:  {ties}")

    # ===========================================================
    # PART D: Sensitivity Analysis
    # ===========================================================
    print(f"\n{'='*80}")
    print("  PART D: Sensitivity Analysis")
    print(f"{'='*80}\n")

    print("  d(sigma)/d(y_ext) measures how much the predicted dispersion changes")
    print("  with the external field strength. Systems with large |d(sigma)/d(y_ext)|")
    print("  AND large TGP-MOND split are the best discriminators.\n")

    sens_data = []

    for i, s in enumerate(systems):
        y_ext = compute_y_ext(s['host_mass'], s['D_kpc'])
        r = all_results[i]

        ds_dy_tgp  = d_sigma_d_yext(s['M_bar_Msun'], s['r_half_pc'], y_ext, 'tgp')
        ds_dy_mond = d_sigma_d_yext(s['M_bar_Msun'], s['r_half_pc'], y_ext, 'mond')

        # Discriminating power: split / sigma_err
        disc_power = r['split'] / r['sig_err']

        sens_data.append(dict(
            name       = s['name'],
            y_ext      = y_ext,
            y_int      = r['y_int'],
            ratio      = r['ratio'],
            ds_dy_tgp  = ds_dy_tgp,
            ds_dy_mond = ds_dy_mond,
            ds_dy_diff = abs(ds_dy_tgp - ds_dy_mond),
            split      = r['split'],
            sig_err    = r['sig_err'],
            disc_power = disc_power,
        ))

    print(f"  {'System':<18s} {'y_ext/y_int':>11s} {'d_sig/dy_T':>10s} {'d_sig/dy_M':>10s} "
          f"{'split':>6s} {'sig_err':>7s} {'disc_pwr':>8s}")
    print(f"  {'-'*18} {'-'*11} {'-'*10} {'-'*10} {'-'*6} {'-'*7} {'-'*8}")

    for sd in sens_data:
        print(f"  {sd['name']:<18s} {sd['ratio']:11.3f} {sd['ds_dy_tgp']:10.4f} "
              f"{sd['ds_dy_mond']:10.4f} {sd['split']:6.2f} {sd['sig_err']:7.1f} "
              f"{sd['disc_power']:8.3f}")

    # Rank by discriminating power
    ranked = sorted(sens_data, key=lambda x: x['disc_power'], reverse=True)

    print(f"\n  Ranked by discriminating power (split / sigma_err):\n")
    for rank, sd in enumerate(ranked, 1):
        label = "*** BEST ***" if rank == 1 else ""
        print(f"    {rank}. {sd['name']:<18s}  disc_power = {sd['disc_power']:.3f}  "
              f"(split = {sd['split']:.2f}, err = {sd['sig_err']:.1f})  {label}")

    # Plot-ready data: split vs y_ext/y_int
    print(f"\n  Plot-ready data: TGP-MOND split vs y_ext/y_int ratio\n")
    print(f"  {'y_ext/y_int':>11s}  {'split[km/s]':>11s}  {'System':<18s}")
    print(f"  {'-'*11}  {'-'*11}  {'-'*18}")
    for sd in sorted(sens_data, key=lambda x: x['ratio']):
        print(f"  {sd['ratio']:11.3f}  {sd['split']:11.2f}  {sd['name']:<18s}")

    print(f"\n  Key insight: discriminating power peaks when y_ext/y_int ~ 1")
    print(f"  (quadrature vs linear differ most when arguments are comparable)")

    # ===========================================================
    # PART E: Future Predictions
    # ===========================================================
    print(f"\n{'='*80}")
    print("  PART E: Future Predictions")
    print(f"{'='*80}\n")

    print("  Question: What precision on sigma is needed to discriminate TGP from MOND")
    print("  at 3-sigma significance for each system?\n")

    print(f"  {'System':<18s} {'split':>6s} {'need_err':>8s} {'curr_err':>8s} "
          f"{'improvement':>11s} {'feasible?':>9s}")
    print(f"  {'-'*18} {'-'*6} {'-'*8} {'-'*8} {'-'*11} {'-'*9}")

    for i, r in enumerate(all_results):
        split = r['split']
        # For 3-sigma discrimination: split / sigma_err_needed = 3
        if split > 0:
            err_needed = split / 3.0
        else:
            err_needed = np.inf

        curr_err = r['sig_err']
        if err_needed < curr_err and err_needed > 0:
            improvement = curr_err / err_needed
            feasible = "Yes" if improvement < 10 else "Difficult"
        elif err_needed >= curr_err:
            improvement = 1.0
            feasible = "ALREADY"
        else:
            improvement = np.inf
            feasible = "No"

        # Check if already discriminating at 3-sigma
        if split / curr_err >= 3.0:
            feasible = "DONE (3sig)"

        print(f"  {r['name']:<18s} {split:6.2f} {err_needed:8.3f} {curr_err:8.1f} "
              f"{improvement:11.1f}x {feasible:>9s}")

    print(f"\n  Upcoming surveys and instruments relevant to EFE discrimination:\n")
    surveys = [
        ("LSST/Rubin Obs.", "2025-2035",
         "Discovery of ~100s new ultra-faint dwarfs; photometric masses"),
        ("DESI", "2024-2029",
         "Spectroscopic velocities for faint dwarfs; sigma to ~1 km/s precision"),
        ("4MOST", "2025-2030",
         "Southern hemisphere spectroscopy; MW satellite kinematics"),
        ("ELT/MOSAIC", "2028+",
         "30m-class spectroscopy; sigma < 0.5 km/s for ultra-faints"),
        ("Roman Space Tel.", "2027+",
         "Deep imaging for new ultra-diffuse galaxies around nearby groups"),
        ("JWST (ongoing)", "2023-2033+",
         "Resolved stellar kinematics in Local Group dwarfs"),
    ]

    for name, timeline, desc in surveys:
        print(f"    {name:<20s} ({timeline})")
        print(f"      {desc}")

    print(f"\n  Priority targets for EFE discrimination (ranked by feasibility):\n")

    priority = []
    for i, r in enumerate(all_results):
        split = r['split']
        if split > 0:
            err_needed = split / 3.0
            improvement = r['sig_err'] / err_needed if err_needed > 0 else np.inf
        else:
            err_needed = np.inf
            improvement = np.inf
        priority.append((r['name'], split, err_needed, improvement, r['sig_err']))

    # Sort by improvement factor needed (lower = more feasible)
    priority.sort(key=lambda x: x[3])

    for rank, (name, split, err_need, improv, curr) in enumerate(priority, 1):
        status = "ALREADY 3-sig" if improv <= 1.0 else f"need {improv:.1f}x improvement"
        print(f"    {rank}. {name:<18s}  split={split:.2f} km/s, "
              f"need err<{err_need:.2f}, current={curr:.1f}  => {status}")

    # ===========================================================
    # Final Summary
    # ===========================================================
    print(f"\n{'='*80}")
    print("  FINAL SUMMARY")
    print(f"{'='*80}\n")

    print(f"  Systems analyzed:      {len(systems)}")
    print(f"  chi2_Newton (total):   {chi2_newt:.2f}  (reduced: {chi2_newt/n_used:.2f})")
    print(f"  chi2_MOND   (total):   {chi2_mond:.2f}  (reduced: {chi2_mond/n_used:.2f})")
    print(f"  chi2_TGP    (total):   {chi2_tgp:.2f}  (reduced: {chi2_tgp/n_used:.2f})")
    print(f"  Delta_chi2 (TGP-MOND): {dchi2:+.2f}")
    print(f"  Statistical winner:    {winner}")
    print(f"  Per-system wins:       TGP={tgp_wins}, MOND={mond_wins}, TIE={ties}")

    best_disc = ranked[0]
    print(f"\n  Best discriminator:    {best_disc['name']}")
    print(f"    y_ext/y_int = {best_disc['ratio']:.3f}")
    print(f"    TGP-MOND split = {best_disc['split']:.2f} km/s")
    print(f"    Discriminating power = {best_disc['disc_power']:.3f}")

    print(f"""
  CONCLUSIONS:

  1. The TGP quadrature EFE prescription consistently predicts HIGHER velocity
     dispersions than MOND's linear prescription, because the external field
     suppresses the MOND-like boost LESS in TGP.

  2. The maximum discriminating power occurs when y_ext/y_int ~ 1, where
     quadrature vs linear addition differ most (sqrt(2) vs 2 at equal arguments).

  3. Systems close to their hosts (high y_ext) with low baryonic mass (low y_int)
     are ideal: Tucana III, Segue 1 are in the strong-EFE regime.

  4. Ultra-diffuse systems (Crater II, Antlia 2, And XIX) have low y_int due
     to their large r_half, making them sensitive to the EFE prescription.

  5. The NGC 1052 satellites (DF2, DF4) are critical: their very low observed
     sigma favors MOND's stronger EFE suppression, UNLESS systematic effects
     (distance, tidal stripping) reduce the baryonic mass estimate.

  6. Next-generation spectroscopy (ELT/MOSAIC, DESI) will achieve the precision
     needed to discriminate at 3-sigma for the best target systems.
""")

    print(sep)
    print("  gs47_efe_extended.py -- analysis complete")
    print(sep)


if __name__ == '__main__':
    main()
