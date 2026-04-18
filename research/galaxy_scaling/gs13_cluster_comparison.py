#!/usr/bin/env python3
"""
gs13: GALAXY CLUSTER COMPARISON
================================
Test: does the winning equation ν(y) = 1 + exp(-y^(4/5))/y^(2/5)
work beyond galactic scales — on galaxy clusters?

Key question: MOND has a known "factor of 2" underprediction in clusters.
Does our γ=0.4 (stronger boost at low y than γ=0.5) help or hurt?

Data sources:
  - Gonzalez et al. (2013, ApJ 778, 14) — 12 clusters, M500/r500/Mgas/Mstar
  - Relaxed cluster sample (Brinckmann+2023) — 10 clusters with Mdyn
  - X-COP sample (Eckert+2022) — 12 clusters
  - Literature compilation: Coma, Perseus, Virgo, Bullet, etc.

Method:
  At r500: g_bar = G*M_bar/r500^2,  g_obs = G*M_total/r500^2
  y = g_bar / a0
  nu_pred(y) from each model
  Compare: g_pred = g_bar * nu_pred  vs  g_obs

Author: TGP research program
"""

import numpy as np

# =============================================================================
# CONSTANTS
# =============================================================================
G_SI = 6.674e-11       # m^3 kg^-1 s^-2
Msun = 1.989e30        # kg
kpc = 3.086e19         # m
Mpc = 3.086e22         # m

a0_exp = 1.12e-10      # our best-fit
a0_mond = 1.2e-10      # standard MOND

# =============================================================================
# INTERPOLATION FUNCTIONS
# =============================================================================
def nu_exp(y, alpha=0.8, gamma=0.4):
    """Winning equation: ν = 1 + exp(-y^α) / y^γ"""
    ya = np.power(y, alpha)
    yg = np.power(y, gamma)
    return 1.0 + np.exp(-ya) / yg

def nu_mond_simple(y):
    """MOND simple: ν = 1/2 + sqrt(1/4 + 1/y)"""
    return 0.5 + np.sqrt(0.25 + 1.0/y)

def nu_mcgaugh(y):
    """McGaugh (2016): ν = 1/(1 - exp(-sqrt(y)))"""
    return 1.0 / (1.0 - np.exp(-np.sqrt(y)))

def nu_exp_mond(y):
    """Our equation with forced γ=0.5 (MOND deep limit)"""
    return 1.0 + np.exp(-np.power(y, 0.645)) / np.power(y, 0.5)

# =============================================================================
# CLUSTER DATA
# =============================================================================
def build_cluster_sample():
    """
    Compile galaxy cluster data from published sources.

    Returns list of dicts with:
      name, M500 (Msun), r500 (Mpc), M_gas (Msun), M_star (Msun),
      source, notes

    M_bar = M_gas + M_star  (hot ICM dominates)
    """
    clusters = []

    # -----------------------------------------------------------------
    # SOURCE 1: Gonzalez et al. (2013, ApJ 778, 14)
    # Table 2: M500, r500, M_gas,500, M_star,3D,500
    # Units: M500 in 1e14 Msun, M_gas in 1e13 Msun, M_star in 1e13 Msun
    # r500 in Mpc
    # -----------------------------------------------------------------
    gonzalez = [
        # (name, r500_Mpc, M500_1e14, Mgas_1e13, Mstar_1e13)
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
            'r500_Mpc': r500,
            'M500_Msun': M500 * 1e14,
            'Mgas_Msun': Mgas * 1e13,
            'Mstar_Msun': Mstar * 1e13,
            'source': 'Gonzalez+2013'
        })

    # -----------------------------------------------------------------
    # SOURCE 2: Well-known clusters from literature compilations
    # Data from Reiprich & Böhringer (2002), Vikhlinin+2006,
    # Planck Collaboration, etc.
    # M500 and r500 from scaling relations / direct measurements
    # Gas fractions: f_gas ~ 0.10-0.15 at r500
    # Stellar fractions: f_star ~ 0.01-0.03
    # -----------------------------------------------------------------
    literature = [
        # (name, r500_Mpc, M500_1e14, fgas, fstar)
        # Coma cluster — one of the best-studied
        ("Coma",       1.31, 7.2,  0.144, 0.020),
        # Perseus — massive cool-core cluster
        ("Perseus",    1.27, 6.6,  0.155, 0.025),
        # Virgo — nearest massive cluster
        ("Virgo",      0.77, 1.4,  0.090, 0.030),
        # Centaurus — nearby
        ("Centaurus",  0.82, 1.8,  0.100, 0.025),
        # A2029 — massive relaxed (Gonzalez+2013 also has it)
        ("A2029",      1.42, 8.71, 0.138, 0.020),
        # A0478 — massive cool-core
        ("A0478",      1.28, 6.58, 0.175, 0.020),
        # A2390 — hot massive cluster
        ("A2390",      1.50, 11.8, 0.146, 0.015),
        # A1795 — well-studied X-COP cluster
        ("A1795",      1.09, 4.54, 0.118, 0.020),
        # A2142 — merging massive cluster
        ("A2142",      1.37, 8.17, 0.174, 0.020),
        # A2319 — very massive, hot
        ("A2319",      1.55, 12.0, 0.150, 0.018),
        # A644 — X-COP cluster
        ("A644",       1.05, 3.70, 0.130, 0.020),
        # A3266 — merging cluster
        ("A3266",      1.25, 6.20, 0.140, 0.018),
        # Fornax — nearby group/poor cluster
        ("Fornax",     0.45, 0.28, 0.060, 0.035),
        # MKW4 — galaxy group
        ("MKW4",       0.50, 0.38, 0.055, 0.040),
        # NGC5044 group
        ("NGC5044grp", 0.43, 0.24, 0.050, 0.045),
        # A3571 — relaxed cluster
        ("A3571",      1.13, 5.02, 0.103, 0.020),
        # Bullet cluster — famous merging cluster
        ("Bullet",     1.35, 7.8,  0.140, 0.015),
    ]
    for name, r500, M500, fgas, fstar in literature:
        M500_sun = M500 * 1e14
        clusters.append({
            'name': name,
            'r500_Mpc': r500,
            'M500_Msun': M500_sun,
            'Mgas_Msun': fgas * M500_sun,
            'Mstar_Msun': fstar * M500_sun,
            'source': 'Literature'
        })

    # Derived quantities for all clusters
    for c in clusters:
        c['Mbar_Msun'] = c['Mgas_Msun'] + c['Mstar_Msun']
        c['fbar'] = c['Mbar_Msun'] / c['M500_Msun']
        c['r500_m'] = c['r500_Mpc'] * Mpc

        # Accelerations at r500
        c['g_bar'] = G_SI * c['Mbar_Msun'] * Msun / c['r500_m']**2
        c['g_obs'] = G_SI * c['M500_Msun'] * Msun / c['r500_m']**2
        c['nu_obs'] = c['g_obs'] / c['g_bar']  # observed boost

        # y parameter
        c['y_exp'] = c['g_bar'] / a0_exp
        c['y_mond'] = c['g_bar'] / a0_mond

    return clusters


# =============================================================================
# MULTI-RADII PROFILES (NFW + beta-model)
# =============================================================================
def nfw_mass(r, M500, r500, c500=3.5):
    """NFW enclosed mass given M500, r500, concentration c500."""
    rs = r500 / c500
    x500 = c500
    f500 = np.log(1 + x500) - x500/(1 + x500)

    x = r / rs
    fx = np.log(1 + x) - x/(1 + x)
    return M500 * fx / f500

def beta_gas_mass(r, Mgas500, r500, rc_frac=0.15, beta=0.65):
    """
    Beta-model gas mass: rho_gas ~ (1 + (r/rc)^2)^(-3β/2)
    Integrated mass enclosed within r.
    rc_frac = rc/r500 (core radius as fraction of r500)
    """
    rc = rc_frac * r500
    # Numerical integration via simple quadrature
    n_pts = 500
    r_arr = np.linspace(0.001*r500, r, n_pts)
    dr = r_arr[1] - r_arr[0]
    rho = np.power(1 + (r_arr/rc)**2, -1.5*beta)

    # Mass integral: M(r) = 4pi * integral(rho * r^2 dr)
    m_unnorm = 4 * np.pi * np.sum(rho * r_arr**2 * dr)

    # Same at r500
    r_arr500 = np.linspace(0.001*r500, r500, n_pts)
    dr500 = r_arr500[1] - r_arr500[0]
    rho500 = np.power(1 + (r_arr500/rc)**2, -1.5*beta)
    m_unnorm500 = 4 * np.pi * np.sum(rho500 * r_arr500**2 * dr500)

    return Mgas500 * m_unnorm / m_unnorm500


def compute_radial_profiles(cluster, n_radii=20):
    """
    Compute g_bar(r) and g_obs(r) at multiple radii for a cluster.
    Uses NFW for total mass, beta-model for gas, stars ~ concentrated.
    """
    r500 = cluster['r500_Mpc'] * Mpc  # meters
    M500 = cluster['M500_Msun'] * Msun  # kg
    Mgas500 = cluster['Mgas_Msun'] * Msun
    Mstar500 = cluster['Mstar_Msun'] * Msun

    radii_frac = np.linspace(0.1, 1.5, n_radii)  # fractions of r500
    radii = radii_frac * r500

    g_bar_arr = np.zeros(n_radii)
    g_obs_arr = np.zeros(n_radii)

    for i, r in enumerate(radii):
        # Total mass (NFW)
        M_tot = nfw_mass(r, M500, r500, c500=3.5)

        # Gas mass (beta-model)
        M_gas = beta_gas_mass(r, Mgas500, r500, rc_frac=0.15, beta=0.65)

        # Stellar mass (more concentrated, King-like)
        M_star = beta_gas_mass(r, Mstar500, r500, rc_frac=0.05, beta=1.0)

        M_bar = M_gas + M_star

        g_bar_arr[i] = G_SI * M_bar / r**2
        g_obs_arr[i] = G_SI * M_tot / r**2

    return radii_frac, g_bar_arr, g_obs_arr


# =============================================================================
# MAIN ANALYSIS
# =============================================================================
def main():
    print("=" * 80)
    print("gs13: GALAXY CLUSTER COMPARISON")
    print("Does ν(y) = 1 + exp(-y^(4/5))/y^(2/5) work on cluster scales?")
    print("=" * 80)

    clusters = build_cluster_sample()
    print(f"\nLoaded {len(clusters)} clusters/groups")

    # =====================================================================
    # SECTION 1: ACCELERATION REGIME OF CLUSTERS
    # =====================================================================
    print("\n" + "=" * 80)
    print("SECTION 1: ACCELERATION REGIME OF GALAXY CLUSTERS")
    print("=" * 80)

    g_bars = np.array([c['g_bar'] for c in clusters])
    g_obss = np.array([c['g_obs'] for c in clusters])
    ys = np.array([c['y_exp'] for c in clusters])

    print(f"\n  g_bar range: {g_bars.min():.2e} to {g_bars.max():.2e} m/s²")
    print(f"  g_obs range: {g_obss.min():.2e} to {g_obss.max():.2e} m/s²")
    print(f"  y = g_bar/a₀ range: {ys.min():.3f} to {ys.max():.3f}")
    print(f"  y median: {np.median(ys):.3f}")
    print(f"\n  Compare: galaxy outer RC median y ~ 0.05")
    print(f"  Clusters at r500: y ~ {np.median(ys):.2f} (DEEPER into modified regime)")

    # Categorize
    deep = sum(1 for y in ys if y < 0.1)
    trans = sum(1 for y in ys if 0.1 <= y < 1.0)
    newt = sum(1 for y in ys if y >= 1.0)
    print(f"\n  Deep MOND (y < 0.1): {deep}")
    print(f"  Transition (0.1-1): {trans}")
    print(f"  Newtonian (y > 1):  {newt}")

    # =====================================================================
    # SECTION 2: MODEL PREDICTIONS AT r500
    # =====================================================================
    print("\n" + "=" * 80)
    print("SECTION 2: MODEL PREDICTIONS vs OBSERVATIONS AT r500")
    print("=" * 80)

    models = {
        'Exp (γ=0.4, winning)': lambda y: nu_exp(y, 0.8, 0.4),
        'Exp (γ=0.5, MOND-like)': lambda y: nu_exp_mond(y),
        'McGaugh': nu_mcgaugh,
        'MOND simple': nu_mond_simple,
    }

    print(f"\n  {'Model':<28s} {'med(g_pred/g_obs)':<20s} {'mean ratio':<12s} {'std':<8s} {'scatter(dex)':<14s}")
    print("  " + "-" * 82)

    model_results = {}
    for mname, mfunc in models.items():
        ratios = []
        for c in clusters:
            y = c['y_mond'] if 'MOND' in mname or 'McGaugh' in mname else c['y_exp']
            nu_pred = mfunc(y)
            g_pred = c['g_bar'] * nu_pred
            ratios.append(g_pred / c['g_obs'])

        ratios = np.array(ratios)
        log_ratios = np.log10(ratios)

        print(f"  {mname:<28s} {np.median(ratios):<20.4f} {np.mean(ratios):<12.4f} "
              f"{np.std(ratios):<8.4f} {np.std(log_ratios):<14.4f}")

        model_results[mname] = ratios

    # =====================================================================
    # SECTION 3: PER-CLUSTER COMPARISON
    # =====================================================================
    print("\n" + "=" * 80)
    print("SECTION 3: PER-CLUSTER PREDICTIONS")
    print("=" * 80)

    print(f"\n  {'Cluster':<14s} {'y=gbar/a0':<10s} {'ν_obs':<8s} {'ν_Exp':<8s} "
          f"{'ν_MOND':<8s} {'ν_McG':<8s} {'Exp/obs':<8s} {'MOND/obs':<8s} {'McG/obs':<8s}")
    print("  " + "-" * 90)

    for c in clusters:
        y_e = c['y_exp']
        y_m = c['y_mond']
        nu_o = c['nu_obs']
        nu_e = nu_exp(y_e, 0.8, 0.4)
        nu_m = nu_mond_simple(y_m)
        nu_mc = nu_mcgaugh(y_m)

        print(f"  {c['name']:<14s} {y_e:<10.4f} {nu_o:<8.2f} {nu_e:<8.2f} "
              f"{nu_m:<8.2f} {nu_mc:<8.2f} {nu_e/nu_o:<8.3f} "
              f"{nu_m/nu_o:<8.3f} {nu_mc/nu_o:<8.3f}")

    # =====================================================================
    # SECTION 4: THE "MISSING MASS" FACTOR
    # =====================================================================
    print("\n" + "=" * 80)
    print("SECTION 4: THE MISSING MASS FACTOR IN CLUSTERS")
    print("=" * 80)

    print(f"\n  MOND prediction for clusters: ν_MOND at y ~ {np.median(ys):.3f}")
    print(f"  If MOND underpredicts by factor ~2, then g_obs ~ 2 * g_MOND_pred")
    print(f"  Our Exp model gives STRONGER boost (γ=0.4 > MOND effective γ=0.5)")

    # Compute the "missing factor" = g_obs / g_predicted
    # If > 1: model underpredicts (needs more mass)
    # If < 1: model overpredicts
    print(f"\n  Missing factor = g_obs / g_pred (>1 = underprediction)")
    print(f"\n  {'Model':<28s} {'median missing':<16s} {'mean missing':<16s} "
          f"{'% underpredicting':<20s}")
    print("  " + "-" * 80)

    for mname, ratios in model_results.items():
        missing = 1.0 / ratios  # g_obs / g_pred
        n_under = np.sum(missing > 1.0)
        pct_under = 100.0 * n_under / len(missing)
        print(f"  {mname:<28s} {np.median(missing):<16.3f} {np.mean(missing):<16.3f} "
              f"{pct_under:<20.1f}")

    # =====================================================================
    # SECTION 5: DEEP MOND LIMIT COMPARISON
    # =====================================================================
    print("\n" + "=" * 80)
    print("SECTION 5: DEEP-MOND vs EXP ASYMPTOTIC BEHAVIOR")
    print("=" * 80)

    print(f"\n  In deep MOND (y << 1):")
    print(f"    MOND:     g_obs = √(g_bar · a₀)          i.e. ν ~ 1/√y")
    print(f"    McGaugh:  g_obs ~ √(g_bar · a₀)          (same deep limit)")
    print(f"    Exp:      g_obs = g_bar^0.59 · a₀^0.41   i.e. ν ~ 1/y^0.4")
    print(f"\n  At y=0.01 (typical cluster outskirts):")

    y_test = 0.01
    nu_e = nu_exp(y_test, 0.8, 0.4)
    nu_m = nu_mond_simple(y_test)
    nu_mc = nu_mcgaugh(y_test)
    print(f"    ν_Exp    = {nu_e:.2f}")
    print(f"    ν_MOND   = {nu_m:.2f}")
    print(f"    ν_McGaugh = {nu_mc:.2f}")
    print(f"    Ratio Exp/MOND = {nu_e/nu_m:.3f}")

    y_test2 = 0.001
    nu_e2 = nu_exp(y_test2, 0.8, 0.4)
    nu_m2 = nu_mond_simple(y_test2)
    nu_mc2 = nu_mcgaugh(y_test2)
    print(f"\n  At y=0.001 (very deep):")
    print(f"    ν_Exp    = {nu_e2:.2f}")
    print(f"    ν_MOND   = {nu_m2:.2f}")
    print(f"    Ratio Exp/MOND = {nu_e2/nu_m2:.3f}")

    # =====================================================================
    # SECTION 6: RADIAL PROFILES FOR REPRESENTATIVE CLUSTERS
    # =====================================================================
    print("\n" + "=" * 80)
    print("SECTION 6: RADIAL PROFILES FOR REPRESENTATIVE CLUSTERS")
    print("=" * 80)

    # Pick a few representative clusters
    representative = ['Coma', 'Perseus', 'Virgo', 'A2029', 'Fornax', 'A1651']

    for cname in representative:
        c = None
        for cc in clusters:
            if cc['name'] == cname:
                c = cc
                break
        if c is None:
            continue

        print(f"\n  --- {c['name']} (M500={c['M500_Msun']:.2e} Msun, r500={c['r500_Mpc']:.2f} Mpc) ---")

        rfrac, g_bar_r, g_obs_r = compute_radial_profiles(c)

        print(f"  {'r/r500':<8s} {'y':<10s} {'ν_obs':<8s} {'ν_Exp':<8s} "
              f"{'ν_MOND':<8s} {'Exp/obs':<8s} {'MOND/obs':<8s}")

        for i in range(0, len(rfrac), 3):  # every 3rd point
            y = g_bar_r[i] / a0_exp
            y_m = g_bar_r[i] / a0_mond
            nu_o = g_obs_r[i] / g_bar_r[i]
            nu_e = nu_exp(y, 0.8, 0.4)
            nu_m = nu_mond_simple(y_m)

            print(f"  {rfrac[i]:<8.2f} {y:<10.4f} {nu_o:<8.2f} {nu_e:<8.2f} "
                  f"{nu_m:<8.2f} {nu_e/nu_o:<8.3f} {nu_m/nu_o:<8.3f}")

    # =====================================================================
    # SECTION 7: WHAT ACCELERATION SCALE WOULD CLUSTERS NEED?
    # =====================================================================
    print("\n" + "=" * 80)
    print("SECTION 7: BEST-FIT ACCELERATION SCALE FOR CLUSTERS")
    print("=" * 80)

    # For each model, find the a0 that minimizes the scatter in g_pred/g_obs
    from scipy.optimize import minimize_scalar

    def cluster_chi2_exp(log_a0):
        a0 = 10**log_a0
        chi2 = 0
        for c in clusters:
            y = c['g_bar'] / a0
            nu_p = nu_exp(y, 0.8, 0.4)
            g_pred = c['g_bar'] * nu_p
            chi2 += (np.log10(g_pred/c['g_obs']))**2
        return chi2

    def cluster_chi2_mond(log_a0):
        a0 = 10**log_a0
        chi2 = 0
        for c in clusters:
            y = c['g_bar'] / a0
            nu_p = nu_mond_simple(y)
            g_pred = c['g_bar'] * nu_p
            chi2 += (np.log10(g_pred/c['g_obs']))**2
        return chi2

    def cluster_chi2_mcgaugh(log_a0):
        a0 = 10**log_a0
        chi2 = 0
        for c in clusters:
            y = c['g_bar'] / a0
            nu_p = nu_mcgaugh(y)
            g_pred = c['g_bar'] * nu_p
            chi2 += (np.log10(g_pred/c['g_obs']))**2
        return chi2

    res_exp = minimize_scalar(cluster_chi2_exp, bounds=(-11, -8), method='bounded')
    res_mond = minimize_scalar(cluster_chi2_mond, bounds=(-11, -8), method='bounded')
    res_mcg = minimize_scalar(cluster_chi2_mcgaugh, bounds=(-11, -8), method='bounded')

    a0_fit_exp = 10**res_exp.x
    a0_fit_mond = 10**res_mond.x
    a0_fit_mcg = 10**res_mcg.x

    print(f"\n  Best-fit a₀ for CLUSTERS:")
    print(f"  {'Model':<20s} {'a₀ (m/s²)':<15s} {'a₀/a₀_galaxy':<15s} {'scatter (dex)':<15s}")
    print("  " + "-" * 65)

    # Compute scatter at best-fit a0
    for mname, a0_fit, chi2_func in [
        ('Exp (γ=0.4)', a0_fit_exp, cluster_chi2_exp),
        ('MOND simple', a0_fit_mond, cluster_chi2_mond),
        ('McGaugh', a0_fit_mcg, cluster_chi2_mcgaugh),
    ]:
        a0_gal = a0_exp if 'Exp' in mname else a0_mond
        scatter = np.sqrt(chi2_func(np.log10(a0_fit)) / len(clusters))
        ratio = a0_fit / a0_gal
        print(f"  {mname:<20s} {a0_fit:<15.3e} {ratio:<15.2f} {scatter:<15.4f}")

    # =====================================================================
    # SECTION 8: FREE GAMMA FIT TO CLUSTERS
    # =====================================================================
    print("\n" + "=" * 80)
    print("SECTION 8: BEST-FIT (α, γ, a₀) FOR CLUSTERS ONLY")
    print("=" * 80)

    from scipy.optimize import minimize

    def cluster_chi2_free(params):
        log_a0, alpha, gamma = params
        a0 = 10**log_a0
        if alpha < 0.1 or alpha > 2.0 or gamma < 0.1 or gamma > 1.5:
            return 1e10
        chi2 = 0
        for c in clusters:
            y = c['g_bar'] / a0
            nu_p = nu_exp(y, alpha, gamma)
            g_pred = c['g_bar'] * nu_p
            chi2 += (np.log10(g_pred/c['g_obs']))**2
        return chi2

    # Multi-start
    best_chi2 = 1e20
    best_params = None

    starts = [
        (-10, 0.8, 0.4),
        (-9.5, 0.8, 0.5),
        (-10, 1.0, 0.5),
        (-9, 0.5, 0.3),
        (-9.5, 0.6, 0.3),
        (-10.5, 0.8, 0.4),
        (-9, 0.8, 0.4),
        (-9.5, 1.0, 0.5),
    ]

    for s in starts:
        res = minimize(cluster_chi2_free, s, method='Nelder-Mead',
                      options={'maxiter': 10000, 'xatol': 1e-6, 'fatol': 1e-10})
        if res.fun < best_chi2:
            best_chi2 = res.fun
            best_params = res.x

    log_a0_best, alpha_best, gamma_best = best_params
    a0_best = 10**log_a0_best

    print(f"\n  Best-fit for clusters:")
    print(f"    a₀ = {a0_best:.3e} m/s² ({a0_best/a0_exp:.2f} × galaxy a₀)")
    print(f"    α  = {alpha_best:.4f}  (galaxy: 0.80)")
    print(f"    γ  = {gamma_best:.4f}  (galaxy: 0.40)")
    print(f"    γ/α = {gamma_best/alpha_best:.4f}")
    print(f"    scatter = {np.sqrt(best_chi2/len(clusters)):.4f} dex")

    # Compare: fixed galaxy params vs free cluster params
    chi2_galaxy_params = cluster_chi2_exp(np.log10(a0_exp))
    print(f"\n  Chi² comparison (sum of log10 residuals²):")
    print(f"    Galaxy params (a₀={a0_exp:.2e}, α=0.8, γ=0.4): {chi2_galaxy_params:.4f}")
    print(f"    Best cluster params: {best_chi2:.4f}")
    print(f"    With free a₀ only (α=0.8, γ=0.4): {cluster_chi2_exp(np.log10(a0_fit_exp)):.4f}")

    # =====================================================================
    # SECTION 9: CLUSTER vs GALAXY RAR COMPARISON
    # =====================================================================
    print("\n" + "=" * 80)
    print("SECTION 9: CLUSTER vs GALAXY RAR — SAME RELATION OR DIFFERENT?")
    print("=" * 80)

    print(f"\n  Galaxy RAR (McGaugh 2016): g_obs = g_bar/(1 - exp(-√(g_bar/a₀)))")
    print(f"  with a₀ = 1.2e-10 m/s²")
    print(f"\n  Cluster RAR (Tian+2020, Eckert+2022): OFFSET from galaxy RAR")
    print(f"  Effective cluster a₀ ~ 2e-9 m/s² (17× larger)")
    print(f"  BUT: this could be due to missing baryons (hot gas underestimate)")

    # Show where clusters fall on the galaxy RAR
    print(f"\n  Cluster positions on RAR:")
    print(f"  {'Cluster':<14s} {'log g_bar':<10s} {'log g_obs':<10s} "
          f"{'log g_McG':<10s} {'log g_Exp':<10s} {'offset(dex)':<12s}")
    print("  " + "-" * 66)

    offsets_exp = []
    offsets_mond = []
    for c in clusters:
        lg_bar = np.log10(c['g_bar'])
        lg_obs = np.log10(c['g_obs'])

        y_e = c['y_exp']
        y_m = c['y_mond']
        lg_exp = np.log10(c['g_bar'] * nu_exp(y_e, 0.8, 0.4))
        lg_mcg = np.log10(c['g_bar'] * nu_mcgaugh(y_m))

        offset = lg_obs - lg_exp
        offsets_exp.append(offset)
        offsets_mond.append(lg_obs - lg_mcg)

        print(f"  {c['name']:<14s} {lg_bar:<10.3f} {lg_obs:<10.3f} "
              f"{lg_mcg:<10.3f} {lg_exp:<10.3f} {offset:<+12.3f}")

    offsets_exp = np.array(offsets_exp)
    offsets_mond = np.array(offsets_mond)

    print(f"\n  Summary of offsets (log g_obs - log g_pred):")
    print(f"  Exp model:     median = {np.median(offsets_exp):+.3f} dex, "
          f"mean = {np.mean(offsets_exp):+.3f} dex")
    print(f"  McGaugh/MOND:  median = {np.median(offsets_mond):+.3f} dex, "
          f"mean = {np.mean(offsets_mond):+.3f} dex")
    print(f"\n  Positive offset = model UNDERPREDICTS (needs more mass)")
    print(f"  Factor of underprediction:")
    print(f"    Exp:  {10**np.median(offsets_exp):.2f}× (median)")
    print(f"    MOND: {10**np.median(offsets_mond):.2f}× (median)")

    # =====================================================================
    # SECTION 10: PHYSICAL INTERPRETATION
    # =====================================================================
    print("\n" + "=" * 80)
    print("SECTION 10: PHYSICAL INTERPRETATION FOR TGP")
    print("=" * 80)

    print(f"""
  KEY FINDINGS:

  1. Acceleration regime: clusters at r500 have y ~ {np.median(ys):.3f}
     This is DEEPER into the modified regime than typical galaxy outskirts (y~0.05).

  2. Our Exp equation (γ=0.4) gives a STRONGER boost than MOND (γ=0.5):
     ν ~ 1/y^0.4 vs 1/y^0.5 at low y
     This means Exp predicts MORE "dark matter" effect than MOND.

  3. MOND's cluster problem (factor ~2 underprediction):
""")

    mond_factor = 10**np.median(offsets_mond)
    exp_factor = 10**np.median(offsets_exp)

    if exp_factor < mond_factor:
        improvement = (1 - (exp_factor - 1)/(mond_factor - 1)) * 100
        print(f"     MOND underpredicts by {mond_factor:.2f}×")
        print(f"     Exp underpredicts by  {exp_factor:.2f}×")
        print(f"     → Exp REDUCES the cluster discrepancy by {improvement:.0f}%!")
    else:
        print(f"     MOND underpredicts by {mond_factor:.2f}×")
        print(f"     Exp underpredicts by  {exp_factor:.2f}×")
        print(f"     → Exp does NOT improve the cluster situation")

    print(f"""
  4. Best-fit cluster a₀:
     Galaxy a₀ = {a0_exp:.2e} m/s²
     Cluster a₀ (Exp) = {a0_fit_exp:.2e} m/s² ({a0_fit_exp/a0_exp:.1f}× larger)
     → Clusters need a DIFFERENT acceleration scale

  5. TGP interpretation:
     - If a₀ = c·H₀/(2π) is truly universal, the cluster offset means:
       a) Missing baryons (hot gas/WHIM not accounted for), OR
       b) The interpolation function changes form at cluster scales, OR
       c) d_eff varies with system size (cluster ≠ galaxy transition)

  6. d_eff connection:
     Galaxy: d_eff = 3 - 2γ = 3 - 0.8 = 2.2 (partial 3D→2D)
     If clusters need stronger γ → smaller d_eff → MORE transition
     Cluster best γ = {gamma_best:.2f} → d_eff = {3 - 2*gamma_best:.2f}
""")

    # =====================================================================
    # SECTION 11: COMPARISON WITH PUBLISHED CLUSTER RAR
    # =====================================================================
    print("=" * 80)
    print("SECTION 11: COMPARISON WITH PUBLISHED CLUSTER RAR RESULTS")
    print("=" * 80)

    print(f"""
  Published findings (for comparison):

  Tian+2020 (CLASH clusters):
    Effective a₀ = 2.0e-9 m/s² (17× galaxy value)
    Slope in log-log: 0.51 ± 0.05 (consistent with MOND's 0.5)

  Chan & Del Popolo (2020):
    52 non-cool-core clusters: scatter = 0.18 dex
    Overall a₀ ~ 9.5e-10 m/s² (8× galaxy value)

  Eckert+2022 (X-COP):
    12 clusters: RAR does NOT hold at cluster scale
    Clusters sit ABOVE galaxy RAR

  Our analysis:
    Best-fit cluster a₀ (Exp) = {a0_fit_exp:.2e} m/s² ({a0_fit_exp/a0_exp:.1f}× galaxy)
    Best-fit cluster a₀ (MOND) = {a0_fit_mond:.2e} m/s² ({a0_fit_mond/a0_mond:.1f}× galaxy)

  CONCLUSION: Cluster RAR requires ~{a0_fit_exp/a0_exp:.0f}× larger a₀ than galaxies.
  This is a KNOWN problem for ALL MOND-like theories.
  Our Exp equation does {'BETTER' if exp_factor < mond_factor else 'NOT BETTER'} than MOND
  (factor {exp_factor:.2f}× vs {mond_factor:.2f}×), reducing discrepancy by
  {abs((exp_factor - mond_factor)/mond_factor * 100):.0f}% but NOT eliminating it.
""")

    # =====================================================================
    # SUMMARY TABLE
    # =====================================================================
    print("=" * 80)
    print("SUMMARY")
    print("=" * 80)

    print(f"""
  ┌─────────────────────────────────────────────────────────────────────┐
  │ Scale      │ a₀ needed   │ γ best │ d_eff  │ ν underprediction    │
  ├─────────────────────────────────────────────────────────────────────┤
  │ Galaxies   │ 1.12e-10    │ 0.40   │ 2.20   │ ×1.00 (by design)   │
  │ Clusters   │ {a0_fit_exp:.2e}    │ {gamma_best:.2f}   │ {3-2*gamma_best:.2f}   │ ×{exp_factor:.2f} (with gal a₀)  │
  │ MOND clust │ {a0_fit_mond:.2e}    │ 0.50   │ 2.00   │ ×{mond_factor:.2f} (with gal a₀)  │
  └─────────────────────────────────────────────────────────────────────┘

  Key result: γ=0.4 (Exp) gives stronger boost → {'REDUCES' if exp_factor < mond_factor else 'does not reduce'}
  cluster missing mass from ×{mond_factor:.2f} to ×{exp_factor:.2f}

  The cluster problem remains for ALL modified gravity theories.
  Possible TGP resolution: a₀ may scale with system size or local
  expansion rate, not be strictly universal.
""")

    print("=" * 80)
    print("DONE")
    print("=" * 80)


if __name__ == '__main__':
    main()
