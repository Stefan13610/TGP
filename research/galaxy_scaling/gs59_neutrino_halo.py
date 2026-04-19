#!/usr/bin/env python3
"""
gs59_neutrino_halo.py - Full Poisson-Boltzmann neutrino halo calculation for TGP

Extends gs55 simple Tremaine-Gunn estimates with self-consistent
Poisson-Boltzmann iteration for sterile neutrino halos in galaxy clusters.

Physics: Fermi-Dirac neutrinos in TGP-enhanced gravitational potential.
References: Sanders (2003), Angus (2009)

The key physics:
  - Neutrinos are virialized in the cluster potential (kT_eff = m * sigma^2)
  - Their density follows Boltzmann: rho(r) = rho_boundary * exp(-[Phi(r)-Phi_boundary]/sigma^2)
  - Capped by Tremaine-Gunn phase-space limit
  - Self-consistent iteration includes neutrino self-gravity
"""

import numpy as np
from scipy.integrate import quad, cumulative_trapezoid

# ============================================================
# CONSTANTS
# ============================================================
G       = 6.674e-11       # m^3/(kg*s^2)
Msun    = 1.989e30        # kg
kpc     = 3.0857e19       # m
hbar    = 1.0546e-34      # J*s
k_B     = 1.381e-23       # J/K
c_light = 2.998e8         # m/s
eV_kg   = 1.783e-36       # kg per eV/c^2
eV_J    = 1.602e-19       # J per eV
T_nu    = 1.95            # K, cosmic neutrino temperature
n_cosmic = 56e6           # /m^3 per flavor (56/cm^3)
g_s     = 2               # spin states
a0      = 1.2e-10         # m/s^2, TGP/MOND acceleration scale

# Neutrino masses to test (eV)
m_s_list = [1.5, 2.0, 5.0, 11.0, 20.0]

# Radii for profile output (kpc)
r_profile_kpc = np.array([10, 50, 100, 200, 500, 1000, 1500])

# ============================================================
# TGP nu function
# ============================================================
def nu_tgp(y):
    """TGP interpolating function nu(y) where y = g_N / a0."""
    y = np.asarray(y, dtype=float)
    result = np.ones_like(y)
    mask = y > 0
    result[mask] = 1.0 + np.exp(-y[mask]**0.8) / y[mask]**0.5714
    return result

# ============================================================
# PART A: Cluster model - beta profile + TGP potential
# ============================================================

def beta_rho0(M_bar, R500, r_c, beta):
    """Compute central density rho_0 for beta model so that M_bar(R500) = M_bar."""
    def integrand(r):
        return 4.0 * np.pi * r**2 / (1.0 + (r / r_c)**2)**(3.0 * beta / 2.0)
    shape_integral, _ = quad(integrand, 0, R500)
    return M_bar / shape_integral

def make_radial_grid(R500_m, n=500):
    """Create logarithmic radial grid from 1 kpc to 3*R500."""
    r_min = 1.0 * kpc
    r_max = 3.0 * R500_m
    return np.logspace(np.log10(r_min), np.log10(r_max), n)

def compute_baryonic_mass_profile(r_grid, rho0, r_c, beta):
    """Compute enclosed baryonic mass M_bar(r) on grid."""
    M_bar = np.zeros_like(r_grid)
    for i, R in enumerate(r_grid):
        def integrand(r):
            return 4.0 * np.pi * r**2 * rho0 / (1.0 + (r / r_c)**2)**(3.0 * beta / 2.0)
        M_bar[i], _ = quad(integrand, 0, R, limit=100)
    return M_bar

def compute_tgp_potential(r_grid, M_enc):
    """
    Compute TGP-enhanced gravitational potential.
    Phi_TGP(r) = -integral_r^inf g_TGP(r') dr'
    g_TGP(r) = nu(g_N/a0) * g_N(r)

    Returns potential referenced to Phi(inf)=0, so Phi < 0 everywhere.
    """
    g_N = G * M_enc / r_grid**2
    y = g_N / a0
    nu_val = nu_tgp(y)
    g_TGP = nu_val * g_N

    n = len(r_grid)
    Phi = np.zeros(n)

    # Integrate from each r to r_max
    for i in range(n - 1):
        Phi[i] = -np.trapezoid(g_TGP[i:], r_grid[i:])

    # Add tail beyond r_max: Phi_tail ~ -G*M/r_max (Newtonian at large r)
    tail = -G * M_enc[-1] / r_grid[-1]
    Phi[-1] = tail
    Phi[:-1] += tail

    return Phi

def compute_tgp_potential_from_total(r_grid, M_bar, M_extra):
    """Compute TGP potential using total enclosed mass = M_bar + M_extra."""
    M_total = M_bar + M_extra
    return compute_tgp_potential(r_grid, M_total)

# ============================================================
# PART B: Neutrino density profiles
# ============================================================

def tremaine_gunn_density(m_s_eV, sigma_ms):
    """Tremaine-Gunn phase-space density limit.
    rho_TG = g_s * m^4 * sigma^3 / (6*pi^2 * hbar^3)
    """
    m_kg = m_s_eV * eV_kg
    return g_s * m_kg**4 * sigma_ms**3 / (6.0 * np.pi**2 * hbar**3)

def boltzmann_density(r_grid, Phi, m_s_eV, sigma_ms, R500_m):
    """
    Isothermal Boltzmann neutrino density with TG cap.

    Neutrinos are virialized: kT_eff = m_s * sigma^2.
    Density referenced to boundary (R500 or grid edge):
      rho(r) = rho_boundary * exp(-(Phi(r) - Phi_boundary) / sigma^2)

    At the boundary (virial radius), rho = rho_cosmic (background).
    The potential well depth Delta_Phi = Phi(center) - Phi(boundary) < 0,
    so the exponential enhances density toward center.

    Capped by Tremaine-Gunn limit.
    """
    m_kg = m_s_eV * eV_kg
    rho_cosmic = n_cosmic * m_kg

    # Reference potential at boundary (use ~ 2*R500 as virial boundary)
    r_boundary = 2.0 * R500_m
    Phi_boundary = np.interp(r_boundary, r_grid, Phi)

    # Delta Phi relative to boundary: negative inside cluster (deeper well)
    Delta_Phi = Phi - Phi_boundary  # < 0 inside cluster

    # Exponent: -Delta_Phi / sigma^2 > 0 inside cluster
    exponent = -Delta_Phi / sigma_ms**2

    # Cap exponent to avoid overflow
    exponent = np.clip(exponent, 0, 700)
    rho_Boltz = rho_cosmic * np.exp(exponent)

    # TG cap
    rho_TG = tremaine_gunn_density(m_s_eV, sigma_ms)
    rho_nu = np.minimum(rho_Boltz, rho_TG)
    return rho_nu

def fermi_density(r_grid, Phi, m_s_eV):
    """
    Degenerate Fermi (maximum packing) density.
    rho_Fermi(r) = g_s * m^4 * v_esc(r)^3 / (6*pi^2*hbar^3)
    """
    m_kg = m_s_eV * eV_kg
    v_esc = np.sqrt(2.0 * np.abs(Phi))
    rho_Fermi = g_s * m_kg**4 * v_esc**3 / (6.0 * np.pi**2 * hbar**3)
    return rho_Fermi

def enclosed_mass_from_density(r_grid, rho):
    """Compute enclosed mass M(r) = 4*pi * integral_0^r rho*r'^2 dr'."""
    integrand = 4.0 * np.pi * r_grid**2 * rho
    M_enc = np.zeros_like(r_grid)
    M_enc[1:] = cumulative_trapezoid(integrand, r_grid)
    return M_enc

def interp_mass_at_R500(r_grid, M_enc, R500_m):
    """Interpolate enclosed mass at R500."""
    return np.interp(R500_m, r_grid, M_enc)

# ============================================================
# PART C: Self-consistent iteration
# ============================================================

def self_consistent_neutrino(r_grid, M_bar, R500_m, m_s_eV, sigma_ms,
                              max_iter=50, tol=0.01):
    """
    Iterate to self-consistency:
    1. Start with Phi = Phi_TGP (baryons only)
    2. Compute rho_nu from Phi
    3. Update total mass: M_total(r) = M_bar(r) + M_nu(r)
    4. Recompute Phi from M_total using TGP nu function
    5. Repeat until |M_nu change| < 1%
    """
    M_nu = np.zeros_like(r_grid)
    M_nu_R500_prev = 0.0

    for iteration in range(max_iter):
        Phi = compute_tgp_potential_from_total(r_grid, M_bar, M_nu)
        rho_nu = boltzmann_density(r_grid, Phi, m_s_eV, sigma_ms, R500_m)
        M_nu = enclosed_mass_from_density(r_grid, rho_nu)
        M_nu_R500 = interp_mass_at_R500(r_grid, M_nu, R500_m)

        if iteration > 0 and M_nu_R500_prev > 0:
            change = abs(M_nu_R500 - M_nu_R500_prev) / M_nu_R500_prev
            if change < tol:
                break
        M_nu_R500_prev = M_nu_R500

    return rho_nu, M_nu, M_nu_R500, iteration + 1

# ============================================================
# NFW reference profile
# ============================================================

def nfw_density(r_grid, M200, r_s):
    """NFW dark matter density profile for reference."""
    rho_crit = 9.5e-27  # kg/m^3 (approximate)
    R200 = (3.0 * M200 / (4.0 * np.pi * 200.0 * rho_crit))**(1.0/3.0)
    c200 = R200 / r_s
    f_c = np.log(1.0 + c200) - c200 / (1.0 + c200)
    rho_s = M200 / (4.0 * np.pi * r_s**3 * f_c)

    x = r_grid / r_s
    rho = rho_s / (x * (1.0 + x)**2)
    return rho

# ============================================================
# Cluster data
# ============================================================

clusters = {
    'Coma':    {'M_bar': 1.0e14*Msun, 'M_obs': 7.0e14*Msun, 'R500': 1200*kpc,
                'sigma': 1000e3, 'r_c': 250*kpc, 'beta': 2.0/3.0},
    'Perseus': {'M_bar': 8.0e13*Msun, 'M_obs': 5.5e14*Msun, 'R500': 1100*kpc,
                'sigma': 900e3,  'r_c': 200*kpc, 'beta': 2.0/3.0},
    'Virgo':   {'M_bar': 4.0e13*Msun, 'M_obs': 3.0e14*Msun, 'R500': 900*kpc,
                'sigma': 700e3,  'r_c': 180*kpc, 'beta': 2.0/3.0},
    'A1689':   {'M_bar': 2.0e14*Msun, 'M_obs': 1.2e15*Msun, 'R500': 1300*kpc,
                'sigma': 1200e3, 'r_c': 300*kpc, 'beta': 2.0/3.0},
    'Bullet':  {'M_bar': 3.4e14*Msun, 'M_obs': 5.5e14*Msun, 'R500': 1000*kpc,
                'sigma': 1000e3, 'r_c': 200*kpc, 'beta': 2.0/3.0},
}

# ============================================================
# MAIN
# ============================================================

def process_cluster(name, params):
    """Process one cluster through Parts A-E."""
    M_bar_total = params['M_bar']
    M_obs       = params['M_obs']
    R500_m      = params['R500']
    sigma       = params['sigma']
    r_c         = params['r_c']
    beta        = params['beta']

    print(f"\n{'='*70}")
    print(f"  CLUSTER: {name}")
    print(f"{'='*70}")
    print(f"  M_bar     = {M_bar_total/Msun:.2e} Msun")
    print(f"  M_obs     = {M_obs/Msun:.2e} Msun")
    print(f"  R500      = {R500_m/kpc:.0f} kpc")
    print(f"  sigma     = {sigma/1e3:.0f} km/s")
    print(f"  r_c       = {r_c/kpc:.0f} kpc")
    print(f"  beta      = {beta:.3f}")

    # --- PART A: Set up cluster model ---
    rho0 = beta_rho0(M_bar_total, R500_m, r_c, beta)
    r_grid = make_radial_grid(R500_m, n=500)
    M_bar = compute_baryonic_mass_profile(r_grid, rho0, r_c, beta)

    # Verify normalization
    M_bar_R500 = interp_mass_at_R500(r_grid, M_bar, R500_m)
    print(f"\n  [Part A] Baryonic mass profile:")
    print(f"    rho_0 (central)   = {rho0:.4e} kg/m^3")
    print(f"    M_bar(R500) check = {M_bar_R500/Msun:.3e} Msun (target: {M_bar_total/Msun:.2e})")

    # TGP potential (baryons only)
    Phi_baryons = compute_tgp_potential(r_grid, M_bar)
    Phi_center = Phi_baryons[0]
    Phi_R500 = np.interp(R500_m, r_grid, Phi_baryons)
    Phi_boundary = np.interp(2.0*R500_m, r_grid, Phi_baryons)
    print(f"    Phi_TGP(center)   = {Phi_center:.4e} m^2/s^2")
    print(f"    Phi_TGP(R500)     = {Phi_R500:.4e} m^2/s^2")
    print(f"    Phi_TGP(2*R500)   = {Phi_boundary:.4e} m^2/s^2")
    print(f"    Delta_Phi(center) = {Phi_center - Phi_boundary:.4e} m^2/s^2")
    print(f"    Delta_Phi/sigma^2 = {(Phi_center - Phi_boundary)/sigma**2:.3f}")

    # TGP dynamical mass (what TGP predicts without neutrinos)
    g_N = G * M_bar / r_grid**2
    y = g_N / a0
    nu_val = nu_tgp(y)
    M_TGP_dyn = nu_val * M_bar  # effective dynamical mass at each r
    M_TGP_R500 = interp_mass_at_R500(r_grid, M_TGP_dyn, R500_m)
    M_deficit = M_obs - M_TGP_R500
    print(f"\n    M_TGP_dyn(R500)   = {M_TGP_R500/Msun:.3e} Msun")
    print(f"    M_deficit         = {M_deficit/Msun:.3e} Msun")
    print(f"    Deficit fraction  = {M_deficit/M_obs*100:.1f}%")

    # --- PART B: Non-self-consistent neutrino profiles ---
    print(f"\n  [Part B] Non-self-consistent neutrino mass at R500:")
    print(f"    {'m_s (eV)':>10} {'M_nu_Boltz':>14} {'M_nu_Fermi':>14} {'rho_TG':>14} {'rho_center':>14}")
    print(f"    {'':>10} {'(Msun)':>14} {'(Msun)':>14} {'(kg/m^3)':>14} {'(kg/m^3)':>14}")

    results_B = {}
    for m_s in m_s_list:
        rho_boltz = boltzmann_density(r_grid, Phi_baryons, m_s, sigma, R500_m)
        rho_fermi = fermi_density(r_grid, Phi_baryons, m_s)
        M_nu_boltz = enclosed_mass_from_density(r_grid, rho_boltz)
        M_nu_fermi = enclosed_mass_from_density(r_grid, rho_fermi)
        M_boltz_R500 = interp_mass_at_R500(r_grid, M_nu_boltz, R500_m)
        M_fermi_R500 = interp_mass_at_R500(r_grid, M_nu_fermi, R500_m)
        rho_TG = tremaine_gunn_density(m_s, sigma)
        print(f"    {m_s:10.1f} {M_boltz_R500/Msun:14.3e} {M_fermi_R500/Msun:14.3e} "
              f"{rho_TG:14.3e} {rho_boltz[0]:14.3e}")
        results_B[m_s] = {'M_boltz': M_boltz_R500, 'M_fermi': M_fermi_R500}

    # --- PART C: Self-consistent iteration ---
    print(f"\n  [Part C] Self-consistent iteration results:")
    print(f"    {'m_s (eV)':>10} {'M_nu(R500)':>14} {'Iterations':>12} {'Deficit filled':>16}")

    results_C = {}
    rho_profiles = {}
    for m_s in m_s_list:
        rho_nu, M_nu, M_nu_R500, n_iter = self_consistent_neutrino(
            r_grid, M_bar, R500_m, m_s, sigma)
        if M_deficit > 0:
            frac = M_nu_R500 / M_deficit * 100.0
        else:
            frac = 0.0  # negative deficit (Bullet overproduces)
        print(f"    {m_s:10.1f} {M_nu_R500/Msun:14.3e} {n_iter:12d} {frac:15.1f}%")
        results_C[m_s] = {'M_nu': M_nu_R500, 'frac': frac, 'n_iter': n_iter}
        rho_profiles[m_s] = rho_nu

    # --- PART D: Compare with TGP deficit ---
    print(f"\n  [Part D] Deficit analysis:")
    print(f"    M_obs             = {M_obs/Msun:.3e} Msun")
    print(f"    M_TGP_dyn(R500)   = {M_TGP_R500/Msun:.3e} Msun")
    print(f"    M_deficit         = {M_deficit/Msun:.3e} Msun")
    if M_deficit > 0:
        print(f"    Deficit fraction  = {M_deficit/M_obs*100:.1f}%")
    else:
        print(f"    NOTE: TGP already EXCEEDS M_obs (no deficit)")

    min_ms_90 = None
    for m_s in m_s_list:
        M_nu_R500 = results_C[m_s]['M_nu']
        M_total = M_TGP_R500 + M_nu_R500
        if M_deficit > 0:
            frac = M_nu_R500 / M_deficit * 100.0
        else:
            frac = 0.0
        print(f"    m_s={m_s:5.1f} eV: M_total = {M_total/Msun:.3e}, "
              f"deficit filled = {frac:.1f}%")
        if M_deficit > 0 and frac >= 90.0 and min_ms_90 is None:
            min_ms_90 = m_s

    if M_deficit <= 0:
        print(f"    => No deficit to fill (TGP alone overshoots)")
        min_ms_90 = 0.0  # sentinel: no neutrinos needed
    elif min_ms_90 is not None:
        print(f"    => Minimum m_s for >90% fill: {min_ms_90} eV")
    else:
        print(f"    => No tested mass fills >90% of deficit")

    # --- PART E: Density profiles (Coma only) ---
    if name == 'Coma':
        print(f"\n  [Part E] Neutrino density profiles rho_nu(r) [kg/m^3]:")
        header = f"    {'r (kpc)':>10}"
        for m_s in m_s_list:
            header += f" {'m='+str(m_s)+'eV':>14}"
        header += f" {'NFW_DM':>14}"
        print(header)

        # NFW reference
        r_s_nfw = 300.0 * kpc
        rho_nfw = nfw_density(r_grid, M_obs, r_s_nfw)

        for r_kpc in r_profile_kpc:
            r_m = r_kpc * kpc
            line = f"    {r_kpc:10.0f}"
            for m_s in m_s_list:
                rho_val = np.interp(r_m, r_grid, rho_profiles[m_s])
                line += f" {rho_val:14.4e}"
            rho_nfw_val = np.interp(r_m, r_grid, rho_nfw)
            line += f" {rho_nfw_val:14.4e}"
            print(line)

    return {
        'r_grid': r_grid,
        'M_bar': M_bar,
        'M_TGP_R500': M_TGP_R500,
        'M_deficit': M_deficit,
        'results_B': results_B,
        'results_C': results_C,
        'rho_profiles': rho_profiles,
        'min_ms_90': min_ms_90
    }


def main():
    print("=" * 70)
    print("  gs59_neutrino_halo.py")
    print("  Full Poisson-Boltzmann neutrino halo calculation for TGP")
    print("  Sterile neutrino halos in galaxy clusters")
    print("=" * 70)
    print(f"\n  Cosmic neutrino density: n = {n_cosmic:.2e} /m^3 = {n_cosmic/1e6:.0f} /cm^3")
    print(f"  Cosmic neutrino T: {T_nu} K")
    print(f"  a0 = {a0:.1e} m/s^2")
    print(f"  Spin states g_s = {g_s}")

    # Check key scales
    for m_s in m_s_list:
        m_kg = m_s * eV_kg
        rho_cosmic = n_cosmic * m_kg
        rho_TG_1000 = tremaine_gunn_density(m_s, 1000e3)
        print(f"\n  m_s = {m_s} eV:")
        print(f"    rho_cosmic     = {rho_cosmic:.3e} kg/m^3")
        print(f"    rho_TG(1000)   = {rho_TG_1000:.3e} kg/m^3")
        print(f"    enhancement    = {rho_TG_1000/rho_cosmic:.3e}")

    # Process all clusters
    all_results = {}
    cluster_order = ['Coma', 'Perseus', 'Virgo', 'A1689', 'Bullet']
    for name in cluster_order:
        all_results[name] = process_cluster(name, clusters[name])

    # ============================================================
    # PART G: Summary
    # ============================================================
    print("\n" + "=" * 70)
    print("  PART G: SUMMARY TABLE")
    print("=" * 70)

    print(f"\n  {'Cluster':>10} {'M_deficit':>14} {'M_nu(2eV)':>14} {'M_nu(11eV)':>14} "
          f"{'fill(2eV)':>12} {'fill(11eV)':>12} {'min_ms_90':>10}")
    print(f"  {'':>10} {'(1e14 Msun)':>14} {'(1e14 Msun)':>14} {'(1e14 Msun)':>14} "
          f"{'(%)':>12} {'(%)':>12} {'(eV)':>10}")
    print("  " + "-" * 90)

    for name in cluster_order:
        res = all_results[name]
        M_def = res['M_deficit']
        rc = res['results_C']
        M_2  = rc[2.0]['M_nu']
        M_11 = rc[11.0]['M_nu']
        f_2  = rc[2.0]['frac']
        f_11 = rc[11.0]['frac']
        ms90 = res['min_ms_90']
        if ms90 is not None and ms90 == 0.0:
            ms90_str = "N/A*"  # no deficit
        elif ms90 is not None:
            ms90_str = f"{ms90}"
        else:
            ms90_str = ">20"
        print(f"  {name:>10} {M_def/Msun/1e14:14.3f} {M_2/Msun/1e14:14.3f} "
              f"{M_11/Msun/1e14:14.3f} {f_2:12.1f} {f_11:12.1f} {ms90_str:>10}")

    print(f"\n  * N/A means TGP alone already exceeds M_obs (no deficit)")

    # Find minimum m_s that works for ALL clusters (excl. those with no deficit)
    clusters_with_deficit = [n for n in cluster_order
                             if all_results[n]['M_deficit'] > 0]
    print(f"\n  Clusters with positive deficit: {', '.join(clusters_with_deficit)}")

    print(f"\n  Global minimum m_s for >90% deficit fill (clusters with deficit):")
    global_min = None
    for m_s in m_s_list:
        all_fill = True
        for name in clusters_with_deficit:
            frac = all_results[name]['results_C'][m_s]['frac']
            if frac < 90.0:
                all_fill = False
                break
        if all_fill:
            global_min = m_s
            print(f"    => m_s = {m_s} eV fills >90% for ALL deficit clusters")
            break

    if global_min is None:
        print(f"    => No tested mass fills >90% for ALL deficit clusters simultaneously")
        for name in clusters_with_deficit:
            frac_20 = all_results[name]['results_C'][20.0]['frac']
            print(f"       {name}: {frac_20:.1f}% filled at 20 eV")

    # Also including ALL clusters
    print(f"\n  Global minimum m_s for >90% deficit fill (ALL clusters incl Bullet):")
    for m_s in m_s_list:
        all_fill = True
        for name in cluster_order:
            res = all_results[name]
            if res['M_deficit'] <= 0:
                continue  # no deficit needed
            frac = res['results_C'][m_s]['frac']
            if frac < 90.0:
                all_fill = False
                break
        if all_fill:
            print(f"    => m_s = {m_s} eV fills >90% for ALL clusters")
            break
    else:
        print(f"    => No tested mass fills >90% for ALL clusters")

    # Comparison with literature
    print(f"\n  Comparison with literature:")
    print(f"    Sanders (2003): ~2 eV neutrinos can partially explain cluster mass")
    print(f"      in MOND. Central densities reach TG limit for m > few eV.")
    print(f"    Angus (2009): 11 eV sterile neutrinos can explain cluster masses")
    print(f"      in MOND with full self-consistent Poisson-Boltzmann calculation.")
    print(f"")
    print(f"    TGP advantage: smaller deficit (~32% vs ~40% for MOND)")
    print(f"    => TGP needs LESS neutrino mass to close the gap.")

    # TGP + 2eV assessment
    print(f"\n  Assessment: Is TGP + 2 eV sterile neutrino sufficient?")
    all_ok_2eV = all(
        all_results[n]['M_deficit'] <= 0 or all_results[n]['results_C'][2.0]['frac'] >= 90.0
        for n in cluster_order
    )
    if all_ok_2eV:
        print(f"    YES - 2 eV sterile neutrinos fill >90% of the deficit")
        print(f"    for ALL clusters. TGP+2eV is viable.")
    else:
        print(f"    NO - 2 eV is insufficient for some clusters.")
        print(f"    Cluster-by-cluster at 2 eV:")
        for name in cluster_order:
            res = all_results[name]
            if res['M_deficit'] <= 0:
                print(f"      {name:>10}: no deficit (TGP exceeds M_obs)")
            else:
                f2 = res['results_C'][2.0]['frac']
                status = "OK" if f2 >= 90 else "INSUFFICIENT"
                print(f"      {name:>10}: {f2:.1f}% [{status}]")

        # Find what mass IS sufficient
        for m_s in m_s_list:
            ok = all(
                all_results[n]['M_deficit'] <= 0 or
                all_results[n]['results_C'][m_s]['frac'] >= 90.0
                for n in cluster_order
            )
            if ok:
                print(f"\n    => Minimum viable m_s for all clusters: {m_s} eV")
                if m_s < 11:
                    print(f"    This is BETTER than Angus (2009) 11 eV requirement in MOND.")
                elif m_s == 11:
                    print(f"    This is comparable to Angus (2009) 11 eV requirement in MOND.")
                else:
                    print(f"    This is worse than Angus (2009) 11 eV requirement in MOND.")
                break
        else:
            print(f"\n    => Even 20 eV is not enough for all clusters.")
            print(f"    Additional physics may be needed (e.g., non-thermal neutrino")
            print(f"    distributions, cluster-specific baryon fractions, or")
            print(f"    higher-order TGP effects).")

    print(f"\n{'='*70}")
    print(f"  gs59 complete.")
    print(f"{'='*70}")


if __name__ == '__main__':
    main()
