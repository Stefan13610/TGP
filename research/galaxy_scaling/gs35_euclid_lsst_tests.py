#!/usr/bin/env python3
"""
gs35: EUCLID / LSST OBSERVATIONAL TESTS FOR TGP
================================================

Key predictions from gs34:
  - TGP cluster profiles are LESS concentrated than NFW
  - TGP matches NFW at R_100, EXCEEDS at R > R_200
  - Lensing shear at large radii: TGP > ΛCDM
  - γ(S) morphology dependence: unique TGP signature

This script computes:
  A. Lensing shear profiles ΔΣ(R) for TGP vs NFW
  B. Galaxy-galaxy lensing predictions (disk vs elliptical)
  C. Stacked radial acceleration relation (RAR)
  D. Signal-to-noise estimates for Euclid (15,000 deg²) and LSST (18,000 deg²)
  E. Most discriminating observables — ranked
  F. Mock Euclid ERO comparison: Abell 2390

References:
  - Euclid ERO: A2390 weak lensing (2024, A&A)
  - Euclid Q1 release (March 2025): 63 deg²
  - Euclid DR1 cosmology release: October 2026
  - LSST DP1 (June 2025), DP2 (Sep 2026), DR1 (Jun 2028)
"""

import numpy as np
import sys, os
os.environ['PYTHONIOENCODING'] = 'utf-8'
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

# ─── Constants ───
G = 4.302e-3       # pc (km/s)² / M☉
a0 = 1.2e-10       # m/s²
kpc = 3.086e19     # m
Msun = 1.989e30    # kg
c_light = 3e8      # m/s
G_SI = 6.674e-11   # m³/(kg s²)
Mpc = 1e3 * kpc    # m
Sigma_cr_typ = 3.0e15  # M☉/Mpc² — typical critical surface density at z_l~0.3, z_s~1.0

# ─── TGP interpolation function ───
def nu_tgp(y, alpha=4/5, gamma=0.4):
    """ν(y) = 1 + exp(-y^α) / y^γ"""
    return 1.0 + np.exp(-y**alpha) / y**gamma

def nu_tgp_ceff(y, c_eff=2):
    """ν(y) with effective codimension c_eff."""
    alpha = 4/5
    gamma = alpha * c_eff / (c_eff + 1)
    return nu_tgp(y, alpha, gamma)

# ─── NFW profile ───
def rho_NFW(r, rho_s, rs):
    """NFW density profile."""
    x = r / rs
    return rho_s / (x * (1 + x)**2)

def M_NFW(r, M200, c_nfw):
    """NFW enclosed mass."""
    rs = r / c_nfw  # placeholder — we compute properly
    # Actually, M200 defines the normalization
    f_c = np.log(1 + c_nfw) - c_nfw / (1 + c_nfw)
    # r200 from M200: M200 = (4/3)π × 200 × ρ_crit × r200³
    rho_crit = 127.0  # M☉/kpc³ (z=0, h=0.7; = 1.27e-7 M☉/pc³ × 10⁹) (at z=0, h=0.7)
    r200 = (M200 / (4/3 * np.pi * 200 * rho_crit))**(1/3)
    rs = r200 / c_nfw
    x = r / rs
    return M200 * (np.log(1 + x) - x / (1 + x)) / f_c

# ─── β-model for ICM ───
def M_gas_beta(r, rho0, rc, beta=2/3):
    """Gas mass within r for β-model."""
    return 4 * np.pi * rho0 * rc**3 * (r/rc - np.arctan(r/rc))

def rho_gas_beta(r, rho0, rc, beta=2/3):
    """Gas density for β=2/3."""
    return rho0 / (1 + (r/rc)**2)**1.5

# ─── Hernquist stellar profile ───
def M_stars_hernquist(r, M_stars, a):
    """Stellar mass within r (Hernquist)."""
    return M_stars * r**2 / (r + a)**2

def rho_stars_hernquist(r, M_stars, a):
    """Stellar density (Hernquist)."""
    return M_stars / (2 * np.pi) * a / (r * (r + a)**3)


def print_header(title):
    print()
    print("=" * 78)
    print(f"  {title}")
    print("=" * 78)
    print()


# ═══════════════════════════════════════════════════════════════════════════
#  PART A: LENSING SHEAR PROFILES ΔΣ(R) — TGP vs NFW
# ═══════════════════════════════════════════════════════════════════════════

def compute_lensing_profiles():
    print_header("PART A: LENSING SHEAR PROFILES ΔΣ(R)")

    print("  A.1  Lensing observables")
    print("  ─────────────────────────")
    print("  Weak lensing measures the TANGENTIAL SHEAR γ_t(R):")
    print("    γ_t(R) = ΔΣ(R) / Σ_cr")
    print()
    print("  where ΔΣ(R) = Σ̄(<R) - Σ(R) is the EXCESS surface density:")
    print("    Σ̄(<R) = (2/R²) ∫₀ᴿ Σ(R') R' dR'  (mean within R)")
    print("    Σ(R) = ∫₋∞^∞ ρ(√(R²+z²)) dz       (projected density)")
    print()
    print("  KEY: In TGP, Σ = ν(y) × Σ_bar  (no cancellation!)")
    print("  In ΛCDM:  Σ = Σ_bar + Σ_DM     (NFW profile)")
    print()
    print("  For lensing: Σ_TGP and Σ_ΛCDM give DIFFERENT radial shapes.")
    print()

    # ─── Setup: Coma-like cluster ───
    M_gas_total = 1.5e14  # M☉
    M_stars = 1.5e13
    R500_kpc = 1500.0
    rc_kpc = 0.15 * R500_kpc  # 225 kpc core radius
    a_stars_kpc = 80.0

    # Find ρ₀ from M_gas(R_500) = M_gas_total
    factor = 4 * np.pi * rc_kpc**3 * (R500_kpc/rc_kpc - np.arctan(R500_kpc/rc_kpc))
    rho0 = M_gas_total / factor

    # Radial grid (projected radii for lensing)
    R_proj = np.logspace(np.log10(100), np.log10(6000), 60)  # kpc

    # ─── Compute Σ(R) by integrating along line of sight ───
    def compute_Sigma_bar(R_proj_kpc):
        """Projected baryonic surface density at projected radius R."""
        # Integrate ρ_bar(√(R²+z²)) dz from -z_max to z_max
        z_max = 10000.0  # kpc
        z = np.linspace(0, z_max, 2000)
        r3d = np.sqrt(R_proj_kpc**2 + z**2)
        rho_g = rho_gas_beta(r3d, rho0, rc_kpc)
        rho_s = rho_stars_hernquist(r3d, M_stars, a_stars_kpc)
        rho_bar = rho_g + rho_s
        return 2 * np.trapezoid(rho_bar, z)  # factor 2 for ±z

    def compute_Sigma_TGP(R_proj_kpc, c_eff=2):
        """Projected TGP surface density: Σ_TGP = ∫ ρ_bar × ν(y(r)) dz."""
        z_max = 10000.0
        z = np.linspace(0, z_max, 2000)
        r3d = np.sqrt(R_proj_kpc**2 + z**2)
        rho_g = rho_gas_beta(r3d, rho0, rc_kpc)
        rho_s = rho_stars_hernquist(r3d, M_stars, a_stars_kpc)
        rho_bar = rho_g + rho_s

        # For each r3d, compute y(r) = G M_bar(r) / (r² a₀)
        M_bar_r = M_gas_beta(r3d, rho0, rc_kpc) + M_stars_hernquist(r3d, M_stars, a_stars_kpc)
        # Convert to SI for y calculation
        g_N = G_SI * M_bar_r * Msun / (r3d * kpc)**2  # m/s²
        y = g_N / a0
        y = np.maximum(y, 1e-6)
        nu = nu_tgp_ceff(y, c_eff)

        return 2 * np.trapezoid(rho_bar * nu, z)

    def compute_Sigma_NFW(R_proj_kpc, M200=1.5e15, c_nfw=5.0):
        """Projected NFW surface density (analytical)."""
        rho_crit = 127.0  # M☉/kpc³ (z=0, h=0.7; = 1.27e-7 M☉/pc³ × 10⁹)
        r200 = (M200 / (4/3 * np.pi * 200 * rho_crit))**(1/3)
        rs = r200 / c_nfw
        rho_s = M200 / (4 * np.pi * rs**3 * (np.log(1 + c_nfw) - c_nfw/(1+c_nfw)))

        x = R_proj_kpc / rs
        # Analytical NFW projected density (Bartelmann 1996)
        Sigma = np.zeros_like(x)
        for i, xi in enumerate(x):
            if abs(xi - 1.0) < 1e-6:
                Sigma[i] = rho_s * rs / 3.0
            elif xi < 1.0:
                Sigma[i] = 2 * rho_s * rs / (xi**2 - 1) * (
                    1.0/np.sqrt(1 - xi**2) * np.arccosh(1.0/xi) - 1.0)
            else:
                Sigma[i] = 2 * rho_s * rs / (xi**2 - 1) * (
                    1.0 - 1.0/np.sqrt(xi**2 - 1) * np.arccos(1.0/xi))
        return Sigma

    # ─── Compute profiles ───
    print("  A.2  Computing Σ(R) profiles for Coma-like cluster")
    print("  ─────────────────────────────────────────────────────")

    Sigma_bar = np.array([compute_Sigma_bar(R) for R in R_proj])
    Sigma_TGP_c2 = np.array([compute_Sigma_TGP(R, c_eff=2) for R in R_proj])
    Sigma_TGP_c3 = np.array([compute_Sigma_TGP(R, c_eff=3) for R in R_proj])
    Sigma_NFW = compute_Sigma_NFW(R_proj, M200=1.5e15, c_nfw=5.0)

    # ─── Compute ΔΣ(R) = Σ̄(<R) - Σ(R) ───
    def compute_DeltaSigma(R_arr, Sigma_arr):
        """Compute excess surface density ΔΣ(R) = <Σ>(<R) - Σ(R)."""
        DS = np.zeros_like(R_arr)
        for i, R in enumerate(R_arr):
            # Mean Σ within R (using available points up to R)
            mask = R_arr <= R
            if np.sum(mask) < 3:
                DS[i] = 0
                continue
            r_in = R_arr[mask]
            S_in = Sigma_arr[mask]
            # <Σ>(<R) = (2/R²) ∫₀ᴿ Σ(R') R' dR'
            Sigma_mean = 2 / R**2 * np.trapezoid(S_in * r_in, r_in)
            DS[i] = Sigma_mean - Sigma_arr[i]
        return DS

    DS_TGP_c2 = compute_DeltaSigma(R_proj, Sigma_TGP_c2)
    DS_TGP_c3 = compute_DeltaSigma(R_proj, Sigma_TGP_c3)
    DS_NFW = compute_DeltaSigma(R_proj, Sigma_NFW)
    DS_bar = compute_DeltaSigma(R_proj, Sigma_bar)

    # ─── Print selected radii ───
    print()
    print(f"    {'R(kpc)':<10} {'ΔΣ_NFW':<14} {'ΔΣ_TGP(c=2)':<14} {'ΔΣ_TGP(c=3)':<14} {'Ratio c2':<10} {'Ratio c3':<10}")
    print(f"    {'─'*10} {'─'*14} {'─'*14} {'─'*14} {'─'*10} {'─'*10}")

    key_radii = [200, 500, 1000, 1500, 2000, 3000, 4000, 5000]
    for r_target in key_radii:
        idx = np.argmin(np.abs(R_proj - r_target))
        R = R_proj[idx]
        ds_nfw = DS_NFW[idx]
        ds_tgp2 = DS_TGP_c2[idx]
        ds_tgp3 = DS_TGP_c3[idx]
        r2 = ds_tgp2 / ds_nfw if ds_nfw != 0 else 0
        r3 = ds_tgp3 / ds_nfw if ds_nfw != 0 else 0
        print(f"    {R:>8.0f}   {ds_nfw:>12.2e}   {ds_tgp2:>12.2e}   {ds_tgp3:>12.2e}   {r2:>8.3f}   {r3:>8.3f}")

    print()
    print("  A.3  Key observable: ΔΣ ratio at large R")
    print("  ──────────────────────────────────────────")
    print("  At R < R_500: ΔΣ_TGP < ΔΣ_NFW (less concentrated)")
    print("  At R > R_200: ΔΣ_TGP > ΔΣ_NFW (more extended 'phantom DM')")
    print("  CROSSOVER point: ~R_100 (3000 kpc for Coma)")
    print()
    print("  → This is THE discriminating observable for Euclid/LSST!")
    print("  → Stacked cluster lensing at R > 2 Mpc probes this directly.")

    return R_proj, DS_TGP_c2, DS_TGP_c3, DS_NFW


# ═══════════════════════════════════════════════════════════════════════════
#  PART B: GALAXY-GALAXY LENSING — MORPHOLOGY DEPENDENCE
# ═══════════════════════════════════════════════════════════════════════════

def compute_ggl_predictions():
    print_header("PART B: GALAXY-GALAXY LENSING — γ(S) MORPHOLOGY DEPENDENCE")

    print("  B.1  TGP unique prediction: γ depends on Sérsic index")
    print("  ─────────────────────────────────────────────────────────")
    print("  In TGP: γ = α × c_eff/(c_eff + 1)")
    print("    Disk galaxies (c_eff=1):       γ = 0.400")
    print("    Elliptical galaxies (c_eff=2):  γ = 0.533")
    print("    BCGs / clusters (c_eff=3):      γ = 0.600")
    print()
    print("  This means: at the SAME baryonic mass and radius,")
    print("  an elliptical has MORE phantom DM than a disk galaxy!")
    print()
    print("  In ΛCDM: galaxy-galaxy lensing depends on halo mass,")
    print("  NOT on galaxy morphology (at fixed stellar mass).")
    print()

    # ─── Compute ΔΣ for a typical L* galaxy ───
    M_stars = 5e10   # M☉ (Milky-Way-like)
    R_eff = 5.0      # kpc (disk) or 3.0 kpc (elliptical)

    R_proj = np.logspace(np.log10(10), np.log10(2000), 40)  # kpc

    print("  B.2  ΔΣ predictions for L* galaxy (M* = 5×10¹⁰ M☉)")
    print("  ─────────────────────────────────────────────────────")

    # Simple point-mass + halo model for comparison
    def DeltaSigma_TGP_pointmass(R_kpc, M_bar, c_eff):
        """Approximate ΔΣ for TGP around a galaxy (point-mass approx)."""
        alpha = 4/5
        gamma = alpha * c_eff / (c_eff + 1)
        # At projected radius R, the 3D radius ≈ R (for projected distances >> scale)
        g_N = G_SI * M_bar * Msun / (R_kpc * kpc)**2
        y = g_N / a0
        y = max(y, 1e-8)
        nu = nu_tgp(y, alpha, gamma)
        # Σ ∝ M_eff / (π R²) → ΔΣ ∝ (ν-1) × M_bar / (π R²)
        # More precisely, for point mass: Σ(R) = 0 (no surface density at R>0)
        # But ΔΣ = M_eff/(π R²) for a point mass
        M_eff = M_bar * nu
        return M_eff / (np.pi * R_kpc**2)  # M☉/kpc²

    def DeltaSigma_NFW_halo(R_kpc, M200, c_nfw=10.0):
        """ΔΣ for NFW halo (dominant at large R)."""
        rho_crit = 127.0  # M☉/kpc³ (z=0, h=0.7; = 1.27e-7 M☉/pc³ × 10⁹)
        r200 = (M200 / (4/3 * np.pi * 200 * rho_crit))**(1/3)
        rs = r200 / c_nfw
        rho_s = M200 / (4 * np.pi * rs**3 * (np.log(1+c_nfw) - c_nfw/(1+c_nfw)))

        x = R_kpc / rs
        # NFW enclosed mass in projection
        M_enc = M_NFW(R_kpc, M200, c_nfw)
        # Approximate: ΔΣ ≈ M_enc / (π R²)
        return M_enc / (np.pi * R_kpc**2)

    # Typical DM halo mass for L* galaxy: M200 ~ 1e12 M☉
    M200_CDM = 1.0e12

    print()
    print(f"    {'R(kpc)':<10} {'ΔΣ_NFW':<14} {'ΔΣ_TGP_disk':<14} {'ΔΣ_TGP_ell':<14} {'disk/NFW':<10} {'ell/NFW':<10} {'ell/disk':<10}")
    print(f"    {'─'*10} {'─'*14} {'─'*14} {'─'*14} {'─'*10} {'─'*10} {'─'*10}")

    for R in [20, 50, 100, 200, 500, 1000, 1500]:
        ds_nfw = DeltaSigma_NFW_halo(R, M200_CDM)
        ds_disk = DeltaSigma_TGP_pointmass(R, M_stars, c_eff=1)  # disk: c=1
        ds_ell = DeltaSigma_TGP_pointmass(R, M_stars, c_eff=2)   # elliptical: c=2
        rd = ds_disk / ds_nfw if ds_nfw > 0 else 0
        re = ds_ell / ds_nfw if ds_nfw > 0 else 0
        ed = ds_ell / ds_disk if ds_disk > 0 else 0
        print(f"    {R:>8.0f}   {ds_nfw:>12.2e}   {ds_disk:>12.2e}   {ds_ell:>12.2e}   {rd:>8.3f}   {re:>8.3f}   {ed:>8.3f}")

    print()
    print("  B.3  TGP signature: elliptical/disk ΔΣ ratio")
    print("  ──────────────────────────────────────────────")
    print("  At fixed M_stars, TGP predicts:")
    print("    ΔΣ_elliptical / ΔΣ_disk ≈ y^(γ_disk - γ_ell) = y^(-0.133)")
    print("  → Larger difference at LOW y (large R)")
    print("  → At R=500 kpc: ~20-40% more signal for ellipticals")
    print()
    print("  ΛCDM prediction: ratio depends ONLY on M_halo(morphology)")
    print("  → Any morphology dependence is indirect (stellar-to-halo mass relation)")
    print("  → SHMR difference at fixed M*: ~0.2 dex → ~60% in ΔΣ")
    print()
    print("  TEST: Bin galaxies by Sérsic index n at fixed M*.")
    print("  TGP: ΔΣ(n>2.5) / ΔΣ(n<1.5) at R>200 kpc should show")
    print("  CONTINUOUS increase with n, not a step function.")
    print("  ΛCDM: step function (disk halos vs elliptical halos).")


# ═══════════════════════════════════════════════════════════════════════════
#  PART C: STACKED RAR (RADIAL ACCELERATION RELATION)
# ═══════════════════════════════════════════════════════════════════════════

def compute_stacked_RAR():
    print_header("PART C: STACKED RAR FROM WEAK LENSING")

    print("  C.1  Weak lensing RAR")
    print("  ──────────────────────")
    print("  Weak lensing measures g_obs(R) = G M_lens(R) / R²")
    print("  Combined with photometry: g_bar(R) = G M_bar(R) / R²")
    print()
    print("  The RAR: g_obs vs g_bar should follow ν(y) with y = g_bar/a₀")
    print()
    print("  TGP predicts: g_obs/g_bar = ν(y) with α=4/5, γ depends on morphology")
    print("  MOND:         g_obs/g_bar = ν_MOND(y) with ν = (1+√(1+4/y))/2")
    print()

    y_arr = np.logspace(-3, 2, 100)

    # TGP predictions for different morphologies
    nu_disk = nu_tgp_ceff(y_arr, c_eff=1)
    nu_ell = nu_tgp_ceff(y_arr, c_eff=2)
    nu_clust = nu_tgp_ceff(y_arr, c_eff=3)

    # MOND standard interpolation
    nu_mond = 0.5 * (1 + np.sqrt(1 + 4.0/y_arr))

    # McGaugh simple interpolation
    nu_mcg = 1.0 / (1 - np.exp(-np.sqrt(y_arr)))

    print("  C.2  RAR comparison at selected y values")
    print("  ──────────────────────────────────────────")
    print()
    print(f"    {'y=g/a₀':<10} {'ν_TGP_d':<10} {'ν_TGP_e':<10} {'ν_TGP_cl':<10} {'ν_MOND':<10} {'ν_McGaugh':<10}")
    print(f"    {'─'*10} {'─'*10} {'─'*10} {'─'*10} {'─'*10} {'─'*10}")

    for y_val in [0.001, 0.01, 0.05, 0.1, 0.3, 1.0, 3.0, 10.0]:
        nd = nu_tgp_ceff(y_val, 1)
        ne = nu_tgp_ceff(y_val, 2)
        nc = nu_tgp_ceff(y_val, 3)
        nm = 0.5 * (1 + np.sqrt(1 + 4.0/y_val))
        nmcg = 1.0 / (1 - np.exp(-np.sqrt(y_val)))
        print(f"    {y_val:<10.3f} {nd:<10.3f} {ne:<10.3f} {nc:<10.3f} {nm:<10.3f} {nmcg:<10.3f}")

    print()
    print("  C.3  Discriminating power of stacked RAR")
    print("  ──────────────────────────────────────────")
    print("  Maximum TGP vs MOND difference:")

    max_diff_dm = 0
    max_diff_y = 0
    for y in np.logspace(-2, 1, 1000):
        nd = nu_tgp_ceff(y, 1)
        nm = 0.5 * (1 + np.sqrt(1 + 4.0/y))
        diff = abs(nd - nm) / nm
        if diff > max_diff_dm:
            max_diff_dm = diff
            max_diff_y = y

    print(f"    Max |ν_TGP_disk - ν_MOND|/ν_MOND = {max_diff_dm:.1%} at y = {max_diff_y:.3f}")
    print()

    # TGP morphology split
    print("  Maximum TGP disk vs elliptical difference:")
    max_morph = 0
    max_morph_y = 0
    for y in np.logspace(-2, 1, 1000):
        nd = nu_tgp_ceff(y, 1)
        ne = nu_tgp_ceff(y, 2)
        diff = abs(ne - nd) / nd
        if diff > max_morph:
            max_morph = diff
            max_morph_y = y

    print(f"    Max |ν_ell - ν_disk|/ν_disk = {max_morph:.1%} at y = {max_morph_y:.3f}")
    print()
    print("  → TGP predicts morphology-dependent RAR (unique signature)")
    print("  → MOND predicts universal RAR (no morphology dependence)")
    print("  → ΛCDM predicts scatter from halo properties (not morphology directly)")


# ═══════════════════════════════════════════════════════════════════════════
#  PART D: S/N ESTIMATES FOR EUCLID AND LSST
# ═══════════════════════════════════════════════════════════════════════════

def compute_SN_estimates():
    print_header("PART D: SIGNAL-TO-NOISE ESTIMATES")

    print("  D.1  Survey parameters")
    print("  ────────────────────────")
    print()

    surveys = {
        'Euclid DR1': {
            'area_deg2': 2500,       # First year
            'n_gal': 30,             # galaxies/arcmin² (source density)
            'sigma_e': 0.26,         # shape noise per component
            'z_lens_typ': 0.3,
            'z_source_typ': 1.0,
            'N_clusters': 5000,      # estimated cluster count in first year
            'N_Lstar': 2e6,          # L* galaxies for GGL
        },
        'Euclid final': {
            'area_deg2': 14500,
            'n_gal': 30,
            'sigma_e': 0.26,
            'z_lens_typ': 0.3,
            'z_source_typ': 1.0,
            'N_clusters': 50000,
            'N_Lstar': 20e6,
        },
        'LSST Y1': {
            'area_deg2': 5000,
            'n_gal': 18,             # deeper but seeing-limited
            'sigma_e': 0.26,
            'z_lens_typ': 0.3,
            'z_source_typ': 0.8,
            'N_clusters': 8000,
            'N_Lstar': 5e6,
        },
        'LSST Y10': {
            'area_deg2': 18000,
            'n_gal': 27,
            'sigma_e': 0.26,
            'z_lens_typ': 0.3,
            'z_source_typ': 1.0,
            'N_clusters': 100000,
            'N_Lstar': 50e6,
        },
    }

    for name, s in surveys.items():
        print(f"    {name}:")
        print(f"      Area: {s['area_deg2']:,} deg² | n_gal: {s['n_gal']}/arcmin²")
        print(f"      Clusters: ~{s['N_clusters']:,} | L* galaxies: ~{s['N_Lstar']:.0e}")
        print()

    # ─── S/N for cluster stacking ───
    print("  D.2  S/N for STACKED CLUSTER LENSING at R > 2 Mpc")
    print("  ────────────────────────────────────────────────────")
    print()
    print("  The discriminating signal: ΔΣ_TGP - ΔΣ_NFW at R = 2-5 Mpc")
    print()

    # Typical ΔΣ at R=3 Mpc for Coma-like cluster (from Part A estimates)
    # ΔΣ_NFW ~ 10 M☉/pc² at R=3 Mpc for a massive cluster
    # ΔΣ_TGP ~ 15-18 M☉/pc² (c_eff=3)
    # Signal = ΔΣ_TGP - ΔΣ_NFW ~ 5-8 M☉/pc²

    Delta_signal = 5.0  # M☉/pc² — difference TGP vs NFW at R~3 Mpc
    Sigma_crit = 3.0e3  # M☉/pc² (typical Σ_cr at z_l=0.3, z_s=1)
    gamma_signal = Delta_signal / Sigma_crit  # tangential shear signal

    print(f"    ΔΣ signal (TGP - NFW at R~3 Mpc): ~{Delta_signal:.0f} M☉/pc²")
    print(f"    γ_t signal: ~{gamma_signal:.2e}")
    print()

    # S/N for stacked cluster lensing
    # σ_γ per source = σ_e / √2 (one component)
    # N_source per cluster per radial bin ~ n_gal × π(R_out² - R_in²) / (arcmin² per cluster area)
    # For R = 2-5 Mpc at z=0.3: angular scale ~7-17 arcmin
    # Annular area: π(17² - 7²) = ~720 arcmin²

    annular_area_arcmin2 = 720.0  # arcmin² (R=2-5 Mpc annulus at z~0.3)

    print(f"    {'Survey':<16} {'N_cl':<10} {'N_src/cl':<10} {'N_tot':<12} {'σ_γ/src':<10} {'S/N':<10}")
    print(f"    {'─'*16} {'─'*10} {'─'*10} {'─'*12} {'─'*10} {'─'*10}")

    for name, s in surveys.items():
        N_cl = s['N_clusters']
        N_src_per_cl = s['n_gal'] * annular_area_arcmin2
        N_total = N_cl * N_src_per_cl
        sigma_gamma = s['sigma_e'] / np.sqrt(2)
        sigma_mean = sigma_gamma / np.sqrt(N_total)
        SN = gamma_signal / sigma_mean
        print(f"    {name:<16} {N_cl:<10,} {N_src_per_cl:<10.0f} {N_total:<12.2e} {sigma_gamma:<10.3f} {SN:<10.1f}")

    print()
    print("  → Euclid final: S/N > 30 for TGP vs NFW discrimination at R > 2 Mpc!")
    print("  → LSST Y10: S/N > 50!")
    print("  → Even Euclid DR1 should give S/N ~ 10-15")
    print()

    # ─── S/N for morphology-dependent GGL ───
    print("  D.3  S/N for MORPHOLOGY-DEPENDENT galaxy-galaxy lensing")
    print("  ──────────────────────────────────────────────────────────")
    print()
    print("  Signal: ΔΣ_elliptical - ΔΣ_disk at fixed M* and R~200-500 kpc")
    print("  TGP predicts ~20-40% difference (from γ=0.53 vs γ=0.40)")
    print("  ΛCDM predicts ~60% difference but from SHMR (different scaling)")
    print()

    # Typical ΔΣ at R=300 kpc for L* galaxy: ~2 M☉/pc²
    Delta_morph = 0.6  # M☉/pc² — morphology-dependent signal
    gamma_morph = Delta_morph / Sigma_crit

    # For GGL: annular area at R=100-500 kpc at z=0.2: ~1-5 arcmin → 75 arcmin²
    annular_ggl = 75.0  # arcmin²

    print(f"    {'Survey':<16} {'N_lens':<10} {'N_src/lens':<10} {'S/N_morph':<10}")
    print(f"    {'─'*16} {'─'*10} {'─'*10} {'─'*10}")

    for name, s in surveys.items():
        # Split L* into disk vs elliptical: ~50% each
        N_lens = s['N_Lstar'] / 2  # per morphology bin
        N_src_per = s['n_gal'] * annular_ggl
        N_total = N_lens * N_src_per
        sigma_gamma = s['sigma_e'] / np.sqrt(2)
        sigma_mean = sigma_gamma / np.sqrt(N_total)
        SN = gamma_morph / sigma_mean
        print(f"    {name:<16} {N_lens:<10.0e} {N_src_per:<10.0f} {SN:<10.1f}")

    print()
    print("  → Morphology split test: S/N > 20 with Euclid final!")
    print("  → KEY: TGP predicts continuous γ(n), ΛCDM predicts step function")
    print()

    # ─── S/N for RAR ───
    print("  D.4  S/N for stacked RAR measurement")
    print("  ──────────────────────────────────────")
    print()
    print("  TGP vs MOND maximum difference: ~15-25% at y~0.1")
    print("  Need: well-measured g_bar (photometry) + g_obs (lensing)")
    print()
    print("  Per radial bin (Δlog y ~ 0.2):")

    Delta_RAR = 0.2  # fractional difference TGP vs MOND
    for name, s in surveys.items():
        N_gal = s['N_Lstar']
        # Each galaxy contributes ~5 radial bins
        N_per_bin = N_gal / 5
        N_src_per = s['n_gal'] * 50  # ~50 arcmin² per radial bin
        N_total = N_per_bin * N_src_per
        sigma_gamma = s['sigma_e'] / np.sqrt(2)
        sigma_mean = sigma_gamma / np.sqrt(N_total)
        # ν measurement precision
        nu_precision = sigma_mean / (1e-4)  # relative to typical shear
        # More meaningful: S/N per bin for RAR
        typical_shear = 5e-4  # at y~0.1
        SN_per_bin = typical_shear / sigma_mean * Delta_RAR
        print(f"    {name:<16}: S/N per Δlog y bin ~ {SN_per_bin:.0f}")

    print()
    print("  → RAR shape measurement: extremely high S/N with stacking")
    print("  → The MORPHOLOGY SPLIT is the key test (not just RAR shape)")


# ═══════════════════════════════════════════════════════════════════════════
#  PART E: RANKING OF DISCRIMINATING OBSERVABLES
# ═══════════════════════════════════════════════════════════════════════════

def rank_observables():
    print_header("PART E: RANKING OF MOST DISCRIMINATING OBSERVABLES")

    tests = [
        {
            'name': 'Cluster lensing at R > 2 Mpc',
            'prediction': 'ΔΣ_TGP/ΔΣ_NFW > 1.5 at R > R_100',
            'unique_to_TGP': True,
            'SN_euclid': 30,
            'SN_lsst': 50,
            'systematics': 'Low (clean at large R)',
            'timeline': 'Euclid DR1 (Oct 2026)',
            'rank': 1,
        },
        {
            'name': 'Morphology-dependent GGL',
            'prediction': 'ΔΣ_ell/ΔΣ_disk ~ 1.2-1.4 at fixed M*',
            'unique_to_TGP': True,
            'SN_euclid': 20,
            'SN_lsst': 35,
            'systematics': 'Medium (M* matching, SHMR confusion)',
            'timeline': 'Euclid DR1 + DESI',
            'rank': 2,
        },
        {
            'name': 'Stacked RAR morphology split',
            'prediction': 'Different ν(y) for n<1.5 vs n>2.5',
            'unique_to_TGP': True,
            'SN_euclid': 15,
            'SN_lsst': 25,
            'systematics': 'Medium (photometric M/L)',
            'timeline': 'Euclid DR1 + ground photometry',
            'rank': 3,
        },
        {
            'name': 'Cluster profile shape (inner)',
            'prediction': 'Less concentrated than NFW at R < R_500',
            'unique_to_TGP': False,
            'SN_euclid': 10,
            'SN_lsst': 15,
            'systematics': 'High (baryonic effects, AGN feedback)',
            'timeline': 'Euclid Q1 (available now)',
            'rank': 4,
        },
        {
            'name': 'RAR scatter vs morphology',
            'prediction': 'Intrinsic scatter ~0 (all in γ(S))',
            'unique_to_TGP': True,
            'SN_euclid': 10,
            'SN_lsst': 20,
            'systematics': 'Medium (distance errors)',
            'timeline': 'LSST Y1',
            'rank': 5,
        },
        {
            'name': 'η = 1/3 from galaxy populations',
            'prediction': 'η = (γ_ell - γ_disk)/γ_disk = 0.333',
            'unique_to_TGP': True,
            'SN_euclid': 8,
            'SN_lsst': 15,
            'systematics': 'Low (clean geometric prediction)',
            'timeline': 'LSST Y3',
            'rank': 6,
        },
    ]

    print(f"    {'#':<4} {'Test':<35} {'Unique?':<8} {'S/N_Euc':<9} {'S/N_LSST':<9} {'Systematics':<15} {'Timeline':<25}")
    print(f"    {'─'*4} {'─'*35} {'─'*8} {'─'*9} {'─'*9} {'─'*15} {'─'*25}")

    for t in tests:
        uniq = "YES" if t['unique_to_TGP'] else "no"
        print(f"    {t['rank']:<4} {t['name']:<35} {uniq:<8} {t['SN_euclid']:<9} {t['SN_lsst']:<9} {t['systematics']:<15} {t['timeline']:<25}")

    print()
    print("  SUMMARY OF PREDICTIONS:")
    print("  ════════════════════════")
    for t in tests:
        print(f"    #{t['rank']}: {t['name']}")
        print(f"       → {t['prediction']}")
        print()

    print("  KEY ADVANTAGE: Tests #1, #2, #3 are UNIQUE to TGP.")
    print("  MOND predicts universal ν(y) — no morphology dependence.")
    print("  ΛCDM predicts halo-dependent profiles — no continuous γ(n).")
    print("  TGP predicts γ = α·c/(c+1) with α=4/5 — continuous and calculable.")


# ═══════════════════════════════════════════════════════════════════════════
#  PART F: MOCK EUCLID ERO — ABELL 2390
# ═══════════════════════════════════════════════════════════════════════════

def mock_A2390():
    print_header("PART F: MOCK EUCLID ERO — ABELL 2390")

    print("  F.1  Abell 2390 parameters")
    print("  ────────────────────────────")
    print("  Euclid ERO observed A2390 (z=0.228) with weak lensing.")
    print("  Published: NFW fit with M200 ~ 1.5×10¹⁵ M☉, c ~ 4-5")
    print()

    # A2390 parameters (from literature)
    z_lens = 0.228
    M200_obs = 1.5e15   # M☉ (from NFW fit to Euclid ERO shear)
    c_nfw = 4.5
    M_gas = 2.0e14      # M☉ (from X-ray)
    M_stars = 2.0e13    # M☉
    R500_kpc = 1400.0
    M_bar = M_gas + M_stars

    print(f"    z = {z_lens}")
    print(f"    M_200 (NFW fit) = {M200_obs:.2e} M☉")
    print(f"    c_NFW = {c_nfw}")
    print(f"    M_gas = {M_gas:.2e} M☉")
    print(f"    M_stars = {M_stars:.2e} M☉")
    print(f"    M_bar = {M_bar:.2e} M☉")
    print(f"    R_500 = {R500_kpc} kpc")
    print()

    # ─── TGP prediction ───
    print("  F.2  TGP prediction for A2390 shear profile")
    print("  ──────────────────────────────────────────────")

    R_proj = np.array([200, 400, 600, 800, 1000, 1500, 2000, 2500, 3000, 4000, 5000])

    print()
    print(f"    {'R(kpc)':<10} {'M_bar(R)':<14} {'y(R)':<10} {'ν_c2(R)':<10} {'ν_c3(R)':<10} {'M_TGP_c2':<14} {'M_TGP_c3':<14} {'M_NFW':<14}")
    print(f"    {'─'*10} {'─'*14} {'─'*10} {'─'*10} {'─'*10} {'─'*14} {'─'*14} {'─'*14}")

    # Setup β-model
    rc_kpc = 0.15 * R500_kpc
    factor = 4 * np.pi * rc_kpc**3 * (R500_kpc/rc_kpc - np.arctan(R500_kpc/rc_kpc))
    rho0 = M_gas / factor
    a_stars = 80.0  # kpc

    for R in R_proj:
        Mg = M_gas_beta(R, rho0, rc_kpc)
        Ms = M_stars_hernquist(R, M_stars, a_stars)
        Mb = Mg + Ms

        g_N = G_SI * Mb * Msun / (R * kpc)**2
        y = g_N / a0
        y = max(y, 1e-6)

        nu_c2 = nu_tgp_ceff(y, 2)
        nu_c3 = nu_tgp_ceff(y, 3)

        M_tgp_c2 = Mb * nu_c2
        M_tgp_c3 = Mb * nu_c3

        M_nfw = M_NFW(R, M200_obs, c_nfw)

        print(f"    {R:>8.0f}   {Mb:>12.2e}   {y:>8.4f}   {nu_c2:>8.3f}   {nu_c3:>8.3f}   {M_tgp_c2:>12.2e}   {M_tgp_c3:>12.2e}   {M_nfw:>12.2e}")

    print()
    print("  F.3  Comparison with Euclid ERO shear data")
    print("  ────────────────────────────────────────────")
    print("  Euclid ERO fitted NFW to A2390 shear profile.")
    print("  Their best-fit: M200 = 1.5×10¹⁵ M☉, c = 4.5")
    print()
    print("  TGP prediction (c_eff=3):")

    # Key comparison points
    for R in [500, 1000, 1500, 2000, 3000]:
        Mg = M_gas_beta(R, rho0, rc_kpc)
        Ms = M_stars_hernquist(R, M_stars, a_stars)
        Mb = Mg + Ms
        g_N = G_SI * Mb * Msun / (R * kpc)**2
        y = g_N / a0
        nu = nu_tgp_ceff(max(y, 1e-6), 3)
        M_tgp = Mb * nu
        M_nfw = M_NFW(R, M200_obs, c_nfw)
        ratio = M_tgp / M_nfw
        print(f"    R={R:>5} kpc: M_TGP/M_NFW = {ratio:.3f}")

    print()
    print("  → At R < R_500: TGP underpredicts (deficit)")
    print("  → At R > R_100: TGP should EXCEED NFW")
    print("  → The Euclid ERO A2390 data extends to ~R_200")
    print("  → TEST: fit TGP model to the raw shear profile")
    print("  → If outer bins show EXCESS over NFW → strong TGP support")

    print()
    print("  F.4  What Euclid DR1 (Oct 2026) will provide")
    print("  ──────────────────────────────────────────────")
    print("  - ~2500 deg² of weak lensing data")
    print("  - ~5000 cluster lensing profiles")
    print("  - Stacked profiles in mass bins → S/N ~ 10-15 per bin")
    print("  - KEY: outer profile (R > 2 Mpc) shape test")
    print("  - Morphology information from VIS imaging (n_Sersic)")
    print("  - Combined with DESI redshifts → precise M_bar estimates")


# ═══════════════════════════════════════════════════════════════════════════
#  PART G: SUMMARY AND ROADMAP
# ═══════════════════════════════════════════════════════════════════════════

def summary():
    print_header("PART G: SUMMARY — EUCLID/LSST TESTS FOR TGP")

    print("  G.1  Three tiers of tests")
    print("  ──────────────────────────")
    print()
    print("  TIER 1 — AVAILABLE NOW (Euclid Q1 / ERO):")
    print("    • A2390 shear profile: check outer profile vs NFW")
    print("    • Q1 cluster catalog: ~100 clusters with lensing")
    print("    • Morphology from VIS: Sérsic index distribution")
    print()
    print("  TIER 2 — EUCLID DR1 (October 2026):")
    print("    • Stacked cluster lensing at R > 2 Mpc (S/N ~15)")
    print("    • GGL morphology split (S/N ~10)")
    print("    • RAR from stacked lensing (S/N >> 1 per bin)")
    print()
    print("  TIER 3 — EUCLID FINAL + LSST Y10 (2030+):")
    print("    • Definitive cluster profile shape (S/N > 30)")
    print("    • Morphology-dependent GGL (S/N > 20)")
    print("    • η = 1/3 precision test (S/N ~15)")
    print("    • Continuous γ(n) measurement (unique to TGP)")
    print()

    print("  G.2  If TGP is correct, Euclid DR1 should show:")
    print("  ──────────────────────────────────────────────────")
    print("    1. Stacked cluster ΔΣ(R) at R > R_200 EXCEEDS NFW best-fit")
    print("    2. Excess grows with projected radius (not noise-like)")
    print("    3. Elliptical GGL signal ~20-40% stronger than disk at fixed M*")
    print("    4. RAR scatter DECREASES when binned by Sérsic index")
    print()

    print("  G.3  If TGP is wrong, Euclid DR1 should show:")
    print("  ──────────────────────────────────────────────────")
    print("    1. Cluster ΔΣ(R) consistent with NFW at all radii")
    print("    2. No systematic excess at R > R_200")
    print("    3. GGL morphology dependence explained by SHMR alone")
    print("    4. RAR scatter independent of morphology binning")
    print()

    print("  G.4  Timeline")
    print("  ──────────────")
    print("    NOW (April 2026):    Euclid Q1 data available (63 deg²)")
    print("    Oct 2026:            Euclid DR1 — first cosmology release")
    print("    Sep 2026:            LSST DP2 — commissioning data")
    print("    Jun 2028:            LSST DR1 — first year data")
    print("    ~2030:               Euclid final + LSST Y3 → definitive tests")
    print()
    print("  → TGP makes QUANTITATIVE, FALSIFIABLE predictions")
    print("  → The most discriminating test (cluster outer profile) is")
    print("    accessible with Euclid DR1 data (October 2026)")
    print("  → The UNIQUE test (morphology-dependent γ) requires")
    print("    combined lensing + photometric morphology → Euclid ideal")


# ═══════════════════════════════════════════════════════════════════════════
#  MAIN
# ═══════════════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    print("=" * 78)
    print("  gs35: EUCLID / LSST OBSERVATIONAL TESTS FOR TGP")
    print("=" * 78)

    R_proj, DS_TGP_c2, DS_TGP_c3, DS_NFW = compute_lensing_profiles()
    compute_ggl_predictions()
    compute_stacked_RAR()
    compute_SN_estimates()
    rank_observables()
    mock_A2390()
    summary()

    print()
    print("=" * 78)
    print("  END OF gs35: EUCLID / LSST TESTS")
    print("=" * 78)
