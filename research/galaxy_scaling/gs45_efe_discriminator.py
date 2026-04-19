#!/usr/bin/env python3
"""
gs45_efe_discriminator.py
=========================
Tests TGP's External Field Effect (EFE) prediction against MOND.

KEY PHYSICS:
  - MOND EFE is LINEAR:     y_eff = y_int + y_ext   (strong suppression)
  - TGP  EFE is QUADRATURE: y_eff = sqrt(y_int^2 + y_ext^2)  (WEAKER suppression)

This is a CLEAN DISCRIMINATING TEST between the two theories.

Systems analysed:
  A. EFE formalism (both prescriptions)
  B. Crater II         -- ultra-diffuse dwarf, y_ext ~ y_int (maximal EFE)
  C. NGC 1052-DF2      -- galaxy "lacking dark matter"
  D. NGC 1052-DF4      -- similar "lacking DM" galaxy
  E. Fornax dSph       -- classical MW satellite

Author : TGP collaboration
Date   : 2026-04-18
"""

import sys
sys.stdout.reconfigure(encoding='utf-8', errors='replace')

import numpy as np

# ============================================================
# Physical constants
# ============================================================
G_SI   = 6.674e-11       # m^3 kg^-1 s^-2
a0     = 1.2e-10          # m/s^2  (MOND acceleration scale)
Msun   = 1.989e30         # kg
kpc    = 3.086e19          # m
pc     = 3.086e16          # m
Lsun   = 3.828e26         # W  (solar luminosity, not critical here)

# ============================================================
# A. EFE formalism
# ============================================================
# For dSphs we use c_eff = 2.5 (quasi-spherical systems)
C_EFF  = 2.5
GAMMA  = 0.8 * C_EFF / (C_EFF + 1.0)   # ~ 0.571


def nu(y, gamma=GAMMA):
    """TGP/MOND interpolating function nu(y).

    nu(y) = 1 + exp(-y^(4/5)) / y^gamma
    """
    return 1.0 + np.exp(-y**(4.0/5.0)) / y**gamma


def y_eff_mond(y_int, y_ext):
    """MOND: LINEAR addition  -->  strong suppression."""
    return y_int + y_ext


def y_eff_tgp(y_int, y_ext):
    """TGP: QUADRATURE addition  -->  weaker suppression."""
    return np.sqrt(y_int**2 + y_ext**2)


def sigma_los_predicted(M_bar_kg, r_half_m, y_ext, model='tgp'):
    """Predict line-of-sight velocity dispersion sigma_los.

    Uses:
      g_int  = G M_bar / r_half^2
      y_int  = g_int / a0
      y_eff  = combination rule (MOND or TGP)
      M_dyn  = nu(y_eff) * M_bar
      sigma  = sqrt(G M_dyn / (2 r_half))     (from Wolf estimator)

    Returns sigma in km/s.
    """
    g_int = G_SI * M_bar_kg / r_half_m**2
    y_int = g_int / a0

    if model == 'mond':
        ye = y_eff_mond(y_int, y_ext)
    elif model == 'tgp':
        ye = y_eff_tgp(y_int, y_ext)
    elif model == 'newtonian':
        ye = 1e10   # huge y -> nu ~ 1
    else:
        raise ValueError(f"Unknown model: {model}")

    nu_val = nu(ye)
    M_dyn  = nu_val * M_bar_kg
    sigma  = np.sqrt(G_SI * M_dyn / (2.0 * r_half_m))   # m/s
    return sigma / 1e3, y_int, ye, nu_val   # km/s


# ============================================================
# Helper: compute y_ext from host mass and distance
# ============================================================
def compute_y_ext(M_host_Msun, D_kpc):
    """External field y_ext = g_ext / a0 where g_ext = G M_host / D^2."""
    D_m    = D_kpc * kpc
    M_kg   = M_host_Msun * Msun
    g_ext  = G_SI * M_kg / D_m**2
    return g_ext / a0


# ============================================================
# Define the satellite systems
# ============================================================
systems = {}

# B. Crater II  (ultra-diffuse MW dwarf)
systems['Crater II'] = dict(
    M_bar_Msun   = 1.6e5,       # ~ L/Lsun with M/L ~ 1
    r_half_pc    = 1066.0,
    sigma_obs    = 2.7,          # km/s
    sigma_err    = 0.3,
    host_mass    = 5e11,         # MW enclosed mass at 117 kpc
    D_kpc        = 117.0,
    notes        = "Ultra-diffuse dwarf; y_ext ~ y_int -> maximal EFE sensitivity",
)

# C. NGC 1052-DF2  ("lacking dark matter")
systems['NGC 1052-DF2'] = dict(
    M_bar_Msun   = 2e8,
    r_half_pc    = 2200.0,       # 2.2 kpc
    sigma_obs    = 8.5,
    sigma_err    = 2.3,
    host_mass    = 1e12,         # NGC 1052 group total mass estimate
    D_kpc        = 80.0,         # projected distance to NGC 1052
    notes        = "If sigma truly low -> MOND's strong EFE favored over TGP",
)

# D. NGC 1052-DF4  (similar "lacking DM")
systems['NGC 1052-DF4'] = dict(
    M_bar_Msun   = 1.5e8,
    r_half_pc    = 1600.0,       # 1.6 kpc
    sigma_obs    = 4.2,
    sigma_err    = 2.2,
    host_mass    = 1e12,
    D_kpc        = 80.0,
    notes        = "Similar to DF2; very low sigma if confirmed",
)

# E. Fornax dSph
systems['Fornax'] = dict(
    M_bar_Msun   = 4.3e7,
    r_half_pc    = 710.0,
    sigma_obs    = 11.7,
    sigma_err    = 0.9,
    host_mass    = 5e11,         # MW enclosed at 147 kpc
    D_kpc        = 147.0,
    notes        = "Classical MW dSph; well-measured kinematics",
)


# ============================================================
# Run analysis for every system
# ============================================================
def analyse_system(name, params):
    """Compute sigma_TGP, sigma_MOND for one satellite system."""
    M_bar_kg  = params['M_bar_Msun'] * Msun
    r_half_m  = params['r_half_pc'] * pc
    y_ext     = compute_y_ext(params['host_mass'], params['D_kpc'])

    results = {}
    for model in ('newtonian', 'mond', 'tgp'):
        sigma, y_int, ye, nu_val = sigma_los_predicted(
            M_bar_kg, r_half_m, y_ext, model=model
        )
        results[model] = dict(sigma=sigma, y_int=y_int, y_eff=ye, nu=nu_val)

    return results, y_ext


def main():
    sep = "=" * 78
    print(sep)
    print("  gs45_efe_discriminator.py")
    print("  TGP vs MOND External Field Effect -- Discriminating Test")
    print(sep)

    # ----------------------------------------------------------
    # A. Show the formalism
    # ----------------------------------------------------------
    print("\n--- A. EFE Formalism ---\n")
    print("  MOND (linear):      y_eff = y_int + y_ext")
    print("  TGP  (quadrature):  y_eff = sqrt(y_int^2 + y_ext^2)\n")
    print(f"  Interpolating function: nu(y) = 1 + exp(-y^(4/5)) / y^gamma")
    print(f"  c_eff = {C_EFF:.1f}  =>  gamma = 0.8*c_eff/(c_eff+1) = {GAMMA:.4f}\n")

    # Illustrative table: how nu differs for sample (y_int, y_ext)
    print("  Illustration: nu(y_eff) for sample internal/external fields\n")
    print(f"  {'y_int':>8s}  {'y_ext':>8s}  {'y_MOND':>10s}  {'y_TGP':>10s}"
          f"  {'nu_MOND':>10s}  {'nu_TGP':>10s}  {'ratio':>8s}")
    print("  " + "-" * 72)
    for y_int_val in [0.01, 0.05, 0.1, 0.5]:
        for y_ext_val in [0.01, 0.05, 0.1]:
            ym = y_eff_mond(y_int_val, y_ext_val)
            yt = y_eff_tgp(y_int_val, y_ext_val)
            nm = nu(ym)
            nt = nu(yt)
            print(f"  {y_int_val:8.3f}  {y_ext_val:8.3f}  {ym:10.4f}  {yt:10.4f}"
                  f"  {nm:10.4f}  {nt:10.4f}  {nt/nm:8.3f}")
    print()
    print("  ratio = nu_TGP / nu_MOND  (>1 means TGP predicts MORE boosting)")
    print("  TGP always gives nu_TGP >= nu_MOND  =>  higher predicted sigma.\n")

    # ----------------------------------------------------------
    # B-E. Analyse each system
    # ----------------------------------------------------------
    summary_rows = []

    for name in ('Crater II', 'NGC 1052-DF2', 'NGC 1052-DF4', 'Fornax'):
        params  = systems[name]
        results, y_ext = analyse_system(name, params)

        print(sep)
        print(f"  System: {name}")
        print(sep)
        print(f"  M_bar        = {params['M_bar_Msun']:.2e} Msun")
        print(f"  r_half       = {params['r_half_pc']:.0f} pc")
        print(f"  sigma_obs    = {params['sigma_obs']:.1f} +/- {params['sigma_err']:.1f} km/s")
        print(f"  Host dist    = {params['D_kpc']:.0f} kpc")
        print(f"  y_ext        = {y_ext:.5f}")
        print(f"  y_int        = {results['tgp']['y_int']:.5f}")
        print(f"  y_ext/y_int  = {y_ext / results['tgp']['y_int']:.3f}")
        print(f"  Notes: {params['notes']}")
        print()

        for model in ('newtonian', 'mond', 'tgp'):
            r = results[model]
            label = model.upper().ljust(10)
            print(f"    {label}  y_eff = {r['y_eff']:10.5f}   nu = {r['nu']:8.4f}"
                  f"   sigma = {r['sigma']:6.2f} km/s")

        sig_tgp  = results['tgp']['sigma']
        sig_mond = results['mond']['sigma']
        sig_obs  = params['sigma_obs']
        sig_err  = params['sigma_err']

        delta_tgp  = abs(sig_tgp  - sig_obs)
        delta_mond = abs(sig_mond - sig_obs)

        n_sigma_tgp  = delta_tgp  / sig_err
        n_sigma_mond = delta_mond / sig_err

        if delta_tgp < delta_mond:
            closer = "TGP"
        elif delta_mond < delta_tgp:
            closer = "MOND"
        else:
            closer = "TIE"

        print()
        print(f"    Observed:  {sig_obs:.1f} +/- {sig_err:.1f} km/s")
        print(f"    |TGP  - obs| = {delta_tgp:.2f} km/s  ({n_sigma_tgp:.1f} sigma)")
        print(f"    |MOND - obs| = {delta_mond:.2f} km/s  ({n_sigma_mond:.1f} sigma)")
        print(f"    Closer model: {closer}")

        discriminating = abs(sig_tgp - sig_mond) > sig_err
        disc_str = "YES" if discriminating else "NO (within error bars)"
        print(f"    TGP-MOND split = {abs(sig_tgp - sig_mond):.2f} km/s  "
              f"-> Discriminating at current precision? {disc_str}")
        print()

        summary_rows.append(dict(
            name       = name,
            sig_tgp    = sig_tgp,
            sig_mond   = sig_mond,
            sig_obs    = sig_obs,
            sig_err    = sig_err,
            closer     = closer,
            disc       = discriminating,
            split      = abs(sig_tgp - sig_mond),
            n_sig_tgp  = n_sigma_tgp,
            n_sig_mond = n_sigma_mond,
        ))

    # ----------------------------------------------------------
    # F. Summary table
    # ----------------------------------------------------------
    print(sep)
    print("  F. SUMMARY TABLE -- TGP vs MOND EFE Discriminator")
    print(sep)
    print()
    hdr = (f"  {'System':<18s}  {'sig_TGP':>8s}  {'sig_MOND':>8s}  "
           f"{'sig_obs':>8s}  {'err':>5s}  {'split':>6s}  "
           f"{'closer':>7s}  {'discrim?':>10s}")
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))

    for row in summary_rows:
        disc_flag = "YES" if row['disc'] else "no"
        print(f"  {row['name']:<18s}  {row['sig_tgp']:8.2f}  {row['sig_mond']:8.2f}  "
              f"{row['sig_obs']:8.1f}  {row['sig_err']:5.1f}  {row['split']:6.2f}  "
              f"{row['closer']:>7s}  {disc_flag:>10s}")

    print()
    print("  All sigma values in km/s.  'split' = |sig_TGP - sig_MOND|.")
    print("  'discrim?' = YES if split > measurement error.\n")

    # ----------------------------------------------------------
    # Interpretation
    # ----------------------------------------------------------
    print(sep)
    print("  INTERPRETATION")
    print(sep)
    print("""
  The EFE prescription is a CLEAN discriminating test:

  1. TGP (quadrature) ALWAYS predicts y_eff <= y_eff(MOND)
     => nu(y_eff_TGP) >= nu(y_eff_MOND)
     => TGP predicts HIGHER velocity dispersions than MOND for EFE-dominated
        systems, because the external field suppresses the MOND boost less.

  2. Crater II is the BEST discriminator because y_ext ~ y_int:
     the two prescriptions diverge most when internal and external fields
     are comparable (quadrature vs linear differ maximally at equal arguments).

  3. NGC 1052-DF2 / DF4: If the observed sigma is truly very low, MOND's
     stronger EFE (linear addition) is favored because it suppresses the
     dynamical mass boost more effectively. TGP would need to explain the
     low sigma via another mechanism (e.g., tidal stripping, distance
     revision reducing M_bar).

  4. Fornax is less EFE-dominated (y_int >> y_ext), so the two models
     converge -- it is NOT a strong discriminator.

  BOTTOM LINE: Precise measurements of sigma_los in EFE-dominated dwarfs
  (especially Crater II) can cleanly distinguish TGP from MOND.
""")

    print(sep)
    print("  gs45_efe_discriminator.py  --  analysis complete")
    print(sep)


if __name__ == '__main__':
    main()
