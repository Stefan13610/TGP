#!/usr/bin/env python3
"""
gs61_sparc_membrane_fit.py
==========================

Pelny fit SPARC (175 galaktyk) modelem membrane z gs9d
+ a0(z=0) = c*H0/(2*pi) z gs46/gs60.

MODEL MEMBRANE (gs9d, geometric mean mechanism):

  H(M) = sqrt(r_S * r_H) = sqrt(G*M / (c*H0))

  Model rotation:
    v_tot^2(r) = v_bar^2(r) + sqrt(kappa * G * M_bar(r) * c * H0)
             = v_bar^2(r) * (1 + r/H(r))  with M(r) from V_bar^2

  kappa = 1 : pure membrane model from gs9d (predicts v^4 = G*M*c*H0)
  kappa = (2*pi*a0_obs)/(c*H0) ~ 0.183 : renormalized to observed MOND

  If kappa=1 holds naturally -> BTFR slope = 4, normalization a0_eff = c*H0
  If observed MOND emerges -> slope = 4 but kappa = 0.18 needed

  a0(z=0) = c*H0/(2*pi) = 1.04e-10  (TGP prediction from gs60)

TESTS:
======
1. Per-galaxy chi2_red for 5 models:
   - Newton (baseline, no modification)
   - MOND simple (nu = 0.5 + sqrt(0.25 + 1/y))
   - TGP standard (nu = 1 + exp(-y^0.8)/y^0.5714)  [gs57 falsified]
   - Form 4 (nu = 1/(1-exp(-y^0.497)))  [gs56 best-fit]
   - MEMBRANE (gs9d + gs60, parameter-free in kappa=1 form)

2. BTFR slope: log(M_bar) vs log(V_flat)
   - MOND predicts: log(M) = 4*log(V) - log(G*a0)
   - Membrane predicts: log(M) = 4*log(V) - log(G*c*H0)  (shifted by 2*pi)
   - Observed: slope ~ 4, intercept consistent with a0 = 1.2e-10

3. Inferred kappa: fit kappa globally so membrane matches SPARC
   - If kappa ~ 1: gs9d prediction correct
   - If kappa ~ 0.18: renormalization needed; 2*pi factor is wrong

4. Mass-scale trends: systematic residuals by M_bar category

OUTPUT: full table + BTFR fit + kappa estimate + falsification criteria.
"""

import numpy as np
import os
import sys
import io
import warnings

warnings.filterwarnings('ignore')
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

# =============================================================================
# Constants
# =============================================================================
G     = 6.67430e-11          # m^3/(kg*s^2)
c     = 2.99792458e8         # m/s
H0    = 2.1844e-18           # s^-1 (67.4 km/s/Mpc)
a0_obs = 1.2e-10             # m/s^2, observed MOND scale
a0_TGP = c * H0 / (2*np.pi)  # m/s^2, TGP prediction: 1.04e-10

kpc_to_m    = 3.0857e19
km_to_m     = 1.0e3
Msun        = 1.989e30        # kg
Lsun        = 3.828e26        # W (irrelevant here, just for reference)

ML_disk_default = 0.5        # M/L for disk at 3.6 micron (SPARC convention)
ML_bul_default  = 0.7        # M/L for bulge at 3.6 micron
M_HI_to_bar     = 1.33       # He correction

# =============================================================================
# Interpolation functions (comparison baselines)
# =============================================================================

def nu_newton(y):
    return np.ones_like(y)

def nu_mond(y):
    y = np.maximum(y, 1e-30)
    return 0.5 + np.sqrt(0.25 + 1.0/y)

def nu_tgp(y):
    """TGP standard: nu(y) = 1 + exp(-y^0.8) / y^0.5714"""
    y = np.maximum(y, 1e-30)
    return 1.0 + np.exp(-y**0.8) / y**0.5714

def nu_form4(y):
    """Form 4 best-fit: nu(y) = 1/(1 - exp(-y^0.497))"""
    y = np.maximum(y, 1e-30)
    val = 1.0 - np.exp(-y**0.497)
    val = np.maximum(val, 1e-30)
    return 1.0 / val

# =============================================================================
# MEMBRANE MODEL (gs9d geometric mean, parametrized by kappa)
# =============================================================================

def v_membrane_squared(Vbar2_kms2, R_kpc, kappa=1.0):
    """
    Compute v_total^2 in (km/s)^2 via membrane model.

    Derivation (gs9d + kappa):
      M_bar(r) = V_bar^2(r) * r / G  [enclosed baryonic mass]
      H(r) = sqrt(kappa * G * M_bar(r) / (c*H0))  [geometric-mean membrane depth]
      v_tot^2(r) = V_bar^2(r) * (1 + r/H(r))
                = V_bar^2(r) + V_bar^2(r)*r/H(r)
                = V_bar^2(r) + sqrt(kappa * G * M_bar(r) * c*H0)
                = V_bar^2(r) + sqrt(kappa * V_bar^2(r) * r * c*H0)

    For kappa=1: v_flat^2 -> sqrt(G*M*c*H0) = sqrt(2pi*G*M*a0_TGP)
    -> BTFR: v^4 = 2pi * G * M * a0_TGP = G * M * c*H0

    For kappa = a0_obs / (2pi*a0_TGP) = a0_obs/(c*H0) -> v^4 = G*M*a0_obs (MOND)
    """
    Vbar2_ms2 = Vbar2_kms2 * km_to_m**2    # (m/s)^2
    R_m       = R_kpc * kpc_to_m           # meters
    # additional term:  sqrt(kappa * V_bar^2 * r * c*H0)
    add_ms2   = np.sqrt(np.maximum(kappa * Vbar2_ms2 * R_m * c * H0, 0.0))
    vtot2_ms2 = Vbar2_ms2 + add_ms2        # (m/s)^2
    return vtot2_ms2 / km_to_m**2          # back to (km/s)^2

def v_nu(nu_func, Vbar2_kms2, R_kpc, a0=a0_obs):
    """Standard nu(y) formulation: v_tot^2 = nu(g_bar/a0) * V_bar^2."""
    Vbar2_ms2 = Vbar2_kms2 * km_to_m**2
    R_m       = R_kpc * kpc_to_m
    gbar      = Vbar2_ms2 / R_m            # m/s^2
    y         = gbar / a0
    return nu_func(y) * Vbar2_kms2

# =============================================================================
# SPARC Table I loader (baryonic masses, Vflat)
# =============================================================================

def load_sparc_table(path):
    """
    Parse SPARC_Lelli2016c.mrt (Lelli+2016 Table 1) using whitespace split.

    Column layout (19 fields):
      [0] Galaxy  [1] T  [2] D  [3] e_D  [4] f_D
      [5] Inc    [6] e_Inc  [7] L[3.6]  [8] e_L[3.6]
      [9] Reff   [10] SBeff  [11] Rdisk [12] SBdisk
      [13] MHI   [14] RHI   [15] Vflat  [16] e_Vflat
      [17] Q     [18] Ref
    Returns list of dicts.
    """
    galaxies = []
    with open(path, 'r', encoding='utf-8', errors='replace') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 18:
                continue
            # Galaxy name field -- first token (may contain dashes etc.)
            # Hubble type (parts[1]) must be integer
            try:
                T = int(parts[1])
                D = float(parts[2])
            except (ValueError, IndexError):
                continue
            try:
                L36    = float(parts[7])
                MHI    = float(parts[13])
                Vflat  = float(parts[15])
                eVflat = float(parts[16])
                Q      = int(parts[17])
            except (ValueError, IndexError):
                continue
            galaxies.append({
                'name': parts[0],
                'L36': L36, 'MHI': MHI,
                'Vflat': Vflat, 'eVflat': eVflat, 'Q': Q,
            })
    return galaxies

# =============================================================================
# Rotmod loader (per-galaxy rotation curves)
# =============================================================================

def load_rotmod(filepath):
    """Load SPARC rotmod. Returns arrays (R_kpc, Vobs, errV, Vgas, Vdisk, Vbul)."""
    R, Vobs, errV, Vgas, Vdisk, Vbul = [], [], [], [], [], []
    with open(filepath, 'r') as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith('#'):
                continue
            parts = s.split()
            if len(parts) < 6:
                continue
            try:
                R.append(float(parts[0]))
                Vobs.append(float(parts[1]))
                errV.append(float(parts[2]))
                Vgas.append(float(parts[3]))
                Vdisk.append(float(parts[4]))
                Vbul.append(float(parts[5]))
            except ValueError:
                continue
    return (np.array(R), np.array(Vobs), np.array(errV),
            np.array(Vgas), np.array(Vdisk), np.array(Vbul))

def vbar2(Vgas, Vdisk, Vbul, ml_disk=ML_disk_default, ml_bul=ML_bul_default):
    return Vgas**2 + ml_disk*Vdisk**2 + ml_bul*Vbul**2

def chi2_reduced(Vobs, Vpred, errV, dof_offset=0):
    err = np.maximum(errV, 1.0)   # avoid div-by-zero; 1 km/s floor
    chi2 = np.sum((Vobs - Vpred)**2 / err**2)
    dof = max(len(Vobs) - dof_offset, 1)
    return chi2 / dof

# =============================================================================
# MAIN
# =============================================================================

def main():
    print("=" * 78)
    print("  gs61_sparc_membrane_fit.py")
    print("  Full SPARC 175-galaxy fit with gs9d membrane + gs60 a0 = cH0/(2pi)")
    print("=" * 78)
    print()
    print(f"Constants:")
    print(f"  H0       = {H0:.4e} s^-1  (67.4 km/s/Mpc)")
    print(f"  a0_obs   = {a0_obs:.3e} m/s^2  (MOND empirical)")
    print(f"  a0_TGP   = c*H0/(2pi) = {a0_TGP:.3e} m/s^2  (TGP prediction gs60)")
    print(f"  c*H0     = {c*H0:.3e} m/s^2  (membrane v^4 normalization)")
    print(f"  ratio    = c*H0 / a0_obs = {(c*H0)/a0_obs:.3f}")
    print(f"  kappa_MOND = a0_obs/(c*H0) = {a0_obs/(c*H0):.4f}  (needed if gs9d is off by 2pi^2)")
    print()

    # Path setup
    here = os.path.dirname(os.path.abspath(__file__))
    data_dir   = os.path.join(here, 'Rotmod_LTG')
    table_path = os.path.join(here, 'SPARC_Lelli2016c.mrt')

    # ---------------------------------------------------------------
    # Load Table I
    # ---------------------------------------------------------------
    print("=" * 78)
    print("  Loading SPARC Table I (Lelli+2016)")
    print("=" * 78)
    table = load_sparc_table(table_path)
    print(f"  Loaded {len(table)} galaxies from Table I")
    # Quality cut: use Q=1 for BTFR (flat velocity reliable)
    table_q1 = [g for g in table if g['Q'] == 1 and g['Vflat'] > 0]
    print(f"  Q=1 (flat RC reliable) with Vflat>0: {len(table_q1)} galaxies")
    print()

    # Build baryonic masses dict {name: {M_bar, Vflat, eVflat}}
    gal_info = {}
    for g in table:
        Mstar_9 = ML_disk_default * g['L36']              # 1e9 Msun (disk-only proxy)
        Mgas_9  = M_HI_to_bar * g['MHI']                  # 1e9 Msun
        Mbar_9  = Mstar_9 + Mgas_9
        gal_info[g['name']] = {
            'L36': g['L36'], 'MHI': g['MHI'],
            'Mbar_9': Mbar_9, 'Vflat': g['Vflat'],
            'eVflat': g['eVflat'], 'Q': g['Q'],
        }

    # ---------------------------------------------------------------
    # Scan Rotmod_LTG folder
    # ---------------------------------------------------------------
    print("=" * 78)
    print("  Loading rotmod files from Rotmod_LTG/")
    print("=" * 78)
    all_files = [f for f in os.listdir(data_dir) if f.endswith('_rotmod.dat')]
    print(f"  Found {len(all_files)} rotmod files")
    print()

    per_galaxy = {}   # name -> dict with R, Vobs, errV, Vbar2, Mbar_9
    for fname in all_files:
        name = fname.replace('_rotmod.dat', '')
        # Some naming mismatches: try variations
        alt_names = [name, name.replace('-', ''), name.replace(' ', '')]
        tbl_name = None
        for an in alt_names:
            if an in gal_info:
                tbl_name = an
                break
        if tbl_name is None:
            # try exact match in table by strip
            for tn in gal_info:
                if tn.replace(' ', '').replace('-', '') == name:
                    tbl_name = tn
                    break
        if tbl_name is None:
            continue
        fpath = os.path.join(data_dir, fname)
        R, Vobs, errV, Vg, Vd, Vb = load_rotmod(fpath)
        if len(R) < 3:
            continue
        Vb2 = vbar2(Vg, Vd, Vb)
        per_galaxy[tbl_name] = {
            'R': R, 'Vobs': Vobs, 'errV': errV,
            'Vgas': Vg, 'Vdisk': Vd, 'Vbul': Vb,
            'Vbar2': Vb2,
            'info': gal_info[tbl_name],
        }
    print(f"  Matched {len(per_galaxy)} galaxies with Table I entries")
    print()

    # ---------------------------------------------------------------
    # Per-galaxy chi2 for 5 models (fixed M/L)
    # ---------------------------------------------------------------
    print("=" * 78)
    print("  Per-galaxy chi2_red (fixed M/L_disk=0.5, M/L_bul=0.7)")
    print("=" * 78)
    print()

    results = {}    # name -> dict of chi2r per model
    for name, d in per_galaxy.items():
        R, Vobs, errV, Vb2 = d['R'], d['Vobs'], d['errV'], d['Vbar2']

        # Model predictions
        Vp_newton   = np.sqrt(np.maximum(Vb2, 0))
        Vp_mond     = np.sqrt(np.maximum(v_nu(nu_mond, Vb2, R, a0_obs), 0))
        Vp_tgp      = np.sqrt(np.maximum(v_nu(nu_tgp, Vb2, R, a0_obs), 0))
        Vp_form4    = np.sqrt(np.maximum(v_nu(nu_form4, Vb2, R, a0_obs), 0))
        Vp_mem_k1   = np.sqrt(np.maximum(v_membrane_squared(Vb2, R, kappa=1.0), 0))
        Vp_mem_kM   = np.sqrt(np.maximum(v_membrane_squared(Vb2, R, kappa=a0_obs/(c*H0)), 0))

        results[name] = {
            'Newton':        chi2_reduced(Vobs, Vp_newton, errV),
            'MOND':          chi2_reduced(Vobs, Vp_mond, errV),
            'TGP_std':       chi2_reduced(Vobs, Vp_tgp, errV),
            'Form4':         chi2_reduced(Vobs, Vp_form4, errV),
            'Memb_k=1':      chi2_reduced(Vobs, Vp_mem_k1, errV),
            'Memb_k=0.183':  chi2_reduced(Vobs, Vp_mem_kM, errV),
            'N':             len(R),
        }

    model_names = ['Newton', 'MOND', 'TGP_std', 'Form4', 'Memb_k=1', 'Memb_k=0.183']

    # Summary stats
    print("  Aggregate chi2_red over all galaxies:")
    print(f"  {'Model':15s}  {'mean':>8s}  {'median':>8s}  {'std':>8s}  {'wins':>5s}")
    print("  " + "-" * 52)
    for m in model_names:
        vals = [results[g][m] for g in results]
        wins = sum(1 for g in results
                   if min(model_names, key=lambda mm: results[g][mm]) == m)
        print(f"  {m:15s}  {np.mean(vals):8.2f}  {np.median(vals):8.2f}  "
              f"{np.std(vals):8.2f}  {wins:5d}")
    print()

    # Showcase galaxies
    showcase = ['DDO154', 'NGC2403', 'NGC3198', 'NGC7331', 'NGC2841', 'UGC02885']
    print("  Per-galaxy chi2_red (showcase):")
    hdr = f"  {'Galaxy':12s} {'N':>4s} "
    for m in model_names:
        hdr += f" {m:>13s}"
    print(hdr)
    print("  " + "-" * len(hdr))
    for g in showcase:
        if g not in results:
            continue
        r = results[g]
        line = f"  {g:12s} {r['N']:4d}  "
        for m in model_names:
            line += f"  {r[m]:12.2f}"
        print(line)
    print()

    # ---------------------------------------------------------------
    # BTFR test
    # ---------------------------------------------------------------
    print("=" * 78)
    print("  BTFR (Baryonic Tully-Fisher Relation)")
    print("  Empirical: log(M_bar) = 4 * log(V_flat) - log(G*a0)")
    print("=" * 78)
    print()

    # Collect Q=1 galaxies with matched data
    btfr_data = []
    for name, d in per_galaxy.items():
        info = d['info']
        if info['Q'] != 1 or info['Vflat'] <= 0:
            continue
        Mbar = info['Mbar_9'] * 1e9 * Msun      # kg
        Vflat_ms = info['Vflat'] * km_to_m       # m/s
        # From rotmod: compute v_flat predictions from membrane model
        R = d['R']
        Vb2 = d['Vbar2']
        # outer 30% of rotmod to define V_flat
        n_outer = max(3, len(R)//3)
        idx_outer = slice(-n_outer, None)

        Vp_mond_sq   = v_nu(nu_mond, Vb2, R, a0_obs)
        Vp_form4_sq  = v_nu(nu_form4, Vb2, R, a0_obs)
        Vp_memk1_sq  = v_membrane_squared(Vb2, R, kappa=1.0)
        Vp_memkM_sq  = v_membrane_squared(Vb2, R, kappa=a0_obs/(c*H0))

        V_mond_flat  = np.sqrt(np.mean(np.maximum(Vp_mond_sq [idx_outer], 0)))
        V_form4_flat = np.sqrt(np.mean(np.maximum(Vp_form4_sq[idx_outer], 0)))
        V_memk1_flat = np.sqrt(np.mean(np.maximum(Vp_memk1_sq[idx_outer], 0)))
        V_memkM_flat = np.sqrt(np.mean(np.maximum(Vp_memkM_sq[idx_outer], 0)))

        btfr_data.append({
            'name': name, 'Mbar': Mbar,
            'Vflat_obs': Vflat_ms,
            'V_mond':  V_mond_flat * km_to_m,
            'V_form4': V_form4_flat * km_to_m,
            'V_memk1': V_memk1_flat * km_to_m,
            'V_memkM': V_memkM_flat * km_to_m,
        })

    print(f"  BTFR sample size (Q=1 with Vflat>0, matched rotmod): {len(btfr_data)}")
    print()

    # Do linear regression log(M) = slope * log(V) + intercept
    def lin_fit(logV, logM):
        A = np.vstack([logV, np.ones_like(logV)]).T
        slope, intercept = np.linalg.lstsq(A, logM, rcond=None)[0]
        resid = logM - (slope*logV + intercept)
        sigma = np.std(resid, ddof=2)
        n = len(logV)
        # standard error on slope
        xmean = np.mean(logV)
        sxx = np.sum((logV - xmean)**2)
        se_slope = sigma / np.sqrt(sxx)
        return slope, intercept, sigma, se_slope

    logM = np.log10([b['Mbar'] / Msun for b in btfr_data])   # log10(Mbar / Msun)
    logV_obs = np.log10([b['Vflat_obs'] / km_to_m for b in btfr_data])  # log10(Vflat / km/s)

    print("  BTFR fit: log10(M_bar/Msun) = slope * log10(V_flat/[km/s]) + intercept")
    print()
    print(f"  {'Source':20s}  {'slope':>8s}  {'se_slope':>9s}  {'intercept':>10s}  "
          f"{'a0_eff':>13s}  {'scatter':>8s}")
    print("  " + "-" * 72)

    for label, Vkey in [
        ('Observed (SPARC)', 'Vflat_obs'),
        ('MOND simple',      'V_mond'),
        ('Form 4',           'V_form4'),
        ('Membrane k=1',     'V_memk1'),
        ('Membrane k=MOND',  'V_memkM'),
    ]:
        logV = np.log10([b[Vkey] / km_to_m for b in btfr_data])
        slope, intercept, sigma, se_slope = lin_fit(logV, logM)
        # a0_eff from BTFR: M = V^4 / (G*a0) -> log M = 4*log V - log(G*a0)
        # intercept_expected = -log(G*a0) = -log G + log(1/a0)
        # If slope is near 4, we can force slope=4 and extract a0_eff
        # from log10(M) - 4*log10(V) = -log10(G*a0)
        diff = logM - 4*logV
        mean_diff = np.mean(diff)
        # mean_diff = log10(M_Msun) - 4*log10(V_km/s)
        # M = V^4/(G*a0); but M in Msun, V in km/s:
        # M_kg = V_kg_m^4 / (G*a0) = (V*1e3)^4 / (G*a0) = 1e12*V^4 / (G*a0)
        # M_Msun = 1e12*V^4 / (G*a0*Msun)
        # log10(M_Msun) = 12 + 4*log10(V) - log10(G*a0*Msun)
        # mean_diff = 12 - log10(G*a0*Msun)
        # a0 = 10^(12 - mean_diff) / (G*Msun)
        log_Ga0_Msun = 12 - mean_diff
        Ga0Msun = 10**log_Ga0_Msun
        a0_eff = Ga0Msun / (G * Msun)
        print(f"  {label:20s}  {slope:8.3f}  +/-{se_slope:7.3f}  {intercept:10.3f}  "
              f"{a0_eff:13.3e}  {sigma:8.3f}")

    print()
    print("  Interpretation:")
    print("    * Observed slope should be ~ 4.0 (MOND BTFR)")
    print("    * Forced slope=4, a0_eff:")
    print(f"      - Observed SPARC: should give ~ 1.2e-10")
    print(f"      - Membrane k=1:   would give ~ c*H0 = {c*H0:.2e}")
    print(f"      - Membrane k=MOND: should give ~ 1.2e-10")
    print()

    # ---------------------------------------------------------------
    # Kappa inference: fit kappa globally to minimize total chi2_red
    # ---------------------------------------------------------------
    print("=" * 78)
    print("  Global kappa fit (membrane model)")
    print("=" * 78)
    print()

    kappa_scan = np.logspace(-2, 0.5, 41)  # 0.01 to 3.16
    total_chi2 = []
    for kappa in kappa_scan:
        t = 0.0
        n_tot = 0
        for name, d in per_galaxy.items():
            R, Vobs, errV, Vb2 = d['R'], d['Vobs'], d['errV'], d['Vbar2']
            Vp_sq = v_membrane_squared(Vb2, R, kappa=kappa)
            Vp = np.sqrt(np.maximum(Vp_sq, 0))
            err = np.maximum(errV, 1.0)
            t += np.sum((Vobs - Vp)**2 / err**2)
            n_tot += len(R)
        total_chi2.append(t / n_tot)

    total_chi2 = np.array(total_chi2)
    idx_min = int(np.argmin(total_chi2))
    kappa_best = kappa_scan[idx_min]
    print(f"  Scanned kappa in [{kappa_scan[0]:.3f}, {kappa_scan[-1]:.3f}]")
    print(f"  Best kappa = {kappa_best:.4f}")
    print(f"  Min chi2_red_total = {total_chi2[idx_min]:.3f}")
    print()
    print(f"  Predicted kappas:")
    print(f"    gs9d geometric-mean (derived):      1.0000")
    print(f"    MOND match (a0_obs/(c*H0)):         {a0_obs/(c*H0):.4f}")
    print(f"    Best fit to SPARC:                  {kappa_best:.4f}")
    print()

    # Full scan output
    print("  Scan table (selected kappa values):")
    show_idx = list(range(0, len(kappa_scan), 5)) + [idx_min]
    show_idx = sorted(set(show_idx))
    print(f"  {'kappa':>8s}  {'chi2_red':>10s}")
    for i in show_idx:
        marker = " <-- best" if i == idx_min else ""
        print(f"  {kappa_scan[i]:8.4f}  {total_chi2[i]:10.3f}{marker}")
    print()

    # ---------------------------------------------------------------
    # Mass-scale systematic trends
    # ---------------------------------------------------------------
    print("=" * 78)
    print("  Mass-scale systematic trends (median chi2_red per category)")
    print("=" * 78)
    print()

    # Bin by Mbar
    cats = [
        ('dwarf',    0,     1e9),
        ('sub-MW',   1e9,   1e10),
        ('MW-like',  1e10,  1e11),
        ('massive',  1e11,  1e13),
    ]
    print(f"  {'Category':12s}  {'N':>4s}  " + "  ".join(f"{m:>11s}" for m in model_names))
    print("  " + "-" * (20 + 13*len(model_names)))
    for lbl, Mlo, Mhi in cats:
        names_in_cat = [n for n in per_galaxy
                        if Mlo*Msun <= per_galaxy[n]['info']['Mbar_9']*1e9*Msun < Mhi*Msun]
        if not names_in_cat:
            continue
        line = f"  {lbl:12s}  {len(names_in_cat):4d}  "
        for m in model_names:
            vals = [results[n][m] for n in names_in_cat]
            line += f"  {np.median(vals):11.2f}"
        print(line)
    print()

    # ---------------------------------------------------------------
    # SUMMARY + FALSIFICATION CRITERIA
    # ---------------------------------------------------------------
    print("=" * 78)
    print("  SUMMARY & FALSIFICATION CRITERIA")
    print("=" * 78)
    print()
    # Compute means over all galaxies
    mem_k1_mean   = np.mean([results[g]['Memb_k=1']     for g in results])
    mem_kM_mean   = np.mean([results[g]['Memb_k=0.183'] for g in results])
    mond_mean     = np.mean([results[g]['MOND']         for g in results])
    form4_mean    = np.mean([results[g]['Form4']        for g in results])
    tgp_mean      = np.mean([results[g]['TGP_std']      for g in results])

    # Compute median BTFR slope and kappa_best
    print(f"  Aggregate results for {len(per_galaxy)} galaxies:")
    print(f"    Mean chi2_red:")
    print(f"      MOND simple:        {mond_mean:8.3f}")
    print(f"      Form 4 (delta=.497): {form4_mean:8.3f}")
    print(f"      TGP standard (gs57): {tgp_mean:8.3f}")
    print(f"      Membrane kappa=1:    {mem_k1_mean:8.3f}")
    print(f"      Membrane kappa_M:    {mem_kM_mean:8.3f}")
    print()
    print(f"    Best-fit global kappa: {kappa_best:.4f}")
    print(f"      (gs9d predicts 1.0;  renormalization suggests {a0_obs/(c*H0):.4f})")
    print()
    print(f"  VERDICTS:")
    if kappa_best > 0.5 and kappa_best < 2.0:
        print(f"    [+] gs9d geometric-mean mechanism VIABLE (kappa ~ 1)")
    else:
        print(f"    [!] gs9d kappa=1 is WRONG by factor {1.0/kappa_best:.2f}")
        print(f"        -> 2pi factor in gs9d argument needs repair")
    if mem_kM_mean < 1.2 * mond_mean:
        print(f"    [+] Membrane (renormalized) matches MOND quality within 20%")
    else:
        print(f"    [!] Membrane worse than MOND by {(mem_kM_mean/mond_mean-1)*100:.0f}%")
        print(f"        -> smooth transition, not hard 3D->2D switch, likely needed")
    print()
    print("  NEXT STEPS (gs62-gs64):")
    print("    gs62: cluster test with Hubble damping (Bullet overshoot dissolution?)")
    print("    gs63: Euclid a0(z) detection design")
    print("    gs64: substrate oscillation derivation from sek05 V(psi)")
    print()
    print("=" * 78)
    print("  gs61 complete.")
    print("=" * 78)


if __name__ == '__main__':
    main()
