#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gs14: SCALING WITH SIZE
========================
Does gamma (or effective a0) scale with system size?

Evidence so far:
  - Massive spirals: prefer gamma=0.40 (gs12)
  - Dwarfs: prefer gamma=0.50 (gs12)
  - Clusters: need gamma>0.50 or a0 ~ 6x larger (gs13)

This script:
  1. Per-galaxy fit with FREE gamma (alpha=0.8 fixed, a0=1.12e-10 fixed)
     -> get best gamma_i for each galaxy
  2. Per-galaxy fit with FREE a0 (alpha=0.8, gamma=0.4 fixed)
     -> get best a0_i for each galaxy
  3. Correlate gamma_i and a0_i with galaxy properties:
     - Vflat (proxy for mass)
     - L[3.6] (luminosity)
     - M_bar (baryon mass from fit)
     - R_last (outermost radius)
     - SB (surface brightness)
     - y_typical (median acceleration parameter)
  4. Add cluster data points to see full scaling
  5. Test: is gamma(M) a smooth function?
"""

import sys
import io
import os
import numpy as np
from scipy.optimize import minimize, minimize_scalar
from scipy.stats import spearmanr, pearsonr
import warnings
warnings.filterwarnings('ignore')

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

G_SI = 6.674e-11
M_sun = 1.989e30
kpc = 3.086e19
Mpc = 3.086e22

a0_ref = 1.12e-10   # galaxy best-fit

# ============================================================
# GALAXY LOADER (from gs12)
# ============================================================
class Galaxy:
    def __init__(self, filepath):
        self.name = os.path.basename(filepath).replace('_rotmod.dat', '')
        lines = []
        self.distance_Mpc = None
        with open(filepath, 'r') as f:
            for line in f:
                ls = line.strip()
                if ls.startswith('#') and 'Distance' in ls:
                    parts = ls.split()
                    for i, p in enumerate(parts):
                        if p == '=' and i+1 < len(parts):
                            try: self.distance_Mpc = float(parts[i+1])
                            except: pass
                elif not ls.startswith('#') and ls:
                    lines.append(ls.split())
        data = np.array(lines, dtype=float)
        self.rad_kpc = data[:, 0]
        self.vobs = data[:, 1]
        self.errv = data[:, 2]
        self.vgas = data[:, 3]
        self.vdisk = data[:, 4]
        self.vbul = data[:, 5]
        self.sbdisk = data[:, 6]
        self.sbbul = data[:, 7]
        self.npts = len(self.rad_kpc)
        self.has_bulge = np.any(self.vbul > 0)
        self.errv_safe = np.maximum(self.errv, np.maximum(3.0, 0.03*np.abs(self.vobs)))
        self.rad_m = self.rad_kpc * kpc
        self.r_last_kpc = self.rad_kpc[-1]

    def compute_gbar(self, Yd, Yb=None):
        if Yb is None: Yb = Yd
        v2 = np.abs(self.vgas)*self.vgas + Yd*self.vdisk**2 + Yb*self.vbul**2
        return v2 * 1e6 / self.rad_m

    def get_vflat(self):
        """Estimate Vflat from outer 1/3 of rotation curve."""
        n = max(1, self.npts // 3)
        return np.median(self.vobs[-n:])

    def get_mbar(self, Yd, Yb=None):
        """Estimate baryonic mass from V^2 * R at outermost point."""
        if Yb is None: Yb = Yd
        v2 = np.abs(self.vgas[-1])*self.vgas[-1] + Yd*self.vdisk[-1]**2 + Yb*self.vbul[-1]**2
        # M_bar ~ v2 * R / G  (in solar masses)
        v2_si = v2 * 1e6  # (km/s)^2 -> (m/s)^2
        R_si = self.rad_m[-1]
        return v2_si * R_si / G_SI / M_sun

    def get_sb_mean(self):
        """Mean disk surface brightness."""
        return np.mean(self.sbdisk[self.sbdisk > 0]) if np.any(self.sbdisk > 0) else 0.0


def load_galaxies(data_dir):
    galaxies = []
    for fname in sorted(os.listdir(data_dir)):
        if fname.endswith('_rotmod.dat'):
            try:
                gal = Galaxy(os.path.join(data_dir, fname))
                if gal.npts >= 5:
                    galaxies.append(gal)
            except: pass
    return galaxies


def load_master_table(filepath):
    """Load SPARC master table (Lelli+2016c.mrt)."""
    props = {}
    with open(filepath, 'r') as f:
        lines = f.readlines()
    # Find last separator
    last_sep = 0
    for i, line in enumerate(lines):
        if line.startswith('---'):
            last_sep = i
    for line in lines[last_sep+1:]:
        parts = line.split()
        if len(parts) >= 17:
            name = parts[0]
            try:
                props[name] = {
                    'T': float(parts[1]),
                    'D': float(parts[2]),
                    'Inc': float(parts[5]),
                    'L36': float(parts[7]),   # 10^9 Lsun
                    'MHI': float(parts[13]),   # 10^9 Msun
                    'Vflat': float(parts[15]),
                    'Q': int(parts[17]),
                }
            except:
                pass
    return props


# ============================================================
# MODELS
# ============================================================
def nu_exp(y, alpha, gamma):
    return 1.0 + np.exp(-np.power(y, alpha)) / np.power(y, gamma)

def nu_mcgaugh(y):
    return 1.0 / (1.0 - np.exp(-np.sqrt(np.maximum(y, 1e-30))))


# ============================================================
# PER-GALAXY FITTERS
# ============================================================
def chi2_galaxy(galaxy, a0, alpha, gamma, Yd, Yb=None):
    """Compute chi2 for a galaxy with given params."""
    if Yb is None: Yb = Yd
    gbar = galaxy.compute_gbar(Yd, Yb)
    mask = gbar > 0
    if np.sum(mask) < 3: return 1e10
    y = gbar[mask] / a0
    y = np.maximum(y, 1e-10)
    try:
        nu = nu_exp(y, alpha, gamma)
    except: return 1e10
    if not np.all(np.isfinite(nu)): return 1e10
    gobs = gbar[mask] * nu
    vp = np.sqrt(np.abs(gobs) * galaxy.rad_m[mask]) / 1e3
    return np.sum(((galaxy.vobs[mask] - vp) / galaxy.errv_safe[mask])**2)


def fit_galaxy_free_gamma(galaxy, a0=a0_ref, alpha=0.8):
    """Fit galaxy with free gamma and free Yd. Return best (gamma, Yd, chi2)."""
    best = {'gamma': 0.4, 'Yd': 0.5, 'Yb': 0.5, 'chi2': 1e10}

    for gamma_init in np.arange(0.15, 0.75, 0.05):
        for yd_init in [0.3, 0.5, 0.8, 1.2]:
            def obj(params):
                gamma = params[0]
                log_Yd = params[1]
                log_Yb = params[2] if galaxy.has_bulge else params[1]
                if gamma < 0.05 or gamma > 1.0: return 1e10
                Yd = 10**log_Yd
                Yb = 10**log_Yb
                if Yd < 0.05 or Yd > 15: return 1e10
                if Yb < 0.05 or Yb > 15: return 1e10
                return chi2_galaxy(galaxy, a0, alpha, gamma, Yd, Yb)

            x0 = [gamma_init, np.log10(yd_init)]
            bounds = [(0.05, 1.0), (-1.3, 1.2)]
            if galaxy.has_bulge:
                x0.append(np.log10(yd_init))
                bounds.append((-1.3, 1.2))

            try:
                res = minimize(obj, x0, method='L-BFGS-B', bounds=bounds,
                             options={'maxiter': 500})
                if res.fun < best['chi2']:
                    best['gamma'] = res.x[0]
                    best['Yd'] = 10**res.x[1]
                    best['Yb'] = 10**(res.x[2] if galaxy.has_bulge else res.x[1])
                    best['chi2'] = res.fun
            except:
                pass

    return best


def fit_galaxy_free_a0(galaxy, alpha=0.8, gamma=0.4):
    """Fit galaxy with free a0 and free Yd. Return best (a0, Yd, chi2)."""
    best = {'a0': a0_ref, 'Yd': 0.5, 'Yb': 0.5, 'chi2': 1e10}

    for a0_init_log in [-11.0, -10.5, -10.0, -9.5, -9.0]:
        for yd_init in [0.3, 0.5, 0.8, 1.2]:
            def obj(params):
                log_a0 = params[0]
                log_Yd = params[1]
                log_Yb = params[2] if galaxy.has_bulge else params[1]
                a0 = 10**log_a0
                Yd = 10**log_Yd
                Yb = 10**log_Yb
                if Yd < 0.05 or Yd > 15: return 1e10
                if Yb < 0.05 or Yb > 15: return 1e10
                return chi2_galaxy(galaxy, a0, alpha, gamma, Yd, Yb)

            x0 = [a0_init_log, np.log10(yd_init)]
            bounds = [(-12.0, -8.0), (-1.3, 1.2)]
            if galaxy.has_bulge:
                x0.append(np.log10(yd_init))
                bounds.append((-1.3, 1.2))

            try:
                res = minimize(obj, x0, method='L-BFGS-B', bounds=bounds,
                             options={'maxiter': 500})
                if res.fun < best['chi2']:
                    best['a0'] = 10**res.x[0]
                    best['Yd'] = 10**res.x[1]
                    best['Yb'] = 10**(res.x[2] if galaxy.has_bulge else res.x[1])
                    best['chi2'] = res.fun
            except:
                pass

    return best


def fit_galaxy_free_alpha_gamma(galaxy, a0=a0_ref):
    """Fit galaxy with free alpha, gamma, Yd. Return best (alpha, gamma, Yd, chi2)."""
    best = {'alpha': 0.8, 'gamma': 0.4, 'Yd': 0.5, 'Yb': 0.5, 'chi2': 1e10}

    for alpha_init in [0.4, 0.6, 0.8, 1.0]:
        for gamma_init in [0.2, 0.35, 0.5, 0.65]:
            for yd_init in [0.3, 0.6, 1.0]:
                def obj(params):
                    alpha = params[0]
                    gamma = params[1]
                    log_Yd = params[2]
                    log_Yb = params[3] if galaxy.has_bulge else params[2]
                    if alpha < 0.1 or alpha > 2.0: return 1e10
                    if gamma < 0.05 or gamma > 1.5: return 1e10
                    Yd = 10**log_Yd
                    Yb = 10**log_Yb
                    if Yd < 0.05 or Yd > 15: return 1e10
                    return chi2_galaxy(galaxy, a0, alpha, gamma, Yd, Yb)

                x0 = [alpha_init, gamma_init, np.log10(yd_init)]
                bounds = [(0.1, 2.0), (0.05, 1.5), (-1.3, 1.2)]
                if galaxy.has_bulge:
                    x0.append(np.log10(yd_init))
                    bounds.append((-1.3, 1.2))

                try:
                    res = minimize(obj, x0, method='L-BFGS-B', bounds=bounds,
                                 options={'maxiter': 500})
                    if res.fun < best['chi2']:
                        best['alpha'] = res.x[0]
                        best['gamma'] = res.x[1]
                        best['Yd'] = 10**res.x[2]
                        best['Yb'] = 10**(res.x[3] if galaxy.has_bulge else res.x[2])
                        best['chi2'] = res.fun
                except:
                    pass

    return best


# ============================================================
# MAIN
# ============================================================
def main():
    print("=" * 80)
    print("gs14: SCALING WITH SYSTEM SIZE")
    print("Does gamma or a0 depend on galaxy mass/size?")
    print("=" * 80)

    data_dir = os.path.join(os.path.dirname(__file__), 'Rotmod_LTG')
    master_file = os.path.join(os.path.dirname(__file__), 'SPARC_Lelli2016c.mrt')

    galaxies = load_galaxies(data_dir)
    props = load_master_table(master_file)
    print(f"\nLoaded {len(galaxies)} galaxies, {len(props)} master table entries")

    # ================================================================
    # PHASE 1: Per-galaxy free gamma fit
    # ================================================================
    print("\n" + "=" * 80)
    print("PHASE 1: FIT EACH GALAXY WITH FREE GAMMA (alpha=0.8, a0=1.12e-10 fixed)")
    print("=" * 80)

    results = []
    for i, gal in enumerate(galaxies):
        if (i+1) % 20 == 0:
            print(f"  Fitting galaxy {i+1}/{len(galaxies)}...", flush=True)

        res_gamma = fit_galaxy_free_gamma(gal)

        # Also get fixed-gamma chi2 for comparison
        chi2_04 = 1e10
        for yd_init in [0.3, 0.5, 0.8, 1.2]:
            def obj04(params):
                Yd = 10**params[0]
                Yb = 10**(params[1] if gal.has_bulge else params[0])
                if Yd < 0.05 or Yd > 15: return 1e10
                return chi2_galaxy(gal, a0_ref, 0.8, 0.4, Yd, Yb)
            x0 = [np.log10(yd_init)]
            bnds = [(-1.3, 1.2)]
            if gal.has_bulge:
                x0.append(np.log10(yd_init))
                bnds.append((-1.3, 1.2))
            try:
                r04 = minimize(obj04, x0, method='L-BFGS-B', bounds=bnds)
                if r04.fun < chi2_04: chi2_04 = r04.fun
            except: pass

        chi2_05 = 1e10
        for yd_init in [0.3, 0.5, 0.8, 1.2]:
            def obj05(params):
                Yd = 10**params[0]
                Yb = 10**(params[1] if gal.has_bulge else params[0])
                if Yd < 0.05 or Yd > 15: return 1e10
                return chi2_galaxy(gal, a0_ref, 0.8, 0.5, Yd, Yb)
            x0 = [np.log10(yd_init)]
            bnds = [(-1.3, 1.2)]
            if gal.has_bulge:
                x0.append(np.log10(yd_init))
                bnds.append((-1.3, 1.2))
            try:
                r05 = minimize(obj05, x0, method='L-BFGS-B', bounds=bnds)
                if r05.fun < chi2_05: chi2_05 = r05.fun
            except: pass

        # Galaxy properties
        vflat = gal.get_vflat()
        mbar = gal.get_mbar(res_gamma['Yd'], res_gamma['Yb'])
        sb = gal.get_sb_mean()

        # y_typical: median y at the outer half of RC
        gbar_outer = gal.compute_gbar(res_gamma['Yd'], res_gamma['Yb'])
        n_outer = max(1, gal.npts // 2)
        gbar_median_outer = np.median(gbar_outer[-n_outer:][gbar_outer[-n_outer:] > 0])
        y_typical = gbar_median_outer / a0_ref if gbar_median_outer > 0 else 0.1

        # Master table properties
        vflat_mt = props.get(gal.name, {}).get('Vflat', 0)
        L36 = props.get(gal.name, {}).get('L36', 0)
        Q = props.get(gal.name, {}).get('Q', 9)
        D = props.get(gal.name, {}).get('D', 0)

        results.append({
            'name': gal.name,
            'npts': gal.npts,
            'gamma_best': res_gamma['gamma'],
            'Yd_best': res_gamma['Yd'],
            'chi2_free': res_gamma['chi2'],
            'chi2_04': chi2_04,
            'chi2_05': chi2_05,
            'dchi2_04': chi2_04 - res_gamma['chi2'],  # >0 if free gamma better
            'dchi2_05': chi2_05 - res_gamma['chi2'],
            'vflat': vflat,
            'vflat_mt': vflat_mt,
            'mbar': mbar,
            'L36': L36,
            'r_last': gal.r_last_kpc,
            'sb_mean': sb,
            'y_typical': y_typical,
            'Q': Q,
            'D': D,
        })

    # ================================================================
    # RESULTS: Distribution of gamma
    # ================================================================
    print("\n" + "=" * 80)
    print("SECTION 1: DISTRIBUTION OF BEST-FIT GAMMA")
    print("=" * 80)

    gammas = np.array([r['gamma_best'] for r in results])
    print(f"\n  N galaxies: {len(gammas)}")
    print(f"  gamma: mean={np.mean(gammas):.3f}, median={np.median(gammas):.3f}, "
          f"std={np.std(gammas):.3f}")
    print(f"  gamma range: {np.min(gammas):.3f} to {np.max(gammas):.3f}")

    # Histogram bins
    bins = [0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95]
    counts, _ = np.histogram(gammas, bins=bins)
    print(f"\n  Histogram:")
    for i in range(len(counts)):
        bar = '#' * counts[i]
        print(f"  {bins[i]:.2f}-{bins[i+1]:.2f}: {counts[i]:3d} {bar}")

    # ================================================================
    # SECTION 2: gamma vs galaxy properties
    # ================================================================
    print("\n" + "=" * 80)
    print("SECTION 2: CORRELATIONS gamma vs GALAXY PROPERTIES")
    print("=" * 80)

    # Only use Q<=2 galaxies for cleaner results
    good = [r for r in results if r['Q'] <= 2]
    print(f"\n  Using Q<=2 galaxies: {len(good)}")

    gammas_g = np.array([r['gamma_best'] for r in good])

    correlations = []
    for prop_name, prop_key, log_it in [
        ('Vflat (measured)', 'vflat', True),
        ('Vflat (master table)', 'vflat_mt', True),
        ('M_bar (from fit)', 'mbar', True),
        ('L[3.6] (luminosity)', 'L36', True),
        ('R_last (kpc)', 'r_last', True),
        ('SB_mean', 'sb_mean', True),
        ('y_typical', 'y_typical', True),
        ('Distance', 'D', True),
    ]:
        vals = np.array([r[prop_key] for r in good])
        mask = (vals > 0) & np.isfinite(vals)
        if np.sum(mask) < 10:
            continue

        x = np.log10(vals[mask]) if log_it else vals[mask]
        y = gammas_g[mask]

        rho, pval = spearmanr(x, y)
        rho_p, pval_p = pearsonr(x, y)

        sig = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else ""
        correlations.append((prop_name, rho, pval, rho_p, pval_p, np.sum(mask)))

        print(f"\n  gamma vs log({prop_name}):")
        print(f"    N={np.sum(mask)}, Spearman rho={rho:+.3f} (p={pval:.4f}) {sig}")
        print(f"    Pearson r={rho_p:+.3f} (p={pval_p:.4f})")

    # ================================================================
    # SECTION 3: Binned gamma by Vflat
    # ================================================================
    print("\n" + "=" * 80)
    print("SECTION 3: BINNED GAMMA BY VELOCITY (MASS PROXY)")
    print("=" * 80)

    vflats = np.array([r['vflat'] for r in good])
    gammas_g = np.array([r['gamma_best'] for r in good])

    # Velocity bins
    v_bins = [(0, 50), (50, 80), (80, 120), (120, 170), (170, 250), (250, 400)]
    print(f"\n  {'Vflat range':<15s} {'N':<5s} {'mean gamma':<12s} {'median gamma':<13s} "
          f"{'std gamma':<10s} {'mean dchi2(0.4)':<16s} {'mean dchi2(0.5)':<16s}")
    print("  " + "-" * 90)

    for vlo, vhi in v_bins:
        mask = (vflats >= vlo) & (vflats < vhi)
        if np.sum(mask) < 2: continue
        g = gammas_g[mask]
        dc04 = np.array([r['dchi2_04'] for r in good])[mask]
        dc05 = np.array([r['dchi2_05'] for r in good])[mask]
        print(f"  {vlo:3d}-{vhi:3d} km/s    {np.sum(mask):<5d} {np.mean(g):<12.3f} "
              f"{np.median(g):<13.3f} {np.std(g):<10.3f} {np.mean(dc04):<16.2f} "
              f"{np.mean(dc05):<16.2f}")

    # ================================================================
    # SECTION 4: Binned gamma by y_typical
    # ================================================================
    print("\n" + "=" * 80)
    print("SECTION 4: BINNED GAMMA BY y_typical (ACCELERATION REGIME)")
    print("=" * 80)

    y_typs = np.array([r['y_typical'] for r in good])

    y_bins = [(0, 0.01), (0.01, 0.03), (0.03, 0.1), (0.1, 0.3), (0.3, 1.0), (1.0, 100)]
    print(f"\n  {'y range':<15s} {'N':<5s} {'mean gamma':<12s} {'median gamma':<13s} "
          f"{'std gamma':<10s} {'regime'}")
    print("  " + "-" * 70)

    for ylo, yhi in y_bins:
        mask = (y_typs >= ylo) & (y_typs < yhi)
        if np.sum(mask) < 2: continue
        g = gammas_g[mask]
        regime = "deep MOND" if yhi <= 0.1 else ("transition" if yhi <= 1.0 else "Newtonian")
        print(f"  {ylo:.2f}-{yhi:.2f}       {np.sum(mask):<5d} {np.mean(g):<12.3f} "
              f"{np.median(g):<13.3f} {np.std(g):<10.3f} {regime}")

    # ================================================================
    # SECTION 5: Top-10 galaxies preferring low vs high gamma
    # ================================================================
    print("\n" + "=" * 80)
    print("SECTION 5: GALAXIES WITH EXTREME GAMMA PREFERENCES")
    print("=" * 80)

    sorted_results = sorted(good, key=lambda r: r['gamma_best'])

    print(f"\n  --- LOWEST GAMMA (most Exp-like, gamma << 0.5) ---")
    print(f"  {'Galaxy':<15s} {'gamma':<7s} {'Vflat':<7s} {'Mbar':<10s} {'y_typ':<8s} "
          f"{'dchi2(0.5)':<12s} {'Q'}")
    for r in sorted_results[:15]:
        print(f"  {r['name']:<15s} {r['gamma_best']:<7.3f} {r['vflat']:<7.0f} "
              f"{r['mbar']:<10.2e} {r['y_typical']:<8.4f} {r['dchi2_05']:<12.1f} {r['Q']}")

    print(f"\n  --- HIGHEST GAMMA (most MOND-like, gamma >> 0.5) ---")
    for r in sorted_results[-15:]:
        print(f"  {r['name']:<15s} {r['gamma_best']:<7.3f} {r['vflat']:<7.0f} "
              f"{r['mbar']:<10.2e} {r['y_typical']:<8.4f} {r['dchi2_05']:<12.1f} {r['Q']}")

    # ================================================================
    # SECTION 6: Per-galaxy free a0 (alpha=0.8, gamma=0.4 fixed)
    # ================================================================
    print("\n" + "=" * 80)
    print("PHASE 2: FIT EACH GALAXY WITH FREE a0 (alpha=0.8, gamma=0.4 fixed)")
    print("=" * 80)

    for i, gal in enumerate(galaxies):
        if (i+1) % 20 == 0:
            print(f"  Fitting galaxy {i+1}/{len(galaxies)}...", flush=True)

        res_a0 = fit_galaxy_free_a0(gal)
        # Add to existing results
        for r in results:
            if r['name'] == gal.name:
                r['a0_best'] = res_a0['a0']
                r['Yd_a0'] = res_a0['Yd']
                r['chi2_a0'] = res_a0['chi2']
                break

    # Correlate a0 with properties
    print("\n" + "=" * 80)
    print("SECTION 6: CORRELATIONS a0_best vs GALAXY PROPERTIES")
    print("=" * 80)

    good = [r for r in results if r['Q'] <= 2 and 'a0_best' in r]
    a0s = np.array([r['a0_best'] for r in good])

    print(f"\n  a0: mean={np.mean(a0s):.3e}, median={np.median(a0s):.3e}, "
          f"std={np.std(a0s):.3e}")
    print(f"  a0/a0_ref: mean={np.mean(a0s)/a0_ref:.2f}, "
          f"median={np.median(a0s)/a0_ref:.2f}")

    for prop_name, prop_key, log_it in [
        ('Vflat', 'vflat', True),
        ('M_bar', 'mbar', True),
        ('L[3.6]', 'L36', True),
        ('R_last', 'r_last', True),
        ('y_typical', 'y_typical', True),
    ]:
        vals = np.array([r[prop_key] for r in good])
        mask = (vals > 0) & np.isfinite(vals)
        if np.sum(mask) < 10: continue

        x = np.log10(vals[mask]) if log_it else vals[mask]
        y = np.log10(a0s[mask])

        rho, pval = spearmanr(x, y)
        sig = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else ""
        print(f"\n  log(a0) vs log({prop_name}):")
        print(f"    N={np.sum(mask)}, Spearman rho={rho:+.3f} (p={pval:.4f}) {sig}")

    # Binned a0 by Vflat
    print(f"\n  Binned a0 by Vflat:")
    vflats = np.array([r['vflat'] for r in good])
    a0s_arr = np.array([r['a0_best'] for r in good])

    print(f"  {'Vflat range':<15s} {'N':<5s} {'mean a0/a0_ref':<16s} {'median a0/a0_ref':<18s}")
    print("  " + "-" * 54)
    for vlo, vhi in v_bins:
        mask = (vflats >= vlo) & (vflats < vhi)
        if np.sum(mask) < 2: continue
        a = a0s_arr[mask] / a0_ref
        print(f"  {vlo:3d}-{vhi:3d} km/s    {np.sum(mask):<5d} {np.mean(a):<16.3f} "
              f"{np.median(a):<18.3f}")

    # ================================================================
    # SECTION 7: FULL SCALE — galaxies + clusters
    # ================================================================
    print("\n" + "=" * 80)
    print("SECTION 7: FULL SCALE — FROM DWARFS TO CLUSTERS")
    print("=" * 80)

    # Galaxy data points (binned)
    print(f"\n  System size scaling:")
    print(f"  {'System':<25s} {'Vflat/sigma':<12s} {'M_bar':<12s} {'gamma_eff':<10s} "
          f"{'a0_eff':<12s} {'a0/a0_ref':<10s}")
    print("  " + "-" * 80)

    for vlo, vhi in v_bins:
        mask_v = np.array([(r['vflat'] >= vlo and r['vflat'] < vhi and r['Q'] <= 2)
                          for r in results])
        if np.sum(mask_v) < 3: continue
        sub = [r for r, m in zip(results, mask_v) if m]
        g_med = np.median([r['gamma_best'] for r in sub])
        a0_med = np.median([r.get('a0_best', a0_ref) for r in sub])
        m_med = np.median([r['mbar'] for r in sub if r['mbar'] > 0])
        vmed = np.median([r['vflat'] for r in sub])
        print(f"  Galaxy ({vlo}-{vhi} km/s)    {vmed:<12.0f} {m_med:<12.2e} "
              f"{g_med:<10.3f} {a0_med:<12.3e} {a0_med/a0_ref:<10.2f}")

    # Add cluster data points (from gs13)
    cluster_points = [
        ("Groups (Fornax etc)", 200, 3e13, 0.37, 6.2e-10),
        ("Poor clusters", 500, 2e14, 0.37, 6.2e-10),
        ("Rich clusters (Coma)", 1000, 7e14, 0.37, 6.2e-10),
        ("Massive clusters", 1500, 1e15, 0.37, 6.2e-10),
    ]
    for name, v, m, g, a0 in cluster_points:
        print(f"  {name:<25s} {v:<12.0f} {m:<12.2e} {g:<10.3f} {a0:<12.3e} "
              f"{a0/a0_ref:<10.2f}")

    # ================================================================
    # SECTION 8: chi2 preference map
    # ================================================================
    print("\n" + "=" * 80)
    print("SECTION 8: chi2 PREFERENCE MAP — gamma=0.4 vs gamma=0.5")
    print("=" * 80)

    good = [r for r in results if r['Q'] <= 2]

    # dchi2 = chi2(gamma=0.4) - chi2(gamma=0.5): negative = 0.4 better
    pref_04 = sum(1 for r in good if r['chi2_04'] < r['chi2_05'])
    pref_05 = sum(1 for r in good if r['chi2_05'] < r['chi2_04'])
    equal = sum(1 for r in good if abs(r['chi2_04'] - r['chi2_05']) < 1)

    print(f"\n  Prefer gamma=0.4: {pref_04} galaxies")
    print(f"  Prefer gamma=0.5: {pref_05} galaxies")
    print(f"  No strong preference (|dchi2|<1): {equal}")

    dchi2_all = np.array([r['chi2_04'] - r['chi2_05'] for r in good])
    vflats_g = np.array([r['vflat'] for r in good])

    rho, pval = spearmanr(vflats_g, dchi2_all)
    print(f"\n  Correlation: dchi2(0.4 vs 0.5) vs Vflat:")
    print(f"    Spearman rho={rho:+.3f} (p={pval:.4f})")
    print(f"    Negative rho = massive galaxies prefer gamma=0.4")

    # Binned
    print(f"\n  {'Vflat range':<15s} {'N':<5s} {'prefer 0.4':<12s} {'prefer 0.5':<12s} "
          f"{'mean dchi2':<12s}")
    print("  " + "-" * 56)
    for vlo, vhi in v_bins:
        mask = (vflats_g >= vlo) & (vflats_g < vhi)
        if np.sum(mask) < 2: continue
        dc = dchi2_all[mask]
        n04 = np.sum(dc < -1)
        n05 = np.sum(dc > 1)
        print(f"  {vlo:3d}-{vhi:3d} km/s    {np.sum(mask):<5d} {n04:<12d} {n05:<12d} "
              f"{np.mean(dc):<12.2f}")

    # ================================================================
    # SECTION 9: PHASE 3 — per-galaxy free (alpha, gamma) for
    # representative sample
    # ================================================================
    print("\n" + "=" * 80)
    print("PHASE 3: FREE (alpha, gamma) FOR 30 REPRESENTATIVE GALAXIES")
    print("=" * 80)

    # Pick ~30 galaxies spanning mass range
    good_sorted = sorted([r for r in results if r['Q'] <= 2 and r['vflat'] > 20],
                        key=lambda r: r['vflat'])
    n_total = len(good_sorted)
    indices = np.linspace(0, n_total-1, min(30, n_total), dtype=int)
    rep_names = [good_sorted[i]['name'] for i in indices]

    print(f"\n  Fitting {len(rep_names)} galaxies with free (alpha, gamma, Yd)...")

    rep_results = []
    for name in rep_names:
        gal = None
        for g in galaxies:
            if g.name == name:
                gal = g
                break
        if gal is None: continue

        res = fit_galaxy_free_alpha_gamma(gal)
        vf = gal.get_vflat()
        mb = gal.get_mbar(res['Yd'], res['Yb'])
        rep_results.append({
            'name': name,
            'alpha': res['alpha'],
            'gamma': res['gamma'],
            'gamma_over_alpha': res['gamma']/res['alpha'] if res['alpha'] > 0 else 0,
            'Yd': res['Yd'],
            'chi2': res['chi2'],
            'vflat': vf,
            'mbar': mb,
        })

    print(f"\n  {'Galaxy':<15s} {'Vflat':<7s} {'alpha':<7s} {'gamma':<7s} "
          f"{'g/a':<6s} {'Yd':<7s} {'chi2/dof':<9s}")
    print("  " + "-" * 65)
    for r in sorted(rep_results, key=lambda x: x['vflat']):
        npts = 0
        for g in galaxies:
            if g.name == r['name']:
                npts = g.npts
                break
        chi2_red = r['chi2'] / max(1, npts - 3)
        print(f"  {r['name']:<15s} {r['vflat']:<7.0f} {r['alpha']:<7.3f} "
              f"{r['gamma']:<7.3f} {r['gamma_over_alpha']:<6.3f} {r['Yd']:<7.3f} "
              f"{chi2_red:<9.3f}")

    # Correlations for representative sample
    if len(rep_results) >= 10:
        alphas_r = np.array([r['alpha'] for r in rep_results])
        gammas_r = np.array([r['gamma'] for r in rep_results])
        vflats_r = np.array([r['vflat'] for r in rep_results])
        ga_ratios = np.array([r['gamma_over_alpha'] for r in rep_results])

        print(f"\n  Correlations for representative sample:")
        rho, pval = spearmanr(np.log10(vflats_r), gammas_r)
        print(f"    gamma vs log(Vflat): rho={rho:+.3f}, p={pval:.4f}")
        rho, pval = spearmanr(np.log10(vflats_r), alphas_r)
        print(f"    alpha vs log(Vflat): rho={rho:+.3f}, p={pval:.4f}")
        rho, pval = spearmanr(np.log10(vflats_r), ga_ratios)
        print(f"    gamma/alpha vs log(Vflat): rho={rho:+.3f}, p={pval:.4f}")

        print(f"\n  gamma/alpha stats:")
        print(f"    mean={np.mean(ga_ratios):.3f}, median={np.median(ga_ratios):.3f}, "
              f"std={np.std(ga_ratios):.3f}")
        print(f"    min={np.min(ga_ratios):.3f}, max={np.max(ga_ratios):.3f}")

    # ================================================================
    # SUMMARY
    # ================================================================
    print("\n" + "=" * 80)
    print("SUMMARY: SCALING WITH SIZE")
    print("=" * 80)

    good = [r for r in results if r['Q'] <= 2]
    gammas_all = np.array([r['gamma_best'] for r in good])
    vflats_all = np.array([r['vflat'] for r in good])

    rho_main, pval_main = spearmanr(vflats_all, gammas_all)

    print(f"""
  1. Per-galaxy free gamma distribution:
     mean = {np.mean(gammas_all):.3f}, median = {np.median(gammas_all):.3f},
     std = {np.std(gammas_all):.3f}

  2. Correlation gamma vs Vflat (mass proxy):
     Spearman rho = {rho_main:+.3f}, p = {pval_main:.4f}
     {'SIGNIFICANT' if pval_main < 0.05 else 'NOT significant'} trend

  3. chi2 preference (gamma=0.4 vs 0.5):
     Prefer 0.4: {pref_04}, Prefer 0.5: {pref_05}
     {'Massive galaxies prefer lower gamma' if rho_main < -0.1 else 'No clear mass dependence'}

  4. Cluster extension (from gs13):
     Clusters need a0 ~ 6x galaxy => gamma_eff needs to be HIGHER
     Or: same gamma but different a0 at cluster scale

  5. Key question: is gamma truly mass-dependent, or is the scatter
     consistent with measurement noise on a UNIVERSAL gamma?
     gamma std = {np.std(gammas_all):.3f} around mean {np.mean(gammas_all):.3f}
     -> {'Large scatter suggests possible mass dependence' if np.std(gammas_all) > 0.1 else 'Small scatter consistent with universal gamma'}
""")

    print("=" * 80)
    print("DONE")
    print("=" * 80)


if __name__ == '__main__':
    main()
