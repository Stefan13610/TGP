#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gs10: External Field Effect & Environment Analysis
===================================================

Goal: Check whether galaxy ENVIRONMENT affects rotation curve fits.
In MOND, the External Field Effect (EFE) from a host galaxy or cluster
modifies internal dynamics. In TGP, g_ext could shift the dimensional
transition scale.

Steps:
  1. Load all 175 SPARC galaxies + master table (distances, luminosities)
  2. Classify environment:
     - Ursa Major cluster (distance method 4 in SPARC)
     - Known groups/pairs from literature
     - Isolated field galaxies
  3. Fit each galaxy with McGaugh + MOND simple (no environment)
  4. Check if residuals correlate with environment proxies
  5. Fit with simple EFE model: g_obs = g_bar * nu((g_bar + g_ext)/a0)
  6. Compare: does adding g_ext improve fits for non-isolated galaxies?

Data: SPARC Rotmod_LTG (175 galaxies) + SPARC_Lelli2016c.mrt (master table)
"""

import sys
import io
import os
import re
import numpy as np
from scipy.optimize import minimize, minimize_scalar, differential_evolution
import warnings
warnings.filterwarnings('ignore')

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

# ============================================================
# CONSTANTS
# ============================================================
G = 6.674e-11
c = 2.998e8
M_sun = 1.989e30
kpc = 3.086e19
Mpc = 3.086e22
H0 = 67.4e3 / (1e6 * kpc)
a0_fiducial = 1.2e-10  # m/s^2

# ============================================================
# INTERPOLATION FUNCTIONS
# ============================================================
def nu_mcgaugh(y):
    """McGaugh: g_obs = g_bar / (1 - exp(-sqrt(y)))"""
    sy = np.sqrt(np.maximum(y, 1e-30))
    return 1.0 / (1.0 - np.exp(-sy))

def nu_mond_simple(y):
    """MOND simple: nu = 1/2 + sqrt(1/4 + 1/y)"""
    return 0.5 + np.sqrt(0.25 + 1.0 / np.maximum(y, 1e-30))

def nu_mcgaugh_efe(y_total, y_bar):
    """McGaugh with EFE: nu evaluated at y_total = (g_bar + g_ext)/a0
    but applied to g_bar only: g_obs = g_bar * nu(y_total)

    Actually the proper MOND EFE is more subtle:
    g_obs = g_bar * nu(y_total) where y_total includes external field
    The external field raises the effective y, REDUCING the MOND boost
    """
    sy = np.sqrt(np.maximum(y_total, 1e-30))
    return 1.0 / (1.0 - np.exp(-sy))

# ============================================================
# GALAXY DATA
# ============================================================
class Galaxy:
    def __init__(self, filepath):
        self.name = os.path.basename(filepath).replace('_rotmod.dat', '')
        self.distance = None
        lines = []
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('# Distance'):
                    self.distance = float(line.split('=')[1].strip().split()[0])
                elif not line.startswith('#') and line:
                    lines.append(line.split())
        data = np.array(lines, dtype=float)
        self.rad_kpc = data[:, 0]
        self.vobs = data[:, 1]
        self.errv = data[:, 2]
        self.vgas = data[:, 3]
        self.vdisk = data[:, 4]
        self.vbul = data[:, 5]
        self.npts = len(self.rad_kpc)
        self.has_bulge = np.any(self.vbul > 0)
        self.errv_safe = np.maximum(self.errv, np.maximum(3.0, 0.03 * np.abs(self.vobs)))
        self.rad_m = self.rad_kpc * kpc

    def compute_gbar(self, Yd, Yb=None):
        if Yb is None:
            Yb = Yd
        v2_bar = np.abs(self.vgas) * self.vgas + Yd * self.vdisk**2 + Yb * self.vbul**2
        return v2_bar * 1e6 / self.rad_m  # m/s^2

    def compute_gobs(self):
        return (self.vobs * 1e3)**2 / self.rad_m


def load_master_table(filepath):
    """Parse SPARC_Lelli2016c.mrt fixed-format table."""
    props = {}
    with open(filepath, 'r') as f:
        lines = f.readlines()

    # Find the LAST separator line — data follows it
    last_sep = 0
    for i, line in enumerate(lines):
        if line.startswith('---'):
            last_sep = i

    for line in lines[last_sep+1:]:
        parts = line.split()
        if len(parts) < 18:
            continue
        try:
            name = parts[0]
            T = int(parts[1])
            D = float(parts[2])
            e_D = float(parts[3])
            f_D = int(parts[4])
            Inc = float(parts[5])
            e_Inc = float(parts[6])
            L36 = float(parts[7])    # 10^9 Lsun
            Reff = float(parts[9])
            SBeff = float(parts[10])
            Rdisk = float(parts[11])
            SBdisk = float(parts[12])
            MHI = float(parts[13])   # 10^9 Msun
            RHI = float(parts[14])
            Vflat = float(parts[15])
            e_Vflat = float(parts[16])
            Q = int(parts[17])

            props[name] = {
                'T': T, 'D': D, 'e_D': e_D, 'f_D': f_D,
                'Inc': Inc, 'e_Inc': e_Inc,
                'L36': L36,
                'Reff': Reff, 'SBeff': SBeff,
                'Rdisk': Rdisk, 'SBdisk': SBdisk,
                'MHI': MHI, 'RHI': RHI,
                'Vflat': Vflat, 'e_Vflat': e_Vflat,
                'Q': Q,
            }
        except (ValueError, IndexError):
            continue

    return props


def classify_environment(props):
    """Classify galaxies by environment using available SPARC data.

    Proxies:
    1. Distance method 4 = Ursa Major cluster member
    2. Known satellite/group members from names and distances
    3. Luminosity as proxy for isolation (dwarfs more likely satellites)
    4. Distance: very nearby galaxies often in Local Group environment
    """
    env = {}

    # Ursa Major cluster members (f_D = 4, D ~ 18 Mpc)
    uma_members = set()

    # Known group/cluster associations from SPARC + literature
    # Local Group and nearby (D < 5 Mpc) — subject to MW/M31 external field
    local_group = set()

    # Virgo-influenced (D ~ 15-20 Mpc, but not UMa)
    virgo_region = set()

    for name, p in props.items():
        env_class = 'field'  # default
        g_ext_estimate = 0.0

        # 1. Ursa Major cluster (f_D = 4)
        if p['f_D'] == 4:
            env_class = 'uma_cluster'
            uma_members.add(name)
            # UMa cluster: sigma ~ 150 km/s, R ~ 1 Mpc
            # g_ext ~ sigma^2 / R ~ (150e3)^2 / (1e6 * kpc) ~ 7.3e-13 m/s^2
            # Or: g_ext from cluster potential at typical position
            # More precise: UMa is a loose group, sigma ~ 148 km/s
            # Using NFW: g_ext ~ 0.01 * a0 at typical member distance
            g_ext_estimate = 0.01 * a0_fiducial

        # 2. Very nearby (D < 4 Mpc) — Local Group environment
        #    MW gravitational field at D ~ 1-4 Mpc:
        #    g_MW ~ G*M_MW / D^2 ~ 6.67e-11 * 1e12 * 2e30 / (D_m)^2
        elif p['D'] < 4.0:
            env_class = 'local_volume'
            D_m = p['D'] * Mpc
            M_MW = 1.0e12 * M_sun
            g_ext_estimate = G * M_MW / D_m**2

        # 3. Nearby (4-10 Mpc) — still potentially influenced
        elif p['D'] < 10.0:
            env_class = 'nearby'
            # Weaker external field
            D_m = p['D'] * Mpc
            M_host = 5e11 * M_sun  # typical nearby massive galaxy
            g_ext_estimate = G * M_host / D_m**2

        # 4. Distant field galaxies (D > 10 Mpc, not UMa)
        else:
            env_class = 'field'
            g_ext_estimate = 0.0  # assume isolated

        # Additional flag: low-luminosity dwarfs are more affected by EFE
        is_dwarf = p['L36'] < 0.5  # < 5e8 Lsun

        env[name] = {
            'class': env_class,
            'g_ext': g_ext_estimate,
            'g_ext_over_a0': g_ext_estimate / a0_fiducial,
            'is_dwarf': is_dwarf,
            'D': p['D'],
            'L36': p['L36'],
            'Vflat': p['Vflat'],
            'Q': p['Q'],
            'T': p['T'],
        }

    return env


# ============================================================
# FITTING
# ============================================================
def fit_galaxy_no_efe(galaxy, nu_func, a0_bounds=(-10.5, -9.5)):
    """Fit galaxy without EFE. Free params: a0, Yd, (Yb)."""

    def objective(params):
        log_a0, log_Yd = params[0], params[1]
        log_Yb = params[2] if galaxy.has_bulge else log_Yd

        a0 = 10**log_a0
        Yd = 10**log_Yd
        Yb = 10**log_Yb

        gbar = galaxy.compute_gbar(Yd, Yb)
        mask = gbar > 0
        if np.sum(mask) < 3:
            return 1e10

        y = gbar[mask] / a0
        y = np.maximum(y, 1e-10)

        try:
            nu = nu_func(y)
        except:
            return 1e10

        if not np.all(np.isfinite(nu)):
            return 1e10

        gobs_pred = gbar[mask] * nu
        v_pred = np.sqrt(np.abs(gobs_pred) * galaxy.rad_m[mask]) / 1e3

        chi2 = np.sum(((galaxy.vobs[mask] - v_pred) / galaxy.errv_safe[mask])**2)
        return chi2

    bounds = [a0_bounds, (-0.7, 1.2)]
    if galaxy.has_bulge:
        bounds.append((-0.7, 1.2))

    try:
        res = differential_evolution(objective, bounds, seed=42, maxiter=150, tol=1e-6, polish=True)
        if not np.isfinite(res.fun):
            return None

        nparams = len(bounds)
        ndof = galaxy.npts - nparams

        return {
            'chi2': res.fun,
            'chi2_red': res.fun / ndof if ndof > 0 else np.inf,
            'ndof': ndof,
            'a0': 10**res.x[0],
            'Yd': 10**res.x[1],
            'Yb': 10**res.x[2] if galaxy.has_bulge else 10**res.x[1],
            'params': res.x,
        }
    except:
        return None


def fit_galaxy_with_efe(galaxy, nu_func, g_ext_range=(0, 0.5e-10)):
    """Fit galaxy WITH external field. Free params: a0, Yd, (Yb), g_ext."""

    def objective(params):
        log_a0, log_Yd = params[0], params[1]
        idx = 2
        if galaxy.has_bulge:
            log_Yb = params[idx]; idx += 1
        else:
            log_Yb = log_Yd
        log_gext = params[idx]

        a0 = 10**log_a0
        Yd = 10**log_Yd
        Yb = 10**log_Yb
        g_ext = 10**log_gext

        gbar = galaxy.compute_gbar(Yd, Yb)
        mask = gbar > 0
        if np.sum(mask) < 3:
            return 1e10

        # EFE: use (gbar + g_ext) as argument to nu, but boost only gbar
        y_total = (gbar[mask] + g_ext) / a0
        y_total = np.maximum(y_total, 1e-10)

        try:
            nu = nu_func(y_total)
        except:
            return 1e10

        if not np.all(np.isfinite(nu)):
            return 1e10

        gobs_pred = gbar[mask] * nu
        v_pred = np.sqrt(np.abs(gobs_pred) * galaxy.rad_m[mask]) / 1e3

        chi2 = np.sum(((galaxy.vobs[mask] - v_pred) / galaxy.errv_safe[mask])**2)
        return chi2

    bounds = [(-10.5, -9.5), (-0.7, 1.2)]
    if galaxy.has_bulge:
        bounds.append((-0.7, 1.2))
    bounds.append((-14.0, -9.0))  # log10(g_ext): 1e-14 to 1e-9

    try:
        res = differential_evolution(objective, bounds, seed=42, maxiter=150, tol=1e-6, polish=True)
        if not np.isfinite(res.fun):
            return None

        nparams = len(bounds)
        ndof = galaxy.npts - nparams

        idx = 2
        if galaxy.has_bulge:
            Yb = 10**res.x[idx]; idx += 1
        else:
            Yb = 10**res.x[1]
        g_ext_fit = 10**res.x[idx]

        return {
            'chi2': res.fun,
            'chi2_red': res.fun / ndof if ndof > 0 else np.inf,
            'ndof': ndof,
            'a0': 10**res.x[0],
            'Yd': 10**res.x[1],
            'Yb': Yb,
            'g_ext': g_ext_fit,
            'g_ext_over_a0': g_ext_fit / (10**res.x[0]),
            'params': res.x,
        }
    except:
        return None


def fit_galaxy_fixed_a0(galaxy, nu_func, a0, g_ext=0.0):
    """Fit with fixed a0 and g_ext, only Yd/Yb free."""

    def objective(params):
        log_Yd = params[0]
        log_Yb = params[1] if galaxy.has_bulge else log_Yd

        Yd = 10**log_Yd
        Yb = 10**log_Yb

        gbar = galaxy.compute_gbar(Yd, Yb)
        mask = gbar > 0
        if np.sum(mask) < 3:
            return 1e10

        y = (gbar[mask] + g_ext) / a0
        y = np.maximum(y, 1e-10)

        try:
            nu = nu_func(y)
        except:
            return 1e10

        if not np.all(np.isfinite(nu)):
            return 1e10

        gobs_pred = gbar[mask] * nu
        v_pred = np.sqrt(np.abs(gobs_pred) * galaxy.rad_m[mask]) / 1e3

        return np.sum(((galaxy.vobs[mask] - v_pred) / galaxy.errv_safe[mask])**2)

    bounds = [(-0.7, 1.2)]
    if galaxy.has_bulge:
        bounds.append((-0.7, 1.2))

    try:
        res = minimize(objective, [0.0] * len(bounds), method='L-BFGS-B', bounds=bounds)
        if not np.isfinite(res.fun):
            return None
        ndof = galaxy.npts - len(bounds)
        return {
            'chi2': res.fun,
            'chi2_red': res.fun / ndof if ndof > 0 else np.inf,
            'ndof': ndof,
            'Yd': 10**res.x[0],
            'Yb': 10**res.x[1] if galaxy.has_bulge else 10**res.x[0],
        }
    except:
        return None


# ============================================================
# MAIN
# ============================================================
def main():
    print("=" * 80)
    print("gs10: ENVIRONMENT & EXTERNAL FIELD EFFECT ANALYSIS")
    print("=" * 80)

    base_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(base_dir, 'Rotmod_LTG')
    master_file = os.path.join(base_dir, 'SPARC_Lelli2016c.mrt')

    # Load master table
    print("\nLoading master table...")
    props = load_master_table(master_file)
    print(f"  Loaded properties for {len(props)} galaxies")

    # Classify environment
    env = classify_environment(props)

    # Count by class
    class_counts = {}
    for name, e in env.items():
        c = e['class']
        class_counts[c] = class_counts.get(c, 0) + 1

    print("\nEnvironment classification:")
    for c, n in sorted(class_counts.items()):
        print(f"  {c}: {n}")

    # Load rotation curves
    print("\nLoading rotation curves...")
    galaxies = []
    for fname in sorted(os.listdir(data_dir)):
        if fname.endswith('_rotmod.dat'):
            try:
                gal = Galaxy(os.path.join(data_dir, fname))
                if gal.npts >= 5 and gal.name in env:
                    galaxies.append(gal)
            except:
                pass
    print(f"  Loaded {len(galaxies)} galaxies with >= 5 points and properties")

    total_pts = sum(g.npts for g in galaxies)
    print(f"  Total data points: {total_pts}")

    # ============================================================
    # PHASE 1: Fit ALL galaxies WITHOUT EFE (baseline)
    # ============================================================
    print("\n" + "=" * 80)
    print("PHASE 1: BASELINE FITS (no external field)")
    print("=" * 80)

    results_no_efe = {}

    for model_name, nu_func in [('McGaugh', nu_mcgaugh), ('MOND_simple', nu_mond_simple)]:
        print(f"\n--- {model_name} ---")

        for i, gal in enumerate(galaxies):
            if (i + 1) % 50 == 0:
                print(f"  Fitting {i+1}/{len(galaxies)}...")

            res = fit_galaxy_no_efe(gal, nu_func)
            if res is not None:
                results_no_efe[(gal.name, model_name)] = res

        fitted = sum(1 for g in galaxies if (g.name, model_name) in results_no_efe)
        print(f"  Fitted: {fitted}/{len(galaxies)}")

        chi2_reds = [results_no_efe[(g.name, model_name)]['chi2_red']
                     for g in galaxies if (g.name, model_name) in results_no_efe]
        print(f"  chi2_red: median={np.median(chi2_reds):.3f}, mean={np.mean(chi2_reds):.3f}")

        a0s = [results_no_efe[(g.name, model_name)]['a0']
               for g in galaxies if (g.name, model_name) in results_no_efe]
        print(f"  a0 (x1e-10): median={np.median(a0s)*1e10:.3f}, std={np.std(a0s)*1e10:.3f}")

    # ============================================================
    # PHASE 2: Residuals vs environment
    # ============================================================
    print("\n" + "=" * 80)
    print("PHASE 2: RESIDUALS vs ENVIRONMENT (McGaugh baseline)")
    print("=" * 80)

    model_name = 'McGaugh'

    # Collect residuals by environment class
    env_residuals = {}
    env_a0s = {}
    env_chi2s = {}

    for gal in galaxies:
        key = (gal.name, model_name)
        if key not in results_no_efe:
            continue

        res = results_no_efe[key]
        e = env[gal.name]
        ec = e['class']

        if ec not in env_residuals:
            env_residuals[ec] = []
            env_a0s[ec] = []
            env_chi2s[ec] = []

        env_a0s[ec].append(res['a0'])
        env_chi2s[ec].append(res['chi2_red'])

    print(f"\n{'Environment':<18} {'N':<5} {'median a0':<14} {'std a0':<12} {'median chi2r':<14} {'mean chi2r':<12}")
    print("-" * 75)

    for ec in sorted(env_residuals.keys()):
        if len(env_a0s.get(ec, [])) == 0:
            continue
        a = np.array(env_a0s[ec])
        c = np.array(env_chi2s[ec])
        print(f"{ec:<18} {len(a):<5} {np.median(a)*1e10:<14.4f} {np.std(a)*1e10:<12.4f} "
              f"{np.median(c):<14.4f} {np.mean(c):<12.4f}")

    # Dwarfs vs non-dwarfs
    print("\n--- Dwarfs vs massive galaxies ---")
    for label, is_dwarf_val in [('Dwarf (L<5e8)', True), ('Massive (L>5e8)', False)]:
        a0s = [results_no_efe[(g.name, model_name)]['a0']
               for g in galaxies
               if (g.name, model_name) in results_no_efe and env[g.name]['is_dwarf'] == is_dwarf_val]
        chi2s = [results_no_efe[(g.name, model_name)]['chi2_red']
                 for g in galaxies
                 if (g.name, model_name) in results_no_efe and env[g.name]['is_dwarf'] == is_dwarf_val]
        if a0s:
            print(f"  {label}: N={len(a0s)}, median a0={np.median(a0s)*1e10:.4f}, "
                  f"std a0={np.std(a0s)*1e10:.4f}, median chi2r={np.median(chi2s):.3f}")

    # Quality flag comparison
    print("\n--- By SPARC quality flag ---")
    for q in [1, 2, 3]:
        a0s = [results_no_efe[(g.name, model_name)]['a0']
               for g in galaxies
               if (g.name, model_name) in results_no_efe and env[g.name].get('Q') == q]
        chi2s = [results_no_efe[(g.name, model_name)]['chi2_red']
                 for g in galaxies
                 if (g.name, model_name) in results_no_efe and env[g.name].get('Q') == q]
        if a0s:
            print(f"  Q={q}: N={len(a0s)}, median a0={np.median(a0s)*1e10:.4f}, "
                  f"std a0={np.std(a0s)*1e10:.4f}, median chi2r={np.median(chi2s):.3f}")

    # ============================================================
    # PHASE 3: Top 20 worst fits — are they satellites?
    # ============================================================
    print("\n" + "=" * 80)
    print("PHASE 3: WORST-FIT GALAXIES — ENVIRONMENT CHECK")
    print("=" * 80)

    all_fits = []
    for gal in galaxies:
        key = (gal.name, model_name)
        if key in results_no_efe:
            res = results_no_efe[key]
            e = env[gal.name]
            all_fits.append({
                'name': gal.name,
                'chi2_red': res['chi2_red'],
                'a0': res['a0'],
                'Yd': res['Yd'],
                'env': e['class'],
                'D': e['D'],
                'L36': e['L36'],
                'Vflat': e['Vflat'],
                'g_ext_est': e['g_ext'],
                'is_dwarf': e['is_dwarf'],
                'Q': e['Q'],
                'npts': gal.npts,
            })

    all_fits.sort(key=lambda x: x['chi2_red'], reverse=True)

    print(f"\n--- 20 WORST-FIT galaxies ---")
    print(f"{'Galaxy':<15} {'chi2r':<8} {'a0(e-10)':<10} {'Yd':<7} {'Env':<15} {'D(Mpc)':<8} {'L(1e9)':<8} {'Q':<3}")
    print("-" * 74)
    for f in all_fits[:20]:
        print(f"{f['name']:<15} {f['chi2_red']:<8.2f} {f['a0']*1e10:<10.3f} {f['Yd']:<7.2f} "
              f"{f['env']:<15} {f['D']:<8.1f} {f['L36']:<8.3f} {f['Q']:<3}")

    print(f"\n--- 20 BEST-FIT galaxies ---")
    for f in all_fits[-20:]:
        print(f"{f['name']:<15} {f['chi2_red']:<8.2f} {f['a0']*1e10:<10.3f} {f['Yd']:<7.2f} "
              f"{f['env']:<15} {f['D']:<8.1f} {f['L36']:<8.3f} {f['Q']:<3}")

    # Environment breakdown for worst vs best
    print("\n--- Environment breakdown: worst 30 vs best 30 ---")
    worst30 = all_fits[:30]
    best30 = all_fits[-30:]

    for label, subset in [('Worst 30', worst30), ('Best 30', best30)]:
        if not subset:
            print(f"  {label}: (empty)")
            continue
        counts = {}
        for f in subset:
            counts[f['env']] = counts.get(f['env'], 0) + 1
        print(f"  {label}: {dict(sorted(counts.items()))}")
        dwarf_frac = sum(1 for f in subset if f['is_dwarf']) / len(subset)
        print(f"    Dwarf fraction: {dwarf_frac:.1%}")
        print(f"    Median D: {np.median([f['D'] for f in subset]):.1f} Mpc")
        print(f"    Median Q: {np.median([f['Q'] for f in subset]):.0f}")

    # ============================================================
    # PHASE 4: FIT WITH EXTERNAL FIELD (free g_ext per galaxy)
    # ============================================================
    print("\n" + "=" * 80)
    print("PHASE 4: FITS WITH FREE EXTERNAL FIELD g_ext")
    print("=" * 80)

    results_with_efe = {}

    print("\nFitting all galaxies with free g_ext (McGaugh + EFE)...")
    for i, gal in enumerate(galaxies):
        if (i + 1) % 50 == 0:
            print(f"  Fitting {i+1}/{len(galaxies)}...")

        res = fit_galaxy_with_efe(gal, nu_mcgaugh)
        if res is not None:
            results_with_efe[gal.name] = res

    print(f"  Fitted: {len(results_with_efe)}/{len(galaxies)}")

    # Compare: chi2 improvement from adding g_ext
    print(f"\n--- chi2 improvement from EFE ---")

    improvements = []
    for gal in galaxies:
        key_no = (gal.name, 'McGaugh')
        if key_no in results_no_efe and gal.name in results_with_efe:
            chi2_no = results_no_efe[key_no]['chi2']
            chi2_efe = results_with_efe[gal.name]['chi2']
            delta = chi2_no - chi2_efe  # positive = EFE better
            ndof_diff = 1  # one extra parameter

            improvements.append({
                'name': gal.name,
                'delta_chi2': delta,
                'chi2_no': chi2_no,
                'chi2_efe': chi2_efe,
                'g_ext': results_with_efe[gal.name]['g_ext'],
                'g_ext_over_a0': results_with_efe[gal.name]['g_ext_over_a0'],
                'env': env[gal.name]['class'],
                'is_dwarf': env[gal.name]['is_dwarf'],
                'D': env[gal.name]['D'],
                'L36': env[gal.name]['L36'],
                'npts': gal.npts,
            })

    improvements.sort(key=lambda x: x['delta_chi2'], reverse=True)

    print(f"\n{'Galaxy':<15} {'dchi2':<8} {'g_ext/a0':<10} {'Env':<15} {'D':<8} {'Dwarf':<6} {'Npts':<5}")
    print("-" * 67)
    print("Top 15 galaxies where EFE helps most:")
    for imp in improvements[:15]:
        print(f"{imp['name']:<15} {imp['delta_chi2']:<8.2f} {imp['g_ext_over_a0']:<10.4f} "
              f"{imp['env']:<15} {imp['D']:<8.1f} {'Y' if imp['is_dwarf'] else 'N':<6} {imp['npts']:<5}")

    print("\nTop 15 galaxies where EFE helps least (or hurts):")
    for imp in improvements[-15:]:
        print(f"{imp['name']:<15} {imp['delta_chi2']:<8.2f} {imp['g_ext_over_a0']:<10.4f} "
              f"{imp['env']:<15} {imp['D']:<8.1f} {'Y' if imp['is_dwarf'] else 'N':<6} {imp['npts']:<5}")

    # Statistics by environment class
    print(f"\n--- Mean delta-chi2 by environment ---")
    for ec in sorted(set(e['env'] for e in improvements)):
        subset = [x for x in improvements if x['env'] == ec]
        if subset:
            deltas = [x['delta_chi2'] for x in subset]
            g_exts = [x['g_ext_over_a0'] for x in subset]
            sig_improve = sum(1 for d in deltas if d > 4.0)  # > 2sigma for 1 dof
            print(f"  {ec:<18}: N={len(subset)}, mean dchi2={np.mean(deltas):.2f}, "
                  f"median g_ext/a0={np.median(g_exts):.4f}, "
                  f"significant improve: {sig_improve}/{len(subset)}")

    # By dwarf status
    print(f"\n--- Mean delta-chi2: dwarfs vs massive ---")
    for label, dw in [('Dwarfs', True), ('Massive', False)]:
        subset = [x for x in improvements if x['is_dwarf'] == dw]
        if subset:
            deltas = [x['delta_chi2'] for x in subset]
            g_exts = [x['g_ext_over_a0'] for x in subset]
            sig = sum(1 for d in deltas if d > 4.0)
            print(f"  {label}: N={len(subset)}, mean dchi2={np.mean(deltas):.2f}, "
                  f"median g_ext/a0={np.median(g_exts):.4f}, significant: {sig}/{len(subset)}")

    # ============================================================
    # PHASE 5: GLOBAL FIT — fixed a0, with/without environment-estimated g_ext
    # ============================================================
    print("\n" + "=" * 80)
    print("PHASE 5: GLOBAL FIT — fixed a0, estimated vs zero g_ext")
    print("=" * 80)

    a0_test = 1.2e-10

    total_chi2_no_efe = 0
    total_chi2_est_efe = 0
    total_ndof = 0
    n_fitted = 0

    for gal in galaxies:
        if gal.name not in env:
            continue

        # Without EFE
        res_no = fit_galaxy_fixed_a0(gal, nu_mcgaugh, a0_test, g_ext=0.0)
        # With estimated EFE
        g_ext_est = env[gal.name]['g_ext']
        res_efe = fit_galaxy_fixed_a0(gal, nu_mcgaugh, a0_test, g_ext=g_ext_est)

        if res_no is not None and res_efe is not None:
            total_chi2_no_efe += res_no['chi2']
            total_chi2_est_efe += res_efe['chi2']
            total_ndof += res_no['ndof']
            n_fitted += 1

    print(f"\n  Galaxies fitted: {n_fitted}")
    if total_ndof > 0:
        print(f"  Total chi2 (no EFE):         {total_chi2_no_efe:.1f}  (chi2/ndof = {total_chi2_no_efe/total_ndof:.4f})")
        print(f"  Total chi2 (estimated EFE):  {total_chi2_est_efe:.1f}  (chi2/ndof = {total_chi2_est_efe/total_ndof:.4f})")
        print(f"  Delta chi2:                  {total_chi2_no_efe - total_chi2_est_efe:.1f}")
    else:
        print(f"  No valid fits obtained.")
    print(f"  (positive = EFE helps)")

    # By environment
    print(f"\n--- Delta chi2 by environment (est. g_ext) ---")
    for ec in sorted(class_counts.keys()):
        chi2_no = 0; chi2_efe = 0; ndof = 0; n = 0
        for gal in galaxies:
            if gal.name not in env or env[gal.name]['class'] != ec:
                continue
            g_ext_est = env[gal.name]['g_ext']
            res_no = fit_galaxy_fixed_a0(gal, nu_mcgaugh, a0_test, g_ext=0.0)
            res_efe = fit_galaxy_fixed_a0(gal, nu_mcgaugh, a0_test, g_ext=g_ext_est)
            if res_no and res_efe:
                chi2_no += res_no['chi2']
                chi2_efe += res_efe['chi2']
                ndof += res_no['ndof']
                n += 1
        if n > 0:
            delta = chi2_no - chi2_efe
            print(f"  {ec:<18}: N={n}, dchi2={delta:+.1f}, "
                  f"chi2/ndof: {chi2_no/ndof:.4f} -> {chi2_efe/ndof:.4f}")

    # ============================================================
    # PHASE 6: CORRELATION — fitted g_ext vs environment proxy
    # ============================================================
    print("\n" + "=" * 80)
    print("PHASE 6: CORRELATION — fitted g_ext vs distance & luminosity")
    print("=" * 80)

    # Collect fitted g_ext values
    gext_data = []
    for gal in galaxies:
        if gal.name in results_with_efe and gal.name in env:
            e = env[gal.name]
            r = results_with_efe[gal.name]
            gext_data.append({
                'name': gal.name,
                'g_ext_fit': r['g_ext'],
                'g_ext_over_a0': r['g_ext_over_a0'],
                'D': e['D'],
                'L36': e['L36'],
                'Vflat': e['Vflat'],
                'env': e['class'],
                'is_dwarf': e['is_dwarf'],
            })

    if gext_data:
        D_arr = np.array([x['D'] for x in gext_data])
        gext_arr = np.array([x['g_ext_fit'] for x in gext_data])
        L_arr = np.array([x['L36'] for x in gext_data])

        # Spearman correlations
        from scipy.stats import spearmanr

        rho_D, p_D = spearmanr(D_arr, gext_arr)
        print(f"\n  Spearman correlation: g_ext_fit vs Distance")
        print(f"    rho = {rho_D:.4f}, p = {p_D:.4e}")
        print(f"    {'SIGNIFICANT' if p_D < 0.05 else 'not significant'}")

        rho_L, p_L = spearmanr(L_arr, gext_arr)
        print(f"\n  Spearman correlation: g_ext_fit vs Luminosity")
        print(f"    rho = {rho_L:.4f}, p = {p_L:.4e}")
        print(f"    {'SIGNIFICANT' if p_L < 0.05 else 'not significant'}")

        # Inverse distance proxy: g_ext ~ 1/D^2 (tidal field from MW)
        invD2 = 1.0 / D_arr**2
        rho_invD2, p_invD2 = spearmanr(invD2, gext_arr)
        print(f"\n  Spearman correlation: g_ext_fit vs 1/D^2 (tidal proxy)")
        print(f"    rho = {rho_invD2:.4f}, p = {p_invD2:.4e}")
        print(f"    {'SIGNIFICANT' if p_invD2 < 0.05 else 'not significant'}")

        # Vflat as mass proxy
        vflat_arr = np.array([x['Vflat'] for x in gext_data])
        mask_v = vflat_arr > 0
        if np.sum(mask_v) > 10:
            rho_V, p_V = spearmanr(vflat_arr[mask_v], gext_arr[mask_v])
            print(f"\n  Spearman correlation: g_ext_fit vs Vflat")
            print(f"    rho = {rho_V:.4f}, p = {p_V:.4e}")
            print(f"    {'SIGNIFICANT' if p_V < 0.05 else 'not significant'}")

        # g_ext distribution by environment
        print(f"\n--- g_ext/a0 distribution by environment ---")
        for ec in sorted(set(x['env'] for x in gext_data)):
            subset = [x['g_ext_over_a0'] for x in gext_data if x['env'] == ec]
            if subset:
                print(f"  {ec:<18}: N={len(subset)}, median={np.median(subset):.4f}, "
                      f"mean={np.mean(subset):.4f}, std={np.std(subset):.4f}")

    # ============================================================
    # SUMMARY
    # ============================================================
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print("""
Key questions answered:
1. Do worst-fit galaxies concentrate in non-isolated environments?
2. Does adding g_ext improve fits preferentially for non-isolated galaxies?
3. Does fitted g_ext correlate with environment proxies (distance, luminosity)?
4. What is the sign and scale of the environmental effect?

If environment matters:
  - Satellites should have systematically different residuals
  - g_ext should correlate with 1/D^2 or cluster membership
  - Dwarfs in dense environments should need larger g_ext
  - This would support EFE as part of the physics, not just a detail
""")

    print("=" * 80)
    print("DONE")
    print("=" * 80)


if __name__ == '__main__':
    main()
