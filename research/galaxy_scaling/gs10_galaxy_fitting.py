#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gs10: HARD NUMERICS — Galaxy-by-galaxy rotation curve fitting
=============================================================

Strategy:
  1. Load ALL 175 SPARC late-type galaxies (individual rotation curves)
  2. For each galaxy, compute g_bar from baryonic components (gas + disk + bulge)
  3. Fit a GENERALIZED interpolation function ν(y) with free parameters
  4. Find optimal equation form by minimizing total χ² across all galaxies
  5. Compare multiple functional forms head-to-head

Functional forms tested:
  F1: ν = 1/(1 - exp(-√y))                         [McGaugh empirical, 1 param: a0]
  F2: ν = 1/2 + √(1/4 + 1/y)                       [MOND simple, 1 param: a0]
  F3: ν = (1 + y^(-n))^(1/n)                        [Generalized MOND, 2 params: a0, n]
  F4: ν = 1 + 1/(√y + y^α)                          [Hybrid family, 2 params: a0, α]
  F5: ν = (1 + (y/β)^(-n))^(1/n)                    [Rescaled gen., 3 params: a0, n, β]
  F6: ν = 1 + exp(-y^α) / y^γ                       [Exponential family, 3 params: a0, α, γ]

Each galaxy also has a free mass-to-light ratio Υ_disk (and Υ_bulge if present).

Data: SPARC Rotmod_LTG (175 galaxies, Lelli+ 2016)
  Columns: Rad(kpc) Vobs(km/s) errV Vgas Vdisk Vbul SBdisk SBbul
  Vgas, Vdisk, Vbul are velocity contributions assuming Υ=1 for disk/bulge
"""

import sys
import io
import os
import numpy as np
from scipy.optimize import minimize, minimize_scalar, differential_evolution
from scipy.stats import chi2 as chi2_dist
import warnings
warnings.filterwarnings('ignore')

# Windows encoding fix
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

# ============================================================
# CONSTANTS
# ============================================================
G = 6.674e-11       # m^3 kg^-1 s^-2
c = 2.998e8          # m/s
M_sun = 1.989e30     # kg
kpc = 3.086e19       # m
H0 = 67.4e3 / (1e6 * kpc)  # 1/s

# ============================================================
# INTERPOLATION FUNCTIONS
# ============================================================
def nu_mcgaugh(y):
    """McGaugh empirical: g = gbar / (1 - exp(-sqrt(gbar/a0)))"""
    sy = np.sqrt(np.abs(y)) * np.sign(y)
    return 1.0 / (1.0 - np.exp(-sy))

def nu_mond_simple(y):
    """MOND simple: ν = 1/2 + sqrt(1/4 + 1/y)"""
    return 0.5 + np.sqrt(0.25 + 1.0/np.abs(y))

def nu_generalized(y, n):
    """Generalized: ν = (1 + y^(-n))^(1/n)"""
    return (1.0 + np.abs(y)**(-n))**(1.0/n)

def nu_hybrid(y, alpha):
    """Hybrid: ν = 1 + 1/(√y + y^α)"""
    sy = np.sqrt(np.abs(y))
    return 1.0 + 1.0/(sy + np.abs(y)**alpha)

def nu_rescaled_gen(y, n, beta):
    """Rescaled generalized: ν = (1 + (y/β)^(-n))^(1/n)"""
    yb = np.abs(y) / beta
    return (1.0 + yb**(-n))**(1.0/n)

def nu_exponential(y, alpha, gamma):
    """Exponential family: ν = 1 + exp(-y^α) / y^γ"""
    ya = np.abs(y)**alpha
    yg = np.abs(y)**gamma
    return 1.0 + np.exp(-ya) / yg

# ============================================================
# GALAXY DATA LOADER
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
        self.rad_kpc = data[:, 0]    # kpc
        self.vobs = data[:, 1]       # km/s
        self.errv = data[:, 2]       # km/s
        self.vgas = data[:, 3]       # km/s (signed — can be negative)
        self.vdisk = data[:, 4]      # km/s (at Υ_disk = 1)
        self.vbul = data[:, 5]       # km/s (at Υ_bulge = 1)
        self.sbdisk = data[:, 6]     # L/pc²
        self.sbbul = data[:, 7]      # L/pc²

        self.npts = len(self.rad_kpc)
        self.has_bulge = np.any(self.vbul > 0)

        # Minimum error floor: 5 km/s or 5% of Vobs
        self.errv_safe = np.maximum(self.errv, np.maximum(5.0, 0.05 * np.abs(self.vobs)))

        # Precompute radii in meters
        self.rad_m = self.rad_kpc * kpc

    def compute_gbar(self, Yd, Yb=None):
        """Compute baryonic acceleration given mass-to-light ratios.

        g_bar = V²_bar / r, where V²_bar = |Vgas|*Vgas + Yd * Vdisk² + Yb * Vbul²

        Vgas is kept with sign (but squared with sign preservation: |Vgas|*Vgas)
        Vdisk, Vbul are always positive (velocity contributions)
        """
        if Yb is None:
            Yb = Yd  # same M/L for disk and bulge if not specified

        # V²_bar = Vgas² (with sign) + Yd * Vdisk² + Yb * Vbul²
        v2_bar = np.abs(self.vgas) * self.vgas + Yd * self.vdisk**2 + Yb * self.vbul**2

        # Convert to acceleration: g = V² / r
        # V in km/s → m/s: ×1e3; r in kpc → m
        v2_bar_si = v2_bar * 1e6  # (km/s)² → (m/s)²

        gbar = v2_bar_si / self.rad_m
        return gbar

    def compute_gobs(self):
        """Observed centripetal acceleration: g_obs = V²_obs / r"""
        return (self.vobs * 1e3)**2 / self.rad_m


def load_all_galaxies(data_dir):
    """Load all SPARC galaxies from directory."""
    galaxies = []
    for fname in sorted(os.listdir(data_dir)):
        if fname.endswith('_rotmod.dat'):
            try:
                gal = Galaxy(os.path.join(data_dir, fname))
                if gal.npts >= 5:  # minimum 5 data points for fitting
                    galaxies.append(gal)
            except Exception as e:
                print(f"  SKIP {fname}: {e}")
    return galaxies


# ============================================================
# FITTING ENGINE
# ============================================================
def chi2_galaxy(params, galaxy, nu_func, nu_nparams):
    """Compute χ² for a single galaxy given interpolation function.

    params = [log10(a0), log10(Yd), (log10(Yb) if bulge), (nu_params...)]
    """
    idx = 0
    log_a0 = params[idx]; idx += 1
    log_Yd = params[idx]; idx += 1
    if galaxy.has_bulge:
        log_Yb = params[idx]; idx += 1
    else:
        log_Yb = log_Yd
    nu_params = params[idx:]

    a0 = 10**log_a0
    Yd = 10**log_Yd
    Yb = 10**log_Yb

    gbar = galaxy.compute_gbar(Yd, Yb)
    gobs = galaxy.compute_gobs()

    # Handle negative gbar (counterrotating gas)
    mask = gbar > 0
    if np.sum(mask) < 3:
        return 1e10

    y = gbar[mask] / a0

    # Ensure y is positive and finite
    y = np.maximum(y, 1e-10)

    try:
        if nu_nparams == 0:
            nu = nu_func(y)
        elif nu_nparams == 1:
            nu = nu_func(y, nu_params[0])
        elif nu_nparams == 2:
            nu = nu_func(y, nu_params[0], nu_params[1])
        else:
            return 1e10
    except (FloatingPointError, OverflowError, ValueError):
        return 1e10

    if not np.all(np.isfinite(nu)):
        return 1e10

    # Predicted velocity: V²_pred = gobs_pred * r = gbar * ν * r
    # But actually: gobs_pred = gbar * ν
    # V_pred = sqrt(gobs_pred * r) = sqrt(gbar * ν * r)
    gobs_pred = gbar[mask] * nu
    v_pred = np.sqrt(np.abs(gobs_pred) * galaxy.rad_m[mask]) / 1e3  # back to km/s

    # χ² on velocities
    residuals = (galaxy.vobs[mask] - v_pred) / galaxy.errv_safe[mask]
    chi2 = np.sum(residuals**2)

    return chi2


def fit_galaxy(galaxy, nu_func, nu_nparams, nu_param_bounds=None, a0_init=-9.9):
    """Fit a single galaxy with given interpolation function.

    Returns: dict with best-fit params, chi2, reduced chi2, etc.
    """
    # Parameter bounds
    bounds = [
        (-11.0, -9.0),    # log10(a0): 1e-11 to 1e-9
        (-0.7, 1.2),      # log10(Yd): 0.2 to 15 (M/L in 3.6μm band)
    ]
    x0 = [a0_init, 0.0]   # a0 ~ 1.2e-10, Yd ~ 1.0

    if galaxy.has_bulge:
        bounds.append((-0.7, 1.2))  # log10(Yb)
        x0.append(0.0)

    if nu_param_bounds is not None:
        bounds.extend(nu_param_bounds)
        # Initial values at midpoint of bounds
        for b in nu_param_bounds:
            x0.append(0.5 * (b[0] + b[1]))

    ndata = galaxy.npts
    nparams = len(x0)
    ndof = ndata - nparams

    if ndof < 1:
        return None

    # Optimize with differential evolution (global) then polish with L-BFGS-B
    try:
        result_de = differential_evolution(
            chi2_galaxy, bounds,
            args=(galaxy, nu_func, nu_nparams),
            seed=42, maxiter=200, tol=1e-6,
            polish=True
        )
        best_chi2 = result_de.fun
        best_params = result_de.x
    except Exception:
        return None

    if not np.isfinite(best_chi2):
        return None

    # Extract parameters
    idx = 0
    log_a0 = best_params[idx]; idx += 1
    log_Yd = best_params[idx]; idx += 1
    log_Yb = best_params[idx] if galaxy.has_bulge else log_Yd
    if galaxy.has_bulge:
        idx += 1
    nu_params_fit = best_params[idx:]

    return {
        'name': galaxy.name,
        'npts': ndata,
        'ndof': ndof,
        'chi2': best_chi2,
        'chi2_red': best_chi2 / ndof if ndof > 0 else np.inf,
        'log_a0': log_a0,
        'a0': 10**log_a0,
        'Yd': 10**log_Yd,
        'Yb': 10**log_Yb,
        'nu_params': nu_params_fit,
        'has_bulge': galaxy.has_bulge,
        'params': best_params,
    }


def fit_galaxy_fixed_a0(galaxy, nu_func, nu_nparams, a0, nu_param_values=None):
    """Fit galaxy with FIXED a0 and nu_params — only Yd (and Yb) free.

    This is for the global fit: fix equation, optimize only mass-to-light.
    """
    bounds = [(-0.7, 1.2)]  # log10(Yd)
    x0 = [0.0]

    if galaxy.has_bulge:
        bounds.append((-0.7, 1.2))
        x0.append(0.0)

    def objective(ml_params):
        full_params = [np.log10(a0)] + list(ml_params)
        if nu_param_values is not None:
            full_params = full_params + list(nu_param_values)
        return chi2_galaxy(full_params, galaxy, nu_func, nu_nparams)

    try:
        result = minimize(objective, x0, method='L-BFGS-B', bounds=bounds)
        best_chi2 = result.fun
        best_ml = result.x
    except Exception:
        return None

    if not np.isfinite(best_chi2):
        return None

    ndata = galaxy.npts
    nparams = len(x0)
    ndof = ndata - nparams

    return {
        'name': galaxy.name,
        'npts': ndata,
        'ndof': ndof,
        'chi2': best_chi2,
        'chi2_red': best_chi2 / ndof if ndof > 0 else np.inf,
        'a0': a0,
        'Yd': 10**best_ml[0],
        'Yb': 10**best_ml[1] if galaxy.has_bulge else 10**best_ml[0],
        'has_bulge': galaxy.has_bulge,
    }


# ============================================================
# GLOBAL FIT: Fix equation params, fit only M/L per galaxy
# ============================================================
def global_chi2(eq_params, galaxies, nu_func, nu_nparams, nu_param_indices):
    """Total χ² across ALL galaxies for given equation parameters.

    eq_params = [log10(a0), (nu_shape_params...)]
    Each galaxy gets its own Yd (and Yb) fitted internally.
    """
    log_a0 = eq_params[0]
    a0 = 10**log_a0
    nu_params = eq_params[1:]

    total_chi2 = 0.0
    for gal in galaxies:
        result = fit_galaxy_fixed_a0(gal, nu_func, nu_nparams, a0,
                                      nu_param_values=nu_params if len(nu_params) > 0 else None)
        if result is None:
            total_chi2 += 1e6
        else:
            total_chi2 += result['chi2']

    return total_chi2


# ============================================================
# MODEL DEFINITIONS
# ============================================================
MODELS = {
    'McGaugh': {
        'func': nu_mcgaugh,
        'nparams': 0,
        'bounds': [],       # no shape params
        'label': 'g/(1-exp(-sqrt(g/a0)))',
    },
    'MOND_simple': {
        'func': nu_mond_simple,
        'nparams': 0,
        'bounds': [],
        'label': '1/2 + sqrt(1/4 + a0/g)',
    },
    'Generalized_n': {
        'func': nu_generalized,
        'nparams': 1,
        'bounds': [(0.3, 3.0)],   # n
        'label': '(1 + (g/a0)^(-n))^(1/n)',
    },
    'Hybrid_alpha': {
        'func': nu_hybrid,
        'nparams': 1,
        'bounds': [(0.5, 2.0)],   # alpha
        'label': '1 + 1/(sqrt(y) + y^alpha)',
    },
    'Rescaled_gen': {
        'func': nu_rescaled_gen,
        'nparams': 2,
        'bounds': [(0.3, 3.0), (0.3, 3.0)],  # n, beta
        'label': '(1 + (y/beta)^(-n))^(1/n)',
    },
    'Exponential': {
        'func': nu_exponential,
        'nparams': 2,
        'bounds': [(0.1, 2.0), (0.1, 1.5)],  # alpha, gamma
        'label': '1 + exp(-y^alpha)/y^gamma',
    },
}


# ============================================================
# MAIN
# ============================================================
def main():
    print("=" * 80)
    print("gs10: HARD NUMERICS — Galaxy-by-galaxy rotation curve fitting")
    print("=" * 80)

    # Load data
    data_dir = os.path.join(os.path.dirname(__file__), 'Rotmod_LTG')
    print(f"\nLoading galaxies from {data_dir}...")
    galaxies = load_all_galaxies(data_dir)
    print(f"Loaded {len(galaxies)} galaxies (>= 5 data points)")

    total_pts = sum(g.npts for g in galaxies)
    print(f"Total data points: {total_pts}")
    print(f"Galaxies with bulge: {sum(1 for g in galaxies if g.has_bulge)}")

    # ============================================================
    # PHASE 1: Galaxy-by-galaxy free fits (all params free per galaxy)
    # ============================================================
    print("\n" + "=" * 80)
    print("PHASE 1: Individual galaxy fits (a0, Yd, Yb, shape params — all free)")
    print("=" * 80)

    phase1_results = {}

    for model_name, model_info in MODELS.items():
        print(f"\n--- Model: {model_name} [{model_info['label']}] ---")

        results = []
        failed = 0

        for i, gal in enumerate(galaxies):
            if (i + 1) % 25 == 0:
                print(f"  fitting galaxy {i+1}/{len(galaxies)}...")

            res = fit_galaxy(
                gal, model_info['func'], model_info['nparams'],
                nu_param_bounds=model_info['bounds'] if model_info['bounds'] else None
            )

            if res is not None:
                results.append(res)
            else:
                failed += 1

        phase1_results[model_name] = results

        if results:
            chi2_reds = [r['chi2_red'] for r in results]
            a0s = [r['a0'] for r in results]
            yds = [r['Yd'] for r in results]

            print(f"  Fitted: {len(results)}/{len(galaxies)} (failed: {failed})")
            print(f"  chi2_red: median={np.median(chi2_reds):.3f}, mean={np.mean(chi2_reds):.3f}")
            print(f"  a0 (×1e-10): median={np.median(a0s)*1e10:.3f}, "
                  f"mean={np.mean(a0s)*1e10:.3f}, std={np.std(a0s)*1e10:.3f}")
            print(f"  Yd: median={np.median(yds):.3f}, mean={np.mean(yds):.3f}")

            # Total chi2
            total_chi2 = sum(r['chi2'] for r in results)
            total_ndof = sum(r['ndof'] for r in results)
            print(f"  Total chi2/ndof = {total_chi2:.1f}/{total_ndof} = {total_chi2/total_ndof:.4f}")

            if model_info['nparams'] > 0:
                shape_params = np.array([r['nu_params'] for r in results])
                for j in range(model_info['nparams']):
                    pname = ['n', 'alpha', 'beta', 'gamma'][j] if j < 4 else f'p{j}'
                    print(f"  {pname}: median={np.median(shape_params[:,j]):.4f}, "
                          f"mean={np.mean(shape_params[:,j]):.4f}, "
                          f"std={np.std(shape_params[:,j]):.4f}")

    # ============================================================
    # PHASE 2: Global fit — fix equation params, fit only M/L per galaxy
    # ============================================================
    print("\n" + "=" * 80)
    print("PHASE 2: Global fit — shared a0 + shape params, per-galaxy Yd/Yb")
    print("=" * 80)

    phase2_results = {}

    for model_name, model_info in MODELS.items():
        print(f"\n--- Model: {model_name} ---")

        # Bounds for global params: [log10(a0), shape_params...]
        global_bounds = [(-10.5, -9.5)]  # log10(a0)
        global_bounds.extend(model_info['bounds'])

        n_global = 1 + model_info['nparams']

        print(f"  Global params: {n_global} (a0 + {model_info['nparams']} shape)")
        print(f"  Per-galaxy params: Yd (+ Yb if bulge)")
        print(f"  Running differential evolution...")

        try:
            result_de = differential_evolution(
                global_chi2, global_bounds,
                args=(galaxies, model_info['func'], model_info['nparams'], None),
                seed=42, maxiter=100, tol=1e-4,
                polish=True, workers=1
            )

            best_global = result_de.x
            best_total_chi2 = result_de.fun

            a0_best = 10**best_global[0]
            shape_best = best_global[1:]

            print(f"  Best a0 = {a0_best:.4e} ({a0_best/1.2e-10:.4f} × 1.2e-10)")
            if len(shape_best) > 0:
                for j, sp in enumerate(shape_best):
                    pname = ['n', 'alpha', 'beta', 'gamma'][j] if j < 4 else f'p{j}'
                    print(f"  Best {pname} = {sp:.6f}")

            # Now get per-galaxy results with best global params
            per_galaxy = []
            for gal in galaxies:
                res = fit_galaxy_fixed_a0(
                    gal, model_info['func'], model_info['nparams'],
                    a0_best,
                    nu_param_values=shape_best if len(shape_best) > 0 else None
                )
                if res is not None:
                    per_galaxy.append(res)

            total_ndof_global = sum(r['ndof'] for r in per_galaxy) - n_global
            print(f"  Total chi2 = {best_total_chi2:.1f}")
            print(f"  Total ndof = {total_ndof_global}")
            print(f"  chi2/ndof = {best_total_chi2/total_ndof_global:.4f}")

            # Yd distribution
            yds = [r['Yd'] for r in per_galaxy]
            chi2_reds = [r['chi2_red'] for r in per_galaxy]
            print(f"  Yd: median={np.median(yds):.3f}, mean={np.mean(yds):.3f}, std={np.std(yds):.3f}")
            print(f"  chi2_red per galaxy: median={np.median(chi2_reds):.3f}, mean={np.mean(chi2_reds):.3f}")

            # Count "good" fits (chi2_red < 2)
            good = sum(1 for cr in chi2_reds if cr < 2.0)
            print(f"  Galaxies with chi2_red < 2: {good}/{len(per_galaxy)}")

            # BIC for model comparison
            n_total = sum(r['npts'] for r in per_galaxy)
            n_free_total = sum(2 if r['has_bulge'] else 1 for r in per_galaxy) + n_global
            bic = best_total_chi2 + n_free_total * np.log(n_total)
            print(f"  BIC = {bic:.1f} (n_free_total = {n_free_total})")

            phase2_results[model_name] = {
                'a0': a0_best,
                'shape_params': shape_best,
                'total_chi2': best_total_chi2,
                'total_ndof': total_ndof_global,
                'chi2_red': best_total_chi2 / total_ndof_global,
                'bic': bic,
                'per_galaxy': per_galaxy,
                'n_good': good,
                'n_galaxies': len(per_galaxy),
            }

        except Exception as e:
            print(f"  FAILED: {e}")
            import traceback
            traceback.print_exc()

    # ============================================================
    # PHASE 3: Model comparison summary
    # ============================================================
    print("\n" + "=" * 80)
    print("PHASE 3: MODEL COMPARISON")
    print("=" * 80)

    print(f"\n{'Model':<20} {'a0 (×1e-10)':<14} {'Shape':<20} {'chi2/ndof':<12} {'BIC':<12} {'Good/N':<10}")
    print("-" * 88)

    # Sort by BIC
    sorted_models = sorted(phase2_results.items(), key=lambda x: x[1]['bic'])

    best_bic = sorted_models[0][1]['bic'] if sorted_models else 0

    for model_name, res in sorted_models:
        shape_str = ', '.join(f'{s:.3f}' for s in res['shape_params']) if len(res['shape_params']) > 0 else '—'
        delta_bic = res['bic'] - best_bic
        print(f"{model_name:<20} {res['a0']*1e10:<14.4f} {shape_str:<20} "
              f"{res['chi2_red']:<12.4f} {delta_bic:<12.1f} "
              f"{res['n_good']}/{res['n_galaxies']}")

    # ============================================================
    # PHASE 4: Worst/best galaxies per model
    # ============================================================
    print("\n" + "=" * 80)
    print("PHASE 4: GALAXY-BY-GALAXY DIAGNOSTICS (best model)")
    print("=" * 80)

    if sorted_models:
        best_model_name = sorted_models[0][0]
        best_res = sorted_models[0][1]
        print(f"\nBest model: {best_model_name}")
        print(f"a0 = {best_res['a0']:.4e}")

        pg = best_res['per_galaxy']
        pg_sorted = sorted(pg, key=lambda x: x['chi2_red'])

        print(f"\n--- 10 BEST-FIT galaxies ---")
        print(f"{'Galaxy':<20} {'Npts':<6} {'chi2_red':<10} {'Yd':<8} {'Yb':<8}")
        for r in pg_sorted[:10]:
            yb_str = f"{r['Yb']:.3f}" if r['has_bulge'] else "—"
            print(f"{r['name']:<20} {r['npts']:<6} {r['chi2_red']:<10.4f} {r['Yd']:<8.3f} {yb_str:<8}")

        print(f"\n--- 10 WORST-FIT galaxies ---")
        for r in pg_sorted[-10:]:
            yb_str = f"{r['Yb']:.3f}" if r['has_bulge'] else "—"
            print(f"{r['name']:<20} {r['npts']:<6} {r['chi2_red']:<10.4f} {r['Yd']:<8.3f} {yb_str:<8}")

        print(f"\n--- Yd distribution ---")
        yds = np.array([r['Yd'] for r in pg])
        percentiles = [5, 16, 25, 50, 75, 84, 95]
        for p in percentiles:
            print(f"  {p}th percentile: Yd = {np.percentile(yds, p):.3f}")

        # Physical reasonableness: Yd in 3.6μm band should be 0.2–1.0 for stars
        phys_ok = sum(1 for yd in yds if 0.2 <= yd <= 1.2)
        print(f"\n  Physically reasonable Yd (0.2–1.2): {phys_ok}/{len(yds)} ({100*phys_ok/len(yds):.0f}%)")

    # ============================================================
    # PHASE 5: Delta-chi2 between models for each galaxy
    # ============================================================
    if len(phase2_results) >= 2:
        print("\n" + "=" * 80)
        print("PHASE 5: PER-GALAXY MODEL PREFERENCE")
        print("=" * 80)

        # Compare top 2 models
        m1_name, m1_res = sorted_models[0]
        m2_name, m2_res = sorted_models[1]

        print(f"\nComparing {m1_name} vs {m2_name}:")

        # Build name-indexed dicts
        m1_dict = {r['name']: r['chi2'] for r in m1_res['per_galaxy']}
        m2_dict = {r['name']: r['chi2'] for r in m2_res['per_galaxy']}

        common = set(m1_dict.keys()) & set(m2_dict.keys())

        prefer_m1 = 0
        prefer_m2 = 0
        delta_list = []

        for name in common:
            delta = m1_dict[name] - m2_dict[name]  # negative = m1 better
            delta_list.append((name, delta))
            if delta < -2:  # significant preference
                prefer_m1 += 1
            elif delta > 2:
                prefer_m2 += 1

        print(f"  Galaxies preferring {m1_name} (Δχ² < -2): {prefer_m1}")
        print(f"  Galaxies preferring {m2_name} (Δχ² > +2): {prefer_m2}")
        print(f"  No strong preference: {len(common) - prefer_m1 - prefer_m2}")

        delta_list.sort(key=lambda x: x[1])
        print(f"\n  Strongest preference for {m1_name}:")
        for name, d in delta_list[:5]:
            print(f"    {name}: Δχ² = {d:.2f}")
        print(f"\n  Strongest preference for {m2_name}:")
        for name, d in delta_list[-5:]:
            print(f"    {name}: Δχ² = {d:.2f}")

    print("\n" + "=" * 80)
    print("DONE — gs10 complete")
    print("=" * 80)


if __name__ == '__main__':
    main()
