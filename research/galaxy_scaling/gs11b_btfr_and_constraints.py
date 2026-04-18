#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gs11b: BTFR verification + constrained model tests
====================================================

Key questions from gs11:
  1. BTFR: does the winning equation match observed M ~ v^3.85?
     Test DIRECTLY on SPARC Vflat data (not asymptotic theory)
  2. Constrained model: gamma = alpha/2 (1 shape param) — how good?
  3. Forced gamma=0.5: how bad is it really? Per-galaxy comparison
  4. Deep MOND: is y<<1 regime even sampled by SPARC data?
"""

import sys
import io
import os
import numpy as np
from scipy.optimize import minimize, differential_evolution
from scipy.stats import spearmanr, linregress
import warnings
warnings.filterwarnings('ignore')

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

G = 6.674e-11
M_sun = 1.989e30
kpc = 3.086e19
a0_obs = 1.2e-10

# ============================================================
# GALAXY LOADER
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
        if Yb is None: Yb = Yd
        v2_bar = np.abs(self.vgas) * self.vgas + Yd * self.vdisk**2 + Yb * self.vbul**2
        return v2_bar * 1e6 / self.rad_m

def load_master_table(filepath):
    props = {}
    with open(filepath, 'r') as f:
        lines = f.readlines()
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
            props[name] = {
                'T': int(parts[1]), 'D': float(parts[2]),
                'L36': float(parts[7]),   # 10^9 Lsun
                'MHI': float(parts[13]),  # 10^9 Msun
                'Vflat': float(parts[15]),
                'e_Vflat': float(parts[16]),
                'Q': int(parts[17]),
            }
        except (ValueError, IndexError):
            continue
    return props


# ============================================================
# MODEL FUNCTIONS
# ============================================================
def nu_exp(y, alpha, gamma):
    return 1.0 + np.exp(-y**alpha) / y**gamma

def nu_exp_constrained(y, alpha):
    """gamma = alpha/2"""
    return 1.0 + np.exp(-y**alpha) / y**(alpha/2)

def nu_mcgaugh(y):
    return 1.0 / (1.0 - np.exp(-np.sqrt(np.maximum(y, 1e-30))))

def nu_mond_simple(y):
    return 0.5 + np.sqrt(0.25 + 1.0/np.maximum(y, 1e-30))


# ============================================================
# FITTING ENGINE
# ============================================================
def fit_galaxy_model(galaxy, nu_func, a0, extra_params=None):
    """Fit Yd (and Yb) for given nu_func and fixed a0."""
    def objective(params):
        log_Yd = params[0]
        log_Yb = params[1] if galaxy.has_bulge else log_Yd
        Yd = 10**log_Yd
        Yb = 10**log_Yb
        gbar = galaxy.compute_gbar(Yd, Yb)
        mask = gbar > 0
        if np.sum(mask) < 3: return 1e10
        y = gbar[mask] / a0
        y = np.maximum(y, 1e-10)
        try:
            if extra_params is not None:
                nu = nu_func(y, *extra_params)
            else:
                nu = nu_func(y)
        except: return 1e10
        if not np.all(np.isfinite(nu)): return 1e10
        gobs_pred = gbar[mask] * nu
        v_pred = np.sqrt(np.abs(gobs_pred) * galaxy.rad_m[mask]) / 1e3
        return np.sum(((galaxy.vobs[mask] - v_pred) / galaxy.errv_safe[mask])**2)

    bounds = [(-0.7, 1.2)]
    if galaxy.has_bulge:
        bounds.append((-0.7, 1.2))

    res = minimize(objective, [0.0]*len(bounds), method='L-BFGS-B', bounds=bounds)
    if not np.isfinite(res.fun): return None
    ndof = galaxy.npts - len(bounds)
    return {
        'chi2': res.fun, 'chi2_red': res.fun/ndof if ndof > 0 else np.inf,
        'ndof': ndof, 'Yd': 10**res.x[0],
        'Yb': 10**res.x[1] if galaxy.has_bulge else 10**res.x[0],
    }


def global_fit_1param(galaxies, alpha_val):
    """Global chi2 for constrained model: gamma=alpha/2, given alpha."""
    a0 = 1.2e-10  # will optimize separately
    total = 0
    for gal in galaxies:
        res = fit_galaxy_model(gal, nu_exp_constrained, a0, extra_params=(alpha_val,))
        if res: total += res['chi2']
        else: total += 1e6
    return total


# ============================================================
# MAIN
# ============================================================
def main():
    print("=" * 80)
    print("gs11b: BTFR VERIFICATION + CONSTRAINED MODEL TESTS")
    print("=" * 80)

    base_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(base_dir, 'Rotmod_LTG')
    master_file = os.path.join(base_dir, 'SPARC_Lelli2016c.mrt')

    # Load
    props = load_master_table(master_file)
    galaxies = []
    for fname in sorted(os.listdir(data_dir)):
        if fname.endswith('_rotmod.dat'):
            try:
                gal = Galaxy(os.path.join(data_dir, fname))
                if gal.npts >= 5:
                    galaxies.append(gal)
            except: pass
    print(f"Loaded {len(galaxies)} galaxies")

    # ============================================================
    # SECTION 1: BTFR DIRECTLY FROM SPARC
    # ============================================================
    print("\n" + "=" * 80)
    print("SECTION 1: BTFR FROM SPARC DATA (Vflat vs M_bar)")
    print("=" * 80)

    # M_bar = Yd * L[3.6] * M_sun/L_sun + 1.33 * M_HI (1.33 for He correction)
    # Use Yd = 0.5 (typical disk M/L at 3.6um) for stellar mass
    Yd_assumed = 0.5

    btfr_data = []
    for gal in galaxies:
        if gal.name not in props:
            continue
        p = props[gal.name]
        if p['Vflat'] <= 0 or p['Q'] > 2:  # need Vflat and good quality
            continue

        M_star = Yd_assumed * p['L36'] * 1e9  # Msun
        M_gas = 1.33 * p['MHI'] * 1e9         # Msun (He correction)
        M_bar = M_star + M_gas

        if M_bar <= 0:
            continue

        btfr_data.append({
            'name': gal.name,
            'M_bar': M_bar,
            'Vflat': p['Vflat'],
            'e_Vflat': p['e_Vflat'],
            'Q': p['Q'],
        })

    print(f"  Galaxies with Vflat > 0 and Q <= 2: {len(btfr_data)}")

    log_M = np.array([np.log10(b['M_bar']) for b in btfr_data])
    log_V = np.array([np.log10(b['Vflat']) for b in btfr_data])

    # Linear regression: log(M) = slope * log(V) + intercept
    slope, intercept, r_value, p_value, std_err = linregress(log_V, log_M)
    print(f"\n  BTFR fit: log(M_bar) = {slope:.3f} * log(Vflat) + {intercept:.3f}")
    print(f"  => M_bar ~ Vflat^{slope:.3f}")
    print(f"  => Vflat^{slope:.2f} ~ M_bar, i.e., BTFR exponent = {slope:.2f}")
    print(f"  R^2 = {r_value**2:.4f}")
    print(f"  Scatter: {np.std(log_M - slope*log_V - intercept):.3f} dex")

    # Also fit V = a * M^b
    slope2, intercept2, r2, p2, se2 = linregress(log_M, log_V)
    print(f"\n  Inverse: log(Vflat) = {slope2:.4f} * log(M_bar) + {intercept2:.4f}")
    print(f"  => Vflat ~ M_bar^{slope2:.4f}")
    print(f"  => M_bar ~ Vflat^{1/slope2:.2f}")

    # Compare with predictions
    print(f"\n  --- Comparison ---")
    print(f"  Observed BTFR: M ~ v^{slope:.2f} (this sample)")
    print(f"  McGaugh+ 2012: M ~ v^3.85 +/- 0.09")
    print(f"  MOND (gamma=0.5): M ~ v^4.00")
    print(f"  Winning eq (gamma=0.41): M ~ v^3.37 (asymptotic)")
    print(f"  Discrepancy winning vs observed: {abs(slope - 3.37)/slope * 100:.1f}%")
    print(f"  Discrepancy MOND vs observed: {abs(slope - 4.0)/slope * 100:.1f}%")

    # KEY: at what y do the SPARC galaxies actually sit?
    print(f"\n  --- y-range actually probed by SPARC ---")
    # For each galaxy, compute typical y = gbar/a0 at outer radii
    y_min_all = []
    y_max_all = []
    y_outer_all = []
    for gal in galaxies:
        gbar = gal.compute_gbar(0.5)
        mask = gbar > 0
        if np.sum(mask) < 2: continue
        y_vals = gbar[mask] / a0_obs
        y_min_all.append(np.min(y_vals))
        y_max_all.append(np.max(y_vals))
        # Outer 3 points
        if len(y_vals) >= 3:
            y_outer_all.append(np.median(y_vals[-3:]))
        else:
            y_outer_all.append(y_vals[-1])

    y_min_all = np.array(y_min_all)
    y_max_all = np.array(y_max_all)
    y_outer_all = np.array(y_outer_all)

    print(f"  y_min per galaxy: median={np.median(y_min_all):.4f}, min={np.min(y_min_all):.4e}")
    print(f"  y_max per galaxy: median={np.median(y_max_all):.2f}, max={np.max(y_max_all):.1f}")
    print(f"  y_outer (RC flat part): median={np.median(y_outer_all):.4f}")
    print(f"  Galaxies with y_outer < 0.1 (deep MOND): {np.sum(y_outer_all < 0.1)}")
    print(f"  Galaxies with y_outer < 0.01 (very deep MOND): {np.sum(y_outer_all < 0.01)}")
    print(f"  Galaxies with y_outer 0.1-1 (transition): {np.sum((y_outer_all >= 0.1) & (y_outer_all < 1.0))}")
    print(f"  Galaxies with y_outer > 1 (Newtonian-ish): {np.sum(y_outer_all >= 1.0)}")

    # ============================================================
    # SECTION 2: CONSTRAINED MODEL (gamma = alpha/2)
    # ============================================================
    print("\n" + "=" * 80)
    print("SECTION 2: CONSTRAINED MODEL gamma = alpha/2")
    print("=" * 80)

    # Global fit: optimize alpha and a0 jointly
    def global_chi2_constrained(params):
        log_a0 = params[0]
        alpha = params[1]
        a0 = 10**log_a0
        total = 0
        for gal in galaxies:
            res = fit_galaxy_model(gal, nu_exp_constrained, a0, extra_params=(alpha,))
            if res: total += res['chi2']
            else: total += 1e6
        return total

    print("  Optimizing (a0, alpha) with gamma=alpha/2...")
    result = differential_evolution(
        global_chi2_constrained,
        [(-10.3, -9.5), (0.3, 1.5)],
        seed=42, maxiter=80, tol=1e-4, polish=True
    )
    a0_con = 10**result.x[0]
    alpha_con = result.x[1]
    gamma_con = alpha_con / 2
    chi2_con = result.fun

    print(f"  Best: a0 = {a0_con:.4e} ({a0_con/1.2e-10:.4f} x 1.2e-10)")
    print(f"  Best: alpha = {alpha_con:.4f}")
    print(f"  Implied: gamma = {gamma_con:.4f}")

    # Compare with unconstrained (from gs10)
    chi2_uncon = 5623.8  # from gs10 Phase 2 Exponential
    a0_uncon = 1.1216e-10
    alpha_uncon = 0.8095
    gamma_uncon = 0.4059

    # Per-galaxy fits for constrained model
    per_gal_con = []
    for gal in galaxies:
        res = fit_galaxy_model(gal, nu_exp_constrained, a0_con, extra_params=(alpha_con,))
        if res:
            per_gal_con.append({'name': gal.name, **res})

    total_ndof_con = sum(r['ndof'] for r in per_gal_con) - 2  # 2 global params
    n_good_con = sum(1 for r in per_gal_con if r['chi2_red'] < 2.0)

    print(f"\n  Constrained (gamma=alpha/2): chi2={chi2_con:.1f}, chi2/ndof={chi2_con/total_ndof_con:.4f}")
    print(f"  Unconstrained:              chi2={chi2_uncon:.1f}")
    print(f"  Delta chi2 = {chi2_con - chi2_uncon:.1f} (constrained worse by this)")
    print(f"  Good fits (chi2_red < 2): {n_good_con}/{len(per_gal_con)}")

    n_total = sum(gal.npts for gal in galaxies)
    n_free_con = sum(2 if gal.has_bulge else 1 for gal in galaxies) + 2  # alpha, a0
    n_free_uncon = sum(2 if gal.has_bulge else 1 for gal in galaxies) + 3  # alpha, gamma, a0
    bic_con = chi2_con + n_free_con * np.log(n_total)
    bic_uncon = chi2_uncon + n_free_uncon * np.log(n_total)
    print(f"\n  BIC constrained:   {bic_con:.1f} (n_free={n_free_con})")
    print(f"  BIC unconstrained: {bic_uncon:.1f} (n_free={n_free_uncon})")
    print(f"  Delta BIC = {bic_con - bic_uncon:.1f}")
    if bic_con < bic_uncon:
        print(f"  => CONSTRAINED MODEL PREFERRED by BIC! (simpler + adequate)")
    else:
        print(f"  => Unconstrained model still preferred (delta={bic_con-bic_uncon:.1f})")

    # ============================================================
    # SECTION 3: FORCED gamma=0.5 (exact MOND deep limit)
    # ============================================================
    print("\n" + "=" * 80)
    print("SECTION 3: FORCED gamma = 0.5 (exact MOND deep limit)")
    print("=" * 80)

    def global_chi2_mond_deep(params):
        log_a0 = params[0]
        alpha = params[1]
        a0 = 10**log_a0
        total = 0
        for gal in galaxies:
            res = fit_galaxy_model(gal, nu_exp, a0, extra_params=(alpha, 0.5))
            if res: total += res['chi2']
            else: total += 1e6
        return total

    print("  Optimizing (a0, alpha) with gamma=0.5 fixed...")
    result05 = differential_evolution(
        global_chi2_mond_deep,
        [(-10.3, -9.5), (0.3, 1.5)],
        seed=42, maxiter=80, tol=1e-4, polish=True
    )
    a0_05 = 10**result05.x[0]
    alpha_05 = result05.x[1]
    chi2_05 = result05.fun

    print(f"  Best: a0 = {a0_05:.4e} ({a0_05/1.2e-10:.4f} x 1.2e-10)")
    print(f"  Best: alpha = {alpha_05:.4f}")
    print(f"  chi2 = {chi2_05:.1f}")

    per_gal_05 = []
    for gal in galaxies:
        res = fit_galaxy_model(gal, nu_exp, a0_05, extra_params=(alpha_05, 0.5))
        if res:
            per_gal_05.append({'name': gal.name, **res})

    total_ndof_05 = sum(r['ndof'] for r in per_gal_05) - 2
    n_good_05 = sum(1 for r in per_gal_05 if r['chi2_red'] < 2.0)

    n_free_05 = sum(2 if gal.has_bulge else 1 for gal in galaxies) + 2
    bic_05 = chi2_05 + n_free_05 * np.log(n_total)

    print(f"\n  Comparison:")
    print(f"  {'Model':<30} {'chi2':<10} {'chi2/ndof':<12} {'BIC':<10} {'Good':<8}")
    print(f"  {'-'*70}")
    print(f"  {'Exp (free alpha,gamma)':<30} {chi2_uncon:<10.1f} {chi2_uncon/3169:<12.4f} {bic_uncon:<10.1f} {'128/171':<8}")
    print(f"  {'Exp (gamma=alpha/2)':<30} {chi2_con:<10.1f} {chi2_con/total_ndof_con:<12.4f} {bic_con:<10.1f} {f'{n_good_con}/171':<8}")
    print(f"  {'Exp (gamma=0.5, MOND)':<30} {chi2_05:<10.1f} {chi2_05/total_ndof_05:<12.4f} {bic_05:<10.1f} {f'{n_good_05}/171':<8}")

    # Per-galaxy comparison: which galaxies suffer most from gamma=0.5?
    print(f"\n--- Per-galaxy: gamma=0.41 vs gamma=0.5 ---")

    con_dict = {r['name']: r['chi2'] for r in per_gal_con}
    g05_dict = {r['name']: r['chi2'] for r in per_gal_05}

    deltas = []
    for name in set(con_dict.keys()) & set(g05_dict.keys()):
        delta = g05_dict[name] - con_dict[name]  # positive = gamma=0.5 worse
        deltas.append((name, delta, con_dict[name], g05_dict[name]))

    deltas.sort(key=lambda x: x[1], reverse=True)

    print(f"\n  Galaxies where gamma=0.5 is MUCH WORSE:")
    print(f"  {'Galaxy':<15} {'dchi2':<10} {'chi2(g=a/2)':<14} {'chi2(g=0.5)':<14}")
    for name, d, c1, c2 in deltas[:15]:
        print(f"  {name:<15} {d:<10.2f} {c1:<14.2f} {c2:<14.2f}")

    print(f"\n  Galaxies where gamma=0.5 is BETTER:")
    for name, d, c1, c2 in deltas[-10:]:
        print(f"  {name:<15} {d:<10.2f} {c1:<14.2f} {c2:<14.2f}")

    prefer_con = sum(1 for _, d, _, _ in deltas if d > 2)
    prefer_05 = sum(1 for _, d, _, _ in deltas if d < -2)
    print(f"\n  Prefer gamma=alpha/2: {prefer_con} galaxies")
    print(f"  Prefer gamma=0.5:     {prefer_05} galaxies")
    print(f"  No strong preference: {len(deltas) - prefer_con - prefer_05}")

    # ============================================================
    # SECTION 4: EFFECTIVE BTFR FROM MODEL (at observed radii)
    # ============================================================
    print("\n" + "=" * 80)
    print("SECTION 4: MODEL-PREDICTED BTFR AT OBSERVED RADII")
    print("=" * 80)

    # For each galaxy with Vflat, compute what Vflat the model predicts
    # using the best-fit Yd and the last few data points

    for model_label, a0_m, alpha_m, gamma_m in [
        ("Constrained (g=a/2)", a0_con, alpha_con, alpha_con/2),
        ("Forced g=0.5", a0_05, alpha_05, 0.5),
        ("McGaugh", 0.765e-10, None, None),
    ]:
        print(f"\n--- {model_label} ---")

        predicted_vflat = []
        observed_vflat = []

        for gal in galaxies:
            if gal.name not in props or props[gal.name]['Vflat'] <= 0:
                continue
            if props[gal.name]['Q'] > 2:
                continue

            # Fit Yd for this model
            if gamma_m is not None:
                res = fit_galaxy_model(gal, nu_exp, a0_m, extra_params=(alpha_m, gamma_m))
            else:
                res = fit_galaxy_model(gal, nu_mcgaugh, a0_m)

            if res is None:
                continue

            Yd = res['Yd']
            Yb = res['Yb']

            # Compute v at outermost data point
            gbar = gal.compute_gbar(Yd, Yb)
            mask = gbar > 0
            if np.sum(mask) < 3:
                continue

            y_vals = gbar[mask] / (a0_m if a0_m else 0.765e-10)
            y_vals = np.maximum(y_vals, 1e-10)

            if gamma_m is not None:
                nu_vals = nu_exp(y_vals, alpha_m, gamma_m)
            else:
                nu_vals = nu_mcgaugh(y_vals)

            gobs_pred = gbar[mask] * nu_vals
            v_pred = np.sqrt(np.abs(gobs_pred) * gal.rad_m[mask]) / 1e3

            # Use last 3 points as "flat" velocity
            v_flat_pred = np.median(v_pred[-3:])
            v_flat_obs = props[gal.name]['Vflat']

            M_star = 0.5 * props[gal.name]['L36'] * 1e9
            M_gas = 1.33 * props[gal.name]['MHI'] * 1e9
            M_bar = M_star + M_gas
            if M_bar <= 0:
                continue

            predicted_vflat.append(v_flat_pred)
            observed_vflat.append(v_flat_obs)

        if len(predicted_vflat) > 10:
            pred = np.array(predicted_vflat)
            obs = np.array(observed_vflat)

            # How well does model reproduce Vflat?
            ratio = pred / obs
            print(f"  N galaxies: {len(pred)}")
            print(f"  V_pred/V_obs: median={np.median(ratio):.3f}, "
                  f"mean={np.mean(ratio):.3f}, std={np.std(ratio):.3f}")

            # Residual scatter
            log_res = np.log10(pred) - np.log10(obs)
            print(f"  log(Vpred/Vobs): scatter = {np.std(log_res):.3f} dex")

    # ============================================================
    # SECTION 5: BTFR SLOPE FROM DIFFERENT MODELS
    # ============================================================
    print("\n" + "=" * 80)
    print("SECTION 5: BTFR SLOPE PREDICTED BY EACH MODEL")
    print("=" * 80)

    # Compute M_bar for each galaxy, then fit BTFR using model-predicted Vflat
    for model_label, a0_m, alpha_m, gamma_m in [
        ("Constrained (g=a/2)", a0_con, alpha_con, alpha_con/2),
        ("Forced g=0.5", a0_05, alpha_05, 0.5),
        ("McGaugh", 0.765e-10, None, None),
        ("MOND simple", 0.783e-10, None, None),
    ]:
        log_M_list = []
        log_V_list = []

        for gal in galaxies:
            if gal.name not in props or props[gal.name]['Vflat'] <= 0:
                continue
            if props[gal.name]['Q'] > 2:
                continue

            if gamma_m is not None:
                res = fit_galaxy_model(gal, nu_exp, a0_m, extra_params=(alpha_m, gamma_m))
            elif 'McGaugh' in model_label:
                res = fit_galaxy_model(gal, nu_mcgaugh, a0_m)
            else:
                res = fit_galaxy_model(gal, nu_mond_simple, a0_m)

            if res is None:
                continue

            Yd = res['Yd']
            # Total baryonic mass with fitted Yd
            M_star = Yd * props[gal.name]['L36'] * 1e9
            M_gas = 1.33 * props[gal.name]['MHI'] * 1e9
            M_bar = M_star + M_gas
            if M_bar <= 0:
                continue

            log_M_list.append(np.log10(M_bar))
            log_V_list.append(np.log10(props[gal.name]['Vflat']))

        if len(log_M_list) > 10:
            log_M_arr = np.array(log_M_list)
            log_V_arr = np.array(log_V_list)
            s, i, r, p, se = linregress(log_V_arr, log_M_arr)
            print(f"  {model_label:<30}: M ~ v^{s:.2f} (r2={r**2:.4f}, scatter={np.std(log_M_arr-s*log_V_arr-i):.3f} dex)")

    # Raw observed BTFR for comparison
    s_raw, i_raw, r_raw, p_raw, se_raw = linregress(log_V, log_M)
    print(f"  {'Raw observed (Yd=0.5 fixed)':<30}: M ~ v^{s_raw:.2f} (r2={r_raw**2:.4f})")

    print("\n" + "=" * 80)
    print("DONE")
    print("=" * 80)


if __name__ == '__main__':
    main()
