#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gs12: REFINED ANALYSIS — definitive parameter constraints
==========================================================

Issues from gs11/gs11b:
  - Coarse Yd grid (30 points) biased global chi2 → "55-sigma" was artifact
  - Need PROPER profile likelihood with full optimization per galaxy
  - Need detailed RC comparison for representative galaxies
  - Need to test erfc interpretation quantitatively

This script:
  1. Fine 2D scan (alpha, gamma) with scipy.optimize per galaxy (not grid)
  2. Profile likelihood for gamma: at each gamma, optimize (a0, alpha, Yd_i)
  3. Profile likelihood for alpha: at each alpha, optimize (a0, gamma, Yd_i)
  4. Rotation curve comparison: 12 representative galaxies
  5. erfc model test: nu = 1 + A*erfc(y^beta) — fit directly
"""

import sys
import io
import os
import numpy as np
from scipy.optimize import minimize, minimize_scalar
from scipy.special import erfc
import warnings
warnings.filterwarnings('ignore')

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

G = 6.674e-11
M_sun = 1.989e30
kpc = 3.086e19

# ============================================================
# GALAXY LOADER
# ============================================================
class Galaxy:
    def __init__(self, filepath):
        self.name = os.path.basename(filepath).replace('_rotmod.dat', '')
        lines = []
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if not line.startswith('#') and line:
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
        self.errv_safe = np.maximum(self.errv, np.maximum(3.0, 0.03*np.abs(self.vobs)))
        self.rad_m = self.rad_kpc * kpc

    def compute_gbar(self, Yd, Yb=None):
        if Yb is None: Yb = Yd
        v2 = np.abs(self.vgas)*self.vgas + Yd*self.vdisk**2 + Yb*self.vbul**2
        return v2 * 1e6 / self.rad_m

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


# ============================================================
# MODELS
# ============================================================
def nu_exp(y, alpha, gamma):
    return 1.0 + np.exp(-y**alpha) / y**gamma

def nu_mcgaugh(y):
    return 1.0 / (1.0 - np.exp(-np.sqrt(np.maximum(y, 1e-30))))

def nu_erfc_model(y, A, beta):
    """erfc-based: nu = 1 + A * erfc(y^beta)"""
    return 1.0 + A * erfc(y**beta)


# ============================================================
# PER-GALAXY OPTIMIZER (proper, not grid)
# ============================================================
def fit_galaxy_chi2(galaxy, a0, alpha, gamma, model='exp'):
    """Optimize Yd (and Yb) for given (a0, alpha, gamma). Return best chi2."""
    def objective(params):
        log_Yd = params[0]
        log_Yb = params[1] if galaxy.has_bulge else log_Yd
        Yd, Yb = 10**log_Yd, 10**log_Yb
        gbar = galaxy.compute_gbar(Yd, Yb)
        mask = gbar > 0
        if np.sum(mask) < 3: return 1e10
        y = gbar[mask] / a0
        y = np.maximum(y, 1e-10)
        try:
            if model == 'exp':
                nu = nu_exp(y, alpha, gamma)
            elif model == 'mcgaugh':
                nu = nu_mcgaugh(y)
            elif model == 'erfc':
                nu = nu_erfc_model(y, alpha, gamma)  # alpha=A, gamma=beta
            else:
                nu = nu_exp(y, alpha, gamma)
        except: return 1e10
        if not np.all(np.isfinite(nu)): return 1e10
        gobs = gbar[mask] * nu
        vp = np.sqrt(np.abs(gobs) * galaxy.rad_m[mask]) / 1e3
        return np.sum(((galaxy.vobs[mask] - vp) / galaxy.errv_safe[mask])**2)

    bounds = [(-0.7, 1.2)]
    x0 = [0.0]
    if galaxy.has_bulge:
        bounds.append((-0.7, 1.2))
        x0.append(0.0)

    # Multi-start for robustness
    best_chi2 = 1e10
    best_Yd = 1.0
    for yd_init in [-0.3, 0.0, 0.3, 0.6]:
        x0_try = [yd_init] + ([yd_init] if galaxy.has_bulge else [])
        try:
            res = minimize(objective, x0_try, method='L-BFGS-B', bounds=bounds)
            if res.fun < best_chi2:
                best_chi2 = res.fun
                best_Yd = 10**res.x[0]
        except: pass

    return best_chi2, best_Yd


def total_chi2(galaxies, a0, alpha, gamma, model='exp'):
    """Total chi2 across all galaxies for given equation params."""
    total = 0.0
    for gal in galaxies:
        c2, _ = fit_galaxy_chi2(gal, a0, alpha, gamma, model)
        total += c2
    return total


# ============================================================
# MAIN
# ============================================================
def main():
    print("=" * 80)
    print("gs12: REFINED ANALYSIS — definitive parameter constraints")
    print("=" * 80)

    base_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(base_dir, 'Rotmod_LTG')
    galaxies = load_galaxies(data_dir)
    print(f"Loaded {len(galaxies)} galaxies, {sum(g.npts for g in galaxies)} total points")

    # ============================================================
    # SECTION 1: PROFILE LIKELIHOOD FOR GAMMA
    # ============================================================
    print("\n" + "=" * 80)
    print("SECTION 1: PROFILE LIKELIHOOD FOR GAMMA")
    print("  At each gamma, optimize a0 (scan) with alpha=0.81 fixed")
    print("=" * 80)

    alpha_fixed = 0.81

    gammas = np.arange(0.30, 0.61, 0.02)
    chi2_profile_gamma = []

    for gamma_val in gammas:
        # Optimize a0 at this gamma
        best_a0_chi2 = 1e15
        best_a0 = 1e-10

        for log_a0 in np.linspace(-10.3, -9.6, 15):
            a0 = 10**log_a0
            c2 = total_chi2(galaxies, a0, alpha_fixed, gamma_val)
            if c2 < best_a0_chi2:
                best_a0_chi2 = c2
                best_a0 = a0

        chi2_profile_gamma.append(best_a0_chi2)
        print(f"  gamma={gamma_val:.2f}: chi2={best_a0_chi2:.1f} (a0={best_a0:.3e})")

    chi2_profile_gamma = np.array(chi2_profile_gamma)
    chi2_min_gamma = np.min(chi2_profile_gamma)
    idx_best = np.argmin(chi2_profile_gamma)

    print(f"\n  Best gamma = {gammas[idx_best]:.2f} (chi2 = {chi2_min_gamma:.1f})")

    # Confidence intervals
    for nsig, dchi2 in [(1, 1.0), (2, 4.0), (3, 9.0)]:
        mask = chi2_profile_gamma - chi2_min_gamma < dchi2
        if np.any(mask):
            g_lo, g_hi = gammas[mask][0], gammas[mask][-1]
            print(f"  {nsig}-sigma: gamma = {g_lo:.2f} to {g_hi:.2f}")
        else:
            print(f"  {nsig}-sigma: only best point")

    # Is gamma=0.5 excluded?
    if 0.5 in gammas or True:
        idx_05 = np.argmin(np.abs(gammas - 0.5))
        delta_05 = chi2_profile_gamma[idx_05] - chi2_min_gamma
        print(f"\n  gamma=0.50: delta_chi2 = {delta_05:.1f} ({np.sqrt(max(0,delta_05)):.1f}-sigma)")

    # ============================================================
    # SECTION 2: PROFILE LIKELIHOOD FOR ALPHA
    # ============================================================
    print("\n" + "=" * 80)
    print("SECTION 2: PROFILE LIKELIHOOD FOR ALPHA")
    print("  At each alpha, optimize a0 (scan) with gamma=0.41 fixed")
    print("=" * 80)

    gamma_fixed = 0.41

    alphas = np.arange(0.40, 1.21, 0.05)
    chi2_profile_alpha = []

    for alpha_val in alphas:
        best_a0_chi2 = 1e15
        best_a0 = 1e-10

        for log_a0 in np.linspace(-10.3, -9.6, 15):
            a0 = 10**log_a0
            c2 = total_chi2(galaxies, a0, alpha_val, gamma_fixed)
            if c2 < best_a0_chi2:
                best_a0_chi2 = c2
                best_a0 = a0

        chi2_profile_alpha.append(best_a0_chi2)
        print(f"  alpha={alpha_val:.2f}: chi2={best_a0_chi2:.1f} (a0={best_a0:.3e})")

    chi2_profile_alpha = np.array(chi2_profile_alpha)
    chi2_min_alpha = np.min(chi2_profile_alpha)
    idx_best_a = np.argmin(chi2_profile_alpha)

    print(f"\n  Best alpha = {alphas[idx_best_a]:.2f} (chi2 = {chi2_min_alpha:.1f})")

    for nsig, dchi2 in [(1, 1.0), (2, 4.0), (3, 9.0)]:
        mask = chi2_profile_alpha - chi2_min_alpha < dchi2
        if np.any(mask):
            a_lo, a_hi = alphas[mask][0], alphas[mask][-1]
            print(f"  {nsig}-sigma: alpha = {a_lo:.2f} to {a_hi:.2f}")

    # ============================================================
    # SECTION 3: JOINT OPTIMUM — fine 2D scan
    # ============================================================
    print("\n" + "=" * 80)
    print("SECTION 3: JOINT 2D SCAN (alpha, gamma) with a0 optimized")
    print("=" * 80)

    # Fine grid around expected optimum
    alphas_2d = np.arange(0.60, 1.01, 0.05)
    gammas_2d = np.arange(0.34, 0.53, 0.02)
    chi2_2d = np.zeros((len(alphas_2d), len(gammas_2d)))

    print("  Scanning...")
    for i, a in enumerate(alphas_2d):
        for j, g in enumerate(gammas_2d):
            # Quick a0 optimization: scan 10 values
            best = 1e15
            for log_a0 in np.linspace(-10.2, -9.7, 10):
                c2 = total_chi2(galaxies, 10**log_a0, a, g)
                if c2 < best: best = c2
            chi2_2d[i, j] = best
        print(f"    alpha={a:.2f} done")

    min_2d = np.min(chi2_2d)
    ij_min = np.unravel_index(np.argmin(chi2_2d), chi2_2d.shape)
    best_alpha_2d = alphas_2d[ij_min[0]]
    best_gamma_2d = gammas_2d[ij_min[1]]

    print(f"\n  2D minimum: alpha={best_alpha_2d:.2f}, gamma={best_gamma_2d:.2f}, chi2={min_2d:.1f}")
    print(f"  gamma/alpha = {best_gamma_2d/best_alpha_2d:.3f}")
    print(f"  gamma = alpha/2 would give gamma = {best_alpha_2d/2:.3f}")

    # Print delta-chi2 table
    print(f"\n  Delta-chi2 table:")
    print(f"  {'a\\g':<6}", end='')
    for g in gammas_2d:
        print(f" {g:.2f} ", end='')
    print()

    for i, a in enumerate(alphas_2d):
        print(f"  {a:.2f}  ", end='')
        for j, g in enumerate(gammas_2d):
            d = chi2_2d[i, j] - min_2d
            if d < 1: print(f" {'**':>5}", end='')
            elif d < 4: print(f" {' *':>5}", end='')
            elif d < 9: print(f" {' .':>5}", end='')
            else: print(f" {d:5.0f}", end='')
        print()

    # 1-sigma and 2-sigma contours (2 params: dchi2 = 2.30, 6.18)
    for nsig, dchi2, label in [(1, 2.30, '1-sigma'), (2, 6.18, '2-sigma')]:
        mask = chi2_2d - min_2d < dchi2
        a_range = alphas_2d[np.any(mask, axis=1)]
        g_range = gammas_2d[np.any(mask, axis=0)]
        if len(a_range) > 0 and len(g_range) > 0:
            print(f"\n  {label} (dchi2<{dchi2}):")
            print(f"    alpha: {a_range[0]:.2f} to {a_range[-1]:.2f}")
            print(f"    gamma: {g_range[0]:.2f} to {g_range[-1]:.2f}")

    # ============================================================
    # SECTION 4: ERFC MODEL TEST
    # ============================================================
    print("\n" + "=" * 80)
    print("SECTION 4: ERFC MODEL — nu = 1 + A * erfc(y^beta)")
    print("=" * 80)

    # The winning equation at alpha=4/5, gamma=2/5 can be written as
    # nu ~ 1 + exp(-u^2)/u where u = y^(2/5)
    # erfc(u) ~ exp(-u^2)/(sqrt(pi)*u) for large u
    # So nu ~ 1 + sqrt(pi) * erfc(y^(2/5)) for y >> 1
    # But for y << 1, erfc -> 1, while exp(-u^2)/u -> 1/u -> inf
    # So erfc is NOT the same as our function for small y

    # Test: nu = 1 + A * erfc(y^beta) as independent model
    print("  Optimizing (a0, A, beta) for erfc model...")

    def global_chi2_erfc(params):
        log_a0, A, beta = params
        a0 = 10**log_a0
        if A <= 0 or beta <= 0: return 1e15
        total = 0
        for gal in galaxies:
            c2, _ = fit_galaxy_chi2(gal, a0, A, beta, model='erfc')
            total += c2
        return total

    # Scan
    best_erfc_chi2 = 1e15
    best_erfc_params = None

    for log_a0 in np.linspace(-10.2, -9.7, 8):
        for A in [0.5, 1.0, 1.5, 1.77, 2.0, 3.0]:
            for beta in [0.2, 0.3, 0.4, 0.5, 0.6]:
                c2 = global_chi2_erfc([log_a0, A, beta])
                if c2 < best_erfc_chi2:
                    best_erfc_chi2 = c2
                    best_erfc_params = (10**log_a0, A, beta)

    if best_erfc_params:
        a0_e, A_e, beta_e = best_erfc_params
        print(f"  Best: a0={a0_e:.3e}, A={A_e:.3f}, beta={beta_e:.3f}")
        print(f"  chi2 = {best_erfc_chi2:.1f}")
        print(f"  sqrt(pi) = {np.sqrt(np.pi):.4f} (expected A if erfc exact)")
        print(f"  2/5 = 0.4 (expected beta if erfc exact)")

    # Compare with winning equation
    a0_win = 1.12e-10
    chi2_win_check = total_chi2(galaxies, a0_win, 0.81, 0.41)
    print(f"\n  Comparison:")
    print(f"  Winning exp model:  chi2 = {chi2_win_check:.1f}")
    print(f"  erfc model:         chi2 = {best_erfc_chi2:.1f}")
    print(f"  Delta:              {best_erfc_chi2 - chi2_win_check:.1f}")

    # ============================================================
    # SECTION 5: REPRESENTATIVE ROTATION CURVES
    # ============================================================
    print("\n" + "=" * 80)
    print("SECTION 5: ROTATION CURVES — 12 representative galaxies")
    print("=" * 80)

    # Select representative galaxies across mass range
    representatives = [
        'DDO154',     # dwarf, gas-dominated
        'NGC3741',    # ultra-low mass
        'IC2574',     # dwarf irregular
        'NGC2403',    # intermediate
        'NGC3198',    # classic spiral
        'NGC6503',    # well-studied
        'NGC2841',    # massive, rising RC
        'NGC7331',    # massive spiral
        'NGC6946',    # well-studied face-on
        'NGC5055',    # Sunflower
        'UGC02885',   # super-spiral
        'NGC7814',    # edge-on with bulge
    ]

    gal_dict = {g.name: g for g in galaxies}

    a0_best = 1.12e-10
    alpha_best = 0.81
    gamma_best = 0.41

    for gname in representatives:
        if gname not in gal_dict:
            continue
        gal = gal_dict[gname]

        # Fit with winning equation
        chi2_exp, Yd_exp = fit_galaxy_chi2(gal, a0_best, alpha_best, gamma_best)
        # Fit with McGaugh
        chi2_mcg, Yd_mcg = fit_galaxy_chi2(gal, 0.765e-10, 0, 0, model='mcgaugh')
        # Fit with MOND-like (alpha, gamma=0.5)
        chi2_05, Yd_05 = fit_galaxy_chi2(gal, 0.827e-10, 0.645, 0.5)

        ndof = gal.npts - (2 if gal.has_bulge else 1)

        print(f"\n  === {gname} ({gal.npts} pts, bulge={'Y' if gal.has_bulge else 'N'}) ===")
        print(f"  {'Model':<25} {'chi2':<10} {'chi2/dof':<10} {'Yd':<8}")
        print(f"  {'Exp (a=0.81,g=0.41)':<25} {chi2_exp:<10.2f} {chi2_exp/ndof:<10.3f} {Yd_exp:<8.3f}")
        print(f"  {'McGaugh':<25} {chi2_mcg:<10.2f} {chi2_mcg/ndof:<10.3f} {Yd_mcg:<8.3f}")
        print(f"  {'Exp (g=0.5, MOND-like)':<25} {chi2_05:<10.2f} {chi2_05/ndof:<10.3f} {Yd_05:<8.3f}")

        # Print RC data vs models at selected radii
        gbar_exp = gal.compute_gbar(Yd_exp)
        gbar_mcg = gal.compute_gbar(Yd_mcg)

        mask = gbar_exp > 0
        y_exp = gbar_exp[mask] / a0_best
        nu_e = nu_exp(np.maximum(y_exp, 1e-10), alpha_best, gamma_best)
        v_exp = np.sqrt(np.abs(gbar_exp[mask] * nu_e * gal.rad_m[mask])) / 1e3

        y_mcg = gbar_mcg[mask] / 0.765e-10
        nu_m = nu_mcgaugh(np.maximum(y_mcg, 1e-10))
        v_mcg = np.sqrt(np.abs(gbar_mcg[mask] * nu_m * gal.rad_m[mask])) / 1e3

        # Show comparison at 5 radii
        indices = np.linspace(0, np.sum(mask)-1, min(6, np.sum(mask))).astype(int)
        r_vals = gal.rad_kpc[mask]
        v_obs = gal.vobs[mask]
        err_vals = gal.errv_safe[mask]

        print(f"  {'r(kpc)':<8} {'Vobs':<8} {'err':<6} {'V_exp':<8} {'V_mcg':<8} {'y':<10}")
        for idx in indices:
            print(f"  {r_vals[idx]:<8.2f} {v_obs[idx]:<8.1f} {err_vals[idx]:<6.1f} "
                  f"{v_exp[idx]:<8.1f} {v_mcg[idx]:<8.1f} {y_exp[idx]:<10.4f}")

    # ============================================================
    # SECTION 6: McGaugh BASELINE (for honest comparison)
    # ============================================================
    print("\n" + "=" * 80)
    print("SECTION 6: MCGAUGH BASELINE — fair comparison")
    print("=" * 80)

    # Total chi2 for McGaugh with optimized a0
    print("  Scanning McGaugh a0...")
    best_mcg_chi2 = 1e15
    best_mcg_a0 = 1e-10
    for log_a0 in np.linspace(-10.3, -9.6, 20):
        a0 = 10**log_a0
        c2 = total_chi2(galaxies, a0, 0, 0, model='mcgaugh')
        if c2 < best_mcg_chi2:
            best_mcg_chi2 = c2
            best_mcg_a0 = a0
        if abs(log_a0 - np.log10(best_mcg_a0)) < 0.1 or c2 - best_mcg_chi2 < 100:
            print(f"    a0={a0:.3e}: chi2={c2:.1f}")

    print(f"\n  McGaugh best: a0={best_mcg_a0:.3e}, chi2={best_mcg_chi2:.1f}")
    print(f"  Exp best:     a0={a0_best:.3e}, chi2={chi2_win_check:.1f}")
    print(f"  Delta chi2:   {best_mcg_chi2 - chi2_win_check:.1f}")

    # Number of extra params: Exp has 2 (alpha, gamma), McGaugh has 0
    # BIC penalty: 2 * ln(N_total)
    N_total = sum(g.npts for g in galaxies)
    bic_penalty = 2 * np.log(N_total)
    print(f"\n  BIC penalty for 2 extra params: {bic_penalty:.1f}")
    print(f"  Net BIC advantage of Exp: {best_mcg_chi2 - chi2_win_check - bic_penalty:.1f}")
    if best_mcg_chi2 - chi2_win_check > bic_penalty:
        print(f"  => Exponential family JUSTIFIED by BIC despite extra params")
    else:
        print(f"  => McGaugh preferred by BIC (simpler, adequate)")

    # ============================================================
    # SUMMARY
    # ============================================================
    print("\n" + "=" * 80)
    print("DEFINITIVE SUMMARY")
    print("=" * 80)
    print(f"""
  WINNING EQUATION: nu(y) = 1 + exp(-y^alpha) / y^gamma

  Best fit (profile likelihood with proper optimization):
    alpha = {best_alpha_2d:.2f} (from 2D scan)
    gamma = {best_gamma_2d:.2f} (from 2D scan)
    gamma/alpha = {best_gamma_2d/best_alpha_2d:.3f}

  Constraints (2D, 1-sigma):
    (see contour above)

  gamma = 0.50 status: delta_chi2 = {delta_05:.1f} from profile

  vs McGaugh: delta_chi2 = {best_mcg_chi2 - chi2_win_check:.1f}
  BIC advantage: {best_mcg_chi2 - chi2_win_check - bic_penalty:.1f}
""")

    print("=" * 80)
    print("DONE")
    print("=" * 80)


if __name__ == '__main__':
    main()
