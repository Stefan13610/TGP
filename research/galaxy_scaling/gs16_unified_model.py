#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gs16: UNIFIED MODEL — GALAXY + CLUSTER
========================================
gs15 showed missing baryons can't close the gap. Key tension:
  - Galaxies: gamma=0.4 wins (SPARC, 171 galaxies, 30sigma vs 0.5)
  - Clusters: gamma=0.4 underpredicts by x1.84 (WORSE than MOND)

Critical observation: galaxies & clusters OVERLAP in y ~ 0.05!
At same y, galaxies want gamma=0.4 but clusters want gamma>0.5.
=> gamma(y) alone CAN'T explain both.

What DIFFERS between galaxies and clusters?
  1. GEOMETRY: galaxies = disks, clusters = spheres
  2. TOTAL MASS: 10^9-10^11 vs 10^13-10^15 Msun
  3. SIZE: 1-50 kpc vs 500-2000 kpc
  4. TRANSITION SCALE: H = sqrt(r_S * r_H) = sqrt(GM/(c*H0))
     - Galaxy (10^11): H ~ 10 kpc, r_out/H ~ 1-5
     - Cluster (10^14): H ~ 300 kpc, r500/H ~ 3-5
     - Cluster (10^15): H ~ 1 Mpc, r500/H ~ 1.5

This script tests 5 models:

  A) gamma(M) — gamma depends on total system mass
  B) a0(R) — acceleration scale depends on system size
  C) Geometry factor — spherical vs disk correction
  D) Second transition — additional boost at very low y
  E) Running gamma(y) with y-dependent floor

All models fitted to BOTH galaxy + cluster data simultaneously.
"""

import sys
import io
import os
import numpy as np
from scipy.optimize import minimize, minimize_scalar
import warnings
warnings.filterwarnings('ignore')

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

G_SI = 6.674e-11
M_sun = 1.989e30
kpc_m = 3.086e19
Mpc_m = 3.086e22
c_light = 3e8
H0 = 70e3 / Mpc_m  # 70 km/s/Mpc in 1/s

a0_ref = 1.12e-10

# ============================================================
# GALAXY LOADER
# ============================================================
class Galaxy:
    def __init__(self, filepath):
        self.name = os.path.basename(filepath).replace('_rotmod.dat', '')
        lines = []
        with open(filepath, 'r') as f:
            for line in f:
                ls = line.strip()
                if not ls.startswith('#') and ls:
                    lines.append(ls.split())
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
        self.rad_m = self.rad_kpc * kpc_m
        self.is_cluster = False
        self.geometry = 'disk'

    def compute_gbar(self, Yd, Yb=None):
        if Yb is None: Yb = Yd
        v2 = np.abs(self.vgas)*self.vgas + Yd*self.vdisk**2 + Yb*self.vbul**2
        return v2 * 1e6 / self.rad_m

    def estimate_mass(self, Yd=0.5):
        v2 = np.abs(self.vgas[-1])*self.vgas[-1] + Yd*self.vdisk[-1]**2
        return v2 * 1e6 * self.rad_m[-1] / G_SI / M_sun


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
# CLUSTER DATA STRUCTURE (matching galaxy interface)
# ============================================================
class ClusterPoint:
    """Single cluster at r500 — treated as 1-point 'rotation curve'."""
    def __init__(self, name, r500_Mpc, M500_1e14, Mbar_frac):
        self.name = name
        self.is_cluster = True
        self.geometry = 'sphere'
        self.npts = 1
        self.has_bulge = False

        self.r500_m = r500_Mpc * Mpc_m
        self.M500 = M500_1e14 * 1e14 * M_sun  # kg
        self.Mbar = Mbar_frac * self.M500
        self.M_total_Msun = M500_1e14 * 1e14

        # "Observed" velocity = sqrt(G*M500/r500)
        self.v_obs = np.sqrt(G_SI * self.M500 / self.r500_m) / 1e3  # km/s
        self.v_err = 0.05 * self.v_obs  # 5% error
        self.g_bar = G_SI * self.Mbar / self.r500_m**2
        self.g_obs = G_SI * self.M500 / self.r500_m**2

    def get_chi2(self, g_pred):
        """Chi2 from velocity comparison."""
        v_pred = np.sqrt(abs(g_pred) * self.r500_m) / 1e3
        return ((v_pred - self.v_obs) / self.v_err)**2


def build_clusters():
    """Build cluster sample from gs13."""
    data = [
        # (name, r500_Mpc, M500_1e14, fbar)
        ("A0122",  0.89, 2.26, 0.112), ("A1651",  1.18, 5.15, 0.143),
        ("A2401",  0.68, 0.95, 0.118), ("A2721",  1.03, 3.46, 0.142),
        ("A2811",  1.04, 3.59, 0.138), ("A2955",  0.68, 0.99, 0.097),
        ("A2984",  0.67, 0.95, 0.152), ("A3112",  1.02, 3.23, 0.154),
        ("A3693",  0.90, 2.26, 0.133), ("A4010",  0.92, 2.41, 0.142),
        ("AS0084", 0.91, 2.37, 0.110), ("AS0296", 0.78, 1.45, 0.095),
        ("Coma",   1.31, 7.2,  0.164), ("Perseus",1.27, 6.6,  0.180),
        ("Virgo",  0.77, 1.4,  0.120), ("Centaurus",0.82,1.8, 0.125),
        ("A2029",  1.42, 8.71, 0.158), ("A0478",  1.28, 6.58, 0.195),
        ("A2390",  1.50, 11.8, 0.161), ("A1795",  1.09, 4.54, 0.138),
        ("A2142",  1.37, 8.17, 0.194), ("A2319",  1.55, 12.0, 0.168),
        ("A644",   1.05, 3.70, 0.150), ("A3266",  1.25, 6.20, 0.158),
        ("Fornax", 0.45, 0.28, 0.095), ("MKW4",   0.50, 0.38, 0.095),
        ("NGC5044g",0.43,0.24, 0.095), ("A3571",  1.13, 5.02, 0.123),
        ("Bullet", 1.35, 7.8,  0.155),
    ]
    return [ClusterPoint(n, r, m, f) for n, r, m, f in data]


# ============================================================
# INTERPOLATION FUNCTIONS
# ============================================================
def nu_exp(y, alpha=0.8, gamma=0.4):
    y = np.maximum(y, 1e-15)
    return 1.0 + np.exp(-np.power(y, alpha)) / np.power(y, gamma)


# ============================================================
# MODEL A: gamma(M) = gamma_0 + delta * log10(M/M_ref)
# ============================================================
def model_A_chi2(params, galaxies, clusters, alpha=0.8):
    gamma_0, delta, log_a0 = params
    a0 = 10**log_a0
    if gamma_0 < 0.1 or gamma_0 > 0.8: return 1e15
    if abs(delta) > 0.3: return 1e15
    M_ref = 1e10  # Msun

    total_chi2 = 0

    # Galaxy part — subsample for speed
    for gal in galaxies:
        M_est = gal.estimate_mass(0.5)
        gamma_eff = gamma_0 + delta * np.log10(max(M_est, 1e6) / M_ref)
        gamma_eff = np.clip(gamma_eff, 0.05, 1.5)

        # Quick fit Yd
        best_c = 1e10
        for Yd in [0.3, 0.5, 0.8, 1.2]:
            gbar = gal.compute_gbar(Yd)
            mask = gbar > 0
            if np.sum(mask) < 3: continue
            y = gbar[mask] / a0
            try:
                nu = nu_exp(y, alpha, gamma_eff)
                gobs = gbar[mask] * nu
                vp = np.sqrt(np.abs(gobs) * gal.rad_m[mask]) / 1e3
                c2 = np.sum(((gal.vobs[mask] - vp) / gal.errv_safe[mask])**2)
                if c2 < best_c: best_c = c2
            except: pass
        total_chi2 += best_c

    # Cluster part
    for cl in clusters:
        M = cl.M_total_Msun
        gamma_eff = gamma_0 + delta * np.log10(max(M, 1e6) / M_ref)
        gamma_eff = np.clip(gamma_eff, 0.05, 1.5)
        y = cl.g_bar / a0
        nu = nu_exp(y, alpha, gamma_eff)
        g_pred = cl.g_bar * nu
        total_chi2 += cl.get_chi2(g_pred)

    return total_chi2


# ============================================================
# MODEL B: a0(R) = a0_gal * (1 + (R/R_c)^n)
# ============================================================
def model_B_chi2(params, galaxies, clusters, alpha=0.8, gamma=0.4):
    log_a0_gal, log_Rc_kpc, n_power = params
    a0_gal = 10**log_a0_gal
    Rc = 10**log_Rc_kpc  # kpc
    if n_power < 0.1 or n_power > 5: return 1e15
    if Rc < 1 or Rc > 1e5: return 1e15

    total_chi2 = 0

    for gal in galaxies:
        R_eff = gal.rad_kpc[-1]  # outermost radius in kpc
        a0_eff = a0_gal * (1.0 + (R_eff / Rc)**n_power)

        best_c = 1e10
        for Yd in [0.3, 0.5, 0.8, 1.2]:
            gbar = gal.compute_gbar(Yd)
            mask = gbar > 0
            if np.sum(mask) < 3: continue
            y = gbar[mask] / a0_eff
            try:
                nu = nu_exp(y, alpha, gamma)
                gobs = gbar[mask] * nu
                vp = np.sqrt(np.abs(gobs) * gal.rad_m[mask]) / 1e3
                c2 = np.sum(((gal.vobs[mask] - vp) / gal.errv_safe[mask])**2)
                if c2 < best_c: best_c = c2
            except: pass
        total_chi2 += best_c

    for cl in clusters:
        R_eff = cl.r500_m / kpc_m  # r500 in kpc
        a0_eff = a0_gal * (1.0 + (R_eff / Rc)**n_power)
        y = cl.g_bar / a0_eff
        nu = nu_exp(y, alpha, gamma)
        g_pred = cl.g_bar * nu
        total_chi2 += cl.get_chi2(g_pred)

    return total_chi2


# ============================================================
# MODEL C: Geometry factor — disk vs sphere
# ============================================================
def model_C_chi2(params, galaxies, clusters, alpha=0.8):
    gamma_disk, gamma_sphere, log_a0 = params
    a0 = 10**log_a0
    if gamma_disk < 0.05 or gamma_disk > 1.0: return 1e15
    if gamma_sphere < 0.05 or gamma_sphere > 1.5: return 1e15

    total_chi2 = 0

    for gal in galaxies:
        gamma = gamma_disk
        best_c = 1e10
        for Yd in [0.3, 0.5, 0.8, 1.2]:
            gbar = gal.compute_gbar(Yd)
            mask = gbar > 0
            if np.sum(mask) < 3: continue
            y = gbar[mask] / a0
            try:
                nu = nu_exp(y, alpha, gamma)
                gobs = gbar[mask] * nu
                vp = np.sqrt(np.abs(gobs) * gal.rad_m[mask]) / 1e3
                c2 = np.sum(((gal.vobs[mask] - vp) / gal.errv_safe[mask])**2)
                if c2 < best_c: best_c = c2
            except: pass
        total_chi2 += best_c

    for cl in clusters:
        gamma = gamma_sphere
        y = cl.g_bar / a0
        nu = nu_exp(y, alpha, gamma)
        g_pred = cl.g_bar * nu
        total_chi2 += cl.get_chi2(g_pred)

    return total_chi2


# ============================================================
# MODEL D: Second transition — nu = 1 + exp(-y^a)/y^g + A/y^g2
# ============================================================
def nu_double(y, alpha, gamma, A2, gamma2, y_break):
    """Two-component: standard Exp + extra boost below y_break."""
    y = np.maximum(y, 1e-15)
    nu1 = 1.0 + np.exp(-np.power(y, alpha)) / np.power(y, gamma)
    # Second component: smooth turn-on below y_break
    switch = 1.0 / (1.0 + np.power(y / y_break, 2))
    nu2 = A2 * switch / np.power(y, gamma2)
    return nu1 + nu2


def model_D_chi2(params, galaxies, clusters, alpha=0.8, gamma=0.4):
    log_a0, A2, gamma2, log_ybreak = params
    a0 = 10**log_a0
    y_break = 10**log_ybreak
    if A2 < 0 or A2 > 10: return 1e15
    if gamma2 < 0 or gamma2 > 1.0: return 1e15
    if y_break < 1e-4 or y_break > 1.0: return 1e15

    total_chi2 = 0

    for gal in galaxies:
        best_c = 1e10
        for Yd in [0.3, 0.5, 0.8, 1.2]:
            gbar = gal.compute_gbar(Yd)
            mask = gbar > 0
            if np.sum(mask) < 3: continue
            y = gbar[mask] / a0
            try:
                nu = nu_double(y, alpha, gamma, A2, gamma2, y_break)
                gobs = gbar[mask] * nu
                vp = np.sqrt(np.abs(gobs) * gal.rad_m[mask]) / 1e3
                c2 = np.sum(((gal.vobs[mask] - vp) / gal.errv_safe[mask])**2)
                if c2 < best_c: best_c = c2
            except: pass
        total_chi2 += best_c

    for cl in clusters:
        y = cl.g_bar / a0
        try:
            nu = nu_double(y, alpha, gamma, A2, gamma2, y_break)
        except: return 1e15
        g_pred = cl.g_bar * nu
        total_chi2 += cl.get_chi2(g_pred)

    return total_chi2


# ============================================================
# MODEL E: TGP transition scale H = sqrt(GM/(c*H0))
# ============================================================
def model_E_chi2(params, galaxies, clusters, alpha=0.8):
    gamma_base, gamma_scale, log_a0 = params
    a0 = 10**log_a0
    if gamma_base < 0.1 or gamma_base > 0.8: return 1e15
    if gamma_scale < -0.3 or gamma_scale > 0.3: return 1e15

    total_chi2 = 0

    for gal in galaxies:
        M = gal.estimate_mass(0.5) * M_sun  # kg
        r_out = gal.rad_m[-1]
        # TGP transition scale
        H_trans = np.sqrt(G_SI * M / (c_light * H0))
        ratio = r_out / H_trans
        # gamma increases with r/H (deeper into transition)
        gamma_eff = gamma_base + gamma_scale * np.log10(max(ratio, 0.01))
        gamma_eff = np.clip(gamma_eff, 0.05, 1.2)

        best_c = 1e10
        for Yd in [0.3, 0.5, 0.8, 1.2]:
            gbar = gal.compute_gbar(Yd)
            mask = gbar > 0
            if np.sum(mask) < 3: continue
            y = gbar[mask] / a0
            try:
                nu = nu_exp(y, alpha, gamma_eff)
                gobs = gbar[mask] * nu
                vp = np.sqrt(np.abs(gobs) * gal.rad_m[mask]) / 1e3
                c2 = np.sum(((gal.vobs[mask] - vp) / gal.errv_safe[mask])**2)
                if c2 < best_c: best_c = c2
            except: pass
        total_chi2 += best_c

    for cl in clusters:
        M = cl.M500
        r_out = cl.r500_m
        H_trans = np.sqrt(G_SI * M / (c_light * H0))
        ratio = r_out / H_trans
        gamma_eff = gamma_base + gamma_scale * np.log10(max(ratio, 0.01))
        gamma_eff = np.clip(gamma_eff, 0.05, 1.2)

        y = cl.g_bar / a0
        nu = nu_exp(y, alpha, gamma_eff)
        g_pred = cl.g_bar * nu
        total_chi2 += cl.get_chi2(g_pred)

    return total_chi2


# ============================================================
# BASELINE: fixed gamma=0.4
# ============================================================
def baseline_chi2(galaxies, clusters, alpha=0.8, gamma=0.4, a0=a0_ref):
    total = 0
    for gal in galaxies:
        best_c = 1e10
        for Yd in [0.3, 0.5, 0.8, 1.2]:
            gbar = gal.compute_gbar(Yd)
            mask = gbar > 0
            if np.sum(mask) < 3: continue
            y = gbar[mask] / a0
            try:
                nu = nu_exp(y, alpha, gamma)
                gobs = gbar[mask] * nu
                vp = np.sqrt(np.abs(gobs) * gal.rad_m[mask]) / 1e3
                c2 = np.sum(((gal.vobs[mask] - vp) / gal.errv_safe[mask])**2)
                if c2 < best_c: best_c = c2
            except: pass
        total += best_c

    for cl in clusters:
        y = cl.g_bar / a0
        nu = nu_exp(y, alpha, gamma)
        g_pred = cl.g_bar * nu
        total += cl.get_chi2(g_pred)

    return total


# ============================================================
# MAIN
# ============================================================
def main():
    print("=" * 80)
    print("gs16: UNIFIED MODEL — GALAXY + CLUSTER")
    print("Can one model explain both scales?")
    print("=" * 80)

    data_dir = os.path.join(os.path.dirname(__file__), 'Rotmod_LTG')
    all_galaxies = load_galaxies(data_dir)
    clusters = build_clusters()

    # Subsample galaxies for speed: every 4th galaxy (sorted by name)
    galaxies = all_galaxies[::4]
    print(f"\nLoaded {len(all_galaxies)} galaxies, using {len(galaxies)} subsample + {len(clusters)} clusters")
    sys.stdout.flush()

    # ================================================================
    # TRANSITION SCALES
    # ================================================================
    print("\n" + "=" * 80)
    print("SECTION 0: TGP TRANSITION SCALES")
    print("=" * 80)

    print(f"\n  H = sqrt(GM / (c*H0)) = sqrt(r_S * r_H)")
    print(f"  r_H = c/H0 = {c_light/H0/Mpc_m:.0f} Mpc")
    print(f"\n  {'System':<20s} {'M (Msun)':<12s} {'R_out (kpc)':<12s} "
          f"{'H (kpc)':<10s} {'R/H':<8s}")
    print("  " + "-" * 62)

    test_systems = [
        ("Dwarf (DDO154)", 1e9, 8),
        ("MW-like", 5e10, 20),
        ("NGC2841", 3e11, 50),
        ("Fornax group", 2.8e13, 450),
        ("Virgo cluster", 1.4e14, 770),
        ("Coma cluster", 7.2e14, 1310),
        ("A2390 massive", 1.18e15, 1500),
    ]
    for name, M, R_kpc in test_systems:
        H_m = np.sqrt(G_SI * M * M_sun / (c_light * H0))
        H_kpc = H_m / kpc_m
        ratio = R_kpc / H_kpc
        print(f"  {name:<20s} {M:<12.2e} {R_kpc:<12.0f} {H_kpc:<10.1f} {ratio:<8.2f}")

    # ================================================================
    # BASELINES
    # ================================================================
    sys.stdout.flush()
    print("\n" + "=" * 80)
    print("SECTION 1: BASELINES")
    print("=" * 80)
    sys.stdout.flush()

    chi2_base_04 = baseline_chi2(galaxies, clusters, 0.8, 0.4, a0_ref)
    chi2_base_05 = baseline_chi2(galaxies, clusters, 0.8, 0.5, a0_ref)

    n_gal = sum(g.npts for g in galaxies)
    n_clust = len(clusters)
    n_total = n_gal + n_clust
    print(f"\n  Total data points: {n_total} ({n_gal} galaxy + {n_clust} cluster)")
    print(f"  Baseline chi2 (gamma=0.4): {chi2_base_04:.1f} (chi2/dof = {chi2_base_04/n_total:.3f})")
    print(f"  Baseline chi2 (gamma=0.5): {chi2_base_05:.1f} (chi2/dof = {chi2_base_05/n_total:.3f})")

    # Decompose
    chi2_gal_04 = baseline_chi2(galaxies, [], 0.8, 0.4, a0_ref)
    chi2_gal_05 = baseline_chi2(galaxies, [], 0.8, 0.5, a0_ref)
    chi2_cl_04 = baseline_chi2([], clusters, 0.8, 0.4, a0_ref)
    chi2_cl_05 = baseline_chi2([], clusters, 0.8, 0.5, a0_ref)

    print(f"\n  Decomposed:")
    print(f"  {'Component':<20s} {'chi2(g=0.4)':<15s} {'chi2(g=0.5)':<15s} {'winner'}")
    print(f"  {'Galaxies':<20s} {chi2_gal_04:<15.1f} {chi2_gal_05:<15.1f} "
          f"{'g=0.4' if chi2_gal_04 < chi2_gal_05 else 'g=0.5'}")
    print(f"  {'Clusters':<20s} {chi2_cl_04:<15.1f} {chi2_cl_05:<15.1f} "
          f"{'g=0.4' if chi2_cl_04 < chi2_cl_05 else 'g=0.5'}")
    print(f"  {'Total':<20s} {chi2_base_04:<15.1f} {chi2_base_05:<15.1f} "
          f"{'g=0.4' if chi2_base_04 < chi2_base_05 else 'g=0.5'}")

    # ================================================================
    # MODEL A: gamma(M)
    # ================================================================
    print("\n" + "=" * 80)
    print("MODEL A: gamma(M) = gamma_0 + delta * log10(M / 10^10 Msun)")
    print("=" * 80)

    sys.stdout.flush()
    best_A = {'chi2': 1e15}
    for g0 in [0.35, 0.45]:
        for d in [-0.05, 0.0, 0.05]:
            for la0 in [-10.0, -9.9]:
                try:
                    res = minimize(model_A_chi2, [g0, d, la0],
                                  args=(galaxies, clusters),
                                  method='Nelder-Mead',
                                  options={'maxiter': 3000, 'xatol': 1e-4})
                    if res.fun < best_A['chi2']:
                        best_A = {'chi2': res.fun, 'params': res.x}
                except: pass

    if 'params' in best_A:
        g0, d, la0 = best_A['params']
        print(f"\n  Best fit: gamma_0 = {g0:.4f}, delta = {d:.4f}, a0 = {10**la0:.3e}")
        print(f"  chi2 = {best_A['chi2']:.1f} (vs baseline {chi2_base_04:.1f})")
        print(f"  dchi2 = {chi2_base_04 - best_A['chi2']:.1f}")
        print(f"\n  Implied gamma at different masses:")
        for M in [1e8, 1e9, 1e10, 1e11, 1e12, 1e13, 1e14, 1e15]:
            ge = g0 + d * np.log10(M / 1e10)
            print(f"    M = {M:.0e}: gamma = {ge:.3f}")

    # ================================================================
    # MODEL B: a0(R)
    # ================================================================
    print("\n" + "=" * 80)
    print("MODEL B: a0(R) = a0_gal * (1 + (R/R_c)^n)")
    print("=" * 80)

    sys.stdout.flush()
    best_B = {'chi2': 1e15}
    for la0 in [-10.0, -9.95]:
        for lrc in [2.5, 3.0]:
            for n in [0.5, 1.0, 2.0]:
                try:
                    res = minimize(model_B_chi2, [la0, lrc, n],
                                  args=(galaxies, clusters),
                                  method='Nelder-Mead',
                                  options={'maxiter': 3000, 'xatol': 1e-4})
                    if res.fun < best_B['chi2']:
                        best_B = {'chi2': res.fun, 'params': res.x}
                except: pass

    if 'params' in best_B:
        la0, lrc, n = best_B['params']
        print(f"\n  Best fit: a0_gal = {10**la0:.3e}, R_c = {10**lrc:.0f} kpc, n = {n:.3f}")
        print(f"  chi2 = {best_B['chi2']:.1f} (vs baseline {chi2_base_04:.1f})")
        print(f"  dchi2 = {chi2_base_04 - best_B['chi2']:.1f}")
        print(f"\n  Implied a0 at different scales:")
        for R in [5, 20, 50, 100, 500, 1000, 1500]:
            a0_eff = 10**la0 * (1 + (R / 10**lrc)**n)
            print(f"    R = {R:5d} kpc: a0 = {a0_eff:.3e} ({a0_eff/a0_ref:.2f}x galaxy)")

    # ================================================================
    # MODEL C: Geometry (disk vs sphere)
    # ================================================================
    print("\n" + "=" * 80)
    print("MODEL C: Separate gamma for DISKS (galaxies) vs SPHERES (clusters)")
    print("=" * 80)

    sys.stdout.flush()
    best_C = {'chi2': 1e15}
    for gd in [0.35, 0.4]:
        for gs in [0.5, 0.6, 0.7]:
            for la0 in [-10.0, -9.9, -9.5]:
                try:
                    res = minimize(model_C_chi2, [gd, gs, la0],
                                  args=(galaxies, clusters),
                                  method='Nelder-Mead',
                                  options={'maxiter': 3000, 'xatol': 1e-4})
                    if res.fun < best_C['chi2']:
                        best_C = {'chi2': res.fun, 'params': res.x}
                except: pass

    if 'params' in best_C:
        gd, gs, la0 = best_C['params']
        print(f"\n  Best fit: gamma_disk = {gd:.4f}, gamma_sphere = {gs:.4f}, a0 = {10**la0:.3e}")
        print(f"  chi2 = {best_C['chi2']:.1f} (vs baseline {chi2_base_04:.1f})")
        print(f"  dchi2 = {chi2_base_04 - best_C['chi2']:.1f}")
        print(f"  d_eff (disk)   = {3 - 2*gd:.2f}")
        print(f"  d_eff (sphere) = {3 - 2*gs:.2f}")

    # ================================================================
    # MODEL D: Second transition
    # ================================================================
    print("\n" + "=" * 80)
    print("MODEL D: Second transition — extra boost at very low y")
    print("nu = 1 + exp(-y^a)/y^g + A2 * switch(y/y_break) / y^g2")
    print("=" * 80)

    sys.stdout.flush()
    best_D = {'chi2': 1e15}
    for la0 in [-10.0, -9.95]:
        for A2 in [0.3, 1.0, 2.0]:
            for g2 in [0.2, 0.4]:
                for lyb in [-1.5, -1.0]:
                    try:
                        res = minimize(model_D_chi2, [la0, A2, g2, lyb],
                                      args=(galaxies, clusters),
                                      method='Nelder-Mead',
                                      options={'maxiter': 3000, 'xatol': 1e-4})
                        if res.fun < best_D['chi2']:
                            best_D = {'chi2': res.fun, 'params': res.x}
                    except: pass

    if 'params' in best_D:
        la0, A2, g2, lyb = best_D['params']
        print(f"\n  Best fit: a0 = {10**la0:.3e}, A2 = {A2:.4f}, gamma2 = {g2:.4f}, "
              f"y_break = {10**lyb:.4f}")
        print(f"  chi2 = {best_D['chi2']:.1f} (vs baseline {chi2_base_04:.1f})")
        print(f"  dchi2 = {chi2_base_04 - best_D['chi2']:.1f}")
        print(f"\n  nu values with double transition:")
        for y_test in [0.001, 0.01, 0.05, 0.1, 0.5, 1.0, 5.0]:
            nu_std = nu_exp(y_test, 0.8, 0.4)
            nu_dbl = nu_double(y_test, 0.8, 0.4, A2, g2, 10**lyb)
            print(f"    y={y_test:<6.3f}: nu_std={nu_std:.3f}, "
                  f"nu_double={nu_dbl:.3f}, ratio={nu_dbl/nu_std:.3f}")

    # ================================================================
    # MODEL E: TGP r/H dependent gamma
    # ================================================================
    print("\n" + "=" * 80)
    print("MODEL E: gamma depends on r/H where H = sqrt(GM/(c*H0))")
    print("gamma(r/H) = gamma_base + gamma_scale * log10(r/H)")
    print("=" * 80)

    sys.stdout.flush()
    best_E = {'chi2': 1e15}
    for gb in [0.35, 0.45]:
        for gsc in [-0.05, 0.05, 0.1]:
            for la0 in [-10.0, -9.9]:
                try:
                    res = minimize(model_E_chi2, [gb, gsc, la0],
                                  args=(galaxies, clusters),
                                  method='Nelder-Mead',
                                  options={'maxiter': 3000, 'xatol': 1e-4})
                    if res.fun < best_E['chi2']:
                        best_E = {'chi2': res.fun, 'params': res.x}
                except: pass

    if 'params' in best_E:
        gb, gsc, la0 = best_E['params']
        print(f"\n  Best fit: gamma_base = {gb:.4f}, gamma_scale = {gsc:.4f}, a0 = {10**la0:.3e}")
        print(f"  chi2 = {best_E['chi2']:.1f} (vs baseline {chi2_base_04:.1f})")
        print(f"  dchi2 = {chi2_base_04 - best_E['chi2']:.1f}")
        print(f"\n  Implied gamma for test systems:")
        for name, M, R_kpc in test_systems:
            H_m = np.sqrt(G_SI * M * M_sun / (c_light * H0))
            H_kpc = H_m / kpc_m
            ratio = R_kpc / H_kpc
            ge = gb + gsc * np.log10(max(ratio, 0.01))
            ge = np.clip(ge, 0.05, 1.2)
            print(f"    {name:<20s}: r/H = {ratio:.2f}, gamma = {ge:.3f}")

    # ================================================================
    # COMPARISON
    # ================================================================
    print("\n" + "=" * 80)
    print("MODEL COMPARISON")
    print("=" * 80)

    results = [
        ("Baseline (g=0.4)", chi2_base_04, 1),  # a0 only
        ("Baseline (g=0.5)", chi2_base_05, 1),
    ]
    if 'params' in best_A:
        results.append(("A: gamma(M)", best_A['chi2'], 3))
    if 'params' in best_B:
        results.append(("B: a0(R)", best_B['chi2'], 3))
    if 'params' in best_C:
        results.append(("C: disk vs sphere", best_C['chi2'], 3))
    if 'params' in best_D:
        results.append(("D: second transition", best_D['chi2'], 4))
    if 'params' in best_E:
        results.append(("E: gamma(r/H)", best_E['chi2'], 3))

    print(f"\n  {'Model':<25s} {'chi2':<12s} {'n_params':<10s} {'dchi2 vs base':<14s} "
          f"{'BIC penalty':<12s} {'net dBIC':<10s}")
    print("  " + "-" * 83)
    for name, chi2, npar in sorted(results, key=lambda x: x[1]):
        dchi2 = chi2_base_04 - chi2
        bic_pen = (npar - 1) * np.log(n_total)
        net_dbic = dchi2 - bic_pen
        print(f"  {name:<25s} {chi2:<12.1f} {npar:<10d} {dchi2:<+14.1f} "
              f"{bic_pen:<12.1f} {net_dbic:<+10.1f}")

    # ================================================================
    # CLUSTER PREDICTIONS WITH BEST MODEL
    # ================================================================
    print("\n" + "=" * 80)
    print("CLUSTER PREDICTIONS WITH BEST UNIFIED MODEL")
    print("=" * 80)

    # Use Model C (geometry) as it's most physically motivated
    if 'params' in best_C:
        gd, gs, la0 = best_C['params']
        a0 = 10**la0
        print(f"\n  Using Model C: gamma_disk={gd:.3f}, gamma_sphere={gs:.3f}, a0={a0:.3e}")
        print(f"\n  {'Cluster':<14s} {'g_obs/g_pred':<12s} {'status'}")
        print("  " + "-" * 40)
        for cl in clusters:
            y = cl.g_bar / a0
            nu = nu_exp(y, 0.8, gs)
            g_pred = cl.g_bar * nu
            ratio = cl.g_obs / g_pred
            status = "OK" if 0.8 <= ratio <= 1.2 else ("close" if 0.6 <= ratio <= 1.4 else "BAD")
            print(f"  {cl.name:<14s} {ratio:<12.3f} {status}")

    # ================================================================
    # PHYSICAL INTERPRETATION
    # ================================================================
    print("\n" + "=" * 80)
    print("PHYSICAL INTERPRETATION FOR TGP")
    print("=" * 80)

    if 'params' in best_C:
        gd, gs, la0 = best_C['params']
        print(f"""
  GEOMETRY HYPOTHESIS:

  Disk galaxies (gamma = {gd:.3f}):
    d_eff = {3-2*gd:.2f} -> {(1-gd/0.5)*100:.0f}% of full 3D->2D transition
    Explanation: disk geometry already constrains mass distribution
    to ~2D -> transition parameter is SMALLER

  Spherical clusters (gamma = {gs:.3f}):
    d_eff = {3-2*gs:.2f} -> deeper transition
    Explanation: full 3D system needs MORE dimensional reduction
    to produce flat "rotation curves" at cluster scale

  This makes TGP sense:
    - In a disk, matter is already partially 2D
    - The substrate deformation is "pre-configured" for 2D
    - Less gamma needed to reach the modified regime
    - In a sphere, the substrate must do MORE work
    - Full 3D->2D transition requires larger gamma

  Testable prediction:
    - Elliptical galaxies (spheroidal) should prefer gamma closer to {gs:.2f}
    - Edge-on vs face-on shouldn't matter (it's 3D geometry, not projection)
    - Galaxy groups (intermediate) should have intermediate gamma
""")

    print("=" * 80)
    print("DONE")
    print("=" * 80)


if __name__ == '__main__':
    main()
