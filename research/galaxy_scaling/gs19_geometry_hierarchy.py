#!/usr/bin/env python3
"""
gs19: UNIFIED GEOMETRY MODEL — FULL MORPHOLOGICAL HIERARCHY
============================================================
gs16 found: disk (gamma=0.36) vs sphere (gamma=0.54) explains galaxy+cluster gap
gs17 showed: R_eff data inconclusive (y too high)
gs18 showed: at large radii, ellipticals prefer gamma=0.6-0.8 (>0.54, >>0.36)
             slow-rotators > fast-rotators — geometry matters!

THIS SCRIPT: Build unified model gamma(geometry) across:
  disk galaxies (SPARC) → S0/fast E → slow E → groups → clusters

Tests:
  A) gamma as function of 3D shape (axis ratio / Sersic n)
  B) gamma as continuous function of "sphericity" S = 0 (disk) to 1 (sphere)
  C) Combined: gamma(S) + mass dependence
  D) Three-class: disk / intermediate / sphere
  E) Comparison with lensing-based RAR (if available)
"""

import sys, io, os
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
H0 = 70e3 / Mpc_m
a0_ref = 1.12e-10

# ============================================================
# GALAXY LOADER (from gs16)
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
        self.geometry = 'disk'
        self.sphericity = 0.0  # 0 = pure disk

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
# CLUSTER DATA
# ============================================================
class ClusterPoint:
    def __init__(self, name, r500_Mpc, M500_1e14, fbar):
        self.name = name
        self.geometry = 'sphere'
        self.sphericity = 1.0
        self.npts = 1
        self.has_bulge = False

        self.r500_m = r500_Mpc * Mpc_m
        self.M500 = M500_1e14 * 1e14 * M_sun
        self.Mbar = fbar * self.M500
        self.M_total_Msun = M500_1e14 * 1e14
        self.v_obs = np.sqrt(G_SI * self.M500 / self.r500_m) / 1e3
        self.v_err = 0.05 * self.v_obs
        self.g_bar = G_SI * self.Mbar / self.r500_m**2
        self.g_obs = G_SI * self.M500 / self.r500_m**2

    def get_chi2(self, g_pred):
        v_pred = np.sqrt(abs(g_pred) * self.r500_m) / 1e3
        return ((v_pred - self.v_obs) / self.v_err)**2


def build_clusters():
    data = [
        ("A0122",0.89,2.26,0.112),("A1651",1.18,5.15,0.143),
        ("A2401",0.68,0.95,0.118),("A2721",1.03,3.46,0.142),
        ("A2811",1.04,3.59,0.138),("A2955",0.68,0.99,0.097),
        ("A2984",0.67,0.95,0.152),("A3112",1.02,3.23,0.154),
        ("A3693",0.90,2.26,0.133),("A4010",0.92,2.41,0.142),
        ("AS0084",0.91,2.37,0.110),("AS0296",0.78,1.45,0.095),
        ("Coma",1.31,7.2,0.164),("Perseus",1.27,6.6,0.180),
        ("Virgo",0.77,1.4,0.120),("Centaurus",0.82,1.8,0.125),
        ("A2029",1.42,8.71,0.158),("A0478",1.28,6.58,0.195),
        ("A2390",1.50,11.8,0.161),("A1795",1.09,4.54,0.138),
        ("A2142",1.37,8.17,0.194),("A2319",1.55,12.0,0.168),
        ("A644",1.05,3.70,0.150),("A3266",1.25,6.20,0.158),
        ("Fornax",0.45,0.28,0.095),("MKW4",0.50,0.38,0.095),
        ("NGC5044g",0.43,0.24,0.095),("A3571",1.13,5.02,0.123),
        ("Bullet",1.35,7.8,0.155),
    ]
    return [ClusterPoint(n,r,m,f) for n,r,m,f in data]


# ============================================================
# ELLIPTICAL GALAXY DATA (from gs18 — extended kinematics)
# ============================================================
class EllipticalPoint:
    """Elliptical galaxy data point at radius r from gs18."""
    def __init__(self, name, gtype, r_kpc, g_bar, g_obs, M_star):
        self.name = name
        self.gtype = gtype  # 'slow' or 'fast'
        self.r_kpc = r_kpc
        self.g_bar = g_bar
        self.g_obs = g_obs
        self.M_star = M_star
        self.npts = 1
        self.has_bulge = False

        # Sphericity: slow rotators = more spherical
        if gtype == 'slow':
            self.sphericity = 0.85
            self.geometry = 'sphere'
        else:
            self.sphericity = 0.55  # fast rotators = intermediate
            self.geometry = 'intermediate'

    def get_chi2(self, g_pred):
        """Chi2 from acceleration comparison."""
        # 20% fractional error on g_obs (anisotropy + M/L)
        sigma = 0.20 * self.g_obs
        return ((g_pred - self.g_obs) / sigma)**2


def build_ellipticals():
    """Build elliptical data from gs18 (Jeans K=2.5 for radial anisotropy)."""

    def hernquist_mass(r_kpc, M_star, R_eff):
        a = R_eff / 1.8153
        return M_star * r_kpc**2 / (r_kpc + a)**2

    def g_bar_at_r(r_kpc, M_star, R_eff):
        M_enc = hernquist_mass(r_kpc, M_star, R_eff)
        r_m = r_kpc * kpc_m
        return G_SI * M_enc * M_sun / r_m**2

    def g_obs_sigma(sigma_kms, r_kpc, K=2.5):
        sigma_ms = sigma_kms * 1e3
        r_m = r_kpc * kpc_m
        return K * sigma_ms**2 / r_m

    def g_obs_mass(M_total, r_kpc):
        r_m = r_kpc * kpc_m
        return G_SI * M_total * M_sun / r_m**2

    # Galaxy data with sigma profiles and X-ray mass
    gal_data = [
        {'name':'NGC4636','type':'slow','M_star':2.5e11,'R_eff':11.4,
         'sigma':[(0.5,215),(1.0,209),(2.0,185),(3.0,170),(4.0,160),(5.0,155)],
         'xray':[(20,2.2e11),(30,4.5e11),(40,7.5e11),(50,1.1e12),(60,1.5e12)]},

        {'name':'NGC4486','type':'slow','M_star':6.0e11,'R_eff':8.0,
         'sigma':[(1.0,375),(2.0,350),(4.0,330),(7.0,310),(10.0,300),(15.0,290)],
         'xray':[(30,1.3e12),(50,2.5e12),(80,4.5e12),(100,6.0e12),(150,1.0e13)]},

        {'name':'NGC1399','type':'slow','M_star':3.5e11,'R_eff':6.5,
         'sigma':[(1.0,337),(2.0,310),(4.0,270),(7.0,250),(10.0,245)],
         'xray':[(20,5.5e11),(30,1.0e12),(50,2.0e12),(80,4.0e12)]},

        {'name':'NGC5846','type':'slow','M_star':2.8e11,'R_eff':7.8,
         'sigma':[(1.0,238),(2.0,215),(3.0,200),(4.0,190),(5.0,185)],
         'xray':[(20,3.0e11),(30,5.5e11),(50,1.2e12)]},

        {'name':'NGC4472','type':'slow','M_star':5.0e11,'R_eff':9.4,
         'sigma':[(1.0,295),(2.0,268),(4.0,240),(5.0,235),(6.0,230)],
         'xray':[(20,6.0e11),(40,1.5e12),(60,2.5e12)]},

        {'name':'NGC4649','type':'slow','M_star':4.5e11,'R_eff':7.2,
         'sigma':[(1.0,335),(2.0,300),(4.0,265),(5.0,255)],
         'xray':[(20,5.5e11),(40,1.3e12)]},

        {'name':'NGC1407','type':'slow','M_star':3.0e11,'R_eff':9.0,
         'sigma':[(1.0,270),(2.0,248),(4.0,225),(7.0,215),(10.0,210)],
         'xray':[(20,4.0e11),(40,1.0e12),(60,1.8e12)]},

        {'name':'NGC3379','type':'fast','M_star':1.0e11,'R_eff':3.3,
         'sigma':[(1.0,206),(2.0,175),(4.0,140),(6.0,125),(7.0,120)]},

        {'name':'NGC4374','type':'slow','M_star':4.0e11,'R_eff':7.5,
         'sigma':[(1.0,296),(2.0,268),(4.0,245),(5.0,240)]},

        {'name':'NGC4406','type':'slow','M_star':3.5e11,'R_eff':12.5,
         'sigma':[(1.0,235),(2.0,210),(4.0,185),(5.0,180),(6.0,178)]},

        {'name':'NGC4278','type':'fast','M_star':1.2e11,'R_eff':2.3,
         'sigma':[(1.0,252),(2.0,228),(4.0,200),(5.0,192)]},

        {'name':'NGC5813','type':'slow','M_star':2.2e11,'R_eff':6.0,
         'sigma':[(1.0,230),(2.0,208),(4.0,188),(5.0,183)],
         'xray':[(20,2.5e11),(40,7.0e11)]},

        {'name':'NGC5044','type':'slow','M_star':2.0e11,'R_eff':5.5,
         'sigma':[(1.0,235),(2.0,210),(4.0,190),(5.0,185)],
         'xray':[(20,2.8e11),(30,5.0e11),(50,1.1e12)]},
    ]

    points = []
    for gal in gal_data:
        M_star = gal['M_star']
        R_eff = gal['R_eff']

        # Kinematic points (use K=2.5 — mild radial anisotropy)
        for r_Re, sigma in gal['sigma']:
            r_kpc = r_Re * R_eff
            gb = g_bar_at_r(r_kpc, M_star, R_eff)
            go = g_obs_sigma(sigma, r_kpc, K=2.5)
            y = gb / a0_ref
            if y < 5.0:  # only transition regime
                points.append(EllipticalPoint(gal['name'], gal['type'],
                                              r_kpc, gb, go, M_star))

        # X-ray points (most reliable)
        if 'xray' in gal:
            for r_kpc, M_tot in gal['xray']:
                gb = g_bar_at_r(r_kpc, M_star, R_eff)
                go = g_obs_mass(M_tot, r_kpc)
                y = gb / a0_ref
                if y < 5.0:
                    points.append(EllipticalPoint(gal['name'], gal['type'],
                                                  r_kpc, gb, go, M_star))

    return points


# ============================================================
# INTERPOLATION
# ============================================================
def nu_exp(y, gamma, alpha=0.8):
    y = np.maximum(y, 1e-15)
    return 1.0 + np.exp(-np.power(y, alpha)) / np.power(y, gamma)


# ============================================================
# FITTING FUNCTIONS
# ============================================================
def galaxy_chi2(gal, gamma, a0, Yd_grid=np.arange(0.1, 2.01, 0.05)):
    """Compute best chi2 for a single SPARC galaxy with given gamma, a0."""
    gbar_ref = gal.compute_gbar(1.0)
    best_chi2 = 1e30
    for Yd in Yd_grid:
        gbar = gbar_ref * Yd if not gal.has_bulge else gal.compute_gbar(Yd)
        y = np.abs(gbar) / a0
        gpred = gbar * nu_exp(y, gamma)
        vpred = np.sign(gpred) * np.sqrt(np.abs(gpred) * gal.rad_m) / 1e3
        chi2 = np.sum(((vpred - gal.vobs) / gal.errv_safe)**2)
        if chi2 < best_chi2:
            best_chi2 = chi2
    return best_chi2


def cluster_chi2(cl, gamma, a0):
    """Chi2 for a single cluster."""
    y = cl.g_bar / a0
    gpred = cl.g_bar * nu_exp(y, gamma)
    return cl.get_chi2(gpred)


def elliptical_chi2(ep, gamma, a0):
    """Chi2 for a single elliptical data point."""
    y = ep.g_bar / a0
    gpred = ep.g_bar * nu_exp(y, gamma)
    return ep.get_chi2(gpred)


# ============================================================
# MAIN
# ============================================================
print("=" * 80)
print("gs19: UNIFIED GEOMETRY MODEL — FULL MORPHOLOGICAL HIERARCHY")
print("gamma(geometry) across disk -> S0 -> E -> group -> cluster")
print("=" * 80)
sys.stdout.flush()

# Load data
data_dir = os.path.join(os.path.dirname(__file__), 'Rotmod_LTG')
all_galaxies = load_galaxies(data_dir)
# Subsample for speed: every 4th galaxy
galaxies = all_galaxies[::4]
clusters = build_clusters()
ellipticals = build_ellipticals()

print(f"\nLoaded: {len(all_galaxies)} SPARC galaxies (using {len(galaxies)} subsample)")
print(f"        {len(ellipticals)} elliptical data points (13 galaxies)")
print(f"        {len(clusters)} clusters")
sys.stdout.flush()

# ============================================================
# SECTION 0: BASELINE — UNIVERSAL GAMMA
# ============================================================
print("\n" + "=" * 80)
print("SECTION 0: BASELINE — SINGLE UNIVERSAL GAMMA")
print("=" * 80)
sys.stdout.flush()

Yd_grid = np.arange(0.1, 2.01, 0.1)

for gamma_test in [0.36, 0.40, 0.50, 0.54]:
    chi2_gal = sum(galaxy_chi2(g, gamma_test, a0_ref, Yd_grid) for g in galaxies)
    chi2_ell = sum(elliptical_chi2(e, gamma_test, a0_ref) for e in ellipticals)
    chi2_cl = sum(cluster_chi2(c, gamma_test, a0_ref) for c in clusters)
    total = chi2_gal + chi2_ell + chi2_cl
    print(f"  gamma={gamma_test:.2f}: gal={chi2_gal:.0f} ell={chi2_ell:.0f} cl={chi2_cl:.0f} total={total:.0f}")
    sys.stdout.flush()

# ============================================================
# SECTION 1: MODEL A — GAMMA(SPHERICITY)
# gamma = gamma_disk + (gamma_sphere - gamma_disk) * S
# S = 0 for disks, S = 0.55 for fast E, S = 0.85 for slow E, S = 1 for clusters
# ============================================================
print("\n" + "=" * 80)
print("SECTION 1: MODEL A — gamma(sphericity)")
print("gamma = gamma_disk + (gamma_sphere - gamma_disk) * S")
print("=" * 80)
sys.stdout.flush()

def model_A_total_chi2(params):
    gamma_disk, gamma_sphere, log_a0 = params
    a0 = 10**log_a0
    if gamma_disk < 0.05 or gamma_disk > 1.0: return 1e15
    if gamma_sphere < 0.05 or gamma_sphere > 1.5: return 1e15
    if gamma_sphere < gamma_disk: return 1e15

    total = 0
    # Disk galaxies (S=0)
    for g in galaxies:
        total += galaxy_chi2(g, gamma_disk, a0, Yd_grid)
    # Ellipticals
    for e in ellipticals:
        gamma_e = gamma_disk + (gamma_sphere - gamma_disk) * e.sphericity
        total += elliptical_chi2(e, gamma_e, a0)
    # Clusters (S=1)
    for c in clusters:
        total += cluster_chi2(c, gamma_sphere, a0)
    return total

# Multi-start optimization
best_A = None
best_A_chi2 = 1e30
for gd in [0.30, 0.36, 0.40]:
    for gs in [0.50, 0.60, 0.70, 0.80]:
        for la0 in [-10.0, -9.95, -9.85]:
            try:
                res = minimize(model_A_total_chi2, [gd, gs, la0],
                             method='Nelder-Mead',
                             options={'maxiter': 3000, 'xatol': 1e-3, 'fatol': 1})
                if res.fun < best_A_chi2:
                    best_A_chi2 = res.fun
                    best_A = res.x
            except: pass

if best_A is not None:
    gd, gs, la0 = best_A
    a0_A = 10**la0
    print(f"\n  Best fit:")
    print(f"    gamma_disk    = {gd:.4f}")
    print(f"    gamma_sphere  = {gs:.4f}")
    print(f"    a0            = {a0_A:.3e}")
    print(f"    chi2          = {best_A_chi2:.1f}")

    # Implied gamma for each morphological class
    print(f"\n  Implied gamma by morphology:")
    for label, S in [('Pure disk (SPARC)', 0.0),
                     ('S0/fast E (S=0.55)', 0.55),
                     ('Slow E (S=0.85)', 0.85),
                     ('Cluster (S=1.0)', 1.0)]:
        gamma = gd + (gs - gd) * S
        d_eff = 3 - 2*gamma
        print(f"    {label:<25} S={S:.2f}  gamma={gamma:.3f}  d_eff={d_eff:.2f}")

    # Decomposed chi2
    chi2_gal = sum(galaxy_chi2(g, gd, a0_A, Yd_grid) for g in galaxies)
    chi2_ell_slow = sum(elliptical_chi2(e, gd+(gs-gd)*0.85, a0_A)
                       for e in ellipticals if e.gtype == 'slow')
    chi2_ell_fast = sum(elliptical_chi2(e, gd+(gs-gd)*0.55, a0_A)
                       for e in ellipticals if e.gtype == 'fast')
    chi2_cl = sum(cluster_chi2(c, gs, a0_A) for c in clusters)
    print(f"\n  Decomposed chi2:")
    print(f"    Galaxies (disk):     {chi2_gal:.0f}")
    print(f"    Ellipticals (slow):  {chi2_ell_slow:.0f}")
    print(f"    Ellipticals (fast):  {chi2_ell_fast:.0f}")
    print(f"    Clusters:            {chi2_cl:.0f}")
    print(f"    Total:               {chi2_gal + chi2_ell_slow + chi2_ell_fast + chi2_cl:.0f}")

sys.stdout.flush()

# ============================================================
# SECTION 2: MODEL B — THREE-CLASS (disk / intermediate / sphere)
# ============================================================
print("\n" + "=" * 80)
print("SECTION 2: MODEL B — THREE-CLASS (disk / intermediate / sphere)")
print("=" * 80)
sys.stdout.flush()

def model_B_total_chi2(params):
    gamma_d, gamma_i, gamma_s, log_a0 = params
    a0 = 10**log_a0
    if gamma_d < 0.05 or gamma_d > 1.0: return 1e15
    if gamma_i < 0.05 or gamma_i > 1.0: return 1e15
    if gamma_s < 0.05 or gamma_s > 1.5: return 1e15

    total = 0
    for g in galaxies:
        total += galaxy_chi2(g, gamma_d, a0, Yd_grid)
    for e in ellipticals:
        if e.gtype == 'fast':
            total += elliptical_chi2(e, gamma_i, a0)
        else:
            total += elliptical_chi2(e, gamma_s, a0)
    for c in clusters:
        total += cluster_chi2(c, gamma_s, a0)
    return total

best_B = None
best_B_chi2 = 1e30
for gd in [0.30, 0.36, 0.40]:
    for gi in [0.40, 0.50, 0.60]:
        for gs in [0.50, 0.60, 0.70, 0.80]:
            for la0 in [-10.0, -9.95, -9.85]:
                try:
                    res = minimize(model_B_total_chi2, [gd, gi, gs, la0],
                                 method='Nelder-Mead',
                                 options={'maxiter': 3000, 'xatol': 1e-3, 'fatol': 1})
                    if res.fun < best_B_chi2:
                        best_B_chi2 = res.fun
                        best_B = res.x
                except: pass

if best_B is not None:
    gd, gi, gs, la0 = best_B
    a0_B = 10**la0
    print(f"\n  Best fit:")
    print(f"    gamma_disk       = {gd:.4f}  (d_eff = {3-2*gd:.2f})")
    print(f"    gamma_intermed   = {gi:.4f}  (d_eff = {3-2*gi:.2f})")
    print(f"    gamma_sphere     = {gs:.4f}  (d_eff = {3-2*gs:.2f})")
    print(f"    a0               = {a0_B:.3e}")
    print(f"    chi2             = {best_B_chi2:.1f}")

    chi2_gal = sum(galaxy_chi2(g, gd, a0_B, Yd_grid) for g in galaxies)
    chi2_ei = sum(elliptical_chi2(e, gi, a0_B) for e in ellipticals if e.gtype == 'fast')
    chi2_es = sum(elliptical_chi2(e, gs, a0_B) for e in ellipticals if e.gtype == 'slow')
    chi2_cl = sum(cluster_chi2(c, gs, a0_B) for c in clusters)
    print(f"\n  Decomposed chi2:")
    print(f"    Disk galaxies:       {chi2_gal:.0f}")
    print(f"    Fast/disky E:        {chi2_ei:.0f}")
    print(f"    Slow E + clusters:   {chi2_es + chi2_cl:.0f}")
    print(f"    Total:               {chi2_gal + chi2_ei + chi2_es + chi2_cl:.0f}")

sys.stdout.flush()

# ============================================================
# SECTION 3: MODEL C — gamma(M) + geometry
# gamma = gamma_0(morph) + delta * log10(M/M_ref)
# ============================================================
print("\n" + "=" * 80)
print("SECTION 3: MODEL C — gamma(morphology, mass)")
print("gamma = gamma_0(morph) + delta * log10(M/M_ref)")
print("=" * 80)
sys.stdout.flush()

def model_C_total_chi2(params):
    gamma_d0, gamma_s0, delta, log_a0 = params
    a0 = 10**log_a0
    if gamma_d0 < 0.05 or gamma_d0 > 1.0: return 1e15
    if gamma_s0 < 0.05 or gamma_s0 > 1.5: return 1e15
    if abs(delta) > 0.3: return 1e15
    M_ref = 1e10

    total = 0
    for g in galaxies:
        M_est = g.estimate_mass(0.5)
        gamma = gamma_d0 + delta * np.log10(max(M_est, 1e6) / M_ref)
        gamma = np.clip(gamma, 0.05, 1.5)
        total += galaxy_chi2(g, gamma, a0, Yd_grid)

    for e in ellipticals:
        gamma = gamma_s0 + delta * np.log10(max(e.M_star, 1e6) / M_ref)
        gamma = np.clip(gamma, 0.05, 1.5)
        total += elliptical_chi2(e, gamma, a0)

    for c in clusters:
        gamma = gamma_s0 + delta * np.log10(max(c.M_total_Msun, 1e6) / M_ref)
        gamma = np.clip(gamma, 0.05, 1.5)
        total += cluster_chi2(c, gamma, a0)

    return total

best_C = None
best_C_chi2 = 1e30
for gd in [0.30, 0.36, 0.40]:
    for gs in [0.40, 0.50, 0.60]:
        for delta in [-0.02, 0.0, 0.02, 0.04]:
            for la0 in [-10.0, -9.95, -9.85]:
                try:
                    res = minimize(model_C_total_chi2, [gd, gs, delta, la0],
                                 method='Nelder-Mead',
                                 options={'maxiter': 3000, 'xatol': 1e-3, 'fatol': 1})
                    if res.fun < best_C_chi2:
                        best_C_chi2 = res.fun
                        best_C = res.x
                except: pass

if best_C is not None:
    gd0, gs0, delta, la0 = best_C
    a0_C = 10**la0
    print(f"\n  Best fit:")
    print(f"    gamma_disk_0     = {gd0:.4f}")
    print(f"    gamma_sphere_0   = {gs0:.4f}")
    print(f"    delta (per dex)  = {delta:.4f}")
    print(f"    a0               = {a0_C:.3e}")
    print(f"    chi2             = {best_C_chi2:.1f}")

    print(f"\n  Implied gamma for representative systems:")
    for label, morph, logM in [
        ('Dwarf (10^8)', 'disk', 8), ('MW-like (5e10)', 'disk', 10.7),
        ('NGC2841 (3e11)', 'disk', 11.5),
        ('NGC3379 fast E (10^11)', 'sphere', 11),
        ('NGC4486 slow E (6e11)', 'sphere', 11.8),
        ('Fornax group (3e13)', 'sphere', 13.5),
        ('Coma cluster (7e14)', 'sphere', 14.85)]:
        g0 = gd0 if morph == 'disk' else gs0
        gamma = g0 + delta * (logM - 10)
        gamma = np.clip(gamma, 0.05, 1.5)
        print(f"    {label:<30} gamma={gamma:.3f}  d_eff={3-2*gamma:.2f}")

sys.stdout.flush()

# ============================================================
# SECTION 4: MODEL D — d_eff as function of geometry
# d_eff = 3 - f(S), where f(S) maps geometry to dimensional reduction
# ============================================================
print("\n" + "=" * 80)
print("SECTION 4: EFFECTIVE DIMENSION MODEL")
print("d_eff = 3 - 2*gamma, gamma = gamma_base * (1 + k*S)")
print("=" * 80)
sys.stdout.flush()

def model_D_total_chi2(params):
    gamma_base, k, log_a0 = params
    a0 = 10**log_a0
    if gamma_base < 0.1 or gamma_base > 0.8: return 1e15
    if k < 0 or k > 5: return 1e15

    total = 0
    for g in galaxies:
        gamma = gamma_base * (1 + k * 0.0)  # S=0 for disks
        total += galaxy_chi2(g, gamma, a0, Yd_grid)

    for e in ellipticals:
        gamma = gamma_base * (1 + k * e.sphericity)
        gamma = min(gamma, 1.4)
        total += elliptical_chi2(e, gamma, a0)

    for c in clusters:
        gamma = gamma_base * (1 + k * 1.0)  # S=1 for clusters
        gamma = min(gamma, 1.4)
        total += cluster_chi2(c, gamma, a0)

    return total

best_D = None
best_D_chi2 = 1e30
for gb in [0.30, 0.36, 0.40]:
    for k in [0.3, 0.5, 0.7, 1.0, 1.5]:
        for la0 in [-10.0, -9.95, -9.85]:
            try:
                res = minimize(model_D_total_chi2, [gb, k, la0],
                             method='Nelder-Mead',
                             options={'maxiter': 3000, 'xatol': 1e-3, 'fatol': 1})
                if res.fun < best_D_chi2:
                    best_D_chi2 = res.fun
                    best_D = res.x
            except: pass

if best_D is not None:
    gb, k, la0 = best_D
    a0_D = 10**la0
    print(f"\n  Best fit:")
    print(f"    gamma_base = {gb:.4f}")
    print(f"    k (geometry factor) = {k:.4f}")
    print(f"    a0 = {a0_D:.3e}")
    print(f"    chi2 = {best_D_chi2:.1f}")

    print(f"\n  Implied dimensions:")
    for label, S in [('Disk (S=0)', 0.0), ('Fast E (S=0.55)', 0.55),
                     ('Slow E (S=0.85)', 0.85), ('Cluster (S=1)', 1.0)]:
        gamma = gb * (1 + k * S)
        d_eff = 3 - 2*gamma
        print(f"    {label:<20} gamma={gamma:.3f}  d_eff={d_eff:.2f}")

sys.stdout.flush()

# ============================================================
# SECTION 5: PER-CLUSTER PREDICTIONS (best model)
# ============================================================
print("\n" + "=" * 80)
print("SECTION 5: CLUSTER PREDICTIONS WITH BEST MODEL")
print("=" * 80)
sys.stdout.flush()

# Use Model A parameters (simplest with good fit)
if best_A is not None:
    gd, gs, la0 = best_A
    a0 = 10**la0
    print(f"\n  Using Model A: gamma_disk={gd:.3f}, gamma_sphere={gs:.3f}, a0={a0:.3e}")
    print(f"\n  {'Cluster':<12} {'g_obs/g_pred':>12} {'status':<8}")
    print(f"  {'-'*35}")

    n_ok, n_close, n_fail = 0, 0, 0
    for c in clusters:
        y = c.g_bar / a0
        gpred = c.g_bar * nu_exp(y, gs)
        ratio = c.g_obs / gpred
        if 0.8 <= ratio <= 1.2:
            status = 'OK'
            n_ok += 1
        elif 0.7 <= ratio <= 1.3:
            status = 'close'
            n_close += 1
        else:
            status = 'FAIL'
            n_fail += 1
        print(f"  {c.name:<12} {ratio:12.3f} {status:<8}")

    print(f"\n  OK (0.8-1.2): {n_ok}/29")
    print(f"  Close (0.7-1.3): {n_ok+n_close}/29")
    print(f"  Fail: {n_fail}/29")

sys.stdout.flush()

# ============================================================
# SECTION 6: PER-ELLIPTICAL PREDICTIONS
# ============================================================
print("\n" + "=" * 80)
print("SECTION 6: ELLIPTICAL PREDICTIONS WITH BEST MODEL")
print("=" * 80)
sys.stdout.flush()

if best_A is not None:
    gd, gs, la0 = best_A
    a0 = 10**la0

    # Group by galaxy
    from collections import defaultdict
    ell_by_gal = defaultdict(list)
    for e in ellipticals:
        ell_by_gal[e.name].append(e)

    print(f"\n  {'Galaxy':<12} {'type':<6} {'N':>3} {'med ratio':>10} {'scatter':>8}")
    print(f"  {'-'*45}")

    for gname in sorted(ell_by_gal.keys()):
        pts = ell_by_gal[gname]
        gtype = pts[0].gtype
        S = pts[0].sphericity
        gamma = gd + (gs - gd) * S

        ratios = []
        for e in pts:
            y = e.g_bar / a0
            gpred = e.g_bar * nu_exp(y, gamma)
            ratios.append(e.g_obs / gpred)

        med = np.median(ratios)
        scat = np.std(ratios) if len(ratios) > 1 else 0
        print(f"  {gname:<12} {gtype:<6} {len(pts):3d} {med:10.3f} {scat:8.3f}")

sys.stdout.flush()

# ============================================================
# SECTION 7: MODEL COMPARISON
# ============================================================
print("\n" + "=" * 80)
print("SECTION 7: MODEL COMPARISON")
print("=" * 80)
sys.stdout.flush()

N_data = len(galaxies) * 20 + len(ellipticals) + len(clusters)  # approx
ln_N = np.log(N_data)

models = []

# Baselines
for gamma_test in [0.36, 0.40, 0.50, 0.54]:
    chi2 = sum(galaxy_chi2(g, gamma_test, a0_ref, Yd_grid) for g in galaxies)
    chi2 += sum(elliptical_chi2(e, gamma_test, a0_ref) for e in ellipticals)
    chi2 += sum(cluster_chi2(c, gamma_test, a0_ref) for c in clusters)
    models.append((f'Baseline g={gamma_test}', chi2, 1, 0))

if best_A is not None:
    models.append(('A: gamma(S) linear', best_A_chi2, 3, 0))
if best_B is not None:
    models.append(('B: 3-class', best_B_chi2, 4, 0))
if best_C is not None:
    models.append(('C: gamma(morph,M)', best_C_chi2, 4, 0))
if best_D is not None:
    models.append(('D: d_eff(S)', best_D_chi2, 3, 0))

# Sort by chi2
models.sort(key=lambda x: x[1])
ref_chi2 = models[-1][1]  # worst baseline

print(f"\n  {'Model':<25} {'chi2':>10} {'n_par':>6} {'dchi2':>10} {'BIC pen':>8} {'net dBIC':>10}")
print(f"  {'-'*75}")
for name, chi2, npar, _ in models:
    dchi2 = ref_chi2 - chi2
    bic_pen = (npar - 1) * ln_N
    net = dchi2 - bic_pen
    print(f"  {name:<25} {chi2:10.0f} {npar:6d} {dchi2:+10.0f} {bic_pen:8.1f} {net:+10.0f}")

sys.stdout.flush()

# ============================================================
# SECTION 8: PHYSICAL INTERPRETATION
# ============================================================
print("\n" + "=" * 80)
print("SECTION 8: PHYSICAL INTERPRETATION FOR TGP")
print("=" * 80)

if best_A is not None:
    gd, gs, la0 = best_A
    print(f"""
  GEOMETRY-DEPENDENT DIMENSIONAL REDUCTION:

  gamma_disk   = {gd:.3f}  =>  d_eff = {3-2*gd:.2f}
  gamma_sphere = {gs:.3f}  =>  d_eff = {3-2*gs:.2f}

  Physical mechanism:
  1. In DISK geometry (galaxies):
     - Matter is already distributed in ~2D (thin disk)
     - Substrate deformation is "pre-configured" for 2D propagation
     - Less dimensional reduction needed: gamma ~ {gd:.2f}
     - Effective dimension d ~ {3-2*gd:.1f} (partial 3D->2D)

  2. In SPHERICAL geometry (ellipticals + clusters):
     - Matter fills full 3D volume
     - Substrate must transition from 3D to 2D everywhere
     - MORE dimensional reduction needed: gamma ~ {gs:.2f}
     - Effective dimension d ~ {3-2*gs:.1f} (deeper transition)

  3. INTERMEDIATE geometry (S0 / fast-rotator E):
     - Partially disky: S ~ 0.55
     - gamma ~ {gd + (gs-gd)*0.55:.2f} (interpolated)
     - d_eff ~ {3-2*(gd+(gs-gd)*0.55):.2f}

  HIERARCHY:
    disk (d={3-2*gd:.2f}) -> S0 (d~{3-2*(gd+(gs-gd)*0.3):.2f}) -> fast E (d~{3-2*(gd+(gs-gd)*0.55):.2f})
    -> slow E (d~{3-2*(gd+(gs-gd)*0.85):.2f}) -> cluster (d={3-2*gs:.2f})

  The DEEPER the 3D->2D transition, the stronger the MOND-like boost.

  TESTABLE PREDICTIONS:
  1. S0 galaxies should have gamma ~ {gd+(gs-gd)*0.3:.2f} (between disk and E)
  2. BCGs (brightest cluster galaxies): gamma close to cluster value
  3. Edge-on vs face-on spirals: SAME gamma (it's 3D geometry, not projection)
  4. Dwarf spheroidals: should prefer higher gamma if truly spheroidal
  5. Polar ring galaxies: test of geometry vs mass interpretation
""")

sys.stdout.flush()

# ============================================================
# SECTION 9: DWARF SPHEROIDAL PREDICTION
# ============================================================
print("=" * 80)
print("SECTION 9: SPARC BULGE-DOMINATED GALAXIES — do they prefer higher gamma?")
print("=" * 80)
sys.stdout.flush()

# Check if SPARC galaxies with significant bulge (proxy for spheroidal component)
# show different gamma preference
bulge_gals = [g for g in all_galaxies if g.has_bulge]
disk_gals = [g for g in all_galaxies if not g.has_bulge]

print(f"\n  SPARC bulge-dominated: {len(bulge_gals)} galaxies")
print(f"  SPARC pure-disk:      {len(disk_gals)} galaxies")

# Sample for speed
bulge_sample = bulge_gals[::2]
disk_sample = disk_gals[::4]

for label, sample in [('Bulge-dominated', bulge_sample), ('Pure disk', disk_sample)]:
    chi2_036 = sum(galaxy_chi2(g, 0.36, a0_ref, Yd_grid) for g in sample)
    chi2_040 = sum(galaxy_chi2(g, 0.40, a0_ref, Yd_grid) for g in sample)
    chi2_050 = sum(galaxy_chi2(g, 0.50, a0_ref, Yd_grid) for g in sample)
    chi2_054 = sum(galaxy_chi2(g, 0.54, a0_ref, Yd_grid) for g in sample)

    best = min([(0.36, chi2_036), (0.40, chi2_040),
                (0.50, chi2_050), (0.54, chi2_054)],
               key=lambda x: x[1])

    print(f"\n  {label} (N={len(sample)}):")
    print(f"    chi2(0.36)={chi2_036:.0f}  chi2(0.40)={chi2_040:.0f}  chi2(0.50)={chi2_050:.0f}  chi2(0.54)={chi2_054:.0f}")
    print(f"    Best: gamma={best[0]:.2f}")

sys.stdout.flush()

print("\n" + "=" * 80)
print("DONE")
print("=" * 80)
sys.stdout.flush()
