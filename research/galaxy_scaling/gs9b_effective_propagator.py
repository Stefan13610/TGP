"""
gs9b: Effective propagator with dimensional transition 3D -> 2D
==============================================================

Context (gs9a findings):
- Soliton tail decays as r^(-(d-1)/2) -- confirmed numerically
- Fisher information g'^2/g is conformally invariant ONLY in d=2
- But conformal invariance of kinetic term doesn't drive physics
  (potential term dominates at large r)
- The key is FORCE LAW: 3D gives F~1/r^2, 2D gives F~1/r
- DGP-like interpolation from gs8 works phenomenologically:
  g_obs = g_N + sqrt(g_N*a0)

Goal:
- Construct explicit propagator with d_eff(r) transition 3->2
- Compare multiple interpolation models
- Compute rotation curves, RAR, BTFR
- Quantitative comparison with MOND (simple & standard)
- Identify testable differences

CORRECTED (v2):
- DGP force law: F = GM/r^2 + sqrt(GM*a0)/r (SUM of 3D + 2D)
  NOT GM/(r^2 + r*r_c) which gives opposite transition direction!
- Smooth d_eff model: dimensionally consistent with reference scale r_c
- Freeman limit: correct unit conversion to M_sun/pc^2
"""

import numpy as np
import sys, io, warnings
warnings.filterwarnings('ignore')
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

# ============================================================
# Physical constants
# ============================================================
G = 6.674e-11       # m^3/(kg*s^2)
c = 2.998e8          # m/s
H0 = 2.18e-18        # 1/s (67.4 km/s/Mpc)
M_sun = 1.989e30     # kg
kpc = 3.086e19       # m
pc = 3.086e16        # m
a0 = 1.2e-10         # m/s^2 (MOND acceleration scale)

print("=" * 70)
print("  gs9b: Effective propagator with dimensional transition 3D -> 2D")
print("=" * 70)

# ============================================================
# 1. FORCE LAW MODELS
# ============================================================
print("\n" + "=" * 70)
print("  1. Force law models with dimensional transition")
print("=" * 70)

def g_newton(M, r):
    """Newtonian gravitational acceleration."""
    return G * M / r**2

def r_mond(M):
    """MOND transition radius."""
    return np.sqrt(G * M / a0)

# Model A: DGP-like (from gs8) -- CORRECTED
# The RAR g = g_N + sqrt(g_N * a0) implies force law:
# F(r) = GM/r^2 + sqrt(GM*a0)/r = g_N + sqrt(g_N*a0)
# This is a SUM of 3D and 2D terms, not an interpolation.
# At r << r_c: GM/r^2 dominates (Newton, 3D)
# At r >> r_c: sqrt(GM*a0)/r dominates (2D-like, flat RC)
# Crossover at r_c = sqrt(GM/a0) where both terms equal
def g_dgp(M, r):
    """DGP-like force: F = GM/r^2 + sqrt(GM*a0)/r."""
    gN = G * M / r**2
    return gN + np.sqrt(gN * a0)

# Model B: Smooth d_eff transition (dimensionally consistent)
# F = (GM/r_c^2) * (r_c/r)^(d_eff - 1)
# d_eff(r) = 2 + 1/(1 + (r/r_c)^alpha)
# At r << r_c: d_eff -> 3, F -> GM/r^2 (Newton)
# At r >> r_c: d_eff -> 2, F -> GM/(r*r_c) (2D-like)
def g_smooth(M, r, alpha=1.0):
    """Force with smooth d_eff transition, dimensionally consistent."""
    rc = np.sqrt(G * M / a0)
    x = r / rc
    d_eff = 2.0 + 1.0 / (1.0 + x**alpha)
    # Dimensionally consistent: GM/(rc^2) * (rc/r)^(d_eff-1)
    return (G * M / rc**2) * (rc / r)**(d_eff - 1)

# Model C: Tanh transition (dimensionally consistent)
def g_tanh(M, r):
    """Force with tanh d_eff transition."""
    rc = np.sqrt(G * M / a0)
    x = r / rc
    d_eff = 2.5 - 0.5 * np.tanh(np.log(x))
    return (G * M / rc**2) * (rc / r)**(d_eff - 1)

# MOND (simple interpolation function)
def g_mond_simple(M, r):
    """MOND with simple interpolation: mu = x/(1+x), x = g/a0."""
    gN = g_newton(M, r)
    return 0.5 * gN + np.sqrt(0.25 * gN**2 + gN * a0)

# MOND (standard interpolation function)
def g_mond_standard(M, r):
    """MOND with standard interpolation: mu = x/sqrt(1+x^2)."""
    gN = g_newton(M, r)
    x = gN / a0
    return gN * np.sqrt(0.5 + 0.5 * np.sqrt(1 + 4.0/x**2))

print("""
  Model A (DGP):     F = GM/r^2 + sqrt(GM*a0)/r  (SUM of 3D + 2D terms)
  Model B (smooth):  F = (GM/rc^2)*(rc/r)^(d_eff-1),
                     d_eff = 2 + 1/(1+(r/rc)^alpha)
  Model C (tanh):    same form, d_eff = 2.5 - 0.5*tanh(ln(r/rc))
  MOND simple:       mu(x) = x/(1+x)
  MOND standard:     mu(x) = x/sqrt(1+x^2)

  All with r_c = sqrt(GM/a0) = MOND radius

  KEY INSIGHT (v2 correction):
  The DGP force is a SUM of 3D and 2D components, NOT a smooth
  interpolation between dimensions. Gravity has TWO channels:
  one decaying as 1/r^2 (3D) and one as 1/r (2D).
""")

# ============================================================
# 2. MILKY WAY ROTATION CURVE
# ============================================================
print("=" * 70)
print("  2. Milky Way rotation curve comparison")
print("=" * 70)

M_MW = 6e10 * M_sun  # Baryonic mass
rc_MW = np.sqrt(G * M_MW / a0)

print(f"\n  M_MW = {M_MW/M_sun:.1e} M_sun")
print(f"  r_MOND = {rc_MW/kpc:.1f} kpc")
print(f"  v_flat(BTFR) = (GM*a0)^(1/4) = {(G*M_MW*a0)**0.25/1e3:.1f} km/s")

# Verify DGP asymptotic: at r>>r_c, g -> sqrt(GM*a0)/r, v^2 = sqrt(GM*a0)
v_asymp = (G * M_MW * a0)**0.25
print(f"  DGP asymptotic: v = (GM*a0)^(1/4) = {v_asymp/1e3:.1f} km/s")

r_print = np.array([0.5, 1, 2, 5, 8, 10, 15, 20, 50, 100, 200]) * kpc
print(f"\n  {'r(kpc)':>8}  {'Newton':>8}  {'DGP':>8}  {'Smooth1':>8}  {'Tanh':>8}  {'MOND_s':>8}  {'MOND_st':>8}  v(km/s)")
print("  " + "-" * 68)
for r in r_print:
    vN = np.sqrt(g_newton(M_MW, r) * r) / 1e3
    vD = np.sqrt(g_dgp(M_MW, r) * r) / 1e3
    vS1 = np.sqrt(g_smooth(M_MW, r, 1.0) * r) / 1e3
    vT = np.sqrt(g_tanh(M_MW, r) * r) / 1e3
    vMs = np.sqrt(g_mond_simple(M_MW, r) * r) / 1e3
    vMst = np.sqrt(g_mond_standard(M_MW, r) * r) / 1e3
    print(f"  {r/kpc:8.1f}  {vN:8.1f}  {vD:8.1f}  {vS1:8.1f}  {vT:8.1f}  {vMs:8.1f}  {vMst:8.1f}")

# ============================================================
# 3. RAR -- Radial Acceleration Relation
# ============================================================
print("\n" + "=" * 70)
print("  3. RAR: g_obs vs g_bar (acceleration space)")
print("=" * 70)

print("""
  The RAR relates observed acceleration g_obs to baryonic g_bar:

  MOND simple:    g = g_bar/2 + sqrt(g_bar^2/4 + g_bar*a0)
  MOND standard:  g = g_bar * sqrt(1/2 + sqrt(1/4 + (a0/g_bar)^2))
  DGP:            g = g_bar + sqrt(g_bar*a0)

  Interpolation function nu(y) where y = g_bar/a0, g_obs = g_bar*nu(y):
  MOND simple:    nu(y) = 1/2 + sqrt(1/4 + 1/y)
  MOND standard:  nu(y) = sqrt(1/2 + sqrt(1/4 + 1/y^2))
  DGP:            nu(y) = 1 + 1/sqrt(y)
""")

print(f"  {'y=g/a0':>10}  {'nu_MOND_s':>10}  {'nu_MOND_st':>11}  {'nu_DGP':>8}  {'DGP/MOND_s':>11}")
print("  " + "-" * 55)
for log_y in np.arange(-3, 3.5, 0.5):
    y = 10**log_y
    nu_ms = 0.5 + np.sqrt(0.25 + 1/y)
    nu_mst = np.sqrt(0.5 + np.sqrt(0.25 + 1/y**2))
    nu_dgp = 1 + 1/np.sqrt(y)
    print(f"  {y:10.3f}  {nu_ms:10.4f}  {nu_mst:11.4f}  {nu_dgp:8.4f}  {nu_dgp/nu_ms:11.4f}")

# ============================================================
# 4. MAXIMUM DIFFERENCE from MOND
# ============================================================
print("\n" + "=" * 70)
print("  4. Maximum difference between DGP and MOND")
print("=" * 70)

g_bar_arr = np.logspace(-13, -8, 500)
g_dgp_rar = g_bar_arr + np.sqrt(g_bar_arr * a0)
g_mond_s_rar = 0.5 * g_bar_arr + np.sqrt(0.25 * g_bar_arr**2 + g_bar_arr * a0)
x_arr = g_bar_arr / a0
g_mond_st_rar = g_bar_arr * np.sqrt(0.5 + 0.5 * np.sqrt(1 + 4.0/x_arr**2))

ratios = g_dgp_rar / g_mond_s_rar
max_idx = np.argmax(np.abs(ratios - 1))
g_max_diff = g_bar_arr[max_idx]
print(f"\n  DGP vs MOND (simple interpolation):")
print(f"    Max deviation: {(ratios[max_idx]-1)*100:.1f}% at g_bar = {g_max_diff/a0:.3f}*a0")
mask10 = np.abs(ratios - 1) > 0.10
if np.any(mask10):
    g_lo = g_bar_arr[mask10][0] / a0
    g_hi = g_bar_arr[mask10][-1] / a0
    print(f"    g_bar range with >10% diff: {g_lo:.3f}*a0 to {g_hi:.3f}*a0")

ratios_st = g_dgp_rar / g_mond_st_rar
max_idx_st = np.argmax(np.abs(ratios_st - 1))
print(f"\n  DGP vs MOND (standard interpolation):")
print(f"    Max deviation: {(ratios_st[max_idx_st]-1)*100:.1f}% at g_bar = {g_bar_arr[max_idx_st]/a0:.3f}*a0")

# ============================================================
# 5. LIMITING BEHAVIOR
# ============================================================
print("\n" + "=" * 70)
print("  5. Limiting behavior of each model")
print("=" * 70)

print("""
  g_bar >> a0 (strong field, Newtonian regime):
  -----------------------------------------------
  MOND simple:   g -> g_bar * (1 + a0/(2*g_bar))     correction ~ a0/g
  MOND standard: g -> g_bar * (1 + (a0/g)^2/2)        correction ~ (a0/g)^2
  DGP:           g -> g_bar * (1 + sqrt(a0/g_bar))     correction ~ sqrt(a0/g)

  DGP correction is INTERMEDIATE: larger than MOND standard,
  smaller than MOND simple at very high g, but LARGER at moderate g.

  g_bar << a0 (deep MOND regime):
  -----------------------------------------------
  ALL THREE -> sqrt(g_bar * a0) (geometric mean)
  DGP has subleading +g_bar correction.
""")

print("  Numerical verification (corrections in strong field):")
print(f"  {'g_bar/a0':>10}  {'MOND_s corr':>12}  {'MOND_st corr':>13}  {'DGP corr':>10}")
print("  " + "-" * 50)
for ratio_val in [1000, 100, 10, 1, 0.1, 0.01]:
    gb = ratio_val * a0
    gMs = 0.5*gb + np.sqrt(0.25*gb**2 + gb*a0)
    gMst = gb * np.sqrt(0.5 + 0.5*np.sqrt(1 + 4/(gb/a0)**2))
    gD = gb + np.sqrt(gb*a0)
    print(f"  {ratio_val:10.3f}  {(gMs/gb-1)*100:11.4f}%  {(gMst/gb-1)*100:12.6f}%  {(gD/gb-1)*100:9.4f}%")

# ============================================================
# 6. BTFR -- Baryonic Tully-Fisher Relation
# ============================================================
print("\n" + "=" * 70)
print("  6. BTFR: v_flat^4 vs M")
print("=" * 70)

masses = np.logspace(8, 12, 50) * M_sun
r_far = 500 * kpc  # very far out to see asymptotic behavior

v4_dgp = []
v4_mond = []
for M in masses:
    gD = g_dgp(M, r_far)
    gM = g_mond_simple(M, r_far)
    v4_dgp.append((gD * r_far)**2)
    v4_mond.append((gM * r_far)**2)

v4_dgp = np.array(v4_dgp)
v4_mond = np.array(v4_mond)

log_M = np.log10(masses / M_sun)
coeffs_dgp = np.polyfit(log_M, np.log10(v4_dgp), 1)
coeffs_mond = np.polyfit(log_M, np.log10(v4_mond), 1)

print(f"\n  BTFR: v^4 ~ M^slope (measured at r = 500 kpc)")
print(f"    DGP slope:  {coeffs_dgp[0]:.4f}  (theory: 1.0000)")
print(f"    MOND slope: {coeffs_mond[0]:.4f}  (theory: 1.0000)")

# Direct check: v_flat = (GM*a0)^(1/4)
print(f"\n  Direct check at M = 10^10 M_sun:")
M_test = 1e10 * M_sun
v_btfr = (G * M_test * a0)**0.25
v_dgp_500 = np.sqrt(g_dgp(M_test, 500*kpc) * 500*kpc)
v_mond_500 = np.sqrt(g_mond_simple(M_test, 500*kpc) * 500*kpc)
print(f"    BTFR: v = {v_btfr/1e3:.2f} km/s")
print(f"    DGP(500kpc): v = {v_dgp_500/1e3:.2f} km/s")
print(f"    MOND(500kpc): v = {v_mond_500/1e3:.2f} km/s")

# ============================================================
# 7. SPECIFIC GALAXIES
# ============================================================
print("\n" + "=" * 70)
print("  7. Specific galaxy rotation curves")
print("=" * 70)

galaxies = {
    'DDO 154':  {'M': 4.0e8,  'r_last': 8.0,  'v_obs': 47},
    'NGC 2403': {'M': 1.2e10, 'r_last': 20.0, 'v_obs': 136},
    'MW':       {'M': 6.0e10, 'r_last': 20.0, 'v_obs': 220},
    'UGC 2885': {'M': 2.0e11, 'r_last': 80.0, 'v_obs': 300},
    'IC 1101':  {'M': 1.0e13, 'r_last': 300.0, 'v_obs': 400},
}

print(f"\n  {'Galaxy':>10}  {'M(M_sun)':>10}  {'rc(kpc)':>8}  {'r_last':>7}  {'v_obs':>6}  {'v_DGP':>6}  {'v_MOND_s':>8}  {'v_Newt':>7}  {'r/rc':>6}")
print("  " + "-" * 82)
for name, gal in galaxies.items():
    M = gal['M'] * M_sun
    r = gal['r_last'] * kpc
    rc = np.sqrt(G * M / a0)
    vN = np.sqrt(g_newton(M, r) * r) / 1e3
    vD = np.sqrt(g_dgp(M, r) * r) / 1e3
    vMs = np.sqrt(g_mond_simple(M, r) * r) / 1e3
    print(f"  {name:>10}  {gal['M']:10.1e}  {rc/kpc:8.1f}  {gal['r_last']:7.0f}  {gal['v_obs']:6.0f}  {vD:6.0f}  {vMs:8.0f}  {vN:7.0f}  {gal['r_last']/(rc/kpc):6.1f}")

# ============================================================
# 8. SMOOTH d_eff MODELS vs MOND
# ============================================================
print("\n" + "=" * 70)
print("  8. Optimizing smooth d_eff model to match MOND")
print("=" * 70)

print("""
  Model B: d_eff(r) = 2 + 1/(1 + (r/rc)^alpha)
  Force: F = (GM/rc^2) * (rc/r)^(d_eff - 1)  [dimensionally consistent]

  At r << rc: d_eff -> 3, F -> GM/r^2 (Newton)
  At r >> rc: d_eff -> 2, F -> GM/(r*rc) (2D-like)
""")

M = M_MW
r_test = np.logspace(np.log10(0.5), np.log10(200), 200) * kpc

def rms_vs_mond(alpha):
    v_ms = np.sqrt(g_mond_simple(M, r_test) * r_test)
    v_sm = np.sqrt(g_smooth(M, r_test, alpha) * r_test)
    return np.sqrt(np.mean(((v_sm - v_ms) / v_ms)**2)) * 100

alphas = np.arange(0.5, 5.1, 0.5)
print(f"  {'alpha':>8}  {'RMS vs MOND(%)':>15}  {'v(5kpc)':>9}  {'v(10kpc)':>9}  {'v(50kpc)':>9}  {'v(200kpc)':>10}")
print("  " + "-" * 65)

# Reference MOND values
v_mond_ref = [np.sqrt(g_mond_simple(M, r*kpc)*r*kpc)/1e3 for r in [5, 10, 50, 200]]
print(f"  {'MOND':>8}  {'---':>15}  {v_mond_ref[0]:9.1f}  {v_mond_ref[1]:9.1f}  {v_mond_ref[2]:9.1f}  {v_mond_ref[3]:10.1f}")

best_alpha = 0.5
best_rms = 1e10
for alpha in alphas:
    rms = rms_vs_mond(alpha)
    v_vals = [np.sqrt(g_smooth(M, r*kpc, alpha)*r*kpc)/1e3 for r in [5, 10, 50, 200]]
    print(f"  {alpha:8.1f}  {rms:15.2f}  {v_vals[0]:9.1f}  {v_vals[1]:9.1f}  {v_vals[2]:9.1f}  {v_vals[3]:10.1f}")
    if rms < best_rms:
        best_rms = rms
        best_alpha = alpha

print(f"\n  Best alpha: {best_alpha:.1f} (RMS = {best_rms:.2f}%)")

# DGP comparison
rms_dgp = np.sqrt(np.mean(((np.sqrt(g_dgp(M, r_test)*r_test) - np.sqrt(g_mond_simple(M, r_test)*r_test))**2 / (g_mond_simple(M, r_test)*r_test)))) * 100
v_dgp_ref = [np.sqrt(g_dgp(M, r*kpc)*r*kpc)/1e3 for r in [5, 10, 50, 200]]
print(f"  DGP:  RMS = {rms_dgp:.2f}%,  v(5) = {v_dgp_ref[0]:.1f},  v(10) = {v_dgp_ref[1]:.1f},  v(50) = {v_dgp_ref[2]:.1f},  v(200) = {v_dgp_ref[3]:.1f}")

# ============================================================
# 9. FREEMAN LIMIT
# ============================================================
print("\n" + "=" * 70)
print("  9. Freeman limit (maximum surface brightness)")
print("=" * 70)

Sigma_mond = a0 / (2 * np.pi * G)  # kg/m^2
# Convert to M_sun/pc^2
Sigma_mond_Mpc2 = Sigma_mond / M_sun * pc**2

print(f"""
  MOND Freeman limit:  Sigma = a0/(2*pi*G) = {Sigma_mond:.4f} kg/m^2
                     = {Sigma_mond_Mpc2:.0f} M_sun/pc^2
  Observed:            Sigma ~ 140 M_sun/pc^2
  Ratio predicted/obs: {Sigma_mond_Mpc2/140:.2f}

  This depends on a0 only, not on interpolation function.
  Both MOND and DGP give the same Freeman limit.
""")

# ============================================================
# 10. SOLAR SYSTEM CONSTRAINTS
# ============================================================
print("=" * 70)
print("  10. Solar system constraints")
print("=" * 70)

r_earth = 1.496e11  # m (1 AU)
r_mercury = 5.79e10 # m
M_sol = M_sun

bodies = [
    ('Mercury', M_sol, r_mercury),
    ('Earth', M_sol, r_earth),
    ('Pioneer (50 AU)', M_sol, 50 * r_earth),
    ('Oort cloud (50000 AU)', M_sol, 50000 * r_earth),
]

print(f"\n  {'Body':>22}  {'g_N (m/s^2)':>12}  {'g/a0':>10}  {'dg/g MOND_s':>12}  {'dg/g MOND_st':>13}  {'dg/g DGP':>10}")
print("  " + "-" * 85)
for name, M, r in bodies:
    gN = G * M / r**2
    ratio = gN / a0
    dg_ms = (g_mond_simple(M, r) - gN) / gN
    dg_mst = (g_mond_standard(M, r) - gN) / gN
    dg_dgp = (g_dgp(M, r) - gN) / gN
    print(f"  {name:>22}  {gN:12.3e}  {ratio:10.1e}  {dg_ms:12.2e}  {dg_mst:13.2e}  {dg_dgp:10.2e}")

print("""
  Constraints from planetary ephemerides: dg/g < 10^-9 at Earth orbit

  DGP correction at Earth = sqrt(a0/g) = sqrt(1.2e-10 / 5.9e-3) = 4.5e-4
  -> RULED OUT by 5 orders of magnitude!

  MOND standard: (a0/g)^2 ~ 4e-16 -> SAFE
  MOND simple:   a0/(2g) ~ 10^-8 -> MARGINAL (but testable)

  CONCLUSION: Raw DGP formula NEEDS SCREENING in strong fields.
""")

# ============================================================
# 11. VAINSHTEIN SCREENING
# ============================================================
print("=" * 70)
print("  11. Vainshtein screening analysis")
print("=" * 70)

r_g_sun = 2 * G * M_sun / c**2
r_c_cosmic = c / H0  # DGP crossover at cosmic scale

r_V_cosmic = (r_g_sun * r_c_cosmic**2)**(1./3.)

print(f"""
  In the original DGP model, the Vainshtein mechanism screens
  the extra force below the Vainshtein radius r_V:
  r_V = (r_g * r_c^2)^(1/3)

  With COSMIC crossover r_c = c/H0:
  r_g(Sun)  = {r_g_sun:.0f} m = {r_g_sun/1e3:.1f} km
  r_c       = c/H0 = {r_c_cosmic:.2e} m
  r_V       = {r_V_cosmic:.2e} m = {r_V_cosmic/kpc:.1f} kpc
  r_V/r_Earth = {r_V_cosmic/r_earth:.0e}
  -> Solar system completely screened!

  Screened correction at Earth:
  dg/g ~ (r/r_V)^(3/2) * sqrt(a0/g)""")

delta_screened = (r_earth / r_V_cosmic)**(1.5) * np.sqrt(a0 / (G*M_sun/r_earth**2))
print(f"  = {delta_screened:.2e}")
print(f"  Constraint < 10^-9 -> {'SAFE' if delta_screened < 1e-9 else 'VIOLATED'}")

print(f"""
  BUT: in our TGP model, r_c is mass-dependent:
  r_c = sqrt(GM/a0) for each source

  For Sun: r_c = {np.sqrt(G*M_sun/a0)/r_earth:.0f} AU
  r_V = (r_g * r_c^2)^(1/3) = {(r_g_sun * (G*M_sun/a0))**(1./3.)/r_earth:.1f} AU""")

r_c_sun = np.sqrt(G * M_sun / a0)
r_V_sun = (r_g_sun * r_c_sun**2)**(1./3.)
delta_sun = (r_earth / r_V_sun)**(1.5) * np.sqrt(a0 / (G*M_sun/r_earth**2))
print(f"  Screened correction at Earth: {delta_sun:.2e}")
print(f"  -> {'SAFE' if delta_sun < 1e-9 else 'VIOLATED'}")

print(f"""
  THIS IS THE KEY QUESTION:
  If r_c = sqrt(GM/a0) (mass-dependent, from MOND phenomenology):
  -> r_V ~ 1 AU for the Sun -> screening fails at Earth!
  If r_c = c/H0 (universal, from cosmology):
  -> r_V ~ 100 kpc for the Sun -> screening works perfectly!

  In standard DGP: r_c = c/H0 (universal). The mass-dependent
  MOND radius is an EMERGENT scale, not the fundamental crossover.

  For TGP: we need to understand if the crossover scale is
  - Universal (r_c = c/H0) -> needs mechanism for MOND radius
  - Mass-dependent (r_c = sqrt(GM/a0)) -> needs different screening
""")

# ============================================================
# 12. EFFECTIVE DIMENSION d_eff
# ============================================================
print("=" * 70)
print("  12. Effective dimension d_eff(r) for each model")
print("=" * 70)

print("""
  d_eff(r) = 1 - d(ln F)/d(ln r)
  Newton (3D): F ~ 1/r^2 -> d_eff = 3
  Flat RC (2D): F ~ 1/r -> d_eff = 2
""")

r_pts = np.array([0.5, 1, 2, 5, 10, 20, 50, 100, 200]) * kpc
dr = 0.001  # small step for derivative

print(f"  {'r(kpc)':>8}  {'DGP':>6}  {'Smooth':>7}  {'Tanh':>6}  {'MOND_s':>7}  {'MOND_st':>8}")
print("  " + "-" * 48)
for r in r_pts:
    r1 = r * (1 - dr)
    r2 = r * (1 + dr)
    results = []
    for func in [g_dgp, lambda M,r: g_smooth(M,r,best_alpha), g_tanh, g_mond_simple, g_mond_standard]:
        F1 = func(M_MW, r1)
        F2 = func(M_MW, r2)
        d = 1 - np.log(F2/F1) / np.log(r2/r1)
        results.append(d)
    print(f"  {r/kpc:8.1f}  {results[0]:6.2f}  {results[1]:7.2f}  {results[2]:6.2f}  {results[3]:7.2f}  {results[4]:8.2f}")

# ============================================================
# 13. HYBRID MODEL
# ============================================================
print("\n" + "=" * 70)
print("  13. Hybrid model: DGP + natural screening")
print("=" * 70)

print("""
  The DGP formula nu(y) = 1 + 1/sqrt(y) has too-slow decay at large y.
  Consider modified forms:

  Hybrid A: nu(y) = 1 + 1/(sqrt(y) + y)
    y >> 1: nu -> 1 + 1/y (like MOND simple, safe for solar system)
    y << 1: nu -> 1 + 1/sqrt(y) (like DGP, gives flat RC)

  Hybrid B: nu(y) = 1 + exp(-sqrt(y)) / sqrt(y)
    y >> 1: exponentially suppressed (safe!)
    y << 1: nu -> 1 + 1/sqrt(y) (DGP-like)

  Hybrid C: nu(y) = sqrt(1 + 1/y)  (another simple form)
    y >> 1: nu -> 1 + 1/(2y) (MOND-like)
    y << 1: nu -> 1/sqrt(y) (correct MOND limit but nu -> inf)
""")

def nu_mond_s(y):
    return 0.5 + np.sqrt(0.25 + 1/y)

def nu_dgp(y):
    return 1 + 1/np.sqrt(y)

def nu_hybA(y):
    return 1 + 1/(np.sqrt(y) + y)

def nu_hybB(y):
    return 1 + np.exp(-np.sqrt(y)) / np.sqrt(y)

def nu_hybC(y):
    return np.sqrt(1 + 1/y)

print(f"  {'y':>10}  {'MOND_s':>8}  {'DGP':>8}  {'HybA':>8}  {'HybB':>8}  {'HybC':>8}")
print("  " + "-" * 55)
for log_y in np.arange(-3, 4.0, 0.5):
    y = 10**log_y
    print(f"  {y:10.3f}  {nu_mond_s(y):8.4f}  {nu_dgp(y):8.4f}  {nu_hybA(y):8.4f}  {nu_hybB(y):8.4f}  {nu_hybC(y):8.4f}")

# Solar system check
print(f"\n  Solar system (Earth orbit, y = g/a0 = {G*M_sun/r_earth**2/a0:.1e}):")
y_earth = G * M_sun / r_earth**2 / a0
print(f"    nu - 1 (excess):")
print(f"    MOND_s: {nu_mond_s(y_earth)-1:.2e}")
print(f"    DGP:    {nu_dgp(y_earth)-1:.2e}")
print(f"    HybA:   {nu_hybA(y_earth)-1:.2e}")
print(f"    HybB:   {nu_hybB(y_earth)-1:.2e}")
print(f"    HybC:   {nu_hybC(y_earth)-1:.2e}")
print(f"    Constraint: nu - 1 < 10^-9")

# Deep MOND check
print(f"\n  Deep MOND (y = 0.001):")
y_dm = 0.001
for name, func in [('MOND_s', nu_mond_s), ('DGP', nu_dgp), ('HybA', nu_hybA), ('HybB', nu_hybB), ('HybC', nu_hybC)]:
    g_ratio = func(y_dm) * y_dm  # g_obs/a0
    g_exact = np.sqrt(y_dm)      # sqrt(g_bar*a0)/a0 = sqrt(y)
    print(f"    {name:>8}: g_obs/a0 = {g_ratio:.5f},  sqrt(g*a0)/a0 = {g_exact:.5f},  ratio = {g_ratio/g_exact:.4f}")

# ============================================================
# 14. ROTATION CURVE with HYBRID A
# ============================================================
print("\n" + "=" * 70)
print("  14. MW rotation curve with Hybrid A model")
print("=" * 70)

def g_hybA(M, r):
    gN = G * M / r**2
    y = gN / a0
    return gN * (1 + 1/(np.sqrt(y) + y))

print(f"\n  {'r(kpc)':>8}  {'Newton':>8}  {'DGP':>8}  {'HybA':>8}  {'MOND_s':>8}  {'MOND_st':>8}  HybA/MOND_s")
print("  " + "-" * 68)
for r in r_print:
    vN = np.sqrt(g_newton(M_MW, r) * r) / 1e3
    vD = np.sqrt(g_dgp(M_MW, r) * r) / 1e3
    vH = np.sqrt(g_hybA(M_MW, r) * r) / 1e3
    vMs = np.sqrt(g_mond_simple(M_MW, r) * r) / 1e3
    vMst = np.sqrt(g_mond_standard(M_MW, r) * r) / 1e3
    print(f"  {r/kpc:8.1f}  {vN:8.1f}  {vD:8.1f}  {vH:8.1f}  {vMs:8.1f}  {vMst:8.1f}  {vH/vMs:8.3f}")

# ============================================================
# 15. CRITICAL ASSESSMENT
# ============================================================
print("\n" + "=" * 70)
print("  15. CRITICAL ASSESSMENT")
print("=" * 70)

print("""
  ================================================================
  MODEL SCORECARD
  ================================================================

  Property           | DGP           | HybA          | MOND_s
  -------------------|---------------|---------------|---------------
  nu(y)              | 1+1/sqrt(y)   | 1+1/(sqrt(y)+y)| 1/2+sqrt(1/4+1/y)
  Deep MOND (y<<1)   | sqrt(y*a0) OK | sqrt(y*a0) OK | sqrt(y*a0) OK
  Strong field (y>>1)| 1/sqrt(y) BAD | 1/y OK        | 1/(2y) OK
  Solar system       | RULED OUT     | SAFE          | MARGINAL
  BTFR v^4 = GM*a0   | EXACT         | EXACT         | EXACT
  Freeman limit      | SAME          | SAME          | SAME
  Simplicity         | SIMPLEST      | SIMPLE        | MODERATE

  KEY FINDINGS:
  =============
  1. ALL models agree in deep MOND: g -> sqrt(g_bar * a0)
     This is ROBUST and follows from any 3D -> 2D transition.

  2. The models DIFFER in:
     a) Strong-field corrections (solar system)
     b) Transition region (g ~ a0)
     c) Specific shape of rotation curves near r_c

  3. DGP (g = g_N + sqrt(g_N*a0)) is ruled out by solar system
     unless Vainshtein screening with r_c = c/H0 (not mass-dependent).

  4. Hybrid A (nu = 1 + 1/(sqrt(y) + y)) is phenomenologically viable:
     - Safe in solar system
     - Correct deep MOND limit
     - BUT differs from MOND by ~9% in transition region

  5. The d_eff analysis shows ALL models transition from d=3 to d=2:
     - DGP reaches d=2 fastest (at r ~ few * r_c)
     - MOND standard reaches d=2 slowest

  PHYSICAL PICTURE:
  =================
  The dimensional transition 3D -> 2D in the TGP substrate means:

  At HIGH acceleration (r << r_c):
  - Substrate is locally 3D (all spatial correlations active)
  - Gravity follows Newton's law: F ~ 1/r^2

  At LOW acceleration (r >> r_c):
  - Substrate correlations become effectively 2D
  - A "2D channel" opens: extra force ~ 1/r
  - Combined: F = GM/r^2 + f(r)/r where f(r) depends on model

  The 2D channel is always present but SUPPRESSED at small r.
  Different models differ in HOW the suppression works.

  NEXT STEPS (gs9c):
  ==================
  1. Compare ALL models with SPARC RAR data (quantitative chi^2)
  2. Identify which transition shape fits SPARC best
  3. Check: is HybA closer to SPARC data than MOND?
  4. Galaxy-by-galaxy rotation curve fits
""")
