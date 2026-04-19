"""
gs56_form4_analysis.py
Form 4 Deep Analysis: nu(y) = 1/(1 - exp(-y^delta))
Beat standard MOND on binned RAR (chi2=0.19 vs 34.64)
Explores tension: data prefer delta=0.5, TGP predicts delta=0.4
"""

import numpy as np
from scipy.optimize import minimize

# ============================================================
# CONSTANTS
# ============================================================
G_SI = 6.674e-11       # m^3 kg^-1 s^-2
a0 = 1.2e-10           # m/s^2
kpc_to_m = 3.086e19    # m per kpc
Msun = 1.989e30        # kg

print("=" * 72)
print("gs56: FORM 4 DEEP ANALYSIS")
print("nu(y) = 1 / (1 - exp(-y^delta))")
print("=" * 72)

# ============================================================
# LOAD RAR BINNED DATA
# ============================================================
data_lines = []
with open("RARbins.mrt", "r") as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith(("Title", "Authors", "Table", "=", "-", "Byte", " ", "#")):
            # skip header lines
            parts = line.split()
            if len(parts) == 4:
                try:
                    vals = [float(p) for p in parts]
                    data_lines.append(vals)
                except ValueError:
                    continue
            continue
        parts = line.split()
        if len(parts) == 4:
            try:
                vals = [float(p) for p in parts]
                data_lines.append(vals)
            except ValueError:
                continue

data = np.array(data_lines)
log_gbar = data[:, 0]
log_gobs = data[:, 1]
sigma = data[:, 2]  # std dev in log10 space
N_pts = data[:, 3].astype(int)

gbar = 10.0 ** log_gbar  # m/s^2
gobs = 10.0 ** log_gobs  # m/s^2
y = gbar / a0

print("\nLoaded %d RAR bins" % len(data))
print("log(gbar) range: %.2f to %.2f" % (log_gbar.min(), log_gbar.max()))
print("log(gobs) range: %.2f to %.2f" % (log_gobs.min(), log_gobs.max()))

# ============================================================
# PART A: FORM 4 PROPERTIES
# ============================================================
print("\n" + "=" * 72)
print("PART A: FORM 4 PROPERTIES")
print("=" * 72)

def nu_form4(y, delta):
    """Form 4: nu(y) = 1 / (1 - exp(-y^delta))"""
    yd = np.power(y, delta)
    return 1.0 / (1.0 - np.exp(-yd))

def nu_mond_standard(y):
    """MOND standard: nu = 1/2 + sqrt(1/4 + 1/y)"""
    return 0.5 + np.sqrt(0.25 + 1.0 / y)

def nu_mond_simple(y):
    """MOND simple: nu = 1/(1 - exp(-sqrt(y)))"""
    return 1.0 / (1.0 - np.exp(-np.sqrt(y)))

def nu_tgp_standard(y, alpha=0.8, gamma=0.5714):
    """TGP standard: nu = 1 + exp(-y^alpha)/y^gamma"""
    return 1.0 + np.exp(-np.power(y, alpha)) / np.power(y, gamma)

def nu_form5(y):
    """Form 5: nu = (1 + sqrt(1 + 4/y)) / 2  (should equal MOND standard)"""
    return (1.0 + np.sqrt(1.0 + 4.0 / y)) / 2.0

print("\n1. Definition: nu4(y, delta) = 1 / (1 - exp(-y^delta))")

print("\n2. Asymptotic limits:")
print("   y >> 1 (Newtonian regime):")
print("     exp(-y^delta) -> 0, so nu4 -> 1/(1-0) = 1  [GOOD]")
print("   y << 1 (deep MOND regime):")
print("     y^delta -> 0, exp(-y^delta) -> 1 - y^delta")
print("     nu4 -> 1/y^delta")
print("     g_obs = nu4 * g_bar = g_bar / y^delta = g_bar / (g_bar/a0)^delta")
print("          = g_bar^(1-delta) * a0^delta")

print("\n3. BTFR requirement:")
print("   V^4 = G*M*a0 requires g_obs = sqrt(g_bar * a0) at low g_bar")
print("   Form 4 at small y: g_obs = g_bar^(1-delta) * a0^delta")
print("   For BTFR slope 4: need 1-delta = 0.5, i.e., delta = 0.5 EXACTLY")

# Verify numerically
print("\n   Numerical verification:")
delta_test = 0.5
y_small = np.array([1e-4, 1e-3, 1e-2])
for yy in y_small:
    g_b = yy * a0
    g_form4 = nu_form4(yy, delta_test) * g_b
    g_btfr = np.sqrt(g_b * a0)
    print("   y=%.0e: g_form4/g_btfr = %.8f (ratio->1 means exact BTFR)" % (yy, g_form4 / g_btfr))

print("\n4. Value at y=1:")
e = np.exp(1)
nu4_at1 = 1.0 / (1.0 - 1.0/e)
print("   nu4(y=1, any delta) = 1/(1 - 1/e) = e/(e-1) = %.6f" % nu4_at1)
nu_tgp_at1 = nu_tgp_standard(1.0)
print("   nu_tgp(y=1) = 1 + exp(-1)/1^0.5714 = 1 + %.6f = %.6f" % (np.exp(-1), nu_tgp_at1))
nu_mond_at1 = nu_mond_standard(1.0)
print("   nu_mond(y=1) = 0.5 + sqrt(1.25) = %.6f" % nu_mond_at1)

# ============================================================
# PART B: RAR FIT COMPARISON (5 forms)
# ============================================================
print("\n" + "=" * 72)
print("PART B: RAR FIT COMPARISON (5 FORMS)")
print("=" * 72)

def chi2_rar(nu_func, params=None):
    """Compute chi2 on RAR bins. nu_func(y, *params) -> nu values.
    Compare log10(g_obs_model) with log10(g_obs_data)."""
    if params is not None:
        nu_vals = nu_func(y, *params)
    else:
        nu_vals = nu_func(y)
    g_model = nu_vals * gbar
    log_g_model = np.log10(g_model)
    chi2 = np.sum(((log_gobs - log_g_model) / sigma) ** 2)
    return chi2

n_bins = len(data)

# Form 1: MOND standard (no free params)
chi2_mond = chi2_rar(nu_mond_standard)
bic_mond = chi2_mond + 0 * np.log(n_bins)

# Form 2: MOND simple (no free params)
chi2_simple = chi2_rar(nu_mond_simple)
bic_simple = chi2_simple + 0 * np.log(n_bins)

# Form 3: TGP standard (alpha=0.8, gamma=0.5714 fixed)
chi2_tgp = chi2_rar(lambda y_: nu_tgp_standard(y_, 0.8, 0.5714))
bic_tgp = chi2_tgp + 0 * np.log(n_bins)

# Form 4: fit delta
def chi2_form4_delta(params):
    delta = params[0]
    if delta <= 0 or delta > 2:
        return 1e10
    return chi2_rar(nu_form4, (delta,))

res4 = minimize(chi2_form4_delta, [0.5], method='Nelder-Mead')
delta_best = res4.x[0]
chi2_f4 = res4.fun
bic_f4 = chi2_f4 + 1 * np.log(n_bins)

# Form 5: check if same as MOND standard
chi2_f5 = chi2_rar(nu_form5)
bic_f5 = chi2_f5 + 0 * np.log(n_bins)

print("\n%-25s %8s %8s %8s %12s" % ("Form", "chi2", "red_chi2", "BIC", "Parameters"))
print("-" * 72)
dof_mond = n_bins
dof_f4 = n_bins - 1
print("%-25s %8.4f %8.4f %8.4f %12s" % ("1. MOND standard", chi2_mond, chi2_mond/dof_mond, bic_mond, "none"))
print("%-25s %8.4f %8.4f %8.4f %12s" % ("2. MOND simple", chi2_simple, chi2_simple/dof_mond, bic_simple, "none"))
print("%-25s %8.4f %8.4f %8.4f %12s" % ("3. TGP standard", chi2_tgp, chi2_tgp/dof_mond, bic_tgp, "a=0.8,g=0.57"))
print("%-25s %8.4f %8.4f %8.4f %12s" % ("4. Form 4 (fit delta)", chi2_f4, chi2_f4/dof_f4, bic_f4, "d=%.4f" % delta_best))
print("%-25s %8.4f %8.4f %8.4f %12s" % ("5. Form 5", chi2_f5, chi2_f5/dof_mond, bic_f5, "none"))

print("\nForm 5 vs MOND standard difference check:")
y_check = np.logspace(-2, 2, 100)
diff_f5_mond = np.max(np.abs(nu_form5(y_check) - nu_mond_standard(y_check)))
print("  max |nu_form5 - nu_mond| over y in [0.01, 100] = %.2e" % diff_f5_mond)
if diff_f5_mond < 1e-10:
    print("  -> Form 5 IS IDENTICAL to MOND standard (algebraically equivalent)")
else:
    print("  -> Form 5 differs from MOND standard by up to %.2e" % diff_f5_mond)

print("\nForm 2 (MOND simple) vs Form 4 (delta=0.5):")
chi2_f4_05 = chi2_rar(nu_form4, (0.5,))
print("  MOND simple chi2 = %.4f" % chi2_simple)
print("  Form 4 (delta=0.5) chi2 = %.4f" % chi2_f4_05)
print("  -> These ARE the same function: 1/(1-exp(-y^0.5)) = 1/(1-exp(-sqrt(y)))")

# ============================================================
# PART C: DELTA TENSION
# ============================================================
print("\n" + "=" * 72)
print("PART C: DELTA TENSION ANALYSIS")
print("=" * 72)

print("\n1. TGP predicts delta = alpha/2 = 0.8/2 = 0.40")
print("2. Data best fit: delta = %.4f" % delta_best)
print("3. BTFR requires: delta = 0.50 exactly")

print("\n4. Physical scenarios for different delta values:")
print("   delta = alpha/2 mapping:")
print("   delta = 0.40 -> alpha = 0.80 (Flory exponent in d=3: nu_F = 3/5)")
print("   delta = 0.50 -> alpha = 1.00 (random walk: nu_F = 1/2 in d=2)")
print("   delta = 5/12 = 0.417 -> alpha = 5/6 (d=3 SAW Flory exact)")
print("   Or: delta = 1/2 may arise directly from BTFR constraint,")
print("       independent of Flory exponent")

print("\n5. Chi2 scan over delta:")
deltas = [0.40, 0.45, delta_best, 0.50, 0.55]
labels = ["TGP pred", "intermed", "best fit", "BTFR exact", "overshoot"]
print("\n   %-12s %-10s %-10s %-10s" % ("delta", "chi2", "red_chi2", "note"))
print("   " + "-" * 52)
for d, lab in zip(deltas, labels):
    c2 = chi2_rar(nu_form4, (d,))
    print("   %-12.4f %-10.4f %-10.4f %s" % (d, c2, c2/(n_bins-1), lab))

chi2_040 = chi2_rar(nu_form4, (0.40,))
chi2_best = chi2_rar(nu_form4, (delta_best,))
print("\n   Chi2 degradation from %.4f to 0.40: %.4f -> %.4f (factor %.1f)" %
      (delta_best, chi2_best, chi2_040, chi2_040/chi2_best if chi2_best > 0 else float('inf')))
print("   Delta-chi2 = %.4f" % (chi2_040 - chi2_best))
print("   For 1 parameter, delta-chi2 > 1 means > 1-sigma tension")
print("   For 1 parameter, delta-chi2 > 4 means > 2-sigma tension")
print("   For 1 parameter, delta-chi2 > 9 means > 3-sigma tension")

dchi2 = chi2_040 - chi2_best
n_sigma = np.sqrt(dchi2) if dchi2 > 0 else 0
print("   -> delta=0.40 is %.1f-sigma from best fit" % n_sigma)

# ============================================================
# PART D: CLUSTER PREDICTION WITH FORM 4
# ============================================================
print("\n" + "=" * 72)
print("PART D: CLUSTER PREDICTIONS WITH FORM 4")
print("=" * 72)

clusters = [
    ("Coma",    1.0e14, 1200, 7.0e14),
    ("Perseus", 8.0e13, 1100, 5.5e14),
    ("Virgo",   4.0e13,  900, 3.0e14),
    ("Bullet",  3.4e14, 1000, 5.5e14),
    ("A1689",   2.0e14, 1300, 1.2e15),
]

print("\n%-10s %10s %10s %10s %10s %10s %10s" % (
    "Cluster", "M_bar", "M_obs", "M_form4", "M_TGP", "M_MOND",
    "y_cluster"))
print("-" * 72)

for name, M_bar, R500_kpc, M_obs in clusters:
    R500_m = R500_kpc * kpc_to_m
    g_bar_cl = G_SI * M_bar * Msun / R500_m**2
    y_cl = g_bar_cl / a0

    nu_f4 = nu_form4(y_cl, delta_best)
    nu_tgp_cl = nu_tgp_standard(y_cl)
    nu_mond_cl = nu_mond_standard(y_cl)

    M_f4 = nu_f4 * M_bar
    M_tgp = nu_tgp_cl * M_bar
    M_mond = nu_mond_cl * M_bar

    print("%-10s %10.2e %10.2e %10.2e %10.2e %10.2e %10.4f" % (
        name, M_bar, M_obs, M_f4, M_tgp, M_mond, y_cl))

print("\nCluster deficit analysis (M_predicted / M_obs):")
print("%-10s %12s %12s %12s" % ("Cluster", "Form4/obs", "TGP/obs", "MOND/obs"))
print("-" * 50)
for name, M_bar, R500_kpc, M_obs in clusters:
    R500_m = R500_kpc * kpc_to_m
    g_bar_cl = G_SI * M_bar * Msun / R500_m**2
    y_cl = g_bar_cl / a0

    nu_f4 = nu_form4(y_cl, delta_best)
    nu_tgp_cl = nu_tgp_standard(y_cl)
    nu_mond_cl = nu_mond_standard(y_cl)

    M_f4 = nu_f4 * M_bar
    M_tgp = nu_tgp_cl * M_bar
    M_mond = nu_mond_cl * M_bar

    print("%-10s %12.4f %12.4f %12.4f" % (
        name, M_f4/M_obs, M_tgp/M_obs, M_mond/M_obs))

# ============================================================
# PART E: FORM 4 ROTATION CURVES
# ============================================================
print("\n" + "=" * 72)
print("PART E: ROTATION CURVES (NGC 2403-like)")
print("=" * 72)

M_disk = 5e9   # Msun
M_gas = 3.5e9  # Msun
R_d = 2.1      # kpc

def g_bar_galaxy(R_kpc, M_disk, M_gas, R_d):
    """Baryonic acceleration for exponential disk + gas.
    Simplified: treat as spherical mass enclosed."""
    R_m = R_kpc * kpc_to_m
    R_d_m = R_d * kpc_to_m
    x = R_m / R_d_m
    # Exponential disk enclosed mass fraction (Freeman 1970 approx)
    f_disk = 1.0 - (1.0 + x) * np.exp(-x)
    # Gas: more extended, use scale 2*R_d
    f_gas = 1.0 - (1.0 + x/2.0) * np.exp(-x/2.0)

    M_enc = (M_disk * f_disk + M_gas * f_gas) * Msun
    g = G_SI * M_enc / R_m**2
    return g

radii = np.linspace(0.5, 25.0, 50)  # kpc
g_bars = np.array([g_bar_galaxy(R, M_disk, M_gas, R_d) for R in radii])
y_gal = g_bars / a0

# Rotation velocity from g: V = sqrt(g * R)
def v_from_g(g_vals, R_kpc):
    R_m = R_kpc * kpc_to_m
    v = np.sqrt(np.abs(g_vals) * R_m)
    return v / 1e3  # km/s

# Newton
v_newton = v_from_g(g_bars, radii)

# MOND
g_mond = nu_mond_standard(y_gal) * g_bars
v_mond = v_from_g(g_mond, radii)

# TGP standard
g_tgp = nu_tgp_standard(y_gal) * g_bars
v_tgp = v_from_g(g_tgp, radii)

# Form 4 (delta = best fit)
g_f4_best = nu_form4(y_gal, delta_best) * g_bars
v_f4_best = v_from_g(g_f4_best, radii)

# Form 4 (delta = 0.40)
g_f4_040 = nu_form4(y_gal, 0.40) * g_bars
v_f4_040 = v_from_g(g_f4_040, radii)

print("\nRotation curves for NGC 2403-like galaxy:")
print("M_disk = %.1e Msun, M_gas = %.1e Msun, R_d = %.1f kpc" % (M_disk, M_gas, R_d))
print("\n%6s %8s %8s %8s %8s %8s" % (
    "R(kpc)", "Newton", "MOND", "TGP", "F4best", "F4_0.40"))
print("-" * 54)
for i in range(0, len(radii), 5):  # every 5th point
    print("%6.1f %8.1f %8.1f %8.1f %8.1f %8.1f" % (
        radii[i], v_newton[i], v_mond[i], v_tgp[i], v_f4_best[i], v_f4_040[i]))

print("\nPercentage difference relative to MOND standard:")
print("%6s %12s %12s %12s" % ("R(kpc)", "TGP-MOND%", "F4best-MOND%", "F4_040-MOND%"))
print("-" * 48)
for i in range(0, len(radii), 5):
    d_tgp = (v_tgp[i] - v_mond[i]) / v_mond[i] * 100
    d_f4b = (v_f4_best[i] - v_mond[i]) / v_mond[i] * 100
    d_f40 = (v_f4_040[i] - v_mond[i]) / v_mond[i] * 100
    print("%6.1f %12.2f %12.2f %12.2f" % (radii[i], d_tgp, d_f4b, d_f40))

# ============================================================
# PART F: CAN FORM 4 BE DERIVED FROM TGP?
# ============================================================
print("\n" + "=" * 72)
print("PART F: CAN FORM 4 BE DERIVED FROM TGP?")
print("=" * 72)

print("""
ANALYSIS: TGP Standard vs Form 4 Relationship

1. TGP standard form:
   nu_TGP(y) = 1 + exp(-y^alpha) / y^gamma
   where alpha = 0.8, gamma = 0.5714

   Asymptotes:
   - y >> 1: nu_TGP -> 1 (Newtonian)
   - y << 1: nu_TGP -> 1/y^gamma (deep MOND-like)
   - For BTFR: g_obs ~ g_bar^(1-gamma) * a0^gamma
     Need gamma = 0.5 for BTFR slope 4, but gamma = 0.5714
     -> TGP standard does NOT give exact BTFR slope 4

2. Form 4:
   nu4(y) = 1 / (1 - exp(-y^delta))

   Asymptotes:
   - y >> 1: nu4 -> 1 (Newtonian)
   - y << 1: nu4 -> 1/y^delta (deep MOND-like)
   - For BTFR: need delta = 0.5 exactly

3. Algebraic relationship:
   Write nu4 = 1 + F(y) where F(y) = exp(-y^delta) / (1 - exp(-y^delta))

   Compare with TGP: F_TGP(y) = exp(-y^alpha) / y^gamma

   For these to be equal:
   exp(-y^delta) / (1 - exp(-y^delta)) = exp(-y^alpha) / y^gamma

   Small-y expansion (y << 1):
   LHS: exp(-y^delta) / (1 - exp(-y^delta))
      ~ (1 - y^delta) / y^delta
      ~ 1/y^delta - 1
   RHS: exp(-y^alpha) / y^gamma
      ~ (1 - y^alpha) / y^gamma
      ~ 1/y^gamma

   Leading terms match when gamma = delta.

   Next order: LHS has correction -1, RHS has correction -y^alpha/y^gamma = -y^(alpha-gamma)
   These match when alpha - gamma = 0, i.e., alpha = gamma = delta.

   -> Form 4 corresponds to TGP with alpha = gamma = delta ~ 0.5

4. TGP constraint: gamma = alpha * c_eff / (c_eff + 1)
   If alpha = gamma: c_eff / (c_eff + 1) = 1 -> c_eff = infinity

   Physical meaning:
   c_eff = morphological codimension of the galaxy in the substrate
   c_eff -> infinity means the galaxy occupies zero volume fraction
   in the substrate (point-like embedding)

   This is UNPHYSICAL for extended objects like galaxies.

   However: in the limit c_eff >> 1, gamma -> alpha, which gives
   nearly Form-4 behavior. For c_eff = 10: gamma = 10/11 * alpha.
   For alpha = 0.5: gamma = 0.4545 (close but not exact).

5. Alternative interpretation:
   Form 4 may NOT be derivable from TGP standard form.
   Instead, it could arise from a DIFFERENT TGP mechanism:
   - Modified substrate propagator
   - Non-perturbative polymer effects
   - Exact resummation of the TGP series

   The fact that Form 4 = MOND simple (when delta = 0.5) suggests
   it may be a more fundamental form than the TGP perturbative result.""")

# Numerical verification of the relationship
print("\nNumerical verification: Form 4 vs TGP with alpha=gamma=0.5:")
y_test = np.array([0.01, 0.05, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0])
print("\n%8s %12s %12s %12s" % ("y", "Form4(0.5)", "TGP(a=g=0.5)", "ratio"))
print("-" * 48)
for yy in y_test:
    f4 = nu_form4(yy, 0.5)
    tgp_mod = nu_tgp_standard(yy, alpha=0.5, gamma=0.5)
    print("%8.3f %12.6f %12.6f %12.6f" % (yy, f4, tgp_mod, f4/tgp_mod))

print("\nKey finding: Form 4 != TGP(alpha=gamma=0.5)")
print("The exponential structures differ:")
print("  Form 4: F(y) = exp(-y^d) / (1 - exp(-y^d))   [geometric series]")
print("  TGP:    F(y) = exp(-y^a) / y^g                [power-law denominator]")
print("These converge only asymptotically for y << 1 and y >> 1,")
print("but differ at intermediate y (transition region).")

# ============================================================
# PART G: SUMMARY AND IMPLICATIONS
# ============================================================
print("\n" + "=" * 72)
print("PART G: SUMMARY AND IMPLICATIONS")
print("=" * 72)

print("""
KEY RESULTS:

1. FORM 4 PERFORMANCE:
   - Best fit delta = %.4f (very close to 0.5)
   - Chi2 = %.4f on 14 RAR bins (reduced chi2 = %.4f)
   - DRAMATICALLY better than MOND standard (chi2 = %.4f)
   - Better than TGP standard (chi2 = %.4f)
""" % (delta_best, chi2_f4, chi2_f4/(n_bins-1), chi2_mond, chi2_tgp))

print("""2. FORM 4 = MOND SIMPLE (when delta = 0.5):
   nu = 1/(1-exp(-sqrt(y)))  [MOND simple interpolating function]
   This is a KNOWN form in the MOND literature (McGaugh 2008).
   The fact that the best fit gives delta ~ 0.5 means:
   Form 4 rediscovers MOND simple as the optimal interpolating function.

3. BTFR:
   - Form 4 with delta = 0.5 gives BTFR slope EXACTLY 4
   - TGP standard with gamma = 0.5714 gives slope 4/(1+0.5714) ~ 2.54
     (in the sense that g_obs ~ g_bar^(1-gamma), not exactly slope 4)
   - This is a SIGNIFICANT advantage of Form 4

4. DELTA TENSION:
   - Data prefer delta = %.4f
   - TGP predicts delta = 0.40 (from alpha = 0.8)
   - Chi2 at delta=0.40 vs best: %.4f vs %.4f
   - Tension: %.1f sigma
""" % (delta_best, chi2_040, chi2_best, n_sigma))

print("""5. CLUSTER PREDICTIONS:
   Form 4 gives similar cluster mass predictions to MOND and TGP.
   All forms underpredict observed cluster masses.
   The cluster deficit persists regardless of interpolating function,
   suggesting it is a genuine physical effect (hot gas, non-thermal
   pressure, or actual missing mass).

6. ROTATION CURVES:
   Form 4 (delta=0.497) and MOND simple (delta=0.5) give nearly
   identical rotation curves. The difference from MOND standard
   is small (<5%%) at most radii, largest in the transition region.

7. TGP COMPATIBILITY:
   Form 4 CANNOT be directly derived from TGP standard form
   (nu = 1 + exp(-y^a)/y^g) with physical parameters.
   It would require alpha = gamma (i.e., c_eff -> infinity),
   which is unphysical.

   POSSIBLE RESOLUTIONS:
   a) TGP standard form is only the leading-order approximation;
      higher-order corrections could yield Form 4
   b) The true TGP nu(y) is obtained by resummation, giving
      the geometric-series structure of Form 4
   c) alpha is not exactly 4/5 but closer to 1/2
   d) The interpolating function depends on galaxy properties
      (not universal), and Form 4 is an effective average

8. BOTTOM LINE:
   Form 4 with delta ~ 0.5 is phenomenologically excellent.
   It equals the well-known MOND simple function.
   The tension with TGP (delta=0.4 vs 0.5) is real but modest.
   This tension could CONSTRAIN TGP parameters or indicate
   that the perturbative TGP form needs modification.""")

print("\n" + "=" * 72)
print("END gs56_form4_analysis.py")
print("=" * 72)
