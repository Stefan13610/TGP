#!/usr/bin/env python3
"""
gs48: ALTERNATIVE nu(y) FUNCTIONAL FORMS FROM MEMBRANE PHYSICS
===============================================================

Investigates whether the TGP chi2/dof ~ 2x MOND gap (from gs37/gs44) is due
to the specific functional form of nu(y) or is structural (from alpha=4/5).

Tests 5 alternative membrane-derived nu(y) forms against:
  - MOND simple interpolation function
  - SPARC RAR data (McGaugh+2016 binned means)

Sections:
  A. Current nu(y) shape analysis (TGP vs MOND)
  B. Alternative membrane-derived forms (5 candidates)
  C. RAR comparison -- fit free parameters to binned RAR data
  D. Asymptotic behavior check (deep-MOND, Newtonian, BTFR)
  E. Verdict -- which form (if any) closes the gap?
"""

import numpy as np
from scipy.optimize import minimize
import sys
import os

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

# --- Constants ---
a0 = 1.2e-10  # m/s^2
ALPHA = 4.0 / 5.0  # Flory exponent, fixed throughout

# ================================================================
#  FORMATTING UTILITIES
# ================================================================

def print_header(title):
    print()
    print("=" * 78)
    print(f"  {title}")
    print("=" * 78)
    print()


def print_subheader(title):
    print(f"  {title}")
    print(f"  {'-' * len(title)}")
    print()


# ================================================================
#  CORE nu(y) FUNCTIONS
# ================================================================

def nu_tgp(y, c_eff=1.5):
    """Current TGP: nu(y) = 1 + exp(-y^alpha) / y^gamma"""
    y = np.asarray(y, dtype=float)
    y = np.maximum(y, 1e-15)
    alpha = ALPHA
    gamma = alpha * c_eff / (c_eff + 1)
    return 1.0 + np.exp(-y**alpha) / y**gamma


def nu_mond_simple(y):
    """MOND simple interpolation: nu = 0.5*(1 + sqrt(1 + 4/y))"""
    y = np.asarray(y, dtype=float)
    y = np.maximum(y, 1e-15)
    return 0.5 * (1.0 + np.sqrt(1.0 + 4.0 / y))


def nu_mond_standard(y):
    """MOND standard interpolation: nu = 1/sqrt(1 - exp(-sqrt(y)))"""
    y = np.asarray(y, dtype=float)
    y = np.maximum(y, 1e-15)
    arg = np.sqrt(y)
    return 1.0 / np.sqrt(1.0 - np.exp(-arg))


# --- Alternative forms ---

def nu_form1(y, beta1, c_eff=1.5):
    """Form 1 (Generalized exponential):
    nu1 = 1 + exp(-y^alpha)/y^gamma * (1 + beta1 * y^alpha)
    """
    y = np.asarray(y, dtype=float)
    y = np.maximum(y, 1e-15)
    alpha = ALPHA
    gamma = alpha * c_eff / (c_eff + 1)
    ya = y**alpha
    return 1.0 + np.exp(-ya) / y**gamma * (1.0 + beta1 * ya)


def nu_form2(y, A1, A2, c_eff=1.5):
    """Form 2 (Double-exponential / two membrane modes):
    nu2 = 1 + A1*exp(-y^alpha)/y^gamma + A2*exp(-y^(2*alpha))/y^(2*gamma)
    """
    y = np.asarray(y, dtype=float)
    y = np.maximum(y, 1e-15)
    alpha = ALPHA
    gamma = alpha * c_eff / (c_eff + 1)
    return (1.0 + A1 * np.exp(-y**alpha) / y**gamma
            + A2 * np.exp(-y**(2*alpha)) / y**(2*gamma))


def nu_form3(y, c_eff=1.5):
    """Form 3 (Modified suppression):
    nu3 = 1 + (1 - exp(-y^(-alpha/2)))^2 / y^gamma
    """
    y = np.asarray(y, dtype=float)
    y = np.maximum(y, 1e-15)
    alpha = ALPHA
    gamma = alpha * c_eff / (c_eff + 1)
    suppression = (1.0 - np.exp(-y**(-alpha / 2)))**2
    return 1.0 + suppression / y**gamma


def nu_form4(y, c_eff=1.5):
    """Form 4 (RAR-matched):
    nu4 = 1 / (1 - exp(-y^(alpha/2)))
    Note: c_eff enters only through interpretation; functional form is fixed.
    """
    y = np.asarray(y, dtype=float)
    y = np.maximum(y, 1e-15)
    alpha = ALPHA
    arg = y**(alpha / 2)
    val = 1.0 - np.exp(-arg)
    val = np.maximum(val, 1e-30)
    return 1.0 / val


def nu_form4_gen(y, delta):
    """Form 4 generalized: nu = 1/(1 - exp(-y^delta))
    delta is free; motivated by alpha/2 = 0.4 but let it vary.
    """
    y = np.asarray(y, dtype=float)
    y = np.maximum(y, 1e-15)
    arg = y**delta
    val = 1.0 - np.exp(-arg)
    val = np.maximum(val, 1e-30)
    return 1.0 / val


def nu_form5(y, kappa, c_eff=1.5):
    """Form 5 (Substrate viscosity):
    nu5 = 1 + exp(-y^alpha) / (y^gamma * (1 + kappa*y^(1-alpha)))
    """
    y = np.asarray(y, dtype=float)
    y = np.maximum(y, 1e-15)
    alpha = ALPHA
    gamma = alpha * c_eff / (c_eff + 1)
    return 1.0 + np.exp(-y**alpha) / (y**gamma * (1.0 + kappa * y**(1.0 - alpha)))


# ================================================================
#  SPARC RAR DATA (McGaugh+2016 binned means from RARbins.mrt)
# ================================================================

# Binned RAR data: log10(gbar), log10(gobs), scatter, N
RAR_BINS = np.array([
    [-11.69, -10.75, 0.18,  75],
    [-11.48, -10.70, 0.16, 191],
    [-11.26, -10.57, 0.15, 308],
    [-11.05, -10.43, 0.14, 323],
    [-10.84, -10.31, 0.13, 283],
    [-10.62, -10.17, 0.12, 229],
    [-10.41, -10.03, 0.14, 242],
    [-10.19,  -9.90, 0.13, 229],
    [ -9.98,  -9.76, 0.14, 204],
    [ -9.76,  -9.60, 0.14, 177],
    [ -9.55,  -9.44, 0.14, 137],
    [ -9.34,  -9.29, 0.11, 104],
    [ -9.12,  -9.11, 0.13,  85],
    [ -8.91,  -8.90, 0.11,  43],
])

log_gbar = RAR_BINS[:, 0]
log_gobs_data = RAR_BINS[:, 1]
log_gobs_sigma = RAR_BINS[:, 2]
N_bins = RAR_BINS[:, 3]

# Convert to linear for nu computation
gbar = 10.0**log_gbar
gobs_data = 10.0**log_gobs_data
y_rar = gbar / a0  # dimensionless acceleration


# ================================================================
#  PART A: CURRENT nu(y) SHAPE ANALYSIS
# ================================================================

def part_a():
    print_header("PART A: CURRENT nu(y) SHAPE ANALYSIS")

    y = np.logspace(-3, 2, 500)

    nu_t = nu_tgp(y, c_eff=1.5)
    nu_m_simple = nu_mond_simple(y)
    nu_m_standard = nu_mond_standard(y)

    # Use standard MOND as primary comparison (matches the user's formula)
    nu_m = nu_m_standard
    delta_nu = nu_t - nu_m
    ratio = nu_t / nu_m

    print_subheader("A.1  nu(y) values at key accelerations")
    print(f"  {'y':>10s}  {'nu_TGP':>10s}  {'nu_MOND_std':>12s}  {'nu_MOND_sim':>12s}  {'delta':>10s}  {'ratio':>8s}")
    print(f"  {'-'*10}  {'-'*10}  {'-'*12}  {'-'*12}  {'-'*10}  {'-'*8}")
    for yv in [0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0]:
        idx = np.argmin(np.abs(y - yv))
        print(f"  {yv:10.3f}  {nu_t[idx]:10.4f}  {nu_m[idx]:12.4f}  {nu_m_simple[idx]:12.4f}  "
              f"{delta_nu[idx]:10.4f}  {ratio[idx]:8.4f}")

    # Logarithmic slope: d(ln nu)/d(ln y)
    print()
    print_subheader("A.2  Logarithmic slopes d(ln nu)/d(ln y)")
    ln_y = np.log(y)
    ln_nu_t = np.log(nu_t)
    ln_nu_m = np.log(nu_m)
    slope_t = np.gradient(ln_nu_t, ln_y)
    slope_m = np.gradient(ln_nu_m, ln_y)

    print(f"  {'y':>10s}  {'slope_TGP':>10s}  {'slope_MOND':>11s}  {'diff':>10s}")
    print(f"  {'-'*10}  {'-'*10}  {'-'*11}  {'-'*10}")
    for yv in [0.001, 0.01, 0.1, 0.3, 0.5, 1.0, 3.0, 10.0, 100.0]:
        idx = np.argmin(np.abs(y - yv))
        print(f"  {yv:10.3f}  {slope_t[idx]:10.4f}  {slope_m[idx]:11.4f}  {slope_t[idx]-slope_m[idx]:10.4f}")

    # Identify where discrepancy is largest
    print()
    print_subheader("A.3  Largest discrepancies")
    idx_max_abs = np.argmax(np.abs(delta_nu))
    idx_max_ratio = np.argmax(np.abs(ratio - 1.0))
    # Focus on intermediate regime y in [0.01, 10]
    mask_mid = (y >= 0.01) & (y <= 10.0)
    y_mid = y[mask_mid]
    delta_mid = delta_nu[mask_mid]
    ratio_mid = ratio[mask_mid]
    idx_mid_max = np.argmax(np.abs(delta_mid))

    print(f"    Global max |delta_nu|:  y = {y[idx_max_abs]:.4f}, delta = {delta_nu[idx_max_abs]:.4f}")
    print(f"    Global max |ratio-1|:   y = {y[idx_max_ratio]:.4f}, ratio = {ratio[idx_max_ratio]:.4f}")
    print(f"    Intermediate (0.01-10) max |delta|: y = {y_mid[idx_mid_max]:.4f}, delta = {delta_mid[idx_mid_max]:.4f}")

    # Deep MOND asymptotic slopes
    print()
    print_subheader("A.4  Asymptotic deep-MOND behavior (y -> 0)")
    y_deep = np.logspace(-3, -1.5, 50)
    nu_t_deep = nu_tgp(y_deep, c_eff=1.5)
    nu_m_deep = nu_m_standard[y < 0.04]  # approximate
    # Fit power law nu ~ C * y^(-p) for y << 1
    from numpy.polynomial import polynomial as P
    log_y_deep = np.log10(y_deep)
    log_nu_t_deep = np.log10(nu_t_deep)
    coeff_t = np.polyfit(log_y_deep, log_nu_t_deep, 1)
    print(f"    TGP:  nu ~ y^({coeff_t[0]:.4f})  [expected gamma = {0.8*1.5/2.5:.4f}]")

    y_deep2 = np.logspace(-3, -1.5, 50)
    nu_m_deep2 = nu_mond_standard(y_deep2)
    log_nu_m_deep2 = np.log10(nu_m_deep2)
    coeff_m = np.polyfit(np.log10(y_deep2), log_nu_m_deep2, 1)
    print(f"    MOND: nu ~ y^({coeff_m[0]:.4f})  [expected -0.5]")


# ================================================================
#  PART B: ALTERNATIVE MEMBRANE-DERIVED FORMS
# ================================================================

def part_b():
    print_header("PART B: ALTERNATIVE MEMBRANE-DERIVED FORMS")

    y = np.logspace(-3, 2, 500)
    nu_m = nu_mond_standard(y)
    nu_t = nu_tgp(y, c_eff=1.5)

    forms = {
        "Current TGP": nu_t,
        "Form 1 (beta1=0.5)": nu_form1(y, beta1=0.5, c_eff=1.5),
        "Form 1 (beta1=1.0)": nu_form1(y, beta1=1.0, c_eff=1.5),
        "Form 2 (A1=1,A2=0.5)": nu_form2(y, A1=1.0, A2=0.5, c_eff=1.5),
        "Form 2 (A1=1,A2=1.0)": nu_form2(y, A1=1.0, A2=1.0, c_eff=1.5),
        "Form 3": nu_form3(y, c_eff=1.5),
        "Form 4 (delta=0.4)": nu_form4_gen(y, delta=0.4),
        "Form 4 (delta=0.5)": nu_form4_gen(y, delta=0.5),
        "Form 5 (kappa=0.5)": nu_form5(y, kappa=0.5, c_eff=1.5),
        "Form 5 (kappa=1.0)": nu_form5(y, kappa=1.0, c_eff=1.5),
        "MOND standard": nu_m,
    }

    print_subheader("B.1  nu(y) comparison at key y values")
    y_check = [0.01, 0.1, 0.3, 1.0, 3.0, 10.0]
    print(f"  {'Form':<26s}", end="")
    for yv in y_check:
        print(f"  {'y='+str(yv):>8s}", end="")
    print()
    print(f"  {'-'*26}", end="")
    for _ in y_check:
        print(f"  {'-'*8}", end="")
    print()

    for name, nu_vals in forms.items():
        print(f"  {name:<26s}", end="")
        for yv in y_check:
            idx = np.argmin(np.abs(y - yv))
            print(f"  {nu_vals[idx]:8.3f}", end="")
        print()

    # RMS deviation from MOND
    print()
    print_subheader("B.2  RMS deviation from MOND standard (log space, y in [0.01, 10])")
    mask = (y >= 0.01) & (y <= 10.0)
    print(f"  {'Form':<26s}  {'RMS(log10)':>10s}  {'max|dlog|':>10s}")
    print(f"  {'-'*26}  {'-'*10}  {'-'*10}")
    for name, nu_vals in forms.items():
        if name == "MOND standard":
            continue
        log_diff = np.log10(nu_vals[mask]) - np.log10(nu_m[mask])
        rms = np.sqrt(np.mean(log_diff**2))
        maxd = np.max(np.abs(log_diff))
        print(f"  {name:<26s}  {rms:10.4f}  {maxd:10.4f}")


# ================================================================
#  PART C: RAR COMPARISON -- FIT TO BINNED DATA
# ================================================================

def chi2_rar(nu_func, params, param_names):
    """Compute chi2 against binned RAR data."""
    try:
        kwargs = {name: val for name, val in zip(param_names, params)}
        nu_vals = nu_func(y_rar, **kwargs)
        gobs_pred = nu_vals * gbar
        log_gobs_pred = np.log10(np.maximum(gobs_pred, 1e-20))
        residuals = (log_gobs_pred - log_gobs_data) / log_gobs_sigma
        return np.sum(residuals**2)
    except (ValueError, RuntimeWarning, FloatingPointError):
        return 1e10


def part_c():
    print_header("PART C: RAR COMPARISON -- FIT FREE PARAMETERS TO BINNED DATA")

    n_data = len(log_gbar)
    print(f"  Using {n_data} binned RAR data points from McGaugh+2016")
    print(f"  a0 = {a0:.2e} m/s^2")
    print()

    # --- MOND baselines ---
    print_subheader("C.1  MOND baselines (no free parameters)")

    for name, func in [("MOND simple", nu_mond_simple), ("MOND standard", nu_mond_standard)]:
        nu_vals = func(y_rar)
        gobs_pred = nu_vals * gbar
        log_pred = np.log10(gobs_pred)
        chi2 = np.sum(((log_pred - log_gobs_data) / log_gobs_sigma)**2)
        print(f"    {name:<20s}: chi2 = {chi2:8.2f}, chi2/dof = {chi2/n_data:6.3f}")

    # --- Current TGP with c_eff free ---
    print()
    print_subheader("C.2  Current TGP (c_eff free)")

    def chi2_tgp_ceff(params):
        c_eff = params[0]
        if c_eff <= 0.1 or c_eff > 10:
            return 1e10
        nu_vals = nu_tgp(y_rar, c_eff=c_eff)
        gobs_pred = nu_vals * gbar
        log_pred = np.log10(np.maximum(gobs_pred, 1e-20))
        return np.sum(((log_pred - log_gobs_data) / log_gobs_sigma)**2)

    best_tgp = None
    best_chi2_tgp = 1e10
    for c0 in [0.5, 1.0, 1.5, 2.0, 3.0, 5.0]:
        res = minimize(chi2_tgp_ceff, [c0], method='Nelder-Mead',
                       options={'xatol': 1e-6, 'fatol': 1e-8, 'maxiter': 5000})
        if res.fun < best_chi2_tgp:
            best_chi2_tgp = res.fun
            best_tgp = res.x

    c_eff_best = best_tgp[0]
    gamma_best = ALPHA * c_eff_best / (c_eff_best + 1)
    print(f"    Best c_eff = {c_eff_best:.4f}")
    print(f"    gamma = {gamma_best:.4f}")
    print(f"    BTFR slope = {2.0/(1.0-gamma_best):.3f}")
    print(f"    chi2 = {best_chi2_tgp:.2f}, chi2/dof = {best_chi2_tgp/(n_data-1):.3f}")

    # --- Form 1: Generalized exponential ---
    print()
    print_subheader("C.3  Form 1: Generalized exponential (beta1, c_eff free)")

    def chi2_form1(params):
        beta1, c_eff = params
        if c_eff <= 0.1 or c_eff > 10 or beta1 < -5 or beta1 > 20:
            return 1e10
        try:
            nu_vals = nu_form1(y_rar, beta1=beta1, c_eff=c_eff)
            if np.any(nu_vals <= 0) or np.any(~np.isfinite(nu_vals)):
                return 1e10
            gobs_pred = nu_vals * gbar
            log_pred = np.log10(np.maximum(gobs_pred, 1e-20))
            return np.sum(((log_pred - log_gobs_data) / log_gobs_sigma)**2)
        except:
            return 1e10

    best_f1 = None
    best_chi2_f1 = 1e10
    for b0 in [0.0, 0.5, 1.0, 2.0, 5.0]:
        for c0 in [0.5, 1.0, 1.5, 2.0, 3.0]:
            res = minimize(chi2_form1, [b0, c0], method='Nelder-Mead',
                           options={'xatol': 1e-6, 'fatol': 1e-8, 'maxiter': 10000})
            if res.fun < best_chi2_f1:
                best_chi2_f1 = res.fun
                best_f1 = res.x

    print(f"    Best beta1 = {best_f1[0]:.4f}, c_eff = {best_f1[1]:.4f}")
    gamma_f1 = ALPHA * best_f1[1] / (best_f1[1] + 1)
    print(f"    gamma = {gamma_f1:.4f}")
    print(f"    BTFR slope = {2.0/(1.0-gamma_f1):.3f}")
    print(f"    chi2 = {best_chi2_f1:.2f}, chi2/dof = {best_chi2_f1/(n_data-2):.3f}")

    # --- Form 2: Double-exponential ---
    print()
    print_subheader("C.4  Form 2: Double-exponential (A1, A2, c_eff free)")

    def chi2_form2(params):
        A1, A2, c_eff = params
        if c_eff <= 0.1 or c_eff > 10 or A1 < 0 or A2 < -10 or A1 > 20 or A2 > 20:
            return 1e10
        try:
            nu_vals = nu_form2(y_rar, A1=A1, A2=A2, c_eff=c_eff)
            if np.any(nu_vals <= 0) or np.any(~np.isfinite(nu_vals)):
                return 1e10
            gobs_pred = nu_vals * gbar
            log_pred = np.log10(np.maximum(gobs_pred, 1e-20))
            return np.sum(((log_pred - log_gobs_data) / log_gobs_sigma)**2)
        except:
            return 1e10

    best_f2 = None
    best_chi2_f2 = 1e10
    for a1 in [0.5, 1.0, 2.0]:
        for a2 in [-1.0, 0.0, 0.5, 1.0, 2.0]:
            for c0 in [1.0, 1.5, 2.0, 3.0]:
                res = minimize(chi2_form2, [a1, a2, c0], method='Nelder-Mead',
                               options={'xatol': 1e-6, 'fatol': 1e-8, 'maxiter': 10000})
                if res.fun < best_chi2_f2:
                    best_chi2_f2 = res.fun
                    best_f2 = res.x

    print(f"    Best A1 = {best_f2[0]:.4f}, A2 = {best_f2[1]:.4f}, c_eff = {best_f2[2]:.4f}")
    gamma_f2 = ALPHA * best_f2[2] / (best_f2[2] + 1)
    print(f"    gamma = {gamma_f2:.4f}")
    print(f"    BTFR slope = {2.0/(1.0-gamma_f2):.3f}")
    print(f"    chi2 = {best_chi2_f2:.2f}, chi2/dof = {best_chi2_f2/(n_data-3):.3f}")

    # --- Form 3: Modified suppression (c_eff free) ---
    print()
    print_subheader("C.5  Form 3: Modified suppression (c_eff free)")

    def chi2_form3(params):
        c_eff = params[0]
        if c_eff <= 0.1 or c_eff > 10:
            return 1e10
        try:
            nu_vals = nu_form3(y_rar, c_eff=c_eff)
            if np.any(nu_vals <= 0) or np.any(~np.isfinite(nu_vals)):
                return 1e10
            gobs_pred = nu_vals * gbar
            log_pred = np.log10(np.maximum(gobs_pred, 1e-20))
            return np.sum(((log_pred - log_gobs_data) / log_gobs_sigma)**2)
        except:
            return 1e10

    best_f3 = None
    best_chi2_f3 = 1e10
    for c0 in [0.5, 1.0, 1.5, 2.0, 3.0, 5.0]:
        res = minimize(chi2_form3, [c0], method='Nelder-Mead',
                       options={'xatol': 1e-6, 'fatol': 1e-8, 'maxiter': 5000})
        if res.fun < best_chi2_f3:
            best_chi2_f3 = res.fun
            best_f3 = res.x

    gamma_f3 = ALPHA * best_f3[0] / (best_f3[0] + 1)
    print(f"    Best c_eff = {best_f3[0]:.4f}")
    print(f"    gamma = {gamma_f3:.4f}")
    print(f"    BTFR slope = {2.0/(1.0-gamma_f3):.3f}")
    print(f"    chi2 = {best_chi2_f3:.2f}, chi2/dof = {best_chi2_f3/(n_data-1):.3f}")

    # --- Form 4: RAR-matched (delta free) ---
    print()
    print_subheader("C.6  Form 4: RAR-matched nu=1/(1-exp(-y^delta)), delta free")

    def chi2_form4(params):
        delta = params[0]
        if delta <= 0.01 or delta > 3.0:
            return 1e10
        try:
            nu_vals = nu_form4_gen(y_rar, delta=delta)
            if np.any(nu_vals <= 0) or np.any(~np.isfinite(nu_vals)):
                return 1e10
            gobs_pred = nu_vals * gbar
            log_pred = np.log10(np.maximum(gobs_pred, 1e-20))
            return np.sum(((log_pred - log_gobs_data) / log_gobs_sigma)**2)
        except:
            return 1e10

    best_f4 = None
    best_chi2_f4 = 1e10
    for d0 in [0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0]:
        res = minimize(chi2_form4, [d0], method='Nelder-Mead',
                       options={'xatol': 1e-8, 'fatol': 1e-10, 'maxiter': 5000})
        if res.fun < best_chi2_f4:
            best_chi2_f4 = res.fun
            best_f4 = res.x

    delta_best = best_f4[0]
    print(f"    Best delta = {delta_best:.6f}")
    print(f"    alpha/2 = {ALPHA/2:.4f}")
    print(f"    Deviation from alpha/2: {abs(delta_best - ALPHA/2):.4f}")
    # Deep MOND: nu ~ 1/(y^delta) so gamma_eff = delta
    print(f"    Effective deep-MOND slope = {delta_best:.4f}")
    print(f"    BTFR slope = {2.0/(1.0-delta_best):.3f}")
    print(f"    chi2 = {best_chi2_f4:.2f}, chi2/dof = {best_chi2_f4/(n_data-1):.3f}")

    # --- Form 5: Substrate viscosity (kappa, c_eff free) ---
    print()
    print_subheader("C.7  Form 5: Substrate viscosity (kappa, c_eff free)")

    def chi2_form5(params):
        kappa, c_eff = params
        if c_eff <= 0.1 or c_eff > 10 or kappa < -5 or kappa > 50:
            return 1e10
        try:
            nu_vals = nu_form5(y_rar, kappa=kappa, c_eff=c_eff)
            if np.any(nu_vals <= 0) or np.any(~np.isfinite(nu_vals)):
                return 1e10
            gobs_pred = nu_vals * gbar
            log_pred = np.log10(np.maximum(gobs_pred, 1e-20))
            return np.sum(((log_pred - log_gobs_data) / log_gobs_sigma)**2)
        except:
            return 1e10

    best_f5 = None
    best_chi2_f5 = 1e10
    for k0 in [0.0, 0.5, 1.0, 2.0, 5.0]:
        for c0 in [0.5, 1.0, 1.5, 2.0, 3.0]:
            res = minimize(chi2_form5, [k0, c0], method='Nelder-Mead',
                           options={'xatol': 1e-6, 'fatol': 1e-8, 'maxiter': 10000})
            if res.fun < best_chi2_f5:
                best_chi2_f5 = res.fun
                best_f5 = res.x

    gamma_f5 = ALPHA * best_f5[1] / (best_f5[1] + 1)
    print(f"    Best kappa = {best_f5[0]:.4f}, c_eff = {best_f5[1]:.4f}")
    print(f"    gamma = {gamma_f5:.4f}")
    print(f"    BTFR slope = {2.0/(1.0-gamma_f5):.3f}")
    print(f"    chi2 = {best_chi2_f5:.2f}, chi2/dof = {best_chi2_f5/(n_data-2):.3f}")

    # --- Summary table ---
    print()
    print_subheader("C.8  SUMMARY: chi2 comparison")

    # Recompute MOND baselines for table
    nu_ms = nu_mond_simple(y_rar)
    chi2_ms = np.sum(((np.log10(nu_ms * gbar) - log_gobs_data) / log_gobs_sigma)**2)
    nu_mst = nu_mond_standard(y_rar)
    chi2_mst = np.sum(((np.log10(nu_mst * gbar) - log_gobs_data) / log_gobs_sigma)**2)

    results = [
        ("MOND simple", chi2_ms, 0, n_data),
        ("MOND standard", chi2_mst, 0, n_data),
        ("TGP current", best_chi2_tgp, 1, n_data),
        ("Form 1 (gen. exp.)", best_chi2_f1, 2, n_data),
        ("Form 2 (double exp.)", best_chi2_f2, 3, n_data),
        ("Form 3 (mod. supp.)", best_chi2_f3, 1, n_data),
        ("Form 4 (RAR-matched)", best_chi2_f4, 1, n_data),
        ("Form 5 (viscosity)", best_chi2_f5, 2, n_data),
    ]

    print(f"  {'Model':<22s}  {'chi2':>8s}  {'n_par':>5s}  {'dof':>5s}  {'chi2/dof':>9s}  {'vs MOND':>8s}")
    print(f"  {'-'*22}  {'-'*8}  {'-'*5}  {'-'*5}  {'-'*9}  {'-'*8}")
    for name, chi2, npar, ndat in results:
        dof = ndat - npar
        ratio_vs_mond = chi2 / chi2_mst if chi2_mst > 0 else float('inf')
        print(f"  {name:<22s}  {chi2:8.2f}  {npar:5d}  {dof:5d}  {chi2/dof:9.3f}  {ratio_vs_mond:8.2f}x")

    # Store results for use by other parts
    return {
        'tgp': (best_tgp, best_chi2_tgp),
        'f1': (best_f1, best_chi2_f1),
        'f2': (best_f2, best_chi2_f2),
        'f3': (best_f3, best_chi2_f3),
        'f4': (best_f4, best_chi2_f4),
        'f5': (best_f5, best_chi2_f5),
    }


# ================================================================
#  PART D: ASYMPTOTIC BEHAVIOR CHECK
# ================================================================

def part_d(fit_results):
    print_header("PART D: ASYMPTOTIC BEHAVIOR CHECK")

    y_deep = np.logspace(-4, -2, 200)
    y_newt = np.logspace(1, 3, 200)

    print_subheader("D.1  Deep-MOND limit (y -> 0): effective power-law slope")

    models = []

    # Current TGP
    c_eff_tgp = fit_results['tgp'][0][0]
    nu_deep = nu_tgp(y_deep, c_eff=c_eff_tgp)
    slope = np.polyfit(np.log10(y_deep), np.log10(nu_deep), 1)[0]
    gamma_eff = -slope
    btfr = 2.0 / (1.0 - gamma_eff)
    models.append(("TGP current", gamma_eff, btfr))

    # Form 1
    b1, c1 = fit_results['f1'][0]
    nu_deep = nu_form1(y_deep, beta1=b1, c_eff=c1)
    slope = np.polyfit(np.log10(y_deep), np.log10(nu_deep), 1)[0]
    gamma_eff_f1 = -slope
    btfr_f1 = 2.0 / (1.0 - gamma_eff_f1)
    models.append(("Form 1 (gen. exp.)", gamma_eff_f1, btfr_f1))

    # Form 2
    a1, a2, c2 = fit_results['f2'][0]
    nu_deep = nu_form2(y_deep, A1=a1, A2=a2, c_eff=c2)
    slope = np.polyfit(np.log10(y_deep), np.log10(nu_deep), 1)[0]
    gamma_eff_f2 = -slope
    btfr_f2 = 2.0 / (1.0 - gamma_eff_f2)
    models.append(("Form 2 (double exp.)", gamma_eff_f2, btfr_f2))

    # Form 3
    c3 = fit_results['f3'][0][0]
    nu_deep = nu_form3(y_deep, c_eff=c3)
    slope = np.polyfit(np.log10(y_deep), np.log10(nu_deep), 1)[0]
    gamma_eff_f3 = -slope
    btfr_f3 = 2.0 / (1.0 - gamma_eff_f3)
    models.append(("Form 3 (mod. supp.)", gamma_eff_f3, btfr_f3))

    # Form 4
    delta4 = fit_results['f4'][0][0]
    nu_deep = nu_form4_gen(y_deep, delta=delta4)
    slope = np.polyfit(np.log10(y_deep), np.log10(nu_deep), 1)[0]
    gamma_eff_f4 = -slope
    btfr_f4 = 2.0 / (1.0 - gamma_eff_f4)
    models.append(("Form 4 (RAR-matched)", gamma_eff_f4, btfr_f4))

    # Form 5
    k5, c5 = fit_results['f5'][0]
    nu_deep = nu_form5(y_deep, kappa=k5, c_eff=c5)
    slope = np.polyfit(np.log10(y_deep), np.log10(nu_deep), 1)[0]
    gamma_eff_f5 = -slope
    btfr_f5 = 2.0 / (1.0 - gamma_eff_f5)
    models.append(("Form 5 (viscosity)", gamma_eff_f5, btfr_f5))

    # MOND
    nu_deep_m = nu_mond_standard(y_deep)
    slope_m = np.polyfit(np.log10(y_deep), np.log10(nu_deep_m), 1)[0]
    gamma_eff_m = -slope_m
    btfr_m = 2.0 / (1.0 - gamma_eff_m)
    models.append(("MOND standard", gamma_eff_m, btfr_m))

    print(f"  {'Model':<22s}  {'gamma_eff':>10s}  {'BTFR slope':>11s}  {'|BTFR-3.85|':>12s}")
    print(f"  {'-'*22}  {'-'*10}  {'-'*11}  {'-'*12}")
    for name, geff, btfr_val in models:
        dev = abs(btfr_val - 3.85)
        print(f"  {name:<22s}  {geff:10.4f}  {btfr_val:11.3f}  {dev:12.3f}")

    print()
    print_subheader("D.2  Newtonian limit (y -> inf): nu should -> 1")

    print(f"  {'Model':<22s}  {'nu(y=10)':>10s}  {'nu(y=100)':>11s}  {'nu(y=1000)':>11s}")
    print(f"  {'-'*22}  {'-'*10}  {'-'*11}  {'-'*11}")

    y_test = np.array([10.0, 100.0, 1000.0])

    c_eff_tgp = fit_results['tgp'][0][0]
    vals = nu_tgp(y_test, c_eff=c_eff_tgp)
    print(f"  {'TGP current':<22s}  {vals[0]:10.6f}  {vals[1]:11.6f}  {vals[2]:11.6f}")

    b1, c1 = fit_results['f1'][0]
    vals = nu_form1(y_test, beta1=b1, c_eff=c1)
    print(f"  {'Form 1':<22s}  {vals[0]:10.6f}  {vals[1]:11.6f}  {vals[2]:11.6f}")

    a1, a2, c2 = fit_results['f2'][0]
    vals = nu_form2(y_test, A1=a1, A2=a2, c_eff=c2)
    print(f"  {'Form 2':<22s}  {vals[0]:10.6f}  {vals[1]:11.6f}  {vals[2]:11.6f}")

    c3 = fit_results['f3'][0][0]
    vals = nu_form3(y_test, c_eff=c3)
    print(f"  {'Form 3':<22s}  {vals[0]:10.6f}  {vals[1]:11.6f}  {vals[2]:11.6f}")

    delta4 = fit_results['f4'][0][0]
    vals = nu_form4_gen(y_test, delta=delta4)
    print(f"  {'Form 4':<22s}  {vals[0]:10.6f}  {vals[1]:11.6f}  {vals[2]:11.6f}")

    k5, c5 = fit_results['f5'][0]
    vals = nu_form5(y_test, kappa=k5, c_eff=c5)
    print(f"  {'Form 5':<22s}  {vals[0]:10.6f}  {vals[1]:11.6f}  {vals[2]:11.6f}")

    vals = nu_mond_standard(y_test)
    print(f"  {'MOND standard':<22s}  {vals[0]:10.6f}  {vals[1]:11.6f}  {vals[2]:11.6f}")

    print()
    print_subheader("D.3  BTFR constraint assessment")
    print("    Observed BTFR slope: 3.85 +/- 0.09 (Lelli+2016)")
    print("    Acceptable range: [3.67, 4.03] (2-sigma)")
    print()
    for name, geff, btfr_val in models:
        in_range = "YES" if 3.67 <= btfr_val <= 4.03 else "NO"
        print(f"    {name:<22s}: BTFR = {btfr_val:.3f}  acceptable: {in_range}")


# ================================================================
#  PART E: VERDICT
# ================================================================

def part_e(fit_results):
    print_header("PART E: VERDICT")

    # Extract chi2 values
    chi2_tgp = fit_results['tgp'][1]
    chi2_f1 = fit_results['f1'][1]
    chi2_f2 = fit_results['f2'][1]
    chi2_f3 = fit_results['f3'][1]
    chi2_f4 = fit_results['f4'][1]
    chi2_f5 = fit_results['f5'][1]

    nu_mst = nu_mond_standard(y_rar)
    chi2_mond = np.sum(((np.log10(nu_mst * gbar) - log_gobs_data) / log_gobs_sigma)**2)

    best_alt_chi2 = min(chi2_f1, chi2_f2, chi2_f3, chi2_f4, chi2_f5)
    best_alt_name = ["Form 1", "Form 2", "Form 3", "Form 4", "Form 5"][
        [chi2_f1, chi2_f2, chi2_f3, chi2_f4, chi2_f5].index(best_alt_chi2)]

    print_subheader("E.1  Key findings")
    print(f"    MOND standard chi2:         {chi2_mond:.2f}")
    print(f"    Current TGP chi2:           {chi2_tgp:.2f}  (ratio: {chi2_tgp/chi2_mond:.2f}x)")
    print(f"    Best alternative chi2:      {best_alt_chi2:.2f}  ({best_alt_name}, ratio: {best_alt_chi2/chi2_mond:.2f}x)")
    print()

    improvement = (chi2_tgp - best_alt_chi2) / chi2_tgp * 100
    gap_closed = (chi2_tgp - best_alt_chi2) / (chi2_tgp - chi2_mond) * 100 if chi2_tgp > chi2_mond else 0

    print(f"    Improvement over current TGP:  {improvement:.1f}%")
    print(f"    Chi2 gap closed vs MOND:       {gap_closed:.1f}%")
    print()

    print_subheader("E.2  Is the gap structural or functional?")
    if best_alt_chi2 < chi2_mond * 1.1:
        print("    CONCLUSION: The chi2 gap is FUNCTIONAL, not structural.")
        print(f"    {best_alt_name} achieves chi2 within 10% of MOND while preserving alpha=4/5.")
        print("    The specific exponential form of TGP's nu(y) is the bottleneck,")
        print("    not the Flory exponent alpha=4/5 itself.")
    elif best_alt_chi2 < chi2_mond * 1.5:
        print("    CONCLUSION: The gap is PARTIALLY structural.")
        print(f"    {best_alt_name} reduces the gap significantly but cannot fully close it.")
        print("    Alpha=4/5 constrains the transition shape, but the functional form matters too.")
    else:
        print("    CONCLUSION: The gap is LARGELY STRUCTURAL.")
        print("    Even with alternative functional forms, alpha=4/5 prevents matching MOND.")
        print("    The Flory exponent constrains the UV-IR transition too rigidly.")

    print()
    print_subheader("E.3  Recommendation for TGP")

    # Check Form 4 specifically since it's the most MOND-like
    delta4 = fit_results['f4'][0][0]
    print(f"    Form 4 best-fit delta = {delta4:.4f}")
    print(f"    TGP prediction alpha/2 = {ALPHA/2:.4f}")
    print(f"    MOND uses sqrt(y) => delta = 0.5")
    print()

    if abs(delta4 - 0.5) < abs(delta4 - ALPHA/2):
        print("    The data PREFER delta ~ 0.5 (MOND) over delta = alpha/2 = 0.4 (TGP).")
    else:
        print("    The data PREFER delta ~ alpha/2 = 0.4 (TGP) over delta = 0.5 (MOND).")
    print()

    print("    RECOMMENDATIONS:")
    print("    1. If chi2 gap is functional: adopt the best-performing alternative form")
    print("       that preserves membrane physics motivation (alpha=4/5 fixed).")
    print("    2. If gap is structural: accept that TGP's rotation curve fits will be")
    print("       systematically ~2x MOND, or consider whether alpha could receive")
    print("       quantum/thermal corrections shifting it toward 1/2.")
    print("    3. Form 4 (nu = 1/(1-exp(-y^delta))) is the most promising template --")
    print("       it structurally matches MOND's form while allowing membrane-motivated delta.")
    print("    4. For rotation curve fitting, the intermediate y regime (0.1 < y < 1)")
    print("       is where all the leverage is -- focus on the transition shape there.")


# ================================================================
#  MAIN
# ================================================================

if __name__ == "__main__":
    print("=" * 78)
    print("  gs48: ALTERNATIVE nu(y) FUNCTIONAL FORMS FROM MEMBRANE PHYSICS")
    print("  Investigating whether functional form or alpha=4/5 drives the chi2 gap")
    print("=" * 78)

    import warnings
    warnings.filterwarnings('ignore', category=RuntimeWarning)

    part_a()
    part_b()
    fit_results = part_c()
    part_d(fit_results)
    part_e(fit_results)

    print()
    print("=" * 78)
    print("  gs48 COMPLETE")
    print("=" * 78)
