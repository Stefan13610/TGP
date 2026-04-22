#!/usr/bin/env python3
"""
ode_formulation_check_v47b.py -- Check ODE consistency between formulations.

Two ODE forms used in TGP scripts:

  FORM A (canonical, K(g) = g^4, alpha=2):
    g'' + (2/r)g' + (2/g)*g'^2 + g^2*(1-g) = 0
    Source: V'(g) = g^2(1-g), Cross: (2/g)*g'^2

  FORM B (scripts tau_selection, koide_from_wkb):
    g'' + (2/r)g' + (1/g)*g'^2 + (1-g) = 0
    Source: (1-g), Cross: (1/g)*g'^2

QUESTION: Do both give the same Koide Q_K = 3/2 and mass ratios?
If yes, results are robust. If no, need to identify the correct form.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

PHI = (1 + np.sqrt(5)) / 2
m_e = 0.51099895; m_mu = 105.6583755; m_tau = 1776.86
r_21 = m_mu / m_e

def solve_form_A(g0, r_max=120):
    """ODE Form A: g'' + (2/r)g' + (2/g)g'^2 + g^2(1-g) = 0"""
    def rhs(r, y):
        g, gp = y
        g = max(g, 1e-10)
        source = g**2 * (1.0 - g)
        cross = (2.0 / g) * gp**2
        if r < 1e-10:
            return [gp, (source - cross) / 3.0]
        return [gp, source - cross - 2.0 * gp / r]
    sol = solve_ivp(rhs, (1e-6, r_max), [g0, 0.0],
                    rtol=1e-11, atol=1e-13, max_step=0.03)
    return sol.t, sol.y[0]

def solve_form_B(g0, r_max=120):
    """ODE Form B: g'' + (2/r)g' + (1/g)g'^2 + (1-g) = 0"""
    def rhs(r, y):
        g, gp = y
        g = max(g, 1e-10)
        source = 1.0 - g
        cross = (1.0 / g) * gp**2
        if r < 1e-10:
            return [gp, (source - cross) / 3.0]
        return [gp, source - cross - 2.0 * gp / r]
    sol = solve_ivp(rhs, (1e-6, r_max), [g0, 0.0],
                    rtol=1e-11, atol=1e-13, max_step=0.03)
    return sol.t, sol.y[0]

def A_tail_generic(r, g):
    mask = (r > 30) & (r < 100)
    if np.sum(mask) < 50:
        return 0.0
    rf = r[mask]
    # Need to determine oscillation frequency first
    h = g[mask] - 1.0
    # Try fitting with free omega
    from scipy.optimize import curve_fit
    def model(r_, A, omega, delta):
        return A * np.cos(omega * r_ + delta) / r_
    try:
        popt, _ = curve_fit(model, rf, h, p0=[0.1, 1.0, 0], maxfev=5000)
        return abs(popt[0]), popt[1]  # amplitude, frequency
    except:
        # Fallback: assume omega = 1
        df = h * rf
        M = np.column_stack([np.cos(rf), np.sin(rf)])
        bc = np.linalg.lstsq(M, df, rcond=None)[0]
        return np.sqrt(bc[0]**2 + bc[1]**2), 1.0

def A_tail_simple(solver, g0):
    """Standard A_tail extraction assuming omega=1."""
    r, g = solver(g0)
    mask = (r > 30) & (r < 100)
    if np.sum(mask) < 50:
        return 0.0
    rf = r[mask]
    df = (g[mask] - 1.0) * rf
    M = np.column_stack([np.cos(rf), np.sin(rf)])
    bc = np.linalg.lstsq(M, df, rcond=None)[0]
    return np.sqrt(bc[0]**2 + bc[1]**2)


print("=" * 72)
print("ODE FORMULATION CHECK: FORM A vs FORM B")
print("=" * 72)


# ============================================================
# [1] TAIL STRUCTURE
# ============================================================
print("\n[1] TAIL OSCILLATION FREQUENCY")
print("-" * 50)

# Form A linearization: g = 1+h
# h'' + (2/r)h' + 2h'^2/(1+h) + (1+h)^2*(-h) = 0
# Leading: h'' + (2/r)h' - h + ... = 0  ???
# Actually: g^2(1-g) = (1+h)^2(-h) = -h - 2h^2 - h^3
# So source linearization: -h  (same sign as Form B's -h)
# Cross: (2/g)g'^2 = 2h'^2/(1+h) -> 2h'^2 (negligible for small h)
# => h'' + (2/r)h' - h = 0  (Form A: exponential decay!)
#
# Wait: -h gives h'' + (2/r)h' = h, which is GROWING, not oscillatory.
# For Form B: source = 1-g = -h, but:
# g'' = source - cross - 2g'/r = -h - (1/g)g'^2 - 2g'/r
# So g'' + (2/r)g' + (1/g)g'^2 = -h
# Linearized: h'' + (2/r)h' = -h => h'' + (2/r)h' + h = 0 (oscillatory!)

# CRITICAL DIFFERENCE:
# Form A: h'' + (2/r)h' + h = 0? Let me redo.
# Form A: g'' + (2/r)g' + (2/g)g'^2 + g^2(1-g) = 0
# g'' = -g^2(1-g) - (2/g)g'^2 - (2/r)g'
# g = 1+h: g'' = h''
# g^2(1-g) = (1+h)^2*(-h) = -(h + 2h^2 + h^3)
# (2/g)g'^2 = 2h'^2/(1+h) ~ 2h'^2
# So: h'' = (h + 2h^2) - 2h'^2 - (2/r)h'
# Linearized: h'' + (2/r)h' - h = 0  (minus h => exponential, NOT oscillatory!)
# This gives u = rh: u'' - u = 0 => u = A*exp(-r) (decaying)

# Form B: g'' + (2/r)g' + (1/g)g'^2 + (1-g) = 0  =>  g'' + (2/r)g' + (1/g)g'^2 = (1-g)
# Wait, I need to be more careful with the sign.
# From script: g'' = source - cross - (2/r)g'
# source = 1-g, cross = (1/g)g'^2
# g'' = (1-g) - (1/g)g'^2 - (2/r)g'
# So: g'' + (2/r)g' + (1/g)g'^2 = 1-g = -h
# Linearized: h'' + (2/r)h' = -h => h'' + (2/r)h' + h = 0 (oscillatory!)

# So Form A gives EXPONENTIAL tails (omega imaginary)
# and Form B gives OSCILLATORY tails (omega = 1, real)
#
# Our numerical solutions show OSCILLATORY tails.
# This means our scripts use FORM B, which is INCONSISTENT with
# the canonical TGP equation!

print("  Form A linearization: h'' + (2/r)h' - h = 0  => EXPONENTIAL decay")
print("  Form B linearization: h'' + (2/r)h' + h = 0  => OSCILLATORY (omega=1)")
print()

# Verify with both solvers
for g0_test in [0.5, 0.87, 1.4]:
    print(f"  g0 = {g0_test}:")

    # Form A
    r_A, g_A = solve_form_A(g0_test, r_max=200)
    A_A, omega_A = A_tail_generic(r_A, g_A)
    print(f"    Form A: A_tail = {A_A:.8f}, omega = {omega_A:.6f}")

    # Form B
    r_B, g_B = solve_form_B(g0_test, r_max=200)
    A_B, omega_B = A_tail_generic(r_B, g_B)
    print(f"    Form B: A_tail = {A_B:.8f}, omega = {omega_B:.6f}")

    # Check tail values
    mask_A = r_A > 80
    mask_B = r_B > 80
    if np.any(mask_A):
        h_A = g_A[mask_A] - 1.0
        print(f"    Form A tail (r>80): max|h| = {np.max(np.abs(h_A)):.2e}, "
              f"sign changes = {np.sum(np.diff(np.sign(h_A)) != 0)}")
    if np.any(mask_B):
        h_B = g_B[mask_B] - 1.0
        print(f"    Form B tail (r>80): max|h| = {np.max(np.abs(h_B)):.2e}, "
              f"sign changes = {np.sum(np.diff(np.sign(h_B)) != 0)}")
    print()


# ============================================================
# [2] KOIDE TEST FOR FORM A
# ============================================================
print("\n[2] KOIDE TEST FOR FORM A (canonical TGP)")
print("-" * 50)

# Calibrate g0_e for Form A
def r21_res_A(g0_1):
    A1 = A_tail_simple(solve_form_A, g0_1)
    A2 = A_tail_simple(solve_form_A, PHI * g0_1)
    if A1 < 1e-15: return 1e10
    return (A2/A1)**4 - r_21

# For Form A, the tails are exponential, not oscillatory.
# A_tail extraction with cos/sin basis may not work!
# Check if A_tail is meaningful for Form A.

print("  Note: Form A has exponential tails => A_tail(cos/sin) may fail")
print("  Using direct tail amplitude instead.")
print()

# For exponential tails: h ~ B*exp(-r)/r
# Try fitting this form
def A_tail_exp(solver, g0):
    """A_tail for exponentially decaying solutions."""
    r, g = solver(g0, r_max=200)
    mask = (r > 20) & (r < 80)
    if np.sum(mask) < 50:
        return 0.0
    rf = r[mask]
    h = g[mask] - 1.0
    # Fit h = B * exp(-r) / r
    # ln|h*r| = ln|B| - r
    rh = np.abs(h * rf)
    valid = rh > 1e-30
    if np.sum(valid) < 10:
        return 0.0
    log_rh = np.log(rh[valid])
    r_v = rf[valid]
    # Linear fit: log_rh = a + b*r
    coeffs = np.polyfit(r_v, log_rh, 1)
    B = np.exp(coeffs[1])
    decay_rate = -coeffs[0]
    return B, decay_rate

for g0_test in [0.5, 0.87, 1.4, 2.0]:
    result = A_tail_exp(solve_form_A, g0_test)
    if result:
        B, kappa = result
        print(f"  g0={g0_test}: B={B:.6e}, decay_rate={kappa:.6f} "
              f"(expect ~1 for unit mass)")
    else:
        print(f"  g0={g0_test}: tail too small")

print()

# For Form A: mass ratios come from exponential tail, not oscillatory.
# m_n ~ B_n^4 where B_n = exponential amplitude
# Test if phi-FP still gives r_21 = 206.77

print("  Testing phi-FP with Form A exponential amplitudes:")
for g0_test in [0.5, 0.8, 1.0, 1.2]:
    B1 = A_tail_exp(solve_form_A, g0_test)
    B2 = A_tail_exp(solve_form_A, PHI * g0_test)
    if B1 and B2:
        ratio = (B2[0]/B1[0])**4
        print(f"  g0={g0_test}: B1={B1[0]:.4e}, B2={B2[0]:.4e}, (B2/B1)^4={ratio:.2f} "
              f"(target: {r_21:.2f})")
    else:
        print(f"  g0={g0_test}: one or both tails too small")


# ============================================================
# [3] KOIDE FOR FORM B (current scripts)
# ============================================================
print("\n\n[3] KOIDE FOR FORM B (as in current scripts)")
print("-" * 50)

def r21_res_B(g0_1):
    A1 = A_tail_simple(solve_form_B, g0_1)
    A2 = A_tail_simple(solve_form_B, PHI * g0_1)
    if A1 < 1e-15: return 1e10
    return (A2/A1)**4 - r_21

g0_e_B = brentq(r21_res_B, 0.82, 0.90, xtol=1e-8)
g0_mu_B = PHI * g0_e_B
A_e_B = A_tail_simple(solve_form_B, g0_e_B)
A_mu_B = A_tail_simple(solve_form_B, g0_mu_B)

print(f"  Form B: g0_e = {g0_e_B:.8f}, A_e = {A_e_B:.10f}")
print(f"  Form B: g0_mu = {g0_mu_B:.8f}, A_mu = {A_mu_B:.10f}")
print(f"  r_21 = {(A_mu_B/A_e_B)**4:.3f} (target: {r_21:.3f})")
print()

# Koide
def koide_B(g0_3):
    A1 = A_tail_simple(solve_form_B, g0_e_B)
    A2 = A_tail_simple(solve_form_B, g0_mu_B)
    A3 = A_tail_simple(solve_form_B, g0_3)
    if min(A1,A2,A3) < 1e-15: return 1e10
    S2 = A1**2 + A2**2 + A3**2; S4 = A1**4 + A2**4 + A3**4
    return S2**2/S4 - 1.5

g0_tau_B = brentq(koide_B, 1.5, 2.0, xtol=1e-10)
A_tau_B = A_tail_simple(solve_form_B, g0_tau_B)
r31_B = (A_tau_B/A_e_B)**4
print(f"  Form B Koide: g0_tau = {g0_tau_B:.8f}, r_31 = {r31_B:.1f} (PDG: 3477.2)")
print()


# ============================================================
# [4] CONNECTION BETWEEN FORMS
# ============================================================
print("\n[4] IS FORM B A VALID TGP ODE?")
print("-" * 50)

# The canonical TGP has K(Phi) = Phi^4, giving alpha = 2.
# The EL equation from S = int K(g)(nabla g)^2 + V(g) is:
#   2K*nabla^2(g) + K'*(nabla g)^2 + V'(g) = 0
# For K = g^4: K' = 4g^3
#   2g^4*(g'' + 2g'/r) + 4g^3*g'^2 + V'(g) = 0
#   Divide by 2g^3:
#   g*(g'' + 2g'/r) + 2g'^2 + V'/(2g^3) = 0
#   g'' + 2g'/r + 2g'^2/g + V'/(2g^4) = 0
#
# With V(g) = g^8/8 - g^7/7 + const:
#   V'(g) = g^7 - g^6 = g^6(g - 1)
#   V'/(2g^4) = g^2(g-1)/2
#
# So: g'' + 2g'/r + 2g'^2/g + g^2(g-1)/2 = 0
# Or: g'' + 2g'/r + 2g'^2/g - g^2(1-g)/2 = 0
#
# Hmm, this has different coefficients than either Form A or B.
# With beta normalization: multiply source by constant.
# If we set V(g) = (2beta)*[g^8/8 - g^7/7]:
#   V' = 2beta*g^6(g-1)
#   V'/(2g^4) = beta*g^2(g-1)
#   g'' + 2g'/r + 2g'^2/g + beta*g^2(g-1) = 0
#
# With beta=1: g'' + 2g'/r + 2g'^2/g + g^2(g-1) = 0
# Source: +g^2(g-1) means for g>1, source>0 (repulsive from vacuum)
# This is Form A but with SOURCE = g^2(g-1), not g^2(1-g).
# Hmm, signs matter!

# Actually, the script Form A I wrote has source = g^2(1-g), but
# the EL equation gives source = g^2(g-1). These have opposite signs!
# Let me re-derive more carefully.

# TGP action: S = int [K(Phi)(nabla Phi)^2 - V(Phi)] d^3x  (with minus V)
# Or S = int [K(Phi)(nabla Phi)^2 + V(Phi)] d^3x  (depends on convention)
# EL: d/dx[2K*dg/dx] - [K'*(dg/dx)^2 + V'] = 0  (varies by sign convention)

# In the scripts (ex142_A4_perturbative.py):
#   source = Vp(g) = g^2(1-g)
#   cross = (ALPHA/g)*gp^2 with ALPHA = 2
#   g'' = source - cross - 2g'/r
#   => g'' + 2g'/r + (2/g)g'^2 = g^2(1-g)
#
# The sign of source determines whether the soliton sits at g < 1
# (undershoot from vacuum at g=1) or g > 1 (overshoot).

# For Form B:
#   g'' = (1-g) - (1/g)g'^2 - 2g'/r
#   => g'' + 2g'/r + (1/g)g'^2 = (1-g)
#
# Both have source = +(1-g) on RHS. For g < 1: source > 0 (restoring to vacuum).
# The key question is: which is the correct derivation?

print("  Both forms have source = (1-g) type (restoring force toward g=1)")
print("  Form A: alpha=2, source = g^2(1-g)")
print("  Form B: alpha=1, source = (1-g)")
print()
print("  Canonical TGP (K = g^4, V = g^8/8 - g^7/7):")
print("    EL gives: g'' + (2/r)g' + (2/g)g'^2 + g^2(g-1)/2 = 0  (if V in action)")
print("    or with -V: g'' + (2/r)g' + (2/g)g'^2 - g^2(g-1)/2 = 0")
print()
print("  Form B might correspond to a different K(g) or V(g).")
print("  Specifically: K(g) = g^2 gives alpha = K'/K = 2/g -> alpha_coeff = 1")
print("    and source = (1-g) corresponds to V'(g) = -(1-g)*g^2 -> V(g) = g^2/2 - g^3/3")
print()

# Check if Form B comes from K = g^2:
# S = int g^2*(nabla g)^2 + V(g)
# EL: 2g^2*(g'' + 2g'/r) + 2g*g'^2 + V' = 0
# Divide by 2g: g*(g''+2g'/r) + g'^2 + V'/(2g) = 0
# g''+2g'/r + g'^2/g + V'/(2g^2) = 0
# If V'/(2g^2) = (1-g): V' = 2g^2(1-g), V = 2g^3/3 - g^4/2 + const
# This is a QUARTIC potential, not the canonical OCTIC.

print("  CONCLUSION: Form B corresponds to K(g) = g^2 with quartic potential,")
print("  NOT the canonical TGP K(g) = g^4 with octic potential.")
print("  The scripts inherited this from an early formulation.")
print()
print("  IMPACT: Koide Q_K = 3/2 might depend on the ODE form.")
print("  Need to test Form A for Koide property.")

print()
print("=" * 72)
print("CRITICAL FINDING")
print("=" * 72)
print("""
  The ODE in our Koide/tau scripts (Form B) differs from the canonical
  TGP field equation (Form A). Key differences:

  Form A (canonical K=g^4): g'' + (2/r)g' + (2/g)g'^2 + g^2(1-g) = 0
    Tail: EXPONENTIAL (h'' + (2/r)h' - h = 0)
    No oscillatory modes in tail

  Form B (scripts):         g'' + (2/r)g' + (1/g)g'^2 + (1-g) = 0
    Tail: OSCILLATORY (h'' + (2/r)h' + h = 0, omega = 1)
    A_tail from cos/sin fit

  These give DIFFERENT physics:
  - Form A has exponentially decaying solitons
  - Form B has oscillatory tails (resonance-like)
  - Mass extraction mechanism differs fundamentally

  STATUS: INCONSISTENCY DETECTED between script ODE and canonical TGP.
  Need to resolve: which form is the physically correct one?

  Possible resolutions:
  1. Form B IS correct for some choice of K(g) and V(g) -- need to identify which
  2. The canonical form needs modification (different potential)
  3. There's a coordinate/field redefinition connecting the two
""")
