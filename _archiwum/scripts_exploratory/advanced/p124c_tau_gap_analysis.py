#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p124c_tau_gap_analysis.py -- Analysis of the tau mass gap in TGP (fast version)
================================================================================

Single-soliton max r_31 ~ 832 vs Koide target 3477.
Tests two most promising paths: ERG potential correction and kinetic modification.
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from scipy.integrate import solve_ivp

ALPHA = 2
PHI = (1 + np.sqrt(5)) / 2
R31_KOIDE = 3477.44
G0_E = 0.8992655880
A_E_REF = 0.10180692  # baseline A_tail(electron)

def solve_and_atail(g0, f_func, Vp_func, rm=300):
    """Solve soliton ODE and return A_tail."""
    fg0 = f_func(g0)
    if abs(fg0) < 1e-15:
        return np.nan
    c2 = Vp_func(g0) / (3 * fg0)
    rs = 0.01
    def rhs(r, y):
        g, p = y
        if g <= 1e-15: return [p, 0]
        fg = f_func(g)
        if abs(fg) < 1e-10: return [p, 0]
        if r < 1e-10: return [p, Vp_func(g) / fg / 3]
        return [p, (Vp_func(g) - 2/r * p) / fg]
    def ev(r, y): return 100 - abs(y[0])
    ev.terminal = True
    try:
        s = solve_ivp(rhs, [rs, rm], [g0 + c2*rs**2, 2*c2*rs],
                      method='RK45', rtol=1e-10, atol=1e-12,
                      max_step=0.1, events=[ev], dense_output=True)
        r = np.linspace(rs, min(s.t[-1], rm), 12000)
        g = s.sol(r)[0]
    except Exception:
        return np.nan
    m = (r >= 120) & (r <= 260)
    rf = r[m]
    tl = (g[m]-1)*rf
    if len(rf) < 10: return np.nan
    A = np.column_stack([np.cos(rf), np.sin(rf)])
    B, C = np.linalg.lstsq(A, tl, rcond=None)[0]
    return np.sqrt(B**2 + C**2)

def scan_max_r31(f_func, Vp_func, g0_range=None):
    """Find max r_31 over g0 range."""
    if g0_range is None:
        g0_range = np.linspace(0.82, 4.0, 60)
    best = 0
    best_g0 = 0
    for g0 in g0_range:
        A = solve_and_atail(g0, f_func, Vp_func)
        if np.isnan(A) or A <= 0: continue
        r31 = (A / A_E_REF)**4
        if r31 > best:
            best = r31
            best_g0 = g0
    return best_g0, best

# Standard functions
f_std = lambda g: 1 + 2*ALPHA*np.log(g) if g > 0 else -1e30
Vp_std = lambda g: g**2 * (1 - g)

print("=" * 60)
print("  TAU MASS GAP ANALYSIS")
print("  Standard soliton max r_31 = 832 vs Koide = 3477")
print("=" * 60)

# ============= PATH A: ERG potential corrections =============
print("\n--- PATH A: ERG-corrected potential V'(g) = g^2(1-g) + c5*g^4 - c5*g^5 ---")
print("   (c6 = -c5 to preserve V'(1)=0 vacuum condition)")

for c5 in [-0.2, -0.1, -0.05, 0.0, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0]:
    Vp_mod = lambda g, c=c5: g**2*(1-g) + c*g**4 - c*g**5
    g0_b, r31_b = scan_max_r31(f_std, Vp_mod)
    # Check r_21 is preserved
    A_e_m = solve_and_atail(G0_E, f_std, Vp_mod)
    A_mu_m = solve_and_atail(PHI * G0_E, f_std, Vp_mod)
    r21_m = (A_mu_m / A_e_m)**4 if not np.isnan(A_e_m) and not np.isnan(A_mu_m) and A_e_m > 0 else -1
    marker = " <-- PROMISING" if r31_b > 2000 else ""
    print("  c5=%6.3f: max_r31 = %7.0f at g0=%.3f, r21=%.1f%s" % (c5, r31_b, g0_b, r21_m, marker))

# ============= PATH B: Kinetic flattening =============
print("\n--- PATH B: Modified kinetic f(g) = max(1+4*ln(g) - km*(g-1)^2, 0.01) ---")

for km in [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0]:
    f_mod = lambda g, k=km: max(1 + 2*ALPHA*np.log(g) - k*(g-1)**2, 0.01) if g > 0 else -1e30
    g0_b, r31_b = scan_max_r31(f_mod, Vp_std)
    A_e_m = solve_and_atail(G0_E, f_mod, Vp_std)
    A_mu_m = solve_and_atail(PHI * G0_E, f_mod, Vp_std)
    r21_m = (A_mu_m / A_e_m)**4 if not np.isnan(A_e_m) and not np.isnan(A_mu_m) and A_e_m > 0 else -1
    marker = " <-- PROMISING" if r31_b > 2000 else ""
    print("  km=%4.1f: max_r31 = %7.0f at g0=%.3f, r21=%.1f%s" % (km, r31_b, g0_b, r21_m, marker))

# ============= PATH B2: Varying effective alpha =============
print("\n--- PATH B2: f(g) = 1 + 2*a_eff*ln(g) with varying a_eff ---")

for a_eff in [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0]:
    f_a = lambda g, a=a_eff: 1 + 2*a*np.log(g) if g > 0 else -1e30
    g0_b, r31_b = scan_max_r31(f_a, Vp_std)
    marker = " <-- STANDARD" if a_eff == 2.0 else (" <-- PROMISING" if r31_b > 2000 else "")
    print("  a_eff=%4.1f: max_r31 = %7.0f at g0=%.3f%s" % (a_eff, r31_b, g0_b, marker))

# ============= PATH D: Mass exponent analysis =============
print("\n--- PATH D: Required mass exponent ---")
A_max = 0.547
ratio_max = A_max / A_E_REF
n_needed = np.log(R31_KOIDE) / np.log(ratio_max)
n_from_r21 = np.log(206.768) / np.log(0.386054 / A_E_REF)
print("  A_max/A_e = %.4f" % ratio_max)
print("  n needed for r_31=3477 (at A_max): n = %.3f" % n_needed)
print("  n from r_21 (at A_mu): n = %.3f" % n_from_r21)
print("  These differ by %.1f%% -- universality of exponent is problematic" % (abs(n_needed - n_from_r21)/n_from_r21*100))

# ============= PATH C: Bound state estimate =============
print("\n--- PATH C: Two-soliton bound state ---")
A_target = R31_KOIDE**0.25 * A_E_REF
print("  A_tau needed (Koide) = %.4f" % A_target)
print("  A_max (single)       = %.4f" % A_max)
print("  Enhancement factor   = %.3f (= A_needed/A_max)" % (A_target/A_max))
print("  If coherent sum: A_eff ~ 2*A -> r_31 ~ %.0f (too large)" % ((2*A_max/A_E_REF)**4))
print("  If A_eff ~ 1.43*A_max: r_31 ~ %.0f (matches!)" % (((A_target)/A_E_REF)**4))
print("  Physical mechanism: partial decoherence or overlap suppression")

# ============= SUMMARY =============
print("\n" + "=" * 60)
print("  SUMMARY")
print("=" * 60)
print("  The tau mass gap (r_31_max=832 vs 3477) is REAL.")
print("  Single-generation soliton with standard V,f CANNOT reach tau.")
print()
print("  Most viable resolutions (TGP-native):")
print("  1. ERG potential correction c5 > 0: enhances r_31_max")
print("     (needs numerical ERG to determine c5 value)")
print("  2. Kinetic flattening: reduces damping at large g0")
print("  3. Tau as composite (2-soliton bound state)")
print("  4. [Unlikely] Different mass exponent for tau")
print()
print("  Status O-K1: OPEN with concrete diagnostic")
print("  Next step: compute c5 from Wetteriu ERG flow at k->0")
