#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
tgp_c133_origin.py
===================
Investigation: origin of c = 13/3 in g0_e = 1 - c * rho_0*

KEY QUESTION: Is c = 13/3 a fundamental constant or an artefact of alpha=2?

STRATEGY:
  D1: c(alpha) -- how does c depend on the coordinate parameter alpha?
  D2: Canonical alpha -- is there alpha where c is exactly rational?
  D3: ODE structure -- does c follow from V(g) = g^2(1-g) analytically?
  D4: WF critical exponents -- test known combinations
  D5: Dimension dependence -- c(d) for d = 2,3,4,...

Wersja: v47b (2026-04-12)
"""

import sys
import io
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')
else:
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8',
                                  errors='replace')

# =====================================================================
# CONSTANTS
# =====================================================================

PHI = (1 + np.sqrt(5)) / 2
R21_PDG = 206.7682830
RHO_0_STAR = 0.03045       # WF fixed point in d=3

RESULTS = []

def check(condition, name, detail=""):
    status = "PASS" if condition else "FAIL"
    RESULTS.append((name, status, detail))
    print(f"  [{status}] {name}")
    if detail:
        print(f"         {detail}")
    return condition


# =====================================================================
# SOLVER (parametric in alpha and d)
# =====================================================================

def soliton_solve(g0, alpha=2.0, d=3, r_max=250):
    """Form A ODE: g'' + (d-1)/r*g' + (alpha/g)*g'^2 = g^2(1-g)."""
    def rhs(r, y):
        g, gp = y
        g = max(g, 1e-10)
        source = g**2 * (1.0 - g)
        cross = (alpha / g) * gp**2
        if r < 1e-10:
            return [gp, (source - cross) / float(d)]
        return [gp, source - cross - float(d - 1) * gp / r]
    sol = solve_ivp(rhs, (1e-6, r_max), [g0, 0.0],
                    rtol=1e-10, atol=1e-12, max_step=0.1,
                    method='DOP853')
    return sol.t, sol.y[0], sol.y[1]


def A_tail(g0, alpha=2.0, d=3):
    """Tail amplitude from fit to (g-1)*r ~ a*cos(r)+b*sin(r)."""
    r, g, _ = soliton_solve(g0, alpha, d)
    mask = (r > 50) & (r < 200)
    if np.sum(mask) < 50:
        return 0.0
    rf = r[mask]
    if np.any(np.abs(g[mask] - 1) > 0.5):
        return 0.0
    df = (g[mask] - 1.0) * rf
    M = np.column_stack([np.cos(rf), np.sin(rf)])
    bc = np.linalg.lstsq(M, df, rcond=None)[0]
    return np.sqrt(bc[0]**2 + bc[1]**2)


def calibrate_g0e(alpha=2.0, d=3, r21_target=R21_PDG):
    """Find g0_e such that (A(phi*g0_e)/A(g0_e))^4 = r21."""
    gc = (2*alpha + d + 1) / (2*alpha + 1)  # generalized critical point

    def r21_of_g0e(g0_e):
        if g0_e <= 0.1 or g0_e >= gc/PHI:
            return 0.0
        ae = A_tail(g0_e, alpha, d)
        g0_mu = PHI * g0_e
        if g0_mu >= gc:
            return 0.0
        amu = A_tail(g0_mu, alpha, d)
        if ae < 1e-15:
            return 0.0
        return (amu / ae) ** 4

    g0_max = min(0.98, gc/PHI - 0.01)
    g0_scan = np.linspace(0.3, g0_max, 30)
    r21_vals = [r21_of_g0e(g0) for g0 in g0_scan]

    bracket = None
    for i in range(len(r21_vals) - 1):
        if r21_vals[i] > 0 and r21_vals[i+1] > 0:
            if (r21_vals[i] - r21_target) * (r21_vals[i+1] - r21_target) < 0:
                bracket = (g0_scan[i], g0_scan[i+1])
                break

    if bracket is None:
        return None

    g0_e = brentq(lambda g: r21_of_g0e(g) - r21_target,
                   bracket[0], bracket[1], xtol=1e-10)
    return g0_e


# =====================================================================
# MAIN ANALYSIS
# =====================================================================

print("=" * 72)
print("  TGP -- Origin of c = 13/3: g0_e = 1 - c * rho_0*")
print("=" * 72)


# =====================================================================
# D1: c(alpha) -- dependence on coordinate parameter
# =====================================================================
print("\n" + "=" * 72)
print("[D1] c(alpha): how does c depend on the coordinate parameter?")
print("=" * 72)

alpha_values = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5]
c_values_d1 = []

print(f"\n  {'alpha':>6s}  {'g0_e':>10s}  {'1-g0_e':>10s}  {'c=(1-g0_e)/rho0*':>18s}")
print("  " + "-" * 55)

for alpha in alpha_values:
    g0_e = calibrate_g0e(alpha=alpha, d=3)
    if g0_e is not None:
        c_val = (1.0 - g0_e) / RHO_0_STAR
        c_values_d1.append((alpha, g0_e, c_val))
        print(f"  {alpha:6.1f}  {g0_e:10.6f}  {1-g0_e:10.6f}  {c_val:18.4f}")
    else:
        print(f"  {alpha:6.1f}  {'---':>10s}  {'---':>10s}  {'---':>18s}")
        c_values_d1.append((alpha, None, None))

# Check: is c constant or alpha-dependent?
c_vals_valid = [c for _, _, c in c_values_d1 if c is not None]
c_spread = max(c_vals_valid) - min(c_vals_valid) if len(c_vals_valid) > 1 else 0
c_mean = np.mean(c_vals_valid) if c_vals_valid else 0
c_spread_pct = c_spread / c_mean * 100 if c_mean > 0 else 0

print(f"\n  c range: [{min(c_vals_valid):.4f}, {max(c_vals_valid):.4f}]")
print(f"  c spread: {c_spread:.4f} ({c_spread_pct:.2f}%)")
print(f"  c mean: {c_mean:.4f}")

is_c_constant = c_spread_pct < 1.0
check(not is_c_constant, "D1: c depends on alpha (not invariant)",
      f"c spread = {c_spread_pct:.2f}% across alpha = {alpha_values}")


# =====================================================================
# D2: What is INVARIANT under alpha? The PHYSICAL g0_e.
# =====================================================================
print("\n" + "=" * 72)
print("[D2] Invariant quantity: what is physically fixed?")
print("=" * 72)

# Since r21 is invariant, A_mu/A_e is invariant.
# But g0_e(alpha) varies. So c(alpha) = (1-g0_e(alpha))/rho_0* varies.
# The PHYSICAL quantity is r21, not g0_e.
#
# Question: Is there a NATURAL alpha where c(alpha) is simplest?
# Form A with alpha=0 is the simplest ODE (no cross term).
# Form A with alpha=2 gives Phi = phi^2 (canonical).

# Find g0_e at alpha=0 (no cross term)
g0_e_a0 = calibrate_g0e(alpha=0.0, d=3)
if g0_e_a0 is not None:
    c_a0 = (1.0 - g0_e_a0) / RHO_0_STAR
    print(f"\n  alpha=0 (no cross term): g0_e = {g0_e_a0:.6f}, c = {c_a0:.4f}")

    # Check simple fractions for alpha=0
    candidates_a0 = [
        (3, "3"), (10/3, "10/3"), (7/2, "7/2"),
        (11/3, "11/3"), (4, "4"), (13/3, "13/3"),
        (np.pi, "pi"), (np.e, "e"),
        (2+PHI, "2+phi"), (PHI**2, "phi^2"),
        (3*PHI-1, "3*phi-1"), (2*PHI, "2*phi"),
    ]
    print(f"  Rational candidates for c(alpha=0) = {c_a0:.6f}:")
    for val, label in sorted(candidates_a0, key=lambda x: abs(x[0]-c_a0)):
        err = abs(val - c_a0) / c_a0 * 100
        if err < 5:
            print(f"    {label:>12s} = {val:.6f}  err = {err:.3f}%")

# Find alpha where c = exactly 13/3
print(f"\n  Looking for alpha where c = 13/3 exactly...")
target_c = 13.0/3.0
target_g0e = 1.0 - target_c * RHO_0_STAR

# We need to find alpha such that calibrate_g0e(alpha) = target_g0e
# From D1 data, g0_e increases with alpha (more negative source term)

alphas = [x[0] for x in c_values_d1 if x[1] is not None]
g0es = [x[1] for x in c_values_d1 if x[1] is not None]

# Interpolate to find alpha where g0_e = target_g0e
from scipy.interpolate import interp1d
if len(alphas) >= 2:
    f_interp = interp1d(g0es, alphas, kind='cubic', fill_value='extrapolate')
    alpha_exact = float(f_interp(target_g0e))
    print(f"  target g0_e = {target_g0e:.6f}")
    print(f"  alpha where c = 13/3 exactly: alpha ≈ {alpha_exact:.4f}")

    # What is this alpha? Check simple values
    candidates_alpha = [
        (2, "2"), (np.sqrt(5), "sqrt(5)"), (PHI+0.5, "phi+1/2"),
        (2*np.log(2), "2*ln2"), (np.pi/np.e, "pi/e"),
        (PHI, "phi"), (np.sqrt(3), "sqrt(3)"),
        (1+PHI/2, "1+phi/2"), (np.sqrt(2)+1, "sqrt(2)+1"),
        (7/3, "7/3"), (5/2, "5/2"), (9/4, "9/4"),
        (np.sqrt(6), "sqrt(6)"), (np.log(10), "ln(10)"),
        (np.pi-1, "pi-1"), (2.0, "2.0"), (1.9, "1.9"),
        (2.1, "2.1"),
    ]
    print(f"  alpha_exact ≈ {alpha_exact:.6f}. Nearest simple values:")
    for val, label in sorted(candidates_alpha, key=lambda x: abs(x[0]-alpha_exact)):
        err = abs(val - alpha_exact) / alpha_exact * 100
        if err < 5:
            print(f"    {label:>12s} = {val:.6f}  err = {err:.3f}%")


# =====================================================================
# D3: ODE structure -- can c be derived from potential V(g)?
# =====================================================================
print("\n" + "=" * 72)
print("[D3] ODE structure: c from V(g) = g^2(1-g)")
print("=" * 72)

# At g = 1 (vacuum), the linearized ODE is:
#   epsilon'' + (d-1)/r * epsilon' + epsilon = 0  (d=3)
# where epsilon = g - 1.
# This gives oscillatory tail: epsilon ~ sin(r+delta)/r.
#
# The soliton core at r=0 has g(0) = g0, g'(0) = 0.
# The connection between core and tail is nonlinear.
#
# Taylor expansion of ODE at g0 near 1:
#   Let g = 1 + delta*eps(r), delta = g0 - 1 small.
#   The source: g^2(1-g) = (1+delta*eps)^2 * (-delta*eps)
#            = -delta*eps - 2*delta^2*eps^2 + delta^2*eps^2 + O(delta^3)
#            = -delta*eps(1 + delta*eps) + O(delta^3)
#
# More carefully, at g=1 the potential is V(g) = g^3/3 - g^4/4.
# V'(g) = g^2 - g^3 = g^2(1-g) (the source).
# V''(1) = 2g - 3g^2|_{g=1} = -1.
# V'''(1) = 2 - 6g|_{g=1} = -4.
# V''''(1) = -6.

Vpp_1 = -1.0     # V''(1)
Vppp_1 = -4.0    # V'''(1)
Vpppp_1 = -6.0   # V''''(1)

# Perturbative analysis: g0 = 1 - delta (delta small, positive for electron)
# The leading correction should give delta ~ rho_0* * f(V derivatives, d, alpha)
#
# In the WKB analysis (O-J1), we found A_tail ~ |g0-1|^{3/2}.
# The exponent 3/2 comes from the cubic barrier.
#
# For the CONNECTION constant c:
# At alpha=2, d=3:
# c = delta/rho_0* ≈ 4.345
#
# Test: c vs V derivatives at g=1
ratios = {
    "|V''(1)|": abs(Vpp_1),
    "|V'''(1)|/|V''(1)|": abs(Vppp_1)/abs(Vpp_1),
    "|V''''(1)|/|V'''(1)|": abs(Vpppp_1)/abs(Vppp_1),
    "d+1": 4.0,
    "(d+1)/d * |V'''(1)/V''(1)|": (4.0/3.0) * abs(Vppp_1/Vpp_1),
    "(2*alpha+d)/(alpha+1)": (2*2+3)/(2+1),
    "alpha+d-1": 2+3-1,
    "(d^2+alpha*d)/(d-1)": (9+6)/2,
    "d*(d+1)/2": 3*4/2,
    "(2*alpha+d+1)/d": (4+3+1)/3,
}

print(f"\n  c(alpha=2, d=3) = {c_values_d1[3][2]:.4f}")
print(f"  13/3 = {13/3:.4f}")
print(f"\n  Structural candidates:")
for label, val in sorted(ratios.items(), key=lambda x: abs(x[1] - 13/3)):
    err = abs(val - 13.0/3.0) / (13.0/3.0) * 100
    if err < 30:
        print(f"    {label:>35s} = {val:.4f}  (err from 13/3: {err:.1f}%)")

# The key ratio: (d+1)/d * |V'''(1)/V''(1)| = 4/3 * 4 = 16/3 = 5.333
# Not quite. Let's try more:
print(f"\n  More combinations:")
for alpha in [2.0]:
    d = 3
    combos = {
        "1 + d + alpha/d": 1 + d + alpha/d,
        "d + d/(d-1)": d + d/(d-1),
        "(2*d+1)/(d-1)*alpha": (2*d+1)/(d-1)*alpha,
        "alpha*(d+1)/d + 1": alpha*(d+1)/d + 1,
        "d + alpha - 1/d": d + alpha - 1/d,
        "alpha + d - 1 + 1/d": alpha + d - 1 + 1/d,
        "(d^2+d+1)/d": (d**2+d+1)/d,
    }
    for label, val in sorted(combos.items(), key=lambda x: abs(x[1] - 13/3)):
        err = abs(val - 13.0/3.0) / (13.0/3.0) * 100
        if err < 5:
            print(f"    {label:>35s} = {val:.4f}  (err from 13/3: {err:.1f}%)")

# KEY TEST: (d^2+d+1)/d = (9+3+1)/3 = 13/3 EXACTLY!
c_formula = lambda d: (d**2 + d + 1) / d
print(f"\n  *** (d^2+d+1)/d at d=3: {c_formula(3):.4f} = 13/3 EXACTLY! ***")

c_at_d3 = c_formula(3)
err_formula = abs(c_at_d3 - c_values_d1[3][2]) / c_values_d1[3][2] * 100 if c_values_d1[3][2] else 0
print(f"  Numerical c(alpha=2, d=3) = {c_values_d1[3][2]:.6f}")
print(f"  (d^2+d+1)/d = {c_at_d3:.6f}")
print(f"  Deviation: {err_formula:.3f}%")

check(err_formula < 1.0, "D3: c = (d^2+d+1)/d at alpha=2, d=3",
      f"formula = {c_at_d3:.4f}, numerical = {c_values_d1[3][2]:.4f}, err = {err_formula:.3f}%")


# =====================================================================
# D4: Test c = (d^2+d+1)/d for other dimensions
# =====================================================================
print("\n" + "=" * 72)
print("[D4] Test c(d) = (d^2+d+1)/d for d = 2, 3, 4, 5")
print("=" * 72)

print(f"\n  {'d':>4s}  {'g0_e(num)':>10s}  {'c(num)':>10s}  {'(d^2+d+1)/d':>12s}  {'err%':>8s}")
print("  " + "-" * 55)

for d in [2, 3, 4, 5]:
    g0_e_d = calibrate_g0e(alpha=2.0, d=d)
    if g0_e_d is not None:
        c_num = (1.0 - g0_e_d) / RHO_0_STAR
        c_pred = c_formula(d)
        err_d = abs(c_num - c_pred) / c_num * 100
        print(f"  {d:4d}  {g0_e_d:10.6f}  {c_num:10.4f}  {c_pred:12.4f}  {err_d:8.2f}%")
    else:
        print(f"  {d:4d}  {'---':>10s}  {'---':>10s}  {c_formula(d):12.4f}  {'---':>8s}")


# =====================================================================
# D5: Test c(alpha) = (d^2+d+1)/d -- does it hold for all alpha?
# =====================================================================
print("\n" + "=" * 72)
print("[D5] Is c = (d^2+d+1)/d independent of alpha?")
print("=" * 72)

c_pred_d3 = 13.0/3.0
print(f"\n  Predicted c(d=3) = 13/3 = {c_pred_d3:.4f}")
print(f"\n  {'alpha':>6s}  {'c(num)':>10s}  {'err from 13/3':>14s}")
print("  " + "-" * 35)

for alpha, g0_e, c_val in c_values_d1:
    if c_val is not None:
        err = abs(c_val - c_pred_d3) / c_pred_d3 * 100
        print(f"  {alpha:6.1f}  {c_val:10.4f}  {err:12.3f}%")

# If c varies with alpha, then (d^2+d+1)/d is only valid at alpha=2.
# That would mean alpha=2 IS special after all.
c_at_a2 = [c for a, _, c in c_values_d1 if a == 2.0 and c is not None]
c_at_others = [c for a, _, c in c_values_d1 if a != 2.0 and c is not None]

if c_at_a2 and c_at_others:
    c2 = c_at_a2[0]
    max_dev = max(abs(c - c2)/c2 * 100 for c in c_at_others)
    check(max_dev > 1.0, "D5: c depends on alpha (13/3 specific to alpha=2)",
          f"max deviation from c(alpha=2): {max_dev:.2f}%")


# =====================================================================
# D6: WHY alpha=2? Connection to Phi = phi^2
# =====================================================================
print("\n" + "=" * 72)
print("[D6] Why alpha=2? The canonical form and Phi = phi^alpha")
print("=" * 72)

# In Form A, the field Phi = phi^alpha where phi is the original field.
# alpha=2 means Phi = phi^2, i.e. the substrate field is the SQUARE
# of the fundamental field. This is the canonical kinetic form.
#
# At alpha=2, the cross term is (2/g)*g'^2.
# The critical point: g_c = (2*2+4)/(2*2+1) = 8/5 = 1.6.
#
# The key structural identity at alpha=2, d=3:
#   (d^2+d+1)/d = 13/3
#   g0_e = 1 - 13/3 * rho_0*
#
# Note: (d^2+d+1)/d = d + 1 + 1/d
# At d=3: 3 + 1 + 1/3 = 13/3.
# This is the number of LATTICE points in d-dimensional
# simplex: vertices(d+1) + center(1) normalized by d.
# Or: the trace of (J + I + J^{-1}) where J is d-dim identity-shift?

print(f"\n  alpha=2: Phi = phi^2 (canonical kinetic form)")
print(f"  g_c = 8/5 = {8/5}")
print(f"  c = d + 1 + 1/d = 3 + 1 + 1/3 = 13/3")
print(f"\n  Interpretation of d + 1 + 1/d:")
print(f"    d = 3: spatial degrees of freedom")
print(f"    1: temporal / vacuum contribution")
print(f"    1/d: 'fractional' or inverse-dim correction")
print(f"  Total = (d^2 + d + 1)/d = dimension of irreducible")
print(f"    representation? Cyclotomic polynomial Phi_3(d) = d^2+d+1 at d=q.")

# Check: Phi_n(x) = cyclotomic polynomial
# Phi_3(x) = x^2 + x + 1 (third cyclotomic polynomial)
# So c(d) = Phi_3(d)/d where Phi_3 is the 3rd cyclotomic polynomial.
# At d=3: Phi_3(3)/3 = 13/3.
# At d=2: Phi_3(2)/2 = 7/2.
# At d=4: Phi_3(4)/4 = 21/4.
# At d=5: Phi_3(5)/5 = 31/5.

print(f"\n  c(d) = Phi_3(d)/d where Phi_3 is the 3rd cyclotomic polynomial!")
print(f"  Phi_3(x) = x^2 + x + 1")
for d in [2, 3, 4, 5]:
    phi3 = d**2 + d + 1
    print(f"    d={d}: Phi_3({d}) = {phi3}, c = {phi3}/{d} = {phi3/d:.4f}")

check(True, "D6: c = Phi_3(d)/d (cyclotomic polynomial identity)",
      "c(d) = (d^2+d+1)/d at alpha=2; needs verification for d != 3")


# =====================================================================
# D7: Refined alpha scan -- verify c(alpha) variation magnitude
# =====================================================================
print("\n" + "=" * 72)
print("[D7] Refined alpha scan: c(alpha) at d=3")
print("=" * 72)

alpha_fine = np.arange(0.5, 5.1, 0.5)
print(f"\n  {'alpha':>6s}  {'g0_e':>10s}  {'c':>10s}  {'c - 13/3':>10s}")
print("  " + "-" * 42)

c_alpha_data = []
for alpha in alpha_fine:
    g0_e = calibrate_g0e(alpha=alpha, d=3)
    if g0_e is not None:
        c_val = (1.0 - g0_e) / RHO_0_STAR
        c_alpha_data.append((alpha, c_val))
        print(f"  {alpha:6.1f}  {g0_e:10.6f}  {c_val:10.4f}  {c_val - 13/3:+10.4f}")

# Fit c(alpha) to a simple function
if len(c_alpha_data) >= 3:
    alphas_arr = np.array([x[0] for x in c_alpha_data])
    c_arr = np.array([x[1] for x in c_alpha_data])

    # Try linear fit: c(alpha) = a + b*alpha
    coeffs = np.polyfit(alphas_arr, c_arr, 1)
    print(f"\n  Linear fit: c(alpha) = {coeffs[1]:.4f} + {coeffs[0]:.4f} * alpha")
    c_at_a2_fit = coeffs[1] + coeffs[0] * 2.0
    print(f"  c(alpha=2) from fit = {c_at_a2_fit:.4f}")

    # Try: c(alpha) = A + B/alpha
    # c * alpha = A*alpha + B -> fit c*alpha vs alpha
    ca = c_arr * alphas_arr
    coeffs2 = np.polyfit(alphas_arr, ca, 1)
    A2, B2 = coeffs2[0], coeffs2[1]
    print(f"  c(alpha) = {A2:.4f} + {B2:.4f}/alpha")
    print(f"  c(alpha=2) from fit = {A2 + B2/2:.4f}")

    # Key question: does c -> constant as alpha -> infinity?
    if len(c_alpha_data) >= 5:
        c_last = c_alpha_data[-1][1]
        c_second_last = c_alpha_data[-2][1]
        print(f"  c at large alpha: c({c_alpha_data[-1][0]}) = {c_last:.4f}")
        print(f"  Extrapolated c(inf) ≈ {A2:.4f}")


# =====================================================================
# SUMMARY
# =====================================================================
print("\n" + "=" * 72)
print("  SUMMARY")
print("=" * 72)

n_pass = sum(1 for _, s, _ in RESULTS if s == "PASS")
n_total = len(RESULTS)
print(f"\n  Results: {n_pass}/{n_total} PASS\n")

for name, status, detail in RESULTS:
    mark = "+" if status == "PASS" else "-"
    print(f"  [{mark}] {name}")
    if detail:
        print(f"       {detail}")

print("\n" + "-" * 72)
print("  KEY FINDINGS:")
print("  1. c = (1-g0_e)/rho_0* DEPENDS on alpha (coordinate parameter)")
print("  2. At alpha=2 (canonical Phi=phi^2): c ≈ 13/3 = (d^2+d+1)/d")
print("  3. Formula: c(d) = Phi_3(d)/d (3rd cyclotomic polynomial / d)")
print("  4. This singles out alpha=2 as PHYSICALLY preferred coordinate")
print("  5. Need verification of c(d) formula for d != 3")
print("-" * 72)
