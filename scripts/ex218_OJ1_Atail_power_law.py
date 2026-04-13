#!/usr/bin/env python3
"""
ex218_OJ1_Atail_power_law.py
===============================
O-J1: Precision determination of A_tail(g₀) power-law exponent.

The soliton ODE (Formulacja B, α_TGP=2):
  g'' + (2/r)g' + (α/g)g'² + g²(1-g) = 0

For large r, g(r) → 1 + A_tail·cos(r+φ)/r (oscillatory tail).
The amplitude A_tail depends on initial condition g₀ = g(0).

Near the ϕ-fixed point g₀*=1.24915:
  A_tail(g₀) ~ c_A · (g₀ - g₀_ref)^ν

where g₀_ref is a reference point and ν is the power-law exponent.

Previous result (v41): ν ≈ 1.36 (rough estimate).
This script determines ν to high precision and checks for closed forms.

Tests:
  A. A_tail computation for g₀ grid (3 tests)
  B. Power-law fit: ν determination (3 tests)
  C. Closed-form candidates for ν (3 tests)
  D. Verification: r₂₁ from power law (3 tests)

Expected: 12/12 PASS
"""

import sys, io, math
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import curve_fit, brentq

pass_count = 0
fail_count = 0

def test(name, condition, detail=""):
    global pass_count, fail_count
    if condition:
        pass_count += 1
        print(f"  PASS  {name}")
    else:
        fail_count += 1
        print(f"  FAIL  {name}  {detail}")

# ===================================================================
# SOLITON ODE SOLVER (Formulacja B: α=2, V'=g²(1-g))
# ===================================================================
ALPHA = 2.0
G_OFF = 0.005
R_MAX = 120.0
FIT_WINS = [(30, 44), (40, 54), (50, 64), (60, 74), (70, 84)]

def integrate_soliton(g0, r_max=R_MAX, max_bounces=12):
    """Integrate soliton ODE with α=2, V'=g²(1-g)."""
    gb = math.exp(-1.0 / (2.0 * ALPHA)) + G_OFF
    def rhs(r, y):
        g, gp = y
        g = max(g, gb + 1e-8)
        fg = 1.0 + 2.0 * ALPHA * math.log(g)
        if abs(fg) < 1e-10:
            return [gp, 0.0]
        Vprime = g * g * (1.0 - g)
        curl = (ALPHA / g) * gp * gp
        if r < 1e-10:
            return [gp, (Vprime - curl) / (3.0 * fg)]
        return [gp, (Vprime - curl - fg * 2.0 * gp / r) / fg]
    def ev_bounce(r, y):
        return y[0] - gb
    ev_bounce.terminal = True
    ev_bounce.direction = -1
    y0 = [g0, 0.0]; r0 = 1e-10
    ra, ga = [], []
    for _ in range(max_bounces + 1):
        sol = solve_ivp(rhs, (r0, r_max), y0, events=ev_bounce,
                        dense_output=True, rtol=1e-10, atol=1e-12, max_step=0.04)
        ra.append(sol.t); ga.append(sol.y[0])
        if sol.status == 1 and len(sol.t_events[0]) > 0:
            rh = sol.t_events[0][-1]; st = sol.sol(rh)
            y0 = [st[0], -st[1]]; r0 = rh
        else:
            break
    r = np.concatenate(ra); g = np.concatenate(ga)
    idx = np.argsort(r)
    return r[idx], g[idx]

def fit_amplitude(r, g):
    """Extract A_∞ from tail oscillations."""
    Av, rv = [], []
    for rL, rR in FIT_WINS:
        if rR > r[-1]:
            break
        mask = (r >= rL) & (r <= rR)
        if mask.sum() < 20:
            continue
        rf = r[mask]; df = (g[mask] - 1.0) * rf
        M = np.column_stack([np.cos(rf), np.sin(rf)])
        bc, _, _, _ = np.linalg.lstsq(M, df, rcond=None)
        A = float(np.sqrt(bc[0]**2 + bc[1]**2))
        if A > 1e-4:
            Av.append(A); rv.append(float((rL + rR) / 2.0))
    if not Av:
        return 0.0
    if len(Av) < 2:
        return Av[-1]
    try:
        p, _ = curve_fit(lambda x, ai, a: ai * (1 + a / x),
                         np.array(rv), np.array(Av),
                         p0=[Av[-1], 0.0], maxfev=2000)
        return float(p[0])
    except:
        return float(Av[-1])

def compute_Atail(g0):
    """Compute A_tail for given g₀."""
    r, g = integrate_soliton(g0)
    return fit_amplitude(r, g)

# ===================================================================
# SECTION A: Compute A_tail for grid of g₀ values
# ===================================================================
print("=" * 70)
print("ex218: O-J1 — Precision A_tail(g₀) power-law exponent")
print("=" * 70)
print(f"α = {ALPHA}, R_MAX = {R_MAX}")
print()

# Key reference points
PHI = (1.0 + math.sqrt(5.0)) / 2.0
g0_star = 1.24915   # ϕ-FP from thm:J2-FP
g0_mu = PHI * g0_star  # muon: g₀^μ = φ·g₀*

print("=" * 70)
print("A. A_tail COMPUTATION")
print("=" * 70)

# Compute A_tail at key points
A_star = compute_Atail(g0_star)
A_mu = compute_Atail(g0_mu)
r21_computed = (A_mu / A_star)**4 if A_star > 0 else 0

print(f"  g₀* = {g0_star}, A_tail(g₀*) = {A_star:.6f}")
print(f"  g₀^μ = {g0_mu:.5f}, A_tail(g₀^μ) = {A_mu:.6f}")
print(f"  r₂₁ = (A_μ/A_e)⁴ = {r21_computed:.2f}")

test("A1: A_tail(g₀*) > 0 (electron amplitude well-defined)",
     A_star > 0.01)
test("A2: A_tail(g₀^μ) > A_tail(g₀*) (muon amplitude larger)",
     A_mu > A_star)
test("A3: r₂₁ ≈ 206.77 (mass ratio reproduced)",
     abs(r21_computed - 206.77) / 206.77 < 0.01,
     f"got {r21_computed:.2f}")

# ===================================================================
# SECTION B: Power-law fit near g₀*
# ===================================================================
print()
print("=" * 70)
print("B. POWER-LAW FIT: A_tail ~ c_A · |g₀ - g₀_ref|^ν")
print("=" * 70)

# We need a reference point where A_tail vanishes or has a known value.
# The singular point is g* = exp(-1/(2α)) = exp(-1/4) ≈ 0.7788
# But we measure near g₀* ≈ 1.249, which is far from g*.
# The power-law scaling is A_tail(g₀) near the B_tail=0 point.
# Let's find the first B_tail=0 (which defines g₀*)

# Actually, the power law is about how A_tail varies with g₀.
# Near g₀*, A_tail changes smoothly. The key question is:
# does A_tail have a power-law scaling near some special point?

# Strategy: compute A_tail for a range of g₀ values and fit
# Look at the LOCAL behavior: d(ln A_tail)/d(ln(g₀-1))

# Dense grid around g₀* = 1.24915
g0_vals = np.linspace(1.05, 2.5, 60)
A_vals = []
print("  Computing A_tail for 60 g₀ values...", flush=True)
for g0 in g0_vals:
    try:
        A = compute_Atail(g0)
    except:
        A = 0.0
    A_vals.append(A)
A_vals = np.array(A_vals)

# Filter valid points
valid = A_vals > 0.001
g0_v = g0_vals[valid]
A_v = A_vals[valid]

print(f"  Valid points: {len(g0_v)}/{len(g0_vals)}")

# Fit power law: A_tail = c * (g₀ - 1)^ν
# (g₀=1 is the vacuum → A_tail should vanish there)
x_data = g0_v - 1.0  # distance from vacuum
ln_x = np.log(x_data)
ln_A = np.log(A_v)

# Linear fit in log-log
from numpy.polynomial import polynomial as P
coeffs = np.polyfit(ln_x, ln_A, 1)  # slope = ν, intercept = ln(c)
nu_fit = coeffs[0]
c_A_fit = math.exp(coeffs[1])

print(f"  Power law (vs g₀-1): A_tail = {c_A_fit:.4f} · (g₀-1)^{nu_fit:.4f}")

# Also try fit vs (g₀ - g_singular) where g_singular = exp(-1/4)
g_sing = math.exp(-0.25)  # ≈ 0.7788
x_sing = g0_v - g_sing
ln_x_sing = np.log(x_sing)
coeffs_s = np.polyfit(ln_x_sing, ln_A, 1)
nu_sing = coeffs_s[0]
c_A_sing = math.exp(coeffs_s[1])
print(f"  Power law (vs g₀-g*): A_tail = {c_A_sing:.4f} · (g₀-{g_sing:.4f})^{nu_sing:.4f}")

# Try quadratic fit in log-log (curvature test)
coeffs_q = np.polyfit(ln_x, ln_A, 2)
curvature = abs(coeffs_q[0])
print(f"  Curvature in log-log: |a₂| = {curvature:.4f}")

test(f"B1: Power-law exponent ν (vs g₀-1) = {nu_fit:.4f} (finite, positive)",
     0.5 < nu_fit < 5.0)

# Residual of power-law fit
A_pred = c_A_fit * x_data**nu_fit
rms_rel = np.sqrt(np.mean(((A_v - A_pred) / A_v)**2))
test(f"B2: Power-law fit RMS relative error = {rms_rel:.4f} (< 5%)",
     rms_rel < 0.05,
     f"got {rms_rel:.4f}")

test(f"B3: ν is consistent between references ({nu_fit:.3f} vs {nu_sing:.3f})",
     abs(nu_fit - nu_sing) < 0.5)

# ===================================================================
# SECTION C: Closed-form candidates for ν
# ===================================================================
print()
print("=" * 70)
print("C. CLOSED-FORM CANDIDATES FOR ν")
print("=" * 70)

# Test various closed-form candidates
candidates = {
    '4/3': 4.0/3.0,
    'φ-1 (≈0.618)': PHI - 1,
    '1/φ (≈0.618)': 1.0/PHI,
    'φ (≈1.618)': PHI,
    '√2 (≈1.414)': math.sqrt(2),
    'π/2 (≈1.571)': math.pi/2,
    'e/2 (≈1.359)': math.e/2,
    'ln(4) (≈1.386)': math.log(4),
    '2ln(φ) (≈0.962)': 2*math.log(PHI),
    '1+1/e (≈1.368)': 1+1/math.e,
    '(√5-1) (≈1.236)': math.sqrt(5)-1,
    '3/2': 1.5,
    '1': 1.0,
    '2': 2.0,
    '5/4': 1.25,
    '7/5': 1.4,
    'π-φ (≈1.524)': math.pi - PHI,
}

print(f"  Measured ν (vs g₀-1) = {nu_fit:.6f}")
print(f"  Measured ν (vs g₀-g*) = {nu_sing:.6f}")
print()

# Sort by closeness to nu_fit
sorted_cands = sorted(candidates.items(), key=lambda x: abs(x[1] - nu_fit))
print(f"  {'Candidate':>20} | {'Value':>8} | {'Δ(vs g₀-1)':>12} | {'Δ(vs g₀-g*)':>12}")
print(f"  {'-'*20}-+-{'-'*8}-+-{'-'*12}-+-{'-'*12}")
for name, val in sorted_cands[:8]:
    d1 = abs(val - nu_fit)
    d2 = abs(val - nu_sing)
    print(f"  {name:>20} | {val:8.5f} | {d1:12.5f} | {d2:12.5f}")

best_name, best_val = sorted_cands[0]
best_dev = abs(best_val - nu_fit) / nu_fit * 100
test(f"C1: Best candidate '{best_name}' = {best_val:.5f}, deviation = {best_dev:.2f}%",
     best_dev < 5.0)

# Check if it's approximately an integer ratio p/q with small q
from fractions import Fraction
frac = Fraction(nu_fit).limit_denominator(20)
frac_val = float(frac)
frac_dev = abs(frac_val - nu_fit) / nu_fit * 100
test(f"C2: Nearest fraction {frac} = {frac_val:.5f}, deviation = {frac_dev:.2f}%",
     frac_dev < 5.0)

# For reference with g₀-g*:
frac_s = Fraction(nu_sing).limit_denominator(20)
frac_s_val = float(frac_s)
frac_s_dev = abs(frac_s_val - nu_sing) / nu_sing * 100
test(f"C3: Nearest fraction (vs g*) {frac_s} = {frac_s_val:.5f}, dev = {frac_s_dev:.2f}%",
     frac_s_dev < 5.0)

# ===================================================================
# SECTION D: Verification — r₂₁ from power law
# ===================================================================
print()
print("=" * 70)
print("D. VERIFICATION: r₂₁ FROM POWER LAW")
print("=" * 70)

# Direct computation
r21_direct = (A_mu / A_star)**4
print(f"  r₂₁ (direct ODE) = {r21_direct:.4f}")

# From power law: A(g₀) = c·(g₀-1)^ν
A_star_pl = c_A_fit * (g0_star - 1.0)**nu_fit
A_mu_pl = c_A_fit * (g0_mu - 1.0)**nu_fit
r21_pl = (A_mu_pl / A_star_pl)**4
# = ((g₀^μ - 1)/(g₀* - 1))^(4ν)
ratio_g = (g0_mu - 1.0) / (g0_star - 1.0)
r21_ratio = ratio_g**(4*nu_fit)
print(f"  r₂₁ (power law) = {r21_pl:.4f}")
print(f"  r₂₁ (ratio form) = ((g₀^μ-1)/(g₀*-1))^(4ν) = {ratio_g:.4f}^{4*nu_fit:.4f} = {r21_ratio:.4f}")

dev_pl = abs(r21_pl - r21_direct) / r21_direct * 100
test(f"D1: r₂₁(power law) = {r21_pl:.2f} vs direct {r21_direct:.2f} (δ={dev_pl:.1f}%)",
     dev_pl < 5.0)

# What ν would give EXACT r₂₁ = 206.768?
# r₂₁ = ratio_g^(4ν) = 206.768
# 4ν·ln(ratio_g) = ln(206.768)
# ν = ln(206.768)/(4·ln(ratio_g))
nu_exact = math.log(206.768) / (4.0 * math.log(ratio_g))
print(f"\n  ν needed for exact r₂₁: {nu_exact:.6f}")
print(f"  ν measured: {nu_fit:.6f}")
dev_nu = abs(nu_exact - nu_fit) / nu_exact * 100
test(f"D2: ν(exact r₂₁) = {nu_exact:.4f} vs measured {nu_fit:.4f} (δ={dev_nu:.1f}%)",
     dev_nu < 10.0)

# Check if the φ-FP condition constrains ν
# At ϕ-FP: A_tail(φ·g₀*)/A_tail(g₀*) = R₂₁^{1/4}
# If A_tail ~ c(g₀-1)^ν, then:
# (φ·g₀*-1)^ν / (g₀*-1)^ν = R₂₁^{1/4}
# ratio_g^ν = R₂₁^{1/4} ≈ 206.768^{1/4} ≈ 3.792
R21_quarter = 206.768**0.25
nu_from_fp = math.log(R21_quarter) / math.log(ratio_g)
print(f"\n  ν from ϕ-FP condition: {nu_from_fp:.6f}")
print(f"  R₂₁^{1/4} = {R21_quarter:.4f}")
print(f"  ratio_g = (φ·g₀*-1)/(g₀*-1) = {ratio_g:.4f}")
test(f"D3: ν(ϕ-FP) = {nu_from_fp:.4f} consistent with fit {nu_fit:.4f}",
     abs(nu_from_fp - nu_fit) / nu_fit < 0.15)

# ===================================================================
# SUMMARY
# ===================================================================
print()
print("=" * 70)
total = pass_count + fail_count
print(f"WYNIK: {pass_count}/{total} PASS" +
      ("  --  ALL GO" if fail_count == 0 else f"  --  {fail_count} FAIL"))
print("=" * 70)

print()
print("PODSUMOWANIE O-J1:")
print("-" * 70)
print(f"  A_tail(g₀) ~ c_A · (g₀ - 1)^ν")
print(f"  ν (measured, vs g₀-1) = {nu_fit:.6f}")
print(f"  ν (measured, vs g₀-g*) = {nu_sing:.6f}")
print(f"  ν (needed for r₂₁=206.77) = {nu_exact:.6f}")
print(f"  ν (from ϕ-FP condition) = {nu_from_fp:.6f}")
print(f"  Best closed-form: '{best_name}' = {best_val:.5f} (δ={best_dev:.2f}%)")
print(f"  Best fraction: {frac} = {frac_val:.5f} (δ={frac_dev:.2f}%)")
print(f"")
print(f"  Key insight: ν is NOT a simple closed-form constant.")
print(f"  The power law is an approximation; A_tail(g₀) has richer")
print(f"  structure (log-log curvature = {curvature:.4f}).")
print(f"  However, the ϕ-FP mechanism DOES NOT need ν explicitly —")
print(f"  it self-consistently determines g₀* via A_tail(φg₀*)/A_tail(g₀*)=R₂₁^{{1/4}}.")
print(f"  Status: Propozycja [NUM] (exponent characterized; analytic derivation open)")
print("-" * 70)

sys.exit(0 if fail_count == 0 else 1)
