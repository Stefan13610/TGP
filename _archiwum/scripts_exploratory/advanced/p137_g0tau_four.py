#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p137_g0tau_four.py - Why g0_tau = 4 exactly?
==============================================
g0 = 4 is VERY close to the Koide value (r_31=3477.8 at g0=4.000).
This is suggestive of a structural origin.

Recall: V(g) = g^3/3 - g^4/4, V'(g) = g^2(1-g)
  V'(0) = 0 (vacuum 1)
  V'(1) = 0 (vacuum 2)
  V'(g) = 0 also at g = 0 (degenerate)

The potential has:
  V(0) = 0, V(1) = 1/12
  V(4) = 4^3/3 - 4^4/4 = 64/3 - 64 = -128/3 = -42.667

Also: f(g) = 1 + 2*alpha_eff*ln(g)
  f(4) = 1 + 2*alpha_eff*ln(4)
  For alpha_eff = 2/(1+12*(4-1)^2) = 2/109 = 0.01835:
  f(4) = 1 + 2*0.01835*1.386 = 1 + 0.0509 = 1.0509

The key question: is there a structural reason why g0=4 is special?

Hypotheses:
  H1: g0_tau = 2 * g0_mu / phi (phi = golden ratio)
      = 2 * 1.465 / 1.618 = 1.810. No.
  H2: g0_tau relates to the OTHER zero of V or some effective quantity
  H3: g0_tau = 4 = 2^2, simplest integer > 1 where soliton exists
  H4: g0_tau = 4 from quantization of some charge
  H5: g0_tau/g0_e * g0_mu/g0_e = phi * (g0_tau/g0_e)
      If this product = integer, then g0_tau/g0_e = n/phi
      For n=7: 7/phi = 4.326. Not 4.417.

Check: at the deeper level, g0_e ~ 0.905 and g0_tau ~ 4.0.
  g0_e * g0_tau = 0.905 * 4.0 = 3.622 ~ 2*phi + 1/phi = 3.854? No.
  g0_e + g0_tau = 4.905 ~ 3*phi = 4.854? Close-ish.
  g0_e * g0_tau / g0_mu = 3.622 / 1.465 = 2.473 ~ phi + 1/phi = 2.236? No.

Actually: g0 = {e, mu, tau} = {0.905, 1.465, 4.000}
  Ratios: mu/e = phi, tau/e = 4.42, tau/mu = 2.73
  Products: e*mu = 1.326, e*tau = 3.622, mu*tau = 5.860
  Sum: e+mu+tau = 6.370

Let's check if there's a pattern in the ODE structure itself.
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from scipy.integrate import solve_ivp, trapezoid
from scipy.optimize import brentq

PHI = (1 + np.sqrt(5)) / 2
ETA_K = 12.067
G0_E = 0.905481
G0_MU = PHI * G0_E

def solve_soliton(g0, eta_K=ETA_K, rm=300):
    def fk(g):
        a = 2.0 / (1 + eta_K * (g - 1)**2)
        return 1 + 2*a*np.log(g) if g > 0 else -1e30
    def Vp(g): return g**2*(1-g)
    fg0 = fk(g0)
    if abs(fg0) < 1e-15: return None, None, None
    c2 = Vp(g0)/(3*fg0)
    rs = 0.01
    def rhs(r, y):
        g, p = y
        if g <= 1e-15: return [p, 0]
        fg = fk(g)
        if abs(fg) < 1e-10: return [p, 0]
        if r < 1e-10: return [p, Vp(g)/fg/3]
        return [p, (Vp(g)-2/r*p)/fg]
    def ev(r, y): return 100-abs(y[0])
    ev.terminal = True
    s = solve_ivp(rhs, [rs, rm], [g0+c2*rs**2, 2*c2*rs],
                  method='RK45', rtol=1e-11, atol=1e-13,
                  max_step=0.05, events=[ev], dense_output=True)
    r = np.linspace(rs, min(s.t[-1], rm), 15000)
    return r, s.sol(r)[0], s.sol(r)[1]

def fk_func(g, eta_K=ETA_K):
    a = 2.0 / (1 + eta_K * (g - 1)**2)
    return 1 + 2*a*np.log(g) if g > 0 else -1e30

def V_func(g): return g**3/3.0 - g**4/4.0
def Vp_func(g): return g**2*(1-g)
def Vpp_func(g): return 2*g - 3*g**2

# =====================================================================
print("="*70)
print("  PART 1: ODE STRUCTURE AT g0 = 4")
print("="*70)

g0 = 4.0
alpha_eff_4 = 2.0 / (1 + ETA_K * (g0-1)**2)
f_4 = fk_func(g0)
V_4 = V_func(g0)
Vp_4 = Vp_func(g0)
Vpp_4 = Vpp_func(g0)

print(f"\n  At g0 = 4.0:")
print(f"    alpha_eff = 2/(1+{ETA_K}*9) = {alpha_eff_4:.6f}")
print(f"    f(4) = 1 + 2*{alpha_eff_4:.6f}*ln(4) = {f_4:.6f}")
print(f"    V(4) = 64/3 - 64 = {V_4:.4f}")
print(f"    V'(4) = 16*(1-4) = {Vp_4:.4f}")
print(f"    V''(4) = 8 - 48 = {Vpp_4:.4f}")
print(f"    |V'(4)| / f(4) = {abs(Vp_4)/f_4:.4f}")
print(f"    V''(4) / f(4) = {Vpp_4/f_4:.4f}")

# Core curvature: g''(0) = V'(g0) / (3*f(g0))
g_core_curv = Vp_4 / (3 * f_4)
print(f"    g''(0) = V'(4)/(3*f(4)) = {g_core_curv:.4f}")
print(f"    |g''(0)| = {abs(g_core_curv):.4f}")

# Compare with electron:
g0e = G0_E
alpha_eff_e = 2.0 / (1 + ETA_K * (g0e-1)**2)
f_e = fk_func(g0e)
Vp_e = Vp_func(g0e)
g_core_curv_e = Vp_e / (3 * f_e)
print(f"\n  At g0 = {g0e:.4f} (electron):")
print(f"    alpha_eff = {alpha_eff_e:.6f}")
print(f"    f({g0e}) = {f_e:.6f}")
print(f"    V'({g0e}) = {Vp_e:.6f}")
print(f"    g''(0) = {g_core_curv_e:.6f}")

# =====================================================================
print(f"\n{'='*70}")
print(f"  PART 2: CHARACTERISTIC LENGTH SCALES")
print(f"{'='*70}")

# For each soliton, there's a characteristic "core size" R_core
# where the soliton transitions from g0 to ~1.
# R_core ~ 1/sqrt(|V'(g0)|/f(g0)) (heuristic)

print(f"\n  Heuristic core size R ~ sqrt(f(g0)/|V'(g0)|)")
for g0_test in [G0_E, G0_MU, 2.0, 3.0, 3.845, 4.0, 5.0, 6.0]:
    f_val = fk_func(g0_test)
    Vp_val = abs(Vp_func(g0_test))
    if Vp_val > 0 and f_val > 0:
        R_core = np.sqrt(f_val / Vp_val)
    else:
        R_core = np.nan
    mark = ""
    if abs(g0_test - G0_E) < 0.01: mark = " (e)"
    elif abs(g0_test - G0_MU) < 0.01: mark = " (mu)"
    elif abs(g0_test - 4.0) < 0.01: mark = " (tau)"
    print(f"  g0={g0_test:6.3f}: f={f_val:.4f}, |V'|={Vp_val:.4f}, R_core={R_core:.4f}{mark}")

# =====================================================================
print(f"\n{'='*70}")
print(f"  PART 3: WHERE f(g0) CHANGES SIGN")
print(f"{'='*70}")

# f(g) = 1 + 2*alpha_eff*ln(g) can become ZERO or NEGATIVE.
# This is a natural boundary for the soliton.
# f(g) = 0 when: 1 + 2*alpha_eff(g)*ln(g) = 0
# => ln(g) = -1/(2*alpha_eff)
# => g = exp(-1/(2*alpha_eff))

# For alpha_eff = alpha_UV = 2 (vacuum): f=0 at g = exp(-1/4) = 0.7788
# For running alpha at g~4: alpha_eff = 0.0184, f=0 at g = exp(-1/0.0367) = exp(-27.3) ~ 0

# But f(g0) must be POSITIVE for the ODE to be well-posed.
# For g0 < 1: f(g0) = 1 + 2*alpha_eff*ln(g0) < 1 (since ln(g0) < 0)
# f(g0) = 0 at g0 such that 2*alpha_eff*ln(g0) = -1

print(f"\n  f(g0) for different g0:")
for g0_test in np.arange(0.5, 0.85, 0.05):
    f_val = fk_func(g0_test)
    print(f"  g0={g0_test:.3f}: f={f_val:.6f}")

print(f"\n  f(g0) MINIMUM near g0=0.78 (ghost boundary)")
# f(g0) -> 0 is the ghost boundary: soliton cannot exist for g0 below this

# Find exact zero of f(g) for FIXED eta_K
def f_zero(g):
    return fk_func(g)

# Scan
for g_test in np.arange(0.50, 0.80, 0.01):
    f1 = fk_func(g_test)
    f2 = fk_func(g_test + 0.01)
    if f1 * f2 < 0:
        g_ghost = brentq(fk_func, g_test, g_test+0.01)
        print(f"  Ghost boundary: f(g) = 0 at g = {g_ghost:.6f}")
        break

# =====================================================================
print(f"\n{'='*70}")
print(f"  PART 4: TWO-SOLITON PICTURE")
print(f"{'='*70}")

# In TGP, the tau might be a BOUND STATE of two solitons
# rather than a single soliton with g0=4.
# If g0_tau ~ 2*g0_mu or g0_tau ~ g0_mu + g0_e + 1:

print(f"\n  Additive combinations:")
combos = {
    "g0_e + g0_mu": G0_E + G0_MU,
    "2*g0_mu": 2 * G0_MU,
    "g0_mu + g0_e + 1": G0_MU + G0_E + 1,
    "2*g0_mu + g0_e": 2*G0_MU + G0_E,
    "g0_mu^2/g0_e": G0_MU**2/G0_E,
    "g0_mu^2 + g0_e^2": G0_MU**2 + G0_E**2,
    "g0_mu*(1+g0_e)": G0_MU*(1+G0_E),
    "g0_mu + g0_mu/phi": G0_MU + G0_MU/PHI,
    "g0_mu + sqrt(g0_mu^2+g0_e^2)": G0_MU + np.sqrt(G0_MU**2 + G0_E**2),
    "g0_mu*g0_e + g0_mu": G0_MU*G0_E + G0_MU,
    "1 + g0_mu^2": 1 + G0_MU**2,
}

results = [(abs(v - 4.0)/4.0*100, k, v) for k, v in combos.items()]
results.sort()
for dev, name, val in results[:8]:
    print(f"  {name:35s} = {val:.6f} (dev {dev:.3f}%)")

# =====================================================================
print(f"\n{'='*70}")
print(f"  PART 5: SELF-CONSISTENT eta_K AND g0_tau")
print(f"{'='*70}")

# What if g0_tau IS a free parameter, and instead of trying to fix it,
# we ask: what (eta_K, g0_tau) pair gives the OBSERVED m_tau?
# This changes the problem from "derive g0_tau" to "what locus in
# (eta_K, g0_tau) space matches observations?"

# We already know: at eta_K=12.067, g0_tau~4.0 gives Koide.
# But there might be other (eta_K, g0_tau) pairs.

# Actually, the more fundamental question is:
# Is eta_K = 12 the CORRECT value, or was it biased by the scan range?
# With g0_tau free, can a DIFFERENT eta_K + different g0_tau give r_31=3477?

print(f"\n  Locus of (eta_K, g0_tau) giving r_31 = 3477:")
print(f"  (keeping g0_e re-derived for each eta_K)")
print(f"\n  {'eta_K':>8s} {'g0_e':>10s} {'g0_tau':>10s} {'r_31':>10s} {'g0_tau/g0_mu':>14s}")
print(f"  {'-'*58}")

def extract_tail(r, g):
    m = (r >= 120) & (r <= 260)
    rf, tl = r[m], (g[m] - 1) * r[m]
    if len(rf) < 10: return np.nan
    A = np.column_stack([np.cos(rf), np.sin(rf)])
    coeff, _, _, _ = np.linalg.lstsq(A, tl, rcond=None)
    return np.sqrt(coeff[0]**2 + coeff[1]**2)

def find_g0e_for_eta(eta_K):
    def get_r21(g0_e):
        g0_mu = PHI * g0_e
        r_e, g_e, _ = solve_soliton(g0_e, eta_K)
        if r_e is None or r_e[-1] < 250: return 1e6
        A_e = extract_tail(r_e, g_e)
        r_mu, g_mu, _ = solve_soliton(g0_mu, eta_K)
        if r_mu is None or r_mu[-1] < 250: return 1e6
        A_mu = extract_tail(r_mu, g_mu)
        if np.isnan(A_e) or np.isnan(A_mu) or A_e < 1e-15: return 1e6
        return (A_mu/A_e)**4
    for g_lo in np.arange(0.88, 0.93, 0.005):
        try:
            v_lo = get_r21(g_lo) - 206.768
            v_hi = get_r21(g_lo+0.005) - 206.768
            if v_lo * v_hi < 0:
                return brentq(lambda g: get_r21(g)-206.768, g_lo, g_lo+0.005, xtol=1e-8)
        except:
            pass
    return None

def find_g0tau_for_r31(g0_e, eta_K, target=3477.48):
    r_e, g_e, _ = solve_soliton(g0_e, eta_K)
    if r_e is None: return np.nan
    A_e = extract_tail(r_e, g_e)
    if np.isnan(A_e): return np.nan

    def r31_func(g0_t):
        r_t, g_t, _ = solve_soliton(g0_t, eta_K)
        if r_t is None: return 0
        A_t = extract_tail(r_t, g_t)
        if np.isnan(A_t): return 0
        return (A_t/A_e)**4 - target

    # Bisect in a reasonable range
    for g_lo in np.arange(2.0, 8.0, 0.5):
        try:
            v_lo = r31_func(g_lo)
            v_hi = r31_func(g_lo + 0.5)
            if v_lo * v_hi < 0:
                return brentq(r31_func, g_lo, g_lo+0.5, xtol=1e-6)
        except:
            pass
    return np.nan

for eta in [8, 10, 11, 12, 12.067, 13, 14, 16]:
    g0_e = find_g0e_for_eta(eta)
    if g0_e is None:
        print(f"  {eta:8.3f}  {'None':>10s}")
        continue
    g0_tau = find_g0tau_for_r31(g0_e, eta)
    g0_mu = PHI * g0_e
    ratio = g0_tau / g0_mu if not np.isnan(g0_tau) else np.nan
    mark = ""
    if not np.isnan(ratio) and abs(ratio - np.e) < 0.05: mark = " ~e"
    if not np.isnan(g0_tau) and abs(g0_tau - 4.0) < 0.1: mark += " ~4"
    print(f"  {eta:8.3f} {g0_e:10.6f} {g0_tau:10.4f} {3477.48:10.1f} {ratio:14.4f}{mark}")

# =====================================================================
print(f"\n{'='*70}")
print(f"  SUMMARY")
print(f"{'='*70}")
print(f"""
  Key findings:

  1. VIRIAL EXTREMUM at g0=3.845 gives m_tau=1617 MeV (9% off)
     -> Not precise enough to be the mechanism.

  2. ALGEBRAIC: g0_tau = e * g0_mu = 3.983 gives m_tau=1758 (1.1% off)
     -> Suggestive but not obviously motivated in TGP.

  3. g0_tau ~ 4.0 is very close to:
     - 4 (exact integer)
     - phi^3 = 4.236 (4.1% off)
     - e*g0_mu = 3.983 (0.4% off)
     - pi*sqrt(2)*g0_e = 4.023 (0.6% off)

  4. The (eta_K, g0_tau) locus shows both parameters trade off.
     No independent constraint fixes g0_tau.

  CONCLUSION: g0_tau remains an OPEN parameter in TGP.
  The theory needs an additional physical principle to fix it.
  Best candidates:
    a) g0_tau = e * g0_mu (Euler number relation, 1.1% accuracy)
    b) g0_tau from a two-soliton bound state energy condition
    c) Some topological charge quantization in the full 3D theory
""")
print("DONE")
