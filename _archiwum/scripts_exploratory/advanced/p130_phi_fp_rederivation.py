#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p130_phi_fp_rederivation.py -- Re-derive phi-FP with running alpha
===================================================================

ULTRA-MINIMAL VERSION for speed.

Key question: for each eta_K, find g0_e such that r_21 = 206.768
(keeping g0_mu = phi*g0_e), then check max r_31.

Author: TGP project, session v42+
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

PHI = (1 + np.sqrt(5)) / 2
R21_PDG = 206.768
R31_PDG = 3477.48

def solve_soliton(g0, eta_K=0.0, p_run=2, c5=0.0, rm=300):
    """ODE with running alpha and optional V correction."""
    def fk(g):
        a = 2.0 / (1 + eta_K * abs(g - 1)**p_run)
        return 1 + 2*a*np.log(g) if g > 0 else -1e30
    def Vp(g):
        return g**2*(1-g) + c5*g**4*(1-g)
    fg0 = fk(g0)
    if abs(fg0) < 1e-15: return None, None
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
                  method='RK45', rtol=1e-10, atol=1e-12,
                  max_step=0.1, events=[ev], dense_output=True)
    r = np.linspace(rs, min(s.t[-1], rm), 15000)
    return r, s.sol(r)[0]

def get_A(r, g):
    m = (r >= 120) & (r <= 260)
    rf, tl = r[m], (g[m] - 1) * r[m]
    if len(rf) < 10: return np.nan, np.nan
    A = np.column_stack([np.cos(rf), np.sin(rf)])
    coeff, _, _, _ = np.linalg.lstsq(A, tl, rcond=None)
    B, C = coeff
    return np.sqrt(B**2 + C**2), np.arctan2(C, B)

def get_r21(g0_e, eta_K=0.0, p_run=2, c5=0.0):
    g0_mu = PHI * g0_e
    r_e, g_e = solve_soliton(g0_e, eta_K, p_run, c5)
    if r_e is None or r_e[-1] < 250: return np.nan
    A_e, _ = get_A(r_e, g_e)
    r_mu, g_mu = solve_soliton(g0_mu, eta_K, p_run, c5)
    if r_mu is None or r_mu[-1] < 250: return np.nan
    A_mu, _ = get_A(r_mu, g_mu)
    if np.isnan(A_e) or np.isnan(A_mu) or A_e < 1e-15: return np.nan
    return (A_mu / A_e)**4

def find_g0e(eta_K, p_run=2, c5=0.0):
    """Find g0_e giving r_21 = 206.768. Fast: 15-pt scan + brentq."""
    def obj(g0_e):
        r21 = get_r21(g0_e, eta_K, p_run, c5)
        return r21 - R21_PDG if not np.isnan(r21) else 1e6
    for g_lo, g_hi in zip(np.linspace(0.70, 0.95, 14), np.linspace(0.72, 0.97, 14)):
        try:
            v_lo, v_hi = obj(g_lo), obj(g_hi)
            if v_lo * v_hi < 0:
                return brentq(obj, g_lo, g_hi, xtol=1e-7)
        except: pass
    return None

def max_r31(g0_e, eta_K, p_run=2, c5=0.0):
    """Max r_31. Coarse 20-pt scan."""
    r_e, g_e = solve_soliton(g0_e, eta_K, p_run, c5)
    if r_e is None or r_e[-1] < 250: return 0, 0
    A_e, d_e = get_A(r_e, g_e)
    if np.isnan(A_e): return 0, 0
    best, bg0 = 0, 0
    for g0 in np.linspace(1.3, 4.0, 25):
        r_t, g_t = solve_soliton(g0, eta_K, p_run, c5)
        if r_t is None or r_t[-1] < 250: continue
        A_t, _ = get_A(r_t, g_t)
        if np.isnan(A_t): continue
        r31 = (A_t / A_e)**4
        if r31 > best: best, bg0 = r31, g0
    return best, bg0

# ================================================================
#  PART A: Standard running alpha (p=2)
# ================================================================
print("=" * 70)
print("  PART A: phi-FP RE-DERIVATION WITH RUNNING alpha (p=2)")
print("=" * 70)
print("  eta_K   g0_e(new)  r_21     max_r31   ratio_to_3477")
print("  " + "-" * 60)

for eta in [0, 2, 5, 8, 12, 18, 25, 35, 50]:
    g0_e = find_g0e(eta)
    if g0_e is None:
        print("  %5.1f   (no phi-FP solution)" % eta)
        continue
    r21 = get_r21(g0_e, eta)
    r31, g0t = max_r31(g0_e, eta)
    ratio = r31 / R31_PDG
    mk = " <<<" if 0.85 < ratio < 1.15 else ""
    print("  %5.1f   %.6f   %7.1f   %8.1f    %.3f%s" % (eta, g0_e, r21, r31, ratio, mk))

# ================================================================
#  PART B: Steeper running (p=4, 6)
# ================================================================
print("\n" + "=" * 70)
print("  PART B: STEEPER RUNNING |g-1|^p")
print("=" * 70)
print("  p   eta    g0_e      r_21     max_r31   ratio")
print("  " + "-" * 55)

for p_run in [4, 6]:
    for eta in [10, 30, 100, 300]:
        g0_e = find_g0e(eta, p_run=p_run)
        if g0_e is None: continue
        r21 = get_r21(g0_e, eta, p_run=p_run)
        r31, g0t = max_r31(g0_e, eta, p_run=p_run)
        ratio = r31 / R31_PDG
        mk = " <<<" if 0.85 < ratio < 1.15 else ""
        print("  %d   %5.0f  %.6f   %7.1f   %8.1f    %.3f%s" % (p_run, eta, g0_e, r21, r31, ratio, mk))

# ================================================================
#  PART C: Potential correction only (c5)
# ================================================================
print("\n" + "=" * 70)
print("  PART C: POTENTIAL CORRECTION c5*g^4*(1-g)")
print("=" * 70)
print("  c5     g0_e      r_21     max_r31   ratio")
print("  " + "-" * 50)

for c5 in [-0.3, -0.1, 0, 0.1, 0.3, 0.5, 1.0, 2.0, 5.0]:
    g0_e = find_g0e(0, c5=c5)
    if g0_e is None:
        print("  %5.2f  (no solution)" % c5)
        continue
    r21 = get_r21(g0_e, 0, c5=c5)
    r31, g0t = max_r31(g0_e, 0, c5=c5)
    ratio = r31 / R31_PDG
    mk = " <<<" if 0.85 < ratio < 1.15 else ""
    print("  %5.2f  %.6f   %7.1f   %8.1f    %.3f%s" % (c5, g0_e, r21, r31, ratio, mk))

# ================================================================
#  PART D: Combined (eta + c5)
# ================================================================
print("\n" + "=" * 70)
print("  PART D: COMBINED eta + c5")
print("=" * 70)
print("  eta   c5     g0_e      r_21     max_r31   ratio")
print("  " + "-" * 55)

for eta, c5 in [(5, 0.5), (5, 1.0), (10, 0.5), (10, 1.0), (15, 0.5), (20, 0.5)]:
    g0_e = find_g0e(eta, c5=c5)
    if g0_e is None: continue
    r21 = get_r21(g0_e, eta, c5=c5)
    r31, g0t = max_r31(g0_e, eta, c5=c5)
    ratio = r31 / R31_PDG
    mk = " <<<" if 0.85 < ratio < 1.15 else ""
    print("  %4.0f  %5.2f  %.6f   %7.1f   %8.1f    %.3f%s" % (eta, c5, g0_e, r21, r31, ratio, mk))

# ================================================================
#  PART E: Best solution - full phase analysis
# ================================================================
print("\n" + "=" * 70)
print("  PART E: SUMMARY")
print("=" * 70)
print()
print("  The key question: can SIMULTANEOUS r_21=206.8 AND r_31=3477")
print("  be achieved with any single-parameter correction?")
print()
print("  If ratio ~ 1.0 for any row above, the answer is YES.")
print("  If all ratios << 1, tau requires qualitatively new physics.")
