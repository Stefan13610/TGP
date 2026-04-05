#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Quick diagnostic: r_31 vs g0_tau at eta_K=12.067, g0_e=0.905481"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from scipy.integrate import solve_ivp

PHI = (1 + np.sqrt(5)) / 2
ETA_K = 12.067
G0_E = 0.905481

def solve_soliton(g0, eta_K=ETA_K, rm=300):
    def fk(g):
        a = 2.0 / (1 + eta_K * (g - 1)**2)
        return 1 + 2*a*np.log(g) if g > 0 else -1e30
    def Vp(g):
        return g**2*(1-g)
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
                  method='RK45', rtol=1e-11, atol=1e-13,
                  max_step=0.05, events=[ev], dense_output=True)
    r = np.linspace(rs, min(s.t[-1], rm), 15000)
    return r, s.sol(r)[0]

def extract_tail(r, g):
    m = (r >= 120) & (r <= 260)
    rf, tl = r[m], (g[m] - 1) * r[m]
    if len(rf) < 10: return np.nan
    A = np.column_stack([np.cos(rf), np.sin(rf)])
    coeff, _, _, _ = np.linalg.lstsq(A, tl, rcond=None)
    return np.sqrt(coeff[0]**2 + coeff[1]**2)

# Get A_e
r_e, g_e = solve_soliton(G0_E)
A_e = extract_tail(r_e, g_e)
print(f"A_e = {A_e:.8f}")
print(f"r_e max = {r_e[-1]:.1f}")

# Scan g0_tau from 1.5 to 6.0
print(f"\n{'g0_tau':>8s} {'A_tau':>12s} {'r_31':>12s} {'r_max':>8s}")
print("-"*45)
for g0_t in np.arange(1.5, 6.5, 0.25):
    r_t, g_t = solve_soliton(g0_t)
    if r_t is None:
        print(f"{g0_t:8.2f}  FAILED")
        continue
    A_t = extract_tail(r_t, g_t)
    if np.isnan(A_t):
        print(f"{g0_t:8.2f}  A=NaN, r_max={r_t[-1]:.0f}")
        continue
    r31 = (A_t / A_e)**4
    print(f"{g0_t:8.2f} {A_t:12.8f} {r31:12.1f} {r_t[-1]:8.0f}")

print("\nDONE")
