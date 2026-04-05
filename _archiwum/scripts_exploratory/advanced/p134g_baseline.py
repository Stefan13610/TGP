#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Quick: r_31 at g0_tau=4.0 with and without running alpha"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from scipy.integrate import solve_ivp
PHI = (1 + np.sqrt(5)) / 2

def solve_soliton(g0, eta_K=0.0, rm=300):
    def fk(g):
        a = 2.0 / (1 + eta_K * (g - 1)**2)
        return 1 + 2*a*np.log(g) if g > 0 else -1e30
    def Vp(g): return g**2*(1-g)
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

def get_A(r, g):
    m = (r >= 120) & (r <= 260)
    rf, tl = r[m], (g[m]-1)*r[m]
    if len(rf) < 10: return np.nan
    M = np.column_stack([np.cos(rf), np.sin(rf)])
    c, _, _, _ = np.linalg.lstsq(M, tl, rcond=None)
    return np.sqrt(c[0]**2 + c[1]**2)

print("Comparison: eta_K=0 (baseline) vs eta_K=12.067 (running)")
print("="*65)

for eta in [0.0, 12.067]:
    print(f"\neta_K = {eta}")
    # Standard g0_e for this eta
    if eta == 0:
        g0_e = 0.89927  # baseline
    else:
        g0_e = 0.905481  # from p131
    g0_mu = PHI * g0_e

    r_e, g_e = solve_soliton(g0_e, eta)
    A_e = get_A(r_e, g_e)
    r_mu, g_mu = solve_soliton(g0_mu, eta)
    A_mu = get_A(r_mu, g_mu)
    r21 = (A_mu/A_e)**4

    print(f"  g0_e={g0_e:.6f}, g0_mu={g0_mu:.6f}")
    print(f"  A_e={A_e:.8f}, A_mu={A_mu:.8f}")
    print(f"  r_21 = {r21:.1f}")

    print(f"\n  g0_tau    A_tau     r_31      m_tau(MeV)")
    for g0_t in [2.0, 3.0, 4.0, 4.5, 5.0]:
        r_t, g_t = solve_soliton(g0_t, eta)
        if r_t is None:
            print(f"  {g0_t:.1f}  FAILED")
            continue
        A_t = get_A(r_t, g_t)
        r31 = (A_t/A_e)**4
        m_tau = 0.51099895 * r31
        print(f"  {g0_t:5.1f}   {A_t:10.6f}  {r31:8.1f}  {m_tau:10.1f}")

print("\n\nKey question: Does running alpha (eta_K=12) qualitatively change")
print("the r_31 vs g0_tau relationship, or just rescale it?")
print("\nDONE")
