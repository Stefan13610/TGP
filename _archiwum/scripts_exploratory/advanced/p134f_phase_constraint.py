#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p134f_phase_constraint.py - What constrains g0_tau?
=====================================================
Check if delta phase spacing selects g0_tau (and thus r_31).
From p131: Delta(e->mu) = 120.01 deg. Does Delta(mu->tau) = 120 deg fix g0_tau?
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

PHI = (1 + np.sqrt(5)) / 2
ETA_K = 12.067
G0_E = 0.905481
G0_MU = PHI * G0_E

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

def extract_tail_full(r, g):
    m = (r >= 120) & (r <= 260)
    rf, tl = r[m], (g[m] - 1) * r[m]
    if len(rf) < 10: return np.nan, np.nan, np.nan, np.nan
    A = np.column_stack([np.cos(rf), np.sin(rf)])
    coeff, _, _, _ = np.linalg.lstsq(A, tl, rcond=None)
    B, C = coeff
    amp = np.sqrt(B**2 + C**2)
    phase = np.arctan2(C, B)  # delta such that xi ~ A*cos(r - delta) = B*cos(r) + C*sin(r)
    return B, C, amp, phase

# Get electron and muon references
r_e, g_e = solve_soliton(G0_E)
B_e, C_e, A_e, d_e = extract_tail_full(r_e, g_e)

r_mu, g_mu = solve_soliton(G0_MU)
B_mu, C_mu, A_mu, d_mu = extract_tail_full(r_mu, g_mu)

print("="*70)
print("  REFERENCE VALUES")
print("="*70)
print(f"  electron: g0={G0_E:.6f}, A={A_e:.8f}, delta={np.degrees(d_e):.2f} deg")
print(f"  muon:     g0={G0_MU:.6f}, A={A_mu:.8f}, delta={np.degrees(d_mu):.2f} deg")
delta_em = d_mu - d_e
# Normalize to [0, 360)
delta_em_deg = np.degrees(delta_em) % 360
print(f"  Delta(e->mu) = {delta_em_deg:.2f} deg")
print(f"  r_21 = (A_mu/A_e)^4 = {(A_mu/A_e)**4:.1f}")

# Scan tau: amplitude, phase, r_31
print(f"\n{'='*70}")
print(f"  TAU SCAN: g0_tau from 1.5 to 8.0")
print(f"{'='*70}")
print(f"  {'g0':>6s} {'A_tau':>10s} {'delta':>10s} {'D(mu->tau)':>12s} {'D(e->tau)':>11s} {'r_31':>10s}")
print(f"  {'-'*65}")

for g0_t in np.arange(1.5, 8.1, 0.25):
    r_t, g_t = solve_soliton(g0_t)
    if r_t is None:
        print(f"  {g0_t:6.2f}  FAILED")
        continue
    B_t, C_t, A_t, d_t = extract_tail_full(r_t, g_t)
    if np.isnan(A_t):
        print(f"  {g0_t:6.2f}  NaN")
        continue
    r31 = (A_t / A_e)**4
    d_mt = np.degrees(d_t - d_mu) % 360
    if d_mt > 180: d_mt -= 360
    d_et = np.degrees(d_t - d_e) % 360
    if d_et > 180: d_et -= 360
    mark = ""
    if abs(d_mt - 120) < 5: mark = " <-- 120?"
    if abs(d_mt + 120) < 5: mark = " <-- -120?"
    if abs(d_et - 240) < 5: mark += " 240!"
    if abs(r31 - 3477.5) < 100: mark += " KOIDE!"
    print(f"  {g0_t:6.2f} {A_t:10.6f} {np.degrees(d_t):10.2f} {d_mt:12.2f} {d_et:11.2f} {r31:10.1f}{mark}")

# Check: is there a g0_tau where Delta(mu->tau) = 120 deg?
print(f"\n{'='*70}")
print(f"  PHASE CONDITION: find g0_tau where Delta(mu->tau) = +-120 deg")
print(f"{'='*70}")

target_phases = [120.0, -120.0, 240.0, -240.0]

for target in target_phases:
    target_rad = np.radians(target)
    print(f"\n  Target: Delta(mu->tau) = {target:.0f} deg")

    prev_val = None
    prev_g0 = None
    found = False

    for g0_t in np.arange(1.5, 10.0, 0.1):
        r_t, g_t = solve_soliton(g0_t)
        if r_t is None: continue
        _, _, A_t, d_t = extract_tail_full(r_t, g_t)
        if np.isnan(A_t): continue
        d_mt = d_t - d_mu
        # Normalize
        val = ((d_mt - target_rad + np.pi) % (2*np.pi)) - np.pi

        if prev_val is not None and prev_val * val < 0 and abs(val) < 1:
            # Bracket found, bisect
            g_lo, g_hi = prev_g0, g0_t
            for _ in range(30):
                g_mid = (g_lo + g_hi) / 2
                r_m, g_m = solve_soliton(g_mid)
                if r_m is None:
                    g_lo = g_mid
                    continue
                _, _, A_m, d_m = extract_tail_full(r_m, g_m)
                if np.isnan(A_m):
                    g_lo = g_mid
                    continue
                v_mid = ((d_m - d_mu - target_rad + np.pi) % (2*np.pi)) - np.pi
                if v_mid * prev_val < 0:
                    g_hi = g_mid
                else:
                    g_lo = g_mid
                    prev_val = v_mid
                if abs(g_hi - g_lo) < 1e-6:
                    break

            g_sol = (g_lo + g_hi) / 2
            r_sol, g_sol_data = solve_soliton(g_sol)
            _, _, A_sol, d_sol = extract_tail_full(r_sol, g_sol_data)
            r31_sol = (A_sol / A_e)**4
            d_mt_sol = np.degrees(d_sol - d_mu) % 360
            if d_mt_sol > 180: d_mt_sol -= 360

            print(f"    FOUND: g0_tau = {g_sol:.6f}")
            print(f"           A_tau = {A_sol:.8f}")
            print(f"           r_31 = {r31_sol:.1f}")
            print(f"           m_tau = {0.51099895 * r31_sol:.1f} MeV (obs: 1776.86)")
            print(f"           Delta(mu->tau) = {d_mt_sol:.4f} deg")
            found = True

        prev_val = val
        prev_g0 = g0_t

    if not found:
        print(f"    Not found in [1.5, 10.0]")

print(f"\n{'='*70}")
print(f"  CONCLUSION")
print(f"{'='*70}")
print(f"  If the phase condition Delta(mu->tau) = 120 deg")
print(f"  uniquely determines g0_tau, then r_31 is a TRUE prediction.")
print(f"  Otherwise, g0_tau must be constrained by another mechanism.")
print("DONE")
