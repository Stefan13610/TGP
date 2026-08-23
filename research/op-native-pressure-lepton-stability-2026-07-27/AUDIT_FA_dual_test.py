#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# DECISIVE TEST: Do the F-S sub-edge modes (0/2/3) survive in the CANONICAL F-A formulation?
# L04 verdict: canonical TGP is alpha=2, K=phi^4 == F-A. F-S(log) is non-canonical.
# If F-A solitons are clean (no sub-edge modes) -> the 0/2/3 structure is an F-S artifact.
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from scipy.integrate import solve_ivp
sys.path.insert(0, '../op-spectral-analysis-Phi-2026-07-03')
from Phase2_bvp_spectrum import soliton_profile, FORM, assemble_and_solve, spectrum_on_background

G0_E, G0_MU, G0_TAU = 1.24915, 2.02117, 3.18912

# ---------------------------------------------------------------
# F-A soliton ODE, derived from same action as F-S:
#   F(u)(u''+2u'/r) + (1/2)F'(u)u'^2 = W'(u)
#   F-A: F=u^4, F'=4u^3, W'=u^7-u^6  =>  u'' + 2u'/r = u^3 - u^2 - 2u'^2/u
# ---------------------------------------------------------------
def fa_soliton(g0, R=60.0, N=4000):
    def rhs(r_, y):
        u, up = y
        u = max(u, 1e-6)
        if r_ < 1e-10:
            return [up, (u**3 - u**2)/3.0]
        return [up, u**3 - u**2 - 2*up**2/u - 2*up/r_]
    def blow(r_, y):
        return abs(y[0]) - 12.0          # detect runaway
    blow.terminal = True
    h = R/N
    r_grid = (np.arange(N)+0.5)*h
    sol = solve_ivp(rhs, [1e-6, R+h], [g0, 0.0], method='DOP853',
                    max_step=0.02, rtol=1e-10, atol=1e-13,
                    events=[blow], dense_output=True)
    end = sol.t[-1]
    runaway = (sol.t_events[0].size > 0)
    sel = (r_grid >= 1e-6) & (r_grid <= end)
    rr = r_grid[sel]
    yy = sol.sol(rr)
    u, up = yy[0], yy[1]
    d2 = np.array([rhs(rr[j], [u[j], up[j]])[1] for j in range(len(rr))])
    return rr, u, up, d2, runaway, float(u[-1]), float(np.min(u)), float(np.max(u))

print("="*80)
print("DECISIVE F-A TEST: canonical formulation soliton spectra")
print("="*80)
print("F-A vacuum edge = +1 (W''(1)=+1, STABLE).  F-S vacuum edge = -1 (W''(1)=-1, tachyonic).")
print()

# First: does F-A even SUPPORT crown solitons (g0>1 -> 1)?
print("STEP 1: Can F-A support crown solitons (g0>1 descending to vacuum=1)?")
print(f"{'g0':>10} {'runaway?':>10} {'u_end':>10} {'u_min':>10} {'u_max':>10} {'pts':>6}")
for tag,g0 in [('e',G0_E),('mu',G0_MU),('tau',G0_TAU)]:
    rr,u,up,d2,ra,uend,umin,umax = fa_soliton(g0)
    print(f"{g0:>10.4f} {str(ra):>10} {uend:>10.4f} {umin:>10.4f} {umax:>10.2f} {len(rr):>6}  ({tag})")

print()
print("STEP 2: Scan a RANGE of g0 in F-A (be fair - maybe crown solitons live elsewhere)")
print(f"{'g0':>8} {'runaway?':>9} {'u_end':>9} {'u_max':>9}")
for g0 in [1.01,1.05,1.1,1.2,1.4,1.7,2.0,2.5,3.0]:
    rr,u,up,d2,ra,uend,umin,umax = fa_soliton(g0)
    print(f"{g0:>8.2f} {str(ra):>9} {uend:>9.4f} {umax:>9.2f}")

print()
print("STEP 3: For any F-A profile that stays bounded, compute spectrum & count sub-edge modes")
print("  (F-A edge=+1, so 'bound/unstable' = modes below +1; 'negative' = below 0)")
print(f"{'g0':>8} {'bounded?':>9} {'lam_min':>10} {'N(<0)':>7} {'N(<+1edge)':>11}")
for g0 in [1.01,1.05,1.1,1.2,1.4,1.7,2.0]:
    rr,u,up,d2,ra,uend,umin,umax = fa_soliton(g0)
    bounded = (not ra) and umax < 11
    if not bounded:
        print(f"{g0:>8.2f} {'NO':>9} {'--':>10} {'--':>7} {'--':>11}")
        continue
    try:
        vals = spectrum_on_background('F-A', rr, u, up, d2, l=0)
        n0 = int(np.sum(vals < -1e-6))
        nedge = int(np.sum(vals < 1.0 - 1e-6))
        print(f"{g0:>8.2f} {'YES':>9} {vals[0]:>10.4f} {n0:>7} {nedge:>11}")
    except Exception as ex:
        print(f"{g0:>8.2f} {'ERR':>9} {str(ex)[:30]}")

print()
print("="*80)
print("STEP 4: Side-by-side reminder of F-S (non-canonical) sub-edge counts")
print("="*80)
for tag,g0 in [('e',G0_E),('mu',G0_MU),('tau',G0_TAU)]:
    r,g,gp,d2,nb,gmin = soliton_profile('F-S', g0, R=60.0, N=4000)
    vals = spectrum_on_background('F-S', r, g, gp, d2, l=0)
    below = vals[vals < -1.0 - 1e-6]
    print(f"  F-S {tag:>3} (g0={g0:.4f}, bounces={nb}): {len(below)} modes below edge(-1) = {np.round(below,3).tolist()}")
print("="*80)
