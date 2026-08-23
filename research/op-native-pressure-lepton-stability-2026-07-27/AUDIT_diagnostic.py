#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# INDEPENDENT AUDIT DIAGNOSTIC - not trusting prior narrative, measuring directly.
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
sys.path.insert(0, '../op-spectral-analysis-Phi-2026-07-03')
from Phase2_bvp_spectrum import soliton_profile, FORM, assemble_and_solve, vacuum_case

G0_E, G0_MU, G0_TAU = 1.24915, 2.02117, 3.18912

def full_spectrum(g0, R, N, l=0, kk=60):
    r, g, gp, d2, nb, gmin = soliton_profile('F-S', g0, R=R, N=N)
    f = FORM['F-S']
    Q = f['W2'](g) - 0.5*f['Fpp'](g)*gp**2 - f['Fp'](g)*(d2+2*gp/r)
    F_nodes = f['F'](g); g_mid = 0.5*(g[:-1]+g[1:]); F_mid = f['F'](g_mid)
    vals = assemble_and_solve(r, F_nodes, F_mid, Q, l=l, k_eigs=kk)
    return vals, nb, len(r)

print("="*78)
print("AUDIT TEST 1: Does N_neg(<0) scale with BOX SIZE R? (artifact test)")
print("  If N_neg grows with R  -> counting discretized CONTINUUM (artifact)")
print("  If N_neg is R-invariant -> physical bound/unstable modes")
print("="*78)
print(f"{'g0':>8} {'R':>6} {'N':>7} {'N(<0)':>7} {'N(<-1edge)':>11} {'lam_min':>10} {'lam0-(-1)':>10}")
for tag,g0 in [('e',G0_E),('mu',G0_MU),('tau',G0_TAU)]:
    for R in (40.0,60.0,80.0):
        N = int(R/60*4000)
        vals,nb,npts = full_spectrum(g0,R,N,l=0,kk=60)
        n0 = int(np.sum(vals < -1e-6))
        nedge = int(np.sum(vals < -1.0 - 1e-6))   # below F-S continuum edge
        print(f"{g0:>8.4f} {R:>6.0f} {npts:>7} {n0:>7} {nedge:>11} {vals[0]:>10.4f} {vals[0]+1:>10.4f}  bounces={nb}")
    print()

print("="*78)
print("AUDIT TEST 2: Vacuum F-S baseline (empty space, NO soliton)")
print("  If vacuum already has many N_neg(<0) -> those are edge artifacts")
print("="*78)
for R in (40.0,60.0,80.0):
    N = int(R/60*4000)
    r, vals = vacuum_case('F-S', -1.0, R=R, N=N, l=0)
    n0 = int(np.sum(vals < -1e-6))
    nedge = int(np.sum(vals < -1.0 - 1e-6))
    print(f"  VACUUM R={R:>4.0f}: N(<0)={n0:>3}  N(<-1edge)={nedge:>3}  lam_min={vals[0]:+.4f}")

print("="*78)
print("AUDIT TEST 3: Bound modes below F-S edge for e/mu/tau (R=60, the '0/2/3' claim)")
print("  Counting eigenvalues below vacuum edge = TRUE bound-state count")
print("="*78)
for tag,g0 in [('e',G0_E),('mu',G0_MU),('tau',G0_TAU)]:
    vals,nb,npts = full_spectrum(g0,60.0,4000,l=0,kk=60)
    below = vals[vals < -1.0 - 1e-6]
    print(f"  {tag:>3} (g0={g0:.4f}, bounces={nb}): {len(below)} modes below edge -> {np.round(below,4).tolist()}")
print("="*78)
