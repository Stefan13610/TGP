#!/usr/bin/env python3
# Quick version - test g0 interpolation
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
sys.path.insert(0, '../op-spectral-analysis-Phi-2026-07-03')

from Phase2_bvp_spectrum import soliton_profile, FORM, assemble_and_solve
from scipy.stats import pearsonr

G0_E, G0_MU, G0_TAU = 1.24915, 2.02117, 3.18912
G_GHOST = np.exp(-1/4)

print("="*80)
print("QUICK INTERPOLATION TEST: λ_min vs Bounce Count")
print("="*80)

# Sample 8 intermediate values
g0_vals = np.linspace(G0_E, G0_TAU, 8)
data = []

print("\nSampling 8 g₀ values from e to τ:\n")
print(f"{'g₀':>10} {'Bounces':>10} {'g_min':>10} {'λ_min':>12} {'N_neg':>8}")
print("-"*50)

for g0 in g0_vals:
    try:
        r, g, gp, d2, bounces, g_min = soliton_profile('F-S', g0, R=60.0, N=4000)
        f = FORM['F-S']
        Q = f['W2'](g) - 0.5*f['Fpp'](g)*gp**2 - f['Fp'](g)*(d2+2*gp/r)
        F_nodes = f['F'](g)
        g_mid = 0.5*(g[:-1]+g[1:])
        F_mid = f['F'](g_mid)
        eigenvals = assemble_and_solve(r, F_nodes, F_mid, Q, l=0, k_eigs=40)
        lambda_min = eigenvals[0]
        n_neg = np.sum(eigenvals < -1e-6)
        print(f"{g0:>10.6f} {bounces:>10d} {g_min:>10.6f} {lambda_min:>+12.6f} {n_neg:>8d}")
        data.append({'g0': g0, 'bounces': bounces, 'lambda': lambda_min, 'n_neg': n_neg, 'g_min': g_min})
    except Exception as e:
        print(f"{g0:>10.6f} ERROR: {str(e)[:40]}")

if len(data) > 1:
    print("\n" + "="*80)
    print("ANALYSIS: Deterministic N_neg by Bounce Count")
    print("="*80)
    
    bounces_arr = np.array([d['bounces'] for d in data])
    lambda_arr = np.array([d['lambda'] for d in data])
    nneg_arr = np.array([d['n_neg'] for d in data])
    
    # Group by bounces
    bounce_map = {}
    for d in data:
        b = d['bounces']
        if b not in bounce_map:
            bounce_map[b] = []
        bounce_map[b].append(d)
    
    print("\nBounce Count → N_neg Mapping:")
    for b in sorted(bounce_map.keys()):
        vals = [d['n_neg'] for d in bounce_map[b]]
        unique = set(vals)
        status = "✓ DETERMINISTIC" if len(unique)==1 else "⚠ VARIES"
        print(f"  bounces={b}: N_neg={vals} {status}")
    
    # Correlation
    try:
        corr, pval = pearsonr(bounces_arr, lambda_arr)
        print(f"\nCorrelation: bounces ↔ λ_min = {corr:+.4f} (p={pval:.2e})")
        if abs(corr) > 0.9:
            print("  → STRONG correlation (bounces determine λ_min)")
        elif abs(corr) > 0.7:
            print("  → MODERATE correlation")
        else:
            print("  → WEAK correlation")
    except:
        pass

print("\n" + "="*80)
