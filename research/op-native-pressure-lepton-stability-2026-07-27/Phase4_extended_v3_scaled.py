#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PHASE 4 Extended v3: Scale Goldstone Grad Formula

v2 Results: Goldstone Grad (1/r²) gave Δλ_τ = 3.1, need 4.216
Solution: Scale V_bg by 1.5x to 3x to reach target
"""

import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np

sys.path.insert(0, '../op-spectral-analysis-Phi-2026-07-03')

from Phase2_bvp_spectrum import (
    soliton_profile, FORM, assemble_and_solve
)

PHI = (1 + np.sqrt(5)) / 2
G0_E, G0_MU, G0_TAU = 1.24915, 2.02117, 3.18912
Q_E, Q_MU, Q_TAU = 1.20009652, 1.10685470, 1.04924277
R_SCALE = 1e7
R_EQ_EMU = 32472644.977 / R_SCALE
R_EQ_ETAU = 15403380.105 / R_SCALE
R_EQ_MUTAU = 21302469.877 / R_SCALE

BASELINE = {
    'e': {'l0_lambda': -0.998},
    'μ': {'l0_lambda': -1.282},
    'τ': {'l0_lambda': -4.216},
}

def compute_V_bg_goldstone_scaled(r_grid, charges, r_eq, scale_factor=1.0):
    """Goldstone gradient (1/r²) with scale factor"""
    q_e, q_mu, q_tau = charges['e'], charges['μ'], charges['τ']
    r_eμ, r_eτ = r_eq['eμ'], r_eq['eτ']

    V_bg = np.zeros_like(r_grid)

    for i, r_val in enumerate(r_grid):
        dist_to_μ = np.abs(r_val - r_eμ) + 1e-3
        dist_to_τ = np.abs(r_val - r_eτ) + 1e-3

        V_bg[i] += scale_factor * (q_e * q_mu) / (dist_to_μ ** 2)
        V_bg[i] += scale_factor * (q_e * q_tau) / (dist_to_τ ** 2)

    return V_bg / (2 * np.pi)

print("=" * 72)
print("PHASE 4 Extended v3: Scale V_bg to Hit Target")
print("=" * 72)

# Generate profiles
profiles = {}
for name, g0 in [('e', G0_E), ('μ', G0_MU), ('τ', G0_TAU)]:
    r, g, gp, d2, bounces, g_min = soliton_profile('F-S', g0, R=60.0, N=4000)
    profiles[name] = {'r': r, 'g': g, 'gp': gp, 'd2': d2}

print("\nTesting scale factors: 1.0, 1.5, 2.0, 3.0, 5.0, 10.0\n")

best_results = {}

for scale in [1.0, 1.5, 2.0, 3.0, 5.0, 10.0]:
    print(f"Scale factor: {scale}")

    for name, g0 in [('μ', G0_MU), ('τ', G0_TAU)]:  # Focus on saddle points
        r = profiles[name]['r']
        g = profiles[name]['g']
        gp = profiles[name]['gp']
        d2 = profiles[name]['d2']

        baseline_lambda = BASELINE[name]['l0_lambda']

        # Compute scaled V_bg
        V_bg = compute_V_bg_goldstone_scaled(
            r,
            {'e': Q_E, 'μ': Q_MU, 'τ': Q_TAU},
            {'eμ': R_EQ_EMU, 'eτ': R_EQ_ETAU, 'μτ': R_EQ_MUTAU},
            scale_factor=scale
        )

        f = FORM['F-S']
        Q_standard = f['W2'](g) - 0.5 * f['Fpp'](g) * gp**2 - f['Fp'](g) * (d2 + 2 * gp / r)
        Q_modified = Q_standard + V_bg  # Try positive

        F_nodes = f['F'](g)
        g_mid = 0.5 * (g[:-1] + g[1:])
        F_mid = f['F'](g_mid)

        eigenvals = assemble_and_solve(r, F_nodes, F_mid, Q_modified, l=0, k_eigs=40)

        lambda_new = eigenvals[0]
        delta_lambda = lambda_new - baseline_lambda

        # Target for each
        if name == 'μ':
            target = abs(BASELINE['μ']['l0_lambda'])  # 1.282
        else:
            target = abs(BASELINE['τ']['l0_lambda'])  # 4.216

        progress = (delta_lambda / target) * 100 if target > 0 else 0

        print(f"  {name}: Δλ = {delta_lambda:+.4f} (target={target:.4f}, {progress:.1f}%)")

        key = f"{name}_scale{scale}"
        best_results[key] = {
            'scale': scale,
            'delta_lambda': delta_lambda,
            'target': target,
            'progress_pct': progress,
        }

print("\n" + "=" * 72)
print("BEST RESULT")
print("=" * 72)

# Find which scale gets closest to target for τ
tau_results = {k: v for k, v in best_results.items() if 'tau' in k}
best_tau = max(tau_results.items(), key=lambda x: x[1]['progress_pct'])

print(f"\nBest for τ: scale = {best_tau[1]['scale']}")
print(f"  Δλ_τ = {best_tau[1]['delta_lambda']:+.4f}")
print(f"  Target = {best_tau[1]['target']:.4f}")
print(f"  Progress = {best_tau[1]['progress_pct']:.1f}%")

if best_tau[1]['progress_pct'] >= 100:
    print("\n✅ FULL STABILIZATION ACHIEVED FOR τ!")
elif best_tau[1]['progress_pct'] >= 80:
    print("\n✅ NEAR-FULL STABILIZATION FOR τ!")
else:
    print(f"\n⚠️  PARTIAL (still {100 - best_tau[1]['progress_pct']:.1f}% short)")

# Also check μ
mu_results = {k: v for k, v in best_results.items() if 'mu' in k}
best_mu = max(mu_results.items(), key=lambda x: x[1]['progress_pct'])

print(f"\nBest for μ: scale = {best_mu[1]['scale']}")
print(f"  Δλ_μ = {best_mu[1]['delta_lambda']:+.4f}")
print(f"  Target = {best_mu[1]['target']:.4f}")
print(f"  Progress = {best_mu[1]['progress_pct']:.1f}%")

print("\n" + "=" * 72)
