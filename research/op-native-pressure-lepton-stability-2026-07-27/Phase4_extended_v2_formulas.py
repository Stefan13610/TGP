#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PHASE 4 Extended v2: Test Different V_bg Formulas

Problem: Exponential decay V_bg doesn't work (Δλ too small, N_neg worsens)
Solution: Try alternative physical formulas for pressure potential

Candidates:
  1. V_bg from 1/r falloff (dipole-like, 2D natural)
  2. V_bg from Goldstone gradient (direct pressure gradient)
  3. V_bg from profile distortion (how much does g change in background?)
  4. V_bg with different coupling scales
"""

import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
from scipy.interpolate import interp1d

sys.path.insert(0, '../op-spectral-analysis-Phi-2026-07-03')

try:
    from Phase2_bvp_spectrum import (
        soliton_profile, spectrum_on_background, Q_of, FORM, G_GHOST, assemble_and_solve
    )
    print("✓ Imported CP-7 functions")
except ImportError as e:
    print(f"✗ Failed: {e}")
    sys.exit(1)

# Constants
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

print("=" * 72)
print("PHASE 4 Extended v2: Alternative V_bg Formulas")
print("=" * 72)

# Generate profiles
profiles = {}
for name, g0 in [('e', G0_E), ('μ', G0_MU), ('τ', G0_TAU)]:
    r, g, gp, d2, bounces, g_min = soliton_profile('F-S', g0, R=60.0, N=4000)
    profiles[name] = {'r': r, 'g': g, 'gp': gp, 'd2': d2}
    print(f"{name}: ✓")

# =====================================================================
# FORMULA 1: V_bg from 1/r falloff (dipole-like, natural in 2D)
# =====================================================================

def compute_V_bg_dipole(r_grid, charges, r_eq):
    """
    V_bg ~ Σ (q_i q_j / r_ij^α) with 1/r falloff (α=1, dipole-like)
    This is more natural for 2D than exponential decay.
    """
    q_e, q_mu, q_tau = charges['e'], charges['μ'], charges['τ']
    r_eμ, r_eτ = r_eq['eμ'], r_eq['eτ']

    V_bg = np.zeros_like(r_grid)

    for i, r_val in enumerate(r_grid):
        dist_to_μ = np.abs(r_val - r_eμ) + 0.5  # +0.5 to avoid divergence at r=r_μ
        dist_to_τ = np.abs(r_val - r_eτ) + 0.5

        # 1/r falloff with coupling constants
        V_bg[i] += (q_e * q_mu) / dist_to_μ
        V_bg[i] += (q_e * q_tau) / dist_to_τ

    return V_bg / (2 * np.pi)  # Normalize by 2π (Goldstone scale)

# =====================================================================
# FORMULA 2: V_bg from second derivative of Goldstone energy
# =====================================================================

def compute_V_bg_goldstone_grad(r_grid, charges, r_eq):
    """
    E_press = -Σ q_i q_j log(r_ij) / (2π)
    dE_press/dr ~ -Σ q_i q_j / r_ij
    d²E_press/dr² ~ Σ q_i q_j / r_ij²

    This is the second derivative (curvature) of pressure energy.
    """
    q_e, q_mu, q_tau = charges['e'], charges['μ'], charges['τ']
    r_eμ, r_eτ = r_eq['eμ'], r_eq['eτ']

    V_bg = np.zeros_like(r_grid)

    for i, r_val in enumerate(r_grid):
        dist_to_μ = np.abs(r_val - r_eμ) + 1e-3
        dist_to_τ = np.abs(r_val - r_eτ) + 1e-3

        # Second derivative of log potential
        V_bg[i] += (q_e * q_mu) / (dist_to_μ ** 2)
        V_bg[i] += (q_e * q_tau) / (dist_to_τ ** 2)

    return V_bg / (2 * np.pi)

# =====================================================================
# FORMULA 3: V_bg from profile distortion
# =====================================================================

def compute_V_bg_profile_distortion(r, g, r_eq, interp_dict, charges):
    """
    Estimate V_bg from how much the profile changes when embedded in background.

    Idea: If soliton profile g(r) is distorted by background ψ_bg,
    the distortion energy goes as dE ~ ∫ (∂²E/∂g²) |Δg| dr

    Δg ~ strength of background × coupling_scale
    """
    q_e, q_mu, q_tau = charges['e'], charges['μ'], charges['τ']
    r_eμ, r_eτ = r_eq['eμ'], r_eq['eτ']

    # Strength of background field at core (rough estimate)
    bg_strength_μ = q_mu / r_eμ if r_eμ > 0.1 else q_mu
    bg_strength_τ = q_tau / r_eτ if r_eτ > 0.1 else q_tau

    # Distortion: how much g changes due to background
    # Rough: Δg ~ (q_i * q_j / r_ij) for coupling
    V_bg = np.zeros_like(r)

    for i, r_val in enumerate(r):
        # Distance to other soliton centers
        dist_to_μ = np.abs(r_val - r_eμ) + 0.1
        dist_to_τ = np.abs(r_val - r_eτ) + 0.1

        # Profile curvature gives d²E/dg²
        # For F-S: W''(g) = 2g - 3g²
        W2_g = 2 * g[i] - 3 * g[i]**2

        # V_bg ~ W''(g) * (distortion amplitude)
        V_bg[i] += W2_g * (q_e * q_mu / (dist_to_μ ** 1.5))
        V_bg[i] += W2_g * (q_e * q_tau / (dist_to_τ ** 1.5))

    return V_bg / (2 * np.pi)

# =====================================================================
# Test all three formulas
# =====================================================================

print("\n" + "=" * 72)
print("Testing Three V_bg Formulas")
print("=" * 72)

formulas = {
    'Dipole (1/r)': compute_V_bg_dipole,
    'Goldstone Grad (1/r²)': compute_V_bg_goldstone_grad,
    'Profile Distortion': compute_V_bg_profile_distortion,
}

results_all = {}

for formula_name, formula_func in formulas.items():
    print(f"\n{'='*72}")
    print(f"Formula: {formula_name}")
    print(f"{'='*72}")

    results = {}

    for name, g0 in [('e', G0_E), ('μ', G0_MU), ('τ', G0_TAU)]:
        r = profiles[name]['r']
        g = profiles[name]['g']
        gp = profiles[name]['gp']
        d2 = profiles[name]['d2']

        baseline_lambda = BASELINE[name]['l0_lambda']

        # Compute V_bg with this formula
        if formula_name == 'Profile Distortion':
            V_bg = formula_func(
                r, g,
                {'eμ': R_EQ_EMU, 'eτ': R_EQ_ETAU, 'μτ': R_EQ_MUTAU},
                {},
                {'e': Q_E, 'μ': Q_MU, 'τ': Q_TAU}
            )
        else:
            V_bg = formula_func(
                r,
                {'e': Q_E, 'μ': Q_MU, 'τ': Q_TAU},
                {'eμ': R_EQ_EMU, 'eτ': R_EQ_ETAU, 'μτ': R_EQ_MUTAU}
            )

        # Try BOTH signs
        for sign_name, V_sign in [('positive', +1), ('negative', -1)]:
            f = FORM['F-S']
            Q_standard = f['W2'](g) - 0.5 * f['Fpp'](g) * gp**2 - f['Fp'](g) * (d2 + 2 * gp / r)
            Q_modified = Q_standard + V_sign * V_bg

            F_nodes = f['F'](g)
            g_mid = 0.5 * (g[:-1] + g[1:])
            F_mid = f['F'](g_mid)

            eigenvals = assemble_and_solve(r, F_nodes, F_mid, Q_modified, l=0, k_eigs=40)

            lambda_new = eigenvals[0]
            delta_lambda = lambda_new - baseline_lambda
            n_neg_new = int(np.sum(eigenvals < -1e-6))

            print(f"\n  {name.upper()} ({sign_name} V_bg):")
            print(f"    λ: {baseline_lambda:+.4f} → {lambda_new:+.4f}")
            print(f"    Δλ: {delta_lambda:+.4f}")
            print(f"    N_neg: {BASELINE[name].get('l0_Nneg', 0)} → {n_neg_new}")

            results[f"{name}_{sign_name}"] = {
                'delta_lambda': delta_lambda,
                'lambda_new': lambda_new,
                'n_neg_new': n_neg_new,
            }

    results_all[formula_name] = results

# =====================================================================
# Summary
# =====================================================================

print("\n" + "=" * 72)
print("SUMMARY: Which Formula Works Best?")
print("=" * 72)

for formula_name in formulas.keys():
    print(f"\n{formula_name}:")
    r = results_all[formula_name]

    # Look for best outcome (largest Δλ_μ or Δλ_τ with Δλ > 0)
    delta_mu_pos = r['μ_positive']['delta_lambda']
    delta_mu_neg = r['μ_negative']['delta_lambda']
    delta_tau_pos = r['τ_positive']['delta_lambda']
    delta_tau_neg = r['τ_negative']['delta_lambda']

    best_delta_mu = max(delta_mu_pos, delta_mu_neg)
    best_delta_tau = max(delta_tau_pos, delta_tau_neg)
    best_sign_mu = "positive" if delta_mu_pos > delta_mu_neg else "negative"
    best_sign_tau = "positive" if delta_tau_pos > delta_tau_neg else "negative"

    print(f"  Best μ: Δλ = {best_delta_mu:+.4f} ({best_sign_mu})")
    print(f"  Best τ: Δλ = {best_delta_tau:+.4f} ({best_sign_tau})")

print("\n" + "=" * 72)
print("Test complete. Check which formula gives largest Δλ > 0.")
print("=" * 72)
