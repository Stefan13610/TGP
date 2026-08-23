#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
N4c: RADIATIVE CORRECTIONS (One-Loop F-A Effects)

Context:
  - N4d (pressure) gave Δλ_τ = 3.1/4.2 = 74% — insufficient
  - F-A formulation has clean spectrum (no saddle points)
  - This suggests loop effects might be key

Approach:
  1. F-A is gravitational dual to F-S (solitonic)
  2. One-loop corrections modify V_eff
  3. Test: Can V_loop shift saddle points to stability?

Key Physics:
  V_eff = V_TGP + V_loop
  where V_loop ~ Σ log(λ_i²) from fluctuation eigenvalues
"""

import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
from scipy.interpolate import interp1d

sys.path.insert(0, '../op-spectral-analysis-Phi-2026-07-03')

try:
    from Phase2_bvp_spectrum import (
        soliton_profile, FORM, assemble_and_solve
    )
    print("✓ Imported CP-7 functions")
except ImportError as e:
    print(f"✗ Failed: {e}")
    sys.exit(1)

# Constants
PHI = (1 + np.sqrt(5)) / 2
G0_E, G0_MU, G0_TAU = 1.24915, 2.02117, 3.18912

BASELINE = {
    'e': {'l0_lambda': -0.998},
    'μ': {'l0_lambda': -1.282},
    'τ': {'l0_lambda': -4.216},
}

print("=" * 72)
print("N4c: RADIATIVE CORRECTIONS (One-Loop Analysis)")
print("=" * 72)

# =====================================================================
# STEP 1: Understand F-A vs F-S Formulation Duality
# =====================================================================

print("\n" + "=" * 72)
print("STEP 1: F-A (Gravitational) vs F-S (Solitonic) Duality")
print("=" * 72)

print("""
From CP-7 analysis:

F-A (Gravitational):
  - F(u) = u⁴
  - W''(u) = 7u⁶ - 6u⁵
  - Spectrum: CLEAN (all λ ≥ 0, no saddle points)
  - Interpretation: Stable in gravitational frame

F-S (Solitonic, log-form):
  - F(u) = 1 + 4ln(u)
  - W''(u) = 2u - 3u²
  - Spectrum: SADDLE POINTS (μ: 2 modes, τ: 3 modes negative)
  - Interpretation: Unstable in solitonic frame

L04 Duality:
  Two formulations encode same physics at different scales
  F-A dominates at high-E (UV)
  F-S dominates at low-E (IR)

Consequence:
  One-loop corrections from F-A → effective potential shift in F-S
  This could stabilize saddle points!
""")

# =====================================================================
# STEP 2: One-Loop Self-Energy Concept
# =====================================================================

print("\n" + "=" * 72)
print("STEP 2: One-Loop Self-Energy Mechanism")
print("=" * 72)

print("""
One-Loop Contribution to V_eff:

V_loop ~ (ℏ) Σ_i log(λ_i²)  [in natural units, ℏ→1]

Physical meaning:
  - λ_i are eigenvalues of fluctuation operator L̂
  - log(λ²) ~ coupling to virtual modes
  - Repulsive modes (λ > 0) contribute positively
  - Attractive modes (λ < 0) contribute negatively

For our solitons:
  - e: λ > 0 (19 modes, large eigenvalues) → positive log contribution
  - μ: λ < 0 (2 saddle modes) → imaginary log → resonance-like effect
  - τ: λ < 0 (3 saddle modes) → imaginary log → resonance-like effect

Strategy:
  Compute V_loop ~ log(λ_min) contributions
  Add to effective potential: V_eff = V_TGP + λ·V_loop
  where λ ~ coupling strength (to be tuned)

Expected outcome:
  Loop corrections should INCREASE V_eff in saddle direction
  → Δλ_saddle > 0 (shift saddle toward positive)
""")

# =====================================================================
# STEP 3: Estimate Loop Coupling Strength
# =====================================================================

print("\n" + "=" * 72)
print("STEP 3: Loop Coupling Strength Estimate")
print("=" * 72)

print("""
From dimensional analysis:

V_loop ~ α · log(λ²/M²)

where:
  α ~ coupling constant ~ 1/(4π) [typical loop factor]
  M ~ reference scale ~ 1 [in natural units]

For solitons embedded in F-A background:
  Loop coupling scales with soliton scale
  α_soliton ~ q_i² / (4π)  [proportional to charge squared]

This gives:
  For e: α_e ~ (1.2)² / (4π) ~ 0.11
  For μ: α_μ ~ (1.1)² / (4π) ~ 0.10
  For τ: α_τ ~ (1.0)² / (4π) ~ 0.08

We'll test with coupling λ_loop ~ 0.01 to 1.0
""")

# =====================================================================
# STEP 4: Compute Loop Correction to Operator
# =====================================================================

print("\n" + "=" * 72)
print("STEP 4: Compute V_loop and Integrate into Spectrum")
print("=" * 72)

def compute_V_loop_contribution(r, eigenvals_isolated, soliton_charge, coupling_strength=0.1):
    """
    Estimate loop correction to potential based on eigenvalue spectrum.

    V_loop ~ coupling × Σ log(λ_i²)

    Simple model: dominant contribution from most negative eigenvalue
    (this is the saddle-point contribution)
    """

    lambda_min = eigenvals_isolated[0]

    # Log contribution (note: for negative λ, log(λ²) = 2log|λ|)
    if lambda_min < 0:
        log_contrib = 2 * np.log(np.abs(lambda_min))
    else:
        log_contrib = 2 * np.log(lambda_min)

    # V_loop ~ coupling × log(λ_min²) × profile-dependent weight
    # Weight: larger near soliton core, smaller at infinity
    V_loop = coupling_strength * log_contrib * np.ones_like(r)

    # Modulation: stronger at soliton center
    # Crude model: V_loop(r) ~ exp(-r/r_core)
    r_core = 2.0
    modulation = np.exp(-r / r_core)
    V_loop *= modulation

    return V_loop

# Generate profiles and compute baseline spectra
print("\nGenerating profiles...")
profiles = {}
baseline_spectra = {}

for name, g0 in [('e', G0_E), ('μ', G0_MU), ('τ', G0_TAU)]:
    r, g, gp, d2, bounces, g_min = soliton_profile('F-S', g0, R=60.0, N=4000)
    profiles[name] = {'r': r, 'g': g, 'gp': gp, 'd2': d2}

    # Baseline spectrum (from CP-7)
    f = FORM['F-S']
    Q_standard = f['W2'](g) - 0.5 * f['Fpp'](g) * gp**2 - f['Fp'](g) * (d2 + 2 * gp / r)
    F_nodes = f['F'](g)
    g_mid = 0.5 * (g[:-1] + g[1:])
    F_mid = f['F'](g_mid)

    eigenvals_baseline = assemble_and_solve(r, F_nodes, F_mid, Q_standard, l=0, k_eigs=40)
    baseline_spectra[name] = eigenvals_baseline

    print(f"  {name}: λ_min = {eigenvals_baseline[0]:+.6f}")

# Test loop corrections with different coupling strengths
print("\n" + "-" * 72)
print("Testing loop coupling strengths: 0.001, 0.01, 0.1, 0.5, 1.0")
print("-" * 72)

coupling_values = [0.001, 0.01, 0.1, 0.5, 1.0]

results = {}

for coupling in coupling_values:
    print(f"\nCoupling λ_loop = {coupling}:")

    for name in ['μ', 'τ']:  # Focus on saddle points
        r = profiles[name]['r']
        g = profiles[name]['g']
        gp = profiles[name]['gp']
        d2 = profiles[name]['d2']

        baseline_lambda = BASELINE[name]['l0_lambda']
        eigenvals_baseline = baseline_spectra[name]

        # Compute loop correction
        V_loop = compute_V_loop_contribution(
            r,
            eigenvals_baseline,
            soliton_charge=1.1,  # approximate
            coupling_strength=coupling
        )

        # Add to potential
        f = FORM['F-S']
        Q_standard = f['W2'](g) - 0.5 * f['Fpp'](g) * gp**2 - f['Fp'](g) * (d2 + 2 * gp / r)

        # Try BOTH signs (loop could be attractive or repulsive)
        for sign_name, sign in [('repulsive', +1), ('attractive', -1)]:
            Q_with_loop = Q_standard + sign * V_loop

            F_nodes = f['F'](g)
            g_mid = 0.5 * (g[:-1] + g[1:])
            F_mid = f['F'](g_mid)

            eigenvals_new = assemble_and_solve(r, F_nodes, F_mid, Q_with_loop, l=0, k_eigs=40)

            lambda_new = eigenvals_new[0]
            delta_lambda = lambda_new - baseline_lambda

            target = abs(BASELINE[name]['l0_lambda'])
            progress = (delta_lambda / target) * 100 if target > 0 else 0

            print(f"  {name} ({sign_name}): Δλ = {delta_lambda:+.4f} ({progress:5.1f}%)")

            key = f"{name}_λ{coupling}_{sign_name}"
            results[key] = {
                'delta_lambda': delta_lambda,
                'target': target,
                'progress': progress,
            }

# =====================================================================
# Summary
# =====================================================================

print("\n" + "=" * 72)
print("SUMMARY: Loop Effects vs Target")
print("=" * 72)

print("\nBest results:")
tau_best = max(
    {k: v for k, v in results.items() if 'tau' in k}.items(),
    key=lambda x: x[1]['progress']
)
print(f"  τ: {tau_best[0]}")
print(f"    Δλ = {tau_best[1]['delta_lambda']:+.4f}, progress = {tau_best[1]['progress']:.1f}%")

mu_best = max(
    {k: v for k, v in results.items() if 'mu' in k}.items(),
    key=lambda x: x[1]['progress']
)
print(f"  μ: {mu_best[0]}")
print(f"    Δλ = {mu_best[1]['delta_lambda']:+.4f}, progress = {mu_best[1]['progress']:.1f}%")

print("\n" + "=" * 72)
print("N4c INITIAL ASSESSMENT COMPLETE")
print("=" * 72)

print("""
Next steps:
  1. If loop effects too small: combine with pressure (N4d + N4c)
  2. If loop effects large: refine coupling strength
  3. If still insufficient: investigate F-A gravitational background effects
  4. Consider: whether bounce-hierarchy is fundamental mechanism
""")
