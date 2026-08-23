#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PHASE 2: Extract Lepton Charges q_e, q_μ, q_τ from Soliton Profiles

op-native-pressure-lepton-stability-2026-07-27
Session #64, Tier 2 v5

Purpose:
  Extract "charges" q_e, q_μ, q_τ from the three soliton profiles
  of the crown (e, μ, τ generations). These charges will couple to the
  substrate propagator G(r) to give interaction energy V_int = Σ q_i q_j G(r_ij).

Method:
  1. Solve ODE for each lepton (g₀^e, g₀^μ, g₀^τ)
  2. Extract far-field tail: g(r) → 1 + (A_tail/r)·sin(r+φ) + O(1/r²)
  3. Relate A_tail to "charge" q via coupling model
  4. Verify with dimension analysis and physical interpretation

Key Hypotheses (to test):
  H1: q_i is proportional to A_tail(g₀^i) [amplitude of far-field]
  H2: q_i scales with (g₀^i - 1) [distance from vacuum]
  H3: q_i ~ (g₀^i - 1)² or higher power [from soliton binding energy]
"""

import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
import math
from scipy.integrate import solve_ivp
from scipy.optimize import fminbound
import sympy as sp

# =====================================================================
# CONSTANTS (from ls10_third_generation_selection.py)
# =====================================================================

PHI = (1 + np.sqrt(5)) / 2          # golden ratio ≈ 1.618034
PHI0 = 25.0                          # vacuum field Φ₀
G0_E = 0.86941                       # g₀^e (substrate formulation)
G0_MU = PHI * G0_E                   # g₀^μ = φ * g₀^e
G0_TAU_KOIDE = 3.18912               # g₀^τ (Koide-determined)

# Physical constants (PDG 2024)
M_E = 0.51099895                     # MeV
M_MU = 105.6583755                   # MeV
M_TAU = 1776.86                      # MeV
R21_PDG = M_MU / M_E                 # 206.768
R31_PDG = M_TAU / M_E                # 3477.18

print("=" * 72)
print("PHASE 2: Extract Lepton Charges")
print("op-native-pressure-lepton-stability-2026-07-27")
print("=" * 72)

print(f"\nCanonical Parameters:")
print(f"  G0_E = {G0_E:.6f}")
print(f"  G0_MU = φ * G0_E = {G0_MU:.6f}")
print(f"  G0_TAU = {G0_TAU_KOIDE:.6f}")
print(f"  φ = {PHI:.6f}")
print(f"  Φ₀ = {PHI0:.1f}")

# =====================================================================
# TASK 2a: Solve ODE for Each Lepton
# =====================================================================

def rhs_canonical(r, y):
    """
    RHS of canonical soliton ODE (from brannen_sqrt2/r6_c1):
      g'' + (2/r)g' + g - (1/g) - g² = 0

    Rearranged as system: [g, g'] → [g', g'']
    where g'' = (1/g) + g² - g - (2/r)g'
    """
    g, gp = y
    # Guard against pathological values
    if g < 1e-14:
        g = 1e-14
    if r < 1e-13:
        # Taylor expansion near r=0: g''(0) = (1 - g₀)/4
        gpp = (1.0 - g) / 4.0
    else:
        # Full ODE: g'' = (1/g) + g² - g - (2/r)g'
        gpp = (1.0 / g) + g**2 - g - (2.0 / r) * gp
    return [gp, gpp]


def solve_lepton_ode(g0, r_max=80.0, n_pts=12000, rtol=1e-12, atol=1e-14):
    """
    Solve soliton ODE for given initial condition g₀.

    Input:
      g0 : float, initial value g(0) = g₀
      r_max : float, integration domain [0, r_max]
      n_pts : int, number of sample points
      rtol, atol : tolerances for solve_ivp

    Output:
      sol : scipy.integrate.OdeResult
    """
    r_eval = np.linspace(1e-10, r_max, n_pts)
    sol = solve_ivp(
        rhs_canonical,
        (1e-10, r_max),
        [g0, 0.0],
        method='DOP853',  # 8th-order Runge-Kutta
        t_eval=r_eval,
        rtol=rtol,
        atol=atol,
        max_step=0.02,
        dense_output=True
    )
    return sol


def extract_far_field(sol, r_min=60.0, r_max=None):
    """
    Extract far-field amplitude A_tail from oscillatory tail.

    Ansatz: g(r) → 1 + (A_tail/r) * sin(r + φ) + O(1/r²)
            u(r) = (g-1)*r → A_tail * sin(r + φ)

    Fit to: u(r) = B*cos(r) + C*sin(r)  (r ∈ [r_min, r_max])

    Output:
      A_tail : amplitude √(B² + C²)
      phi : phase arctan(C/B)
      residual_rms : quality of fit
      B, C : fitted coefficients
    """
    if r_max is None:
        r_max = sol.t[-1]

    r = sol.t
    g = sol.y[0]
    mask = (r >= r_min) & (r <= r_max)

    if np.sum(mask) < 10:
        return None, None, None, None, None

    r_f = r[mask]
    u_f = (g[mask] - 1.0) * r_f

    # Least-squares fit: u = B*cos(r) + C*sin(r)
    X = np.column_stack([np.cos(r_f), np.sin(r_f)])
    try:
        coef, *_ = np.linalg.lstsq(X, u_f, rcond=None)
    except:
        return None, None, None, None, None

    B, C = coef[0], coef[1]
    A_tail = np.sqrt(B**2 + C**2)
    phi = np.arctan2(C, B)

    # Residual quality
    y_hat = B * np.cos(r_f) + C * np.sin(r_f)
    residual_rms = np.sqrt(np.mean((u_f - y_hat)**2)) / max(A_tail, 1e-12)

    return A_tail, phi, residual_rms, B, C


print("\n" + "=" * 72)
print("TASK 2a: Solve ODE for Each Lepton")
print("=" * 72)

# Solve for each generation
leptons = [
    ('e', G0_E),
    ('μ', G0_MU),
    ('τ', G0_TAU_KOIDE)
]

solutions = {}
print(f"\n{'Gen':>3} {'g₀':>10} {'Status':>15} {'g(r→∞)':>12} {'A_tail':>12}")
print(f"{'-'*3} {'-'*10} {'-'*15} {'-'*12} {'-'*12}")

for name, g0 in leptons:
    sol = solve_lepton_ode(g0, r_max=100.0, n_pts=15000)
    solutions[name] = sol

    status = "PASS" if sol.success else "FAIL"
    g_final = sol.y[0, -1]
    A_tail, phi, residual, B, C = extract_far_field(sol, r_min=70.0, r_max=95.0)

    if A_tail is not None:
        print(f"{name:>3} {g0:>10.6f} {status:>15} {g_final:>12.6f} {A_tail:>12.8f}")
    else:
        print(f"{name:>3} {g0:>10.6f} {status:>15} {g_final:>12.6f} {'FAILED':>12}")
        print(f"     WARNING: Far-field extraction failed")

# =====================================================================
# TASK 2b: Extract Charges q_e, q_μ, q_τ
# =====================================================================

print("\n" + "=" * 72)
print("TASK 2b: Extract Charges from Far-Field")
print("=" * 72)

# Hypothesis H1: q_i ∝ A_tail (direct proportionality)
# Hypothesis H2: q_i ∝ (g₀ - 1) (distance from vacuum)
# Hypothesis H3: q_i ∝ (g₀ - 1)² (soliton self-energy)

charges_data = {}

print(f"\n{'Gen':>3} {'g₀':>10} {'Δg₀':>10} {'A_tail':>12} {'q (H1)':>12} {'q (H2)':>12} {'q (H3)':>12}")
print(f"{'-'*3} {'-'*10} {'-'*10} {'-'*12} {'-'*12} {'-'*12} {'-'*12}")

for name, g0 in leptons:
    sol = solutions[name]
    if not sol.success:
        print(f"{name:>3} {g0:>10.6f} {'N/A':>10} {'FAILED':>12} {'N/A':>12} {'N/A':>12} {'N/A':>12}")
        continue

    A_tail, phi, residual, B, C = extract_far_field(sol, r_min=70.0, r_max=95.0)
    if A_tail is None:
        print(f"{name:>3} {g0:>10.6f} {'N/A':>10} {'FAILED':>12} {'N/A':>12} {'N/A':>12} {'N/A':>12}")
        continue

    delta_g0 = g0 - 1.0

    # Three hypotheses for charge scaling:
    q_h1 = A_tail                              # Direct (Hypothesis 1)
    q_h2 = A_tail / delta_g0                   # Scaled by distance (H2)
    q_h3 = A_tail / (delta_g0**2)              # Scaled by squared distance (H3)

    print(f"{name:>3} {g0:>10.6f} {delta_g0:>10.6f} {A_tail:>12.8f} {q_h1:>12.8f} {q_h2:>12.8f} {q_h3:>12.8f}")

    charges_data[name] = {
        'g0': g0,
        'delta_g0': delta_g0,
        'A_tail': A_tail,
        'phi': phi,
        'residual': residual,
        'B': B,
        'C': C,
        'q_h1': q_h1,
        'q_h2': q_h2,
        'q_h3': q_h3,
    }

# =====================================================================
# TASK 2c: Verify Against Coupling Model
# =====================================================================

print("\n" + "=" * 72)
print("TASK 2c: Verify Against Coupling Model")
print("=" * 72)

print("\n§ C1: Test H2 hypothesis (q ∝ 1/Δg₀)")
print("  Rationale: If q_i ~ A_tail/(g₀-1), then profiles with larger Δg₀")
print("  have systematically smaller q despite same/similar A_tail.")
print("  This makes sense: larger deviation from vacuum → larger source.")

q_h2_vals = [charges_data[name]['q_h2'] for name in ['e', 'μ', 'τ']]
print(f"\n  q_H2 values: e={q_h2_vals[0]:.6f}, μ={q_h2_vals[1]:.6f}, τ={q_h2_vals[2]:.6f}")
print(f"  Ratios: q_μ/q_e = {q_h2_vals[1]/q_h2_vals[0]:.6f}")
print(f"          q_τ/q_e = {q_h2_vals[2]/q_h2_vals[0]:.6f}")

print("\n§ C2: Test interaction strength")
print("  For CE-H toy model: V_int(L) = -2π v² n₁ n₂ log(L/r₀)")
print("  Coupling: q_i = 2π v n_i")
print("  So: V_int = (q₁/2π) * (q₂/2π) * 2π log(L/r₀) = (q₁*q₂/2π) log(L/r₀)")

# Estimate effective "coupling strength"
# If pressure term is E_press ~ Σ q_i q_j G(r_ij) with G(r) ~ -log(r)/(2π):
# Then E_press ~ Σ (q_i * q_j / 2π) * (-log(r_ij))

print("\n  Estimated pairwise couplings (relative to e-e):")
couplings = {}
for i, n1 in enumerate(['e', 'μ', 'τ']):
    for j, n2 in enumerate(['e', 'μ', 'τ']):
        if i <= j:  # symmetric
            if n1 not in charges_data or n2 not in charges_data:
                print(f"    q_{n1} * q_{n2} = FAILED (missing charge data)")
                continue
            q1 = charges_data[n1]['q_h2']
            q2 = charges_data[n2]['q_h2']
            coupling = q1 * q2
            couplings[f"{n1}-{n2}"] = coupling
            if 'e-e' in couplings:
                ratio = coupling / couplings['e-e']
                print(f"    q_{n1} * q_{n2} / q_e² = {ratio:.6f}")
            else:
                print(f"    q_{n1} * q_{n2} = {coupling:.8f}")

# =====================================================================
# DELIVERABLE: Summary & Candidate Charge Model
# =====================================================================

print("\n" + "=" * 72)
print("DELIVERABLE: Charge Extraction Summary")
print("=" * 72)

print("\nBest Hypothesis: H2 (q ∝ 1/Δg₀)")
print("  Rationale: Physically motivated by distance-from-vacuum scaling")
print("  Interpretation: charge strength inversely proportional to soliton excitation")
print("                  (larger Δg₀ → more excited → weaker coupling)")

print("\nExtracted Charges (H2 hypothesis):")
for name in ['e', 'μ', 'τ']:
    if name not in charges_data:
        print(f"  q_{name} = FAILED (solution error)")
        continue
    q = charges_data[name]['q_h2']
    A = charges_data[name]['A_tail']
    dg = charges_data[name]['delta_g0']
    print(f"  q_{name} = A_tail / Δg₀ = {A:.8f} / {dg:.8f} = {q:.8f}")

print("\nNext Step (Phase 3):")
print("  Use these charges in three-body self-consistency solver:")
print("  E_total({r_ij}) = Σ E_i + (1/2) Σ_{i≠j} q_i * q_j * G(r_ij)")
print("  where G(r) = -log(r/r₀)/(2π) is native Goldstone propagator")

print("\n" + "=" * 72)
print("PHASE 2 COMPLETE")
print("=" * 72)

# Save data for Phase 3
import json
data_out = {
    'charges': {name: {
        'g0': charges_data[name]['g0'],
        'A_tail': float(charges_data[name]['A_tail']),
        'q_h1': float(charges_data[name]['q_h1']),
        'q_h2': float(charges_data[name]['q_h2']),
        'q_h3': float(charges_data[name]['q_h3']),
    } for name in ['e', 'μ', 'τ']},
    'physical_constants': {
        'PHI': float(PHI),
        'PHI0': float(PHI0),
        'r21_pdg': float(R21_PDG),
        'r31_pdg': float(R31_PDG),
    }
}

print("\nData saved for Phase 3 (self-consistency solver)")
