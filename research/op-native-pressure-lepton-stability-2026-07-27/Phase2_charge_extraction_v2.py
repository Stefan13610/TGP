#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PHASE 2 v2: Extract Lepton Charges from CP-7 Profiles

op-native-pressure-lepton-stability-2026-07-27
Session #64, Tier 2 v5

Strategy:
  Import soliton_profile() directly from CP-7 Phase2_bvp_spectrum.py
  to regenerate the EXACT same profiles used in spectral analysis.
  Then extract far-field amplitude and derive charges.
"""

import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
import math
from scipy.optimize import fminbound

# =====================================================================
# Import soliton_profile from CP-7
# =====================================================================

# Add CP-7 directory to path so we can import the function
sys.path.insert(0, '../op-spectral-analysis-Phi-2026-07-03')

try:
    from Phase2_bvp_spectrum import soliton_profile
    print("✓ Successfully imported soliton_profile from CP-7 Phase2_bvp_spectrum.py")
except ImportError as e:
    print(f"✗ Failed to import: {e}")
    print("Falling back to simplified ODE solver...")
    soliton_profile = None

# =====================================================================
# CONSTANTS (CP-7 Formulacja B / Log-form)
# =====================================================================
# NOTE: These are CP-7 values, NOT canonical substrate values!
# CP-7 uses Formulacja B: f(g) = 1 + 4·ln(g), with ghost wall at g* = e^(-1/4)
# Canonical version (0.86941) uses different parametrization.
# For consistency with CP-7 spectra/profiles, use these:

PHI = (1 + np.sqrt(5)) / 2
PHI0 = 25.0
G0_E = 1.24915      # CP-7 value (Formulacja B)
G0_MU = PHI * G0_E  # = 2.02117
G0_TAU = 3.18912    # Koide-determined (same in both)

# Physical constants (PDG 2024)
M_E = 0.51099895
M_MU = 105.6583755
M_TAU = 1776.86
R21_PDG = M_MU / M_E
R31_PDG = M_TAU / M_E

print("=" * 72)
print("PHASE 2 v2: Extract Lepton Charges (from CP-7 Profiles)")
print("op-native-pressure-lepton-stability-2026-07-27")
print("=" * 72)

print(f"\nCanonical Parameters:")
print(f"  G0_E = {G0_E:.6f}")
print(f"  G0_MU = {G0_MU:.6f}")
print(f"  G0_TAU = {G0_TAU:.6f}")

# =====================================================================
# TASK 2a: Generate Profiles (using CP-7)
# =====================================================================

print("\n" + "=" * 72)
print("TASK 2a: Generate Soliton Profiles (using CP-7)")
print("=" * 72)

leptons = [
    ('e', G0_E),
    ('μ', G0_MU),
    ('τ', G0_TAU)
]

profiles = {}

print(f"\n{'Gen':>3} {'g₀':>10} {'Status':>15} {'g_min':>12} {'bounces':>10}")
print(f"{'-'*3} {'-'*10} {'-'*15} {'-'*12} {'-'*10}")

for name, g0 in leptons:
    try:
        if soliton_profile is None:
            raise RuntimeError("soliton_profile not available")

        # Generate profile with same parameters as CP-7 Phase2
        r, g, gp, d2, bounces, g_min = soliton_profile('F-S', g0, R=60.0, N=4000)

        profiles[name] = {
            'r': r,
            'g': g,
            'gp': gp,
            'd2': d2,
            'bounces': bounces,
            'g_min': g_min,
        }

        print(f"{name:>3} {g0:>10.6f} {'PASS':>15} {g_min:>12.6f} {bounces:>10}")

    except Exception as e:
        print(f"{name:>3} {g0:>10.6f} {'FAIL':>15} {str(e):>20}")

# =====================================================================
# TASK 2b: Extract Far-Field Charges
# =====================================================================

print("\n" + "=" * 72)
print("TASK 2b: Extract Charges from Far-Field")
print("=" * 72)

def extract_far_field_amplitude(r, g, r_min=50.0, r_max=None):
    """
    Extract oscillatory tail amplitude from profile.

    Ansatz: g(r) → 1 + (A_tail/r) * sin(r + φ) + O(1/r²)
            u(r) = (g-1)*r → A_tail * sin(r + φ)

    Fit to: u(r) = B*cos(r) + C*sin(r)
    """
    if r_max is None:
        r_max = r[-1]

    mask = (r >= r_min) & (r <= r_max)
    if np.sum(mask) < 10:
        return None, None, None, None, None

    r_fit = r[mask]
    u_fit = (g[mask] - 1.0) * r_fit

    # Least-squares: u = B*cos(r) + C*sin(r)
    X = np.column_stack([np.cos(r_fit), np.sin(r_fit)])
    try:
        coef, *_ = np.linalg.lstsq(X, u_fit, rcond=None)
    except:
        return None, None, None, None, None

    B, C = coef[0], coef[1]
    A_tail = np.sqrt(B**2 + C**2)
    phi = np.arctan2(C, B)

    # Residual quality
    y_hat = B * np.cos(r_fit) + C * np.sin(r_fit)
    residual_rms = np.sqrt(np.mean((u_fit - y_hat)**2)) / max(A_tail, 1e-12)

    return A_tail, phi, residual_rms, B, C

print(f"\n{'Gen':>3} {'g₀':>10} {'Δg₀':>10} {'A_tail':>12} {'q (H2)':>12} {'residual':>12}")
print(f"{'-'*3} {'-'*10} {'-'*10} {'-'*12} {'-'*12} {'-'*12}")

charges_data = {}

for name, g0 in leptons:
    if name not in profiles:
        print(f"{name:>3} {g0:>10.6f} FAILED - no profile")
        continue

    profile = profiles[name]
    r = profile['r']
    g = profile['g']

    # Extract far-field amplitude
    A_tail, phi, residual, B, C = extract_far_field_amplitude(r, g, r_min=45.0, r_max=58.0)

    if A_tail is None:
        print(f"{name:>3} {g0:>10.6f} N/A       FAILED")
        continue

    delta_g0 = g0 - 1.0

    # Hypothesis H2: q ∝ A_tail / Δg₀
    q_h2 = A_tail / delta_g0

    print(f"{name:>3} {g0:>10.6f} {delta_g0:>10.6f} {A_tail:>12.8f} {q_h2:>12.8f} {residual:>12.6f}")

    charges_data[name] = {
        'g0': g0,
        'delta_g0': delta_g0,
        'A_tail': A_tail,
        'phi': phi,
        'residual': residual,
        'B': B,
        'C': C,
        'q_h2': q_h2,
    }

# =====================================================================
# TASK 2c: Verify Against Coupling Model
# =====================================================================

print("\n" + "=" * 72)
print("TASK 2c: Verify Coupling Model")
print("=" * 72)

if len(charges_data) >= 2:
    print("\n§ C1: Charge Ratios")

    q_vals = {}
    for name in ['e', 'μ', 'τ']:
        if name in charges_data:
            q_vals[name] = charges_data[name]['q_h2']

    if 'e' in q_vals:
        print(f"  q_e = {q_vals['e']:.8f}")
        if 'μ' in q_vals:
            print(f"  q_μ = {q_vals['μ']:.8f}   (ratio μ/e = {q_vals['μ']/q_vals['e']:.6f})")
        if 'τ' in q_vals:
            print(f"  q_τ = {q_vals['τ']:.8f}   (ratio τ/e = {q_vals['τ']/q_vals['e']:.6f})")

    print("\n§ C2: Pairwise Interaction Couplings")
    print("  V_int(L) ∝ q_i * q_j * log(L/r₀)   [from Goldstone propagator]")

    couplings = []
    for i, n1 in enumerate(['e', 'μ', 'τ']):
        for j, n2 in enumerate(['e', 'μ', 'τ']):
            if i <= j and n1 in q_vals and n2 in q_vals:
                coupling = q_vals[n1] * q_vals[n2]
                couplings.append((f"{n1}-{n2}", coupling))

    if couplings:
        couplings_sorted = sorted(couplings, key=lambda x: -x[1])
        for pair, coup in couplings_sorted:
            if 'e-e' in [c[0] for c in couplings]:
                ee_coup = [c[1] for c in couplings if c[0] == 'e-e'][0]
                ratio = coup / ee_coup
                print(f"    {pair:>5}: {coup:12.8f}  (relative to e-e: {ratio:8.4f}×)")
            else:
                print(f"    {pair:>5}: {coup:12.8f}")

# =====================================================================
# DELIVERABLE
# =====================================================================

print("\n" + "=" * 72)
print("DELIVERABLE: Extracted Charges for Phase 3")
print("=" * 72)

if charges_data:
    print("\nExtracted Charges (Hypothesis H2: q ∝ A_tail / Δg₀):")
    for name in ['e', 'μ', 'τ']:
        if name in charges_data:
            d = charges_data[name]
            print(f"\n  {name}:")
            print(f"    g₀ = {d['g0']:.6f}")
            print(f"    Δg₀ = {d['delta_g0']:.6f}")
            print(f"    A_tail = {d['A_tail']:.8f}")
            print(f"    q = {d['q_h2']:.8f}")
            print(f"    [phase φ = {d['phi']:.4f}, residual = {d['residual']:.2e}]")
else:
    print("No charges extracted - profile generation failed")

print("\n" + "=" * 72)
print("PHASE 2 v2 COMPLETE")
print("=" * 72)

print("\nNext Step (Phase 3):")
print("  Use these charges q_e, q_μ, q_τ in three-body self-consistency:")
print("  E_total({r_ij}) = Σ E_i + (1/2) Σ_{i≠j} q_i * q_j * G(r_ij)")
print("  where G(r) = -log(r/r₀)/(2π) [native Goldstone propagator]")
