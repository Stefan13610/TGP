#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
γ.1 P3.5+ — Verify λ.1 P2.3 hypothesis (10/3)·e² ≈ 8π·g̃_T-Λ connection

Hypothesis discovered in Phase 3:
  - 8π · 0.98 = 24.6300
  - (10/3)·e² = 24.6302
  → g̃_T-Λ = 0.98003 ≈ (10/3)·e²/(8π)

Test: jest to numerical coincidence czy algebraic identity?
"""

import sympy as sp
import numpy as np


def main():
    print("=" * 78)
    print("  γ.1 P3.5+ — λ.1 P2.3 connection: (10/3)·e² ≈ 8π · g̃")
    print("=" * 78)
    print()

    # Sympy exact
    e_sq = sp.exp(2)
    eight_pi = 8 * sp.pi

    lambda1_phi_eff = sp.Rational(10, 3) * e_sq
    t_lambda_phi_eff_pure = eight_pi  # g̃=1
    t_lambda_g_tilde_required = lambda1_phi_eff / eight_pi

    print("Algebraic forms:")
    print(f"  λ.1 P2.3 hypothesis:  Φ_eff = (10/3)·e² = {lambda1_phi_eff}")
    print(f"  T-Λ pure structural:  Φ_eff = 8π = {eight_pi}")
    print()
    print(f"  Required g̃ for λ.1 ↔ T-Λ:  g̃ = (10/3)·e²/(8π) = {t_lambda_g_tilde_required}")
    print()

    # Numerical
    lambda1_num = float(lambda1_phi_eff)
    eight_pi_num = float(eight_pi)
    g_tilde_num = float(t_lambda_g_tilde_required)

    print(f"Numerical evaluation:")
    print(f"  (10/3)·e²        = {lambda1_num:.10f}")
    print(f"  8π               = {eight_pi_num:.10f}")
    print(f"  g̃ = (10·e²)/(24π) = {g_tilde_num:.10f}")
    print()

    # Compare to T-Λ reported g̃ = 0.98
    g_tilde_T_lambda = 0.98
    drift = abs(g_tilde_num - g_tilde_T_lambda) / g_tilde_T_lambda * 100
    print(f"Compare to T-Λ reported g̃:")
    print(f"  T-Λ closure quoted: g̃ = 0.98 (full M_Pl convention, exact match)")
    print(f"  λ.1↔T-Λ implied:    g̃ = {g_tilde_num:.6f}")
    print(f"  drift: {drift:.4f}%")
    print()

    # Check if g̃ exact = (10·e²)/(24π) is meaningful
    print("Hypothesis: g̃ structural = 5·e²/(12π) = (10/3)·e²/(8π)")
    print()
    g_tilde_alt = sp.Rational(5, 12) * e_sq / sp.pi
    print(f"  5·e²/(12π) = {sp.simplify(g_tilde_alt)} = {float(g_tilde_alt):.10f}")
    print(f"  Same as (10/3)·e²/(8π)? {sp.simplify(g_tilde_alt - t_lambda_g_tilde_required) == 0}")
    print()

    # Implication: pure T-Λ Ω_Λ_TGP and λ.1-corrected
    print("Implications:")
    print(f"  Ω_Λ pure T-Λ (g̃=1):     2π/9 = {float(2*sp.pi/9):.6f}")
    print(f"  Ω_Λ T-Λ (g̃=0.98):       2π·0.98/9 = {float(2*sp.pi*0.98/9):.6f}")
    print(f"  Ω_Λ z λ.1 hypothesis:    (10/3)·e²/36 = {float(lambda1_phi_eff/36):.6f}")
    print(f"  Ω_Λ Planck observed:     0.6847 ± 0.0073")
    print()

    # Test: what does g̃ = 5e²/(12π) mean?
    print("Algebraic interpretation g̃ = 5·e²/(12π):")
    print("  - π pochodzi z geometric ρ_crit = 3M_Pl²H₀²/(8π)")
    print("  - e² pochodzi z (10/3)·e² hipotezy λ.1 P2.3 (Brannen e_Euler² match)")
    print("  - 5/12 = 5/(12) — purely numerical factor")
    print()
    print("  Tested in λ.1: 'Φ_eff = (10/3)·e² ≈ 24.63' to 0.12% z cosmological 24.66")
    print("  γ.1 adds: (10/3)·e² = 8π·(5e²/12π) is CONSISTENT z T-Λ structural Φ_eff = 8π")
    print("           pod g̃ = 5e²/12π ≈ 0.98003, which is T-Λ's g̃ ≈ 0.98 fit")
    print()

    # Sigma comparison vs g̃=0.98003 vs T-Λ g̃=0.98
    omega_L_planck = 0.6847
    sigma_planck = 0.0073

    omega_L_lambda1 = float(lambda1_phi_eff / 36)
    omega_L_T_lambda_g098 = float(2 * sp.pi * 0.98 / 9)
    omega_L_T_lambda_pure = float(2 * sp.pi / 9)

    print(f"Ω_Λ values vs Planck 0.6847 ± 0.0073:")
    candidates = [
        ("T-Λ pure (g̃=1)", omega_L_T_lambda_pure),
        ("T-Λ g̃=0.98", omega_L_T_lambda_g098),
        ("λ.1 (10/3)·e²/36", omega_L_lambda1),
        ("Brannen 24.783/36", 24.783/36),
        ("Old cosmological 24.66/36", 24.66/36),
    ]
    for name, val in candidates:
        sigma = (val - omega_L_planck) / sigma_planck
        print(f"  {name:<35}: Ω_Λ = {val:.6f}, deviation = {sigma:+.2f}σ")
    print()

    # α_s under each
    g0_e = 0.86941
    print(f"α_s under each Φ_eff (formula 27·g₀^e/(8·Φ_eff)):")
    phi_eff_candidates = [
        ("T-Λ pure 8π", float(eight_pi)),
        ("T-Λ g̃=0.98 (8π·0.98)", float(eight_pi) * 0.98),
        ("λ.1 (10/3)·e²", lambda1_num),
        ("Brannen 24.783", 24.783),
        ("Old cosmological 24.66", 24.66),
    ]
    pdg_alpha_s = 0.1180
    pdg_sigma = 0.0009
    for name, phi in phi_eff_candidates:
        alpha = 27 * g0_e / (8 * phi)
        sig = (alpha - pdg_alpha_s) / pdg_sigma
        print(f"  {name:<35}: Φ_eff = {phi:.4f}, α_s = {alpha:.6f}, dev = {sig:+.2f}σ")
    print()

    # Summary
    print("=" * 78)
    print("  Final synthesis γ.1 H5 (Φ_eff = 8π, T-Λ structural)")
    print("=" * 78)
    print()
    print("  Pure structural prediction (g̃=1):")
    print(f"    Φ_eff = 8π = {eight_pi_num:.4f}")
    print(f"    Ω_Λ_TGP = 2π/9 = {float(2*sp.pi/9):.4f}")
    print(f"    α_s = 27·g₀^e/(8·8π) = {27*g0_e/(8*eight_pi_num):.6f}")
    print(f"    Tensions: Ω_Λ +1.84σ Planck, α_s -1.39σ PDG")
    print()
    print("  T-Λ-corrected (g̃ = 0.98 fit):")
    print(f"    g̃ ≈ 0.98 = (10/3)·e²/(8π) = 5·e²/(12π)")
    print(f"    Φ_eff_eff ≈ {eight_pi_num*0.98:.4f} ≈ (10/3)·e² = {lambda1_num:.4f}")
    print(f"    Ω_Λ_eff = 0.98·(2π/9) = {2*np.pi*0.98/9:.4f}")
    print(f"    α_s under = {27*g0_e/(8*eight_pi_num*0.98):.6f}")
    print()
    print("  λ.1 P2.3 hypothesis (10/3)·e² ≈ Φ_eff:")
    print(f"    Identical to T-Λ-corrected within 0.001%")
    print(f"    (10/3)·e² = 8π · (5e²/12π) = 8π·0.98003")
    print()
    print("  KLUCZOWE: λ.1 i T-Λ formalisms są ALGEBRAICZNIE consistent")
    print("            via g̃ = 5e²/(12π).")
    print()
    print("=" * 78)


if __name__ == "__main__":
    main()
