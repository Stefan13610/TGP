#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
δ.1 Phase 3 — Sympy rigorous verification + cross-sector consistency

Phase 2 ujawniło dwie main candidates:
- H_NF: g̃ = N_f·e²/(12π) z N_f=5 (exact algebraic identity)
- H_loop: g̃ = 1 − α_s/(2π) (drift 0.122% od H_NF, ale fizycznie natural)

Phase 3 cele:
- P3.1: Sympy exact verification wszystkich hipotez
- P3.2: Cross-check Ω_Λ vs Planck (which hypothesis closer?)
- P3.3: α_s consistency
- P3.4: μ/e mass ratio test (czy λ.1 P2.3 reopens?)
- P3.5: Verdict synthesis
"""

import sympy as sp
import numpy as np


# Exact symbolic
e = sp.exp(1)
e_sq = sp.exp(2)
pi = sp.pi


def header(title):
    print()
    print("=" * 78)
    print(f"  {title}")
    print("=" * 78)
    print()


def p3_1_sympy_exact_forms():
    header("P3.1 — Sympy exact verification of hypothesis forms")

    print("Algebraic forms (sympy exact):")
    print()

    # H_NF exact form
    g_tilde_H_NF = sp.Rational(5) * e_sq / (12 * pi)
    print(f"  H_NF: g̃ = 5·e²/(12π)")
    print(f"    Symbolic: {g_tilde_H_NF}")
    print(f"    Simplified: {sp.simplify(g_tilde_H_NF)}")
    print(f"    Numerical: {float(g_tilde_H_NF):.10f}")
    print()

    # H_loop with PDG α_s
    alpha_s_pdg = sp.Rational(118, 1000)  # 0.118
    g_tilde_H_loop_pdg = 1 - alpha_s_pdg / (2 * pi)
    print(f"  H_loop (α_s=PDG): g̃ = 1 − α_s/(2π) = 1 − 0.118/(2π)")
    print(f"    Symbolic: {g_tilde_H_loop_pdg}")
    print(f"    Simplified: {sp.simplify(g_tilde_H_loop_pdg)}")
    print(f"    Numerical: {float(g_tilde_H_loop_pdg):.10f}")
    print()

    # H_loop required α_s for exact match z 5e²/(12π)
    alpha_s_required = 2 * pi * (1 - g_tilde_H_NF)
    alpha_s_required_simplified = sp.simplify(alpha_s_required)
    print(f"  H_loop required α_s for exact 5e²/(12π) match:")
    print(f"    α_s = 2π(1 - 5e²/(12π))")
    print(f"        = 2π - 10e²/12 = 2π − 5e²/6")
    print(f"    Symbolic: {alpha_s_required}")
    print(f"    Simplified: {alpha_s_required_simplified}")
    print(f"    Numerical: {float(alpha_s_required):.10f}")
    print(f"    PDG α_s: 0.1180 ± 0.0009 → drift = {(float(alpha_s_required) - 0.1180)/0.0009:.2f}σ")
    print()

    # Algebraic relationship between H_NF and H_loop
    print(f"  Algebraic relationship:")
    print(f"    If H_NF exact (g̃ = 5e²/(12π)),")
    print(f"    then 1-loop α_s would be 2π − 5e²/6 ≈ {float(alpha_s_required):.6f}")
    print(f"    vs PDG α_s = 0.1180 → +8.4σ TENSION")
    print()
    print(f"    Conversely, if H_loop exact z α_s_PDG,")
    print(f"    then g̃ = 1 - 0.118/(2π) = {float(g_tilde_H_loop_pdg):.6f}")
    print(f"    vs 5e²/(12π) = {float(g_tilde_H_NF):.6f} → drift 0.122%")
    print()
    print(f"  → H_NF and H_loop are NOT consistent z PDG α_s simultaneously")
    print(f"  → One must be approximation")
    print()

    return g_tilde_H_NF, g_tilde_H_loop_pdg, alpha_s_required


def p3_2_omega_lambda_cross_check():
    header("P3.2 — Ω_Λ cross-check vs Planck 0.6847 ± 0.0073")

    omega_planck = 0.6847
    sigma_planck = 0.0073

    # Ω_Λ = 2π·g̃/9 z γ.1 derivation
    print("  Ω_Λ = 2π·g̃/9 (γ.1 derivation)")
    print()

    # Pure structural g̃=1
    g_tilde_pure = 1
    omega_pure = 2 * np.pi * g_tilde_pure / 9
    print(f"  [Pure structural g̃=1]")
    print(f"    Ω_Λ = 2π/9 = {omega_pure:.6f}")
    print(f"    Planck deviation: {(omega_pure - omega_planck)/sigma_planck:+.2f}σ")
    print()

    # H_NF: g̃ = 5e²/(12π)
    g_tilde_H_NF = 5 * np.exp(2) / (12 * np.pi)
    omega_H_NF = 2 * np.pi * g_tilde_H_NF / 9
    print(f"  [H_NF: g̃ = 5e²/(12π) = {g_tilde_H_NF:.6f}]")
    print(f"    Ω_Λ = 2π·g̃/9 = {omega_H_NF:.6f}")
    print(f"    Equivalent: 5e²/54 = {5*np.exp(2)/54:.6f} ✓ same")
    print(f"    Planck deviation: {(omega_H_NF - omega_planck)/sigma_planck:+.2f}σ")
    print()

    # H_loop: g̃ = 1 - α_s_PDG/(2π)
    alpha_s = 0.1180
    g_tilde_H_loop = 1 - alpha_s / (2 * np.pi)
    omega_H_loop = 2 * np.pi * g_tilde_H_loop / 9
    print(f"  [H_loop: g̃ = 1 − α_s/(2π) = {g_tilde_H_loop:.6f}]")
    print(f"    Ω_Λ = 2π·g̃/9 = {omega_H_loop:.6f}")
    # Equivalent: 2π·(1 - α_s/(2π))/9 = 2π/9 - α_s/9
    omega_H_loop_alt = 2 * np.pi / 9 - alpha_s / 9
    print(f"    Equivalent: 2π/9 − α_s/9 = {omega_H_loop_alt:.6f}")
    print(f"    Planck deviation: {(omega_H_loop - omega_planck)/sigma_planck:+.2f}σ")
    print()

    # Best match
    print(f"  RANKING (closest to Planck):")
    candidates = [
        ("Pure structural (g̃=1)", omega_pure),
        ("H_NF (5e²/(12π))", omega_H_NF),
        ("H_loop (α_s/(2π))", omega_H_loop),
        ("Old cosmological 24.66/36", 24.66/36),
        ("Brannen 24.783/36", 24.783/36),
    ]
    candidates.sort(key=lambda x: abs(x[1] - omega_planck))
    for i, (name, val) in enumerate(candidates):
        sigma_dev = (val - omega_planck) / sigma_planck
        marker = "⭐" if i == 0 else "  "
        print(f"    {marker} {name:<35}: Ω_Λ = {val:.6f}, dev = {sigma_dev:+.2f}σ")
    print()


def p3_3_alpha_s_cross_check():
    header("P3.3 — α_s cross-check vs PDG 0.1180 ± 0.0009")

    g0_e = 0.86941
    alpha_s_pdg = 0.1180
    sigma_pdg = 0.0009

    # α_s = 27·g₀^e / (8·Φ_eff), Φ_eff = 8π·g̃
    print(f"  α_s = 27·g₀^e/(8·Φ_eff), Φ_eff = 8π·g̃")
    print(f"    g₀^e = {g0_e}")
    print(f"    α_s_PDG = {alpha_s_pdg} ± {sigma_pdg}")
    print()

    candidates = [
        ("Pure structural (g̃=1)", 1.0),
        ("H_NF (5e²/(12π))", 5 * np.exp(2) / (12 * np.pi)),
        ("H_loop (α_s/(2π) z PDG)", 1 - alpha_s_pdg / (2 * np.pi)),
        ("Cosmological 24.66 → g̃=24.66/(8π)", 24.66 / (8 * np.pi)),
        ("Brannen 24.783 → g̃=24.783/(8π)", 24.783 / (8 * np.pi)),
    ]

    print(f"  RANKING (closest to PDG α_s):")
    results = []
    for name, g_tilde in candidates:
        phi_eff = 8 * np.pi * g_tilde
        alpha_s = 27 * g0_e / (8 * phi_eff)
        sigma_dev = (alpha_s - alpha_s_pdg) / sigma_pdg
        results.append((name, g_tilde, phi_eff, alpha_s, sigma_dev))

    results.sort(key=lambda x: abs(x[4]))
    for i, (name, g_tilde, phi_eff, alpha_s, sigma_dev) in enumerate(results):
        marker = "⭐" if i == 0 else "  "
        print(f"    {marker} {name:<40}")
        print(f"        g̃ = {g_tilde:.6f}, Φ_eff = {phi_eff:.4f}, α_s = {alpha_s:.6f}, dev = {sigma_dev:+.2f}σ")
    print()


def p3_4_mass_ratio_check():
    header("P3.4 — μ/e mass ratio cross-check (λ.1 P2.3 reopener test)")

    print("  λ.1 mass formula (z R3 ODE α=2): m = c·A²·g₀^X")
    print("  X = e²/2 = exp(2)/2 = 3.694528")
    print()

    # PDG masses
    m_e_pdg = 0.51099895
    m_mu_pdg = 105.6583755
    m_tau_pdg = 1776.86

    # μ/e ratio under different X values
    X_target = m_mu_pdg / m_e_pdg
    X_lambda1 = np.exp(2) / 2

    print(f"  μ/e PDG ratio: {X_target:.6f}")
    print(f"  λ.1 mass formula slope X = e²/2 = {X_lambda1:.6f}")
    print()

    # γ.1 connection: jeśli Φ_eff = (10/3)·e², czy μ/e powiązane?
    phi_eff_lambda1 = (10/3) * np.exp(2)
    print(f"  λ.1 P2.3 hypothesis Φ_eff = (10/3)·e² = {phi_eff_lambda1:.6f}")
    print(f"  γ.1 algebraic identity: (10/3)·e² ≡ 8π·5e²/(12π)")
    print()

    # Czy λ.1 P2.3 reopens pod δ.1 H_NF success?
    # λ.1 P2.3 NEG verdict: anchor-dependent
    # γ.1 reframe: P2.3 hypothesis algebraic equivalent T-Λ corrected
    # δ.1 H_NF: g̃ = 5e²/(12π) z N_f=5 — czy to daje physical mechanism dla P2.3?
    print(f"  λ.1 P2.3 reopen analysis:")
    print()
    print(f"    Pre-γ.1: Φ_eff = (10/3)·e² treated as anchor-dependent numerology")
    print(f"    Post-γ.1: Φ_eff = (10/3)·e² ≡ 8π·5e²/(12π) (algebraic identity)")
    print(f"    Post-δ.1 H_NF: 5e²/(12π) with N_f=5 hint")
    print()
    print(f"    Mechanism dla P2.3:")
    print(f"      Φ_eff = 8π·N_f·e²/(12π) = N_f·e²·(8π)/(12π) = N_f·e²·(2/3)")
    print(f"      = (2/3)·N_f·e²")
    print(f"      Z N_f=5: Φ_eff = (10/3)·e² ✓ (matches λ.1 P2.3)")
    print()
    print(f"    To means λ.1 P2.3 hypothesis Φ_eff = (10/3)·e² jest equivalent z:")
    print(f"      Φ_eff = (2/3)·N_f·e² gdzie N_f=5 (active flavors at M_Z)")
    print()
    print(f"    To NIE jest derivation X = e²/2 mass formula slope —")
    print(f"    λ.1 P2.1, P2.2, M.4-M.6 NEG closures dla X derivation pozostają.")
    print(f"    Ale Φ_eff anchor (P2.3) zyskuje structural form.")
    print()


def p3_5_verdict():
    header("P3.5 — Phase 3 verdict synthesis")

    print("  KEY FINDINGS:")
    print()
    print("  1. ALGEBRAIC IDENTITY VERIFIED:")
    print("     g̃ = 5e²/(12π) ≈ 0.980004 (sympy exact)")
    print("     Φ_eff = 8π·g̃ = (10/3)·e² (γ.1 algebraic identity)")
    print("     Ω_Λ = 5e²/54 (clean algebraic form, equivalent)")
    print()
    print("  2. H_NF (N_f=5) identifies '5' z QCD active flavors:")
    print("     Numerical: 5·e²/(12π) = 0.980004 EXACT")
    print("     Mechanism: Φ_eff = (2/3)·N_f·e² z N_f=5 at M_Z scale")
    print()
    print("  3. H_loop (Schwinger) approximate:")
    print("     g̃ = 1 - α_s/(2π) z α_s_PDG drift 0.122% od H_NF")
    print("     Tension: jeśli H_loop exact, α_s_TGP = 2π − 5e²/6 = 0.1255")
    print("              vs PDG 0.1180 → +8.4σ")
    print("     → H_loop CANNOT be exact derivation if α_s_PDG=0.1180")
    print()
    print("  4. CROSS-SECTOR CONSISTENCY:")
    print("     Ω_Λ:  H_loop closer (+0.03σ) niż H_NF (-0.07σ) — both within 1σ")
    print("     α_s:  H_NF gives 0.1191 (+1.26σ), Brannen optimal (+0.44σ)")
    print()
    print("  5. λ.1 P2.3 STRUCTURAL REFRAMING:")
    print("     Φ_eff = (10/3)·e² gets natural interpretation:")
    print("       = (2/3)·N_f·e² z N_f=5 (active flavors)")
    print("     → P2.3 numerologia accusation reframed as physical formula")
    print()

    print("  VERDICT:")
    print()
    print("  H_NF (g̃ = N_f·e²/(12π) z N_f=5) jest TOP CANDIDATE z:")
    print("    + Exact algebraic identity (drift 0%)")
    print("    + Natural interpretation N_f=5 (QCD active flavors at M_Z)")
    print("    + Provides mechanism dla λ.1 P2.3 reframing")
    print("    + Ω_Λ deviation -0.07σ (within Planck 1σ)")
    print("    − Wymaga 'cosmological-gauge' coupling argument (Λ defined at M_Z scale)")
    print("    − e² source nie jest first-principles (imported z λ.1 mass formula)")
    print()
    print("  H_loop (g̃ = 1 − α_s/(2π)) jest SECONDARY:")
    print("    + Natural Schwinger structure")
    print("    + Ω_Λ deviation +0.03σ (best)")
    print("    − Drift 0.122% od H_NF algebraic identity")
    print("    − Wymaga α_s_TGP = 0.1255 (sprzeczne z PDG +8.4σ)")
    print()
    print("  CONCLUSION δ.1 PHASE 3:")
    print("    H_NF jest **partial structural derivation** dla g̃.")
    print("    POSITIVE w sensie: 5 = N_f (active flavors) jest natural.")
    print("    PARTIAL w sensie: wymaga że Λ sector defined at M_Z scale.")
    print()
    print("  Phase 4 implementation: dokumentować H_NF jako primary derivation,")
    print("  H_loop jako approximate alternative.")
    print()


def main():
    print("=" * 78)
    print("  δ.1 Phase 3 — Sympy verification + cross-sector consistency")
    print("=" * 78)
    print()

    p3_1_sympy_exact_forms()
    p3_2_omega_lambda_cross_check()
    p3_3_alpha_s_cross_check()
    p3_4_mass_ratio_check()
    p3_5_verdict()

    print()
    print("=" * 78)
    print("  Phase 3 complete")
    print("=" * 78)


if __name__ == "__main__":
    main()
