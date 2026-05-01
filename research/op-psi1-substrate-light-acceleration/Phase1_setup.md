---
title: "ψ.1.Phase1 setup — L₅ coupling structural derivation (5 sub-tests)"
date: 2026-05-01
cycle: ψ.1.Phase1
parent: "[[program.md]]"
status: setup
tags:
  - TGP
  - psi1
  - phase1
  - setup
---

# ψ.1.Phase1 — L₅ coupling structural derivation

## Sub-tests (5)

### T1.1 — L₅ candidate scan + φ.1 X→λX scale-invariance
- Construct 4 candidate dim-6 EFT operators coupling (∂lnX) to F²:
  - L₅_a: (∂lnX)²·F² (CANONICAL — scalar, derivative-only)
  - L₅_b: (∂lnX)²·F·F̃ (parity-odd alternative)
  - L₅_c: (□lnX)·F² (single-derivative alternative — fails by parts)
  - L₅_d: lnX·F² (massless dilaton — VIOLATES φ.1 scale-symmetry)
- Verify φ.1 invariance: only L₅_a survives (and L₅_b for axion-like, but breaks parity in scalar c)
- PASS if L₅_a uniquely scale-invariant scalar coupling

### T1.2 — UV matching β_g sign (3 independent channels)
- **Channel A:** AS NGFP fixed-point Wilson coefficient analysis (Reuter+ 2002)
- **Channel B:** Heavy-mode integration: substrate-coupled fermion loop one-loop estimate β_g ~ Σ Q_f²·m_f²/(16π²Λ²)
- **Channel C:** Cosmological consistency: BBN photon spectrum at z~10⁹ requires |β_g·(∂lnX)²/Λ²| < 10⁻⁴
- 3-channel sign agreement → β_g sign definitive
- PASS if all 3 channels yield consistent sign (anticipated β_g < 0 → CLOCK ACCELERATES light)

### T1.3 — Effective scalar c shift derivation
- Photon EOM with L₅: ∂_ν[(1 + β_g(∂lnX)²/Λ²) F^νμ] = 0
- Plane-wave ansatz: ω² = c_local² · k², solve for c_local
- Derive: c_local = c₀/√(1 + β_g(∂lnX)²/Λ²) ≈ c₀(1 - β_g(∂lnX)²/(2Λ²))
- PASS if scalar c shift formula consistent with EOM

### T1.4 — ω.1 EOM source maximization F·F̃ via E∥B parallel + null controls
- Use τ.3 inheritance: ω.1 EOM □(lnX) = -(g/f_X²) E·B sources (∂lnX) ≠ 0 only when F·F̃ ≠ 0
- F·F̃ = -4 E·B parallel maximization (cos θ = 1)
- Configurations:
  - E∥B parallel: max source ✓
  - E⊥B perpendicular: F·F̃ = 0, NO source — null CONTROL ✓
  - Pure E (B=0): F·F̃ = 0, null CONTROL ✓
  - Pure B (E=0): F·F̃ = 0, null CONTROL ✓
- PASS if all configurations match expected (∂lnX)² ≠ 0 vs = 0

### T1.5 — Viability gate Λ ≲ 10 GeV via Sagnac LIGO-class today
- Schwinger-class lab E·B parallel: |E·B| ~ 10¹⁸ V·T/m
- (∂lnX)² ~ ((g/f_X²) E·B · L_eff)² for Yukawa range L_eff ~ 1 µm at m_X = 100 MeV
- Δc/c ~ |β_g|·(∂lnX)²/Λ² ~ |β_g|·10⁻¹² at Λ = 100 MeV
- Sagnac fazowy: Δφ = ω·L·Δc/c₀ ~ 10¹⁵·0.1·10⁻¹²/c₀ ~ 6×10⁻⁷ rad
- LIGO-class Sagnac sensitivity: 10⁻¹¹ rad → SNR ~10⁴
- PASS if Λ ≲ 10 GeV (window 1000× wider than τ.3 100 MeV bound)

## Phase verdict criterion

5/5 PASS → ψ.1.Phase 2 forward.
4/5 → re-examine alt-couplings.
≤3/5 → ψ.1 STALL.

## Cross-cycle inheritance

- φ.1: substrate scale-symmetry X→λX axiom (L₅_a derivative-only form)
- ω.1: photon-substrate axion coupling EOM □(lnX) = -(g/f_X²) E·B
- σ.1: leading helicity-dependent dispersion (NO scalar c at leading)
- τ.3: L₄ analogous mass channel — ψ.1 is L₅ analog for photon
