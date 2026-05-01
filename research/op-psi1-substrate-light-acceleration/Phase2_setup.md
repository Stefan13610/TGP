---
title: "ψ.1.Phase2 setup — sympy LOCK Δc/c + Sagnac/TOF engineering (7 sub-tests)"
date: 2026-05-01
cycle: ψ.1.Phase2
parent: "[[program.md]]"
status: setup
tags:
  - TGP
  - psi1
  - phase2
  - sympy
  - setup
---

# ψ.1.Phase2 — sympy LOCK Δc/c + lab Sagnac/TOF engineering

## Sub-tests (7)

### T2.1 — sympy LOCK Δc/c formula
- Symbolic derivation: from L_em + L₅, derive c_local
- Target formula: $\Delta c/c_0 = -(\beta_g/(2\Lambda^2))\,(\partial \ln X)^2$
- Verify sympy diff(target − derived) = 0 EXACT
- PASS if sympy LOCK exact

### T2.2 — F² = 2(B²-E²) vs F·F̃ = -4 E·B distinction
- ψ.1 uses F² (in L₅) but ω.1 source uses F·F̃ (∂lnX from E·B)
- Lab E∥B: maximizes F·F̃ → sources (∂lnX) ≠ 0; F² = 2(B²-E²) independent
- Cross-coupling: (∂lnX) sourced by E·B then dotted into F² in L₅
- Verify both expressions sympy-LOCKED + numerical at Schwinger-class
- PASS if both formulas consistent

### T2.3 — Yukawa Greens for substrate gradient
- Substrate EOM with mass: (□ + m_X²)(lnX) = source
- Greens function G(r) = -e^(-m_X·r)/(4πr)
- Two regimes:
  - Light substrate (m_X⁻¹ >> L_field): ∂lnX ~ ρ·L/(4π) — field-region dominated
  - Heavy substrate (m_X⁻¹ << L_field): ∂lnX ~ ρ/(4π m_X²) — localized
- Verify Greens function structure
- PASS if both regimes consistent + scaling laws derived

### T2.4 — Λ-cutoff scan
- Scan Λ ∈ {M_Pl, TeV, GeV, 100 MeV, 10 MeV, 1 MeV}
- Compute Δc/c, Sagnac Δφ, TOF Δt for each
- Compare to current sensitivities (Sagnac LIGO 10⁻¹¹ rad, TOF as 10⁻¹⁸ s)
- Identify viability window
- PASS if Λ ≲ 10 GeV detectable today via Sagnac

### T2.5 — Sagnac fazowy SNR + TOF zs SNR numerical at Schwinger-class
- Sagnac with L = 10 cm, λ = 1064 nm:
  - signal Δφ = 6×10⁻⁷ rad at Λ = 100 MeV
  - LIGO-class noise floor: 10⁻¹¹ rad after integration
  - SNR ~ 10⁴ (4σ in single shot)
- TOF with L = 10 cm:
  - signal Δt = 0.3 zs
  - attoclock 2030+ noise floor: 1 zs
  - SNR ~ 10⁻³ today, ~10⁰ at 2030+
- PASS if Sagnac wykonalne dziś + TOF frontier 2030+

### T2.6 — Cross-coupling z σ.1 (sub-leading scalar vs leading helicity-dependent)
- Build consistency matrix:
  | Order | Channel | Status |
  |-------|---------|--------|
  | O((∂lnX)¹) | σ.1 helicity birefringence | leading, axion-like F·F̃ |
  | O((∂lnX)²) scalar | **ψ.1 L₅ scalar c shift** | sub-leading SOURCEABLE LAB |
  | O((∂lnX)²) helicity | sub-leading axion correction | suppressed |
  | O((∂lnX)≥3) | higher EFT | further suppressed |
- ψ.1 IS the sub-leading σ.1 scalar channel — NO TENSION
- PASS if consistency matrix complete + no leading scalar c(X) tension

### T2.7 — 4 alt-L₅ couplings cross-falsification
- 4 candidate L₅ forms (a/b/c/d) tested via E∥B vs E⊥B + sign-flip + pure-field:
  | Form | E∥B | E⊥B | Sign-flip E·B | Pure E | Pure B |
  |------|-----|-----|---------------|--------|--------|
  | L₅_a (∂lnX)²·F² (CANONICAL) | signal scalar | null | same | null | null |
  | L₅_b (∂lnX)²·F·F̃ | helicity-signal | null | FLIPS | null | null |
  | L₅_c (□lnX)·F² | reduces by parts | parts-equiv | parts-equiv | parts | parts |
  | L₅_d (∂lnX)²·(E²-B²)² | F² inherent + extra | extra null | same | extra E | extra B |
- Two-axis differential test (helicity × parity × field-config) uniquely identifies form
- PASS if 4 forms experimentally distinguishable

## Phase verdict criterion

7/7 PASS → ψ.1.Phase 3 forward.
6/7 → minor patch + Phase 3 advance.
≤5/7 → ψ.1 STALL.
