---
title: "ψ.1.Phase3 setup — lab predictions + 4-channel convergence (6 sub-tests)"
date: 2026-05-01
cycle: ψ.1.Phase3
parent: "[[program.md]]"
status: setup
tags:
  - TGP
  - psi1
  - phase3
  - predictions
  - setup
---

# ψ.1.Phase3 — Lab predictions + 4-channel convergence

## Sub-tests (6)

### T3.1 — Sagnac fazowy lab E∥B prediction (WYKONALNE DZIŚ)
- Setup: Mach-Zehnder lub Sagnac fazowy interferometer
- Beam #1 (signal): pass through E∥B Schwinger-class region (or ELI-NP class), L = 10 cm
- Beam #2 (reference): identical geometry, E⊥B (null source) for chopper
- Predicted Δφ ~ 6×10⁻⁷ rad (β_g = 1, Λ = 100 MeV), 4σ in single shot at LIGO-class 10⁻¹¹ rad
- Falsification: null at LIGO-class precision z Schwinger-class E∥B → β_g ~ 0 OR Λ > 100 MeV
- PASS if SNR ≥ 4σ wykonalne dziś (≥1σ minimum)

### T3.2 — TOF dual-arm zs-precision frontier 2030+
- Setup: 2 synchronized photon paths (E∥B vs E⊥B), L = 10 cm
- TOF measurement via attoclock 2030+ technology (1 zs precision)
- Predicted Δt ~ 0.3 zs (β_g = 1, Λ = 100 MeV)
- SNR ~10⁰ at 2030+, ~10² at 2035+ frontier
- Alternative: nonlinear phase-shift mixing in BBO crystals (higher harmonic generation timing)
- PASS if 2030+ TOF wykonalne ≥1σ

### T3.3 — Cosmological scalar c shift residual NULL (consistent z Webb/Murphy)
- Cosmological FRW: Hubble flow H₀ ~ 10⁻¹⁸ s⁻¹ → (∂lnX)_cosmo ~ H₀ ~ 10⁻¹⁸/c
- Δc/c_cosmo ~ |β_g| · 10⁻³⁶/Λ² ~ |β_g| · 10⁻⁶⁰ (Λ = 100 MeV) → << any sensitivity
- Webb/Murphy 1e-7 NULL on Δα/α: dependence Δα/α + Δ(c/c) z multiplet not separable; ψ.1 contributes Δc/c << 10⁻¹⁵
- Status: NULL at all current cosmological sensitivities — consistent
- PASS if cosmological residual << current sensitivity (consistent, no false positive)

### T3.4 — Magnetar atmosphere FRB time-of-flight (NOVEL)
- FRB pulses passing through magnetar plasma: B ~ 10¹⁵ G, locally E∥B from twist
- Path-integrated (∂lnX)² along line of sight → cumulative Δt
- For SGR 1935+2154 FRB 200428 (galactic magnetar burst):
  - DM = 332.7 pc·cm⁻³ (cumulative dispersion measure)
  - if 1% path through E∥B regions B ~10¹⁵ G E ~10¹⁰ V/m, Δt ~ μs at 1.4 GHz
  - signature: frequency-dependence Δt ∝ 1/ω² CONTRAST with standard plasma DM ∝ 1/ω²
  - ψ.1 predicts ω-INDEPENDENT shift (scalar c) — DISCRIMINATOR vs DM
- Falsification: FRB array (CHIME, ASKAP) timing residual ω-independent at sensitivity 1 ms
- PASS if ω-independent FRB shift detectable at CHIME/ASKAP precision

### T3.5 — 4 alt-L₅ couplings cross-channel falsification pattern
- Joint signature pattern (Sagnac × TOF × cosmo × FRB) identifies which L₅ form realized:
  | Form | Sagnac | TOF | Cosmo | FRB |
  |------|--------|-----|-------|-----|
  | L₅_a (∂lnX)²·F² | scalar shift | scalar Δt | NULL | ω-independent |
  | L₅_b (∂lnX)²·F·F̃ | helicity Δφ | helicity Δt | helical PMF | ω-, helicity- correlated |
  | L₅_c (□lnX)·F² | reduces — equivalent to L₅_a | parts-equivalent | parts | parts |
  | L₅_d (∂lnX)²·(E²-B²)² | quadratic | quadratic | suppressed | quadratic ω scaling |
- Joint pattern uniquely identifies L₅ realized
- PASS if 4 forms cross-distinguished

### T3.6 — 4-channel ψ.1 convergence
- ✓ Channel 1 LAB Sagnac: SNR ~10⁴ wykonalne dziś, Λ ≤ 10 GeV reachable
- ✓ Channel 2 LAB TOF: SNR ~10⁰ at 2030+, Λ extension to ~100 GeV
- ✓ Channel 3 COSMO: NULL consistent (B too small)
- ✓ Channel 4 ASTRO FRB: ω-independent residual at CHIME/ASKAP 1 ms precision
- 4-channel target: ≥3 channels active
- PASS if 4/4 convergence

## Phase verdict criterion

6/6 PASS → ψ.1 program END.
5/6 → final patch + END projection.
≤4/6 → ψ.1 STALL revisit.

## Master/phase ledger projection

ψ.1.Phase1 5 + ψ.1.Phase2 7 + ψ.1.Phase3 6 = **18 sub-tests**
- Master ledger: 751 → 769
- Phase ledger: 503 → 521
