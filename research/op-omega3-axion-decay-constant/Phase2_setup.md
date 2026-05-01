---
title: "ω.3.Phase2 setup — sympy LOCK + cosmological + F-cluster"
date: 2026-05-01
cycle: ω.3.Phase2
status: SETUP
parent: "[[program.md]]"
predecessor: "[[Phase1_results.md]]"
tags:
  - TGP
  - omega3
  - phase2
  - sympy-lock
  - f_a
  - setup
---

# ω.3.Phase2 setup

**Score gate:** ≥6/7 PASS = Phase 3 ENABLED.

## Sub-tests (7)

- **O2.1** **f_a sympy LOCK rigorous** = (N_A·2π²·M_GUT)/(536/75) z **closed-form
  sympy** (Rational, pi, no float).

- **O2.2** **f_a numerical reproduction** = 4.847·10¹⁷ GeV via UV.2 M_TGP × 1/E_TGP.
  - **Drift gate**: < 1% vs UV.2 K-LOCK trace.

- **O2.3** **g_aγ numerical via direct + via inversion**:
  - direct: g_aγ = α_em·E_TGP/(2π·f_a)
  - inversion: g_aγ = g_axion(ω.2)/f_a
  - sympy diff = 0 EXACT requirement.

- **O2.4** **Cosmological misalignment / ALP DM consistency**:
  - TGP ALP m_a free → no specific Ω_a h² band prediction
  - if pre-inflation PQ: misalignment Ω_a h² ∝ θ_i²·(f_a)^(7/6) for QCD-axion
    formula reference; for ALP, mechanism depends on m_a (free)
  - **Gate**: f_a in super-GUT regime CONSISTENT z anthropic θ_i misalignment scenario.

- **O2.5** **Isocurvature pre-inflation PQ constraint**:
  - if PQ broken before inflation + θ_i ~ O(1): Δ_iso ~ H_inf/(2π·f_a·θ_i)
  - Planck PR4 limit: Δ_iso < 10⁻¹¹ at k = 0.05 Mpc⁻¹
  - For TGP f_a = 4.85·10¹⁷ + θ_i = 1: H_inf < 2π·f_a·10⁻¹¹ ≈ 3·10⁷ GeV
  - **Gate**: TGP f_a allows H_inf up to ~10⁷ GeV (compatible z standard slow-roll inflation r ~ few × 10⁻¹⁰).

- **O2.6** **F-cluster post-ω.3 preservation**:
  - F4 (α₀ = 1069833/264500) untouched ✓ (ω.3 nie ingeruje)
  - F5 (g̃ = 0.9803) untouched ✓
  - F6 (κ = √(32π) DERIVED post-χ.1) untouched ✓
  - XS1 (√α₀ = κ_TGP drift 0.042%) untouched ✓
  - **Gate**: max F-cluster drift < 0.1% post-ω.3.

- **O2.7** **4-channel ω.3 cascade self-consistency**:
  - Channel 1: g* (UV.1 NGFP) = 71/100 ✓
  - Channel 2: N_A (ξ.1 photon-ring) = 500/57 ✓
  - Channel 3: K_struct = N_A·2π² (UV.2) ≈ 173.15 ✓
  - Channel 4: E_TGP (ω.2) = 536/75 ✓
  - All 4 anchors flow into f_a = (g*/...)·(N_A·2π²·M_GUT)/(E_TGP); sympy diff = 0 across cascade.

## Inputs (LOCKED)

```
g* (UV.1)            = 71/100
N_A (xi.1)           = 500/57
E_TGP (omega.2)      = 536/75
K_struct (UV.2)      = N_A·2pi^2 ≈ 173.15
M_TGP (UV.2)         = K·M_GUT ≈ 3.4630·10^18 GeV
M_GUT (input)        = 2.0·10^16 GeV (SM 2-loop)
g_axion (omega.2)    = alpha_em·E_TGP/(2pi) ≈ 8.30·10^-3
alpha_em             = 1/137.036
```

## Strategy

Phase 2 LOCKS f_a = M_TGP/E_TGP sympy + numerical, validates cosmological band
(misalignment + isocurvature), confirms F-cluster preservation, and verifies
4-channel cascade self-consistency. Phase 3 will then expand to predictions +
forward-gates.

## Cross-references

- [[program.md]]
- [[Phase1_results.md]]
