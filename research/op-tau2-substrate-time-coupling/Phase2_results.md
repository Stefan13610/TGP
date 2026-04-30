---
title: "τ.2.Phase2 results — proper time + sympy LOCK 7/7 PASS"
date: 2026-04-30
cycle: τ.2.Phase2
status: PASS
parent: "[[program.md]]"
tags:
  - TGP
  - tau2
  - phase2
  - sympy
  - results
---

# τ.2.Phase2 results — 7/7 FULL CASCADE

## Sub-test outcomes

| ID | Test | Result |
|----|------|--------|
| **T2.1** | Proper time τ = ∫√g_00^eff dt sympy LOCK | ✅ PASS |
| **T2.2** | Substrate-induced time dilation O(∂ ln X) = 0 | ✅ PASS |
| **T2.3** | Polarization-Zeeman cross-coupling z σ.1 | ✅ PASS |
| **T2.4** | Cs hyperfine (9.2 GHz) vs optical Sr/Yb (430 THz) protection | ✅ PASS |
| **T2.5** | Clock-rate ratio R(X1)/R(X2) = 1 sympy EXACT | ✅ PASS |
| **T2.6** | ω.1+σ.1 consistency (scalar c NULL → scalar ω NULL) | ✅ PASS |
| **T2.7** | Alt-couplings m_e ∝ X^α FALSIFIED | ✅ PASS |

**Score: 7/7 → Phase 3 forward**

## Key sympy LOCK results

### T2.1: Proper time formal definition
- Effective metric: g_00^eff = 1 + ε² where ε = (∂ ln X)/Λ.
- Proper time: τ = ∫√(g_00^eff) dt.
- Series expansion: dτ/dt = 1 + ε²/2 + O(ε⁴).
- **Linear coefficient = 0** (sympy verified) ⟹ NO O(∂ ln X) substrate time-dilation.

### T2.2: O(∂ ln X) = 0 confirmed
δτ/τ = 0 at leading O(∂ ln X). First nonzero correction at O(ε²) = O((∂ ln X / Λ)²) ~ Λ-suppressed dim-6 EFT.

### T2.5: Clock-rate ratio sympy EXACT
For ω_atom ~ m_e c² α_em²/(ℏ n²) under X→λX:
- m_e → m_e λ^α_m
- α_em → α_em λ^α_a
- ω_atom → ω_atom λ^(α_m + 2α_a)

**Demand R(X1)/R(X2) = 1 ⟹ α_m + 2α_a = 0** (sympy LOCKED constraint).

φ.1 already forces α_m = 0 ⟹ α_a = 0. Joint protection theorem.

### T2.3: Polarization-Zeeman cross-coupling (NOVEL)
σ.1 LINEAR birefringence Δv_φ = ±gn/(2k) ⟹ photon helicity drift L → polarization rotation θ = gnL/2 ⟹ atomic σ+/σ- selection rule sees ROTATED drive ⟹ differential AC Stark on Zeeman levels m_J=±1.

This is a NOVEL τ.2 + σ.1 + ω.1 cross-prediction (Faraday-like substrate-induced rotation seen via atomic Zeeman).

### T2.6: ω.1 + σ.1 consistency
- ω.1: α_em invariant (axion topological coupling, NOT dilaton-like)
- σ.1: scalar c(X) NULL at leading O(n) (only polarization-dependent split)
- ⟹ τ.2: scalar atomic frequency NULL at leading

Three independent protections from three independent cycles — robust against single-channel violation.

### T2.7: Alt-couplings cross-channel falsified
- m_e ∝ X^α (α≠0): FALSIFIED by Webb/Murphy NULL 1e-7
- ℏ ∝ X^β (β≠0): FALSIFIED by canonical [x,p]=iℏ universality + lab QM
- α_em ∝ X^γ (γ≠0): FALSIFIED by Webb/Murphy NULL 1e-7
- hyperfine ∝ X^δ (δ≠0): FALSIFIED by Hg/Sr 5e-17 NULL

⟹ X-invariant canonical UNIQUE consistent option.

## Phase verdict

**τ.2.Phase 2 PASS (FULL CASCADE 7/7) → Phase 3 forward**

Proper-time formalism sympy-LOCKED. Clock-rate ratio R(X1)/R(X2) = 1 EXACT under φ.1 X→λX gauge. Polarization-Zeeman cross-coupling identified as novel signature linking τ.2 to σ.1 birefringence.
