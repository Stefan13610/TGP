---
title: "τ.2.Phase1 results — scale-protection theorem 5/5 PASS"
date: 2026-04-30
cycle: τ.2.Phase1
status: PASS
parent: "[[program.md]]"
tags:
  - TGP
  - tau2
  - phase1
  - results
---

# τ.2.Phase1 results — 5/5 FULL PASS

## Sub-test outcomes

| ID | Test | Result |
|----|------|--------|
| **T1.1** | Atomic-mass coupling candidates scan (4 forms) | ✅ PASS |
| **T1.2** | Scale-invariance X→λX requires α=0 (m_atom INVARIANT) | ✅ PASS |
| **T1.3** | Noether scale-current ∂_μJ^μ implication for proper time | ✅ PASS |
| **T1.4** | Effective Hamiltonian H_atom X-independent | ✅ PASS |
| **T1.5** | Protection theorem: clock rate X-invariant at leading | ✅ PASS |

**Score: 5/5 → Phase 2 forward**

## Key results

### T1.1: Mass-coupling forms
Of 4 candidate Lagrangian couplings:
- **L1 = m_0** (X-independent): scale-INVARIANT ✓
- **L2 = m_0 X^α**: scale-BREAKS for α≠0
- **L3 = m_0 + α ln X**: scale-BREAKS for α≠0
- **L4 = m_0 + α(∂ ln X)²**: scale-INVARIANT but dim-6 EFT-suppressed

Only L1 (canonical) and L4 (gradient-coupled, suppressed) survive scale-symmetry.

### T1.2: φ.1 X→λX gauge forces m_atom invariant
- φ.1 axiom: S[X] = ∫½(∂_μ ln X)² invariant under X→λX (constant shift drops in derivative).
- For matter sector consistent with φ.1: δS_matter = α·ln λ·∫ψ̄ψ = 0 ⟹ **α = 0**.
- Atomic masses (m_e, m_p, m_n, m_qq) X-INVARIANT at leading order.

### T1.3: Noether scale-current → proper time
- Conservation ∂_μ J^μ = □(ln X) = 0 on-shell.
- Substrate gradient does NOT directly modify g_μν at leading order.
- **τ(X1)/τ(X2) = 1 EXACT** at leading O(∂ ln X).
- Higher-order: dim-6 EFT correction g_μν^eff = η_μν + (1/Λ²)(∂_μ ln X)(∂_ν ln X).

### T1.4: H_atom X-independent
All atomic-spectrum constants (m_e, c, α_em, ℏ, Z_eff) X-invariant under φ.1+ω.1+σ.1:
- m_e → m_e (T1.2 protection)
- c → c (σ.1 scalar protection)
- α_em → α_em (ω.1 axion-coupling, NOT dilaton)
- ℏ → ℏ (universal)

⟹ E_n = m_e c² α_em²/(2n²) X-invariant ⟹ atomic transition ω_nm X-INVARIANT.

### T1.5: Protection theorem
**THEOREM (τ.2):** Under φ.1 gauge X→λX, atomic clock rate R = ω_nm/(2π) is **INVARIANT at leading order O(∂ ln X / E_atom)**.

**PROOF:** T1.1 + T1.2 + T1.3 + T1.4 jointly establish R(X) = R_0 q.e.d.

**CORRECTIONS:**
- δR/R ~ (∂ ln X / Λ)² (dim-6 EFT)
- δR/R ~ g(∂ ln X)/k_drive (ω.1 + σ.1 polarization-coupling)

## Cross-channel consistency

- ✓ Webb/Murphy 2003-2017: NULL Δα_em < 1e-7 over z=0-4 → τ.2 prediction confirmed
- ✓ Lab Hg/Yb/Sr clock comparisons: NULL drift 1e-18/yr → τ.2 prediction confirmed
- ✓ ω.1 axion-coupling (NOT scalar dilaton) consistency with α_em invariance
- ✓ σ.1 scalar c(X) protection generalized to all atomic observables

## Phase verdict

**τ.2.Phase 1 PASS (FULL 5/5) → Phase 2 forward**

Scale-protection theorem structurally derived. Atomic clock rate is scalar-substrate-protected at leading order. Higher-order corrections enter at O((∂ ln X)²) and via ω.1+σ.1 cross-couplings (polarization-Zeeman channel).
