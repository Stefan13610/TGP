---
title: "op-T34-normalization-amendment-2026-05-09 — Audit + amendment of T3.4 ξ_eff normalization chain"
date: 2026-05-09
type: research-cycle
priority: P1_CRITICAL
parent: "[[../op-sigma-3PN-radiative-2026-05-09/Phase2_adversarial_verification.md]]"
target: "Resolve ξ_eff factor-4 normalization gap; restore framework GR amplitude consistency"
classification: AMENDMENT_CYCLE
status: ACTIVE — Phase 1
predecessor:
  - "[[../op-sigma-3PN-radiative-2026-05-09/Phase2_results.md]] (Phase 2 finding)"
  - "[[../op-sigma-3PN-radiative-2026-05-09/Phase2_adversarial_verification.md]] (adversarial confirmation)"
related:
  - "[[../op7/OP7_T3_results.md]] (T3.4 source z gaps)"
  - "[[../op7/op7_t3_4_xi_coupling.py]] (Gap 1 line 132, Gap 2 line 140)"
  - "[[../op-c0-derivation-from-substrate-2026-05-09/Phase1_sympy.py]] (line 65: ξ = 4π·G·Φ_0² — third inconsistency)"
tags:
  - amendment-cycle
  - T34-normalization
  - factor-4-gap
  - xi-eff-correction
  - chain-audit
---

# op-T34-normalization-amendment-2026-05-09

## §0 — Mission

**Resolve T3.4 ξ_eff normalization gap z explicit first-principles re-derivation.**

Phase 2 of `op-sigma-3PN-radiative` + adversarial verification ujawniły że:
- **Three different ξ_eff values** appear w TGP cycle chain
- All three CONFLICT z each other
- Resolution requires clean re-derivation z full factor tracking

## §1 — Three inconsistent ξ_eff values

| Source | Stated ξ_eff | Implication for h^σ/h^GR (z c_0 = 4π LOCK) |
|---|---|---|
| **OP-7 T3.4 text** ([[../op7/OP7_T3_results.md]] §2 T3.4) | ξ = G·Φ_0² | ratio = 1/4 = 0.25 (TGP underpredicts) |
| **c_0 cycle Phase 1 sympy** ([[../op-c0-derivation-from-substrate-2026-05-09/Phase1_sympy.py]] line 65) | ξ = 4π·G·Φ_0² | ratio = π ≈ 3.14 (TGP overpredicts) |
| **Adversarial corrected** (per Phase 2 audit) | ξ = 4·G·Φ_0² | ratio = 1.0 EXACT (GR match) |

**Wniosek:** Każda z trzech wartości daje fundamentalnie inną fizyczną predykcję.
Tylko jedna może być poprawna. Cycle musi rozstrzygnąć **niezależnie**.

## §2 — Strategy: clean re-derivation z first principles

### §2.1 — Pre-declared methodology

NO inheritance z prior cycles. Build fresh:

1. **Standard textbook references** (NIE TGP cycles):
   - Misner-Thorne-Wheeler "Gravitation" (1973) §36.3-36.10
   - Maggiore "Gravitational Waves Vol I" (2008) §3.1-3.4
   - Wald "General Relativity" (1984) §11.2

2. **Path A Lagrangian as STARTING POINT** (OP-7 T3.1 form):
   ```
   L_σ = -(1/4)(∂_μσ_ab)(∂^μσ^ab) - (1/2)m_σ²σ_ab σ^ab - (ξ_eff/2)σ_ab T^{ab,TT}
   ```

3. **Variational principle z explicit factor tracking** dla EOM derivation
4. **Standard retarded Green function** dla massless wave equation
5. **Standard PN identity** ∫T^ij = (1/2)·d²Q/dt² (Maggiore Eq. 3.81)
6. **TT-projection** at observer z standard projector
7. **Matching** to GR h^TT_GR = (2G/c⁴r)·d²Q^M_TT/dt² (MTW Eq. 36.22)

### §2.2 — Required ξ_eff value derivation

Find ξ_eff such that h_TGP^σ = h_GR EXACTLY. Document each factor explicitly.

### §2.3 — Verify cycle LOCKS preserved

After derivation, check:
- c_0 = 4π z cycle #1 — does this LOCK survive?
- κ_σ = 1/(3π) z cycle #2 — independent of T3.4 amendment
- c_0·κ_σ = 4/3 z joint — should remain valid (phase-only)
- β_ppE = 0 — should remain valid

### §2.4 — Cross-check Phase 2 finding

Confirm post-amendment ratio h^σ/h^GR z corrected ξ_eff.

## §3 — Phase plan

| Phase | Goal | Deliverable |
|---|---|---|
| 0 | Setup + scope | This README |
| 1 | First-principles re-derivation z factor tracking | Phase1_sympy.py + Phase1_results.md |
| 2 | Cross-check c_0 LOCK + joint LOCK preservation | Phase2_sympy.py |
| FINAL | Amendment notice + classification | Phase_FINAL_close.md |

**Estimated effort:** 1-2 sesji.

## §4 — Falsifier

If first-principles re-derivation gives ξ_eff inconsistent z c_0 = 4π LOCK:
- Cycle classifies STRUCTURAL_NO_GO
- Both T3.4 + cycle #1 require deeper revision
- Pivot do nonlinear δΦ self-coupling (Route B) lub framework extension

If first-principles re-derivation gives ξ_eff = 4·G·Φ_0² consistent z c_0 = 4π:
- T3.4 ξ_eff = G·Φ_0² CONFIRMED to be off by factor 4
- Amendment notice issued
- Phase 2 σ-3PN cycle upgraded to STRUCTURAL_DERIVED post-amendment

## §5 — Anti-pattern compliance

- NIE inherit factors z prior cycles (clean derivation)
- NIE multi-candidate fit (specific factor predicted from first principles)
- NIE drift hardening (text books only as reference)
- Honest reporting: report whatever ξ_eff first-principles gives, NOT predetermined value

## §6 — Cross-references

- [[../op-sigma-3PN-radiative-2026-05-09/Phase2_adversarial_verification.md]] — predecessor finding
- [[../op7/op7_t3_4_xi_coupling.py]] — T3.4 source z gaps
- [[../op7/OP7_T3_results.md]] — T3.4 results (NEEDS amendment notice)
- [[../op-c0-derivation-from-substrate-2026-05-09/Phase1_sympy.py]] line 65 — c_0 cycle ξ form
- [[../op-c0-derivation-from-substrate-2026-05-09/Phase_FINAL_close.md]] — c_0 LOCK
- [[../op-kappa-sigma-2body-PN-2026-05-09/Phase_FINAL_close.md]] — κ_σ LOCK (independent)
- [[../op-h-TT-calibration-2026-05-09/Phase_FINAL_close.md]] — predecessor adversarial protocol
