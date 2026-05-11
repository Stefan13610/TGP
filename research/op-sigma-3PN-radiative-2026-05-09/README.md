---
title: "op-sigma-3PN-radiative-2026-05-09 — Route A escape: σ-coupling at 3PN+ as TT source"
date: 2026-05-09
type: research-cycle
priority: P1_CRITICAL
parent: "[[../op-h-TT-calibration-2026-05-09/Phase_FINAL_close.md]]"
target: "Resolve R5 risk via σ_ij providing radiative TT modes at 3PN+ order"
classification: ESCAPE_ROUTE_CYCLE
status: ACTIVE — Phase 0
predecessor:
  - "[[../op-scalar-mode-LIGO-bound-2026-05-09/]] (DOWNGRADED to STRUCTURAL_CONDITIONAL)"
  - "[[../op-h-TT-calibration-2026-05-09/]] (8/8 PASS adversarial verification)"
  - "[[../op-emergent-metric-from-interaction-2026-05-09/]] (CLOSED, 57/57 PASS)"
related:
  - "[[../op7/OP7_T3_results.md]] (T3.4 ξ_eff = G·Φ_0² LOCK)"
  - "[[../closure_2026-04-26/sigma_ab_pathB/]] (σ_ab Path B audit)"
  - "[[../op-ppE-mapping/Phase1.5_G_SPA_lock.md]] (G_SPA = 48 LOCK)"
tags:
  - escape-route-A
  - 3PN-radiative
  - sigma-coupling
  - second-order-Phi
  - R5-risk-resolution
  - multi-session
---

# op-sigma-3PN-radiative-2026-05-09

## §0 — Mission

Resolve R5 risk z [[../op-scalar-mode-LIGO-bound-2026-05-09/]] (downgraded
2026-05-09 to STRUCTURAL_CONDITIONAL): TGP linearized δg_eff^ij = δ^ij·b_1·δΦ
NIE produkuje h_+, h_× at observer (TT-projection of δ^ij·X = 0 IDENTICALLY).

**Hipoteza Route A (Phase 1 cycle #3 §5.1):** σ_ij at 3PN+ order provides
PROPER radiative TT contribution at 1/r far-field, NIE just 1/r² near-field
z linearized order.

**Mechanism (proposed):**
1. Linear δΦ z binary source (1/r quadrupole, scalar)
2. σ_ij = (∂_iΦ)(∂_jΦ) near-field tensor (1/r² order)
3. σ_ij as **EFFECTIVE TENSOR SOURCE** for 2nd-order δΦ generation
4. 2nd-order δΦ from σ source at far-field: 1/r decay z TENSOR character
5. δg_eff^TT^observer ≠ 0 z 2nd-order δΦ_TT contribution

If mechanism works: TGP framework reproduces GR-like h_+, h_× at 3PN+
(NIE at linearized).

## §1 — Critical context

### §1.1 — Phase 1 cycle calibration finding

[[../op-h-TT-calibration-2026-05-09/Phase2_sympy.py]]: TT-projection of
δ^ij·X = 0 IDENTICALLY (rigorous sympy via projector identity).

⟹ Linear order TGP δg_eff has NO TT modes. Escape route MUST be at
higher order (nonlinear δΦ products giving tensor effective source).

### §1.2 — σ_ij Path B audit

[[../closure_2026-04-26/sigma_ab_pathB/results.md]]: σ_ab is
**operatorowo derived** from ŝ-EOM, NIE independent DOF. M² = 2 m_s².

W RADIATION context, σ_ab ma EOM:
```
□ σ_ab + 2m_s²·σ_ab = source z gradient products
```

This means σ_ab CAN propagate jako wave (M²= 2m_s² > 0 mass term),
ALE at LIGO band: m_s ≪ ω, effectively massless propagation.

If σ_ab propagates jako wave z 1/r far-field (radiative), to provides
h_+, h_× modes at observer.

### §1.3 — PN counting

| Order | Linear δΦ | σ near-field | 2nd-order δΦ (z σ source) |
|---|---|---|---|
| 0 | M_total const (no rad) | 0 (vacuum) | 0 |
| 1 | dipole (vanishing COM) | 0 | 0 |
| 2 (quadrupole) | 1/r leading | 1/r² | ? |
| 3 | 1/r tail | 1/r³ | 1/r ← potential RADIATIVE |
| 4 | 1/r 4PN | 1/r⁴ | 1/r² |

**Key prediction:** σ-induced 2nd-order δΦ becomes RADIATIVE (1/r) at
3PN+ order.

## §2 — Strategy

### §2.1 — Phase plan

| Phase | Goal | Deliverable |
|---|---|---|
| 0 | Setup + scope | This README |
| 1 | Linear δΦ + σ_ij near-field structure (sympy) | Phase1_sympy.py |
| 2 | σ_ij volume integral + tensor moments | Phase2_sympy.py |
| 3 | 2nd-order δΦ from σ source + far-field expansion | Phase3_sympy.py |
| 4 | TT projection of 2nd-order δΦ_2 → h_+, h_× amplitudes | Phase4_sympy.py |
| 5 | Compare z LIGO observed amplitudes (calibration) | Phase5_sympy.py |
| FINAL | Classification: DERIVED / CONDITIONAL / NO_GO | Phase_FINAL_close.md |

### §2.2 — Time budget

**Multi-session estimated: 3-5 sesji.**

**Single-session realistic scope:** Phase 0 + Phase 1 (linear δΦ + σ
near-field structure).

Phase 2-5 deferred do multi-session continuation.

### §2.3 — Falsifier

If 2nd-order δΦ z σ source NIE gives 1/r tensor radiation at 3PN+:
- Cycle classifies STRUCTURAL_NO_GO
- Pivot do Route B (nonlinear δΦ self-coupling) or Route C/D

If 2nd-order δΦ DOES give 1/r tensor radiation:
- Cycle DERIVED
- Magnitude vs LIGO bound (5% scalar polarization) jest separate question

## §3 — Probability assessment

| Outcome | Probability |
|---|---|
| Pełen DERIVED (σ at 3PN+ resolves R5) | 30-40% |
| STRUCTURAL CONDITIONAL (mechanism works but magnitude calibration uncertain) | 30-40% |
| STRUCTURAL_NO_GO (σ doesn't provide radiative TT even at 3PN+) | 20-30% |
| EARLY_HALT (multi-session technical heaviness) | 5-10% |

## §4 — Anti-pattern compliance

- NIE multi-candidate fit: pre-declared mechanism (σ second-order)
- NIE drift hardening: PN counting structural, not fitted
- Honest acknowledgment: Phase 1 cycle #3 verdict was incorrect, this cycle
  attempts to address fundamental issue

## §5 — Cross-references

- [[../op-h-TT-calibration-2026-05-09/Phase_FINAL_close.md]] — Phase 3 cycle #3 error verification
- [[../op-scalar-mode-LIGO-bound-2026-05-09/Phase_FINAL_close.md]] — DOWNGRADED status
- [[../closure_2026-04-26/sigma_ab_pathB/results.md]] — σ_ab Path B audit
- [[../op7/OP7_T3_results.md]] — ξ_eff = G·Φ_0² LOCK
- [[../op-emergent-metric-from-interaction-2026-05-09/Phase4_results.md]] — Phase 4 emergent-metric reference
