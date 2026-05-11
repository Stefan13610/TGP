---
title: "op-h-TT-calibration-2026-05-09 — quantitative h_TT^TGP/h_TT^GR ratio z linearized emergent-metric"
date: 2026-05-09
type: research-cycle
priority: P2_MEDIUM
parent: "[[../op-scalar-mode-LIGO-bound-2026-05-09/Phase_FINAL_close.md]]"
target: "Resolve ~4√π factor mismatch w h_TT^TGP / h_TT^GR amplitude"
classification: QUANTITATIVE_CALIBRATION_CYCLE → ADVERSARIAL_VERIFICATION_CYCLE
status: 🟠 CLOSED — STRUCTURAL_CONDITIONAL_HALT (adversarial finding triggered amendment cascade)
folder_status: closed-conditional-adversarial
sympy_total: "16/16 PASS (Phase 1+2)"
close_date: 2026-05-09
phase_final_close: "[[./Phase_FINAL_close.md]]"
key_finding: "Identified + rigorously verified subtle Phase 3 cycle #3 error: TT-projection of δ^ij·X = 0 IDENTICALLY ⟹ TGP linearized has NO h_+, h_× at observer"
amendment_triggered: "Cycle #3 downgrade STRUCTURAL DERIVED → STRUCTURAL_CONDITIONAL; TGP_FOUNDATIONS §3.6.10.4 amended; PREDICTIONS_REGISTRY amended"
predecessor:
  - "[[../op-scalar-mode-LIGO-bound-2026-05-09/]] (CLOSED, Phase 3 multipole 9/9)"
  - "[[../op-c0-derivation-from-substrate-2026-05-09/]] (c_0 ≈ 4π heuristic)"
  - "[[../op-emergent-metric-from-interaction-2026-05-09/]] (closed STRUCTURAL DERIVED)"
  - "[[../op7/OP7_T3_results.md]] (T3.4 ξ/G ≈ 1.06 GW150914 calibration)"
related:
  - "[[../op-ppE-mapping/Phase1.5_G_SPA_lock.md]] (G_SPA = 48 LOCK, OP-7 framework)"
  - "[[../closure_2026-04-26/sigma_ab_pathB/]] (σ_ab Path B, ξ_eff = G·Φ_0²·O(1))"
tags:
  - h-TT-calibration
  - quantitative-amplitude
  - 4-sqrt-pi-factor
  - GW150914-anchor
---

# op-h-TT-calibration-2026-05-09

## §0 — Mission

Quantitatively resolve **h_TT^TGP / h_TT^GR amplitude ratio** z linearized
emergent-metric framework.

**Phase 3 cycle #3 finding:**
```
h_+^TGP ~ b_1·q/(K_1·Φ_0·c⁴·r) · d²(Q_xx - Q_yy)/dt²
h_+^GR  ~ (2·G/c⁴·r) · d²(Q_xx - Q_yy)/dt²
```

Z M9.1''-canonical b_1 = -4, q² = 4π G K_1 Φ_0² (Phase 5 LOCK):
```
b_1·q/(K_1·Φ_0) = -8√(π·G/K_1)        (dla K_1 = 1: -8√(π G))
Ratio TGP/GR = 8√(π G) / (2G) = 4·√(π/G) ≈ 4√π ≈ 7.09 (dla G=1)
```

**Discrepancy z observation: GW150914 amplitude ~ 1e-21 matches GR.**
TGP predicts ~ 7e-21 — too large.

OR matches GR via O(1) calibration factor (analog OP-7 T3.4 ξ/G ≈ 1.06).

## §1 — Hypothesis: where does 4√π factor come from + how does it resolve?

### §1.1 — Hypothesis A: Phase 5 q² normalization context-specific

Phase 5 derived q² = 4π G K_1 Φ_0² z **STATIC** source Newton coupling
(matching G_eff = G_Newton in 1PN limit).

W RADIATION context (1.5PN+), q effectively might be DIFFERENT — analog
do how σ-coupling enters at 2PN-orbital not 1PN.

**Possible resolution:** q_radiation = q_static / (4√π) ≈ q_static / 7.09
gives h_TT^TGP / h_TT^GR ≈ 1.

### §1.2 — Hypothesis B: b_1 normalization ambiguity

b_1 = f'(1) where f = 1/A. Z M9.1'' A_M911 = ψ/(4-3ψ):
- A(1) = 1, A'(1) = 4 ⟹ a_1 = 4
- 1/A Taylor: B_inv has b_1 = -a_1 = -4 (PPN constraint)

Sprawdzić: czy "b_1" w h_TT formula odnosi się do A czy 1/A Taylor coefs?
Phase 3 used b_1 = f'(1) = -4. Może powinno być b_1 = A'(1) = +4?

**Possible resolution:** factor of 4 reduction if proper normalization.

### §1.3 — Hypothesis C: GW150914 ξ/G calibration jest already correct

OP-7 T3.4 derived ξ = G·Φ_0² z GW150914 numerical match z ξ/G ≈ 1.06.
This means the EFFECTIVE coupling is properly normalized, but Phase 5
LOCK (q² = 4π G K_1 Φ_0²) used DIFFERENT normalization derived from
Newton limit.

**Possible resolution:** unify normalization conventions z OP-7 T3.4.
Effective q_radiation might be smaller than naive Phase 5 q_static
by factor sqrt(ξ_eff/G) ≈ 1.03 (from GW150914 calibration).

But this is much smaller than 7×. So Hypothesis C alone doesn't fix.

### §1.4 — Hypothesis D: missing factor of 1/(2√π) z proper TT projection

In linearized GR derivation, h_TT involves spatial TT-projection operator
P_TT^ij,kl which has trace structure z normalization 1/2 (standard).

W TGP linearized, similar TT projection but z SCALAR field source δΦ Y_2m,
projection might include additional factor 1/(2√π) z spherical harmonic
normalization.

**Possible resolution:** factor 1/(2√π) reduction → ratio 4√π/(2√π)·√π/π... 
hmm not clean. Worth checking.

## §2 — Phase plan

### Phase 0: Setup (this README)

### Phase 1: Re-derive h_TT^TGP carefully
- Explicit linearized GW from binary in TGP framework
- Compare term-by-term z standard GR quadrupole derivation
- Identify exact factor source

### Phase 2: Test hypotheses A-D
- Sympy verification each hypothesis
- Numerical comparison

### Phase 3: Cross-check z OP-7 T3.4 GW150914 calibration
- Use ξ/G ≈ 1.06 jako anchor
- Determine proper normalization

### Phase FINAL: Verdict
- DERIVED: ratio resolved cleanly, consistent z GR
- STRUCTURAL_CONDITIONAL: identifies framework adjustment needed
- STRUCTURAL_NO_GO: 4√π factor genuine, framework prediction differs from GR

## §3 — Probability assessment

| Outcome | Probability |
|---|---|
| Pełen DERIVED (ratio = 1 + small calibration) | 35-50% |
| STRUCTURAL CONDITIONAL (factor needs adjustment, plausible) | 30-40% |
| STRUCTURAL_NO_GO (TGP genuinely predicts factor 7× GR) | 15-25% |

## §4 — Time budget

**Single-session realistic scope:** Phase 0-1 (re-derivation + hypothesis A test).
**Multi-session:** 2-3 sesji dla full Phase 1-3 + FINAL close.

## §5 — Cross-references

- [[../op-scalar-mode-LIGO-bound-2026-05-09/Phase3_results.md]] — 4√π identification
- [[../op-c0-derivation-from-substrate-2026-05-09/Phase_FINAL_close.md]] — c_0 = 4π heuristic
- [[../op7/OP7_T3_results.md]] — ξ = G·Φ_0² LOCK + GW150914 calibration
- [[../op-emergent-metric-from-interaction-2026-05-09/Phase5_sympy.py]] — Phase 5 q² LOCK
- [[../op-ppE-mapping/Phase1.5_G_SPA_lock.md]] — SPA chain LOCK
- [[../../TGP_FOUNDATIONS.md]] §3.6.10 — joint cycles status
