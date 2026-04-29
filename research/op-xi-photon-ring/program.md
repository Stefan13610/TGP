---
title: "ξ.1 — ξ-factor photon-ring resolution (3-phase mini-cycle)"
date: 2026-04-29
cycle: ξ.1
status: ACTIVE (Phase 1 PRE-EXECUTION)
parent: "[[../../INDEX.md]]"
siblings:
  - "[[../op-cross-sector-charge/program.md]]"
  - "[[../op-bh-alpha-threshold/program.md]]"
  - "[[../op-sc-alpha-origin/program.md]]"
tags:
  - TGP
  - xi-factor
  - photon-ring
  - heat-kernel
  - a2-coefficient
  - falsification
---

# ξ.1 — ξ-factor photon-ring resolution

> **Cel programu:** Zamknąć **single unresolved ξ factor** w XS.1
> PARTIALLY DERIVED — różnicę między **target_shift = 0.114** (F4 sympy
> rational LOCKED) a **target_shift = (1/2)(1 − 3/3.88) = 0.1134** (Phase
> 2 strict form). Promote XS.1 z PARTIALLY DERIVED → **DERIVED** (best
> case) lub map ξ-factor na UV-route selection (fallback).

---

## Motywacja

XS.1 zamknęło substrate-action derivation √α₀ = κ_TGP z **6 LOCKED
inputs** + ZERO free parameters, ale pozostawiło **jeden O(1) ξ factor**
unresolved:

```
F4 rational  : 1069833/264500 ≈ 4.04472   (|Δ|/κ_TGP² = 0.084%)
Phase 2 strict: (1/2)(1−3/3.88)/0.168² = 4.0179 (|Δ|/κ_TGP² = 0.747%)
Split        : 1.05%
```

Hipoteza ξ.1: split jest albo (A) computational artifact (F4 = a₂-derived
true value, strict = bare approximation), albo (B) UV-completion-dependent
(różne UV completion routes give różne ξ frames), albo (C) calibration
convention (both consistent within EFT 1% precision).

UV7 (PREDICTIONS_REGISTRY): currently `STRUCTURAL-POSTULATE`:
> Δ_target = 0.114 in heat-kernel a₂ frame; a₂ ⊃ (1/2)V''²; ξ_geom=1, α(α−1)=2

ξ.1 audit czy a₂ first-principles derivation z TGP inputs (F1+F2+F3+F4)
selects unique frame (A vs B) lub deklaruje degeneracy (C).

---

## 3-fazowy plan

### Phase 1 — EFT photon-ring frame audit (5 sub-tests)

**Cel:** Sprawdzić, czy fundamental inputy a₂ derivation są LOCKED w
istniejących closurach (ξ_geom=1, α(α−1)=2, ψ_ph=1.168, F4 rational
provenance, target_shift strict provenance).

| # | Test | Falsification rule |
|---|------|--------------------|
| ξ1.1 | ξ_geom = 1 derivation z Phase 2.B.2 (M9.1″ vacuum) | jeżeli ξ_geom ma free parameter > 1 wariancji → frame nie jest closed |
| ξ1.2 | α(α−1) = 2 derivation z Phase 1.A.1 (Theorem alpha2) | jeżeli α≠2 dopuszczone → a₂ structure ambiguous |
| ξ1.3 | ψ_ph − 1 = 0.168 z M9.2-D photon-ring | jeżeli ψ_ph drift > 0.5% → photon-ring frame nie zamknięty |
| ξ1.4 | F4 rational 1069833/264500 provenance audit | jeżeli F4 nie z arithmetic identity Δ_target/[(ψ_ph−1)²·ξ_geom] → trace inconsistency |
| ξ1.5 | Phase 2 strict (1/2)(1 − 3/3.88) / 0.168² provenance audit | jeżeli strict nie z direct geometric photon-ring shift → strict not bare-form |

**Verdict gate:** 5/5 PASS → all inputs LOCKED, proceed Phase 2 derivation.
≤ 4/5 → audit gap pierwsze, Phase 2 deferred.

### Phase 2 — Heat-kernel a₂ first-principles derivation (7 sub-tests)

**Cel:** Wyprowadzić Δ_target z a₂ heat-kernel coefficient pod TGP
substrate inputs (F1 single-Φ, F2 K_geo·φ⁴, F3 φ_eq=1, F4 α₀ anchor).
Decide czy F4 rational lub Phase 2 strict (lub trzecia forma) jest
true a₂-derived value.

| # | Test | Falsification rule |
|---|------|--------------------|
| ξ2.1 | Heat-kernel a₂ structure: a₂ ⊃ (1/2)V''(Φ)² + Ricci terms | a₂ structure deviates z Birrell-Davies/Avramidi |
| ξ2.2 | TGP V''(Φ) evaluation pod F2 (K_geo·φ⁴): V'' = … | V'' depends on free parameter not LOCKED |
| ξ2.3 | M9.1″ FRW background: Ricci R suppressed 10⁻¹²² dex | Ricci contribution NOT subdominant |
| ξ2.4 | a₂ → Δ_target conversion: Δ_target = ξ_geom · α(α−1) · (a₂-ratio) | conversion has free O(1) parameter |
| ξ2.5 | Frame A test: a₂ derivation yields F4 sympy rational 1069833/264500 | jeżeli a₂ NOT 4.0447 sympy → Frame A FALSIFIED |
| ξ2.6 | Frame B test: a₂ derivation yields Phase 2 strict 0.1134 | jeżeli a₂ NOT 4.0179 → Frame B FALSIFIED |
| ξ2.7 | Classification: DERIVED (Frame A or B unique) / PARTIALLY DERIVED (both consistent within UV uncertainty) / STRUCTURAL HINT (both possible, no a₂ derivation) | symbolic + numeric chain produces no decision |

**Verdict gate:** ≥ 6/7 PASS → ξ.1 DERIVED or PARTIALLY DERIVED, proceed Phase 3.
≤ 5/7 → ξ.1 STRUCTURAL HINT, Phase 3 still proceeds (limited registry update).

### Phase 3 — Cross-sector falsification + UV-route map (7 sub-tests)

**Cel:** Wygenerować predictions XI1–XI3 (Frame A) lub XS7–XS9 (Frame B
UV-route), zaktualizować precision gates dla XS1–XS6, sharpen
ngEHT/LnH₉ falsification triggers.

| # | Test | New prediction |
|---|------|----------------|
| ξ3.1 | If Frame A: F4 LOCKED-DERIVATIVE; tighten ngEHT precision gate 5% → 0.5% | XI1: combined ngEHT × SC ≤0.5% achievable post-ξ.1 |
| ξ3.2 | If Frame B: register UV-route selection pointer (which AS/string/LQG/CDT prefers Frame B) | XS7: AS NGFP / LQG / CDT pojedyncza UV preferowana przez Frame B |
| ξ3.3 | Cross-sector consistency cascade: F4 + F5 + F6 status update | XI2: F-cluster locked-derivative under ξ.1 |
| ξ3.4 | Photon-ring lever recalculation: minimum ngEHT mass-sample size to reject identity at 5σ | XI3: ngEHT 2030+ minimum 2-3 sources for ξ-tightened identity |
| ξ3.5 | Lepton orthogonality re-audit (XS3): does ξ-derivation reveal latent κ_TGP factor in Koide? | XS3 strengthen / weaken |
| ξ3.6 | RG stability of ξ across IR → UV (one-loop) | XI4: ξ-factor RG-invariant lub flow direction noted |
| ξ3.7 | PREDICTIONS_REGISTRY entries XI1–XI3 (Frame A) or XS7–XS9 (Frame B) + XS1 precision update | 3-6 new pre-registered predictions |

**Verdict gate:** 7/7 PASS → ξ.1 program END, 3-6 new predictions, XS.1 promoted (best case) or UV-route registered (fallback).

---

## Cumulative target

19 sub-tests across 3 phases (5 + 7 + 7), parallel structure to XS.1
(19 = 5+7+7), BH.1 (19 = 5+7+7), SC.1 (17 = 4+6+7).

Master ledger target: 336 → **355** (+19 z ξ.1) at program END.

---

## Predecessors

- [`../op-cross-sector-charge/Phase2_results.md`](../op-cross-sector-charge/Phase2_results.md) — XS.1 PARTIALLY DERIVED, ξ unresolved
- [`../op-cross-sector-charge/Phase3_results.md`](../op-cross-sector-charge/Phase3_results.md) — XS5 LOCKED-derivative (F4 sub-percent)
- [`../op-bh-alpha-threshold/Phase2_results.md`](../op-bh-alpha-threshold/Phase2_results.md) — Phase 2 strict form (1/2)(1−3/3.88)
- [`../op-phase3-uv-completion/Phase3_E_results.md`](../op-phase3-uv-completion/Phase3_E_results.md) — UV7 STRUCTURAL-POSTULATE (Δ_target a₂ frame)
- [`../op-phase2-quantum-gravity/Phase2_B_results.md`](../op-phase2-quantum-gravity/Phase2_B_results.md) — F4 rational provenance (Phase 2.B.3)
- [`../closure_2026-04-26/`](../closure_2026-04-26/) — F4 sympy rational LOCKED

## Decision gate after each phase

| Phase | PASS gate | Action on PASS | Action on FAIL |
|-------|-----------|-----------------|------------------|
| Phase 1 | 5/5 | proceed Phase 2 derivation | abandon (premise inputs not LOCKED) |
| Phase 2 | ≥6/7 | proceed Phase 3 (DERIVED or PARTIALLY) | Phase 3 with STRUCTURAL HINT classification |
| Phase 3 | ≥6/7 | ξ.1 program END, 3-6 new predictions | partial registration |

## Open issues (long-term, NOT blocking ξ.1)

After ξ.1 closure, remaining long-term track:
- Full 2-loop a₂ derivation z renormalizable substrate-action
- UV completion route selection (which AS/string/LQG/CDT actually picks ξ-frame)
- Phase 0 minimal axiomatic kernel for photon-ring shift

Te pozostają w `op-uv-renormalizability-research/` long-term track.
