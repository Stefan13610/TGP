---
title: "UV.1 — AS NGFP first-principles closure for N_A normalization (3-phase mini-cycle)"
date: 2026-04-29
cycle: UV.1
status: ACTIVE (Phase 1 PRE-EXECUTION)
parent: "[[../../INDEX.md]]"
siblings:
  - "[[../op-xi-photon-ring/program.md]]"
  - "[[../op-cross-sector-charge/program.md]]"
  - "[[../op-bh-alpha-threshold/program.md]]"
  - "[[../op-sc-alpha-origin/program.md]]"
tags:
  - TGP
  - uv-completion
  - asymptotic-safety
  - NGFP
  - Reuter
  - heat-kernel
  - a2-coefficient
  - photon-ring
  - falsification
---

# UV.1 — AS NGFP first-principles closure for N_A normalization

> **Cel programu:** Domknąć **single residual N_A = 8.7719 algebraic
> provenance** zostawiony przez ξ.1.Phase3 z asymptotic-safety NGFP
> first-principles (Reuter 1998 / Eichhorn 2018). ξ.1.Phase3 wykazał,
> że AS NGFP jest **najbliższym UV-route match dla N_A** (Δ 0.068%), ale
> illustrative — nie derived. UV.1 attempts: derive N_A directly z
> Litim invariant g*·λ* = 0.1349 + η_N* = -2 + heat-kernel a₂ scaling
> pod NGFP RG flow. Promote UV7 → DERIVED, ξ.1 → DERIVED (best case),
> XS.1 → DERIVED (cascade).

---

## Motywacja

ξ.1.Phase3 zostawiło **dokładnie jedną zewnętrzną zaślepkę:**
- N_A = 500/57 = 8.7719 (TGP a₂-derived from Frame A target_shift = 0.114)
- najbliższy natural integer: 9 (Δ 2.6%)
- najbliższy UV-route match (illustrative): AS NGFP ≈ 8.7660 (Δ 0.068%)

ξ.1 nie wyprowadziło N_A **algebraicznie** — fit AS NGFP był heurystyczny
(`9.0 * (1 − 0.026)`), nie derived z g*, λ*, η_N*. UV.1 attempts to close
this gap by:

1. **Phase 1 audit:** wszystkie 5 fundamental AS NGFP constants LOCKED w
   istniejących closurach (Phase 3.A KEYSTONE).
2. **Phase 2 derivation:** wyprowadzić N_A = 8.7719 lub blisko (Δ < 0.1%)
   z {g*, λ*, η_N*, Litim invariant 0.1349}.
3. **Phase 3 predictions:** UV1–UV3 falsifiable predictions tied to AS RG
   flow (LIGO O5 dispersion, ngEHT photon-ring sharpening, RG cross-scale).

---

## 3-fazowy plan

### Phase 1 — AS NGFP foundational audit (5 sub-tests)

**Cel:** Sprawdzić, czy fundamental inputs do N_A derivation są LOCKED
w istniejących Phase 3.A / Phase 3.E closurach.

| # | Test | Falsification rule |
|---|------|--------------------|
| UV1.1 | Litim invariant g*·λ* = 0.1349 LOCKED (Phase 3.A) | drift > 5% vs Reuter 1998 → NGFP foundation broken |
| UV1.2 | η_N* = −2 LOCKED (NGFP anomalous dimension; Reuter 1998) | η_N* drift > 10% → NGFP fixed-point ambiguous |
| UV1.3 | Scale separation m_Φ/Λ_EFT = 60.93 dex > 50 dex gate (Phase 3.A) | dex < 50 → EFT-NGFP bridge broken |
| UV1.4 | T-FP IR consistency 12/12 POSITIVE (Phase 3.A.5) | <12/12 → IR-UV bridge inconsistent |
| UV1.5 | Heat-kernel a₂ → α₀ reproducibility 0.0% drift (Phase 3.E.4) | drift > 0.1% → a₂ frame inconsistent |

**Verdict gate:** 5/5 PASS → all foundational AS constants LOCKED, proceed
Phase 2 derivation. ≤ 4/5 → audit gap pierwsze, Phase 2 deferred.

### Phase 2 — N_A first-principles derivation z NGFP scaling (7 sub-tests)

**Cel:** Wyprowadzić N_A = 8.7719 (lub blisko, Δ < 0.5%) z AS NGFP
parameters (g*, λ*, η_N*, Litim invariant).

| # | Test | Falsification rule |
|---|------|--------------------|
| UV2.1 | a₂ vacuum form on FRW under NGFP (cross-link ξ.1.Phase2) | a₂ structure drifts od Birrell-Davies pod NGFP |
| UV2.2 | V''(Φ_eq=1) = 2β z β UV-fixed by NGFP RG flow | V''(1) free parameter pod NGFP |
| UV2.3 | NGFP scaling of heat-kernel renormalization (heat-kernel under η_N* = −2) | a₂ does NOT pick up (1 + η_N*/2)·corrections |
| UV2.4 | N_A natural-integer attempt: N_A = 9·(1 − Δ) gdzie Δ z NGFP (Δ < 5%) | Δ z NGFP > 5% → not natural integer base |
| UV2.5 | N_A algebraic attempt: N_A = f(g*, λ*) closed-form match Δ < 0.5% | best fit Δ > 0.5% → N_A not first-principles z NGFP |
| UV2.6 | 2-loop residual interpretation: 0.068% within 2-loop EFT band? | residual > 1% → 2-loop hypothesis falsified |
| UV2.7 | Classification: DERIVED (sympy-clean) / PARTIALLY (Δ < 0.5% close) / STRUCTURAL (Δ > 0.5%) | chain inconsistent → no classification |

**Verdict gate:** ≥ 6/7 PASS → UV.1 DERIVED or PARTIALLY DERIVED, proceed
Phase 3. ≤ 5/7 → STRUCTURAL HINT, Phase 3 still proceeds (limited registry).

### Phase 3 — UV-cascade predictions + status promotions (6 sub-tests)

**Cel:** 5–6 new predictions UV1–UV6 (or UV-A through UV-F), cascade
status promotions through ξ.1, XS.1, UV7.

| # | Test | New prediction |
|---|------|----------------|
| UV3.1 | N_A LOCKED at NGFP-derived value (post-Phase 2) | UV1: N_A = derived value (Δ from AS NGFP closed-form) |
| UV3.2 | LIGO O5 2027+ GW dispersion bound on NGFP RG scaling | UV2: graviton dispersion δ < 10⁻¹⁵ from NGFP η_N* |
| UV3.3 | ngEHT 2030+ photon-ring sharpening post-UV.1: N_A precision band | UV3: ngEHT × NGFP combined ≤0.3% by 2030+ |
| UV3.4 | ξ-factor RG-invariance verified across UV scales (k_IR → k_UV) | UV4: ξ.1.Phase3 RG-invariance confirmed at AS endpoint |
| UV3.5 | Status cascade: ξ.1 promotion, XS.1 promotion, UV7 promotion | UV5: ξ.1 → DERIVED (or PARTIALLY refined²), XS.1 → DERIVED |
| UV3.6 | Long-term track handoff: UV1–UV7 reduced to UV2–UV7 (UV1 closed) | UV6: research-track items 7 → 6 |

**Verdict gate:** ≥ 5/6 PASS → UV.1 program END, 5–6 new predictions, ξ.1
i XS.1 promoted (best case) or refined (fallback).

---

## Cumulative target

18 sub-tests across 3 phases (5 + 7 + 6), nieco mniej niż ξ.1 (5+7+7=19),
SC.1 (4+6+7=17), BH.1 (5+7+7=19), XS.1 (5+7+7=19). UV.1 jest "tightest"
pod kątem expected outcomes — Phase 2 może osiągnąć tylko PARTIALLY
DERIVED (rezerwujemy realistyczność co do AS NGFP first-principles power).

Master ledger target: 355 → **373** (+18 z UV.1) at program END.

---

## Predecessors

- [`../op-xi-photon-ring/Phase3_results.md`](../op-xi-photon-ring/Phase3_results.md) — ξ.1 program END (N_A = 8.7719 left open)
- [`../op-xi-photon-ring/Phase2_results.md`](../op-xi-photon-ring/Phase2_results.md) — heat-kernel a₂ derivation
- [`../op-phase3-uv-completion/Phase3_A_results.md`](../op-phase3-uv-completion/Phase3_A_results.md) — AS NGFP KEYSTONE 12/12
- [`../op-phase3-uv-completion/Phase3_E_results.md`](../op-phase3-uv-completion/Phase3_E_results.md) — UV7 STRUCTURAL-DERIVED (Δ_target a₂ frame)
- [`../op-phase3-uv-completion/Phase3_R_final_results.md`](../op-phase3-uv-completion/Phase3_R_final_results.md) — 4-of-4 UV synthesis
- [`../op-cross-sector-charge/Phase3_results.md`](../op-cross-sector-charge/Phase3_results.md) — XS.1 PARTIALLY DERIVED (refined)

## Decision gate after each phase

| Phase | PASS gate | Action on PASS | Action on FAIL |
|-------|-----------|-----------------|------------------|
| Phase 1 | 5/5 | proceed Phase 2 derivation | abandon (NGFP foundation not LOCKED) |
| Phase 2 | ≥6/7 | proceed Phase 3 (DERIVED or PARTIALLY) | Phase 3 with STRUCTURAL HINT |
| Phase 3 | ≥5/6 | UV.1 program END, 5–6 new predictions | partial registration |

## Open issues (long-term, NOT blocking UV.1)

After UV.1 closure, remaining long-term track:
- Full UV-complete renormalizability proof (UV completeness theorem)
- String / LQG / CDT alternative routes (only AS prioritized w UV.1)
- Phase 0 minimal axiomatic kernel for NGFP-substrate bridge

Te pozostają w `op-uv-renormalizability-research/` long-term track
(UV2–UV7 supplement, redukcja z 7 do 6 items po UV.1 closure).

---

## Cross-references

- AS NGFP framework: Reuter 1998 PRD 57 (FRG fixed point); Eichhorn 2018 reviews
- Heat-kernel: Birrell-Davies 1982 Eq. 6.45; Avramidi 2000
- TGP closures: closure_2026-04-26 (F4 rational), Phase 2.B.3 (α₀ chain)
- ξ.1 N_A target: 500/57 = 8.7719 (Frame A 1-loop a₂-corrected)
