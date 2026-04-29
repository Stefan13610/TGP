---
title: "ε.1 program — ε_ph structural decomposition (photon-ring scale)"
date: 2026-04-29
type: program
status: ACTIVE
predecessor: "[[../op-uv-as-ngfp/Phase3_results.md]]"
tags:
  - TGP
  - epsilon-photon-ring
  - structural-decomposition
  - heat-kernel
  - F4-anchor
  - cross-sector
---

# ε.1 program — ε_ph structural decomposition

> **Cel:** Zamknąć ε_ph (M9.1″ refined photon-ring scale = 0.16788) jako
> **structural decomposition ε_ph = ψ_ph − 1** z M9.1″ null geodesic
> ψ_ph = 4/(3 + 0.4250) = 1.16788, oraz **falsyfikować** alternatywne
> identity candidates (1.168/(2π), normalized ratios, κ_TGP-based).
> Status target: PARTIALLY DERIVED → DERIVED via cross-sector consistency.

---

## Strategiczny kontekst

ε_ph appears w trzech kluczowych miejscach:

1. **F4 chain α₀ = target_shift_F4 / ε_ph²** = 0.114 / 0.16788² = 4.04489
   (drift 0.0038% z F4 sympy 1069833/264500) — heat-kernel reproducibility
   (UV.1.Phase1.UV1.5)
2. **ξ.1 frame audit** (ξ_geom · α(α−1) · ε_ph² / 2 / N_ref) — geometric
   chain dla target_shift derivation
3. **M9.1″ photon-ring scale** — refined null-geodesic anchor (gradient-from
   coarse 0.168 M9.2-D do precise 0.16788 closes 36× drift improvement)

**Open question:** Is ε_ph = 0.16788 an **algebraic identity** (closed-form
ratio z fundamental constants) lub **structural footprint** ψ_ph − 1
(M9.1″ geodesic-derived)?

**Pre-cycle audit (UV.1 era):**
- Candidate "1.168/(2π) = 0.18589" → drift +10.7% **REJECTED**
- Candidate "(ψ_ph−1)/ψ_ph = 0.14376" → drift −14.3% **REJECTED**
- Candidate "1/(2π·κ_TGP) ~ 1.99" → dimensional fail **REJECTED**
- Candidate "ε_ph = ψ_ph − 1" → tautological (decomposition, NOT identity)
- Candidate "ε_ph = 1/(2β² + δ_F4)" → TBD audit

→ **Hypothesis:** ε_ph is **structural decomposition** ψ_ph − 1, where ψ_ph
itself emerges from M9.1″ null geodesic computation. Identity candidates
all fail at >10% drift level. F4 chain provides **algebraic anchor**
(target_shift = 0.114 = 57/500), so α₀ closure via target_shift/ε_ph²
fixes ε_ph implicitly through α₀ = 1069833/264500.

---

## 3-phase plan

### **Phase 1: ε_ph numerical landscape audit (5 sub-tests)**

**Cel:** Verify ε_ph = 4197/25000 = 0.16788 sympy-exact, confirm M9.1″ refined
provenance, falsify identity candidates 2-5.

**Sub-tests:**
- E1.1: ε_ph = 0.16788 sympy-exact rational (4197/25000)
- E1.2: ψ_ph = 4/(3 + 0.4250) = 1.16788 z M9.1″ null geodesic
- E1.3: ε_ph² = 0.028224 (= 441/15625) sympy exact
- E1.4: F4 chain consistency α₀ = target_shift_F4 / ε_ph² = 4.04489 (drift 0.0038%)
- E1.5: Refinement gain M9.2-D (0.168) → M9.1″ (0.16788) closes 36×

### **Phase 2: ε_ph = ψ_ph − 1 structural decomposition (7 sub-tests)**

**Cel:** Pokazać że ε_ph **structurally** = ψ_ph − 1 z M9.1″ photon-ring
geodesic; falsify all 5 alternative identity candidates; show
α₀ = target_shift/ε_ph² constraint LOCKS ε_ph implicitly.

**Sub-tests:**
- E2.1: ε_ph = ψ_ph − 1 structural identity (sympy exact)
- E2.2: Identity candidate falsifications (1.168/(2π), (ψ_ph−1)/ψ_ph,
  1/(2πκ_TGP), 1/(2β²+δ_F4), φ-FP-derived)
- E2.3: F4 chain locks ε_ph implicitly (α₀ = 1069833/264500 + target_shift = 57/500
  ⟹ ε_ph² = 441/15625 sympy exact)
- E2.4: Heat-kernel a₂ frame consistency (ε_ph² entering ξ_geom · α(α−1) · ε_ph²/2)
- E2.5: M9.1″ vs M9.2-D refinement audit (0.527% F4-strict split = a₂ EFT 1-loop)
- E2.6: NGFP RG-stability of ε_ph² ratio under common β-rescaling
- E2.7: Classification — STRUCTURAL DECOMPOSITION (not algebraic identity)

### **Phase 3: predictions E1-E6 + classification (6 sub-tests)**

**Cel:** Generate 6 falsifiable predictions about ε_ph structural footprint
across TGP sectors; classification PARTIALLY DERIVED.

**Sub-tests:**
- E3.1: E1: ngEHT 2030+ ε_ph precision target 0.05% (independent F4 anchor check)
- E3.2: E2: Cross-sector ε_ph appearance (BH/SC/UV/F4) max drift < 0.5%
- E3.3: E3: Identity-falsification roadmap (5 candidates, all > 5% drift)
- E3.4: E4: ε_ph² as definitional ratio (RG-invariant, LISA 2035+)
- E3.5: E5: ε_ph closure z F4 chain unique (no degenerate solutions)
- E3.6: E6: 5-channel falsification convergence

---

## Verdict gates

- **Phase 1 ≥ 4/5 PASS** → Phase 2 proceeds
- **Phase 2 ≥ 6/7 PASS** → ε_ph DECOMPOSITION confirmed, Phase 3 proceeds
- **Phase 3 ≥ 5/6 PASS** → ε.1 program END, classification PARTIALLY DERIVED

**Total:** 18 sub-tests across 3 phases, ledger 373 → 391.

---

## Status taxonomy expectations

- **ε_ph = ψ_ph − 1** structural decomposition: **DERIVED (refined)**
- Identity candidates (5): **FALSIFIED**
- F4 chain implicit lock: **LOCKED-derivative** (already DERIVED via UV.1)
- Heat-kernel a₂ frame: **STRUCTURAL** (already locked via ξ.1)

---

## Materiał wykonawczy

- **Phase 1:** [`Phase1_setup.md`](Phase1_setup.md) + [`phase1_eps_audit.py`](phase1_eps_audit.py) + [`Phase1_results.md`](Phase1_results.md)
- **Phase 2:** [`Phase2_setup.md`](Phase2_setup.md) + [`phase2_eps_decomposition.py`](phase2_eps_decomposition.py) + [`Phase2_results.md`](Phase2_results.md)
- **Phase 3:** [`Phase3_setup.md`](Phase3_setup.md) + [`phase3_eps_predictions.py`](phase3_eps_predictions.py) + [`Phase3_results.md`](Phase3_results.md)

---

## Cross-references

- [`../op-uv-as-ngfp/Phase3_results.md`](../op-uv-as-ngfp/Phase3_results.md) — UV.1 program END (predecessor)
- [`../op-xi-photon-ring/Phase1_results.md`](../op-xi-photon-ring/Phase1_results.md) — ξ.1 frame audit (ε_ph in ξ-factor)
- [`../op-phase3-uv-completion/Phase3_E_results.md`](../op-phase3-uv-completion/Phase3_E_results.md) — F4 chain α₀ closure
- [`../../INDEX.md`](../../INDEX.md) — master ledger 373 → 391
