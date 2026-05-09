---
title: "op-S07-alternative-f-psi-derivation-2026-05-09 — alternative f(ψ) post-M9.1'' falsification"
date: 2026-05-09
type: research-cycle
priority: P1_CRITICAL
parent: "[[../../audyt/S07_M911_derivation/README.md]]"
target: "S07 audit realization — find alternative f(ψ) ansatz post-(4-3ψ)/ψ falsification"
classification: TGP_BLOCKER_RESOLUTION
status: ACTIVE — Phase 0
folder_status: active
predecessor_cycles:
  - "[[../op-Phi-vacuum-scale-2026-05-09/]] (matter sector dual-V CLOSED)"
  - "[[../op-V-canonical-consistency-audit-2026-05-09/]]"
  - "[[../op-MAG-Phase5-V-reference-clarification-2026-05-09/]]"
  - "[[../op-dual-V-structure-clarification-2026-05-09/]]"
  - "[[../op-Phase5-MAG-erratum-2026-05-09/]]"
  - "[[../op-Phi0-spatial-variation-predictions-2026-05-09/]]"
falsification_source:
  - "[[../op-ppE-mapping/Phase1.5_G_SPA_lock.md]] (G_SPA=48 sympy LOCK 5/5)"
  - "[[../op-GWTC3-reanalysis/Phase2_RERUN_2026-05-09_corrected_beta.md]] (5.02σ RULED OUT)"
g0_reference:
  - "[[../op-g0-r3-from-canonical-projection/]] (M9.1'' canonical V_M911 LOCK)"
tags:
  - S07-realization
  - alternative-f-psi
  - M911-falsification-resolution
  - gravity-sector-blocker
  - first-principles-derivation
  - 1PN-exact-GR
  - 2PN-bounded
  - GWTC-3-compliant
---

# op-S07-alternative-f-psi-derivation-2026-05-09

## Mission

Realize **S07 audit** ([[../../audyt/S07_M911_derivation/README.md]]) post-2026-05-09
M9.1'' specific (4−3ψ)/ψ form **5σ falsification** by GWTC-3.

**Goal:** find alternative ansatz `f(ψ)` (gdzie g_tt = −c₀² f(ψ)) który:
1. Preserves passing tests: 1PN exact GR, mass spectrum (V-independent A_tail), κ
2. Satisfies GWTC-3 constraint: |β_ppE^(b=-1)| ≤ 0.78 (1σ) at η=1/4
3. Wynika z **first-principles** derivation (NIE empirical ansatz)

**Falsifier (within scope):** jeśli ŻADEN alternative f(ψ) z naturally
parametrized family nie satisfies all constraints + first-principles
derivation, M9.1'' approach (gravity emergence z f(ψ) form na M9.1''-style
structure) jest FALSIFIED. TGP wymaga radically different gravity emergence
mechanism (post-S07 fallback).

## Status pre-cycle (2026-05-09)

### Falsification context

| Quantity | M9.1'' value | GWTC-3 constraint | σ-tension |
|---|---|---|---|
| f(ψ) ansatz | (4−3ψ)/ψ | — | — |
| Δα_3 (g_tt @ U³) | -5/6 | — | — |
| G_SPA (SPA chain prefactor) | 48 sympy-exact | — | — |
| β_ppE^TGP^(b=-1) at η=1/4 | -15/4 ≈ -3.75 | \|β\| ≤ 0.78 (1σ) | **5.02σ RULED OUT** |
| BF_TGP/GR (GWTC-3 combined ~90 BBH) | 3.55·10⁻⁶ | — | OVERWHELMING GR |

### Status sectorów post-falsification

| Sektor | Status | Independent of f(ψ)? |
|---|---|---|
| Gravity (V_M9.1'') | **(4-3ψ)/ψ FALSIFIED 5σ** | NO — gravity sector blocker |
| Matter (V_orig) | ✅ OK (dual-V framework) | YES (separate sector) |
| EFT Phi_0 | ✅ scale-dependent free param | YES |
| Spin (op-SPIN-SU2) | ✅ V-independent (RP² topology) | YES |
| Particle masses (P4) | ✅ A_tail mechanism | YES (V-independent) |

## Phase plan

### Phase 0: Balance sheet + hard constraints C1-C10 formalization

**Deliverables:**
- ✅ This README (Phase 0 setup) — DONE
- Phase0_balance.md — full balance sheet
- Phase0_constraints_C1_C10_sympy.py — C1-C10 sympy formalization
- Gate: ≥8/10 constraints sympy-formalized

### Phase 1: Reconnaissance — enumerate candidate f(ψ) families

**Plausible parametrizations:**
- **Polynomial:** f(ψ) = a₀ + a₁ψ + a₂ψ² + ...
- **Rational:** f(ψ) = P(ψ)/Q(ψ) (rational; (4-3ψ)/ψ był specjalny przypadek)
- **Exponential:** f(ψ) = exp(-aψ) lub f(ψ) = e^(-2U(ψ))
- **Hybrid:** kombinacje
- **Constraint-natural:** start from C1-C10 algebraically forced form

**Deliverables:**
- Phase1_setup.md
- Phase1_candidate_enumeration_sympy.py
- Phase1_results.md
- Gate: ≥5 candidate families enumerated z preliminary C2 (PPN exact) check

### Phase 2: Constraint testing per candidate (sympy heavy)

For each surviving candidate from Phase 1:
- **C1**: EOM compatibility (α=2 vacuum)
- **C2**: 1PN expansion → γ_PPN = β_PPN = 1 EXACT (sympy)
- **C3**: 2PN expansion → Δα_3 numeric
- **C4**: Compute G_SPA(f) → check |Δα_3·G_SPA| ≤ 8.32 at η=1/4
- **C5**: Newton limit κ = 3/(4·Φ_0) preserved
- **C6**: Mass spectrum (A_tail, P22-style verification)
- **C7**: Vacuum stability m_sp² > 0
- **C8**: Hyperbolicity natural BH cutoff (ψ_h horizon)
- **C9**: √(-g_eff) consistency
- **C10**: Dual-V matter sector independence

**Deliverables:**
- Phase2_setup.md (per-candidate test plan)
- Phase2_constraint_testing_sympy.py (per candidate)
- Phase2_results.md
- Gate: ≥1 candidate satisfies all C1-C10

### Phase 3: First-principles derivation (if Phase 2 PASS)

If candidate satisfies all C1-C10, derive z fundamentu:
- **Path A**: H_Γ coarse-graining → continuum action → forced f(ψ)
- **Path B**: variational uniqueness given C1+C2+C8 → forced f(ψ)
- **Path C**: forced by C1-C10 algebraically → unique f(ψ)

**Deliverables:**
- Phase3_setup.md
- Phase3_derivation_sympy.py (Path A/B/C)
- Phase3_results.md
- Gate: at least one Path A/B/C closure

### Phase 4: Framework integration (if Phase 3 PASS)

- Replace V_M9.1'' formula in sek08a, sek08c
- Update G.0 P22 mass spectrum (verify still works)
- Update PPN P23 (verify γ=β=1 still exact)
- Update PREDICTIONS_REGISTRY M911-P1/P2/P3 z new ansatz
- Update TGP_FOUNDATIONS.md §3.5 + dual-V annotation

### Phase FINAL: ABSOLUTE BINDING gate

Classify outcome:
- **DERIVED**: alternative f(ψ) found + first-principles derived + all C1-C10
- **STRUCTURAL_CONDITIONAL**: alternative satisfies C1-C10 but only empirical
- **EARLY_HALT**: NO alternative satisfies C1-C10 → M9.1'' approach FALSIFIED structurally

## Hard constraints C1-C10 (binding)

| # | Constraint | Source | Sympy testable? |
|---|---|---|---|
| C1 | α=2 vacuum Φ-EOM preserved | T-D-uniqueness (TGP foundation) | YES |
| C2 | 1PN exact GR: γ_PPN=β_PPN=1 EXACT | G.0 P23 sympy LOCK 5/5 | YES (PPN expansion) |
| C3 | GWTC-3 2PN: \|β_ppE^(b=-1)\| ≤ 0.78 (1σ) at η=1/4 | observational 2026-05-09 | YES (SPA chain) |
| C4 | \|Δα_3·G_SPA\| ≤ 8.32 (η=1/4) | derived z C3 | YES |
| C5 | Newton limit: κ = 3/(4·Φ_0) preserved | G.0 P32 INVARIANT | YES |
| C6 | Mass spectrum: m_μ/m_e≈207, m_τ/m_e≈3477 | G.0 P22 (V-independent A_tail) | YES (analytical) |
| C7 | Vacuum stability: m_sp² > 0 (no tachyon) | G.0 P21 | YES |
| C8 | Hyperbolicity natural BH cutoff: ψ_h horizon | structural (M9.1'' philosophy) | YES |
| C9 | Volume element √(-g_eff) consistent z f(ψ) | sek08a | YES |
| C10 | Dual-V matter sector independence: V_orig nie affected | dual-V clarification 2026-05-09 | YES (sectoral) |

## Methodological constraints (CALIBRATION_PROTOCOL binding)

Per [[../../meta/CALIBRATION_PROTOCOL.md]]:

| Anti-pattern | Forbidden | Detection |
|---|---|---|
| 1 | Multi-candidate fit z minimum drift selection | Treat candidates equally; no post-hoc selection |
| 2 | Constructed criterion post-hoc | All criteria pre-declared in C1-C10 |
| 3 | Drift hardening fitted corrections | NO empirical fudge factors |
| 4 | Algebraic re-arrangement masquerading as derivation | First-principles MUST trace to H_Γ or variational |
| 5 | Definitional tautology | NO circular "this satisfies because we defined it that way" |
| 6 | Sympy-rationalization "DERIVED" without first-principles | Path A/B/C MUST produce constructive proof |

**Honest reporting MANDATORY.** EARLY_HALT z honest acknowledgment > forced
DERIVED.

## Sympy verification standards

Per CALIBRATION_PROTOCOL:
- check() helper z PASS/FAIL counter
- ≥5 sympy tests per candidate f(ψ)
- W końcowym raporcie X/Y PASS prominently
- PYTHONIOENCODING=utf-8 dla Polish characters
- 4-level verification dla key claims (sympy + hand + numerical + alternative route)

## Time budget (multi-session expected)

- **Session 1 (current):** Phase 0 + Phase 1 reconnaissance
- **Session 2-3:** Phase 2 (constraint testing, sympy heavy)
- **Session 4:** Phase 3 (derivation if candidate found)
- **Session 5:** Phase 4 + Phase FINAL (integration + close)

**Honest expectation:** może wymagać iterowania z user'em (jak op-Phi-vacuum-scale).

## Probability assessment (post-handoff)

| Outcome | Probability |
|---|---|
| Alternative f(ψ) found + first-principles derived | 25-35% |
| Alternative satisfies C1-C10 but only empirical | 30-40% |
| NO alternative → M9.1'' approach FALSIFIED structurally | 20-30% |
| Multi-candidate (several work) → fundamental ambiguity | 10-15% |

## Open questions (carryover od previous sessions)

- Czy K(ψ)=ψ⁴ (T-D-uniqueness) pozostaje invariant pod f(ψ) update?
  (P22 mass spectrum używa A_tail mechanism który jest V-independent.)
- Czy √(-g) form forsuje jednoznacznie f(ψ)? (G.0 P21 used (4-3ψ)/ψ assumption.)
- Czy hyperboliczna struktura (BH horizon) jest **wymaganą cechą** f(ψ),
  czy tylko **opcjonalną**?

## Cross-references

- [[../../audyt/S07_M911_derivation/README.md]] — S07 audit issue source
- [[../../audyt/EXTERNAL_REVIEW_2026-05-06.md]] §EXT-3 — external review
- [[../../core/sek08a_akcja_zunifikowana/sek08a_akcja_zunifikowana.tex]]
- [[../../core/sek08c_metryka_z_substratu/sek08c_metryka_z_substratu.tex]]
- [[../op-ppE-mapping/Phase1.5_G_SPA_lock.md]] — falsification source (G_SPA=48)
- [[../op-GWTC3-reanalysis/Phase2_RERUN_2026-05-09_corrected_beta.md]] — 5σ ruled out
- [[../op-g0-r3-from-canonical-projection/]] — G.0 closure (V_M911 LOCK)
- [[../../TGP_FOUNDATIONS.md]] §3.5 — dual-V structure
- [[../../HANDOFF_NEXT_SESSION_S07_alternative_f_psi.md]] — handoff origin
