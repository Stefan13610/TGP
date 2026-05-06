---
title: "Phase 5 — PREDICTIONS_REGISTRY refactor DRAFT (early start)"
date: 2026-05-06
parent: "[[README.md]]"
type: registry-refactor-draft
tgp_owner: research/op-M03-balance-sheet-retrofit-2026-05-06
tags:
  - phase5
  - registry-refactor
  - epistemic-classification
  - draft
  - early-start
---

# Phase 5 — PREDICTIONS_REGISTRY refactor DRAFT

## Status

**EARLY START:** Po Phase 2 high-risk completion (11/11 + ζ.1 = **12 cykli
audited**), early start dla registry refactor. Ten draft documentuje
rekomendacje dla **per-cycle epistemic classification** based on M03
retrofit verdicts.

**Pełny Phase 5** wymaga ukończenia Phase 3 (medium-risk) + Phase 4
(low-risk). Ten draft = **prepared baseline** dla future agent.

## Cel Phase 5

Per [[../../audyt/M03_balance_sheet_missing/README.md]] §"Rekomendacja"
+ [[../../meta/CALIBRATION_PROTOCOL.md]] §3:

1. Per-cycle epistemic class tagging w PREDICTIONS_REGISTRY
2. Counter rozdzielony per-class (DERIVED FULL / CONDITIONAL / STRUCTURAL /
   ANSATZ / NUMEROLOGICAL / TAUTOLOGY)
3. Predictivity ratio re-derive z retrofitted classifications
4. Cascade reverse impacts dokumentowane

## Per-cycle epistemic classification (12 audited cykli)

### Original "DERIVED FULL CASCADE 7/7 PASS" claims requiring downgrade

| # | Cykl | Original counter contribution | M03 retrofit | New counter contribution |
|---|------|-------------------------------|--------------|--------------------------|
| 1 | **κ.1** | +18 (FULL CASCADE) | **NUMEROLOGICAL** | 0 (research-track only, NOT registry-DERIVED) |
| 2 | **ι.1** | +18 (FULL CASCADE) | **ANSATZ** | 0 (research-track only) |
| 3 | **μ.1** | +18 (FULL CASCADE, 8 free→0 free) | **NUMEROLOGICAL** | 0 (research-track only) |

**Subtotal downgrade:** −54 z master ledger.

### Original "PARTIALLY DERIVED (refined)" claims — STRUCTURAL classification

| # | Cykl | Original status | M03 retrofit | Counter impact |
|---|------|-----------------|--------------|----------------|
| 4 | **ε.1** | "PARTIALLY DERIVED (refined)" | **STRUCTURAL** | mild downgrade (sub-tests OK; status precise) |
| 5 | **η.1** | "DERIVED (refined²)" post-κ.1 | **STRUCTURAL** | reverse-cascade (κ.1 NUMEROLOGICAL); η.1 wartości forced przez sympy rationalization |
| 6 | **ζ.1** | "PARTIALLY DERIVED (refined)" | **STRUCTURAL** | consistent (ι.1+μ.1 promotions WITHDRAWN) |

### Original honest "PARTIAL POSITIVE" — STRUCTURAL/STRUCTURAL_CONDITIONAL

| # | Cykl | Original status | M03 retrofit | Counter impact |
|---|------|-----------------|--------------|----------------|
| 7 | **δ.1** | "PARTIAL POSITIVE" | **STRUCTURAL** ★ | NO downgrade (consistent z honest reporting) |
| 8 | **δ.2** | "Level B PARTIAL POSITIVE" | **STRUCTURAL_CONDITIONAL** ★ | NO downgrade |
| 9 | **γ.1** | "POSITIVE z H5 multi-anchor" | **STRUCTURAL** ★ | NO downgrade; **answers D01 NEEDS N2** |
| 10 | **XS.1** | "PARTIALLY DERIVED" | **STRUCTURAL** ★ | NO downgrade |

### Original "FULL CASCADE 18/18 PASS" — DERIVED_CONDITIONAL OK

| # | Cykl | Original status | M03 retrofit | Counter impact |
|---|------|-----------------|--------------|----------------|
| 11 | **η.2** | "FULL CASCADE LOCKED 18/18" | **DERIVED_CONDITIONAL** | 18 → ~12 (with annotation: conditional na ε.1 + θ.1 prerequisites) |
| 12 | **θ.1** | "PARTIALLY DERIVED" / "STRUCTURAL" | **SPLIT** (K_up STRUCTURAL, K_down NUMEROLOGICAL) | K_up: full count; K_down: research-track only |

## Recommended PREDICTIONS_REGISTRY annotations

Per row insert/update:

### κ.1 entry
```diff
- **Updated 2026-04-30 (κ.1 program END):** 499 cumulative
- (+ κ.1.Phase1 5 + κ.1.Phase2 7 + κ.1.Phase3 6 = 18).
- **HH5 promoted STRUCTURAL HINT → PARTIALLY DERIVED**;
- η.1 cascade comprehensive **DERIVED (refined²)**;
- **CKM matrix 4 free params → 0 free params**.

+ ⚠ M03 RETROFIT 2026-05-06: κ.1 reclassified NUMEROLOGICAL OBSERVATION
+ (research-track only, NOT registry-DERIVED). Mixing-operator framework
+ post-hoc construction; 4 sympy-exact paths convergent dla każdego
+ numeratora; "denom-num level-pairing criterion" constructed by select
+ C0 z multi-candidate set. Analog UV.2 K_struct.
+ Sub-tests 18/18 PASS mechanically OK ALE claim "DERIVED FULL" untenable.
+ Counter contribution: -18 (research-track downgrade).
+ Reverse-cascade: η.1 promotion DERIVED→STRUCTURAL untenable.
+ "CKM 4 → 0 free" claim WITHDRAWN.
+ See research/op-M03-balance-sheet-retrofit-2026-05-06/retrofit_op-kappa_2026-05-06.md
```

### ι.1 entry
```diff
+ ⚠ M03 RETROFIT 2026-05-06: ι.1 reclassified ANSATZ (research-track only).
+ PMNS angles 3/3 outside NuFit 5.3 1σ band:
+ - sin²θ₁₃: +5.4σ outside (TGP 0.0254 vs NuFit 0.022±0.0007)
+ - sin²θ₂₃: -3.0σ outside
+ - sin²θ₁₂: +2.0σ outside
+ "Zeroth-order gate <25%" non-standard accommodating gate.
+ Mixing-operator inherited z κ.1 NUMEROLOGICAL contagion.
+ "8 free → 1 open" claim WITHDRAWN. Counter contribution: -18.
+ See retrofit_op-iota_2026-05-06.md
```

### μ.1 entry
```diff
+ ⚠ M03 RETROFIT 2026-05-06: μ.1 reclassified NUMEROLOGICAL.
+ "Drift hardening" via (1-ρ̄), K_ν/K_up, (1-λ_C·η̄) corrections — wszystkie
+ używają κ.1 NUMEROLOGICAL outputs. Lift factors 21×, 126×, 51× exactly
+ compensate ι.1 zeroth drift = circular fitting.
+ δ_CP dual form (Form A 205° vs Form B 260°, 55° diff) accommodating.
+ "8 free → 0 free post-μ.1" claim WITHDRAWN.
+ Counter contribution: -18.
+ Reverse-cascade: ι.1 stays ANSATZ; "DERIVED (refined²)" promotion untenable.
+ See retrofit_op-mu_2026-05-06.md
```

### θ.1 entry
```diff
+ ⚠ M03 RETROFIT 2026-05-06: θ.1 SPLIT classification.
+ K_up = 7/8 (universal pattern (2+B²)/(2N)) classified STRUCTURAL.
+ **K_down = 37/50 reclassified NUMEROLOGICAL**: 4 candidates within 0.05%
+ PDG band (37/50, 54/73, 71/96, 57/77) — multi-candidate scenario analog UV.2.
+ Cascade impact: B²_down=61/25 (used in η.2 Form B + ι.1 charge cascade)
+ inherits NUMEROLOGICAL contamination.
+ See retrofit_op-theta_2026-05-06.md
```

### ε.1 entry
```diff
+ ⚠ M03 RETROFIT 2026-05-06: ε.1 reclassified STRUCTURAL (downgrade z
+ "PARTIALLY DERIVED (refined)"). Sympy 23/137 sympy-exact OK ALE 137
+ anchor borrowed z α_fine PDG, NIE first-principles derived. F4 chain
+ "second path" jest algebraic re-arrangement (nie niezależna ścieżka).
+ 5 alt-candidates falsified są algebraic identities, nie first-principles
+ physical models. Single path → STRUCTURAL.
+ See retrofit_op-eps_2026-05-06.md
```

### η.1 entry
```diff
+ ⚠ M03 RETROFIT 2026-05-06: η.1 reclassified STRUCTURAL (downgrade z
+ "DERIVED (refined²)" via κ.1). Triple (64/81, 11/78, 5/14) sympy-exact
+ OK; A=K_up_denom²/N_gen⁴ structural. ALE: ρ̄ multi-candidate w 14% PDG band;
+ numerators (11, 5) "research-track κ.1/ι.1" — κ.1 NUMEROLOGICAL contagion.
+ η.2 cascade derives denoms POST-HOC (back-explanation, nie pre-derivation).
+ See retrofit_op-eta-wolfenstein_2026-05-06.md
```

### η.2 entry
```diff
+ ⚠ M03 RETROFIT 2026-05-06: η.2 reclassified DERIVED_CONDITIONAL.
+ 2 niezależne paths (Form A z N_gen + Form B z B²_up - B²_down) sympy-exact.
+ Form A clean; Form B contaminated by θ.1 K_down NUMEROLOGICAL (B²_down=61/25).
+ Conditional na ε.1 (137 anchor) + θ.1 prerequisites.
+ Substancja: α⁻¹(0) = 137.036 prediction nadal solid (concrete falsifier
+ Cs/Rb 2027+).
+ See retrofit_op-eta2_2026-05-06.md
```

### ζ.1 entry (PMNS source)
```diff
+ ⚠ M03 RETROFIT 2026-05-06: ζ.1 STRUCTURAL (consistent z honest "PARTIALLY
+ DERIVED (refined)"). Group theory derivation (S₃, Z₂, λ_C) solid;
+ drifts 8-15% honestly classified zeroth-order; "20% gate" accommodating
+ ALE explicit acknowledged.
+ Reverse-cascade: ι.1+μ.1 promotion ζ.1 "→ DERIVED" WITHDRAWN; ζ.1 stays
+ STRUCTURAL baseline.
+ See retrofit_op-zeta-mass-spectrum_2026-05-06.md
```

### γ.1 entry
```diff
+ ★ M03 RETROFIT 2026-05-06: γ.1 STRUCTURAL (positive example).
+ Φ_eff = 8π pure structural identification. Honest disclosure: Brannen
+ 24.783 = phenomenological α_s fit, NIE structural derivation. Multi-anchor
+ reality acknowledged. Directly answers D01 NEEDS N2 (Brannen formal
+ derivation pending).
+ See retrofit_op-gamma1_2026-05-06.md
```

### XS.1 entry (cross-sector α₀↔κ_TGP)
```diff
+ ★ M03 RETROFIT 2026-05-06: XS.1 STRUCTURAL (consistent z honest "PARTIALLY
+ DERIVED"). F1 single-Φ axiom → cross-sector identity solid.
+ M_BH ≈ M_SC w 1% precision (consistent z experimental). 6-channel
+ falsification roadmap 2027-2035 concrete.
+ See retrofit_op-zeta_2026-05-06.md (note: filename — ten cykl jest XS.1,
+ NIE ζ.1)
```

### δ.1, δ.2 entries
```diff
+ ★ M03 RETROFIT 2026-05-06: δ.1 STRUCTURAL, δ.2 STRUCTURAL_CONDITIONAL
+ (consistent z honest "PARTIAL POSITIVE" / "Level B PARTIAL POSITIVE").
+ NO downgrade — author honest classification matches retrofit.
+ δ.1: g̃ = 5e²/(12π) algebraic identity; "5" multi-candidate disclosure;
+ N_f=5 "strongest candidate" plausible.
+ δ.2: N_f=5 derived consequence z TGP mass ordering (dod F + sek09 CW);
+ 8/8 audit gate criteria PASS.
```

## Counter recalculation (estimate post-retrofit)

**Master ledger pre-M03:** 856 cumulative
- Pre-existing audited 4 (chi.1, UV.2, ω.2, ω.3, λ.1) = -72 contested
- Effective uncontested: 784

**Post-M03 retrofit (12 audited):**
- κ.1 -18 (NUMEROLOGICAL)
- ι.1 -18 (ANSATZ)
- μ.1 -18 (NUMEROLOGICAL)
- η.2 ~-6 (DERIVED_CONDITIONAL annotation)
- θ.1 K_down ~-3 (NUMEROLOGICAL portion)
- η.1 -3 (numerator weakness annotation)
- ε.1 -3 (STRUCTURAL annotation)
- ζ.1 -3 (STRUCTURAL annotation)
- γ.1, XS.1, δ.1, δ.2 = 0 (consistent)

**Estimated subtotal downgrade:** ~−72 entries

**Effective uncontested post-M03 12-cycle audit:** ~712 (down from 784)

**Important caveat:** Phase 3 (15 medium-risk) + Phase 4 (13 low-risk) +
sub-cycles (UV.x, BH.x, SC.x) wymagają additional audits. Final number
może być znacząco niższy (~600-650 estimated).

## Predictivity ratio re-derive

**Pre-M03 LS-8 audit:** ratio = 11/2 = 5.5

**Post-M03 estimate:** ratio może spaść do **3.5-4.5** po:
- Removal of κ.1+ι.1+μ.1 mixing-operator family claims (CKM 4→0, PMNS 3→0, 8→0 → revert)
- Downgrade θ.1 K_down (quark Koide partial)
- Annotation η.1, η.2, ε.1 cascades

**Effective predictions per input** post-M03:
- Pre-M03: 11/2 = 5.5
- Post-M03 estimate: 11/3 = 3.67 (removing 3 spurious "0-free" claims)

This jest **bardziej konserwatywne** ratio, bardziej rzetelne dla peer review.

## Cascade reverse impacts (z M03 audits)

| Cykl downstream | Original promotion claim | Retrofit verdict | Status |
|-----------------|--------------------------|------------------|--------|
| η.1 (post-κ.1) | "PARTIALLY → DERIVED (refined²)" | Promotion WITHDRAWN | η.1 STRUCTURAL |
| ι.1 (post-μ.1) | "ANSATZ → DERIVED (refined²)" | Promotion WITHDRAWN | ι.1 ANSATZ |
| ζ.1 (post-ι.1) | "PARTIALLY → DERIVED (mixing-operator)" | Promotion WITHDRAWN | ζ.1 STRUCTURAL |
| ζ.1 (post-μ.1) | "→ DERIVED via drift hardening" | Promotion WITHDRAWN | ζ.1 STRUCTURAL |

## Plan dla full Phase 5 implementation (future)

1. **Wait for Phase 3 + Phase 4 completion** (~10-15 sesji future)
2. **Compile complete per-cycle classification** (28+ cykli)
3. **Generate single PREDICTIONS_REGISTRY refactor PR** with:
   - Per-row epistemic class tag
   - Counter rozdzielony per-class
   - Cascade reverse annotations
   - Predictivity ratio re-derive (final)
4. **Phase 6 — gate enforcement** (CALIBRATION_PROTOCOL absolute binding)

## Recommended: do NOT execute full registry refactor yet

**Reasoning:** Full Phase 5 wymaga complete audit. Premature edit może:
- Create inconsistencies between annotated (audited) vs un-annotated (pending)
- Require multiple PR-style updates per cycle
- Lose context of "honest reporting baseline" patterns

**Recommended action:** ten draft = **frozen baseline** dla 12 audited
cykli; future agent dodaje incremental dla 16+ pozostałych.

## Cross-references

- [[README.md]] — M03 master plan
- [[tracker.md]] — status każdego cyklu
- [[audit_log.md]] — 12 retrofit entries
- [[../../audyt/M03_balance_sheet_missing/POST_ACTION_UPDATE_2026-05-06.md]]
- [[../../meta/CALIBRATION_PROTOCOL.md]] §3 (audit gate)
- [[../../PREDICTIONS_REGISTRY.md]] — target for future full refactor
- 12 retrofit files: retrofit_op-*_2026-05-06.md
