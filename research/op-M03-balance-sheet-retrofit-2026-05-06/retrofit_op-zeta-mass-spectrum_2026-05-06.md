---
title: "Phase 0 balance sheet retrofit — op-zeta-mass-spectrum (real ζ.1)"
date: 2026-05-06
parent: "[[README.md]]"
type: balance-sheet-retrofit
cycle_audited: ζ.1 (op-zeta-mass-spectrum)
cycle_path: "[[../op-zeta-mass-spectrum/README.md]]"
auditor: Claudian
classification: STRUCTURAL
tgp_owner: research/op-M03-balance-sheet-retrofit-2026-05-06
tags:
  - phase0
  - balance-sheet-retrofit
  - zeta1
  - PMNS
  - lambda-C-source
  - iota-cascade-source
---

# Phase 0 balance sheet retrofit — ζ.1 (op-zeta-mass-spectrum)

## Metadata cyklu

- **Cykl:** [[../op-zeta-mass-spectrum/README.md]]
- **Data closure:** 2026-04-29 (ζ.1.Phase2 PASS 7/7)
- **Data retrofit:** 2026-05-06
- **Auditor:** Claudian
- **Klasyfikacja końcowa:** **STRUCTURAL** (consistent z author "PARTIALLY DERIVED")

## 1. Co cykl twierdzi że robi

Z [[../op-zeta-mass-spectrum/Phase2_results.md]]:

> PMNS angles all closed z group structure: sin²θ₁₂ = 1/3 (S₃ democratic),
> sin²θ₂₃ = 1/2 (Z₂ atmospheric + K(ν)=1/2), sin²θ₁₃ = λ_C²/2 (cross-sector
> Cabibbo lock); 5 alternative parameterizations FALSIFIED.
> Classification: PARTIALLY DERIVED (refined).

Główne claims:

- **C1:** sin²θ₁₂ = 1/3 (S₃ ⊂ GL(3,𝔽₂) democratic; drift 8.58%)
- **C2:** sin²θ₂₃ = 1/2 (Z₂ atmospheric + K_ν = 1/2; drift 12.59%)
- **C3:** sin²θ₁₃ = λ_C²/2 (cross-sector Cabibbo; drift 15.57%)
- **C4:** λ_C = 165/167 = 0.22550 (GL form factor)
- **C5:** "PARTIALLY DERIVED (refined); max drift 15.6% < 20% gate"

## 2. Phase 0 balance sheet (CALIBRATION_PROTOCOL §2)

### 2.1 External inputs

```
- NuFit 5.3 sin²θ₁₂ = 0.307 ± 0.013     [4.2% band]
- NuFit 5.3 sin²θ₁₃ = 0.022 ± 0.0007    [3.2% band]
- NuFit 5.3 sin²θ₂₃ = 0.572 ± 0.024     [4.2% band]
- PDG λ_C = 0.22500 ± 0.00067           [0.30% band]
```

### 2.2 Structural axioms (TGP-internal)

```
- N_gen = 3                               [why_n3 + R3 ODE]
- K_ν = 1/2 (Majorana chirality)          [θ.1 K-taxonomy — STRUCTURAL]
- S₃ ⊂ GL(3,𝔽₂) democratic permutation    [structural group theory — pre-derived]
- Z₂ atmospheric reflection μ-τ swap      [structural symmetry]
- λ_C = 165/167 (GL form factor)          [STATUS: GL(3,𝔽₂) derivation — structural]
```

### 2.3 Derived outputs

```
- sin²θ₁₂ = 1/3 (S₃ democratic)
- sin²θ₂₃ = 1/2 (Z₂ + K_ν)
- sin²θ₁₃ = λ_C²/2 = 0.0254
- λ_C = 0.22550 (165/167 GL form factor)
```

## 2.4 Tautology test

#### sin²θ₁₂ = 1/3 (S₃ democratic)

S₃ democratic permutation eigenstate (1,1,1)/√3 daje |U_e2|² = 1/3.
**Structural argument** z group theory.

**Werdykt:** ✅ **PASS** (S₃ axiom pre-derived).

#### sin²θ₂₃ = 1/2

Z₂ atmospheric + K_ν = 1/2 — structural argument.

**Werdykt:** ✅ **PASS** (Z₂ + K_ν pre-derived).

#### sin²θ₁₃ = λ_C²/2

```
sin²θ₁₃ = λ_C²/2 = (0.2255)²/2 = 0.02543
```

Cross-sector Cabibbo lock — λ_C inherited z GL(3,𝔽₂) form factor 165/167.
λ_C = 165/167 jest **structural** (group theory), not fitted.

**Werdykt:** ✅ **PASS** (λ_C structural).

### 2.5 Falsifiability test

**Aktualny problem (jak ι.1):**

| Angle | TGP | NuFit 5.3 (1σ) | Tension |
|-------|-----|----------------|---------|
| sin²θ₁₂ | 0.333 | [0.294, 0.320] | +2.0σ outside |
| sin²θ₁₃ | 0.0254 | [0.0213, 0.0227] | +5.4σ outside (!) |
| sin²θ₂₃ | 0.500 | [0.548, 0.596] | -3.0σ outside |

**Wszystkie 3 angles już *outside* 1σ.**

ALE: ζ.1 explicit klasyfikuje jako "**PARTIALLY DERIVED (refined)**" — NIE
"DERIVED FULL". Author używa "20% gate dla zeroth-order" — accommodating
ale **honest** disclosure że to zeroth-order, "1-loop RG corrections" needed.

To kontrastuje z ι.1 który claims "DERIVED FULL CASCADE 7/7" przy tych
samych drifts.

**Werdykt falsifiability:** ⚠ **PARTIAL** — drifty już outside 1σ, ale
honest "PARTIALLY DERIVED" classification.

### 2.6 Independent-path cross-validation

**Path 1:** S₃ democratic (group theory) → sin²θ₁₂ = 1/3
**Path 2:** Z₂ + K_ν (chirality) → sin²θ₂₃ = 1/2
**Path 3:** λ_C cross-sector → sin²θ₁₃ = λ_C²/2

**Convergence:** 3 paths są **strukturalne**, każdy z innego źródła
(group theory, chirality, GL form factor) — **independent paths**.

**Werdykt independent-path:** ✅ **PASS** — 3 niezależne strukturalne
arguments.

## 3. Audit gate checklist

```
☑ Phase 0 balance sheet exists
☑ Tautology test PASS (S₃, Z₂, λ_C wszystkie pre-derived structural)
⚠ Falsifiability test PARTIAL (drifty outside 1σ, honest "20% gate" dla zeroth-order)
☑ Independent-path PASS (3 niezależne paths)
☑ Alt-scan: 5 candidates FALSIFIED at 18-100% drift
☑ NIE post-hoc structural motivations (group theory pre-derived)
☑ NIE circular anchor
☑ NIE inheriting drift > parent × 5×
```

**7 ☑ + 1 ⚠ z 8 ⇒ STRUCTURAL.**

## 4. Klasyfikacja końcowa

| Klasa | Spełnia? |
|-------|----------|
| DERIVED FULL | NO — drifts 8-15% outside 1σ wymagają 1-loop corrections |
| DERIVED CONDITIONAL | NO — wymaga 1-loop completion (μ.1 attempted via fitted corrections — NUMEROLOGICAL) |
| **STRUCTURAL** | **YES** — group theory + λ_C all structural; honest "PARTIALLY DERIVED (refined)" |
| ANSATZ | partial — zeroth-order hypothesis ALE strukturalna baza solid |
| NUMEROLOGICAL | NO — group theory derivation legit |
| TAUTOLOGY | NO |

**Final verdict:** **STRUCTURAL** (consistent z author "PARTIALLY DERIVED (refined)";
5-th positive example honest reporting).

## 5. Comparison ze status oryginalnym

| Element | Original claim | Retrofit verdict |
|---------|----------------|------------------|
| Status | "PARTIALLY DERIVED (refined); max drift 15.6% < 20% gate" | **STRUCTURAL** (consistent) |
| Honest reporting | "1-loop RG corrections needed" + "20% gate zeroth-order" | ★ **POSITIVE EXAMPLE** (5-ty) |
| ι.1 promotion ζ.1 → DERIVED | via mixing-operator framework | **REVERSED** — ι.1 jest ANSATZ, μ.1 NUMEROLOGICAL; ζ.1 stays STRUCTURAL |

## 6. Recommended action

- ☑ **NO-OP** — ζ.1 honestly classified
- ☑ **MINOR ANNOTATION** w PREDICTIONS_REGISTRY: "STRUCTURAL"
- ☐ CRITIQUE — niepotrzebne (positive cycle)
- ☑ **CASCADE_AUDIT impact** ↔ **ι.1, μ.1 reverse**:
  - ι.1 cascade promotion ζ.1 "PARTIALLY DERIVED → DERIVED" via mixing-operator: **WITHDRAWN**
  - μ.1 cascade promotion ζ.1 "via drift hardening": **WITHDRAWN**
  - ζ.1 stays as **honest baseline** — drifty 8-15% explicit zeroth-order
- ☐ CORE_IMPACT — brak

## 7. Notes

**ζ.1 jest 5-tym POSITIVE EXAMPLE honest reporting.**

| Cykl | Honest classification | M03 verdict |
|------|----------------------|-------------|
| δ.1 | "PARTIAL POSITIVE" | STRUCTURAL ✓ |
| δ.2 | "Level B PARTIAL POSITIVE" | STRUCTURAL_CONDITIONAL ✓ |
| γ.1 | "POSITIVE z H5 multi-anchor" | STRUCTURAL ✓ |
| XS.1 | "PARTIALLY DERIVED" | STRUCTURAL ✓ |
| **ζ.1** | **"PARTIALLY DERIVED (refined); zeroth-order"** | **STRUCTURAL ✓** |

**5 z 12 cykli audytowanych są honest reports** (42%) — stable baseline.

**Critical insight: ζ.1 jest BASELINE for ι.1+μ.1 family.**

ζ.1 honestly says: "drifty 8-15%, zeroth-order, 1-loop corrections needed".

ι.1 took ζ.1 outputs i claimed "DERIVED FULL CASCADE" via mixing-operator
framework (post-hoc construction).

μ.1 took ι.1 + κ.1 outputs i claimed "drift hardening" via fitted corrections.

**ζ.1 → ι.1 → μ.1 cascade chain**: honest baseline → over-claiming via
post-hoc framework → over-claiming via fitted hardening.

**Wniosek:** ζ.1 baseline jest SOLID; ι.1 + μ.1 over-claiming są **not
legitimate promotions** of ζ.1.

**Pozytywne strony ζ.1:**
- Group-theoretic derivation 3 PMNS angles (S₃, Z₂, λ_C structural)
- Honest "20% gate zeroth-order" disclosure
- 5 algebraic alts falsified
- λ_C structural via GL form factor 165/167

**Krytyczne strony ζ.1:**
- Drifty 8-15% już outside NuFit 1σ band — even honestly classified
  zeroth-order, claim "PARTIALLY DERIVED" wymaga reservation
- 1-loop RG corrections "needed" — open problem (NIE legitimately closed
  by μ.1 fitted hardening)

## 8. Cross-references

- [[../op-zeta-mass-spectrum/README.md]] — cykl audited
- [[../op-zeta-mass-spectrum/Phase2_results.md]] — group-theoretic derivation
- [[retrofit_op-iota_2026-05-06.md]] — ι.1 (cascade downstream, ANSATZ)
- [[retrofit_op-mu_2026-05-06.md]] — μ.1 (cascade downstream, NUMEROLOGICAL)
- [[retrofit_op-theta_2026-05-06.md]] — θ.1 (K_ν source)
- [[../op-cross-sector-charge]] — XS.1 (NIE ζ.1 — separate cycle)
- [[README.md]] — M03 master plan
- [[audit_log.md]] — running log
- [[tracker.md]] — status updated do DONE_STRUCTURAL
- [[../../meta/CALIBRATION_PROTOCOL.md]]
