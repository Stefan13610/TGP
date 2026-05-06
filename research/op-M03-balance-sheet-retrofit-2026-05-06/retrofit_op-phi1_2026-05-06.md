---
title: "Phase 0 balance sheet retrofit — op-phi1-substrate-action-variational (φ.1)"
date: 2026-05-06
parent: "[[README.md]]"
type: balance-sheet-retrofit
cycle_audited: op-phi1-substrate-action-variational
cycle_path: "[[../op-phi1-substrate-action-variational/README.md]]"
auditor: Claudian
classification: STRUCTURAL
tgp_owner: research/op-M03-balance-sheet-retrofit-2026-05-06
tags:
  - phase0
  - balance-sheet-retrofit
  - retrospective
  - phi1
  - phase3-medium-risk
  - positive-example
related:
  - "[[../op-phi1-substrate-action-variational/Phase3_results.md]]"
  - "[[../op-upsilon1-closure-cross-family/Phase3_results.md]]"
---

# Phase 0 balance sheet retrofit — φ.1 (op-phi1-substrate-action-variational)

## Metadata cyklu

- **Cykl:** [[../op-phi1-substrate-action-variational/README.md]]
- **Data oryginalnego closure:** 2026-04-30 (Phase 3 PASS, 18/18)
- **Data retrofit:** 2026-05-06
- **Auditor:** Claudian (M03 Phase 3, medium-risk #4 — referenced jako AXIOM przez ω.1)
- **Klasyfikacja końcowa:** **STRUCTURAL** ★

## 1. Co cykl twierdzi że robi

Z [[../op-phi1-substrate-action-variational/Phase3_results.md]] verdict
(18/18 PASS, "FULL CONVERGENCE 4/4"):

> "Substrate-action variational principle S[X] = ∫½(∂_μ ln X)² d⁴x
> AXIOM-LIFTED. Triple closure π.1 + τ.1 + υ.1 sympy-LOCK reproduced
> from single Lagrangian. RG-stability of closure factor = experimental
> signature. 5 alt-actions falsified. **Universal closure law promoted
> DERIVED → AXIOM.**"

Główne claims:

- **C1**: ℒ = ½(∂_μ ln X)² **unique** action density (5 alt-actions FALSIFIED)
- **C2**: closure(X) = (X_ref/X_obs)^(1/N_gen) → derived z EL + boundary at L/N_gen
- **C3**: X → λX **Noether scale-symmetry**, J^μ = ∂^μ(ln X)
- **C4**: 6 isotopes (⁷Be, ³⁷Ar, ⁵¹Cr, ⁷¹Ga, ⁹⁸Mo, ¹³⁷Cs) z single substrate-action
- **C5**: RG-stability R_TGP invariant — experimental falsifier (FRIB 2030+)
- **C6**: Universal closure law **AXIOM-LIFTED** z DERIVED status

## 2. Phase 0 balance sheet (CALIBRATION_PROTOCOL §2)

### 2.1 External inputs

```
- BEST 2022 ⁷¹Ga R_obs = 0.8084 ± 0.0295            [arXiv:2205.05421]
- SAGE-Ar 1996 R_obs = 0.79 ± 0.10 (³⁷Ar)           [Abdurashitov+ 1996]
- N_gen = 3 charged-lepton generations              [PDG, K-taxonomy TGP-internal]
- Z proton numbers (atomic isotopes)                [periodic table]
```

### 2.2 Structural axioms (TGP-internal LOCKED)

```
- υ.1 universal closure (X_ref/X_obs)^(1/N_gen)     [op-upsilon1, DERIVED 2026-04-29]
- π.1 NME M ∝ (76/A)^(1/3)                          [op-pi1, DERIVED]
- τ.1 overlap f = (Z_a/Z_t)^(1/3)                    [op-tau1, DERIVED]
- N_gen = 3 (cascade primality)                      [empirical + K-taxonomy + υ.1]
- X → λX gauge symmetry of substrate field           [φ.1 builds on this]
- Boundary sampling at x = L/N_gen                   [N_gen-cascade subdivision postulate]
```

### 2.3 Derived outputs

```
- O1: ℒ = ½(∂ ln X)²                                 (unique w klasie X^a(∂X)^b)
- O2: EL: □(ln X) = 0 → ln X linear in x             (sympy direct)
- O3: closure = (X_ref/X_obs)^(1/N_gen)              (z EL + boundary L/N_gen)
- O4: J^μ = ∂^μ(ln X), ∂_μ J^μ = 0                   (Noether scale-current)
- O5: R_TGP RG-stable across energy scales           (scale-current conservation)
- O6: 6 isotope predictions z single ℒ                (lab 2030+ + post-confirms)
```

### 2.4 Tautology test (CRITICAL)

**O1 (Lagrangian uniqueness):** w klasie ∫X^a(∂X)^b z dim-4 +
scale-invariance + Noether → only a=-2, b=2 gives closure-form EL.
NIE tautology w sensie redukcji do tożsamości; **structural uniqueness
within stated EFT class** (dim-4, scale-symmetric, X→λX gauge).
**PASS** (within class).

**O2 (EL → linear ln X):** standard Euler-Lagrange from L=½(∂ln X)²
gives □(ln X) = 0; in 1D static: ln X linear w x. Mathematical identity
**PASS** (correct derivation).

**O3 (closure derivation):** ln X linear z EL + boundary X(0)=X_ref,
X(L)=X_obs → ln X(x) = ln X_ref + (x/L)(ln X_obs - ln X_ref). Sampling
przy x = L/N_gen daje closure = (X_ref/X_obs)^(1/N_gen). To jest
**mathematical identity** dla wybranego sampling point.

**Krytyczne pytanie:** dlaczego sampling przy x = L/N_gen? Odpowiedź
program.md: "N_gen=3 cascade subdivision: boundary X(0)=X_ref,
X(L)=X_obs, intermediate sampling at x=L/N_gen". To jest **postulate**
N_gen-cascade subdivision. Sampling point NIE jest derived from L lub
EL — jest external boundary condition consistent z υ.1 form. φ.1 więc
**reprodukuje υ.1** z L+EL+postulated-sampling, nie *deriviuje*
1/N_gen exponent z first principles.

**Niemniej:** L+EL+sampling-postulate jest **mathematically clean
unified description** kompatybilna z N_gen=3 cascade primality
triple-locked (empirical + K-taxonomy + υ.1). To NIE jest tautology
(L jest substantive); jest **structural lift z DERIVED → AXIOM-level
unification** consistent w sensie variational principle.

**O5 (RG-stability):** Noether scale-current J^μ = ∂^μ(ln X) z
∂_μ J^μ = □(ln X) = 0 — direct mathematical consequence chosen L. RG
-invariance closure factor jest **prediction** (testable experimental
falsifier).

**Werdykt tautology test:** PASS dla derivations (correct math); O3
**relies on N_gen-sampling postulate** which itself is consistent
z υ.1 form but NOT first-principles derived w φ.1.

### 2.5 Falsifiability test (CRITICAL)

**C5 (RG-stability):** "ANY R_TGP RG-running observed → φ.1 falsified"
— concrete experimental falsifier. FRIB 2030+, KamLAND-Zen + CUPID-LSM
joint, BEST extensions. **PASS**.

**6 isotope predictions explicit:**
| Isotope | R_TGP | Status |
|---------|-------|--------|
| ⁷¹Ga | 0.7751 | BEST 2022 1.13σ POST-CONFIRMED |
| ³⁷Ar | 0.7619 | SAGE-Ar 1996 0.28σ POST-CONFIRMED |
| ⁵¹Cr | 0.7693 | GALLEX/SAGE post-confirmed |
| ⁷Be | 0.6531 | lab 2030+ |
| ⁹⁸Mo | 0.7793 | FRIB 2030+ |
| ¹³⁷Cs | 0.8014 | CUPID-LSM 2030+ |

Each value falsifiable z explicit experimental horizon. **PASS**.

**4th gen / sterile constraints:** "direct evidence of 4th chiral
generation would invalidate φ.1 N_gen=3 axiom" — concrete falsifier
EWPM bounds. **PASS**.

**Werdykt falsifiability test:** PASS — multiple concrete experimental
falsifiers across 4 channels.

### 2.6 Independent-path cross-validation (CRITICAL for DERIVED)

**O3 closure form — 1 main path:** L → EL → linear ln X → boundary
sampling at L/N_gen → (X_ref/X_obs)^(1/N_gen). **Single mathematical
path**.

**Cross-confirmation:** φ.1 closure form **reproduces** υ.1 sympy-LOCK
exact. To NIE jest niezależna derivation z różnych domen — to jest
**re-derivation** tej samej form z variational principle (lift z
DERIVED → AXIOM-level w sensie unified Lagrangian description).

**Alt-action falsification (5 alts):** ∫½(∂X)², ∫½(∂X)²-V, ∫X(∂X)²,
∫(∂ln X)⁴, general X^a(∂X)^b → **structural falsification map**
within class dim-4 + scale-inv. To jest **conditional uniqueness**
(within class, NOT against arbitrary actions). PASS w stated class.

**N_gen=3 cascade primality:** triple-locked z 3 sources:
- (a) empirical: 3 charged-lepton generations PDG
- (b) K-taxonomy: K_tot = K_e·K_μ·K_τ factorization (TGP-internal)
- (c) υ.1: 1/N_gen exponent (TGP-internal, derived earlier)

3 paths but **2 z 3 są TGP-internal** — to NIE jest 3 niezależne
domains. **PARTIAL multi-source** convergence.

**4-channel "FULL CONVERGENCE" (F3.5):** 1) Lagrangian uniqueness
(in class), 2) Noether scale-symmetry (mathematical consequence chosen
L), 3) υ.1 reproduction (POST-CONFIRMED, czyli reproduces pre-existing),
4) RG-stability prediction (LIVE 2030+). Channels 1+2+3 są
**internal consistency checks** dla wybranej L (analog ω.1 pattern).
Channel 4 = single LIVE prediction.

**Identical pattern do ω.1:** "FULL CONVERGENCE" framing dla 3
internal mathematical checks + 1 LIVE prediction = **borderline
promotional**.

**Werdykt independent-path:** STRUCTURAL — uniqueness w stated EFT
class + sympy-exact reproduction υ.1 form + multi-isotope consistency.
**NOT** multi-domain niezależna DERIVATION; primary path 1 + internal
re-derivation.

## 3. Audit gate checklist

```
☑ Phase 0 balance sheet exists (this file)
☑ Tautology test PASS (correct derivations); N_gen-sampling postulate flagged
☑ Falsifiability test PASS (R_TGP RG-stability + 6 isotope predictions + 4th gen exclusion)
☐ Independent-path cross-validation PARTIAL — 1 main mathematical path; uniqueness conditional na EFT class
☑ Alt-scan ≥4 candidates (5 alt-actions FALSIFIED w stated class)
☑ NIE used post-hoc structural motivations (N_gen=3 z empirical PDG; X→λX standard scale-symmetry)
☑ NIE circular anchor (φ.1 reproduces υ.1; nie self-reference)
☑ NIE inheriting drift > parent × 5× (uses υ.1 form direct, no drift accumulation)
```

**7/8 ☑ + 1 ☐** — independent-path PARTIAL → max status STRUCTURAL.

"FULL CONVERGENCE" framing borderline (analog ω.1) — flagged for
Phase 5 registry annotation.

## 4. Klasyfikacja końcowa

| Klasa | Spełnia? |
|-------|----------|
| DERIVED FULL | NO — single mathematical path; 1/N_gen sampling jest postulate, nie derived |
| DERIVED CONDITIONAL | partial — closure form sympy-exact z L+EL+sampling; conditional na N_gen-cascade subdivision postulate |
| **STRUCTURAL** | **YES** — Lagrangian uniqueness w stated EFT class + sympy-exact closure reproduction + 6 isotope consistency map |
| ANSATZ | NO — concrete falsifiers + structural derivation w stated class |
| NUMEROLOGICAL | NO — no multi-candidate fit z minimum drift; alt-actions FALSIFIED przez structural arguments |
| TAUTOLOGY | NO — outputs są substantive physics |

**Final verdict:** **STRUCTURAL** ★

**Strukturalne cechy positive example:**
- **Honest "POST-CONFIRMED" vs "DERIVED" usage:** F3.1 calls ⁷¹Ga
  1.13σ "POST-CONFIRMED via φ.1 variational principle" — accurate
  (φ.1 reproduces υ.1 form, observational status low σ).
- **LIVE 2030+ explicit** dla RG-stability i 4 isotopes — no over-claim.
- **Alt-action falsification map** w stated class substantive (5 alts +
  general X^a(∂X)^b form).
- **N_gen=3 triple-locked** explicit acknowledged jako empirical +
  K-taxonomy + υ.1; NIE pretends to derive z first principles.

**"AXIOM-LIFTED" framing:** rhetorical promotion z DERIVED → AXIOM-level
via variational principle. **Substantive within Lagrangian unification
sense** (single L description), ale **NIE means N_gen=3 derived** —
N_gen=3 stays as triple-locked input. Phase 5 registry annotation:
"AXIOM-lifted = unified Lagrangian description; N_gen=3 stays cascade-
primality input".

**"FULL CONVERGENCE 4/4" framing:** identical borderline pattern do
ω.1 — 3 internal consistency checks + 1 LIVE prediction. Phase 5
re-characterize: "1 unified L description passes 3 internal mathematical
consistency tests + RG-stability LIVE 2030+ falsifier".

**Phase 6 gate compliance:**
1. ✓ Phase0_balance.md exists (this file)
2. ✓ "AXIOM-LIFTED" framing flagged jako Lagrangian-unification, NIE z
   first-principles N_gen=3 derivation
3. ✓ Brak constructed criterion (no multi-candidate selection)
4. ⚠ "FULL CONVERGENCE" framing borderline promotional — flagged Phase 5
5. ✓ Brak sympy-rationalization-as-DERIVED (sympy LOCK jest direct EL of L)

## 5. Comparison ze status oryginalnym

| Element | Original claim | Retrofit verdict |
|---------|----------------|------------------|
| Status YAML | `status: PASS`, "FULL CONVERGENCE 4/4", 18/18 PASS, "AXIOM-LIFTED" | STRUCTURAL — Lagrangian uniqueness w stated class + closure reproduction + 6 isotope consistency. "AXIOM-LIFTED" = unified description, NOT N_gen=3 derivation. |
| Counter | "Cumulative ledger: 661 + 18 = 679" | Stays w sensie sub-tests substantive; "FULL CONVERGENCE" downgraded do "structural-form + LIVE prediction" |
| Sub-tests | 5+7+6 = 18/18 PASS | Substantive (Phase 1 form-scan rigorous, Phase 2 sympy LOCK + alt falsification, Phase 3 isotope predictions + RG-stability) |
| Independence | "Triple closure π.1+τ.1+υ.1 reproduced; 4-channel FULL CONVERGENCE" | Re-characterized: π.1+τ.1+υ.1 są **same closure form** z 3 isotope sectors; 4 channels = 3 internal + 1 LIVE |

## 6. Recommended action

- [x] **NO-OP klasy** — STRUCTURAL classification z honest reporting OK
- [x] **DOWNGRADE annotation** w Phase 5 PREDICTIONS_REGISTRY:
      "AXIOM-LIFTED" → "Lagrangian-unification AXIOM-level; N_gen=3 stays
      cascade-primality triple-locked input; closure form reproducible z
      L+EL+sampling-postulate"
- [x] **"FULL CONVERGENCE 4/4" annotation:** "1 unified Lagrangian
      passes 3 internal consistency checks + 1 LIVE 2030+ RG-stability
      prediction"
- [ ] CRITIQUE — nie wymaga
- [ ] CASCADE_AUDIT — **flag dla downstream**: ω.1 uses φ.1 jako AXIOM.
      ω.1 retrofit już acknowledges this; NIE wymaga additional cascade.
- [ ] CORE_IMPACT — none

## 7. Notes

**Strukturalna pozycja φ.1 w hierarchii TGP:**

φ.1 jest **lifting cycle** — promuje υ.1 (closure form, DERIVED) do
AXIOM-level przez unified Lagrangian description. To NIE jest nowa
fizyka, jest **mathematical consolidation**.

Pattern: L = ½(∂ ln X)² → EL → linear ln X → boundary sampling →
closure form. Single L description **reproduces** triple closure
(π.1 NME + τ.1 overlap + υ.1 universal). To jest substantive sense
"AXIOM-LIFTED".

**Comparison z ω.1 + BH.1 patterns:**

| Cykl | "FULL CONVERGENCE" framing | Real structure |
|------|----------------------------|----------------|
| BH.1 | n=2 z 3 niezależnych constraints (Z₂ + MICR-2 + non-overkill) | Multi-domain DERIVED |
| ω.1 | 4 channels = 3 internal mathematical checks + 1 LIVE PARTIAL | Structural |
| **φ.1** | **4 channels = 3 internal checks + 1 LIVE prediction** | **Structural (analog ω.1)** |

φ.1 fits ω.1 pattern: substantive Lagrangian uniqueness w stated class +
"FULL CONVERGENCE" framing borderline promotional. Phase 5 annotation
recommended.

**N_gen=3 sampling postulate vs derivation:**

closure form (X_ref/X_obs)^(1/N_gen) wymaga sampling przy x = L/N_gen.
φ.1 program treats this as **N_gen-cascade subdivision** — natural
boundary condition consistent z empirical N_gen=3. ALE samo N_gen=3
NIE pochodzi z φ.1; jest input from K-taxonomy + empirical PDG +
υ.1 (which derived 1/N_gen exponent).

**Triple-locked N_gen=3** jest **multi-source consistency check**, NOT
single-cycle derivation. Phase 5 registry annotation should clarify:
- φ.1 derives closure-form **given** N_gen
- N_gen=3 is empirical-cascade-primality input

This is honest accounting; cycle nie over-claim N_gen-derivation.

**ω.1 dependency:** ω.1 substrate-EM coupling uses φ.1 substrate-action
**directly as AXIOM** (½f_X²(∂ln X)² term). φ.1 STRUCTURAL classification
nie problematizuje ω.1 — substrate-action remains valid Lagrangian-
unification level even if N_gen=3 is input not derived. ω.1 retrofit
already acknowledges this dependency.

## 8. Cross-references

- [[../op-phi1-substrate-action-variational/README.md]] — cykl audited
- [[../op-phi1-substrate-action-variational/Phase3_results.md]] — main claims source
- [[../op-phi1-substrate-action-variational/program.md]] — 3-phase plan + Lagrangian hypothesis
- [[../op-upsilon1-closure-cross-family/Phase3_results.md]] — υ.1 closure DERIVED (parent)
- [[../op-pi1-bb0nu-nme-isotope/]] — π.1 NME (input)
- [[../op-tau1-closure-overlap-coulomb/]] — τ.1 overlap (input)
- [[../op-omega1-substrate-em-coupling/]] — ω.1 uses φ.1 as AXIOM
- [[retrofit_op-omega1_2026-05-06.md]] — ω.1 retrofit (consistent classification)
- [[README.md]] — M03 master plan
- [[audit_log.md]] — appended 2026-05-06 (Phase 3 #4)
- [[tracker.md]] — status updated to DONE_STRUCTURAL
- [[../../meta/CALIBRATION_PROTOCOL.md]] — protocol source
- [[../../meta/CALIBRATION_GATE_ENFORCEMENT.md]] — Phase 6 gate
