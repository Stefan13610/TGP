---
title: "Phase 0 balance sheet retrofit — op-omega1-substrate-em-coupling (ω.1)"
date: 2026-05-06
parent: "[[README.md]]"
type: balance-sheet-retrofit
cycle_audited: op-omega1-substrate-em-coupling
cycle_path: "[[../op-omega1-substrate-em-coupling/README.md]]"
auditor: Claudian
classification: STRUCTURAL
tgp_owner: research/op-M03-balance-sheet-retrofit-2026-05-06
tags:
  - phase0
  - balance-sheet-retrofit
  - retrospective
  - omega1
  - phase3-medium-risk
  - positive-example
related:
  - "[[../op-omega1-substrate-em-coupling/Phase3_results.md]]"
  - "[[../op-phi1-substrate-action-variational/Phase3_results.md]]"
---

# Phase 0 balance sheet retrofit — ω.1 (op-omega1-substrate-em-coupling)

## Metadata cyklu

- **Cykl:** [[../op-omega1-substrate-em-coupling/README.md]]
- **Data oryginalnego closure:** 2026-04-30 (Phase 3) — z honest downgrade 2026-05-01
- **Data retrofit:** 2026-05-06
- **Auditor:** Claudian (M03 Phase 3, medium-risk #2)
- **Klasyfikacja końcowa:** **STRUCTURAL** ★ (positive example, self-corrected 2026-05-01)

## 1. Co cykl twierdzi że robi

Z [[../op-omega1-substrate-em-coupling/Phase3_results.md]] verdict
(18/18 PASS, "FULL CONVERGENCE"):

> "Substrate ↔ EM axion-like topological coupling
> ℒ = -¼F² + ½f_X²(∂ln X)² + (g/4)(ln X)FF̃ generuje 4 falsifiable
> experimental channels. 3 alt-couplings cross-channel FALSIFIED.
> ω.1 closes G3 gap z em_from_substrate."

Główne claims:

- **C1**: ℒ_ω.1 axion-like form jest **structurally unique** non-trivial coupling (gauge inv + scale inv + parity-odd channel)
- **C2**: Modified Maxwell + substrate EOMs sympy LOCK
- **C3**: 3 alt-couplings cross-channel FALSIFIED (minimal/dilaton/gradient)
- **C4**: g LOCK — **4 candidates** (κ_TGP, α_em, 1/(2π), η_chir) compatible z PVLAS bound
- **C5**: β = 0.34±0.09° "LIVE PARTIAL candidate" (downgraded 2026-05-01 z "POST-CONFIRM")
- **C6**: F·F̃ ∝ E·B, NOT B² alone (corrected 2026-05-01)

## 2. Phase 0 balance sheet (CALIBRATION_PROTOCOL §2)

### 2.1 External inputs

```
- PVLAS-IV / OSQAR-II bound g/f_a < 6.6·10⁻¹¹ GeV⁻¹      [Della Valle+ 2016]
- Planck PR4 + ACT 2024 β = 0.34 ± 0.09° (~3.8σ hint)    [Eskilt 2024]
- Webb/Murphy 2003-17 NULL on α_em(z) variation          [no Δα/α observed]
- Magnetar B ~ 10¹⁵ G + E_rot ~ 10¹⁰⁻¹² V/m              [pulsar pole geometry]
- Substrate scale candidate M_TGP ~ 10¹⁶-10¹⁹ GeV       [GUT/Planckian, 3 orders range]
- α_em ≈ 1/137.036                                       [CODATA]
```

### 2.2 Structural axioms (TGP-internal LOCKED)

```
- φ.1 substrate-action ½f_X²(∂ ln X)²              [op-phi1 Phase 3, AXIOM]
- Scale invariance X → λX (ln X → ln X + ln λ)     [φ.1 Noether]
- Gauge invariance U(1)                            [standard]
- F·F̃ = total divergence (topological density)     [standard, parity-odd]
- κ_TGP = 2.012 (cross-sector charge candidate)    [TGP-SC v2, XS.1]
```

### 2.3 Derived outputs

```
- O1: ℒ_ω.1 = -¼F² + ½f_X²(∂ln X)² + (g/4)(ln X)F·F̃   (Phase 1 form-scan)
- O2: Modified Maxwell ∂_ν F^νμ = g F̃^μν ∂_ν(ln X)    (sympy EL, Phase 2)
- O3: Modified substrate □(ln X) = (g/4f_X²) F·F̃        (sympy EL, Phase 2)
- O4: Δχ(η) = (g/2) ∫ ∂(ln X)/∂η dη                     (Phase 3, CMB rotation)
- O5: g ∈ {κ_TGP, α_em, 1/(2π), η_chir} — 4 candidates  (Phase 1 W1.4, NOT locked)
- O6: f_a ~ M_TGP ~ 10¹⁶-10¹⁹ GeV — 3 orders range      (Phase 3 W3.1)
```

### 2.4 Tautology test (CRITICAL)

**O1 (axion form unique):** Argument: minimal coupling A^μ ∂_μ(ln X) =
total div in Lorenz gauge (gauge-trivial); F·F̃ = ε^μνρσ F_μν F_ρσ is
unique parity-odd dim-4 operator coupling F to ln X scale-invariant.
Substytucja: alternative forms (dilaton (ln X)F², gradient (∂ln X)²F²)
violate scale-symmetry preservation OR introduce mass-scale OR are dim-8.
Output **NIE kasuje się tożsamościowo** — uniqueness comes z structural
constraints (gauge + scale + parity + dim ≤ 4), explicit physical
content. **PASS** (within EFT dim-4 truncation).

**O2/O3 (modified EOMs):** Direct Euler-Lagrange of ℒ_ω.1 (sympy LOCK
verified Phase 2 W2.1+W2.2). Pure derivation — **PASS**.

**O5 (g LOCK):** 4 candidates z theoretical compatibility z PVLAS bound.
**No selection criterion** — author honestly leaves open. NIE tautology
(each candidate testable separately), ale also NIE single-value LOCK.
**PARTIAL** — claim "g LOCK candidates" is honest list, not closure.

**Werdykt tautology test:** PASS dla form + EOMs (głównych derivations);
PARTIAL dla g specific value (multi-candidate honest disclosure).

### 2.5 Falsifiability test (CRITICAL)

**C1 (axion form):** falsifier per W3.5 — 3 alt-couplings cross-channel
FALSIFIED: dilaton z Webb/Murphy 2003-17 NULL on α_em(z); minimal
gauge-trivial; gradient dim-8 cosmologically irrelevant. Each alt has
distinct experimental signature, NULL detection forces axion as unique
viable form **w obrębie EFT dim-4 + scale-symmetric class**. **PASS**
(within stated class).

**C5 (β CMB birefringence):** Planck PR4 + ACT 2024 β = 0.34±0.09° at
~3.8σ. Author **honestly downgraded** 2026-05-01 (ψ.1.v2 critique cycle):
"POST-CONFIRM" → "LIVE PARTIAL candidate". Falsifier: SO 2027+,
LiteBIRD 2029+ corroboration to upgrade or null. Theoretical TGP range
β ~ 0–1.5° (g·O(1)·ln(1100)) — band wider than Planck error 0.09° but
NOT order-of-magnitude wider; falsifiable z 5σ requirement. **PASS**.

**C2/C3 (E·B magnetar):** F·F̃ ∝ E·B prediction (corrected 2026-05-01),
testable via FAST/SKA 2030+ pulsar timing anomalies near magnetic
poles. **PASS**.

**Quasar W3.4:** Δχ(z) ∝ ln(1+z) distinct from Faraday (1+z)⁻²·∫B_∥ —
SKA 2030+ z > 4 polarimetry. **PASS**.

**Werdykt falsifiability test:** PASS — każdy główny channel ma
konkretny experimental falsifier 2027-2030+.

### 2.6 Independent-path cross-validation (CRITICAL for DERIVED)

**O1 form uniqueness — 1 structural path** (gauge + scale + parity + dim≤4
EFT class). **NOT** multi-path derivation, ale uniqueness argument
substantive **w obrębie stated class**. STRUCTURAL grade.

**Cross-channel falsification map (W3.5):** dilaton/minimal/gradient
FALSIFIED z 3 niezależnych domains (Webb/Murphy α_em null,
gauge-triviality, dim-8 suppression). To jest **falsification cross-
validation**, NOT independent derivation paths. STRUCTURAL grade.

**4-channel "FULL CONVERGENCE" (W3.6):** 1) gauge inv = mathematical
identity for chosen form; 2) scale inv preservation = derivation step;
3) EOM sympy = direct EL of L_ω.1; 4) CMB obs ~3.8σ hint = LIVE PARTIAL.
Channels 1+2+3 są **internal consistency checks** dla wybranej form
(NOT 3 niezależne paths derivujące unikalność). Channel 4 = single
observational hint (LIVE PARTIAL). **"FULL CONVERGENCE"** framing jest
**promotional** — properly characterized: "1 form passes 3 mathematical
consistency tests + 1 observational LIVE PARTIAL".

**g specific value — 0 paths** (4 candidates open).

**Werdykt independent-path:** STRUCTURAL (uniqueness in narrow class +
falsification cross-validation), NOT DERIVED (no multi-path convergent
derivation z różnych domen).

## 3. Audit gate checklist

```
☑ Phase 0 balance sheet exists (this file)
☑ Tautology test PASS (form + EOMs); PARTIAL (g value)
☑ Falsifiability test PASS (β + E·B + quasar Δχ + alt-coupling falsification map)
☐ Independent-path cross-validation PARTIAL — 1 structural path (uniqueness in EFT class), NOT multi-path
☑ Alt-scan ≥4 candidates with discrimination (4 coupling forms; 3 FALSIFIED z niezależnych obs/structural arguments)
☑ NIE used post-hoc structural motivations (gauge+scale+parity standard symmetries)
☑ NIE circular anchor (φ.1 substrate-action z niezależnego cyklu)
☑ NIE inheriting drift > parent × 5× (uses φ.1 as axiom, no drift)
```

**7/8 ☑ + 1 ☐** — independent-path PARTIAL → max status STRUCTURAL.
g specific value NOT locked → status STRUCTURAL even jeśli wszystkie
inne PASS-uje.

**"FULL CONVERGENCE" framing oryginalny** = over-promotional dla
mathematical consistency checks (gauge+scale+EOM są 3 *internal* checks,
NIE 3 *independent paths*).

## 4. Klasyfikacja końcowa

| Klasa | Spełnia? |
|-------|----------|
| DERIVED FULL | NO — independent-path PARTIAL; g multi-candidate; β LIVE PARTIAL only |
| DERIVED CONDITIONAL | NO — wymaga ≥2 niezależne paths convergent; aktualnie 1 structural path |
| **STRUCTURAL** | **YES** — form uniqueness w EFT class dim-4 + cross-sector falsification map robust; β prediction concrete falsifier (SO/LiteBIRD 2027-2029) |
| ANSATZ | NO — concrete falsifiers + structural uniqueness argument |
| NUMEROLOGICAL | NO — no multi-candidate fit z minimum drift selection (g list is honest open set, NOT fitted) |
| TAUTOLOGY | NO — outputs są substantive physics |

**Final verdict:** **STRUCTURAL** ★

**Strukturalne cechy positive example:**
- **Self-correction 2026-05-01** (ψ.1.v2 critique): "POST-CONFIRM" → "LIVE PARTIAL candidate" — explicit honest downgrade with rationale.
- **Mathematical correction 2026-05-01**: "B² sourcing" → "F·F̃ ∝ E·B" — physics correctness over rhetorical convenience.
- **g multi-candidate honest disclosure**: 4 candidates listed without artificial selection criterion (avoiding κ.1-style constructed criterion antipattern).
- **Cross-channel falsification map** (W3.5): 3 alt-couplings FALSIFIED z niezależnych domains — substantive scientific content, NOT mere consistency.
- "FULL CONVERGENCE" framing pozostaje **promotional** dla 4-channel structure
  (3 internal mathematical checks + 1 LIVE PARTIAL observational hint), ale
  NIE compromituje classification — wszystkie sub-tests substancjalne.

**Phase 6 gate compliance:**
1. ✓ Phase0_balance.md exists (this file)
2. ✓ Brak status promotion (Phase 1 PASS → Phase 2 PASS → Phase 3 explicit downgrade 2026-05-01, opposite direction)
3. ✓ Brak constructed criterion (g list honest open, NIE selected)
4. ⚠ "FULL CONVERGENCE" framing borderline promotional — flagged for Phase 5 registry annotation
5. ✓ Brak sympy-rationalization-as-DERIVED (EOMs sympy direct from L_ω.1, structural derivation)

## 5. Comparison ze status oryginalnym

| Element | Original claim | Retrofit verdict |
|---------|----------------|------------------|
| Status YAML | `status: PASS`, "FULL CONVERGENCE 4/4", 18/18 PASS | STRUCTURAL — form unique w narrow EFT class; cross-channel falsification map robust; observational status LIVE PARTIAL only |
| Counter | "Cumulative ledger: 679 + 18 = 697 post-ω.1" | Stays w sensie sub-tests substancjalne; "FULL CONVERGENCE" annotation downgraded do "structural-form + observational LIVE PARTIAL" |
| Sub-tests | 5+7+6 = 18/18 PASS | Sub-tests substantive (NOT mechanical): Phase 1 form-scan rigorous, Phase 2 sympy LOCK direct EL, Phase 3 cross-channel falsification map z niezależnych obs |
| Independence | "4 channels FULL CONVERGENCE" | Re-characterized: 3 internal mathematical checks + 1 LIVE PARTIAL observational hint; NOT 4 niezależne derivation paths |

**Drift discipline:** uses φ.1 substrate-action as direct axiom (no
drift); β prediction band 0–1.5° vs Planck 0.34±0.09° = order O(1)
match within natural range — falsifiable not over-fitted.

## 6. Recommended action

- [x] **NO-OP klasy** — STRUCTURAL classification z honest reporting OK
- [x] **DOWNGRADE annotation** w Phase 5 PREDICTIONS_REGISTRY refactor:
  "FULL CONVERGENCE" → "structural-form-locked + observational LIVE
  PARTIAL". Counter entry stays jako STRUCTURAL with explicit experimental
  falsifier (β SO/LiteBIRD 2027-2029).
- [ ] CRITIQUE — nie wymaga (cycle ALREADY self-corrected 2026-05-01)
- [ ] CASCADE_AUDIT — nie wymaga (downstream σ.1, τ.2, ζ.1 są research-
  track candidates, NOT promoted DERIVED)
- [ ] CORE_IMPACT — none (cycle is research-level, not core LaTeX)

## 7. Notes

**Honest self-correction pattern (2026-05-01) jest exemplary:**

Przed: "PARTIAL POST-CONFIRM" + "B² sourcing" claims.

Po ψ.1.v2 critique cycle: "LIVE PARTIAL candidate" (3.8σ ≠ confirmed) +
"F·F̃ ∝ E·B" (rotation-induced parallel E required, not pure B alone).

To jest dokładnie wzorzec self-correction discipline z
[[../../meta/CALIBRATION_PROTOCOL.md]] §4 ("mark-as-unproven, NIE
rollback"). Sub-tests PASS preserved; tylko **interpretacja statusu**
downgraded. ω.1 dostarcza pre-Phase-6 example tej discipline w
naturalnej formie (przed binding date).

**g multi-candidate vs κ.1 constructed criterion — kluczowe odróżnienie:**

ω.1 g ∈ {κ_TGP, α_em, 1/(2π), η_chir} — list is **open set** of theoretically
compatible candidates without **selection criterion**. PASS per Phase 6
("brak constructed criterion").

κ.1 ρ̄_num=11 z 4 sympy-exact paths + post-hoc "denom-num pairing" rule
**constructed by select** C0 = NUMEROLOGICAL.

Identical surface pattern (4 candidates), opposite Phase 6 verdict due
to **selection criterion presence**. ω.1 stops at "any of these is
compatible, pending future cycle"; κ.1 went further to "this specific
combination wins by this rule" without falsifier for the rule.

**"FULL CONVERGENCE" framing borderline:**

Author counted: 1) gauge inv (manifest), 2) scale inv preservation
(W1.3 derivation), 3) EOM sympy LOCK (W2.1+W2.2), 4) CMB obs (~3.8σ
LIVE PARTIAL). Channels 1+2+3 są mathematical/derivation-internal dla
chosen L_ω.1 form. Channel 4 = single observational hint. Properly
characterized: "1 axion-like form passes 3 internal consistency checks
+ 1 observational LIVE PARTIAL". Phase 5 registry refactor should
annotate accordingly.

## 8. Cross-references

- [[../op-omega1-substrate-em-coupling/README.md]] — cykl audited
- [[../op-omega1-substrate-em-coupling/Phase3_results.md]] — main claims source
- [[../op-omega1-substrate-em-coupling/Phase2_results.md]] — sympy LOCK derivations
- [[../op-omega1-substrate-em-coupling/program.md]] — 3-phase plan + form-scan
- [[../op-phi1-substrate-action-variational/Phase3_results.md]] — φ.1 axiom source
- [[../op-cross-sector-charge/]] — XS.1 κ_TGP candidate source
- [[README.md]] — M03 master plan
- [[audit_log.md]] — appended 2026-05-06 (Phase 3 #2)
- [[tracker.md]] — status updated to DONE_STRUCTURAL
- [[../../meta/CALIBRATION_PROTOCOL.md]] — protocol source
- [[../../meta/CALIBRATION_GATE_ENFORCEMENT.md]] — Phase 6 gate
