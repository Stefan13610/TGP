---
title: "Phase 0 balance sheet retrofit — op-eta-wolfenstein (η.1)"
date: 2026-05-06
parent: "[[README.md]]"
type: balance-sheet-retrofit
cycle_audited: η.1 (op-eta-wolfenstein)
cycle_path: "[[../op-eta-wolfenstein/README.md]]"
auditor: Claudian
classification: STRUCTURAL
tgp_owner: research/op-M03-balance-sheet-retrofit-2026-05-06
tags:
  - phase0
  - balance-sheet-retrofit
  - eta1
  - wolfenstein
  - CKM
  - high-risk-resolved
---

# Phase 0 balance sheet retrofit — η.1 (op-eta-wolfenstein)

## Metadata cyklu

- **Cykl:** [[../op-eta-wolfenstein/README.md]]
- **Data oryginalnego closure:** 2026-04-29 (η.1.Phase3 PASS 6/6, program END 18/18)
- **Data retrofit:** 2026-05-06
- **Auditor:** Claudian
- **Klasyfikacja końcowa:** **STRUCTURAL** (downgrade z claim "PARTIALLY DERIVED (refined)")

## 1. Co cykl twierdzi że robi

Z [[../op-eta-wolfenstein/Phase2_results.md]] werdykt:

> Triple (A, ρ̄, η̄) = (64/81, 11/78, 5/14) LOCKED z drift < 0.05% per
> component; 5/5 alternatives FALSIFIED at 0.5% threshold; cross-sector
> denom-family analysis documented (η̄ denom 14 shares prime 7 z K_up = 7/8
> numerator); V_ub_TGP refined cascade drift 8.93%; J_TGP = 2.93·10⁻⁵ drift 4.51%;
> classification PARTIALLY DERIVED (refined).

Główne claims:

- **C1:** Wolfenstein triple `(A, ρ̄, η̄) = (64/81, 11/78, 5/14)` sympy-exact
- **C2:** A = 64/81 = `K_up_denom²/N_gen⁴` (cross-sector cascade z θ.1 i N_gen)
- **C3:** Cross-sector denom-prime sharing (78 ↔ K_lepton denom 3; 14 ↔ K_up num 7)
- **C4:** 5/5 alternative triples FALSIFIED at 0.5% threshold
- **C5:** V_ub_TGP = 0.00348 (drift 8.93% vs PDG 0.00382)

## 2. Phase 0 balance sheet (CALIBRATION_PROTOCOL §2)

### 2.1 External inputs (PDG)

```
- PDG A = 0.790 ± 0.012        [3 sig figs, 1.5% band]
- PDG ρ̄ = 0.141 ± 0.020        [3 sig figs, 14% band ⚠ WIDE]
- PDG η̄ = 0.357 ± 0.014        [3 sig figs, 4% band]
- PDG λ_C = 0.22500 ± 0.00067  [4 sig figs, 0.30% band]
- PDG |V_ub| = 0.00382 ± 0.00010 [3 sig figs, 2.6% band]
- PDG sin(2β) = 0.699 ± 0.017   [3 sig figs, 2.4% band]
- PDG |V_td/V_ts| = 0.205 ± 0.006 [3 sig figs, 2.9% band]
```

### 2.2 Structural axioms (TGP-internal LOCKED)

```
- N_gen = 3                                  [why_n3 + R3 ODE]
- K_up_denom = 8 (z K_up = 7/8)              [θ.1 STRUCTURAL — retrofit_op-theta_2026-05-06.md]
- K_up_num = 7                                [θ.1 STRUCTURAL]
- K_lepton_num = 2, K_lepton_denom = 3       [θ.1 K-taxonomy]
- λ_C = 0.22550                               [ζ.1: GL form factor 165/167; STATUS: NIE-AUDYTOWANE]
```

**Numerator anchors (research-track open):**
- Numerator A = 64 = K_up_denom² ✓ (forced)
- Numerator ρ̄ = 11 — **research-track κ.1/ι.1** (NIE pre-derived; OPEN)
- Numerator η̄ = 5 — **research-track κ.1/ι.1** (NIE pre-derived; OPEN)

### 2.3 Derived outputs (the cycle claims)

```
- Output 1: A = 64/81 = 0.79012             (PDG 0.790, drift 0.016%)
- Output 2: ρ̄ = 11/78 = 0.14103             (PDG 0.141, drift 0.018%)
- Output 3: η̄ = 5/14 = 0.35714              (PDG 0.357, drift 0.040%)
- Output 4: V_ub = 0.00348                   (PDG 0.00382, drift 8.93%)
- Output 5: J = 2.93·10⁻⁵                    (drift 4.51%)
- Output 6: sin(2β) = 0.7090                 (drift 1.43%)
- Output 7: |V_td/V_ts| = 0.2098             (drift 2.33%)
```

### 2.4 Tautology test (CRITICAL)

#### Output 1: A = 64/81

```
A = 64/81
  = K_up_denom² / N_gen⁴ (z η.2 B2.1 + cross-sector cascade)
  = 8² / 3⁴ ✓ structural
```

**Tautology check:** A = 64/81 wykorzystuje pre-derived axiomy
(K_up_denom z θ.1 STRUCTURAL, N_gen z R3 ODE) — niezależna informacja.

**Werdykt A:** ✅ **PASS** (forced by structural axioms).

⚠ **Caveat:** η.2 derives `81 = N_gen⁴` PO η.1 ustaleniu A = 64/81. Czyli
**reverse engineering:** η.1 fitted (64/81), η.2 explains why (cross-sector
cascade). To **research-track derivation** post-hoc, nie pre-derivation.

#### Output 2: ρ̄ = 11/78

```
Numerator 11: research-track κ.1/ι.1 (Phase3 H4)
Denominator 78 = 2·N_gen·B²_up_num = 2·3·13 (z η.2 B2.2)

Tautology check:
- 78 forced by cross-sector cascade (η.2 derived) ✓
- 11 NOT derived w η.1 ani η.2 — "research-track κ.1/ι.1" OPEN
- ρ̄ = 11/78 jest **partially derived** (denom forced, numerator open)
```

**Werdykt ρ̄:** ⚠ **PARTIAL** — denom forced, but **numerator 11 open**.

⚠ **PDG band ρ̄ jest 14%** (PDG 0.141 ± 0.020). TGP claim 11/78 = 0.141
mieści się w paśmie 0.13% (drift 0.018%). **Jakikolwiek ułamek 0.121-0.161
pasuje** — to jest **multi-candidate scenario** podobny do θ.1 K_down!

Alternatywne ułamki w paśmie:
- 11/78 = 0.14103 (TGP best, drift 0.018%)
- 9/64 = 0.14063 (drift 0.31%)
- 13/92 = 0.14130 (drift 0.21%)
- 10/71 = 0.14085 (drift 0.13%)

Wszystkie w 14% PDG band → **tylko 11/78 selected by minimum drift**.

#### Output 3: η̄ = 5/14

```
Numerator 5: research-track κ.1/ι.1
Denominator 14 = K_up_num · K_lepton_num = 7·2 (z η.2 B2.3)

Tautology check:
- 14 forced by cross-sector cascade ✓
- 5 NOT derived — "research-track κ.1/ι.1" OPEN
- η̄ = 5/14 partially derived
```

**Werdykt η̄:** ⚠ **PARTIAL** — denom forced, **numerator 5 open**.

#### Werdykt overall tautology test

⚠ **PASS QUALIFIED** — denoms forced przez cross-sector cascade (η.2),
ale **numerators (11, 5) są research-track κ.1/ι.1 dependency**.
Output ρ̄ 11/78 ma multi-candidate scenario w paśmie 14% PDG.

### 2.5 Falsifiability test (CRITICAL)

**5 alt-triples falsified at 0.5% threshold (T2.5):**

| Alternative | max drift | FALSIFIED |
|---|---|---|
| C1: (4/5, 1/7, 5/14) | 1.32% | YES |
| C2: (3/4, 1/7, 1/π) | 10.84% | YES |
| C3: (8/9, 1/7, 5/14) | 12.52% | YES |
| C4: Wolfenstein-2002 | 39.01% | YES |
| C5: PDG-2012 | 7.09% | YES |

**Krytyczne:** 0.5% threshold jest **mniejsze niż PDG band 14% dla ρ̄**!
To znaczy że alt-triples z drift 0.5%-14% **NIE są falsified** ani przez
TGP threshold ani przez PDG.

**Realistyczne alt-scan w PDG band:**
- (64/81, **11/78**, 5/14) — TGP best, drift 0.04%
- (64/81, **9/64**, 5/14) — drift ~0.3% (within 14% PDG band, NOT falsified by PDG)
- (64/81, **10/71**, 5/14) — drift ~0.13% (within band)
- (4/5, 11/78, 5/14) — drift A 1.27% (BORDERLINE, but A drift 1.27% < PDG A band 1.5%!)

**Werdykt falsifiability test:** ⚠ **PARTIAL** —
- η̄ falsifiable (4% PDG band, TGP 0.04% drift, 100× tighter)
- A falsifiable (1.5% PDG band, TGP 0.016%, 90× tighter)
- **ρ̄ NIE falsifiable** (14% PDG band, multi-candidate scenario)

**Future falsifier:** Belle II 2027+ + LHCb Run 4 2030+ projected
σ(ρ̄) ~3% — **may** discriminate; aktualnie nie.

### 2.6 Independent-path cross-validation

**Path 1:** Sympy rationalization (find p/q matching PDG, denom < 100)
**Path 2:** Cross-sector cascade z θ.1 + η.2 — ale η.2 derive denomy
*post-hoc* dla η.1 wartości (reverse engineering).

**Werdykt independent-path:** ⚠ **PARTIAL** —
- Realnie 1 niezależna fizyczna ścieżka (sympy rationalization z PDG)
- η.2 cross-sector cascade jest *back-explanation*, nie pre-derivation

## 3. Audit gate checklist

```
☑ Phase 0 balance sheet exists
⚠ Tautology test PARTIAL (denoms forced, numerators open research-track)
⚠ Falsifiability test PARTIAL (ρ̄ multi-candidate w 14% PDG band)
⚠ Independent-path PARTIAL (η.2 to back-explanation)
⚠ Alt-scan 5 algebraic candidates (not 4 first-principles)
⚠ Post-hoc structural motivations (η.2 derives denoms post-hoc)
☑ NIE circular anchor (numerators are open, not circular)
☑ NIE inheriting drift > parent × 5×
```

**5 ⚠ z 8 ⇒ status max STRUCTURAL.**

## 4. Klasyfikacja końcowa

| Klasa | Spełnia? | Uzasadnienie |
|-------|----------|--------------|
| DERIVED FULL | NO | numerators (11, 5) research-track open |
| DERIVED CONDITIONAL | NO | conditional na κ.1/ι.1, ale te HIGHER risk niż η.1 |
| **STRUCTURAL** | **YES** | denoms forced by cross-sector cascade; sympy-exact; 5/5 alt-falsified |
| ANSATZ | NO | jest stronger niż ansatz (concrete Belle II + LHCb predictions) |
| NUMEROLOGICAL | partial | ρ̄ 11/78 ma multi-candidate scenario, ale denominator forced |
| TAUTOLOGY | NO | output ma niezależną informację (denoms cross-sector) |

**Final verdict:** **STRUCTURAL** (downgrade z "PARTIALLY DERIVED (refined)").

**Caveat ρ̄:** Z perspektywy precyzji obecnego PDG (band 14%), ρ̄ = 11/78
jest *nieodróżnialne* od wielu sąsiednich ułamków. Status `partially
NUMEROLOGICAL`. Po Belle II 2027+ + LHCb Run 4 2030+ z σ(ρ̄) ~3%, to
zostanie discriminated → potential promotion do STRUCTURAL+.

## 5. Comparison ze status oryginalnym

| Element | Original claim | Retrofit verdict |
|---------|----------------|------------------|
| Status | "PARTIALLY DERIVED (refined)" | **STRUCTURAL** (downgrade) |
| Counter PREDICTIONS_REGISTRY | +18 | sub-tests mechanical OK; status downgrade w Phase 5 |
| Sub-tests | 18/18 PASS | mechanically PASS; substancjalnie: numerators open, ρ̄ multi-candidate |
| Independence | claim cross-sector locks | weryfikacja: η.2 cascade derives denoms *post-hoc* |

## 6. Recommended action

- ☐ NO-OP — wymaga downgrade
- ☑ **DOWNGRADE w PREDICTIONS_REGISTRY**: status `PARTIALLY DERIVED → STRUCTURAL`
  w Phase 5 registry refactor; with annotation "ρ̄ multi-candidate w PDG 14% band"
- ☐ CRITIQUE — niepotrzebne (no canceling tautology); ale dodać note o
  numerators 11, 5 research-track open
- ☑ **CASCADE_AUDIT**:
  - **κ.1, ι.1** — numerators 11, 5 source — wymagają audyt (κ.1 już w high-risk queue)
  - **η.2** — derives η.1 denoms; back-explanation dependency dokumentowana
- ☐ CORE_IMPACT — brak

## 7. Notes

**Pozytywne strony η.1:**
- Sympy-exact triple (64/81, 11/78, 5/14)
- Concrete Belle II 2027+ + LHCb Run 4 2030+ falsifiers (4/4 LIVE)
- Cross-sector denom-prime sharing (suggestive structural pattern)
- 5/5 algebraic alt-falsified

**Krytyczne strony η.1:**
- **ρ̄ PDG band 14%** — multi-candidate scenario, numerator 11 selected
  z minimum drift (analog UV.2 K_struct + θ.1 K_down)
- **Numerators 11, 5 są "research-track κ.1/ι.1"** — czyli OPEN; cycle
  claims rely on future cycle resolution
- Cross-sector cascade derivation η.2 jest *back-explanation*, nie
  pre-derivation; η.1 fitted first, η.2 explained why

**Pattern:** η.1 reprezentuje "**structural retrofit**" pattern — TGP wartości
fitted przez sympy rationalization (η.1 Phase 1), potem cross-sector
explanation back-fitted (η.2 Phase 2). Ten flow jest legitimate, ale
**audytowo distinct** od first-principles derivation.

## 8. Cross-references

- [[../op-eta-wolfenstein/README.md]] — cykl audited
- [[../op-eta-wolfenstein/Phase2_results.md]] — main derivation
- [[../op-eta-wolfenstein/Phase3_results.md]] — predictions
- [[retrofit_op-eta2_2026-05-06.md]] — η.2 retrofit (downstream cykl)
- [[retrofit_op-theta_2026-05-06.md]] — θ.1 retrofit (K_up_denom prerequisite)
- [[retrofit_op-eps_2026-05-06.md]] — ε.1 retrofit (137 anchor)
- [[README.md]] — M03 master plan
- [[audit_log.md]] — running log (added 2026-05-06)
- [[tracker.md]] — status updated do DONE_STRUCTURAL
- [[../../meta/CALIBRATION_PROTOCOL.md]]
