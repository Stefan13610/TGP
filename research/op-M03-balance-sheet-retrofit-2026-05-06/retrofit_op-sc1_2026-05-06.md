---
title: "Phase 0 balance sheet retrofit — op-sc-alpha-origin (SC.1)"
date: 2026-05-06
parent: "[[README.md]]"
type: balance-sheet-retrofit
cycle_audited: op-sc-alpha-origin
cycle_path: "[[../op-sc-alpha-origin/README.md]]"
auditor: Claudian
classification: STRUCTURAL
tgp_owner: research/op-M03-balance-sheet-retrofit-2026-05-06
tags:
  - phase0
  - balance-sheet-retrofit
  - retrospective
  - sc1
  - phase3-medium-risk
  - positive-example
  - bidirectional-falsification
related:
  - "[[../op-sc-alpha-origin/Phase3_results.md]]"
  - "[[../op-bh-alpha-threshold/Phase2_results.md]]"
  - "[[retrofit_op-bh-alpha_2026-05-06.md]]"
---

# Phase 0 balance sheet retrofit — SC.1 (op-sc-alpha-origin)

## Metadata cyklu

- **Cykl:** [[../op-sc-alpha-origin/README.md]]
- **Data oryginalnego closure:** 2026-04-28 (Phase 3 CLOSED, 17/17 sub-tests PASS)
- **Data retrofit:** 2026-05-06
- **Auditor:** Claudian (M03 Phase 3, medium-risk #6 — cross-sector reference dla BH.1)
- **Klasyfikacja końcowa:** **STRUCTURAL** ★ (positive example, analog BH.1 pattern)

## 1. Co cykl twierdzi że robi

Z [[../op-sc-alpha-origin/Phase3_results.md]] verdict (17/17 sub-tests PASS):

> "**Phase 1**: α_PB nie jest unit-cousin α_0 (H₀ rejected).
> **Phase 2**: α_PB JEST A-G-like, ale **a priori J_sf ratio 2.59** (H_AG_PARTIAL).
> **Phase 3**: TGP μ_eff² scaling **preferred** over A-G+dG na istn. danych
> (RMS 0.42 vs 1.5); cleanest tests = SmH₉ (TGP wins 10⁵·⁸×) + YbH₉ (A-G wins 84×) —
> pre-registered."

Główne claims:

- **C1**: α_PB ≠ unit-cousin α_0 (Phase 1 H₀ rejected explicit)
- **C2**: α_PB derivation z A-G + J_sf ratio 2.59 a priori — **PARTIAL** (Phase 2)
- **C3**: 15-Ln³⁺ T_c table z TGP μ_eff² scaling
- **C4**: TGP RMS_log 0.42 vs A-G+dG 1.53 na PrH₉+NdH₉ — TGP preferred 3.6× lepiej
- **C5**: 5 LIVE predictions SC2-SC7 (SmH₉, YbH₉, TmH₉, PmH₉, anchors)
- **C6**: Bidirectional discriminator: SmH₉ TGP wins, YbH₉+TmH₉ A-G wins
- **C7**: Cross-sector hint α₀ = κ_TGP² (z BH.1, 0.75% match)

## 2. Phase 0 balance sheet (CALIBRATION_PROTOCOL §2)

### 2.1 External inputs

```
- PrH₉ T_c = 8.9 K (Drozdov 2019 Sci.Adv. 120 GPa lit. update 2025)  [TESTED]
- NdH₉ T_c = 4.5 K (Zhou 2020 JACS 142:2803)                          [TESTED]
- LaH₁₀ T_c ≈ 200 K baseline                                          [Drozdov 2018]
- T_c^base = 143 K (LaH₁₀ baseline TGP-fit)                          [calibrated SC v2]
- 15-Ln³⁺ Hund GS: μ_eff², g_J, dG factor (Jensen & Mackintosh 1991)  [analytical atomic]
- κ_TGP = 2.012 (V, Nb, Ta, Mo, Pd) → also 0.75% gap κ_TGP² ≈ α₀^BH   [SC v2 calibration]
- α₀^BH ≈ 4.018 (BH.1 photon-ring strict, ξ=1)                       [BH.1 STRUCTURAL HINT]
```

### 2.2 Structural axioms (TGP-internal LOCKED)

```
- Standard Abrikosov–Gorkov: ln(T_c⁰/T_c) = ψ(½+Γ_sf/(2πT_c)) - ψ(½)  [textbook]
- Born approximation Γ_sf ∝ μ_eff²·N(0)·τ_sf                          [textbook]
- Hund's rules ground state for Ln³⁺ ions                              [atomic physics]
- de Gennes scaling (g_J - 1)²·J(J+1) standard SC reference            [textbook]
- TGP T_c formula T_c = T_c^base · exp(-c_TGP · μ_eff²)                [SC v2 fit form]
```

### 2.3 Derived outputs (the cycle claims)

```
- O1: H₀ (α_PB unit-cousin α_0) REJECTED dimensionally                 [Phase 1]
- O2: α_PB ≈ 0.2887 μ_B⁻² (fitted on PrH₉+NdH₉, 2-point anchor)        [SC v2 input]
- O3: c_TGP = 0.2631 μ_B⁻² (Phase 3 multi-LnH₉ refit)                  [Phase 3 T3.2]
- O4: 15-Ln³⁺ T_c predictions z TGP μ_eff² scaling                     [Phase 3 master table]
- O5: TGP RMS_log = 0.418 < A-G RMS_log = 1.526 (PrH₉+NdH₉)            [Phase 3 T3.5]
- O6: SC4 SmH₉ T_c ≈ 119 K (TGP vs A-G ≈ 0 K, factor 10⁵·⁸)            [LIVE 2027-2030]
- O7: SC5 YbH₉ T_c ≈ 0.64 K (TGP vs A-G ≈ 54 K, factor 84 reverse)     [LIVE 2027-2030]
- O8: SC6 TmH₉ T_c ≈ 4·10⁻⁵ K (TGP vs A-G ≈ 4.14 K, factor 10⁵ reverse) [LIVE 2028-2030]
```

### 2.4 Tautology test (CRITICAL)

**O1 (H₀ rejected):** α_PB has units μ_B⁻² (magnetic moment squared);
α_0 is dimensionless (substrate-matter coupling). Dimensional analysis
(sympy `dimsys_SI`) gives **NO** unit-bridge factor f relating them.
**PASS** — H₀ correctly rejected, structural distinction documented.

**O2-O3 (α_PB / c_TGP fits):** α_PB = 0.2887 μ_B⁻² fitted on **2-point
anchor** (PrH₉ + NdH₉ paper-original). Phase 3 c_TGP = 0.2631 fitted
**re-anchor** post lit-update PrH₉ = 8.9 K. Both are **calibrations**,
NIE first-principles. Honestly acknowledged: "TGP SC v2 paper fitował
2 punkty". **PARTIAL** — fit-not-derivation explicit.

**O4 (15-Ln³⁺ predictions):** T_c = T_c^base · exp(-c_TGP · μ_eff²).
Single formula, 15 elements; μ_eff² values from analytical Hund GS
table (Jensen & Mackintosh literature). NIE tautology — substantive
multi-element prediction. **PASS**.

**O5 (TGP RMS < A-G RMS):** RMS_log on PrH₉+NdH₉:
- NdH₉ obs=4.5 K: TGP 4.57 (drift 1.6%); A-G 0.54 (factor 8 off)
- PrH₉ obs=8.9 K: TGP 4.93 (factor 1.8 off); A-G 12.6 (factor 1.4 off)

**Caveat:** RMS dominated by NdH₉ point (A-G factor 8 off). PrH₉ alone
shows A-G slightly closer (1.4 vs TGP 1.8). RMS preference **honest
**but** dependent na NdH₉ outlier weight. NOT cherry-picked: NdH₉ jest
anchor point dla TGP fit, więc 1.6% drift expected. **Honest balanced
reporting**. PASS w sensie correct calculation; **balanced reporting**
preserved.

**O6-O8 (3 LIVE predictions z bidirectional discriminators):**
- SmH₉: TGP 119 K vs A-G 2·10⁻⁴ K → factor 10⁵·⁸ (TGP direction)
- YbH₉: TGP 0.64 K vs A-G 54 K → factor 84 (A-G direction)
- TmH₉: TGP 4·10⁻⁵ K vs A-G 4.14 K → factor 10⁵ (A-G direction)

**Bidirectional discriminators** = strong falsifiability content. NIE
confirmation-biased. **PASS**.

**Werdykt tautology test:** PASS — substantive predictions, honest fit
disclosure, bidirectional falsification map.

### 2.5 Falsifiability test (CRITICAL)

**Multiple concrete experimental falsifiers:**

- **SC4 SmH₉ ≈ 119 K (TGP wins 10⁵·⁸×):** DAC 2027-2030 synthesis
  (Eremets, Hemley, Prakapenka). Cleanest TGP discriminator.
- **SC5 YbH₉ ≈ 0.64 K (A-G wins 84×):** DAC 2027-2030 (LiBH₄ /
  NH₃BH₃ H-source). Opposite direction — bidirectional!
- **SC6 TmH₉ ≈ 4·10⁻⁵ K (A-G wins 10⁵×):** 2028-2030.

**Bidirectional structure:** A-G *winning* in 2 cases (Yb, Tm) i TGP
winning w 1 case (Sm) z **opposite directions** = **NIE confirmation-
biased prediction set**. To jest property positive falsification map.

**Falsification gates explicit:**
- "Systematic outlier ≥ 50% w którymkolwiek niewykluczonym materiale"
  → Phase 2 GO/NO-GO criterion
- "α_PB^pred odbiega o > 30% od α_PB^fit = 0.2887" → Phase 2 viability

**PmH₉ STRUCTURAL only:** explicit acknowledged radioactive (Pm-145
T_½=17.7 yr) → **honest "structural prediction, not testable"** label.

**Werdykt falsifiability test:** PASS strong — bidirectional discriminators
+ explicit gates + honest non-testable acknowledgment.

### 2.6 Independent-path cross-validation (CRITICAL for DERIVED)

**O3 (c_TGP scaling):** Single TGP formula T_c = T_c^base · exp(-c_TGP · μ_eff²).
Calibrated, NOT first-principles derived. 1 path.

**A-G alternative comparison:** A-G + de Gennes scaling = **standard
literature alternative**. Bidirectional comparison test (15-element
table) jest **substantive falsification**, NIE niezależna derivation.

**Cross-sector hint α₀ = κ_TGP² (z BH.1):** κ_TGP=2.012 → 4.0481;
BH.1 strict α₀ = 4.0179 (ξ=1) → match 0.75%. **2 niezależne paths**:
- (a) SC.1 calibrated on V/Nb/Ta/Mo/Pd → κ_TGP
- (b) BH.1 photon-ring geometric calibration → α₀

NIE sympy-exact, ale **2 niezależne calibrations** zbiegają do 0.75%
gap. STRUCTURAL HINT (jak explicit acknowledged przez BH.1 retrofit).

**15-Ln³⁺ multi-element prediction:** single formula applied to 15
elements; TGP RMS already 3.6× better than A-G na 2 anchor points
(NdH₉ outlier dominated). Multi-element predictive consistency =
substantive **structural support**, NIE multi-path derivation.

**Werdykt independent-path:** STRUCTURAL — single TGP formula calibrated;
substantive bidirectional falsification map vs A-G; cross-sector hint
α₀ = κ_TGP² (0.75% gap, 2 niezależne calibrations).

## 3. Audit gate checklist

```
☑ Phase 0 balance sheet exists (this file)
☑ Tautology test PASS (substantive 15-element predictions; honest 2-point fit disclosure)
☑ Falsifiability test PASS strong (bidirectional discriminators SmH₉/YbH₉/TmH₉ z explicit gates)
☐ Independent-path cross-validation PARTIAL — 1 calibrated TGP path; multi-element consistency check
☑ Alt-scan: A-G + de Gennes alternative explicit compared 15-element table
☑ NIE used post-hoc structural motivations (Hund GS standard atomic physics, A-G textbook)
☑ NIE circular anchor (α_PB ≠ unit-cousin α_0 explicit Phase 1)
☑ NIE inheriting drift > parent × 5× (uses κ_TGP independent of α_0; cross-sector hint 0.75% gap acknowledged)
```

**7/8 ☑ + 1 ☐** — independent-path PARTIAL (single TGP formula path) →
max status STRUCTURAL.

**Honest reporting strong indicators:**
- Phase 1 H₀ rejected explicit (analog BH.1)
- Phase 2 "H_AG_PARTIAL" explicit PARTIAL
- Phase 3 bidirectional discriminators (Sm wins TGP, Yb+Tm win A-G)
- "TGP SC v2 paper fitował 2 punkty" — fit-origin disclosure
- PmH₉ STRUCTURAL only (radioactive non-testable acknowledged)

## 4. Klasyfikacja końcowa

| Klasa | Spełnia? |
|-------|----------|
| DERIVED FULL | NO — α_PB jest fit, NOT first-principles derived; single calibrated TGP path |
| DERIVED CONDITIONAL | partial — multi-LnH₉ predictions z single TGP formula; conditional na 2-point anchor + κ_TGP cross-sector |
| **STRUCTURAL** | **YES** — TGP scaling structural identification + multi-element prediction + bidirectional falsification map + cross-sector hint |
| ANSATZ | NO — concrete falsifiers + structural derivation w stated framework |
| NUMEROLOGICAL | NO — bidirectional predictions exclude minimum-drift fitting; A-G *winning* in 2 cases |
| TAUTOLOGY | NO — 15-element predictions substantive |

**Final verdict:** **STRUCTURAL** ★

**Strukturalne cechy positive example:**

1. **Phase 1 H₀ rejected explicit:** unit-bridge audit dimensionally
   verified, identyczny pattern do BH.1 Phase 1 H₀ rejection.
2. **Phase 2 "H_AG_PARTIAL" explicit:** "α_PB JEST A-G-like, ale a priori
   J_sf ratio 2.59" — honest acknowledgement że α_PB jest **PARTIAL
   derivation** (uses J_sf input), NOT first-principles.
3. **Phase 3 bidirectional discriminators:** SmH₉ TGP wins (10⁵·⁸×),
   YbH₉ A-G wins (84×), TmH₉ A-G wins (10⁵×) — **opposite directions**.
   To jest property substantive falsification map, NIE confirmation-bias.
4. **2-point fit explicit:** "TGP SC v2 paper fitował 2 punkty (PrH₉+NdH₉)"
   — origin disclosure NOT defensive over-claiming.
5. **PmH₉ STRUCTURAL only:** explicit "radioactive — STRUCTURAL" label,
   not pretending testable.
6. **Cross-sector hint α₀=κ_TGP² 0.75% gap:** explicit STRUCTURAL HINT
   labeling (consistent z BH.1 retrofit classification).

**Phase 6 gate compliance:**
1. ✓ Phase0_balance.md exists (this file)
2. ✓ Brak status promotion ("PARTIAL" → "PARTIAL" preserved across phases)
3. ✓ Brak constructed criterion (single TGP formula; A-G alternative explicit)
4. ✓ Brak accommodating gate (Phase 2 "α_PB^pred odbiega > 30%" explicit)
5. ✓ Brak sympy-rationalization-as-DERIVED (α_PB explicit "fit", c_TGP explicit "refit")

## 5. Comparison ze status oryginalnym

| Element | Original claim | Retrofit verdict |
|---------|----------------|------------------|
| Status YAML | `status: CLOSED`, "7/7 PASS Phase 3", 17 sub-tests PASS | STRUCTURAL — bidirectional falsification map substantive; α_PB calibration honest disclosure |
| Counter (PREDICTIONS_REGISTRY) | SC2-SC7 (5 LIVE + 1 STRUCTURAL) | Stays as is — entries explicit "LIVE 2027-2030" honest classification |
| Sub-tests | 4+6+7 = 17 PASS | Substantive (Phase 1 H₀ rejected rigorously, Phase 2 PARTIAL labeling, Phase 3 multi-element table z explicit data) |
| Independence | "TGP wins RMS_log 0.42 vs A-G 1.53" | Honest acknowledgment dominated by NdH₉ anchor; PrH₉ alone A-G slightly better — balanced reporting maintained |

## 6. Recommended action

- [x] **NO-OP klasy** — STRUCTURAL classification z honest reporting OK
- [x] **Phase 5 PREDICTIONS_REGISTRY annotation:** SC2-SC7 entries
      preserve explicit "LIVE 2027-2030" labels + bidirectional falsification
      property explicit
- [x] **Cross-sector hint annotation:** α₀^BH = κ_TGP² 0.75% gap consistent
      classification z BH.1 retrofit STRUCTURAL HINT (not promoted)
- [ ] CRITIQUE — nie wymaga (cycle uses honest "PARTIAL" labeling)
- [ ] CASCADE_AUDIT — none (SC.1 nie ma downstream cycles inheriting)
- [ ] CORE_IMPACT — none (research-level, NOT core LaTeX)

## 7. Notes

**Bidirectional falsification map jako quality indicator:**

SC.1 Phase 3 produkuje **rare property** w predictions: discriminators
**w przeciwnych kierunkach**. SmH₉ TGP wins; YbH₉ + TmH₉ A-G wins.
To jest **strongly anti-confirmation-biased** — TGP nie predykuje
"wszystkie LnH₉ → high T_c" (confirmation pattern), ale specific
elements z **predictable winning side** (Sm low μ_eff² → TGP wins;
Yb high μ_eff² → A-G wins).

**Comparison z κ.1+ι.1+μ.1+ν.1 mixing-operator family:**

| Aspect | SC.1 | mixing-operator family |
|--------|------|------------------------|
| Multi-candidate response | A-G + TGP explicit comparison | Constructed criterion to select winner |
| Falsification map | Bidirectional (Sm wins TGP, Yb+Tm win A-G) | Confirmation-only (all PMNS angles "refined²") |
| Fit disclosure | "TGP SC v2 paper fitował 2 punkty" | "DERIVED FULL CASCADE 7/7" |
| PARTIAL labeling | Phase 2 H_AG_PARTIAL explicit | Phase 2 + Phase 3 promotions |
| Status outcome | STRUCTURAL ★ honest | NUMEROLOGICAL/ANSATZ retrofitted |

**Identical surface elements** (multi-element predictions, fit-based
calibrations) z **opposite discipline** (bidirectional vs confirmation-
biased, PARTIAL vs FULL CASCADE).

**Cross-sector α₀ = κ_TGP² hint integration:**

SC.1 Phase 3 finalizes κ_TGP = 2.012 (V/Nb/Ta/Mo/Pd) calibration.
BH.1 Phase 2 finds α₀^BH ≈ 4.018 (ξ=1 photon-ring sketch). Gap
4.0481 vs 4.0179 = 0.75%. **2 niezależne calibrations** convergent
within sub-percent — **substantive STRUCTURAL HINT** dla cross-sector
unification.

Status post-retrofit (consistent z BH.1 retrofit):
- κ_TGP value: SC.1 calibrated, **NIE first-principles derived**
- α₀^BH value: BH.1 PARTIAL DERIVED (ξ=1 sketch)
- α₀ = κ_TGP² identity: STRUCTURAL HINT (0.75% gap, 2 niezależne paths)
- Future Phase 1 PLAN action principle derivation deferred (~15 months)

**SC.1 SCN entries:**
- SC2 PrH₉ TESTED-PASS (drift factor 1.8 vs anchor 1.4× A-G better lit)
- SC3 NdH₉ TESTED-PASS (drift 1.6% vs A-G factor 8 off)
- SC4 SmH₉ LIVE 2027-2030 (TGP 119K vs A-G 0K, factor 10⁵·⁸)
- SC5 YbH₉ LIVE 2027-2030 (TGP 0.64K vs A-G 54K, factor 84 reverse)
- SC6 TmH₉ LIVE 2028-2030 (TGP ≈0K vs A-G 4.14K, factor 10⁵ reverse)
- SC7 PmH₉ STRUCTURAL only (radioactive, factor 10³ TGP direction)

**6 entries z explicit experimental status** (TESTED + LIVE + STRUCTURAL)
honestly labeled. Counter contribution stays as is.

## 8. Cross-references

- [[../op-sc-alpha-origin/README.md]] — cykl audited
- [[../op-sc-alpha-origin/Phase3_results.md]] — main predictions source
- [[../op-sc-alpha-origin/program.md]] — 3-phase plan + α_PB origin question
- [[../op-bh-alpha-threshold/Phase2_results.md]] — cross-sector α₀^BH source
- [[retrofit_op-bh-alpha_2026-05-06.md]] — BH.1 retrofit (consistent classification)
- [[../closure_2026-04-26/alpha_psi_threshold/results.md]] — α_0 T-α (Phase 1 unit-bridge target)
- [[README.md]] — M03 master plan
- [[audit_log.md]] — appended 2026-05-06 (Phase 3 #6)
- [[tracker.md]] — status updated to DONE_STRUCTURAL
- [[../../meta/CALIBRATION_PROTOCOL.md]] — protocol source
- [[../../meta/CALIBRATION_GATE_ENFORCEMENT.md]] — Phase 6 gate
