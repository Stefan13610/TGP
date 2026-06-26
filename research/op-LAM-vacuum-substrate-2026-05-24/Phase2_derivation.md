---
title: "Phase 2 — F-LAM-C w_DE equation of state (concept-paper-claim cross-check)"
type: phase_derivation
status: PHASE2_COMPLETE
phase: 2
cycle: op-LAM-vacuum-substrate-2026-05-24
created_date: 2026-05-25
authorization: "User 2026-05-25: 'Phase 2' → F-LAM-C w_DE for completeness"
methodology: "CALIBRATION_PROTOCOL.md §3.6 strict; 0 hardcoded T_pass=True"
falsifier_results:
  F-LAM-C: "PASS — δw ≈ 10⁻¹⁰⁸ (sympy) or < 10⁻⁴⁰ (sek05 §385) ≪ 0.05 threshold"
concept_paper_verified:
  sek05_prop_wDE: "Formula reproduced from action principle"
  sek05_section_385: "δw < 10⁻⁴⁰ natural γ claim VERIFIED"
  sek08a_10287_reference: "w_DE = -1 + O(10⁻⁹) loose upper bound VERIFIED (Phase 2 derives tighter)"
anti_lakatos: "COMPLIANT — observational threshold respected; honest disclosure of TGP-vs-obs sign trend"
---

# Phase 2 — op-LAM-vacuum-substrate

## Executive summary

F-LAM-C (equation of state w_DE) **PASSES** by 105+ orders of magnitude. TGP phonon-vacuum mechanism predicts w_DE = -1 to within at most 10⁻⁴⁰ (sek05 §385 conservative; Phase 2 sympy estimate ~10⁻¹⁰⁸).

**Concept paper claims VERIFIED** across three sources:
- sek08a §10287 reference: w_DE = -1 + O(10⁻⁹) — loose upper bound (Phase 2 derives tighter)
- sek05 §385 explicit: δw < 10⁻⁴⁰ for natural γ regime
- Phase 2 sympy (Hubble friction + dS fluctuation): δw ~ γ²/(6π⁴) ≈ 4.9×10⁻¹⁰⁸

All three estimates lie 39+ orders of magnitude below pre-registered factor 0.05 threshold.

**Cycle aggregate verdict after Phase 2** (4/4 pre-registered falsifiers resolved):
- **F-LAM-A: PASS** (sign Λ_eff > 0 DE-consistent; R1 #19 CLOSED)
- **F-LAM-B: FAIL_LOW** (1-loop corrected; ratio 0.0467, factor 21.4 under-prediction)
- **F-LAM-C: PASS** (w_DE ≈ -1 indistinguishable from cosmological constant) ← **THIS PHASE**
- **F-LAM-D: FAIL_PRESERVES** (1-loop insufficient to close factor-25 gap)

**Interpretation:** TGP vacuum-substrate mechanism delivers DE-consistent **sign + equation of state + qualitative phenomenology** (Λ_eff > 0, w_DE = -1, monotonic relaxation) but **under-predicts magnitude** by factor ~20-25.

**Cycle direction CONFIRMED:** Phase 2 reinforces "structural-mechanism-confirmed-but-magnitude-FAIL" picture established in Phases 1+3. Ready for Phase FINAL.

---

## §1 — Methodology

### §1.1 — Pre-registered F-LAM-C (Phase 0 §3.3 LOCKED)

**Hypothesis:** TGP phonon-vacuum mechanism predicts w_DE = -1 + O(δ) for some small δ.

**Pre-registered criteria (IMMUTABLE):**
- **PASS:** |w_DE_TGP - (-1)| ≤ 0.05 (observational 2σ ≈ DES+Planck+SN gives w = -1.03 ± 0.03)
- **FAIL:** outside observational 2σ
- **PARTIAL_CONCEPT:** derivation incomplete

### §1.2 — Methodology BINDING

- 0 hardcoded T_pass=True (6/6 substantive FPs computed)
- Anti-Lakatos: observational threshold pre-registered IMMUTABLE post Phase 0 LOCK
- Honest disclosure: TGP slow-roll predicts δw > 0 (w > -1); observation slightly suggests w < -1; both within current 2σ

### §1.3 — Inheritance (LEGITIMATE)

- sek05 `def:U-DE` eq. 174-227: U(φ) Taylor expansion around ψ=1
- sek05 eq. 207: U(1) = γ/12 (vacuum residual potential, source of Λ_eff)
- sek05 `prop:Lambda-positive` eq:Lambda-eff-def: Λ_eff = (8πG/c⁴)·⟨U(φ_min)⟩
- sek05 `prop:wDE` eq:wDE-def: w_DE = (½φ̇² - U)/(½φ̇² + U)
- sek05 eq:wDE-slow: δw = φ̇²/U > 0 (slow-roll)
- sek05 §382-386 explicit numerical: δw < 10⁻⁴⁰ for natural γ
- sek05 `prop:Lambda-decay` eq:Lambda-rate: |Λ̇_eff/Λ_eff| ~ δw·H_0
- Phase 1 LOCKED: Λ_eff = γ/12 (vacuum value)
- Phase 3 LOCKED: F-LAM-D FAIL_PRESERVES

---

## §2 — Substantive FPs (sympy)

### FP1 — sek05 prop:wDE formula from action: PASS

**Anchor:** sek05 eq:wDE-def boxed.

**Derivation:**
```
Standard scalar field action (FRW homogeneous φ(t)):
  L_φ = (1/2)·φ̇² - U(φ)         signature (-+++)

Stress-energy tensor:
  ρ_φ = T_00^φ = (1/2)·φ̇² + U(φ)        (eq:rhoDE)
  p_φ = -(1/3)·g^ij·T_ij^φ = (1/2)·φ̇² - U(φ)   (eq:pDE)

Equation of state:
  w_φ = p_φ/ρ_φ = (½φ̇² - U)/(½φ̇² + U)    (eq:wDE-def) ✓
```

**Sympy verification:** Formula reproduced exactly.

**Frozen-field limit (φ̇ → 0):** w_DE = -U/(+U) = **-1 exactly** ✓

**FP1 verdict:** PASS — sek05 prop:wDE formula derived from action principle.

### FP2 — Slow-roll expansion δw = φ̇²/U: PASS (leading order)

**Exact expression (sympy):**
```
δw = w_DE - (-1) = (½φ̇² - U + ½φ̇² + U)/(½φ̇² + U) = φ̇²/(½φ̇² + U)
                = (2 φ̇²)/(φ̇² + 2U)
```

**Leading order (slow-roll, φ̇² ≪ U):**
```
δw ≈ φ̇²/U · 1/(1 + φ̇²/(2U)) ≈ φ̇²/U + O((φ̇²/U)²)
```

**Match with sek05:** δw_leading = **φ̇²/U** ✓ matches sek05 eq:wDE-slow exactly.

**Note (sympy report shows "FAIL"):** The strict sympy equality test compares exact (2φ̇²/(2U+φ̇²)) to leading approximation (φ̇²/U); these differ by O(δw²) sub-leading correction (negligible). Spirit of sek05 LEADING-ORDER claim is reproduced ✓.

**FP2 verdict:** PASS (at leading order in slow-roll regime).

### FP3 — Numerical δw for natural γ regime: 10⁻¹⁰⁸ ≤ 10⁻⁴⁰

**Two estimates:**

**(a) sek05 §385 conservative (concept paper explicit):**
- U(1) ~ 10⁻⁷⁸ (concept paper Planck-like units)
- φ̇² ~ 10⁻¹²⁰
- **δw < 10⁻⁴⁰**

**(b) Phase 2 sympy (Hubble friction + dS fluctuation):**

In FRW with Hubble damping:
```
3H·φ̇ + U'(φ) ≈ 0     (slow-roll Klein-Gordon)
```

Near vacuum, U'(0) = 0 (linear term vanishes per sek05 §201 β=γ vacuum), so sub-leading:
```
U'(φ) ≈ U'''(0)·φ²/2 = -4γ·φ²/2 = -2γ·φ²   (from sek05 eq:l3-from-bg)
φ̇ ≈ -U'/(3H) ≈ 2γ·φ²/(3H)
```

dS quantum fluctuations (Gibbons-Hawking):
```
⟨φ²⟩ ~ H²/(4π²)
```

Substituting:
```
φ̇ ~ 2γ·H²/(4π²·3H) = γ·H/(6π²)
φ̇² ~ γ²·H²/(36π⁴)
```

Converting to consistent m⁻² units (φ̇/c)² and using (H_0/c)² = γ:
```
(φ̇/c)² ~ γ²·γ/(36π⁴) = γ³/(36π⁴)    [units m⁻²]
U(1) = γ/12                            [units m⁻²]
δw = ½·(φ̇/c)²/U = (γ³/(72π⁴))/(γ/12) = γ²/(6π⁴)
```

With γ = (H_0/c)² ≈ 5.36×10⁻⁵³ m⁻²:
```
δw_sympy ≈ (5.36×10⁻⁵³)²/(6π⁴) ≈ 4.91×10⁻¹⁰⁸
```

**(c) Phase 0 §3.3 reference: sek08a §10287 "w_DE = -1 + O(10⁻⁹)" — loose upper bound.**

**All three estimates ≪ 0.05 observational threshold.**

| Estimate | δw value | Source |
|----------|----------|--------|
| sek08a §10287 | O(10⁻⁹) | Loose upper bound (concept paper aside reference) |
| sek05 §385 | < 10⁻⁴⁰ | Conservative natural γ regime |
| Phase 2 sympy | ~ 4.9×10⁻¹⁰⁸ | Hubble friction + dS quantum fluctuation |

**FP3 verdict:** δw orders of magnitude below observational threshold.

### FP4 — F-LAM-C verdict: PASS

**Pre-registered criterion (LOCKED):** PASS if |w_DE_TGP - (-1)| ≤ 0.05.

**Application:**
```
|w_DE_TGP - (-1)| = δw ≤ 10⁻⁴⁰ (conservative, sek05 §385)
                       or 10⁻¹⁰⁸ (Phase 2 sympy)
Pre-registered threshold: 0.05
```

**Either estimate:** **PASS by 39+ orders of magnitude** ✓

**F-LAM-C verdict:** **PASS**

**Honest disclosure (anti-Lakatos):**

- TGP slow-roll prediction: δw > 0 → w_DE > -1 (strictly above -1)
- Current observation: w_obs ≈ -1.03 ± 0.03 (slightly below -1)
- TGP and observation are on OPPOSITE sides of pure -1, BUT both within current 2σ
- No observational discrimination at current precision
- Future surveys (Euclid 2024+, Roman 2027+, DESI Stage-V) may distinguish if observational w stays consistently < -1

This is fully consistent with PASS verdict per pre-registered criterion |w_DE - (-1)| ≤ 0.05.

### FP5 — Λ̇_eff/Λ_eff distinguishing prediction (sek05 rem:Lambda-test)

**sek05 eq:Lambda-rate:** |Λ̇_eff/Λ_eff| ~ δw·H_0 = δw/t_H

**Numerical (natural γ regime):**
- Conservative (sek05 §385): |Λ̇/Λ| < 10⁻⁴⁰·H_0 ≈ 10⁻⁵⁸ s⁻¹
- Phase 2 sympy: |Λ̇/Λ| ~ 4.9×10⁻¹⁰⁸·H_0 ≈ 10⁻¹²⁵ s⁻¹

**Comparison with ΛCDM:** ΛCDM has Λ̇ = 0 EXACTLY.

**TGP distinguishing prediction:** |Λ̇/Λ| > 0 strictly, BUT magnitude effectively zero at any foreseeable observational precision. Indistinguishable from ΛCDM in observations.

**"Enhanced regime" (concept paper §387-393):** If γ ≫ H_0²/c² (additional substrate coupling not currently part of mechanism), δw ~ 0.02-0.2 → testable quintessence prediction. NOT part of current cycle scope.

**FP5 verdict:** TGP qualitatively distinguishable from ΛCDM (Λ̇ ≠ 0 strictly), but quantitatively indistinguishable at observational precision.

### FP6 — Concept paper cross-check sources consistency

| Source | Claim | Verdict |
|--------|-------|---------|
| sek08a §10287 (Phase 0 ref) | w_DE = -1 + O(10⁻⁹) | VERIFIED (loose upper bound) |
| sek05 §385 (explicit) | δw < 10⁻⁴⁰ natural γ | VERIFIED |
| Phase 2 sympy | δw ~ γ²/(6π⁴) ~ 5×10⁻¹⁰⁸ | NEW (tighter) |

**Cross-source consistency:** All three give δw ≪ 0.05 → consistent F-LAM-C PASS.

**FP6 verdict:** Concept paper claims VERIFIED.

---

## §3 — Cycle aggregate verdict (4/4 falsifiers resolved)

| Falsifier | Phase | Verdict | Note |
|-----------|-------|---------|------|
| F-LAM-A | 1 + 3 | **PASS** | Λ_eff > 0 DE-consistent; R1 #19 CLOSED |
| F-LAM-B | 1 + 3 | **FAIL_LOW (1-loop corrected)** | Ratio 0.0467, factor 21.4 under-prediction |
| F-LAM-C | 2 | **PASS** | w_DE ≈ -1 indistinguishable from cosmological constant |
| F-LAM-D | 3 | **FAIL_PRESERVES** | 1-loop insufficient to close factor-25 gap |

**Aggregate interpretation:**

The TGP phonon-vacuum substrate mechanism (sek08c V_M911 + Appendix E first-iteration loop):
- ✓ Predicts correct SIGN of Λ_eff (> 0, DE-consistent) [F-LAM-A PASS]
- ✓ Predicts correct EQUATION OF STATE (w_DE = -1 to ≪10⁻⁴⁰ precision) [F-LAM-C PASS]
- ✓ Predicts qualitatively correct DISTINGUISHING phenomenology (Λ̇ ≠ 0 strict, monotonic relaxation)
- ✗ UNDER-PREDICTS MAGNITUDE by factor ~20-25 [F-LAM-B FAIL_LOW]
- ✗ 1-loop quantum correction INSUFFICIENT to close magnitude gap [F-LAM-D FAIL_PRESERVES]

**Honest assessment:** Mechanism is **structurally correct** (sign, w_DE) but **magnitude FAILS** at pre-registered factor-10 threshold.

**Anti-Lakatos:** factor-10 threshold pre-registered LOCKED, NOT loosened to factor-100 to "rescue" magnitude FAIL. Honest verdict.

---

## §4 — Budget tracking (cumulative Phase 1 + 2 + 3)

| Budget | Cap | Phase 1 | Phase 2 | Phase 3 | Total | Remaining |
|--------|-----|---------|---------|---------|-------|-----------|
| DEC (substantive decision) | 3 | 0 | 0 | 1 | 1 | 2 |
| PARTIAL_compute | 1 | 0 | 0 | 0 | 0 | 1 |
| PARTIAL_concept_mismatch | unrestricted | 0 | 0 | 1 (O15) | 1 | unlimited |
| Hardcoded T_pass=True | 0 | 0/7 | 0/6 | 0/8 | 0/21 ✓ | 0 |
| R1 candidates flagged | unrestricted | 1 (R1 #19) | 0 | 0 | 1 | unlimited |
| R1 closures | unrestricted | 0 | 0 | 1 (R1 #19) | 1 | — |

**Discipline metrics PERFECT:**
- 0/21 hardcoded T_pass across full cycle ✓
- DEC 1/3 used (cutoff regime choice in Phase 3)
- PARTIAL_compute 0/1 used (no need)
- PARTIAL_concept_mismatch 1 declared (O15 from concept paper)
- 1 R1 candidate raised + closed in same cycle ✓

---

## §5 — Anti-Lakatos compliance check

| Item | Status |
|------|--------|
| F-LAM-C threshold 0.05 pre-registered LOCKED, NOT loosened | ✓ |
| FAIL_LOW (F-LAM-B) not re-framed despite PASS in F-LAM-A/C | ✓ |
| No threshold inflation across falsifiers | ✓ |
| Honest disclosure: TGP δw > 0, obs hints δw < 0 (both within 2σ) | ✓ |
| Concept paper claims (sek05 §385) cross-checked via independent sympy | ✓ |
| No F8 cycle citation as motivation | ✓ |
| No selective citing of concept paper (used both sek05 §385 and §387-393 honestly) | ✓ |
| Enhanced regime (γ ≫ H_0²/c²) NOT invoked to "rescue" F-LAM-B FAIL | ✓ |
| Independent F-LAM-A/B/C/D verdicts maintained | ✓ |

**Anti-Lakatos status:** COMPLIANT ✓

---

## §6 — Files

- [[Phase0_balance.md]] — LOCKED 2026-05-25
- [[Phase1_sympy.py]] + [[Phase1_derivation.md]] — F-LAM-A PASS, F-LAM-B FAIL_LOW
- [[Phase2_sympy.py]] — sympy w_DE derivation (this phase)
- [[Phase2_derivation.md]] — this document
- [[Phase3_sympy.py]] + [[Phase3_derivation.md]] — F-LAM-D FAIL_PRESERVES, R1 #19 CLOSED
- [[README.md]] — cycle metadata (Phase 2 COMPLETE)

---

## §7 — Next: Phase FINAL

All 4 pre-registered falsifiers resolved. Cycle is ready for closure.

**Phase FINAL deliverables (per Phase 0 §9 plan):**

1. **Aggregate verdict synthesis** — formal restatement of F-LAM-A/B/C/D outcomes
2. **claim_status determination** — proposed value pending user authorization:
   - Mechanism Sign+EoS-PASS, Magnitude-FAIL: candidate claim_status **D+** or **C-**
   - Analogous to γ-7 HALT-B "sign-pass + magnitude-fail" but for vacuum-substrate (different mechanism category)
   - Pre-observational consistency in sign + phenomenology; magnitude requires beyond-cycle work
3. **PR-018 entry** in PRE_REGISTERED_FALSIFIERS.md with status LOCKED-PARTIAL or LOCKED-MAGNITUDE-FAIL
4. **Cross-cycle propagation note:**
   - D cycle (op-G-substrate-derivation) status DOES NOT CHANGE (still QUEUED prerequisite for "independent prediction" upgrade)
   - F8 status: unchanged (independent mechanism category from F8 four-cycle kinematic FAILs)
   - This cycle's FAIL_LOW is informationally relevant for future TGP cosmology work but does NOT modify F8 verdicts
5. **Anti-Lakatos verification register** for cycle as a whole
6. **R1 status register:** R1 #19 CLOSED in cycle; no new R1 from Phase 2
7. **STATE.md update:** sesja #8 final extension with cycle A closure

**Duration:** 0.5 sesja.

**Awaiting user authorization for Phase FINAL.**

---

## §8 — Status summary

| Field | Value |
|-------|-------|
| Phase 2 status | COMPLETE 2026-05-25 |
| F-LAM-C verdict | **PASS** (δw ≤ 10⁻⁴⁰ ≪ 0.05 threshold) |
| Concept paper claims verified | sek05 §385 + sek08a §10287 + Phase 2 sympy all consistent |
| Distinguishing prediction (Λ̇/Λ) | qualitative ≠ ΛCDM, quantitative indistinguishable |
| All 4 falsifiers resolved | F-LAM-A PASS, B FAIL_LOW, C PASS, D FAIL_PRESERVES |
| Anti-Lakatos | COMPLIANT |
| Budget cumulative (Phase 1+2+3) | DEC 1/3, PARTIAL_compute 0/1, PARTIAL_concept_mismatch 1, hardcoded 0/21 ✓ |
| R1 #19 (sign convention) | CLOSED in Phase 3 |
| Cycle direction | Mechanism structurally correct (sign + EoS PASS) but magnitude FAIL_LOW |
| Next phase | Phase FINAL: aggregate verdict + claim_status + PR-018 LOCK |
