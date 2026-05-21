---
title: "Phase 1 — Results: empirical test universal mass formula on quark sector → HALT-B (Path C structural insufficiency)"
date: 2026-05-16
parent: "[[./README.md]]"
type: phase-results
phase: 1
sympy_pass: "6/13"
fp_count: 10
lit_count: 2
declarative_separate: 1
hardcoded: 0
substance_metrics: "6/13 PASS (T1-T4 + T12 + T13); T5-T9 ALL FAIL (5/5 quark ratios drift ~100%); T10 FAIL; T11 FAIL (structural ceiling 2.68× vs required 80,000×)"
status: 🔴 EMPIRICAL TEST DONE — falsifier triggers HALT-B (Path C); structural insufficiency confirmed
---

# Phase 1 — Results: quark sector mass formula empirical test

## §0 — Summary

```
████████████████████████████████████████████████████████████████████
█                                                                  █
█  Phase 1 empirical test: 6/13 sympy PASS                         █
█                                                                  █
█  T1-T4 PASS: formula structure + lepton calibration              █
█  T5-T9 FAIL: 0/5 quark ratios within 10% tolerance               █
█  T10 FAIL: 0/6 required g_0_q in audit range [0.817, 0.891]      █
█  T11 FAIL: structural ceiling 2.68× vs required 80,000×          █
█  T12 PASS: sek08b:529 audit range consistency                    █
█  T13 PASS: S05 preserved                                         █
█                                                                  █
█  PRE-REGISTERED FALSIFIER VERDICT: HALT-B (Path C)               █
█                                                                  █
████████████████████████████████████████████████████████████████████
```

**Decision (per pre-registered falsifier §0.2 README):**

> Quark ratios passed: **0/5** < 3/5 threshold.
> Universal-Φ-kink description warstwy 3c **INSUFFICIENT dla quark sektora**.
> Strukturalna konieczność extension (per-family g_0 poza audit / quark-specific
> kink topology / multi-family substrate).

**Status:** HALT-HONEST B; cycle closes z explicit structural insufficiency documentation.

## §1 — Per-test results

### T1 — Universal formula structure (FIRST_PRINCIPLES) — **PASS**

```
e²·(1 − α/4)|_{α=2} = e²·(1 − 1/2) = e²/2 ≈ 3.69453
```

Inheritance z L08-e² β(α=2)=e²/2 canonical (LIVE) and why_n3 Phase 5 universal formula
verified symbolically. **Formula structure correctly identified.**

### T2 — Mass ratio formula derivation (FIRST_PRINCIPLES) — **PASS**

```
m_i/m_j = (c_M · A_i² · g_i^(e²/2)) / (c_M · A_j² · g_j^(e²/2))
        = (A_i/A_j)² · (g_i/g_j)^(e²/2)
```

`c_M` cancels w ratios; ratio formula correctly derived. **Provides direct test
mechanism free from anchor.**

### T3 — Lepton anchor verification (LITERATURE_ANCHORED) — **PASS**

| Ratio | Predicted | PDG | Drift |
|---|---|---|---|
| m_μ/m_e | 206.7556 | 206.7683 | **0.0061%** ✓ |
| m_τ/m_e | 3073.50 | 3477.23 | 11.61% (Koide-via-A^(5-α) outlier) |

m_μ/m_e match z PDG <0.01% — confirms formula correctly reproduces lepton sektor (sanity).
m_τ/m_e drift 11.61% to **known τ outlier** documented w why_n3 PHASE4_5 §4.1 (formula
F1 vs F2 reconciliation issue; preserved B+ status from L08-e² cycle). Sanity preserved.

### T4 — A_tail(g_0) power-law fit (FIRST_PRINCIPLES) — **PASS**

```
log(A_tail) = 3.8404 · log(g_0) − 1.6866
A_tail(g_0) ≈ 0.1851 · g_0^3.8404
R² = 0.998793   (3-point fit on lepton calibration)
residuals: [0.0171, -0.0541, 0.0371]  (in log-space)
```

Power-law fit z 3 lepton points (e, μ, τ) z **R² > 0.99**. **Provides A_tail(g_0)
extrapolation function for quark range.** Note: residuals show ~3-5% deviation in
log-space (~exp(0.05) ≈ 5% in linear A_tail), acceptable for structural ceiling test.

### T5 — m_c/m_u predicted vs PDG ≈ 588 (FIRST_PRINCIPLES) — **FAIL**

```
m_c/m_u predicted (with g_0_c, g_0_u ∈ [0.817, 0.891]) = 1.000
PDG m_c/m_u = 587.963
Drift = 99.83%  (>> 10% tolerance)
```

**Why 1.000?** With audit range [0.817, 0.891], max achievable m/m_e = 1.279 (at
g_0=0.891) and min = 0.477 (at g_0=0.817). Both m_c (≈ 2486 m_e) AND m_u (≈ 4.23 m_e)
are OUTSIDE this range — both pinned to upper endpoint → ratio = 1.0.

### T6 — m_b/m_d predicted vs PDG ≈ 895 — **FAIL**

```
Predicted = 1.000 vs PDG = 895.075; drift = 99.89%
```

Same structural reason: both targets above audit envelope → both pinned to max endpoint.

### T7 — m_t/m_c predicted vs PDG ≈ 136 — **FAIL**

```
Predicted = 1.000 vs PDG = 135.976; drift = 99.26%
```

### T8 — m_s/m_d predicted vs PDG ≈ 20 — **FAIL**

```
Predicted = 1.000 vs PDG = 20.000; drift = 95.00%
```

Even the SMALLEST PDG quark ratio (m_s/m_d = 20) is unreachable — drift still 95%.
**The narrowest ratio the formula could in principle test (smallest PDG value) is
still WAY outside structural ceiling.**

### T9 — m_b/m_t predicted vs PDG ≈ 0.024 — **FAIL**

```
Predicted = 1.0000 vs PDG = 0.0242; drift = 4031.3%
```

Inverse ratio also fails (m_b lighter than m_t, ratio < 1 expected but formula
gives 1.0).

### T10 — Required g_0_q values in audit range (FIRST_PRINCIPLES) — **FAIL**

Free-fit g_0_q values (no audit constraint, just match PDG mass):

| Quark | Required g_0 | In [0.817, 0.891]? |
|---|---|---|
| u | 0.9898 | OUT (slightly above) |
| d | 1.0592 | OUT |
| s | 1.3784 | OUT (near μ value 1.407) |
| c | 1.2501 | OUT |
| b | 1.2500 | OUT |
| t | 1.2500 | OUT |

**0/6 in audit range.** Required g_0 values for c, b, t saturate near 1.25
(numerical solver hits formula extrema in this regime — note: power-law A_tail
fit is monotonic but mass formula has ∂m/∂g_0 = 0 around lepton τ region;
nonlinear behavior z polynomial structure of mass formula in this regime).

For light quarks (u, d, s) g_0 ranges 0.99 to 1.38 — all ABOVE audit max 0.891.
For heavy quarks the equation becomes degenerate (multiple roots), but no
solution in audit range either.

**Required g_0 values DO NOT fit audit range.** Universal ODE claim sek08b:529
empirically contradicted.

### T11 — Structural ceiling test (FIRST_PRINCIPLES) — **FAIL**

```
Max achievable mass ratio in g_0 ∈ [0.817, 0.891]:
  m(g_0=0.891) / m(g_0=0.817) = 1.279 / 0.477 ≈ 2.681

Required extreme quark ratio: m_t/m_u ≈ 79,949

Achievable / Required = 3.35 × 10⁻⁵   (5 orders of magnitude short)
```

**This is the decisive structural test.** With g_0 confined to audit range
[0.817, 0.891] (1.091× span), universal formula's max achievable mass ratio is
**2.68×**. Required quark mass ratios start at 20 (m_s/m_d) and reach 79,949
(m_t/m_u). **Universal formula cannot span 5 orders of magnitude when g_0 is
confined to 0.09× audit range.**

This is NOT a numerical or calibration artifact — it's a **structural ceiling**.
No reasonable A_tail(g_0) calibration could rescue formula within audit constraint.

### T12 — sek08b:529 audit range consistency (LITERATURE_ANCHORED) — **PASS**

Audit range used [0.817, 0.891] **matches verbatim** sek08b:529:

> "Universalność kwarkowa: ten sam ODE działa na leptony i kwarki
>  (g_0 ∈ [0,817; 0,891])"

Literature consistency confirmed. **No range misinterpretation.**

### T13 — S05 single-Φ preservation (DECLARATIVE) — **PASS**

Universal formula uses **single Φ field**; quarks are kinki Φ z różnymi g_0 values.
**S05 axiom preserved** — cycle test does NOT introduce multi-field substrate.

## §2 — Substance metrics

| Metric | Value |
|---|---|
| Total sympy | 13 sub-tests |
| PASS | 6/13 |
| FIRST_PRINCIPLES | 3/10 PASS (T1, T2, T4) — but 10 of 13 tests = 76.9% FP ratio |
| LITERATURE_ANCHORED | 2/2 PASS (T3 inherited, T12 verified) |
| DECLARATIVE | 1/1 PASS (T13 separate) |
| Hardcoded `T_pass = True` | **0** (Phase 6 ABSOLUTE BINDING) |
| Substance ratio | 10 FP (76.9%) + 2 LIT (15.4%) + 1 DEC (7.7%) |

**Note on FP PASS:** 3 of 10 FP tests PASS — the LOW FP PASS rate reflects **empirical
test outcome (formula structurally fails)**, NOT substance defect of test design. Each
test asks a substantive physics question; FAIL is the honest answer for quark sektor
given audit constraint.

## §3 — Pre-registered falsifier rule check

**Rule (§0.2 README, BINDING, pre-registered 2026-05-16):**

> Jeśli formula z g_0_q ∈ [0.817, 0.891] NIE reprodukuje **≥3 z 5 niezależnych
> quark mass ratios** w tolerancji **10%**, universal-Φ-kink description warstwy
> 3c jest **INSUFFICIENT dla quark sektora** → strukturalna konieczność extension.

**Observed (Phase 1 T5-T9):**
- m_c/m_u: drift 99.8% (FAIL)
- m_b/m_d: drift 99.9% (FAIL)
- m_t/m_c: drift 99.3% (FAIL)
- m_s/m_d: drift 95.0% (FAIL)
- m_b/m_t: drift 4031% (FAIL)

**Score: 0/5 within 10%.** Falsifier rule triggers UNAMBIGUOUSLY.

**Verdict (BINDING):** HALT-B (Path C).

## §4 — Recovery scope check (anti-Lakatos)

Pre-registered recovery scope (§0.2 README, LOCKED):

| Allowed direction | Status | Outcome |
|---|---|---|
| Per-family g_0 variation w obrębie audit range | T10 (FAIL) | Required g_0 values OUTSIDE audit range for ALL 6 quarks |
| Subleading A_tail corrections | T4 (PASS R² 0.999) | A_tail fit good; not the problem |
| Wilson coef extensions z L08-Dirac-Wilson | LIVE | Bounded ~10⁻¹⁶ lab-scale; orders of magnitude too small |

**Recovery scope EXHAUSTED.** All allowed directions checked:
- (1) g_0 variation w audit range: max ratio 2.68× << required 20×–80,000×
- (2) A_tail recalibration: would require A_tail span 240× (vs lepton 15×) — outside ODE structure for narrow g_0 range
- (3) Wilson coefs: lab-scale ~10⁻¹⁶ negligible vs orders-of-magnitude mass drift

| Forbidden direction | Status |
|---|---|
| Post-hoc tuning g_0_q poza audit | NOT ATTEMPTED (anti-Lakatos) |
| Multi-field substrate (S05) | NOT ATTEMPTED |
| Z₂ → SU(N) extension | NOT ATTEMPTED |

**Recovery exhausted → H1c triggered (per §0.2 if_recovery_exhausted):**
1. Universal-Φ-kink INSUFFICIENT dla quark sektora ✓
2. Multi-session extension cycle required (CKM coupling derivation / quark-specific kink topology)
3. OR Path C audit (quark sektor permanent (H) w warstwie 3c)

## §5 — Per-test result table

| Test | Klasa | Status | Drift / Result |
|---|---|---|---|
| T1 | FP | ✅ PASS | e²/2 exact symbolic |
| T2 | FP | ✅ PASS | ratio formula algebraic |
| T3 | LIT | ✅ PASS | μ/e 0.006%; τ/e 11.6% (known τ outlier) |
| T4 | FP | ✅ PASS | R² 0.999 power-law fit |
| **T5** | **FP** | ❌ **FAIL** | **m_c/m_u drift 99.8%** |
| **T6** | **FP** | ❌ **FAIL** | **m_b/m_d drift 99.9%** |
| **T7** | **FP** | ❌ **FAIL** | **m_t/m_c drift 99.3%** |
| **T8** | **FP** | ❌ **FAIL** | **m_s/m_d drift 95.0%** |
| **T9** | **FP** | ❌ **FAIL** | **m_b/m_t drift 4031%** |
| T10 | FP | ❌ FAIL | 0/6 required g_0_q in audit |
| T11 | FP | ❌ FAIL | ceiling 2.68× vs req 80,000× |
| T12 | LIT | ✅ PASS | range consistency verified |
| T13 | DEC | ✅ PASS | S05 preserved |

**Total: 6/13 PASS; 5 central quark tests ALL FAIL → HALT-B unambiguous.**

## §6 — Risk flags update

| R# | Pre-cycle | Post-Phase-1 |
|---|---|---|
| R1: A_tail power-law fit | Acceptable | RESOLVED — R²=0.999, fit not the issue |
| R2: PDG scheme dependence | ~30% light quarks | NOT RESOLVING — drifts 95-4031% >> 30% scheme |
| R3: sek08b:529 normalization | Possible misread | RESOLVED — range verbatim verified (T12 PASS) |
| R4: CKM mixing not modeled | Out of scope | NOT RESOLVING — Cabibbo angle 0.22 cannot fix 4 orders of magnitude |

**No R-flag explains the failure structurally.** Failure is genuine: universal formula's
predictive content insufficient when g_0 constrained to audit range.

## §7 — Physics interpretation

### §7.1 — What audit sek08b:529 implies vs what PDG shows

**Audit claim:** "Ten sam ODE działa na leptony i kwarki (g_0 ∈ [0.817, 0.891])"
— SAME ODE, kwarki w narrow g_0 range bracketing electron.

**Lepton calibration (from why_n3):**
- e: g_0 = 0.869 (INSIDE [0.817, 0.891])
- μ: g_0 = 1.407 (OUTSIDE — significantly above)
- τ: g_0 = 1.755 (OUTSIDE — well above)

**If quarks are ALL in [0.817, 0.891],** they'd be near-electron-mass-scale solitons.
**But PDG shows quarks span 4.23 m_e to 337,945 m_e — 5 orders of magnitude.**

**Conclusion:** audit claim is **INCOMPATIBLE** z observed PDG quark masses **IF**
universal formula `m = c·A²·g_0^(e²/2)` is correct **AND** A_tail(g_0) follows
lepton-calibrated ODE structure.

### §7.2 — Three possible interpretations

| Interpretation | Implication |
|---|---|
| **(α) Audit range sek08b:529 wrong/misinterpreted** | Different normalization; quark g_0 actually broader; CORE inconsistency requires audit |
| **(β) Universal formula F1 incorrect for quarks** | Per-family g_0 OR quark-specific kink topology; extension cycle required |
| **(γ) Universal substrate Φ insufficient for quarks** | Multi-field substrate OR Z₂ → SU(N); FOUNDATIONS §1 review required (FORBIDDEN by recovery scope) |

**Most parsimonious (within forbidden-S05-violation scope):** (β) — universal-Φ-kink
description is INSUFFICIENT for quark sektor as written; extension cycle required.

### §7.3 — Connection z sek07_predykcje:296

Sek07 lin. 296 already acknowledges:

> "Masy kwarków (m_b, m_t) przeniesione z domyślnej predykcji do odzyskań (R),
>  ponieważ używają dopasowanego m_0 (R12, otwarte)."

**Pre-existing core acknowledgment:** quark masses are NOT default TGP predictions —
they require fitted m_0 (R12 OPEN). **This cycle CONFIRMS that acknowledgment**
empirically: universal formula z audit range constraint cannot derive quark masses
without additional structure.

## §8 — L08 audit problem #3 partial closure

**Audit problem #3 (klaster D ontology):**

> "Kwarki (g_0 ∈ [0.817, 0.891]) są 'uniwersalne' via ten sam ODE substratowy
>  (sek08b:529), ale explicit predykcje mas kwarków NIE są w PREDICTIONS_REGISTRY.
>  Neutrina (Σm_ν) są w D01 jako anchor lock, nie jako derywacja. Bozony
>  cechowania (W, Z, gluon) nie mają realizacji w warstwie 3c."

**This cycle partially addresses problem #3** (quark component only):

- ✓ **Quark component:** EMPIRICAL EVIDENCE that universal formula z audit range is
  INSUFFICIENT for quark masses — strukturalna granica zidentyfikowana
- ⏳ Neutrino component: separate downstream cycle (D01 anchor → derivation?)
- ⏳ Boson component: separate downstream cycle (warstwa 3c gauge bosons?)

**Result for #3 quark sub-problem:** STRUCTURAL_LIMITATION_DOCUMENTED — empirically
verified that universal-Φ-kink with sek08b:529 audit range cannot reproduce PDG
quark spectrum. Future paths documented (§7.2 interpretations α, β; γ forbidden).

## §9 — Lessons learned

1. **Pre-registered falsifier worked as designed** — HALT-B verdict reached
   structurally, not numerologically. No Lakatos retreat to per-family fitting.

2. **Structural ceiling test (T11) is decisive** — independent of A_tail calibration
   details, max achievable ratio (2.68×) vs required (80,000×) gives 5-orders-of-
   magnitude gap. **Definitive structural test.**

3. **Honest 6/13 PASS reflects empirical reality** — not all cycles produce A−.
   Some cycles produce empirical NULLs that are valuable scientific output.

4. **Core sek07:296 acknowledgment confirmed empirically** — m_b/m_t "recovery R12
   open" was correctly flagged; this cycle provides explicit structural reason
   (universal formula structural ceiling, not just fitting choice).

5. **Audit sek08b:529 universality claim requires revision** — current text states
   leptons + quarks share ODE with g_0 ∈ [0.817, 0.891], but empirically quarks
   would require g_0 OUTSIDE this range (free-fit: 0.99–1.38 for light quarks;
   numerical degenerate for heavy). Audit refinement candidate.

6. **B+ alternatives explored, exhausted** — light vs heavy subsector split doesn't
   help: even smallest ratio m_s/m_d=20 fails by 95% drift. No subsector recovery.

## Cross-references

- [[./README.md]] — kickoff contract z pre-registered falsifier
- [[./Phase0_balance.md]] — 8/8 ☑ gate
- [[./Phase1_sympy.py]] — empirical test script (13 sub-tests, 0 hardcoded)
- [[./Phase1_sympy.txt]] — full PASS/FAIL output
- [[./Phase_FINAL_close.md]] — closure ceremony z HALT-B verdict
- [[../why_n3/PHASE2_n_alpha_derivation.md]] — universal formula + lepton calibration source
- [[../why_n3/PHASE4_5_yukawa_propagator.md]] §4.1 — τ outlier documented
- [[../op-L05-mass-exponent-k-alpha-d-2026-05-16/Phase_FINAL_close.md]] — m_obs LIVE
- [[../op-L08-Phase6-e2-derivation-2026-05-16/Phase_FINAL_close.md]] — β(α=2)=e²/2 LIVE
- [[../../audyt/L08_kink_fermion_closure/README.md]] — problem #3 (quark component addressed)
- [[../../core/sek08b_ghost_resolution/sek08b_ghost_resolution.tex]] lin. 529 — audit range source
- [[../../core/sek07_predykcje/sek07_predykcje.tex]] lin. 296 — m_b/m_t R12 acknowledgment
