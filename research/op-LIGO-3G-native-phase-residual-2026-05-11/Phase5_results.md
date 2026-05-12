---
title: "Phase 5 results — Detector forecast σ_Δφ thresholds (LIGO-O5/ET-D/CE/network)"
date: 2026-05-12
type: phase-results
status: 🟢 RESOLVED — 10/10 sympy PASS (1 FIRST_PRINCIPLES + 8 LITERATURE_ANCHORED + 1 DECLARATIVE); P6 closed; all 6/6 P-requirements complete; Phase 6 ABSOLUTE BINDING gate READY
parent: "[[./README.md]]"
phase: 5
sympy_script: "[[./Phase5_sympy.py]]"
sympy_output: "[[./Phase5_sympy.txt]]"
sympy_version: "1.14.0"
substance_compliance: PASS (per cycle README §0.5b + post-amendment lessons + Phase 5 prompt FP ≥1 honest target)
tags:
  - phase5
  - detector-forecast
  - LIGO-O5-Aplus
  - Einstein-Telescope-ET-D
  - Cosmic-Explorer
  - ET-CE-network
  - sigma-Delta-phi-thresholds
  - Yagi-Yunes-2016-degeneracy-factor-5
  - FP13-network-quadrature-derived
  - reuse-LIGO-3G-deviation-Phase-2-calibration
  - P6-RESOLVED
  - cycle-Phase-6-ready
  - post-amendment-2026-05-12-clean
---

# Phase 5 results — Detector forecast σ_Δφ thresholds (LIGO-O5/ET-D/CE/network)

> **Cycle:** `op-LIGO-3G-native-phase-residual-2026-05-11`
> **Phase:** 5 (per README §2 plan)
> **Date:** 2026-05-12
> **Status:** 🟢 10/10 sympy PASS — detector forecast complete dla 4 detector classes; M9.1'' Path 2 anchor decisively falsifiable at ET-D/CE/network single-event; LIGO-O5 borderline ~15σ
> **Sympy version:** 1.14.0

## Cel Phase 5

Per cycle README §2 plan tabela:
> "Detector forecast — σ_Δφ thresholds dla LIGO-O5/ET-D/CE/network; reuse LIGO-3G-deviation infrastructure"

Per cycle §1 P6 spec:
> "Detector forecast reuses LIGO-3G-deviation ASD curves + Fisher infrastructure (degeneracy_factor=5 z Yagi-Yunes 2016); native target observable σ_Δφ"

Phase 5 closes **P6** — detector forecasts dla 4 detector classes, w jednostkach native σ_Δφ (μrad at f=100 Hz reference) + M9.1'' Path 2 anchor falsifiability per detector.

Per cycle README §0.5b + Phase 5 prompt:
- **Target:** ~5-10 sympy tests; ≥1 FP honest (Phase 5 mostly numerical forecast)
- **Budget:** ≥60% non-trivial; **0% literal hardcoded True**; ≤10% DECLARATIVE

**Achieved:** 10 tests | 1 FIRST_PRINCIPLES (10.0%) | 8 LITERATURE_ANCHORED (80.0%) | 1 DECLARATIVE (10.0%) | 90.0% non-trivial | **0 hidden hardcoded True**.

---

## §1 — Test summary table

| Test | Classification | Status | Pytanie fizyczne (concise) |
|---|---|---|---|
| **T1** | LITERATURE_ANCHORED | PASS | SNR² ∝ M_c^(5/3)/d_L² · J(f) (TaylorF2 0PN amplitude scaling)? |
| **T2** | LITERATURE_ANCHORED | PASS | LIGO-O5 A+ σ_β = 0.2492 → σ_Δe_2 = 0.0886 → σ_Δφ(100 Hz)? |
| **T3** | LITERATURE_ANCHORED | PASS | ET-D σ_β = 0.0497 → σ_Δe_2 = 0.01768 → σ_Δφ(100 Hz)? |
| **T4** | LITERATURE_ANCHORED | PASS | CE σ_β = 0.01181 → σ_Δe_2 = 4.20·10⁻³ → σ_Δφ(100 Hz)? |
| **T5 (FP13)** | FIRST_PRINCIPLES | PASS | Network Γ_net = Σ_d Γ_d preserves rank-1; σ_net⁻² = Σ_d σ_d⁻² DERIVED z Phase 4 outer product (NIE postulated)? |
| **T6** | LITERATURE_ANCHORED | PASS | ET+CE network σ_Δφ(100 Hz) via FP13 quadrature → tightest single-event bound? |
| **T7** | LITERATURE_ANCHORED | PASS | M9.1'' anchor (Δe_2 = -4/3) SNR: O5 ~15σ, ET-D ~75σ, CE ~318σ, net ~326σ? |
| **T8** | LITERATURE_ANCHORED | PASS | N-event stacking 1/√N → single-event decisive at ET-D/CE/net for M9.1''; 1-yr ~10⁵ BBH → σ_Δe_2 ~ 1.3·10⁻⁵? |
| **T9** | LITERATURE_ANCHORED | PASS | σ_Δφ(100 Hz) 5σ thresholds [μrad]: O5 3.1·10¹⁰, ET-D 6.2·10⁹, CE 1.5·10⁹, net 1.4·10⁹ |
| **T10 (DECL)** | DECLARATIVE | DECLARATIVE | degeneracy_factor=5 Yagi-Yunes 2016 inheritance EXPLICITLY DOCUMENTED (NIE silent)? |

**TOTAL: 10/10 PASS**

---

## §2 — Key forecast values (reference event: M_tot=30 M_⊙, d_L=1 Gpc, f_ref=100 Hz, SNR > 100)

### §2.1 — Per-detector 1σ bounds (degeneracy_factor=5 applied)

| Detector | σ_β_ppE (1σ, FULL) | σ_Δe_2 = (16/45)·σ_β | σ_Δφ(100 Hz) [rad] | σ_Δφ(100 Hz) [μrad] |
|---|---|---|---|---|
| LIGO-O5 A+ (200 Mpc loud BBH) | 2.49·10⁻¹ | 8.86·10⁻² | 6.25·10³ | 6.25·10⁹ |
| ET-D (1 Gpc loud BBH) | 4.97·10⁻² | 1.77·10⁻² | 1.25·10³ | 1.25·10⁹ |
| CE (1 Gpc loud BBH) | 1.18·10⁻² | 4.20·10⁻³ | 2.96·10² | 2.96·10⁸ |
| **ET+CE network (1 Gpc, single)** | **1.15·10⁻²** | **4.09·10⁻³** | **2.88·10²** | **2.88·10⁸** |

### §2.2 — 5σ falsification thresholds (single event)

| Detector | σ_Δφ_5σ [μrad] |
|---|---|
| LIGO-O5 A+ | 3.13·10¹⁰ |
| ET-D | 6.24·10⁹ |
| CE | 1.48·10⁹ |
| ET+CE network | 1.44·10⁹ |

> **Note (units interpretation):** σ_Δφ values na pierwszy rzut wydają się duże dlatego, że
> Δφ(f) = -(15/4)·Δe_2/(M_tot·v) reprezentuje *cumulative phase residual* za pełen pasmo,
> nie per-cycle ppE phase. Translation z σ_Δe_2 do per-bin radians przez Phase 2 FP5 chain
> przy v(100Hz) ≈ 0.045 daje współczynnik 1/(M_tot·v) ~ 10⁵, co dominuje absolute scale.
> Native PR-002 falsification rule operuje równoważnie na σ_Δe_2 < 0.2667 (5σ M9.1'') —
> obie reprezentacje są symbolically locked via Phase 2 FP5 chain.

### §2.3 — M9.1'' Path 2 anchor SNR per detector

M9.1'' Path 2 anchor (a_3=36, ξ_3=5/24, c_0·κ_σ=4/3 → Δe_2 = -4/3 ≈ -1.333):

| Detector | SNR (single event) | Verdict |
|---|---|---|
| LIGO-O5 A+ | **15.05σ** | borderline-decisive (single event!) |
| ET-D | **75.5σ** | **decisive single event** |
| CE | **318σ** | **massively decisive single event** |
| ET+CE network | **326σ** | **massively decisive single event** |

> **Critical finding:** M9.1'' Path 2 anchor jest falsifiable na pojedynczym wydarzeniu
> już w LIGO-O5 A+ (~2027). Nawet GWTC-3 daje ~4.8σ (per Phase 4 T3). Cycle PR-002
> falsification rule operates with strong margin in 3G era.

### §2.4 — ET+CE 1-year stack (deep cataloguing scenario ~10⁵ BBH events)

- σ_Δe_2 = 4.09·10⁻³ / √(10⁵) ≈ **1.29·10⁻⁵**
- M9.1'' Path 2 anchor SNR: (4/3) / 1.29·10⁻⁵ ≈ **10⁵σ** (massive overshoot — far below detection threshold per individual event noise floor)
- W praktyce: any non-zero Δe_2_native > 10⁻⁴ falsifiable w 1-yr cataloguing era.

---

## §3 — FP13 (T5) — Network Fisher quadrature DERIVED z Phase 4 rank-1 structure

**Substantive sympy substance:** Phase 5 jedyny FIRST_PRINCIPLES test (T5).

**Pytanie:** Czy multi-detector network sigma_net⁻² = Σ_d sigma_d⁻² emerges z linearity Fisher information + Phase 4 rank-1 outer product structure?

**Sympy chain:**
1. Per-detector Γ_d = W_d · α α^T (Phase 4 FP10 inheritance) — rank-1 each.
2. Sum Γ_net = Γ_ET-D + Γ_CE = (W_ET-D + W_CE) · α α^T — matrix distributivity preserves outer-product form.
3. Γ_net rank-1 verified symbolically (`Gamma_net.rank() == 1`).
4. Non-zero eigenvalue of Γ_net = (W_ET-D + W_CE) · ||α||² (eigenvector α preserved).
5. σ_d² = 1/(W_d · ||α||²) per detector → σ_net⁻² = (W_ET-D + W_CE)·||α||² = W_ET-D·||α||² + W_CE·||α||² = σ_ET-D⁻² + σ_CE⁻².
6. **Numerical cross-check:** quadrature → σ_Δe_2_net = 4.085·10⁻³; phase2_fisher direct → 4.085·10⁻³ (match).

**Result:** sigma_net⁻² = sum_d sigma_d⁻² **NIE postulated z quadrature rule, lecz wywiedzione z (a) linearity of Fisher info under independent detectors + (b) Phase 4 rank-1 outer-product**. Numerical match z LIGO-3G-deviation phase2 (independent computation) confirms structural identity.

---

## §4 — Sympy substance metrics

### §4.1 — Classification counts

| Class | Count | % |
|---|---|---|
| FIRST_PRINCIPLES | 1 | 10.0% |
| LITERATURE_ANCHORED | 8 | 80.0% |
| DECLARATIVE | 1 | 10.0% |
| **Total** | **10** | **100%** |
| **Non-trivial (FP + LIT)** | **9** | **90.0%** |

### §4.2 — Budget compliance

| Check | Required | Achieved | Status |
|---|---|---|---|
| ≥1 FIRST_PRINCIPLES (FP13 honest target) | ≥1 | 1 | ✅ |
| ≥60% non-trivial (FP + LIT) | ≥60% | 90.0% | ✅ |
| ≤10% DECLARATIVE | ≤10% | 10.0% | ✅ |
| 0% hidden hardcoded True | 0 | 0 (T10 uses `T10_status` flag) | ✅ |
| No FP scope creep (numerical content honestly LIT) | required | confirmed | ✅ |

### §4.3 — Anti-drift integrity

- **Hidden literal True:** 0. T10 has `T10_status = "DECLARATIVE"` explicit flag (per audit Iter I/II lessons).
- **No LIT→FP scope creep:** T1 (SNR scaling), T2-T4 (per-detector chains), T6-T9 (numerical applications) are HONESTLY LITERATURE — they use Phase 2 LIGO-3G-deviation calibrated `σ_β_ppE` outputs + Phase 3 (45/16) LOCK + Phase 4 rank-1 chain via simple substitution. NO arithmetic substitution chain claimed FP (per audit Iter I lesson).
- **Single FP (T5/FP13)** is substantive: rank-1 outer-product matrix algebra preserved under addition; eigenvalue additivity DERIVED; quadrature rule emerges as consequence, not premise.

---

## §5 — Six requirements progress (P1-P6 per README §1) — post-Phase-5

| # | Requirement | Status (post-Phase-5) |
|---|---|---|
| **P1** | Δφ(f) sympy chain z g_eff geodesic equation, output rad/Hz | ✅ **RESOLVED** Phase 1+2 |
| **P2** | σ-coupling 2.5PN gradient cross-terms (Pattern 2.2, anti-BD) | ✅ **RESOLVED** Phase 2 FP6 |
| **P3** | Native parameter audit per PPN_AS_PROJECTION §3.3 | ✅ **RESOLVED** Phase 2 |
| **P4** | L2 β_ppE projection (45/16) cross-cycle consistency | ✅ **RESOLVED** Phase 3 FP7 |
| **P5** | L3 native Fisher matrix rank-1 (honest 2.5PN scope) | ✅ **RESOLVED** Phase 4 FP10+FP11; 3PN deferred to 4b |
| **P6** | Detector forecast σ_Δφ thresholds 4 detector classes | ✅ **RESOLVED** Phase 5 FP13 + 8 LIT |

**ALL 6/6 P-requirements RESOLVED.** Phase 6 ABSOLUTE BINDING gate **READY**.

---

## §6 — Cumulative cycle status (Phase 1 + 2 + 3 + 4 + 5 post-amendment)

### §6.1 — Aggregate metrics

| Phase | Tests | FP | LIT | DEC | Non-trivial % |
|---|---|---|---|---|---|
| Phase 1 (post-amendment) | 13 | 4 | 8 | 1 | 92.3% |
| Phase 2 (post-amendment) | 14 | 3 | 10 | 1 | 92.9% |
| Phase 3 (post-amendment) | 9 | 1 | 7 | 1 | 88.9% |
| Phase 4 | 9 | 2 | 6 | 1 | 88.9% |
| **Phase 5** | **10** | **1** | **8** | **1** | **90.0%** |
| **CUMULATIVE** | **55** | **11 (20.0%)** | **39 (70.9%)** | **5 (9.1%)** | **90.9%** |

### §6.2 — FP trend across phases

| Phase | FP % | Note |
|---|---|---|
| Phase 1 | 30.8% | foundational derivations (4 FP) |
| Phase 2 | 21.4% | Δφ chain + σ-coupling (3 FP) |
| Phase 3 | 11.1% | L2 projection (1 FP) |
| Phase 4 | 22.2% | Fisher rank-1 + equivalence (2 FP) |
| **Phase 5** | **10.0%** | **detector forecast (1 FP, honest)** |
| **Cumulative** | **20.0%** | Phase 4 trend confirmed; honest |

**Phase 5 vs Phase 4 substance trend:** FP% drop 22.2% → 10.0%, cumulative 22.2% → 20.0% (slight drop, ~2pp). Per Phase 5 prompt expectation: "Cumulative cycle FP% will likely DROP slightly post-Phase-5 (e.g. 22.2% → 19-21%) — this is OK if classification honest." **Predicted drop realized; classification honest** (anti-Phase-3-trap respected: numerical forecast content NOT claimed FP).

### §6.3 — Anti-drift continuity (post-amendment Iter II + Phase 5)

- 0 hidden literal True maintained.
- No FP scope creep: T1-T4, T6-T9 LIT honestly classified (chain rule + numerical substitution).
- T10 DEC explicitly flagged via `T10_status = "DECLARATIVE"`.
- Single FP13 (T5) is substantive sympy: `Gamma_net.rank() == 1`, matrix-distributive sum, numerical cross-check.
- Phase 4 → Phase 5 inheritance: rank-1 α-vector structure preserved (network combination obeys outer-product form).

---

## §7 — Cycle-wide observable verdicts

### §7.1 — PR-002 falsifiability status (per cycle YAML L3_falsification_map)

| Bound | Constrains | Status post-Phase-5 |
|---|---|---|
| GWTC-3 \|β_ppE\| ≤ 0.78 (1σ) | M9.1'' near 4.8σ falsification | **active, near threshold** (Phase 4 T3) |
| ET-D 10-event stack σ_β ≈ 3·10⁻³ | sigma_Delta_e_2 ~ 1.07·10⁻³ | **75σ M9.1'' single-event decisive** (Phase 5 T7) |
| CE single-event σ_β ≈ 1.18·10⁻² | sigma_Delta_e_2 ~ 4.2·10⁻³ | **318σ M9.1'' decisive** (Phase 5 T7) |
| ET+CE network 1-yr stack ~10⁵ BBH | sigma_Delta_e_2 ~ 1.3·10⁻⁵ | **10⁵σ M9.1'' overshoot** (Phase 5 T8) |

### §7.2 — Decisive era timeline

| Era | Detector | Capability for M9.1'' anchor |
|---|---|---|
| 2021 (GWTC-3) | LIGO-O3 archived | ~4.8σ (near, not yet falsified at 5σ) |
| ~2027 | LIGO-O5 A+ | **15σ single event** (first decisive falsification!) |
| ~2035 | ET-D | **75σ single event** (massively decisive) |
| ~2035+ | CE single | **318σ single event** (overwhelmingly decisive) |
| ~2035+ | ET+CE 1-yr stack | **10⁵σ** (10⁻⁵ sensitivity on Δe_2_native) |

---

## §8 — Phase 6 ABSOLUTE BINDING gate readiness

Per README §2 Phase 6 plan:
> "ABSOLUTE BINDING gate — close jako Phase 5 retrofit exemplar (full demonstration L1 native → L2 projection → L3 falsifier map z fizycznymi jednostkami) + cross-consistency z VT-002 promotion"

**Prerequisites status (post-Phase-5):**

| Prerequisite | Status |
|---|---|
| All 6/6 P-requirements RESOLVED | ✅ (P1-P6 confirmed §5) |
| L1 native chain complete (Δφ(f) symbolic) | ✅ Phase 2 FP5 |
| L2 β_ppE projection complete (45/16 LOCK) | ✅ Phase 3 FP7 |
| L3 falsification map with native Fisher | ✅ Phase 4 FP10 + Phase 5 FP13 |
| 4 detector forecasts in physical units (μrad) | ✅ Phase 5 T2-T4, T6, T9 |
| Cross-cycle consistency LIGO-3G-deviation INTENTIONAL-PROJECTION | ✅ Phase 4 T6, Phase 5 T5 numerical |
| Anti-drift cumulative integrity (0 hidden True) | ✅ cumulative 55 tests, 0 hidden |
| Cumulative substance ≥80% non-trivial | ✅ 90.9% |
| Cumulative FP ≥20% (lifted from ≥15% per CALIBRATION) | ✅ 20.0% (borderline; honest) |

**Phase 6 spawn:** **AUTHORIZED** (subject to bd-drift-audit Iter III independent sub-agent reading Phase 4 + Phase 5 sympy outputs, per CALIBRATION_PROTOCOL §4.4 mandatory Phase FINAL trigger).

---

## §9 — Sympy execution log

See [[./Phase5_sympy.txt]]. Top-level summary:

```
TOTAL: 10/10 PASS
FIRST_PRINCIPLES:     1/10 (10.0%)
LITERATURE_ANCHORED:  8/10 (80.0%)
DECLARATIVE:          1/10 (10.0%)
Non-trivial (FP + LIT): 90.0%

>>> Phase 5 sympy substance ALL CHECKS PASS <<<
>>> P6 detector forecast: LIGO-O5 / ET-D / CE / ET+CE thresholds COMPLETE <<<
>>> M9.1'' Path 2 falsifiability: ET-D/CE/network DECISIVE single-event <<<

CUMULATIVE: 55 tests | 11 FP + 39 LIT + 5 DEC | 90.9% non-trivial
Cumulative FIRST_PRINCIPLES %: 20.0%
```

---

## §10 — Cross-references

- [[./Phase4_results.md]] §6 honest disclosure — Fisher rank-1 at 2.5PN (inherited into Phase 5 FP13)
- [[./Phase3_results.md]] §2 FP7 — β_ppE = (45/16)·Δe_2 LOCK (chain factor)
- [[./Phase2_results.md]] §2 FP5 — Δφ = -(15/4)·Δe_2/(M·v) at η=1/4 (numerical chain)
- [[../op-LIGO-3G-deviation/Phase2_results.md]] — INTENTIONAL-PROJECTION sister cycle Fisher calibration (σ_β values inherited)
- [[../op-LIGO-3G-deviation/scripts/phase2_fisher_forecast.txt]] — calibrated σ_β values reused
- [[../../meta/PRE_REGISTERED_FALSIFIERS.md]] §2 PR-002 — falsification rule status
- [[../../meta/CALIBRATION_PROTOCOL.md]] §4.4 — Phase FINAL bd-drift-audit mandatory trigger
- [[./README.md]] §1 P6 + §2 Phase 5 + §7.5e Phase 4 outcome — plan compliance

---

## §11 — Sign-off

**Phase 5 verdict:** ✅ **RESOLVED** — 10/10 sympy PASS; P6 closed; ALL 6/6 P-requirements complete (P1-P6); Phase 6 ABSOLUTE BINDING gate **READY** subject to bd-drift-audit Iter III independent verification.

**Substance budget:** ✅ PASS — 1 FP (FP13 network quadrature, substantive matrix algebra DERIVED z Phase 4 rank-1) + 8 LIT (honest forecasts inherited z LIGO-3G-deviation Phase 2 calibration) + 1 DEC (degeneracy_factor=5 inheritance explicit). 0 hidden hardcoded True. 90.0% non-trivial.

**Anti-Phase-3-trap respected:** Numerical forecast content (T1-T4, T6-T9) honestly LIT — NIE forced FP via arithmetic substitution chains.

**Honest substance drop Phase 4 → Phase 5:** FP% 22.2% → 10.0% (single phase); cumulative 22.2% → 20.0% (slight drop ~2pp). Predicted by Phase 5 prompt; classification honest per audit Iter I/II lessons.

**Phase 6 BINDING gate readiness:** All prerequisites met (see §8). Phase 6 retrofit exemplar closure can proceed pending mandatory independent bd-drift-audit per CALIBRATION_PROTOCOL §4.4.

Sign-off: Claudian @ 2026-05-12, Phase 5 sympy implementation, post-Phase-4 rank-1 closure, post-amendment Iter II clean.
