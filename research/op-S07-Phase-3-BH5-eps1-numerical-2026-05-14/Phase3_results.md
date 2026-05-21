---
title: "Phase 3 Results — numerical projections + family discriminability matrix + 4-way M9.1'' anchor"
date: 2026-05-14
type: phase-results
parent: "[[./README.md]]"
phase: 3
phase_focus: "Numerical projections per family + LIGO-O5/CE/ngEHT/LISA σ family discriminability + cross-cycle 4-way anchor matrix"
status: 🟢 CLOSED-PHASE-3 — Phase 3 complete; ready dla Phase FINAL closure ceremony
predecessors:
  - "[[./Phase1_results.md]] (BH5 KEY DERIVATION)"
  - "[[./Phase2_results.md]] (ε.1 KEY DERIVATION + cross-channel ratio invariant)"
  - "[[./Phase3_setup.md]] (10 tests scoped; family matrix + 4-way anchor)"
  - "[[./Phase3_sympy.py]] + [[./Phase3_sympy.txt]] (executed; 10/10 PASS)"
sympy_total: "34/34 PASS cumulative (Phase 1: 12 + Phase 2: 12 + Phase 3: 10)"
substance_metrics: "28 FP (82.4%) + 6 LIT (17.6%) + 0 DEC counted; 6 DEC structural separate (DEC-1..6); 0 hardcoded"
key_finding: "Family discriminability matrix: LIGO-O5 discriminates poly/quad/trans for poly-vs-others (5σ); CE 2030+ promotes quad-vs-trans to 5σ; ngEHT alone INSUFFICIENT at recovery boundary; 4-way M9.1'' anchor consistency PASSED"
tags:
  - phase-3-results
  - numerical-projections
  - family-discriminability-matrix
  - 4-way-anchor-passed
  - LISA-projection
---

# Phase 3 Results — numerical projections + family discriminability matrix

> **Phase 3 closure:** 10/10 sympy PASS, 8 FP (80.0%), 0 hardcoded. **Cumulative cycle: 34/34 PASS, 28 FP (82.4%), 0 hardcoded.** Family discriminability matrix per detector reveals: LIGO-O5 stack100 discriminates poly/quad/trans for 2/3 family pairs (poly-quad 6.4σ, poly-trans 5.5σ, quad-trans 0.88σ); Cosmic Explorer 2030+ promotes ALL 3 pairs to 5σ. **4-way M9.1'' anchor matrix at α=-4 effective PASSED** (BH5 + ε.1 + S07-reset α/3 + emergent-metric c_0·κ_σ=4/3 simultaneously consistent).

## §1 — Executive summary

**Sympy substance metrics (Phase 3):** 10/10 PASS, 8 FP (80.0%), 2 LIT (T8 LIGO PSD + T9 Sgr A* literature anchors), 0 hardcoded.

**Cumulative cycle metrics post-Phase-3:**

| Metric | Phase 1 | Phase 2 | Phase 3 | **Cumulative** |
|---|---:|---:|---:|---:|
| Sympy PASS | 12/12 | 12/12 | 10/10 | **34/34** |
| FIRST_PRINCIPLES | 10 (83.3%) | 10 (83.3%) | 8 (80.0%) | **28 (82.4%)** |
| LITERATURE_ANCHORED | 2 | 2 | 2 | 6 (17.6%) |
| Hardcoded T_pass=True | 0 | 0 | 0 | **0 preserved** |
| DEC structural separate | 2 | 2 | 2 | 6 (DEC-1..6) |

**Per-restart-era benchmark:** S07-reset cumulative 27/27, 81.5% FP; inflation cumulative 41/41, 80.5% FP. **This cycle 34/34, 82.4% FP — comparable; clean execution preserved.**

## §2 — Family discriminability matrix (per detector per family pair)

**Fiducial values (post-PE z S07-reset PR-010 recovery boundaries):**
- α = ±0.832 (1σ GWTC-3 boundary; transcendental channel)
- β_q = ±0.4 (1σ derived; quadratic channel pre-bounded)
- κ_geom = 1.0 (M9.1'' anchor upper bound; conservative for discriminability)
- Δψ_ringdown = 0.20 (op-bh-alpha-threshold T3.2)

**BH5 channel signal magnitudes at boundary fiducials:**

| Family | δω_QNM/ω_GR (BH5) | Channel value |
|---|---|---|
| Polynomial | 0 EXACT | null channel |
| Quadratic (β_q=+0.4) | κ_geom·0.4·0.04 = **1.6%** | β_q-linear |
| Transcendental (α=+0.832) | κ_geom·0.832²·0.04/2 = **1.38%** | α²-quadratic |

**ε.1 quad channel signal magnitudes at boundary fiducials:**

| Family | δε_ph²/ε_ph²_GR (quad) | Channel value |
|---|---|---|
| Polynomial | 0 EXACT | null channel |
| Quadratic (β_q=+0.4) | 0.4/9 ≈ **4.44%** | β_q-linear |
| Transcendental (α=+0.832) | 0.832²/18 ≈ **3.85%** | α²/2 |

**BH5 family discriminability matrix per detector (T1-T2-T5):**

| Detector | σ_BH5 | poly-quad | poly-trans | quad-trans | Conclusion |
|---|---|---|---|---|---|
| **LIGO-O5 stack100** | 0.25% | **6.4σ** ✅ | **5.5σ** ✅ | 0.88σ ❌ | **2/3 pairs discriminable**; quad-vs-trans NIE |
| **Cosmic Explorer stack100** | 0.025% | **64σ** ✅ | **55σ** ✅ | **8.8σ** ✅ | **ALL 3 pairs discriminable** (×10 vs LIGO-O5) |
| **LISA 2035+ EMRI** | 0.1% | **16σ** ✅ | **14σ** ✅ | 2.2σ ❌ | poly vs other 2 ✅; quad-vs-trans needs CE |

**ε.1 channel discriminability ngEHT (T3):** ALL signals <5σ at recovery boundary fiducials (4.44%, 3.85%, 0.59% vs σ_stack=6.3%) → **ngEHT alone INSUFFICIENT** for family discrimination at boundary. Combined z BH5 cross-channel (T4) marginal improvement.

**Cross-channel coupled BH5+ε.1 joint SNR (T4, transcendental family):**
- BH5 alone: 5.5σ
- ε.1 alone: 0.61σ
- Joint (Pythagorean): √(5.5² + 0.61²) ≈ **5.53σ** (small improvement; BH5-dominated)

## §3 — 4-way M9.1'' anchor matrix at α=-4 effective (T6 — KEY CROSS-CYCLE TEST)

**Independent anchors derived in 4 different cycles, all evaluated at M9.1'' effective parameter values:**

| # | Anchor | Source | Value at α=-4 / d²f=8 | Status |
|---|---|---|---|---|
| 1 | **BH5 trans channel** | This cycle Phase 1 | δω/ω = κ_geom·8/2·0.04 = κ_geom·0.16 ∈ [8%, 16%] for κ_geom∈[0.5, 1.0] | ✅ matches op-bh-alpha-threshold T3.2 LIVE 8-16% |
| 2 | **ε.1 quad channel** | This cycle Phase 2 | δε_ph²/ε_ph²_GR (quad) = 4/9 ≈ 44.4% | ✅ family-discriminator (distinct z linear total +14.6%) |
| 3 | **S07-reset Δe_2 = α/3** | S07-reset Phase 2 T10 | Δe_2 (α=-4) = **-4/3 EXACT** | ✅ S07-reset Phase 2 LOCK |
| 4 | **c_0·κ_σ = 4/3 EXACT** | emergent-metric Phase 4 Path 2 anchor | c_0·κ_σ = 4π·1/(3π) = **4/3 EXACT** | ✅ structural |

**S07-reset constraint at α=-4 (consistent z c_0·κ_σ=4/3 LOCK):**
```
-4·ξ_3 + 4 - a_3/8 + (4/3) = α/3 = -4/3
→ -4·ξ_3 - a_3/8 = -4 - 4/3 - 4/3 = -20/3   (1-param {ξ_3, a_3} family)
```

**4-way anchor consistency: ✅ PASSED** — wszystkie 4 independently derived w 4 different cycles, simultaneously consistent at M9.1'' anchor point. Cross-cycle framework coherence demonstrated.

## §4 — Three-layer L1/L2/L3 final summary (per PPN_AS_PROJECTION §3.1 BINDING)

### §4.1 — L1 Native predictions (PRIMARY)

**Native observables:**
1. δω_QNM/ω_GR (BH5 ringdown frequency relative shift) per family
2. δε_ph²/ε_ph²_GR (ε.1 photon ring quadrant relative shift) per family

**Symbolic predictions per family (Phase 1+2 derivation):**

| Family | δω_QNM/ω_GR | δε_ph²/ε_ph²_GR (quad) |
|---|---|---|
| Polynomial | 0 (null channel) | 0 (null channel) |
| Quadratic | κ_geom · β_q · (Δψ)² | β_q/9 |
| Transcendental | κ_geom · α² · (Δψ)² / 2 | α²/18 |

**Cross-channel ratio invariant (transcendental):**
```
BH5(trans) / ε.1(trans) = 9 · κ_geom · (Δψ_ringdown)²    [α CANCELS]
```

### §4.2 — L2 PPN/ppE projection consistency map

**Cross-channel coupling matrix per family** (preserved from Phase 2 §2.2):

| Family | ppE inspiral β_ppE | BH5 ringdown δω/ω | ε.1 quad δε_ph² | Coupling |
|---|---|---|---|---|
| Polynomial | (15/16)·α | 0 | 0 | inspiral-only (orthogonal) |
| Quadratic | (15/16)·α | κ_geom·β_q·(Δψ)² | β_q/9 | inspiral via α; ringdown+ε.1 via β_q (independent) |
| Transcendental | (15/16)·α | κ_geom·α²·(Δψ)²/2 | α²/18 | all 3 couple via shared α |

### §4.3 — L3 Falsification map (post-Phase-3 NUMERICAL)

| Bound | Constrains | Window | Status |
|---|---|---|---|
| op-bh-alpha-threshold T3.2 BH5 LIVE δf/f∈[8%, 16%] at M9.1'' | M9.1'' anchor | 8-16% | ✅ **PASSED** Phase 1 T7 + Phase 3 T6 |
| op-eht +14.6% photon ring shadow shift M9.1'' (total) | linear-dominated; quad component this cycle = 4/9 | quad-only family-discriminator | ✅ honest scope annotation |
| LIGO-O5 A+ ~2027 stack100 σ_BH5=0.25% | poly-quad 6.4σ; poly-trans 5.5σ ✅; quad-trans 0.88σ | 2/3 pairs 5σ-discriminable | **PENDING-DATA** ~2027 |
| Cosmic Explorer ~2030 stack100 σ_BH5=0.025% | ALL 3 pairs 8.8σ-64σ | ALL 5σ-discriminable | **PENDING-DATA** ~2030 |
| LISA 2035+ EMRI σ_BH5=0.1% | poly vs other 2 14σ-16σ; quad-trans 2.2σ | 2/3 pairs 5σ; CE remains needed | **PENDING-DATA** ~2035 |
| ngEHT ~2030 10-SMBH stack σ_eps=6.3% | ALL ε.1 signals <5σ at boundary | ngEHT alone INSUFFICIENT | **PENDING-DATA + needs CE coupling** |
| Cross-channel ratio invariant 9·κ_geom·(Δψ)² (trans) | κ_geom geometric, NIE α | future joint BH5+ε.1 measurement | **PENDING-DATA pre-observational test** |
| 4-way M9.1'' anchor matrix | BH5 + ε.1 + S07-reset α/3 + c_0·κ_σ=4/3 | all 4 simultaneously consistent | ✅ **PASSED** Phase 3 T6 |

## §5 — Substance per test detail (Phase 3)

**Anti-pattern compliance:** 0/10 hardcoded T_pass=True; 100% non-trivial substance.

| # | Test | Class | Symbolic check |
|---|---|---|---|
| T1 | BH5 family discriminability LIGO-O5 (poly-quad 6.4σ, poly-trans 5.5σ, quad-trans 0.88σ) | **FP** | `sp.Abs(bh5_poly - bh5_quad) / sigma` numerical evaluation per fiducial |
| T2 | BH5 Cosmic Explorer (all 3 pairs 8.8σ-64σ) | **FP** | same with σ_CE=0.025% |
| T3 | ε.1 ngEHT (all <5σ at boundary; INSUFFICIENT alone) | **FP** | per-pair sigma calculation |
| T4 | Cross-channel coupled BH5+ε.1 joint SNR (5.53σ vs 5.5σ alone) | **FP** | Pythagorean joint SNR |
| T5 | LISA 2035+ EMRI (poly-quad 16σ, poly-trans 14σ, quad-trans 2.2σ) | **FP** | same with σ_LISA=0.1% |
| T6 | **4-way M9.1'' anchor matrix consistency** | **FP** | symbolic verification per anchor + cross-anchor constraint |
| T7 | Recovery region α∈[-0.832, 0.832] PASSING per channel | **FP** | β_ppE=0.78 boundary, BH5 ≤1.4%, ε.1 ≤3.85% within all detectors |
| T8 | LIGO-O5 PSD-coupled SNR_BH5 at M=30 M⊙ (BH5 LIVE) | **LIT** | SNR=20, σ=2.5% per BH5 LIVE T3.2 |
| T9 | Photon orbit numerical Sgr A* θ_shadow≈50 μas | **LIT** | M=4.1·10⁶ M⊙, D=8.2 kpc → θ ≈ 50 μas literature |
| T10 | Pre-observational 5σ thresholds per detector | **FP** | LIGO-O5=1.25%, CE=0.125%, LISA=0.5%, ngEHT=31.6% |

## §6 — Anti-Lakatos PR-012 compliance Phase 3

| Sub-check | Status |
|---|---|
| Recovery scope α∈[-0.832, 0.832] preserved | ✅ PASS (T7 verifies) |
| β_q channel pre-bounded [-0.4, 0.4] | ✅ PASS (T1+T3 use boundary fiducials) |
| Brak post-hoc rule revision | ✅ PASS |
| Brak H1c/H1d backstops | ✅ PASS (only H1a recovery success vs H1b structural) |
| S05 violation NIE | ✅ PASS (DEC-5+DEC-6) |
| Direct Φ-quantum exchange NIE | ✅ PASS (Phase 1+2 T9 inheritance) |

**Anti-Lakatos: 6/6 sub-checks PASS przez Phase 3.**

## §7 — Verdict draft

**Phase 3 verdict: H1a CONFIRMED (pre-observational discriminability established)**

Combined Phase 1+2+3 demonstrates:
1. **Symbolic family marker mapping:** BH5 + ε.1 channels per S07 family verified analytically
2. **Cross-channel coupling matrix:** poly orthogonal; quad independent; trans coupled via shared α
3. **Cross-channel ratio invariant** (NEW Phase 2): pre-observational discriminator α-independent
4. **Family discriminability matrix:** LIGO-O5+CE+LISA+ngEHT sensitivity envelopes mapped per family pair
5. **4-way M9.1'' anchor matrix consistency PASSED:** cross-cycle coherence demonstrated
6. **Recovery region α∈[-0.832, 0.832] PASSING per channel:** anti-Lakatos preserved

**Per PR-012 LOCKED decision rule:** family discrimination via combined LIGO-O5+CE+ngEHT detectors achievable pre-observationally; pending observational verification 2027+ (LIGO-O5 first decisive era for poly-vs-others; CE 2030+ for full family).

## §8 — Phase FINAL entry gate

**Phase 3 → Phase FINAL transition criteria:**
1. ✅ Phase 3 sympy 10/10 PASS, 8 FP (80.0% > 70% Phase 3 target)
2. ✅ Cumulative 34/34 PASS, 28 FP (82.4% > 75% binding)
3. ✅ Three-layer L1/L2/L3 final per Phase 3 §4
4. ✅ Cross-cycle inheritance preserved 9/9 + 4-way M9.1'' anchor PASSED
5. ✅ Anti-Lakatos PR-012 6/6 sub-checks PASS (Phase 1+2+3 all)
6. ✅ ASK-RULE 4/4 Triggers PASS preserved
7. ✅ Substantive Phase 3 findings: discriminability matrix + 4-way anchor consistency

**Phase FINAL entry: AUTHORIZED.**

## §9 — Sign-off

**Phase 3 closed:** 2026-05-14 sesja P3-numerical (Claudian).

**Status:** 🟢 **CLOSED-PHASE-3.** **Cumulative:** 34/34 PASS, 28 FP (82.4%), 6 LIT, 6 DEC separate, 0 hardcoded.

**Substantywne Phase 3 findings:**
1. **Discriminability matrix per detector** (LIGO-O5/CE/ngEHT/LISA per family pair)
2. **Cross-channel coupled BH5+ε.1 joint SNR** (BH5-dominated; ε.1 marginal contribution)
3. **4-way M9.1'' anchor consistency PASSED** (BH5 + ε.1 + S07-reset α/3 + c_0·κ_σ=4/3 simultaneously consistent)
4. **Recovery region α∈[-0.832, 0.832] PASSING per channel per detector** (anti-Lakatos preserved)

**Cross-references:**
- [[./README.md]] (BINDING contract; PR-012)
- [[./Phase1_results.md]] (BH5 KEY DERIVATION)
- [[./Phase2_results.md]] (ε.1 KEY DERIVATION + cross-channel ratio invariant)
- [[./Phase3_setup.md]] (sympy substance plan)
- [[./Phase3_sympy.py]] + [[./Phase3_sympy.txt]] (10/10 PASS, 8 FP)
- Phase FINAL closure ceremony: [[./Phase_FINAL_close.md]] (combined same session per Opcja A heroic)
