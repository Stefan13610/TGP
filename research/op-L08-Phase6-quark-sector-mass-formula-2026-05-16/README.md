---
title: "op-L08-Phase6-quark-sector-mass-formula — test why_n3 universal mass formula on quark sektor (PDG 5 niezależnych mass ratios)"
date: 2026-05-16
type: research-cycle
folder_status: parking
parent: "[[../../TGP_FOUNDATIONS.md]]"

# ============== KICKOFF CONTRACT (BINDING post-2026-05-10) ==============
contract:
  # --- L1: Native (MANDATORY) ---
  L1_native:
    output_observable: "Quark mass ratios m_c/m_u, m_b/m_d, m_t/m_c, m_s/m_d, m_b/m_t [dimensionless]; ich predykcja z universal mass formula why_n3 Phase 5 (m_obs = c_M · A_tail² · g_0^(e²/2) for α=2) z g_0_q values constrained do audit range sek08b:529 [0.817, 0.891] gdzie A_tail(g_0) is ODE-derived (calibration: same R3 ODE jak dla leptonów)"
    measurement_instrument: "PDG 2024 quark masses (MS-bar 2 GeV for light; m(m) for c, b; pole for t): m_u=2.16(7) MeV, m_d=4.67(7) MeV, m_s=93.4(8.6) MeV, m_c=1.270(2) GeV, m_b=4.18(3) GeV, m_t=172.69(30) GeV; inheriting lepton anchor PDG m_e=0.5110 MeV, m_μ/m_e=206.7682, m_τ/m_e=3477.23"
    native_coefs_constrained:
      - "g_0_q per quark family (u, d, s, c, b, t) — values testowane vs audit sek08b:529 range"
      - "A_tail(g_0) functional form — ODE-derived, calibrated z 3 lepton points (e, μ, τ)"
      - "Universal exponent e²/2 ≈ 3.69 (L08-e² β(α=2)=e²/2 canonical, LIVE)"
      - "Mass anchor c_M (cancels w ratios)"
    falsification_rule: "Jeśli why_n3 Phase 5 universal mass formula m_obs = c_M · A_tail² · g_0^(e²/2) z g_0_q values constrained do audit range [0.817, 0.891] (sek08b:529) NIE reprodukuje ≥3 z 5 niezależnych quark mass ratios w tolerancji 10% (z A_tail(g_0) calibrated z lepton ODE structure), universal-Φ-kink description warstwy 3c jest INSUFFICIENT dla quark sektora → strukturalna konieczność extension (per-family g_0 OR quark-specific kink topology OR multi-family substrate)."
    pre_registration_date: "2026-05-16"

  # --- L2: Cross-framework reduction (OPTIONAL — last stage) ---
  L2_framework_reduction:
    target_frameworks:
      - "Standard Model Yukawa couplings (y_q · v Higgs)"
      - "CKM mixing-induced effective mass spectrum"
    reduction_type: "not-attempted"
    validation_transfer: ""
    failure_disposition: "L1-stands"

  # --- L3: Falsification map (consistency) ---
  L3_falsification_map:
    - { bound: "PDG m_c/m_u = 588 (light/charm)", constrains: "(A_c/A_u)² · (g_0_c/g_0_u)^(e²/2) within audit range", window: "10% tolerance", status: "pending Phase 1" }
    - { bound: "PDG m_b/m_d = 894 (light/bottom)", constrains: "(A_b/A_d)² · (g_0_b/g_0_d)^(e²/2)", window: "10%", status: "pending Phase 1" }
    - { bound: "PDG m_t/m_c = 136 (charm/top)", constrains: "(A_t/A_c)² · (g_0_t/g_0_c)^(e²/2)", window: "10%", status: "pending Phase 1" }
    - { bound: "PDG m_s/m_d = 20 (down/strange)", constrains: "(A_s/A_d)² · (g_0_s/g_0_d)^(e²/2)", window: "10%", status: "pending Phase 1" }
    - { bound: "PDG m_b/m_t = 0.024 (bottom/top)", constrains: "(A_b/A_t)² · (g_0_b/g_0_t)^(e²/2)", window: "10%", status: "pending Phase 1" }
    - { bound: "Lepton anchor m_μ/m_e = 206.77, m_τ/m_e = 3477", constrains: "Formula sanity (inheritance LIVE z L05+why_n3)", window: "0.1%", status: "inherited LIVE — sanity check only" }

# ============== END KICKOFF CONTRACT ==============

tgp_status:
  level: L1
  kind: derivation
  output_type: observable
  core_compatibility: review-only
  may_edit_core: false
  has_needs_file: false
  has_findings_file: false
  exports_findings: false
  open_bridges:
    - "L08 (kink fermion closure) — problem #3 (quarks/neutrinos/bozony) partial component"
    - "L08-Dirac-Wilson — Wilson framework reusable jako boundary on δa_e-like quark corrections"
  depends_on:
    - "research/why_n3 Phase 1-5 closed 2026-05-01 (universal mass formula)"
    - "op-L05-mass-exponent-k-alpha-d-2026-05-16 A− (m_obs ≠ M_full, k_obs(α=2,d=3)=3)"
    - "op-L08-Phase6-e2-derivation-2026-05-16 B+ (β(α=2)=e²/2 canonical)"
    - "op-L08-Phase6-FR-antisymmetry-2026-05-16 A− (Pauli + Fock antisym)"
    - "op-L08-Phase6-Clifford-emergence-2026-05-16 A− (Cl(1,3) algebra)"
    - "op-L08-Phase6-Dirac-propagator-2026-05-16 A− (S_F^TGP)"
    - "core/sek08b_ghost_resolution lin. 529: g_0_q ∈ [0.817, 0.891] universal ODE"
  impacts:
    - "audyt/L08_kink_fermion_closure problem #3 (quarks) — partial decision"
    - "TGP_FOUNDATIONS §4 warstwa 3c — partial-(D) status for matter sektor"
    - "core/sek07_predykcje lin. 296 m_b/m_t recovery R12 — informed (potential closure or strengthening)"
    - "PREDICTIONS_REGISTRY — potential PR-014 (quark mass formula falsifier)"
  source_of_status:
    - "audyt/L08_kink_fermion_closure/README.md problem #3 OPEN (multi-session)"
    - "core/sek08b_ghost_resolution lin. 529 universal ODE claim"

predecessors:
  - "[[../op-L05-mass-exponent-k-alpha-d-2026-05-16/]] (A−; m_obs ≠ M_full; k_obs(α=2,d=3)=3)"
  - "[[../op-L08-Phase6-FR-antisymmetry-2026-05-16/]] (A−; antisym Fock; γ_exchange=π)"
  - "[[../op-L08-Phase6-Clifford-emergence-2026-05-16/]] (A−; {γ^μ,γ^ν}=2η^μν)"
  - "[[../op-L08-Phase6-Dirac-propagator-2026-05-16/]] (A−; S_F^TGP standard form)"
  - "[[../op-L08-Phase6-Dirac-precision-Wilson-coefs-2026-05-16/]] (A−; Wilson framework reusable)"
  - "[[../op-L08-Phase6-e2-derivation-2026-05-16/]] (B+; β(α=2)=e²/2)"
  - "[[../why_n3/]] Phase 1-5 closed (universal mass formula 2026-05-01)"

related:
  - "[[../../audyt/L08_kink_fermion_closure/README.md]] (problem #3)"
  - "[[../../audyt/L05_mass_exponent_drift/README.md]] (CLOSED-RESOLVED via L05 cycle)"
  - "[[../../meta/CYCLE_KICKOFF_TEMPLATE.md]] §1-2 (BINDING contract)"
  - "[[../../meta/CALIBRATION_PROTOCOL.md]] (Phase 0 8/8 gate)"
  - "[[../../meta/PRE_REGISTERED_FALSIFIERS.md]] (PR-014 candidate)"
  - "[[../../core/sek08b_ghost_resolution/sek08b_ghost_resolution.tex]] lin. 529 (quark ODE universality claim)"
  - "[[../../core/sek07_predykcje/sek07_predykcje.tex]] lin. 296 (m_b/m_t recovery R12)"

classification: DERIVATION — L08 audit problem #3 quark component
priority: high (P2 OPEN klaster D ontology; partial closure dla quark sub-problem of #3)
goal: "Empirical test why_n3 Phase 5 universal mass formula m_obs = c_M · A_tail² · g_0^(e²/2) na quark sektorze z g_0_q ∈ [0.817, 0.891] audit range sek08b:529 constraint. Five independent quark mass ratios (m_c/m_u, m_b/m_d, m_t/m_c, m_s/m_d, m_b/m_t) compared to PDG 2024 z 10% tolerance threshold. Decision trichotomous: A− (≥3/5 ratios within 10%), B+ (light or heavy subsector partial), HALT-B (Path C structural insufficiency → universal-Φ-kink description requires extension)."
estimated_effort: "~1 sesja (Phase 0 + Phase 1 sympy + Phase FINAL)"
target_window: "Phase 1: sympy 13 sub-tests (T1-T13); honest verdict per actual ratio reproduction outcomes; pre-registered tolerance 10% BINDING."

six_requirements_target:
  - "P1: Universal formula structure verification — m = c_M · A_tail² · g_0^(e²/2) for α=2 canonical (T1 FP)"
  - "P2: Mass ratio formula derivation — m_i/m_j = (A_i/A_j)² · (g_i/g_j)^(e²/2) (T2 FP)"
  - "P3: Lepton anchor sanity preservation — m_μ/m_e and m_τ/m_e match PDG <0.1% (T3 LIT inheritance)"
  - "P4: Five PDG quark mass ratios predicted within 10% tolerance — m_c/m_u, m_b/m_d, m_t/m_c, m_s/m_d, m_b/m_t (T5-T9 FP, central tests)"
  - "P5: Audit range cross-check — required g_0_q values consistent z sek08b:529 [0.817, 0.891] OR fail (T10, T12)"
  - "P6: S05 single-Φ axiom preservation throughout (T13 DEC)"

risk_flags:
  - "R1: A_tail(g_0) ODE-derived form may not be representable as simple power-law fit z 3 lepton points — multi-branch ODE solutions possible (NOTE: r3_alpha2_full_closure.py provides actual ODE data; this cycle uses analytical extrapolation z 3 calibration points jako proxy)"
  - "R2: PDG quark masses are renormalization-scheme dependent (MS-bar at 2 GeV for light; m(m) for c, b; pole for t) — universal formula doesn't specify scale; this cycle uses bare PDG values acknowledging scheme variation ~30% for light quarks"
  - "R3: Sek08b:529 audit range may correspond to different normalization than R3 g_0 (e.g., g_0_quark may refer to chiral subgroup parameter) — clarification deferred do separate ontology audit cycle"
  - "R4: CKM mixing not modeled — universal formula gives mass eigenstates; CKM rotation z weak eigenstates may add subleading corrections (~Cabibbo-angle scale 0.22) — out of scope this cycle"

phase_plan:
  Phase_0: "Balance sheet + 8/8 ☑ gate + scope; lepton calibration A_tail(g_0) form (power-law fit z 3 points)"
  Phase_1: "Sympy T1-T13: formula structure, mass ratio formula, lepton verification, A_tail calibration, 5 quark ratios prediction, constraint analysis, audit cross-check"
  Phase_FINAL: "Honest verdict (A−/B+/HALT-B per data); L08 audit problem #3 partial closure note; downstream impact"

tags:
  - L08
  - quark-sector
  - mass-formula
  - universal-Phi-kink
  - audit-sek08b-529
  - PDG-comparison
  - decision-trichotomous
  - claim-status-pending
  - audit-closure-candidate
  - cycle-scaffold-2026-05-16
---

# op-L08-Phase6-quark-sector-mass-formula-2026-05-16

> **Cel:** Empiryczny test why_n3 Phase 5 universal mass formula
> `m_obs = c_M · A_tail² · g_0^(e²/2)` (α=2 canonical) na **quark sektorze**
> z g_0_q values constrained do **audit sek08b:529 range [0.817, 0.891]** (universal
> ODE claim dla leptonów i kwarków). Pięć niezależnych quark mass ratios PDG 2024
> testowanych w tolerancji 10%. **Decision trichotomous: A−/B+/HALT-B (Path C).**

## §0 — Cel + native-first contract

[CITE: `meta/CYCLE_KICKOFF_TEMPLATE.md` §1; `meta/PPN_AS_PROJECTION.md` §3.1; `audyt/L08_kink_fermion_closure/README.md` problem #3]

### §0.1 — Native observable target

**Co fizycznie liczymy:**

- 5 niezależnych quark mass ratios (m_c/m_u, m_b/m_d, m_t/m_c, m_s/m_d, m_b/m_t) [dimensionless]
- predykowane z universal formula m = c_M · A_tail² · g_0^(e²/2)
- z g_0_q values w audit range [0.817, 0.891] (sek08b:529)
- A_tail(g_0) ODE-derived (kalibrowana z 3 lepton points: e, μ, τ)

**Instrument:** PDG 2024 quark masses (per kierunki MS-bar/pole) + lepton anchor sanity check.

### §0.2 — Pre-registered falsification rule

**Decision rule WRITTEN BEFORE any calculation (2026-05-16):**

> Jeśli why_n3 Phase 5 universal mass formula z g_0_q values constrained do audit range
> [0.817, 0.891] NIE reprodukuje **≥3 z 5 niezależnych quark mass ratios** w tolerancji
> **10%** (z A_tail(g_0) calibrated z lepton ODE structure), **universal-Φ-kink description
> warstwy 3c jest INSUFFICIENT dla quark sektora** → strukturalna konieczność extension
> (per-family g_0 poza audit range / quark-specific kink topology / multi-family substrate).

```
pre_registration_date: 2026-05-16
pre_registration_hash: <auto-set by git commit SHA>
recovery_scope:
  allowed_directions:
    - "Per-family g_0 variation w obrębie audit range [0.817, 0.891]"
    - "Subleading A_tail corrections (analog L05 reconciliation; max ~lepton calibration spread)"
    - "Wilson coef extensions z L08-Dirac-Wilson framework (precision — bounded 10⁻¹⁶ lab-scale)"
  forbidden_directions:
    - "Post-hoc tuning g_0_q per quark POZA audit range (NIE family-systematic, Lakatos)"
    - "Multi-field substrate (S05 violation)"
    - "Substrate Z₂ → SU(N) extension (FOUNDATIONS §1 violation)"
  if_recovery_exhausted:
    - "H1c: Universal-Φ-kink INSUFFICIENT dla quark sektora"
    - "Multi-session extension cycle (CKM coupling derivation OR quark-specific kink topology)"
    - "OR Path C audit (quark sektor permanent (H) w warstwie 3c)"
```

### §0.3 — TGP-native check (mandatory, pre-Phase-1)

Q1-Q8 checklist per `meta/CYCLE_LIFECYCLE.md` §Phase 0 README template.

- [x] **Q1 (Pattern coverage):** Reviewed `meta/TGP_NATIVE_COMPUTATIONAL_PATTERNS.md`.
      Pattern 2.7 (asymptotic matching) — A_tail(g_0) lepton calibration; Pattern 2.6 (Derrick) inherited L05
- [x] **Q2 (Red flags):** NONE — no BD-form, no m_Φ usage, no Φ-quantum carrier framing
- [x] **Q3 (Inherited LOCKs):** L05 m_obs≠M_full LIVE; L08-e² β=e²/2 LIVE; lepton calibration why_n3 LIVE
- [x] **Q4 (Standard-physics tools):** sympy algebra + power-law fit — universal; native part: TGP-specific universal formula + audit range constraint
- [x] **Q5 (m_Φ usage):** N/A — pracujemy z classical soliton mass formula, no Φ propagator
- [x] **Q6 (GR limit framing):** N/A — flat space classical mass formula
- [x] **Q7 (ASK-RULE self-check):** Form-meaning split for "kwarki w TGP":
      KWARKI w warstwa 3c są **kinki Φ z różnymi g_0_quark** (NIE SM-elementarne
      fields w QCD-style); CKM mixing emergesowy z overlap kink amplitudes
      (NIE postulowany via Yukawa-style coupling) — declarative, ASK-RULE Trigger A documented
- [x] **Q8 (BD-drift audit plan):** Manual self-audit in Phase FINAL (no Φ propagator usage, low risk)

### §0.4 — Pre-flight methodology read confirmation

**BINDING per `meta/CYCLE_KICKOFF_TEMPLATE.md` §2.6:**

- [x] Przeczytano [[../../meta/PPN_AS_PROJECTION.md]] §3.1 — three-layer L1/L2/L3 (inherited z 5 sister cycles 2026-05-16)
- [x] Przeczytano [[../../meta/TGP_NATIVE_COMPUTATIONAL_PATTERNS.md]] §1-§4 — anti-BD-drift (inherited)
- [x] Przeczytano [[../../meta/CYCLE_KICKOFF_TEMPLATE.md]] §1-§2 — kickoff contract
- [x] Przeczytano [[../../meta/CALIBRATION_PROTOCOL.md]] — 8/8 gate template
- [x] Przeczytano [[../../audyt/L08_kink_fermion_closure/README.md]] — problem #3 statement
- [x] Przeczytano [[../why_n3/PHASE2_n_alpha_derivation.md]] — universal formula details + lepton calibration data
- [x] Przeczytano [[../op-L05-mass-exponent-k-alpha-d-2026-05-16/Phase_FINAL_close.md]] — m_obs ≠ M_full LIVE inheritance

**Sign-off:** Claudian (theoretical physics agent) @ 2026-05-16 sesja Q-quark

### §0.5 — Sympy substance plan

Per `meta/AUDIT_2026-05-11_sympy_substance.md` §4 — substance-first lessons:

- [x] Każdy test sympy ma **explicit pytanie fizyczne** które weryfikuje
- [x] **≥75% testów** to non-trivial symbolic manipulation
- [x] **0 hardcoded T_pass=True** (Phase 6 ABSOLUTE BINDING)
- [x] PASS/FAIL determined empirically z explicit tolerance threshold

**Plan testów Phase 1 (13 tests target):**

| Test | Klasa | Pytanie fizyczne | Pre-decyzja PASS |
|---|---|---|---|
| T1 | **FIRST_PRINCIPLES** | Universal formula structure: m = c_M · A_tail² · g_0^(e²/2) for α=2 canonical | Exponent symbolic check |
| T2 | **FIRST_PRINCIPLES** | Mass ratio formula: m_i/m_j = (A_i/A_j)² · (g_i/g_j)^(e²/2) | Algebraic equivalence |
| T3 | **LITERATURE_ANCHORED** | Lepton verification: m_μ/m_e + m_τ/m_e match PDG <0.1% (inherited why_n3 LIVE) | Inherited PASS, sanity check |
| T4 | **FIRST_PRINCIPLES** | A_tail(g_0) power-law fit z 3 lepton points; residual quality | Fit r² > 0.99 expected |
| T5 | **FIRST_PRINCIPLES** | m_c/m_u predicted with g_0_c, g_0_u ∈ audit range [0.817, 0.891]; tolerance 10% vs PDG 588 | central test #1 |
| T6 | **FIRST_PRINCIPLES** | m_b/m_d predicted similarly; 10% vs PDG 894 | central test #2 |
| T7 | **FIRST_PRINCIPLES** | m_t/m_c predicted; 10% vs PDG 136 | central test #3 |
| T8 | **FIRST_PRINCIPLES** | m_s/m_d predicted; 10% vs PDG 20 | central test #4 |
| T9 | **FIRST_PRINCIPLES** | m_b/m_t predicted; 10% vs PDG 0.024 | central test #5 |
| T10 | **FIRST_PRINCIPLES** | Required g_0_q values (assuming A_tail calibration) to match PDG; check all 6 ∈ audit range | Constraint compatibility |
| T11 | **FIRST_PRINCIPLES** | Maximum achievable mass ratio max(m)/min(m) with g_0 ∈ [0.817, 0.891]; compare to required m_t/m_u ≈ 80,000 | Structural ceiling test |
| T12 | **LITERATURE_ANCHORED** | Cross-check sek08b:529 explicit range vs cycle assumption | Literature consistency |
| T13 | **DECLARATIVE** | S05 single-Φ preservation (single ODE, single substrate Φ) | separate from PASS count |

**Target:** sympy 13 sub-tests; **honest PASS count determined by empirical outcome** (NO predetermined pass count); decision trichotomous per pre-registered falsification rule §0.2.

**Substance ratio:** 10 FP (76.9%) + 2 LIT (15.4%) + 1 DEC (7.7%, separate). Above 75% FP threshold.

---

## §1 — Phase 0: balance sheet + 8/8 gate

[Patrz `Phase0_balance.md` w tym folderze]

## §2 — Phase 1: empirical test

[Patrz `Phase1_sympy.py` + `Phase1_sympy.txt` + `Phase1_results.md` w tym folderze]

## §FINAL — Honest verdict + L08 audit problem #3 partial closure

[Patrz `Phase_FINAL_close.md` w tym folderze; verdict A−/B+/HALT-B per actual data]

---

## Status

🟡 **PARKING — scaffold opened 2026-05-16 sesja Q-quark**. Pre-flight methodology read confirmation:
**complete** (§0.4). Authorization: user kickoff prompt (sesja Q-quark, 13. cykl 2026-05-16 — first post-Wilson cycle).

Phase 0 commit gate:
1. README.md (this file) z BINDING contract — **DONE**
2. PR-014 candidate documented w Phase_FINAL (pending verdict outcome)

This session deliverables target:
- README.md (this file) z BINDING contract — **DONE**
- Phase0_balance.md — **PLANNED**
- Phase1_sympy.py — **PLANNED** (10 FP + 2 LIT + 1 DEC; 0 hardcoded)
- Phase1_sympy.txt — **PLANNED** (full output)
- Phase1_results.md — **PLANNED**
- Phase_FINAL_close.md — **PLANNED** (honest verdict per data)

---

**Cycle scaffolded:** 2026-05-16 (Claudian, theoretical physics expert role; opens after
12 sister cycles closed 2026-05-16; addresses L08 audit problem #3 quark component
of klaster D ontology).

**Cross-references:**
- [[../../meta/CYCLE_KICKOFF_TEMPLATE.md]] §1-§2 (BINDING contract)
- [[../../meta/CALIBRATION_PROTOCOL.md]] (8/8 gate template)
- [[../../audyt/L08_kink_fermion_closure/README.md]] (problem #3 statement)
- [[../why_n3/PHASE2_n_alpha_derivation.md]] (universal formula + lepton calibration)
- [[../op-L05-mass-exponent-k-alpha-d-2026-05-16/Phase_FINAL_close.md]] (m_obs ≠ M_full LIVE)
- [[../op-L08-Phase6-e2-derivation-2026-05-16/Phase_FINAL_close.md]] (β(α=2)=e²/2 canonical)
- [[../../core/sek08b_ghost_resolution/sek08b_ghost_resolution.tex]] lin. 529 (universal ODE claim)
