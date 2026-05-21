---
title: "op-L08-Phase6-hadron-topology-confinement — composite hadron N-M≡0 mod 3 rule z compact U(1) J_phase winding"
date: 2026-05-16
type: research-cycle
folder_status: parking
parent: "[[../../TGP_FOUNDATIONS.md]]"

# ============== KICKOFF CONTRACT (BINDING post-2026-05-10) ==============
contract:
  # --- L1: Native (MANDATORY) ---
  L1_native:
    output_observable: "Hadron composition classification: dla każdej konfiguracji N quark + M antiquark, czy n_total = Σ n_i ∈ ℤ (allowed jako izolowalna cząstka) czy n_total ∉ ℤ (forbidden). Conkretnie: predicate ISOLABLE(config) ∈ {True, False} dla 14+ PDG hadron'ów + 4+ konfiguracji forbidden. Rule statement: N_q - N_q̄ ≡ 0 (mod 3) gdzie quark winding ∈ {±1/3, ±2/3}."
    measurement_instrument: "PDG 2024 hadron compositions (baryons: p, n, Δ⁺⁺, Δ⁻, Λ, Σ⁺, Σ⁻, Ξ⁰; mesons: π⁺, π⁻, π⁰, K⁺, K⁰, J/ψ); LHCb 2015 pentaquarks P_c(4380), P_c(4450); BESIII tetraquark Z_c(3900); LHCb 2021 T_cc(3875)"
    native_coefs_constrained:
      - "n_total ∈ ℤ constraint dla izolowanej cząstki z compact U(1) (dodatekO thm:winding_quant)"
      - "Quark winding assignments n_q ∈ {+2/3, -1/3} dla {u-type, d-type} z SM electric charges (input PDG)"
      - "Composition rule N_q - N_q̄ ≡ 0 (mod 3) wymuszona przez fractional winding sum integer constraint"
    falsification_rule: "Jeśli compact U(1) winding quantization z dodatekO thm:winding_quant + SM quark fractional winding (±1/3, ±2/3) NIE klasyfikuje correctly wszystkich 14+ obserwowanych PDG hadron'ów jako n_total ∈ ℤ (allowed) AND wszystkich 4+ forbidden konfiguracji (izolowany kwark, di-kwark, 4-kwark) jako n_total ∉ ℤ (forbidden), strukturalna mechanism INSUFFICIENT → wymaga albo (a) substrate extension SU(N>1) (S05 violation, FORBIDDEN by recovery scope) lub (b) revisiting quark winding assignment (acceptable). Threshold: ≥14/14 PDG hadrons must classify correctly AND ≥4/4 forbidden configs must classify correctly. Less = STRUCTURAL_INSUFFICIENCY."
    pre_registration_date: "2026-05-16"

  # --- L2: Cross-framework reduction (OPTIONAL — last stage) ---
  L2_framework_reduction:
    target_frameworks:
      - "QCD SU(3)_color confinement (Wilson 1974, color singlet rule)"
      - "Lattice QCD hadron spectrum"
    reduction_type: "structural-equivalence"
    validation_transfer: "TGP compact U(1) confinement rule matches QCD SU(3) color singlet rule N_q - N_q̄ ≡ 0 (mod 3) for hadron composition; differs in MECHANISM (topological vs energetic) but agrees in composition predictions"
    failure_disposition: "L1-stands"

  # --- L3: Falsification map (consistency) ---
  L3_falsification_map:
    - { bound: "PDG observed hadrons (all)", constrains: "composition rule N-M ≡ 0 mod 3", window: "100% match required", status: "pending Phase 1" }
    - { bound: "LHCb 2015 pentaquark P_c(4380), P_c(4450) uudcc̄", constrains: "4q+1q̄ allowed (N-M=3)", window: "structural", status: "pending Phase 1" }
    - { bound: "BESIII 2013 Z_c(3900) (tetraquark cc̄uū-like)", constrains: "2q+2q̄ allowed (N-M=0)", window: "structural", status: "pending Phase 1" }
    - { bound: "LHCb 2021 T_cc(3875) (tetraquark ccūd̄)", constrains: "structure consistency", window: "structural", status: "pending Phase 1" }
    - { bound: "NEVER OBSERVED: isolated free quark", constrains: "N-M ≡ 1 mod 3 forbidden", window: "non-observation", status: "consistency LIVE" }
    - { bound: "NEVER OBSERVED: isolated di-quark", constrains: "N-M ≡ 2 mod 3 forbidden", window: "non-observation", status: "consistency LIVE" }

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
    - "L08 audit problem #3 (quarks/neutrinos/bosons) — quark confinement sub-component"
    - "Audit sek08b:529 universalność kwarkowa — potential reframe via N-M rule"
  depends_on:
    - "core/formalizm/dodatekO_u1_formalizacja.tex thm:winding_quant (n[γ] ∈ ℤ)"
    - "core/formalizm/dodatekE_pi1_formal.tex (kink stany niosą kwantyzację)"
    - "research/op-lambda1-e2-amplitude-emergence/phase1L5_amplitude_phase_separation.md (J_amp vs J_phase split)"
    - "research/why_n3/PHASE3_RP2_defect_quantization.md (RP² topology → spin-1/2)"
    - "research/op-L08-Phase6-FR-antisymmetry-2026-05-16 (antisym Fock)"
    - "research/op-L08-Phase6-Clifford-emergence-2026-05-16 (Cl(1,3) algebra)"
    - "research/exploration_neutrino_g0_2026-05-16/topology_playground.py (motivating exploration)"
  impacts:
    - "audyt/L08_kink_fermion_closure problem #3 quark sub-component — derivation candidate (replacing HALT-B mass formula approach)"
    - "TGP_FOUNDATIONS §4 warstwa 3c — strukturalna closure dla quark confinement"
    - "PREDICTIONS_REGISTRY — hadron composition rule (PR-015 candidate)"
    - "audyt/sek07_predykcje lin. 296 m_b/m_t recovery R12 — informed (confinement explains free-quark non-observation)"
  source_of_status:
    - "audyt/L08_kink_fermion_closure problem #3 OPEN"
    - "research/exploration_neutrino_g0_2026-05-16/topology_playground.txt (motivating data: 13/13 hadrons classified correctly w playground)"

predecessors:
  - "[[../op-L08-Phase6-quark-sector-mass-formula-2026-05-16/]] (HALT-B; superseded approach — mass-formula direct fails; this cycle uses topology instead)"
  - "[[../op-L08-Phase6-FR-antisymmetry-2026-05-16/]] (A−; antisym Fock; γ_exchange=π)"
  - "[[../op-L08-Phase6-Clifford-emergence-2026-05-16/]] (A−; {γ^μ,γ^ν}=2η^μν)"
  - "[[../op-L08-Phase6-Dirac-propagator-2026-05-16/]] (A−; S_F^TGP)"
  - "[[../op-lambda1-e2-amplitude-emergence/]] (J_amp vs J_phase split)"
  - "[[../exploration_neutrino_g0_2026-05-16/]] (motivating playground)"

related:
  - "[[../../audyt/L08_kink_fermion_closure/README.md]] problem #3"
  - "[[../../meta/CYCLE_KICKOFF_TEMPLATE.md]] §1-2 BINDING contract"
  - "[[../../meta/CALIBRATION_PROTOCOL.md]] 8/8 gate"
  - "[[../../meta/PRE_REGISTERED_FALSIFIERS.md]] PR-015 candidate"
  - "[[../../core/formalizm/dodatekO_u1_formalizacja.tex]] U(1) winding quantization"
  - "[[../../core/formalizm/dodatekE_pi1_formal.tex]] kink quantization"

classification: DERIVATION — L08 audit problem #3 quark confinement structural
priority: high (P2 OPEN klaster D ontology; addresses quark sub-component of #3 after mass formula HALT-B)
goal: "Wyprowadzić structural mechanism confinement kwarków z compact U(1) J_phase winding quantization (dodatekO thm:winding_quant) + SM fractional quark charges. Show that hadron composition rule N_q - N_q̄ ≡ 0 (mod 3) jest WYMUSZONA strukturalnie przez integer-winding requirement. Test rule against 14+ observed PDG hadrons + 4+ forbidden configurations. Decision A-/B+ per match rate; pre-registered threshold 100%."
estimated_effort: "~1 sesja (Phase 0 + Phase 1 sympy + Phase FINAL)"
target_window: "Phase 1: 13-15 sympy sub-tests (T1 winding theorem, T2 quark assignments, T3 rule derivation, T4-T9 hadron classifications, T10-T11 LIT exotic states, T12 general theorem, T13 S05 DEC). Goal: A- if 14/14 PASS + structural derivation clean; B+ if partial."

six_requirements_target:
  - "P1: Winding quantization theorem symbolic verification (T1 FP, from dodatekO)"
  - "P2: Quark fractional winding assignment from SM charges (T2 FP)"
  - "P3: Composition rule derivation: n_total ∈ ℤ ⟺ N_q-N_q̄ ≡ 0 mod 3 (T3 FP, T12 general)"
  - "P4: 14+ PDG hadron classifications correct (T4-T8 FP central tests)"
  - "P5: 4+ forbidden configurations correctly excluded (T6 FP)"
  - "P6: S05 single-Φ preserved (T13 DEC)"

risk_flags:
  - "R1: 1/3 fractional charge for quarks taken from SM, NOT derived from TGP foundations — partial closure flag (B+ if not addressed)"
  - "R2: Quantitative confinement σ ≈ 1 GeV/fm NOT derived (mapping, not derivation) — out of scope; topological mechanism is qualitatively different from QCD energetic mechanism"
  - "R3: Exotic states (pentaquark, tetraquark, dibaryon) — structural composition rule allowed, but BINDING/STABILITY is separate question requiring dynamics"
  - "R4: 'Composition rule' jest topology-level; phenomenological stability requires kink-kink interaction analysis (Yukawa + Berry connections) — separate cycle"

phase_plan:
  Phase_0: "Balance sheet + 8/8 ☑ gate; scope: topological confinement mechanism only (not quantitative σ)"
  Phase_1: "Sympy T1-T13: winding theorem, quark assignments, rule derivation, hadron classifications, forbidden configs, exotic states LIT, S05 DEC"
  Phase_FINAL: "Honest verdict A-/B+ per data; L08 problem #3 quark partial closure; downstream impact"

tags:
  - L08
  - hadron-topology
  - confinement
  - compact-U(1)
  - winding-quantization
  - composition-rule
  - audit-closure-candidate
  - quark-component-of-problem-3
  - cycle-scaffold-2026-05-16
---

# op-L08-Phase6-hadron-topology-confinement-2026-05-16

> **Cel:** Wyprowadzić strukturalny mechanism confinement kwarków z compact U(1)
> J_phase winding quantization (dodatekO thm:winding_quant) + SM fractional quark
> charges. Pokazać że hadron composition rule **N_q - N_q̄ ≡ 0 (mod 3)** jest
> WYMUSZONA strukturalnie przez integer-winding requirement (n_total ∈ ℤ).
> Test reguły vs 14+ obserwowanych PDG hadronów + 4+ forbidden konfiguracji.
> **Decision A-/B+ per dane Phase 1.**

## §0 — Cel + native-first contract

[CITE: `meta/CYCLE_KICKOFF_TEMPLATE.md` §1; `dodatekO_u1_formalizacja.tex` thm:winding_quant]

### §0.1 — Native observable target

**Co fizycznie liczymy:**

- Klasyfikacja każdej kombinacji N quark + M antiquark jako:
  - **ALLOWED** (n_total = Σn_i ∈ ℤ → izolowalna cząstka w compact U(1))
  - **FORBIDDEN** (n_total ∉ ℤ → strukturalnie zakazana)
- Test dla **14+ obserwowanych PDG hadronów** (baryons + mesons + exotic)
- Test dla **4+ FORBIDDEN konfiguracji** (izolowany kwark, di-kwark, 4-kwark, 5-kwark)
- Derivation reguły ogólnej: **N_q - N_q̄ ≡ 0 (mod 3)** z compactness θ ∈ [0, 2π)

**Instrument:** PDG 2024 baryon/meson compositions; LHCb 2015 pentaquark; BESIII tetraquark.

### §0.2 — Pre-registered falsification rule

**Decision rule WRITTEN BEFORE any calculation (2026-05-16):**

> Jeśli compact U(1) winding quantization (dodatekO thm:winding_quant) + SM
> quark fractional winding (n_u-type = +2/3, n_d-type = -1/3) NIE klasyfikuje
> **correctly wszystkich ≥14 obserwowanych PDG hadronów** jako n_total ∈ ℤ
> (allowed) **AND** wszystkich ≥4 forbidden konfiguracji (izolowany kwark,
> di-kwark, etc.) jako n_total ∉ ℤ (forbidden) — strukturalny mechanism
> **INSUFFICIENT** → wymaga substrate extension lub revisiting quark
> winding assignment.

```
pre_registration_date: 2026-05-16
recovery_scope:
  allowed_directions:
    - "SM electroweak fractional quark charges as INPUT (not derived this cycle)"
    - "Composite-state allowance dla n_total ∈ ℤ sum constraint"
  forbidden_directions:
    - "Post-hoc tuning quark winding assignments per hadron (Lakatos)"
    - "Substrate extension SU(N>1) → S05 violation (FOUNDATIONS §1)"
    - "Multi-Φ-field substrate → S05 violation"
  if_recovery_exhausted:
    - "H1c: Compact U(1) insufficient; substrate extension required (forbidden)"
    - "H2c: Quark fractional charges incorrect (revisit electroweak inputs)"
    - "H3c: Topological mechanism partial; quantitative dynamics needed (σ derivation)"
```

### §0.3 — TGP-native check (mandatory, pre-Phase-1)

- [x] **Q1 (Pattern coverage):** Reviewed `meta/TGP_NATIVE_COMPUTATIONAL_PATTERNS.md`. Topology + compactness arguments są TGP-native (S05 + dodatekO)
- [x] **Q2 (Red flags):** NONE — no BD-form, no m_Φ usage, no fixed scalar mass
- [x] **Q3 (Inherited LOCKs):** dodatekO thm:winding_quant; lambda1 J_amp/J_phase split
- [x] **Q4 (Standard-physics tools):** Topology (compact U(1)) jest universal; native-relevance: compact substrate Φ = |Φ|exp(iθ) z S05
- [x] **Q5 (m_Φ usage):** N/A — pracujemy z winding number topology
- [x] **Q6 (GR limit framing):** N/A — flat space topology
- [x] **Q7 (ASK-RULE self-check):** Form-meaning explicit: kwarki w TGP są kinki Φ z fractional J_phase winding; confinement jest TOPOLOGICAL (cannot create isolated fractional winding in compact U(1)), NIE energetic (jak QCD gluon flux). Te dwie ramy są phenomenologically equivalent ale strukturalnie różne. ASK-RULE Trigger A documented.
- [x] **Q8 (BD-drift audit plan):** Manual self-audit in Phase FINAL (low risk — no Φ propagator, topology only)

### §0.4 — Pre-flight methodology read confirmation

**BINDING per `meta/CYCLE_KICKOFF_TEMPLATE.md` §2.6:**

- [x] Przeczytano [[../../meta/PPN_AS_PROJECTION.md]] §3.1 (inherited)
- [x] Przeczytano [[../../meta/TGP_NATIVE_COMPUTATIONAL_PATTERNS.md]] §1-§4 (inherited)
- [x] Przeczytano [[../../meta/CYCLE_KICKOFF_TEMPLATE.md]] §1-§2
- [x] Przeczytano [[../../meta/CALIBRATION_PROTOCOL.md]]
- [x] Przeczytano [[../../audyt/L08_kink_fermion_closure/README.md]] problem #3
- [x] Przeczytano [[../op-lambda1-e2-amplitude-emergence/phase1L5_amplitude_phase_separation.md]] J_amp/J_phase
- [x] Przeczytano [[../exploration_neutrino_g0_2026-05-16/topology_playground.py]] motivating exploration

**Sign-off:** Claudian (theoretical physics agent) @ 2026-05-16 sesja R-topology

### §0.5 — Sympy substance plan

- [x] Każdy test sympy ma **explicit pytanie fizyczne**
- [x] **≥75% testów** to non-trivial symbolic manipulation
- [x] **0 hardcoded T_pass=True** (Phase 6 ABSOLUTE BINDING)
- [x] PASS/FAIL determined empirically z explicit criteria

**Plan testów Phase 1 (13 tests):**

| Test | Klasa | Pytanie fizyczne |
|---|---|---|
| T1 | **FIRST_PRINCIPLES** | Winding quantization theorem: n[γ] = (1/2π)∮dθ ∈ ℤ z kompaktowości θ ∈ [0, 2π) |
| T2 | **FIRST_PRINCIPLES** | Quark winding assignment: n_q = q_q/e_0; u-type → +2/3, d-type → -1/3 |
| T3 | **FIRST_PRINCIPLES** | Composition rule derivation: n_total ∈ ℤ ⟺ Σ n_i = Σ(±1/3 lub ±2/3) ∈ ℤ ⟺ N_q - N_q̄ ≡ 0 (mod 3) |
| T4 | **FIRST_PRINCIPLES** | Baryon classification (8 baryons): p, n, Δ⁺⁺, Δ⁻, Λ, Σ⁺, Σ⁻, Ξ⁰ — wszystkie integer? |
| T5 | **FIRST_PRINCIPLES** | Meson classification (6 mesons): π⁺, π⁻, π⁰, K⁺, K⁰, J/ψ — wszystkie integer? |
| T6 | **FIRST_PRINCIPLES** | Forbidden configurations: isolated u, isolated d, ud diquark, uuud 4-quark, uuudd 5-quark — wszystkie NON-integer? |
| T7 | **FIRST_PRINCIPLES** | Pentaquark allowed: 4q+1q̄ (N-M=3) — composition uudcc̄ |
| T8 | **FIRST_PRINCIPLES** | Tetraquark cases: 2q+2q̄ allowed, 3q+1q̄ forbidden, 4q+0q̄ forbidden |
| T9 | **FIRST_PRINCIPLES** | Dibaryon allowed: 6q+0q̄ (N-M=6 ≡ 0 mod 3) — H-dibaryon uuddss |
| T10 | **LITERATURE_ANCHORED** | LHCb 2015 P_c(4380), P_c(4450) — uudcc̄ composition; consistency check |
| T11 | **LITERATURE_ANCHORED** | BESIII Z_c(3900), LHCb T_cc(3875) — tetraquark observations; consistency |
| T12 | **FIRST_PRINCIPLES** | General theorem: for composite of arbitrary u-type/d-type quarks: integer n_total ⟺ N_q - N_q̄ ≡ 0 (mod 3) symbolically |
| T13 | **DECLARATIVE** | S05 single-Φ preservation: compact U(1) θ jest phase pojedynczego complex Φ; brak multi-field substrate |

**Target:** 12-13 sympy PASS (T1-T12); 1 DEC (T13). **Honest PASS count determined by data.**

**Substance ratio:** 10 FP (76.9%) + 2 LIT (15.4%) + 1 DEC (7.7%). Above 75% FP threshold.

---

## §1 — Phase 0: balance sheet + 8/8 gate

[Patrz `Phase0_balance.md`]

## §2 — Phase 1: structural derivation + classification

[Patrz `Phase1_sympy.py` + `Phase1_sympy.txt` + `Phase1_results.md`]

## §FINAL — Verdict + L08 audit problem #3 quark sub-closure

[Patrz `Phase_FINAL_close.md`]

---

## Status

🟡 **PARKING — scaffold opened 2026-05-16 sesja R-topology**. Pre-flight methodology read confirmation: **complete** (§0.4).

**Cycle position:** 14. cykl sesji 2026-05-16 (post-Q-quark HALT-B; quark sub-problem #3 approach #2: topology zamiast mass-formula).

Phase 0 commit gate:
1. README.md (this file) z BINDING contract — **DONE**
2. PR-015 candidate documented w Phase_FINAL (pending verdict)

This session deliverables target:
- README.md — **DONE**
- Phase0_balance.md — **PLANNED**
- Phase1_sympy.py — **PLANNED** (10 FP + 2 LIT + 1 DEC; 0 hardcoded)
- Phase1_sympy.txt — **PLANNED**
- Phase1_results.md — **PLANNED**
- Phase_FINAL_close.md — **PLANNED**

---

**Cycle scaffolded:** 2026-05-16 (Claudian, theoretical physics expert role; opens after 13 sister cycles closed 2026-05-16. Addresses L08 problem #3 quark sub-component via TOPOLOGICAL mechanism — different approach than HALT-B'ed mass-formula cycle, motivated by exploration_neutrino_g0_2026-05-16 playground).

**Cross-references:**
- [[../../meta/CYCLE_KICKOFF_TEMPLATE.md]] §1-§2 BINDING contract
- [[../../meta/CALIBRATION_PROTOCOL.md]] 8/8 gate
- [[../../audyt/L08_kink_fermion_closure/README.md]] problem #3 source
- [[../op-L08-Phase6-quark-sector-mass-formula-2026-05-16/Phase_FINAL_close.md]] HALT-B predecessor (mass-formula approach)
- [[../exploration_neutrino_g0_2026-05-16/topology_playground.py]] motivating exploration
- [[../../core/formalizm/dodatekO_u1_formalizacja.tex]] thm:winding_quant LIVE source
