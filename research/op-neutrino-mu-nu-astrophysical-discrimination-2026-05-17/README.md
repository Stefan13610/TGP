---
title: "op-neutrino-mu-nu-astrophysical-discrimination — survey wszystkich astrofizycznych bounds dla dyskryminacji μ_ν^TGP scenarios A vs B"
date: 2026-05-17
type: research-cycle
folder_status: active
parent: "[[../../TGP_FOUNDATIONS.md]]"

contract:
  L1_native:
    output_observable: "Joint statistical discrimination scenario A (m_X-scale, μ_ν^TGP_A = 3.55·10⁻¹² μ_B) vs scenario B (SM-like Lee-Shrock, μ_ν^TGP_B = 3.2·10⁻²⁰ μ_B) z comprehensive survey wszystkich istniejących astrofizycznych bounds dla μ_ν. Per-bound tension level w σ jednostkach przy użyciu joint uncertainty propagation (TGP CI z m_X anchor + bound systematic). Konkretnie: czy któryś bound (TRGB, SN1987A, BBN, Solar RSFP, BH accretion, globular clusters) odrzuca scenario A z >2σ, czy wszystkie pozostają compatible przy honest CI?"
    measurement_instrument: "Astrophysical observations: TRGB photometry (Hubble + ground), SN1987A neutrino burst (Kamiokande-II + IMB + Baksan), CMB/light-element BBN inference (Planck + LUNA), solar neutrino spectra (Super-K + SNO + Borexino), BH accretion disk emission (X-ray observatories), globular cluster RGB-tip photometry (HST + JWST). Comparison: sympy log-space joint σ_tension."
    native_coefs_constrained:
      - "Scenario A: μ_ν^TGP = 3.55·10⁻¹² μ_B z m_X uncertainty CI [1.28, 3.55]·10⁻¹² μ_B (cycle 3 + cycle 4 propagation)"
      - "Scenario B: μ_ν^TGP = 3.2·10⁻²⁰ μ_B z v_H scale (cycle 6 Lee-Shrock)"
      - "Per-bound systematic log-σ (stellar opacity, plasmon emissivity, BBN inputs, RSFP magnetic field assumption)"
    falsification_rule: "A- DISCRIMINATION: jakikolwiek bound wyklucza scenario A z joint σ_tension > 2.0 → strukturalnie preferred scenario B (SM-like). A- BOTH CONSISTENT: wszystkie bounds compatible z obu scenarios przy joint CI (max σ_tension ≤ 2.0) → dual-scenario stoi; dyskryminacja wymaga XLZD/DARWIN ~2030+. B+ PARTIAL: mieszane wyniki — niektóre bounds problematyczne (1 < σ_tension ≤ 2), nie wszystkie discriminative. HALT-B: problemy methodology lub interpretation w odczycie bounds."
    pre_registration_date: "2026-05-17"

  L2_framework_reduction:
    target_frameworks:
      - "Standard plasmon decay γ* → νν̄ z μ_ν vertex (Raffelt 1990 formalism)"
      - "SN1987A neutrino cooling argument (Magill+2018; Lattimer-Cooperstein 1988)"
      - "BBN N_eff shift z magnetic moment-induced thermalization (Cyburt+2016)"
      - "Akhmedov 1988 RSFP w solar magnetic field"
    reduction_type: "consistency-check"
    failure_disposition: "L1-stands"

  L3_falsification_map:
    - { bound: "TRGB Capozzi-Raffelt 2020 μ_ν < 1.2·10⁻¹² μ_B (2σ)", constrains: "scenarios A, B", window: "TRGB stellar evolution", status: "pending Phase 1" }
    - { bound: "SN1987A Magill+2018 μ_ν < 1.3·10⁻¹² μ_B (95% CL)", constrains: "scenarios A, B", window: "neutrino burst cooling", status: "pending Phase 1" }
    - { bound: "Arceo-Diaz+2015 (ωCen globular) μ_ν < 2.2·10⁻¹² μ_B (95% CL)", constrains: "scenarios A, B", window: "globular RGB-tip", status: "pending Phase 1" }
    - { bound: "Viaux+2013 (M5 globular) μ_ν < 4.5·10⁻¹² μ_B (95% CL)", constrains: "scenarios A, B", window: "M5 RGB-tip", status: "pending Phase 1" }
    - { bound: "BBN N_eff Cyburt+2016 μ_ν < ~10⁻¹⁰ μ_B (cosmological)", constrains: "scenarios A, B", window: "primordial nucleosynthesis", status: "pending Phase 1" }
    - { bound: "Solar RSFP non-observation Super-K/Borexino μ_ν < ~7·10⁻¹¹ μ_B", constrains: "scenarios A, B", window: "Akhmedov mechanism", status: "pending Phase 1" }
    - { bound: "BH accretion disk Latimer-Burrows 2007 μ_ν < ~10⁻¹⁰ μ_B", constrains: "scenarios A, B", window: "hot dense plasma", status: "pending Phase 1" }

tgp_status:
  level: L1
  kind: empirical-discrimination
  output_type: joint-statistical-discrimination
  core_compatibility: review-only
  may_edit_core: false
  open_bridges:
    - "Cycle 6 dual-scenario validation"
    - "Per-bound systematic uncertainty propagation"
    - "Joint-bound discrimination methodology"

predecessors:
  - "[[../op-WZ-emergence-quantitative-loop-2026-05-17/]] (cycle 6, dual-scenario source)"
  - "[[../op-neutrino-red-giant-tension-analysis-2026-05-17/]] (cycle 4, joint CI methodology to replicate)"
  - "[[../op-neutrino-L_kink-bracketing-2026-05-17/]] (cycle 3, scenario A prediction source)"

classification: EMPIRICAL-DISCRIMINATION-SURVEY (NIE first-principles derivation) — comprehensive bound survey + joint CI per-bound tension assessment
priority: high (sesja 2026-05-17 cycle 6 wykazał DUAL_SCENARIO; cycle 4 sprawdziło tylko 1 z ~6 bounds)
goal: "Określić strukturalnie czy comprehensive astrofizyczny bound survey może zdyskryminować scenario A vs scenario B beyond singleton TRGB, w zgodzie z cycle 4 joint-CI methodology. Decision tree: A- DISCRIMINATION (jeden bound wyklucza A z >2σ); A- BOTH CONSISTENT (wszystkie compatible przy joint CI); B+ PARTIAL (mieszane); HALT-B (methodology problems)."
estimated_effort: "~2-3h (comprehensive survey + 7 bounds × 2 scenarios joint CI computation)"

six_requirements_target:
  - "P1: Comprehensive bound survey z literature (T1 LIT)"
  - "P2: TRGB bound reproduction z cycle 4 (T2 FP — consistency check)"
  - "P3: SN1987A constraint per scenario (T3 FP)"
  - "P4: BBN N_eff constraint per scenario (T4 FP)"
  - "P5: Solar ν RSFP non-observation per scenario (T5 FP)"
  - "P6: BH accretion disk constraint per scenario (T6 FP)"

risk_flags:
  - "R1: Per-bound systematic uncertainties różnią się (TRGB stellar opacity ~0.3 dex; SN1987A neutrino emissivity ~0.5 dex; BBN N_eff ~0.2 dex; Solar B_⊙ ~factor 2). Wszystkie muszą być honestly incorporated"
  - "R2: Naive bound combination może overestimate discrimination — correlated systematics (plasmon physics shared TRGB+SN+globular; cycle 4 lesson)"
  - "R3: Scenario A jest na granicy najtighter bounds (TRGB 1.2·10⁻¹², SN1987A ~1.3·10⁻¹²); subjective methodology choices mogą tip verdict w obu kierunkach"
  - "R4: Scenario B jest tak głęboko pod current bounds (10⁻²⁰ vs 10⁻¹² najtighter) że trivialnie compatible — to NIE jest discrimination per se, tylko consistency"
  - "R5: BH accretion bound jest model-dependent (Latimer-Burrows assumptions); może być cherry-pickable w obu kierunkach — używać CONSERVATIVE bound"
  - "R6: Pre-registered falsifier MUST być stosowany bez post-hoc adjustment per CALIBRATION_PROTOCOL §1"

phase_plan:
  Phase_0: "Balance + 8/8 gate; scope: discrimination survey, NIE re-derivation prediction; cycle 4 methodology replication"
  Phase_1: "Sympy T1-T8: bound survey, per-bound joint CI tension dla A i B, joint discrimination verdict"
  Phase_FINAL: "Verdict per decision tree + downstream implications dla PR-016 dual-scenario status"

tags:
  - neutrino
  - mu-nu
  - astrophysics
  - discrimination
  - dual-scenario
  - bound-survey
  - sesja-2026-05-17-cycle-7
  - empirical-discrimination
---

# op-neutrino-mu-nu-astrophysical-discrimination-2026-05-17

> **Cel:** Określić strukturalnie czy comprehensive astrofizyczny bound survey
> (TRGB + SN1987A + BBN + Solar RSFP + BH accretion + globular clusters)
> dyskryminuje TGP μ_ν^TGP **scenario A (3.55·10⁻¹² μ_B)** vs **scenario B
> (3.2·10⁻²⁰ μ_B)** ustanowione w cycle 6, czy oba pozostają consistent przy
> honest joint CI methodology (cycle 4 protocol replication).

## §0 — Cel + native-first contract

### §0.1 — Native observable target

**Pytanie kluczowe:** Czy któryś tighter astrofizyczny bound (poza TRGB cycle 4)
**wyklucza scenario A z joint σ_tension > 2.0**, czy dual-scenario musi
zostać dla XLZD/DARWIN ~2030+ experimental discrimination?

**Inputs:**
- TGP scenario A central: μ_ν^TGP_A = 3.55·10⁻¹² μ_B (cycle 3 spinor B, m_X = 60 MeV)
- TGP scenario A CI: [1.28·10⁻¹², 3.55·10⁻¹²] μ_B (cycle 4 m_X uncertainty propagation 60→100 MeV)
- TGP scenario A geomean: 2.13·10⁻¹² μ_B (log-fair central)
- TGP scenario A log-σ: 0.22 dex (cycle 4 value)
- TGP scenario B central: μ_ν^TGP_B = 3.2·10⁻²⁰ μ_B (cycle 6 Lee-Shrock z v_H = 246 GeV)
- TGP scenario B log-σ: ~0.3 dex (v_H precise, m_e precise, m_ν 0.05-0.2 eV uncertainty)

**Astrophysical bounds (Phase 1 will detail):**
- TRGB Capozzi-Raffelt 2020: 1.2·10⁻¹² μ_B (TRGB 2σ), stellar log-σ ~0.30 dex
- SN1987A Magill+2018: ~1.3·10⁻¹² μ_B (95% CL), emissivity log-σ ~0.45 dex
- ωCen Arceo-Diaz+2015: 2.2·10⁻¹² μ_B (95% CL), stellar log-σ ~0.30 dex
- M5 Viaux+2013: 4.5·10⁻¹² μ_B (95% CL), stellar log-σ ~0.30 dex
- BBN N_eff: ~10⁻¹⁰ μ_B (cosmological, weaker), log-σ ~0.20 dex
- Solar RSFP: ~7·10⁻¹¹ μ_B (Super-K + Borexino 90% CL), B_⊙ log-σ ~0.30 dex
- BH accretion: ~10⁻¹⁰ μ_B (Latimer-Burrows 2007 conservative), log-σ ~0.50 dex

**Output:**
- Per-bound joint σ_tension per scenario (A, B) using cycle 4 log-space methodology
- Discrimination verdict per decision tree (§0.2)
- PR-016 dual-scenario status update (preserved or discriminated)

### §0.2 — Pre-registered falsification rule

```
A- DISCRIMINATION (jeden bound z joint σ > 2.0):
  → Scenario A excluded by tighter bound
  → SM-like scenario B preferred strukturalnie
  → PR-016 promoted to SINGLE-SCENARIO (B)
  → XLZD/DARWIN null result expected
  → Cycle 3 prediction effectively retracted (preserved jako derived-but-excluded)

A- BOTH CONSISTENT (wszystkie bounds z joint σ ≤ 2.0):
  → Dual-scenario stoi po pełnym survey
  → Empirical discrimination requires next-gen experiments (XLZD/DARWIN ~2030+)
  → PR-016 dual-scenario STRENGTHENED (passed comprehensive survey)

B+ PARTIAL (1 < max σ ≤ 2.0; mixed):
  → Niektóre bounds problematyczne ale nie definitive
  → Honest flag: scenario A "stressed" by accumulated tension
  → Dual-scenario preserved z honest stress disclosure

HALT-B (methodology/interpretation problems):
  → Bound reading ambiguous, systematics ill-defined
  → No verdict; flag dla refined-systematics future cycle
```

```
pre_registration_date: 2026-05-17
recovery_scope:
  allowed:
    - "Joint uncertainty propagation z BOTH TGP and bound (cycle 4 protocol replication)"
    - "Per-bound systematic uncertainty z literature anchors (stellar models, neutrino emissivity, BBN inputs)"
    - "Conservative bound reading (when range available: use most conservative compatible z 2σ)"
    - "Decision tree applied AS-IS without post-hoc threshold adjustment"
  forbidden:
    - "Cherry-picking weakest bound to favor scenario A"
    - "Cherry-picking tightest bound to favor scenario B"
    - "Hardcoded T_pass=True"
    - "Threshold modification po obejrzeniu verdict (per CALIBRATION_PROTOCOL §1)"
    - "Naive point-vs-point comparison without joint CI (cycle 4 anti-pattern)"
    - "Bound combination z assumed-independent likelihoods when systematics correlated"
```

### §0.3 — TGP-native check Q1-Q8

- [x] **Q1:** Discrimination survey jest standard empirical methodology — nie wymaga TGP-native re-derivation
- [x] **Q2:** No m_Φ usage; m_X anchor inherited z L06 (per cycle 3/4)
- [x] **Q3:** Inherited z cycle 3 (μ_ν^TGP_A formula) + cycle 6 (μ_ν^TGP_B formula); cycle 4 methodology (joint CI)
- [x] **Q4:** Standard log-space combined-σ z stellar/BBN/RSFP literature
- [x] **Q5:** N/A — empirical comparison cycle
- [x] **Q6:** N/A — no new TGP equations
- [x] **Q7:** ASK-RULE — "discrimination" = ≥1 bound z joint σ > 2.0 (statistical disagreement); "consistency" = all bounds joint σ ≤ 2.0; thresholds pre-registered
- [x] **Q8:** Manual audit w Phase FINAL §6 (substance breakdown)

### §0.4 — Pre-flight methodology read confirmation

- [x] [[../../meta/PPN_AS_PROJECTION.md]] §3.1 — three-layer L1/L2/L3 (this cycle: L1 native observable = joint σ per scenario; L3 = falsification map z 7 bounds)
- [x] [[../../meta/TGP_NATIVE_COMPUTATIONAL_PATTERNS.md]] §1-§4 — anti-BD-drift (N/A: empirical comparison cycle, no BD-form formulas)
- [x] [[../../meta/M9_RESTRUCTURE_NOTE.md]] §1.4 + §3 — M9.1'' inheritance (N/A: μ_ν sektor, no M9.1'' Taylor coefs)
- [x] [[../../meta/CYCLE_KICKOFF_TEMPLATE.md]] §1-§2 — kickoff contract (this README zgodny z template L1/L2/L3 sections)
- [x] [[../../meta/CALIBRATION_PROTOCOL.md]] — anti-overclaim (this is empirical-discrimination cycle, max claim A-)

**Pre-flight predecessor read:**
- [x] [[../op-WZ-emergence-quantitative-loop-2026-05-17/Phase_FINAL_close.md]] §1.2 (dual-scenario establishment) + §3 (cycle 3 dual-scenario revision)
- [x] [[../op-neutrino-red-giant-tension-analysis-2026-05-17/Phase_FINAL_close.md]] §1 (joint CI 0.67σ methodology — REPLICATE)
- [x] [[../op-neutrino-red-giant-tension-analysis-2026-05-17/Phase1_sympy.py]] T6 (log-space combined σ implementation)
- [x] [[../op-neutrino-L_kink-bracketing-2026-05-17/Phase_FINAL_close.md]] §1 (μ_ν^TGP_A scenario A details)

**Sign-off:** Claudian @ 2026-05-17 (7th cycle sesji, kontynuacja cycle 6 dual-scenario)

### §0.5 — Sympy substance plan

| Test | Klasa | Pytanie |
|---|---|---|
| **T1** | **LIT** | Comprehensive astrofizycznych bounds survey z sources (7 bounds: TRGB + SN1987A + 2 globular + BBN + Solar + BH disk) |
| **T2** | **FP** | TRGB Capozzi-Raffelt 2020 joint CI (cycle 4 reproduction; consistency check methodology) per scenario A, B |
| **T3** | **FP** | SN1987A Magill+2018 joint CI tension per scenario A, B |
| **T4** | **FP** | BBN N_eff Cyburt+2016 joint CI tension per scenario A, B |
| **T5** | **FP** | Solar ν RSFP non-observation joint CI tension per scenario A, B |
| **T6** | **FP** | BH accretion disk Latimer-Burrows 2007 joint CI tension per scenario A, B + globular ωCen/M5 supplementary check |
| **T7** | **FP** | Joint statistical assessment — discrimination verdict z pre-registered decision tree |
| **T8** | **DEC** | S05 preservation; no new free parameters; comparison-only cycle |

**Substance ratio:** 6 FP + 1 LIT + 1 DEC = 75% FP ✓. Hardcoded T_pass=True target: 0.

## §1 — Status

🟢 **ACTIVE — opened 2026-05-17 (7th cycle sesji 2026-05-17, kontynuacja po cycle 6)**

## §2 — Cross-references

### Predecessors (LIVE inheritance):
- [[../op-WZ-emergence-quantitative-loop-2026-05-17/]] (cycle 6, dual-scenario)
- [[../op-neutrino-red-giant-tension-analysis-2026-05-17/]] (cycle 4, joint CI methodology)
- [[../op-neutrino-L_kink-bracketing-2026-05-17/]] (cycle 3, scenario A formula)

### Downstream targets:
- [[../../STATE.md]] — sesja 2026-05-17 z 7 cykli (post-cycle-7 entry)
- [[../../PREDICTIONS_REGISTRY.md]] — PR-016 status post-discrimination
- [[../../audyt/L08_kink_fermion_closure/README.md]] — problem #3 neutrino sub-component post-discrimination

### Literature (T1 LIT):
- Capozzi F., Raffelt G. 2020, *Astrophysical neutrino electromagnetic properties revised*, arXiv:2007.03694 (TRGB best bound)
- Magill G., Plestid R., Pospelov M., Tsai Y.-D. 2018, *Dipole portal to heavy neutral leptons*, Phys. Rev. D 98, 115015 (SN1987A updated)
- Raffelt G.G. 1990, *Astrophysical methods*, Phys Rep 198, 1 (SN1987A classical + plasmon decay foundation)
- Arceo-Diaz S., et al. 2015, *Constraint on the magnetic dipole moment of neutrinos by the tip-RGB luminosity in ωCen*, Astropart. Phys. 70, 1 (ωCen globular)
- Viaux N. et al. 2013, *Particle-physics constraints from M5 RGB tip-luminosity*, A&A 558, A12 (M5 globular)
- Cyburt R.H., Fields B.D., Olive K.A., Yeh T.-H. 2016, *Big-bang nucleosynthesis after Planck*, Rev. Mod. Phys. 88, 015004 (BBN N_eff)
- Akhmedov E.K. 1988, *Resonance enhancement of the neutrino spin precession in solar matter*, Phys. Lett. B 213, 64 (RSFP)
- Latimer A., Burrows A. 2007, *Neutrino magnetic moments and their effects on supernova neutrino signals*, ApJ 661, 320 (BH/dense plasma)
- Borexino Collaboration 2017, *Limiting neutrino magnetic moments with Borexino Phase-II solar data*, Phys. Rev. D 96, 091103 (Solar limit)

---

**Scaffolded:** 2026-05-17 Claudian (7th cycle sesji 2026-05-17 dual-scenario discrimination).
