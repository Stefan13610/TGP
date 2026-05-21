---
title: "op-neutrino-red-giant-tension-analysis — does TGP μ_ν prediction create real astrophysical tension?"
date: 2026-05-17
type: research-cycle
folder_status: active
parent: "[[../../TGP_FOUNDATIONS.md]]"

contract:
  L1_native:
    output_observable: "Sensitivity + uncertainty analysis dla TGP prediction μ_ν^TGP ≈ 3.5·10⁻¹² μ_B (z cycle 3 L_kink bracketing) vs red-giant astrophysical bound (Capozzi-Raffelt 2020: μ_ν < 3·10⁻¹² μ_B at 95% CL). Konkretnie: (a) TGP prediction uncertainty z propagacji m_X anchor uncertainty (L06: factor 1.7 z target), (b) suppression factor sensitivity (heurystyczny placeholder vs alternative loop forms), (c) statistical assessment tension level (1σ/2σ/3σ)."
    measurement_instrument: "Sympy numerical propagation z calibrated inputs; sensitivity scan m_X w [60, 150] MeV range; suppression factor scan (L_kink/λ_C)^n dla n ∈ [1,3]"
    native_coefs_constrained:
      - "TGP prediction μ_ν^TGP_central = 3.55·10⁻¹² μ_B (cycle 3 spinor B)"
      - "m_X uncertainty ~factor 1.7 (L06 anchor → derived target 100 MeV)"
      - "Suppression factor power n: heurystyczny n=2 (placeholder); n=1 or n=3 alternative bounds"
    falsification_rule: "TENSION REAL jeśli TGP_central > Capozzi-Raffelt bound z statistical significance >2σ, considering joint uncertainty. NO TENSION jeśli TGP_central w 1σ range of Capozzi-Raffelt OR uncertainty propagation gives confidence interval overlapping bound. MARGINAL jeśli 1σ < tension < 2σ (within scope of model uncertainty)."
    pre_registration_date: "2026-05-17"

  L2_framework_reduction:
    target_frameworks:
      - "Standard stellar evolution z plasmon decay → ν+ν̄"
      - "Tip-of-the-red-giant-branch (TRGB) brightness as tracer"
      - "Solar neutrino magnetic moment constraints"
    reduction_type: "consistency-check"
    failure_disposition: "L1-stands"

  L3_falsification_map:
    - { bound: "Capozzi-Raffelt 2020 μ_ν < 1.2·10⁻¹² μ_B (TRGB best-fit 2σ)", constrains: "TGP μ_ν^TGP", window: "experimental astrophysics", status: "pending Phase 1" }
    - { bound: "Raffelt 1990 classical bound μ_ν < 3·10⁻¹² μ_B", constrains: "older but conservative", window: "experimental", status: "pending Phase 1" }
    - { bound: "Viaux+2013 globular cluster M5 μ_ν < 4.5·10⁻¹² μ_B", constrains: "M5 bright RGB", window: "experimental", status: "pending Phase 1" }
    - { bound: "XENONnT 2022 μ_ν < 6.3·10⁻¹² μ_B", constrains: "lab", window: "experimental", status: "PASS z cycle 3" }

tgp_status:
  level: L1
  kind: tension-analysis
  output_type: sensitivity-assessment
  core_compatibility: review-only
  may_edit_core: false
  open_bridges:
    - "Cycle 3 L_kink bracketing prediction validation"
    - "m_X anchor sensitivity (L06 NUMERICAL ANCHOR status)"

predecessors:
  - "[[../op-neutrino-L_kink-bracketing-2026-05-17/]] (cycle 3, μ_ν prediction source)"
  - "[[../op-neutrino-RP2-wake-extension-2026-05-17/]] (cycle 2, spinor channel)"
  - "[[../op-neutrino-omega-motion-wake-2026-05-17/]] (cycle 1, β-task)"
  - "[[../op-L06-axion-mass-derivation-2026-05-16/]] (m_X anchor source)"

classification: TENSION-ANALYSIS (NIE first-principles derivation) — sensitivity + uncertainty propagation
priority: medium-high (sesja 2026-05-17 cycle 3 wykazał early tension warning)
goal: "Określić czy TGP μ_ν^TGP ≈ 3.5·10⁻¹² μ_B jest w **realnej** tension z astrofizycznymi bounds dla red giants (Capozzi-Raffelt 2020 best), czy mieści się w joint uncertainty z m_X anchor + heuristic suppression. Decision: TENSION REAL (>2σ) → TGP revision needed; MARGINAL (1-2σ) → flag dla follow-up; NO TENSION (<1σ) → TGP consistent z bounds."
estimated_effort: "~1.5h (compact sensitivity cycle)"

six_requirements_target:
  - "P1: Best red-giant bound extracted z literature (T1 LIT)"
  - "P2: TGP prediction central + 1σ/2σ range (T3 FP)"
  - "P3: m_X anchor sensitivity scan (T4 FP)"
  - "P4: Suppression factor sensitivity (T5 FP)"
  - "P5: Tension statistical assessment (T6 FP)"
  - "P6: Falsifiability re-assessment post-tension (T7 FP)"

risk_flags:
  - "R1: Astrophysical bound itself ma systematic uncertainties (stellar model + opacity tables); 2σ range może być >2× nominal"
  - "R2: Heuristic suppression (L_kink/λ_C_ν)² jest placeholder; alternative forms (e.g., G_F·m_ν²) mogą dawać OOM different"
  - "R3: m_X anchor z L06 jest factor 1.7 z target 100 MeV — sensitivity scan needed"
  - "R4: Cycle 4 tension verdict zależy od interpretation choices — honest CI propagation required"

phase_plan:
  Phase_0: "Balance + 8/8 gate; scope: sensitivity analysis, NIE re-derivation of mechanism"
  Phase_1: "Sympy T1-T8: bound extraction, uncertainty propagation, sensitivity scans, tension assessment"
  Phase_FINAL: "Verdict tension status + downstream implications"

tags:
  - neutrino
  - red-giant
  - astrophysics
  - tension-analysis
  - sensitivity
  - mu-nu
  - sesja-2026-05-17-cycle-4
  - empirical-fit
---

# op-neutrino-red-giant-tension-analysis-2026-05-17

> **Cel:** Określić strukturalnie czy TGP prediction μ_ν^TGP ≈ 3.5·10⁻¹² μ_B
> (cycle 3) jest w **realnej** tension z red-giant astrofizycznymi bounds,
> czy mieści się w joint uncertainty propagacji.

## §0 — Cel + native-first contract

### §0.1 — Native observable target

**Pytanie kluczowe:** Czy TGP prediction tworzy **falsifying tension** z best
red-giant data, czy mieści się w uncertainty?

**Inputs:**
- TGP central: μ_ν^TGP_central = 3.55·10⁻¹² μ_B (cycle 3 spinor B scenario)
- TGP uncertainty: propagated z (a) m_X anchor [60, 150] MeV range; (b) suppression power n ∈ [1, 3]
- Red-giant bounds:
  - Capozzi-Raffelt 2020: μ_ν < 1.2·10⁻¹² μ_B (TRGB best-fit, 2σ)
  - Raffelt 1990 classical: μ_ν < 3·10⁻¹² μ_B (conservative)
  - Viaux+2013 (M5 globular): μ_ν < 4.5·10⁻¹² μ_B (95% CL)

**Output:**
- Tension level: REAL (>2σ) / MARGINAL (1-2σ) / NONE (<1σ)
- TGP prediction CI propagated (joint uncertainty z m_X + n)
- Empirical commitment update post-tension

### §0.2 — Pre-registered falsification rule

```
TENSION_REAL (>2σ above bound):
  → TGP scenario B + heuristic n=2 ruled out at 2σ
  → Revision needed: alternative L_kink scenario OR rigorous suppression form
  → Possibly downgrade cycle 3 from B+ constraining → HALT-B requires revision

TENSION_MARGINAL (1-2σ):
  → Within model uncertainty; honestly flagged
  → Falsifiability strengthened — next-gen bound tightening will resolve

NO_TENSION (<1σ):
  → TGP prediction passes red-giant test ✓
  → Cycle 3 prediction CONFIRMED structurally
```

```
pre_registration_date: 2026-05-17
recovery_scope:
  allowed:
    - "Uncertainty propagation z established inputs (m_X anchor uncertainty, suppression power n range)"
    - "Statistical tension assessment z standard methodology"
    - "Bound comparison z multiple red-giant sources"
  forbidden:
    - "Post-hoc tuning suppression form aby dopasować TGP poniżej bound"
    - "Hardcoded T_pass=True"
    - "Cherry-picking weakest bound"
```

### §0.3 — TGP-native check

- [x] Q1: Sensitivity analysis jest standard practice w physics
- [x] Q2: No m_Φ usage
- [x] Q3: Inherited z cycle 3 (μ_ν central; m_X anchor)
- [x] Q4: Standard error propagation
- [x] Q5: N/A
- [x] Q6: N/A
- [x] Q7: ASK-RULE — "tension" = statistical disagreement z bound, NOT theoretical issue
- [x] Q8: Manual audit Phase FINAL

### §0.4 — Pre-flight read confirmation

- [x] [[../op-neutrino-L_kink-bracketing-2026-05-17/Phase_FINAL_close.md]] §3.1 (TGP position)
- [x] [[../op-neutrino-L_kink-bracketing-2026-05-17/Phase_FINAL_close.md]] §3.2 (falsifiability)
- [x] Cycle 3 Phase1_sympy.py — μ_ν computation
- [x] [[../../meta/CYCLE_KICKOFF_TEMPLATE.md]], [[../../meta/CALIBRATION_PROTOCOL.md]]

### §0.5 — Sympy substance plan

| Test | Klasa | Pytanie |
|---|---|---|
| T1 | **LIT** | Best red-giant bound — Capozzi-Raffelt 2020 μ_ν < 1.2·10⁻¹² μ_B (2σ TRGB) |
| T2 | **FP** | Statistical interpretation: 1σ ≈ 0.6·10⁻¹², 2σ ≈ 1.2·10⁻¹² |
| T3 | **FP** | TGP central + uncertainty propagation z m_X scan |
| T4 | **FP** | m_X anchor sensitivity [60, 150] MeV |
| T5 | **FP** | Suppression power n ∈ [1, 3] sensitivity |
| T6 | **FP** | Tension level: σ_tension = (TGP - bound) / σ_combined |
| T7 | **FP** | Falsifiability re-assessment: next-gen sensitivity required |
| T8 | **DEC** | S05 preservation; no new free parameters |

**Substance ratio:** 6 FP + 1 LIT + 1 DEC = 75% FP ✓. Hardcoded: 0.

## §1 — Status

🟢 **ACTIVE — opened 2026-05-17 (4th cycle sesji)**

## §2 — Cross-references

- [[../op-neutrino-L_kink-bracketing-2026-05-17/]] (predecessor, prediction source)
- [[../op-neutrino-RP2-wake-extension-2026-05-17/]] (spinor mechanism)
- [[../op-neutrino-omega-motion-wake-2026-05-17/]] (β-task)
- [[../op-L06-axion-mass-derivation-2026-05-16/]] (m_X anchor)
- [[../../audyt/NUMERICAL_ANCHORS_REGISTRY.md]] (anchor framework)

---

**Scaffolded:** 2026-05-17 Claudian (4th cycle sesji β-task-resolution).
