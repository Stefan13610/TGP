---
title: "op-L01-N5-retrofit-native-EW — first-principles SU(2)×U(1) gauge anomaly + EW phase transition z TGP framework"
date: 2026-05-13
type: research-cycle
folder_status: parking
parent: "[[../../TGP_FOUNDATIONS.md]]"

contract:
  L1_native:
    output_observable: "EW phase transition critical temperature T_EW ~ 159 GeV; sphaleron rate suppression dla baryogenesis BAU; gauge boson trace contributions ⟨T^μ_μ_EW⟩ [GeV⁴]"
    measurement_instrument: "Lattice EW phase transition simulations (Kajantie+1996; D'Onofrio+2014); Planck 2018 BAU η_B = 6.13·10⁻¹⁰; future precision EWPO at FCC-ee"
    native_coefs_constrained:
      - "β_SU(2) 1-loop coefficient b_2 = 19/6 dla SM EW (NIE general SU(2))"
      - "β_U(1) hypercharge b_1 = -41/10"
      - "Φ_eq(t_EW) profile at T = T_EW ~ 159 GeV (z ~ 10¹⁵)"
    falsification_rule: "Jeśli FCC-ee precision EWPO (S, T, U parameters lub W boson mass) measurement shows |Δ_TGP/Δ_SM| > 0.1% z 5σ confidence beyond SM expectations, TGP EW gauge sector w g_eff[{Φ_i}] background insufficient → wymaga (a) direct Φ-W±/Z portal vertex (S05 challenged) lub (b) revised emergent-metric Phase 1 ansatz {A,B,C} dla EW symmetry breaking."
    pre_registration_date: "2026-05-13"

  L2_framework_reduction:
    target_frameworks:
      - "SM-EW-1-loop-curved-background"
      - "EW phase transition Friedmann radiation-era"
      - "Sphaleron rate calculation lattice"
    reduction_type: "not-attempted"
    failure_disposition: "L1-stands"

  L3_falsification_map:
    - { bound: "Planck 2018 BAU η_B = 6.13·10⁻¹⁰", constrains: "Sphaleron rate EW epoch", window: "1σ consistent z SM crossover scenario", status: "inherited PASS" }
    - { bound: "PDG 2024 sin²θ_W(M_Z) = 0.23121 ± 0.00004", constrains: "g, g' running consistency", window: "structural", status: "structural" }
    - { bound: "ATLAS+CMS m_W = 80.369 ± 0.013 GeV", constrains: "EW symmetry breaking scale", window: "1σ consistent", status: "inherited PASS" }

tgp_status:
  level: T2
  kind: retrofit
  output_type: observable
  core_compatibility: review-only
  may_edit_core: false
  has_needs_file: false
  has_findings_file: false
  exports_findings: false
  open_bridges:
    - "Sphaleron rate non-perturbative calculation (lattice deferred)"
    - "BAU baryogenesis explicit mechanism (beyond scope)"
  depends_on:
    - "op-L01-rho-stress-energy-bridge-2026-05-04"
    - "op-emergent-metric-from-interaction-2026-05-09"
    - "op-L01-N4-retrofit-native-Higgs-2026-05-13 (sister)"
  impacts:
    - "L01 NEEDS §N5 retrofit replaces C LITERATURE_ANCHORED predecessor"
  source_of_status:
    - "RESEARCH_RESTART_2026-05-11 §1.3 N5-EW candidate"

predecessors:
  - "[[../op-L01-N5-EW-gauge-anomaly-2026-05-11/]] (CLOSED-DOWNGRADED C — LITERATURE_ANCHORED)"

related:
  - "[[../op-L01-N1-retrofit-native-EM-2026-05-13/]]"
  - "[[../op-L01-N2-retrofit-native-QCD-2026-05-13/]]"
  - "[[../op-L01-N4-retrofit-native-Higgs-2026-05-13/]]"

classification: RETROFIT — C → A− target
goal: "First-principles SU(2)×U(1) gauge anomaly mechanism w TGP framework: β_SU(2), β_U(1) running w g_eff[{Φ_i}]; T_EW phase transition consistent z Φ_eq(t_EW); BAU consistency."
estimated_effort: "~1.5 sesji"
target_window: "Phase 1: β_2 = 19/6·g³/(16π²), β_1 = -41/10·g'³/(16π²); sin²θ_W(M_Z) = 0.23121; T_EW ~ 159 GeV; m_W = 80.369 GeV; sphaleron rate structural"

six_requirements_target:
  - "P1: β_SU(2) i β_U(1) 1-loop coefficients explicit symbolic"
  - "P2: Weinberg angle sin²θ_W = g'²/(g²+g'²) z explicit derivation"
  - "P3: EW phase transition T_EW ~ 159 GeV consistent z V_eff(H, T)"
  - "P4: Planck BAU η_B preserved (sphaleron freeze-out)"
  - "P5: Φ_eq(t_EW) cosmology constraint"
  - "P6: S05 single-Φ; no direct Φ-W/Z portal"

risk_flags:
  - "R1: Non-perturbative sphaleron rate wymaga lattice (only structural here)"
  - "R2: BAU explicit mechanism (CP-violation source) beyond TGP scope"
  - "R3: EW phase transition crossover (NIE first-order) — SM preserved by TGP"

phase_plan:
  Phase_0: "Balance + scope"
  Phase_1: "First-principles β_2, β_1, sin²θ_W, T_EW + Φ_eq EW epoch"
  Phase_FINAL: "Closure + L3 falsification map"

tags:
  - L01
  - L01-N5-retrofit
  - EW-gauge-anomaly
  - SU2-U1-running
  - sphaleron-BAU
  - retrofit-C-to-A
  - cycle-scaffold-2026-05-13
---

# op-L01-N5-retrofit-native-EW-2026-05-13

> **Cel:** first-principles SU(2)×U(1) gauge anomaly + EW phase transition.
> Replace C LITERATURE_ANCHORED predecessor z ≥75% FP.

## §0 — Cel + contract

### §0.1 — Native observable

- T_EW critical temperature ~159 GeV
- BAU η_B = 6.13·10⁻¹⁰ consistency
- Weinberg angle sin²θ_W(M_Z) = 0.23121

### §0.2 — Pre-registered rule

```
pre_registration_date: 2026-05-13
recovery_scope:
  allowed: ["β-function precision refinement (g, g')", "epoch Φ_eq calibration"]
  forbidden: ["direct Φ-W/Z portal (S05)", "BAU post-hoc tuning"]
```

### §0.3 — Q1-Q8 OK

### §0.4 — Pre-flight read

- [x] PPN_AS_PROJECTION §3.1
- [x] TGP_NATIVE_COMPUTATIONAL_PATTERNS §1-§4
- [x] M9_RESTRUCTURE_NOTE §1.4 + §3
- [x] CYCLE_KICKOFF_TEMPLATE §1-§2

**Sign-off:** Claudian @ 2026-05-13

### §0.5 — Sympy plan

| Test | Klasa | Pytanie |
|---|---|---|
| T1 | FP | β_SU(2)^SM b_2 = 19/6 i β_U(1)^SM b_1 = -41/10 (sign + value symbolic) |
| T2 | FP | Weinberg angle sin²θ_W = g'²/(g²+g'²) (geometric definition) |
| T3 | FP | M_W² = g² v²/4; M_Z² = (g² + g'²) v²/4 (gauge boson masses z VEV) |
| T4 | FP | β-function running: g and g' evolve oppositely — g shrinks z μ, g' grows |
| T5 | FP | Sphaleron rate ~ exp(-E_sph/T) z E_sph ~ 4π v/g — exponential suppression structural |
| T6 | FP | g_eff[{Φ_i}] Pattern 2.1 → EW sektor couples to Φ via metric (NIE direct vertex) |
| T7 | LIT | sin²θ_W(M_Z) = 0.23121 ± 0.00004 (PDG 2024) |
| T8 | LIT | m_W = 80.369 ± 0.013 GeV (PDG 2024) + η_B = 6.13·10⁻¹⁰ (Planck) |
| T9 | DECLARATIVE | S05 + scope |

**Target:** 6 FP + 2 LIT + 1 DEC.

## Status

🟡 **PARKING — scaffold 2026-05-13**.
