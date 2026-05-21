---
title: "op-neutrino-L_kink-bracketing — L_kink range estimate + μ_ν^TGP quantitative bracketing"
date: 2026-05-17
type: research-cycle
folder_status: active
parent: "[[../../TGP_FOUNDATIONS.md]]"

contract:
  L1_native:
    output_observable: "Range estimate L_kink dla TGP n=0 neutrino soliton z TGP-native scales (NIE fundamental Compton-only). Resulting quantitative μ_ν^TGP w obu channels (scalar β-task + spinor RP² extension). Conkretnie: (a) L_kink lower bound z core-scale arguments (m_X substrate); (b) L_kink upper bound z asymptotic Compton tail; (c) intermediate scenario z g_0-weighted scaling; (d) μ_ν^TGP range [μ_min, μ_max] propagated do każdego scenariusza."
    measurement_instrument: "Sympy dimensional analysis z calibrated TGP inputs: m_X ≈ 60 MeV (L06 numerical anchor), g_0_ν ≈ 0.22 (lekkie ν z why_n3 playground), A_tail_ν ≈ 6.1·10⁻⁴, m_ν = 0.1 eV PDG ν₃ scale"
    native_coefs_constrained:
      - "L_kink tail ~ ℏ/(m_ν·c) = λ_C ≈ 2 mm (lub macroscopic Yukawa decay)"
      - "L_kink core ~ ℏ/(m_X·c) ≈ 3.3 fm (z L06 substrate anchor 60 MeV)"
      - "L_kink intermediate ~ ℏ/(g_0_ν·M_eff)·A_tail factor"
      - "μ_scalar(L_kink) z β-task: δθ_wake ~ e·B·v·L_kink²/c² → Larmor-like μ"
      - "μ_spinor(L_kink) z RP² extension: μ ~ e·β·ℏ/(4m_eff·L_kink/λ_C)"
    falsification_rule: "Jeśli pełny range [μ_min, μ_max] **wszystkich** 3 scenarios całkowicie powyżej XENONnT bound (6.3·10⁻¹² μ_B) → TGP mechanism candidate RULED OUT. Jeśli pełny range całkowicie poniżej SM Dirac (3·10⁻²⁰ μ_B) → mechanism existuje strukturalnie ale empirically irrelevant. PASS scenario: w 3 scenariuszach co najmniej jeden ma μ_ν^TGP w testable window 10⁻²⁰ < μ < 10⁻¹² μ_B (falsifiable by next-gen)."
    pre_registration_date: "2026-05-17"

  L2_framework_reduction:
    target_frameworks:
      - "Standard QFT Compton wavelength λ_C = ℏ/(mc)"
      - "'t Hooft-Polyakov monopole: core radius ~ 1/(eM_W) vs tail ~ 1/M_H"
      - "TGP m_eff(ψ) = c_M·A_tail²·g_0^(e²/2) from why_n3 PHASE2"
    reduction_type: "scale-bracketing consistency"
    failure_disposition: "L1-stands (range estimate, NIE pinned-precision)"

  L3_falsification_map:
    - { bound: "XENONnT μ_ν < 6.3·10⁻¹² μ_B (2022)", constrains: "upper limit μ_ν^TGP", window: "experimental", status: "pending Phase 1 T7" }
    - { bound: "GEMMA μ_ν < 2.9·10⁻¹¹ μ_B (2012)", constrains: "consistency", window: "experimental", status: "pending T7" }
    - { bound: "Red giant cooling μ_ν < 3·10⁻¹² μ_B", constrains: "astrophysics", window: "experimental", status: "pending T7" }
    - { bound: "SM Dirac μ_ν^SM ≈ 3·10⁻²⁰ μ_B", constrains: "reference (lower)", window: "consistency", status: "pending T7" }
    - { bound: "XLZD/DARWIN target ~10⁻¹² μ_B (2030+)", constrains: "future falsifiability", window: "projected", status: "discussion" }

tgp_status:
  level: L1
  kind: quantitative-bracketing
  output_type: range-estimate
  core_compatibility: review-only
  may_edit_core: false
  open_bridges:
    - "L_kink determination strict (deferred dedicated cycle z first-principles soliton ODE)"
    - "W/Z sector quantitative loop integration (problem #3 boson — still OPEN)"

predecessors:
  - "[[../op-neutrino-omega-motion-wake-2026-05-17/]] (β-task A-)"
  - "[[../op-neutrino-RP2-wake-extension-2026-05-17/]] (β REFINED A-)"
  - "[[../op-L06-axion-mass-derivation-2026-05-16/]] (m_X ≈ 60 MeV numerical anchor)"
  - "[[../why_n3/PHASE2_n_alpha_derivation.md]] (m_eff = c_M·A_tail²·g_0^(e²/2))"
  - "[[../exploration_neutrino_g0_2026-05-16/]] (g_0_ν, A_tail_ν calibration)"

related:
  - "[[../../audyt/NUMERICAL_ANCHORS_REGISTRY.md]] (L06 anchor m_X)"
  - "[[../../meta/CYCLE_KICKOFF_TEMPLATE.md]] §1-2"

classification: QUANTITATIVE-BRACKETING (NIE first-principles derivation) — range estimate
priority: medium (continues sesja 2026-05-17 neutrino line; bracketing complements structural cycles 1+2)
goal: "Zwrócić quantitative range μ_ν^TGP per 3 scenarios L_kink, używając TGP-native inputs (m_X anchor, g_0_ν calibration). Decision: B+ partial (honest bracketing, NIE pinned precision) jeśli range testable; HALT jeśli całkowicie ruled out lub całkowicie below SM."
estimated_effort: "~1.5h (compact bracketing cycle)"

six_requirements_target:
  - "P1: L_kink^tail z Compton wavelength explicit (T1 FP)"
  - "P2: L_kink^core z m_X substrate anchor (T2 LIT)"
  - "P3: L_kink^intermediate z g_0-weighted (T3 FP)"
  - "P4: μ_scalar(L_kink) computed (T5 FP)"
  - "P5: μ_spinor(L_kink) computed (T6 FP)"
  - "P6: Falsifiability window verified (T7 FP)"

risk_flags:
  - "R1: m_X = 60 MeV jest NUMERICAL ANCHOR (L06, NIE structural derivation) — bracketing inherits anchor status, honestly documented"
  - "R2: 'Effective mass m_eff' dla neutrino w μ formuła — czy używamy m_ν (PDG) czy m_eff_TGP? Honest disclaimer: użyjemy m_ν = 0.1 eV jako observable mass"
  - "R3: Conversion δθ_wake → μ_ν wymaga loop integration (W/Z sector OPEN); używamy heuristic Larmor-like mapping"
  - "R4: Three scenarios mogą dawać overlap range — to jest OK, oznacza honest uncertainty propagation"

phase_plan:
  Phase_0: "Balance + 8/8 gate; scope: range bracketing, NOT first-principles L_kink derivation"
  Phase_1: "Sympy T1-T8: L_kink scenarios + μ_ν computations + falsifiability check"
  Phase_FINAL: "Verdict B+ partial (honest bracketing); empirical commitments"

tags:
  - neutrino
  - L-kink
  - bracketing
  - quantitative-estimate
  - mu-nu-range
  - sesja-2026-05-17-cycle-3
  - falsifiability-window
---

# op-neutrino-L_kink-bracketing-2026-05-17

> **Cel:** Range estimate L_kink dla TGP neutrino soliton z TGP-native inputs;
> quantitative μ_ν^TGP bracketing dla obu kanałów (scalar β-task + spinor RP² extension);
> falsifiability window dla testable scenarios.

## §0 — Cel + native-first contract

### §0.1 — Native observable target

**Trzy strukturalne podejścia do L_kink:**

| Approach | L_kink formula | Numerical | Interpretacja |
|---|---|---|---|
| **Compton tail** | λ_C = ℏ/(m_ν·c) | ~2 mm dla m_ν=0.1 eV | Asymptotic Yukawa decay (large r) |
| **Substrate core** | L_X = ℏc/m_X (m_X=60 MeV) | ~3.3 fm | TGP-native substrate scale (L06 anchor) |
| **g_0-weighted** | L_eff = λ_C·g_0_ν / g_0_e | ~0.5 mm | Soliton-specific calibration (interpolation) |

**μ_ν^TGP computation dla każdego scenariusza:**

Scalar channel (β-task):
$$\mu_{scalar}(L_{kink}) \sim \mu_B \cdot \alpha \cdot \frac{v}{c} \cdot \frac{L_{kink}^2 \cdot m_e c}{ℏ}$$

Spinor channel (RP² extension):
$$\mu_{spinor}(L_{kink}) \sim \mu_B \cdot \frac{e \cdot \beta}{4} \cdot \frac{\hbar}{m_{eff}} \cdot \text{(W/Z-like suppression)}$$

### §0.2 — Pre-registered falsification rule

**Decision tree post-Phase 1:**

- **B+ PASS (range testable):** co najmniej 1 scenario μ_ν^TGP w 10⁻²⁰ < μ < 10⁻¹² μ_B
- **B+ PARTIAL (range overlap):** wszystkie scenarios w testable window, ale range zbyt szeroki dla pinning prediction
- **HALT (above bound):** wszystkie scenarios > 6.3·10⁻¹² μ_B (XENONnT) — TGP mechanism RULED OUT
- **HALT (below SM):** wszystkie scenarios < 3·10⁻²⁰ μ_B (SM Dirac) — mechanism strukturalnie istnieje ale empirically dead

```
pre_registration_date: 2026-05-17
recovery_scope:
  allowed:
    - "Three L_kink scenarios bracketing (range estimate)"
    - "Heuristic δθ_wake → μ conversion (Larmor-like)"
    - "Honest classification as NUMERICAL/PARTIAL (NOT first-principles derivation)"
  forbidden:
    - "Post-hoc tuning L_kink aby fit specific empirical bound"
    - "Hardcoded μ_ν values"
    - "Claiming structural pinning of L_kink z bracketing only"
```

### §0.3 — TGP-native check

- [x] Q1: Bracketing approach z TGP-native inputs (m_X anchor, g_0_ν calibration)
- [x] Q2: No m_Φ usage (m_X ≠ m_Φ; m_X jest substrate anchor)
- [x] Q3: Inherited: β-task source S; RP² extension spinor channel; L06 m_X; why_n3 PHASE2 m_eff
- [x] Q4: 't Hooft-Polyakov core/tail analogy reference
- [x] Q5: N/A
- [x] Q6: N/A (flat space)
- [x] Q7: ASK-RULE — "bracketing" explicit (NOT precision derivation); honest disclaimer
- [x] Q8: Manual audit Phase FINAL

### §0.4 — Pre-flight read confirmation

- [x] [[../op-neutrino-omega-motion-wake-2026-05-17/Phase_FINAL_close.md]] (β-task)
- [x] [[../op-neutrino-RP2-wake-extension-2026-05-17/Phase_FINAL_close.md]] (RP² ext)
- [x] [[../op-L06-axion-mass-derivation-2026-05-16/Phase_FINAL_close.md]] (m_X anchor)
- [x] [[../../audyt/NUMERICAL_ANCHORS_REGISTRY.md]] §1 anchor #2
- [x] [[../why_n3/PHASE2_n_alpha_derivation.md]] (m_eff formula)
- [x] [[../exploration_neutrino_g0_2026-05-16/playground.py]] (g_0_ν, A_tail_ν)
- [x] [[../../meta/CYCLE_KICKOFF_TEMPLATE.md]]
- [x] [[../../meta/CALIBRATION_PROTOCOL.md]]

### §0.5 — Sympy substance plan

| Test | Klasa | Pytanie |
|---|---|---|
| T1 | **FP** | Compton wavelength λ_C = ℏ/(m_ν·c) dla m_ν=0.1 eV — explicit numerical |
| T2 | **LIT** | Substrate L_X = ℏc/m_X z L06 anchor m_X≈60 MeV — explicit numerical |
| T3 | **FP** | g_0-weighted L_eff = λ_C·g_0_ν/g_0_e — interpolation scale |
| T4 | **FP** | Range bracket [L_min, L_max] z 3 scenarios |
| T5 | **FP** | μ_scalar(L_kink) z β-task formula — 3 scenarios |
| T6 | **FP** | μ_spinor(L_kink) z RP² formula — 3 scenarios |
| T7 | **FP** | Falsifiability check vs XENONnT/GEMMA/SM Dirac — które scenarios pass? |
| T8 | **DEC** | S05 preservation: no new free parameters introduced; inputs są TGP-native (anchor + calibration) |

**Substance ratio:** 6 FP + 1 LIT + 1 DEC = 75% FP ✓. Hardcoded T_pass=True: 0.

## §1 — Status

🟢 **ACTIVE — opened 2026-05-17 (3rd cycle sesji 2026-05-17)**

## §2 — Cross-references

- [[../op-neutrino-omega-motion-wake-2026-05-17/]] — β-task A-
- [[../op-neutrino-RP2-wake-extension-2026-05-17/]] — RP² A-
- [[../op-L06-axion-mass-derivation-2026-05-16/]] — m_X anchor source
- [[../../audyt/NUMERICAL_ANCHORS_REGISTRY.md]] — anchor classification

---

**Scaffolded:** 2026-05-17 Claudian (3rd cycle sesji β-task-resolution).
