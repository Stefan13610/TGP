---
title: "op-L06-axion-mass-derivation — Próba derywacji m_X strukturalnie + honest confirmation FREE-parameter status (L06 audit closure)"
date: 2026-05-16
type: research-cycle
folder_status: parking
parent: "[[../../TGP_FOUNDATIONS.md]]"

# ============== KICKOFF CONTRACT (BINDING post-2026-05-10) ==============
contract:
  L1_native:
    output_observable: "Wartość liczbowa m_X [MeV] z derywacji strukturalnej (jeśli udana) lub explicit obstruction proofs dla 5 candidate paths (A-E); klasyfikacja statusu m_X (DERIVED / FREE PARAMETER / NUMERICAL ANCHOR) [verdict]; cross-cycle consistency check ψ.1 (100 MeV) vs τ.3 (0.83 MeV) [boolean]"
    measurement_instrument: "Symbolic enumeration TGP available scales (Φ₀, v, γ, M_Pl, H_0, α, g_ω.1, f_X); dimensional analysis dla candidate combinations; explicit OOM check vs target 100 MeV; sympy symbolic verification dla każdego path"
    native_coefs_constrained:
      - "Path A: m_X² ?= |U''(φ_min)| substrate breathing mode → expected mismatch (γ-scale)"
      - "Path B: m_X ?= g·f_X (τ.3 derivation) → produces 0.83 MeV ≠ 100 MeV (cross-cycle inconsistency)"
      - "Path C: m_X ?= √(α·H_0·M_Pl) lub other dimensional combination → expected OOM mismatch"
      - "Path D: m_X from radiative Coleman-Weinberg analog → expected sub-MeV scale, NOT 100 MeV"
      - "Path E: m_X = FREE PARAMETER (audit § A.7 option 2 — ω.3 endorsed) → honest acknowledgment"
    falsification_rule: "Jeśli któryś z paths A-D daje m_X ≈ 100 MeV (consistency within 10%) z first-principles + brak free fitting parametru, wówczas TGP m_X jest DERIVED i audit L06 jest CLOSED-FULL (A−). Jeśli wszystkie paths A-D zawodzą z explicit structural obstructions, m_X jest CONFIRMED FREE PARAMETER (B+ partial, audit § A.7 option 2 strukturalnie verified). Jeśli paths są niewystarczające (gap w analizie), HALT-B."
    pre_registration_date: "2026-05-16"

  L2_framework_reduction:
    target_frameworks:
      - "Peccei-Quinn (1977) axion mass m_a·f_a = m_π·f_π from QCD instanton — DOES NOT APPLY to TGP ALP (no color anomaly)"
      - "Goldstone-Nambu (1960-61) — Z₂ Goldstone is massless w pure Z₂-symmetric H_Γ"
      - "Coleman-Weinberg (1973) radiative mass generation — applicable for sub-leading effects"
      - "Hierarchy problem (general): m_h vs M_Pl analog dla m_X vs available TGP scales"
    reduction_type: "consistency-mapping"
    validation_transfer: "TGP ALP m_a behaviour analog do QCD ALP m_a = FREE; nie analog do QCD axion m_a·f_a = const"
    failure_disposition: "L1-stands"

  L3_falsification_map:
    - { bound: "Audit L06 § A.7 option 2 + ω.3 endorsement", constrains: "m_X = FREE PARAMETER post-ω.3", window: "structural; closure pending forward ω.4 derivation attempt", status: "this cycle attempts forward derivation" }
    - { bound: "ψ.1 Phase1_setup.md:50-56 — m_X = 100 MeV input", constrains: "Yukawa range L_eff = 1.97 μm if m_X = 100 MeV", window: "phenomenological; SNR optimization target", status: "input remains phenomenological post-this-cycle if Path E confirmed" }
    - { bound: "τ.3 B7_greens_function_results.md — m_X = g·f_X = 0.83 MeV", constrains: "B7 KEY PHYSICS heavy-regime universal lab; m_X·L_lab ~ 4·10⁹", window: "structural; bulk signal ZERO", status: "cross-cycle inconsistency z ψ.1 (100 MeV) — both phenomenological" }
    - { bound: "ω.3 phase1_structural_setup.txt:71 — 'TGP m_a status: FREE PARAMETER'", constrains: "ω.3 explicitly defers m_X to ω.4 cycle", window: "structural; ALP classification preserved", status: "this is the ω.4 cycle attempt" }
    - { bound: "Goldstone theorem dla exactly-symmetric Z₂ substrate", constrains: "m_X = 0 strukturalnie if no explicit Z₂-breaking term in H_Γ", window: "structural; requires verification że TGP H_Γ jest Z₂-exact", status: "T-test in Phase 1" }

tgp_status:
  level: L1
  kind: derivation
  output_type: structural-verdict-or-derivation
  core_compatibility: review-only
  may_edit_core: false
  has_needs_file: false
  has_findings_file: false
  exports_findings: false
  open_bridges:
    - "Forward ω.4+ structural derivation attempt (this cycle is the attempt)"
    - "Cross-cycle m_X harmonization (ψ.1: 100 MeV; τ.3: 0.83 MeV; ω.3: FREE) — separate housekeeping"
    - "Coleman-Weinberg radiative mass detailed calculation (sub-MeV regime; if Path D needs refinement)"
  depends_on:
    - "audyt/L06_axion_mass_locked/README.md (audit issue source)"
    - "research/op-omega3-axion-decay-constant/phase1_structural_setup.py + Phase2_results.md (m_a FREE classification)"
    - "research/op-tau3-substrate-clock-acceleration/B7_greens_function_results.md (m_X = g·f_X = 0.83 MeV)"
    - "research/op-psi1-substrate-light-acceleration/Phase1_setup.md:50-56 (m_X = 100 MeV phenomenological)"
    - "research/op-L07-zero-sum-Z2-derivation-2026-05-16/ (Z₂ substrate symmetry derived earlier today)"
    - "core/sek05_ciemna_energia/sek05_ciemna_energia.tex eq.U-phi-explicit (V''(1) = -γ)"
    - "research/closure_2026-04-26/Lambda_from_Phi0/ (γ = M_Pl²·H_0² T-Λ closure)"
  impacts:
    - "audyt/L06_axion_mass_locked — Path 2 closure (structural derivation attempt completed)"
    - "ψ.1 phenomenological m_X = 100 MeV — annotation 'remains FREE post-ω.4' if Path E confirmed"
    - "τ.3 m_X = 0.83 MeV — annotation 'remains derived-from-coupling, but coupling free' if Path E confirmed"
    - "ω.3 m_a FREE status — STRENGTHENED if Path E confirmed"
    - "PREDICTIONS_REGISTRY — TT13/TT14/WW7-WW12 conditional-on-m_X annotation needed if Path E confirmed"

predecessors:
  - "[[../op-L07-zero-sum-Z2-derivation-2026-05-16/]] (Z₂ substrate symmetry derivation, used dla Goldstone theorem application)"
  - "[[../op-omega3-axion-decay-constant/]] (m_a FREE classification — this cycle is forward ω.4)"
  - "[[../op-tau3-substrate-clock-acceleration/]] (m_X = g·f_X derivation channel)"
  - "[[../op-psi1-substrate-light-acceleration/]] (m_X = 100 MeV phenomenological input)"

related:
  - "[[../../audyt/L06_axion_mass_locked/README.md]] (this cycle addresses)"
  - "[[../../audyt/L07_zero_sum_axiom/README.md]] (Z₂ structure inherited)"
  - "[[../../meta/CYCLE_KICKOFF_TEMPLATE.md]]"
  - "[[../../STATE.md]]"

classification: DERIVATION — L06 audit Path 2 (forward structural derivation attempt of m_X)
priority: high (P2 OPEN; ω.3 explicit forward-gate; cross-cycle inconsistency needs resolution)
goal: "Forward structural derivation attempt of m_X (TGP ALP axion mass) — tested via 5 candidate paths: (A) substrate breathing mode V''(1); (B) g·f_X coupling product; (C) dimensional analysis on TGP available scales; (D) Coleman-Weinberg radiative mass; (E) honest FREE-parameter acknowledgment. Pre-registered B+ outcome likely (Path E confirmation) — analog to L08 e²-derivation cycle (2026-05-16 partial closure). Resolution of cross-cycle m_X inconsistency (ψ.1 100 MeV vs τ.3 0.83 MeV)."
estimated_effort: "~1 sesja (Phase 0 + Phase 1 symbolic + Phase FINAL closure; honest partial expectation pre-registered)"
target_window: "Phase 1: 4 path tests (A, B, C, D) z structural obstruction documentation; Path E confirmation via cross-cycle audit; ~10-12 tests target"

six_requirements_target:
  - "P1: Path A: substrate breathing mode |V''(1)| = γ scale check → expected M_Pl·√H_0 ~ 10⁻³ eV ≠ 100 MeV (OOM mismatch ~8)"
  - "P2: Path B: m_X = g·f_X = 0.83 MeV (τ.3 inherited) ≠ 100 MeV (ψ.1 used) → cross-cycle inconsistency explicit"
  - "P3: Path C: dimensional combinations of (Φ₀, γ, M_Pl, H_0, α) for m_X ≈ 100 MeV → exhaustive enumeration shows no natural combination"
  - "P4: Path D: Coleman-Weinberg radiative mass from substrate-photon coupling → estimated sub-MeV scale ≠ 100 MeV"
  - "P5: Path E: m_X = FREE PARAMETER honest acknowledgment z audit § A.7 option 2 + ω.3 endorsement → structurally confirmed"
  - "P6: Cross-cycle harmonization recommendation: ψ.1 and τ.3 m_X values are BOTH phenomenological choices for different SNR scenarios (NIE sprzeczne, NIE derived)"

risk_flags:
  - "R1: Wishful thinking — pressure to 'find' m_X derivation. Mitigation: pre-registered B+ acceptable, audit option 2 endorsed"
  - "R2: Coleman-Weinberg calculation deep field-theoretic computation; cycle scope = OOM estimate, NIE full loop integral"
  - "R3: Cross-cycle inconsistency (ψ.1 100 MeV vs τ.3 0.83 MeV) — could be perceived as 'failure'; honest framing as BOTH phenomenological mitigates"
  - "R4: Path E acknowledgment status: NIE 'cycle failure' — analog do L08 e² CLOSED-NEGATIVE 2026-05-01 + numerical anchor reinforced 2026-05-16"
  - "R5: ω.3 declares m_X FREE — does forward ω.4 derivation attempt have value? YES: explicit obstruction proofs strengthen FREE status structurally"

phase_plan:
  Phase_0: "Balance sheet + 6/6 gate + scope (4 derivation paths + 1 honest acknowledgment; pre-registered B+ expectation)"
  Phase_1: "First-principles symbolic: 4 path tests + cross-cycle consistency + audit endorsement verification"
  Phase_FINAL: "Closure ceremony z honest B+ verdict (analog L08 e²-derivation 2026-05-16); audit L06 Path 2 partial closure; cross-cycle harmonization recommendations"

tags:
  - L06
  - L06-Phase1
  - axion-mass
  - m_X
  - phenomenology
  - FREE-parameter
  - structural-derivation-attempt
  - audit-L06-Path2
  - omega4-cycle
  - cycle-scaffold-2026-05-16
---

# op-L06-axion-mass-derivation-2026-05-16

> **Cel:** Forward structural derivation attempt of TGP axion mass m_X (audit L06 Path 2,
> ω.3 forward-gate ω.4+). 4 candidate paths (A breathing mode, B coupling product,
> C dimensional, D radiative) z honest acknowledgment Path E (FREE parameter) jako
> pre-registered B+ outcome. Cross-cycle harmonization ψ.1 (100 MeV) vs τ.3 (0.83 MeV).

## §0 — Cel + native-first contract

[CITE: `audyt/L06_axion_mass_locked/README.md`; `op-omega3-axion-decay-constant/phase1_structural_setup.txt:71`
"TGP m_a status: FREE PARAMETER (open omega.4+ cycle target)"]

### §0.1 — Native observable target

**Co fizycznie liczymy:**

- `m_X_path_A_result` [MeV] — substrate breathing mode |V''(1)| estimate
- `m_X_path_B_result` [MeV] — g·f_X product (τ.3 inherited)
- `m_X_path_C_results` [list MeV] — dimensional combinations of TGP scales
- `m_X_path_D_estimate` [MeV] — Coleman-Weinberg OOM estimate
- `m_X_status` ∈ {DERIVED-VALUE, FREE-PARAMETER, NUMERICAL-ANCHOR, HALT}
- `cross_cycle_consistency_status` ∈ {CONSISTENT, INCONSISTENT-BUT-PHENOMENOLOGICAL, SPRZECZNOŚĆ}

**Instrument:** Symbolic dimensional analysis; OOM comparison vs 100 MeV target;
sympy symbolic verification dla każdego path; literature anchor checks (PQ, Goldstone, CW).

### §0.2 — Pre-registered falsification rule

**Decision rule WRITTEN BEFORE any calculation (2026-05-16):**

> Jeśli któryś z paths A-D daje m_X ≈ 100 MeV (consistency within 10%) z first-principles
> + brak free fitting parametru → TGP m_X jest DERIVED, audit L06 jest CLOSED-FULL (A−).
>
> Jeśli wszystkie paths A-D zawodzą z explicit structural obstructions, m_X CONFIRMED
> FREE PARAMETER (B+ partial, audit § A.7 option 2 strukturalnie verified).
>
> Jeśli paths niewystarczające lub HALT zachodzi (np. some path needs deeper computation
> outside cycle scope), HALT-B z explicit obstruction note.

```
pre_registration_date: 2026-05-16
recovery_scope:
  allowed_directions:
    - "Subleading corrections to Coleman-Weinberg OOM estimate (within order of magnitude)"
    - "Alternative dimensional combinations of TGP scales (exhaustive enumeration)"
    - "Cross-cycle harmonization: m_X as 'phenomenological choice' rather than 'derived value'"
  forbidden_directions:
    - "Free fitting parameter introduced to force m_X = 100 MeV (would not be derivation)"
    - "Post-hoc redefinition of 'derivation' to accept order-of-magnitude only"
    - "Acceptance of 100 MeV as 'natural' without OOM check vs Φ₀, γ, M_Pl, H_0"
  if_recovery_exhausted: "Honest B+ verdict: m_X = FREE PARAMETER structurally confirmed; audit § A.7 option 2 endorsed; cross-cycle ψ.1/τ.3 inconsistency dispositioned as phenomenological choice diversity"
```

### §0.3 — TGP-native check (mandatory)

- [x] **Q1 (Pattern coverage):** Pattern 2.9 (no-derivation / numerical anchor / free param)
      relevant; analog L08 e² numerical anchor (CLOSED-NEGATIVE 2026-05-01 + reinforced 2026-05-16)
- [x] **Q2 (Red flags):** WIELKIE potential red flag — pressure to find derivation. Cycle
      explicitly mitigates z pre-registered B+ acceptable + audit option 2 endorsement
- [x] **Q3 (Inherited LOCKs):** Z₂ substrate symmetry derived (L07 cycle today); ω.3 m_a
      FREE classification LIVE; τ.3 g·f_X = 0.83 MeV LIVE; T-Λ closure γ = M_Pl²·H_0² LIVE
- [x] **Q4 (Standard-physics tools):** Peccei-Quinn (1977), Goldstone-Nambu (1960-61),
      Coleman-Weinberg (1973) standard tools; native-relevance: TGP ALP regime → PQ
      relation NIE applicable (no QCD anomaly)
- [x] **Q5 (m_Φ usage):** m_X jest TARGET observable; this cycle attempts derivation
- [x] **Q6 (GR limit framing):** N/A — particle physics question, NOT cosmology
- [x] **Q7 (ASK-RULE self-check):** methodology cited; gaps in R-flags; honest partial
      expectation pre-registered (B+ acceptable)
- [x] **Q8 (BD-drift audit plan):** Phase FINAL flag jeśli Path E confirmation introduces
      Brans-Dicke-style "free parameter wisi w teorii"

### §0.4 — Pre-flight methodology read confirmation

**BINDING per `meta/CYCLE_KICKOFF_TEMPLATE.md` §2.6:**

- [x] Przeczytano [[../../audyt/L06_axion_mass_locked/README.md]]
- [x] Przeczytano [[../op-omega3-axion-decay-constant/phase1_structural_setup.txt]]:71
      ("TGP m_a status: FREE PARAMETER (open omega.4+ cycle target)")
- [x] Przeczytano [[../op-tau3-substrate-clock-acceleration/B7_greens_function_results.md]]:97
      (m_X = g·f_X = 0.83 MeV)
- [x] Przeczytano [[../op-psi1-substrate-light-acceleration/Phase1_setup.md]]:50-56
      (m_X = 100 MeV phenomenological input)
- [x] Przeczytano [[../op-L07-zero-sum-Z2-derivation-2026-05-16/]] (Z₂ structure)
- [x] Przeczytano [[../../meta/CYCLE_KICKOFF_TEMPLATE.md]] §1-§2

**Sign-off:** Claudian (theoretical physics agent) @ 2026-05-16

### §0.5 — Sympy substance plan

**Plan testów Phase 1 (target 11 tests, ≥9 first-principles):**

| Test | Klasa | Pytanie fizyczne |
|---|---|---|
| T1 | **FIRST_PRINCIPLES** | Path A: V''(1) = -γ < 0 z sek05 eq.U-phi-explicit; tachyonic at φ=1 vacuum point |
| T2 | **FIRST_PRINCIPLES** | Path A scale: √\|V''\| = √γ = M_Pl·√H_0; explicit numerical = ~10⁻⁵ eV² → mass ~5·10⁻³ eV = 5 meV ≠ 100 MeV (OOM mismatch 11) |
| T3 | **FIRST_PRINCIPLES** | Path B: m_X = g·f_X = 8.3·10⁻³ × 100 MeV = 0.83 MeV (τ.3 derivation); ≠ 100 MeV ψ.1 input → CROSS-CYCLE INCONSISTENCY explicit |
| T4 | **FIRST_PRINCIPLES** | Path C dimensional enumeration: 6 combinations of (M_Pl, H_0, Φ₀, α, α_s) tested; żadna nie daje ~100 MeV naturally |
| T5 | **FIRST_PRINCIPLES** | Path C analog QCD: Λ_QCD = 217 MeV ≈ 2× target; ALE TGP axion is ALP (no QCD anomaly) → mismatch reasonable structurally |
| T6 | **FIRST_PRINCIPLES** | Path D Coleman-Weinberg: m_X²_CW ~ (g²/16π²)·Λ_UV² with g~10⁻³, Λ_UV~M_Pl → m_X_CW ~ M_Pl/10³ × 1/12 ≈ 10²⁵ eV — TOO BIG; restricted scale Λ_UV ~ Λ_QCD → m_X_CW ~ 10 keV ≠ 100 MeV |
| T7 | **FIRST_PRINCIPLES** | Goldstone theorem: Z₂ exactly preserved at H_Γ level (L07 derivation) ⇒ axion is Goldstone (massless) IF no explicit breaking; m_X > 0 wymaga explicit breaking term |
| T8 | **FIRST_PRINCIPLES** | TGP H_Γ Z₂-breaking source check: explicit term δH_Z₂ ABSENT z fundamental TGP Lagrangian (single-Φ axiom S05) → m_X = 0 strukturalnie jeśli pure Z₂-symmetric → conflict z observed m_X > 0 in ψ.1/τ.3 |
| T9 | **FIRST_PRINCIPLES** | Resolution: m_X > 0 wymaga emergent breaking (not fundamental). Sources: substrate-photon coupling g_ω.1 (introduces explicit electromagnetic Z₂-breaking through F·F̃ ω.1 EOM); → m_X² ~ g_ω.1²·(electromagnetic scale)²; sub-MeV estimate |
| T10 | **FIRST_PRINCIPLES** | Path E: ω.3 endorsement check — "TGP m_a status: FREE PARAMETER" (phase1_structural_setup.txt:71); audit § A.7 option 2 endorsed; this cycle confirms strukturalnie |
| T11 | **LITERATURE_ANCHORED** | Comparison: QCD axion m_a·f_a = m_π·f_π fixed by PQ instanton; TGP ALP m_a free (no QCD anomaly); analog Coleman-Weinberg pseudo-Goldstone scenarios |
| T12 | **DECLARATIVE** | S05 single-Φ + Z₂ substrate preserved; no new free parameters introduced; m_X status declassified from "locked 100 MeV" to "FREE phenomenological choice" w ψ.1/τ.3 |

**Target:** 11/11 PASS sympy (T1-T11) + 1 structural declaration (T12 separate).

**Ratio:** 10 FIRST_PRINCIPLES (90.9%) + 1 LITERATURE_ANCHORED (9.1%) + 1 DECLARATIVE separate.

---

## §1 — Phase 0: balance sheet

[Patrz `Phase0_balance.md`]

## §2 — Phase 1: native derivation

[Patrz `Phase1_sympy.py` + `Phase1_results.md`]

## §FINAL — Closure

[Patrz `Phase_FINAL_close.md`]

---

## Status

🟢 **ACTIVE — opened 2026-05-16** per user authorization "L06 axion-mass cycle".

Cycle scope: focused 4-path derivation attempt + 1 honest acknowledgment (Path E); audit
L06 Path 2 partial closure (ω.4 forward-gate); cross-cycle ψ.1/τ.3 harmonization. Full
phenomenology rewrite ψ.1/τ.3 m_X annotations DEFERRED do separate housekeeping cycle.

This session deliverables:
- README.md (this file) z BINDING contract — **DONE**
- Phase0_balance.md — **PLANNED**
- Phase1_sympy.py — **PLANNED** (11-test first-principles symbolic)
- Phase1_sympy.txt — **PLANNED** (sympy output)
- Phase1_results.md — **PLANNED**
- Phase_FINAL_close.md — **PLANNED if Phase 1 PASS or HALT-acceptable**

Pre-registered partial outcome: **B+ expected** (Path E confirmation; m_X = FREE PARAMETER
strukturalnie verified) — analog do L08-e²-derivation cycle (2026-05-16) honest partial.

---

**Cross-references:**
- [[../../audyt/L06_axion_mass_locked/README.md]] (cycle addresses this audit issue)
- [[../op-omega3-axion-decay-constant/phase1_structural_setup.txt]] (m_a FREE classification source)
- [[../op-tau3-substrate-clock-acceleration/B7_greens_function_results.md]] (m_X = g·f_X channel)
- [[../op-psi1-substrate-light-acceleration/Phase1_setup.md]] (m_X = 100 MeV phenomenological)
- [[../op-L07-zero-sum-Z2-derivation-2026-05-16/]] (Z₂ structure inherited)
- [[../../meta/CYCLE_KICKOFF_TEMPLATE.md]]
- [[../../STATE.md]]
