---
title: "op-S07-reset-alternative-f-psi — alternative f(ψ) metric structures post M9.1'' GWTC-3 5.02σ rejection [CLOSED-RESOLVED 2026-05-13 sesja P-FINAL; A−]"
date: 2026-05-11
last_updated: 2026-05-13 sesja P-FINAL (Phase FINAL closure ceremony; A−)
type: research-cycle
folder_status: closed-resolved   # CLOSED 2026-05-13 sesja P-FINAL per user authorization "Opcja A (recommended): Phase FINAL closure ceremony z claim_status A−"; cumulative 27/27 sympy PASS (22 FP / 5 LIT / 4 DEC separate); WIP slot 1/5 FREED
claim_status: A-MINUS   # STRUCTURAL_DERIVED_NATIVE z L2 not-fully-FP-attempted; H1a TENTATIVE pending observational LIGO-O5 A+ ~2027
parent: "[[../../TGP_FOUNDATIONS.md]]"

# ============== KICKOFF CONTRACT (BINDING post-2026-05-10) ==============
contract:
  L1_native:
    output_observable: "β_ppE^(b=-1) phase coefficient projection [dimensionless] dla alternative f(ψ) family member; LIGO-O5 A+ ~2027 SNR threshold dla single-event detection; Δφ(f) phase residual [radians] inherited z LIGO-3G-native methodology"
    measurement_instrument: "LIGO-Virgo-KAGRA O4-O5 + future LIGO-3G (Cosmic Explorer + Einstein Telescope) GW phase precision; GWTC-3 90 BBH posteriors current bounds"
    native_coefs_constrained:
      - "f(ψ) functional family parameters (specific to candidate; e.g., f(ψ) = 1 + α·(ψ-1) lub polynomial)"
      - "Δe_2_native = -4·ξ_3 + 4 - a_3/8 + c_0·κ_σ (z LIGO-3G-native A− cycle inheritance)"
    falsification_rule: "Jeśli wszystkie f(ψ) z S07 freedom family give β_ppE^(b=-1) outside GWTC-3 1σ window |β_ppE| ≤ 0.78 OR z LIGO-O5 A+ 5σ single-event excluded, S07 freedom INSUFFICIENT do escape M9.1'' falsification → framework wymaga architecture revision lub acceptance M9.1'' jako framework-level falsification (H1b verdict)."
    pre_registration_date: "2026-05-13"

  L2_framework_reduction:
    target_frameworks:
      - "ppE basis 2.5PN phase"
      - "Bayesian GWTC-3 90 BBH combined analysis"
      - "GR limit recovery low-velocity"
    reduction_type: "not-attempted"
    failure_disposition: "L1-stands"

  L3_falsification_map:
    - { bound: "GWTC-3 |β_ppE| ≤ 0.78 (1σ)", constrains: "f(ψ) family β_ppE prediction", window: "alternative family member must fall within", status: "pending Phase 1 enumeration" }
    - { bound: "LIGO-O5 A+ ~2027 single-event SNR threshold", constrains: "Δφ(f) phase residual", window: "inherited z LIGO-3G-native PR-002 LOCKED", status: "pending observational" }
    - { bound: "GR low-velocity limit recovery", constrains: "f(ψ) → 1 as ψ → ψ_0 cosmological", window: "structural requirement", status: "pending Phase 1" }

# ============== END KICKOFF CONTRACT ==============

tgp_status:
  level: T2
  kind: recovery
  output_type: observable
  core_compatibility: review-only
  may_edit_core: false
  has_needs_file: false
  has_findings_file: false
  exports_findings: false
  open_bridges:
    - "Full Bayesian fit alternative f(ψ) families do GWTC-3 (multi-session deferred)"
    - "BH5 QNM ringdown + ε.1 photon ring channels (separate cycles)"
  depends_on:
    - "op-newton-momentum/M9_1_pp_P1_results.md (M9.1'' canonical — FALSIFIED)"
    - "op-ppE-mapping/Phase1.5_G_SPA_lock.md (β = -15/4 sympy LOCK)"
    - "op-GWTC3-reanalysis/Phase2_RERUN_2026-05-09_corrected_beta.md (5.02σ TGP rejection)"
    - "op-emergent-metric-from-interaction-2026-05-09 (Phase 4 zero-β region {A,B,C} family)"
    - "op-LIGO-3G-native-phase-residual-2026-05-11 (A− native methodology)"
    - "core/sek07_solver/sek07_solver.tex (S07 alternative f(ψ) freedom)"
  impacts:
    - "M9.1'' Path 2 anchor recovery enumeration (per M9_RESTRUCTURE_NOTE §2)"
  source_of_status:
    - "RESEARCH_RESTART_2026-05-11 §5.2 priority 4 (deferred — reactivated 2026-05-13)"

predecessors:
  - "[[../op-newton-momentum/M9_1_pp_P1_results.md]] (M9.1'' kanoniczna metryka — FALSIFIED)"
  - "[[../op-ppE-mapping/Phase1.5_G_SPA_lock.md]] (β = -15/4 sympy LOCK)"
  - "[[../op-GWTC3-reanalysis/Phase2_RERUN_2026-05-09_corrected_beta.md]] (5.02σ TGP rejection)"
  - "[[../op-emergent-metric-from-interaction-2026-05-09/]] (Phase 4 zero-β region {A,B,C})"
  - "[[../op-LIGO-3G-native-phase-residual-2026-05-11/]] (A−)"

related:
  - "[[../../core/sek07_solver/sek07_solver.tex]] (S07 freedom)"

classification: RECOVERY z post-falsification reset
priority: medium (long-term, multi-session)
goal: "Wyprowadzić alternative f(ψ) metric structures S07-compatible że (a) recover GR limit low-velocity, (b) compatible z GWTC-3 within 1σ, (c) preserve cross-cycle N1-N5 + Q2 closures, (d) preserve S05. Brak Lakatos OR-clauses — jeśli wszystkie alternatywy fail → H1b verdict (TGP framework requires architecture revision)."
estimated_effort: "~5-8 sesji (mathematically intensive multi-session)"
target_window: "Phase 1: S07 alternative f(ψ) family enumeration (operational classes z GR-limit constraint). Phase 2: β_ppE predictions dla każdej alternative (corrected SPA chain z Phase 1.5). Phase 3: GWTC-3 Bayesian compatibility test. Phase 4: cross-cycle consistency."

six_requirements_target:
  - "P1: S07 alternative f(ψ) family enumeration (operational classes, GR-limit recovery)"
  - "P2: β_ppE prediction dla każdej alternative (corrected Phase 1.5 SPA chain)"
  - "P3: GR limit recovery verification dla każdej alternative (ψ → ψ_0 cosmological → GR)"
  - "P4: GWTC-3 Bayesian compatibility test (|β_ppE| ≤ 0.78 at 1σ)"
  - "P5: Cross-cycle consistency z L01 N1-N5 + Q2 + cluster closures preserved"
  - "P6: S05 single-Φ axiom preserved bezwarunkowo"

risk_flags:
  - "R1: Alternative f(ψ) może wymagać post-hoc tuning — BLOCKED anti-Lakatos protocol"
  - "R2: GR limit nie recovery dla some alternatives → strukturalnie ruled out (per P3)"
  - "R3: Bayesian GWTC-3 compatibility nontrivial dla each alternative"
  - "R4: BH5 (QNM ringdown) + ε.1 (photon ring) channels musi pozostać LIVE"
  - "R5: M911-P2 multi-coefficient ratios WITHDRAWN — full SPA chain re-derivation needed"
  - "R6: S05 nie violated by alternative f(ψ) (each must preserve single-Φ)"

phase_plan:
  Phase_0: "Balance sheet + pre-flight methodology read confirmation (this README + Phase0_balance.md)"
  Phase_1: "S07 alternative f(ψ) family enumeration (operational classes z GR-limit constraint, NIE Lakatos OR-clause backstop)"
  Phase_2: "β_ppE predictions dla każdej alternative (corrected SPA chain z Phase 1.5)"
  Phase_3: "GR limit + GWTC-3 Bayesian compatibility test dla każdej alternative"
  Phase_4: "Selection: which alternative passes wszystkie tests + cross-cycle consistency"
  Phase_FINAL: "Closure z H1a/H1b verdict; brak H1c/H1d backstops"

tags:
  - S07-reset
  - alternative-f-psi
  - M9.1''-falsified-followup
  - ppE-beta-correction
  - GWTC-3-compatibility
  - SPA-chain-rederivation
  - reactivated-binding-2026-05-13
  - multi-session
---

# op-S07-reset-alternative-f-psi — CLOSED-RESOLVED 2026-05-13 sesja P-FINAL

> **Status:** Scaffold HALTED 2026-05-11 → REACTIVATED 2026-05-13 z BINDING template →
> Phase 1 12/12 PASS → Phase 2 15/15 PASS → **CLOSED-RESOLVED 2026-05-13 sesja P-FINAL**
> per closure ceremony [[./Phase_FINAL_close.md]]. **claim_status: A−** (STRUCTURAL_DERIVED_NATIVE
> z L2 not-fully-FP-attempted). **H1a TENTATIVE** pending observational LIGO-O5 A+ ~2027.

## Closure summary (post 2026-05-13 sesja P-FINAL)

- **Cumulative sympy:** 27/27 PASS (Phase 1: 12/12 + Phase 2: 15/15)
- **Substance:** 22 FP (81.5%) / 5 LIT (18.5%) / 4 DEC separate; 0 hardcoded
- **6/6 P-requirements RESOLVED**
- **PR-010 status:** LOCKED-PENDING-DATA (LIGO-O5 A+ ~2027 first decisive era)
- **WIP slot 1/5 FREED**

**KEY FINDINGS:**
1. **Phase 1 linear scaling:** β_ppE^poly(α) = (15/16)·α; recovery region α ∈ [-0.832, 0.832]
2. **Phase 2 Bayesian α-mapping:** α_ML(GWTC-3) ≈ 0; LIGO-O5 A+ projection σ_α^O5 ≈ 0.266 (×3.13)
3. **Phase 2 family distinguishability:** d²f/dψ²(ψ_0) = {0, 2β_q, α²} dla {poly, quad, trans}
4. **Phase 2 cross-cycle:** Δe_2_native(α) = α/3 z M9.1'' anchor consistency exact (α=-4 → -4/3)

Patrz [[./Phase_FINAL_close.md]] pełną closure ceremony.

---

## Original cycle activation context (preserved historical)

> **Activation context:** Scaffold HALTED 2026-05-11; **REACTIVATED 2026-05-13** z BINDING template
> rewrite per `meta/RESEARCH_RESTART_2026-05-11.md` §1.2 reactivation procedure.

## §0 — Cel + native-first contract

### §0.1 — Native observable

- β_ppE^(b=-1) projection dla S07 alternative f(ψ) families [dimensionless]
- Δφ(f) phase residual [radians] inherited z LIGO-3G-native methodology
- GR limit recovery verification

### §0.2 — Pre-registered rule (anti-Lakatos)

```
pre_registration_date: 2026-05-13
recovery_scope:
  allowed_directions:
    - "f(ψ) family enumeration WITHIN S07 freedom (polynomial / transcendental / piecewise)"
    - "GR-limit recovery constraint (mandatory)"
  forbidden_directions:
    - "Post-hoc tuning specific f(ψ) form post-data (Lakatos R1)"
    - "OR-clause H1c, H1d alternatives without pre-bounded scope"
    - "Direct Φ-photon vertex (S05 violation)"
  if_recovery_exhausted: "H1b verdict — TGP framework wymaga architecture revision OR acceptance M9.1'' jako framework-level falsification"
```

### §0.3 — Q1-Q8 OK

- [x] Q1-Q8 per CYCLE_LIFECYCLE.md template

### §0.4 — Pre-flight methodology read confirmation

- [x] Przeczytano [[../../meta/PPN_AS_PROJECTION.md]] §3.1
- [x] Przeczytano [[../../meta/TGP_NATIVE_COMPUTATIONAL_PATTERNS.md]] §1-§4
- [x] Przeczytano [[../../meta/M9_RESTRUCTURE_NOTE.md]] §1.4 + §3
- [x] Przeczytano [[../../meta/CYCLE_KICKOFF_TEMPLATE.md]] §1-§2

**Sign-off:** Claudian @ 2026-05-13 (BINDING template rewrite)

### §0.5 — Sympy substance plan (multi-session)

Phase 0 (this session): structural Phase 0 balance only — Phase 1 sympy planowany dla
przyszłej sesji z dedicated WIP slot (multi-session 5-8 sesji estymata).

Plan Phase 1 (deferred):
- T1-T3 FP: f(ψ) family classification z GR-limit constraint (ψ → ψ_0 → f → 1)
- T4-T6 FP: β_ppE^(b=-1) symbolic derivation dla 3-4 alternative forms
- T7-T8 LIT: GWTC-3 1σ window |β_ppE| ≤ 0.78
- T9 DEC: anti-Lakatos commitment + S05 preservation

**Target Phase 1:** ≥6 FP + 2 LIT + 1 DEC.

## §1 — Phase 0 status

🟡 **PARKING — reactivated 2026-05-13** z BINDING template. Pre-flight: complete. Validator
status: PENDING (run po commit).

Phase 0 commit gate:
1. Validator PASS — **PENDING run**
2. PR-010 entry — **DONE**
3. User authorization "active" + WIP slot — **PENDING**

**Phase 1 sympy NOT executed w tej sesji** — long-term multi-session work (5-8 sesji estymata).
Tej sesji deliverable: clean BINDING template rewrite, validator PASS, PR-010 entry.

---

**Cycle reactivated:** 2026-05-13 (Claudian, BINDING template rewrite per RESEARCH_RESTART §1.2).

**Predecessor halt notice preserved:** 2026-05-11 halt due to 0/5 BINDING fields; rewrite
addresses all gaps (contract:: block + output_type + §0.4 + §0.5 + PR-010 entry).
