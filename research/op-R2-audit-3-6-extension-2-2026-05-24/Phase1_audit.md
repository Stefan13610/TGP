---
title: "Phase 1 — Audit 3 items + draft §3.6.11-13 sub-rules"
type: phase_audit
status: LOCKED
phase: 1
parent_cycle: op-R2-audit-3-6-extension-2-2026-05-24
execution_date: 2026-05-24
---

# Phase 1 — Audit 3 items + draft §3.6.11-13 sub-rules

**Status:** LOCKED 2026-05-24.

---

## §1 — Item 1 audit: PARTIAL category refinement

### Pattern analysis

**Two distinct PARTIAL categories observed across cycles:**

#### Category A: PARTIAL_compute
- Computation completed but result OFF threshold by small margin
- Example: γ-1 retry Phase 2 T_P2_4 sign convention (pre-reg +2π, computed -2π; 99% magnitude match, literal threshold +2π)
- Example: γ-1 Phase 3 T_P3_3 (R²_log vs R²_exp diff 0.0127 vs threshold 0.02)
- **Compute path WAS executed; result vs threshold = partial agreement**

#### Category B: PARTIAL_concept_mismatch
- TGP-native does NOT have direct observational counterpart
- Computation conceptually inapplicable
- Example: γ-3 Phase 4 T_P4_1 (Ω_m — TGP nie ma GR-equivalent ρ_critical)
- Example: γ-3 Phase 4 T_P4_2 (CMB T = 2.725 K — observational input, not TGP-derived)
- **NIE compute fail; rather TGP framework structurally NIE applicable to specific test**

### Why current cycle 1/2/7 "1 PARTIAL allowed" is insufficient

- Budget treats both categories same → 2 concept-mismatch PARTIALS w γ-3 Phase 4 violated budget
- But concept-mismatch PARTIALs honestly declare TGP scope limits, NIE methodology failures
- Forcing strict 1 PARTIAL budget when honest assessment yields 2+ concept-mismatches incentivizes:
  - Hiding mismatches behind FAIL or PASS (anti-Lakatos hazard)
  - OR artificial single-PARTIAL packaging
  - Both undermine honest reporting

### Proposed §3.6.11 sub-rule

**§3.6.11 — PARTIAL taxonomy + budget refinement (BINDING)**

> **Pre-registration must classify each FP's possible PARTIAL outcome jako jeden z dwóch:**
>
> **(a) `PARTIAL_compute`** — Computation executed; result diverges od threshold by quantifiable margin. Examples: sign convention mismatch z 99% magnitude match; R² difference 0.013 vs threshold 0.02.
>
> **(b) `PARTIAL_concept_mismatch`** — TGP-native framework structurally lacks direct counterpart for observational test. Examples: TGP nie ma GR-equivalent (Ω_m); observational input nie TGP-derived (T_CMB).
>
> **Budget:**
> - `PARTIAL_compute`: max 1 per cycle (preserves cycle 1/2/7 numerical discipline)
> - `PARTIAL_concept_mismatch`: NO hard limit; each must explicitly document why TGP-native cannot test (concept paper § / declared limit / scope statement)
>
> **Pre-registration requirement:** Each FP musi specify PASS / PARTIAL_compute / PARTIAL_concept_mismatch / DEFERRED / FAIL outcomes z explicit thresholds.

### Severity assessment

**LOW** — Methodology refinement; closes budget ambiguity. Does NOT invalidate γ-3 verdicts (Phase 4 PARTIALs were honestly declared).

### Item 1 verdict: **CLOSED → §3.6.11 drafted**

---

## §2 — Item 2 audit: Concept paper qualitative claims methodology

### Pattern analysis

**γ-3 Phase 5 specifically:**

Concept paper §5 LOCKED 2026-05-21 claim:
> "Akceleracja ekspansji | Więcej Volume → więcej frontier surface → więcej creation → positive feedback | NATIVE"

Phase 5 sympy execution revealed:
- Statement A: "Creation rate growth ∝ R²" — true (V̇ = 4πR²c)
- Statement B: "Spatial expansion R̈ > 0" — false in linear frontier-at-c model
- §5 claim CONFLATED A and B

**Implication:** Concept paper LOCKED claims do NOT automatically equal technical truths. Some §5 claims are:
- (a) Rigorously derived facts
- (b) Plausible structural arguments (need verification)
- (c) Qualitative/intuitive intuitions (potentially flawed)

LOCKED status protects against post-hoc modification but does NOT certify rigor.

### Why current methodology gap exists

- Anti-Lakatos LOCK protects against ex-post threshold modification
- But pre-LOCK rigor of qualitative claims was NOT explicitly required
- Downstream cycles (γ-1, γ-2, γ-3) treated §5 claims jako pre-registered foundation
- Phase 5 sympy odkryło że one z 6 phenomena (acceleration) had QUALITATIVE pre-registration

### Proposed §3.6.12 sub-rule

**§3.6.12 — Concept paper claim rigor classification (BINDING)**

> **Każdy concept paper claim referenced jako pre-registration foundation MUST be classified pre-LOCK as jeden z trzech:**
>
> **(I) `DERIVED`** — Rigorously derived from minimal axioms; reproducible w sympy/symbolic computation. Concept paper sekcja musi cite derivation OR link to research cycle z derivation.
>
> **(II) `STRUCTURAL_PLAUSIBLE`** — Plausible structural argument; not yet rigorously derived but reasonable inference from axioms. Concept paper sekcja musi explicitly flag jako "do verification w future cycle".
>
> **(III) `QUALITATIVE`** — Intuitive/qualitative claim; NIE technical derivation. Concept paper sekcja musi explicitly flag jako QUALITATIVE z note: "downstream cycles SHOULD NOT rely on this claim as pre-registration without independent rigor."
>
> **Pre-LOCK audit:** PRZED LOCKing concept paper (or any meta document) z pre-registered claims, all referenced claims musi mieć (I/II/III) classification explicit.
>
> **Existing LOCKED documents:** Retroactive R2 audit for classifications; identified flaws (e.g., §5 acceleration conflation) flagged for future revision PRZED new cycle pre-registration depends on them. NIE retroactive verdict modification.

### Severity assessment

**MEDIUM** — Identifies structural methodology weakness; clarifies concept paper claim status. Resolution partly depends on Item 3 (c(Φ) may rehabilitate §5 acceleration claim).

### Item 2 verdict: **CLOSED → §3.6.12 drafted**

---

## §3 — Item 3 audit: Fundamental constants identification (CRITICAL)

### Pattern analysis

**γ-3 Phase 3 specifically:**

```python
# Phase3_sympy.py:
t_sym = sp.Symbol('t', positive=True)
c_sym_val = sp.Symbol('c', positive=True)
R_late = c_sym_val * t_sym
# Implicitly: c = standard speed of light constant
```

**Concept paper basis (LOCKED 2026-05-21):**
- §1.1: "TGP = Teoria Generowanej Przestrzeni" (space emergent z Phi)
- §3.4: "NIE wprowadza GR metric explicitly (przestrzeń jest emergent z Phi, nie pre-existing manifold)"

**Implication:** Space + signal propagation are emergent properties. c is NOT a fundamental constant a priori — it is property of how Phi-substrate excitations propagate.

**Possible mechanisms** (concept paper consistent):
- (a) σ-mode group velocity in background Φ: v_g(k, m_σ(Φ))
- (b) Emergent metric g_μν(Φ) → effective c(Φ) at cosmological scales
- (c) Frontier saturation depends on local Phi → effective frontier velocity ≠ c_0

**Phase 0 §3.6.8 BINDING (drafted 2026-05-23):**
> "Implicit symmetries: U(1) phase symmetry preserved; Z₂ + RP² preserved; Cosmological principle — TO BE TESTED, not assumed"

**What's missing:** Fundamental constants are NOT enumerated. c, ℏ, G are used implicitly z standard values w numerical computations, but their TGP-native status is unaudited.

### Why this is HIGH severity

W TGP-native cosmology, c MAY emerge z Phi configuration. If true:
- γ-3 Phase 3 R(t) = c·t derivation potentially invalidates (c was assumed constant)
- F-γ-3 PASS_TARGET MAY be downgraded OR re-validated z different mechanism
- F8 LITERAL FAIL MAY flip to PASS if c(Φ) gives non-linear R(t) z R̈ > 0
- Concept paper §5 acceleration claim MAY be rehabilitated via c(Φ) mechanism

This is **most consequential** of three R1 flags.

### Proposed §3.6.13 sub-rule

**§3.6.13 — Fundamental constants identification BINDING (HIGH priority)**

> **Każdy cycle pre-registration MUST explicitly enumerate fundamental constants i numerical scales used w derivations, z classification:**
>
> **(α) `TGP_FUNDAMENTAL`** — Constant appears w minimal TGP Lagrangian/axioms as parameter (e.g., m_σ derived z λ, v). Justified by axiom set.
>
> **(β) `EMERGENT_FROM_PHI`** — Constant emerges z Phi-substrate dynamics; may vary z background Phi configuration. Examples: signal speed c (per §1.1 ontology), effective metric g_μν.
>
> **(γ) `OBSERVATIONAL_ANCHOR`** — External observational value used as input (e.g., t_universe, T_CMB). Must be honestly declared as external.
>
> **(δ) `APPROXIMATION_LIMIT`** — Constant approximation valid only w specific regime. Must declare regime explicitly + flag for full derivation w future cycles.
>
> **Pre-registration requirement:**
> - Each derivation MUST list every constant used + classification
> - If any constant classified (β) (emergent), cycle MUST verify whether constant actually behaves as constant w investigated regime, OR derive functional form explicitly
> - If classified (δ) (approximation), regime of validity MUST be stated PRZED derivation
>
> **Specific TGP-native constants requiring classification:**
> - **c (signal speed):** (β) EMERGENT_FROM_PHI per §1.1 + §3.4 ontology. Cycles must justify c=c_0 OR derive c(Φ).
> - **m_σ (sigma mode mass):** (α) TGP_FUNDAMENTAL via m_σ² = 2λv² (function of TGP parameters).
> - **G_Newton (gravitational coupling):** Currently NOT in TGP Lagrangian → undefined w current TGP scope. Per AX-DECL declared limit.
> - **t_universe (cosmological age):** (γ) OBSERVATIONAL_ANCHOR (stellar ages).
> - **T_CMB (2.725 K):** (γ) OBSERVATIONAL_ANCHOR per γ-3 Phase 4 honest declaration.

### Severity assessment

**HIGH** — Phase 0 implicit assumption audit gap potentially invalidates Phase 3-5 verdicts under c(Φ) framework. Application to γ-3' new cycle:
- Phase 0 musi enumerate c jako (β) EMERGENT_FROM_PHI
- Cycle musi derive c(Φ_frontier) functional form PRZED R(t) numerical
- F-γ-3 + F8 verdicts pod c(Φ) reanalysis — may give different results

### Item 3 verdict: **CLOSED → §3.6.13 drafted (HIGH priority)**

---

## §4 — Cross-rule consistency check

| Rule | Conflict z §3.6.1-10? | Conflict z §3.6.11/12/13? |
|------|-----------------------|---------------------------|
| §3.6.11 PARTIAL taxonomy | NO (refines PARTIAL category w cycle 1/2/7) | N/A |
| §3.6.12 Concept paper rigor | NO (orthogonal to §3.6.1-10; meta-document audit) | NO |
| §3.6.13 Constants identification | NO (extends §3.6.8 implicit assumptions BINDING) | NO |

**No conflicts.** Three sub-rules complement existing §3.6.1-10.

---

## §5 — Summary

| Item | R1 flag origin | Sub-rule | Severity | Status |
|------|---------------|----------|----------|--------|
| 1 | γ-3 Phase 4 | §3.6.11 PARTIAL taxonomy | LOW | **CLOSED → drafted** |
| 2 | γ-3 Phase 5 | §3.6.12 Concept paper rigor | MEDIUM | **CLOSED → drafted** |
| 3 | User obs 2026-05-24 | §3.6.13 Constants identification | **HIGH** | **CLOSED → drafted** |

**Audit result:** 3/3 items CLOSED z proposed sub-rules. Ready dla Phase 4 propagation.

---

## §6 — Phase 1 status: CLOSED

**Phase 1 verdict:** **All 3 items CLOSED with §3.6.11-13 sub-rules drafted.**

**Phase 4 ready:** propagate §3.6.11-13 do CALIBRATION_PROTOCOL.md + annotate γ-3 cycle + STATE.md.

---

**END OF PHASE 1 — Audit closure LOCKED 2026-05-24**
