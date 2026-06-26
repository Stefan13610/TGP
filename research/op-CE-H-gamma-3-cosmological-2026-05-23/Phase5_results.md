---
title: "Phase 5 results — F8 Acceleration emergence (POSITIVE PREDICTION)"
type: phase_results
status: LOCKED
phase: 5
parent_cycle: op-CE-H-gamma-3-cosmological-2026-05-23
execution_date: 2026-05-23
substantive_FP: 4
hardcoded_FP: 0
PASS_count: 1
FAIL_count: 3
PARTIAL_count: 0
DEFERRED_count: 0
F8_verdict: FAIL_LITERAL
R1_flag: candidate_2
---

# Phase 5 results — F8 Acceleration emergence (POSITIVE PREDICTION)

**Status:** LOCKED 2026-05-23.
**Execution:** Phase5_sympy.py.
**Outcome:** F8 = **FAIL LITERAL** (3 substantive FP FAIL); pre-registered severity POSITIVE (NIE KILLER).
**Critical finding:** Concept paper §5 "positive feedback → acceleration" claim **CONFLATED** "creation rate growth" z "spatial expansion acceleration".

---

## §1 — Fixed-point results

### T_P5_1 — ä derivation from R(t) = c·t — **FAIL**

**Computed:**
- R(t) = c·t (Phase 3 model)
- ȧ = c, **ä = 0**, q = 0
- Observed: q_0 ≈ -0.55 (acceleration ä > 0)

**Verdict:** **FAIL** — TGP linear gives ä = 0, observed ä > 0.

### T_P5_2 — Positive feedback analysis — **FAIL**

**Computed:**
- dV/dt = 4πR²c (grows quadratically in R)
- But R̈ = 0 (linear R(t) → no acceleration)
- "Creation rate growth ∝ R²" ≠ "spatial expansion R̈ > 0"

**Verdict:** **FAIL** — Concept paper §5 claim challenged.

**Critical observation:** Concept paper §5 "Akceleracja ekspansji | Więcej Volume → więcej frontier surface → więcej creation → positive feedback | NATIVE" appears to CONFLATE two distinct notions:
- Creation rate growth ∝ R² ✓ (true)
- Spatial expansion R̈ > 0 ✗ (NIE follows w linear frontier-at-c model)

### T_P5_3 — w_eff vs -1 — **FAIL**

**Computed:**
- TGP linear R = c·t maps to Milne (curvature-dominated empty universe)
- Effective w_eff = **-1/3**
- Pre-registered: w_DE ∈ [-1.2, -0.8] (±0.2 from -1)
- -1/3 ≈ -0.333 OUTSIDE [-1.2, -0.8]

**Verdict:** **FAIL** — TGP w_eff = -1/3 NIE matches observed w_DE ≈ -1.

### T_P5_4 — F8 honest verdict — **PASS** (meta-declaration)

**Honest declaration provided:**
- LITERAL: F8 = FAIL (3/4 sub-tests FAIL)
- STRUCTURAL: TGP-native non-linear R(t) could PASS (NIE derived w γ-3 scope)
- DEFERRED reinterpretation: SN + CMB + BAO w TGP framework

**Verdict:** **PASS** (honest meta-declaration counts).

---

## §2 — F8 PRIMARY VERDICT

**Pre-registered (concept paper §7 F8):**
> "TGP MUSI **naturalnie** dać accelerating expansion late universe (w_DE ≈ -1 ± 0.2) jako konsekwencję frontier creation, **bez** ad-hoc tuning."
>
> "Tolerancja: w_DE ∈ [-1.2, -0.8]"
>
> "Severity: POSITIVE — jeśli wymaga ad-hoc parameters by uzyskać accelerating expansion, strukturalna porażka argumentu 'naturalności'."

**Phase 5 verdict:**

$$\boxed{F8 = \text{FAIL LITERAL}}$$

**Reasoning:**
- TGP-native simplest model (R = c·t linear) gives ä = 0 (NOT ä > 0)
- w_eff = -1/3 (OUTSIDE [-1.2, -0.8])
- "Naturalnie" claim is FALSIFIED for current TGP model
- NIE ad-hoc tuning attempted (anti-Lakatos preserved)

**Severity per pre-registration:** POSITIVE — "strukturalna porażka argumentu naturalności" but NIE HALT-B trigger.

**Implications:**
- γ-3 cycle claim_status: cannot be A+ (F8 FAIL); likely **A- or B+** lub HALT-B (user decision)
- Concept paper §5 acceleration claim → REQUIRES UPDATE
- TGP-native cosmology limits identified explicitly

---

## §3 — Cycle 1/2/7 compliance

| Aspect | Status |
|--------|--------|
| Substantive FP | 4/4 (1 PASS + 3 FAIL) |
| Hardcoded T_pass | 0 ✓ |
| Anti-Lakatos | ✓ (no ad-hoc rescue attempted) |
| §3.6 BINDING | ✓ |
| 0 hardcoded FP T_pass=True | ✓ |

---

## §4 — R1 flag CANDIDATE #2

### Pattern detected:

**Concept paper §5 strukturalna siła argumentu** claimed "6 phenomena emerge from one principle". Phase 5 challenges acceleration phenomenon — TGP-native simplest model does NOT give ä > 0.

**Methodology observation:** Concept paper §5 made claim about "positive feedback gives acceleration" without rigorous derivation. Phase 5 sympy execution reveals this was QUALITATIVE/INTUITIVE claim, not technically derived.

**This is methodology pattern:** Concept paper claims may need **R2 audit** dla rigorous derivation BEFORE pre-registration LOCK.

### Implications dla §3.6 extension

Possible **§3.6.11 sub-rule candidate:** "Concept paper qualitative claims require explicit pre-registration LOCK audit BEFORE downstream cycle pre-registration."

**Defer to:** Future R2 audit (post γ-3 closure).

---

## §5 — Anti-Lakatos summary

| Check | Status |
|-------|--------|
| Modify F8 threshold ex post? | NO ✓ ([-1.2, -0.8] STAYS) |
| Rename F8 to avoid FAIL? | NO ✓ (F8 LITERAL FAIL declared) |
| Add ad-hoc acceleration mechanism? | NO ✓ |
| Reinterpret R = c·t to non-linear? | NO ✓ (linear model pre-registered Phase 3) |
| Hide FAIL behind PARTIAL? | NO ✓ (3 FAIL declared) |
| Hidden falsification? | NO ✓ (F8 FAIL openly declared) |

**Anti-Lakatos LOCK: PRESERVED.**

---

## §6 — γ-3 cycle implications

### Current verdict tracking:

| Falsifier | Pre-registered Severity | Phase 5 status |
|-----------|-------------------------|----------------|
| F-γ-3 (F4) Hubble | PRIMARY KILLER | **PASS_TARGET** ✓ |
| F5 Ω_m | SECONDARY KILLER | PARTIAL (Phase 4) |
| F6 CMB | HARD CONSTRAINT | PARTIAL (Phase 4) |
| F7 BBN | HARD CONSTRAINT | DEFERRED (Phase 4) |
| F8 Acceleration | POSITIVE PREDICTION | **FAIL LITERAL** ⚠ |
| F9 Local creation | NULL CONSISTENCY | TBD Phase 6 |

### Concept paper §5 strukturalna siła argumentu

**Original claim:** "6 niezależnych observed phenomena wypada z jednej zasady. Jeśli choć jedno z tych przewidywań fail, CE-H jest falsyfikowane."

**Phase 5 reading:**
- F-γ-3 PASS_TARGET: phenomenon (Hubble H_0) ✓ derived
- F8 FAIL: phenomenon (acceleration) NIE derived naturally
- → §5 statement implies CE-H falsified

**Phase 5 disposition (honest):** Strict reading suggests HALT-B trigger; §7 severity POSITIVE suggests less severe. **AMBIGUITY** — user decision required Phase FINAL.

---

## §7 — Cross-references dla update

**Concept paper §5 update candidates (TBD post γ-3 close):**
- "Akceleracja ekspansji | NATIVE (do verification)" → revise to "REQUIRES non-linear R(t) extension; NIE z linear frontier-at-c"
- "w_DE ≈ -1 | NATIVE (do verification)" → revise based on Phase 5 finding

**Concept paper §7 F8 update:** Current pre-registered status PRE-COMPUTATION; Phase 5 = FAIL LITERAL. Update post γ-3 close.

---

## §8 — Phase 5 status: CLOSED

- ✅ 4/4 substantive FP completed
- ✅ 0 hardcoded T_pass
- ✅ F8 FAIL declared honestly
- ✅ Anti-Lakatos preserved (no rescue attempted)
- ⚠ R1 flag CANDIDATE #2 declared (concept paper qualitative claim methodology)
- ✅ Concept paper §5 conflation identified honestly

**Phase 5 verdict:** **F8 LITERAL FAIL + R1 flag CANDIDATE #2 + Concept paper §5 needs update.**

**Phase 6 ready:** F9 (NULL CONSISTENCY) + F-γ-4 (SPECULATIVE).

---

**END OF PHASE 5 RESULTS — LOCKED 2026-05-23**
