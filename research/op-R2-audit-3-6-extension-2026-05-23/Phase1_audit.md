---
title: "Phase 1 — Audit 4 items + draft §3.6 extension (op-R2-audit-3-6-extension-2026-05-23)"
type: phase_audit
status: COMPLETE
date_completed: 2026-05-23
phase: 1
parent_cycle: op-R2-audit-3-6-extension-2026-05-23
items_audited: 4
---

# Phase 1 — Audit 4 items + draft §3.6 extension

**Status:** COMPLETE 2026-05-23.

---

## §0 — Verdict summary

| Item | Verdict | Drafted §3.6 subsection |
|------|---------|------------------------|
| EXT-§3.6-1 (sign conventions) | ✅ CLOSED | §3.6.6 |
| EXT-§3.6-2 (fit DoF equalization) | ✅ CLOSED | §3.6.7 |
| EXT-§3.6-3 (implicit assumption enumeration) | ✅ CLOSED | §3.6.8 |
| EXT-§3.6-4 (numerical precision validation) | ✅ CLOSED | §3.6.9 |

**Aggregate Phase 1:** 4 CLOSED + 0 DEFERRED + 0 ESCALATED → **R2_PASS preliminary**.

---

## §1 — EXT-§3.6-1: Sign convention explicit derivation

### §1.1 Pattern evidence

**Source:** γ-1 Phase 2 T_P2_4 HONEST FAIL.

**Specifics:**
- Pre-registered (Phase 0 §8.6): $V_{int}(L)/L_z \approx +2\pi v^2 n_1 n_2 \log(L/r_0)$ z POSITIVE slope
- Native result: $-2\pi v^2 n_1 n_2 \log(L/r_0)$ z NEGATIVE slope
- Magnitude match w 99% (|slope| = 6.2191 vs 2π = 6.2832)
- Sign convention reversed

**Root cause:** Phase 0 §8 derived form via Green function method correctly, BUT skipped sign verification z physics (same-sign vortices repel ↔ V_int decreases with L).

### §1.2 Methodology gap

Current §3.6.1-3.6.5 BINDING requires:
- (a) Analytical derivation
- (b) Symbolic computation
- (c) Phase 0 documentation
- §3.6.2 forbids heuristic shortcuts dla DERIVED VALUES

**Gap:** §3.6 NIE explicit addresses SIGNS as derived values. Sign was implicit aspect of derivation.

### §1.3 Proposed §3.6.6 draft

**§3.6.6 — Sign convention derivation (BINDING)**

For each pre-registered formula involving signs (slope coefficients, attraction/repulsion, growth/decay rates), Phase 0 analytical pre-derivation MUST include explicit sign verification via:

(a) **Physical principle:** State which physical principle determines sign (e.g., "like-charges repel in 2D Coulomb", "kink-antikink attract via Casimir-like long-range force")

(b) **Limiting case verification:** Verify sign in known limit (e.g., "at L→r_0, V_int should be high for like-sign because charges close" → positive in standard convention)

(c) **Convention statement:** Explicit convention chosen (e.g., "V_int = energy increase when bringing charges together; positive convention dla repulsion")

**Forbidden shortcut:** Assuming sign from form mathematical analogy without physical principle derivation.

**Example (γ-1 Phase 2 lesson):**
- ❌ "V_int ~ 2π log(L/r_0) → assume positive coefficient"
- ✅ "Same-sign 2D Coulomb charges REPEL → V_int(L) high at small L, low at large L → coefficient -2π for log(L/r_0)"

### §1.4 Verdict EXT-§3.6-1: CLOSED

Per decision matrix §3 README: "Sign requirement formalized + Phase 0 template slot defined" → CLOSED.

---

## §2 — EXT-§3.6-2: Fit parameter DoF equalization

### §2.1 Pattern evidence

**Source:** γ-1 Phase 3 T_P3_3 LITERAL FAIL.

**Specifics:**
- Log fit: V = A + B·log(L) — **2 parameters**, R²_log = 0.9998
- Exp+offset fit: V = C·exp(-mL) + D — **3 parameters**, R²_exp = 0.9871
- Difference 0.0127 < pre-registered 0.02 threshold
- 3-parameter exp+offset can mimic any smooth decreasing function

**Root cause:** Pre-registered comparison used asymmetric parameter counts; extra DoF biased R² toward exp fit.

### §2.2 Methodology gap

Current §3.6 NIE addresses parameter count asymmetry in fit comparisons.

### §2.3 Proposed §3.6.7 draft

**§3.6.7 — Fit parameter DoF equalization (BINDING)**

For each pre-registered FP-class falsifier involving comparison of fit forms (e.g., R²_log vs R²_exp), the comparison MUST use:

(a) **Equal parameter count:** Both fit forms with same number of free parameters (e.g., 2-param log vs 2-param exp, NIE 2-param log vs 3-param exp+offset)

(b) **OR adjusted R² (Akaike/BIC):** Information criteria accounting for parameter count if fits MUST have different DoF dla theoretical reasons

(c) **Phase 0 explicit declaration:** Number of parameters per fit form documented; jeśli unequal, justification + AIC/BIC adjustment specified

**Forbidden shortcut:** Comparing R² of fits with asymmetric parameter counts without adjustment.

**Example (γ-1 Phase 3 lesson):**
- ❌ "R²_log (2-param) vs R²_exp+offset (3-param), threshold R²_log > R²_exp + 0.02"
- ✅ "R²_log (2-param: A, B) vs R²_exp (2-param: C, m), threshold R²_log > R²_exp + 0.02"
- ✅ "ΔAIC_log < ΔAIC_exp by margin 4 (significantly preferring log)"

### §2.4 Verdict EXT-§3.6-2: CLOSED

Per decision matrix §3 README: "DoF equalization protocol formalized + comparison rule explicit" → CLOSED.

---

## §3 — EXT-§3.6-3: Implicit assumption enumeration

### §3.1 Pattern evidence

**Source:** FFS pre-screening T7 (revealed Phase 4 2026-05-20).

**Specifics:**
- Pre-screening T7: σ = π·v² (implicit assumption q=1 effective)
- Phase 4 strict Nielsen-Olesen: σ = π·q²·v² (z q=1/3 fractional)
- Implicit q=1 absorbed silently w convenient pre-screening formula
- Phase 4 reveal: factor 10 quantitative gap depending on interpretation

**Root cause:** Pre-screening formula assumed q=1 for "effective" purposes; assumption was implicit, NIE enumerated explicit.

### §3.2 Methodology gap

Current §3.6 NIE addresses implicit assumption enumeration. Forbids heuristic shortcuts BUT some assumptions are by nature implicit (background, normalization, limit choices).

### §3.3 Proposed §3.6.8 draft

**§3.6.8 — Implicit assumption enumeration (BINDING)**

For each pre-registered formula z derived value, Phase 0 analytical pre-derivation MUST enumerate explicit:

(a) **Background assumptions:** What is assumed about field configuration vs general case (e.g., "spherically symmetric ansatz", "linearized about VEV")

(b) **Normalization conventions:** Field units, dimensional choices (e.g., "Phi has dimensions of mass", "ρ scaled by v")

(c) **Limit choices:** Which limits taken (e.g., "small-coupling λ << 1", "non-relativistic", "static")

(d) **Effective parameter substitutions:** Any parameter set to convenient value (e.g., "q=1 effective in pre-screening", "v=1 numerical normalization")

(e) **Implicit symmetries:** Symmetries imposed for tractability (e.g., "z-translation invariant")

**Phase 0 BALANCE SHEET MUST include explicit "assumptions enumeration" subsection per FP-class falsifier.**

**Forbidden shortcut:** Using convenient pre-screening or simplified formulas without explicit assumption documentation.

**Example (FFS T7 q=1 implicit lesson):**
- ❌ "σ = π·v² (string tension formula)" — implicit q=1
- ✅ "σ = π·q²·v² (Nielsen-Olesen z winding q); pre-screening uses q=1 effective for order-of-magnitude estimate; strict q=1/3 evaluated separately in full cycle"

### §3.4 Verdict EXT-§3.6-3: CLOSED

Per decision matrix §3 README: "Implicit assumption checklist drafted + Phase 0 slot defined" → CLOSED.

---

## §4 — EXT-§3.6-4: Numerical precision validation

### §4.1 Pattern evidence

**Source:** CE-H Phase 3 T_P3_2 HONEST FAIL (also covered by current §3.6, but pattern shows additional aspect).

**Specifics:**
- Pre-registered: m=1.0 expected decay rate
- Analytical correct: m·√2 ≈ 1.4142
- Numerical fitted: 1.40 (match w 1% analytical)
- 40% off pre-registered threshold

**Root cause:** Pre-registration used heuristic m=1.0 without verifying numerical precision via independent analytical computation.

### §4.2 Methodology gap

Current §3.6.1 requires "analytical derivation" but doesn't explicit specify precision standard (5%, 10%, 1%).

### §4.3 Proposed §3.6.9 draft

**§3.6.9 — Numerical precision validation (BINDING)**

Pre-registered numerical thresholds dla FP-class falsifiers MUST be validated against analytical computation z **explicit precision standard**:

(a) **5% accuracy standard:** Analytical pre-derivation must produce expected value w **±5% accuracy** of subsequent sympy execution (unless explicit looser tolerance justified physically)

(b) **Precision documentation:** Phase 0 documents:
   - Analytical expected value
   - Analytical derivation steps (verifying each)
   - Expected precision basis (5% default; relaxation justified)
   - Sympy verification plan

(c) **Mismatch handling:** Jeśli analytical pre-derivation gives value differing from sympy by > 5% (beyond expected numerical noise), Phase 0 MUST be RE-AUDITED before LOCK (anti-Lakatos discipline still preserved: NIE modification post-LOCK)

**Forbidden shortcut:** Pre-registering numerical thresholds based on heuristic intuition without analytical precision validation.

**Example (CE-H T_P3_2 lesson):**
- ❌ "Expected m = 1.0 (heuristic intuition)" → fitted 1.40, 40% off
- ✅ "Expected m·√2 ≈ 1.4142 (analytical derivation z tanh tail expansion); threshold tolerance 10% → fitted 1.40 within 1% of analytical → PASS"

### §4.4 Verdict EXT-§3.6-4: CLOSED

Per decision matrix §3 README: "Precision validation requirement formalized + standard 5% accuracy" → CLOSED.

---

## §5 — Aggregate Phase 1 verdict

| Item | Verdict | Justification |
|------|---------|---------------|
| EXT-§3.6-1 | ✅ CLOSED | Sign convention §3.6.6 drafted z 3-step requirement |
| EXT-§3.6-2 | ✅ CLOSED | DoF equalization §3.6.7 drafted z equal-param OR AIC/BIC |
| EXT-§3.6-3 | ✅ CLOSED | Implicit assumption §3.6.8 drafted z 5-category enumeration |
| EXT-§3.6-4 | ✅ CLOSED | Numerical precision §3.6.9 drafted z 5% standard |

**4 CLOSED + 0 DEFERRED + 0 ESCALATED → R2_PASS preliminary.**

---

## §6 — Anti-Lakatos check

- ✅ Each item verdict per decision matrix §3 LOCKED
- ✅ No threshold modifications ex post
- ✅ No new items added mid-cycle (forbidden #1)
- ✅ Drafted addendum text references specific patterns explicit (NIE generic rules)
- ✅ Examples cite specific cycles + tests (concrete, NIE Lakatos abstraction)
- ✅ Self-irony cascade acknowledged Phase 0 §8

---

## §7 — Implications dla Phase 4

**Phase 4 task:** Edit CALIBRATION_PROTOCOL.md adding §3.6.6-3.6.9 BINDING (4 sub-rules drafted §1.3, §2.3, §3.3, §4.3).

**Verification:** Drafted text reviewed Phase 1; ready for actual edit Phase 4.

---

## §8 — Pattern resolution mapping

| Pattern instance | Cycle | Mitigation rule |
|------------------|-------|------------------|
| 1 (T_P3_2 m vs m·√2) | CE-H Poziom β | §3.6.9 (numerical precision validation) |
| 2 (T_P4_3 q=1 implicit) | FFS Phase 4 | §3.6.8 (implicit assumption enumeration) |
| 3 (T_P2_4 sign) | γ-1 Phase 2 | §3.6.6 (sign convention derivation) |
| 4 (T_P3_3 DoF) | γ-1 Phase 3 | §3.6.7 (fit DoF equalization) |

**Each pattern instance addressed by specific §3.6 sub-rule.** Coverage complete dla observed patterns.

**Future-proofing:** New patterns (5th, 6th instance) may emerge and trigger further extensions; methodology evolution legitimate per Phase 0 §8.

---

**END OF PHASE 1 — Audit 4 items + draft §3.6 extension LOCKED 2026-05-23**

**Aggregate Phase 1:** 4 CLOSED + 0 DEFERRED + 0 ESCALATED. Drafted §3.6.6-3.6.9 ready for Phase 4 propagation.
