---
title: "op-CE-H-3D-native-interaction-retry-2026-05-23 — Poziom γ-1 retry z §3.6 extension compliance"
type: research_cycle
status: PRE_REGISTERED_LOCKED
folder_status: active
pre_registration_date: 2026-05-23
parent_cycle: op-CE-H-3D-native-interaction-2026-05-22 (γ-1 A- conditional; original)
parent_audit: op-R2-audit-3-6-extension-2026-05-23 (R2_PASS; §3.6 extension BINDING source)
predecessor_concept_paper: meta/TGP_GENERATED_SPACE_COSMOLOGY_2026-05-21.md
hypothesis_codes:
  - CE-H structural feature (R3 3/3 accepted)
  - F-γ-1 CRUCIAL TEST (LOCKED 2026-05-21, original; retry z corrected pre-reg)
  - F-γ-2 SECONDARY (LOCKED 2026-05-21; conditional na F-γ-1 PASS)
test_level: NATIVE_3D
test_scope: defect-defect interaction at large separation w 3D — retry z §3.6 extension compliance
methodology: native equations FIRST + §3.6.6-3.6.10 BINDING enforcement
claim_status_target: A | A- (depending on F-γ-1 retry verdict)
authorization_chain:
  - "2026-05-23: User 'ok zgadzam się z rekomendacjami działaj' — Path A authorization"
discipline:
  - anti-Lakatos LOCKED
  - strict cycle 1/2/7 (0 hardcoded FP T_pass=True)
  - max 1 DEC budget per cycle
  - R1+R2+R3 BINDING
  - §3.6 EXTENSION BINDING (2026-05-23, this is first cycle post-§3.6.6-10 BINDING)
falsifiers_pre_registered: F-γ-1 + F-γ-2 (same form as original; pre-reg numerical values corrected per §3.6 extension)
---

# op-CE-H-3D-native-interaction-retry — Poziom γ-1 retry BINDING contract

**Pre-registration date:** 2026-05-23
**Status:** LOCKED — żadnych modyfikacji ex post bez HALT-B.

---

## §0 — Origin i scope (retry rationale)

Niniejszy cykl jest **retry oryginalnego Poziom γ-1** ([[../op-CE-H-3D-native-interaction-2026-05-22/]] A- conditional, F-γ-1 LITERAL FAIL + SUBSTANTIVE PARTIAL).

**Retry trigger:** R2 §3.6 extension audit ([[../op-R2-audit-3-6-extension-2026-05-23/]] R2_PASS 2026-05-23) propagated CALIBRATION_PROTOCOL §3.6.6-3.6.10 BINDING addressing 4 pattern aspects (sign, DoF, implicit, precision) which caused original γ-1 LITERAL FAILs.

**Anti-Lakatos discipline preserved:**
- Original γ-1 cycle LOCKED z A- conditional verdict (immutable per audit trail invariant)
- Retry is **NEW cycle z NEW pre-registration** (NIE retroactive modification)
- Phase 0 z §3.6 extension compliance — CORRECTED pre-registration applies forward only

**Cel:** Test F-γ-1 CRUCIAL z §3.6 extension applied, expecting CLEAN PASS (both literal + substantive).

---

## §1 — Pre-registered falsifiers (LOCKED 2026-05-23, z §3.6 extension applied)

### F-γ-1 — Native 3D long-range interaction (CRUCIAL TEST)

**Pre-registracja (same as original cycle; substance unchanged):**

V_int(L) przy L → ∞ MUSI mieć formę power-law lub logarithmic lub Coulomb-like — CLEARLY distinguishable from pure exponential.

**Pre-registered numerical threshold (CORRECTED per §3.6.6-9 BINDING):**

(a) **Sign convention (§3.6.6):** Slope of log fit = **-2π** (per physical principle: same-sign 2D Coulomb charges REPEL → V_int decreases with L → coefficient NEGATIVE for log(L/r_0) form)

(b) **Fit DoF equalization (§3.6.7):** Compare **2-param log** (A + B·log(L)) vs **2-param exp** (C·exp(-mL)) — equal parameter counts; NIE 3-param exp+offset

(c) **PASS criteria:**
  - R²_log ≥ 0.95 (same as original)
  - R²_log > R²_exp + 0.02 z EQUAL-PARAM comparison (corrected per §3.6.7)
  - Slope magnitude |B| ≈ 2π within 5% (precision §3.6.9)

(d) **Substantive criterion (§3.6.6 verification):** Slope SIGN = negative (physical principle: repulsion)

### F-γ-2 — Self-consistency closure z native bg (SECONDARY)

Same as original; CONDITIONAL na F-γ-1 PASS.

### F-γ-3, F-γ-4 — NOT activated (cosmological scope)

Pozostają PENDING_POZIOM_GAMMA_3+.

---

## §2 — 10 forbidden post-hoc moves (inherited)

Same as original γ-1 + extension:
1. Modyfikacja F-γ-1 corrected thresholds ex post
2. Renaming "exponential" to avoid FAIL
3. Adding fields by rescue
4. Cherry-picking parameter range
5. Re-defining "long-range" by include mid-distance
6. Switching defect type mid-cycle (DEC reserved)
7. Hardcoding FP T_pass=True
8. Using DEC budget > 1
9. Introducing new axioms by rescue
10. Bypassing §3.6 extension (§3.6.6-10 BINDING)

---

## §3 — L1/L2/L3 falsification map

### L1 (cycle-local)
- F-γ-1 z §3.6 extension applied (corrected sign + equal-param fits + precision)
- F-γ-2 self-consistency (conditional)

### L2 (framework targets)
- CE-H Poziom β A−→A clean upgrade trajectory
- FFS C6 PARTIAL → RESOLVED disposition

### L3 (cross-cycle propagation, conditional na PASS)
- §3.6 extension first practical PASS validation
- Path B Poziom γ-2 + γ-3 sequence enabled

---

## §4 — Phase plan (5 substantive faz + Phase 0 + FINAL)

### Phase 0 — Balance sheet + §3.6 extension compliance (BINDING)
**Scope:** Explicit compliance z §3.6.6 (sign), §3.6.7 (DoF), §3.6.8 (implicit), §3.6.9 (precision). Pre-registered slope = -2π.
**Deliverable:** `Phase0_balance.md`

### Phase 1 — 3D vortex ansatz + EL + mass spectrum
**Scope:** Same as original γ-1 Phase 1 (unchanged substance); reuse 4/4 PASS result.
**Deliverable:** `Phase1_results.md` (annotation referencing original Phase 1 sympy)

### Phase 2 — Two-vortex V_int(L) z corrected pre-reg sign
**Scope:** Same numerical computation (reuse original V_int data); CORRECTED pre-reg sign verification.
**Deliverable:** `Phase2_results.md` (annotation z corrected pre-reg)

### Phase 3 — Differential test z 2-param equal-param fits
**Scope:** NEW sympy z 2-param exp fit (NIE exp+offset). Recompute R²_exp z fair comparison.
**Deliverable:** `Phase3_sympy.py` + `Phase3_results.md`

### Phase 4 — Conditional F-γ-2 self-consistency
**Scope:** ONLY IF F-γ-1 PASS clean. Self-consistency closure z native log bg.
**Deliverable:** `Phase4_sympy.py` + `Phase4_results.md`

### Phase FINAL — Closure
**Scope:** F-γ-1 + F-γ-2 verdicts, claim_status, cross-cycle propagation.
**Deliverable:** `Phase_FINAL_close.md`

**Total estimated:** 1-2 dni (most substance reused; key new work is corrected Phase 3 + Phase 4 conditional).

---

## §5 — Discipline declarations

### Strict cycle 1/2/7 inherited
- 0 hardcoded FP T_pass=True
- Substantive FP tests compute-then-compare

### Max DEC budget = 1
- Same DEC choice as original (vortex line default; point defect fallback)
- Reserved 0/1 currently

### R1+R2+R3 BINDING (CALIBRATION §6)

### §3.6 EXTENSION BINDING (this cycle = first practical application)

### Anti-BD-drift LOCK

### Anti-Lakatos LOCK
- Pre-registration LOCKED 2026-05-23 before any sympy
- §3.6 extension applied **forward only** (original γ-1 preserved at A-)

---

## §6 — Reuse vs new sympy clarification

**Reused (substance unchanged):**
- Phase 1: EL equation + mass spectrum (analytical, dimensional — same result)
- Phase 2: V_int(L) numerical data (same physics, same integration)

**New (corrected pre-reg or different fit):**
- Phase 2 pre-registered sign verification (corrected to -2π per §3.6.6)
- Phase 3: 2-param exp fit (NIE 3-param) per §3.6.7
- Phase 4: F-γ-2 (originally NOT executed; now activated if F-γ-1 PASS)

**Anti-Lakatos clarity:** Original γ-1 phase results IMMUTABLE per audit trail invariant. Retry uses same NUMERICAL DATA (substance) z CORRECTED PRE-REGISTRATION (methodology).

---

## §7 — Risk register

| ID | Risk | Mitigation | Severity |
|----|------|-----------|----------|
| R1 | Phase 3 R²_exp z 2-param fit higher than expected | Pre-registered threshold absolute (R²_log > R²_exp + 0.02); honest report if fail | LOW (analysis suggests pure exp fits poorly) |
| R2 | New methodology gap revealed in this cycle | R1 flag dla future R2 audit (5th instance scenario; methodology evolution) | LOW |
| R3 | Phase 4 F-γ-2 reveals additional issue | Pre-registered honest scenario; A- conditional acceptable | MEDIUM |
| R4 | C6 over-claim post clean PASS | Explicit C6 disposition: RESOLVED pending γ-2 + γ-3 (NIE full closure z γ-1 alone) | MEDIUM |
| R5 | Retry pattern (retry of retry of retry) — slippery slope | Anti-Lakatos LOCK: NO retroactive modifications; each retry is NEW cycle | LOW |

---

## §8 — Authorization gate

- ✅ 2026-05-23: γ-1 retry authorized
- ⏳ Phase 1-3 execution (batch authorization expected)
- ⏳ Phase 4 conditional na F-γ-1 PASS
- ⏳ Phase FINAL

**Single-session execution intended** (most substance reused; key new work limited).

---

**END OF README — γ-1 retry BINDING contract LOCKED 2026-05-23**
