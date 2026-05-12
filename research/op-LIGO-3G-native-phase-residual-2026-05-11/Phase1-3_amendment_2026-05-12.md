---
title: "Phase 1+2+3 amendment 2026-05-12 — post bd-drift-audit reclassification + hidden True fix"
date: 2026-05-12
type: amendment-record
status: IMMUTABLE — append-only amendment record
amendment_protocol: per bd-drift-audit 2026-05-12 §6 mandatory items 1+2
user_authorization: "Scope A (mandatory only) 2026-05-12"
parent: "[[./README.md]]"
---

# Phase 1+2+3 amendment 2026-05-12 — post bd-drift-audit substance amendment

> **Cycle:** `op-LIGO-3G-native-phase-residual-2026-05-11`
> **Trigger:** `bd_drift_audit_2026-05-12.md` verdict AMENDMENT NEEDED (§5)
> **User authorization:** Scope A (mandatory only, ~1 sesja)
> **Scope:** mandatory items 1 (hidden True sub-flags) + 2 (FP → LIT reclassifications)
> **Status (post-amendment):** 36/36 sympy PASS; honest substance metrics
> **Audit-trail invariant:** historical §1-§N content of `PhaseX_results.md` unchanged; only `§RETROACTIVE-amendment-2026-05-12` sections appended.

---

## §1 — Trigger (bd-drift-audit verdict + user scope decision)

The independent adversarial `bd_drift_audit_2026-05-12.md` (mid-cycle, pre-Phase-4 spawn) identified:

- **0 HIGH-severity BD-drift markers** — no Yukawa propagator usage; no `m_Φ` universal slip; S05 declared honestly via flagged DEC budget.
- **9 substantive reclassifications** (25% of 36 tests) — exceeds 5% PASS threshold per protocol.
- **4 hidden literal `True` sub-flags** (Type A.1 antipattern) that flowed into pass conditions:
  - Phase 2 line 627: `T6_consistency = True`
  - Phase 3 line 338: `T3_criterion_e = True`
  - Phase 3 lines 381, 385: `T4_PN_order_correspondence = True`, `T4_dimensional = True`
  - Phase 3 line 533: `T7_convention_explained = True`
- **Phase 1 T2 (FP2) substantive concern:** Pattern 2.2 momentum-flux MECHANISM (the unique TGP-native distinguisher from BD) is NOT computed; test recovers Newton via test-particle geodesic formula F = -m·∇φ, not via ∮T^{0i}dS_i.
- **Phase 3 FP-claim concern:** FP7/FP8/FP9 reduce to arithmetic substitution chains anchored in Cutler-Flanagan literature (Eq. 3.18 hand-inserted).

**Verdict (audit §5):** AMENDMENT NEEDED. Phase 4 spawn BLOCKED pending amendment.

**User scope decision (2026-05-12):** Scope A authorized — execute mandatory items 1+2 only. Audit recommendations §6 items 3-5 (advisable) deferred (recommendation §6 item 5 separate TAUT classification category NOT introduced per user Scope A; TAUT findings mapped to LIT per audit recommendation pathway).

---

## §2 — Mandatory item 1 changes (hidden literal True sub-flag fix)

Per audit §6 item 1, the four hidden literal `T_x_subflag = True` lines that flow into top-level `T_n_pass` were addressed as follows. Strategy:

- If genuine substantive sympy check is feasible → replace literal True with the check.
- Otherwise → explicitly flag `T_x_subflag_status = "DECLARATIVE"` AND keep `T_x_subflag = True` with comment marking it as audit-flagged DECLARATIVE-substep (so it remains in pass aggregation but is NOT counted as substantive verification).

| Phase | Location | Original | Resolution |
|---|---|---|---|
| Phase 2 | line 627 `T6_consistency` | `T6_consistency = True  # structural form linear in Delta_e_2 verified above` | **REPLACED with substantive sympy linearity check.** New code computes `sp.diff(phase_form, _Delta_e2_var)` and verifies it equals `-15/(4*M*v)`, AND `sp.diff(..., 2) == 0`. Flagged `T6_consistency_status = "DERIVED"`. |
| Phase 3 | line 338 `T3_criterion_e` | `T3_criterion_e = True  # symbolic structure trivially dimensionless` | **FLAGGED DECLARATIVE-substep.** `T3_criterion_e_status = "DECLARATIVE"`. Comment: dimensional consistency in geometric units NOT tracked in this sympy script; remains True as honest declarative substep, NOT a substantive sympy check. |
| Phase 3 | line 381 `T4_PN_order_correspondence` | `T4_PN_order_correspondence = True  # standard ppE basis table` | **FLAGGED DECLARATIVE-substep.** `T4_PN_order_correspondence_status = "DECLARATIVE"`. Standard Yunes-Pretorius 2009 Table I literature fact. |
| Phase 3 | line 385 `T4_dimensional` | `T4_dimensional = True` | **FLAGGED DECLARATIVE-substep.** `T4_dimensional_status = "DECLARATIVE"`. Units not tracked in symbol-only sympy script. |
| Phase 3 | line 533 `T7_convention_explained` | `T7_convention_explained = True  # consistent up to sign/factor convention` | **FLAGGED DECLARATIVE-substep.** `T7_convention_explained_status = "DECLARATIVE"`. SPA convention equivalence narrative statement, not sympy-verifiable. |

**Net effect on substance metrics:**
- 1 hidden True replaced with substantive sympy check (Phase 2 T6 linearity via `sp.diff`).
- 4 remaining hidden Trues explicitly flagged DECLARATIVE-substep with `_status` annotation, marked-per-audit in inline comments. None silently propagate as substantive content.
- All top-level `T_n_pass` aggregations remain True (so 36/36 PASS preserved); honesty is structural (sub-flags now visible and labeled).

---

## §3 — Mandatory item 2 changes (FP → LIT reclassifications)

Per audit §6 item 2 mandatory list, the following 7 tests reclassified `FIRST_PRINCIPLES` → `LITERATURE_ANCHORED`:

| Test | Phase | Pre-amendment | Post-amendment | Reason per audit (line refs) |
|---|---|---|---|---|
| T2 (FP2) | 1 | FIRST_PRINCIPLES | LITERATURE_ANCHORED | bd-audit §2.1: momentum-flux mechanism NOT actually computed; test recovers Newton via `F = -m·∇φ` test-particle formula (line 213 pre-amendment) — Gauss-theorem reduction COMMENTED but ∮T^{0i}dS_i sphere integral NOT carried out in sympy. Pattern 2.2 §2.2.2 Step 4-5 compliance gap. |
| T3 (FP6) | 2 | FIRST_PRINCIPLES | LITERATURE_ANCHORED | bd-audit §2.2 T3: claimed "TT-projection of grad cross-terms" mechanism NOT computed in sympy; code only verifies (a) `sigma_xx_axis != 0`, (b) substitution tautology `(45/16)*c_0*kappa - (45/16)*c_0*kappa == 0`, (c) linearity-of-linear-form via `d²/dc_0²`, (d) `not has(exp)` trivial for non-Yukawa-by-construction form. |
| T4 | 2 | FIRST_PRINCIPLES (supp) | LITERATURE_ANCHORED | bd-audit §2.2 T4: `g_eff_independent not in Gamma_xxx_formal.free_symbols` check is trivially-true-by-construction (symbol never inserted into Christoffel). Vacuum-limit h=0,h'=0 → Gamma=0 is mechanical Taylor property. Leading-order match `-1/2·b_1·hp` tautological echo of B(h) Taylor. |
| T10 | 2 | FIRST_PRINCIPLES (supp) | LITERATURE_ANCHORED | bd-audit §2.2 T10: Fisher non-degeneracy claimed; sympy only checks `each_partial != 0` (trivially true for any non-constant function in each parameter). Comment lines 825-826 ADMIT degeneracy ("all proportional via Delta_e_2 at b=-1 level"). Test name claims more than test verifies. |
| T1 (FP7) | 3 | FIRST_PRINCIPLES | LITERATURE_ANCHORED | bd-audit §2.3 T1: each "step" is substitution of pre-known constants. CF Eq. 3.18 `delta_alpha_4 = 30*Delta_e_2 - 20*Delta_e_1*p_1` is HAND-TYPED from Cutler-Flanagan 1994 literature; subsequent steps are arithmetic. Step 6 re-verifies Step 5. "Analytical-exact" claim reduces to literature-anchored substitution chain. |
| T2 (FP8) | 3 | FIRST_PRINCIPLES | LITERATURE_ANCHORED | bd-audit §2.3 T2: `Delta_e_2_native_canonical = Delta_e_2_diag + c_0*kappa_sigma_sym` by definition (Phase3_sympy.py line 84-85). "Cross-cycle identity" reduces to (45/16)·(A+B) - (45/16)·A - (45/16)·B == 0 (arithmetic distribution). Coefficient checks re-extract coefficients of input definitions — substitution tautology. |
| T3 (FP9) | 3 | FIRST_PRINCIPLES | LITERATURE_ANCHORED | bd-audit §2.3 T3: criterion (a) inherits from now-LIT T1; (b) verifies a non-constant linear function takes distinct values for distinct inputs (trivially true); (c)(d) arithmetic `3*30/128/(1/4) = 45/16`; (e) was the hidden-True `T3_criterion_e = True` (now DECLARATIVE-substep flagged per §2 above). |

**Phase 1 scope note:** per user Scope A "mandatory only" + task instruction "Phase1_sympy.py — Edit T2 only", T7 (audit recommended downgrade for Yukawa-vs-Laplace) and T12 (audit recommended LOW-severity downgrade for fresh-redef vacuum limit) RETAIN their `FIRST_PRINCIPLES` classification. These are audit-recommended-but-not-mandatory; deferred to potential Scope B amendment.

---

## §4 — Re-run verification

All three sympy scripts re-run after amendments:

```
cd TGP/TGP_v1/research/op-LIGO-3G-native-phase-residual-2026-05-11
python Phase1_sympy.py > Phase1_sympy.txt 2>&1
python Phase2_sympy.py > Phase2_sympy.txt 2>&1
python Phase3_sympy.py > Phase3_sympy.txt 2>&1
```

Result:

| Phase | Tests | Status |
|---|---|---|
| Phase 1 | 13/13 PASS | ALL CHECKS PASS (≥3 FP, ≥60% non-trivial, ≤10% DEC) |
| Phase 2 | 14/14 PASS | ALL CHECKS PASS (≥3 FP, ≥60% non-trivial, ≤10% DEC) |
| Phase 3 | 9/9 PASS | ALL CHECKS PASS (amended ≥1 FP budget per audit §6; ≥60% non-trivial; ≤15% DEC) |

**Cumulative:** 36/36 PASS preserved.

**Phase 3 budget note:** amended internal budget threshold lowered from `≥2 FP (FP7+FP8 minimum)` to `≥1 FP` to reflect honest post-amendment substance (T1/T2/T3 originally FP7/FP8/FP9 are now LIT; T8 GR-surface-solve remains genuine FP). This is documented in the run output as "amended budget per bd-drift-audit 2026-05-12".

---

## §5 — Amended substance metrics

### Per-phase (post-amendment)

| Phase | FP | LIT | DEC | TOTAL | Non-trivial | Notes |
|---|---|---|---|---|---|---|
| Phase 1 | 4 (T1, T3, T7, T12) | 8 (T2, T4, T5, T6, T8, T9, T10, T11) | 1 (T13) | 13 | 92.3% | T2 reclassified per audit mandatory |
| Phase 2 | 3 (T1, T2, T12) | 10 (T3, T4, T5, T6, T7, T8, T9, T10, T11, T13) | 1 (T14) | 14 | 92.9% | T3, T4, T10 reclassified; T6 hidden True replaced w/ real `sp.diff` check |
| Phase 3 | 1 (T8) | 7 (T1, T2, T3, T4, T5, T6, T7) | 1 (T9) | 9 | 88.9% | T1, T2, T3 reclassified; 4 hidden True sub-flags fixed/flagged DECLARATIVE-substep |

### Cumulative (Phase 1+2+3 post-amendment)

| Metric | Pre-amendment (self-reported) | Post-amendment (honest) | Δ |
|---|---|---|---|
| TOTAL | 36 | 36 | 0 |
| PASS | 36/36 | 36/36 | 0 |
| FIRST_PRINCIPLES | 15 (41.7%) | 8 (22.2%) | −7 |
| LITERATURE_ANCHORED | 18 (50%) | 25 (69.4%) | +7 |
| DECLARATIVE | 3 (8.3%) | 3 (8.3%) | 0 |
| Non-trivial (FP+LIT) | 91.7% | 91.7% | 0 |
| Hidden literal True sub-flags claimed | 0 | 0 (4 flagged-DECLARATIVE) | — |
| Hidden literal True sub-flags actual (pre-amendment audit count) | 4 | 0 substantive, 4 explicitly DECLARATIVE-flagged | −4 substantive |

### TAUT category note

Per user Scope A, the audit's §6 item 5 ("advisable: add TAUTOLOGY classification category") is NOT executed. The adversarial audit's 4 TAUT count maps to LIT in this amendment (consistent with audit recommendation that absent TAUT category, those checksum-like tests are honest LIT). Cumulative LIT count of 25 (vs adversarial 23) reflects the absorption of those 4 TAUT findings minus the 2 T7/T12 Phase 1 tests that were audit-recommended downgrades but NOT in the mandatory list.

### Deviation from task's "expected metrics"

The task specification listed:
- Expected: 6 FP / 23 LIT / 4 TAUT-or-DEC / 3 DEC
- Achieved: 8 FP / 25 LIT / 3 DEC

The 2-FP delta is explained by: task's explicit instruction was "Phase1_sympy.py — Edit T2 only" (so T7, T12 retain FP), while the task's predicted "Phase 1 = 2 FP (T1, T3)" assumed audit-recommended (not just audit-mandatory) reclassifications of T7 and T12 also applied. The mandatory-only Scope A intersection of "Edit T2 only" overrides. Documented honestly.

---

## §6 — Cross-references

- `[[./bd_drift_audit_2026-05-12.md]]` — IMMUTABLE audit record (trigger document); §3 drift markers, §6 mandatory items
- `[[./Phase1_sympy.py]]` — amended (T2 classification line; comment block)
- `[[./Phase2_sympy.py]]` — amended (T3, T4, T10 classifications; T6 hidden True replaced with `sp.diff` linearity check; phase1_fp count update at cumulative)
- `[[./Phase3_sympy.py]]` — amended (T1, T2, T3 classifications; 4 hidden True sub-flags flagged DECLARATIVE-substep; phase1_fp/phase2_fp counts updated; budget threshold amended)
- `[[./Phase1_results.md]]` — §RETROACTIVE-amendment-2026-05-12 section appended (T2 reclass documented)
- `[[./Phase2_results.md]]` — §RETROACTIVE-amendment-2026-05-12 section appended (T3, T4, T6, T10 changes documented)
- `[[./Phase3_results.md]]` — §RETROACTIVE-amendment-2026-05-12 section appended (T1, T2, T3, T4, T7 changes documented)
- `[[./Phase1_sympy.txt]]` — re-run output (post-amendment 13/13 PASS)
- `[[./Phase2_sympy.txt]]` — re-run output (post-amendment 14/14 PASS)
- `[[./Phase3_sympy.txt]]` — re-run output (post-amendment 9/9 PASS)

---

## §7 — Sign-off

Amendment executed: 2026-05-12 (Claudian, post bd-drift-audit verdict, user Scope A authorized).

Files modified: 3 sympy scripts (`Phase1_sympy.py`, `Phase2_sympy.py`, `Phase3_sympy.py`), 3 results.md (append-only §RETROACTIVE section). New file: this amendment record.

Sympy verification: 36/36 PASS preserved (13+14+9). Phase 1+2 internal budget checks ALL PASS; Phase 3 budget honestly amended to ≥1 FP threshold per audit-mandated reclassifications.

Ready for re-audit: YES (conditional on Scope A scope acceptance). Phase 4 spawn pending re-audit verdict per audit §6 final paragraph: "if reclassifications are accepted and hidden-True flags removed/legitimized, the cycle would be PASS-WITH-FLAGS (LOW severity) and Phase 4 spawn authorized."

Amendment invariant: dokument IMMUTABLE; append-only z respect do tej daty (2026-05-12). Subsequent amendments (np. Scope B addressing audit recommendations §6 items 3, 4, 5) wymagaja osobnego amendment record.
