---
title: "Phase 3 — F-INT-C PR-020 LOCK candidate definition + threshold derivation"
type: phase_PR_candidate
phase: 3
status: PHASE3_COMPLETE
cycle: op-S07-emergent-metric-integration-2026-06-01
created_date: 2026-06-01
authorization: "User 2026-06-01: 'działaj Phase 3'"
audit_target: F-INT-C (Phase 0 §3 LOCKED)
f_int_c_verdict: "PASS_PARTIAL_HEURISTIC (PR-020 LOCK candidate fully specified; numerical anchors c_0=4π and κ_σ=1/(3π) HEURISTIC; rigorous pinning deferred to O1/O2 future cycles)"
sympy_phase3: "10/10 PASS (Phase3_sympy.py inheritance verification + threshold derivation)"
methodology_binding: "CALIBRATION_PROTOCOL.md §3.6 BINDING; Phase 0 §4 AUDIT + §4.5 predecessor verdict invariance LOCK"
anti_lakatos: "COMPLIANT — PR-020 LOCK candidate inherits LOCKED predecessor results (emergent-metric Phase 3+4, c_0/κ_σ joint LOCK); 0 hardcoded T_pass=True; no predecessor verdicts modified"
discipline:
  hardcoded_T_pass: "0/10 ✓"
  dec_used: "0/3 (cumulative across Phase 0-3)"
  partial_compute_used: "0/1"
  partial_concept_mismatch: "0"
documentation_observation: "c_0-derivation Phase FINAL §3.3 says 'β_ppE ≈ 0.08' but Phase 3 sympy confirms β_ppE^new at GW150914 calibration is 0.225 (the 0.08 is the c_0·κ_σ deviation from 4/3, NOT the resulting β_ppE value). This is INFORMATIONAL documentation observation — does NOT modify predecessor verdict; flagged as minor cleanup opportunity (add to R1 #21 or future Phase FINAL annotation pass)."
---

# Phase 3 — F-INT-C PR-020 LOCK candidate

## §0 — Executive verdict

**F-INT-C: PASS_PARTIAL_HEURISTIC**

PR-020 LOCK candidate is fully specified with all required attributes: (a) native observable + units, (b) TGP predicted value range, (c) cycle reference, (d) measurement instrument + timeline. Numerical anchors are HEURISTIC (c_0 = 4π geometric, κ_σ = 1/(3π) heuristic, joint c_0·κ_σ = 4/3 EXACT via clean π cancellation) per LOCKED predecessor cycles op-c0-derivation + op-kappa-sigma. Rigorous pinning is DEFERRED to future cycles O1 (κ_σ Hadamard 2-body PN) and O2 (c_0 covariant Path A→B), per F-INT-D inventory Phase 2 §2.1.2.

PR-020 is therefore classified **LOCKED-PR020-CONDITIONAL** in PRE_REGISTERED_FALSIFIERS.md format (analog to PR-018 STRUCTURAL_PARTIAL classification).

**Phase 3 sympy verification:** 10/10 PASS ([[Phase3_sympy.py]], [[Phase3_sympy.txt]]) — confirms all PR-020 numerical claims compute-then-compare against LOCKED predecessors with 0/10 hardcoded T_pass=True.

---

## §1 — PR-020 native observable definition

### §1.1 — Observable

**β_ppE^new at 2.5PN (b = −1) inspiral phase for BBH events at η = 1/4 (equal mass)**

Where β_ppE^new is the post-Einsteinian (ppE) framework parameter at 2.5PN order (b = −1 phase term) for the TGP emergent-metric Phase 4 Path 2 σ-coupling recovery. Equivalently, the inspiral phase residual relative to GR baseline:

```
Δφ_GW(f) = (3/(128η)) · (πMf)^(b) · β_ppE      (per ppE framework convention)
```

Native form per [[../op-emergent-metric-from-interaction-2026-05-09/Phase3_results.md]] LOCK:

```
β_ppE^new |_{η=1/4} = (45/16) · [Δe_2 + c_0 · κ_σ]

Δe_2 = -a_1·ξ_3 - 3 - 4·a_2/a_1² + 4·b_2/a_1² - 8·a_3/a_1³ + 16·a_2²/a_1⁴
```

with M9.1''-canonical PN coefs (Path 2 preserved):
- a_1 = 4, a_2 = 12, b_2 = 4, a_3 = 36, ξ_3 = 5/24
- Δe_2 = −4/3 (Phase 3 sympy FP2 verified)
- → β_ppE^M9.1'' = (45/16)·(−4/3) = **−15/4 = −3.75** (the FALSIFIED reference, Phase 3 sympy FP1)

With σ-coupling addition (Phase 4 Path 2 STRUCTURAL_PREFERRED):

```
β_ppE^new = -15/4 + (45/16) · c_0 · κ_σ
```

### §1.2 — Units

Dimensionless. β_ppE^new is a dimensionless residual in the ppE framework (phase deviation per GR template at given PN order).

---

## §2 — TGP predicted value range

### §2.1 — Geometric target (clean π cancellation)

Per joint LOCK of [[../op-c0-derivation-from-substrate-2026-05-09/Phase_FINAL_close.md]] + [[../op-kappa-sigma-2body-PN-2026-05-09/Phase_FINAL_close.md]]:

- c_0 = 4π (geometric, from OP-7 T3.4 ξ_eff = 4π·G·Φ_0² → Path A→B conversion)
- κ_σ = 1/(3π) (heuristic, from 1/π orbital averaging × 1/3 σ-trace structure)
- **c_0 · κ_σ = 4π · 1/(3π) = 4/3 EXACT** (clean π cancellation, two independent sources, Phase 3 sympy FP4)

At geometric target:

```
β_ppE^new(geometric) = -15/4 + (45/16)·(4/3) = -15/4 + 15/4 = 0 EXACT
```

(Phase 3 sympy FP3 verified.)

### §2.2 — GW150914 calibrated estimate

Per c_0-derivation Phase FINAL §1: "c_0 (z GW150914 ξ/G ≈ 1.06 calibration): c_0 ≈ 4π·1.06 ≈ 13.32".

At GW150914 calibration:
- c_0·κ_σ = 4π·1.06 · 1/(3π) = 4·1.06/3 ≈ **1.4133**
- β_ppE^new(GW150914) = (15/4)·(c_0·κ_σ − 4/3) = (15/4)·0.06 = **0.225** ≈ **+0.22**

(Phase 3 sympy FP5 verified.)

**Documentation observation:** c_0-derivation Phase FINAL §3.3 states "β_ppE ≈ 0.08 within GWTC-3 bound 0.78". This is INACCURATE — the 0.08 is the deviation of c_0·κ_σ from 4/3 (= 1.413 − 1.333), NOT the resulting β_ppE^new value. Phase 3 sympy FP5 confirms β_ppE^new ≈ 0.225 at GW150914 calibration. **This is INFORMATIONAL only**, does NOT modify predecessor verdict per Phase 0 §4.5 LOCK; flagged as minor doc cleanup candidate.

### §2.3 — Pre-registered TGP value range

**PR-020 TGP value:**
- **Geometric central value:** β_ppE^new = 0 EXACT (at c_0·κ_σ = 4/3)
- **Heuristic uncertainty:** ± 0.225 (from GW150914 calibration spread; reflects gap between geometric prediction and observation-fit calibration)
- **Range:** β_ppE^new ∈ [−0.225, +0.225] under current heuristic c_0/κ_σ pinning

**CONDITIONAL caveat:** Range tightens if future O1/O2 cycles deliver rigorous c_0/κ_σ that constrain calibration spread.

---

## §3 — Threshold derivation (current + future tightening)

### §3.1 — Current bound: GWTC-3

Per [[../op-GWTC3-reanalysis/Phase2_RERUN_2026-05-09_corrected_beta.md]] + [[../op-emergent-metric-from-interaction-2026-05-09/Phase4_results.md]] §1:

- **GWTC-3 1σ bound: |β_ppE| ≤ 0.78** (combined ~90 BBH posterior)
- Implied c_0·κ_σ window: **[1.0560, 1.6107]** (width 0.555)
- Phase 3 sympy FP6 verified.

**Geometric target 4/3 = 1.333 is INSIDE GWTC-3 window** (FP7 verified). **GW150914 calibrated 1.413 also INSIDE** (FP8 verified). **Recovery framework is currently observationally consistent.**

### §3.2 — Future tightening: ET-D / CE / LISA projections

Per Phase 0 §2 mandatory reading inheritance + standard literature projections for 3G GW detectors:

| Instrument | Timeline | Projected sensitivity tightening | Implied \|β_ppE\| bound | Implied c_0·κ_σ window |
|------------|----------|----------------------------------|--------------------------|-------------------------|
| Current LIGO-Virgo (O4) | ~2024-2027 | ~1× | 0.78 | [1.056, 1.611] |
| LIGO O5 + KAGRA | ~2028-2030 | ~2× tighter | ~0.39 | [1.195, 1.473] |
| Einstein Telescope D (ET-D) | ~2035+ | ~10× tighter (BBH at z<5) | ~0.078 | **[1.306, 1.361]** (width 0.056) |
| Cosmic Explorer (CE) | ~2035+ | similar to ET-D for BBH | ~0.078 | similar |
| LISA | ~2035+ | SMBH inspirals, complementary | varied by mass range | requires separate ppE mapping |

(Phase 3 sympy FP9 verified ET-D implied window.)

### §3.3 — Critical falsification gate at ET-D precision

At ET-D bound |β_ppE| ≤ 0.078:

- **Geometric target c_0·κ_σ = 4/3 = 1.333: INSIDE** (FP9 verified)
- **GW150914 calibrated c_0·κ_σ = 1.413: OUTSIDE ET-D window** (FP9 verified)

**This is the heart of PR-020:**
- IF future GW data confirms β_ppE^new ≈ 0 to ET-D precision (|β_ppE| ≲ 0.078), recovery framework is **VALIDATED** at the geometric clean-π-cancellation prediction
- IF future GW data places β_ppE^new at ~0.22 (the GW150914 calibration central), recovery framework would face **TENSION** between geometric prediction and observation
- IF future GW data excludes β_ppE^new = 0 at >5σ for any value, recovery framework is **FALSIFIED**

---

## §4 — Falsification criteria (LOCKED)

PR-020 pre-registered decision rules (immutable post-LOCK per Phase 0 §6 forbidden moves):

| Verdict | Condition | Interpretation |
|---------|-----------|----------------|
| **SOFT_PASS** (current state) | GWTC-3 \|β_ppE\| ≤ 0.78 includes 0 | Recovery consistent within current observational precision |
| **PASS_NARROW_GEOMETRIC** | Future bound (e.g., ET-D) narrows to \|β_ppE\| ≲ 0.078 AND TGP value confirmed at β_ppE^new ≈ 0 | Geometric clean-π-cancellation prediction validated; rigorous c_0/κ_σ pinning vindicated post O1/O2 cycles |
| **PASS_NARROW_CALIBRATED** | Future bound narrows AND TGP value confirmed near 0.22 (GW150914-style calibration regime) | Recovery framework consistent but rigorous c_0/κ_σ require non-geometric values; O1/O2 cycles must re-pin |
| **TENSION** | Future bound narrows to 0.078 ≲ \|β_ppE\| ≲ 0.78 AND TGP value near 0.22 | Calibrated regime survives; geometric regime falsified; framework requires explanation of c_0·κ_σ ≠ 4/3 EXACT |
| **HARD_FAIL** | Future bound excludes 0 at >5σ AND TGP value pinned outside compatible range | **Recovery framework FALSIFIED at ET-D/CE/LISA precision** |

**Anti-Lakatos LOCK:** thresholds 0.78 (current), 0.078 (ET-D projected) IMMUTABLE post Phase 3 LOCK. No post-hoc loosening allowed.

---

## §5 — Conditional caveats (PR-020-CONDITIONAL)

PR-020 is classified **LOCKED-PR020-CONDITIONAL** with explicit dependencies on O1/O2 future cycles:

### §5.1 — Heuristic numerical anchors

- **c_0 = 4π** is HEURISTIC (geometric clean from OP-7 T3.4 Path A→B conversion; multi-session covariant derivation deferred per emergent-metric Phase FINAL §7 O2)
- **κ_σ = 1/(3π)** is HEURISTIC (1/π × 1/3 structural plausibility argument; explicit Hadamard 2-body PN regularization deferred per emergent-metric Phase FINAL §7 O1)
- **Joint c_0·κ_σ = 4/3 EXACT** is rigorous in the heuristic regime (clean π cancellation, two independent sources LOCKED) but ANCHORED in heuristic individual values

### §5.2 — What rigorous pinning would do

If O1 (κ_σ Hadamard rigorous) and O2 (c_0 covariant rigorous) deliver values that maintain joint c_0·κ_σ = 4/3 EXACT structurally, PR-020 promotes to **LOCKED-PR020-RIGOROUS**. If joint product shifts (e.g., to 4/3 × correction factor), PR-020 TGP value range re-pins; PR-020 falsification criteria thresholds remain LOCKED (anti-Lakatos).

### §5.3 — Independence from heuristic-vs-rigorous status

Phase 3 sympy verifies: **all PR-020 falsification thresholds are inherited from observational instruments (GWTC-3 0.78, ET-D 0.078), NOT from heuristic c_0/κ_σ values.** Therefore, threshold structure is robust to future rigorous pinning; only the TGP central value range may shift.

This is the precise meaning of PASS_PARTIAL_HEURISTIC per Phase 0 §3: "Observational target identified but numerical value depends on heuristic c_0/κ_σ; pre-register as PR-020-CONDITIONAL with explicit caveat".

---

## §6 — Compatibility check with PR-010 (Phase 0 §8 R-INT-7)

PR-010 LOCKED (per [[../../meta/PRE_REGISTERED_FALSIFIERS.md]]): α ∈ [−0.832, 0.832] for S07 polynomial family.

PR-010 and PR-020 are **DIFFERENT parameterizations of the SAME observable category** (PN coefficients in inspiral phase):
- **PR-010:** S07 polynomial family α (single-coefficient parametrization of M9.1''-class f(ψ) polynomial expansion)
- **PR-020:** emergent-metric 3-function ansatz parameters {a_n, b_n, c_n, c_0·κ_σ}

Per Phase 3 sympy FP10:
- No direct contradiction between PR-010 and PR-020
- PR-010 is bounded at GWTC-3 1σ (|α| ≤ 0.832, equivalent to |β_ppE| ≲ 0.78 in mapping)
- PR-020 specifies **DIFFERENT precision regime** (ET-D/CE/LISA ~2035+ projections) and **DIFFERENT parameterization** (3-function emergent-metric vs S07 polynomial)

**Compatibility verdict (FP10 PASS):** PR-010 and PR-020 coexist as complementary pre-registrations. PR-010 bounds the FALSIFIED S07-polynomial framework remnant; PR-020 bounds the RECOVERY emergent-metric framework prediction at future precision.

Per Phase 0 §3 F-INT-C: "PR-020 should target DIFFERENT observable or DIFFERENT precision regime than PR-010" — satisfied.

---

## §7 — Anti-Lakatos verification (Phase 3)

| Check | Status |
|-------|--------|
| PR-020 LOCK candidate inherits LOCKED predecessors (emergent-metric Phase 3+4, c_0/κ_σ joint LOCK)? | YES ✓ |
| No predecessor verdicts modified? | YES ✓ (Phase 0 §4.5 LOCK preserved) |
| 0 hardcoded T_pass=True in Phase 3 sympy? | YES ✓ (10/10 compute-then-compare) |
| Heuristic c_0/κ_σ NOT auto-promoted to rigorous? | YES — explicitly LOCKED-PR020-CONDITIONAL §5 ✓ |
| Thresholds (0.78, 0.078) inherited from observational instruments, NOT TGP fit? | YES — GWTC-3 from LIGO collaboration; ET-D from instrument projection ✓ |
| PR-020 NOT framed as F8 work? | YES — gravity sector framework, NOT F8 mechanism ✓ |
| Cycle A/B/D verdicts unchanged? | YES ✓ |
| F8 cycle verdicts unchanged? | YES ✓ |
| PR-010 unchanged? | YES (PR-010 LOCKED preserved; PR-020 complementary parameterization) ✓ |
| Falsification criteria pre-LOCKED IMMUTABLE? | YES — §4 explicit ✓ |
| Documentation observation (c_0-derivation §3.3 β_ppE typo) modified predecessor? | NO — INFORMATIONAL flag only; predecessor verdict UNCHANGED ✓ |
| Anticipated PR-020-CONDITIONAL outcome matches Phase 0 §13 prediction? | YES ✓ |

**Anti-Lakatos status:** COMPLIANT ✓

---

## §8 — PR-020 LOCK candidate full text (for PRE_REGISTERED_FALSIFIERS.md append)

The following is the proposed PR-020 entry for [[../../meta/PRE_REGISTERED_FALSIFIERS.md]] append (Phase FINAL of THIS cycle will commit):

---

### PR-020 (LOCKED 2026-06-01): TGP emergent-metric Phase 4 Path 2 β_ppE^new 2.5PN binary inspiral falsifier

- **Cycle:** [[../research/op-S07-emergent-metric-integration-2026-06-01/]] (Phase 0 + Phase 1 audit + Phase 2 supersession + Phase 3 PR-020 candidate definition + Phase FINAL, sesja #10 sprint 2026-06-01)
- **Pre-registration date:** 2026-06-01 (Phase 0 LOCKED at user "ok działaj z op-S07-emergent-metric-integration-cycle")
- **Pre-registration commit:** <git SHA to be inscribed at PR-020 LOCK commit; PRE_REGISTERED_FALSIFIERS.md + Phase_FINAL_close.md + STATE.md + S07 cycle status annotation update scheduled as single PR-020 activation commit>
- **Native observable:** β_ppE^new at 2.5PN (b = −1) inspiral phase residual for BBH events at η = 1/4 (equal-mass binary), per ppE framework convention. Native form: β_ppE^new |_{η=1/4} = (45/16)·[Δe_2 + c_0·κ_σ], with Δe_2 = −a_1·ξ_3 − 3 − 4a_2/a_1² + 4b_2/a_1² − 8a_3/a_1³ + 16a_2²/a_1⁴. At M9.1''-canonical PN coefs (Path 2 preserved: a_1=4, a_2=12, b_2=4, a_3=36, ξ_3=5/24): β_ppE^new = −15/4 + (45/16)·c_0·κ_σ.
- **TGP predicted value (geometric central):** β_ppE^new = **0 EXACT** at joint c_0·κ_σ = 4/3 (clean π cancellation from c_0 = 4π geometric × κ_σ = 1/(3π) heuristic, two independent LOCKED predecessor cycles)
- **TGP predicted value (GW150914-calibrated):** β_ppE^new ≈ +0.225 (with c_0 = 4π·1.06 from GW150914 ξ/G calibration, c_0·κ_σ ≈ 1.413)
- **TGP value range (current heuristic uncertainty):** β_ppE^new ∈ [−0.225, +0.225] under current c_0/κ_σ heuristic pinning
- **Decision rules (LOCKED, verbatim Phase 3 §4):**
  - **SOFT_PASS (current state):** GWTC-3 1σ \|β_ppE\| ≤ 0.78 includes 0 — recovery consistent within current observational precision
  - **PASS_NARROW_GEOMETRIC:** future bound (e.g., ET-D ~2035 \|β_ppE\| ≲ 0.078) narrows AND TGP value confirmed at β_ppE^new ≈ 0 — geometric clean-π-cancellation prediction validated; rigorous c_0/κ_σ pinning vindicated post O1/O2 future cycles
  - **PASS_NARROW_CALIBRATED:** future bound narrows AND TGP value confirmed near 0.22 (GW150914-style calibration regime) — recovery framework consistent but rigorous c_0/κ_σ require non-geometric values; O1/O2 cycles must re-pin
  - **TENSION:** future bound narrows to 0.078 ≲ \|β_ppE\| ≲ 0.78 AND TGP value near 0.22 — calibrated regime survives; geometric regime falsified; framework requires explanation of c_0·κ_σ ≠ 4/3 EXACT
  - **HARD_FAIL:** future bound excludes 0 at >5σ AND TGP value pinned outside compatible range — recovery framework FALSIFIED at ET-D/CE/LISA precision
- **Falsification target:** TGP emergent-metric Phase 4 Path 2 σ-coupling recovery framework (alternative to M9.1''-canonical 2.5PN prediction FALSIFIED at 5.02σ by GWTC-3 per M911-P1). PR-020 pre-registers the RECOVERY prediction at zero-β geometric target.
- **Confidence threshold:** Phase 3 LOCKED current GWTC-3 1σ (\|β_ppE\| ≤ 0.78); future ET-D/CE/LISA projection \|β_ppE\| ≲ 0.078 (factor 10 tightening) — both inherited from observational instrument projections, NOT TGP fit; IMMUTABLE post-LOCK.
- **Recovery scope (LOCKED, anti-Lakatos per Phase 0 §6 forbidden moves):**
  ```yaml
  allowed_directions:
    - "Future O1 cycle (op-kappa-sigma-Hadamard-rigorous-...) — rigorous κ_σ Hadamard 2-body PN regularization, 3-5 sesji estimated; may re-pin κ_σ value (impacts c_0·κ_σ joint product)"
    - "Future O2 cycle (op-c0-covariant-PathA-PathB-rigorous-...) — rigorous c_0 covariant Path A→B derivation, 3-5 sesji; may re-pin c_0 value"
    - "Future O3 research program (mechanism v for P6 R5 risk) — framework extension to resolve LIGO scalar mode amplitude with m_Φ ~ M_Pl regime; multi-session"
    - "Future GW data updates (ET-D / CE / LISA / 3G observation runs) — observational tightening per LOCKED thresholds"
    - "Joint c_0·κ_σ structural identity preservation — clean π cancellation rigorously re-derived post-O1+O2"
  forbidden_directions:
    - "Loosen GWTC-3 1σ threshold \|β_ppE\| ≤ 0.78 to accommodate larger TGP value (threshold inheritance LOCKED)"
    - "Loosen ET-D projected threshold \|β_ppE\| ≲ 0.078 (projection inheritance LOCKED)"
    - "Re-frame HARD_FAIL as 'partial success' or 'methodology issue' if triggered (verdict tampering)"
    - "Cite F8 cycles (γ-3/3'/5/7) as motivation for any future PR-020-related cycle"
    - "Cite cycle A (LAM PR-018) or cycle D (G PR-019) as evidence for or against PR-020"
    - "Modify predecessor verdicts (emergent-metric STRUCTURAL DERIVED, c_0/κ_σ STRUCTURAL DERIVED heuristic, S07 STRUCTURAL_CONDITIONAL_HALT, σ-3PN STRUCTURAL DERIVED) based on PR-020 outcome"
    - "Auto-promote heuristic c_0/κ_σ to rigorous DERIVED without O1/O2 cycle execution"
    - "Frame PR-020 as a 'predicted observation' rather than pre-registered falsifier (epistemic framing distinction)"
    - "Post-hoc adjust TGP value range to fit any newly-released GW data without explicit re-pinning rationale"
    - "Modify PR-010 (S07 polynomial α bound) based on PR-020 outcome (different parametrization)"
    - "Frame PR-020 outcome as F8 rescue OR as cycle A upgrade trigger (orthogonal scopes)"
  if_recovery_exhausted: "Cycle CLOSED-RESOLVED with HARD_FAIL verdict. TGP emergent-metric Phase 4 Path 2 σ-coupling recovery framework FALSIFIED at ET-D/CE/LISA precision. Recovery framework path forward: (a) re-evaluate Phase 4 Path 1 (3PN parameter tuning, alternative not yet excluded but Path 2 preferred for SU(2) consistency); (b) framework extension via mechanism v (O3); (c) accept TGP gravity sector framework requires further substantive revision beyond emergent-metric scope. F8 status UNCHANGED. PREDICTIONS_REGISTRY M911-P* status DOWNGRADE to PERMANENT-FALSIFIED (no recovery path validated)."
  ```
- **Status:** **LOCKED-PR020-CONDITIONAL** (PR-020 LOCK candidate fully specified at Phase 3 sympy 10/10 PASS; F-INT-C PASS_PARTIAL_HEURISTIC verdict; numerical anchors c_0=4π, κ_σ=1/(3π) HEURISTIC; rigorous pinning DEFERRED to O1/O2 future cycles; threshold structure ROBUST to rigorous pinning regardless of c_0/κ_σ re-evaluation)
- **Phase 3 sympy verification:** [[../research/op-S07-emergent-metric-integration-2026-06-01/Phase3_sympy.py]] **10/10 PASS** (compute-then-compare against LOCKED predecessor results; 0/10 hardcoded T_pass=True). Key verified values:
  - β_ppE^M9.1''-canonical = −15/4 = −3.75 (FALSIFIED reference)
  - β_ppE^new(c_0·κ_σ = 4/3 geometric) = 0 EXACT
  - β_ppE^new(GW150914 calibrated) = +0.225
  - GWTC-3 1σ window on c_0·κ_σ: [1.0560, 1.6107] (width 0.555); both geometric 4/3 and calibrated 1.413 INSIDE
  - ET-D projected window: [1.306, 1.361] (width 0.056); geometric INSIDE, GW150914-calibrated OUTSIDE — falsification gate active at ET-D precision
  - PR-010 / PR-020 cross-parameterization compatibility: NO direct contradiction (different parameterizations of PN observables; PR-010 S07-polynomial α; PR-020 emergent-metric c_0·κ_σ; both within GWTC-3 1σ)
- **Notes:**
  - **Independence from cycle D PR-019** (γ-derivation HONEST_NEGATIVE) explicitly declared (Phase 0 §1.4 + Phase 2 §1.4): PR-020 is gravity-sector framework falsifier; cycle D's γ epistemic status is orthogonal foundational scale question. Whether γ is calibrated or rigorously derived does NOT modify PR-020 structure or thresholds.
  - **Independence from cycle A PR-018** (Λ_eff STRUCTURAL_PARTIAL C+) explicitly declared: PR-020 is gravity sector inspiral phase; PR-018 is vacuum cosmological constant. Different observable categories.
  - **Independence from F8 cycles** explicitly declared: PR-020 is gravity sector framework; F8 cycles tested dark energy acceleration mechanisms (kinematic, c-variation, clumping, vacuum stress-energy). Different observable categories.
  - **Independence from PR-017 PSR** (NS surface O(U³) gravitational redshift): different observable (orbital drift / NS surface vs BBH inspiral phase).
  - **Compatibility with PR-010** (S07 polynomial α bound): NO direct contradiction; PR-020 is emergent-metric parameterization while PR-010 is S07-polynomial parameterization — both within GWTC-3 1σ for current observations; PR-020 specifies future-tightening regime (ET-D/CE/LISA ~2035+).
  - **R1 #21 reference:** TGP_FOUNDATIONS §3.6 documentation drift (3 annotation-level gaps) noted in Phase 1 audit; cleanup deferred to Phase FINAL of this cycle.
  - **Documentation observation:** c_0-derivation Phase FINAL §3.3 states "β_ppE ≈ 0.08" at GW150914 calibration; Phase 3 sympy FP5 confirms β_ppE^new = +0.225 (the 0.08 is the c_0·κ_σ deviation from 4/3, not the β_ppE value). INFORMATIONAL flag; does NOT modify predecessor verdict per Phase 0 §4.5 LOCK.
  - **Cross-cycle inheritance LOCKs (LEGITIMATE only):** emergent-metric Phase 3+4 LOCK (β_ppE^new formula + GWTC-3 window); c_0-derivation Phase 1 LOCK (c_0 = 4π heuristic); κ_σ Phase 1 LOCK (κ_σ = 1/(3π) heuristic); joint c_0·κ_σ = 4/3 EXACT (clean π cancellation); GWTC-3 from LIGO/Virgo collaboration observational anchor; ET-D/CE/LISA projections from instrument literature (LITERATURE_ANCHORED).
  - **Substance ceiling:** LOCKED-PR020-CONDITIONAL per pre-registered heuristic uncertainty disclosure. NOT promoted to LOCKED-PR020-RIGOROUS pending O1+O2 future cycles. Analog PR-018 STRUCTURAL_PARTIAL classification (sign + structure PASS, magnitude FAIL_LOW with rigorous pinning needed) — both pre-register honest classification matching actual derivation status.
  - **Estimated 1 sesja FINAL LOCK (this entry; Phase 3 sub-phase within sesja #10 sprint).**

---

## §9 — F-INT-C verdict

**PASS_PARTIAL_HEURISTIC** per Phase 0 §3 F-INT-C acceptance criterion:

> "PASS_PARTIAL_HEURISTIC: Observational target identified but numerical value depends on heuristic c_0/κ_σ (4π / 1/(3π)); pre-register as PR-020-CONDITIONAL with explicit caveat (PASS_HEURISTIC contingent on O1/O2 future cycle rigorous pinning)"

**Justification chain:**
1. PR-020 native observable identified: β_ppE^new at 2.5PN inspiral phase for BBH (§1)
2. TGP predicted value range fully specified: geometric 0 ± 0.225 heuristic (§2)
3. Threshold inherited from observational instruments: GWTC-3 0.78 (current), ET-D 0.078 (~2035) (§3)
4. Falsification criteria pre-LOCKED with 5 verdict categories (SOFT_PASS / PASS_NARROW_GEOMETRIC / PASS_NARROW_CALIBRATED / TENSION / HARD_FAIL) (§4)
5. Conditional caveats explicit: heuristic c_0/κ_σ pinning, rigorous deferred to O1/O2 (§5)
6. Compatibility with PR-010 verified (§6)
7. Phase 3 sympy 10/10 PASS compute-then-compare against LOCKED predecessors
8. Anti-Lakatos 12/12 checks COMPLIANT (§7)

**This is PASS, not FAIL_NO_LOCKBOX:** PR-020 IS a falsifiable lockbox — at ET-D/CE/LISA ~2035+ precision, geometric (β=0) vs GW150914-calibrated (β≈0.22) are distinguishable.

**This is PARTIAL_HEURISTIC, not PASS_PR020_LOCK_FULL:** Numerical anchors (c_0, κ_σ individually) are heuristic; promotion to LOCK_RIGOROUS requires O1+O2 future cycles to validate joint c_0·κ_σ = 4/3 EXACT structurally beyond clean-π-cancellation argument.

---

## §10 — Phase 3 statistics

```
Phase 3 sympy:                  10/10 PASS ✓
Hardcoded T_pass=True:           0/10
DEC used:                        0/3 (cumulative across Phase 0-3)
PARTIAL_compute used:            0/1 (cumulative)
PARTIAL_concept_mismatch:        0 (cumulative)
R1 raised in Phase 3:            0 (R1 #21 from Phase 1 unchanged)
Anti-Lakatos checks:            12/12 COMPLIANT ✓

PR-020 LOCK candidate:           FULLY SPECIFIED
Status:                          LOCKED-PR020-CONDITIONAL
Heuristic dependencies:          c_0 = 4π, κ_σ = 1/(3π); joint c_0·κ_σ = 4/3 EXACT
Rigorous pinning:                deferred to O1 + O2 future cycles
Threshold structure robust:      YES (observational anchors, not TGP fit)
Documentation observation:       c_0-derivation Phase FINAL §3.3 minor typo (β_ppE ≈ 0.08 vs actual 0.225); INFORMATIONAL only

Aggregate Phase 0 → Phase 3:
  F-INT-A (Phase 1):  PASS_WITH_ANNOTATIONS
  F-INT-B (Phase 2):  PASS_FULL_SUPERSESSION
  F-INT-D (Phase 2):  PASS_INVENTORY (4 outstanding items)
  F-INT-C (Phase 3):  PASS_PARTIAL_HEURISTIC ⭐
  → ALL 4 FALSIFIERS RESOLVED. Cycle ready for Phase FINAL.
```

---

## §11 — Recommendation: proceed to Phase FINAL

Per Phase 0 §10 Phase plan: **Phase FINAL = aggregate verdict + claim_status + PR-020 LOCK entry append to PRE_REGISTERED_FALSIFIERS.md + S07 status flip (per F-INT-B PASS_FULL_SUPERSESSION) + folder status updates + STATE.md sesja closure.**

Estymacja: 0.5 sesji.

**Phase FINAL deliverables anticipated:**
1. `Phase_FINAL_close.md` — aggregate closure ceremony
2. `meta/PRE_REGISTERED_FALSIFIERS.md` append: PR-020 LOCK entry (per §8 above)
3. `research/op-S07-alternative-f-psi-derivation-2026-05-09/README.md` + `Phase_FINAL_close.md` — supersession annotation update (per F-INT-B PASS_FULL_SUPERSESSION; classification annotation NOT verdict modification)
4. `TGP_FOUNDATIONS.md` §3.6.9 + §3.6.10.6 — annotation cleanups CL-1 + CL-2 (per F-INT-D inventory + R1 #21)
5. `research/op-S07-emergent-metric-integration-2026-06-01/README.md` — folder_status flip ACTIVE → CLOSED-RESOLVED
6. `STATE.md` — sesja #10 closure entry

**Cumulative cycle status after Phase FINAL:**
- claim_status: **CLOSED-RESOLVED INTEGRATION_COMPLETE** (4/4 falsifiers PASS or PASS-with-qualification; PR-020 LOCKED-CONDITIONAL; S07 supersession annotation applied; concept paper integration confirmed substantively complete with annotation cleanups)
- New R1: R1 #21 (TGP_FOUNDATIONS §3.6 doc drift) — addressed by Phase FINAL CL-1+CL-2
- Predecessor verdicts: ALL PRESERVED per §4.5 LOCK
- F8 status: UNCHANGED
- Publication-readiness: OUT OF SCOPE for this cycle; separate user decision

**Await user authorization** for Phase FINAL.

---

## §12 — Cross-references

- [[Phase0_balance.md]] — pre-registration LOCKED 2026-06-01
- [[Phase1_audit.md]] — F-INT-A PASS_WITH_ANNOTATIONS
- [[Phase2_supersession.md]] — F-INT-B PASS_FULL_SUPERSESSION + F-INT-D PASS_INVENTORY
- [[Phase3_sympy.py]] + [[Phase3_sympy.txt]] — 10/10 sympy verification
- [[../op-emergent-metric-from-interaction-2026-05-09/Phase3_results.md]] — β_ppE^new formula LOCK source
- [[../op-emergent-metric-from-interaction-2026-05-09/Phase4_results.md]] — GWTC-3 compliance window LOCK source
- [[../op-c0-derivation-from-substrate-2026-05-09/Phase_FINAL_close.md]] §1 — c_0 = 4π geometric heuristic LOCK
- [[../op-kappa-sigma-2body-PN-2026-05-09/Phase_FINAL_close.md]] — κ_σ = 1/(3π) heuristic LOCK
- [[../../meta/PRE_REGISTERED_FALSIFIERS.md]] — PR-020 LOCK candidate target (Phase FINAL will append)
- [[../../meta/PRE_REGISTERED_FALSIFIERS.md]] PR-010 entry — S07 polynomial α bound (compatibility check §6)
- [[../../STATE.md]] — sesja #10 entry (Phase 3 entry will follow)
