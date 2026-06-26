---
title: "Phase 2 — F-PSR-B observational comparison + F-PSR-C cross-system"
type: phase_results
status: PHASE2_COMPLETE_AWAITING_USER_REVIEW
phase: 2
cycle: op-PSR-orbital-drift-2026-05-24
parent_motivation: "Phase 1 LOCKED Δz_poly(α, U) magnitude derivation"
created_date: 2026-05-24
authorization: "User 2026-05-24: 'działaj Phase 2'"
methodology_binding: "CALIBRATION_PROTOCOL.md §3.6 BINDING"
anti_lakatos_lock: PRESERVED
substantive_fp_total_phase2: 5
substantive_fp_pass_phase2: 5
substantive_fp_fail_phase2: 0
cumulative_fp_total: 11
cumulative_fp_pass: 9
cumulative_findings: 2  # Phase 1 gauge mismatches
hardcoded_T_pass_count: 0
partial_compute_count: 0
partial_compute_cumulative: "0/1"
dec_budget_used: "0/3"
falsifier_status:
  F-PSR-A: "PASS (Phase 1; magnitude derivation procedure)"
  F-PSR-B: "PASS_CONSISTENT (Phase 2; S07 α ⊆ NICER allowed range)"
  F-PSR-C: "PASS (Phase 2; cross-system independent)"
---

# Phase 2 — F-PSR-B + F-PSR-C observational verdict

## §1 — Executive summary

Phase 2 of B cycle tested TGP polynomial α-family against NICER NS measurements (J0030+0451, J0740+6620) and supplementary spectroscopic NS redshift (EXO 0748-676 Cottam+2002).

**Result:**
- **F-PSR-B: PASS_CONSISTENT** — S07-LOCKED α ∈ [-0.832, 0.832] is fully consistent with NICER bounds for both systems
- **F-PSR-C: PASS** — cross-system check independent; same α-range allowed by both systems

**Cumulative B-cycle (Phase 1 + Phase 2): 9 PASS + 2 gauge findings (Phase 1) = 11 substantive items, 0 hardcoded T_pass, 0 DEC, 0 PARTIAL_compute.**

---

## §2 — Substantive FP results (sympy executed)

**Script:** [[Phase2_sympy.py]] (executed; full output preserved).

| FP ID | Description | Status |
|-------|-------------|--------|
| T1 | F-PSR-B: S07 α range ⊆ NICER bounds (both systems) | **PASS** ✓ |
| T2 | FAIL_TINY check: |Δz_TGP_max|/σ_z vs 0.1 threshold | **PASS** ✓ (ratios 0.109-0.115, just ABOVE FAIL_TINY threshold) |
| T3 | F-PSR-C: cross-system α-range consistency | **PASS** ✓ |
| T4 | Cottam (EXO 0748-676) supplementary consistency | **PASS** ✓ |
| T5 | Cumulative B-cycle verdict synthesis | **PASS** ✓ |

**Summary:** 5/5 PASS Phase 2.

---

## §3 — NICER observational input (LITERATURE_ANCHORED)

### §3.1 — J0030+0451 (Miller et al. 2019)

| Quantity | Nominal | Asymmetric error |
|----------|---------|-------------------|
| M | 1.44 M_sun | +0.15 / -0.14 |
| R | 13.02 km | +1.24 / -1.06 |
| U = GM/(c²R) | **0.1633** | ±0.0219 (13.4% precision) |
| z_GR = 1/√(1-2U) - 1 | **0.2187** | ±0.0396 (18.1% precision) |

### §3.2 — J0740+6620 (Miller et al. 2021 NICER+XMM; Cromartie 2020 Shapiro)

| Quantity | Nominal | Asymmetric error |
|----------|---------|-------------------|
| M | 2.072 M_sun | ±0.066 |
| R | 13.7 km | +2.6 / -1.5 |
| U = GM/(c²R) | **0.2234** | ±0.0342 (15.3% precision) |
| z_GR = 1/√(1-2U) - 1 | **0.3444** | ±0.0830 (24.1% precision) |

### §3.3 — Cottam et al. 2002 (supplementary; disputed)

| Quantity | Value | Note |
|----------|-------|------|
| z (Fe XXVI/XXV) | 0.35 ± 0.05 | EXO 0748-676 XMM spectroscopy; **disputed by later analyses** |
| Inferred U (GR) | 0.226 | |

**Cottam treated as supplementary only**; not relied upon for primary verdict.

---

## §4 — F-PSR-B observational verdict

### §4.1 — Methodology

Per Phase 1 LOCKED:
$$\Delta z_{poly}(\alpha, U) = \alpha \cdot \left(\frac{U^2}{8} + \frac{5U^3}{16} + \frac{11U^4}{16}\right) + \mathcal{O}(U^5)$$

Per Phase 0 acceptance criterion: |Δz_TGP(α)| ≤ σ_z_obs.

Solve for α_max: $|\alpha|_{max} = \sigma_{z}^{obs} / |\Delta z_{poly}(\alpha=1, U)|$

### §4.2 — α-allowed range per system

| System | U | σ_z | |Δz_poly(α=1)| | |α|_max (NICER) |
|--------|---|-----|----------------|-----------------|
| J0030+0451 | 0.163 | 0.040 | 0.0052 | **7.64** |
| J0740+6620 | 0.223 | 0.083 | 0.0114 | **7.26** |
| Cottam EXO 0748 (sup.) | 0.226 | 0.050 | 0.0117 | 4.26 |

**S07-LOCKED:** |α| ≤ 0.832

**Comparison:**
- NICER J0030 |α|_max = 7.64 ≫ 0.832 (S07)
- NICER J0740 |α|_max = 7.26 ≫ 0.832 (S07)
- Cottam |α|_max = 4.26 > 0.832 (S07)

**All three observational bounds ALLOW the S07-LOCKED α range** with substantial margin.

### §4.3 — F-PSR-B verdict: **PASS_CONSISTENT**

S07 polynomial α-range is fully consistent with NS surface observations. No constraint tighter than S07.

### §4.4 — FAIL_TINY borderline check

Per Phase 0: FAIL_TINY if |Δz_TGP_max|/σ_z < 0.1.

| System | |Δz_TGP_max at α=0.832| | σ_z | Ratio |
|--------|------------------------|-----|-------|
| J0030 | 0.00432 | 0.0396 | **0.109** |
| J0740 | 0.00951 | 0.0830 | **0.115** |

Both ratios are **just above 0.1 threshold** → **NOT FAIL_TINY**.

Effective interpretation: signal is **marginal but not negligible** relative to current precision. Not a strong falsification test, but not "unobservable" either. With modest precision improvements (factor ~3), this would become a real test.

---

## §5 — F-PSR-C cross-system independence verdict

### §5.1 — Cross-system check methodology

Hypothesis: same α governs both J0030 and J0740 (linear-in-α scaling assumption).

Combined NICER bound (intersection of allowed ranges): **|α| ≤ 7.26** (limited by tighter system).

Predicted signals at S07 edge α = 0.832:
- J0030: Δz_pred = 0.0043 (vs σ_z = 0.040) → consistent
- J0740: Δz_pred = 0.0095 (vs σ_z = 0.083) → consistent
- Signal ratio J0740/J0030 = 2.20 (different U → different signal scale; expected)

### §5.2 — F-PSR-C verdict: **PASS**

Both systems consistent with same α range. No system-dependent α required. Linear-in-α scaling assumption preserved.

---

## §6 — Cumulative B-cycle metrics

| Metric | Phase 1 | Phase 2 | Cumulative |
|--------|---------|---------|------------|
| Substantive FPs | 6 | 5 | **11** |
| PASS | 4 | 5 | **9** |
| Gauge findings | 2 | 0 | 2 |
| FAIL | 0 | 0 | 0 |
| Hardcoded T_pass=True | 0 | 0 | **0 ✓** |
| PARTIAL_compute | 0 | 0 | **0/1 budget** |
| DEC | 0 | 0 | **0/3 budget** |
| Anti-Lakatos compliance | ✓ | ✓ | **COMPLIANT** |

**Verdict status:**
- F-PSR-A: **PASS** (Phase 1)
- F-PSR-B: **PASS_CONSISTENT** (Phase 2)
- F-PSR-C: **PASS** (Phase 2)

---

## §7 — Physical interpretation

### §7.1 — What this cycle established

TGP polynomial α-family (post-M9.1'' GWTC-3 falsification) is:
1. **CONSISTENT** with current NS surface observations (NICER + spectroscopy)
2. **NOT tightly constrained** by these — NICER systematic precision (15-24%) is wider than LIGO ppE (which gives the dominant α ≤ 0.832 bound)
3. **Marginal-signal regime**: predicted deviation (0.4-1.0%) is at the edge of 10% detection threshold

### §7.2 — What this cycle did NOT establish

This cycle did NOT:
- Falsify TGP polynomial family (signal too small relative to precision)
- Strongly constrain α beyond S07-LOCKED LIGO ppE bounds
- Test BH shadow / photon ring observables (those in PR-012 scope)
- Test binary pulsar orbital dynamics (U_orbit ~ 10⁻⁶ → O(U³) ~ 10⁻¹⁸, far below current precision)

### §7.3 — Future-test target

Future instruments could turn this into strong test:
- **NICER-Plus / eXTP**: σ_R / R → 1-2% → σ_z → 2-3% → |α|_max ~ 0.2-0.3 (tighter than S07)
- **SKA pulsar timing**: orbital observables to 10⁻⁹ level (would be sensitive to O(U³) for some Galactic Center pulsars near Sgr A*)
- **NS-NS merger spectroscopy** (post-detection era): direct surface redshift measurements

**Recommendation:** Register cycle B observable as **future-test target** in PR registry (PR-017 candidate); claim_status pending pre-observational pattern.

---

## §8 — Anti-Lakatos compliance check

| Item | Status |
|------|--------|
| F8 FAILs cited as motivation? | NO ✓ |
| F8_FORENSIC cited as positive evidence? | NO ✓ |
| E1/E2 explorations cited as predictions? | NO ✓ |
| Factor-10 threshold from γ-7? | NO ✓ (σ_obs used) |
| Phase 0 thresholds modified post-hoc? | NO ✓ |
| New falsifiers added in Phase 2? | NO ✓ |
| Threshold lowered to make PASS? | NO ✓ |
| 0 hardcoded T_pass=True | ✓ verified |
| Independent observational input only? | ✓ (Miller, Riley, Cromartie, Cottam — LITERATURE_ANCHORED) |

**Anti-Lakatos status:** COMPLIANT ✓

---

## §9 — Phase 2 status + Phase FINAL recommendation

| Field | Value |
|-------|-------|
| Phase 2 status | COMPLETE (awaiting user review) |
| F-PSR-A | PASS (Phase 1) |
| F-PSR-B | PASS_CONSISTENT (Phase 2) |
| F-PSR-C | PASS (Phase 2) |
| Cumulative FPs | 11 (9 PASS + 2 gauge findings) |
| Hardcoded T_pass | 0 |
| Anti-Lakatos | COMPLIANT |
| Recommended next | **Phase FINAL closure ceremony** + PR-017 registry entry |

**Anticipated claim_status (informational, not pre-registered):**

Per Phase 0 §3 (F-PSR-B): "FAIL_TINY (signal below precision): |Δ_TGP| < 0.1·σ_obs → consistent BUT not a strong test (registered as future-test target)".

Our T2 ratios (0.109, 0.115) are **just above** FAIL_TINY threshold (0.1), so technically NOT FAIL_TINY. But signal-to-precision ratio of ~0.11 indicates **weak falsification test**.

Recommended claim_status: **B+** (positive consistency, weak observational discrimination, future-test target registered).

This is analogous to PR-012 PR-010 "pre-observational pattern" — TGP consistent with current data; future improvements may discriminate.

---

## §10 — Files

- [[Phase2_sympy.py]] — sympy + numerical script (executed)
- [[Phase2_results.md]] — this report

**Awaiting user authorization** for Phase FINAL (aggregate closure + PR-017 entry).
