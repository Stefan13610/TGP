---
title: "Phase 0 — Balance sheet + 8/8 ☑ gate + scope (BINDING)"
date: 2026-05-16
parent: "[[./README.md]]"
type: phase-balance
phase: 0
status: 🟢 BALANCE OK — proceed to Phase 1 (8/8 ☑ verified)
---

# Phase 0 — Balance sheet

Per `meta/CALIBRATION_PROTOCOL.md` §2 (BINDING template for any DERIVED claim aspiration).

## §1 — External inputs (PDG 2024, observational)

| Symbol | Wartość | Precision | Źródło |
|---|---|---|---|
| m_u | 2.16 MeV | ±0.07 (MS-bar at 2 GeV) | PDG 2024 |
| m_d | 4.67 MeV | ±0.07 (MS-bar at 2 GeV) | PDG 2024 |
| m_s | 93.4 MeV | ±8.6 (MS-bar at 2 GeV) | PDG 2024 |
| m_c | 1.270 GeV | ±0.002 (m_c(m_c)) | PDG 2024 |
| m_b | 4.18 GeV | ±0.03 (m_b(m_b)) | PDG 2024 |
| m_t | 172.69 GeV | ±0.30 (pole) | PDG 2024 |
| m_e | 0.5110 MeV | ±0.0000003 (exact) | PDG 2024 |
| m_μ/m_e | 206.7682830 | 8e-7 | PDG 2024 (lepton ratio) |
| m_τ/m_e | 3477.23 | inherited why_n3 | PDG 2024 |

**Derived test ratios (5 niezależnych):**
| Ratio | Wartość | Źródło |
|---|---|---|
| m_c/m_u | 587.96 | PDG arithmetic |
| m_b/m_d | 894.97 | PDG arithmetic |
| m_t/m_c | 135.98 | PDG arithmetic |
| m_s/m_d | 20.00 | PDG arithmetic |
| m_b/m_t | 0.02421 | PDG arithmetic |

## §2 — Structural axioms (TGP-internal LOCKED)

| Anchor | Wartość | Independent LOCK source |
|---|---|---|
| **Universal mass formula** | m_obs = c_M · A_tail² · g_0^(e²/2) | why_n3 Phase 5 closure 2026-05-01 (CLOSED) |
| **β(α=2) = e²/2** | 3.6945 | op-L08-Phase6-e2-derivation-2026-05-16 B+ (LIVE) |
| **k_obs(α=2, d=3) = 3** | 3 | op-L05-mass-exponent-k-alpha-d-2026-05-16 A− (LIVE) |
| **m_obs ≠ M_full** distinction | operational | op-L05 A− (LIVE) |
| **Lepton calibration** | (e: g=0.869, A=0.110); (μ: g=1.407, A=0.650); (τ: g=1.755, A=1.666) | why_n3 PHASE1_psi_g0 + PHASE2_n_alpha (CLOSED 2026-05-01) |
| **Audit quark g_0 range** | [0.817, 0.891] | core/sek08b_ghost_resolution lin. 529 |
| **e_Euler** | 2.71828... | mathematical constant |
| **α=2 canonical** | 2 | op-L04-ODE-canonicalization-2026-05-04 LOCKED |

## §3 — Derived outputs (this cycle aspires to test)

| Output | Status | Phase 1 source |
|---|---|---|
| 5 quark mass ratios predicted z universal formula | EMPIRICAL TEST target | T5-T9 |
| Required g_0_q values (back-solved from PDG) | DERIVED-TARGET | T10 |
| Max achievable ratio with g_0 ∈ audit range | DERIVED-TARGET (structural ceiling) | T11 |
| Decision trichotomous (A−/B+/HALT-B) | EMPIRICAL VERDICT | Phase FINAL |

## §4 — Tautology test (CRITICAL per CALIBRATION_PROTOCOL §2.4)

**Question:** Are the output mass ratios definitionally reducible to PDG identity?

**Sympy substitution:** 
- m_i = c_M · A_tail(g_0_i)² · g_0_i^(e²/2) — A_tail and g_0 are INDEPENDENT parameters (g_0 = primary; A_tail = ODE-derived but treated jako empirical calibration here)
- c_M cancels w ratios; e² is mathematical constant; structure has TWO free per-quark parameters (g_0_i, A_tail_i)
- **6 quarks × 2 parameters = 12 free numbers** to fit 6 masses + cross-check audit range constraint
- But A_tail(g_0) is FUNCTIONAL (ODE-derived; calibrated z 3 lepton points → fixed function of g_0)
- So effectively **6 g_0_q parameters to fit 6 masses + 1 constraint** (audit range [0.817, 0.891])
- IF audit range is honored, the 6 g_0_q are restricted to span 1.09× ratio — structural test

**Verdict:** **NOT tautological** — formula's predictive content survives substitution. Audit range constraint is genuine BOUND (not free parameter). PASS tautology test.

## §5 — Falsifiability test (CRITICAL per CALIBRATION_PROTOCOL §2.5)

**Question:** Does an experimental value of axiom exist that would falsify the match?

**Pre-registered falsification rule (§0.2 README):**
> Jeśli formula z g_0_q ∈ [0.817, 0.891] NIE reprodukuje ≥3 z 5 ratios w 10% → INSUFFICIENT.

**Falsifier specificity:**
- 5 niezależnych ratios + 10% tolerance threshold
- 6 g_0_q parameters restricted do narrow range [0.817, 0.891] (1.09× span)
- A_tail(g_0) functional form fixed z lepton calibration

**Question:** Could a NUMEROLOGICAL fitting "save" the formula?

**Analysis:** With g_0_q restricted to 1.09× range:
- (g_0_max/g_0_min)^(e²/2) ≈ 1.091^3.6945 ≈ 1.378 (max enhancement)
- A_tail(g_0=0.817) ≈ 0.087 (extrapolated lepton power-law) → A_tail(g_0=0.891) ≈ 0.121
- (A_max/A_min)² ≈ (0.121/0.087)² ≈ 1.93
- **Max achievable mass ratio (within audit range): 1.378 × 1.93 ≈ 2.66**
- vs PDG quark ratios: 20 (smallest) to 80,000 (m_t/m_u)

**Verdict:** **FALSIFIABLE structurally** — max achievable ratio (2.66) << required ratios (20+). PASS falsifiability test (formula CAN fail by huge margin, not numerologically rescuable within audit constraint).

## §6 — Independent-path cross-validation

**Question:** Is there independent path from axioms to output predicting same result?

**Path A (this cycle):** Universal formula z lepton-calibrated A_tail(g_0), audit g_0 range constraint → predict 5 ratios.

**Path B (independent check):** Structural ceiling test (T11) — analytical max-min ratio within audit range. If max achievable ratio < 20 (smallest PDG ratio), formula DEFINITIVELY fails regardless of specific calibration.

**Path C (sanity):** Lepton subsector verification (T3) — formula reproduces PDG m_μ/m_e and m_τ/m_e (inheritance from why_n3 Phase 5). Independent confirmation formula structure correct.

**Verdict:** **3-path validation** structure (Path A predictions; Path B ceiling; Path C lepton sanity). Convergent: jeśli max ratio < 20 → ALL 5 quark predictions fail → HALT-B inevitable structurally.

## §7 — Inherited LOCKs (z poprzednich cykli)

| Inheritance | Status | Źródło |
|---|---|---|
| Universal mass formula m = c·A²·g_0^(e²/2) | LIVE (CLOSED 2026-05-01) | why_n3 Phase 5 |
| m_obs ≠ M_full operational distinction | LIVE | L05 A− 2026-05-16 |
| k_obs(α=2, d=3) = 3 | LIVE | L05 A− 2026-05-16 |
| β(α=2) = e²/2 canonical | LIVE | L08-e² B+ 2026-05-16 |
| Antisym Fock (γ_exchange=π) | LIVE | L08-FR A− 2026-05-16 |
| Cl(1,3) algebra (γ^μ) | LIVE | L08-Clifford A− 2026-05-16 |
| S_F^TGP standard form | LIVE | L08-Dirac A− 2026-05-16 |
| Wilson framework dla precision corrections | LIVE | L08-Dirac-Wilson A− 2026-05-16 |
| g_0_q ∈ [0.817, 0.891] universal ODE claim | LITERATURE LOCK | sek08b:529 |
| Lepton calibration (e, μ, τ) g_0/A_tail values | LIVE (CLOSED) | why_n3 PHASE1 + PHASE2 |
| PDG 2024 quark masses | LITERATURE | PDG external |
| S05 single-Φ axiom | LIVE | TGP_FOUNDATIONS §1 |

## §8 — 8/8 ☑ gate checklist (Phase 6 ABSOLUTE BINDING per CALIBRATION_PROTOCOL)

Per `meta/CALIBRATION_PROTOCOL.md` §3 audit gate:

```
☑ Phase 0 balance sheet exists (this file)
☑ Tautology test PASS (§4) — formula's predictive content survives substitution
☑ Falsifiability test PASS (§5) — falsifier triggers within structural ceiling
☑ Independent-path cross-validation (§6) — 3-path validation
☑ Alt-scan documented (single formula candidate; falsifier defines failure mode)
☑ NIE post-hoc structural motivations — audit range source is sek08b:529 (pre-existing)
☑ NIE circular anchor — A_tail calibrated z lepton ODE data independent z quark PDG
☑ NIE inheriting drift > parent × 5× — inherits LIVE A− cycles (L05, L08-Dirac etc.)
```

**Gate verdict:** 8/8 ☑ — proceed to Phase 1 sympy.

## §9 — 6/6 P-requirements gate (target)

| P# | Requirement | Verification source (Phase 1) |
|---|---|---|
| P1 | Universal formula structure m=c·A²·g_0^(e²/2) at α=2 | T1 (FP) |
| P2 | Mass ratio formula m_i/m_j = (A_i/A_j)²·(g_i/g_j)^(e²/2) | T2 (FP) |
| P3 | Lepton anchor sanity preserved (μ/e + τ/e match PDG <0.1%) | T3 (LIT inheritance) |
| P4 | 5 PDG quark ratios test outcome | T5-T9 (FP central tests) |
| P5 | Audit range cross-check + structural ceiling | T10-T12 (FP+LIT) |
| P6 | S05 single-Φ preserved | T13 (DEC, separate) |

## §10 — R-flags (risks)

Per README §risk_flags:

| R# | Risk | Mitigation |
|---|---|---|
| R1 | A_tail(g_0) ODE-derived vs power-law fit | Use 3-point power-law fit z lepton calibration; document fit residual; flag uncertainty |
| R2 | PDG quark masses scheme-dependent (~30% light) | Use bare PDG values; document scheme variation; tolerance 10% accommodates ~half |
| R3 | sek08b:529 normalization vs R3 g_0 | Treat audit range as primary constraint; flag if discordant z R3 lepton numerology |
| R4 | CKM mixing not modeled | Out of scope; flag as separate cycle component |

## §11 — Scope (in/out)

**IN:**
- Empirical test universal mass formula on 6 PDG quark masses (5 niezależnych ratios)
- A_tail(g_0) power-law calibration z 3 lepton points (sympy fit)
- Audit range constraint [0.817, 0.891] enforced
- Decision trichotomous A−/B+/HALT-B per pre-registered falsifier

**OUT:**
- CKM mixing derivation (separate cycle scope)
- Neutrino sector (separate cycle scope, also problem #3 component)
- Boson sector (W, Z, gluons) (separate cycle scope)
- Direct ODE solving dla quark g_0 (uses analytical extrapolation jako proxy; full ODE multi-session)
- Multi-loop QED Wilson corrections (L08-Dirac-Wilson framework documented, not deepened)

## §12 — Native vs Standard physics split

Per CYCLE_KICKOFF_TEMPLATE §0.3:

- **Native part:** Universal formula struct m = c·A²·g_0^(e²/2) (TGP why_n3 specific); g_0_q ∈ [0.817, 0.891] audit range (TGP-specific universal ODE claim)
- **Standard part:** Sympy algebra + power-law fit + PDG comparison (universal)
- **Mode:** **Native-with-mapping** — universal formula z TGP why_n3 + PDG comparison; L2 reduction to SM Yukawa nicht attempted (out of scope)

## §13 — Validator gate compatibility

Validator `tooling/validate_kickoff.py` checks:
- `contract::L1_native::output_observable` non-empty ✓ (5 quark mass ratios)
- `contract::L1_native::measurement_instrument` non-empty ✓ (PDG 2024 quark masses)
- `contract::L1_native::falsification_rule` pre-registered ✓ (≥3/5 ratios @ 10%)
- `pre_registration_date` matches cycle date ✓ (2026-05-16)
- 6/6 P-requirements declared ✓

**Status:** Phase 0 balance OK; **8/8 ☑** verified; proceed to Phase 1 sympy.

## Cross-references

- [[./README.md]] — kickoff contract (full)
- [[./Phase1_sympy.py]] — empirical test script (next deliverable)
- [[../../meta/CYCLE_KICKOFF_TEMPLATE.md]] §2.6 (Phase 0 template)
- [[../../meta/CALIBRATION_PROTOCOL.md]] §2 (8/8 gate)
- [[../../audyt/L08_kink_fermion_closure/README.md]] (problem #3 statement)
- [[../why_n3/PHASE2_n_alpha_derivation.md]] (universal formula + lepton data)
- [[../op-L05-mass-exponent-k-alpha-d-2026-05-16/Phase_FINAL_close.md]] (LIVE inheritance)
- [[../op-L08-Phase6-e2-derivation-2026-05-16/Phase_FINAL_close.md]] (LIVE inheritance)
- [[../../core/sek08b_ghost_resolution/sek08b_ghost_resolution.tex]] lin. 529 (audit range source)
