---
title: "FINDINGS — op-L01-N1-EM-trace-anomaly-TGP-2026-05-11 (eksportowalne wyniki)"
date: 2026-05-11
parent: "[[./README.md]]"
type: findings
total_findings: 25
total_sympy_pass: "16/16 (Phase 1 + Phase 2)"
six_requirements: "6/6 RESOLVED"
classification: STRUCTURAL_DERIVED
tags:
  - findings
  - eksportowalne
  - L01-N1-closure
  - cross-cycle-N1-derived
---

# FINDINGS — eksportowalne wyniki cyklu

## §1 — Phase 1: 1-loop QED on g_eff[{Φ_i}]

| ID | Finding | Source |
|---|---|---|
| **F1.1** | β-function 1-loop QED `β(α) = (2/(3π))·α²` (single Dirac fermion, Capper-Duff-Halpern 1974) | Phase 1 §1.4 + sympy T1 |
| **F1.2** | Trace anomaly prefactor `β(α)/(2α) = α/(3π) ≈ 7.74·10⁻⁴` (sympy LOCK) | Phase 1 §1.4 + sympy T2 |
| **F1.3** | **L01 ADDENDUM §3.2 Q3 typo identified:** `α²/(3π) ≈ 7.7·10⁻⁷` → correct value `α/(3π) ≈ 7.74·10⁻⁴` (factor ~1000 typo) | Phase 1 §1.4 |
| **F1.4** | Trace anomaly explicit form: `T^μ_μ_EM,1-loop = (α/(3π))·F² + curvature × F² mixing + Riegert localized terms` | Phase 1 §1.3 |
| **F1.5** | Riegert localization σ_eff = function(ψ) — S05 single-Φ preserved bezwarunkowo | Phase 1 §2.2 + sympy T7 |
| **F1.6** | Generic 3-funkcyjny ansatz {A(ψ), B(ψ), C(ψ)} per emergent-metric Phase 1 — R1 (M9.1'' contamination) eliminowany strukturalnie | Phase 1 §1.1 + sympy T6 |
| **F1.7** | Photon dispersion 1-loop preserved canonical (ω² = c²k²) po field redefinition | Phase 1 sympy T8 (R3 partial) |
| **F1.8** | Phase 2 disjointness target: curvature × F² operators TGP-reduced **NIE produkują** ψ.1.v3 dim-6 EFT form (anticipated, verified Phase 2) | Phase 1 §3 |

## §2 — Phase 2: TGP-specific reduction + Theorem 2.1 (Disjointness)

| ID | Finding | Source |
|---|---|---|
| **F2.1** | R[g_eff] reduction w 1PN limit (b_1 = -a_1, γ=1): `R = a_1·∇²h`; pure-2-derivative ψ scalar | Phase 2 §1.2 + sympy T1 |
| **F2.2** | R^{μν}[g_eff] reduction: `R^{ij} = -a_1·∂^i∂^j h` (tensor 2-derivative; scalar part identically zero with γ=1) | Phase 2 §1.3 + sympy T2 |
| **F2.3** | T^μ_μ_EM,1-loop_TGP zawiera 6 distinct operator classes: `(α/(3π))·F²`, `(∂²ψ)·F²`, `(∂_μ∂_ν ψ)·F^{μρ}F^ν_ρ`, `σ_ab·F²`, `□F²`, `Riegert local with σ_eff = function(ψ)` | Phase 2 §1.4 + sympy T3 |
| **F2.4** | **Theorem 2.1 (Disjointness):** trace anomaly TGP-reduced operator classes są strukturalnie DISJOINT od ψ.1.v3 canonical basis B = {L₅'_a, L₅'_b}; różnice ∂ψ leg counting + tensor structure | Phase 2 §2 + sympy T4 |
| **F2.5** | **Q1 closure dedicated derivation:** ψ.1.v3 *NIE pokrywa* quantum trace anomaly EM (operator class disjoint), N1 dedicated cycle strukturalnie *konstruktywnie* potwierdzony | F2.4 + Phase 2 §2.2 |
| **F2.6** | GW170817 c_GW=c_EM full check: Δc/c ≈ (α/(3π))·R_cosmo/m_e² ≈ **10⁻⁸⁰** ≪ 9·10⁻²² bound (~58 OOM margin) | Phase 2 §3.1 + sympy T5 |
| **F2.7** | c_GW dispersion **NIE modyfikowane** przez 1-loop QED (no graviton in QED loop); c_GW = c structurally preserved | Phase 2 §3.2 + sympy T6 |
| **F2.8** | R5 (perturbative breakdown) documented: cycle valid B ≪ B_QED ≈ 4.4·10⁹ T; B ≳ B_QED extreme magnetar surface defer non-perturbative analysis | Phase 2 sympy T7 |
| **F2.9** | S05 single-Φ preserved verified: all reduced operators contain only Φ source field | Phase 2 sympy T8 |
| **F2.10** | ρ_EM_quantum[{Φ_i}] = native source field dla Φ-EOM, computable z Phase 1+2 framework w obecności emergent-metric ansatz {A, B, C} | Phase 2 §1.4 |

## §3 — Phase 3: Phenomenology

| ID | Finding | Source |
|---|---|---|
| **F3.1** | Lab regime ρ_EM_quantum (B=1 T): \|ρ\| ~ 7·10⁻¹⁵ kg/m³ inside electromagnet volume; niewykrywalne acceleration | Phase 3 §1.1 |
| **F3.2** | Magnetar typical (B=10¹⁰ T): ρ_EM_quantum/ρ_NS ~ 1.7·10⁻¹² | Phase 3 §2.1 |
| **F3.3** | Magnetar extreme (B=10¹¹ T, B/B_QED ~23): ratio ~1.7·10⁻¹⁰ z **corrected** α/(3π) prefactor (NIE 10⁻¹² jak L01 ADDENDUM cytuje z typo) | Phase 3 §2.2 |
| **F3.4** | L01 ADDENDUM §3.2 Q3 typo correction propagated: `α²/(3π) ≈ 7.7·10⁻⁷` → `α/(3π) ≈ 7.74·10⁻⁴`; 2-3 OOM correction w magnetar estimates | Phase 3 §2.3 |
| **F3.5** | η_TGP_EM_quantum (Pt vs Ti) = **0 strukturalnie** z universal coupling structure (S05) | Phase 3 §3.1 |
| **F3.6** | MICROSCOPE bound automatic PASS: η_TGP = 1.32·10⁻²⁶ (B9 baseline preserved) ≪ 10⁻¹⁵, ~11 OOM margin | Phase 3 §3.2 |
| **F3.7** | Eöt-Wash + LLR Nordtvedt: automatic PASS z universal coupling | Phase 3 §3.3-§3.4 |
| **F3.8** | **R6 (QEP violations from Asorey-2015 type non-local Lagrangians) closed strukturalnie** — TGP universal coupling z S05 immune do tej klasy QEP violation mechanism | Phase 3 §4 |
| **F3.9** | TT10 magnetar test mechanism decoupling preserved: trace anomaly ≪ surface ρ_NS dla all magnetar regimes; **τ.3 L4 mechanism testowany czysto** | Phase 3 §2.4 |

## §4 — Phase 4: Three-layer L1/L2/L3 closure

| ID | Finding | Source |
|---|---|---|
| **F4.1** | L1 native predictions: ρ_EM_quantum jest **natywnym source field** dla Φ-EOM (NIE PPN/ppE artifact); 5 native obserwabli enumerated (lab, magnetar, solar system, Schwinger lab, cosmology) | Phase 4 §1 |
| **F4.2** | L2 projection: quantum trace anomaly **NIE wprowadza nowych** PPN/ppE/cosmology free parameters; modyfikuje precision (deferred precision items, NIE new free coefs) | Phase 4 §2 |
| **F4.3** | L3 falsification map: 7 observational bounds enumerated z status; 4 future tests dla concrete falsifier roles (MICROSCOPE-2, CE GW170817-analog, Schwinger lab, magnetar atomic spectroscopy) | Phase 4 §3 |
| **F4.4** | Native parameter audit: 6 derived items + 6 forced strukturalnie + 0 disjoint sektory (Theorem 2.1) + 3 deferred precision; **NIE rośnie liczba swobodnych parametrów** od N1 closure | Phase 4 §4 |
| **F4.5** | **6/6 P-requirements RESOLVED** (P1-P6 all PASS); STRUCTURAL_DERIVED classification | Phase 4 §5 |

## §5 — Cross-cycle (synthesis findings)

| ID | Finding | Source |
|---|---|---|
| **F5.1** | **L01 N1 closure source:** ten cykl daje *konstruktywną* derywację N1 (zamiast deferred placeholder); zamyka N1 w klasie STRUCTURAL_DERIVED | Phase 4 §6-§7 |
| **F5.2** | **Cross-cycle convergence updated diagnostic** — pięć niezależnych diagnoz (L01, τ.3, ψ.1, Q2, op-L01-N1) zbieżnych na separable sector structure jako *strukturalna własność TGP*, NIE post-hoc tuning | Phase 2 §6 |
| **F5.3** | **L01 ADDENDUM §3.2 Q3 typo correction propagation list** (Phase 4 §8): L01 ADDENDUM, L01 NEEDS, L01 README, τ.3 ADDENDUM, PREDICTIONS_REGISTRY | Phase 3 §7 + Phase 4 §8 |

## §6 — Risk register (final status)

| Risk | Status | Source |
|---|---|---|
| **R1** (M9.1'' contamination) | **closed strukturalnie** — generic ansatz {A, B, C} per emergent-metric, NIE M9.1'' specific | Phase 1 §4.1 + sympy T6 |
| **R2** (operator class re-overlap z ψ.1.v3) | **closed konstruktywnie** — Theorem 2.1 explicit | Phase 2 §2 + sympy T4 |
| **R3** (GW170817 c-violation) | **closed numerycznie** — Δc/c ~ 10⁻⁸⁰ ≪ bound 10⁻²², ~58 OOM margin | Phase 2 §3 + sympy T5+T6 |
| **R4** (S05 violation z quantum loops) | **closed strukturalnie** — Riegert σ_eff = function(ψ) | Phase 1 §2.2 + sympy T7 + T8 |
| **R5** (perturbative breakdown w extreme magnetar) | **honestly documented** — cycle valid B ≪ B_QED; deferred non-perturbative cycle | Phase 2 §3.4 + Phase 3 §5 + sympy T7 |
| **R6** (QEP violations from Asorey-2015) | **closed strukturalnie** — TGP universal coupling z S05 immune | Phase 3 §4 |

## §7 — Sympy LOCK summary

| Phase | Tests | PASS | Notes |
|---|---|---|---|
| Phase 1 | 8 | 8/8 | β-function, prefactor, conformal invariance, anomaly form, dimensional, R1+R4 guards, dispersion |
| Phase 2 | 8 | 8/8 | R reduction, R^μν reduction, operator enumeration, disjointness, GW170817, c_GW, B_QED limit, S05 |
| **Total** | **16** | **16/16** | **100% PASS** |

## §8 — Predictions for PREDICTIONS_REGISTRY (planned entries Phase_FINAL)

| Entry ID | Prediction | Test | Status |
|---|---|---|---|
| **M911-EM-quantum-trace-anomaly** | Quantum trace anomaly EM contributes ρ_EM_quantum z prefactor `α/(3π) ≈ 7.74·10⁻⁴` × F²/c² | Lab (B=1 T): undetectable; Schwinger-class lab (2030+) future | TESTED-PASS (lab undetectable as predicted) |
| **M911-EM-quantum-magnetar** | Magnetar B=10¹¹ T: ρ_EM_quantum/ρ_NS ~ 10⁻¹⁰ (corrected) | TT10 magnetar X-ray timing — decoupled (testuje L4); deferred non-perturbative for B ≳ B_QED | LIVE — deferred precision |
| **M911-EM-quantum-MICROSCOPE** | η_TGP_EM_quantum (Pt vs Ti) = 0 strukturalnie | MICROSCOPE 2017+: η ≤ 10⁻¹⁵; future MICROSCOPE-2 η ≤ 10⁻¹⁸ | TESTED-PASS |
| **M911-EM-quantum-GW170817** | Δc/c od trace anomaly ~ 10⁻⁸⁰ | GW170817 |c_GW/c| ≤ 9·10⁻²² | TESTED-PASS (~58 OOM margin) |

---

**Findings exported.** 25 distinct findings + 6/6 risks + 16/16 sympy + 6/6 P-requirements.
