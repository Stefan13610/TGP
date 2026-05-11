---
title: "Phase 4 — Three-layer L1/L2/L3 closure + native parameter audit + 6/6 P-requirements verify"
date: 2026-05-11
parent: "[[./README.md]]"
type: phase-results
phase: 4
status: 🟢 RESOLVED — three-layer presentation complete + 6/6 P-requirements PASS
sub_needs_resolved: [N0.12]
six_requirements_status: "6/6 RESOLVED (P1-P6)"
sympy_total: "24/24 PASS (Phase 1 + Phase 2 + Phase 3)"
predecessor: "[[./Phase3_results.md]]"
tags:
  - phase4
  - three-layer-presentation
  - native-observables-first
  - native-parameter-audit
  - six-requirements-verify
  - cosmology-mandatory
---

# Phase 4 — Three-layer L1/L2/L3 closure

## §0 — Methodology binding

Per [[../../meta/PPN_AS_PROJECTION.md]] §3.1 mandatory three-layer presentation +
§4 cosmology generalization (binding 2026-05-10+).

## §1 — L1: Native predictions (primary)

### §1.1 — Native source field: ρ_QCD(T)

**Native obiekt** (NIE PPN/ppE artifact):

```
ρ_QCD(T)[{Φ_i}] = ρ_QCD_vacuum + ρ_QCD_thermal(T)

  ρ_QCD_vacuum  = -(b_0/8) · ⟨α_s G²/π⟩ / c_0² ≈ 2.8·10¹⁸ kg/m³
                  (constant SVZ gluon condensate; substrate-decoupled per Q2 F1)
  ρ_QCD_thermal(T) = -Δ(T) · T⁴ / c_0²
                  Δ(T) = (ε(T) - 3p(T))/T⁴ z lattice EoS profile (HotQCD)
                  Peak Δ_max ≈ 4 near T_c = 156 ± 9 MeV (crossover)
                  Δ(T → 0) → 0 strukturalnie (Q2 F1 verification)
                  Δ(T → ∞) → 0 (Stefan-Boltzmann conformal)
```

z prefactor `b_0/8 = 9/8` dla N_f=3 light quarks (low-T) lub `b_0/8 = 7/8` dla
N_f=6 (high-T, sympy LOCK Phase 1 T1).

### §1.2 — Native L1 observables

| Observable | Native form | Constraint na native coefs |
|---|---|---|
| **(L1.a)** Vacuum gluon condensate `⟨α_s G²/π⟩_0` | SVZ-1979 + lattice 2018+ | external input z uncertainty band [0.005, 0.020] GeV⁴ |
| **(L1.b)** Crossover temperature T_c | HotQCD lattice consensus | 156 ± 9 MeV (chiral susceptibility peak) |
| **(L1.c)** Hubble H(T) w QCD epoce | Friedmann z transient ρ_QCD(T) | t_H(T_c) ~ 10⁻⁵ s (standard cosmology consistent) |
| **(L1.d)** ⁴He primordial Y_p (BBN) | Y_p = 0.247 ± 0.001 | TGP = ΛCDM (post-QCD non-source) |
| **(L1.e)** D/H (BBN) | D/H = 2.5·10⁻⁵ ± 0.1·10⁻⁵ | TGP = ΛCDM |
| **(L1.f)** ω_Λ (today) | ω_Λ = 0.6889 (Q2 F1 + T-Λ) | g̃ ≈ 0.98 z Q2 derivation |
| **(L1.g)** Stochastic GW background QCD epoch | weak (smooth crossover) | NIE dominantny dla PTA |

### §1.3 — Strukturalna konsekwencja

**Native level:** ρ_QCD(T) jest *transient phase-transition source* dla Φ-EOM w
hot QCD epoch (z ~ 10¹², t ~ 10⁻⁵ s). W IR limit (T<<Λ_QCD) reduces to constant
substrate-decoupled vacuum condensate.

**Cross-cycle Q2 F1 KONSTRUKTYWNIE VERIFIED dla QCD sektora** (Phase 2 §3) +
**Q2 closure + tej cyklu BBN/CMB/PTA bound checks all PASS automatic.**

## §2 — L2: Projection charts

### §2.1 — Cosmology projection (FRW, perturbations)

| Cosmology observable | Native counterpart | TGP prediction (post-N2) |
|---|---|---|
| **w_eff(z)** | g_eff[Φ̄(t)] z V(Φ_eq)=const | -1 EXACT (matter-decoupling preserved) |
| **H(z) w QCD epoce** | Friedmann z ρ_QCD_thermal source | standard ΛCDM + lattice EoS |
| **H(z~10⁹) BBN** | post-QCD non-source ρ_QCD ≈ 0 | standard BBN H(z) |
| **CMB ω_b, ω_m** | matter-decoupling per Q2 F1 | preserved exactly |
| **Σ(z, k), μ(z, k)** | modified gravity charts | unchanged od emergent-metric prediction |

### §2.2 — PPN/ppE projection (gravity sector)

QCD trace anomaly **NIE wprowadza nowych** PPN/ppE parametrów — analogous do N1
cycle architecture.

| PPN/ppE parameter | Trace anomaly contribution |
|---|---|
| γ_PPN | unchanged (matter-coupling structure preserved) |
| β_PPN | unchanged |
| β_ppE (2.5PN GW) | unchanged (per emergent-metric Phase 4) |
| α_ppE (amplitude) | unchanged |

### §2.3 — Stochastic GW spectrum chart

| GW band | TGP prediction | Observation |
|---|---|---|
| LIGO (kHz) | h_TT^σ = h_TT^GR exact | BBH detections preserved |
| LISA (mHz) | crossover-suppressed QCD signal | not yet operational; TGP NIE predicts dominant QCD |
| **PTA (nHz)** | **No first-order QCD signal** | NANOGrav 15-yr SMBHB consensus preserved |

## §3 — L3: Falsification map

### §3.1 — Per-observation falsifier mapping

| Observational bound | Constrains native item | Status |
|---|---|---|
| **PDG 2024 ⁴He Y_p = 0.245 ± 0.003** | H(z~10⁹) BBN preservation | PASS within 0.55σ |
| **PDG 2024 D/H = 2.527·10⁻⁵ ± 0.030·10⁻⁵** | analog | PASS within 0.26σ |
| **Planck 2018 ω_b = 0.02237 ± 0.00015** | matter-decoupling Q2 F1 | PASS exactly |
| **Planck 2018 ω_m = 0.1430 ± 0.0011** | analog | PASS exactly |
| **Planck 2018 Ω_Λ = 0.6889 ± 0.0056** | g̃ via T-Λ + Q2 | PASS (g̃=0.98 derived) |
| **T-Λ ratio = 1.020 (closure_2026-04-26)** | matter sector decoupled bezwarunkowo | PASS empirically |
| **NANOGrav 15-yr stochastic GW** | crossover smoothness (lattice 2+1 flavor) | compatible (SMBHB) |
| **HotQCD T_c = 156 ± 9 MeV** | lattice consensus | input (NOT prediction) |
| **SVZ ⟨α_s G²/π⟩_0 = 0.012 GeV⁴** | sum rules + lattice | input z uncertainty band |
| **PDG 2024 Λ_QCD = 217 ± 8 MeV** | dimensional transmutation OOM-correct | input |

### §3.2 — Concrete falsifier scenarios

**Co BY zafalsyfikowało N2 closure verdict:**

1. **BBN observation Y_p > 0.250 lub Y_p < 0.240** (5σ deviation): falsyfikuje
   matter-decoupling assumption → falsyfikuje Q2 F1 + TGP framework.
2. **Discovery first-order QCD phase transition w lattice** (sprzeczne z
   2+1 flavor consensus): falsyfikuje crossover assumption → predicted
   dominant QCD PTA signal byłby observable.
3. **CMB ω_b, ω_m measurements szyfują** (post-Planck 2018 5σ shift): mogłyby
   sygnalizować problem z matter-decoupling assumption.
4. **Detection ratio T-Λ deviating dramatically od 1.020** (e.g., > 50% off):
   falsyfikuje Q2 F1 mechanism + cycle base assumption.

### §3.3 — Future tests (concrete falsifiers planned)

| Future test | Target precision | Falsifier role |
|---|---|---|
| CMB-S4 (2030+) | ω_b ~ 10⁻⁴ precision | Higher precision matter-decoupling test |
| LiteBIRD CMB | improvement on r, n_s | Inflation epoch (NOT direct N2) |
| LISA stochastic GW (2035+) | mHz band | EW phase transition (separate from QCD N2) |
| SKA PTA (2030+) | improved nHz band | SMBHB vs new physics discrimination |
| Lattice QCD precision | sub-percent T_c, Δ(T) | refinement (NOT new constraint) |

## §4 — Native parameter audit

```
Independent native source structures fixed by op-L01-N2-QCD cycle:

Constrained / DERIVED w tym cyklu:
  β_QCD(g) = -(b_0/(16π²))·g³ z b_0 = (11/3)N_c - (2/3)N_f         [sympy T1]
  Trace anomaly form T^μ_μ_QCD = (β(g)/(2g))·G²                    [sympy T4]
  Riegert decomposition w g_eff[{Φ_i}] z σ_eff = function(ψ)       [Phase 1 §5]
  Asymptotic freedom (β < 0)                                       [sympy T2]
  Q2 F1 konstruktywna verification dla QCD sektora                 [Phase 2 §3]
  BBN/CMB/PTA bounds preserved                                     [Phase 3]

Forced from substrate symmetry / Q2 mechanism (NIE wprowadzane w tym cyklu):
  ρ_QCD_vacuum substrate-decoupled od bare Λ (Q2 F1 mechanism)
  ρ_QCD_thermal(T → 0) → 0 (lattice IR limit + Q2 F1)
  S05 single-Φ axiom (Riegert σ_eff = function(ψ))
  c_GW = c_EM (no graviton in QCD loops; same as N1 §3.2)
  WEP universality (universal coupling structure; same as N1 R6 closure)
  
External lattice/sum-rule inputs (NOT new free params, deferred precision):
  Λ_QCD = 217 ± 8 MeV (PDG 2024 MS-bar N_f=5)
  ⟨α_s G²/π⟩_0 = 0.012 GeV⁴ ± O(50%) (SVZ + lattice)
  ⟨q̄q⟩ chiral condensate ≈ -(250 MeV)³ (lattice)
  T_c = 156 ± 9 MeV (HotQCD 2+1 flavor)
  α_s(M_Z) = 0.1179 ± 0.0009 (PDG 2024 world avg)
  Quark masses m_u, m_d, m_s, m_c, m_b, m_t (PDG 2024)
  HotQCD/Wuppertal-Budapest lattice EoS Δ(T) profile

Disjoint sectors (verified, structural separation):
  ψ.1.v3 dim-6 EFT basis (operator class) — analog do N1 Theorem 2.1
  N1 cycle EM trace anomaly (separate gauge sector)
  N4 Higgs sector (deferred)
  N5 EW gauge anomaly (deferred)

Total native parameter count for L01-N2 sektora:
  Constrained:   1 prefactor (β_QCD) + 1 trace anomaly + Riegert structure = 3 derived items
  Forced:         5 strukturalne tożsamości (Q2 F1 + S05 + GW170817 + WEP + IR limit)
  External:       7 lattice/sum-rule inputs (NIE new free params)
  Disjoint:      verified vs ψ.1, N1, future N4/N5
```

**Konkluzja audit:** N2 closure cycle:
- **NIE wprowadza nowych swobodnych parametrów** w TGP framework (analogous do N1).
- **Modyfikuje precyzję** native predictions w QCD epoce (cosmology integration);
  external lattice inputs traktowane jako *deferred precision*, NIE *new free*.
- **Konstruktywnie potwierdza** Q2 F1 mechanism dla QCD sektora.
- **Strukturalnie consistent** z all observational bounds (BBN/CMB/PTA + T-Λ).

## §5 — Six P-requirements verification

| # | Requirement | Status | Evidence |
|---|---|---|---|
| **P1** | ρ_QCD(T) explicit formula z β_QCD(g) + Λ_QCD lattice input w `g_eff[{Φ_i}]` (NIE M9.1''!) | ✅ **PASS** | Phase 1 §1, sympy T7 (R1 guard); generic ansatz {A, B, C} preserved |
| **P2** | Cosmology integration — Friedmann eq z transient ρ_QCD(T) source w QCD epoce | ✅ **PASS** | Phase 2 §2; lattice EoS profile + H(T_c) = 10⁻²² GeV consistent |
| **P3** | BBN ⁴He, D/H constraint preserved (~1% precision) | ✅ **PASS** | Phase 3 §1; Y_p 0.55σ, D/H 0.26σ deviations |
| **P4** | CMB ω_b/ω_m preservation (matter-decoupling cross-check z Q2) | ✅ **PASS** | Phase 3 §2; Planck 2018 preserved exactly + T-Λ ratio 1.020 |
| **P5** | PTA stochastic GW background compatibility (NANOGrav 15-yr 2023) | ✅ **PASS** | Phase 3 §3; crossover lattice consensus → no first-order signal |
| **P6** | Q2 F1 closure (matter vacua decoupled od bare Λ) verified konstruktywnie z dedicated QCD derivation | ✅ **PASS** | Phase 2 §3; Q2 F1 konstruktywnie verified; T-Λ ratio 1.020 preserved |

**6/6 RESOLVED** — STRUCTURAL_DERIVED classification verified.

## §6 — Cumulative cycle status

```
op-L01-N2-QCD-trace-anomaly-cosmology-2026-05-11:
  Phase 0 (balance sheet):                     6/6 gate PASS         ✅ DONE
  Phase 1 (formal QCD trace anomaly):          8/8 sympy PASS         ✅ DONE
  Phase 2 (cosmology + Q2 F1 verify):          8/8 sympy PASS         ✅ DONE
  Phase 3 (BBN/CMB/PTA bounds):                8/8 sympy PASS         ✅ DONE
  Phase 4 (three-layer closure):               6/6 P-requirements RESOLVED ← TUTAJ

Cumulative sympy: 24/24 PASS (Phase 1+2+3)
Six requirements: 6/6 RESOLVED
Risks: R1-R7 all addressed:
  R1 (M9.1'' contamination):                CLOSED Phase 1 (generic ansatz)
  R2 (non-perturbative regime):             HONESTLY DOCUMENTED Phase 1 (lattice external)
  R3 (BBN incompatibility):                 CLOSED Phase 3 §1 (suppression 10⁻⁶⁸)
  R4 (S05 violation):                       CLOSED Phase 1 (composite operator + Riegert)
  R5 (Q2 inconsistency):                    CLOSED Phase 2 §3 (konstruktywna verification)
  R6 (PTA false positive):                  CLOSED Phase 3 §3 (crossover consensus)
  R7 (crossover not phase transition):      CLOSED Phase 2 §4 (lattice 2+1 flavor)
```

**Status post-Phase-4:** **STRUCTURAL_DERIVED**.

## §7 — Cross-cycle propagation list

Phase_FINAL_close.md will execute:

1. **L01 NEEDS.md §N2 status update:** IN_PROGRESS → **CLOSED** z linkiem do
   Phase_FINAL_close.md.
2. **L01 README.md status update:** N2 closure note (analog do N1 closure).
3. **L01 NEEDS Q2 + T.2 entry update** z konstruktywna Q2 F1 verification.
4. **Q2 cycle FINDINGS update:** N2 cycle uses F1+F4+F7 z konstruktywna potwierdzeniem.
5. **PREDICTIONS_REGISTRY.md** — new entries dla M911-QCD-trace-anomaly variants.
6. **closure_2026-04-26 T-Λ** cross-link — N2 dostarcza konstruktywne potwierdzenie
   że matter sector vacua NIE zaburzają T-Λ ratio.

## §8 — Cross-references

- [[./README.md]]
- [[./Phase0_balance.md]] (6/6 gate PASS)
- [[./Phase1_results.md]] (β_QCD LOCK, 8/8)
- [[./Phase2_results.md]] (cosmology integration + Q2 verify, 8/8)
- [[./Phase3_results.md]] (BBN/CMB/PTA bounds, 8/8)
- [[../../meta/PPN_AS_PROJECTION.md]] §3.1, §3.3, §4 (binding methodology)
- [[../op-L01-rho-stress-energy-bridge-2026-05-04/]] (parent cycle)
- [[../op-L01-N1-EM-trace-anomaly-TGP-2026-05-11/]] (sister cycle architecture)
- [[../op-Q2-vacuum-budget-2026-05-10/]] (parent-mechanism)

---

**Phase 4 close:** Three-layer L1/L2/L3 presentation complete. Native parameter
audit complete. **6/6 P-requirements RESOLVED.**

Next: Phase_FINAL_close.md (sign-off + 7/7 risks resolved + STRUCTURAL_DERIVED
classification) + FINDINGS.md + cross-cycle propagation updates.
