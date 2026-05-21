---
title: "Phase 0 — Balance: μ_ν^TGP astrofizyczna dyskryminacja scenarios A/B"
date: 2026-05-17
parent: "[[./README.md]]"
type: phase-zero-balance
phase: 0
status: 🟢 ACTIVE
sympy_substance_ratio: "6 FP / 1 LIT / 1 DEC = 75% FP"
hardcoded_T_pass: 0
---

# Phase 0 — Balance discrimination cyklu

## §1 — Inputs

### §1.1 — TGP predictions (LIVE z cycles 3, 4, 6)

| Element | Source | Value | Status |
|---|---|---|---|
| μ_ν^TGP_A central (scenario A) | Cycle 3 spinor B | 3.55·10⁻¹² μ_B | LIVE |
| μ_ν^TGP_A CI (m_X propagation) | Cycle 4 T6 | [1.28, 3.55]·10⁻¹² μ_B | LIVE |
| μ_ν^TGP_A geomean (log-fair central) | Cycle 4 T6 | 2.13·10⁻¹² μ_B | LIVE |
| μ_ν^TGP_A log-σ (m_X anchor) | Cycle 4 T6 | 0.22 dex | LIVE |
| μ_ν^TGP_B central (scenario B) | Cycle 6 T5 Lee-Shrock | 3.2·10⁻²⁰ μ_B | LIVE |
| μ_ν^TGP_B log-σ (v_H + m_e + m_ν) | Cycle 6 propagation | ~0.3 dex (m_ν dominant) | LIVE |
| m_X anchor (L06) | L06 numerical anchor | 60 MeV (factor 1.7 target 100 MeV) | NUMERICAL ANCHOR |
| Suppression power n (heuristic) | Cycle 3 placeholder | n = 2 | OPEN heuristic |

### §1.2 — Astrophysical bounds (LIT z literature)

| Bound | μ_max (μ_B) | CL | Stellar/bound log-σ | Source |
|---|---|---|---|---|
| TRGB Capozzi-Raffelt 2020 | 1.2·10⁻¹² | 2σ (95%) | 0.30 dex | arXiv:2007.03694 |
| SN1987A Magill+2018 (updated) | ~1.3·10⁻¹² | 95% CL | 0.45 dex | Phys. Rev. D 98, 115015 |
| ωCen Arceo-Diaz+2015 | 2.2·10⁻¹² | 95% CL | 0.30 dex | Astropart. Phys. 70, 1 |
| M5 Viaux+2013 | 4.5·10⁻¹² | 95% CL | 0.30 dex | A&A 558, A12 |
| Raffelt 1990 (classical SN/RG) | 3·10⁻¹² | 95% CL | 0.50 dex | Phys Rep 198, 1 (conservative anchor) |
| BBN N_eff Cyburt+2016 | ~10⁻¹⁰ | cosmological | 0.20 dex | Rev. Mod. Phys. 88, 015004 |
| Solar RSFP Borexino+SuperK | ~7·10⁻¹¹ | 90% CL | 0.30 dex (B_⊙) | Phys. Rev. D 96, 091103 |
| BH accretion Latimer-Burrows 2007 | ~10⁻¹⁰ | model-dep | 0.50 dex | ApJ 661, 320 |

### §1.3 — Methodology inputs (cycle 4 REPLICATE)

| Element | Source | Status |
|---|---|---|
| Log-space combined σ formula | Cycle 4 T6 | LIVE methodology |
| log_diff = log10(TGP_geomean) - log10(bound) | Cycle 4 T6 | LIVE |
| combined_σ = sqrt(log_σ_TGP² + log_σ_bound²) | Cycle 4 T6 | LIVE |
| Joint tension = log_diff / combined_σ | Cycle 4 T6 | LIVE |

## §2 — Outputs

### §2.1 — Phase 1 deliverables

- **T1 LIT:** Comprehensive bound survey z sources + per-bound systematic uncertainty
- **T2-T6 FP:** Per-bound joint σ_tension dla scenario A i scenario B (5 bounds × 2 scenarios = 10 σ values)
- **T7 FP:** Joint statistical discrimination verdict per pre-registered decision tree
- **T8 DEC:** S05 preservation declaration

### §2.2 — Phase FINAL deliverable

- Discrimination verdict: A- DISCRIMINATION / A- BOTH CONSISTENT / B+ PARTIAL / HALT-B
- Per-bound table z joint σ values
- PR-016 status post-discrimination
- Downstream implications (problem #3 neutrino, XLZD/DARWIN forecasting)

## §3 — Risk register

| Risk | Severity | Mitigation |
|---|---|---|
| **R1** Per-bound systematic uncertainties (TRGB stellar opacity ~0.3 dex; SN1987A emissivity ~0.5 dex; BBN N_eff ~0.2 dex; Solar B_⊙ ~factor 2; BH disk ~0.5 dex) | medium | Use literature-anchored systematic log-σ per bound (Phase 1 T1) |
| **R2** Correlated systematics (plasmon physics shared TRGB+SN+globular) — naive combination overestimates discrimination | high | Treat bounds INDIVIDUALLY (no naive multiplication); honest disclosure correlated systematics |
| **R3** Scenario A na granicy tightest bounds — methodology choices could tip verdict w obu kierunkach | high | Use cycle 4 protocol EXACTLY (joint CI log-space combined-σ); pre-registered thresholds |
| **R4** Scenario B trivially compatible (10⁻²⁰ << 10⁻¹²) — not real discrimination | medium | Document honestly: scenario B passes all bounds trivially; primary discrimination question is whether scenario A fails any |
| **R5** BH accretion bound model-dependent — could be cherry-picked | low | Use Latimer-Burrows 2007 conservative ~10⁻¹⁰ μ_B; flag model-dependence honestly |
| **R6** Post-hoc threshold adjustment temptation | high | Pre-registered §0.2 thresholds strictly enforced; no re-tuning |
| **R7** Cycle 4 NO TENSION verdict mogłaby fail to extend gdy multiple bounds aggregate | medium | Examine each bound z TRZECH possible thresholds (>2σ, 1-2σ, <1σ); decision tree handles aggregate honestly |

## §4 — Methodology binding

**This cycle REPLICATES cycle 4 joint CI methodology EXACTLY:**

```python
# Per cycle 4 Phase1_sympy.py T6:
log_TGP_mid = log10(mu_TGP_geomean)         # scenario-specific
log_bound = log10(mu_max_2sigma)            # per-bound
log_diff = log_TGP_mid - log_bound
combined_log_sigma = sqrt(TGP_log_σ² + bound_log_σ²)
sigma_tension = log_diff / combined_log_sigma
```

**Decision per σ_tension:**
- σ > 2: TENSION REAL (discrimination)
- 1 < σ ≤ 2: TENSION MARGINAL (flag)
- σ ≤ 1: NO TENSION (consistent)

**Sign convention:** σ_tension > 0 means TGP above bound (constrained); σ_tension < 0 means TGP below bound (compatible, no constraint).

## §5 — 8/8 gate checklist

- [x] **G1:** Phase 0 balance sheet exists (this file)
- [x] **G2:** README z BINDING contract + pre-registered falsifier (§0.2)
- [x] **G3:** Sympy substance plan: 6 FP + 1 LIT + 1 DEC = 75% FP
- [x] **G4:** Zero hardcoded T_pass=True target (Phase 6 ABSOLUTE BINDING)
- [x] **G5:** Pre-flight read confirmation §0.4 (5 methodology files + 4 predecessor)
- [x] **G6:** TGP-native check Q1-Q8 (§0.3)
- [x] **G7:** P-requirements target 6/6 (P1-P6 declared w README)
- [x] **G8:** Risk register §3 z 7 risks honestly identified

**Verdict:** 🟢 **8/8 PASS.** Cycle ready dla Phase 1 sympy execution.

## §6 — Scope discipline

**IN-SCOPE:**
- Comprehensive bound survey z established literature
- Per-bound joint CI tension computation (cycle 4 methodology)
- Pre-registered decision tree application
- Honest verdict report

**OUT-OF-SCOPE (deferred to future cycles):**
- New TGP-native derivation of μ_ν (problem #3 boson sub-component, multi-session)
- Refinement of suppression power n (rigorous QED loop, post-W/Z sektor closure)
- L_X structural derivation (HALT-B w cycle 5)
- XLZD/DARWIN numerical sensitivity refinement (separate observational cycle)
- Bayesian combined likelihood z correlated systematics (requires dedicated stats cycle)

---

**Sign-off:** Claudian @ 2026-05-17 (7th cycle sesji 2026-05-17, post-cycle-6 dual-scenario discrimination).
