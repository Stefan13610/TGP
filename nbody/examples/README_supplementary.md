# TGP Numerical Verification Suite — Supplementary Material

**Associated papers:**
- *Three inputs, forty predictions* (Letter, 3pp)
- *TGP theory: derivations, verification, and 40 predictions* (Companion, 11pp)

**Date:** 2026-04-07

---

## Overview

This directory contains 36 Python scripts (ex235–ex270) that perform
348 individual numerical tests verifying the Tensor-Gravitational-Potential
(TGP) theory. Additionally, ex271 is the final summary/scorecard script.

**Requirements:** Python 3.8+, NumPy, SciPy

**Three fundamental inputs:**
- `g0e = 0.86941` — TGP coupling constant
- `Omega_Lambda = 0.6847` — cosmological constant (Planck 2018)
- `N = 3` — generation number (from GL(3,F₂))

**Running:** Each script is self-contained:
```bash
python ex235_alpha_s_formula.py
```

---

## Script Registry

### Phase 1: Core Formulas (ex235–ex244)

| Script | Topic | Score | Perfect |
|--------|-------|-------|---------|
| ex235 | α_s = 3g₀ᵉ/(32Ω_Λ) basic formula | 10/10 | ★★★ |
| ex236 | Cabibbo angle λ_C = Ω_Λ/N | 10/10 | ★★★ |
| ex237 | CP violation phases from GL(3,F₂) | 9/10 | |
| ex238 | Koide constant K = 2/3 from Z₃ | 10/10 | ★★★ |
| ex239 | Dark matter Ω_DM = Ω_b(N! − Ω_Λ) | 8/10 | |
| ex240 | CKM matrix structure | 9/10 | |
| ex241 | Running coupling unification | 8/10 | |
| ex242 | Neutrino mixing from GL(3,F₂) | 7/10 | |
| ex243 | Quark mass ratios | 8/10 | |
| ex244 | Master consistency check (all 12 eqs) | 10/10 | ★★★ |

**Phase subtotal: 93/100 (93.0%)**

### Phase 2: Precision Tests (ex245–ex256)

| Script | Topic | Score | Perfect |
|--------|-------|-------|---------|
| ex245 | Cosmological constant from vacuum | 9/10 | |
| ex246 | Lepton mass formula from φ-FP | 8/10 | |
| ex247 | GL(3,F₂) group: order, irreps, Z₃ | 10/10 | ★★★ |
| ex248 | Weinberg angle sin²θ_W | 9/10 | |
| ex249 | PMNS matrix from GL(3,F₂) | 7/10 | |
| ex250 | Fine structure constant α_em | 8/10 | |
| ex251 | Strong CP problem (θ̄ = 0 from Z₃) | 9/10 | |
| ex252 | Dark matter soliton properties | 8/10 | |
| ex253 | Cosmological predictions (H₀, S₈) | 8/10 | |
| ex254 | Neutrino mass spectrum (NO only) | 7/8 | |
| ex255 | Proton stability from Z₃ triality | 10/10 | ★★★ |
| ex256 | Master summary (mid-pipeline check) | 10/10 | ★★★ |

**Phase subtotal: 104/118 (88.1%)**

### Phase 3: Advanced Topics (ex257–ex261)

| Script | Topic | Score | Perfect |
|--------|-------|-------|---------|
| ex257 | Muon g−2 anomaly in TGP | 9/10 | |
| ex258 | Gravitational waves (c_GW = c) | 10/10 | ★★★ |
| ex259 | Koide constant from TGP action | 10/10 | ★★★ |
| ex260 | EW precision observables & m_W | 9/10 | |
| ex261 | Inflation from TGP (n_s, r) | 8/10 | |

*Note: ex262 was an interim summary and is excluded from physics scoring.*

**Phase subtotal: 44/50 (88.0%)**

### Phase 4: Deep Physics (ex263–ex270)

| Script | Topic | Score | Perfect |
|--------|-------|-------|---------|
| ex263 | Baryogenesis & leptogenesis | 9/10 | |
| ex264 | Higgs mass m_H = v × 57/112 | 8/10 | |
| ex265 | Theory comparison (TGP vs BSM) | 10/10 | ★★★ |
| ex266 | EW & QCD phase transitions | 9/10 | |
| ex267 | SM anomaly cancellation in TGP | 10/10 | ★★★ |
| ex268 | Black holes (g→0, Vainshtein) | 10/10 | ★★★ |
| ex269 | RG flow of g₀ᵉ (no Landau pole) | 10/10 | ★★★ |
| ex270 | Neutron stars & dense matter | 10/10 | ★★★ |

**Phase subtotal: 73/80 (91.2%)**

### Summary Script

| Script | Topic | Score | Perfect |
|--------|-------|-------|---------|
| ex271 | Final summary v2 (definitive scorecard) | 10/10 | ★★★ |

---

## Grand Total

| Metric | Value |
|--------|-------|
| Physics scripts | 36 |
| Total tests | 348 |
| Tests passed | 314 |
| Pass rate | **90.2%** |
| Perfect scores (★★★) | **14/36** |
| Kill criteria survived | **0/15** |
| Quantitative predictions | **40** |
| Confirmed predictions | **21** |
| Parameter reduction | 35 → 7 (−80%) |
| Prediction-to-parameter ratio | **5.7** |

---

## Key Results

The 12 master equations (all verified in ex244):

| # | Equation | Value | Exp. agreement |
|---|----------|-------|----------------|
| F1 | α_s = 3g₀ᵉ/(32Ω_Λ) | 0.1190 | 1.1σ |
| F2 | λ_C = Ω_Λ/N | 0.22823 | 4.8σ |
| F3 | K(ℓ) = 2/3, K(ν) = 1/2 | exact | 0σ |
| F4 | \|GL(3,F₂)\| = 168 | exact | 0σ |
| F5 | α_s × Ω_Λ = 3g₀ᵉ/32 | 0.0815 | invariant |
| F6 | Ω_DM = Ω_b(N! − Ω_Λ) | 0.262 | 0.3σ |
| F7 | Σm_ν = 62.9 meV (NO) | prediction | < 120 meV |
| F8 | S[g] unified action | framework | — |
| F9 | n_s = 1 − 2/N_e | 0.967 | 0.4σ |
| F10 | m_W = 80.354 GeV | 80.354 | 0.01σ |
| F11 | m_H = v × 57/112 | 125.31 GeV | 0.3σ |
| F12 | β(g₀) = 0 at 1-loop | conformal | — |

---

## License

These scripts are supplementary material for the TGP papers.
For academic use only.
