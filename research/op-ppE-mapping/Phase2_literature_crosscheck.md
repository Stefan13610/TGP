---
title: "Phase 2 — Literature cross-check β_ppE^TGP^(b=-1) vs catalog modyfikacji GR"
date: 2026-05-07
parent: "[[README.md]]"
type: phase2-results
tgp_owner: research/op-ppE-mapping
tags:
  - phase2
  - literature
  - cross-check
  - dCS
  - sGB
  - Einstein-Aether
  - Brans-Dicke
  - ppE-catalog
related:
  - "[[README.md]]"
  - "[[Phase1_results.md]]"
  - "[[../../audyt/T01_LIGO3G_falsifier/PPN_TO_PPE_MAPPING.md]]"
---

# Phase 2 — Literature cross-check

## §0 — TL;DR

| Test | Wynik |
|------|-------|
| Single-coefficient `b_ppE^TGP = -1` jest TGP-distinguishing? | **NIE** — dzielony z dCS, sGB, Einstein-Æther |
| Multi-coefficient ratio pattern jest TGP-distinguishing? | **TAK** — 0 free parameters w TGP vs ≥1 free w innych |
| `|β_ppE^TGP|` magnitude w basenie LIGO-O3 GWTC-3 bounds? | **TAK** — borderline (TGP ~ 8·10⁻², bound ~ 10⁻¹) |
| Sygnał już w GWTC-3 reanalysis? | **MOŻLIWE** — re-analysis recommended (T01 NEEDS Q3 / FINDINGS Tier 5) |

## §1 — Catalog modyfikacji GR z 2PN-phase deviation (b_ppE = -1)

Z Yunes, Yagi, Pretorius 2016 (arXiv:1603.08955) i Yagi-Yunes 2016
(arXiv:1602.04674), modyfikacje GR dające phase deviation przy
b_ppE = -1 (2PN-phase):

### 1.1 Dynamical Chern-Simons (dCS) gravity

**Theory:** Quadratic curvature term `(α_dCS/4) ϑ R̃R` gdzie ϑ
jest pseudoscalarnym dilaton field i R̃R = ε^μνρσ R_μν^αβ R_ρσαβ.

**β_ppE^(dCS)^(b=-1) form:** Yagi-Yunes 2016 §5.2:
```
β_ppE^(dCS) = (1549225 / 11812864) · ζ_dCS · χ_dCS(η, χ)
```
gdzie ζ_dCS = α_dCS²/M⁴ (free coupling, dimensionless) i χ_dCS jest
mass+spin-zależnym prefactorem.

**Free parameters:** 1 (ζ_dCS).

**Constraint z GWTC-3:** ζ_dCS < O(1) (Yagi-Yunes 2016) → very weak.

**Distinction od TGP:**
- dCS β_ppE proportional do `ζ_dCS` (free coupling).
- TGP β_ppE proportional do `(5/6)` (FORCED by α=2 + hyperbolic).
- Distinguishing: ratio z 3PN-phase coefficient. dCS w 3PN ma
  dodatkowy prefactor zależny od ζ_dCS²; TGP ma fixed ratio -23/10.

### 1.2 Scalar-Gauss-Bonnet (sGB) gravity

**Theory:** `(α_sGB/4) f(ϑ) G_GB` z Gauss-Bonnet topological
invariant + scalar field ϑ.

**β_ppE^(sGB)^(b=-1) form:** Endlich-Khoury-Solomon 2017
(arXiv:1709.10018), Yagi-Yunes 2016:
```
β_ppE^(sGB) = -(5/7168) · ζ_sGB · ξ(η)
```
z ζ_sGB = 16π α_sGB²/M⁴ (free coupling).

**Free parameters:** 1 (ζ_sGB).

**Distinction:** podobny problem jak dCS — ratio 3PN/2PN zależy
od ζ_sGB i postaci f(ϑ).

### 1.3 Einstein-Æther theory

**Theory:** Vector-tensor with timelike vector `u^μ` z Lorentz
violation parameters c_1, c_2, c_3, c_4.

**β_ppE^(EÆ)^(b=-1) form:** Yagi-Eling-Foster 2014
(arXiv:1311.7144):
```
β_ppE^(EÆ) = -3/(128 η) · (1 - c_T/c_s) · F_EAE(c_1, c_2, c_3, c_4, η)
```
z 4 free parameters c_i.

**Free parameters:** 4 (c_1, c_2, c_3, c_4).

**KEY: TGP F1 (single Φ) wyklucza dodatkowy vector field** —
Einstein-Æther jest strukturalnie incompatible z TGP foundations.

### 1.4 Brans-Dicke higher-PN

**Theory:** `S_BD = (1/16π) ∫ d⁴x √(-g) (φ R - ω/φ ∂φ²)` z
1 free parameter ω.

**β_ppE^(BD)^(b=-1) form:** lower-PN dipole (b=-7) jest dominant;
2PN phase appears jako sub-leading. Yagi-Yunes 2016 daje:
```
β_ppE^(BD)^(2PN) ≈ -(5/3584) · (1/(2ω+3)) · χ_BD(η)
```

**Free parameters:** 1 (ω).

**KEY: M9.1'' jest *strukturalnie* skalar-tensorowa** ale z **0 free
parameters** (α=2 LOCKED, hyperbolic LOCKED). To jest *bardziej
restryktywna* niż BD — restorable do BD limit by α → 0 i odpowiedni
ansatz f(ψ).

### 1.5 Tensor: porównanie

| Theory | Free params | b_ppE = -1 | Ratio do 3PN | Distinguishing od TGP |
|--------|-------------|------------|---------------|-----------------------|
| **TGP M9.1''** | **0** | yes (5/6 forced) | -23/10 (forced) | reference |
| dCS | 1 (ζ_dCS) | yes | depends ζ | ratio z 3PN różna |
| sGB | 1 (ζ_sGB) | yes | depends ζ | ratio z 3PN różna |
| Einstein-Æther | 4 (c_i) | yes | depends c_i | strukturalnie różna (vector field) |
| Brans-Dicke | 1 (ω) | yes (sub-leading) | depends ω | dominant b=-7 (TGP brak dipole) |

**Konkluzja:** **multi-coefficient ratio pattern** jest TGP-distinguishing
signature. Single coefficient detection (b=-1) NIE jest distinguishing
— wymaga 2-3 coefficients dla unique identification.

## §2 — Comparison z aktualnymi LIGO bounds

### 2.1 GWTC-3 ToGR ppE constraints (Abbott et al. 2021, 2023)

LIGO Tests of GR (TGR) papers giving ppE bounds:

| Reference | ppE order | Bound (90% CL) | Sample |
|-----------|-----------|----------------|--------|
| Abbott et al. 2021 (PRD 103:122002, GWTC-2) | 0PN-3.5PN phase | δφ_n < O(10⁻¹) per coef | 13 BBH events |
| Abbott et al. 2023 (Phys. Rev. X 13:041039, GWTC-3 ToGR) | refined | similar OOM, slight improvements | 35 events |

**Specific 2PN-phase (b=-1) bound dla GW150914-like:** ~10⁻¹
(Yunes-Yagi-Pretorius 2016 single-event Fisher).

**TGP β_ppE^TGP^(b=-1) ≈ 8 · 10⁻²** (Phase 1 §6.3 lock; OOM 5.5·10⁻²
do 1.2·10⁻¹).

**Stosunek bound/prediction:**
- Single-event LIGO-O3: bound ~10⁻¹, TGP ~8·10⁻². **TGP jest
  borderline** (margin ~1.3×, *NIE* sfalsyfikowane).
- Stack 100 BBH O3+O4: bound ~10⁻², TGP ~8·10⁻². **DECISIVE**.

### 2.2 Side-cykl rekomendowany: GWTC-3 reanalysis dla M911-P1

T01 NEEDS Q3 / FINDINGS Tier 5: re-analysis aktualnych GWTC-3 ppE
constraints w 2PN-phase coefficient z TGP-specific β_ppE^TGP =
-5/64 prior (zamiast generic ppE prior).

**Output spodziewany:**
- Bayes factor TGP vs GR z aktualnych GWTC-3 events.
- Możliwość: TGP już *wykryte* w O3 reanalysis (z low significance
  ~1-2σ); lub TGP *bordering* falsification.

**Estymowana praca:** ~1 sesja Python (bilby + ppE prior wstawienie).

## §3 — Założenia A1–A6 walidacja

Z [[../../audyt/T01_LIGO3G_falsifier/PPN_TO_PPE_MAPPING.md]] §5:

| ID | Założenie | Status w Phase 1+2 |
|----|-----------|---------------------|
| A1 | M9.1'' two-body Lagrangian istnieje | **VALIDATED** w Phase 1 §4 (linear superposition weak-field) |
| A2 | Quadrupole formula struktura zachowana, deviation O((5/6)U³) wchodzi przez metric perturbation | **PRELIMINARY** — założone w Phase 1 §6.2 (G_SPA = 1 OOM); wymaga retarded Green's function calc w M9.1'' (Phase 1.5 future) |
| A3 | dE/dt luminosity modyfikuje 2PN+ orbital | **PRELIMINARY** — propagates from A2 |
| A4 | SPA stacjonarna w M9.1'' (adiabatic inspiral) | **VALIDATED** — adiabatic condition spójna w cały band ET (f < f_ISCO) |
| A5 | TGP nie wprowadza nowych radiacyjnych DOF | **VALIDATED** w PREDICTIONS_REGISTRY GW1, GW2 (3 DOF, c_T = c_s) |
| A6 | Konwencja PN counting | **VALIDATED** w [[../../audyt/T01_LIGO3G_falsifier/CONVENTION_DECISION.md]] (PHASE primary) |

**4/6 VALIDATED, 2/6 PRELIMINARY** (A2, A3 wymagają tighter lock).

Phase 1.5 (future cycle continuation) zamknie A2/A3 do precyzji 5%.

## §4 — Output Phase 2 → Phase 3

Phase 3 ([[Phase3_paper_ready.md]]) syntezuje Phase 1 + Phase 2:
- Tabela porównawcza β_ppE^TGP vs catalog (5 theories).
- Multi-coefficient TGP-signature plot (do ngEHT/LISA presentation).
- Falsifier statement liczbowy zaktualizowany dla T01 wpisu.
