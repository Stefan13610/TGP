---
title: "Phase 2 — Bayes factor TGP vs GR (Laplace approximation, GWTC-3)"
date: 2026-05-07
parent: "[[README.md]]"
type: phase2-results
tgp_owner: research/op-GWTC3-reanalysis
tags:
  - phase2
  - Bayes-factor
  - Laplace
  - TGP-vs-GR
  - GWTC-3
  - reconciliation
related:
  - "[[README.md]]"
  - "[[Phase0_balance.md]]"
  - "[[Phase1_GWTC3_bounds.md]]"
  - "[[scripts/bayes_factor_tgp_vs_gr.py]]"
  - "[[scripts/bayes_factor_tgp_vs_gr.txt]]"
---

# Phase 2 — Bayes factor TGP vs GR

## §0 — TL;DR

```
TGP M911-P1 prediction in LIGO units:
   δφ̂_4_TGP_central = -0.018  (z β_TGP^(b=-1) = -5/64 / 4.336)
   OOM: δφ̂_4_TGP ∈ [-0.0277, -0.0127]

vs GWTC-3 combined posterior (~90 BBH):
   δφ̂_4_obs = +0.050 ± 0.182 (1σ)

Deviation: TGP is at 0.37σ from observed centroid → CONSISTENT within 1σ.

Bayes factor: BF_TGP/GR ≈ 0.97 → INCONCLUSIVE
   (TGP and GR equally well-supported by current data)

Detection power: TGP signal magnitude / GW uncertainty = 0.10
   → ~10⁵ BBH events needed for 5σ confirmation in current paradigm
   → multi-coefficient marginalized Fisher; single-coefficient is much faster
```

## §1 — TGP prediction in LIGO fractional units

Z [[Phase1_GWTC3_bounds.md]] §4 conversion:

```
β_ppE^(b=-1) = 4.336 · δφ̂_4    (η=1/4)

TGP central:   β_TGP = -5/64 = -0.0781    →   δφ̂_4_TGP = -0.0180
TGP OOM low:   |β| = 0.055               →   |δφ̂_4| = 0.0127
TGP OOM high:  |β| = 0.120               →   |δφ̂_4| = 0.0277
```

Sign convention: TGP daje *deeper* potential well niż GR (Δα_3 = -5/6 < 0),
więc phase advance → δφ̂_4 < 0.

## §2 — TGP overlay vs GWTC-3 published posteriors

| Scenariusz | δφ̂_4_obs | σ_obs (1σ) | δφ̂_4_TGP | (TGP - obs)/σ | Within range |
|------------|------------|--------------|------------|------------------|--------------|
| GWTC-3 combined (~90 BBH) | +0.050 | 0.182 | -0.018 | -0.37 | **YES (1σ)** |
| GW150914-like single | +0.200 | 0.608 | -0.018 | -0.36 | YES (1σ) |
| GWTC-2 combined | +0.100 | 0.274 | -0.018 | -0.43 | YES (1σ) |

**Konkluzja §2:** TGP M911-P1 prediction jest **consistent z GWTC-3
within 1σ** — żadne **tension**, żadne **rejection**.

## §3 — Bayes factor (Laplace approximation)

Standard Laplace approximation: dla gaussian posterior z mean θ_obs
i width σ_obs, likelihood ratio:

```
BF_TGP/GR = exp[-(δφ_TGP - δφ_obs)² / (2σ²)] / exp[-(0 - δφ_obs)² / (2σ²)]
          = exp[(δφ_obs² - (δφ_obs - δφ_TGP)²) / (2σ²)]
```

| Scenariusz | BF_TGP/GR | log10(BF) | Verdict |
|------------|------------|-----------|---------|
| GWTC-3 combined (~90 BBH) | 0.969 | -0.014 | **INCONCLUSIVE** |
| GW150914-like single | 0.990 | -0.004 | INCONCLUSIVE |
| GWTC-2 combined | 0.974 | -0.011 | INCONCLUSIVE |

**Jeffreys interpretation:**
- |log10(BF)| < 0.5 → **INCONCLUSIVE** (no preference between models)
- 0.5 < log10(BF) < 1 → weak preference
- log10(BF) > 1 → strong preference (>10× preferred)

**Wszystkie scenariusze są w zakresie INCONCLUSIVE** — GWTC-3 nie
preferuje ani TGP, ani GR z aktualnych danych. Ta jest expected:
TGP signal (δφ̂_4 ~ -0.018) jest ~10× *poniżej* current GW
posterior width (~0.18 1σ).

## §4 — Detection power: ile zdarzeń potrzeba dla 5σ?

```
Current TGP signal-to-noise ratio: |δφ̂_4_TGP| / σ_GWTC3 = 0.018 / 0.182 = 0.10

Naive scaling: σ_combined ∝ 1/√N for stack of N events.
Want 5σ detection: 0.018 / σ_required = 5 → σ_required = 0.0036

Current GWTC-3 ~90 events: σ ≈ 0.182
For σ = 0.0036: factor (0.182/0.0036)² = 2560× more events
→ N_events_needed = 90 · 2560 = ~230,000 BBH events
```

## §5 — Reconciliation z Path A Phase 2-3

### 5.1 Tension w wnioskach

[[../op-LIGO-3G-deviation/Phase3_falsifier_thresholds.md]] §2
zawiera (calibrated):
- LIGO-O5 single ~3·10⁻² → TGP/β_5σ ≈ 2.6 → "first decisive 2027+"
- ET-D single ~10⁻² → TGP/β_5σ ≈ 7.8 → "decisive >5σ"
- CE single ~3·10⁻³ → TGP/β_5σ ≈ 26 → "decisive >50σ"

Ten Phase 2 (GWTC-3 reanalysis) zawiera:
- Current GWTC-3 (~90 BBH stack) → TGP signal/σ = 0.10 → INCONCLUSIVE
- Need ~230,000 BBH events for 5σ → faaar od ET-D + CE projected ~10⁵/yr

**To nie zgadza się** — GWTC-3 reanalysis jest *factor ~50× mniej
optimistyczna* niż Path A.

### 5.2 Pochodzenie tension: jednoznaczny vs marginalized Fisher

| Path A (calibrated, optimistic) | Phase 2 (GWTC-3, conservative) |
|------------------------------------|-----------------------------------|
| **Single-coefficient Fisher** dla `β_ppE^(b=-1)` z degeneracy_factor=5 dla intrinsic params (M_chirp, η, χ_eff) | **Multi-coefficient marginalized** Fisher z LIGO ToGR style — fitting wszystkich 8 PN coefficients (0PN-3.5PN) jednocześnie |
| `dΨ/dβ_ppE = u^(-1)` (single derivative) | All `dΨ/dα̂_n` derivatives correlated → bounds significantly weaker per coefficient |
| Bound β_5σ ~10⁻³ (ET-D single, calibrated) | Bound δφ̂_4 ~0.18 (1σ GWTC-3 combined) → β_5σ ~3.9 absolute |

**To są DWA RÓŻNE tests**, nie konflikt:

- **Path A:** "Pure ppE-deviation" prior — assume TGP signal is present
  in single coefficient, fit absolute β.
- **Phase 2 / GWTC-3 ToGR:** "Generic PN deviation" prior —
  marginalize nad wszystkimi 8 PN coefficients simultaneously.

### 5.3 Honest implication

**TGP M9.1'' jest distinguishable** z 3G era (ET-D + CE) **TYLKO**
w specific scenariuszach:
1. **Theory-priored Bayes inference** (single-coefficient β_TGP prior z
   TGP-specific value -5/64): lock predyktywności wymaga TGP-specific
   prior, nie generic ppE.
2. **Multi-coefficient ratio test** (M911-P2): testuje pattern
   {-23/10, -38/23, +337/228}, nie absolute coefficients. To jest
   **STRONGER** test, bo nie wymaga wykrycia każdego coefficient
   individually.

**TGP NIE jest distinguishable** w generic ppE marginalized analysis
(jak LIGO ToGR papers) — bounds są zbyt słabe.

### 5.4 Konkluzja Reconciliation

Path A Phase 3 detection thresholds **POZOSTAJĄ VALID** dla
TGP-specific analysis (single-coefficient β prior, lub multi-coefficient
ratio test). ALE: pre-publication signal w generic GWTC-3 ToGR papers
jest **NIE-DOSTĘPNY** — required jest dedicated TGP analysis pipeline.

## §6 — Bayes factor z OOM uncertainty (G_SPA)

TGP β ma 30% OOM uncertainty z G_SPA. Marginalizacja:

| Scenariusz | BF central | BF (best in window) | BF (worst) |
|------------|-------------|----------------------|--------------|
| GWTC-3 combined | 0.969 | 0.979 | 0.948 |
| GW150914-like single | 0.990 | 0.993 | 0.984 |
| GWTC-2 combined | 0.974 | 0.982 | 0.959 |

**Wszystkie pozostają INCONCLUSIVE.** OOM uncertainty G_SPA NIE zmienia
verdict.

## §7 — Output → Phase 3

Phase 3 ([[Phase3_verdict.md]]) syntezuje:
- Verdict: TGP consistent z GWTC-3 within 1σ; nie sfalsyfikowane,
  nie potwierdzone.
- Reconciliation: Path A optimistic detection thresholds wymagają
  *dedicated* TGP-specific analysis (nie generic ppE).
- Recommendations dla autora: jaki workflow dla pre-publication
  signal jest *actually* viable.
