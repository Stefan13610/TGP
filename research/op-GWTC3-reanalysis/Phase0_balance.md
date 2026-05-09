---
title: "Phase 0 balance sheet — op-GWTC3-reanalysis"
date: 2026-05-07
parent: "[[README.md]]"
type: phase0-balance
tgp_owner: research/op-GWTC3-reanalysis
tags:
  - phase0
  - balance-sheet
  - M03-gate
related:
  - "[[README.md]]"
  - "[[../../meta/CALIBRATION_PROTOCOL.md]]"
---

# Phase 0 balance sheet — op-GWTC3-reanalysis

## §1 — Inputy

### 1.1 Files read

| Plik | Wkład |
|------|-------|
| [[../op-ppE-mapping/Phase1_results.md]] | β_ppE^TGP^(b=-1) = -5/64 ≈ -7.81·10⁻² LOCKED |
| [[../op-ppE-mapping/Phase3_paper_ready.md]] | Multi-coefficient ratios (M911-P2 cross-check) |
| [[../op-LIGO-3G-deviation/Phase2_results.md]] | Fisher framework + degeneracy_factor=5 |
| [[../op-LIGO-3G-deviation/Phase3_falsifier_thresholds.md]] | Calibrated thresholds; LIGO-O3 ~10⁻¹ baseline |
| [[../../audyt/T01_LIGO3G_falsifier/FINDINGS.md]] §2.5 | Tier 5 recommendation |

### 1.2 Pre-existing locked anchors

| Anchor | Wartość | Status |
|--------|---------|--------|
| β_ppE^TGP^(b=-1) central | -5/64 ≈ -7.81·10⁻² | LOCKED (op-ppE-mapping Phase 1) |
| OOM window | [5.5·10⁻², 1.2·10⁻¹] | LOCKED |
| η (equal-mass reference) | 1/4 | LOCKED |
| G_SPA central | 1.0 | PRELIMINARY (~30% precision) |
| degeneracy_factor (Fisher) | 5 | LOCKED literature (Yagi-Yunes 2016) |

### 1.3 Literature inputs (kluczowe)

| Reference | Rola |
|-----------|------|
| Abbott et al. 2021 PRD 103:122002 (GWTC-2 ToGR) | Single-event + combined δφ̂_n bounds |
| Abbott et al. 2023 PRX 13:041039 (GWTC-3 ToGR) | Updated GWTC-3 bounds (~90 events) |
| Yunes & Pretorius 2009 PRD 80:122003 | ppE definitions + b_ppE = -1 mapping |
| Mishra et al. 2016 PRD 93:084054 | 2PN GR coefficient φ̂_4^GR (for fractional ↔ absolute) |
| Will 2014 LRR 17:4 | Fisher matrix conventions |

## §2 — Outputs

| Output | Plik | Epistemic class |
|--------|------|------------------|
| GWTC-3 bounds compilation (TGR Table III/IV) | [[Phase1_GWTC3_bounds.md]] | LITERATURE-COMPILATION |
| Unit conversion δφ̂_n ↔ β_ppE | [[Phase2_Bayes_factor.md]] §1 | NUMERICAL |
| TGP prediction overlay vs bounds | [[Phase2_Bayes_factor.md]] §2 | NUMERICAL |
| Bayes factor TGP vs GR (Laplace approx) | [[Phase2_Bayes_factor.md]] §3 | NUMERICAL-APPROXIMATE |
| Verdict + recommendation | [[Phase3_verdict.md]] | SYNTHESIS |

## §3 — Predyktywność ratio

| Metric | Wartość |
|--------|---------|
| N_locked_inputs | 1 (β_TGP) + 4 (δφ̂_n bounds, GWTC-3 + GWTC-2) = 5 |
| N_outputs | 1 conversion + 1 overlay + 1 Bayes factor + 1 verdict = 4 |
| Ratio | 4/5 = 0.8 (literature-side analysis, conservative) |

## §4 — Założenia

| ID | Założenie | Walidacja |
|----|-----------|-----------|
| G1 | Published δφ̂_n posteriors są approximately gaussian (allowing Laplace approximation) | YES — LIGO ToGR papers explicit gaussian-like centroid + 90% CL |
| G2 | TGP β_TGP^(b=-1) jest *single-coefficient* deviation (no spread w prior) | YES — sympy LOCK |
| G3 | Conversion δφ̂_n (fractional, LIGO convention) ↔ β_ppE (absolute, ppE convention) | NEED careful: φ̂_4^GR (2PN GR) coefficient z Mishra 2016 jest reference |
| G4 | Bayes factor Laplace approximation (gaussian peaks) jest sufficient dla ~10% precision | YES dla ~OOM verdict |
| G5 | GWTC-3 published bounds reprezentują full O1-O3 BBH catalog | YES (90 events post-cleaning) |

## §5 — Honest scope

**Co ten cykl POTRAFI:**
- Literature-grade reanalysis: published bounds + TGP prior overlay
- Bayes factor estimation (Laplace approximation, ~10% precision)
- Verdict z 4 opcji: (a) consistent + room, (b) tentative signal,
  (c) excluded, (d) within 1σ best-fit
- ~1 sesja Python work (numpy + scipy)

**Co ten cykl NIE robi:**
- Full reanalysis z LIGO H5 strain data (wymaga ~50+ godz + bilby/pycbc + GraceDB access)
- Per-event posterior re-derivation
- Multi-coefficient simultaneous Bayes (op-LIGO-3G-deviation Phase 4 future)
- New strain data analysis (only published summary statistics)

**Wartość mimo limitations:**
- Dostarcza **first-order verdict** TGP vs GR z aktualnych GWTC-3
- Identifikuje czy TGP jest sfalsyfikowane lub *signal-pre-publication*
- ~24-48 godz roadmap dla full reanalysis (jeśli pozytywny verdict)

## §6 — Phase 0 sign-off

- ✓ Inputy zinwentaryzowane.
- ✓ Outputy planowane z explicit epistemic class.
- ✓ Założenia G1-G5 znane i będą walidowane w Phase 1-2.
- ✓ Honest scope: literature-grade, NOT production-grade.
- ✓ M03 gate enforcement compliant.

**Phase 0 SIGNED 2026-05-07.** → Phase 1 GWTC-3 bounds compilation.
