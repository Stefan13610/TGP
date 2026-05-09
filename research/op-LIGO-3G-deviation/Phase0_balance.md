---
title: "Phase 0 balance sheet — op-LIGO-3G-deviation"
date: 2026-05-07
parent: "[[README.md]]"
type: phase0-balance
tgp_owner: research/op-LIGO-3G-deviation
tags:
  - phase0
  - balance-sheet
  - M03-gate
  - op-LIGO-3G-deviation
related:
  - "[[README.md]]"
  - "[[../../meta/CALIBRATION_PROTOCOL.md]]"
---

# Phase 0 balance sheet — op-LIGO-3G-deviation

## §1 — Inputy

| Plik | Wkład |
|------|-------|
| [[../op-ppE-mapping/Phase1_results.md]] | β_ppE^TGP^(b=-1) = -5/64 ≈ -7.81·10⁻² LOCKED (sympy 14/14) |
| [[../op-ppE-mapping/Phase3_paper_ready.md]] | Multi-coefficient ratios {-23/10, -38/23, +337/228} |
| [[../../audyt/T01_LIGO3G_falsifier/CYCLE_KICKOFF_op-LIGO-3G-deviation.md]] | Phase 1-3 setup recipe, gates G1-G6 |
| [[../../audyt/T01_LIGO3G_falsifier/SENSITIVITY_BACK_OF_ENVELOPE.md]] | OOM bounds z literatury (Chamberlain-Yunes 2017, Maggiore 2020) |

## §2 — Outputs

| Output | Plik | Epistemic class |
|--------|------|------------------|
| TaylorF2 + ppE waveform model | [[Phase1_waveform_setup.md]] | NUMERICAL |
| Fisher matrix single-event ET-D, CE, LIGO-O5 | [[Phase2_results.md]] §2 | NUMERICAL |
| Stacked SNR forecasts | [[Phase2_results.md]] §3 | NUMERICAL |
| Degeneracy analysis (M_chirp, η, χ_eff × β) | [[Phase2_results.md]] §4 | NUMERICAL |
| Falsifier thresholds tabela | [[Phase3_falsifier_thresholds.md]] | SYNTHESIS |

## §3 — Predyktywność ratio

| Metric | Wartość |
|--------|---------|
| N_locked_inputs | 1 (β_ppE^TGP from Path B) + 4 (ASD curves z literatury) = 5 |
| N_outputs | 4 detector × 2 (single + stack) = 8 SNR thresholds + 1 degeneracy + 1 multi-coef extension = 10 |
| Ratio | 10/5 = 2.0 (forecasting cycle, multi-detector) |

## §4 — Założenia

| Założenie | Walidacja |
|-----------|-----------|
| TaylorF2 baseline reprodukuje GR Fisher analyses do 5% | LIGO public software validated |
| ppE single-coefficient deviation OK w `inspiral` (f < f_ISCO) | tak, weak-field assumption |
| Sky-averaging (4π/5) i inclination-averaging | standard w forecasting |
| ASD curves z public sources (T1800042, ET-D Maggiore, CE Reitze) | citations |

## §5 — Phase 0 sign-off

- ✓ Inputy zinwentaryzowane (β_ppE^TGP locked).
- ✓ Outputy planowane.
- ✓ Predyktywność ratio 2.0 (forecasting).
- ✓ M03 gate enforcement compliant.

**Phase 0 SIGNED 2026-05-07.** → Phase 1 setup.
