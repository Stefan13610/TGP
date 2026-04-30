---
title: "τ.3.Phase3 setup — predictions + 4-channel convergence (6 tests)"
date: 2026-04-30
cycle: τ.3.Phase3
status: SETUP
parent: "[[program.md]]"
tags:
  - TGP
  - tau3
  - phase3
  - setup
  - predictions
---

# τ.3.Phase3 setup

6 sub-tests deriving falsifiable lab predictions across 4 cross-channels.

## Sub-tests

- **T3.1** Sr/Yb 1e-18/yr lab differential E∥B vs E⊥B prediction:
  current sensitivity sets Λ > X bound depending on null-result.
- **T3.2** ELI-NP / HERMES 2030+ frontier (10²² W/cm² laser, B ~ 10⁹ G):
  E·B vastly larger → δω/ω boost ~ (E·B_new/E·B_old)² → improved Λ reach.
- **T3.3** Cosmological residual: primordial B fields ~ 1 nG induce sub-leading
  α_em(z) scatter beyond Webb 1e-7 — possible τ.3 signature in next-gen QSO.
- **T3.4** Magnetar atmosphere F·F̃ regional clock acceleration: Chandra/NICER
  atomic spectroscopy in pole regions where E∥B → δω/ω ≠ 0, observable as 
  spectral line shift correlated z magnetic pole geometry.
- **T3.5** 4 alt-L4-couplings cross-channel falsification: predictions for
  β F·F̃ (linear), γ ln(F·F̃), δ (E²-B²), η (E·B)² couplings — distinguishable
  via signature scaling (B² vs B⁴ vs E·B sign).
- **T3.6** 4-channel τ.3 convergence: lab + frontier + cosmo + magnetar
  jointly constrain Λ and α_g sign.

## Pass criteria

- 6/6 PASS = τ.3 program END
- ≥5/6 PASS = τ.3 conditional close
- <5/6 = Phase 3 retry / τ.3 weak close

## Output

`phase3_tau3_predictions.py` runs all 6 sub-tests sequentially.
