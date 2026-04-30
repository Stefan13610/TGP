---
title: "τ.2.Phase3 setup — predictions + 4-channel convergence (6 tests)"
date: 2026-04-30
cycle: τ.2.Phase3
status: SETUP
parent: "[[program.md]]"
tags:
  - TGP
  - tau2
  - phase3
  - predictions
  - falsification
---

# τ.2.Phase3 setup

6 sub-tests deriving experimental predictions + 4-channel convergence.

## Sub-tests

- **T3.1** Atomic clock cosmological drift NULL (Webb/Murphy 2003-2017 1e-7 precision)
- **T3.2** Lab Hg/Yb/Sr clock comparison precision (current ~1e-18/yr, future 1e-22/yr 2035+)
- **T3.3** Strong-gradient residuals: magnetar atomic spectroscopy + lab E·B clock
- **T3.4** Polarization-Zeeman cross-coupling z σ.1 (circularly polarized drive)
- **T3.5** Alt-clock-couplings cross-channel falsification (4 forms: m·X^α, ℏ·X^β, α_em·X^γ, hyperfine·X^δ)
- **T3.6** 4-channel τ.2 convergence

## Pass criteria

- 6/6 = FULL CONVERGENCE → τ.2 program END
- ≥5/6 = τ.2 PASS
- <5/6 = retry

## Output

`phase3_tau2_predictions.py` runs all 6 prediction tests.
