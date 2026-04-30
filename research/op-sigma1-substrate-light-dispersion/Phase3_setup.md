---
title: "σ.1.Phase3 setup — predictions + 4-channel convergence (6 tests)"
date: 2026-04-30
cycle: σ.1.Phase3
status: SETUP
parent: "[[program.md]]"
tags:
  - TGP
  - sigma1
  - phase3
  - predictions
  - falsification
---

# σ.1.Phase3 setup

6 sub-tests deriving experimental predictions + 4-channel convergence.

## Sub-tests

- **W3.1** Lab phase-velocity measurement (Mach-Zehnder + B field)
- **W3.2** Pulsar polarized dispersion residuals (PSR J0030/J0740, NICER + IXPE)
- **W3.3** CMB E/B-mode chirality (Planck PR4 + ACT 2024 + SO 2027+ + LiteBIRD 2029+)
- **W3.4** Atomic clock ratio gradient sensitivity (Hg/Yb/Sr clock comparisons)
- **W3.5** Alt-dispersion cross-channel falsification (scalar c(X), tensor c(X), Lorentz-violating)
- **W3.6** 4-channel σ.1 convergence

## Pass criteria

- 6/6 = FULL CONVERGENCE → σ.1 program END
- ≥5/6 = σ.1 PASS
- <5/6 = retry

## Output

`phase3_sigma1_predictions.py` runs all 6 prediction tests.
