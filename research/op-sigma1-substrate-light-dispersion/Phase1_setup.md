---
title: "σ.1.Phase1 setup — dispersion relation (5 tests)"
date: 2026-04-30
cycle: σ.1.Phase1
status: SETUP
parent: "[[program.md]]"
tags:
  - TGP
  - sigma1
  - phase1
  - setup
---

# σ.1.Phase1 setup

5 sub-tests deriving plane-wave dispersion relation from modified Maxwell ω.1.

## Sub-tests

- **W1.1** Plane-wave Fourier decomposition A_μ = a_μ e^{ik·x}
- **W1.2** Dispersion relation derivation: substituting do `∂_ν F^νμ = g F̃^μν ∂_ν(ln X)`
- **W1.3** Polarization eigen-modes: helicity basis e_±^μ = (e_1 ± i e_2)/√2
- **W1.4** Slowly-varying gradient validity: |∂(ln X)|/k² ≪ 1 (WKB)
- **W1.5** Gauge structure preservation under A_μ → A_μ + ∂_μ Λ

## Pass criteria

- 5/5 PASS = Phase 2 forward
- ≥4/5 PASS = Phase 2 forward z conditional
- <4/5 = Phase 1 retry / σ.1 abort

## Output

`phase1_sigma1_dispersion.py` runs all 5 sub-tests sequentially.
