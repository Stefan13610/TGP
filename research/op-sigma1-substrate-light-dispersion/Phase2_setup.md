---
title: "σ.1.Phase2 setup — phase/group velocity LOCK (7 tests)"
date: 2026-04-30
cycle: σ.1.Phase2
status: SETUP
parent: "[[program.md]]"
tags:
  - TGP
  - sigma1
  - phase2
  - sympy
  - velocity
---

# σ.1.Phase2 setup

7 sub-tests sympy-LOCKING phase/group velocity + effective optical metric.

## Sub-tests

- **W2.1** Phase velocity v_φ_± = ω_±/k via sympy LOCK
- **W2.2** Group velocity v_g_± = ∂ω_±/∂k
- **W2.3** Polarization-averaged c_eff = (v_+ + v_-)/2 = 1 (no scalar c(X))
- **W2.4** Birefringence Δv_φ = v_+ − v_− = g·n_∥/k (where n_∥ ≡ k̂·∇ln X)
- **W2.5** Helicity-dependent optical cone Q^(±)_μν(p) k^μ k^ν = 0 (effective dispersion geometry, NOT full classical metric)
- **W2.6** Scale-invariance preservation X→λX (Δv invariant)
- **W2.7** ω.1 EOM back-reaction consistency

## Pass criteria

- 7/7 = FULL CASCADE → Phase 3 forward
- ≥5/7 = Phase 3 forward
- <5/7 = Phase 2 retry / σ.1 abort

## Output

`phase2_sigma1_velocity.py` runs sympy verification on all 7 sub-tests.
