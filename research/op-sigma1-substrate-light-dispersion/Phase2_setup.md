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
- **W2.4** Birefringence Δv = v_+ - v_- = g·∂(ln X)/k
- **W2.5** Effective optical metric g_μν^opt = η_μν + δg_μν(∂ ln X)
- **W2.6** Scale-invariance preservation X→λX (Δv invariant)
- **W2.7** ω.1 EOM back-reaction consistency

## Pass criteria

- 7/7 = FULL CASCADE → Phase 3 forward
- ≥5/7 = Phase 3 forward
- <5/7 = Phase 2 retry / σ.1 abort

## Output

`phase2_sigma1_velocity.py` runs sympy verification on all 7 sub-tests.
