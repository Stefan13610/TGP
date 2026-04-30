---
title: "τ.2.Phase2 setup — proper-time + sympy LOCK (7 tests)"
date: 2026-04-30
cycle: τ.2.Phase2
status: SETUP
parent: "[[program.md]]"
tags:
  - TGP
  - tau2
  - phase2
  - sympy
  - proper-time
---

# τ.2.Phase2 setup

7 sub-tests sympy-LOCKING proper time + clock-rate ratio + cross-coupling z σ.1.

## Sub-tests

- **T2.1** Proper time τ = ∫√(g_00^eff) dt sympy LOCK z g_00^eff = 1 + ε(∂ ln X)²
- **T2.2** Substrate-induced time dilation δτ/τ at leading O(∂ ln X) = 0
- **T2.3** Polarization-Zeeman cross-coupling: drive ω_drive split via σ.1 → atomic Zeeman level shift
- **T2.4** Cs hyperfine (~9.2 GHz) vs optical Sr/Yb (~430 THz) relativistic invariance protection
- **T2.5** Clock-rate ratio R(X1)/R(X2) = 1 sympy EXACT under X→λX
- **T2.6** ω.1 + σ.1 consistency: scalar c(X) NULL → scalar atomic frequency NULL
- **T2.7** Alt-couplings m_e ∝ X^α breaks scale-symmetry FALSIFIED

## Pass criteria

- 7/7 = FULL CASCADE → Phase 3 forward
- ≥5/7 = Phase 3 forward
- <5/7 = Phase 2 retry / τ.2 abort

## Output

`phase2_tau2_propertime.py` runs sympy verification on all 7 sub-tests.
