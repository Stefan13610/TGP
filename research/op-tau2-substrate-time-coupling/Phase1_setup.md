---
title: "τ.2.Phase1 setup — scale-protection theorem (5 tests)"
date: 2026-04-30
cycle: τ.2.Phase1
status: SETUP
parent: "[[program.md]]"
tags:
  - TGP
  - tau2
  - phase1
  - setup
---

# τ.2.Phase1 setup

5 sub-tests deriving scale-symmetry protection of atomic masses + clock rates.

## Sub-tests

- **T1.1** Atomic-mass coupling candidates scan (4 forms: m·X^α, m+α·ln X, m·exp(α ln X), m + α(∂ ln X)²)
- **T1.2** Scale-invariance X→λX requires α=0 (m_atom INVARIANT)
- **T1.3** Noether scale-current ∂_μ J^μ = 0 implication for proper time
- **T1.4** Effective Hamiltonian H_atom = m c² + p²/(2m) + V(α_em) X-independent
- **T1.5** Protection theorem: clock rate ∝ |E_n - E_m|/ℏ X-invariant at leading

## Pass criteria

- 5/5 PASS = Phase 2 forward
- ≥4/5 PASS = Phase 2 conditional
- <4/5 = Phase 1 retry / τ.2 abort

## Output

`phase1_tau2_protection.py` runs all 5 sub-tests sequentially.
