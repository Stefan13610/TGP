---
title: "τ.3.Phase1 setup — L4 coupling structural derivation (5 tests)"
date: 2026-04-30
cycle: τ.3.Phase1
status: SETUP
parent: "[[program.md]]"
tags:
  - TGP
  - tau3
  - phase1
  - setup
  - L4
---

# τ.3.Phase1 setup

5 sub-tests deriving structural properties of L4 gradient-coupled mass and viability.

## Sub-tests

- **T1.1** L4 candidate forms scan (sign-of-α_g, magnitude, alt variants):
  L4_a = m_0 + α_g(∂ ln X)², L4_b = m_0 + β F·F̃, L4_c = m_0 + γ ln(F·F̃), 
  L4_d = m_0 + η (E·B)².
- **T1.2** UV matching candidates: AS NGFP fixed-point Wilson coefficient sign;
  decoupling theorem heavy-mode integration; cosmological consistency of α_g sign.
- **T1.3** Effective mass derivation: m_e_eff(X) = m_e^(0) + (α_g/Λ²)(∂_μ ln X)(∂^μ ln X)
  expanded to leading order in (∂ ln X / Λ).
- **T1.4** Substrate gradient source: ω.1 EOM gives □(ln X) = -(g/f_X²) E·B;
  parallel E∥B in lab → max source; perpendicular → zero source.
- **T1.5** Viability gate: detectable iff Λ ≲ 100 MeV (1e-18/yr Sr threshold);
  upper Λ bound from null current Sr/Yb differential clock test.

## Pass criteria

- 5/5 PASS = Phase 2 forward
- ≥4/5 PASS = Phase 2 conditional
- <4/5 = Phase 1 retry / τ.3 abort

## Output

`phase1_tau3_l4coupling.py` runs all 5 sub-tests sequentially.
