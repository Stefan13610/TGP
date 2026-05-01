---
title: "UV.2.Phase2 setup — sympy LOCK + numerical M_TGP reproduction + UV-IR cascade"
date: 2026-05-01
cycle: UV.2.Phase2
status: SETUP
parent: "[[program.md]]"
predecessor: "[[Phase1_results.md]]"
tags:
  - TGP
  - UV2
  - phase2
  - sympy-lock
  - M_TGP
  - setup
---

# UV.2.Phase2 setup

**Score gate:** ≥6/7 PASS = Phase 2 forward.

## Sub-tests (7)

- **U2.1** TGP dim-less invariant ledger: enumeruj wszystkie LOCKED dim-less
  invariants {g*, λ*, η_N*, N_A, c_χ, α₀, g̃, κ, F4, F5, XS1} i ich struktualne
  algebraic combinations (level-1 products/quotients, level-2 with π factors).
- **U2.2** K_struct sympy form derivation: w 4-cand scan (Phase 1) wybór K_a
  = N_A · 2π² uniqueness w drift < 1% partial-lock band.
- **U2.3** M_Pl/M_GUT prediction: sympy `M_Pl/M_GUT = 2π²·N_A^(3/2)/√g*`
  vs PDG·M_GUT obs. **Drift target:** < 1%.
- **U2.4** M_TGP numerical reproduction: M_TGP_UV.2 = K·M_GUT vs χ.1.
  **Drift target:** < 0.5% (within M_GUT 2-loop SM uncertainty).
- **U2.5** G_N → M_Pl chain post-UV.2: G_N predicted z M_GUT only;
  drift vs CODATA. **Drift target:** < 0.5% (within M_GUT unc.).
- **U2.6** UV-IR cascade self-consistency: AS NGFP UV cutoff Λ_AS → M_TGP IR
  scale matching consistent z η_N* = -2 marginal IR (sub-percent residual).
- **U2.7** F-cluster post-UV.2: F4/F5/F6 untouched (algebraic + EFT-scaling
  + κ-DERIVED) ↔ XS1 √α₀ = κ_TGP preserved drift < 1%.

## Inputs from Phase 1 (LOCKED)

```
K_struct LOCK   = N_A · 2π² ≈ 173.15
M_Pl/M_GUT pred = 2π²·N_A^(3/2)/√g* ≈ 608.62
M_GUT (input)   = 2.00·10¹⁶ GeV (SM 2-loop, ~10-30% theory unc.)
g*              = 71/100 (UV.1 NGFP)
N_A             = 500/57 (ξ.1)
α₀              = 1069833/264500 (F4)
κ_TGP           = 2.012 (XS1)
g̃               = 0.9803 (F5)
M_Pl PDG        = 1.220890·10¹⁹ GeV
M_TGP χ.1       = 3.4734·10¹⁸ GeV (chi.1 anchor)
G_N CODATA      = 6.67430·10⁻¹¹ m³ kg⁻¹ s⁻²
```

## Strategy

**Sympy-rational LOCK**: K = (500/57) · 2π² formal-rational form (z transcendental
2π² factor); reproduces target K_target = M_TGP/M_GUT obs ≈ 173.67 at drift 0.30%.

**Joint-system propagation**: K LOCK + χ.1 G_N relation + M_GUT anchor →
predicts M_TGP, M_Pl, G_N(SI) wszystkie w jednym dim-less framework
(no M_Pl PDG input).

**Falsifier**: jeśli M_GUT 2-loop SM running drifts > 30% post-threshold corrections,
UV.2 K_struct hypothesis falsified (M_GUT no longer reliable single anchor).

## Cross-references

- [[program.md]]
- [[Phase1_results.md]]
