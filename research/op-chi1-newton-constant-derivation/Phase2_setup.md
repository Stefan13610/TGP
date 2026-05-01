---
title: "χ.1.Phase2 setup — sympy LOCK + JOINT M_TGP anchor + numerical κ reproduction"
date: 2026-05-01
cycle: χ.1.Phase2
status: SETUP
parent: "[[program.md]]"
predecessor: "[[Phase1_results.md]]"
tags:
  - TGP
  - chi1
  - phase2
  - sympy-lock
  - M_TGP-joint
  - newton-constant
  - setup
---

# χ.1.Phase2 setup

**Score gate:** ≥6/7 PASS = Phase 2 forward.

## Sub-tests (7)

- **X2.1** AS NGFP RG-flow integration: solve `dg/d ln k = (η_N* + 2)g + ...`
  z η_N*=-2 → `G(k)·k² = g*` marginal limit; threshold matching przy
  `k = M_TGP`: **G_N · M_TGP² = g*** sympy-exact (UV.1 anchor; ξ_grav=1 baseline)
- **X2.2** ξ_grav structural form scan (4 candidates):
  - (a) ξ_grav = 1 (trivial / marginal limit)
  - (b) ξ_grav = N_A = 500/57 (ξ.1 photon-ring inheritance)
  - (c) ξ_grav = (X_ref/X_obs)^(2/N_gen) (υ.1 closure law)
  - (d) ξ_grav = κ_TGP² = 4.048 (XS1 cross-sector — disfavored)
- **X2.3** **JOINT M_TGP lock** z 3 orthogonal-anchor constraints:
  - (1) GUT scale matching: M_TGP ≈ M_GUT ~ 2·10¹⁶ GeV (gauge unification baseline)
  - (2) M_Pl natural relation: M_TGP = M_Pl · √(g*/ξ_grav) (χ.1 hypothesis inversion)
  - (3) Phase 2.E.3 g̃ entropy scaling — orthogonal cross-check
- **X2.4** Numerical κ reproduction: plug values → predict κ = √(32π G_N)
  in M_Pl=1 units. **Drift target vs F6 anchor 10.0265: < 0.1%**
- **X2.5** G_N → M_Pl chain: M_Pl = G_N^(-1/2); predict M_Pl ± drift vs
  PDG 1.220890·10¹⁹ GeV. **Drift target: < 0.1%**
- **X2.6** Quantum-gravity self-consistency: 1-loop graviton self-energy
  match z F5 (g̃ = 0.9803); 2-loop FRG (UV-research-track) preserves G(k)
  marginal scaling
- **X2.7** Cross-sector consistency: F4 (α₀) ↔ F6 (κ post-χ.1) — verify
  √α₀ = κ_TGP identity (XS1) doesn't depend on G_N reinterpretation

## Inputs from Phase 1 (LOCKED)

```
Hypothesis: G_N = g* / (M_TGP² · ξ_grav)
c_χ        = √3                        (canonical kinetic match)
g*         = 0.71                       (UV.1 AS NGFP, η_N* = -2 marginal IR)
N_A        = 500/57 ≈ 8.7719           (ξ.1 photon-ring sympy-exact)
F6 κ       = √(32π) ≈ 10.0265           (G_N=1 natural units, F6 STRUCTURAL)
M_Pl       = 1.220890·10¹⁹ GeV         (PDG)
M_GUT      ≈ 2·10¹⁶ GeV                 (gauge unification)
g̃         = 0.9803                      (F5 EFT entropy scaling)
α₀         = 1069833/264500 ≈ 4.0447   (F4 algebraic anchor)
κ_TGP      = 2.012                      (XS1 cross-sector)
```

## Strategy

**Joint M_TGP–ξ_grav lock**: hypothesis `M_TGP² = g*·M_Pl²/ξ_grav`
provides 1 equation for 2 unknowns. Phase 2 picks ξ_grav from structural
inheritance (X2.2 winner) → uniquely fixes M_TGP via M_Pl PDG anchor.

**Cross-check loop**: predicted M_TGP must satisfy (a) within partial-lock
band [10¹⁶, 10¹⁹] GeV, (b) consistency z gauge-unification GUT scale,
(c) entropy scaling g̃ self-consistency.

**Falsifier**: jeśli żaden ξ_grav z 4 candidates nie daje M_Pl = PDG within
0.1% drift, χ.1 hypothesis form falsified (need additional structural factor).

## Cross-references

- [[program.md]]
- [[Phase1_results.md]]
