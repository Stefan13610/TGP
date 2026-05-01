---
title: "UV.2.Phase1 setup — structural setup + alt-K_struct falsification"
date: 2026-05-01
cycle: UV.2.Phase1
status: SETUP
parent: "[[program.md]]"
tags:
  - TGP
  - UV2
  - phase1
  - M_TGP
  - K_struct
  - setup
---

# UV.2.Phase1 setup

**Score gate:** ≥4/5 PASS = Phase 1 forward.

## Sub-tests (5)

- **U1.1** M_TGP–M_GUT separation hypothesis: substrate scale > gauge unification
  scale by structural factor; cyrkularność z M_Pl jest broken przy M_GUT independent.
  - Postulować: M_TGP = K_struct · M_GUT z K_struct ∈ (10, 1000) (sub-Planckian).
  - Verify: M_TGP/M_Pl = √(g*/N_A) ≈ 0.2845 z χ.1 → M_GUT/M_Pl = √(g*/N_A)/K_struct.
- **U1.2** K_struct candidates 4-form scan:
  - (a) **K_a = N_A · 2π²** (S³-volume × ξ.1 photon-ring inheritance)
  - (b) K_b = N_A² · √(2π) (square-N_A correction)
  - (c) K_c = (4π) · N_A · κ_TGP (XS1 cross-sector dim-less)
  - (d) K_d = α₀ · 4π² · √N_A (F4 algebraic anchor)
  - Compute K vs target K_target = M_TGP_χ.1/M_GUT ≈ 173.7.
  - Selection: drift < 1% partial-lock band.
- **U1.3** M_Pl/M_GUT cross-prediction: M_Pl/M_GUT = K_struct/√(g*/N_A)
  = 2π²·N_A^(3/2)/√g* (under K_a).
  - vs PDG M_Pl/M_GUT obs = 1.221·10¹⁹ / 2·10¹⁶ ≈ 610.5.
  - **Drift target:** < 1%.
- **U1.4** Joint-system consistency: M_TGP–M_Pl–M_GUT triple-lock
  pod {χ.1, UV.2} dual-derivation.
  - χ.1 fixes M_TGP/M_Pl = √(g*/N_A) (UV.1 + ξ.1 inheritance).
  - UV.2 fixes M_TGP/M_GUT = K_struct (UV.2 hypothesis).
  - → M_Pl/M_GUT = K_struct/√(g*/N_A) (derived).
- **U1.5** Alt-anchor ansatz falsification (5 alts):
  - (i) M_TGP = M_Pl postulate — circular (χ.1 tautology, no derivation gain)
  - (ii) M_TGP = M_GUT direct — fails (M_TGP/M_GUT = 1; observed 174)
  - (iii) M_TGP = √(M_Pl·M_GUT) geometric mean — fails (≈ 4.94·10¹⁷, drift 7×)
  - (iv) M_TGP = M_GUT/g* — fails (≈ 2.82·10¹⁶, drift 100×)
  - (v) **M_TGP = K_struct·M_GUT z K = N_A·2π² UV.2 winner** (drift 0.30%)

## Inputs from χ.1 (LOCKED)

```
g*         = 71/100 = 0.71               (UV.1 AS NGFP)
N_A        = 500/57 ≈ 8.7719            (ξ.1 photon-ring)
M_Pl PDG   = 1.220890·10¹⁹ GeV
M_TGP_χ.1  = 3.4734·10¹⁸ GeV (joint-lock z M_Pl PDG)
M_GUT (SM 2-loop) = 2.0·10¹⁶ GeV (theoretical ~10-30% unc.)
α₀         = 1069833/264500 ≈ 4.0447    (F4)
κ_TGP      = 2.012                       (XS1)
```

## Strategy

**Anchor swap**: zamienić M_Pl PDG anchor (χ.1) na M_GUT gauge-unification anchor.
M_GUT pochodzi z α_1 = α_2 = α_3 RG-running (sin²θ_W consistency) — independent
od grav physics. Cyrkularność broken: M_TGP derived z {g*, N_A, M_GUT} bez M_Pl input.

**Falsifier**: jeśli żaden K candidate z 4-scan nie daje drift < 1% vs M_TGP_χ.1,
UV.2 hypothesis form falsified (need different M-anchor lub additional structural factor).

## Cross-references

- [[program.md]]
- [[../op-chi1-newton-constant-derivation/Phase3_results.md]]
- [[../op-uv-as-ngfp/Phase3_results.md]]
- [[../op-xi-photon-ring/Phase3_results.md]]
