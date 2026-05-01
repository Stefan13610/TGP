---
title: "χ.1.Phase1 setup — structural setup + alt-G ansatz falsification"
date: 2026-05-01
cycle: χ.1.Phase1
status: SETUP
parent: "[[program.md]]"
tags:
  - TGP
  - chi1
  - phase1
  - newton-constant
  - emergent-gravity
  - structural-setup
  - setup
---

# χ.1.Phase1 setup

**Score gate:** ≥4/5 PASS = Phase 1 forward.

## Sub-tests (5)

- **X1.1** Substrate-emergent metric ansatz: derive `g_eff_μν = η_μν + κ·h_μν^TT
  + (c_χ/3)·η_μν·ln X` jako unique kanoniczna dekompozycja preserwująca:
  - (a) gauge inv `δh_μν = ∂_μξ_ν + ∂_νξ_μ`
  - (b) Single-Φ axiom F1 (1 trace mode)
  - (c) Phase 2.A spectrum (2 TT + 1 scalar)
- **X1.2** Coupling assignment: postulate `G_N = g*/(M_TGP²·ξ_grav)`;
  derive `c_χ = 1` z scale-symmetry X→λX (analog Noether current φ.1)
- **X1.3** Stueckelberg matching: φ.1 EL eq `□(ln X) = 0` ↔ linear graviton
  trace eq `□h_b = κ·T^μ_μ`; consistency vacuum (T^μ_μ = 0) → `h_b = c_χ·ln X`
  exact mod boundary
- **X1.4** F-cluster ledger consistency: F4 (α₀ = 1069833/264500) ↔ F5 (g̃ =
  0.9803) ↔ F6 (κ = 10.0265) self-consistency pod χ.1 framework — verify
  post-χ.1 reinterpretacja F6 NIE narusza F4/F5/XS5
- **X1.5** **Alt-G ansatz falsification (5 alts):**
  - (i) `G_N ∝ M_TGP^(-2)` plain (no g*) — fail: ignores AS RG-flow
  - (ii) `G_N = 1/M_Pl²` direct (M_Pl postulat) — circular; no derywacja
  - (iii) `G_N = α_em / M_TGP²` — fail z XS3 lepton-orthogonal
  - (iv) `G_N = κ_TGP² / M_TGP²` — fail z XS1/XS3 cross-sector orthogonality
  - (v) `G_N = g*/(M_TGP²·ξ_grav)` — TEST candidate (χ.1 hypothesis)

## Inputs from prior cycles (LOCKED)

```
g*        = 0.71              (UV.1 AS NGFP)
λ*        = 0.19              (UV.1 AS NGFP)
η_N*      = -2                (UV.1 AS NGFP marginal limit)
N_A       = 500/57            (ξ.1 photon-ring sympy-exact)
F4 α₀     = 1069833/264500    (≈ 4.04472)
F5 g̃     = 0.9803             (Phase 2.E.3 entropy scaling)
F6 κ      = √(32πG_N) ≈ 10.0265 (G_N=1 natural units; STRUCTURAL)
M_TGP     ∈ [10¹⁶, 10¹⁹] GeV  (joint-lock partial; central ~10¹⁶ GUT scale)
G_N       = 6.67430e-11 m³ kg⁻¹ s⁻² (CODATA 2022)
M_Pl      = 1.220890e19 GeV   (PDG)
```

## Cross-references

- [[program.md]]
- [[../op-phi1-substrate-action-variational/Phase3_results.md]] — φ.1 AXIOM
- [[../op-uv-as-ngfp/Phase3_results.md]] — UV.1 g*, η_N*
- [[../op-xi-photon-ring/Phase3_results.md]] — ξ.1 ξ-factor, N_A
- [[../op-phase2-quantum-gravity/Phase2_A_results.md]] — F6 STRUCTURAL
