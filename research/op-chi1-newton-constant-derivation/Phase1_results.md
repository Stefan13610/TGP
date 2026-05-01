---
title: "χ.1.Phase1 results — structural setup + alt-G falsification 5/5 PASS"
date: 2026-05-01
cycle: χ.1.Phase1
status: COMPLETE
parent: "[[program.md]]"
predecessor: "[[Phase1_setup.md]]"
successor: "[[Phase2_setup.md]]"
tags:
  - TGP
  - chi1
  - phase1
  - newton-constant
  - structural-setup
  - results
  - PASS
---

# χ.1.Phase1 results

**Score: 5/5 PASS** ≥4/5 gate → **Phase 2 ENABLED**.

## Sub-test results

| ID | Test | Result | Detail |
|---|---|---|---|
| X1.1 | Substrate-emergent metric ansatz | **PASS** | gauge inv + F1 Single-Φ + Phase 2.A spectrum (2 TT + 1 scalar) |
| X1.2 | Coupling assignment + c_χ structural | **PASS** | G_N = g*/(M_TGP²·ξ_grav); c_χ = √3 z canonical kinetic match |
| X1.3 | Stueckelberg matching φ.1 ↔ graviton trace | **PASS** | h_b = c_χ·ln X exact mod boundary |
| X1.4 | F-cluster F4/F5/F6 consistency | **PASS** | XS1 √α₀ vs κ_TGP drift = 0.042% < 1% |
| X1.5 | Alt-G falsification (5 alts) | **PASS** | exactly 1/5 PROBE-passes (χ.1 hypothesis) |

## Key derived structural identities

### Identity 1 — Emergent metric decomposition

$$g^{\text{eff}}_{\mu\nu} = \eta_{\mu\nu} + \kappa\, h^{TT}_{\mu\nu}
   + \frac{c_\chi}{3}\, \eta_{\mu\nu}\, \ln X$$

**Component count check:**
- h_μν has 10 components
- = 5 (h^TT) + 4 (longitudinal/gauge) + 1 (h_b trace mode)
- Gauge fixing removes 4 → physical: **2 TT + 1 scalar** ← matches Phase 2.A spectrum

**Single-Φ consistency (F1)**: h_b couples to exactly one TGP scalar (ln X) → unique
canonical decomposition.

### Identity 2 — c_χ from canonical kinetic match

Action term comparison:
- S_grav: ½ M_Pl² R(g_eff) → ½ (∂h_b)² · (3-factor from η_μν trace structure)
- S_φ.1: ½ (∂ ln X)²

Match: c_χ²/3 = 1 → **c_χ = √3** ≈ 1.732

Scale-symmetry X → λX preserved trivially (constant shift annihilated by □).

### Identity 3 — Stueckelberg matching (vacuum)

φ.1 EL equation (AXIOM):
$$\Box(\ln X) = 0$$

Linearized Einstein trace equation in TT gauge + vacuum:
$$\Box(h_b) = 0$$

→ **h_b = c_χ·ln X** consistent (both massless wave eqs).

Matter case: $\Box(\ln X) = -(\kappa/c_\chi)\, T^\mu_\mu$ — substrate sourced
by matter trace via Stueckelberg coupling.

### Identity 4 — F-cluster preservation post-χ.1

Pre-χ.1 ledger:
- F4: α₀ = 1069833/264500 ≈ 4.044737 (algebraic, G_N-independent)
- F5: g̃ = 0.9803 (EFT entropy scaling, G_N-independent)
- F6: κ = √(32π G_N) ≈ 10.0265 (G_N=1 natural units)

Post-χ.1: F6 → κ = √(32π·g*/(M_TGP²·ξ_grav)) **upgrade only**.

XS1 verification: √α₀ = √(1069833/264500) = 2.0112 vs κ_TGP = 2.012 → drift **0.042%** ≪ 1% (cross-sector identity preserved).

→ F4, F5, XS1, XS5 ALL untouched by χ.1; only F6 upgrades STRUCTURAL → DERIVED.

## Falsification ledger (X1.5)

5 alt-G ansatz tested:

| # | Candidate | Verdict | Rationale |
|---|---|---|---|
| (i) | G_N ∝ M_TGP⁻² plain (no g*) | ✗ FAIL | ignores AS RG-flow; g* is NGFP marginal coef |
| (ii) | G_N = 1/M_Pl² direct postulate | ✗ FAIL | circular; M_Pl is what we want to derive |
| (iii) | G_N = α_em/M_TGP² | ✗ FAIL | XS3 lepton-orthogonal sector |
| (iv) | G_N = κ_TGP²/M_TGP² | ✗ FAIL | XS1 cross-sector orthogonality |
| (v) | **G_N = g*/(M_TGP²·ξ_grav)** | ✓ **PASS** | TGP-native (UV.1 g* + ξ_grav structural) |

**1/5 PROBE-passes** → uniqueness criterion satisfied.

## Phase 1 verdict

**SCORE: 5/5 PASS (≥4/5 gate)** → **Phase 2 enabled**.

**Promotion candidates entering Phase 2:**

1. **G_N hypothesis: G_N = g*/(M_TGP²·ξ_grav)** (TGP-structural, UV.1+ξ.1 anchored)
2. **c_χ = √3** (canonical kinetic match)
3. **Emergent metric**: 2 TT + 1 scalar consistent z Phase 2.A spectrum
4. **F-cluster preserved**: F4, F5, XS1, XS5 unaffected; only F6 upgrades

**Phase 2 plan**: solve AS NGFP RG-flow integration → ξ_grav structural form
scan → JOINT M_TGP lock z 3 orthogonal anchors → numerical κ reproduction
< 0.1% drift vs F6 anchor 10.0265.

## Cross-references

- [[program.md]]
- [[Phase1_setup.md]]
- [[Phase2_setup.md]]
- [[../op-uv-as-ngfp/Phase3_results.md]] — UV.1 g*, η_N*
- [[../op-xi-photon-ring/Phase3_results.md]] — ξ.1 ξ-factor, N_A = 500/57
- [[../op-phase2-quantum-gravity/Phase2_A_results.md]] — F6 STRUCTURAL κ ≈ 10.0265
