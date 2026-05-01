---
title: "χ.1.Phase3 setup — predictions + 4-channel convergence (CODATA G_N + PDG M_Pl + LISA + Cavendish)"
date: 2026-05-01
cycle: χ.1.Phase3
status: SETUP
parent: "[[program.md]]"
predecessor: "[[Phase2_results.md]]"
tags:
  - TGP
  - chi1
  - phase3
  - predictions
  - convergence
  - CODATA
  - PDG
  - LISA
  - setup
---

# χ.1.Phase3 setup

**Score gate:** ≥5/6 PASS = χ.1 program END (FULL CONVERGENCE).

## Sub-tests (6)

- **X3.1 — G_N(0) prediction vs CODATA 2022**
  - CODATA 2022: G_N = 6.67430·10⁻¹¹ m³ kg⁻¹ s⁻² (uncertainty 1.5·10⁻⁵, rel.).
  - χ.1 prediction: G_N_χ.1 = g*/(M_TGP²·N_A) — converted z natural GeV⁻²
    via ℏc·G_N·c² = G_N(SI) (dimensional restoration).
  - **Drift gate:** < 5·10⁻⁵ (within experimental band).
  - **Caveat:** χ.1 currently anchored M_TGP via M_Pl PDG. **Independent test**
    after M_TGP derived from non-Pl source (UV.2 mini-cycle).

- **X3.2 — M_Pl prediction vs PDG**
  - PDG M_Pl = 1.220890·10¹⁹ GeV.
  - χ.1: M_Pl_χ.1 = G_N^(-1/2) = M_TGP·√(N_A/g*) tautologically reproduces PDG
    (M_TGP defined z M_Pl). **Consistency check** (not independent prediction).
  - **Drift gate:** < 10⁻⁴ (mechanical reproducibility).

- **X3.3 — G_eff(z) cosmological evolution**
  - sek08 §6109 framework: G_eff(z) = G_N/ψ(z) z ψ(z) = X(z)/X_0.
  - φ.1 EL eq w FRW: ψ(z) bounded soft.
  - **Range:** ψ(z=2) ∈ [0.95, 1.05] → ΔG_eff/G_eff < 5%.
  - **Falsifier:** DESI DR3 2027+ + LSST 2030+ f σ_8(z) growth-rate consistency.

- **X3.4 — LISA EMRI G-running test**
  - η_N* = −2 NGFP marginal IR signature → ξ-factor running < 0.5% across
    LISA chirp band 0.1–100 mHz.
  - **Falsifier:** observed running > 0.5% → falsifies η_N* = −2 → χ.1 ξ_grav form.
  - **Forward gate** (LIVE 2035+).

- **X3.5 — Lab Cavendish-type G_N precision**
  - F1 Single-Φ → composition-independent G_N (Equivalence Principle preserved).
  - **Prediction:** G_N drift < 10⁻⁶ across labs (BIPM 2030+ projected precision).
  - **Falsifier:** composition-dependent G_N > 10⁻⁶ → falsifies F1+χ.1.
  - **Forward gate** (LIVE 2030+).

- **X3.6 — 4-channel χ.1 convergence summary**
  - (1) UV-running anchor: g* = 0.71 (UV.1 NGFP)
  - (2) F6 κ reproduction: κ = √(32π) ≈ 10.0265 (X2.4 drift 0.0001%)
  - (3) F-cluster consistency: F4 ↔ F5 ↔ XS1 (drift < 0.5%)
  - (4) Observational: CODATA G_N + PDG M_Pl (X3.1 + X3.2)
  - **Convergence gate:** all 4 channels within own targets.

## Inputs from Phase 2 (LOCKED)

```
G_N_χ.1   = 6.7088·10⁻³⁹ GeV⁻² (sympy-exact: g*/(M_TGP²·N_A))
M_TGP_χ.1 = 3.4734·10¹⁸ GeV
M_Pl_χ.1  = 1.2209·10¹⁹ GeV (matches PDG by construction)
κ_χ.1     = √(32π) = 10.026513 (matches F6 0.0001%)
g*        = 0.71  (UV.1 AS NGFP)
N_A       = 500/57 (ξ.1 photon-ring)
```

## CODATA / PDG benchmarks (2022/2024)

```
G_N (CODATA 2022, SI)    = 6.67430·10⁻¹¹ m³ kg⁻¹ s⁻² ± 1.5·10⁻⁵ rel.
M_Pl (PDG 2024)          = 1.220890·10¹⁹ GeV
ℏc (PDG)                 = 1.97327·10⁻¹⁶ GeV·m
c (exact)                = 2.99792458·10⁸ m/s
GeV→kg                   = 1.78266·10⁻²⁷ kg/GeV
```

**Conversion chain GeV⁻² → SI:**
$$G_N^{\text{SI}} = G_N^{\text{nat}} \cdot \frac{(\hbar c)^3}{(\text{GeV→kg})^2 \cdot c^4 \cdot \text{GeV}^2}$$

equivalent: $G_N^{\text{SI}} = G_N^{\text{GeV}^{-2}} \cdot (\hbar c)^3 / [\text{GeV→kg}]^2 / c^4$ — with χ.1
M_TGP = 0.2845 M_Pl this should reproduce CODATA at PDG-anchor precision.

## Strategy

**X3.1 + X3.2 are tautological** by χ.1 construction (M_TGP defined z M_Pl PDG
anchor). They serve as **dimensional-conversion verification** — confirming
that natural-units G_N → SI conversion path doesn't have arithmetic errors.

**True independent predictions** = X3.3 (cosmological), X3.4 (LISA), X3.5 (lab).
These are LIVE forward gates — chi.1 program END is conditional on
**Phase 2 KEYSTONE-confirmed** structural form holding through 2027–2035+.

**Phase 3 dual-mode**: tests X3.1 + X3.2 verify the structural form numerically
zamykając Phase 2 KEYSTONE; X3.3–X3.6 plant forward-gate flags.

## Cross-references

- [[program.md]]
- [[Phase1_results.md]]
- [[Phase2_results.md]]
