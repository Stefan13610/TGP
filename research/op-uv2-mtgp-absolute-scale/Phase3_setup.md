---
title: "UV.2.Phase3 setup — predictions + 4-channel convergence"
date: 2026-05-01
cycle: UV.2.Phase3
status: SETUP
parent: "[[program.md]]"
predecessor: "[[Phase2_results.md]]"
tags:
  - TGP
  - UV2
  - phase3
  - predictions
  - convergence
  - setup
---

# UV.2.Phase3 setup

**Score gate:** ≥5/6 PASS = UV.2 program END (FULL CONVERGENCE).

## Sub-tests (6)

- **U3.1** **M_GUT structural prediction**: M_GUT_UV.2 = M_TGP_χ.1/(N_A·2π²) ≈ 2.005·10¹⁶ GeV.
  - Vs SM 2-loop gauge unification central 2.0·10¹⁶ GeV: drift 0.30%.
  - Vs SM theoretical band [1.0, 2.5]·10¹⁶ GeV: well within band.
  - **Drift gate:** < 5% vs SM RG central value.

- **U3.2** **M_Pl reproduction post-UV.2 (no PDG anchor)**:
  M_Pl_UV.2 = M_GUT · 2π²·N_A^(3/2)/√g* — z M_GUT independent.
  - **Drift gate:** < 1% vs PDG.

- **U3.3** **G_N(SI) reproduction post-UV.2 (no PDG)**:
  G_N = g*/(M_TGP²·N_A) z M_TGP_UV.2 = K·M_GUT.
  - **Drift gate:** < 1% vs CODATA 2022.

- **U3.4** **f_a axion decay constant prediction (post-UV.2 → ω.3 enabling)**:
  f_a ~ M_TGP_UV.2 / (some structural factor) — band derivation.
  - QCD axion classical band: f_a ∈ [10⁹, 10¹²] GeV (Pecci-Quinn).
  - TGP ω.1 sets g_axion/f_a ~ α_em·E_TGP/(2π) = 8.30·10⁻³ → f_a relation.
  - **Forward gate** for ω.3 mini-cycle.

- **U3.5** **Gauge-grav unification structural**: M_GUT structurally promoted
  od observational → STRUCTURAL anchor (TGP-side derived).
  - SM RG-running α_1=α_2=α_3 at M_GUT cross-check (2-loop residual).
  - **Falsifier**: SM 2-loop M_GUT push outside [1.5, 2.5]·10¹⁶ post-threshold-corrections falsifies UV.2 K-LOCK.

- **U3.6** **4-channel UV.2 convergence summary**:
  - (1) UV-anchor (g* = 0.71): UV.1 NGFP exact
  - (2) Photon-ring (N_A = 500/57): ξ.1 exact
  - (3) M_GUT independent: SM 2-loop gauge unification ~2·10¹⁶ GeV
  - (4) M_Pl reproduction: PDG 1.221·10¹⁹ GeV (drift < 1%)

## Inputs from Phase 2 (LOCKED)

```
K_struct LOCK     = N_A · 2π² ≈ 173.15
M_TGP_UV.2        = 3.4630·10¹⁸ GeV (z M_GUT, drift 0.30% vs χ.1)
M_Pl_UV.2         = 1.2172·10¹⁹ GeV (drift 0.30% vs PDG)
G_N_UV.2 (SI)     = 6.7145·10⁻¹¹ m³ kg⁻¹ s⁻² (drift 0.60% vs CODATA)
g* = 71/100       UV.1
N_A = 500/57      ξ.1
M_GUT (input)     = 2.0·10¹⁶ GeV (SM 2-loop)
```

## Strategy

**True independent prediction**: M_GUT (gauge-unification scale) jest measured
zewnętrznie z SM RG running α_1=α_2=α_3, niezwiązane z grav physics. UV.2
przewiduje M_TGP, M_Pl, G_N z M_GUT + dim-less alone — to JEST niezależna
predykcja (nie tautological consistency check).

**Forward-gate**: f_a derivation post-UV.2 enables ω.3 axion decay constant
mini-cycle.

## Cross-references

- [[program.md]]
- [[Phase1_results.md]]
- [[Phase2_results.md]]
