---
title: "τ.1 program — closure overlap exponent 1/N_gen Coulomb-readjustment lift (cross-sector hardening)"
date: 2026-04-30
cycle: τ.1
status: PRE-PHASE1
parent: "[[../../INDEX.md]]"
predecessor: "[[../op-rho1-71Ge-cross-section/Phase3_results.md]]"
tags:
  - TGP
  - tau1
  - closure
  - overlap
  - coulomb
  - cross-sector
  - program
---

# τ.1 program — closure overlap exponent 1/N_gen Coulomb-readjustment lift

## Cel

Lift z **STRUCTURAL HINT → DERIVED** post-ρ.1 closure form:

```
f_overlap = (Z_a / Z_t)^{1/N_gen}    # N_gen = 3 (number of fermion generations)
```

ρ.1 zarejestrował tę formę jako STRUCTURAL HINT z f_overlap² = (Z_Ga/Z_Ge)^{2/3}
= (31/32)^{2/3} = 0.9791 dla ⁷¹Ga(ν_e, e⁻)⁷¹Ge electron-capture-like process.
τ.1 testuje uniwersalność tej formy across **3+ different EC/β-capture
nuclear sectors** (²H, ⁷Be, ⁵¹Cr, ³⁷Ar, ⁷¹Ga, ⁹⁸Mo) + **derives 1/N_gen
exponent** z TGP B²-cascade structure.

## Kontekst

Coulomb-readjustment proton-overlap correction enters bound-state ν-capture
cross-section jako multiplicative factor on |M_GT|²:

```
σ_corrected = σ_Bahcall · |⟨φ_p^t | φ_p^a⟩|² 
            ≈ σ_Bahcall · (Z_a / Z_t)^{α}
```

Standard atomic-physics derivations (Bahcall 1962, Haxton 2013) suggest
α ∈ [1, 2] z full quantum-mechanical overlap calculation. TGP-specific
ansatz post-ρ.1: **α = 2/N_gen = 2/3** (entered jako f_overlap²).

## Hypothesis τ.1

**1/N_gen exponent derives z TGP B²-cascade closure across N_gen=3 fermion
generations:**

```
f_overlap = (Z_a/Z_t)^{1/N_gen}
```

Substrate-action interpretation:
- N_gen=3 dla TGP (Standard Model 3 generations, locked z B²-cascade
  B²_lep · B²_ν · B²_up · B²_down primality structure)
- Geometric mean across N_gen channels: (Z_a/Z_t)^{1/3} as N_gen-th root
- Squared dla cross-section quadratic w |M|²: f_overlap² = (Z_a/Z_t)^{2/3}

**Cross-sector universal test:** apply f_overlap² do other EC-class
reactions w literaturze + measure Δ²/N_dof against published Bahcall
+ Haxton + Engfer/Bahcall ν-capture cross-sections.

## Reference data — EC/ν-capture target nuclei

| Reaction | Z_a | Z_t | (Z_a/Z_t)^{2/3} | Δ% | Status |
|---|---:|---:|---:|---:|---|
| ²H(ν_e,e⁻)pp | 1 | 1 | 1.000 | 0% | trivial |
| ⁷Be EC → ⁷Li | 4 | 3 | 1.211 | +21.1% | solar ν |
| ³⁷Ar(ν_e,e⁻)³⁷K¹ | 17 | 18 | 0.962 | −3.8% | radiochemical |
| ⁵¹Cr EC → ⁵¹V | 24 | 23 | 1.029 | +2.9% | GALLEX/SAGE source |
| **⁷¹Ga(ν_e,e⁻)⁷¹Ge** | **31** | **32** | **0.979** | **−2.1%** | **ρ.1 anchor ★** |
| ⁹⁸Mo(ν_e,e⁻)⁹⁸Tc | 42 | 43 | 0.985 | −1.5% | proposed FRIB |
| ¹³⁷Cs EC → ¹³⁷Xe | 55 | 54 | 1.012 | +1.2% | radiometric |

(¹) ³⁷Ar EC → ³⁷Cl, but for ν_e capture on ³⁷Cl→³⁷Ar Z_a=17, Z_t=18 (analog).

## Sub-test plan (5+7+6 = 18 sub-tests)

### Phase 1 — alt-power landscape + viability gate (5 sub-tests)

- P1.1 — literature audit: Bahcall 1962/1997 + Haxton 2013 + Engfer 2010
  Coulomb-overlap formulas dla EC/ν-capture
- P1.2 — alt-power scan: α ∈ {1/4, 1/3, 1/2, 2/3, 1, 4/3, 2} sympy-exact
  applied do ⁷¹Ga ρ.1 anchor
- P1.3 — best α selection on minimal-prime denom + TGP-cascade fit
- P1.4 — viability across cross-sector targets (4-7 nuclei)
- P1.5 — viability gate: 1/N_gen=1/3 retained as TGP-native ansatz

### Phase 2 — 1/N_gen derivation + sympy LOCK (7 sub-tests)

- P2.1 — substrate-action: N_gen-fold geometric mean across cascade
- P2.2 — sympy-exact form lock: f_overlap = (Z_a/Z_t)^{1/3}
- P2.3 — cross-sector predictions: ⁷Be, ³⁷Ar, ⁵¹Cr, ⁹⁸Mo ν-capture
- P2.4 — combined cross-section recompute z f_overlap² universal
- P2.5 — orthogonality check vs Coulomb correction F(Z, E_e) (separate
  factor, not double-counted)
- P2.6 — Δ²/N_dof goodness-of-fit across N≥3 nuclei
- P2.7 — promotions: f_overlap STRUCTURAL HINT → DERIVED, N_gen=3
  closure-anchor LOCKED, alt α∈{1/4, 1/2, 1} FALSIFIED

### Phase 3 — cross-sector predictions + 4-channel convergence (6 sub-tests)

- P3.1 — ⁷¹Ga(ν_e,e⁻)⁷¹Ge POST-CONFIRM ρ.1 anchor
- P3.2 — ⁷Be solar ν capture: Borexino-II precision 2030+
- P3.3 — ⁹⁸Mo(ν_e,e⁻)⁹⁸Tc FRIB / iThemba 2030+ predictions
- P3.4 — ⁵¹Cr / ³⁷Ar source-calibration cross-checks (historical reanalysis)
- P3.5 — orthogonal: pp solar ν (Z=1=Z trivial limit) — no closure shift
- P3.6 — 4-channel τ.1 falsification convergence

## PASS gates

- Phase 1: ≥4/5 PASS
- Phase 2: ≥6/7 PASS
- Phase 3: ≥5/6 PASS = τ.1 program END
- 6/6 Phase 3 = FULL CONVERGENCE

## Cross-references

- [[../op-rho1-71Ge-cross-section/Phase3_results.md]] — ρ.1 ⁷¹Ga anchor
- [[../op-pi1-bb0nu-nme-isotope/Phase3_results.md]] — π.1 1/A^{1/3} closure
  family (NME isotope-cross)
- [[../../INDEX.md]]
- [[../../PREDICTIONS_REGISTRY.md]]
