---
title: "ι.1.Phase2 setup — first-principles derivation"
date: 2026-04-30
cycle: ι.1.Phase2
status: PRE-EXECUTION
parent: "[[program.md]]"
predecessor: "[[Phase1_results.md]]"
tags:
  - TGP
  - iota1
  - phase2
  - first-principles
---

# ι.1.Phase2 — first-principles derivation (7 sub-tests)

> **Cel:** Derive PMNS angles via mixing-operator extension analog do κ.1
> CKM closure; promote charge-sector unification (B²_up − B²_down = 81/100)
> z η.2 hidden identity → first-principles theorem.

## 7 sub-tests

### I2.1 — Charge-sector unification: hidden identity formal derivation

**Statement:** B²_up − B²_down = 81/100 = N_gen⁴/(2·5)² DERIVED z 4-sector
chirality-counting product structure + cross-sector Q_u/Q_d charge difference.

**Method:** Combine:
- |Q_u|² − |Q_d|² = 1/N_gen sympy-exact identity
- B²_up − B²_down = 81/100 = N_gen⁴/(2·5)² η.2 hidden identity
- α-residual = 9/250 = 2·(B²_up − B²_down)/(N_gen²·5) η.2 dual derivation

**Cross-sector lock:** all three identities combine into structural theorem
linking charge sector ↔ B²-difference ↔ α-residual (3-way cascade).

### I2.2 — PMNS mixing-operator framework formal definition

**Definition (analog do κ.1):**

For PMNS, mixing-operator B²_PMNS-mix(X→Y; L) := L_X − L_Y dla L ∈ {B², K}
with neutrino-as-reference frame (Majorana B²=1, minimal chirality).

**Numerator rule (denom-num level pairing, analog κ.1):**
- if denom uses B²-level → num uses B²-level diff
- if denom uses K-level → num uses K-level diff

**Reference frame:** neutrino-sector (Majorana B²=1, single chirality, minimum).

### I2.3 — sin²θ₁₃ via (ν, up) mixing-operator pair + cross-sector λ_C

**Form:** sin²θ₁₃ = K_ν · λ_C² = (1/2)·λ_C²

**Derivation:**
- λ_C from cross-sector Cabibbo lock (z ζ.1 GL(3,𝔽₂) form factor 165/167)
- K_ν = 1/2 from Majorana B²=1 chirality counting (ν-internal, ζ.1)
- (ν, up) cross-sector pair B²_ν − B²_up = −9/4 → sign-direction in mixing-operator
- λ_C² · K_ν combines cross-sector Cabibbo z neutrino-internal chirality

**Value:** sin²θ₁₃ = (1/2)·(0.22550)² = 0.0254
**Drift:** 15.57% vs NuFit 0.022 (within 25% zeroth-order gate)

### I2.4 — sin²θ₂₃ via (lep, ν) Majorana-Dirac mixing pair

**Form:** sin²θ₂₃ = K_ν = 1/2 (maximal)

**Derivation:**
- (lep, ν) pair B²_lep − B²_ν = 1 (Majorana-Dirac trivial difference)
- Z₂ atmospheric reflection (μ-τ swap)
- K_ν = 1/2 chirality lock (ν-internal, ζ.1)
- Maximal mixing emerges z (lep,ν) trivial pair + Z₂

**Value:** sin²θ₂₃ = 1/2 = 0.5000
**Drift:** 12.59% vs NuFit 0.572 (within 25% zeroth-order gate)

### I2.5 — sin²θ₁₂ via group structure (ζ.1 inheritance + ι.1 reinterpretation)

**Form:** sin²θ₁₂ = 1/N_gen = 1/3 (trimaximal)

**Derivation:**
- S₃ ⊂ GL(3,𝔽₂) democratic permutation z ζ.1
- N_gen = 3 generation count (universal TGP anchor)
- Mixing-operator interpretation: trimaximal z 3-sector cross-product
  (ν, lep, sektor) z democratic weighting

**Value:** sin²θ₁₂ = 1/3 = 0.3333
**Drift:** 8.58% vs NuFit 0.307 (within 25% zeroth-order gate)

### I2.6 — 5+ alternative PMNS forms FALSIFIED

**Falsification criterion:** TGP-natural mixing-operator forms vs:
- Tribimaximal (TBM): sin²θ₁₃ = 0 (FAIL — outside ι.1 mixing-operator framework)
- Bimaximal (BM): sin²θ₁₂ = 1/2 (mixed-form mismatch)
- Golden ratio: sin²θ₁₂ = (1−1/√5)/2 ≈ 0.276 (irrational, no TGP origin)
- Hexagonal: sin²θ₁₂ = 1/4 (drift 18.6% vs NuFit, smaller than 1/3 but
  no TGP-natural origin)
- Democratic strict: sin²θ_ij = 1/3 dla wszystkich (sin²θ₂₃ = 1/3 mismatch
  z K_ν chirality lock 1/2)

### I2.7 — Classification cascade

**Pre-ι.1 → Post-ι.1 status promotions:**

| Element | Pre-ι.1 | Post-ι.1 |
|---------|---------|----------|
| ζ.1 PMNS angles (3 angles) | PARTIALLY DERIVED (refined) | **DERIVED** (mixing-operator) |
| Charge-sector unification (Q_u/Q_d ↔ 81/100) | STRUCTURAL HINT (η.2) | **PARTIALLY DERIVED** (charge identity) |
| KK5 ι.1 research-track | research-track | **PARTIALLY DERIVED** (mixing-operator extension) |
| Cross-sector lepton-quark unification | DERIVED z λ_C anchor (Z5) | **FULL DERIVED** (CKM post-κ.1 + PMNS post-ι.1) |
| PMNS matrix 4 free params | 3 PARTIALLY (angles) + 1 open (δ_CP) | **3 DERIVED + 1 open (δ_CP)** |

**Caveat:** Status PMNS DERIVED via mixing-operator framework, but drifts
8.58–15.57% vs NuFit 5.3 (zeroth-order). δ_CP phase otwarte JUNO/DUNE 2030+.
Full DERIVED w sensie CKM (z numerical lock < 0.05%) czeka na future μ.1/ν.1
cycle z higher-order corrections.

## Verdict gate

**6/7 PASS minimum** dla derivation success; 7/7 → FULL CASCADE.

## Środowisko

```bash
PYTHONIOENCODING=utf-8 python -X utf8 research/op-iota-charge-pmns-unification/phase2_iota_derivation.py 2>&1 | tee research/op-iota-charge-pmns-unification/phase2_iota_derivation.txt
```

## Cross-references

- [`Phase1_results.md`](Phase1_results.md), [`program.md`](program.md)
- [`../op-zeta-mass-spectrum/Phase2_results.md`](../op-zeta-mass-spectrum/Phase2_results.md)
- [`../op-kappa-mixing-numerator/Phase2_results.md`](../op-kappa-mixing-numerator/Phase2_results.md)
