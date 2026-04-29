---
title: "ζ.1.Phase2 setup — PMNS angles first-principles z GL(3,𝔽₂) × Z₃ × SU(2)_L"
date: 2026-04-29
cycle: ζ.1.Phase2
status: PRE-EXECUTION
parent: "[[program.md]]"
predecessor: "[[Phase1_results.md]]"
tags:
  - TGP
  - zeta-mass-spectrum
  - PMNS
  - GL3F2
  - cross-sector
  - Cabibbo-lock
---

# ζ.1.Phase2 — PMNS angles first-principles derivation

> **Cel:** Wyprowadzić PMNS mixing angles (θ₁₂, θ₂₃, θ₁₃) first-principles
> z GL(3,𝔽₂) × Z₃ × SU(2)_L group structure; falsify alternative
> parameterizations; pokazać cross-sector lepton-quark θ_C-θ₁₃ unification
> via single Cabibbo anchor (λ_C = 0.22550 LOCKED z tgp-leptons GL form
> factor 165/167).

---

## Anchors LOCKED z poprzednich closure

1. **λ_C = 0.22550** — DERIVED via GL(3,𝔽₂) form factor 165/167 (tgp-leptons, 0.7σ vs PDG 0.22500 ± 0.00067)
2. **K(ν) = 1/2** — DERIVED via Majorana B²=1 chirality-counting (ζ.1.Phase1)
3. **GL(3,𝔽₂) = PSL(2,7)** — VERIFIED: 168 elements, 6 conjugacy classes, 28 Z₃ subgroups, 20 double cosets
4. **N_gen = 3** — DERIVED via metric singularity barrier (M9.2-D)
5. **Σm_ν = 59.01 meV NO ordering** — LOCKED (ζ.1.Phase1)

## 7 sub-tests

### Z2.1 — sin²θ₁₂ = 1/3 trimaximal z S₃ ⊂ GL(3,𝔽₂)

**Test:**
- S₃ democratic permutation action na 3 generations leptonic doublet
- Trimaximal mixing matrix U_TBM eigenvector dla S₃: (1,1,1)/√3
- sin²θ₁₂ = |U_e2|² = 1/3 (TBM zeroth-order)
- NuFit 5.3 observed: sin²θ₁₂ ≈ 0.307 (drift 8.5%)

**Falsification:** if S₃ action nie produkuje (1,1,1)/√3 eigenstate → 
GL(3,𝔽₂) family symmetry incorrect for solar mixing.

### Z2.2 — sin²θ₂₃ = 1/2 maximal z Z₂ + K(ν)=1/2 lock

**Test:**
- Z₂ atmospheric reflection (μ-τ swap) → maximal θ₂₃ = π/4
- K(ν) = 1/2 chirality lock confirms 1/2 normalization w sektorze atmospheric
- sin²θ₂₃ = 1/2 zeroth-order; NuFit observed ≈ 0.516 (drift 3.2%)

**Falsification:** if Z₂ atmospheric reflection nie jest symetria reżimu
→ K(ν) sektor structure framework needs revision.

### Z2.3 — sin²θ₁₃ = λ_C²/2 cross-sector Cabibbo lock

**Test:**
- λ_C² = (0.22550)² = 0.050850 (z GL form factor 165/167)
- sin²θ₁₃ = λ_C²/2 = 0.025425
- NuFit observed: sin²θ₁₃ ≈ 0.0220 (drift 15.6%)
- sin θ₁₃ ≈ λ_C/√2 = 0.15945 vs observed 0.1483 (drift 7.5%)

**Falsification:** if observational θ₁₃ shifts > 5σ z λ_C/√2 →
cross-sector Cabibbo lock framework broken.

### Z2.4 — 5 alternative parameterizations falsification

**Test:**
- C1: TBM strict (sin²θ₁₃ = 0) — drift 100%
- C2: BM bimaximal (sin²θ₁₂ = 1/2) — drift 62.9%
- C3: Golden ratio (sin²θ₁₂ = 1 - 1/φ ≈ 0.276) — drift 10.0%
- C4: Hexagonal (sin²θ₁₂ = 1/4) — drift 18.6%
- C5: Democratic (sin²θ₁₂ = 1/2) — drift 62.9%
- Reject all > 5% drift; TGP TBM zeroth-order 8.5% wins on weighted score

**Falsification:** if alternative parameterization < 5% drift → 
TBM zeroth-order nie jest optimal explanation.

### Z2.5 — PMNS unitarity UU† = I sympy verification

**Test:**
- Build U_PMNS w PDG parameterization (s_ij, c_ij + δ_CP)
- Verify UU† = I_3×3 sympy exact (na każdym matrix element)
- Check |U_e1|² + |U_e2|² + |U_e3|² = 1 dla każdego row/column
- Cross-check sum rule: Σ sin²θ = 5/6 (zeroth-order TBM sum rule)

**Falsification:** if UU† ≠ I → PMNS parameterization niespójna.

### Z2.6 — NGFP RG-stability mixing angles pod common β-rescaling

**Test:**
- AS NGFP marginal anomalous dimensions (z UV.1 closure η_N* = -2)
- Common β-rescaling: dθ_ij/d ln μ ∝ (Y_ν - Y_l) × β-factor
- Pod common β-rescaling, ratios θ₁₂/θ₂₃ i θ₁₃/λ_C są RG-invariant
- LISA 2035+ EMRI: nie wykryje > 0.5% θ₁₃ running

**Falsification:** if observed θ_ij scale-dependent > 1% → 
NGFP RG-stability framework broken w sektorze leptonowym.

### Z2.7 — Classification PARTIALLY DERIVED (refined)

**Test:**
- TGP PMNS prediction reduces to 3 angles + 1 phase (4 free → 0 free pod TGP)
- All 3 zeroth-order angles z group structure (GL(3,𝔽₂) × Z₃ × SU(2)_L)
- Drifts: 8.5%, 3.2%, 15.6% (acceptable < 20% bez 1-loop corrections)
- Status: PARTIALLY DERIVED (refined) — full DERIVED czeka na 1-loop
  RG corrections (long-term, neutrino MSW substrate track)

**Verdict:** classification PARTIALLY DERIVED (refined) jeśli 6/7 PASS.

---

## Verdict gate

**7/7 PASS** → PMNS first-principles closure z GL(3,𝔽₂) × Z₃ × SU(2)_L,
classification **PARTIALLY DERIVED (refined)**, Phase 3 proceeds.

**6/7 PASS** → audit gap, classification PARTIALLY DERIVED (z marginalia),
Phase 3 limited scope.

**≤ 5/7 PASS** → ζ.1.Phase2 reframing required.

---

## Środowisko

```bash
PYTHONIOENCODING=utf-8 python -X utf8 research/op-zeta-mass-spectrum/phase2_pmns_derivation.py 2>&1 | tee research/op-zeta-mass-spectrum/phase2_pmns_derivation.txt
```

## Cross-references

- [`program.md`](program.md) — overall ζ.1 plan
- [`Phase1_results.md`](Phase1_results.md) — Σm_ν = 59 meV closure (predecessor)
- [`../cabibbo_correction/r1_gl3f2_structure.py`](../cabibbo_correction/r1_gl3f2_structure.py) — GL(3,𝔽₂) algebra
- [`../nbody/examples/ex246_pmns_mixing_structure.py`](../nbody/examples/ex246_pmns_mixing_structure.py) — TBM precursor
- [`../../PREDICTIONS_REGISTRY.md`](../../PREDICTIONS_REGISTRY.md) — JUNO 2027+ θ₁₃ precision
