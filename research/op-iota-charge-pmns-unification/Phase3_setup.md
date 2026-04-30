---
title: "ι.1.Phase3 setup — predictions + falsification convergence"
date: 2026-04-30
cycle: ι.1.Phase3
status: PRE-EXECUTION
parent: "[[program.md]]"
predecessor: "[[Phase2_results.md]]"
tags:
  - TGP
  - iota1
  - phase3
  - predictions
  - falsification
---

# ι.1.Phase3 — predictions + falsification convergence (6 sub-tests)

> **Cel:** Lock ι.1 PMNS mixing-operator predictions z post-Phase 2 framework;
> register II1-II6 falsification entries; close ι.1 program z full cross-sector
> lepton-quark unification verdict.

## 6 sub-tests

### I3.1 (II1) — JUNO 2027+ sin²θ₁₃ ultra-sharp window

**Prediction:** sin²θ₁₃ = K_ν · λ_C² = (1/2)·(0.22550)² = **0.0254**
- TGP framework window post-ι.1: **0.024–0.027** (zeroth-order ±10%)
- NuFit 5.3 current: 0.022 (drift 15.57%)
- JUNO 2027+ projected sensitivity: ≤ 1% precision on sin²θ₁₃

**Falsification gate:** JUNO 2027+ central value outside TGP [0.024, 0.027]
window → ι.1 mixing-operator framework FALSIFIED dla θ₁₃.

### I3.2 (II2) — DUNE/T2HK 2030+ sin²θ₂₃ hardened window

**Prediction:** sin²θ₂₃ = K_ν = **1/2 = 0.500** (maximal mixing)
- TGP framework window post-ι.1: **0.48–0.52** (Majorana-Dirac chirality lock)
- NuFit 5.3 current: 0.572 (2nd octant), drift 12.59%
- DUNE/T2HK 2030+ projected sensitivity: ≤ 2% on sin²θ₂₃ + octant resolution

**Falsification gate:** DUNE central value outside TGP [0.48, 0.52] window OR
strong octant preference (1st vs 2nd) z δ > 5σ → ι.1 K_ν chirality lock
FALSIFIED.

### I3.3 (II3) — PMNS matrix 4 free → 3 DERIVED + 1 open (δ_CP) closure

**Statement:** Post-ι.1, PMNS matrix params:
- sin²θ₁₃ = (1/2)·λ_C² **DERIVED** (mixing-operator (ν,up))
- sin²θ₂₃ = 1/2 **DERIVED** (Majorana-Dirac chirality lock)
- sin²θ₁₂ = 1/3 **DERIVED** (S₃ ⊂ GL(3,𝔽₂))
- δ_CP **OPEN** (μ.1/ν.1 cycle target, JUNO/DUNE 2030+)

**Falsification gate:** ≥1 z 3 angles inconsistent z TGP framework window
post-experimental closure 2030+ → PMNS derivation status PARTIAL/FAIL.

### I3.4 (II4) — Cross-sector lepton-quark unification full closure

**Statement:** Combined CKM (κ.1) + PMNS (ι.1) post-derivation:
- CKM: 4 free → 0 free (κ.1 closure)
- PMNS: 4 free → 3 DERIVED + 1 open (ι.1 closure)
- **Combined: 8 free → 1 open (δ_CP)** = 7/8 fundamental mixing params DERIVED

**Anchor:** cross-sector λ_C couples both sektors via Cabibbo lock z ζ.1
GL(3,𝔽₂) form factor 165/167; 4-sector chirality-counting B²-cross-product
unifies CKM Wolfenstein i PMNS angles w jeden mixing-operator framework.

**Falsification gate:** quantum coherence test 2027+ shows independent
CKM/PMNS structure → unification falsified.

### I3.5 (II5) — Future research-track hint (μ.1/ν.1 cycle)

**Hint:** δ_CP phase derivation + residual PMNS drift hardening (15.57% →
< 5% target) require:
- Higher-order mixing-operator corrections (analog do κ.1 K-level extension)
- Cross-sector phase coupling (CKM δ phase + PMNS δ_CP phase ↔ shared topology)
- Possible 5-sector extension (sterile neutrino slot at B²_sterile = ?)

**Status:** logged dla μ.1/ν.1 cycle (Wolfenstein + PMNS post-DUNE 2030+
hardening campaign).

### I3.6 (II6) — N-channel ι.1 falsification convergence

**Setup:** count independent falsification channels post-ι.1:
- C1: JUNO 2027+ sin²θ₁₃ window violation → ι.1 fail
- C2: DUNE 2030+ sin²θ₂₃ + octant violation → ι.1 fail
- C3: T2HK 2030+ overlap consistency check
- C4: 0νββ Majorana phase test (KamLAND-Zen, NEXT) → tests Majorana B²=1 ref
- C5: Cosmological Σm_ν tension → tests neutrino-as-reference framework
- C6: Sterile neutrino exclusion (DUNE near + STEREO/PROSPECT) → tests 4-sector
  closure (no 5th sektor)

**Convergence target:** ≥ 5 independent channels w 2027–2032 timeframe.

## Verdict gate

**5/6 PASS minimum** dla program END; 6/6 → FULL CONVERGENCE.

## Środowisko

```bash
PYTHONIOENCODING=utf-8 python -X utf8 research/op-iota-charge-pmns-unification/phase3_iota_predictions.py 2>&1 | tee research/op-iota-charge-pmns-unification/phase3_iota_predictions.txt
```

## Cross-references

- [[Phase2_results.md]], [[program.md]]
- [[../../PREDICTIONS_REGISTRY.md]] — II1-II6 entries (post-ι.1 registration)
- [[../op-kappa-mixing-numerator/Phase3_results.md]] — κ.1 closure baseline
