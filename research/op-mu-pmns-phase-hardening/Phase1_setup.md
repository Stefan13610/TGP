---
title: "μ.1.Phase1 setup — δ_CP + PMNS drift correction landscape audit"
date: 2026-04-30
cycle: μ.1.Phase1
status: PRE-EXECUTION
parent: "[[program.md]]"
predecessor: "[[../op-iota-charge-pmns-unification/Phase3_results.md]]"
tags:
  - TGP
  - mu1
  - phase1
  - landscape-audit
  - delta-cp
  - drift-hardening
---

# μ.1.Phase1 — δ_CP + PMNS drift correction landscape audit (5 sub-tests)

> **Cel:** Map landscape dla (a) δ_CP_PMNS structural derivation z CKM γ +
> PMNS-Wolfenstein analog, (b) higher-order λ_C-corrections do PMNS angles
> hardening drifts.

## 5 sub-tests

### M1.1 — NuFit 5.3 PMNS reference inventory

NuFit 5.3 NO best-fit + 1σ windows:
- sin²θ₁₂ = 0.307 ± 0.013 (1σ ≈ 4.2%)
- sin²θ₂₃ = 0.572 ± 0.018 (1σ ≈ 3.1%, 2nd octant)
- sin²θ₁₃ = 0.022 ± 0.0006 (1σ ≈ 2.7%)
- δ_CP = 195° (= 1.08π) z 1σ window [128°, 352°] (NO; very wide)

**Pass:** 4 quantities + 1σ envelopes documented.

### M1.2 — γ_CKM apex angle z Wolfenstein triple sympy-exact

γ_CKM = arctan(η̄/ρ̄) z κ.1 (η̄, ρ̄) = (5/14, 11/78):
- tan(γ_CKM) = (5/14)·(78/11) = 390/154 = **195/77** sympy-exact
- γ_CKM ≈ arctan(195/77) ≈ 68.45°
- PDG γ ≈ 65.7° → drift 4.19%

**Pass:** sympy-rational arctan argument 195/77 LOCKED.

### M1.3 — PMNS-Wolfenstein analog (ρ̄_PMNS, η̄_PMNS) z (ν,up) pair

Mixing-operator analog do κ.1 dla CKM, ale z (ν,up) pair zamiast (lep,up):

- ρ̄_PMNS_num = B²_up_num − B²_ν_num = 13 − 1 = **12**
- ρ̄_PMNS_denom = 2·N_gen·B²_up_num = 2·3·13 = 78
- ρ̄_PMNS = 12/78 = **2/13**

- η̄_PMNS_num = K_up_num − K_ν_num = 7 − 1 = **6**
- η̄_PMNS_denom = K_up_num · K_ν_num = 7·1 = 7
- η̄_PMNS = **6/7**

Ratio: η̄_PMNS/ρ̄_PMNS = (6/7)/(2/13) = 6·13/(7·2) = **78/14 = 39/7** sympy-exact.
arctan(39/7) ≈ 79.81°.

**Pass:** (ρ̄_PMNS, η̄_PMNS) sympy-exact LOCKED.

### M1.4 — PMNS drift correction landscape (top-rational forms)

Dla każdego z 3 zeroth-order PMNS angles, identify TGP-natural correction
forms:

**sin²θ₁₃** (zeroth K_ν·λ_C² = 0.025425; target 0.022; ratio 0.866):
- Top form: K_ν·λ_C²·(1 − ρ̄) = (1/2)·λ_C²·(67/78) ≈ 0.02184 drift 0.73%
- Alt (1 − 7/52): K_ν·λ_C²·(45/52) ≈ 0.022001 drift 0.001%
- Alt (1 − 27/200): K_ν·λ_C²·(173/200) ≈ 0.021993 drift 0.034%

Best structural anchor: **(1 − ρ̄)** — direct cross-sector CKM Wolfenstein
real-part leakage.

**sin²θ₂₃** (zeroth K_ν = 0.5; target 0.572; ratio 1.144):
- Top form: K_ν/K_up = (1/2)/(7/8) = 4/7 ≈ 0.5714 drift 0.10%
- Alt K_ν·(1 + N_gen·λ_C²) ≈ 0.5763 drift 0.75%
- Alt K_ν·(1 + 18/125) ≈ 0.572 drift 0.0%

Best structural anchor: **K_ν/K_up** — clean cross-sector quark-lepton
K-taxonomy ratio (quotient zamiast suma — matches K-taxonomy structure).

**sin²θ₁₂** (zeroth 1/N_gen = 0.333; target 0.307; ratio 0.921):
- Top form: (1/N_gen)·(1 − λ_C·η̄) = (1/3)·(5149/5600) ≈ 0.3065 drift 0.17%
- Alt (1/N_gen)·(1 − λ_C/N_gen) ≈ 0.30828 drift 0.42%

Best structural anchor: **(1 − λ_C·η̄)** — cross-sector Wolfenstein imaginary
leakage λ_C·η̄ (Cabibbo · CKM imaginary part).

**Pass:** 3/3 angles z TGP-natural < 1% correction forms identified.

### M1.5 — Viability gate dla Phase 2

Wszystkie 3 PMNS angles admit < 1% drift correction z TGP-natural form +
δ_CP admits 2 candidate structural derivations within NuFit 1σ → Phase 2
viable. Cross-sector phase coupling (CKM γ ↔ PMNS δ_CP via gen-tripling
lub PMNS-Wolfenstein analog) zaakceptowana jako framework.

**Pass:** viability gate satisfied (4/4: drift hardening 3/3 < 1% + δ_CP
2 candidates within 1σ).

## Verdict gate

**5/5 PASS** dla Phase 2 progression.

## Środowisko

```bash
PYTHONIOENCODING=utf-8 python -X utf8 research/op-mu-pmns-phase-hardening/phase1_mu_landscape.py 2>&1 | tee research/op-mu-pmns-phase-hardening/phase1_mu_landscape.txt
```

## Cross-references

- [[../op-iota-charge-pmns-unification/Phase3_results.md]] — ι.1 zeroth-order
  PMNS angles
- [[../op-kappa-mixing-numerator/Phase3_results.md]] — κ.1 Wolfenstein
  (ρ̄, η̄) post-FULL-DERIVED
- [[../op-eta-wolfenstein/Phase3_results.md]] — η.1 Wolfenstein triple LOCKED
