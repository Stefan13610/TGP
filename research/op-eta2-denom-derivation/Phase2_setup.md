---
title: "η.2.Phase2 setup — Wolfenstein denom + residual derivation (7 sub-tests)"
date: 2026-04-30
cycle: η.2.Phase2
status: PRE-EXECUTION
parent: "[[program.md]]"
predecessor: "[[Phase1_results.md]]"
tags:
  - TGP
  - eta2
  - denom-derivation
  - residual-cascade
  - chirality-counting
  - falsification
---

# η.2.Phase2 — Wolfenstein denom-derivation + α-residual cascade

> **Cel:** Rigorous first-principles derivation of (81, 78, 14) Wolfenstein
> denoms z 4-sector B²-cross-product + residual 0.036 derivation z B²-cascade
> form lub N_gen-direct form. Falsify 5 alternative cross-product hypotheses.

## 7 sub-tests

### B2.1 — A_TGP denom 81 = N_gen⁴ uniqueness check

**Test:** Verify 81 = 3⁴ uniquely among small-denom B²-cross-products:
- N_gen⁴ = 81 ✓
- (B²_lepton)⁵ = 32, (B²_lepton + 1)⁴ = 81 ✓ (both match!)
- N_gen³ · K_lepton_denom = 27·3 = 81 ✓ (also!)
- ...

**Falsification gate:** if ≥ 3 distinct simple forms give 81, derivation
NOT unique → STRUCTURAL HINT only. If 1-2 forms match AND share TGP-natural
structure → DERIVED candidate.

### B2.2 — ρ̄_TGP denom 78 = 2·N_gen·B²_up_num uniqueness

**Test:** 78 = 2·3·13 vs alternative cross-products:
- 2·N_gen·B²_up_num = 2·3·13 ✓ (B²_up = 13/4 numerator)
- 78 = K_lepton_denom · 26 = 3·26 (no B² connection)
- 78 = K_ν_denom · 39 (39 = 3·13)
- ...

**Falsification gate:** unique if only 1 form combines TGP-primes (2, 3, 5, 7)
+ B² values do exact 78.

### B2.3 — η̄_TGP denom 14 = K_up_num·K_lepton_num uniqueness

**Test:** 14 = 2·7 vs alternatives:
- K_up_num · K_lepton_num = 7·2 = 14 ✓
- K_ν_denom · 7 = 14 ✓
- (B²_lepton + B²_ν) · K_up_num = 3·7 = 21 (no)
- ...

**Falsification gate:** uniqueness via prime-7 + prime-2 cross-sector double-link
(K_up to η̄ + K_lepton/K_ν to η̄).

### B2.4 — Residual 0.036 derivation z B²-cascade

**Test:** Verify residual = α⁻¹(0) − 137 = 9/250 emerges z:

Form A: residual = N_gen²/(2·5³) = 9/250 [direct, ratio drift 0.0025%]
Form B: residual = 2·(B²_up − B²_down)/(N_gen²·5) [B²-cascade form]
        = 2·(81/100)/(9·5) = 9/250 sympy-exact

**Identity check:** Form A ≡ Form B sympy?
9/250 = 9/250 (yes, both reduce do same rational).

**Drift check:** vs measured residual 0.0359990840:
9/250 − 0.0359990840 = 0.0000009160, drift = 9.16·10⁻⁷ / 0.036 ≈ 2.5·10⁻³%
= 0.0025% (well within structural tolerance 0.05%).

### B2.5 — 5 alternative residual derivations FALSIFY

**Test:** Falsify alternative cross-products at 0.5% threshold:
- C1: residual ?= 1/(N_gen²·something) clean
- C2: residual ?= λ_C² · k (Cabibbo-cascade form)
- C3: residual ?= K_up · K_down · (something)
- C4: residual ?= ε_ph² scaling
- C5: residual ?= η̄·ρ̄·A scaling

**Falsification gate:** ≥ 4/5 z C1-C5 drift > 0.5% → 9/250 form unique.

### B2.6 — Cross-sector denom-prime cascade SYMPY check

**Test:** Construct unified Wolfenstein triple form z B²-cascade:
```
A    = (B²_up_num − ?)·... / N_gen⁴
ρ̄   = ? / (2·N_gen·B²_up_num)
η̄   = ? / (K_up_num · K_lepton_num)
```

**Goal:** sympy-derive numerators (64, 11, 5) z some 4-sector cross-product
combinations OR confirm że tylko denoms są derivable, numerators stay
empirical (more honest classification).

### B2.7 — Classification cascade ACTIVATION test

**Test:** Based on B2.1-B2.6 outcomes, declare:
- IF 3/3 denom uniqueness PASS + residual derivation PASS + 5/5 alternatives
  FALSIFY → **η.1 promotion PARTIALLY DERIVED → DERIVED** (Wolfenstein triple
  denoms structurally locked); **α.1 residual STRUCTURAL HINT → PARTIALLY
  DERIVED** (9/250 = N_gen²/(2·5³) sympy-exact).
- IF 2/3 denoms unique + residual derivation PASS → η.1 stays PARTIALLY DERIVED
  (refined²); α.1 residual upgrades do PARTIALLY DERIVED.
- IF ≤ 1/3 unique OR residual derivation FAILS → no upgrade, both stay current
  classification (research-track honest).

## Verdict gate

- **7/7 PASS** → Phase 3 predictions + η.2 program END declaration.
- **6/7 PASS** → Phase 3 z minor caveat.
- **≤ 5/7 PASS** → η.2.Phase2 reframing required.

## Środowisko

```bash
PYTHONIOENCODING=utf-8 python -X utf8 research/op-eta2-denom-derivation/phase2_denom_derivation.py 2>&1 | tee research/op-eta2-denom-derivation/phase2_denom_derivation.txt
```

## Cross-references

- [`Phase1_results.md`](Phase1_results.md), [`program.md`](program.md)
- [`../op-eta-wolfenstein/Phase3_results.md`](../op-eta-wolfenstein/Phase3_results.md)
- [`../op-theta-quark-koide/Phase3_results.md`](../op-theta-quark-koide/Phase3_results.md)
- [`../op-alpha-fine-structure/Phase3_results.md`](../op-alpha-fine-structure/Phase3_results.md)
