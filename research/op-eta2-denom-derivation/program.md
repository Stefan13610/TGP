---
title: "η.2 program — Wolfenstein denom-derivation + α-residual cascade extension"
date: 2026-04-30
cycle: η.2
status: ACTIVE
parent: "[[../../INDEX.md]]"
predecessor: "[[../op-eta-wolfenstein/Phase3_results.md]]"
sibling_anchors:
  - "[[../op-theta-quark-koide/Phase3_results.md]]"
  - "[[../op-zeta-mass-spectrum/Phase3_results.md]]"
  - "[[../op-alpha-fine-structure/Phase3_results.md]]"
tags:
  - TGP
  - eta2
  - wolfenstein-denom
  - chirality-counting
  - alpha-residual
  - cross-sector-cascade
---

# η.2 — Wolfenstein denom-derivation + α-residual cascade extension

> **Cel:** Rigorous derivation of Wolfenstein triple denominators (81, 78, 14)
> z 4-sector chirality-counting B² cross-product (lepton/neutrino/up/down) +
> check if α-residual 0.036 emerges z 4-sector cross-product OR F4-extension.
>
> **Inheritance:**
> - η.1 (Phase3 6/6): Wolfenstein triple (A, ρ̄, η̄) = (64/81, 11/78, 5/14) LOCKED, classification PARTIALLY DERIVED (refined)
> - θ.1 (Phase3 6/6): K-taxonomy 4-sector, K = (2+B²)/(2N) for N=3, B²-values (lepton 2, ν 1, up 13/4, down 61/25)
> - ζ.1 (Phase3 6/6): Single Cabibbo anchor λ_C = 0.22550 governs CKM + PMNS
> - α.1 (Phase3 6/6): 137 sympy-LOCKED via F4 chain ψ_ph = 4/(3+17/40) = 160/137; residual 0.036 STRUCTURAL HINT
>
> **Hypothesis to test:** denoms (81, 78, 14) i residual 0.036 są derivable z
> 4-sector B² cross-product struktury (chirality-counting universal pattern).
> If TRUE → η.1 Wolfenstein triple promoted PARTIALLY DERIVED → **DERIVED**;
> α.1 residual promoted STRUCTURAL HINT → **PARTIALLY DERIVED**.
> If FALSE (denom-pattern lacks cross-sector derivation) → η.1 / α.1
> classifications stay PARTIALLY DERIVED + STRUCTURAL HINT (research-track open).

## Master ledger trajectory

- Pre-η.2: **463** (post-α.1)
- Post-η.2 (target): **481** (5 + 7 + 6 = +18)

## 3-phase plan

| Phase | Sub-tests | Goal | Verdict gate |
|---|---|---|---|
| **Phase 1** | 5 | Cross-sector denom landscape audit (81, 78, 14, 8, 50, 3, 2, 137, 250 inventory) + B²-cross-product candidate ranking + prime factor inheritance graph | **5/5 PASS** → Phase 2 viable |
| **Phase 2** | 7 | First-principles derivation attempt (81, 78, 14) z 4-sector B² cross-product + residual 0.036 emergence test + 5 alternative cross-product hypotheses falsification | **7/7 PASS** → η.2 promotion cascade ACTIVATED, **6/7 PASS** → minor caveat, **<6 PASS** → reframing required |
| **Phase 3** | 6 | 6 predictions HH1-HH6: sharper Wolfenstein bounds, residual 0.036 derivation check (η.2-refined), cross-sector denom-pattern uniqueness test, 4-sector B²-cascade verification, lepton-quark-α prime-cascade unification, 5-channel η.2 falsification convergence | **6/6 PASS** → **η.2 program END** |

**Total:** 5 + 7 + 6 = **18 sub-tests**.

## Why now?

After η.1 + α.1 closures, two cross-sector hints remain unresolved:

1. **η.1.H4 (STRUCTURAL hint):** denom-prime sharing prime-3 (A=81, ρ̄=78, K_lepton denom) + prime-7 (η̄=14, K_up=7/8); derivation OPEN.

2. **α.1.A5 (STRUCTURAL HINT):** residual 0.036 = α⁻¹(0) − 137; best fit 9/250 drift 0.0025% but soft denom 250=2·5³ lacks cross-sector TGP meaning; research-track η.2/β.1 future cycle.

η.2 is the natural next cycle to **stress-test** both hints simultaneously
via 4-sector B²-cascade unification framework.

## Hypothesis structure

### H_η2_A: 4-sector B²-cross-product → Wolfenstein denoms

```
B²_lepton = 2     (Dirac, 2 chiralities)
B²_ν      = 1     (Majorana, 1 chirality)
B²_up     = 13/4  (Dirac + QCD/color asymmetry)
B²_down   = 61/25 (Dirac + QCD effective)

Hypothesis: denoms (81, 78, 14) = polynomial(B²_i) for some
            chirality-counting cross-product OR sum/LCM/product structure.

Examples to test:
  - 81 = 3⁴ → 4 sectors × 3 generations? Or (B²_lepton)²·(N_gen)²·...?
  - 78 = 2·3·13 → 13 from B²_up = 13/4 numerator?
  - 14 = 2·7 → prime-7 unique do K_up = 7/8 numerator
```

### H_η2_B: residual 0.036 z 4-sector cross-product

```
α⁻¹(0) − 137 = 0.0359990...
9/250 fit = 0.036 (drift 0.0025%) — soft denom

Hypothesis: residual = f(B²_lepton, B²_ν, B²_up, B²_down) / G
            where G is some cross-sector geometric form factor.

Examples to test:
  - residual = (B²_lepton − B²_ν)/(...) z F4 chain?
  - residual = (B²_up − B²_down)/(some power-of-N)?
  - residual ?= cross-sector η̄ + ρ̄ correction term?
```

## Falsification gate

- **Phase 2 5/5 alternatives FALSIFIED** (drift > 1% across cross-products) → 
  TGP-best derivation < 0.05% drift uniquely promoted DERIVED.
- **Phase 2 ≥3/5 alternatives PASS** (<1% drift) → derivation NOT unique → 
  STRUCTURAL HINT only (research-track stays open).

## Cross-references

- [`../op-eta-wolfenstein/Phase3_results.md`](../op-eta-wolfenstein/Phase3_results.md) — η.1 triple LOCKED + H4 hint
- [`../op-theta-quark-koide/Phase3_results.md`](../op-theta-quark-koide/Phase3_results.md) — K-taxonomy 4-sector
- [`../op-alpha-fine-structure/Phase3_results.md`](../op-alpha-fine-structure/Phase3_results.md) — α.1 + residual A5 hint
- [`../../PREDICTIONS_REGISTRY.md`](../../PREDICTIONS_REGISTRY.md) — H4, A5 entries
- [`../../INDEX.md`](../../INDEX.md) — master ledger 463 → 481 target

## Decyzja

η.2 ma być honest stress-test 2 cross-sector hints (H4 + A5) jednocześnie.
Outcome może być pozytywny (cascade DERIVED z 4-sector B²) lub negatywny
(denom-pattern strukturalnie incompatible z B²-cascade) — oba scenariusze
są walidne i będą logged transparently.
