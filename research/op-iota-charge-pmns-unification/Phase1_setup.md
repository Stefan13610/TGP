---
title: "ι.1.Phase1 setup — charge-sector + PMNS landscape audit"
date: 2026-04-30
cycle: ι.1.Phase1
status: PRE-EXECUTION
parent: "[[program.md]]"
predecessor: "[[program.md]]"
tags:
  - TGP
  - iota1
  - phase1
  - landscape-audit
---

# ι.1.Phase1 — charge-sector + PMNS landscape audit (5 sub-tests)

> **Cel:** Audit cross-sector pair B²-differences + PMNS angle anchors + hidden
> identity B²_up − B²_down = 81/100 + charge-sector ratio Q_u/Q_d connection.

## 5 sub-tests / sub-questions

### I1.1 — Hidden identity (B²_up − B²_down) = 81/100 sympy re-confirmation

**Test:** B²_up − B²_down = 13/4 − 61/25 = (325 − 244)/100 = **81/100** sympy-exact;
verify identity 81/100 = N_gen⁴/(2·5)² = 81/100 ✓.

**Cross-sector positioning:** charge-sector B²-difference unique w 4-sector landscape.

**Falsification gate:** identity drift > 0.1% sympy → η.2 hidden identity broken.

### I1.2 — Cross-sector pair B²-differences inventory

**Test:** Compute wszystkie 6 pair B²-differences w 4-sector taxonomy:
- (ν, up): B²_ν − B²_up = 1 − 13/4 = **−9/4**
- (ν, down): B²_ν − B²_down = 1 − 61/25 = **−36/25**
- (lep, up): B²_lep − B²_up = 2 − 13/4 = **−5/4** = −QCD_up
- (lep, down): B²_lep − B²_down = 2 − 61/25 = **−11/25** = −QCD_down
- (lep, ν): B²_lep − B²_ν = 2 − 1 = **1** = Majorana-Dirac trivial
- (up, down): B²_up − B²_down = **81/100** = charge-sector identity

**Falsification gate:** ≥1 z 6 pairs not sympy-rational → mixing-operator framework broken.

### I1.3 — PMNS angle sympy values + drift analysis

**Test:** Z ζ.1 results:
- sin²θ₁₂ = 1/3 trimaximal, drift 8.6% vs NuFit 5.3 0.307
- sin²θ₂₃ = 1/2 maximal, drift 12.6% vs NuFit 0.572 (2nd octant)
- sin²θ₁₃ = λ_C²/2 ≈ 0.0254 z λ_C = 0.22550, drift 15.6% vs NuFit 0.022

**Status pre-ι.1:** PARTIALLY DERIVED (refined). Window dla ι.1 promotion: drifts
może być reduced via mixing-operator interpretation z (ν,up)/(lep,ν) pair.

### I1.4 — Charge-ratio Q_u/Q_d ↔ B²-difference structural cross-check

**Test:** Q_u = +2/3, Q_d = −1/3, ratio Q_u/Q_d = −2 (sign-charge ratio).
Connect z B²-difference via:
- |Q_u|² + |Q_d|² = 4/9 + 1/9 = 5/9
- |Q_u|² − |Q_d|² = 3/9 = 1/3
- (B²_up − B²_down)/something z (|Q_u|², |Q_d|²) algebra

**Cel:** sprawdzić czy 81/100 = N_gen⁴/(2·5)² ma reinterpretation w terminach charge.

### I1.5 — Viability gate dla Phase 2

**Criteria:** 4/5 above PASS + cross-sector pair coverage ≥ 5/6 + PMNS drift gap
za sin²θ z mixing-operator pair-form ≤ current ζ.1 drift → Phase 2 viable.

**Fallback:** jeśli viability fails → ι.1 reframing, no Phase 2.

## Środowisko

```bash
PYTHONIOENCODING=utf-8 python -X utf8 research/op-iota-charge-pmns-unification/phase1_iota_landscape.py 2>&1 | tee research/op-iota-charge-pmns-unification/phase1_iota_landscape.txt
```

## Cross-references

- [`program.md`](program.md), [`../op-kappa-mixing-numerator/Phase2_results.md`](../op-kappa-mixing-numerator/Phase2_results.md)
- [`../op-zeta-mass-spectrum/Phase2_results.md`](../op-zeta-mass-spectrum/Phase2_results.md)
- [`../op-eta2-denom-derivation/Phase2_results.md`](../op-eta2-denom-derivation/Phase2_results.md)
