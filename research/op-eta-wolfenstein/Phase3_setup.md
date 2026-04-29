---
title: "η.1.Phase3 setup — 6 predictions H1-H6 + cross-sector cascade"
date: 2026-04-29
cycle: η.1.Phase3
status: PRE-EXECUTION
parent: "[[program.md]]"
predecessor: "[[Phase2_results.md]]"
tags:
  - TGP
  - eta-wolfenstein
  - CKM
  - predictions
  - cross-sector
  - falsification-roadmap
---

# η.1.Phase3 — 6 predictions H1-H6 + cross-sector cascade

> **Cel:** Generate 6 falsifiable predictions across CKM Wolfenstein
> sektor + cross-sector cascade; triple LOCKED (64/81, 11/78, 5/14);
> classification PARTIALLY DERIVED (refined) confirmed; η.1 program END.

---

## 6 sub-tests / predictions

### T3.1 (H1) — Belle II 2027+ |V_ub| refined window

**Prediction:**
- TGP V_ub_η.1 = (64/81) · λ_C³ · √((11/78)² + (5/14)²) = 0.003479
- PDG 2024: |V_ub| = 0.00382 ± 0.00010 (drift 8.93%)
- Belle II projected (2027+): σ(|V_ub|_excl) ~ 1.5%, σ(|V_ub|_incl) ~ 2%
- **Falsification gate:** if Belle II measures |V_ub| outside [3.40, 4.00]·10⁻³
  at >5σ, η.1 triple lock broken
- **Confirmation gate:** if |V_ub| stays w window, Wolfenstein triple
  (64/81, 11/78, 5/14) cross-sector locked
- **Cross-sector signature:** if |V_ub| measured at ~3.50·10⁻³, η̄ = 5/14
  promoted DERIVED

### T3.2 (H2) — LHCb Run 4 (2030+) Jarlskog J refined window

**Prediction:**
- TGP J_η.1 = (64/81)² · λ_C⁶ · (5/14) = 2.932·10⁻⁵
- PDG 2024: J = (3.07 ± 0.10)·10⁻⁵ (drift 4.51%)
- LHCb Run 4 projected: σ(J) ~ 1% absolute
- **Falsification gate:** if LHCb J outside [2.85, 3.30]·10⁻⁵ at >5σ,
  Wolfenstein λ⁶ cascade z η.1 triple broken
- **Confirmation gate:** if J ≈ 2.93·10⁻⁵ confirmed, A² · η̄ structural
  product LOCKED

### T3.3 (H3) — Unitarity triangle β angle prediction

**Prediction:**
- TGP β = arg(- V_cd · V_cb* / V_td · V_tb*) z (A_TGP, λ_C, ρ̄_TGP, η̄_TGP)
- PDG 2024: sin(2β) = 0.699 ± 0.017 → β ≈ 22.2° ± 0.7°
- TGP β predicted z η.1 triple
- **Falsification gate:** if Belle II + LHCb measure sin(2β) outside
  [0.65, 0.75] at >5σ, triangle β closure broken
- **Confirmation gate:** if sin(2β) ≈ 0.69-0.71, η.1 triangle apex
  geometry confirmed

### T3.4 (H4) — Cross-sector A ↔ K_up cross-link (denom-3/7 hint)

**Prediction:**
- A_TGP = 64/81 (denom 81 = 3⁴)
- η̄_TGP = 5/14 (denom 14 = 2·7)
- K_up = 7/8 (numerator 7, denom 8 = 2³)
- K_lepton = 2/3 (denom 3)
- **Cross-sector denom-prime sharing:** prime 3 between A & K_lepton;
  prime 7 between η̄ denom & K_up numerator
- **Open hypothesis:** chirality-counting B²-extension might predict
  Wolfenstein denoms via 4-sector cross-product
- **Falsification gate:** if rigorous derivation excludes denom-3/7 sharing,
  cross-sector cascade hint broken
- **Confirmation gate:** if future η.2 cycle derives 81 = 3⁴ z B²_lepton
  cascade & 14 = 2·7 z B²_up · 4/7 form, cross-sector framework promoted DERIVED

### T3.5 (H5) — V_td/V_ts cross-check refined (B-B̄ mixing)

**Prediction:**
- TGP V_td_η.1 = A · λ_C³ · √((1-ρ̄)² + η̄²)  z η.1 triple
- TGP V_ts_η.1 = A · λ_C² (Wolfenstein O(λ⁴))
- Ratio |V_td/V_ts| = λ_C · √((1-ρ̄)² + η̄²) / 1 ≈ 0.225 · √((1-0.141)² + 0.357²)
  = 0.225 · √(0.737 + 0.127) = 0.225 · √0.865 ≈ 0.209
- PDG 2024: |V_td/V_ts| = 0.205 ± 0.006 (z B_s/B_d mixing) (drift ~2%)
- LHCb Run 4: σ(|V_td/V_ts|) ~ 1%
- **Falsification gate:** if |V_td/V_ts| outside [0.195, 0.220] at >5σ,
  Wolfenstein triple cascade broken
- **Confirmation gate:** if ratio ≈ 0.205-0.215 measured, B-B̄ mixing
  cross-sector lock z (ρ̄, η̄) confirmed

### T3.6 (H6) — 4-channel η.1 falsification convergence

**Prediction:**
- Belle II 2027+: |V_ub| (H1)
- LHCb Run 4 2030+: Jarlskog J (H2)
- Belle II + LHCb (combined): sin(2β) (H3)
- LHCb Run 4 2030+: |V_td/V_ts| (H5)
- **Convergence:** ≥ 3 z 4 channels muszą converge within 5σ TGP
  predictions for η.1 triple lock stabilization
- **Falsification gate:** if ≥ 2 z 4 channels reject η.1 > 5σ, triple
  classification PARTIALLY DERIVED → reverts to STRUCTURAL

---

## Verdict gate

**6/6 PASS** → η.1 program END, classification PARTIALLY DERIVED (refined),
ledger 427 → 445.

**5/6 PASS** → η.1 program END z minor gap.

**≤ 4/6 PASS** → η.1.Phase3 reframing required.

---

## Środowisko

```bash
PYTHONIOENCODING=utf-8 python -X utf8 research/op-eta-wolfenstein/phase3_eta_predictions.py 2>&1 | tee research/op-eta-wolfenstein/phase3_eta_predictions.txt
```

## Cross-references

- [`program.md`](program.md) — overall η.1 plan
- [`Phase1_results.md`](Phase1_results.md) — top-5 rationals ranked
- [`Phase2_results.md`](Phase2_results.md) — triple LOCKED
- [`../../PREDICTIONS_REGISTRY.md`](../../PREDICTIONS_REGISTRY.md) — H1-H6 entries (LIVE 2027+)
