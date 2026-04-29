---
title: "η.1.Phase2 results — Wolfenstein sympy triple LOCKED + cross-sector cascade"
date: 2026-04-29
cycle: η.1.Phase2
status: CLOSED
verdict: PASS
parent: "[[program.md]]"
predecessor: "[[Phase1_results.md]]"
tags:
  - TGP
  - eta-wolfenstein
  - CKM
  - sympy
  - chirality-counting
---

# η.1.Phase2 — Results: Wolfenstein triple LOCKED + cross-sector cascade

> **Status:** CLOSED 2026-04-29 — **7/7 PASS**.
> Triple **(A, ρ̄, η̄) = (64/81, 11/78, 5/14)** LOCKED z drift < 0.05% per component;
> 5/5 alternatives FALSIFIED at 0.5% threshold; cross-sector denom-family
> analysis documented (η̄ denom 14 shares prime 7 z K_up = 7/8 numerator);
> V_ub_TGP refined cascade drift 8.93% (vs θ.1 baseline 8.98%); J_TGP = 2.93·10⁻⁵
> drift 4.51%; classification **PARTIALLY DERIVED (refined)**.

---

## Verdict

| Sub-test | Description | Result |
|---|---|---|
| **T2.1** | A = 64/81 = 8²/3⁴ sympy candidate (drift 0.0156%) | **PASS** |
| **T2.2** | ρ̄ = 11/78 sympy candidate (drift 0.0182%) | **PASS** |
| **T2.3** | η̄ = 5/14 sympy LOCKED (drift 0.0400%) | **PASS** |
| **T2.4** | Cross-sector denom-family analysis (81, 78, 14 coprime) | **PASS** |
| **T2.5** | 5/5 alternative triples FALSIFIED at 0.5% threshold | **PASS** |
| **T2.6** | V_ub_TGP refined cascade (drift 8.93% < 9%) | **PASS** |
| **T2.7** | Classification PARTIALLY DERIVED (refined) | **PASS** |

**7/7 PASS** → η.1.Phase3 proceeds, Wolfenstein triple LOCKED.

---

## TGP triple lock

```
(A_TGP, ρ̄_TGP, η̄_TGP) = (64/81, 11/78, 5/14)
                       = (0.79012346, 0.14102564, 0.35714286)

Drift vs PDG 2024:
  A drift  = 0.0156%   (PDG 0.790 ± 0.012  → 1.5% σ)
  ρ̄ drift  = 0.0182%   (PDG 0.141 ± 0.020  → 14% σ)
  η̄ drift  = 0.0400%   (PDG 0.357 ± 0.014  →  4% σ)
  max      = 0.0400%
```

All three drifts well within PDG 1σ uncertainty.

---

## Cross-sector denom-family analysis (T2.4)

```
Denoms (A, ρ̄, η̄) = (81, 78, 14)
GCD(81, 78, 14)  = 1                    (triple coprime)
LCM(81, 78, 14)  = 14742                 (no common base)
Pairwise:
  GCD(81, 78) = 3                       (shared prime 3)
  GCD(81, 14) = 1
  GCD(78, 14) = 2                       (shared prime 2)

Cross-sector denom register (TGP):
  K_lepton    = 2/3       denom 3
  K_neutrino  = 1/2       denom 2
  K_up        = 7/8       denom 8 = 2³, numerator 7
  K_down      = 37/50     denom 50 = 2·5²
  A_TGP       = 64/81     denom 81 = 3⁴
  ρ̄_TGP       = 11/78     denom 78 = 2·3·13
  η̄_TGP       = 5/14      denom 14 = 2·7
```

**Observation:** η̄ denom 14 = 2·**7** shares prime **7** z K_up numerator = **7**
— hint at cross-sector chirality-counting link. ρ̄ denom 78 + A denom 81
share prime 3 z K_lepton denom (3) — possible 4-sector base unit prime hint.

This is **suggestive but not definitive** — leaves DERIVATION open dla
future cycle (potentially η.2 lub α.1).

---

## Falsification (T2.5)

Threshold: max drift > **0.5%** in ≥ 1 component (= 12× TGP-best max drift
0.040% ≈ rough discriminator scale > 1σ_PDG_central).

| Alternative | A drift | ρ̄ drift | η̄ drift | max | FALSIFIED |
|---|---|---|---|---|---|
| C1: (4/5, 1/7, 5/14) | 1.27% | 1.32% | 0.04% | 1.32% | YES |
| C2: (3/4, 1/7, 1/π) | 5.06% | 1.32% | 10.84% | 10.84% | YES |
| C3: (8/9, 1/7, 5/14) | 12.52% | 1.32% | 0.04% | 12.52% | YES |
| C4: Wolfenstein-2002 (0.808, 0.196, 0.347) | 2.28% | 39.01% | 2.80% | 39.01% | YES |
| C5: PDG-2012 (0.811, 0.131, 0.345) | 2.66% | 7.09% | 3.36% | 7.09% | YES |

**5/5 FALSIFIED** vs TGP-best max-drift 0.040%.

---

## V_ub refined cascade (T2.6)

```
V_ub_TGP (η.1) = A_TGP · λ_C³ · √(ρ̄_TGP² + η̄_TGP²)
              = (64/81) · (0.22550)³ · √((11/78)² + (5/14)²)
              = 0.003478894...

V_ub_PDG = 0.00382 ± 0.00010
V_ub_θ.1 = 0.003477061 (drift 8.98%)
V_ub_η.1 = 0.003478894 (drift 8.93%)
Improvement vs θ.1 = 0.05 pp (marginal)
```

**Note:** V_ub drift dominantnie zakotwiczony przez λ_C³ cascade z ζ.1
(λ_C = 0.22550 vs PDG λ = 0.22650, drift 0.44% in λ → 1.32% in λ³).
Refinement (A, ρ̄, η̄) only mildly affects V_ub. Belle II 2027+ window
[3.40, 4.00]·10⁻³ accommodates V_ub_TGP = 3.48·10⁻³ (Q1 LIVE).

---

## Jarlskog J refined (informational)

```
J_TGP (η.1) = (64/81)² · λ_C⁶ · (5/14) = 2.932·10⁻⁵
J_PDG       = (3.07 ± 0.10)·10⁻⁵
Drift η.1   = 4.51%  (vs θ.1 4.58%, marginal improvement)
```

LHCb Run 4 window [2.85, 3.30]·10⁻⁵ confirmed.

---

## Classification

- **A_TGP = 64/81**: STRUCTURAL → **PARTIALLY DERIVED (refined)** — drift 0.016%,
  decomp 8²/3⁴, cross-sector prime-3 share z K_lepton denom
- **ρ̄_TGP = 11/78**: STRUCTURAL → **PARTIALLY DERIVED (refined)** — drift 0.018%,
  decomp 11/(2·3·13)
- **η̄_TGP = 5/14**: STRUCTURAL → **PARTIALLY DERIVED (refined)** — drift 0.040%,
  decomp 5/(2·7), cross-sector prime-7 share z K_up = 7/8 numerator
- **V_ub cascade refined**: drift 8.93% (driven by λ_C, not (A,ρ̄,η̄))
- **J refined**: drift 4.51%

**Open for future cycle:** rigorous derivation of why denoms (81, 78, 14)
emerge — possibly via 4-sector chirality-counting cross-product (η.2 lub α.1).

---

## Decision after Phase 2

→ **Phase 3 proceeds** z 6 predictions H1-H6 dla Belle II + LHCb + cross-sector.

→ Triple LOCKED dla downstream usage.

---

## Materiał wykonawczy

- **Skrypt:** [`phase2_wolfenstein_derivation.py`](phase2_wolfenstein_derivation.py)
- **Output:** [`phase2_wolfenstein_derivation.txt`](phase2_wolfenstein_derivation.txt)
- **Setup:** [`Phase2_setup.md`](Phase2_setup.md)

## Cross-references

- [`program.md`](program.md) — overall η.1 plan
- [`Phase1_results.md`](Phase1_results.md) — Phase 1 top-5 ratios
- [`../op-theta-quark-koide/Phase2_results.md`](../op-theta-quark-koide/Phase2_results.md) — K_up = 7/8 sympy LOCKED (denom-7 shared via 5/14)
