---
title: "υ.1.Phase1 results — structural homology π.1 ↔ τ.1 (5/5 PASS)"
date: 2026-04-30
cycle: υ.1.Phase1
status: PASS
parent: "[[program.md]]"
predecessor: "[[Phase1_setup.md]]"
tags:
  - TGP
  - upsilon1
  - phase1
  - homology
  - PASS
---

# υ.1.Phase1 results — structural homology π.1 ↔ τ.1

**Score: 5/5 PASS → Phase 2.**

> **Headline:** Formal algebraic homology potwierdzona — π.1 `(A_anchor/A_iso)^{1/3}`
> i τ.1 `(Z_a/Z_t)^{1/3}` mają identyczną strukturę `(X_ref/X_obs)^{1/N_gen}`;
> alt-N_gen scan wyklucza wszystkie N_gen∈{1,2,4,5,6} → tylko N_gen=3 LOCKED;
> oba rooting (3D-volume π.1 + N_gen-cascade τ.1) niezależnie zbiegają w 1/3.

## Sub-test results

### U1.1 — formal algebraic homology ✓ PASS

| family | closure form | substrate X |
|---|---|---:|
| π.1 NME | `(A_anchor / A_iso)^{1/3}` | mass-number A |
| τ.1 overlap | `(Z_a / Z_t)^{1/3}` | charge Z |
| υ.1 universal | `(X_ref / X_obs)^{1/N_gen}`, N_gen=3 | substrate count X |

Identyczna struktura algebraiczna; różnica tylko w `X`.

### U1.2 — substrate-action root scan ✓ PASS

- **π.1 root**: V_nuc ∝ R³, A ∝ V → R ∝ A^{1/3}; surface-density 1/R ∝ A^{-1/3}.
  Exp **1/3 z 3D-volume geometrycznego**.
- **τ.1 root**: f_overlap = geometric mean GM(x_1,x_2,x_3) = (x_1·x_2·x_3)^{1/3}
  across N_gen=3 generations. Exp **1/3 z cascade primality**.

**Konwergencja**: oba routy independent dają 1/3 — projekcje wspólnej
substrate-action invariance pod X-dimensional gauge (3D-volume ≡ N_gen=3).

### U1.3 — alt-N_gen falsification ✓ PASS

| N_gen | 1/N_gen | matches π.1 | matches τ.1 |
|---:|---:|:---:|:---:|
| 1 | 1.000 | ✗ | ✗ |
| 2 | 0.500 | ✗ | ✗ |
| **3** | **0.333** | **✓** | **✓** |
| 4 | 0.250 | ✗ | ✗ |
| 5 | 0.200 | ✗ | ✗ |
| 6 | 0.167 | ✗ | ✗ |

**Unikalny N_gen=3 LOCKED** simultaneously dla obu rodzin.

### U1.4 — cross-family viability ✓ PASS

- π.1 specjalizacja: closure(76, 136) = **0.8237** → M(¹³⁶Xe)/M(⁷⁶Ge) ~17.5% reduction
- τ.1 specjalizacja: closure(31, 32) = **0.9895** → f² = 0.9791, R_TGP = **0.7751** (vs BEST 0.8084±0.0295, 1.13σ ✓)

**Zero free parameters** — N_gen=3 LOCKED; X_ref/X_obs z danych.

### U1.5 — Phase 1 gate ✓ PASS

5/5 PASS — υ.1.Phase1 z FULL CASCADE → proceed Phase 2.

## Cross-references

- [[program.md]]
- [[Phase1_setup.md]]
- [[../op-pi1-bb0nu-nme-isotope/Phase3_results.md]]
- [[../op-tau1-closure-overlap-coulomb/Phase3_results.md]]
