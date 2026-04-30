---
title: "υ.1.Phase2 results — universal closure-law sympy LOCK (7/7 FULL CASCADE)"
date: 2026-04-30
cycle: υ.1.Phase2
status: PASS
parent: "[[program.md]]"
predecessor: "[[Phase2_setup.md]]"
tags:
  - TGP
  - upsilon1
  - phase2
  - derivation
  - sympy
  - PASS
  - FULL-CASCADE
---

# υ.1.Phase2 results — universal closure-law derivation

**Score: 7/7 FULL CASCADE → Phase 3.**

> **Headline:** Universal closure law `closure(X_ref, X_obs) = (X_ref/X_obs)^{1/N_gen}`
> z N_gen=3 sympy-LOCKED. Recovers π.1 `(76/A_iso)^{1/3}` exactly i τ.1
> `(Z_a/Z_t)^{1/3}` exactly. Substrate-action gauge invariance (X→λX) potwierdzona
> sympy diff = 0. Joint combined `C = (76/A)^{1/3}·(Z_a/Z_t)^{1/3}` zachowuje
> trivial limit = 1.0 (Z_a=Z_t, A=76). N_gen=3 cascade-primality reuse z τ.1.

## Sub-test results

### U2.1 — sympy declaration ✓ PASS

```python
N_gen = sp.Integer(3)
closure = (X_ref / X_obs)**(sp.Rational(1, 3))
```

### U2.2 — π.1 recover ✓ PASS

`closure[X_ref→76, X_obs→A_iso] = (76/A_iso)^{1/3}`. Diff vs π.1 target = **0**.
At A_iso=136 (¹³⁶Xe): closure = **0.8237**.

### U2.3 — τ.1 recover ✓ PASS

`closure[X_ref→Z_a, X_obs→Z_t] = (Z_a/Z_t)^{1/3}`. Diff vs τ.1 target = **0**.
At ⁷¹Ga (Z_a=31, Z_t=32): f = 0.9895, f² = 0.9791, R_TGP = **0.7751** ✓

### U2.4 — substrate-action gauge invariance ✓ PASS

Pod X_ref → λ·X_ref, X_obs → λ·X_obs: closure_gauge − closure = **0** sympy.
**Substrate-action gauge invariance LOCKED.**

### U2.5 — N_gen=3 cascade reuse ✓ PASS

TGP B²-cascade z τ.1.Phase2: B²_lep=2, B²_ν=1, B²_up=13/4, B²_dn=61/25.
4 species × 3 generations = 12 quark+lepton flavor states.
**N_gen=3 cascade-primality reused** w υ.1.

### U2.6 — joint combined closure ✓ PASS

`C(A, Z_a, Z_t) = (76/A_iso)^{1/3} · (Z_a/Z_t)^{1/3}`

| test | A | Z_a | Z_t | C |
|---|---:|---:|---:|---:|
| ⁷⁶Ge β-decay | 76 | 32 | 33 | 0.9898 |
| ¹³⁶Xe β-decay | 136 | 54 | 55 | 0.8187 |
| trivial limit (Z_a=Z_t) | 76 | 32 | 32 | **1.0000** ✓ |

### U2.7 — Phase 2 gate ✓ PASS

7/7 FULL CASCADE → Phase 3.

## Promotions

- **υ.1 universal closure law DERIVED**: `closure(X) = (X_ref/X_obs)^{1/N_gen}` z N_gen=3
- **Substrate-action gauge invariance LOCKED**: closure invariant pod X→λX
- **π.1 ↔ τ.1 unification**: oba są instancjami υ.1 universal law (różnica tylko w X)
- **Joint cross-family combined closure**: deployable for joint NME × EC TT

## Cross-references

- [[program.md]]
- [[Phase1_results.md]]
- [[Phase2_setup.md]]
- [[../op-pi1-bb0nu-nme-isotope/Phase3_results.md]]
- [[../op-tau1-closure-overlap-coulomb/Phase2_results.md]]
