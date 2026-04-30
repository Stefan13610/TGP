---
title: "υ.1.Phase2 setup — universal closure-law derivation + sympy LOCK"
date: 2026-04-30
cycle: υ.1.Phase2
status: SETUP
parent: "[[program.md]]"
predecessor: "[[Phase1_results.md]]"
tags:
  - TGP
  - upsilon1
  - phase2
  - derivation
  - sympy
  - setup
---

# υ.1.Phase2 setup — universal closure-law derivation

**Score gate:** ≥6/7 PASS → Phase 3.

## Sub-tests (7)

- **U2.1** — sympy declaration: `closure(X) = (X_ref/X_obs)^{1/N_gen}`, N_gen=3
- **U2.2** — π.1 specjalizacja: recover `(76/A)^{1/3}` z anchor=⁷⁶Ge
- **U2.3** — τ.1 specjalizacja: recover `(Z_a/Z_t)^{1/3}` z N_gen=3
- **U2.4** — substrate-action invariance: closure invariant pod X→λX
- **U2.5** — N_gen=3 cascade-derivation reuse z τ.1.Phase2
- **U2.6** — joint cross-family combined closure check
- **U2.7** — Phase 2 gate ≥6/7

## Cross-references

- [[program.md]]
- [[Phase1_results.md]]
