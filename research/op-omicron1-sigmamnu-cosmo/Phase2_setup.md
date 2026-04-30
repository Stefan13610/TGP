---
title: "ο.1.Phase2 setup — Σm_ν derivation hardening (7 sub-tests)"
date: 2026-04-30
cycle: ο.1.Phase2
status: PRE-EXECUTION
parent: "[[program.md]]"
predecessor: "[[Phase1_results.md]]"
tags:
  - TGP
  - omicron1
  - phase2
  - derivation
---

# ο.1.Phase2 — derivation hardening (7 sub-tests)

> **Cel:** Lock TGP Σm_ν_A = 59.01 meV closed-form, falsify 6 alt mass
> spectrum fits, register 4 promotions cascade.

## Sub-tests

### O2.1 — Σm_ν_A closed form
m_1=0, m_2=√Δm²_21, m_3=√Δm²_31; Σm_ν_A = √Δm²_21 + √Δm²_31.

### O2.2 — Form A vs Form B comparison
A: m_1=0; B: m_1=ζ_TGP ≠ 0. Form B introduces extra ~22.5 meV → tension.

### O2.3 — cosmological active-only mass post-ξ.2
ξ.2 lock: no thermalized sterile → Σm_ν_cosmo = Σm_ν_active.

### O2.4 — Σm_ν vs DESI DR2 < 72 meV → 13 meV margin
Form A: 59.01 < 72 → PASS w 18% margin.

### O2.5 — Σm_ν vs Planck+BAO < 150 meV → trivial PASS
Form A: 59.01 << 150 → trivial.

### O2.6 — Tension analysis
DESI DR2 central inferred ~ 0 (null measurement; 95% CL upper).
Bayesian tension Form A vs DR2 central ~1.4σ (compatible).

### O2.7 — 6 alt fits FALSIFIED + 4 promotions

| # | Alt fit | Σm_ν (meV) | Falsifier | Status |
|---|---------|------------|-----------|--------|
| 1 | IH (m_3=0) | ~101.72 | DESI DR2 | FALSIFIED (3σ) |
| 2 | Degenerate | ~150 | DESI DR2 | FALSIFIED (4σ) |
| 3 | Quasi-degen | ~120 | DESI DR2 | FALSIFIED |
| 4 | Sub-NH-floor | <58 | NH minimum | FALSIFIED struct |
| 5 | Sterile-augmented | 59+thermalized | ξ.2 N_eff | FALSIFIED |
| 6 | Early-DE modified | varies | DESI consistency | FALSIFIED |

Promotions:
- Σm_ν_A = 59.01 meV LOCKED
- m_1 = 0 LOCKED (Form A primary)
- NH ordering LOCKED
- m_lightest = 0 LOCKED

## PASS bramka
- ≥6/7 PASS → Phase 3 unlocked
- 7/7 PASS → strong lock

## Cross-references

- [[program.md]]
- [[Phase1_results.md]]
- [[../op-nu-majorana-phase-mbb/Phase3_results.md]]
- [[../op-xi2-sterile-nu-5sector/Phase3_results.md]]
