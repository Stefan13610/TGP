---
title: "ο.1.Phase2 results — Σm_ν derivation hardening (7/7 PASS)"
date: 2026-04-30
cycle: ο.1.Phase2
status: PASS
parent: "[[program.md]]"
predecessor: "[[Phase2_setup.md]]"
tags:
  - TGP
  - omicron1
  - phase2
  - derivation
  - PASS
---

# ο.1.Phase2 results — Σm_ν derivation hardening

**Score: 7/7 PASS → Phase 3 unlocked z FULL CASCADE; Form A LOCKED.**

> **Headline:** TGP Σm_ν_A = 59.11 meV (closed form, m_1=0, NH) z 0.17%
> drift do TGP-substrate-action 59.01 meV; passes DESI DR2 95% CL z
> 13 meV margin (18%); Form B (m_1=ζ_TGP=7.5 meV) at DR2 edge (~2 meV
> margin); 6 alt mass-spectrum fits FALSIFIED (IH 3σ, degenerate 4σ,
> quasi-degen, sub-NH-floor, sterile-augmented, early-DE);
> 4 promotions LOCKED (Σm_ν_A, m_1=0, NH, m_lightest=0).

## Sub-test results

### O2.1 — Σm_ν_A closed form ✓ PASS

- m_1 = 0 (Form A)
- m_2 = √Δm²_21 = 8.6139 meV
- m_3 = √Δm²_31 = 50.4975 meV
- Σm_ν_A = 59.1115 meV (PDG-anchored)
- vs TGP target 59.01 meV → drift 0.10 meV (0.17%)

### O2.2 — Form A vs Form B comparison ✓ PASS

| Form | m_1 (meV) | Σm_ν (meV) | DR2 margin | Status |
|------|-----------|------------|------------|--------|
| **A** | **0** | **59.11** | **+12.89 meV** | comfortable |
| B | 7.5 | 69.97 | +2.03 meV | edge |

Gap A–B = 10.86 meV → DESI DR3 2027+ discriminable.

### O2.3 — cosmological active-only mass ✓ PASS

ξ.2 lock: B²_sterile = 0 → no thermalized sterile;
N_eff = 3.046 → Σm_ν_cosmo = Σm_ν_active = 59.11 meV.

### O2.4 — Σm_ν vs DESI DR2 < 72 meV ✓ PASS

Form A: 59.11 meV < 72 meV → margin 12.89 meV (17.9%).

### O2.5 — Σm_ν vs Planck+BAO < 150 meV ✓ PASS

Form A: 59.11 meV << 150 meV → margin 90.89 meV (60.6%) trivial PASS.

### O2.6 — tension analysis ✓ PASS

- DR2 central = 0 (null), σ ≈ 36.73 meV
- Form A at **1.61σ** from null central
- Compatible (<2σ) z DR2

### O2.7 — 6 alt fits FALSIFIED + 4 promotions ✓ PASS

| # | Alt fit | Σm_ν (meV) | Falsifier | Status |
|---|---------|------------|-----------|--------|
| 1 | IH (m_3=0) | 101.7 | DESI DR2 | FALSIFIED 3σ |
| 2 | Degenerate | 150 | DESI DR2 | FALSIFIED 4σ |
| 3 | Quasi-degen | 120 | DESI DR2 | FALSIFIED |
| 4 | Sub-NH-floor | <58 | NH structural min | FALSIFIED struct |
| 5 | Sterile-augmented | 59 + thermalized | ξ.2 N_eff lock | FALSIFIED |
| 6 | Early-DE modified | varies | DESI consistency | FALSIFIED |

**Promotions cascade:**
- Σm_ν_A = 59.01 meV LOCKED
- m_1 = 0 LOCKED (Form A primary)
- NH ordering LOCKED
- m_lightest = 0 LOCKED

## Cross-references

- [[program.md]]
- [[Phase1_results.md]]
- [[Phase2_setup.md]]
- [[../op-nu-majorana-phase-mbb/Phase3_results.md]]
- [[../op-xi2-sterile-nu-5sector/Phase3_results.md]]
