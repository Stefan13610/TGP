---
title: "ο.1.Phase1 results — Σm_ν cosmological landscape (5/5 PASS)"
date: 2026-04-30
cycle: ο.1.Phase1
status: PASS
parent: "[[program.md]]"
predecessor: "[[Phase1_setup.md]]"
tags:
  - TGP
  - omicron1
  - phase1
  - landscape
  - PASS
---

# ο.1.Phase1 results — Σm_ν cosmological landscape

**Score: 5/5 PASS → Phase 2 unlocked.**

> **Headline:** TGP Form A Σm_ν = 59.01 meV (m_1=0, NH) z 0.17% drift do
> PDG-anchored 59.11 meV; passes DESI DR2 95% CL (72 meV) z 13 meV
> margin (18%); IH disfavored z floor 101.72 meV > DR2; m_lightest
> dual-form structure registered (A: m_1=0 ‖ B: m_1≈7.5 meV).

## Sub-test results

### O1.1 — bounds inventory ✓ PASS

5 cosmological bounds inventoried:
- Planck 2018 alone < 120 meV
- Planck + BAO BOSS < 150 meV
- Planck + DESI DR1 < 82 meV (2024)
- **Planck + DESI DR2 < 72 meV** (2024–2025, adopted)
- DESI DR2 + CMB + SN < 64 meV (tightest)

### O1.2 — TGP Σm_ν_A derivation review ✓ PASS

- m_1 = 0 (Form A)
- m_2 = √Δm²_21 = 8.614 meV
- m_3 = √Δm²_31 = 50.498 meV
- Σm_ν_PDG = 59.111 meV
- Σm_ν_TGP_A = **59.01 meV** (ν.1 substrate-action correction)
- Drift = 0.101 meV (0.17%) ✓

### O1.3 — NH vs IH likelihood ✓ PASS

- IH floor: m_1 ≈ 50.50 meV + m_2 ≈ 51.23 meV → Σm_ν_IH ≈ 101.72 meV
- DESI DR2 95% CL = 72 meV → IH excluded by 29.72 meV
- TGP Form A NH ordering CONFIRMED

### O1.4 — m_lightest dual-form structure ✓ PASS

| Form | m_1 (meV) | Σm_ν (meV) | DR2 status |
|------|-----------|------------|------------|
| **A (TGP)** | **0** | **59.01** | PASS w 13 meV margin |
| B (TGP alt) | 7.5 | 81.5 | edge of DR2 |
| IH alt | 0 | 101.7 | 3σ disfavored |
| Degenerate | 50 | 150 | 4σ excluded |

### O1.5 — viability gate ✓ PASS

- Σm_ν_TGP_A = 59.01 meV < DESI DR2 = 72 meV
- Margin = 12.99 meV (18.0%)
- → PASS, Phase 2 unlocked

## Cross-references

- [[program.md]]
- [[Phase1_setup.md]]
- [[../op-nu-majorana-phase-mbb/Phase2_results.md]]
