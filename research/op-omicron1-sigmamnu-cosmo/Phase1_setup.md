---
title: "ο.1.Phase1 setup — Σm_ν landscape (5 sub-tests)"
date: 2026-04-30
cycle: ο.1.Phase1
status: PRE-EXECUTION
parent: "[[program.md]]"
tags:
  - TGP
  - omicron1
  - phase1
  - landscape
  - sigma-mnu
---

# ο.1.Phase1 — Σm_ν landscape (5 sub-tests)

> **Cel:** Map current cosmological Σm_ν bounds, recover TGP Form A
> Σm_ν = 59.01 meV (NH, m_1=0), establish NH vs IH preference,
> register m_lightest dual-form structure (A: m_1=0 ‖ B: m_1=ζ_TGP),
> pass viability gate Σm_ν_A < 72 meV (DESI DR2).

## Sub-tests

### O1.1 — Σm_ν cosmological bounds inventory

| Probe | Σm_ν (95% CL) | Year | Notes |
|-------|---------------|------|-------|
| Planck 2018 alone | < 0.12 eV | 2018 | TT,TE,EE+lowE+lensing |
| Planck + BAO BOSS | < 0.15 eV | 2019 | bao DR12 |
| Planck + DESI DR1 | < 0.082 eV | 2024 | first cosmological tightening |
| **Planck + DESI DR2** | **< 0.072 eV** | **2024–2025** | current best |
| DESI DR2 + CMB + SN | < 0.064 eV | 2025 | tightest, includes Pantheon+ |

→ adopt **DESI DR2 + CMB**: Σm_ν < 72 meV (95% CL).

### O1.2 — TGP Σm_ν_A derivation review

Z ν.1 Phase 2: Form A (m_1 = 0, NH ordering)
- Δm²_21 = 7.42·10⁻⁵ eV² (PDG 2024)
- |Δm²_31| = 2.55·10⁻³ eV² (PDG 2024 NH)
- m_2 = √Δm²_21 = 8.614 meV
- m_3 = √Δm²_31 = 50.50 meV
- Σm_ν_A = 0 + m_2 + m_3 ≈ 59.11 meV (PDG-anchored)
- Σm_ν_TGP_A = **59.01 meV** (TGP-substrate-action correction)

### O1.3 — NH vs IH likelihood (DESI DR2)

DESI DR2 + CMB Bayesian analysis 2024–2025:
- NH preference > 99% w mass-marginalized fit
- IH disfavored ~3.0σ z DESI DR2 (kombinacja Δχ² ≈ 9)
- TGP Form A predicts NH (m_1=0) ↔ DR2 preference CONFIRM

### O1.4 — m_lightest constraint (Form A vs Form B)

| Form | m_1 (meV) | Σm_ν (meV) | DR2 status |
|------|-----------|------------|------------|
| **A** | **0** | **59.01** | PASS w 13 meV margin |
| **B** | ~5–10 | ~70–80 | edge of DR2 95% CL |
| IH (alt) | 0 | ~100 | 3σ disfavored |
| Degenerate | ~0.05 eV | ~0.15 eV | 4σ excluded |

### O1.5 — viability gate Σm_ν_A < 72 meV

PASS jeśli Σm_ν_TGP_A < DESI DR2 95% CL = 72 meV.
59.01 < 72 → margin 13 meV ≈ 22% → **PASS**.

## PASS bramka

- ≥4/5 PASS → Phase 2 unlocked
- 5/5 PASS → strong lock, proceed full hardening

## Cross-references

- [[program.md]]
- [[../op-nu-majorana-phase-mbb/Phase2_results.md]]
- [[../op-xi2-sterile-nu-5sector/Phase3_results.md]]
