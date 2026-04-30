---
title: "ο.1.Phase3 setup — Σm_ν predictions + falsification (6 sub-tests)"
date: 2026-04-30
cycle: ο.1.Phase3
status: PRE-EXECUTION
parent: "[[program.md]]"
predecessor: "[[Phase2_results.md]]"
tags:
  - TGP
  - omicron1
  - phase3
  - predictions
  - falsification
---

# ο.1.Phase3 — predictions + falsification (6 sub-tests)

> **Cel:** Register 7-channel cosmological Σm_ν falsification convergence
> z DESI DR2 (2024–2025 confirmed), DESI DR3 2027+, Simons Observatory
> 2025+, CMB-S4 2030+, KATRIN 2030+ direct mass, Euclid + Roman 2027–2030+.

## Sub-tests

### O3.1 — DESI DR3 2027+ ~50 meV projection
Form A 59.01 meV vs ~50 meV: tension ~9 meV → 2σ.

### O3.2 — Simons Observatory 2025+ projection
~40 meV CMB lensing + delensing target.

### O3.3 — CMB-S4 2030+ ~20 meV ultimate
Form A 59.01 meV decisive 5σ falsification ALBO confirmation.

### O3.4 — KATRIN 2030+ direct mass m_β ~ 0.2 eV
Orthogonal beta-decay endpoint; Form A m_β ≈ √Σm_i²·|U_ei|² ~ 9 meV → null.

### O3.5 — Euclid + Roman 2027–2030+ weak lensing
~30 meV galaxy clustering + lensing combined.

### O3.6 — 7-channel ο.1 cosmological convergence

| # | Channel | Observable | Action | Date |
|---|---------|------------|--------|------|
| 1 | DESI DR2 | Σm_ν<72 meV | Form A 59.01 PASS | 2024–2025 ✓ |
| 2 | DESI DR3 | Σm_ν<~50 meV | Form A 2σ tension | 2027+ |
| 3 | Simons Obs | Σm_ν<~40 meV | Form A edge | 2025+ |
| 4 | CMB-S4 | Σm_ν<~20 meV | Form A 5σ decisive | 2030+ |
| 5 | KATRIN | m_β<0.2 eV | Form A null | 2030+ |
| 6 | Euclid+Roman | weak lens+gal | Form A consistent | 2027–2030+ |
| 7 | LiteBIRD | CMB pol | r-Σm_ν cross | 2030+ |

## PASS bramka
- ≥5/6 PASS → ο.1 program END
- 6/6 PASS → FULL CONVERGENCE

## Środowisko

```bash
PYTHONIOENCODING=utf-8 python -X utf8 research/op-omicron1-sigmamnu-cosmo/phase3_sigmamnu_predictions.py 2>&1 | tee research/op-omicron1-sigmamnu-cosmo/phase3_sigmamnu_predictions.txt
```

## Cross-references

- [[program.md]]
- [[Phase2_results.md]]
- [[../op-xi2-sterile-nu-5sector/Phase3_results.md]]
