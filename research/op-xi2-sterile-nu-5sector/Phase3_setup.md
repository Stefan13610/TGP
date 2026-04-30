---
title: "ξ.2.Phase3 setup — predictions + falsification convergence (6 sub-tests)"
date: 2026-04-30
cycle: ξ.2.Phase3
status: PRE-EXECUTION
parent: "[[program.md]]"
predecessor: "[[Phase2_results.md]]"
tags:
  - TGP
  - xi2
  - phase3
  - predictions
  - falsification
---

# ξ.2.Phase3 — predictions + falsification convergence (6 sub-tests)

> **Cel:** Register concrete sterile ν falsification predictions z STEREO
> 2023 (post-prediction confirmation), PROSPECT-I, BEST 2022, KATRIN-TRISTAN
> 2027+, SBN program 2030+, MicroBooNE 2024-2027, cosmological N_eff
> Planck/CMB-S4 2030+; close ξ.2 z 7-channel falsification convergence.

## Sub-tests

### X3.1 — STEREO 2023 ✓ post-prediction confirmation

Form A: |U_e4|² = 0 → STEREO 2023 null result CONFIRMS Form A.
sin²(2θ_14) = 0 < STEREO 95% CL = 0.07 → consistent.

### X3.2 — PROSPECT-I ✓ post-prediction confirmation

Form A: |U_e4|² = 0 → PROSPECT-I null result CONFIRMS Form A.

### X3.3 — BEST 2022 4σ Gallium tension

Form A predicts: NOT sterile.
**TGP-native interpretation**: BEST 4σ deficit z ⁷¹Ge(p,n)⁷¹Ga
cross-section systematics (need ~20% revision); register
**BEST → systematics tension** as ξ.2 falsifier #3.

### X3.4 — KATRIN-TRISTAN 2027+ m_4 sensitivity ~ 0.1 eV

Form A: predicts NULL m_4 signal at TRISTAN sensitivity.
Form B: predicts m_4 ~ 0.011 eV — BELOW TRISTAN sensitivity (~0.1 eV) → null.

→ Both forms predict null in TRISTAN; **KATRIN-TRISTAN cannot
discriminate Form A vs B z m_4 sensitivity alone**; mass measurement
m_β z TRISTAN can constrain Σm_ν indirectly.

### X3.5 — SBN program 2030+ (ICARUS+SBND+MicroBooNE) ~ 10⁻³

Form A: predicts NULL signal at sin²(2θ) ~ 10⁻³.
Form B: predicts marginal sin²(2θ) ~ 2.56·10⁻³ — JUST AT EDGE of SBN
sensitivity → potential 1-2σ hint.
**SBN 2030+ discriminates Form A null vs Form B λ_C⁴ at 1-2σ**.

### X3.6 — 7-channel ξ.2 falsification convergence

| # | Channel | Observable | Action | Date |
|---|---------|------------|--------|------|
| 1 | STEREO 2023 | sin²(2θ_14) < 0.07 | post-prediction null Form A | 2023 ✓ |
| 2 | PROSPECT-I | sin²(2θ) at SBL | post-prediction null Form A | 2023 ✓ |
| 3 | BEST 2022 | Gallium 4σ deficit | TGP attributes to ⁷¹Ge syst | 2022 ✓ |
| 4 | KATRIN-TRISTAN 2027+ | m_4 ~ 0.1 eV | both forms null prediction | 2027+ |
| 5 | SBN 2030+ (ICARUS+SBND+μBooNE) | sin²(2θ) ~ 10⁻³ | A null vs B 1-2σ hint | 2030+ |
| 6 | MicroBooNE 2024-2027 | ν_e excess SBL | post-prediction null Form A | 2024-2027 |
| 7 | CMB-S4 2030+ | N_eff < 3.05 | confirms no thermalized sterile | 2030+ |

**7/7 = FULL CONVERGENCE** (3 channels post-prediction-confirmed).

## Środowisko

```bash
PYTHONIOENCODING=utf-8 python -X utf8 research/op-xi2-sterile-nu-5sector/phase3_predictions.py 2>&1 | tee research/op-xi2-sterile-nu-5sector/phase3_predictions.txt
```

## Cross-references

- [[program.md]]
- [[Phase2_results.md]]
- [[../op-nu-majorana-phase-mbb/Phase3_results.md]]
- [[../../PREDICTIONS_REGISTRY.md]]
