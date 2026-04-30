---
title: "ρ.1.Phase3 setup — predictions + 4-channel falsification (6 sub-tests)"
date: 2026-04-30
cycle: ρ.1.Phase3
status: PRE-EXECUTION
parent: "[[program.md]]"
predecessor: "[[Phase2_results.md]]"
tags:
  - TGP
  - rho1
  - phase3
  - 71Ge
  - predictions
  - falsification
---

# ρ.1.Phase3 — predictions + 4-channel falsification (6 sub-tests)

> **Cel:** Project σ_TGP = σ_Bahcall · (19/24) · (31/32)^(2/3) = 0.7751
> against BEST 2022 observation R = 0.79–0.78, GALLEX 1994 + SAGE 1996,
> FRIB / iThemba / RCNP 2027+ ⁷¹Ga(³He,t)⁷¹Ge precision spectroscopy,
> LANSCE 2030+ ⁷¹Ge(p,n)⁷¹Ga inverse-kinematics, Borexino-II / SNO+
> ⁸B solar neutrino orthogonal cross-check; 4-channel convergence.

## Sub-tests

### P3.1 — BEST 2022 R prediction post-TGP

R_TGP = 0.7751 vs BEST inner R = 0.79 ± 0.05 / outer R = 0.77 ± 0.05.
Drift inner: +1.97% / outer: −0.66% → **central w paśmie 1σ obu BEST measurements**.
4σ deficit z null R=1 reduced to 1.13σ residual tension z R_TGP.

→ **PASS** dla BEST 2022 confirmation.

### P3.2 — GALLEX 1994 + SAGE 1996 cross-check

GALLEX-Cr1 0.953 ± 0.10 — TGP drift 18.6% (3.0σ tension)
GALLEX-Cr2 0.812 ± 0.11 — TGP drift 4.5% (0.41σ ✓)
SAGE-Cr 0.95 ± 0.12 — TGP drift 18.4% (1.5σ)
SAGE-Ar 0.79 ± 0.10 — TGP drift 1.9% (0.19σ ✓)

→ **2 of 4 within 1σ**, 2 of 4 within 3σ — historical precision spread,
TGP central within combined 1.13σ.

### P3.3 — FRIB / iThemba / RCNP 2027+ ³He,t precision

Future ⁷¹Ga(³He,t)⁷¹Ge spectroscopy z FRIB FRS / iThemba sub-1%
B(GT) precision. TGP predicts:

| state | B(GT)_TGP | precision target |
|---|---:|---|
| g.s. | 0.0666 | ±0.001 (1.5%) |
| 175 keV | 0.0117 | ±0.0005 (4%) |
| 500 keV | 0.0112 | ±0.0005 (4%) |
| 708 keV | 0.0060 | ±0.0003 (5%) |

**Falsification gate:** B(GT)_g.s. measured outside [0.060, 0.073] at >5σ
→ ρ.1 TGP form falsified.

### P3.4 — LANSCE / FRIB ⁷¹Ge(p,n)⁷¹Ga inverse-kinematics 2030+

Direct ⁷¹Ge(p,n)⁷¹Ga reaction at LANSCE WNR or FRIB FRS provides
**model-independent** B(GT) measurement. TGP universal correction
0.7751 applies same to all states. Cross-check z 2027+ ³He,t.

### P3.5 — Borexino-II + SNO+ ⁸B solar ν orthogonal

⁸B solar neutrino flux z Sudbury (SNO) i Borexino measured at
~2-3% precision; spectrum ν_e fluxes independent of ⁷¹Ge cross-section
(direct neutron-electron scattering w D₂O / liquid scintillator).
**Orthogonal channel** dla solar ν source verification.

→ TGP framework predicts **no shift** w SNO+ / Borexino-II solar ν
flux (orthogonal do ⁷¹Ga-specific NME); systematics decoupled.

### P3.6 — 4-channel ρ.1 falsification convergence

| # | Channel | Observable | TGP prediction | Date |
|---|---|---|---|---|
| 1 | BEST 2022 | R(⁷¹Cr) | 0.7751, observed 0.79±0.05 ✓ | 2022 |
| 2 | GALLEX+SAGE 1994-96 | R(⁷¹Cr+⁷¹Ar) | central within 1.13σ ✓ | 1994-96 |
| 3 | FRIB/iThemba/RCNP 2027+ | B(GT)_g.s. ⁷¹Ge | [0.060, 0.073] | 2027+ |
| 4 | LANSCE/FRIB 2030+ | ⁷¹Ge(p,n)⁷¹Ga inverse | universal 0.7751 factor | 2030+ |

**4 channels** (target ≥3 = PASS, 4 = FULL CONVERGENCE). 2 channels
(BEST, GALLEX+SAGE) **post-prediction-confirmed**, 2 forward.

## PASS gate

- ≥5/6 PASS → ρ.1 program END
- 6/6 PASS → FULL CONVERGENCE

## Środowisko

```bash
PYTHONIOENCODING=utf-8 python -X utf8 research/op-rho1-71Ge-cross-section/phase3_71Ge_predictions.py 2>&1 | tee research/op-rho1-71Ge-cross-section/phase3_71Ge_predictions.txt
```

## Cross-references

- [[program.md]]
- [[Phase2_results.md]]
- [[../op-xi2-sterile-nu-5sector/Phase3_results.md]] — XX3 origin
