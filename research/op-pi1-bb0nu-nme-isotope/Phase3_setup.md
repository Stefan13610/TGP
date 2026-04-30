---
title: "π.1.Phase3 setup — predictions + falsification (6 sub-tests)"
date: 2026-04-30
cycle: π.1.Phase3
status: PRE-EXECUTION
parent: "[[program.md]]"
predecessor: "[[Phase2_results.md]]"
tags:
  - TGP
  - pi1
  - phase3
  - 0nubb
  - predictions
  - falsification
---

# π.1.Phase3 — predictions + falsification (6 sub-tests)

> **Cel:** Project T_{1/2} Form A/B cascade against KZ-900 (2027+),
> LEGEND-1000 (2030+), nEXO (2030+), NEXT-HD (2030+), CUPID (2030+),
> SNO+ (2026+); register cross-isotope T_{1/2} ratio falsification
> as decisive π.1 closure; 7-channel falsification convergence.

## Sub-tests

### P3.1 — KamLAND-Zen 800 + KZ-900 2027+ Xe-136

KZ-800 current limit: T_{1/2} > 2.3×10²⁶ yr.
KZ-900 2027+ projected sensitivity: T_{1/2} ~ 5×10²⁷ yr.
Form A predicts T(Xe-136, QRPA) = 6.18×10²⁹ yr → 124× below target.
Form B predicts T(Xe-136, QRPA) = 1.47×10²⁹ yr → 29× below target.
→ **Both Forms safe** at KZ-900.

### P3.2 — LEGEND-1000 2030+ Ge-76

Sensitivity ~1.3×10²⁸ yr.
Form A T(Ge-76, EDF) = 2.08×10³⁰ yr → 160× below target.
Form B T(Ge-76, EDF) = 4.95×10²⁹ yr → 38× below target.
→ Both Forms safe; Ge-76 NME-cleanest isotope (σ_NME=0.465).

### P3.3 — nEXO 2030+ Xe-136 + NEXT-HD 2030+ Xe-136

nEXO projected: T_{1/2} ~ 5.7×10²⁷ yr (90% CL).
NEXT-HD projected: T_{1/2} ~ 4×10²⁷ yr (HP gas).
Form A m_ββ → m_ββ ~ 0.5 meV ↔ T ~ 1.7×10²⁹ yr safe.
**Discrimination Form A vs B: m_ββ separation 1.665 meV vs nEXO
target σ ≈ 0.5 meV → 3.33σ.**

### P3.4 — CUPID 2030+ Te-130 + SNO+ Te-130

CUPID-1T 2030+ sensitivity ~ 2×10²⁷ yr.
SNO+ Phase II: ~ 10²⁶ yr.
Form A T(Te-130, EDF) = 2.81×10²⁹ yr → 140× below CUPID target.
Form B T(Te-130, EDF) = 6.69×10²⁸ yr → 33× below.
→ Both Forms safe; Te-130 useful for cross-isotope.

### P3.5 — Cross-isotope T_{1/2} ratio falsification

Decisive 3-way Ge-76 / Te-130 / Xe-136 consistency check:

| Pair | TGP-EDF | Lit-mean range | Falsification trigger |
|------|---------|----------------|------------------------|
| T(Ge)/T(Xe) | 2.49 | 2.22–5.15 | exp ratio outside range falsifies NME methodology |
| T(Te)/T(Xe) | 0.71 | 0.66–0.92 | tightest cross-iso test |

If KZ-900 + LEGEND-1000 + CUPID detect 0νββ z mismatched ratios →
**either NME methods all wrong OR new physics (TGP m_ββ scenario shifts)**.

### P3.6 — 7-channel π.1 falsification convergence

| # | Channel | Observable | Action | Date |
|---|---------|------------|--------|------|
| 1 | KZ-800 / KZ-900 | T(Xe-136) | both forms safe | 2024-2027+ |
| 2 | LEGEND-1000 | T(Ge-76) | both forms safe | 2030+ |
| 3 | nEXO | T(Xe-136) | 3.33σ A vs B | 2030+ |
| 4 | NEXT-HD | T(Xe-136 HP) | 3σ A vs B | 2030+ |
| 5 | CUPID | T(Te-130) | both forms safe | 2030+ |
| 6 | SNO+ | T(Te-130) | mid-sensitivity | 2026+ |
| 7 | Theory NME-cross | T-ratio Ge/Te/Xe | NME systematic | 2026+ |

## PASS bramka

- ≥5/6 PASS → π.1 program END
- 6/6 PASS → FULL CONVERGENCE

## Środowisko

```bash
PYTHONIOENCODING=utf-8 python -X utf8 research/op-pi1-bb0nu-nme-isotope/phase3_nme_predictions.py 2>&1 | tee research/op-pi1-bb0nu-nme-isotope/phase3_nme_predictions.txt
```

## Cross-references

- [[program.md]]
- [[Phase2_results.md]]
- [[../op-nu-majorana-phase-mbb/Phase3_results.md]]
