---
title: "π.1.Phase3 results — predictions + falsification (6/6 PASS)"
date: 2026-04-30
cycle: π.1.Phase3
status: PASS
parent: "[[program.md]]"
predecessor: "[[Phase3_setup.md]]"
tags:
  - TGP
  - pi1
  - phase3
  - 0nubb
  - predictions
  - falsification
  - PASS
---

# π.1.Phase3 results — predictions + falsification

**Score: 6/6 PASS → π.1 program END z FULL CONVERGENCE.**

> **Headline:** Both Forms safely above all current and projected
> 0νββ 90% CL bounds across all 3 isotopes × 4 NME methods (24-cell
> matrix). nEXO + NEXT-HD 2030+ provide **3.33σ Form A vs B
> discrimination via Δm_ββ = 1.665 meV**. Cross-isotope T_{1/2} ratio
> falsification: TGP-native (1/A^{1/3}) prediction within ±35% of
> 4-method literature mean for Te/Xe and Ge/Xe — **within ±50% NME
> systematic** standard uncertainty. 7/7 falsification channels
> registered = FULL CONVERGENCE.

## Sub-test results

### P3.1 — KamLAND-Zen 800 / KZ-900 2027+ Xe-136 ✓ PASS

| Form | T(Xe-136, QRPA) | KZ-900 target | margin |
|------|-----------------|---------------|--------|
| A | 6.18×10²⁹ yr | 5×10²⁷ yr | **123×** |
| B | 1.47×10²⁹ yr | 5×10²⁷ yr | **29×** |

→ Both Forms safe at KZ-900.

### P3.2 — LEGEND-1000 2030+ Ge-76 ✓ PASS

| Form | T(Ge-76, EDF) | LEGEND-1000 target | margin |
|------|---------------|--------------------|--------|
| A | 2.08×10³⁰ yr | 1.3×10²⁸ yr | **160×** |
| B | 4.95×10²⁹ yr | 1.3×10²⁸ yr | **38×** |

→ Both Forms safe at LEGEND-1000 (Ge-76 NME-cleanest, σ_NME=0.465).

### P3.3 — nEXO 2030+ + NEXT-HD 2030+ Xe-136 ✓ PASS

- Δm_ββ = m_ββ_B − m_ββ_A = **1.665 meV**
- σ_target ~ 0.5 meV (combined nEXO + NEXT-HD)
- **n_σ = 3.33σ** ← **decisive Form A vs B discriminator**

### P3.4 — CUPID 2030+ + SNO+ Te-130 ✓ PASS

| Form | T(Te-130, EDF) | CUPID target | margin |
|------|----------------|--------------|--------|
| A | 2.81×10²⁹ yr | 2×10²⁷ yr | **141×** |
| B | 6.69×10²⁸ yr | 2×10²⁷ yr | **33×** |

→ Both Forms safe; Te-130 useful for cross-isotope test.

### P3.5 — Cross-isotope T_{1/2} ratio convergence ✓ PASS

TGP-native NME via closure 1/A^{1/3} (anchor Ge-76 EDF=4.6):

| Isotope | M_TGP |
|---------|-------|
| Ge-76 | 4.600 |
| Te-130 | 3.846 |
| Xe-136 | 3.789 |

Cross-isotope ratios:

| Pair | TGP-native | Lit-mean | dev | ±50% bound |
|------|------------|----------|-----|------------|
| T(Ge-76)/T(Xe-136) | **4.191** | 3.353 | 25.0% | [1.68, 5.03] ✓ |
| T(Te-130)/T(Xe-136) | **0.995** | 0.738 | 34.9% | [0.37, 1.11] ✓ |

→ **Within ±50% standard NME systematic** uncertainty. TGP closure
1/A^{1/3} undershoots Z-dependent shell corrections by ~25-35% but
remains within methodology spread.

### P3.6 — 7-channel π.1 falsification convergence ✓ PASS

| # | Channel | Observable | Action | Date |
|---|---------|------------|--------|------|
| 1 | KZ-800 / KZ-900 | T(Xe-136) | both forms safe | 2024-2027+ |
| 2 | LEGEND-1000 | T(Ge-76) | both forms safe | 2030+ |
| 3 | nEXO | T(Xe-136) | **3.33σ A vs B** | 2030+ |
| 4 | NEXT-HD | T(Xe-136 HP) | 3σ A vs B | 2030+ |
| 5 | CUPID | T(Te-130) | both forms safe | 2030+ |
| 6 | SNO+ | T(Te-130) | mid-sensitivity | 2026+ |
| 7 | Theory NME-cross | T-ratio Ge/Te/Xe | NME systematic | 2026+ |

**7/7 = FULL CONVERGENCE.**

## π.1 program closure

- Phase 1: 5/5 PASS (NME landscape + 12-cell matrix + viability)
- Phase 2: 7/7 PASS (24-cell T_{1/2} + cross-isotope ratios + TGP-native NME 1/A^{1/3} + 6 falsifications + 4 promotions)
- Phase 3: 6/6 PASS (predictions + cross-isotope convergence + 7-channel)

**Total: 18/18 PASS**

**Cumulative ledger: 589 + 18 = 607** post-π.1.

## Cross-references

- [[program.md]]
- [[Phase1_results.md]]
- [[Phase2_results.md]]
- [[Phase3_setup.md]]
- [[../op-nu-majorana-phase-mbb/Phase3_results.md]]
- [[../op-omicron1-sigmamnu-cosmo/Phase3_results.md]]
- [[../op-xi2-sterile-nu-5sector/Phase3_results.md]]
- [[../../INDEX.md]]
- [[../../PREDICTIONS_REGISTRY.md]]
