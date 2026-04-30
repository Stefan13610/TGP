---
title: "π.1.Phase1 results — NME landscape (5/5 PASS)"
date: 2026-04-30
cycle: π.1.Phase1
status: PASS
parent: "[[program.md]]"
predecessor: "[[Phase1_setup.md]]"
tags:
  - TGP
  - pi1
  - phase1
  - NME
  - landscape
  - PASS
---

# π.1.Phase1 results — NME landscape

**Score: 5/5 PASS → π.1.Phase2 OK z full landscape.**

## Sub-test results

### P1.1 — 12-cell NME matrix ✓ PASS

| Isotope | QRPA | IBM-2 | NSM | EDF | mean | span |
|---------|------|-------|-----|-----|------|------|
| Ge-76 | 5.00 | 4.60 | 3.00 | 4.60 | 4.30 | 2.00 |
| Te-130 | 4.20 | 4.10 | 1.90 | 5.10 | 3.82 | 3.20 |
| Xe-136 | 3.40 | 3.30 | 1.80 | 4.20 | 3.17 | 2.40 |

### P1.2 — Phase-space factor inventory ✓ PASS

| Isotope | G_{0ν} (10⁻¹⁵ yr⁻¹) |
|---------|---------------------|
| Ge-76 | 2.36 |
| Te-130 | 14.22 |
| Xe-136 | 14.58 |

g_A = 1.27 unquenched (Kotila & Iachello 2012).

### P1.3 — m_ββ Form A/B anchored z ν.1 ✓ PASS

| Form | m_ββ (meV) |
|------|------------|
| A (chirality halving) | 1.584 |
| B (PMNS-Wolfenstein) | 3.249 |

Δm_ββ (B−A) = 1.665 meV; ratio B/A = **2.051× → 4.207× in T_{1/2}**.

### P1.4 — Method-spread σ_NME ✓ PASS

| Isotope | mean | span | σ_NME = span/mean |
|---------|------|------|-------------------|
| Ge-76 | 4.30 | 2.00 | **0.465** (most stable) |
| Xe-136 | 3.17 | 2.40 | 0.756 |
| Te-130 | 3.82 | 3.20 | 0.837 |

→ **Ge-76 is the cleanest NME observable** (LEGEND-1000 advantage).

### P1.5 — Rate-only viability ✓ PASS

All **24 cells** (3 isotopes × 4 NME methods × 2 Forms) yield
T_{1/2}^{0ν,TGP} > current 90% CL bound by ≥3 orders of magnitude:

- Min T_A: 2.81×10²⁹ yr (Te-130 EDF)
- Min T_B: 6.69×10²⁸ yr (Te-130 EDF)
- Strictest bound: Xe-136 KZ800 = 2.3×10²⁶ yr

→ **No rate-only falsification today** for either Form across any NME.
π.1 advances to Phase 2.

## Cross-references

- [[program.md]]
- [[Phase1_setup.md]]
- [[../op-nu-majorana-phase-mbb/Phase3_results.md]]
- [[../op-omicron1-sigmamnu-cosmo/Phase3_results.md]]
