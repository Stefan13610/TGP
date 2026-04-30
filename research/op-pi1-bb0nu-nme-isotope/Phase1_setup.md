---
title: "π.1.Phase1 setup — NME landscape (5 sub-tests)"
date: 2026-04-30
cycle: π.1.Phase1
status: PRE-EXECUTION
parent: "[[program.md]]"
predecessor: "[[../op-omicron1-sigmamnu-cosmo/Phase3_results.md]]"
tags:
  - TGP
  - pi1
  - phase1
  - NME
  - landscape
---

# π.1.Phase1 — NME landscape (5 sub-tests)

> **Cel:** Inventory NME literature (QRPA, IBM-2, NSM, EDF) and PSF
> values for Ge-76, Te-130, Xe-136; anchor m_ββ Form A/B z ν.1;
> compute method-spread σ; verify rate-only viability.

## Sub-tests

### P1.1 — NME inventory (QRPA/IBM-2/NSM/EDF × Ge-76/Te-130/Xe-136)

12-cell matrix; report mean ± span per isotope.

### P1.2 — Phase-space factor G_{0ν} inventory

Kotila & Iachello 2012 (10⁻¹⁵ yr⁻¹):
- Ge-76: 2.36
- Te-130: 14.22
- Xe-136: 14.58

g_A = 1.27 unquenched (literature default).

### P1.3 — m_ββ Form A=1.584 ‖ Form B=3.249 meV from ν.1

Anchored z [[../op-nu-majorana-phase-mbb/Phase3_results.md]].

### P1.4 — Method-spread σ_NME

For each isotope, compute σ = (max − min) / mean. Identifies which
isotope is most NME-stable.

### P1.5 — Rate-only viability gate

T_{1/2}^{0ν,TGP-A} > T_{1/2}^{exp,90%CL} and T_{1/2}^{0ν,TGP-B} >
T_{1/2}^{exp,90%CL} for all 3 isotopes × 4 NME methods. If yes, no rate
falsification yet → π.1 advances.

Current 90% CL bounds:
- Ge-76 (LEGEND-200 + Gerda): T_{1/2} > 1.8×10²⁶ yr
- Te-130 (CUORE): T_{1/2} > 2.2×10²⁵ yr
- Xe-136 (KamLAND-Zen 800): T_{1/2} > 2.3×10²⁶ yr

## PASS bramka

- ≥4/5 PASS → π.1.Phase2 OK
- 5/5 PASS → strong launch

## Środowisko

```bash
PYTHONIOENCODING=utf-8 python -X utf8 research/op-pi1-bb0nu-nme-isotope/phase1_nme_landscape.py 2>&1 | tee research/op-pi1-bb0nu-nme-isotope/phase1_nme_landscape.txt
```

## Cross-references

- [[program.md]]
- [[../op-nu-majorana-phase-mbb/Phase3_results.md]]
