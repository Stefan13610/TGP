---
title: "Phase 0 balance — cluster mass deficit literature inventory + 6/6 gate"
date: 2026-05-11
parent: "[[./README.md]]"
type: phase-balance
phase: 0
status: 🟡 in progress (scaffold; needs literature deep-dive Phase 1)
tags:
  - phase0
  - balance-sheet
  - cluster-mass-deficit
  - literature-inventory
---

# Phase 0 — Balance sheet (scaffold version)

## §1 — Sub-needs inventory (initial, ~12 items)

### Theoretical sub-needs
- **N0.1**: Cluster-scale ROFM extension formalism (g_eff[{Φ_i}] aggregation
  od galaxy → galaxy_group → cluster)
- **N0.2**: Multi-source interaction enhancement mechanism w high-density limit
- **N0.3**: Transition scale ~Mpc emergent czy explicit (analog kpc galactic)
- **N0.4**: NFW-like vs cuspless profile prediction z TGP

### Numerical/observational sub-needs
- **N0.5**: Cluster sample dataset (Coma, A1689, Bullet, A2744, +large sample)
- **N0.6**: v_circ profile vs r prediction w TGP
- **N0.7**: Bullet Cluster lensing-vs-X-ray offset reproduction
- **N0.8**: Mass-temperature relation M-T_X cluster preserve

### BSM addition sub-needs (if H1a fails)
- **N0.9**: Sterile ν mass + mixing analysis (~2 eV, sin²2θ ~ 10⁻³)
- **N0.10**: Planck 2018 N_eff = 3.046 ± 0.18 bound check
- **N0.11**: BBN ⁴He Y_p sensitivity do additional ν
- **N0.12**: CMB-S4 N_eff future precision (±0.04) prediction

## §2 — Literature anchors needed (Phase 1 prep)

### Cluster mass deficit observational
- Coma cluster: Andernach 1991, Geller & Beers 1982; modern X-ray + lensing
- A1689: Limousin et al. 2007; A2744: Merten et al. 2011
- Bullet Cluster: Clowe et al. 2006 (collisionless DM evidence)
- General cluster sample: Vikhlinin et al. 2009 (M-T relation)

### MOND-like alternative resolutions
- Sanders 1999, 2003: MOND clusters problem
- Angus, Famaey, Diaferio 2010: MOND + sterile ν cluster fit
- Famaey, McGaugh 2012: MOND review (cluster scale issue explicit)

### TGP/emergent gravity cluster
- Verlinde 2017: emergent gravity cluster predictions (different framework)
- TGP galaxy_scaling: gs13-gs55 cluster results

### Sterile neutrino bounds
- Planck 2018 X. Cosmological parameters
- de Salas et al. 2018: global ν oscillation fit
- KATRIN 2024: m_ν direct measurement bounds

## §3 — 6/6 gate

| Gate | Status |
|---|---|
| **G1**: Predecessor cycles enumerated | ✅ N3 + galaxy_scaling cycles |
| **G2**: Sub-needs initial list | ✅ 12 items (§1) |
| **G3**: Literature currency check | 🟡 partial (anchors identified, full review Phase 1) |
| **G4**: Methodology binding declared | ✅ Native-first + S05 + ROFM extension |
| **G5**: Risk register R1-R6 declared | ✅ YAML frontmatter |
| **G6**: Probability assessment recorded | ✅ README §Probability |

**6/6 PASS (scaffold).** Phase 1 will deepen literature review + start formal
derivation.

## §4 — Methodology constraints (binding)

1. Native-first methodology
2. S05 single-Φ axiom preserved bezwarunkowo
3. §5.1 BD/Horndeski demarcation
4. **ROFM extension from N3 SPARC architecture** (recursive multi-source aggregation)
5. **NIE post-hoc fit** — TGP prediction must be derived strukturalnie z g_eff
6. Cluster scale ~Mpc transition czytane jako emergent property, NIE explicit cutoff
7. Sympy LOCK dla każdego analytic step (Phase 1+)
8. Honest STRUCTURAL_NO_GO jeśli ROFM extension napotka strukturalną przeszkodę

## §5 — Next session targets

**Phase 1 (next 1-2 sessions):**
1. Deep literature review (Vikhlinin 2009, Limousin 2007, Sanders 1999, Famaey 2012)
2. ROFM aggregation formalism od galaxy → cluster
3. g_eff[{Φ_i}] multi-source interaction explicit
4. Initial sympy LOCK 4-8 tests dla cluster mechanism

## §6 — Cross-references

- [[./README.md]]
- [[../op-L01-N3-SPARC-rho-consistency-2026-05-11/Phase_FINAL_close.md]] §5 (cluster outside N3 scope)
- [[../galaxy_scaling/]] gs13, gs21, gs36, gs47, gs55
- [[../op-emergent-metric-from-interaction-2026-05-09/]]
- Andernach 1991; Vikhlinin 2009; Clowe 2006; Famaey-McGaugh 2012; Planck 2018

---

**Phase 0 scaffold ready.** Cycle status: 🟡 OPEN — awaiting Phase 1 in
subsequent session.
