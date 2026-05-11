---
title: "NEEDS — op-cluster-mass-deficit-resolution-2026-05-11 (otwarte items + future tests)"
date: 2026-05-11
parent: "[[./README.md]]"
type: needs
status: 🟢 cycle CLOSED (H1b adopted); otwarte items są empirical future tests
classification: STRUCTURAL_DERIVED (H1b)
tags:
  - needs
  - deferred-precision
  - future-tests
  - sterile-nu-direct-detection
---

# NEEDS — otwarte items (post-cluster-cycle closure)

## §0 — Status

Cykl cluster mass deficit zamknięty 2026-05-11 jako **STRUCTURAL_DERIVED**
(H1b: TGP + sterile ν 2 eV) z **24/24 sympy PASS** + **6/6 P-requirements
RESOLVED** + **6.4σ combined multi-experiment falsifiability post-2030+**.

Poniższe items są **future empirical tests** + deferred precision, NIE blokery.

## §1 — Future empirical tests (falsifiable)

### §1.1 — JSNS² (J-PARC 2024-2027) — EARLIEST PROBE

**Status:** active experiment 2024+

**H1b boundary prediction:**
- sin²2θ_14 ~ 10⁻³ at Δm² ≈ 4 eV²
- JSNS² sensitivity at 10⁻³ — **exactly at boundary**

**Falsifiability:**
- Detection at sin²2θ ≈ 10⁻³: **H1b PARTIALLY CONFIRMED**
- Null at sin²2θ < 10⁻³: **H1b sin²2θ parameter wymaga refinement**

### §1.2 — CMB-S4 (DOE+NSF 2030+) — N_eff + Σm_ν

**Status:** scheduled deployment 2030+

**H1b predictions:**
- ΔN_eff = 0.05 ± 0.01 → 1.25σ marginal detection
- Σm_ν NH = 0.062 eV → 1.55σ detection
- Σm_ν IH = 0.102 eV → 2.55σ detection

**Falsifiability:**
- Detection at >2σ ΔN_eff: H1b CONFIRMED
- Null at <0.5σ: **H1b falsified ~3σ**

### §1.3 — Project 8 (CRES, MIT+UW 2030+) — m_β + sterile kink

**Status:** Phase II 2025+; final 2030+

**H1b prediction:**
- m_β sterile contribution: 0.032 eV
- Spectral kink at E_endpoint - 2 eV
- Combined m_β: 0.051 eV (active 0.04 + sterile 0.032)

**Falsifiability:**
- Detection of sterile kink + Δm_β > 10 meV: **H1b CONFIRMED (3σ)**
- Null: H1b sterile mass refinement needed

### §1.4 — Euclid (ESA 2023-2030 survey) — M-T_X precision

**Status:** active survey 2024+

**H1b prediction:**
- f_sterile_ν = 4.73 ± 0.71 (cluster mean; CV 15%)
- M-T_X residual scatter ≤ 0.05 dex (improved nad current 0.13)

**Falsifiability:**
- Global f_sterile constant within CV < 5%: **H1b CONFIRMED (5σ)**
- f_sterile cluster-mass-dependent: H1b sektor refinement

### §1.5 — Athena (ESA L-mission 2035+) — ICM mapping

**Status:** approved L-mission; launch 2035

**H1b prediction:**
- Sterile ν NFW spatial profile r_s ≈ 300 kpc
- Spatial resolution 124 elements within r_s for Coma
- M-T_X scatter < 0.05 dex post-Athena

**Falsifiability:**
- NFW spatial profile match: H1b CONFIRMED quantitatively
- Cored profile (warm DM signature): H1b warm DM correction needed

### §1.6 — PROSPECT-II + JSNS² + Daya Bay-II (2024-2027) combined

**Status:** all running

Combined factor 5-10× improvement nad current sterile ν bounds at Δm² ~ 4 eV².

## §2 — Deferred precision items

### §2.1 — Larger cluster sample (R1 partial)

**Status:** 🟡 deferred — 10 clusters analyzed; Euclid future expansion

**Phase 2+** future:
- Stacking analysis of ~10⁵ Euclid clusters
- Mass-redshift evolution z (z = 0 to z = 1)
- Sub-cluster scale (groups, infall regions)

### §2.2 — Full hydrodynamic N-body simulation (R2 partial)

**Status:** 🟡 deferred — analytical analysis only Phase 1-4

**Future cycle candidate:** `op-cluster-N-body-TGP-sterile-nu-simulation`
(estymata: ~4-6 sesji; requires external simulation code GADGET/AREPO)

### §2.3 — Cluster mass function dN/dM evolution

**Status:** 🟡 deferred — Phase 1-4 focus on individual clusters

**Future work:** TGP + sterile ν dN/dM(z) prediction; cross-comparison z observed
cluster mass function evolution (Vikhlinin 2009, Allen et al. 2011, eROSITA 2024+).

### §2.4 — Warm DM free-streaming corrections

**Status:** 🟡 deferred — Phase 2 assumed NFW; full warm DM modification needed

**Sterile ν m_ν = 2 eV** jest marginally warm (λ_fs ~ 100 kpc - 1 Mpc):
- Cluster cores: cored profile (warm DM signature)
- Cluster outskirts: NFW-like
- Phase 2+ refinement: cored NFW hybrid profile

## §3 — Cross-cycle propagation tasks

### §3.1 — L01 NEEDS (cluster scope explicit)

| Item | Status |
|---|---|
| L01 README: cross-link cluster cycle CLOSED | ⬜ next session |
| L01 NEEDS: cluster scope (separate od L01 N1-N5) noted | ⬜ next session |

### §3.2 — N3 SPARC cycle (galactic baseline preserved)

| Item | Status |
|---|---|
| N3 Phase_FINAL_close: cluster cycle CLOSED cross-link | ⬜ next session |
| N3 FINDINGS: cluster scope explicit (preserved galactic) | ⬜ next session |

### §3.3 — Q2 cycle (sterile ν compatibility)

| Item | Status |
|---|---|
| Q2 NEEDS: sterile ν matter sektor (NIE substrate vacuum) noted | ⬜ next session |

### §3.4 — PREDICTIONS_REGISTRY (5 new M911-cluster-* entries)

| Item | Status |
|---|---|
| M911-cluster-mass-deficit-H1b | ⬜ next session |
| M911-cluster-sterile-nu-CMB-S4 | ⬜ next session |
| M911-cluster-sterile-nu-Project-8 | ⬜ next session |
| M911-cluster-sterile-nu-JSNS2 | ⬜ next session |
| M911-cluster-sterile-nu-Euclid-Athena | ⬜ next session |

## §4 — Future cycles (deferred)

### §4.1 — op-cluster-N-body-TGP-sterile-nu-simulation

**Scope:** Full hydrodynamic N-body simulation w TGP + sterile ν framework
(estymata: ~4-6 sesji; requires external simulation code).

### §4.2 — op-warm-DM-cosmological-imprint

**Scope:** Sterile ν 2 eV cosmological imprint w P(k), Ly-α forest, dwarf
galaxy abundance (cross-validation z Lyman-α + 21cm).

## §5 — Cross-references

- [[./README.md]]
- [[./FINDINGS.md]]
- [[./Phase_FINAL_close.md]]
- [[../op-L01-N3-SPARC-rho-consistency-2026-05-11/]]
- [[../../PREDICTIONS_REGISTRY.md]] (M911-cluster-* entries planned)

---

**Cycle CLOSED konstruktywnie 2026-05-11.** Otwarte items są future empirical
tests + deferred precision + cross-cycle propagation tasks.
