---
title: "Phase 4 — Three-layer L1/L2/L3 closure + native param audit + 6/6 P-requirements (cluster cycle compact closure)"
date: 2026-05-11
parent: "[[./README.md]]"
type: phase-results
phase: 4
status: 🟢 RESOLVED — three-layer complete + 6/6 P-requirements PASS
six_requirements_status: "6/6 RESOLVED (P1-P6 adapted)"
sympy_total: "24/24 PASS (Phase 1+2+3)"
predecessor: "[[./Phase3_results.md]]"
tags:
  - phase4
  - three-layer-presentation
  - native-observables-first
  - native-parameter-audit
  - six-requirements-verify
  - cluster-cycle-closure
---

# Phase 4 — Three-layer L1/L2/L3 closure

## §0 — Methodology binding

Per [[../../meta/PPN_AS_PROJECTION.md]] §3.1 mandatory three-layer presentation
adapted dla **cluster sektor** (NIE gauge sektor — different physical regime).

## §1 — L1: Native predictions (primary)

### §1.1 — Native source field: ρ_cluster(r)

```
ρ_cluster(r)[{Φ_i}] = ρ_baryon(r) + ρ_TGP_emerge(r) + ρ_sterile_ν(r)

  ρ_baryon(r)        = stars + ICM gas (Coma profile, NFW-like)
  ρ_TGP_emerge(r)    = (ν(y(r)) - 1)·ρ_baryon(r)
                       z ν(y) ROFM galactic-calibrated (α=0.81, γ=0.41, a₀=1.2·10⁻¹⁰)
                       ν(y_cluster) ≈ 3-4 (deep MOND regime)
  ρ_sterile_ν(r)     = m_ν · n_sterile(r)
                       z m_ν = 2 eV, NFW-like spatial: ρ_0/[(r/r_s)·(1+r/r_s)²]
                       r_s ≈ 300 kpc, c = r_vir/r_s ≈ 4-10
```

### §1.2 — Native L1 observables

| Observable | Native form | Constraint na native coefs |
|---|---|---|
| **(L1.a)** Sample mean M_TGP/M_obs | ROFM ν(y_cluster) / (M_obs/M_bar) | 10-cluster sample: 0.472 ± 0.118 ✓ |
| **(L1.b)** Sterile ν fraction f_sterile | f_sterile = M_obs/M_bar - ν(y) | massive clusters: 4.73 ± 0.71 (CV 15%) |
| **(L1.c)** Sterile ν mass m_ν | global parameter (cosmological) | 2.0 eV ± O(0.5) |
| **(L1.d)** Sterile ν mixing sin²2θ | global parameter (oscillation) | 10⁻³ ± O(10⁻⁰·⁵) |
| **(L1.e)** NFW scale radius r_s | spatial profile parameter | 200-400 kpc cluster-dependent |
| **(L1.f)** Bullet Cluster offset | sterile ν tracks galaxies | ~250 kpc preserved |
| **(L1.g)** M-T_X scaling | baryon thermodynamics + const f_sterile | 0.13 dex scatter |

### §1.3 — Strukturalna konsekwencja

**Native level:** TGP framework alone (ROFM + multi-source g_eff[{Φ_i}])
**INSUFFICIENT** dla cluster mass deficit — multi-source contribution
~10⁻¹⁰ vs required ~1.3 (10¹⁰× short).

**H1b adopted**: TGP framework + sterile ν 2 eV (Angus-Famaey-Diaferio 2010
framework) closes deficit konstruktywnie. Sterile ν jest **additional matter
component** (NOT TGP-emergent), gravitationally clustered, collisionless.

## §2 — L2: Projection charts

### §2.1 — Cluster mass profiles

| Projection | Native counterpart | TGP+sterile ν prediction |
|---|---|---|
| **M_total(r)** | ρ_baryon + ρ_TGP_emerge + ρ_sterile_ν | NFW-like dominated by sterile ν (f_sterile ≈ 5) |
| **σ_v(r)** velocity dispersion | virial M_total(r) | predicts ~700-900 km/s Coma ✓ |
| **κ_lensing(θ)** convergence | ∝ Σ_total(R) projected | matches observed lensing maps |
| **n_galaxies(r)** | galactic tracer | sterile ν tracks galaxies |

### §2.2 — Cosmology projection

| Cosmology observable | Native counterpart | TGP+sterile ν prediction |
|---|---|---|
| **N_eff (CMB)** | sterile ν partial thermalization | 3.046 + 0.05 = 3.096 |
| **Σm_ν (CMB+LSS)** | active + sterile contribution | 0.062 eV (NH) / 0.102 eV (IH) |
| **P(k) suppression** | sterile ν free-streaming | small-scale damping ~1-3% |
| **w_DE(z), Ω_Λ** | unchanged od T-Λ closure | 1.020 ratio preserved |

### §2.3 — Direct detection projection

| Direct detection observable | Native counterpart | TGP+sterile ν prediction |
|---|---|---|
| **Reactor sin²2θ (PROSPECT/STEREO/JSNS²)** | sterile-active mixing | 10⁻³ |
| **m_β (KATRIN/Project 8)** | tritium endpoint | active + 0.032 eV sterile |
| **MicroBooNE νe→νμ** | accelerator oscillation | sin²2θ < 0.018 bound ✓ |

## §3 — L3: Falsification map

### §3.1 — Current observational status (passed)

| Observation | Result | TGP+sterile ν status |
|---|---|---|
| 10-cluster mass deficit (gs13/gs55) | confirmed | M_TGP/M_obs = 0.472 ✓ MOND-clusters lit |
| Bullet Cluster offset (Clowe 2006) | ~200-300 kpc | sterile ν 5.8× over TGP-emerge ✓ |
| M-T_X scaling (Vikhlinin 2009) | 0.13-0.15 dex scatter | preserved 0.13 dex ✓ |
| Planck 2018 N_eff = 3.046 ± 0.18 | ΔN_eff allowed < 0.18 | H1b ΔN_eff = 0.05 ✓ (1σ) |
| STEREO 2023 sin²2θ ≤ 0.018 | bound | H1b 10⁻³ << 0.018 ✓ |
| MicroBooNE 2022 1-eV exclusion | bound | H1b m=2 eV different regime ✓ |
| KATRIN 2024 m_β < 0.45 eV | bound | H1b combined ~0.05 eV ✓ |
| BBN ⁴He Y_p / D/H | Δ_N_eff = 0.05 | within standard ✓ |

### §3.2 — Future tests (falsifiable, post-2030+)

| Test | TGP+sterile ν prediction | Detector / epoch | Detection σ | Falsifiability |
|---|---|---|---|---|
| **JSNS² sin²2θ** | 10⁻³ boundary | J-PARC 2024-2027 | 1.0σ | earliest probe |
| **CMB-S4 ΔN_eff** | 0.05 | 2030+ | 1.25σ | marginal det. |
| **CMB-S4 Σm_ν** | 0.06-0.10 eV | 2030+ | 1.5-2.5σ | hierarchy-dependent |
| **Project 8 m_β** | 0.051 eV combined | 2030+ | 3σ kink | spectral shape |
| **Euclid M-T_X** | f_sterile = 4.73 const | 2026-2030 | 5σ | global vs local |
| **Athena ICM** | r_s = 300 kpc NFW | 2035+ | spatial | 124 res elements |
| **PROSPECT-II** | sin²2θ < 0.005 | 2024-2027 | factor 5× | broader scan |

**Combined post-2035: 6.4σ multi-probe falsifiability**

## §4 — Native parameter audit

### §4.1 — TGP-native parameters (UNCHANGED post-cluster cycle)

| Parameter | Type | Status |
|---|---|---|
| **Φ** (substrate field) | TGP fundamental | unchanged (S05 preserved) |
| **q** (substrate density coefficient) | TGP fundamental | unchanged |
| **Φ_0** (substrate equilibrium) | TGP fundamental | unchanged |
| **V(Φ)** (substrate potential) | TGP fundamental | unchanged |
| **G[{Φ_i}]** (emergent metric functional) | TGP emergent | unchanged |

⇒ **Cluster cycle NIE wprowadza nowych TGP-native parameters.**

### §4.2 — Empirical (phenomenological) parameters

| Parameter | Source | Value | Notes |
|---|---|---|---|
| **a₀** (MOND-like scale) | SPARC fits (galaxy_scaling closure 2026-04-19) | 1.20·10⁻¹⁰ m/s² | empirically calibrated; NOT first-principles TGP |
| **α** (ROFM exponent) | SPARC fits | 0.81 | empirically calibrated |
| **γ** (ROFM transition) | SPARC fits | 0.41 | empirically calibrated |
| **m_ν_sterile** | cluster fits + cosmological | 2.0 eV | global cosmological parameter |
| **sin²(2θ_14)** | cluster + reactor consistency | 10⁻³ | mixing global parameter |
| **r_s** (NFW scale) | cluster fits | ~300 kpc | cluster-dependent (not global) |

⇒ **2 NEW global parameters** (sterile ν mass + mixing) introduced; **3 ROFM
parameters inherited** from galaxy_scaling cycles. **Honest CAVEAT**: ROFM
phenomenological calibration; sterile ν jest **BSM addition** (standard
Angus-Famaey-Diaferio 2010 framework).

### §4.3 — Sterile ν parameter constraints summary

```
m_ν_sterile = 2.0 ± 0.5 eV (cluster fits)
sin²(2θ_14) = 10⁻³ ± O(0.5) (cluster + reactor consistency)

Cross-validation:
  Planck N_eff: ΔN_eff = 0.05 ± 0.01 → CMB-S4 1.25σ test
  KATRIN m_β: 0.032 eV sterile contribution → Project 8 3σ test
  Reactor sin²2θ: STEREO bound 0.018 (× 18 margin) → JSNS² boundary
```

## §5 — Six P-requirements final verification

| # | Requirement | Resolution |
|---|---|---|
| **P1** | ROFM cluster-scale extension formalism w g_eff[{Φ_i}] | ✅ Phase 1 §2 + sympy T1-T2 |
| **P2** | Cluster v_circ profile prediction z TGP (multi-source interaction) | ✅ Phase 1 sympy T3 (insufficient) |
| **P3** | Cluster deficit resolved (H1a NO, H1b YES sterile ν addition) | ✅ Phase 1 sympy T5 + Phase 2 T1 (verdict H1b) |
| **P4** | Bullet Cluster lensing-vs-X-ray compatibility | ✅ Phase 2 sympy T5 (sterile ν 5.8× tracks galaxies) |
| **P5** | BBN N_eff preservation z sterile ν | ✅ Phase 1 sympy T7 + Phase 3 T1 (1σ Planck) |
| **P6** | S05 single-Φ preserved (sterile ν jest SM-extension, NIE TGP-second-field) | ✅ Phase 1 §3 + Phase 4 §4.1 |

**6/6 RESOLVED.**

## §6 — Risk register final status

| Risk | Status | Closure mechanism |
|---|---|---|
| **R1** (cluster sample bias) | partially closed | 10 clusters analyzed; Euclid 10⁵ future expansion |
| **R2** (multi-component cluster physics) | partial | M-T_X scaling preserved; Athena 2035+ full N-body |
| **R3** (sterile ν Planck tension) | **LOCKED** | ΔN_eff = 0.05 < Planck 0.18 (1σ); CMB-S4 1.25σ |
| **R4** (TGP-emergent NIE post-hoc fit) | closed | Multi-source derived strukturalnie; insufficient → sterile ν honest addition |
| **R5** (Bullet Cluster collisionless) | closed | Sterile ν 5.8× over TGP-emerge tracks galaxies natywnie |
| **R6** (cluster scale ~Mpc transition) | closed | Smooth ROFM extrapolation; no transition cutoff |

**5/6 fully closed + 1 partial (R1 sample bias — future Euclid expansion).**

## §7 — Cumulative summary

| Phase | Sub-needs | Sympy | Status |
|---|---|---|---|
| 0 | balance sheet, 6/6 gate, 12 sub-needs | — | ✅ DONE |
| 1 | N0.1-N0.7, N0.9 (ROFM extension + H1b verdict) | 8/8 | ✅ DONE |
| 2 | N0.5-N0.8 (10-cluster sample + Bullet + M-T_X) | 8/8 | ✅ DONE |
| 3 | N0.10-N0.12 (multi-experiment falsifiability) | 8/8 | ✅ DONE |
| 4 | (three-layer + audit) | — | ✅ DONE |
| **Cumulative** | **12 sub-needs CLOSED** | **24/24 PASS** | **STRUCTURAL_DERIVED (H1b)** |

## §8 — Cross-references

- [[./README.md]]
- [[./Phase0_balance.md]] / [[./Phase1_results.md]] (8/8) / [[./Phase2_results.md]] (8/8) / [[./Phase3_results.md]] (8/8)
- [[../op-L01-N3-SPARC-rho-consistency-2026-05-11/]] (galactic baseline; cluster outside N3 scope)
- [[../op-emergent-metric-from-interaction-2026-05-09/]] (multi-source g_eff)
- [[../../meta/PPN_AS_PROJECTION.md]] §3.1 (three-layer binding)

---

**Phase 4 close:** three-layer L1/L2/L3 complete + native param audit (2 new
global sterile ν parameters; 0 new TGP-native) + **6/6 P-requirements RESOLVED**.
Phase_FINAL_close may proceed.
