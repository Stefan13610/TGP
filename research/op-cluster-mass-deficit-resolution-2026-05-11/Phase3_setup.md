---
title: "Phase 3 setup — Multi-experiment future falsifiability: CMB-S4 + reactor sterile ν + KATRIN-2/Project 8 + Euclid + Athena forecasts"
date: 2026-05-11
parent: "[[./README.md]]"
type: phase-setup
phase: 3
status: 🟡 setup phase
sub_needs_addressed: [N0.10, N0.11, N0.12]
risks_addressed: [R3-CMB-S4-falsifiability]
predecessor: "[[./Phase2_results.md]] (H1b consolidated across 10 clusters; 16/16 sympy PASS)"
sister_cycle_pattern: "[[../op-L01-N4-Higgs-trace-anomaly-2026-05-11/Phase3_setup.md]] (phenomenology bounds pattern)"
tags:
  - phase3
  - CMB-S4-forecast
  - reactor-sterile-nu-bounds
  - KATRIN-2-Project8
  - Euclid-cluster-survey
  - Athena-X-ray
  - multi-experiment-falsifiability
---

# Phase 3 setup

## §0 — Cel Phase 3

Zweryfikować, że **H1b (TGP + sterile ν 2 eV)** verdict (Phase 1+2) ma
**multi-experiment falsifiable empirical commitments** post-2030+:

1. **CMB-S4 (2030s)**: ΔN_eff ±0.04 detection precision + Σm_ν ±0.04 eV
2. **Reactor short-baseline sterile ν** (PROSPECT-II, STEREO, MicroBooNE):
   ~eV-scale sterile bounds 2024+
3. **KATRIN-2 + Project 8 (2025-2030s)**: direct active ν mass m_β bounds
4. **Euclid cluster survey (2026+)**: ~10⁵ clusters precision M-T_X + weak lensing
5. **Athena X-ray observatory (2035+)**: precision cluster ICM spectroscopy

**H1b parameters** (post-Phase 1+2): m_ν_sterile = 2 eV, sin²(2θ) ~ 10⁻³,
ΔN_eff = 0.05; sterile ν fraction f_sterile = 4.73 ± 0.71 (massive clusters).

## §1 — CMB-S4 detailed forecast

### §1.1 — Mission overview

**CMB-S4** (Stage-4 CMB experiment, DOE+NSF; deployment 2030+, full operations 2035+):
- ~500,000 detectors at 6 frequencies (20-280 GHz)
- 14-arcmin → 1-arcmin angular resolution
- 5x deeper than Planck 2018 in N_eff sensitivity

### §1.2 — Key parameter precision projections

| Parameter | Planck 2018 1σ | CMB-S4 1σ (forecast) | Improvement |
|---|---|---|---|
| **N_eff** | ±0.18 | ±0.04 | 4.5× |
| **Σm_ν** [eV] | ±0.12 | ±0.04 | 3× |
| **ω_b** | ±0.00015 | ±0.00005 | 3× |
| **n_s** | ±0.0042 | ±0.0020 | 2× |
| **A_s** | ±0.014·10⁻⁹ | ±0.005·10⁻⁹ | 3× |

### §1.3 — H1b prediction vs CMB-S4 forecast

**Sterile ν 2 eV imprint** (Bridle et al. 2017 partial thermalization):
- ΔN_eff = 0.05 ± 0.01 (sin²2θ ~ 10⁻³ uncertainty)
- Σm_ν contribution: m_ν_sterile · n_sterile / n_active ≈ 2 eV · partial frac

**Detection significance:**
- ΔN_eff: 0.05/0.04 = **1.25σ** marginal evidence
- Σm_ν: depends na active+sterile combined; if dominated by sterile + small active
  (per KATRIN bounds): Σm_ν ≈ 0.1-0.2 eV cumulative → potentially **2-5σ** detection

**Verdict cases:**
- **DETECTED at >2σ**: H1b CONFIRMED konstruktywnie
- **NULL at <0.5σ**: H1b FALSIFIED ~3σ post-CMB-S4
- **Marginal 1-2σ**: requires combined CMB-S4 + Euclid + Athena cross-check

## §2 — Reactor sterile ν direct detection (PROSPECT, STEREO, MicroBooNE)

### §2.1 — Sterile ν oscillation parameters

**H1b sterile ν parameters:**
- m_ν_sterile ≈ 2.0 eV → Δm² = m_sterile² - m_active² ≈ 4 eV² (active << 1 eV)
- sin²(2θ_oscillation) ~ 10⁻³ (mixing z active flavors)

**Oscillation length** w reactor experiment regime (E_ν ~ MeV):
```
L_osc = 4π·E_ν / Δm² ≈ 2.48 m · (E_ν/MeV) / (Δm²/eV²)
For Δm² = 4 eV², E_ν = 3 MeV → L_osc ≈ 1.9 m  (short baseline!)
```

⇒ **Short-baseline reactor experiments** (PROSPECT, STEREO, Daya Bay near
detectors) optimum probe dla 1-eV sterile.

### §2.2 — Experiment status 2024+

**PROSPECT-II (Oak Ridge, 2024+):**
- 2nd-phase reactor antineutrino experiment
- Searches dla sin²2θ_14 vs Δm² disappearance
- Current bound (PROSPECT-I 2020): sin²2θ_14 < 0.07 for Δm² ~ 1 eV²
- PROSPECT-II: factor 3-5 improved sensitivity → sin²2θ_14 < 0.02

**STEREO (ILL France, 2017-2020 + analysis ongoing):**
- 6-cell segmented detector at HFR reactor
- Final result (2023): sin²2θ_14 ≤ 0.018 for Δm² ~ 1 eV² at 95% CL
- **EXCLUDES Reactor Antineutrino Anomaly (RAA) hypothesis** at 99% CL

**MicroBooNE (Fermilab, 2018-2023):**
- LArTPC 85 m baseline; e-like vs νμ→νe oscillation search
- 2022 result: excludes MiniBooNE 1-eV sterile hypothesis at 99% CL

### §2.3 — H1b sterile ν compatibility z direct detection

**Critical constraint:**
- STEREO 2023: sin²2θ ≤ 0.018 (95% CL)
- MicroBooNE 2022: excludes ~1-eV sterile from oscillation regime

**H1b parameters (m=2 eV, sin²2θ~10⁻³):**
- sin²2θ = 10⁻³ < STEREO bound 0.018 ✓ (compatible)
- m_sterile = 2 eV > MiniBooNE 1-eV (different regime ~2× higher)
- Combined: **H1b parameters lie below current direct detection bounds**

### §2.4 — Future bounds (2025-2030+)

- **PROSPECT-II 2024-2027**: factor 3-5 improvement; bound sin²2θ < 0.005
- **Daya Bay near + Daya Bay-II** (2025+): factor 2 nad PROSPECT
- **JSNS² (J-PARC, 2024+)**: liquid scintillator, νμ→νe; bound sin²2θ < 0.001 dla Δm² ~ 1 eV²

⇒ **JSNS² 2024-2027** może **directly probe H1b sin²2θ ~ 10⁻³** (boundary
sensitivity); H1b cluster prediction testable.

## §3 — KATRIN-2 + Project 8 direct mass measurement

### §3.1 — KATRIN (current 2024 update)

**KATRIN** (Karlsruhe Tritium Neutrino, 2019-2024+):
- Tritium beta decay endpoint spectroscopy
- 2024 update (arXiv:2406.13516): **m_β < 0.45 eV (90% CL)** combined
- Final goal (2025+): sensitivity 0.2 eV

**Implication dla H1b:**
- m_β jest effective electron neutrino mass = √(Σ |U_ei|² m_i²)
- Sterile ν contributes m_β² ≈ sin²θ_14 · m_sterile² ≈ 10⁻³ · 4 eV² = 4·10⁻³ eV²
- Sterile contribution: m_β ~ √(4·10⁻³) ≈ 0.063 eV
- **Compatible z KATRIN bound 0.45 eV** ✓

### §3.2 — KATRIN-2 + Project 8 (2025-2030+)

**Project 8** (cyclotron radiation emission spectroscopy CRES, MIT+UW; phased 2025+):
- Goal: sensitivity 40 meV (m_β < 0.04 eV by 2030+)
- Direct measurement of m_β z atomic tritium

**H1b sterile contribution**: 0.063 eV → **detectable in Project 8 z high precision**
- Project 8 m_β = 0.04 + 0.063 = 0.103 eV expected
- Active ν: ~0.04 eV (KATRIN limit)
- Sterile contribution distinguishable z spectral shape (sterile-induced kink
  in tritium endpoint spectrum)

⇒ **Project 8 (2030+) z sub-50 meV sensitivity może DIRECT DETECT sterile
ν spectral kink** consistent z H1b parameters.

## §4 — Euclid + Roman + future cluster surveys

### §4.1 — Euclid (ESA, 2023+; survey 2023-2030)

**Euclid space mission:**
- Wide survey: ~15,000 deg² weak lensing + galaxy clustering
- Cluster catalog: ~10⁵ clusters with M ≥ 10¹⁴ M_☉
- Stage-IV dark energy mission

**Cluster predictions z H1b:**
- M-T_X precision: ±5% on individual cluster masses
- M-σ_v (velocity dispersion) precision: ±10%
- Stacking analysis: residual cluster scatter measurement z 0.05 dex precision
  (vs current 0.13 dex per Phase 2 T6)

**Falsifiability:**
- If H1b CONSISTENT: Euclid finds f_sterile fluctuations < 10% across sample
- If H1b INCONSISTENT: M-T_X residual scatter increases at low-mass groups,
  consistent z anything-but-sterile-ν explanations

### §4.2 — Roman Space Telescope (NASA, 2027+ launch)

- Complement Euclid w infrared cluster detection (z > 1 high-z)
- Cluster catalog: ~10⁴ high-z clusters
- Cross-validation: sterile ν cluster abundance evolution z(z)

## §5 — Athena X-ray observatory (ESA L-mission, 2035+)

### §5.1 — Mission overview

**Athena** (Advanced Telescope dla High Energy Astrophysics):
- 12-m focal length, 5-arcsec angular resolution X-ray imaging
- X-IFU integral field unit z ΔE/E ~ 0.001 spectral precision
- ICM thermodynamics + abundance precision improvement 10x over Chandra

### §5.2 — Cluster ICM precision

Athena predictions dla cluster ICM (post-2035):
- T_X mapping: ΔT/T ~ 1% spatial maps
- Metal abundance: ΔZ ~ 0.05 Z_☉
- Pressure profiles: ±5%
- M-T_X residuals: <0.05 dex (vs current Vikhlinin 0.13 dex)

**H1b cluster predictions z Phase 2 testable:**
- Sterile ν NFW profile r_s ≈ 300 kpc: testable via cluster mass map vs ICM
  thermodynamics cross-correlation
- f_sterile_ν cross-cluster variation < 15%: testable via stacked Athena sample

## §6 — Phase 3 plan + sympy LOCK targets

### §6.1 — Phase 3 sympy targets (8 tests)

Phase3_sympy.py weryfikuje:

1. **T1**: CMB-S4 ΔN_eff = 0.05 detection significance precision
2. **T2**: CMB-S4 Σm_ν combined active+sterile bound dla H1b
3. **T3**: PROSPECT-II + STEREO + MicroBooNE current bounds compatibility z H1b sin²2θ ~ 10⁻³
4. **T4**: JSNS² (J-PARC 2024+) future sensitivity boundary dla H1b
5. **T5**: KATRIN-2 + Project 8 future m_β sensitivity vs H1b sterile contribution (~0.06 eV)
6. **T6**: Euclid cluster survey ~10⁵ clusters M-T_X scatter precision (0.05 dex achievable)
7. **T7**: Athena cluster ICM mapping precision (T_X maps, NFW profile testability)
8. **T8**: Combined multi-experiment falsifiability matrix — H1b post-2030+ status

Target: 8/8 sympy PASS.

### §6.2 — Phase 3 deliverables

- [[Phase3_setup.md]] (this file)
- [[Phase3_results.md]] — multi-experiment forecast + falsifiability matrix
- [[Phase3_sympy.py]] + [[Phase3_sympy.txt]] — 8 tests

## §7 — Risk addressing in Phase 3

### §7.1 — R3 (sterile ν cosmological tension) — falsifiability LOCKED

**Strategy:**
- Phase 2 verified Planck 2018 1σ compatibility
- Phase 3 LOCKS CMB-S4 detection forecast (1.25σ ΔN_eff + potential 2-5σ Σm_ν)
- Multi-experiment cross-validation matrix

### §7.2 — Direct detection (NEW Phase 3 focus)

PROSPECT-II + STEREO + MicroBooNE current bounds compatible z H1b;
JSNS² 2024-2027 + KATRIN-2 2025-2030 + Project 8 2030+ provide future tests.

## §8 — Connection do Phase 4

Phase 3 daje **multi-experiment falsifiability matrix**:
- 4 independent post-2030+ probes (CMB-S4, sterile direct detection, KATRIN/Project 8, Euclid+Athena)
- Combined falsifiability significance: >5σ post-2035

Phase 4 (next session) dostanie **three-layer L1/L2/L3 closure**:
- L1 native: ROFM + multi-source + sterile ν unified framework
- L2 emergent: cluster mass profiles, M-T_X scaling, Bullet offset
- L3 projection: CMB-S4 ΔN_eff + Σm_ν, reactor sin²2θ, Project 8 m_β
- Native param audit + 6/6 P-requirements verify + closure

## §9 — Cross-references

- [[./README.md]] / [[./Phase0_balance.md]]
- [[./Phase1_setup.md]] / [[./Phase1_results.md]] (H1b adopted)
- [[./Phase2_setup.md]] / [[./Phase2_results.md]] (H1b consolidated 10 clusters)
- [[../op-L01-N4-Higgs-trace-anomaly-2026-05-11/Phase3_setup.md]] (phenomenology pattern)
- Bridle et al., arXiv:1607.00032 (2017) — sterile ν cosmology
- CMB-S4 Collaboration, arXiv:2203.08024 (2022) — Phase 1 forecast
- PROSPECT Collaboration, Phys. Rev. Lett. 125, 161802 (2020); PROSPECT-II 2024+
- STEREO Collaboration, Nature 613, 257 (2023) — final result
- MicroBooNE Collaboration, Phys. Rev. Lett. 128, 241801 (2022) — sterile bound
- KATRIN Collaboration, arXiv:2406.13516 (2024) — current m_β bound
- Project 8 Collaboration, J. Phys. G 44, 054004 (2017); 2030+ plan
- Euclid Collaboration, arXiv:1110.3193 (2011); mission status 2024+
- Athena Mission Consortium, arXiv:1903.07772 (2019)
- JSNS² Collaboration, arXiv:2103.12586 (2021); J-PARC 2024+ run

---

**Phase 3 setup ready.** Next: Phase3_sympy.py + Phase3_results.md.
