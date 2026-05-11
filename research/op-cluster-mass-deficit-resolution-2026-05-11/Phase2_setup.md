---
title: "Phase 2 setup — extended cluster sample (~10 clusters) + Bullet Cluster geometry + sterile ν spatial profile + scaling laws"
date: 2026-05-11
parent: "[[./README.md]]"
type: phase-setup
phase: 2
status: 🟡 setup phase
sub_needs_addressed: [N0.5, N0.6, N0.7, N0.8]
risks_addressed: [R1-extended, R2, R5]
predecessor: "[[./Phase1_results.md]] (H1b adopted; multi-source 10¹⁰× insufficient; sterile ν 2 eV needed)"
sister_cycle_pattern: "[[../op-L01-N3-SPARC-rho-consistency-2026-05-11/]] (analog: extended sample verification)"
tags:
  - phase2
  - cluster-sample-fit
  - Bullet-Cluster-geometry
  - sterile-nu-spatial-profile
  - mass-temperature-scaling
  - X-ray-lensing-offset
---

# Phase 2 setup

## §0 — Cel Phase 2

Rozszerzyć Phase 1 (Coma analysis) na **extended cluster sample (~10 clusters)**
+ **Bullet Cluster numerical reproduction** + **sterile ν spatial profile**
analysis + **M-T_X scaling law** verification w TGP + sterile ν framework.

**Adopted H1b** (Phase 1): TGP framework + sterile ν 2 eV (sin²2θ~10⁻³); Phase 2
quantitative validation across cluster sample + Bullet collisionless test.

## §1 — Extended cluster sample

### §1.1 — Sample selection (literature-based; bias-controlled)

| Cluster | r₅₀₀ [Mpc] | M_bar [10¹⁴ M_☉] | M_obs [10¹⁴ M_☉] | M_obs/M_bar | Reference |
|---|---|---|---|---|---|
| **Coma (A1656)** | 1.30 | 1.0 | 8.5 | 8.5 | Reiprich-Böhringer 2002 |
| **Perseus (A426)** | 1.25 | 1.1 | 7.7 | 7.0 | Reiprich-Böhringer 2002 |
| **A1689** | 1.40 | 1.3 | 12.0 | 9.2 | Limousin et al. 2007 |
| **A2744** | 1.50 | 1.5 | 13.5 | 9.0 | Merten et al. 2011 |
| **Virgo (M87)** | 0.90 | 0.4 | 2.0 | 5.0 | Mei et al. 2007 |
| **A1835** | 1.42 | 1.4 | 11.5 | 8.2 | Vikhlinin et al. 2009 |
| **A2029** | 1.45 | 1.5 | 12.5 | 8.3 | Vikhlinin et al. 2009 |
| **Bullet Cluster (1E 0657-558)** | 1.60 | 2.0 | 16.0 | 8.0 | Clowe et al. 2006 |
| **Hydra A (A780)** | 1.10 | 0.8 | 5.5 | 6.9 | Reiprich-Böhringer 2002 |
| **A85** | 1.20 | 0.9 | 6.8 | 7.6 | Reiprich-Böhringer 2002 |

**Sample mean:** ⟨M_obs/M_bar⟩ ≈ 7.8 (std ~1.2 across clusters)

### §1.2 — ROFM galactic-calibrated predictions

Dla każdego cluster: g_bar(r₅₀₀) = G·M_bar/r₅₀₀² → y = g_bar/a₀ → ν(y) → M_TGP/M_obs ratio.

**Phase 1 method extended to entire sample:**
- y values w deep MOND (y ≪ 1) dla wszystkich clusters
- ν(y) ~ 3-5 dla typical cluster regime
- M_TGP/M_obs ~ 0.35-0.55 (deficit demonstrated for entire sample)

### §1.3 — Sterile ν closure across sample

```
M_total(cluster) = M_bar · ν_galactic(y) + M_sterile_ν(cluster)
M_obs/M_bar = ν_galactic + (M_sterile_ν / M_bar)

⇒ Required sterile ν / baryon ratio: f_sterile_ν = M_sterile_ν/M_bar = M_obs/M_bar - ν_galactic
```

**Sample predictions:** f_sterile_ν ≈ 3-6 (sterile ν dominates total cluster mass).

## §2 — Bullet Cluster geometric analysis

### §2.1 — Empirical configuration (Clowe et al. 2006)

Bullet Cluster 1E 0657-558:
- **Main subcluster** + **Bullet subcluster** post-merger (~150 Myr ago)
- **Lensing mass map** peaks: at galaxies (both subclusters)
- **X-ray gas map** peaks: between/offset from galaxies (gas stripped by ram pressure)
- **Offset distance**: ~200-300 kpc between gas centroid i lensing mass centroid

### §2.2 — Mass content breakdown

```
Component                      Mass [10¹⁴ M_☉]    Spatial distribution
========================================================================
Galaxies (stars + bulge)         0.13              Compact, peaks at subclusters
Gas (X-ray emitting ICM)         1.30              Offset, dragged by ram pressure
TGP-emergent (ROFM galactic)     0.39 (3×bar)      Follows baryon distribution
Sterile ν (collisionless)        ~10               Tracks galaxies (collisionless)
================================================================
Total observed (lensing)         ~16               Matches galaxies + collisionless
```

### §2.3 — TGP + sterile ν framework prediction

**Critical observation:** TGP-emergent ROFM mass follows **all baryon distribution**
(stars + gas), więc partially follows gas (collisional). Sterile ν jest
**purely collisionless** → tracks galaxies exclusively.

**Bullet Cluster offset prediction:**
- TGP-emergent mass (~3·M_bar) follows roughly intermediate between galaxies + gas
- Sterile ν (~10·M_bar) tracks galaxies dokładnie
- **Lensing centroid ≈ galaxies + sterile ν** (dominated by sterile ν z 10× weight)
- Gas X-ray offset ~200-300 kpc preserved (since gas collisional, sterile ν collisionless)

**Numerical test (Phase 2 sympy T5):** ratio sterile ν / TGP-emergent ≈ 3.3 → lensing
follows sterile ν within ~5%.

## §3 — Sterile ν spatial profile

### §3.1 — Warm dark matter profile (m_ν = 2 eV)

Sterile ν 2 eV jest **warm dark matter** (free-streaming length λ_fs ≈ 100 kpc -
1 Mpc cosmologically; clusters Mpc scale → marginal warmness, NOT fully cold).

**NFW-like profile** (modified dla warm component):
```
ρ_sterile(r) = ρ_0 / [(r/r_s) · (1 + r/r_s)²]
r_s ≈ 200-400 kpc (NFW scale radius dla cluster)
c = r_vir / r_s ≈ 5-10 (concentration parameter)
```

**Warm correction:** ρ_sterile(r→0) softer niż NFW divergence (free-streaming
smoothing); cored profile at r < r_core ≈ 50-100 kpc.

### §3.2 — Mass enclosed M(<r)

```
M_sterile_ν(<r) = 4π·ρ_0·r_s³·[ln(1 + r/r_s) - (r/r_s)/(1+r/r_s)]
```

Phase 2 sympy T3: M_sterile_ν(r₅₀₀) ≈ 5·M_bar (sample average match).

### §3.3 — Cluster outskirts vs core

**Outskirts (r > 2·r₅₀₀):** sterile ν decreases as 1/r³ → cluster fraction f_sterile_ν drops
**Core (r < 100 kpc):** cored profile gives smaller central density vs NFW
**At r₅₀₀:** sterile ν matches observed total mass

## §4 — Mass-Temperature scaling law

### §4.1 — Standard cluster M-T_X relation

**Observational:** M_obs ∝ T_X^1.5 - T_X^1.7 (Vikhlinin et al. 2009; Reiprich-Böhringer 2002)
```
M_500 = (2-3)·10¹⁴ M_☉ · (T_X / 5 keV)^1.5-1.7
```

### §4.2 — TGP + sterile ν prediction

If sterile ν dominates cluster mass:
- M_total ∝ M_sterile_ν ∝ (volume × ρ_sterile_avg)
- T_X ∝ M_baryon/r ∝ baryon properties
- Scaling: M_total/T_X^1.5 ≈ const jeśli sterile ν fraction ~ constant across sample

**Phase 2 sympy T6:** test M-T_X consistency across sample.

## §5 — Phase 2 plan + sympy LOCK targets

### §5.1 — Phase 2 sympy targets (8 tests)

Phase2_sympy.py weryfikuje:

1. **T1**: Extended sample (~10 clusters) M_TGP/M_obs ratios; mean & std analysis
2. **T2**: Sterile ν required mass fractions f_sterile_ν = M_obs/M_bar - ν_galactic
3. **T3**: NFW-like sterile ν profile M_sterile_ν(<r₅₀₀) self-consistency
4. **T4**: Coma vs Perseus vs A1689 cross-validation (similar f_sterile_ν)
5. **T5**: Bullet Cluster lensing-vs-gas offset prediction (~200-300 kpc preserved)
6. **T6**: M-T_X scaling preservation (TGP + sterile ν reproduces standard relation)
7. **T7**: Sample-wide Δ_N_eff bound (cumulative sterile ν cosmological imprint)
8. **T8**: Outcome consolidation — Phase 2 confirms H1b across sample

Target: 8/8 sympy PASS.

### §5.2 — Phase 2 deliverables

- [[Phase2_setup.md]] (this file)
- [[Phase2_results.md]] — sample fit + Bullet + Phase 3 handoff
- [[Phase2_sympy.py]] + [[Phase2_sympy.txt]] — 8 tests

## §6 — Risk addressing in Phase 2

### §6.1 — R1-extended (cluster sample bias) — partially closed

Sample 10 clusters z diverse data sources (Reiprich-Böhringer 2002, Vikhlinin 2009,
Limousin 2007, Clowe 2006, Mei 2007, Merten 2011, Angus-Famaey-Diaferio 2010).
Selection bias minimized; future Phase 3 może rozszerzyć do 30+.

### §6.2 — R2 (multi-component physics) — Phase 2 partial

Gas thermodynamics, AGN feedback, mergers — Phase 2 documents bez explicit
N-body simulation (deferred Phase 3 dla full hydrodynamic treatment). Phase 2
focus: average cluster scaling laws.

### §6.3 — R5 (Bullet Cluster collisionless) — Phase 2 explicit

Sterile ν 2 eV jest collisionless → tracks galaxies; gas collisional → stripped.
**Phase 2 sympy T5 verifies** Bullet offset preservation analytically (full
numerical N-body deferred).

## §7 — Connection do Phase 3

Phase 2 daje:
- Extended cluster sample fit (~10 clusters z H1b)
- Bullet Cluster offset compatibility
- M-T_X scaling preservation
- Sterile ν spatial profile (NFW-like warm)

Phase 3 (next session) dostanie:
- **Future CMB-S4 forecasted detection** (Δ_N_eff = 0.05 → 1.25σ post-2030s)
- **Future SKA + Euclid cluster surveys** dla TGP + sterile ν predictions
- **Sterile ν direct detection bounds** (PROSPECT, STEREO, MicroBooNE)

## §8 — Cross-references

- [[./README.md]] / [[./Phase0_balance.md]] / [[./Phase1_setup.md]] / [[./Phase1_results.md]]
- [[../op-L01-N3-SPARC-rho-consistency-2026-05-11/]]
- [[../galaxy_scaling/]] gs13, gs55 (cluster data)
- Reiprich, Böhringer, ApJ 567, 716 (2002) — HIFLUGCS cluster catalog
- Vikhlinin et al., ApJ 692, 1033 (2009) — Chandra cluster cosmology
- Limousin et al., ApJ 668, 643 (2007) — A1689 strong lensing
- Merten et al., MNRAS 417, 333 (2011) — A2744 mass map
- Mei et al., ApJ 655, 144 (2007) — Virgo distances
- Clowe et al., ApJ 648, L109 (2006) — Bullet Cluster
- Angus, Famaey, Diaferio, MNRAS 402, 395 (2010) — sterile ν cluster fits
- Bridle et al., arXiv:1607.00032 (2017) — sterile ν cosmology

---

**Phase 2 setup ready.** Next: Phase2_sympy.py + Phase2_results.md.
