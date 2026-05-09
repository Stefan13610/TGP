---
title: "Phase 1 — Waveform model setup (TaylorF2 + ppE deviation)"
date: 2026-05-07
parent: "[[README.md]]"
type: phase1-results
tgp_owner: research/op-LIGO-3G-deviation
tags:
  - phase1
  - waveform
  - TaylorF2
  - ppE
  - SPA
related:
  - "[[README.md]]"
  - "[[Phase0_balance.md]]"
  - "[[scripts/phase2_fisher_forecast.py]]"
  - "[[../op-ppE-mapping/Phase1_results.md]]"
---

# Phase 1 — Waveform model setup

## §1 — Decyzja model: TaylorF2 + ppE deviation

Z [[../../audyt/T01_LIGO3G_falsifier/CYCLE_KICKOFF_op-LIGO-3G-deviation.md]] §2.1:

**Wybrany model:** **TaylorF2 GR baseline + ppE single-coefficient
deviation** (perturbative inspiral).

**Argumenty:**
- M9.1'' two-body Lagrangian (z `op-ppE-mapping/` Phase 1) jest
  perturbative w (5/6) U³ — full NR-TGP nie istnieje.
- Dla inspiral signal (f < f_ISCO) U ≪ 1 i perturbative TaylorF2
  jest wystarczające.
- TaylorF2 jest standard dla LIGO O3/O4 ToGR analyses (LSC papers).

**Cutoff:** f_max = f_ISCO = 4400/(M_total/M_⊙) Hz; wykluczać
merger/ringdown gdzie U → 1/2.

## §2 — Phase formula (TaylorF2 + ppE)

```
Ψ(f) = 2π f t_c − Φ_c − π/4
     + (3/(128 η)) · v⁻⁵ · [1 + α_2 v² + α_3 v³ + α_4 v⁴ + α_5 v⁵ (1+ln v) + α_6 v⁶ + ...]
     + β_ppE · v^b_ppE                                                       ← TGP deviation
```

gdzie:
- v = (π M f)^(1/3) (PN parameter, M = total mass in geometric units)
- α_n: PN coefficients (Buonanno et al. 2009 TaylorF2)
- β_ppE^TGP^(b=-1) = -5/64 ≈ -7.81·10⁻² (LOCKED z [[../op-ppE-mapping/Phase1_results.md]])

## §3 — Amplitude (sky-averaged)

```
|h(f)| = (1/d_L) · sqrt(5/24) · (G·M_chirp)^(5/6) · (πc²)^(-2/3) · f^(-7/6) · sqrt(4/5)
```

(sky-averaging factor 4/5 dla average inclination i orientation).

## §4 — ASD curves (analytical fits do design sensitivities)

### LIGO-O5 (A+ design)
Reference: LIGO-T1800042. Fit z dominant components:
- seismic ~ f^(-4)
- thermal Brownian ~ flat
- shot noise ~ f^(2)

```python
ASD(f) = sqrt( (3e-23)² · [0.07·(f/215)^(-4.14) + 5e-3/(1+(f/215)²) + 0.5·(1+(f/215)²)] )
```

### Einstein Telescope (ET-D)
Reference: Hild 2010 / Maggiore 2020. Cryogenic + xylophone.

```python
ASD(f) = sqrt( (5e-25)² · [2.0·(f/200)^(-4.5) + 0.05/(1+(f/200)²) + 0.4·(1+(f/200)^1.6)] )
       + low_freq_extension below 10 Hz
```

### Cosmic Explorer (CE)
Reference: Reitze 2019, Hall+Evans 2019. 40km arms.

```python
ASD(f) = sqrt( (3e-25)² · [1.0·(f/100)^(-4.5) + 0.02/(1+(f/100)²) + 0.3·(1+(f/100)²)] )
       + low_freq_extension below 7 Hz
```

**Caveat:** te są analityczne fits OOM-precision (factor ~3-5x off
od production-grade noise curves). Dla precise SNR estimates produkcyjne
analyses (LSC pipeline) wymagają full ASD interpolation. **Phase 2
output będzie referencjonalny — należy mapować na literature
SNR (Maggiore 2020 ET-D, Reitze 2019 CE).**

## §5 — Output → Phase 2

Skrypt: [[scripts/phase2_fisher_forecast.py]]

Outputs:
- SNR_inner_product(asd, M_c, η, d_L, f_lo, f_hi)
- σ_β = 1/sqrt(F_ββ) z degeneracy_factor = 5 (Yagi-Yunes 2016 for
  realistic Fisher z M_chirp / η / χ_eff covariance)
- β_5σ = 5 · σ_β

Wyniki w [[Phase2_results.md]].
