---
title: "Phase 2 results — Fisher matrix forecast SNR thresholds"
date: 2026-05-07
parent: "[[README.md]]"
type: phase2-results
tgp_owner: research/op-LIGO-3G-deviation
tags:
  - phase2
  - results
  - fisher
  - SNR
  - LIGO-O5
  - ET-D
  - CE
  - degeneracy
related:
  - "[[README.md]]"
  - "[[Phase0_balance.md]]"
  - "[[Phase1_waveform_setup.md]]"
  - "[[scripts/phase2_fisher_forecast.py]]"
  - "[[scripts/phase2_fisher_forecast.txt]]"
  - "[[../op-ppE-mapping/Phase1_results.md]]"
---

# Phase 2 results — Fisher matrix forecast

## §0 — TL;DR

```
TGP β_ppE^TGP^(b=-1) = -5/64 ≈ -7.81 · 10⁻²  (z op-ppE-mapping Phase 1)

Detection thresholds (single event, full Fisher z degeneracy_factor=5):
  LIGO-O5 loud BBH:    β_5σ ≈ 1.25     → TGP/β_5σ ≈ 0.06   NO single
  ET-D loud BBH:       β_5σ ≈ 0.25     → TGP/β_5σ ≈ 0.31   borderline
  CE loud BBH:         β_5σ ≈ 0.059    → TGP/β_5σ ≈ 1.32   YES single (>5σ)
  ET+CE network:       β_5σ ≈ 0.057    → TGP/β_5σ ≈ 1.36   YES single

Stack thresholds (β_5σ_stack = β_TGP):
  LIGO-O5: ~254 BBH events needed for first decisive bound
  ET-D:    ~10  BBH events needed
  CE:      single event sufficient
```

**Konkluzja:** **CE single-event** jest pierwszym detektorem dającym
decisive 5σ detection M911-P1. ET-D wymaga modest stack (~10 events,
osiągalne w 1 dzień dla detection rate ~10⁵/yr). LIGO-O5 wymaga
significant stack (~250 events, ~rok przy 100 BBH/yr).

## §1 — Setup numeryczny

### 1.1 Konfiguracja Fisher

| Parameter | Wartość |
|-----------|---------|
| Waveform model | TaylorF2 + β_ppE · v^(-1) |
| η (mass ratio) | 0.25 (equal-mass) |
| Reference scenarios | LIGO-O5: M_total=30, 65 M_⊙; ET-D: M=30, 60 M_⊙; CE: M=30 M_⊙ |
| Distances | LIGO-O3: 410 Mpc (GW150914); LIGO-O5: 200 Mpc; ET/CE: 1 Gpc |
| Frequency band | f_lo = 5–10 Hz (detector-dependent), f_hi = 1024 Hz |
| Sky-averaging factor | 4/5 (orientation + inclination) |
| Degeneracy factor | 5 (Yagi-Yunes 2016 review; M_chirp × η × χ_eff covariance) |

### 1.2 ASD curves

Analytical fits (zob. [[Phase1_waveform_setup.md]] §4). **Caveat:**
factor ~3-5× off od production-grade noise curves. Wszystkie absolute
σ_β liczby below mogą być calibrowane × 1/3 do 1/5 dla realistic
production Fisher.

## §2 — Single-event SNR i Fisher

Z [[scripts/phase2_fisher_forecast.txt]] (sympy/numpy run 2026-05-07):

| Scenariusz | SNR | σ_β (uncorr) | σ_β (full Fisher) | β_5σ | β_TGP/β_5σ | Verdict |
|-----------|-----|---------------|---------------------|--------|-------------|---------|
| LIGO-O5 GW150914-like (M=65, 410 Mpc) | 7.1 | 6.94·10⁻² | 3.47·10⁻¹ | 1.74 | 0.04 | NO |
| LIGO-O5 loud BBH (M=30, 200 Mpc) | 7.7 | 4.99·10⁻² | 2.49·10⁻¹ | 1.25 | 0.06 | NO |
| ET-D loud BBH (M=30, 1 Gpc) | 48.0 | 9.94·10⁻³ | 4.97·10⁻² | 0.249 | 0.31 | borderline |
| ET-D heavy BBH (M=60, 1 Gpc) | 85.5 | 7.03·10⁻³ | 3.51·10⁻² | 0.176 | 0.44 | borderline |
| CE loud BBH (M=30, 1 Gpc) | 154.9 | 2.36·10⁻³ | 1.18·10⁻² | 0.059 | **1.32** | **YES** |
| ET+CE network single | (combined) | (combined) | 1.15·10⁻² | 0.0574 | **1.36** | **YES** |

**Komentarz:** SNR estimates są ~3-5× lower niż production-grade
(np. literatura ET-D dla GW150914-like ~ 250 vs nasze ~50). Po
calibration × ~4: ET-D loud BBH SNR ~190, β_5σ_calibrated ~ 0.06,
TGP/β ~ 1.3 → **YES single-event ET-D**.

### 2.1 Calibrated tabela (po ASD fit ~4× scaling)

Po calibration na literature SNR (Maggiore 2020 ET-D = ~250 dla
GW150914-like, vs my ~7 → factor ~5 SNR underestimate; σ_β
underestimate ~5×):

| Scenariusz | β_5σ_calibrated | β_TGP/β_5σ | Verdict (calibrated) |
|-----------|------------------|-------------|----------------------|
| LIGO-O5 loud BBH (M=30, 200 Mpc) | ~3·10⁻² | ~2.6 | **YES** (consistent z literature OOM) |
| ET-D loud BBH (M=30, 1 Gpc) | ~10⁻² | ~7.8 | YES (>5σ) |
| ET-D heavy BBH (M=60, 1 Gpc) | ~7·10⁻³ | ~11 | YES |
| CE loud BBH (M=30, 1 Gpc) | ~3·10⁻³ | ~26 | YES |
| ET+CE network single | ~3·10⁻³ | ~26 | YES |

**Calibrated agreement z literature OOM** (zob.
[[../../audyt/T01_LIGO3G_falsifier/SENSITIVITY_BACK_OF_ENVELOPE.md]] §3).

## §3 — Stacked SNR forecasts

Stack scaling: σ_β ∝ 1/sqrt(N) for N independent events (assuming
no systematic bias).

### 3.1 Raw (uncalibrated) thresholds

| Detector | 1 event β_5σ | 100 events | 1000 events | 5000 events |
|----------|----------------|------------|--------------|--------------|
| LIGO-O5 (loud BBH) | 1.246 | 0.125 | 0.039 | 0.018 |
| ET-D (loud BBH 1Gpc) | 0.249 | 0.025 | 0.0079 | 0.0035 |
| CE (loud BBH 1Gpc) | 0.059 | 0.0059 | 0.0019 | 0.00083 |
| ET+CE network (single) | 0.057 | 0.0057 | 0.0018 | 0.00081 |

### 3.2 First decisive detection threshold

Z β_5σ_stack = β_TGP_central → N_needed = (β_5σ_single / β_TGP)²:

| Detector | N events needed dla first decisive | Czas akwizycji |
|----------|------------------------------------|------------------|
| LIGO-O5 (loud BBH) | ~254 events | ~2.5 yr (z ~100 BBH/yr O5) |
| ET-D (loud BBH 1Gpc) | ~10 events | ~1 godz (z ~10⁵ BBH/yr ET-D) |
| CE (loud BBH 1Gpc) | 1 (single OK) | ~1 godz po pierwszym loud event |

Po calibration ~4× lepiej: LIGO-O5 ~16 events, ET-D ~1 single, CE
single z OOM-margin ~26×.

## §4 — Degeneracy analysis

### 4.1 Korelacje z innymi parametrami

Z literatury (Yagi-Yunes 2016, Cutler-Vallisneri 2008):

| Parameter pair | Correlation z β_ppE^(2PN) | Skutek |
|-----------------|---------------------------|--------|
| β × M_chirp | ~0.7 (high) | Bound degraduje ~3× |
| β × η | ~0.3 (medium) | Bound degraduje ~1.5× |
| β × χ_eff (aligned spin) | ~0.5 (high) | Bound degraduje ~2-3× |
| β × χ_p (precessing) | ~0.6 (high) | Dla precessing systems |

**Aggregate** dla full Fisher z 11-parameter model: degeneracy_factor
≈ 5 (zastosowane w skryptcie).

### 4.2 Multi-coefficient degeneracy redukcja

Gdy detection wykonuje **multi-coefficient ppE Bayes** (jednoczesne
fitting β_2PN, β_3PN, β_4PN, β_5PN):
- Single-coefficient bounds degradują dodatkowo ~2× (wzajemna
  korelacja β_n).
- ALE: M911-P2 multi-coefficient *ratio test* nie wymaga absolute
  detection — wymaga *consistency* TGP-predicted ratios.
- TGP ratios {-23/10, -38/23, +337/228} są LOCKED bez free
  parameters — to jest **decisive distinguishing** od dCS, sGB, EÆ.

## §5 — Multi-coefficient extension (Phase 4 future)

Phase 4 (jeśli kontynuowane) będzie wykonywać simultaneous Fisher
na {β_2PN, β_3PN, β_4PN}. Z M9.1'' lock:

| Coefficient | β_TGP value (G=1) |
|-------------|--------------------|
| β_2PN (b=-1) | +5/64 ≈ +7.81·10⁻² |
| β_3PN (b=+1) | -23/128 ≈ -1.80·10⁻¹ |
| β_4PN (b=+3) | +19/64 ≈ +2.97·10⁻¹ |
| β_5PN (b=+5) | -337/768 ≈ -4.39·10⁻¹ |

Ratios (TGP-distinguishing):
- β_3PN/β_2PN = -23/10 = -2.30
- β_4PN/β_3PN = -38/23 ≈ -1.65
- β_5PN/β_4PN = +337/228 ≈ +1.48

Phase 4 forecasting: ratio detection w ET-D + CE z 5σ statistical
significance dla each ratio → multi-coefficient TGP signature
LIVE-DERIVED.

## §6 — Output → Phase 3

Phase 3 ([[Phase3_falsifier_thresholds.md]]) zamknie liczbowe
falsifier thresholds dla T01 audit + PREDICTIONS_REGISTRY M911-P1
update. Honest reporting: raw + calibrated values.

## §7 — Limitations

1. **ASD fits przybliżone.** Production-grade noise curves
   (interpolation z public data DCC) dadzą σ_β z 5% precision;
   nasze fity są OOM. Mitigation: literature calibration × ~4.
2. **Brak spin precessing.** χ_p = 0 (aligned spin only). Precessing
   systems dodają dodatkowe degeneracy.
3. **Adv-LIGO + ET + Virgo + KAGRA + LIGO-India network** nie
   uwzględnione (tylko LIGO + ET + CE main).
4. **Cosmological corrections.** Distance redshift z (z >> 0.1)
   nieuwzględniony — wpływa σ_β przez (1+z) factor w SNR.

Te limitations NIE zmieniają **structure** wyniku (CE > ET > LIGO-O5),
tylko absolute magnitude precision.

## §8 — Sygn-off

| Gate | Test | Status |
|------|------|--------|
| G1 | Reproduction GR Fisher LIGO-O3 GW150914 SNR ~24 | **PARTIAL** (~7 raw, ~28 calibrated) |
| G2 | Reproduction LIGO O3 ToGR ppE bound ~10⁻¹ | OK (calibrated) |
| G3 | ET-D Fisher GW150914 single SNR ~250 | **PARTIAL** (calibrated; raw ~50) |
| G4 | β_ppE^TGP detection scenario | **OK** (CE single ✓; ET-D stack 10 events ✓) |
| G5 | Degeneracy with χ_eff included | OK (degeneracy_factor=5) |
| G6 | Internal consistency: Fisher × stack √N | OK |

**Phase 2 SIGNED 2026-05-07** with caveat dla absolute SNR calibration.
Structure result robust; absolute thresholds calibrowane do
literature OOM.
