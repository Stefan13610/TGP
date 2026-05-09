---
title: "Cycle kickoff brief — op-LIGO-3G-deviation (Path A; Fisher analysis ET-D + CE)"
date: 2026-05-07
parent: "[[README.md]]"
type: cycle-kickoff
tgp_owner: audyt/T01_LIGO3G_falsifier
target_cycle: research/op-LIGO-3G-deviation/
tags:
  - kickoff
  - cycle-brief
  - LIGO-3G
  - Einstein-Telescope
  - Cosmic-Explorer
  - Fisher
  - SNR
  - ppE
  - 2PN-phase
  - M911
  - T01
  - EXT-5
  - path-A
related:
  - "[[README.md]]"
  - "[[NEEDS.md]]"
  - "[[CONVENTION_DECISION.md]]"
  - "[[SENSITIVITY_BACK_OF_ENVELOPE.md]]"
  - "[[CYCLE_KICKOFF_op-ppE-mapping.md]]"
  - "[[FALSIFIER_STATEMENT_DRAFT.md]]"
---

# Cycle kickoff — `research/op-LIGO-3G-deviation/`

> **Cel pliku.** Dostarczyć **gotowy do skopiowania** brief
> uruchomieniowy dla cyklu `research/op-LIGO-3G-deviation/`, który
> ma zamknąć ścieżkę A z [[README.md]] (Fisher analysis ET-D + CE
> dla detekcji M911-P1). Brief jest *gotowy do startu* po zamknięciu
> [[CYCLE_KICKOFF_op-ppE-mapping.md]] Phase 1.4.
>
> **Polityka:** ten plik **przygotowuje** cykl, ale **nie tworzy go**.
> Również: ten cykl **DEPENDS ON** `op-ppE-mapping/` Phase 1.4
> (β_ppE^TGP^(b=−1) liczbowo zamknięty).

## §1 — Charakterystyka cyklu

| Pole | Wartość |
|------|---------|
| **Folder docelowy** | `research/op-LIGO-3G-deviation/` |
| **Klasa** | numerical-Fisher + statistical-falsifier |
| **Domyka** | T01 N1 (waveform mapping), N3 (Fisher), N4 (SNR thresholds), N5 (higher-PN), N7 (validity regime), N11 (basin spójność) |
| **Otwiera bridge** | nie (terminal Phase A; może się rozszerzyć w paper-track Path D) |
| **Estymata pracy** | ~3–6 sesji numerycznych (1 sesja = 4–6 godzin pracy + numerical) |
| **Output główny** | `Phase2_Fisher_results.md` — SNR thresholds + Fisher matrix + degeneracy analysis |
| **Output secondary** | `Phase3_falsifier_thresholds.md` — locked threshold dla [[FALSIFIER_STATEMENT_DRAFT.md]] |
| **Dependencies** | **HARD:** `op-ppE-mapping/` Phase 1.4 zamknięty (β_ppE^TGP). **SOFT:** S07 closure (awans z "test ansatzu" do "test derywowanej predykcji"). |

## §2 — Cele i kryteria sukcesu

### 2.1 Phase 1 (waveform model setup)

**Cel:** Wybór i implementacja waveform modelu z β_ppE^TGP^(b=−1)
deviation (rozstrzyga T01 problem #5 "dualizm modeli waveform").

**Kryteria PASS:**
- [P1.1] Wybór: **GR-base TaylorF2 + ppE single-coefficient deviation**
  (perturbative inspiral, no NR-TGP). Argumenty: M9.1'' two-body
  Lagrangian (z `op-ppE-mapping/` Phase 1.2) jest perturbative w
  (5/6) U³; full NR-TGP nie istnieje i jego budowa byłaby ~2-letnim
  cyklem; dla inspiral signal (f < f_ISCO) perturbative jest wystarczające.
- [P1.2] Implementacja w **bilby** lub **pyCBC** ppE plugin —
  reuse standardowych Bayes inference frameworks dla LIGO data.
- [P1.3] Inspiral cutoff: f_max = f_ISCO = c³/(6√6 π G M_total) lub
  niżej; wykluczać merger/ringdown (gdzie U → 1/2 i perturbative
  nie trzyma).

**Kryterium FAIL** (wymaga rewizji):
- TaylorF2 z ppE NIE reprodukuje M9.1'' two-body Lagrangian dynamics
  → wybór waveform modelu wymaga zaawansowanego (Yunes-style EOB-ppE
  hybrid).

### 2.2 Phase 2 (Fisher matrix analysis)

**Cel:** Single-event SNR threshold dla 5σ detekcji β_ppE^TGP^(b=−1)
w ET-D i CE.

**Kryteria PASS:**
- [P2.1] Fisher matrix dla 11-parameter model (M_chirp, η, χ_eff, χ_p,
  d_L, ι, ψ, t_c, φ_c, RA, Dec) + β_ppE^TGP^(b=−1).
- [P2.2] Sky-position averaging i inclination averaging (standardowo
  dla forecast).
- [P2.3] SNR threshold dla 5σ detekcji β_ppE^TGP^(b=−1) (z liczbowym
  β_ppE^TGP z `op-ppE-mapping/` Phase 1.4) — output *liczbowy* dla:
  - LIGO-O5 (A+ design): SNR_threshold_O5
  - ET-D single (Maggiore 2020): SNR_threshold_ET
  - CE single (Reitze 2019): SNR_threshold_CE
  - ET+CE network: SNR_threshold_ETCE
- [P2.4] Population study: dla każdej ET/CE detection rate (Maggiore
  2020 §2; ~10⁵ BBH/yr ET-D), policzyć N_events potrzebnych dla
  cumulative 5σ stack detection.
- [P2.5] Degeneracy analysis: **β_ppE × M_chirp**, **β_ppE × χ_eff**
  korelacje. Standard w literaturze: 30–80% bound degradation z spin.

**Kryterium FAIL:**
- TGP β_ppE^TGP_locked < ET-D single-event Fisher bound (5σ): sygnał,
  że trzeba większego stack lub wyższych masa (chirp); honest report.
- Degeneracy z spin > 90%: coefficient nie identyfikowalny single-event;
  wymaga multi-event stacking.

### 2.3 Phase 3 (falsifier threshold lock)

**Cel:** Zamknąć liczbowo *threshold* dla wpisu falsifier w
`PREDICTIONS_REGISTRY.md`.

**Output (numerical):**
- `β_ppE^TGP_threshold_5σ_ET-D = X.X · 10^(−Y)` (single event GW150914-like)
- `β_ppE^TGP_threshold_5σ_ET-D_stack_100 = X.X · 10^(−Y)`
- `β_ppE^TGP_threshold_5σ_ET+CE_stack_5000 = X.X · 10^(−Y)`
- Kompletna tabela: M_chirp ∈ {10, 20, 30, 50} M_⊙ × {single, stack 30, stack 100, stack 1000}.

### 2.4 Phase 4 (multi-coefficient extension, opcjonalne)

**Cel:** Rozszerzyć detekcję na M911-P2 (β_ppE_3PN_phase = 23/12 → liczbowo)
i M911-P3 (β_ppE_4PN_phase). Daje multi-coefficient TGP signature
(zob. [[PPN_TO_PPE_MAPPING.md]] §3 last paragraph).

**Wymaga:** `op-ppE-mapping/` Phase 1.4 mieć policzone *wszystkie*
β_(N-PN_phase) dla N = 2, 3, 4 (nie tylko 2PN).

## §3 — Setup numeryczny (gotowy do skopiowania jako Phase 0)

### 3.1 Tools i environment

```bash
# Python 3.11+, dependencies:
pip install bilby pycbc lalsimulation numpy scipy emcee dynesty
# ppE plugin (modyfikowany TaylorF2):
pip install ppEmcee  # (or local fork; ppE often custom-implemented)
```

### 3.2 Reference noise curves (download)

| Detector | Source | URL |
|----------|--------|-----|
| LIGO-O5 (A+) | LIGO Scientific Collaboration 2020 | DCC: T1800042 / aLIGOAPlusDesignSensitivityT1800042 |
| ET-D | ET design study (Hild 2010) | gwic.ligo.org / et-design |
| Cosmic Explorer | Reitze 2019 | cosmicexplorer.org/data |

Cykl powinien zacommitować copy noise curves do
`research/op-LIGO-3G-deviation/data/` z citation.

### 3.3 Sample populations (BBH/BNS)

Z LIGO O3 catalog (GWTC-3) i synthetic populations (Madau-style):
- **GW150914-clone:** M_total = 65 M_⊙, η = 0.247, χ_eff = -0.02, d_L = 410 Mpc
- **GW170817-clone:** M_total = 2.74 M_⊙, η = 0.249, χ_eff = 0.005, d_L = 40 Mpc
- **Loud-BBH ET archetype:** M_total = 30 M_⊙, η = 0.25, χ_eff = 0,
  d_L variable (200 Mpc do 4 Gpc)
- **EMRI archetype** (LISA cross-check): M_1 = 10⁶ M_⊙, M_2 = 10 M_⊙

### 3.4 Skrypty struktura

```
research/op-LIGO-3G-deviation/
├── README.md
├── Phase0_balance.md
├── Phase1_waveform_setup.md           — ppE TaylorF2 implementation
├── Phase2_fisher_setup.py              — Fisher matrix script
├── Phase2_fisher_runs/                 — output dla każdego scenariusza
│   ├── ligo-o5_loudBBH_30Msun_200Mpc/
│   ├── et-d_loudBBH_30Msun_1Gpc/
│   ├── ce_loudBBH_30Msun_1Gpc/
│   └── etce-network_stack/
├── Phase2_results.md                   — tabela SNR thresholds
├── Phase3_falsifier_thresholds.md      — locked numbers dla rejestru
├── Phase4_multicoef.md                 — opcjonalne (M911-P2, P3)
├── data/                                — noise curves
├── scripts/                             — auxiliary
└── NEEDS.md                             — open follow-ups
```

## §4 — Wymagania techniczne

### 4.1 Skill set (subagent / autor)

- LIGO data analysis (LSC pipeline familiarity).
- Bayesian inference (bilby / pycbc).
- ppE framework (Yunes–Pretorius operational).
- Fisher matrix forecasting (Cutler-Flanagan + Vallisneri 2008
  on Fisher pitfalls).
- Numerical Python proficiency.

### 4.2 Estymowane sesje

| Phase | Praca | Szacunkowo godzin |
|-------|-------|-------------------|
| Phase 1 (waveform setup) | TaylorF2 + ppE plugin | 4–6 h |
| Phase 2 (Fisher single events) | run + analysis ×4 detectors | 6–10 h |
| Phase 2 (degeneracy & spin) | extended Fisher | 3–4 h |
| Phase 2 (population study) | synthetic samples | 3–5 h |
| Phase 3 (threshold lock + reporting) | tabularization | 2–3 h |
| Phase 4 (multi-coefficient, optional) | extended Fisher | 4–6 h |
| **Razem** | | **~22–34 h, ~4–6 sesji** |

## §5 — Walidacja i gates

### 5.1 Pre-commit gate (Phase 0 balance check)

Zgodnie z **CALIBRATION_PROTOCOL § ABSOLUTE BINDING gate**:

`Phase0_balance.md` zawierający:
- inputs: noise curves, ppE plugin code, β_ppE^TGP value (locked
  z `op-ppE-mapping/` Phase 1.4), M9_1_pp_P1 baseline.
- outputs: SNR thresholds, falsifier thresholds.
- predyktywność ratio: N_outputs (Phase 2 + Phase 3 + Phase 4)
  / N_locked_inputs.
- epistemic class: NUMERICAL (Fisher matrix → forecasting; nie
  DERIVED w sensie analitycznym).

### 5.2 Output validation gates

| Gate | Test | Pass criterion |
|------|------|----------------|
| G1 | Reproduction GR Fisher LIGO-O3 GW150914 | within 5% literature |
| G2 | Reproduction LIGO O3 ToGR ppE bound | within 30% (Abbott et al. 2021) |
| G3 | ET-D Fisher GW150914 single | within 5% Maggiore 2020 |
| G4 | β_ppE^TGP detection scenario | clean 5σ for at least one configuration |
| G5 | Degeneracy with χ_eff | bound degradation < 80% (limit pełnego mascowania) |
| G6 | Internal consistency: Fisher × stack √N | reproduces from single events |

## §6 — Output → T01 closure path

Po zamknięciu `op-LIGO-3G-deviation/` Phase 3:

1. **Update** [[FALSIFIER_STATEMENT_DRAFT.md]] §1 tabela detekcji:
   zaktualizować "OOM" do liczb z Phase 3.
2. **Update** [[FALSIFIER_STATEMENT_DRAFT.md]] §1 falsifier clause:
   zaktualizować threshold "[β_th]" liczbą + scenariusze
   M_chirp / SNR / N_events.
3. **Update** [[SENSITIVITY_BACK_OF_ENVELOPE.md]]: oznaczyć jako
   superseded (lub zachować jako preliminary OOM, dodać
   "actually-locked" tabelę).
4. **Update** [[NEEDS.md]] N1, N3, N4, N5, N7, N11: zamknąć status.
5. **Update** [[README.md]] sekcja "Postęp domknięcia" tabela:
   "Path A SENSITIVITY preview" → "Path A EXECUTED, locked".
6. **(Autor)** Publikować draft paper "M9.1'' 2PN-phase deviation:
   ET/CE forecast" (Path D).

T01 status awansuje wtedy z **CLOSED-PARTIAL** do **CLOSED-EXECUTED
(wszystkie 4 ścieżki A/B/C/D wykonane lub w trakcie)**.

## §7 — Ryzyka i fall-back

| Ryzyko | Prawdopodobieństwo | Impact | Mitigation |
|--------|--------------------|--------| -----------|
| ppE plugin nie obsługuje custom b_ppE = -1 (bilby/pycbc default może być limited) | medium (40%) | medium (custom code 1–2h dodatkowe) | implement local fork (~50 lines Python) |
| TGP β_ppE^TGP poniżej ET single-event bound | medium (30%, OOM ~10⁻¹ vs bound ~10⁻³ → margin 100×, ale Fisher precise może zmienić) | low (stack rozwiązuje) | move to stack analysis Phase 2.4 |
| Degeneracy z χ_eff całkowita | low (10%) | high (TGP nie identyfikowalne single-event) | use BNS (low spin) or precise event selection |
| Phase 1.4 dependency NIE zamknięta | medium (zależy od op-ppE-mapping timeline) | very high (cykl czeka) | sekwencyjne uruchomienie; nie startować Path A przed B |
| LIGO-O3 reanalysis (Q3 z [[NEEDS.md]]) sfalsyfikuje TGP już teraz | very low (OOM safe margin) | terminal (TGP M9.1'' upada) | honest report; możliwy pivot do S07 ścieżki C (przyznanie statusu (P)) |

## §8 — Powiązanie z innymi cyklami

- **`op-ppE-mapping/`** (Path B) — **HARD dependency** Phase 1.4.
  Bez β_ppE^TGP locked, Path A nie ma co testować.
- **S07 (M9.1'' derivation)** — **SOFT dependency** Phase 4 (paper).
  Po S07 closure, paper Path D awansuje z "test ansatzu" do "test
  derywowanej teorii".
- **GW1, GW2, GW6, BH5** (PREDICTIONS_REGISTRY) — **orthogonal
  channels**; TGP-pattern multi-coefficient + BH5 QNM ringdown
  daje multi-prong falsifier-set.
- **LISA EMRI** — przyszły rozszerzający bridge (M911-P_LISA).
  EMRI U → 0.3 — strong field, multi-coefficient pattern detektywny
  na pojedynczym zdarzeniu.

## Cross-references

- [[README.md]] — diagnoza T01 problemu #1 (sensitivity), #5 (waveform model)
- [[NEEDS.md]] — N1, N3, N4, N5, N7, N11 (closed by this cycle)
- [[CONVENTION_DECISION.md]] — adopted phase-PN convention
- [[SENSITIVITY_BACK_OF_ENVELOPE.md]] — preview tego cyklu
- [[CYCLE_KICKOFF_op-ppE-mapping.md]] — predecessor (HARD dependency)
- [[FALSIFIER_STATEMENT_DRAFT.md]] — target update po Phase 3
- [[../S07_M911_derivation/README.md]] — meta-bloker (test ansatzu vs test fundamentu)
- [[../../meta/CALIBRATION_PROTOCOL.md]] — Phase 0 gate enforcement
- [[../../PREDICTIONS_REGISTRY.md]] — target rejestru BH5, ε.1, GW1-6 cross-check

## Bibliografia operacyjna

- LIGO Scientific Collaboration, Phys. Rev. D **103**, 122002 (2021), arXiv:2010.14529 — GWTC-2 ppE constraints baseline
- LIGO Scientific Collaboration, *GWTC-3 ppE updates*, arXiv:2112.06861 — current bounds
- Maggiore et al., JCAP **03**, 050 (2020), arXiv:1912.02622 — ET science case
- Reitze et al., Bull. Am. Astron. Soc. **51**, 035 (2019) — CE
- Chamberlain & Yunes, Phys. Rev. D **96**, 084039 (2017), arXiv:1704.08268 — 3G implications baseline
- Vallisneri, Phys. Rev. D **77**, 042001 (2008) — Fisher matrix pitfalls (REQUIRED reading)
- Hild et al., Class. Quantum Grav. **27**, 015003 (2010) — ET-D design
- Ashton et al., Astrophys. J. Suppl. **241**, 27 (2019) — bilby framework
- Biwer et al., PASP **131**, 024503 (2019) — pycbc framework
