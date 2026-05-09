---
title: "op-GWTC3-reanalysis — TGP M911-P1 Bayes factor analysis na public LIGO/Virgo GWTC-3 data"
date: 2026-05-07
parent: "[[../README.md]]"
type: research-cycle
tgp_owner: research/op-GWTC3-reanalysis
status: active
tags:
  - research
  - cycle
  - GWTC-3
  - Bayes-factor
  - TGP
  - M911
  - 2PN-phase
  - ppE
  - reanalysis
  - tier-5
  - T01
related:
  - "[[../op-ppE-mapping/Phase1_results.md]]"
  - "[[../op-LIGO-3G-deviation/Phase3_falsifier_thresholds.md]]"
  - "[[../../audyt/T01_LIGO3G_falsifier/FALSIFIER_STATEMENT_DRAFT.md]]"
  - "[[../../audyt/T01_LIGO3G_falsifier/FINDINGS.md]]"
  - "[[../../PREDICTIONS_REGISTRY.md]]"
tgp_status:
  folder_status: paused
  level: research
  kind: research-cycle
  core_compatibility: side-channel
  may_edit_core: false
  exports_findings: true
  has_needs_file: false
  has_findings_file: false
  open_bridges: []
  depends_on: ["op-ppE-mapping/Phase1_results.md (β_ppE^TGP locked)", "op-LIGO-3G-deviation/Phase2_results.md (Fisher framework)"]
  impacts: ["audyt/T01_LIGO3G_falsifier (FINDINGS Tier 5 closure)", "PREDICTIONS_REGISTRY (M911-P1 status pre-confirmation possible)"]
  source_of_status:
    - "[[../../audyt/T01_LIGO3G_falsifier/FINDINGS.md]] §2.5 Tier 5"
    - "[[../../audyt/T01_LIGO3G_falsifier/SESSION_REPORT_2026-05-07_C-B-A-D.md]] §10"
  promoted_to_core: null
  polluted_74394a8: false
  pre_existing_findings: false
  pre_existing_needs: false
  last_yaml_update: 2026-05-07
---

# op-GWTC3-reanalysis — TGP M911-P1 reanalysis na GWTC-3

## Cel cyklu

**Side-channel pre-publication confirmation** dla T01 M911-P1
predykcji `β_ppE^TGP^(b=-1) = -5/64 ≈ -7.81·10⁻²` (LOCKED z
[[../op-ppE-mapping/Phase1_results.md]]).

Aktualne LIGO/Virgo/KAGRA Tests of GR (ToGR) papers (Abbott et al.
2021 PRD 103:122002 GWTC-2; Abbott et al. 2023 PRX 13:041039 GWTC-3)
podają bounds dla **fractional 2PN-phase deviation** δφ_4 (lub
δφ̂_4 w niektórych notacjach). Z calibrated Fisher analysis w Path A
([[../op-LIGO-3G-deviation/Phase3_falsifier_thresholds.md]] §2):
LIGO-O3 stack 100 BBH ma β_5σ ~10⁻², a TGP β_TGP ~7.8·10⁻² → stosunek
~7.8 → **decisive detection możliwa już z aktualnych public data**.

**Cykl pyta:** Czy GWTC-3 ToGR public bounds wykluczają lub
*sygnalizują* M9.1'' deviation?

## Approach

**Honest scope:** Pełna reanalysis z LIGO H5 strain data + bilby/
pycbc Bayesian inference jest poza zakresem 1-sesyjnego cyklu
(wymaga ~50+ godz pracy + access do danych). Ten cykl wykonuje
**literature-grade reanalysis**:

1. Zebrać public **published 2PN-phase bounds** z Abbott et al. 2023
   PRX 13:041039 (Tables III, IV) + Abbott et al. 2021 PRD 103:122002.
2. Konwersja jednostek: LIGO `δφ̂_n` (fractional) ↔ ppE `β_ppE^(b=-1)`
   (absolute) via SPA chain prefactor.
3. Compute residual TGP-prediction vs published combined GWTC-3
   posterior.
4. Estimate Bayes factor TGP vs GR using Laplace approximation
   on published 2PN-phase posteriors.
5. Report verdict: TGP (a) consistent + room for detection, (b)
   tentative signal at low significance, (c) rejected/excluded,
   or (d) within 1σ of best-fit.

## Status

**Aktywny od 2026-05-07.** Phase 0 + Phase 1 + Phase 2 + Phase 3
EXECUTED w sesji uruchomienia.

## Phases

| Phase | Plik | Status |
|-------|------|--------|
| 0 | [[Phase0_balance.md]] | EXECUTED 2026-05-07 |
| 1 | [[Phase1_GWTC3_bounds.md]] | EXECUTED 2026-05-07 |
| 2 | [[Phase2_Bayes_factor.md]] | EXECUTED 2026-05-07 |
| 3 | [[Phase3_verdict.md]] | EXECUTED 2026-05-07 |

## Cross-references

- [[../../audyt/T01_LIGO3G_falsifier/FINDINGS.md]] §2.5 Tier 5 — recommendation
- [[../../audyt/T01_LIGO3G_falsifier/SESSION_REPORT_2026-05-07_C-B-A-D.md]] §10 — pierwszy z trzech recommended
- [[../op-ppE-mapping/Phase1_results.md]] — β_ppE^TGP input
- [[../op-LIGO-3G-deviation/Phase3_falsifier_thresholds.md]] — Fisher framework
- [[../../PREDICTIONS_REGISTRY.md]] M911-P1 — target rejestru update

## Bibliografia kluczowa

- Abbott et al. (LIGO/Virgo) 2021, *Phys. Rev. D* **103**, 122002
  (GWTC-2 ToGR; arXiv:2010.14529) — primary bounds source
- Abbott et al. (LIGO/Virgo/KAGRA) 2023, *Phys. Rev. X* **13**,
  041039 (GWTC-3 ToGR; arXiv:2112.06861) — primary bounds source
- Yunes & Pretorius 2009, *Phys. Rev. D* **80**, 122003 — ppE
  framework
- Cutler & Flanagan 1994, *Phys. Rev. D* **49**, 2658 — SPA
- Mishra, Iyer, Sundararajan 2016, *Phys. Rev. D* **93**, 084054 —
  3PN waveform reference
