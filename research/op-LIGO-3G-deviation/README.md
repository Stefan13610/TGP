---
title: "op-LIGO-3G-deviation — Fisher matrix forecasting M911-P1 dla ET-D + CE"
date: 2026-05-07
parent: "[[../README.md]]"
type: research-cycle
tgp_owner: research/op-LIGO-3G-deviation
status: active
tags:
  - research
  - cycle
  - LIGO-3G
  - Einstein-Telescope
  - Cosmic-Explorer
  - Fisher
  - ppE
  - 2PN-phase
  - M911
  - T01
  - EXT-5
  - path-A
related:
  - "[[../../audyt/T01_LIGO3G_falsifier/CYCLE_KICKOFF_op-LIGO-3G-deviation.md]]"
  - "[[../../audyt/T01_LIGO3G_falsifier/SENSITIVITY_BACK_OF_ENVELOPE.md]]"
  - "[[../op-ppE-mapping/Phase1_results.md]]"
  - "[[../../audyt/T01_LIGO3G_falsifier/FALSIFIER_STATEMENT_DRAFT.md]]"
  - "[[../../PREDICTIONS_REGISTRY.md]]"
tgp_status:
  folder_status: paused
  level: research
  kind: research-cycle
  core_compatibility: extension
  may_edit_core: false
  exports_findings: true
  has_needs_file: false
  has_findings_file: false
  open_bridges: []
  depends_on: ["op-ppE-mapping/Phase1_results.md (β_ppE^TGP locked)"]
  impacts: ["audyt/T01_LIGO3G_falsifier (FALSIFIER thresholds locked)", "PREDICTIONS_REGISTRY (M911-P1 LIVE thresholds)"]
  source_of_status:
    - "[[../../audyt/T01_LIGO3G_falsifier/CYCLE_KICKOFF_op-LIGO-3G-deviation.md]]"
  promoted_to_core: null
  polluted_74394a8: false
  pre_existing_findings: false
  pre_existing_needs: false
  last_yaml_update: 2026-05-07
---

# op-LIGO-3G-deviation — Fisher matrix forecasting M911-P1

## Cel cyklu

Wyznaczyć **Fisher matrix forecasts** dla detekcji M911-P1
(β_ppE^TGP^(b=-1) = -5/64 z [[../op-ppE-mapping/Phase1_results.md]])
w 4 detektorach 3G era:

- LIGO-O5 (A+ design, ~2027)
- Einstein Telescope (ET-D, ~2035)
- Cosmic Explorer (CE, ~2035+)
- ET+CE network

Konkretnie:
- SNR thresholds dla 5σ detekcji β_ppE^(b=-1) single-event,
- Stacked SNR dla N events (BBH/yr rate w każdym detector),
- Degeneracy z M_chirp, η, χ_eff (parameter covariance),
- Falsifier liczbowy thresholds dla
  [[../../audyt/T01_LIGO3G_falsifier/FALSIFIER_STATEMENT_DRAFT.md]]
  i [[../../PREDICTIONS_REGISTRY.md]] M911-P1.

**Dependency:** β_ppE^TGP value z `op-ppE-mapping` Phase 1.4 (LOCKED).

## Status

**Aktywny od 2026-05-07.** Phase 0 + Phase 1 + Phase 2 EXECUTED w
sesji uruchomienia. Pełna orchestracja patrz `Phase0_balance.md` →
`Phase2_results.md`.

## Phases

| Phase | Plik | Status |
|-------|------|--------|
| 0 | [[Phase0_balance.md]] | EXECUTED 2026-05-07 |
| 1 | [[Phase1_waveform_setup.md]] | EXECUTED 2026-05-07 |
| 2 | [[Phase2_results.md]] | EXECUTED 2026-05-07 |
| 3 | [[Phase3_falsifier_thresholds.md]] | EXECUTED 2026-05-07 |

## Cross-references

- [[../../audyt/T01_LIGO3G_falsifier/CYCLE_KICKOFF_op-LIGO-3G-deviation.md]] — szczegółowy kickoff brief
- [[../op-ppE-mapping/Phase1_results.md]] — β_ppE^TGP input source
- [[../../audyt/T01_LIGO3G_falsifier/SENSITIVITY_BACK_OF_ENVELOPE.md]] — preview audytu
