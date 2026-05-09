---
title: "op-ppE-mapping — Mapowanie M9.1'' (5/6) U³ deviation na ppE phase parameter"
date: 2026-05-07
parent: "[[../README.md]]"
type: research-cycle
tgp_owner: research/op-ppE-mapping
status: active
tags:
  - research
  - cycle
  - ppE
  - 2PN-phase
  - M911
  - SPA
  - inspiral
  - T01
  - EXT-5
  - path-B
related:
  - "[[../../audyt/T01_LIGO3G_falsifier/CYCLE_KICKOFF_op-ppE-mapping.md]]"
  - "[[../../audyt/T01_LIGO3G_falsifier/PPN_TO_PPE_MAPPING.md]]"
  - "[[../../audyt/T01_LIGO3G_falsifier/CONVENTION_DECISION.md]]"
  - "[[../../audyt/T01_LIGO3G_falsifier/FALSIFIER_STATEMENT_DRAFT.md]]"
  - "[[../op-newton-momentum/M9_1_pp_P1_results.md]]"
  - "[[../op-newton-momentum/M9_1_pp_setup.md]]"
  - "[[../../PREDICTIONS_REGISTRY.md]]"
tgp_status:
  folder_status: paused
  level: research
  kind: research-cycle
  core_compatibility: extension
  may_edit_core: false
  exports_findings: true
  has_needs_file: true
  has_findings_file: false
  open_bridges: []
  depends_on: ["op-newton-momentum (M9.1'' P1 sympy LOCK 5/5)", "T01 audit kickoff"]
  impacts: ["audyt/T01_LIGO3G_falsifier (FALSIFIER_STATEMENT_DRAFT [β_th] lock)", "PREDICTIONS_REGISTRY (M911-P1, M911-P2 entry update)", "research/op-LIGO-3G-deviation (Phase 1.4 input)"]
  source_of_status:
    - "[[../../audyt/T01_LIGO3G_falsifier/CYCLE_KICKOFF_op-ppE-mapping.md]]"
  promoted_to_core: null
  polluted_74394a8: false
  pre_existing_findings: false
  pre_existing_needs: false
  last_yaml_update: 2026-05-07
---

# op-ppE-mapping — Mapowanie M9.1'' (5/6) U³ deviation na ppE phase parameter

## Cel cyklu

Wyprowadzić **liczbowo** współczynniki ppE Yunes–Pretorius dla
fazy waveformu inspiralu BBH w M9.1'' kanonicznej metryce
hiperbolicznej `g_tt = -c²·(4-3ψ)/ψ`. Konkretnie:

- β_ppE^TGP^(b=−1) (2PN-phase coefficient) liczbowo dla equal-mass
  BBH (η = 1/4),
- β_ppE^TGP^(b=+1) (3PN-phase coefficient),
- β_ppE^TGP^(b=+3) (4PN-phase coefficient),
- multi-coefficient TGP-distinguishing pattern {β_2PN, β_3PN, β_4PN}.

**Wejście:** M9.1'' analytic c_n=2..7 z [[../op-newton-momentum/M9_1_pp_P1_results.md]] §3.2 (sympy LOCK 5/5).

**Wyjście:** β_ppE^TGP locked → update [[../../audyt/T01_LIGO3G_falsifier/FALSIFIER_STATEMENT_DRAFT.md]] §1 + [[../../PREDICTIONS_REGISTRY.md]] M911-P1 entry.

## Status

**Aktywny od 2026-05-07.** Phase 0 + Phase 1 PASS w sesji uruchomienia
(2026-05-07). Pełna orchestracja patrz `Phase0_balance.md` →
`Phase1_results.md`.

## Phases

| Phase | Plik | Status |
|-------|------|--------|
| 0 | [[Phase0_balance.md]] | EXECUTED 2026-05-07 |
| 1 | [[Phase1_results.md]] | EXECUTED 2026-05-07 |
| 2 | [[Phase2_literature_crosscheck.md]] | EXECUTED 2026-05-07 |
| 3 | [[Phase3_paper_ready.md]] | EXECUTED 2026-05-07 |

## Cross-references

- [[../../audyt/T01_LIGO3G_falsifier/CYCLE_KICKOFF_op-ppE-mapping.md]] — szczegółowy kickoff brief
- [[../../audyt/T01_LIGO3G_falsifier/PPN_TO_PPE_MAPPING.md]] — analytical preview audytu
- [[../../audyt/T01_LIGO3G_falsifier/CONVENTION_DECISION.md]] — PHASE convention
- [[../op-newton-momentum/M9_1_pp_P1_results.md]] — wejście (c_n=2..7)
- [[../op-newton-momentum/M9_1_pp_setup.md]] — M9.1'' setup
