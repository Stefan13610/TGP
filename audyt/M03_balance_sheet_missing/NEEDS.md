---
title: "NEEDS — M03 brak balance sheet pre-74394a8 cykli"
date: 2026-05-04
parent: "[[README.md]]"
type: needs
tgp_owner: audyt/M03_balance_sheet_missing
tags:
  - needs
  - balance-sheet
  - retrofit
  - methodology
---

# NEEDS — M03 brak balance sheet dla 27+ pre-74394a8 cykli

## Otwarte luki

| ID | Luka | Typ | Source | Kandydat dostawcy |
|----|------|-----|--------|-------------------|
| N1 | Balance sheet retrofit dla high-risk mini-cycles (ε.1, ζ.1, θ.1, η.1, η.2, α.1, κ.1, ι.1, μ.1) | audit-method | brak | nowy cykl `op-balance-sheet-retrofit-all` Phase 1+2 |
| N2 | Balance sheet retrofit dla medium-risk cykli (ψ.1, σ.1, τ.1, τ.2, ω.1, BH.1, SC.1, XS.1) | audit-method | brak | cykl Phase 2 |
| N3 | Balance sheet retrofit dla low-risk cykli (UV.1, υ.1, ξ.2, ο.1, π.1, ρ.1, φ.1) | audit-method | brak | cykl Phase 2 |
| N4 | PREDICTIONS_REGISTRY refactor z post-retrofit klasyfikacjami | editorial | rekomendacja Phase 3 | cykl Phase 3 |
| N5 | Counter master 856 → split na DERIVED FULL / CONDITIONAL / NUMEROLOGICAL / STRUCTURAL | editorial | brak | cykl Phase 3 |
| N6 | Re-derive predictivity ratio z post-retrofit klasyfikacjami | analytical | LS-8 audit | cykl Phase 3 |
| N7 | CALIBRATION_PROTOCOL lock jako absolute binding gate dla future cycles | policy | [[../../meta/CALIBRATION_PROTOCOL.md]] | cykl Phase 4 (policy decision) |

## Pytania otwarte

- **Q1**: Czy mini-cycles (ε.1 ... μ.1) używają wspólnej struktury
  „mixing-operator framework" z cascade locks? Jeśli tak, są podejrzane
  o ten sam pattern co χ.1/UV.2.
  - Source: PREDICTIONS_REGISTRY head 80 (status cascade activated)
  - Wpływa na: priority sorting w Phase 1

- **Q2**: Ile z 784 effective uncontested entries faktycznie przejdzie
  tautology + falsifiability test?
  - Source: brak (otwarte pytanie)
  - Wpływa na: faktyczna predyktywność TGP

- **Q3**: Czy retrospect rollback (Option A) czy forward-patch markery
  (Option B) jest lepsze dla pre-74394a8 cykli?
  - Source: powiązane z M02 decyzją
  - Wpływa na: integralność rejestru

- **Q4**: Czy Calibration Protocol jest *naprawdę* binding od 2026-05-04+,
  czy tylko aspirational? Jakie są konsekwencje gdy subagent commitnie
  cykl bez balance sheet?
  - Source: CALIBRATION_PROTOCOL policy
  - Wpływa na: prevention futurystycznego status creep

## Blokery

| ID | Bloker | Czeka na | Source |
|----|--------|----------|--------|
| B1 | Retrofit 27+ cykli to ~60 dni pracy — wymaga dedicated theorist time | autor + dedicated audytor | rekomendacja M03 |
| B2 | M02 decision (rollback vs forward-patch) musi być *przed* M03 retrofit żeby uniknąć duplikacji | [[../M02_ledger_pollution]] | rekomendacja M03 |
| B3 | Phase 4 policy lock wymaga decyzji: czy CALIBRATION_PROTOCOL ma „veto power" nad subagentem? | autor decision | [[../../meta/CALIBRATION_PROTOCOL.md]] |

## Closed needs

| ID | Co | Domknięte przez | Data |
|----|----|-----------------|------|
| (puste) | | | |
