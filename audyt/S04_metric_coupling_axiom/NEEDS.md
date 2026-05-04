---
title: "NEEDS — S04 ax:metric-coupling vs L_mat"
date: 2026-05-04
parent: "[[README.md]]"
type: needs
tgp_owner: audyt/S04_metric_coupling_axiom
tags:
  - needs
  - axiom
  - L_mat
  - critical
---

# NEEDS — S04 ax:metric-coupling vs L_mat

## Otwarte luki

| ID | Luka | Typ | Source | Kandydat dostawcy |
|----|------|-----|--------|-------------------|
| N1 | Formal kowariantna derywacja `L_mat = -(q/Φ_0)·φ·ρ` z czysto metrycznego sprzęgania `L[ψ_m, g_eff]` | derivation | `sek08a` lin. 74–83 — brak formalnej derywacji | nowy cykl `op-A4-formal-derivation` Phase 1 |
| N2 | Test piątej siły TGP: gradient dφ/dr · ρ wokół solitona vs MICROSCOPE 10⁻¹⁵ | numerical | brak (audit § A.4) | cykl A4 Phase 2 |
| N3 | Mapowanie TGP `q/Φ_0` na Brans–Dicke ω_BD; sprawdzenie ω_BD > 40 000 (Cassini) | numerical | brak | cykl A4 Phase 3 |
| N4 | Test kompozycji M9.2 WEP: dwa źródła o różnym `q` lub niejednorodnym ρ | numerical | audit § B.9 explicit pending | cykl A4 Phase 2 (lub B6) |
| N5 | Formal sprawdzenie spójności Option-2 dla materii nierelatywistycznej (T^μ_μ ≠ ρ_rest) | derivation | audit § N.2 | cykl A4 Phase 1 |
| N6 | Cassini Shapiro re-test z `dφ/dr·ρ` korekta dyspersji sygnału | numerical | brak | cykl A4 Phase 3 |

## Pytania otwarte

- **Q1**: Czy `ax:metric-coupling` jest aksjomatem niezależnym od
  `ax:zrodlo`, czy konsekwencją? (Option-2 zakłada, że są kompatybilne;
  formalnie nie wykazane.)
  - Source: audit § N.1
  - Wpływa na: aksjomatyka §1 FOUNDATIONS

- **Q2**: Czy istnieje formal `L_mat[ψ_m, g_eff]` z którego po wariacji
  i całkowaniu (dla danego pola materii ψ_m, np. Diraca) wychodzi
  *dokładnie* `-(q/Φ_0)·φ·ρ` w sek08a, czy jest to interpretacyjny
  „handwave"?
  - Source: audit § N.2
  - Wpływa na: legitimacy całej akcji zunifikowanej

- **Q3**: Co właściwie znaczy `ρ` w L_mat dla pola Diraca (T^μ_μ = m·ψ̄ψ)
  vs pola elektromagnetycznego (T^μ_μ = 0)? Czy fotony nie sprzęgają
  z TGP?
  - Source: powiązane z [[../L01_rho_operational]]
  - Wpływa na: predykcje cosmologiczne (ψ.1, σ.1 — varying-c)

## Blokery

| ID | Bloker | Czeka na | Source |
|----|--------|----------|--------|
| B1 | A4 Phase 1 wymaga decyzji o `ρ` definicji (Q3) | [[../L01_rho_operational]] | audit § C.3 |
| B2 | Phase 2 test fifth-force wymaga zewnętrznych danych MICROSCOPE | brak | publiczne dane |
| B3 | Reform L_mat może wymagać re-run M9.x (S02 dependency) | [[../S02_volume_element_M9]] | rekomendacja S04 |

## Closed needs

| ID | Co | Domknięte przez | Data |
|----|----|-----------------|------|
| (puste) | | | |
