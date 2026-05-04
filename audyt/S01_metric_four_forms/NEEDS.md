---
title: "NEEDS — S01 cztery sprzeczne formy metryki"
date: 2026-05-04
parent: "[[README.md]]"
type: needs
tgp_owner: audyt/S01_metric_four_forms
tags:
  - needs
  - metric
  - sek08c
  - critical
---

# NEEDS — S01 cztery sprzeczne formy metryki

## Otwarte luki

| ID | Luka | Typ | Source | Kandydat dostawcy |
|----|------|-----|--------|-------------------|
| N1 | sek08c nie zawiera kanonicznej derywacji M9.1'' (forma IV) — tylko nagłówkowy komentarz | analytical-bridge | `core/sek08c_metryka_z_substratu/sek08c_metryka_z_substratu.tex` lin. 64–129 | `research/op-newton-momentum/M9_1_pp_*` (źródło derywacji) |
| N2 | `thm:antipodal-uniqueness` deklaruje jedyność, ale jest błędne dla O(ε²) (β_metric=3) | lemma | `sek08c.tex` lin. 345–371 | autor (decyzja: wycofać czy uwarunkować) |
| N3 | Linearyzacja `sek08c.tex` lin. 183–191 daje β=2 — sprzeczna z M9.1'' β=1 | analytical-bridge | `sek08c.tex` lin. 183–191 | sek08c-rewrite cycle |
| N4 | Forma (II) eksponencjalna nie ma niezależnego wyprowadzenia z budżetu | derivation | `sek08c.tex` lin. 208–211 | (do wycofania jako OBSOLETE) |
| N5 | Brak formalnego twierdzenia, że M9.1'' jest *jedyną* metryką dającą β=γ=1 dla pojedynczego pola Φ z `eq:budget-def` | uniqueness-theorem | brak (audit B8) | `op-metric-uniqueness-M911` (proponowany cykl) |

## Pytania otwarte

- **Q1**: Czy zachować formę (I) jako historyczną referencję, czy usunąć
  z `sek08c.tex` całkowicie?
  - Source: audit § M.2 zaleca markery inline; rekomendacja S01 sugeruje
    rewrite z dodatkiem historycznym.
  - Wpływa na: czytelność manuskryptu, główny `main.tex` compile.

- **Q2**: Czy `eq:budget-def` (`B = N_B · s₀`) faktycznie wymusza M9.1''
  (forma IV), czy pozwala na inne hiperboliczne formy `g_tt = -c²·F(ψ)/ψ`
  z `F(ψ_0)=1, F'(1) = -1, F'(4/3)=0`?
  - Source: brak formalnego twierdzenia jedyności w sek08c
  - Wpływa na: legitimacy formy (IV) jako „kanonicznej"

- **Q3**: Czy „master formula" β = f''(1)/f'(1)² + 2c₂/f'(1) jest
  poprawnym kryterium PPN w teoriach z niestandardowym kinetic
  prefactorem K(φ), czy jest to ad-hoc sprzęganie metryki+kinetyki?
  - Source: audit § A.3, § M.2
  - Wpływa na: status β_PPN=1 (S03)

## Blokery

| ID | Bloker | Czeka na | Source |
|----|--------|----------|--------|
| B1 | sek08c-rewrite nie może iść bez decyzji w Q1 i Q2 | autor + nowy cykl `op-metric-uniqueness-M911` | rekomendacja S01 |
| B2 | S02 (volume element re-run) wymaga klarownej formy (IV) w sek08c | N1 | [[../S02_volume_element_M9]] |
| B3 | S03 (β_PPN konwencja) wymaga decyzji Q3 | autor lub formalna derywacja Lagrangianu | [[../S03_beta_PPN_convention]] |

## Closed needs

| ID | Co | Domknięte przez | Data |
|----|----|-----------------|------|
| (puste) | | | |
