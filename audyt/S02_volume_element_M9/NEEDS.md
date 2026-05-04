---
title: "NEEDS — S02 √(-g) niespójny z M9.1''"
date: 2026-05-04
parent: "[[README.md]]"
type: needs
tgp_owner: audyt/S02_volume_element_M9
tags:
  - needs
  - volume-element
  - M9
  - critical
---

# NEEDS — S02 `√(-g)` niespójny z M9.1''

## Otwarte luki

| ID | Luka | Typ | Source | Kandydat dostawcy |
|----|------|-----|--------|-------------------|
| N1 | κ = 3/(4Φ_0) wyprowadzone z formy (I) `√(-g)=c·ψ` — ponownie obliczyć dla (IV) | numerical | `sek08a` lin. 451+ | cykl B6 M9.x re-run |
| N2 | Φ-EOM `eq:field-eq-reproduced` używa volume z (I); wyprowadzić dla (IV) | analytical-bridge | `sek08a` lin. 80–86, `sek08_formalizm` | cykl B6 |
| N3 | M9.2 m_field PRE-AUDIT — ponownie ekstrahować z linearyzacji w `√(-g)=c·ψ/(4-3ψ)` | numerical | `op-newton-momentum/M9_2_results.md` | cykl B6 |
| N4 | M9.3 kwadrupol Peters-Mathews PRE-AUDIT; re-fit GW170817 `c_GW=c` | numerical | `op-newton-momentum/M9_3_results.md` | cykl B6 |
| N5 | c₂ PPN coefficient (= -1) potwierdzony tylko z (I); re-derive z (IV) | numerical | `op-newton-momentum/M9_1_pp_P3_results.md` | cykl B6 |
| N6 | M9.2 WEP MICROSCOPE 10⁻¹⁵ — re-test z poprawną metryką | numerical | `op-newton-momentum/M9_2_results.md` § B9 audit | cykl B6 |

## Pytania otwarte

- **Q1**: Czy zmiana `√(-g) = c·ψ → c·ψ/(4-3ψ)` zachowuje znak β_PPN=1 w
  master formula, czy modyfikuje master formula (S03)?
  - Source: audit § A.2 + § A.3
  - Wpływa na: S03 (β_PPN konwencja)

- **Q2**: Dla `c_GW = c` (GW170817) czy korekta volume O(U) wpada w
  bound `< 3·10⁻¹⁵`?
  - Source: README key predictions table
  - Wpływa na: jeden z 5 falsyfikatorów TGP

- **Q3**: Czy m_field przeliczone z (IV) zachowuje `m_field/m_g ≈
  3.48·10⁻²` (wyprowadzane jako bezpieczne dla MICROSCOPE)?
  - Source: M9_2_results.md lin. 90–113
  - Wpływa na: pred WEP test

## Blokery

| ID | Bloker | Czeka na | Source |
|----|--------|----------|--------|
| B1 | Cały cykl B6 czeka na S01 zamknięcie strukturalnie | [[../S01_metric_four_forms]] N1 | rekomendacja S02 |
| B2 | Re-fit GW170817 wymaga zewnętrznych danych LIGO/Virgo (publiczne, OK) | brak | — |
| B3 | β_PPN sequencing: B6 musi być zsynchronizowane z S03 decision Q3 | [[../S03_beta_PPN_convention]] | audit § M.2 |

## Closed needs

| ID | Co | Domknięte przez | Data |
|----|----|-----------------|------|
| (puste) | | | |
