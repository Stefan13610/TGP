---
title: "NEEDS — S03 β_PPN konwencja"
date: 2026-05-04
parent: "[[README.md]]"
type: needs
tgp_owner: audyt/S03_beta_PPN_convention
tags:
  - needs
  - PPN
  - critical
---

# NEEDS — S03 β_PPN konwencja

## Otwarte luki

| ID | Luka | Typ | Source | Kandydat dostawcy |
|----|------|-----|--------|-------------------|
| N1 | Brak Lagrangianu, z którego f(ψ) = (4-3ψ)/ψ wynika wariacyjnie | derivation | audit § B.8 (cyrkularność „potrójnej motywacji") | nowy cykl `op-M911-lagrangean-derivation` |
| N2 | c₂ = -1 derived z α=2 (forma I); brak niezależnej weryfikacji dla (IV) | numerical | audit § A.3 + § M.2 | cykl B6 M9.x re-run |
| N3 | Linearyzacja `sek08c.tex` lin. 183–191 (β=2) sprzeczna z M9.1'' β=1 | analytical-bridge | `sek08c.tex` lin. 183–191 | sek08c-rewrite (S01) |
| N4 | Brak explicit caveat w README/FOUNDATIONS że β=1 jest „master formula choice" | editorial | README.md, TGP_FOUNDATIONS.md | krótki edit |
| N5 | Niezależny eksperymentalny test c₂ (np. z dyspersji GW wyższych rzędów) | data-bridge | brak | LIGO 3G post-2030 |

## Pytania otwarte

- **Q1**: Czy przyjąć Ścieżkę A (B8 Lagrangean derivation), B (explicit
  acknowledgment), czy C (wycofanie M9.1'')?
  - Source: rekomendacja S03
  - Wpływa na: status M9.1'' jako kanonicznej, S01, S02, główny artykuł

- **Q2**: Czy „master formula" jest standardową konwencją PPN dla teorii
  z niestandardowym K(φ), czy jest specyficzna dla TGP?
  - Source: literatura PPN (Will 1993; Damour-Esposito-Farèse 1992)
  - Wpływa na: legitimacy konwencji w eyes of external referees

- **Q3**: Jeśli f(ψ) wynika z Lagrangianu, jaki jest minimalny wymóg
  na akcję `S[Φ, g]` z którego wynikają oba: (a) `g_tt = -c²(4-3ψ)/ψ`,
  (b) c₂ = -1?
  - Source: B8 audytu
  - Wpływa na: czy istnieje czysto skalarna teoria realizująca M9.1''

## Blokery

| ID | Bloker | Czeka na | Source |
|----|--------|----------|--------|
| B1 | Bez Q1 decyzji nie idzie sek08c-rewrite (S01) | autor | audit § M.2 |
| B2 | Cykl B8 wymaga niezależnego fizyka teoretyka (nie subagenta) | autor decision | rekomendacja S03 Ścieżka A |
| B3 | Wszystkie predykcje PPN (γ=1, β=1, α₁₋₃, ζ₁₋₄, ξ_PPN) wiszą na konwencji | OP-7 T6 12/12 PASS | audit § A.3 |

## Closed needs

| ID | Co | Domknięte przez | Data |
|----|----|-----------------|------|
| (puste) | | | |
