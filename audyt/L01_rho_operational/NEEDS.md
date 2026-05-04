---
title: "NEEDS — L01 ρ operational definition"
date: 2026-05-04
parent: "[[README.md]]"
type: needs
tgp_owner: audyt/L01_rho_operational
tags:
  - needs
  - rho
  - ontology
---

# NEEDS — L01 ρ operational definition

## Otwarte luki

| ID | Luka | Typ | Source | Kandydat dostawcy |
|----|------|-----|--------|-------------------|
| N1 | Formal `ρ = -T^μ_μ/c²` w sek08a + mapping na sektory SM | derivation | `sek08a` lin. 79–83 (brak), audit § C.3 | nowy cykl `op-rho-stress-energy-bridge` Phase 1 |
| N2 | Treatment fotonów (T^μ_μ = 0) — czy `ρ_EM = 0` znaczy brak sprzężenia z Φ czy tylko brak skalarnego sprzęgania | analytical-bridge | brak | cykl Phase 1 |
| N3 | Trace anomaly w QFT na curved background — stabilność `ρ` pod renormalizacją | derivation | brak | cykl Phase 2 |
| N4 | SPARC N-body re-derivation z explicit T^μ_μ → ρ mapping | numerical | `nbody/` skrypty | cykl Phase 3 |
| N5 | Spójność `ρ` z `ax:zrodlo` w sek08_formalizm | analytical | `rem:materia-hierarchia` | cykl Phase 1 |

## Pytania otwarte

- **Q1**: Czy fotony (T^μ_μ = 0) sprzęgają z TGP przez g_eff (bez `φ·ρ`
  termu) czy w ogóle nie?
  - Source: brak explicit decision
  - Wpływa na: ψ.1, σ.1 (varying-c), GW170817 c_GW=c

- **Q2**: Dla pola Diraca z masą m: `T^μ_μ = m·ψ̄ψ`. Czy `ρ = m·|ψ_D|²/c²`
  jest poprawne czy potrzebny additional renorm?
  - Source: audit § C.3 + standard QFT
  - Wpływa na: predykcje masy leptonów, kinki (warstwa 3c)

- **Q3**: SPARC fits używają `ρ = ρ_baryon` (HI + stars + bulge) — czy
  jest to non-relativistic limit `T^μ_μ = -ρ_rest·c²` poprawnie zachowany?
  - Source: nbody scripts
  - Wpływa na: SPARC galaxy rotation predictions

## Blokery

| ID | Bloker | Czeka na | Source |
|----|--------|----------|--------|
| B1 | Cykl L01 musi iść równolegle z S04 (formal derivation L_mat) | [[../S04_metric_coupling_axiom]] | rekomendacja L01 |
| B2 | Phase 3 SPARC re-derivation może dotknąć fits, które są historycznie LOCKED | autor decision | nbody scripts |

## Closed needs

| ID | Co | Domknięte przez | Data |
|----|----|-----------------|------|
| (puste) | | | |
