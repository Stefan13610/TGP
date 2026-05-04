---
title: "NEEDS — D01 drifting parameters"
date: 2026-05-04
parent: "[[README.md]]"
type: needs
tgp_owner: audyt/D01_drifting_numbers
tags:
  - needs
  - drift
  - anchor-lock
  - propagation
---

# NEEDS — D01 drifting parameters

## Otwarte luki

| ID | Luka | Typ | Source | Kandydat dostawcy |
|----|------|-----|--------|-------------------|
| N1 | α_s = 0.1184 globalna propagacja przez 15+ plików (LP-5, LP-7, sek00 sync, companion) | editorial | audit § B.3 (B3-v2 pending) | nowy cykl `op-anchor-lock-propagation` Phase 2 |
| N2 | Φ₀ = 24.783 globalna propagacja przez 15+ plików (C10-v2 pending) | editorial | audit § C.10 | cykl Phase 2 |
| N3 | m_H sync sek09, dodatekU (TGP-side 125.31; reference 125.20 PDG 2024) | editorial | audit § C.9 | cykl Phase 2 |
| N4 | Σm_ν = 59.01 meV propagacja do README key predictions table | editorial | audit § B.4 | cykl Phase 2 |
| N5 | g₀^e = 0.86941 lock po L04 decyzji + propagacja | editorial + decision | LP-6 N-1 | po [[../L04_ODE_dualism_alpha]] |
| N6 | k vs p = 4 vs 3 unification po L04+L05 | editorial | LP-4 + R3 CORRECTIONS | po [[../L05_mass_exponent_drift]] |
| N7 | Re-run wszystkich verification scripts z locked wartościami | numerical | brak | cykl Phase 3 |
| N8 | Re-derive predictivity ratio dla *jednej* trajektorii | analytical | LS-8 audit (was for selected trajectory) | cykl Phase 3 |

## Pytania otwarte

- **Q1**: Czy `m_H` TGP-side = 125.31 (z formal `v × 57/112`) jest *predykcją*
  czy *kalibracją*? Jeśli predykcja, dlaczego deklaruje wartość 0.7σ
  od PDG zamiast PDG-target?
  - Source: README claim
  - Wpływa na: status m_H jako DERIVED

- **Q2**: Czy spread Φ₀ (24.66 / 24.783 / 25.0) reprezentuje fizyczną
  niepewność (np. systematic w F11 mass formula) czy artefakt
  niedokończonej propagacji?
  - Source: ROADMAP_v3:65, 122, 896
  - Wpływa na: status Φ₀ jako kalibracja vs derived

- **Q3**: Czy LP-5 (α_s=0.1174) i LP-7 (α_s=0.1171) używają tego samego
  Φ₀ co LP-1/sek00? Jeśli nie, drift α_s jest artefaktem różnych Φ₀.
  - Source: LP-5, LP-7 vs sek00
  - Wpływa na: czy drift α_s redukuje się do drift Φ₀

## Blokery

| ID | Bloker | Czeka na | Source |
|----|--------|----------|--------|
| B1 | g₀^e lock i k/p resolution czeka na L04 | [[../L04_ODE_dualism_alpha]] | rekomendacja D01 |
| B2 | Phase 2 propagation potrzebuje stabilnej trajektorii Φ₀ → α_s → m_W → m_H | rekomendacja Phase 1 | audit § C.10 |
| B3 | Re-derivation predictivity ratio może obniżyć liczbę z 5.5 do 4–4.5 | LS-8 re-run | rekomendacja Phase 3 |

## Closed needs

| ID | Co | Domknięte przez | Data |
|----|----|-----------------|------|
| (puste) | | | |
