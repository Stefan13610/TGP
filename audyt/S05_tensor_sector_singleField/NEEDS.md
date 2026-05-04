---
title: "NEEDS — S05 σ_ab vs single-Φ axiom"
date: 2026-05-04
parent: "[[README.md]]"
type: needs
tgp_owner: audyt/S05_tensor_sector_singleField
tags:
  - needs
  - axiom
  - sigma_ab
  - tensor
  - critical
---

# NEEDS — S05 σ_ab vs single-Φ axiom

## Otwarte luki

| ID | Luka | Typ | Source | Kandydat dostawcy |
|----|------|-----|--------|-------------------|
| N1 | Formal derywacja effective action `Γ_eff[σ_ab]` z `S_TGP[Φ]` przez integration-out ŝ | derivation | brak (audit § E.4) | nowy cykl `op-OP7-axiom-decision` Ścieżka A |
| N2 | M² = 2 m_s² (spektralny gap σ_ab) jest *postulatem* w T3-extended Bethe-Salpeter, nie wyprowadzony z S[Φ] | derivation | T3-extended results | cykl OP7-axiom Ścieżka A |
| N3 | T3 EOM dla σ_ab — czy wynika z `δS/δŝ` po projekcji, czy traktuje σ_ab jak niezależne pole? | analytical-bridge | OP-7 T3 results 44/47 PASS | cykl OP7-axiom |
| N4 | Linearyzacja jednoskalarnej Φ-EOM daje TYLKO scalar modes — pokazać, że TT modes wychodzą z `Γ_eff[σ_ab]` | numerical | M9.3 results | cykl OP7-axiom Ścieżka A |
| N5 | Falsyfikator "no breathing mode 1000+ events" — explicit predykcja gdy σ_ab to kompozyt, nie pole | data-bridge | README falsyfikatory | LIGO Run 5+ post-2030 |
| N6 | Decyzja autorska: Ścieżka A (kompozyt), B (dwupolowa), czy C (tylko skalar)? | decision | rekomendacja S05 | autor |

## Pytania otwarte

- **Q1**: Czy formal `Γ_eff[σ_ab]` w sense Wilsonowskim jest możliwy
  do wyprowadzenia z `S[Φ]` z `Φ = ŝ²` substratem v2 GL?
  - Source: closure_2026-04-26 Path B PRIMARY claim
  - Wpływa na: cała grawitacja TGP, smoking gun w GW

- **Q2**: Jeśli Ścieżka B (dwupolowa), czy TGP nie kolapsuje do
  Horndeski/Brans–Dicke? Co odróżnia TGP fenomenologicznie od scalar-tensor?
  - Source: rem:not-scalar-tensor sek08a lin. 129–151
  - Wpływa na: unique falsifiable predictions

- **Q3**: Jeśli Ścieżka C (tylko skalar), jak wytłumaczyć obserwowane
  GW150914 h_+, h_× (LIGO measure)? Czy TGP w GR-limit *trywialnie*
  daje TT polaryzacje, czy explicit nie?
  - Source: T5 fit GW150914
  - Wpływa na: status TGP vs LIGO data

- **Q4**: Path B PRIMARY (closure_2026-04-26) deklaruje σ_ab kompozytem
  z box-of-product algebra — czy to faktycznie jest formal derivation
  M² = 2m_s² czy strukturalny argument?
  - Source: closure_2026-04-26/correction_to_OP7_T3.md
  - Wpływa na: T3-extended legitimacy

## Blokery

| ID | Bloker | Czeka na | Source |
|----|--------|----------|--------|
| B1 | Cały cykl OP7-axiom-decision wymaga decyzji autora (Ścieżka A/B/C) | autor | rekomendacja S05 |
| B2 | Ścieżka A (preferowana) wymaga zaawansowanej fizyki (functional integration-out) | dedicated theorist | brak |
| B3 | Ścieżka B wymaga przepisania FOUNDATIONS §1 (binding aksjomatyka) | autor decision + review | rekomendacja S05 |
| B4 | Status OP-7 jako CLOSED 96.9% PASS jest do reklasyfikacji jako CLOSED-CONDITIONAL na Ścieżkę | rekomendacja S05 + decyzja autora | audit § B.6 |

## Closed needs

| ID | Co | Domknięte przez | Data |
|----|----|-----------------|------|
| (puste) | | | |
