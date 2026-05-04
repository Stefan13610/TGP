---
title: "NEEDS — L04 ODE dualism α=1 vs α=2"
date: 2026-05-04
parent: "[[README.md]]"
type: needs
tgp_owner: audyt/L04_ODE_dualism_alpha
tags:
  - needs
  - ODE
  - alpha
  - decision
---

# NEEDS — L04 ODE dualism α=1 vs α=2

## Otwarte luki

| ID | Luka | Typ | Source | Kandydat dostawcy |
|----|------|-----|--------|-------------------|
| N1 | **Decyzja autorska**: K=g² (α=1) czy K=g⁴ (α=2) jako kanoniczna formulacja TGP | decision | LP-6 § N-1 | autor |
| N2 | Globalny lock wybranej formulacji w core/sek08* (~30 lokalizacji) | editorial | brak | po decyzji N1 |
| N3 | Sync `m_obs = c·A^(5-α)` formula (R3) z decyzją — czy k=4 czy p=3? | analytical | `why_n3/CORRECTIONS_2026-05-01.md` | po decyzji N1 |
| N4 | README claim „α=2 algebraic theorem" — albo prawdziwie używać K=g⁴, albo wycofać | editorial | README lin. ~102 | po decyzji N1 |
| N5 | Druga formulacja (po wyborze) → dodatek historyczny z explicit „equivalent in weak field, different in strong field" | editorial | brak | po decyzji N1 |
| N6 | Quark sector: K=g² daje cross-sector Koide ✓ — czy K=g⁴ to też daje? | numerical | LP-5 13/13 PASS | po decyzji N1 lub w cyklu L04 Phase 2 |

## Pytania otwarte

- **Q1**: Jakie są **fizyczne** argumenty za α=2 (poza „algebraic theorem")?
  - Source: README, FOUNDATIONS §3 (twierdzenie alfa2 w core/sek08b)
  - Wpływa na: wybór

- **Q2**: Jeśli K=g² jest preferowana (ghost-free + stabilna + reprodukuje
  τ), dlaczego README forsuje K=g⁴?
  - Source: README vs LP-6 § N-1 explicit recommendation
  - Wpływa na: integralność manuskryptu

- **Q3**: Czy istnieje *trzecia* formulacja, która łączy zalety obu?
  Np. `K(φ) = φ²(1+φ²)/2` — daje α=2 limit dla małych φ, α=4 dla
  dużych?
  - Source: open speculation
  - Wpływa na: czy decyzja jest binarna

## Blokery

| ID | Bloker | Czeka na | Source |
|----|--------|----------|--------|
| B1 | **Decyzja N1 jest fundamentalna i blokuje wszystkie inne luki w L04, L05** | autor | rekomendacja L04 |
| B2 | Bez decyzji N1, L05 (mass exponent) i L03 (spektralna analiza) nie idą | [[../L05_mass_exponent_drift]], [[../L03_K_phi_stability]] | rekomendacja L04 |
| B3 | D01 (drift g₀^e) zależy od wyboru formulacji ODE | [[../D01_drifting_numbers]] | rekomendacja L04 |

## Closed needs

| ID | Co | Domknięte przez | Data |
|----|----|-----------------|------|
| (puste) | | | |
