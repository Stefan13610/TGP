---
title: "R7: UV Completion — unifikacja sprzężeń przy M_Pl"
date: 2026-05-03
tgp_status:
  folder_status: "needs-bridge"
  level: L1
  kind: derivation
  core_compatibility: "unknown"
  last_reviewed_against_core: "unknown"
  may_edit_core: false
  exports_findings: true
  has_needs_file: true
  has_findings_file: true
  open_bridges: []
  depends_on: []
  impacts: []
  source_of_status:
    - "README.md H1: 'R7: UV Completion — unifikacja sprzężeń przy M_Pl'"
  promoted_to_core: null
  polluted_74394a8: false
  pre_existing_findings: false
  pre_existing_needs: false
  last_yaml_update: "2026-05-03"
---

# R7: UV Completion — unifikacja sprzężeń przy M_Pl

## Problem

TGP nie ma explicit UV completion. Wpływ dynamicznego Φ na running
stałych sprzężenia α₁, α₂, α₃ nie jest obliczony.

## Obecny status

- F12: β(g₀) = 0 at 1-loop (konforemność) — ZAMKNIĘTE
- Ale: Modyfikacja propagatorów przez K(Φ) zmienia β-functions
- Unifikacja α₁ = α₂ = α₃ przy M_Pl: NIE SPRAWDZONE

## Plan ataku

### Krok 1: β-funkcje z K(Φ)
- Standardowe β-funkcje SM: β_i = b_i·α_i²/(2π)
- Modyfikacja TGP: K(Φ) = Φ² zmienia kinetyczny coupling → threshold corrections
- Obliczenie: δβ_i od dynamicznego Φ

### Krok 2: Running do M_Pl
- Rozwiązanie RGE: dα_i/d(ln μ) = β_i(α; Φ)
- Sprawdzenie: czy krzywe się spotykają

### Krok 3: Porównanie
- SM: brak unifikacji (miss ~2–3% przy M_GUT)
- MSSM: unifikacja przy ~10^16 GeV
- TGP: sprawdzić

## Kryterium zamknięcia

Wykres α_i(μ) z modyfikacjami TGP. Odpowiedź: czy unifikacja zachodzi.

## Priorytet

**BONUSOWY** — nie blokuje publikacji. Ale pozytywny wynik znacząco
wzmocniłby teorię.

## Referencje rdzenia

- `sek09_cechowanie.tex` (sector cechowania)
- `dodatekV_su3_formalizacja.tex`
- `dodatekU_su2_formalizacja.tex`
- `dodatekO_u1_formalizacja.tex`

## Status

- [ ] Krok 1: Obliczenie δβ_i
- [ ] Krok 2: Numeryczne running
- [ ] Krok 3: Analiza wyników
