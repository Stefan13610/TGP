# R1: Korekcja Cabibbo (Ω_Λ/N)²

## Problem

Kąt Cabibbo — największe pojedyncze napięcie w TGP:

```
TGP (zerowy rząd):  λ_C = Ω_Λ/N = 0.6847/3 = 0.22823
PDG:                sin(θ_C) = 0.22483 ± 0.00005
Napięcie:           4.8σ
```

Kill criterion K2 przeżywa (< 30% off), ale 4.8σ wymaga wyjaśnienia.

## Obecna derywacja

- Mieszanie międzygeneracyjne CKM mediowane przez Z₃ ⊂ GL(3,𝔽₂)
- Amplitude ∝ Ω_Λ/N (tree-level)
- Brak korekcji wyższego rzędu

## Plan ataku

### Krok 1: Struktura algebraiczna CKM z GL(3,𝔽₂)
- GL(3,𝔽₂) ma 168 elementów, Z₃ podgrupę cykliczną
- Macierz CKM: V_ij = ⟨i|Z₃|j⟩ — obliczyć elementy macierzowe
- Sprawdzić: V_us = Ω_Λ/N + δV_us, gdzie δV_us ~ (Ω_Λ/N)²

### Krok 2: Sub-leading corrections
- Oblicz δλ_C = c₂·(Ω_Λ/N)² + c₃·(Ω_Λ/N)³
- Ustal c₂ z struktury grupy (nie fitowaniem!)
- Cel: |λ_C(TGP) − λ_C(PDG)| < 2σ

### Krok 3: Konsekwencje
- Nowe predykcje: V_cb, V_ub z korekcjami
- Sprawdzenie unitarności CKM

## Kryterium zamknięcia

λ_C(TGP, z korekcją) − λ_C(PDG) < 2σ

## Pliki do scalenia z rdzeniem (po zamknięciu)

- `tgp_companion.tex` §F2: nowy paragraf z korekcją
- `scripts/cabibbo_correction_verify.py`: skrypt weryfikacyjny

## Referencje rdzenia

- `tgp_companion.tex`, linie 329–342, 791–792
- `tgp_letter.tex`, linia 135
- `nbody/examples/ex247_*.py` (obecna walidacja λ_C)
- `nbody/examples/ex249_*.py` (δ_CKM, β_UT)

## Status

- [ ] Krok 1: Algebraiczna struktura CKM z GL(3,𝔽₂)
- [ ] Krok 2: Obliczenie (Ω_Λ/N)² korekcji
- [ ] Krok 3: Weryfikacja numeryczna
- [ ] Krok 4: Scalenie z rdzeniem
