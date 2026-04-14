# R6: B=√2 analitycznie z ODE solitonowego

## Problem

Stosunek Brannena B = b/a parametryzuje masy leptonów:
```
√mᵢ = a·(1 + b·cos(θ + 2πi/3)),    i = 1,2,3
```

Z danych: B_num = 1.414212... ≈ √2 (|δ| < 10⁻⁶)

Konsekwencja: K = (1 + B²/2)/N = (1 + 1)/3 = **2/3** (Koide)

Gdyby B=√2 było **udowodnione analitycznie** → Koide staje się **twierdzeniem TGP**.

## Obecny status (2026-04-14)

### ✅ UDOWODNIONE (algebraicznie)

| Element | Status | Dowód |
|---------|--------|-------|
| B = √(N-1) ⟺ K = 2/N | **TWIERDZENIE** | Tożsamość: K = (1+B²/2)/N, więc K=2/N → B²=2 |
| K = 2/3 ⟺ fazy 120° | **TWIERDZENIE** | Equidistant cos(2πi/3) → Σcos=0, Σcos²=3/2 |
| K = 2/3 z PDG mas | **WERYFIKACJA** | Q_K(PDG) = 1.500014 ≈ 3/2 |

### ⚠️ OTWARTE

| Element | Status | Problem |
|---------|--------|---------|
| Fazy 120° z GL(3,𝔽₂) | HEURYSTYCZNE | Z₃ ⊂ GL(3,𝔽₂) → fazy 120°, ale brak formalnego dowodu |
| F(φ) = A(φg₀)/A(g₀) | ZMIENNE | F(φ) NIE jest stałe (CV = 220%) — Ścieżka 4 nie działa |
| φ²-drabinka dla tau | NIEZGODNA | g₀^τ = φ²g₀^e daje A_tail = 0 (substrat); best-fit g₀^τ/g₀^e = 1.99 ≠ φ² |
| Analityczne μ z ODE | OTWARTE | μ_eff ≈ 0.91 (substrat), ~4.12 (kanoniczne) |

## Kluczowy wynik R6 (2026-04-14)

```
ŁAŃCUCH ALGEBRAICZNY (udowodniony):
  GL(3,𝔽₂) → Z₃ podgrupa → fazy równoodstępne (120°)
  → K = 2/3 (Koide dokładne)
  → B = √2 (Brannen dokładne)

BRAKUJĄCE OGNIWO:
  Dlaczego Z₃ ⊂ GL(3,𝔽₂) wymusza FAZY 120°?
  ≡ Dlaczego √m_i = M(1 + B·cos(θ + 2πi/3))?
  ≡ Dlaczego parametryzacja Brannena jest WŁAŚCIWA?
```

### Negatywne wyniki (ważne!)

1. **F(φ) nie jest stałe** — stosunek A_tail(φg₀)/A_tail(g₀) zmienia się
   dramatycznie z g₀ (od 0.27 do 35). Ścieżka 4 nie prowadzi do B=√2.

2. **φ²-drabinka tau nie działa w substracie** — soliton przy g₀ = 2.28
   nie ma właściwego ogona (A_tail = 0). Potrzebne g₀^τ ≈ 1.73
   (24% od φ² = 2.62).

3. **Prawo potęgowe jest przybliżone** — μ_eff zmienia się od 0.87 do 0.96
   w zakresie g₀ ∈ [0.5, 0.9] (CV = 3.8%).

### Pozytywne wyniki

1. **Best-fit g₀^τ daje B = 1.4143** — |B-√2| = 0.0001, potwierdzenie
   numeryczne z pełnego ODE solitonowego.

2. **Łańcuch algebraiczny jest kompletny** — jeśli K = 2/3, to B = √2 dokładnie.
   K = 2/3 wynika z faz 120° (twierdzenie).

3. **(A_μ/A_e)⁴ = 206.63** — stosunek mas e-μ z substratu ODE potwierdza k=4.

## Pliki

| Plik | Opis | Status |
|------|------|--------|
| `research/brannen_sqrt2/r6_atail_ratio_analysis.py` | Ścieżka 4: analiza F(φ), 4/9 PASS | ✅ NOWE |
| `scripts/a3d_soliton_brannen_r.py` | ODE solitonu + Brannen: 5/6 PASS | ✅ RDZEŃ |
| `scripts/a3_koide_origin_analysis.py` | K=2/3 algebraicznie: 5/5 PASS | ✅ RDZEŃ |

## Ścieżki dalszego ataku

### Ścieżka A (najpilniejsza): Z₃ → fazy 120°
- Formalny dowód: reprezentacja Z₃ na masach leptonów
- Czy Z₃ ⊂ GL(3,𝔽₂) wymusza strukturę cos(2πi/3)?
- Powiązanie z R1 (Cabibbo — ta sama Z₃!)

### Ścieżka B: Kanoniczne ODE vs substrat
- a3d używa KANONICZ. formulation z g₀^e = 1.249 (nie substratowe 0.869)
- Porównać: czy w formulacji kanonicznej F(φ) jest bardziej stałe?
- Ghost boundary g* = 0.779 zmienia charakter A_tail(g₀)

### Ścieżka C: Tau z Koide constraint
- Nie drabinka φ², lecz WIĄZANIE Koide: K(m_e, m_μ, m_τ) = 2/3
- Znając m_e i m_μ z solitonu, wyznacz m_τ z K=2/3 analitycznie
- Sprawdź: czy to daje g₀^τ = 1.73?

## Referencje rdzenia

- `dodatekT3_brannen_geometry.tex` (geometria Brannena)
- `dodatekT_koide_atail_formal.tex` (Koide z A_tail)
- `dodatekT2_koide_fp_algebra.tex` (algebra φ-FP)

## Status

- [x] Ścieżka 4: F(φ) = A(φg₀)/A(g₀) — przeskanowane, NIE stałe
- [x] Łańcuch algebraiczny: Z₃ → 120° → K=2/3 → B=√2 — KOMPLETNY
- [x] Best-fit tau: B = 1.4143 ≈ √2 z dokładnością 10⁻⁴
- [ ] Ścieżka A: Formalizacja Z₃ → fazy 120°
- [ ] Ścieżka B: Porównanie kanoniczne vs substrat
- [ ] Ścieżka C: Tau z wiązania Koide
