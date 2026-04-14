# R6: B=√2 analitycznie z ODE solitonowego

## Problem

Stosunek Brannena B = b/a parametryzuje masy leptonów:
```
√mᵢ = a·(1 + b·cos(θ + 2πi/3)),    i = 1,2,3
```

Z danych: B_num = 1.414212... ≈ √2 (|δ| < 10⁻⁶)

Konsekwencja: K = (1 + B²/2)/N = (1 + 1)/3 = **2/3** (Koide)

Gdyby B=√2 było **udowodnione analitycznie** → Koide staje się **twierdzeniem TGP**.

## Obecny status

- a3d_soliton_brannen_r.py: 5/6 PASS, T5 = EXPECTED FAIL
- B wynika z profili solitonowych: g(r; g₀^e), g(r; g₀^μ), g(r; g₀^τ)
- g₀^μ = φ·g₀^e, g₀^τ = φ²·g₀^e (φ-FP scaling)
- √m ∝ A_tail² → B = b/a zależy od stosunków A_tail

## Plan ataku

### Ścieżka 1: Analiza ODE (główna)
- ODE solitonowy: g'' + (2/r)g' = V'(g), g(0)=g₀, g'(0)=0, g(∞)=1
- z K(g) = g^{2α}, α=2
- Rozwinięcie asymptotyczne: g(r) ~ 1 - A(g₀)·e^{-μr}/r
- Obliczenie A(g₀) analitycznie (metoda dopasowania asymptotyk)
- Cel: A(φ·g₀)/A(g₀) = f(φ) → B = F(f(φ))

### Ścieżka 2: Symetria ODE
- Analiza Lie punktowych symetrii ODE: g'' + (2/r)g' - V'(g) = 0
- Czy ODE ma ukrytą symetrię dyskretną wymuszającą B=√2?
- Grupa symetrii + warunki brzegowe → wiązanie na B

### Ścieżka 3: Perturbacyjna
- Niech B = √2 + ε
- Ekspandować masy: m_i(ε) = m_i(0) + ε·δm_i + ...
- Sprawdzić: czy ε=0 jest **punktem stacjonarnym** jakiegoś funkcjonału
- Np.: ∂/∂B [Σ√m / (Σm)^{1/2}] = 0 przy B=√2?

### Ścieżka 4: Z φ-FP dokładnie
- φ-FP: g₀^(n+1) = φ·g₀^(n) → A^(n+1) = F(φ)·A^(n) + korekcje
- Jeśli F(φ) jest dokładnie obliczalne → B wynika algebraicznie
- Klucz: czy A_tail(φ·g₀)/A_tail(g₀) jest stałe (niezależne od g₀)?

## Kryterium zamknięcia

**Twierdzenie:** "Dla ODE solitonowego z α=2, d=3, φ-FP: B = √2 dokładnie"

## Pliki do scalenia z rdzeniem

- Rozszerzenie `dodatekT3_brannen_geometry.tex`
- Nowy akapit w `tgp_companion.tex` (sekcja Koide)

## Referencje rdzenia

- `dodatekT3_brannen_geometry.tex` (geometria Brannena)
- `dodatekT_koide_atail_formal.tex` (Koide z A_tail)
- `dodatekT2_koide_fp_algebra.tex` (algebra φ-FP)
- `scripts/a3d_soliton_brannen_r.py` (5/6 PASS)
- `scripts/a3_koide_origin_analysis.py` (5/5 PASS)

## Status

- [ ] Ścieżka 1: Analityczne A(g₀) z ODE
- [ ] Ścieżka 4: F(φ) = A(φ·g₀)/A(g₀) — obliczenie
- [ ] Ścieżka 2: Analiza Lie symetrii
- [ ] Ścieżka 3: ε-ekspansja wokół B=√2
- [ ] Połączenie w dowód B=√2
