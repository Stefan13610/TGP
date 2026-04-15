# R3: Dlaczego N=3 generacji? (T-OP3)

## Problem

**NAJFUNDAMENTALNIEJSZE otwarte pytanie TGP.**

GL(3,F₂) z |GL|=168 **zakłada** N=3. Nie wyprowadza go z fizyki.
"Dlaczego 3 generacje?" to otwarte pytanie całej fizyki cząstek, nie tylko TGP.

## Obecny status (2026-04-15)

### ✅ GŁÓWNY WYNIK: N=3 z bariery metrycznej

| Element | Status | Wynik |
|---------|--------|-------|
| Singularność metryczna g₀_crit | **ZWERYFIKOWANE** | g(r) → 0 w rdzeniu solitonu |
| g₀_crit(1D) = 4/3 analitycznie | **TWIERDZENIE** | Z prawa zachowania F(g₀) = 0 |
| g₀_crit(2D) ≈ 1.7324 (≠ √3) | **NUMERYCZNE** | Δ = 3.9×10⁻⁴ od √3 |
| g₀_crit(3D) ≈ 2.2062 (substrat) | **NUMERYCZNE** | Bisection, tol < 10⁻⁸ |
| g₀_crit(3D, canonical) ≈ 1.546 | **NUMERYCZNE** | K=g⁴ daje NIŻSZĄ barierę! |
| τ(Koide) < g₀_crit(sub) | **POTWIERDZONE** | g₀^τ = 1.729 < 2.206 ✓ |
| τ(φ²) > g₀_crit(sub) | **POTWIERDZONE** | g₀^τ(φ²) = 2.276 > 2.206 ✗ |
| 4. generacja > g₀_crit | **POTWIERDZONE** | g₀^(4) = 2.354 > 2.206 ✗ |
| dm/dg₀ → ∞ przy barierze | **POTWIERDZONE** | Masa dywerguje — twardy limit |

### ✅ NOWE: N=3 BEZ KOIDE (z α + bariera)

```
Kluczowy wynik (r3_n3_from_barrier.py):

Z skalowania m ~ A^{2α} i eksperymentalnych mas:
  α = ln(m_μ/m_e) / (2·ln(A_μ/A_e)) = 2.35

  τ:   g₀ = 1 + A_e·(3477)^{1/4.7} = 1.742 < g₀_crit = 2.206  ✓
  4th: g₀ = 1 + A_e·(58477)^{1/4.7} = 2.354 > g₀_crit = 2.206  ✗

  N = 3 wynika z SAMEJ bariery + skalowania mas!
  Koide jest KONSYSTENTNE ale NIE POTRZEBNE jako input.
```

### ⚠️ NAPIĘCIA I OTWARTE PYTANIA

| Element | Problem |
|---------|---------|
| g₀_crit zależy od α | Substrate (α=1): g₀_crit=2.206 vs Canonical (α=2): g₀_crit=1.546 |
| Canonical daje N≤2 | g₀_crit(can)=1.546 — τ NIE mieści się pod barierą! |
| φ²-drabinka daje N=2 | τ NIE siedzi na prostej drabince φ² |
| Analityczne g₀_crit(3D)? | Brak zamkniętej formy. Best: 48/π² ≈ g₀_crit² (Δ=0.004) |

## Hipoteza auto-przestrzeni

### Mechanizm

```
1. Soliton z g₀ > 1 ma ROZCIĄGNIĘTY rdzeń: g_ij = g₀·δ_ij
2. Oscylacja ODE kołysze g(r) PONIŻEJ vacuum (g < 1)
3. Głębokość dołka rośnie z g₀ (nieliniowo)
4. Przy g₀ = g₀_crit: dołek sięga g = 0
5. g = 0 to SINGULARNOŚĆ METRYKI: przestrzeń zanika
6. Solitony z g₀ > g₀_crit NIE MOGĄ ISTNIEĆ

Interpretacja fizyczna:
"Materia generuje przestrzeń — także WEWNĄTRZ siebie.
 Cięższe cząstki rozciągają rdzeń, który oscylując
 kurczy się do zera — tworzy DZIURĘ w przestrzeni.
 To jest fizyczna granica istnienia cząstki."
```

### Profil solitonu przy krytyczności (substrat, d=3)

```
g₀ = 1.5:  g_min = 0.851  (depth 0.649) — bezpieczny
g₀ = 1.8:  g_min = 0.700  (depth 1.100) — bezpieczny
g₀ = 2.0:  g_min = 0.543  (depth 1.457) — bezpieczny
g₀ = 2.1:  g_min = 0.422  (depth 1.678) — bezpieczny
g₀ = 2.2:  g_min = 0.147  (depth 2.053) — MARGINALNY
g₀ ≈ 2.21: g_min → 0      (depth 2.21)  — SINGULARNOŚĆ!
```

### Dywergencja masy przy barierze

```
dm/dg₀ przy g₀ bliskim g₀_crit:
  g₀ = 2.16: dm/dg₀ =  30,768
  g₀ = 2.18: dm/dg₀ =  48,515
  g₀ = 2.20: dm/dg₀ = 101,427
  g₀ = 2.204: dm/dg₀ = 217,693  → ∞

dm/dg₀ DYWERGUJE → masa rośnie bez granic przy barierze.
To jest TWARDY LIMIT na masę cząstki.
```

## g₀_crit(d) — zależność od wymiaru

### Twierdzenie (1D, dokładne)

```
ODE: g'' + (1/g)(g')² = 1 - g    (brak tłumienia)

Prawo zachowania: q = (g')² spełnia dq/dg + 2q/g = 2(1-g)
Rozwiązanie:      g²(g')² = 2g³/3 - g⁴/2 + C

BC: q(g₀) = 0 → C = g₀⁴/2 - 2g₀³/3
g_min = 0 iff F(g₀) = 0 → g₀ = 4/3  ■

Dowód: F(x) = 2x³/3 - x⁴/2 ma zera przy x = 0 i x = 4/3.
```

### Wartości numeryczne (substrat)

```
d=1: g₀_crit = 1.3333  = 4/3 (dokładne)
d=2: g₀_crit = 1.7324  (≠ √3 = 1.7321, Δ = 3.9×10⁻⁴)
d=3: g₀_crit = 2.2062  (przypadek fizyczny)
d=4: g₀_crit = 2.7645
d=5: g₀_crit = 3.4187
d=6: g₀_crit = 4.1811

Best fit: g₀_crit ≈ 0.218·d^{1.46} + 1.12  (max err 1.6%)
```

### g₀_crit ZALEŻY OD α! (nowe)

```
ODE: g'' + (2α-1)/g · (g')² + ((d-1)/r)g' = (1-g)/g²

d=3:
  α=1 (substrat, K=g²): g₀_crit = 2.206
  α=2 (kanoniczny, K=g⁴): g₀_crit = 1.546

Wyższe α → silniejsze sprzężenie → NIŻSZA bariera!
To fundamentalnie wpływa na liczbę generacji.
```

## Ile generacji? (analiza wg ODE)

### Substrat (α=1): N=3 ✓

```
φ-drabinka (g₀ spacing):
  e:   g₀ = 0.869 < 2.206 ✓
  μ:   g₀ = 1.407 < 2.206 ✓
  τ:   g₀ = 2.276 > 2.206 ✗  → N = 2 (z φ-drabiną)

Mass scaling α=2.35 (z experiment):
  e:   g₀ = 0.869 < 2.206 ✓
  μ:   g₀ = 1.407 < 2.206 ✓
  τ:   g₀ = 1.742 < 2.206 ✓
  4th: g₀ = 2.354 > 2.206 ✗  → N = 3 ✓

Koide (K=2/3):
  e:   g₀ = 0.869 < 2.206 ✓
  μ:   g₀ = 1.407 < 2.206 ✓
  τ:   g₀ = 1.729 < 2.206 ✓
  4th: impossible             → N = 3 ✓
```

### Kanoniczny (α=2): N≤2 ⚠️

```
g₀_crit(can) = 1.546
  e:   g₀ = 0.869 < 1.546 ✓
  μ:   g₀ = 1.407 < 1.546 ✓  (margin tylko 10%!)
  τ:   g₀ > 1.546            ✗ → N ≤ 2

NAPIĘCIE: kanoniczny ODE ma za niską barierę dla 3 generacji.
Możliwe rozwiązania:
  (a) Fizyczny ODE to substrat, nie kanoniczny
  (b) Wartości g₀ różnią się między ODE (ten sam fizyczny lepton
      ma INNY g₀ w substracie vs kanonicznym)
  (c) Poprawka do prostego skalowania A^{2α}
```

## Obecne heurystyki (żadna nie jest dowodem)

| Argument | Status | Uwagi |
|----------|--------|-------|
| **BARIERA + SKALOWANIE** | **MECHANIZM** | α=2.35, N=3 bez Koide |
| **AUTO-PRZESTRZEŃ** | **MECHANIZM** | g₀_crit z singularności metryki |
| **1D TWIERDZENIE** | **DOWÓD** | g₀_crit(1D) = 4/3 z prawa zachowania |
| d=3 → k=4 → WKB: 3 stany | HEURYSTYKA | Brak dowodu formalnego |
| |GL(3,F₂)| = 168 | TAUTOLOGIA | Zakłada N=3 |
| N_ν = 2.984 ± 0.008 (LEP) | EKSPERYMENT | Potwierdza 3, nie wyjaśnia |
| 4. generacja zakazana dynamicznie | NUMERYCZNE | H8: PASS |
| Anomaly cancellation | POWIĄZANIE | N_gen = N_color |

## Ścieżki ataku

### Ścieżka 6 (najpilniejsza): Auto-przestrzeń → N=3
- ✅ Analityczne g₀_crit(1D) = 4/3 — DONE
- ✅ g₀_crit(2D) ≈ 1.7324 (nie √3) — DONE
- ✅ g₀_crit zależy od α — DONE
- ✅ N=3 z bariery + α=2.35 — DONE
- ⬜ Rozwiązać napięcie substrat vs kanoniczny
- ⬜ Analityczne g₀_crit(3D) — brak zamkniętej formy

### Ścieżka 5: Z solitonowego WKB
- ⬜ WKB: N_bound = ∫√(2|V|) dr / π
- Trudność: solitony mają ~63 węzłów niezależnie od g₀ (sekcja 6 skryptu)
- Zliczanie węzłów NIE daje generacji

### Ścieżka 2: Dynamiczna (stabilność)
- ✅ g₀^(4) > g₀_crit (dla obu ODE)
- ✅ dm/dg₀ → ∞ przy barierze

### Ścieżka 4: Algebraiczna (eliminacja)
- GL(2,F₂) = S₃ (6 el.): zbyt mało
- GL(4,F₂) (20160 el.): daje 4 generacje — zakazane

## Pliki

| Plik | Opis | Status |
|------|------|--------|
| `r3_metric_singularity.py` | g₀_crit = 2.206, bariera metryczna | ✅ |
| `r3_self_space_stability.py` | Krajobraz A_tail(g₀), hipoteza auto-przestrzeni | ✅ |
| `r3_g0crit_analytical.py` | g₀_crit(1D) = 4/3 exact, dimension scan | ✅ |
| `r3_n3_from_barrier.py` | **N=3 z bariery + α, canonical vs substrate** | ✅ NOWE |

## Kryterium zamknięcia

Twierdzenie: "W teorii solitonów z K=Φ², d=3, istnieją dokładnie 3 stabilne
sektory masowe z powodu singularności metrycznej przy g₀_crit i skalowania mas."

Status: **CZĘŚCIOWO UDOWODNIONE** — mechanizm działa, ale zależy od wyboru K(g).

## Checklist

- [x] Singularność metryczna g₀_crit — ZWERYFIKOWANE
- [x] g₀_crit(1D) = 4/3 — TWIERDZENIE
- [x] g₀_crit(2D) ≈ 1.7324 (≠ √3) — ZWERYFIKOWANE
- [x] g₀_crit zależy od α — POTWIERDZONE
- [x] τ(Koide) < g₀_crit(sub) — POTWIERDZONE
- [x] 4. generacja > g₀_crit — POTWIERDZONE
- [x] dm/dg₀ → ∞ przy barierze — POTWIERDZONE
- [x] N=3 z α=2.35 + bariera (bez Koide) — NOWE
- [ ] Rozwiązać napięcie substrat vs kanoniczny
- [ ] Analityczne g₀_crit(3D)
- [ ] Wyprowadzić Koide z teorii solitonów
- [ ] Formalizacja dowodu
