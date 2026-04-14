# R3: Dlaczego N=3 generacji? (T-OP3)

## Problem

**NAJFUNDAMENTALNIEJSZE otwarte pytanie TGP.**

GL(3,𝔽₂) z |GL|=168 **zakłada** N=3. Nie wyprowadza go z fizyki.
"Dlaczego 3 generacje?" to otwarte pytanie całej fizyki cząstek, nie tylko TGP.

## Obecny status (2026-04-14)

### ✅ NOWY WYNIK: Hipoteza auto-przestrzeni

| Element | Status | Wynik |
|---------|--------|-------|
| Singularność metryczna g₀_crit = 2.206 | **ZWERYFIKOWANE** | g(r) → 0 w rdzeniu solitonu |
| τ(Koide) poniżej bariery | **ZWERYFIKOWANE** | g₀^τ = 1.729 < 2.206 ✓ |
| τ(φ²) powyżej bariery | **POTWIERDZONE** | g₀^τ(φ²) = 2.276 > 2.206 ✗ |
| 4. generacja powyżej bariery | **POTWIERDZONE** | g₀^(4) = 3.683 > 2.206 ✗ |
| Mechanizm auto-przestrzeni | **MECHANIZM** | Soliton generuje przestrzeń → g → 0 → singularność |

### ⚠️ NAPIĘCIA

| Element | Problem |
|---------|---------|
| φ²-drabinka τ > g₀_crit | τ NIE siedzi na prostej drabince φ²! |
| N=2 na φ-drabince | Tylko e, μ mieszczą się pod barierą z φ-drabiną |
| N=3 z Koide | Z wiązaniem K=2/3 trzy generacje poniżej bariery |

## Kluczowy wynik: Hipoteza auto-przestrzeni

### Mechanizm

```
1. Soliton z g₀ > 1 ma ROZCIĄGNIĘTY rdzeń: g_ij = g₀·δ_ij
2. Oscylacja ODE kołysze g(r) PONIŻEJ vacuum (g < 1)
3. Głębokość dołka rośnie z g₀ (nieliniowo)
4. Przy g₀ = g₀_crit ≈ 2.206: dołek sięga g = 0
5. g = 0 to SINGULARNOŚĆ METRYKI: przestrzeń zanika
6. Solitony z g₀ > g₀_crit NIE MOGĄ ISTNIEĆ

Interpretacja fizyczna:
"Materia generuje przestrzeń — także WEWNĄTRZ siebie.
 Cięższe cząstki rozciągają rdzeń, który oscylując
 kurczy się do zera — tworzy DZIURĘ w przestrzeni.
 To jest fizyczna granica istnienia cząstki."
```

### Profil solitonu przy krytyczności

```
g₀ = 1.5:  g_min = 0.851  (głębokość 0.65)  — bezpieczny
g₀ = 1.8:  g_min = 0.700  (głębokość 1.10)  — bezpieczny
g₀ = 2.0:  g_min = 0.543  (głębokość 1.46)  — bezpieczny
g₀ = 2.1:  g_min = 0.422  (głębokość 1.68)  — bezpieczny
g₀ = 2.2:  g_min = 0.147  (głębokość 2.05)  — MARGINALNY
g₀ = 2.21: g_min → 0      (głębokość 2.21)  — SINGULARNOŚĆ!
```

### Implikacja dla τ

```
φ²-drabinka:     g₀^τ = φ²·g₀^e = 2.276 > g₀_crit = 2.206  ✗
Koide-wiązanie:  g₀^τ = 1.729           < g₀_crit = 2.206  ✓

Wniosek: τ NIE siedzi na prostej drabince φ²!
Masa τ jest OGRANICZONA BARIERĄ, a jej dokładna wartość
wynika z wiązania Koide K(m_e, m_μ, m_τ) = 2/3.

To wyjaśnia też dlaczego φ²-drabinka daje złe wyniki
w analizie R6 (A_tail = 0 dla g₀ = 2.28).
```

### Ile generacji?

```
Z drabinką φ:                    Z Koide:
  e:   g₀ = 0.869 < 2.206 ✓      e:   0.869 < 2.206 ✓
  μ:   g₀ = 1.407 < 2.206 ✓      μ:   1.407 < 2.206 ✓
  τ:   g₀ = 2.276 > 2.206 ✗      τ:   1.729 < 2.206 ✓
  4th: g₀ = 3.683 > 2.206 ✗      4th: impossible   ✗

N = 2 (prosta drabinka)           N = 3 (z wiązaniem Koide)
```

## Obecne heurystyki (żadna nie jest dowodem)

| Argument | Status | Problem |
|----------|--------|---------|
| **AUTO-PRZESTRZEŃ (NOWE)** | **MECHANIZM** | g₀_crit = 2.206, N=3 z Koide |
| d=3 → k=4 → WKB: 3 stany związane | HEURYSTYKA | "Bariera duchowa" — brak formalnego dowodu |
| |GL(3,𝔽₂)| = 168 = (2N+1)·2^N·N! | TAUTOLOGIA | Zakłada N=3 |
| N_ν = 2.984 ± 0.008 (LEP) | EKSPERYMENT | Potwierdza 3, nie wyjaśnia |
| 4. generacja dynamicznie zakazana | NUMERYCZNE | H8: PASS |
| Anomaly cancellation: N_gen = N_color | POWIĄZANIE | Nie dowodzi |

## Ścieżki ataku

### Ścieżka 6 (NOWA, najpilniejsza): Auto-przestrzeń → N=3
- Analityczne g₀_crit z V(g) = g³/3 - g⁴/4 - 1/12
- Sprawdzić: czy g₀_crit = 2φ-1 = √5? (bliskie: 2.236 vs 2.206, Δ = 1.4%)
- Dlaczego N=3 a nie N=2? → bo τ jest PONIŻEJ bariery (wiązanie Koide)
- Sformalizować: g₀_crit jest cechą potencjału, nie parametrem

### Ścieżka 5: Z solitonowego WKB
- WKB: liczba stanów związanych N_bound = ∫√(2|V|) dr / π
- Sprawdzić: czy N_bound = 3 **dokładnie** dla α=2, d=3

### Ścieżka 2: Dynamiczna (stabilność)
- g₀^(4) = φ³ · g₀^e ≈ 3.68 > g₀_crit ✓ (potwierdzone)
- NOWE: g₀^τ(φ²) = 2.28 TEŻ > g₀_crit! → φ²-drabinka nie działa dla τ

### Ścieżka 4: Algebraiczna (eliminacja)
- GL(2,𝔽₂) = S₃ (6 elementów): zbyt mało
- GL(4,𝔽₂) (20160 elementów): daje 4 generacje — zakazane przez barierę

## Pliki

| Plik | Opis | Status |
|------|------|--------|
| `r3_metric_singularity.py` | g₀_crit = 2.206, bariera metryczna | ✅ NOWE |
| `r3_self_space_stability.py` | Krajobraz A_tail(g₀), hipoteza auto-przestrzeni | ✅ NOWE |

## Kryterium zamknięcia

Twierdzenie: "W teorii solitonów z K=Φ², d=3, istnieją dokładnie 3 stabilne
sektory masowe z powodu singularności metrycznej przy g₀_crit i wiązania Koide."

## Status

- [x] Singularność metryczna g₀_crit = 2.206 — ZWERYFIKOWANE
- [x] τ(Koide) < g₀_crit < τ(φ²) — POTWIERDZONE
- [x] 4. generacja powyżej bariery — POTWIERDZONE
- [x] Mechanizm fizyczny: auto-generowana przestrzeń — SFORMUŁOWANE
- [ ] Analityczne g₀_crit z potencjału V(g)
- [ ] Dlaczego dokładnie 3 z Koide + bariera?
- [ ] Formalizacja dowodu
