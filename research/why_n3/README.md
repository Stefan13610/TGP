# R3: Dlaczego N=3 generacji? (T-OP3)

## Problem

**NAJFUNDAMENTALNIEJSZE otwarte pytanie TGP.**

GL(3,𝔽₂) z |GL|=168 **zakłada** N=3. Nie wyprowadza go z fizyki.
"Dlaczego 3 generacje?" to otwarte pytanie całej fizyki cząstek, nie tylko TGP.

## Obecne heurystyki (żadna nie jest dowodem)

| Argument | Status | Problem |
|----------|--------|---------|
| d=3 → k=4 → WKB: 3 stany związane | HEURYSTYKA | "Bariera duchowa" — brak formalnego dowodu |
| |GL(3,𝔽₂)| = 168 = (2N+1)·2^N·N! | TAUTOLOGIA | Zakłada N=3 |
| N_ν = 2.984 ± 0.008 (LEP) | EKSPERYMENT | Potwierdza 3, nie wyjaśnia |
| 4. generacja dynamicznie zakazana (g₀ > g_crit = 8/5) | NUMERYCZNE | H8: PASS — ale nie odpowiada "dlaczego" |
| Anomaly cancellation: N_gen = N_color | POWIĄZANIE | Nie dowodzi |

## Ścieżki ataku

### Ścieżka 1: Topologiczna
- Czy π₁(przestrzeni konfiguracji solitonów w d=3) wymusza 3 klasy homotopii?
- Solitony w R³ z warunkami brzegowymi g→1 → π₂(S²) = Z
- Ale dlaczego Z → Z₃?

### Ścieżka 2: Dynamiczna (stabilność)
- Udowodnić: g₀^(N+1) > g_crit = 8/5 dla N ≥ 4
- g₀^(4) = φ³ · g₀^e ≈ 4.23 > 1.6 ✓ (numerycznie)
- Ale: Dlaczego g₀^(3) = φ² · g₀^e ≈ 2.27 < 1.6? (NIE < 1.6! To > 1.6!)
- Problem: g₀^τ = 2.27 > g_crit = 1.6, a τ istnieje. Sprawdzić subtelności.

### Ścieżka 3: Informacyjna
- Entropia Shanona na substracie Z₂ w d=3
- Maksymalna informacja mutualnej między sektorami
- Czy 3 sektory maksymalizują jakiś funkcjonał informacyjny?

### Ścieżka 4: Algebraiczna (eliminacja)
- GL(2,𝔽₂) = S₃ (6 elementów): zbyt mało — nie daje pełnej struktury CKM
- GL(4,𝔽₂) (20160 elementów): daje 4 generacje — eksperymentalnie zakazane
- GL(3,𝔽₂) jest "Goldilocksem" — ale dlaczego akurat 𝔽₂?

### Ścieżka 5: Z solitonowego WKB (najkonkretniejsza)
- ODE: g'' + (2/r)g' = V'(g) z V odpowiednim dla α=2
- WKB: liczba stanów związanych N_bound = ∫√(2|V|) dr / π
- Sprawdzić: czy N_bound = 3 **dokładnie** dla α=2, d=3, K(Φ)=Φ²

## Kryterium zamknięcia

Twierdzenie: "W teorii solitonów z K=Φ², d=3, istnieją dokładnie 3 stabilne
sektory masowe" (lub równoważne).

## Pliki do scalenia z rdzeniem

- Nowy paragraf w `sek04_stale.tex` lub `sek08_formalizm.tex`
- Nowy `dodatek_N3_proof.tex`

## Referencje rdzenia

- `dodatekK_wkb_atail.tex` (analiza WKB)
- `dodatekF_hierarchia_mas.tex` (hierarchia mas)
- `nbody/examples/ex240_*.py` (N_gen z K + oscylacji)

## Status

- [ ] Ścieżka 5: WKB — policzyć N_bound analitycznie
- [ ] Ścieżka 2: Sprawdzić g₀^(N) vs g_crit formalnie
- [ ] Ścieżka 4: Algebraiczne wykluczenie GL(N≠3, 𝔽₂)
- [ ] Ścieżka 1: Topologia przestrzeni konfiguracji
- [ ] Ścieżka 3: Entropia — eksploracyjna
