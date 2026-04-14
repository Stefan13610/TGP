# R5: Prawo skalowania m ∝ A_tail⁴

## Problem

Formuła masowa m_n = c_M · A_tail⁴ jest **fundamentem sektora leptonowego**.
Daje r₂₁ = (A_μ/A_e)⁴ = 206.768 z dokładnością 0.0001% wobec PDG (206.77).

Ale potęga k=4 jest **postulatem**. Dlaczego nie k=3 lub k=5?

## Obecne uzasadnienie

| Argument | Typ | Status |
|----------|-----|--------|
| Wymiarowy: w d=3, zbieżność ogona wymaga n > 2d/(d-1)=3, więc k=4 | HEURYSTYKA | ex188: OK |
| E²=0 (wiriał): E_bind ∝ A² | NIEKOMPLETNE | Działa, ale nie daje k=4 |
| E³=0: E_perturbacyjne | OBALONY | E³ = -0.647 ≠ 0 |
| prop:K-exponent (stopień V → wykładnik) | OBALONY | ex150: DISPROVED |

## Plan ataku

### Ścieżka 1: Twierdzenie wiriałowe (najprostsza)
- Soliton z K(Φ)=Φ^α w d=3
- Derrick's theorem: (d-2)·E_kin = d·E_pot (dla skalowania r → λr)
- Z α=2: E_kin = ∫Φ²(∇g)²d³x, E_pot = ∫P(g)d³x
- Sprawdzić: czy E_bind ∝ A_tail^{2α} = A⁴ wynika z α=2

### Ścieżka 2: Z funkcjonału działania
- S[g] = ∫[K(g)·(∇g)² + P(g)]d³x
- Masa jako ∂S/∂(parametr deformacji)
- Parametr: g₀ (amplituda centralna)
- m(g₀) = dS/dg₀ — jakie potęgi A_tail pojawiają się?

### Ścieżka 3: Analiza asymptotyczna ogona
- Ogon solitonu: g(r) ~ 1 - A·exp(-μr)/r^ν dla r → ∞
- A = A_tail = amplituda ogona
- Masa jako moment zerowy: m ∝ ∫(g-1)d³x
- Obliczenie: m = ∫(A·e^{-μr}/r^ν)·4πr²dr = 4πA·∫r^{2-ν}·e^{-μr}dr
- Dla ν=1 (Yukawa): m ∝ A/μ² — liniowe w A!
- Ale: masa **fizyczna** = energia wiązania, nie moment zerowy
- Trzeba obliczyć E_bind[A] dokładnie

### Ścieżka 4: Relacja między A_tail a g₀
- Z ODE: A_tail = f(g₀) — jaką ma postać?
- Jeśli A ∝ (g₀ - 1)^{1/4} → m ∝ A⁴ ∝ (g₀-1)
- Sprawdzić numerycznie skalowanie A(g₀)

## Kryterium zamknięcia

**Twierdzenie:** "Dla solitonu z K(Φ)=Φ^α, α=2, d=3: masa ∝ A_tail^{2α} = A⁴"

## Pliki do scalenia z rdzeniem

- Rozszerzenie `dodatekJ_ogon_masy.tex`
- Nowy akapit w `dodatekF_hierarchia_mas.tex`
- Skrypt dowodowy → `scripts/`

## Referencje rdzenia

- `dodatekJ_ogon_masy.tex` (teoria ogona)
- `dodatekF_hierarchia_mas.tex` (hierarchia)
- `dodatekK_wkb_atail.tex` (WKB + A_tail)
- `scripts/ex188_A4_dimensional_argument.py`
- `scripts/ex150_*.py` (obalenie prop:K-exponent)
- `scripts/ex152_*.py` (energie solitonów)

## Status

- [ ] Ścieżka 1: Wiriał z α=2 — obliczenie E_bind(A)
- [ ] Ścieżka 4: Numeryczne skalowanie A(g₀) → sprawdzenie A ∝ (g₀-1)^{1/4}
- [ ] Ścieżka 2: Masa z dS/dg₀
- [ ] Ścieżka 3: Asymptotyka ogona
- [ ] Połączenie ścieżek w dowód
