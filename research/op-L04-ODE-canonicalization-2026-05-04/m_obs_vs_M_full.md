---
title: "L04 — analiza fizyczna: m_obs vs M_full distinction"
date: 2026-05-04
parent: "[[README.md]]"
type: physical-analysis
tgp_owner: research/op-L04-ODE-canonicalization-2026-05-04
tags:
  - L04
  - m_obs
  - M_full
  - mass-distinction
  - physical-analysis
related:
  - "[[../why_n3/CORRECTIONS_2026-05-01.md]]"
  - "[[../why_n3/PHASE2_n_alpha_derivation.md]]"
---

# Analiza fizyczna: `m_obs` vs `M_full` distinction

## 1. Sedno problemu

Pierwotny audit L04 ([[../../audyt/L04_ODE_dualism_alpha]]) wskazywał
na **pozorny dualizm**: dwie formulacje ODE (α=1 K=g² i α=2 K=g⁴)
dające inne `g₀^e` i inne wzory masy `M ∝ A^k`.

**Insight użytkownika 2026-05-01** rozwiązał ten dualizm:

> „m = c·A_tail⁴ to masa OBSERWOWALNA, nie pełna masa cząstki. Do
> bariery potrzebujemy pełną masę, nie obserwowalną."

Distinction `m_obs ≠ M_full` zmienia całą interpretację. Dualizm α=1
vs α=2 NIE jest dualizmem teorii — jest *artefaktem mieszania dwóch
różnych mas*.

## 2. Dwa rodzaje masy w TGP

### 2.1 `M_full` — pełna energia wewnętrzna (mass-internal)

```
M_full(g₀, α) = K + V_eff
            = ∫_0^∞ [½·K(g)·(g')² + V(g)] · 4πr² dr
```

gdzie:
- `K(g) = g^(2α)` — kinetic prefactor
- `V(g) = (1/3)g³ − (1/4)g⁴` (z β=γ=1 vacuum condition)
- ODE: `g'' + (2/r)g' + α(g')²/g = V'(g)/K(g)`

**Cechy strukturalne:**

- Definiowane jako **całka po całym profilu solitonu** (not asymptotic)
- Zachowane w sensie energii dla regularnych rozwiązań ODE
- **Bariera g₀_crit operuje na M_full** — to *strukturalna własność
  ODE*, niezależna od α: dla `g₀ > g₀_crit(α)` rozwiązanie ODE wchodzi
  w singularność `g_min → 0` w skończonym `r_max`
- Może być **ujemne** dla solitonów typu „excess" w false vacuum

**Z r3_observable_vs_full_mass.txt (Sekcja 1):**

| α | M_full(e) | M_full(μ) | M_full(τ) |
|---|-----------|-----------|-----------|
| 0.5 | -0.0034 | -0.0356 | -0.1302 |
| 1.0 | -0.0032 | -0.0507 | -0.2214 |
| 2.0 | -0.0027 | -0.0588 | -0.1771 |

`M_full < 0` dla wszystkich α w fizycznym zakresie (e/μ/τ). Wartości
są O(0.01-0.3) i nie są w prostej proporcji do PDG mas (m_e:m_μ:m_τ =
1:206:3477).

### 2.2 `m_obs` — masa obserwowalna (mass-asymptotic)

```
m_obs(g₀, α) = c_M · A_tail²(g₀, α) · g₀^[e²(1−α/4)]
```

gdzie:
- `A_tail` = amplituda asymptotycznego rozwiązania `δg(r) ~ A_tail·exp(-m·r)/r`
- `c_M` = uniwersalna stała kalibracyjna (pojedynczy parametr dla całej
  drabiny generacji)

**Cechy fizyczne:**

- Definiowane przez **asymptotyczną projekcję** ogona oscylacyjnego
  (analog ADM mass w GR)
- **Empirycznie zweryfikowane** dla α∈[0.25, 4.0] z residuum <0.1%
  w fizycznym zakresie [α=1, α=2]
- **Spektakularna zgodność** z PDG dla α=2 (TGP-canonical):
  - `m_μ/m_e` = 206.77 (PDG 206.7682), **diff −0.001%** ✓
  - `m_τ/m_e` = 3477.4 (z Koide K=2/3, PDG 3477.23), **diff −0.085%** ✓
- **n(α) = e²·(1−α/4)** — zaskakujące pojawienie się klasycznej stałej e²

**Wartości n(α):**

| α | n(α) = e²(1−α/4) | Specjalna wartość |
|---|-------------------|---------------------|
| 0 | e² ≈ 7.389 | (formal limit) |
| 1 | 3e²/4 ≈ 5.542 | substratowa α=1 (R3 oryginalne) |
| 2 | e²/2 ≈ 3.695 | **TGP-canonical (K=g⁴)** |
| 3 | e²/4 ≈ 1.847 | (X = e²/4 unit) |
| 4 | 0 | **Hobart-Derrick balance** (m=c·A² czyste) |

## 3. Dlaczego `m_obs ≠ M_full`

### 3.1 Analogie w innych teoriach

W innych teoriach pola istnieje analogiczna distinction:

| Teoria | Mass-internal | Mass-asymptotic |
|--------|--------------|------------------|
| **GR** | masa Komara `M_K = (1/4πG)∫_S T^t_t dS` | masa ADM `M_A = (1/16πG)∫(∂_iH_ii − ∂_jH_ij)dS_j` |
| **QFT** | bare mass `m_0` (UV scale) | renormalized mass `m_phys` (IR scale, on-shell) |
| **EM** | energia pola `E_field = ½∫E²dV` | całkowity ładunek `Q = ∮E·dS` |
| **Solitony NLσM** | full energy `E_NL = ∫T^00 d³x` | topological charge `Q_top` |
| **TGP** | **`M_full` (strukturalne)** | **`m_obs` (asymptotyczne)** |

### 3.2 Mechanizm różnicy w TGP

`m_obs` jest projekcją asymptotyczną — *wagą z dystansu*. Soliton wokół
otoczki `Φ_eq` daje fluktuacje `δΦ ~ A_tail·exp(-m_sp·r)/r`. Obserwator
w nieskończoności mierzy `m_obs` przez prawo Newtona/grawitacji
asymptotycznej. **Mass-tail coupling** zależy od kinetic prefactor:

```
m_obs ~ A² · (kinetic dressing) · (core renormalization)
```

Wave-function renormalization Z(g₀, α) = g₀^n(α):

- Dla α=4 (Hobart-Derrick): K(g)=g⁸ — kinetic „payouje" maksymalnie,
  brak core dressing → `m_obs ~ A²` (klasyczne)
- Dla α=0: K(g)=1 stałe — kinetic nie zna geometrii, *maximal* core
  dressing → `m_obs ~ A²·g₀^e²` (dziedziczy całą stałą e²)
- Pomiędzy: liniowa interpolacja `n(α) = e²(1−α/4)`

### 3.3 Bariera `g₀_crit` — własność strukturalna ODE

Mechanizm bariery (R3 PHASE 1, [[../why_n3/r3_phase1_psi_g0_identification.txt]]):

- Dla `g₀ > g₀_crit(α, d)` rozwiązanie ODE solitonu blow-uje (g_min → 0)
- Bariera istnieje **dla każdego α** — to topologiczna własność klasy
  ODE, nie konkretnego α
- **Wartość g₀_crit zależy od α:**
  - α=1: g₀_crit = 2.2062
  - α=2: g₀_crit = 1.8744 (= 4/3·1.40554, Lorentzian horizon M9.1''
    via PHASE1 identification)

**Niezależność od mass formula:**

> Bariera operuje strukturalnie na profilu solitonu (M_full), niezależnie
> od mass formula A^p details. Mass formula jest osobnym problemem —
> observable mass = c·A^p(α).
>
> [[../why_n3/r3_observable_vs_full_mass.txt]] Sekcja 4

## 4. Co tłumaczy „dualizm" α=1 vs α=2

### 4.1 LP-4 (k=4) — działa tylko dla α=1

[[../../meta/PLAN_DOMKNIECIA_MASTER.md]] LP-4 zamknął się 9/9 PASS na
**k=4** z trzema argumentami:

1. **Zero-mode**: E_kin = E_pot virial → E₂ = 0
2. **Konwergencja wymiarowa**: `k = 2(d-1)/(d-2) = 4` w d=3
3. **Dyskryminacja**: k=3 → r₂₁=55, k=4 → 207, k=5 → 784

**Krytyczne zauważenie:** LP-4 9/9 PASS używał **substratowego ODE
K=g²** (`scripts/lp4_mass_exponent_verification.py`):

> Substrate r₂₁ = 206.74 (δ = 0.013%) z k_eff = 4.000 dokładnie...
> Formulacja substratowa K=g² jest numerycznie stabilna i daje dokładne
> wyniki. Kanoniczna K=g⁴ jest niestabilna dla g₀ > 1.3.
>
> [[../../meta/PLAN_DOMKNIECIA_MASTER.md]] LP-4 finding

W świetle Phase 2 analysis: `k=4` to **specjalny przypadek α=1** w
formule `m_obs = c·A²·g₀^[e²(1−α/4)]`:

```
α = 1: m_obs = c · A² · g₀^(3e²/4) ≈ c · A² · g₀^5.542
```

Empiryczne `(A_μ/A_e)^4 = 207` dla α=1 wynika z faktu, że dla `g₀^μ = φ·g₀^e`:
- `(g₀_μ/g₀_e)^5.542 = φ^5.542 ≈ 14.30`
- `(A_μ/A_e)^2 ≈ 14.37` (numerycznie)
- Ich iloczyn ≈ 206

To znaczy: **pseudo-uniwersalność** `m ∝ A^4` jest *artefaktem α=1*.
Dla α=2 (TGP-canonical) `(A_μ/A_e)^3.69` musi być pomnożone przez
`(g₀_μ/g₀_e)^3.69 ≈ φ^3.69` ≈ 5.91², żeby dać 207. R5 K² (rozumiane
jako m=c·A^4) **fail spektakularnie dla α=2** (+490% drift).

### 4.2 Argument konwergencyjny `k = 2(d-1)/(d-2) = 4`

LP-4 argument konwergencyjny **nie jest universal**. Jest sformulowany
dla:

- ODE: `g'' + (2/r)g' + α(g')²/g = V'(g)/K(g)`
- Asymptotyka tail: `δg ~ A_tail · exp(-m·r)/r` (Yukawa form)
- Skończoność energii: `E = ∫(½K(g)(g')² + V(g))r²dr < ∞`

**Dla α=1 (K=g²)**: warunek skończoności energii daje k = 4. Ale
**dla α=2 (K=g⁴)** ten sam warunek daje **inną liczbę** k (ze skanu
empirycznego k_eff dla α=2 jest ≈3, nie 4).

Pełny argument konwergencyjny jako *wielowymiarowa* relacja
**k(α, d)** wymaga osobnej analizy (otwarty problem,
[[../../audyt/L05_mass_exponent_drift]]).

## 5. Dlaczego Phase 2 (z e²) jest fundamentalna

### 5.1 R5 ↔ Phase 2 bridge — analytical theorem

[[../mass_scaling_k4/R5_PHASE2_ANALYTICAL_BRIDGE_2026-05-02.md]] **Twierdzenie**:

> R5 mass formula `m = c·K²` (z K~A² universal) jest **równoważna**
> Phase 2 universal `m_obs = c_M·A²·g₀^[e²(1−α/4)]` **wtedy i tylko
> wtedy gdy α = 1**.

**Dowód** (sekcja 2):

Empirical scaling Phase 2: `m_obs ~ A^(5−α)`, czyli `c_M A² g₀^n(α) ~ A^(5−α)`.

Slope condition Phase 2: `slope = log(g₀)/log(A) = (3−α)/n(α)`

Slope condition R5: `c A^4 = c_M A² g₀^n(α)` ⇒ `slope_R5_req = 2/n(α)`

Equivalence: `(3−α)/n(α) = 2/n(α)` ⇔ `3−α = 2` ⇔ **α=1** ∎

### 5.2 Co to znaczy

R5 K² mass formula (i implicit LP-4 k=4) **NIE jest fundamentalnym
mechanizmem niezależnym od Phase 2** — jest *strukturalną konsekwencją*
Phase 2 universal mass formula dla specyficznego substratu α=1.

Phase 2 jest **fundamental**, R5 K² to **derivative**.

Dla TGP-canonical α=2:
- R5 K² ratio = 1221 (PDG mass ratio = 207, **mismatch +490%**)
- Phase 2: m_μ/m_e = 206.77, **diff −0.001%**

To jednoznacznie wskazuje, że **α=2 + Phase 2 universal** jest
poprawną kanoniczną formulacją TGP.

## 6. Podsumowanie

| Wielkość | α=1 (substratowa K=g²) | α=2 (TGP-canonical K=g⁴) |
|----------|-------------------------|----------------------------|
| Status w sek08 | NIE wybrana strukturalnie | **WYBRANA** przez `thm:D-uniqueness` |
| Phase 2 mass formula | Działa (n(1) = 3e²/4) | Działa (n(2) = e²/2) |
| LP-4 k=4 | k=4 ✓ specjalny przypadek | k=3 (nie k=4) |
| R5 K² formula | Spójna z Phase 2 | **Niespójna** (R5 ≡ Phase 2 IFF α=1) |
| Mass spectrum vs PDG | −0.013% (substrate ODE) | **−0.001%** (Phase 2 + Koide) |
| Bariera g₀_crit | 2.2062 | 1.8744 (= 4/3·1.40554, Lorentzian horizon) |
| g₀^τ pod barierą | 1.7293 (margin +0.477) | 1.755 (margin +0.119) |
| 4-ta gen zakazana | g₀^4=2.798 > 2.21 ✓ | g₀^4=2.840 > 1.874 ✓ |

**Werdykt:** TGP jest **jednoznacznie kanoniczna z α=2**. „Dualizm"
α=1 vs α=2 jest pozorny — α=1 to *specjalny przypadek* universal
Phase 2 mass formula, która jest fundamentalna dla każdego α.

## 7. Co pozostaje otwarte

- **X = e²/4 RG derivation** — Phase 6 Q5 R⁵-bridge NEGATIVE; e² to
  *empirical discovery*, awaiting analitycznej derywacji z RG flow lub
  Hobart-Derrick balance
- **k(α, d) generalization** — pełna analiza warunku skończoności
  energii dla różnych α w d=3
- **Dlaczego α=4 jest Hobart-Derrick balance** — n(4) = 0 numerycznie
  ścisłe; topologiczne lub wariacyjne uzasadnienie pożądane

Patrz [[NEEDS.md]] dla pełnej listy.

## Cross-references

- [[README.md]] (cykl L04)
- [[canonical_form_evidence.md]] (3 niezależne dowody α=2)
- [[ODE_class_taxonomy.md]] (klasy operatorów C1-C3)
- [[mass_formula_unification.md]] (Phase 2 universal vs LP-4/LP-6)
- [[../why_n3/PHASE2_n_alpha_derivation.md]]
- [[../why_n3/CORRECTIONS_2026-05-01.md]]
- [[../mass_scaling_k4/R5_PHASE2_ANALYTICAL_BRIDGE_2026-05-02.md]]
