# R6: B=√2 analitycznie z ODE solitonowego

## Problem

Stosunek Brannena B = b/a parametryzuje masy leptonów:
```
√mᵢ = a·(1 + b·cos(θ + 2πi/3)),    i = 1,2,3
```

Z danych: B_num = 1.414212... ≈ √2 (|δ| < 10⁻⁶)

Konsekwencja: K = (1 + B²/2)/N = (1 + 1)/3 = **2/3** (Koide)

Gdyby B=√2 było **udowodnione analitycznie** → Koide staje się **twierdzeniem TGP**.

## Obecny status (2026-04-16)

### ✅ UDOWODNIONE (algebraicznie)

| Element | Status | Dowód |
|---------|--------|-------|
| B = √(N-1) ⟺ K = 2/N | **TWIERDZENIE** | Tożsamość: K = (1+B²/2)/N, więc K=2/N → B²=2 |
| K = 2/3 ⟺ fazy 120° | **TWIERDZENIE** | Equidistant cos(2πi/3) → Σcos=0, Σcos²=3/2 |
| K = 2/3 z PDG mas | **WERYFIKACJA** | Q_K(PDG) = 1.500014 ≈ 3/2 |
| Fazy 120° TRYWIALNE | **TWIERDZENIE** | Fourier na Z₃: DOWOLNE 3 liczby → fazy 120° (DFT) |
| B = √2 ↔ K = 2/3 ↔ CV(√m) = 1 | **TWIERDZENIE** | Łańcuch tożsamości algebraicznych |

### ✅ NUMERYCZNE (2026-04-15)

| Element | Wynik | Dokładność |
|---------|-------|------------|
| r₂₁ = (A_μ/A_e)⁴ | 206.55 | 0.10% od PDG |
| g₀^τ(Koide) → r₃₁ | 3474.15 | 0.09% od PDG |
| B(Koide) | 1.41421356 | |K-2/3| = 6×10⁻¹⁰ |
| eta asymmetry c₁ | 0.72538 | stałe do 10⁻⁵ |

### ⚠️ OTWARTE

| Element | Status | Problem |
|---------|--------|---------|
| Dlaczego CV(√m) = 1? | **OTWARTE** | = dlaczego K=2/3? = dlaczego g₀^τ = 1.729? |
| Co wyznacza g₀^τ? | **KLUCZOWE** | g₀^τ ≠ φ²·g₀^e (bariery), ≠ 2·g₀^e (1% off) |
| g₀^τ/g₀^μ ≈ √(3/2) | **DO ZBADANIA** | 3/2 = K⁻¹ — przypadek? |

## Kluczowe wyniki R6

### Łańcuch dowodowy

```
Level 0: Fazy 120° TRYWIALNE (Fourier na Z₃)         ✅ UDOWODNIONE
Level 1: B = √2 ↔ K = 2/3 ↔ CV(√m) = 1             ✅ TOŻSAMOŚCI
Level 2: K(PDG) = 0.666660 ≈ 2/3                     ✅ EMPIRYCZNE
Level 3: DLACZEGO K = 2/3 z dynamiki solitonu?        ⚠️ OTWARTE
```

### Struktura η(δ) — korekcja nieliniowa (2026-04-15)

Definicja: `η(δ) = A_tail(1+δ)/|δ|`  (δ = g₀ - 1, ze znakiem)

```
Kluczowe właściwości:
  η → 1     dla |δ| → 0   (teoria perturbacji: A ~ |δ|)
  η_def < 1  (deficit: g₀ < 1, solitony deficytowe)
  η_exc > 1  (excess: g₀ > 1, solitony nadmiarowe)

Asymetria:
  (η_exc(δ) - η_def(δ))/δ → c₁ = 0.72538  (STAŁE!)
  Jedyna korekcja nieliniowa: z termu (1/g)g'² w ODE
  
Fizyczne wartości:
  η_e  = 0.954 (deficit, δ_e = -0.131)
  η_μ  = 1.161 (excess, δ_μ = +0.407)
  η_μ/η_e = 1.217 → boost faktor (η_μ/η_e)⁴ = 2.195
  
Konsekwencja:
  r₂₁ = (δ_μ/δ_e)⁴ × (η_μ/η_e)⁴ = 94.1 × 2.195 = 206.6
  Liniowe przybliżenie daje 94 → η korekta daje prawidłowe 207!
```

### Tau constraint (2026-04-15)

```
PROBLEM: φ²·g₀^e = 2.276 > g₀_crit = 2.250 (ZABLOKOWANE!)
         Phi-drabinka dla tau nie działa.

Koide daje: g₀^τ = 1.729
  g₀^τ/g₀^e = 1.989 ≈ 2      (diff 0.55%)
  g₀^τ/g₀^μ = 1.229 ≈ √(3/2)  (diff 0.37%)

Hipoteza g₀^τ = 2·g₀^e: K = 0.673 (1% off od 2/3)
Hipoteza g₀^τ = √(3/2)·g₀^μ: g₀^τ = 1.723 (0.35% off)

WNIOSEK: g₀^τ prawie = 2·g₀^e LUB √(3/2)·g₀^μ,
ale żadna algebraiczna forma nie daje DOKŁADNIE K = 2/3.
```

### Pełny łańcuch derywacji

```
WHAT WE CAN DERIVE:
  1. g_ij = g·δ_ij → substrate ODE (α=1)       [PROVEN]
  2. g₀^e = 0.86941 ← Compton matching          [INPUT: λ_C]
  3. g₀^μ = φ·g₀^e ← φ-drabinka                 [ASSUMPTION]
  4. r₂₁ = 206.55   ← ODE + (2)+(3)             [DERIVED: 0.10%]
  5. g₀^τ = 1.729   ← ???                        [OPEN GAP]
  6. K = 2/3         ← zależy od (5)             [≡ step 5]
  7. B = √2          ← algebraicznie z (6)       [PROVEN if K=2/3]
  8. N = 3           ← bariera + (6)             [PROVEN if K=2/3]

THE GAP: Step 5. What determines g₀^τ?
  → K = 2/3 jest WIĄZANIEM, nie derywacją (na razie)
```

### Negatywne wyniki (ważne!)

1. **F(φ) nie jest stałe** — CV = 220%, Ścieżka 4 nie prowadzi do B=√2
2. **φ²-drabinka tau nie działa** — g₀^τ = 2.28 > g₀_crit = 2.25
3. **c_M = E/A⁴ nie jest stałe** — zmienia się od -20 do -1200
4. **r₂₁ NIE jest uniwersalne** — zależy silnie od g₀^e (CV = 460%)
5. **g₀^τ = 2·g₀^e daje K = 0.673, nie 2/3** — 1% off
6. **Fazy δᵢ ogona NIE są Z₃-symetryczne** (2026-04-16, `r6_tail_phase_z3.py`)
   - δ_e ≈ 176.71°, δ_μ ≈ 5.64°, δ_τ ≈ 3.76°
   - Różnica e↔μ to 171° (nie 120°), μ↔τ to 1.88°
   - Interpretacja: δ ma **nieciągłość przy g₀=1** (deficit vs excess soliton),
     więc Koide K=2/3 **NIE pochodzi** od Z₃-działania na fazy ogona.
   - Test 3/6 PASS (Koide PASS, Z₃/Brannen FAIL)
7. **K=2/3 NIE jest ekstremum żadnego lokalnego funkcjonału S[g₀^τ]**
   (2026-04-16, `r6_koide_variational.py`)
   - Przetestowano 10 fizycznych funkcjonałów: Shannon H, E_solitonic,
     L2-mass, Fisher info, free energy, CV(√m), Sum log m, K_Koide sama.
   - **Wszystkie mają ekstrema na brzegu** skanu (g₀=1.40 lub 2.15), nie w interior.
   - g₀^τ(Koide) = 1.72932 **nie jest punktem krytycznym** żadnego z nich.
   - g₀^τ/g₀^e zmienia się od 2.00 do 3.40 przy zachowaniu K=2/3
     (CV = 22% — brak drabinki algebraicznej).
   - Test 2/4 PASS (Koide inversion PASS, stałe ratio FAIL).
   - **Wniosek**: K=2/3 **nie pochodzi z lokalnej zasady wariacyjnej**
     na parametrach g₀.
8. **c₁ ≠ 1 - ln(3)/4** (2026-04-16, `r6_c1_high_precision.py` + `r6_c1_richardson.py`)
   - Pomysł: `4·(1-c₁) = ln(3)` ⇔ `exp(4·(1-c₁)) = 3` — asymetria deficit/excess
     ODE kodowałaby N=3 generacji poprzez entropię log-3.
   - Pomiar wysoko-precyzyjny (DOP853, rtol=1e-13, δ ∈ [0.002, 0.05]):
     c₁(δ) = 0.7252580 + (-4.9×10⁻²)·δ² + O(δ⁴), residuum fit < 10⁻⁷
   - Ekstrapolowane c₁ = 0.72525802
   - `1 - ln(3)/4 = 0.72534693` — **diff = 8.9×10⁻⁵**, czyli 100× większa
     niż precyzja ekstrapolacji. **REJECTED**.
   - Wniosek: c₁ nie ma prostej postaci zamkniętej wśród testowanych
     kandydatów klasycznych (ln 3, π, e, φ w prostych kombinacjach).
   - Pozostaje: c₁ ≈ 0.72526 jest **subtelnym artefaktem dynamiki ODE** α=1,
     być może analityczny w oryginalnym języku funkcji specjalnych (Bessel?).

## Pliki

| Plik | Opis | Status |
|------|------|--------|
| `r6_fourier_z3_proof.py` | Fourier na Z₃ + Koide: 9/9 PASS | ✅ RDZEŃ |
| `r6_atail_ratio_analysis.py` | Ścieżka 4: F(φ), 4/9 PASS | ✅ NEGATYWNY |
| `r6_eta_koide_attack.py` | η(δ) korekcja + c₁ = 0.725 | ✅ NOWE |
| `r6_koide_from_ode.py` | K jako f(g₀^e): g₀^τ inversion | ✅ NOWE |
| `r6_tau_constraint.py` | Hipotezy na g₀^τ: test 2·g₀^e, √(3/2)·g₀^μ | ✅ NOWE |
| `r6_tail_phase_z3.py` | Test Z₃ na fazach δᵢ ogona ODE | ✅ NEGATYWNY |
| `r6_koide_variational.py` | Test zasady wariacyjnej K=2/3 | ✅ NEGATYWNY |
| `r6_c1_closed_form_test.py` | Skan zamkniętych form c₁ | ✅ NOWE |
| `r6_c1_high_precision.py` | Wysoko-precyzyjny pomiar c₁ (DOP853) | ✅ NOWE |
| `r6_c1_richardson.py` | Richardson ekstrapolacja c₁ → δ=0 | ✅ NEGATYWNY |

## Ścieżki dalszego ataku

### Ścieżka A: g₀^τ/g₀^μ = √(3/2)?
- 3/2 = K⁻¹ — jest związek kauzalny?
- Czy soliton ODE wymusza √(K⁻¹) jako krok tau?
- Powiązanie z R₃₁: g₀ symetria lustra

### Ścieżka B: Entropia / topologia
- Czy g₀^τ minimalizuje jakiś funkcjonał? → **WYKLUCZONE lokalne** (`r6_koide_variational.py`, 2026-04-16)
- Pozostaje: funkcjonał **nielokalny** (np. over tail-wrap number)
- Topological charge constraint → otwarte

### Ścieżka C: Z₃ ⊂ GL(3,𝔽₂)
- Formalizacja: Z₃ na masach → K = 2/3 jako ZASADA SYMETRII
- Nie potrzeba derywacji z ODE — Koide z GRUPY
- Ale dlaczego GL(3,𝔽₂) a nie inna grupa?

### Ścieżka D (NOWA, 2026-04-16): Deficit–excess asymmetry c₁ = 0.72538

Po wykluczeniu:
- (R6.6) Z₃ na fazy δᵢ ogona
- (R6.7) Lokalnej zasady wariacyjnej S[g₀] z K=2/3 jako ekstremum

**Pozostaje znacząca STAŁA numeryczna:** asymetria deficit/excess
`c₁ = (η_exc(δ) - η_def(δ))/δ → 0.72538` (stała do 10⁻⁵).

**Hipoteza D1**: c₁ jest **jedyną fizyczną skalą** TGP na której ODE różnicuje
rozgałęzienia deficit (e) vs excess (μ,τ). Jeśli K=2/3 da się wyrazić przez c₁,
mamy derywację Koide z asymetrii substrate.

**Hipoteza D2** (TESTOWANE 2026-04-16): c₁ = 1 - ln(3)/4?
   - Pomysł: `4·(1-c₁) = ln(3)` ⇔ c₁ wiąże N=3 przez Shannon entropię.
   - Pomiar wysoko-precyzyjny: c₁ = 0.72525802 (Richardson extrap., 10⁻⁷)
   - `1 - ln(3)/4 = 0.72534693` — diff 8.9×10⁻⁵ → **REJECTED** (100× nad precyzję)
   - Żaden z klasycznych kandydatów (ln, π, e, φ) nie pasuje.

**Plan ZREWIDOWANY**:
1. Sprawdzić czy c₁ = limit całki z funkcji Bessela (perturbacyjna struktura ODE)
2. Rozwiązać perturbacyjnie: g = 1 + δ·f, gdzie f spełnia liniowy ODE z Besselem
3. Jeśli c₁ = ∫... (Bessel integral) to mimo że nie jest "ładne", to ma interpretację

### Ścieżka E (NOWA): Tail winding number jako constraint

Skoro Z₃ nie działa na **fazach ogona**, może działa na **liczbie nawinięć** (winding).
Tail-wrap number n(g₀) = integer, bo ogon sin(r+δ)/r jest okresowy.
Może: n_e + n_μ + n_τ = 0 mod 3 **wymusza** K=2/3?

## Status checklist

- [x] Ścieżka 4: F(φ) = A(φg₀)/A(g₀) — NIE stałe
- [x] Łańcuch algebraiczny: Z₃ → 120° → K=2/3 → B=√2 — KOMPLETNY
- [x] Best-fit tau: B = 1.4142 z |K-2/3| = 6×10⁻¹⁰
- [x] η(δ) korekcja: c₁ = 0.72538, perturbacyjnie zrozumiana
- [x] K(g₀^e) skan: K zależy od g₀^e i g₀^τ
- [x] g₀^τ candidates: 2·g₀^e najlepszy (0.55%), ale nie dokładny
- [x] r₂₁ nie jest uniwersalne — zależy od g₀^e
- [x] Z₃ na fazy ogona δᵢ — WYKLUCZONE (negatyw, 2026-04-16)
- [x] Lokalna zasada wariacyjna dla K=2/3 — WYKLUCZONE (2026-04-16)
- [ ] **Ścieżka D**: c₁ = 0.72538 analitycznie (czy ln(3)/4?)
- [ ] **Ścieżka E**: tail winding number n(g₀) → Z₃ constraint
- [ ] Derywacja g₀^τ z ODE (the gap!)
- [ ] Związek g₀^τ/g₀^μ ≈ √(3/2) z K
