---
title: "FAZA 6 (FINAL): X = e²/4 vs α-em — separate sectors, hipoteza odrzucona"
date: 2026-05-01
type: phase-results
phase: 6
status: CLOSED-NEGATIVE — α-em bridge odrzucone po sprawdzeniu TGP-equations
parent: "[[tgp_emergent_dirac_propagator.md]]"
related:
  - "[[r3_phase6_alpha_em_connection.py]]"
  - "[[PHASE2_n_alpha_derivation.md]]"
  - "[[../../core/formalizm/dodatekO_u1_formalizacja.tex]]"
tags:
  - TGP
  - R3
  - emergent-dirac
  - alpha-fine-structure
  - phase6-final
  - separate-sectors
  - amplitude-vs-phase
---

# FAZA 6 (FINAL) — X = e²/4 vs α-em: separate sectors

> **Status:** CLOSED-NEGATIVE. Po sprawdzeniu rzeczywistych równań TGP
> w `core/formalizm/dodatekO_u1_formalizacja.tex`, hipoteza "X jako
> TGP-modified α-em" została **odrzucona w sposób honest**:
> X należy do **amplitude sector** (R3, grawitacja), α_em do **phase sector**
> (U(1), EM). To są **dwa oddzielne sprzężenia** w TGP (`J_amp` vs `J_phase`).

---

## 1. Korekta poprzedniej Fazy 6

Wcześniejsza wersja (`PHASE6_Q5_R5_bridge_first_attempt.md`) próbowała
wyprowadzić X = e²/4 przez **5-wymiarową strukturę Kaluza-Klein** —
co było **zewnętrzną teorią importowaną**, nie naturalnym rozszerzeniem TGP.

User słusznie zauważył:
- "R⁵" w R3 było tylko **wewnętrznym aliasem** dla mass formula sector
  w `r3_atail_bridge.py`
- Standardowe 5D KK ma poważne problemy z wymaganiami wstępnymi
  (kompaktyfikacja, higher modes, hierarchy, Lorentz invariance)
- TGP **explicit nie ma** extra dim w aksjomatach (`TGP_FOUNDATIONS.md`
  §1: "jedno fundamentalne pole skalarne Φ z Z₂")

Stara hipoteza R⁵-bridge **wycofana** jako conceptual import bez
ugruntowania w TGP-formalizmie.

---

## 2. Nowa hipoteza (insight uzytkownika 2026-05-01)

> "X = e²/4 przypomina trochę stałą struktury subtelnej. Skoro w TGP
>  c, ℏ, ε₀ nie są stałe to trochę tłumaczyłoby odchyły. Chodzi o
>  zapis stałej struktury subtelnej jeżeli podstawimy c=1 itd."

### 2.1 Strukturalne podobieństwo

Stała struktury subtelnej w **Heaviside-Lorentz natural units**
(ℏ = c = ε₀ = 1):

```
α_HL = e_charge² / (4π)
```

R3 mass formula coefficient:

```
X_R3 = e_Euler² / 4 = 1.84726
```

**Formalna różnica: tylko czynnik 1/π** (4π → 4).

### 2.2 Semantyczna pułapka

Litera "e" oznacza **dwie różne rzeczy**:
- W α: `e_charge` — ładunek elementarny (bezwymiarowy w naturalnych
  jednostkach: e_charge ≈ 0.30282 w HL natural)
- W X: `e_Euler` — liczba Eulera ≈ 2.71828

To są **niezwiązane standardowo**. Pytanie czy w TGP istnieje **unifying
calibration** gdzie się łączą.

---

## 3. Numeryczne testy — co odkryto

### 3.1 X i α to **różne wielkości numerycznie**

```
X_R3   = 1.84726
α_SI   = 0.007297 (1/137.036)
Ratio  = X / α_SI = 253.14
```

Różnica **4 rzędów wielkości** wyklucza prostą interpretację
"X to α w innych jednostkach".

### 3.2 Test K = X/α na clean integers

| Kandydat | Wartość | diff |
|----------|---------|------|
| 253 | 253.000 | -0.056% |
| 256 (2^8) | 256.000 | +1.13% |
| 254 | 254.000 | +0.34% |
| 252.5 | 252.500 | -0.25% |

253 jest **najbliżej**, ale 253 = 11·23 (prime factorization), bez
oczywistego fizycznego znaczenia. **Brak czystej clean matchu.**

### 3.3 Niefizyczna kalibracja gdy zmusimy X = α

Jeśli założymy X = α_HL_TGP w TGP-natural units, dostajemy:

```
e_charge_TGP = √(4π·X) = √23.21 = 4.818
```

vs standard `e_charge_HL ≈ 0.303` — **15.9× większy**. To **niefizyczne**
dla ładunku elektronu.

---

## 4. Trzy hipotezy (niedomknięte)

### H1: Numerologia
X = e²/4 ma podobną strukturę algebraiczną do α_HL, ale to
**przypadkowa zbieżność formy**. Nie ma fundamental connection.

### H2: Deep analogy (najbardziej obiecująca)
Obie formuły pochodzą z **wave-function renormalization** lub
**partition function evaluation**. W TGP, z varying c/ℏ/ε₀,
"fundamental coupling unit" mógłby mieć inną strukturę niż α_SI.

**Możliwy mechanizm:** w TGP gdzie `c_loc(ψ) = c·√A(ψ)` (Faza 1+3),
i `ℏ_loc`, `ε₀_loc` mogą varies, **lokalna α_TGP(ψ)** może mieć formę:

```
α_TGP(ψ) = e²/(4π·ε₀_loc·ℏ_loc·c_loc)
        = α_SI · [ε₀·ℏ·c] / [ε₀_loc·ℏ_loc·c_loc(ψ)]
```

Jeśli ratio `[ε₀·ℏ·c] / [ε₀_loc·ℏ_loc·c_loc]` osiąga taki konkretny
wzór że α_TGP ↔ X_R3, mielibyśmy connection. **Ale to wymagałoby
ekstremalnej modyfikacji** (factor ~250×) — **niefizyczne** dla
generic TGP.

### H3: TGP-α emerguje jako different coupling
W TGP może istnieć **inne fundamental coupling** α_TGP = e²/4 gdzie
"e" ma inne fizyczne znaczenie (np. exp(action), topological winding,
RG fixed point). To NIE jest standardowe α_em — to **TGP-specific
fundamental constant** o STRUKTURALNIE podobnej formie do α_em.

W tym wypadku R3 mass formula używa α_TGP jako **kanonicznego coupling
unit**, a α_em w naszym świecie jest **innym (znacznie mniejszym)
parameter**.

---

## 5. TGP varying c/ℏ/ε₀ — sprawdzenie

User słusznie zauważył: w TGP `c, ℏ, ε₀` nie są stałe (z M9.1''
metryki, c_loc varies z ψ). Sprawdzono dla generacji:

| gen | ψ | A(ψ) | c_loc/c | 1/A |
|-----|------|--------|---------|-----|
| vac | 1.000 | 1.000 | 1.000 | 1.000 |
| e | 0.950 | 1.211 | 1.100 | 0.826 |
| μ | 1.155 | 0.463 | 0.681 | 2.159 |
| τ | 1.288 | 0.106 | 0.325 | 9.469 |

**Lokalna prędkość światła dla τ jest 0.325·c** — dramatycznie różna
od vacuum! Jeśli ε₀_loc i ℏ_loc compensują c_loc (przez Maxwell
relation `ε₀·μ₀·c² = 1`), α_TGP(ψ) **może być stałe** lokalnie.

Ale to **nie tłumaczy** dlaczego X = e²/4 zamiast e²/(4π) — czynnik π
musi pochodzić z innego mechanizmu.

---

## 6. Co byłoby potrzebne żeby zamknąć

### 6.1 Najczystsza ścieżka: **dedykowany cykl** R3 ↔ α-em

TGP ma już cykl `op-alpha-fine-structure` (status DERIVED refined²).
Powinien istnieć **bridge cycle** który łączy:

- R3 mass formula (z X = e²/4)
- α_em derivation w op-alpha-fine-structure
- Wspólne fundamental constants (jeśli istnieją)

Taki cykl mógłby:
1. Spróbować wyprowadzić **obie** wartości (X i α_em) z **jednego**
   strukturalnego anchora
2. Sprawdzić czy "e_Euler" i "e_charge" mają wspólne pochodzenie w TGP
3. Pokazać czy 4π → 4 ma natural origin w TGP-Z₂ symmetry lub
   radial measure

### 6.2 Konkretne pytania do osobnego cyklu

1. Czy `op-alpha-fine-structure` przewiduje `α_em = e_Euler²/(4π·K)`
   gdzie K to TGP-strukturalny factor?

2. Czy istnieje **bezwymiarowa stała w TGP** o wartości ~253 (= X/α_SI)?

3. Czy TGP-modified Coulomb law w substracie (ψ ≠ 1) zmienia formalną
   strukturę α?

4. Czy "e_Euler" w R3 ma fundamental meaning lub tylko statistical fit?

---

## 7. HONEST framing dla user'a publication agent

**To co publication agent powinien wiedzieć:**

### Co jest pewne (verified)
- R3 mass formula `m_obs = c_M · A² · g₀^[e²(1-α/4)]` matche PDG
  μ/e na 0.014%
- Numeryczny X = 1.847 wygrywa wśród alternatywnych prostych form
  (e²/4, 37/20, (3+e·φ)/4) ale przewaga ~20%
- Spin-1/2 emerguje strukturalnie z RP² topology
- m_0 = 0 (zero bare mass) potwierdzone

### Co jest sugestywne ale niedowiedzione
- X = e²/4 ma podobną strukturę algebraiczną do α_HL
- Może istnieć TGP-fundamental coupling α_TGP strukturalnie analogiczny
  do α_em ale z różną kalibracją
- To może łączyć R3 z innym sektorem TGP (op-alpha-fine-structure)

### Co jest **niewspierane** w obecnym stanie
- "X = e²/4 jest fundamentalna stała derywowana z TGP"
- "X to α_em w TGP-natural units"
- "5D Kaluza-Klein R⁵-bridge derywuje X" (wycofane jako external import)

### Honest framing
W publication: **"X = e²/4 jest empirycznym odkryciem wymagającym
dedykowanego cyklu derywacji. Strukturalne podobieństwo do α_HL =
e²/(4π) jest sugestywne ale wymaga wspólnego anchora w TGP-formalizmie."**

---

## 8. Open problems (po Fazie 6 redo)

1. **Czy istnieje TGP-cykl łączący R3 z op-alpha-fine-structure?**
   (Najczystsza ścieżka derywacji X)

2. **Czy 4π → 4 ma natural origin w TGP-Z₂ symmetry?**
   (Może radial measure dla soliton w 3D z Z₂-projection daje 4 zamiast 4π)

3. **Co to znaczy "e" w R3 fundamentally?**
   (Liczba Eulera, ładunek, oba, lub trzecia interpretacja?)

4. **Czy α_em jest constant w TGP?**
   (Ze varying c_loc(ψ), local α_TGP może varies — testowalne przez
   atomic clocks w polach grawitacyjnych)

---

## 9. Pliki Fazy 6 (REDO)

| Plik | Zawartość |
|------|-----------|
| `r3_phase6_alpha_em_connection.py` | Numeryczna analiza X vs α_em |
| `r3_phase6_alpha_em_connection.txt` | Output (negative explicit derivation) |
| `PHASE6_alpha_em_connection.md` | Ten dokument |
| `r3_phase6_r5_bridge.py/.txt` | Wcześniejsza próba R⁵-bridge (DEPRECATED) |
| `PHASE6_Q5_R5_bridge_first_attempt.md` | Wcześniejszy dokument (DEPRECATED) |

---

## 10. Quick check w TGP-equations (przed inwestowaniem w pełen cykl)

Użytkownik 2026-05-01: *"przed miesiącami pracy najpierw sprawdziłbym czy
w równaniach na c lub ε₀ nie pojawia się π w mianowniku, i policzył to
na równaniach TGP. Może okaże się, że w ramach TGP przyjęcie 1 dla tych
zmiennych jest po prostu niemożliwe."*

To było **kluczowe pragmatyczne pytanie** — 15-minutowy check zaoszczędził
miesiące potencjalnie bezowocnej pracy.

### 10.1 Co znalazłem w `core/formalizm/dodatekO_u1_formalizacja.tex`

TGP **MA** explicit derywację α_em z substratu (sekcja O.4, linie 300-403):

```
α_em^(sub) = e²·μ₀·c₀ / (4π·ℏ₀)
           = Φ₀ / (8π·J_phase)
           ≈ 24.66 / (8π · 92.5)
           ≈ 1/94.09  (Planck scale)
```

**Czynnik 8π = 2·(4π) explicit w mianowniku** — czyli `4π` rzeczywiście
jest naturalne w TGP-EM derivation.

### 10.2 Co potwierdziło Twoją intuicję

**Pozytywne ustalenia:**
- ✓ TGP rzeczywiście ma `4π` (właściwie `8π`) w α_em derivation
- ✓ `c, ℏ, ε₀, μ₀` są emergent — **NIE stałe fundamentalne**
  - `μ₀ = ℏ₀² / (2J·v²·a_sub²·e²)` — funkcja ładunku!
  - `ε₀ = 1/(c₀²·μ₀)` — emergent
- ✓ Substytucja `c = ε₀ = 1` strukturalnie modyfikuje teorię
  (te wielkości są inter-tied przez `Φ₀, J, a_sub`)

### 10.3 Co odrzuciło hipotezę α-bridge

**Negatywne ustalenia (decydujące):**

- ✗ **Numerycznie X ≠ α**:
  - X_R3 = 1.847
  - α_em^(sub) ≈ 0.0106 (Planck) lub 0.0073 (SI)
  - Różnica **2 rzędy wielkości** — nie da się skleić bez
    ekstremalnej redefinicji ładunku (15× standard)

- ✗ **Czynnik strukturalny: 8π vs 4 to nie kwestia "varying constants"** —
  to **różne mianowniki**:
  - α_em^TGP ma 8π (z full sphere + Heaviside)
  - X_R3 ma 4 (z 3D radial measure + Hobart-Derrick)

- ✗ **Liczba Eulera (e_Euler) NIE pojawia się** w TGP-α-derivation.
  TGP używa: `Φ₀, J, π, 2π` (z winding number `θ ∈ [0, 2π)`).
  Brak natural source dla `e_Euler ≈ 2.718` w U(1) sector.

### 10.4 Strukturalne wyjaśnienie: separate sectors w TGP

`dodatekO_u1_formalizacja.tex:405-421` (uwaga "Dwa parametry sprzężenia"):

> "TGP ma jeden substrat Γ = (V, E), ale dwa sektory sprzężeń:
> - **J_amp**: sprzężenie amplitudowe (|ψ_i|) — daje pole Φ, geometrię, grawitację.
> - **J_phase**: sprzężenie fazowe (θ_i) — daje pole A_μ, elektromagnetyzm."

R3 mass formula używa **J_amp sector** (amplituda → R3 soliton, mass).
α_em pochodzi z **J_phase sector** (faza → photon, EM coupling).

To są **dwa oddzielne sprzężenia w tym samym substracie**. Mogą być różne
nawet jeśli mikroskopicznie zaczęły się jako `J_amp = J_phase = J` (RG flow
może je rozdzielić).

**Brak natural reason dla X (amplitude) ↔ α_em (phase) bridge w TGP.**

---

## 11. Wnioski

1. **User insight był strukturalnie ostry** — `e²/4` vs `e²/(4π)` to
   genuine formalne podobieństwo. **TGP rzeczywiście ma 4π w α_em** —
   to nie było wymyślone.

2. **Pragmatyczny check zaoszczędził miesiące** — dedykowany cykl
   R3↔α-em okazałby się bezowocny, bo sektory są **strukturalnie
   oddzielne** w TGP-formalizmie.

3. **Numerycznie X ≠ α** — różnica 2 rzędów wielkości + brak common
   `e_Euler` source w U(1) sector wyklucza fundamental connection.

4. **Honest classification dla publication agent**:

   > X = e²/4 to **EMPIRICAL FIT w R3 amplitude sector** z `e_Euler`
   > **statystycznym anchor** (lepszym niż 37/20, (3+e·φ)/4 o ~20%).
   > Strukturalne podobieństwo do α_HL = e²/(4π) jest **formalne, nie
   > fizyczne**: oba sektory (amplitude vs phase) są oddzielne w TGP.
   > Derywacja X pozostaje OPEN, ale ścieżka **w obrębie R3 amplitude
   > sector** (RG flow R3 ODE samo, bez U(1)), nie przez α-em bridge.

5. **Co NIE jest podtrzymywane** w obecnym stanie TGP:
   - "X = e²/4 jest fundamentalna stała derywowana z TGP"
   - "X to TGP-modified α-em"
   - "R⁵-bridge derywuje X" (wycofane jako external import)

---

## 12. Otwarte ścieżki **w obrębie R3 amplitude sector**

Bez α-em bridge, nadal można szukać derywacji X = e²/4 w:

1. **R3 ODE jako effective theory** — 1-loop renormalization w samym
   R3 ODE (nie 5D, nie U(1) — czysto wewnątrz amplitude sector)

2. **Hobart-Derrick balance** — α=4 jest punktem gdzie n=0; może
   `e²` pojawia się przez specyficzny scaling argument przy α=4

3. **Wave-function renormalization Z_φ** dla R3 amplitude field —
   może AS NGFP `g* = 0.71, λ* = 0.19, η_N* = -2` (UV.1) daje
   anomalous dimension `e²/4·(4-α)` w odpowiednim limicie

4. **Statystyczna interpretacja** — X = 1.847 jako "mean field
   coupling" R3 amplitude; może `e²/4` jest po prostu **statistical
   fit**, a "true" formula ma inny strukturalny charakter

Żadna z tych ścieżek nie wymaga α-em bridge. Wszystkie pozostają OPEN.

---

## 13. Pliki Fazy 6 (FINAL)

| Plik | Status |
|------|--------|
| `r3_phase6_alpha_em_connection.py` | Numeryczna analiza X vs α-em (negative) |
| `r3_phase6_alpha_em_connection.txt` | Output skryptu |
| `PHASE6_alpha_em_connection.md` | Ten dokument (CLOSED-NEGATIVE) |
| `r3_phase6_r5_bridge.py/.txt` | DEPRECATED (R⁵ external import) |
| `PHASE6_Q5_R5_bridge_first_attempt.md` | DEPRECATED |

---

## 14. Status końcowy R3 (po wszystkich Fazach 1-6)

```
┌────────────────────────────────────────────────────────────────────┐
│  R3 + Emergent Dirac Program — STATUS FINAL 2026-05-01            │
│                                                                    │
│  ✅ Faza 1: ψ↔g₀ identification (CLOSED)                          │
│  ✅ Faza 2: Mass formula m = c·A²·g₀^[e²(1-α/4)] (CLOSED)         │
│  ✅ Faza 3: RP² → spin-1/2 (CLOSED)                               │
│  ✅ Faza 4: Yukawa coupling, m_0=0 (CLOSED)                       │
│  ✅ Faza 5: Full propagator, vacuum=Dirac (CLOSED)                │
│  ⚠️  Faza 6: α-em bridge (CLOSED-NEGATIVE)                        │
│         X = e²/4 nie wynika z α-em (separate sectors)             │
│         Pozostaje EMPIRICAL fit w R3 amplitude sector             │
│                                                                    │
│  Mass ratios:                                                      │
│    m_μ/m_e: -0.001% PDG  ✓                                        │
│    m_τ/m_e: -0.085% PDG  ✓ (R3 native formula)                    │
│                                                                    │
│  Open problems (poza scope):                                       │
│    - Analytical derivation X = e²/4 (w R3 amplitude sector)       │
│    - Tau drift -11.57% w Faza 5 propagator (reconciliation        │
│      A^(5-α) vs A²·g₀^(e²/2))                                     │
│    - Pełna 4D dynamika, gauge structure                           │
└────────────────────────────────────────────────────────────────────┘
```

---

**Autor:** Faza 6 FINAL — po pragmatycznym check'u uzytkownika 2026-05-01.
**Data:** 2026-05-01.
**Status:** CLOSED-NEGATIVE. α-em bridge odrzucone.
**Wnioski:** R3 amplitude sector i U(1) phase sector są oddzielne w TGP;
X = e²/4 pozostaje empirical fit w R3, nie α-related.
