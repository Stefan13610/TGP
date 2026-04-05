# Sesja v38 — Plan i analiza (2026-03-29)

## Stan wejściowy po v37

### Zamknięte w v37

| ID | Wynik |
|----|-------|
| O-L4 | Mechanizm α*≠2 wyjaśniony: gałąź fizyczna F(α), α_min≈2.54 |
| O-L5 | α*₃=2.929524 (σ=0, R_MAX=120→300, Δ=2.3×10⁻⁸) |
| ex92 | Profil F(α)+G₃(α) na [1.0,5.0], 116/120 punktów |
| ex93 | Precyzyjne α*₃, G₃ monotoniczne na [2.70,3.30] |
| ex94 | Weryfikacja R_MAX=300, analiza algebraiczna S₃ |

### Otwarte priorytety v38

| Priorytet | ID | Zadanie |
|-----------|-----|---------|
| 🔴 Wysoki | O-L6 | Lepsza precyzja G₃ → definitywne potwierdzenie/obalenie S₃=8 |
| 🟡 Średni | O-L4b | Analityczna formuła α_min: dF/dα=0 → numeryczna weryfikacja |
| 🟢 Niski | O-L7 | Algebraiczne badanie α*τμ=2.8534 (G₃/F=R32) |

---

## Analiza wstępna O-L6 (pre-session)

### Problem precyzji G₃

Z ex92–ex94: G₃(α*₃) = 3495.2 vs R31=3477.2 → diff = +18 jednostek (+0.52%).

**Źródło szumu**: A_tail(φ²·z₀) przy φ²·z₀≈3.3 ekstrapolowana z okien [20,36]→[50,66].
Krzywa A(r) maleje jako ~A_∞·(1 + a/r + b/r²) — ta forma jest dopasowywana z 4 okien.
Przy dużym g₀ (φ²·z₀≈3.3) soliton jest szybkozmienný na małych r, ale wyraźnie oscylacyjny dla r>10.

**Hipoteza**: okna [20,36]→[50,66] są za małe dla dużego g₀. Rozszerzenie do [60,80]→[100,120] może dać lepszą ekstrapolację.

### Plan ex95

```python
# Rozszerzone okna dla A_tail przy dużych g₀
FIT_WINS_EXTENDED = [(20,36),(30,46),(40,56),(50,66),(60,76),(70,86),(80,96),(90,106)]
R_MAX = 120  # wystarczy, bo okna sięgają do r=106

# Test: porównaj A_tail(φ²·z₀) ze standardowymi i rozszerzonymi oknami
# Dla α∈[2.90, 2.95], krok 0.002
```

### Analiza wstępna O-L4b

Z ex92, dane dla gałęzi fizycznej:

| α | F(α) |
|---|------|
| 2.345 | 218.6 |
| 2.432 | ~206.8 (α*₁) |
| 2.479 | 201.0 |
| 2.613 | 204.8 |
| 2.636 | ~206.8 (α*₂) |
| 2.748 | 217.7 |

Minimum F ≈ (201.0+204.8)/2 = 202.9 przy α≈(2.479+2.613)/2 = 2.546.

**Warunek analityczny**: dF/dα = 0 przy α=α_min.

F(α) = (A_tail(φ·z₀(α), α) / A_tail(z₀(α), α))^4

dF/dα = 0 ↔ d/dα [ln A_tail(φ·z₀(α), α) − ln A_tail(z₀(α), α)] = 0

To wymaga znajomości:
- ∂A_tail/∂g₀ (jak amplituda zmienia się ze startowym g₀)
- ∂A_tail/∂α (jak amplituda zmienia się z α)
- dz₀/dα (jak punkt kwantowania zmienia się z α)

Numeryczne: liczymy dF/dα z ex93 danych i szukamy analitycznego wyrażenia.

---

## Plan skryptów v38

### ex95_G3_precision.py

```python
# Cel: lepsza precyzja G₃ przez rozszerzone okna fit
# Porównanie standardowych FIT_WINS vs rozszerzonych
# Skan α∈[2.90, 2.95], krok=0.002
# R_MAX=150 (trochę wyższy niż 120 dla stabilności w r>100)
```

### ex96_F_derivative.py (opcjonalne po ex95)

```python
# Cel: numeryczna dF/dα na gęstej siatce
# Znajdowanie α_min z wysoką precyzją
# Test kandydatów algebraicznych: czy α_min = (α*₁+α*₂)/2 = 2.539?
```

---

## Priorytety v38

```
1. ex95_G3_precision — rozszerzone okna, weryfikacja S₃=8
2. Ocena wyniku ex95 — definitywne potwierdzenie/obalenie S₃=8
3. ex96_F_derivative — (jeśli czas) numeryczna analityka O-L4b
```

---

## Wyniki ex95 (in-session, 2026-03-29)

### Parametry ex95

```
R_MAX=150, α*₃=2.929524 (ex93/ex94)
FAZA 1: diagnostyka okien przy α*₃
FAZA 2: skan G₃(α) na [2.90,2.95], krok=0.001, Pool(6)
FAZA 3: brentq z oknami std/ext/far
Czas: 893s
```

---

### ❌ KRYTYCZNA KOREKTA: α*₃≈2.929 jest artefaktem numerycznym

**Profil A_tau per-okno przy α*₃, φ²·z₀=3.307:**

| Okno [rL,rR] | A_tau | Uwaga |
|--------------|-------|-------|
| [20,36] | **2.838** | **ZA NISKO** — soliton nie w asymptotyce |
| [30,46] | 3.012 | skok! |
| [40,56] | 3.010 | plateau |
| [50,66] | 2.981 | |
| [60,76] | 2.985 | |
| [80,96] | 3.006 | |
| [120,136] | **2.994** | **wartość prawdziwa** |

**Asymptotyczna A_tau ≈ 2.994** (okna [60,136] dają stabilne ~2.993–3.006).

**Porównanie Atau i G₃ z różnych zestawów okien przy α*₃=2.9295:**

| Okna | Atau | G₃ | G₃−R31 |
|------|------|-----|--------|
| std [20,66] | **2.652** | **3410** | **−67** ← ZANIŻONE |
| ext [20,106] | 2.910 | 4953 | +1476 |
| far [60,136] | 2.895 | 4859 | +1382 |

**Wyjaśnienie artefaktu**: okno [20,36] jest anomalnie niskie (2.838) bo soliton z g₀=φ²·z₀≈3.307 ma dużą strefę bliską (r<30) gdzie A(r) jeszcze nie jest asymptotyczna. Standardowe okna [20,66] włączają ten niski punkt → ekstrapolacja biasuje Atau w dół → G3_std jest fałszywie niska.

**Rzeczywista G₃ przy α∈[2.90, 2.95]:**
- G₃_ext ≈ 4700–5000 — **powyżej R31=3477 wszędzie**
- G₃_far ≈ 4724–5008 — **powyżej R31=3477 wszędzie**
- Brak zera G₃=R31 w α∈[2.90,2.95] przy poprawnych oknach

**Wniosek**: α*₃≈2.929 i S₃≈8 z ex93/ex94 są **artefaktami** metody ekstrapolacji ze standardowymi oknami. **S₃=8 NIE jest potwierdzone.**

---

### Co mówi ex95 o prawdziwym sektorze tau?

Prawdziwa G₃ (A_tau≈2.994) daje G₃≈5550 w okolicach α≈2.93 — znacznie powyżej R31=3477. Prawdziwe zero G₃=R31 musi być w innym zakresie α.

Z ex92 (std, coarse): G₃ spada poniżej R31 dopiero przy α≈3.13 lub α≈3.82 — ale te punkty też były ze std oknami (niepoprawne). **Potrzebny ex96**: skan G₃ z poprawnymi oknami (tylko [60,136]) na szerokim zakresie α∈[2.5, 5.0] żeby znaleźć prawdziwe zero G₃=R31.

---

### Status O-L5 po ex95

**KOREKTA**: poprzednie wyniki ex93/ex94 były oparte na artefakcie.

| Wielkość | ex93/ex94 (std, błędne) | ex95 (poprawna metoda) |
|----------|------------------------|------------------------|
| α*₃ | 2.929524 | **NIEOKREŚLONA** |
| G₃(α*₃) | 3495 | ~5550 (powyżej R31) |
| S₃=8 | 383 ppm (artefakt) | **NIEOKREŚLONA** |

**O-L5: powrót do stanu "Program"** — prawdziwe α*₃ wymaga ex96 z poprawnymi oknami.

---

### Plan v38 po ex95

```
✅ ex96: skan G₃(α) z oknami FAR [60,136] na α∈[2.5, 5.0] — WYKONANY
⏳ ex97: test zbieżności A_tau z R_MAX=300, okna [100,280]
⏳ O-L4b: numeryczna dF/dα → α_min analityczny
```

---

## Wyniki ex96 (in-session, 2026-03-28)

### Parametry ex96

```
R_MAX=150, okna WYŁĄCZNIE FAR [(60,76)→(120,136)]
FAZA 1: coarse skan α∈[2.5, 5.0], krok=0.05, Pool(6)
FAZA 2: brentq G₃=R31 (jeśli zmiana znaku)
FAZA 3: dense skan wokół minimów |G₃−R31|, krok=0.005
FAZA 4: weryfikacja F=R21 (α*₁,₂) z FAR oknami
Czas: 597.4 s
```

---

### ❌❌ WYNIK FUNDAMENTALNY: G₃_far > R31 wszędzie na α∈[2.5, 5.0]

**Coarse skan (krok=0.05):**

| α | F_far | G3_far | G3−R31 |
|---|-------|--------|--------|
| 2.50 | 198.13 | **4667.4** | **+1190** |
| 2.75 | 206.48 | 5026.7 | +1549 |
| 3.00 | 233.18 | 5818.8 | +2342 |
| 3.50 | 304.69 | 7934.5 | +4457 |
| 4.00 | 384.09 | 9907.4 | +6430 |
| 5.00 | 671.81 | 13177.4 | +9700 |

**Min G₃_far = 4667.4** przy α=2.50 — nadal 1190 jednostek (34%) powyżej R31=3477.

**Brak zmiany znaku G₃−R31 na całym α∈[2.5, 5.0].**

Dense skan α∈[2.4, 2.6] potwierdza: G₃_far monotonicznie rośnie od 4535 (α=2.40) do 4847 (α=2.60).

---

### Zaskakujące: F_far ma inne zera niż F_std

**α*₁ z FAR**: brentq w [2.43,2.44] → **α*₁_far ≈ 2.434** (blisko α*₁_std=2.432)
**α*₂ z FAR**: brentq w [2.75,2.80] → **α*₂_far = 2.753825** (α*₂_std=2.636 — różnica +0.118!)

| Okna | α*₁ | α*₂ | S₂=α*₁+α*₂ |
|------|-----|-----|------------|
| STD (ex93) | 2.43184 | 2.63558 | 5.06742 |
| FAR (ex96) | ~2.434 | 2.75382 | ~5.188 |
| Różnica | ~+0.002 | **+0.118** | **+0.121** |

**Interpretacja**: Near-field bias w oknie [20,36] zaniża Amu (φ·z₀≈2.04) bardziej na gałęzi prawej F(α) (α>2.54, F rośnie) niż na lewej. Skutek: α*₂_std jest systematycznie zaниżone o ~0.12.

---

### Diagnoza sytuacji po ex96

**Opcja A — α*₃ nie istnieje** (G₃_far > R31 wszędzie):
- Formuła G₃=(A_tail(φ²·z₀)/A_tail(z₀))^4=R31 nie ma rozwiązania na fizycznej gałęzi
- τ-sektor wymaga innej parametryzacji (inny mnożnik zamiast φ², inna formuła)

**Opcja B — A_tail(φ²·z₀) jest nadal zawyżone** (plateau r∈[60,136] ≠ A_∞):
- Per-okno: A_tau≈2.994 na [60,136] — stabilne (~0.22%)
- Ale czy to prawdziwa asymptota? Soliton z g₀≈3.31 potrzebuje bardzo długiej `r` do pełnej asymptotyki
- Przy R_MAX=150, „plateau" [60,136] może być nadal postresonansowym tranzytem
- Potrzeba R_MAX=300+ z oknami [150,280]: jeśli A_tau_far(r→∞)≪2.994, to G₃_true < 4667

**Opcja C — near-field bias dotyczy też Ae i/lub Amu** (ale w przeciwnym kierunku):
- F_far daje α*₂≈2.754 zamiast 2.636 — zmiana znaku A_mu przy STD vs FAR
- Systematyczne zniekształcenie całej relacji F, G₃

---

### Plan ex97

```python
# ex97_amplitude_convergence.py
# Cel: test zbieżności A_tau(g₀=φ²·z₀) z R_MAX
# Dla α=2.929 (ex93 α*₃): oblicz per-okno A_tau dla okien [20,36] aż do [250,266]
# R_MAX = 300
# Oczekiwanie: jeśli A_tau plateau ~2.994 się utrzymuje → Opcja A (α*₃ nie istnieje)
#              jeśli A_tau powoli spada poniżej 1.8 → G₃ może być < R31
```

**Wniosek wstępny**: α*₃ = NIEOKREŚLONE; S₃=8 = NIEOKREŚLONE. Potrzebna diagnostyka konwergencji A_tau do R_MAX=300.

---

## Wyniki ex97 (in-session, 2026-03-28)

### Parametry ex97

```
R_MAX=300, okna od [20,36] do [270,286] (26 okien, szerokość 16)
Trzy testowe α: 2.929524, 2.500, 3.500
Czas: 40.7 s (bez Poola — tylko 3 punkty)
```

---

### ✅ DEFINITYWNE POTWIERDZENIE: Plateau A_tau jest prawdziwą asymptotą

**α=2.929524, φ²·z₀=3.307:**

| Zakres r | A_tau | G₃ | G₃−R31 |
|----------|-------|----|--------|
| std [20,66] | 2.757 (fit) | 3989 | +512 |
| mid [60,136] | **2.994** (plateau) | 5546 | +2069 |
| far [120,220] | **2.996** (plateau) | 5558 | +2080 |
| very far [200,286] | **2.996** (plateau) | 5564 | +2087 |

**Trend A_tau na [150,286]: dA/dr = +0.000009/r → ZERO (plateau potwierdzone)**

Plateau ≈ 2.994–2.998 utrzymuje się **nieprzerwanie od r=60 do r=286**. Nie ma powolnego zaniku ku mniejszym wartościom.

**Analogicznie:**
- α=2.500: plateau A_tau ≈ 2.727 (okna [60,286]), G₃ ≈ 4671 (+1194 nad R31)
- α=3.500: plateau A_tau ≈ 3.312 (okna [60,286]), G₃ ≈ 7982 (+4505 nad R31)

---

### ❌ DEFINITYWNY WYNIK: G₃_true > R31 WSZĘDZIE

| α | A_tau (FAR, trwałe) | G₃_true | R31 | Nadwyżka |
|---|---------------------|---------|-----|----------|
| 2.500 | 2.727 | 4671 | 3477 | **+34%** |
| 2.929 | 2.996 | 5558 | 3477 | **+60%** |
| 3.500 | 3.312 | 7982 | 3477 | **+130%** |

**Formuła G₃ = (A_tail(φ²·z₀) / A_tail(z₀))^4 = R31 NIE MA rozwiązania na fizycznej gałęzi α∈[2.5, 5.0].**

---

### Interpretacja fizyczna

**Dlaczego STD dawało fałszywe zero?**

Soliton τ (g₀=φ²·z₀≈3.31) ma dużą "strefę przejściową" (near-field) dla r<50. Okno [20,36] mierzy A=2.838 — soliton jest jeszcze w tranzycie. STD ekstrapolacja: A_∞_std ≈ 2.652 (drastycznie zaniżone). G₃_std ≈ 3410 ≈ R31 — przypadkowe trafienie.

Ae (z₀≈1.26, małe g₀) jest odporne na near-field bias → Ae_std ≈ Ae_far ≈ 0.347. Cały błąd pochodzi z Atau.

**Wnioski dla TGP:**
1. Model τ-solitonu z G₃=(A_tau(φ²·z₀)/A_e(z₀))^4=R31 i mnożnikiem φ² **nie ma rozwiązania**
2. τ-lepton wymaga innej parametryzacji (inny mnożnik, inna formuła amplitudowa)
3. α*₂ z FAR okien = 2.754 (nie 2.636) — zmiana S₂ o +0.12
4. S₃=α*₁+α*₂+α*₃=8 jest hipotezą bez podstawy numerycznej

---

### Status otwarty po ex97

```
✅ ex98: analiza mnożnika m* i α*₁,₂_far — WYKONANE
⏳ O-L4b: α_min z FAR oknami (może się różnić od STD)
⏳ Nowe O-L8: Skąd pochodzi sektor τ w TGP?
```

---

## Wyniki ex98 (in-session, 2026-03-29)

### Parametry ex98

```
R_MAX=150, FAR okna, brentq F_far=R₂₁ i G₃=R₃₁
FAZA 1: α*₁_far, α*₂_far precyzyjnie
FAZA 2: profil A(g₀) vs g₀/z₀ dla 4 wartości α, krok=0.1
FAZA 3: m*(α) — brentq G₃(m·z₀)=R₃₁ dla 13 α ∈[2.43,3.50]
FAZA 4: analiza algebraiczna m*(α)
FAZA 5: skan G₃_far/F_far = R₃₂
Czas: 853 s
```

---

### FAZA 1: Precyzyjne α*₁,₂ z oknami FAR

| Wielkość | FAR (ex98) | STD (ex93) | Δ |
|----------|-----------|-----------|---|
| α*₁ | **2.436011** | 2.431838 | +0.004 |
| α*₂ | **2.753825** | 2.635577 | **+0.118** |
| S₂=α*₁+α*₂ | **5.189836** | 5.067415 | **+0.122** |
| G₃(α*₁_far) | 4588 | — | >>R₃₁=3477 |
| G₃(α*₂_far) | 5034 | — | >>R₃₁=3477 |

α*₂ przesuwa się o +0.118 — duży efekt near-field biasu dla gałęzi prawej F(α).

---

### FAZA 3: m*(α) — mnożnik dla G₃=R₃₁

| α | m* | m*/φ² | G₃(m*·z₀) |
|---|---|---|---|
| 2.43 | 2.4956 | 0.9533 | 3477.2 |
| 2.50 | 2.4895 | 0.9509 | 3477.2 |
| 2.75 | 2.4440 | 0.9334 | 3477.2 |
| 3.00 | 2.3964 | 0.9153 | 3477.2 |
| 3.50 | 2.2787 | 0.8704 | 3477.2 |

**m* zależy od α** (cv=2.51%, zakres [2.279, 2.496]) — **brak uniwersalnego mnożnika**.

**Fit liniowy**: m\*(α) = −0.199·α + 2.989

Przy **α=α_TGP=2**: m\*(2) = 2.989 − 0.398 = **2.591 ≈ φ²=2.618** (−1%)

→ Hipoteza: φ²-mnożnik jest motywowany przez α=α_TGP=2, ale fizyczne α*≠2.

---

### FAZA 5: G₃_far/F_far = R₃₂

| Zakres α | G₃/F | vs R₃₂=16.82 |
|----------|-------|--------------|
| 2.40–3.50 | 20–25 | +20%–50% |
| 4.90–5.00 | 19.4–19.6 | +15%–17% |

**Brak zmiany znaku G₃/F−R₃₂ na [2.4,5.0]** — minimum G₃/F=19.43 przy α=4.95, nadal +15% nad R₃₂.

---

### Synteza: Co wiemy o sektorze τ po ex95–ex98?

**Tabela standardowych warunków z FAR oknami:**

| Warunek | Status FAR | Uwaga |
|---------|-----------|-------|
| F=(A(φ·z₀)/A(z₀))^4=R₂₁ | ✅ MA rozwiązanie | α*₁=2.436, α*₂=2.754 |
| G₃=(A(φ²·z₀)/A(z₀))^4=R₃₁ | ❌ BRAK rozwiązania | G₃_far>R₃₁ wszędzie |
| G₃/F=(A(φ²·z₀)/A(φ·z₀))^4=R₃₂ | ❌ BRAK rozwiązania | G₃/F>R₃₂ wszędzie |
| G₃(m·z₀)=R₃₁, m=const | ❌ m nie jest stały | m*(α)∈[2.28,2.50], zależy od α |

**Wniosek**: Hierarchia φ (e→μ→τ) z mnożnikami (1, φ, φ²) nie działa dla sektora τ z poprawnymi amplitudami. Model TGP w obecnej postaci nie przewiduje masy τ przez G₃=R₃₁.

---

### O-L8 (nowe): Skąd pochodzi sektor τ?

Kandydaci do dalszego badania:
1. **α_TGP=2 jako punkt kanoniczny**: m*(2)≈φ², czyli warunek τ jest spełniony przy α=2 z przeskalowanym mnożnikiem
2. **Inna formuła amplitudowa**: może τ = (A(z₀)/A_ref)^4=R₃₁ gdzie A_ref to amplituda przy innym warunku
3. **F=R₃₁**: czy istnieje α* takie że F_far=(A(φ·z₀)/A(z₀))^4=R₃₁? (wymaga dużych α, F rośnie monot.)
4. **Energia solitonu**: może masa leptonu = energia solitonu, nie amplituda
5. **Model wymaga rozszerzenia**: τ może wymagać nowego wkładu (np. SU(3) lub innego sektora)

---

## Wyniki ex99 (in-session, 2026-03-29)

### Parametry ex99

```
R_MAX=150, FAR okna, skan F_far(α) na α∈[5,15], krok=0.5
FAZA 1: coarse F, G₃ na [5,15]
FAZA 2: brentq F_far=R₃₁
FAZA 3: dense skan wokół α*₃_F
FAZA 4: analiza algebraiczna
Czas: 2191 s
```

---

### ⚠️ Dzikie oscylacje F_far i G₃_far dla α>5

| α | F_far | G₃_far | Uwaga |
|---|-------|--------|-------|
| 5.0 | 672 | 13177 | gładkie |
| 7.0 | 1928 | **781** | G₃ PONIŻEJ R₃₁! |
| 7.5 | 2518 | **13640** | G₃ POWYŻEJ R₃₁ ← |
| 8.0 | 3306 | **422927** | |
| 8.5 | 5332 | **3311950** | |
| 13.0 | **58** | 653006 | F prawie zero! |
| 13.5 | **38409** | 1671002 | skok o czynnik 660! |
| 14.5 | **10** | 55281 | F prawie zero! |

**Oscylacje F i G₃ przy α>8 są rezonansami fazowymi** — okna FAR [60,136] trafiają na zero/maksimum oscylacji A(r)·r, co daje skoki o rzędy wielkości. Nie są fizycznymi przejściami.

---

### F_far = R₃₁ przy α≈8.074

brentq znalazło kilka "skrzyżowań", ale:
- **α*₃_F ≈ 8.074**: F_far=3467.6 (vs R₃₁=3477.2, diff −9.6, −0.28%)
  - W gładkim regionie [7.5, 8.5], G₃_far=422927 >> R₃₁
  - α*₃_F/α*₁_far = 8.074/2.436 = 3.314 ≈ √11?
  - S₃_F = 2.436+2.754+8.074 = **13.264** — brak czystej formuły
  - α*₃_F ≈ **8** = 4·α_TGP: bliskość do 8 intryująca, różnica 0.074 (~9200 ppm)
- Pozostałe "skrzyżowania" (α≈12.38, 13.28, 14.42): F_far≈0 przy tych α (rezonans fazowy) — ARTEFAKTY

### G₃_far = R₃₁ przy α≈7.1–7.4 (nowe!)

Z coarse scan (niezaplanowane odkrycie):
- α=7.0: G₃_far=781 (PONIŻEJ R₃₁=3477) → phi² multiplier DZIAŁA przy dużym α!
- α=7.5: G₃_far=13640 (POWYŻEJ R₃₁)
- Przejście G₃_far=R₃₁ między α=7.0 i α=7.5
- Ale: F_far≈1928–2518 >> R₂₁=206.77 → warunek elektronu-muona niezachowany
- Czy to fizyczne czy rezonans? Wymaga testu R_MAX=300 (stabilność)

---

### Synteza po ex95–ex99

**Tabela skrzyżowań z FAR oknami:**

| Warunek | α gdzie spełniony | F lub G₃ tam | Fizyczny? |
|---------|------------------|--------------|-----------|
| F=R₂₁ (muon) | α*₁=2.436, α*₂=2.754 | ✅ | ✅ tak |
| G₃=R₃₁ (tau, φ²) | **brak na [2.5,5.0]** | G₃>>R₃₁ | ❌ |
| G₃=R₃₁ (tau, φ²) | **α≈7.1–7.4** | F≠R₂₁ | ❓ (rezonans?) |
| F=R₃₁ | **α≈8.074** | G₃>>R₃₁ | ❓ (rezonans?) |
| G₃/F=R₃₂ | **brak na [2.4,5.0]** | zawsze >16.82 | ❌ |

**Wniosek**: W fizycznym zakresie α∈[2.4,5.0] — gdzie spełniony jest warunek muona F=R₂₁ — warunek taona G₃=R₃₁ NIE jest spełniony. Sektor τ w modelu TGP wymaga albo innego zakresu α, albo nowej fizyki.

---

### Status O-L5 po v38

**O-L5 = OTWARTY FUNDAMENTALNIE**

| Wynik | Status |
|-------|--------|
| α*₃≈2.929 (ex93/ex94) | ❌ ARTEFAKT STD okien |
| G₃_far=R₃₁ na [2.5,5.0] | ❌ BRAK rozwiązania |
| m*(α) = const | ❌ m* zależy od α (cv=2.5%) |
| G₃_far=R₃₁ przy α≈7.1 | ❓ potencjalne, ale F≠R₂₁ |
| F_far=R₃₁ przy α≈8.074 | ❓ bliskie α=8=4·α_TGP, ale niefizyczne |

**Priorytet v39**: ex100 — weryfikacja G₃_far=R₃₁ przy α∈[7.0,7.5] z R_MAX=300 (test stabilności vs rezonans)
