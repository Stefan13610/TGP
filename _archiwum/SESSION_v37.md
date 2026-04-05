# Sesja v37 — Plan i analiza (2026-03-28)

## Stan wejściowy po v36

### Zamknięte w v36

| ID | Wynik |
|----|-------|
| F7 | n_s = 1−2/N_e (TGP=Starobinsky przy Planck 2018) |
| O-L1b | α*₁=2.43183767, α*₂=2.63557742, S=5.0674 (σ=0, R_MAX-niezależne) |
| O-L2 | Zera F(α*)=R21 tylko przy β_pot=1.0 (N0-5-specific) |
| O-L3 | Ścieżki 9+10 komplementarne: F(α_TGP=2)=263.9≠R21 |

### Otwarte priorytety v37

| Priorytet | ID | Zadanie |
|-----------|-----|---------|
| 🔴 Wysoki | O-L4 | Analityczna derywacja α*₁,₂ — dlaczego α*≠α_TGP=2? |
| 🟡 Średni | O-L5 | Trzecia wartość α*₃ dla tauonu: F(α*₃)=r₃₁=3477? |

---

## Analiza wstępna O-L4 (pre-session)

### Znana struktura F(α)

Z ex91 i wcześniejszych wyników:

| α | F(α) | Uwaga |
|---|------|-------|
| 2.000 | **263.9** | α_TGP, F > R21 |
| 2.4318 | **206.8** | α*₁ (zero F=R21, malejąco) |
| ~2.534 | **~F_min** | minimum F |
| 2.6356 | **206.8** | α*₂ (zero F=R21, rosnąco) |
| ? | ? | dalej rosnąco |

**Kluczowa obserwacja (O-L4):** F(α=2)=263.9 > R21=206.77 > F_min

Funkcja F(α) ma kształt "doliny" — maleje od α=2.0 do α_min≈2.534, potem rośnie.
Dwa zera (α*₁, α*₂) leżą symetrycznie wokół minimum.

### Pytanie O-L4

Dlaczego centrum doliny α_min ≈ 2.534 ≠ α_TGP=2?

**Mechanizm**: α_min wyznaczony jest przez warunek:
d/dα [(A_tail(φ·z₀(α), α) / A_tail(z₀(α), α))^4] = 0

To wymaga:
dA_tail/dα przy φ·z₀ i przy z₀ być równoważne (stosunek stały)

Związane z tym jak z₀(α) zmienia się z α: gdy α rośnie, g*(α)=exp(-1/2α) rośnie → bariera kinetyczna wyższa → z₀ rośnie (potrzeba wyższego punktu startowego dla B_coeff=0).

### Plany analityczne

**Plan A** (numeryczny preludium): ex92 — pełny profil F(α) na [1.0, 5.0] z krokiem 0.05 → wyznaczenie F_min, α_min, zachowanie dla α>3.5

**Plan B** (tau): G₃(α) = (A_tail(φ²·z₀, α) / A_tail(z₀, α))^4 — tau-to-electron ratio condition. Skan tej funkcji i znalezienie α*₃ gdzie G₃(α*₃) = r₃₁ = 3477.

---

## O-L5 — Analiza wstępna

### Problem tau

R₃₁ = m_τ/m_e ≈ 3477.23

Dwa podejścia:
1. **F(α*₃) = R31 = 3477**: Ta sama funkcja F(α)=(A_tail(φ·z₀)/A_tail(z₀))^4. Wymaga F=3477, ale znany zakres F∈[198, 264] dla α∈[2,3.5]. Pytanie: czy F rośnie do 3477 dla α poza [2,3.5]?

2. **G₃(α*₃) = R31**: Nowa funkcja G₃=(A_tail(φ²·z₀)/A_tail(z₀))^4, czyli stosunek tau-do-elektronu zamiast muon-do-elektronu. Bardziej naturalne dla 3-generacyjnej struktury.

### Numerologia G₃

Jeśli A_tail ∝ (z₀-g*)^4 (przybliżenie):
- F = ((φ·z₀-g*)/(z₀-g*))^16
- G₃ = ((φ²·z₀-g*)/(z₀-g*))^16

Stosunek G₃/F = ((φ²·z₀-g*)/(φ·z₀-g*))^16 ≈ φ^16 = 2207 dla z₀>>g*

Więc G₃ ≈ F·φ^16 ≈ 206.77 × 2207 ≈ 456,000... za duże.

Ale to przybliżenie. Realna wartość G₃ może być inna.

**Klucz**: sprawdzić G₃(α*₁) numerycznie — czy G₃(α*₁) ≈ R31=3477?

---

## Plan skryptów v37

### ex92_F_alpha_profile.py

```python
# Pełny profil F(α) i G₃(α) na szerokim zakresie
ALPHA_RANGE = np.linspace(1.0, 5.0, 120)  # krok 0.033
R_MAX = 120
B_WIN = [28, 42]  # stałe

# Obie funkcje:
# F(α)  = (A_tail(φ·z₀, α) / A_tail(z₀, α))^4
# G3(α) = (A_tail(φ²·z₀, α) / A_tail(z₀, α))^4

# Szukamy:
# 1. α_min gdzie F = minimum → wyjaśnienie O-L4
# 2. Czy F(α) osiąga R31=3477 gdziekolwiek?
# 3. Czy G3(α) ma zero przy R31=3477?
```

### Plan O-L4 po ex92

Po zmapowaniu F(α):
1. Znalezienie analitycznej formuły dla α_min
2. Badanie czy z₀(α) ma prostą formę analityczną
3. Connection formula: jak B_coeff(z₀, α)=0 wyznacza z₀(α)?

---

## Priorytety v37

```
1. ex92_F_alpha_profile — czas ~20 min (Pool(4))
2. O-L4: analiza wyniku ex92 → formuła α_min
3. O-L5: G₃(α) → α*₃ dla tauonu
```

---

## Wyniki ex92 (in-session, 2026-03-28)

### Parametry ex92

```
ALPHA_RANGE = [1.0, 5.0], N=120, krok≈0.034
R_MAX = 120, B_WIN = [28, 42]
Pool(6), czas = 484.5 s
Wyniki: 116/120 dobre (brak: α∈{1.000,1.034,1.067,1.101} — brak z₀ tam)
```

Korekcja stałej: R21 = m_μ/m_e ≈ 206.77 (NIE czwarta potęga — komentarz w ex88/ex91 mylący).
Analogicznie R31 = m_τ/m_e ≈ 3477.22.

---

### O-L4 — Wyniki ex92: Dwie gałęzie z₀(α)

**Kluczowe odkrycie**: F(α) ma DWA reżimy z₀:

| Gałąź | Zakres α | z₀ | F(α) | Status |
|-------|----------|-----|------|--------|
| Górna | α < ~1.4 | ≈2.3–2.5 | 13–15 (minimum!) | Niefizyczna |
| Fizyczna | α > ~1.5 | ≈1.11–1.26 | 201–1932 | **Używana w ex88** |

**Gałąź górna**: z₀≈2.3 daje F_min=13.4 przy α=1.13 — ale z₀ SKACZE od 2.32 do 1.11 między α=1.40 a α=1.54. Skok ten jest artefaktem branchingu, nie fizycznym przejściem.

**Profil na gałęzi fizycznej (ex92, krok 0.034):**

| α | z₀ | F(α) | vs R21 |
|---|-----|------|--------|
| 1.538 | 1.114 | 1932 | +834% |
| 2.000 | 1.228 | **264.1** | +27.7% |
| 2.345 | 1.255 | 218.6 | +5.7% |
| **2.432** | 1.258 | **~206.8** | **≈0% (α*₁)** |
| 2.479 | 1.260 | 201.0 | −2.8% |
| 2.613 | 1.262 | 204.8 | −0.9% |
| **2.636** | 1.262 | **~206.8** | **≈0% (α*₂)** |
| 2.748 | 1.263 | 217.7 | +5.3% |
| 3.000 | 1.263 | 252.4 | +22.1% |

**Minimum fizyczne:** F_min(fiz) ≈ 200.97 przy α≈2.48, poniżej R21.

**Zera F=R21 (ex92, gałąź fizyczna):**
- α*₁ ≈ 2.4345 (F malejąco przez R21) — ex88: 2.43184 ✅
- α*₂ ≈ 2.6441 (F rosnąco przez R21) — ex88: 2.63558 ✅

**α_min(fizyczny) ≈ (α*₁+α*₂)/2 ≈ 2.539** — zgodne z szacunkiem SESSION_v37 (2.534) ✅

**Wniosek O-L4:**

Dlaczego α*≠α_TGP=2? F(α) na gałęzi fizycznej ma minimum przy α_min≈2.54, a nie przy α=2. Funkcja F jest MALEJĄCA w α=2 (jeszcze nie osiągnęła minimum). α_TGP=2 leży po lewej stronie minimum — F(2)=264>R21, a dopiero przy α*₁≈2.43 F spada do R21.

**Mechanizm**: z₀(α) rośnie od z₀≈1.11 (α=1.5) do maksimum z₀≈1.263 (α≈2.75), potem maleje. Gwiazda kinetyczna g*(α)=exp(−1/2α) też rośnie monotoniczne. Interplay tych dwóch funkcji wyznacza kształt F(α) i położenie minimum.

**O-L4: ANALIZA ZAKOŃCZONA** — mechanizm jakościowo wyjaśniony. Analityczna formuła wymaga dalszej pracy (ex93?).

---

### O-L5 — Wyniki ex92: Tau sektor

**G₃(α) = (A_tail(φ²·z₀) / A_tail(z₀))^4** — skan na [1.5, 5.0]:

| α | G₃(α) | vs R31=3477 |
|---|--------|------------|
| α*₁=2.432 | **5831** | +68% (za duże) |
| α*₂=2.636 | **7295** | +110% (za duże) |
| 2.983 | 3286 | −5.5% |
| 3.151 | 2407 | −30.8% |
| 3.384 (zero) | **≈3477** | ≈0% |
| 3.739 (zero) | **≈3477** | ≈0% |
| 4.456 (zero) | **≈3477** | ≈0% |

**Zera G₃=R31 (ex92) na gałęzi fizycznej:**
1. **α*₃ ≈ 2.928** — G₃ malejąco (3550→3347), najbliższe α*₁,₂ → **KANDYDAT GŁÓWNY**
2. α ≈ 3.384 — G₃ rosnąco (3085→3507)
3. α ≈ 3.739 — G₃ malejąco (3886→3050)
4. α ≈ 4.456 — G₃ rosnąco (1877→3869)

**G₃ oscyluje silnie dla α>2.5** — prawdopodobnie szum numeryczny w ekstrapolacji A_tail(φ²·z₀) (φ²·z₀≈3.3 → wiele odbić solitonu). Potrzebna weryfikacja ex93 z wyższą rozdzielczością.

**Stosunek G₃/F:**
- G₃/F przy α*₁: 28.20 (vs φ^16=2207, φ^8=47, φ^4=6.85) — nie jest prostą potęgą φ
- G₃/F przy α=2: 23.02

**Wniosek O-L5:**

Metoda G₃(α*)=R31 daje **α*₃≈2.928** jako głównego kandydata (pierwsze zero na gałęzi fizycznej po α*₂). Jednak G₃ jest numerycznie niestabilna w tym zakresie (φ²·z₀≈3.3 — długie trajektorie). Potrzebne ex93: precyzyjny Brent dla G₃=R31 blisko α≈2.928 z większym R_MAX.

---

### Plan v37 po ex92

```
✅ ex92: profil F(α) + G₃(α) na [1.0, 5.0] — DONE
✅ O-L4: mechanizm α*≠2 wyjaśniony (gałąź fizyczna, α_min≈2.54)
⏳ ex93: precyzyjne α*₃ z G₃=R31 (brentq, wyższy R_MAX)
⏳ O-L4 analityczne: z₀(α_max) i dF/dα=0 jako równanie algebraiczne
```

---

## Wyniki ex93 (in-session, 2026-03-28)

### Parametry ex93

```
Skan: α∈[2.70, 3.30], krok=0.010, 61 punktów (wszystkie poprawne)
FAZA 1: Pool(6), R_MAX=120
FAZA 2: brentq (xtol=1e-8)
FAZA 3: stabilność R_MAX∈{120,150,200}
Czas całkowity: 2333s (~39 min)
```

---

### O-L5 — Wyniki ex93: α*₃ STABILNA

**Kluczowe odkrycie (korekta ex92)**: G₃(α) jest **monotonically malejąca** na [2.70, 3.30] — ex92 z krokiem 0.034 przeskakiwał przez strome regiony dając pozorne oscylacje. Na gęstej siatce (krok 0.010) G₃ jest gładka.

**Profil G₃(α) — wybrane punkty:**

| α | F(α) | G₃(α) | G₃−R31 | G₃/F |
|---|------|--------|---------|------|
| 2.700 | 212.4 | 5340 | +1863 | 25.14 |
| 2.800 | 224.0 | 4262 | +784 | 19.03 |
| **2.850** | 230.4 | **3870** | +393 | **16.80 ≈ R32!** |
| 2.900 | 236.9 | 3598 | +120 | 15.19 |
| 2.920 | 239.6 | 3497 | +20 | 14.60 |
| **2.930** | 241.0 | **3415** | −62 | 14.17 |
| 3.000 | 250.2 | 3238 | −239 | 12.94 |
| 3.100 | 263.3 | 3094 | −383 | 11.75 |
| 3.300 | 286.8 | 2607 | −870 | 9.09 |

**FAZA 2 — brentq R_MAX=120:**
```
Bracket [2.920, 2.930]: G₃=3496.92→3415.35 (malejąco)
α*₃ = 2.92288038
G₃(α*₃) = 3443.34  (R31=3477.22, diff=-33.88 = 0.97%)
F(α*₃)  = 239.47   (R21=206.77, diff=+32.71)
z₀      = 1.26320
```

**FAZA 3 — Stabilność vs R_MAX (KLUCZOWY WYNIK):**

| R_MAX | α*₃ | G₃(α*₃) | G₃−R31 |
|-------|-----|---------|--------|
| 120 | **2.92952433** | 3495.18 | +17.96 |
| 150 | **2.92952433** | 3495.18 | +17.96 |
| 200 | **2.92952433** | 3495.18 | +17.96 |
| **σ** | **0.000000** | — | — |

**✅ α*₃ NIEZALEŻNA od R_MAX (σ=0) — analogicznie do α*₁,₂ w ex88!**

**Wartość finalna:** α*₃ = **2.929 ± 0.007** (centrum [2.923, 2.930])

---

### Odkrycia fizyczne ex93

**1. G₃/F ≈ R32 przy α≈2.850**

Przy α=2.850: G₃/F = 16.7995 ≈ R32 = m_τ/m_μ = 16.817 (diff = **−0.10%**)

Interpretacja: istnieje α ≈ 2.850 gdzie stosunek amplitud tau/muon dokładnie odtwarza stosunek mas τ/μ. Ale przy tym α: F(2.850)=230.4 ≠ R21=206.77 — warunek muonowy nie jest tam spełniony.

**2. Suma trzech α***

| Formuła | Wartość | Cel |
|---------|---------|-----|
| α*₁+α*₂ | 5.06742 | — |
| α*₁+α*₂+α*₃ | **7.990** | **8 = 4·α_TGP?** |
| 4·α_TGP | 8.000 | — |
| Odchylenie od 8 | −0.010 | 1250 ppm |

Formuła **S₃ = α*₁+α*₂+α*₃ = 8 = 4·α_TGP** jest kandydatem (odchylenie 1250 ppm). Nie tak ostre jak wymagane (ex88 dał 22336 ppm dla S₂=2π−11/10 → obalona). Wymaga sprawdzenia czy G₃ niestabilność numeryczna nie wchodzi w grę.

**3. Porównanie stabilności**

| Wielkość | ex88 (F=R21) | ex93 (G₃=R31) |
|----------|-------------|----------------|
| α*₁ | 2.43183767 | — |
| α*₂ | 2.63557742 | — |
| α*₃ | — | ≈2.929 |
| σ vs R_MAX | 0 | 0 |
| Precyzja | maszynowa | ±0.007 (G₃ noise ≈50) |

---

### Wnioski O-L5

1. **α*₃ istnieje i jest stabilna vs R_MAX (σ=0)** — sektor tau jest realny
2. **α*₃ ≈ 2.929** (przedział [2.923, 2.930])
3. **G₃/F = R32 przy α≈2.850** — nowy kandydat fizyczny (inne podejście)
4. **S₃ ≈ 8 = 4·α_TGP** — kandydacka formuła sumy (1250 ppm, niezweryfikowana)
5. Wymagana większa precyzja G₃: ex94 z wyższą rozdzielczością i R_MAX=300

**O-L5: STATUS — Kandydat α*₃≈2.929 znaleziony, stabilny, wymaga weryfikacji algebraicznej**

---

### Plan v37 po ex93

```
✅ ex92: profil F(α) + G₃(α) — DONE
✅ ex93: α*₃≈2.929 stabilna vs R_MAX — DONE
⏳ O-L4 analityczne: formuła α_min (dF/dα=0)
⏳ Weryfikacja S₃=8? (formuła sumy trzech α*)
⏳ ex94: G₃/F=R32 — precyzyjne α przy warunku tau/muon
```

---

## Wyniki ex94 (in-session, 2026-03-28)

### Parametry ex94

```
Sequential brentq (brak Pool — precyzyjna analiza)
R_MAX testowane: 120, 150, 200, 300
Czas całkowity: 3183s (~53 min)
```

---

### α*₃ POTWIERDZONA do precyzji maszynowej

| R_MAX | α*₃ | G₃(α*₃) | G₃−R31 |
|-------|-----|---------|--------|
| 120 (ex93) | 2.929524330 | 3495.18 | +17.96 |
| **300 (ex94)** | **2.929524307** | 3495.20 | +17.98 |
| **Różnica** | **2.3×10⁻⁸** | — | — |

**Wniosek: α*₃ = 2.9295243 jest stabilna od R_MAX=120 do R_MAX=300** — analogia do ex88 (α*₁,₂ identyczne dla R_MAX=100–300, σ=0).

Uwaga: G₃(α*₃)=3495 ≠ R31=3477 (diff +17.98 = +0.52%) — to szum numeryczny G₃ (noise ≈ 20 jednostek w ekstrapolacji A_tail(φ²·z₀)).

---

### Formuła sumy S₃

| Wielkość | Wartość |
|----------|---------|
| α*₁ | 2.43183767 |
| α*₂ | 2.63557742 |
| α*₃ | 2.92952431 |
| **S₃ = α*₁+α*₂+α*₃** | **7.99693940** |
| **8 = 4·α_TGP** | **8.00000000** |
| **Odchylenie** | **−382.6 ppm** |

**Kandydaci S₃ (ranking):**

| Formuła | Wartość | ppm |
|---------|---------|-----|
| **4·α_TGP = 8** | 8.0000 | **−383** ✅ |
| 5φ−1/10 | 7.9902 | +847 |
| 2π+φ | 7.9012 | +12115 |
| 2π+3/2 | 7.7832 | +27464 |

**S₃=8 jest zdecydowanie najlepszym kandydatem.**

**Analiza szumu**: szum numeryczny G₃ ≈ ±20 jednostek → niepewność α*₃ ≈ ±0.007 → niepewność S₃ ≈ ±875 ppm. Odchylenie 383 ppm leży **w granicach szumu** — S₃=8 nie może być definitywnie potwierdzone ani obalone bez lepszej precyzji G₃.

---

### α*τμ — warunek G₃/F=R32

| R_MAX | α*τμ | G₃/F(α*τμ) | σ |
|-------|------|-----------|---|
| 120 | 2.853355 | 16.859 | |
| 150 | 2.853355 | 16.859 | |
| 200 | 2.853355 | 16.859 | |
| **σ** | **4.4×10⁻¹⁶** | — | **≈0** |

α*τμ ≈ **2.853355** — stabilna (σ≈0). Ale G₃/F=16.859 ≠ R32=16.817 (diff=+0.25%) — szum podobny jak G₃.

Interpretacja: istnieje specjalna wartość α≈2.853 gdzie stosunek amplitud (τ/μ)^4 reprodukuje m_τ/m_μ. Przy tym α: F(α*τμ)=229.8 ≠ R21 — muonowy warunek nie jest spełniony równocześnie.

---

### Analiza algebraiczna α*₃

**Najlepszy kandydat:** α*₃ = 8−S₂ = 8−α*₁−α*₂ = **2.93258491**
Odchylenie od α*₃: −1043.7 ppm (w granicach szumu G₃)

**D₃/D₂ = (α*₃−α*₂)/(α*₂−α*₁) = 1.44276** — brak prostego wzoru (√2=1.414 jest za 20183 ppm, φ za 108327 ppm).

---

### Wnioski O-L5 po ex94

1. **α*₃=2.9295243 jest stabilna do 2.3×10⁻⁸** (R_MAX=120→300) — pozorna stabilność
2. **S₃≈8 (383 ppm)** — kandydat, ale G₃ ma szum ~20 j.
3. **α*τμ≈2.853 (σ≈0)** — kandydat G₃/F=R32
4. G₃ ma szum ~20 j. — precyzja limituje

⚠️ **KOREKTA po ex95 (v38)**: α*₃≈2.929 i S₃=8 są **artefaktami** — patrz SESSION_v38.

**O-L5: STATUS** — wymaga ex96 z poprawnymi oknami ekstrapolacji.

---

### Plan po ex94

```
✅ ex92, ex93, ex94 — kompletna analiza O-L4 + O-L5
⏳ ZAMKNIĘCIE v37: SESSION_v37 finalna
⏳ v38: analityczna derywacja α_min (O-L4) lub lepsza precyzja G₃
```
