# Phase 0 — LOCK: op-fluctuation-extended-nbody-2026-08-31

**Status: PHASE0-LOCKED 2026-08-31. Zero obliczeń przed zamknięciem tego pliku.**
**Autoryzacja:** user 2026-08-31: „zacomituj i zajmij się
research/op-substrate-fluctuation-channel-2026-08-23/" → realizacja
zaległych research-NEEDS **N1 + N2** tamtego cyklu ([[../op-substrate-fluctuation-channel-2026-08-23/NEEDS.md]])
jako cykl-następca z własnym LOCKiem (wymóg metodologiczny: każda nowa
numeryka = własny cykl z Phase 0).
**Rejestr:** structural-emergence (inwentarz własności kanału na poziomie 0),
NIE empirical-novelty. Zero claimów obserwacyjnych/grawitacyjnych.

---

## 1. Kontekst i delimitacja

- **Poziom:** 0 (warstwa Gaussowska substratu, jak cykl-rodzic). Maszyneria
  dziedziczona z [[../op-substrate-fluctuation-channel-2026-08-23/Phase0_balance.md]]
  §3.1 (siatka L³ periodyczna, G_m przez FFT, defekt kotwiczący M-D)
  oraz Amendment A1 (estymator przypadku krytycznego: stała B̂ z fitu
  G = A·d^p + B na defekcie punktowym; propagator „connected").
- **Wynik rodzica (punkt startu):** F_fl(d) = ½ln(1−(G(d)/G₀)²) < 0
  uniwersalnie; zasięg 2μ; na krytyczności −1/d² dla defektów PUNKTOWYCH.
  Zastrzeżenie rodzica: −1/d² ≠ Newton −1/d; pytania o obiekty rozciągłe
  (N1) i N-ciałowość (N2) jawnie odroczone.
- **Zakazy delimitacyjne:** (a) zakaz dublowania `op-bloch-chain-stability-2026-08-31`
  (poziom 1, sektor dynamiczny, RÓWNOLEGŁA SESJA — nie dotykać jej katalogu);
  (b) zakaz dublowania Q1/Q2 `op-bath-two-sectors` (zamknięty);
  (c) poziom 0 wyłącznie.

## 2. Pytania binarne

> **QE (z N1):** czy dla DWÓCH OBIEKTÓW ROZCIĄGŁYCH (kule dyskretne
> promienia R z zamrożonym polem — model rdzenia solitonu) istnieje
> reżim (R, okno d), w którym kanał fluktuacyjny na krytyczności daje
> lokalny wykładnik potęgowy **p = −1.0 ± 0.15** (Newtonowskie −1/d)
> — czy też wykładnik dalekiego pola pozostaje −2 (R wchodzi tylko
> w amplitudę), a przy kontakcie jedynie stromieje?

> **QN (z N2):** czy kanał fluktuacyjny jest NIEADDYTYWNY parowo dla
> trzech defektów (ΔF₃ = F_fl(ABC) − Σ_par F_fl ≠ 0) z **uniwersalnym
> znakiem** członu trójciałowego — i jaki to znak (osłabia czy wzmacnia
> przyciąganie parowe)? Kontrast: kanał klasyczny (źródłowy) jest
> w teorii Gaussowskiej addytywny DOKŁADNIE (superpozycja) — do
> weryfikacji numerycznej jako kontrola.

## 3. Modele (ZAMKNIĘTE — zero zmian po starcie obliczeń)

### 3.1 QE — dwa klastry kotwiczące

Klaster A/B = zbiór węzłów kuli dyskretnej {x : |x − x_c| ≤ R},
R ∈ {1, 2, 3} (n_sites = 7, 33, 123). Rozdzielenie d = |x_B − x_A|
wzdłuż osi siatki (center-to-center). Pełna macierz kowariancji
Σ(d) rozmiaru (n_A+n_B)² z wpisami G(x_i − x_j); część fluktuacyjna
(v-niezależna, jak u rodzica — rozkład brzegowy Gaussa na węzłach
zamrożonych):

```
F_fl(d; R) = ½ [ ln det Σ(d) − ln det Σ_A − ln det Σ_B ]
```

(Σ_A/Σ_B = bloki własne; równoważnie ½ ln det(I − Σ_A^{-1}K Σ_B^{-1}Kᵀ),
K = blok mieszany). Nierówność Fischera ⟹ F_fl ≤ 0 zawsze — do użycia
jako kontrola maszynerii, nie jako wynik.

**Przypadek krytyczny (m=0):** wpisy macierzy z propagatora connected
G_c(r) = G(r) − B̂, gdzie B̂ z fitu S-typu rodzica (G = A·d^p + B,
defekt punktowy, okno jak u rodzica L=128: d ∈ [8,24]) — dziedziczenie
Amendment A1 wprost. Przypadki masywne: pełne G, bez stałej.

### 3.2 QN — trzy defekty punktowe

Konfiguracje na siatce (defekty punktowe, pinning):
- **T (trójkąt prostokątny):** (0,0,0), (d,0,0), (0,d,0) — odległości
  d, d, d√2;
- **C (kolinearna):** (0,0,0), (d,0,0), (2d,0,0).

```
ΔF₃ = F_fl(ABC) − [F_fl(AB) + F_fl(BC) + F_fl(AC)]
```

wszystko z log-det (3×3 vs 2×2). Kontrola klasyczna: U_S(ABC) dla
źródeł q=1 — superpozycja dokładna (odchył < 1e−12).

## 4. Kryteria (PRE-REJESTROWANE)

### QE (werdykt = odpowiedź binarna na N1; obie odpowiedzi informatywne)

- **QE-0 (tożsamość maszynerii):** dla R=0 (klaster 1-węzłowy) F_fl
  z log-det ≡ ½ln(1−(G(d)/G₀)²) co do < 1e−12 (ta sama wielkość,
  inna ścieżka kodu). FAIL ⟹ INCONCLUSIVE całości.
- **QE-1 (uniwersalność przeżywa rozciągłość):** F_fl(d;R) < 0
  i monotonicznie rosnące z d dla wszystkich R ∈ {1,2,3} we wszystkich
  oknach §5 (m=0.2 i m=0-connected).
- **QE-2 (sygnatura zasięgu przeżywa rozciągłość):** m=0.2, L=128,
  R ∈ {1,2}: fit ln|F_fl| = a − κ_D·d − 2ln d w oknie §5 daje
  κ_D = 2μ(0.2) ± 10%, R² > 0.99 (μ(m) = arccosh(1+m²/2)).
- **QE-3 (pytanie N1 — wykładnik dalekiego pola):** m=0-connected,
  L=128: slope ln|F_fl| vs ln d w oknach dalekich §5 dla R ∈ {1,2,3};
  raport slope(R). Interpretacja pre-rejestrowana: jeśli wszystkie
  slope = −2.0 ± 0.25 ⟹ „R wchodzi tylko w amplitudę" (odpowiedź NIE
  w dalekim polu); jeśli którykolwiek slope = −1.0 ± 0.15 z R² ≥ 0.99
  ⟹ odpowiedź TAK.
- **QE-4 (pytanie N1 — skan bliski kontaktu):** m=0-connected, L=128:
  lokalny wykładnik p_loc(d) (3-punktowe sąsiednie ilorazy log-log)
  na całym zakresie d od d_min = 2R+2 do końca okna dalekiego.
  **Odpowiedź TAK na N1** ⟺ istnieje ciągły przebieg ≥ 3 kolejnych d
  z |p_loc − (−1.0)| ≤ 0.15. **NIE** ⟺ brak takiego przebiegu.
  Raportować pełny profil p_loc(d) dla każdego R (bez selekcji).
- **QE-5 (amplituda, deskryptywnie):** fit A(R) w |F_fl| ≈ A(R)/d²
  (okna dalekie); porównanie z przewidywaniem pojemnościowym
  (monopol: A ~ (C_A C_B)², C_R ~ pojemność kuli) — TYLKO raport,
  bez kryterium PASS/FAIL.
- **QE-kontrole:** (a) QE-0; (b) F_fl ≤ 0 wszędzie (Fischer);
  (c) dryf slope QE-3 między L=96 a L=128 dla R=2 < 0.1 (okno wspólne);
  (d) G_m > 0 (masywne). Pad którejkolwiek ⟹ INCONCLUSIVE (maszyneria).

### QN (werdykt = QN-1 ∧ QN-2 rozstrzygane łącznie; QN-3 kontrola)

- **QN-1 (sympy exact):** rozwinięcie ΔF₃ w potęgach g_ij = G(d_ij)/G₀
  do rzędu 4 włącznie: wyprowadzić człon wiodący (oczekiwany kandydat:
  g₁₂g₁₃g₂₃, rząd 3) i jego znak SYMBOLICZNIE. Werdykt: człon wiodący
  istnieje i ma określony znak dla g_ij > 0 (TAK/NIE).
- **QN-2 (numerycznie, uniwersalność znaku):** znak ΔF₃ ten sam dla
  WSZYSTKICH punktów: m ∈ {0.2, 0-connected}, konfiguracje {T, C},
  d ∈ {4, 6, 8, 10, 12}, L=128 (20 punktów). Raport |ΔF₃|/|Σ_par F_fl|
  dla każdego punktu. Osiągalny FAIL: znak zależny od konfiguracji/d.
- **QN-3 (kontrola klasyczna):** addytywność kanału źródłowego
  U_S(ABC) − Σ_par U_S = 0 co do < 1e−12 (superpozycja Gaussowska).
  Pad ⟹ INCONCLUSIVE (błąd maszynerii, nie fizyka).
- **QN-4 (zgodność z analityką):** dla m=0.2, konfiguracja T, d=8:
  |ΔF₃_num − ΔF₃_analit(rząd≤4)| / |ΔF₃_num| < 0.05 (małe g_ij ⟹
  rozwinięcie powinno działać). Osiągalny FAIL: rozjazd rzędu 1.

## 5. Okna dopasowań (ZAMKNIĘTE przed obliczeniami)

| Test | L | m | R | okno d (oś siatki, center-to-center) |
|---|---|---|---|---|
| QE-2 fit κ_D | 128 | 0.2 | 1 | d ∈ [8, 18] |
| QE-2 fit κ_D | 128 | 0.2 | 2 | d ∈ [10, 20] |
| QE-3 slope daleki | 128 | 0-conn | 1 | d ∈ [10, 28] |
| QE-3 slope daleki | 128 | 0-conn | 2 | d ∈ [14, 28] |
| QE-3 slope daleki | 128 | 0-conn | 3 | d ∈ [16, 28] |
| QE-4 skan p_loc | 128 | 0-conn | 1/2/3 | d ∈ [2R+2, 28] |
| QE-kontrola (c) | 96 vs 128 | 0-conn | 2 | d ∈ [14, 24] |
| QN-2 | 128 | 0.2 i 0-conn | — | d ∈ {4,6,8,10,12} |

Zakaz przesuwania okien po obejrzeniu danych. Fit z R² < 0.99 w oknie
⟹ INCONCLUSIVE dla tego punktu, NIE dobieranie okna. B̂ wyznaczane
NIEZALEŻNIE na defekcie punktowym (okno rodzica [8,24], L=128) przed
policzeniem jakiegokolwiek F_fl klastrowego; jedno B̂ na L (nie per-R).

## 6. Drzewo decyzyjne

- **QE: NIE (slope −2, brak przebiegu p_loc ≈ −1):** kanał fluktuacyjny
  dwóch ciał NIE produkuje Newtonowskiego −1/d przez samą rozciągłość;
  amplituda rośnie z R, wykładnik nie. Do NEEDS: ścieżki pozostałe dla
  −1/d to N-ciałowość zbiorowa / ośrodek o skończonej gęstości defektów
  (ekranowanie) / poziom 1 — user-gate.
- **QE: TAK:** istnieje reżim −1/d — natychmiastowy user-gate (wynik
  zmieniałby bilans programu „most do grawitacji"); ZERO propagacji
  do rdzenia w tym cyklu.
- **QN: nieaddytywność z uniwersalnym znakiem:** inwentarz poziomu 0
  uzupełniony o własność N-ciałową kanału (kontrast z addytywnym
  kanałem klasycznym — zgłosić znak: osłabia/wzmacnia). Deskryptywnie
  odnotować relację do op-nbody-additivity (sektor solitonowy,
  addytywny z poprawką ~exp) — bez przenoszenia między poziomami.
- **Dowolne INCONCLUSIVE:** zgłosić wprost; żadnych wniosków z tej gałęzi.

## 7. Forbidden moves

1. Zmiana kryteriów §4, okien §5, modeli §3 po starcie obliczeń.
2. Edycje rdzenia .tex — zakaz (wszystko przez NEEDS, user-gate).
3. Dotykanie katalogu `op-bloch-chain-stability-2026-08-31` (równoległa
   sesja) i wspólnych plików koordynacyjnych (STATE.md — dopiero po
   idle tamtej sesji; NEEDS rodzica już zaktualizowane przed LOCKiem).
4. Claimy grawitacyjne/obserwacyjne — zakaz (rejestr structural-emergence).
5. Każdy skrypt musi mieć testy zdolne dać FAIL; SUMMARY cytuje liczby.
6. Wynik negatywny → wprost w Phase_FINAL_close.
7. Commit — autoryzowany przez usera („zacomituj", 2026-08-31); push
   tylko za jawną zgodą.

## 8. Fazy i deliverables

- **Phase 1 (analitycznie, sympy):** (a) rozwinięcie ΔF₃ do rzędu 4
  (QN-1); (b) daleko-polowe rozwinięcie F_fl dwóch klastrów w obrazie
  pojemnościowym (monopol) — forma amplitudy A(R) i wykładnik;
  → `Phase1_analytic.py` + `Phase1_output.txt`.
- **Phase 2 (QE numerycznie):** G przez FFT (L=96/128), klastry R∈{1,2,3},
  log-det; QE-0..QE-5 wg §4–§5; → `Phase2_extended.py` + `Phase2_output.txt`.
- **Phase 3 (QN numerycznie):** trójki T/C, ΔF₃, kontrola klasyczna;
  QN-2..QN-4; → `Phase3_nbody.py` + `Phase3_output.txt`.
- **Zamknięcie:** `Phase_FINAL_close.md`, `NEEDS.md`, `README.md`;
  wpis STATE.md po zwolnieniu pliku przez równoległą sesję.

*LOCK zamknięty 2026-08-31 przed napisaniem jakiegokolwiek kodu.*
