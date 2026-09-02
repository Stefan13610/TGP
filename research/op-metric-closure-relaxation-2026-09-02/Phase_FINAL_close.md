---
title: "Phase_FINAL_close — zamknięcie: Q-PASS-NUCLEATION wg litery (sol×C-BAR: nukleacja DOLNA zbieżna siatka×dt, N_det=1±0) przy pełnym ROZJEŹDZIE PRIMARY↔C-BAR: wariant metryczny w(ψ) załamuje się 10/10 (biegun ψ=4/3 PRZYCIĄGA, bo U_b(g_ceil)=−0.0219<0 — kanoniczne U nie znika na granicy metryki, w przeciwieństwie do korpusowego V_M9.1''); P1 3/3 PASS (po korekcie macierzy P1c)"
date: 2026-09-02
type: phase-final-close
tgp_owner: research/op-metric-closure-relaxation-2026-09-02
status: CLOSED
verdict: "Q: Q-PASS-NUCLEATION wg litery LOCKa §3 — para sol×C-BAR: NUCLEATION-DN na obu siatkach (h=0.025/0.0125) i obu dt/2, t₀=2.0 wszędzie, N_det=1±0 (zgodność ±1 w 4/4); obiekt = kula rdzeniowa r∈[0,~10] (rdzeń solitonu inwertuje w dół na podłogę QB-2, min g→0.549), objętość zewnętrzna wspina się do studni barierowej g=1.2541=g_ceil+0.0994. DESKRYPTYWNIE: to kreacja POJEDYNCZEGO obiektu podprogowego, nie kaskada mnożenia (przyrost N: 0→1, dalej stały). ROZJAZD PRIMARY↔C-BAR totalny: PRIMARY (w(ψ)=ψ/(4−3ψ) z eq:vol-element-M911) = BREAKDOWN 10/10 (gen t=8.72/8.77 identycznie dla 3 podłóg — zerowa czułość na g_floor; sol/lat t≤0.06 — starty częściowo poza domeną metryki: sol g₀=1.4651, lat g_max=1.4734>g_ceil); przyczyna strukturalna (analityczna, MD §2): U_b(g_ceil)=U(g_ceil)−U(1)=−0.0219<0 ⟹ w·U_b→−∞ przy ψ→4/3 — biegun metryczny PRZYCIĄGA jednorodny dryf zamiast go odpychać („dynamiczna bariera" LOCKa §0 działa tylko na dodatni człon gradientowy); kanoniczne U=g⁷/7−g⁸/8 NIE znika w ψ=4/3, w przeciwieństwie do korpusowego V_M9.1''=−γψ²(4−3ψ)²/12 (podwójne zero), którego iloczyn w·V=−γψ³(4−3ψ)/12 jest skończony. lat×C-BAR: STATIONARY jednorodne g≡0.5354=g_floor−δ (całość na podłodze; wkład Q-FAIL-owy pary). P1: P1a 6/6 PASS (dryf 0.0 ≤1e−10, oba warianty, 3 geometrie), P1b 2/2 PASS (BREAKDOWN t=2.750/3.130 — dokładna reprodukcja poprzednika 2.75/3.13), P1c 4/4 PASS po korekcie macierzy (Phase_correction_note_p1c_matrix.md: lat×PRIM wyłączony z KONTROLI — start poza domeną metryki, fakt z npz; pierwotny FAIL zachowany). Phase 3 widmowe NIE WYKONANE (warunek Q-PASS-STATIC niespełniony); charakterystyka kaskady bez progów wykonana (Phase3_output.txt). Sanity g_max vs g_ceil: struktury C-BAR żyją TUŻ NAD granicą (g_ceil+0.0994, studnia bariery κ=100: κδ²=|U′| ⟹ δ=0.0994 — samozgodne); w PRIMARY brak stanów końcowych."
anti_lakatos_lock: PRESERVED
tags: [metric-closure, two-sided-closure, nucleation, metric-pole-attracts, primary-cbar-split, upward-runaway, closed]
related:
  - "[[Phase0_balance.md]]"
  - "[[Phase_method_decisions.md]]"
  - "[[Phase_correction_note_p1c_matrix.md]]"
  - "[[README.md]]"
  - "[[NEEDS.md]]"
  - "[[../op-metametric-boundary-2026-09-01/Phase_FINAL_close.md]]"
  - "[[../op-blocked-soliton-bang-2026-07-04/README.md]]"
  - "[[../../core/sek08a_akcja_zunifikowana/sek08a_akcja_zunifikowana.tex]]"
---

# Phase FINAL — zamknięcie cyklu op-metric-closure-relaxation

**Status: CLOSED-EXECUTED (2026-09-02, jedna sesja: LOCK → method_decisions
→ Phase 1 → Phase 2 → Phase 3-kaskada → zamknięcie).** Kryteria LOCKa
(`Phase0_balance.md` §2–3, §6) stosowane DOSŁOWNIE; zero zmian
progów/detektorów/seedów po starcie (jedna korekta macierzy kontroli — §4).

---

## 0. Werdykty

| Pytanie | Werdykt | Jedno zdanie |
|---|---|---|
| P1a (próżnia, oba warianty) | **PASS 6/6** | dryf 0.0 (gate ≤1e−10) w PRIMARY i C-BAR, 3 geometrie — referencja próżniowa U−U(1) (MD §1) daje dokładną stacjonarność |
| P1b (ciągłość z poprzednikiem) | **PASS 2/2** | bez domknięć: BREAKDOWN t=2.750 (h=0.025) / 3.130 (h=0.0125) — identyczne z poprzednikiem (2.75/3.13) |
| P1c (sektor stabilny, kontrola detektorów) | **PASS 4/4** (po korekcie macierzy — §4) | wszystkie → próżnia STATIONARY, ZERO alarmów OBU detektorów (w tym przy zasianych N_seed_dn=1, N_seed_up=1) |
| **Q** (stan metametryczny / kaskada przy obustronnym domknięciu) | **Q-PASS-NUCLEATION** (litera §3) | sol×C-BAR: nukleacja DOLNA zbieżna (obie siatki × dt i dt/2; t₀=2.0; N_det=1±0) |
| Phase 3 (widmo) | **NIE WYKONANE** | LOCK dopuszcza tylko przy Q-PASS-STATIC; wykonano wariant kaskadowy (deskryptywny, bez progów) |

**Zastrzeżenie deskryptywne (obowiązkowe, litera „rozjazd raportować"):**
pozytyw nukleacyjny niesie wyłącznie wariant kontrolny **C-BAR (goły
sufit C²)**. Wariant **PRIMARY (fizyka metryki M9.1'')** załamał się
we WSZYSTKICH 10 biegach — patrz §3b. „Kreacja na granicy metryki" ma
nośnik numeryczny w klasie „sufit", NIE w klasie „czynnik metryczny
z kanonicznym U".

## 1. Wejścia i domknięcia (rejestr; MD §1, §5)

- **g_ceil = √(4/3) = 1.1547005** [INPUT dokumentacyjny; CYTAT sek08a:
  „√(−g_eff) = c₀ψ/(4−3ψ) (eq:vol-element-M911)"; metryka
  „ds² = −c₀²(4−3ψ)/ψ dt² + ψ/(4−3ψ)δᵢⱼdxⁱdxʲ (eq:metric-M911-canonical)"];
  g_thr_up = 1.0773503.
- g_floor PRIMARY = 0.5458938 (QB-2 0.298); wrażliwość genezy
  {0.4438468, 0.5753260}; κ=100 FROZEN; seed=20260902; amp=1e−3;
  g₀_μ=1.4650974; β=γ=1; tło 2π z npz (READ-ONLY, mtime 2026-08-31
  21:41 niezmieniony po wszystkich odczytach).
- PRIMARY: E=∫w(ψ)[½K|∇g|²+U_b]dx, U_b=U−U(1)+V_fl (referencja
  próżniowa WYMUSZONA literą P1a — bez niej δE/δg|₁=w′(1)·2·U(1)=1/7≠0;
  MD §1); jawna wariacja z członem w′(ψ)·2g·[…] — MD §2; regularizacja
  bieguna TYLKO powyżej g_ceil−1e−6 (płaska kontynuacja W̃, W̃′=0 —
  dokładny gradient energii zregularyzowanej; MD §3).
- C-BAR: w≡1 + V_ceil=(κ/3)(g−g_ceil)³ dla g>g_ceil.

## 2. Phase 1 (Phase1_output.txt)

- **P1a 6/6:** ‖g−1‖∞ = 0.0 przez t=10 (PRIM/CBAR × radial h=0.0125,
  3D 2π N=32, 3D 4π N=48).
- **P1b 2/2:** bez domknięć BREAKDOWN t=2.750/3.130 ∈ [2.7,3.2] —
  ciągłość z poprzednikiem DOKŁADNA (2.75/3.13).
- **P1c 4/4** (lat CBAR t=19, gen PRIM t=7, gen CBAR t=7, sol CBAR t=16;
  wszystkie ‖g−1‖∞ ≤ 9e−9, zero alarmów; detektory zdolne do FAIL
  i czyste). Pre-korekta: lat×PRIM = BREAKDOWN t=0.04
  (`Phase1_output_pre_correction.txt`; §4).

## 3. Phase 2 — RACHUNEK CENTRALNY (Phase2_output.txt)

### 3a. Tabela relaksacji (14 biegów głównych + 2 kontrole dt/2)

| Start | Wariant | Siatka | Status | t | Stan końcowy / ostatni zdrowy | g_max vs g_ceil |
|---|---|---|---|---|---|---|
| geneza f1 | PRIM | N=48 | BREAKDOWN | 8.72 | t=8: g∈[1.0028,1.0245], N_dn=N_up=0 | (brak stanu) |
| geneza f1 | PRIM | N=64 | BREAKDOWN | 8.77 | j.w. (identyczna trajektoria) | (brak stanu) |
| geneza f2 | PRIM | N=48/64 | BREAKDOWN | 8.72/8.77 | IDENTYCZNIE jak f1 | — |
| geneza f3 | PRIM | N=48/64 | BREAKDOWN | 8.72/8.77 | IDENTYCZNIE jak f1 | — |
| soliton | PRIM | h=0.025 | BREAKDOWN | 0.06 | start poza domeną (g₀=1.4651>g_ceil) | — |
| soliton | PRIM | h=0.0125 | BREAKDOWN | 0.06 | j.w. | — |
| sieć 2π | PRIM | N=32 | BREAKDOWN | 0.04 | start poza domeną (g_max=1.4734, 6.6% obj.>g_ceil) | — |
| sieć 2π | PRIM | N=48 | BREAKDOWN | 0.04 | j.w. (g_max=1.4688) | — |
| **soliton** | **C-BAR** | h=0.025 | **NUCLEATION-DN** | t₀=2.0 | g∈[0.5495,1.2541], N_dn=1/seed0 | **+0.0994** |
| **soliton** | **C-BAR** | h=0.0125 | **NUCLEATION-DN** | t₀=2.0 | g∈[0.5495,1.2541], N_dn=1/seed0 | **+0.0994** |
| soliton | C-BAR | dt/2 ×2 | NUCLEATION-DN | t₀=2.0 | N_det=1 (zgodność ±1: 4/4) | +0.0994 |
| sieć 2π | C-BAR | N=32 | STATIONARY | 18 | **jednorodne g≡0.5354** (=g_floor−δ) | −0.6193 |
| sieć 2π | C-BAR | N=48 | STATIONARY | 19 | j.w.; zbieżność podsiatki 3.8e−9 | −0.6193 |

INCOMPLETE: 0 biegów. Nukleacji GÓRNEJ nie potwierdzono nigdzie
(w PRIMARY regiony górne pojawiały się tuż przed załamaniem — okno
10 j.cz. niepotwierdzone; N_up=1 w chwili załamania gen/lat/sol).

### 3b. ROZJAZD PRIMARY↔C-BAR (deskryptywnie, litera drzewa §6)

**Totalny: PRIMARY 10/10 BREAKDOWN (załamania NIE-nukleacyjne),
C-BAR 6/6 zdarzeń regularnych (4× nukleacja zbieżna + 2× stacjonarność).**
Struktura rozjazdu:

1. **Biegun metryczny PRZYCIĄGA zamiast odpychać.** Fakt analityczny
   (MD §2, potwierdzony numerycznie): U_b(g_ceil) = U(g_ceil) − U(1)
   = −0.004039 − 0.017857 = **−0.021896 < 0**, więc gęstość
   w(ψ)·U_b → **−∞** przy ψ→4/3. „Nieskończony koszt objętościowy"
   z LOCKa §0 mnoży UJEMNĄ gęstość kanonicznego U — dla jednorodnego
   dryfu w górę biegun jest studnią bez dna, nie barierą; bariera
   działa wyłącznie na (dodatni) człon gradientowy w·½K|∇g|².
   Geneza: pole niemal jednorodnie dryfuje w górę (t=8:
   g∈[1.0028,1.0245], zero struktury nad progami) i wpada w biegun
   t≈8.7 — identycznie dla trzech podłóg (zerowa czułość na g_floor,
   jak u poprzednika) i zgodnie między siatkami (8.72 vs 8.77).
2. **Kanoniczne U nie jest potencjałem pary metrycznej.** Korpusowe
   V_M9.1''(ψ) = −γψ²(4−3ψ)²/12 ma PODWÓJNE ZERO w ψ=4/3 i iloczyn
   w·V = −γψ³(4−3ψ)/12 jest skończony na granicy; kanoniczne
   U = g⁷/7 − g⁸/8 zera tam nie ma. Zalockowany funkcjonał PRIMARY
   (w × kanoniczne U) łączy elementy dwóch niekompatybilnych domknięć —
   to jest strukturalna PRZYCZYNA załamań, nie artefakt siatki
   (obserwacja deskryptywna; żadna zmiana modelu nie została wykonana).
3. **Starty sol/lat leżą częściowo POZA domeną metryki** (ψ>4/3:
   sol g₀=1.4651, lat g_max≈1.47; fakt z danych READ-ONLY) — w PRIMARY
   załamanie natychmiastowe (t≤0.06) w strefie regularizacji
   (W̃_clip≈1.9e5 przy zalockowanym dt=0.01). Obserwacja LOCKa §0
   („tła bloch g_max=1.141–1.143 tuż pod g_ceil") dotyczyła tła 1D;
   tła 3D granicę PRZEKRACZAJĄ.
4. **C-BAR pokazuje, co robi obustronne domknięcie, gdy sufit trzyma:**
   soliton → inwersja rdzenia W DÓŁ na podłogę QB-2 (nukleacja obiektu
   podprogowego) + wspinaczka objętości zewnętrznej do studni
   barierowej g = g_ceil+0.0994 (samozgodne: κδ²=|U′(g*)| ⟹ δ=0.0994);
   sieć 2π → kolaps całości NA PODŁOGĘ (jednorodne g≡0.5354=g_floor−δ,
   δ=0.011 — zgodne z MD poprzednika §2). Sanity pre-rejestrowane:
   struktury końcowe C-BAR żyją TUŻ NAD granicą metryczną (nie tuż pod,
   jak tła bloch) — w studni kary, głębokość ~κ.

## 4. Korekty / incydenty / higiena (anti-Lakatos)

- ✓ LOCK + method_decisions zamknięte PRZED kodem; progi/detektory/
  seedy/siatki nietknięte po starcie; zakaz strojenia g_ceil/g_floor/κ
  dotrzymany; C-BAR nieusunięty; INCONCLUSIVE-owe załamania NIE
  reinterpretowane jako pozytyw.
- ✓ **Korekta 1** (`Phase_correction_note_p1c_matrix.md`, zapisana PRZED
  użyciem skorygowanego wyniku): macierz KONTROLI P1c włączała
  lat×PRIMARY wbrew własnej zamrożonej zasadzie domenowej MD §6 —
  na fałszywej przesłance faktograficznej (obserwacja LOCKa o g_max
  tła dotyczyła tła 1D; tło 3D ma g_max=1.4734>g_ceil, zweryfikowane
  w npz READ-ONLY). Pierwotny output zachowany
  (`Phase1_output_pre_correction.txt`: lat×PRIM BREAKDOWN t=0.04).
  Macierz P2 NIETKNIĘTA (lat/sol×PRIMARY biegły w Phase 2 wg litery
  LOCKa — wyniki w tabeli §3a).
- ✓ Korekta deskryptywnej uwagi końcowej w `Phase3_cascade.py`
  (geometria obiektu opisana z danych: kula rdzeniowa, nie powłoka) —
  tekst opisowy, zero progów/kryteriów; output przegenerowany.
- ✓ Referencja próżniowa U−U(1) w funkcjonale PRIMARY: decyzja
  ZAMROŻONA w MD §1 PRZED obliczeniami (wymuszona literą P1a;
  konwencja linii), nie post-hoc.
- ✓ Tła npz READ-ONLY: mtime 2026-08-31 21:41 niezmieniony po
  wszystkich odczytach (weryfikacja w każdym starcie `lat` i na końcu
  sesji). Rdzeń `.tex` NIETKNIĘTY; STATE.md nieedytowane; git
  nieużywany; katalogi innych cykli tylko odczyt; pełne ścieżki bez
  `cd`; `ls` po każdym zapisie.
- ✓ Nazewnictwo: LOCK §5 przewidywał `Phase3_spectrum.py` (warunkowy —
  warunek Q-PASS-STATIC NIE zaszedł ⟹ nie utworzony); wariant
  kaskadowy zrealizowany jako `Phase3_cascade.py` → `Phase3_output.txt`
  (LOCK nie specyfikował nazwy pliku wariantu kaskadowego).
- Środowisko: CPython 3.14.2, numpy 2.4.3, scipy 1.17.1, sympy 1.14.0.

## 5. Phase 3 — charakterystyka kaskady (Phase3_output.txt; bez progów)

- N_obj(t) detektora dolnego: 0 → 1 (t=2.0) → stały (brak mnożenia);
  przyrosty na j.cz.: pojedyncza kreacja. Detektor górny: N_up
  1 → 0 (t=1: rdzeń chwilowo pod progiem po pierwszym zrzucie) →
  1 (t=3+, objętość zewnętrzna nad progiem); na h=0.025 w t=12
  przejściowo 2 regiony górne (pierścień r∈[19.2,21.5] + zewnętrze) —
  na h=0.0125 pojedynczy region (różnica siatkowa raportowana wprost;
  detekcyjny werdykt DN od niej nie zależy).
- Rozmiary obiektu dolnego (j. komórek, h=0.025 / h=0.0125):
  t₀: 44/89; t₀+5: 185/370; t₀+10: 398/813 — obiekt ROŚNIE
  (poszerzająca się kula rdzeniowa), fizyczny zasięg zgodny między
  siatkami (r∈[0,9.94] vs [0,10.16]).
- Energia: monotoniczny spadek (0.138 → −3060 do t=12) — sanity flow;
  stan w t=12 wciąż ewoluuje (bieg zakończony potwierdzeniem okna
  detektora, litera LOCKa).
- **Deskryptywnie:** „kaskada" w sensie mnożenia obiektów NIE wystąpiła;
  wystąpiła pojedyncza, zbieżna siatka×dt kreacja obiektu podprogowego
  (inwersja rdzenia solitonu na podłogę QB-2) przy jednoczesnym
  osadzeniu tła na granicy górnej (studnia bariery). Pre-rejestrowany
  pozytyw autora zaszedł w minimalnej krotności (N_det=1).

## 6. Odczyt (deskryptywnie, bez claimów poza klasą zbadaną)

1. **Hipoteza „kreacja na granicy metryki" dostaje nośnik numeryczny
   wyłącznie w klasie SUFIT (C-BAR):** przy twardym górnym domknięciu
   sektor tachionowy nie ucieka, lecz TWORZY obiekt (rdzeń → podłoga)
   i osadza tło na granicy (g_ceil+δ_κ). To pierwszy pozytyw
   nukleacyjny linii (poprzednik: 18/18 załamań bez nukleacji).
2. **Czynnik metryczny w(ψ) z KANONICZNYM U nie domyka:** biegun
   przyciąga (U_b(4/3)<0). Naturalny kandydat domknięcia metrycznego
   pozostaje w korpusie: para (w, V_M9.1'') z podwójnym zerem V
   w ψ=4/3 — poza zalockowanym modelem tego cyklu (user-gate, NEEDS).
3. Wynik genezowy: w pudle L=4π (k_min=0.5, pasmo tachionowe częściowo
   w pudle) struktura przestrzenna NADAL nie zdążyła powstać przed
   ucieczką jednorodną (t=8: brak regionów nad progami) — w PRIMARY;
   klasa C-BAR dla genezy nie była zalockowana (macierz LOCKa §3).

## 7. Pliki cyklu

`Phase0_balance.md` (LOCK) · `Phase_method_decisions.md` (FROZEN) ·
`Phase1_gate.py` → `Phase1_output.txt`
(+ `Phase1_output_pre_correction.txt`) ·
`Phase_correction_note_p1c_matrix.md` ·
`Phase2_two_sided_relax.py` → `Phase2_output.txt`
+ `Phase2_relaxed_states.npz` (16 stanów) + `Phase2_results/`
(json/npz per bieg) + `Phase2_results_batch*.log` ·
`Phase3_cascade.py` → `Phase3_output.txt` (`Phase3_spectrum.py`
NIE UTWORZONY — warunek LOCKa niespełniony) · `NEEDS.md` (user-gated) ·
`README.md` (log).

## 8. Mapowanie na drzewo decyzyjne LOCKa §6

**Q-PASS-NUCLEATION** → NEEDS: hipoteza „kreacja na granicy metryki"
z nośnikiem numerycznym (z zastrzeżeniem klasy C-BAR); PILNE: geneza
Γ+s_i (poziom 0) + reinterpretacja Q-FAIL-i kanonicznych jako
relaksacji (user-gate, dopiski core). Rozjazd PRIMARY↔C-BAR
raportowany deskryptywnie (§3b) — dosłownie wg drzewa.
