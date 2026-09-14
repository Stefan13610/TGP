---
title: "Phase0_balance — LOCK: klasa dynamiki pary metrycznej (w, V_M9.1'') — czy prawo zachowania substratu (A: Cahn–Hilliard) i bezwładność z lapse M9.1'' (B: dynamika 2. rzędu) tłumią relaksację i podtrzymują/wytwarzają strukturę w strefie słabego pola (ψ<2/3)?"
date: 2026-09-02
type: phase0-lock
tgp_owner: research/op-dynamics-class-M911-2026-09-02
status: PHASE0-LOCKED
computations_performed: ZERO
authorization: "User 2026-09-02: „Ok A + B działaj ;)" (wybór wariantów A+B z rekomendacji po zamknięciu op-metric-pair-M911; hipoteza usera: kreacja w przedmetrycznym stanie zaburzonego substratu przy słabym ψ)"
anti_lakatos_lock: ACTIVE
related:
  - "[[README.md]]"
  - "[[../op-metric-pair-M911-2026-09-02/Phase_FINAL_close.md]]"
  - "[[../op-metric-pair-M911-2026-09-02/NEEDS.md]]"
  - "[[../../core/sek08a_akcja_zunifikowana/sek08a_akcja_zunifikowana.tex]]"
---

# Phase 0 — LOCK cyklu `op-dynamics-class-M911`

**ZERO obliczeń wykonanych przed zapisaniem tego dokumentu.**

## 0. Pytanie i dziedziczona diagnoza

Poprzednik (op-metric-pair-M911): para korpusowa (w, V_M9.1'', K=ψ⁴)
jest SAMODOMKNIĘTA (Q-A-PASS), ale w klasie relaksacyjnej gradient flow
(dynamika 1. rzędu, NIEzachowana, bez bezwładności) wszystko spływa do
próżni ψ≡1 w czasie ~1/𝒰″(1)=O(10) (Q-B-FAIL). Ustalone analitycznie
(Phase 1 poprzednika): strefa spinodalna 𝒰″=γψ(3ψ−2)<0 istnieje
WYŁĄCZNIE przy słabym polu ψ<2/3; zalockowane starty poprzednika jej
nie próbkowały (ψ_min=0.837). Szybkość relaksacji była własnością KLASY
dynamiki: (a) parametr niezachowany (Allen–Cahn: mod k=0 gaśnie w tempie
𝒰″(1)=1), (b) zero bezwładności, (c) zero wagi metrycznej czasu (lapse).

**Pytania binarne (obydwa w tej samej zamrożonej parze (w, V_M9.1'', K)):**
- **Q-CONS (wariant A, zachowanie substratu):** czy przy dynamice
  ZACHOWANEJ (Cahn–Hilliard: ∂ψ/∂t=∇²μ, μ=δE/δψ, ∫ψ=const) struktura
  nukleuje lub TRWA (w tym: czy start słabopolowy ψ_min<2/3 przetrwa
  jako obiekt), czy wszystko dyfunduje do jednorodnej średniej?
  Obowiązkowo deskryptywnie: ILOKROTNIE wolniejsza jest relaksacja
  względem poprzednika (czasy t_dev: pierwsze t z dev<1e−3).
- **Q-INER (wariant B, bezwładność+lapse):** czy przy dynamice
  2. rzędu z wagą czasową z metryki M9.1'' (B(ψ)=w²K=ψ⁶/(4−3ψ)²,
  zachowawcza, H=const) powstają/trwają zlokalizowane struktury
  (oscylony; pre-rejestrowany POZYTYW autora dla startu słabopolowego)
  lub nukleacja — czy pole dyspersuje do quasi-jednorodności?

Uwaga interpretacyjna zalockowana: FAIL obu wariantów = wynik ważny —
mówi, że quietyzm pary metrycznej NIE jest artefaktem klasy
przetłumionej i kieruje hipotezę kreacji do genezy Γ+s_i / sprzężenia
z materią (jak drzewo poprzednika). PASS któregokolwiek = pierwszy
nośnik „tłumienia relaksacji" wewnątrz czystego sektora.

## 1. Model ZAMKNIĘTY (formy dziedziczone, dynamiki zalockowane TERAZ)

**Para (FROZEN, cytaty w MD poprzednika §1 — dziedziczone dosłownie):**
w(ψ)=ψ/(4−3ψ), V=−γψ²(4−3ψ)²/12, K=K_geo ψ⁴; K_geo=γ=1; odczyt B
kinetyki (𝒦=ψ⁴, w·g^ij≡1); 𝒰=w·V=γ(ψ⁴/4−ψ³/3), 𝒰′=γψ²(ψ−1);
E[ψ]=∫[½𝒦|∇ψ|²+𝒰]dx; dziedzina ψ∈(0,4/3). ZAKAZ modyfikacji form.

**Wariant A (Q-CONS):** ∂ψ/∂t = ∇²μ, μ = δE/δψ, mobilność M≡1
[FROZEN — izolacja czystego efektu prawa zachowania; alternatywy
M∈{w,ψ⁴} odnotowane, NIErealizowane]. Zachowanie ∫ψ dx (gate).
E maleje (dE/dt=−∫|∇μ|² — sanity).

**Wariant B (Q-INER):** z akcji M9.1'' (człon czasowy w·K·|g^tt|ψ̇²,
|g^tt|=ψ/(c₀²(4−3ψ)), c₀=1): gęstość kinetyczna czasowa ½B(ψ)ψ̇²,
**B(ψ)=w²K=ψ⁶/(4−3ψ)²**; EOM: **B ψ̈ + ½B′ ψ̇² = −δE/δψ**;
H=∫[½Bψ̇²+½𝒦|∇ψ|²+𝒰]dx zachowane (gate dryfu). Start: ψ̇≡0.
Bez tłumienia (η=0 — granica przetłumiona już zbadana przez poprzednika).

**Starty (zalockowane 3, wszystkie 3D periodyczne):**
- (i) **geneza:** ψ=1+szum, L=4π, N∈{32,48}, seed=20260904, amp=1e−3,
  konstrukcja pasmowa dziedziczona (|n_i|≤8, te same współczynniki
  na obu siatkach);
- (ii) **dip słabopolowy (NOWY, testuje strefę spinodalną):**
  ψ = 1 − 0.45·exp(−r²/(2σ²)), σ=1.5, centrum pudła, L=4π, N∈{32,48};
  ψ_min=0.55 < 2/3 (wewnątrz strefy 𝒰″<0);
- (iii) **sieć 2π:** klucze `2pi__A1.0__N{32,48}` z npz READ-ONLY
  poprzednika op-3d-canonical-lattice; procedura skalowania FROZEN
  dziedziczona (ψ_raw=g², ψ₀=1+s(ψ_raw−1), s=0.30/(maxψ_raw−1)
  per siatka, ψ_max=1.30).

**Macierz:** 3 starty × 2 siatki × 2 warianty = 12 biegów głównych;
dt/2 warunkowo „przy zdarzeniach" (nukleacja / kandydat
STRUCTURE/OSCILLON / BOUNDARY / BREAKDOWN) — obie siatki pary.

**Detektory (dziedziczone DOSŁOWNIE):** dolny ψ<5/6, górny ψ>7/6;
ndimage.label + sklejanie periodyczne; N_seed z t=0; nukleacja =
N>N_seed utrzymane ≥10 j.cz.; zbieżność: kierunek + N_det ±1 na
(2 siatki × dt, dt/2). UWAGA pre-rejestrowana: dip ma N_seed_dn=1
(0.55<5/6) — trwanie zasianego obiektu NIE jest nukleacją (od tego
jest kryterium STRUCTURE/OSCILLON).

**Obsługa granic (dziedziczona):** zero podłóg/barier; pas ψ>4/3−1e−6
= BREAKDOWN-BOUNDARY; min ψ<1e−6 = BREAKDOWN-BOUNDARY-LOWER;
niefinityczność = BREAKDOWN. **Pre-rejestracja sztywności (wariant B):**
B→0 przy ψ→0 (zanik bezwładności — lokalna częstość |𝒰″|/B ~ 32/ψ⁵);
zejście pola pod ψ_stiff=0.12 flagowane STIFF; następujące po nim
załamanie klasyfikowane INCONCLUSIVE-STIFF (sztywność integratora,
NIE werdykt fizyczny) — rozstrzyga bieg dt/2.

**Rejestr WEJŚĆ:** seed=20260904, amp=1e−3; dip (0.45, σ=1.5);
skalowanie sieci per-siatka do 1.30; M≡1; B=w²K; dt_A=0.01,
dt_B=0.0025 (kontrole dt/2); t_max A=200, B=100; stacjonarność A:
‖ψ̇‖∞≤1e−8; okno oceny B: [t_max−20, t_max], obłożenie ≥80% próbek;
gate energii B: ε_H=|H(t)−H(0)|/max(|H_rel(0)|,0.01) ≤ 0.02;
progi detektorów 5/6, 7/6; dev strukturalny 0.05; podsiatka 5e−3;
ψ_stiff=0.12.

## 2. Fazy i kryteria (zalockowane)

### Phase 1 — bramki maszynerii (dowolny FAIL ⟹ STOP)
- **G1 (A):** próżnia ψ≡1 zostaje (t=10, dryf ≤1e−10, N=32, L∈{2π,4π});
  zachowanie masy na starcie dip N=32, t=10: |Δψ̄|/|ψ̄| ≤ 1e−13;
  E niemalejąco NIE rośnie (wzrost między próbkami ≤1e−12).
- **G2 (B):** próżnia zostaje (jw., dryf ≤1e−10; ψ̇≡0);
  test fali liniowej: ψ=1+1e−6·cos(kx), k=1 (mod m=2 w L=4π), t=20:
  częstość zmierzona vs analityczna ω=√(𝒦(1)k²+𝒰″(1))/√B(1)=√2 —
  zgodność ≤1e−3 względnie; dryf H ≤1e−3 względem energii fali;
  pochodna B′(ψ)=12ψ⁵(2−ψ)/(4−3ψ)³ vs sympy — zgodność 1e−12
  w {0.5, 1, 7/6, 1.3} (osiągalny FAIL implementacji).
- **G3 (detektory):** zasiany dip (min 0.6) i bump (max 1.3) wykryte
  1±0, czysta próżnia zero alarmów (3D N=48, L=4π).

### Phase 2 — wariant A (Q-CONS): 6 biegów + warunkowe dt/2
Do ‖ψ̇‖∞≤1e−8 / nukleacji / pasa / t_max=200. **Werdykty (litera):**
- **Q-CONS-PASS-NUCLEATION:** detektor, zbieżnie (2 siatki × dt, dt/2, ±1).
- **Q-CONS-PASS-STRUCTURE:** stan końcowy (STATIONARY lub TMAX)
  z dev=½(max−min) ≥0.05 na OBU siatkach, zbieżność podsiatkowa
  (16³) ≤5e−3, obiekt w detektorach (N_dn+N_up≥1) we WSZYSTKICH
  próbkach ostatnich 20 j.cz., potwierdzone w dt/2 na obu siatkach.
- **Q-CONS-FAIL:** wszystkie pary jednorodne (dev<0.05) na końcu —
  dyfuzja do średniej (raport: średnia vs 1; czasy t_dev vs poprzednik).
- **Q-CONS-INCONCLUSIVE:** pozostałe (w tym BOUNDARY jako kategoria
  deskryptywna; niezbieżność).
Deskryptywnie OBOWIĄZKOWO: tabela t_dev(dev<1e−3) A vs poprzednik
(7.0/16.0/16.0) — ilościowa odpowiedź „czy zachowanie tłumi relaksację".

### Phase 3 — wariant B (Q-INER): 6 biegów + warunkowe dt/2
Do nukleacji / pasa / załamania / t_max=100 (bez stopu stacjonarności —
dynamika zachowawcza). Kryteria na skalarach okna [80,100] per siatka
(porównanie POLOWE między siatkami NIE jest wymagane — dekoherencja
fazowa dynamiki falowej, pre-rejestrowane odstępstwo od konwencji
linii). **Werdykty (litera):**
- **Q-INER-PASS-NUCLEATION:** detektor, zbieżnie (jw.).
- **Q-INER-PASS-OSCILLON:** w oknie [80,100]: min dev ≥0.05 ORAZ
  obłożenie detektorowe (N_dn+N_up≥1) ≥80% próbek, na OBU siatkach,
  przy ε_H≤0.02 — potwierdzone w dt/2 na obu siatkach.
- **Q-INER-FAIL-DISPERSE:** max dev w oknie <0.05 na obu siatkach
  (pole rozproszone; dla genezy amp 1e−3 oznacza to brak wzmocnienia —
  raport wprost).
- **Q-INER-INCONCLUSIVE:** pozostałe (STIFF/DRIFT/BOUNDARY jako
  kategorie deskryptywne — NIE pozytyw).
Deskryptywnie OBOWIĄZKOWO: trajektoria ψ_min(t) dla dipa (czy strefa
spinodalna pogłębia zaburzenie), N_obj(t) obu detektorów, ε_H(t).

### Phase 4 — zamknięcie
`Phase_FINAL_close.md` (wzorzec linii), `NEEDS.md` (user-gated,
drzewo §5), `README.md`. Widma NIE liczymy (poza zakresem locka;
przy PASS-STRUCTURE/OSCILLON — charakterystyka deskryptywna obiektu:
profil, rozmiar, kontrast, okres oscylacji z serii dev(t)).

## 3. Forbidden moves
Dziedziczone z linii + (a) formy (w,V,K,𝒰) i odczyt B NIETYKANE;
(b) zakaz podłóg/barier/tłumienia η; (c) detektory i progi (5/6, 7/6,
0.05, okna, 80%, ε_H, ψ_stiff) niezmienialne po pierwszym biegu;
(d) M≡1 i B=w²K niezmienialne (to JEST przedmiot testu); (e) rdzeń
.tex/STATE/git nietykane; katalogi innych cykli tylko odczyt (npz
READ-ONLY, mtime weryfikowany); (f) INCONCLUSIVE/STIFF/BOUNDARY ≠
pozytyw; (g) rejestr WEJŚĆ flagowany w każdym wyniku.

## 4. Deliverables
`Phase_method_decisions.md` (FROZEN), `M911_common.py` (model+detektory
+starty), `Phase1_gate.py`+output, `Phase2_conserved.py`+output+npz,
`Phase3_inertial.py`+output+npz, `Phase_FINAL_close.md`, `NEEDS.md`,
`README.md`.

## 5. Drzewo decyzyjne
```text
P1 FAIL → STOP (maszyneria)
Q-CONS-PASS-* lub Q-INER-PASS-* → NEEDS: „tłumienie relaksacji /
   trwanie struktury" ma nośnik w czystym sektorze; PILNE: charakterystyka
   obiektu + user-gate: który mechanizm (zachowanie vs bezwładność)
   nośny dla hipotezy słabopolowej usera; kandydat: sprzężenie z materią
   na tym wariancie dynamiki
Q-CONS-FAIL i Q-INER-FAIL → quietyzm pary potwierdzony we WSZYSTKICH
   trzech klasach dynamiki (przetłumiona/zachowana/bezwładna) —
   hipoteza kreacji ostatecznie kierowana do genezy Γ+s_i (poziom 0)
   / sprzężenia z materią (osobny LOCK); deskryptywnie: o ile
   zachowanie/bezwładność SPOWALNIA relaksację (to też odpowiedź)
Q-*-INCONCLUSIVE → NEEDS: wniosek metodologiczny (integrator sztywny
   przy B→0 = granica techniczna klasy B) + ewentualny re-lock dt
mieszane (jeden PASS, drugi FAIL) → NEEDS: identyfikacja mechanizmu
   (user-gate na cykl następczy z tym mechanizmem)
```

---

**LOCK ZAMKNIĘTY 2026-09-02. Zmiany poniżej tej linii po starcie
obliczeń = forbidden move.**
