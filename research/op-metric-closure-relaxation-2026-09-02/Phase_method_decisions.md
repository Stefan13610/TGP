---
title: "Phase_method_decisions — funkcjonał PRIMARY z czynnikiem metrycznym w(ψ)=ψ/(4−3ψ) (cytat eq:vol-element-M911), jawna wariacja δE/δg z członem w′; gęstość energii względem próżni (wymuszone literą P1a); regularizacja bieguna wyłącznie powyżej g_ceil−1e−6 (płaska kontynuacja = dokładny gradient energii zregularyzowanej); wariant C-BAR κ=100; detektor górny (1+g_ceil)/2; siatki i starty zalockowane; decyzje ZAMROŻONE przed startem obliczeń"
date: 2026-09-02
type: method-decisions
tgp_owner: research/op-metric-closure-relaxation-2026-09-02
status: FROZEN-PRE-COMPUTE
computations_performed: ZERO
related:
  - "[[Phase0_balance.md]]"
  - "[[../op-metametric-boundary-2026-09-01/Phase_method_decisions.md]]"
  - "[[../../core/sek08a_akcja_zunifikowana/sek08a_akcja_zunifikowana.tex]]"
  - "[[../op-3d-canonical-lattice-2026-08-31/Phase_method_decisions.md]]"
---

# Decyzje metodyczne (zamrożone PRZED jakimkolwiek obliczeniem cyklu)

Kryteria/progi/punkty LOCKa (`Phase0_balance.md`) NIETKNIĘTE. Jedyne
czynności przed zamrożeniem: odczyt LOCKa i dokumentów dziedziczonych,
cytat z korpusu (sek08a), weryfikacja środowiska importem (CPython
3.14.2, numpy 2.4.3, scipy 1.17.1, sympy 1.14.0 — identyczne
z poprzednikiem), weryfikacja mtime npz tła 2π (2026-08-31 21:41 —
zgodny z zapisem poprzednika, READ-ONLY) oraz arytmetyka
transkrypcyjna rejestru wejść (pierwiastki progów; tabela §5).
Analiza wariacyjna §2 wykonana symbolicznie na papierze PRZED kodem —
jest częścią zamrożenia, nie obliczeniem cyklu.

**Rejestr WEJŚĆ (flagowany w każdym zależnym wyniku):**
- **g_ceil = √(4/3) = 1.1547005** [INPUT dokumentacyjny z korpusu,
  M9.1'' canonical — CYTAT sek08a (`rem:hyp-unified-action-M911-canonical`,
  linie ~338–360):
  > „Element objętościowy: √(−g_eff) = c₀ ψ/(4−3ψ) (M9.1'' canonical,
  > `eq:vol-element-M911`, sek08c body; A2 audit closure 2026-05-02)."
  > „Metryka efektywna: M9.1'' hiperboliczna ds² = −c₀²(4−3ψ)/ψ dt²
  > + ψ/(4−3ψ) δᵢⱼ dxⁱdxʲ (`eq:metric-M911-canonical`, sek08c body)."
  Biegun w ψ=4/3; w zmiennej g (ψ=g², bo φ=(Φ/Φ₀)^{1/2}):
  g_ceil=√(4/3)];
- g_floor: PRIMARY 0.5458938 (próg QB-2 0.298); wrażliwość TYLKO dla
  biegu genezowego: {0.4438468, 0.5753260} (progi {0.197, 0.331})
  [INPUT, op-substrate-fluctuation-channel Phase 3, mapowanie
  g_floor=√(Φ_c/Φ_vac) dziedziczone z MD poprzednika §1];
- κ = 100 (kara C², FROZEN, dziedziczone; ta sama wartość dla podłogi
  i dla bariery C-BAR);
- seed = 20260902, amplituda szumu 1e−3 [LOCK §2];
- g₀_μ = φ·0.90548 = 1.4650974 [INPUT, definicja jak MD op-3d-canonical];
- β = γ = 1 [INPUT];
- tło 2π: `../op-3d-canonical-lattice-2026-08-31/Phase2_backgrounds3d.npz`
  (klucze `2pi__A1.0__N{32,48}`; READ-ONLY, mtime weryfikowany po użyciu);
- kotwica siatek radialnych h ∈ {0.025, 0.0125}, R=60 [dziedziczone].

## 1. Funkcjonał (LOCK §2) i forma czynnika metrycznego

**PRIMARY:** E[g] = ∫ w(ψ)·[½K(g)|∇g|² + U_b(g)] dx, ψ = g²,
**w(ψ) = ψ/(4−3ψ)** (wyłącznie forma z eq:vol-element-M911 — forbidden
move 8: zakaz modyfikacji wykładników), K = g⁴,
**U_b(g) = U(g) − U(1) + V_fl(g)**, U = βg⁷/7 − γg⁸/8,
V_fl = (κ/3)(g_floor−g)³ dla g<g_floor (dziedziczona kara C²).

**Decyzja FROZEN (referencja próżniowa):** gęstość potencjalna wchodzi
do nawiasu jako U(g) − U(1) (energia WZGLĘDEM próżni; U(1)=1/56).
Uzasadnienie (pre-compute, analityczne): wariacja (§2) zawiera człon
w′(ψ)·2g·[…]; przy surowym U mielibyśmy δE/δg|_{g≡1} = w′(1)·2·U(1)
= 8/56 = 1/7 ≠ 0 — próżnia analitycznie NIEstacjonarna i P1a (wymóg
LOCKa §3: „próżnia zostaje w g≡1, oba warianty, dryf ≤1e−10") byłoby
niespełnialne z definicji, w OBU sektorach (także P1c). Litera LOCKa
P1a wymusza zatem referencję próżniową; jest ona ponadto konwencją
całej linii (P1b/P1c poprzednika: E_rel względem referencji,
ε(2π) = E[g]−E[1]). Czynnik w(ψ) pozostaje nietknięty.

**C-BAR (kontrolny, nieusuwalny):** w ≡ 1;
U_b = U(g) − U(1) + V_fl(g) + V_ceil(g),
**V_ceil = (κ/3)(g−g_ceil)³ dla g>g_ceil**, κ=100 (C²: V_ceil′ =
κ(g−g_ceil)², V_ceil″ = 2κ(g−g_ceil) — obie znikają w g_ceil;
struktura wariacyjna zachowana, jak podłoga). Stała −U(1) w C-BAR
nie wpływa na dynamikę (w≡1); trzymana dla jednolitości.

## 2. Wyprowadzenie gradient flow PRIMARY (jawna wariacja)

Definiuję W(g) := w(g²) = g²/(4−3g²). Wtedy
E[g] = ∫ [½·𝒦(g)|∇g|² + 𝒰(g)] dx, **𝒦(g) = W(g)·K(g) = g⁶/(4−3g²)**,
**𝒰(g) = W(g)·U_b(g)**.

Wariacja standardowa funkcjonału tej postaci:
**δE/δg = −∇·(𝒦∇g) + ½𝒦′(g)|∇g|² + 𝒰′(g)**, gdzie
- W′(g) = w′(ψ)·2g = [4/(4−3ψ)²]·2g = **8g/(4−3g²)²**
  (w′(ψ) = 4/(4−3ψ)² z ilorazu),
- 𝒦′ = W′K + WK′  (człon w′·2g·½K|∇g|² siedzi w ½𝒦′|∇g|²),
- 𝒰′ = W′U_b + WU_b′  (człon w′·2g·U_b jawnie obecny).

Gradient flow (schemat dziedziczony FROZEN): ∂g/∂t = −δE/δg.
Sprawdzenie stacjonarności próżni: 𝒰′(1) = W′(1)·U_b(1) + W(1)·U_b′(1)
= 8·0 + 1·0 = 0 (U_b(1)=0 dzięki referencji próżniowej; U_b′(1)=U′(1)=0
z β=γ). Tachioniczność zachowana: 𝒰″(1) = U″(1) = −1
(człony W″U_b i 2W′U_b′ znikają w g=1).

**C-BAR:** 𝒦 = K, 𝒰 = U_b — silnik identyczny z dziedziczonym
(Flow{Radial,3D} z podmianą U_tot).

**Dyskretyzacja:** dokładnie dziedziczona struktura strumieniowa
(gradient dyskretnej energii): t_flux = 𝒦(g_mid)Δg/h,
t_quad = ¼𝒦′(g_mid)(Δg/h)², człon lokalny h·𝒰′(g) — podmiana
(K,K′,U_tot′) → (𝒦,𝒦′,𝒰′) zachowuje dokładność gradientu dyskretnego
E_h (𝒦 na midpointach jak K u poprzednika). Radialnie miara r²dr,
Neumann; 3D periodycznie.

**Krok czasowy (dziedziczony FROZEN):** stabilizowany semi-implicit
Euler (I − dt·A_t·L)(g^{n+1}−g^n) = dt·rhs(g^n), A_t = 1.05·𝒦(max g^n)
(adaptowany; w PRIMARY 𝒦 zamiast K — konsekwentnie), L = laplasjan
(3D FFT / radialnie solve_banded). **dt = 0.01, kontrola dt/2 = 0.005;
t_max = 200; stacjonarność ‖rhs‖∞ ≤ 1e−8 próbkowana co Δt=1.**
Niefinityczność ⟹ BREAKDOWN (obsługa błędów dziedziczona z korekty 1
poprzednika; klasyfikacja: załamanie nie-nukleacyjne → INCONCLUSIVE).
E(t) raportowane deskryptywnie (monotonia = sanity).

## 3. Numeryczna obsługa bieguna w (FROZEN; LOCK: TYLKO powyżej g_ceil−1e−6)

w(ψ) ma biegun w ψ=4/3 (to jest fizyczne domknięcie). Regularizacja
WYŁĄCZNIE powyżej **g_c := g_ceil − 1e−6** (litera LOCKa):

**W̃(g) := w(ψ̃), ψ̃ := min(g², ψ_c), ψ_c := g_c²;
W̃′(g) := 8g/(4−3g²)² dla g² < ψ_c, W̃′(g) := 0 dla g² ≥ ψ_c.**

- Dla g < g_c formuły DOKŁADNE (zero regularizacji w domenie metryki).
- Płaska kontynuacja (W̃ stałe, W̃′=0 ponad progiem) czyni flow
  DOKŁADNYM gradient flow dobrze określonej energii zregularyzowanej
  Ẽ (W̃ ciągłe; kink C⁰ w g_c — jak kary); nie fabrykuje sił
  z ekstrapolacji poza domenę metryki (dla ψ>4/3 metryka M9.1''
  nie istnieje: g_tt zmienia znak).
- Wartość na progu: W̃_max = w(ψ_c) ≈ 1.92e5 (policzona w skrypcie
  z formuły) — skończona, duża: „nieskończony koszt objętościowy"
  reprezentowany saturacją.
- Warunek na ψ̃=min(g²,ψ_c) (nie na g) obsługuje też patologiczne
  g<−g_c bez zmiany znaku W̃ (pole nie powinno tam wejść; wejście
  = kandydat BREAKDOWN).

## 4. Detektory nukleacji (FROZEN; dolny dziedziczony + górny z LOCKa)

- **Dolny (dziedziczony dosłownie):** maska m_dn = (g < g_thr_dn),
  g_thr_dn = (1+g_floor)/2; `scipy.ndimage.label` (3D 6-spójność,
  sklejanie periodyczne union-find przez pary ścian; radialnie 1D bez
  periodyczności); próbkowanie co Δt=1; N_seed_dn = liczba obiektów
  w t=0; NUKLEACJA_dn: pierwsze t₀ z N_dn(t) > N_seed_dn dla wszystkich
  t ∈ {t₀,…,t₀+10} (≥10 j.cz.).
- **Górny (zalockowany TERAZ, te same reguły):** maska m_up =
  (g > g_thr_up), **g_thr_up = (1+g_ceil)/2 = 1.0773503**; N_seed_up
  z t=0; NUKLEACJA_up analogicznie.
- Oba detektory biegną RÓWNOLEGLE w każdym biegu; zdarzenie nukleacji =
  pierwszy z detektorów, który potwierdzi okno (raportowany kierunek).
- **Zbieżność (Q-PASS-NUCLEATION):** nukleacja (ten sam kierunek)
  na OBU siatkach pary ORAZ w biegach dt/2 na obu siatkach; liczba
  obiektów w chwili detekcji zgodna ±1 między czterema biegami danego
  (start × wariant × podłoga).
- Detektory NIEZMIENIALNE po pierwszym biegu (forbidden move 9).

## 5. Siatki, starty, tabela progów (FROZEN)

| wielkość | wartość |
|---|---|
| g_ceil | 1.1547005 |
| g_thr_up | 1.0773503 |
| g_floor f1/f2/f3 | 0.4438468 / **0.5458938 (PRIMARY)** / 0.5753260 |
| g_thr_dn f1/f2/f3 | 0.7219234 / 0.7729469 / 0.7876630 |

**Siatki:** radialnie h ∈ {0.025, 0.0125}, R=60; sieć 2π: L=2π,
N ∈ {32,48}; **geneza: L=4π, N ∈ {48,64}** (k_min=0.5 — pasmo
tachionowe częściowo w pudle). Porównania siatkowe: radialnie
interpolacja fine→coarse; 3D wspólna podsiatka 16³ (sieć: stride 2/3;
geneza: stride 3/4). Zalockowanych siatek NIE zmniejszamy.

**Starty:**
(i) geneza: g₀ = 1 + f, szum pasmowy |n_i| ≤ 8 (k = n/2 w L=4π,
    tj. k ≤ 4): rng(20260902) losuje standard_normal((17,17,17,2)),
    hermityzacja C_sym = (C + conj(C[::-1,::-1,::-1]))/2, wbudowanie
    TYCH SAMYCH współczynników w siatki N=48 i N=64 (identyczne pole
    ciągłe), normalizacja max|f| = 1e−3 ze zbudowanego N=64;
(ii) soliton radialny μ: profil ODE kanoniczny (DOP853, rtol=1e−12,
    atol=1e−14, max_step=0.02 — verbatim dziedziczone), g₀_μ=1.4650974,
    dense-output na obu siatkach radialnych; UWAGA pre-rejestrowana:
    rdzeń profilu (g do 1.465) leży POWYŻEJ g_ceil — start częściowo
    poza domeną metryki PRIMARY; bieg wykonujemy wg litery (wynik
    raportowany, przewidywany kandydat BREAKDOWN/INCONCLUSIVE
    w PRIMARY), C-BAR w pełni określony;
(iii) sieć 2π: klucze `2pi__A1.0__N{32,48}` z npz (READ-ONLY, mtime
    weryfikowany po każdym odczycie).

**Macierz P2 (LOCK §3): 14 biegów głównych dt=0.01:**
- geneza × PRIMARY × {f1,f2,f3} × N∈{48,64} — 6;
- soliton × {PRIMARY, C-BAR} × f2 × h∈{0.025,0.0125} — 4;
- sieć 2π × {PRIMARY, C-BAR} × f2 × N∈{32,48} — 4.
Biegi dt/2: warunkowo, gdy para (start×wariant×podłoga) ma nukleację
na obu siatkach (wymóg zbieżności werdyktu).

## 6. Phase 1 — specyfikacja (FROZEN)

- **P1a:** start DOKŁADNIE g≡1, bez zaburzeń, t=10, gate
  max_t‖g−1‖∞ ≤ 1e−10; oba warianty (PRIMARY f2, C-BAR f2);
  geometrie: radialna h=0.0125, 3D L=2π N=32, 3D L=4π N=48.
- **P1b (reprodukcja BREAKDOWN bez domknięć):** w≡1, BEZ V_fl i BEZ
  V_ceil (czyste U tach., silnik dziedziczony), start solitonowy,
  h ∈ {0.025, 0.0125}; PASS ⟺ BREAKDOWN z t_break ∈ [2.7, 3.2]
  na obu siatkach (poprzednik: 2.75 / 3.13) — osiągalny FAIL
  ciągłości.
- **P1c (sektor stabilny m²=+γ, nieusuwalny):** U_stab = −βg⁷/7+γg⁸/8,
  referencja próżniowa analogicznie (U_stab − U_stab(1)); pełny silnik
  z domknięciami i OBOMA detektorami; macierz (FROZEN): sieć 2π N=32
  × {PRIMARY, C-BAR}, geneza N=48 × {PRIMARY, C-BAR}, soliton h=0.025
  × C-BAR (5 biegów; soliton w PRIMARY wyłączony z KONTROLI, bo start
  leży częściowo poza domeną metryki PRIMARY — kontrola wymaga startu
  w domenie ważności modelu; jego zachowanie w sektorze tachionowym
  i tak raportuje P2). PASS ⟺ wszystkie: STATIONARY, ‖g−1‖∞ ≤ 1e−3,
  ZERO alarmów OBU detektorów.
- Dowolny FAIL ⟹ STOP (litera LOCKa).

## 7. Phase 3 (warunkowe przy Q-PASS-STATIC) — Hessian E konsekwentnie

Druga wariacja E wokół stanu zrelaksowanego g₀ (funkcjonał §1–§3,
z W̃): δ²E[φ] = ∫ [𝒦|∇φ|² + 2𝒦′φ∇g₀·∇φ + (½𝒦″|∇g₀|² + 𝒰″)φ²] dx.
Implementacja: DOKŁADNY Hessian dyskretnej energii E_h (druga pochodna
analityczna wyrażenia strumieniowego per ściana + per komórka;
macierz rzadka symetryczna). Problem własny uogólniony
H_h φ = ω² M φ, **M = diag(W̃(g₀)K(g₀))** (waga kinetyczna czasowa
z czynnikiem w — konsekwentnie z akcją ważoną w; C-BAR: M=diag(K(g₀))).
eigsh (najmniejsze), v0 deterministyczne rng(20260902). Mody zerowe
(translacyjne/podłogowe) identyfikowane PRZED interpretacją (overlap
z ∂g₀/∂x_i; lokalizacja w strefach kar). Zbieżność siatkowa
≤ 0.05·max(|ω²_min|,0.1); **Q3-PASS: ω²_min ≥ −1e−3**.
Przy Q-PASS-NUCLEATION zamiast spektrum: charakterystyka kaskady BEZ
progów — N_obj(t) obu detektorów, przyrosty/j.cz., rozkłady rozmiarów
w t₀, t₀+5, t₀+10 (dziedziczona forma deskryptywna).

## 8. Higiena wykonania

Pełne ścieżki bez `cd`; `ls` po każdym zapisie; runy >10 min w tle
(proces ≤ ~55 min; genezowe 64³ z checkpointami npz co 10 j.cz.,
wznawialne `--resume`); wyniki per-bieg json/npz w `Phase2_results/`;
`--verdict` składa całość do `Phase2_output.txt`
+ `Phase2_relaxed_states.npz`. INCOMPLETE raportowane z przyczyną.
Rdzeń `.tex`, STATE.md, git, katalogi innych cykli — NIETYKANE
(npz READ-ONLY, mtime po każdym odczycie). INCONCLUSIVE ≠ pozytyw.

**FROZEN 2026-09-02, przed uruchomieniem jakiegokolwiek skryptu cyklu.**
