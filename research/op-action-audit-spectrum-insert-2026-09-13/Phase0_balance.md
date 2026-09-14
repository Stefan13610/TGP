---
title: "Phase0_balance — LOCK: audyt jednej akcji (pęd kanoniczny, Hamiltonian, dyspersja próżni) + operacyjna definicja i rachunek ΔE_insert we właściwej parze (w, V_M9.1'') — realizacja P0.1+P0.3 audytu zewnętrznego"
date: 2026-09-13
type: phase0-lock
tgp_owner: research/op-action-audit-spectrum-insert-2026-09-13
status: PHASE0-LOCKED
computations_performed: ZERO
authorization: "User 2026-09-13: wybór ścieżki «Audyt analityczny P0.1+P0.3» (spośród opcji: P0.1+P0.3 / dopisek core N4 / LOCK N2 materia / LOCK N3 dynamika 2. rzędu) po analizie audytu zewnętrznego TGP_analiza_i_priorytety.pdf"
anti_lakatos_lock: ACTIVE
related:
  - "[[README.md]]"
  - "[[HANDOFF_PROMPT.md]]"
  - "[[../op-metric-pair-M911-2026-09-02/Phase_FINAL_close.md]]"
  - "[[../op-metric-pair-M911-2026-09-02/NEEDS.md]]"
  - "[[../op-metametric-boundary-2026-09-01/Phase_FINAL_close.md]]"
  - "[[../../core/sek08a_akcja_zunifikowana/sek08a_akcja_zunifikowana.tex]]"
  - "[[../../core/sek08c_metryka_z_substratu/sek08c_metryka_z_substratu.tex]]"
---

# Phase 0 — LOCK cyklu `op-action-audit-spectrum-insert`

**ZERO obliczeń wykonanych przed zapisaniem tego dokumentu.**

---

## 0. Pytanie i kontekst

Audyt zewnętrzny (TGP_analiza_i_priorytety.pdf, snapshot commita
5c543d76 z 2026-09-02) wskazał dwa najtańsze testy rozstrzygające,
których korpus dotąd nie wykonał:

- **P0.1:** „wyprowadzić z jednej zapisanej akcji pęd kanoniczny,
  Hamiltonian, statykę i dyspersję próżni, po podstawieniu zależnej
  od ψ metryki" — poprzednik (op-metric-pair-M911) wykonał krajobraz
  statyczny (Q-A-PASS), ale Phase 4 (widmo) nie zaszła (warunkowa),
  a sektor CZASOWY akcji (współczynnik M(ψ) przy ψ̇²) nigdy nie był
  audytowany we właściwej parze. Audyt (rozdz. 2) raportuje, że
  w STAREJ hybrydzie literalne równanie falowe dawało tachion
  (ω²=k²−1) — pytanie, czy właściwa para to naprawia, jest OTWARTE.
- **P0.3:** Q1-POS poprzednika (op-metametric-boundary) porównywał
  energię tego samego profilu względem RÓŻNYCH teł (mieszane znaki
  {−,+,−}); audyt (rozdz. 5): to nie jest poprawnie zdefiniowany koszt
  wstawienia obiektu. Potrzebna definicja ΔE_insert na WSPÓLNYM tle,
  z kontrolą pudła i siatki, zanim jakikolwiek claim o „zerokosztowej
  kreacji" (hipoteza metametryczna, NEEDS N5 poprzednika).

**Pytania binarne:**
- **Q-D1 (dyspersja z jednej akcji, analitycznie — RACHUNEK CENTRALNY
  części pierwszej):** czy z jednej literalnej akcji
  (eq:S-TGP-unified-M911-canonical, para w, V_M9.1'', K=K_geoψ⁴,
  odczyt PRIMARY B) wynika: M(ψ)>0 na dziedzinie, 𝒰″(ψ*)>0
  oraz ω²(k)>0 dla wszystkich k (próżnia bez tachionu i bez ducha
  w dynamice zachowawczej)?
- **Q-D2 (koszt wstawienia, lekka numeryka radialna):** czy operacyjnie
  zdefiniowany ΔE_insert(A;R,h) (§1) jest zbieżny (h, R) i jaki ma znak
  na rodzinie amplitud A — czy istnieje kanał ΔE_insert≤0 (kreacja
  opłacalna), czy koszt jest wszędzie dodatni (bariera)?

Uwaga interpretacyjna zalockowana: Q-D1-FAIL (tachion/duch we
właściwej parze) też jest wynikiem ważnym — potwierdza problem znaków
audytu na poziomie rdzenia i wymaga decyzji usera o procedurze
wariacji/sygnaturze (user-gate, bez samodzielnych napraw). Q-D2
z kosztem wszędzie dodatnim NIE obala hipotezy metametrycznej —
mówi, że jej nośnik musi leżeć poza czystym sektorem grawitacyjnym
(spójnie z Q-B-FAIL poprzednika).

## 1. Model ZAMKNIĘTY

- **Formy z korpusu (implementator MUSI odtworzyć z cytatem
  w method_decisions, forbidden move a):** w(ψ)=ψ/(4−3ψ)
  [eq:vol-element-M911, sek08c]; V_M9.1''(ψ)=−γψ²(4−3ψ)²/12
  [eq:V-M911]; K(ψ)=K_geo·ψ⁴ [eq:K-coupling-unified]; pełna metryka
  ds²=−c₀²(4−3ψ)/ψ dt² + ψ/(4−3ψ) δᵢⱼdxⁱdxʲ oraz √−g_eff=c₀ψ/(4−3ψ)
  [sek08a ~310–380, sek08c rem:vol-element-M911]. K_geo=γ=c₀=1
  [LOCK, bezwymiarowo].
- **Odczyt kinetyki:** PRIMARY = odczyt B (rozstrzygnięty CYTATEM
  w op-metric-pair-M911, MD §2: w·g_eff^{ij}≡1 ⟹ 𝒦=K_geoψ⁴).
  DZIEDZICZONY — nie rozstrzygać ponownie; odczyt A odnotować
  w wynikach równolegle (tylko tam, gdzie różnicuje).
- **Sektor czasowy (przedmiot audytu, NIE zakładać z góry):**
  M(ψ) wyprowadzić ŚCIŚLE z literalnej akcji: człon
  √−g_eff·½K·|g_eff^{tt}|ψ̇², symbolicznie (sympy), z cytatem formy
  metryki. Wynik M(ψ) jest OUTPUTEM Phase 1, nie wejściem.
- **Funkcjonał statyczny odniesienia (gate spójności):**
  E[ψ]=∫[½ψ⁴|∇ψ|²+𝒰(ψ)]dx, 𝒰=w·V=γ(ψ⁴/4−ψ³/3) — DOKŁADNIE PRIMARY
  poprzednika; redukcja statyczna Hamiltonianu MUSI go odtworzyć
  symbolicznie (tożsamość sympy, nie float).
- **ΔE_insert (definicja operacyjna, FROZEN):**
  ΔE_insert(A;R,h) = E[ψ_A^relax] − E[ψ≡1] na IDENTYCZNYM pudle
  radialnym [0,R], siatce h i warunkach brzegowych (jednorodny
  Neumann w 0 i R), gdzie ψ_A^relax = wynik relaksacji gradient flow
  funkcjonału E z WIĘZEM ψ(0)=A (pin centrum; jedyny więz; bez
  podłóg/barier — dziedzina (0,4/3), pas ψ>4/3−1e−6 = klasyfikacja
  BREAKDOWN-BOUNDARY jak u poprzednika). Tło wspólne = próżnia ψ≡1
  w TYM SAMYM pudle (żadnych porównań między różnymi tłami — sedno
  korekty P0.3). Rodzina amplitud A ZALOCKOWANA:
  {0.50, 0.70, 5/6, 7/6, 1.25, 1.30}; pudła R∈{60,120};
  siatki h∈{0.025,0.0125}; dt=0.01, stacjonarność ‖ψ̇‖∞≤1e−8,
  t_max=200; seed niepotrzebny (starty deterministyczne: gauss
  σ=5.0 wokół 0 z ψ(∞)=1, amplituda dopasowana do A — jak bump
  poprzednika [INPUT-MD]).
- **Rozróżnienie zalockowane:** ΔE_insert to ENERGIA stanu końcowego
  względem próżni; BARIERA kreacji (ścieżka minimalna) jest POZA
  zakresem tego cyklu — zakaz wnioskowania o barierze z ΔE_insert.
- **Ograniczenie zakresu (odnotowane):** w czystym sektorze jedynym
  jednorodnym stanem stacjonarnym jest ψ*=1 (Q-A poprzednika), więc
  „rodzina ośrodków n" audytu redukuje się tu do rodziny amplitud A
  na tle próżni; rozszerzenie na ośrodki n≠1 wymaga sprzężenia
  z materią (NEEDS N2 poprzednika) — poza tym LOCKiem.
- **Rejestr WEJŚĆ:** K_geo=γ=c₀=1; rodzina A; R, h, dt, progi;
  σ_bump=5.0; brak seeda.

## 2. Fazy i kryteria (zalockowane)

### Phase 1 — kanonika z jednej akcji (sympy, zero numeryki)
- P1a: z literalnej akcji: gęstość Lagranżjanu L(ψ,ψ̇,∇ψ)
  = ½M(ψ)ψ̇² − ½𝒦(ψ)|∇ψ|² − 𝒰(ψ) z jawnymi M,𝒦,𝒰 (symbolicznie,
  z cytatami form); pęd π=∂L/∂ψ̇; Hamiltonian H[π,ψ]; równanie ruchu
  Eulera–Lagrange'a z członami M′,𝒦′.
- P1b (gate spójności statyki): δH/δψ przy π=0 ≡ δE/δψ funkcjonału
  PRIMARY poprzednika — tożsamość sympy (simplify=0). FAIL ⟹ STOP
  i raport (oznaczałoby, że akcja i funkcjonał relaksacji poprzednika
  NIE pochodzą z jednego zapisu — dokładnie rozjazd, przed którym
  ostrzega audyt rozdz. 2).
- P1c (gate implementacji): wartości M,𝒦,𝒰 w {0.5, 1, 7/6, 1.3}
  sympy vs float — zgodność 1e−12 (osiągalny FAIL).
- Jeżeli literalna akcja WYMAGA traktowania g_eff jako zmiennej
  niezależnej albo wariacji euklidesowej (kaweat audytu) — STOP,
  raport wprost, user-gate (zakaz samodzielnego wyboru procedury).

### Phase 2 — Q-D1: dyspersja próżni (sympy + kontrola float)
- Linearyzacja równania ruchu wokół ψ*=1 (ψ=1+ε·e^{i(kx−ωt)}):
  ω²(k) = [𝒦(1)k² + 𝒰″(1)] / M(1) — wyprowadzić, nie postulować;
  raportować m²=𝒰″(1)/M(1), c_s²=𝒦(1)/M(1) oraz znaki M(ψ), 𝒦(ψ)
  na całej dziedzinie (0,4/3).
- **Q-D1-PASS:** M(1)>0 ∧ 𝒰″(1)>0 ∧ ω²(k)>0 ∀k≥0 ORAZ M(ψ)>0,
  𝒦(ψ)>0 na (0,4/3). **Q-D1-FAIL:** przeciwnie (tachion: ω²(k)<0
  dla pewnych k; duch: M≤0 gdziekolwiek na dziedzinie) — raport
  wprost z miejscem załamania znaku.
- Odczyt A (𝒦_A=ψ⁵/(4−3ψ)): dyspersję policzyć równolegle,
  odnotować, NIE wpływa na werdykt (PRIMARY=B).

### Phase 3 — Q-D2: ΔE_insert (lekka numeryka radialna)
- P3a (bramka maszynerii): próżnia ψ≡1 bez więzu zostaje (dryf
  ≤1e−10, t=10, obie siatki); relaksacja z pinem A=1 daje
  ΔE_insert=0±1e−10 (kontrola zerowa). FAIL ⟹ STOP.
- P3b: pełna macierz A×R×h (6×2×2=24 relaksacje radialne, tanie);
  werdykt per A: zbieżność h (|Δ|/max(|ΔE|,1e−6)≤5e−3) i R
  (|ΔE(120)−ΔE(60)| raportowane; brak członu objętościowego =
  warunek poprawności definicji).
- **Werdykty (litera):** **Q-D2-COST** (ΔE_insert(A)>0 dla wszystkich
  A≠1, zbieżnie) / **Q-D2-CHANNEL** (istnieje A z ΔE_insert≤0,
  zbieżnie — kanał opłacalny; raportować profil ψ_A^relax:
  rozciągłość, ψ_max/ψ_min) / **Q-D2-INCONCLUSIVE** (brak zbieżności
  h lub jawna zależność objętościowa od R — wtedy definicja wymaga
  poprawki, raport bez werdyktu znaku; BREAKDOWN-BOUNDARY = osobna
  kategoria deskryptywna).
- Deskryptywnie obowiązkowo: ΔE_insert(A) jako tabela + monotonia;
  odniesienie do liczb Q1-POS poprzednika (−0.179 / +16156.6) —
  wyjaśnienie różnicy definicyjnej (wspólne tło vs różne tła).

## 3. Forbidden moves
Dziedziczone z op-metric-pair-M911 (§3) + (a) formy wyłącznie
z sek08a/sek08c z cytatem, zakaz modyfikacji; (b) M(ψ) jest wynikiem
wyprowadzenia, nie założeniem; (c) zakaz podłóg/barier; (d) rodzina A,
progi i więz pinu niezmienialne po pierwszym biegu; (e) rdzeń
.tex/STATE.md/git NIETYKANE (commit/STATE robi sesja główna);
(f) katalogi innych cykli tylko odczyt; (g) INCONCLUSIVE ≠ pozytyw;
(h) zakaz wnioskowania o barierze kreacji z ΔE_insert.

## 4. Deliverables
`Phase_method_decisions.md` (FROZEN, cytaty form + jawna wariacja),
`Phase1_canonical.py`+output, `Phase2_dispersion.py`+output,
`Phase3_insert_cost.py`+output+json/npz per bieg,
`Phase_FINAL_close.md`, `NEEDS.md`, `README.md`.

## 5. Drzewo decyzyjne
```text
P1b FAIL (akcja ≠ funkcjonał statyczny) → STOP; user-gate: rozjazd
   rdzeniowy akcja↔relaksacja (dokładnie ostrzeżenie audytu rozdz. 2)
Q-D1-PASS → operator kinetyczny (M,𝒦) ZALOCKOWANY dla przyszłego cyklu
   dynamiki 2. rzędu (NEEDS N3 poprzednika staje się dobrze postawiony);
   m², c_s² = punkt odniesienia dla P1.3 (nośnik propagacji) i P2 (sondy)
Q-D1-FAIL → user-gate CORE: problem znaków audytu potwierdzony we
   właściwej parze (procedura wariacji / sygnatura / konwencja) —
   ZERO samowolnych napraw
Q-D2-COST → „zerokosztowa kreacja" bez nośnika w czystym sektorze
   (spójne z Q-B-FAIL poprzednika) → wzmacnia gałąź N2 (materia)
Q-D2-CHANNEL → istnieje kanał opłacalny → kandydat startów dla cyklu
   dynamiki 2. rzędu (z Q-D1) LUB dla genezy poziomu 0
Q-D2-INCONCLUSIVE → NEEDS metodologiczny: poprawka definicji
   (więz, brzeg, skala pudła) — bez claimów o znaku
```

---

**LOCK ZAMKNIĘTY 2026-09-13. Zmiany poniżej tej linii po starcie
obliczeń = forbidden move.**
