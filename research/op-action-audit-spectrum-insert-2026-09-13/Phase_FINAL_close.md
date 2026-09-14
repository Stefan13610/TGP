---
title: "Phase_FINAL_close — zamknięcie: Q-D1-PASS (z jednej akcji: M(ψ)=ψ⁶/(4−3ψ)²>0 WYPROWADZONE z członu √−g·½K·|g^tt|ψ̇², 𝒦=ψ⁴>0 na (0,4/3), 𝒰″(1)=1>0, ω²(k)=k²+1>0 ∀k — próżnia bez tachionu i ducha we właściwej parze; m²=c_s²=1; kaweat: literalna kontrakcja ZNAKOWANYM g^tt daje ω²=k²−1 tachionu audytu + napięcie znaku statyka↔R3 ODE, user-gate) + Q-D2-INCONCLUSIVE (ΔE_insert(A;R,h): P3a 8/8 PASS, znak DODATNI wszędzie przy skończonym h i |ΔE(120)−ΔE(60)|=0 dokładnie — zero członu objętościowego, ALE brak zbieżności h: ΔE∝h (punktowy pin ma zerową pojemność w 3D; dokładne minima Newtona: stosunek 0.4988) — definicja wymaga więzu skończonej skali; A=1.30 → BREAKDOWN-BOUNDARY-LOWER zbieżnie)"
date: 2026-09-13
type: phase-final-close
tgp_owner: research/op-action-audit-spectrum-insert-2026-09-13
status: CLOSED
verdict: "Q-D1: PASS wg litery (sympy: z literalnej akcji eq:S-TGP-unified-M911-canonical z metryką eq:metric-M911-canonical podstawioną jako funkcja ψ: L=½Mψ̇²−½𝒦|∇ψ|²−𝒰 z M(ψ)=K_geoψ⁶/(c₀(4−3ψ)²) WYNIKIEM członu √−g·½K·|g^tt|ψ̇² [litera LOCKa §1], 𝒦=K_geoψ⁴ [odczyt B DZIEDZICZONY], 𝒰=wV=γ(ψ⁴/4−ψ³/3); π=Mψ̇; H=∫[π²/2M+½𝒦|∇ψ|²+𝒰]; P1b tożsamość δH/δψ|π=0≡δE_PRIMARY/δψ simplify=0 PASS; P1c 12/12 ≤5.7e−14 po correction note 1 (harness+fma, progi/punkty nietknięte); ω²(k)=k²+1 wyprowadzone linearyzacją, m²=1, c_s²=1, M>0 i 𝒦>0 na (0,4/3), odczyt A w ψ*=1 identyczny (𝒦_A(1)=1). Deskryptywnie: kontrakcja znakowanym g^tt daje ω²=k²−1 (dokładnie tachion audytu rozdz. 2) i to ONA odtwarza R3 ODE (1−ψ)/ψ², podczas gdy statyka formy kanonicznej daje (ψ−1)/ψ² — napięcie znaku rdzenia raportowane wprost, user-gate w NEEDS. Q-D2: INCONCLUSIVE wg litery (brak zbieżności h): P3a 8/8 PASS (dryf 0.0, ΔE(pin A=1)=0.0 dokładnie); P3b 24/24 wykonane: ΔE_insert>0 dla wszystkich A≠1 przy każdym skończonym h, |ΔE(R=120)−ΔE(R=60)|=0 dokładnie (zero członu objętościowego — sedno korekty P0.3 działa), ALE conv_h∈{0.94–1.00}≫5e−3 dla A∈{0.5,0.7,5/6,7/6,1.25}: ΔE∝h→0 (zerowa pojemność punktowego więzu w 3D; potwierdzone dokładnymi minimami warunkowymi Newtona: ΔE_exact(h/2)/ΔE_exact(h)=0.4988, diag deskryptywna); A=1.30: BREAKDOWN-BOUNDARY-LOWER zbieżnie na obu siatkach (pin przy 1.3 wypycha węzeł sąsiedni poniżej dolnego krańca dziedziny t≈2.5–2.7, kategoria deskryptywna). ZAKAZ wnioskowania o barierze kreacji z ΔE_insert dotrzymany; INCONCLUSIVE ≠ pozytyw."
anti_lakatos_lock: PRESERVED
tags: [action-audit, canonical-momentum, hamiltonian, dispersion, no-tachyon-canonical, sign-tension-R3, insert-cost, zero-capacity-pin, inconclusive-qd2, qd1-pass, closed]
related:
  - "[[Phase0_balance.md]]"
  - "[[Phase_method_decisions.md]]"
  - "[[Phase_correction_note_P1c_gate.md]]"
  - "[[README.md]]"
  - "[[NEEDS.md]]"
  - "[[../op-metric-pair-M911-2026-09-02/Phase_FINAL_close.md]]"
  - "[[../op-metametric-boundary-2026-09-01/Phase_FINAL_close.md]]"
  - "[[../../core/sek08a_akcja_zunifikowana/sek08a_akcja_zunifikowana.tex]]"
  - "[[../../core/sek08c_metryka_z_substratu/sek08c_metryka_z_substratu.tex]]"
---

# Phase FINAL — zamknięcie cyklu op-action-audit-spectrum-insert

**Status: CLOSED-EXECUTED (2026-09-13, jedna sesja: LOCK → method_decisions
→ Phase 1 → Phase 2 → Phase 3 → zamknięcie).** Kryteria LOCKa
(`Phase0_balance.md` §1–2, §5) stosowane DOSŁOWNIE; zero zmian
kryteriów/progów/form/rodziny A/więzu po starcie; jedna korekta
implementacyjna gate'u (correction note 1 — §5).

---

## 0. Werdykty

| Pytanie | Werdykt | Jedno zdanie |
|---|---|---|
| **Q-D1** (dyspersja z jednej akcji — RACHUNEK CENTRALNY cz. 1) | **Q-D1-PASS** (litera §2 Phase 2) | M(1)=1>0, 𝒰″(1)=1>0, ω²(k)=k²+1>0 ∀k≥0, M>0 i 𝒦>0 na całej (0,4/3) — próżnia właściwej pary bez tachionu i bez ducha (w zamrożonym odczycie \|g^tt\| członu czasowego) |
| P1b (gate akcja↔funkcjonał statyczny) | **PASS** | δH/δψ\|_{π=0} ≡ δE_PRIMARY/δψ — tożsamość sympy simplify=0 (H przy π=0 JEST funkcjonałem PRIMARY poprzednika) |
| P1c (gate implementacji) | **PASS 12/12** | ≤5.7e−14 (próg 1e−12) po correction note 1 (harness porównuje w identycznym wejściu binarnym; referencja float M z fma; progi/punkty NIETKNIĘTE) |
| P3a (bramka maszynerii) | **PASS 8/8** | dryf próżni 0.0 ≤1e−10 (4 pudła), pin A=1 ⟹ ΔE=0.0 dokładnie ±1e−10 (4 pudła) |
| **Q-D2** (ΔE_insert na wspólnym tle) | **Q-D2-INCONCLUSIVE** (litera §2 Phase 3: „brak zbieżności h") | znak DODATNI wszędzie przy skończonym h i zero zależności od R, ALE ΔE∝h→0 (punktowy pin ma zerową pojemność w 3D) — conv_h≈1≫5e−3; definicja wymaga więzu skończonej skali; BEZ werdyktu znaku |
| A=1.30 | **BREAKDOWN-BOUNDARY-LOWER** (kategoria deskryptywna, zbieżna na obu siatkach) | pin ψ(0)=1.30 wypycha węzeł sąsiedni poniżej dolnego krańca dziedziny (ψ₁→<0) w t≈2.5–2.7 |

**Uwagi interpretacyjne zalockowane (LOCK §0, stosowane dosłownie):**
Q-D2-INCONCLUSIVE ≠ pozytyw; ZAKAZ wnioskowania o barierze kreacji
z ΔE_insert; Q-D2 nie obala hipotezy metametrycznej — mówi o definicji.

## 1. Wejścia i formy (rejestr; MD — rejestr WEJŚĆ, §1–3)

- **Formy z korpusu (CYTATY w MD §1, forbidden move a dotrzymany):**
  w(ψ)=ψ/(4−3ψ) [eq:vol-element-M911, sek08c ~510–519];
  V_M9.1''(ψ)=−γψ²(4−3ψ)²/12 [eq:V-M911, sek08a ~977–981];
  K(ψ)=K_geoψ⁴ [eq:K-coupling-unified ~169–172]; metryka
  ds²=−c₀²(4−3ψ)/ψ dt²+ψ/(4−3ψ)δᵢⱼdxⁱdxʲ [eq:metric-M911-canonical,
  sek08c ~420–424] ⟹ g^tt=−ψ/(c₀²(4−3ψ)), g^ij=(4−3ψ)/ψ δ^ij;
  √−g=c₀ψ/(4−3ψ). K_geo=γ=c₀=1 [LOCK].
- **Odczyt kinetyki przestrzennej:** PRIMARY = odczyt B DZIEDZICZONY
  (MD §2 poprzednika op-metric-pair-M911, nie rozstrzygany ponownie):
  w·g^ij≡1 ⟹ 𝒦=K_geoψ⁴; odczyt A (𝒦_A=ψ⁵/(4−3ψ)) odnotowany
  równolegle.
- **Człon czasowy (litera LOCKa §1):** √−g·½K·|g^tt|·ψ̇² ⟹ M jest
  WYNIKIEM Phase 1. Kaweat znaku odnotowany w MD §3 PRZED
  obliczeniami i policzony równolegle (§2 niżej).
- Rodzina A={0.50,0.70,5/6,7/6,1.25,1.30}; R∈{60,120};
  h∈{0.025,0.0125}; dt=0.01, t_max=200, stacjonarność 1e−8;
  start gauss σ=5.0 [INPUT-MD], ψ₀(0)=A dokładnie; pas 4/3−1e−6;
  progi P3a 1e−10, conv_h 5e−3; brak seeda (deterministycznie);
  zero podłóg/barier.

## 2. Phase 1–2 — kanonika i dyspersja (Phase1_output.txt, Phase2_output.txt)

**Z jednej literalnej akcji (sympy, wszystkie tożsamości simplify=0):**

- **M(ψ) = K_geo ψ⁶/(c₀(4−3ψ)²)** [WYNIK członu √−g·½K·|g^tt|ψ̇²;
  M(1)=1; M>0 na (0,4/3); M′=12ψ⁵(2−ψ)/(4−3ψ)³];
- **𝒦(ψ) = K_geo ψ⁴** [odczyt B; 𝒦(1)=1; 𝒦>0 na (0,4/3); 𝒦′=4ψ³];
- **𝒰(ψ) = w·V = γ(ψ⁴/4−ψ³/3)** [tożsamość wielomianowa sympy;
  𝒰′=γψ²(ψ−1); 𝒰″=γψ(3ψ−2), 𝒰″(1)=γ=1>0];
- **π = M(ψ)ψ̇**;
- **H[π,ψ] = ∫[π²/(2M) + ½𝒦|∇ψ|² + 𝒰]d³x** — przy π=0 dokładnie
  funkcjonał PRIMARY poprzednika (gate P1b: tożsamość simplify=0);
- **EOM:** Mψ̈ + ½M′ψ̇² = ∇·(𝒦∇ψ) − ½𝒦′|∇ψ|² − 𝒰′ (= −δE_PRIMARY/δψ);
- **Dyspersja (wyprowadzona linearyzacją ψ=1+εe^{i(kx−ωt)}, nie
  postulowana):** ω²(k) = [𝒦(1)k²+𝒰″(1)]/M(1) = **k² + 1**;
  **m² = 1, c_s² = 1**. Zgodność z formułą LOCKa simplify=0.
- **Odczyt A równolegle:** 𝒦_A=ψ⁵/(4−3ψ)>0 na (0,4/3); 𝒦_A(1)=1 ⟹
  dyspersja w próżni IDENTYCZNA (ω²=k²+1, m²=c_s²=1); różnicuje
  poza próżnią: c_s²(ψ): B=(4−3ψ)²/ψ², A=(4−3ψ)/ψ (obie dodatnie).
- **Q-D1-PASS** — wszystkie warunki litery spełnione.

**Deskryptywnie (MD §3 kaweat; NIE werdyktotwórcze; user-gate):**
literalna kontrakcja +½Kg^{μν}∂ψ∂ψ ze ZNAKOWANYM g^tt<0 (sygnatura
(−,+,+,+) jak zapisana w eq:metric-M911-canonical) daje
L_sgn = −½Mψ̇²+½𝒦|∇ψ|²−𝒰 (duch w członie czasowym) i **ω² = k²−1 —
dokładnie tachion audytu (rozdz. 2)**. Napięcie znaku rdzenia
(zweryfikowane sympy, obie redukcje simplify=0): statyka formy
kanonicznej |g^tt| ⟹ ψ''+2/rψ'+2ψ'²/ψ = **(ψ−1)/ψ²** (zaniki
Yukawy e^{−r}/r wokół próżni), a R3 ODE korpusu [eq:R3-ODE, dowód
prop:V-M911-canonical] ma **(1−ψ)/ψ²** (ogony oscylacyjne sin(r)/r) —
tę drugą odtwarza wyłącznie wariacja znakowana, która w sektorze
czasowym daje ducha/tachion. Jedna rzeczywista wariacja Lorentzowska
NIE daje jednocześnie zdrowej dyspersji k²+1 i R3 ODE w zapisanej
postaci — dokładnie rozjazd, przed którym ostrzega audyt rozdz. 2.
User-gate w [[NEEDS.md]] (N2); ZERO samowolnych napraw rdzenia.

## 3. Phase 3 — Q-D2: ΔE_insert (Phase3_output.txt, Phase3_results/)

**P3a: PASS 8/8** — próżnia bez więzu: dryf 0.0 (t=10, 4 pudła);
pin A=1: ΔE=0.0 dokładnie (4 pudła; start A=1 jest dokładnie ψ≡1,
𝒰′(1)=0 w arytmetyce float — zera dokładne).

**P3b: tabela ΔE_insert(A;R,h) = E[ψ_A^relax]−E[ψ≡1] (wspólne
tło/pudło/siatka/brzeg; 24/24 biegi):**

| A | ΔE(h=0.025) | ΔE(h=0.0125) | conv_h | ≤5e−3 | \|ΔE(R120)−ΔE(R60)\| | status |
|---|---|---|---|---|---|---|
| 0.50 | +2.6057e−3 | +1.3002e−3 | 1.00e+0 | NIE | 0.0 (dokładnie) | TMAX (plateau E) |
| 0.70 | +1.5810e−3 | +7.8955e−4 | 1.00e+0 | NIE | 0.0 | TMAX (plateau E) |
| 5/6 | +6.8937e−4 | +3.4471e−4 | 1.00e+0 | NIE | 0.0 | TMAX (plateau E) |
| 7/6 | +1.8963e−3 | +9.5637e−4 | 9.83e−1 | NIE | 0.0 | TMAX (plateau E) |
| 1.25 | +6.9203e−3 | +3.5666e−3 | 9.40e−1 | NIE | 0.0 | TMAX (plateau E) |
| 1.30 | +1.0417e−1 | +1.1191e−1 | 6.92e−2 | NIE | 0.0 | BREAKDOWN-BOUNDARY-LOWER (t≈2.5–2.7) |

(wartości identyczne dla R=60 i R=120 co do wszystkich cyfr —
profil zlokalizowany, ogon e^{−r}: ψ(r=10)−1 ≈ 2.5e−9)

- **Zależność od R: ZERO dokładnie** — brak członu objętościowego;
  operacyjna definicja na wspólnym tle (korekta P0.3) działa w tym
  aspekcie poprawnie (kontrast z +16156.6 poprzednika, gdzie koszt
  dominował człon 𝒰(1)·V(R) między RÓŻNYMI tłami).
- **Zbieżność h: BRAK dla wszystkich A** — ΔE maleje ~2× przy h→h/2.
  Przyczyna strukturalna (diagnostyka deskryptywna
  `Phase3_diag_output.txt`): **punktowy więz ψ(0)=A ma zerową
  pojemność w 3D** — minimizer warunkowy zapada się do rdzenia skali
  siatki (r_half=0, jeden węzeł; ψ(h)≈0.96–0.98 już przy A=0.5–0.7),
  a ΔE ∝ h → 0⁺. Dokładne dyskretne minima warunkowe (Newton,
  ‖dH‖≤1e−11, niezależnie od flow): stosunek ΔE(h/2)/ΔE(h)=0.4988
  (A=0.70 i A=1.25) — czyste skalowanie ~h.
- **Znak (deskryptywnie, bez werdyktu — litera INCONCLUSIVE):**
  przy każdym skończonym h ΔE_insert > 0 dla wszystkich A≠1 (spójne
  z Q-A-PASS: próżnia jest jedynym minimum, każde odkształcenie
  kosztuje); monotonia: rośnie z |A−1| po obu stronach, asymetrycznie
  (górna strona droższa: A=7/6 vs 5/6: 9.6e−4 vs 3.4e−4).
- **A=1.30 (kategoria deskryptywna, zbieżna na obu siatkach i obu R):**
  BREAKDOWN-BOUNDARY-LOWER — pin przy 1.30 z członem ½𝒦′|∇ψ|²
  (𝒦′=4ψ³>0) na kolanie skali siatki wypycha węzeł sąsiedni W DÓŁ
  poza dolny kraniec dziedziny (ψ₁: 1.3 → <0 w t≈2.5–2.7); NIE jest
  to werdykt pozytywny ani wejście w pas górny (pas 4/3−1e−6
  nieaktywowany w całym cyklu).
- **Odniesienie definicyjne do Q1-POS poprzednika (−0.179/+16156.6;
  NIE reprodukcja):** tamte liczby porównywały ten sam profil względem
  RÓŻNYCH teł (próżnia vs stan pusty; mieszane znaki {−,+,−})
  w sektorze kanonicznym tachionowym (U′=g⁶(1−g), próżnia na
  maksimum); tu różnica energii DWÓCH stanów na TYM SAMYM tle
  w parze właściwej (próżnia = minimum, Q-A-PASS). Ujemne wartości
  tamtej konstrukcji nie są sprzeczne z dodatnim znakiem tutaj —
  to inna wielkość, inny sektor; człon objętościowy tu kasuje się
  z definicji.

## 4. Mapowanie na drzewo decyzyjne LOCKa §5

- **Q-D1-PASS** → „operator kinetyczny (M,𝒦) ZALOCKOWANY dla
  przyszłego cyklu dynamiki 2. rzędu (NEEDS N3 poprzednika staje się
  dobrze postawiony); m², c_s² = punkt odniesienia dla P1.3 (nośnik
  propagacji) i P2 (sondy)" — dosłownie; M=ψ⁶/(4−3ψ)², 𝒦=ψ⁴,
  𝒰=γ(ψ⁴/4−ψ³/3), m²=c_s²=1. Zastrzeżenie deskryptywne: kaweat
  znaku (§2) czyni wybór konwencji |g^tt| jawnym warunkiem tego
  locka — stąd równoległy user-gate N2.
- **Q-D2-INCONCLUSIVE** → „NEEDS metodologiczny: poprawka definicji
  (więz, brzeg, skala pudła) — bez claimów o znaku" — dosłownie;
  konkretnie: więz skończonej skali (N1 w NEEDS). BREAKDOWN-BOUNDARY
  = osobna kategoria deskryptywna (A=1.30, dolny kraniec).

## 5. Korekty / incydenty / higiena (anti-Lakatos)

- ✓ LOCK przeczytany W CAŁOŚCI przed wszystkim; method_decisions
  FROZEN (cytaty form, jawny człon czasowy i wariacja, schemat pinu,
  pas graniczny) zamknięte PRZED jakimkolwiek kodem cyklu; zero zmian
  kryteriów/progów/form/rodziny A/więzu po starcie.
- ✓ **Korekta 1** (`Phase_correction_note_P1c_gate.md`, zapisana PRZED
  użyciem wyniku; pierwotny output zachowany:
  `Phase1_output_pre_correction1.txt`): harness P1c porównywał sympy
  w dokładnym 13/10 z float w double(1.3) (sam wkład reprezentacji
  wejścia 1.36e−12 > progu przy M′(1.3)≈3.1e4) + naiwna forma float
  M traciła 2.1e−12 przez zaokrąglenie 3ψ wzmocnione 1/(4−3ψ)²;
  korekta: identyczne wejście binarne + math.fma; progi/punkty/formy
  NIETKNIĘTE; po korekcie 12/12 ≤5.7e−14.
- ✓ Kosmetyka Phase 2: literalny `%s` (brak argumentu formatu)
  w jednej linii deskryptywnej odczytu A — poprawiony i przebieg
  powtórzony przed użyciem; zero wpływu na wartości (deterministyczne).
- ✓ **Incydent (zero wpływu na wyniki):** status TMAX zamiast
  STATIONARY w biegach P3b — zdiagnozowany jako strukturalna własność
  zamrożonego schematu „krok→projekcja" (tożsamość punktu stałego:
  rhs_i=0 dla i≥2, rhs₁=−d·A_t·a₀/m₁; zmierzone 76.6 = przewidziane);
  energia ma plateau do 1e−12 przez ostatnie ≥5 j.cz.; diagnostyka
  deskryptywna (`Phase3_diag_constrained_min.py` →
  `Phase3_diag_output.txt`): dokładne minima warunkowe Newtona
  potwierdzają znak i skalowanie ~h; przesunięcie wartości punktu
  stałego flow względem dokładnego minimum: +4.4% (A=0.70) do ~2×
  (A=1.25) — raportowane; schemat NIE był zmieniany po starcie
  (żaden bieg nie był reinterpretowany jako STATIONARY); werdykt
  INCONCLUSIVE nie zależy od tego przesunięcia (skalowanie ~h
  wspólne dla flow i minimów dokładnych).
- ✓ Zakaz podłóg/barier dotrzymany (jedyna obsługa granic: pas
  klasyfikacyjny; BREAKDOWN-BOUNDARY-LOWER klasyfikowany, nie
  korygowany); rodzina A/więz/progi niezmienione po pierwszym biegu;
  M(ψ) wyprowadzone, nie założone; INCONCLUSIVE nie reinterpretowane;
  zakaz wnioskowania o barierze kreacji dotrzymany.
- ✓ Rdzeń `.tex` NIETKNIĘTY; STATE.md nieedytowane; git nieużywany;
  katalogi innych cykli tylko odczyt; pełne ścieżki bez `cd`;
  `ls` po każdym zapisie (zero artefaktów zagnieżdżonych ścieżek);
  sandbox bez /dev/null i heredoc (wszystko przez pliki skryptów).
- Środowisko: CPython 3.14.2, numpy 2.4.3, scipy 1.17.1, sympy 1.14.0
  (identyczne z oboma poprzednikami).

## 6. Odczyt (deskryptywnie, bez claimów poza klasą zbadaną)

1. **P0.1 audytu wykonany:** z jednej zapisanej akcji wyprowadzono
   pęd kanoniczny, Hamiltonian, statykę i dyspersję próżni.
   W zamrożonym odczycie |g^tt| właściwa para (w, V_M9.1'', K=ψ⁴)
   NAPRAWIA tachion starej hybrydy: ω²=k²+1 (m²=+1, bez ducha,
   M>0 wszędzie na dziedzinie). Gate P1b domyka spójność akcja↔
   funkcjonał relaksacyjny poprzednika (H|_{π=0} = E_PRIMARY
   dokładnie).
2. **Ale audytowy problem znaków NIE znika — przenosi się na poziom
   konwencji wariacji:** literalna kontrakcja znakowana odtwarza
   R3 ODE (fundament spektrum mas), lecz kosztem tachionu/ducha;
   odczyt |g^tt| daje zdrową dynamikę, lecz jego statyka to
   (ψ−1)/ψ², nie R3. Rozstrzygnięcie (konwencja/sygnatura/procedura
   wariacji w rdzeniu) jest user-gated (N2) — dokładnie kaweat
   audytu rozdz. 2, teraz z precyzyjną lokalizacją.
3. **P0.3 audytu wykonany połowicznie z zyskiem metodologicznym:**
   definicja ΔE_insert na wspólnym tle eliminuje człon objętościowy
   DOKŁADNIE (główna wada Q1-POS usunięta), ale punktowy więz ψ(0)=A
   ma zerową pojemność w 3D — ΔE∝h→0, więc wielkość w tej postaci
   nie ma granicy kontinuum innej niż 0⁺. Deskryptywnie: przy każdym
   skończonym h koszt jest dodatni i rosnący z |A−1| (spójne
   z samodomkniętym krajobrazem Q-A-PASS); żadnego kanału ΔE≤0 nie
   zaobserwowano w klasie zbadanej. Poprawka definicji (więz
   skończonej skali) — user-gated (N1).

## 7. Pliki cyklu

`Phase0_balance.md` (LOCK) · `HANDOFF_PROMPT.md` ·
`Phase_method_decisions.md` (FROZEN) ·
`Phase_correction_note_P1c_gate.md` ·
`Phase1_canonical.py` → `Phase1_output.txt`
(+ `Phase1_output_pre_correction1.txt` — pierwotny, zachowany) ·
`Phase2_dispersion.py` → `Phase2_output.txt` ·
`Phase3_insert_cost.py` → `Phase3_output.txt` + `Phase3_results/`
(json+npz per bieg, `gates_P3a.json`, `verdict.json`)
+ `Phase3_results_batch.log` ·
`Phase3_diag_constrained_min.py` → `Phase3_diag_output.txt`
(diagnostyka deskryptywna) · `NEEDS.md` (user-gated) ·
`README.md` (log).
