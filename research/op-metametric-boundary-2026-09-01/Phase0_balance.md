---
title: "Phase0_balance — LOCK: granica metametryczna jako właściwy punkt odniesienia stabilności — równanie stanu kreacji, podłoga substratowa z QB-2, stabilność w stanie zrelaksowanym; załamanie przez nukleację = pre-rejestrowany POZYTYW"
date: 2026-09-01
type: phase0-lock
tgp_owner: research/op-metametric-boundary-2026-09-01
status: PHASE0-LOCKED
computations_performed: ZERO
authorization: "User 2026-09-01: hipoteza wskazana w README op-blocked-soliton-bang (»koszty tworzenia nowych solitonów są zerowe w takim stanie«; »wielki wybuch trwa na granicy metryki«); dyspozycja: »zakłada, że to powinna być granica układu; obliczenia mogą być trudne, ale jeżeli załamią się ze względu na generowanie obiektów, to w sumie będzie wynik pozytywny. I tak zapisz nowy cykl«."
anti_lakatos_lock: ACTIVE
related:
  - "[[README.md]]"
  - "[[../op-blocked-soliton-bang-2026-07-04/README.md]]"
  - "[[../op-blocked-soliton-bang-2026-07-04/POST_CLOSE_BRAINSTORM_UPDATE.md]]"
  - "[[../op-substrate-fluctuation-channel-2026-08-23/Phase_FINAL_close.md]]"
  - "[[../op-3d-canonical-lattice-2026-08-31/Phase_FINAL_close.md]]"
  - "[[../op-3d-canonical-lattice-2026-08-31/Phase2_backgrounds3d.npz]]"
---

# Phase 0 — LOCK cyklu `op-metametric-boundary`

**ZERO obliczeń wykonanych przed zapisaniem tego dokumentu. Kryteria,
progi, kontrole, detektor nukleacji i forbidden moves zamknięte PRZED
kodem.**

---

## 0. Hipoteza dziedziczona i re-rama pytania o stabilność

**Hipoteza autora (op-blocked-soliton-bang, „Kluczowa luka"):** właściwą
granicą układu jest stan METAMETRYCZNY — stan, w którym koszt kreacji
nowych solitonów/obiektów ΔE_create → 0 (kreacja darmowa). Jeśli taka
granica istnieje, układ relaksuje DO NIEJ, i to ona (nie próżnia g=1,
nie arbitralna sieć d) jest właściwym punktem odniesienia stabilności.
Wszystkie dotychczasowe werdykty Q-FAIL (bloch-chain, symplectic,
3d-canonical) mierzyły stabilność w zespole KANONICZNYM (ustalona liczba
obiektów) względem konfiguracji, które NIE są stanem zrelaksowanym —
mogły więc mierzyć nie „śmierć solitonu", lecz „relaksację ku granicy,
której nikt nie wyznaczył".

**Fakt konwencyjny (do formalizacji w Phase 1):** U(g)=g⁷/7−g⁸/8 ma
U(1)=+1/56 > U(0)=0 i U″(1)=−1 — próżnia na maksimum potencjału; stan
pusty g=0 energetycznie niżej. Obserwowane ucieczki g→0 we wszystkich
cyklach = relaksacja w tę stronę. W kontinuum granicy nie widać;
kandydat na granicę: poziom 0 — **spinodala rozrzedzenia
Φ_c/Φ_vac = 0.298 (skan 0.197–0.331)** z op-substrate-fluctuation-channel
QB-2.

## 1. Trzy pytania binarne

- **Q1 (kontinuum, równanie stanu kreacji):** czy w kontinuum
  kanonicznym istnieje skończona granica metametryczna — tj. czy koszt
  kreacji μ przechodzi przez 0 na jakiejś skończonej konfiguracji?
  (Spodziewany negatyw = formalizacja intuicji „układ relaksuje bez
  końca, bo granicy nie ma w kontinuum" — pełnoprawny wynik.)
- **Q2 (podłoga substratowa):** czy z podłogą pola wyprowadzoną
  z QB-2 (nie ad-hoc!) relaksacja zatrzymuje się na stanie NIESTAŁYM
  (strukturalnym) — kandydacie na stan metametryczny? Trzy rozłączne
  wyjścia (§3).
- **Q3 (warunkowe):** czy stan zrelaksowany z Q2 jest stabilny:
  ω²_min ≥ 0 na tym tle (pierwszy pomiar stabilności względem
  WŁAŚCIWEGO punktu odniesienia)?

## 2. Model ZAMKNIĘTY

- **Sektor:** kanoniczny tachionowy (K=g⁴, U=g⁷/7−g⁸/8), jak w cyklach
  3D; spektra weight-1 bez regularyzacji.
- **Podłoga substratowa (Q2/Q3):** g_floor wyprowadzone z progu QB-2
  przez mapowanie pole↔Φ z dodatekB (φ=(Φ/Φ₀)^{1/2} — dokładne
  mapowanie i wartość FROZEN w `Phase_method_decisions.md` PRZED
  startem, z cytatem równań źródłowych); **skan wrażliwości po całym
  zakresie QB-2: Φ_c/Φ_vac ∈ {0.197, 0.298, 0.331}** (3 wartości
  zalockowane). Implementacja podłogi: kara/projekcja zachowująca
  strukturę wariacyjną (wybór FROZEN; wymóg: gate energii ewolucji
  dalej osiągalny).
- **Konfiguracje startowe (zalockowane):** (i) izolowany soliton μ
  radialny (kotwica z 3d-canonical: λ_min=−1.646589 [INPUT]);
  (ii) sieć 2π z `Phase2_backgrounds3d.npz` (READ-ONLY); (iii) próżnia
  g≡1 + małe zaburzenie losowe (ziarno zalockowane: seed=20260901,
  amplituda 1e−3) — start „genezowy".
- **Rejestr WEJŚĆ:** g₀_μ=φ·0.90548, próg QB-2 (3 wartości), seed,
  β=γ=1 — flagować.

## 3. Fazy i kryteria (zalockowane)

### Phase 1 — Q1: równanie stanu kreacji (kontinuum; sympy + tania kwadratura radialna)
- **P1a (analitycznie, sympy):** struktura U: U(0)=0 < U(1)=1/56,
  U″(1)=−1; wniosek formalny o kierunku relaksacji; brak minimum
  lokalnego U w (0,1).
- **P1b (kwadratura radialna):** ΔE_create(soliton μ | próżnia)
  = E[g_sol]−E[g≡1] na R=60 (h z tabeli kotwicy) — znak i wartość;
  analogicznie względem stanu pustego.
- **P1c (sieć):** ε(2π) = (E_cell − E_cell[g≡1])/V dla tła 2π
  (z npz, bez relaksowania) — znak.
- **Werdykt Q1:** granica metametryczna w kontinuum ISTNIEJE ⟺
  istnieje w policzonym zbiorze konfiguracja z μ=0 oddzielająca
  reżim opłacalnej kreacji od kosztownej. Spodziewane NIE — wtedy
  Q1-NEG (formalizacja relaksacji bez granicy) i przejście do Q2.

### Phase 2 — Q2: relaksacja z podłogą substratową (RACHUNEK CENTRALNY)
- **P2a (gate maszynerii):** ewolucja/relaksacja z podłogą na próżni
  bez zaburzeń: MUSI zostać w g≡1 (osiągalny FAIL); gate energii
  ≤1e−6 tam, gdzie podłoga nieaktywna.
- **P2b (relaksacja):** tłumiony gradient flow (schemat FROZEN) dla
  3 konfiguracji startowych × 3 wartości podłogi × 2 siatki
  (radialnie h z tabeli kotwicy; 3D N∈{32,48} dla sieci 2π);
  do stacjonarności ‖ġ‖∞ ≤ 1e−8 LUB do detekcji nukleacji LUB do
  t_max=200 (INCOMPLETE raportowany).
- **DETEKTOR NUKLEACJI (pre-rejestrowany):** obiekt = spójny region
  g < g_thr := (1+g_floor)/2, odseparowany od zasianych (etykietowanie
  spójnych komponent co Δt=1); **NUKLEACJA = wzrost liczby obiektów
  ponad liczbę zasianych, utrzymany ≥10 jednostek czasu**.
- **KRYTERIUM AUTORA (pre-rejestrowane jako POZYTYW):** jeśli rachunek
  „załamuje się" przez GENEROWANIE OBIEKTÓW — tj. nukleacja wg
  detektora, ZBIEŻNA (obecna na obu siatkach i przy dt/2; liczba
  obiektów w chwili detekcji zgodna ±1) — to werdykt
  **Q2-PASS-NUCLEATION: reżim bezkosztowej kreacji potwierdzony**
  (hipoteza „wielki wybuch trwa na granicy metryki" dostaje pierwszy
  nośnik numeryczny). Załamanie NIE-nukleacyjne (NaN, brak zbieżności
  detektora między siatkami) = INCONCLUSIVE, nie pozytyw.
- **Pozostałe wyjścia:** **Q2-PASS-STATIC** — relaksacja staje na
  stanie NIESTAŁYM (‖g−const‖∞ ≥ 0.05, stacjonarnym, zbieżnym
  siatkowo ≤5e−3) = stan metametryczny istnieje jako konfiguracja;
  **Q2-FAIL** — wszystkie starty relaksują do stanu jednorodnego
  (g≡1 lub g≡g_floor) na wszystkich podłogach — granica substratowa
  nie generuje struktury.
- **P2c (kontrola sektora stabilnego, nieusuwalna):** identyczna
  procedura w m²=+γ: relaksacja do próżni, ZERO nukleacji — osiągalny
  FAIL detektora (fałszywe alarmy dyskwalifikują maszynerię).

### Phase 3 — Q3 (TYLKO przy Q2-PASS-STATIC): stabilność w stanie metametrycznym
- Widmo ω²_min wokół stanu zrelaksowanego (maszyneria weight-1 jak
  w cyklach 3D; dla tła radialnego — operator radialny; dla sieci —
  Bloch Γ/X/M/R); mody zerowe (translacje/podłoga) zidentyfikowane
  PRZED interpretacją; zbieżność ≤0.05·max(|ω²_min|,0.1).
- **Q3-PASS:** ω²_min ≥ −tol (tol=1e−3) zbieżnie — pierwszy stabilny
  obiekt sektora tachionowego, zmierzony względem właściwej granicy.
- **Q3-FAIL:** ω²_min < 0 zbieżnie — nawet stan metametryczny
  niestabilny.
- Przy Q2-PASS-NUCLEATION Phase 3 NIE ma tła statycznego — zamiast
  tego deskryptywna charakteryzacja kaskady (tempo kreacji, rozkład
  rozmiarów) BEZ progów (zero skanowania do celu).

## 4. Forbidden moves
1. Zero zmian kryteriów/progów/detektora/seeda po starcie; korekty
   wyłącznie dla udokumentowanego błędu implementacji (correction_note
   PRZED użyciem; pierwotne outputy zachowane).
2. Detektor nukleacji NIE może być zmieniany po pierwszym biegu;
   kontrola P2c (zero fałszywych alarmów) nieusuwalna.
3. Podłoga WYŁĄCZNIE z derywacji QB-2/dodatekB (3 zalockowane
   wartości); zakaz strojenia g_floor do wyniku.
4. Tła z npz READ-ONLY; kotwica [INPUT], niestrojona.
5. Wyniki negatywne wprost; INCONCLUSIVE ≠ pozytyw (litera §3);
   zakaz reinterpretacji Q2-FAIL jako „prawie nukleacja".
6. Rdzeń `.tex` NIETKNIĘTY; STATE.md nieedytowany; git nieużywany;
   katalogi innych cykli tylko do odczytu.
7. Decyzje implementacyjne FROZEN przed startem faz; higiena ścieżek
   (bez `cd`, weryfikacja materializacji).

## 5. Deliverables
`Phase_method_decisions.md` (w tym mapowanie g_floor z cytatem),
`Phase1_creation_cost.py`+output, `Phase2_floor_relax.py`+output
(+ npz stanów zrelaksowanych), `Phase3_metametric_spectrum.py`+output
(warunkowy), `Phase_FINAL_close.md` (werdykty Q1/Q2/Q3 wg litery),
`NEEDS.md` (user-gated), `README.md` (log faz).

## 6. Drzewo decyzyjne
```text
Q1-POS (granica w kontinuum istnieje) → charakteryzacja + Q2 nadal
   (podłoga może przesuwać granicę)
Q1-NEG → formalizacja „relaksacji bez granicy w kontinuum" → Q2
P2a FAIL → STOP (maszyneria podłogi nieważna)
P2c FAIL (fałszywe nukleacje) → STOP (detektor nieważny)
Q2-PASS-NUCLEATION → NEEDS: hipoteza granicy metametrycznej dostaje
   nośnik numeryczny; kandydat: program genezy Γ+s_i (post-close
   blocked-soliton-bang) jako osobny cykl; re-interpretacja Q-FAIL-i
   1D/3D jako relaksacji (dopisek user-gate)
Q2-PASS-STATIC → Phase 3; przy Q3-PASS: PILNE — pierwszy stabilny
   obiekt sektora tachionowego (konsekwencje dla znaku W: user-gate)
Q2-FAIL → NEEDS: podłoga substratowa nie generuje struktury w klasie
   zbadanej; granica metametryczna wymaga poziomu 0 wprost (Γ+s_i)
Q3-FAIL → NEEDS: nawet stan zrelaksowany niestabilny — hipoteza
   granicy jako punktu odniesienia osłabiona (raport wprost)
```

---

**LOCK ZAMKNIĘTY 2026-09-01. Zmiany poniżej tej linii po starcie obliczeń
= forbidden move (poza datami realizacji faz).**
