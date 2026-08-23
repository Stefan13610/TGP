---
title: "Phase0_balance — DRAFT: B2-prime proba #4 — pelny P2: uniwersalnosc w znaku kretu na SIATCE SKOSNEJ (tlo asymetryczne w krecie, symetryczne w pozycjach; L=256)"
date: 2026-07-05
type: phase0-draft
tgp_owner: research/op-asymmetric-lattice-2026-07-05
status: LOCKED
anti_lakatos_lock: ACTIVE (autoryzacja autora: polecenie realizacji cyklu, 2026-07-12; lock przed kodem; wykonawca: osobny agent, nowa sesja — spelnione)
tags: [gravity-bridge, branch-B2prime, oblique-lattice, vortex-lensing, circulation-channel, universality, P2, L256]
related:
  - "[[../op-shortwave-lattice-2026-07-05/Phase_FINAL_close.md]]"
  - "[[../op-lattice-background-2026-07-05/Phase_FINAL_close.md]]"
  - "[[../op-vortex-refraction-2026-07-05/Phase_FINAL_close.md]]"
---

# Phase 0 — DRAFT: B2-prime proba #4 — pelny P2 na siatce skosnej

## 0. Hipoteza robocza i miejsce w programie

Stan lancucha (wszystko zmierzone, nic do wiary):
- P1 (POLICZALNOSC) ROZSTRZYGNIETE POZYTYWNIE (proba #3 6/6 +
  proba #2 5/6): eikonal substratu przewiduje kat ugiecia w calej
  swojej domenie (b_eff >= 1.5 lambda), ratio 0.79-1.57; znak KU
  defektowi w 19/19 waznych runow; granica stosowalnosci zmierzona
  (b_eff ~ lambda -> ratio 0.40-0.45 niezaleznie od L).
- P2 (UNIWERSALNOSC W ZNAKU KRETU) OTWARTE: w probach #2 i #3
  alpha_odd = 0 Z KONSTRUKCJI (odbicie osi wiazki + sprzezenie
  zespolone bylo dokladna symetria tla — pre-rejestrowane
  i potwierdzone do 4e-04 / 4.9e-11 deg). Pomiar kanalu
  cyrkulacyjnego wymaga tla, na ktorym ZADNA symetria nie wiaze
  alpha(+b) z alpha(-b).
- Regula decyzyjna proby #3 wskazuje ten op jako pierwszy
  z trzech kandydatow ("pelny P2 na tle asymetrycznym").

Kluczowa obserwacja konstrukcyjna (fundament tego cyklu):
odbicie wzgledem osi wiazki zlozone ze sprzezeniem zespolonym
ZAWSZE zachowuje ladunki wirow (obie operacje z osobna odwracaja
kret; zlozenie go zachowuje), a odwzorowuje pozycje y -> 2y_osi-y.
Zatem alpha_odd jest wiazane symetria wtedy i tylko wtedy, gdy
POZYCJE wirow sa lustrzane wzgledem osi wiazki Z ZACHOWANIEM
LADUNKOW. Mozna wiec zbudowac tlo, ktore:
1. ma pozycje wirow (prawie) lustrzane -> eikonal (slepy na kret,
   zalezny tylko od |psi|) przewiduje Delta_eik ~ 0,
2. ma wzor LADUNKOW maksymalnie lamiacy lustro -> sprzezenie
   pakietu z cyrkulacja tla (czlon ~ i grad(theta) . grad(u),
   pierwszy rzad w grad theta) NIE jest niczym skasowane,
3. ma sile netto = 0 DOKLADNIE na kazdym wirze (kasowanie
   z SYMETRII INWERSYJNEJ otoczenia, nie z D4).

To jest **eksperyment zerowy typu Eotvosa dla kanalu
geometrycznego**: Delta_meas = alpha(+b) - alpha(-b) rozny od
Delta_eik ~ 0 = bezposrednia detekcja kanalu sprzezonego z kretem
(analogia grawitomagnetyzmu); Delta_meas zgodne z Delta_eik =
slepota na kret ZMIERZONA tam, gdzie nic jej nie wymusza.

Tlo: **SIATKA SKOSNA** (szachownica scinana o L/4 = 64 na rzad):
wiry w wezlach siatki Bravais a1 = (128, 0), a2 = (64, 64),
ladunek (-1)^(n1+n2), 8 wirow na torusie 256 x 256:

```text
  y=-64:  +1 (-64,-64)   -1 (+64,-64)     <- rzad wiazki (rozprasza +1)
  y=0:    -1 (0,0)       +1 (-128,0)
  y=+64:  +1 (+64,+64)   -1 (-64,+64)
  y=-128: +1 (0,-128)    -1 (+128,-128)
```

Wlasnosci (do weryfikacji w Phase0_check):
- POZYCJE: dokladnie lustrzane wzgledem y=-64 (rzad y=0 i y=-128
  w odleglosci 64, oba z soczewka w x=0 w oknie pomiaru;
  x = -128 == +128 na torusie),
- LADUNKI: lustro+sprzezenie odwraca ladunki rzedow posrednich
  (y=0 i y=-128) -> ZADNA symetria nie wiaze alpha(+b) z alpha(-b),
- SILA: otoczenie kazdego wiru jest inwersyjnie symetryczne
  z ZACHOWANIEM ladunkow (parzystosc (-1)^(n1+n2) jest funkcja
  parzysta wektora) -> F_netto = 0 DOKLADNIE (kasowanie parami
  v/-v), takze na poziomie pelnego pola (inwersja wokol kazdego
  rdzenia jest symetria konfiguracji; grid ja respektuje —
  nominaly na wezlach),
- n_tot = 0; gestosc 2x proby #3 (najblizszy sasiad 90.5).

Pytania:

> (P2-a) UNIWERSALNOSC IZOTROPOWA: czy na tle BEZ symetrii
>      lustrzanej kat ugiecia w KAZDYM punkcie +/-b zgadza sie ze
>      SLEPYM NA KRET eikonalem substratu (znak zgodny z predykcja,
>      pasmo [0.5, 2.0])? To test policzalnosci superpozycji
>      8 soczewek bez oslony symetrii.
> (P2-b) KANAL CYRKULACYJNY (RDZEN): czy asymetria
>      Delta_meas = alpha(+b) - alpha(-b) jest w calosci wyjasniona
>      eikonalem (Delta_eik ~ 0 z konstrukcji pozycji), czy zostaje
>      skladowa sprzezona ze znakiem kretu tla?

G5 PASS (a+b) = P2 rozstrzygniete pozytywnie: kanal geometryczny
SLEPY NA ZNAK KRETU rowniez tam, gdzie nic tego nie wymusza —
druga wlasnosc grawitacyjna (uniwersalnosc) zmierzona, nie
zalozona. G5b FAIL przy G5a PASS = ODKRYCIE kanalu cyrkulacyjnego
(wielkosc zmierzona) — nie falsyfikacja mostu, lecz nowy kanal do
charakteryzacji (osobny op).

## 1. Predykcje algebraiczne do zalockowania (PRZED kodem)

### 1.1 Rownanie #1 — wspolczynnik zalamania (bez zmian)

```text
n^2(x,y) = (omega^2 - V''(|psi_bg|)) / (omega^2 - V''(s*))
A_t = V'''(s*) * s* * xi^2 = 2.5690;  m_TV = 0.9374
k0/lambda/v_g:  omega=1.1: 0.8140 / 7.72 / 0.3700
                omega=1.5: 1.6561 / 3.79 / 0.5520
n^2 zalezy WYLACZNIE od |psi_bg| — eikonal slepy na kret
Z DEFINICJI; kanal cyrkulacyjny (jesli istnieje) siedzi w czlonie
~ 2 i kappa grad(theta_bg).grad(u) rownania pakietu, ktorego
eikonal NIE zawiera. Ten cykl mierzy jego skutek katowy.
```

### 1.2 Rownanie #2 — bilans sil i zlamanie lustra (fundament tla)

```text
sila: otoczenie wiru w siatce Bravais o ladunku (-1)^(n1+n2):
  sasiad w v = n1 a1 + n2 a2 ma ladunek q(v) = (-1)^(n1+n2);
  q(-v) = q(v) -> wklady par v/-v kasuja sie DOKLADNIE (inwersja);
  Phase0_check weryfikuje numerycznie (pary inwersyjne do
  zbieznosci; UWAGA: brak kasowania per powloka-kwadrat jak w D4 —
  kasowanie jest parami, sumowac W PARACH v/-v).
lustro (y -> -128-y wzgl. y=-64) + sprzezenie: pozycje przechodza
  na siebie, ladunki rzedow y=0 i y=-128 sie ODWRACAJA
  -> tlo NIE jest odwzorowane na siebie: alpha_odd niewiazana.
symetrie POZOSTALE (pre-rejestracja testow implementacji G6a):
  S-inv: inwersja wokol rdzenia A (r -> 2A - r) jest DOKLADNA
    symetria tla z zachowaniem ladunkow -> run(+b, pakiet
    prawobiezny) = run(-b, pakiet LEWOBIEZNY) co do zaokraglen.
  S-conj: sprzezenie zespolone = translacja tla o (128, 0)
    -> run na tle sprzezonym z pakietem przesunietym o (128,0)
    = sprzezenie runu bazowego BITOWO.
niestabilnosc siodlowa: tempo NIEZNANE dla siatki skosnej (mody
  scinajace moga byc miekksze niz w szachownicy!) — RYZYKO #1
  cyklu; rozstrzyga bramka G2, nie zalozenie.
```

### 1.3 Rownanie #3 — zaklad eikonalny i skala efektu

Formula ogonowa skalibrowana na ray tracingu proby #3
(t_deg(d) = (180/pi)*(pi/2)*A_t/((omega^2-V''(s*))*d^2);
kontrola: omega=1.1, d=20: formula 1.74 vs trace 1.79 — zgodna
do 3%): t_deg(d) = 698/d^2 (omega=1.1); 169/d^2 (omega=1.5).

```text
orientacyjnie (baza: trace proby #3 dla wkladu rzedu wiazki):
  omega=1.1: alpha(+/-14) ~ 3.6-3.9; (+/-20) ~ 1.6-2.0;
             (+/-24) ~ 1.0-1.5 (wklady rzedow posrednich
             z x=0: t(64-b) i t(64+b) po obu stronach)
  omega=1.5: alpha(+/-8) ~ 3.0; (+/-12) ~ 1.2
Delta_eik = alpha_pred(+b) - alpha_pred(-b): z lustrzanosci
  POZYCJI oczekiwane ~ 0 (roznice tylko przez slad wzoru ladunkow
  w |psi| — czlon (grad theta)^2 w rownaniu amplitudy, drugi
  rzad; orientacyjnie |Delta_eik| < 0.1 deg). WIAZACO liczy
  Phase0_check: ray tracing na zrelaksowanym tlu skosnym,
  OSOBNO dla +b i -b, z czuloscia plaszczyzn (wzorzec D6).
sigma_alpha (z proby #3, ten sam aparat): 0.20 (1.1, b~14) do
  0.014 (1.5, b=16) -> 3*sigma_Delta ~ 0.83 deg (najgorzej,
  1.1/14) i ~0.15 deg (najlepiej, 1.5/8). Najczystsze zero:
  seria 1.5.
```

### 1.4 Budzet czasu (KRYTERIUM, wzorzec prob #2/#3)

```text
geometria: rozpraszacz +1 w (-64,-64); start x0 = -100;
  prog wyzwolenia ~x = -54; przelot do 5. probki:
  omega=1.1: ~tau=140; omega=1.5: ~tau=99
bramka G2: tlo bez pakietu, II rzad, okno tau=600 (12000 krokow);
  kryterium: max przemieszczenie kazdego z 8 rdzeni <= 1.0 na
  tau in [0, 450]; (450, 600] raport poza kryterium
kontrast (zmierzone): szachownica L=256 — 0.0000 przez tau=600;
  para — anihilacja 107.5
```

## 2. Model i parametry (propozycja do LOCKa)

```text
substrat: bez zmian (kappa=0.5, a=0.5, b=1.6, c=1.0; psi in C;
  N=512, L=256, dx=0.5, periodycznie, eps=0.30, S_GUARD=10 flaga,
  seed=20260704)
tlo: SIATKA SKOSNA 8 wirow (sekcja 0); ansatz: suma CZTERECH par
  wspolliniowych poziomych theta_h_pair (konstrukcja zwalidowana):
  rzad y=-64:  pair(+1@(-64,-64), -1@(+64,-64))
  rzad y=0:    pair(+1@(-128,0),  -1@(0,0))
  rzad y=+64:  pair(+1@(+64,+64), -1@(-64,+64))
  rzad y=-128: pair(+1@(0,-128),  -1@(+128,-128))
  (obrazy k in [-2,2]; amplituda: iloczyn 8 profili rdzeni;
  diagnostyka szwu OBOWIAZKOWA — 4 pary to nowa kombinacja;
  relaksacja 2000 krokow dt=0.02; residuum raportowane)
sektor falowy, referencja blizniacza, pakiet, okno pomiaru,
  b_eff, alpha_pred(b_eff): BEZ ZMIAN vs proba #3 (D1-D6 tamtego
  LOCKa; nominaly detektora: 4+4 rdzenie)
NOWE — pakiet LEWOBIEZNY (wylacznie G6a-S-inv): start
  x0 = +? (inwersyjny obraz startu prawobieznego wzgledem A:
  (-28, -64-b)), kierunek -x, okno lustrzane; procedura
  doprecyzowana w sekcji 10 przed kodem
skan glowny (KRYTERIALNY, zamkniety) — PARY +/-b:
  omega=1.1 x b in {+/-14, +/-20, +/-24}   (6 runow)
  omega=1.5 x b in {+/-8, +/-12}           (4 runy)
  wszystkie b_nom >= 1.5*lambda + 0.9 (domena eikonalu:
  1.1: prog 12.48; 1.5: prog 6.59); odleglosc wiazki od soczewek
  rzedow posrednich rowniez >= 1.5 lambda w kazdym punkcie
symetrie (G6a): S-inv: (1.1, -14) lewobiezny vs (1.1, +14);
  S-conj: tlo sprzezone + pakiet przesuniety o (128,0): bitowo
liniowosc (G6b): u0 -> u0/2 przy (1.1, +14)
obserwacje POZA kryteriami: (1.1, +/-28) — strefa najsilniejszych
  wkladow wielosoczewkowych (blizej rzedow posrednich)
metrologia: jak proba #3 (D8): dryft sekularny energii <= 1e-4;
  gamma_core <= 0.01 (okna r<5 wokol WSZYSTKICH 8 rdzeni);
  wrap-guard 0.9L; determinizm bitowy w Phase 1; detektor rdzeni
  4+4 (klastrowanie D14, nominaly z sekcji 0)
budzet krokow: relaksacja 2000; bramka G2 12000 (tau=600);
  runy falowe 6000 (tau=300); brak 5 probek = run niewazny
koszt: jak proba #3 (~4x krok proby #2); bramka ~15 min,
  run ~3-6 min, calosc ~1.5-2 h
```

## 3. Scenariusze

1. `background_check + LIFETIME GATE (G2)`: relaksacja siatki
   skosnej; szew, kret (n_tot=0, 8 rdzeni), residuum; ewolucja II
   rzedu BEZ pakietu, tau=600. BRAMKA: FAIL -> STOP.
   (RYZYKO #1: mody scinajace siatki skosnej.)
2. `dispersion_TV (G3)`: kalibracja k in [0.5, 1.7],
   m in {22, 30, 41, 55, 68} (powtorka pomiaru proby #3 na tym
   samym gridzie — cytowanie zabronione, pomiar tani).
3. `deflection(+/-b, omega) (G4, G5a, G5b)`: skan glowny 5 par;
   b_eff, b_eff/lambda per run; alpha vs alpha_pred(b_eff);
   Delta_meas vs Delta_eik per para.
4. `symmetry (G6a)`: S-inv (lewobiezny) i S-conj (sprzezenie
   + translacja).
5. `linearity (G6b)`: u0/2 przy (1.1, +14).
6. obserwacje poza kryteriami: (1.1, +/-28).

## 4. Kryteria (propozycja do LOCKa)

- **G1 techniczne**: determinizm bitowy; brak NaN/S_GUARD; dryft
  sekularny energii <= 1e-4; gamma_core <= 0.01; wrap-guard;
  kret tla zachowany (8 rdzeni, n_tot = 0) w kazdym runie.
- **G2 BRAMKA — czas zycia siatki skosnej**: ewolucja II rzedu
  bez pakietu, tau = 600; max przemieszczenie kazdego z 8 rdzeni
  <= 1.0 dla tau <= 450. FALSYFIKATOR: mody scinajace siatki
  skosnej uciekaja szybciej niz 3x przelot -> tlo niezdatne;
  STOP, deliverable = czasy zycia (skosna vs szachownica).
  Bramki nie wolno oslabic po pierwszym runie.
- **G3 kalibracja**: blad <= 5% dla k in [0.5, 1.7] (pomiar,
  nie cytat).
- **G4 istnienie**: |alpha(1.1, +14)| > 3*sigma ORAZ |alpha|
  maleje z |b| PO OBU STRONACH serii 1.1 (|14| > |20| > |24|).
- **G5a — policzalnosc bez oslony symetrii (P2-a)**: dla
  WSZYSTKICH DZIESIECIU punktow: sign(alpha) = sign(alpha_pred)
  ORAZ alpha/alpha_pred(b_eff) in [0.5, 2.0].
  FALSYFIKATOR: superpozycja soczewek nieczytelna z substratu
  na tle asymetrycznym.
- **G5b RDZEN — zero kanalu cyrkulacyjnego (P2-b)**: dla KAZDEJ
  z 5 par: |Delta_meas - Delta_eik| <= 3*sigma_Delta, gdzie
  Delta_meas = alpha(+b) - alpha(-b), Delta_eik =
  alpha_pred(b_eff+) - alpha_pred(b_eff-) (kazda predykcja przy
  wlasnym b_eff — pochlania artefakt detektora ~0.5, errata #4
  proby #2), sigma_Delta = sqrt(sigma+^2 + sigma-^2).
  FALSYFIKATOR (= odkrycie): skladowa katowa sprzezona z kretem
  tla. Regula wzmacniajaca (pre-rejestrowana): FAIL w >= 2 parach
  ze spojnym znakiem nadwyzki = "kanal cyrkulacyjny zmierzony";
  FAIL w 1 parze = "anomalia pojedynczej pary" (raport, bez
  deklaracji odkrycia).
- **G6 symetrie i liniowosc**: (a) S-inv: |alpha(+14,prawy) -
  alpha(-14,lewy)| <= 0.05*|alpha(+14)| ORAZ S-conj bitowo /
  w tej samej tolerancji; (b) u0/2: roznica <= 10%.

## 5. Regula decyzyjna (propozycja)

```text
G2 FAIL -> STOP przed pomiarem; deliverable: czasy zycia
  (skosna vs szachownica); kandydaci: mniejsze scinanie (s=32),
  L=512.
G5a PASS + G5b PASS (+ G4, G6 PASS)
  -> P2 ROZSTRZYGNIETE POZYTYWNIE: kanal geometryczny slepy na
     znak kretu ZMIERZONY na tle bez oslony symetrii; most
     geometryczny ma obie wlasnosci grawitacji (uniwersalny znak
     + uniwersalnosc zrodla). Nastepne opy (osobne locki):
     (i) skalowanie z |n|, (ii) granica newtonowska 2D + G_eff.
G5a PASS + G5b FAIL (wg reguly wzmacniajacej)
  -> ODKRYCIE: kanal cyrkulacyjny istnieje; deliverable:
     wielkosc Delta na 5 parach; osobny op charakteryzacji
     (skalowanie z b, omega, konfiguracja kretu). P2 w brzmieniu
     "slepota na kret" rozstrzygniete NEGATYWNIE — bez zmiany
     statusu P1.
G5a FAIL (wielkosc przy dobrym znaku) -> policzalnosc superpozycji
     zalamuje sie bez oslony symetrii; powrot do analizy.
G5a FAIL (znak) -> anomalia kierunkowa; powrot do analizy.
G3/G4 FAIL -> stop przed interpretacja.
G6a FAIL -> blad implementacji: errata techniczna (dzien runow),
     dopiero potem ewaluacja G5.
```

## 6. Forbidden moves (propozycja)

1. Zero pinningu i zero wiezow; stabilnosc siatki skosnej ma sie
   utrzymac SAMA (rozstrzyga G2).
2. BRAMKI G2 nie wolno oslabic, skrocic ani zreinterpretowac.
3. Zadnych zmian progow/parametrow/procedur po pierwszym runie
   (wyjatki pre-rejestrowane: jedno obnizenie dt; jeden kontrolny
   run dx=0.25).
4. Zestaw kryterialny (5 par) ZAMKNIETY; obserwacje +/-28
   niedolaczalne; zadnego wykluczania par post hoc.
5. Delta_meas NIE wolno raportowac bez Delta_eik obok (surowa
   asymetria zawiera wklad eikonalny i artefakt b_eff detektora).
6. Ewentualnego G5b FAIL nie wolno interpretowac jako
   frame-dragging/Lense-Thirring GR; tau != czas; zero deklaracji
   obserwacyjnych. Nie promowac do core/ bez osobnej autoryzacji.
7. Wynik S-inv/S-conj (G6a) nie jest pomiarem P2 (testy
   implementacji; jedyny pomiar P2 to G5b na parach +/-b
   prawobieznych).

## 7. Czego ten cykl NIE robi

- NIE mierzy skalowania ugiecia z |n| ani granicy newtonowskiej.
- NIE charakteryzuje kanalu cyrkulacyjnego (jesli wykryty —
  osobny op ze skalowaniem).
- NIE powtarza obserwacji superpozycji za soczewka (D16 proby #3
  wykonana) ani rezimu dyfrakcyjnego.
- Zero GR; tau != czas.

## 8. Deliverables

- `Phase0_balance.md` (po autoryzacji: LOCKED + sekcja 10
  pre-code PRZED kodem; wzorzec: sekcja 10 proby #3, D1-D16 +
  nowe D dla: par inwersyjnych bilansu sil, detektora 4+4,
  pakietu lewobieznego, definicji Delta/sigma_Delta)
- `Phase0_check.py`/`Phase0_output.txt` (bilans sil W PARACH
  INWERSYJNYCH; weryfikacja zlamania lustra przez ladunki
  i lustrzanosci pozycji; tlo zrelaksowane + szew 4 par;
  ray tracing WIAZACY alpha_pred(+/-b) i Delta_eik per para
  + czulosc plaszczyzn; tabela projektowa b/lambda i odleglosci
  od soczewek posrednich; budzet czasu)
- `Phase1_gate_and_calibration.py`/`Phase1_output.txt` (G2
  BRAMKA tau=600 na 8 rdzeniach + G3 + background_check)
- `Phase2_deflection.py`/`Phase2_output.txt` (G4, G5a, G5b, G6b
  + obserwacje +/-28)
- `Phase3_symmetry.py`/`Phase3_output.txt` (G6a: S-inv, S-conj)
- `Phase_FINAL_close.md`, `README.md`

## 9. Status

**LOCKED (2026-07-12).** Autoryzacja autora: polecenie realizacji
cyklu (2026-07-12). Wykonawca: OSOBNY AGENT w nowej sesji
(spelnione). Obowiazuje pelny anty-Lakatos. Konstrukcja tla
rozdziela kanaly: pozycje soczewek lustrzane (eikonal:
Delta_eik ~ 0), wzor kretu lamie lustro (kanal cyrkulacyjny
niczym nieskasowany), sila netto = 0 z inwersji (stabilnosc
mierzona bramka G2 — ryzyko modow scinajacych spisane z gory).
Zestaw 5 par +/-b w calosci w domenie eikonalu; kryterium G5b
z pre-rejestrowana regula odrozniajaca "zero" od "odkrycia"
i od "anomalii pojedynczej pary". Doprecyzowania pre-code:
sekcja 10 (spisane PRZED pierwsza linijka kodu; wzorzec:
sekcja 10 LOCKa proby #3, D1-D16, ze zmianami wskazanymi
w draftach sekcji 9: bilans sil PARAMI INWERSYJNYMI v/-v
zamiast powlok D4; detektor 4+4 z nominalami z sekcji 0;
procedura pakietu lewobieznego dla S-inv; definicje
Delta_meas/Delta_eik/sigma_Delta per para — D16-D19).

## 10. Doprecyzowania pre-code (spisane przed pierwsza linijka kodu)

Zadne z ponizszych nie zmienia progow G1-G6, modelu ani predykcji —
to ujednoznacznienie procedur. D1-D15 przeniesione WPROST z sekcji 10
proby #3 (op-shortwave-lattice) z adaptacja do tla 8 wirow siatki
skosnej (detektor 4+4, gamma_core na 8 rdzeniach, brak runu
superpozycji); D16-D19 NOWE (zmiany wskazane w sekcji 9 draftu).

- **D1 (integrator):** identycznie jak proby #2/#3 — pelny sektor
  zespolony II rzedu, leapfrog pozycyjny `psi_{n+1} = 2 psi_n -
  psi_{n-1} + dt^2 * RHS(psi_n)`, `RHS(psi) = kappa*Lap(psi) -
  psi*(a - b|psi| + c|psi|^2)`; dt = 0.05. `psi_tau` do obserwabli
  = roznica CENTRALNA `(psi_{n+1} - psi_{n-1})/(2 dt)`.
- **D2 (pakiet PRAWOBIEZNY, start analityczny):** `psi(0) = psi_bg
  + e^{i theta_bg} * u0 * env * cos(k0 (x - x0))`, `psi(-dt) =
  psi_bg + e^{i theta_bg} * u0 * env * cos(k0 (x - x0) + omega dt)`;
  env = exp(-(x-x0)^2/(2 sx^2) - (y-y0)^2/(2 sy^2)), sx = sy = 8,
  u0 = 1e-3, x0 = -100, y0 = -64 + b_nom (b_nom ZE ZNAKIEM: b<0 =
  pakiet po stronie -y rzedu wiazki); k0 = sqrt((omega^2 -
  V''(s*)) / kappa).
- **D3 (referencja blizniacza):** identycznie: delta_psi = psi -
  psi_ref, psi_ref ewoluuje TYM SAMYM integratorem z psi_bg bez
  pakietu, w lockstep. Pozycje rdzeni do b_eff, okna i gamma_core
  czytane z psi_ref (detektor 4+4, D14). To NIE jest pinning.
- **D4 (okno i wyzwolenie):** wyzwolenie CENTROIDEM:
  x_c(|delta_psi|^2) > x1 + 10, x1 = BIEZACA (z psi_ref) pozycja x
  rdzenia A (+1 rzedu y=-64, nominal (-64,-64)). W chwili
  wyzwolenia progi ZAMRAZANE: okno = {x1_frozen + 10 < x <
  x2_frozen - 10}, x2 = pozycja x rdzenia -1 rzedu y=-64 (nominal
  (+64,-64), z psi_ref). Potem 5 probek co Delta_tau = 4:
  P = -Re sum(conj(delta_psi_tau) * grad(delta_psi)) dx^2 w oknie;
  usrednienie WEKTOROWE (Px, Py); `alpha := -sign(b_nom) *
  atan2(Py, Px)` [deg] — alpha > 0 = KU rdzeniowi A dla OBU znakow
  b; sigma_alpha = odch. std probek katowych.
- **D5 (b_eff):** min po tau (probki co 10 krokow, do wyzwolenia)
  odleglosci periodycznej centroidu |delta_psi|^2 od rdzenia A
  (pozycja z psi_ref). alpha_pred(b_eff) = ray tracing przy
  y_start = -64 + sign(b_nom) * b_eff na tym samym zrelaksowanym
  tle. b_eff/lambda raportowane per run; zadnego wykluczania
  punktow post hoc.
- **D6 (eikonal na rzeczywistym tle):** n^2 = (omega^2 -
  V''(|psi_bg|)) / (omega^2 - V''(s*)) na PELNYM zrelaksowanym tle
  siatki skosnej (8 soczewek + obrazy periodyczne); bilinearna
  interpolacja n^2 i gradientu; ray tracing hamiltonowski RK4,
  dsigma = 0.02, start x = -100, y = -64 + sign(b)*b_eff, kierunek
  +x, |p(start)| = n(start). Kat czytany w plaszczyznie x_eval =
  centroid pakietu w 3. probce; CZULOSC raportowana na krancach
  zamrozonego okna. Phase0_check: czulosc plaszczyzn x_eval in
  {-54, -27, 0, +27, +54} dla +b i -b OSOBNO na kazdym punkcie.
  Przechwyt: r < 0.6 od ktoregokolwiek z 8 rdzeni -> NaN (raport).
- **D7 (G3, dyspersja):** identycznie: tlo jednorodne s*,
  zaburzenie rzeczywiste u0*cos(k x), m in {22, 30, 41, 55, 68}
  (k in [0.540, 1.669]), okno tau = 100, zliczanie przejsc przez
  zero; prog 5%; dyspersja sieciowa analityczna obok.
- **D8 (G1, metryki):** dryft sekularny energii <= 1e-4 (ostatnia
  vs pierwsza cwiartka); gamma_core = slope fitu ln RMS |delta_psi|
  w SUMIE okien r < 5 wokol WSZYSTKICH OSMIU rdzeni (pozycje
  z psi_ref), fit po tranzycie, prog 0.01. Kret tla:
  n_tot(psi_ref, koniec) = 0 ORAZ 8 rdzeni obecnych (4 klastry na
  znak, D14). Wrap-guard: |przemieszczenie centroidu od startu| <
  0.9 L = 230.4. S_GUARD = 10 flaga. Determinizm bitowy w Phase 1
  (relaksacja + run falowy). Zapisy co 10 krokow; seed = 20260704.
- **D9 (budzet krokow):** relaksacja tla 2000 krokow dt = 0.02;
  bramka G2: 12000 krokow (tau = 600) bez pakietu; runy falowe:
  limit 6000 krokow (tau = 300). BEZ wyjatku superpozycji (D16
  proby #3 nie wystepuje w tym cyklu). Brak 5 probek = run
  niewazny do kryterium, przyczyna zapisana.
- **D10 (G6b, liniowosc):** |alpha(u0) - alpha(u0/2)| <= 0.1 *
  max(|alpha(u0)|, |alpha(u0/2)|) przy (omega = 1.1, b = +14).
- **D11 (G6a — pre-rejestracja procedur symetrii):**
  (a) **S-inv:** tlo TO SAMO psi_bg (inwersja wokol rdzenia A,
  r -> 2A - r, jest dokladna symetria konfiguracji z zachowaniem
  ladunkow — sekcja 1.2); run testowy = (omega=1.1, b_nom=-14)
  z pakietem LEWOBIEZNYM wg D17; baseline = (1.1, +14) prawobiezny
  (powtorka runu Phase 2). Kryterium: |alpha(+14, prawy) -
  alpha(-14, lewy)| <= 0.05 * |alpha(+14)|. Oczekiwane
  zaokraglenia, nie bitowosc (ansatz nie jest bitowo inwersyjnie
  symetryczny: obciecie sumy obrazow k in [-2,2]); diagnostyka:
  max|psi_bg - I[psi_bg]| raportowana (I = odbicie indeksow
  (i,j) -> (256-i, 256-j) mod 512 = inwersja wokol nominalnego A).
  (b) **S-conj:** tlo sprzezone + pakiet przesuniety o (128, 0),
  zaimplementowane w UKLADZIE PRZESUNIETYM: tlo_c :=
  np.roll(conj(psi_bg), N/2, axis=0) (sprzezenie tla = translacja
  o (128,0), sekcja 1.2; roll o N/2 komorek = dokladna operacja
  siatkowa), pakiet standardowy D2 w (-100, -64+14) ukladu
  przesunietego = (+28, -64+14) ukladu bazowego; procedura
  pomiarowa identyczna z D4. Oczekiwanie: max|tlo_c - psi_bg| = 0
  -> run bitowo identyczny z baseline; w przeciwnym razie rownosc
  alpha w tolerancji 5%|alpha|; diagnostyka bitowa (max|tlo_c -
  psi_bg|, roznice probek) raportowana. Wynik G6a nie jest
  pomiarem P2 (forbidden move #7).
- **D12 (domena eikonalu i geometria):** wszystkie 10 punktow
  kryterialnych spelnia b_nom >= 1.5*lambda + 0.9 ORAZ odleglosc
  wiazki od soczewek rzedow posrednich (64 - |b_nom| do blizszej,
  64 + |b_nom| do dalszej; soczewki w x=0 rzedow y=0 i y=-128) >=
  1.5*lambda (tabela projektowa w Phase0_check). Obserwacje
  (1.1, +/-28) POZA kryteriami, niedolaczalne (forbidden move #4);
  zadnego wykluczania par post hoc, nawet gdyby b_eff/lambda < 1.5.
- **D13 (narzedzia):** zero `cd` w powloce; sciezki od korzenia
  vaulta; kazdy skrypt pisze output do `PhaseX_output.txt`
  w katalogu opa (python -u, przekierowanie >).
- **D14 (detektor rdzeni 4+4):** winding plakietowy W; komorki
  z W*sign > 0.5 klastrowane po odleglosci periodycznej (promien
  <= 8 komorek od komorki bazowej klastra); wymagane DOKLADNIE
  4 klastry na znak (inaczej flaga core_lost); pozycja klastra =
  centroid deficytu Phi (wagi max(0, Phi*/2 - Phi)) w oknie 13x13
  wokol komorki bazowej; przypisanie 4 klastrow do 4 nominalow
  danego znaku po najmniejszej sumie odleglosci periodycznych
  (pelny przeglad 24 permutacji). Nominaly (sekcja 0):
  +1: (-64,-64), (-128,0), (+64,+64), (0,-128);
  -1: (+64,-64), (0,0), (-64,+64), (+128,-128).
  Klucz (+1, 0) = rdzen A (rozpraszajacy); (-1, 0) = rdzen x2
  rzedu wiazki.
- **D15 (bramka G2 — procedura):** start z zrelaksowanego tla,
  psi(-dt) = psi(0) (zerowa predkosc), dynamika II rzedu D1,
  12000 krokow (tau = 600) BEZ pakietu; pozycje 8 rdzeni (D14) co
  10 krokow; kryterium: max po rdzeniach i po tau <= 450
  odleglosci periodycznej od pozycji poczatkowej <= 1.0;
  (450, 600] raportowane poza kryterium. gamma_s: fit liniowy
  ln(max_disp) vs tau na odcinku max_disp in [0.05, 1.0] (jesli
  osiagalny). Bramki nie wolno oslabic, skrocic ani
  zreinterpretowac; G2 FAIL -> STOP, zero runow ugiecia.
- **D16 (NOWE — bilans sil PARAMI INWERSYJNYMI, Phase0_check):**
  sasiedzi rdzenia w v = n1 a1 + n2 a2 (a1 = (128,0),
  a2 = (64,64)), ladunek q(v) = (-1)^(n1+n2); wklad pary {v, -v}:
  q(v)*[F(v) + F(-v)] = 0 DOKLADNIE (q(-v) = q(v) — parzystosc
  jest funkcja parzysta; F nieparzysta w v). Sumowac WYLACZNIE
  parami (reprezentant: n1 > 0 lub (n1 = 0 i n2 > 0)), grupowanie
  w powloki indeksowe K = max(|n1|, |n2|) = 1..8; raport per para
  (max residuum) i per powloka. Kontrola negatywna: obciecie
  POLPLASZCZYZNOWE (niezamkniete na inwersje) NIE kasuje sie —
  kasowanie pochodzi z par v/-v, nie z symetrii D4. NIE stosowac
  petli powlokowej proby #3 bez parowania (pulapka z README).
- **D17 (NOWE — pakiet LEWOBIEZNY, wylacznie G6a-S-inv):** start
  (x0L, y0L) = (-28, -64 - b) przy b = 14 (obraz inwersyjny startu
  prawobieznego (−100, −64+b) wzgledem nominalnego A = (-64,-64));
  psi(0) = psi_bg + e^{i theta_bg} * u0 * env * cos(k0 (x0L - x)),
  psi(-dt) = psi_bg + e^{i theta_bg} * u0 * env * cos(k0 (x0L - x)
  + omega dt) — faza +omega*dt; pakiet biegnie w -x; env jak D2
  wokol (x0L, y0L). Centroid x liczony PERIODYCZNIE (srednia
  kolowa — pakiet przekracza szew x = +/-128). Wyzwolenie: s :=
  (x1 - x_c) mod L in (10, 60) (obraz warunku x_c > x1 + 10; x1 =
  pozycja x rdzenia A z psi_ref); w chwili wyzwolenia okno
  ZAMRAZANE: 10 < (x1_frozen - x) mod L < ((x1_frozen - x2_frozen)
  mod L) - 10 (lustro okna D4 wzgledem inwersji wokol A); 5 probek
  co Delta_tau = 4, usrednienie wektorowe; `alpha := -sign(b_nom)
  * atan2(Py, -Px)` [deg] (alpha > 0 = KU rdzeniowi A; skladowa
  podluzna -Px > 0 dla ruchu w -x); b_eff jak D5 (odleglosc
  periodyczna od rdzenia A). Pakiet lewobiezny NIE wystepuje
  w zadnym runie kryterialnym (forbidden move #7).
- **D18 (NOWE — definicje G5b, per para):** Delta_meas :=
  alpha(+b) - alpha(-b) (oba runy PRAWOBIEZNE, Phase 2);
  Delta_eik := alpha_pred(b_eff(+b)) - alpha_pred(b_eff(-b)) —
  kazda predykcja przy WLASNYM b_eff swojego runu i wlasnej
  plaszczyznie ewaluacji (centroid 3. probki swojego runu;
  artefakt detektora ~0.5 w b_eff sie odejmuje — errata #4 proby
  #2); sigma_Delta := sqrt(sigma_+^2 + sigma_-^2). Kryterium G5b:
  |Delta_meas - Delta_eik| <= 3*sigma_Delta dla KAZDEJ z 5 par;
  regula wzmacniajaca wg sekcji 4 (FAIL w >= 2 parach ze spojnym
  znakiem nadwyzki = kanal cyrkulacyjny zmierzony; FAIL w 1 parze
  = anomalia pojedynczej pary). Delta_meas NIGDY nie jest
  raportowana bez Delta_eik obok (forbidden move #5).
- **D19 (NOWE — kolejnosc realizacji = protokol):** Phase0_check
  (rachunki wiazace) -> Phase 1 (background_check + BRAMKA G2 +
  G3; G2 FAIL -> STOP calego cyklu, zero runow ugiecia) ->
  Phase 2 (skan 5 par + G6b + obserwacje +/-28) -> Phase 3 (G6a)
  -> FINAL close. Interpretacja kazdej fazy wylacznie po
  komplecie jej wynikow; erraty w dniu runow, przed interpretacja.
