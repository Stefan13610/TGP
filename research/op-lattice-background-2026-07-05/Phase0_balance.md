---
title: "Phase0_balance — DRAFT: B2-prime na siatce szachownicowej — ilosciowa refrakcja na defekcie na tle o zerowej sile netto"
date: 2026-07-05
type: phase0-draft
tgp_owner: research/op-lattice-background-2026-07-05
status: LOCKED
anti_lakatos_lock: ACTIVE (autoryzacja autora: polecenie realizacji cyklu, 2026-07-05; lock przed kodem; wykonawca: osobny agent, nowa sesja — spelnione)
tags: [gravity-bridge, branch-B2prime, lattice-background, vortex-lensing, refraction, eikonal, lens-superposition]
related:
  - "[[../op-vortex-refraction-2026-07-05/Phase_FINAL_close.md]]"
  - "[[../op-stationary-background-2026-07-05/Phase_FINAL_close.md]]"
  - "[[../op-goldstone-mediator-2026-07-04/Phase_FINAL_close.md]]"
---

# Phase 0 — DRAFT: B2-prime na siatce szachownicowej

## 0. Hipoteza robocza i miejsce w programie

Stan lancucha (wszystko zmierzone, nic do wiary):
- B1: materia = DEFEKTY fazowe (wiry) w wygenerowanej przestrzeni
  Phi*; kanal fazowy dlugozasiegowy, ale ladunkowy.
- op-stationary-background: sektor falowy wokol wiru CZYSTY
  (lambda_full = +0.0008); wir = soczewka SKUPIAJACA, n(r) > 1,
  ogon ~ A_t/r^2, SLEPA na znak kretu.
- op-vortex-refraction (B2-prime, proba #1): pomiar NIEWYKONANY —
  para (+1,-1) pod dynamika II rzedu kolabuje INERCYJNIE i anihiluje
  przy tau = 107.5 (< przelot pakietu ~120); kolaps fizyczny
  (kontrola dt/2: 0.0074). Kwazi-stacjonarnosc gradientowa NIE
  przenosi sie na dynamike falowa. P1/P2 OTWARTE. Aparat pomiarowy
  (referencja blizniacza, b_eff, eikonal na rzeczywistym tle,
  kalibracja dyspersji 1.21%) przetestowany i sprawny.

Ten cykl usuwa JEDYNA zidentyfikowana przeszkode — tlo o krotkim
czasie zycia — konstrukcja, nie sila: **szachownica 4 wirow**
(+1,-1,-1,+1) w rozstawie L/2. Na torusie obrazy periodyczne czynia
z niej nieskonczona siec szachownicowa: kazdy wir lezy w punkcie
o symetrii D4 swojego otoczenia, wiec **sila netto znika DOKLADNIE**
(kasowanie z symetrii, nie z odleglosci). Konfiguracja jest siodlem
(mody scinajace), ale dyskretyzacja deterministyczna zachowuje
symetrie do precyzji maszynowej — niestabilnosc rosnie z szumu
zaokraglen ~1e-16, wiec czas ucieczki ~ 37/gamma_s e-foldingow.
Czy to wystarcza, ROZSTRZYGA pomiar (bramka G2, nowa — lekcja
z proby #1), nie zalozenie.

Pytania:

> (P1) POLICZALNOSC: czy kat ugiecia pakietu amplitudowego na
>      defekcie zgadza sie z eikonalem z substratu (znak + wielkosc,
>      pasmo [0.5, 2.0]) — teraz w geometrii WIELOSOCZEWKOWEJ
>      (eikonal liczony na pelnym tle siatki: de facto test
>      superpozycji soczewek, punkt (ii) z reguly decyzyjnej
>      op-stationary)?
> (P2-czesciowo) SYMETRIA: czy pomiar respektuje dokladne symetrie
>      tla (alpha(+b) = alpha(-b); wariant lustrzany rownowazny)?
>      UWAGA spisana z gory: na szachownicy czesc cyrkulacyjna
>      alpha_odd znika Z KONSTRUKCJI (odbicie y wzgledem osi wiazki
>      + sprzezenie zespolone jest dokladna symetria tla, a wklady
>      cyrkulacyjne pary wirow rzedu kasuja sie parami). Pelny,
>      nietrywialny test P2 (uniwersalnosc w znaku kretu na tle
>      asymetrycznym) = OSOBNY op.

G5 PASS = pierwszy ilosciowo potwierdzony kanal o znaku grawitacji
w programie (uniwersalnie przyciagajaca, policzalna geometria
propagacji wokol materii-defektow, wraz z superpozycja soczewek).

## 1. Predykcje algebraiczne do zalockowania (PRZED kodem)

### 1.1 Rownanie #1 — wspolczynnik zalamania (bez zmian vs proba #1)

```text
n^2(x,y) = (omega^2 - V''(|psi_bg|)) / (omega^2 - V''(s*))
ogon pojedynczego defektu: n^2 - 1 -> A_t / (r^2 (omega^2 - V''(s*))),
  A_t = V'''(s*) * s* * xi^2 = 2.5690
kanal amplitudowy propaguje dla omega > m_TV = 0.9374;
k0(omega) = sqrt((omega^2 - V''(s*)) / kappa);  v_g = kappa*k0/omega
```

### 1.2 Rownanie #2 — bilans sil na siatce (fundament tla)

```text
sila na wir siatki = suma po sasiadach i obrazach periodycznych;
szachownica o stalej L/2 na torusie L x L: otoczenie kazdego wiru
ma symetrie D4 (najblizsi sasiedzi w +/-x, +/-y: ladunek przeciwny,
przekatne: rownoimienny; obrazy zachowuja wzor) -> F_netto = 0
DOKLADNIE. Phase0_check weryfikuje numerycznie (suma par obrazow
do zbieznosci) oraz cytuje kontrast: para z proby #1 miala
F != 0 i anihilacje tau = 107.5.
niestabilnosc siodlowa: tempo gamma_s NIEZNANE przed pomiarem;
ucieczka z szumu 1e-16 wymaga ~ln(1e16)/gamma_s = 37/gamma_s;
bramka G2 mierzy to wprost.
```

### 1.3 Rownanie #3 — zaklad eikonalny

Orientacyjnie (tlo PARY, proba #1, Phase0_output.txt; alpha > 0 =
KU wirowi): omega=1.1: +14.6/+5.0/+2.7 deg dla b = 8/12/16;
omega=1.3: +5.8/+2.1/+1.1. Zaklad WIAZACY dla G5 liczy Phase0_check
TEGO cyklu: ray tracing na rzeczywistym, zrelaksowanym tle SIATKI
(wszystkie 4 soczewki + obrazy), z ewaluacja kata w plaszczyznie
pomiaru (D-doprecyzowania po locku). Roznice vs tabela pary
raportowane (drugi wir rzedu w x=+32 tez skupia — geometria
dwusoczewkowa na osi przelotu).

UWAGA falowa (bez zmian): lambda(1.1) = 7.72 ~ b — rezim
dyfrakcyjny; stad pasmo G5 szerokie [0.5, 2.0]; b=6 i omega=1.0
poza kryteriami (obserwacje).

### 1.4 Symetrie tla (pre-rejestracja rozstrzygniec G6)

```text
S1 (odbicie osi wiazki): y -> -64-y (wzgledem y=-32) zlozona ze
   sprzezeniem zespolonym odwzorowuje tlo na siebie
   -> alpha(+b) = alpha(-b) DOKLADNIE (do zaokraglen)
   -> alpha_odd = 0 z konstrukcji (raportowane jako test spojnosci,
      NIE jako pomiar uniwersalnosci cyrkulacyjnej)
S2 (lustro ladunkowe): negacja fazy = translacja o (L/2, 0)
   zlozona z symetria -> wariant lustrzany rownowazny; dodatkowo
   pelne sprzezenie zespolone wiaze runy bitowo (D11 proby #1)
```

### 1.5 Budzet czasu (korekta po probie #1 — teraz KRYTERIUM, nie zalozenie)

```text
przelot: start x0=-60, wyzwolenie pomiaru ~x=-20..-15
  -> tau_przelotu ~ 122 (omega=1.1) / 92 (omega=1.3)
bramka G2: tlo bez pakietu, dynamika II rzedu, okno tau = 400
  (>= 3x przelot); kryterium: max przemieszczenie kazdego rdzenia
  <= 1.0 na tau in [0, 300] (pozycje z pola co 10 krokow)
kontrast (zmierzony, proba #1): para — przemieszczenie 14 juz
  przy tau=50, anihilacja 107.5
```

## 2. Model i parametry (propozycja do LOCKa)

```text
substrat: bez zmian (kappa=0.5, a=0.5, b=1.6, c=1.0; psi in C;
  N=256, L=128, dx=0.5, periodycznie, eps=0.30, S_GUARD=10 flaga,
  seed=20260704)
tlo: SZACHOWNICA 4 wirow w rozstawie L/2:
  +1 w (-32,-32)   -1 w (+32,-32)
  -1 w (-32,+32)   +1 w (+32,+32)
  ansatz: suma DWOCH PAR WSPOLLINIOWYCH POZIOMYCH (konstrukcja
  dokladna D1 z B1: rzad y=-32: theta_h_pair((-32,-32),(+32,-32));
  rzad y=+32: theta_h_pair((+32,+32),(-32,+32))); BEZ transpozycji
  -> BEZ artefaktu punktu C; amplituda: iloczyn 4 profili rdzeni;
  diagnostyka szwu obowiazkowa (max skok fazy poza rdzeniami << pi;
  szew y: ogony dipolowe ~exp(-pi) ~ 0.04 rad — raport);
  relaksacja 2000 krokow (dt=0.02)
sektor falowy: PELNY zespolony II rzedu (jak proba #1):
  d^2psi/dtau^2 = kappa*Lap(psi) - psi*(a - b|psi| + c|psi|^2)
  leapfrog pozycyjny dt=0.05; psi_tau = roznica centralna
referencja blizniacza (sprawdzona w probie #1): psi_ref bez pakietu
  w lockstep tym samym integratorem; delta_psi = psi - psi_ref;
  pozycje rdzeni z psi_ref
pakiet: delta_psi(0) = e^{i theta_bg} * u0 * env * cos(k0 (x - x0));
  start analityczny (D2 proby #1); u0 = 1e-3, sx = sy = 8,
  x0 = -60, y0 = -32 + b; wir rozpraszajacy: +1 w (-32,-32)
okno pomiaru (lekcja proby #1 — prog NIE podaza za rdzeniem):
  maska x in (x1+10, x2-10), gdzie x1, x2 = pozycje x wirow rzedu
  y=-32 z psi_ref W CHWILI WYZWOLENIA (potem ZAMROZONA);
  wyzwolenie: centroid |delta_psi|^2 przekracza x1+10;
  5 probek co Delta_tau = 4; usrednienie WEKTOROWE (Px, Py);
  alpha := -sign(b) * atan2(Py, Px) [deg] (alpha > 0 = KU wirowi)
b_eff: min odleglosc centroidu |delta_psi|^2 od rdzenia
  rozpraszajacego (psi_ref), do chwili wyzwolenia
alpha_pred(b_eff): ray tracing na TYM SAMYM zrelaksowanym tle
  siatki, przy b = b_eff; kat promienia ewaluowany w x rownym
  pozycji centroidu pakietu w 3. (srodkowej) probce; czulosc na
  plaszczyzne ewaluacji (krance okna) raportowana
skan glowny: omega in {1.1, 1.3} x b in {8, 12, 16}
symetrie (G6): (omega=1.1, b=8) w wariantach {b=-8} i {tlo
  lustrzane, b=+8}
liniowosc (G6): u0 -> u0/2 przy (omega=1.1, b=8)
obserwacje poza kryteriami: b=6, omega=1.0, profil |delta_psi|
  za druga soczewka (superpozycja jakosciowo)
metrologia: jak proba #1 (D8): dryft sekularny energii <= 1e-4;
  gamma_core <= 0.01 (fit po tranzycie); wrap-guard 0.9L;
  determinizm bitowy w Phase 1; zapisy co 10 krokow
```

## 3. Scenariusze

1. `background_check + LIFETIME GATE (G2)`: relaksacja siatki;
   diagnostyka szwu, kret (n_tot=0, 4 rdzenie), residuum; potem
   ewolucja II rzedu BEZ pakietu, tau=400: pozycje rdzeni co 10
   krokow, max przemieszczenie, tempo ucieczki gamma_s (fit, jesli
   widoczne). TO JEST BRAMKA: FAIL -> STOP.
2. `dispersion_TV (G3)`: kalibracja jak proba #1 (tam: 1.21%).
3. `deflection(b, omega) (G4, G5)`: skan glowny; b_eff;
   alpha vs alpha_pred(b_eff).
4. `symmetry (G6a)`: b=-8 i tlo lustrzane.
5. `linearity (G6b)`: u0/2.
6. obserwacje poza kryteriami: b=6, omega=1.0, superpozycja.

## 4. Kryteria (propozycja do LOCKa)

- **G1 techniczne**: determinizm bitowy; brak NaN/S_GUARD; dryft
  sekularny energii <= 1e-4; gamma_core <= 0.01; wrap-guard;
  kret tla zachowany (4 rdzenie, n_tot = 0) w kazdym runie.
- **G2 BRAMKA — czas zycia tla pod dynamika pomiarowa** (NOWA,
  lekcja proby #1): ewolucja II rzedu bez pakietu, tau = 400;
  max przemieszczenie kazdego rdzenia <= 1.0 dla tau <= 300.
  FALSYFIKATOR: siodlo szachownicy ucieka z szumu numerycznego
  szybciej niz 3x przelot -> tlo niezdatne; STOP, deliverable =
  zmierzone czasy zycia (para vs siatka). Zadnych runow ugiecia
  przy G2 FAIL.
- **G3 kalibracja**: dyspersja kanalu amplitudowego na tle s*:
  blad <= 5% dla k in [0.5, 1.5].
- **G4 istnienie refrakcji**: |alpha(b=8, omega=1.1)| > 3*sigma
  ORAZ |alpha| maleje z b (omega=1.1) i z omega (b=8).
- **G5 RDZEN — policzalnosc (P1, geometria wielosoczewkowa)**:
  dla WSZYSTKICH szesciu punktow (b in {8,12,16} x omega in
  {1.1,1.3}): znak KU wirowi ORAZ alpha/alpha_pred(b_eff)
  in [0.5, 2.0]. FALSYFIKATOR: geometria propagacji nieczytelna
  z substratu w rezimie falowym.
- **G6 symetrie i liniowosc**: (a) |alpha(+8) - alpha(-8)| <=
  0.05 * |alpha(+8)| ORAZ wariant lustrzany zgodny w tej samej
  tolerancji (pre-rejestracja #1.4: oczekiwane roznice ~zaokraglen;
  wieksza roznica = blad implementacji, nie fizyka); (b) alpha(u0)
  vs alpha(u0/2): roznica <= 10%.

## 5. Regula decyzyjna (propozycja)

```text
G2 FAIL -> STOP przed pomiarem; deliverable: czasy zycia tla
  (siatka vs para); nastepny kandydat (osobny lock): L=256
  (dla pary: czas kolapsu ~liniowy w separacji) albo pomiar
  w ukladzie wspolporuszajacym sie.
G5 PASS (+ G4, G6 PASS)
  -> MOST GEOMETRYCZNY POTWIERDZONY ILOSCIOWO w geometrii
     wielosoczewkowej: policzalna, skupiajaca geometria propagacji
     wokol materii-defektow + superpozycja soczewek dziala.
     Nastepne opy (osobne locki): (i) pelny P2 na tle asymetrycznym
     (czesc cyrkulacyjna nieskasowana symetria), (ii) skalowanie
     z |n| (ugiecie ~ n^2?), (iii) granica newtonowska 2D i G_eff.
G5 FAIL (wielkosc poza pasmem przy dobrym znaku)
  -> eikonal nie opisuje rezimu falowego (lambda ~ b): raport
     rozjazdu; powtorka krotkofalowa (wyzsze omega / wieksze L)
     jako osobny op.
G5 FAIL (zly znak) -> falsyfikacja kanalu geometrycznego wokol
     defektow; program wraca do analizy.
G3/G4 FAIL -> stop przed interpretacja.
G6a FAIL -> blad implementacji: najpierw errata techniczna (dzien
     runow), dopiero potem ewaluacja G5.
```

## 6. Forbidden moves (propozycja)

1. Zero pinningu i zero wiezow: symetria szachownicy ma sie
   utrzymac SAMA (dynamicznie); zadnego przypinania, re-centrowania
   ani projekcji symetryzujacej po starcie.
2. BRAMKI G2 nie wolno oslabic, skrocic ani zreinterpretowac po
   pierwszym runie; przy G2 FAIL zadnych runow ugiecia "na probe".
3. Zadnych zmian progow/parametrow/procedur po pierwszym runie
   (wyjatki pre-rejestrowane: jedno obnizenie dt; jeden kontrolny
   run dx=0.25).
4. Punkty poza-kryterialne (b=6, omega=1.0) niedolaczalne do G5.
5. alpha_odd = 0 z symetrii NIE wolno raportowac jako "pomiar
   uniwersalnosci cyrkulacyjnej" (to wlasnosc konstrukcji tla);
   pelny P2 = osobny op na tle asymetrycznym.
6. Nie interpretowac wyniku jako GR/soczewkowania obserwacyjnego;
   tau != czas. Nie promowac do core/ bez osobnej autoryzacji.

## 7. Czego ten cykl NIE robi

- NIE testuje uniwersalnosci w znaku kretu (P2) nietrywialnie —
  symetria szachownicy kasuje czesc cyrkulacyjna z konstrukcji;
  G6a testuje wylacznie spojnosc implementacji z ta symetria.
- Nie mierzy skalowania z |n|; nie wyprowadza metryki; zero GR.
- Nie rozstrzyga kanalu dylatacyjnego (galaz C — osobna powtorka
  z procedura zachowujaca symetrie, wg zamkniecia
  op-scalar-sector-phistar).

## 8. Deliverables

- `Phase0_balance.md` (po autoryzacji: LOCKED + doprecyzowania
  pre-code w sekcji 10 PRZED kodem; wzorzec: sekcja 10 proby #1)
- `Phase0_check.py`/`Phase0_output.txt` (bilans sil na siatce
  z obrazami; n^2 na zrelaksowanym tle siatki; ray tracing
  alpha_pred(b_grid, omega_grid) + czulosc plaszczyzny ewaluacji;
  k0/v_g; budzet czasu)
- `Phase1_gate_and_calibration.py`/`Phase1_output.txt` (G2 BRAMKA
  + G3 + background_check + czesc G1)
- `Phase2_deflection.py`/`Phase2_output.txt` (G4, G5, G6b
  + obserwacje)
- `Phase3_symmetry.py`/`Phase3_output.txt` (G6a)
- `Phase_FINAL_close.md`, `README.md`

## 9. Status

**LOCKED (2026-07-05).** Model (tlo szachownicowe z dokladnym
kasowaniem sil z symetrii), predykcje (w tym pre-rejestracja
alpha_odd = 0 z konstrukcji i geometrii dwusoczewkowej), nowa
bramka G2 (czas zycia tla pod dynamika pomiarowa — bezposrednia
lekcja z proby #1), kryteria i regula decyzyjna spisane przed
kodem. Autoryzacja autora: polecenie realizacji cyklu
(2026-07-05). Wykonawca: osobny agent, nowa sesja (spelnione).
Obowiazuje pelny anty-Lakatos. Doprecyzowania pre-code: sekcja 10
(spisane PRZED pierwsza linijka kodu; wzorzec: sekcja 10 LOCKa
op-vortex-refraction, D1-D13).

## 10. Doprecyzowania pre-code (spisane przed pierwsza linijka kodu)

Zadne z ponizszych nie zmienia progow G1-G6, modelu ani predykcji —
to ujednoznacznienie procedur. Wiekszosc przeniesiona WPROST
z sekcji 10 proby #1 (op-vortex-refraction, D1-D13), z adaptacja
do tla siatkowego; D14-D16 nowe (detektor 2+2, procedura bramki G2,
obserwacja superpozycji):

- **D1 (integrator):** identycznie jak proba #1 — pelny sektor
  zespolony II rzedu, leapfrog pozycyjny `psi_{n+1} = 2 psi_n -
  psi_{n-1} + dt^2 * RHS(psi_n)`, `RHS(psi) = kappa*Lap(psi) -
  psi*(a - b|psi| + c|psi|^2)`; dt = 0.05. `psi_tau` do obserwabli
  = roznica CENTRALNA `(psi_{n+1} - psi_{n-1})/(2 dt)`.
- **D2 (pakiet, start analityczny):** `psi(0) = psi_bg + e^{i
  theta_bg} * u0 * env * cos(k0 (x - x0))`, `psi(-dt) = psi_bg +
  e^{i theta_bg} * u0 * env * cos(k0 (x - x0) + omega dt)` (pakiet
  prawobiezny); env = exp(-(x-x0)^2/(2 sx^2) - (y-y0)^2/(2 sy^2)),
  sx = sy = 8, u0 = 1e-3, x0 = -60, y0 = -32 + b;
  k0 = sqrt((omega^2 - V''(s*)) / kappa).
- **D3 (referencja blizniacza):** identycznie jak proba #1:
  `delta_psi(tau) = psi(tau) - psi_ref(tau)`, psi_ref ewoluuje TYM
  SAMYM integratorem z psi_bg bez pakietu (start: psi_ref(-dt) =
  psi_ref(0) = psi_bg), w jednej petli (lockstep). Pozycje rdzeni
  do b_eff, okna i gamma_core czytane z psi_ref (detektor D14).
  To NIE jest pinning (zero ingerencji w dynamike).
- **D4 (okno i wyzwolenie — adaptacja do siatki):** wyzwolenie
  CENTROIDEM (lekcja proby #1): x_c(|delta_psi|^2) > x1 + 10,
  gdzie x1 = BIEZACA (z psi_ref) pozycja x rdzenia +1 rzedu y=-32.
  W chwili wyzwolenia progi ZAMRAZANE: okno = {x1_frozen + 10 < x
  < x2_frozen - 10}, x2 = pozycja x rdzenia -1 rzedu y=-32
  (z psi_ref, w chwili wyzwolenia). Potem 5 probek co Delta_tau=4:
  `P = -Re sum(conj(delta_psi_tau) * grad(delta_psi)) dx^2`
  w zamrozonym oknie; usrednienie WEKTOROWE (Px, Py); `alpha :=
  -sign(b) * atan2(Py, Px)` [deg] (alpha > 0 = KU wirowi);
  sigma_alpha = odch. std probek katowych.
- **D5 (b_eff):** min po tau (probki co 10 krokow, do chwili
  wyzwolenia) odleglosci periodycznej centroidu |delta_psi|^2 od
  rdzenia rozpraszajacego (+1, nominal (-32,-32); pozycja z psi_ref).
  alpha_pred(b_eff) = ray tracing przy b = b_eff na tym samym
  zrelaksowanym tle siatki (bez interpolacji tabel).
- **D6 (eikonal na rzeczywistym tle siatki):** n^2(x,y) = (omega^2
  - V''(|psi_bg|)) / (omega^2 - V''(s*)) na PELNYM zrelaksowanym
  tle siatki (4 soczewki + obrazy periodyczne w polu), bilinearna
  interpolacja n^2 i gradientu; ray tracing hamiltonowski:
  dr/dsigma = p, dp/dsigma = grad(n^2)/2, |p(start)| = n(start),
  RK4, dsigma = 0.02; start x = -60, y = -32 + sign(b)*b_eff,
  kierunek +x. NOWE vs proba #1: kat promienia czytany
  w plaszczyznie x_eval = pozycja x centroidu pakietu w 3.
  (srodkowej) probce; CZULOSC raportowana: kat rowniez na krancach
  zamrozonego okna (x1_frozen+10 i x2_frozen-10). Przechwyt:
  r < 0.6 od ktoregokolwiek z 4 rdzeni -> NaN (raport).
- **D7 (G3, dyspersja kanalu amplitudowego):** identycznie jak
  proba #1: tlo jednorodne psi = s*; zaburzenie RZECZYWISTE
  u0*cos(k x), u_tau(0) = 0; k = 2 pi m / L, m in {11, 15, 20, 25,
  30} (k in [0.54, 1.47]); omega_meas ze zliczania przejsc przez
  zero Re(psi - s*) w punkcie probki, okno tau = 100; omega_teor =
  sqrt(kappa k^2 + V''(s*)); dyspersja sieciowa obok (poza
  kryterium).
- **D8 (G1, metryki):** jak proba #1, z rozszerzeniem na 4 rdzenie:
  dryft sekularny energii = |<H> ostatnia cwiartka - <H> pierwsza
  cwiartka| / |H(0)|, prog 1e-4; oscylacja leapfroga raportowana
  osobno. gamma_core = slope fitu ln RMS |delta_psi| w SUMIE okien
  r < 5 wokol WSZYSTKICH CZTERECH rdzeni (pozycje z psi_ref), fit
  po tranzycie (od wyzwolenia), prog 0.01. Kret tla: n_tot(psi_ref,
  koniec) = 0 ORAZ 4 rdzenie obecne (2 klastry na znak, D14).
  Wrap-guard: |przemieszczenie centroidu od startu| < 0.9 L.
  S_GUARD = 10 flaga (nie clamp). Determinizm bitowy w Phase 1
  (relaksacja + run falowy). Zapisy co 10 krokow; seed = 20260704.
- **D9 (budzet krokow):** relaksacja tla: 2000 krokow dt = 0.02.
  Bramka G2: 8000 krokow (tau = 400) bez pakietu. Runy falowe:
  limit 6000 krokow (tau = 300) dla omega = 1.1/1.3; 8000 (tau =
  400) dla obserwacji omega = 1.0; run konczy sie po 5 probkach
  (plus ewentualna kontynuacja obserwacyjna D16) albo na guardzie/
  limicie (brak 5 probek = raport, run niewazny do kryterium,
  przyczyna zapisana).
- **D10 (G6b, liniowosc):** |alpha(u0) - alpha(u0/2)| <= 0.1 *
  max(|alpha(u0)|, |alpha(u0/2)|) przy (omega = 1.1, b = 8).
- **D11 (PRE-REJESTRACJA: warianty symetryczne G6a):** (a) tlo
  lustrzane = negacja fazy ansatzu (theta -> -theta); wszystkie
  wspolczynniki rzeczywiste, pakiet w lokalnej ramce fazy SWOJEGO
  tla -> run lustrzany jest DOKLADNIE sprzezeniem zespolonym runu
  zwyklego -> oczekiwanie: alpha identyczne do zaokraglen (test
  IMPLEMENTACJI). (b) b = -8: odbicie y -> -64 - y (wzgledem osi
  wiazki y = -32) zlozone ze sprzezeniem odwzorowuje tlo siatki
  na siebie (siatka i grid dokladnie symetryczne) -> oczekiwanie:
  alpha(-8) = alpha(+8) do zaokraglen. Obie rownosci = G6a;
  wieksza roznica = blad implementacji (sekcja 1.4), nie fizyka;
  alpha_odd NIE jest raportowana jako pomiar P2 (forbidden move #5).
- **D12 (uwaga dyfrakcyjna, powtorzona):** lambda(1.1) = 7.72,
  lambda(1.3) = 4.93 vs b in {8, 12, 16} — rezim falowy; pasmo G5
  [0.5, 2.0] szerokie z tego powodu; zadnych poprawek
  dyfrakcyjnych po fakcie.
- **D13 (narzedzia):** zero `cd` w powloce; sciezki od korzenia
  vaulta; kazdy skrypt pisze output do `PhaseX_output.txt`
  w katalogu opa.
- **D14 (detektor rdzeni 2+2, NOWY):** winding plakietowy W;
  komorki z W*sign > 0.5 klastrowane po odleglosci periodycznej
  (promien <= 8 komorek od komorki bazowej klastra); wymagane
  DOKLADNIE 2 klastry na znak (inaczej flaga core_lost); pozycja
  klastra = centroid deficytu Phi (wagi max(0, Phi*/2 - Phi))
  w oknie 13x13 wokol komorki bazowej; przypisanie klastrow do
  nominalow (+1: (-32,-32), (+32,+32); -1: (+32,-32), (-32,+32))
  po najmniejszej sumie odleglosci periodycznych.
- **D15 (bramka G2 — procedura, NOWA):** start z zrelaksowanego
  tla, psi(-dt) = psi(0) (zerowa predkosc), dynamika II rzedu D1,
  8000 krokow (tau = 400) BEZ pakietu; pozycje 4 rdzeni (D14) co
  10 krokow; kryterium: max po rdzeniach i po tau <= 300
  odleglosci periodycznej od pozycji poczatkowej <= 1.0;
  przedzial tau in (300, 400] raportowany poza kryterium (zapas).
  gamma_s: fit liniowy ln(max_disp) vs tau na odcinku max_disp in
  [0.05, 1.0] (jesli osiagalny); przy braku wzrostu raport "brak
  widocznej ucieczki w oknie". Bramki nie wolno oslabic, skrocic
  ani zreinterpretowac po pierwszym runie; G2 FAIL -> STOP, zero
  runow ugiecia.
- **D16 (obserwacja superpozycji, poza kryteriami):** run (omega =
  1.1, b = 8) po zebraniu 5 probek kontynuowany do x_c > x2_frozen
  + 10 albo limitu D9; raport jakosciowy: centroid y, szerokosc
  poprzeczna pakietu i kat wektorowy P w oknie x > x2_frozen + 10
  (za druga soczewka rzedu). Bez statusu kryterialnego.
