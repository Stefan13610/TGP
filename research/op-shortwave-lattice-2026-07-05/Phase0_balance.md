---
title: "Phase0_balance — DRAFT: B2-prime proba #3 — refrakcja krotkofalowa na siatce szachownicowej L=256 (kazdy punkt w domenie eikonalu, b_eff >= 1.5 lambda)"
date: 2026-07-05
type: phase0-draft
tgp_owner: research/op-shortwave-lattice-2026-07-05
status: LOCKED
anti_lakatos_lock: ACTIVE (autoryzacja autora: polecenie realizacji cyklu, 2026-07-05; lock przed kodem; wykonawca: osobny agent, nowa sesja — spelnione)
tags: [gravity-bridge, branch-B2prime, lattice-background, vortex-lensing, refraction, eikonal, shortwave, L256]
related:
  - "[[../op-lattice-background-2026-07-05/Phase_FINAL_close.md]]"
  - "[[../op-vortex-refraction-2026-07-05/Phase_FINAL_close.md]]"
  - "[[../op-stationary-background-2026-07-05/Phase_FINAL_close.md]]"
---

# Phase 0 — DRAFT: B2-prime proba #3 — refrakcja krotkofalowa (L=256)

## 0. Hipoteza robocza i miejsce w programie

Stan lancucha (wszystko zmierzone, nic do wiary):
- proba #1 (op-vortex-refraction): tlo pary anihiluje (tau=107.5)
  — 0/9 runow waznych; aparat pomiarowy sprawny.
- proba #2 (op-lattice-background): tlo szachownicowe ROZWIAZALO
  problem tla (bramka G2: przemieszczenie rdzeni 0.0000 przez
  tau<=400); 9/9 runow waznych; znak KU wirowi WSZEDZIE; symetrie
  dokladne (G6a do 4e-04 deg / bitowo); liniowosc dokladna.
  **G5 FAIL formalnie: 5/6 punktow w pasmie [0.5, 2.0]**; wypadl
  wylacznie punkt najglebiej dyfrakcyjny (omega=1.1, b=8):
  b_eff = 7.17 ~ lambda = 7.72, ratio 0.445 przy dobrym znaku —
  dokladnie rezim zastrzezony z gory (D12 proby #2).
- wniosek pomiarowy #2: granica stosowalnosci eikonalu lezy przy
  b_eff ~ lambda; dla b_eff >= 1.5 lambda ratio = 0.74-1.46.

Ten cykl realizuje sciezke reguly decyzyjnej proby #2 ("G5 FAIL
wielkosc przy dobrym znaku -> powtorka krotkofalowa jako osobny
op"): **te same pytanie P1, ale KAZDY punkt kryterialny
zaprojektowany w domenie stosowalnosci eikonalu** (b_nom >=
1.5 lambda + 0.9 marginesu na dip b_eff), dwiema dzwigniami:
1. **wyzsze omega**: dodana seria omega = 1.5 (lambda = 3.79),
2. **L = 256** (N = 512): wieksze b dostepne dla omega = 1.1 bez
   nakladania soczewek (rozstaw siatki L/2 = 128); przy okazji
   czystszy test wielosoczewkowy (slabsze ogony sasiadow).

Pytanie:

> (P1-domkniecie) POLICZALNOSC W DOMENIE: czy kat ugiecia pakietu
>      amplitudowego zgadza sie z eikonalem substratu (znak +
>      wielkosc, pasmo [0.5, 2.0]) we WSZYSTKICH punktach
>      zaprojektowanych tak, by b_eff >= 1.5 lambda — na stabilnym
>      tle szachownicowym L=256?

G5 PASS tutaj + 5/6 z proby #2 = P1 rozstrzygniete: policzalna,
uniwersalnie skupiajaca geometria propagacji wokol materii-defektow
w calej domenie stosowalnosci eikonalu. G5 FAIL przy dobrym znaku
= eikonal nieczytelny NAWET w swojej domenie -> program wraca do
analizy (bez kolejnej powtorki bez nowej idei).

UWAGA spisana z gory: punkt (1.1, b=8) — ten, ktory oblal G5
w probie #2 — wraca tu WYLACZNIE jako obserwacja poza kryteriami
(test: czy ratio ~0.45 jest wlasnoscia rezimu falowego, a nie
geometrii L=128).

## 1. Predykcje algebraiczne do zalockowania (PRZED kodem)

### 1.1 Rownanie #1 — wspolczynnik zalamania (bez zmian)

```text
n^2(x,y) = (omega^2 - V''(|psi_bg|)) / (omega^2 - V''(s*))
ogon defektu: n^2 - 1 -> A_t / (r^2 (omega^2 - V''(s*))),
  A_t = V'''(s*) * s* * xi^2 = 2.5690
kanal amplitudowy: omega > m_TV = 0.9374;
k0(omega) = sqrt((omega^2 - V''(s*)) / kappa);  v_g = kappa*k0/omega
k0/lambda/v_g:  omega=1.1: 0.8140 / 7.72 / 0.3700
                omega=1.3: 1.2738 / 4.93 / 0.4899   (obserwacja)
                omega=1.5: 1.6561 / 3.79 / 0.5520
                omega=1.7: 2.0057 / 3.13 / 0.5899   (obserwacja)
```

### 1.2 Rownanie #2 — projekt punktow kryterialnych (fundament cyklu)

```text
kryterium projektowe: b_nom >= 1.5*lambda + 0.9
  (0.9 = najwiekszy zmierzony dip b_eff w probie #2: b_nom - b_eff
   in [0.32, 0.84]; margines pokrywa go z zapasem)
seria omega=1.1 (lambda=7.72; 1.5*lambda+0.9 = 12.5):
  b in {14, 16, 20}   -> b_eff/lambda oczekiwane ~1.70/1.96/2.48
seria omega=1.5 (lambda=3.79; 1.5*lambda+0.9 = 6.6):
  b in {8, 12, 16}    -> b_eff/lambda oczekiwane ~1.9/3.0/4.0
zestaw kryterialny G5 = te 6 punktow, ZAMKNIETY;
b_eff/lambda raportowane per run; ZADNEGO wykluczania punktow
post hoc, nawet gdyby dip okazal sie wiekszy niz margines.
```

### 1.3 Rownanie #3 — zaklad eikonalny

Orientacyjnie (tlo siatki L=128, proba #2 Phase0, srodek okna;
alpha > 0 = KU wirowi): omega=1.1: b=12: +5.45, b=16: +2.90;
omega=1.3: b=8: +5.41, b=12: +2.17, b=16: +1.17. Ekstrapolacje
orientacyjne dla tego cyklu: omega=1.1: b=14 ~ +3.9, b=20 ~ +1.8;
omega=1.5 (skalowanie ogona ~1/(omega^2 - V''(s*)), czynnik
0.8113/1.3713 = 0.59 wzgledem omega=1.3): b=8 ~ +3.2, b=12 ~ +1.3,
b=16 ~ +0.7. Zaklad WIAZACY liczy Phase0_check TEGO cyklu: ray
tracing na rzeczywistym, zrelaksowanym tle siatki L=256, z ewaluacja
kata w plaszczyznie pomiaru i czuloscia na krance okna (wzorzec:
D6 proby #2). Roznice vs tabela L=128 raportowane (slabsze ogony
sasiadow przy rozstawie 128).

### 1.4 Dyspersja sieciowa przy wysokim k (pre-rejestracja ryzyka)

```text
k0(1.5) = 1.656 lezy POZA pasmem kalibracji proby #2 ([0.5, 1.5])
-> pasmo G3 ROZSZERZONE do [0.5, 1.7]; oszacowanie analityczne
   bledu sieciowego (dx=0.5): k=1.66: ~-1.9%; k=2.01: ~-2.9%
   (omega_lat = sqrt(kappa*(2-2cos(k dx))/dx^2 + V''(s*)))
   — wewnatrz progu 5%, ale MIERZONE, nie zakladane (G3).
omega=1.7 (k0=2.006, k0*dx ~ 1.0) POZA kryteriami (obserwacja).
```

### 1.5 Symetrie tla (jak proba #2, pre-rejestracja G6a)

```text
S1 (odbicie osi wiazki y -> -128-y wzgledem y=-64) + sprzezenie:
   alpha(+b) = alpha(-b) DOKLADNIE (do zaokraglen)
S2 (lustro ladunkowe theta -> -theta): run = sprzezenie zespolone
   -> bitowa rownosc (proba #2: potwierdzone co do bitu)
alpha_odd = 0 z konstrukcji — NIE raportowac jako pomiar P2.
```

### 1.6 Budzet czasu (KRYTERIUM, wzorzec proby #2)

```text
geometria: rdzen rozpraszajacy +1 w (-64,-64); start x0 = -100
  (4.5 sx przed rdzeniem); prog wyzwolenia ~x = -54
przelot do 5. probki: omega=1.1: ~tau=140; omega=1.5: ~tau=99
bramka G2: tlo bez pakietu, dynamika II rzedu, okno tau = 600
  (12000 krokow); kryterium: max przemieszczenie kazdego rdzenia
  <= 1.0 na tau in [0, 450] (>= 3x przelot 140; pozycje z pola
  co 10 krokow); przedzial (450, 600] raportowany poza kryterium
kontrast (zmierzone): siatka L=128 — 0.0000 przez tau=400;
  para L=128 — przemieszczenie 14 przy tau=50, anihilacja 107.5
```

## 2. Model i parametry (propozycja do LOCKa)

```text
substrat: bez zmian (kappa=0.5, a=0.5, b=1.6, c=1.0; psi in C;
  eps=0.30, S_GUARD=10 flaga, seed=20260704)
  NOWE: N=512, L=256, dx=0.5 (bez zmian dx!), periodycznie
tlo: SZACHOWNICA 4 wirow w rozstawie L/2 = 128:
  +1 w (-64,-64)   -1 w (+64,-64)
  -1 w (-64,+64)   +1 w (+64,+64)
  ansatz: suma DWOCH PAR WSPOLLINIOWYCH POZIOMYCH (konstrukcja
  dokladna, wzorzec proby #2: rzad y=-64: theta_h_pair((-64,-64),
  (+64,-64)); rzad y=+64: theta_h_pair((+64,+64),(-64,+64));
  obrazy k in [-2,2] przy L=256); amplituda: iloczyn 4 profili
  rdzeni; diagnostyka szwu obowiazkowa (max skok fazy poza
  rdzeniami << pi); relaksacja 2000 krokow (dt=0.02); residuum
  raportowane (proba #2 przy L=128: 1.6e-06)
sektor falowy: PELNY zespolony II rzedu (bez zmian):
  d^2psi/dtau^2 = kappa*Lap(psi) - psi*(a - b|psi| + c|psi|^2)
  leapfrog pozycyjny dt=0.05; psi_tau = roznica centralna
referencja blizniacza: psi_ref bez pakietu w lockstep (bez zmian)
pakiet: delta_psi(0) = e^{i theta_bg} * u0 * env * cos(k0 (x - x0));
  start analityczny; u0 = 1e-3, sx = sy = 8, x0 = -100,
  y0 = -64 + b; wir rozpraszajacy: +1 w (-64,-64)
okno pomiaru (bez zmian koncepcyjnych): wyzwolenie CENTROIDEM
  (x_c > x1+10, x1 = biezaca pozycja x rdzenia +1 rzedu y=-64
  z psi_ref); w chwili wyzwolenia progi ZAMRAZANE:
  okno = {x1+10 < x < x2-10} (x2 = rdzen -1 rzedu y=-64);
  5 probek co Delta_tau = 4; usrednienie WEKTOROWE (Px, Py);
  alpha := -sign(b) * atan2(Py, Px) [deg] (alpha > 0 = KU wirowi)
b_eff: min odleglosc periodyczna centroidu |delta_psi|^2 od rdzenia
  rozpraszajacego (psi_ref), do chwili wyzwolenia; b_eff/lambda
  raportowane per run
alpha_pred(b_eff): ray tracing na TYM SAMYM zrelaksowanym tle
  siatki L=256, przy b = b_eff; kat w plaszczyznie x centroidu
  3. probki; czulosc na krancach okna raportowana
skan glowny (KRYTERIALNY, zamkniety): omega=1.1 x b {14, 16, 20};
  omega=1.5 x b {8, 12, 16}
symetrie (G6a): (omega=1.1, b=14) w wariantach {b=-14} i {tlo
  lustrzane, b=+14}
liniowosc (G6b): u0 -> u0/2 przy (omega=1.1, b=14)
obserwacje POZA kryteriami:
  (1.1, b=8)  — punkt oblany w probie #2: czy ratio ~0.45
                przenosi sie na L=256 (wlasnosc rezimu, nie siatki)
  (1.3, b=8)  — sprzeg z proba #2 (tam ratio 0.741 przy L=128)
  (1.7, b=8)  — ultra-krotkofalowo (dyspersja sieciowa ~3%)
  superpozycja: (1.5, b=8) kontynuowany za druga soczewke,
                dedykowany run z limitem 10000 krokow (tau=500)
metrologia: jak proba #2 (D8): dryft sekularny energii <= 1e-4;
  gamma_core <= 0.01 (fit po tranzycie, okna r<5 wokol 4 rdzeni);
  wrap-guard 0.9L = 230; determinizm bitowy w Phase 1; zapisy co
  10 krokow; detektor rdzeni 2+2 z klastrowaniem (D14 proby #2,
  nominaly przeskalowane do +/-64)
budzet krokow: relaksacja 2000 (dt=0.02); bramka G2: 12000
  (tau=600); runy falowe: 6000 (tau=300); wyjatek pre-rejestrowany:
  run obserwacyjny superpozycji 10000 (tau=500); brak 5 probek =
  run niewazny, przyczyna zapisana
koszt (informacyjnie): N=512 -> ~4x koszt kroku proby #2;
  bramka ~15 min, kazdy run ugiecia ~5-8 min, calosc ~1.5-2 h
```

## 3. Scenariusze

1. `background_check + LIFETIME GATE (G2)`: relaksacja siatki
   L=256; szew, kret (n_tot=0, 4 rdzenie), residuum; ewolucja II
   rzedu BEZ pakietu, tau=600: pozycje rdzeni co 10 krokow, max
   przemieszczenie, gamma_s (fit, jesli widoczny). BRAMKA:
   FAIL -> STOP.
2. `dispersion_TV (G3)`: kalibracja w pasmie ROZSZERZONYM
   k in [0.5, 1.7]: m in {22, 30, 41, 55, 68} przy L=256
   (k = {0.540, 0.736, 1.006, 1.350, 1.669}).
3. `deflection(b, omega) (G4, G5)`: skan glowny 6 punktow; b_eff;
   b_eff/lambda; alpha vs alpha_pred(b_eff).
4. `symmetry (G6a)`: b=-14 i tlo lustrzane (+14).
5. `linearity (G6b)`: u0/2 przy (1.1, 14).
6. obserwacje poza kryteriami: (1.1,8), (1.3,8), (1.7,8),
   superpozycja (1.5,8).

## 4. Kryteria (propozycja do LOCKa)

- **G1 techniczne**: determinizm bitowy; brak NaN/S_GUARD; dryft
  sekularny energii <= 1e-4; gamma_core <= 0.01; wrap-guard;
  kret tla zachowany (4 rdzenie, n_tot = 0) w kazdym runie.
- **G2 BRAMKA — czas zycia tla pod dynamika pomiarowa**: ewolucja
  II rzedu bez pakietu, tau = 600; max przemieszczenie kazdego
  rdzenia <= 1.0 dla tau <= 450 (>= 3x przelot 140).
  FALSYFIKATOR: siodlo szachownicy przy L=256 ucieka szybciej niz
  3x przelot -> tlo niezdatne; STOP, deliverable = czasy zycia
  (L=256 vs L=128 vs para). Zadnych runow ugiecia przy G2 FAIL;
  bramki nie wolno oslabic po pierwszym runie.
- **G3 kalibracja (pasmo rozszerzone)**: dyspersja kanalu
  amplitudowego na tle s*: blad <= 5% dla k in [0.5, 1.7].
  FALSYFIKATOR: siatka dx=0.5 nie niesie rzetelnie k0(1.5) ->
  stop przed interpretacja serii 1.5 (seria 1.1 ocenialna osobno
  TYLKO jesli tak rozstrzygnie errata w dniu runow — domyslnie
  stop calego G5).
- **G4 istnienie refrakcji**: |alpha(b=14, omega=1.1)| > 3*sigma
  ORAZ |alpha| maleje z b w serii omega=1.1 (14 > 16 > 20 co do
  |alpha|) ORAZ |alpha(1.1, 16)| > |alpha(1.5, 16)| (maleje
  z omega przy wspolnym b=16).
- **G5 RDZEN — policzalnosc w domenie eikonalu (P1)**: dla
  WSZYSTKICH SZESCIU punktow (1.1 x {14,16,20}; 1.5 x {8,12,16}):
  znak KU wirowi ORAZ alpha/alpha_pred(b_eff) in [0.5, 2.0].
  FALSYFIKATOR: eikonal substratu nie opisuje ugiecia nawet
  w swojej domenie stosowalnosci.
- **G6 symetrie i liniowosc**: (a) |alpha(+14) - alpha(-14)| <=
  0.05 * |alpha(+14)| ORAZ wariant lustrzany zgodny w tej samej
  tolerancji (pre-rejestracja #1.5); (b) alpha(u0) vs alpha(u0/2):
  roznica <= 10%.

## 5. Regula decyzyjna (propozycja)

```text
G2 FAIL -> STOP przed pomiarem; deliverable: czasy zycia
  (L=256 vs L=128); kandydat: pomiar w ukladzie wspolporuszajacym.
G5 PASS (+ G4, G6 PASS)
  -> P1 ROZSTRZYGNIETE POZYTYWNIE (lacznie z 5/6 proby #2):
     policzalna, uniwersalnie skupiajaca geometria propagacji
     wokol materii-defektow, z jawnie zmierzona granica
     stosowalnosci eikonalu (b_eff ~ lambda). Nastepne opy
     (osobne locki): (i) pelny P2 na tle asymetrycznym,
     (ii) skalowanie z |n| (ugiecie ~ n^2?), (iii) granica
     newtonowska 2D i G_eff z parametrow substratu.
G5 FAIL (wielkosc poza pasmem przy dobrym znaku)
  -> eikonal nieczytelny w swojej domenie: powazny problem
     policzalnosci kanalu; program WRACA DO ANALIZY (zadnej
     kolejnej powtorki bez nowej idei teoretycznej); raport
     pelnego rozjazdu wszystkich 12 zmierzonych punktow (#2 + #3).
G5 FAIL (zly znak w ktorymkolwiek punkcie) -> falsyfikacja kanalu
     geometrycznego wokol defektow.
G3 FAIL -> stop przed interpretacja (patrz G3: domyslnie stop
     calego G5, errata w dniu runow moze zawezic do serii 1.5 —
     wylacznie na podstawie tego, KTORE k oblalo).
G4 FAIL -> stop przed interpretacja.
G6a FAIL -> blad implementacji: errata techniczna (dzien runow),
     dopiero potem ewaluacja G5.
```

## 6. Forbidden moves (propozycja)

1. Zero pinningu i zero wiezow; symetria szachownicy utrzymuje sie
   SAMA (zwalidowane w probie #2 przy L=128; tu MIERZONE ponownie).
2. BRAMKI G2 nie wolno oslabic, skrocic ani zreinterpretowac po
   pierwszym runie; przy G2 FAIL zadnych runow ugiecia.
3. Zadnych zmian progow/parametrow/procedur po pierwszym runie
   (wyjatki pre-rejestrowane: jedno obnizenie dt; jeden kontrolny
   run dx=0.25).
4. Zestaw kryterialny G5 (6 punktow) ZAMKNIETY: zadnego dolaczania
   obserwacji (1.1,8), (1.3,8), (1.7,8) ani wykluczania punktow
   po fakcie (rowniez gdyby b_eff/lambda wypadl ponizej 1.5 —
   raport, ale punkt zostaje kryterialny).
5. alpha_odd = 0 z symetrii NIE raportowac jako pomiaru P2.
6. Nie interpretowac wyniku jako GR/soczewkowania obserwacyjnego;
   tau != czas. Nie promowac do core/ bez osobnej autoryzacji.
7. Limit 10000 krokow dotyczy WYLACZNIE dedykowanego runu
   obserwacyjnego superpozycji (nie wolno nim ratowac runow
   kryterialnych bez 5 probek).

## 7. Czego ten cykl NIE robi

- NIE testuje uniwersalnosci w znaku kretu (P2) nietrywialnie —
  jak proba #2 (symetria kasuje czesc cyrkulacyjna z konstrukcji);
  pelny P2 = osobny op na tle asymetrycznym.
- NIE bada rezimu dyfrakcyjnego b_eff < 1.5 lambda kryterialnie
  (dwie obserwacje: 1.1/8, 1.3/8 — wylacznie raport).
- Nie mierzy skalowania z |n|; nie wyprowadza metryki; zero GR.

## 8. Deliverables

- `Phase0_balance.md` (po autoryzacji: LOCKED + doprecyzowania
  pre-code w sekcji 10 PRZED kodem; wzorzec: sekcja 10 proby #2,
  D1-D16 — wiekszosc przenosi sie wprost po przeskalowaniu
  geometrii)
- `Phase0_check.py`/`Phase0_output.txt` (bilans sil na siatce
  L=256 — powloki D4 CENTROWANE NA WIRZE, wzorzec erraty #1
  proby #2; dyspersja sieciowa analitycznie przy k0(1.5)/k0(1.7);
  n^2 na zrelaksowanym tle; ray tracing WIAZACY alpha_pred(b, omega)
  + czulosc plaszczyzn; tabela projektowa b_eff/lambda; budzet
  czasu i kosztu)
- `Phase1_gate_and_calibration.py`/`Phase1_output.txt` (G2 BRAMKA
  tau=600 + G3 pasmo [0.5, 1.7] + background_check + czesc G1)
- `Phase2_deflection.py`/`Phase2_output.txt` (G4, G5, G6b
  + obserwacje 1.1/8, 1.3/8, 1.7/8, superpozycja)
- `Phase3_symmetry.py`/`Phase3_output.txt` (G6a)
- `Phase_FINAL_close.md`, `README.md`

## 9. Status

**LOCKED (2026-07-05).** Zestaw kryterialny zaprojektowany
W CALOSCI w domenie stosowalnosci eikonalu (b_nom >= 1.5 lambda
+ 0.9; margines z pomiaru dipu b_eff w probie #2), pasmo
kalibracji rozszerzone do k0(1.5) z pre-rejestrowanym oszacowaniem
bledu sieciowego, bramka G2 przeskalowana do nowego przelotu
(tau<=450, okno 600), obserwacje kontrolne laczace proby #2 i #3
spisane przed kodem. Autoryzacja autora: polecenie realizacji
cyklu (2026-07-05). Wykonawca: osobny agent, nowa sesja
(spelnione). Obowiazuje pelny anty-Lakatos. Doprecyzowania
pre-code: sekcja 10 (spisane PRZED pierwsza linijka kodu; wzorzec:
sekcja 10 LOCKa proby #2, D1-D16; zmiany czysto geometryczne).

## 10. Doprecyzowania pre-code (spisane przed pierwsza linijka kodu)

Zadne z ponizszych nie zmienia progow G1-G6, modelu ani predykcji —
to ujednoznacznienie procedur. Przeniesione WPROST z sekcji 10
proby #2 (op-lattice-background, D1-D16) z adaptacja czysto
geometryczna do N=512, L=256 (nominaly rdzeni +/-64, x0=-100,
prog wyzwolenia ~-54, m-lista G3, limity krokow wg sekcji 2):

- **D1 (integrator):** identycznie jak proba #2 — pelny sektor
  zespolony II rzedu, leapfrog pozycyjny `psi_{n+1} = 2 psi_n -
  psi_{n-1} + dt^2 * RHS(psi_n)`, `RHS(psi) = kappa*Lap(psi) -
  psi*(a - b|psi| + c|psi|^2)`; dt = 0.05. `psi_tau` do obserwabli
  = roznica CENTRALNA `(psi_{n+1} - psi_{n-1})/(2 dt)`.
- **D2 (pakiet, start analityczny):** `psi(0) = psi_bg + e^{i
  theta_bg} * u0 * env * cos(k0 (x - x0))`, `psi(-dt) = psi_bg +
  e^{i theta_bg} * u0 * env * cos(k0 (x - x0) + omega dt)` (pakiet
  prawobiezny); env = exp(-(x-x0)^2/(2 sx^2) - (y-y0)^2/(2 sy^2)),
  sx = sy = 8, u0 = 1e-3, x0 = -100, y0 = -64 + b;
  k0 = sqrt((omega^2 - V''(s*)) / kappa).
- **D3 (referencja blizniacza):** identycznie jak proba #2:
  `delta_psi(tau) = psi(tau) - psi_ref(tau)`, psi_ref ewoluuje TYM
  SAMYM integratorem z psi_bg bez pakietu (start: psi_ref(-dt) =
  psi_ref(0) = psi_bg), w jednej petli (lockstep). Pozycje rdzeni
  do b_eff, okna i gamma_core czytane z psi_ref (detektor D14).
  To NIE jest pinning (zero ingerencji w dynamike).
- **D4 (okno i wyzwolenie):** wyzwolenie CENTROIDEM:
  x_c(|delta_psi|^2) > x1 + 10, gdzie x1 = BIEZACA (z psi_ref)
  pozycja x rdzenia +1 rzedu y=-64. W chwili wyzwolenia progi
  ZAMRAZANE: okno = {x1_frozen + 10 < x < x2_frozen - 10}, x2 =
  pozycja x rdzenia -1 rzedu y=-64 (z psi_ref, w chwili
  wyzwolenia). Potem 5 probek co Delta_tau = 4: `P = -Re
  sum(conj(delta_psi_tau) * grad(delta_psi)) dx^2` w zamrozonym
  oknie; usrednienie WEKTOROWE (Px, Py); `alpha := -sign(b) *
  atan2(Py, Px)` [deg] (alpha > 0 = KU wirowi); sigma_alpha =
  odch. std probek katowych.
- **D5 (b_eff):** min po tau (probki co 10 krokow, do chwili
  wyzwolenia) odleglosci periodycznej centroidu |delta_psi|^2 od
  rdzenia rozpraszajacego (+1, nominal (-64,-64); pozycja
  z psi_ref). alpha_pred(b_eff) = ray tracing przy b = b_eff na
  tym samym zrelaksowanym tle siatki L=256 (bez interpolacji
  tabel). b_eff/lambda raportowane per run (sekcja 1.2); zadnego
  wykluczania punktow post hoc.
- **D6 (eikonal na rzeczywistym tle siatki):** n^2(x,y) = (omega^2
  - V''(|psi_bg|)) / (omega^2 - V''(s*)) na PELNYM zrelaksowanym
  tle siatki L=256 (4 soczewki + obrazy periodyczne w polu),
  bilinearna interpolacja n^2 i gradientu; ray tracing
  hamiltonowski: dr/dsigma = p, dp/dsigma = grad(n^2)/2,
  |p(start)| = n(start), RK4, dsigma = 0.02; start x = -100,
  y = -64 + sign(b)*b_eff, kierunek +x. Kat promienia czytany
  w plaszczyznie x_eval = pozycja x centroidu pakietu w 3.
  (srodkowej) probce; CZULOSC raportowana: kat rowniez na krancach
  zamrozonego okna (x1_frozen+10 i x2_frozen-10). Przechwyt:
  r < 0.6 od ktoregokolwiek z 4 rdzeni -> NaN (raport).
- **D7 (G3, dyspersja kanalu amplitudowego, pasmo ROZSZERZONE):**
  tlo jednorodne psi = s*; zaburzenie RZECZYWISTE u0*cos(k x),
  u_tau(0) = 0; k = 2 pi m / L, m in {22, 30, 41, 55, 68}
  (k = {0.540, 0.736, 1.006, 1.350, 1.669} — pasmo [0.5, 1.7]);
  omega_meas ze zliczania przejsc przez zero Re(psi - s*)
  w punkcie probki, okno tau = 100; omega_teor = sqrt(kappa k^2
  + V''(s*)); dyspersja sieciowa analityczna raportowana obok
  (pre-rejestracja 1.4: k=1.66: ~-1.9%, k=2.01: ~-2.9%).
- **D8 (G1, metryki):** jak proba #2: dryft sekularny energii =
  |<H> ostatnia cwiartka - <H> pierwsza cwiartka| / |H(0)|, prog
  1e-4; oscylacja leapfroga raportowana osobno. gamma_core =
  slope fitu ln RMS |delta_psi| w SUMIE okien r < 5 wokol
  WSZYSTKICH CZTERECH rdzeni (pozycje z psi_ref), fit po tranzycie
  (od wyzwolenia), prog 0.01. Kret tla: n_tot(psi_ref, koniec) = 0
  ORAZ 4 rdzenie obecne (2 klastry na znak, D14). Wrap-guard:
  |przemieszczenie centroidu od startu| < 0.9 L = 230.4. S_GUARD =
  10 flaga (nie clamp). Determinizm bitowy w Phase 1 (relaksacja
  + run falowy). Zapisy co 10 krokow; seed = 20260704.
- **D9 (budzet krokow):** relaksacja tla: 2000 krokow dt = 0.02.
  Bramka G2: 12000 krokow (tau = 600) bez pakietu. Runy falowe:
  limit 6000 krokow (tau = 300). Wyjatek pre-rejestrowany
  (sekcja 2): DEDYKOWANY run obserwacyjny superpozycji (1.5, b=8)
  z limitem 10000 krokow (tau = 500). Run konczy sie po 5 probkach
  (plus kontynuacja obserwacyjna D16 w runie dedykowanym) albo na
  guardzie/limicie (brak 5 probek = raport, run niewazny do
  kryterium, przyczyna zapisana; limitu 10000 NIE wolno uzyc do
  ratowania runow kryterialnych — forbidden move #7).
- **D10 (G6b, liniowosc):** |alpha(u0) - alpha(u0/2)| <= 0.1 *
  max(|alpha(u0)|, |alpha(u0/2)|) przy (omega = 1.1, b = 14).
- **D11 (PRE-REJESTRACJA: warianty symetryczne G6a):** (a) tlo
  lustrzane = negacja fazy ansatzu (theta -> -theta); wszystkie
  wspolczynniki rzeczywiste, pakiet w lokalnej ramce fazy SWOJEGO
  tla -> run lustrzany jest DOKLADNIE sprzezeniem zespolonym runu
  zwyklego -> oczekiwanie: alpha identyczne do zaokraglen (test
  IMPLEMENTACJI). (b) b = -14: odbicie y -> -128 - y (wzgledem osi
  wiazki y = -64) zlozone ze sprzezeniem odwzorowuje tlo siatki
  na siebie (siatka i grid dokladnie symetryczne) -> oczekiwanie:
  alpha(-14) = alpha(+14) do zaokraglen. Obie rownosci = G6a;
  wieksza roznica = blad implementacji (sekcja 1.5), nie fizyka;
  alpha_odd NIE jest raportowana jako pomiar P2 (forbidden
  move #5).
- **D12 (domena eikonalu, powtorzona):** zestaw kryterialny G5
  zaprojektowany w calosci w domenie b_nom >= 1.5 lambda + 0.9
  (sekcja 1.2); obserwacje (1.1,8), (1.3,8), (1.7,8) POZA
  kryteriami (rezim dyfrakcyjny / ultra-krotkofalowy); zadnych
  poprawek dyfrakcyjnych po fakcie; zadnego wykluczania punktow
  kryterialnych, nawet gdyby b_eff/lambda < 1.5.
- **D13 (narzedzia):** zero `cd` w powloce; sciezki od korzenia
  vaulta; kazdy skrypt pisze output do `PhaseX_output.txt`
  w katalogu opa.
- **D14 (detektor rdzeni 2+2):** winding plakietowy W; komorki
  z W*sign > 0.5 klastrowane po odleglosci periodycznej (promien
  <= 8 komorek od komorki bazowej klastra); wymagane DOKLADNIE
  2 klastry na znak (inaczej flaga core_lost); pozycja klastra =
  centroid deficytu Phi (wagi max(0, Phi*/2 - Phi)) w oknie 13x13
  wokol komorki bazowej; przypisanie klastrow do nominalow
  (+1: (-64,-64), (+64,+64); -1: (+64,-64), (-64,+64)) po
  najmniejszej sumie odleglosci periodycznych.
- **D15 (bramka G2 — procedura):** start z zrelaksowanego tla,
  psi(-dt) = psi(0) (zerowa predkosc), dynamika II rzedu D1,
  12000 krokow (tau = 600) BEZ pakietu; pozycje 4 rdzeni (D14) co
  10 krokow; kryterium: max po rdzeniach i po tau <= 450
  odleglosci periodycznej od pozycji poczatkowej <= 1.0;
  przedzial tau in (450, 600] raportowany poza kryterium (zapas).
  gamma_s: fit liniowy ln(max_disp) vs tau na odcinku max_disp in
  [0.05, 1.0] (jesli osiagalny); przy braku wzrostu raport "brak
  widocznej ucieczki w oknie". Bramki nie wolno oslabic, skrocic
  ani zreinterpretowac po pierwszym runie; G2 FAIL -> STOP, zero
  runow ugiecia.
- **D16 (obserwacja superpozycji, poza kryteriami):** DEDYKOWANY
  run (omega = 1.5, b = 8, limit 10000 krokow) po zebraniu
  5 probek kontynuowany do x_c > x2_frozen + 10 albo limitu;
  raport jakosciowy: centroid y, szerokosc poprzeczna pakietu
  i kat wektorowy P w oknie x > x2_frozen + 10 (za druga soczewka
  rzedu). Bez statusu kryterialnego; wynik kryterialny punktu
  (1.5, 8) czytany z osobnego runu standardowego D9 (6000 krokow).
