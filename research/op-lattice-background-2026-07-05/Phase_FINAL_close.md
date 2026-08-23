---
title: "Phase_FINAL_close — zamkniecie: B2-prime proba #2 — refrakcja na defekcie na tle szachownicowym"
date: 2026-07-05
type: phase-final-close
tgp_owner: research/op-lattice-background-2026-07-05
status: CLOSED
verdict: "G1 PASS; G2 BRAMKA PASS (tlo zyje >= tau=400, przemieszczenie 0.0000 — problem tla ROZWIAZANY konstrukcja); G3 PASS (1.21%); G4 PASS (pierwsze WAZNE pomiary refrakcji w programie; znak KU wirowi we wszystkich runach); G5 FAIL formalnie — 5/6 punktow w pasmie [0.5,2.0], punkt (omega=1.1, b=8) ratio 0.445 przy dobrym znaku (rezim gleboko dyfrakcyjny, lambda ~ b_eff; zastrzezenie D12) -> wg reguly: raport rozjazdu, powtorka krotkofalowa jako osobny op; G6a PASS (symetrie do 4e-04 deg / bitowo); G6b PASS (liniowosc dokladna)"
anti_lakatos_lock: PRESERVED
tags: [gravity-bridge, branch-B2prime, lattice-background, vortex-lensing, refraction, eikonal, lens-superposition, closed]
related:
  - "[[Phase0_balance.md]]"
  - "[[README.md]]"
  - "[[../op-vortex-refraction-2026-07-05/Phase_FINAL_close.md]]"
  - "[[../op-stationary-background-2026-07-05/Phase_FINAL_close.md]]"
---

# Phase FINAL — zamkniecie cyklu

## 0. Werdykt

**G2 BRAMKA PASS — przeszkoda, ktora zabila probe #1, jest USUNIETA
KONSTRUKCJA: szachownica 4 wirow pod pelna dynamika II rzedu nie
przemieszcza rdzeni w ogole (max 0.0000 przez tau <= 400; kontrast:
para — 14 przy tau=50, anihilacja 107.5). G4 PASS: pierwsze WAZNE
ilosciowe pomiary kanalu geometrycznego w programie — wszystkie runy
daja ugiecie KU wirowi (znak grawitacji), malejace z b i omega,
liniowe (G6b: roznica 0.000 deg), respektujace dokladne symetrie tla
(G6a: 3.8e-04 deg / bitowo 0). G5 FAIL FORMALNIE: 5 z 6 punktow
kryterialnych w pasmie [0.5, 2.0] wzgledem eikonalu substratu;
punkt (omega=1.1, b=8) — najglebiej dyfrakcyjny (b_eff = 7.17 ~
lambda = 7.72) — daje ratio 0.445, tuz pod pasmem, przy poprawnym
znaku. Wg reguly decyzyjnej: rozjazd raportowany, powtorka
krotkofalowa jako osobny op.**

Glowne USTALENIA cyklu:

1. **Problem tla pomiarowego jest ROZWIAZANY** (konstrukcja, nie
   sila): kasowanie sil z symetrii D4 na torusie dziala rowniez pod
   bezdyssypacyjna dynamika II rzedu; niestabilnosc siodlowa NIE
   ucieka z szumu zaokraglen w oknie 3x przelot (dryft energii runu
   bramki: 1.7e-16). Program ma pierwsze stabilne tlo defektowe
   do pomiarow falowych.
2. **Kanal geometryczny ma ZMIERZONY znak grawitacji**: 9/9 waznych
   runow ugiecia (6 kryterialnych + 3 dodatkowe) — pakiet uginany
   KU defektowi, niezaleznie od b, omega, amplitudy i wariantu
   symetrycznego. To pierwszy komplet waznych pomiarow refrakcji
   na defekcie w programie (proba #1: 0/9 waznych).
3. **Eikonal substratu jest ilosciowo czytelny poza strefa
   dyfrakcyjna**: b_eff >= 11 -> ratio 1.19-1.46; b_eff ~ 7.5 przy
   omega=1.3 (lambda=4.93) -> 0.741; b_eff ~ 7.2 przy omega=1.1
   (lambda=7.72) -> 0.445 (poza pasmem). Granica stosowalnosci
   lezy tam, gdzie zapisano ja z gory (D12): b_eff ~ lambda.

## 1. Co testowano

Po probie #1 (tlo pary anihilowalo przed pomiarem): te same pytania
P1 (policzalnosc kata ugiecia z eikonalu substratu, pasmo [0.5,
2.0]) — teraz w geometrii WIELOSOCZEWKOWEJ (eikonal na pelnym tle
siatki) — oraz czesciowo P2 (zgodnosc z dokladnymi symetriami tla).
Model bez strojenia (kappa=0.5, a=0.5, b=1.6, c=1.0; N=256, L=128);
tlo: szachownica 4 wirow (+1,-1/-1,+1) w rozstawie L/2 — kazdy wir
w punkcie o symetrii D4 otoczenia -> sila netto = 0 DOKLADNIE;
NOWA bramka G2: czas zycia tla pod dynamika pomiarowa PRZED
jakimkolwiek runem z pakietem.

## 2. Przebieg

- **LOCK**: autoryzacja autora (polecenie realizacji cyklu,
  2026-07-05); doprecyzowania pre-code D1-D16 dopisane do sekcji 10
  PRZED pierwsza linijka kodu (D14: detektor 2+2 z klastrowaniem;
  D15: procedura bramki G2; D16: obserwacja superpozycji).
- **Phase 0** (`Phase0_output.txt`): bilans sil — powloki D4
  centrowane na wirze kasuja sie DOKLADNIE (0.0 kazda powloka;
  kontrast pary: |F| = 1.7e-02); tlo zrelaksowane: residuum
  1.6e-06 (para proby #1: 8.1e-03 — siatka jest STACJONARNA,
  nie quasi-); zaklad eikonalny wiazacy na tle siatki + czulosc
  plaszczyzny ewaluacji (0.4-3.8 deg); budzet czasu.
- **Phase 1** (`Phase1_output.txt`): szew ansatzu 0.124 rad << pi
  (BEZ artefaktu punktu C); **BRAMKA G2 PASS**: max przemieszczenie
  rdzeni 0.0000 (tau <= 300 i tau <= 400), n_tot = 0, 4 rdzenie
  przez cale okno, dryft energii 1.7e-16, gamma_s: brak widocznej
  ucieczki; **G3 PASS**: dyspersja max |blad| 1.21% (prog 5%);
  determinizm bitowy (relaksacja + fala). G1(czesc) PASS.
- **Phase 2** (`Phase2_output.txt`): 9 runow, WSZYSTKIE WAZNE
  (5 probek kazdy; zadnych flag; E_dryft <= 1.3e-13; gamma_core
  <= 0.004; kret zachowany). G4 PASS; G5 5/6 (FAIL formalnie);
  G6b PASS (roznica 0.000 deg przy u0/2).
- **Phase 3** (`Phase3_output.txt`): G6a PASS — b=-8 zgodne do
  3.8e-04 deg; tlo lustrzane: pole po relaksacji ROWNE sprzezeniu
  (max roznica 0.0e+00), run bitowo identyczny (pre-rejestracja
  D11 potwierdzona co do bitu). alpha_odd = -1.9e-04 deg = 0
  z konstrukcji (NIE raportowane jako pomiar P2, forbidden move #5).

## 3. Wyniki liczbowe (kluczowe)

```text
BRAMKA G2 (tlo bez pakietu, dynamika II rzedu, tau=400):
  max przemieszczenie rdzenia: 0.0000 (wszystkie 4 rdzenie)
  dryft sekularny energii: 1.67e-16;  kontrast pary: anihilacja 107.5
G3 (dyspersja, tlo s*): max |blad| = 1.21% (prog 5%)  PASS
skan glowny (alpha zmierzone vs alpha_pred(b_eff) [deg]; ratio):
  omega=1.1 b= 8: +9.005+/-0.300  vs +20.235  ratio 0.445  POZA
  omega=1.1 b=12: +6.657+/-0.170  vs  +5.602  ratio 1.188  OK
  omega=1.1 b=16: +3.565+/-0.132  vs  +2.609  ratio 1.366  OK
  omega=1.3 b= 8: +4.590+/-0.191  vs  +6.197  ratio 0.741  OK
  omega=1.3 b=12: +3.110+/-0.225  vs  +2.125  ratio 1.464  OK
  omega=1.3 b=16: +1.497+/-0.078  vs  +1.045  ratio 1.432  OK
G4: |alpha(8,1.1)| = 9.01 > 3 sigma = 0.90; monotonia w b i omega: TAK
G6a: |alpha(+8)-alpha(-8)| = 3.8e-04 deg; lustro: 0.0 (bitowo)
G6b: |alpha(u0)-alpha(u0/2)| = 0.000 deg (prog 0.90)
obserwacje poza kryteriami:
  b=6, omega=1.1: alpha = +8.67 (pred 12.25, ratio 0.71) — silne pole
  omega=1.0, b=8: alpha = +12.19; eikonal przy b_eff=6.37: PRZECHWYT
  superpozycja (D16): centroid delta_psi nie osiagnal x2+10 w limicie
    D9 (tau=300) — obserwacja niewykonana (raport, bez przedluzania)
```

## 4. Interpretacja w granicach modelu

CO WYNIK ZNACZY:
- **tlo szachownicowe jest pierwszym dzialajacym tlem pomiarowym**
  dla sektora falowego wokol defektow: zero dryfu, zero anihilacji,
  wszystkie pomiary wazne — lekcja proby #1 (dryf gradientowy nie
  przewiduje dryfu falowego) przekuta w dzialajaca konstrukcje,
- **znak grawitacji POTWIERDZONY pomiarem**: defekt-materia ugina
  fale amplitudowe KU sobie we wszystkich zmierzonych konfiguracjach
  (rowniez na tle wielosoczewkowym); ugiecie liniowe w amplitudzie
  i zgodne z dokladnymi symetriami tla,
- **policzalnosc (P1) potwierdzona w rezimie polfalowym**
  (b_eff >~ 1.5 lambda): 5/6 punktow w pasmie; eikonal z samego
  substratu przewiduje wielkosc ugiecia do czynnika ~1.5,
- **granica stosowalnosci eikonalu zmierzona**: przy b_eff ~ lambda
  (punkt 1.1/8: lambda/b_eff = 1.08) fala usrednia soczewke po
  swojej szerokosci i ugiecie jest ~2.2x mniejsze od promienia
  geometrycznego — dokladnie zastrzezony z gory rezim dyfrakcyjny
  (D12); formalnie G5 FAIL wg zamknietego zestawu szesciu punktow.

CZEGO WYNIK NIE ZNACZY:
- NIE potwierdza P1 w brzmieniu LOCKa (pasmo wymagane we WSZYSTKICH
  szesciu punktach) — brak promocji do "mostu potwierdzonego",
- NIE testuje uniwersalnosci w znaku kretu (P2) nietrywialnie:
  alpha_odd = 0 z konstrukcji tla (pre-rejestrowane; pelny P2 =
  osobny op na tle asymetrycznym),
- zadnych deklaracji GR/soczewkowania obserwacyjnego; tau != czas;
  bez promocji do core/ (osobna autoryzacja).

## 5. Errata / transparencja

1. **Phase0_check, bilans sil — poprawka metody PRZED runami
   pomiarowymi**: pierwotne sumowanie obrazow z obcieciem
   pudelkowym (nie centrowanym na wirze) zostawialo czlon
   powierzchniowy ~K^-2 (|F| ~ 1.9e-05 przy K=8), pozornie
   sprzeczny z teza o dokladnym kasowaniu; po zgrupowaniu w powloki
   centrowane na wirze (unie orbit D4) kazda powloka = 0.0
   DOKLADNIE. Obie wersje raportowane w `Phase0_output.txt`;
   zadne progi nie ulegly zmianie.
2. **Punkt (1.1, b=8) — analiza rozjazdu (poza kryterium)**:
   alpha_pred jest bardzo stromy w b przy b ~ 7 (pred(8.0) = +14.4,
   pred(7.17) = +20.2): korekta b_eff (centroid wciagany 0.84 ku
   rdzeniowi) podnosi przewidywanie o ~40%. Wzgledem b_nom = 8
   ratio wynosiloby ~0.62-0.65 (w pasmie) — raportowane WYLACZNIE
   jako analiza; D5 wiaze ratio przy b_eff i tak pozostaje
   (G5 FAIL formalny). Czulosc plaszczyzny ewaluacji NIE tlumaczy
   rozjazdu (19.6-20.2 deg na wszystkich plaszczyznach).
3. **Obserwacja D16 (superpozycja za druga soczewka) niewykonana**:
   centroid |delta_psi|^2 nie przekroczyl x2+10 do limitu D9
   (tau=300) — ogon pakietu w oknie spowalnia centroid ponizej v_g.
   Zgodnie z protokolem limitu NIE przedluzono; obserwacja
   odnotowana jako niewykonana.
4. **b_eff w wariancie S1 (b=-8)**: 7.651 vs 7.165 (baseline) —
   roznica 0.49 = 2 x 0.25 przesuniecia detektora plakietowego
   wzgledem osi odbicia (pozycja rdzenia -31.75, os -32.00);
   alpha zgodne do 3.8e-04 deg — wlasnosc detektora, nie dynamiki.
5. Zadnych zmian progow, parametrow ani procedur po pierwszym
   runie; zadnych powtorek runow kryterialnych.

## 6. Decyzja i nastepny krok

Wg reguly LOCKa — sciezka "G5 FAIL (wielkosc poza pasmem przy
dobrym znaku)": **raport rozjazdu (sekcje 3-5) + powtorka
KROTKOFALOWA jako osobny op** (za osobna autoryzacja). Draft
gotowy: [[../op-shortwave-lattice-2026-07-05/Phase0_balance.md]]
(seria omega=1.5 + seria 1.1 na b {14,16,20}, L=256, wszystkie
punkty b_eff >= 1.5 lambda). Zmierzone przeslanki do nowego locka:

1. **Rezim krotkofalowy na TYM SAMYM tle** (tlo juz zwalidowane
   bramka G2): omega wyzsze (k0 rosnie; uwaga na dyspersje sieciowa
   — przy k=1.47 blad juz 1.21%, monitorowac prog G3) i/lub
   L = 256 z b in {12, ..., 24}, tak by kazdy punkt mial
   b_eff >= 1.5 lambda. Sam punkt 1.1/8 przechodzi w rezim
   polfalowy przy L=256 (wieksze b przy tej samej geometrii).
2. **Pelny P2 na tle asymetrycznym** (czesc cyrkulacyjna
   nieskasowana symetria) — niezaleznie od punktu 1.
3. (dalsze, wg reguly dla przyszlego PASS: skalowanie z |n|,
   granica newtonowska 2D i G_eff.)

Bilans mostu do grawitacji po osmiu opach: kanal geometryczny wokol
materii-defektow ma (a) stabilne tlo pomiarowe (NOWE), (b) zmierzony
znak uniwersalnie przyciagajacy (NOWE — pierwsze wazne pomiary),
(c) policzalnosc z substratu potwierdzona w 5/6 punktow, z granica
stosowalnosci eikonalu zmierzona na b_eff ~ lambda. Formalnie P1
pozostaje OTWARTE do rozstrzygniecia krotkofalowego.

## 7. Pliki cyklu

- [[Phase0_balance.md]] — LOCK (sekcja 10: D1-D16)
- `Phase0_check.py` / `Phase0_output.txt` — bilans sil (0.0 per
  powloka D4), tlo zrelaksowane (residuum 1.6e-06), n^2, eikonal
  wiazacy + czulosc plaszczyzn, budzet czasu
- `Phase1_gate_and_calibration.py` / `Phase1_output.txt` —
  **BRAMKA G2 PASS** (0.0000 / tau=400), G3 PASS (1.21%),
  background_check, determinizm
- `Phase2_deflection.py` / `Phase2_output.txt` — 9/9 runow waznych;
  G4 PASS, G5 5/6 (FAIL formalny), G6b PASS + obserwacje
- `Phase3_symmetry.py` / `Phase3_output.txt` — G6a PASS
  (S1: 3.8e-04 deg; S2: bitowo 0)
- [[README.md]] — przeglad opa
- ten plik — zamkniecie
