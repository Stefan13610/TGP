---
title: "op-lattice-background (B2-prime, proba #2) — refrakcja na defekcie na tle szachownicowym o zerowej sile netto"
date: 2026-07-05
type: research-op
status: CLOSED
verdict: "G2 BRAMKA PASS (tlo zyje: 0.0000 przemieszczenia przez tau=400); G3 PASS (1.21%); G4 PASS (9/9 runow waznych, znak KU wirowi wszedzie); G5 FAIL formalnie (5/6 punktow w pasmie [0.5,2.0]; punkt 1.1/8 ratio 0.445 w rezimie b_eff ~ lambda); G6a PASS (do 4e-04 deg / bitowo); G6b PASS"
tags: [gravity-bridge, branch-B2prime, lattice-background, vortex-lensing, refraction, eikonal, lens-superposition, closed]
related:
  - "[[Phase0_balance.md]]"
  - "[[Phase_FINAL_close.md]]"
  - "[[../op-vortex-refraction-2026-07-05/Phase_FINAL_close.md]]"
  - "[[../op-stationary-background-2026-07-05/Phase_FINAL_close.md]]"
---

# op-lattice-background-2026-07-05 (B2-prime, proba #2)

**Status: CLOSED (2026-07-05). LOCK: autoryzacja autora (polecenie
realizacji cyklu); D1-D16 spisane przed kodem; anty-Lakatos
zachowany. Werdykt i pelny raport: [[Phase_FINAL_close.md]].**

Pytanie: **czy ugiecie pakietu amplitudowego na wirze zgadza sie
ilosciowo z eikonalem substratu (P1; pasmo [0.5, 2.0]) — na tle,
ktore zyje dluzej niz pakiet leci?** Proba #1
([[../op-vortex-refraction-2026-07-05/Phase_FINAL_close.md]]) padla
wylacznie na tle: para (+1,-1) anihiluje pod dynamika II rzedu przy
tau = 107.5. To tlo (szachownica 4 wirow, sila netto = 0 DOKLADNIE
z symetrii D4 otoczenia kazdego wiru na torusie) usuwa przyczyne
konstrukcja, nie sila.

**WYNIK:**
- **BRAMKA G2 PASS — problem tla ROZWIAZANY**: szachownica pod
  pelna dynamika II rzedu nie rusza rdzeni W OGOLE (max 0.0000
  przez tau <= 400; dryft energii 1.7e-16). Program ma pierwsze
  stabilne tlo pomiarowe dla sektora falowego.
- **G4 PASS — pierwsze WAZNE pomiary refrakcji w programie**:
  9/9 runow waznych (proba #1: 0/9); ugiecie KU wirowi we
  wszystkich konfiguracjach; monotonia w b i omega; liniowosc
  dokladna (G6b: 0.000 deg); symetrie dokladne (G6a: 3.8e-04 deg,
  lustro bitowo 0 — pre-rejestracja D11 potwierdzona co do bitu).
- **G5 FAIL formalnie, 5/6 w pasmie**: ratio alpha/alpha_pred(b_eff)
  = 1.19-1.46 (b_eff >= 11) i 0.741 (1.3/8); punkt (1.1, b=8) —
  najglebiej dyfrakcyjny, b_eff = 7.17 ~ lambda = 7.72 — ratio
  0.445, dobry znak. Granica stosowalnosci eikonalu zmierzona:
  b_eff ~ lambda (zastrzezenie D12 zapisane z gory).
- Wg reguly decyzyjnej: raport rozjazdu w zamknieciu + **powtorka
  krotkofalowa** (wyzsze omega / L=256) jako osobny op (osobna
  autoryzacja) — draft gotowy:
  [[../op-shortwave-lattice-2026-07-05/Phase0_balance.md]];
  pelny P2 (tlo asymetryczne) niezaleznie.

Szczegoly, kryteria G1-G6 (w tym NOWA bramka G2: czas zycia tla),
regula decyzyjna: [[Phase0_balance.md]]. Pliki: `Phase0_check.py`/
`Phase0_output.txt` (bilans sil 0.0 per powloka D4; eikonal
wiazacy), `Phase1_gate_and_calibration.py`/`Phase1_output.txt`
(BRAMKA G2 + G3), `Phase2_deflection.py`/`Phase2_output.txt`
(skan glowny), `Phase3_symmetry.py`/`Phase3_output.txt` (G6a),
[[Phase_FINAL_close.md]] (zamkniecie).

## Notatki wykonawcze dla nowego agenta (doswiadczenie sesji #B1-#B2prime-1)

**Protokol (bezwzgledny):**
1. NAJPIERW autoryzacja autora -> status LOCKED -> doprecyzowania
   pre-code do sekcji 10 -> DOPIERO kod. Wzorzec doprecyzowan:
   sekcja 10 LOCKa `../op-vortex-refraction-2026-07-05/`
   (D1-D13) — wiekszosc przenosi sie wprost.
2. Zadnych zmian progow/procedur po pierwszym runie; czulosci
   wylacznie "poza kryterium"; erraty w dniu runow, przed
   interpretacja. BRAMKI G2 nie wolno oslabic po pierwszym runie;
   przy G2 FAIL zadnych runow ugiecia.
3. Kazda faza: skrypt + output.txt w katalogu opa; na koncu
   Phase_FINAL_close.md + aktualizacja tego README.
4. KOLEJNOSC FAZ JEST CZESCIA PROTOKOLU: bramka G2 (czas zycia)
   PRZED jakimkolwiek runem z pakietem.

**Znane pulapki techniczne (wszystkie przecwiczone w probie #1 —
nie odkrywaj ich ponownie):**
- **Dryf gradientowy NIE przewiduje dryfu falowego** (glowna lekcja
  proby #1): kwazi-stacjonarnosc pod przeplywem gradientowym
  (predkosc ~ sila) nie przenosi sie na dynamike II rzedu
  (przyspieszenie ~ sila; bez dyssypacji predkosc sie AKUMULUJE).
  Czas zycia tla mierzy sie POD DYNAMIKA POMIAROWA (bramka G2),
  nigdy nie ekstrapoluje z relaksacji.
- **Prog okna pomiaru NIE moze podazac za rdzeniem** (proba #1:
  prog "x_rdzen+10" z pozycja biezaca uciekal przed pakietem —
  wyzwolenie nigdy nie nastapilo). Okno zamrazac w chwili
  wyzwolenia; wyzwolenie centroidem.
- **Referencja blizniacza (D3 proby #1) dziala bardzo dobrze**:
  psi_ref bez pakietu w lockstep tym samym integratorem;
  delta_psi = psi - psi_ref; E_dryft ~ 1.5e-07. Uzywac.
- **Ansatz siatki**: suma DWOCH par poziomych theta_h_pair (rzad
  y=-32: +1/-1; rzad y=+32: -1/+1) — konstrukcja dokladna D1 z B1,
  BEZ transpozycji, wiec BEZ artefaktu punktu C (nowosc vs para
  diagonalna). Kod wzorcowy: `theta_h_pair` w
  `../op-vortex-refraction-2026-07-05/Phase2_deflection.py`
  (tam tez: caly lockstep, b_eff, pomiar P, detektor rdzeni).
  Po zbudowaniu ansatzu ZAWSZE diagnostyka szwu: max skok fazy na
  linkach poza rdzeniami << pi (szew y: ogony dipolowe ~0.04 rad
  oczekiwane — raport).
- **Leapfrog**: psi_tau roznica CENTRALNA; dryft energii = trend
  sekularny (srednie cwiartek), oscylacja O((dt*omega)^2) osobno.
- **Pomiar kata**: usredniac WEKTOROWO (Px, Py), potem kat;
  konwencja alpha := -sign(b)*atan2(Py,Px), alpha > 0 = KU wirowi.
- **Detektor wirow**: winding plakietowy + centroid deficytu Phi;
  UWAGA: na siatce sa 4 rdzenie — detektor z proby #1 bierze
  cells[0] per znak; rozszerzyc na klastrowanie 2+2 (pozycje
  wszystkich czterech, przypisanie do nominalow po najblizszej
  odleglosci periodycznej). Pozycje zawsze z pola.
- **Eikonal**: ray tracing na RZECZYWISTYM zrelaksowanym tle
  (bilinear n^2 + gradient), alpha_pred przy b_eff; NOWE: kat
  promienia ewaluowac w plaszczyznie pomiaru (3. probka), czulosc
  na krance okna raportowac — za pierwsza soczewka pole nie jest
  jeszcze plaskie (druga soczewka rzedu w x=+32).
- **Symetrie siatki**: alpha(+b) = alpha(-b) i rownowaznosc lustra
  sa DOKLADNE z konstrukcji (odbicie osi wiazki + sprzezenie;
  translacja o L/2) — to testy IMPLEMENTACJI (G6a), nie fizyki;
  alpha_odd = 0 z konstrukcji, NIE raportowac jako pomiar P2.
- **Narzedzia**: nie uzywac `cd` w bashu; sciezki pelne od korzenia
  vaulta; outputy przez `> plik.txt`; python + numpy (sprawdzone:
  3.14 / 2.4).

**Kluczowe liczby (cytowane):**
- substrat (op-stationary Phase0): s* = 1.174166, Phi* = 1.378665,
  V''(s*) = 0.878665, m_TV = 0.9374, xi^2 = 0.5690, A_t = 2.5690
- kanal amplitudowy: k0/v_g/lambda (omega=1.1): 0.814/0.370/7.72;
  (omega=1.3): 1.274/0.490/4.93
- eikonal orientacyjny (tlo PARY, proba #1 Phase0): omega=1.1:
  +14.6/+5.0/+2.7 deg (b=8/12/16); omega=1.3: +5.8/+2.1/+1.1
  (zaklad wiazacy: Phase0_check tego opa, na tle siatki)
- kalibracja dyspersji (proba #1, identyczny sektor): 1.21% <= 5%
- KONTRAST dla bramki G2: para (+1,-1) pod dynamika II rzedu —
  przemieszczenie 14 przy tau=50, anihilacja tau=107.5
  (`../op-vortex-refraction-2026-07-05/Phase2b_output.txt`)
