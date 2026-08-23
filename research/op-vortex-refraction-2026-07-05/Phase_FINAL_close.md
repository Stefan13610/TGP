---
title: "Phase_FINAL_close — zamkniecie: B2-prime — refrakcja na defekcie (wirze)"
date: 2026-07-05
type: phase-final-close
tgp_owner: research/op-vortex-refraction-2026-07-05
status: CLOSED
verdict: "G1 FAIL (czesc Phase2, strukturalne); G2 PASS (1.21%); G3 FAIL strukturalne — WSZYSTKIE runy ugiecia niewazne (D9): para wirow ANIHILUJE pod dynamika II rzedu przy tau=107.5, w oknie przelotu pakietu -> STOP wg reguly; G4/G5/G6 NIEROZSTRZYGNIETE (pomiar niewykonalny na tym tle); Phase3 nieuruchomiona"
anti_lakatos_lock: PRESERVED
tags: [gravity-bridge, branch-B2prime, vortex-lensing, refraction, background-collapse, closed]
related:
  - "[[Phase0_balance.md]]"
  - "[[README.md]]"
  - "[[../op-stationary-background-2026-07-05/Phase_FINAL_close.md]]"
  - "[[../op-phi-metric-refraction-2026-07-04/Phase_FINAL_close.md]]"
---

# Phase FINAL — zamkniecie cyklu

## 0. Werdykt

**G2 PASS. G3 FAIL — strukturalne: zaden z dziewieciu runow ugiecia
nie wytworzyl waznego pomiaru (D9: brak 5 probek), bo tlo — para
wirow (+1,-1) — pod PELNA dynamika II rzedu KOLABUJE DO ANIHILACJI
w tau = 107.5, wewnatrz okna przelotu pakietu. G1 (czesc Phase2)
FAIL z tej samej przyczyny (flagi core_lost/no5). STOP wg reguly
decyzyjnej ("G2/G3 FAIL -> stop przed interpretacja"): G4, G5, G6
NIEROZSTRZYGNIETE; Phase 3 nieuruchomiona.**

Glowne USTALENIE cyklu (przetrwa werdykt formalny):

1. **Kwazi-stacjonarnosc pary jest wlasnoscia PRZEPLYWU
   GRADIENTOWEGO, nie substratu falowego.** Cytowany dryf 0.015/tau
   (op-stationary, Phase2b) byl mierzony pod dynamika dyssypatywna
   (predkosc ~ sila). Pod dynamika II rzedu (bezdyssypacyjna) sila
   przyciagania AKUMULUJE predkosc: rdzen przemieszcza sie o 14 do
   tau=50, o 40 do tau=100, anihilacja tau=107.5 (Phase2b, sekcja 1).
   Kontrola dt=0.025: roznica pozycji 0.0074 — kolaps FIZYCZNY,
   nie artefakt integratora. Pre-rejestrowany budzet dryfu #1.4
   ("~3 na przelot tau~200") jest SFALSYFIKOWANY 15-krotnie.
   Argument "separacji skal 25x" porownywal predkosci — dla dynamiki
   inercyjnej to porownanie nie przenosi sie na okno pomiaru.
2. **Strona teoretyczna P1 jest gotowa i spojna:** eikonal na
   RZECZYWISTYM tle pary (Phase0) odtwarza tabele cytowana
   w punktach kryterialnych do ~11% (b=8/12, omega=1.1: -5.8..+4.1%;
   omega=1.3: +10.9/-2.6%; b=16 dodane); kalibracja dyspersji G2:
   blad <= 1.21% (prog 5%). Brakuje wylacznie STABILNEGO tla
   do pomiaru.

## 1. Co testowano

Po op-stationary (sektor falowy wokol wiru czysty, defekt =
soczewka skupiajaca): pierwszy ILOSCIOWY pomiar kanalu
geometrycznego — (P1) zgodnosc kata ugiecia pakietu z eikonalem
substratu (pasmo [0.5, 2.0]); (P2) uniwersalnosc w znaku kretu.
Model bez strojenia (kappa=0.5, a=0.5, b=1.6, c=1.0; N=256, L=128);
tlo: para (+1,-1) szachownicowa; pakiet amplitudowy w lokalnej
ramce fazy tla; korekta dryfu pre-rejestrowana przez b_eff
i referencje blizniacza (D3/D5).

## 2. Przebieg

- **LOCK**: autoryzacja autora (polecenie realizacji cyklu,
  2026-07-05); doprecyzowania pre-code D1-D13 dopisane do sekcji 10
  PRZED pierwsza linijka kodu (m.in. D3: referencja blizniacza
  psi_ref w lockstep; D5: b_eff + alpha_pred ray tracingiem przy
  b_eff; D11: pre-rejestracja dokladnej symetrii sprzezenia dla
  tla lustrzanego).
- **Phase 0** (`Phase0_output.txt`): k0/v_g/lambda; n^2(r) na tle
  pary (zgodne z ogonem A_t/r^2 dla r >= 3); ray tracing na
  rzeczywistym tle vs tabela cytowana (roznice: kryterialne <= 11%,
  b=6/omega=1.0 poza pasmem slabego pola — zgodnie z zastrzezeniem);
  budzet dryfu wg zalozenia #1.4 (~1.2-2.3 na przelot).
- **Phase 1** (`Phase1_output.txt`): **G2 PASS** (max blad 1.21%),
  G1(czesc) PASS: determinizm bitowy (fala + relaksacja), dryft
  sekularny H <= 1.25e-10, szew ansatzu czysty (1.58 = znany
  artefakt C), n_tot = 0. Dryf koncowki relaksacji 0.021/tau
  (cytowane 0.015/tau — spojne, zjazd przyspiesza).
- **Phase 2** (`Phase2_output.txt`): 9 runow (6 kryterialnych + G6
  + 2 obserwacje). WSZYSTKIE niewazne wg D9: pakiet nigdy nie
  przekroczyl okna transmisji (prog x_rdzen+10 "ucieka" wraz
  z kolabujacym rdzeniem), flagi core_lost — detektor windingu
  traci pare po anihilacji. E_dryft = 1.5e-07 (dynamika zdrowa,
  zero flag energetycznych) — problem jest TOPOLOGICZNO-DYNAMICZNY,
  nie numeryczny.
- **Phase 2b** (`Phase2b_output.txt`, diagnostyka w dniu runow,
  przed interpretacja; precedens Phase2b op-stationary):
  separacja(tau): 87 -> 8 na tau in [10, 100], anihilacja 107.5;
  kontrola dt/2: 0.0074; obserwacja poza kryteriami (okno stale):
  klasyfikator kierunku NIEROZSTRZYGNIETY (+0.43 deg — probki
  zbierane w chwili anihilacji soczewki).
- **Phase 3**: NIEURUCHOMIONA — regula decyzyjna nakazuje stop
  po G3 FAIL. D11 pozostaje pre-rejestrowanym przewidywaniem
  (dokladna symetria psi -> psi* wiaze warianty lustrzane).

## 3. Wyniki liczbowe (kluczowe)

```text
G2 (dyspersja, tlo s*): max |blad| = 1.21% (prog 5%)  PASS
eikonal, tlo rzeczywiste vs cytowany (omega=1.1): b=6: +39.1/+43.0;
  b=8: +14.6/+14.0; b=12: +5.0/+5.3; b=16: +2.7 (nowy punkt)
kolaps tla pod dynamika II rzedu (Phase2b):
  tau:        10    30    50    70    90    100   107.5
  separacja:  86.8  76.0  60.2  41.1  19.8  8.2   ANIHILACJA
  przemieszczenie rdzenia: 14.1 (tau=50), 40.2 (tau=100)
  kontrola integratora: |pos(dt=0.05) - pos(dt=0.025)| = 0.0074
budzet #1.4 (zalozenie LOCKa): ~3 na przelot  ->  pomiar: ~40  (x13)
runy Phase 2: 9/9 niewaznych (brak 5 probek; core_lost)
klasyfikator kierunku (poza kryteriami, tlo kolabujace): +0.43 deg
  = nierozstrzygniety (soczewka w ruchu 0.2-0.3 vs v_g = 0.37
  i anihilujaca w oknie probek)
```

## 4. Interpretacja w granicach modelu

CO WYNIK ZNACZY:
- **para wirow na torusie NIE jest tlem pomiarowym dla sektora
  falowego**: bez dyssypacji przyciaganie topologiczne prowadzi do
  inercyjnego kolapsu w czasie krotszym niz przelot pakietu;
  zjazd, ktory pod przeplywem gradientowym byl "wolny" (0.015/tau),
  pod dynamika falowa jest przyspieszajacym spadkiem swobodnym,
- **falsyfikacja dotyczy TLA, nie kanalu geometrycznego**: eikonal
  z substratu jest policzalny i wewnetrznie spojny (Phase0),
  sektor falowy skalibrowany (G2); pytania P1/P2 pozostaja OTWARTE
  — pomiar wymaga tla o czasie zycia >> przelotu pakietu,
- kolaps jest trzecim wcieleniem tej samej przeszkody z op-stationary
  (G3 FAIL): scisla stacjonarnosc pary na skonczonym torusie jest
  niedostepna; nowe ustalenie: rowniez KWAZI-stacjonarnosc jest
  nieprzenoszalna na dynamike II rzedu przy tej separacji (L=128).

CZEGO WYNIK NIE ZNACZY:
- NIE falsyfikuje kanalu geometrycznego wokol defektow (brak
  waznego pomiaru alpha; regula: stop przed interpretacja G4),
- NIE potwierdza ani nie obala uniwersalnosci (G5 niewykonane),
- zadnych deklaracji GR/soczewkowania obserwacyjnego; tau != czas.

## 5. Errata / transparencja

1. **Budzet dryfu #1.4 sfalsyfikowany pomiarem** (zalozenie ~3,
   pomiar ~40 na przelot): zalozenie pochodzilo z ekstrapolacji
   dryfu gradientowego na dynamike falowa — blad metodologiczny
   LOCKa wykryty pierwszym runem i udokumentowany (Phase2b).
   Zgodnie z protokolem NIE zmieniano progow ani procedur; runy
   pozostaja niewazne, wynik raportowany jako rozstrzygniecie
   strukturalne.
2. **Okno "x > x_rdzenia+10" z progiem podazajacym** (D4) okazalo
   sie nieosiagalne na kolabujacym tle (prog ucieka przed pakietem).
   Diagnostyczna wersja z oknem STALYM (Phase2b, poza kryteriami)
   tez nie daje interpretowalnego kata — problem lezy w tle,
   nie w definicji okna.
3. Dryf koncowki relaksacji zmierzony w Phase1: 0.021/tau vs
   cytowane 0.015/tau — spojne (zjazd przyspieszajacy; op-stationary
   mierzyl srednia wczesniejszego odcinka).
4. `Phase2_output.txt` (runy niewazne) zachowany w calosci —
   zadnych powtorek po diagnozie; Phase2b to wylacznie diagnostyka
   + obserwacja poza kryteriami.

## 6. Decyzja i nastepny krok

Wg reguly LOCKa (G3 FAIL -> stop przed interpretacja): cykl
zamkniety bez rozstrzygniecia P1/P2. Zmierzone przeslanki definiuja
nastepny op (za osobna autoryzacja) — **B2-prime na tle o dlugim
czasie zycia**; draft gotowy:
[[../op-lattice-background-2026-07-05/Phase0_balance.md]]
(szachownica 4 wirow, sila netto = 0 z symetrii, bramka G2 czasu
zycia tla). Kandydaci rozwazani do nowego locka:

1. **L=256 przy tej samej fizyce**: sila przyciagania ~1/d, czas
   kolapsu rosnie ~liniowo z separacja; do przeliczenia w Phase0
   nowego opa, czy okno przelotu miesci sie z zapasem >= 3x.
2. **Pomiar w ukladzie wspolporuszajacym sie rdzenia** (alpha
   wzgledem chwilowej pozycji soczewki, okno skrocone tau < 50,
   start pakietu blizej): wymaga pre-rejestracji nowej obserwabli.
3. **Konfiguracje o zerowej sile netto z symetrii** (np. siatka 4
   wirow szachownicowa): do rozstrzygniecia w locku, czy wiez
   symetrii nie narusza zakazu pinningu.
4. Wariant hybrydowy: slaba dyssypacja tla (tlumienie tylko modow
   dlugofalowych tla, nie pakietu) — wymaga starannej analizy
   forbidden moves (ingerencja w dynamike!).

Bilans mostu do grawitacji po siedmiu opach: kanal geometryczny
wokol defektow ma policzalna strone teoretyczna (eikonal z substratu,
n>1, ogon 1/r^2, slepy na znak kretu) i czysty sektor falowy;
pomiar ilosciowy czeka na tlo, ktore zyje dluzej niz pakiet leci.

## 7. Pliki cyklu

- [[Phase0_balance.md]] — LOCK (sekcja 10: D1-D13)
- `Phase0_check.py` / `Phase0_output.txt` — k0/v_g, n^2(r) na tle
  pary, eikonal rzeczywisty vs cytowany, budzet dryfu
- `Phase1_calibration.py` / `Phase1_output.txt` — G2 PASS,
  background_check, determinizm
- `Phase2_deflection.py` / `Phase2_output.txt` — 9 runow niewaznych
  (kolaps tla; G3 FAIL strukturalne)
- `Phase2b_background_lifetime.py` / `Phase2b_output.txt` —
  diagnostyka kolapsu (anihilacja tau=107.5; kontrola dt) +
  obserwacja poza kryteriami
- [[README.md]] — przeglad opa
- ten plik — zamkniecie
