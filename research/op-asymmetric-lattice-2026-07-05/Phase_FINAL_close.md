---
title: "Phase_FINAL_close — B2-prime proba #4 (op-asymmetric-lattice): STOP na bramce G2 — siatka skosna dynamicznie niezdatna (niestabilny mod scinajacy, czas zycia tau=227)"
date: 2026-07-12
type: phase-final-close
tgp_owner: research/op-asymmetric-lattice-2026-07-05
status: CLOSED (STOP na bramce G2, wg reguly decyzyjnej LOCKa)
anti_lakatos_lock: RELEASED (cykl zamkniety; lock dotrzymany w calosci)
tags: [gravity-bridge, branch-B2prime, oblique-lattice, shear-instability, gate-stop, P2, L256]
related:
  - "[[Phase0_balance.md]]"
  - "[[README.md]]"
  - "[[../op-shortwave-lattice-2026-07-05/Phase_FINAL_close.md]]"
---

# FINAL close — B2-prime proba #4: STOP na bramce G2

## 1. Werdykt (sciezka reguly decyzyjnej LOCKa)

**G2 FAIL -> STOP CALEGO CYKLU przed pomiarem** (sekcja 5 LOCKa,
pierwsza galaz reguly decyzyjnej; forbidden move #2 dotrzymany:
bramka nieoslabiona, zero runow ugiecia).

```text
G1 (czesc Phase 1, do momentu bramki): PASS czesciowy
    determinizm bitowy relaksacji: True; H_mono: True; n_tot = 0;
    detektor 4+4 widzi 8 rdzeni w ansatzu i po relaksacji;
    dryft sekularny energii w runie bramki: 2.12e-08 <= 1e-4
G2 (BRAMKA czasu zycia siatki skosnej): FAIL
    kryterium: max przemieszczenie KAZDEGO z 8 rdzeni <= 1.0
    dla tau <= 450; zmierzone: 63.31 (wszystkie 8 rdzeni),
    flaga core_lost przy tau = 227.0
G3: POMINIETE (wg skryptu zalockowanego: G2 FAIL -> stop przed G3)
G4, G5a, G5b, G6a, G6b: NIE MIERZONE (zero runow z pakietem)
```

To jest scenariusz pre-rejestrowany jako RYZYKO #1 cyklu (sekcja
1.2 LOCKa: "tempo NIEZNANE dla siatki skosnej (mody scinajace moga
byc miekksze niz w szachownicy!) — rozstrzyga bramka G2, nie
zalozenie"). Bramka zadzialala dokladnie tak, jak zaprojektowano.

## 2. Deliverable: czasy zycia (pomiar bramki)

```text
tlo                          przemieszczenie / czas zycia
---------------------------------------------------------------
siatka SKOSNA L=256 (8 wirow, scinanie 64/rzad, ten cykl):
    tau= 60: max disp = 14.36;  tau=120: 41.29;  tau=180: 59.98
    UTRATA RDZENI (anihilacja par w rzedach) przy tau = 227.0
    gamma_s (fit ln dmax vs tau, odcinek [0.05, 1.0]) = +0.284/tau
    (prog disp=1.0 przekroczony juz przy tau ~ 25-30;
     przelot pomiarowy wymagal tau ~ 140 przy omega=1.1)
szachownica L=256 (proba #3, cytowane): 0.0000 przez tau = 600
para +1/-1 (proba #1, cytowane):        anihilacja tau = 107.5
```

Charakter modu (z pozycji rdzeni bramki): WSZYSTKIE 8 rdzeni
przemieszcza sie o te sama odleglosc — podsiec +1 sunie w +x,
podsiec -1 w -x (translacyjny mod scinajacy podsieci wzdluz
rzedow). Przy przemieszczeniu ~63 kazdy wir +1 dochodzi do
pozycji wira -1 swojego rzedu (rozstaw 128) -> anihilacje
w rzedach -> core_lost (tau = 227.0).

Zapowiedz w Phase 0 (spisana przed bramka): juz relaksacja
gradientowa (2000 krokow) przesunela podsieci o +/-0.69 w x
(residuum 5.68e-03 vs 9.6e-07 na szachownicy proby #3) — sila
netto w punkcie symetrycznym jest DOKLADNIE zero (kasowanie
parami inwersyjnymi v/-v, zweryfikowane w Phase0_check do zera
maszynowego), ale punkt symetryczny jest SIODLEM: antysymetryczne
przesuniecie podsieci obniza energie. Zero sily nie implikuje
stabilnosci — dokladnie dlatego bramka G2 byla warunkiem wstepnym.

## 3. Pozostale wyniki liczbowe (Phase 0, wiazace — wykonane przed bramka)

- Tabela projektowa: wszystkie 10 punktow kryterialnych w domenie
  eikonalu (b_nom >= 1.5 lambda + 0.9; odleglosci od soczewek
  posrednich 36-56 >> 1.5 lambda). Projekt punktow byl poprawny.
- Bilans sil parami inwersyjnymi: 0.000e+00 per para, per powloka
  indeksowa K=1..8 i w sumie; kontrola negatywna (polplaszczyzna
  bez partnerow -v): |F| = 7.68e-03 — kasowanie jest wlasnoscia
  par v/-v, nie orbit D4; kontrola pudelkowa K=8: artefakt
  powierzchniowy 1.23e-02 (identyczny na wszystkich 8 rdzeniach —
  otoczenia rownowazne z inwersji).
- Lustro: pozycje 8/8 lustrzane wzgledem y=-64; ladunki NIEZGODNE
  w 4/8 (oba rzedy posrednie) -> konstrukcyjny cel tla osiagniety
  (alpha_odd niewiazana zadna symetria) — ale tlo nie zyje.
- S-inv (inwersja wokol A): 8/8 pozycji+ladunkow; S-conj
  (sprzezenie = translacja o (128,0)): 8/8.
- Ray tracing wiazacy (na tle po relaksacji 2000): Delta_eik
  w x=0: -1.08 deg (1.1, |b|=14), -1.08 (|b|=20), -1.11 (|b|=24),
  -0.23 (1.5, |b|=8), -0.25 (1.5, |b|=12). UWAGA: niezerowe
  glownie przez dryft relaksacyjny podsieci (pozycje ZRELAKSOWANE
  utracily scisla lustrzanosc) — kolejny, niezalezny symptom
  niestabilnosci scinania, spisany przed bramka.

## 4. Interpretacja w granicach modelu

**CO ZNACZY:**
- Siatka skosna (szachownica scinana o L/4 = 64 na rzad) przy
  L=256 jest dynamicznie NIEZDATNA jako tlo pomiarowe: ma
  niestabilny miekki mod scinajacy podsieci (gamma_s ~ 0.28/tau,
  zycie tau ~ 227 << wymagane 450), mimo dokladnie zerowej sily
  netto na kazdym wirze. Konfiguracja jest siodlem, nie minimum.
- Program bramkowy dziala: niestabilnosc zostala wykryta PRZED
  jakimkolwiek runem pomiarowym; zadne dane ugiecia nie powstaly,
  wiec zadne nie wymagaja uniewaznienia.
- Pomiar pelnego P2 (kanal cyrkulacyjny) WYMAGA innego tla:
  konstrukcja "pozycje lustrzane + ladunki lamiace lustro"
  pozostaje poprawnym pomyslem eksperymentalnym (Phase0_check
  potwierdzil wlasnosci symetryczne konfiguracji), ale ta
  konkretna realizacja geometryczna nie utrzymuje sie sama
  (forbidden move #1: zero pinningu — respektowane).

**CZEGO NIE ZNACZY:**
- NIE falsyfikuje P2 ani mostu geometrycznego: kanal cyrkulacyjny
  nie zostal zmierzony (ani wykluczony) — pytanie P2-b pozostaje
  OTWARTE, dokladnie tak jak przed cyklem.
- NIE podwaza P1 (rozstrzygnietego w probach #2/#3 na stabilnej
  szachownicy) ani zadnego wczesniejszego wyniku.
- NIE jest wynikiem "negatywnym dla substratu": niestabilnosc
  siodlowa konkretnej konfiguracji wirow to wlasnosc dynamiki
  tla, nie kanalu propagacji.
- Zero deklaracji GR/obserwacyjnych; tau != czas (forbidden
  move #6).

## 5. Errata / transparencja (dzien runow, przed interpretacja)

1. `Phase2_deflection.py` i `Phase3_symmetry.py` zostaly NAPISANE
   rownolegle z trwajaca bramka (przygotowanie), ale NIGDY NIE
   URUCHOMIONE — brak `Phase2_output.txt`/`Phase3_output.txt`
   jest zgodny z protokolem (G2 FAIL -> zero runow ugiecia).
   Pliki pozostawione w katalogu jako zapis pre-rejestrowanej
   implementacji (w tym pakiet lewobiezny D17 i definicje D18).
2. Zadnych zmian progow, parametrow ani procedur po pierwszym
   runie; bramka wykonana w pelnym oknie tau=600 zgodnie z D15
   (run przerwal sie sam na core_lost przy tau=227 — flaga
   detektora, nie skrocenie bramki).
3. Diagnostyka D11 w Phase0_check (max|psi - I[psi]| = 2.35)
   jest zdominowana przez stala/gladka roznice fazy (gauge)
   i dryft relaksacyjny; nie byla kryterium i nie zostala uzyta
   do zadnej decyzji (testy G6a sie nie odbyly).
4. Oszacowanie orientacyjne draftu "|Delta_eik| < 0.1 deg"
   okazalo sie bledne na tle ZRELAKSOWANYM (pomiar wiazacy:
   do -1.11 deg) — rozjazd spisany w Phase0_output przed bramka;
   bez wplywu na werdykt (G5b nie bylo mierzone; procedura D18
   i tak porownywala Delta_meas z wiazacym Delta_eik, nie z 0).

## 6. Decyzja i nastepny krok

Cykl ZAMKNIETY ze STOP. Zgodnie z regula decyzyjna LOCKa,
kandydaci na nastepny op (kazdy za OSOBNA autoryzacja, osobny
lock; zadnej kontynuacji w tym cyklu):

1. **siatka skosna o mniejszym scinaniu (s = 32)** — blizej
   szachownicy; mod scinajacy moze odzyskac stabilnosc, a ladunki
   rzedow posrednich nadal lamia lustro (do sprawdzenia w drafcie:
   czy s=32 zachowuje lustrzanosc pozycji wzgledem osi wiazki);
2. **L = 512** — wieksze separacje (slabsze sprzezenie miedzy
   rzedami; tempo modu scinajacego ~1/d^2 — do oszacowania
   w Phase 0 nastepnego draftu);
3. alternatywnie: inna konstrukcja tla "pozycje lustrzane +
   ladunki lamiace lustro" o sztywnym szkielecie (bez pinningu,
   forbidden move #1 pozostaje).

Deliverable czasow zycia (sekcja 2) jest wkladem pomiarowym tego
cyklu do projektowania kolejnych tel.

## 7. Pliki cyklu

- `Phase0_balance.md` — LOCKED 2026-07-12 (sekcja 10: D1-D19
  spisane przed kodem)
- `Phase0_check.py` / `Phase0_output.txt` — rachunki wiazace
  (wykonane)
- `Phase1_gate_and_calibration.py` / `Phase1_output.txt` —
  bramka G2: FAIL (wykonane; G3 pominiete wg reguly)
- `Phase2_deflection.py`, `Phase3_symmetry.py` — przygotowane,
  NIE uruchomione (errata #1)
- `Phase_FINAL_close.md` — ten plik
- `README.md` — zaktualizowany (status STOP)
