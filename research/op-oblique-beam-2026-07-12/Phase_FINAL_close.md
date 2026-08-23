---
title: "Phase_FINAL_close — B2-prime proba #5 (op-oblique-beam): KANAL CYRKULACYJNY ZMIERZONY (G5b FAIL wariant (i) przy G6a PASS); P2 w brzmieniu 'slepota na kret' rozstrzygniete NEGATYWNIE; policzalnosc czysto eikonalowa zalamana na geometrii skosnej jako skutek tego samego zjawiska"
date: 2026-07-14
type: phase-final
tgp_owner: research/op-oblique-beam-2026-07-12
status: CLOSED
verdict: "G1 PASS; G2 PASS (0.0000/tau=600); G3 PASS (1.12%); G4 FAIL (monotonia zlamana przez tlumienie alpha na galezi +b serii 1.5 — ten sam kanal odd); G5a FAIL (3/8, wszystkie w serii 1.5); G5b FAIL 3/4 par wg reguly wzmacniajacej WARIANT (i) -> KANAL CYRKULACYJNY ZMIERZONY; G6a PASS (S-inv 0.056 deg = 1.3% baseline; S-conj bitowo); G6b PASS (0.000)"
tags: [gravity-bridge, branch-B2prime, oblique-beam, checkerboard, circulation-channel, discovery, P2-negative, L256, closed]
related:
  - "[[Phase0_balance.md]]"
  - "[[README.md]]"
  - "[[../op-shortwave-lattice-2026-07-05/Phase_FINAL_close.md]]"
  - "[[../op-asymmetric-lattice-2026-07-05/Phase_FINAL_close.md]]"
---

# FINAL close — proba #5: wiazka skosna na szachownicy L=256

## 1. Werdykt i sciezka reguly decyzyjnej

**KANAL CYRKULACYJNY ZMIERZONY** (regula wzmacniajaca G5b,
wariant (i), pre-rejestrowana w LOCKu): FAIL w 3 z 4 par
kryterialnych ze SPOJNYM (ujemnym) znakiem nadwyzki, kazda
nadwyzka wieksza od rozrzutu plaszczyzn Delta_eik, przy G6a PASS
(sygnal jest wlasnoscia pola, nie bledem ramki skosnej).

**P2 w brzmieniu "kanal geometryczny slepy na znak kretu"
rozstrzygniete NEGATYWNIE.** Status P1 (policzalnosc eikonalna
na geometriach chronionych lustrem, proby #2+#3) BEZ ZMIAN.

Sciezka reguly decyzyjnej (LOCK, sekcja 5) — uczciwe
rozstrzygniecie zbiegu wynikow:

- G2 PASS -> cykl wolno bylo kontynuowac (0.0000, powtorka #3).
- G3 PASS -> kanal skosny rzetelny (max 1.12% na modach (2q,q)).
- G6a PASS -> zadnej erraty implementacyjnej; ewaluacja G5 wazna.
- G5b FAIL wg wariantu (i) -> "kanal cyrkulacyjny ZMIERZONY"
  (galaz reguly: odkrycie; osobny op charakteryzacji).
- G4 FAIL i G5a FAIL wystapily ROWNOCZESNIE z G5b FAIL. Litera
  reguly przewiduje dla samodzielnego G5a FAIL "powrot do
  analizy", a kombinacji G5a FAIL + G5b FAIL nie enumeruje.
  Rozstrzygniecie (spisane tu, przed zamknieciem, z pelna
  transparencja): wszystkie trzy FAIL-e maja JEDNA wspolna
  przyczyne — zmierzona skladowa katowa sprzezona ze strona
  przelotu wzgledem kretu (kanal odd). Nadwyzka odd (do -3.25 deg
  przy sigma 0.08) przesuwa pojedyncze punkty +b/-b poza pasmo
  [0.5, 2.0] (G5a) i lamie monotonie |alpha| na galezi +b serii
  1.5 (G4: |alpha(1.5,+8)| = 0.26 < |alpha(1.5,+12)| = 0.61 —
  tlumienie az do odwrocenia znaku). "Policzalnosc czysto
  eikonalowa" zalamana na geometrii skosnej NIE jest wiec drugim,
  niezaleznym zjawiskiem — to ten sam kanal cyrkulacyjny, ktorego
  eikonal (slepy na kret Z DEFINICJI) nie zawiera. Werdykt
  odkrycia opiera sie wylacznie na pre-rejestrowanej regule
  wzmacniajacej G5b; zadne kryterium nie zostalo zmienione,
  zadna para nie zostala wykluczona.

Obserwacja POST HOC (poza kryteriami, oznaczona jako taka,
niekryterialna): czesc PARZYSTA pomiaru,
[alpha(+b)+alpha(-b)]/2 vs [pred(+)+pred(-)]/2, lezy w pasmie
[0.5, 2.0] we WSZYSTKICH 4 parach: 1.1/14: 4.93/3.58 = 1.38;
1.1/20: 1.91/1.46 = 1.30; 1.5/8: 1.96/3.43 = 0.57; 1.5/12:
1.25/1.19 = 1.06. Eikonal opisuje czesc parzysta; zawodzi
wylacznie czesc nieparzysta (odd) — spójne z obrazem "eikonal +
brakujacy czlon cyrkulacyjny".

## 2. Wyniki liczbowe (komplet)

### 2.1 Fazy techniczne

```text
G1: PASS we wszystkich fazach (flagi czyste; E-dryft <= 2e-15;
    gamma_core < 0 po tranzycie; n_tot = 0, 4 rdzenie zawsze;
    determinizm bitowy relaksacji i fali: True)
G2 BRAMKA: max przemieszczenie 0.0000 przez tau <= 600 (per
    rdzen 0.0000; E-dryft bramki 0.00e+00) — powtorka proby #3
    zmierzona, nie cytowana. PASS
G3 (mody skosne (2q,q), pasmo [0.5, 1.7]): bledy {-0.02, -0.08,
    -0.33, -0.67, -1.12}% dla |k| {0.549, 0.713, 1.098, 1.372,
    1.646}; analityka sieciowa przewiduje pomiar do <0.05 p.p.
    PASS (prog 5%)
tlo: szew 0.1244 rad; residuum 9.648e-07; rdzen A (-63.75,-63.75)
    — bitowo jak proba #3
```

### 2.2 Skan kryterialny (8 punktow; alpha > 0 = KU A)

```text
punkt        alpha+/-sigma [deg]   pred(b_eff)  ratio   b_eff/lam  czulosc[xi=10,130]
1.1 b=+14:   +4.391 +/- 0.214      +3.425       1.282   1.77       [+3.240, +3.936]
1.1 b=-14:   +5.473 +/- 0.186      +3.726       1.469   1.74       [+3.524, +3.877]
1.1 b=+20:   +1.386 +/- 0.114      +1.430       0.969   2.59       [+1.334, +1.759]
1.1 b=-20:   +2.425 +/- 0.087      +1.491       1.626   2.55       [+1.387, +1.450]
1.5 b= +8:   +0.258 +/- 0.085      +2.240       0.115   2.33       [+2.157, +2.351]  POZA
1.5 b= -8:   +3.671 +/- 0.055      +4.610       0.796   1.83       [+4.514, +4.628]
1.5 b=+12:   -0.609 +/- 0.076      +0.951      -0.641   3.39       [+0.891, +1.048]  POZA (ZNAK)
1.5 b=-12:   +3.116 +/- 0.025      +1.423       2.189   2.90       [+1.349, +1.434]  POZA
-> G5a: FAIL (3/8; wszystkie trzy w serii 1.5 — najmniejsze
   sigma, najsilniejszy kanal odd wzgledem wielkosci alpha)
```

### 2.3 G5b per para (RDZEN pomiaru; D17)

```text
para     Delta_meas  Delta_eik  reszta   3*sigma_Delta  rozrzut_plaszczyzn  werdykt
1.1/14:  -1.082      -0.301     -0.781   0.852          0.356               OK (zero)
1.1/20:  -1.039      -0.061     -0.978   0.430          0.370               FAIL
1.5/8:   -3.414      -2.369     -1.044   0.302          0.091               FAIL
1.5/12:  -3.725      -0.473     -3.252   0.241          0.086               FAIL
(Delta_eik przy WLASNYM b_eff kazdego runu i wlasnej plaszczyznie
 xi_eval; rozrzut = max-min roznicy pred(+)-pred(-) po wspolnych
 plaszczyznach {xi_eval+, xi_eval-, 10, 130})
regula wzmacniajaca: 3 pary FAIL, znak nadwyzki SPOJNY (ujemny),
 kazda nadwyzka > rozrzut plaszczyzn -> WARIANT (i): ODKRYCIE
obserwacja poza kryteriami (niedolaczalna, D12):
 para 1.1/24: Delta_meas = -0.968, Delta_eik = -0.014, reszta
 = -0.953 vs 3*sigma_Delta = 0.238, rozrzut 0.420 — ten sam znak
 (raport; NIE wchodzi do reguly)
```

### 2.4 Symetrie i liniowosc

```text
G6a-S-inv  (pakiet wsteczny -u, b=-14, vs baseline +14, +u):
    +4.44750 vs +4.39129 -> |d| = 0.0562 deg = 1.3% baseline
    <= prog 0.2196 (5%). PASS (zaokraglenia, nie bitowosc —
    pre-rejestrowane D11a; patrz errata E3)
G6a-S-conj (tlo lustrzane theta -> -theta): max|psi_bg_m -
    conj(psi_bg)| = 0.000e+00; alpha bitowo rowne (+4.39129,
    identyczne probki). PASS
G6b (u0 -> u0/2): roznica 0.000 deg (probki identyczne do 2
    miejsc, alpha do 3 miejsc). PASS — pomiar w rezimie liniowym
```

### 2.5 Sygnatura kanalu (opisowo, w jezyku substratu)

1. **Znak**: alpha(+b) < alpha(-b) we wszystkich 5 zmierzonych
   parach (4 kryterialne + obs 1.1/24). Przy u = (2,1)/sqrt(5)
   i n = (-1,2)/sqrt(5): przelot po stronie +eta rdzenia A (+1)
   daje slabsze skupienie KU A, po stronie -eta silniejsze.
2. **Sygnatura w b_eff (przed jakimkolwiek oknem pomiarowym)**:
   strona +b jest ODPYCHANA (b_eff > b_nom: 8 -> 8.83,
   12 -> 12.87; przy 1.1: 14 -> 13.68, 20 -> 19.96 — dip
   mniejszy niz w probie #3), strona -b PRZYCIAGANA (8 -> 6.94,
   12 -> 10.99, 14 -> 13.42, 20 -> 19.71). Asymetria dipu
   b_eff(+) - b_eff(-) rosnie z omega: ~0.25 (1.1) vs ~1.9 (1.5).
   Czesc tego efektu jest juz POCHLONIETA w Delta_eik (predykcje
   przy wlasnym b_eff, D17) — reszta i tak przekracza 3 sigma.
3. **Skalowanie**: nadwyzka |reszty| rosnie z omega (seria 1.5 >>
   1.1 wzgledem alpha) i z |b| w serii 1.5 (-1.04 przy b=8,
   -3.25 przy b=12) — charakteryzacja ilosciowa = osobny op.
4. **S-inv potwierdza geometryczna nature**: run WSTECZNY po tej
   samej stronie wzgledem cyrkulacji (b=-14, kierunek -u) daje
   alpha rowna baseline (+4.45 vs +4.39), podczas gdy run PRZEDNI
   po przeciwnej stronie (b=-14, +u) daje +5.47. Znak kanalu
   podaza za strona przelotu wzgledem KRETU, nie za znakiem
   wspolrzednej eta — dokladnie zachowanie skladowej sprzezonej
   z cyrkulacja.

## 3. Interpretacja

### CO TO ZNACZY (zmierzone)

- Na tle o niezerowym lokalnym krecie istnieje skladowa kata
  ugiecia pakietu amplitudowego sprzezona ze STRONA przelotu
  wzgledem cyrkulacji defektu — kanal, ktorego eikonal substratu
  n^2(|psi|) nie zawiera z definicji. Wielkosc: reszty -0.78 do
  -3.25 deg na parach (przy 3 sigma_Delta 0.24-0.85 i rozrzucie
  plaszczyzn 0.09-0.37).
- Kanal jest niewidoczny na geometriach z lustrem przez os wiazki
  (proby #2/#3: alpha_odd = 0 z konstrukcji) — dlatego P1 tam
  przeszlo i jego status sie NIE zmienia. Ujawnia sie dopiero bez
  oslony symetrii.
- Czesc parzysta ugiecia pozostaje zgodna z eikonalem (post hoc,
  sekcja 1) — obraz "eikonal + brakujacy czlon sprzezony
  z gradientem fazy tla" (kandydat w rownaniu pakietu: czlon
  ~ 2 i kappa grad(theta_bg).grad(u), spisany juz w LOCKu proby
  #4, sekcja 1.1).

### CZEGO TO NIE ZNACZY (forbidden moves #6, #7)

- To NIE jest deklaracja frame-draggingu / Lense-Thirringa / GR;
  tau != czas; zero deklaracji obserwacyjnych. Analogia
  grawitomagnetyczna pozostaje ANALOGIA jezykowa.
- To NIE falsyfikuje mostu geometrycznego ani P1: policzalnosc
  eikonalna w domenie stosowalnosci na tlach chronionych lustrem
  jest zmierzona (11/12 punktow prob #2+#3) i nietknieta.
- Wyniki S-inv/S-conj NIE sa pomiarem P2 (testy implementacji).
- NIE zmierzylismy tu praw skalowania kanalu (b, omega, |n|,
  konfiguracja kretu) — to zadanie osobnego opa charakteryzacji.
- Obserwacja czesci parzystej i sygnatura b_eff sa POST HOC /
  poza kryteriami — raportowane dla pelni obrazu, nie jako
  kryteria.

## 4. Erraty i transparencja (dzien runow: 2026-07-13/14)

- **E1 (przerwa API po LOCKu):** sesja wykonawcy zostala przerwana
  bledem API (limit wydatkow) 2026-07-13, PO zmianie statusu
  frontmatter na LOCKED, a PRZED dopisaniem sekcji 10 i przed
  pierwsza linijka kodu. Kontynuacja po odblokowaniu: najpierw
  dokonczono LOCK (sekcja 9 + sekcja 10 pre-code), potem dopiero
  kod. Zero wplywu na pomiary (zadnego kodu ani runow przed
  kompletnym LOCKiem).
- **E2 (przestoje miedzy fazami):** trzy przestoje organizacyjne
  (po Phase 1, po Phase 2, po Phase 3) — kazda faza rekonstruuje
  tlo deterministycznie od zera (determinizm bitowy potwierdzony
  w G1), wiec przestoje nie maja wplywu na wyniki.
- **E3 (S-inv na poziomie 0.056 deg, nie bitowo):**
  pre-rejestrowane w D11a: ansatz z obcieciem obrazow k in
  [-2,2] nie jest bitowo inwersyjnie symetryczny (diagnostyka
  max|psi_bg - I[psi_bg]| = 2.348 — wartosc zdominowana przez
  otoczenia rdzeni i wrap fazy przy inwersji indeksow); do tego
  artefakt detektora w b_eff (13.962 wstecz vs 13.683 przod,
  wzorzec ~0.25-0.5 z prob #2/#3). Kryterium 5% spelnione
  z 4x zapasem; roznica 0.056 deg = 0.3 sigma pomiaru.
- **Transparencja:** zadnych zmian progow/parametrow/procedur po
  pierwszym runie; maski okien zamrozone (forbidden move #8);
  zestaw 4 par ZAMKNIETY — obserwacje +/-24 nie weszly do zadnego
  kryterium; zadna para nie zostala wykluczona; Delta_meas nigdzie
  nie jest raportowana bez Delta_eik i rozrzutu plaszczyzn
  (forbidden move #5).

## 5. Nastepne kroki (kazdy za osobna autoryzacja, osobny lock)

1. **op charakteryzacji kanalu cyrkulacyjnego** (deliverable tego
   cyklu jako wejscie): skalowanie nadwyzki z b i omega, zaleznosc
   od konfiguracji kretu (np. defekt |n|=2, pary wspolbiezne),
   test kandydata teoretycznego ~ grad(theta).grad(u) przeciw
   zmierzonym 5 parom; projekt eikonalu ROZSZERZONEGO o czlon
   cyrkulacyjny i powtorka predykcji.
2. Dopiero po charakteryzacji: powrot do P2 w brzmieniu
   ilosciowym ("kanal cyrkulacyjny POLICZALNY?") zamiast zerowego
   ("slepota").

## 6. Pliki cyklu

- `Phase0_balance.md` (LOCKED 2026-07-12; sekcja 10 D1-D18
  pre-code)
- `Phase0_check.py` / `Phase0_output.txt` (geometria pasa,
  dyspersja skosna, tlo, ray tracing wiazacy, budzet)
- `Phase1_gate_and_calibration.py` / `Phase1_output.txt` (G2
  BRAMKA 0.0000/600 PASS; G3 skosne 1.12% PASS)
- `Phase2_deflection.py` / `Phase2_output.txt` (4 pary + G6b +
  obs +/-24; G4 FAIL, G5a FAIL, G5b FAIL wariant (i), G6b PASS)
- `Phase3_symmetry.py` / `Phase3_output.txt` (G6a PASS: S-inv
  0.056 deg, S-conj bitowo)
- `Phase_FINAL_close.md`, `README.md`
