---
title: "Phase0_balance — DRAFT: B2-prime proba #5 — pelny P2: wiazka SKOSNA (atan(1/2)) na zwalidowanej szachownicy L=256 (lustro lamie geometria wiazki, nie tlo)"
date: 2026-07-12
type: phase0-draft
tgp_owner: research/op-oblique-beam-2026-07-12
status: LOCKED
anti_lakatos_lock: ACTIVE (autoryzacja autora: polecenie realizacji cyklu, 2026-07-12; lock przed kodem; wykonawca: osobny agent, nowa sesja — spelnione)
tags: [gravity-bridge, branch-B2prime, oblique-beam, checkerboard, vortex-lensing, circulation-channel, universality, P2, L256]
related:
  - "[[../op-asymmetric-lattice-2026-07-05/Phase_FINAL_close.md]]"
  - "[[../op-shortwave-lattice-2026-07-05/Phase_FINAL_close.md]]"
  - "[[../op-lattice-background-2026-07-05/Phase_FINAL_close.md]]"
---

# Phase 0 — DRAFT: B2-prime proba #5 — P2 wiazka skosna na szachownicy

## 0. Hipoteza robocza i miejsce w programie

Stan lancucha (wszystko zmierzone, nic do wiary):
- P1 ZAMKNIETE (proby #2+#3): 11/12 punktow w pasmie [0.5, 2.0],
  19/19 runow znak KU defektowi; granica eikonalu przy
  b_eff ~ lambda.
- P2 OTWARTE. Proba #4 (op-asymmetric-lattice): STOP na bramce G2
  — siatka skosna dynamicznie niezdatna (mod scinajacy podsieci,
  gamma_s = +0.284/tau, utrata rdzeni tau=227). DWIE lekcje
  pomiarowe: (a) zero sily netto NIE implikuje stabilnosci —
  stabilnosc szachownicy to sztywnosc D4; (b) dryft podsieci
  podczas relaksacji zepsul tez zalozenie Delta_eik ~ 0
  (wiazaco do -1.11 deg).
- Wniosek konstrukcyjny: lustro nalezy zlamac GEOMETRIA WIAZKI,
  nie tla. Tlo wraca na SZACHOWNICE L=256 — jedyne tlo
  z pomiarem 0.0000 przemieszczenia przez tau=600 (proba #3)
  i residuum 9.6e-07.

Konstrukcja: pakiet biegnie pod katem phi = atan(1/2) = 26.565 deg
do osi siatki, u = (2,1)/sqrt(5), n = (-1,2)/sqrt(5). Szachownica
ma grupe luster D4 (osie i przekatne przez rdzenie); kierunek
(2,1) NIE lezy w zadnym z nich -> odbicie wzgledem linii wiazki
(+ sprzezenie) NIE jest symetria ukladu -> alpha_odd niewiazana
konstrukcja. Rozne od prob #2/#3 wylacznie ORIENTACJA wiazki;
tlo, siatka, integrator, metrologia — bez zmian.

Roznica wzgledem proby #4 (spisana z gory): Delta_eik NIE jest
~0 — sasiedzi soczewki rozpraszajacej leza asymetrycznie
wzgledem linii wiazki (odleglosci prostopadle rozne dla +b i -b),
wiec eikonal SLEPY NA KRET przewiduje niezerowa asymetrie.
Obserwabla kanalu cyrkulacyjnego = RESZTA Delta_meas - Delta_eik.
Eksperyment pozostaje zerowy w reszcie, nie w surowej roznicy.

Pytania:

> (P2-a) POLICZALNOSC BEZ OSLONY SYMETRII: czy kat ugiecia
>      w KAZDYM punkcie +/-b (geometria bez zadnego lustra przez
>      os wiazki) zgadza sie z eikonalem substratu
>      (znak zgodny z predykcja, pasmo [0.5, 2.0])?
> (P2-b) KANAL CYRKULACYJNY (RDZEN): czy Delta_meas = alpha(+b)
>      - alpha(-b) jest w calosci wyjasniona przez Delta_eik
>      (slepa na kret geometria soczewek), czy zostaje skladowa
>      sprzezona ze znakiem kretu tla?

G5 PASS (a+b) = P2 rozstrzygniete pozytywnie (slepota na kret
zmierzona tam, gdzie nic jej nie wymusza). G5b FAIL wg reguly
wzmacniajacej = ODKRYCIE kanalu cyrkulacyjnego (osobny op
charakteryzacji; bez deklaracji GR).

## 1. Predykcje algebraiczne do zalockowania (PRZED kodem)

### 1.1 Rownanie #1 — kanal amplitudowy i dyspersja skosna

```text
n^2(x,y) = (omega^2 - V''(|psi_bg|)) / (omega^2 - V''(s*));
A_t = 2.5690; m_TV = 0.9374; eikonal slepy na kret Z DEFINICJI.
k0/lambda/v_g:  omega=1.1: 0.8140 / 7.72 / 0.3700
                omega=1.5: 1.6561 / 3.79 / 0.5520
dyspersja sieciowa dla k = k0*(2,1)/sqrt(5) (kazda skladowa
osobno w omega_lat^2 = kappa*[(2-2cos(kx dx)) + (2-2cos(ky dx))]
/dx^2 + V''): czynnik anizotropii (kx^4+ky^4)/k^4 = 17/25 = 0.68
-> blad przy k0(1.5): ~ -1.2% (axialnie: -1.72%, proba #3);
   przy k0(1.1): ~ -0.13%. Oszacowanie analityczne; MIERZONE w G3
   na modach skosnych k = (2 pi/L)*(2q, q), q in {10,13,20,25,30}
   (k = 0.549, 0.713, 1.098, 1.372, 1.647; pasmo [0.5, 1.7]).
```

### 1.2 Rownanie #2 — geometria wiazki skosnej (fundament cyklu)

```text
rozpraszacz: A = +1 w (-64,-64) (pozycja z psi_ref);
wspolrzedne wiazki (wzgledem A, skladowe wrapowane do
[-L/2, L/2)): xi = (2 dx' + dy')/sqrt(5) [wzdluz wiazki],
eta = (-dx' + 2 dy')/sqrt(5) [poprzecznie];
start pakietu: xi = -36 (4.5 sx przed rdzeniem), eta = b.
obrazy soczewek w pasie pomiarowym (xi, eta) wzgledem A:
  (i,j)=(0,1) [-1]: xi= 57.2, eta=+114.5
  (i,j)=(1,0) [-1]: xi=114.5, eta= -57.2
  (i,j)=(1,1) [+1]: xi=171.7, eta= +57.2
  (i,j)=(2,1) [+1]: xi=286.2, eta=  0.0  <- wiazka trafia w NIA;
     okno pomiaru konczy sie na xi=130 << 286 (poza cyklem)
asymetria srodowiska: najblizsze soczewki poza rzedem leza po
  PRZECIWNYCH stronach linii wiazki w roznych odleglosciach
  (57.2 i 114.5) -> Delta_eik != 0, POLICZALNA (slepa na kret);
  zadne lustro przez os wiazki nie istnieje -> alpha_odd wolna.
symetrie POZOSTALE (testy implementacji G6a):
  S-inv: szachownica jest inwersyjnie symetryczna wokol kazdego
    rdzenia -> run(+b, kierunek +u) = run(-b, kierunek -u)
    co do zaokraglen (pakiet WSTECZNY, wylacznie G6a).
  S-conj: tlo lustrzane (theta -> -theta) -> run bitowo rowny
    sprzezeniu runu bazowego (wzorzec D11 prob #2/#3).
```

### 1.3 Rownanie #3 — zaklad eikonalny i skala efektu

Formula ogonowa (kalibracja: trace proby #3, blad 3%):
t_deg(d) = 698/d^2 (omega=1.1); 169/d^2 (omega=1.5).

```text
orientacyjnie: wklad rozpraszacza jak w probie #3 (alpha(14) ~
  3.8; alpha(20) ~ 1.8 przy 1.1; alpha(8) ~ 3.1, alpha(12) ~ 1.2
  przy 1.5); wklady sasiadow pasa: soczewka (1,0) w eta=-57.2:
  odleglosc od wiazki 57.2+b (przy +b) vs 57.2-b (przy -b);
  soczewka (0,1) w eta=+114.5: 114.5-b vs 114.5+b ->
  Delta_eik ~ 2[t(57-b) - t(57+b)] - 2[t(114-b) - t(114+b)]
  orientacyjnie (1.1): b=14: ~ +0.4 deg; b=20: ~ +0.7 deg
                (1.5): b=8:  ~ +0.05 deg; b=12: ~ +0.09 deg
  (znaki i wielkosci ORIENTACYJNE — geometria skosna z oknem
  skonczonym; WIAZACO liczy Phase0_check: ray tracing na
  zrelaksowanym tle, OSOBNO +b i -b, plaszczyzny xi_eval,
  czulosc na krancach okna RAPORTOWANA OBOWIAZKOWO obok Delta_eik)
sigma (proba #3, ten sam aparat): 0.20 (1.1) do 0.014 (1.5,16);
  3 sigma_Delta ~ 0.83 (1.1/14) do ~0.15 (1.5). Najczystsze zero:
  seria 1.5; seria 1.1 niesie wiekszy Delta_eik (test dodatni
  policzalnosci asymetrii + slabszy test zera).
```

### 1.4 Budzet czasu (KRYTERIUM)

```text
start xi=-36; prog wyzwolenia xi=+10; przelot do 5. probki:
  omega=1.1: ~tau=140; omega=1.5: ~tau=99 (jak proba #3)
bramka G2: tlo bez pakietu, II rzad, okno tau=600; kryterium:
  max przemieszczenie kazdego z 4 rdzeni <= 1.0 na tau in [0,450]
  (tlo IDENTYCZNE z proba #3 — oczekiwanie: powtorka 0.0000,
  przy identycznym kodzie relaksacji mozliwa bitowa rownosc;
  MIERZONE ponownie, nie cytowane)
```

## 2. Model i parametry (propozycja do LOCKa)

```text
substrat: bez zmian (kappa=0.5, a=0.5, b=1.6, c=1.0; psi in C;
  N=512, L=256, dx=0.5, periodycznie, eps=0.30, S_GUARD=10 flaga,
  seed=20260704)
tlo: SZACHOWNICA 4 wirow +/-64 — DOKLADNIE jak proba #3
  (2 pary poziome theta_h_pair, amplituda = iloczyn 4 profili,
  relaksacja 2000 dt=0.02, diagnostyka szwu; detektor 2+2,
  nominaly +1: (-64,-64),(+64,+64); -1: (+64,-64),(-64,+64))
sektor falowy, referencja blizniacza: bez zmian (D1, D3 proby #3)
pakiet SKOSNY: u = (2,1)/sqrt(5); n = (-1,2)/sqrt(5);
  r0 = A_nom + (-36)*u + b*n; env izotropowa exp(-|r-r0|^2/(2*8^2))
  (niezmiennicza na obrot — bez zmian ksztaltu);
  delta_psi(0) = e^{i theta_bg} * u0 * env * cos(k0 * u.(r-r0));
  delta_psi(-dt): faza + omega dt (pakiet biezacy w +u); u0=1e-3
wspolrzedne pomiaru: xi, eta wzgledem BIEZACEJ pozycji A
  z psi_ref (skladowe wrapowane); wyzwolenie CENTROIDEM:
  xi_c(|delta_psi|^2) > +10; w chwili wyzwolenia okno ZAMRAZANE:
  pas {+10 < xi < +130} MINUS dyski r < 12 wokol kazdego z 4
  rdzeni (w pasie leza obrazy soczewek (0,1) i (1,0) — wylaczenie
  pre-rejestrowane; wklad wagowany |delta_psi|^2 i tak maly);
  5 probek co Delta_tau = 4; P = -Re sum(conj(delta_psi_tau)
  grad(delta_psi)) dx^2 w oknie; usrednienie WEKTOROWE;
  alpha := -sign(b) * atan2(P.n, P.u) [deg] (alpha > 0 = KU A)
b_eff: min odleglosc periodyczna centroidu od A (do wyzwolenia);
  b_eff/lambda raportowane per run
alpha_pred: ray tracing na TYM SAMYM zrelaksowanym tle; start
  r0(b_eff) = A + (-36)u + sign(b) b_eff n, kierunek u; kat czytany
  w plaszczyznie xi_eval = xi centroidu 3. probki; czulosc:
  xi = +10 i +130; przechwyt r < 0.6
skan glowny (KRYTERIALNY, zamkniety) — PARY +/-b:
  omega=1.1 x b in {+/-14, +/-20}; omega=1.5 x b in {+/-8, +/-12}
  (8 runow; wszystkie b >= 1.5 lambda + 0.9; odleglosc wiazki od
  soczewek pasa >= 33 >= 1.5 lambda w kazdym punkcie)
Delta per para (D18 proby #4): Delta_meas = alpha(+b)-alpha(-b);
  Delta_eik = alpha_pred(b_eff+) - alpha_pred(b_eff-) (kazda przy
  WLASNYM b_eff); sigma_Delta = sqrt(sigma+^2 + sigma-^2)
symetrie (G6a): S-inv: (1.1,-14) pakiet WSTECZNY (-u, start
  xi=+36 po drugiej stronie) vs (1.1,+14); S-conj: tlo lustrzane
liniowosc (G6b): u0/2 przy (1.1,+14)
obserwacje POZA kryteriami: (1.1, +/-24)
metrologia: jak proba #3 (D8): E-dryft <= 1e-4; gamma_core <= 0.01
  (okna r<5 wokol 4 rdzeni); wrap-guard 0.9L; determinizm bitowy
budzet krokow: relaksacja 2000; bramka 12000 (tau=600);
  runy falowe 6000 (tau=300)
koszt: bramka ~15 min; 11 runow Phase2 + 3 runy Phase3 ~ 1-1.5 h
```

## 3. Scenariusze

1. `background_check + LIFETIME GATE (G2)`: identycznie jak
   proba #3 (to samo tlo); MIERZONE ponownie. FAIL -> STOP.
2. `dispersion_TV (G3, mody SKOSNE)`: zaburzenie cos(k.r),
   k = (2 pi/L)(2q, q), q in {10,13,20,25,30}; blad <= 5%.
3. `deflection(+/-b, omega) (G4, G5a, G5b)`: 4 pary; b_eff;
   alpha vs alpha_pred; Delta_meas vs Delta_eik per para.
4. `symmetry (G6a)`: S-inv (pakiet wsteczny), S-conj (lustro).
5. `linearity (G6b)`: u0/2 przy (1.1,+14).
6. obserwacje poza kryteriami: (1.1, +/-24).

## 4. Kryteria (propozycja do LOCKa)

- **G1 techniczne**: jak proba #3 (determinizm bitowy; NaN/
  S_GUARD; E-dryft <= 1e-4; gamma_core <= 0.01; wrap-guard;
  kret tla: 4 rdzenie, n_tot = 0 w kazdym runie).
- **G2 BRAMKA**: max przemieszczenie kazdego rdzenia <= 1.0 dla
  tau <= 450 (okno 600). Oczekiwana powtorka wyniku proby #3;
  MIERZONE. FAIL -> STOP (bez oslabiania).
- **G3 kalibracja SKOSNA**: blad <= 5% dla k in [0.5, 1.7] na
  modach (2q, q). FALSYFIKATOR: anizotropia siatki psuje kanal
  skosny -> stop przed interpretacja.
- **G4 istnienie**: |alpha(1.1,+14)| > 3*sigma ORAZ |alpha|
  maleje z |b| PO OBU STRONACH obu serii (|14|>|20| przy 1.1;
  |8|>|12| przy 1.5).
- **G5a — policzalnosc bez oslony symetrii (P2-a)**: dla
  WSZYSTKICH OSMIU punktow: sign(alpha) = sign(alpha_pred) ORAZ
  alpha/alpha_pred(b_eff) in [0.5, 2.0].
- **G5b RDZEN — zero kanalu cyrkulacyjnego (P2-b)**: dla KAZDEJ
  z 4 par: |Delta_meas - Delta_eik| <= 3*sigma_Delta.
  Regula wzmacniajaca (pre-rejestrowana): (i) FAIL w >= 2 parach
  ze SPOJNYM znakiem nadwyzki ORAZ nadwyzka wieksza od rozrzutu
  plaszczyzn Delta_eik (czulosc D6) = "kanal cyrkulacyjny
  ZMIERZONY"; (ii) FAIL w 1 parze lub nadwyzka wewnatrz rozrzutu
  plaszczyzn = "anomalia/nierozstrzygniete" (raport, bez
  deklaracji odkrycia).
- **G6 symetrie i liniowosc**: (a) S-inv: |alpha(+14,+u) -
  alpha(-14,-u)| <= 0.05*|alpha(+14)|; S-conj: bitowo / ta sama
  tolerancja; (b) u0/2: roznica <= 10%.

## 5. Regula decyzyjna (propozycja)

```text
G2 FAIL -> STOP (nieoczekiwane: to samo tlo co proba #3;
  wskazywaloby na blad srodowiska — errata techniczna).
G5a PASS + G5b PASS (+ G4, G6) -> P2 ROZSTRZYGNIETE POZYTYWNIE:
  kanal geometryczny slepy na znak kretu bez oslony symetrii;
  most ma zmierzone: znak + policzalnosc (P1) + uniwersalnosc
  (P2). Nastepne opy (osobne locki): (i) skalowanie z |n|,
  (ii) granica newtonowska 2D + G_eff.
G5a PASS + G5b FAIL (regula wzmacniajaca, wariant i)
  -> ODKRYCIE kanalu cyrkulacyjnego; deliverable: nadwyzki na
  4 parach; osobny op charakteryzacji (skalowanie z b, omega).
G5b wariant (ii) -> nierozstrzygniete; kandydat: powtorka
  z mniejszym sigma (dluzsze usrednianie — osobny lock).
G5a FAIL (wielkosc) -> policzalnosc geometrii skosnej zalamana;
  powrot do analizy (dyspersja skosna? okno pasowe? — raport).
G5a FAIL (znak) -> anomalia kierunkowa; powrot do analizy.
G3/G4 FAIL -> stop przed interpretacja.
G6a FAIL -> blad implementacji ramki skosnej: errata techniczna
  (dzien runow), dopiero potem ewaluacja G5.
```

## 6. Forbidden moves (propozycja)

1. Zero pinningu i wiezow; tlo bez zmian wzgledem proby #3.
2. Bramki G2 nie wolno oslabic; przy FAIL zero runow ugiecia.
3. Zadnych zmian progow/parametrow/procedur po pierwszym runie
   (wyjatki pre-rejestrowane: jedno obnizenie dt; jeden kontrolny
   run dx=0.25).
4. Zestaw kryterialny (4 pary) ZAMKNIETY; obserwacje +/-24
   niedolaczalne; zadnego wykluczania par post hoc.
5. Delta_meas NIGDY bez Delta_eik i bez rozrzutu plaszczyzn obok.
6. Ewentualnego odkrycia nie interpretowac jako frame-dragging/
   Lense-Thirring GR; tau != czas; bez promocji do core/.
7. Wynik S-inv/S-conj nie jest pomiarem P2 (testy implementacji).
8. Maski okna (pas xi + dyski wykluczenia r<12) ZAMROZONE
   w chwili wyzwolenia; zadnych zmian ksztaltu okna po pierwszym
   runie.

## 7. Czego ten cykl NIE robi

- NIE mierzy skalowania z |n| ani granicy newtonowskiej.
- NIE charakteryzuje kanalu cyrkulacyjnego (jesli wykryty).
- NIE bada rezimu dyfrakcyjnego (wszystkie punkty w domenie).
- Zero GR; tau != czas.

## 8. Deliverables

- `Phase0_balance.md` (LOCKED + sekcja 10 pre-code PRZED kodem;
  wzorzec: sekcja 10 proby #3 D1-D16 + D17/D18 proby #4, ze
  zmianami: ramka skosna xi/eta (definicje, wrapowanie, maski),
  pakiet wsteczny dla S-inv, mody skosne G3, dyski wykluczenia)
- `Phase0_check.py`/`Phase0_output.txt` (geometria pasa: tabela
  obrazow soczewek (xi, eta) i odleglosci od wiazki dla kazdego
  +/-b; dyspersja skosna analitycznie; tlo zrelaksowane;
  ray tracing WIAZACY alpha_pred(+/-b) + Delta_eik + czulosc
  plaszczyzn; tabela projektowa b/lambda; budzet)
- `Phase1_gate_and_calibration.py`/`Phase1_output.txt` (G2 +
  G3 skosne + background_check)
- `Phase2_deflection.py`/`Phase2_output.txt` (G4, G5a, G5b, G6b
  + obserwacje +/-24)
- `Phase3_symmetry.py`/`Phase3_output.txt` (G6a: S-inv, S-conj)
- `Phase_FINAL_close.md`, `README.md`

## 9. Status

**LOCKED (2026-07-12).** Autoryzacja autora: polecenie realizacji
cyklu (2026-07-12). Wykonawca: OSOBNY AGENT w nowej sesji
(spelnione). Obowiazuje pelny anty-Lakatos. Tlo wraca na jedyna
zwalidowana konstrukcje (szachownica L=256: 0.0000 przemieszczenia
przez tau=600, residuum 9.6e-07 — MIERZONE ponownie w G2, nie
cytowane); lustro lamie ORIENTACJA WIAZKI (atan(1/2) poza
wszystkimi lustrami D4 — alpha_odd niewiazana konstrukcja);
Delta_eik != 0 policzalna i wiazana ray tracingiem z obowiazkowa
czuloscia plaszczyzn; kryterium G5b na reszcie Delta_meas -
Delta_eik z pre-rejestrowana regula wzmacniajaca odrozniajaca
odkrycie od anomalii. Zestaw kryterialny 4 pary +/-b ZAMKNIETY,
w calosci w domenie eikonalu. Doprecyzowania pre-code: sekcja 10
(spisana PRZED pierwsza linijka kodu; wzorzec: sekcja 10 proby #3
D1-D16 + D17/D18 proby #4; nowe: pelna definicja ramki xi/eta
z wrapowaniem wzgledem A, maska pasa z dyskami wykluczenia r<12
zamrazana w wyzwoleniu, pakiet wsteczny dla S-inv, mody skosne
(2q,q) w G3, plaszczyzny xi_eval w ray tracingu).

## 10. Doprecyzowania pre-code (spisane przed pierwsza linijka kodu)

Zadne z ponizszych nie zmienia progow G1-G6, modelu ani predykcji —
to ujednoznacznienie procedur. D1, D3, D5, D7-D10, D13-D15
przeniesione WPROST z sekcji 10 proby #3 (op-shortwave-lattice,
D1-D16; identyczne tlo, grid i metrologia), D16 (pakiet
niestandardowy) i D17 (definicje Delta) wg wzorca D17/D18 proby #4
(op-asymmetric-lattice); nowe doprecyzowania wymagane w sekcji 9
draftu oznaczone [NOWE].

- **D1 (integrator):** identycznie jak proby #2/#3 — pelny sektor
  zespolony II rzedu, leapfrog pozycyjny `psi_{n+1} = 2 psi_n -
  psi_{n-1} + dt^2 * RHS(psi_n)`, `RHS(psi) = kappa*Lap(psi) -
  psi*(a - b|psi| + c|psi|^2)`; dt = 0.05. `psi_tau` do obserwabli
  = roznica CENTRALNA `(psi_{n+1} - psi_{n-1})/(2 dt)`.
- **D2 (pakiet SKOSNY, start analityczny) [NOWE]:** wersory:
  u = (2,1)/sqrt(5), n = (-1,2)/sqrt(5); A_nom = (-64,-64);
  r0 = A_nom + (-36)*u + b_nom*n (b_nom ZE ZNAKIEM; b<0 = pakiet
  po stronie -eta linii wiazki); env IZOTROPOWA:
  env = exp(-|r - r0|^2 / (2 * 8^2)) (niezmiennicza na obrot —
  bez zmiany ksztaltu vs proba #3); s := u.(r - r0) = skladowa
  wzdluz wiazki; `psi(0) = psi_bg + e^{i theta_bg} * u0 * env *
  cos(k0 s)`, `psi(-dt) = psi_bg + e^{i theta_bg} * u0 * env *
  cos(k0 s + omega dt)` (pakiet biezacy w +u — dokladny analog
  D2 proby #3 po obrocie); u0 = 1e-3;
  k0 = sqrt((omega^2 - V''(s*)) / kappa).
- **D3 (referencja blizniacza):** identycznie: delta_psi = psi -
  psi_ref, psi_ref ewoluuje TYM SAMYM integratorem z psi_bg bez
  pakietu, w lockstep. Pozycje rdzeni do ramki, b_eff, okna
  i gamma_core czytane z psi_ref (detektor 2+2, D14). To NIE
  jest pinning.
- **D4 (ramka skosna, wyzwolenie, maska okna) [NOWE]:** ramka
  (xi, eta) wzgledem BIEZACEJ pozycji rdzenia A (+1, nominal
  (-64,-64)) z psi_ref: dla punktu (x, y):
  dx' = (x - Ax + L/2) mod L - L/2, dy' = (y - Ay + L/2) mod L -
  L/2 (skladowe wrapowane do [-L/2, L/2));
  xi = (2*dx' + dy')/sqrt(5), eta = (-dx' + 2*dy')/sqrt(5).
  Maska pasa na torusie liczona z TYCH SAMYCH wrapowanych
  skladowych co centroid. Centroid (x_c, y_c) z |delta_psi|^2
  jak w probie #3 (srednia arytmetyczna; zasieg cyklu nie
  przecina szwu dla zadnego runu — wrap-guard aktywny);
  xi_c, eta_c = ramka zastosowana do (x_c, y_c). WYZWOLENIE:
  xi_c > +10. W chwili wyzwolenia maska okna ZAMRAZANA:
  okno = pas {10 < xi < 130} MINUS dyski odleglosci periodycznej
  < 12 od KAZDEGO z 4 rdzeni (pozycje rdzeni z psi_ref w chwili
  wyzwolenia; ramka wzgledem A rowniez zamrozona). Zadnych zmian
  ksztaltu okna po pierwszym runie (forbidden move #8). Potem
  5 probek co Delta_tau = 4: P = -Re sum(conj(delta_psi_tau) *
  grad(delta_psi)) dx^2 w zamrozonym oknie; usrednienie
  WEKTOROWE (Px, Py); rzuty: P.u = (2*Px + Py)/sqrt(5),
  P.n = (-Px + 2*Py)/sqrt(5); `alpha := -sign(b_nom) *
  atan2(P.n, P.u)` [deg] — alpha > 0 = KU rdzeniowi A dla OBU
  znakow b; sigma_alpha = odch. std probek katowych.
- **D5 (b_eff):** min po tau (probki co 10 krokow, do chwili
  wyzwolenia) odleglosci periodycznej centroidu |delta_psi|^2 od
  rdzenia A (pozycja z psi_ref). alpha_pred(b_eff) = ray tracing
  D6 przy b = b_eff (start z sign(b_nom)*b_eff). b_eff/lambda
  raportowane per run; zadnego wykluczania punktow post hoc.
- **D6 (eikonal na rzeczywistym tle, plaszczyzny xi) [NOWE]:**
  n^2 = (omega^2 - V''(|psi_bg|)) / (omega^2 - V''(s*)) na PELNYM
  zrelaksowanym tle szachownicy L=256; bilinearna interpolacja
  n^2 i gradientu; ray tracing hamiltonowski: dr/dsigma = p,
  dp/dsigma = grad(n^2)/2, RK4, dsigma = 0.02; start
  r_start = A + (-36)*u + sign(b_nom)*b_eff*n (A = pozycja
  rdzenia A ze zrelaksowanego tla, detektor D14 — tlo sie nie
  przemieszcza, wiec rowna biezacej z psi_ref), kierunek
  poczatkowy u, |p(start)| = n(r_start). Kat czytany
  w plaszczyznie xi_eval = xi centroidu pakietu w 3. (srodkowej)
  probce: xi promienia (ramka D4 wzgledem tego samego A) rosnie
  monotonicznie; przy przekroczeniu plaszczyzny interpolacja
  liniowa p; alpha_pred = -sign(b_nom) * atan2(p.n, p.u) [deg]
  (p.u = (2 px + py)/sqrt(5), p.n = (-px + 2 py)/sqrt(5)).
  CZULOSC raportowana OBOWIAZKOWO: kat rowniez na krancach okna
  xi = +10 i xi = +130 (rozrzut plaszczyzn = odniesienie reguly
  wzmacniajacej G5b). Phase0_check: czulosc plaszczyzn xi_eval
  in {10, 40, 70, 100, 130} dla +b i -b OSOBNO na kazdym punkcie.
  Przechwyt: odleglosc periodyczna < 0.6 od ktoregokolwiek
  z 4 rdzeni -> NaN (raport).
- **D7 (G3, dyspersja na MODACH SKOSNYCH) [NOWE]:** tlo
  jednorodne psi = s*; zaburzenie RZECZYWISTE 1e-3 * cos(k.r),
  k = (2 pi / L) * (2q, q), q in {10, 13, 20, 25, 30}
  (|k| = {0.549, 0.713, 1.098, 1.372, 1.647} — pasmo [0.5, 1.7];
  mod periodyczny na torusie z konstrukcji); u_tau(0) = 0;
  omega_meas ze zliczania przejsc przez zero Re(psi - s*)
  w punkcie probki (N/2, 10), okno tau = 100; omega_teor =
  sqrt(kappa |k|^2 + V''(s*)); dyspersja sieciowa analityczna
  raportowana obok: omega_lat = sqrt(kappa*[(2-2cos(kx dx)) +
  (2-2cos(ky dx))]/dx^2 + V''(s*)) (anizotropia
  (kx^4+ky^4)/|k|^4 = 17/25); prog G3: 5%; determinizm bitowy
  pierwszego modu.
- **D8 (G1, metryki):** jak proba #3: dryft sekularny energii =
  |<H> ostatnia cwiartka - <H> pierwsza cwiartka| / |H(0)|, prog
  1e-4; oscylacja leapfroga raportowana osobno. gamma_core =
  slope fitu ln RMS |delta_psi| w SUMIE okien r < 5 wokol
  WSZYSTKICH CZTERECH rdzeni (pozycje z psi_ref), fit po
  tranzycie (od wyzwolenia), prog 0.01. Kret tla: n_tot(psi_ref,
  koniec) = 0 ORAZ 4 rdzenie obecne (2 klastry na znak, D14).
  Wrap-guard: |przemieszczenie centroidu od startu| < 0.9 L =
  230.4. S_GUARD = 10 flaga (nie clamp). Determinizm bitowy
  w Phase 1 (relaksacja + run falowy). Zapisy co 10 krokow;
  seed = 20260704.
- **D9 (budzet krokow):** relaksacja tla: 2000 krokow dt = 0.02.
  Bramka G2: 12000 krokow (tau = 600) bez pakietu. Runy falowe:
  limit 6000 krokow (tau = 300). BEZ wyjatku superpozycji (D16
  proby #3 nie wystepuje). Brak 5 probek = run niewazny do
  kryterium, przyczyna zapisana.
- **D10 (G6b, liniowosc):** |alpha(u0) - alpha(u0/2)| <= 0.1 *
  max(|alpha(u0)|, |alpha(u0/2)|) przy (omega = 1.1, b = +14).
- **D11 (G6a — pre-rejestracja procedur symetrii) [NOWE]:**
  (a) **S-inv:** tlo TO SAMO psi_bg (inwersja wokol rdzenia A,
  r -> 2A - r, jest dokladna symetria szachownicy z zachowaniem
  ladunkow — inwersja w 2D zachowuje kret); run testowy =
  (omega=1.1, b_nom=-14) z pakietem WSTECZNYM wg D16; baseline =
  (1.1, +14) pakiet standardowy D2 (powtorka runu Phase 2).
  Kryterium: |alpha(+14, +u) - alpha(-14, -u)| <= 0.05 *
  |alpha(+14)|. Oczekiwane zaokraglenia, nie bitowosc (ansatz
  z obcieciem obrazow k in [-2,2] nie jest bitowo inwersyjnie
  symetryczny; diagnostyka max|psi_bg - I[psi_bg]| raportowana,
  I = inwersja indeksow wokol nominalnego A).
  (b) **S-conj:** tlo lustrzane theta -> -theta (negacja fazy
  ansatzu przed relaksacja; wzorzec D11 prob #2/#3): wszystkie
  wspolczynniki rzeczywiste, pakiet w lokalnej ramce fazy SWOJEGO
  tla -> run lustrzany jest DOKLADNIE sprzezeniem zespolonym runu
  bazowego -> oczekiwanie: bitowo (kontrola max|psi_bg_m -
  conj(psi_bg)|), w przeciwnym razie rownosc alpha w tolerancji
  5%|alpha|; detektor na tle lustrzanym: znaki nominalow
  odwrocone. Wynik G6a nie jest pomiarem P2 (forbidden move #7).
- **D12 (domena eikonalu i geometria pasa):** wszystkie 8 punktow
  kryterialnych spelnia b_nom >= 1.5*lambda + 0.9 (progi: 1.1:
  12.48; 1.5: 6.59) ORAZ odleglosc linii wiazki (eta = b) od
  obrazow soczewek pasa ((1,0): eta=-57.2; (0,1): eta=+114.5) >=
  33 >= 1.5*lambda w kazdym punkcie (tabela w Phase0_check).
  Obraz (2,1) na osi wiazki lezy przy xi=286.2 — poza oknem
  (130) i poza limitem probek. Obserwacje (1.1, +/-24) POZA
  kryteriami, niedolaczalne (forbidden move #4); zadnego
  wykluczania par post hoc, nawet gdyby b_eff/lambda < 1.5.
- **D13 (narzedzia):** zero `cd` w powloce; sciezki od korzenia
  vaulta; kazdy skrypt pisze output przez `python -u ... >` do
  `PhaseX_output.txt` w katalogu opa; dlugie obliczenia w tle
  i CZEKAC na wynik przed interpretacja.
- **D14 (detektor rdzeni 2+2):** identycznie jak proba #3:
  winding plakietowy W; komorki z W*sign > 0.5 klastrowane po
  odleglosci periodycznej (promien <= 8 komorek od komorki
  bazowej); wymagane DOKLADNIE 2 klastry na znak (inaczej flaga
  core_lost); pozycja klastra = centroid deficytu Phi (wagi
  max(0, Phi*/2 - Phi)) w oknie 13x13; przypisanie do nominalow
  (+1: (-64,-64), (+64,+64); -1: (+64,-64), (-64,+64)) po
  najmniejszej sumie odleglosci periodycznych. Klucz (+1, 0) =
  rdzen A (rozpraszajacy).
- **D15 (bramka G2 — procedura):** start z zrelaksowanego tla,
  psi(-dt) = psi(0) (zerowa predkosc), dynamika II rzedu D1,
  12000 krokow (tau = 600) BEZ pakietu; pozycje 4 rdzeni (D14) co
  10 krokow; kryterium: max po rdzeniach i po tau <= 450
  odleglosci periodycznej od pozycji poczatkowej <= 1.0;
  (450, 600] raportowane poza kryterium. gamma_s: fit liniowy
  ln(max_disp) vs tau na odcinku max_disp in [0.05, 1.0] (jesli
  osiagalny). Oczekiwana powtorka 0.0000 proby #3 (to samo tlo)
  — MIERZONE, nie cytowane. Bramki nie wolno oslabic, skrocic
  ani zreinterpretowac; G2 FAIL -> STOP calego cyklu (errata
  techniczna — wskazywaloby na blad srodowiska).
- **D16 (pakiet WSTECZNY, wylacznie G6a-S-inv) [NOWE]:** obraz
  inwersyjny pakietu D2 wzgledem A_nom: dla runu (1.1, b_nom=-14)
  start r0L = A_nom + (+36)*u + b_nom*n = A_nom + 36u - 14n
  (= dokladna inwersja startu baseline r0 = A_nom - 36u + 14n
  wokol A_nom); env izotropowa jak D2 wokol r0L;
  s := u.(r - r0L); `psi(0) = psi_bg + e^{i theta_bg} * u0 * env *
  cos(k0 s)`, `psi(-dt) = psi_bg + e^{i theta_bg} * u0 * env *
  cos(k0 s - omega dt)` — faza -omega*dt: delta ~ cos(k0 s +
  omega tau) biegnie w -u (dla pakietu D2 w +u faza jest
  +omega*dt). Wyzwolenie: xi_c < -10; w chwili wyzwolenia maska
  ZAMRAZANA: pas {-130 < xi < -10} MINUS dyski odleglosci
  periodycznej < 12 od kazdego z 4 rdzeni (lustro inwersyjne
  maski D4). 5 probek co Delta_tau = 4, usrednienie wektorowe;
  `alpha := -sign(b_nom) * atan2(P.n, -P.u)` [deg] (skladowa
  podluzna -P.u > 0 dla ruchu w -u; alpha > 0 = KU rdzeniowi A —
  inwersja odwraca P, wiec definicja odwzorowuje sie na D4);
  b_eff jak D5 (min odleglosc periodyczna centroidu od A do
  wyzwolenia). Pakiet wsteczny NIE wystepuje w zadnym runie
  kryterialnym (forbidden move #7).
- **D17 (definicje G5b, per para) [NOWE]:** Delta_meas :=
  alpha(+b) - alpha(-b) (oba runy pakietem standardowym D2,
  Phase 2); Delta_eik := alpha_pred(b_eff(+b)) -
  alpha_pred(b_eff(-b)) — kazda predykcja przy WLASNYM b_eff
  swojego runu i wlasnej plaszczyznie ewaluacji (xi centroidu
  3. probki swojego runu; artefakt detektora ~0.25-0.5 w b_eff
  sie odejmuje — wzorzec D18 proby #4); sigma_Delta :=
  sqrt(sigma_+^2 + sigma_-^2). Kryterium G5b: |Delta_meas -
  Delta_eik| <= 3*sigma_Delta dla KAZDEJ z 4 par. Rozrzut
  plaszczyzn Delta_eik (odniesienie reguly wzmacniajacej) :=
  max - min po plaszczyznach {xi_eval, 10, 130} roznicy
  alpha_pred(+) - alpha_pred(-) liczonej plaszczyzna po
  plaszczyznie. Regula wzmacniajaca (sekcja 4): (i) FAIL w >= 2
  parach ze SPOJNYM znakiem nadwyzki ORAZ nadwyzka wieksza od
  rozrzutu plaszczyzn = odkrycie; (ii) FAIL w 1 parze lub
  nadwyzka wewnatrz rozrzutu plaszczyzn = anomalia/
  nierozstrzygniete. Delta_meas NIGDY nie jest raportowana bez
  Delta_eik I rozrzutu plaszczyzn obok (forbidden move #5).
- **D18 (kolejnosc realizacji = protokol):** Phase0_check
  (rachunki wiazace) -> Phase 1 (background_check + BRAMKA G2 +
  G3 skosne; G2 FAIL -> STOP calego cyklu, zero runow ugiecia) ->
  Phase 2 (skan 4 par + G6b + obserwacje +/-24) -> Phase 3 (G6a:
  S-inv, S-conj) -> FINAL close. Interpretacja kazdej fazy
  wylacznie po komplecie jej wynikow; erraty w dniu runow,
  przed interpretacja.
