---
title: "Phase0_balance — kanal dylatacyjny (amplitudowy) na tle Phi*: czy substrat ma kanal UNIWERSALNIE przyciagajacy? (most do grawitacji, galaz C)"
date: 2026-07-05
type: phase0-balance
tgp_owner: research/op-scalar-sector-phistar-2026-07-05
status: LOCKED
anti_lakatos_lock: LOCKED (2026-07-05, autoryzacja autora — wybor drogi "sektor skalarny na tle Phi*" po zamknieciu B1 i B2; lock przed pierwsza linijka kodu)
tags: [gravity-bridge, branch-C, dilaton, amplitude-channel, universality, vortex, quadrupole]
related:
  - "[[../op-goldstone-mediator-2026-07-04/Phase_FINAL_close.md]]"
  - "[[../op-phi-metric-refraction-2026-07-04/Phase_FINAL_close.md]]"
  - "[[../op-lock-interaction-gravity-2026-07-04/Phase_FINAL_close.md]]"
---

# Phase 0 — kanal dylatacyjny na tle Phi* (galaz C)

## 0. Hipoteza robocza

Bilans po B1/B2: kanal fazowy U(1) ma dlugi zasieg (2D-Newton,
C2 = 2*pi*kappa*Phi*), ale jest LADUNKOWY; kanal geometryczny czeka
na tlo stacjonarne. Do grawitacji brakuje UNIWERSALNOSCI ZNAKU.
Standardowy argument teorii pola: wymiana SKALARNA przyciaga
ladunki jednoimienne. Substrat ma naturalny sektor skalarny:
fluktuacje amplitudy (dylatacyjne) `delta|psi|` wokol `s*`.

> Czy kanal amplitudowy na tle Phi* jest UNIWERSALNIE przyciagajacy
> (znak niezalezny od ladunkow fazowych zrodel) i policzalny
> z substratu? Jesli tak — substrat MA kanal o dobrym znaku
> grawitacji, a pytaniem programu staje sie ZASIEG (lustrzane
> odbicie wyniku B1: dobry zasieg, zly znak).

Zrodla = wiry (jedyne trwale obiekty na tle Phi*); ich rdzenie sa
studniami amplitudy sprzezonymi do sektora dylatacyjnego.

## 1. Predykcje algebraiczne (PRZED kodem)

### 1.1 Rownanie #1 — ogon amplitudowy wiru (algebraiczny, nie Yukawa)

Zrodlem `delta|psi|` jest gradient fazy `|grad theta|^2 ~ 1/r^2`:

```text
f(r) = s* (1 - xi^2 / r^2) + O(r^-4),   xi^2 = kappa / V''(s*) = 0.5690
```

Zgodnosc lancucha (JUZ zmierzona w B1, cytowana jako kotwica):
f_ODE(8)/s* = 0.9905 vs 1 - xi^2/64 = 0.9911.
Wniosek strukturalny zapisany z gory: kanal amplitudowy miedzy wirami
NIE musi byc swobodna Yukawa (m_a = sqrt(V''(s*)/kappa) = 1.3256) —
ogony sa algebraiczne, bo sektor jest ZRODLOWANY faza.

### 1.2 Rownanie #2 — dekompozycja kanalow (RDZEN konstrukCJI)

Dwie konfiguracje neutralne o IDENTYCZNYCH pozycjach 4 wirow,
rozniace sie wylacznie znakami (grupy w x = -/+X, pary wewnatrz
grupy wzdluz y, odleglosc d):

```text
pozycje: L1=(-X,+d/2), L2=(-X,-d/2), R1=(+X,+d/2), R2=(+X,-d/2)
S ("same"):     n = (+1, +1, -1, -1)
O ("opposite"): n = (+1, -1, -1, +1)
```

Sektor fazowy (harmoniczny, scisle na torusie, z parzystej symetrii
G_per): pary wewnatrzgrupowe kasuja sie miedzy S i O; pary
miedzygrupowe [L1R1,L2R2: G(2X,0)] i [L1R2,L2R1: G(2X,d)] kasuja
sie w SUMIE:

```text
E_phase[S] + E_phase[O] = const(X), NIEZALEZNE od d   (dokladnie)
=>  E_sym(d) := (E_S + E_O)/2  mierzy WYLACZNIE kanal slepy na znak
    (dylatacyjny + sprzezenia krzyzowe amplituda-faza parzyste w n)
E_anty(d) := (E_O - E_S)/2 = 4*pi*J*[G(d) - G((2X,d))]  (znane
    dokladnie z sieciowej funkcji Greena) — wbudowana KONTROLA:
    wspolczynnik -> 4*pi*J = 8.662 (2x C2 z B1)
```

### 1.3 Rownanie #3 — dwie formy konkurencyjne dla kanalu E_sym

```text
(P) zrodlowana algebraiczna:  E_sym(d) = C - B * ln(d/xi) / d^2
    B_pred: policzone NUMERYCZNIE w Phase0_check z profilu ODE
    i ansatzu zlozeniowego (zero parametrow swobodnych)
(Y) swobodna Yukawa:          E_sym(d) = C - A * e^(-m*d),
    m in [0.8, 1.2] * m_a = [1.06, 1.59]  (skala masy dylatonu)
```

OCZEKIWANIE TEORETYCZNE (jawne): kanal istnieje, jest PRZYCIAGAJACY
(B > 0 lub A > 0) i wygrywa (P) — algebraiczny, bo zrodlowany faza.
Rozmiar: |E_sym(6) - E_sym(12)| rzedu 0.1-0.3 (oszacowanie
kappa*Phi**xi^2-owe; ilosciowo lockuje Phase0_check przed pomiarem).

### 1.4 Zastrzezenie interpretacyjne (pre-code)

Uniwersalne przyciaganie + krotki zasieg NIE jest jeszcze grawitacja
(brak zasiegu) — bedzie to wynik strukturalny: substrat MA kanal
o dobrym znaku; pytanie programu przesuwa sie na mechanizm zasiegu.
Odpychanie albo brak kanalu = falsyfikacja tej drogi.

## 2. Model i parametry (LOCKED)

```text
model, potencjal, siatka, dynamika: IDENTYCZNE z B1 (zero strojenia):
  kappa=0.50, a=0.50, b=1.60, c=1.00; psi in C;
  dpsi/dtau = kappa*Lap(psi) - psi*(a - b|psi| + c|psi|^2)
  N=256, L=128 (dx=0.5), periodycznie, dt=0.02, eps=0.30,
  S_GUARD=10 (flaga), seed=20260704
ansatz 4-wirowy: periodyzacja sinh (D1 z B1) + amplituda prof(r_i)
X = 32; skan d0 in {5, 6, 8, 10, 12} x {S, O}
tau_relax: dokladnie 200 krokow (D3 z B1), punkt kontrolny 150
pozycje wirow i d_meas: detektor windingu + centroid deficytu Phi
  (D2 z B1); pozycje wszystkich 4 wirow mierzone w chwili pomiaru E
kontrola dx=0.25 (OBOWIAZKOWA jak w B1): L=64, N=256, X=16,
  d0 in {5, 6, 8}
metrologia: D10 z B1 (H_Gamma co 10 krokow, tol. 1e-10; determinizm
  w Phase 1; brak flagi contact — pole wypelnia pudlo)
```

## 3. Procedura pomiarowa i fit (doprecyzowania pre-code)

- **P1 (fit laczny):** poniewaz podczas relaksacji d_meas[S] rosnie
  a d_meas[O] maleje (fazowe odpychanie/przyciaganie), E_sym NIE jest
  liczone punktowo, tylko przez LACZNY fit LSQ na wszystkich runach:
  `E_cfg = K + Cphi * gph_cfg + amp(d_cfg)`, gdzie
  `gph_cfg = -2*pi * sum_{i<j} n_i n_j G_lat(r_ij)` policzone z
  SIECIOWEJ funkcji Greena (FFT, Laplasjan 5-pkt — jak Phase0 B1)
  na ZMIERZONYCH pozycjach wirow; amp(d) wg formy (P) [k=3: K, Cphi,
  B] albo (Y) [k=4: K, Cphi, A, m skan siatkowy w pasmie];
  d_cfg = srednia d_meas obu grup. AIC = n*ln(RSS/n)+2k;
  Delta_AIC >= 10 rozstrzyga forme (inaczej BRAK — raport).
- **P2 (rezyduum bezmodelowe do G4):** E_res_cfg = E_cfg - K_hat -
  Cphi_hat * gph_cfg z fitu (P); E_res_sym(d0) = srednia rezyduow
  S i O dla danego d0 — empiryczny przebieg kanalu slepego na znak;
  G4 ocenia MONOTONIE tego przebiegu (Spearman = +1 na 5 punktach
  skanu) i rozstep > 10 * floor. FLOOR (doprecyzowane pre-code,
  przed pierwsza linijka kodu): caly fit laczny powtorzony na danych
  z kroku 150; floor = max_d0 |E_res_sym_200(d0) - E_res_sym_150(d0)|
  — stabilnosc wyekstrahowanego kanalu na transjent relaksacji
  (symetryczny dryft rdzeni absorbuje stala K i nie moze maskowac
  floora; roznice polozen absorbuje baza gph na zmierzonych
  pozycjach).
- **P3 (G3):** Cphi_hat / (4*pi*J = 8.6624)... UWAGA: Cphi mnozy
  gph na wszystkich parach, wiec odpowiada 2*pi*J w konwencji
  pairwise; kryterium: Cphi_hat / (2*pi*J = 4.3312) in [0.8, 1.2].
- **P4 (B_pred):** Phase0_check liczy E_sym_ansatz(d) na siatce
  (ansatz zlozeniowy, faza kasowana ta sama algebra) i fituje (P)
  -> B_pred zalockowane PRZED pomiarem; G5 pasmo B/B_pred in
  [0.5, 2.0] (analog pasma eikonalnego B2).
- **P5 (obserwacja, poza kryteriami):** dynamika pary rownoimiennej
  na malym d (geometria G5 z B1, d_ss=4): czy odpychanie slabnie
  ponizej czysto fazowego przewidywania (przyciagajaca poprawka
  dylatacyjna) — raport jakosciowy.

## 4. Kryteria (LOCKED)

- **G1 techniczne:** jak B1 (determinizm w Phase 1, brak NaN/runaway,
  H_Gamma nierosnace, kret zachowany n_tot=0 i |n_i|=1 wszedzie).
- **G2 ogon wiru:** profil 2D po relaksacji:
  |f(r)/[s*(1-xi^2/r^2)] - 1| <= 1% dla r in [6, 10].
  FALSYFIKATOR: ogon nie-algebraiczny -> rownanie #1 pada, forma (P)
  bez podstawy.
- **G3 dekompozycja dziala:** Cphi_hat/(2*pi*J) in [0.8, 1.2].
  FALSYFIKATOR: konstrukcja kwadrupolowa nie kasuje fazy jak w
  algebrze -> pomiar E_sym niewazny, stop przed interpretacja G4/G5.
- **G4 RDZEN — uniwersalne przyciaganie:** przebieg kanalu slepego
  na znak (P2) SCISLE ROSNACY na skanie (przyciaganie; Spearman +1)
  ORAZ rozstep > 10 * floor. FALSYFIKATOR: kanal pusty (rozstep
  ponizej floor) albo malejacy (odpychanie) -> uniwersalnosc znaku
  nie istnieje takze w sektorze skalarnym.
- **G5 forma i policzalnosc (klasyfikator + pasmo):** Delta_AIC >= 10
  rozstrzyga (P) vs (Y); przy (P): B/B_pred in [0.5, 2.0]; przy (Y):
  m w pasmie z #1.3. PASS gdy forma rozstrzygnieta i pasmo trzyma.
- **G6 kontrola dx=0.25:** znak i monotonia E_res_sym na kontroli
  zgodne z siatka glowna (3 punkty; raport ilosciowy).

## 5. Regula decyzyjna (LOCKED)

```text
G4 PASS + G5 (P) w pasmie
  -> substrat MA kanal uniwersalnie przyciagajacy (dylatacyjny),
     algebraiczny, policzalny. Program: mechanizm ZASIEGU
     (nastepne opy, za autoryzacja): (a) kanal dylatacyjny na tle
     niejednorodnym/wielkoskalowym, (b) polaczenie kanalow
     (zasieg fazowy x uniwersalnosc skalarna), (c) powrot B2 po
     stabilizacji obiektu.
G4 PASS + G5 (Y) w pasmie
  -> kanal uniwersalny, masywny (swobodny dylaton): jw., mechanizm
     zasiegu pilniejszy.
G4 PASS + G5 BRAK/poza pasmem
  -> kanal istnieje i przyciaga, ale forma/skala nieczytelna:
     raport; powtorka na gestszym skanie d jako osobny op.
G4 FAIL (odpychanie/pusty)
  -> takze sektor skalarny nie daje uniwersalnosci -> obie
     "proste" drogi znaku wyczerpane; priorytet: stabilizacja
     obiektu i B2-stacjonarne.
G2/G3 FAIL -> stop przed interpretacja (fundament pomiaru).
```

## 6. Forbidden moves

1. Zero pinningu; kret czytany z pola; neutralnosc calkowita (n_tot=0).
2. Zadnych zmian progow/parametrow/procedur po pierwszym runie
   (wyjatki techniczne jak B1: jedno obnizenie dt; kontrola dx=0.25
   jest czescia planu, nie wyjatkiem).
3. Tylko formy (P) i (Y); zadnej trzeciej formy po obejrzeniu danych.
4. Nie interpretowac krotkozasiegowego przyciagania jako grawitacji;
   tau != czas; zadnych deklaracji GR/MOND.
5. Wynik "odpychanie/pusty kanal" raportowac jako pelnoprawna
   falsyfikacje drogi C.
6. Nie promowac do core/ bez osobnej autoryzacji.

## 7. Czego ten cykl NIE robi

- Nie liczy zasiegu grawitacji ani nie testuje kanalu na tlach
  niejednorodnych (nastepne opy).
- Nie rozstrzyga B2 (czeka na tlo stacjonarne).
- Nie wiaze skali z jednostkami fizycznymi.

## 8. Deliverables

- `Phase0_balance.md` (ten plik, LOCKED przy utworzeniu — autoryzacja
  autora 2026-07-05: wybor drogi "sektor skalarny na tle Phi*")
- `Phase0_check.py`/`Phase0_output.txt` (algebra, ogon ODE,
  E_sym_ansatz(d), B_pred, G_lat)
- `Phase1_vortex_tail.py`/`Phase1_output.txt` (G2 + czesc G1)
- `Phase2_channel_decomposition.py`/`Phase2_output.txt` (G3-G6 + P5)
- `Phase_FINAL_close.md`, `README.md`
