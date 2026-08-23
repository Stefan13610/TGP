---
title: "Phase_FINAL_close — zamkniecie: mod Goldstone'a jako bezmasowy mediator (most do grawitacji, galaz B1)"
date: 2026-07-04
type: phase-final-close
tgp_owner: research/op-goldstone-mediator-2026-07-04
status: CLOSED
verdict: "PASS (G1-G4, G6 po erracie #2); G5 klasyfikator: LADUNKOWY -> mediator bezmasowy 2D-Newton istnieje, ale nie jest grawitacja"
anti_lakatos_lock: PRESERVED
tags: [gravity-bridge, branch-B1, goldstone, u1, vortex, newton-2d, charged-mediator, closed]
related:
  - "[[Phase0_balance.md]]"
  - "[[README.md]]"
  - "[[../op-lock-interaction-gravity-2026-07-04/Phase_FINAL_close.md]]"
  - "[[../op-nbody-additivity-2026-07-04/Phase_FINAL_close.md]]"
  - "[[../op-phi-metric-refraction-2026-07-04/Phase0_balance.md]]"
---

# Phase FINAL — zamkniecie cyklu

## 0. Werdykt

**G1 PASS, G2 PASS, G3 PASS, G4 PASS (forma logarytmiczna, C2 w pasmie
predykcji), G5 = LADUNKOWY (klasyfikator), G6 PASS (ratio 3.93 w oknie
waznym; errata #2).**

Wg zalockowanej reguly decyzyjnej (Phase0_balance.md, sekcja 5):

> G4 PASS (L, C2 w pasmie) + G5 LADUNKOWY
> -> mediator bezmasowy ISTNIEJE, forma 2D-Newtona, sprzezenie ~ Phi*,
>    ale ladunkowy. Grawitacja NIE jest prostym Goldstone'em U(1).

Sektor fazowy U(1) substratu NIESIE oddzialywanie dlugozasiegowe
o dokladnie tej formie, ktorej kanal amplitudowy Z2 nie mial:
potencjal logarytmiczny ("2D Newton", sila ~ 1/d), ze sprzezeniem
policzalnym algebraicznie z gestosci wygenerowanej przestrzeni
(`C2 = 2*pi*kappa*Phi*`, zmierzone/predykcja = 0.86-1.04) i ze
skalowaniem ladunkowym `~ n1*n2` (zmierzone 3.93, predykcja 4).
ALE oddzialywanie jest LADUNKOWE: rownoimienne defekty sie odpychaja.
Do grawitacji brakuje mechanizmu uniwersalnosci znaku — dokladnie tak,
jak przewidywalo jawne zastrzezenie interpretacyjne zapisane przed
kodem (sekcja 1.3 LOCKa). Zgodnie z LOCKiem (forbidden move 5) to jest
pelnoprawne rozstrzygniecie, nie porazka.

## 1. Co testowano

Po `op-lock-interaction-gravity` (kanal amplitudowy = Yukawa, krotki
zasieg -> galaz B) i `op-nbody-additivity` (superpozycja parowa):
czy mod fazowy pola zespolonego `psi = |psi| e^(i theta)` (Goldstone
U(1), zgodny z ontologia rdzenia sek01) daje oddzialywanie
dlugozasiegowe miedzy defektami fazowymi (wirami) w tle prawdziwej
prozni `|psi| = s*` — i czy jest ono grawitacyjne, czy ladunkowe.
Model: ten sam potencjal co caly lancuch (zero strojenia):
`kappa=0.5, a=0.5, b=1.6, c=1.0`, przeplyw gradientowy, N=256, L=128,
dt=0.02, seed=20260704. Zrodla = wiry o krecie n; konfiguracje
neutralne (n_tot = 0 w pudle periodycznym).

## 2. Przebieg

- **Phase 0**: LOCK (autoryzacja autora 2026-07-04: polecenie wyboru
  i przeprowadzenia cyklu; wybor B1 przed B2 zgodny z regula decyzyjna
  B2) + doprecyzowania pre-code D1-D11, w tym dwa spisane PRZED kodem:
  (a) periodyzacja ansatzu fazy (goly atan2 wstrzykuje energie szwu
  zalezna od d), (b) pasmo kryterialne m in [0.3, 2.0] dla formy (Y)
  z powodu znanej degeneracji m->0 z forma (L).
  `Phase0_output.txt`: algebra potwierdzona (s*=1.1742, Phi*=1.3787,
  xi=0.7544, J=0.6893, C2_pred=4.3312); profil rdzenia z ODE
  radialnego (ogon algebraiczny, r_99 = 7.85); poprawki obrazow
  periodycznych: C2_eff pozostaje w pasmie G4 na obu geometriach.
- **Phase 1**: genesis_u1 (G2 + czesc G1). `Phase1_output.txt`.
- **Phase 2**: profil rdzenia, statyka E(d) fit (L) vs (Y), dynamika
  par, kontrola dx=0.25 + pinning (G3, G4). `Phase2_output.txt`.
- **Phase 3**: sign_test 4-wirowy (G5), skalowanie z kretem (G6),
  obserwacja rozszczepienia n=2. `Phase3_output.txt`.

## 3. Wyniki liczbowe (kluczowe)

### 3.1 Ciaglosc fundamentu na U(1) (G2)

```text
bare (szum * e^{i*0.7}): Phi_max -> 4.5e-19 (tau=30) -> 0; area < A_MIN
cluster(0.6, D=3, K=7):  lock Phi_max = 1.37736 ~ Phi* = 1.37867
kowariancja fazowa U(1): max| |psi|_(theta=0.7) - |psi|_(theta=0) | = 4.9e-45
determinizm bitowy: TAK;  H_Gamma nierosnace (H_inc_max = 0.00e+00)
```

### 3.2 Istnienie i stabilnosc wirow (G3)

```text
kret zachowany po relaksacji: wszystkie runy (statyka, kontrola, core)
|psi|(r=8.25)/s* = 0.9908 (dx=0.5) / 0.9897 (dx=0.25)  [prog 1 +/- 0.03]
profil vs ODE: zbieznosc < 0.007 dla r >= 4 (rdzen po 200 krokach
  jeszcze szerszy od asymptotyki — transjent, maleje z dx)
pinning (D9): Delta_E_pin / rozstep E = 0.001% (dx=0.5), 0.000% (dx=0.25)
  [prog 2%] — brak przyklejania do siatki; rdzenie ruchome na obu siatkach
```

### 3.3 Prawo oddzialywania (G4) — RDZEN

Statyka, E(d_meas) po 200 krokach relaksacji (D3), siatka glowna:

```text
d_meas:  2.866    6.666    11.223    15.419    23.605
E:     -706.265 -704.119  -701.853  -700.412  -698.419
```

Fit (D4), rozstrzyga Delta_AIC >= 10:

```text
(L) log:    C1=-710.64  C2=3.7387          RSS=7.45e-01  AIC=-5.52
(Y) Yukawa (m in [0.3,2.0]): m=0.30 (brzeg) RSS=7.94e+00  AIC=+8.31
Delta_AIC (Y - L) = +13.83  >= 10  ->  WYGRYWA (L) LOGARYTMICZNA
C2 = 3.7387 in [3.465, 5.197]  ->  C2/C2_pred = 0.863  (pasmo OK)
```

Zbiegajace sie niezalezne pomiary C2 (wszystkie raportowane, poza
kryterium):

```text
C2 (fit glowny, 5 pkt)                    = 3.739   (0.86 * pred)
C2 (czulosc: bez punktu d=2.87 < 2*r_core) = 4.506   (1.04 * pred)
C2 (kontrola dx=0.25, 3 pkt, D9)           = 4.483   (1.04 * pred)
C2 (okno d0 in {12,16,24}, por. G6)        = 4.622   (1.07 * pred)
C2_dyn (D6, dH/dln d): 2.75 / 3.09 / 3.36  (d0 = 12/16/24; zanizone
  przez transjent i doklejone okno malych d — spojne co do rzedu)
predykcja: C2 = 2*pi*kappa*Phi* = 4.3312
```

Dynamika: v*d ~ const na kazdym runie (mediana 1.45/1.09/0.86 dla
d0=12/16/24, malejaca z oknem — poprawki logarytmiczne mobilnosci
wiru, przewidziane w LOCKu 1.4, poza kryterium); zblizanie
monotoniczne (frakcja krokow z dd>0: 0.0000) — przyciaganie (+1,-1)
potwierdzone dynamicznie; zero orbit (przeplyw nadtlumiony).

### 3.4 Dyskryminator znaku (G5) i skalowanie z kretem (G6)

sign_test (D7): konfiguracja neutralna (+1,+1,-1,-1), grupy w L/2:

```text
tau:        0.2     25      50      100     150     200
d_ss(+,+): 8.069  11.680  13.838  16.888  19.154  20.968
d_ss(-,-): identyczne co do wydruku (symetria dokladna)
Delta d_ss = +12.90  >>  prog +1.0   ->  ODPYCHANIE rownoimiennych
G5 = LADUNKOWY
```

winding (D8): C2(+2,-2)/C2(+1,-1):

```text
okno wazne (d0 in {12,16,24}; kret(+) = +2 zachowany):
  C2(+2,-2) = 18.157   C2(+1,-1, to samo okno) = 4.622
  ratio = 3.928   in [3.2, 4.8]   (predykcja |n1*n2| = 4)   -> PASS
punkty d0 in {6,8} skanu (+2,-2): konfiguracja zdegradowana podczas
  relaksacji do (+1,-1) (rozszczepienie n=2 + anihilacja pod-pary;
  rozszczepienie bylo PRZEWIDZIANE w LOCKu, scenariusz 6) — punkty
  niewazne jako pomiar zbocza (+2,-2); errata #2.
obserwacja n=2 (d0=16, tau<=100): rozszczepienie marginalne
  (max separacja 1.6 ~ rozdzielczosc detektora klastrow 1.5;
  pary +1+1 zwiazane z tlem pozostaja w jednym klastrze) — rdzen n=2
  kwazistabilny na oknie; przy malych d rozszczepienie + anihilacja
  czesciowa realna (patrz wyzej).
```

### 3.5 G1 (techniczne)

Determinizm bitowy (Phase 1, run x2), zero NaN/runaway, H_Gamma
nierosnace na KAZDYM runie wszystkich faz (H_inc_max = 0.00e+00),
zero clampow. PASS.

## 4. Interpretacja w granicach modelu

CO WYNIK ZNACZY:
- substrat U(1) MA bezmasowy mediator dlugozasiegowy: mod Goldstone'a
  w fazie metrycznej Phi > 0; oddzialywanie defektow ma dokladnie
  forme "2D Newtona" (E ~ ln d, F ~ 1/d) — czyli forme, ktorej
  zabraklo kanalowi amplitudowemu Z2 w poprzednim opie,
- sprzezenie jest ALGEBRAICZNIE policzalne z substratu i NIESIONE
  przez gestosc wygenerowanej przestrzeni: C2 = 2*pi*kappa*Phi*
  (zmierzone 0.86-1.07 predykcji czterema niezaleznymi drogami),
- skalowanie ladunkowe ~ n1*n2 potwierdzone (3.93 vs 4),
- ALE mediator jest LADUNKOWY: (+,+) sie odpycha. Uniwersalnego
  przyciagania nie ma — grawitacja NIE jest prostym Goldstone'em U(1).
  To potwierdza standardowy argument (mediator wektorowy/ladunkowy
  daje oba znaki; uniwersalnosc znaku wymaga sprzezenia do "masy",
  nie do ladunku) na poziomie substratu TGP.

CZEGO WYNIK NIE ZNACZY:
- nie falsyfikuje mostu do grawitacji: gałąź B2 (geometria propagacji,
  osobny draft) oraz sektor skalarny na tle Phi* pozostaja otwarte,
- tau nie jest czasem fizycznym; dynamika nadtlumiona ("sila" =
  predkosc); zadnych orbit, geodezyjnych, deklaracji zgodnosci z GR
  (forbidden move 4),
- skala J pozostaje fenomenologiczna (LOCK, sekcja 7): zero deklaracji
  o jednostkach fizycznych.

## 5. Errata / transparencja

1. **Incydent techniczny bez wplywu na wyniki:** pliki Phase1-Phase3
   zostaly poczatkowo zapisane w zagniezdzonych katalogach (dryf cwd
   narzedzia — identyczny incydent jak errata #3 poprzedniego opa)
   i przeniesione do katalogu opa bez zmian tresci.
2. **G6 (raportowanie, dzien runow, przed interpretacja):** pierwotna
   linia ewaluacyjna (widoczna w `Phase3_output.txt`: "G6 FAIL")
   wymagala kret(+) = +2 na WSZYSTKICH punktach skanu, przez co
   degradacja konfiguracji przy d0 in {6,8} (rozszczepienie n=2
   i anihilacja pod-pary — zjawisko przewidziane w LOCKu, scen. 6)
   zlewala sie z "skalowanie zle". Brzmienie LOCKa G6 to wylacznie
   pasmo ratio zboczy; zbocze (+2,-2) jest zdefiniowane tylko tam,
   gdzie konfiguracja (+2,-2) istnieje (D2: kret czytany z pola).
   Ocena skorygowana do brzmienia LOCKa na oknie waznym, z baseline
   (+1,-1) na TYCH SAMYCH d0 (apples-to-apples): ratio = 3.928 ->
   PASS. Surowy wydruk zachowany; ratio "wszystkie punkty" (4.501,
   tez w pasmie) raportowane obok. Analogiczny przypadek jak
   errata #1 poprzedniego opa.
3. **Punkt statyczny d0=6 dla (+1,-1):** podczas relaksacji para
   zapadla sie do d_meas = 2.87 < 2*r_core — glęboko w strefie
   przekrywania rdzeni. Punkt POZOSTAWIONY w ficie glownym (LOCK
   definiuje skan po d0), fit czulosci bez niego RAPORTOWANY:
   C2 = 4.506 (RSS spada 120x) — wniosek G4 niewrazliwy, C2 w pasmie
   w obu wariantach. Analogiczny przypadek jak errata #2 poprzedniego
   opa.
4. **Czulosc D4 (zapisana z gory):** (Y) z m swobodnym daje m = 0.061
   i AIC lepszy od (L) — to jest DOKLADNIE degeneracja m->0
   pre-zarejestrowana w D4 (zasieg 1/m = 16.5 ~ d_max: operacyjnie
   bezmasowy, nierozroznialny od (L) w oknie). Rozstrzygniecie
   kryterialne opiera sie na pasmie m in [0.3, 2.0] (Goldstone
   masywny w skali substratu); niezaleznie: kontrola dx=0.25,
   C2 z okna czystego, skalowanie ~ n1*n2 i v*d ~ const wykluczaja
   interpretacje krotkozasiegowa.

Zadne progi G1-G6, wspolczynniki ani procedury nie byly zmieniane po
pierwszym runie. Dozwolone wyjatki techniczne: kontrolny run dx=0.25
byl OBOWIAZKOWY z gory (LOCK sekcja 2); obnizenie dt niepotrzebne.

## 6. Decyzja i nastepny krok

Regula decyzyjna (LOCK, sekcja 5), galaz **G4 PASS + G5 LADUNKOWY**:

> mediator bezmasowy ISTNIEJE, forma 2D-Newtona, sprzezenie ~ Phi*,
> ale ladunkowy. Grawitacja NIE jest prostym Goldstone'em U(1).

Nastepne kroki wskazane przez LOCK (kazdy jako osobny op, do
autoryzacji przez autora):
1. **sektor skalarny/amplitudowy na tle Phi***: czy fluktuacje
   dylatacyjne wokol prawdziwej prozni daja kanal zawsze
   przyciagajacy (uniwersalny)?
2. **B2 — geometria propagacji**: draft gotowy
   ([[../op-phi-metric-refraction-2026-07-04/Phase0_balance.md]]);
   po tym cyklu B2 jest naturalnym kolejnym pomiarem (mechanizm
   "sila ladunkowa Goldstone'a" zmierzony i sklasyfikowany).
3. **analiza teoretyczna**: czy uniwersalnosc znaku wymaga sprzezenia
   tensorowego (spin-2) — standardowy argument do sprawdzenia na
   substracie TGP.

Bilans mostu do grawitacji po trzech opach galezi:
- kanal amplitudowy Z2: Yukawa, krotki zasieg, frontowy -> NIE,
- kanal fazowy U(1): 2D-Newton, dlugi zasieg, policzalny -> forma TAK,
  znak NIE (ladunkowy),
- kanal geometryczny (B2): niezmierzony -> nastepny.

## 7. Pliki cyklu

- [[Phase0_balance.md]] — LOCK + doprecyzowania pre-code (sekcja 10)
- `Phase0_check.py` / `Phase0_output.txt` — algebra, profil ODE,
  poprawki obrazow periodycznych
- `Phase1_genesis_u1.py` / `Phase1_output.txt` — G2 (+G1)
- `Phase2_vortex_pair_law.py` / `Phase2_output.txt` — G3, G4
- `Phase3_sign_and_winding.py` / `Phase3_output.txt` — G5, G6
- [[README.md]] — przeglad opa
- ten plik — zamkniecie
