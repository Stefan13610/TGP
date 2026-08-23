---
title: "Phase_FINAL_close — zamkniecie: Phi jako metryka efektywna propagacji (most do grawitacji, galaz B2)"
date: 2026-07-04
type: phase-final-close
tgp_owner: research/op-phi-metric-refraction-2026-07-04
status: CLOSED
verdict: "G1-G2 PASS; G3/G4/G6 NIEAPLIKOWALNE (pre-rejestrowana klauzula D6: tlo zamrozone tachionowe); G5 empirycznie nierozstrzygniete — B2 na tle zamrozonym niemierzalne, wymaga tla stacjonarnego"
anti_lakatos_lock: PRESERVED
tags: [gravity-bridge, branch-B2, effective-metric, refraction, eikonal, wall-tachyon, closed]
related:
  - "[[Phase0_balance.md]]"
  - "[[README.md]]"
  - "[[../op-goldstone-mediator-2026-07-04/Phase_FINAL_close.md]]"
  - "[[../op-lock-interaction-gravity-2026-07-04/Phase_FINAL_close.md]]"
  - "[[../op-wall-dynamics-2026-07-03/README.md]]"
---

# Phase FINAL — zamkniecie cyklu

## 0. Werdykt

**G1 PASS, G2 PASS. G3/G4/G6 NIEAPLIKOWALNE wg pre-rejestrowanej
klauzuli D6(d). G5 empirycznie NIEROZSTRZYGNIETE (strona teoretyczna:
raport eikonalny D7, zalockowany przed pomiarem).**

Glowny WYNIK cyklu (odkryty rachunkiem przed kodem, pre-zarejestrowany
w D6, potwierdzony rachunkowo w Phase 0 i empirycznie w Phase 2):

> **Zamrozone tlo zlockowanego obiektu jest TACHIONOWE w scianie.**
> Potencjal lancucha ma `V''(s) < 0` dla `s in (0.190, 0.874)`
> (min -0.353 przy s = 0.533); sciana obiektu przechodzi przez to
> pasmo, a studnia przyciagajaca zawsze wiaze stan:
> `lambda_min = -0.0152`, `gamma = 0.123` (Phase 0, iteracja potegowa);
> empirycznie `gamma_meas = 0.125` (zgodnosc 1.3%), mod w 99.99%
> zlokalizowany na scianie. Fizycznie to zamrozony mod WZROSTU frontu:
> obiekt w tym substracie rosnie, wiec jego "zatrzymana" sciana nie
> jest stanem rownowagi. Amplifikacja e^(gamma*tau_lotu) ~ 1e8
> kontaminuje kazdy pomiar ugiecia przed oknem pomiarowym (0/12
> waznych runow; S_GUARD osiagany przed pierwsza probka alpha).

Minimalna realizacja B2 ("sam sektor falowy na zamrozonym tle,
zero nowych sprzezen") jest wiec NIEMIERZALNA — nie dlatego, ze
refrakcja nie istnieje, lecz dlatego, ze samo zamrozenie tla generuje
niestabilnosc. To pelnoprawne rozstrzygniecie techniczno-fizyczne
z jasnym wskazaniem brakujacej struktury: **B2 wymaga tla
STACJONARNEGO** (stabilizacja obiektu — watek op-wall-dynamics)
albo sektora propagacji bez sciany w pasmie V'' < 0.

## 1. Co testowano

Galaz B2: "grawitacja" = refrakcja malych zaburzen w gradiencie Phi.
Minimalne rozszerzenie (jedyna nowa struktura, zadeklarowana wprost):
sektor falowy `d^2u/dtau^2 = kappa*Lap(u) - V''(s_bg)*u` na zamrozonym
tle zlockowanego obiektu (cluster 0.6, D=3, K=7 -> R_obj = 12).
Dwa pytania: (P1) policzalnosc kata z eikonalu, (P2) znak (ku/od
obiektu). Model bez strojenia: `kappa=0.5, a=0.5, b=1.6, c=1.0`,
N=256, L=128, dt=0.05 (leapfrog), seed=20260704.

## 2. Przebieg

- **Phase 0**: LOCK (autoryzacja: polecenie kontynuacji po zamknieciu
  B1) + doprecyzowania D1-D11, w tym PRE-REJESTRACJA ryzyka tachionu
  sciany (D6: rachunek V''<0 + procedura + kryterium waznosci +
  werdykt warunkowy — wszystko przed kodem). `Phase0_output.txt`:
  masy (m_FV=0.7071, m_TV=0.9374), n_in(1.1)=0.6831 (przyklad LOCKa:
  0.68 OK), tlo R_obj=12.000, **lambda_min = -0.01523 (tachion)**,
  ray tracing alpha_pred(b, omega) zalockowany przed pomiarem.
- **Phase 1**: kalibracja dyspersji (G2 + G1). `Phase1_output.txt`.
- **Phase 2**: tachyon_check + skan deflection 12 runow + bariera +
  liniowosc + tlo zywe (G3-G6 wg D6). `Phase2_output.txt`.

## 3. Wyniki liczbowe (kluczowe)

### 3.1 Kalibracja sektora falowego (G2, G1)

```text
dyspersja omega(k), k in [0.49, 1.47], tla FV i TV (10 pomiarow):
  max |blad| vs omega^2 = kappa*k^2 + V'' :  1.52%   (prog 5%)  PASS
  bledy ida za dyspersja SIECIOWA (kappa*(2-2cos(k dx))/dx^2) — jakosc
  dyskretyzacji, nie fizyki; H_u: dryft sekularny max 1.8e-05 (prog
  1e-4), oscylacja leapfrog O((dt*omega)^2) ~ 5e-4 raportowana (D4);
  determinizm bitowy: TAK
```

### 3.2 Tachion sciany (D6) — RDZEN faktyczny cyklu

```text
Phase 0 (rachunek):   lambda_min(L = -kappa*Lap + V''(s_bg)) = -0.01523
                      gamma = sqrt(-lambda_min) = 0.1234
                      mod: <r> = 12.24 ~ R_obj; 99.99% normy w scianie
Phase 2 (empiria):    gamma_meas = 0.1250  (zgodnosc 1.3%)
                      wzrost z szumu 1e-8 do S_GUARD=10 w tau = 195
prognoza kontaminacji: e^(gamma * tau_lotu ~ 155) ~ 2e8  — potwierdzona:
  WSZYSTKIE 12 runow deflection: runaway PRZED pierwsza probka alpha
  (0/12 waznych wg D6c); bariera (omega=0.8): to samo
tlo zywe (D11, obserwacja): wzrost przenosi sie na RUCHOMA sciane
  (gamma_wall w masce statycznej: -0.017, ale runaway rowniez) —
  mieszana dynamika (tlo I rzedu + fala II rzedu) nie usuwa pasma
  V''<0; obiekt dodatkowo rosnie w oknie pomiaru (R: 12 -> ~36)
```

### 3.3 Strona teoretyczna P1 (eikonal D7, zalockowana przed pomiarem)

```text
alpha_pred [deg] (>0 = OD obiektu):     b=8      b=12     b=16     b=20
  omega=1.0:                          +71.1    +24.8    -92.4    -0.24
  omega=1.1:                          +61.4    +19.0    -65.6    -0.17
  omega=1.3:                          +45.0    +11.1    -13.7    -0.10
```

Struktura znaku BOGATSZA niz naiwne oczekiwanie LOCKa (sekcja 1.2
"n < 1 wewnatrz -> zawsze OD"): pasmo `V''(s) < 0` w scianie tworzy
PIERSCIEN o `n^2 > 1` (skupiajacy). Promienie przecinajace wnetrze
(b <= 12) uginaja sie OD obiektu; promienie muskajace pierscien
(b ~ 16 ~ R_obj + polowa grubosci sciany) sa silnie przyciagane KU
obiektowi (az do zakrzywienia -92 deg — rezim soczewki krawedziowej);
b = 20 praktycznie nieugiete. Uwaga o tej mozliwosci zostala zapisana
w D7 PRZED uruchomieniem ray tracingu. To obserwacja teoretyczna
(nie pomiar falowy!) — ale pokazuje, ze "znak refrakcji" w tym
substracie nie jest jedna liczba, tylko funkcja parametru zderzenia.

## 4. Interpretacja w granicach modelu

CO WYNIK ZNACZY:
- minimalna realizacja B2 (fala na zamrozonym tle) jest WEWNETRZNIE
  sprzeczna: zamrozenie nie-stacjonarnego obiektu robi z modu wzrostu
  frontu tachion — kazdy dluzszy pomiar na takim tle jest zdominowany
  przez artefakt zamrozenia, nie przez fizyke propagacji,
- to jest POZYTYWNA wiedza pomiarowa: (i) ilosciowa (lambda_min,
  gamma, lokalizacja — zgodnosc rachunku z empiria 1.3%),
  (ii) strukturalna: geometria propagacji w TGP wymaga najpierw
  STACJONARNEGO obiektu (stabilizacja — op-wall-dynamics) albo
  pytania zadanego tak, by sciana nie lezala w pasmie V'' < 0,
- eikonal (strona teoretyczna P1) dziala i daje falsyfikowalne,
  nietrywialne predykcje (znak zalezny od b, pierscien skupiajacy) —
  gotowe do konfrontacji z pomiarem falowym, gdy tlo bedzie
  stacjonarne.

CZEGO WYNIK NIE ZNACZY:
- NIE rozstrzyga P2 (grawitacyjnosc kierunku) — pomiar falowy sie nie
  odbyl; deklarowanie "refrakcja odpychajaca" na podstawie samego
  eikonalu byloby zlamaniem forbidden move 4,
- NIE falsyfikuje B2 jako mechanizmu — falsyfikuje jedna REALIZACJE
  (tlo zamrozone z rosnacym obiektem),
- zadnych wnioskow o metryce g_mu_nu / GR / soczewkowaniu
  obserwacyjnym (forbidden moves 4-5).

## 5. Errata / transparencja

1. **Phase 1, niezgodnosc kodu z D1 (naprawiona przed interpretacja):**
   pierwsza wersja liczyla u_tau roznica WSTECZNA zamiast centralnej
   zadeklarowanej w D1, co zawyzalo oscylacje H_u (O(dt*omega))
   i przez aliasing zawyzalo miare dryftu; dodatkowo dryft (kryterium
   G1) byl blednie wliczony do oceny G2. Po naprawie do brzmienia D1:
   dryft max 1.8e-05, G1/G2 PASS. Pierwsza wersja outputu zastapiona;
   zmiana czysto implementacyjna (zgodnosc z LOCKiem), zero zmian
   progow.
2. **Kontaminacja w runach deflection** objawila sie flaga S_GUARD
   (runaway) PRZED pierwsza probka alpha, wiec ratio kontaminacji
   z D6(c) nie bylo liczone (NaN w tabeli) — waznosc rozstrzygala
   flaga runaway, co jest przypadkiem silniejszym (kontaminacja
   100%), w pelni w ramach D6(c).
3. Werdykt "NIEAPLIKOWALNE" dla G3/G4/G6 pochodzi w calosci
   z klauzuli D6(d) spisanej PRZED pierwsza linijka kodu — zadne
   progi ani procedury nie byly zmieniane po runach.

## 6. Decyzja i nastepny krok

Pre-rejestrowany werdykt D6(d) + regula decyzyjna LOCKa (galaz
"G4 FAIL / nieczytelnosc"): **B2 w tej postaci ODLOZONE** z precyzyjnym
wskazaniem brakujacej struktury:

1. **Tlo stacjonarne** — warunek wstepny powrotu do B2: stabilizacja
   zlockowanego obiektu (watek [[../op-wall-dynamics-2026-07-03/README.md]])
   albo konfiguracja bez sciany w pasmie V''<0. Po stabilizacji:
   powtorka TEGO SAMEGO protokolu (eikonal D7 juz zalockowany,
   predykcje w tabeli 3.3 zostaja jako zaklad).
2. **Sektor skalarny/amplitudowy na tle Phi*** (z reguly decyzyjnej
   B1) — niezalezna droga do uniwersalnego przyciagania, bez sciany.
3. Sprzezenie `kappa(Phi)/c(Phi)` (aksjomat rdzenia) — dopiero PO
   pomiarze na tle stacjonarnym (kolejnosc z reguly decyzyjnej LOCKa).

Bilans mostu do grawitacji po czterech opach galezi:
- kanal amplitudowy Z2: Yukawa, krotki zasieg -> NIE,
- kanal fazowy U(1): 2D-Newton, dlugozasiegowy, policzalny,
  ale LADUNKOWY -> forma TAK, znak NIE,
- kanal geometryczny B2: na tle zamrozonym NIEMIERZALNY (tachion
  sciany, zmierzony ilosciowo) -> wymaga tla stacjonarnego,
- najblizsze otwarte drogi: stabilizacja obiektu, sektor skalarny
  na tle Phi*.

## 7. Pliki cyklu

- [[Phase0_balance.md]] — LOCK + doprecyzowania pre-code (sekcja 10,
  w tym pre-rejestracja tachionu D6)
- `Phase0_check.py` / `Phase0_output.txt` — masy, n(omega), tlo,
  lambda_min, eikonal alpha_pred
- `Phase1_dispersion.py` / `Phase1_output.txt` — G2 (+G1)
- `Phase2_deflection.py` / `Phase2_output.txt` — tachyon_check,
  skan deflection, bariera, liniowosc, tlo zywe (G3-G6 wg D6)
- [[README.md]] — przeglad opa
- ten plik — zamkniecie
