---
title: "Phase0_balance — DRAFT: oddzialywanie zlockowanych obiektow jako kandydat na most do grawitacji"
date: 2026-07-04
type: phase0-draft
tgp_owner: research/op-lock-interaction-gravity-2026-07-04
status: LOCKED
anti_lakatos_lock: LOCKED (2026-07-04, autoryzacja autora — polecenie przeprowadzenia cyklu; lock przed pierwsza linijka kodu)
tags: [gravity-bridge, lock-interaction, yukawa, force-law, substrate]
related:
  - "[[../op-bare-substrate-genesis-2026-07-04/Phase_FINAL_close.md]]"
  - "[[../op-bare-substrate-genesis-2026-07-04/Phase0_balance.md]]"
  - "[[../op-wall-dynamics-2026-07-03/README.md]]"
  - "[[../../core/sek01_ontologia/sek01_ontologia.tex]]"
---

# Phase 0 — DRAFT: oddzialywanie zlockowanych obiektow (most do grawitacji, krok 1)

## 0. Hipoteza robocza

Hipoteza autora (2026-07-04, po PASS `op-bare-substrate-genesis`):

> Skoro substrat potrafi tworzyc zlockowane obiekty `Phi > 0`, grawitacja
> powinna byc oddzialywaniem pol miedzy nimi — moze nie prostym przy
> wielu obiektach, ale dla prostych przypadkow POLICZALNYM, z szansa na
> wyprowadzenie algebraicznych rownan oddzialywania.

Ten cykl testuje wersje minimalna: **dwa zlockowane obiekty w tym samym
substracie. Czy oddzialuja? Jakim prawem? Czy to prawo jest algebraicznie
wyprowadzalne z substratu i czy ma cechy grawitacji (dlugi zasieg,
skalowanie z "masa"), czy raczej cechy sily powierzchniowej (krotki
zasieg, skalowanie z geometria frontu)?**

## 1. Predykcje algebraiczne zapisane PRZED kodem (rdzen anty-Lakatos)

Zeby wynik byl rozstrzygajacy, DWIE konkurencyjne formy prawa sily sa
zalockowane z gory, wraz z teoretycznym oczekiwaniem:

### 1.1 Rownanie #1 — ogon pola w falszywej prozni

Wokol zlockowanego obiektu pole `s` zanika w `Phi = 0` z masa ekranowania:

```text
m0 = sqrt(V''(0) / kappa) = sqrt(a / kappa) = sqrt(0.5/0.5) = 1.0
s(r) ~ K0(m0 * r)  ~  e^(-m0 r) / sqrt(r)    (2D, zmodyfikowany Bessel)
```

### 1.2 Rownanie #2 — dwie zalockowane formy oddzialywania

Energia oddzialywania dwoch obiektow przedzielonych szczelina `g`
(mierzalna statycznie: `E_int(g) = H[s1+s2] - H[s1] - H[s2]`):

```text
(Y)  Yukawa / substratowa:   E_int(g) = -A * e^(-m*g)     (krotki zasieg)
(P)  Potegowa / gravity-like: E_int(g) = -B * g^(-p)      (dlugi zasieg)
```

TEORETYCZNE OCZEKIWANIE (zapisane wprost, PRZED kodem, zeby wyniku nie
dalo sie reinterpretowac po fakcie): pole `s` z symetria Z2 NIE MA modu
bezmasowego, wiec substrat powinien dac forme (Y) z `m ~ m0 = 1.0`.
Jesli tak wyjdzie — to jest PELNOPRAWNY wynik: bezposrednie oddzialywanie
lockow jest krotkozasiegowe i NIE jest kandydatem na grawitacje; most
wymaga innego mechanizmu (galaz B reguly decyzyjnej). Forma (P) bylaby
niespodzianka i mocnym kandydatem na most (galaz A).

### 1.3 Rownanie #3 — prawo frontu (kalibracja)

Front pojedynczej domeny w 2D (Allen-Cahn):

```text
v(R) = v_inf - c / R      (przewidywanie: c ~ kappa = 0.5)
```

To rownanie MUSI byc dopasowalne na kontroli — inaczej pomiar oddzialywania
(definiowany jako nadwyzka nad wzrost swobodny) nie ma sensu.

### 1.4 Uwaga o dynamice nadtlumionej

Przeplyw gradientowy jest nadtlumiony: "sila" objawia sie jako PREDKOSC
zblizania, nie przyspieszenie (`v ~ -dE_int/dg`, mobilnosc efektywna
fitowana). Ten cykl mierzy STATYKE prawa oddzialywania (forme `E_int(g)`
i skalowania), nie dynamike orbit. Inercja/geodezyjne = poza zakresem.

## 2. Model i parametry (do zalockowania)

Substrat i potencjal IDENTYCZNE jak w op-bare-substrate-genesis
(ciaglosc lancucha: zero strojenia miedzy cyklami):

```text
kappa = 0.50, a = 0.50, b = 1.60, c = 1.00, Phi = s^2
dynamika: ds/dtau = kappa*Lap(s) - V'(s), jawny Euler
grid: N = 256, L = 128.0 (dx = 0.5), brzegi periodyczne, Laplasjan 5-pkt
dt = 0.02
eps = 0.30, A_min = 4 wezly / N^2, S_GUARD = 10 (flaga, nie clamp)
przygotowanie obiektu: cluster(A0=0.6, D=3, K=7) ewoluowany do zadanego
  promienia R_obj (promien = polowa szerokosci Phi>eps wzdluz osi x)
para: s_pair = s_obj1(przesuniety o -X0) + s_obj2(przesuniety o +X0)
  (blad superpozycji ogonow ~ e^(-g) — raportowany, dopuszczalny dla g>=2)
szczelina g = odleglosc miedzy przecieciami Phi=eps na osi laczacej
skan szczelin: g0 in {2, 3, 4, 5, 6, 8}
skan rozmiarow (przy g0=4): R_obj in {8, 12, 16}
okno pomiarowe: od tau=0 do momentu g < 1.0 albo boundary_contact
seed: 20260704 (cykl deterministyczny; seed dla ewentualnych wariantow)
```

## 3. Scenariusze

1. `control_front(R)`: pojedyncza domena, pomiar `R(tau)` -> fit rownania
   #3 (`v_inf`, `c`). Kalibracja odniesienia dla "wzrostu swobodnego".
2. `pair_static(g0)`: statyczny pomiar `E_int(g0)` przez kompozycje pol
   (bez ewolucji) dla wszystkich g0 -> fit form (Y) vs (P).
3. `pair_dynamic(g0)`: ewolucja pary; `Delta_v(g) = v_zblizania(g) -
   2*v_control(R(tau))` (kontrola dopasowana krzywizna/promieniem).
   Niezalezne potwierdzenie prawa sily z dynamiki.
4. `mass_scaling(R1, R2)`: pary rozmiarow {8x8, 12x12, 16x16} przy g0=4.
   Dyskryminator: sila masowa vs frontowa.
5. `universality`: dwa obiekty o TYM SAMYM R_obj przygotowane ROZNA
   historia (A0=0.6 dluzej vs A0=1.0 krocej) -> czy oddzialywanie zalezy
   tylko od stanu, nie od historii.

## 4. Obserwable

```text
- R(tau), v(R) na kontroli; parametry v_inf, c fitu #3 + R^2
- E_int(g) statycznie; fit (Y): A, m; fit (P): B, p; R^2 i AIC obu form
- g(tau), Delta_v(g) dynamicznie; spojnosc znaku i ksztaltu z E_int
- ratio_mass = Delta_v(16x16) / Delta_v(8x8) przy g0=4
- ratio_hist = Delta_v(historia1) / Delta_v(historia2) przy tych samych R,g
- H_Gamma(tau) nierosnace; finite/runaway/boundary_contact flagi
```

## 5. Kryteria (propozycja do LOCKa)

### G1 — techniczne
Determinizm, brak NaN/runaway, H_Gamma nierosnace. PASS/FAIL.

### G2 — kalibracja frontu (rownanie #3)
Fit `v(R) = v_inf - c/R` na kontroli z `R^2 >= 0.98` i `c in [0.2, 1.0]`.
FALSYFIKATOR: front nie spelnia prawa krzywizny -> odniesienie
"wzrost swobodny" niezdefiniowane, pomiar oddzialywania niewykonalny.

### G3 — istnienie oddzialywania
`E_int(g) < 0` (przyciaganie) i `|E_int|` malejace z g, ORAZ
`Delta_v(g) > 0` dla g <= 4 w dynamice (zblizanie ponad wzrost swobodny).
FALSYFIKATOR: brak mierzalnego oddzialywania -> locki sa "gluche" na
siebie -> hipoteza mostu przez bezposrednie oddzialywanie falsyfikowana
w najmocniejszy sposob.

### G4 — RDZEN: prawo jest algebraiczne i rozstrzygalne
Fit (Y) vs (P) na `E_int(g)`: forma wygrywa, gdy `Delta_AIC >= 10`.
Jesli wygrywa (Y): `m_fit in [0.8, 1.2]` potwierdza wyprowadzenie z
substratu (predykcja m0 = 1.0 z rownania #1) — prawo oddzialywania
JEST algebraicznie policzalne z modelu.
FALSYFIKATOR: zadna forma nie wygrywa / m_fit poza pasmem -> prawo
nieczytelne w tym modelu, "policzalnosc prostych przypadkow" pada.

### G5 — dyskryminator: sila masowa czy frontowa
Przy stalym g0=4:
```text
ratio_mass >= 2.0  -> skalowanie MASOWE (bulk, gravity-like)
ratio_mass <= 1.5  -> skalowanie FRONTOWE (powierzchniowe)
posrednie          -> MIXED, raport wprost
```
To kryterium NIE ma zdefiniowanego "PASS" — jest klasyfikatorem galezi
reguly decyzyjnej (obie odpowiedzi sa informacja).

### G6 — uniwersalnosc stanu
`|ratio_hist - 1| <= 0.10`: oddzialywanie zalezy od stanu obiektow, nie
od ich historii.
FALSYFIKATOR: historia przygotowania zmienia sile -> "obiekt" nie jest
dobrze zdefiniowany przez stan pola -> pojecie oddzialywania miedzy
lockami zle postawione.

## 6. Regula decyzyjna cyklu

```text
G3 FAIL (brak oddzialywania)
  -> bezposrednie oddzialywanie lockow NIE istnieje -> most przez sile
     falsyfikowany; przejdz od razu do drogi geometrycznej (galaz B2).

G3+G4 PASS, forma (P) potegowa, G5 masowe
  -> GALAZ A: mocny kandydat na most. Nastepny cykl: granica
     "newtonowska" (UWAGA: w 2D Newton = potencjal logarytmiczny,
     sila ~ 1/r — porownywac z ta forma!), uklady 3-cialowe,
     stala sprzezenia G_eff z parametrow substratu.

G3+G4 PASS, forma (Y) Yukawa (oczekiwana), G5 frontowe
  -> GALAZ B: bezposrednia sila miedzy lockami jest krotkozasiegowa
     i powierzchniowa -> NIE jest grawitacja. Fundament stoi, mechanizm
     do wymiany. Dwie drogi (osobne opy):
     B1: bezmasowy mediator — pole zespolone s in C, U(1), mod fazowy
         (Goldstone) daje oddzialywanie dlugozasiegowe; lock z
         op-bare-substrate-genesis powtorzyc dla |s|.
     B2: droga geometryczna — Phi jako metryka efektywna: predkosc
         propagacji malych zaburzen zalezna od Phi(x); "grawitacja" =
         refrakcja trajektorii w gradiencie Phi, nie sila miedzy cialami.

Wynik mieszany raportowac wprost, bez strojenia progow.
```

Uwaga strukturalna: galaz B NIE jest porazka programu — jest pomiarem,
ktory z trzech mechanizmow (sila bezposrednia / mediator bezmasowy /
geometria propagacji) niesie grawitacje w tej ontologii. Kazda galaz ma
zdefiniowany nastepny krok.

## 7. Forbidden moves (propozycja)

1. Nie przypinac obiektow (zero external potentials, zero trzymania
   pozycji). Wszystko czytane z jednej dynamiki pola.
2. Nie zmieniac progow G1-G6, wspolczynnikow a,b,c,kappa ani procedur po
   pierwszym runie (wyjatki techniczne jak poprzednio: jedno obnizenie dt;
   jeden kontrolny run na drobniejszej siatce dla przypadkow granicznych).
3. Nie dofitowywac TRZECIEJ formy prawa sily po obejrzeniu danych.
   Rozstrzygniecie tylko (Y) vs (P) — inne formy = osobny, nowy cykl.
4. Nie interpretowac nadtlumionej dynamiki jako orbit/geodezyjnych.
5. Nie mylic tau z czasem fizycznym; nie deklarowac zgodnosci z GR/MOND.
6. Wynik "Yukawa, frontowe" raportowac jako pelnoprawne rozstrzygniecie
   (galaz B), nie jako porazke do ukrycia.
7. Nie promowac do core/ bez osobnej autoryzacji.

## 8. Czego ten cykl NIE robi

- Nie liczy metryki g_mu_nu, geodezyjnych ani ukladow wielocialowych.
- Nie stabilizuje obiektow do skonczonego rozmiaru (obiekty rosna; pomiar
  w zdefiniowanym oknie zanim fronty sie zetkna) — stabilizacja to
  osobny watek (op-wall-dynamics).
- Nie rozstrzyga galezi B1/B2 — tylko wskazuje galaz.
- Nie dotyka core/ ani rejestrow predykcji.

## 9. Deliverables

- `Phase0_balance.md` (ten plik; po autoryzacji: status LOCKED)
- `Phase0_check.py` / `Phase0_output.txt` (weryfikacja rownan #1 i #3
  pre-code: m0, ogon K0, oszacowanie v_inf)
- `Phase1_front_calibration.py` / `Phase1_output.txt` (G2)
- `Phase2_pair_interaction.py` / `Phase2_output.txt` (G3, G4)
- `Phase3_mass_universality.py` / `Phase3_output.txt` (G5, G6)
- `Phase_FINAL_close.md`, `README.md`

## 10. Status

**LOCKED (2026-07-04).** Model, predykcje (w tym jawne oczekiwanie Yukawa
dla Z2!), kryteria i regula decyzyjna spisane przed kodem. Autoryzacja
autora udzielona 2026-07-04 (polecenie przeprowadzenia cyklu). Obowiazuje
pelny anty-Lakatos. Doprecyzowania metrologiczne pre-code: sekcja 11
(spisane PRZED pierwsza linijka kodu).

## 11. Doprecyzowania pre-code (spisane przed pierwsza linijka kodu)

Zadne z ponizszych nie zmienia progow G1-G6 ani modelu — to wylacznie
ujednoznacznienie procedur pomiarowych, zeby po runach nie bylo swobody
interpretacyjnej:

- **D1 (kompozycja pary):** przesuniecia obiektow o CALKOWITA liczbe
  wezlow (np.roll, periodycznie); szczelina `g` mierzona z kompozycji
  (liniowa interpolacja przeciec `Phi = eps` na osi y=0), raportowana
  jako `g_meas` zamiast nominalnego `g0`. Zero interpolacji pola.
- **D2 (fit Y vs P):** fit w przestrzeni log: `y = ln(-E_int)`;
  (Y): regresja liniowa `y ~ alpha - m*g`; (P): `y ~ beta - p*ln(g)`.
  RSS liczone w `y`; `AIC = n*ln(RSS/n) + 2k`, `k = 2` dla obu form;
  `Delta_AIC` wg G4. Punkty z `E_int >= 0` dyskwalifikuja fit (G3).
- **D3 (predkosci):** roznice centralne na zapisach co 10 krokow,
  baza `Delta_tau = 2.0`; `v_zblizania = -dg/dtau`.
- **D4 (kontrola dla Delta_v):** kontrola pierwotna = run POJEDYNCZEGO
  obiektu z TEGO SAMEGO snapshotu (identyczne R0, siatka, dt);
  `Delta_v(tau) = v_zblizania(tau) - 2*v_control(tau)`, parowanie po tau
  (rowne R poczatkowe => dopasowanie krzywizna do 1. rzedu). Prawo #3
  z Phase 1 sluzy jako odniesienie wtorne, nie pierwotne.
- **D5 (agregaty):** `Delta_v_bar` (do ratio G5/G6) = srednia `Delta_v`
  interpolowanego w `g in {2.0, 2.5, 3.0, 3.5}`. Warunek dynamiczny G3:
  `Delta_v > 0` w `g in {2.0, 3.0, 4.0}` na runie `g0 = 6`.
- **D6 (boundary_contact):** `Phi > eps` przy odleglosci Chebysheva
  `> 0.45*L` od centrum siatki.
- **D7 (historie G6):** hist1 = `cluster(A0=0.6, D=3, K=7)`,
  hist2 = `cluster(A0=1.0, D=3, K=7)`; obie ewoluowane do tego samego
  `R_obj = 8` (rozny czas ewolucji).
- **D8 (promien obiektu):** polowa szerokosci `Phi > eps` wzdluz osi x
  przez centrum, interpolacja subgridowa (spojnie z sekcja 2).
- **D9 (okno fitu #3):** kalibracja frontu fitowana w `R in [6, 25]`
  (odciecie transjentu formowania locku i strefy brzegowej).
- **D10 (monotonicznosc H):** `H_Gamma` proba co 10 krokow (kadencja
  zapisow), tolerancja `1e-10 * max(1, |H0|)` jak w poprzednim opie.
