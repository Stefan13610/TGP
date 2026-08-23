---
title: "Phase0_balance — addytywnosc parowa oddzialywania lockow w ukladach N-cialowych"
date: 2026-07-04
type: phase0-lock
tgp_owner: research/op-nbody-additivity-2026-07-04
status: LOCKED
anti_lakatos_lock: LOCKED (2026-07-04, autoryzacja autora — polecenie "dzialaj z takim testem"; lock przed pierwsza linijka kodu)
tags: [gravity-bridge, nbody, additivity, superposition, yukawa, fermat-point]
related:
  - "[[../op-lock-interaction-gravity-2026-07-04/Phase_FINAL_close.md]]"
  - "[[../op-bare-substrate-genesis-2026-07-04/Phase_FINAL_close.md]]"
  - "[[../../core/sek08_formalizm/sek08_formalizm.tex]]"
---

# Phase 0 — LOCK: addytywnosc parowa w ukladach N-cialowych

## 0. Hipoteza i po co to rdzeniowi

Po `op-lock-interaction-gravity` (PASS, galaz B): oddzialywanie 2-cialowe
lockow jest Yukawa `E_2(g) ~ -A e^(-m0 g)`, `m0 = sqrt(a/kappa) = 1`.
Newtonowska granica rdzenia TGP (tw. Newtona, sek08_formalizm) zaklada
**liniowa superpozycje** zrodel (`delta_Phi` od wielu mas sie dodaje).
Genesis pokazal zarazem, ze nakladanie ZARODKOW jest silnie nieliniowe
(kolektywny lock, przewaga >= 7x). Pytanie tego cyklu:

> Czy dla ZALOCKOWANYCH obiektow (nasycone domeny, cienkie sciany)
> oddzialywanie N-cialowe jest suma parowa? Gdzie i jakim prawem
> addytywnosc sie zalamuje — i czy prog zalamania jest algebraicznie
> przewidywalny z substratu?

Hipoteza (pre-code): addytywnosc trzyma w rezimie rozseparowanym
(g >= 3), a wiodaca nieaddytywnosc jest czlonem 3-cialowym z kubicznej
czesci potencjalu (`V'''(0) = -2b < 0`), zlokalizowanym w punkcie
Fermata trojki obiektow — czyli WYPROWADZALNYM, nie dowolnym.

## 1. Predykcje algebraiczne zapisane PRZED kodem

Ogon obiektu (surface r=R): `s_i(x) ~ eta * exp(-m0*(|x-c_i|-R))`,
`eta ~ sqrt(eps) = 0.548`, `m0 = 1`, `s_bar = 0.4258`.

### Rownanie #1 — wiodacy czlon 3-cialowy

Rozwiniecie `V(s1+s2+s3) - suma czesci parowych` zostawia czysto
3-cialowe czlony; najnizszy: `-2b * s1*s2*s3` (z `-(b/3)|s|^3`).

```text
E_3body ~ -2b * eta^3 * J3(geometria),
J3 = calka exp(-m0*(r1s+r2s+r3s)) d^2x,  r_is = max(0, |x-c_i|-R)
zdominowana przez punkt FERMATA ukladu trzech powierzchni.
```

### Rownanie #2 — skalowanie nieaddytywnosci (trojkat rownoboczny)

Dla trojkata o szczelinie g (d = 2R+g): odleglosc srodek-powierzchnia
`rho_c = d/sqrt(3) - R`, wiec `J3 ~ exp(-3*m0*rho_c)`, a suma parowa
`~ 3A*exp(-m0*g)`. Stad:

```text
delta_3(g) = [E_N - suma_par E_2] / |suma_par E_2|
           ~ -C * exp(-mu*g),   mu = (sqrt(3) - 1) * m0 ~ 0.73
ZNAK: delta_3 < 0 (nadaddytywne przyciaganie, bo V'''(0) = -2b < 0)
RZAD WIELKOSCI: |delta_3(2)| ~ 10^-3..10^-2 (prefaktor niepewny,
raportowany, poza kryterium)
```

### Rownanie #3 — dyskryminator geometryczny (lancuch)

Dla wspolliniowego lancucha 3 obiektow punkt spotkania trzech ogonow
lezy za obiektem srodkowym: `J3 ~ exp(-m0*(2R+2g))` — zaniedbywalne.

```text
|delta_chain(g)| << |delta_tri(g)|; predykcja: ponizej podlogi 1e-4.
```

Nieaddytywnosc jest wiec GEOMETRYCZNA (Fermat), nie "od liczby cial".

## 2. Model i parametry (LOCKED — ciaglosc lancucha, zero strojenia)

```text
kappa=0.50, a=0.50, b=1.60, c=1.00, Phi=s^2; ds/dtau = kappa*Lap(s)-V'(s)
grid N=256, L=128 (dx=0.5), periodyczne, Laplasjan 5-pkt, dt=0.02
eps=0.30, S_GUARD=10 (flaga), seed=20260704
obiekt: cluster(A0=0.6, D=3, K=7) ewoluowany do R_obj=8 (jak poprzednio)
kompozycje statyczne: SUROWA superpozycja, przesuniecia o pelne wezly
  (D1 poprzedniego opa); E_2 liczone dla DOKLADNIE tych samych ofsetow
  calkowitoliczbowych, ktore wystepuja w konfiguracji N (eliminacja
  bledu pozycjonowania z konstrukcji); E_2(Delta) = E_2(-Delta)
konfiguracje: chain3(g), triangle3(g), hex6(g); d = 2*R_act + g
  triangle: wierzcholki co 120 stopni, circumradius d/sqrt(3)
  hex: wierzcholki co 60 stopni (katy +/-30,+/-90,+/-150), circumradius d
skan statyczny: g in {2.0, 2.5, 3.0, 4.0}
skan G2 (kontrola ciaglosci): pary osiowe g in {2, 2.5, 3, 4, 5, 6}
podloga numeryczna: |delta| < 1e-4 raportowane jako "ponizej podlogi",
  wykluczone z fitow G4 (zapisane pre-code)
dynamika: pair(g0=4) vs chain3(g0=4); kontrola = pojedynczy obiekt
  (D4 poprzedniego opa); Delta_v_bar = srednia Delta_v w g in
  {2.0,2.5,3.0,3.5} na czesci wspolnej nosnikow (D5); predkosci:
  roznice centralne, baza tau=2.0, zapis co 10 krokow (D3, D10)
obserwacja (poza kryteriami): hex6(g=3) dynamicznie vs pair(g0=3) —
  czas domkniecia szczeliny, symetria, mostki Phi>eps
okno pomiarowe dynamiki: do g < 1.0 albo boundary_contact (Chebyshev
  0.45L od centrum siatki, D6)
```

## 3. Obserwable

```text
- E_2(g_meas) osiowe + refit Yukawy (ciaglosc z poprzednim cyklem)
- E_N(konfiguracja, g), suma parowa, delta_N(g) dla chain3/tri/hex6
- znak delta_tri, monotonicznosc |delta_tri|(g), zbocze mu z fitu
  ln|delta_tri| vs g (punkty powyzej podlogi)
- ratio_coop = Delta_v_bar(szczelina chain3) / Delta_v_bar(para)
- H_Gamma nierosnace; determinizm; flagi finite/runaway/contact
```

## 4. Kryteria (LOCKED)

### G1 — techniczne
Determinizm (jeden run x2 bitowo), brak NaN/runaway, H_Gamma nierosnace
na kazdym runie ewolucyjnym. PASS/FAIL.

### G2 — ciaglosc prawa dwucialowego
Refit Yukawy na parach osiowych (g_meas z kompozycji):
`m_fit in [0.8, 1.2]` i `R^2(log) >= 0.97`.
FALSYFIKATOR: prawo 2-cialowe niereprodukowalne -> baza addytywnosci
niezdefiniowana, cykl nie moze orzekac.

### G3 — istnienie rezimu addytywnego
Dla g >= 3: `|delta_N(g)| <= 0.02` dla WSZYSTKICH konfiguracji
{chain3, triangle3, hex6}.
FALSYFIKATOR: nieaddytywnosc w slabym polu -> zalozenie superpozycji
zrodel w newtonowskiej granicy rdzenia BEZ poparcia w substracie;
przed B1/B2 konieczny osobny program wielocialowy.

### G4 — RDZEN: zalamanie addytywnosci jest algebraiczne (Fermat)
Na punktach |delta_tri| powyzej podlogi 1e-4:
(a) `delta_tri < 0` na calym skanie (znak z rownania #1);
(b) `|delta_tri(g)|` scisle rosnace przy malejacym g;
(c) fit `ln|delta_tri| = ln C - mu*g`: `mu in [0.4, 1.1]`
    (predykcja punktowa: (sqrt(3)-1)*m0 ~ 0.73);
(d) dyskryminator geometryczny: `|delta_chain(g)| < |delta_tri(g)|`
    dla kazdego wspolnego g (predykcja: chain ponizej podlogi).
Wymagane (a)+(b)+(d); (c) w pasmie. Jesli punktow powyzej podlogi < 3,
(c) raportowac jako NIEROZSTRZYGALNE (za slaby sygnal), G4 oceniac
z (a)+(b)+(d) na dostepnych punktach; przy ZERO punktow powyzej podlogi
G4 = NIEROZSTRZYGNIETE, raportowac wprost (addytywnosc "za dobra" na
pomiar mechanizmu zalamania).
FALSYFIKATOR: znak dodatni / brak monotonicznosci / mu poza pasmem /
chain >= tri -> mechanizm zalamania NIE jest ogonowo-geometryczny ->
"przewidywalnosc progu" pada.

### G5 — klasyfikator dynamiczny (bez PASS/FAIL)
`ratio_coop in [0.9, 1.1]` -> LOKALNY (sciana-ogon; addytywny w pedzie)
`ratio_coop >= 1.3`        -> KOLEKTYWNY (wzmocnienie ponadparowe)
posrednie                  -> MIXED, raport wprost.

## 5. Regula decyzyjna

```text
G3 PASS + G4 PASS/(a,b,d) -> addytywnosc parowa uzasadniona w rezimie
  rozseparowanym, zalamanie przewidywalne (Fermat). Zalozenie
  superpozycji w rdzeniu POPARTE w toy-skali. Nastepny krok: B1/B2
  bez modyfikacji; nieaddytywnosc parametryzowac jako poprawke
  ~exp(-mu*g) w przyszlych opach wielocialowych.
G3 PASS + G4 FAIL -> addytywnosc jest, mechanizm zalamania nieczytelny:
  raport; ekstrapolacje N-cialowe tylko empiryczne.
G3 FAIL -> superpozycja pada w slabym polu: STOP przed B1/B2,
  osobny program wielocialowy (mapa delta(N, g)) przed mostem.
G5: informacja dla rownan pedu (lokalna odpowiedz sciany vs dzialanie
  kolektywne); nie zmienia galezi.
```

## 6. Forbidden moves

1. Zero relaksacji konfiguracji statycznych przed pomiarem E (surowa
   superpozycja — to ona jest testowana).
2. Nie zmieniac progow G1-G5, parametrow substratu ani procedur po
   pierwszym runie (wyjatki techniczne: jedno obnizenie dt; jeden
   kontrolny run na drobniejszej siatce).
3. Nie dofitowywac innych form zalamania niz exp(-mu*g) po obejrzeniu
   danych (inne formy = nowy cykl).
4. Nie interpretowac dynamiki nadtlumionej jako orbit; tau != czas.
5. Podloge 1e-4 stosowac tak, jak zapisano — nie przesuwac po runie.
6. Nie promowac do core/ bez osobnej autoryzacji.

## 7. Czego cykl NIE robi

- Nie testuje N > 6 ani konfiguracji przypadkowych (mapa delta(N,g)
  to osobny op, jesli G3 FAIL).
- Nie rozstrzyga B1/B2; dostarcza warunek waznosci superpozycji dla
  granicy newtonowskiej rdzenia (sek08: delta_Phi ~ suma zrodel).
- Nie liczy sektora pedowego (dynamika nadtlumiona — pomiar statyki
  i predkosci scian).

## 8. Deliverables

- `Phase0_balance.md` (ten plik), `Phase0_check.py`/`Phase0_output.txt`
  (numeryczna weryfikacja rownan #1-#3: J3, mu_pred, rzad C)
- `Phase1_static_additivity.py` / `Phase1_output.txt` (G1-G4)
- `Phase2_dynamics_cooperativity.py` / `Phase2_output.txt` (G5 + hex obs.)
- `Phase_FINAL_close.md`, `README.md`
