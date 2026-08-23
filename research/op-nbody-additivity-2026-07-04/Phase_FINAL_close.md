---
title: "Phase_FINAL_close — zamkniecie: addytywnosc parowa oddzialywania lockow w ukladach N-cialowych"
date: 2026-07-04
type: phase-final-close
tgp_owner: research/op-nbody-additivity-2026-07-04
status: CLOSED
verdict: "PASS (G1-G4); G5 klasyfikator: LOKALNY — superpozycja parowa POPARTA w rezimie rozseparowanym"
anti_lakatos_lock: PRESERVED
tags: [gravity-bridge, nbody, additivity, superposition, fermat-point, closed]
related:
  - "[[Phase0_balance.md]]"
  - "[[README.md]]"
  - "[[../op-lock-interaction-gravity-2026-07-04/Phase_FINAL_close.md]]"
  - "[[../../core/sek08_formalizm/sek08_formalizm.tex]]"
---

# Phase FINAL — zamkniecie cyklu

## 0. Werdykt

**G1 PASS, G2 PASS, G3 PASS, G4 PASS (mechanizm Fermata potwierdzony
w calosci: znak, monotonicznosc, zbocze, geometria), G5 = LOKALNY.**

Wg zalockowanej reguly decyzyjnej: **addytywnosc parowa oddzialywania
zlockowanych obiektow jest uzasadniona w rezimie rozseparowanym
(g >= 3: |delta| <= 3.3e-3), a jej zalamanie jest algebraicznie
przewidywalne** — niesione przez czlon 3-cialowy `-2b*s1*s2*s3`
zlokalizowany w punkcie Fermata, o skalowaniu `~exp(-mu*g)`.
Zalozenie liniowej superpozycji zrodel w newtonowskiej granicy rdzenia
TGP (sek08, tw. Newtona) jest POPARTE w toy-skali. Galaz B1/B2 mostu
do grawitacji moze isc bez modyfikacji.

## 1. Co testowano

Po `op-lock-interaction-gravity` (2-cialowe prawo Yukawa, m0=1) —
pytanie: czy E_int ukladow N-cialowych jest suma parowa, gdzie sie
zalamuje i czy prog zalamania jest wyprowadzalny z substratu. Substrat
i obiekt identyczne jak w poprzednich cyklach (zero strojenia).
Konfiguracje: chain3 / triangle3 / hex6, g in {2, 2.5, 3, 4}; E_2
liczone dla dokladnie tych samych calkowitoliczbowych ofsetow co w
konfiguracji N (delta izoluje czysto >=3-cialowe czlony: czlony
kwadratowe, w tym gradientowe, kasuja sie tozsamosciowo).

## 2. Przebieg

- **Phase 0**: LOCK (autoryzacja autora 2026-07-04) + weryfikacja
  rownan #1-#3 na modelu ogonowym: mu_pred(model, skonczone R) = 0.59
  vs asymptotyczne (sqrt(3)-1)*m0 = 0.73; |delta_chain/delta_tri| ~
  1.6e-5; znak ujemny. `Phase0_output.txt` PASS.
- **Phase 1**: statyka — G1-G4. `Phase1_output.txt`.
- **Phase 2**: dynamika — G5 + obserwacja hex6. `Phase2_output.txt`.

## 3. Wyniki liczbowe (kluczowe)

### 3.1 Ciaglosc prawa 2-cialowego (G2)

```text
refit Yukawy (pary osiowe, g_meas z kompozycji):
  m = 0.813 in [0.8, 1.2],  R^2(log) = 0.979 >= 0.97  -> PASS
(por. m = 0.883/0.984 z poprzedniego cyklu — spojnie z m0 = 1)
```

### 3.2 Addytywnosc i jej zalamanie (G3, G4 — RDZEN)

```text
delta_N(g) = [E_N - suma_par E_2] / |suma_par E_2|

g     chain3        triangle3     hex6
2.0   -1.3e-09      -6.56e-03     -1.8e-06
2.5   -9.4e-10      -3.93e-03     -8.9e-07
3.0   -6.1e-10      -3.30e-03     -5.6e-07
4.0   -2.2e-10      -1.56e-03     -1.8e-07

G3: max|delta| przy g>=3: 3.3e-3 <= 0.02                    -> PASS
G4a: znak delta_tri < 0 wszedzie (predykcja V'''(0) = -2b)  -> PASS
G4b: |delta_tri| scisle rosnie przy malejacym g             -> PASS
G4c: mu_fit = 0.689 in [0.4, 1.1]
     (predykcja Fermata: 0.73 asympt., 0.59 model skoncz. R) -> PASS
G4d: |delta_chain| < |delta_tri| o ~6 rzedow                -> PASS
```

Nieaddytywnosc jest GEOMETRYCZNA, nie "od liczby cial": trojkat
(punkt Fermata miedzy trzema powierzchniami) niesie sygnal 6 rzedow
wielkosci wiekszy niz lancuch o tej samej liczbie obiektow; hexagon
(rozwarte trojki) lezy 3-4 rzedy ponizej trojkata. Wszystkie trzy
przewidziane pre-code z jednego czlonu `-2b*s1*s2*s3`.

### 3.3 Dynamika (G5) i obserwacja hex6

```text
Delta_v(g=3) przy g0=4:  para = +0.32935,  chain3 = +0.32935
ratio_coop = 1.0000  ->  LOKALNY (sciana-ogon, addytywny)
hex6(g=3) vs pair(g=3): tau(g<=1.5) = 1.80 vs 1.80 (stosunek 1.0000)
```

Zgodnosc do rozdzielczosci wydruku jest przewidziana: wplyw trzeciego
ciala na cudza szczeline tlumiony ~exp(-(2R+2g)) ~ 1e-11. Sciana
odpowiada na LOKALNY ogon po drugiej stronie wlasnej szczeliny —
brak "dzialania kolektywnego" w predkosciach scian.

### 3.4 G1 (techniczne)

Determinizm bitowy (prep x2), zero NaN/runaway, H_Gamma nierosnace
na wszystkich runach ewolucyjnych, zero boundary_contact. PASS.

## 4. Interpretacja w granicach modelu (kontekst rownan pedu TGP)

CO WYNIK ZNACZY:
- dla ZALOCKOWANYCH obiektow (nasycone domeny, cienkie sciany)
  superpozycja parowa dziala w slabym nakladaniu z dokladnoscia
  <=0.3% (g>=3) — inaczej niz dla ZARODKOW (genesis: kolektywnosc),
  bo locki sa nasycone, a nieliniowosc widzi tylko iloczyny ogonow,
- wiodaca poprawka nieaddytywna jest w pelni scharakteryzowana:
  znak (nadaddytywne przyciaganie), nosnik (punkt Fermata), skalowanie
  (`~e^(-mu g)`, mu ~ 0.69) — mozna ja parametryzowac jako poprawke
  wyzszego rzedu w przyszlych opach wielocialowych,
- dla granicy newtonowskiej rdzenia (delta_Phi ~ suma zrodel,
  sek08_formalizm): toy-substrat NIE podwaza superpozycji w rezimie
  rozseparowanym; jej zalamanie wchodzi wykladniczo dopiero przy
  geometriach z blisko zbiegajacymi sie trzema powierzchniami,
- G5 LOKALNY: w dynamice nadtlumionej "sila" na sciane jest lokalnym
  odczytem ogona — pedu nie przenosi zadne dzialanie kolektywne;
  spojne z obrazem, w ktorym rownania pedu buduje sie z pola, nie
  z par czastek.

CZEGO WYNIK NIE ZNACZY:
- nie testowano N > 6, konfiguracji przypadkowych ani rezimu g < 2
  (tam kompozycja przestaje byc superpozycja — obiekty sie mergu-ja),
- tau != czas fizyczny; zero wnioskow o orbitach/geodezyjnych,
- wynik NIE jest dowodem poprawnosci granicy newtonowskiej rdzenia —
  jest usunieciem konkretnego zastrzezenia (superpozycja) w toy-skali.

## 5. Errata / transparencja

1. **G2, tabela par osiowych:** zaokraglenie polozen do pelnych wezlow
   dalo IDENTYCZNA kompozycje dla g0=2.5 i g0=3.0 (ten sam ofset
   19 komorek) — fit zawiera wiec 5 roznych punktow (jeden zdublowany),
   w tym punkt g0=2 z g_meas=1.01 (poza deklarowana waznoscia
   superpozycji, jak w poprzednim cyklu). Kryterium ocenione tak, jak
   zalockowano; m=0.813 lezy w pasmie przy dolnej krawedzi — spojne
   z wynikiem poprzedniego cyklu (0.883; czulosc 0.984).
2. **Phase 2, identycznosc wydrukow:** roznice para/chain3 i
   para/hex6 leza ponizej rozdzielczosci wydruku (przewidziana
   supresja ~1e-11), a nie sa artefaktem wspolnego kodu — konfiguracje
   komponowane niezaleznie (3 i 6 obiektow), mierzone w roznych
   miejscach siatki.

Zadne progi, parametry ani procedury nie byly zmieniane po pierwszym
runie. Wyjatki techniczne (dt, drobna siatka) nie byly potrzebne.
Podloga 1e-4 nie byla przesuwana (wszystkie punkty trojkata nad nia).

## 6. Decyzja i nastepny krok

Regula decyzyjna: G3 PASS + G4 PASS -> **superpozycja parowa
uzasadniona; zalamanie przewidywalne (Fermat)**. Program wraca na
glowna os mostu do grawitacji — drafty przygotowane
([[../op-goldstone-mediator-2026-07-04/Phase0_balance.md]],
[[../op-phi-metric-refraction-2026-07-04/Phase0_balance.md]]),
do autoryzacji autora jako osobne opy:

- **B1**: pole zespolone `s in C`, U(1); mod Goldstone'a jako bezmasowy
  mediator dlugozasiegowy (rdzen juz pisze `s_i -> psi_i = |psi|e^(i
  theta)` w sek01_ontologia — B1 jest testem tej struktury),
- **B2**: `Phi` jako metryka efektywna propagacji malych zaburzen
  (droga geometryczna; odpowiada mechanizmowi `g_00 = -(1-2 delta_Phi/
  Phi_0)` + geodezyjne z sek08_formalizm).

W przyszlych opach wielocialowych nieaddytywnosc parametryzowac jako
`delta ~ -C e^(-mu g)`, mu ~ 0.7, C ~ 3e-2 (z tego cyklu).

## 7. Pliki cyklu

- [[Phase0_balance.md]] — LOCK (predykcje #1-#3, kryteria G1-G5)
- `Phase0_check.py` / `Phase0_output.txt` — weryfikacja rownan pre-code
- `Phase1_static_additivity.py` / `Phase1_output.txt` — G1-G4
- `Phase2_dynamics_cooperativity.py` / `Phase2_output.txt` — G5 + hex6
- [[README.md]] — przeglad opa
- ten plik — zamkniecie
