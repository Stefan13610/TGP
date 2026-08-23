---
title: "Phase_FINAL_close — zamkniecie: oddzialywanie zlockowanych obiektow (most do grawitacji, krok 1)"
date: 2026-07-04
type: phase-final-close
tgp_owner: research/op-lock-interaction-gravity-2026-07-04
status: CLOSED
verdict: "PASS (G1-G4, G6); G5 klasyfikator: FRONTOWE -> GALAZ B"
anti_lakatos_lock: PRESERVED
tags: [gravity-bridge, lock-interaction, yukawa, force-law, branch-B, closed]
related:
  - "[[Phase0_balance.md]]"
  - "[[README.md]]"
  - "[[../op-bare-substrate-genesis-2026-07-04/Phase_FINAL_close.md]]"
  - "[[../op-wall-dynamics-2026-07-03/README.md]]"
---

# Phase FINAL — zamkniecie cyklu

## 0. Werdykt

**G1 PASS, G2 PASS, G3 PASS, G4 PASS (forma Yukawa, m_fit w pasmie
predykcji), G5 = FRONTOWE (klasyfikator), G6 PASS.**

Wg zalockowanej reguly decyzyjnej (Phase0_balance.md, sekcja 6):

> G3+G4 PASS, forma (Y) Yukawa (oczekiwana), G5 frontowe -> **GALAZ B**

Bezposrednie oddzialywanie miedzy zlockowanymi obiektami ISTNIEJE, jest
przyciagajace, algebraicznie wyprowadzalne z substratu — ale jest
**krotkozasiegowe (Yukawa, m ~ m0 = 1) i skaluje sie z frontem, nie z
"masa"**. NIE jest kandydatem na grawitacje. Zgodnie z LOCKiem to
pelnoprawne rozstrzygniecie, nie porazka: most do grawitacji przechodzi
na mechanizmy B1 (bezmasowy mediator, pole zespolone/Goldstone) albo
B2 (Phi jako metryka efektywna propagacji — droga geometryczna).

Wynik jest DOKLADNIE zgodny z oczekiwaniem teoretycznym zapisanym przed
kodem (sekcja 1.2 LOCKa): pole Z2 nie ma modu bezmasowego -> Yukawa.

## 1. Co testowano

Hipoteza autora (po PASS `op-bare-substrate-genesis`): grawitacja jako
oddzialywanie pol miedzy zlockowanymi obiektami, policzalne dla prostych
przypadkow. Wersja minimalna: dwa zlockowane obiekty w tym samym
substracie — czy oddzialuja, jakim prawem, czy prawo ma cechy grawitacji
(dlugi zasieg, skalowanie masowe) czy sily powierzchniowej (krotki
zasieg, skalowanie frontowe). Model identyczny z poprzednim cyklem
(zero strojenia): `kappa=0.5, a=0.5, b=1.6, c=1.0, Phi=s^2`, przeplyw
gradientowy, N=256, L=128, dt=0.02.

## 2. Przebieg

- **Phase 0**: LOCK (autoryzacja autora 2026-07-04) + doprecyzowania
  metrologiczne D1-D10 (sekcja 11, przed pierwsza linijka kodu).
  Weryfikacja rachunkow: `Phase0_output.txt` — m0=1.0, ogon K0
  (zbocze -0.989), front 1D v_planar=0.1599 vs teoria DeltaV/int(f'^2)
  =0.1608 (0.6%). PASS.
- **Phase 1**: kalibracja frontu (G2). `Phase1_output.txt`.
- **Phase 2**: statyka E_int(g) + dynamika pary (G3, G4).
  `Phase2_output.txt`.
- **Phase 3**: skalowanie masowe (G5) + uniwersalnosc historii (G6).
  `Phase3_output.txt`.

## 3. Wyniki liczbowe (kluczowe)

### 3.1 Kalibracja frontu (G2) — rownanie #3

```text
fit v(R) = v_inf - c/R  (787 punktow, R in [6, 25]):
  v_inf = 0.15993   (front 1D z Phase 0: 0.15986 — zgodnosc < 0.1%)
  c     = 0.4898    (predykcja LOCK: c ~ kappa = 0.5)
  R^2   = 0.9932    (prog >= 0.98)
```

### 3.2 Prawo oddzialywania (G3, G4) — RDZEN

Statycznie, E_int(g) = H[s1+s2] - 2H[obj], R_obj = 8:

```text
g_meas:  1.010    2.711    3.886    4.939    5.956    7.964
E_int:  -0.534   -0.217   -0.0816  -0.0295  -0.0105  -0.00129
```

Przyciaganie na calym skanie, |E_int| scisle malejace. Fity (D2):

```text
(Y) Yukawa:  A=1.940, m=0.883, R^2(log)=0.984, AIC=-12.50
(P) potega:  B=1.317, p=2.691, R^2(log)=0.816, AIC=+2.10
Delta_AIC (P-Y) = +14.60  >= 10  ->  WYGRYWA YUKAWA
m_fit = 0.883 in [0.8, 1.2]  ->  zgodne z predykcja substratu m0 = 1.0
[czulosc, bez punktu g_meas=1.01 (superpozycja wazna wg LOCKa dla g>=2):
 m = 0.984, R^2(log)=0.998, Delta_AIC=+15.62 — jeszcze blizej m0]
```

Dynamicznie (kontrola D4 — ten sam snapshot, pojedynczy): Delta_v > 0
na WSZYSTKICH zmierzonych punktach; przekroj (run g0=8):

```text
g:        3.0       3.5       4.0       5.0       6.0       7.0
Delta_v: +0.333    +0.175    +0.101    +0.033    +0.011    +0.004
fit informacyjny: m_dyn = 1.147 (R^2=0.992), Delta_AIC(P-Y)=+74.6
```

Statyka i dynamika NIEZALEZNIE daja to samo prawo: eksponencjalny zanik
z masa ~ 1 = sqrt(a/kappa). Prawo oddzialywania jest ALGEBRAICZNIE
WYPROWADZONE z substratu (rownanie #1 LOCKa), nie dofitowane.

### 3.3 Dyskryminator masowy/frontowy (G5) i uniwersalnosc (G6)

```text
Delta_v(g=3) przy g0=4:  8x8: +0.329   12x12: +0.359   16x16: +0.376
ratio_mass = 1.142  (12/8: 1.091)   <= 1.5  ->  FRONTOWE (powierzchniowe)
ratio_hist = 0.989, |ratio-1| = 0.011 <= 0.10  ->  G6 PASS (wersja mocna)
```

Sila rosnie z rozmiarem zaledwie ~sqrt-owo/slabiej (2x wiekszy obiekt ->
+14% sily), czyli skalowanie POWIERZCHNIOWE, nie objetosciowe/masowe.
Oddzialywanie zalezy wylacznie od stanu pola (rozne historie przygotowania
przy tym samym R daja identyczna sile w granicach 1.1%).

### 3.4 G1 (techniczne)

Determinizm bitowy (Phase 1, run x2), zero NaN/runaway, H_Gamma
nierosnace na KAZDYM runie (H_inc_max = 0.00e+00), zero clampow,
zero boundary_contact. PASS.

## 4. Interpretacja w granicach modelu

CO WYNIK ZNACZY:
- zlockowane obiekty NIE sa "gluche" — oddzialuja przyciagajaco przez
  ogony pola w falszywej prozni,
- prawo oddzialywania jest policzalne algebraicznie z substratu:
  E_int(g) ~ -A e^(-m0 g), m0 = sqrt(V''(0)/kappa) — hipoteza
  "policzalnosci prostych przypadkow" POTWIERDZONA,
- ale to oddzialywanie ma zly profil na grawitacje: zasieg ~ 1/m0 = 1
  (krotki), skalowanie frontowe — jest analogiem sily powierzchniowej /
  oddzialywania kink-antykink, nie sily Newtona.

CZEGO WYNIK NIE ZNACZY:
- nie falsyfikuje mostu do grawitacji jako takiego — falsyfikuje JEDEN
  z trzech mechanizmow (sila bezposrednia); mediator bezmasowy (B1)
  i droga geometryczna (B2) pozostaja otwarte i sa nastepnym krokiem,
- dynamika jest nadtlumiona: "sila" = predkosc zblizania, zadnych
  orbit/geodezyjnych/inercji (forbidden move 4),
- tau nie jest czasem fizycznym; zero deklaracji zgodnosci z GR/MOND
  (forbidden move 5).

## 5. Errata / transparencja

1. **G3-dynamika (raportowanie, dzien runow, przed interpretacja):**
   pierwotna linia ewaluacji sprawdzala dodatniosc Delta_v przez
   interpolacje w punktach {2,3,4} (doprecyzowanie D5); punkt g=2.0
   lezy poza nosnikiem estymatora roznic centralnych (obciecie
   +/-tau=2.0 przy przyspieszajacych frontach tuz przed zetknieciem),
   przez co "niemierzalne w g=2.0" zlewalo sie z "niedodatnie".
   Skorygowano do brzmienia LOCKa G3 ("Delta_v(g) > 0 dla g <= 4"):
   wymagane istnienie i dodatniosc WSZYSTKICH zmierzonych punktow
   g <= 4 (17 punktow, wszystkie dodatnie, pokrycie do g=2.62).
   Analogiczny przypadek jak errata #1 poprzedniego opa.
2. **Punkt statyczny g0=2:** po kompozycji g_meas = 1.01 < 2 (ogony
   podnosza Phi w szczelinie) — poza deklarowana w LOCKu waznoscia
   superpozycji (g >= 2). Punkt pozostawiony w ficie glownym (LOCK
   definiuje skan po g0), a fit czulosci bez niego RAPORTOWANY:
   m = 0.984 — wniosek G4 niewrazliwy (Delta_AIC rosnie do +15.6).
3. Incydent techniczny bez wplywu na wyniki: pliki Phase1/Phase2
   zostaly poczatkowo zapisane w zagniezdzonym katalogu (dryf cwd
   narzedzia) i przeniesione do katalogu opa bez zmian tresci.

Zadne progi G1-G6, wspolczynniki ani procedury dynamiki nie byly
zmieniane po pierwszym runie. Dozwolone wyjatki techniczne (obnizenie
dt, kontrolny run na drobnej siatce) NIE byly potrzebne.

## 6. Decyzja i nastepny krok

Regula decyzyjna (LOCK, sekcja 6): **GALAZ B** — fundament stoi,
mechanizm mostu do wymiany. Dwie drogi, kazda jako osobny op
(drafty przygotowane: [[../op-goldstone-mediator-2026-07-04/Phase0_balance.md]]
i [[../op-phi-metric-refraction-2026-07-04/Phase0_balance.md]];
do autoryzacji przez autora):

- **B1 — bezmasowy mediator:** pole zespolone `s in C` z symetria U(1);
  mod fazowy (Goldstone) daje oddzialywanie dlugozasiegowe; lock z
  op-bare-substrate-genesis powtorzyc dla `|s|`.
- **B2 — droga geometryczna:** `Phi` jako metryka efektywna — predkosc
  propagacji malych zaburzen zalezna od `Phi(x)`; "grawitacja" =
  refrakcja trajektorii w gradiencie `Phi`, nie sila miedzy cialami.

Uwaga strukturalna z LOCKa podtrzymana: ten cykl jest POMIAREM, ktory
z trzech mechanizmow niesie grawitacje w tej ontologii — mechanizm
"sila bezposrednia" zostal zmierzony i odrzucony jako nosnik grawitacji,
przy jednoczesnym potwierdzeniu policzalnosci substratu.

## 7. Pliki cyklu

- [[Phase0_balance.md]] — LOCK + doprecyzowania pre-code (sekcja 11)
- `Phase0_check.py` / `Phase0_output.txt` — weryfikacja rownan #1 i #3
- `Phase1_front_calibration.py` / `Phase1_output.txt` — G2 (+G1)
- `Phase2_pair_interaction.py` / `Phase2_output.txt` — G3, G4
- `Phase3_mass_universality.py` / `Phase3_output.txt` — G5, G6
- [[README.md]] — przeglad opa
- ten plik — zamkniecie
