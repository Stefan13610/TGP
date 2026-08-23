---
title: "Phase0_balance — DRAFT: bezmasowy mediator (pole zespolone U(1), mod Goldstone'a) — most do grawitacji, galaz B1"
date: 2026-07-04
type: phase0-draft
tgp_owner: research/op-goldstone-mediator-2026-07-04
status: LOCKED
anti_lakatos_lock: LOCKED (2026-07-04, autoryzacja autora — polecenie przeprowadzenia cyklu; lock przed pierwsza linijka kodu)
tags: [gravity-bridge, branch-B1, goldstone, u1, vortex, long-range, newton-2d]
related:
  - "[[../op-lock-interaction-gravity-2026-07-04/Phase_FINAL_close.md]]"
  - "[[../op-nbody-additivity-2026-07-04/Phase_FINAL_close.md]]"
  - "[[../../core/sek01_ontologia/sek01_ontologia.tex]]"
---

# Phase 0 — DRAFT: mod Goldstone'a jako bezmasowy mediator (galaz B1)

## 0. Hipoteza robocza

Po `op-lock-interaction-gravity` (galaz B: kanal amplitudowy Z2 jest
Yukawa, krotkozasiegowy) i `op-nbody-additivity` (superpozycja parowa
poparta): jesli substratem jest pole ZESPOLONE `psi = |psi| e^(i theta)`
z symetria U(1) — jak pisze rdzen w sek01_ontologia
(`s_i -> psi_i = |psi_i| e^(i theta_i)`) — to w fazie metrycznej
(`Phi = |psi|^2 > 0`) mod fazowy `theta` jest bezmasowy (Goldstone)
i MOZE niesc oddzialywanie dlugozasiegowe:

> Czy mod Goldstone'a daje w 2D oddzialywanie logarytmiczne
> (sila ~ 1/r — dokladnie forma "2D Newton" z reguly decyzyjnej
> poprzedniego opa), z policzalna stala sprzezenia z substratu?
> I czy to oddzialywanie moze byc grawitacja (uniwersalne
> przyciaganie), czy jest ladunkowe (oba znaki)?

WAZNA reinterpretacja (zapisac wprost przy locku): nosnikami "zrodel"
sa tu DEFEKTY fazowe (wiry o krecie n) w spojnym tle prawdziwej prozni
`|psi| = s*`, a nie locki w falszywej prozni. Sens fizyczny: locki z
genesis rosna i wypelniaja przestrzen faza metryczna — defekty zyja
W wygenerowanej przestrzeni. Ogon amplitudowy pozostaje Yukawa
(zmierzony); B1 pyta o sektor fazowy.

## 1. Predykcje algebraiczne do zalockowania (PRZED kodem)

### 1.1 Rownanie #1 — profil rdzenia wiru

```text
xi = sqrt(kappa / V''(s*)) = sqrt(0.5/0.87869) = 0.754
|psi|(r) -> s* dla r >> xi;  |psi|(0) = 0 (rdzen = lokalne Phi=0)
```

### 1.2 Rownanie #2 — prawo oddzialywania pary wirow (RDZEN)

Czesc fazowa energii: `0.5*kappa*|psi|^2*|grad theta|^2` -> efektywna
sztywnosc `J = kappa * s*^2 = kappa * Phi*` (sprzezenie NIESIONE przez
gestosc wygenerowanej przestrzeni Phi* = 1.3787 — zapisac wprost):

```text
J = kappa * Phi* = 0.5 * 1.3787 = 0.6893
E_int(d) = -2*pi*J * n1*n2 * ln(d/xi) + const
para (+1,-1): E rosnie z d -> PRZYCIAGANIE, sila F = 2*pi*J/d ~ 1/d
   (to jest DOKLADNIE "2D Newton": potencjal logarytmiczny)
predykcja ilosciowa zbocza: C2 = 2*pi*J = 4.331
```

Dwie zalockowane formy konkurencyjne (rozstrzyga Delta_AIC >= 10):

```text
(L) logarytmiczna:  E_int(d) = C1 + C2*ln(d)      (dlugi zasieg)
(Y) Yukawa:         E_int(d) = C1 - A*e^(-m*d)    (krotki zasieg)
```

OCZEKIWANIE TEORETYCZNE (jawne, pre-code): wygrywa (L) z
`C2 in [0.8, 1.2] * 2*pi*J = [3.46, 5.20]`. Gdyby wygralo (Y) —
Goldstone jest masywny (artefakt siatki/pinning rdzenia albo model
nie robi tego, co obiecuje symetria) — galaz B1 falsyfikowana.

### 1.3 Rownanie #3 — znak i skalowanie z kretem

```text
n1*n2 = -1 (wir-antywir):  przyciaganie
n1*n2 = +1 (rownoimienne):  ODPYCHANIE   <- kluczowy dyskryminator
|n1*n2| = 4 ((+2,-2)):      zbocze x4 wzgledem (+1,-1)
```

JAWNE zastrzezenie interpretacyjne (pre-code): grawitacja jest
UNIWERSALNIE przyciagajaca; mediator U(1) daje oba znaki. Wynik
"(L) + oba znaki" oznacza: istnieje bezmasowy mediator dlugozasiegowy
o formie 2D-Newtona, ale LADUNKOWY — do grawitacji brakuje mechanizmu
uniwersalnosci znaku. To jest pelnoprawne rozstrzygniecie (patrz
regula decyzyjna), nie porazka.

### 1.4 Dynamika nadtlumiona

`v(d) ~ mobilnosc * 2*pi*J/d` -> `v*d ~ const` (poprawki logarytmiczne
mobilnosci wiru raportowane, poza kryterium). Zadnych orbit.

## 2. Model i parametry (propozycja do LOCKa)

```text
pole: psi(x) in C; H = int[ 0.5*kappa*|grad psi|^2 + V(|psi|) ]
V(m) = 0.5*a*m^2 - (b/3)*m^3 + 0.25*c*m^4  (m = |psi| >= 0)
kappa=0.50, a=0.50, b=1.60, c=1.00 (ciaglosc lancucha, zero strojenia)
dynamika: dpsi/dtau = kappa*Lap(psi) - V'(|psi|)*psi/|psi| (grad. flow)
grid: N=256, L=128 (dx=0.5), periodyczne, Laplasjan 5-pkt, dt=0.02
eps=0.30, S_GUARD=10 (flaga), seed=20260704
UWAGA techniczna do LOCKa: xi=0.75 ~ 1.5*dx — rdzen na granicy
  rozdzielczosci; obowiazkowy jeden kontrolny run dx=0.25 (N=512 albo
  L=64) dla przypadku granicznego (wyjatek techniczny jak w poprzednich
  opach, zapisany z gory)
neutralnosc: w pudle periodycznym calkowity kret = 0 (wszystkie
  konfiguracje neutralne; "pojedynczy wir" = para o d_ref = L/2)
przygotowanie pary (n1,n2): ansatz theta(x) = n1*atan2(y-y1, x-x1)
  + n2*atan2(y-y2, x-x2), |psi| = s* * prof(r1) * prof(r2),
  prof(r) = r/sqrt(r^2 + 2*xi^2); potem krotka relaksacja tau_relax
  (do zalockowania, np. 200 krokow) przed pomiarem statycznym
skan statyczny: d in {6, 8, 12, 16, 24} (d << L/2; poprawki obrazow
  periodycznych oszacowac w Phase0_check, raportowac)
dynamika: d0 in {12, 16, 24}; okno do d < 3 albo anihilacji
metrologia: przejac D1-D10 z op-lock-interaction (zapisy co 10 krokow,
  roznice centralne tau=2.0, podlogi raportowane) — potwierdzic przy locku
```

## 3. Scenariusze

1. `genesis_u1`: powtorka kluczowego wyniku genesis dla |psi| przy
   theta=const: bare zanika, cluster lockuje (ciaglosc fundamentu).
2. `vortex_core`: profil rdzenia po relaksacji vs rownanie #1.
3. `pair_static(d)`: E(d) pary (+1,-1); fit (L) vs (Y).
4. `pair_dynamic(d0)`: d(tau); test v*d ~ const; spojnosc znaku.
5. `sign_test`: konfiguracja neutralna 4-wirowa (+1,+1,-1,-1), pary
   rownoimienne blisko, przeciwne daleko: czy (+1,+1) sie ODPYCHA.
6. `winding`: para (+2,-2): zbocze/4 vs (+1,-1); obserwacja
   stabilnosci rdzenia n=2 (przewidywanie: rozszczepienie na 2x n=1).

## 4. Kryteria (propozycja do LOCKa)

- **G1 techniczne**: determinizm, brak NaN/runaway, H_Gamma nierosnace
  (przeplyw gradientowy), flagi kontaktu. PASS/FAIL.
- **G2 ciaglosc fundamentu**: genesis_u1 odtwarza jakosciowo genesis
  (bare: Phi->0; cluster: lock z Phi ~ Phi*). FALSYFIKATOR: U(1) psuje
  sam lock -> B1 zle postawione na tym substracie.
- **G3 istnienie i stabilnosc wirow**: po relaksacji rdzen zachowuje
  kret (n niezmienne), |psi| -> s* dla r >> xi, rdzen nie "przykleja"
  sie do siatki w sposob zmieniajacy klasyfikacje (kontrola dx=0.25).
- **G4 RDZEN — prawo dlugozasiegowe**: fit (L) vs (Y) na E(d):
  `Delta_AIC >= 10` dla (L) ORAZ `C2 in [3.46, 5.20]`
  (predykcja 2*pi*kappa*Phi* = 4.33 — sprzezenie z gestosci
  wygenerowanej przestrzeni). FALSYFIKATOR: (Y) wygrywa albo C2 poza
  pasmem -> mediator nie jest bezmasowy/policzalny.
- **G5 dyskryminator znaku (klasyfikator, bez PASS/FAIL)**:
  (+1,+1) odpycha i (+1,-1) przyciaga -> LADUNKOWY;
  oba przyciagaja -> UNIWERSALNY (niespodzianka, gravity-like);
  raport wprost.
- **G6 skalowanie z kretem**: zbocze (+2,-2) / zbocze (+1,-1)
  in [3.2, 4.8] (predykcja 4).

## 5. Regula decyzyjna (propozycja)

```text
G4 PASS (L, C2 w pasmie) + G5 LADUNKOWY
  -> mediator bezmasowy ISTNIEJE, forma 2D-Newtona, sprzezenie ~ Phi*,
     ale ladunkowy. Grawitacja NIE jest prostym Goldstone'em U(1).
     Nastepne kroki (osobne opy): (i) sektor skalarny/amplitudowy na
     tle Phi* (dylatacyjny, zawsze przyciagajacy?), (ii) B2 —
     geometria propagacji; (iii) analiza: czy uniwersalnosc znaku
     wymaga tensora (spin-2), jak w standardowym argumencie.
G4 PASS + G5 UNIWERSALNY
  -> mocny kandydat na most: nastepny op — G_eff z parametrow
     substratu, uklady wielocialowe (superpozycja juz poparta),
     odpowiednik precesji/soczewkowania w 2D.
G4 FAIL (Yukawa albo C2 poza pasmem)
  -> Goldstone masywny/nieczytelny: jesli kontrola dx=0.25 wskazuje
     pinning rdzenia -> jeden powtorzony run (wyjatek techniczny);
     jesli nie -> B1 falsyfikowane, priorytet dla B2.
G2/G3 FAIL -> raport, stop przed interpretacja (fundament).
```

## 6. Forbidden moves (propozycja)

1. Nie przypinac wirow (zero external pinning); kret czytany z pola.
2. Nie zmieniac progow/parametrow/procedur po pierwszym runie
   (wyjatki: jedno obnizenie dt; jeden kontrolny run dx=0.25 —
   zapisany z gory jako obowiazkowy dla przypadkow granicznych).
3. Nie dofitowywac trzeciej formy prawa (tylko (L) vs (Y)).
4. Nie interpretowac dynamiki nadtlumionej jako orbit; tau != czas.
5. Wynik "LADUNKOWY" raportowac jako rozstrzygniecie, nie porazke.
6. Neutralnosc kretu obowiazkowa; zadnych konfiguracji z n_tot != 0.
7. Nie promowac do core/ bez osobnej autoryzacji.

## 7. Czego ten cykl NIE robi

- Nie liczy metryki, geodezyjnych, spin-2 ani ukladow N>4 wirow.
- Nie rozstrzyga B2 (osobny draft: op-phi-metric-refraction).
- Nie laczy skali J z jednostkami fizycznymi (q rdzenia pozostaje
  fenomenologiczne — por. rem. o statusie q w sek08).

## 8. Deliverables

- `Phase0_balance.md` (po autoryzacji: LOCKED + ew. poprawki pre-code)
- `Phase0_check.py`/`Phase0_output.txt` (xi, J, C2_pred, profil rdzenia
  1D, oszacowanie poprawek obrazow periodycznych w E(d))
- `Phase1_genesis_u1.py` / `Phase1_output.txt` (G2)
- `Phase2_vortex_pair_law.py` / `Phase2_output.txt` (G3, G4)
- `Phase3_sign_and_winding.py` / `Phase3_output.txt` (G5, G6)
- `Phase_FINAL_close.md`, `README.md`

## 9. Status

**LOCKED (2026-07-04).** Model, predykcje (w tym jawne zastrzezenie
o ladunkowosci U(1)!), kryteria i regula decyzyjna spisane przed kodem.
Autoryzacja autora udzielona 2026-07-04 (polecenie przeprowadzenia
cyklu; wybor galezi B1 przed B2 zgodny z regula decyzyjna B2, ktora
w galezi FAIL sama wskazuje "priorytet B1"). Obowiazuje pelny
anty-Lakatos. Doprecyzowania metrologiczne pre-code: sekcja 10
(spisane PRZED pierwsza linijka kodu).

## 10. Doprecyzowania pre-code (spisane przed pierwsza linijka kodu)

Zadne z ponizszych nie zmienia progow G1-G6, modelu, parametrow ani
predykcji — to ujednoznacznienie procedur pomiarowych i poprawna
implementacja zalockowanego ansatzu w pudle periodycznym:

- **D1 (periodyzacja ansatzu fazy):** goly `atan2` NIE jest periodyczny
  — na szwie y = +/-L/2 zostawia skok fazy zalezny od d
  (dla d=24: ~0.74 rad), czyli sztuczna energie szwu O(10) niszczaca
  pomiar C2. Implementacja zalockowanego ansatzu w pudle periodycznym:
  `theta(z) = sum_i n_i * sum_{k=-2..+2} arg sinh(pi*(z - z_i - k*L)/L)`
  (forma zamknieta nieskonczonej sumy obrazow atan2 w y + jawna suma
  obrazow w x). Dla konfiguracji neutralnych `e^{i theta}` jest scisle
  periodyczne w y i periodyczne w x z bledem ~e^(-4pi) ~ 3e-6.
  Rownowazne ansatzowi z sekcji 2 w granicy L -> inf. Amplituda wg
  sekcji 2, odleglosci r_i w konwencji najblizszego obrazu.
  Diagnostyka szwu (max skok fazy na brzegach) raportowana.
- **D2 (pozycje wirow i d_meas):** kret czytany z pola: winding
  plakietowy z owinietych roznic faz (suma po plakietce / 2pi);
  rdzenie = plakietki |n| >= 0.5 grupowane po sasiedztwie (Chebyshev
  <= 3 wezly, torus); pozycja rdzenia = centroid deficytu Phi
  (waga max(0, 0.5*Phi* - Phi)) w oknie +/-6 wezlow wokol plakietek;
  `d_meas` = odleglosc min-image miedzy centroidami przeciwnych znakow,
  mierzona W TEJ SAMEJ chwili co E.
- **D3 (tau_relax):** dokladnie 200 krokow (tau = 4.0) przed pomiarem
  statycznym; E i d_meas z tego samego kroku; punkt kontrolny po 150
  krokach raportowany (stabilnosc C2 wzgledem transjentu; poza
  kryterium).
- **D4 (fit (L) vs (Y)):** fit w przestrzeni E (nie log).
  (L): LSQ `E ~ C1 + C2*ln(d_meas)`, k = 2.
  (Y): `E ~ C1 - A*e^(-m*d_meas)`: skan siatkowy po m + LSQ (C1, A),
  k = 3. KRYTERIALNIE m in [0.3, 2.0]: hipoteza (Y) znaczy fizycznie
  MASYWNY Goldstone o masie skali substratu (skale: sqrt(V''(0)/kappa)
  = 1.0, sqrt(V''(s*)/kappa) = 1.33). Znana degeneracja zapisana
  Z GORY: dla m*d_max <~ 1 forma (Y) redukuje sie do (L) (granica
  bezmasowa: -A*e^(-m*d) ~ (A*m)*d - A, a przy malych m mimikruje
  ln(d)) — fit z m swobodnym raportowany WYLACZNIE jako czulosc, poza
  kryterium. AIC = n*ln(RSS/n) + 2k; Delta_AIC = AIC_Y - AIC_L;
  (L) wygrywa gdy Delta_AIC >= +10, (Y) gdy <= -10, inaczej BRAK.
- **D5 (predkosci):** zapisy co 10 krokow; roznice centralne z baza
  Delta_tau = 2.0 (jak D3 poprzedniego opa); v_zbl = -dd/dtau;
  test v*d ~ const: mediana i rozstep na oknie d in [4, d0-2];
  poza kryterium (sekcja 1.4).
- **D6 (C2 z drogi dynamicznej):** w runach dynamicznych rejestrowane
  (d(tau), H(tau)); zbocze dH/dln(d) na oknie d in [5, d0-3]
  raportowane jako NIEZALEZNE potwierdzenie C2 (poza kryterium G4;
  po drodze kwazistatycznej cala dyssypowana energia pochodzi
  z E_int).
- **D7 (geometria G5):** +1 w (-32,-4) i (-32,+4); -1 w (+32,-4)
  i (+32,+4): pary rownoimienne d_ss = 8, grupy przeciwnych znakow
  w odleglosci L/2 (sily miedzygrupowe znosza sie z symetrii obrazow).
  Okno tau = 200. Klasyfikacja: Delta d_ss >= +1.0 w OBU grupach ->
  ODPYCHANIE rownoimiennych (LADUNKOWY); <= -1.0 w obu ->
  przyciaganie (UNIWERSALNY); inaczej "nierozstrzygniete" wprost.
- **D8 (G6):** scan (+2,-2) na identycznej siatce d i identyczna
  procedura D3/D4; ratio = C2(+2,-2)/C2(+1,-1); baseline (+1,-1)
  policzony ta sama sciezka kodu w tym samym pliku Phase 3 (identyczna
  procedura, zero zaleznosci od artefaktow Phase 2). Rozszczepienie
  rdzenia n=2: osobny dlugi run obserwacyjny (poza kryterium).
- **D9 (kontrola dx=0.25):** L=64, N=256, dt=0.02 (stabilnosc:
  dt < dx^2/(4*kappa) = 0.03125); profil rdzenia + statyka
  d in {6, 8, 12}; raport C2_ctrl (3 punkty, poza kryterium G4,
  sluzy G3) + test pinningu: energia pary przy przesunieciu rdzeni
  o pol wezla (Delta_E_pin) raportowana dla dx=0.5 i dx=0.25;
  klasyfikacja "pinning istotny" gdy Delta_E_pin > 2% rozstepu E(d)
  na skanie albo rdzenie nieruchome w dynamice.
- **D10 (metrologia ogolna):** przejeta z op-lock-interaction:
  H_Gamma probkowane co 10 krokow, tolerancja monotonicznosci
  1e-10 * max(1, |H0|); S_GUARD = 10 (flaga, nie clamp); determinizm
  bitowy sprawdzany w Phase 1 (run x2); flaga boundary_contact
  stosowalna tylko w scenariuszu genesis (w scenariuszach wirowych
  |psi| ~ s* wypelnia pudlo z definicji — flaga bezprzedmiotowa).
- **D11 (genesis_u1):** na siatce opa (N=256, L=128):
  bare = szum U(-0.05, +0.05) * e^(i*0.7), seed 20260704, 6000 krokow,
  kryterium zaniku: metric_area(Phi > eps) < 4/N^2;
  cluster = cluster(A0=0.6, D=3, K=7) * e^(i*0.7) ewoluowany do
  R >= 8; lock: Phi_max in [0.8, 1.2] * Phi*. Kowariancja fazowa
  U(1): |psi| z runu theta0 = 0.7 vs identyczny run rzeczywisty
  (theta0 = 0) — max roznica raportowana (oczekiwana ~ eps maszynowe).
