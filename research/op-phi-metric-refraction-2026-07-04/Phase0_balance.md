---
title: "Phase0_balance — DRAFT: Phi jako metryka efektywna propagacji (refrakcja) — most do grawitacji, galaz B2"
date: 2026-07-04
type: phase0-draft
tgp_owner: research/op-phi-metric-refraction-2026-07-04
status: LOCKED
anti_lakatos_lock: LOCKED (2026-07-04, autoryzacja autora — polecenie kontynuacji programu po zamknieciu B1; lock przed pierwsza linijka kodu)
tags: [gravity-bridge, branch-B2, effective-metric, refraction, eikonal, wave-sector]
related:
  - "[[../op-lock-interaction-gravity-2026-07-04/Phase_FINAL_close.md]]"
  - "[[../../core/sek01_ontologia/sek01_ontologia.tex]]"
  - "[[../../core/sek08_formalizm/sek08_formalizm.tex]]"
---

# Phase 0 — DRAFT: refrakcja malych zaburzen na tle Phi(x) (galaz B2)

## 0. Hipoteza robocza i uczciwe postawienie problemu

Galaz B2 z reguly decyzyjnej `op-lock-interaction-gravity`: "grawitacja"
= refrakcja trajektorii malych zaburzen w gradiencie `Phi`, nie sila
miedzy cialami. Rdzen (sek01, aksjomat c(Phi); sek08c metryka z
substratu) POSTULUJE, ze predkosc propagacji zalezy od `Phi` tak, iz
regiony gestsze uginaja promienie KU sobie (przyciaganie, soczewka).

Ten cykl pyta o wersje minimalna i falsyfikowalna:

> Czy MINIMALNY substrat (ten sam potencjal, bez zadnych nowych
> sprzezen), uzupelniony wylacznie o sektor falowy (druga pochodna
> po tau), daje refrakcje zaburzen na tle zlockowanego obiektu —
> w ktora STRONE ugina, i czy kat ugiecia jest policzalny (eikonal)?

DWA jawnie oddzielone pytania (zapisac wprost przy locku):
(P1) POLICZALNOSC: czy ugiecie zgadza sie z rachunkiem eikonalnym
     z samego substratu (znak + wielkosc)?
(P2) GRAWITACYJNOSC: czy znak jest przyciagajacy (ku obiektowi)?

## 1. Predykcje algebraiczne do zalockowania (PRZED kodem)

### 1.1 Rownanie #1 — sektor falowy i dyspersja

Minimalne rozszerzenie (JEDYNA nowa struktura tego cyklu, zapisana
wprost): dynamika falowa malych zaburzen `u` na ZAMROZONYM tle `s_bg`:

```text
d^2u/dtau^2 = kappa*Lap(u) - V''(s_bg(x))*u      (linearyzacja)
dyspersja jednorodna: omega^2 = kappa*k^2 + V''(s_bg)
V''(0)   = a                — proznia falszywa: m_FV = sqrt(0.5) = 0.707
V''(s*)  = 0.87869          — proznia prawdziwa: m_TV = 0.937
```

Tlo zamrozone = uczciwa realizacja "propagacji NA tle Phi(x)" bez
mieszania z ruchem frontow (fronty rosna w przeplywie gradientowym —
osobny sektor, juz zmierzony).

### 1.2 Rownanie #2 — wspolczynnik zalamania i ZNAK ugiecia (RDZEN)

Eikonal dla `omega` propagujacej w obu fazach (`omega > m_TV`):

```text
n^2(x) = k^2(x)/k_ref^2 = (omega^2 - V''(s_bg(x))) / (omega^2 - V''(0))
wewnatrz obiektu: V''(s*) > V''(0)  =>  n < 1  (osrodek "rzadszy")
```

OCZEKIWANIE TEORETYCZNE (jawne, pre-code — analog uczciwej Yukawy
z poprzedniego opa): promienie uginaja sie OD obiektu (refrakcja
ODPYCHAJACA). Czyli: minimalny substrat ma policzalna geometrie
propagacji (P1: TAK), ale o znaku antygrawitacyjnym (P2: NIE).
Jesli tak wyjdzie — to PELNOPRAWNY wynik: "Phi jako metryka" wymaga
sprzezenia typu c(Phi)/kappa(Phi) (aksjomat rdzenia) jako DODATKOWEJ
struktury, ktorej minimalny potencjal nie generuje. Wynik przeciwny
(ugiecie KU obiektowi) bylby niespodzianka i mocnym kandydatem na most.

Przyklad liczbowy (do weryfikacji w Phase0_check): omega = 1.1:
`n_in = sqrt(0.331/0.710) = 0.68` — kontrast silny, kat mierzalny.

### 1.3 Rownanie #3 — kat ugiecia (eikonal, do porownania z pomiarem)

```text
alpha_pred(b, omega) = calka eikonalna po promieniu (ray tracing
numeryczny w Phase0_check na tym samym tle s_bg, ktore dostanie
pakiet) — zadnych wzorow dopasowywanych po fakcie
```

### 1.4 Zakres czestosci

```text
omega < m_FV:          brak propagacji (zaniedbac)
m_FV < omega < m_TV:   obiekt = bariera (calkowite odbicie/tunelowanie)
                       — scenariusz kontrolny, obserwacja
omega > m_TV:          refrakcja — glowny pomiar; skan omega
                       in {1.0, 1.1, 1.3} (kontrast n maleje z omega
                       => predykcja: |alpha| maleje z omega)
```

## 2. Model i parametry (propozycja do LOCKa)

```text
tlo: s_bg = zlockowany obiekt z przeplywu gradientowego (cluster
  A0=0.6, D=3, K=7 -> R_obj = 12), ZAMROZONE (ciaglosc: identyczna
  preparatyka jak poprzednie opy)
sektor falowy: u(x, tau), rownanie #1, integrator leapfrog
  (symplektyczny); dt = 0.05 (Courant: dt < dx/sqrt(kappa) = 0.707)
grid: N=256, L=128 (dx=0.5), periodyczne; UWAGA: pakiet nie moze
  okrazyc pudla w oknie pomiaru — okno tau ograniczone z gory
pakiet: gaussowski, nosna k0 = sqrt((omega^2 - a)/kappa) w prozni
  falszywej, szerokosci sigma_x, sigma_y (do zalockowania, np. 8 i 12),
  start w odleglosci ~40 od obiektu, parametry zderzenia
  b in {8, 12, 16, 20} (b > R_obj czesciowo, skan przekrywa brzeg)
obserwable ugiecia: sredni pęd poprzeczny pakietu
  <p_y> = int( u_tau * du/dy ) (pseudo-pęd falowy) oraz przesuniecie
  centroidu |u|^2 po przejsciu; kat alpha = atan(delta_p_y / p_x)
energia falowa: H_u = int[ 0.5*u_tau^2 + 0.5*kappa*|grad u|^2
  + 0.5*V''(s_bg)*u^2 ] — ZACHOWANA (sektor konserwatywny!)
seed = 20260704; amplituda pakietu u0 << s_bar (liniowosc; test
  zbieznosci u0 -> u0/2 zapisany z gory)
metrologia D1-D10 z op-lock-interaction przejac tam, gdzie stosowalne
```

## 3. Scenariusze

1. `dispersion_check`: jednorodne tla (FV i TV): pomiar omega(k)
   vs rownanie #1 (kalibracja sektora falowego).
2. `barrier_regime`: omega in (m_FV, m_TV): obiekt jako bariera —
   obserwacja (odbicie), kontrola spojnosci.
3. `deflection(b, omega)`: przelot pakietu, pomiar alpha(b, omega);
   porownanie ZNAKU i WIELKOSCI z alpha_pred (rownanie #3).
4. `linearity`: u0 -> u0/2: alpha niezmienne (liniowosc potwierdzona).
5. `frozen_vs_live` (obserwacja, poza kryteriami): jeden run z tlem
   zywym (pelna nieliniowa ewolucja s+u drugiego rzedu) — czy wniosek
   o znaku przezywa odmrozenie tla.

## 4. Kryteria (propozycja do LOCKa)

- **G1 techniczne**: determinizm; brak NaN; zachowanie H_u w sektorze
  falowym (dryft wzgledny <= 1e-4 na okno pomiaru); pakiet nie okraza
  pudla w oknie. PASS/FAIL.
- **G2 kalibracja dyspersji**: omega(k) na obu tlach zgodne z
  rownaniem #1 z bledem <= 5% w zakresie pomiarowym. FALSYFIKATOR:
  sektor falowy zle zdyskretyzowany -> pomiar ugiecia niewykonalny.
- **G3 istnienie refrakcji**: |alpha(b=12, omega=1.1)| przekracza
  3x niepewnosc pomiarowa (bootstrap po oknach czasowych) i alpha
  maleje z b oraz z omega (jakosciowo, zgodnie z #1.4).
- **G4 RDZEN — policzalnosc (P1)**: znak alpha zgodny z eikonalem
  ORAZ `alpha/alpha_pred in [0.5, 2.0]` dla wszystkich (b, omega)
  z G3. FALSYFIKATOR: geometria propagacji nieczytelna z substratu.
- **G5 klasyfikator kierunku (P2, bez PASS/FAIL)**:
  `OD obiektu` (oczekiwane) -> minimalny substrat = refrakcja
  antygrawitacyjna; `KU obiektowi` -> gravity-like (niespodzianka);
  raport wprost.
- **G6 liniowosc**: alpha(u0) vs alpha(u0/2): roznica <= 10%.

## 5. Regula decyzyjna (propozycja)

```text
G4 PASS + G5 "OD obiektu" (oczekiwane)
  -> P1 TAK, P2 NIE: propagacja na tle Phi jest policzalna, ale
     minimalny substrat ugina OD gestosci. Grawitacja przez geometrie
     WYMAGA sprzezenia kappa(Phi)/c(Phi) (tak, jak postuluje aksjomat
     rdzenia) — to jest zmierzone wskazanie, CO trzeba dodac.
     Nastepny op (za autoryzacja): substrat z kappa(Phi) — minimalna
     realizacja aksjomatu c(Phi) i test, czy wtedy alpha zmienia znak
     i odtwarza soczewke; rownolegle B1 (sektor fazowy) niezalezny.
G4 PASS + G5 "KU obiektowi"
  -> mocny kandydat na most geometryczny: nastepny op — soczewkowanie
     ilosciowe, zaleznosc od "masy" (R_obj), granica newtonowska 2D.
G4 FAIL -> geometria propagacji nieczytelna: raport; sprawdzic
     dyskretyzacje (jeden kontrolny run dx=0.25 — wyjatek techniczny);
     jesli utrzymane -> B2 w tej postaci odlozone, priorytet B1.
G2 FAIL -> stop przed interpretacja.
```

## 6. Forbidden moves (propozycja)

1. Tlo zamrozone w scenariuszach kryterialnych; "zywe" tlo tylko jako
   obserwacja 5 (bez wplywu na G1-G6).
2. Nie zmieniac progow/parametrow/procedur po pierwszym runie
   (wyjatki: jedno obnizenie dt; jeden kontrolny run dx=0.25).
3. Nie dodawac zadnych sprzezen kappa(Phi)/c(Phi) w TYM cyklu —
   to jest dokladnie hipoteza nastepnego opa, nie tego.
4. Nie interpretowac pseudo-pedu falowego jako pedu czastki; zadnych
   deklaracji o geodezyjnych/GR/soczewkowaniu obserwacyjnym.
5. Wynik "OD obiektu" raportowac jako pelnoprawne rozstrzygniecie
   (wskazanie brakujacej struktury), nie jako porazke.
6. Nie promowac do core/ bez osobnej autoryzacji.

## 7. Czego ten cykl NIE robi

- Nie wyprowadza metryki g_mu_nu ani rownan Einsteina (sek08c robi to
  formalnie; tu tylko toy-test kierunku refrakcji).
- Nie testuje kappa(Phi) — dopiero nastepny op, jesli G5 = "OD".
- Nie liczy trajektorii czastek (brak czastek w modelu) — wylacznie
  pakiety falowe.

## 8. Deliverables

- `Phase0_balance.md` (po autoryzacji: LOCKED + ew. poprawki pre-code)
- `Phase0_check.py`/`Phase0_output.txt` (m_FV, m_TV, n(omega),
  eikonal/ray-tracing alpha_pred(b, omega) na tle modelowym)
- `Phase1_dispersion.py` / `Phase1_output.txt` (G2 + G1)
- `Phase2_deflection.py` / `Phase2_output.txt` (G3, G4, G5, G6)
- `Phase_FINAL_close.md`, `README.md`

## 9. Status

**LOCKED (2026-07-04).** Model (w tym jawnie zadeklarowane MINIMALNE
rozszerzenie o sektor falowy), predykcje (w tym uczciwe oczekiwanie
refrakcji ODPYCHAJACEJ!), kryteria i regula decyzyjna spisane przed
kodem. Autoryzacja autora: polecenie kontynuacji programu po
zamknieciu B1 (ktorego regula decyzyjna wskazuje B2 jako nastepny
pomiar). Obowiazuje pelny anty-Lakatos. Doprecyzowania pre-code:
sekcja 10 (spisane PRZED pierwsza linijka kodu).

## 10. Doprecyzowania pre-code (spisane przed pierwsza linijka kodu)

Zadne z ponizszych nie zmienia progow G1-G6, modelu ani predykcji —
to ujednoznacznienie procedur ORAZ jawna pre-rejestracja jednego
ryzyka fizycznego wykrytego rachunkiem przed kodem (D6):

- **D1 (integrator):** leapfrog w formie pozycyjnej:
  `u_{n+1} = 2u_n - u_{n-1} + dt^2 (kappa*Lap(u_n) - V''(s_bg)*u_n)`;
  start: `u_0` = pakiet, `u_{-1}` z analitycznego `u_tau(0)`;
  `u_tau` do obserwabli = roznica centralna `(u_{n+1}-u_{n-1})/(2dt)`.
  Stabilnosc: dt=0.05 << 2/omega_max ~ 0.49.
- **D2 (pakiet):** `u(x,y,0) = u0 * exp(-(x-x0)^2/(2*sx^2)
  - (y-b)^2/(2*sy^2)) * cos(k0*(x-x0))`,
  `u_tau(x,y,0) = +omega * u0 * exp(...) * sin(k0*(x-x0))`
  (pakiet prawobiezny); sx=8, sy=12 (wartosci z sekcji 2);
  u0 = 1e-3; x0 = -(R_obj + 40) zaokraglone do wezla; y0 = b.
- **D3 (obserwable i chwila pomiaru):** pseudo-ped
  `P = -int(u_tau * grad u) dx^2`; KRYTERIALNIE P liczone w
  polprzestrzeni transmisji `x > R_obj + 5` (wyklucza okolice
  sciany; por. D6), P globalne raportowane obok.
  `alpha = atan(P_y / P_x)` (P_y(0)=0 z symetrii pakietu, wiec
  Delta P_y = P_y). Pomiar: 5 probek co Delta_tau=4, poczynajac od
  chwili, gdy centroid |u|^2 w oknie transmisji przekroczy
  x = R_obj + 20; alpha = srednia, sigma_alpha = odchylenie std
  probek (bootstrap po oknach czasowych z G3). Centroid |u|^2
  w oknie transmisji raportowany jako obserwabla wtorna.
- **D4 (G1, dryft H_u):** leapfrog daje OGRANICZONA oscylacje energii
  O((dt*omega)^2) ~ 1e-3 bez dryfu sekularnego; "dryft" z G1 =
  |<H_u> ostatnia cwiartka okna - <H_u> pierwsza cwiartka| / H_u(0)
  (trend sekularny), prog 1e-4; amplituda oscylacji raportowana
  osobno (poza kryterium).
- **D5 (G2, dyspersja):** fale plaskie `u = cos(k*x)`, `u_tau = 0`,
  k = 2*pi*m/L, m in {10, 15, 20, 25, 30}; omega z zliczania przejsc
  przez zero probki u(punkt) na oknie tau=100; tla: FV (s_bg=0)
  i TV (s_bg=s*); prog 5% vs rownanie #1; dyspersja sieciowa
  (kappa*(2-2cos(k*dx))/dx^2 zamiast kappa*k^2) raportowana.
- **D6 (PRE-REJESTRACJA: tachion tla zamrozonego):** rachunek przed
  kodem: `V''(s) = a - 2bs + 3cs^2 < 0` dla `s in (0.190, 0.877)`
  (minimum -0.353 przy s = b/(3c) = 0.533). Sciana zamrozonego
  obiektu przechodzi przez to pasmo, a studnia przyciagajaca w 1D/2D
  ZAWSZE wiaze stan -> oczekiwane `lambda_min(L) < 0` dla
  `L = -kappa*Lap + V''(s_bg)` -> tlo zamrozone jest TACHIONOWE
  w scianie (fizycznie: zamrozony mod wzrostu frontu; w zywym
  przeplywie to po prostu rosniecie obiektu, nie niestabilnosc).
  Procedura zapisana Z GORY:
  (a) Phase0_check liczy lambda_min (odwrotna/przesunieta iteracja
      potegowa) i gamma = sqrt(-lambda_min), oraz lokalizacje modu;
  (b) scenariusz `tachyon_check` (tlo zamrozone + szum 1e-8, seed
      20260704, BEZ pakietu): pomiar tempa wzrostu gamma_meas vs
      gamma — empiryczne potwierdzenie;
  (c) kryterium waznosci pomiaru ugiecia (KONTAMINACJA): jesli
      w runie deflection ||u|| w oknie sciany (r < R_obj+5) rosnie
      wykladniczo (fit gamma > 0.02) I przekracza 10% ||u|| w oknie
      transmisji w chwili pomiaru -> pomiar alpha NIEWAZNY;
  (d) jesli wszystkie runy deflection niewazne wg (c): G3/G4/G6
      NIEAPLIKOWALNE (nie FAIL) — werdykt cyklu: "geometria
      propagacji na TLE ZAMROZONYM niemierzalna, bo zamrozenie
      generuje tachion sciany; B2 wymaga tla stacjonarnego
      (op-wall-dynamics) albo sektora bez sciany" + raport
      eikonalny (strona teoretyczna P1). To NIE jest zmiana progow:
      to pre-rejestrowany warunek stosowalnosci pomiaru.
  UWAGA: rownanie liniowe => kontaminacji NIE wykryje G6 (wszystko
  skaluje sie z u0) — dlatego kryterium (c) jest strukturalne.
- **D7 (eikonal, rownanie #3):** ray tracing w formie hamiltonowskiej
  `dr/dsigma = p`, `dp/dsigma = grad(n^2)/2`, `|p(start)| = n(start)`,
  RK4, dsigma=0.05, n^2 z bilinearnej interpolacji na siatce tla;
  start x=-60, kierunek +x; alpha_pred = kat p na wyjsciu (x>+60).
  Ta sama funkcja tla co pakiet. UWAGA zapisana z gory: w scianie
  V''<0 => n^2 > 1 (pierscien skupiajacy) — znak wypadkowy moze
  zalezec od b; predykcja kierunku = wynik ray tracingu, nie
  uproszczony argument "n<1 w srodku" z sekcji 1.2.
- **D8 (wrap guard):** okno pomiaru konczy sie zanim centroid pakietu
  przebedzie 0.9*L od startu; flaga wrap raportowana.
- **D9 (tlo):** preparatyka identyczna jak lancuch: cluster(A0=0.6,
  D=3, K=7), przeplyw gradientowy do R_obj >= 12 (half-width
  Phi > eps=0.30 wzdluz osi x, interpolacja subgridowa), ZAMROZONE;
  determinizm bitowy preparatyki i sektora falowego w Phase 1.
- **D10 (metrologia):** zapisy co 10 krokow; S_GUARD=10 flaga
  (na u — sygnalizuje tachion, nie clamp); seed=20260704;
  barrier_regime: omega=0.8 (obserwacja, poza kryteriami);
  skan glowny omega in {1.0, 1.1, 1.3}, b in {8, 12, 16, 20}.
- **D11 (frozen_vs_live, scenariusz 5):** "tlo zywe" = s_bg
  ewoluuje rownolegle przeplywem gradientowym (ds/dtau = kappa*Lap s
  - V'(s)), V''(s_bg(tau)) w rownaniu falowym aktualizowane co krok;
  wylacznie obserwacja (bez wplywu na G1-G6).
