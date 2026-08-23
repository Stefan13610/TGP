---
title: "Phase0_balance — DRAFT: B2-prime — ilosciowa refrakcja fal na defekcie (wirze) — most geometryczny wokol materii-defektow"
date: 2026-07-05
type: phase0-draft
tgp_owner: research/op-vortex-refraction-2026-07-05
status: LOCKED
anti_lakatos_lock: ACTIVE (autoryzacja autora: polecenie realizacji cyklu, 2026-07-05; lock przed kodem; wykonawca: osobny agent, nowa sesja)
tags: [gravity-bridge, branch-B2prime, vortex-lensing, effective-metric, refraction, eikonal, universality]
related:
  - "[[../op-stationary-background-2026-07-05/Phase_FINAL_close.md]]"
  - "[[../op-phi-metric-refraction-2026-07-04/Phase_FINAL_close.md]]"
  - "[[../op-goldstone-mediator-2026-07-04/Phase_FINAL_close.md]]"
---

# Phase 0 — DRAFT: B2-prime — refrakcja na defekcie (wirze)

## 0. Hipoteza robocza i miejsce w programie

Stan lancucha (wszystko zmierzone, nic do wiary):
- B1: materia = DEFEKTY fazowe w wygenerowanej przestrzeni Phi*;
  kanal fazowy dlugozasiegowy, ale ladunkowy.
- B2: refrakcja na duzym obiekcie niemierzalna (tachion sciany,
  gamma = 0.125 — wlasnosc SCIANY rosnacego obiektu).
- op-stationary-background: wokol wiru sektor falowy CZYSTY
  (lambda_full = +0.0008; gamma_rdzenie = 0.0003; nawet toy-sektor
  +0.0198); deficyt amplitudy defektu to soczewka SKUPIAJACA:
  n(r) > 1, ogon algebraiczny n^2 - 1 ~ 2.5690 / (r^2 (omega^2 -
  V''(s*))), SLEPA na znak kretu; demo pakietu: kierunek KU wirowi;
  tlo quasi-stacjonarne (dryf pary 0.015/tau, 25x wolniej od pakietu).

Ten cykl pyta ILOSCIOWO:

> (P1) POLICZALNOSC: czy kat ugiecia pakietu amplitudowego na wirze
>      zgadza sie z eikonalem z samego substratu (znak + wielkosc,
>      pasmo [0.5, 2.0])?
> (P2) UNIWERSALNOSC: czy ugiecie nie zalezy od ZNAKU kretu defektu
>      (test grawitacyjnosci: materia-defekt ugina fale tak samo
>      niezaleznie od "ladunku")?

G4+G5 PASS = pierwszy ilosciowo potwierdzony kanal o znaku
grawitacji: uniwersalnie przyciagajaca geometria propagacji wokol
materii-defektow.

## 1. Predykcje algebraiczne do zalockowania (PRZED kodem)

### 1.1 Rownanie #1 — wspolczynnik zalamania wokol defektu

```text
n^2(r) = (omega^2 - V''(f(r))) / (omega^2 - V''(s*))
f(r): profil wiru z ODE radialnego (kod wzorcowy: Phase0_check
  poprzednich opow); ogon: n^2 - 1 -> A_t / (r^2 (omega^2 - V''(s*))),
  A_t = V'''(s*) * s* * xi^2 = 2.5690
kanal amplitudowy propaguje dla omega > m_TV = 0.9374;
k0(omega) = sqrt((omega^2 - V''(s*)) / kappa);  v_g = kappa*k0/omega
```

### 1.2 Rownanie #2 — zaklad eikonalny (JUZ policzony, cytowany)

Z `op-stationary-background/Phase0_output.txt` (ray tracing na
profilu ODE; alpha > 0 = KU wirowi):

```text
omega=1.1:  b=6: +43.0   b=8: +14.0   b=12: +5.30   [deg]
omega=1.3:  b=6: +12.4   b=8: +5.24   b=12: +2.11
omega=1.0:  strefa przechwytu przy b <~ 6 (petle) — poza kryteriami
```

Phase0_check TEGO cyklu przelicza te tabele na rzeczywistym tle
(para, nie izolowany profil) + dodaje b=16; roznice vs tabela
cytowana raportowane (oczekiwane < kilka %, drugi wir w odleglosci
~90). UWAGA falowa zapisana z gory: dlugosc fali lambda = 2*pi/k0
= 7.7 (omega=1.1) jest porownywalna z b — rezim dyfrakcyjny;
dlatego pasmo G4 jest szerokie [0.5, 2.0], a punkty b=6 (silne pole)
i omega=1.0 (przechwyt) sa POZA kryterium (raport).

### 1.3 Rownanie #3 — dekompozycja parzysta/nieparzysta (P2)

```text
alpha_even = [alpha(b) + alpha(-b)] / 2      (czesc uniwersalna,
   soczewkowa — porownywana z eikonalem; n(r) nie zna znaku kretu)
alpha_odd  = [alpha(b) - alpha(-b)] / 2      (czesc cyrkulacyjna,
   typu Iordanskiego/Magnusa — sprzezona ze znakiem kretu;
   RAPORTOWANA, poza kryterium)
test lustrzany: zamiana znakow wirow tla (+1 <-> -1) przy tym samym
   b musi dac alpha_even identyczne w granicach niepewnosci
```

### 1.4 Dryf tla (korekta pre-rejestrowana)

Tlo = para (+1,-1) szachownicowa: zjazd zmierzony w
op-stationary-background: ~0.015/tau (przyspieszajacy; 1.21 na
tau=80 w drugiej polowie relaksu x4). Obowiazkowa korekta:
pozycja rdzenia mierzona z pola (detektor windingu) W OKNIE
przelotu pakietu; parametr zderzenia b_eff = najmniejsza odleglosc
centroidu pakietu od rdzenia (obie wielkosci z pomiaru, nie
z nominalu). Budzet okna: przy locie pakietu ~tau=200 rdzen
przesuwa sie o ~3 — b_eff czyni pomiar poprawnym z definicji.

## 2. Model i parametry (propozycja do LOCKa)

```text
substrat: jak B1 (kappa=0.5, a=0.5, b=1.6, c=1.0; psi in C;
  N=256, L=128, dx=0.5, periodycznie, eps=0.30, S_GUARD=10 flaga,
  seed=20260704)
tlo: para (+1 w (-32,-32), -1 w (+32,+32)); ansatz ZLOZONY
  (para pozioma D1 + para pionowa transponowana przez punkt
  posredni C — kod wzorcowy: op-stationary-background/
  Phase2_vortex_background.py, funkcje theta_h_pair/theta_v_pair;
  UWAGA: goly iloczyn sinh dla pary diagonalnej ma skok pi na
  szwie — NIE uzywac); relaksacja 2000 krokow (dt=0.02)
sektor falowy: PELNY zespolony II rzedu:
  d^2psi/dtau^2 = kappa*Lap(psi) - psi*(a - b|psi| + c|psi|^2)
  leapfrog pozycyjny dt=0.05; u_tau = roznica centralna (D1 z B2)
pakiet (kanal amplitudowy, w lokalnej ramce fazy tla):
  delta_psi = e^{i theta_bg} * u0 * env * cos(k0 (x - x0));
  predkosc poczatkowa analitycznie (jak D2 z B2);
  u0 = 1e-3, sx = 8, sy = 8; start x0 = -60, y0 = y_wir + b
skan glowny: omega in {1.1, 1.3} x b in {8, 12, 16}
  (b=6 i omega=1.0: obserwacje poza kryteriami)
uniwersalnosc (G5): (omega=1.1, b=8) w czterech wariantach:
  {wir +1, b > 0}, {wir +1, b < 0}, {tlo lustrzane -1/+1, b > 0},
  {tlo lustrzane, b < 0} -> alpha_even/alpha_odd wg #1.3
obserwable: pseudo-ped delta_psi w oknie transmisji (jak D3 B2):
  P = -Re< d(delta_psi)/dtau, grad(delta_psi) >, okno x > x_rdzenia+10;
  5 probek co Delta_tau=4 od przekroczenia okna przez centroid;
  usrednienie WEKTOROWE (Px, Py), potem kat; sigma z probek
kontrola dryfu rdzenia: gamma_core (okno r < 5) <= 0.01 (cytowany
  wynik: 0.0003); pozycje rdzeni zapisywane co 10 krokow
metrologia: determinizm bitowy w Phase 1; dryft sekularny energii
  <= 1e-4 (D4 z B2); wrap-guard (0.9 L)
```

## 3. Scenariusze

1. `background_check`: relaksacja tla; residuum, dryf, kret,
   diagnostyka szwu (max skok fazy poza rdzeniami << pi).
2. `dispersion_TV`: kanal amplitudowy na tle jednorodnym s*
   (fale plaskie, omega(k) vs kappa*k^2 + V''(s*)) — kalibracja G2.
3. `deflection(b, omega)`: skan glowny; alpha vs alpha_pred; b_eff.
4. `universality`: cztery warianty G5 (#1.3).
5. `linearity`: u0 -> u0/2 przy (omega=1.1, b=8).
6. obserwacje poza kryteriami: b=6 (silne pole), omega=1.0
   (przechwyt), alpha_odd (efekt cyrkulacyjny).

## 4. Kryteria (propozycja do LOCKa)

- **G1 techniczne**: determinizm; brak NaN/S_GUARD; dryft sekularny
  energii <= 1e-4; gamma_core <= 0.01 w kazdym runie; wrap-guard;
  kret tla zachowany. PASS/FAIL.
- **G2 kalibracja**: dyspersja kanalu amplitudowego na tle s*:
  blad <= 5% dla k in [0.5, 1.5]. FALSYFIKATOR: sektor zle
  zdyskretyzowany.
- **G3 istnienie refrakcji**: |alpha(b=8, omega=1.1)| > 3*sigma
  (probki czasowe) ORAZ |alpha| maleje z b (przy omega=1.1)
  i z omega (przy b=8).
- **G4 RDZEN — policzalnosc (P1)**: dla WSZYSTKICH szesciu punktow
  kryterialnych (b in {8,12,16} x omega in {1.1,1.3}): znak KU
  wirowi ORAZ alpha/alpha_pred(b_eff) in [0.5, 2.0].
  FALSYFIKATOR: geometria propagacji nieczytelna z substratu
  w rezimie falowym.
- **G5 RDZEN — uniwersalnosc (P2)**: |alpha_even(wir +1) -
  alpha_even(tlo lustrzane)| <= 0.2 * |srednia| przy (b=8,
  omega=1.1); alpha_odd raportowane osobno (poza kryterium).
  FALSYFIKATOR: ugiecie zalezne od znaku kretu -> kanal
  geometryczny tez "ladunkowy".
- **G6 liniowosc**: alpha(u0) vs alpha(u0/2): roznica <= 10%.

## 5. Regula decyzyjna (propozycja)

```text
G4 PASS + G5 PASS
  -> MOST GEOMETRYCZNY POTWIERDZONY ILOSCIOWO: uniwersalnie
     przyciagajaca, policzalna geometria propagacji wokol
     materii-defektow. Nastepne opy (osobne locki): (i) skalowanie
     z "masa" defektu (|n| = 2: A_t ~ n^2 -> ugiecie x4?),
     (ii) uklady wielu defektow (superpozycja soczewek),
     (iii) granica newtonowska 2D i G_eff z parametrow substratu.
G4 PASS + G5 FAIL
  -> geometria policzalna, ale znakowo-zalezna (dominacja czesci
     cyrkulacyjnej) — niespodzianka; analiza rozdzialu soczewka/
     cyrkulacja jako osobny op.
G4 FAIL (wielkosc poza pasmem przy dobrym znaku)
  -> eikonal nie opisuje rezimu falowego (lambda ~ b, dyfrakcja):
     raport ilosciowy rozjazdu; powtorka w rezimie krotkofalowym
     (wyzsze omega / wieksze b / wieksze L) jako osobny op.
G4 FAIL (zly znak) -> falsyfikacja kanalu geometrycznego wokol
     defektow; program wraca do analizy.
G2/G3 FAIL -> stop przed interpretacja.
```

## 6. Forbidden moves (propozycja)

1. Zero pinningu; kret czytany z pola; korekta dryfu WYLACZNIE
   przez b_eff i pozycje mierzone (procedura #1.4), zadnego
   "przypinania" tla.
2. Zadnych zmian progow/parametrow/procedur po pierwszym runie
   (wyjatki: jedno obnizenie dt; jeden kontrolny run dx=0.25).
3. Punkty poza-kryterialne (b=6, omega=1.0, alpha_odd) nie moga
   byc dolaczone do G4 po fakcie; kryterialny zestaw jest zamkniety.
4. Nie interpretowac alpha_odd jako grawitacji; nie deklarowac
   GR/soczewkowania obserwacyjnego; tau != czas.
5. Wynik "G5 FAIL" raportowac jako pelnoprawne rozstrzygniecie.
6. Nie promowac do core/ bez osobnej autoryzacji.

## 7. Czego ten cykl NIE robi

- Nie mierzy skalowania z |n| ani ukladow wielu defektow (nastepne).
- Nie rozstrzyga kanalu dylatacyjnego (galaz C, osobna procedura).
- Nie wyprowadza metryki g_mu_nu; zero deklaracji GR.

## 8. Deliverables

- `Phase0_balance.md` (po autoryzacji: LOCKED + doprecyzowania
  pre-code w sekcji 10 PRZED kodem)
- `Phase0_check.py`/`Phase0_output.txt` (n^2(r) na tle pary,
  alpha_pred(b_grid, omega_grid) + porownanie z tabela cytowana,
  k0/v_g, budzet dryfu)
- `Phase1_calibration.py`/`Phase1_output.txt` (G2 + background_check
  + czesc G1)
- `Phase2_deflection.py`/`Phase2_output.txt` (G3, G4, G6 + obserwacje)
- `Phase3_universality.py`/`Phase3_output.txt` (G5 + alpha_odd)
- `Phase_FINAL_close.md`, `README.md`

## 9. Status

**LOCKED (2026-07-05).** Model, predykcje (w tym cytowany zaklad
eikonalny i jawne zastrzezenie dyfrakcyjne), kryteria i regula
decyzyjna spisane przed kodem. Autoryzacja autora: polecenie
realizacji cyklu (2026-07-05). Wykonawca: osobny agent, nowa
sesja (spelnione). Obowiazuje pelny anty-Lakatos. Doprecyzowania
pre-code: sekcja 10 (spisane PRZED pierwsza linijka kodu).

## 10. Doprecyzowania pre-code (spisane przed pierwsza linijka kodu)

Zadne z ponizszych nie zmienia progow G1-G6, modelu ani predykcji —
to ujednoznacznienie procedur oraz jawna pre-rejestracja dwoch
obserwacji strukturalnych (D11, D12) wykrytych rachunkiem/analiza
przed kodem:

- **D1 (integrator):** pelny sektor zespolony II rzedu, leapfrog
  pozycyjny: `psi_{n+1} = 2 psi_n - psi_{n-1} + dt^2 * RHS(psi_n)`,
  `RHS(psi) = kappa*Lap(psi) - psi*(a - b|psi| + c|psi|^2)`;
  dt=0.05 << 2/omega_max ~ 0.49 (omega_max^2 = 8 kappa/dx^2
  + max V''). `psi_tau` do obserwabli = roznica CENTRALNA
  `(psi_{n+1} - psi_{n-1})/(2 dt)` (errata #1 B2).
- **D2 (pakiet, start analityczny):** `psi(0) = psi_bg + e^{i
  theta_bg} * u0 * env * cos(k0 (x - x0))`, `psi(-dt) = psi_bg +
  e^{i theta_bg} * u0 * env * cos(k0 (x - x0) + omega dt)` (pakiet
  prawobiezny — wzorzec: demo G5 op-stationary); env = exp(-(x-x0)^2
  /(2 sx^2) - (y-y0)^2/(2 sy^2)), sx = sy = 8, u0 = 1e-3,
  x0 = -60, y0 = y_wir(-32) + b; k0 = sqrt((omega^2 - V''(s*))
  / kappa).
- **D3 (referencja blizniacza — formalizacja korekty dryfu #1.4):**
  `delta_psi(tau) = psi(tau) - psi_ref(tau)`, gdzie psi_ref
  ewoluuje TYM SAMYM integratorem z psi_bg bez pakietu (start:
  psi_ref(-dt) = psi_ref(0) = psi_bg), w jednej petli (lockstep).
  To odejmuje wlasny ruch tla (zjazd pary + relaksacja artefaktu C)
  z obserwabli dokladnie w rzedzie liniowym; pozycje rdzeni do b_eff
  i okna czytane z psi_ref (detektor windingu + centroid deficytu,
  D2 z B1). To NIE jest pinning (zero ingerencji w dynamike).
- **D4 (obserwabla i chwila pomiaru):** pseudo-ped `P = -Re
  sum(conj(delta_psi_tau) * grad(delta_psi)) dx^2` w oknie
  transmisji `x > x_rdzen(+1 przy (-32,-32); z psi_ref) + 10`
  (prog x zamrozony w chwili wyzwolenia); wyzwolenie: frakcja
  |delta_psi|^2 w oknie > 0.5; potem 5 probek co Delta_tau = 4;
  usrednienie WEKTOROWE (Px, Py) -> alpha; sigma_alpha = odch. std
  probek katowych. Konwencja znaku: `alpha := -sign(b) *
  atan2(Py, Px)` [deg] — alpha > 0 = KU wirowi dla obu znakow b.
- **D5 (b_eff):** b_eff = min po tau (probki co 10 krokow, do chwili
  wyzwolenia pomiaru) odleglosci centroidu |delta_psi|^2 (globalnego)
  od rdzenia rozpraszajacego (pozycja z psi_ref, z uwzglednieniem
  okresowosci). alpha_pred(b_eff) = ray tracing PRZY b = b_eff na tym
  samym zrelaksowanym tle co run (bez interpolacji tabeli).
- **D6 (eikonal na rzeczywistym tle):** n^2(x,y) = (omega^2 -
  V''(|psi_bg|)) / (omega^2 - V''(s*)) z bilinearna interpolacja;
  ray tracing hamiltonowski (D7 z B2): dr/dsigma = p, dp/dsigma =
  grad(n^2)/2, |p(start)| = n(start), RK4, dsigma = 0.02; start
  x = -60 (28 przed rdzeniem), y = -32 + b, kierunek +x; wyjscie
  x > +40 (72 za rdzeniem; drugi wir w odleglosci >= 58, wklad
  n^2-1 < 0.3%). Phase0_check raportuje roznice wzgledem tabeli
  cytowanej (izolowany profil ODE, start/meta -60/+60 wzgledem
  rdzenia) — oczekiwane kilka % (krotszy odcinek naplywu).
- **D7 (G2, dyspersja kanalu amplitudowego):** tlo jednorodne
  psi = s*; zaburzenie RZECZYWISTE u0*cos(k x), u_tau(0) = 0 (pelna
  dynamika II rzedu; pole pozostaje rzeczywiste — czysty kanal
  amplitudowy); k = 2 pi m / L, m in {11, 15, 20, 25, 30}
  (k in [0.54, 1.47] — wewnatrz pasma kryterialnego [0.5, 1.5]);
  omega_meas ze zliczania przejsc przez zero Re(psi - s*) w punkcie
  probki, okno tau = 100; omega_teor = sqrt(kappa k^2 + V''(s*));
  dyspersja sieciowa raportowana obok (poza kryterium).
- **D8 (G1, metryki):** dryft sekularny energii = |<H> ostatnia
  cwiartka - <H> pierwsza cwiartka| / |H(0)| (H = kinetyczna +
  gradientowa + potencjalna pelnego psi), prog 1e-4; oscylacja
  leapfroga O((dt omega)^2) raportowana osobno. gamma_core = slope
  fitu ln RMS |delta_psi| w oknach r < 5 wokol OBU rdzeni (pozycje
  z psi_ref), fit na DRUGIEJ POLOWIE runu (po tranzycie pakietu;
  przejsciowy wzrost przy przelocie jest oczekiwany i nie jest
  niestabilnoscia), prog 0.01. Kret tla: n_tot(psi_ref, koniec) = 0
  i obecne oba rdzenie. Wrap-guard: |przemieszczenie centroidu od
  startu| < 0.9 L. S_GUARD = 10 flaga (nie clamp). Determinizm
  bitowy: powtorka pelnej relaksacji tla + jednego runu falowego
  w Phase 1. Zapisy co 10 krokow; seed = 20260704.
- **D9 (budzet krokow):** relaksacja tla: 2000 krokow dt=0.02
  (przeplyw gradientowy, jak op-stationary). Runy falowe: limit
  6000 krokow (tau = 300) dla omega = 1.1/1.3; 8000 (tau = 400)
  dla obserwacji omega = 1.0 (v_g = 0.246); run konczy sie po
  5 probkach albo na guardzie/limicie (brak 5 probek = raport,
  run niewazny do kryterium, przyczyna zapisana).
- **D10 (G6, liniowosc):** |alpha(u0) - alpha(u0/2)| <= 0.1 *
  max(|alpha(u0)|, |alpha(u0/2)|) przy (omega = 1.1, b = 8).
- **D11 (PRE-REJESTRACJA: symetria lustrzana):** tlo lustrzane =
  negacja fazy ansatzu (theta -> -theta; amplituda identyczna);
  poniewaz wszystkie wspolczynniki rownan sa rzeczywiste,
  a pakiet wstrzykiwany jest w lokalnej ramce fazy SWOJEGO tla,
  caly run lustrzany jest DOKLADNIE sprzezeniem zespolonym runu
  zwyklego -> oczekiwanie: alpha(lustro, b) = alpha(wir +1, b)
  z dokladnoscia zaokraglen. G5 w brzmieniu LOCKa testuje wiec
  (a) te symetrie implementacyjnie, (b) fizyczna uniwersalnosc
  niesie czesc parzysta alpha_even po znaku b (roznica alpha(+b)
  vs alpha(-b) = czesc cyrkulacyjna alpha_odd, raportowana).
  Oczekiwanie zapisane PRZED kodem; dane rozstrzygaja.
- **D12 (uwaga dyfrakcyjna, powtorzona):** lambda(1.1) = 7.72,
  lambda(1.3) = 4.93 vs b in {8,12,16} — rezim falowy; pasmo G4
  [0.5, 2.0] jest szerokie z tego powodu (sekcja 1.2); zadnych
  poprawek dyfrakcyjnych po fakcie.
- **D13 (narzedzia):** zero `cd` w powloce; sciezki pelne od
  korzenia vaulta; kazdy skrypt pisze output do
  `PhaseX_output.txt` w katalogu opa.
