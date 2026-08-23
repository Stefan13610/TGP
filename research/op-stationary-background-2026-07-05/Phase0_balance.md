---
title: "Phase0_balance — tlo stacjonarne dla sektora falowego: babel krytyczny vs wir (odblokowanie B2)"
date: 2026-07-05
type: phase0-balance
tgp_owner: research/op-stationary-background-2026-07-05
status: LOCKED
anti_lakatos_lock: LOCKED (2026-07-05, autoryzacja autora — wybor opcji 2 "stabilizacja obiektu" po zamknieciu galezi C; lock przed pierwsza linijka kodu)
tags: [gravity-bridge, stationary-background, critical-bubble, vortex, wave-sector, tachyon-resolution]
related:
  - "[[../op-phi-metric-refraction-2026-07-04/Phase_FINAL_close.md]]"
  - "[[../op-goldstone-mediator-2026-07-04/Phase_FINAL_close.md]]"
  - "[[../op-lock-interaction-gravity-2026-07-04/Phase_FINAL_close.md]]"
---

# Phase 0 — tlo stacjonarne dla sektora falowego

## 0. Hipoteza robocza

B2 padlo, bo zamrozone tlo rosnacego obiektu ma tachion sciany
(gamma = 0.125, zmierzone). Pytanie tego cyklu:

> Czy substrat (bez zadnych nowych sprzezen) POSIADA tlo stacjonarne
> i falowo stabilne, na ktorym pomiar propagacji (B2-prime) jest
> wazny? Dwaj kandydaci z lancucha: (K1) babel krytyczny R_c
> (stacjonarny przez balans naciag/cisnienie) i (K2) wir
> (stacjonarny przez TOPOLOGIE — zgodnie z reinterpretacja B1:
> "defekty zyja W wygenerowanej przestrzeni").

Hipoteza: K1 jest siodlem (mod wzrostu = ten sam tachion co B2),
K2 jest stabilny w PELNEJ linearyzacji zespolonej — a tachion B2 byl
artefaktem zlozenia: (a) niestacjonarnego tla, (b) toy-sektora
RZECZYWISTEGO, ktory gubi czlony fazowe stabilizujace wir.

## 1. Predykcje algebraiczne (PRZED kodem)

### 1.1 Rownanie #1 — babel krytyczny (K1)

Z zalockowanego prawa frontu B (op-lock-interaction, G2 PASS):
`v(R) = v_inf - c/R`, `v_inf = 0.15993`, `c = 0.4898` (CYTOWANE,
zero strojenia):

```text
R_c_pred = c / v_inf = 3.0626    (pasmo testu: [0.8, 1.2]*R_c_pred)
lambda_min(Hessian na bablu) < 0  (mod wzrostu; sciana babla = ta
sama sciana, ktorej tachion zmierzylo B2: gamma ~ 0.12)
=> K1 stacjonarny, ale NIE stabilny -> odrzucony jako tlo
```

### 1.2 Rownanie #2 — dwa operatory na tle wirowym (RDZEN)

Tlo: para wirow (+1,-1) w konfiguracji szachownicowej
(p i p + (L/2, L/2)) — sily kasuja sie z symetrii obrazow.

```text
L_toy  = -kappa*Lap + V''(|psi_bg|)     (sektor B2: rzeczywiste u)
L_full = Hessian H wokol psi_bg (pelna linearyzacja zespolona):
  A(u) = -kappa*Lap(u) + e^{i theta}[ V''(f)*alpha + i*(V'(f)/f)*beta ]
  f = |psi_bg|, theta = arg(psi_bg),
  alpha = Re(e^{-i theta} u), beta = Im(e^{-i theta} u)
```

OCZEKIWANIE (jawne, pre-code):
- L_toy MA mod ujemny (pierscien V''(f)<0 przy r in ~(0.5, 2.6)
  wokol rdzenia wiaze stan w 2D) — czyli toy-sektor B2 jest wadliwy
  NAWET na tle stacjonarnym; diagnostyka, bez PASS/FAIL.
- L_full NIE ma szybkiego modu ujemnego (wir = minimum w klasie
  topologicznej; czlony fazowe V'(f)/f i sprzezenia znosza wiazanie).
  Jedyny dopuszczalny slad: powolne siodlo translacyjne pary
  (skala 2*pi*J/d^2 ~ 1e-3 rozmyta po normie modu -> |lambda| <~ 1e-4).

### 1.3 Rownanie #3 — sektor falowy i predykcja eikonalna (na przyszlosc)

Sektor falowy = JEDYNA nowa struktura B2 (d^2/dtau^2) zastosowana
do substratu B1 (pole zespolone): `d^2psi/dtau^2 = kappa*Lap(psi)
- V'(|psi|)*psi/|psi|`. Tlo wirowe jest jego dokladnym rozwiazaniem
statycznym. Fale amplitudowe (omega > m_TV = 0.9374) widza:

```text
n^2(r) = (omega^2 - V''(f(r))) / (omega^2 - V''(s*))
ogon:  n^2 - 1 = A_t / (r^2 (omega^2 - V''(s*))),
       A_t = V'''(s*) * s* * xi^2 = 3.8450 * 1.174166 * 0.56896
           = 2.5690   (V'''(s*) = -2b + 6c*s*)
n > 1 wszedzie (deficyt Phi jest osrodkiem GESTSZYM) i n(r) NIE
zalezy od znaku kretu -> refrakcja KU wirowi, UNIWERSALNA,
o ogonie ALGEBRAICZNYM 1/r^2 — przeciwienstwo obiektu B2!
alpha_pred(b, omega): ray tracing w Phase0_check (zaklad na nastepny
op pomiarowy). Fale fazowe (Goldstone, bezmasowe): c = sqrt(kappa)
niezalezne od f w wiodacym rzedzie — refrakcja amplitudowa ich nie
dotyczy; rozprazanie na cyrkulacji (typu Iordanskiego) jest
znakowo-nieparzyste i POZA zakresem tego cyklu.
```

## 2. Model i parametry (LOCKED)

```text
substrat: jak B1 (kappa=0.5, a=0.5, b=1.6, c=1.0; N=256, L=128,
  dx=0.5, periodycznie, seed=20260704, eps=0.30, S_GUARD=10 flaga)
relaksacja tla: przeplyw gradientowy dt=0.02 (jak B1)
dynamika falowa: leapfrog dt=0.05 (jak B2, D1: pozycyjny; Courant OK)
K1 (babel, pole rzeczywiste jak op-lock-interaction):
  ansatz s(r) = s* * 0.5*(1 - tanh((r - R)/1.5)); bisekcja R w [1, 6]
  do Delta_R <= 0.25/2; klasyfikacja po 4000 krokach przeplywu:
  Phi_max < eps -> kurczy sie; area(Phi>eps) > 1.5x start -> rosnie
K2 (wir-para): ansatz sinh-periodyzowany (D1 B1), wiry w
  (-32,-32,+1) i (+32,+32,-1); relaksacja 2000 krokow; pomiary
  stacjonarnosci w krokach 1000 i 2000
widma: przesunieta iteracja potegowa (shift 20), 4000 iteracji
  (L_toy — rzeczywisty; L_full — zespolony, iloczyn Re<u,v>);
  ROZDZIELCZOSC metody dla L_full ograniczona kontinuum Goldstone'a
  od 0: raportujemy lambda_min z deklarowana rozdzielczoscia 1e-3
noise-run (empiryczny test tachionu, analog D6b z B2):
  psi = psi_bg + 1e-6 * szum zespolony (seed), tau = 150 (3000
  krokow); gamma_fit = nachylenie ln||psi - psi_bg|| w oknach:
  rdzenie (r < 5 od kazdego rdzenia) i globalnie, okno tau in [50,150]
pakiet demonstracyjny (obserwacja): kanal amplitudowy,
  omega = 1.1, k0 = sqrt((omega^2 - V''(s*))/kappa) = 0.8141,
  sx=8, sy=10, u0=1e-3, start (x0=-60, y=-26) -> przelot kolo wiru
  (-32,-32) z b=6; znak ugiecia z pseudo-pedu poprzecznego delta-psi
  w oknie transmisji x > -10 (jak D3 B2)
```

## 3. Kryteria (LOCKED)

- **G1 techniczne:** determinizm bitowy (relaksacja tla x2), brak
  NaN/runaway, H_Gamma nierosnace w relaksacjach; w dynamice II rzedu
  dryft sekularny energii <= 1e-4 (definicja D4 z B2), oscylacja
  leapfrog raportowana.
- **G2 (K1 — babel krytyczny):** bisekcja rozstrzyga R_c z
  Delta_R <= 0.25; `R_c in [0.8, 1.2] * 3.0626`; lambda_min na
  zamrozonym bablu ~R_c: < 0 (mod wzrostu). PASS = pomiar
  rozstrzygniety (K1 sklasyfikowany jako siodlo). FALSYFIKATOR
  hipotezy pomocniczej: R_c poza pasmem -> prawo frontu nie
  ekstrapoluje sie do malych R (raport).
- **G3 (K2 — stacjonarnosc):** po 2000 krokach relaksacji:
  residuum `max|kappa*Lap(psi) - V'(|psi|)psi/|psi||` <= 1e-4;
  dryf pozycji rdzeni miedzy krokiem 1000 a 2000 <= 0.02;
  kret zachowany (+1/-1, n_tot=0). PASS/FAIL.
- **G4 (RDZEN — brak tachionu na tle wirowym):**
  (a) widmowo: lambda_min(L_full) >= -1e-3 (granica rozdzielczosci,
      zapisana z gory);
  (b) empirycznie: gamma_fit(noise-run, okno rdzeni) <= 0.01 ORAZ
      gamma_fit(globalnie) <= 0.01 ORAZ zero flag S_GUARD/NaN.
  Obie nogi wymagane. lambda_min(L_toy) RAPORTOWANE jako diagnostyka
  toy-sektora (oczekiwanie < -0.01; bez PASS/FAIL).
  FALSYFIKATOR: tachion takze na tle wirowym -> substrat nie ma
  stacjonarnego tla falowego w minimalnym zestawie.
- **G5 (klasyfikator kierunku, bez PASS/FAIL):** tabela
  alpha_pred(b, omega) z ray tracingu (Phase0_check) + znak ugiecia
  pakietu demonstracyjnego: `KU wirowi` (oczekiwane, uniwersalne) /
  `OD wiru` / `nierozstrzygniete`; pomiar ILOSCIOWY alpha = nastepny
  op z wlasnym LOCKiem (zaklad: tabela z Phase0_check).

## 4. Regula decyzyjna (LOCKED)

```text
G3 + G4 PASS
  -> TLO STACJONARNE ISTNIEJE: wir (stabilizacja topologia, nie sila).
     B2-prime odblokowane. Nastepny op (za autoryzacja): ilosciowy
     pomiar refrakcji na tle wirowym wg protokolu B2 (G2-G6 tamtego
     locka) z sektorem PELNYM zespolonym; alpha_pred z tego cyklu
     jako zalockowany zaklad. Jesli przy tym lambda_min(L_toy) < 0:
     toy-sektor rzeczywisty zdyskwalifikowany takze na tle
     stacjonarnym — B2-prime MUSI uzywac sektora pelnego (zapis
     do protokolu nastepnego opa).
G4 FAIL
  -> substrat bez stacjonarnego tla falowego w minimalnym zestawie;
     stop; analiza teoretyczna zrodla niestabilnosci przed dalszymi
     pomiarami propagacji.
G3 FAIL -> konfiguracja szachownicowa nie jest operacyjnie
     stacjonarna: raport (dryf/residuum), jeden dozwolony wariant
     techniczny: dluzsze relaksowanie x4 (zapisany z gory); jesli
     utrzymane -> stop.
G2: klasyfikacja K1 niezalezna (nie blokuje G3/G4).
```

## 5. Forbidden moves

1. Zero pinningu; kret czytany z pola; n_tot = 0.
2. Zadnych zmian progow/procedur po pierwszym runie (wyjatki
   zapisane z gory: jedno obnizenie dt; relaksacja x4 przy G3).
3. Nie deklarowac "soczewkowania grawitacyjnego" z demonstracji
   znaku (G5 to klasyfikator; ilosciowo — nastepny op).
4. tau != czas; zadnych deklaracji GR/MOND.
5. Nie promowac do core/ bez osobnej autoryzacji.

## 6. Czego ten cykl NIE robi

- Nie mierzy ilosciowo alpha(b, omega) (nastepny op, wlasny lock).
- Nie rozstrzyga uniwersalnosci znaku sily miedzy defektami
  (galaz C — osobna procedura).
- Nie dotyka watku solitonow korony (op-wall-dynamics 2026-07-03
  to INNY problem: stabilnosc mu/tau, CP-7).

## 7. Deliverables

- `Phase0_balance.md` (LOCKED przy utworzeniu)
- `Phase0_check.py`/`Phase0_output.txt` (R_c_pred, A_t, n^2(r),
  alpha_pred(b, omega), lambda_thr)
- `Phase1_critical_bubble.py`/`Phase1_output.txt` (G2 + czesc G1)
- `Phase2_vortex_background.py`/`Phase2_output.txt` (G3, G4 + G5-demo)
- `Phase_FINAL_close.md`, `README.md`
