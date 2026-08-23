---
title: "Phase0_balance — LOCK: blocked soliton bang toy model"
date: 2026-07-04
type: phase0-lock
tgp_owner: research/op-blocked-soliton-bang-2026-07-04
status: LOCKED
anti_lakatos_lock: PRESERVED
tags: [blocked-relaxation, soliton, substrate, toy-model, gravity-bridge]
related:
  - "[[../op-wall-dynamics-2026-07-03/README.md]]"
  - "[[../op-spectral-analysis-Phi-2026-07-03/README.md]]"
  - "[[../nbody/tgp_pde_solver.py]]"
  - "[[../galaxy_scaling/gs4_n_soliton_monte_carlo.py]]"
---

# Phase 0 — LOCK: Blokowana relaksacja N solitonow

## 0. Hipoteza robocza

Ten cykl testuje hipoteze autora z rozmowy 2026-07-04:

> Substrat nie musi bezposrednio produkowac metryki/grawitacji. Wystarczy,
> ze dopuszcza N lokalnych zaburzen typu soliton z wlasnym profilem Phi,
> ktore pojedynczo moga sie relaksowac, ale w ukladzie wielu obiektow
> wzajemnie blokuja pelna relaksacje. Grawitacja jest kandydacko energia
> / relacja wynikajaca z tej nierelaksowalnosci.

Wersja minimalna NIE twierdzi, ze policzona zostala grawitacja. Testuje
tylko, czy mechanizm "single decays, many block" ma stabilny rdzen
numeryczny i czy generuje nietrywialna energie interakcyjna.

## 1. Model zalockowany przed kodem

Uzywamy pola skalarnego `phi(x,y)` na siatce 2D jako toy-modelu substratu.
Proznia to `phi = 1`. Solitony sa lokalnymi zrodlami gaussowskimi z
pozycja `x_i` i amplituda `a_i`.

Ewolucja pola:

```text
dphi/dt = kappa * Lap(phi)
          - mobility(overlap) * m2 * (phi - 1)
          - lambda4 * (phi - 1)^3
          + source({x_i, a_i})
```

gdzie:

```text
overlap(x) = sum_i exp(-|x-x_i|^2 / (2 block_sigma^2))
mobility   = 1 / (1 + block_strength * max(overlap - 1, 0)^2)
```

Interpretacja:

- `m2*(phi-1)` jest lokalnym powrotem substratu do prozni.
- `mobility < 1` tlumi tylko kanal relaksacji do prozni w obszarach,
  gdzie profile wielu solitonow nachodza na siebie.
- Dyfuzja i samooddzialywanie nie sa blokowane; blokada nie jest zwyklym
  zamrozeniem calego solvera.

Ewolucja amplitud solitonow:

```text
da_i/dt = amp_rate * a_i * (overlap_i - decay_threshold)
```

gdzie `overlap_i` liczy wklad innych solitonow w promieniu `block_sigma`.
Samotny soliton (`overlap_i = 0`) traci amplitude. Soliton w grupie moze
sie utrzymac, jezeli sasiadujace profile dostarczaja wystarczajacej
blokady.

Regula kreacji w Phase 1 jest wylaczona. Phase 1 testuje tylko
utrzymanie/zanik juz zadanych solitonow. Regula kreacji moze byc dodana
w pozniejszym cyklu po osobnym locku, bo latwo wprowadzic artefakt
"perpetuum mobile" przez sam algorytm.

## 2. Parametry Phase 1 (LOCKED)

Jednostki bezwymiarowe:

- grid: `N = 64`, `L = 24`
- seed smoke: `12345`
- field steps: `steps = 1200`
- `dt = 0.025`
- `kappa = 0.18`
- `m2 = 0.55`
- `lambda4 = 0.12`
- `source_sigma = 0.75`
- `block_sigma = 2.40`
- `source_strength = 0.55`
- `block_strength = 18.0`
- `amp_rate = 0.035`
- `decay_threshold = 0.72`
- amplitude floor for "alive": `alive_threshold = 0.15`
- initial amplitude: `1.0`

Scenariusze smoke:

1. `single`: jeden soliton w centrum.
2. `cluster`: osiem solitonow w zwartej konfiguracji pierscien + centrum.
3. `pair_scan`: dwa solitony o separacjach `d = {2,3,4,6,8}` z amplitudami
   zamrozonymi, zeby zmierzyc energie interakcyjna.

## 3. Obserwable

Raportujemy:

- `N_alive(t)` — liczba solitonow z `a_i >= alive_threshold`.
- `amp_mean`, `amp_min`, `amp_max`.
- `E_positive = int [0.5*kappa|grad phi|^2 + 0.5*m2(phi-1)^2 + 0.25*lambda4(phi-1)^4] dA`.
- `E_functional = E_positive - int source*(phi-1) dA`.
- `blocked_fraction = mean(mobility < 0.5)`.
- `mean_overlap_at_solitons`.
- `V_pair(d) = E_functional(two at d) - 2 E_functional(single fixed)`.

## 4. Kryteria sukcesu i falsyfikacji

### C1 — deterministycznosc

Ten sam seed i parametry daja identyczne wyniki w jednym uruchomieniu.
Wykonanie skryptu bez wyjatku i bez NaN/Inf = PASS techniczny.

### C2 — samotny soliton ma kanal zaniku

Po `steps = 1200` scenariusz `single` ma:

```text
N_alive <= 0  albo  amp_mean < 0.15
```

Jesli samotny soliton pozostaje stabilny przy tych parametrach, mechanizm
"stabilnosc jako wlasnosc kolektywna" nie zostal pokazany w Phase 1.

### C3 — wiele solitonow blokuje relaksacje

Po `steps = 1200` scenariusz `cluster` ma:

```text
N_alive >= 4  oraz  amp_mean >= 0.35  oraz  blocked_fraction > 0.01
```

Jesli cluster zanika tak samo jak single, wynik jest negatywny dla toy
mechanizmu v0.

### C4 — energia interakcyjna jest nietrywialna

Dla `pair_scan` wymagamy:

```text
max(V_pair) - min(V_pair) > 1e-3
```

To NIE wymaga prawa `1/r`. Wymaga tylko, zeby overlap wielu profili
generowal mierzalna energie interakcyjna.

## 5. Forbidden moves

1. Nie zmieniac progow C1-C4 po pierwszym uruchomieniu Phase 1.
2. Nie stroic parametrow pod pozytywny wynik bez nowego Phase0 lock.
3. Nie interpretowac wyniku jako dowodu GR, MOND ani kosmologii.
4. Nie promowac wyniku do `core/` bez osobnej autoryzacji.
5. Raportowac wyniki negatywne wprost.
6. Nie ukrywac eksplozji numerycznej przez clamp jako sukcesu. Clamp moze
   chronic przed overflow, ale jesli `max(phi)` dobija do limitu, run
   jest oznaczany jako runaway.

## 6. Deliverables

- `README.md`
- `Phase1_blocked_relaxation_toy.py`
- `Phase1_output.txt`
- `Phase2_long_run_scan.py`
- `Phase2_output.txt`
- `Phase_FINAL_close.md`

## 7. Czego ten cykl NIE robi

- Nie rozstrzyga rzeczywistej stabilnosci solitonow TGP.
- Nie liczy pelnego pola 3D.
- Nie uzywa metryki jako inputu.
- Nie zawiera jeszcze dynamicznej kreacji nowych solitonow.
- Nie dotyka `core/` ani rejestrow predykcji.

## 8. Deklaracja LOCK

Model, parametry, scenariusze i kryteria C1-C4 zapisano przed wykonaniem
Phase 1. Wynik dodatni, mieszany albo negatywny bedzie raportowany bez
zmiany progow.
