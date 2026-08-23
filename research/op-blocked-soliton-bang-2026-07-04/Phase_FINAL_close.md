---
title: "Phase_FINAL_close — op-blocked-soliton-bang"
date: 2026-07-04
type: final-close
status: CLOSED-EXECUTED
anti_lakatos_lock: PRESERVED
---

# Phase FINAL — close

## 1. Co zostalo wykonane

Utworzono izolowany research-op:

```text
research/op-blocked-soliton-bang-2026-07-04/
```

Zawiera:

- `Phase0_balance.md` — lock hipotezy, modelu i kryteriow.
- `Phase1_blocked_relaxation_toy.py` — minimalny toy 2D.
- `Phase1_output.txt` — smoke run.
- `Phase2_long_run_scan.py` — skan czasu i blokowania.
- `Phase2_output.txt` — wynik skanu.
- `README.md` — werdykt cyklu.

Nie edytowano `core/`, rejestrow predykcji ani planu Cursor.

## 2. Werdykt wzgledem Phase0

Phase 1:

```text
C1 PASS  deterministic finite execution
C2 FAIL  single nie spadl ponizej progu przy steps=1200
C3 PASS  cluster przezywa i ma strefe blokady pola
C4 PASS  energia pary jest nietrywialna
```

Formalny wynik Phase1: **mixed, 3/4 PASS**.

Phase 2 pokazal, ze C2 FAIL jest przede wszystkim skala czasu:

```text
single @ 1200 steps: amp_mean = 0.469429, C2 FAIL
single @ 2400 steps: amp_mean = 0.220364, C2 FAIL
single @ 3600 steps: amp_mean = 0.103445, C2 PASS
cluster @ 3600 steps: N_alive = 8/8, C3 PASS
```

Wersja fizyczna hipotezy "samotny soliton jest realnie niestabilny" jest
zgodna z Phase2: soliton samotny zanika, tylko wolniej niz zalockowany
czas smoke-runu.

## 3. Co jest najciekawsze

Energia pary reaguje na kanal field-blocking:

```text
block_strength=0      V_span=0.0148794478
block_strength=6      V_span=0.0765983378
block_strength=18     V_span=0.2479631809
block_strength=36     V_span=0.3864610308
```

To jest najwazniejszy sygnal mostu: oddzialywanie nie jest tylko
ksiegowaniem amplitud solitonow. Zalezy od blokady relaksacji pola.

## 4. Czego NIE wolno nadinterpretowac

Ten cykl nie pokazuje jeszcze:

- genezy zaburzen z golego substratu relacyjnego `Gamma`;
- samorzutnego powstawania solitonow z dynamiki `s_i`;
- prawa `1/r`;
- metryki GR;
- MOND;
- realnej kosmologii;
- dynamicznej nukleacji nowych solitonow;
- zerowego kosztu kreacji w stanie blokowanym.

To jest tylko pierwszy pozytywny toy-test, ze mechanizm "single decays,
many block" nie jest natychmiast sprzeczny numerycznie.

## 5. Post-close korekta: model byl za pozny

Po zamknieciu cyklu autor doprecyzowal, ze fizycznie wazniejszy jest
etap przedsolitonowy:

```text
1. goly substrat relacyjny Gamma;
2. zaburzenie s_i != 0, jeszcze nie pelnoprawny soliton;
3. lokalne Phi_i = s_i^2 > 0 generuje zalazek metrycznosci;
4. w tej metrycznosci rozpad/kreacja zaburzen sa bezkosztowe;
5. kaskada zaburzen;
6. lock: dopiero teraz mowimy o solitonach, a kreacja/relaksacja zamyka sie.
```

Ten research-op zaczal od punktu 6. Zatem jego wynik pozostaje wazny jako
test poznej fazy locku, ale NIE jest testem genezy substratowej.

Szczegoly post-close zapisano w `POST_CLOSE_BRAINSTORM_UPDATE.md`.

## 6. Najwazniejszy next step

Nastepny research-op powinien dotyczyc poprawionego problemu:

```text
research/op-bare-substrate-genesis-2026-07-04/
```

Model powinien startowac od zmiennych substratowych `s_i` na grafie
`Gamma`, a nie od listy gotowych solitonow. Minimalna dynamika:

```text
ds_i/dtau = - dH_Gamma/ds_i + seed/noise + nonlinear metric feedback
Phi_i = s_i^2
metric_enabled_i = Phi_i > epsilon
```

gdzie `tau` jest parametrem selekcji/relaksacji, nie czasem fizycznym
przed powstaniem fazy metrycznej.

Osobny modul powinien mierzyc to, co autor nazwal kluczowe:

> w stanie blokowanym koszt tworzenia nowych solitonow powinien byc zerowy.

Minimalny test po korekcie:

```text
1. Start: Gamma + male zaburzenie s_i, bez gotowych solitonow.
2. Policz lokalny koszt utworzenia kolejnego zaburzenia DeltaE_create(x).
3. Sprawdz, czy DeltaE_create -> 0 w obszarze lokalnej metrycznosci.
4. Pozwol na kaskade tylko tam, gdzie koszt jest numerycznie zerowy.
5. Sprawdz, czy kaskada konczy sie lockiem: brak dalszej kreacji/relaksacji.
```

To oddzieli fizyczna "zerokosztowa kreacje" od algorytmicznego spawn-rule.

## 7. Status koncowy

Status: **CLOSED-EXECUTED, toy-positive / theory-open**.

Najkrotszy werdykt:

> Substratowy model blokowanej relaksacji potrafi numerycznie rozroznic
> samotny niestabilny soliton od stabilizowanego klastra i generuje
> mierzalna energie interakcyjna. Po korekcie autora: jest to test
> poznej fazy locku, nie test genezy z golego substratu. Wlasciwy
> nastepny test musi zaczynac od `Gamma/s_i`, bez wrzucania solitonow.
