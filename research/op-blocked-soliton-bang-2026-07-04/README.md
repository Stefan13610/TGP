---
title: "op-blocked-soliton-bang — blokowana relaksacja N solitonow"
type: research_cycle
status: CLOSED-EXECUTED
phase: FINAL
folder_status: closed-mixed-positive-toy
claim_status: "Toy model wykonany: Phase1 C1/C3/C4 PASS, C2 FAIL przy zalockowanych 1200 krokach; Phase2 pokazuje, ze samotny soliton przekracza prog zaniku przy 3600 krokach, a klaster pozostaje blokowany. POST-CLOSE BRAINSTORM: ten op testuje tylko pozna faze locku gotowych solitonow; NIE testuje wlasciwej genezy z golego substratu relacyjnego Gamma/s_i."
created_date: 2026-07-04
closed_date: 2026-07-04
anti_lakatos_lock: PRESERVED
related:
  - "[[../op-wall-dynamics-2026-07-03/README.md]]"
  - "[[../op-spectral-analysis-Phi-2026-07-03/README.md]]"
  - "[[../nbody/tgp_pde_solver.py]]"
  - "[[POST_CLOSE_BRAINSTORM_UPDATE.md]]"
---

# op-blocked-soliton-bang

## Cel

Pierwszy falsyfikowalny toy-test hipotezy:

```text
samotny soliton -> realnie niestabilny / relaksuje sie
N solitonow     -> overlap profili Phi blokuje relaksacje
blokada         -> mierzalna energia interakcyjna, kandydat mostu do grawitacji
```

Metryka nie jest inputem. W tym cyklu obserwujemy tylko pole `phi`,
amplitudy solitonow, blokade relaksacji i energie interakcyjna.

## Metoda

- **Phase 0** ([[Phase0_balance.md]]): zamknieto model, parametry,
  kryteria C1-C4 i forbidden moves przed uruchomieniem kodu.
- **Phase 1** ([[Phase1_blocked_relaxation_toy.py]] ->
  [[Phase1_output.txt]]): minimalny model 2D `64x64`, `steps=1200`.
- **Phase 2** ([[Phase2_long_run_scan.py]] ->
  [[Phase2_output.txt]]): skan czasu, progu zaniku i sily blokowania.

## Wyniki Phase 1

Zalockowane kryteria:

| Kryterium | Wynik | Szczegol |
|---|---:|---|
| C1 deterministic finite execution | PASS | brak NaN/Inf, brak runaway |
| C2 single soliton decay channel | FAIL | po 1200 krokach `amp_mean=0.469429`, nadal powyzej `0.15` |
| C3 cluster blocked survival | PASS | `N_alive=8/8`, `amp_mean=2.5`, `blocked=0.133789` |
| C4 nontrivial pair interaction | PASS | `V_pair_span=0.2479631809` |

Phase 1 werdykt: **3/4 PASS, 1/4 FAIL**. Zgodnie z lockiem C2
pozostaje formalnie negatywne dla czasu 1200.

## Wyniki Phase 2

Phase2 nie przepisuje C1-C4; diagnozuje, co znaczy C2 FAIL.

Najwazniejsze liczby:

- `single` przy 1200: `amp_mean=0.469429`, C2 FAIL.
- `single` przy 2400: `amp_mean=0.220364`, C2 FAIL.
- `single` przy 3600: `amp_mean=0.103445`, `N_alive=0/1`, C2 PASS.
- `cluster` przy 1200/2400/3600: `N_alive=8/8`, C3 PASS.

Interpretacja: samotny soliton faktycznie ma kanal zaniku, ale zalockowany
czas Phase1 byl zbyt krotki, zeby przekroczyc prog `alive_threshold`.

Skan `block_strength` pokazuje rozdzial dwoch efektow:

- przy `block_strength=0` klaster nadal ma amplitudy, ale C3 FAIL, bo brak
  blokady pola (`blocked_fraction=0`);
- przy `block_strength >= 6` klaster przechodzi C3;
- `V_pair_span` rosnie z `block_strength`:
  `0.014879` -> `0.076598` -> `0.247963` -> `0.386461`.

To oznacza, ze energia pary jest czula na fizyczny kanal blokowanej
relaksacji pola, a nie tylko na arbitralna regule amplitud.

## Werdykt

Toy-mechanizm **ma nogi numerycznie** w wersji minimalnej:

1. samotny soliton jest niestabilny w dluzszym czasie;
2. klaster pozostaje dluzywotny;
3. pole ma niezerowa strefe blokady;
4. energia interakcyjna zalezy od separacji i sily blokowania.

Ale wynik jest tylko **toy-positive / theory-open**. Nie jest to jeszcze
dowod grawitacji, MOND, kosmologii ani pelnej metryki.

## Post-close korekta interpretacji

Po zamknieciu cyklu autor doprecyzowal, ze wlasciwy problem genezy jest
wczesniejszy niz ten toy-model. Ten op zaczyna od juz zadanych solitonow
i pyta, czy klaster blokuje relaksacje. Poprawny program genezy powinien
zaczynac od **golego substratu relacyjnego**:

```text
Gamma + lokalne zmienne s_i
-> male zaburzenie s_i != 0
-> Phi_i = s_i^2 > 0 wlacza lokalna metrycznosc
-> w tej lokalnej metryce kreacja/rozpad zaburzen sa bezkosztowe
-> kaskada zaburzen
-> lock: dopiero wtedy powstaja stabilne solitony
```

W tej poprawionej interpretacji niniejszy cykl dotyczy tylko ostatniego
odcinka:

```text
gotowe solitony -> blokowana relaksacja -> energia interakcyjna
```

Nie testuje on, czy solitony powstaja samorzutnie z rownania substratu.
Szczegoly zapisano w [[POST_CLOSE_BRAINSTORM_UPDATE.md]].

## Kluczowa luka: zerokosztowa kreacja

Autor doprecyzowal w trakcie cyklu:

> wazny jest mechanizm, w ktorym koszty tworzenia nowych solitonow sa
> zerowe w takim stanie.

Ten cykl tego nie liczyl. Phase0 celowo wylaczyl kreacje nowych solitonow,
zeby nie ukryc artefaktu algorytmicznego. Nastepny cykl powinien dodac
osobny, zalockowany kanal:

```text
stan blokowany N solitonow
-> lokalny koszt nukleacji DeltaE_create(x)
-> test, czy DeltaE_create -> 0 na granicy/metrycznej frontiere
-> kreacja nowych solitonow bez naruszania energii w fazie stabilnej
```

Po post-close korekcie jest to jeszcze mocniejsze: nastepny test nie
powinien dodawac probe-solitonu na tle gotowego klastra jako glowny
model. Powinien zaczac od `s_i` na `Gamma` i sprawdzic, czy zaburzenia
same przechodza przez faze bezkosztowej kreacji, kaskady i locku.
Bez tego most do "wielki wybuch trwa na granicy metryki" pozostaje
hipoteza robocza, nie wynikiem numerycznym.

## Pliki

- [[Phase0_balance.md]]
- [[Phase1_blocked_relaxation_toy.py]]
- [[Phase1_output.txt]]
- [[Phase2_long_run_scan.py]]
- [[Phase2_output.txt]]
- [[Phase_FINAL_close.md]]
- [[POST_CLOSE_BRAINSTORM_UPDATE.md]]
