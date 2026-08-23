---
title: "op-phi-metric-refraction — Phi jako metryka efektywna propagacji (most do grawitacji, galaz B2)"
date: 2026-07-04
type: research-op
status: CLOSED
verdict: "G1-G2 PASS; G3/G4/G6 NIEAPLIKOWALNE (tachion sciany tla zamrozonego, pre-rejestrowany w D6 i zmierzony: gamma=0.125); B2 wymaga tla stacjonarnego"
tags: [gravity-bridge, branch-B2, effective-metric, refraction, eikonal, wall-tachyon, closed]
related:
  - "[[Phase0_balance.md]]"
  - "[[Phase_FINAL_close.md]]"
  - "[[../op-lock-interaction-gravity-2026-07-04/Phase_FINAL_close.md]]"
  - "[[../op-goldstone-mediator-2026-07-04/README.md]]"
---

# op-phi-metric-refraction-2026-07-04 (galaz B2)

**Status: CLOSED (2026-07-04). LOCK autoryzowany i zachowany
(anty-Lakatos PRESERVED).** Pelne zamkniecie: [[Phase_FINAL_close.md]].

**Wynik jednym zdaniem:** zamrozone tlo zlockowanego obiektu jest
TACHIONOWE w scianie (pasmo V''(s)<0; lambda_min=-0.0152,
gamma_meas=0.125 — zgodnosc rachunku z empiria 1.3%), wiec pomiar
refrakcji na tle zamrozonym jest niemierzalny (0/12 waznych runow wg
pre-rejestrowanej klauzuli D6); dyspersja skalibrowana (G2 PASS,
blad <=1.5%), a zalockowany eikonal przewiduje znak ugiecia ZALEZNY
od b (OD obiektu dla b<=12, KU obiektowi przy b~16 przez pierscien
n^2>1 w scianie) — predykcje czekaja na tlo STACJONARNE
(op-wall-dynamics) jako warunek powrotu do B2.

Kontynuacja galezi B z `op-lock-interaction-gravity`. Pytanie B2:
**czy male zaburzenia propagujace sie na tle zlockowanego obiektu
Phi(x) doznaja refrakcji — w ktora strone i czy kat jest policzalny
(eikonal) z samego substratu?**

## Konstrukcja rozstrzygniecia

- JEDYNE rozszerzenie modelu (zadeklarowane wprost): sektor falowy
  malych zaburzen `d^2u/dtau^2 = kappa*Lap(u) - V''(s_bg)*u` na
  ZAMROZONYM tle; sektor konserwatywny (H_u zachowane — kryterium G1).
- Dwa jawnie oddzielone pytania: P1 policzalnosc (kat vs eikonal,
  pasmo [0.5, 2.0]) i P2 grawitacyjnosc (znak: ku obiektowi / od).
- Uczciwe oczekiwanie pre-code (analog zalockowanej Yukawy z
  poprzedniego opa): poniewaz `V''(s*) > V''(0)`, obiekt jest osrodkiem
  "rzadszym" (n < 1) -> refrakcja OD obiektu (antygrawitacyjna).
  Taki wynik = zmierzone wskazanie, ze "Phi jako metryka" wymaga
  sprzezenia `kappa(Phi)/c(Phi)` (aksjomat rdzenia) jako dodatkowej
  struktury — definiuje nastepny op. Ugiecie KU obiektowi =
  niespodzianka i mocny kandydat na most geometryczny.
- Kalibracja: dyspersja `omega^2 = kappa*k^2 + V''` na obu tlach
  (m_FV = 0.707, m_TV = 0.937); skan omega in {1.0, 1.1, 1.3},
  parametry zderzenia b in {8, 12, 16, 20}; test liniowosci u0 -> u0/2.

Szczegoly, kryteria G1-G6 i pelna regula decyzyjna:
[[Phase0_balance.md]]. Projekt siostrzany (galaz B1):
[[../op-goldstone-mediator-2026-07-04/Phase0_balance.md]].
