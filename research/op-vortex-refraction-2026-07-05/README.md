---
title: "op-vortex-refraction (B2-prime) — ilosciowa refrakcja fal na defekcie: most geometryczny wokol materii-defektow"
date: 2026-07-05
type: research-op
status: CLOSED
verdict: "G2 PASS; G3 FAIL strukturalne (tlo — para wirow — anihiluje pod dynamika II rzedu przy tau=107.5, w oknie przelotu pakietu; 9/9 runow niewaznych) -> STOP; G4/G5/G6 NIEROZSTRZYGNIETE; P1/P2 otwarte — pomiar wymaga tla o dlugim czasie zycia"
tags: [gravity-bridge, branch-B2prime, vortex-lensing, refraction, eikonal, universality, background-collapse, closed]
related:
  - "[[Phase0_balance.md]]"
  - "[[Phase_FINAL_close.md]]"
  - "[[../op-stationary-background-2026-07-05/Phase_FINAL_close.md]]"
  - "[[../op-phi-metric-refraction-2026-07-04/Phase_FINAL_close.md]]"
---

# op-vortex-refraction-2026-07-05 (B2-prime)

**Status: CLOSED (2026-07-05). Werdykt i pelne uzasadnienie:
[[Phase_FINAL_close.md]]. LOCK zachowany (anty-Lakatos PRESERVED).**

Pytanie brzmialo: **czy ugiecie pakietu amplitudowego na wirze
zgadza sie ilosciowo z eikonalem substratu (P1) i czy jest
UNIWERSALNE w znaku kretu (P2)?**

## Wynik w jednym akapicie

Pomiar okazal sie NIEWYKONALNY na zaplanowanym tle: para wirow
(+1,-1), kwazi-stacjonarna pod przeplywem gradientowym (dryf
0.015/tau), pod PELNA dynamika II rzedu kolabuje INERCYJNIE
(bez dyssypacji sila akumuluje predkosc) i anihiluje przy
tau = 107.5 — wewnatrz okna przelotu pakietu. Kolaps jest fizyczny
(kontrola dt/2: roznica 0.0074). Pre-rejestrowany budzet dryfu
(~3 na przelot) sfalsyfikowany ~13x. Wg reguly decyzyjnej: STOP
przed interpretacja; P1/P2 pozostaja OTWARTE. Strona teoretyczna
jest gotowa: eikonal na rzeczywistym tle zgodny z tabela cytowana
(<= 11% w punktach kryterialnych, b=16 dodane), kalibracja
dyspersji G2 PASS (1.21%).

## Pliki

- [[Phase0_balance.md]] — LOCK + doprecyzowania pre-code D1-D13
- `Phase0_check.py` / `Phase0_output.txt` — k0/v_g, n^2(r) na tle
  pary, eikonal na rzeczywistym tle vs cytowany, budzet dryfu
- `Phase1_calibration.py` / `Phase1_output.txt` — G2 (dyspersja)
  + background_check + determinizm
- `Phase2_deflection.py` / `Phase2_output.txt` — skan glowny:
  9/9 runow niewaznych (kolaps tla)
- `Phase2b_background_lifetime.py` / `Phase2b_output.txt` —
  diagnostyka kolapsu + obserwacja poza kryteriami
- [[Phase_FINAL_close.md]] — zamkniecie, errata, nastepne kroki

## Notatki wykonawcze dla nastepnych sesji (doswiadczenie tego cyklu)

Wszystkie pulapki z poprzednich opow potwierdzone (ansatz zlozony,
roznica centralna, usrednianie wektorowe, detektor windingu, zero
`cd`). NOWE doswiadczenie:

- **Dryf gradientowy NIE przewiduje dryfu falowego.** Kwazi-
  stacjonarnosc zmierzona przeplywem gradientowym (predkosc ~ sila)
  nie przenosi sie na dynamike II rzedu (przyspieszenie ~ sila):
  para o separacji 90 na L=128 anihiluje w tau~108. Kazdy przyszly
  op z tlem defektowym MUSI zmierzyc czas zycia tla POD DYNAMIKA
  POMIAROWA w Phase 1 (test: ewolucja II rzedu bez pakietu,
  pozycje rdzeni vs tau) ZANIM zaplanuje okno pomiaru.
- **Prog okna podazajacy za rdzeniem** (x > x_rdzen+10 z pozycja
  biezaca) na ruchomym tle ucieka przed pakietem — wyzwolenie
  nigdy nie nastepuje. Definiowac okna wzgledem pozycji rdzenia
  w chwili najwiekszego zblizenia pakietu albo w ukladzie
  wspolporuszajacym sie.
- **Referencja blizniacza (psi_ref w lockstep, D3)** dziala
  technicznie bardzo dobrze (E_dryft 1.5e-07, czyste delta_psi)
  — zalecana w kazdym przyszlym pomiarze pakietowym.
- Kandydaci na tlo nastepnego opa: L=256 (czas kolapsu ~liniowy
  w separacji — przeliczyc w Phase0!), uklad wspolporuszajacy sie,
  konfiguracje o zerowej sile netto (analiza zakazu pinningu
  wymagana). Szczegoly: [[Phase_FINAL_close.md]], sekcja 6.
  **Draft nastepnego opa gotowy:**
  [[../op-lattice-background-2026-07-05/Phase0_balance.md]]
  (DRAFT-FOR-LOCK, czeka na autoryzacje).

**Kluczowe liczby zakladu (bez zmian, zrodlo:
`../op-stationary-background-2026-07-05/Phase0_output.txt`):**
s* = 1.174166, Phi* = 1.378665, V''(s*) = 0.878665, m_TV = 0.9374,
xi^2 = 0.5690, A_t = 2.5690; alpha_pred(omega=1.1): +43.0 / +14.0 /
+5.30 deg (b = 6/8/12); NOWE (to tlo, Phase0 tego opa): +39.1 /
+14.6 / +5.0 / +2.7 deg (b = 6/8/12/16).
