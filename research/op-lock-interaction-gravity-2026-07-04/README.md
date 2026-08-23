---
title: "op-lock-interaction-gravity — oddzialywanie zlockowanych obiektow (most do grawitacji, krok 1)"
date: 2026-07-04
type: research-op
status: CLOSED
verdict: "PASS (G1-G4, G6); G5: FRONTOWE -> GALAZ B (Yukawa, nie-grawitacja)"
tags: [gravity-bridge, lock-interaction, yukawa, force-law, branch-B, closed]
related:
  - "[[Phase0_balance.md]]"
  - "[[Phase_FINAL_close.md]]"
  - "[[../op-bare-substrate-genesis-2026-07-04/Phase_FINAL_close.md]]"
---

# op-lock-interaction-gravity-2026-07-04

**Status: CLOSED (2026-07-04). LOCK zalozony i zachowany (anty-Lakatos
PRESERVED). Werdykt: G1-G4 PASS, G6 PASS, G5 = FRONTOWE -> GALAZ B.**

Kontynuacja po PASS `op-bare-substrate-genesis`. Pytanie cyklu: **czy dwa
zlockowane obiekty oddzialuja, jakim prawem, i czy to prawo ma cechy
grawitacji (dlugi zasieg, skalowanie z masa), czy sily powierzchniowej
(krotki zasieg, skalowanie z frontem)?**

## Rozstrzygniecie

- Oddzialywanie ISTNIEJE: przyciagajace, potwierdzone statycznie
  (`E_int(g) < 0`, monotoniczne) i dynamicznie (`Delta_v > 0`).
- Prawo: **Yukawa** `E_int(g) ~ -A e^(-m g)` wygrywa z forma potegowa
  z `Delta_AIC = +14.6` (prog 10); `m_fit = 0.883` (czulosc: 0.984)
  w pasmie zalockowanej predykcji substratu `m0 = sqrt(a/kappa) = 1.0`
  — prawo WYPROWADZONE algebraicznie z modelu, nie dofitowane.
  Niezalezne potwierdzenie z dynamiki: `m_dyn = 1.15`.
- Skalowanie: `ratio_mass = 1.14 <= 1.5` -> **FRONTOWE** (powierzchniowe,
  nie masowe/bulk).
- Uniwersalnosc: `ratio_hist = 0.989` — sila zalezy od stanu, nie od
  historii przygotowania (G6 PASS, wersja mocna).
- Kalibracja: `v(R) = v_inf - c/R`, `v_inf = 0.1599` (zgodnosc z frontem
  1D < 0.1%), `c = 0.490 ~ kappa`, `R^2 = 0.993`.

Wynik jest DOKLADNIE zgodny z oczekiwaniem teoretycznym zapisanym przed
kodem (Z2 bez modu bezmasowego -> Yukawa). Wg zalockowanej reguly
decyzyjnej: **bezposrednia sila miedzy lockami NIE jest grawitacja**
(krotki zasieg + skalowanie frontowe). Most przechodzi na — do wyboru
i autoryzacji autora jako osobne opy:

- **B1**: pole zespolone `s in C`, U(1), mod Goldstone'a = bezmasowy
  mediator dlugozasiegowy;
- **B2**: `Phi` jako metryka efektywna propagacji — "grawitacja" jako
  refrakcja w gradiencie `Phi` (droga geometryczna).

Szczegoly, liczby, errata i pelna regula: [[Phase_FINAL_close.md]];
LOCK i kryteria G1-G6: [[Phase0_balance.md]].
