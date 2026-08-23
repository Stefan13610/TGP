---
title: "op-nbody-additivity — addytywnosc parowa oddzialywania lockow (test superpozycji dla rdzenia)"
date: 2026-07-04
type: research-op
status: CLOSED
verdict: "PASS (G1-G4); G5: LOKALNY — superpozycja parowa poparta, zalamanie przewidywalne (Fermat)"
tags: [gravity-bridge, nbody, additivity, superposition, fermat-point, closed]
related:
  - "[[Phase0_balance.md]]"
  - "[[Phase_FINAL_close.md]]"
  - "[[../op-lock-interaction-gravity-2026-07-04/Phase_FINAL_close.md]]"
---

# op-nbody-additivity-2026-07-04

**Status: CLOSED (2026-07-04). LOCK zalozony i zachowany. Werdykt:
G1-G4 PASS, G5 = LOKALNY.**

Test pomostowy miedzy zmierzonym 2-cialowym prawem Yukawy
(`op-lock-interaction-gravity`) a zalozeniem liniowej superpozycji
zrodel w newtonowskiej granicy rdzenia TGP (sek08, tw. Newtona).
Pytanie: **czy E_int ukladow N zlockowanych obiektow jest suma parowa;
gdzie i jakim prawem sie zalamuje?**

## Rozstrzygniecie

- **Addytywnosc trzyma** w rezimie rozseparowanym: |delta| <= 3.3e-3
  dla g >= 3 (prog LOCKa: 0.02) — chain3, triangle3, hex6.
- **Zalamanie jest algebraiczne (mechanizm Fermata):** wiodaca poprawka
  to czlon 3-cialowy `-2b*s1*s2*s3` zlokalizowany w punkcie Fermata:
  znak ujemny (nadaddytywne przyciaganie) — PASS; |delta_tri| rosnie
  przy malejacym g — PASS; zbocze `mu_fit = 0.689` w pasmie [0.4, 1.1]
  (predykcja (sqrt(3)-1)*m0 = 0.73 asympt., 0.59 dla skonczonego R)
  — PASS; dyskryminator geometryczny: |delta_chain| ~ 1e-9-1e-10,
  o ~6 rzedow ponizej trojkata — PASS.
- **Dynamika (G5): LOKALNY** — `ratio_coop = 1.0000`; sciana odpowiada
  wylacznie na ogon po drugiej stronie wlasnej szczeliny (supresja
  trzeciego ciala ~e^(-(2R+2g)) ~ 1e-11, potwierdzona obserwacja hex6).
- Ciaglosc lancucha: refit Yukawy m = 0.813, R^2(log) = 0.979 (G2).

Konsekwencja dla rdzenia: zalozenie superpozycji zrodel w granicy
newtonowskiej NIE jest podwazone przez substrat w rezimie
rozseparowanym; poprawka nieaddytywna parametryzowalna jako
`delta ~ -C e^(-mu g)`, mu ~ 0.7, C ~ 3e-2. Program wraca na os
mostu: galezie B1 (Goldstone/U(1)) i B2 (Phi jako metryka efektywna)
— osobne opy, do autoryzacji.

Szczegoly, liczby i errata: [[Phase_FINAL_close.md]];
LOCK, predykcje i kryteria: [[Phase0_balance.md]].
