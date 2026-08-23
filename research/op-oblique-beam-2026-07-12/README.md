---
title: "op-oblique-beam (B2-prime, proba #5) — pelny P2: wiazka skosna atan(1/2) na zwalidowanej szachownicy L=256"
date: 2026-07-12
type: research-op
status: CLOSED (2026-07-14)
verdict: "KANAL CYRKULACYJNY ZMIERZONY (G5b FAIL 3/4 par, regula wzmacniajaca wariant (i), G6a PASS); P2 w brzmieniu 'slepota na kret' NEGATYWNIE; P1 bez zmian"
tags: [gravity-bridge, branch-B2prime, oblique-beam, checkerboard, circulation-channel, discovery, P2-negative, L256, closed]
related:
  - "[[Phase0_balance.md]]"
  - "[[Phase_FINAL_close.md]]"
  - "[[../op-asymmetric-lattice-2026-07-05/Phase_FINAL_close.md]]"
  - "[[../op-shortwave-lattice-2026-07-05/Phase_FINAL_close.md]]"
---

# op-oblique-beam-2026-07-12 (B2-prime, proba #5)

**Status: CLOSED (2026-07-14). Werdykt: G1 PASS; G2 BRAMKA PASS
(0.0000 przez tau=600 — zmierzona powtorka proby #3); G3 PASS
(mody skosne (2q,q), max 1.12%); G4 FAIL; G5a FAIL (3/8); G5b
FAIL 3/4 par wg pre-rejestrowanej reguly wzmacniajacej WARIANT
(i); G6a PASS (S-inv 0.056 deg = 1.3%; S-conj bitowo); G6b PASS.
-> KANAL CYRKULACYJNY ZMIERZONY. P2 w brzmieniu "kanal
geometryczny slepy na znak kretu" rozstrzygniete NEGATYWNIE;
status P1 bez zmian.**

Pytanie bylo: czy kanal geometryczny jest slepy na znak kretu
tam, gdzie zadna symetria tego nie wymusza? Odpowiedz: **NIE** —
na wiazce skosnej atan(1/2) (kierunek (2,1) poza wszystkimi
lustrami D4 szachownicy) reszta Delta_meas - Delta_eik jest
spojnie ujemna na wszystkich parach, przekracza 3 sigma_Delta
w 3/4 par kryterialnych i rozrzut plaszczyzn w kazdej z nich;
potwierdzona niezaleznie w obserwacji 1.1/24 i w sygnaturze
b_eff (+b odpychane, -b przyciagane, asymetria rosnie z omega).
G6a PASS wyklucza blad ramki: run wsteczny (-u) po tej samej
stronie wzgledem cyrkulacji odtwarza baseline (4.45 vs 4.39),
podczas gdy run przedni po przeciwnej stronie daje 5.47 — znak
kanalu podaza za strona przelotu wzgledem KRETU. G4/G5a FAIL to
skutek tego samego zjawiska (tlumienie alpha na galezi +b serii
1.5 az do odwrocenia znaku przy (1.5,+12): ratio -0.641); czesc
PARZYSTA pomiaru pozostaje w pasmie eikonalu we wszystkich 4
parach (post hoc). Zero deklaracji GR (forbidden move #6).

Szczegoly, komplet liczb, erraty: [[Phase_FINAL_close.md]];
LOCK i kryteria: [[Phase0_balance.md]].

## Wyniki kluczowe (skrot)

```text
G5b per para (Delta_meas / Delta_eik / reszta / 3sigma_Delta /
              rozrzut plaszczyzn):
  1.1/14: -1.082 / -0.301 / -0.781 / 0.852 / 0.356  OK (zero)
  1.1/20: -1.039 / -0.061 / -0.978 / 0.430 / 0.370  FAIL
  1.5/8:  -3.414 / -2.369 / -1.044 / 0.302 / 0.091  FAIL
  1.5/12: -3.725 / -0.473 / -3.252 / 0.241 / 0.086  FAIL
  (obs 1.1/24: -0.968 / -0.014 / -0.953 / 0.238 / 0.420 —
   ten sam znak, poza kryteriami)
G5a: 1.1: OK OK OK OK; 1.5: POZA(0.115) OK POZA(-0.641 ZNAK)
     POZA(2.189)
b_eff: +8->8.83, -8->6.94; +12->12.87, -12->10.99;
       +14->13.68, -14->13.42; +20->19.96, -20->19.71
G6a: S-inv |4.44750 - 4.39129| = 0.056 <= 0.220; S-conj 0.00e+00
```

## Notatki wykonawcze (doswiadczenie sesji #B2prime-5)

- Protokol utrzymany: LOCK + sekcja 10 (D1-D18) przed kodem;
  bramka G2 przed pierwszym pakietem; zero zmian progow po
  pierwszym runie; zestaw 4 par zamkniety; erraty w dniu runow
  (E1: przerwa API po LOCKu przed sekcja 10/kodem — zero wplywu;
  E2: przestoje miedzy fazami — determinizm bitowy; E3: S-inv
  na poziomie 0.056 deg, pre-rejestrowane D11a).
- Ramka skosna (xi/eta wrapowane wzgledem A z psi_ref) dziala:
  maska pasa + dyski r<12 zamrazane w wyzwoleniu; ray tracing
  z plaszczyznami xi_eval; wszystko zwalidowane G6a.
- Pakiet wsteczny D16 (faza -omega*dt, wyzwolenie xi_c < -10,
  okno lustrzane) — wzorzec dla przyszlych opow z pakietami
  w dowolnym kierunku.
- Delta_eik przy wlasnym b_eff POCHLANIA czesc kanalu odd (dryft
  poprzeczny przed wyzwoleniem zmienia b_eff!) — reszta i tak
  przekracza progi; op charakteryzacji powinien rozdzielic
  "odd w b_eff" od "odd w kacie" jawnie.
- Sygnal odd rosnie z omega (1.5 >> 1.1) i z |b| (1.5/12 > 1.5/8
  w reszcie) — projekt skanu charakteryzacji od tego zaczac.
