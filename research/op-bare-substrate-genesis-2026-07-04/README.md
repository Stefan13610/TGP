---
title: "op-bare-substrate-genesis — geneza samopodtrzymujacych struktur z golego substratu"
date: 2026-07-04
type: research-op
status: CLOSED
verdict: PASS (G1-G6)
tags: [substrate-genesis, soliton, lock, nucleation, gravity-bridge, closed]
related:
  - "[[Phase0_balance.md]]"
  - "[[Phase_FINAL_close.md]]"
  - "[[../op-blocked-soliton-bang-2026-07-04/README.md]]"
  - "[[../op-lock-interaction-gravity-2026-07-04/Phase0_balance.md]]"
---

# op-bare-substrate-genesis-2026-07-04

**Pytanie cyklu**: czy goly substrat (`Phi = 0`, faza niemetryczna) potrafi
SAM — z jednej dynamiki pola, bez gotowych solitonow i bez osobnej reguly
przezycia — wygenerowac lokalna, samopodtrzymujaca sie strukture `Phi > 0`
przez kolektywny lock nakladajacych sie zaburzen?

**Odpowiedz**: TAK. `G1-G6 PASS`. Szczegoly: [[Phase_FINAL_close.md]].

## Wynik w trzech liniach

1. **Single decays / many lock**: pojedynczy seed potrzebuje amplitudy
   A_c in (4.0, 5.0) ~ 4x powyzej prawdziwej prozni; klaster K=7 takich
   samych seedow lockuje juz przy A0_c ~ 0.46-0.57 (przewaga kolektywna
   >= 7-8.7x). Prog jest ostra bifurkacja (szerokosc 0.0094), potwierdzona
   na drobniejszej siatce (bez pinningu).
2. **Koszt kreacji**: DeltaE_create = -0.0146 na GRANICY regionu
   metrycznego (samoczynna), +0.0144 w martwym substracie, +0.0228 we
   wnetrzu locku. Kreacja jest "za darmo" tylko na froncie metryki.
3. **Zastrzezenie**: lock to rosnaca domena prawdziwej prozni, nie
   skonczony soliton; `tau` to parametr relaksacji, nie czas; zero
   twierdzen o GR.

## Model (LOCK)

Siatka 2D 128x128 (dx=0.5, tor), `s(x)` w R, `Phi = s^2`,
`ds/dtau = kappa*Lap(s) - V'(s)`,
`V(s) = 0.5*a*s^2 - (b/3)|s|^3 + 0.25*c*s^4`, a=0.5, b=1.6, c=1.0,
kappa=0.5. `s=0` metastabilne (bariera |s|~0.426), prawdziwa proznia
s*~1.174, V(s*)~-0.044. Zero clampow, zero regul amplitud, zero spawnow.

## Pliki

| Plik | Rola |
|---|---|
| [[Phase0_balance.md]] | LOCK modelu, kryteria G1-G6, poprawki pre-code (sekcja 11) |
| `Phase0_check.py` -> `Phase0_output.txt` | weryfikacja rachunkow analitycznych (pre-code) |
| `Phase1_substrate_genesis.py` -> `Phase1_output.txt` | bare / single / cluster; G1-G4 |
| `Phase2_threshold_and_create_cost.py` -> `Phase2_output.txt` | skan progu, pinning N=256, DeltaE_create; G5-G6 |
| [[Phase_FINAL_close.md]] | zamkniecie, errata, decyzja |

Reprodukcja (deterministyczna, seed 20260704):

```bash
python Phase0_check.py
python Phase1_substrate_genesis.py
python Phase2_threshold_and_create_cost.py
```

## Nastepny cykl

[[../op-lock-interaction-gravity-2026-07-04/Phase0_balance.md]] — most do
grawitacji: oddzialywanie pol miedzy zlockowanymi obiektami, z zalockowana
predykcja algebraiczna i ostrym dyskryminatorem "sila frontowa vs sila
masowa". Status: DRAFT-FOR-LOCK.
