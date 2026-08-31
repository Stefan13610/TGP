---
title: "op-bloch-chain-stability — następca Q1: samouzgodniony łańcuch periodyczny + Bloch w sektorze tachionowym (ω²(n) modu runaway)"
date: 2026-08-31
type: research-cycle
status: CLOSED — Q-FAIL (PRIMARY wg rulingu; strict-literal: Q-INCONCLUSIVE)
tgp_owner: research/op-bloch-chain-stability-2026-08-31
authorization: "User 2026-08-31: 'ok działaj z krokami 1,2,3' (krok 3)"
related:
  - "[[Phase0_balance.md]]"
  - "[[Phase_method_decisions.md]]"
  - "[[Phase3_verdict_ruling.md]]"
  - "[[Phase_FINAL_close.md]]"
  - "[[NEEDS.md]]"
  - "[[../op-bath-two-sectors-2026-08-23/Phase_FINAL_close.md]]"
  - "[[../op-nonlinear-charge-constraint-2026-07-03/README.md]]"
---

# op-bloch-chain-stability (2026-08-31)

**Status: CLOSED (2026-08-31) — werdykt Q-FAIL** (odczyt PRIMARY wg
[[Phase3_verdict_ruling.md]]; strict-literal: Q-INCONCLUSIVE — oba
raportowane w [[Phase_FINAL_close.md]]): samouzgodniony łańcuch
periodyczny w sektorze tachionowym NIE stabilizuje — ω²_min(d) =
−1.222…−1.230 < 0 zbieżnie dla wszystkich istniejących teł
(d ∈ {3π, 4π, 6π}; dla d ∈ {π, 2π} tło niestałe nie istnieje),
ucieczka nieliniowa t_esc = 4.58–4.92 ≤ 2t*_ref przy obu znakach;
kontrole P1a/P3b/P3c czyste.

Następca pytania Q1 z `op-bath-two-sectors` (Q1-INCONCLUSIVE — komórka
zero-flux z obciętym tłem niezdolna rozstrzygnąć), zbudowany DOKŁADNIE wg
NEEDS N2 poprzednika: tło SAMOUZGODNIONE (relaksacja przed linearyzacją)
+ periodyczne BC + analiza pasmowa Blocha — konstrukcja, która w Q2
poprzednika dała czystą zbieżność, przeniesiona na sektor dynamiczny.

**Pytanie binarne:** czy istnieje separacja d łańcucha periodycznego
w sektorze tachionowym, przy której ω²_min(d) > 0 (najniższa gałąź
pasmowa fluktuacji wokół tła samouzgodnionego)?

Kryteria, punkty skanu, kontrole i forbidden moves: [[Phase0_balance.md]]
(LOCK zamknięty i zacommitowany PRZED jakimkolwiek obliczeniem).

## Log faz

- 2026-08-31: Phase 0 LOCK zapisany. Obliczenia: NIEROZPOCZĘTE.
- 2026-08-31: `Phase_method_decisions.md` FROZEN przed startem (forma EL
  + rozjazd z #63: cykl w akcji kanonicznej K=g⁴, U=g⁷/7−g⁸/8; f_ε → K_ε
  tylko w ewolucji; Newton na pół-komórce; Bloch węzłowy z fazą e^{ikd};
  solver scipy.linalg.eigh; reguła Goldstone'a; metryki błędu).
- 2026-08-31: **Phase 1 PASS** (`Phase1_output.txt`): P1a tach/stab
  maxerr 8.3e−5 / 3.2e−5, ratio 4.000/4.000 (rząd 2); P1b |ΔE|/E ≤ 6.9e−15,
  dt-konsystencja ≤ 1.3e−13.
- 2026-08-31: **Phase 2** (`Phase2_output.txt` + addendum): istnienie dla
  d ∈ {3π, 4π, 6π} (‖R‖∞ ≤ 5.7e−12, ‖g_N−g_2N‖∞ ≤ 4.4e−5); d ∈ {π, 2π}:
  kolaps do próżni (strukturalny — okres orbit ≥ 2π, potwierdzone S3).
  Incydent: start zalockowany dla 6π dał tło 2-garbne → DODANY start
  S3single (jednogarbne, g_min=0.1675); nic nie usunięto.
- 2026-08-31: **Phase 3** (`Phase3_output.txt`): ω²_min(d) ZBIEŻNE 4/4:
  3π: −1.222191; 4π: −1.229340; 6π(2-garb): −1.222209; 6π(1-garb):
  −1.229829 — wszystkie ujemne, argmin k=0, mod amplitudowy
  (nie-translacyjny; Goldstone λ~−1e−5 zidentyfikowany osobno);
  P3b PASS 5/5 (maxerr ≤ 2.7e−5); P3c PASS (+1.88310/+1.88311, kotwica
  poprzednika odtworzona); cross-check u-formy ≤ 8.8e−6.
  `Phase3_verdict_ruling.md` (kwantyfikator d przy istnieniu częściowym)
  zapisany PRZED Phase 4.
- 2026-08-31: **Phase 4** (`Phase4_output.txt`): ucieczka we WSZYSTKICH
  16 biegach t_esc = 4.58–4.92 ≤ 2t*_ref = 7.24 (oba znaki; kontrole
  ε=0.1 i dt/2 stabilne); σ_fit vs √|ω²_min| odchył ≤ 2.6%; gate energii
  ≤ 1.6e−7 PASS.
- 2026-08-31: **CLOSED — Q-FAIL** ([[Phase_FINAL_close.md]]);
  [[NEEDS.md]] (user-gated: Limitations 1D + podniesiony priorytet 3D —
  argument węzłowy Hilla czyni negatyw 1D strukturalnym).
