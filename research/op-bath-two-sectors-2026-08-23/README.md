---
title: "op-bath-two-sectors — runaway w kąpieli na fazach zmierzonych (Q1) + derywacyjny test dwóch sektorów znaku W (Q2)"
type: research_cycle
status: PHASE0-LOCKED
phase: "0 (LOCK zamknięty; zero obliczeń)"
folder_status: locked-awaiting-execution
claim_status: "PHASE0-LOCKED 2026-08-23 — następca op-lattice-bath-runaway (CLOSED-GATE-STOP) po rozstrzygnięciu N1/N2. Q1: czy mod runaway dostaje ω²>0 przy skończonej gęstości n, z ogonami/fazami ZMIERZONYMI z układów M-P (potęgowy, EL Noty kanonicznej) i M-L (logarytmiczny eq:J-ode, nośnik Δ=120.01°). Q2 (hipoteza autora): czy znak tachionowy efektywnego W EMERGUJE z gęstości w akcji STABILNEJ (eq:field-eq-reproduced) — test 'pojedynczy obiekt vs kąpiel' jako dwóch sektorów jednej akcji. Kryteria Q1/Q2-PASS/FAIL/INCONCLUSIVE, kontrole negatywne (Yukawa bez cos; próżnia w komórce; d=∞→Yukawa m²=+γ) zalockowane PRZED kodem."
created_date: 2026-08-23
closed_date: null
authorization: "User 2026-08-23: 'ok, przygotuj prompt' — realizacja Phase 1–3: nowa sesja, osobny agent (prompt: HANDOFF_PROMPT.md)"
anti_lakatos_lock: ACTIVE
related:
  - "[[Phase0_balance.md]]"
  - "[[HANDOFF_PROMPT.md]]"
  - "[[../op-lattice-bath-runaway-2026-08-23/README.md]]"
  - "[[../op-lattice-bath-runaway-2026-08-23/ANALIZA_N1_pochodzenie-faz_2026-08-23.md]]"
  - "[[../op-lattice-bath-runaway-2026-08-23/ANALIZA_N2_znak-W-z-akcji_2026-08-23.md]]"
---

# op-bath-two-sectors (PHASE0-LOCKED, zero obliczeń)

## Dwa pytania binarne

> **Q1:** czy mod runaway solitonu dostaje **ω² > 0** przy skończonej
> gęstości sąsiadów n — z fazami **zmierzonymi**, nie dokumentowanymi?
>
> **Q2:** czy znak tachionowy W **emerguje z gęstości** w akcji stabilnej —
> tj. czy rozdwojenie znaku z N2 to dwa sektory jednej akcji
> (izolowany obiekt ↔ kąpiel), a nie konflikt dwóch akcji?

## Struktura (kryteria: [[Phase0_balance.md]])

| Faza | Treść | Gate/werdykt |
|---|---|---|
| **1** | (A, δ) zmierzone per gatunek × {M-P, M-L}; drabina d\* par; kontrola Yukawa | R²≥0.999, R-stabilność, P1c czysta |
| **2** | **Q1:** baseline #63 V3 → skan 5 gęstości × 3 konfiguracje → ω²_min + ewolucja ± | Q1-PASS / Q1-FAIL / Q1-INCONCLUSIVE |
| **3** | **Q2 (niezależna):** tła ψ_n w sektorze stabilnym, d∈{∞,8,6,4,2} → znak ω²_min(d) | Q2-PASS / Q2-FAIL / Q2-INCONCLUSIVE |

Możliwe kombinacje wyników i ich znaczenie — drzewo decyzyjne LOCKa §6.
Szczególnie: **Q2-PASS** dałby derywacyjną kotwicę maszynerii 2 (koniec
dylematu „którą akcję wybrać"); **Q1-PASS** dałby pierwszy policzalny
mechanizm stabilności leptonów zgodny z ontologią.

## Wejście dla agenta

Gotowy prompt: [[HANDOFF_PROMPT.md]] (skopiować do nowej sesji w całości).
Kod do reuse: `N1_provenance_check.py` (układy+fit), `Phase3_nonlinear_evolution.py`
(#63, dynamika+gate), `RETRO_oscillating_tail_lock.py` (model drabiny).

## Status

- **Phase 0: 🟢 LOCKED (2026-08-23). Zero obliczeń.**
- Phase 1–3: oczekują (nowa sesja / osobny agent).
