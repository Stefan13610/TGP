---
title: "op-bath-two-sectors — runaway w kąpieli na fazach zmierzonych (Q1) + derywacyjny test dwóch sektorów znaku W (Q2)"
type: research_cycle
status: CLOSED
phase: "FINAL (Phase 1–3 wykonane 2026-08-23; zamknięcie 2026-08-29)"
folder_status: closed
claim_status: "CLOSED 2026-08-29 — Q1-INCONCLUSIVE, Q2-FAIL. Phase 1: gate PASS wg rulingu (Phase1_gate_ruling.md; strict-reading z opcjonalnym τ: FAIL); (A,δ) zmierzone [120,260]: M-P e(0.0936,−75.34°) μ(0.6159,+97.20°), M-L e(0.0959,−81.43°) μ(0.3637,+38.58°) — Δ_ML(e→μ)=120.01° reprodukuje N1; drabiny 2π±5% PASS 6/6; P1c czysta. Q1-INCONCLUSIVE: P2a baseline PASS (λ_min=−1.3896 vs #63 −1.389, t*=3.62, gate ≤1.6e−8), ale ω²_min NIEZBIEŻNE 0/15 punktów (znak podąża za rozmiarem komórki R≶π); komparator BEZ kąpieli daje niemal identyczne spektra i breakdowny — wariant komórkowy z obciętym tłem niezdolny rozstrzygnąć Q1 (N3 pozostaje OTWARTE). Q2-FAIL: kontrola d=∞ PASS (m²=γ±0.00%), ω²_min(d)=+1.34/+1.57/+1.88/+2.47 dla d=8/6/4/2 — zbieżne, dodatnie, ROSNĄCE z gęstością; odpowiedź monotoniczna. Hipoteza dwóch sektorów OBALONA w klasie zbadanej; znak W = otwarty problem AKSJOMATYCZNY (decyzja ontologiczna autora). NEEDS user-gated."
created_date: 2026-08-23
closed_date: 2026-08-29
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

- Phase 0: 🟢 LOCKED (2026-08-23).
- Phase 1: ✅ wykonana (2026-08-23) — gate PASS wg rulingu
  ([[Phase1_gate_ruling.md]]); Δ_ML(e→μ)=120.01°, drabiny 2π±5% PASS 6/6.
- Phase 2 (Q1): ⚠️ **Q1-INCONCLUSIVE** — niezbieżność siatka×komórka 0/15;
  komparator bez kąpieli pada identycznie (artefakt komórki dominuje).
- Phase 3 (Q2): ❌ **Q2-FAIL** — ω²_min(d) dodatnie i rosnące z gęstością;
  hipoteza dwóch sektorów obalona w klasie zbadanej.
- **CLOSED 2026-08-29.** Werdykty i liczby: [[Phase_FINAL_close.md]];
  propozycje user-gated: [[NEEDS.md]].
