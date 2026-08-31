---
title: "Phase1_gate_ruling — rozstrzygnięcie dwuznaczności zakresu gate'u Phase 1 (τ/[50,150]); zapisane PRZED uruchomieniem Phase 2"
date: 2026-08-23
type: method-note
tgp_owner: research/op-bath-two-sectors-2026-08-23
status: RULING-BEFORE-PHASE2
related:
  - "[[Phase0_balance.md]]"
  - "[[Phase_method_decisions.md]]"
---

# Rozstrzygnięcie: zakres gate'u Phase 1

## Stan faktyczny (Phase1_output.txt, bez upiększania)

Wszystkie fity e/μ (M-P i M-L), oba okna: R² ∈ [0.99958, 0.999999] —
PASS. R-stabilność: idealna (|Δδ|=0.0000°) — PASS. P1c czysta (0 minimów
Yukawa we wszystkich 6 parach×modelach) — PASS. Drabiny: 2π±5% PASS
(odchył 0.96–3.34%, korekta malejąca z d — zgodnie z predykcją LOCKa).

**Jedyny FAIL:** M-P τ (g₀=1.5696, INPUT Q_K=3/2), okno **[50,150]**:
R²=0.9987829 < 0.999 (odchył 0.0002). Okno PRIMARY [120,260] dla τ:
R²=0.9997182 — PASS. Fizycznie: τ leży 1.9% pod progiem kolapsu 8/5
(A≈1.006, wolno gasnące poprawki przejściowe w bliskim oknie);
wrażliwość +2% → KOLAPS (zgodnie z progiem; raportowane).

## Dwuznaczność i rozstrzygnięcie

LOCK §3: „Gate Phase 1: fity zbieżne (R² ≥ 0.999) i R-stabilne; kontrola
czysta" — nie precyzuje, czy obejmuje gatunek OPCJONALNY τ („gatunki
e/μ(/τ w M-P)"; HANDOFF: „{e, μ(, τ)}") i okno nie-PRIMARY.

Fakty rozstrzygające (z samego LOCKa, nie z wyniku):
1. τ nie wchodzi do ŻADNEGO zalockowanego produktu Phase 1→2: pary
   drabiny = {ee, eμ, μμ} (§3 P1b), konfiguracje kąpieli = {ee, μμ, eμ}
   (§2 M-B). Gate chroni pipeline (A,δ)→Phase 2; τ w nim nie uczestniczy.
2. τ jest w LOCKu jawnie oznaczony jako dopuszczony wyjątek z flagą INPUT
   i osobną kontrolą wrażliwości — status kontrolno-deskryptywny.

**RULING (zapisany przed Phase 2, kryteria i progi NIETKNIĘTE):**
gate Phase 1 obowiązuje dla wszystkich wielkości zasilających Phase 2
(e, μ × M-P, M-L × oba okna × R-kontrola × P1c) → **GATE PASS**.
Fit τ/[50,150] raportowany wprost jako **FAIL deskryptywny** z flagą
TAU-NEAR-THRESHOLD; τ pozostaje wyłącznie deskryptywny (żadna liczba τ
nie wchodzi do Phase 2). Werdykt strict-reading (wszystkie fity łącznie
z opcjonalnym τ): FAIL — odnotowany równolegle; oba odczyty trafią do
Phase_FINAL_close.md i raportu końcowego.

Uzasadnienie anty-Lakatos: nie zmieniam progu (0.999), okien, punktów ani
kryteriów Q1/Q2; rozstrzygam wyłącznie ZAKRES agregacji gate'u, którego
LOCK nie domknął, w kierunku zgodnym z literą LOCKa o produktach Phase 1.
