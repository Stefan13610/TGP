---
title: "op-lattice-bath-runaway — audyt maszynerii 2 (bramka) + test runaway solitonu w kąpieli sąsiadów"
type: research_cycle
status: CLOSED-GATE-STOP
phase: "A (bramka) — FAIL na A3; Phase 1–3 NIEURUCHOMIONE"
folder_status: closed
claim_status: "CLOSED-GATE-STOP 2026-08-23 — bramka A: A1 PASS (g0_crit=8/5, |Δ|=3.7e−11; formuła na 4 parach (α,d)), A2 PASS (ω_tail=1.00000 dla α∈{1,2,3}, kontrola negatywna czysta), A3 FAIL MERYTORYCZNY (fazy ogona z dodatekH lin. 1126–1129 niereprodukowalne: δ_e=−75.50° vs −81.4±2°, δ_μ=+88.48° vs +38.6±2°, Δ(e→μ)=163.98° vs 120±1°; wykluczone konwencja/okno/wariant układu w 3 iteracjach diagnostycznych; A_μ dodatekH sprzeczne z własnymi skryptami rdzenia), A5 ROZSTRZYGNIĘTE (różnica = znak źródła W→−W; dwa różne układy), A4/A6 wykonane. Reguła bramki → STOP całego cyklu przed Phase 1. Rachunek centralny (ω²(n) w kąpieli) NIEPOLICZONY. Podtest amplitudy A≈|g0−1| PASS (0.993–1.011). Q_K=3/2 = INPUT (flagowane)."
created_date: 2026-08-23
closed_date: 2026-08-23
authorization: "User 2026-08-23: 'tak zrób też audyt i start z fazą 0' — realizacja Phase A: agent-implementator 2026-08-23"
anti_lakatos_lock: PRESERVED
related:
  - "[[Phase0_balance.md]]"
  - "[[PhaseA_report.md]]"
  - "[[Phase_FINAL_close.md]]"
  - "[[NEEDS.md]]"
  - "[[../op-native-pressure-lepton-stability-2026-07-27/README.md]]"
  - "[[../op-native-pressure-lepton-stability-2026-07-27/ANALIZA_retrospektywa_oscylacyjny-lock_2026-08-16.md]]"
  - "[[../op-nonlinear-charge-constraint-2026-07-03/README.md]]"
---

# op-lattice-bath-runaway (CLOSED-GATE-STOP 2026-08-23)

## Pytanie binarne

> Czy mod runaway solitonu dostaje **ω² > 0** przy jakiejś skończonej gęstości
> sąsiadów n — w konfiguracji, którą ontologia TGP uznaje za fizyczną,
> a która nigdy w korpusie nie została policzona?

**Status pytania: NIEROZSTRZYGNIĘTE** — cykl zatrzymany na bramce
(audyt maszynerii 2), zanim rachunek centralny został uruchomiony.

## Wynik faz

| Faza | Treść | Werdykt |
|---|---|---|
| **A** | Audyt maszynerii 2 (A1–A6) | **A1 PASS · A2 PASS · A3 FAIL · A5 rozstrzygnięte · A4/A6 wykonane → BRAMKA ZAMKNIĘTA (STOP)** |
| **1** | (κ, φ, A) per gatunek; d\* par | NIEURUCHOMIONA (reguła bramki) |
| **2** | RACHUNEK CENTRALNY (baseline #63 V3 → skan n → ω²_min) | NIEURUCHOMIONA |
| **3** | ω²(n) deskryptywnie | nie dotyczy |

Sedno A3 (pełna diagnoza: [[PhaseA_report.md]] §3): fazy ogona per
gatunek — dokładnie to wejście, które Phase 1→2 miały brać z rdzenia —
są w dodatekH niereprodukowalne z udokumentowanego ODE (Δ(e→μ)=2π/3
nie odtwarza się: zmierzono 163.98°), przy czym amplitudy A_tail
z własnych skryptów rdzenia zgadzają się z niezależną re-implementacją
do 5 cyfr (błąd implementacji wykluczony).

Co przeżyło audyt: próg kolapsu 8/5 (1e−10), ω=1 uniwersalne,
oscylacyjność ogona, A≈|g₀−1| (±1.1%). Pada warstwa FAZOWA.

## Dokumenty

- [[Phase0_balance.md]] — LOCK (nienaruszony).
- [[PhaseA_report.md]] — pełny raport bramki + tabela A4 (8 skryptów).
- [[Phase_FINAL_close.md]] — zamknięcie wg drzewa decyzyjnego LOCKa §6.
- [[NEEDS.md]] — user-gated: N1 pochodzenie faz w dodatekH, N2 znak W
  z akcji, N3 ewentualny nowy cykl na fazach zmierzonych, N4 rysa τ(g₀=4).
