---
title: "Phase2_correction_note — błąd implementacyjny WYŁĄCZNIE deskryptywnej diagnostyki „symetria czworek ±λ" (suma zamiast różnicy posortowanych widm); zero wpływu na bramki Gate A/C1/C2/C3 i ich progi"
date: 2026-08-31
type: correction-note
tgp_owner: research/op-symplectic-Jspectrum-2026-08-31
status: APPLIED-BEFORE-USE
related:
  - "[[Phase0_balance.md]]"
  - "[[Phase_method_decisions.md]]"
  - "[[Phase2_gate_nls.py]]"
---

# Correction note — Phase 2 (zapisana PRZED użyciem poprawionego wyniku)

Pierwszy bieg `Phase2_gate_nls.py` (output pierwotny ZACHOWANY:
`Phase2_output_run1_original.txt`): WSZYSTKIE bramki zalockowane PASS
(Gate A kotwica |δ|=1.19e−7; C1 max Re λ=1.5e−7 ≤ tol; C2 λ=2.915
zbieżnie; C3 maxerr=6.8e−4, ratio 4.00). Błąd dotyczy wyłącznie
pomocniczej, NIE-kryterialnej linii diagnostycznej „symetria czworek
±λ": wzór `max|sort(λ) + sort(−λ)|` — dla widma symetrycznego multiset
{−λ} = {λ}, więc wyrażenie redukuje się do `2·max|λ|` (stąd wartości
5.2e3/2.1e4 ≈ 2·|λ|_UV) i niczego nie testuje. Poprawny test symetrii:
`max_i min_j |λ_i + λ_j|` (odległość widma od jego odbicia λ→−λ).

Żadne kryterium/próg/bramka nie ulega zmianie; korekta wyłącznie tej
linii diagnostycznej + rerun (output poprawiony: `Phase2_output.txt`;
pierwotny zachowany obok). Ta sama poprawna forma diagnostyki wchodzi
do Phase 3.
