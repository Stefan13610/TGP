---
title: "Phase_correction_note — korekta 1: obsługa załamania niefinitycznego w silniku gradient flow (Phase 2); zero zmian kryteriów/progów/detektora"
date: 2026-09-01
type: correction-note
tgp_owner: research/op-metametric-boundary-2026-09-01
status: DOCUMENTED-BEFORE-USE
related:
  - "[[Phase0_balance.md]]"
  - "[[Phase_method_decisions.md]]"
  - "[[Phase2_floor_relax.py]]"
---

# Korekta 1 (Phase 2): przechwycenie załamania niefinitycznego

**Objaw (pierwotny output zachowany:
`Phase2_results/job_tach_sol_f2_h025_pre_correction.txt`):** pierwszy bieg
P2b `tach_sol_f2_h025` zakończył się nieobsłużonym wyjątkiem
`ValueError: array must not contain infs or NaNs` (solve_banded na
niefinitycznym rhs) w kroku 275 (t ≈ 2.75), zamiast klasyfikacją biegu.

**Diagnoza (przed użyciem jakiegokolwiek wyniku P2b):** przyczyna
FIZYCZNA to ucieczka g → +∞ rdzenia solitonu (mod niestabilny siodła,
wzrost ~e^{1.65t}; U(g) = g⁷/7 − g⁸/8 → −∞ dla g → ∞ — znana własność
zalockowanego modelu, odnotowana jawnie w MD op-3d-canonical §1:
„runaway kanoniczny może iść w g→duże, U nieograniczone z dołu").
Trajektoria zweryfikowana diagnostycznie: monotoniczny wzrost max g
w r≈0 od t≈1, blowup skończono-czasowy t≈2.7. Podłoga substratowa
(g_floor < 1) z konstrukcji LOCKa ogranicza pole WYŁĄCZNIE od dołu —
nie dotyka tego kanału.

**Błąd implementacji (znaleziony):** silnik `run_flow` sprawdzał
finityczność stanu PO kroku, ale sam krok (solve_banded/FFT) dostawał
już niefinityczny rhs i rzucał wyjątek — załamanie nie było
klasyfikowane zgodnie z literą LOCKa („załamanie NIE-nukleacyjne (NaN…)
= INCONCLUSIVE, nie pozytyw").

**Korekta (wyłącznie obsługa błędu):**
1. kontrola finityczności rhs w `step` + try/except wokół kroku;
2. klasyfikacja biegu `BREAKDOWN` z czasem t_break = k·dt;
3. zachowanie OSTATNIEGO finitycznego stanu (g_prev) jako stanu
   raportowanego w npz + diagnostyka kierunku załamania
   (min/max g przed załamaniem);
4. `np.seterr` na ignore dla over/invalid w silniku (NaN propaguje do
   detekcji zamiast crashować proces).

Kryteria, progi, detektor, seed, siatki, schemat flow — NIETKNIĘTE.
Diagnostyka trajektorii wykonana pomocniczo (wydruk min/max g) nie
zasila żadnej bramki. Pierwotny output FAIL niezamazany.
