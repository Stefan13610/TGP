---
title: "op-3d-canonical-lattice — CLOSED: Q-FAIL — sieć sc solitonów kanonicznych μ istnieje tylko przy d=2π i jest tachionowa (ω²_min=−1.674350, potwierdzone nieliniowo); 3D pogłębia niestabilność względem 1D (−1.22)"
date: 2026-08-31
type: research-cycle
status: CLOSED
tgp_owner: research/op-3d-canonical-lattice-2026-08-31
authorization: "User 2026-08-31: opcja 'b' z NEEDS op-3d-lattice-bath-stability (jawna zmiana bramki: kotwice #63 przestają obowiązywać w 3D; kotwica kanoniczna mierzona w cyklu) + realizacja"
related:
  - "[[Phase0_balance.md]]"
  - "[[Phase_FINAL_close.md]]"
  - "[[NEEDS.md]]"
  - "[[../op-3d-lattice-bath-stability-2026-08-31/Phase_FINAL_close.md]]"
---

# op-3d-canonical-lattice (2026-08-31)

**Status: CLOSED. Werdykt: Q-FAIL (PRIMARY po tłach istniejących;
ruling zapisany przed Phase 4).**

Re-lock pytania ω²(n) w 3D po GATE-FAIL-STOP poprzednika: cały cykl
w akcji kanonicznej (K=g⁴, U=g⁷/7−g⁸/8; spektra bez regularyzacji).
Most radialny→kartezjański PRZESZEDŁ (kotwica kanoniczna −1.646589
reprodukowana w 3D z odchyłem 0.33% przy h=0.3; diagnoza poprzednika
potwierdzona: bariera rozdzielczości była własnością ściany f_ε —
kieszeń kanoniczna V(r) jest szeroka: FWHM=1.20, w rdzeniu).

**Wynik centralny:** sieć sc istnieje wyłącznie przy d=2π
(π/d*₁: ucieczka g→0; 3π/4π: kolaps do próżni) i jest TACHIONOWA:
ω²_min(2π) = −1.674350 (N=48, zbieżnie, po odjęciu 3 modów
translacyjnych; argmin Γ; per k: Γ/X/M/R = −1.674/−1.664/−1.654/−1.645).
Nieliniowo (superkomórka 2×2×2, oba tła, 8 biegów): ucieczka
t_esc = 3.98–4.20 ≤ 2·t*_izo = 9.42 przy obu znakach; gate energii
≤2.7e−8. Kontrole: P3b PASS 10/10 (po korekcie 1: krotność 8 vs k=10
w ARPACK — udokumentowana, poprawiona przed użyciem), P3c PASS (+1.065).
Trend vs 1D: −1.674 < −1.222 — trzeci wymiar POGŁĘBIA niestabilność;
hipoteza stabilizacji gęstością zamknięta także w 3D sc.

Kryteria: [[Phase0_balance.md]] (LOCK); wynik: [[Phase_FINAL_close.md]];
dalsze kroki (user-gated): [[NEEDS.md]].

## Log faz

- 2026-08-31: Phase 0 LOCK zapisany. Obliczenia: NIEROZPOCZĘTE.
- 2026-08-31: Phase_method_decisions.md FROZEN (delta vs poprzednik;
  g₀_μ zamrożone formułą φ·0.90548=1.4650974).
- 2026-08-31: Phase 1 PASS — kotwica kanoniczna λ_min(w1)=−1.646589
  (tabela zbieżności h∈{0.05,0.025,0.0125}); kieszeń V(r): FWHM=1.20
  (anty-pułapka OK); t*_ref=4.336; most: −1.652 przy h=0.3 (0.33%);
  t*_izo(3D)=4.710. P1a/P1d dziedziczone z cytatem.
- 2026-08-31: Phase 2 — istnienie TYLKO d=2π (oba starty; ‖R‖≤6.7e−10,
  dgrid≤4.6e−3). Phase2_backgrounds3d.npz zapisany.
- 2026-08-31: Phase 3 — ω²_min(2π)=−1.674350 (ZBIEŻNE), translacje
  zidentyfikowane (λ~O(h²), coverage 1.0); P3c PASS; P3b FAIL w punkcie
  d=4π N=48 → diagnostyka FD-exact + correction note (krotność 8 vs
  k=10) → po korekcie PASS 10/10.
- 2026-08-31: Phase3_verdict_ruling.md zapisany PRZED Phase 4.
- 2026-08-31: Phase 4 — 8/8 biegów UCIECZKA ≤2t*; incydent: limit czasu
  tła ubił proces po komplecie A0.7 → kontynuacja A1.0 w 4 osobnych
  procesach (Phase4_continue_A10.py, delta czysto implementacyjna).
- 2026-08-31: **Q-FAIL**; Phase_FINAL_close.md + NEEDS.md zapisane.
  Cykl CLOSED.
