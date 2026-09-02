---
title: "op-metric-closure-relaxation — relaksacja z obustronnym domknięciem: podłoga QB-2 + granica metryczna ψ=4/3 (M9.1''); test stanu metametrycznego i kaskady nukleacji"
date: 2026-09-02
type: research-cycle
status: CLOSED
tgp_owner: research/op-metric-closure-relaxation-2026-09-02
authorization: "User 2026-09-02: „ok, działaj" (sekwencja N1→re-lock Q2 z NEEDS op-metametric-boundary)"
related:
  - "[[Phase0_balance.md]]"
  - "[[../op-metametric-boundary-2026-09-01/Phase_FINAL_close.md]]"
  - "[[../op-blocked-soliton-bang-2026-07-04/README.md]]"
---

# op-metric-closure-relaxation (2026-09-02)

**Status: CLOSED (2026-09-02). Werdykt: Q-PASS-NUCLEATION (litera) —
nukleacja dolna zbieżna w parze sol×C-BAR (N_det=1±0); ROZJAZD
PRIMARY↔C-BAR totalny: wariant metryczny w(ψ) załamuje się 10/10,
bo biegun ψ=4/3 PRZYCIĄGA kanoniczne U (U−U(1)<0 na granicy) —
szczegóły [[Phase_FINAL_close.md]].**

Następca op-metametric-boundary (Q1-POS / Q2-INCONCLUSIVE przez nowy
kanał ucieczki g→+∞). **N1 rozstrzygnięte dokumentacyjnie:** korpus MA
naturalne górne domknięcie — biegun metryki efektywnej M9.1''
(√−g_eff=c₀ψ/(4−3ψ), g_tt→0 przy ψ=4/3) ⟹ **g_ceil=√(4/3)=1.15470** —
dosłownie „granica metryki" z hipotezy autora. Obserwacja: tła bloch
miały g_max=1.141–1.143, tuż pod granicą.

**Pytanie:** czy z obustronnym domknięciem (podłoga QB-2 + czynnik
metryczny w akcji, wariant kontrolny z gołą barierą dla separacji
efektów) i większym pudłem (L=4π — pasmo tachionowe próżni częściowo
w pudle) relaksacja osiąga stan metametryczny (statyczny strukturalny)
albo kaskadę nukleacji (pre-rejestrowany POZYTYW; detektor dziedziczony
+ nowy detektor obiektów górnych przy granicy metrycznej)?

Kryteria: [[Phase0_balance.md]].

## Log faz

- 2026-09-02: Phase 0 LOCK zapisany. Obliczenia: NIEROZPOCZĘTE.
- 2026-09-02: `Phase_method_decisions.md` FROZEN (wariacja PRIMARY
  z członem w′, referencja próżniowa U−U(1) wymuszona literą P1a,
  regularizacja bieguna tylko >g_ceil−1e−6, detektor górny 1.0773503).
- 2026-09-02: Phase 1 — P1a 6/6 PASS (dryf 0.0), P1b 2/2 PASS
  (BREAKDOWN t=2.750/3.130 = reprodukcja poprzednika), P1c 4/4 PASS
  po korekcie macierzy kontroli
  ([[Phase_correction_note_p1c_matrix.md]]; tło 3D ma g_max=1.4734 >
  g_ceil — lat×PRIM poza domeną kontroli; pierwotny output zachowany).
- 2026-09-02: Phase 2 — 14 biegów + 2 dt/2: PRIMARY 10/10 BREAKDOWN
  (geneza t=8.72/8.77 identycznie dla 3 podłóg; sol/lat t≤0.06 —
  starty poza domeną); sol×C-BAR NUCLEATION-DN zbieżna (t₀=2.0,
  N_det=1 w 4/4); lat×C-BAR STATIONARY jednorodne g≡0.5354 (podłoga).
  **Werdykt: Q-PASS-NUCLEATION.**
- 2026-09-02: Phase 3 widmowe NIE WYKONANE (warunek Q-PASS-STATIC
  niespełniony); charakterystyka kaskady bez progów →
  `Phase3_output.txt` (pojedyncza kreacja obiektu podprogowego —
  inwersja rdzenia na podłogę QB-2; tło w studni bariery
  g_ceil+0.0994).
- 2026-09-02: `Phase_FINAL_close.md`, `NEEDS.md` (user-gated) —
  cykl ZAMKNIĘTY.
