---
title: "op-dynamics-class-M911 — klasa dynamiki pary metrycznej: zachowanie substratu (Cahn–Hilliard) i bezwładność z lapse M9.1'' vs quietyzm relaksacyjny; strefa słabego pola ψ<2/3"
folder_status: closed
date: 2026-09-02
type: research-cycle
status: CLOSED
verdict: "Phase 1 PASS komplet; Q-CONS-FAIL (dyfuzja 6/6 do jednorodnej średniej, ~2–3× wolniej — zachowanie substratu nie ratuje struktury); Q-INER-INCONCLUSIVE wg litery przy zbieżnych kategoriach: geneza DYSPERSJA, sieć→GÓRNA granica 4/3 (czasy identyczne siatka×dt), dip słabopolowy→kolaps do DOLNEJ granicy przedmetrycznej ψ→0 (t≈5.55 zbieżnie; finał STIFF pre-rejestrowany). Klasa bezwładna jest granico-szukająca — nie relaksuje."
tgp_owner: research/op-dynamics-class-M911-2026-09-02
authorization: "User 2026-09-02: „Ok A + B działaj ;)" (rekomendacja A+B po zamknięciu op-metric-pair-M911; hipoteza usera: kreacja w przedmetrycznym stanie zaburzonego substratu przy słabym ψ na granicy rozlewającego się solitonu)"
related:
  - "[[Phase0_balance.md]]"
  - "[[../op-metric-pair-M911-2026-09-02/Phase_FINAL_close.md]]"
  - "[[../op-metric-pair-M911-2026-09-02/NEEDS.md]]"
---

# op-dynamics-class-M911 (2026-09-02)

Następca op-metric-pair-M911 (Q-A-PASS: para samodomknięta;
Q-B-FAIL: relaksacja gradient flow spływa do próżni w t≈7–16).
Pytanie tego cyklu: czy szybki quietyzm był własnością KLASY dynamiki
(przetłumiona, niezachowana, bez lapse), a nie pary (w, V_M9.1'')?

**Warianty:** A — dynamika ZACHOWANA (Cahn–Hilliard, ∫ψ=const, M≡1);
B — dynamika 2. rzędu BEZWŁADNA z wagą czasową metryki M9.1''
(B(ψ)=w²K=ψ⁶/(4−3ψ)², H zachowane). Starty: geneza z szumu
(seed=20260904), NOWY dip słabopolowy ψ_min=0.55<2/3 (strefa
spinodalna 𝒰″<0 — hipoteza usera), sieć 2π przeskalowana (1.30).
Detektory i progi dziedziczone (dolny 5/6, górny 7/6).

**Pytania:** Q-CONS (A) — nukleacja / trwała struktura / dyfuzja do
średniej + ILOKROTNE spowolnienie vs poprzednik; Q-INER (B) —
nukleacja / oscylony (okno [80,100], dev≥0.05, obłożenie ≥80%) /
dyspersja. Kryteria: [[Phase0_balance.md]].

## Log faz

- 2026-09-02: Phase 0 LOCK zapisany (autoryzacja usera „Ok A + B
  działaj"). Implementator = ta sama sesja.
- 2026-09-02: `Phase_method_decisions.md` FROZEN przed kodem
  (B=w²K z |g^tt| M9.1''; analiza stabilności dt_B; ψ_stiff=0.12
  pre-rejestrowane). `M911_common.py` — wspólny rdzeń.
- 2026-09-02: Phase 1 — **PASS komplet** (G1: dryf 0.0, masa 1.1e−16,
  E monotone; G2: ω vs dyspersja dyskretna 1.4e−11, dryf H 2.9e−5,
  B′ 1.1e−14; G3: detektory 3/3). Incydent pre-compute bez wpływu:
  SyntaxError (em-dash w źródle ASCII) poprawiony przed pierwszym
  przebiegiem bramek. `Phase1_output.txt`.
- 2026-09-02: Phase 2 (wariant A) — **Q-CONS-FAIL** (6/6 STATIONARY
  w średniej: gen 0.999999 t=22, dip 0.987947 t=49, sieć 0.935 t=12;
  masa zachowana; zdarzeń zero ⟹ dt/2 niewymagane). `Phase2_output.txt`,
  `Phase2_relaxed_states.npz`.
- 2026-09-02: Phase 3 (wariant B) — **Q-INER-INCONCLUSIVE** wg litery;
  kategorie zbieżne: geneza TMAX-DYSPERSJA (dev~6e−4, ε_H≤1.5e−10);
  sieć BREAKDOWN-BOUNDARY górna (4.6450/4.6450 N32, 4.1100/4.1087 N48 —
  dt vs dt/2); dip BREAKDOWN-BOUNDARY-LOWER (5.5725/5.5500, dt i dt/2
  identycznie; finał STIFF pre-rejestrowany). 6 biegów + 4 dt/2.
  `Phase3_output.txt`, `Phase3_final_states.npz`.
- 2026-09-02: Zamknięcie: `Phase_FINAL_close.md` + `NEEDS.md`
  (user-gated: N1 fizyka ψ→0 poziomu 0/UV; N2 materia na dynamice
  bezwładnej; N3 re-lock techniczny; N4 dopisek core; N5 status
  hipotezy słabopolowej). npz tła READ-ONLY niezmieniony. **CLOSED.**
