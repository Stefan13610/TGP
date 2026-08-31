---
title: "op-fluctuation-extended-nbody — obiekty rozciągłe (czy −1/d?) + N-ciałowość kanału fluktuacyjnego, poziom 0"
type: research_cycle
status: CLOSED-EXECUTED
phase: "FINAL (LOCK → Phase 1–3 → close, jedna sesja)"
folder_status: closed
claim_status: "CLOSED-EXECUTED 2026-08-31. QE (N1 rodzica): NIE — wykładnik dalekiego pola kanału fluktuacyjnego pozostaje −2 dla kul R∈{1,2,3} (slopes −2.07/−2.09/−2.15, R²≥0.9999); przy kontakcie STROMIEJE (do ~−3), nigdy nie łagodnieje do −1 (najdłuższy przebieg p_loc=−1±0.15: 0); R wchodzi tylko w amplitudę wg obrazu pojemnościowego A~½C_R²(4π)⁻² (zgodność 93–96%); uniwersalność znaku i zasięg 2μ przeżywają rozciągłość. QN (N2 rodzica): TAK — ΔF₃ = g₁₂g₁₃g₂₃ + O(g⁴) (sympy exact), znak DODATNI (osłabia przyciąganie parowe) we wszystkich 20 punktach; |ΔF₃|/|ΣF_par| = 1.1–4.0% na krytyczności, zanik ~d⁻³; kanał klasyczny addytywny dokładnie (<1e−12). Poziom 0; zero claimów obserwacyjnych."
created_date: 2026-08-31
closed_date: 2026-08-31
authorization: "User 2026-08-31: „zacomituj i zajmij się research/op-substrate-fluctuation-channel-2026-08-23/" → realizacja NEEDS N1+N2 rodzica jako cykl-następca"
anti_lakatos_lock: ACTIVE
related:
  - "[[Phase0_balance.md]]"
  - "[[Phase_FINAL_close.md]]"
  - "[[../op-substrate-fluctuation-channel-2026-08-23/README.md]]"
  - "[[../op-substrate-fluctuation-channel-2026-08-23/NEEDS.md]]"
---

# op-fluctuation-extended-nbody (2026-08-31)

Cykl-następca `op-substrate-fluctuation-channel` (poziom 0): wykonanie jego
research-NEEDS **N1** (obiekty rozciągłe — czy istnieje reżim Newtonowskiego
−1/d) i **N2** (N-ciałowość log-det). Maszyneria dziedziczona: siatka L³,
G przez FFT, pinning, Amendment A1 (propagator connected na krytyczności).

Pytania, modele, kryteria, okna i forbidden moves: [[Phase0_balance.md]]
(LOCK zamknięty przed kodem).

## Status
- Phase 0: 🟢 LOCKED (2026-08-31).
- Phase 1: 🟢 5/5 PASS (sympy: ΔF₃ do rzędu 4; tożsamość monopolowa
  det = 1−c²C_AC_B; wykładnik −2 niezależnie od pojemności).
  Incydent P1-5a (bug składniowy, crash nie FAIL) — udokumentowany w close.
- Phase 2: 🟢 9/9 PASS → **QE: NIE** (wykładnik zostaje −2; R tylko
  w amplitudzie; profile p_loc bez selekcji w outputach).
- Phase 3: 🟢 3/3 PASS → **QN: TAK** (nieaddytywny, znak dodatni
  uniwersalny, kontrast z dokładnie addytywnym kanałem klasycznym).
- Zamknięcie: [[Phase_FINAL_close.md]] · NEEDS N1–N4 user-gated ([[NEEDS.md]]).
