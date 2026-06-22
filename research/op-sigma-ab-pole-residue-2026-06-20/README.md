---
title: "op-sigma-ab-pole-residue — czy framework dostarcza warunek pole-residue ustalający C_σ (κ_E) jako predykcję? WERDYKT: NIE (decydujący negatyw)"
date: 2026-06-20
type: cycle_readme
status: 🟢 CLOSED (2026-06-20) — WERDYKT NEGATYWNY: brak pole-residue; C_σ pozostaje wolnym parametrem; κ_E = wolny parametr (opcja b)
cycle: op-sigma-ab-pole-residue-2026-06-20
parent_cycles:
  - "[[../op-CG4-substrate-closure-2026-06-20/Phase3_renorm.md]] (C_σ UV-czuły = wolny parametr)"
  - "[[../closure_2026-04-26/sigma_ab_pathB/results.md]] (Path B; Bethe-Salpeter OTWARTE §5)"
anti_lakatos_lock: ACTIVATED
authorization: "User 2026-06-20: 'działaj z a osobny cykl zobaczmy co z tego wyjdzie'"
tags: [sigma-ab, pole-residue, bound-state, Bethe-Salpeter, spin-2, kappa-E, free-parameter, anti-Lakatos]
---

# op-sigma-ab-pole-residue

> **Pytanie:** op-CG4 Phase 3 dowiódł, że $C_\sigma$ jest UV-czuły = wolny parametr. Jedyna droga do
> predykcji $\kappa_E$: izolowany biegun (bound state) σ_ab z residuum on-shell ustalającym $C_\sigma$.
> **Czy framework go dostarcza?**
>
> **WERDYKT (value-blind, wyliczony): NIE.** σ_ab to **kontinuum** (cięcie bąbla atan, próg $4m_s^2$),
> nie izolowany biegun; kontakt φ⁴ (substrat M0) ma **zerową projekcję na falę spin-2** (L=2) ⟹ brak
> bound-state ⟹ brak residuum ⟹ **$C_\sigma$ pozostaje wolny**. „$M^2{=}2m_s^2$" to coeff OPE, NIE masa bieguna.

## Wynik (Phase 1)

| Kryterium | Werdykt | Podstawa (sympy) |
|---|---|---|
| **C-POLE** izolowany biegun | ❌ FAIL | $\Pi(p)=\arctan(p/2m)/(4\pi p)$ — cięcie od $p^2=-4m^2$, brak prostego bieguna |
| **C-KERNEL** jądro w spin-2 | ❌ FAIL | kontakt: $\int_{-1}^1 P_2\,dx=0$; $(k\cdot k')\sim x$: $\int xP_2=0$; tylko $\ge2$-pochodne ($x^2$): $4/15\neq0$ |
| **C-RESIDUE** residuum ustala $C_\sigma$ | ❌ FAIL | brak bieguna ⟹ brak residuum on-shell |
| **C-MATCH** biegun przy $2m_s^2$? | USTALENIE | $2m_s^2$ = coeff OPE/heredity, nie pozycja bieguna |

**Reguła agregatu (Phase0_lock §3):** C-KERNEL FAIL → C-POLE FAIL → C-RESIDUE FAIL ⟹ **framework NIE dostarcza
pole-residue**; $\kappa_E$ = **genuine wolny parametr** sektora radiacyjnego ⟹ **opcja (b)** jest uczciwym domknięciem.

## Znaczenie

Cykl (a) **zamyka** pytanie postawione po op-CG4: nie ma ukrytego warunku normalizacyjnego, który zamieniłby
$\kappa_E$ w predykcję. Fizyczny powód jest czysty i ogólny: **interakcja kontaktowa (s-wave) nie wiąże w kanale
d-wave (spin-2)**, więc tensorowy nośnik GW nie ma izolowanego stanu związanego — tylko kontinuum. To czyni
$\kappa_E$ jawnie wolnym parametrem EFT sektora radiacyjnego (potwierdza i domyka op-CG4 Phase 3).

## Residual (jawny)
Interakcja $\ge2$-pochodne (kwadrupolowa, np. wariant gradient-bond v2) MA niezerowe L=2 — ALE: (i) to operator
**irrelevant** (wyższy wymiar), nie część zwalidowanego substratu M0; (ii) nawet gdyby obecny, wiązanie
niegwarantowane, a pozycja bieguna BS **nie** równa się generycznie $2m_s^2$ ⟹ nie ratuje predykcji czysto.
Negatyw jest **robust** dla ustalonego frameworku.

## Pliki
- [[Phase0_lock.md]] — kryteria, reguła agregatu, forbidden moves
- [[Phase0_balance.md]] — gate
- [[Phase1_pole.py]] / [[Phase1_output.txt]] — rachunek (spektrum + projekcja L=2 + warunek BS)
- [[Phase_FINAL_close.md]] — werdykt
