---
title: "Phase 0 — pre-registration LOCK: czy Φ jest zmienną pomocniczą (non-radiative) w LIVE TGP?"
date: 2026-06-13
type: phase-balance
phase: 0
status: 🔒 LOCKED 2026-06-13 (before any computation)
cycle: op-phi-radiative-dof-audit-2026-06-13
scoping: "[[../../meta/SCOPING_op-phi-radiative-dof-audit_2026-06-13.md]]"
authorization: "User 2026-06-13: 'działaj z audytem'"
parent_verdict: "[[../op-PSR-Pdot-energy-balance-2026-06-13/]] PR-025 TRIGGERED (pending ratification)"
tags:
  - pre-registration
  - Dirac-constraint-analysis
  - fast-kill
  - structural-audit
---

# Phase 0 — pre-registration LOCK

## §1 — Pytanie

Czy w LOCKED inventory LIVE TGP istnieje struktura (wiąz, symetria cechowania,
modyfikacja kinetyki) czyniąca mod oddechowy δΦ NIE-radiacyjnym (auxiliary),
tak by P_total = κ_E·P_GR zamiast κ_E·P_GR + (1/6)P_GR?

Cykl czysto strukturalny: **zero dostępu do danych obserwacyjnych** (R_obs
pulsarów nietykane). Rozstrzyga wyłącznie status PR-025: „both branches"
vs „exhaustive over LIVE".

## §2 — Wejścia LOCKED

| Input | Wartość | Źródło |
|---|---|---|
| Zlinearyzowany L sektora Φ | L_Φ = (K₁/2)[(δΦ̇)² − (∇δΦ)²] − (m²/2)δΦ² − (q/Φ₀)ρδΦ, K₁ > 0 | Phase 5 `op-emergent-metric-from-interaction` §2 (LOCKED) |
| Sektor σ (Path A) | L_σ = −¼(∂σ_ab)² − ½m_σ²σ² − (ξ_eff/2)σ_abT^{ab,TT} | `op-sigma-3PN-radiative` Phase 3 §1.1 (LOCKED) |
| Symetrie LIVE działające na pola | Z₂ (dyskretna, Φ→−Φ); S05 U(1) (faza θ; moduł Φ inwariantny) | FOUNDATIONS §1, §5.1 |
| Wymóg Lorentza | substrat Lorentz-invariant (źródło α_i ≡ 0 strukturalnych) | FOUNDATIONS §3.6.2 L2 |
| Statyka | G̃(ω=0, k) = 1/k² (Newton 1/r; q² = 4πG z PR-025 T2a LOCKED) | PR-025 |

## §3 — Falsyfikatory LOCKED (reguły IMMUTABLE)

| ID | Test | Reguła decyzyjna |
|---|---|---|
| **F-AUX-A** | Analiza więzów Diraca (mode-wise Hessian + algorytm więzów) dla L_Φ; walidacja metody na EM (A₀ musi wyjść auxiliary) | Hessian osobliwy w sektorze δΦ ⟹ constraint istnieje ⟹ **DERIVED (Φ auxiliary)**; Hessian regularny ⟹ DOF_Φ = 1 ⟹ **NEGATIVE** |
| **F-AUX-B** | Czy istnieje Lorentz-invariantny propagator z G̃(0,k)=1/k² i bez bieguna radiacyjnego? | Wymuszenie F(s)=1/s na s>0 ⟹ biegun w s=0 ⟹ **NEGATIVE** (statyka⇔radiacja nierozdzielne); znaleziona kauzalna alternatywa ⟹ DERIVED-PATH (osobny PR) |
| **F-AUX-C** | Inwentarz symetrii: czy Z₂ lub S05 U(1) generuje wiąz pierwszej klasy na δΦ? | Z₂ dyskretna (brak generatora) ∧ U(1) działa tylko na θ (moduł inwariantny) ⟹ **NEGATIVE** |
| **F-AUX-D** | Agregat: werdykt WYLICZONY z flag A∧B∧C | wszystkie NEGATIVE ⟹ **HONEST_NEGATIVE** („Φ auxiliary" niewyprowadzalne w LIVE) ⟹ PR-025 upgrade **EXHAUSTIVE-OVER-LIVE**; jakiekolwiek DERIVED ⟹ PR-025 flagowany do re-run |

## §4 — Pre-derywacja (scoping §2, R1 #18): oczekiwane HONEST_NEGATIVE we wszystkich trzech. Rachunek nadrzędny nad oczekiwaniem.

## §5 — Forbidden moves

1. Zakaz modyfikacji L poza LOCKED inventory (modyfikacja = nowy PR, nie rescue).
2. Zakaz dostępu do danych obserwacyjnych.
3. Zakaz reinterpretacji PR-025 przed F-AUX-D.
4. Zakaz hardcoded T_pass; werdykt linii FINAL wyliczany z flag (lekcja GST #19).
5. Ghost-risk nie-FP masywnego tensora σ: poza zakresem — wolno tylko zarejestrować w DOUBTS.
6. Budżet nowych stałych: 0.
