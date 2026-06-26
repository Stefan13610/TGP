---
title: "Phase 1 results — Dirac DOF audit: Φ-auxiliary niewyprowadzalne w LIVE"
date: 2026-06-13
type: phase-results
phase: 1
status: 🟢 EXECUTED — 13/13 sympy PASS — F-AUX-D = HONEST_NEGATIVE (computed from flags)
parent: "[[./README.md]]"
predecessor: "[[./Phase0_balance.md]] (LOCKED before computation)"
sympy_script: "[[./Phase1_sympy.py]]"
sympy_output: "[[./Phase1_sympy.txt]]"
tags:
  - Dirac-constraint-analysis
  - radiative-dof
  - HONEST-NEGATIVE
---

# Phase 1 — wyniki

## §0 — Executive summary

**13/13 PASS, 0 hardcoded, werdykt WYLICZONY z flag (D.1).**

| Falsyfikator | Flaga | Treść |
|---|---|---|
| F-AUX-A | **NEGATIVE** | Hessian kinetyczny sektora Φ = K₁ > 0 regularny ⟹ zero więzów pierwotnych ⟹ DOF_Φ = 1 (propagujący). Walidacja metody na EM: A₀ poprawnie wykryte jako auxiliary (rank 2/3, nullspace = A₀). Struktura więzów k-niezależna ⟹ żaden multipol (w tym kwadrupol) nie jest wyróżniony. |
| F-AUX-B | **NEGATIVE** | Lorentz-invariancja + lock statyczny G̃(0,k)=1/k² ⟹ F(s)=1/s tożsamościowo ⟹ biegun na powłoce radiacyjnej s=k²−ω²=0 z residuum 1 ≠ 0. Kontrpróba F=1/(s+a): zabija biegun, ale statyka przechodzi w Yukawę ⟹ niszczy q²=4πG (PR-025 T2a). **Statyka i radiacja są nierozdzielne.** |
| F-AUX-C | **NEGATIVE** | S05 U(1) działa wyłącznie na fazę θ (∂\|Φ\|/∂α = 0 — moduł inwariantny); Z₂ nie ma ciągłego rozszerzenia na rzeczywistym Φ (Im(e^{iε}x)=x·sinε=0 ⟹ ε∈{0,π}) ⟹ brak generatora Noether ⟹ brak więzu pierwszej klasy na mod oddechowy. |
| **F-AUX-D** | **HONEST_NEGATIVE** | „Φ jako zmienna pomocnicza" NIE jest wyprowadzalne w LIVE TGP. **PR-025 upgrade: „both branches" → EXHAUSTIVE-OVER-LIVE.** |

## §1 — Sens fizyczny (jednym akapitem)

OTW utrzymuje potencjał newtonowski jako nieradiacyjny, bo dyfeomorfizmy
czynią komponenty 00/0i więzami — promieniuje tylko TT. LIVE TGP nie posiada
żadnej struktury pełniącej tę rolę dla δΦ: kinetyka jest hiperboliczna bez
więzów (A), ten sam biegun propagatora niesie Newtona i falę (B), a inwentarz
symetrii nie zawiera niczego, co działałoby na moduł Φ (C). Mod oddechowy Φ
jest fizycznym kanałem radiacyjnym z konieczności strukturalnej — kanał
(1/6)P_GR z PR-025 nie ma w LIVE żadnej drogi anihilacji.

## §2 — Zgodność z pre-rejestracją

Pre-derywacja scoping §2 / Phase 0 §4 zrealizowana we wszystkich trzech
punktach. Zero modyfikacji progów, zero nowych stałych, zero dostępu do danych
obserwacyjnych (cykl czysto strukturalny).

## §3 — Honest notes

- B.1 używa kroku „F analityczna na s>0 ⟹ F(s)=1/s wszędzie" (twierdzenie
  o identyczności); sympy weryfikuje algebrę, analityczność jest założeniem
  klasy teorii lokalnych Lorentz-invariantnych — nielokalne G̃ poza zakresem
  (W-AUX-2).
- INFO.1: wszystkie 6 komponent σ_ab propaguje bez więzów — nie-Fierz-Pauli
  masywny tensor generycznie niesie ghost (mod skalarny złej normy). POZA
  zakresem werdyktu (Phase 0 §5.5), zarejestrowane jako W-AUX-3.
