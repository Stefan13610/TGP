---
title: "Phase 0 BALANCE (gate) — op-c0-derivation-from-substrate. Cykl DOTYKA statusów epistemicznych (kandydat do downgrade cyklu 2026-05-09 + §3.6.8) ⟹ gate wymagany; lean jawny; budżet stałych 0."
date: 2026-06-22
type: phase_balance
status: 🔒 GATE-OPEN (pre-Phase1)
phase: 0
cycle: op-c0-derivation-from-substrate-2026-06-22
anti_lakatos_lock: ACTIVE
---

# Phase 0 — BALANCE GATE

> Cel gate'u: zanim cykl wytworzy jakąkolwiek liczbę zmieniającą status epistemiczny
> (downgrade cyklu 2026-05-09 lub korekta FOUNDATIONS §3.6.8), jawnie zbilansować, co stoi
> na szali, jaki jest lean i co jest chronione.

## §1 — Czy cykl dotyka statusów epistemicznych?
**TAK.** Werdykt może:
- (a) **downgradować** `op-c0-derivation-from-substrate-2026-05-09` (STRUCTURAL DERIVED → SUPERSEDED),
- (b) **korygować** FOUNDATIONS §3.6.8 („c₀ framework-derivable, deferred") oraz §3.6.10.1–3,
- (c) **rozszerzać** bilans parametrów (jeśli c₀ wolne: budżet 3 → potencjalnie 4, LUB c₀ i κ_E to
  ten sam parametr → budżet pozostaje 3 ale okno recovery = strojenie).

⟹ **Gate WYMAGANY** (reguła domowa: balance przed dotknięciem statusów). Żadna edycja poza
folderem tego cyklu bez osobnej autoryzacji usera (Phase0_lock §6.5).

## §2 — Co na szali (dwustronnie)
| Wynik | Zysk dla TGP | Koszt dla TGP |
|---|---|---|
| c₀ DERIVED (scheme-indep.) | sektor grawitacyjny odzyskuje ≥1 falsyfikowalną predykcję; okno 4/3 testowalne | wymaga jawnej tożsamości/struktury — wysoki próg dowodu |
| c₀ = wolny parametr UV | uczciwość param-counting; spójność z #33/#34/#35 | sektor GW NIE-PREDYKTYWNY; „smoking gun" breathing-mode i c_GW=c pozostają (R2) ale amplituda/faza = fit |
| UNDERDETERMINED | status zachowany; jasny gap | brak postępu netto (ale lokalizacja gap = wartość) |

## §3 — Lean jawny (powtórzony z LOCK §4)
Pesymistyczny: **c₀ = ten sam wolny parametr UV co C_σ** (wiersz 1). Powód do leanu:
c₀ i C_σ są współczynnikami tego samego operatora σ_ab~∂ŝ∂ŝ; F1 dowodzi UV-czułości jego
normalizacji p²; brak znanej tożsamości anulującej. Pozytyw wymaga obalenia tego leanu.

## §4 — Chronione (NIE wchodzą do bilansu downgrade)
- R1 m_σ²=2m_s² (masa, OPE), R2 modi/symetria, R5 struktura σ_ab, R3 matching ξ_eff (warunek), R4 α/Φ.
- GW2/GW3/GW4/GW5/GW6 w registry — niezależne od c₀.

## §5 — Budżet stałych
**0 nowych stałych.** Cykl NIE wprowadza parametrów; ustala status istniejącego (c₀).

## §6 — Warunek przejścia do Phase 1
Spełniony: LOCK §1–§7 zapisany; gate otwarty; lean jawny; zakres i forbidden moves ustalone.
**Phase 1 = analityczna rekoncyliacja** (czy c₀ dziedziczy `D` z F1) + walidacja silnika (reprodukcja −16/35).
