---
title: "Phase 0 BALANCE (gate) — op-sigma-status-propagation-audit"
date: 2026-06-22
type: phase_balance_gate
status: 🔒 GATE PASSED
phase: 0
cycle: op-sigma-status-propagation-audit-2026-06-20
---

# Phase 0 — BALANCE GATE (dotyka statusów rdzenia)

Obowiązkowy gate, bo audyt dotyka **statusów** rdzenia i dokumentów głównych (nie liczb).

## §1 — Co audyt MOŻE zmienić
- Status / framing sektora tensorowego GW (κ_E / C_σ / amplituda) → spójność z F (#33/#34).
- Globalną liczbę parametrów w narracji (2 → 3, lub zakresowanie do sektora skalarnego).
- Annotacje statusu LIVE sektora radiacyjnego (FOUNDATIONS CL-*).

## §2 — Czego audyt NIE MOŻE zmienić
- Żadnej liczby fizycznej (κ_E=5/6, −16/35, m_σ²=2m_s², 12/N_e², r, n_s, …).
- GW4 (LOCKED), GW1/GW2/GW3/GW5/GW6 (niezależne od κ_E).
- Treści dowodów #33/#34 (tylko propagacja ich wniosku).
- Parametryzacji skalarnej (Φ₀, a_Γ vs γ, q — pre-existing, poza zakresem).

## §3 — Symetria edycji (anti-confirmation-bias)
Każda poprawka musi być **neutralna value-blind**: podnosi spójność (κ_E=wolny parametr),
NIE poprawia ani NIE pogarsza „atrakcyjności" teorii selektywnie. Dodanie C_σ jako 3. parametru
**obniża** parsimony ratio (6 → ~4 jeśli liczone globalnie) — audyt zgłasza to wprost,
nie ukrywa (forbidden: cherry-picking korzystnego liczenia).

## §4 — Budżet stałych
Nowych stałych: **0**. Nowych etykiet LaTeX: 0 (wykorzystujemy istniejące rem:sigma-Csigma-free,
rem:sigma-params jako cele cross-ref).

## §5 — Build gate
Po edycjach core: `pdflatex main.tex` musi dać exit 0; pre-existing dangling refs (#32:
ax:substrat, ssec:disformal, app:A-aksjomaty) NIE są regresją tego cyklu.

## §6 — Higiena
Usunąć __pycache__ / zagnieżdżone błędne katalogi jeśli powstaną (brak skryptów w tym cyklu —
audyt czysto tekstowy).

**GATE: PASSED** — audyt może przejść do Phase 1 i (po autoryzacji) do edycji.
