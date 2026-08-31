---
title: "Phase1_correction_note — dwa udokumentowane błędy implementacyjne skryptu Phase1_analytic.py (brak założenia g₀>0 → nieredukowalne (g₀²)^{7/2}; zbędny znak − przy konwencji euler_equations w P1b); zero zmian kryteriów"
date: 2026-08-31
type: correction-note
tgp_owner: research/op-symplectic-Jspectrum-2026-08-31
status: APPLIED-BEFORE-USE
related:
  - "[[Phase0_balance.md]]"
  - "[[Phase_method_decisions.md]]"
  - "[[Phase1_analytic.py]]"
---

# Correction note — Phase 1 (zapisana PRZED użyciem poprawionego wyniku)

Pierwszy bieg `Phase1_analytic.py` (output pierwotny ZACHOWANY:
`Phase1_output_run1_original.txt`) dał FAIL w 5 testach: P1b (2)
i P1c E_p/E_q rzędu 0/1 (3). Wzorzec FAIL-ów (przechodzą wszystkie
tożsamości niezależne od U(ρ) oraz P1c-CENTRALNA i Wirtinger; padają
wyłącznie porównania zawierające U(ρ) po stronie euler_equations)
wskazywał artefakt implementacyjny, nie błąd derywacji. Diagnoza
(sesja 2026-08-31, reprodukowalna):

1. **Brak założenia pozytywności pola.** `g_0 = Function('g_0')(x)`
   bez `positive=True` ⟹ sympy NIE redukuje `sqrt(g_0²)⁷ = (g_0²)^{7/2}`
   do `g_0⁷` (test: `simplify((g_0²)^{7/2} − g_0⁷)` = 0 tylko przy
   `positive=True`). U = βρ⁷/7 − γρ⁸/8 ma nieparzyste potęgi ρ ⟹
   wszystkie porównania z U po stronie ρ=√(p²+q²) fałszywie niezerowe.
   Fizycznie g_d > 0 na wszystkich tłach (g_min ≥ 0.167) — założenie
   poprawne. KOREKTA: `positive=True` dla pól tła w P1b/P1c.

2. **Zbędny znak − w P1b.** Konwencja sympy: `euler_equations(L,[f],[x])`
   zwraca lhs = ∂L/∂f − d/dx ∂L/∂f′, czyli DOKŁADNIE δE/δf dla E=∫L dx
   (test: L=½f′² ⟹ lhs = −f″ = δ∫½f′²/δf ✓). Skrypt w P1b brał
   `dEdg = −lhs` (błędny znak); w P1c znak był już poprawny
   (`E_p = eqs[0].lhs`) — pozostała po redakcji martwa podwójna
   przypiska E_p do usunięcia. KOREKTA: `dEdg = +lhs`, czyszczenie
   martwego kodu.

Żadne kryterium, próg, tolerancja ani lista teł NIE ulega zmianie;
korekta dotyczy wyłącznie założeń symbolicznych i znaku odczytu
konwencji bibliotecznej. Testowane tożsamości pozostają identyczne.
Pierwotny output zachowany obok poprawionego (`Phase1_output.txt`).
