---
title: "Phase 2 — F-FM-V: creation velocity — TIEBREAKER_DERIVED (CONDITIONAL): v_c = 2c/3 EXACT (B-k4); collapsed prediction log₁₀G = 2.025"
type: phase_result
status: PHASE2_COMPLETE
phase: 2
cycle: op-frontier-microphysics-2026-06-11
created_date: 2026-06-11
authorization: "User 2026-06-11: 'FM Phase 2'"
sympy_script: "[[./Phase2_sympy.py]]"
sympy_output: "[[./Phase2_sympy.txt]]"
sympy: "8/8 PASS; 0 hardcoded T_pass"
falsifier_resolved: "F-FM-V = TIEBREAKER_DERIVED (CONDITIONAL on A-ii exact + A-iv monochromatic entry): v_c = 2c/3 EXACT ⇒ ε = 2/3 ⇒ B-k4; B-k3 EXCLUDED unconditionally (K1 + K4a). PR-022 condition (i) CONDITIONALLY satisfied. NO PR-022 (conditions ii-iv outstanding)"
anti_lakatos_lock: PRESERVED
---

# Phase 2 — z jaką prędkością wchodzi materia kreowana? (tiebreaker B-k3 vs B-k4)

## §0 — Verdict at a glance

**F-FM-V = TIEBREAKER_DERIVED (CONDITIONAL): v_c = 2c/3 EXACT ⇒ ε = 2/3 ⇒ B-k4.**
Zbiór dwupunktowy FCR kolabuje. Predykcja bezparametrowa:

```
log₁₀G = p·log₁₀(1091),  p = 2/3 EXACT (EdS)  ⇒  log₁₀G = 2.025
```

**PASS_BAND (krawędź: 0.025 dex), 0.97 dex PONIŻEJ obserwowanego 3.0** — kryteria value-blind
wybrały punkt DALSZY od obserwacji, dokładnie w kierunku pre-flagowanym w Phase 0 §7.
Raportowane bez korekt.

## §1 — Wykluczenie B-k3 (bezwarunkowe; K1 + K4a)

1. **K1 (masywność przy zdarzeniu księgowym):** ledger event = pojawienie się dM > 0.
   Z Phase 1 FP6: stabilna masywna materia tylko przy |Φ| > Φ₀/√3 (wewnątrz ściany);
   przy wejściu do bulku m_eff → m_Φ = √(2λ)Φ₀ > 0. Element masywny ⇒ |v| < c ŚCIŚLE
   ⇒ **dokładne v_c = c nie ma fizycznego nośnika** (newtonowski zapis (1/2)c² nie odpowiada
   energii kinetycznej żadnego elementu).
2. **K4(a) (audyt radiation-first):** jedyna droga ratunku B-k3 — materia rodzi się bezmasowo
   (wchodzi przy c), kondensuje do solitonów później w bulku. Ale kondensacja = zdarzenie
   dM > 0 we wnętrzu bulku = **kreacja w bulku — zakazana przez A-i LOCKED-claim** („bulk
   saturated, blocked"). Pojawienie się masy jest ograniczone do warstwy ściany (grubość
   O(δ), δ/ct → 0; głębokość krawędzi stabilności x* ≈ 0.659δ) ⇒ przy wejściu element JEST
   masywny ⇒ **escape hatch CLOSED; B-k3 EXCLUDED**.

Wykluczenie czysto semantyczno-strukturalne (LOCKED źródła), zero odniesień do wartości wzrostu.

## §2 — Pozytywna selekcja v_c = 2c/3 (dwie niezależne zbieżne linie)

**Linia 1 — K2 (mean-flow boundary value):** ciągłość bez źródeł + jednorodność ⇒ ∇·u = 2/t;
ansatz izotropowy u_r = f(t)x + g(t)/x² ma człon g dywergencyjnie zerowy, zabity regularnością
w x = 0 ⇒ **u = (2/3)x/t UNIKALNE**. Materia na powierzchni wejścia x = ct⁻ to w całości
materia świeżo kreowana — jej średnia prędkość MUSI być wartością brzegową przepływu:
**v_c = u(ct, t) = 2c/3 EXACT**. Warunkowo: **A-ii exact** (jednorodność — COR-1, Phase 3)
+ **A-iv** (wejście monochromatyczne, zerowa dyspersja — NOWE założenie, zadeklarowane).

> **Uczciwa nota ewaluacyjna K2 (§3.6 anti-goalpost):** pierwotne sformułowanie K2 („mismatch
> wymaga warstwy dragu") samo w sobie NIE wiąże — prędkości pekuliarne zanikają adiabatycznie
> v_pec ∝ 1/a (FP7: wykładniki trajektorii {2/3, 1/3}; moda 1/3 ⇒ v ∝ t^(−2/3) = a⁻¹), więc
> mismatch mógłby zanikać bez dragu. K2 wiąże przez argument wartości brzegowej średniego
> przepływu (powyżej). To samo kryterium, uczciwie wyliczone — mechanizm udokumentowany,
> żadne nowe kryterium nie zostało dodane.

**Linia 2 — K3 (zerowa siła substratowa, model-independent):** z EQ-1 energia solitonu jest
funkcjonałem lokalnego tła: E_sol = E(⟨Φ⟩(x)) ⇒ siła F = −E′(⟨Φ⟩)·∇⟨Φ⟩. W nasyconym bulku
⟨Φ⟩ = Φ₀ = const ⇒ ∇⟨Φ⟩ = 0 ⇒ **F_substrat = 0 EXACT dla DOWOLNEGO E** (zero inputu
modelowego). Wobec tego tło Eulera musi domykać się bez źródeł: residual przepływu
u = (2/3)x/t względem grawitacji g = −(ε/3)x/t²:

```
res = ∂u/∂t + u·u′ − g = (3ε−2)/9 · x/t²  ⇒  res = 0 ⟺ ε = 2/3 (jedyne zero)
```

— odtwarza LOCKED formułę Δ_bulk = |res|/|u·u′| = |3ε−2|/4 (Δ(3/2) = 5/8 ✓ FCR). **Niezależna,
wymuszona selekcja ε = 2/3.** Korolarium: COR-2/C-2 — „wymagana siła substratowa O(1)H²x"
okazuje się być ZEREM, a balans zachodzi tożsamościowo przy ε = 2/3 (formalne zaksięgowanie
w Phase 3).

**NEW_POINT nie zrealizowany** (forbidden move #7 — sprawdzone): średnia prędkość wejścia
= wartość przepływu = dokładnie LOCKED punkt B-k4.

## §3 — Domknięcie księgowe: marginalność ⟺ krytyczność (FP5)

Marginalność (LOCKED FCR P3) przy v_c = 2c/3: M = 2c³t/(9G) (ξ = 2/9) ⇒
ρ̄ = M/((4π/3)(ct)³) = **1/(6πGt²) = 3H_m²/(8πG) EXACT** przy H_m = 2/(3t) — sektor materii
jest **dokładnie krytyczny względem własnego wyprowadzonego przepływu** (ta sama tożsamość,
którą standardowa kosmologia definiuje ρ_crit); pętla ε = 4πGρ̄t² = 2/3 domyka się. Punkt
B-k4 = dokładny sektor EdS zanurzony w przestrzeni generowanej R = ct.

## §4 — Predykcja po kolapsie (mechanicznie; bands LOCKED)

| Punkt | p | log₁₀G | Band | Status |
|---|---|---|---|---|
| **B-k4 SELECTED** | **2/3 EXACT (EdS)** | **2.025** | PASS_BAND (krawędź 0.025 dex) | **predykcja po kolapsie** |
| B-k3 excluded | (√55−1)/6 | 3.249 | PASS_BAND | for the record |

Odległość od obserwowanego log₁₀G = 3.0: **0.97 dex PONIŻEJ** (czynnik ~9.4 w G). Uczciwie:
PASS jest pasmowo-mechaniczny; fizyczna rozbieżność wymaga otwartej dyskusji na etapie PR
(NIE ratujemy — zapisujemy).

## §5 — Noty strukturalne (FP7)

- **Numerologia Schwarzschilda rozpuszczona:** przy M = 2c³t/9G: 2GM/c² = (4/9)ct < R = ct —
  wszechświat NIE jest w warunku horyzontu; identyfikacja „zero-energy universe" z FCR B-k3
  nie jest zrealizowana (żaden LOCKED werdykt tego nie twierdził — nota porządkowa).
- **M3 (wall-energetics): NON-DISCRIMINATING** — podział energii kinetyczna/wewnętrzna wolny
  na tym poziomie; pełny budżet ściany = Phase 3 A2 (z zadeklarowanym napięciem R-FM-3:
  ΔV = const vs ρ̄c² ∝ t⁻²).

## §6 — Status warunków PR-022 (po Phase 2)

| Warunek | Status |
|---|---|
| (i) tiebreaker derived | **CONDITIONALLY SATISFIED** (v_c = 2c/3; warunki: A-ii, A-iv) |
| (ii) A-ii homogeneity | OPEN → Phase 3 (COR-1) |
| (iii) C-2 substrate balance | **PRE-RESOLVED** (F_sub = 0 forced; formalne zaksięgowanie Phase 3) |
| (iv) A2 Φ→matter bridge | OPEN → Phase 3 (COR-3; napięcie R-FM-3 zadeklarowane) |

**NO PR-022** (forbidden move #6). Predykcja-kandydat po kolapsie: **log₁₀G = 2.025
bezparametrowo** — jeśli Phase 3 domknie (ii)-(iv).

## §7 — Anti-Lakatos (Phase 2): COMPLIANT ✓

0/8 hardcoded; B-k3 wykluczony semantycznie (LOCKED A-i + masywność), nie wartością wzrostu;
selekcja idzie PRZECIW obserwacji (2.025 vs 3.0 — kierunek pre-flagowany Phase 0 §7, zapisany
PRZED wyprowadzeniem); K2 uczciwie osłabione i przewiązane z dokumentacją mechanizmu (bez
nowych kryteriów); A-iv zadeklarowane jako nowe założenie (nie ukryte); G_obs wyłącznie
porównawczo (FP8 audit); 0 predecessor verdicts modified; 0 nowych stałych; oba punkty
raportowane z odległościami.
