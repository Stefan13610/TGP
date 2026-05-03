---
title: "λ.1 Phase 2 — Sympy LOCK + analytical derivation"
date: 2026-05-01
phase: 2
parent: "[[Phase1_results.md]]"
status: ACTIVE
tags:
  - TGP
  - lambda1
  - phase2
  - analytical-derivation
---

# λ.1 Phase 2 — Sympy LOCK + analytical derivation

## Cel Phase 2

Po Phase 1 (3.0/5 PASS, gate passed), zadanie Phase 2 to **explicit
analytical derivation** mechanizmu który produkuje:

```
n(α) = (e²/4) · (4-α)
Φ_eff = (10/3) · e²
```

**z TGP-fundamentu**, nie tylko jako empirical fit.

---

## Status wejściowy z Phase 1

**Co JUŻ wiemy:**
- L1.5: phase compact U(1) wyklucza e_Euler; amplitude continuous ℝ pozwala
- L1.3: exp pojawia się trywialnie w partition function Z = ∫Dφ·exp(-βH)
- L1.2: γ_φ ∝ (4-α) plausible z Hobart-Derrick scaling
- L1.1: 10/3 = (gravity DOF in 4D)/N_gen counting argument

**Co BRAKUJE:**
- Konkretny mechanizm wykładnik = e²/4 (nie tylko (4-α) factor)
- Explicit derivation Φ_eff = (10/3)·e² z stat-mech
- Wyjaśnienie czemu wartość **2** w `e²` (nie e, nie e³)

---

## Sub-tasks Phase 2 (4)

### P2.1 — Explicit R3 partition function calculation

**Cel:** Pochodzić wave-function renormalization Z_φ z explicit
evaluation log det O dla R3 soliton background.

**Metoda:**
- Linearizacja akcji wokół R3 soliton: S = S_sol + (1/2)·δg·O·δg
- Operator O = -∇² + V''(g_sol(r))
- Heat kernel expansion + ζ-function regularization
- Compute log det O = -ζ_O'(0) jako funkcja g₀
- Sprawdzić czy (1/2)·log det O ~ (e²/4)·(4-α)·log g₀

**PASS criterion:** Show że log det O ma postać proportional do log(g₀)
ze współczynnikiem zalezne od α; bonus jeśli wartość matche e²/4.

**Plik:** `phase2_partition_function_explicit.py`

---

### P2.2 — Multi-loop / semiclassical analysis

**Cel:** Sprawdzić czy e_Euler pojawia się naturalnie w:
- Multi-loop resummed series (e w factorial growth)
- Semiclassical instanton (e^(-S/ℏ) factors)
- Continuous-limit dyskretnego procesu (user's M4)

**Metoda:**
- Standard scalar 3D effective field theory wokół vacuum
- Borel-resummed perturbation series
- Semiclassical action S_inst dla R3-like soliton
- Sprawdzić exp(-S_inst) factor structure

**PASS criterion:** Identify konkretny mechanizm produkujący e_Euler
(nawet jeśli nie e²/4 dokładnie).

**Plik:** `phase2_multiloop_semiclassical.py`

---

### P2.3 — Φ_eff derivation z substrate stat-mech

**Cel:** Explicit derivation Φ_eff = (10/3)·e² z partition function
substrate.

**Metoda:**
- Hamilton H_Γ = -J·Σ(φ·φ)² substrate (z sek10)
- Mean-field decomposition φ_i = √Φ₀ + δφ_i
- Self-consistent equation dla Φ₀ (mean-field gap)
- Bare Φ_0 = 115, screened Φ_eff = 24.66 — pochodzić explicit
- Sprawdzić czy 10/3 = ratio z TGP-stat-mech parameters

**PASS criterion:** Pokazać explicit formula Φ_eff = K · e² gdzie K
wynika z TGP-substrate calculation (nie post-hoc rationalization).

**Plik:** `phase2_phi_eff_derivation.py`

---

### P2.4 — Sympy LOCK + final synthesis

**Cel:** Jeśli P2.1, P2.2 lub P2.3 znajdą mechanizm, sympy LOCK
finalnej formuły. Inaczej — honest negative.

**Metoda:**
- Symbolic verification: empirical formula = derived formula (sympy diff)
- Cross-check: partition function approach + multi-loop approach +
  Φ_eff approach — czy są konsystentne?
- Final synthesis: jeśli pozytywne, λ.1 promotion plan

**PASS criterion:** Sympy verification znajdzie mechanizm który
**reproduces** n(α)=(e²/4)·(4-α) lub Φ_eff=(10/3)·e² (lub oba).

**Plik:** `phase2_sympy_lock_synthesis.py`

---

## Score gate Phase 2

```
≥ 3/4 PASS → Phase 3 forward (predictions + cross-channel)
< 3/4 PASS → λ.1 program END z postmortem
```

**Kill criterion:** Jeśli żaden z P2.1-P2.4 nie produkuje konkretnego
mechanizmu, λ.1 closes z negative result. To jest **honest** outcome —
wszystkie hint'y zostają jako "empirical fit, no fundamental derivation
found w obecnym TGP-formalism".

---

## Order of execution

1. **P2.1** najpierw (najbardziej technicznie konkretny, dobrze
   zdefiniowany problem)
2. **P2.3** drugi (Φ_eff derivation — może być najszybszy)
3. **P2.2** trzeci (najszerszy zakres, multiple sub-mechanisms)
4. **P2.4** ostatni (synteza wszystkich wyników)

---

**Status:** ACTIVE — startujemy 2026-05-01.
