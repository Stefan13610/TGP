---
title: "NEEDS — op-L01-N1-EM-trace-anomaly-TGP-2026-05-11 (residual sub-needs deferred)"
date: 2026-05-11
parent: "[[./README.md]]"
type: needs
status: 🟢 cycle CLOSED — residual deferred items NIE blokują N1 closure
deferred_items_count: 4
all_blockers_count: 0
tags:
  - needs
  - cycle-closed
  - deferred-precision
  - non-perturbative-extension
---

# NEEDS — residual sub-needs po cycle close

## §0 — Status

🟢 **cycle CLOSED** — STRUCTURAL_DERIVED. Wszystkie six P-requirements RESOLVED.
Residual items poniżej **NIE blokują** N1 closure; są *deferred precision* lub
*future cycle extensions*.

## §1 — Deferred items (NOT blockers)

### N_open.1: Wilson coefficients γ_1, γ_2, γ_3, γ_4 numerical pinning

**Status:** OPEN — deferred precision.

**Problem:** Phase 1+2 establish operator classes structurally:
```
T_anomaly_TGP = (α/(3π))·F² + γ_1·(∂²ψ)·F² + γ_2·(∂_μ∂_ν ψ)·F^{μρ}F^ν_ρ
              + γ_3·σ_ab·F² + γ_4·□F² + Riegert local
```

Numerical values γ_i wymagają explicit Schwinger-DeWitt heat kernel coefficients
w obecności emergent-metric ansatz {A(ψ), B(ψ), C(ψ)} expansion — to jest
multi-loop technical work.

**Co potrzebne:**
1. Explicit Schwinger-DeWitt b_n coefficients dla diagonal Phase 1 ansatz.
2. Verification że γ_i są O(α/(3π)) w magnitude (renormalization-fixed, NIE free).
3. Sub-leading corrections O(α²) at 2-loop level.

**Kandydat dostawcy:** `op-trace-anomaly-precision-extension-TGPxxx` (~3-5 sesji
multi-loop QED w background metryki).

**Typ:** derivation (multi-loop QED, technical).

**Nie blokuje:** N1 closure jest w klasie STRUCTURAL_DERIVED z renormalization-fixed
*structure*. Numerical pinning γ_i to *precision refinement*.

### N_open.2: Magnetar regime B ≳ B_QED non-perturbative analysis

**Status:** OPEN — deferred to non-perturbative cycle.

**Problem:** Phase 3 estimates dla B = 10¹¹ T używają linear extrapolation z 1-loop
formula. B/B_QED ≈ 23 jest deeply non-perturbative regime — Schwinger pair
production back-reaction, Heisenberg-Euler 4-loop+ corrections, vacuum
birefringence.

**Co potrzebne:**
1. Heisenberg-Euler effective Lagrangian z full QED corrections w extreme B.
2. Schwinger pair production rate w magnetar atmosphere conditions.
3. Photon splitting + vacuum birefringence + cyclotron absorption.
4. Cross-check że ρ_EM_quantum/ρ_NS pozostaje ≪ 1 w extreme regime.

**Kandydat dostawcy:** `op-EM-trace-anomaly-Schwinger-extension-TGPxxx` (~4-6
sesji non-perturbative QED w magnetar regime).

**Typ:** derivation (non-perturbative QED + magnetar phenomenology).

**Nie blokuje:** N1 closure status preserved; cycle scope explicit B ≪ B_QED
documented (R5).

### N_open.3: Schwinger-class lab predictions (2030+)

**Status:** OPEN — pending experimental capability.

**Problem:** future lab fields E ~ 10¹⁵-10¹⁸ V/m + B ~ 10⁹ T macroscopic (XCELS,
ELI, MTW@PW, future facilities) byłyby concrete falsifier dla N1 trace anomaly.

**Co potrzebne:**
1. Detailed lab-scale prediction `δΦ` od `ρ_EM_quantum` w E∥B configuration.
2. Atomic clock differential measurement protocol (Schwinger Δclock TT z τ.3 +
   ω.1 cycles).
3. Sensitivity threshold — z jaką precision lab-Schwinger byłby N1 falsifier.

**Kandydat dostawcy:** synthesis cycle z τ.3 + ω.1 + N1 — `op-Schwinger-lab-roadmap-TGPxxx`
(~2-4 sesji synthesis + lab feasibility).

**Typ:** synthesis (forecast + experimental design).

**Nie blokuje:** N1 closure preserved; future test enabling.

### N_open.4: Cross-cycle update tasks (cosmetic)

**Status:** OPEN — propagation tasks.

**Problem:** Phase_FINAL_close.md identyfikuje propagation list:
1. L01 ADDENDUM §3.2 Q3 typo correction
2. L01 NEEDS.md §T.1 (N1) status update CLOSED
3. L01 README.md Q1 status update
4. τ.3 ADDENDUM §2 numerical update
5. PREDICTIONS_REGISTRY.md new entries

**Co potrzebne:** dedicated propagation session — short edits.

**Estymata:** ~30 min — 1 h.

**Nie blokuje:** N1 closure status preserved; updates są cosmetic / cleanup.

## §2 — Pytania otwarte (Q-style)

(Brak nowych Q po Phase 4 closure — wszystkie identyfikowane risks R1-R6 closed
lub honestly documented.)

## §3 — Cross-references

- [[./README.md]] §"Risk flags" R1-R6
- [[./Phase_FINAL_close.md]] §"Cross-cycle propagation list"
- [[./FINDINGS.md]] §6 risk register final status
- [[../op-L01-rho-stress-energy-bridge-2026-05-04/NEEDS.md]] §N1 (parent — *zamknięty*
  przez ten cykl)

---

**NEEDS list:** 4 deferred items, 0 blockers. Cycle CLOSED.
