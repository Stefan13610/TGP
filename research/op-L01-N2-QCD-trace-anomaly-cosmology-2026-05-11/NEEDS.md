---
title: "NEEDS — op-L01-N2-QCD-trace-anomaly-cosmology-2026-05-11 (residual deferred)"
date: 2026-05-11
parent: "[[./README.md]]"
type: needs
status: 🟢 cycle CLOSED — all sub-needs N0.1-N0.12 resolved; residual = deferred precision items NOT blocking
deferred_items_count: 4
all_blockers_count: 0
tags:
  - needs
  - cycle-closed
  - deferred-precision
  - extension-cycles
---

# NEEDS — residual sub-needs po cycle close

## §0 — Status

🟢 **cycle CLOSED 2026-05-11** — STRUCTURAL_DERIVED. Wszystkie 12 sub-needs
N0.1-N0.12 resolved. 6/6 P-requirements RESOLVED. Residual items poniżej **NIE
blokują** N2 closure; są *deferred precision* items (lattice external) lub *extension cycles* dla SM sektorów outside QCD.

## §1 — Deferred items (NOT blockers)

### N_open.1: Lattice QCD precision refinement

**Status:** OPEN — deferred precision (analogous do N1 Wilson γ_i deferred precision).

**Problem:** Phase 1+2+3 use lattice QCD inputs (Λ_QCD, ⟨α_s G²/π⟩, T_c, EoS Δ(T)
profile) z documented uncertainty bands. Multi-loop β_QCD + sub-percent lattice
precision could refine numerical predictions.

**Co potrzebne:**
1. Multi-loop β_QCD (b_1, b_2, b_3 coefficients): Λ_QCD precision do percent.
2. Sub-percent lattice EoS Δ(T) profile (HotQCD post-2024 expected).
3. Higher-order Wilson coefficients dla curvature × G² terms (analogous N1 γ_i).

**Kandydat dostawcy:** `op-QCD-trace-anomaly-precision-extension-TGPxxx`
(~3-5 sesji multi-loop QED + lattice precision integration).

**Nie blokuje:** N2 closure jest STRUCTURAL_DERIVED z renormalization-fixed structure.

### N_open.2: Higgs sector extension (L01 N4)

**Status:** OPEN — deferred do separate cycle.

**Problem:** N1 + N2 cycles closed dla EM (gauge U(1)) i QCD (gauge SU(3)) sektorów.
Higgs sektor (scalar field z SSB) wymaga separate extension cycle.

**Co potrzebne:** 1-loop Higgs sector trace anomaly w `g_eff[{Φ_i}]` z SSB
cancellation framework + quantum fluctuations h(x) wokół vacuum.

**Kandydat dostawcy:** `op-Higgs-trace-anomaly-extension-TGPxxx` (~3-5 sesji).

**Nie blokuje:** L01 NEEDS N4 osobny.

### N_open.3: EW gauge anomaly extension (L01 N5)

**Status:** OPEN — deferred do separate cycle (depends on N2 architecture extension).

**Problem:** SU(2)×U(1) electroweak gauge anomaly + EW phase transition cosmology.

**Co potrzebne:** Extension N2 architecture do non-Abelian SU(2)×U(1); EW phase
transition (T_EW ~ 100 GeV); LISA stochastic GW background prediction (jeśli
first-order EW transition).

**Kandydat dostawcy:** `op-EW-trace-anomaly-extension-TGPxxx` (~3-5 sesji).

### N_open.4: Cross-cycle propagation tasks (cosmetic)

**Status:** OPEN — propagation tasks po cycle close.

**Items dla Phase_FINAL §7:**
1. L01 NEEDS.md §N2 status update IN_PROGRESS → CLOSED
2. L01 README.md post-N2 closure note + cross-cycle convergence diagnostic update (5-fold → 6-fold)
3. L01 NEEDS §T.2 (N2 three-layer) update z konstruktywnymi wynikami
4. Q2 cycle FINDINGS reciprocal verification reference
5. PREDICTIONS_REGISTRY.md new entries M911-QCD-* (4 entries)
6. closure_2026-04-26 T-Λ cross-link

**Estymata:** ~30 min — 1 h.

**Nie blokuje:** N2 closure status preserved; updates są cosmetic.

## §2 — Pytania otwarte (Q-style)

(Brak nowych Q po Phase 4 closure — wszystkie identyfikowane risks R1-R7 closed
lub honestly documented.)

## §3 — Cross-references

- [[./README.md]] §"Risk flags" R1-R7
- [[./Phase_FINAL_close.md]] §7 cross-cycle propagation list
- [[./FINDINGS.md]] §6 risk register final status
- [[../op-L01-rho-stress-energy-bridge-2026-05-04/NEEDS.md]] §N2 (parent — *zamknięty*
  przez ten cykl)

---

**NEEDS list:** 4 deferred items, 0 blockers. **Cycle CLOSED 2026-05-11.**
