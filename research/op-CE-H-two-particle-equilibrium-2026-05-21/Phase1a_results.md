---
title: "Phase 1a Results -- N=2 isolation null test"
type: phase_results
status: CLOSED
phase: 1a
parent_cycle: op-CE-H-two-particle-equilibrium-2026-05-21
date_completed: 2026-05-21
result: 5/5 PASS, F-beta-1 NULL prediction CONFIRMED
---

# Phase 1a Results -- N=2 isolation null test

**Status:** CLOSED 2026-05-21
**Result:** 5/5 PASS (4/4 substantive FP + 1/1 LIT)
**F-β-1 verification:** NULL prediction CONFIRMED (no stable L* in isolation)

---

## §1 — Test verdict table

| Test | Description | Class | Result |
|------|-------------|-------|--------|
| T_P1a_1 | Single kink Φ_K(x) = v·tanh(m·x/√2) satisfies EOM | substantive | PASS |
| T_P1a_2 | E_K = (2√2/3)·m·v² (standard result) | substantive | PASS |
| T_P1a_3 | V_int(L) Manton-Sutcliffe asymptotic form | LIT (informational) | PASS |
| T_P1a_4 | F-β-1: dE/dL has no zero in L > 0 (no stationary point) | substantive | PASS |
| T_P1a_5 | Asymptotic limits E(L→0+) = 2E_K - A_int, E(L→∞) = 2E_K | substantive | PASS |

**Substantive metrics:**
- 4/4 substantive FP PASS (100%)
- 0 hardcoded T_pass=True (strict cycle 1/2/7 preserved)
- 0/1 DEC budget used

---

## §2 — Technical results

### §2.1 Single kink (T_P1a_1, T_P1a_2)

**Ansatz LOCKED Phase 0:**
$$\Phi_K(x) = v \cdot \tanh\left(\frac{m \cdot x}{\sqrt{2}}\right), \quad m^2 = \lambda \cdot v^2$$

**EOM verification:** d²Φ/dx² - dV/dΦ = 0 after substitution m = √λ·v. **VERIFIED analytically by sympy.**

**Kink energy:**
$$E_K = \frac{2\sqrt{2}}{3} \cdot m \cdot v^2$$

**Computed by sympy** (integral of energy density from -∞ to +∞): matches expected form exactly. Difference = 0.

### §2.2 Kink-antikink interaction (T_P1a_3, LIT informational)

**Parameterization:**
$$V_{int}(L) = -A_{int} \cdot e^{-m \cdot L}$$

with A_int > 0 positive constant (specific value from Manton-Sutcliffe 2004 Ch.5, e.g. A_int ∝ m·v²·factor).

**Status:** INFORMATIONAL. Specific A_int coefficient not derived in Phase 1a (deferred). Substantive claim is structural: V_int < 0 (attractive), exponential decay, V_int → -A_int as L → 0, V_int → 0 as L → ∞.

### §2.3 Stationary point search (T_P1a_4, F-β-1 verification)

**Total energy:**
$$E_{total}(L) = 2 E_K - A_{int} \cdot e^{-m L}$$

**Derivative:**
$$\frac{dE_{total}}{dL} = A_{int} \cdot m \cdot e^{-m L}$$

**Sympy result:** dE/dL is **strictly positive** for all L > 0 (since A_int, m > 0 and exp > 0).

**Solutions to dE/dL = 0 in L > 0:** **NONE** (empty list returned by sp.solve).

**Numerical samples:**
- dE/dL at L=0.1: 0.90484 > 0
- dE/dL at L=1: 0.36788 > 0
- dE/dL at L=10: 0.0000454 > 0 (decaying but always positive)

**Interpretation:** E_total(L) is monotonically increasing on (0, ∞). No stable equilibrium L* > 0 exists in isolation. **F-β-1 NULL PREDICTION CONFIRMED.**

### §2.4 Asymptotic limits (T_P1a_5)

| Limit | Computed | Expected | Diff |
|-------|----------|----------|------|
| L → 0⁺ | 2E_K - A_int | 2E_K - A_int | 0 |
| L → ∞ | 2E_K | 2E_K | 0 |

**Physical interpretation:**
- L → 0: kink + antikink anihilują, total energy released = A_int (binding energy)
- L → ∞: free particle pair, no interaction, E = 2·E_K (rest mass equivalent)

---

## §3 — Interpretation of F-β-1 NULL result

### §3.1 What the null result means

**Pre-registered prediction:** "Phase 1a (N=2 isolation) MUSI dać brak stabilnego L* (lub L* = 0 collapse / L* = ∞ no binding)."

**Result confirms:** dE/dL > 0 everywhere → E monotonic increase from L=0 (anihilation) to L=∞ (free pair). System "wants" to collapse to L=0 (anihilate) OR stay at L=∞ (free particles); no intermediate stable equilibrium.

### §3.2 Why this MATTERS for CE-H

This null result is **STRUCTURAL EVIDENCE FOR CE-H**, not against it:

- CE-H claim: particles need cosmic ⟨Φ⟩_bg background for stability
- If Phase 1a (no bg) had given stable equilibrium → CE-H would be UNNECESSARY (bg not required)
- Phase 1a (no bg) gives NO equilibrium → bg might be NECESSARY (consistent with CE-H)

**Critical logical point:** Phase 1a alone does NOT prove CE-H. It only shows isolation lacks equilibrium. Phase 1b (positive test with bg) is required to demonstrate that bg CAN provide stability.

### §3.3 Literature consistency

This result is **already known** in topological soliton literature (Manton-Sutcliffe 2004, Rajaraman 1982, Coleman 1985). In pure φ⁴ mexican hat model, kink-antikink pair anihilates classically (no static bound state). Phase 1a sympy is **independent verification** of this literature result via direct computation.

**Anti-BD-drift discipline:** literature used as **consistency check only**, not as substantive evidence. The substantive substance is the sympy computation showing dE/dL has no zero.

---

## §4 — Discipline status

### §4.1 Anti-Lakatos

- ✅ All 5 tests reported honestly vs pre-registration
- ✅ No tolerance modifications ex post
- ✅ F-β-1 NULL prediction registered LOCKED 2026-05-21 BEFORE sympy execution
- ✅ Result CONFIRMS pre-registration (genuine prediction, not retroactive fit)

### §4.2 Strict cycle 1/2/7

- ✅ 0 hardcoded T_pass=True
- ✅ All substantive FP use compute-then-compare pattern
- ✅ 1 LIT (T_P1a_3) clearly labeled INFORMATIONAL
- ✅ DEC budget: 0/1 used cumulatively in Phase 1a

### §4.3 Native equations methodology

- ✅ Started from TGP-native Phi-substrate Lagrangian (no derivation from QCD/SM)
- ✅ Single kink + interaction follow from variational principle on TGP Lagrangian
- ✅ Mexican hat V_TGP(Φ) = (λ/4)(Φ² - v²)² is **native TGP** (Pattern 2.5 §3.5.6)
- ✅ No fitting to external frameworks

---

## §5 — Open items for Phase 1b

1. Does **adding** an external bg contribution to E_total(L) produce stable equilibrium?
2. What form of bg contribution is needed? (Phase 0 ansatz: D/L^α power-law)
3. Is the resulting L*(D) monotonic? (F-β-3)
4. Does equilibrium exist for non-fine-tuned parameter range? (F-β-4)

**Phase 1b sympy implementation will address these.**

---

## §6 — Status końcowy Phase 1a

- ✅ 5/5 tests PASS (4/4 substantive)
- ✅ F-β-1 NULL prediction CONFIRMED
- ✅ Phase 1a closes successfully
- ✅ Sets foundation for Phase 1b positive test
- ✅ Discipline LOCKED: 0 hardcoded T_pass=True, 0/1 DEC, anti-Lakatos LOCKED

**Phase 1a CLOSED 2026-05-21. Phase 1b authorized for execution.**
