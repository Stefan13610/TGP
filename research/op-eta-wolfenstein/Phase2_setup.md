---
title: "η.1.Phase2 setup — Wolfenstein sympy triple lock + cross-sector cascade (7 sub-tests)"
date: 2026-04-29
cycle: η.1.Phase2
status: PRE-EXECUTION
parent: "[[program.md]]"
predecessor: "[[Phase1_results.md]]"
tags:
  - TGP
  - eta-wolfenstein
  - CKM
  - sympy
  - chirality-counting
---

# η.1.Phase2 — Wolfenstein sympy triple lock + cross-sector cascade (7 sub-tests)

> **Cel:** Lock (A_TGP, ρ̄_TGP, η̄_TGP) jako rational triple z Phase-1
> top-rankings; verify cross-sector denom-family; compute V_ub_TGP
> refined drift; falsify 5 alternatives; classification PARTIALLY DERIVED.

---

## 7 sub-tests / podtesty

### T2.1 — A = 64/81 sympy candidate (first non-trivial)

**Cel:** Test A = 64/81 = 0.79012 jako structural anchor (drift 0.016% z PDG 0.790).
Note: 79/100 trivially exact ale empirical, nie structural; 64/81 = 8²/3⁴
może mieć geometric interpretation.

**Verdict criterion:** drift < 0.05% AND 64/81 admits structural decomposition
(8²/3⁴ = 2⁶/3⁴ — power-of-2 over power-of-3 hint cross-sector?).

### T2.2 — ρ̄ = 11/78 sympy candidate

**Cel:** Test ρ̄ = 11/78 jako structural anchor (drift 0.018% z PDG 0.141).
Note: 78 = 2·3·13. Coprime z 14 (η̄ denom).

**Verdict criterion:** drift < 0.05% AND ρ̄ rationality cross-checked.

### T2.3 — η̄ = 5/14 sympy lock (structural anchor)

**Cel:** Confirm η̄ = 5/14 LOCKED (Phase 1 already showed drift 0.04%, rank #1).
Cross-sector denom-14 analysis: 14 = 2·7. Hint: B²_up = 13/4 = (numerator 13,
B²_lepton = 2 numerator-shared denom 4 ↔ K_up denom 8 = 2³). Czy denom 14
appears anywhere indicating structural cross-link?

**Verdict criterion:** drift < 0.05% (already 0.040%); 5/14 confirmed sympy
exact.

### T2.4 — Cross-sector denom-family analysis

**Cel:** Examine GCD/LCM struktur denom (A, ρ̄, η̄) = (81, 78, 14):
- LCM(81, 78, 14) = ?
- GCD pairwise = ?
- czy istnieje shared base unit?
- Compare z K-taxonomy denoms (8, 50, 1, 2/3 base 6=2N=2·3, 4 lepton)

Alternative: A = 79/100 (trivial PDG fit), denom 100 = 4·25.

**Verdict criterion:** Document structure honestly — pass if cross-sector
analysis yields coherent observation (positive lub negative) z documented
denominators.

### T2.5 — 5 alternative Wolfenstein triples FALSIFIED

**Cel:** Show that alternative triples (A, ρ̄, η̄) drift > 1% in some component:
- C1: (4/5, 1/7, 5/14) = (0.800, 0.1429, 0.3571) — 4/5 drifts A 1.27%
- C2: (3/4, 1/7, 1/π) = (0.750, 0.1429, 0.3183) — drifts in all 3
- C3: (√(5/8), 7/50, 5/14) = (0.7906, 0.14, 0.3571) — root form, irrational A
- C4: (Wolfenstein 2002 stale: 0.808, 0.196, 0.347) — historical
- C5: PDG-2012: A=0.811, ρ̄=0.131, η̄=0.345 — older PDG

Compute drift max(A_drift, ρ̄_drift, η̄_drift) for each — pass if all > 1%.

**Verdict criterion:** All 5 alternatives FALSIFIED z drift > 1% in ≥ 1
component AND TGP-best (64/81, 11/78, 5/14) max-drift < 0.1% across all 3.

### T2.6 — V_ub_TGP refined cascade

**Cel:** Compute V_ub_TGP = A_TGP · λ_C³ · √(ρ̄_TGP² + η̄_TGP²)
z (A_TGP, ρ̄_TGP, η̄_TGP) = (64/81, 11/78, 5/14) i λ_C = 0.22550.
Compare z PDG |V_ub| = 0.00382.

**Verdict criterion:** Drift V_ub_TGP vs PDG < 9% (improvement on θ.1 8.98%).

### T2.7 — Classification PARTIALLY DERIVED (refined)

**Cel:** Verify that:
- (A_TGP, ρ̄_TGP, η̄_TGP) all locked rationals z drift < 0.1%
- 5 alternatives FALSIFIED
- Cross-sector denom-family analysis documented
- V_ub_TGP refined drift improved vs θ.1 baseline
→ classification: **STRUCTURAL → PARTIALLY DERIVED (refined)**.

**Verdict criterion:** All previous 6 tests PASS AND meta-criteria above.

---

## Verdict gate

**7/7 PASS** → η.1.Phase2 CLOSED, K-Wolfenstein triple LOCKED, Phase 3 proceeds.

**6/7 PASS** → η.1.Phase2 closes z minor gap, Phase 3 proceeds.

**≤ 5/7 PASS** → η.1.Phase2 reframing required.

---

## Środowisko

```bash
PYTHONIOENCODING=utf-8 python -X utf8 research/op-eta-wolfenstein/phase2_wolfenstein_derivation.py 2>&1 | tee research/op-eta-wolfenstein/phase2_wolfenstein_derivation.txt
```

## Cross-references

- [`program.md`](program.md) — overall η.1 plan
- [`Phase1_results.md`](Phase1_results.md) — top-5 rationals ranked
- [`../op-theta-quark-koide/Phase2_results.md`](../op-theta-quark-koide/Phase2_results.md) — K_up=7/8 sympy LOCKED z denom 8 (cross-sector hint)
