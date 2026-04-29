---
title: "η.1.Phase1 setup — Wolfenstein (A, ρ̄, η̄) numerical landscape audit (5 sub-tests)"
date: 2026-04-29
cycle: η.1.Phase1
status: PRE-EXECUTION
parent: "[[program.md]]"
predecessor: "[[../op-theta-quark-koide/Phase3_results.md]]"
tags:
  - TGP
  - eta-wolfenstein
  - CKM
  - Wolfenstein
  - audit
---

# η.1.Phase1 — Wolfenstein (A, ρ̄, η̄) numerical landscape audit (5 sub-tests)

> **Cel:** Verify PDG 2024 Wolfenstein values; rank top-5 rational candidates
> per parameter; audit cross-sector cascade z λ_C lock; unitarity triangle
> apex closure; Jarlskog J consistency.

---

## 5 sub-tests / podtesty

### T1.1 — A = 0.790 ± 0.012 PDG audit + top-5 rational candidates

**Cel:** Sympy 30-digit precision verify A = 0.790 PDG; rank top-5 rational
candidates (denom ≤ 100) z drift < 2%.

**Method:** Stern-Brocot rational-approximation search; sort by drift_pct.

**Pass criterion:** ≥ 1 candidate z drift < 1%; top-5 ranked.

### T1.2 — ρ̄ = 0.141 ± 0.020 PDG audit + top-5 rational candidates

**Cel:** Sympy 30-digit verify ρ̄ = 0.141 PDG; rank top-5 rational candidates
(denom ≤ 100) z drift < 2%.

**Pass criterion:** ≥ 1 candidate z drift < 2%; top-5 ranked.

### T1.3 — η̄ = 0.357 ± 0.014 PDG audit + top-5 rational candidates

**Cel:** Sympy 30-digit verify η̄ = 0.357 PDG; rank top-5 rational candidates
(denom ≤ 100); test 5/14 hypothesis explicitly.

**Pass criterion:** 5/14 ranks #1 z drift < 0.1%; top-5 listed.

### T1.4 — Unitarity triangle apex (ρ̄ + iη̄) closure α + β + γ = π

**Cel:** Verify CKM unitarity triangle z PDG (A, λ, ρ̄, η̄):
- α = arg(-V_td·V_tb* / V_ud·V_ub*)
- β = arg(-V_cd·V_cb* / V_td·V_tb*)
- γ = arg(-V_ud·V_ub* / V_cd·V_cb*)
- Σ = α + β + γ should equal π exact (CKM unitarity)

**Method:** Compute V_td, V_tb*, V_ud, V_ub*, V_cd, V_cb* z Wolfenstein
parameterization (4-th order); arg via sp.atan2; sum vs π.

**Pass criterion:** |Σ - π| < 10⁻³ rad (numerical closure).

### T1.5 — Jarlskog J = A²·λ_C⁶·η̄ z PDG values vs PDG J

**Cel:** Compute J z PDG (A, λ_C, η̄); compare z PDG-2024 measured
J = (3.07 ± 0.10)·10⁻⁵.

**Method:**
- J_TGP = A² · λ_C⁶ · η̄ (Wolfenstein definition)
- Drift_pct = |J_TGP - J_PDG| / J_PDG · 100

**Pass criterion:** Drift < 5%.

---

## Verdict gate

**5/5 PASS** → Phase 2 proceeds z best-candidate triple (A, ρ̄, η̄).

**4/5 PASS** → Phase 2 proceeds z minor gap (note in classification).

**≤ 3/5 PASS** → η.1.Phase1 reframing required.

---

## Środowisko

```bash
PYTHONIOENCODING=utf-8 python -X utf8 research/op-eta-wolfenstein/phase1_wolfenstein_audit.py 2>&1 | tee research/op-eta-wolfenstein/phase1_wolfenstein_audit.txt
```

## Cross-references

- [`program.md`](program.md) — overall η.1 plan
- [`../op-theta-quark-koide/Phase3_results.md`](../op-theta-quark-koide/Phase3_results.md) — θ.1 predecessor (V_ub gap)
- [`../../PREDICTIONS_REGISTRY.md`](../../PREDICTIONS_REGISTRY.md) — H1-H6 future entries
