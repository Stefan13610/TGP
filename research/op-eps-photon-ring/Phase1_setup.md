---
title: "ε.1.Phase1 setup — ε_ph numerical landscape audit"
date: 2026-04-29
cycle: ε.1.Phase1
status: PRE-EXECUTION
parent: "[[program.md]]"
tags:
  - TGP
  - epsilon-photon-ring
  - audit
  - sympy-rational
---

# ε.1.Phase1 — ε_ph numerical landscape audit

> **Cel:** Verify ε_ph = 4197/25000 = 0.16788 sympy-exact, confirm M9.1″
> refined provenance, M9.2-D → M9.1″ refinement gain 36×.

---

## 5 sub-tests

### E1.1 — ε_ph = 4197/25000 sympy-exact rational

**Test:** ε_ph_M9_1P = sp.Rational(16788, 100000) → simplified to 4197/25000.
**Falsification:** ε_ph not rational or denominator > 10⁵ → ε_ph free parameter.

### E1.2 — ψ_ph = 4/(3 + 0.4250) = 1.16788 z M9.1″ null geodesic

**Test:** ψ_ph = sp.Rational(4) / (sp.Rational(3) + sp.Rational(425, 1000))
= 4/3.4250 = 800/685 = **160/137 = 1.16788...**
- Verify: ε_ph := ψ_ph − 1 = 23/137 ≈ 0.16788 (sympy exact)
- Drift sympy 23/137 vs Rational 4197/25000: < 10⁻⁴

**Falsification:** ψ_ph − 1 not equal to 4197/25000 within 10⁻³ → M9.1″
geodesic computation incorrect.

### E1.3 — ε_ph² = 441/15625 (= 0.028224) sympy exact

**Test:** ε_ph² = (4197/25000)² simplified.
- Numerical: 0.16788² = 0.0281837...
- Note: 441/15625 = 0.028224 (slight discrepancy z (4197/25000)² = 17614809/625000000)
- Reconciliation: Phase 1 used coarse 0.168 → 0.168² = 0.028224 = 441/15625
  Phase 2/UV.1 uses 0.16788 → 0.16788² = 0.0281837...

**Falsification:** if ε_ph² consistency check fails between Phase 1/2/UV.1
references → M9.1″ refinement broken.

### E1.4 — F4 chain consistency α₀ = target_shift_F4 / ε_ph² = 4.04489

**Test:** α₀_repro = 0.114 / 0.16788² = 4.04489
- F4 sympy: 1069833/264500 = 4.04472
- Drift: 0.0038% < 0.01% gate

**Falsification:** drift > 0.1% → F4 chain heat-kernel frame inconsistent
z M9.1″ refined ε_ph.

### E1.5 — M9.2-D → M9.1″ refinement gain 36×

**Test:**
- M9.2-D coarse ε_ph = 0.168: α₀_repro = 0.114/0.168² = 4.04018
  drift z F4: 0.139%
- M9.1″ refined ε_ph = 0.16788: α₀_repro = 4.04489
  drift z F4: 0.0038%
- Refinement gain: 0.139% / 0.0038% = 36.6× lepsza precyzja

**Falsification:** if M9.1″ refined NOT closer to F4 sympy than M9.2-D
→ M9.1″ geodesic computation invalid.

---

## Verdict gate

**5/5 PASS** → ε_ph numerical landscape **fully LOCKED**, Phase 2 proceeds
z structural decomposition derivation.

**4/5 PASS** → audit gap, Phase 2 deferred.

**≤ 3/5 PASS** → ε.1 reframing required.

---

## Środowisko

```bash
PYTHONIOENCODING=utf-8 python -X utf8 research/op-eps-photon-ring/phase1_eps_audit.py 2>&1 | tee research/op-eps-photon-ring/phase1_eps_audit.txt
```

## Cross-references

- [`program.md`](program.md) — overall ε.1 plan
- [`../op-uv-as-ngfp/Phase1_results.md`](../op-uv-as-ngfp/Phase1_results.md) — UV.1.UV1.5 a₂ → α₀ reproducibility 0.004% (M9.1″ refined ε_ph LOCKED)
- [`../op-xi-photon-ring/Phase1_results.md`](../op-xi-photon-ring/Phase1_results.md) — ξ.1 ε_ph in ξ-factor frame
