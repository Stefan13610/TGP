---
title: "α.1.Phase1 setup — α_QED numerical landscape audit"
date: 2026-04-30
cycle: α.1.Phase1
status: PRE-EXECUTION
parent: "[[program.md]]"
predecessor: "[[../op-eta-wolfenstein/Phase3_results.md]]"
tags:
  - TGP
  - alpha-fine-structure
  - prime-137
  - audit
---

# α.1.Phase1 — α_QED numerical landscape audit

> **Cel:** Audit liczbowy α_QED⁻¹ ≈ 137.035999084 + cross-sector
> prime-137/prime-7 mapping. Test, czy zeroth-order anchor α_QED⁻¹ ≈ 137
> (z ψ_ph = 160/137 prime denominator) jest strukturalnie spójny.

## 5 sub-tests

### A1.1 — CODATA 2022 + PDG 2024 α_QED reference values

**Test:**
- α_QED⁻¹(0)  = 137.035999084 ± 0.000000021 (CODATA 2022, atomic recoil
  Cs / Rb cross-check 81 ppt precision)
- α_QED⁻¹(M_Z) = 127.952 ± 0.030 (PDG 2024, vacuum polarization running)
- α_QED(0) = 7.297352569·10⁻³

**Gate:** values internally consistent; running α(0) → α(M_Z)
amounts to 7.6% increase via SM vacuum polarization.

### A1.2 — Top-N rational candidates dla α_QED⁻¹ (denom ≤ 200)

**Test:** find best p/q approximations of 137.035999084:
- Brute-force scan p/q dla q ∈ [1, 200]
- Rank by absolute drift |p/q − target|/target
- Identify top-5 unique rationals

**Expected outcomes:**
- 137/1 (trivial) drift 0.0263%
- Other small-denom rationals will appear z drift 0.001–0.1%
- Identify if any "clean" structural rational (like 5/14 dla η̄) emerges

**Gate:** zeroth-order anchor 137/1 OR a clean small-denom rational
fits within structural-anchor band (drift < 0.05%).

### A1.3 — Structural form: 137 = ψ_ph · 137 − ε_ph · 137 sympy exact

**Test:**
- ψ_ph = 160/137, ε_ph = 23/137 (z ε.1 LOCKED)
- ψ_ph · 137 = 160 (sympy exact)
- ε_ph · 137 = 23 (sympy exact)
- (ψ_ph − ε_ph) · 137 = (160 − 23) = 137 (sympy exact)
- α_QED⁻¹(0) − 137 = 0.035999084 (residual)

**Residual cascade probes:**
- 23/640 = 0.0359375 (drift 0.18% vs 0.036)
- 23·(ε_ph/14.7) = ?
- η̄·κ_TGP·something = ?
- ε_ph² · scale where scale = ?

**Gate:** 137 emerges sympy-exact z F4 chain (already DERIVED z ε.1);
residual 0.036 to be addressed Phase 2.

### A1.4 — Cross-sector prime-137 mapping

**Test:** which TGP rationals contain 137 or its factors (137 prime)?
Inventory:
- ψ_ph = 160/137 (z ε.1) — DERIVED
- ε_ph = 23/137 (z ε.1) — DERIVED
- α_QED⁻¹ ≈ 137.036 (PDG)
- Other TGP rationals: A=64/81 (3⁴), ρ̄=11/78 (2·3·13), η̄=5/14 (2·7),
  K_up=7/8 (2³), K_down=37/50, K_lepton=2/3, K_ν=1/2, N_A=500/57 (5³·8/3·19),
  λ_C=165/167 (5·11·3/167)

**Cross-sector primes summary table:**
- 137 — appears w ψ_ph, ε_ph; 137 itself prime
- 167 — appears w λ_C; 167 prime
- 7 — appears w η̄ (=5/14), K_up (=7/8 numerator)
- 3 — appears w A (=81=3⁴), ρ̄ (=78=2·3·13), K_lepton (=2/3 denom)
- 2 — appears w η̄, ρ̄, K_lepton, K_up, K_ν

**Gate:** 137 unique: only ψ_ph, ε_ph share prime-137. No cross-sector
overlap z innych cycles → 137 is **isolated structural prime** w TGP
landscape, anchoring photon-ring sektor.

### A1.5 — α_QED⁻¹ RG-running consistency (α(0) → α(M_Z))

**Test:**
- α(0) = 1/137.035999084 = 7.29735e-3
- α(M_Z) = 1/127.952 = 7.81570e-3
- Running fraction = (α(M_Z)/α(0) − 1) = 7.10% (vacuum polarization
  z lepton + hadron loops)
- TGP test: czy ratio α(0)/α(M_Z) = 137.036/127.952 = 1.07097 admits
  TGP rational anchor?

**Possible anchors:**
- 1.07097 ≈ 1 + ε_ph · 0.4225 (= 1 + 23/137 · 0.4225 ≈ 1 + 0.0709) — interesting!
  Note 0.4225 = α₀/9.575? Or 0.4225 ≈ target_shift_F4 = 57/500 = 0.114? No.
- target_shift_F4 = 0.114; 0.0709 / 0.114 ≈ 0.622
- α₀ = 4.045; ε_ph² · α₀ = (23/137)² · 4.045 = 0.0282·4.045 = 0.114 (= target_shift!)
- So ratio − 1 ≈ ε_ph² · α₀ / 1.61 — likely coincidence

**Gate:** TGP form factors locked under common β-rescaling RG-flow
(NGFP marginal a₂); test if α_QED form anchored by 137 prime is
RG-invariant under common γ_an scale change z UV.1 framework.

## Verdict gate

**5/5 PASS** → Phase 2 first-principles derivation viable.
**4/5 PASS** → revise threshold or document as partial-pass.
**≤ 3/5 PASS** → STRUCTURAL HINT only, terminate research-track.

## Środowisko

```bash
PYTHONIOENCODING=utf-8 python -X utf8 research/op-alpha-fine-structure/phase1_alpha_audit.py 2>&1 | tee research/op-alpha-fine-structure/phase1_alpha_audit.txt
```

## Cross-references

- [`program.md`](program.md) — overall α.1 plan
- [`../op-eps-photon-ring/Phase3_results.md`](../op-eps-photon-ring/Phase3_results.md) — ε.1.E7 prime-137 hint logged
- [`../op-eta-wolfenstein/Phase3_results.md`](../op-eta-wolfenstein/Phase3_results.md) — η.1.H4 prime-7 hint
- [`../../PREDICTIONS_REGISTRY.md`](../../PREDICTIONS_REGISTRY.md) — E7 entry
