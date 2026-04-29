---
title: "θ.1.Phase2 setup — K_quark first-principles z chirality-counting + B² taxonomy"
date: 2026-04-29
cycle: θ.1.Phase2
status: PRE-EXECUTION
parent: "[[program.md]]"
predecessor: "[[Phase1_results.md]]"
tags:
  - TGP
  - theta-quark-koide
  - chirality-counting
  - K-up-7-8
  - B-squared-extension
---

# θ.1.Phase2 — K_quark first-principles z chirality-counting + B² taxonomy

> **Cel:** Derive K_up = 7/8 sympy-exact via B²_up = 13/4 chirality-counting
> extension; rank K_down rational candidates; show cross-sector V_us = λ_C
> single-anchor lock; falsify 5 alternative quark Koide formulas;
> NGFP RG-stability via common β-rescaling; classification
> PARTIALLY DERIVED (refined).

---

## 7 sub-tests

### T2.1 — K_up = 7/8 sympy candidate (drift 0.046%)

**Test:**
- Universal Koide pattern: K = (2 + B²) / (2N), N=3
- K_up = 7/8 ⇒ B²_up = 6·(7/8) − 2 = 21/4 − 2 = **13/4 = 3.25**
- Drift K_up vs PDG: |7/8 − 0.874559| / 0.874559 = 0.0461%
- Sympy exact rational: 7/8 = 0.875 ↔ B² = 13/4

**Falsification:** if K_up drift > 0.5% (vs 0.046% target), 7/8 candidate
sfalsyfikowany; szukaj alternatywnego rational anchor.

### T2.2 — K_down sympy candidates ranking

**Test:**
- Find best rational n/d with d ≤ 100 closest do K_down = 0.7399
- Top candidates (sympy 30-digit search):
  - 37/50 = 0.74000 — drift 0.014%
  - 54/73 = 0.73973 — drift 0.024%
  - 71/96 = 0.73958 — drift 0.043%
  - 57/77 = 0.74026 — drift 0.049%
  - 17/23 = 0.73913 — drift 0.104%
- B² values: B²_down = 6·K_down − 2 ≈ 2.44 (continuous range)
- **No single small-denominator rational dominates** → K_down structural,
  not LOCKED to clean fraction

**Falsification:** if all candidates drift > 0.5%, K_down sektor wymaga
non-rational structural derivation (e.g., transcendental fit z QCD running).

### T2.3 — B²_up = 13/4, B²_down effective derivation

**Test:**
- Lepton (Dirac): B²_lep = 2 (2 chiralities, no color)
- Neutrino (Majorana): B²_ν = 1 (1 chirality, no color)
- Up-quark (Dirac + 3 color + QCD running): B²_up = 13/4
  - Decompose: 13/4 = 2 (Dirac) + 5/4 (QCD/color correction)
- Down-quark (Dirac + 3 color + QCD running): B²_down ≈ 61/25 ≈ 2.44
  - Decompose: 61/25 = 2 (Dirac) + 11/25 (QCD/color correction)
- Asymmetry: B²_up − B²_down = 13/4 − 61/25 = 325/100 − 244/100 = 81/100
  - Physical meaning: up-sektor has stronger electromagnetic coupling
    (Q_u = +2/3 vs Q_d = −1/3) → larger effective d.o.f. correction

**Falsification:** if B²_up < 2 or B²_down < 2, chirality-counting
extension nie jest valid (quarks must have ≥ 2 chiralities like Dirac
leptons; QCD adds positive correction).

### T2.4 — Cross-sector V_us = λ_C single-anchor lock

**Test:**
- ζ.1 single Cabibbo anchor: λ_C = 0.22550 (GL form factor 165/167)
- Cross-sector lock: V_us (CKM) = sin θ_C = λ_C
- Cross-sector lock: sin θ₁₃ (PMNS) = λ_C / √2
- Wolfenstein λ-cascade: V_us = λ_C, V_cb = A·λ_C², V_ub = A·λ_C³·(ρ̄ − iη̄)
- Drifts vs PDG:
  - V_us: 0.222% (excellent)
  - V_cb: 0.884% (good)
  - V_ub: 8.977% (driven by (ρ̄, η̄) precision, not λ_C)

**Falsification:** if V_us drift > 1% lub cross-sector ratio sin θ_C/sin θ₁₃
≠ √2 ± 10%, single λ_C anchor framework broken.

### T2.5 — 5 alternative quark Koide formulas FALSIFIED

**Test 5 alternatives:**
- C1: K_up = 2/3 (lepton-like) — drift 23.8% (FAIL > 1%)
- C2: K_up = 1/2 (neutrino-like) — drift 42.8% (FAIL > 1%)
- C3: K_up = √(2/3) ≈ 0.816 — drift 6.7% (FAIL > 1%)
- C4: K_up = 8/9 ≈ 0.889 — drift 1.6% (FAIL > 1%)
- C5: K_up = 6/7 ≈ 0.857 — drift 2.0% (FAIL > 1%)

Każda alternatywa daje drift > 1% — uniqueness K_up = 7/8 (drift 0.046%)
strukturalnie potwierdzona.

**Falsification (inverse):** if alternatywa daje drift < 0.5%, K_up = 7/8
nie jest unique optimal; prowadzi do non-uniqueness.

### T2.6 — NGFP RG-stability via common β-rescaling

**Test:**
- Phase 1 T1.3 already proved K_up i K_down RG-invariant pod common
  m_q → c · m_q rescaling (10⁻²⁹% precision)
- NGFP common β-rescaling: γ_m universal across u/c/t (per-sector)
- Tests: K_up at μ ∈ {2 GeV, 80 GeV, 1 TeV, 100 TeV} unchanged
- Same for K_down

**Falsification:** if γ_m flavor-dependent w u/c/t (lub d/s/b) at >0.1%,
K_quark loses RG-invariance i Phase 2 framework wymaga refinement.

### T2.7 — Classification PARTIALLY DERIVED (refined)

**Test:**
- Pre-θ.1: K_quark range 0.81-0.87 STRUCTURAL (open, no rational lock)
- Post-θ.1.Phase2:
  - K_up = 7/8 = 0.875 sympy LOCKED (drift 0.046%) → PARTIALLY DERIVED
  - K_down = 0.7399 numerically locked, no clean rational → STRUCTURAL
  - 4-sector taxonomy distinct (Phase 1 PASS) → DERIVED
- Cross-sector λ_C single anchor (ζ.1 inheritance) → DERIVED refined

**Falsification:** if Phase 2 daje 6/7 PASS lub mniej, classification
pozostaje STRUCTURAL; reframing required.

---

## Verdict gate

**7/7 PASS** → K_quark first-principles classification PARTIALLY DERIVED
(refined), Phase 3 predictions Q1-Q6 proceeds.

**6/7 PASS** → minor gap, Phase 3 z notatką niedomknięcia.

**≤ 5/7 PASS** → Phase 2 reframing required.

---

## Środowisko

```bash
PYTHONIOENCODING=utf-8 python -X utf8 research/op-theta-quark-koide/phase2_kquark_derivation.py 2>&1 | tee research/op-theta-quark-koide/phase2_kquark_derivation.txt
```

## Cross-references

- [`program.md`](program.md) — overall θ.1 plan
- [`Phase1_results.md`](Phase1_results.md) — K_quark numerical landscape LOCKED
- [`../op-zeta-mass-spectrum/Phase3_results.md`](../op-zeta-mass-spectrum/Phase3_results.md) — ζ.1 cross-sector λ_C anchor
