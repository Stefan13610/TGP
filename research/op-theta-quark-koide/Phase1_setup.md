---
title: "θ.1.Phase1 setup — K_quark numerical landscape audit (5 sub-tests)"
date: 2026-04-29
cycle: θ.1.Phase1
status: PRE-EXECUTION
parent: "[[program.md]]"
tags:
  - TGP
  - theta-quark-koide
  - quark-masses
  - K-taxonomy
  - Wolfenstein
---

# θ.1.Phase1 — K_quark numerical landscape audit

> **Cel:** Verify TGP predictions K_up = 0.8746 i K_down = 0.7398 z PDG
> MS-bar quark masses @ M_Z; audit RG-invariance over 6 orders of magnitude;
> 4-sector K-taxonomy distinct values; Wolfenstein λ-cascade vs PDG.

---

## 5 sub-tests

### T1.1 — K_up = 0.8746 sympy z PDG MS-bar @ M_Z

**Test:**
- Koide K_up = (m_u + m_c + m_t) / (√m_u + √m_c + √m_t)²
- PDG 2024 MS-bar @ M_Z (MeV):
  - m_u = 1.46
  - m_c = 810
  - m_t = 170478
- Numerical: K_up ≈ 0.8746 (literature value, RG-invariant)
- Sympy high-precision evaluation z 30-digit Float

**Falsification:** if K_up drift > 1% z PDG anchors, MS-bar quark mass
inputs nieprawidłowe lub Koide formula nie ma sensu w sektorze quark.

### T1.2 — K_down = 0.7398 sympy z PDG MS-bar @ M_Z

**Test:**
- Koide K_down = (m_d + m_s + m_b) / (√m_d + √m_s + √m_b)²
- PDG 2024 MS-bar @ M_Z (MeV):
  - m_d = 3.17
  - m_s = 63.1
  - m_b = 3089
- Numerical: K_down ≈ 0.7398 (literature value, RG-invariant)
- Sympy high-precision evaluation z 30-digit Float

**Falsification:** if K_down drift > 1% z PDG anchors, MS-bar quark mass
inputs nieprawidłowe lub Koide formula nie ma sensu w sektorze quark.

### T1.3 — K_quark RG-invariance over 6 orders of magnitude

**Test:**
- Common β-rescaling theorem: K invariant pod (m_i → c · m_i) all i
- QCD running m_q(μ) ~ m_q(μ₀) · [α_s(μ)/α_s(μ₀)]^(γ_m/β_0)
- Universal anomalous dimension γ_m equal across u/c/t (and d/s/b) → 
  common multiplicative factor → K_quark RG-invariant
- Verify: K_up at μ=2 GeV vs μ=M_Z vs μ=10 TeV → identical to 10⁻⁵
- Same dla K_down

**Falsification:** if K_quark drift > 0.1% over μ ∈ [2 GeV, 10 TeV] →
common β-rescaling broken (would require flavor-dependent γ_m).

### T1.4 — 4-sector K-taxonomy distinct values

**Test:**
- K_lepton = 2/3 ≈ 0.6667 (Dirac, B²=2)
- K_neutrino = 1/2 = 0.5000 (Majorana, B²=1)
- K_up ≈ 0.8746 (PDG MS-bar @ M_Z)
- K_down ≈ 0.7398 (PDG MS-bar @ M_Z)
- All 4 distinct (no two values within 1%)
- K_up > K_down > K_lepton > K_neutrino (hierarchy)

**Falsification:** if any 2 K values coincide w 1%, sektor-separation
broken (np. K_quark = K_lepton would unify Dirac taxonomies).

### T1.5 — Wolfenstein λ-cascade vs PDG

**Test:**
- Wolfenstein parameterization (PDG 2024):
  - λ = 0.22650 ± 0.00048
  - A = 0.790 ± 0.017
  - ρ̄ = 0.141 ± 0.022
  - η̄ = 0.357 ± 0.011
- TGP single Cabibbo anchor: λ_C = 0.22550 (GL form factor 165/167)
- Cross-sector ζ.1 inheritance: λ_C = sin θ_C = V_us
- λ_C vs PDG λ: drift |0.22550 - 0.22650| / 0.22650 ≈ 0.44%
- V_us = λ_C (TGP)
- V_cb = A · λ_C² (Wolfenstein)
- V_ub ≈ A · λ_C³ (magnitude)

**Falsification:** if cross-sector λ_C-Wolfenstein drift > 1% → ζ.1
single-anchor lock broken; Cabibbo angle nie governs CKM cascade.

---

## Verdict gate

**5/5 PASS** → K_quark numerical landscape LOCKED, Phase 2 proceeds.

**4/5 PASS** → audit gap, Phase 2 deferred lub limited scope.

**≤ 3/5 PASS** → θ.1 reframing required.

---

## Środowisko

```bash
PYTHONIOENCODING=utf-8 python -X utf8 research/op-theta-quark-koide/phase1_kquark_audit.py 2>&1 | tee research/op-theta-quark-koide/phase1_kquark_audit.txt
```

## Cross-references

- [`program.md`](program.md) — overall θ.1 plan
- [`../op-zeta-mass-spectrum/Phase3_results.md`](../op-zeta-mass-spectrum/Phase3_results.md) — ζ.1 program END (predecessor)
- [`../../PREDICTIONS_REGISTRY.md`](../../PREDICTIONS_REGISTRY.md) — C2 Cabibbo entry
