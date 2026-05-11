---
title: "op-kappa-sigma-2body-PN-2026-05-09 — derivacja κ_σ(η) geometric factor"
date: 2026-05-09
type: research-cycle
priority: P1_CRITICAL
parent: "[[../op-emergent-metric-from-interaction-2026-05-09/Phase_FINAL_close.md]]"
target: "Numerical κ_σ(η=1/4) z 2-body PN binding-energy modification"
classification: NUMERICAL_DERIVATION_CYCLE
status: 🟢 CLOSED — STRUCTURAL DERIVED (heuristic numerical)
folder_status: closed-resolved-heuristic
sympy_total: "7/7 PASS (Phase 1)"
close_date: 2026-05-09
phase_final_close: "[[./Phase_FINAL_close.md]]"
twin_cycle: "[[../op-c0-derivation-from-substrate-2026-05-09/]] (closed heuristic 5/5 PASS)"
joint_result: "c_0·κ_σ = 4/3 EXACT (clean π cancellation, Phase 4 zero-β target reproduced)"
predecessor:
  - "[[../op-emergent-metric-from-interaction-2026-05-09/]] (CLOSED, 57/57 PASS)"
  - "[[../op-c0-derivation-from-substrate-2026-05-09/]] (Phase 1 STRUCTURAL DERIVED, twin cycle)"
related:
  - "[[../op-ppE-mapping/Phase1.5_G_SPA_lock.md]] (G_SPA = 48 sympy-exact)"
  - "[[../closure_2026-04-26/sigma_ab_pathB/]] (σ_ab Path B audit)"
tags:
  - kappa-sigma-derivation
  - 2body-PN
  - sigma-cross-coupling
  - binding-energy-modification
  - numerical-pinning
---

# op-kappa-sigma-2body-PN-2026-05-09

## §0 — Mission

Derive **κ_σ(η=1/4)** geometric factor z 2-body PN binding-energy
modification, gdzie σ-cross-coupling wnosi *strukturalnie nowy* wkład
do orbital E_orb przy 2PN-orbital level.

**Wynik:** Phase 3 emergent-metric LOCK'uje strukturalną formę:
```
Δe_2^σ = c_0 · κ_σ(η)         (linear in c_0, 2PN-orbital)
β_ppE^σ = (45/16)·c_0·κ_σ(η)  (SPA chain → 2.5PN-phase coefficient)
```

Cykl wyprowadza **κ_σ(η)** numerically dla η=1/4 (equal-mass binary).

## §1 — Critical context

### §1.1 — σ_ij^cross structure

Z [[../op-emergent-metric-from-interaction-2026-05-09/Phase3_sympy.py]]:
```
σ_ij^cross = (∂_iΦ_1)(∂_jΦ_2) + (∂_iΦ_2)(∂_jΦ_1) - (2/3)δ_ij(∇Φ_1·∇Φ_2)
```

Anisotropy along binary axis. Magnitude ~ G²M_1M_2/(r_1²r_2²).

### §1.2 — Coupling do orbital E_orb

W effective Lagrangian dla 2-body system:
```
L_orb = (η v²/2) + (M/r) + corrections
```

σ-coupling adds correction at 2PN-orbital (v⁴ ~ U²):
```
L_orb^σ = (η v²/2) + (M/r) + ... + ΔL_σ(c_0, σ_ij^cross, v_i, v_j)
```

Gdzie ΔL_σ ∝ c_0 · σ_ij^cross · v^i v^j (z g_eff_ij·v^iv^j coupling).

### §1.3 — Δe_2^σ derivation (target)

Δe_2 to 2PN-orbital binding-energy correction. From 2-body Lagrangian:
```
ΔE_orb^σ ~ c_0 · ⟨σ_ij^cross · v_i v_j⟩_orbit / (Φ_0² c²)
        ~ c_0 · κ_σ(η) · (η·v²) · v⁴

⟹ Δe_2^σ = c_0 · κ_σ(η)        (Phase 3 emergent-metric LOCK)
```

**κ_σ(η) jest geometric integral:** orbital averaging σ_cross over circular orbit.

## §2 — Strategy

### Step 1: Compute σ_ij^cross at orbital position

Test particle on circular orbit z partner mass at separation r_12.
Compute σ_ij^cross at test particle position:
```
σ_ij^cross|_at_test = M_partner · geometric_function(r_12, position)
```

### Step 2: Contract z orbital velocity

At circular orbit v² = M/r_12 (Newton). Contract σ_ij with v^i v^j:
```
σ_ij^cross v^i v^j ~ M_partner / r_12² · v² · angular_factor
                 ~ M·v⁴/r_12² (at equal mass)
```

### Step 3: Average over orbital phase

Circular orbit: average over θ ∈ [0, 2π]. Anisotropic σ_cross gives
non-trivial angular average.

### Step 4: Identify κ_σ

Match z Δe_2 = c_0 · κ_σ. Read off κ_σ as numerical coefficient.

## §3 — Hypothesis: κ_σ ≈ 1/(3π)?

Z cycle #1 Phase 1 cross-check, **predicted κ_σ ≈ 1/(3π) ≈ 0.106**
(via c_0·κ_σ = 4/3, c_0 ≈ 4π).

**1/(3π)** jest "naturalna" wartość: typowy multipole integral factor 1/π
+ 1/3 z anisotropy averaging. Plausible structural form.

Cykl ma to sprawdzić DERIVATION-WISE, nie a priori postulować.

## §4 — Phase plan

### Phase 0: Setup (this README + scope)

### Phase 1: σ_ij^cross orbital averaging
- Compute σ_cross at orbital probe position (sympy)
- Contract z circular orbital velocity
- Average over orbital phase
- Numerical κ_σ(η=1/4)

### Phase 2: Cross-check c_0·κ_σ z Phase 4 target
- Compute c_0·κ_σ z derived κ_σ + cycle #1 c_0 estimate
- Compare to 4/3 target
- Honest reporting deviation

### Phase 3: Sensitivity analysis + falsifier

### Phase FINAL: classification

## §5 — Probability assessment

| Outcome | Pre-cycle | Notes |
|---|---|---|
| Pełen DERIVED z κ_σ value | 35-50% | If orbital averaging tractable single-session |
| STRUCTURAL CONDITIONAL | 30-40% | If multi-session 2-body Lagrangian derivation needed |
| STRUCTURAL_NO_GO | 10-20% | If σ-coupling structure incompatible z Phase 4 |
| EARLY_HALT | 5-10% | If technical heaviness exceeds session |

## §6 — Time budget

**Estimated:** 3-5 sesji multi-session.
**Single-session realistic scope:** Phase 0-1 (orbital averaging analytical),
numerical κ_σ estimate.

## §7 — Cross-references

- [[../op-c0-derivation-from-substrate-2026-05-09/]] — joint closure dependent
- [[../op-emergent-metric-from-interaction-2026-05-09/Phase3_sympy.py]] — σ_cross derivation
- [[../op-ppE-mapping/Phase1.5_G_SPA_lock.md]] — G_SPA = 48 SPA chain reference
- [[../closure_2026-04-26/sigma_ab_pathB/results.md]] — σ_ab vanishes in static spherical
