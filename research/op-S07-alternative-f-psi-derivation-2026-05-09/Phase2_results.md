---
title: "Phase 2 results — constraint testing summary"
date: 2026-05-09
parent: "[[README.md]]"
type: phase2-results
phase: 2
sympy_status: "13/13 baseline + 11/11 obstruction = 24/24 PASS"
---

# Phase 2 results — constraint testing + structural obstruction

## §0 — Executive summary

Phase 2 ukończone w SCOPED-mode (per Phase2_setup §2 decision):

| Sub-phase | Sympy | Status |
|---|---|---|
| Baseline reproduction (M9.1'') | 13/13 | ✅ PASS |
| Structural obstruction analysis | 11/11 | ✅ PASS |
| **Total** | **24/24** | **✅ PASS** |

**Key findings:**
1. M9.1''-class framework consistency confirmed (sympy verified)
2. R3 ODE is **f-independent** in M9.1''-class (anti-podal + universal U_eff)
3. Newton matching forces c_1 = -2/b_1 (f-dependent through Taylor coef)
4. F1 (GR-clone) identified as priority candidate for Phase 3
5. Full alternative derivation requires matter-coupled EOM (multi-session)

## §1 — Phase 2 baseline (M9.1'' reproduction)

**Script:** [[Phase2_baseline_M911_reproduction_sympy.py]]
**Output:** [[Phase2_baseline_M911_reproduction_sympy.txt]]

**Verified:**
- f_M911 * h_M911 = 1 (anti-podal)
- f(1) = 1, h(1) = 1
- V_M911(1) = -1/12 (gamma=1)
- U_eff_M911 = V_M911 * h_M911 = ψ⁴/4 - ψ³/3 (G.0 universal form)
- α_n^M911 (from Phase 1.5 LOCK 7/7): -2, +2, -7/3 at n=1, 2, 3
- Δα_3 = -5/6 (5σ falsification source)
- β_ppE^M911 = 15/4 = 3.75 (vs GWTC-3 bound 0.78 — violation 4.81x)

**Conclusion:** framework infrastructure consistent, Phase 1.5 LOCK reproduces.

## §2 — Phase 2 structural obstruction

**Script:** [[Phase2_structural_obstruction_sympy.py]]
**Output:** [[Phase2_structural_obstruction_sympy.txt]]

### §2.1 — Key structural insight

**M9.1''-class assumptions:**
- (i) K(ψ) = ψ⁴ (T-D-uniqueness, C1)
- (ii) Anti-podal f·h = 1 (M9.1'' structural feature, C9 default)
- (iii) Static EOM = R3 ODE (G.0 universal projection)
- (iv) Vacuum at ψ=1 (C7 baseline, m_sp² = +γ > 0)

**Result:** under these constraints, **V_grav = U_eff · f UNIQUELY determined by f**:

```
V_grav(ψ) = (ψ⁴/4 - ψ³/3) · f(ψ)
```

This is the structural counterpart to G.0 P21 V_M911 derivation —
generalized to any anti-podal f.

### §2.2 — Vacuum stability (C7) is universal

For ANY f(ψ) in M9.1''-class, linearized EOM around ψ=1 gives:
```
ε'' + (2/r) ε' = -ε + O(ε²)
m_sp² = +1 (in γ=1 units), INDEPENDENT of f
```

⟹ C7 automatically satisfied for entire M9.1''-class.

### §2.3 — Newton matching algebraic constraint

For 1PN matching α_1 = -2:
```
α_1 = b_1 · c_1 = -2  ⟹  c_1 = -2/b_1
```

For 1PN β_PPN = 1 (α_2 = +2):
```
α_2 = b_1·c_2 + b_2·c_1²/2 = +2
⟹ c_2 = 2/b_1 - 2·b_2/b_1³
```

### §2.4 — Δα_3 = 0 strategy: structural test

For Δα_3 = 0 (i.e., α_3 = -3/2 = GR exact):
```
α_3 = b_1·c_3 + b_2·c_1·c_2 + b_3·c_1³/6 = -3/2
```

Solving for c_3:
```
c_3 = (-3/2 - b_2·c_1·c_2 - b_3·c_1³/6) / b_1
```

**KEY**: c_3 is also determined by EOM nonlinear structure. Consistency of
these TWO equations (Newton-matched c_3 vs EOM c_3) is the structural test.

### §2.5 — M9.1'' verification

For M9.1'' (b_1=-4, b_2=8, b_3=-24):
- Newton: c_1 = 1/2 ✓
- 1PN: c_2 = -1/4 ✓
- α_3 = -7/3 (Phase 1.5 LOCK) gives c_3 = 5/24
- For α_3 = -3/2 (GR match), would require c_3 = 0
- **Discrepancy: c_3_actual = 5/24 ≠ 0 = c_3_required** → Δα_3 = -5/6 ≠ 0 ✓

### §2.6 — F1 (GR-clone) preliminary

For F1 (b_1=-2, b_2=4, b_3=-9):
- Newton: c_1 = 1 (from -2/-2)
- 1PN: c_2 = 0 (from 2/-2 - 2·4/(-8) = -1 + 1 = 0)
- For α_3 = -3/2 (GR match): c_3_required = 0
- **Whether EOM gives c_3 = 0 is the Phase 3 question**

## §3 — Phase 2 status

✅ Framework consistency: confirmed sympy-rigorously  
✅ R3 ODE f-independence: identified  
✅ Algebraic constraints C2: derived (Newton + 1PN)  
🔄 Δα_3 = 0 STRATEGY: requires EOM c_3 derivation per candidate (Phase 3)  
🔄 F1 priority: c_3_F1_required = 0; check if EOM gives this

## §4 — Cross-references

- [[Phase0_balance.md]] — constraint inventory
- [[Phase1_results.md]] (TODO if needed)
- [[Phase2_baseline_M911_reproduction_sympy.py]]
- [[Phase2_structural_obstruction_sympy.py]]
- [[Phase3_F1_attempt_sympy.py]] — Phase 3 F1 derivation attempt
- [[../op-g0-r3-from-canonical-projection/phase2_P21_vacuum_uniqueness.py]] — G.0 V_M911 reference
