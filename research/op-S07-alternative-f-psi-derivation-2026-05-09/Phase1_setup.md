---
title: "Phase 1 — Setup: candidate f(ψ) family enumeration"
date: 2026-05-09
parent: "[[README.md]]"
type: phase1-setup
phase: 1
---

# Phase 1 — Setup: candidate f(ψ) family enumeration

## §0 — Goal

Enumerate ≥5 candidate parametrizations of f(ψ) (gravity sector ansatz)
that satisfy:
- f(1) = 1 (vacuum at infinity, GR matching at r → ∞)
- Symbolic structure admitting alpha_1 = -2 condition (1PN matching)
- Φ-EOM solvability (analytical or numerical)

## §1 — Enumeration plan

### Family P (polynomial)

| Subfamily | Form | Free params | Notes |
|---|---|---|---|
| P_2 | f(ψ) = 1 + b_1(ψ-1) + b_2(ψ-1)² | b_1, b_2 | Min-degree; loses BH horizon |
| P_3 | f(ψ) = 1 + b_1(ψ-1) + b_2(ψ-1)² + b_3(ψ-1)³ | 3 | Degree-3, more freedom |
| P_4 | similar to ψ⁴ | 4 | Higher degree |

### Family R (rational)

| Subfamily | Form | Free params | Notes |
|---|---|---|---|
| R_(1,1) | f(ψ) = (a₀+a₁ψ)/(b₀+b₁ψ) | 3 (modulo scale) | M9.1'' (4-3ψ)/ψ FALSIFIED |
| R_(2,1) | f(ψ) = (a₀+a₁ψ+a₂ψ²)/(b₀+b₁ψ) | 4 | More freedom |
| R_(2,2) | f(ψ) = (a₀+...)/(b₀+...) | 5 | Quad/quad |

### Family E (exponential)

| Subfamily | Form | Free params | Notes |
|---|---|---|---|
| E_lin | f(ψ) = exp(-2(ψ-1)) | 0 (after match) | Pure GR-style isotropic |
| E_quad | f(ψ) = exp(-2(ψ-1) + b(ψ-1)²) | 1 | 2PN tunable |
| E_general | f(ψ) = exp(g(ψ)) for polynomial g | N | General |

### Family GR-exact-isotropic (special case)

| Subfamily | Form | Free params | Notes |
|---|---|---|---|
| GR_iso | f(ψ(U)) = ((1-U/2)/(1+U/2))² (forced by ψ(U) ansatz) | 0 | Pure GR clone |

This is the "trivial pass" candidate — f(ψ) chosen such that g_tt is
EXACTLY the GR isotropic Schwarzschild form. Then by construction Δα_3 = 0
at all orders. Question: does this satisfy other constraints (mass spectrum,
hyperbolicity, etc.)?

### Family C (constraint-natural / variational)

f(ψ) forced by C1-C10 algebraic system (constraint solver, NOT empirical
fit). Phase 3 territory.

## §2 — Phase 1 sympy script structure

```
For each candidate family:
  1. Define f_candidate(psi) symbolically (with free params if any)
  2. Check f_candidate(1) = 1 (anchor)
  3. Compute b_1 = f'(1), b_2 = f''(1), b_3 = f'''(1)
  4. Identify constraints on free params for:
     - alpha_1 = -2 (with c_1 from Phi-EOM)
     - alpha_2 = 2 (with c_1, c_2 from Phi-EOM)
     - alpha_3 = -3/2 (Δα_3 = 0 strategy, optional)
  5. Optional: derive natural h(psi) from f-h coupling assumption
  6. Note: full constraint testing in Phase 2 (with Phi-EOM coupling).
```

## §3 — Phase 1 gate criteria

| # | Criterion |
|---|---|
| 1 | ≥5 candidate families enumerated |
| 2 | For each: f(1) = 1 verified (sympy) |
| 3 | For each: b_1, b_2, b_3 computed (sympy) |
| 4 | Algebraic constraints for C2 (1PN) identified |
| 5 | Anti-pattern check: NO post-hoc selection |
| 6 | Honest reporting of which candidates fail prelim |

## §4 — Cross-references

- [[Phase0_balance.md]] — constraint inventory
- [[Phase0_constraints_C1_C10_sympy.py]] — formalization
- [[../op-g0-r3-from-canonical-projection/phase2_P21_vacuum_uniqueness.py]] — G.0 reference
