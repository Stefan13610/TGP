---
title: "Phase 2 — Setup: constraint testing per candidate"
date: 2026-05-09
parent: "[[README.md]]"
type: phase2-setup
phase: 2
---

# Phase 2 — Setup: constraint testing per candidate

## §0 — Goal

For each Phase 1 candidate, test compatibility with hard constraints
C1-C10 via sympy-rigorous derivation of Phi-EOM coupling and PN
expansion of f(ψ(U)).

## §1 — Critical observation post-Phase 1

The full Phi-EOM derivation requires:
1. Choice of h(ψ) ansatz (independent of f, OR coupled via f·h=1)
2. K(ψ) = ψ^4 (T-D-uniqueness, preserved)
3. V_grav(ψ) derived from R3 ODE projection (G.0-style)
4. Static spherical solution ψ(r) = 1 + ε(r) where ε(U) Taylor coefs
   come from EOM

**Then α_n = (Taylor coef of f(ψ(U)) at U^n)**, computed by composing
the Taylor series.

## §2 — Phase 2 scope decision

Full Phi-EOM derivation per candidate is **multi-session work**.
Pragmatic Phase 2 scope:

1. **Reproduce M9.1'' baseline**: derive Phi-EOM solution + α_n^M911
   sympy-rigorously to validate framework (target: α_3 = -7/3, Δα_3 = -5/6).

2. **F1 GR-exact-clone test**: attempt to derive Phi-EOM with:
   - f_F1(ψ) = (3-ψ)²/(1+ψ)²
   - h_F1(ψ) = 1/f_F1 = (1+ψ)²/(3-ψ)²  (anti-podal preserve)
   - K(ψ) = ψ⁴
   - V_grav derived from R3 ODE projection
   Goal: check if Δα_3 = 0 emerges (would be "GR clone" verification).

3. **Structural-obstruction analysis**: identify if ANY alternative f(ψ)
   in the (M9.1''-class) family can satisfy:
   - C1 (α=2 vacuum)
   - C2 (1PN exact GR)
   - C3 (|Δα_3·G_SPA| ≤ 8.32)
   - C9 (sqrt(-g) consistency)
   - R3 ODE projection (G.0 universal V_grav)
   simultaneously.

4. **Honest verdict**: classify outcome:
   - DERIVED: F1 (or other) passes all C1-C10 with first-principles derivation
   - STRUCTURAL_CONDITIONAL: passes some C, fails others (deferred to Phase 3)
   - EARLY_HALT: structural obstruction proven (M9.1''-class FALSIFIED)

## §3 — Anti-pattern compliance

- NIE post-hoc selection: F1 chosen as "most natural GR-mimic", justified
  ex ante (canonical isotropic Schwarzschild form).
- NIE drift hardening: V_grav derivation strictly z R3 ODE projection
  (NIE empirical fit to avoid Δα_3 ≠ 0).
- NIE constructed criterion: C1-C10 pre-declared in Phase 0.

## §4 — Phase 2 deliverables

- Phase2_baseline_M911_reproduction_sympy.py
- Phase2_F1_GR_clone_test_sympy.py
- Phase2_structural_obstruction_sympy.py
- Phase2_results.md

## §5 — Gate criteria

| # | Criterion |
|---|---|
| 1 | M9.1'' baseline α_n reproduced sympy (4-LEVEL VERIFICATION) |
| 2 | F1 attempted: Δα_3 derived OR structural obstruction shown |
| 3 | If F1 negative: structural argument extends to broader class |
| 4 | Honest verdict z classification |
| 5 | All anti-patterns checked |

## §6 — Cross-references

- [[../op-g0-r3-from-canonical-projection/phase2_P21_vacuum_uniqueness.py]] — G.0 V_M911 derivation reference
- [[../op-ppE-mapping/Phase1.5_G_SPA_lock.md]] — α_n^M911 LOCK 7/7
