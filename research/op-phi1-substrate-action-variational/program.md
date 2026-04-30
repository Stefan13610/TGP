---
title: "φ.1 program — closure substrate-action variational principle"
date: 2026-04-30
cycle: φ.1
status: ACTIVE
parent: "[[../op-upsilon1-closure-cross-family/Phase3_results.md]]"
tags:
  - TGP
  - phi1
  - substrate-action
  - variational-principle
  - lagrangian
  - axiom-lift
---

# φ.1 — closure substrate-action variational principle

**Goal:** Lift υ.1 universal closure law DERIVED → AXIOM-level via
explicit Lagrangian/action-functional derivation z X→λX gauge invariance
+ N_gen=3 cascade primality jako geometric/topological invariant w
substrate-action density.

## Hypothesis

The universal closure law from υ.1:
$$\text{closure}(X_{\text{ref}}, X_{\text{obs}}) = \left(\frac{X_{\text{ref}}}{X_{\text{obs}}}\right)^{1/N_{\text{gen}}}$$
with N_gen=3, is **not a postulate** but the unique extremum of a
substrate-action functional with:

1. **Action density:** $\mathcal{L} = \tfrac{1}{2} (\partial_\mu \ln X)(\partial^\mu \ln X)$
2. **Scale invariance:** X → λX is a Noether symmetry (∂_μ ln X invariant)
3. **Euler-Lagrange equation:** $\Box (\ln X) = 0$ (massless scalar in log-coord)
4. **N_gen=3 cascade subdivision:** boundary X(0) = X_ref, X(L) = X_obs, with
   intermediate sampling at fractional point x = L/N_gen yields
   $X(L/N_{\text{gen}}) = X_{\text{ref}} \cdot (X_{\text{obs}}/X_{\text{ref}})^{1/N_{\text{gen}}}$
5. **Closure factor:** $X_{\text{ref}} / X(L/N_{\text{gen}}) = (X_{\text{ref}}/X_{\text{obs}})^{1/N_{\text{gen}}}$

This **exactly reproduces** υ.1 sympy-DERIVED form, and lifts it to AXIOM-level
via single variational principle.

## Reference data

| Source | Closure form | υ.1 status | φ.1 status |
|---|---|:---:|:---:|
| π.1 NME | `M ∝ (76/A)^{1/3}` | DERIVED | AXIOM-LIFT |
| τ.1 overlap | `f = (Z_a/Z_t)^{1/3}` | DERIVED | AXIOM-LIFT |
| υ.1 universal | `closure(X) = (X_ref/X_obs)^{1/N_gen}` | DERIVED | AXIOM-LIFT |
| φ.1 Lagrangian | `L = ½(∂ ln X)²` | — | **NEW** |

## Phase plan (5 + 7 + 6 = 18 sub-tests)

### Phase 1 — Lagrangian ansatz scan + EL equation derivation (5 sub-tests)

- **F1.1** — Lagrangian candidate scan (4 candidates: log-derivative, power-law, mass-term, geometric)
- **F1.2** — Euler-Lagrange equation derivation z best ansatz → recover closure
- **F1.3** — Scale invariance X→λX as Noether symmetry → conserved current
- **F1.4** — N_gen=3 fractional sampling z 4-species cascade subdivision
- **F1.5** — gate ≥4/5 PASS

### Phase 2 — sympy LOCK + structural derivation (7 sub-tests)

- **F2.1** — sympy LOCK: L = ½(∂ ln X)² → EL → closure(X_ref, X_obs)
- **F2.2** — sympy LOCK: Noether current J^μ = ∂^μ(ln X) z X→λX symmetry
- **F2.3** — Alt-Lagrangian falsification scan (4 alt forms)
- **F2.4** — N_gen=3 z fermion-generation cascade primality embedded
- **F2.5** — Cross-family combined Lagrangian (π.1 + τ.1 + υ.1 unified)
- **F2.6** — Boundary-condition uniqueness (fixed X_ref, X_obs → unique solution)
- **F2.7** — gate ≥6/7 PASS

### Phase 3 — predictions + 4-channel convergence (6 sub-tests)

- **F3.1** — π.1 ↔ τ.1 ↔ υ.1 triple closure z φ.1 Lagrangian (POST-CONFIRM)
- **F3.2** — Cross-sector predictions (³⁷Ar, ⁷Be, ⁵¹Cr, ⁷¹Ga, ⁹⁸Mo, ¹³⁷Cs)
- **F3.3** — alt-action variational falsification (5 alt actions)
- **F3.4** — Ward identity / scale-current conservation experimental signature
- **F3.5** — 4-channel φ.1 convergence
- **F3.6** — Future N_gen extension hint (4-gen, 5-gen sterile sector)

## Cross-references

- [[../op-upsilon1-closure-cross-family/Phase3_results.md]] — υ.1 closure law DERIVED
- [[../op-pi1-bb0nu-nme-isotope/Phase3_results.md]] — π.1 NME 1/A^{1/3}
- [[../op-tau1-closure-overlap-coulomb/Phase3_results.md]] — τ.1 overlap 1/N_gen
- [[../../INDEX.md]]
- [[../../PREDICTIONS_REGISTRY.md]]
