---
title: "ω.1.Phase1 results — coupling-form scan + structural invariance (5/5 FULL PASS)"
date: 2026-04-30
cycle: ω.1.Phase1
status: PASS
parent: "[[program.md]]"
predecessor: "[[Phase1_setup.md]]"
tags:
  - TGP
  - omega1
  - phase1
  - coupling-scan
  - PASS
  - FULL-PASS
---

# ω.1.Phase1 results

**Score: 5/5 FULL PASS → Phase 2 forward.**

> **Headline:** Wśród 4 kandydujących form coupling EM ↔ substrate **jedynie
> axion-like topological coupling** $\tfrac{g}{4}(\ln X) F_{\mu\nu}\tilde F^{\mu\nu}$
> spełnia jednocześnie 4 kryteria strukturalne: gauge invariance, scale
> invariance pod X→λX, non-triviality, dim-4 EFT relevance. Naive minimal
> coupling $A^\mu J_\mu$ z $J_\mu = \partial_\mu \ln X$ jest gauge-trivial
> (znika w Lorenz gauge). Dilaton $(\ln X)F^2$ łamie skali, gradient
> $(\partial \ln X)^2 F^2$ jest dim-8 EFT-irrelevant.

## Sub-test results

### W1.1 — Coupling-form scan (4 candidates) ✓ PASS

| Candidate | Form | Gauge | Scale | NonTriv | Dim-4 | Pass |
|---|---|:---:|:---:|:---:|:---:|:---:|
| minimal A·J | `e A^μ ∂_μ(ln X)` | Y | Y | **n** (gauge-trivial) | Y | ✗ |
| dilaton | `(g/4)(ln X)F²` | Y | **n** | Y | Y | ✗ |
| **axion** | **`(g/4)(ln X)FF̃`** | **Y** | **Y** | **Y** | **Y** | **✓** |
| gradient | `(g/4)(∂ ln X)²F²` | Y | Y | Y | **n** (dim-8) | ✗ |

**Unique winner**: axion-like topological coupling.

### W1.2 — Gauge invariance ✓ PASS

Identity: $F_{\mu\nu}\tilde F^{\mu\nu} = 4\,\partial_\mu(A_\nu \tilde F^{\mu\nu})$ — total divergence.

Pod $A_\mu \to A_\mu + \partial_\mu \Lambda$: $\delta F_{\mu\nu} = 0$ → $\delta(F\tilde F) = 0$.
Manifestly gauge-invariant. Integration by parts daje Chern-Simons-like formę
z $\partial_\mu(\ln X)$ jako background gauge field.

### W1.3 — Scale invariance X → λX ✓ PASS

Pod $X \to \lambda X$: $\ln X \to \ln X + \ln \lambda$, $\partial_\mu \ln X$ unchanged.

Action shift: $\delta S = (g/4) \ln\lambda \int F\tilde F \,d^4x$.

Dla finite-energy U(1) (F → 0 at ∞): $\int F\tilde F \,d^4x = 0$ → **scale-invariant exactly**.

Dla topologically non-trivial uplift (non-abelian): $\delta S$ trivial mod 2π
(Witten quantization).

### W1.4 — g LOCK candidates ✓ PASS

[g] = 0 (dimensionless). Kandydaci z TGP framework:

| Constant | Value | Source |
|---|---:|---|
| κ_TGP | 2.0120 | XS.1 cross-sector SC charge |
| α₀ = κ²_TGP | 4.0180 | BH photon-ring |
| α_em | 7.297×10⁻³ | fine structure |
| η_chir = 19/24 | 0.7917 | chirality factor |
| 1/(2π) | 0.1592 | axion EFT norm |

Phase 2: g LOCK via cross-sector identity. Phase 3: g/f_a bound z PVLAS/CMB.

### W1.5 — Phase 1 gate ✓ PASS

5/5 FULL PASS.

## Promotions post-Phase 1

- **Axion-like topological coupling form UNIQUE-LOCKED** (4 alts FAIL: minimal trivial, dilaton breaks scale, gradient dim-8)
- **Lagrangian skeleton confirmed**: $\mathcal{L}_{\omega.1} = -\tfrac{1}{4}F^2 + \tfrac{1}{2}f_X^2(\partial \ln X)^2 + \tfrac{g}{4}(\ln X)F\tilde F$
- **Gauge + scale invariance** simultaneously preserved
- **g dimensionless**, LOCK candidates from TGP framework identified
- **Phase 2 forward**: sympy LOCK + modified EOM derivation

## Cross-references

- [[program.md]]
- [[Phase1_setup.md]]
- [[../op-phi1-substrate-action-variational/Phase3_results.md]]
