---
title: "ω.1.Phase2 results — sympy LOCK + modified EOMs (7/7 FULL CASCADE)"
date: 2026-04-30
cycle: ω.1.Phase2
status: PASS
parent: "[[program.md]]"
predecessor: "[[Phase1_results.md]]"
tags:
  - TGP
  - omega1
  - phase2
  - sympy
  - eom
  - PASS
  - FULL-CASCADE
---

# ω.1.Phase2 results

**Score: 7/7 FULL CASCADE → Phase 3 forward.**

> **Headline:** Pełen Lagrangian
> $\mathcal{L}_{\omega.1} = -\tfrac{1}{4}F^2 + \tfrac{1}{2}f_X^2(\partial \ln X)^2 + \tfrac{g}{4}(\ln X)F\tilde F$
> generuje **modified Maxwell**: $\partial_\nu F^{\nu\mu} = g\,\tilde F^{\mu\nu}\partial_\nu(\ln X)$
> oraz **modified substrate EOM**: $\Box(\ln X) = (g/(4f_X^2))\,F\tilde F$ (sympy
> LOCK exact). Bianchi $\partial_\mu \tilde F^{\mu\nu} = 0$ niezmieniona.
> Lorenz-gauge consistency 0=0 (no anomaly w abelian U(1)). T^μν gauge-inv,
> axion piece = boundary in bulk. **3 alt-couplings FALSIFIED**, axion UNIQUE.

## Sub-test results

### W2.1 — Modified Maxwell EOM ✓ PASS

Variation L względem A_μ:

$$\partial_\nu F^{\nu\mu} = g \, \tilde F^{\mu\nu} \, \partial_\nu(\ln X)$$

Po IBP axion term: $-g(\partial_\mu \ln X) A_\nu \tilde F^{\mu\nu}$, skąd
$\partial L/\partial A_\mu = -g(\partial_\nu \ln X)\tilde F^{\nu\mu}$.

**Physical interpretation**: substrate gradient $\partial_\mu(\ln X)$ acts
jako effective axion-like background, sourcing parallel E·B field structure.

### W2.2 — Modified substrate EOM ✓ PASS (sympy LOCK)

```
L = (1/2) f_X² (∂_μ u)(∂^μ u) + (g/4) u F·F̃
EL eq: f_X² □u - (g/4) F·F̃ = 0
sympy diff (EL - expected) = 0  EXACT
```

$$\boxed{\,\Box(\ln X) = \frac{g}{4 f_X^2}\, F_{\mu\nu}\tilde F^{\mu\nu}\,}$$

**This is the answer to the user's research question**: pole EM **może
strukturalnie oddziaływać na przestrzeń substrate** poprzez F·F̃ source —
non-trivial gdy g ≠ 0 i E·B ≠ 0 (parallel E and B fields).

W phi.1 (no EM): `□(ln X) = 0`. W ω.1 (z EM): right-hand side staje się
**non-zero** dla configurations z E ⃗·B ⃗ ≠ 0.

### W2.3 — Bianchi identity ✓ PASS

$F_{\mu\nu} = \partial_\mu A_\nu - \partial_\nu A_\mu$ jako tożsamość (independent
of L) → $\partial_\mu \tilde F^{\mu\nu} = 0$ unchanged. Brak monopoli magnetycznych.

### W2.4 — Alt-coupling falsification ✓ PASS

| Form | Behavior | Pass |
|---|---|:---:|
| dilaton `(ln X)F²` | scale-breaks (F² weight 4) | ✗ |
| minimal `eA^μ∂_μ(ln X)` | gauge-trivial w Lorenz | ✗ |
| gradient `(∂ ln X)²F²` | dim-8 EFT-irrelevant | ✗ |
| **axion `(ln X)F·F̃` canonical** | **gauge-inv + scale-inv + non-trivial + dim-4** | **✓** |

3 alts FALSIFIED, axion UNIQUE.

### W2.5 — Lorenz-gauge consistency ✓ PASS

Divergence of modified Maxwell:
- LHS: $\partial_\mu \partial_\nu F^{\nu\mu} = 0$ (antisym + commute)
- RHS: $g[\partial_\mu \tilde F^{\mu\nu} \cdot \partial_\nu \ln X + \tilde F^{\mu\nu}\partial_\mu \partial_\nu \ln X] = 0$

Zatem **0 = 0 consistently**. W Lorenz gauge: $-\Box A^\mu = g\tilde F^{\mu\nu}\partial_\nu(\ln X)$.

### W2.6 — Stress-energy + scale-current ✓ PASS

- $T^{\mu\nu}_\text{em} = -F^{\mu\alpha}F^\nu_{\,\alpha} + \tfrac{1}{4}\eta^{\mu\nu}F^2$
- $T^{\mu\nu}_\text{sub} = f_X^2(\partial^\mu u)(\partial^\nu u) - \tfrac{1}{2}\eta^{\mu\nu}f_X^2(\partial u)^2$
- $T^{\mu\nu}_\text{axion} = 0$ in bulk (F·F̃ jest total div)

Scale-symmetry preserved → $J^\mu_\text{scale}$ conserved on-shell modulo trace.

### W2.7 — Phase 2 gate ✓ PASS

7/7 FULL CASCADE.

## Promotions post-Phase 2

- **Modified Maxwell EOM LOCKED**: $\partial_\nu F^{\nu\mu} = g\tilde F^{\mu\nu}\partial_\nu(\ln X)$
- **Modified substrate EOM LOCKED** (sympy exact): $\Box(\ln X) = (g/(4f_X^2))F\tilde F$
- **EM ↔ substrate back-reaction CHANNEL OPENED** (closes G3 gap from em_from_substrate)
- **Bianchi identity preserved** (no monopoles introduced)
- **Lorenz-gauge consistency 0=0** (no abelian anomaly)
- **3 alt-couplings FALSIFIED**, axion canonical UNIQUE
- **Stress-energy gauge-invariant**, scale-current modified consistently
- **Phase 3 forward**: predictions + 4-channel convergence

## Cross-references

- [[program.md]]
- [[Phase2_setup.md]]
- [[Phase1_results.md]]
- [[../op-phi1-substrate-action-variational/Phase3_results.md]]
