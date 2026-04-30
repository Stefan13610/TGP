---
title: "σ.1 program — substrate-light dispersion (c-mechanism)"
date: 2026-04-30
cycle: σ.1
status: ACTIVE
parent: "[[../op-omega1-substrate-em-coupling/Phase3_results.md]]"
tags:
  - TGP
  - sigma1
  - c-mechanism
  - dispersion
  - substrate-light
---

# σ.1 program

> **Goal:** From modified Maxwell ω.1 EOM derive **dispersion relation
> ω²(k, ∂ ln X)** for plane-wave EM in slowly-varying substrate gradient,
> giving **polarization-dependent phase/group velocity** v_±(X) =
> 1 ± (g/(2k))·n̂·∂(ln X) i emergent **effective optical metric** g_μν^opt.
> Closes user's c-mechanism question: prędkość światła zależy od substrate
> gradient w sposób polarization-dependent (axion-photon mixing).

## Foundation (post-ω.1)

ω.1 LOCKED Lagrangian:
$$\mathcal{L}_{\omega.1} = -\tfrac{1}{4}F^2 + \tfrac{1}{2}f_X^2(\partial \ln X)^2 + \tfrac{g}{4}(\ln X)F\tilde F$$

Modified Maxwell EOM:
$$\partial_\nu F^{\nu\mu} = g\,\tilde F^{\mu\nu}\,\partial_\nu(\ln X)$$

σ.1 task: solve this for plane-wave A_μ ∝ e^{ik·x} in background ∂(ln X) ≠ 0.

## Critical insight — birefringence not scalar c(X)

In ω.1 axion-like coupling, scalar speed-of-light variation `c(X)` is
**ruled out** (would require dilaton `(ln X)F²` form, FALSIFIED in ω.1
Phase 1 by scale-invariance + Webb/Murphy α_em NULL).

What ω.1 DOES allow: **chirality-dependent phase velocity** for left/right
circular polarizations:
$$v_\pm(k, \partial \ln X) = 1 \pm \frac{g\,n̂\cdot\partial(\ln X)}{2k}$$

Average velocity (polarization-averaged): `c_eff = 1` at this order
(consistent z Webb/Murphy NULL in α_em).

Difference: birefringence Δv = v_+ - v_- = g·n̂·∂(ln X)/k → integrates
over path do polarization rotation Δχ (już zaobserwowane w CMB Planck PR4).

## Phase plan (5+7+6 = 18 sub-tests)

### Phase 1 — dispersion relation (5 tests)

- **W1.1** Plane-wave ansatz w modified Maxwell EOM
- **W1.2** Dispersion relation ω²(k, ∂ ln X) derivation
- **W1.3** Polarization eigen-modes (L/R circular)
- **W1.4** Slowly-varying gradient validity (WKB)
- **W1.5** Gauge structure preservation

### Phase 2 — phase/group velocity LOCK (7 tests)

- **W2.1** Phase velocity v_φ ± δv (sympy LOCK)
- **W2.2** Group velocity v_g = ∂ω/∂k
- **W2.3** Polarization-averaged effective c_eff
- **W2.4** Birefringence Δv = v_+ - v_-
- **W2.5** Effective optical metric g_μν^opt
- **W2.6** Scale-invariance preservation X→λX
- **W2.7** ω.1 EOM consistency (substrate back-reaction)

### Phase 3 — predictions + 4-channel convergence (6 tests)

- **W3.1** Lab phase-velocity measurement (Mach-Zehnder + B field)
- **W3.2** Pulsar polarized dispersion residuals
- **W3.3** CMB E/B-mode chirality (already ω.1 partial)
- **W3.4** Atomic clock ratio gradient sensitivity
- **W3.5** Alt-dispersion cross-channel falsification
- **W3.6** 4-channel σ.1 convergence

## Cross-references

- [[../op-omega1-substrate-em-coupling/Phase3_results.md]] — direct predecessor
- [[../op-phi1-substrate-action-variational/Phase3_results.md]] — substrate action
- [[../../INDEX.md]]
- [[../../PREDICTIONS_REGISTRY.md]]
