---
title: "σ.1 program — substrate-light dispersion (c-mechanism)"
date: 2026-04-30
cycle: σ.1
status: ACTIVE
parent: "[[../op-omega1-substrate-em-coupling/Phase3_results.md]]"
notation_patch: "2026-05-01 dimensional + sign + framing critique applied"
tags:
  - TGP
  - sigma1
  - c-mechanism
  - dispersion
  - substrate-light
---

# σ.1 program

> **NOTE 2026-05-01 (σ.1 critique patch):** Notation harmonized after technical critique:
> (i) dispersion written in dimensionally-consistent form `ω_±² = k² ± g·k·n_∥`
> with `n_∥ ≡ k̂·∇ln X` (covariant: `p_μ = g ∂_μ ln X`, special case `ω² = k² ± k p_∥`),
> NOT `ω² = k² ± g(n·k̂)` (which would mix dim 1/L² and 1/L·1/L = 1/L²
> only if `(n·k̂)` is read as `n_∥` — the explicit form removes ambiguity);
> (ii) group-velocity sign corrected: `v_g,± ≈ 1 + (gn_∥)²/(8k²)` (positive,
> earlier draft erroneously wrote `1 − …`); (iii) "effective optical metric LOCKED"
> softened to **helicity-dependent optical cone / effective dispersion geometry**
> — the WKB result is a leading null-cone for static spatial gradient `n_i ∥ k_i`,
> NOT a full classical pseudo-Riemannian metric (full CFJ-class structure
> `(k²)² + p²k² − (p·k)² = 0` is more general and beyond the present derivation).
> All physical conclusions (linear phase birefringence, no scalar c(X), Δχ rotation)
> are unchanged.

> **Goal:** From modified Maxwell ω.1 EOM derive **dispersion relation
> ω²(k, ∂ ln X)** for plane-wave EM in slowly-varying substrate gradient,
> giving **polarization-dependent phase/group velocity** v_±(X) =
> 1 ± (g/(2k))·n_∥ (where `n_∥ = k̂·∇ln X`) i emergent **helicity-dependent
> optical cone** (Carroll-Field-Jackiw-class effective dispersion geometry).
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
$$v_{\varphi,\pm}(k, n_\parallel) = 1 \pm \frac{g\,n_\parallel}{2k},\qquad n_\parallel \equiv \hat k\cdot\nabla\ln X$$

Average velocity (polarization-averaged): `c_eff = 1` at this order
(consistent z Webb/Murphy NULL in α_em).

Difference: birefringence Δv_φ = v_+ − v_− = g·n_∥/k → integrates
over path do polarization rotation Δχ = (g/2)∫n_∥ ds = (g/2)∫d(ln X)
(LIVE PARTIAL hint w CMB Planck PR4 + ACT 2024, ~3.8σ candidate).

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
