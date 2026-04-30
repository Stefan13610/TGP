---
title: "σ.1.Phase1 results — dispersion relation (5/5 FULL PASS)"
date: 2026-04-30
cycle: σ.1.Phase1
status: PASS
parent: "[[program.md]]"
predecessor: "[[Phase1_setup.md]]"
tags:
  - TGP
  - sigma1
  - phase1
  - dispersion
  - PASS
  - FULL-PASS
---

# σ.1.Phase1 results

**Score: 5/5 FULL PASS → Phase 2 forward.**

> **Headline:** Plane-wave EM in modified Maxwell ω.1 z substrate gradient
> n_μ = ∂_μ(ln X) generuje **dispersion relation**
> $\omega_\pm^2(k, n) = k^2 \pm g\,(n\!\cdot\!k)$
> z polarization eigen-modes = circular L/R helicities. WKB validity
> |n|/k << 1 satisfied across cosmologies + labs (n/k ~ 1e-28 dla CMB).
> Gauge structure preserved → 2 transverse physical d.o.f.

## Sub-test results

### W1.1 — Plane-wave ansatz ✓ PASS

A_μ = a_μ e^{ik·x} → F_{μν} = i(k_μ a_ν - k_ν a_μ), F̃^{μν} = i ε^{μνρσ} k_ρ a_σ.
Background ∂_μ(ln X) = n_μ (constant gradient). LHS modified Maxwell w Lorenz
gauge: -k² a^μ. RHS: i g ε^{μνρσ} n_ν k_ρ a_σ. Equation = 4×4 eigenval problem.

### W1.2 — Dispersion relation ✓ PASS

W helicity basis e_± = (e_x ± i e_y)/√2 (transverse, k = ẑ direction):

$$\boxed{\omega_\pm^2(k, n) = k^2 \pm g\,(n\cdot \hat k)}$$

WKB leading: ω_± ≈ k ± g(n·k̂)/2. Canonical axion-photon dispersion z parity-odd
splitting at O(g/k).

### W1.3 — Polarization eigen-modes (L/R circular) ✓ PASS

Eigenvectors of M(k,n):
- e_+ (right circular): λ_+ = -k² + g(n·k)
- e_- (left circular):  λ_- = -k² - g(n·k)

Linear polarization NIE jest eigen-mode; rotuje jako superpozycja e_+ + e_-
z różnymi fazami → axion-induced birefringence Δχ = (g/2) n·L.

### W1.4 — WKB gradient validity ✓ PASS

|∂(ln X)|/k << 1 satisfied:
- CMB photon: n/k ~ 1e-28 (cosmological H_0 vs k ~ 1e2 m^-1)
- Optical photon: n/k ~ 1e-33
- Lab gradient: n/k ~ 1e-37 (B^2 sourced via box(ln X))

WKB regime fully justified across all astrophysical + lab scales.

### W1.5 — Gauge structure preservation ✓ PASS

A_μ → A_μ + ∂_μ Λ leaves F_{μν}, F̃^{μν} invariant → modified Maxwell EOM
gauge-invariant. Lorenz / Coulomb gauges yield identical dispersion.
Residual gauge freedom a → a + α k preserves all eqs → 2 transverse d.o.f.

## Promotions post-Phase 1

- **σ.1 dispersion relation LOCKED**: ω_±² = k² ± g(n·k) z polarization eigen-modes L/R
- **WKB validity established** dla cosmologies + labs (n/k ≤ 1e-28)
- **Gauge invariance preserved** → 2 physical transverse polarizations
- **Phase 2 forward**: phase/group velocity sympy LOCK + optical metric

## Cross-references

- [[program.md]]
- [[Phase1_setup.md]]
- [[../op-omega1-substrate-em-coupling/Phase2_results.md]] — modified Maxwell parent
- [[../op-omega1-substrate-em-coupling/Phase3_results.md]] — ω.1 closure
