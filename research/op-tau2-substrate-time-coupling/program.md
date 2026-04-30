---
title: "τ.2 program — substrate-time coupling (clock rates pod ln X gradient)"
date: 2026-04-30
cycle: τ.2
status: ACTIVE
parent: "[[../op-sigma1-substrate-light-dispersion/Phase3_results.md]]"
tags:
  - TGP
  - tau2
  - clock
  - substrate-time
  - scale-invariance
---

# τ.2 program

> **Goal:** Determine how atomic clock rates depend na substrate gradient
> ∂_μ(ln X). Derive **scale-symmetry protection theorem**: φ.1 X→λX gauge
> wymaga atomic masses X-invariant → leading-order clock-rate drift = 0
> (consistent z Webb/Murphy α_em NULL). Identify **new τ.2 signatures**:
> polarization-dependent Zeeman shifts (connects σ.1 birefringence),
> strong-gradient residuals (magnetar/lab E·B), cosmological clock
> comparison NULL.

## Foundation (post-σ.1)

φ.1 substrate-action: $\mathcal{L}_{\phi.1} = \tfrac{1}{2}f_X^2(\partial_\mu \ln X)(\partial^\mu \ln X)$
→ X→λX gauge symmetry, Noether scale-current J^μ = ∂^μ(ln X).

ω.1 axion-like EM coupling: $\mathcal{L}_{\omega.1} \supset (g/4)(\ln X)F\tilde F$
→ modified Maxwell + modified substrate EOMs.

σ.1 dispersion result: photon polarization-dependent v_φ_± = 1 ± gn/(2k);
NO scalar c(X) at leading O(gn/k); polarization-averaged c_eff = 1.

τ.2 task: extend do MASSIVE matter (atoms, clocks).

## Critical scale-protection insight

**Scale invariance X → λX wymaga**:
- Massless scalars: protected (kinetic ½(∂ ln X)²)
- Massive fermions/scalars: m must be X-INDEPENDENT lub scale-breaking
- Atomic masses (m_e, m_p) → X-INVARIANT z φ.1 axiom
- Atomic transition frequencies ∝ m_e c² α_em² → X-INVARIANT
- **Clock rates protected at leading order O(∂(ln X))**

This generalizes σ.1 scalar c(X) protection do all scalar observables.

## Possible NEW τ.2 signatures (beyond protection)

1. **Polarization-Zeeman cross-coupling**: σ.1 birefringence affects
   Zeeman levels in atomic clocks driven by circularly polarized light;
   ground-state hyperfine splitting in B field → δω ∝ g·∂(ln X)/k_drive.

2. **Strong-gradient regime**: magnetar pole regions z E·B ≠ 0 source
   `□(ln X) ≠ 0` (ω.1 W3.2) → local atomic clocks experience gradient.

3. **Gravitational-time analog**: substrate gradient acts as effective
   gravitational potential? Open question — possibly higher-order.

## Phase plan (5+7+6 = 18 sub-tests)

### Phase 1 — scale-protection theorem (5 tests)

- **T1.1** Atomic-mass coupling candidates scan (4 forms)
- **T1.2** Scale-invariance X→λX requires m_atom INVARIANT
- **T1.3** Noether scale-current implication for proper time
- **T1.4** Effective Hamiltonian H_atom X-independent
- **T1.5** Protection theorem: clock rate ∝ frequency standard X-invariant

### Phase 2 — proper-time formalism + sympy LOCK (7 tests)

- **T2.1** Proper time τ = ∫√(g_00^eff) dt sympy LOCK
- **T2.2** Substrate-induced time dilation analog δτ/τ at leading order
- **T2.3** Polarization-Zeeman cross-coupling z σ.1
- **T2.4** Cs hyperfine vs optical Sr/Yb relativistic protection
- **T2.5** Clock-rate ratio R(X1)/R(X2) sympy LOCK
- **T2.6** ω.1 + σ.1 consistency cross-check
- **T2.7** Alt-couplings (m_e ∝ X^α, ℏ ∝ X^β) FALSIFIED

### Phase 3 — predictions + 4-channel convergence (6 tests)

- **T3.1** Atomic clock cosmological drift NULL (Webb/Murphy)
- **T3.2** Lab Hg/Yb/Sr clock comparison precision
- **T3.3** Strong-gradient residuals (magnetar, lab E·B)
- **T3.4** Polarization-Zeeman cross-coupling z σ.1
- **T3.5** Alt-clock-coupling cross-channel falsification
- **T3.6** 4-channel τ.2 convergence

## Cross-references

- [[../op-sigma1-substrate-light-dispersion/Phase3_results.md]] — direct predecessor
- [[../op-omega1-substrate-em-coupling/Phase3_results.md]] — EM coupling parent
- [[../op-phi1-substrate-action-variational/Phase3_results.md]] — substrate action axiom
- [[../../INDEX.md]]
- [[../../PREDICTIONS_REGISTRY.md]]
