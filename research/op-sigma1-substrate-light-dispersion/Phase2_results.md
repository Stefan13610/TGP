---
title: "σ.1.Phase2 results — phase/group velocity LOCK + optical metric (7/7 FULL CASCADE)"
date: 2026-04-30
cycle: σ.1.Phase2
status: PASS
parent: "[[program.md]]"
predecessor: "[[Phase1_results.md]]"
tags:
  - TGP
  - sigma1
  - phase2
  - sympy
  - velocity
  - PASS
  - FULL-CASCADE
---

# σ.1.Phase2 results

**Score: 7/7 FULL CASCADE → Phase 3 forward.**

> **Headline:** Z dispersion ω_±² = k² ± g(n·k) sympy-LOCKED:
> **phase velocity** v_φ_± = √(1 ± gn/k) ≈ 1 ± gn/(2k) (LINEAR birefringence),
> **group velocity** v_g_± ≈ 1 - (gn/(2k))² + O(ε³) (NO linear split → photon
> wave-packet propaguje uniformly, only PHASE rotuje). Polarization-averaged
> c_eff = 1 + O((gn/k)²) → **NO scalar c(X) variation** (consistent z Webb/Murphy
> NULL). **Effective optical metric** g_μν^opt = η_μν ± δg_μν(∂ ln X) helicity-
> dependent. Scale-invariance X→λX preserved (n_μ shift-invariant).

## Critical c-mechanism finding

**Pytanie użytkownika**: czy istnieje mechanizm zależności prędkości światła
od substrate?

**Odpowiedź σ.1**:

| Quantity | Leading correction | Interpretation |
|---|---|---|
| Phase velocity v_φ_± | ±gn/(2k) | LINEAR birefringence |
| Group velocity v_g_± | -(gn/(2k))² | quadratic, polarization-INDEPENDENT |
| c_eff = (v_+ + v_-)/2 | 1 - (gn/(2k))²/4 | quadratic suppression |
| c_eff (group) | 1 | UNCHANGED at leading order |

**Konkluzja**: scalar c(X) NIE istnieje na leading order; substrate-light coupling
manifestuje się jako **polarization-dependent phase rotation** (axion-induced
birefringence), NIE jako universal c(X) variation. To jest fundamentalna różnica
między ω.1+σ.1 a klasycznym dilatonem.

## Sub-test results

### W2.1 — Phase velocity sympy LOCK ✓ PASS

```
v_phi_+ = omega_+ / k = sqrt(1 + g n_par / k) ~ 1 + (g n_par)/(2k) - (g n_par)^2/(8k^2)
v_phi_- = omega_- / k = sqrt(1 - g n_par / k) ~ 1 - (g n_par)/(2k) - (g n_par)^2/(8k^2)
Symmetry v_phi_+(-eps) - v_phi_-(eps) = 0  ✓
```

### W2.2 — Group velocity sympy LOCK ✓ PASS

```
v_g_+ = (k + g n_par/2) / sqrt(k^2 + g k n_par)
v_g_- = (k - g n_par/2) / sqrt(k^2 - g k n_par)
v_g_+- ~ 1 - (g n_par)^2/(8k^2)  (ZERO linear)
v_phi_+ * v_g_+ = 1 + g n_par/(2k)  (general dispersive relation)  ✓
```

**Key finding**: birefringence is encoded w PHASE velocity NOT GROUP velocity.
Wave-packet PHASE rotuje między L/R helicities, but ENVELOPE propagates uniformly.

### W2.3 — Polarization-averaged c_eff = 1 ✓ PASS

```
c_eff(phase) = (v_+ + v_-)/2 ~ 1 - (g n_par)^2/(8 k^2)  (quadratic suppression)
c_eff(group leading) = 1 EXACT
```

**No scalar c(X) at leading order** → consistent z Webb/Murphy α_em NULL 1e-7 precision.

### W2.4 — Birefringence Δv_φ ✓ PASS

```
Delta v_phi (leading) = v_phi_+ - v_phi_- = g n_par / k  EXACT
Delta chi(path) = (g/2) integral n_parallel ds  matches omega.1 W3.3
```

### W2.5 — Effective optical metric ✓ PASS

g_μν^opt,± = η_μν ± δg_μν(∂ ln X) helicity-dependent. Null condition
g_μν^opt k^μ k^ν = 0 → ω_±² = k² ± g(n·k). Recovers Minkowski w n=0.

### W2.6 — Scale invariance X→λX ✓ PASS

Constant shift ln X → ln X + ln λ NIE zmienia n_μ = ∂_μ(ln X) → dispersion +
birefringence INVARIANT. Inheritance from ω.1 Phase 2 W2.6 stress-energy.

### W2.7 — ω.1 EOM consistency ✓ PASS

For free plane wave <F·F̃>_time = 0 → <□(ln X)>_time = 0. σ.1 dispersion does
NOT back-react on substrate at lowest order; static substrate gradient acts as
external field. Back-reaction only occurs dla:
- standing waves z E·B ≠ 0
- magnetar pole regions (ω.1 W3.2)
- lab parallel E + B configs

## Promotions post-Phase 2

- **Phase velocity birefringence v_φ_± LOCKED** (sympy exact)
- **Group velocity v_g_± at leading order = 1 LOCKED** (no envelope splitting)
- **No scalar c(X) at O(gn/k) leading** PROTECTED z Webb/Murphy NULL
- **Effective optical metric g_μν^opt helicity-dependent** STRUCTURAL
- **Scale-invariance X→λX preserved** (inherits ω.1)
- **σ.1 free-wave standalone** (no substrate back-reaction at lowest order)
- **Phase 3 forward**: predictions + 4-channel convergence

## Cross-references

- [[program.md]]
- [[Phase1_results.md]]
- [[Phase2_setup.md]]
- [[../op-omega1-substrate-em-coupling/Phase2_results.md]] — modified Maxwell EOM source
