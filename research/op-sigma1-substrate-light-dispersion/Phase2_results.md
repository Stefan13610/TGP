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

> **NOTE 2026-05-01 (σ.1 critique patch):** Three corrections vs. earlier draft:
> (i) dispersion written in dimensionally-explicit form `ω_±² = k² ± g·k·n_∥`
> (with `n_∥ ≡ k̂·∇ln X`), replacing ambiguous `ω² = k² ± g(n·k)`;
> (ii) **group velocity sign FLIPPED to POSITIVE** — direct sympy + analytic
> `dω/dk` give `v_g,± = 1 + (gn_∥)²/(8k²) + O(ε³)` (NOT minus); earlier draft
> printed wrong sign on the leading O(ε²) correction (the underlying sympy series
> result is correct after the corrected print);
> (iii) "effective optical metric LOCKED" softened to **helicity-dependent optical
> cone / effective dispersion geometry** — the WKB result is a leading null-cone
> for the special case of static spatial gradient `n_i ∥ k_i`, NOT a full
> classical pseudo-Riemannian metric (full CFJ structure
> `(k²)² + p²k² − (p·k)² = 0` is more general and beyond this derivation).
> Physical conclusions (no scalar c(X), linear v_φ birefringence, Δχ) unchanged.

> **Headline:** Z dispersion ω_±² = k² ± g·k·n_∥ sympy-LOCKED:
> **phase velocity** v_φ_± = √(1 ± gn_∥/k) ≈ 1 ± gn_∥/(2k) (LINEAR birefringence),
> **group velocity** v_g_± = 1 + (gn_∥)²/(8k²) + O(ε³) (NO linear split → photon
> wave-packet propaguje uniformly, only PHASE rotuje; sign POSITIVE, sympy-verified).
> Polarization-averaged
> c_eff = 1 − O((gn_∥/k)²)/8 → **NO scalar c(X) variation at leading order**
> (consistent z Webb/Murphy NULL). **Helicity-dependent optical cone**
> Q^(±)_μν(p) k^μ k^ν = 0 (with p_μ = g ∂_μ ln X) — leading WKB null-cone, NOT
> full classical metric. Scale-invariance X→λX preserved (n_μ shift-invariant).

## Critical c-mechanism finding

**Pytanie użytkownika**: czy istnieje mechanizm zależności prędkości światła
od substrate?

**Odpowiedź σ.1**:

| Quantity | Leading correction | Interpretation |
|---|---|---|
| Phase velocity v_φ_± | ±g·n_∥/(2k) | LINEAR birefringence |
| Group velocity v_g_± | +(g·n_∥)²/(8k²) | quadratic POSITIVE, polarization-INDEPENDENT (sympy-verified) |
| c_eff(phase) = (v_+ + v_-)/2 | 1 − (g·n_∥)²/(8k²) | quadratic suppression in PHASE avg |
| c_eff(group) = (v_g+ + v_g−)/2 | 1 + (g·n_∥)²/(8k²) | UNCHANGED at leading order in linear, slight super-luminal at O(ε²) |

**Konkluzja**: scalar c(X) NIE istnieje na leading order; substrate-light coupling
manifestuje się jako **polarization-dependent phase rotation** (axion-induced
birefringence), NIE jako universal c(X) variation. To jest fundamentalna różnica
między ω.1+σ.1 a klasycznym dilatonem.

## Sub-test results

### W2.1 — Phase velocity sympy LOCK ✓ PASS

```
Dispersion: omega_+-^2 = k^2 +- g k n_par,  n_par = k_hat . grad(ln X)
v_phi_+ = omega_+ / k = sqrt(1 + g n_par / k) ~ 1 + (g n_par)/(2k) - (g n_par)^2/(8k^2)
v_phi_- = omega_- / k = sqrt(1 - g n_par / k) ~ 1 - (g n_par)/(2k) - (g n_par)^2/(8k^2)
Symmetry v_phi_+(-eps) - v_phi_-(eps) = 0  ✓
```

### W2.2 — Group velocity sympy LOCK ✓ PASS  *(sign-corrected 2026-05-01)*

```
v_g_+ = (k + g n_par/2) / sqrt(k^2 + g k n_par)
v_g_- = (k - g n_par/2) / sqrt(k^2 - g k n_par)
sympy series at n_par=0 to O(n_par^2):
  v_g_+ to O(n_par^2):  1 + g^2 n_par^2 / (8 k^2)   [POSITIVE]
  v_g_- to O(n_par^2):  1 + g^2 n_par^2 / (8 k^2)   [POSITIVE]
  Coeff(v_g_+) - g^2/(8k^2) = 0   ✓
  Coeff(v_g_-) - g^2/(8k^2) = 0   ✓
  v_g_+- ~ 1 + (g n_par)^2/(8 k^2)  (ZERO linear, POSITIVE quadratic)
Direct check: v_g = d omega/dk = (2k +- g n_par)/(2 omega) = 1 + (g n_par)^2/(8 k^2) + O(eps^3)
v_phi_+ * v_g_+ = 1 + g n_par/(2k)  (general dispersive relation)  ✓
```

**Sign clarification**: earlier draft printed `v_g ~ 1 − (gn)²/(8k²)`; the
actual sympy series + analytic `dω/dk` both give POSITIVE sign. The cancellation
of LINEAR birefringence in `v_g` (Δv_g = 0 at O(ε)) is unchanged — only the
sub-leading O(ε²) sign was wrong.

**Key finding**: birefringence is encoded w PHASE velocity NOT GROUP velocity.
Wave-packet PHASE rotuje między L/R helicities, but ENVELOPE propagates uniformly.

### W2.3 — Polarization-averaged c_eff = 1 ✓ PASS

```
c_eff(phase) = (v_+ + v_-)/2 ~ 1 - (g n_par)^2/(8 k^2)  (quadratic suppression in PHASE)
c_eff(group leading O(eps)) = 1 EXACT
c_eff(group O(eps^2)) ~ 1 + (g n_par)^2 / (8 k^2)  (slight super-luminal envelope at O(ε²))
```

**No scalar c(X) at leading order** → consistent z Webb/Murphy α_em NULL 1e-7 precision.

### W2.4 — Birefringence Δv_φ ✓ PASS

```
Delta v_phi (leading) = v_phi_+ - v_phi_- = g n_par / k  EXACT
Delta chi(path) = (g/2) integral n_parallel ds  matches omega.1 W3.3
```

### W2.5 — Helicity-dependent optical cone ✓ PASS  *(language softened 2026-05-01)*

Modified Maxwell for helicity ± can be cast as a **helicity-dependent optical
cone / effective dispersion geometry**:
$$Q^{(\pm)}_{\mu\nu}(p)\, k^\mu k^\nu = 0,\qquad p_\mu \equiv g\,\partial_\mu \ln X$$
Special-case (static, `n_i ∥ k_i`) WKB inversion gives σ.1 dispersion
`ω_±² = k² ± g·k·n_∥`. Recovers Minkowski (`Q^(±) → η`) w n=0.

> **Important framing note**: this is **NOT** a full classical pseudo-Riemannian
> metric. It is a leading-order null-cone / effective dispersion geometry
> derived in WKB for a single special configuration. The full Carroll-Field-
> Jackiw-class dispersion is `(k²)² + p²k² − (p·k)² = 0` for arbitrary `p_μ`,
> which is more general than what σ.1 derives here. The σ.1 program LOCKS the
> WKB special case + framing; the full structure is left for future work.

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

- **Phase velocity birefringence v_φ_± LOCKED** (sympy exact): `v_φ,± = 1 ± g·n_∥/(2k) − (g·n_∥)²/(8k²) + …`
- **Group velocity v_g_± at leading O(ε) = 1 LOCKED** (no linear envelope splitting); **at O(ε²): `v_g,± = 1 + (g·n_∥)²/(8k²)` POSITIVE** (sympy + analytic, sign-corrected 2026-05-01)
- **No scalar c(X) at O(g·n_∥/k) leading** PROTECTED z Webb/Murphy NULL
- **Helicity-dependent optical cone / effective dispersion geometry** (NOT full classical metric — WKB special case)
- **Scale-invariance X→λX preserved** (inherits ω.1)
- **σ.1 free-wave standalone** (no substrate back-reaction at lowest order)
- **Phase 3 forward**: predictions + 4-channel convergence

## Cross-references

- [[program.md]]
- [[Phase1_results.md]]
- [[Phase2_setup.md]]
- [[../op-omega1-substrate-em-coupling/Phase2_results.md]] — modified Maxwell EOM source
