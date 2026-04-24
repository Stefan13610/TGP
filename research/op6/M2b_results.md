# OP-6 / M2-b — Numerical results: envelope stiffness of H₁

**Status:** CLOSED, 2026-04-23. Outcome: **FAIL** (K_eff(Φ) is not ∝ Φ).
**Inputs:** `M2a_analytical_sketch.md`, `M2b_numerical_plan.md`.
**Script:** `m2b_envelope_stiffness_1d.py`.
**Raw data:** `m2b_scan_results.txt`.

> **⚠ Corrective note (2026-04-23, see `M2c_H3_reality_check.md`).**
> §6–§7 below recommend pivoting to "H₃ (quadratic bond)". That
> recommendation is **algebraically wrong**: the bilinear bond
> `−J Σ (φ_iφ_j)² = −J Σ Φ_iΦ_j` gives α(φ) = 1, not α(φ) = 2.
> The correct pivot target is the GL functional H_GL, i.e.
> `Σ J(φ_iφ_j)² (φ_j − φ_i)²`. The numerical verdict in §3 and the
> closure of M1-B in §6 are unaffected; only the recommended
> replacement axiom changes. Read "H₃" in §6–§7 as "H_GL" and
> consult M2c for the full correction.

---

## 1. What was tested

The only mechanism left alive by M2-a (§2) for deriving the kinetic
functional F_kin = (K_eff/2) Φ (∇Φ)² — equivalently α = 2 in the φ
variable — from the canonical H₁ is **envelope coarse-graining at
fixed Φ = ⟨ŝ²⟩**.

The direct test: for H₁ on a 1D lattice with continuous ŝ, measure
the effective kinetic stiffness K_eff(⟨Φ⟩) from the Ornstein-Zernike
expansion of the Φ-Φ structure factor,

    1 / S(k_hat) ≈ 1/χ_Φ + K_eff · k_hat² ,   k_hat² = 2(1 − cos k) .

Then scan control parameters (β, m₀², λ₀) to vary ⟨Φ⟩ and fit
K_eff(⟨Φ⟩) = C · ⟨Φ⟩^p. The TGP target is p = +1 (⇔ K_eff ∝ Φ ⇔
α = 2 in φ). The Gaussian (λ₀ = 0) sub-limit predicts p = −1
analytically (§3 below), serving as a validation of the method.

## 2. Method validation (Part B, Gaussian sub-limit)

For H₁ with λ₀ = 0 the partition function is exactly Gaussian and
one can compute K_eff analytically:

    cosh(ν) = m₀²/(2J),
    ⟨Φ⟩_exact = 1 / (2β J sinh(ν)),
    K_eff_exact = (β²J²/2) tanh(ν).

Hence K_eff · ⟨Φ⟩ = (β J / 4) · sech(ν) = β / (2 m₀²), which varies
slowly with m₀². Across the scan range (J=1, m₀² ∈ {2.02, 2.05,
2.10, 2.20}, β ∈ {0.5, 1, 2, 4}), the analytical log-log slope is
p_analytical = −1.744.

The MC measurement gives p_MC = −1.706 ± 0.676 (n = 10 valid
points), consistent with the analytical value within statistical
noise. Individual K_eff measurements are very noisy (typical
jackknife error ~comparable to the value), because ξ_Φ ≤ 3
lattice units and the OZ fit has only 3–4 modes of curvature to
work with, but the *exponent* is determined by the log-range
spanned by ⟨Φ⟩ (~1 decade) and is much less sensitive to
per-point noise. **Method validated.**

## 3. Result (Part C, Full H₁ with λ₀ = 1)

Scan: β ∈ {0.3, 0.5, 0.7, 1, 1.5, 2, 3}, m₀² ∈ {−1, 0, 1, 2, 3},
λ₀ = 1, J = 1, N = 1024, n_therm = 3000, n_measure = 8000 sweeps,
jackknife n_jk = 8.

Aggregate log-log fit:

    K_eff(⟨Φ⟩) = 29.7 · ⟨Φ⟩^(−0.28 ± 0.43)   (n = 20 valid points)

- **p = −0.28 ± 0.43.**
- **z-score against TGP target p = +1: (−0.28 − 1.00)/0.43 = −3.0 σ.**
- **TGP target p = +1 is rejected at 3σ.**

Qualitatively: K_eff is approximately Φ-independent (p ≈ 0),
**not** linearly growing with Φ (p = +1), in the 1D envelope
coarse-graining of H₁.

## 4. Interpretation

The numerical result in 1D is consistent with the analytical
narrowing of M2-a and inconsistent with the TGP target. Specifically:

- Mechanisms (i)–(iv) were analytically ruled out in M2-a (LPA,
  LPA' with Z_k(ρ), derivative expansion, φ²(∇φ)² operator).
- Mechanism (v) — envelope coarse-graining at fixed Φ — now fails
  in the direct 1D MC test as well. K_eff is neither ∝ Φ (TGP
  target) nor even positively correlated with ⟨Φ⟩ in any regime
  probed.
- A sharp analytical indicator (M2-a §2.v.b) also disfavoured the
  WF result: Δ_ε ≠ 2 Δ_σ in Ising universality.

The 1D test is preliminary (not at the 3D WF fixed point), but
combined with the analytical arguments it is strong evidence
that the M1-B programme — derive α = 2 from H₁ by RG — does
not close.

## 5. Caveats

- 1D has no Wilson–Fisher fixed point. The 1D measurement is a
  **surrogate** for the full WF answer, not a substitute. It could
  in principle happen that K_eff(Φ) changes sign of slope when
  one moves from d = 1 (disordered massive theory) to d = 3 at WF.
  This is the escape hatch for M1-B.
- The OZ fit on small ξ_Φ is noisy; individual K_eff values have
  large jackknife errors. The conclusion rests on the *exponent*
  of the log-log fit, which is more robust than per-point K_eff.
- Our Gaussian analytical prediction for the slope in this
  parameter range is p ≈ −1.74, not exactly −1, because we are
  not in the strict critical (ν → 0) limit. The TGP target p = +1
  is far enough from any value produced by our MC to reject it
  regardless of these finite-size/range caveats.

## 6. Consequence for OP-6

1. **M1-B is numerically disfavoured.** The hypothesis that H₁
   (bilinear bond, core-paper axiom) produces K(Φ) ∝ Φ under
   envelope coarse-graining is not supported by either the
   analytical sketch (M2-a) or the 1D MC measurement (this).
2. The **next action** (*corrected 2026-04-23 per M2c*) is not
   "M1-A = pivot to bilinear H₃", because the bilinear bond
   `−J Σ (φ_iφ_j)² = −J Σ Φ_iΦ_j` gives α(φ) = 1, not α(φ) = 2.
   The correct pivot target is **M1-A′ = axiom = H_GL**:
   adopt as a substrate-level axiom the GL functional
   `F = Σ_⟨ij⟩ J (φ_iφ_j)² (φ_j − φ_i)² + on-site terms`,
   equivalently the continuum form
   `F_kin^{sub}[φ] = (K_geo/2) ∫ φ⁴ (∇φ)² d^dx`.
   This is exactly what `prop:substrate-action` in dodatekB
   actually uses. The cost is that the axiom is an
   *effective-theory* GL functional, not a microscopic spin
   Hamiltonian. See `M2c_H3_reality_check.md` §3 for the full
   comparison of pivot options.
3. M3 (3D FRG/MC confirmation) would be the escape-hatch check
   to see if the 1D failure is 1D-specific. It is not strictly
   required before pivoting to M1-A′, but it would harden the
   case for v2.

## 7. Recommended workshop edits (non-blocking)

*(Corrected 2026-04-23 per M2c: "H₃" was the wrong pivot target;
the correct one is H_GL.)*

After v2 of the core paper is actually planned, these are the
edits to make:

1. Core paper tex: axiom `eq:H-Gamma` becomes **H_GL**, phrased
   either as a microscopic GL bond
   `Σ_⟨ij⟩ J (φ_iφ_j)² (φ_j − φ_i)² + on-site terms` or directly
   as the continuum functional
   `F_kin^{sub}[φ] = (K_geo/2) ∫ φ⁴ (∇φ)² d^dx`. Theorem
   `thm:alpha2` is unchanged (it is a classification inside
   (C1)–(C3)); but the substrate Hamiltonian now supplies
   (C1)–(C3) directly, instead of via the gap that v1 left open.
2. dodatekB `prop:substrate-action` is unchanged in statement
   and proof — it already starts from the H_GL form.
3. `M1_H_Gamma_canonical.md` §5: recommend M1-A′ (pivot to H_GL)
   instead of M1-B. Distinguish H_GL from the bilinear "H₃"
   inventoried in the original M1 §1.
4. `TGP_v1/README.md` etc.: document the pivot and the H_GL/H₃
   distinction.

These are all v2 content. For **v1 (current release)**, the
correct stance remains: OP-6 is an open problem (as already
stated in KNOWN_ISSUES.md and the published paper).

## 8. Exit status

- M2-b: **closed**, verdict **FAIL** for M1-B.
- M2-a: closed (no change).
- M1: closed with the refinement that M1-B is disfavoured;
  M1-A is the recommended path for any future v2 of the core
  paper.
- M3: not required to advance; can be added later to harden the
  case in 3D.
- M4–M5: contingent on M1-A; content carries over if the axiom
  is switched to H₃ because the effective theory for Φ is then
  derived from H₃ directly as in M1 §4.2.
- v2 core: now becomes an M1-A-style update (change axiom),
  not an M1-B-style derivation. This is conceptually clean and
  honest about scope.
