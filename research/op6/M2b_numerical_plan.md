# OP-6 / M2-b — Numerical plan for envelope stiffness of H₁

**Status:** working document, 2026-04-23.
**Predecessor:** `M2a_analytical_sketch.md`.
**Input:** canonical H₁ (M1), the analytical narrowing of M2-a.
**Question:** in envelope coarse-graining of H₁, is K_eff(Φ) ∝ Φ
(equivalently α=2 in φ), or not?

---

## 1. Why we cannot use an off-the-shelf FRG

The analytical sketch (M2-a, §2.v) narrows the remaining path to
**envelope coarse-graining**: integrate out the fast ŝ fluctuations
*at fixed profile of the Z₂-even density* Φ(x) = ⟨ŝ(x)²⟩.

This is not the same computation as standard LPA' for the WF
fixed point in 3D Ising:

- Standard LPA': Γ_k[σ] with σ = ⟨ŝ⟩ (Z₂-odd "magnetisation").
  Ansatz Γ_k = ∫ [(Z_k(σ²)/2)(∇σ)² + U_k(σ²)].
- Envelope: Γ_k[Φ] with Φ = ⟨ŝ²⟩ (Z₂-even "density").
  Ansatz Γ_k = ∫ [(K_k(Φ)/2)(∇Φ)² + V_k(Φ)].

In the symmetric phase σ = 0 but Φ ≠ 0, so these are distinct
effective actions. There is no direct conversion between Z_k(σ²)
and K_k(Φ). A proper FRG for Γ_k[Φ] requires either

(a) a two-field flow Γ_k[σ, Φ] with a Hubbard–Stratonovich
    constraint Φ = ŝ², then Legendre-transforming,
(b) a composite-operator RG treating Φ = ŝ² as an independent
    field with its own propagator,
(c) direct computation of the Φ-Φ 4-point function from the
    underlying ŝ theory.

(a)–(c) are all genuine research projects. None is a one-script
implementation.

**Conclusion.** The tractable M2-b deliverable is not a custom
FRG — it is a **direct Monte-Carlo measurement** of K_eff(⟨Φ⟩)
from H₁ configurations, via the Φ-Φ correlator.

## 2. What we actually measure

Assume the effective action has the Ornstein–Zernike form near
the symmetric point ⟨Φ⟩ = Φ₀:

    F_eff[Φ] = ∫ d^d x [ (K_eff(Φ)/2) (∇Φ)² + V_eff(Φ) ]

Expanding around Φ₀, the Gaussian fluctuations of δΦ = Φ − Φ₀
satisfy (in momentum space)

    ⟨δΦ(−q) δΦ(q)⟩ = 1 / [ V_eff''(Φ₀) + K_eff(Φ₀) q² ]

In OZ-normalised form

    G_Φ(q) = χ_Φ / (1 + ξ_Φ² q²)

where χ_Φ = 1/V_eff'', ξ_Φ² = K_eff/V_eff'' = K_eff · χ_Φ. This
gives the key identity

    **K_eff(⟨Φ⟩) = ξ_Φ² / χ_Φ**                              (*)

with

    χ_Φ = Σ_x C_Φ(x),
    M_2 = Σ_x x² C_Φ(x),
    ξ_Φ² = M_2 / (2 d χ_Φ),

where C_Φ(x) = ⟨Φ(0)Φ(x)⟩ − ⟨Φ⟩² is the connected two-point
function of Φ and d is the spatial dimension.

Formula (*) converts MC measurements of the Φ-Φ correlator into
an estimate of K_eff.

To extract the **Φ-dependence** of K_eff, we scan control
parameters (β, m₀², λ₀) to vary ⟨Φ⟩, giving a set of pairs
(⟨Φ⟩, K_eff(⟨Φ⟩)). A fit to

    K_eff(Φ) = C · Φ^p

then reports the exponent p. The four regimes of interest:

| p | α(φ) | Interpretation |
|---|------|----------------|
| −1 | 0 | Gaussian (§M2-a.1.3): K(ŝ)=const, K(Φ)=const/Φ |
|  0 | 1/2 | K constant in Φ, "flat stiffness" |
| +1 | 2 | **TGP target: K(Φ) ∝ Φ ⇔ K(φ) ∝ φ⁴** |
|  · |  · | anything else: non-canonical |

The M2-a §1.3 algebra predicts p = −1 (Gaussian sub-limit).
Mechanism (v) is viable only if we find **p ≈ +1** in the
interacting regime.

## 3. Implementation: `m2b_envelope_stiffness_1d.py`

A 1D lattice simulation is chosen as the first test because:

1. No phase transition in 1D, so the correlator decays
   exponentially with a finite, measurable ξ_Φ for all (β, m²,
   λ). No critical slowing-down, no finite-size effects on
   ξ that would dominate.
2. The question "is K_eff(Φ) ∝ Φ as an algebraic functional?"
   does not require being at a fixed point. If the answer is NO
   in 1D, it does not directly prove NO in 3D, but it is
   strong preliminary evidence that envelope coarse-graining
   does **not** produce such a kinetic dependence.
3. 1D MC is fast and allows comfortable scans over (β, m², λ).

### 3.1 Hamiltonian

    H₁_1D = Σ_i [ (m₀²/2) ŝ_i² + (λ₀/4) ŝ_i⁴ ] − J Σ_⟨ij⟩ ŝ_i ŝ_j

Continuous ŝ_i ∈ ℝ; periodic BC; Z₂ symmetry ŝ ↦ −ŝ.

### 3.2 MC scheme

- Single-site Metropolis with Gaussian proposal ŝ' = ŝ + η,
  η ~ N(0, σ_prop²).
- Sweep = N local updates. Thermalise for n_therm sweeps.
- Measure every sweep for n_sweep sweeps.
- Observable set: Φ_i = ŝ_i²; C_Φ(x) via FFT; χ_Φ and ξ_Φ via (*).

### 3.3 Parameter scan

Fix J = 1, λ₀ = 1; vary (β, m₀²) to produce a range of
⟨Φ⟩ values spanning roughly one decade. Typical ranges:

    β ∈ {0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0}
    m₀² ∈ {0.5, 1.0, 2.0}
    λ₀ = 1.0

For each (β, m₀²) combination, record (⟨Φ⟩, K_eff). Combine all
into a single scatter plot and fit power-law K_eff = C · ⟨Φ⟩^p.

### 3.4 Deliverables

1. Table of (β, m₀², ⟨Φ⟩, χ_Φ, ξ_Φ², K_eff).
2. Log-log plot of K_eff vs ⟨Φ⟩ (PNG).
3. Best-fit exponent p with uncertainty.
4. Verdict:
   - p ∈ (0.7, 1.3) → M2-b passes, M1-B vindicated, proceed to
     M3 (3D FRG/MC).
   - p ∈ (−0.3, 0.3) → K_eff is Φ-independent; α(Φ) = 1/2,
     α(φ) = 1. Partial: need 3D confirmation before pivot.
   - p ∈ (−1.3, −0.7) → Gaussian-like. M1-B fails; pivot to
     M1-A (change axiom to H₃).
   - Anything else → weird, investigate.

### 3.5 Sanity checks

- **Free limit (J = 0):** ŝ_i are independent, Φ_i =
  ŝ_i² are independent, C_Φ(x) = 0 for x ≠ 0. Should give
  ξ_Φ = 0 → K_eff = 0. A check that the code does not
  spuriously report K_eff ≠ 0 when it shouldn't.
- **Pure Gaussian (λ₀ = 0, m₀² > 0, J small):** C_Φ(x) =
  2 G_ŝ(x)² (Wick contraction); K_eff computable analytically.

## 4. Relation to M2-a mechanism (v)

M2-a §2.v concluded that mechanism (v) is the only remaining
candidate analytically but that its WF behaviour cannot be
derived analytically without FRG. This 1D MC is a first-pass
surrogate: if the answer in 1D is clearly inconsistent with
K_eff ∝ Φ (as M2-a expects via the Δ_ε ≠ 2Δ_σ mismatch), we
have strong evidence against mechanism (v). A full 3D test
with continuous φ⁴ on the lattice remains as M3.

## 5. If this test fails

If K_eff(Φ) is clearly NOT ∝ Φ in the MC scan, we have
numerical evidence that H₁ does not produce α=2 via envelope
coarse-graining. Combined with M2-a's elimination of
mechanisms (i)–(iv), this closes M1-B: **α=2 is not a
consequence of the published core-paper axiom**.

The correct follow-up is M1-A: change the substrate
Hamiltonian axiom in a (future) v2 of the core paper from H₁
to a form with a quartic bond (H₃ / H₄). The change is
conceptually clean — it makes the axiom match what the
derivation actually needs — and it is honest about the
scope of the claim.
