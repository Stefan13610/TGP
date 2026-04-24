# M3 — Numerical computation of C_β / C_γ at the WF fixed point

**Status:** Plan only. Not yet executed.
**Target:** compute the ratio `C_β / C_γ` at the 3D Ising
Wilson–Fisher fixed point of the Migdal–Kadanoff RG, starting from
the v2 substrate Hamiltonian on a `d=3, z=6` hypercubic lattice.

Closes **Route 3 of OP-2** by testing the prediction
`C_β / C_γ ≈ 1.74` needed for `β = γ` at criticality (dodatek B
eq:B-bg-ratio).

---

## 1. Starting point

From dodatek B §app:B-mapa:

    β = (λ₀ / a_sub²) · C_β,
    γ = (λ₀² / (v² a_sub²)) · C_γ.

At the WF fixed point `r* = −2.251, u* = 3.917` (dodatek B
eq:B-alpha-eff):

    β / γ |_* = (|r*| / u*) · (C_β / C_γ) = 0.575 · (C_β / C_γ).

Vacuum condition `β = γ` ⇒ `C_β / C_γ ≈ 1.74`.

## 2. What needs extending in the existing MK-RG code

The MK-RG pipeline used in dodatek B (reference:
`tooling/scripts/substrate/renormalization_substrate.py`) currently
tracks:

- `r_ℓ` (quadratic coupling),
- `u_ℓ` (quartic coupling ŝ⁴).

**M3 extends to track also:**

- `B_ℓ` — coefficient of the Φ³ = ŝ⁶ operator (cubic in Φ),
- `Γ_ℓ` — coefficient of the Φ⁴ = ŝ⁸ operator (quartic in Φ),
- `G_ℓ` — coefficient of the GL gradient bond `Φ_i Φ_j (ΔΦ)²` to
  match M2's induced kinetic operator.

Then:

    C_β = B_* / v*²,          C_γ = Γ_* / v*²,
    where v*² = |r*| / u* = 0.575.

## 3. Algorithm

1. Initialise `(r₀, u₀, J₀, B₀=0, Γ₀=0, G_0 = J₀ A)`.
2. Apply one Migdal–Kadanoff step (bond-moving + decimation) on the
   `d=3, z=6` hypercubic lattice. Generate the induced
   operator mixing matrix.
3. Read off the flow of `(r, u, B, Γ, G)` under the RG step.
4. Iterate until fixed point `(r*, u*, B*, Γ*, G*)`.
5. Compute `C_β / C_γ = (B*/v*²) / (Γ*/v*²) = B* / Γ*`.

## 4. Operator basis

Effective action truncated at 4th order in Φ and 2nd order in
gradients:

    F_eff[Φ] = ∫ d³x [
          (r/2) Φ
        + (u/4) Φ²
        + (B/3) Φ³
        + (Γ/4) Φ⁴
        + (K_geo/2) Φ_i Φ_j (∇Φ)²    (GL gradient)
        + higher                      (checked irrelevant)
    ]

where the last gradient term couples to the scale `G`.
Operator-mixing under MK step: this 5-operator truncation should
close under one RG iteration up to irrelevant `Φ⁵` and
`Φ^n (∇Φ)²` (n ≥ 3) terms.

## 5. Validation against op6 M3-a

The M3-a code (`research/op6/m3a_block_rg_1d.py`) validates that
block-RG gives the correct trivial sub-limits (free field,
Gaussian). M3 here will run the same validation protocol in 3D
MK-RG for the Gaussian and Ising sub-limits:

- Gaussian (u=0, B=0, Γ=0): r flows trivially, C_β = 0, C_γ = 0.
  Sanity check.
- Pure Ising (B=0, Γ=0, only r, u): must recover the published
  MK-RG fixed point `r* = −2.251, u* = 3.917`.

## 6. Expected sources of the cubic operator

At the free-theory starting point `B₀ = 0`, so a cubic operator
must be *generated* by the RG flow. Candidate generation
mechanisms:

1. **Operator mixing from ŝ⁶ contractions.** Two `u · ŝ⁴` vertices
   contracted on one internal ŝ line generate a `ŝ⁶ = Φ³` vertex
   at the integrated-out scale. This is the leading 1-loop source
   of `B_ℓ`.
2. **GL bond induced cubic.** Expanding `Φ_i Φ_j (ΔΦ)²` around the
   saddle `Φ = Φ₀` gives a quadratic + cubic + quartic local
   series in the fluctuation `δΦ`; the cubic piece contributes to
   `B`.
3. **Wave-function mixing.** If the order parameter gets
   multiplicative renormalisation `Φ → Z_Φ · Φ`, then the effective
   `(r/2) Z_Φ Φ + ...` generates cubic after shift to
   `φ = Φ / Φ₀`.

All three should be coded and their relative magnitudes compared
with Route A of M2.

## 7. Deliverables

- **Extended Python script** `mk_rg_bgamma.py` analogous to the
  existing dodatek B code, with operator-mixing matrix.
- **Scan table** of `(r_ℓ, u_ℓ, B_ℓ, Γ_ℓ)` flow.
- **Fixed-point value** `C_β / C_γ = 1.74 ± δ` with error bar `δ`
  from truncation uncertainty.
- **Results markdown** `M3_results.md` documenting outcome and
  comparison with M2 Route A.

## 8. Cross-checks against the conformal bootstrap

At the 3D Ising WF fixed point, scaling dimensions are known to
high precision:

- Δ_ε = 1.41262528(29) — ε = φ² operator,
- Δ_{ε'} = 3.8303(18) — ε' = φ⁴ operator.

The eigenvalues of the MK-RG flow matrix linearised around
`(r*, u*, B*, Γ*)` should give, after re-identification:

    y_ε = d − Δ_ε = 1.587,
    y_{ε'} = d − Δ_{ε'} = −0.830  (irrelevant).

Getting the MK-RG eigenvalues to agree with bootstrap within the
known MK systematic (typically 5–10%) is a sanity check.

## 9. Timeline

- 1 session: extend the MK-RG code with cubic + quartic Φ
  operators.
- 1 session: scan + fit fixed-point ratio.
- 1 session: write up results + compare with M2.

## 10. Success / failure criteria

| Outcome | Interpretation |
|---|---|
| `C_β / C_γ` converges to `1.74 ± 0.09` (within MK systematic) | **OP-2 Route 3 closed.** Combined with M2 Route 2, β = γ at all scales is established. |
| `C_β / C_γ` converges to a value incompatible with 1.74 (≥ 3σ) | β = γ is NOT protected at criticality; `thm:beta-eq-gamma-triple` Route 3 fails. Opens a new OP: why does Route 1 (exact at vacuum) not extend? |
| `C_β / C_γ` does not converge under truncation | M3 needs higher-order truncation. Extend to Φ⁵ and re-run. |

## 11. Parallelism with op6 M3-a

This is the analogue at the *potential* level of what M3-a was at
the *bond* level. Same methodology (block-RG), same validation
pattern (Gaussian sub-limit, run on physical limit, statistical
error from repeated runs). Difference: 3D hypercubic lattice
instead of 1D chain, 5-operator truncation instead of 2-operator.

## 12. Next step

Code stub of `mk_rg_bgamma.py` in this folder. First target:
Gaussian sub-limit validation. Then extend to full v2 GL bond.
