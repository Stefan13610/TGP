# M2 — Effective potential U_eff(φ) from H_Γ (analytical sketch)

**Status:** Plan only. Not yet executed.
**Target:** derive `U_eff(φ) = (β/3) φ³ − (γ/4) φ⁴ + R(φ)` with
`R(φ)` composed only of RG-irrelevant operators at the 3D Ising
Wilson–Fisher fixed point, starting from the v2 substrate
Hamiltonian `H_Γ` on a `d=3, z=6` hypercubic lattice.

Closes **OP-1** (form uniqueness) and **Route 2 of OP-2**
(all-scale protection via universality). Supplies the structural
input for **M4** (OP-4, value of γ).

---

## 1. Starting point

    H_Γ[ŝ] = Σ_i [π_i² / 2μ + (m₀²/2) ŝ_i² + (λ₀/4) ŝ_i⁴]
             + J Σ_⟨ij⟩ A_ij · ŝ_i² ŝ_j² (ŝ_j² − ŝ_i²)²        (v2)

with the classical partition function

    Z = ∫ Dŝ e^{−β H_Γ[ŝ]},     β = 1/T.

Order parameter: `Φ_i = ŝ_i²` (Z₂-even), `φ_i := Φ_i / Φ₀`.

## 2. Outline of the derivation

We want the 1PI generating functional `Γ[φ]` (effective action).
Two independent routes. Run both; they must agree.

### Route A — Hubbard–Stratonovich on the composite Φ

**Step A.1.** Linearise the local quartic using an auxiliary field
σ_i conjugate to Φ_i = ŝ_i²:

    exp[−(λ₀/4) ŝ_i⁴]
      = ∫ Dσ_i · exp[−(1/λ₀) σ_i² − σ_i ŝ_i²]         (up to const).

This trades the quartic term for a Gaussian integral over σ
coupled linearly to Φ.

**Step A.2.** Do the same for the bond term. The v2 bond
`ŝ_i² ŝ_j² (ŝ_j² − ŝ_i²)²` contains factors of Φ_i and Φ_j up to
the 4th power; write

    ŝ_i² ŝ_j² (ŝ_j² − ŝ_i²)² = Φ_i Φ_j (Φ_j − Φ_i)²

and use a bilinear H-S decoupling on the factor `(Φ_j − Φ_i)²`
with conjugate field ξ_{⟨ij⟩}:

    exp[−β J A_ij Φ_i Φ_j (Φ_j − Φ_i)²]
      = ∫ Dξ · exp[−ξ_{ij}² / (β J A_ij Φ_i Φ_j) + i ξ_{ij} (Φ_j − Φ_i)].

Notice the dependence of the ξ-field Gaussian weight on `Φ_i Φ_j`:
this is the source of the `φ⁴` structure in the effective kinetic
term (`prop:substrate-action`) and must reproduce it.

**Step A.3.** Gaussian integral over ŝ. After A.1 + A.2 the
remaining ŝ-dependence is bilinear (kinetic + `m₀² ŝ²` +
`σ ŝ²`), plus a Z₂-odd source `ξ ŝ²`:

    Z = ∫ Dσ Dξ Dŝ · exp[−(1/λ₀) σ² − ξ-weights + S_quad[ŝ; σ, ξ]].

Performing the ŝ Gaussian integral gives a `(1/2) Tr ln M[σ, ξ]`
correction to the σ, ξ action, where `M` is the ŝ-field Hessian.

**Step A.4.** Legendre transform. Define the effective action
`Γ[Φ, ⟨ξ⟩]` by Legendre-transforming over σ and then projecting
onto the subspace where the source for σ is just `−Φ/Φ₀` (the
order parameter). Expanding `Γ` around the broken-Z₂ saddle
`Φ = Φ₀, ξ = 0`:

    Γ[Φ] = Γ[Φ₀] + ∫ d^dx [ (K_geo/2) φ⁴ (∇φ)² + U_eff(φ) + higher ].

The kinetic part must match `prop:substrate-action` (consistency
check: non-trivial validation that the derivation is correct).

**Step A.5.** Local piece of `Γ`. Collect all local terms in
`φ = Φ/Φ₀`:

    U_eff(φ) = c₀ + c₁ φ + c₂ φ² + c₃ φ³ + c₄ φ⁴ + c_k φ^k (k ≥ 5).

Z₂ on ŝ makes `Φ = ŝ²` non-negative and `φ ≥ 0`; this does *not*
forbid odd powers of φ (since φ is not Z₂-odd — it is already
Z₂-even). Instead, constraints on `c_k` come from:

- **c₀:** absorbed in the vacuum energy → redefine `Γ[Φ₀] = 0`.
- **c₁:** must vanish at the saddle `φ = 1` (choice of Φ₀).
  Combined with the normalisation `U(1) = −γ/12` (implied by
  β = γ vacuum), this fixes two relations among `c_k`.
- **c₂:** gives the mass of fluctuations of φ around vacuum.
  Physical interpretation: `m_sp² = γ` (dodatek B L533), which is
  a *prediction* of the derivation, not an input.
- **c₃, c₄:** the target coefficients β/3, −γ/4.
- **c_k (k ≥ 5):** RG irrelevance to check against 3D Ising WF
  scaling dimensions (see Route B and §4).

**Expected outcome of Route A:**

    c₃ = β / 3 = (1/3) (λ₀ / a_sub²) C_β,
    c₄ = −γ / 4 = −(1/4) (λ₀² / (v² a_sub²)) C_γ,
    c_k (k ≥ 5) → 0 in the IR (RG irrelevance).

The constants `C_β, C_γ` will come out as explicit loop integrals
over the ŝ Gaussian propagator; plug them into M4 to close OP-4.

### Route B — Wilsonian RG (momentum shell)

**Step B.1.** Decompose `ŝ = ŝ_slow + ŝ_fast` in momentum shells.
Integrate out `ŝ_fast` perturbatively in `λ₀` to one loop.

**Step B.2.** Read off the resulting effective action for
`ŝ_slow`, re-express in the `Φ_slow = ŝ_slow²` variable.

**Step B.3.** Identify the induced cubic and quartic vertices in
Φ_slow. Cross-check against Route A.

**Step B.4.** Match to MK-RG (§M3) by lattice-regularising the
momentum shell at scales `a_sub, 2 a_sub, 4 a_sub, …`

Route B is more technical but less prone to H-S ambiguities.
Routes A and B *must* yield the same `β, γ`; discrepancy means an
error in one of them.

## 3. Operator enumeration for the local piece

At the 3D Ising WF fixed point, primary operators and their
scaling dimensions (conformal bootstrap):

| Operator | Δ (3D Ising) | RG class |
|---|---|---|
| ε = φ² | 1.41262… | relevant |
| ε' = φ⁴ | 3.830 | IRRELEVANT |
| ε'' = φ⁶ | ≈ 6.9 | irrelevant |
| T_{μν} (stress tensor) | 3 (exact) | marginal |
| φ | 0.51814… | relevant |
| φ³ (descendant of φ) | ≥ 3.52 | irrelevant at subleading |

**Classification of c_k:**

- `c_3 φ³`: descendant of φ at Δ = 1.5, relevant.
- `c_4 φ⁴`: ε' operator at Δ = 3.83, IRRELEVANT at IR.

Wait — this collides with Route 2 of `thm:beta-eq-gamma-triple`
which asserts β = γ is protected by Z₂ symmetry. **Resolution:**
in the 3D Ising universality class the *relevant* deformation is
`φ² = ε`, and the combination `(β/3) φ³ − (γ/4) φ⁴` is the
*lowest* operator cluster that carries the non-trivial
interaction content. The cubic and quartic are not independent
operators at the WF fixed point — they are tied by the single
self-coupling `u*` of the Wilson–Fisher flow. This is the
content of "β = γ" Route 2 and M2 must reproduce it *as an
identity*, not as a dynamical fine-tuning.

## 4. Irrelevance of higher operators

Sextic φ⁶ at WF has Δ ≈ 6.9 ≫ 3 = d ⇒ irrelevant.
Octic φ⁸ even more so. Therefore at the IR fixed point only the
cubic + quartic structure survives, matching (♦).

**This is the structural content of OP-1 closure:** uniqueness of
U(φ) up to cubic+quartic is dictated by the 3D Ising WF fixed
point, modulo irrelevant corrections of calculable size.

## 5. Consistency checks

Successful M2 must reproduce:

1. **Kinetic match:** local part of `Γ[φ]` yields
   `F_kin^geo[φ] = ∫ (K_geo/2) φ⁴ (∇φ)²` (prop:substrate-action).
2. **β = γ at vacuum:** `c_3 / c_4 = −β/3 / (−γ/4) · (1/something)`
   must collapse to `β = γ` after evaluating at the WF fixed point.
3. **Dynamical coefficient:** `W(1) = 7β/3 − 2γ = γ/3` when β=γ
   (sek08 L1483). This requires the "time-like" part of `Γ[φ]` to
   carry the same β, γ with the `7/3, 2` ratios.
4. **Substrate → TGP map:**
   `β = (λ₀ / a_sub²) C_β` and `γ = (λ₀² / (v² a_sub²)) C_γ`
   must emerge as loop-integral results, not be put in by hand.
5. **Second-sound mass:** `m_sp² = γ` (dodatek B L533) as a
   derived output.

If all five check out, M2 closes.

## 6. Tools

- sympy for symbolic expansions.
- Mathematica (optional) for the `Tr ln M` expansion.
- For numerical cross-checks, reuse the 1D MC code in
  `research/op6/m2b_envelope_stiffness_1d.py` (sentinels for the
  operator-mixing matrix).

## 7. Known risks

- **Route-A ambiguity:** multiple equivalent H-S decouplings.
  Mitigation: choose the *collective-field* decoupling
  (`σ ↔ Φ`) which has a direct physical interpretation.
- **Saddle-point expansion blowing up:** the GL bond's `(Φ_j − Φ_i)²`
  factor vanishes at the saddle, so naive perturbation may have
  flat directions. Mitigation: separate the "bond gradient"
  channel from the "bond amplitude" channel.
- **Mismatch with Route B:** if present, diagnose by checking the
  anomalous dimension of `ŝ²` via the two routes against the
  3D Ising bootstrap value `η_ε ≈ 0.03627`.

## 8. Deliverables

- **sympy worksheet** doing the A.1 → A.5 computation on the
  lattice, truncated at 1 loop and `O(φ⁴)` in the expansion.
- **LaTeX writeup** suitable for inclusion as appendix in
  dodatek B-II (or a new dodatek O).
- **Table of induced operators** with their WF scaling dimensions
  and RG class.

## 9. Next step

Execute Step A.1 on paper first (a couple of lines), then code
A.1-A.3 in sympy. If A.3 (Gaussian integral over ŝ) goes through,
A.4-A.5 are straightforward.
