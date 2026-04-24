# OP-1 / OP-2 / OP-4 — derivation of the self-interference potential U(φ)

This folder is the clean workspace for a joint attack on the three
open problems that together form the **potential cluster** of TGP
core:

- **OP-1** — uniqueness of the form
  `U(φ) = (β/3) φ³ − (γ/4) φ⁴`.
- **OP-2** — protection of the vacuum relation `β = γ` at all scales
  (not only at `φ_vac = 1`).
- **OP-4** — derivation of the quartic coefficient `γ` (equivalent to
  the TGP coupling `g₀ᵉ = 0.86941`) from `H_Γ` alone.

After the v2 axiom pivot (2026-04-24, OP-6 axiomatic closure) the
derivation chain

```
  H_Γ  →  F_kin^geo[φ]  =  ∫ (K_geo/2) φ⁴ (∇φ)² dx
```

is **rigorous** at the kinetic level (`prop:substrate-action`,
`thm:alpha2`). The **potential** side

```
  H_Γ  →  U_eff(φ)  =  (β/3) φ³ − (γ/4) φ⁴ + (irrelevant)
```

is currently **postulated** as the minimal cubic+quartic form
compatible with Z₂ symmetry of the composite `Φ = ŝ²`. The aim of
this workspace is to promote that postulate to a theorem, and hence
to close OP-1, OP-2, OP-4 jointly.

## Programme

- **M1** — inventory and canonical form of `U(φ)` across the TGP
  corpus (core paper, companion papers, dodatek B, dodatek N,
  dodatek Q). Pin down which parameters are definitional, which are
  external inputs, and which OPs are genuinely coupled.
  *Status: not started.* See `M1_potential_inventory.md`.
- **M2** — analytical derivation of the effective potential
  `U_eff(φ)` from `H_Γ` (v2 GL bond) by Hubbard-Stratonovich
  transformation on the composite field `Φ = ŝ²`, followed by
  saddle-point expansion. Expected output:
  - cubic term `β φ³` from the local `λ₀ ŝ⁴` vertex projected onto
    the Z₂-even channel,
  - quartic term `γ φ⁴` from two-loop contribution + GL bond
    contraction,
  - higher terms irrelevant at the 3D Ising WF fixed point
    (Δ_ε = 1.413, sextic dimension ≥ 3).
  *Status: not started.* See `M2_effective_potential_sketch.md`.
- **M3** — numerical closure of `C_β / C_γ ≈ 1.74` via extended
  MK-RG on the v2 GL Hamiltonian. Leverages existing MK-RG
  pipeline from dodatek B §app:B-MK (where `β/γ |_* = 0.575 · C_β/C_γ`
  is already derived analytically). Delivers:
  - cubic-operator flow equation for `β(ℓ)`,
  - quartic-operator flow equation for `γ(ℓ)`,
  - fixed-point ratio `C_β / C_γ` to ≤ 5% accuracy,
  - test of the prediction `C_β / C_γ ≈ 1.74`.
  *Status: not started.* See `M3_CbCg_numerical_plan.md`.
- **M4** — derivation of `γ` as a function of `(J, λ₀, m₀², d, z)`
  from the M2 effective action, and confrontation with the
  empirical TGP calibration `g₀ᵉ = 0.86941`.
  *Status: blocked on M2 + M3.*

## Relation to existing work

- `axioms/substrat/dodatekB_substrat.tex` §app:B-MK already contains
  `β/γ = 0.575 · C_β/C_γ` at the Wilson-Fisher fixed point of MK-RG,
  and Theorem `thm:beta-eq-gamma-triple` with three partial routes
  to `β = γ`. M3 is the numerical completion of Route 3 of that
  theorem; M2 is the analytical deepening of Route 2.
- `partial_proofs/most_gamma_phi/dodatekQ2_most_gamma_phi_lematy.tex`
  collects coarse-graining lemmas that are candidate input to M2.
- `core/formalizm/dodatekQ_coarse_graining_formal.tex` holds the
  formal skeleton that M2 will instantiate.
- `core/formalizm/dodatekN_erg_renormalizacja.tex` hosts the
  ergodic + RG machinery we will call out.
- `research/op6/` is the predecessor — OP-6 closed at the axiom
  level (v2 pivot). This workspace takes the next step down the
  derivation chain, from kinetic to potential.

## Non-goals

- Changing any core-paper axiom. `H_Γ` v2 is frozen; this is a
  derivation effort *inside* the v2 axiomatic.
- Re-deriving the kinetic side `F_kin^geo`. That is
  `prop:substrate-action` in dodatek B; this workspace cites it
  as a theorem.
- Touching the companion papers (TGP-QM, TGP-Leptons, TGP-SC).
  Those already work at the coarse-grained `Φ` level and are
  transparent to any internal refinement of `U(φ)`.
- Re-opening `β = γ` at the vacuum. Route 1 of
  `thm:beta-eq-gamma-triple` is exact and definitional; this
  workspace extends it to all scales, not replaces it.

## Success criteria

This workspace closes when:

1. **M2 delivers** `U_eff(φ) = (β/3) φ³ − (γ/4) φ⁴ + higher` with
   all higher terms argued irrelevant at the WF fixed point (closes
   **OP-1**).
2. **M3 delivers** `C_β / C_γ = 1.74 ± δ` from MK-RG with δ ≤ 5%,
   with MC cross-check at `L ≥ 32` (closes **OP-2** Route 3 and
   hence all-scale protection via WF universality).
3. **M4 delivers** `γ (J, λ₀, m₀²)` at the WF fixed point, matched
   to `g₀ᵉ = 0.86941` within the known TGP calibration error
   (closes **OP-4**).

A negative result (e.g. `C_β / C_γ` converges to a value
incompatible with 1.74, or higher operators are not irrelevant) is
equally acceptable as an outcome and would be documented as the
analogue of M2-c/M3-a in the op6 programme.

## Files in this folder

- `README.md` — this file.
- `plan_of_attack.md` — top-level strategy and risk assessment.
- `M1_potential_inventory.md` — catalog of U(φ) across TGP corpus.
  *Stub; to be completed.*
- `M2_effective_potential_sketch.md` — analytical plan for
  H_Γ → U_eff via Hubbard-Stratonovich + saddle-point.
  *Stub; to be completed.*
- `M3_CbCg_numerical_plan.md` — numerical plan for C_β / C_γ via
  extended MK-RG. *Stub; to be completed.*
