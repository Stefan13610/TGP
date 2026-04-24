# Plan of attack — OP-1 / OP-2 / OP-4 potential cluster

**Date opened:** 2026-04-24.
**Successor to:** OP-6 (`research/op6/`), which closed at the axiom
level via the v2 pivot.

---

## 1. Target

Derive the self-interference potential

    U_eff(φ) = (β/3) φ³ − (γ/4) φ⁴ + (irrelevant)                (⋆)

from the v2 substrate Hamiltonian

    H_Γ = Σ_i [ π²_i / 2μ + (m₀²/2) ŝ_i² + (λ₀/4) ŝ_i⁴ ]
          + J Σ_⟨ij⟩ A_ij · ŝ_i² ŝ_j² (ŝ_j² − ŝ_i²)²            (v2)

with Z₂ symmetry `ŝ → −ŝ` and composite field `Φ = ŝ²`,
`φ = Φ / Φ₀`. Equation (⋆) currently lives as a **postulate**
in the core paper; the aim is to make it a **theorem** and in the
same step fix `γ` from microscopic parameters.

This simultaneously closes OP-1 (form), OP-2 (all-scale β = γ) and
OP-4 (value of γ).

## 2. Strategy

The approach is a **two-line derivation**:

- **Line A — analytical (M2):** Hubbard-Stratonovich on Φ = ŝ²
  trades the quartic `ŝ⁴` vertex for a Gaussian integral over an
  auxiliary field σ. Integrating out ŝ and expanding the resulting
  `ln det` around the Z₂-broken saddle `⟨ŝ⟩ = v, Φ = v²` gives an
  effective action `Γ_eff[φ]` whose local part is (⋆). The cubic
  coefficient `β` arises from the anomalous dimension of the
  composite `ŝ²`; the quartic `γ` arises from the one-loop
  fluctuation determinant + tree-level GL bond contribution.
- **Line B — numerical (M3):** extend the existing MK-RG pipeline
  from dodatek B to track the flow of the cubic operator
  (coefficient `B_ℓ`, with `C_β := B_* / v*²`) and the quartic
  operator (coefficient `Γ_ℓ`, with `C_γ := Γ_* / v*²`). Read off
  `C_β / C_γ` at the WF fixed point and test against the prediction
  `C_β / C_γ ≈ 1.74` required for `β = γ` at criticality.

Lines A and B must agree within their error bars. Agreement closes
the cluster; disagreement indicates either an operator we missed in
Line A (open OP-1) or a non-WF fixed point (harder).

## 3. Milestone structure (M1 → M4)

| M | Title | Depends on | Output |
|---|---|---|---|
| M1 | Potential inventory | — | `M1_potential_inventory.md` + parameter dependency graph |
| M2 | Effective potential (analytical) | M1 | `M2_effective_potential_sketch.md` + sympy/Mathematica worksheet |
| M3 | C_β / C_γ (numerical MK-RG) | M1 | `M3_CbCg_numerical_plan.md` + Python pipeline analogous to `m3a_block_rg_1d.py` |
| M4 | γ from microscopic (J, λ₀, m₀²) | M2 + M3 | γ(J, λ₀, m₀²) + calibration check against g₀ᵉ = 0.86941 |

Workflow is M1 → (M2 ∥ M3) → M4. M2 and M3 are independent and
should be run in parallel once M1 clarifies the operator basis.

## 4. Risk register

| Risk | Probability | Mitigation |
|---|---|---|
| `U_eff` has a non-trivial sextic or higher term that is NOT irrelevant at WF | low | check sextic anomalous dimension against bootstrap data (Δ_{ε²} ≈ 3.83 in 3D Ising → irrelevant) |
| `C_β / C_γ` converges but to a value ≠ 1.74 | medium | would indicate `β ≠ γ` at criticality; would force re-opening of OP-2 Route 1 vs Route 3 reconciliation |
| MK-RG truncation underestimates the cubic coefficient flow | medium | cross-check against FRG (Polchinski) truncation in `core/formalizm/` |
| H-S decoupling introduces spurious cubic term proportional to `⟨ŝ⟩ ≠ 0` choice | low | restrict to Z₂-symmetric phase expansion around `⟨ŝ⟩ = 0`, then Legendre-transform to 1PI |
| The WF fixed point is not the relevant IR fixed point for the v2 GL bond | medium | verify by scaling-dimension analysis analogous to M3-c of op6; the v2 bond `ŝ_i² ŝ_j² (Δŝ²)²` is Z₂-even in ŝ and flows to 3D Ising universality class |
| Numerical MK-RG is run on `H_Γ` v2 but the effective potential at low ℓ picks up contributions from modes outside the truncation | high | document truncation level explicitly; follow op6 M2-b pattern of validating against exact sub-limits |

## 5. Methodological anchors

- **Line A uses:** Hubbard-Stratonovich (`Φ = ŝ²`), Legendre
  transform to 1PI, saddle-point expansion around the broken-Z₂
  solution of `U_eff`. Reference: Zinn-Justin §5, 2002 ed.
- **Line B uses:** Migdal-Kadanoff RG on the `d = 3, z = 6`
  hypercubic lattice, extended with explicit cubic and quartic
  operators. Reference: dodatek B §app:B-MK (current code), with
  extension to operator mixing matrix.
- **Consistency:** Both lines must sit inside the 3D Ising WF
  universality class. Scaling-dimension checks use the conformal
  bootstrap data already used in op6 M3-c.

## 6. Deliverables

At closure:

1. **Paper patch:** promote core-paper postulate on `U(φ)` to
   Theorem; append a short appendix (dodatek B-II or
   dodatek O) documenting the M2 derivation.
2. **Core paper table:** OP-1, OP-2, OP-4 rows move from `HY` to
   `TH`.
3. **KNOWN_ISSUES.md** entry documenting the closure (analogous to
   the 2026-04-24 v2 pivot entry for OP-6).
4. **Zenodo v2.1.0 release** of tgp-core-paper with updated theorem
   layer.

## 7. Timeline (soft)

- M1: 1 session (inventory + dependency graph).
- M2: 2-3 sessions (H-S + saddle-point + operator analysis).
- M3: 2-3 sessions (extend MK-RG code + scan).
- M4: 1 session (assemble γ and calibration check).

Not a commitment, just a sanity check.

## 8. Fallbacks

- If M2 stalls on operator enumeration, retreat to numerical-first
  approach: run M3 pipeline with a more complete operator basis and
  read off which operators survive to the fixed point, then
  reverse-engineer Line A.
- If M3 gives `C_β / C_γ ≠ 1.74`, investigate whether the discrepancy
  is fine-structure (e.g. operator mixing) or fundamental
  (e.g. wrong universality class). The latter would force a new
  axiom-level discussion; the former is a quantitative refinement
  inside the current v2 axiomatic.
