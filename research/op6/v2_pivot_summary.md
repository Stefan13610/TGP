# OP-6 v2 pivot — summary of the axiom change

**Date:** 2026-04-24.
**Status:** Axiom-level change committed to tex source. v2 release
pending companion-paper audit.
**Master record:** `KNOWN_ISSUES.md` 2026-04-24 entry.

---

## 1. One-sentence summary

The substrate Hamiltonian in `eq:H-Gamma` of the core paper has
been changed from a **bilinear Z₂-odd bond** (v1) to a
**Ginzburg–Landau gradient bond in Φ = ŝ²** (v2), closing OP-6
at the axiom level after the v1 derivation path was ruled out
on three independent grounds.

## 2. The axiom change

### v1 (published, Zenodo 10.5281/zenodo.19670324)

    H_Γ = Σ_i [π²/2μ + (m₀²/2) ŝ_i² + (λ₀/4) ŝ_i⁴]
          − J Σ_⟨ij⟩ A_ij · ŝ_i ŝ_j            ← bilinear bond

### v2 (tex source, release pending)

    H_Γ = Σ_i [π²/2μ + (m₀²/2) ŝ_i² + (λ₀/4) ŝ_i⁴]
          + J Σ_⟨ij⟩ A_ij · ŝ_i² ŝ_j² · (ŝ_j² − ŝ_i²)²   ← GL bond

with J > 0, λ₀ > 0, A_ij ∈ {0, 1}.

### Continuum form (Prop. `prop:substrate-action` in dodatek B)

    F_kin^geo[φ]  =  ∫ d^dx  (K_geo/2) · φ⁴ · (∇φ)²,
    K_geo        =  2d · J · a_sub^{2−d} .

## 3. Why we pivoted

The v1 derivation path claimed that Migdal–Kadanoff coarse-graining
of the bilinear bond generates the geometric coupling
`K_ij = J(φ_iφ_j)²`, thereby supplying (C1)–(C3) for Theorem
`thm:alpha2`. Between 2026-04-22 and 2026-04-24 this claim was
tested on three independent grounds:

| Probe | Technique | Outcome |
|---|---|---|
| M2-a | Analytical RG (4 mechanisms) | All fail |
| M2-b | 1D envelope-RG MC | p=−0.28±0.43, **3σ reject** of p=+1 |
| M3-a | 1D block-RG MC, B∈{1..16} | p=−0.52±0.28, **5.5σ reject** of p=+1 uniformly |
| M3-c | Scaling dimensions, 3D Ising WF (conformal bootstrap) | K^{(1)}/K^{(0)} ~ k^{1.41} → 0 |

Minimal Z₂-preserving extensions of v1 H_Γ (tricritical,
higher-derivative, O(N), long-range, two-field, gauged Z₂) were
also evaluated and shown not to escape the obstruction (M3-c §6).

## 4. What is preserved

- **Ontological chain "nothingness + Z₂ → reality".** The
  microscopic amplitude ŝ_i remains Z₂-odd; the Z₂-even composite
  Φ_i = ŝ_i² drives the geometric bond.
- **All four axioms** (ax:N0, ax:substrate, ax:source, ax:P3).
  Only the bond term inside `eq:H-Gamma` changed.
- **On-site structure.** π²/2μ kinetic, m₀² mass, λ₀ φ⁴
  self-interaction — all unchanged.
- **Proposition 2** (Z₂ symmetry of H_Γ). The new GL bond is even
  under ŝ → −ŝ (the bond is built from ŝ², which is Z₂-even), so
  the proposition still holds.
- **Proposition `prop:substrate-action`** in dodatek B. It has
  used the GL form since 2026-03-24; the pivot simply elevates
  this object to the axiom level.
- **Theorem `thm:alpha2`**. Its proof is axiom-independent (it is
  a classification inside (C1)–(C3)), so the proof is unchanged.

## 5. What strengthens under v2

- **"No substrate ⇒ no geometry" as property of the axiom.** The
  GL bond vanishes on any configuration with ŝ_i = 0. In v1, a
  bilinear bond is zero when *both* neighbours are zero, but was
  not explicitly tied to Φ. In v2, the bond is quadratic in
  Φ = ŝ², so empty substrate literally cannot support metric
  bonds. This is ontologically stronger than v1.
- **α = 2 becomes algebraic consequence of axiom.** Via
  Prop. `prop:substrate-action`, K(φ) = K_geo φ⁴ follows in one
  step from the GL bond; the Euler–Lagrange equation gives
  α(φ) = 2 directly. No RG derivation required.
- **`thm:alpha2` acquires a new role.** In v1 it was the *only*
  supplier of α = 2 (through an axiomatic classification without
  a substrate-level derivation). In v2 it becomes a **consistency
  check**: among all local Φ-covariant kinetic operators, the GL
  form is the unique one satisfying (C1)–(C3). So v2 H_Γ is the
  *minimal* choice consistent with substrate-metric geometry.

## 6. What weakens under v2

Honest accounting:

- **The bond is no longer bilinear in ŝ.** It is 8th-order in ŝ
  (4th-order in Φ × (Φ_j − Φ_i)²). This is less "elegant" than
  v1 but more ontologically principled (§5). Similar bonds appear
  in standard Ginzburg–Landau theory for interface stiffness.
- **The "Z₂-odd substrate magnet" image is less direct.** In v1,
  H_Γ was a scalar φ⁴-with-bilinear-bond Ising-like theory. In v2,
  the *on-site* part is unchanged (still Z₂-odd φ⁴), but the bond
  couples the Z₂-even composite Φ. One can still read the theory
  as "Z₂-odd ŝ whose square Φ drives geometry", but it is not a
  pure scalar-field theory in ŝ.
- **Mean-field analysis of the v2 bond is more involved.** The
  eff. potential for Φ acquires contributions from the bond
  fluctuations that differ from the bilinear case. For v2 this
  is not an obstacle (Prop. `prop:substrate-action` handles the
  continuum limit explicitly), but mean-field intuition from v1
  does not transfer directly.

## 7. Files changed

| File | Change |
|---|---|
| `tgp-core-paper/paper/tgp_core.tex` | `eq:H-Gamma` bond term: bilinear → GL. Added 5-line axiom annotation. Extended `rem:alpha2-closure` (7 → ~30 lines) documenting OP-6 closure. Rewrote MK coarse-graining paragraph in §3 to point at Prop. `prop:substrate-action`. |
| `tgp-core-paper/KNOWN_ISSUES.md` | New top section "2026-04-24 — OP-6 closed via axiom pivot". |
| `TGP_v1/research/op6/M1_H_Gamma_canonical.md` | §5 recommendation reversed (M1-B → M1-A′). §6–§8 retained verbatim with a superseded-stamp. |
| `TGP_v1/research/op6/README.md` | M3-a/M3-c entries added; non-goals updated. |

New files:
- `TGP_v1/research/op6/M3c_scaling_dimensions.md` (analytical)
- `TGP_v1/research/op6/M3a_block_rg_1d_plan.md`
- `TGP_v1/research/op6/m3a_block_rg_1d.py`
- `TGP_v1/research/op6/m3a_block_rg_results.txt`
- `TGP_v1/research/op6/m3a_run.log`
- `TGP_v1/research/op6/M3a_results.md`
- `TGP_v1/research/op6/v2_pivot_summary.md` (this file)

Dodatek B was **not** modified (it already used the GL form).

## 8. What remains to do before v2 Zenodo release

1. **Companion-paper audit.** Confirm that `tgp-qm-paper`,
   `tgp-leptons-paper`, `tgp-sc-paper` do not rely on the v1
   bilinear bond form. Expected: none does (all use coarse-
   grained Φ dynamics, not the microscopic bond).
2. **Regenerate PDF + diff.** Produce the updated paper PDF and
   a v1→v2 diff for transparency.
3. **Update the paper abstract and §0 summary** to state the
   axiom clearly.
4. **Update `TGP_v1/README.md`** to reflect v2 status.
5. **Tag a new Zenodo release** with a DOI distinct from v1.
6. **Keep v1 accessible** on Zenodo unchanged; this file + the
   KNOWN_ISSUES note are the reconciliation record.

None of these are blocking for any internal work; they are the
deliverables for an external v2 release.

## 9. Non-goals

- **No retroactive claim that v1 was "correct" in a different
  sense.** v1 had a real gap (OP-6) and we closed it by pivot.
- **No re-opening of the M1-B derivation.** Three independent
  probes agree. Further effort on H₁ → H_GL is a dead end.
- **No changes to (C1)–(C3)**, `thm:alpha2`, or
  `prop:substrate-action` beyond the metadata. The mathematical
  content is unchanged; only its placement in the argument
  changed (from "open derivation" to "axiom").
