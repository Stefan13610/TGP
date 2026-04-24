# OP-6 rework — rigorous continuum limit of H_Γ

This folder is the clean workspace for a fresh attempt at **OP-6**
(the rigorous continuum limit of the TGP substrate Hamiltonian).
It supersedes the earlier `research/continuum_limit/` attempt whose
central synthesis script `a1_alpha2_frg_synthesis.py` has been
withdrawn; see `tgp-core-paper/KNOWN_ISSUES.md` for the record.

## Programme

- **M1** — consolidate all circulating H_Γ definitions into one
  canonical Hamiltonian and pin variable conventions.
  *Status: done.* See `M1_H_Gamma_canonical.md`. Canonical choice
  is H₁ (bilinear bond, core-paper axiom) with the plan to derive
  H₃ (quartic bond) as an effective form via RG (M1-B).
- **M2** — derive the effective kinetic functional
  F_kin[φ] = (K_eff/2) φ⁴ (∇φ)² from the canonical H_Γ by RG
  coarse-graining, and hence α = 2 via Euler–Lagrange.
  *Status: closed with negative result.*
  - **M2-a** (analytical sketch): mechanisms (i)–(iv) in FRG / LPA'
    / derivative expansion / extended operator basis all fail to
    produce K(φ) ∝ φ⁴ from H₁. Only "envelope" coarse-graining
    (mechanism v) survives as an open possibility.
    See `M2a_analytical_sketch.md`.
  - **M2-b** (1D MC numerical test): envelope coarse-graining of
    H₁ in 1D gives K_eff(Φ) ∝ Φ^(−0.28 ± 0.43), **rejecting**
    the TGP target K_eff ∝ Φ at 3σ. Method validated on the
    Gaussian sub-limit against exact analytical prediction.
    See `M2b_results.md`.
  - **M2-c** (reality check, 2026-04-23; closed 2026-04-24):
    audit of what "H₃" actually is. The bilinear bond
    `−J Σ(φ_iφ_j)² = −J Σ Φ_iΦ_j` (i.e. the "H₃" of M1 §1) gives
    α(φ) = 1, NOT α(φ) = 2. The object that actually gives
    α(φ) = 2 is the GL functional H_GL = `Σ J(φ_iφ_j)²(φ_j − φ_i)²`,
    which is what dodatekB's `prop:substrate-action` actually uses.
    Analytical claim confirmed numerically by Part E of the MC
    scan: H₃ bilinear gives p_Φ = −2.49 ± 0.63 (n = 12),
    rejecting the α=2 target (p_Φ = +1) at **5.5σ**. See
    `M2c_H3_reality_check.md`.
  - Consequence: M1-B (derive α=2 from H₁) is numerically
    disfavoured. The pivot target is **not** "H₃ bilinear" (that
    gives α=1); it is **H_GL** (GL effective functional with
    density-dependent stiffness). M1-A is therefore more invasive
    than originally framed, because H_GL is an effective-theory
    axiom, not a microscopic spin Hamiltonian. Two live options
    for v2: (1) adopt H_GL as the axiom (M1-A′); (2) close
    M1-B by showing H₁ → H_GL via RG (harder than M2-b's
    question). Both remain open.
- **M3** — FRG / lattice MC hardening of M2-b's 1D failure.
  - **M3-c** (analytical scaling dimensions, 2026-04-24; closed):
    at the 3D Ising Wilson–Fisher fixed point, scaling-dimension
    data (Δ_ε = 1.41) force `K^{(1)}/K^{(0)} ~ k^{1.41} → 0` in the
    IR. So `K_*(Φ) → const`, not `K_*(Φ) ∝ Φ`. Minimal Z₂-preserving
    extensions (tricritical, higher-derivative, O(N), long-range,
    two-field, gauged Z₂) evaluated — none escape the obstruction.
    Prior on M1-B lowered from < 50% to < 5%. See
    `M3c_scaling_dimensions.md`.
  - **M3-a** (1D block-RG numerical test, 2026-04-24; closed):
    direct block-RG on H₁ in 1D with block sizes B ∈ {1, 2, 4, 8,
    16}. Exponent `p_B` is **RG-invariant** at ≈ −0.52 (drift =
    0.0σ across four octaves), TGP target p = +1 rejected at
    **5.5σ** at every block size. `C_B ∝ 1/B` confirms block-RG
    is working. See `M3a_results.md`.
  - **M3-b** (3D cluster MC): demoted from "blocking" to
    "hardening-only". Not required to advance OP-6.
  *Status: M3-c + M3-a closed; M1-B rejected on three independent
  grounds. M3-b optional.*
- **M4** — effective action for Φ under the M1-A′ pivot
  (continuum form of H_GL is already the direct content of
  `prop:substrate-action`; M4 becomes documentation of the
  axiom switch rather than a fresh derivation).
  *Status: not started.*
- **M5** — Γ-convergence (or Banach contraction) of the block-spin
  map on a simplified model where it can actually be proved.
  *Status: not started.*

## Relation to other files

- `../continuum_limit/cg_strong_numerical.py` remains a legitimate
  numerical probe of block-averaging. We will extend it (or write
  a successor in this folder) for M3.
- `../continuum_limit/a1_alpha2_frg_synthesis.py` is WITHDRAWN.
  Do not cite or re-use its tests.
- `axioms/substrat/dodatekB_substrat.tex` `prop:substrate-action`
  is the formal target: once M2 closes, that proposition becomes a
  theorem instead of a postulate.

## Non-goals

- Re-litigating whether TGP as a physical theory is correct.
- Changing the published core paper axioms before M2 has actually
  closed. Any axiom change is a v2 decision. (M2 has now closed
  with a negative result; M3-a + M3-c have closed the "can block-RG
  save M1-B?" question also negatively, so M1-B is now rejected on
  three independent grounds. The v2 release itself is still a
  separate deliverable and will require a decision between M1-A,
  M1-A′, or a Z₂-preserving minimal extension of H₁ — see
  `M3a_results.md` §7.)

## Files in this folder

- `README.md` — this file.
- `M1_H_Gamma_canonical.md` — inventory of H_Γ variants, variable
  conventions, recommended canonical choice (H₁) and M2 sub-plan.
  *Partially superseded by M2c (see in-file corrective note).*
- `M2a_analytical_sketch.md` — analytical ruling-out of the
  standard RG mechanisms for deriving α=2 from H₁.
- `M2b_numerical_plan.md` — plan for the numerical surrogate of
  the envelope FRG (direct MC measurement of K_eff(Φ)).
- `M2b_results.md` — numerical results: M1-B disfavoured.
  *Pivot target corrected from "H₃" to "H_GL" per M2c; see
  in-file corrective note.*
- `M2c_H3_reality_check.md` — correction: the bilinear "H₃" of
  M1's inventory gives α=1, not α=2. The actual α=2-generating
  object is the GL functional H_GL. Identifies two live v2
  options.
- `m2b_envelope_stiffness_1d.py` — the 1D MC implementation.
  *Part E adds a numerical check that bilinear H₃ gives K(Φ)=const
  (p_Φ ≈ 0), distinct from H_GL.*
- `m2b_scan_results.txt` — raw parameter-scan table from
  `m2b_envelope_stiffness_1d.py`.
- `M3c_scaling_dimensions.md` — analytical scaling-dimension
  argument that at 3D Ising WF, `K_*(Φ) → const` in the IR; closes
  the "can block-RG save M1-B?" question negatively.
- `M3a_block_rg_1d_plan.md` — plan for numerical block-RG test.
- `m3a_block_rg_1d.py` — 1D block-RG implementation.
- `m3a_block_rg_results.txt` — raw data from the block-RG scan.
- `m3a_run.log` — full run log.
- `M3a_results.md` — numerical result: `p_B` is RG-invariant at
  ≈ −0.52 across B ∈ {1, …, 16}; `+1 rejected at 5.5σ`.
