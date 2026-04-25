# M7 — GL-bond eigenvalue at the M3 fixed point (verdict)

**Date:** 2026-04-25.
**Test designation:** Track-B sketch of P3.2 of
`external_review_2026-04-25/review_response_plan.md`.
**Verdict:** **STRONGLY IRRELEVANT.** The diagonal projection of
the J_GL response at the M3 fixed point gives **|λ_GL| ≈ 0.07**,
two orders of magnitude below the relevance threshold (|λ| > 1).
Bare J_GL decays under MK; the v2 GL bond operator, linearised at
the 3D Ising WF fixed point, **cannot drive closure of OP-2b on
its own**.

This rules in favour of the third option of M6 §7: M6's
closure-in-principle at J_GL ≈ +5.89 was indeed an artefact of
treating J_GL as a fixed parameter at unphysically large bare
value. Under proper RG flow, J_GL flows to zero in the IR.

Setup: `M7_glbond_eigenvalue.md`.
Code:  `mk_rg_glbond_eig.py`.
Raw:   `mk_rg_glbond_eig_results.txt`.

---

## 1. Numerical result

### 1.1 M3 fixed point (n_ops=8, n_quad=1200, s_max=10)

```
r* = -2.45694    u* = +3.01611    B* = -5.12163    G* = +9.00600
B*/Γ* = -0.5687     v*² = 0.81461       1/v*² = +1.2276
K_new at M3 FP = +5.377
```

(Reproduces M3 baseline ✓.)

### 1.2 Eigenvalue estimates

Three independent extractions of the J_GL eigenvalue from one MK
step at the M3 FP:

| Method | λ_GL | Stability |
|---|---|---|
| **L² norm projection onto O_GL** | **+0.0704** | extremely stable: identical to 4 digits across n_outer ∈ {20,30,40,50} and n_max ∈ {8,12,16} |
| Monomial-basis extraction (2,6) coeff | +0.028 to +0.094 | varies with basis truncation n_max |
| Numerical-perturbation finite-diff | +0.003 to +0.009 | unreliable: breakdown of linear approximation at large s |
| Canonical scaling 4/K_new^4 | +0.0048 | tree-level baseline |

The **L² norm projection is the most trustworthy** (no basis
truncation error, no finite-difference cancellation). All
estimates agree on the qualitative answer:

```
|λ_GL| < 0.1   ≪   relevance threshold = 1.
```

### 1.3 Operator-mixing diagnostic

The shape ratio c(4,4) / c(2,6) of the linear response on the
degree-8 bond subspace:

```
GL-bond shape (1 : -2 : 1)   →  c(4,4)/c(2,6) = -2.
Numerical (n_outer=30, n_max=12): -3.52.
Numerical (n_outer=50, n_max=16): -2.08.
```

The ratio approaches −2 as the basis is enlarged (n_max=16), but
remains slightly off. **Conclusion:** the J_GL → J_GL response is
*almost* pure-GL at the M3 FP, with sub-leading mixing into the
"(s_1+s_3)^8 cross" operator. The diagonal eigenvalue computed via
the L² projection is therefore a faithful single-number summary
of the GL relevance.

### 1.4 Bilateral-symmetry sanity

`|c(2,6) − c(6,2)| / max ≈ 1×10⁻⁸` across all robustness configurations.
The numerical 2D quadrature respects bilateral symmetry to machine
precision. ✓

## 2. Decision against M7 §6 criteria

| Criterion | Threshold | Result | Verdict |
|---|---|---|---|
| Relevant | |λ_GL| > 1 | 0.07 | NO |
| Slow irrelevant | 0.5 < |λ_GL| < 1 | 0.07 | NO |
| **Strongly irrelevant** | |λ_GL| < 0.5 | 0.07 | **YES** |

J_GL is not just irrelevant but **strongly irrelevant**: bare J_GL
decays by a factor of ~14 per MK step (1 / 0.07) at the M3 fixed
point. After ~3 MK steps starting from O(1) bare J_GL, the
coupling is at the level of the floating-point noise.

## 3. Comparison with M3, M4, M5, M6

| Test | Channel | Result | Closure J/μ/η |
|---|---|---|---|
| M3 | bilinear bond only | gap 1.796 | — |
| M4 | + H-S Jacobian μ ln(s² + ε²) | NEGATIVE | none |
| M5 | + Z_Φ via η-deformation | NEGATIVE | none |
| **M6 Track A** | + GL bond, J_GL fixed | closure-in-principle | J_GL ≈ +5.89 (outside perturbative regime) |
| **M7 (this)** | GL eigenvalue at M3 FP | **|λ_GL| ≈ 0.07** | **strongly irrelevant** |

The picture is now consistent:

- All three single-channel candidates (H-S, Z_Φ, GL) **fail** to
  close OP-2b at the level of single-operator additions to
  single-site MK-RG.
- M6's closure-in-principle for GL was, as suspected, an artefact
  of large bare J_GL. Under proper RG flow, J_GL is irrelevant
  and decays to zero.
- The OP-2b gap at the M3 FP is **genuinely** open at the level
  of single-site MK-RG; closing it requires either:
  - an interacting Jacobian incorporating multiple irrelevant
    operators with strong off-diagonal mixing (NPRG-level
    treatment), or
  - a fundamentally different channel not captured by the M3 §6
    enumeration.

## 4. What this proves (and does not prove)

**It does prove:**

- The GL-bond operator's diagonal eigenvalue at the M3 FP is
  |λ_GL| ≈ 0.07, far below the relevance threshold.
- M6's closure-in-principle finding cannot be lifted to a true
  closure under proper RG flow: bare J_GL flows to zero in the IR.
- The v2 GL bond, treated as a perturbation of the 3D Ising
  WF FP at the level of single-site MK-RG, is **not** the missing
  physics for OP-2b.
- The M3 → M3 + GL bond extension does not produce a relevant
  direction in coupling space at the linearised level.

**It does not prove:**

- That J_GL is irrelevant at *every* fixed point of the extended
  flow. If the extended flow has a *different* (non-Wilson-Fisher)
  fixed point, J_GL there might be relevant. This M7 only addresses
  the M3 FP.
- That the OP-2b gap is unconquerable. Off-diagonal mixing in the
  full operator basis (J_GL ↔ "(s_1+s_3)^8 cross" ↔ ...) could
  produce a relevant eigenvector that is *not* aligned with O_GL
  individually. M7 reports only the diagonal element.
- That single-site MK-RG is the right framework. NPRG (Wetterich)
  resolves these ambiguities and is the gold-standard cross-check.

## 5. Resolution of M6's caveat

M6 left open: "Track A's first-order J_GL truncation is not
trustworthy near J_GL ≈ 5.89; full Track B (J_GL flowing) or
NPRG required to confirm."

M7 (Track B sketch) provides a partial Track-B answer at the
linearised level: **J_GL is irrelevant** at the M3 FP, hence the
"closure" found in M6 is not realised under proper RG flow. The
caveat in M6 §6 ("closure-in-principle only") is thus **resolved
in the negative**: the GL bond does not, as a single operator
flowing under MK, close OP-2b.

## 6. What's next

With all three single-channel candidates exhausted:

1. **NPRG (Wetterich) cross-check (P3.4).** The non-perturbative
   gold standard. Confirms or denies M3-M7 at the level of full
   operator basis with all non-perturbative resummations. This
   is now the recommended next step for OP-2b.

2. **Full operator-basis Jacobian at the M3 FP (extended Track B).**
   Compute the full N×N matrix of eigenvalues for the operator
   sector {bond_GL, bond_8_cross, bond_K^2, ...}. If any eigenvector
   has |λ| > 1 with non-negligible projection on O_GL, GL is
   *indirectly* relevant via mixing.

3. **Document the OP-2b gap as a genuine open problem.** M3-M7
   together constitute strong evidence that single-channel single-
   site MK-RG cannot close OP-2b. Update the v2 paper to reflect
   this: OP-2b is "open at the level of single-site MK-RG; closure
   requires NPRG or operator-mixing analysis."

The recommendation is to do **(3) first** (paper update) and
**(1) next** (NPRG). Option (2) is a follow-up if NPRG itself is
inconclusive.

## 7. Files

- `mk_rg_glbond_eig.py` — M7 implementation (eigenvalue extraction
  via 2D F(s_1, s_3) projection).
- `mk_rg_glbond_eig_results.txt` — full output:
    - validation: M3 FP reproduced.
    - L² norm projection: λ_GL_norm = 0.07043 (stable across
      all robustness configurations).
    - monomial-basis extraction: λ_GL_avg ∈ [0.028, 0.094]
      (basis-dependent).
    - numerical-perturbation cross-check: unreliable due to
      large-s breakdown.
    - canonical baseline: 4/K_new⁴ = 0.0048.
- `M7_glbond_eigenvalue.md` — analytical setup.
- `M7_results.md` — this verdict.

## 8. Bottom line

**Track-B (linearised, J_GL eigenvalue at M3 FP):** the GL coupling
is **strongly irrelevant** (|λ_GL| ≈ 0.07). M6's closure at
J_GL ≈ +5.89 is **not** realised under proper RG flow.

**Status of OP-2b after M3, M4, M5, M6, M7:**

- M3 baseline: gap −1.80.
- M4 (H-S Jacobian): no closure for any μ.
- M5 (Z_Φ via η-deformation): no closure for any η.
- M6 (GL bond, fixed J_GL Track A): closure at unphysical
  J_GL ≈ +5.89 (outside perturbative regime).
- **M7 (GL eigenvalue, Track B): J_GL is strongly irrelevant
  (|λ_GL| ≈ 0.07).**

**Conclusion:** the OP-2b gap is **open** at the level of
single-channel single-site MK-RG. All three single-channel
candidates from the M3 §6 list have been falsified. Closure
requires either NPRG (P3.4) or operator-mixing analysis at the
extended Jacobian level.
