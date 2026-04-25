# M8 — NPRG (Wetterich) cross-check of OP-2b at the 3D Ising WF FP

**Date:** 2026-04-25.
**Test designation:** P3.4 of `external_review_2026-04-25/review_response_plan.md`.
**Goal:** Compute `B*/Γ*` at the 3D Ising WF fixed point using
non-perturbative RG (Wetterich) at the level of LPA with a Litim
regulator, and compare with the MK-RG result `B*/Γ* = −0.5687`
(M3, M7).
**Decision:** if NPRG-LPA also gives `B*/Γ* < 0`, the M3-M7 verdict
that OP-2b is open at the level of single-channel scalar Z₂ theory
is confirmed in the gold-standard non-perturbative scheme. If NPRG
gives the v1-vacuum value `+1.74` (i.e. `+1` for `B*/Γ*` after the
v² rescaling), MK-RG was misleading and OP-2b closes at NPRG level.

Companion artefacts: `nprg_lpa_3d.py`, `nprg_lpa_3d_results.txt`,
`M8_results.md`.

---

## 1. Why NPRG, why now

After M3 (MK-RG baseline, `B*/Γ* = −0.57`), M4 (H-S Jacobian,
negative), M5 (`Z_Φ` via η-deformation, negative), M6 (GL bond
Track A, closure-in-principle at unphysical `J_GL ≈ +5.89`) and
M7 (GL eigenvalue at M3 FP, `|λ_GL| ≈ 0.07` — strongly irrelevant),
all three single-channel candidates from M3 §6 are exhausted. The
question that remains is **methodological**: is the negative
`B*/Γ* ≈ −0.57` a feature of the Z₂ Ising WF FP, or an artefact of
the Migdal–Kadanoff scheme (real-space block + decimation)?

The Wetterich exact RG equation
```
∂_t Γ_k[Φ] = (1/2) Tr [ ∂_t R_k · (Γ_k^{(2)}[Φ] + R_k)^{-1} ]
```
is non-perturbative and exact. Truncating to a local potential
ansatz (LPA) gives the "gold standard" cross-check: any sign of
`B*/Γ*` consistent across MK and Wetterich is a feature of 3D
Z₂-Ising universality, not of the scheme.

## 2. LPA flow equation in d=3 with Litim regulator

LPA ansatz (η = 0, Z = 1):
```
Γ_k[Φ] = ∫ d³x [ (1/2)(∂Φ)² + V_k(Φ) ].
```
Litim regulator: `R_k(q²) = (k² − q²) θ(k² − q²)`. With this choice,
the Wetterich-loop integral reduces to one over a 3-sphere of radius
`k`, and the flow of `V_k(Φ)` becomes (standard result, Litim 2001;
Berges–Tetradis–Wetterich 2002):
```
∂_t V_k(Φ) = (k⁵ / (6π²)) · 1 / (k² + V_k''(Φ))                 (1)
```
In dimensionless variables `ρ̃ = Φ²/(2 k^{d−2}) = Φ²/(2k)` (d=3),
`v(ρ̃) = V_k(Φ)/k^d = V_k/k³`, the flow becomes (with `t = ln(k/Λ)`):
```
∂_t v(ρ̃) = −3 v(ρ̃) + ρ̃ W(ρ̃) + (1/(6π²)) / (1 + m²(ρ̃))         (2)
```
where
```
W(ρ̃) = v'(ρ̃),    m²(ρ̃) = W(ρ̃) + 2 ρ̃ W'(ρ̃) = ∂²V/∂Φ² (in dim. units).
```

The dimensionless **fixed-point equation** is `∂_t v_* = 0`:
```
3 v_*(ρ̃) = ρ̃ W_*(ρ̃) + (1/(6π²)) / (1 + m_*²(ρ̃))                (3)
```

## 3. Polynomial truncation

Expand around `ρ̃ = 0`:
```
v_*(ρ̃) = Σ_{k=1}^N a_k ρ̃^k.                                   (4)
```
(The constant term `a_0` is fixed by `(3) at ρ̃=0`: `3 a_0 =
1/(6π² (1+a_1))`. We discard it as a non-universal vacuum offset.)

From (4):
```
W_*(ρ̃) = Σ_{k=1}^N k a_k ρ̃^{k−1},
m_*²(ρ̃) = Σ_{k=1}^N k(2k−1) a_k ρ̃^{k−1}                       (5)
        = a_1 + 6 a_2 ρ̃ + 15 a_3 ρ̃² + 28 a_4 ρ̃³ + 45 a_5 ρ̃⁴ + ...
```

Define `[x]_k` ≡ coefficient of `ρ̃^k` in `x`. Order-by-order, eq. (3)
gives, for `k = 1, …, N`:
```
(3 − k) a_k = (1/(6π²)) · [1/(1 + m_*²(ρ̃))]_k.                 (6)
```
At `k=3` (canonical sextic), the LHS vanishes; the equation becomes
a constraint on `a_1, a_2, a_3, a_4`. This is the "marginal
direction" — Wilson–Fisher's fixed point sits where the constraint
is solved by a non-trivial set of lower coefficients.

The system (6) is N coupled algebraic equations in N unknowns
`(a_1, …, a_N)`, solved by Newton iteration.

## 4. Mapping to MK-RG observables

MK-RG convention: `V(ŝ) = Σ_k c_{2k}/(2k) · ŝ^{2k}` with `(r, u, B, Γ)
= (c_2, c_4, c_6, c_8)`. Under the substitution `ŝ² = 2 ρ̃`:
```
V(ŝ) = Σ_k (c_{2k}/(2k)) · (2ρ̃)^k = Σ_k a_k ρ̃^k,
       a_k = c_{2k} · 2^{k−1} / k.                             (7)
```
Inverting:
```
r  = a_1,
u  = a_2,
B  = (3/4) a_3,                                                (8)
Γ  = (1/2) a_4.
```
The convention-invariant ratio:
```
B / Γ  = (3/2) · (a_3 / a_4).                                  (9)
```
MK-RG (M3, N_ops=8): `B*/Γ* = −0.5687`, hence
```
a_3 / a_4   =   (2/3) · (−0.5687)   =   −0.3791                (target)
```
This is the number NPRG-LPA must reproduce (in sign and order of
magnitude) to confirm the M3-M7 verdict.

## 5. Validation: critical exponent ν

Linearise (2) around the FP `v_*`:
```
δv(ρ̃, t) = e^{−λ t} g(ρ̃),        ν = 1/λ_relevant.            (10)
```
The relevant eigenvalue corresponds to deformations of `a_1`
(thermal direction). At LPA with Litim regulator the literature
value is `ν_LPA ≈ 0.6496` (Litim 2001, also Berges-Tetradis-Wetterich
2002 §A). Our polynomial truncation is required to reproduce
`ν ∈ [0.62, 0.66]` for the FP solution to be considered correct;
this is the validation gate before reading off `a_3/a_4`.

## 6. Polynomial-vs-MK-RG comparison: what is fair

MK-RG's `B*/Γ* = −0.5687` came from N_ops = 8 polynomial truncation
in the `(r, u, B, Γ, …, c_{16})` basis. The natural NPRG-LPA
counterpart is polynomial truncation of `v_*(ρ̃)` to degree N
(i.e. up to `a_N ρ̃^N`). Since MK-RG used 4 even powers up to ŝ⁸
(equivalent to ρ̃⁴), the matching NPRG truncation is **N = 4** (with
robustness scan to N = 6, 8, 10, 12 for convergence assessment).

Both schemes are polynomial truncations; the difference is the
loop kernel:
- MK-RG: bond-move (K → 4K) + decimation (Gauss–Legendre quadrature
  of `exp(−V(ŝ_2) + h ŝ_2)` in source-shifted variable).
- NPRG-LPA: Wetterich-Litim, `1/(1 + W + 2ρ̃ W')` integrated over
  3-sphere.

If the two schemes agree on the **sign** of `B*/Γ*` at matched
polynomial order, then `B*/Γ* < 0` is a feature of 3D Z₂-Ising and
not of either scheme. If they disagree, MK-RG is misleading and the
v1 prediction is restored.

## 7. Decision matrix

| Outcome at N=8 | sign of `a_3/a_4` | verdict on OP-2b |
|---|---|---|
| `a_3/a_4 ∈ [−0.5, −0.2]` | negative, ~MK-RG | **OP-2b confirmed open**, M3-M7 universality-class feature |
| `a_3/a_4 ∈ [−0.1, +0.1]` | near-zero | inconclusive; need higher truncation or LPA' |
| `a_3/a_4 ∈ [+0.5, +1.0]` | positive, ~v1 vacuum | MK-RG was wrong, OP-2b closes at NPRG |
| `a_3/a_4 ∉ [−1, +1]` | extreme | numerical instability or wrong FP branch |

Also: `ν_LPA ∈ [0.62, 0.66]` is required as a validation gate for
the FP solution.

## 8. Caveats and what M8 does NOT do

1. **LPA only**, no `Z_k` flow (anomalous dimension `η = 0`). In
   3D Ising the actual `η ≈ 0.036` is small, and we expect `B*/Γ*`
   to be insensitive to `η` at the percent level. LPA' could be a
   future refinement (M9 if needed).
2. **Symmetric-phase polynomial expansion** around `ρ̃ = 0`. The
   broken-phase expansion around the FP minimum `ρ̃_0 > 0` is an
   alternative scheme; literature suggests both converge to the same
   FP function. We use symmetric expansion to match MK-RG's
   small-ŝ Taylor.
3. **No GL-bond operator** in NPRG. Justified: M7 showed J_GL is
   strongly irrelevant at the M3 FP, hence J_GL = 0 at the NPRG FP
   too; the OP-2b ratio at the FP is independent of J_GL.
4. **Litim regulator only.** Different regulators give different LPA
   `ν` values (`ν ≈ 0.65` vs `ν ≈ 0.69` for sharp cutoff). For
   the **ratio** `a_3/a_4`, regulator-dependence is expected to be
   smaller. The Litim choice is most economical (closed-form loop
   integral).

## 9. Implementation outline (`nprg_lpa_3d.py`)

```python
def fp_residuals(a, N):
    # a = (a_1, ..., a_N), N coupled FP equations
    # Build m²(ρ̃) as polynomial in ρ̃
    # Compute 1/(1+m²(ρ̃)) via formal polynomial inversion to order N
    # Form residuals R_k = (3-k)*a_k - c * coeff_k_inv  (k=1..N)
    return R

def solve_fp(N, a_init=None):
    # Newton/scipy.optimize.fsolve from Gaussian initial guess
    # a_init: from N-1 truncation, padded with zeros
    return a_star, ν, eigenvalues

def main():
    # Validation: N=2 should give a_1 ≈ -0.18, a_2 ≈ ... + ν_LPA ≈ 0.65
    # Production: N=4, 6, 8, 10, 12
    # Output: a_3/a_4, hence B*/Γ*
    # Robustness: vary regulator scale (sanity)
```

## 10. Files

- `M8_NPRG_setup.md` — this file (analytical setup).
- `nprg_lpa_3d.py` — polynomial-truncation LPA FP solver.
- `nprg_lpa_3d_results.txt` — raw output: FP coefficients, ν,
  `B*/Γ*` per truncation, robustness scan.
- `M8_results.md` — verdict.

## 11. Bottom line of M8 (to be written)

If `(a_3/a_4)_NPRG ≈ −0.38` (matching MK-RG): **OP-2b confirmed
open at the level of single-component scalar Z₂ field theory.** The
gap is a feature of the WF FP, not of the MK scheme. Closure
requires multi-component or non-scalar physics (e.g. coupling to
the conjugate momentum sector, OP-7 tensor sector, or genuinely
new substrate dynamics).

If `(a_3/a_4)_NPRG ≈ +0.67` (matching v1 vacuum): MK-RG was
misleading; OP-2b closes structurally; revisit M3-M7 verdicts.

The result of §10 will write itself in `M8_results.md` once the
solver is run.
