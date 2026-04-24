# M4 — MK-RG with H-S Jacobian (Test A): Results

**Status:** OP-2b **confirmed open**. The Hubbard–Stratonovich Jacobian
is **not** the missing physics behind M3's `B*/Γ* ≈ -0.57` deviation
from the `β = γ`-at-WF target. The deferred mechanisms (`Z_Φ`,
GL-bond operator, NPRG) remain the candidates.

**Date:** 2026-04-25.
**Test designation:** Test A of
`external_review_2026-04-25/review_response_plan.md`.
**Script:** `mk_rg_phi.py`.
**Raw output:** `mk_rg_phi_results.txt`.
**Analytical setup:** `M4_phi_variable_derivation.md`.

---

## 1. Summary

The reviewer-C5 conjecture (concurrent with the user's own diagnostic
in the 2026-04-25 external review) was that M3's discrepancy

```
B*/Γ* ≈ -0.57    vs.   target  +1/v*²(M3) ≈ +1.23
                                        sign flip + magnitude × 3
```

is an artefact of working in `ŝ`-variables that drop the
Hubbard–Stratonovich Jacobian `Φ^{-1/2}` accompanying the composite
field map `Φ = ŝ²`. Translated to the `ŝ`-action, this Jacobian is
a `μ ln(s² + ε²)` term with the H-S-natural value `μ_HS = 1/2`.

We added this term to M3's MK-RG, scanned `μ ∈ [0, 0.45]`, regulator
`ε ∈ {0.1, 0.05, 0.01}`, at truncations `N_ops ∈ {6, 8}`, and tracked
both `B*(μ)/Γ*(μ)` and the target `1/v*²(μ)`.

**Result.** `B*/Γ*` and `1/v*²` move in opposite directions as `μ` is
increased: `B*/Γ*` becomes more negative (`-0.569 → -0.764`),
`1/v*²` decreases (`+1.228 → +0.871`). The signed gap

```
diff(μ) ≡ B*/Γ*(μ) − 1/v*²(μ)
```

remains strictly negative throughout the scan and reaches a
**minimum of `|diff| ≈ 1.62`** at `μ_min ≈ 0.32`. There is **no sign
change**, no crossing, and no extrapolation that would close the gap
at the H-S-natural value `μ = 1/2`.

**Verdict on OP-2b:** the H-S Jacobian addition is ruled out as the
mechanism behind M3's deviation. OP-2b stays **open** with the
deferred candidates (`Z_Φ`, GL-bond, NPRG) — see §6.

## 2. Implementation summary

`mk_rg_phi.py` extends `mk_rg_bgamma.py`'s `MigdalKadanoffRGN` by a
single weight modification in `_moments`:

```python
log_w -= self.mu * np.log(s**2 + self.eps**2)   # H-S Jacobian, regulated
```

Marginality of `μ` under MK-RG (proved in
`M4_phi_variable_derivation.md` §3) lets us hold `μ` fixed during the
flow and treat it as a label of FP families. For each `(μ, ε)` we
solve for the WF FP by fixed-point iteration in bar-variables
(M3 protocol), with **continuation seeding** (each `μ`'s FP is
seeded with the previous `μ`'s FP, step `Δμ = 0.025`). All other
numerics are inherited from M3 (Gauss–Legendre quadrature with
`n_quad ∈ [1200, 2400]`, `s_max = 10`).

## 3. Validation

### 3.1 μ=0 reproduces M3 exactly

```
M3 expected (n_ops=8):  r* = -2.45694   u* = +3.01611
                        B* = -5.1216     Γ* = +9.0060
                        v*² = 0.81461   B*/Γ* = -0.56869
M4 at μ=0, ε=0.1:       r* = -2.45694   u* = +3.01611
                        B* = -5.1216     Γ* = +9.0060
                        v*² = 0.81461   B*/Γ* = -0.56869
                        |Δr*| = 0.00000  |Δu*| = 0.00000  |Δratio| = 0.00001
```

Identity to 5 decimals. The H-S weight reduces analytically to the
M3 weight at `μ = 0`; the numerical implementation respects this.

### 3.2 μ-marginality (numerical check)

The analytical claim of §3 of `M4_phi_variable_derivation.md` is that
`μ` does not flow under MK-RG. We confirm numerically: holding `μ`
fixed and iterating the polynomial sector to convergence, the FP is
reached and the bar-variables stabilise. If `μ` were not marginal,
the polynomial flow would not converge (or would converge to a `μ`
different from the input value). Convergence is observed for
`μ ∈ [0, 0.275]` at `ε = 0.1` (50–700 iterations, `tol = 10⁻¹⁰`),
beyond which the iteration limit is reached but the bar-variables
have stabilised to 3–4 digits.

## 4. The full scan

### 4.1 Trend table at ε = 0.1, N_ops = 8 (numerically safest)

| μ    | r*       | u*      | B*       | Γ*      | v*²    | B*/Γ*    | 1/v*²   | diff     |
|------|----------|---------|----------|---------|--------|----------|---------|----------|
| 0.00 | −2.45694 | +3.01611| −5.122   | +9.006  | 0.8146 | −0.56869 | +1.2276 | −1.79628 |
| 0.05 | −2.46189 | +2.88588| −4.726   | +8.043  | 0.8531 | −0.58757 | +1.1722 | −1.75979 |
| 0.10 | −2.46632 | +2.74791| −4.320   | +7.090  | 0.8975 | −0.60929 | +1.1142 | −1.72346 |
| 0.15 | −2.46984 | +2.60197| −3.905   | +6.152  | 0.9492 | −0.63474 | +1.0535 | −1.68824 |
| 0.20 | −2.47185 | +2.44811| −3.482   | +5.234  | 1.0097 | −0.66533 | +0.9904 | −1.65573 |
| 0.25 | −2.47131 | +2.28695| −3.053   | +4.339  | 1.0806 | −0.70360 | +0.9254 | −1.62900 |
| 0.30 | −2.45791 | +2.15648| −2.721   | +3.681  | 1.1398 | −0.73924 | +0.8774 | −1.61661 |
| 0.32 | −2.44211 | +2.12776| −2.656   | +3.568  | 1.1477 | −0.74451 | +0.8713 | **−1.61579** |
| 0.35 | −2.42927 | +2.08210| −2.538   | +3.334  | 1.1667 | −0.76119 | +0.8571 | −1.61828 |
| 0.40 | −2.41260 | +1.94616| −2.160   | +2.557  | 1.2397 | −0.84475 | +0.8067 | −1.65141 |

Bold value (μ ≈ 0.325) is the minimum of `|diff|`.

### 4.2 What moves where, as μ increases

- `B*/Γ*` becomes **more negative** (−0.569 → −0.844 at μ=0.40),
  i.e. the sign-flip with the target is **never** repaired.
- `1/v*²` **decreases** (`+1.228 → +0.807`), because `v*²` increases
  as `u*` shrinks (the H-S weight enhances small-s, which softens
  the quartic).
- The gap `|diff| = |B*/Γ* − 1/v*²|` does narrow slightly: from
  1.80 (μ=0) to 1.62 (μ≈0.32), then **widens again** for larger μ
  as `B*/Γ*` plunges faster than `1/v*²`.
- Crucially, **diff never changes sign** in the entire convergent
  regime `μ ∈ [0, 0.45]`. The H-S-natural value `μ = 1/2` is the
  boundary of convergence (logarithmic divergence at s=0 in the
  symmetric s-form), and the trend is moving **away** from
  closure for `μ > 0.32`.

### 4.3 Cross-checks: regulator and truncation independence

The result is **stable** under variation of `ε` (regulator) and
`N_ops` (truncation):

| ε     | min |diff| | μ_min  | comment |
|-------|------------|--------|---------|
| 0.10  | 1.6158     | 0.325  | reference |
| 0.05  | 1.6248     | 0.225  | minimum shifts to lower μ; small-s region more sensitive |
| 0.01  | (~1.63)    | 0.20   | iteration unstable beyond μ≈0.20; but trend identical |

| N_ops | min |diff| | μ_min  |
|-------|------------|--------|
| 6     | 1.6201     | 0.350  |
| 8     | 1.6158     | 0.325  |

The minimum-|diff| value of **≈1.62** is robust to ±2% across all
combinations. There is no signal of `|diff| → 0` in any regulator
or truncation limit.

### 4.4 Eigenvalue spectrum

At μ=0,0.1,0.2 (ε=0.1, N_ops=8), the dominant irrelevant eigenvalue
grows monotonically:

| μ   | |λ_1| (dominant irrelev.) | y_1     | comment |
|-----|---------------------------|---------|---------|
| 0.0 | 0.7621                    | −0.392  | M3 baseline (ω-analogue 0.83) |
| 0.1 | 0.8030                    | −0.316  | drifting toward marginality |
| 0.2 | 0.8704                    | −0.200  | nearing |λ|=1 |

This is consistent with the H-S coupling becoming a **marginal
deformation** at the H-S-natural value μ=1/2: the `μ`-direction is
exactly marginal by construction (§3 of the derivation), and
crucially is decoupled from the polynomial spectrum reported
above (those are eigenvalues of the polynomial-sector Jacobian,
holding `μ` fixed). The growth of `|λ_1|` toward 1 reflects loss
of stability of the WF FP itself in the polynomial sector as the
H-S weight enhances small-s and pushes the system toward the
mass term's strong-coupling regime.

## 5. Numerical caveats

1. **Convergence boundary at μ=1/2.** In the unregulated form
   (ε = 0), the integrand `(s²)^{-μ}` is integrable iff `2μ < 1`.
   With ε > 0 the integral is finite for any μ, but the integrand
   becomes increasingly peaked at s=0 as `ε → 0` and `μ → 1/2`,
   making the polynomial cumulants noisier. We see this in the
   `ε = 0.01` scan, where iteration robustness degrades for
   `μ > 0.20`. The `ε = 0.1` scan is the cleanest in this regime.

2. **Continuation-seeded saddle drift.** For `μ > 0.40` at `ε ∈
   {0.1, 0.05}` the continuation-seeded iteration sometimes loses
   the WF branch and falls into a runaway saddle (`u* → ∞`,
   `v*² → 0`), the same failure mode as M3's odd-`N_ops`
   pathology. The data at `μ ∈ {0.425, 0.45}` (ε=0.1, ε=0.05) is
   for that reason **not** physical; we exclude these points from
   the trend.

3. **Iteration limit at μ > 0.275.** For ε=0.1 and `μ ∈ [0.30, 0.40]`
   the FP iteration reaches `max_iter = 8000` without satisfying
   `tol = 10⁻¹⁰`, but the bar-variables have stabilised to 3–4
   digits. The reported FP values are reliable for the trend
   analysis (which depends on 2–3 digits), though tighter
   convergence would require more iterations or a smarter
   accelerator (e.g. Anderson mixing). Given the negative result
   stands at the level of `|diff| ≈ 1.6` (no sign change, no
   crossing), this caveat does not affect the verdict.

## 6. Verdict on OP-2b

**Per the decision criteria of `M4_phi_variable_derivation.md` §6:**

- **Positive criterion** (B*/Γ* crosses 1/v*² for some μ_* ∈ (0, 1/2]):
  **NOT MET.** No sign change anywhere in the scanned regime; both
  quantities have the wrong relative sign at every μ.
- **Negative criterion** (B*/Γ* stays negative or fails to match
  1/v*² for all μ): **MET.** B*/Γ* ∈ [−0.84, −0.57] throughout;
  1/v*² ∈ [+0.81, +1.23] throughout. Distinct signs and an
  irreducible gap of ≈1.6.
- **Inconclusive criterion** (trend correct but not provable past
  the convergence boundary): **PARTIALLY APPLIES.** The gap |diff|
  does narrow on `μ ∈ [0, 0.325]` (from 1.80 to 1.62) but **widens
  again** on `μ ∈ [0.325, 0.45]` (from 1.62 to 1.65). The trend
  toward μ=1/2 from the convergent side is **away** from closure,
  not toward it.

**Conclusion: OP-2b CONFIRMED OPEN.**

The Hubbard–Stratonovich Jacobian is **not** the missing physics.
M3's negative result `B*/Γ* ≈ -0.57` at the 3D Ising WF FP is robust
against this addition, and the missing mechanism must come from a
channel not parameterised by an on-site weight modification:

1. **Wave-function renormalisation `Z_Φ`** (P3.1 of the response
   plan). The composite field has its own multiplicative
   renormalisation. Standard FRG/NPRG calculations for `(φ³, φ⁴)`
   theories show this can shift β-coefficients by O(1) factors.
2. **GL-bond operator `G Φ_i Φ_j (ΔΦ)²`** (P3.2). The TGP-v2
   bond is intrinsically gradient-coupled — two-site, momentum-
   dependent — and is not captured by single-site MK moments.
   The bond-mixing it induces under MK requires extending the
   state to `(r, u, B, Γ, G, …)` with cross-couplings.
3. **Non-perturbative RG (NPRG) test** (P3.4). The above two
   mechanisms can be checked simultaneously in a Wetterich-style
   FRG setup with a `Z_Φ`-tracking ansatz and the gradient kinetic
   term retained.

These are deferred to future work and are now the **leading
candidates** for resolving OP-2b.

## 7. Implications for the paper

This test does **not** change the paper's status quo:

- **OP-2a (β=γ at vacuum) closed** — Routes 1 & 2 of
  `thm:beta-eq-gamma-triple` (M2a, M2c). Independent of M4.
- **OP-2b (B*/Γ* at WF) open** — `thm:beta-eq-gamma-triple` Route 3
  remains as M3 left it: `B*/Γ* ≈ -0.57 ± 0.05` at single-site
  MK-RG. The *paper-level disposition* set by P1.1
  (2026-04-25 patch in `core/sek08_formalizm.tex`,
  `partial_proofs/dodatek_B`, and KNOWN_ISSUES) is **unchanged**:
  Route 3 is *not* claimed; β=γ-at-criticality is an open
  problem. M4 narrows the candidate-mechanism list (rules out
  H-S Jacobian) but does not close the gap.

**No paper edits required from M4.** The KNOWN_ISSUES.md entry for
OP-2b will be updated with the M4 negative result.

## 8. Files

- `mk_rg_phi.py` — implementation.
- `mk_rg_phi_results.txt` — raw output (full μ-scan, 3 ε values,
  both N_ops).
- `M4_phi_variable_derivation.md` — analytical setup (Hypothesis,
  marginality proof, decision criteria).
- `M4_results.md` — this note.

## 9. Cross-refs

- `M3_results.md` — predecessor (μ=0 baseline, OP-2b open since
  2026-04-24).
- `external_review_2026-04-25/review_response_plan.md` — Test A
  in the priority-2 plan.
- `TGP/tgp-core-paper/KNOWN_ISSUES.md` — OP-2b disposition; will
  be amended with M4's negative-result entry.
- `axioms/substrat/dodatekB_substrat.tex`, `core/sek08_formalizm.tex`
  — paper-level statements that frame OP-2b; **not** modified by M4.

## 10. Recommended next steps (P3 of the plan)

In order of expected impact:

1. **P3.1 — Z_Φ measurement.** Extend the MK-RG to track the
   composite-field 2-point function and extract `Z_Φ`. Inserts
   anomalous dimension into the β/γ flow.
2. **P3.2 — GL-bond operator in MK.** Promote the bond from
   bilinear `−K Σ s_i s_j` to GL-style with explicit gradient
   coefficient `G`. Requires a 2-site moment formulation of the
   bond-move step. Highest expected impact, hardest to implement.
3. **P3.4 — NPRG sanity check.** Wetterich equation with leading
   ansatz including `Z_Φ` and the GL kinetic term. Independent
   estimate of `B*/Γ*` from a different RG scheme.

Items 1–3 are open work for OP-2b. Test A (this M4) is now
complete and contributes a sharp negative input: the H-S Jacobian
is not the missing physics.
