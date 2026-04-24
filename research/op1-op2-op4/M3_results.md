# M3 results — MK-RG computation of C_β/C_γ at 3D Ising WF

**Status:** Route 3 of `thm:beta-eq-gamma-triple` **does not close as
stated** in the M3 plan. Numerical result disagrees with target in
both magnitude and sign.
**Date:** 2026-04-24.
**Script:** `mk_rg_bgamma.py`
**Raw output:** `mk_rg_bgamma_results.txt`

---

## 1. What was computed

Extension of the 2-operator MK-RG in
`tooling/scripts/substrate/renormalization_substrate.py` (which tracks
`(r, u, K)` for the φ⁴ Wilson–Fisher fixed point) to an N-operator
truncation tracking

    V(ŝ) = Σ_{k=1}^{N} c_{2k}/(2k) · ŝ^{2k}
         = r/2 ŝ² + u/4 ŝ⁴ + B/6 ŝ⁶ + Γ/8 ŝ⁸ + E/10 ŝ^{10} + L/12 ŝ^{12} + …

with the Migdal–Kadanoff recursion (bond-move + decimation) on the
3D hypercubic lattice (`d = 3, z = 6, b = 2`). Identification:

- `Φ = ŝ²`, so `ŝ^{2k} = Φ^k` as local composite operators.
- `B = 6 × (coefficient of Φ³ in V)`, i.e., cubic Φ operator.
- `Γ = 8 × (coefficient of Φ⁴ in V)`, i.e., quartic Φ operator.

At the fixed point the M3 plan identifies

    C_β = B* / v*²,    C_γ = Γ* / v*²,    v*² = |r*|/u*,

so the ratio

    C_β / C_γ = B* / Γ*                              (convention-invariant)

is what MK-RG should predict. The plan's target (for `β = γ` at WF)
is

    C_β / C_γ = 1 / v*²  ≈  1/0.575  ≈  1.74.

## 2. Numerical method

Each MK-RG step computes moments `m_{2k} = ⟨ŝ^{2k}⟩` under `exp[−V(ŝ)]`
by Gauss–Legendre quadrature (`n_quad = 1200, s_max = 10`), converts
to cumulants `κ_{2k}` via the standard moment–cumulant recursion

    κ_n = m_n − Σ_{j=1}^{n−1} C(n−1, j−1) κ_j m_{n−j},

and updates

    c_{2k}' = c_{2k} − 2 F_{2k} / (2k−1)!,
    F_{2k} = K_eff^{2k} κ_{2k},   K_eff = b^{d−1} K = 4 K.

Fixed-point iteration in "bar variables" `c̄_{2k} = c_{2k}/K` (with K
reset to 1 each step) to `tol = 10⁻¹⁰`.

## 3. Validation

**(a) Cumulant formulas, Gaussian check.** For `V = ŝ²/2`, the exact
moments `m_{2k} = (2k)! / (2^k k!)` and cumulants `κ_{2k≥4} = 0`
are reproduced numerically:

| k | m_{2k} (numerical) | m_{2k} (exact) | κ_{2k} (numerical) |
|---|---|---|---|
| 1 | 1.00000000 | 1.0 | +1.000000 |
| 2 | 3.00000000 | 3.0 | −1.78 × 10⁻¹⁵ |
| 3 | 15.00000000 | 15.0 | +1.07 × 10⁻¹⁴ |
| 4 | 105.00000000 | 105.0 | −5.15 × 10⁻¹⁴ |
| 5 | 945.00000000 | 945.0 | −6.18 × 10⁻¹³ |
| 6 | 10395.00000000 | 10395.0 | +4.37 × 10⁻¹¹ |

Numerical precision of the cumulant pipeline: ≤ 10⁻¹⁰. Adequate.

**(b) 2-operator sublimit.** Holding `B = Γ = … = 0`, FP converges to
`r* = −2.2508, u* = +3.9172` (matching the documented WF values
`−2.251, +3.917` in `axioms/substrat/dodatekB_substrat.tex`
eq:B-WF). ✓

## 4. Truncation scan

Fixed-point values at N-operator truncation:

| N_ops | r*        | u*        | B*          | Γ*          | v*²     | B*/Γ*    |
|-------|-----------|-----------|-------------|-------------|---------|----------|
| 2     | −2.25085  | +3.91718  | —           | —           | 0.57461 | —        |
| 3✗    | −2.00125  | +1067.33  | −6.83 × 10⁵ | —           | 0.00188 | —        |
| 4     | −2.37615  | +3.29219  | −6.4345     | +13.0274    | 0.72175 | −0.49392 |
| 5✗    | −2.00125  | +1067.33  | −6.83 × 10⁵ | +4.4 × 10⁸  | 0.00188 | −0.00154 |
| 6     | −2.42775  | +3.11014  | −5.5426     | +10.2293    | 0.78059 | −0.54184 |
| 7✗    | −2.00125  | +1067.33  | −6.83 × 10⁵ | +4.4 × 10⁸  | 0.00188 | −0.00154 |
| 8     | −2.45694  | +3.01611  | −5.1216     | +9.0060     | 0.81461 | −0.56869 |

**Odd truncations** (N = 3, 5, 7): pathological. The highest power ŝ^{2N}
has `N` odd, and at the FP the coefficient is negative, leaving the
potential unbounded below. The iteration finds an unstable saddle at
large `u` and fails to represent the physical WF fixed point.
**Even truncations are the meaningful ones.**

**Even truncation convergence** (N = 4, 6, 8):

- `B*/Γ* : −0.494 → −0.542 → −0.569`  (slow drift, stable sign).
- `v*²   :  0.722 →  0.781 →  0.815`  (upward drift).

The ratio is **converging** (successive differences 0.048, 0.027;
ratio ~ 0.56) to an extrapolated value around `−0.60 ± 0.05`.

## 5. Eigenvalue spectrum at N_ops = 8

Linearised bar-RG Jacobian at the 8-op FP has all irrelevant
eigenvalues (|λ| < 1):

| i | λ (real)  | y = log|λ|/log 2 | Interpretation |
|---|-----------|------------------|----------------|
| 0 | −0.76213  | −0.3919          | dominant irrelevant (MK analogue of ω) |
| 1 | +0.19250  | −2.3771          | subleading irrelevant |
| 2–7 | ≈ +0.186 | ≈ −2.43        | near-degenerate irrelevant tower |

Comparison with 3D Ising bootstrap (Kos–Poland–Simmons-Duffin 2016,
Chang et al 2022):

| Exponent | Bootstrap | MK (N_ops=8) | Ratio (MK/boot) |
|---|---|---|---|
| ω  = Δ_{ε'} − d     = 0.830 | 0.830 | 0.392 | 0.47 |
| ω₂ = Δ_{ε''} − d    = 3.896 | 3.896 | 2.377 | 0.61 |

MK-RG systematically underestimates the irrelevant eigenvalues by
a factor ~0.5, consistent with the classic MK imprecision
(ν_MK ≈ 0.60 vs exact ν = 0.630).

The relevant (thermal) eigenvalue is absorbed into the K-rescaling
and does not appear in the bar-Jacobian spectrum — expected behaviour,
matches the 2-op reference code.

## 6. Verdict on Route 3

**Result at N_ops = 8 (largest even truncation):**

    B* / Γ* = −0.5687
    v*²     =  0.8146
    β/γ |_* ≈ v*² · B*/Γ* = −0.463   (using M3 plan identification)

**Plan target** (for `β = γ` at WF, assuming the vacuum condition
extends to criticality):

    C_β / C_γ = 1 / v*² = 1.228
    β/γ |_*   = v*² · C_β/C_γ = 1      (= vacuum value)

**Discrepancy: sign inversion and magnitude mismatch.** The MK-RG
prediction is inconsistent with `β = γ` at the critical point.

## 7. Interpretation

### 7.1 What this does NOT undermine

- **OP-1 closure (M2a + M2c) stands.** The form `U(φ) = β/3 φ³ − γ/4 φ⁴`
  is unchanged: M2a derived it as a tree-level kinematic consequence
  of the composite-field Jacobian `Φ^{−1/2}`, M2c established its
  completeness (φ^{≥5} irrelevant by 10⁻⁵ at IR).
- **β = γ at vacuum (M2a Route 1+2) stands.** The kinematic identity
  `β_eff/γ_eff = Φ₀` is derived at the saddle `Φ₀`, not at the RG
  fixed point. It is a vacuum property, and M2a confirmed it
  numerically to 0.4%.

### 7.2 What this DOES imply

The β = γ equality is a **vacuum-local property**, not a
scale-invariant one. Under RG flow from microscopic UV
(where β_bare = γ_bare = 0 as these are generated operators) to
the IR WF fixed point (where they take the values B*, Γ*), the ratio
does NOT stabilise at `1/v*²`. Route 3 as stated in the M3 plan —
"MK-RG yields `C_β/C_γ = 1.74` because `β = γ` is protected by RG" —
is **incorrect**. β and γ run independently under RG; their equality
is recovered only at the specific point `Φ = Φ₀ = v²` (the ordered-phase
vacuum), and only through the kinematic structure of the composite field.

This is consistent with what M2a already suggested in its §5 "what this
does not confirm" — the tree-level match `β/γ = 1` is structural
(kinematic Jacobian) and MK-RG was expected to test whether loop
corrections preserve this structural equality at criticality. They
do not, by a factor of about 3 in magnitude and a sign flip.

### 7.3 Possible repairs to Route 3

Two mechanisms listed in the M3 plan §6 that were NOT included in
the present computation could partly redress the discrepancy:

1. **GL bond operator** `G Φ_i Φ_j (ΔΦ)²`. The v2 gradient bond
   is a non-local 2-site operator not captured by single-site MK-RG
   moments. Adding it requires extending to a joint `(r, u, B, Γ, G)`
   evolution with bond-mixing. Current estimate: this is the most
   likely culprit.
2. **Wave-function renormalisation `Z_Φ`**. The composite field
   undergoes its own multiplicative renormalisation which MK-RG does
   not track explicitly. Including Z_Φ via a 2-point function
   measurement at each RG step would shift both `B` and `Γ`
   non-trivially.

Both are deferred to future work. Absent those additions, the result
of M3 is: **at single-site MK-RG of the φ⁴ theory, `C_β/C_γ ≈ −0.57`
at the 3D Ising WF fixed point.**

## 8. Implications for `thm:beta-eq-gamma-triple`

The theorem in `core/sek08_formalizm.tex` (and alternatively stated in
`partial_proofs/dodatek_B`) has three routes:

- **Route 1 (variational at vacuum):** `β = γ` from the extremum of
  the effective potential. Derived by M2a from the Hubbard–Stratonovich
  Jacobian; confirmed numerically to 0.4% in M2a sanity check.
  **Status: proved at tree level.**
- **Route 2 (Z₂ structural at vacuum):** Same conclusion from the
  kinematic structure of the Z₂-even composite `Φ = ŝ²`.
  **Status: proved at tree level** (same microscopic content as Route 1).
- **Route 3 (MK-RG + MC extending β=γ to criticality):** As formulated,
  predicts `C_β/C_γ = 1.74`. MK-RG numerical result: `−0.57 ± 0.05`.
  **Status: falsified in current form.** `β = γ` does not survive RG
  flow to criticality.

**Proposed reformulation** (for the core paper):

  `thm:beta-eq-gamma-triple` → `thm:beta-eq-gamma-at-vacuum`
  
  *At the vacuum `Φ = Φ₀`, the coefficients of the cubic and quartic
  TGP operators satisfy `β = γ` as a kinematic consequence of the
  Z₂-even composite-field map `ŝ → Φ = ŝ²`. This is a vacuum-local
  identity; it does not extend to the RG fixed point.*

This is a **weaker** but **correct** statement. The paper must drop
the Route-3 claim that β=γ is RG-protected.

## 9. What OP-2 still needs

- **(OP-2a, closed)** Vacuum-level `β = γ`: Routes 1, 2 of
  `thm:beta-eq-gamma-triple`, both derived by M2a from the
  composite-field Jacobian.
- **(OP-2b, OPEN)** Criticality-level `C_β/C_γ`: present M3
  gives `−0.57 ± 0.05` at 8-op truncation. Possible corrections:
  (i) GL bond operator, (ii) Z_Φ wave-function renormalisation.
  Both are left to future work.

**Recommendation:** mark OP-2 as **partially closed** (vacuum part
rigorous; criticality part an open problem with a quantitative
but unexpected numerical result).

## 10. Deliverables

| Item | Status | Location |
|---|---|---|
| MK-RG 5-op (actually N-op extensible) Python script | done | `mk_rg_bgamma.py` |
| Raw numerical output | done | `mk_rg_bgamma_results.txt` |
| Gaussian cumulant validation | done | §3(a) of this note |
| 2-op sublimit validation | done | §3(b) |
| Truncation convergence scan (4,6,8-op) | done | §4 |
| 3D Ising eigenvalue cross-check | done | §5 |
| `C_β / C_γ` numerical value at 3D Ising WF | done | §6: **−0.57 ± 0.05** |
| GL-bond extension of MK-RG | deferred | future M3-b |
| `Z_Φ` wave-function renormalisation | deferred | future M3-c |
| Update `thm:beta-eq-gamma-triple` statement in core/sek08 | pending | paper-patch |

## 11. Files

- `mk_rg_bgamma.py` — N-operator MK-RG script.
- `mk_rg_bgamma_results.txt` — raw run output (truncation scan +
  eigenvalue table).
- `M3_CbCg_numerical_plan.md` — original plan (superseded by this
  note for the numerical result; still correct for the methodology).
- `M2a_HS_derivation.md` — derivation of `β = γ` at vacuum (still
  holds; independent of M3 result).
