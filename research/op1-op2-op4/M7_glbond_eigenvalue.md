# M7 — GL-bond eigenvalue at the M3 fixed point (analytical setup)

**Date:** 2026-04-25.
**Test designation:** Track-B sketch of P3.2 of
`external_review_2026-04-25/review_response_plan.md`.
**Scope:** analytical setup supporting `mk_rg_glbond_eig.py`.
Verdict in `M7_results.md`.

---

## 1. Why M7 is needed

M6 (Track A) found a closure of B*/Γ* = 1/v*² at J_GL ≈ +5.89, but
the value lies **outside** the perturbative regime (|J_GL| ≲ 1)
where the first-order-in-J_GL expansion of M6 is trustworthy. The
central caveat from M6 §3:

> "Track A's first-order J_GL truncation is grossly violated near
>  the crossing; J_GL² and higher terms are O(1) corrections to the
>  answer."

The deeper question that M6 leaves open: under proper MK flow with
J_GL flowing alongside (r, u, B, Γ, …), does the GL coupling **grow**
under coarse-graining (relevant) or **shrink** to zero (irrelevant)?

If λ_GL > 1 (relevant): a tiny bare J_GL grows to finite IR value,
plausibly driving closure even from a small starting point.
If |λ_GL| < 1 (irrelevant): J_GL decays to zero in the IR, and the
v2 bond cannot — at the level of single-site MK linearised at the
M3 FP — close OP-2b regardless of bare value.

M7 computes λ_GL directly.

## 2. Definition of λ_GL

Consider the extended state (cb_2, cb_4, …, cb_{2N}, J_GL_bar) where
cb_{2k} = c_{2k}/K (M3's bar variables) and J_GL_bar = J_GL/K
(natural choice if J_GL transforms identically to the on-site
couplings under K-rescaling; alternative conventions are noted below).

Linearise the MK flow at the M3 FP (cb_M3*, J_GL_bar = 0). The
J_GL eigenvalue is

```
λ_GL = ∂J_GL_bar_post / ∂J_GL_bar_pre |_{cb_M3*, J_GL_pre = 0}.
```

Equivalently, the diagonal J_GL component of the (N+1)×(N+1)
linearised Jacobian at the extended FP — assuming the off-diagonal
mixing with the c-sector is small (justified at first order in J_GL,
where the c-sector flow is unchanged from M3 plus an O(J_GL)
correction; this O(J_GL) correction shifts the FP location but not
the J_GL eigenvalue at leading order).

Verdict criteria:

| |λ_GL| | Status of GL bond at M3 FP |
|---|---|
| > 1 | RELEVANT — grows under MK; can drive closure from small bare value |
| ≈ 1 | MARGINAL — slow flow; depends on subleading corrections |
| < 1 | IRRELEVANT — decays to zero under MK; cannot close OP-2b on its own |

## 3. Computational recipe

### 3.1 The bond Hamiltonian decomposition

After bond move + decimation through one eliminated site `s_2`
between two surviving neighbours `s_1, s_3`, the partition contribution
is

```
exp[F(s_1, s_3)]
  = ∫ds_2 exp[-V(s_2) + K_eff·s_2·(s_1+s_3) - J_GL_eff·(bond_12 + bond_23)]
```

with K_eff = b^{d-1}·K = 4K, J_GL_eff = b^{d-1}·J_GL = 4 J_GL,
bond_12 = s_1² s_2^6 - 2 s_1^4 s_2^4 + s_1^6 s_2², and similarly
for bond_23.

The new effective Hamiltonian satisfies S_eff = -F + (on-site terms).
Splitting F into on-site + cross pieces:

```
F(s_1, s_3) = a(s_1) + a(s_3) + b(s_1, s_3) + const
```

with `b(s_1, 0) = b(0, s_3) = 0`. Then:
- V_new(s) = V(s) - 2·a(s)  (each surviving site neighbours 2 eliminated sites)
- H_bond_new(s_1, s_3) = -b(s_1, s_3) (the new bond Hamiltonian)

Identifying the new K and J_GL:

```
H_bond_new(s_1, s_3) = -K_new·s_1·s_3
                       + J_GL_new·[s_1²s_3^6 - 2 s_1^4 s_3^4 + s_1^6 s_3²]
                       + (other higher-bond operators)
```

so

```
K_new     = +[coefficient of s_1·s_3  in F_cross(s_1, s_3)]
J_GL_new  = -[coefficient of s_1²s_3^6 in F_cross(s_1, s_3)]
```

(F_cross differs from F by s_1-only and s_3-only terms; for monomials
with both a, b > 0 the coefficients agree.)

### 3.2 First-order-in-J_GL evaluation

To first order in J_GL_eff:

```
F(s_1, s_3) = F_M3(s_1+s_3) - J_GL_eff·⟨bond_12 + bond_23⟩(s_1, s_3)
            + O(J_GL_eff²),
```

where ⟨·⟩ is the average under the source-shifted M3 single-site
measure exp[-V(s_2) + h·s_2] with h = K_eff·(s_1+s_3). Defining
M_n(h) = ⟨s_2^n⟩(h):

```
⟨bond_12 + bond_23⟩(s_1, s_3)
   = (s_1² + s_3²)·M_6(h)
     - 2(s_1^4 + s_3^4)·M_4(h)
     + (s_1^6 + s_3^6)·M_2(h),
   h = K_eff·(s_1+s_3).
```

Crucially, M_n(h) is **even** in h for even n (since V is even in
s_2), so the small-h Taylor expansion only involves even powers
of h.

### 3.3 Extraction of J_GL_post (linear response)

Define
```
ΔF/J_GL_pre(s_1, s_3)
   = -4·[(s_1²+s_3²) M_6(h) - 2(s_1^4+s_3^4) M_4(h)
                                + (s_1^6+s_3^6) M_2(h)],
```
(the factor −4 absorbs J_GL_eff = 4 J_GL_pre).

`J_GL_post` at first order is then

```
J_GL_post = -[coeff of s_1²s_3^6 in ΔF] · J_GL_pre + O(J_GL_pre²)
                     ≡ λ_pre · J_GL_pre + O(J_GL_pre²).
```

Bilateral symmetry (s_1 ↔ s_3) gives a consistency check:

```
[coeff of s_1²s_3^6 in ΔF]  =  [coeff of s_1^6 s_3² in ΔF].
```

### 3.4 Bar-rescaling

In M3's bar convention K_bar = 1 each step. Therefore

```
K_pre_bar = 1,        K_post = K_new (numerically computed at M3 FP),
J_GL_pre_bar = J_GL_pre,
J_GL_post_bar = J_GL_post / K_new.
```

The eigenvalue is

```
λ_GL = ∂J_GL_post_bar / ∂J_GL_pre_bar = λ_pre / K_new.
```

### 3.5 Operator-mixing caveat

A subtlety: ⟨bond_12 + bond_23⟩, when expanded in monomials
{s_1^a s_3^b}, does **not** have the pure GL-bond shape
(1 : −2 : 1) in (a, b) ∈ {(2,6), (4,4), (6,2)}. The analytic
structure (see §4) gives

```
[coeff at (2,6)]  =  μ_6^{(6)}·K_eff^6/45 - μ_4^{(4)}·K_eff^4/2
                     + μ_2^{(2)}·K_eff²/2,
[coeff at (4,4)]  =  μ_6^{(6)}·K_eff^6/24 - μ_4^{(4)}·K_eff^4/6,
[coeff at (6,2)]  =  [coeff at (2,6)]            (bilateral symmetry),
```

where μ_n^{(k)} = ∂^k M_n/∂h^k|_{h=0}. The (2,6) coefficient
defines J_GL_post via the convention `J_GL ↔ s_i² s_j^6`. The
(4,4) coefficient of the linear response then reports a **shape**
ratio c_(4,4)/c_(2,6) which equals −2 only if the linear response
is purely GL-shaped.

If c_(4,4)/c_(2,6) ≠ −2, the GL bond, under one MK step at the M3 FP,
generates **other** degree-8 bond operators (such as `s_i^4 s_j^4`
"pure-quartic-cross" or the symmetrised `(s_i s_j)^4` operator,
distinct from O_GL). For the purposes of computing λ_GL we use the
(2,6) extraction: this is the linear response of the J_GL coupling
defined by `J_GL · s_i² s_j^6 + (sym)`. The shape mismatch is
reported as an operator-mixing diagnostic.

We also report a **norm-projection** estimate
λ_GL_norm = ⟨ΔF/J_GL_pre, O_GL⟩ / ⟨O_GL, O_GL⟩ · 1/K_new,
which is the projection onto the GL direction in the L² metric on
[−s_max, s_max]² (uniform Gauss-Legendre). Agreement of the two
estimates within a few % indicates the GL extraction is robust.

## 4. Analytic eigenvalue (sanity check)

Direct expansion of M_n(h) = Σ_k μ_n^{(k)}/k! · h^k (only even k for
even V) and substitution into
⟨bond_12 + bond_23⟩|_{(2,6)-monomial} gives, at the M3 FP and K=1:

```
λ_pre = -[coeff of s_1²s_3^6 in ΔF/J_GL_pre]
      =  4·[μ_6^{(6)} K_eff^6/45 - μ_4^{(4)} K_eff^4/2 + μ_2^{(2)} K_eff^2/2]
      with K_eff = 4 (so K_eff^2 = 16, K_eff^4 = 256, K_eff^6 = 4096)
      = (16384/45) μ_6^{(6)} - 512 μ_4^{(4)} + 32 μ_2^{(2)}.
```

The Taylor coefficients μ_2^{(2)}, μ_4^{(4)}, μ_6^{(6)} are computed
from the M3 FP cumulants:

```
μ_2^{(2)} = κ_4 + 2 κ_2²
μ_4^{(4)} = (lengthier; computed numerically by polynomial fit)
μ_6^{(6)} = (idem)
```

In practice we **compute λ_pre numerically** from the 2D F(s_1, s_3)
projection — the analytic formula serves only as a cross-check.

## 5. Canonical-scaling baseline

Tree-level (Gaussian) prediction for λ_GL:

- Bond move: J_GL → b^{d-1}·J_GL = 4·J_GL.
- s-rescaling under bar convention (K_bar = 1): each s rescales by
  √K_new (set by FP condition). The bond_GL has 8 powers of s, so
  bond_GL → K_new^4·bond_GL.
- For J_GL · bond_GL to remain in the action: J_GL_bar →
  J_GL_bar / K_new^4 (canonical part).
- Net canonical: λ_GL_can = b^{d-1} / K_new^4 = 4 / K_new^4.

At the 3D Ising MK FP, K_new ≈ ?? (to be computed in the script).
For order-of-magnitude: if K_new ~ O(1) at the FP, then
λ_GL_can ~ O(1) — could be relevant or marginal. **The canonical
prediction depends on K_new and is computed in the script for
direct comparison with the numerical eigenvalue.**

## 6. Decision criteria

| Outcome | Verdict |
|---|---|
| |λ_GL| > 1 | GL bond is **relevant**: small bare J_GL grows under MK. Track-A's closure-in-principle finding becomes a serious candidate (re-test at extended FP). |
| 0.5 < |λ_GL| < 1 | GL bond is **irrelevant but slow**: bare J_GL decays, but slowly; finite-J_GL effects could persist over RG depth. |
| |λ_GL| < 0.5 | GL bond is **strongly irrelevant**: bare J_GL decays fast; M3 + GL ≈ M3 in IR. Track-A closure was an artefact of perturbatively-large J_GL. |
| Negative real or complex λ_GL with |·| > 1 | Anomalous; warrants double-check (could indicate operator-mixing into a relevant direction, requiring full Jacobian). |

## 7. Implementation outline

`mk_rg_glbond_eig.py`:

1. Find M3 FP cb_M3* (8-op truncation, M3 baseline).
2. Compute K_new at cb_M3* via one M3 RG step.
3. Build outer 1D Gauss-Legendre grid (n_outer ≈ 30) on [−s_max, s_max].
4. For each outer pair (s_1[i], s_3[j]), compute h = K_eff(s_1+s_3),
   and via vectorised inner Gauss-Legendre quadrature (n_quad = 1200)
   compute M_2(h), M_4(h), M_6(h).
5. Form the 2D linear-response array
   `dF[i,j] = -4·[(s1²+s3²)·M_6 - 2(s1^4+s3^4)·M_4 + (s1^6+s3^6)·M_2]`.
6. Project dF onto monomials s_1^a s_3^b with a, b ≥ 0 even,
   a + b ≤ 12 (or ≤ 16 for full operator basis), via L²-weighted lstsq.
   - Read off λ_pre = -[coeff(2,6)] and bilateral-symmetry check from
     [coeff(6,2)] (should agree).
   - Read off shape diagnostic c_(4,4)/c_(2,6) (expected −2 if pure GL).
7. λ_GL = λ_pre / K_new.
8. Output:
   - Numerical λ_GL.
   - Bilateral consistency: |coeff(2,6) − coeff(6,2)| / |coeff(2,6)|.
   - Shape ratio c_(4,4) / c_(2,6) (compared to GL = −2).
   - Norm-projection λ_GL_norm.
   - Canonical baseline 4/K_new^4.
   - Verdict against §6 criteria.

## 8. Files

- `mk_rg_glbond_eig.py` — implementation (this M7).
- `mk_rg_glbond_eig_results.txt` — raw output.
- `M7_results.md` — verdict.
