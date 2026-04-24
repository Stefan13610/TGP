# M6 — GL-bond operator in MK-RG (analytical setup)

**Date:** 2026-04-25.
**Test designation:** Test C (P3.2) of
`external_review_2026-04-25/review_response_plan.md`.
**Scope:** analytical setup supporting `mk_rg_glbond.py`.
Verdict in `M6_results.md`.

---

## 1. Hypothesis under test

After M4 (H-S Jacobian, negative) and M5 (Z_Φ via η-deformation,
negative), the only remaining single-channel candidate from the
M3 §6 list is the **GL-bond operator** — the v2 axiom-level bond:

```
H_Γ = ... + J_GL Σ_{<ij>} A_ij Φ_i Φ_j (Φ_j − Φ_i)²            (v2)
```

In `ŝ`-variables (`Φ = ŝ²`):

```
J_GL ŝ_i² ŝ_j² (ŝ_j² − ŝ_i²)² = J_GL [ŝ_i² ŝ_j^6 − 2 ŝ_i^4 ŝ_j^4 + ŝ_i^6 ŝ_j²].
```

This is a degree-8, two-site, momentum-dependent operator. It is
**not** captured by single-site MK-RG moments — it is intrinsically
**non-local** in a way that the H-S Jacobian (on-site weight) and
Z_Φ (field-strength rescaling) are not.

**Question:** does the GL bond, when added to M3's bilinear bond as
a perturbation, shift the WF FP enough to close OP-2b?

If yes → the v2 axiom-level bond is the missing physics.
If no → all three single-channel candidates are exhausted; OP-2b
is genuinely open at the level of single-operator additions to
single-site MK.

## 2. Approach: J_GL as on-site-weight perturbation

The cleanest tractable test in the M3 framework is to add the GL
bond as a perturbation to the bilinear bond, keeping M3's
decimation-as-source structure intact. The trick is that the GL
bond, when one site is integrated out, contributes additional
polynomial terms to the on-site `V'(s)` at the surviving site.

In M3's decimation step, integrating out `s_2` with neighbours
`s_1, s_3` gives:

```
F(h) = ln ∫ ds_2 exp[−V(s_2) + K_eff h s_2],    h = K_eff (s_1+s_3)
```

(Note: bond-move K → K_eff = b^{d-1} K = 4K is already absorbed
into the weight `K_eff h s_2`.)

With the GL bond added, the integrand picks up
`exp[−J_GL · (bond_12 + bond_23)]` where

```
bond_12 = s_1² s_2² (s_2² − s_1²)² = s_1² s_2^6 − 2 s_1^4 s_2^4 + s_1^6 s_2².
```

and similarly for `bond_23`. **Crucially, this is even in s_2 and
polynomial of degree 8**, so it contributes to the moments
`m_n = ⟨s_2^n⟩` of the integrand at order `J_GL`. Specifically, to
first order in `J_GL`:

```
F_GL(s_1, s_3) ≈ F_M3(h) − J_GL [⟨bond_12⟩ + ⟨bond_23⟩] + O(J_GL²)
```

where `⟨·⟩` is the M3 single-site average (under the
weight `exp[-V(s_2) + K_eff h s_2]`, which depends on h, hence on
s_1+s_3 — but the bond expectations also depend explicitly on
s_1 and s_3 through their polynomial coefficients).

**Simplification.** Since `bond_12 = s_1² s_2^6 − 2 s_1^4 s_2^4 +
s_1^6 s_2²` is linear in three powers of s_1 with separated
coefficients (which are 6th, 4th, 2nd moments of s_2), and bond_23
has the same structure with `s_3` replacing `s_1`, the J_GL
correction to the effective potential is:

```
ΔV_GL'(s_1, s_3) = +J_GL [s_1² ⟨s_2^6⟩(h) − 2 s_1^4 ⟨s_2^4⟩(h) + s_1^6 ⟨s_2²⟩(h)
                        + s_3² ⟨s_2^6⟩(h) − 2 s_3^4 ⟨s_2^4⟩(h) + s_3^6 ⟨s_2²⟩(h)]
```

where `⟨s_2^n⟩(h)` is the n-th moment of `s_2` in M3's source-h
measure.

For the M3-style decomposition `V'(s) = V(s) − 2[F(s) − F_0]` to
work, we extract the **on-site contribution** by setting `s_1 =
s, s_3 = 0` (or symmetrically), and take half:

```
ΔV_on-site_GL(s) ≈ +J_GL [s² ⟨s²^6⟩(K_eff·s) − 2 s^4 ⟨s²^4⟩(K_eff·s)
                        + s^6 ⟨s²^2⟩(K_eff·s)]
```

This is now a function of `s` alone, polynomial in `s` (after
expanding `⟨s_2^n⟩(h)` as a power series in h = K_eff·s), and can be
distributed across the c_{2k} updates.

## 3. Two-track implementation

We pursue **two complementary tests**:

### Track A: J_GL as fixed parameter, scan

Treat `J_GL` as a fixed parameter throughout the flow (analogous
to μ in M4, η in M5). At each MK step, add `ΔV_on-site_GL` to V'.
Find FP, report `B*/Γ*` and `1/v*²` as a function of J_GL.
Closure criterion: B*/Γ* − 1/v*² changes sign at some
physically reasonable `|J_GL|`.

### Track B: J_GL flowing alongside (r, u, B, Γ, ...)

Track J_GL as an additional coupling. The flow equation is
schematically:

```
J_GL' = b^{d-1} J_GL · (rescaling factor from operator dimensions)
       + O(J_GL · couplings)  (operator mixing)
```

The leading-order (canonical) scaling: bond_12 contains ŝ²ŝ^6 +
ŝ^4ŝ^4 + ŝ^6ŝ², total degree 8. Under field rescaling
`s → b^{(d-2)/2} s = b^{1/2} s` (d=3, eta=0), each ŝ² picks up b,
so `bond_12 → b^4 · bond_12`. With bond move `b^{d-1} = 4` for
the operator coefficient: J_GL → J_GL · b^{d-1} · b^{-(d_op - d)}
where `d_op` is the operator's scaling dimension... let's just
measure it numerically.

**For this first pass, we focus on Track A (the J_GL scan). Track
B is a fall-through extension.**

## 4. Implementation outline

`mk_rg_glbond.py` extends `mk_rg_bgamma.py`'s `MigdalKadanoffRGN`:

```python
class MigdalKadanoffGL(MigdalKadanoffRGN):
    def __init__(self, n_ops, J_GL=0.0, ...):
        super().__init__(...)
        self.J_GL = float(J_GL)

    def rg_step(self, couplings, K):
        # Standard M3 update first.
        new_c, K_new = super().rg_step(couplings, K)
        # GL correction (Track A): add J_GL contribution to c_{2,4,6,8,...}.
        if self.J_GL != 0.0:
            # See section 5 below for the explicit correction formula.
            ...
        return new_c, K_new
```

The correction depends on which c_{2k} the GL bond modifies.

## 5. Explicit GL contribution to c_{2k} update

Working in the decimation integral, set `s_3 = 0` (extract on-site
piece). The integrand is:

```
exp[−V(s_2) + K_eff·s·s_2] · exp[−J_GL · bond_12(s, s_2)] + (s ↔ 0 piece)
```

To first order in J_GL:

```
F(s, 0) ≈ F_M3(K_eff·s) − J_GL · ⟨bond_12(s, s_2)⟩
       = F_M3(K_eff·s) − J_GL · [s² ⟨s_2^6⟩(K_eff·s)
                                 − 2 s^4 ⟨s_2^4⟩(K_eff·s)
                                 + s^6 ⟨s_2^2⟩(K_eff·s)]
```

Where `⟨s_2^n⟩(h)` is the cumulant-generating-function-derived
n-th moment. By the same generating-function logic as M3:

```
⟨s_2^n⟩(h) = Σ_k (binomial coefficients) · κ_{2j} · h^{n-2j}
```

— i.e., a polynomial in h of degree n. Plugging h = K_eff·s and
expanding, we get a polynomial in s. The s^{2k} coefficient of
`F(s, 0)` then receives a J_GL contribution that we compute as:

```
Δc_{2k}(J_GL) = (some combinatorial factor) × κ_{2j} × K_eff^{2j}
              × J_GL × (powers of s combined)
```

A clean implementation: numerically compute `⟨bond_12(s, s_2)⟩` at
each Gauss-Legendre node `s` using the same single-site quadrature
as M3 (no extra integral needed), then read off the polynomial
coefficients. This avoids writing out the cumulant expansion
analytically.

## 6. Decision criteria

Define `R(J_GL) ≡ B*(J_GL)/Γ*(J_GL)` and `T(J_GL) ≡ 1/v*²(J_GL)`.

1. **Closure (positive):** `R(J_GL_*) = T(J_GL_*)` for some J_GL_*
   in the convergent regime (likely |J_GL_*| ≲ 1 at the M3 FP
   scale). The GL bond closes OP-2b. **Falsifiability:** the
   bootstrap-derived value of J_GL (if computable) should match
   J_GL_*.
2. **Closure-in-principle:** crossing exists but at unphysically
   large J_GL. GL bond contributes but is not dominant (need
   higher-order J_GL² or NPRG).
3. **No closure:** no crossing in the convergent regime. **All
   three single-channel candidates exhausted**; OP-2b is genuinely
   open and requires either NPRG (P3.4) or a fundamentally new
   mechanism.

## 7. Caveats

- **First-order in J_GL.** Track A truncates at O(J_GL). For large
  |J_GL|, the J_GL² and higher terms will matter. We restrict to
  small/moderate |J_GL| where this is justified.
- **On-site projection.** Setting s_3 = 0 to extract the on-site
  contribution loses the genuine "bond" piece of the GL operator.
  This is a partial test: it captures how the GL bond renormalises
  the on-site V (which is what M3's c_{2k} track), but not the
  separate bond flow. Track B (full bond flow) would address this.
- **Bond move.** The GL bond after bond move would normally pick
  up `b^{d-1}` parallel-bond combination factor. Since we're
  treating the GL contribution at first order in J_GL, this
  factor is implicitly absorbed by using K_eff = b^{d-1} K in the
  decimation integral.

## 8. Files

- `mk_rg_glbond.py` — implementation (Track A: J_GL fixed, scan).
- `mk_rg_glbond_results.txt` — raw output.
- `M6_results.md` — verdict on GL bond closing OP-2b.
