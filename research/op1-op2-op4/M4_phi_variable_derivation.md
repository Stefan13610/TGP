# M4 — MK-RG in Φ-variables (analytical setup)

**Date:** 2026-04-25.
**Test designation:** "Test A" of `external_review_2026-04-25/review_response_plan.md`.
**Scope:** analytical derivation supporting the implementation in
`mk_rg_phi.py`. Conclusion in `M4_results.md`.

---

## 1. Hypothesis under test

M3 (`mk_rg_bgamma.py`, `M3_results.md`) tracked the operator basis

```
V(ŝ) = Σ_{k≥1} c_{2k}/(2k) · ŝ^{2k}
```

in the **Z₂-symmetric ŝ-variables** of the substrate, found a 3D-Ising
WF fixed point with

```
B*/Γ* = −0.57 ± 0.05    (target for β=γ at vacuum: +1/v*² ≈ +1.23)
```

The user's conjecture (concurrent with reviewer C5 in the
2026-04-25 external review): the sign-mismatch and magnitude error
arise because **a non-polynomial Z₂-even operator was dropped** from
the M3 basis, namely the Hubbard–Stratonovich Jacobian

```
Φ^{−1/2}     (Φ ≥ 0, the TGP-natural variable)
   ⇕  composite-field map  Φ = ŝ²
(s² + ε²)^{−μ}   (regulated)  in ŝ-variables, with μ_HS = 1/2.
```

In the Φ-action this Jacobian is a `(1/2) · ln Φ` term in `V_Φ`; in
the ŝ-action it is `μ ln(s²) = 2μ ln|s|` with **μ_HS = 1/2**.

**Question:** does adding `μ ln(s² + ε²)` to the on-site V move the
WF fixed point of the polynomial sector, in particular does
`B*(μ)/Γ*(μ)` cross `1/v*²(μ)` somewhere in `μ ∈ [0, 1/2]`?

If yes → β=γ is recovered at the H-S-natural μ; OP-2b closes.
If no → M3's −0.57 is robust against the Jacobian; OP-2b fundamentally
open and one of the deferred mechanisms (Z_Φ, GL bond) must be
responsible.

## 2. Why the H-S Jacobian was missing from M3

M3 starts from the Z₂-symmetric s-action

```
S_s = Σ_i V_s(s_i) − K Σ_{<ij>} s_i s_j,
V_s(s) = Σ_k c_{2k}/(2k) s^{2k}    (no log term).
```

This generates the partition function `Z_s = ∫ Π_i ds_i e^{−S_s}`.

The change of variables Φ_i = s_i² (with **two branches** s = ±√Φ
collapsed by Z₂) gives

```
Z_s = ∫ Π_i dΦ_i · Φ_i^{−1/2} · e^{−S_Φ}    (formal; the Φ ≥ 0 measure
                                              hides a sum over Z₂ branches)
S_Φ = Σ_i V_Φ(Φ_i) − K Σ_{<ij>} ŝ_i ŝ_j,
V_Φ(Φ) = V_s(√Φ).
```

i.e., the **integration measure** in Φ-variables carries the Jacobian
`Φ^{-1/2}`. Equivalently, V_Φ_eff = V_Φ + (1/2) ln Φ.

In the **TGP starting point**, the Φ-action is given *axiomatically*
(GL-substrate bond + on-site potential) — there is no a-priori reason
to include the H-S Jacobian. M3's s-variable formulation is consistent
with **the s-action being primary**; the user's conjecture is
equivalent to asking what happens if we instead take the **Φ-action
as primary** and the s-formulation is derived (with the inverse
Jacobian Φ^{1/2} = |s| in the ŝ-action measure, i.e. `−ln|s|` in V_s,
i.e. **μ = −1/2** in our conventions).

Both signs of μ are physically interesting:
- **μ = +1/2** (H-S Jacobian explicit in s-action): test of "Φ is the
  TGP-natural variable, integrate over it; the s-formulation of MK
  must compensate".
- **μ = −1/2** (s-action primary, Φ derived): probably equivalent to
  M3 by definition.

## 3. Marginality of μ under single-site MK-RG

**Claim:** Under bond-move + decimation + bar-rescaling with the
bilinear bond `−K Σ s_i s_j`, the operator `μ ln(s² + ε²)` is
**exactly marginal** (μ' = μ at every step, no flow).

**Proof.**

(i) **Decimation.** Integrating out a site `s_2` with neighbours
`s_1, s_3` and the modified weight:

```
e^{−V'(s_1, s_3)} = ∫ ds_2 (s_2² + ε²)^{−μ} e^{−V_poly(s_2) + K_eff s_2 (s_1+s_3)}
                 ≡ Z(h),    h ≡ K_eff (s_1 + s_3).
```

`Z(h)` is an analytic function of `h` (for ε > 0 the integrand is
smooth; for ε = 0 with μ < 1/2 the singularity at s_2 = 0 is
integrable). Z₂ symmetry of the integrand under s_2 → −s_2 (both
`(s² + ε²)^{−μ}` and `e^{−V_poly}` are even) makes `Z(h) = Z(−h)`,
hence

```
log Z(h) = log Z(0) + Σ_{n≥2 even} κ_n h^n / n!,
```

where `κ_n` are cumulants of the **modified** weight (μ-dependent).
**No `ln(h)` term is generated** — the decimated effective potential
is polynomial in `(s_1 + s_3)`. Consequently, after the standard MK
manipulation `V'(s) = V(s) − 2[F(s) − F_0]`, the log coefficient is
unchanged: μ' = μ.

(ii) **Bond-move.** `K → K · b^{d−1} = 4K`. Acts only on bonds, not
on the on-site potential. μ' = μ.

(iii) **Bar-rescaling** (M3 convention `c̄_{2k} = c_{2k}/K^k` after
each step). The log term `μ ln(s² + ε²)` under field rescaling
`s → s · K^{−1/2}` becomes `μ ln(s²/K + ε²/K) = μ ln(s² + ε²) − μ ln K`
(treating ε as rescaled equivalently). The additive constant `−μ ln K`
is absorbed into the field-independent `F_0`, leaving μ unchanged.
μ' = μ.

**Conclusion:** μ labels a **continuous family of MK-RG fixed
points**. Each μ gives a (potentially) different polynomial WF FP
location `(c_2*(μ), c_4*(μ), c_6*(μ), …)`. The empirical question:
how does `B*/Γ*` depend on μ, and does it cross `1/v*²` at some
μ ∈ (0, 1/2]?

## 4. Effect of μ on polynomial cumulants

Although μ doesn't flow, it **changes the cumulants** that drive
the polynomial flow. Concretely:

```
m_n(c, μ) = ⟨s^n⟩ = ∫ds s^n (s²+ε²)^{−μ} e^{−V_poly(s,c)} / Z(c, μ),
κ_n(c, μ) = (moment-to-cumulant recursion).
```

For the same polynomial couplings `c`, larger μ shifts the moment
weight toward the tail (because the small-|s| region is suppressed
by `(s²)^{−μ}` for large s² · ε^{-2}? no — wait: `(s²+ε²)^{−μ}` is
**larger** for small `s² + ε²`, which **enhances** the small-|s|
region by a factor `~|s|^{−2μ}` for `|s| > ε`). So increasing μ:
- **enhances** small-|s| weight (where V_poly is small)
- mildly **suppresses** large-|s| tails
- net effect: shifts `⟨s^{2k}⟩` distributions toward smaller values

This in turn rescales the κ_{2k} pattern, hence the FP location.

## 5. Convergence and regulator

For Z₂-symmetric quadrature the integrand `(s²)^{−μ} e^{−V_poly}`
near `s=0` behaves as `|s|^{−2μ}`, integrable iff `2μ < 1`, i.e.
**μ < 1/2** in the unregulated case (ε = 0).

The H-S-natural value `μ = 1/2` is therefore the **boundary** of
convergence: a logarithmic divergence at `s=0` in the symmetric s-form,
even though it is perfectly finite in the asymmetric `Φ ≥ 0` form
(`∫_0^∞ dΦ Φ^{-1/2} e^{−V_Φ}` converges).

For numerical work we use the regulator ε > 0 and scan μ ∈
{0, 0.1, 0.2, 0.3, 0.4} with ε ∈ {0.1, 0.01}. We extrapolate the
trend toward μ = 1/2 to assess whether `B*/Γ*` is approaching
`1/v*²` (the β=γ target).

The regulator ε corresponds physically to a UV cutoff on `Φ` from
below: `Φ ≥ ε²`, i.e. the substrate has a minimum non-empty
vacancy density. This is consistent with the TGP picture (`Φ ≥ 0`,
with `Φ = 0` only at N₀ defects) and provides a natural physical
regulator if the test gives a μ-trend.

## 6. Decision criteria for OP-2b

Let `R(μ) ≡ B*(μ)/Γ*(μ)` and `T(μ) ≡ 1/v*²(μ)` (the β=γ target).

1. **OP-2b restored (positive result):** `R(μ_*) = T(μ_*)` for some
   μ_* ∈ (0, 0.5], with sign(R) flipping from − (M3) to + as μ
   increases. β=γ at WF would then be a property of the H-S-natural
   FP at μ ≈ 1/2.

2. **OP-2b confirmed open (negative result):** `R(μ)` stays negative
   or does not match `T(μ)` for all `μ ∈ [0, 0.5)`. M3's result is
   robust; the missing physics must come from another channel:
   wave-function renormalisation `Z_Φ` (P3.1), or GL bond operator
   (P3.2), or NPRG (P3.4).

3. **Inconclusive (mixed result):** `R(μ)` trends correctly but does
   not cross `T(μ)` in the convergent regime; extrapolation to
   μ = 1/2 is suggestive but not provable without going past the
   logarithmic divergence.

In all cases this test **rules in or out** the H-S Jacobian as the
mechanism behind M3's deviation; the result feeds into the disposition
of OP-2b in `KNOWN_ISSUES.md`.

## 7. Files

- `mk_rg_phi.py` — implementation (extension of `mk_rg_bgamma.py`).
- `mk_rg_phi_results.txt` — raw output.
- `M4_results.md` — final interpretation and OP-2b verdict.
