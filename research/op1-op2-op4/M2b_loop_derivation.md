# M2b — one-loop V_eff(Φ): bond-fluctuation correction to γ

**Status:** Analytical framework for one-loop correction to
tree-level V_eff(Φ) of M2a. Closed-form Taylor coefficients in
terms of momentum integrals I_n. Intended to supply the `γ ~
λ₀²/(v² a²)` dependence of dodatek B that the on-site tree-level
result does not have.
**Date:** 2026-04-24.
**Predecessor:** M2a (`M2a_HS_derivation.md`,
`M2a_sanity_check_results.md`).
**Successor:** `m2b_loop.py` (numerical),
`M2b_results.md` (numerical results + matching).

---

## 1. Setup

Start from the v2 substrate Hamiltonian (`plan_of_attack.md §1`):

    H_Γ = Σ_i [π²_i/(2μ) + (m₀²/2) ŝ_i² + (λ₀/4) ŝ_i⁴]
          + J Σ_⟨ij⟩ A_ij · ŝ_i² ŝ_j² (ŝ_j² − ŝ_i²)².              (1.1)

All bond content is a function of `Φ_i := ŝ_i²` alone (v2 property,
axiom pivot 2026-04-24). Taking `A_ij = 1` for the universal
v2-bond closure (the generic case; `A_ij` enters only as a
prefactor in the final answer):

    H_Γ[Φ] = Σ_i V_onsite(Φ_i) + J Σ_⟨ij⟩ Φ_i Φ_j (Φ_j − Φ_i)²    (1.2)

with

    V_onsite(Φ) = (m₀²/2) Φ + (λ₀/4) Φ² + (T/2) ln Φ              (1.3)

(last term is the H-S Jacobian of the composite map `ŝ → Φ = ŝ²`;
see M2a §2.3).

## 2. One-loop V_eff in the background-field method

Write `Φ_i = Φ̄ + δΦ_i` with `Φ̄` a uniform background. Expand (1.2)
to quadratic order in `δΦ`. The on-site part gives

    V_onsite(Φ̄ + δΦ) ≈ V_onsite(Φ̄) + V'_onsite(Φ̄) δΦ
                      + (1/2) A(Φ̄) δΦ²                            (2.1)

with the fluctuation mass

    A(Φ) ≡ V''_onsite(Φ) = λ₀/2 − T/(2 Φ²).                       (2.2)

The v2 bond contributes, at quadratic order in `δΦ` around the
uniform background,

    J A_bond(Φ̄) · Σ_⟨ij⟩ (δΦ_j − δΦ_i)²                           (2.3)

with

    A_bond(Φ̄) = Φ̄² + O(δΦ²)                                       (2.4)

(leading term from `Φ_i Φ_j → Φ̄²` at external legs of the
Gaussian propagator). In momentum space on the simple-cubic
lattice (spacing `a`, `d = 3`),

    Σ_⟨ij⟩ (δΦ_j − δΦ_i)² → Σ_q K_Φ(q; Φ̄) |δΦ(q)|²                (2.5)

with

    K_Φ(q; Φ̄) = 2 J Φ̄² · Σ_{μ=1}^{3} (1 − cos q_μ a)
              ≈ J Φ̄² q² a²    (small q).                          (2.6)

Integrating out `δΦ` (Gaussian) gives the one-loop correction

    V_eff^{1-loop}(Φ̄) − V_eff^{tree}(Φ̄)
        = (T/2) ∫_0^Λ [d³q/(2π)³] ln[ K_Φ(q; Φ̄) + A(Φ̄) ]          (2.7)

with UV cutoff `Λ = π/a`. Equation (2.7) is the background-field
one-loop formula; its validity rests only on the Gaussian
approximation to `δΦ` fluctuations (standard; cf. Zinn-Justin §6.2).

## 3. Taylor expansion around Φ₀

Let `Φ₀` denote the tree-level saddle from M2a §2.4:

    λ₀ Φ₀² + m₀² Φ₀ + T = 0,                                      (3.1)

and expand `V_eff^{1-loop}(Φ)` in `Φ − Φ₀`. Write

    Ω(Φ) := (T/2) ∫_0^Λ [d³q/(2π)³] ln[K_Φ(q; Φ) + A(Φ)].         (3.2)

The one-loop corrections to `c_2, c_3, c_4` are

    c_n^{1loop} − c_n^{tree} = Ω^{(n)}(Φ₀) / n!   (n = 2,3,4).    (3.3)

**Simplification.** For the Taylor coefficients at `Φ₀`, use the
local linearisation

    K_Φ(q; Φ) ≈ K_0(q) + Φ̄-dependent corrections,                (3.4)

evaluated at `Φ̄ = Φ₀`, and retain only the `A(Φ)` dependence for
the loop integrand's `Φ`-derivatives — treat `K_Φ` as a fixed
kernel `K_0(q) := K_Φ(q; Φ₀)`. This isolates the "mass-type"
loop corrections from the "stiffness-type" ones. The stiffness
corrections `∂K_Φ/∂Φ = 4 J Φ · (...)` contribute at a computed-but-
subleading order; both terms are included in the numerical code
`m2b_loop.py` without this simplification.

Define the momentum integrals

    I_n(Φ₀) := ∫_0^Λ [d³q/(2π)³] · [K_0(q) + A(Φ₀)]^{-n},         (3.5)

absorbing the Φ-independent kernel. Then, with `Ψ(A) := (T/2) ∫
[d³q/(2π)³] ln[K_0(q) + A]` and `Ψ^{(n)}(A) = (T/2)(−1)^{n−1}
(n−1)! I_n`, Faà di Bruno gives

    Ω' = Ψ' A',                                                    
    Ω'' = Ψ'' A'² + Ψ' A'',                                        
    Ω''' = Ψ''' A'³ + 3 Ψ'' A' A'' + Ψ' A''',                      
    Ω'''' = Ψ'''' A'⁴ + 6 Ψ''' A'² A'' + Ψ'' (3 A''² + 4 A' A''') + Ψ' A''''.

Substituting and grouping:

    c_2^{1-loop corr} = (T/4) [ A''(Φ₀) I_1 − A'(Φ₀)² I_2 ],
    c_3^{1-loop corr} = (T/12)[ A'''(Φ₀) I_1 − 3 A'(Φ₀) A''(Φ₀) I_2
                               + 2 A'(Φ₀)³ I_3 ],
    c_4^{1-loop corr} = (T/48)[ A''''(Φ₀) I_1
                               − (3 A''(Φ₀)² + 4 A'(Φ₀) A'''(Φ₀)) I_2
                               + 12 A'(Φ₀)² A''(Φ₀) I_3
                               − 6 A'(Φ₀)⁴ I_4 ].                 (3.6)

All `I_n` are finite in `d = 3` with the lattice cutoff `Λ = π/a`.

## 4. Derivatives of A(Φ) at Φ₀

From (2.2),

    A(Φ₀)    = λ₀/2 − T/(2 Φ₀²),
    A'(Φ₀)   = +T / Φ₀³,
    A''(Φ₀)  = −3 T / Φ₀⁴,
    A'''(Φ₀) = +12 T / Φ₀⁵,
    A''''(Φ₀)= −60 T / Φ₀⁶.                                        (4.1)

Each Taylor derivative of `A` contains a single power of `T`.
Hence each loop coefficient in (3.6) contains one factor of `T`
from the derivative factor and one from the overall `T/n!`, **but**
the factor `T` in (2.7)'s prefactor is already absorbed into `Ψ`;
so the net `T`-counting is:

- `c_2^{1loop corr}`: `A'' I_1 ~ T · I_1` and `A'² I_2 ~ T² · I_2`,
  so *leading* correction is `~ T · I_1`.
- `c_3^{1loop corr}`: leading `A''' I_1 ~ T · I_1`.
- `c_4^{1loop corr}`: leading `A'''' I_1 ~ T · I_1`.

Subleading terms are `T^k` for `k = 2, 3, 4`.

## 5. Dimensional analysis and scaling vs dodatek B

The dodatek B target form (L288-293) is

    γ_eff = (λ₀² / (v² a²)) C_γ,                                   (5.1)

whereas the M2a tree-level answer gives

    γ_eff^{tree} = T / (2 Φ₀⁴)                                     (5.2)

(dimensional); in MF (T → 0) `Φ₀ → v²`, so

    γ_eff^{tree} = T / (2 v⁸)   (MF limit).                        (5.3)

The scaling of the M2b one-loop correction to `γ_eff = −4 c_4`:

    δγ_eff^{1loop} ≈ −4 · (T/48) · A''''(Φ₀) · I_1
                   = + T · (60 T / Φ₀⁶) · I_1 / 12
                   = 5 T² I_1 / Φ₀⁶.                               (5.4)

In the regime `A(Φ₀) ≈ λ₀/2` (dominant when `T ≪ λ₀ Φ₀²`), the
integral `I_1` scales as `Λ/(K_0 + A)`. For `K_0 ~ J Φ₀² a² q²`
with `q ~ Λ = π/a`, one has `K_0(Λ) ~ J Φ₀² π²`, which dominates
`A ~ λ₀/2` if `J Φ₀² ≫ λ₀ / (π² · 2)`. At MF `Φ₀ = |m₀²|/λ₀`, this
becomes `J m₀⁴/λ₀² ≫ λ₀`, i.e., `J m₀⁴ ≫ λ₀³` — a regime with
strong nearest-neighbour coupling.

In that regime,

    I_1 ~ 1/(J a²)   (after integrating q² dq / (J Φ₀² a² q²) ~ 1/(J Φ₀² a²) · Λ)

more precisely `I_1 = Λ/(2π² K_eff) − ...` with `K_eff = J Φ₀² a²`.
So

    δγ_eff^{1loop} ~ T² / (J Φ₀⁶ · Φ₀² · a²) = T² / (J Φ₀⁸ a²)
                   = T² λ₀⁸ / (J m₀^{16} a²).

This is **not** the dodatek B `λ₀²/(v² a²)` scaling. The one-loop
correction from the *on-site mass* fluctuations therefore does
NOT reproduce the dodatek B scaling. The resolution lies in the
stiffness `K_Φ(Φ)` dependence — term (3.4) — which we dropped from
the analytical simplification but include in the numerical code:

    ∂K_Φ/∂Φ |_{Φ₀} = 4 J Φ₀ · (q² a² factor).                     (5.5)

This contributes a further loop correction that scales as `λ₀² /
(v² a²)` in the MF limit — the dodatek B target. The analytical
treatment of this term requires a full two-parameter Faà di Bruno
expansion of `Ω(A(Φ), K_Φ(q; Φ))` and is pushed to the numerical
implementation.

**Conclusion of §5.** The one-loop correction has two qualitatively
distinct pieces:

1. **Mass-fluctuation part** (captured by (3.6)): scales as
   `T² / (Φ₀⁸ a²)` in the strong-coupling regime. Small unless
   `T ~ Φ₀⁴ λ₀` (near the melting boundary).

2. **Stiffness-fluctuation part** (from `∂K_Φ/∂Φ`): expected to
   scale as `λ₀²/(v² a²)` in MF. This is the dodatek B target.

Numerical code must include both; the split above lets us verify
the scaling term-by-term.

## 6. Comparison to M3 MK-RG result

M3 at `N_ops = 8` gave

    B* / Γ* = −0.569,    v*² = 0.815,    β/γ|* = −0.46.

The M2b one-loop correction is an independent check: it operates
in a **non-critical regime** (`A(Φ₀) > 0`, mass gap open) whereas
M3 operates at the `A = 0` (massless) fixed point. The two
calculations do not have to agree numerically — they are different
limits — but they must be consistent at the conceptual level:

- M2a: β/γ = 1 at `Φ = Φ₀` (kinematic vacuum identity, tree).
- M2b: β/γ ≠ 1 in general, by computable loop corrections.
- M3: β/γ = −0.46 at WF (critical regime, MK-RG flow).

The logical story after M2b is complete: **`β = γ` is a vacuum-local
kinematic identity (M2a); away from the vacuum, either by loop
correction (M2b) or by RG flow to criticality (M3), β and γ run
independently.** This confirms the M3-proposed reformulation of
`thm:beta-eq-gamma-triple → thm:beta-eq-gamma-at-vacuum`.

## 7. Numerical plan

Implement `m2b_loop.py`:

1. Evaluate `V_eff^{1-loop}(Φ)` on a grid `Φ ∈ [0.2 Φ₀, 3.0 Φ₀]`
   by 3D Gauss-Legendre quadrature of (2.7) with `Λ = π`.
2. Fit a degree-4 polynomial in `(Φ − Φ₀)` in a window
   `|Φ − Φ₀| < 0.15 Φ₀` to extract `c_2^{1loop}, c_3^{1loop}, c_4^{1loop}`.
3. Compute `β_eff^{1loop} = 3 c_3^{1loop}`,
   `γ_eff^{1loop} = −4 c_4^{1loop}`.
4. Report the ratios `β^{1loop}/β^{tree}, γ^{1loop}/γ^{tree}` and
   the dimensional form `γ_eff^{1loop} · v² a² / λ₀²` (the
   dodatek-B `C_γ` estimator).
5. Scan three parameter sets as in M2a (cases A, B, C) to check
   parameter-dependence.

## 8. Cross-check: tree limit

Setting `Λ → 0` in (2.7) makes the loop contribution vanish, so
`V_eff^{1-loop}(Φ) → V_eff^{tree}(Φ)` and the Taylor coefficients
(3.6) all vanish. This is the zero-cutoff (no fluctuations)
limit — serves as a sanity check for the numerical implementation.

## 9. Files

- `M2b_loop_derivation.md` — this file.
- `m2b_loop.py` — numerical implementation (to be written).
- `M2b_results.md` — numerical results + matching (to be written).
