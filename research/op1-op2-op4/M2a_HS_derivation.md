# M2a — H-S derivation of V_eff(Φ): first pass

**Status:** First-pass calculation of tree-level V_eff(Φ)
completed. Key structural result identified.
**Date:** 2026-04-24.
**Method:** collective-field / composite-field (Hubbard–Stratonovich)
integration of the microscopic ŝ degrees of freedom.

---

## 1. The operational target

Dodatek B §app:B-MC (L229-236) defines the object we are deriving:

    V_eff(Φ) := −T ln P(Φ),                                      (1.1)

where `P(Φ)` is the Boltzmann distribution of the block-averaged
composite `Φ = ⟨ŝ²⟩_block`, and the Taylor expansion around its
saddle `Φ₀` gives

    V_eff(Φ) = const + c_2 (Φ−Φ₀)² + c_3 (Φ−Φ₀)³ + c_4 (Φ−Φ₀)⁴ + …

with the identifications (dodatek B L235)

    β_eff = 3 c_3,        γ_eff = −4 c_4.                         (1.2)

The TGP form `U(φ) = (β/3) φ³ − (γ/4) φ⁴` with `φ = Φ/Φ₀`
(dimensionless) arises **provided** `c_1 = c_2 = 0` at the
saddle — which for `c_1` is automatic (saddle) and for `c_2`
requires criticality (massless composite).

## 2. Derivation by composite-field trick

### 2.1 Setup

All ŝ-dependence in `H_Γ` (v2 form) goes through `ŝ² = Φ`.
Replacing `ŝ_i² → Φ_i` everywhere (allowed since `H_Γ` is a
functional of `ŝ²` alone) gives

    H_Γ[Φ] = Σ_i [(m₀²/2) Φ_i + (λ₀/4) Φ_i²]
             + J Σ_⟨ij⟩ A_ij Φ_i Φ_j (Φ_j − Φ_i)²                (2.1)

and the partition function

    Z = ∫Dŝ DΦ · Π_i δ(Φ_i − ŝ_i²) · e^{−β_T H_Γ[Φ]}.             (2.2)

### 2.2 ŝ integration = Jacobian

The per-site ŝ integral (with Φ_i ≥ 0) is the Jacobian of
`ŝ → ŝ² = Φ`:

    ∫dŝ_i · δ(Φ_i − ŝ_i²) = 1 / √Φ_i = Φ_i^{−1/2}.                (2.3)

So the marginal distribution of Φ is

    P(Φ) ∝ Π_i Φ_i^{−1/2} · e^{−β_T H_Γ[Φ]}.                      (2.4)

### 2.3 V_eff(Φ)

Taking `−T ln` of (2.4) gives (per site)

    v_eff(Φ) = (m₀²/2) Φ + (λ₀/4) Φ² + (T/2) ln Φ + const.        (2.5)

This is the **exact tree-level** local V_eff(Φ). All tree-level
information is contained in three ingredients:
- the bare quadratic `(m₀²/2) Φ` (mass),
- the bare quartic `(λ₀/4) Φ²` (self-interaction),
- the **composite-field Jacobian** `(T/2) ln Φ`.

The GL bond does not contribute to the *local* piece — it is
gradient-type and produces `F_kin^geo` (already closed by
`prop:substrate-action`).

### 2.4 Saddle

    v_eff′(Φ) = m₀²/2 + (λ₀/2) Φ + T/(2Φ),
    v_eff′(Φ₀) = 0
      ⇒ λ₀ Φ₀² + m₀² Φ₀ + T = 0
      ⇒ Φ₀ = [ |m₀²| + √(m₀⁴ − 4 λ₀ T) ] / (2 λ₀)                 (2.6)

(valid for `m₀² < 0`, `λ₀ > 0`, `m₀⁴ > 4 λ₀ T`).

Mean-field limit `T → 0`: `Φ₀ → |m₀²|/λ₀ = v²` ✓ (matches
dodatek B L300 `v² = (Jz − m₀²)/λ₀` up to the on-site-only
approximation).

### 2.5 Taylor coefficients at saddle

    v_eff″(Φ₀)   = λ₀/2 − T/(2 Φ₀²),
    v_eff‴(Φ₀)   = T/Φ₀³,
    v_eff⁽⁴⁾(Φ₀) = −3T/Φ₀⁴.

Hence

    c_2 = v_eff″(Φ₀)/2    = (λ₀ Φ₀² − T) / (4 Φ₀²),
    c_3 = v_eff‴(Φ₀)/6    = T / (6 Φ₀³),
    c_4 = v_eff⁽⁴⁾(Φ₀)/24 = −T / (8 Φ₀⁴).                        (2.7)

### 2.6 β_eff, γ_eff via dodatek B definitions

Using (1.2):

    β_eff = 3 c_3 = T / (2 Φ₀³),
    γ_eff = −4 c_4 = T / (2 Φ₀⁴).                                 (2.8)

**Ratio in natural units** (`Φ₀` has dimension of Φ):

    β_eff / γ_eff = Φ₀.                                           (2.9)

## 3. Key result: β = γ ⇔ `Φ₀ = 1` (natural normalisation)

Recall the TGP definition `φ = Φ / Φ₀` — i.e., `φ = 1` at the
vacuum. In that natural dimensionless normalisation the
coefficients of `U(φ) = (β/3) φ³ − (γ/4) φ⁴` are obtained from
the dimensional coefficients `β_eff, γ_eff` by

    β = β_eff · Φ₀³,      γ = γ_eff · Φ₀⁴,                        (3.1)

so that

    β = T/2,            γ = T/2,            β/γ = 1.              (3.2)

**This is the derivation of Route 1 (and, operationally, Route 2)
of `thm:beta-eq-gamma-triple` from microscopic tree-level
dynamics.** The identity `β = γ` at the vacuum is **not** a
fine-tuning — it is the automatic consequence of:

1. the Jacobian `Φ^{−1/2}` of the composite field (a kinematic,
   purely geometric fact about `Φ = ŝ²`);
2. the Taylor expansion of `ln Φ` (mathematical identity);
3. the dimensionless `φ = Φ / Φ₀` normalisation (TGP convention).

All three are structural. No dynamical information beyond the
existence of the saddle is used.

## 4. Comparison to dodatek B substrate map

Dodatek B L288-293 states

    β = (λ₀ / a_sub²) · C_β,
    γ = (λ₀² / (v² a_sub²)) · C_γ.                                (4.1)

Our tree-level result (3.2) gives `β = γ = T/2` (dimensionless
units). To reconcile:

- `β_tree / γ_tree = 1` ≡ `C_β / C_γ · (a_sub² · v² / λ₀) = 1`
  ≡ `C_β / C_γ = λ₀ / (a_sub² v²)` at tree level.

This is **not** `C_β / C_γ ≈ 1.74` (the MK-RG prediction at WF).
The discrepancy `1 vs 1.74` is the loop-correction content that
separates the bare tree-level answer from the renormalised WF
answer. **M2b (loop correction) and M3 (MK-RG flow) are both
needed to bridge this gap.**

## 5. Interpretation

The H-S derivation delivers three distinct things:

| Finding | Status | Implication |
|---|---|---|
| V_eff has a cubic term (c_3 ≠ 0) | structural (from `(T/2) ln Φ`) | Cubic term in `U(φ)` is **automatic** for any Z₂-symmetric microscopic theory in the composite-field formulation |
| `β_eff · Φ₀⁴ = γ_eff · Φ₀⁴` (i.e. `β = γ` in `φ_vac = 1` units) | structural | **Route 1 + Route 2 of thm:beta-eq-gamma-triple unified into one derivation** |
| `β_tree = γ_tree = T/2` | tree-level prediction | Loop corrections will modify this; matching dodatek B's `C_β/C_γ ≈ 1.74` is M2b + M3 work |
| Mass term `c_2 ≠ 0` in general | bare result | `c_2 = 0` requires criticality (tuning `T = λ₀ Φ₀²`). This is the critical-point condition for the WF fixed point. |

## 6. Immediate cross-check

The derivation makes the sharp prediction

    β_eff / γ_eff = Φ₀  (raw, dimensional units)

or equivalently

    β / γ = 1     (in `φ = Φ/Φ₀` normalisation)

at tree level for **any** `(m₀², λ₀, T)` compatible with the
existence of a saddle `Φ₀`. This can be verified numerically in a
1-hour Python script sampling `P(Φ)` from the single-site
Boltzmann distribution and fitting a polynomial around the peak.

## 7. What M2a does **not** deliver

1. **Loop corrections:** the tree-level `β = γ = T/2` is not the
   loop-renormalised value. Corrections from bond fluctuations
   around uniform `Φ` will move the effective `γ` off from
   `T/2` and introduce the `v²` factor of (4.1). Deferred to
   **M2b**.
2. **Irrelevance of φ⁵, φ⁶:** the higher Taylor coefficients
   `c_5, c_6, …` of `v_eff` are

   ```
   c_5 = v^{(5)}(Φ₀)/120 = + T / (5 · 2 · Φ₀⁵) · (4!/5!) = + T / (600 Φ₀⁵) · …
   ```

   (general pattern: `c_k = (−1)^k · (k−1)! · T / (2 · k! · Φ₀^k)` for `k ≥ 3`,
   falling like `Φ₀^{−k}`). At the WF fixed point, these get
   dressed by their respective scaling dimensions `Δ_{φ^k}`; need
   to check against conformal-bootstrap data (`Δ_{ε²} ≈ 3.83`
   already used in op6 M3-c). Deferred to **M2b**.
3. **Value of γ from `(J, λ₀, m₀²)` alone:** `β = γ = T/2` has
   the wrong dependence (T instead of λ₀²/v²) — a signature of
   the missing loops. Deferred to **M2b + M3** convergence.
4. **Dynamical coefficient `W(1) = 7β/3 − 2γ`** (sek08 L1483):
   sits at the non-stationary level and is not touched by
   tree-level static M2a.

## 8. Status of OP-1, OP-2, OP-4 after M2a

| OP | Was | Now | Path forward |
|---|---|---|---|
| OP-1 (form) | postulated | **Cubic + quartic structure + β=γ emerges automatically from composite-field H-S at tree level.** Higher-operator irrelevance: M2b pending. | M2b (loop + scaling dim) |
| OP-2 (all-scale β=γ) | only vacuum (Route 1) | **Route 2 derived: β = γ in dimensionless units is automatic consequence of the Φ-Jacobian.** RG all-scale: M3 pending. | M3 (MK-RG of C_β/C_γ) |
| OP-4 (value of γ) | empirical only | **Tree-level γ = T/2 in dimensionless units.** Needs loops to recover λ₀² dependence. | M2b + M3 convergence |

M2a is thus a **partial but structural** advance: it replaces the
"minimum cubic+quartic" postulate with a **derived** cubic+quartic
form + β=γ identity, at tree level. The gap to the full WF value
of γ is quantified (C_β/C_γ : 1 → 1.74) and localised (loops +
RG flow).

## 9. Code sketch — numerical confirmation of (2.9)

A 1D single-site sanity check (~30 lines, numpy):

```python
import numpy as np

# Parameters (dimensionless units: J=1)
m0sq = -2.0   # in units of J; m0^2 < 0 for ordered phase
lambda0 = 4.0
T = 1.0
N_samples = 10**6

# Single-site Boltzmann: P(s) ∝ exp(-beta*H(s)), H = (m0sq/2)s^2 + (lambda0/4)s^4
def U_s(s): return 0.5*m0sq*s**2 + 0.25*lambda0*s**4

# Metropolis over s, then Phi = s**2, histogram Phi to get P(Phi), fit around peak
# Predict: Phi_0 = (|m0sq| + sqrt(m0sq**2 - 4*lambda0*T)) / (2*lambda0)
# Predict: c_3/c_4 ratio in polynomial fit gives exactly T/(2*Phi_0^3) / (T/(2*Phi_0^4)) = Phi_0
# Equivalently: beta_eff/gamma_eff = Phi_0.

# Expected: beta_eff * Phi_0^3 = gamma_eff * Phi_0^4 = T/2 (should be ~0.5 here)
```

Priority: optional pre-M2b cross-check. ~1 hour wall-clock.

## 10. Decision point

Two parallel next tracks:

- **M2b — loop correction:** perturbative expansion of `v_eff(Φ)`
  beyond tree level, via bond-fluctuation expansion around uniform
  Φ. Expected to produce the `v²` dependence in `γ` (dodatek B
  eq:gamma-from-substrate).
- **M3 — MK-RG extension:** tracks operator flow and measures
  `C_β/C_γ` at WF. Should come out to ≈ 1.74 if the programme
  is internally consistent.

**Recommendation:** start M3 next — it's numerical, self-contained,
and gives an immediate sanity check on the whole framework. M2b can
follow once M3 gives the target.
