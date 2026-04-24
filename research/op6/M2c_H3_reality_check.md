# OP-6 / M2-c — Reality check: what does "H₃" actually give?

**Status:** CLOSED, 2026-04-24. Outcome: **H₃ bilinear rejected
at 5.5σ as α=2 target** (Part E numerical confirmation).
**Predecessor:** `M2b_results.md`.
**Purpose:** Before accepting the "M1-A pivot = replace H₁ with H₃" as
the recommended fallback, we verify that H₃ *as stated in M1* actually
produces α=2. It does not. The α=2 result of `dodatekB_substrat.tex`
uses a **different** object than the H₃ inventoried in M1. This
document corrects the record.

---

## 1. Two different "quartic-bond" forms

In lattice scalar field theory there are two natural ways to involve
a density-dependent kinetic coupling, and they are frequently
conflated.

### Form A: quartic (bilinear in Φ) bond

    H₃_bilinear   =   on-site   −   J  Σ_⟨ij⟩  (φ_i φ_j)²
                  =   on-site   −   J  Σ_⟨ij⟩  Φ_i Φ_j ,   Φ = φ²

This is the form *literally written* in `M1_H_Gamma_canonical.md` §1
under "H₃". In the Φ variable it is the standard Ising bilinear
bond, with continuum kinetic term

    F_kin[Φ]  =  (J a^{2−d} / 2) ∫ (∇Φ)² d^d x .

That is K(Φ) = const. Converting to the φ variable via Φ = φ²,
(∇Φ)² = 4 φ² (∇φ)², one obtains

    F_kin[φ]  =  2 J a^{2−d} ∫ φ² (∇φ)² d^d x ,
    K(φ)     =  4 J a^{2−d} · φ² ,
    α(φ)    =  K'·φ/(2K)  =  **1** .

So form A gives α = 1, **not 2**.

### Form B: GL functional with density-dependent stiffness

    F_kin^{GL}[φ]   =   Σ_⟨ij⟩  K_{ij} (φ_j − φ_i)² ,
    K_{ij}          =   J (φ_i φ_j)²

This is what `dodatekB_substrat.tex:530–641` (`prop:substrate-action`)
actually uses, and it is a **Ginzburg–Landau functional**, not a
bond-term Hamiltonian. Leading expansion:

    (φ_j − φ_i)²   ≈   a² (ê·∇φ)²  +  O(a³) ,
    K_{ij}         ≈   J φ^4  +  O(a) ,

so summing over bonds and integrating,

    F_kin^{GL}[φ]  =  2 d J a^{2−d} ∫ φ⁴ (∇φ)² / 2 · d^d x / d
                   =  (K_geo / 2) ∫ φ⁴ (∇φ)² d^d x ,
    K(φ)           =  K_geo · φ⁴ ,
    α(φ)          =  **2** . ✓

This is the correct derivation of the TGP kinetic form.

## 2. Why this matters

The pivot target recommended in `M2b_results.md` was stated as
"change axiom to H₃". Because the M1 inventory identifies H₃ with
form A, the recommendation as written is **wrong**: form A gives
α=1, not α=2.

The correct pivot target is form B, which we now designate **H_GL**:

    H_GL  =  on-site terms (as in H₁)
            +  Σ_⟨ij⟩  J (φ_i φ_j)² (φ_j − φ_i)²

where φ_i = √(Φ_i / Φ₀) and Φ_i = ⟨ŝ_i²⟩.

H_GL is a **GL effective functional** with density-dependent
stiffness, not a microscopic spin Hamiltonian. Adopting it as
a core-paper axiom is conceptually more invasive than adopting
form A would have been:

- H₁ (bilinear ŝ bond): natural microscopic spin Hamiltonian.
- H₃ form A: natural spin Hamiltonian (Ising in Φ). Wrong α.
- H_GL form B: *effective-theory* axiom with a built-in density
  stiffness; not a microscopic spin bond.

So the clean axiomatic identification "core paper's substrate
Hamiltonian = microscopic Ising-like system" becomes strained
under the H_GL pivot. The axiom is now "substrate effective
kinetic functional is GL with density-dependent stiffness",
which is more of a postulate about the EFT structure than a
microscopic Hamiltonian.

## 3. Consequences for the OP-6 programme

### 3.1 M1 inventory is incomplete

The inventory in `M1_H_Gamma_canonical.md` §1 lists H₃ and H₄
(H₄ ≡ H₃ without on-site terms) both as "quartic-bond forms".
Neither is H_GL. The correct list of relevant microscopic /
effective-theory objects is four:

| Label | Form | Bond structure | K(φ) | α(φ) |
|-------|------|----------------|------|------|
| H₁ | bilinear ŝ bond | −J ŝ_i ŝ_j | const · 1 | 0 (in ŝ) |
| H₃_bilinear | bilinear φ² bond | −J (φ_iφ_j)² | const · φ² | **1** |
| H_GL | GL stiffness | Σ J(φ_iφ_j)²(φ_j−φ_i)² | K_geo · φ⁴ | **2** |
| — | quartic relative bond | −J (φ_j − φ_i)⁴ | const · (∇φ)²-anomalous | — |

Only H_GL gives α = 2.

### 3.2 M2-b verdict unchanged

M2-b's negative result (K_eff ∝ Φ rejected at 3σ from 1D MC of H₁)
still holds. H₁ does not generate α=2. Nothing in this document
changes that conclusion.

### 3.3 Pivot target corrected

The recommended pivot is not "replace H₁ with H₃ bilinear" (that
gives α=1). The recommended pivot — if one is willing to adopt a
GL effective-theory axiom — is **replace H₁ with H_GL**:

    Axiom (v2, M1-A'):  F_kin^{substrate}[φ] =
        (K_geo/2) ∫ φ⁴ (∇φ)² d^dx,
    or microscopically:
        F = Σ_⟨ij⟩ J (φ_iφ_j)² (φ_j - φ_i)²
                   + on-site GL terms.

This is the *direct* axiomatization of what prop:substrate-action
uses, with no derivation gap. Conceptually, it promotes the GL
kinetic functional from a "derived object" to an "axiom about
the effective theory".

### 3.4 Open question for any v2

Is H_GL natural as a substrate-level axiom? Arguments for:

- It is the unique local, Z₂-even, derivative-expansion-leading
  kinetic functional in class (C1)-(C3) with density-dependent
  stiffness. (Theorem thm:alpha2 in the core paper is effectively
  a classification inside (C1)-(C3).)
- It is exactly what dodatekB's proof already uses.

Arguments against:

- It is not a microscopic bond term; it is an effective GL
  functional. TGP's "substrate is a spin-like lattice system"
  motivation is weakened.
- The bond structure Σ K_{ij}(φ_j-φ_i)² is phenomenologically
  natural in GL theory but is not a standard starting point in
  RG: normally one proves K flows to a certain density form, not
  assumes it.

### 3.5 A new potential route: H₁ → H_GL via RG?

An alternative that preserves the "microscopic = H₁" narrative:
ask whether H₁ produces H_GL as an *effective* functional under
block coarse-graining. This is distinct from M2-b (which asked
whether H₁ produces K_eff(Φ) ∝ Φ under envelope coarse-graining).

Concretely: at scale L_B, the effective kinetic functional derived
from H₁ might have the structure Σ K_B(φ_i, φ_j)·(φ_j − φ_i)²,
where K_B is generated by the RG flow. If K_B → K_geo · φ^4 in the
IR, then H₁ is consistent with the GL axiom at the effective level.

This is a *stronger* claim than M2-b's and correspondingly harder to
test. It is a legitimate M3-class problem.

## 4. Numerical verification: form A does NOT give α=2

Part E of `m2b_envelope_stiffness_1d.py` (added 2026-04-24) runs
form A (bilinear Φ-bond `−J Σ (s_is_j)²`) with the same OZ
framework as Part C, scanning (β, m², λ) subject to the stability
constraint λ > 4J.

**Scan:** β ∈ {0.2, 0.3, 0.5, 0.7, 1.0}, (m², λ) ∈
{(0.5, 6), (1.0, 6), (2.0, 6), (3.0, 8), (4.0, 10)} with J = 1,
N = 1024, n_therm = 3000, n_measure = 8000 sweeps.

**Result:**

    K_eff(⟨Φ⟩) = 9.24 · ⟨Φ⟩^(−2.49 ± 0.63)   (n = 12 valid points)

- Test of "H₃ bilinear ⇒ α=2" (target p_Φ = +1):
  z = (−2.49 − 1.00)/0.63 = **−5.5 σ** → rejected.
- Test of bare-bond prediction (target p_Φ = 0):
  z = −2.49/0.63 = −3.9 σ → rejected too.

### 4.1 Interpretation

The **primary finding** is the 5.5σ rejection of p_Φ = +1: the
measured K_eff(⟨Φ⟩) for H₃-bilinear is very far from linear in Φ.
This *directly confirms* the analytical claim of §1: form A is
not the α=2 target. Combined with the 3σ rejection of p_Φ = +1
from Part C (H₁ envelope coarse-graining), both microscopic
candidates in M1's original "H₁ vs H₃" dichotomy are numerically
incompatible with α(φ) = 2. The only remaining α=2-generating
object is H_GL.

The **secondary finding** is that the measured slope p_Φ ≈ −2.5
differs significantly from the bare-bond prediction p_Φ = 0. Two
non-exclusive explanations:

1. **Non-Gaussian fluctuation renormalization.** The MC is in the
   strong-coupling regime λ > 4J (forced by stability), so the
   effective K_eff is heavily dressed by loop corrections beyond
   the tree-level "−J Φ_iΦ_j → const·(∂Φ)²" algebra. The dressed
   K_eff(⟨Φ⟩) resembles a Gaussian-like 1/Φ² behaviour (as in
   Part B of the H₁ scan), consistent with the composite-operator
   character of the Φ = s² correlator we measure.

2. **Short correlation length.** Many OZ fits fail (NaN) because
   ξ_Φ < 1 lattice unit at the strong-coupling parameters
   required. The 12 valid points are concentrated in a narrow
   range of ⟨Φ⟩ (0.17–0.77), amplifying the sensitivity of the
   slope to small numerical effects.

Neither of these secondary issues affects the primary conclusion:
**p_Φ is very far from +1, so form A (H₃ bilinear) does not
generate α=2 by envelope coarse-graining.**

### 4.2 Comparison table (numerical)

| Hamiltonian | Bond | Measured p_Φ | σ from p_Φ = +1 |
|-------------|------|--------------|-----------------|
| H₁ (bilinear ŝ bond)     | −J ŝ_iŝ_j            | −0.28 ± 0.43 | 3.0 σ (Part C) |
| H₃ bilinear (bilinear Φ bond) | −J (s_is_j)² = −J Φ_iΦ_j | −2.49 ± 0.63 | 5.5 σ (Part E) |
| H_GL | Σ J(s_is_j)²(s_j−s_i)² | +1 (by algebraic construction) | — |

The "target" column highlights what α = 2 requires: p_Φ ≈ +1.
Neither H₁ nor H₃-bilinear is close; only H_GL attains it, and
that by construction.

See `m2b_scan_results.txt` for the raw table of (β, m², λ, ⟨Φ⟩,
χ, ξ², K, K_err) for all 25 Part-E scan points.

## 5. What this adjusts in the written record

The following documents need a corrective note:

1. `M1_H_Gamma_canonical.md` §1 — add a new entry "H_GL" for the
   GL-stiffness form, distinguishing it from "H₃ bilinear". §5 and
   §6 should reference H_GL for α=2, not H₃.
2. `M2b_results.md` §6-7 — correct "recommended pivot: axiom = H₃"
   to "recommended pivot: axiom = H_GL (GL effective functional)".
   Note that the ultimate axiomatic cost is slightly higher than
   the M1-B path but is still the only route compatible with α=2.
3. `README.md` — reflect the distinction and include H_GL in the
   file listing.
4. `cg_strong_numerical.py` (core paper) — the "CG-4 = K ∝ Φ"
   result was based on the wrong K ~ ξ identification. That
   result does not establish K(Φ) ∝ Φ for H₃ bilinear; form A
   analytically gives K(Φ) = const. This should be documented
   in KNOWN_ISSUES (but is secondary: cg_strong_numerical.py is
   now labelled as "exploratory" in the core-paper retraction).

## 6. Bottom line

The v2-core pivot is conceptually *clean* only if the axiom is H_GL
(the GL effective kinetic functional), because that is what the α=2
proof actually uses. A "pivot to H₃ bilinear" would be algebraically
wrong **and** is now numerically ruled out (5.5σ, Part E).

OP-6 therefore reduces to two live options, both non-trivial:

- **Option 1 (H_GL axiom):** accept that TGP's substrate is
  specified by a GL effective kinetic functional, not a microscopic
  spin Hamiltonian. Clean at the axiom level, costly at the
  interpretive level.
- **Option 2 (M1-B revisited):** prove that H₁ generates H_GL as an
  effective functional under RG. Clean interpretively, but this is
  a *harder* problem than M2-b (we'd need to show not only that the
  effective kinetic stiffness depends on Φ, but that it has the
  specific functional form K ∝ φ⁴ · (φ_j-φ_i)²). Currently unknown.

The honest v1 stance — OP-6 open — remains correct. v2 is still a
pending deliverable contingent on either Option 1 or Option 2
closing.

### Current numerical status

| Hamiltonian  | Measured p_Φ        | Verdict for α=2 target (p=+1) |
|--------------|---------------------|-------------------------------|
| H₁           | −0.28 ± 0.43 (n=20) | rejected at 3.0 σ (M2-b Part C) |
| H₃ bilinear  | −2.49 ± 0.63 (n=12) | rejected at 5.5 σ (M2-c Part E) |
| H_GL         | +1 exactly          | constructed to satisfy it       |

This closes M2. The only microscopic/effective object in the
current inventory compatible with α=2 is H_GL, which is an
effective-theory axiom rather than a spin-Hamiltonian axiom.
