# OP-6 / M1 — Canonical H_Γ and the variable-convention audit

**Status:** superseded in part by `M2c_H3_reality_check.md`, 2026-04-23.
**Author:** internal review after withdrawal of `a1_alpha2_frg_synthesis.py`
(see `tgp-core-paper/KNOWN_ISSUES.md`).

**Purpose:** Before any re-derivation of α = 2 from the discrete
substrate, we consolidate *one* canonical Hamiltonian in *one*
canonical variable set, list every currently circulating variant,
and state precisely what M2 has to derive.

> **⚠ Corrective note (2026-04-23, see M2c).** §1, §4.2, §5, §6, §8
> of this document treat "H₃ = −J Σ (φ_iφ_j)²" as the α=2 target.
> That identification is **algebraically wrong**: the bilinear
> Ising-like bond `−J Σ (φ_iφ_j)² = −J Σ Φ_iΦ_j` gives a
> Φ-independent stiffness K(Φ) = const, equivalently K(φ) ∝ φ²,
> equivalently **α(φ) = 1**, not 2. The object that actually gives
> α(φ) = 2 is the Ginzburg–Landau functional
> `F = Σ_⟨ij⟩ J(φ_iφ_j)² (φ_j − φ_i)²` (here denoted **H_GL**),
> which is what `prop:substrate-action` in dodatekB_substrat.tex
> actually uses. Accordingly the §5 "M1-A = pivot to H₃" decision
> is void; the corrected pivot target is H_GL, and the full
> analysis is in `M2c_H3_reality_check.md`. The rest of this
> document is retained for historical traceability.

---

## 1. The four Hamiltonians currently in circulation

Inventory from a systematic search of `TGP_v1/`. All four have
been used at some point as "H_Γ" in proofs / scripts / LaTeX.

### H₁ — *published-core-paper form* (linear bond)

```
H_Γ = Σ_i [ π̂_i² / (2μ)  +  (m₀²/2) ŝ_i²  +  (λ₀/4) ŝ_i⁴ ]
    − J Σ_{⟨ij⟩} A_{ij} ŝ_i ŝ_j
```

- Source: `core/sek01_ontologia/sek01_ontologia.tex:737–743`,
  `eq:H-Gamma`, Axiom 0 of core paper.
- Also in `tgp-core-paper/paper/tgp_core.tex` eq:H-Gamma (published).
- Structure: Z₂ scalar on a graph with a **linear** bond term
  `−J ŝ_i ŝ_j` (not squared).

### H₂ — *static renormalization form* (linear bond, no kinetic)

```
H_Γ = Σ_i [ (m₀²/2) ŝ_i²  +  (λ₀/4) ŝ_i⁴ ]  −  J Σ_{⟨ij⟩} ŝ_i ŝ_j
```

- Source: `axioms/substrat/dodatekB_substrat.tex:19–23`, `eq:B-H`.
- Used for RG/MC/thermodynamic studies.
- Identical to H₁ with `π²/(2μ)` dropped (static sector).
- Consistent with H₁.

### H₃ — *geometric-coupling form* (quadratic bond)

```
H_Γ^(geo) = [on-site terms]  −  J Σ_{⟨ij⟩} (φ_i φ_j)²
```

- Source: `axioms/substrat/dodatekB_substrat.tex:530–600`,
  `prop:substrate-action`.
- In this model the bond is `(φ_i φ_j)²`, not `φ_i φ_j`.
- *This is what the derivation of α = 2 in dodatekB uses*.
- **Not equal to H₁ or H₂.**

### H₄ — *minimal coupling form* (letter)

```
H_Γ = −J Σ_{⟨ij⟩} (φ_i φ_j)²
```

- Source: `tgp_letter.tex:203`.
- Bond-only version of H₃; typically used as a shorthand.
- Equivalent to H₃ when on-site terms are suppressed / absorbed
  in the fixed-point analysis.

---

## 2. The inconsistency

The published core paper declares H₁ as axiom, then claims in
Theorem `thm:alpha2` and in row 16/21 of the results table that
α = 2 follows algebraically from the substrate.

The actual derivation of α = 2 lives in `prop:substrate-action`
(dodatekB) and uses H₃, **not H₁**. That proposition starts from
the *geometric coupling* K_{ij} = J(φ_i φ_j)² and obtains the
continuum functional F_kin = (K_geo/2) φ⁴ (∇φ)². This cannot be
derived from H₁ by simple coarse-graining — H₁ has a linear bond
and in the naïve continuum limit produces a *constant* kinetic
prefactor (see §4 below), not a φ⁴ one.

The two Hamiltonians are therefore **different microscopic models**.
The published axiom is H₁; the working derivation of α = 2 requires H₃.
This is the source of the circular reference we hit in
`a1_alpha2_frg_synthesis.py`.

**Summary of the inconsistency (corrected 2026-04-23 per M2c):**

| Object | Bond / Form | Continuum K(φ) | α (in φ-var.) |
|-------|------|----------------|--------------|
| H₁ / H₂ | `−J ŝ_i ŝ_j` | const (in ŝ); ∝ 1/Φ | 0 (in ŝ); −½ (in Φ) |
| H₃ / H₄ (bilinear) | `−J (φ_iφ_j)² = −J Φ_iΦ_j` | ∝ φ² | **1**, not 2 |
| **H_GL** | `Σ J(φ_iφ_j)²(φ_j−φ_i)²` (GL functional) | K_geo · φ⁴ | **2** ✓ |

The core paper effectively uses α = 2 from the H_GL line of
`prop:substrate-action` while declaring H₁ (bilinear ŝ bond) as
axiom. H_GL is not the same as H₃: the former is a GL effective
functional with density-dependent stiffness, the latter is a
bilinear Ising-like bond in Φ. The original version of this
summary incorrectly equated the two.

---

## 3. Variable-convention pinning

Let us fix one convention for the rest of the OP-6 programme.

- `ŝ_i` — substrate operator on site `i ∈ V`. Dimensionless,
  Z₂-odd: `ŝ_i ↦ −ŝ_i`.
- `⟨ŝ_i⟩ = σ_i` — mean-field amplitude (real, Z₂-odd).
- `Φ_i ≡ ⟨ŝ_i²⟩` — local mean-squared amplitude. Z₂-even.
  This is the block-precursor of the continuum field Φ.
- `Φ(x)` — continuum field after block averaging:
  `Φ(x) = (1/N_B) Σ_{i ∈ B(x)} ⟨ŝ_i²⟩`.
- `Φ₀` — vacuum value of Φ. Single TGP scale.
- `φ ≡ Φ/Φ₀` — dimensionless macroscopic field. Positive.
- `g ≡ √φ = √(Φ/Φ₀)` — "amplitude" variable used in the
  soliton code (`f_kin(g) = 1 + 2α ln g`).

Relation: `Φ = Φ₀ φ = Φ₀ g²`.

### The α convention

Given a kinetic functional
F[ψ] = (1/2) ∫ K(ψ) (∇ψ)² d^d x,
the corresponding Euler–Lagrange equation is
∇² ψ + (K'(ψ) / (2 K(ψ))) (∇ψ)² = 0.
We call **α(ψ) ≡ K'(ψ) · ψ / (2 K(ψ))** the TGP kinetic exponent
*in the ψ variable*.

For `K(ψ) = ψ^n`: α(ψ) = n/2. Thus:

- K(φ) = φ⁴ → α(φ) = 2.
- K(g) = g² → α(g) = 1.
- K(Φ) = Φ → α(Φ) = 1/2.
- K(Φ) = 1/Φ → α(Φ) = −1/2.

**α is variable-dependent**. The core paper statement
"α = 2" is implicitly in the φ variable.

---

## 4. Change-of-variable audit

We now check what `a1_alpha2_frg_synthesis.py` actually computed.

### 4.1. From H₁ / H₂ (linear bond)

Coarse-graining `−J Σ_{⟨ij⟩} ŝ_i ŝ_j` to leading order in the
lattice spacing `a`:

- ŝ_i → ŝ(x), ŝ_j = ŝ(x + a ê) → ŝ(x) + a ê·∇ŝ + (a²/2)(ê·∇)²ŝ + …
- ŝ_i ŝ_j = ŝ² + a ŝ ê·∇ŝ + (a²/2) ŝ (ê·∇)² ŝ + …
- Σ_{⟨ij⟩} ŝ_i ŝ_j / (# bonds) ~ ŝ² + (a²/2d) Σ_μ ŝ ∂²_μ ŝ

The continuum of `−J Σ_{⟨ij⟩} ŝ_i ŝ_j` thus produces
- A mass-like term `const · ŝ²`
- A kinetic term **of the form `(J/2) a^{2−d} (∇ŝ)²`** after
  integration by parts — coefficient **constant in ŝ**.

⇒ K(ŝ) = const. In the variable Φ = ŝ²:
- (∇ŝ)² = (∇Φ)² / (4Φ)
- K(ŝ) (∇ŝ)² = const · (∇Φ)² / (4Φ)
- ⇒ **K(Φ) = const / Φ** ⇒ α(Φ) = −1/2.

This is exactly the `K_1(Φ) = Z_0/(4Φ)` result from the
script. It is algebraically correct and corresponds to a
*free scalar field*, not α = 2 in any variable.

(If instead we keep the variable ŝ: K(ŝ) = const ⇒ α(ŝ) = 0.)

### 4.2. From H₃ / H₄ (bilinear-in-Φ bond) — **CORRECTED**

**Note (2026-04-23, per M2c):** the algebra below was written as if
coarse-graining `−J Σ (φ_iφ_j)²` produced `(K_geo/2) φ⁴ (∇φ)²`. It
does not. In the Φ variable (Φ = φ²) this bond *is* the Ising
bilinear `−J Σ Φ_iΦ_j`, so the standard lattice-gradient expansion
yields K(Φ) = const, i.e. K(φ) = 4·const · φ², i.e. **α(φ) = 1**.

Correct leading-order expansion:

- (φ_iφ_j)² = Φ_iΦ_j with Φ ≡ φ²
- Using Φ_iΦ_j = (Φ_i² + Φ_j²)/2 − (Φ_j − Φ_i)²/2:
  −J Σ_⟨ij⟩ Φ_iΦ_j = (on-site · Φ²) + (J/2) Σ_⟨ij⟩ (Φ_j − Φ_i)²
- The bond piece Σ_⟨ij⟩ (Φ_j − Φ_i)² ≈ a² Σ_⟨ij⟩ (ê·∇Φ)²
  →  d · a^{2−d} ∫ (∇Φ)² d^dx.

So F_kin[Φ] = (J a^{2−d}/2) ∫ (∇Φ)² d^dx up to lattice factors, i.e.

  **K(Φ) = const**   ⇒  α(Φ) = 1/2,
  **K(φ) ∝ φ²**      ⇒  **α(φ) = 1**.

The object that actually gives α(φ) = 2 is the GL functional H_GL
(§4.3 below), not the bilinear bond written here as H₃. The
original (incorrect) claim, preserved for diff-trace:

> ~~⇒ K(φ) = K_geo · φ⁴ ⇒ α(φ) = 2.~~  **Wrong for this bond.**

### 4.3. From H_GL (GL functional with density-dependent stiffness)

This is what `prop:substrate-action` in dodatekB actually uses.
Define

    F_kin^{GL}[φ] = Σ_⟨ij⟩ K_{ij} (φ_j − φ_i)²,  K_{ij} = J (φ_i φ_j)².

Leading expansion:

- (φ_j − φ_i)² ≈ a² (ê·∇φ)²
- K_{ij} ≈ J φ⁴ + O(a)
- Σ_⟨ij⟩ → (d · N) bonds per dimension; Σ → a^{−d} ∫ d^dx

hence

    F_kin^{GL}[φ] = (K_geo/2) ∫ φ⁴ (∇φ)² d^dx,
    K_geo = 2d J a^{2−d},
    K(φ) = K_geo · φ⁴,
    **α(φ) = 2**. ✔

H_GL is an effective GL functional, not a microscopic bond
Hamiltonian. See `M2c_H3_reality_check.md` §1 for details.

### 4.3. Reconciling the soliton code

The soliton code uses `f_kin(g) = 1 + 2α ln g`, i.e. *logarithmic*
rather than power-law kinetic prefactor. Neither H₁ nor H₃ gives
this by leading-order coarse-graining. A log correction would have
to come from:

- loop-level / RG-induced anomalous dimension η > 0,
- a Weyl-anomaly-like contribution from the effective metric
  `g_ij ∝ Φ^p · δ_ij`, or
- an explicit resummation of higher-order bond expansions.

The honest position is that `f_kin(g) = 1 + 2α ln g` is **not
currently derived from H_Γ at all** — it is postulated in the
effective field theory, motivated by the φ⁴ pure-power form
(which gives α=2 at leading order).

---

## 5. M1 decision: which Hamiltonian is canonical

> **Corrective note (2026-04-23, M2c):** Option M1-A below was
> written under the mistaken premise that the bilinear "H₃" bond
> gives α=2. It does not. The truly α=2-generating axiom is
> **H_GL** (GL functional). Read "M1-A" in what follows as a
> *historical* option that is algebraically void; the corrected
> pivot candidate is **M1-A′ = promote H_GL to axiom**. See M2c §3.
> Option M1-B is unchanged.

Two clean ways forward:

### Option M1-A. Declare H₃ / H₄ canonical  ~~[void, see M2c]~~

- Replace `eq:H-Gamma` in `core/sek01_ontologia/sek01_ontologia.tex`
  and in the published core paper with
  `H_Γ = [on-site] − J Σ (ŝ_iŝ_j)²`.
- dodatekB `prop:substrate-action` becomes a direct theorem from
  the axiom, no gap.  **(Algebraically incorrect — this bond gives
  α=1, not α=2.)**
- **Cost:** a substantive axiomatic change to the published paper.
  Triggers a v2 release. Also: Z₂ symmetry `ŝ ↦ −ŝ` of the
  quadratic bond is less restrictive (the bond is already
  Z₂-even) — the symmetry structure of TGP changes qualitatively.

### Option M1-A′. Declare H_GL canonical  [corrected, per M2c]

- Adopt as a substrate-level axiom:
  `F_kin^{sub}[φ] = (K_geo/2) ∫ φ⁴ (∇φ)² d^dx`, equivalently
  the microscopic GL functional
  `F = Σ_⟨ij⟩ J (φ_iφ_j)² (φ_j − φ_i)² + on-site terms`.
- `prop:substrate-action` becomes a direct theorem from the axiom.
- **Cost:** the axiom is an *effective-theory* statement, not a
  microscopic spin Hamiltonian. Conceptually more invasive than
  "M1-A" as originally imagined. See M2c §3 for trade-offs.

### Option M1-B. Keep H₁ canonical; derive H₃ as effective

- Keep `eq:H-Gamma` as published.
- Treat `H₃` as an *effective* Hamiltonian produced by integrating
  out short-wavelength modes of H₁. This is the natural RG picture.
- Required derivation: under block-averaging to scale `L_B`, the
  bilinear bond of H₁ renormalises into a quartic bond that
  *dominates* the IR. This is a non-trivial claim that needs a
  genuine RG computation.
- **Cost:** the derivation is hard and may fail. But if it succeeds,
  the core paper needs no axiomatic change — only a strengthening
  of the connection between H₁ and `prop:substrate-action`.

### Recommendation (SUPERSEDED 2026-04-24)

> **Update 2026-04-24:** this recommendation is **reversed**. The
> M1-B path was shown to fail on three independent grounds (M2-b
> 1D envelope 3σ reject, M3-a 1D block-RG 5.5σ reject, M3-c 3D
> Ising WF scaling-dimension argument). The now-authoritative
> recommendation is **M1-A′ (declare H_GL canonical)**. The v2
> of the core paper (tgp-core-paper/paper/tgp_core.tex)
> implements this pivot as of 2026-04-24. See
> `M3a_results.md` §7 and `KNOWN_ISSUES.md` 2026-04-24 entry.

~~**Option M1-B** is correct. Reasons:~~

~~1. The Z₂ structure of H₁ is the foundation of the entire TGP
   axiom system; changing it is more invasive than changing the
   substrate Hamiltonian.~~
~~2. If the quartic bond cannot be derived from the bilinear bond
   by RG, then the α = 2 result is *not* a consequence of the
   stated TGP axioms — it is a consequence of a different
   microscopic model. We should find this out, not mask it.~~
~~3. In the RG picture, the appearance of a dominant φ⁴·(∇φ)²
   term at the IR fixed point is a *specific and testable*
   claim, which is exactly what M3 (numerical FRG) should verify.~~

**Post-closure recommendation (2026-04-24): Option M1-A′.** Reasons:

1. **M1-B is numerically and analytically ruled out** — the
   quartic bond cannot emerge from the bilinear bond under RG
   at 3D Ising WF. Three independent probes (M2-b, M3-a, M3-c)
   agree. Continuing M1-B would be unscientific.
2. **The ontological chain "nothingness + Z₂ → reality" is
   preserved** by M1-A′. The microscopic variable ŝ_i stays
   Z₂-odd; Φ_i = ŝ_i² stays Z₂-even; the bond
   `J ŝ_i² ŝ_j² (ŝ_j² − ŝ_i²)²` vanishes on empty-substrate
   configurations (ŝ = 0), so "no substrate → no geometry"
   becomes a *property of the axiom*, not a derived consequence.
3. **The GL bond is actually more ontologically principled than
   the bilinear bond** — in v1 (H₁), substrate bonds were
   non-zero even when the field was zero (bond =
   `−J·0·0 = 0` formally, but activation was not tied to Φ).
   In v2 (H_GL), the bond is quadratic in Φ, so empty space
   cannot propagate.
4. `prop:substrate-action` in dodatekB has used the GL bond since
   2026-03-24. The pivot elevates the object already used in the
   formal appendix to the axiom level, closing the v1 gap.

**Canonical H_Γ (v2):**

    H_Γ = Σ_i [π²/2μ + (m₀²/2) ŝ² + (λ₀/4) ŝ⁴]
          + J Σ_⟨ij⟩ A_ij · ŝ_i² ŝ_j² · (ŝ_j² − ŝ_i²)²

with J > 0, λ₀ > 0, A_ij ∈ {0, 1}. Continuum form (see
Prop. `prop:substrate-action`): `F_kin^geo[φ] =
(K_geo/2) ∫ φ⁴ (∇φ)² d^dx`.

---

## 6. What M2 must derive

> **SUPERSEDED 2026-04-24.** §§6–8 below describe the M2–M5
> programme *under the now-obsolete M1-B recommendation* (derive
> α=2 from bilinear H₁). That programme ran and **failed** on
> three independent grounds; the contingency of §8 ("if M2-a
> fails") has now triggered. The executed and closed record is
> in `M2a_analytical_sketch.md`, `M2b_results.md`,
> `M2c_H3_reality_check.md`, `M3a_results.md`,
> `M3c_scaling_dimensions.md`. The v2 resolution is **M1-A′**
> (§5 above), and the doc-update list of §7 is correspondingly
> replaced by the edits actually made to
> `tgp-core-paper/paper/tgp_core.tex`, `KNOWN_ISSUES.md`, this
> file's §5, and `TGP_v1/research/op6/v2_pivot_summary.md`.

§§6–8 are retained verbatim for historical traceability.

---

With H₁ fixed as canonical, the rigorous programme for α = 2 is:

**(M2-a) Effective bond coupling.**
Show that block-averaging H₁ to scale L_B generates an effective
Hamiltonian `H_eff[ŝ_B]` whose bond term is dominated by a
quartic piece `−J_eff Σ (ŝ_B,i ŝ_B,j)²` at L_B → ∞. Equivalently:
the bilinear bond is irrelevant in the IR, and the leading
relevant bond operator on the symmetric side of the fixed point
is the quartic one.

This is an explicit RG claim. It can be checked perturbatively
(LPA' expansion around the Wilson–Fisher fixed point, with the
bond-operator basis enlarged to include both ŝ_iŝ_j and
(ŝ_iŝ_j)²). It can also be checked numerically on the lattice
(cg_strong_numerical.py is the kind of script that can be
extended to measure this, by reading off the effective bond
spectrum after block averaging).

**(M2-b) Continuum form.**
From the effective quartic bond, derive the continuum kinetic
functional F_kin = (K_eff/2) φ⁴ (∇φ)² as in `prop:substrate-action`.
This is standard coarse-graining once the quartic bond is given.

**(M2-c) α = 2 via Euler–Lagrange.**
From F_kin above, obtain the E–L operator `∇²φ + 2(∇φ)²/φ = 0`,
i.e. α(φ) = 2. This is algebraic (it is the content of Theorem
`thm:alpha2` in the core paper).

Together (M2-a) + (M2-b) + (M2-c) give α = 2 *as a derived*
statement from H₁, not as a postulate. This replaces the
withdrawn A1–A5 chain.

---

## 7. Documentation follow-ups (non-blocking for M2)

These are edits we will do *if and when* the M2–M5 programme
closes successfully (→ v2 of the core paper):

1. `core/sek01_ontologia/sek01_ontologia.tex` — add a remark
   immediately after `eq:H-Gamma` noting that the quartic bond
   in `prop:substrate-action` (dodatekB) is an *effective* coupling
   obtained by RG, not an independent axiom. Cross-reference
   M2-a.

2. `core/sek00_summary/sek00_summary.tex:106–116` — rewrite the
   "Nota o dualizmie α": the two forms K(φ)=φ⁴ and K_sub(g)=g²
   are not independent "dualities" but different coarse-grainings
   of the same H₁.

3. `axioms/substrat/dodatekB_substrat.tex:530` — restate
   `prop:substrate-action` as a theorem conditional on the
   effective quartic bond, cross-referencing M2-a.

4. `tgp_letter.tex:203` — clarify that H_Γ = −J Σ (φ_iφ_j)² is
   shorthand for the effective bond of H₁, not the fundamental
   Hamiltonian.

5. Published core paper (would require v2) — add an explicit
   derivation paragraph around `thm:alpha2` linking (C1)–(C3) to
   H₁ via the RG-derived quartic bond.

Until M2 closes, these remain *wishlist* edits. If M2 fails,
a different set of edits applies (see §8).

---

## 8. Contingency: if M2-a fails

If it turns out that block-averaging H₁ does *not* produce a
dominant quartic bond at the IR, then α = 2 is not a consequence
of H₁ at all. In that case:

- `thm:alpha2` remains true as an axiomatic classification
  under (C1)–(C3), but (C1)–(C3) are not satisfied by H₁.
- The core paper would then need a fresh axiom (or a modified
  H_Γ), and v2 would document the new situation honestly.
- This contingency must be communicated clearly; it is the
  same class of failure as the withdrawn A1–A5 chain, only
  honest about scope.

A decision point at the end of M2 is therefore required.

---

## 9. M1 exit criteria

M1 is closed when:

- [x] All H_Γ definitions in `TGP_v1/` are inventoried.
- [x] The inconsistency H₁ vs H₃ is stated explicitly.
- [x] Variable conventions (ŝ, Φ, φ, g, α(ψ)) are fixed.
- [x] The two algebraic continuum limits (§4.1, §4.2) are
      written out, so the "α=0 vs α=2" confusion is resolved:
      the 0 comes from H₁-leading-order, the 2 comes from H₃.
- [x] The canonical H for OP-6 is chosen (**H₁**).
- [x] M2's three sub-tasks (a/b/c) are stated precisely.

Proceed to M2.
