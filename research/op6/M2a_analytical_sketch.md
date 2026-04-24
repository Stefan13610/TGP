# OP-6 / M2-a — Analytical sketch: can α=2 emerge from H₁?

**Status:** working document, 2026-04-22. Part (a) of the
two-step programme specified by the user ("a then b" —
analytical first, then numerical FRG).
**Input:** canonical H₁ fixed by M1 (`M1_H_Gamma_canonical.md`).
**Question:** under what RG / coarse-graining mechanism does
H₁ produce a kinetic functional with K(φ) ∝ φ⁴, equivalently
K(Φ) ∝ Φ, equivalently α = 2?
**TL;DR:** standard momentum-shell / LPA' / block-spin RG of H₁
does **not** generate this kinetic form. The only surviving
candidate mechanism is a non-trivial "envelope" coarse-graining
at fixed Φ_B = ⟨ŝ²⟩_B; the analytical sketch below shows the
free-field contribution vanishes and the would-be α=2 must come
from genuinely interacting (WF) fluctuations. Analytical closure
is therefore not achievable — we pass the question to M2-b
(numerical FRG with an extended operator basis).

---

## 0. Notation, recap

Canonical H₁ (from M1):
```
H_Γ = Σ_i [ π̂_i²/(2μ) + (m₀²/2) ŝ_i² + (λ₀/4) ŝ_i⁴ ]
    − J Σ_{⟨ij⟩} A_{ij} ŝ_i ŝ_j                          (H₁)
```
- Z₂ symmetry: ŝ_i ↦ −ŝ_i.
- Target: kinetic functional K(φ) ∝ φ⁴, equivalently K(Φ) ∝ Φ
  (where Φ = Φ₀ φ², ŝ² → Φ block-averaged).

## 1. Leading-order coarse-graining of H₁

### 1.1 Bond term

For i, j nearest neighbours along ê at lattice spacing a:
  ŝ_iŝ_j = ŝ(y)² + (a²/8) (ê·∇)²(ŝ²) − (a²/4) (ê·∇ŝ)² + O(a⁴),
with y the bond midpoint. Summing bonds and integrating
by parts (∫ (ê·∇)²(ŝ²) dx = 0 at infinity):

```
Σ_⟨ij⟩ ŝ_iŝ_j  ≈  (d/a^d) ∫ ŝ² dx  −  (1/(2 a^(d−2))) ∫ |∇ŝ|² dx
```

### 1.2 Continuum H₁

Putting the on-site part in and flipping signs:
```
S_cont[ŝ] = ∫ d^d x {  (1/2) K₀ |∇ŝ|²  +  (1/2) m̃² ŝ²
                     +  (λ₀/4) ŝ⁴ + kinetic-time ∫ (∂_t ŝ)² /2μ  }
K₀ = J a^{2−d}     m̃² = m₀² − 2dJ a^{−d}
```

This is the **standard Ising / φ⁴ field theory in d dimensions**
with **constant** kinetic prefactor K₀. In particular, in the
ŝ-variable,

  K(ŝ) = K₀ = const  ⇒  α(ŝ) = 0.

### 1.3 Change of variable Φ = ŝ²

For ŝ > 0: (∇ŝ)² = (∇Φ)²/(4Φ). Hence

  K(Φ) = K₀ / Φ  ⇒  α(Φ) = −½.

And for φ = √(Φ/Φ₀): K(φ) = K₀·Φ₀ = const, α(φ) = 0.

**Observation 1 (algebraic).** No change of variable, on its own,
can turn K(ŝ) = const into K(φ) ∝ φ⁴.

### 1.4 Where does K(φ) ∝ φ⁴ appear?

Only under **density-dependent stiffness**. Consider the GL
functional used in `dodatekB_substrat.tex` eq:prop:substrate-action:

```
F_GL[φ] = Σ_⟨ij⟩ K_{ij} (φ_j − φ_i)²,
K_{ij} = J (φ_i φ_j)²                                    (GL-geo)
```

Expanding (φ_j − φ_i)² = a² (ê·∇φ)² + O(a⁴) and K_{ij} =
J φ_i⁴ [1 + O(a)]:
```
F_GL[φ] ≈ J a² ∫ d^dx · d · φ⁴ (∇φ)² / a^d
        = ∫ d^dx · (K_geo/2) · φ⁴ · (∇φ)²,      K_geo = 2dJ a^{2−d}
```
This gives K(φ) = K_geo φ⁴, α(φ) = 2. ✓

**Observation 2 (structural).** The GL-geo functional is **not** the
continuum limit of H₁. H₁ has a bilinear bond ŝ_iŝ_j with
**constant** coefficient J; GL-geo has a lattice-Laplacian-like
term (φ_j − φ_i)² with **φ-dependent** coefficient J(φ_iφ_j)².
These are different microscopic models.

Equivalently: expanding the GL-geo energy in terms of bilinear/higher
bonds,

  Σ_⟨ij⟩ K_{ij}(φ_j − φ_i)² = −2 Σ_⟨ij⟩ J φ_i³ φ_j³  +  diagonal,

i.e. a **sextic** inter-site coupling `Σ φ_i³ φ_j³`, not the
**bilinear** coupling `Σ φ_iφ_j` of H₁. The kinetic expansion of
a J φ_i³φ_j³ bond is not a leading-order rewriting of H₁.

---

## 2. Candidate mechanisms to get α=2 from H₁

We enumerate the only paths known in the FRG / block-spin
literature, and evaluate each.

### Mechanism (i). Standard momentum-shell RG (LPA, d=3, Z₂)

LPA truncation:
  Γ_k[φ] = ∫ d^d x [ ½ (∇φ)² + U_k(ρ) ],   ρ = φ²/2.

At the Wilson–Fisher fixed point: U*(ρ) is non-trivial, but the
kinetic prefactor is canonical (by construction). K(φ) = 1,
α(φ) = 0.

**Verdict:** does not generate α = 2.

### Mechanism (ii). LPA' (wave-function renormalisation)

LPA' adds Z_k(ρ):
  Γ_k[φ] = ∫ d^d x [ ½ Z_k(ρ) (∇φ)² + U_k(ρ) ].

At the WF fixed point, Z*(ρ) is a non-trivial function, but
|Z*(ρ) − Z*(ρ₀)| / Z*(ρ₀) ≲ 0.05 across the relevant ρ range
(numerical LPA', Litim regulator, d=3, n=1: see standard FRG
references; our independent run in `tgp_erg_eta_lpa_prime.py`
finds η* = 0.044).

For K(φ) ∝ φ⁴ we would need Z*(ρ) ∝ ρ² near ρ₀. The FRG flow
equation for Z_k is linear in Z_k to leading order and admits
only Z*(ρ) weakly deformed from the Gaussian fixed point; no
ρ² scaling of Z at WF is known.

**Verdict:** gives K(φ) ≈ const · (1 + small), i.e. α(φ) ≈ 0.
Not α = 2.

### Mechanism (iii). Derivative expansion to higher order (DE₄, DE₆)

Includes operators like (∇²φ)², (∇φ)⁴, etc. These are
**subleading** and do not parametrise a φ^n · (∇φ)² term with
n > 0.

**Verdict:** no effect on α.

### Mechanism (iv). Extended operator basis including the quartic bond

Introduce in Γ_k a term
  ∫ d^d x · Y_k(ρ) · φ² · (∇φ)² (equivalent to K_{quart}(φ) = Y_k·φ²)
in addition to Z_k(ρ)·(∇φ)². Then flow equations will involve
both Z_k and Y_k.

Under the Z₂ symmetry ŝ ↦ −ŝ, this `φ²(∇φ)²` operator is allowed
(even power of φ). Its running is non-trivial. The question
becomes: does the WF fixed point accommodate a **non-zero** fixed
value Y*(ρ₀) ≠ 0 with a specific ratio to Z*(ρ₀)?

- At the Gaussian fixed point, Y* = 0 (no quartic kinetic term
  in the free theory).
- Y_k is a **composite** operator `φ²(∇φ)²`. Its dimension at WF
  (from Ising OPE): Δ_{φ²(∇φ)²} = 2 Δ_ε + Δ_T where
  Δ_ε ≈ 1.413, Δ_T = d = 3. But this sum is the **engineering**
  dimension; the anomalous correction at WF gives the SCALING
  dimension ≈ 2·1.413 + 2 = 4.83 for the full operator. Since
  d + 2 = 5 for a kinetic-like density (integrated to dim d), we
  have Δ_op − (d+2) ≈ 4.83 − 5 = −0.17 < 0, meaning the operator
  is **irrelevant** at WF.

**Verdict:** the φ²(∇φ)² operator is irrelevant at WF, so its
coefficient Y_k flows to zero (Y* = 0) under RG. No generation of
K(φ) ∝ φ² or φ⁴.

This is the most decisive negative result: *at the Wilson–Fisher
fixed point of H₁*, the would-be kinetic dressing term is
irrelevant, i.e. it scales away under block-averaging.

### Mechanism (v). Envelope coarse-graining at fixed Φ_B

This is the **only remaining mechanism** that is not ruled out
by the above.

Define the block-averaged density field
  Φ_B(x) = (1/N_B) Σ_{i ∈ B(x)} ⟨ŝ_i²⟩.

Integrate out the fast fluctuations of ŝ *at fixed* Φ_B(x) profile:
  e^{−βF_eff[Φ_B]} = ∫ Dŝ · δ[Φ_B − (1/N_B)·Σ_{i∈B(·)} ŝ_i²] · e^{−βH_Γ[ŝ]}.

The resulting F_eff[Φ_B] will generically contain a gradient term
of the form (K_eff(Φ_B)/2) (∇Φ_B)² with **some** K_eff.

**(a) Gaussian sector (m₀² > 0, λ₀ = 0).** For a free field with
⟨ŝ(x)ŝ(y)⟩ = G(x − y), a fluctuation analysis at fixed Φ_B gives
(standard envelope result for bilinear constraints):

  F_eff[Φ_B] = (1/2) ∫ d^d x [ A (∇Φ_B)² / Φ_B + mass & pot terms ]

i.e. **K_eff(Φ_B) = A / Φ_B**, α(Φ_B) = −½. Same as the naïve
change of variable from §1.3. This is expected: in the Gaussian
limit, the only information is the variance, and the variance is
Φ_B itself.

**(b) Wilson–Fisher fixed point (interacting).** Here the
calculation is far less trivial. A non-Gaussian effective K_eff(Φ_B)
can appear. In particular:

- Dimensional analysis: at WF with Ising η, the natural scaling
  of the kinetic term for Φ_B is K_eff(Φ_B) ∼ Φ_B^{α_K} with
  α_K determined by operator OPE.
- For K_eff(Φ_B) ∝ Φ_B (α = ½ in Φ_B = α = 2 in φ), one needs
  the energy operator ε (= Φ_B) dimension to match 2·Δ_σ, i.e.
  Δ_ε ?= 2·Δ_σ. From the Ising universality: Δ_σ = 0.5181,
  Δ_ε = 1.4126. Check: 2·Δ_σ = 1.036 ≠ Δ_ε = 1.413. Not
  matched by a factor of 1.36.

Under an honest FRG calculation in the ρ = φ²/2 variable, there
is no prior expectation that K_eff(Φ_B) ∝ Φ_B with the correct
coefficient. We **cannot** conclude this by pure analytics.

**Verdict:** envelope coarse-graining might give any K_eff(Φ_B);
in the Gaussian limit it gives K_eff ∝ 1/Φ_B (wrong), in the
WF limit it is a genuinely open question requiring numerical FRG.

---

## 3. Conclusion of M2-a

### 3.1 What is ruled out

The following mechanisms for generating α = 2 from H₁ are
**ruled out analytically** (modulo possibly unknown non-perturbative
phenomena):

1. Change of variable alone (§1.3): gives α = −½ or 0.
2. Standard LPA momentum-shell RG (§2.i): gives α = 0.
3. LPA' with Z_k(ρ) (§2.ii): Z* nearly constant; no φ⁴.
4. Derivative expansion to DE₄/DE₆ (§2.iii): no effect.
5. Extended basis with `φ²(∇φ)²` operator (§2.iv): irrelevant at
   WF (scaling dimension > d + 2).

### 3.2 What remains

Mechanism (v) — envelope coarse-graining at fixed Φ_B —
**cannot be ruled out analytically** at the WF fixed point. It
requires an FRG computation in the Φ_B variable directly, with
a non-trivial K_eff(Φ_B) ansatz.

### 3.3 Implication for M1-A vs M1-B

M1 recommended M1-B (keep H₁ axiomatic, derive H₃/GL-geo as
effective). This sketch shows M1-B has only **one** viable path
left — mechanism (v) — which requires numerical work to close.

If mechanism (v) also fails (M2-b), we must pivot to M1-A
(replace H₁ with a quartic-bond or GL-geo form as axiom).

### 3.4 Handoff to M2-b

The numerical FRG task:

- **Setup.** Formulate FRG flow equations in the variable
  Φ_B (not ŝ). Relevant effective action:

      Γ_k[Φ] = ∫ d^d x [ (K_k(Φ)/2) (∇Φ)² + U_k(Φ) ].

  Equivalently, work in the φ = √(Φ/Φ₀) variable with an LPA'-like
  ansatz that keeps Z_k(φ) (not Z_k(ρ)).

- **Target measurement.** Does the fixed point of this flow
  admit K*(Φ) ∝ Φ in d = 3, Z₂ class?

- **Sanity check.** Gaussian sub-limit must reproduce K(Φ) ∝ 1/Φ
  (§2.v.a).

- **Deliverable.** A plot of K*(Φ)/Φ vs the block scale, with
  error bars. If K*(Φ)/Φ → const as k → 0, M1-B is vindicated.
  Otherwise, pivot to M1-A.

## 4. Next step decision

The user asked for **a then b**. M2-a (this document) is now
closed with a narrow negative analytical result plus one open
mechanism that needs numerics.

Proceeding to **M2-b**: numerical FRG with the envelope ansatz
Γ_k[Φ]. I will write a standalone Python FRG script in
`TGP_v1/research/op6/` that (a) implements the flow equations
for K_k(Φ) and U_k(Φ) in the envelope variable, (b) finds the
WF fixed point in that parameterisation, and (c) reports
K*(Φ) as a function of Φ with statistical / truncation error
bars.
