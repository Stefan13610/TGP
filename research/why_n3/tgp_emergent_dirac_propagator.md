# Emergent Dirac Propagator from RP² Topological Defects in TGP

## Status

This note sketches a possible derivation of the effective Dirac propagator in TGP from the topology of the substrate defect sector.

The claim is not that the full Dirac propagator follows from topology alone. Topology can provide the spinorial character of the excitation, while the propagator itself also requires local dynamics, an effective metric, tetrads, a mass scale, and quantization.

The proposed structure is:

```text
relational substrate Γ
→ real Z₂-odd microscopic amplitude s_i
→ scalar geometric order parameter Φ
→ projective orientational defect sector RP²
→ minimal Q_eff = 1/2 defect
→ spinorial quantum state
→ effective Dirac operator on the TGP metric
→ Dirac propagator as the Green function of that operator
```

---

## 1. Microscopic substrate and coarse-grained geometry

The microscopic substrate is not assumed to be a smooth manifold. It is a relational graph:

```text
Γ = (V, E)
```

with a real microscopic amplitude/operator on the nodes:

```text
s_i ∈ R
```

The microscopic variable is Z₂-odd:

```text
s_i → -s_i
```

The canonical effective scalar field is not `s_i` itself, but a coarse-grained quadratic order parameter:

```text
Φ(x) ∝ <s²>_cg
```

or, in dimensionless form:

```text
ψ = Φ / Φ₀
```

Since Φ is quadratic in `s`, the effective geometric scalar does not retain the sign of the microscopic amplitude. Thus the sign information is physically projected out at the scalar-geometric level.

This naturally suggests that any residual orientational sector associated with the microscopic substrate should be projective:

```text
n ~ -n
```

so that the local orientational variable belongs to:

```text
n ∈ RP² = S² / Z₂
```

This is the key topological input for the fermionic sector.

---

## 2. Effective geometric background

The current canonical TGP metric is the M9.1'' metric:

```text
ds² = -c₀² (4 - 3ψ)/ψ dt² + ψ/(4 - 3ψ) δ_ij dx^i dx^j
```

Define:

```text
A(ψ) = (4 - 3ψ)/ψ
B(ψ) = ψ/(4 - 3ψ)
```

Then:

```text
ds² = -c₀² A(ψ) dt² + B(ψ) δ_ij dx^i dx^j
```

and, importantly:

```text
A(ψ) B(ψ) = 1
```

The Lorentzian working domain is:

```text
0 < ψ < 4/3
```

The local speed of light is:

```text
c_loc = c₀ sqrt(A(ψ))
```

A natural diagonal tetrad is:

```text
e^0_t = c₀ sqrt(A)
e^a_i = sqrt(B) δ^a_i
```

with inverse tetrad:

```text
e_0^t = 1 / (c₀ sqrt(A))
e_a^i = 1 / sqrt(B) δ_a^i
```

Since `B = 1/A`, this gives:

```text
e_a^i = sqrt(A) δ_a^i
```

This will later determine the effective Dirac operator.

---

## 3. Projective defect sector

The scalar field Φ describes the coarse-grained local capacity of the substrate to support geometry. However, the fermionic sector is associated not merely with the scalar amplitude, but with a projective orientational defect:

```text
n : S² → RP²
```

Here, `S²` is a sphere enclosing the localized defect core.

Because:

```text
RP² = S² / Z₂
```

antipodal orientations are physically identified:

```text
n ≡ -n
```

This allows a minimal projective hedgehog-like defect whose effective topological charge can be half of the corresponding lifted charge in the covering space.

A useful way to express this is to introduce a local lift:

```text
ñ : S² → S²
```

with the identification:

```text
ñ ~ -ñ
```

Then the effective projective charge may be written schematically as:

```text
Q_eff = 1/2 deg(ñ)
```

For the minimal non-trivial projective defect:

```text
Q_eff = 1/2
```

This object is not a conventional SU(2) skyrmion and not an ordinary vortex. It is better understood as a kink-hedgehog defect in the RP² sector of the substrate.

---

## 4. From RP² topology to spinorial behavior

The crucial point is that the spinorial character does not need to be fundamental. It may arise from the topology of the configuration space of the defect.

Let:

```text
C_defect
```

denote the configuration space of the minimal RP² defect, modulo smooth deformations that preserve the defect sector.

Because the defect lives in a projective target space, the configuration space can contain a non-trivial Z₂ loop. Physically, this loop corresponds to a 2π rotation of the defect.

The proposed topological structure is:

```text
π₁(C_defect) = Z₂
```

Thus there are two classes of loops:

```text
trivial loop
non-trivial loop
```

The non-trivial loop is associated with a 2π rotation of the defect. Quantization may choose the non-trivial one-dimensional representation of this Z₂ structure:

```text
Ψ[non-trivial loop] = -Ψ
```

Therefore:

```text
Ψ(2π) = -Ψ(0)
Ψ(4π) = +Ψ(0)
```

This is the defining transformation behavior of a spinor.

Thus, in TGP, the spinor is not introduced as a fundamental field on a pre-existing spacetime. Instead, it appears as the quantum state of a topological defect of the substrate.

The logical chain is:

```text
s_i is Z₂-odd
→ scalar geometry sees s_i²
→ residual orientation is projective
→ n ∈ RP²
→ minimal defect has Q_eff = 1/2
→ configuration space has a non-trivial Z₂ loop
→ quantization assigns phase -1 to the 2π loop
→ the defect behaves as a spin-1/2 object
```

---

## 5. Effective collective-coordinate description

At low energy, the localized defect can be described by collective coordinates.

The most important collective coordinate is its position:

```text
X^μ(τ)
```

There may also be internal orientational/topological collective coordinates:

```text
q(τ)
```

The field configuration is then approximated as a family:

```text
n(x; X, q)
Φ(x; X, q)
```

Substituting this family into the effective action gives a defect action of the schematic form:

```text
S_defect = ∫ dτ [
    -M_core(ψ) c_loc²
    + rotational / internal kinetic terms
    + topological Berry phase
]
```

The topological term is responsible for the spinorial sign:

```text
exp(i S_top[2π]) = -1
```

equivalently:

```text
S_top[2π] = π mod 2π
```

This term encodes the non-trivial topology of the projective defect configuration space.

The mass of the effective particle is interpreted as the core energy of the defect expressed in local geometric units:

```text
m_eff(ψ) c_loc² = E_core(ψ)
```

Since:

```text
c_loc² = c₀² A(ψ)
```

one obtains schematically:

```text
m_eff(ψ) = E_core(ψ) / [c₀² A(ψ)]
```

In the simplest low-energy approximation, this can be written as:

```text
m_eff(ψ) = m₀ F(ψ)
```

where `F(ψ)` should ultimately be derived from the defect core energy and its coupling to the local substrate state.

---

## 6. Emergent spinor field

After quantization, the defect state carries spin-1/2. The low-energy quantum state can therefore be represented by a spinor.

For a relativistic effective theory, the natural long-wavelength field is a Dirac spinor:

```text
Ψ(x)
```

This does not mean that Ψ is fundamental. It is an effective field describing the propagation of the quantized RP² defect.

The corresponding effective Dirac equation on the emergent TGP geometry is:

```text
D_TGP Ψ = 0
```

with:

```text
D_TGP = i γ^a e_a^μ (∂_μ + Ω_μ) - m_eff(ψ)
```

where:

```text
e_a^μ     = inverse tetrad of the TGP metric
Ω_μ       = spin connection
m_eff(ψ) = effective defect mass
```

The spin connection is:

```text
Ω_μ = 1/4 ω_μ^{bc} γ_bc
```

At this stage, `Spin(3,1)` enters only emergently: not as a fundamental target space of the microscopic substrate, but as the structure required to describe a spinorial excitation propagating on the emergent Lorentzian geometry.

---

## 7. Dirac operator on the M9.1'' metric

For the M9.1'' metric:

```text
A(ψ) = (4 - 3ψ)/ψ
B(ψ) = ψ/(4 - 3ψ) = 1/A(ψ)
```

The inverse tetrad gives:

```text
e_0^t = 1 / (c₀ sqrt(A))
e_a^i = sqrt(A) δ_a^i
```

Therefore the effective Dirac operator is:

```text
D_TGP =
i γ^0 [1 / (c₀ sqrt(A))] (∂_t + Ω_t)
+
i γ^i sqrt(A) (∂_i + Ω_i)
-
m_eff(ψ)
```

In a local patch where ψ varies slowly, the spin connection can be neglected at leading order:

```text
Ω_μ ≈ 0
```

Then:

```text
D_TGP ≈
i γ^0 [1 / (c₀ sqrt(A))] ∂_t
+
i γ^i sqrt(A) ∂_i
-
m_eff(ψ)
```

This is the local effective Dirac operator for a spinorial defect propagating on the TGP background.

---

## 8. Local propagator

In a locally homogeneous region where ψ is approximately constant, use the plane-wave replacement:

```text
∂_t → -iE
∂_i → ip_i
```

The local momentum-space operator becomes:

```text
D_TGP(p; ψ) =
γ^0 E / [c₀ sqrt(A)]
-
γ^i sqrt(A) p_i
-
m_eff(ψ)
```

The corresponding local propagator is the inverse Green function:

```text
S_TGP(p; ψ) = i D_TGP^{-1}(p; ψ)
```

Explicitly:

```text
S_TGP(p; ψ) =
i [
    γ^0 E / [c₀ sqrt(A)]
    -
    γ^i sqrt(A) p_i
    +
    m_eff(ψ)
]
/
[
    E² / (c₀² A)
    -
    A |p|²
    -
    m_eff²(ψ)
    +
    iε
]
```

where:

```text
A(ψ) = (4 - 3ψ)/ψ
```

This is the proposed local TGP-modified Dirac propagator.

---

## 9. Standard Dirac limit

For the vacuum normalization:

```text
ψ = 1
```

we have:

```text
A(1) = 1
B(1) = 1
```

The metric reduces locally to the standard Minkowski form:

```text
ds² = -c₀² dt² + δ_ij dx^i dx^j
```

The local propagator becomes:

```text
S_TGP(p; 1) =
i [
    γ^0 E / c₀
    -
    γ^i p_i
    +
    m
]
/
[
    E² / c₀²
    -
    |p|²
    -
    m²
    +
    iε
]
```

This is the standard Dirac propagator, up to conventions for units and gamma-matrix signs.

Thus the ordinary Dirac propagator appears as the flat-vacuum limit of the effective propagator of the topological defect.

---

## 10. Interpretation

The proposed interpretation is:

```text
Dirac spinors are not fundamental fields placed on spacetime.
They are effective quantum states of minimal topological defects of the relational substrate.
```

The topological part explains the spinorial character:

```text
RP² defect topology
→ non-trivial Z₂ loop in configuration space
→ minus sign under 2π rotation
→ spin-1/2 behavior
```

The geometric/dynamical part explains the propagator:

```text
TGP metric
→ tetrads
→ spin connection
→ effective Dirac operator
→ Green function / propagator
```

Therefore the full chain is:

```text
Γ, s_i, Z₂
→ Φ and projective orientation sector
→ n ∈ RP²
→ minimal Q_eff = 1/2 defect
→ spinorial quantization
→ effective Dirac equation on M9.1''
→ TGP-modified Dirac propagator
```

---

## 11. What is derived and what is assumed

### Derived / motivated from topology

The following elements are topological or quasi-topological:

```text
n ∈ RP²
Q_eff = 1/2
π₁(C_defect) = Z₂
Ψ(2π) = -Ψ
spin-1/2 character
```

### Requires geometry and dynamics

The following elements are not determined by topology alone:

```text
m_eff(ψ)
the exact defect core profile
the stabilization radius of the defect
the full spin connection terms
the iε prescription
the quantum measure
```

The propagator therefore should not be described as coming from topology alone.

A more precise statement is:

```text
Topology provides the spinorial defect sector.
The TGP metric and local defect dynamics provide the Dirac operator.
The propagator is the Green function of this emergent Dirac operator.
```

---

## 12. Conservative claim

A conservative and technically safer statement is:

```text
In TGP, the Dirac propagator can be interpreted as the low-energy Green function of a quantized RP² topological defect of the relational substrate. The projective topology of the defect sector supplies the spinorial transformation law, while the emergent M9.1'' metric supplies the tetrads and local causal structure required for the effective Dirac operator.
```

---

## 13. Stronger claim

A stronger version, requiring further proof of defect stability and quantization, is:

```text
The fermionic sector of TGP emerges from minimal Q_eff = 1/2 kink-hedgehog defects in the RP² orientational sector of the substrate. Quantization of the non-trivial Z₂ loop in the defect configuration space yields spin-1/2 states. Their long-wavelength propagation on the emergent TGP geometry is governed by a Dirac operator, whose inverse is the effective Dirac propagator.
```

---

## 14. Main missing derivation

The central missing derivation is:

```text
s_i ∈ R with Z₂ symmetry
→ projective orientational sector RP²
→ minimal stable Q_eff = 1/2 defect
→ configuration-space loop corresponding to 2π rotation
→ phase -1 under this loop
```

This is the mathematical core needed to make the argument robust.

Once this is established, the emergence of the Dirac operator is the natural low-energy relativistic description of the resulting spin-1/2 excitation on the TGP metric.

---

## 15. Summary

The proposed result can be summarized as:

```text
Topology does not determine the full Dirac propagator by itself.

However, in TGP, the RP² topology of the projective defect sector can provide the spinorial nature of a minimal substrate defect. Once such a defect is quantized as a spin-1/2 excitation, its long-wavelength propagation on the emergent M9.1'' metric is described by an effective Dirac operator. The Dirac propagator is then the Green function of that operator.
```

In compact form:

```text
RP² topology
→ spinorial defect
→ effective Dirac field
→ TGP Dirac operator
→ Dirac propagator
```

---

## 16. Integracja z R3 mass spectrum (dodane 2026-05-01)

> Dodatek po sesji 2026-05-01: integracja niniejszej propozycji RP²-defect
> z R3 mass spectrum closure dla TGP-canonical α=2.

### 16.1 Stan R3 po sesji 2026-05-01

R3 (`why_n3/`) zamknął mass spectrum dla TGP-canonical α=2 (z `K(φ)=K_geo·φ⁴`):

| Wielkość | Wynik | PDG | Diff |
|----------|-------|-----|------|
| m_μ/m_e | 206.56 | 206.77 | −0.099% ✓ |
| m_τ/m_e | 3474.28 | 3477.23 | −0.085% ✓ |
| m_τ/m_μ | 16.820 | 16.817 | +0.015% ✓ |
| g₀^τ (Koide K=2/3) | 1.755 | input | margin +0.119 do bariery |
| g₀⁴ (φ-drabinka) | 2.840 | — | > g₀_crit = 1.874 (4. zakazana) |

**Mass formula:** `m_obs = c · A_tail^(5−α) = c · A_tail³` dla α=2.

**Empiryczny rozkład:** `m_obs = c · A_tail² · g₀^n(α)` z liniowym
`n(α) = -1.851·α + 7.394` (numerycznie diff < 0.003 dla α∈[0.5, 2.0]).

### 16.2 Identyfikacja: m_eff(ψ) z RP²-defect ↔ m_obs z R3

Sekcja 5 powyżej proponuje `m_eff(ψ) = m₀·F(ψ)` jako mass solitonowego defektu
w niskoenergetycznym przybliżeniu, ale **nie wyprowadza F(ψ)**.

R3 dostarcza explicit empiryczny F(ψ) przez identyfikację:

```text
m_eff(g₀, α) = c · A_tail²(g₀, α) · g₀^n(α)
             = c · A_tail²(g₀, α=2) · g₀^(7.394 - 1.851·α)
```

dla α = 2 (TGP-canonical):

```text
n(2) = 7.394 - 1.851·2 = 3.692
m_eff(g₀, α=2) = c · A_tail²(g₀) · g₀^3.692
```

To jest **konkretna formuła zamykająca** Sekcję 5 dla M9.1'' background.

### 16.3 Identyfikacja: ψ (M9.1'' parameter) ↔ g₀ (R3 soliton parameter)

W M9.1'' metryka zależy od `ψ = Φ/Φ₀`. W R3 soliton ma parametr `g₀ = g(r=0)`
(centralna wartość pola). Czy te dwie zmienne są tożsame?

**Wstępna odpowiedź:** TAK z normalizacja, ale **wymaga formalnej weryfikacji**:

```text
ψ_local(r) ≡ g(r)/Φ₀ ?  (jeśli soliton sam jest "lokalnym ψ-bryłą")
g₀ ≡ ψ(r=0) = wartość pola w centrum solitonu
```

Wtedy soliton w R3 to **lokalna konfiguracja ψ** w przestrzeni z M9.1'' metryką
poza solitonem (gdzie ψ→1). Bariera g₀_crit = 1.874 < 4/3·η dla η > 1.4 — to
jest gęsto blisko **maximum dynamiki Lorentzowskiej** ψ < 4/3 (sekcja 2 powyżej):

```text
0 < ψ_local < 4/3   (Lorentzian domain)
g₀_crit(α=2) = 1.874 > 4/3 ?  TAK
```

Hmm, to **PROBLEM**: g₀_crit > 4/3, czyli soliton wchodzi w region poza
Lorentzian domain M9.1''. Albo:
- (a) Identyfikacja `ψ ↔ g` jest błędna (różne normalizacje, scale factor)
- (b) Bariera g_min → 0 (singularność metryki M9.1'' przy ψ=0) jest
  dokładnie tym samym co bariera R3 (singularność g(r) = 0 w solitonie)
- (c) Soliton explore both Lorentzian (g<4/3) i Euclidean-like (g>4/3) regions

**Open problem.** Wymaga osobnego studium relacji R3 ↔ M9.1''.

### 16.4 Spinorialny defekt ↔ R3 soliton — czy to ten sam obiekt?

Sekcja 3 powyżej definiuje **RP² hedgehog defect** w `n(x)`. R3 ma
**radialny scalar soliton** w `g(r)`. Czy to ten sam obiekt?

**Możliwa identyfikacja:**

```text
n(x) = (∇g)/|∇g|     (gradient direction = orientation field)
```

Ale `∇g` na sferycznie symetrycznym solitonie wskazuje radialnie, czyli
n(x) = x̂ (radial unit vector). To jest **trywialny n=x̂ na S²**, czyli
**hedgehog defect z deg = 1**.

Jeśli RP² zniża to do `Q_eff = 1/2`, dostajemy fermion. **Ale czy R3 soliton
naturally istnieje w sektorze RP²?** To wymaga że pole substratu nie jest
pojedynczym Φ (skalarne), tylko ma additional orientational DOF.

**Otwarte pytanie:** Czy fundamentalny TGP-substrate ma:
- (Hipoteza A) tylko Φ jako amplituda, plus sublattice graph topology
  daje effectivnie n ∈ RP²
- (Hipoteza B) Φ jako amplituda + niezależny n ∈ RP² jako "internal compass"
- (Hipoteza C) zespolona Φ z globalną U(1)/Z₂ ↦ RP² po projekcji

Każda z tych hipotez ma różne konsekwencje dla SU(2)_L (przyszłej derywacji).

### 16.5 Drabina mas i n(α) ↔ wave-function renormalization

Empiryczne odkrycie `n(α) = -1.851·α + 7.394` jest **uderzająco liniowe**.
Liniowość sugeruje analytyczną pochodną:

```text
dn/dα = -1.851
```

To jest **na granicy** -e/√2 = -1.922 lub -2/√(φ²) = -1.236... — żaden ładny
wzór nie pasuje natychmiast. **Linear regression by sufficient large alpha**
może to być artefakt średniej; potrzeba derywacji teoretycznej.

**Hipoteza:** n(α) odzwierciedla **wave-function renormalization** Z(α) w
QFT-podobnym efektywnym opisie:

```text
Z_ψ(α) = m_bare/m_renormalized = g₀^n(α)
```

W standardowym QFT, Z to czynnik "ubierania" gołej masy przez jednopętlowe
poprawki. W TGP, "ubieranie" pochodzi z **wnętrza solitonu** (core
contribution K_core ~ A² · g₀^{2α}-średnia), gdzie α-zależny prefactor
generuje α-zależny n.

Liniowy n(α) sugeruje **pierwszy rząd perturbacji** w odpowiednim parametrze
rozwijanym wokół α=4 (gdzie n(4)=0, czyli mass formula jest tylko A²
bez core dressing).

**Jeśli to się sprawdza analitycznie**, dawałoby to wyrazistą interpretację:

```text
α = kinetic prefactor exponent (z K(φ) = φ^{2α})
α = 4  →  bare scalar mass formula m ~ A² (Hobart-Derrick balance)
α < 4  →  core dressing dodaje renormalizację g₀^n(α)
α = 2 (TGP-canonical) →  n=3.69, m_obs ~ A²·g₀^3.69
```

To zsuwa **R3 mass spectrum** z **emergent fermion propagator** poprzez
identyfikację: każdy soliton (e, μ, τ) jest topologicznym defektem RP²,
którego wave-function renormalization jest determined przez kinetic
prefactor TGP-canonical α=2.

### 16.6 Trzy generacje jako trzy minima w RP²-defect spectrum

R3 N=3 z bariery: trzy g₀_i (e, μ, τ) mieszczą się pod g₀_crit, czwarta nie.

**Identyfikacja w RP² defect picture:**

Trzy solitony z różnymi `g₀` reprezentują **trzy stany związane** spin-1/2
fermionu z różnymi masami `m_eff(g₀)`. Bariera `g₀_crit` to **maximum mass
gap** — powyżej ten fermion staje się niestabilny (decay przez vacuum
tunneling do mniejszego g₀).

To jest **bezpośrednio falsyfikowalne**: jeśli LHC odkryłoby 4. generację
leptonu z masą `m_4` mieszczącą się w extrapolowanej φ-drabince, R3+RP²
falsyfikowałyby. Aktualne LEP `N_ν = 2.984 ± 0.008` jest **zgodne** z N=3
ale tylko dla **lekkich** (< m_Z/2) neutrinów; ciężkie neutrinos byłyby
testem.

### 16.7 Konkretne kolejne kroki

W świetle integracji R3 ↔ RP²-defect, propozycja roadmap dla cyklu
emergent-Dirac:

#### Faza 1: Formalna identyfikacja ψ (M9.1'') ↔ g₀ (R3)
- Sprawdzić relację skalowania (czy g₀=ψ czy g₀=ψ^k z jakimś k?)
- Rozstrzygnąć paradoks g₀_crit > 4/3
- Zamknąć: gdzie soliton "żyje" w przestrzeni M9.1'' — poniżej, na, czy
  ponad ψ=4/3 horizon?

#### Faza 2: Derywacja n(α) z field theory
- Pochodzić z wave-function renormalization Z(α) w odpowiednim limicie
- Sprawdzić: czy n(4) = 0 jest naturalnym anchor (Hobart-Derrick balance)
- Sprawdzić: czy n(α) ma analytyczną postać (np. n = (4-α)·φ_const)

#### Faza 3: RP² defect quantization
- Pełna analiza π₁(C_defect) = Z₂ dla R3-solitonu
- Berry phase calculation w R3 background
- Pokazanie spin-1/2 z Q_eff=1/2

#### Faza 4: Yukawa coupling do scalar Φ
- Pokazać że soliton R3 effectively couplem do quantum spinor `ψ` przez
  Yukawa
- Mass term `m_eff·ψ̄ψ` z emergent solitonu
- Derywacja Z_ψ z core dressing

#### Faza 5: Path integral i propagator
- Pełny S = S_TGP[Φ] + S_ψ[ψ̄, ψ] + S_int
- Obliczenie ⟨ψ̄(x)ψ(0)⟩ w background solitonu
- **Verify:** propagator pokrywa się z eq. powyżej (Sekcja 8)

### 16.8 Status meta i priority

**Status:** PROPOSED research direction, **NIE** zamknięta derywacja.

**Empiryczne anchors gotowe:**
- m_obs = c·A²·g₀^n(α), n(α) = -1.851α + 7.394 (numerycznie zweryfikowane)
- N=3 z bariery (zweryfikowane dla α=1,2)
- Mass ratios PDG-match <0.1% (zweryfikowane dla α=2)

**Teoretyczne brakuje:**
- Formalna identyfikacja ψ ↔ g₀
- Derywacja n(α) z first principles
- Sublattice/RP² topology dla poziomu 0 substrate
- Spinor quantization (π₁ analiza dla R3 defect specifically)

**Priority dla TGP-program:** **HIGH long-term**, **MEDIUM short-term**.
Bez tego TGP nie domyka warstwy 3 (fermiony jako emergentne); ale zamknięcie
warstwy 3 wymaga 5-10 lat pracy.

Krótszerminowo: warstwa 3a (fundamental fermion field na metryce) jest
**prostszą drogą** i zachowuje większość predyktywnej mocy R3 mass spectrum,
**bez** wymagania emergent-Dirac. Decyzja w gestii autora.

---

## 17. Pliki referencyjne (added 2026-05-01)

### W tym folderze (R3)

- `r3_alpha2_full_closure.py` — m_obs = c·A^(5−α) dla α=2, 6/6 PASS
- `r3_p_alpha_analytical.py` — odkrycie n(α) liniowy fit (diff < 0.003)
- `r3_observable_vs_full_mass.py` — separacja m_obs vs M_full
- `CORRECTIONS_2026-05-01.md` — dokument korygujący po audycie

### W TGP_v1

- `core/sek08_formalizm/sek08_formalizm.tex` — warstwa 3c hipoteza
- `TGP_FOUNDATIONS.md` §4.3 — warstwy 3a/3b/3c distinction
- `axioms/substrat/dodatekB_substrat.tex` — poziom 0 graph + Z₂ structure
- `meta/AUDYT_TGP_2026-05-01.md` — audyt sytuacji metryki M9.1''

### Literatura zewnętrzna

- **Skyrme model:** Witten 1983 Nucl. Phys. B 223, 433 (current algebra,
  baryons, quark confinement)
- **Lattice fermions:** Susskind 1977 Phys. Rev. D 16, 3031
- **RP² defects:** Mermin 1979 "The topological theory of defects in ordered
  media" Rev. Mod. Phys. 51, 591
- **Emergent gravity + fermions:** Volovik 2003 "The Universe in a Helium
  Droplet" Oxford
- **Composite fermion topology:** Read & Green 2000 Phys. Rev. B 61, 10267

---

**Autor dodatków:** Sesja 2026-05-01 (po rezolucji R3 + diagnoza n(α) liniowy fit).
**Status końcowy:** dokument zachowany jako **research direction** (PROPOSED),
NIE jako derywacja. Klucz: empiryczne anchors z R3 są gotowe; brakuje
formalnej topologii poziomu 0 dla pełnego closure.
