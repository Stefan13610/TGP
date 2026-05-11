---
title: "Phase 1 setup — 1-loop QED effective action S_eff[A_μ; g_eff[{Φ_i}]] na curved background (Birrell-Davies framework adapted)"
date: 2026-05-11
parent: "[[./README.md]]"
type: phase-setup
phase: 1
status: 🟡 setup phase
sub_needs_addressed: [N0.1, N0.2, N0.3]
risks_addressed: [R1, R4]
tags:
  - phase1
  - 1-loop-QED
  - birrell-davies
  - riegert-action
  - effective-action
  - emergent-metric-input
  - S05-preservation
---

# Phase 1 setup

## §0 — Cel Phase 1

Wyprowadzić formalnie **1-loop QED effective action** S_eff[A_μ; g_eff[{Φ_i}]]
na background `g_eff = G[{Φ_i}, σ_ab, Φ̄]` (per
[[../op-emergent-metric-from-interaction-2026-05-09/Phase1_results.md]]) i z niej
**trace anomaly stress-energy tensor** T^μ_μ_EM,1-loop.

Sub-needs: N0.1 (setup) + N0.2 (sympy LOCK β(α)) + N0.3 (Riegert decomposition).

Risks addressed: R1 (M9.1'' contamination guard) + R4 (S05 preservation).

## §1 — Setup: action principles

### §1.1 — Classical EM on g_eff[{Φ_i}] background

W TGP, foton w obecności Φ-konfiguracji widzi metrykę `g_eff^μν[{Φ_i}]`. Klasyczna
akcja EM kanonicznie sprzężona do g_eff:

```
S_EM[A_μ; g_eff] = -1/4 · ∫ d⁴x · √(-g_eff) · g_eff^μα · g_eff^νβ · F_μν · F_αβ
```

z F_μν = ∂_μ A_ν - ∂_ν A_μ (kowariantna pochodna w Abelian gauge group U(1) jest
=∂, bo A jest 1-formą, nie tensorem covariant — Maxwell 2-forma F = dA jest
zawsze antysymetryczna).

**Klasyczny stress-energy tensor** (już z L01 photon_treatment.md §1):
```
T^μ_ν[A] = -F^{μλ} F_{νλ} + (1/4) δ^μ_ν F^{αβ} F_{αβ}
T^μ_μ = 0          (conformal invariance massless gauge w 4D)
```

### §1.2 — Coupling do fermionów (Dirac sektor) na g_eff

Pełna QED akcja:
```
S_QED[A, ψ_D; g_eff] = ∫ d⁴x √(-g_eff) [
   -1/4 · F^μν F_μν
   + ψ̄_D · (i γ^μ_eff D_μ - m) · ψ_D
]
```

z `γ^μ_eff = e^μ_a γ^a` (vierbein-projected gamma matrices), `D_μ = ∂_μ + i e A_μ
+ ω_μ` (gauge + spin connection).

**S05 check (R4 guard):** ten action zawiera *jedno* fundamentalne pole substratu
Φ (poprzez g_eff = G[{Φ_i}]); A_μ i ψ_D są emergent SM fields *na* tym
backgroundzie, NIE oddzielne fundamental fields. 1-loop integration out fermion
loop ψ_D produkuje *effective action dla A_μ*, ale nie dodaje second propagating
field.

### §1.3 — Effective action 1-loop QED (Birrell-Davies §6.1)

Path-integrating out fermion loop:
```
Z_eff[A_μ; g_eff] = ∫ Dψ̄_D Dψ_D · exp(i S_QED)
                 = exp(i S_eff[A_μ; g_eff])
```

z

```
S_eff[A_μ; g_eff] = S_EM_classical[A_μ; g_eff] + S_eff,fermion-loop[A_μ; g_eff]
```

**Fermion 1-loop contribution** (Birrell-Davies §6.1, Schwinger heat kernel):
```
S_eff,fermion-loop = -i Tr ln(i γ^μ_eff D_μ - m)
                  = (i/2) Tr ln(-D̂² + m² + (1/4) R + (i/2) σ^{μν} F_{μν})
```

z `D̂² = D_μ D^μ`, `σ^{μν} = (i/2)[γ^μ, γ^ν]`, R - skalar Ricciego dla g_eff,
F_μν - field strength.

Heat-kernel expansion:
```
Tr ln(...) = -∫₀^∞ ds/s · Tr exp(-s [-D̂² + m² + (R/4) + (i/2)σ^{μν} F_{μν}])
```

z standardowych Schwinger-DeWitt coefficients (Birrell-Davies 6.41-6.44).

## §2 — Trace anomaly: derivation

### §2.1 — Origin: scaling violation z renormalization

Klasyczna akcja S_EM jest *invariant* pod konformalną transformacją:
```
g_μν → Ω²(x) g_μν,    A_μ → A_μ
```

Klasycznie ⇒ T^μ_μ = 0 (conformal invariance).

Po **renormalizacji** (cutoff Λ_UV → ∞ + counterterms), coupling α(μ) running
breaks scale invariance. Dla 1-loop QED:

```
β(α) = μ · dα/dμ = α²/(3π) + O(α³)         [1-loop, m_e ≪ μ]
```

To jest **Capper-Duff-Halpern (1974)** wynik. Sympy LOCK (Phase1_sympy.py) target.

### §2.2 — Trace anomaly formula (Capper-Duff-Halpern + Birrell-Davies §6.3)

W obecności renormalizacji, expectation value stress-energy tensor 1-loop:

```
⟨T^μ_μ_EM⟩_1-loop = (β(α)/(2α)) · F_{αβ} F^{αβ}
                  + a₁ · □ R
                  + a₂ · R²
                  + a₃ · R^{αβ} R_{αβ}
                  + a₄ · R^{αβγδ} R_{αβγδ}
                  + (curvature × F²) cross-terms
                  + non-local Riegert terms
```

W *konwencjach standardowych* (Birrell-Davies 6.135 + Duff 1994 review):

**Dominant term** (independent of curvature):
```
T_anomaly,F² = (β(α)/(2α)) · F² = (α/(3π)/2) · F²
              = (α/(6π)) · F_{μν} F^{μν}
```

z β(α) = α²/(3π), więc β/(2α) = α/(6π) ≈ (1/137)/(18.85) ≈ 3.87·10⁻⁴.

**Wait — let me re-derive carefully.** Classical β-funkcja QED 1-loop:
```
β(α) = (2α²)/(3π) · N_f         [N_f = number of charged fermion species]
```

Dla single Dirac fermion (electron-only), N_f = 1, daje β(α) = α²/(3π) (jak L01
NEEDS §T.1).

Ratio:
```
β(α)/(2α) = α/(6π) ≈ (1/137)/(6π) ≈ 3.87·10⁻⁴
```

L01 ADDENDUM §3.2 cytuje "β(α)/(2α) ≈ 7.7·10⁻⁷" — różnica rzędu factor 500.
**Investigate:**

L01 ADDENDUM Q3 estimate cytuje:
```
β(α)/(2α) ≈ α/(3π) ≈ 7.7·10⁻⁷
```

Hmm — ale `α/(3π) = (1/137)/(9.42) ≈ 7.7·10⁻⁴`, nie 10⁻⁷. Hmm.

Nope: α/(3π) ≈ (7.297·10⁻³)/(9.425) ≈ 7.74·10⁻⁴ — nie 7.7·10⁻⁷. Wygląda jak
*typo* w L01 ADDENDUM (off by 1000 factor — sprawdzić czy to było α²/(3π) zamiast
α/(3π)).

**Let me recompute from formula:**
```
α ≈ 1/137 ≈ 7.297·10⁻³
α² ≈ 5.33·10⁻⁵
α²/(3π) ≈ 5.33·10⁻⁵ / 9.425 ≈ 5.65·10⁻⁶
β(α)/(2α) = (α²/(3π))/(2α) = α/(6π) ≈ 7.297·10⁻³ / 18.85 ≈ 3.87·10⁻⁴
```

Wait — in L01 ADDENDUM §3.2 Q3 (line 270-272 NEEDS):
```
ρ_EM_quantum(x) ≈ -[β(α)/(2α)] · F_μν F^μν / c_0² + O(R²) corrections
β(α)/(2α) = α/(3π) ≈ 7.7·10⁻⁷ (1-loop QED)
```

This formulation says `β(α)/(2α) = α/(3π)`. To znaczy że L01 użyło konwencji:
```
β(α) = (2α²)/(3π)         [wtedy β/(2α) = α/(3π)]
```

A "7.7·10⁻⁷" to numerical typo. Powinno być `α/(3π) = 7.297·10⁻³ / 9.425 = 7.7·10⁻⁴`.

**Conclusion:**

**β(α) = (2α²)/(3π)** is the **canonical convention** for **single Dirac fermion**
1-loop QED β-function (Capper-Duff-Halpern 1974, eq. 4.5):
```
β(α) ≡ μ dα/dμ = (2/(3π)) · α² + O(α³)         (CDH-1974, e=charge → α=e²/(4π))
```

(Note: some textbooks use convention `β(α) = α²/(3π)` z N_f=1 already factored
in; sympy LOCK Phase 1 will document the canonical form explicitly.)

W **canonical form** (CDH 1974 + Duff 1994):
```
β(α)/(2α) = α/(3π)         [single Dirac fermion 1-loop]
```

Numeric:
```
α/(3π) = (1/137.036) / (9.4248) = 7.7397·10⁻⁴
```

**To jest correct value.** L01 ADDENDUM "7.7·10⁻⁷" jest **typo** — nightmare
calc error w ADDENDUM (3 zer mniej).

**Phase 1 sympy LOCK target N0.2:** verify β(α) = (2/(3π))·α² i β(α)/(2α) = α/(3π)
≈ 7.74·10⁻⁴.

**Action item dla Phase 4 closure:** odnotować typo w L01 ADDENDUM §3.2 Q3 i
**propagate correction**. Numerical estimates Phase 3 use *correct* α/(3π) ≈
7.74·10⁻⁴.

### §2.3 — Riegert action (Riegert 1984)

Trace anomaly w 4D ma postać:
```
⟨T^μ_μ⟩ = c_R · E_4 + c_W · W²_{μνρσ} + c_F · F²_{μν}
                                            + (renormalization scheme-dep terms)
```

z:
- E_4 = R²_{μνρσ} - 4 R²_{μν} + R² (Euler density 4D)
- W²_{μνρσ} = R²_{μνρσ} - 2 R²_{μν} + (1/3) R² (Weyl tensor squared)
- F² = F^{μν} F_{μν} (electromagnetic)
- c_R, c_W, c_F coefficients from heat-kernel + β-functions

Dla **QED 1-loop, single Dirac fermion** (Duff 1994, table):
```
c_F = -β(α)/(2α) = -α/(3π)         (sign by signature convention)
c_R = -(N_s + 11 N_f/2 + 62 N_v) / (360 (4π)²)         (here N_f=1 photon)
c_W = +(N_s + (7/2) N_f + 13 N_v) / (120 (4π)²)
```

(N_s=0 scalars, N_f=1 Dirac, N_v=0 vector dla pure-QED z single fermion + photon
loop. Note N_v=1 jeśli photon loop tracked — ale photon jest classical background
in our setup; quantum loops dotyczą tylko fermion loop. Phase 1 sympy LOCK
clarifies konwencja.)

**Riegert non-local action** (Riegert 1984, eq. 5):
```
Γ_anomaly = (1/4) ∫ d⁴x · √(-g_eff) ·
   [ (c_W · W² + c_R · E_4 + c_F · F²) · (1/Δ_4) · (Q_anomaly density) ]
```

z `Δ_4` four-derivative GJMS-like operator (specific to conformal anomaly
non-locality).

**Lokalizowana wersja** (Riegert 1984 §3 + Duff 1994):
```
Γ_anomaly,local = ∫ d⁴x · √(-g_eff) · [
   α₁ · (R - □σ)·F²        +
   α₂ · (R^{μν} F_{μρ} F^ρ_ν)  +
   α₃ · (□F²)              +
   higher-derivative terms
]
```

z auxilary scalar σ wprowadzonym w lokalizacji non-local action.

**Ostrzeżenie strukturalne (R4 guard):** Riegert local form wprowadza auxilary
scalar σ. **MUSIMY zweryfikować że σ NIE jest niezależnym propagującym polem**
w naszej derivacji — w Riegert 1984 σ jest *gauge-fixing parameter* (dilaton
mode of local conformal transformation), nie second fundamental field. W naszym
TGP setup, σ identifikuje się z `lnψ` lub `lnX` (Φ-related), więc S05 preserved
bezwarunkowo.

**Phase 1 sympy LOCK N0.3 target:** verify, że Riegert localization w naszym
g_eff[{Φ_i}] background daje σ ↦ (function of Φ), NIE second fundamental field.

### §2.4 — Stress-energy tensor 1-loop (final form)

Variation Riegert action względem g_eff^μν daje stress-energy tensor 1-loop. Po
usunięciu klasycznego (T^μ_μ = 0) part, **trace** quantum part:

```
⟨T^μ_μ_EM⟩_1-loop = (α/(3π)) · F^{μν} F_{μν}
                  + b₁ · R[g_eff] · F²
                  + b₂ · R^{μν}[g_eff] · F_{μρ} F^ρ_ν
                  + b₃ · □F²
                  + (Riegert non-local terms decomposed via σ ≡ f(Φ))
```

gdzie:
- (α/(3π))·F² jest **dominant pure-photon dim-4 anomaly** — to jest *L01 N1
  target*, *disjoint* od ψ.1.v3 dim-6 EFT.
- b₁, b₂, b₃ ze Schwinger-DeWitt coefficients (CDH 1974 + Birrell-Davies §6.4).
- σ ≡ f(Φ) — **MUSI** być verified w Phase 2 (R4 guard), żeby NIE wprowadzić
  niezależnego propagującego pola.

## §3 — Phase 1 plan + sympy LOCK targets

### §3.1 — Phase 1 sympy targets

Phase1_sympy.py będzie weryfikować:

1. **β-funkcja 1-loop QED:** `β(α) = (2/(3π)) α²` (CDH 1974)
2. **β/(2α) ratio:** `β/(2α) = α/(3π) ≈ 7.74·10⁻⁴` (numerical for α=1/137)
3. **Conformal invariance classical EM:** `T^μ_μ_classical = 0` w 4D
4. **Trace anomaly EM dim-4:** `T^μ_μ_anomaly_F² = (α/(3π))·F²`
5. **Schwinger-DeWitt coefs canonicality:** `b₁, b₂, b₃` formulas (Birrell-Davies
   6.41-6.44)
6. **R1 guard:** verify że nasza derivation NIE używa M9.1'' specific f(ψ) — 
   ansatz {A(ψ), B(ψ), C(ψ)} 3-funkcyjny per emergent-metric Phase 1.
7. **R4 guard (S05 preservation):** verify że Riegert local σ identifies z funkcją
   Φ, NIE jest second fundamental field.
8. **Dimensional analysis check:** [T^μ_μ_EM,1-loop] = [E·L⁻³] = energy density
   ⇒ [ρ_EM_quantum] = [M·L⁻³] (po podziale przez c_0²).

Target: 8/8 sympy PASS.

### §3.2 — Phase 1 deliverables

- [[Phase1_setup.md]] (this file)
- [[Phase1_results.md]] — trace anomaly explicit form, β-function LOCK,
  Riegert decomposition w g_eff[{Φ_i}], R1+R4 guards verified
- [[Phase1_sympy.py]] — sympy script (8 tests)
- [[Phase1_sympy.txt]] — sympy output

## §4 — Risk addressing in Phase 1

### §4.1 — R1 (M9.1'' contamination)

**Strategy:** **explicit** ansatz {A(ψ), B(ψ), C(ψ)} per
[[../op-emergent-metric-from-interaction-2026-05-09/Phase1_results.md]]
"Strukturalna decyzja: ansatz dla g_eff". Nigdy nie wstawiamy specific f_M911(ψ)
= (4-3ψ)/ψ. Phase 1 derivation pozostaje w generalnym ansatzu.

**Verification (sympy test 6):** sprawdzić, że żaden term w T^μ_μ_EM,1-loop nie
zawiera explicit (4-3ψ)/ψ ani factor (4-3ψ).

### §4.2 — R4 (S05 preservation, single-Φ)

**Strategy:** quantum loops integrate out *fermion* (matter sector) na fixed
background g_eff[{Φ_i}]. To NIE wprowadza second fundamental field — Φ pozostaje
jedynym fundamental field substratu, A_μ i ψ_D są emergent SM fields.

**Verification (sympy test 7):** Riegert local σ identifies via g_eff functional
dependence: `σ = -(1/2) ln(det(g_eff)/det(η))` (standard conformal mode
extraction). W naszym Phase 1 ansatz:
```
σ_eff = -(1/2) ln(A·B³ · (1 + σ_ab C/(B·Φ_0² c²))) ≈ -(1/2)[ln A + 3 ln B + (σ-corr)]
       = function of (ψ, ∂Φ)
```

Czyli σ_eff jest *funkcjonałem Φ*, NIE niezależnym fundamental field. S05
preserved.

## §5 — Connection do Phase 2

Phase 1 daje **forma generalna** T^μ_μ_EM,1-loop w obecności g_eff. Phase 2 robi
**TGP-specific reduction**:

1. Substytucja g_eff = G[{Φ_i}, σ_ab, Φ̄] explicit do R, R_μν, etc.
2. Reduction R[g_eff] do funkcji {∂Φ_i + cross-terms ∂Φ_i·∂Φ_j + σ_ab + Φ̄}.
3. Disjointness check: pokaż, że produkt operator structure jest pure-photon dim-4
   + (∂Φ × F²) cross-coupling, NIE dim-6 (∂lnX)·(∂lnX)·F·F z ψ.1.v3 basis.
4. GW170817 c_GW=c_EM dispersion linearization check.

## §6 — Cross-references

- [[./README.md]] §"Centralna hipoteza H1"
- [[./Phase0_balance.md]] §3 NEEDS list, §4 gate criteria
- [[../op-emergent-metric-from-interaction-2026-05-09/Phase1_results.md]] (g_eff ansatz)
- [[../op-L01-rho-stress-energy-bridge-2026-05-04/photon_treatment.md]] §1, §5
- Birrell-Davies, "Quantum Fields in Curved Space" (CUP 1982) ch. 6
- Riegert, Phys. Lett. B 134, 56 (1984)
- Capper, Duff, Halpern, Phys. Rev. D 10, 461 (1974)
- Duff, Class. Quantum Grav. 11, 1387 (1994) — review

---

**Phase 1 setup ready.** Next: Phase1_sympy.py + Phase1_results.md.
