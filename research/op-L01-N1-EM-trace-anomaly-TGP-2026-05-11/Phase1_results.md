---
title: "Phase 1 results — 1-loop QED effective action S_eff[A_μ; g_eff[{Φ_i}]] + trace anomaly + sympy LOCK 8/8"
date: 2026-05-11
parent: "[[./README.md]]"
type: phase-results
phase: 1
status: 🟢 RESOLVED — 8/8 sympy PASS
sub_needs_resolved: [N0.1, N0.2, N0.3]
risks_addressed: [R1, R4, R3-partial]
sympy_script: "[[./Phase1_sympy.py]]"
sympy_output: "[[./Phase1_sympy.txt]]"
predecessor: "[[./Phase1_setup.md]]"
tags:
  - phase1
  - 1-loop-QED-derivation
  - trace-anomaly-explicit
  - capper-duff-halpern
  - riegert-localization
  - S05-preserved
  - L01-typo-correction
---

# Phase 1 results

## §0 — Executive summary

**8/8 sympy PASS.** Phase 1 establishes:

1. **β-funkcja 1-loop QED:** β(α) = (2/(3π))·α² (CDH 1974) — sympy LOCK ✓.
2. **Trace anomaly prefactor:** β(α)/(2α) = α/(3π) ≈ **7.74·10⁻⁴** (NIE 7.7·10⁻⁷;
   typo w L01 ADDENDUM §3.2 propagated to Phase 4 closure).
3. **Trace anomaly explicit form** w obecności g_eff[{Φ_i}]:
   `T^μ_μ_EM,1-loop = (α/(3π))·F² + b₁·R[g_eff]·F² + b₂·R^μν[g_eff]·F_μρ F^ρ_ν + b₃·□F² + Riegert localized terms`.
4. **R1 guard:** generic 3-funkcyjny ansatz {A(ψ), B(ψ), C(ψ)} per emergent-metric
   Phase 1; **NIE** używamy M9.1'' specific f_M911 = (4-3ψ)/ψ. ✓
5. **R4 guard (S05 preservation):** Riegert local conformal mode σ_eff identifies
   z funkcją Φ (via g_eff functional), NIE jest second fundamental field. ✓
6. **R3 partial (GW170817 c_GW=c_EM):** dispersion linearization 1-loop preserves
   ω² = c²k²; pełny test Phase 2/3.
7. **L01 ADDENDUM §3.2 typo identyfikowany** — needs correction propagation in
   Phase 4 closure.

| Check | Result |
|---|---|
| T1: β(α) = (2/(3π))·α² | ✅ PASS |
| T2: β/(2α) = α/(3π) ≈ 7.74·10⁻⁴ (z typo correction) | ✅ PASS |
| T3: T^μ_μ_classical = 0 (conformal invariance) | ✅ PASS |
| T4: T^μ_μ_anomaly = (α/(3π))·F² | ✅ PASS |
| T5: dimensional analysis | ✅ PASS |
| T6: R1 guard (generic ansatz, NIE M9.1'') | ✅ PASS |
| T7: R4 guard (S05, σ_eff = funkcja Φ) | ✅ PASS |
| T8: R3 partial (GW170817 dispersion) | ✅ PASS |
| **TOTAL** | **8/8 PASS** |

## §1 — Effective action 1-loop QED na g_eff[{Φ_i}]

### §1.1 — Pełna QED akcja na curved background

W TGP, A_μ i ψ_D są emergent SM fields propagujące na background `g_eff[{Φ_i}]`
(per [[../op-emergent-metric-from-interaction-2026-05-09/Phase1_results.md]]
Phase 1 ansatz). Klasyczna QED akcja:

```
S_QED[A_μ, ψ_D; g_eff] = ∫ d⁴x · √(-g_eff) · [
   -¼ g_eff^{μα} g_eff^{νβ} F_μν F_αβ
   + ψ̄_D · (i γ^μ_eff D_μ - m_e) · ψ_D
]
```

z `F_μν = ∂_μ A_ν - ∂_ν A_μ`, `D_μ = ∂_μ + i e A_μ + ω_μ` (gauge + spin
connection na g_eff), `γ^μ_eff = e^μ_a γ^a` (vierbein-projected).

**S05 verification (Phase 1 setup §1.2):** Φ jest jedynym fundamental field
substratu; A_μ, ψ_D są emergent. ✓

### §1.2 — Path-integrating out fermion loop

Effective action po integrowaniu out fermion loop (Birrell-Davies §6.1 +
Schwinger heat kernel):

```
S_eff[A_μ; g_eff] = S_EM_classical[A_μ; g_eff] + S_eff,fermion-loop[A_μ; g_eff]
```

z

```
S_eff,fermion-loop[A_μ; g_eff] = -i Tr ln(i γ^μ_eff D_μ - m_e)
                                = (i/2) Tr ln(-D̂² + m_e² + R/4 + (i/2) σ^{μν} F_{μν})
```

gdzie `D̂² = D_μ D^μ`, `σ^{μν} = (i/2)[γ^μ, γ^ν]`, R = R[g_eff] skalar Ricciego.

Heat-kernel expansion (Birrell-Davies 6.41-6.44, Schwinger-DeWitt coefficients):

```
S_eff,fermion-loop ≈ -∫₀^∞ ds/s · Σ_n (s)^(n-2) · a_n[g_eff, F]
```

z koefficientami a_n zawierającymi:
- a_0 = 1 (cosmological-constant divergence; absorbed by renormalization)
- a_1 = R/12 (Einstein-Hilbert UV divergence; absorbed by Newton G renormalization)
- a_2 = standard Gauss-Bonnet + Weyl² + R² + F² + R·F² + ... combinations

**Trace anomaly** wynika z scale-running coefficient α(μ):

```
β(α) = μ · dα/dμ = (2/(3π)) · α² + O(α³)         [single Dirac fermion, 1-loop]
```

(sympy LOCK T1, see Phase1_sympy.py).

### §1.3 — Trace anomaly stress-energy tensor

Variation `S_eff` względem g_eff^μν daje stress-energy tensor. Klasyczna część
T^μ_μ = 0 (conformal invariance), kwantowa część jest *trace anomaly*:

```
⟨T^μ_μ_EM⟩_1-loop = (β(α)/(2α)) · F^{αβ} F_{αβ}
                  + a · □ R[g_eff]
                  + b · R²[g_eff]
                  + c · R^{αβ}[g_eff] R_{αβ}[g_eff]
                  + d · R^{αβγδ}[g_eff] R_{αβγδ}
                  + e · R[g_eff] · F²
                  + f · R^{μν}[g_eff] · F_{μρ} F^ρ_ν
                  + g · □F²
                  + h · Riegert non-local terms
```

W Phase 1 fokusujemy się na **dominant (α/(3π))·F²** term (sympy T4 LOCK) +
**curvature × F² cross-coupling** (b₁·R·F², b₂·R^μν F_μρ F^ρ_ν, b₃·□F²) — to są
te terms które zredukują się do TGP-native ∂Φ + F² coupling w Phase 2.

Curvature-only terms (a, b, c, d) to **conformal anomaly cosmological/quantum-gravity**
contributions; są separate sektorem (closure_2026-04-26 T-Λ + Q2 cycle territory),
NIE direct dla N1.

### §1.4 — Numerical value trace anomaly prefactor

**Sympy LOCK T2:**

```
β(α)/(2α) = α/(3π)         [single Dirac fermion, 1-loop QED]
```

Numerical (α = 1/137.036, fine-structure CODATA):
```
α/(3π) = (7.297 · 10⁻³) / (9.425) = 7.74 · 10⁻⁴
```

**⚠ WAŻNA KOREKTA:** L01 ADDENDUM §3.2 Q3 (line 270 w NEEDS.md) cytuje:

> β(α)/(2α) = α/(3π) ≈ 7.7·10⁻⁷ (1-loop QED)

To jest **typo numerical** — correct value to **7.74·10⁻⁴**, nie 7.7·10⁻⁷
(off by 1000 factor; prawdopodobnie zamieszanie z α²/(3π)·m_e factor lub
dodatkowy α factor z Schwinger limit przeliczenia).

**Konsekwencja dla L01 ADDENDUM Q3 estimate:**

L01 ADDENDUM §3.2 cytuje:
```
ρ_EM_quantum/ρ_NS ~ 10⁻¹² w typowym magnetar
```

z calculation: B²/(2μ_0) ~ 4·10²⁸ J/m³ × 7.7·10⁻⁷ → 3·10⁵ kg/m³, vs ρ_NS ~
4·10¹⁷ kg/m³ → ratio 10⁻¹².

Z **correct** prefactor 7.74·10⁻⁴ (3 OOM więcej):
```
ρ_EM_quantum_corrected ~ 4·10²⁸ J/m³ × 7.74·10⁻⁴ / c_0² 
                        = 4·10²⁸ × 7.74·10⁻⁴ / (9·10¹⁶) kg/m³
                        ≈ 3·10⁸ kg/m³ 
ρ_EM_quantum/ρ_NS ~ 3·10⁸ / 4·10¹⁷ = 7.5·10⁻¹⁰
```

Wciąż **bardzo małe** (10⁻⁹ ÷ 10⁻¹⁰ vs 10⁻¹²), ale 3 OOM większe niż w L01
ADDENDUM. **TT10 magnetar test status preserved** (ρ_EM_quantum nadal
zaniedbywalne dla magnetar surface), ale Phase 3 zrobi pełny estimate z
corrected prefactor.

**Action item Phase 4:** propagate correction do L01 ADDENDUM §3.2 Q3 +
[[../../PREDICTIONS_REGISTRY.md]] entries.

## §2 — Riegert action localization w g_eff[{Φ_i}]

### §2.1 — Riegert non-local form

Riegert (1984) wprowadził local action z auxiliary scalar σ:

```
Γ_anomaly,local = ∫ d⁴x √(-g_eff) · [
   ½ (□σ)² · (Q_anomaly factor)
   - σ · (c_W W² + c_R E_4 + c_F F²)
   + (gauge-fixing terms)
]
```

z `σ` zdefiniowanym przez:
```
σ_eff(x) = -(1/2) ln(det(g_eff(x)) / det(η))         [conformal mode of g_eff]
```

### §2.2 — σ_eff w Phase 1 ansatz

Per emergent-metric Phase 1 ansatz:
```
g_eff^00 = -A(ψ),  g_eff^ij = δ^ij·B(ψ) + σ^ij·C(ψ)/(Φ_0² c²),  g_eff^0i = 0
```

determinant w statycznym limit (dominanting diagonal part):
```
det(g_eff) ≈ -A(ψ) · B³(ψ) · (1 + tr(σ_ab) C/(B Φ_0² c²))
```

(traceless σ_ab daje pełen trace 0, więc cross-correction zaczyna od kwadratowych
terms σ^ij σ_ij; dla naszej Phase 1 leading-order pomijamy te cross-terms i
wracamy do nich w Phase 2).

Zatem:
```
σ_eff(x) ≈ -(1/2) ln(A(ψ) · B³(ψ))
         = -(1/2)·[ln A(ψ) + 3 ln B(ψ)]
         = function of ψ = Φ/Φ_0
         = function of Φ alone
```

**Sympy LOCK T7:** σ_eff zawiera tylko ψ jako free symbol (przez A(ψ), B(ψ));
NIE jest niezależnym fundamental field.

**S05 preservation strukturalne:** Riegert localization wprowadza σ_eff jako
*derived quantity* z Φ-konfiguracji, NIE jako new fundamental degree of freedom.
Quantum 1-loop integration nie zmienia tego — fermion loop dostaje "absorbed"
do effective vertex structure A_μ, σ_eff pozostaje zdeterminowane przez Φ.

### §2.3 — Coefficients c_F, c_R, c_W

Dla 1-loop QED z **single Dirac fermion** + **classical photon background**
(N_s=0 scalars, N_f=1 Dirac fermion, N_v=0 vector loop wewnętrzny — photon jest
classical zewnętrzny):

```
c_F = -β(α)/(2α) = -α/(3π)              [trace anomaly EM coefficient]
c_R = -(11 N_f/2) / (360·(4π)²) · m_e²  [conformal anomaly Type A, Euler density]
c_W = +((7/2) N_f) / (120·(4π)²) · m_e² [conformal anomaly Type B, Weyl²]
```

(Duff 1994 review, table; signs by signature convention -+++.)

**Dla our N1 cycle** dominant term jest c_F coefficient — to daje
(α/(3π))·F² dim-4 anomaly term.

c_R, c_W dotyczą *gravitational* sector trace anomaly; oddzielone od N1
(Q2 vacuum budget territory).

## §3 — Disjointness preview (Phase 2 target)

Operator structure w T^μ_μ_EM,1-loop:

| Operator | Dim | Origin | ψ.1.v3 basis match? |
|---|---|---|---|
| `(α/(3π)) · F²` | 4 | classical β-function 1-loop | **NO** (pure-photon, no ∂lnX) |
| `b₁ · R[g_eff] · F²` | 6 | curvature × F² cross-coupling | partial overlap candidate |
| `b₂ · R^{μν} F_μρ F^ρ_ν` | 6 | tensor-curvature × F² | partial overlap candidate |
| `b₃ · □F²` | 6 | derivative correction | NO (no ∂Φ explicit) |

**Phase 2 target:** w obecności g_eff[{Φ_i}], R[g_eff] i R^{μν}[g_eff] redukują
się do kombinacji `(∂Φ_i, ∂Φ_i·∂Φ_j, σ_ab, Φ̄)`. Question: czy resulting
`R·F²` + `R^{μν} F_{μρ} F^ρ_ν` daje kontrybucję w postaci `(∂lnX)·(∂lnX)·F·F`
(które byłoby ψ.1.v3 sektorem)?

Anticipate: **NIE**, bo:
1. R[g_eff] ~ ∂²ψ + (∂ψ)² (2-derivative scalar curvature) — *struktura inna*
   niż (∂lnX)·(∂lnX) (2 explicit ∂lnX legs).
2. R^{μν}[g_eff] ~ ∂_μ ∂_ν ψ + ∂_μ ψ ∂_ν ψ — same.
3. ψ.1.v3 basis L₅'_a = (∂_μ lnX)(∂_ν lnX)·F^{μρ} F^ν_ρ wymaga *explicite*
   2-tensor (∂lnX) outer product — to jest *różne* od R^{μν} (which contains
   ∂²ψ pieces).

Confirmation Phase 2 explicit reduction.

## §4 — R-guard verification (sympy LOCK)

### §4.1 — R1 guard: M9.1'' contamination — PASS (T6)

W Phase 1 derivation **ani razu** nie podstawiamy specific f_M911(ψ) = (4-3ψ)/ψ.
Generic 3-funkcyjny ansatz {A(ψ), B(ψ), C(ψ)} z emergent-metric Phase 1 jest
*untainted* przez M9.1'' falsification.

**Sympy T6:** A(ψ), B(ψ), C(ψ) deklarowane jako `sympy.Function`, a nie jako
specific algebraic expressions. Ich derivative content w R, R_μν jest *generalny*.

### §4.2 — R4 guard: S05 single-Φ — PASS (T7)

Riegert localization σ_eff *redukuje się* do funkcji ψ = Φ/Φ_0 przez determinant
metryki, NIE jest niezależnym fundamental field.

**Strukturalna konsekwencja:** quantum 1-loop QED na g_eff[{Φ_i}] background
NIE wprowadza second fundamental field. Φ pozostaje single substrate field
(S05 preserved).

### §4.3 — R3 partial: GW170817 dispersion — PASS (T8)

1-loop QED renormalizes `α(μ)` running, ale:
- Photon dispersion na background g_eff[{Φ_i}] po field redefinition A_μ → A_μ/√Z
  pozostaje canoniczna: ω² = c²·k² + O(F²/m_e²)·sub-leading.
- Quantum corrections NIE wprowadzają propagującego scalar mode — są Lagrangian
  modyfikacje (β-function running), nie new dynamical DOF.

Pełny test (anisotropic c, scalar mode amplitude) — Phase 2/3.

## §5 — Phase 1 → Phase 2 handoff

### §5.1 — Co Phase 1 dało

1. **Formal effective action S_eff[A_μ; g_eff]** w postaci Birrell-Davies +
   Riegert localization.
2. **Trace anomaly explicit:** dominant `(α/(3π))·F²` + curvature × F²
   cross-coupling + non-local Riegert (auxiliary σ_eff = funkcja Φ).
3. **Numerical prefactor LOCK:** 7.74·10⁻⁴, NIE 7.7·10⁻⁷ (typo correction).
4. **R1, R4 guards** verified — generic ansatz + S05 preserved.
5. **R3 partial guard** — dispersion linearization preserves ω² = c²k².

### §5.2 — Co Phase 2 musi dostać

1. **Reduction R[g_eff], R^{μν}[g_eff] do TGP-native form** — explicit funkcje
   {∂Φ_i, ∂Φ_i·∂Φ_j, σ_ab, Φ̄}.
2. **Disjointness check vs ψ.1.v3 basis** — pokazać explicit że *resulting*
   operator structure NIE odtwarza dim-6 (∂lnX)·(∂lnX)·F·F.
3. **GW170817 c_GW=c_EM full check** — dispersion analysis w obecności curvature
   × F² mixing terms.
4. **TGP-specific ρ_EM_quantum form** dla {Φ_i} configuration:
   ```
   ρ_EM_quantum[{Φ_i}] = -[1/c_0²] · [ (α/(3π))·F² + (curvature × F²) [reduced] + Riegert [reduced]]
   ```

## §6 — Findings (exportable)

| ID | Finding | Source |
|---|---|---|
| **F1.1** | β-function 1-loop QED `β(α) = (2/(3π))·α²` (single Dirac fermion) | sympy T1 |
| **F1.2** | Trace anomaly prefactor `β(α)/(2α) = α/(3π) ≈ 7.74·10⁻⁴` | sympy T2 + numerical |
| **F1.3** | L01 ADDENDUM §3.2 Q3 cytuje "7.7·10⁻⁷" — **typo**, correct value **7.74·10⁻⁴** (factor 1000 off) | Phase 1 §1.4 |
| **F1.4** | Trace anomaly explicit form: `T^μ_μ_EM,1-loop = (α/(3π))·F² + curvature × F² + Riegert non-local` | Birrell-Davies §6.3 + Phase 1 §1.3 |
| **F1.5** | Riegert localization σ_eff = funkcja Φ — S05 preserved bezwarunkowo | Phase 1 §2.2 + sympy T7 |
| **F1.6** | Generic 3-funkcyjny ansatz {A, B, C} per emergent-metric — R1 guard | Phase 1 §1.1 + sympy T6 |
| **F1.7** | GW170817 dispersion preserved 1-loop, partial check (full Phase 2/3) | sympy T8 |
| **F1.8** | Phase 2 disjointness target: `R·F²`, `R^{μν} F_{μρ} F^ρ_ν` w TGP-reduction NIE produkują `(∂lnX)·(∂lnX)·F·F` ψ.1.v3 form | Phase 1 §3 |

## §7 — Open from Phase 1 (deferred)

- (R3 full): GW170817 dispersion w obecności curvature × F² mixing — Phase 2/3.
- (Cross-cycle): magnetar regime ρ_EM_quantum estimate z **corrected** 7.74·10⁻⁴
  prefactor — Phase 3.
- (PREDICTIONS_REGISTRY): entry M911-EM-quantum z corrected numerical — Phase 4.
- (L01 ADDENDUM correction): typo propagation — Phase 4.

## §8 — Cross-references

- [[./README.md]]
- [[./Phase0_balance.md]] §3 NEEDS list, §4 6/6 gate
- [[./Phase1_setup.md]] §1-§4
- [[./Phase1_sympy.py]] / [[./Phase1_sympy.txt]] (8/8 PASS)
- [[../op-emergent-metric-from-interaction-2026-05-09/Phase1_results.md]] (g_eff ansatz {A, B, C})
- [[../op-L01-rho-stress-energy-bridge-2026-05-04/photon_treatment.md]] §5
- [[../op-L01-rho-stress-energy-bridge-2026-05-04/ADDENDUM_2026-05-10_native_observables_first.md]] §3.2 Q3 (typo source)
- [[../op-psi1-substrate-light-acceleration/ADDENDUM_2026-05-10_native_observables_first.md]] §3 (Q1 disjointness target)
- Birrell, Davies, "Quantum Fields in Curved Space" (CUP 1982) §6.1-6.4
- Riegert, Phys. Lett. B 134, 56 (1984)
- Capper, Duff, Halpern, Phys. Rev. D 10, 461 (1974)
- Duff, Class. Quantum Grav. 11, 1387 (1994) — review

---

**Phase 1 close:** 8/8 sympy PASS. Phase 2 may proceed.
