---
title: "Phase 1 setup — SU(3) Yang-Mills + 6 quark flavors action on g_eff[{Φ_i}]; non-perturbative QCD trace anomaly derivation framework"
date: 2026-05-11
parent: "[[./README.md]]"
type: phase-setup
phase: 1
status: 🟡 setup phase
sub_needs_addressed: [N0.1, N0.2, N0.3, N0.4, N0.11]
risks_addressed: [R1, R2-honest, R4]
predecessor: "[[./Phase0_balance.md]] (6/6 gate PASS)"
sister_cycle_architecture: "[[../op-L01-N1-EM-trace-anomaly-TGP-2026-05-11/Phase1_setup.md]] (Birrell-Davies+Riegert framework, adapted for non-Abelian)"
tags:
  - phase1
  - QCD-non-perturbative
  - SU3-Yang-Mills
  - asymptotic-freedom
  - gluon-condensate
  - SVZ-sum-rules
  - dimensional-transmutation
---

# Phase 1 setup

## §0 — Cel Phase 1

Wyprowadzić formalnie **QCD trace anomaly** w obecności g_eff[{Φ_i}] background:
1. SU(3) YM + 6 quark flavors action.
2. β_QCD(g) 1-loop derivation z asymptotic freedom (b_0 = 7 dla N_c=3, N_f=6).
3. Trace anomaly explicit: T^μ_μ_QCD,1-loop = (β(g)/(2g))·Tr(G_μν G^μν).
4. Non-perturbative regime: gluon condensate ⟨α_s G²/π⟩ ≈ 0.012 GeV⁴ (SVZ + lattice
   external input).
5. Dimensional transmutation: Λ_QCD = 217 MeV (PDG 2024).
6. Riegert-like localization w g_eff[{Φ_i}] z σ_eff = function(ψ) — analog do N1.

Sub-needs: N0.1, N0.2, N0.3, N0.4, N0.11.
Risks addressed: R1 (M9.1'' contamination guard), R2 (non-perturbative regime
honestly documented), R4 (S05 single-Φ preservation).

## §1 — Setup: action principles

### §1.1 — Classical SU(3) Yang-Mills + 6 quarks on g_eff[{Φ_i}]

W TGP, kwarki + gluon na background `g_eff[{Φ_i}]`. Klasyczna QCD akcja
(non-Abelian extension N1 cycle §1.1):

```
S_QCD[A^a_μ, q_f, q̄_f; g_eff] = ∫ d⁴x · √(-g_eff) · [
   -¼ g_eff^{μα} g_eff^{νβ} G^a_μν G^a_αβ
   + Σ_{f=1}^{6} q̄_f · (i γ^μ_eff D_μ - m_f) · q_f
]
```

z:
- `G^a_μν = ∂_μ A^a_ν - ∂_ν A^a_μ + g f^{abc} A^b_μ A^c_ν` (non-Abelian field strength)
- `f^{abc}` SU(3) structure constants
- `D_μ q = (∂_μ + i g A^a_μ T^a) q`, T^a = λ^a/2 fundamental SU(3) generators
- a, b, c = 1, ..., 8 (gluon color indices)
- `f` = u, d, s, c, b, t quark flavors (PDG 2024 masses):
  - m_u ≈ 2.16 MeV, m_d ≈ 4.67 MeV, m_s ≈ 93.4 MeV (light, ≪ Λ_QCD)
  - m_c ≈ 1.27 GeV, m_b ≈ 4.18 GeV, m_t ≈ 172.5 GeV (heavy, > Λ_QCD)

**Differences vs N1 (QED):**

| Aspect | QED (N1) | QCD (N2) |
|---|---|---|
| Gauge group | U(1) Abelian | SU(3) non-Abelian |
| Field strength | F = dA | G = dA + g[A,A] (non-linear) |
| β-function sign | **+** (charge running UP, Landau pole UV) | **−** (asymptotic freedom; running UP IR) |
| Coupling structure | linear photon-fermion | self-coupling gluon-gluon + quark-gluon |
| Confinement | None (free fields) | **YES** (no isolated quarks at low T) |
| Vacuum condensate | None classical | **⟨α_s G²/π⟩, ⟨q̄q⟩** non-perturbative |

**S05 check (R4 guard, analog do N1 §1.2):** A^a_μ, q_f są emergent SM fields *na*
fixed background g_eff[{Φ_i}]; Φ pozostaje single fundamental field. 1-loop
integration produce effective action dla A^a_μ, ale NIE wprowadza second
fundamental field.

### §1.2 — 1-loop effective action via path-integration out fermions

Analogous do N1 §1.3 (Birrell-Davies §6.1 + Schwinger heat kernel). Path-integrating
out quark loops:

```
Z_eff[A^a_μ; g_eff] = ∫ Dq̄ Dq · exp(i S_QCD)
                   = exp(i S_eff^{f}[A^a_μ; g_eff])
```

Plus integrating out gluon loops (ghost loops via Faddeev-Popov + gluon
self-interactions for non-Abelian):

```
S_eff[A^a_μ; g_eff] = S_YM_classical + S_eff^{f-loop} + S_eff^{gluon-loop} + S_eff^{ghost-loop}
```

Dla 1-loop renormalizacji β-funkcji wystarczą:
- N_f Dirac fermion contributions (analog do N1 z N_f=1; tutaj efektywny N_f w
  zależności od skali — dla T ≫ m_t wszystkie 6 active, dla T ≪ m_b efektywnie 4
  light, etc.)
- gluon self-energy (gluon-gluon-gluon + 4-gluon vertex)
- ghost-gluon contribution

### §1.3 — β-function 1-loop QCD (canonical derivation)

Standard result (Gross-Wilczek 1973, Politzer 1973 — Nobel 2004):

```
β_QCD(g) = μ · dg/dμ = -(b_0/(16π²)) · g³ + O(g⁵)
```

z **b_0 coefficient** (1-loop, Yang-Mills + N_f flavors):

```
b_0 = (11/3) N_c - (2/3) N_f
```

Dla **QCD physical** (N_c=3, N_f=6 high-T regime):
```
b_0 = (11/3)(3) - (2/3)(6) = 11 - 4 = 7
```

**Asymptotic freedom verification:** b_0 > 0 ⇒ β_QCD < 0 ⇒ coupling g(μ) **maleje**
gdy μ rośnie (UV-free). To jest *opposite* do QED (β_QED > 0, charge maleje w IR,
Landau pole w UV).

**For α_s = g²/(4π) convention:**
```
β(α_s) = dα_s/dlnμ = 2g·dg/dlnμ / (4π) = -(b_0/(8π²))·g⁴/(4π) = -(b_0/(2π))·α_s²
```

### §1.4 — Non-perturbative regime + Λ_QCD dimensional transmutation

QCD coupling g(μ) **diverges** w IR przy μ → Λ_QCD przez integral 1-loop running:

```
1/g²(μ) = (b_0/(8π²)) · ln(μ/Λ_QCD)
```

⇒ `g(μ) → ∞` for μ → Λ_QCD = scale przy którym przegląd perturbacyjny załamuje się.

**Dimensional transmutation:** Λ_QCD jest *strong-coupling scale* — emerges z
dimensionless coupling g(μ_ref) + reference scale μ_ref:

```
Λ_QCD = μ_ref · exp(-(8π²)/(b_0 · g²(μ_ref)))
       ≈ 217 MeV  [PDG 2024, MS-bar scheme, N_f = 5 flavors near m_b]
```

**Below Λ_QCD:** non-perturbative regime — confinement, chiral symmetry breaking,
gluon condensate.

**Gluon condensate** (Shifman-Vainshtein-Zakharov 1979 sum rules):
```
⟨α_s/π · G^a_μν G^aμν⟩_0 ≈ 0.012 GeV⁴   [SVZ-1979 + lattice 2018+]
```

Z modern lattice QCD updates this can vary in range 0.005-0.02 GeV⁴ (renormalization
scheme dependence). **Tego cyklu approach: traktujemy ⟨α_s G²/π⟩ jako external lattice
input z documented uncertainty band.**

### §1.5 — Trace anomaly explicit form

Klasyczna SU(3) YM jest *conformally invariant* w 4D (analog do EM):
```
T^μ_ν^{cl} = G^a^{μλ} G^a_νλ - ¼ δ^μ_ν G^a_αβ G^aαβ + (quark stress-energy)
T^μ_μ^{cl}|_pure-gauge = 0    (conformal invariance)
T^μ_μ^{cl}|_quark = m_f q̄_f q_f    (mass term, vanishes for massless)
```

**Quantum trace anomaly** (Collins-Duncan-Joglekar 1977, Crewther 1972):
```
T^μ_μ^{QCD,quantum} = (β(g)/(2g)) · G^a_μν G^aμν + Σ_f (1+γ_m) m_f q̄_f q_f
```

**Pure-gauge anomaly** (analog do EM trace anomaly, dim-4 pure-gluonic operator):
```
T^μ_μ_anomaly,gluonic = (β(g)/(2g)) · G²
                      = -(b_0/(16π²)) · g²/2 · G²    [1-loop substitution]
                      = -(b_0/(32π²)) · g² · G²
                      = -(b_0/(8π)) · α_s · G²
                      ≈ -(7/(8π)) · ⟨α_s G²⟩         [N_f=6 high-T]
                      ≈ -(7/8) · ⟨α_s G²/π⟩          [redistribution]
```

Numerical (z SVZ vacuum):
```
T^μ_μ_anomaly_vacuum ≈ -(7/8) · 0.012 GeV⁴ ≈ -0.0105 GeV⁴   [low-energy QCD vacuum]
```

W energy units to jest przeniesione na ρ:
```
ρ_QCD = -⟨T^μ_μ⟩/c_0² ≈ 0.0105 GeV⁴ / c_0² 
                       ≈ 1.87 · 10¹⁷ kg/m³           [non-relativistic equivalent]
```

To jest **gluon condensate gravitational density** — porównywalna z neutron-star
density! Ale to jest w *quantum vacuum*, NIE w macroscopic medium. Per Q2 cycle F1,
ten condensate jest *substrate-decoupled* od bare Λ — *transient phase-transition source*
w wczesnym Wszechświecie (Phase 2 cosmology integration).

### §1.6 — Riegert-like localization w g_eff[{Φ_i}]

Analogous do N1 §2 (Riegert 1984 framework). W obecności curved background:

```
Γ_anomaly_QCD = ∫ d⁴x · √(-g_eff) · [
   c_W · W² + c_R · E_4 + c_G · G²    (gauge sector trace anomaly)
   + curvature × G² mixing terms     (analog R·F²)
   + Riegert non-local with σ_eff = function(ψ)
]
```

z:
- c_W, c_R z curvature-only Type A/B trace anomaly (Duff 1994 — for SU(3)+6 quarks)
- c_G = -β(g)/(2g) = (b_0/(32π²))·g² = (b_0 α_s)/(8π) (gauge trace anomaly coefficient)

**S05 verification (R4 guard, analog do N1 §2.2):** Riegert local σ_eff =
-(1/2) ln(det g_eff) reduces do funkcji ψ przez Phase 1 ansatz {A(ψ), B(ψ), C(ψ)}.
Single-Φ axiom preserved.

## §2 — Phase 1 plan + sympy LOCK targets

### §2.1 — Phase 1 sympy targets (analog do N1; 8 tests)

Phase1_sympy.py będzie weryfikować:

1. **β-funkcja 1-loop QCD:** `β_QCD(g) = -(b_0/(16π²))·g³` z `b_0 = (11/3)N_c - (2/3)N_f`
   for N_c=3, N_f=6 → b_0 = 7
2. **Asymptotic freedom sign:** β_QCD < 0 (vs β_QED > 0) — **opposite to QED**
3. **β(α_s) convention:** `β(α_s) = -(b_0/(2π))·α_s²`
4. **Trace anomaly form:** `T^μ_μ_QCD = (β(g)/(2g))·G²` at 1-loop substitution
5. **Λ_QCD dimensional transmutation:** consistency `Λ_QCD = μ·exp(-1/(2β_0·α_s(μ)))`
   z PDG 2024 value ~ 217 MeV
6. **Gluon condensate dimensional:** `[⟨α_s G²/π⟩] = [GeV⁴]`; SVZ value ~ 0.012 GeV⁴
   ⇒ ρ_QCD_vacuum ≈ 0.0105 GeV⁴ ≈ 1.87·10¹⁷ kg/m³ equivalent
7. **R1 guard:** generic ansatz {A, B, C} (NIE M9.1''); R4 (S05 preservation) —
   gluon condensate jest composite operator z YM stress-energy + Riegert σ_eff =
   function(ψ)
8. **Non-perturbative regime documentation (R2 honest):** verify że tego cyklu
   nasza derivation jest *consistent z* lattice QCD external inputs (Λ_QCD, ⟨α_s G²/π⟩,
   T_c=156 MeV) ale *NIE wyprowadza* ich z first principles — honest accounting

Target: 8/8 sympy PASS.

### §2.2 — Phase 1 deliverables

- [[Phase1_setup.md]] (this file)
- [[Phase1_results.md]] — full derivation + Riegert decomposition + R-guards
- [[Phase1_sympy.py]] + [[Phase1_sympy.txt]] — 8 tests

## §3 — Risk addressing in Phase 1

### §3.1 — R1 (M9.1'' contamination)

**Strategy:** **identical do N1** — explicit ansatz {A(ψ), B(ψ), C(ψ)} per
emergent-metric Phase 1. Nigdy nie wstawiamy specific f_M911(ψ) = (4-3ψ)/ψ.

### §3.2 — R2 (non-perturbative regime — honest documentation)

**Strategy:** lattice QCD + SVZ sum rules są *external*. Tego cyklu nasza derivation
*uses* te wartości jako input, ale NIE wyprowadza ich z first principles.

**Honest documentation w Phase 1 results:**
- β-function 1-loop *jest perturbative*, valid dla μ ≫ Λ_QCD.
- Gluon condensate ⟨α_s G²/π⟩ jest *empirical/lattice/SVZ-input*; ten cykl traktuje
  ją jako parameter z documented uncertainty band 0.005-0.02 GeV⁴.
- Dimensional analysis ρ_QCD ~ Λ_QCD⁴ jest *consistent* z lattice EoS thermal profiles
  (Phase 2).

**NIE blokuje cycle SUCCESS** — analogous do N1 Wilson γ_i numerical pinning
(deferred precision, NIE new free param).

### §3.3 — R4 (S05 single-Φ preservation)

**Strategy:** identical do N1 §4.2.
- Quantum loops integrate out *fermions* (quarks) i *gluons* na fixed background.
  NIE wprowadza second fundamental field.
- Gluon condensate ⟨α_s G²/π⟩ jest *composite operator* z YM stress-energy tensor —
  funkcjonalna od A^a_μ na backgroundzie g_eff[{Φ_i}], NIE niezależny field.
- Riegert local σ_eff = -(1/2) ln(det g_eff) = function(ψ) — analog do N1.

## §4 — Connection do Phase 2

Phase 1 daje **renormalized T^μ_μ_QCD form** + gluon condensate vacuum value.

Phase 2 (multi-session continuation) dostanie **thermal field theory ρ_QCD(T)**:
- HotQCD lattice EoS interaction measure profile (peak near T_c ~ 156 MeV)
- Friedmann equation modyfikacja w QCD epoce (z~10¹²)
- Reduction ρ_QCD(T<<Λ_QCD) → 0 strukturalnie (Q2 F1 konstruktywna verification)

## §5 — Cross-references

- [[./README.md]] §"Centralna hipoteza H1"
- [[./Phase0_balance.md]] §3 NEEDS, §4 6/6 gate
- [[./NEEDS.md]] N0.1-N0.4, N0.11
- [[../op-L01-N1-EM-trace-anomaly-TGP-2026-05-11/Phase1_setup.md]] (sister architecture)
- [[../op-emergent-metric-from-interaction-2026-05-09/Phase1_results.md]] (g_eff ansatz)
- [[../op-Q2-vacuum-budget-2026-05-10/Phase_FINAL_close.md]] §2.3 (transient sources)
- Collins, Duncan, Joglekar, Phys. Rev. D 16, 438 (1977) — canonical CDJ
- Shifman, Vainshtein, Zakharov, Nucl. Phys. B 147, 385 (1979) — SVZ sum rules
- Gross, Wilczek; Politzer (1973) — asymptotic freedom Nobel-Prize-caliber
- Birrell, Davies (CUP 1982) ch. 6 — QFT on curved background
- Riegert, Phys. Lett. B 134, 56 (1984) — conformal anomaly action

---

**Phase 1 setup ready.** Next: Phase1_sympy.py + Phase1_results.md.
