---
title: "Phase 1 setup — Higgs Lagrangian + SSB + classical T^μ_μ_vac=0 + 1-loop quantum h(x) trace anomaly"
date: 2026-05-11
parent: "[[./README.md]]"
type: phase-setup
phase: 1
status: 🟡 setup phase
sub_needs_addressed: [N0.1, N0.2, N0.3, N0.4, N0.5, N0.6, N0.11]
risks_addressed: [R1, R2, R3-partial, R4]
predecessor: "[[./Phase0_balance.md]] (6/6 gate PASS)"
sister_cycle_architecture: "[[../op-L01-N1-EM-trace-anomaly-TGP-2026-05-11/Phase1_setup.md]] (Birrell-Davies+Riegert) + [[../op-L01-N2-QCD-trace-anomaly-cosmology-2026-05-11/Phase1_setup.md]] (β-function + non-perturbative + crossover)"
tags:
  - phase1
  - Higgs-Lagrangian
  - SSB-cancellation
  - 1-loop-fluctuations
  - beta-lambda
  - gamma-m
  - hierarchy-problem
---

# Phase 1 setup

## §0 — Cel Phase 1

Wyprowadzić formalnie **Higgs trace anomaly** w obecności g_eff[{Φ_i}] background:
1. Higgs Lagrangian + SSB potential
2. Classical T^μ_μ + SSB cancellation: T_vac = 0 (post-renormalization)
3. 1-loop quantum fluctuations h(x) wokół vacuum
4. Trace anomaly form: γ_m·m_H²·h² + (β_λ/4)·h⁴
5. β_λ at 1-loop: self-coupling + Yukawa + gauge contributions
6. γ_m anomalous dimension at 1-loop (Yukawa dominant)
7. Riegert-like decomposition w g_eff[{Φ_i}]
8. R1, R2, R3 partial, R4 guards

Sub-needs: N0.1-N0.6 + N0.11.

## §1 — Setup: Higgs Lagrangian + SSB

### §1.1 — Standard Model Higgs

Higgs doublet H w SU(2) gauge group; po unitary gauge fixing daje single physical
Higgs scalar h:

```
L_H = (∂_μ h)·(∂^μ h)/2 - V(h)
V(h) = -μ²·h²/2 + (λ/4)·h⁴            (SSB form)
```

z μ² > 0, λ > 0 (SSB). Vacuum minimum: dV/dh = -μ²·h + λ·h³ = 0 →
h_vac² = μ²/λ ≡ v² → **v = √(μ²/λ) ≈ 246.22 GeV** (PDG).

**Higgs mass:** m_H² = ∂²V/∂h² |_h=v = -μ² + 3λv² = -μ² + 3μ² = 2μ² → 
**m_H = √(2μ²) = √(2λ)·v** ≈ 125.25 GeV (PDG 2024).

**Self-coupling:** λ = m_H²/(2v²) ≈ (125.25)²/(2·246.22²) ≈ 0.129.

### §1.2 — Stress-energy tensor

Dla scalar field z L = (1/2)(∂h)² - V(h):
```
T^μ_ν = ∂^μ h · ∂_ν h - δ^μ_ν · L
T^μ_μ = (∂h)² - 4·L = (∂h)² - 4·[(1/2)(∂h)² - V(h)]
      = -(∂h)² + 4·V(h)
```

### §1.3 — Classical vacuum: T^μ_μ_vac w SSB

Przy h = v (vacuum minimum), ∂h = 0:
```
T^μ_μ_vac_bare = 4·V(v) = 4·[-μ²v²/2 + λv⁴/4] = -2μ²v² + λv⁴
```

Z relacji v² = μ²/λ → μ² = λv² → 2μ²v² = 2λv⁴, więc:
```
T^μ_μ_vac_bare = -2λv⁴ + λv⁴ = -λv⁴
```

⇒ **Bare vacuum trace** `T^μ_μ_vac_bare = -λv⁴ ≠ 0` (non-zero!).

**KLUCZOWA OBSERWACJA:** Standard SM convention jest *renormalization subtraction*
of vacuum energy:
```
V_renormalized(h) = V(h) - V(v) = V(h) - V_min
```

Po tej subtrakcji, V_renormalized(v) = 0, więc:
```
T^μ_μ_vac_renormalized = -(∂h)² + 4·V_renormalized(h) → 0   przy h=v
```

To jest standard convention SM (tzw. "vacuum energy renormalization"). W kontekście
Q2 F1: bare V(v) = -λv⁴ jest **substrate-decoupled** od bare Λ_TGP (Q2 F1 mechanism);
post-renormalization vacuum energy effectively zero w observed Λ.

**L01 NEEDS §T.4 cytuje:**
> "co dla Higgs SSB: 2m²v² = λv⁴ ⇒ T = 0 w klasycznej próżni"

To jest convention-mixed statement. **Properly stated:** post-renormalization
T^μ_μ_vac = 0; bare T^μ_μ_vac = -λv⁴ jest z Q2 F1 substrate-decoupled. Phase 1
sympy LOCK będzie verified obie wartości.

### §1.4 — 1-loop quantum fluctuations h(x) wokół vacuum

Po SSB, fluctuations h(x) wokół v:
```
h(x) = v + δh(x)         z ⟨δh⟩ = 0
```

Quantum 1-loop effective potential (Coleman-Weinberg 1973 framework):
```
V_eff_1-loop(h) = V_classical(h) + (1/64π²)·∑_i n_i·m_i⁴(h)·[ln(m_i²(h)/μ_R²) - C_i]
```

z:
- n_i = degrees of freedom (boson, fermion, gauge)
- m_i²(h) = field-dependent masses (Higgs self-mass, Yukawa fermion masses, gauge boson masses)
- μ_R = renormalization scale (MS-bar)
- C_i = scheme constants (3/2 for scalar/fermion, 5/6 dla gauge)

**Field-dependent masses:**
- Higgs self: m_H_eff²(h) = -μ² + 3λh² → przy h=v: m_H² = 2μ² = 2λv² ✓
- Top Yukawa: m_t_eff(h) = y_t·h/√2 → przy h=v: m_t = y_t·v/√2 ≈ 173 GeV ✓
- W boson: m_W_eff(h) = g·h/2 → przy h=v: m_W ≈ 80.4 GeV
- Z boson: m_Z_eff(h) = √(g²+g'²)·h/2 → przy h=v: m_Z ≈ 91.2 GeV

### §1.5 — Trace anomaly explicit (Coleman-Weinberg formalism)

Trace anomaly w 1-loop QFT (analog do N1 §2):
```
T^μ_μ_quantum = (∑_i β_i / coupling_i) · [field-strength terms]
              + (γ_m_field) · [mass terms]
              + curvature × field² mixing + Riegert non-local
```

Dla Higgs sektor explicit (per CW-1973 + Sirlin-1980):
```
T^μ_μ_Higgs_1-loop ≈ (1+γ_m)·m_H²·h² + (β_λ/4)·h⁴·c_λ + curvature × h² + Riegert
```

z:
- γ_m = anomalous Higgs mass dimension (z Yukawa loops, dominant top)
- β_λ = Higgs self-coupling running
- c_λ = scheme-dependent coefficient

### §1.6 — β_λ at 1-loop (Higgs self-coupling running)

Standard SM 1-loop β_λ (Wetterich 1981, Sher 1989):
```
β_λ = (1/16π²) · [
   24λ²                                                 (Higgs self)
   - 6·y_t⁴                                             (top Yukawa, dominant negative)
   + (3/8)·(2g⁴ + (g² + g'²)²)                           (gauge)
   + 12λ·y_t² - 9λ·g² - 3λ·g'²                          (mixing terms)
   + smaller Yukawa contributions
]
```

**Numerical at electroweak scale (MS-bar, μ_R = m_t):**
```
β_λ ≈ (1/16π²) · [24·(0.129)² - 6·(0.99)⁴ + (3/8)·(0.66⁴ + (0.66² + 0.36²)²) + ...]
     ≈ (1/16π²) · [0.40 - 5.76 + 0.085 + ...]
     ≈ -(5.27)/(16π²)
     ≈ -3.3·10⁻²
```

⇒ **β_λ < 0 at electroweak scale** — Higgs self-coupling **decreases** z increasing
energy (top Yukawa dominant negative contribution). To prowadzi do **stability bound**
(λ → 0 lub negative przy Λ_instability ~ 10⁹⁻¹⁰ GeV per CW analysis).

### §1.7 — γ_m at 1-loop (Higgs mass anomalous dimension)

Standard SM 1-loop γ_m (Higgs mass):
```
γ_m = (1/16π²) · [12λ - 6·y_t² + (gauge corrections)]
```

**Numerical at EW scale:**
```
γ_m ≈ (1/16π²) · [12·0.129 - 6·(0.99)² + ...]
     ≈ (1/16π²) · [1.55 - 5.88 + ...]
     ≈ -(4.3)/(16π²)
     ≈ -2.7·10⁻²
```

⇒ **γ_m < 0** — Higgs mass running negative (top Yukawa dominant).

### §1.8 — Riegert localization w g_eff[{Φ_i}]

Standard Riegert decomposition (analog do N1 §2.3 + N2 §5.1):
```
S_anomaly_Higgs = ∫ d⁴x √(-g_eff) · [
   c_W·W²[g_eff] + c_R·E_4[g_eff]               (curvature-only)
   + γ_m·m_H²·h²                                 (mass-anomalous-dim term)
   + (β_λ/4)·h⁴                                  (self-coupling running)
   + d_1·R[g_eff]·h² + d_2·R^μν·∂_μh·∂_νh        (curvature × field)
   + Riegert non-local with σ_eff = function(ψ)
]
```

z σ_eff = function(ψ) (Riegert auxiliary identifies z Φ — S05 preserved).

## §2 — Risk addressing in Phase 1

### §2.1 — R1 (M9.1'' contamination) — same as N1+N2

Generic ansatz {A(ψ), B(ψ), C(ψ)} per emergent-metric Phase 1; nigdy NIE
M9.1'' specific f(ψ).

### §2.2 — R2 (renormalization scheme) — honestly documented

Tego cyklu derivation use **MS-bar scheme** (modified minimal subtraction):
- Standard SM convention dla 1-loop calculations
- Coleman-Weinberg + Sirlin framework
- Numerical β_λ, γ_m at μ_R = m_t (top mass scale, common reference)
- Document explicit że **scheme choice affects numerical values** (NIE structural form)

### §2.3 — R3 (hierarchy problem) — partial Phase 1

Quadratic divergence δm_H² ~ Λ_UV² jest standard hierarchy problem w SM. W TGP
framework:

**Hipoteza H4 (R3 partial):** Q2 F1 + S05 mechanism może *strukturalnie* protect
m_H przed quadratic Λ_UV² destabilization, analogous do how Q2 F1 protects bare
Λ_TGP od matter vacuum catastrophe (10⁷⁷ OOM).

**Strategy Phase 1:** verify że Q2 F1 mechanism (substrate-vacuum identification)
**zapewnia** strukturalną renormalization condition która kasuje hierarchy.

**Phase 2 follow-up:** explicit verification w cosmological context (EW phase
transition).

**HONEST CAVEAT:** tego cyklu NIE *rozwiąże* hierarchy problem fully (to byłaby
revolution); może *jedynie pokazać consistency* z Q2 F1 mechanism. R3 honestly
documented jako *deferred precision* item.

### §2.4 — R4 (S05) — same as N1+N2

h(x) jest emergent SM scalar field na background g_eff[{Φ_i}]; quantum loops
integrate out fermion (Yukawa) + gauge bosons (W, Z) + Higgs self-loops; **NIE
wprowadzają** second fundamental field.

Riegert σ_eff = function(ψ) — single-Φ axiom preserved.

## §3 — Phase 1 sympy LOCK targets (8 tests)

Phase1_sympy.py weryfikuje:

1. **SSB minimum:** v² = μ²/λ; m_H² = 2μ² = 2λv² (analytic LOCK)
2. **Bare vacuum trace:** T^μ_μ_vac_bare = -λv⁴ ≠ 0
3. **Renormalized vacuum trace:** post-subtraction T_vac_renorm = 0
4. **PDG 2024 numerical:** m_H = 125.25 GeV, v = 246.22 GeV → λ ≈ 0.129
5. **β_λ at 1-loop EW scale:** β_λ ≈ -3.3·10⁻² (top Yukawa dominant negative)
6. **γ_m at 1-loop EW scale:** γ_m ≈ -2.7·10⁻² (top Yukawa dominant)
7. **R1 + R4 guards** (generic ansatz {A, B, C}; S05 σ_eff = function(ψ))
8. **R3 partial** + R2 honest documentation

## §4 — Connection do Phase 2

Phase 1 daje **renormalized Higgs trace anomaly form**.

Phase 2 (multi-session continuation) dostanie:
1. EW phase transition T_EW ~ 100 GeV — lattice consensus crossover (NOT first-order)
   dla m_H=125 GeV (per arXiv:2405.01191 2024)
2. Friedmann z transient ρ_Higgs(T) source w EW epoch
3. Reduction ρ_Higgs(T<<v) → 0 strukturalnie (Q2 F1 verification)

## §5 — Cross-references

- [[./README.md]]
- [[./Phase0_balance.md]]
- [[../op-L01-N1-EM-trace-anomaly-TGP-2026-05-11/Phase1_setup.md]] (sister architecture)
- [[../op-L01-N2-QCD-trace-anomaly-cosmology-2026-05-11/Phase1_setup.md]] (sister architecture)
- Coleman, Weinberg, Phys. Rev. D 7, 1888 (1973) — effective potential 1-loop
- Sirlin, Phys. Rev. D 22, 971 (1980) — SM radiative corrections
- Buchbinder, Odintsov, Shapiro 1992 — 1-loop scalar curved background
- PDG 2024 — m_H, v, λ, y_t, gauge couplings

---

**Phase 1 setup ready.** Next: Phase1_sympy.py + Phase1_results.md.
