---
title: "Phase 2 setup — thermal field theory ρ_QCD(T)/T⁴ profile + Friedmann equation modyfikacja w QCD epoce + Q2 F1 konstruktywna verification"
date: 2026-05-11
parent: "[[./README.md]]"
type: phase-setup
phase: 2
status: 🟡 setup phase
sub_needs_addressed: [N0.5, N0.6, N0.7]
risks_addressed: [R3-partial, R5, R7]
predecessor: "[[./Phase1_results.md]] (8/8 sympy PASS)"
sister_cycle_architecture: "[[../op-L01-N1-EM-trace-anomaly-TGP-2026-05-11/Phase2_results.md]] (g_eff curvature reduction analog)"
tags:
  - phase2
  - thermal-field-theory
  - lattice-QCD-EoS
  - interaction-measure
  - Friedmann-modification
  - Q2-reduction-verification
  - crossover-not-phase-transition
---

# Phase 2 setup

## §0 — Cel Phase 2

Wyprowadzić **thermal-dependent ρ_QCD(T)** z lattice QCD equation-of-state +
zintegrować z **Friedmann equation w QCD epoce** (z~10¹², T~T_c~156 MeV) +
**konstruktywnie potwierdzić Q2 F1** (ρ_QCD(T<<Λ_QCD) → 0 strukturalnie, substrate-
decoupled od bare Λ).

Sub-needs: N0.5 (thermal profile), N0.6 (Friedmann modification), N0.7 (Q2 reduction).
Risks: R3 (BBN preview), R5 (Q2 konstruktywna verification), R7 (crossover smoothness).

## §1 — Setup: thermal field theory dla QCD trace anomaly

### §1.1 — Equilibrium thermodynamics — interaction measure

W termalnej QCD (perfect fluid limit), stress-energy tensor w rest frame:
```
T^μ_ν^{thermal} = diag(ε(T), -p(T), -p(T), -p(T))
T^μ_μ = ε(T) - 3p(T)            (rest-frame trace)
```

z `ε(T)` energy density, `p(T)` pressure.

**Definicja "interaction measure" (canonical termodynamics)**:
```
Δ(T) ≡ (ε - 3p) / T⁴
```

To jest *natural dimensionless quantity* — measures **departure od conformality**:
- Free radiation (massless, conformal): p = ε/3 ⇒ Δ = 0
- Strong interaction (non-conformal): p < ε/3 ⇒ Δ > 0 (positive)

W TGP per L01:
```
ρ_QCD(T) = -T^μ_μ_QCD(T) / c_0²
        = -(ε - 3p) / c_0²
        = -Δ(T)·T⁴ / c_0²       [w naturlanych jednostkach]
```

**Wage:** dla radiation-only case (T<<Λ_QCD, gluons confined w hadrons), Δ → 0
strukturalnie → ρ_QCD(T<<Λ_QCD) = 0. To jest **Q2 F1 konstruktywna verification**.

### §1.2 — Lattice QCD equation-of-state (HotQCD + Wuppertal-Budapest)

**HotQCD collaboration (Bazavov et al. 2014+) lattice results** dla 2+1 flavor:

| T regime | Δ(T) profile | Physical interpretation |
|---|---|---|
| **T < 100 MeV** (hadronic) | Δ(T) → 0 | Confinement; only vacuum condensate (constant) |
| **T = T_c ≈ 156 MeV** | Δ_max ≈ 4 (peak) | Crossover transition; gluons becoming free |
| **T = 200 MeV** | Δ ≈ 3 | Mixed phase, hadrons + free gluons |
| **T = 400 MeV** | Δ ≈ 1 | Free QGP, asymptotic freedom kicks in |
| **T → ∞** | Δ → 0 (logarithmic) | Stefan-Boltzmann free QGP limit |

**T_c definition:** chiral susceptibility peak temperature, lattice consensus 156 ± 9 MeV
(HotQCD 2018+, Wuppertal-Budapest 2020+).

**Crossover, NOT phase transition:** dla 2+1 flavors lattice pokazuje Δ(T) jest
*smooth function* z peak (NOT discontinuity); thermodynamic analytical continuation
przy T=T_c. **R7 risk addressed:** ρ_QCD(T) jest smooth, NIE first-order z latent heat.

### §1.3 — Two contributions: vacuum + thermal

**Important separation:**

```
ρ_QCD(T) = ρ_QCD_vacuum + ρ_QCD_thermal(T)
```

z:
- **ρ_QCD_vacuum** = constant ≈ -⟨α_s G²/π⟩ × 9/8 ≈ -0.0135 GeV⁴ ≈ 2.8·10¹⁸ kg/m³
  (Phase 1 §2.2; SVZ-1979 input). To jest **substrate-decoupled** od bare Λ per Q2 F1.
- **ρ_QCD_thermal(T)** = -Δ_lattice(T) × T⁴ / c_0² (lattice EoS profile);
  dynamic w hot epoch. **To jest source dla Φ-EOM w QCD epoce.**

**Konsekwencje w T-limits:**

| T regime | ρ_QCD_vacuum | ρ_QCD_thermal(T) | Total ρ_QCD |
|---|---|---|---|
| T → 0 (today, T_CMB ~ 2.7 K = 0.23 meV) | 2.8·10¹⁸ kg/m³ | 0 (vanishingly small) | constant gluon condensate |
| T = T_c ≈ 156 MeV | 2.8·10¹⁸ kg/m³ | -Δ_max·T_c⁴/c² ≈ 1.04·10²⁵ kg/m³ | thermal dominates |
| T → ∞ (early universe T >> T_c) | 2.8·10¹⁸ kg/m³ | → 0 (Stefan-Boltzmann conformal) | vacuum dominates (small) |

**Per Q2 cycle F1 substrate-decoupling:** ρ_QCD_vacuum NIE additive do bare Λ —
absorbed structurally przez single-Φ axiom + substrate-vacuum identification.
ρ_QCD_thermal(T) jest **transient phase-transition source** w hot QCD epoch.

## §2 — Friedmann equation modyfikacja w QCD epoce

### §2.1 — Standard FRW Friedmann (radiation-dominated era)

W standard cosmology:
```
H²(T) = (8π G_N / 3) · ρ_total(T)
ρ_total(T) ≈ ρ_radiation(T) = (π²/30) · g_*(T) · T⁴
```

z `g_*(T)` effective relativistic DOF:
- T >> 1 GeV: g_* = 106.75 (all SM particles relativistic)
- T = T_c ~ 156 MeV: g_* ≈ 47 (post-QCD: photons + neutrinos + light leptons +
  light quarks + gluons, gradually decoupling)
- T_c ≈ 156 MeV (transition phase, lattice precise number)
- T << 100 MeV: g_* ≈ 17 (hadronic phase, post-confinement)
- T ~ MeV (BBN): g_* ≈ 10.75 (photons + 3 neutrino species)

### §2.2 — TGP modyfikacja w QCD epoce

W TGP per emergent-metric + L01 framework:
```
□Φ̄ = -V'(Φ̄) - (q/Φ_0)·⟨φ·ρ⟩_matter         [sek08a eq:Phi-EOM]
```

z source `⟨φ·ρ⟩_matter` zawiera **complete stress-energy tensor**:
```
ρ_total_TGP(T) = ρ_radiation(T) + ρ_QCD_thermal(T) + ρ_other_matter(T)
```

Przy T~T_c, **interaction measure peak** Δ_max·T_c⁴ daje:
```
ρ_QCD_thermal(T_c) ≈ 4 · (0.156 GeV)⁴ ≈ 4 · 5.93·10⁻⁴ GeV⁴ 
                  ≈ 2.37·10⁻³ GeV⁴
```

**Compare to ρ_radiation(T_c):**
```
ρ_radiation(T_c) = (π²/30) · 47 · (0.156)⁴ ≈ 0.0917 GeV⁴
```

⇒ **ρ_QCD_thermal/ρ_radiation ≈ 2.6%** at T=T_c. **Dominant term jest ρ_radiation;
ρ_QCD_thermal jest small but non-negligible correction.**

### §2.3 — Modyfikacja H(T) w QCD epoce

W ramach standard ΛCDM cosmology, **interaction measure** is *included* w lattice
EoS calculations of g_*(T). Modyfikacja H(T) w QCD epoce jest:

```
H_TGP(T) = H_standard(T) · (1 + δ(T))
δ(T) ≈ ρ_QCD_thermal/ρ_radiation × O(1)  ~ few percent at T_c, → 0 elsewhere
```

**Tego cyklu key finding:** TGP modyfikacja H(T) w QCD epoce nie wprowadza *new*
contribution beyond what standard cosmology + lattice QCD already include. Standard
BBN/CMB analysis correctly accounts for QCD interaction measure in g_*(T) —
TGP framework jest *consistent* z tym.

### §2.4 — IR limit: ρ_QCD(T<<Λ_QCD) → 0 (Q2 F1 verification)

**Strukturalna konsekwencja:**

W limicie T → 0 (today, T_CMB ~ 2.7 K = 0.23 meV << Λ_QCD = 217 MeV):

1. **Thermal gluons confined w hadrons** (no free gluons w macroscopic universe).
2. **Δ(T<<Λ_QCD) → 0** lattice consensus (HotQCD: Δ vanishes exponentially below T~100 MeV).
3. **ρ_QCD_thermal(T → 0) → 0** strukturalnie.
4. **ρ_QCD_vacuum** (SVZ gluon condensate, ~2.8·10¹⁸ kg/m³) — **substrate-decoupled**
   per Q2 F1 mechanism.

⇒ **Q2 F1 konstruktywnie verified dla QCD sektora:**

> ρ_QCD(today) ≡ ρ_QCD_vacuum (constant gluon condensate) + ρ_QCD_thermal(T~0) (≈0)
> ρ_QCD_vacuum NIE additive do bare Λ_TGP (Q2 F1 mechanism)
> ⇒ ρ_QCD(today) NIE contribuuje do ρ_vac_TGP = M_Pl²·H₀²·g̃/12

**To zachowuje T-Λ ratio empirical 1.020 (2026-04-26 closure 7/7 PASS).**

## §3 — Phase 2 plan + sympy LOCK targets

### §3.1 — Phase 2 sympy targets (8 tests)

Phase2_sympy.py będzie weryfikować:

1. **Dimensional analysis ρ_QCD(T):** `[ρ_QCD(T)] = [GeV⁴]` (energy density);
   konwersja na mass density via `ρ_QCD/c_0² = [kg/m³]`.
2. **Interaction measure peak:** `Δ_max ≈ 4` near `T_c = 156 ± 9 MeV` (HotQCD lattice value).
3. **Stefan-Boltzmann limit:** at `T >> T_c`, `p → ε/3`, `Δ → 0` (asymptotic freedom);
   verify że Δ(T) jest *finite peak* nie diverging.
4. **IR limit:** at `T << Λ_QCD`, `Δ(T) → 0` exponentially; only vacuum condensate
   pozostaje (Q2 F1 verification basis).
5. **Friedmann equation w QCD epoce:** `H²(T) = (8πG_N/3)·ρ_total(T)` z
   `ρ_total = ρ_radiation + ρ_QCD_thermal`; ratio `ρ_QCD_thermal/ρ_radiation ≈ 2.6%`
   at T=T_c.
6. **Q2 F1 konstruktywna verification:** ρ_QCD_vacuum (substrate-decoupled) +
   ρ_QCD_thermal(T → 0) → 0 → today's Λ NIE dostaje QCD additive contribution.
7. **R7 verification:** 2+1 flavor lattice consensus: crossover (smooth Δ(T)),
   NIE first-order phase transition (no discontinuity).
8. **R5 cross-check Q2 F1:** explicit że ρ_QCD(today) NIE rozbije T-Λ ratio 1.020
   empirical match.

Target: 8/8 sympy PASS.

### §3.2 — Phase 2 deliverables

- [[Phase2_setup.md]] (this file)
- [[Phase2_results.md]] — full cosmology integration + Q2 F1 verification
- [[Phase2_sympy.py]] + [[Phase2_sympy.txt]] — 8 tests

## §4 — Risk addressing in Phase 2

### §4.1 — R3 (BBN setup, full check Phase 3)

**Strategy:** Phase 2 sets up framework dla BBN bound check w Phase 3:
- BBN era T ~ 1 MeV << Λ_QCD = 217 MeV ⇒ deeply w hadronic regime
- ρ_QCD_thermal(T~1 MeV) ≈ 0 strukturalnie
- ⇒ H(z~10⁹) standard ΛCDM preserved
- ⁴He Y_p = 0.245 ± 0.003 (PDG 2024) preserved automatic

Phase 3 will explicit numerical verification.

### §4.2 — R5 (Q2 inconsistency check)

**Strategy:** **konstruktywna** verification w §1.3 + §2.4:
- ρ_QCD_vacuum + ρ_QCD_thermal(T→0) → ρ_QCD_vacuum (constant)
- Q2 F1 mechanism (single-Φ + substrate-vacuum) decoupling preserved bezwarunkowo
- T-Λ ratio 1.020 empirical match preserved

**Gdyby tego cyklu derivation pokazałaby ρ_QCD(today) > 0 jako contribution do bare
Λ**, to byłoby falsification Q2 F1 i T-Λ closure. Lattice consensus + dimensional
analysis preserve obie.

### §4.3 — R7 (crossover not phase transition)

**Strategy:** **lattice consensus 2+1 flavor crossover** (HotQCD 2018+, Wuppertal-
Budapest 2020+) — smooth Δ(T) profile bez discontinuity. Per Phase 1 §1.4
description i Phase 2 §1.2 lattice values, ρ_QCD(T) jest smooth function.

**Konsekwencja:** brak strong stochastic GW background z first-order phase
transition; preserve compatibility z PTA NANOGrav 15-yr consensus (SMBHB origin).

## §5 — Connection do Phase 3

Phase 2 daje **thermal ρ_QCD(T) profile + Friedmann modyfikacja + Q2 reduction
verified**.

Phase 3 (next session) dostanie **observational bounds**:
1. BBN ⁴He, D/H predictions z H(z~10⁹) preserved
2. CMB ω_b, ω_m preservation (matter-decoupling)
3. PTA NANOGrav 15-yr compatibility (NO first-order signal)

## §6 — Cross-references

- [[./README.md]]
- [[./Phase0_balance.md]] §3 NEEDS, §4 6/6 gate
- [[./Phase1_results.md]] (β_QCD LOCK + gluon condensate)
- [[./NEEDS.md]] N0.5, N0.6, N0.7
- [[../op-Q2-vacuum-budget-2026-05-10/Phase_FINAL_close.md]] §2.3 (transient sources)
- [[../op-Q2-vacuum-budget-2026-05-10/FINDINGS.md]] F1, F2, F3 (matter-decoupling target)
- [[../closure_2026-04-26/Lambda_from_Phi0/results.md]] (T-Λ ratio 1.020 baseline)
- [[../op-emergent-metric-from-interaction-2026-05-09/Phase1_results.md]] (g_eff)
- HotQCD collaboration (Bazavov et al. 2014+) — lattice EoS Δ(T)
- Wuppertal-Budapest collaboration (Borsanyi et al. 2020+) — lattice EoS update
- PDG 2024 — α_s, Λ_QCD, T_c, BBN parameters

---

**Phase 2 setup ready.** Next: Phase2_sympy.py + Phase2_results.md.
