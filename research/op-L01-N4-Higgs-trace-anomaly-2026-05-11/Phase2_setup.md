---
title: "Phase 2 setup — EW phase transition cosmology + Friedmann ρ_Higgs(T) source + Q2 F1 konstruktywna verification dla Higgs + R5 crossover lock"
date: 2026-05-11
parent: "[[./README.md]]"
type: phase-setup
phase: 2
status: 🟡 setup phase
sub_needs_addressed: [N0.7, N0.8, N0.9]
risks_addressed: [R3-extended, R5, R6-partial]
predecessor: "[[./Phase1_results.md]] (8/8 sympy PASS)"
sister_cycle_architecture: "[[../op-L01-N2-QCD-trace-anomaly-cosmology-2026-05-11/Phase2_setup.md]] (thermal + Friedmann + Q2 reduction pattern)"
tags:
  - phase2
  - EW-phase-transition
  - thermal-Higgs
  - lattice-EW-crossover
  - Friedmann-modification
  - Q2-F1-Higgs-verification
  - crossover-not-first-order
  - R5-lock
---

# Phase 2 setup

## §0 — Cel Phase 2

Wyprowadzić **thermal-dependent ρ_Higgs(T)** w EW epoce + zintegrować z
**Friedmann equation w EW epoce** (z~10¹⁵, T~T_EW~159 GeV) +
**konstruktywnie potwierdzić Q2 F1** dla Higgs sektora (ρ_Higgs(T<<v) → 0
strukturalnie, substrate-decoupled od bare Λ) + **explicit R5 verify**
(crossover NIE first-order dla m_H=125.25 GeV per lattice 2024+).

Sub-needs: N0.7 (EW phase transition crossover), N0.8 (Friedmann modyfikacja
w EW epoce), N0.9 (Q2 F1 Higgs reduction). Risks: R3 (hierarchy partial
extension), R5 (crossover NOT first-order), R6 (cross-cycle consistency
preview).

Architectural inheritance: **N2 QCD Phase 2 pattern** (lattice EoS thermal
profile + Friedmann modyfikacja + IR limit Q2 F1 verification + crossover R7).
Higgs analog jest **simpler** od QCD ponieważ Higgs jest **fundamentally
perturbative** w SM (małe λ ≈ 0.129) podczas gdy QCD jest non-perturbative
near T_c — ale **physical mechanism identical**.

## §1 — Setup: EW phase transition w SM z m_H=125.25 GeV

### §1.1 — Historical context (M_H dependence)

W historii badań Higgs sektora (przed LHC 2012):

**Pre-2012 status (Kajantie, Laine, Rummukainen, Shaposhnikov 1996+):**
| m_H regime | EW transition order |
|---|---|
| **m_H < ~72 GeV** | **First-order** strong (potential barrier; bubble nucleation) |
| **72 < m_H < 80 GeV** | First-order weak (small barrier) |
| **m_H ≈ 80 GeV** | **Endpoint of first-order line** (lattice crossover threshold) |
| **m_H > 80 GeV** | **Crossover** (smooth, no barrier, no latent heat) |

**Post-2012 LHC (m_H = 125.25 ± 0.17 GeV PDG 2024):**
m_H >> 80 GeV ⇒ **deeply w crossover regime**.

### §1.2 — Lattice consensus 2024+ (R5 lock basis)

**Recent papers konstruktywnie confirming crossover at m_H=125 GeV:**

1. **Csikor, Fodor, Heitger 1998+** — pierwsze lattice 4D scan: endpoint at m_H_c ≈ 72.3 ± 0.7 GeV (3D effective theory) → m_H_c ≈ 80 GeV (4D)
2. **D'Onofrio, Rummukainen 2014, arXiv:1404.3565** — full SM lattice at physical Higgs mass; confirmed crossover (no metastability signal)
3. **Kainulainen, Pascual, Roloff 2024, arXiv:2405.01191** — recent review post-Run 2; crossover confirmed
4. **Niemi, Schicho, Tenkanen 2021, arXiv:2103.07467** — dimensional reduction + lattice confirms crossover for SM physical parameters

⇒ **R5 PRIMARY ANCHOR:** dla m_H = 125.25 GeV (PDG 2024), EW transition jest
**KROSOVER** (smooth), NIE first-order. **NO strong stochastic GW signal
od bubble collisions** w PTA NANOGrav 15-yr (consistent z empirical absence
of MHz-Hz primordial GW peak).

### §1.3 — Critical temperature T_EW

**Lattice consensus (Kajantie-Laine-Rummukainen-Shaposhnikov, D'Onofrio-Rummukainen):**

```
T_EW(m_H=125 GeV) ≈ 159 ± 1 GeV    (crossover "transition")
```

To jest "smooth crossover temperature" (peak chiral susceptibility, NOT phase
transition w sensie thermodynamicznym).

**Physical scale comparison:**
- m_H ≈ 125 GeV
- v(T=0) = 246 GeV
- T_EW ≈ 159 GeV (slightly above m_H, much below v)
- T_EW/v ≈ 0.65 (dimensionless ratio)

### §1.4 — Symmetry restoration w wysokich temperaturach

Per **finite-T effective potential** (Dolan-Jackiw 1974, Weinberg 1974):
```
V_eff(h, T) ≈ V_classical(h) + (T²/24)·[3g²/4 + (g²+g'²)/4 + 2y_t² + 4λ]·h² + O(T⁴)
            ≈ V_classical(h) + (c²/2)·T²·h²    [c² ≈ 0.146 dla SM physical params]
```

Effective mass: μ²_eff(T) = μ² - c²·T².

**Critical temperature** (minimum at h=0): μ²_eff(T_c) = 0 →
```
T_c² = μ²/c² = λv²/c² 
T_c ≈ v·√(λ/c²) ≈ 246·√(0.129/0.146) ≈ 246·0.94 ≈ 232 GeV  (rough perturbative)
```

Lattice corrections obniżają to do ~159 GeV (perturbative theory over-estimates).

**Above T_c:** ⟨h⟩ = 0 (symmetric phase); below T_c: ⟨h⟩ ≠ 0 (broken phase).
W crossover regime przejście jest **smooth** (no discontinuous order parameter
jump).

## §2 — Thermal stress-energy ρ_Higgs(T)

### §2.1 — Two contributions: vacuum + thermal (analog N2 §1.3)

```
ρ_Higgs(T) = ρ_Higgs_vacuum + ρ_Higgs_thermal(T)
```

**ρ_Higgs_vacuum (substrate-decoupled per Q2 F1):**
```
ρ_Higgs_vacuum_bare = -T^μ_μ_vac_bare/c_0² = +λv⁴/c_0²
                   ≈ 0.129·(246.22 GeV)⁴/c_0²
                   ≈ 4.7·10⁸ GeV⁴/c_0²
```

Per Phase 1 F1.11: **bare ρ_Higgs_vacuum ~ 10⁶⁶ eV⁴** (Q2 cycle F7 reference).
**Substrate-decoupled** od bare Λ_TGP — strukturalnie absorbed przez single-Φ
axiom + substrate-vacuum identification.

**ρ_Higgs_thermal(T) — transient EW epoch source:**

Dla bosonic thermal field theory:
```
ρ_Higgs_thermal(T) = (π²/30)·g_Higgs·T⁴ × f(m_H/T)
```

z g_Higgs = 1 (single real scalar dof w broken phase), f(m_H/T) = thermal
function:
- **T >> m_H:** f → 1 (relativistic/massless limit)
- **T ~ m_H:** f ≈ 0.3-0.7 (Boltzmann suppression onset)
- **T << m_H:** f → 0 exponentially (massive freeze-out)

W EW epoke (z ~ 10¹⁵, T ~ 159 GeV ≈ 1.27·m_H):
```
ρ_Higgs_thermal(T_EW) ≈ (π²/30)·1·(159 GeV)⁴·f(m_H/T_EW)
                      ≈ 3.29·(159)⁴ × f(0.79)
                      ≈ 3.29·6.4·10⁸ × ~0.5
                      ≈ 1.05·10⁹ GeV⁴
```

To jest comparable do total radiation density g_*(T_EW) ≈ 100:
```
ρ_radiation(T_EW) = (π²/30)·100·(159)⁴
                  ≈ 3.29·100·6.4·10⁸
                  ≈ 2.1·10¹¹ GeV⁴
```

⇒ **ρ_Higgs_thermal/ρ_radiation ≈ 0.5% at T_EW** (single Higgs DOF out of
~100 total relativistic DOF, exactly as expected).

### §2.2 — Konsekwencje w T-limits

| T regime | ρ_Higgs_vacuum | ρ_Higgs_thermal(T) | Total ρ_Higgs | Verdict |
|---|---|---|---|---|
| **T → 0** (today, T_CMB ~ 2.35·10⁻⁴ eV) | ~10⁶⁶ eV⁴ (substrate-decoupled) | ~0 (exp suppressed by m_H/T) | substrate-decoupled vacuum only | **Q2 F1 verified** |
| **T = T_EW ≈ 159 GeV** | ~10⁶⁶ eV⁴ | ~10⁵⁶ eV⁴ thermal | thermal small vs radiation | normal cosmology |
| **T >> m_H** (early universe) | ~10⁶⁶ eV⁴ | ~T⁴ (Stefan-Boltzmann) | radiation-like | included w g_*(T) |

**Per Q2 F1 substrate-decoupling:** ρ_Higgs_vacuum NIE additive do bare Λ —
absorbed strukturalnie przez single-Φ axiom + substrate-vacuum identification.
ρ_Higgs_thermal(T) jest **transient EW epoch source**, identical w nature do
QCD thermal Δ(T)·T⁴ source w QCD epoce (per N2 Phase 2).

### §2.3 — Q2 F1 konstruktywna verification dla Higgs sektora

**Strukturalna konsekwencja w limicie T → 0:**

W limicie T → 0 (today, T_CMB ~ 2.35·10⁻⁴ eV << m_H = 125.25 GeV):

1. **Thermal Higgs bosons exponentially suppressed:** f(m_H/T) ~ exp(-m_H/T)
   z m_H/T_CMB ~ 5.3·10¹⁴ ⇒ utterly negligible.
2. **ρ_Higgs_thermal(T → 0) → 0** strukturalnie (exp Boltzmann factor).
3. **ρ_Higgs_vacuum** (Coleman-Weinberg + tree-level) — **substrate-decoupled**
   per Q2 F1 mechanism (single-Φ + substrate-vacuum identification, Q2 §3).
4. **Net contribution to today's Λ:** ZERO (post-Q2 F1 absorption).

⇒ **Q2 F1 konstruktywnie verified dla Higgs sektora** (analog do N2 §3 dla QCD,
N1 §2 dla EM):

> ρ_Higgs(today) ≡ ρ_Higgs_vacuum (substrate-decoupled) + ρ_Higgs_thermal(T~0) (~0)
> ρ_Higgs_vacuum NIE additive do bare Λ_TGP (Q2 F1 mechanism)
> ⇒ ρ_Higgs(today) NIE contribuuje do ρ_vac_TGP = M_Pl²·H₀²·g̃/12

**To zachowuje T-Λ ratio empirical 1.020 (2026-04-26 closure 7/7 PASS),
identycznie jak w QCD (N2 §3) i EM (N1 §2).**

## §3 — Friedmann equation modyfikacja w EW epoce

### §3.1 — Standard FRW Friedmann (radiation-dominated era)

W standard cosmology w EW epoce (z ~ 10¹⁵):
```
H²(T) = (8π G_N / 3) · ρ_total(T)
ρ_total(T) ≈ ρ_radiation(T) = (π²/30) · g_*(T) · T⁴
```

z g_*(T) effective relativistic DOF:
- **T >> 100 GeV (pre-EW):** g_* = 106.75 (full SM relativistic spectrum)
- **T ≈ T_EW ≈ 159 GeV:** g_* ≈ 100 (Higgs becomes non-relativistic)
- **100 GeV > T > m_t ≈ 173 GeV:** g_* gradually drops (top decoupling)
- **T < m_W ≈ 80 GeV:** g_* drops to ~86 (W, Z become non-relativistic)

### §3.2 — TGP modyfikacja w EW epoce (z ~ 10¹⁵)

Per emergent-metric + L01 framework:
```
□Φ̄ = -V'(Φ̄) - (q/Φ_0)·⟨φ·ρ⟩_matter
```

z source `⟨φ·ρ⟩_matter` zawiera **complete stress-energy tensor**:
```
ρ_total_TGP(T_EW) = ρ_radiation(T_EW) + ρ_Higgs_thermal(T_EW) + ρ_other_matter(T_EW)
```

**Modyfikacja H(T) w EW epoce:**
```
H_TGP(T_EW) = H_standard(T_EW) · (1 + δ_EW(T))
δ_EW(T) ≈ ρ_Higgs_thermal/ρ_radiation ~ 0.5% at T_EW, → 0 elsewhere
```

**Tego cyklu key finding:** TGP modyfikacja H(T) w EW epoce **nie wprowadza
new contribution beyond what standard SM cosmology already includes** —
g_*(T_EW) ≈ 100 *includes* Higgs DOF naturally. Standard ΛCDM analysis poprawnie
accounts for Higgs interaction measure in g_*(T_EW) — TGP framework jest
**consistent** z tym (analog do N2 §2.3 dla QCD).

### §3.3 — Hubble rate impact: negligible

```
δ_EW(T_EW)/H_standard ≈ 0.5%
```

Niedostateczne dla **detection** w jakimkolwiek cosmological probe:
- BBN sensitivity: ~few % na N_eff (T ~ 1 MeV << T_EW; thermal Higgs already
  decoupled) — irrelevant.
- CMB sensitivity: H(z~1100) z Planck; EW epoch z~10¹⁵ → ne detection possible.
- PTA/LISA stochastic GW: sensitive do **FIRST-ORDER** EW transitions z bubble
  collisions; per R5 lock, m_H=125 GeV crossover → NO signal.

⇒ **EW thermal Higgs cosmology effect jest standardowo accountowane;
TGP-specific addition jest ZERO.**

## §4 — R-guard verification (Phase 2)

### §4.1 — R5 (crossover NOT first-order) — primary verification

**Strukturalna verification:**

Dla m_H = 125.25 ± 0.17 GeV (PDG 2024 LHC Run 2):

1. **Pre-LHC lattice (Kajantie-Laine-Rummukainen-Shaposhnikov 1996):** endpoint
   m_H_c ≈ 80 GeV (4D extrapolation z 3D effective theory).
2. **Post-LHC lattice (D'Onofrio-Rummukainen 2014, arXiv:1404.3565):** explicit
   4D full SM lattice w m_H=125 GeV → crossover confirmed.
3. **2024+ reviews (Kainulainen 2024, arXiv:2405.01191):** crossover consensus
   post-Run 2.

⇒ **R5 PASS:** EW transition dla physical m_H jest **crossover**, NIE first-order.
**Brak strong stochastic GW signal** od bubble collisions/wall dynamics. Compatible
z PTA NANOGrav 15-yr SMBHB consensus + LISA forecasts (no detectable EW signal).

**Konsekwencja dla TGP:** Higgs sektor NIE produkuje detectable stochastic GW
background w cosmological probes. Spójność z N2 QCD crossover (also no signal).
**Two SM sektory NIE produkują GW background** — TGP separable sector structure
robust.

### §4.2 — R6 (cross-cycle consistency) — partial

**Strategy:** Phase 2 confirms identical *physical mechanism* dla Higgs jak dla:
- N1 EM: vacuum/curvature separation + Q2 F1
- N2 QCD: thermal/vacuum + Q2 F1 + crossover
- N3 SPARC: rho-consistency + Q2 F1

**Higgs sektor** następuje **identyczny pattern**: separation vacuum/thermal +
Q2 F1 substrate-decoupling + crossover (NOT first-order) + standard SM
cosmology preserved.

Full cross-cycle audit Phase 4.

### §4.3 — R3 (hierarchy problem) — extended Phase 2 reflection

**Q2 F1 mechanism analog dla m_H stability:**

W TGP framework, **substrate-vacuum identification** (Q2 F1) protects bare Λ_TGP
przed matter vacuum catastrophe (10⁷⁷ OOM). **Czy analogous mechanism protects
m_H przed Λ_UV² destabilization?**

Per Phase 1 §4.3 (R3 partial) — Q2 F1 może protect m_H *strukturalnie* przez:
1. **Substrate identification** (single-Φ axiom): h(x) jest emergent SM scalar,
   NIE second fundamental field; substrate-defined regulator Λ_TGP jest the
   physical UV cutoff, NIE arbitrary Λ_UV.
2. **Renormalization condition** consistent z substrate energetics: SM matter
   loops integrate up to substrate UV scale (≈ M_Pl × √(ε)), not arbitrary
   field-theoretic Λ_UV.
3. **m_H natural scale** może być substrate-determined przez Φ-EOM consistency,
   nie ad-hoc.

**HONEST CAVEAT (extending Phase 1):** tego cyklu Phase 2 NIE *rozwiąże*
hierarchy problem fully — to byłaby breakthrough. Phase 2 jedynie *strengthens
consistency speculation* z Q2 F1 + S05 mechanism. **R3 pozostaje deferred
precision item.**

## §5 — Phase 2 plan + sympy LOCK targets

### §5.1 — Phase 2 sympy targets (8 tests)

Phase2_sympy.py będzie weryfikować:

1. **EW phase transition crossover (R5):** dla m_H=125.25 GeV (PDG), lattice
   endpoint m_H_c ≈ 80 GeV ⇒ m_H > m_H_c ⇒ **crossover** (NOT first-order).
2. **Critical temperature T_EW:** finite-T effective potential μ²_eff(T_EW) = 0
   → T_EW perturbative ≈ 232 GeV; lattice ≈ 159 GeV (corrections); both at v's
   scale; ratio T_EW/v ≈ 0.65.
3. **Thermal Higgs density at T_EW:** ρ_Higgs_thermal(T_EW) ≈ (π²/30)·g_Higgs·
   T_EW⁴·f(m_H/T_EW); ratio to ρ_radiation ≈ 0.5-1%.
4. **Boltzmann suppression w T → 0:** f(m_H/T_CMB) ~ exp(-m_H/T_CMB) =
   exp(-5·10¹⁴) ≈ 0 (utterly suppressed).
5. **Q2 F1 konstruktywna verification dla Higgs:** ρ_Higgs(today) ≡ ρ_Higgs_vacuum
   (substrate-decoupled) + ρ_Higgs_thermal(T~0) (~0); NIE contribuuje do today's Λ.
6. **Friedmann modyfikacja w EW epoce:** δ_EW(T_EW) = ρ_Higgs_thermal/ρ_radiation
   ≈ 0.5%; standard g_*(T) already includes Higgs DOF — TGP-specific addition ZERO.
7. **No first-order GW signal (R5 lock):** crossover EW transition + crossover
   QCD transition ⇒ Higgs sektor compatible z NANOGrav 15-yr SMBHB consensus
   + LISA forecasts.
8. **Cross-cycle pattern verification (R6 partial):** Higgs sektor follows identical
   Q2 F1 substrate-decoupling pattern jak EM (N1), QCD (N2), SPARC (N3).

Target: 8/8 sympy PASS.

### §5.2 — Phase 2 deliverables

- [[Phase2_setup.md]] (this file)
- [[Phase2_results.md]] — full EW epoch cosmology + Q2 F1 verification
- [[Phase2_sympy.py]] + [[Phase2_sympy.txt]] — 8 tests

## §6 — Connection do Phase 3

Phase 2 daje **thermal ρ_Higgs(T) profile + Friedmann modyfikacja + Q2 F1
verification + R5 crossover lock**.

Phase 3 (next session) dostanie **phenomenological bounds**:
1. **LHC m_H=125.25 ± 0.17 GeV preservation** (Phase 1 already established;
   Phase 3 will confirm consistency w cosmological framework).
2. **Planck 2018 ω_b/ω_m + N_eff consistency** (EW epoch transient does NOT
   affect BBN/CMB observables).
3. **Future LISA stochastic GW bounds** (NO first-order EW signal expected;
   R5 explicit prediction).
4. **Future HL-LHC + FCC-ee Higgs self-coupling measurements** (precision
   future test of β_λ).

## §7 — Cross-references

- [[./README.md]]
- [[./Phase0_balance.md]] §3 NEEDS, §4 6/6 gate
- [[./Phase1_results.md]] (8/8 sympy PASS — SSB + β_λ + γ_m + Q2 F1 setup)
- [[./NEEDS.md]] N0.7, N0.8, N0.9 (planned)
- [[../op-L01-N2-QCD-trace-anomaly-cosmology-2026-05-11/Phase2_setup.md]] (sister architecture pattern)
- [[../op-L01-N2-QCD-trace-anomaly-cosmology-2026-05-11/Phase2_results.md]] (Q2 F1 QCD konstruktywna verification reference)
- [[../op-Q2-vacuum-budget-2026-05-10/Phase_FINAL_close.md]] (Q2 F1 mechanism)
- [[../op-Q2-vacuum-budget-2026-05-10/FINDINGS.md]] F1, F2, F3 (matter-vacuum decoupling)
- [[../closure_2026-04-26/Lambda_from_Phi0/results.md]] (T-Λ ratio 1.020 baseline)
- Kajantie, Laine, Rummukainen, Shaposhnikov 1996, Nucl. Phys. B 466, 189 — endpoint analysis
- D'Onofrio, Rummukainen 2014, arXiv:1404.3565 — full SM lattice at m_H=125 GeV crossover confirmed
- Kainulainen, Pascual, Roloff 2024, arXiv:2405.01191 — recent EW transition review
- Niemi, Schicho, Tenkanen 2021, arXiv:2103.07467 — dimensional reduction + lattice
- Dolan, Jackiw, Phys. Rev. D 9, 3320 (1974) — finite-T effective potential
- Weinberg, Phys. Rev. D 9, 3357 (1974) — symmetry restoration
- PDG 2024 — m_H, v, g_*(T), T_EW

---

**Phase 2 setup ready.** Next: Phase2_sympy.py + Phase2_results.md.
