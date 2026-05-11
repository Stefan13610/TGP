---
title: "NEEDS — op-Q2-vacuum-budget-2026-05-10 (residual open problems)"
date: 2026-05-10
parent: "[[README.md]]"
type: needs
tgp_owner: research/op-Q2-vacuum-budget-2026-05-10
tags:
  - needs
  - Q2
  - vacuum-budget
  - residual-open
---

# NEEDS — Q2 vacuum budget cycle (residual open problems)

> Q2 cycle zamknął *strukturalnie* native answer dla L01 NEEDS Q2 (matter sector
> vacua decoupled od bare Λ). Pozostają **residual open problems** związane z
> formalną derywacją mechanizmu i phase transition dynamics.

## Otwarte luki

### N1: Formalna kowariantna derywacja "renormalization scheme TGP"

**Status:** OPEN.

**Problem:** Phase FINAL §2.2 argument A2 stwierdza że w TGP renormalization scheme:

```
⟨T^μ_μ⟩_vacuum_TGP = -c_0²·V(Φ_eq)    [substrate-defined]
```

zamiast standardowego QFT zero-point sum. To jest **strukturalna konsekwencja**
single-Φ axiom + substrate-vacuum identification, ale brakuje **formalnej
kowariantnej derywacji** "dlaczego matter sector zero-point energies *strukturalnie*
nie wchodzą do Λ counterterm".

**Co potrzebne:**

1. RG flow analysis substrate Lagrangianu (w TGP-native formalism)
2. Demonstracja że single-Φ axiom forbids independent Λ counterterm
3. Explicit calculation jak `⟨ρ_QCD⟩`, `⟨ρ_Higgs⟩` etc. wchodzą do Φ-EOM
   *jako source*, ale *nie* jako addytywny term do V(Φ_eq)
4. Cross-check w finite-temperature QFT formalism

**Source:** [[Phase_FINAL_close.md]] §5.1

**Kandydat dostawcy:** dedicated cycle `op-TGP-renormalization-scheme-formal`
(estymata: ~3-4 tygodnie, formal QFT + RG flow)

**Typ:** derivation (formal renormalization theory)

**Priorytet:** ŚREDNI — argument strukturalny + empirical PASS (T-Λ 2% match)
wystarczają operacyjnie; formal derivation byłby major theoretical achievement
ale nie jest blocker dla cosmological predictions.

### N2: Phase transition dynamics (transient ρ_QCD(T), ρ_EW(T))

**Status:** OPEN.

**Problem:** Phase FINAL §2.3 stwierdza że SM matter sector vacua są **transient
sources** during phase transitions:

```
QCD epoch (T~Λ_QCD): ρ_QCD(T) ~ Λ_QCD⁴ ~ 10⁵⁵ eV⁴, transient
EW epoch (T~Λ_EW):   ρ_EW(T) ~ Λ_EW⁴ ~ 10⁴¹ eV⁴, transient
```

ale konkretne **profile czasowe** wymagają lattice QCD + thermal field theory
inputs.

**Co potrzebne:**

1. Lattice QCD computation `Λ_QCD(T)` profile across confinement transition
2. Thermal Higgs effective potential V_eff(φ, T) for SSB transition
3. EW phase transition dynamics (potencjalnie first-order w SM extensions)
4. Implikacje dla `H(z)` w odpowiednich epochs
5. Implikacje dla stochastic GW background z phase transitions (LISA target)

**Source:** [[Phase_FINAL_close.md]] §2.3, §5.2

**Kandydat dostawcy:** L01 NEEDS N2 + N5 → planowany cycle
`op-QCD-trace-anomaly-cosmology` (estymata: ~4-6 tygodni)

**Typ:** derivation (non-perturbative QCD + thermal QFT + cosmology)

### N3: Higgs vacuum quantum fluctuations

**Status:** ✅ **CLOSED 2026-05-11** — STRUCTURAL_DERIVED via dedicated cycle.

**Closure source:** [[../op-L01-N4-Higgs-trace-anomaly-2026-05-11/Phase_FINAL_close.md]]
(24/24 sympy PASS, 6/6 P-requirements RESOLVED).

**Resolution (post-N4 2026-05-11):**

1. **`T^μ_μ_Higgs_quantum` explicit form** (sympy LOCK Phase 1):
   `T^μ_μ_quantum = (1+γ_m)·m_H²·h² + (β_λ/4)·h⁴ + curvature·h² + Riegert local`
2. **β_λ ≈ -0.033** (top Yukawa dominant), **γ_m ≈ -0.027 ÷ -0.035** (1-loop EW scale)
3. **Q2 native answer DLA Higgs (verified konstruktywnie):** 1-loop Higgs
   *strukturalnie* substrate-decoupled, NIE residual contribution to bare Λ_TGP.
   OOM gap 55.3 (bare ρ_Higgs_vac vs Λ_obs) absorbed via single-Φ + substrate-
   vacuum identification (Q2 F1 mechanism).
4. **Quantitative bound:** Friedmann modyfikacja w EW epoce TGP-specific
   addition = ZERO (Higgs DOF already in standard g_*(T_EW)≈100); cosmological
   transient ρ_Higgs_thermal(T_EW)/ρ_radiation ≈ 0.83%; Boltzmann today ~exp(-5·10¹⁴).
   g̃ shift = 0 (T-Λ ratio 1.020 ± 0.02 preserved bezwarunkowo).
5. **R5 LOCK:** m_H=125.25 GeV > endpoint 80 GeV → EW crossover (NOT first-order);
   LISA Ω_GW^EW=0 falsifiable post-2035

**Q2.Q1 partial answer (post-N4):** Q2 F1 (matter-vacuum decoupling) **konstruktywnie
verified dla FOUR independent SM sektory** (N1 EM, N2 QCD, N3 SPARC, N4 Higgs)
via four independent diagnostic mechanisms. Strukturalna własność TGP, NIE
post-hoc tuning.

**Typ:** derivation (1-loop Higgs sector w curved background + EW cosmology
+ phenomenology bounds)

### N3.bis: Cluster mass deficit sterile ν compatibility (post-2026-05-11)

**Status:** ✅ **CLOSED 2026-05-11** (compatibility note dla Q2 framework).

**Source:** [[../op-cluster-mass-deficit-resolution-2026-05-11/Phase_FINAL_close.md]]

**Compatibility note:**
- Cluster cycle adopted **sterile ν 2 eV** (Angus-Famaey-Diaferio 2010 framework)
  jako matter component addressing ~35% cluster mass deficit
- **Sterile ν jest MATTER SEKTOR** (subject to Q2 F1 substrate-decoupling) —
  NIE substrate vacuum component
- ω_sterile ~ 0.0005 (cosmological abundance) → negligible w Ω_Λ ≈ 0.69 budget
- **T-Λ ratio 1.020 preserved bezwarunkowo** post-sterile-ν
- ΔN_eff = 0.05 < Planck 2018 ±0.18 (1σ) ✓
- **Q2 F1 mechanism preserved** — sterile ν cluster contribution gravitational
  (clustered), NIE substrate vacuum additive

### N4: Inflation prehistory of Φ_eq

**Status:** OPEN — beyond Q2 scope.

**Problem:** Q2 cycle zakłada że Φ_eq = H₀ jest *obecna wartość*. Pytanie:

- Czy Φ_eq = H(t) zawsze (równa Hubble parameter w każdej epoce)?
- Czy Φ_eq jest *stała* od inflation onwards?
- Jaka była dynamika Φ_eq podczas inflation, reheating, BBN?

To jest powiązane z OP-3 postulate (a_Γ = 1/Φ₀) i wymaga dedicated inflation
analysis.

**Source:** T-Λ closure §7.1.5; [[Phase_FINAL_close.md]] §5

**Kandydat dostawcy:** dedicated cycle `op-inflation-substrate-genesis`
(estymata: ~6-8 tygodni, inflation + substrate dynamics)

**Typ:** derivation (early-universe substrate dynamics)

**Priorytet:** NISKI dla Q2 (nie blokuje obecnego closure); WYSOKI dla long-term
TGP cosmology completeness.

## Pytania otwarte

- **Q2.Q1:** Czy strukturalna decoupling SM-vacuum / substrate-vacuum (F1) jest
  *uniwersalna* (zachowana dla wszystkich BSM extensions), czy zależy od konkretnego
  field content? — wymaga RG flow analysis (N1).

- **Q2.Q2:** Jak Q2 verdict propaguje na **modified gravity charts** (Σ, μ
  cosmological)? — Q2 mówi o background `ρ_vac`; perturbacje wymagają oddzielnego
  analysis.

- **Q2.Q3:** Czy istnieje **direct experimental probe** dla matter-decoupling
  argumentu, oddzielny od T-Λ ratio test? Możliwe kandydaty:
  - QCD phase transition GW signatures (LISA, PTA)
  - CMB B-mode anomalies during pre-recombination phase transitions
  - Λ-evolution constraints (DESI-II, Euclid)

## Blokery

| ID | Bloker | Czeka na | Source |
|----|--------|----------|--------|
| B1 | N1 (formal renormalization derivation) | dedicated formal QFT cycle | RG flow z H_Γ |
| B2 | N2 (phase transition transients) | lattice QCD + thermal field theory | dedicated cycle |
| B3 | N3 (Higgs quantum) | ~~część N1 cycle / dedicated Higgs cycle~~ **CLOSED 2026-05-11** | [[../op-L01-N4-Higgs-trace-anomaly-2026-05-11/]] (STRUCTURAL_DERIVED, 24/24 sympy PASS, Q2 F1 verified dla Higgs sektora konstruktywnie) |
| B4 | N4 (inflation prehistory) | dedicated inflation cycle | early-universe substrate |

## Closed needs (po Q2 cyklu 2026-05-10)

| ID | Co | Domknięte przez | Data |
|----|----|-----------------|------|
| L01 NEEDS Q2 (matter vacua additive?) | structural decoupling argument (F1+F2+F3+F4); empirical PASS via T-Λ ratio (F8) | Phase_FINAL_close.md §2 + §6 | 2026-05-10 |
| Native parameter audit Λ sector | 0 free params (strongest possible lock) | Phase_FINAL_close.md §4 | 2026-05-10 |
| Three-layer presentation methodology compliance | 6/6 mandatory layers verified | Phase_FINAL_close.md §6.3 | 2026-05-10 |

## Cross-references

- [[README.md]] — cykl Q2 indeks
- [[Phase_FINAL_close.md]] — pełna analiza
- [[FINDINGS.md]] — eksportowalne wyniki F1-F8
- [[../op-L01-rho-stress-energy-bridge-2026-05-04/NEEDS.md]] §Q2 — closed by Q2 cycle
- [[../closure_2026-04-26/Lambda_from_Phi0]] — parent T-Λ closure
