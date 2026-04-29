---
title: "ζ.1.Phase1 results — Σm_ν = 59.0 meV neutrino mass-spectrum LOCKED"
date: 2026-04-29
cycle: ζ.1.Phase1
status: CLOSED
verdict: PASS
predecessor: "[[../op-eps-photon-ring/Phase3_results.md]]"
parent: "[[program.md]]"
tags:
  - TGP
  - zeta-mass-spectrum
  - neutrino-masses
  - Koide-K-half
  - Majorana
---

# ζ.1.Phase1 — Results: Σm_ν = 59.0 meV neutrino mass-spectrum LOCKED

> **Status:** CLOSED 2026-04-29 — **5/5 PASS**.
> Σm_ν = **59.01 meV** (NO ordering) LOCKED via K(ν) = 1/2 (Majorana, B²=1)
> + NuFit 5.3 Δm² observational inputs (drift 0.99% z TGP target 59.6 meV);
> normal ordering EXCLUSIVE; inverted ordering **structurally FORBIDDEN**
> by K=1/2 incompatibility; within DESI DR2 0.072 eV bound (margin +22%);
> falsifiable by DESI DR3 2027+ (TGP 59 meV vs projected 40 meV bound).

---

## Verdict

| Sub-test | Description | Result |
|---|---|---|
| **Z1.1** | K(ν) = 1/2 sympy exact (Majorana, B²=1, chirality-counting) | **PASS** |
| **Z1.2** | NuFit 5.3 Δm²₂₁ = 7.53e-5, \|Δm²₃₁\| = 2.453e-3 eV² inputs | **PASS** |
| **Z1.3** | Σm_ν = 59.01 meV (drift 0.99% from 59.6 meV target) | **PASS** |
| **Z1.4** | DESI DR2 0.072 eV bound vs TGP 0.059 eV (margin +22%) | **PASS** |
| **Z1.5** | Inverted ordering FORBIDDEN by K(ν)=1/2 (no real root) | **PASS** |

**5/5 PASS** → Σm_ν mass-spectrum LOCKED; Phase 2 PMNS proceeds.

---

## Z1.1 — K(ν) = 1/2 sympy exact

```
Koide formula  K = (2 + B²) / (2N)     for N=3 generations
Lepton (Dirac, B²=2)        K_lep = (2+2)/6 = 4/6 = 2/3 = 0.66667
Neutrino (Majorana, B²=1)   K_ν   = (2+1)/6 = 3/6 = 1/2 = 0.50000
K_lep matches PDG (2/3) to 10⁻⁵ level                     True
K_ν = 1/2 sympy exact rational                            True
```

**Verdict:** PASS — chirality-counting w Koide formuła wynika z group
theory (Dirac 2 chiralities → B²=2; Majorana single chirality → B²=1).
K_lep = 2/3 i K_ν = 1/2 sympy exact rationals; obie zawarte w
DERIVED status taxonomy.

## Z1.2 — Δm² observational inputs (NuFit 5.3)

```
Δm²₂₁ (KamLAND + solar)                       = 7.53 · 10⁻⁵ eV²
|Δm²₃₁| (atmospheric + reactor + accelerator) = 2.453 · 10⁻³ eV²
Hierarchy ratio                                = 32.6
```

**Verdict:** PASS — NuFit 5.3 Δm² inputs sane (positive, hierarchy ratio
~30); both observational (not first-principles z TGP substrate); inputs
do TGP framework dla mass-spectrum closure.

## Z1.3 — Σm_ν = 59.01 meV via K(ν) = 1/2 closure

```
m_1 (lightest, NO)                            = 0.76 meV
m_2 = √(m₁² + Δm²₂₁)                          = 8.71 meV
m_3 = √(m₁² + |Δm²₃₁|)                        = 49.53 meV
Σm_ν = m₁ + m₂ + m₃                           = 59.01 meV
TGP target                                    = 59.6 meV
Drift |Σ − target|/target                     = 0.99%
K closure at m₁                                = 0.50000  (target 0.5)
```

**Verdict:** PASS — bisection na K(ν) = 1/2 daje m₁ = 0.76 meV (lightest)
i Σm_ν = 59.01 meV. Drift 0.99% z TGP target 59.6 meV reflects minor
NuFit 5.3 vs original-quoted Δm² values (acceptable < 5% gate).
Normal ordering self-consistent z K = 1/2 closure.

## Z1.4 — DESI DR2/DR3 bound vs TGP

```
TGP Σm_ν                                       = 59.01 meV = 0.0590 eV
DESI DR2 (current 95% CL)                      < 72 meV
DESI DR3 (projected 2027+ 95% CL)              < 40 meV
Margin DR2 (bound − TGP)/TGP                   = +22.0%
Margin DR3 (bound − TGP)/TGP                   = −32.2%
Within DR2 bound?                              True
Could be falsified by DR3?                     True (LIVE)
```

**Verdict:** PASS — TGP Σm_ν = 59.0 meV mieści się w DESI DR2 bound 72 meV
z 22% margin. **DR3 (2027+) projected bound 40 meV** mogłaby falsyfikować
TGP — to **feature** (LIVE prediction), nie bug. Odpowiada predicted
falsification gate dla ζ.1.

## Z1.5 — Inverted ordering FORBIDDEN by K(ν) = 1/2

```
IO bisection: m_3 (lightest) ∈ [10⁻⁶, 1] eV
IO Koide residual at m₃ = 1 μeV               = −0.00222
IO Koide residual at m₃ = 1 eV                = −0.16667
No sign change in residual → NO IO root z K=1/2 closure
```

**Verdict:** PASS — IO ordering ma K range ≈ [0.498, 0.333] (od m₃→0
do m₃→degenerate); **K = 0.5 NIE leży w tym range** (close ale poza).
IO **strukturalnie FORBIDDEN** by K(ν) = 1/2. Confirms TGP normal-ordering
exclusivity z chirality-counting framework.

---

## Synthesis

ζ.1.Phase1 zamyka 5-sub-test neutrino mass-spectrum audit:

1. **K(ν) = 1/2** sympy exact via Majorana B²=1 chirality-counting
2. **NuFit 5.3 Δm² inputs** sanity check (Δm²₂₁=7.53e-5, |Δm²₃₁|=2.453e-3)
3. **Σm_ν = 59.01 meV** via K=1/2 closure bisection na m₁=0.76 meV
4. **DESI DR2 within bound** (+22% margin); **DR3 falsifiable** (-32%)
5. **Inverted ordering FORBIDDEN** by K=1/2 incompatibility

**Conclusion:** TGP K(ν) = 1/2 + NuFit Δm² → unique normal-ordering closure
z Σm_ν ≈ 59 meV; phase 2 może proceed dla PMNS angles z structural
GL(3,𝔽₂) × Z₃ × SU(2)_L derivation.

---

## What ζ.1.Phase1 closes

- ✅ K(ν) = 1/2 sympy exact LOCKED (Majorana B²=1)
- ✅ NuFit 5.3 Δm² inputs (observational)
- ✅ Σm_ν = 59.01 meV NO ordering closure (drift 0.99%)
- ✅ DESI DR2 bound 0.072 eV margin +22% (LIVE)
- ✅ Inverted ordering FORBIDDEN structurally

## What ζ.1.Phase1 does NOT close

- ❌ PMNS angle derivations (Phase 2)
- ❌ Cross-sector lepton-quark θ_C-θ₁₃ unification (Phase 3)
- ❌ Δm² first-principles derivation z TGP substrate (long-term, neutrino
  MSW track)

---

## Materiał wykonawczy

- **Skrypt:** [`phase1_neutrino_audit.py`](phase1_neutrino_audit.py)
- **Output:** [`phase1_neutrino_audit.txt`](phase1_neutrino_audit.txt)
- **Setup:** [`Phase1_setup.md`](Phase1_setup.md)

## Cross-references

- [`program.md`](program.md) — overall ζ.1 plan
- [`../../PREDICTIONS_REGISTRY.md`](../../PREDICTIONS_REGISTRY.md) — C3 Σm_ν entry (LIVE 2027+)
- [`../op-eps-photon-ring/Phase3_results.md`](../op-eps-photon-ring/Phase3_results.md) — ε.1 program END

## Decyzja po Phase 1

**ζ.1.Phase1 CLOSED** with 5/5 PASS.

→ **Proceed Phase 2** (PMNS first-principles z GL(3,𝔽₂) × Z₃ × SU(2)_L, 7 sub-tests).
   Master ledger update: 391 → 396 (+5 z Phase 1).
