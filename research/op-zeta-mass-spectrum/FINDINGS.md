---
title: "FINDINGS — op-zeta-mass-spectrum"
date: 2026-05-03
parent: "[[README.md]]"
type: findings
tgp_owner: research/op-zeta-mass-spectrum
source_session: S5 (auto-extraction; cite-only per AGENT_PROTOCOL §3)
tags:
  - findings
---

# FINDINGS — op-zeta-mass-spectrum

> **Sesja 5 auto-generated** (2026-05-03). Ekstrakcja cite-only z istniejących
> plików folderu. **Żadne treści nie są wymyślone** — każdy item ma `source:`
> cytujący plik. Manual review w Sesji 6 (RESEARCH_BUS broadcasts).

## Phase results — frontmatter verdicts

| Plik | Cycle | Status | Verdict | Score | Program |
|------|-------|--------|---------|-------|---------|
| `Phase1_results.md` | `ζ.1.Phase1` | `CLOSED` | `PASS` | — | — |
| `Phase2_results.md` | `ζ.1.Phase2` | `CLOSED` | `PASS` | — | — |
| `Phase3_results.md` | `ζ.1.Phase3` | `CLOSED` | `PASS` | — | `END` |

## Numerical PASS counts (cited from .txt outputs)

- `phase1_neutrino_audit.txt` — 5/5 PASS — context: …nu)         : 1/2 = 0.50000   Sigma m_nu    : 59.6 meV target   goal          : 5/5 PASS -> proceed Phase 2 PMNS  =============================================
- `phase1_neutrino_audit.txt` — 5/5 PASS — context: …Z1.1: PASS   Z1.2: PASS   Z1.3: PASS   Z1.4: PASS   Z1.5: PASS    Cumulative: 5/5 PASS    -> Phase 1 CLOSED with 5/5 PASS   -> Sigma m_nu = 59.01 meV LOCKED vi
- `phase1_neutrino_audit.txt` — 5/5 PASS — context: …ASS   Z1.4: PASS   Z1.5: PASS    Cumulative: 5/5 PASS    -> Phase 1 CLOSED with 5/5 PASS   -> Sigma m_nu = 59.01 meV LOCKED via K(nu)=1/2 + NuFit 5.3   -> Norm
- `phase2_pmns_derivation.txt` — 7/7 PASS — context: …S   Z2.3: PASS   Z2.4: PASS   Z2.5: PASS   Z2.6: PASS   Z2.7: PASS  Cumulative: 7/7 PASS Verdict: 7/7 PASS → PMNS first-principles LOCKED          classificati
- `phase2_pmns_derivation.txt` — 7/7 PASS — context: …2.4: PASS   Z2.5: PASS   Z2.6: PASS   Z2.7: PASS  Cumulative: 7/7 PASS Verdict: 7/7 PASS → PMNS first-principles LOCKED          classification PARTIALLY DERIV
- `phase3_zeta_predictions.txt` — 6/6 PASS — context: …3 (Z3): PASS   Z3.4 (Z4): PASS   Z3.5 (Z5): PASS   Z3.6 (Z6): PASS  Cumulative: 6/6 PASS Verdict: 6/6 PASS → ζ.1 program END          classification PARTIALLY 
- `phase3_zeta_predictions.txt` — 6/6 PASS — context: …4 (Z4): PASS   Z3.5 (Z5): PASS   Z3.6 (Z6): PASS  Cumulative: 6/6 PASS Verdict: 6/6 PASS → ζ.1 program END          classification PARTIALLY DERIVED (refined)…

## TL;DR / Wynik / Verdict sections (cytaty)

### `Phase1_results.md` — Verdict

> ## Verdict
> 
> | Sub-test | Description | Result |
> |---|---|---|
> | **Z1.1** | K(ν) = 1/2 sympy exact (Majorana, B²=1, chirality-counting) | **PASS** |
> | **Z1.2** | NuFit 5.3 Δm²₂₁ = 7.53e-5, \|Δm²₃₁\| = 2.453e-3 eV² inputs | **PASS** |
> | **Z1.3** | Σm_ν = 59.01 meV (drift 0.99% from 59.6 meV target) | **PASS** |
> | **Z1.4** | DESI DR2 0.072 eV bound vs TGP 0.059 eV (margin +22%) | **PASS** |
> | **Z1.5** | Inverted ordering FORBIDDEN by K(ν)=1/2 (no real root) | **PASS** |
> 
> **5/5 PASS** → Σm_ν mass-spectrum LOCKED; Phase 2 PMNS proceeds.
> 
> ---

### `Phase2_results.md` — Verdict

> ## Verdict
> 
> | Sub-test | Description | Result |
> |---|---|---|
> | **Z2.1** | sin²θ₁₂ = 1/3 trimaximal (S₃ ⊂ GL(3,𝔽₂); drift 8.58%) | **PASS** |
> | **Z2.2** | sin²θ₂₃ = 1/2 maximal (Z₂ + K(ν)=1/2; drift 12.6%) | **PASS** |
> | **Z2.3** | sin²θ₁₃ = λ_C²/2 cross-sector Cabibbo (drift 15.6%) | **PASS** |
> | **Z2.4** | 5 alternative parameterizations FALSIFIED (>5% drift) | **PASS** |
> | **Z2.5** | PMNS unitarity UU†=I sympy + sum rule 5/6 (drift 4.7%) | **PASS** |
> | **Z2.6** | NGFP RG-stability marginal factor (1+η_N*/2)=0; LISA margin | **PASS** |
> | **Z2.7** | Classification PARTIALLY DERIVED (refined);...(truncated)

### `Phase3_results.md` — Verdict

> ## Verdict
> 
> | Sub-test | Description | Result |
> |---|---|---|
> | **Z3.1 (Z1)** | DESI DR3 2027+ Σm_ν falsification (margin DR2 +22%, DR3 −32%) | **PASS** |
> | **Z3.2 (Z2)** | JUNO 2027+ θ₁₃ precision (TGP sin²2θ₁₃ = 0.099, drift 16.5%) | **PASS** |
> | **Z3.3 (Z3)** | DUNE/T2HK 2030+ θ₂₃ octant resolution (drift 12.6%, octant LIVE) | **PASS** |
> | **Z3.4 (Z4)** | Cross-sector K-taxonomy distinct (lepton=2/3, ν=1/2, quark 0.81-0.87) | **PASS** |
> | **Z3.5 (Z5)** | Lepton-quark θ_C-θ₁₃ unification (TGP √2 ratio, drift 6.8%) | **PASS** |
> | **Z3.6 (Z6)** | 4-channel falsification convergence (4/4 ≥ 3/4 ...(truncated)

### `Phase3_results.md` — Status taxonomy (post-ζ.1)

> ### Status taxonomy (post-ζ.1)
> 
> - **K(ν) = 1/2** Majorana B²=1: **DERIVED** (chirality-counting sympy exact)
> - **Σm_ν = 59.01 meV NO ordering**: **STRUCTURAL** (K=1/2 + observational Δm²)
>   → **PARTIALLY DERIVED** post-DESI DR3 confirmation
> - **PMNS angles** (θ₁₂=trimaximal, θ₂₃=maximal, θ₁₃=Cabibbo-lock):
>   **PARTIALLY DERIVED (refined)** — 3 z 4 free parameters closed
> - **Cross-sector λ_C** unification: **DERIVED** (GL form factor 165/167)
> - **Cross-sector K-taxonomy**: **DERIVED** (chirality-counting per sektor)

---

## Cross-references

- [[README.md]] — opis folderu + YAML status
- [[NEEDS.md]] — otwarte luki tego folderu
- [[meta/research/RESEARCH_BUS.md]] — broadcast tych findings
- [[meta/research/FOLDER_STATUS_INDEX.md]] — globalna mapa