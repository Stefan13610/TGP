---
title: "ν.1.Phase1 setup — m_ββ landscape audit + Majorana phase candidates"
date: 2026-04-30
cycle: ν.1.Phase1
status: PRE-EXECUTION
parent: "[[program.md]]"
tags:
  - TGP
  - nu1
  - phase1
  - landscape
  - mbb
---

# ν.1.Phase1 — m_ββ landscape audit + Majorana phase candidates (5 sub-tests)

> **Cel:** Map 0νββ effective mass landscape (current limits + future
> sensitivity), enumerate 4 Majorana phase candidates, register TGP NO
> inputs (ζ.1 masses + μ.1 angles), audit drift sources, viability gate.

## Sub-tests

### N1.1 — m_ββ experimental landscape

Current limits:
- KamLAND-Zen 2024 (Xe-136): m_ββ < 36–122 meV (90% CL, NME-dependent)
- CUORE (Te-130): m_ββ < 75–350 meV
- GERDA + LEGEND-200 (Ge-76): m_ββ < 79–180 meV

Future sensitivity:
- KamLAND-Zen 2027+: ~5 meV (Xe-136 LS)
- LEGEND-1000 2030+: ~10 meV → 3 meV (Ge-76, ton-scale)
- nEXO 2030+: ~5–10 meV (Xe-136 TPC, 5 t)
- NEXT-HD 2030+: ~10 meV → 3 meV (Xe-136 HP)

**Convergence target:** ≥3 z 4 next-gen reach m_ββ < 5 meV.

### N1.2 — Majorana phase candidates (4 forms)

| Form | (α₂₁, α₃₁) | Anchor |
|------|-----------|--------|
| A (chirality-halving) | (π/2, 9π/26) | Majorana B²_ν=1 vs Dirac B²_lep=2 + (ν,up) pair |
| B (PMNS-Wolfenstein) | (11π/13, 12π/7 mod 2π) | μ.1 (ρ̄_PMNS, η̄_PMNS) = (2/13, 6/7) |
| C (maximal CP) | (0, π) | reference: pure Dirac no Majorana phase |
| D (democratic) | (2π/3, 4π/3) | S₃ democratic permutation |

### N1.3 — TGP NO ordering inputs

Z ζ.1.Phase2 (m₁ via K(ν)=1/2 Majorana bisection):
- m₁ = 0.76 meV
- m₂ = 8.71 meV (z Δm²₂₁ = 7.42·10⁻⁵ eV²)
- m₃ = 49.53 meV (z Δm²₃₁ = 2.515·10⁻³ eV²)
- Σm_ν = 59.01 meV (DESI DR3 2027+ falsifiable)

Z μ.1.Phase2 (PMNS angles refined²):
- sin²θ₁₂ = 5149/16800 ≈ 0.30649 (drift 0.17%)
- sin²θ₂₃ = 4/7 ≈ 0.57143 (drift 0.10%)
- sin²θ₁₃ = 13627867/624000000 ≈ 0.02184 (drift 0.73%)
- δ_CP Form A = N_gen·arctan(195/77) ≈ 205.36°
- δ_CP Form B = π + arctan(39/7) ≈ 259.82°

### N1.4 — Drift sources audit

| Source | Drift impact |
|--------|-------------|
| NME uncertainty | factor ~3 (NME = 1–6 dla Xe-136) |
| m₁ tension z DESI | ζ.1 m₁=0.76 meV vs DESI DR2 limit Σ<0.072 eV margin +22% |
| sin²θ₁₂ μ.1 | 0.17% drift, negligible |
| sin²θ₁₃ μ.1 | 0.73% drift, negligible |
| δ_CP Form A vs B | structural ambiguity, both within NuFit 1σ |

### N1.5 — Viability gate

Wszystkie 4 candidate forms muszą produkować m_ββ_TGP w obecnych KamLAND-Zen
limits [0, 122] meV. Gate: **4/4 candidate forms compatible** → Phase 2 viable.

## Środowisko

```bash
PYTHONIOENCODING=utf-8 python -X utf8 research/op-nu-majorana-phase-mbb/phase1_nu_landscape.py 2>&1 | tee research/op-nu-majorana-phase-mbb/phase1_nu_landscape.txt
```

## Cross-references

- [[program.md]]
- [[../op-mu-pmns-phase-hardening/Phase3_results.md]]
- [[../op-zeta-mass-spectrum/Phase3_results.md]]
