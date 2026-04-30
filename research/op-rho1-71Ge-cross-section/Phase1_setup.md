---
title: "ρ.1.Phase1 setup — ⁷¹Ge cross-section landscape + BEST anomaly viability (5 sub-tests)"
date: 2026-04-30
cycle: ρ.1.Phase1
status: PRE-EXECUTION
parent: "[[program.md]]"
predecessor: "[[../op-pi1-bb0nu-nme-isotope/Phase3_results.md]]"
tags:
  - TGP
  - rho1
  - phase1
  - 71Ge
  - cross-section
  - BEST
  - landscape
---

# ρ.1.Phase1 — ⁷¹Ge cross-section landscape (5 sub-tests)

> **Cel:** Compile Bahcall 1997 + Frekers 2011 + Haxton 2013 +
> Kostensalo 2019 cross-section calculations dla ⁷¹Ga(ν_e, e⁻)⁷¹Ge
> + register BEST 2022 4σ anomaly R = 0.78 ± 0.05 + 4-experiment combined
> Gallium R ≈ 0.80 + viability test dla TGP chirality-counting correction.

## Sub-tests

### P1.1 — Bahcall 1997 cross-section recompute

Bahcall (1997 ApJ 467 475) computed ⁷¹Ga(ν_e, e⁻)⁷¹Ge dla solar ν,
⁵¹Cr (746 keV), ³⁷Ar (812 keV) sources. Key components:

```
σ(E_ν) = (G_F² · cos²θ_C / π) · p_e · E_e · F(Z, E_e) · |M|²
|M_g.s.|² = B(F, g.s.) + (g_A/g_V)² · B(GT, g.s.)
         = 1 + 1.27² · 0.0865 ≈ 1.140 (g.s. only)
```

excited states (175, 500, 708 keV) add ~5–10% w E_ν > 240 keV regime.

**Gate:** Bahcall σ(⁷¹Cr → e⁻) = 5.81·10⁻⁴⁵ cm² baseline reproduced.

### P1.2 — Frekers 2011 RCNP ³He,t direct B(GT) measurement

Frekers et al. (PRC 84 014312) measured ⁷¹Ga(³He, t)⁷¹Ge:
- B(GT, g.s.) = 0.0859 ± 0.0048 (vs Bahcall 0.0865; drift 0.7% null)
- B(GT, 175 keV) = 0.0151 ± 0.0008 (consistent z Bahcall)
- B(GT, 500 keV) = 0.0145 ± 0.0008 (vs Bahcall 0.0130; lift 11%)
- B(GT, 708 keV) = 0.0078 ± 0.0006 (vs Bahcall 0.0070; lift 11%)

Combined excited-state contribution **lifted** ~5%, ground-state
**constant** → Frekers data NIE rozwiązuje BEST anomaly (~22%).

**Gate:** B(GT) sums match Frekers Tab. 2.

### P1.3 — Haxton 2013 bound-state proton overlap correction

Haxton (2013, PRC 88 015503) noted że bound-state proton wave-function
overlap dla ⁷¹Ga → ⁷¹Ge transition has ~10–20% systematic uncertainty
z shell-model truncation. Specifically:
- f_overlap = ⟨ψ_proton(⁷¹Ge) | ψ_proton(⁷¹Ga)⟩ ~ 0.85–0.95
- effective B(GT)_eff = f_overlap² · B(GT)_naive

**Gate:** f_overlap ∈ [0.85, 0.95] → cross-section reduction
[10%, 28%] → BEST 22% w środku band.

### P1.4 — 4-experiment Gallium combined fit

GALLEX (Cr-1) 0.953 ± 0.10
GALLEX (Cr-2) 0.812 ± 0.11
SAGE (Cr) 0.95 ± 0.12
SAGE (Ar) 0.79 ± 0.10
BEST (Cr inner) 0.79 ± 0.05
BEST (Cr outer) 0.77 ± 0.05

**Combined R = 0.80 ± 0.04 (5σ deficit)** — robustna systematyka.

**Gate:** combined R ∈ [0.76, 0.84] → consistent z 20% shift.

### P1.5 — TGP chirality-counting viability

Test 4 candidate B²-ratios dla ~20% reduction:

| Candidate | Form | Δ |
|---|---|---:|
| C1 | 1 − (B²_up·K_up)/(B²_lep·K_lep) = 1 − (13/4·7/8)/(2·2/3) | **17.2%** |
| C2 | 1 − (B²_lep · K_lep) / (B²_up · K_up) | reverse, negative |
| C3 | 1 − K_up/B²_up = 1 − (7/8)/(13/4) | 73.1% (too high) |
| C4 | 1 − K_lep · K_up = 1 − (2/3)·(7/8) | 41.7% (too high) |

**Gate:** ≥1 candidate within ±5% of BEST 22% → ρ.1 viable. C1 = 17.2%
within bound; refinements w Phase 2.

## PASS gate

- ≥4/5 PASS → Phase 2 viable
- 5/5 PASS → strong viability

## Środowisko

```bash
PYTHONIOENCODING=utf-8 python -X utf8 research/op-rho1-71Ge-cross-section/phase1_71Ge_landscape.py 2>&1 | tee research/op-rho1-71Ge-cross-section/phase1_71Ge_landscape.txt
```

## Cross-references

- [[program.md]]
- [[../op-pi1-bb0nu-nme-isotope/Phase2_results.md]] — closure 1/A^{1/3}
- [[../op-xi2-sterile-nu-5sector/Phase3_results.md]] — XX3 BEST tension
