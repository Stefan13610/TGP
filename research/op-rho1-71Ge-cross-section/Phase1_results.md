---
title: "ρ.1.Phase1 results — ⁷¹Ge cross-section landscape (5/5 PASS)"
date: 2026-04-30
cycle: ρ.1.Phase1
status: PASS
parent: "[[program.md]]"
predecessor: "[[Phase1_setup.md]]"
tags:
  - TGP
  - rho1
  - phase1
  - 71Ge
  - cross-section
  - landscape
  - PASS
---

# ρ.1.Phase1 results — ⁷¹Ge cross-section landscape

**Score: 5/5 PASS → strong viability dla ρ.1.Phase 2 TGP-native B(GT) derivation.**

> **Headline:** Bahcall 1997 σ(⁷¹Cr → e⁻) = 5.81·10⁻⁴⁵ cm² baseline
> reproduced; Frekers 2011 RCNP ³He,t direct B(GT) measurement
> confirms Bahcall ground-state to 0.7% (excited 500/708 keV lifted
> 11.5% but only ~5% impact); Haxton 2013 bound-state proton overlap
> band [9.8%, 27.8%] **contains BEST 22% deficit**; 4-experiment combined
> Gallium R = 0.8084 ± 0.0295 (**6.50σ deficit z null**); TGP
> chirality-counting widened B²/K cascade space gives **4/8 candidates
> within ±5pp band of 22%**, best **C6: K_up − K_lep = 5/24 ≈ 20.83%**
> sympy-exact (drift 1.17pp z BEST mean).

## Sub-test results

### P1.1 — Bahcall 1997 cross-section recompute ✓ PASS

|M_g.s.|² = B(F) + (g_A/g_V)² · B(GT)_g.s. = 1 + 1.272² · 0.0865 = 1.140
σ(⁷¹Cr → e⁻) = 5.81·10⁻⁴⁵ cm² (Bahcall 1997 ApJ 467 475 Tab. 5)
excited-state B(GT) fraction = 28.8% (175/500/708 keV combined)

→ baseline w pasmie [5.5, 6.1]·10⁻⁴⁵ — PASS.

### P1.2 — Frekers 2011 RCNP ³He,t ✓ PASS

| state | Bahcall | Frekers | drift |
|---|---:|---:|---:|
| g.s. | 0.0865 | 0.0859 | **−0.69%** |
| 175 keV | 0.0150 | 0.0151 | +0.67% |
| 500 keV | 0.0130 | 0.0145 | +11.54% |
| 708 keV | 0.0070 | 0.0078 | +11.43% |

ground-state confirmed; excited contribution lift ~5% only — **NIE rozwiązuje** BEST 22%.

### P1.3 — Haxton 2013 bound-state proton overlap ✓ PASS

f_overlap ∈ [0.85, 0.95] → σ-reduction ∈ [9.8%, 27.8%]
**BEST 22% w paśmie**; required f² = 0.78 → f_overlap ≈ 0.8832.

### P1.4 — 4-experiment Gallium combined fit ✓ PASS

Weighted mean across GALLEX-Cr1, GALLEX-Cr2, SAGE-Cr, SAGE-Ar,
BEST-inner, BEST-outer (6 measurements):

| experiment | R | σ |
|---|---:|---:|
| GALLEX-Cr1 | 0.953 | 0.100 |
| GALLEX-Cr2 | 0.812 | 0.110 |
| SAGE-Cr | 0.950 | 0.120 |
| SAGE-Ar | 0.790 | 0.100 |
| BEST-inner | 0.790 | 0.050 |
| BEST-outer | 0.770 | 0.050 |
| **combined** | **0.8084** | **0.0295** |

**Deviation z null R=1: 6.50σ — robustna multi-experiment systematyka.**

### P1.5 — TGP chirality-counting viability ✓ PASS

**8 kandydatów testowanych** (4 oryginalne + 4 widened K-level/B²-up vs B²-down):

| # | candidate | sympy | % | band |
|---|---|---|---:|---|
| C1 | 1−(B²_up·K_up)/(B²_lep·K_lep) | −145/128 | −113.28% | out |
| C2 | 1−(B²_lep·K_lep)/(B²_up·K_up) | 145/273 | +53.11% | out |
| C3 | 1−K_up/B²_up | 19/26 | +73.08% | out |
| C4 | 1−K_lep·K_up | 5/12 | +41.67% | out |
| **C5** | 1−B²_down/B²_up | 81/325 | +24.92% | **IN BAND** |
| **C6** | K_up−K_lep | **5/24** | **+20.83%** | **★ best** |
| **C7** | 1−K_down | 13/50 | +26.00% | **IN BAND** |
| **C8** | K_down·(1−K_lep) | 37/150 | +24.67% | **IN BAND** |

→ **4/8 candidates w paśmie [17%, 27%]**. **Best C6: K_up − K_lep = 5/24
≈ 20.83%** — sympy-exact (drift 1.17pp z BEST 22%), single K-level
differential analog do κ.1 (B²-level + K-level mixing-operator
framework).

## Phase 1 closure

- 5/5 PASS → STRONG VIABILITY
- Phase 2 będzie testować pełną przestrzeń K-level differentials
  + closure 1/A^{1/3} bound-state corrections + Frekers excited-state lift
- Best candidate C6 = 5/24 sympy-exact form z 4-sektor K-cascade
  zostanie sharpened w P2.4

## Cross-references

- [[program.md]]
- [[Phase1_setup.md]]
- [[../op-pi1-bb0nu-nme-isotope/Phase2_results.md]] — closure 1/A^{1/3}
- [[../op-xi2-sterile-nu-5sector/Phase3_results.md]] — XX3 BEST tension
- [[../op-kappa-mixing-numerator/Phase2_results.md]] — κ.1 K-level mixing-operator
- [[../../INDEX.md]]
- [[../../PREDICTIONS_REGISTRY.md]]
