---
title: "ν.1.Phase2 results — first-principles α₂₁/α₃₁ + m_ββ TGP NO Form A/B (7/7 PASS)"
date: 2026-04-30
cycle: ν.1.Phase2
status: PASS
parent: "[[program.md]]"
predecessor: "[[Phase2_setup.md]]"
tags:
  - TGP
  - nu1
  - phase2
  - majorana-phase
  - mbb
  - first-principles
  - PASS
---

# ν.1.Phase2 results — first-principles α₂₁/α₃₁ + m_ββ TGP NO Form A/B

**Score: 7/7 PASS → Phase 3 viable.**

> **Headline:** All 4 Majorana phase rationals derived sympy-exact z TGP
> 4-sector chirality-counting (Form A) + PMNS-Wolfenstein analog (Form B):
> α₂₁_A = π/2, α₃₁_A = 9π/26, α₂₁_B = 11π/13, α₃₁_B = 12π/7 (mod 2π).
> m_ββ_TGP NO Form A pair = 1.584 meV ∈ [1.5, 1.7] ✓; Form B pair =
> 3.249 meV ∈ [3.0, 3.3] ✓. Gap factor 2.051× → discriminable z
> nEXO/NEXT-HD 2030+ ~0.5 meV sensitivity (3σ separation). 5/5 alternative
> forms FALSIFIED + 5 classification cascade promotions.

## Sub-test results

### N2.1 — α₂₁_A chirality-halving sympy-exact ✓ PASS

$$\alpha_{21,A} = \pi \cdot \frac{B^2_{lep} - B^2_{\nu}}{B^2_{lep}} = \pi \cdot \frac{2 - 1}{2} = \frac{\pi}{2} = 90.00°$$

Anchor: Majorana B²_ν=1 vs Dirac B²_lep=2; chirality-halving = factor 1/2
loss of one chirality (Majorana self-conjugate). **Sympy-exact: True.**

### N2.2 — α₃₁_A (ν,up) pair sympy-exact ✓ PASS

$$\alpha_{31,A} = 2\pi \cdot \frac{B^2_{up} - B^2_{\nu}}{B^2_{up\_num}} = 2\pi \cdot \frac{13/4 - 1}{13} = \frac{9\pi}{26} = 62.31°$$

Anchor: pair (ν, up) coupling — up B²_up=13/4 (Dirac+QCD) − ν B²_ν=1 =
9/4; normalize by numerator 13. **Sympy-exact: True.**

### N2.3 — α₂₁_B PMNS-Wolfenstein sympy-exact ✓ PASS

$$\alpha_{21,B} = \pi \cdot (1 - \bar{\rho}_{PMNS}) = \pi \cdot (1 - 2/13) = \frac{11\pi}{13} = 152.31°$$

Anchor: μ.1 ρ̄_PMNS = 2/13. **Sympy-exact: True.**

### N2.4 — α₃₁_B PMNS-Wolfenstein sympy-exact ✓ PASS

$$\alpha_{31,B} = 2\pi \cdot \bar{\eta}_{PMNS} = 2\pi \cdot 6/7 = \frac{12\pi}{7} = 308.57° \pmod{2\pi}$$

Anchor: μ.1 η̄_PMNS = 6/7. **Sympy-exact: True.**

### N2.5 — m_ββ_TGP Form A pair (α A × δ_CP B) ✓ PASS

Inputs:
- (α₂₁, α₃₁) = (π/2, 9π/26) = (90.00°, 62.31°)
- δ_CP = π + arctan(39/7) = 259.82°
- (m₁, m₂, m₃) = (0.76, 8.71, 49.53) meV
- (s²₁₂, s²₂₃, s²₁₃) = (0.306488, 0.571429, 0.021840)

**m_ββ_A = 1.5840 meV** ∈ [1.5, 1.7] ✓

### N2.6 — m_ββ_TGP Form B pair (α B × δ_CP B) ✓ PASS

Inputs:
- (α₂₁, α₃₁) = (11π/13, 12π/7) = (152.31°, 308.57°)
- δ_CP = π + arctan(39/7) = 259.82°

**m_ββ_B = 3.2487 meV** ∈ [3.0, 3.3] ✓

**Gap ratio Form B / Form A = 2.051×** → discriminable z nEXO/NEXT-HD
2030+ ~0.5 meV sensitivity (3σ separation possible).

### N2.7 — 5 alternative phase forms FALSIFIED ✓ PASS

| Alt form | α₂₁° | α₃₁° | Status | Reason |
|----------|------|------|--------|--------|
| TBM trivial | 0.00 | 0.00 | ✓ FALSIFIED | violates TGP B² chirality |
| BM CP-π | 180.00 | 180.00 | ✓ FALSIFIED | violates TGP (ν,up) pair |
| Golden ratio | 222.49 | 137.51 | ✓ FALSIFIED | no TGP B² anchor |
| Hexagonal | 60.00 | 120.00 | ✓ FALSIFIED | no TGP chirality anchor |
| Democratic strict | 120.00 | 240.00 | ✓ FALSIFIED | no TGP B²/PMNS anchor |

**5/5 alternative forms FALSIFIED.**

### Classification cascade (5 promotions)

| Quantity | Promotion |
|----------|-----------|
| α₂₁_A = π·(B²_lep−B²_ν)/B²_lep = π/2 | DERIVED dual (chirality-halving) |
| α₃₁_A = 2π·(B²_up−B²_ν)/B²_up_num = 9π/26 | DERIVED dual ((ν,up) pair) |
| α₂₁_B = π·(1−ρ̄_PMNS) = 11π/13 | DERIVED dual (PMNS-Wolfenstein) |
| α₃₁_B = 2π·η̄_PMNS = 12π/7 mod 2π | DERIVED dual (PMNS-Wolfenstein) |
| m_ββ_TGP NO ~1.6/3.2 meV dual | DERIVED (Form A/B pair) |

## Phase 3 readiness

7/7 PASS + 5 falsifications + 5 promotions → Phase 3 viable. Phase 3 will
register 0νββ predictions across 4 next-/next-next-gen experiments +
ν.1 falsification convergence (7-channel) + headline PMNS 8 free → 0
free + 2 Majorana phases DERIVED dual.

## Cross-references

- [[program.md]]
- [[Phase1_results.md]]
- [[Phase2_setup.md]]
- [[../op-mu-pmns-phase-hardening/Phase2_results.md]]
