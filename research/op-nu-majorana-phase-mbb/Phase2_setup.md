---
title: "ν.1.Phase2 setup — first-principles α₂₁/α₃₁ + m_ββ TGP NO Form A/B (7 sub-tests)"
date: 2026-04-30
cycle: ν.1.Phase2
status: PRE-EXECUTION
parent: "[[program.md]]"
predecessor: "[[Phase1_results.md]]"
tags:
  - TGP
  - nu1
  - phase2
  - majorana-phase
  - mbb
  - first-principles
---

# ν.1.Phase2 — first-principles α₂₁/α₃₁ + m_ββ TGP NO Form A/B (7 sub-tests)

> **Cel:** Sympy-exact derivacja Majorana phases (α₂₁_A, α₃₁_A) z TGP
> 4-sector chirality-halving + (α₂₁_B, α₃₁_B) z PMNS-Wolfenstein analog;
> compute m_ββ_TGP NO ordering dla obu form pair; falsify 5 alternative
> Majorana phase parametrizations + classification cascade (5 promotions).

## Sub-tests

### N2.1 — α₂₁_A chirality-halving sympy-exact

$$\alpha_{21,A} = \pi \cdot \frac{B^2_{lep} - B^2_{\nu}}{B^2_{lep}} = \pi \cdot \frac{2 - 1}{2} = \frac{\pi}{2}$$

Anchor: Majorana B²_ν=1 (single-chirality counting) vs Dirac B²_lep=2
(chirality-doubled). Half-period factor 1/2 odzwierciedla loss of one
chirality channel (Majorana = self-conjugate).

### N2.2 — α₃₁_A (ν,up) pair sympy-exact

$$\alpha_{31,A} = 2\pi \cdot \frac{B^2_{up} - B^2_{\nu}}{B^2_{up\_num}} = 2\pi \cdot \frac{13/4 - 1}{13} = 2\pi \cdot \frac{9/4}{13} = \frac{9\pi}{26}$$

Anchor: pair (ν, up) coupling — up B²_up=13/4 (Dirac+QCD) minus ν
B²_ν=1 (Majorana-naked) gives 9/4; normalization B²_up_num=13 (numerator
of B²_up=13/4 Rational).

### N2.3 — α₂₁_B PMNS-Wolfenstein sympy-exact

$$\alpha_{21,B} = \pi \cdot (1 - \bar{\rho}_{PMNS}) = \pi \cdot (1 - 2/13) = \frac{11\pi}{13}$$

Anchor: μ.1 ρ̄_PMNS = 2/13 (PMNS-Wolfenstein analog, μ.1.Phase2.M2);
"1 − ρ̄" gives complement-from-unit-circle struktura (analog do
δ_CP = π + arctan(η̄/ρ̄)).

### N2.4 — α₃₁_B PMNS-Wolfenstein sympy-exact

$$\alpha_{31,B} = 2\pi \cdot \bar{\eta}_{PMNS} = 2\pi \cdot \frac{6}{7} = \frac{12\pi}{7} \pmod{2\pi}$$

Anchor: μ.1 η̄_PMNS = 6/7 (PMNS-Wolfenstein analog, μ.1.Phase2.M3);
2π · η̄ daje full-period z η̄ ∈ (0, 1) → α₃₁ ∈ (0, 2π).

### N2.5 — m_ββ_TGP Form A pair (α A × δ_CP B)

m_ββ_A = m_ββ(α₂₁=π/2, α₃₁=9π/26, δ_CP=π+arctan(39/7))

Predicted: ~1.5–1.7 meV (z Phase 1 = 1.584 meV).

### N2.6 — m_ββ_TGP Form B pair (α B × δ_CP B)

m_ββ_B = m_ββ(α₂₁=11π/13, α₃₁=12π/7, δ_CP=π+arctan(39/7))

Predicted: ~3.0–3.3 meV (z Phase 1 = 3.249 meV).

### N2.7 — 5 alternative phase forms FALSIFIED + classification cascade

| Alt form | (α₂₁, α₃₁) | Anchor | Falsification reason |
|----------|-----------|--------|---------------------|
| TBM trivial | (0, 0) | tri-bimaximal no-CP | violates TGP B² chirality structure |
| BM CP-π | (π, π) | bimaximal CP-conserving | violates TGP (ν,up) pair structure |
| Golden ratio | (2π/φ, 2π/φ²) | golden ratio symmetry | no TGP B² anchor |
| Hexagonal | (π/3, 2π/3) | hexagonal A₄ | no TGP chirality anchor |
| Democratic strict | (2π/3, 4π/3) | S₃ strict (z form D) | no TGP B²/PMNS anchor (different from D as candidate) |

Classification cascade (5 promotions):
- α₂₁_A: π · (B²_lep−B²_ν)/B²_lep — promoted to **DERIVED dual** (chirality-halving)
- α₃₁_A: 2π · (B²_up−B²_ν)/B²_up_num — promoted to **DERIVED dual** ((ν,up) pair)
- α₂₁_B: π · (1−ρ̄_PMNS) — promoted to **DERIVED dual** (PMNS-Wolfenstein)
- α₃₁_B: 2π · η̄_PMNS — promoted to **DERIVED dual** (PMNS-Wolfenstein)
- m_ββ_TGP NO: dual prediction ~1.6/3.2 meV — promoted to **DERIVED**

## Środowisko

```bash
PYTHONIOENCODING=utf-8 python -X utf8 research/op-nu-majorana-phase-mbb/phase2_majorana_derivation.py 2>&1 | tee research/op-nu-majorana-phase-mbb/phase2_majorana_derivation.txt
```

## Cross-references

- [[program.md]]
- [[Phase1_results.md]]
- [[../op-mu-pmns-phase-hardening/Phase2_results.md]]
- [[../op-zeta-mass-spectrum/Phase3_results.md]]
