---
title: "Phase 5 results — Singletu (N6) + no-signaling (N7) → Bell -cos(α-β) + CHSH 2√2"
date: 2026-05-09
type: phase-results
status: STRUCTURAL_DERIVED
parent: "[[./README.md]]"
phase: 5
needs_resolved: [N6, N7]
sympy_verification: 8/8 PASS
tags:
  - phase5
  - bell-correlation
  - singletu
  - no-signaling
  - CHSH
  - tsirelson-bound
---

# Phase 5 — Singletu + Bell + no-signaling

## Status: **STRUCTURAL DERIVED** — wszystkie six magnetism requirements satisfied

**Sympy:** 8/8 PASS

## Cel Phase 5

Per [[./README.md]] Phase 5 + [[./NEEDS.md]] N6, N7:
- **N6** (Phase 5): Singletu konstrukcja w sklejenie picture
- **N7** (Phase 5): No-signaling formal proof
- **Cel**: derivation E(α,β) = -cos(α-β) (Bell correlation)

## Setup

2-particle Hilbert space: (C²) ⊗ (C²) = 4-dim

Singletu state:
```
|Ψ⟩ = (1/√2)(|+−⟩ − |−+⟩)
```

W TGP-natywnej interpretacji:
- Każdy soliton: 2-state bifurcation (N17 zanik/ekspansja)
- Singletu = **wspólna geometria sklejenia** 2 solitonów (antisymmetric tensor)
- Pomiar = external field selects axis (N19 induced SU(2))

## Kluczowe wyniki

### S2-S3: Singletu właściwości

- ✓ |Ψ⟩ properly normalized: ⟨Ψ|Ψ⟩ = 1
- ✓ Antysymmetric: P₁₂|Ψ⟩ = -|Ψ⟩ (swap operator daje -1)
- ✓ Total spin: J²|Ψ⟩ = 0 (singletu = j=0 representation)
- ✓ J_z|Ψ⟩ = 0 (m=0 component)

### S4-S5: Bell correlation E(α,β) = -cos(α-β)

Setup pomiarów (axes w plane xz):
```
â = (sin α, 0, cos α)
b̂ = (sin β, 0, cos β)
```

Joint observable: (σ·â) ⊗ (σ·b̂)

Correlation:
```
E(α, β) = ⟨Ψ| (σ·â) ⊗ (σ·b̂) |Ψ⟩
```

**Sympy verified: E(α, β) = -cos(α - β)** ✓

To jest QM signature singletu — Bell correlation.

### S6: No-signaling

Marginal P(+|a) na particle 1 niezależnie od axis b:
- ✓ P(+|a) = 1/2 (rotation invariant)
- ✓ Reduced density matrix: ρ₁ = Tr₂(|Ψ⟩⟨Ψ|) = (1/2)I (maximally mixed)
- **Wniosek:** żadna preferencja kierunku na particle 1, no-signaling **formal**

### S7: CHSH violation (Tsirelson bound)

Standard test: a=0, a'=π/2, b=π/4, b'=3π/4

```
S = E(a,b) - E(a,b') + E(a',b) + E(a',b')
  = -cos(0-π/4) + cos(0-3π/4) - cos(π/2-π/4) - cos(π/2-3π/4)
  = -√2/2 + (-√2/2) + (-√2/2) + (-√2/2)
  = -2√2
|S| = 2√2 ≈ 2.828   (Tsirelson bound)
```

- ✓ Sympy verified: |S| = 2√2 exactly
- Local realism predicts |S| ≤ 2 → **violated**
- QM saturates bound: |S| = 2√2

## Six requirements: final check

| # | Requirement | Status | Source |
|---|---|---|---|
| 1 | Dwa wyniki ±1 | ✓ | Bifurcation (N17) |
| 2 | Inna struktura niż wektor R³ | ✓ | 2-spinor (z N18) |
| 3 | 720° symmetry | ✓ | U(2π) = -I (N18, N19, N21) |
| 4 | cos²(θ/2) projekcja | ✓ | 2-spinor Born rule (N18) |
| 5 | **Singletu -cos(α-β)** | ✓ | **Phase 5 (this)** |
| 6 | Brak tabeli predefined values | ✓ | CHSH violation 2√2 |

**Score: 6/6 SATISFIED.**

## Konwergencja całego cyklu

| Phase / Need | Result | Sympy |
|---|---|---|
| **N16** Dynamic equilibrium | Derrick rozbrojony strukturalnie | 3/3 |
| **N17** Bifurcation 2-state | Saddle/branches w fazoplane V(φ) | 7/7 |
| **N18** SU(2) z 2-state | Fundamental rep, Born rule, 720° | 7/7 |
| **N19** External SO(3) embedding | Induced SU(2) z lean direction | 11/11 |
| **N21** Horizon multipole l=1 | Independent SU(2) z M9.1'' | 11/11 |
| **Phase 5** Singletu + Bell | -cos(α-β), no-signaling, 2√2 | 8/8 |

**Total cycle: 47/47 sympy PASS**

## What Phase 5 delivers

✓ **N6 RESOLVED** — singletu konstrukcja z 2-particle Hilbert (antysymetryczny tensor product)
✓ **N7 RESOLVED** — no-signaling **formally** proven (reduced density matrix = (1/2)I)
✓ **Bell correlation -cos(α-β)** derived analytically
✓ **CHSH = 2√2** Tsirelson bound saturated
✓ **All six magnetism requirements** satisfied

## Limitations / honest acknowledgments

### L1: TGP-natywna interpretacja "wspólnej geometrii sklejenia"

Antysymetryczny tensor product jest formalny mathematical object. **Konkretne TGP-natywne mapping** (jak 2 solitony fizycznie wytwarzają singletu state) wymaga:
- Kink-fermion physical pair production mechanism
- Conservation laws (total angular momentum) z TGP framework
- Topological/coupling anchor

To **NIE jest** rozwiązane w Phase 5 — jest **interpretowane**.

### L2: Compatibility, not new prediction

Phase 5 reproduces standard QM Bell physics. To jest **compatibility check** — pokazanie że TGP framework jest **konsystentny** z standard Bell physics.

**NIE jest to** new prediction beyond QM. Jest to demonstracja że TGP-derived 2-spinor + 2-particle Hilbert space dają standardowe QM correlations.

### L3: Fizyczna realność singletu w TGP solitonach

Singletu z 2 solitonów (e.g., positronium-like, electron pair) wymaga produkcji w pair process. TGP framework nie specifies pair production mechanism (out of scope).

## Probability re-update

| Outcome | Pre-Phase-5 | Post-Phase-5 |
|---|---|---|
| Pełen DERIVED | 50-60% | **55-65%** ↑ |
| STRUCTURAL CONDITIONAL | 30-35% | **25-30%** ↓ |
| STRUCTURAL_NO_GO | <10% | <5% |

**Reasoning:** Phase 5 closes wszystkie six requirements. Cykl ma teraz **kompletny analytical foundation** dla TGP-natywnej spinor structure.

## Cross-references

### Within cycle
- [[./Phase5_singletu_Bell_sympy.py]] — sympy 8/8 PASS
- [[./Phase1_N16_results.md]] — dynamic equilibrium foundation
- [[./Phase1_N17_results.md]] — bifurcation 2-state
- [[./Phase1_N18_results.md]] — first SU(2) emergence
- [[./Phase2_N21_results.md]] — second SU(2) (horizon)
- [[./Phase3_N19_results.md]] — third SU(2) (external embedding)

### External
- Standard Bell test literature (Bell 1964, CHSH 1969)
- Tsirelson bound (Cirel'son 1980)

## Status

**Phase 5 RESOLVED z STRUCTURAL DERIVED.** Six requirements **all satisfied**.

**Cycle sympy total:** 47/47 PASS.

**Next:** Phase 6 ABSOLUTE BINDING gate dla cycle close.
