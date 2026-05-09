---
title: "Phase 0 balance sheet — op-MAG-Lorentz-A-mu-coupling"
date: 2026-05-09
type: phase0-balance
status: WIP
parent: "[[./README.md]]"
phase: 0
tags:
  - phase0
  - balance-sheet
  - gate-criteria
  - lorentz-force
  - A-mu-coupling
---

# Phase 0 — Balance sheet

## Cel cyklu

Per [[./README.md]]:
1. EM-strength Lorentz F = qv × B w TGP framework (compatibility + derivation)
2. Spinor amplification mechanism (gravitomagnetic ~10⁻⁴⁴ → EM-strength)
3. Maxwell-equivalent equations (M2: ∇·B=0, M3: Faraday) z TGP
4. Hierarchy problem (czemu α >> G_N)

## External inputs

- **PDG α_em⁻¹** = 137.035999084(21) [9 sig figs]
- **PDG μ_B** = 9.2740100783(28) × 10⁻²⁴ J/T (Bohr magneton)
- **CODATA G_N** = 6.67430(15) × 10⁻¹¹ m³/(kg·s²)
- **PDG m_e** = 0.51099895000(15) MeV
- **PDG e** (elementary charge) = 1.602176634 × 10⁻¹⁹ C (exact)
- **Standard QED:** Pauli equation (z Dirac NR limit), minimum coupling p → p - qA

## Structural axioms (TGP-internal LOCKED)

- **S05 single-Φ axiom:** zachowane
- **Stage 2 photon = standard A_μ:** binding (z op-Phi-decomposition-photon)
- **N17/N18 spinor structure:** DERIVED (z op-SPIN-SU2-substrate-derivation)
- **N2d gravitomagnetic Biot-Savart:** B_g = (3/2)(G m₂/c²)(v₂ × r̂)/r² (z op-MAG-resonance-formalization)
- **M4 g_e = 2 leading order:** DERIVED (z op-MAG-resonance-formalization)
- **M9.1''(Φ̄):** existing TGP result (used dla compatibility)

## Derived outputs (cycle claims)

| Output | Status target |
|---|---|
| Pauli equation Lorentz F = qE + qv×B emergence | DERIVED (compatibility z standard QED) |
| ∇·B = 0 (M2) | DERIVED (Bianchi identity, A_μ structure) |
| Faraday ∇×E = -∂B/∂t (M3) | DERIVED (z A_μ Maxwell sector) |
| Spinor amplification (gravitomagnetic → EM) | OPEN (key novel question) |
| Hierarchy α/G_N ratio explanation | OPEN |
| Magnetic moment μ_B = (q/2m)·g·S | DERIVED (extension M4) |

## 8/8 gate criteria check

### (1) Concrete claim z falsifiable consequences
✓ Pauli equation Lorentz force structure jest **falsifiable** (cyclotron motion measurements, magnetic moment μ_B)

### (2) Phase 0 balance sheet
✓ Niniejszy plik

### (3) NEEDS document
TODO — patrz [[./NEEDS.md]] (do utworzenia)

### (4) Anchor independence
✓ Anchory są:
- Standard QED Pauli equation (external math)
- N17/N18 spinor (TGP-internal LOCKED)
- N2d gravitomagnetic Biot-Savart (TGP-internal LOCKED)
- Stage 2 A_μ (TGP-internal LOCKED)

Brak self-reference.

### (5) Falsifiability test (1σ)
✓ Standard tests:
- Cyclotron frequency ω_c = qB/m (precyzja ~10⁻⁹)
- Magnetic moment electron g_e/2 (12 sig figs)
- Lorentz invariance tests (CPT)

### (6) Tautology check
TODO — przeprowadzić w Phase 6

### (7) Sympy verification
TODO — Phase 1+

### (8) Cross-validate independent path
✓ Mamy mocne TGP-internal anchors plus standard QED comparison

**Status: 6/8 ready, 2/8 (NEEDS doc + sympy) TODO.**

## Falsifiability — empirical tests

Niezależne tests dla cycle claims:

1. **Cyclotron motion:** electron w B field → helix z ω_c = qB/m. Standard QED reprodukuje precyzyjnie. TGP musi.

2. **Magnetic moment electron:** g_e ≈ 2 (M4 done) + α/π anomaly (out of scope, separate cycle)

3. **Lorentz transformations:** all moving observers see same physics. TGP musi preservować.

4. **Faraday induction:** dB/dt → induced E. Standard Maxwell. TGP compatibility.

5. **Hall effect:** carriers in B field → transverse voltage. Lorentz force test.

6. **Aharonov-Bohm:** phase shift z A bez B (locally). QED gauge structure. TGP test.

## Probability assessment (subiektywna pre-cycle)

Per scoped README:
| Outcome | Prob |
|---------|------|
| Pełen DERIVED (TGP-natywne EM-strength) | 15-25% |
| STRUCTURAL CONDITIONAL (compatibility + spinor amplification) | 40-50% |
| MOSTLY STANDARD (TGP nie adds beyond compatibility) | 25-35% |
| EARLY_HALT (mechanism nie weryfikuje się) | 5-10% |

**Most realistic outcome:** STRUCTURAL CONDITIONAL — TGP framework jest **compatible** z Lorentz/Maxwell, spinor amplification mechanism partially clarified, ale hierarchy α/G_N pozostaje open.

## Strategic position

Cykl jest **medium ambition**:
- Most requirements (M1-M3) są standard QED — TGP musi tylko być **compatible**
- Spinor amplification jest **novel TGP test** — main payoff
- Hierarchy α/G_N może łączyć się z op-Phi-vacuum-scale (FUTURE)

**Plan time:** 3-6 weeks (per scoped README).

## Plan Phase 1+

| Phase | Cel | NEEDS |
|---|---|---|
| **Phase 1** | Pauli Lorentz F = qv × B derivation w TGP context | L1, L2 |
| **Phase 2** | Maxwell sector (∇·B=0, Faraday) z A_μ structure | L3, L4 |
| **Phase 3** | Spinor amplification mechanism (gravitomagnetic → EM) | L5, L6 |
| **Phase 4** | Hierarchy α/G_N analysis | L7, L8 |
| **Phase 5** | Aharonov-Bohm + advanced gauge tests | L9 |
| **Phase 6** | ABSOLUTE BINDING gate | — |

## Status

**Phase 0 IN PROGRESS.** Cel: complete NEEDS document, then Phase 1.
