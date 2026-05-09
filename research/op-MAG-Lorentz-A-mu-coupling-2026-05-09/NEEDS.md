---
title: "NEEDS — op-MAG-Lorentz-A-mu-coupling"
date: 2026-05-09
type: needs-list
status: WIP
parent: "[[./README.md]]"
tags:
  - needs
  - phase0
  - lorentz
  - A-mu-coupling
---

# NEEDS — op-MAG-Lorentz-A-mu-coupling

## Format

Każdy NEED ma: ID (Lx), Priority, Phase, Status, Description, Resolution path.

## CRITICAL needs

### L1: Pauli equation Lorentz force structure
**Priority:** CRITICAL
**Phase:** 1
**Status:** OPEN

**Description:** Pokazanie że Pauli equation z TGP spinor (z N18) + standard A_μ
coupling (z Stage 2) daje:
- F = qE + qv × B (classical limit)
- Cyclotron motion ω_c = qB/m
- Magnetic moment μ_B z Pauli term σ·B

**Resolution path:**
- Pauli H = (p - qA)²/(2m) + V - q·σ·B/(2m)
- Heisenberg EOM: dp/dt = -∇H w external A
- Result: F = qE + qv × B

**Bez tego:** Phase 1 fail — niemożliwe pokazać EM-strength Lorentz w TGP.

### L2: Sympy verification Pauli derivation
**Priority:** CRITICAL
**Phase:** 1
**Status:** OPEN

**Description:** Symbolic derivation Lorentz force z Pauli Hamiltonian.

**Resolution path:**
- Sympy: build Pauli H, compute commutators [r, H], [p, H]
- Show v = [r, H]/iℏ, F = [p, H]/iℏ
- Recover qE + qv × B + magnetic moment correction

### L3: ∇·B = 0 (M2) z A_μ structure
**Priority:** CRITICAL
**Phase:** 2
**Status:** OPEN

**Description:** Bianchi identity ∇·(∇×A) = 0 → ∇·B = 0 trivially.
Ale w TGP context: czy A_μ struktura **wymagana** jest TGP-natywnie, czy tylko external?

**Resolution path:**
- Stage 2 result: photon = A_μ jest TGP-compatible
- Vector potential w TGP: jak relates do scalar Φ?
- Sympy: explicitly verify ∇·(∇×A) = 0

### L4: Faraday ∇×E = -∂B/∂t (M3)
**Priority:** CRITICAL
**Phase:** 2
**Status:** OPEN

**Description:** Faraday law z A_μ:
  E = -∇φ - ∂A/∂t
  B = ∇×A
  → ∇×E = -∇×∇φ - ∇×∂A/∂t = 0 - ∂(∇×A)/∂t = -∂B/∂t  ✓

Trivial z A_μ structure. Verify w TGP context.

**Resolution path:**
- Sympy: derive Faraday z A_μ definitions
- Confirm consistency z Stage 2

## IMPORTANT needs (novel TGP content)

### L5: Spinor amplification mechanism (KEY NOVEL)
**Priority:** IMPORTANT
**Phase:** 3
**Status:** OPEN

**Description:** Czy spinor S coupling jednocześnie z:
- Gravitomagnetic B_g (z N2d, ~10⁻⁴⁴ × EM)
- Standard EM A_μ (z Stage 2)

generuje effective coupling który "amplifies" gravitomagnetic do EM-strength?

**Hypothesis (z scoped README):**
- Geometric/topological enhancement z N18 SU(2) double cover structure
- Effective coupling constant różny od G_N
- Connection przez α (fine structure constant)

**Test:**
- Compute total magnetic moment μ z BOTH couplings
- μ_total = μ_EM (z A_μ) + μ_grav (z B_g)?
- μ_EM ~ μ_B (Bohr magneton)
- μ_grav ~ G m → ~10⁻⁴⁴ × μ_B
- W single spinor S, oba contributions add coherently?

**Resolution path:**
- Sympy: extended Pauli H z BOTH B and B_g terms
- Compute g-factor and magnetic moment
- Identify possible amplification mechanism

### L6: μ_B reproduction
**Priority:** IMPORTANT
**Phase:** 3
**Status:** OPEN (extends M4)

**Description:** Pełen μ_B = eℏ/(2m_e) z TGP framework:
- M4 done: g_e = 2 leading
- L6: μ = (eℏ/2m_e)·σ z spinor + A_μ coupling

**Resolution path:**
- Pauli H term: -μ·B = -(q/2m)·σ·B (already derived w M4)
- Numerical match z μ_B PDG value
- Trivial after M4

### L7: Hierarchy α/G_N analysis
**Priority:** IMPORTANT (BUT challenging)
**Phase:** 4
**Status:** OPEN

**Description:** Dlaczego α ≈ 1/137 vs G_N m_e²/(ℏ c) ≈ 10⁻⁴³?
- Ratio α/(G m²) ~ 10⁴⁴ (Planck mass squared problem)
- Czy TGP daje resolution lub tylko inherits standard hierarchy problem?

**Resolution path:**
- Inspect coupling constants w TGP framework
- Φ_0 dependence (cross-cycle z op-Phi-vacuum-scale)
- Possibly STRUCTURAL CONDITIONAL outcome (open question)

### L8: Hall effect / cyclotron precision tests
**Priority:** IMPORTANT
**Phase:** 4-5
**Status:** OPEN

**Description:** TGP musi reprodukować precyzyjnie:
- Cyclotron ω_c = qB/m (10⁻⁹ precision in atomic clocks)
- Hall coefficient
- Lorentz invariance

**Resolution path:**
- Derive z Pauli equation
- Show identity z standard QED prediction
- TGP nie odbiega na poziomie 10⁻⁶

## NICE-TO-HAVE needs

### L9: Aharonov-Bohm effect
**Priority:** NICE-TO-HAVE
**Phase:** 5
**Status:** OPEN

**Description:** Phase shift z A bez B (locally) — gauge structure test.

**Resolution path:**
- Standard QED reproduction
- TGP compatibility check (no contradiction)

### L10: Connection do Stage 2 photon ontology
**Priority:** NICE-TO-HAVE
**Phase:** 1-2
**Status:** OPEN

**Description:** Stage 2 ustaliło photon = standard A_μ. Niniejszy cykl używa A_μ
jako external. Ale **dlaczego** A_μ struktura jest TGP-natywna (NIE δΦ-mode)?

**Resolution path:**
- Reference Stage 2 derivation
- Confirm consistency z TGP single-Φ axiom (S05)
- Show A_μ jest **complementary** field NIE wbrew S05

## Summary by phase

| Phase | Critical | Important | Nice |
|-------|---------|-----------|------|
| 1 (Pauli Lorentz) | L1, L2 | — | L10 |
| 2 (Maxwell) | L3, L4 | — | — |
| 3 (Amplification) | — | L5, L6 | — |
| 4 (Hierarchy) | — | L7 | — |
| 5 (Advanced tests) | — | L8 | L9 |

## Decision criteria

**Phase 0 → Phase 1:**
- ☑ Phase 0 balance sheet
- ☑ NEEDS document (niniejszy)
- ☐ User confirmation start

**Phase 1 → Phase 2:**
- L1, L2 RESOLVED z DERIVED
- Pauli Lorentz force confirmed compatible

**Phase 2 → Phase 3:**
- L3, L4 RESOLVED z DERIVED
- Maxwell sector compatibility verified

**Phase 3 → Phase 4 (lub close):**
- L5 (amplification) — **KEY decision point**:
  - Mechanism wyłania się TGP-natywnie → continue Phase 4 (hierarchy)
  - No clear mechanism → STRUCTURAL CONDITIONAL classify, close cycle
  - Ambiguous → continue z honest acknowledgment

## Honest reporting note

Cycle ma **dwie warstwy**:
1. **Compatibility check** (M1-M3) — likely DERIVED easy
2. **Novel test** (L5 amplification, L7 hierarchy) — UNCERTAIN

Cycle outcome zależy głównie od L5 result.

**Status NEEDS:** WIP, początkowa wersja 2026-05-09.
