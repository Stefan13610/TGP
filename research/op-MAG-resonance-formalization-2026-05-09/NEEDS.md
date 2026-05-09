---
title: "NEEDS — op-MAG-resonance-formalization"
date: 2026-05-09
type: needs-list
status: WIP
parent: "[[./README.md]]"
tags:
  - needs
  - phase0
  - magnetism
  - resonance
---

# NEEDS — explicit wymagania dla cyklu op-MAG-resonance-formalization

## CRITICAL needs (5)

### N1: Definicja "częstotliwości ω" dla TGP soliton
**Priority:** CRITICAL — **FOUNDATIONAL BLOCKER**
**Phase:** 1 (immediate, primary task)
**Status:** OPEN

**Description:** Centralny concept "resonance" wymaga precyzyjnej
definicji "częstotliwości" dla solitonu. Możliwe sources:
- (a) Linear oscillation frequency around meta-stable saddle: ω_lin ~ √γ
- (b) Bifurcation rate (zanik vs ekspansja): ω_bif
- (c) Internal "blinking" frequency (z working hypothesis): ω_blink
- (d) Compton-like frequency m·c²/ℏ (związane z mass)

**Resolution path:**
- Analytical analysis każdej opcji
- Test consistency z QM gyromagnetic results
- Test consistency z observed cyclotron resonance
- Wybór operacjonalna definicji

**Bez tego:** cały cykl nie może procedować — "resonance" jest pustym
sloganem bez ω.

### N2: Coupling Hamiltonian formalization
**Priority:** CRITICAL
**Phase:** 1
**Status:** OPEN

**Description:** Sprzężenie soliton ↔ δΦ field musi mieć konkretną
mathematical structure:
```
H_coupling = ∫ d³x F[δΦ_soliton(x), δΦ_field(x)]
```

Pytania:
- Linear coupling F ~ δΦ_1 · δΦ_2?
- Quadratic / nonlinear F?
- Gauge / Lorentz invariance?

**Resolution path:**
- Variational principle z TGP action
- Effective field theory derivation
- Sympy + analytical work

### N3: Lorentz force F = qv × B derivation
**Priority:** CRITICAL — **MAIN TEST cyklu**
**Phase:** 2
**Status:** OPEN

**Description:** Pokazanie że limit slow-motion daje
F = q(E + v × B) z δΦ-resonance Hamiltonianu.

**Resolution path:**
- Static limit: gradient force from H_coupling
- Slow-motion expansion: leading correction
- Comparison z QED gold standard
- Sympy verification dla simple test cases (point charge, uniform B)

**Decyzja kluczowa:** jeśli derivation NIE daje F = qv × B z dokładnością
~10⁻³, cykl pivots do STRUCTURAL_NO_GO M1.

### N4: Maxwell-equivalent equations (∇·B=0, Faraday)
**Priority:** CRITICAL
**Phase:** 3
**Status:** OPEN

**Description:** Pokazanie że TGP-natywne formalism daje (lub jest
spójny z):
- ∇·B = 0 (no monopoles)
- ∇×E = -∂B/∂t (Faraday)

**Resolution path:**
- Derivation z structure δΦ-coupling
- Bianchi-like identities w TGP
- Compatibility z M9.1'' i Stage 2 photon = A_μ

### N5: g_e ≈ 2 z TGP dynamics
**Priority:** CRITICAL
**Phase:** 4
**Status:** OPEN

**Description:** Reprodukcja anomalous magnetic moment:
- Leading order: g_e = 2 (Dirac result)
- Anomaly: g_e - 2 ≈ α/π (Schwinger)

**Resolution path:**
- Coupling spinor (z op-SPIN-SU2 N18) z δΦ field
- Static + dynamic ("blinking") contributions
- Radiative corrections z Φ̄ vacuum fluctuations

## IMPORTANT needs (4)

### N6: Mach inertia formula testing
**Priority:** IMPORTANT (najambitniej)
**Phase:** 5
**Status:** OPEN

**Description:** Quantitative formula:
```
m_eff ~ ⟨δΦ²_background⟩ · ∫ δΦ²_soliton d³x
```

Test z m_e ~ 511 keV przy reasonable cosmological background.

**Resolution path:**
- Background δΦ amplitude estimate (z FRW solutions, EXT-1 context)
- Soliton intensity integral (od op-SPIN-SU2 jeśli orientation
  resolved)
- Numerical comparison

### N7: Resonance selectivity test
**Priority:** IMPORTANT
**Phase:** 1-2
**Status:** OPEN

**Description:** Pokazanie że frequency-matching daje observed
selectivity:
- Electrons (same m, same ω) coupling jednolity z B
- Nuclei (różne m, różne ω) coupling słabszy lub inny
- Photons (massless?) coupling dynamiczny

### N8: Lorentz invariance preservation
**Priority:** IMPORTANT
**Phase:** 2-3
**Status:** OPEN

**Description:** Pełen framework musi być Lorentz invariant
(modulo M9.1'' background corrections). Cyclotron motion w boost'd
frame should give consistent result.

### N9: Photon emission/absorption mechanism
**Priority:** IMPORTANT
**Phase:** 4
**Status:** OPEN

**Description:** Soliton coupling do A_μ field przez δΦ-resonance.
Photon emission jako relaxation procces w resonance Hamiltonianie.

## NICE-TO-HAVE needs (3)

### N10: Connection do L01 ρ_EM audit
**Priority:** NICE-TO-HAVE
**Phase:** 3-4
**Status:** OPEN

**Description:** Czy niniejszy cykl daje natywny ρ_EM mechanism który
może rozwiązać L01 audit OPEN?

### N11: Gauge invariance verification
**Priority:** NICE-TO-HAVE
**Phase:** 3
**Status:** OPEN

**Description:** Standard EM ma U(1) gauge symmetry. Czy TGP
δΦ-resonance ma analogiczną symmetry?

### N12: Quantum Hall, BCS, etc. predictions
**Priority:** NICE-TO-HAVE
**Phase:** 6
**Status:** OPEN

**Description:** Może niniejszy framework dać natywne predictions
dla Quantum Hall, superconductivity (BCS), etc.

## Summary by phase

| Phase | Critical needs | Important | Nice |
|-------|---------------|-----------|------|
| 1 | **N1** (BLOCKER), N2 | N7 | — |
| 2 | **N3** (MAIN TEST) | N7, N8 | — |
| 3 | N4 | N8 | N10, N11 |
| 4 | N5 | N9 | — |
| 5 | — | N6 | — |
| 6 | — | — | N12 |

**N1 jest FOUNDATIONAL BLOCKER** — bez precyzyjnej definicji "ω"
cały cykl jest pustym wezwaniem. Phase 1 musi to rozstrzygnąć
przed dalszą pracą.

## Decision criteria — Phase 0 → Phase 1

**Wymagane przed startem Phase 1:**
- ☑ Phase 0 balance sheet completed (gate criteria 7-8/8)
- ☑ NEEDS document complete (niniejszy plik)
- ☐ User confirmation: Phase 1 priority "resonance ω definition"
- ☐ Acknowledgment długoterminowej skali (~6-12 miesięcy)

**Critical decision points w cyklu:**

1. **Po N1** (Phase 1): czy ω jest dobrze zdefiniowane?
   - TAK → Phase 2 derivation Lorentz force
   - NIE → cykl pivot lub STRUCTURAL_NO_GO

2. **Po N3** (Phase 2): czy F = qv × B jest reproducted?
   - TAK → Phase 3+
   - PARTIAL → STRUCTURAL CONDITIONAL
   - NIE → STRUCTURAL_NO_GO M1

3. **Po N5** (Phase 4): czy g_e ≈ 2 z TGP dynamics?
   - TAK z dokładnością ~10⁻³ → DERIVED M4
   - NIE → STRUCTURAL CONDITIONAL z gap

## Honest reporting note

Lista NEEDS jest **wstępna** (Phase 0). Może rosnąć w Phase 1+ jak
ujawnia się więcej technicznych wymagań. Każdy nowy NEED dodawany
incrementalnie z timestamp i justification.

**Status NEEDS document:** WIP, początkowa wersja 2026-05-09.

**Cykl scope:** to jest **najambitniejszy cykl** w TGP framework do
tej pory. Pełen DERIVED probability ~15-25%; STRUCTURAL CONDITIONAL
probability ~35-45%. Honest expectation: significant partial success
najprawdopodobniejszy outcome.
