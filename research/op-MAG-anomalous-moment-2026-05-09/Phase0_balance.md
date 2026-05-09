---
title: "Phase 0 balance sheet — op-MAG-anomalous-moment"
date: 2026-05-09
type: phase0-balance
status: WIP
parent: "[[./README.md]]"
phase: 0
tags:
  - phase0
  - balance-sheet
  - anomalous-moment
  - alpha-over-pi
  - phase-amplification
---

# Phase 0 — Balance sheet

## Cel cyklu (z scoped README)

Hipoteza autora:
> "α/π — może wynikać ze wzmacniania się fazy na granicy połączenia
> (wygenerowane pole przy przełączaniu dozwolonej topologii, może
> wzmacniać pewne obszary graniczne tak jakby były mierzone 2 razy)"

**Cel:** test czy g_e - 2 ≈ α/π (Schwinger 1948) emerguje z TGP-natywnej phase amplification at topology boundary.

## External inputs

- **Schwinger 1948:** g_e/2 - 1 = α/(2π) ≈ 1.16×10⁻³ (lowest-order QED radiative)
- **PDG α⁻¹** = 137.035999084(21)
- **PDG g_e/2** = 1.00115965218059(13) (12 sig figs)
- α/π ≈ 0.00232282
- TGP: bifurcation 2-state z N17, SU(2) double cover (z N18, N19, N21)

## Structural axioms (TGP-internal LOCKED)

- **S05 single-Φ axiom:** preserved
- **N17 saddle topology:** V(φ) = γ[φ³/3 - φ⁴/4], saddle at φ=1, separatrix
- **N18 SU(2) lift:** 720° symmetry, double cover signature
- **N19 lean direction:** Bloch sphere structure
- **N21 horizon ψ=4/3:** 2-sphere topology

## Derived outputs (cycle claims)

| Output | Status target |
|---|---|
| Boundary phase formalization | OPEN (Phase 1) |
| Berry / geometric phase computation | OPEN (Phase 1) |
| "Double counting" mechanism | OPEN (Phase 1) |
| α/π emergence (TGP-natywnie OR z input α) | OPEN (Phase 1 KEY decision point) |
| Quantitative reproduction | OPEN (Phase 1+) |

## 8/8 gate criteria check

✓ (1) Concrete claim z falsifiable consequences
✓ (2) Phase 0 balance (this)
TODO (3) NEEDS document → see [[./NEEDS.md]]
✓ (4) Anchor independence (mostly TGP-internal LOCKED)
✓ (5) Falsifiability test (PDG g_e/2 12 sig figs)
TODO (6) Tautology check (Phase 6)
TODO (7) Sympy verification (Phase 1)
✓ (8) Cross-validate independent path (vs standard QED Schwinger calculation)

**Status:** 5/8 ready, 3/8 pending Phase 1.

## Strategic position - QUICK SCAN approach

**To jest najbardziej speculative cykl** z 3 follow-ups. User-rekomendowana strategy:
> "Phase 0-1 quick scan (cheap), if early signals pozytywne — kontynuuj.
>  If negative — close honest."

**Quick scan goal:** Phase 1 analytical test, jeśli mechanism nie weryfikuje się elementarnie → EARLY_HALT z honest acknowledgment.

## Probability assessment (subiektywna pre-Phase-1)

| Outcome | Prob |
|---------|------|
| Pełen DERIVED (α/π z TGP first principles) | 5-10% (very ambitious — α native required) |
| STRUCTURAL DERIVED (α input, ratio α/π structural) | 20-30% |
| ANSATZ (mathematical model OK, no derivation) | 30-40% |
| EARLY_HALT (mechanism nie weryfikuje się) | 30-40% |

**Realistic outcome:** likely EARLY_HALT lub ANSATZ. α emergence z TGP byłaby major theoretical achievement, mało prawdopodobne w quick scan.

## Plan Phase 0-1

### Phase 0 (current)
- ✓ Balance sheet (this)
- TODO: NEEDS list

### Phase 1 (quick scan, 1 sympy session)
- A1: Boundary topology formalization (N17 saddle separatrix)
- A2: Berry phase / geometric phase calc
- A3: "Double counting" mechanism test
- A4: α/π emergence check
- A5: Quick verdict (continue OR EARLY_HALT)

### Phase 2-6 (only if Phase 1 pozytywne)

## Status

**Phase 0 IN PROGRESS.** NEEDS file next.
