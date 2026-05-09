---
title: "NEEDS — op-MAG-anomalous-moment"
date: 2026-05-09
type: needs-list
status: WIP
parent: "[[./README.md]]"
tags:
  - needs
  - phase0
  - anomalous-moment
---

# NEEDS — op-MAG-anomalous-moment

## CRITICAL needs

### A1: Geometric definition "boundary" w N17 topology
**Priority:** CRITICAL
**Phase:** 1
**Status:** OPEN

**Description:** Co dokładnie znaczy "granica połączenia dozwolonej topologii" w TGP?

Kandydaci:
- Saddle separatrix V(φ)=γ[φ³/3 - φ⁴/4] (z N17): granica między basenami zanik/ekspansja
- Horizon ψ=4/3 (z N21): boundary między timelike i spacelike regions
- Bifurcation phase plane: stable/unstable manifolds of saddle

**Resolution path:**
- Sympy: phase-plane analysis V(φ), separatrix curve
- Identify "boundary measure" (length, area, etc.)

### A2: Boundary phase amplification mechanism
**Priority:** CRITICAL
**Phase:** 1
**Status:** OPEN

**Description:** Jaki mechanism amplifies phase na boundary?

Kandydaci:
- Berry phase / geometric phase z adiabatic transport
- Topological winding number
- "Double cover" structure (z N18) — boundary "measured 2x"

**Resolution path:**
- Sympy: compute Berry phase dla loop wokół saddle
- Compare z standard SU(2) Berry phase Ω/2

### A3: α emergence — TGP-natywne czy input?
**Priority:** CRITICAL (decision point)
**Phase:** 1
**Status:** OPEN

**Description:** Czy α (fine structure constant) wyłania się z TGP, czy musi być input?

Standard QED: α jest input parameter (set by experiment).
TGP framework: czy α derivable z S05 axiom + struktur?

**Decision tree:**
- α native → DERIVED possible (most ambitious)
- α input, ratio α/π structural → STRUCTURAL DERIVED
- Neither → EARLY_HALT

### A4: Quantitative α/π reproduction
**Priority:** CRITICAL
**Phase:** 1
**Status:** OPEN

**Description:** Czy formalized hypothesis daje exactly α/π?

**Resolution path:**
- Sympy: derive ratio (boundary phase amp) / (phase volume)
- Compare with α/π ≈ 0.00232282

## IMPORTANT needs

### A5: Higher-order corrections (α², α³)
**Priority:** IMPORTANT
**Phase:** 2 (if Phase 1 pozytywne)
**Status:** OPEN

**Description:** Schwinger α/π jest leading order. Higher orders:
- α/(2π) - 0.328478... (α/π)² + 1.18124... (α/π)³ - ...
- TGP musiałby reprodukować całe expansion

### A6: Cross-check vs full QED calculation
**Priority:** IMPORTANT
**Phase:** 2
**Status:** OPEN

## Decision criteria — QUICK SCAN

**Phase 1 → Phase 2:**
- A1, A2 give konkretny mathematical structure
- A4 daje numerical match z α/π (jeśli α input) lub better
- Mechanism strukturalnie sound

**Phase 1 → EARLY_HALT:**
- A1/A2 nie give sensible boundary structure
- A4 nie reproduces α/π even with α input
- Mechanism remains too vague to formalize

## Status

WIP, początkowa wersja 2026-05-09 (quick scan mode).
