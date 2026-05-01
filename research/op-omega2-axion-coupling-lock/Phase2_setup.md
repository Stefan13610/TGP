---
title: "ω.2.Phase2 setup — sympy LOCK + cross-channel matching"
date: 2026-05-01
cycle: ω.2.Phase2
status: SETUP
parent: "[[program.md]]"
predecessor: "[[Phase1_results.md]]"
tags:
  - TGP
  - omega2
  - phase2
  - sympy-lock
  - cross-channel
  - setup
---

# ω.2.Phase2 setup

**Score gate:** ≥6/7 PASS = Phase 2 forward.

## Sub-tests (7)

- **W2.2.1** CMB β constraint: g·Δ(ln X)_cosmo = 2β = 0.01186 rad — solve dla wszystkich 5 candidates
- **W2.2.2** Δ(ln X) cosmological z φ.1 EL eq w FRW background
- **W2.2.3** PVLAS-V 2030+ projected sensitivity check (g/f_a < 1·10⁻¹¹ GeV⁻¹)
- **W2.2.4** Magnetar E·B sensitivity (FAST/SKA 2030+ pulsar-timing anomaly)
- **W2.2.5** Quasar Δχ/ln(1+z) = g/2 coefficient lock — SKA 2030+ z>4
- **W2.2.6** Combined-channel χ² ranking → unique winner
- **W2.2.7** UV-IR matching consistency (anomaly one-loop universal)

## Strategy

Phase 1 zidentyfikował **g_anomaly = α_em·E_TGP/(2π) = α_em·(536/75)/(2π)**
jako jedyną strukturalnie wyprowadzoną postać z TGP B²-chirality counting.
Phase 2 weryfikuje czy ta forma rzeczywiście zwycięża cross-channel
falsification vs 4 historic LOCK candidates (κ_TGP, α_em alone, 1/(2π), η_chir).

**Falsifier nadrzędny:** jeśli kombinowane dane (CMB + magnetar + quasar)
preferują INNE g niż g_anomaly z ≥3σ, to TGP B²-anomaly ansatz jest
falsyfikowany strukturalnie.

## Cross-references

- [[program.md]]
- [[Phase1_results.md]]
