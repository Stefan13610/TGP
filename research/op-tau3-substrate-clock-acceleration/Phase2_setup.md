---
title: "τ.3.Phase2 setup — sympy LOCK δω/ω + lab E·B engineering chain (7 tests)"
date: 2026-04-30
cycle: τ.3.Phase2
status: SETUP
parent: "[[program.md]]"
tags:
  - TGP
  - tau3
  - phase2
  - setup
  - sympy
  - engineering
---

# τ.3.Phase2 setup

7 sub-tests deriving sympy-LOCKED clock-rate-shift formula and lab engineering chain.

## Sub-tests

- **T2.1** sympy LOCK δω/ω formula: m_e_eff = m_e^(0) + (α_g/Λ²)(∂ ln X)²,
  E_n ∝ m_e c² α_em², so δω/ω = δm_e/m_e = (α_g/(Λ² m_e^(0)))(∂ ln X)².
- **T2.2** Lab E·B engineering: F_μν F̃^μν = -4 E·B (parallel maximizes), 
  ω.1 EOM □(ln X) = (g/(4f_X²))F·F̃ = -(g/f_X²) E·B.
- **T2.3** Substrate gradient Green's function: (□ + m_X²)(ln X) = source,
  with m_X = f_X·g effective substrate mass; Yukawa-like Green's function.
- **T2.4** Λ-cutoff regime scan: numerical δω/ω for Λ ∈ {M_Pl, TeV, GeV, 
  100 MeV, 10 MeV, 1 MeV} at Schwinger-class lab fields.
- **T2.5** Clock-rate shift Sr/Yb numerical at E ~ 10¹⁵ V/m, B ~ 10⁶ G ∥ E.
- **T2.6** Cross-coupling z τ.2: L4 enters as sub-leading τ.2 correction
  (consistent with τ.2 R(X) protected at LEADING O(∂ ln X), L4 at O((∂ ln X)²)).
- **T2.7** 4 alt-L4-couplings cross-falsification: (β F·F̃, γ ln(F·F̃), 
  δ (E²-B²), η (E·B)²) — only α_g(∂ ln X)² is φ.1 X→λX scale-invariant.

## Pass criteria

- 7/7 PASS = Phase 3 forward
- ≥6/7 PASS = Phase 3 conditional
- <6/7 = Phase 2 retry / τ.3 abort

## Output

`phase2_tau3_engineering.py` runs all 7 sub-tests sequentially.
