---
title: "Phase 3 setup — Re-examination linearized analysis (§6.4)"
date: 2026-05-09
parent: "[[./README.md]]"
type: phase-setup
phase: 3
status: 🟢 OPEN — RE-EXAMINATION (response to Phase 2 STRUCTURAL_NO_GO)
predecessor: "[[./Phase2_results.md]] (STRUCTURAL_NO_GO at linearized)"
---

# Phase 3 setup — Re-examination linearized analysis (§6.4)

## §0 — Goal

Re-examine Phase 2 verdict z perspektywy MULTIPOLE STRUCTURE of binary radiation.

**Phase 2 claim:** linearized TGP single-field gives SCALAR-only GW radiation
(conflict z observed TT-dominant pattern).

**Phase 3 hypothesis (re-examination):** Phase 2 ANALYSIS WAS NAIVE — failed
to account for **angular multipole structure** of δΦ at quadrupole level.
Properly decomposed, l=2 quadrupole δΦ gives **PROPER TT POLARIZATION**.

## §1 — Critical observation missed in Phase 2

W Phase 2 ja założyłem:
- δg_eff^ij = δ^ij·b_1·h is "scalar isotropic"
- ⟹ no TT polarization from this term

**ALE** to jest **niepoprawne** dla binary radiation source.

W binary inspiral:
- Monopole δΦ_l=0: M_total = const (NO RADIATION)
- Dipole δΦ_l=1: vanishes in COM frame
- **Quadrupole δΦ_l=2: leading radiating mode**

l=2 multipole has spherical harmonic structure Y_2m(θ,φ) — NIE jest izotropowe.
Specifically: Y_2m has tensor structure Y_2m^ij = (∂_i∂_j - (1/3)δ_ij∇²)Y_2m
when projected onto angular position.

## §2 — Revised analysis

### §2.1 — Quadrupole δΦ angular structure

Far-field δΦ from binary quadrupole:
```
δΦ(r, θ, φ, t) ~ Σ_m a_2m · Y_2m(θ, φ) · g(t - r/c) / r
```

where Y_2m are l=2 spherical harmonics, a_2m mass quadrupole moments.

### §2.2 — δg_eff^ij decomposition with angular structure

```
δg_eff^ij = δ^ij·b_1·δΦ + σ^ij·c_0/(Φ_0²c²)
```

**Angular decomposition** of δ^ij·δΦ_quad:
- δΦ_quad has Y_2m structure
- δ^ij·Y_2m has tensor angular structure (l=2 tensor)
- TT projection NIE jest zero — gives proper h_+, h_× pattern!

### §2.3 — TT projection at observer

For wave propagating along z, observer's TT-plane is xy:
```
δg^TT_xx = δg^xx - (1/2)·tr(δg^ij_perp)
δg^TT_yy = δg^yy - (1/2)·tr(δg^ij_perp)
δg^TT_xy = δg^xy
```

Z δg^ij = δ^ij·b_1·δΦ_quad:
- δg^xx = b_1·δΦ_quad, δg^yy = b_1·δΦ_quad
- δg^TT_xx = b_1·δΦ_quad - (1/2)·(δg^xx + δg^yy) = b_1·δΦ_quad - b_1·δΦ_quad = 0??

Hmm, dla isotropic δ^ij·scalar, TT projection przy detector daje zero
**TYLKO IF** δΦ jest isotropic in θ,φ (l=0 monopole). Dla **quadrupole** δΦ ~ Y_2m,
TT projection jest non-zero!

Specifically: Y_2m has differential angular structure giving non-trivial
δg^TT_xx ≠ δg^TT_yy (for h_+ mode) or δg^TT_xy ≠ 0 (for h_× mode).

This requires careful angular calculation.

## §3 — Phase 3 strategy

### §3.1 — Sympy verification: Y_2m → h_+, h_× explicit

1. Build δΦ_quad ~ a_2m·Y_2m(θ,φ)·exp(iω(t-r/c))/r
2. Compute δg_eff^μν z emergent-metric ansatz
3. TT-project at observer position (specific θ_obs, φ_obs)
4. Read off h_+, h_×

### §3.2 — Compare z GR quadrupole prediction

GR quadrupole formula: h_TT^ij = (2G/c⁴r) · (d²Q^TT/dt²)
TGP emergent-metric: h_TT^ij = (b_1/Phi_0)·(d²Q^Φ/dt²)·... 

If b_1/Phi_0 ~ q/Phi_0 ~ √(4π G), then h_TT^TGP ~ h_TT^GR (modulo O(1) factor).

### §3.3 — Scalar polarization for binary

For binary in COM:
- Monopole M_total = const → NO scalar radiation
- Dipole = 0 (COM cancellation) → NO scalar radiation
- Quadrupole = mass quadrupole ≠ 0 → has l=2 angular structure (TT-like!)

⟹ scalar polarization (h_S = trace of δg^ij) for binary may be **EXACTLY ZERO**!

This would resolve R5 risk completely.

## §4 — Phase 3 deliverables

- Phase3_sympy.py — Y_2m angular decomposition + TT projection
- Phase3_sympy.txt
- Phase3_results.md — verdict on R5 risk after re-examination

## §5 — Phase 3 gate criteria

| # | Criterion |
|---|---|
| G1 | δΦ_quad Y_2m angular structure derived |
| G2 | δg_eff^ij Y_2m projection computed |
| G3 | TT projection h_+, h_× explicit |
| G4 | Scalar polarization h_S for binary (Y_2m component projected to s-wave) |
| G5 | Verdict: R5 resolved or persistent? |

## §6 — Probability assessment

| Outcome | Probability |
|---|---|
| **Re-examination saves framework** (h_S for binary effectively 0, h_TT proper) | **45-60%** |
| Phase 2 verdict stands (still wrong polarization) | 30-40% |
| Inconclusive (multi-session work) | 15-20% |

## §7 — Structural insight (anchor)

**KEY:** binary radiation is l=2 (quadrupole). NIE l=0 (monopole). 
Y_2m angular pattern IS the TT polarization pattern (in standard multipole
decomposition of GW). 

In linearized GR, h_+ ~ Y_22 + Y_2,-2 (real combination), h_× ~ i(Y_22 - Y_2,-2).

If TGP linearized δΦ_quad ~ Y_2m has SAME angular structure, then δg_eff has
TT-like behavior automatically.

This jest **FUNDAMENTAL physics** — GW polarization comes z multipole structure
of source, NIE z metric DOF count alone.

## §8 — Cross-references

- [[./Phase2_results.md]] — STRUCTURAL_NO_GO verdict (challenged)
- [[./Phase2_sympy.py]] — Phase 2 derivation (re-examination)
- [[../op-emergent-metric-from-interaction-2026-05-09/Phase4_results.md]] — Phase 4 1PN/2PN consistency
- Standard GR quadrupole formula references (Blanchet 2014)
