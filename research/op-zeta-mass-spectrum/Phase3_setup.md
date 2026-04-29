---
title: "ζ.1.Phase3 setup — 6 predictions Z1-Z6 + cross-sector K-taxonomy unification"
date: 2026-04-29
cycle: ζ.1.Phase3
status: PRE-EXECUTION
parent: "[[program.md]]"
predecessor: "[[Phase2_results.md]]"
tags:
  - TGP
  - zeta-mass-spectrum
  - PMNS
  - predictions
  - cross-sector
  - falsification-roadmap
---

# ζ.1.Phase3 — 6 predictions Z1-Z6 + cross-sector K-taxonomy unification

> **Cel:** Generate 6 falsifiable predictions across mass-spectrum + PMNS
> sektorów; cross-sector K-taxonomy lepton-quark-neutrino audit;
> classification PARTIALLY DERIVED (refined) confirmed; ζ.1 program END.

---

## 6 sub-tests / predictions

### Z3.1 (Z1) — DESI DR2 → DR3 Σm_ν bound tightening

**Prediction:**
- TGP Σm_ν = 59.01 meV LOCKED
- DESI DR2 (current 95% CL) bound: < 72 meV → margin +22%
- DESI DR3 (2027+ projected) bound: < 40 meV → margin -32%
- **Falsification gate:** if DR3 sets Σm_ν < 0.040 eV at 95% CL, TGP
  K(ν)=1/2 framework sfalsyfikowana; ζ.1 reopens
- **Confirmation gate:** if DR3 detects Σm_ν ≈ 0.059 eV via Λ→ΛCDM
  consistency, TGP K(ν)=1/2 promoted PARTIALLY DERIVED → DERIVED

### Z3.2 (Z2) — JUNO 2027+ θ₁₃ precision target

**Prediction:**
- TGP sin²θ₁₃ = λ_C²/2 = 0.0254251 (drift 15.6% vs current NuFit)
- TGP sin θ₁₃ = λ_C/√2 = 0.15945 (drift 7.5%)
- JUNO 2027+ projected precision: σ(sin²2θ₁₃) ~ 0.5% absolute
- **Falsification gate:** if JUNO measures sin²2θ₁₃ outside TGP ± 1σ
  (~0.085 ± 0.003), λ_C²/2 framework broken
- **Confirmation gate:** if JUNO confirms sin²2θ₁₃ ≈ 0.083-0.087,
  cross-sector Cabibbo lock promoted

### Z3.3 (Z3) — T2K/NOνA + DUNE/T2HK θ₂₃ octant resolution

**Prediction:**
- TGP sin²θ₂₃ = 1/2 (zeroth-order, maximal); current NuFit prefers
  2nd octant 0.572 (drift 12.6%)
- T2K/NOνA tension on octant unresolved as of 2026
- DUNE 2030+ + T2HK 2030+ projected: > 5σ octant determination
- **Falsification gate:** if DUNE/T2HK confirm 2nd octant > 5σ + drift
  > 20%, Z₂ atmospheric framework needs 1-loop refinement
- **Confirmation gate:** if maximal mixing confirmed within 5%, K(ν)=1/2
  + Z₂ atmospheric promoted

### Z3.4 (Z4) — Cross-sector K-taxonomy: K_lepton=2/3, K_quark≠2/3, K_ν=1/2

**Prediction:**
- K_lepton = 2/3 (Dirac, B²=2; matches PDG to 10⁻⁵)
- K_neutrino = 1/2 (Majorana, B²=1; sympy exact)
- K_quark ≈ 0.81-0.87 (RG-invariant range; NOT 2/3)
- **Cross-sector taxonomy:** chirality-counting (B²) governs Koide
  per sector → universal pattern K = (2 + B²)/(2N)
- **Falsification gate:** if quark sektor confirms K_quark = 2/3,
  Dirac/Majorana taxonomy collapses
- **Test:** existing PDG quark mass data + RG running confirm K_quark
  range 0.81-0.87 (not 2/3 nor 1/2)

### Z3.5 (Z5) — Lepton-quark θ_C-θ₁₃ unification (single Cabibbo anchor)

**Prediction:**
- Quark CKM:    V_us = sin θ_C = λ_C = 0.22550 (TGP, 0.22% drift vs PDG 0.22500)
- Lepton PMNS:  sin θ₁₃ = λ_C/√2 = 0.15945 (TGP, 7.5% drift vs NuFit)
- **Single GL(3,𝔽₂) form factor 165/167** governs both sectors
- Cross-sector ratio: sin θ_C / sin θ₁₃ = √2 EXACT (TGP)
- Observed ratio: 0.22500 / 0.14832 = 1.5170 (vs √2 = 1.4142, drift 7.3%)
- **Falsification gate:** if cross-sector ratio drifts > 20%, lepton-quark
  Cabibbo unification framework broken

### Z3.6 (Z6) — 4-channel falsification convergence

**Prediction:**
- DESI DR3 2027+: Σm_ν < 40 meV bound (Z1)
- JUNO 2027+: σ(sin²2θ₁₃) ~ 0.5% (Z2)
- DUNE/T2HK 2030+: θ₂₃ octant > 5σ (Z3)
- μ→eγ MEG-II 2027+: BR < 6×10⁻¹⁴ (Z6 cross-sector check)
- **Convergence:** ≥3 z 4 channels must converge within 5σ of TGP
  predictions dla classification stabilization
- **Falsification gate:** if ≥2 z 4 channels reject TGP > 5σ, ζ.1
  classification PARTIALLY DERIVED → reverts to STRUCTURAL

---

## Verdict gate

**6/6 PASS** → ζ.1 program END, classification PARTIALLY DERIVED (refined),
status cascade ACTIVATED jeśli applicable.

**5/6 PASS** → ζ.1 program END z minor gap, partial cascade.

**≤ 4/6 PASS** → ζ.1.Phase3 reframing required.

---

## Środowisko

```bash
PYTHONIOENCODING=utf-8 python -X utf8 research/op-zeta-mass-spectrum/phase3_zeta_predictions.py 2>&1 | tee research/op-zeta-mass-spectrum/phase3_zeta_predictions.txt
```

## Cross-references

- [`program.md`](program.md) — overall ζ.1 plan
- [`Phase1_results.md`](Phase1_results.md) — Σm_ν = 59 meV closure
- [`Phase2_results.md`](Phase2_results.md) — PMNS first-principles closure
- [`../../PREDICTIONS_REGISTRY.md`](../../PREDICTIONS_REGISTRY.md) — Z1-Z6 entries (LIVE 2027+)
