---
title: "ζ.1 program — neutrino mass-spectrum + PMNS first-principles"
date: 2026-04-29
type: program
status: ACTIVE
predecessor: "[[../op-eps-photon-ring/Phase3_results.md]]"
tags:
  - TGP
  - zeta-mass-spectrum
  - neutrino-masses
  - PMNS
  - Koide-K-half
  - cross-sector
---

# ζ.1 program — neutrino mass-spectrum + PMNS first-principles

> **Cel:** Zamknąć neutrino mass-spectrum (Σm_ν = 59.6 meV) jako TGP
> structural prediction via Koide K(ν) = 1/2 (Majorana, B²=1); zamknąć
> PMNS angles (θ₁₂, θ₂₃, θ₁₃) first-principles z GL(3,𝔽₂) × Z₃ × SU(2)_L
> structure; cross-sector mass-spectrum unification (lepton+quark+ν)
> via existing TGP locks: r_21=206.77 (φ-FP), r_31=3477 (Koide K=2/3),
> N_gen=3 (metric singularity barrier), θ_C=0.22550 (GL form factor).
> Status target: PARTIALLY DERIVED → DERIVED via DESI 2027+ + JUNO 2027+
> falsification windows.

---

## Strategiczny kontekst

Post-ε.1 program END (2026-04-29, ledger 391), TGP closes:
- **Photon-ring scale**: ε_ph = 23/137 prime-137 decomposition (PARTIALLY
  DERIVED refined)
- **UV completion**: AS NGFP unique best UV-route, N_A = 500/57 sympy-exact
  (DERIVED refined²)
- **Cross-sector identity**: √α₀ = κ_TGP (DERIVED refined²)
- **ξ-factor frame**: F4 = 1-loop a₂-corrected, Frame B = tree-level
  (DERIVED refined²)

**Open question dla ζ.1:** Mass-spectrum first-principles closure:
1. **Lepton sector**: r_21 = 206.77, r_31 = 3477 LOCKED via φ-FP (tgp-leptons),
   K_lepton = 2/3 sympy exact via thm:koide-theorem
2. **Quark sector**: K_quark = 0.8746/0.7398 RG-invariant (NOT 2/3, separate)
3. **Neutrino sector**: K(ν) = 1/2 (Majorana, B²=1) LOCKED;
   Σm_ν = 59.6 meV from Δm²_21 + Δm²_31 + K=1/2 (NOT first-principles
   z TGP — observational Δm² inputs)
4. **PMNS angles**: STRUCTURAL z GL(3,𝔽₂) × SU(2)_L interplay; sin²θ₁₂≈1/3,
   sin²θ₂₃≈1/2, sin²θ₁₃≈λ_C²/2 — but full unitary closure pending
5. **Cabibbo angle θ_C**: LOCKED w tgp-core (GL form factor 165/167,
   λ_C = 0.22550 vs PDG 0.22500 ± 0.00067, 0.7σ)

**Pre-cycle audit:**
- TESTED-PASS: r_21 (0.001%), r_31 (0.036%), θ_C (0.7σ), N_gen (exact)
- LIVE: Σm_ν = 59.6 meV (DESI DR2 2027+), PMNS θ₁₃ (JUNO 2027+)
- STRUCTURAL: PMNS angles (need explicit GL × SU(2)_L unitary closure)
- OPEN: NEUTRINO MSW substrate derivation (orthogonal track)

→ **Hypothesis:** Mass-spectrum **structurally closed** w obrębie TGP
single-Φ + Z₂ + Koide K_lepton/K_quark/K_neutrino taxonomy + GL(3,𝔽₂) ×
Z₃ × SU(2)_L mixing. ζ.1 audyts robustness pod current observational
inputs i generuje predictions dla DESI 2027+ + JUNO 2027+ falsification.

---

## 3-phase plan

### **Phase 1: Σm_ν = 59.6 meV neutrino mass-spectrum robustness audit (5 sub-tests)**

**Cel:** Verify Σm_ν = 59.6 meV TGP prediction (K(ν)=1/2 + observational
Δm²) robustness pod DESI DR2 bound 0.12 eV; falsify alternative neutrino
mass orderings; confirm normal ordering exclusivity.

**Sub-tests:**
- Z1.1: K(ν) = 1/2 LOCKED (Majorana, B²=1, 2 chiralities collapse)
- Z1.2: Δm²_21 = 7.53e-5 eV², |Δm²_31| = 2.453e-3 eV² (NuFit 5.3 inputs)
- Z1.3: Σm_ν = m1 + m2 + m3 = 59.6 meV (NO ordering, K=1/2 closure)
- Z1.4: DESI DR2 0.12 eV bound vs TGP 0.0596 eV (50% margin)
- Z1.5: Inverted ordering FORBIDDEN by K(ν) = 1/2 incompatibility

### **Phase 2: PMNS angles first-principles z GL(3,𝔽₂) × SU(2)_L (7 sub-tests)**

**Cel:** Derive PMNS mixing angles (θ₁₂, θ₂₃, θ₁₃) z GL(3,𝔽₂) × Z₃ × SU(2)_L
group structure; falsify alternative parameterizations; show cross-sector
lepton-quark θ₁₃ = λ_C²/2 unification.

**Sub-tests:**
- Z2.1: sin²θ₁₂ = 1/3 trimaximal (S₃ permutation symmetry)
- Z2.2: sin²θ₂₃ = 1/2 maximal (Z₂ atmospheric reflection)
- Z2.3: sin²θ₁₃ = λ_C²/2 ≈ 0.026 (cross-sector Cabibbo lock)
- Z2.4: 5 alternative parameterizations falsification (TBM, BM, golden ratio, etc.)
- Z2.5: PMNS unitarity preserved (UU† = I sympy)
- Z2.6: NGFP RG-stability of mixing angles (common β-rescaling)
- Z2.7: Classification PARTIALLY DERIVED (refined)

### **Phase 3: predictions Z1-Z6 + cross-sector unification (6 sub-tests)**

**Cel:** Generate 6 falsifiable predictions about mass-spectrum + PMNS
across TGP sectors; cross-sector lepton-quark-neutrino unification audit;
classification PARTIALLY DERIVED.

**Sub-tests:**
- Z3.1: Z1: DESI DR2 2027+ Σm_ν bound tightening
- Z3.2: Z2: JUNO 2027+ θ₁₃ precision (target <0.1σ z TGP λ_C²/2)
- Z3.3: Z3: T2K/NOνA 2027+ θ₂₃ octant resolution
- Z3.4: Z4: Cross-sector K-taxonomy: K_lepton=2/3, K_quark≠2/3, K_ν=1/2
- Z3.5: Z5: Lepton-quark θ_C-θ₁₃ unification (single Cabibbo anchor)
- Z3.6: Z6: 4-channel falsification convergence (DESI + JUNO + T2K/NOνA + μ→eγ)

---

## Verdict gates

- **Phase 1 ≥ 4/5 PASS** → Phase 2 proceeds
- **Phase 2 ≥ 6/7 PASS** → PMNS first-principles confirmed, Phase 3 proceeds
- **Phase 3 ≥ 5/6 PASS** → ζ.1 program END, classification PARTIALLY DERIVED

**Total:** 18 sub-tests across 3 phases, ledger 391 → 409.

---

## Status taxonomy expectations

- **K(ν) = 1/2** (Majorana, B²=1): **DERIVED** (Koide chirality-counting)
- **Σm_ν = 59.6 meV**: **STRUCTURAL** (K=1/2 + observational Δm²)
- **PMNS angles (θ₁₂=trimaximal, θ₂₃=maximal, θ₁₃=Cabibbo-lock)**:
  **PARTIALLY DERIVED (refined)** post-Phase 2
- **Cross-sector K-taxonomy**: **DERIVED** (chirality-counting per sector)

---

## Materiał wykonawczy

- **Phase 1:** [`Phase1_setup.md`](Phase1_setup.md) + [`phase1_neutrino_audit.py`](phase1_neutrino_audit.py) + [`Phase1_results.md`](Phase1_results.md)
- **Phase 2:** [`Phase2_setup.md`](Phase2_setup.md) + [`phase2_pmns_derivation.py`](phase2_pmns_derivation.py) + [`Phase2_results.md`](Phase2_results.md)
- **Phase 3:** [`Phase3_setup.md`](Phase3_setup.md) + [`phase3_zeta_predictions.py`](phase3_zeta_predictions.py) + [`Phase3_results.md`](Phase3_results.md)

---

## Cross-references

- [`../op-eps-photon-ring/Phase3_results.md`](../op-eps-photon-ring/Phase3_results.md) — ε.1 program END (predecessor)
- [`../../PREDICTIONS_REGISTRY.md`](../../PREDICTIONS_REGISTRY.md) — L1-L5, C2-C3 (existing lepton/Cabibbo locks)
- [`../../INDEX.md`](../../INDEX.md) — master ledger 391 → 409
