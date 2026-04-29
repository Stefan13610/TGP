---
title: "θ.1 program — quark-sektor Koide K_quark first-principles + K-taxonomy completion"
date: 2026-04-29
type: program
status: ACTIVE
predecessor: "[[../op-zeta-mass-spectrum/Phase3_results.md]]"
tags:
  - TGP
  - theta-quark-koide
  - quark-masses
  - CKM
  - Wolfenstein
  - K-taxonomy
  - cross-sector
---

# θ.1 program — quark-sektor Koide K_quark first-principles + K-taxonomy completion

> **Cel:** Zamknąć quark-sektor Koide K_up = 0.8746 i K_down = 0.7398
> jako TGP structural predictions via chirality-counting taxonomy
> (B²_up, B²_down vs lepton B²=2 / neutrino B²=1); zamknąć Wolfenstein
> λ-cascade (V_us, V_cb, V_ub) z single Cabibbo anchor λ_C = 0.22550
> (cross-sector PMNS-CKM lock z ζ.1); 4-sector K-taxonomy completion
> (K_lepton=2/3, K_ν=1/2, K_up=7/8 candidate, K_down ≈ 0.74). Status target:
> PARTIALLY DERIVED (refined) → DERIVED via Belle II 2027+ + LHCb Run 4
> + EIC 2030+ falsification windows.

---

## Strategiczny kontekst

Post-ζ.1 program END (2026-04-29, ledger 409), TGP closes:
- **Neutrino mass-spectrum**: Σm_ν = 59.01 meV LOCKED via K(ν)=1/2
  (PARTIALLY DERIVED refined)
- **PMNS angles**: sin²θ₁₂=1/3, sin²θ₂₃=1/2, sin²θ₁₃=λ_C²/2 first-principles
  closure z GL(3,𝔽₂) × Z₃ × SU(2)_L (PARTIALLY DERIVED refined)
- **Single Cabibbo anchor**: λ_C = 0.22550 governs both V_us (CKM) i sin θ₁₃
  (PMNS) — cross-sector lepton-quark unification (DERIVED refined)
- **K-taxonomy partial**: K_lepton=2/3 (Dirac, B²=2), K_ν=1/2 (Majorana, B²=1)
  LOCKED; K_quark range 0.81-0.87 (NOT 2/3, separate framework) — open

**Open question dla θ.1:** Quark-sektor first-principles closure:
1. **K_up sympy candidate**: K_up = 0.8746 vs 7/8 = 0.875 → drift 0.046%
   (strong candidate; B²_up = 13/4?)
2. **K_down sympy candidate**: K_down = 0.7398 vs 5/2 / (3+1) ?
   (B²_down = 5/2?, lub inny rational anchor)
3. **K-taxonomy completion**: Universal pattern K = (2 + B²)/(2N)
   for N=3 — quark sektor breaks Dirac/Majorana taxonomy via QCD running
4. **Wolfenstein λ-cascade**: V_us = λ_C, V_cb = A·λ_C², V_ub = A·λ_C³·(ρ̄-iη̄)
   — cross-sector single-anchor lock
5. **CKM unitarity triangle**: J = A²·λ_C⁶·η̄ Jarlskog invariant; cross-check
   z LHCb Run 4 + Belle II 2027+ precision

**Pre-cycle audit:**
- TESTED-PASS: K_lepton=2/3 (10⁻⁵), K(ν)=1/2 (Majorana exact),
  λ_C=0.22550 (0.22%)
- LIVE: K_up RG-invariant 0.8746, K_down RG-invariant 0.7398 (PDG MS-bar @ M_Z)
- STRUCTURAL: K_quark sympy closure (need B²_up, B²_down derivation)
- OPEN: NGFP RG-stability of K_quark over μ ∈ [50 MeV, 100 TeV]

→ **Hypothesis:** Quark-sektor K_up/K_down **structurally closed** via
chirality-counting B² extension + QCD running corrections; Wolfenstein
λ-cascade single-anchor lock z λ_C = 0.22550 (ζ.1 cross-sector inheritance).
θ.1 audyts robustness pod current PDG inputs i generuje predictions
dla Belle II + LHCb + EIC falsification.

---

## 3-phase plan

### **Phase 1: K_quark numerical landscape audit (5 sub-tests)**

**Cel:** Verify K_up=0.8746 i K_down=0.7398 TGP predictions z PDG MS-bar
quark masses @ M_Z; audit RG-invariance over 6 orders of magnitude;
4-sector K-taxonomy distinct values; Wolfenstein λ-cascade vs PDG.

**Sub-tests:**
- T1.1: K_up = 0.8746 sympy z PDG MS-bar @ M_Z (m_u, m_c, m_t)
- T1.2: K_down = 0.7398 sympy z PDG MS-bar @ M_Z (m_d, m_s, m_b)
- T1.3: K_quark RG-invariance (μ ∈ [50 MeV, 100 TeV] common β-rescaling)
- T1.4: 4-sector K-taxonomy distinct values (K_lepton=2/3, K_ν=1/2,
  K_up=0.875, K_down≈0.74) — universal sector-separation
- T1.5: Wolfenstein λ-cascade vs PDG (V_us, V_cb, V_ub)

### **Phase 2: K_quark first-principles z chirality-counting + B² taxonomy (7 sub-tests)**

**Cel:** Derive K_up/K_down z chirality-counting B² extension (B²_up=13/4?,
B²_down=5/2?); falsify alternative formulas; show cross-sector
λ_C single-anchor lock; classification PARTIALLY DERIVED (refined).

**Sub-tests:**
- T2.1: K_up = 7/8 sympy candidate (drift 0.046% vs PDG)
- T2.2: K_down sympy candidates ranking (closest rational anchor)
- T2.3: B²_up, B²_down effective chirality-counting derivation
- T2.4: Cross-sector V_us = λ_C single-anchor lock (Wolfenstein λ ↔ Cabibbo)
- T2.5: 5 alternative quark Koide formulas FALSIFIED
- T2.6: NGFP RG-stability of K_quark (common β-rescaling preservation)
- T2.7: Classification PARTIALLY DERIVED (refined)

### **Phase 3: predictions Q1-Q6 + K-taxonomy completion (6 sub-tests)**

**Cel:** Generate 6 falsifiable predictions about quark sektor + CKM
across TGP windows; 4-sector K-taxonomy unification audit;
classification PARTIALLY DERIVED.

**Sub-tests:**
- T3.1: Q1: Belle II 2027+ V_ub precision (target <1% absolute)
- T3.2: Q2: LHCb Run 4 unitarity triangle Jarlskog invariant J
- T3.3: Q3: EIC 2030+ proton mass-radius (cross-check QCD running)
- T3.4: Q4: Lepton-quark cross-sector λ_C anchor (PMNS-CKM unification)
- T3.5: Q5: K-taxonomy 4-sector completion (K_lepton, K_ν, K_up, K_down)
- T3.6: Q6: 4-channel falsification convergence (Belle II + LHCb + EIC + μ→eγ)

---

## Verdict gates

- **Phase 1 ≥ 4/5 PASS** → Phase 2 proceeds
- **Phase 2 ≥ 6/7 PASS** → K_quark first-principles confirmed, Phase 3 proceeds
- **Phase 3 ≥ 5/6 PASS** → θ.1 program END, classification PARTIALLY DERIVED

**Total:** 18 sub-tests across 3 phases, ledger 409 → 427.

---

## Status taxonomy expectations

- **K_up = 0.8746** (RG-invariant): **STRUCTURAL** → **PARTIALLY DERIVED (refined)** post-Phase 2
- **K_down = 0.7398** (RG-invariant): **STRUCTURAL** → **PARTIALLY DERIVED (refined)** post-Phase 2
- **K_up = 7/8 sympy candidate**: **STRUCTURAL** (drift 0.046%, B²_up extension)
- **Wolfenstein λ-cascade single anchor**: **DERIVED** (λ_C = 0.22550 cross-sector lock)
- **4-sector K-taxonomy**: **DERIVED** (chirality-counting per sector taxonomy)

---

## Materiał wykonawczy

- **Phase 1:** [`Phase1_setup.md`](Phase1_setup.md) + [`phase1_kquark_audit.py`](phase1_kquark_audit.py) + [`Phase1_results.md`](Phase1_results.md)
- **Phase 2:** [`Phase2_setup.md`](Phase2_setup.md) + [`phase2_kquark_derivation.py`](phase2_kquark_derivation.py) + [`Phase2_results.md`](Phase2_results.md)
- **Phase 3:** [`Phase3_setup.md`](Phase3_setup.md) + [`phase3_theta_predictions.py`](phase3_theta_predictions.py) + [`Phase3_results.md`](Phase3_results.md)

---

## Cross-references

- [`../op-zeta-mass-spectrum/Phase3_results.md`](../op-zeta-mass-spectrum/Phase3_results.md) — ζ.1 program END (predecessor)
- [`../../PREDICTIONS_REGISTRY.md`](../../PREDICTIONS_REGISTRY.md) — Q1-Q6 entries (LIVE 2027+)
- [`../../INDEX.md`](../../INDEX.md) — master ledger 409 → 427
