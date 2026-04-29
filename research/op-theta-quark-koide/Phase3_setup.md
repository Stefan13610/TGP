---
title: "θ.1.Phase3 setup — 6 predictions Q1-Q6 + K-taxonomy 4-sector completion"
date: 2026-04-29
cycle: θ.1.Phase3
status: PRE-EXECUTION
parent: "[[program.md]]"
predecessor: "[[Phase2_results.md]]"
tags:
  - TGP
  - theta-quark-koide
  - CKM
  - predictions
  - cross-sector
  - falsification-roadmap
---

# θ.1.Phase3 — 6 predictions Q1-Q6 + K-taxonomy 4-sector completion

> **Cel:** Generate 6 falsifiable predictions across quark + CKM sektor;
> 4-sector K-taxonomy completion (lepton-ν-up-down); classification
> PARTIALLY DERIVED (refined) confirmed; θ.1 program END.

---

## 6 sub-tests / predictions

### T3.1 (Q1) — Belle II 2027+ V_ub precision target

**Prediction:**
- TGP V_ub = A · λ_C³ · √(ρ̄² + η̄²) = 0.00348 (current)
- PDG 2024: |V_ub| = 0.00382 ± 0.00010 (drift 8.98%)
- Belle II projected (2027+): σ(|V_ub|_excl) ~ 1.5%, σ(|V_ub|_incl) ~ 2%
- **Falsification gate:** if Belle II measures |V_ub| outside [3.40, 4.00] · 10⁻³
  at >5σ, λ_C³ cascade z (ρ̄, η̄) PDG framework broken
- **Confirmation gate:** if Belle II resolves |V_ub| ≈ 3.5-3.8 · 10⁻³,
  TGP λ_C³ cascade locked with refined (ρ̄, η̄)

### T3.2 (Q2) — LHCb Run 4 (2030+) Jarlskog invariant J

**Prediction:**
- TGP J = A² · λ_C⁶ · η̄ = 0.0009² · 0.357 · ... ≈ 2.95 · 10⁻⁵
- PDG 2024: J = (3.07 ± 0.10) · 10⁻⁵
- TGP drift vs PDG: ≈ 4%
- LHCb Run 4 projected: σ(J) ~ 1% absolute
- **Falsification gate:** if LHCb J outside [2.85, 3.30] · 10⁻⁵ at >5σ,
  Wolfenstein λ⁶ cascade broken
- **Confirmation gate:** if J ≈ 3.0 · 10⁻⁵ confirmed, CP-violation
  CKM area cross-sector lock

### T3.3 (Q3) — EIC 2030+ proton mass-radius cross-check

**Prediction:**
- TGP no direct proton mass-radius prediction — but K_up = 7/8 sympy lock
  implies QCD running γ_m universal across u/c/t to 10⁻⁵
- EIC 2030+ projected: proton mass-radius < r_M > with ~1% precision
- Cross-check: if EIC measures m_q running consistent z γ_m universal,
  K_up = 7/8 RG-stability confirmed empirically
- **Falsification gate:** if EIC reveals flavor-dependent γ_m at >5%,
  K_up RG-invariance broken (would require K_up(μ) running)
- **Confirmation gate:** if universal γ_m confirmed within 1%, K_up = 7/8
  sympy lock promoted PARTIALLY DERIVED → DERIVED

### T3.4 (Q4) — Lepton-quark cross-sector λ_C anchor (PMNS-CKM unification)

**Prediction:**
- TGP single Cabibbo λ_C = 0.22550 governs:
  - V_us (CKM) = λ_C exact
  - sin θ₁₃ (PMNS) = λ_C / √2 (ζ.1 lock)
  - Cross-sector ratio sin θ_C / sin θ₁₃ = √2 EXACT (TGP)
- Observed ratio (2026): 0.22500 / 0.14832 = 1.517 (drift 7.3%)
- **Falsification gate:** if JUNO 2027+ refines sin θ₁₃ such that
  ratio drifts > 20% from √2, single λ_C anchor framework broken
- **Confirmation gate:** if ratio converges to √2 ± 1% via JUNO precision,
  PMNS-CKM lepton-quark unification promoted DERIVED

### T3.5 (Q5) — K-taxonomy 4-sector completion

**Prediction:**
- TGP universal pattern: K = (2 + B²)/(2N), N=3
  - K_lepton = 2/3 (Dirac, B²=2) — DERIVED
  - K_neutrino = 1/2 (Majorana, B²=1) — DERIVED (ζ.1)
  - K_up = 7/8 (Dirac+QCD, B²_up = 13/4) — PARTIALLY DERIVED (θ.1)
  - K_down ≈ 0.7399 (Dirac+QCD, B²_down ≈ 2.44) — STRUCTURAL refined
- **4-sector taxonomy completed** with 3 LOCKED + 1 STRUCTURAL
- **Falsification gate:** if any sektor disconfirms K = (2+B²)/(2N) form
  by >5%, universal taxonomy broken
- **Confirmation gate:** if K_down B² = 11/25 + 2 = 61/25 sympy verified
  via further QCD-running corrections, all 4 sectors LOCKED → DERIVED

### T3.6 (Q6) — 4-channel falsification convergence (Belle II + LHCb + EIC + cross-PMNS)

**Prediction:**
- Belle II 2027+: |V_ub| precision (Q1)
- LHCb Run 4 2030+: Jarlskog J (Q2)
- EIC 2030+: proton mass-radius (Q3 cross-check QCD γ_m)
- JUNO 2027+: PMNS sin θ₁₃ ↔ CKM V_us cross-sector lock (Q4)
- **Convergence:** ≥3 z 4 channels must converge within 5σ of TGP
  predictions for classification stabilization
- **Falsification gate:** if ≥2 z 4 channels reject TGP > 5σ, θ.1
  classification PARTIALLY DERIVED → reverts to STRUCTURAL

---

## Verdict gate

**6/6 PASS** → θ.1 program END, classification PARTIALLY DERIVED (refined),
status cascade ACTIVATED dla K_up = 7/8 sympy LOCKED.

**5/6 PASS** → θ.1 program END z minor gap, partial cascade.

**≤ 4/6 PASS** → θ.1.Phase3 reframing required.

---

## Środowisko

```bash
PYTHONIOENCODING=utf-8 python -X utf8 research/op-theta-quark-koide/phase3_theta_predictions.py 2>&1 | tee research/op-theta-quark-koide/phase3_theta_predictions.txt
```

## Cross-references

- [`program.md`](program.md) — overall θ.1 plan
- [`Phase1_results.md`](Phase1_results.md) — K_quark numerical landscape LOCKED
- [`Phase2_results.md`](Phase2_results.md) — K_up = 7/8 sympy LOCKED
- [`../../PREDICTIONS_REGISTRY.md`](../../PREDICTIONS_REGISTRY.md) — Q1-Q6 entries (LIVE 2027+)
