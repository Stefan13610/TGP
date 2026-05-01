---
title: "ω.3.Phase3 setup — predictions + 4-channel convergence + program END"
date: 2026-05-01
cycle: ω.3.Phase3
status: SETUP
parent: "[[program.md]]"
predecessor: "[[Phase2_results.md]]"
tags:
  - TGP
  - omega3
  - phase3
  - predictions
  - convergence
  - setup
---

# ω.3.Phase3 setup

**Score gate:** ≥5/6 PASS = ω.3 program END (FULL CONVERGENCE).

## Sub-tests (6)

- **O3.1** **PVLAS-V 2030+ Light-Shining-Through-Walls NULL forecast**:
  Detection probability P_LSW ∝ (g_aγ·B·L)⁴; for TGP g_aγ = 1.71·10⁻²⁰ +
  PVLAS-V projected sensitivity g_aγ ≥ 6.6·10⁻¹¹: TGP signal 1.6·10⁻³⁸ ratio
  → **NULL forecast at all current + 2030+ LSW experiments**.

- **O3.2** **ALPS-II / IAXO 2030+ haloscope/helioscope NULL**:
  - ALPS-II 2030+ projected g_aγ ≥ 2·10⁻¹¹ GeV⁻¹
  - IAXO 2030+ projected g_aγ ≥ 1·10⁻¹² GeV⁻¹ (solar axion)
  - TGP g_aγ = 1.71·10⁻²⁰ → **8 OOM below IAXO** → all helioscopes NULL.
  - **NULL forecast** structurally LOCKED.

- **O3.3** **ADMX/HAYSTAC haloscope band NULL** (ALP regime, m_a free):
  - ADMX/HAYSTAC scan QCD axion band [μeV, meV] = [10⁻⁶, 10⁻³] eV
  - TGP m_a is FREE PARAMETER (ALP, no QCD anomaly N) → if m_a outside
    band, NULL; if m_a within band, depends on coupling
  - For TGP g_aγ = 1.71·10⁻²⁰: 5 OOM below ADMX-Run3 sensitivity even if
    m_a happens to lie in band → **structural NULL forecast**.

- **O3.4** **CMB birefringence cross-check β cosmological**:
  TGP-canonical β = (g_aγ/2)·∫(∂(ln X)/∂η) dη — already locked z ω.1+ω.2;
  ω.3 reproduces consistent prediction post-Planck PR4+ACT 2024 ~3.8σ
  LIVE PARTIAL candidate (downgraded 2026-05-01: hint, not confirmed signal —
  awaits SO/LiteBIRD 2027+ corroboration).

- **O3.5** **ALP fuzzy DM cosmology forecast**:
  - m_a free parameter in TGP ALP regime → no specific TGP DM mass
    prediction
  - if future ω.4 cycle locks m_a (e.g., m_a ~ √(g·f_a·H_inf) or
    similar substrate-action mechanism), DM cosmological band falsifiable
  - **Forward-gate ω.4+** (research-track): m_a structural derivation post-ω.3.

- **O3.6** **4-channel ω.3 convergence summary**:
  - (1) UV-anchor (g* = 0.71): UV.1 NGFP exact
  - (2) Photon-ring (N_A = 500/57): ξ.1 exact
  - (3) Substrate-scale (M_TGP = K·M_GUT): UV.2 LOCKED
  - (4) Triangle anomaly (E_TGP = 536/75): ω.2 LOCKED
  - All 4 channels flow into f_a = (N_A·2π²·M_GUT)/E_TGP = 3125·π²·M_GUT/1273

## Inputs from Phase 2 (LOCKED)

```
f_a sympy-LOCK   = 3125·pi^2·M_GUT/1273
f_a numerical    = 4.8456·10^17 GeV (drift 0.30% vs chi.1)
g_a-gamma        = 1.7129·10^-20 GeV^-1
H_inf max        = 3.04·10^7 GeV (theta_i=1, isocurvature)
r max            = 3.17·10^-23 (tensor-scalar ratio)
F-cluster max drift = 0.083%
```

## Strategy

**Phase 3** delivers the predictive NULL forecasts at all current axion-photon
experiments (PVLAS, ALPS, IAXO, CAST, ADMX/HAYSTAC), reproduces the CMB
birefringence prediction post-ω.3 (consistent with ω.1+ω.2 LIVE PARTIAL
candidate), and confirms the 4-channel ω.3 convergence z UV.1 + ξ.1 + UV.2 +
ω.2 cascade. Forward-gate dla ω.4+ structural m_a derivation noted explicitly.

## Cross-references

- [[program.md]]
- [[Phase1_results.md]]
- [[Phase2_results.md]]
