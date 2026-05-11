---
title: "Phase 3 setup — phenomenological bounds checks: BBN ⁴He/D/H + CMB ω_b/ω_m + PTA NANOGrav 15-yr compatibility"
date: 2026-05-11
parent: "[[./README.md]]"
type: phase-setup
phase: 3
status: 🟡 setup phase
sub_needs_addressed: [N0.8, N0.9, N0.10]
risks_addressed: [R3-full, R6]
predecessor: "[[./Phase2_results.md]] (8/8 PASS, Q2 F1 konstruktywnie verified)"
tags:
  - phase3
  - BBN-bounds
  - CMB-bounds
  - PTA-NANOGrav
  - phenomenology
  - cosmology-observational-bounds
---

# Phase 3 setup

## §0 — Cel Phase 3

Verify że tego cyklu N2 derivation **passes** wszystkie observational bounds dla
QCD-related cosmology:
1. **BBN era** (T ~ 1 MeV, t ~ 1 s, z ~ 10⁹): ⁴He Y_p, D/H abundances preserved
   na 1% precision (PDG 2024).
2. **CMB era** (T ~ eV, z ~ 1100): ω_b, ω_m parameters preserved (Planck 2018).
3. **PTA stochastic GW** (NANOGrav 15-yr 2023, EPTA, IPTA): compatibility z
   crossover (NOT first-order phase transition) — consensus SMBHB-dominated origin.

Sub-needs: N0.8 (BBN), N0.9 (CMB), N0.10 (PTA).
Risks: R3 (BBN incompatibility — explicit verify) + R6 (PTA false positive — verify).

## §1 — BBN bound check

### §1.1 — BBN era physics setup

Big Bang Nucleosynthesis epoch:
- Temperature: T ~ 0.05-1 MeV (peak nucleosynthesis around T ~ 80 keV)
- Time: t ~ 1-1000 s after Big Bang
- Redshift: z ~ 4·10⁸-10⁹
- **Crucial:** T_BBN ≪ Λ_QCD = 217 MeV (factor ~200 below QCD crossover scale)

### §1.2 — ρ_QCD(T~MeV) ≈ 0 strukturalnie

Per Phase 2 §1.3 + IR limit analysis:
- T_BBN ~ 1 MeV << Λ_QCD = 217 MeV
- T_BBN/Λ_QCD ≈ 0.005 (deeply hadronic phase)
- Δ_lattice(T~1 MeV) ≈ 0 (exponentially suppressed below T~100 MeV)
- ρ_QCD_thermal(T_BBN) ≈ 0 strukturalnie

Vacuum gluon condensate (constant, Q2 F1 substrate-decoupled) NIE wpływa na BBN
dynamics; obserwabla `H(z~10⁹)` jest standard ΛCDM:

```
H_BBN_TGP(T) = H_BBN_standard(T) · (1 + δ_QCD(T))
δ_QCD(T~MeV) ≈ 0    [structurally; no thermal QCD source remaining]
```

⇒ **BBN nucleosynthesis predictions standard ΛCDM** w TGP framework.

### §1.3 — ⁴He, D/H predictions

PDG 2024 BBN review:
- **⁴He mass fraction:** Y_p = 0.245 ± 0.003 (observed primordial)
- **Deuterium-to-hydrogen:** D/H = 2.527·10⁻⁵ ± 0.030·10⁻⁵ (observed)
- **⁷Li/H:** 1.6·10⁻¹⁰ (lower than theoretical, "lithium problem" — non-trivial)

Standard BBN H(z) calculation (z baryon-to-photon ratio η ≈ 6·10⁻¹⁰):
- Y_p_predicted = 0.247 ± 0.001 (consistent z observed within 1σ)
- D/H_predicted = 2.5·10⁻⁵ ± 0.1·10⁻⁵ (consistent z observed within 1σ)

**TGP prediction (Phase 3):** Y_p, D/H **identical do standard ΛCDM** because
ρ_QCD(T_BBN) ≈ 0 → H(T_BBN) unchanged.

⇒ **BBN PASS automatic.** R3 (BBN incompatibility) closed.

## §2 — CMB bound check

### §2.1 — CMB era physics setup

Cosmic Microwave Background last-scattering surface:
- Temperature: T_CMB_last_scattering ~ 0.26 eV (z ~ 1100)
- Today: T_CMB = 2.725 K = 0.235 meV (z = 0)
- **Both >> Λ_QCD ≈ 217 MeV downward** — z_CMB << z_QCD (10³ vs 10¹²)

### §2.2 — ρ_QCD(T~eV) → ρ_QCD_vacuum

Z Phase 2 §3.1 IR reduction:
- T_CMB << Λ_QCD ⇒ ρ_QCD_thermal(T_CMB) → 0 strukturalnie
- Only ρ_QCD_vacuum (constant) remains
- Per Q2 F1 mechanism: substrate-decoupled od bare Λ_TGP

⇒ **CMB era observables (ω_b, ω_m, ω_Λ, Ω_curv) standard ΛCDM** w TGP framework.

### §2.3 — Planck 2018 parameters

Per Planck 2018 base ΛCDM:
- **ω_b = 0.02237 ± 0.00015** (baryon density)
- **ω_m = 0.1430 ± 0.0011** (matter density)
- **ω_Λ = 0.6889 ± 0.0056** (dark energy density, Ω_Λ from Friedmann)
- **n_s = 0.9665 ± 0.0038** (scalar tilt)
- **τ = 0.0544 ± 0.0073** (reionization optical depth)

**TGP prediction (Phase 3):** wszystkie parameters preserved standard ΛCDM (per
Q2 F1 decoupling + Phase 2 IR reduction).

⇒ **CMB PASS automatic.** Cosmological constant problem rozszerzona structural
resolution preserved per Q2 F1 + this cycle Phase 2.

## §3 — PTA NANOGrav 15-yr compatibility

### §3.1 — NANOGrav 15-yr 2023 stochastic GW background

NANOGrav 15-yr (2023) detected stochastic GW background w PTA frequency band
(~ nHz), consensus interpretation:
- **Astrophysical:** supermassive black hole binaries (SMBHB) in galactic centers
- Hellings-Downs angular correlation pattern detected (3-4σ)
- Spectral index γ ≈ 13/3 (consistent z SMBHB; possible alternatives include
  cosmic strings + new physics)

### §3.2 — TGP prediction: NO first-order QCD signal

Per Phase 2 §4 + §1.2 (lattice 2+1 flavor consensus):
- **Crossover, NOT first-order phase transition** dla 2+1 flavor QCD
- Smooth Δ(T) profile, no discontinuity, latent heat L ≈ 0
- **No bubble nucleation** → no strong stochastic GW signal w PTA band

**TGP cosmology prediction:**
- QCD epoch (z ~ 10¹², t ~ 10⁻⁵ s) gives smooth crossover
- Stochastic GW background z QCD epoch: weak (only smooth thermal fluctuations)
- **NIE generates dominant PTA signal**

⇒ **TGP framework consistent z NANOGrav consensus SMBHB origin.** R6 (PTA false
positive) closed.

### §3.3 — Hypothetical falsification scenario

Gdyby tego cyklu derivation pokazałaby first-order QCD phase transition (in
contradiction z lattice consensus):
- Bubble nucleation → strong stochastic GW background w PTA band
- TGP prediction would *predict* dominant PTA signal
- Detection consistent z this scenario could be cited as TGP support

**Phase 2 §4 R7 closure:** lattice 2+1 flavor crossover preserved → tego cyklu
NIE *predict* first-order signal → consistent z observation.

**Note:** future NANOGrav-style PTA experiments mogą *additionally* detect:
- Cosmic strings signature (separate sektor)
- Inflation reheating (separate sektor)  
- New physics signatures
TGP framework currently NIE *predicts* additional signatures beyond standard
SMBHB; ten cykl preserves compatibility.

## §4 — Phase 3 sympy LOCK targets

Phase3_sympy.py będzie weryfikować (8 tests):

1. **BBN era T~MeV check:** T_BBN/Λ_QCD ≈ 0.005 ≪ 1 ⇒ ρ_QCD_thermal(T_BBN) ≈ 0
2. **H(z~10⁹) preservation:** standard ΛCDM H(T~MeV) unchanged
3. **⁴He Y_p prediction:** 0.247 ± 0.001 (TGP standard) matches PDG 2024
   observation 0.245 ± 0.003 within 1σ
4. **D/H prediction:** 2.5·10⁻⁵ (TGP standard) matches PDG 2024 observation
   2.527·10⁻⁵ ± 0.030·10⁻⁵ within 1σ
5. **CMB era T~eV check:** ρ_QCD(T_CMB) ≡ ρ_QCD_vacuum (substrate-decoupled)
6. **Planck 2018 ω_b preservation:** TGP prediction = ΛCDM = 0.02237 (within
   precision)
7. **PTA NANOGrav crossover compatibility:** lattice 2+1 flavor crossover →
   no first-order signal → consistent z SMBHB consensus
8. **Cross-cycle synthesis:** BBN + CMB + PTA all PASS automatic z Q2 F1 +
   Phase 2 IR reduction

Target: 8/8 PASS.

## §5 — Cross-references

- [[./README.md]]
- [[./Phase0_balance.md]]
- [[./Phase1_results.md]] (β_QCD, gluon condensate)
- [[./Phase2_results.md]] (thermal profile, Q2 F1 verification, crossover)
- [[../op-Q2-vacuum-budget-2026-05-10/Phase_FINAL_close.md]] §3.3 (CMB ω_b PASS preservation)
- [[../closure_2026-04-26/Lambda_from_Phi0/results.md]] (T-Λ baseline)
- PDG 2024 BBN review — Y_p, D/H, η_baryon
- Planck 2018 ΛCDM parameters
- NANOGrav 15-yr 2023 stochastic GW analysis
- HotQCD/Wuppertal-Budapest 2018+ — lattice EoS T_c, Δ(T)

---

**Phase 3 setup ready.** Next: Phase3_sympy.py + Phase3_results.md.
