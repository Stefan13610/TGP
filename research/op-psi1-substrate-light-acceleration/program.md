---
title: "ψ.1 program — substrate-engineered light acceleration via L₅ + ω.1 chain"
date: 2026-05-01
cycle: ψ.1
parent: "[[INDEX.md]]"
status: ACTIVE
projection:
  ledger_master: "751 → 769"
  ledger_phase: "503 → 521"
  pattern: "5+7+6 = 18 sub-tests"
tags:
  - TGP
  - psi1
  - substrate-light-acceleration
  - L5
  - omega1
  - sigma1-subleading
  - program
---

# ψ.1 — Substrate-engineered light acceleration via L₅

## Thesis

Sub-leading σ.1 channel: scalar (helicity-INDEPENDENT) modification of effective light speed $c_{local}(\mathbf{x})$ in regions of substrate gradient $\partial(\ln X) \neq 0$, via dimension-6 EFT operator $L_5$ scale-invariant pod φ.1 X→λX. Mechanism is **sourceable in lab** through the **same** ω.1 EOM channel as τ.3 (E∥B parallel field → $\Box(\ln X) = -(g/f_X^2)\,\mathbf{E}\cdot\mathbf{B} \neq 0$ → $\partial(\ln X) \neq 0$), and produces measurable scalar shift $\Delta c/c$ that is **NOT helicity-dependent** (contrast to σ.1 leading birefringence).

User's question: **"czy istnieje coś, co może przyśpieszyć światło w sensie skalarnym, lokalnie zmieniając taktowanie substratu?"** ψ.1 answers: **TAK** strukturalnie, jeśli β_g < 0 (UV matching).

## L₅ candidate Lagrangian

Canonical form (φ.1 scale-invariant, dim-6 EFT-suppressed):

$$L_5 \;=\; -\frac{1}{4}\,\frac{\beta_g}{\Lambda^2}\,(\partial_\mu \ln X)(\partial^\mu \ln X) \cdot F_{\nu\rho}F^{\nu\rho}$$

Total photon kinetic term:

$$L_{em} + L_5 \;=\; -\frac{1}{4}\Big[1 + \frac{\beta_g}{\Lambda^2}(\partial \ln X)^2\Big] F_{\nu\rho}F^{\nu\rho}$$

→ Effective dielectric tensor + permeability shift → **scalar c modification**:

$$\boxed{\frac{\Delta c}{c_0} \;=\; -\frac{1}{2}\,\frac{\beta_g}{\Lambda^2}\,(\partial \ln X)^2}$$

For β_g < 0: $\Delta c/c > 0$ → **light ACCELERATES locally**. Sign determined by UV matching (analogous to τ.3 α_g sign analysis).

## Why this is NOT in σ.1 leading

σ.1 leading dispersion $\omega^2_\pm = k^2 \pm g(n\cdot\hat{k})$ is **helicity-dependent** (axion-like F·F̃ coupling) and **explicitly excludes** scalar c(X) at leading order (Webb/Murphy 1e-7 NULL protected). ψ.1 enters at:
- **higher order** in $(\partial \ln X)$ (quadratic, not linear like σ.1)
- **scalar Lorentz-singlet** (F², not F·F̃)
- **EFT-suppressed** by Λ² (not f_X² as in ω.1)

The same hierarchy as τ.2 ↔ τ.3:
- τ.2: leading clock-rate scale-protection (X-invariant canonical UNIQUE)
- τ.3: sub-leading L₄ (∂lnX)² · m_e (sourceable lab)
- σ.1: leading scalar c(X) PROTECTED + helicity-dependent birefringence
- **ψ.1: sub-leading L₅ (∂lnX)² · F² (scalar c shift, sourceable lab)**

NO TENSION z σ.1 — complementary regime, **structurally analogous to τ.3 promotion of τ.2**.

## Engineering chain (parallel to τ.3)

```
lab E∥B (Schwinger-class lub ELI-NP 2030+)
    ↓ ω.1 EOM
□(ln X) = -(g/f_X²) E·B parallel max
    ↓ Yukawa Greens
∂_μ(ln X) ≠ 0 in field region
    ↓ L₅
Δc/c = -(β_g/2Λ²)(∂lnX)²
    ↓ accumulation along path L
phase shift Δφ = ω·L·(Δc/c)/c₀
time-of-flight delay Δt = L·(Δc/c)/c₀
```

## Lab measurability (key advantage over σ.1 birefringence)

For Schwinger-class lab E∥B, $(\partial\ln X)^2 \sim 10^{-30}$ (from τ.3 Phase 2 T2.5 inheritance), Λ = 100 MeV:

$$|\Delta c/c| \;\sim\; \frac{|\beta_g|}{2 \cdot (100\,\text{MeV})^2} \cdot 10^{-30} \;\sim\; |\beta_g| \cdot 10^{-12}$$

For propagation L = 10 cm, photon $\lambda$ = 1064 nm ($\omega = 1.77 \times 10^{15}$ rad/s):

| Observable | Predicted signal | Sensitivity (today) | SNR (β_g = 1) |
|------------|------------------|---------------------|---------------|
| **Sagnac phase Δφ** | $\sim 6 \times 10^{-7}$ rad | $10^{-11}$ rad (LIGO-class) | **~10⁴** ★ |
| **TOF Δt** | $\sim 3 \times 10^{-22}$ s = 0.3 zs | $10^{-18}$ s (attosec) | $10^{-3}$ → 2030+ |
| **Spectral shift** | $\Delta\omega/\omega \sim 10^{-12}$ | $10^{-18}$/yr (Sr/Yb) | **~10⁶** if integration ✓ |

**Key result:** Sagnac fazowy interferometr is **WYKONALNY DZIŚ** dla β_g ~ O(1), Λ ≲ 100 MeV — w przeciwieństwie do σ.1 birefringence która wymaga 2030+ frontier.

## Λ-cutoff dependency table

| Λ | Λ/m_e | Δc/c numerical | Sagnac Δφ (L=10cm, λ=1µm) | Status |
|---|-------|----------------|---------------------------|--------|
| M_Pl | 2.4×10²² | ~10⁻⁵⁰ | 10⁻⁴² | undetectable |
| TeV | 1.96×10⁶ | ~10⁻²⁸ | 10⁻²⁰ | undetectable |
| GeV | 1.96×10³ | ~10⁻²⁰ | 10⁻¹² | borderline 2035+ |
| **100 MeV** | **196** | **~10⁻¹²** | **~10⁻⁴ rad** | **DETECTABLE Sagnac dziś** |
| 10 MeV | 19.6 | ~10⁻⁸ | ~1 rad | overshoot — already constrained |
| 1 MeV | 1.96 | ~10⁻⁴ | huge | excluded by atomic spec |

**Sagnac LIGO-class 10⁻¹¹ rad threshold ⟹ Λ ≲ ~10 GeV detectable now (much wider window than τ.3 100 MeV).**

## 18 sub-test plan

### Phase 1 — L₅ coupling structural derivation (5 tests)

- **T1.1** L₅ candidate scan + φ.1 X→λX scale-invariance check
- **T1.2** UV matching β_g sign (3 channels: AS NGFP + heavy-mode loop + cosmological consistency)
- **T1.3** Effective scalar c shift derivation $\Delta c/c = -(\beta_g/2\Lambda^2)(\partial\ln X)^2$
- **T1.4** ω.1 EOM source maximization F·F̃ via E∥B parallel + null controls
- **T1.5** Viability gate Λ ≲ 10 GeV detectable via Sagnac LIGO-class today

### Phase 2 — sympy LOCK + Sagnac/TOF engineering (7 tests)

- **T2.1** sympy LOCK $\Delta c/c$ formula
- **T2.2** F² = 2(B²-E²) vs F·F̃ = -4 E·B distinction (E∥B max for F·F̃ source, F² independent)
- **T2.3** Yukawa Greens for substrate gradient (heavy/light substrate regimes)
- **T2.4** Λ-cutoff scan
- **T2.5** Sagnac fazowy SNR + TOF zs SNR numerical at Schwinger-class
- **T2.6** Cross-coupling z σ.1 (sub-leading scalar vs leading helicity-dependent consistency matrix)
- **T2.7** 4 alt-L₅ couplings cross-falsification

### Phase 3 — Lab predictions + 4-channel convergence (6 tests)

- **T3.1** Sagnac fazowy lab E∥B prediction (wykonalne dziś, SNR ~10⁴)
- **T3.2** TOF dual-arm zs-precision frontier 2030+
- **T3.3** Cosmological scalar c shift residual (analogous to TT9, B too small)
- **T3.4** Magnetar atmosphere FRB time-of-flight (NOVEL — FRB delays correlated z line-of-sight B integral)
- **T3.5** 4 alt-L₅ cross-channel falsification pattern
- **T3.6** 4-channel ψ.1 convergence

## Prediction projection (TT13-TT18)

| ID | Prediction | Channel | Status target |
|----|-----------|---------|---------------|
| **TT13** | ★ NOVEL Sagnac fazowy E∥B substrate light acceleration (LAB-ENGINEERED, SNR ~10⁴ today) | lab | LIVE NOVEL |
| **TT14** | ★ TOF dual-arm zs-precision (frontier 2030+ attoclock chain) | lab frontier | LIVE FRONTIER |
| **TT15** | Cosmological scalar c shift residual NULL (B too small, consistent z Webb/Murphy) | cosmo | STRUCTURAL |
| **TT16** | ★ NOVEL FRB time-of-flight delay correlated z B²·DM line-of-sight | astrophysical | LIVE NOVEL |
| **TT17** | 4 alt-L₅ couplings cross-channel falsification (F², F·F̃, (E²-B²)², (E·B)²) | discriminator | STRUCTURAL |
| **TT18** | 4-channel ψ.1 convergence (lab Sagnac + TOF + cosmo + FRB) | meta | LIVE NOVEL |

## Cross-cycle position

φ.1 → ω.1 → σ.1 → τ.2 → τ.3 → **ψ.1** (6th cycle in chain):
- τ.3 promoted τ.2 leading protection by adding sub-leading L₄ (mass channel)
- ψ.1 promotes σ.1 leading birefringence by adding sub-leading L₅ (light-speed channel)
- Together: τ.3 ⊕ ψ.1 = full lab-engineering substrate response of matter (mass) AND light (c)

This is the **second lab-engineering predictive TGP cycle** — τ.3 was first (atomic clocks), ψ.1 is second (interferometric phase / TOF).

## Phase verdict (anticipated)

If 18/18 PASS: ψ.1 program END, ledger 751 → 769. **First TGP cycle providing concrete LIGO-class Sagnac experiment WYKONALNY DZIŚ** that falsifies or measures (β_g, Λ) — eksperyment dostępny w 2026, nie 2030+.

