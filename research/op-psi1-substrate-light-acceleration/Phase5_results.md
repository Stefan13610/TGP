---
title: "ψ.1.v2.Phase5 results — eikonal + dispersion + corrected Sagnac 5/5 PASS"
date: 2026-05-01
cycle: ψ.1.v2.Phase5
status: PASS
parent: "[[program.md]]"
tags:
  - TGP
  - psi1
  - phase5
  - correction
  - eikonal
  - sagnac
  - results
---

# ψ.1.v2.Phase5 results — 5/5 FULL CASCADE (correction phase)

## Sub-test outcomes

| ID | Test | Result |
|----|------|--------|
| **T5.1** | Eikonal dispersion $g^{\mu\nu}_{eff} k_\mu k_\nu = 0$ sympy LOCK | ✅ PASS |
| **T5.2** | Anisotropic $c_{local}(\theta)$ sympy LOCK | ✅ PASS |
| **T5.3** | Sagnac chopper differential (E∥B vs E⊥B), realistic SNR | ✅ PASS |
| **T5.4** | Yukawa Greens dla tensor source | ✅ PASS |
| **T5.5** | 4 alt-tensor-L₅' falsification matrix | ✅ PASS |

**Score: 5/5 → Phase 6 forward**

## Key results

### T5.1: Eikonal dispersion relation LOCKED

Z $g^{\mu\nu}_{eff} = \eta^{\mu\nu} + (\xi/\Lambda^2) n^\mu n^\nu$ i konwencji
sygnaturalnej (+,-,-,-):

$$\boxed{g^{\mu\nu}_{eff}\, k_\mu k_\nu = (\omega^2 - |\vec{k}|^2) + \frac{\xi}{\Lambda^2}(n\cdot k)^2 = 0}$$

gdzie $n\cdot k = n^\mu k_\mu = n^0\omega - \vec{n}\cdot\vec{k}$.

sympy diff vs. expected = 0 (po fix Minkowski signature contraction —
n·k = n^μ k_μ z metryką, nie Euclidean dot).

### T5.2: Anisotropic $c_{local}(\theta)$ LOCKED

Statyczny gradient $n^\mu = (0, \vec{n})$, $\vec{n} = \vec\nabla\ln X$, $|\vec{n}| = n_{mag}$,
$\theta = \angle(\vec{k}, \vec{n})$:

$$\boxed{c_{eff}^2(\theta) = 1 - \frac{\xi\, n_{mag}^2}{\Lambda^2}\cos^2\theta}$$

sympy verified:
- $c_{eff}^2(\theta=0) = 1 - \xi n^2/\Lambda^2$ → **parallel: maksymalne spowolnienie**
- $c_{eff}^2(\theta=\pi/2) = 1$ → **perpendicular: BRAK efektu**

Z konwencji Phase 4 ($\xi = +2|\beta_g| > 0$) wszędzie $c_{eff} \leq c_0$ (subluminal).

### T5.3: Sagnac chopper differential — REALISTIC SNR (≪ v1)

**v2 chopper protocol:**
- Configuration A: pętla interferometru w płaszczyźnie zawierającej $\vec\nabla\ln X$
- Configuration B: pętla prostopadła do $\vec\nabla\ln X$
- Differential: $\Delta\phi_{A-B}$ (eliminuje common-mode noise)

**Symbolic result:**
$$\Delta\phi_{A-B} = \frac{4\pi L_{arm}}{\lambda_\gamma} \cdot |\beta_g| \cdot \frac{|\vec\nabla\ln X|^2}{\Lambda^2}$$

**Numerical estimates:**

| Scenario | $|\vec\nabla\ln X|$ | $L_{arm}$ | $\lambda_\gamma$ | $\Delta\phi$ | SNR (1 month, shot-floor 6.2×10⁻¹³ rad) |
|---|---|---|---|---|---|
| Lab-scale realistic | $10^{-3}$ m⁻¹ | 1 m | 1 μm | $5\times 10^{-36}$ rad | $8\times 10^{-24}$ ❌ |
| Magnetar-extreme | 1 m⁻¹ | 1 m | 1 μm | $5\times 10^{-30}$ rad | $8\times 10^{-18}$ ❌ |

**Crucial result:** realistic SNR ≪ v1's claimed 3×10⁴.

ψ.1 lab-Sagnac signal jest **sub-detectable** w obecnych setupach.
**v1 SNR ~3×10⁴ był artefaktem fałszywej skalarnej Δc/c.** Tensor anizotropia
daje znacznie mniejszy efekt, bo:
1. Tylko składowa $\vec{k}\parallel\vec\nabla\ln X$ czuje modyfikację (cos²θ)
2. Konkretny lab gradient $|\vec\nabla\ln X|$ jest tiny (no natural amplification)
3. Wymaga external source (cosmological, astrophysical) — nie tabletop

**Realistic detection channel:** TOF differential dla **astrophysical** sources
(np. magnetar FRB) gdzie cumulative phase shift przez gigaparsek może być
mierzalny w sub-leading O((∂lnX)²) (poprawka do leading σ.1 helicity-dispersion).

### T5.4: Yukawa Greens function for tensor source

Linearized EOM dla efektywnego potencjału generowanego przez tensor source:
$$\Phi_{eff}(r) \sim \frac{|\beta_g|}{\Lambda^2}\cdot \frac{e^{-\Lambda r}}{4\pi r}\cdot (\text{tensor projector})$$

sympy LOCK: $(\nabla^2 - \Lambda^2) G_{Yukawa} = 0$ for $r>0$ ✓.

**Implication:** żaden long-range tensor force — exponential cutoff przy $\Lambda^{-1}$ (sub-fm dla $\Lambda$=100 TeV → no astrophysical/cosmological tensor contamination z tego operator).

### T5.5: 4 alt-tensor falsification matrix

| Operator | $c_{eff}$ aniso | Sagnac chopper | TOF directional | vacuum Cherenkov |
|----------|:---:|:---:|:---:|:---:|
| **L₅'_a** $(\partial\ln X)(\partial\ln X)F\cdot F$ | YES (cos²θ) | **YES (A-B differential)** | YES | safe (subluminal) |
| L₅'_b $(\partial\ln X)(\partial\ln X)F\cdot\tilde F$ | birefringence (helicity-dep) | NO (helicity-locked) | YES (helicity-split) | safe |
| L₅'_c $(\partial\partial\ln X)F\cdot F$ | YES (reduces to L₅'_a) | YES (subset) | YES | safe |
| L₅'_d $\Box(\ln X)F^2$ | NO (scalar — **v1 mistake**) | NO | NO | n/a (no Δc) |

**L₅'_a uniquely** produces Sagnac A-B differential signal:
- L₅'_b helicity-locked, NOT directional
- L₅'_c reduces algebraicznie do L₅'_a (przez parts integration)
- L₅'_d ZERO — to BYŁO ψ.1.v1 fałszywe założenie

**Cross-cycle separation:**
- σ.1 cykl (helicity-dispersion) ⟵ L₅'_b territory (parity-odd)
- ψ.1.v2 cykl (directional anisotropy) ⟵ L₅'_a territory (parity-even, tensor)
- v1 cykl ⟵ L₅'_d territory (gives ZERO Δc, formally NEGATIVE)

## Phase verdict

**ψ.1.v2.Phase 5 PASS (FULL CASCADE 5/5) → Phase 6 forward**

Strukturalne wyniki Phase 5:
- Eikonal dispersion $g^{\mu\nu}_{eff}k_\mu k_\nu = 0$ sympy LOCK z proper Minkowski contractions
- Anisotropic $c_{eff}^2(\theta) = 1 - (\xi/\Lambda^2)n^2\cos^2\theta$ LOCKED
- **Sagnac chopper differential** SNR realistic (~10⁻²⁴ lab, 10⁻¹⁸ ekstremalny) — ≪ v1's fałszywe 3×10⁴
- Yukawa Greens: tensor source ma exp(-Λr) cutoff (no long-range force)
- L₅'_a UNIQUELY produces directional Sagnac signal (reszta kandydatów wyklucza się)

**Cross-cycle implication:** ψ.1.v2 dostarcza **anisotropowy sub-leading kanał**
σ.1 birefringence — TT19-TT23 będą skorygowane predykcje (NIE leading-order
detectable lab effect, tylko astrophysical TOF i ω-INDEPENDENT FRB consistency
z σ.1 + ψ.1 superposition).

**Realistic detection landscape:**
- Lab Sagnac: sub-detection (artifact'em v1 było)
- Magnetar FRB: sub-leading TOF anisotropy on top of σ.1 leading birefringence
- Cosmological: NULL (no large-scale anisotropy from local ∂lnX)
- BBN/CMB: NULL (anisotropic averages to zero on cosmological scales)
