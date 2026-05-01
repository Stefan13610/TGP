---
title: "ω.3.Phase2 results — sympy LOCK + cosmological + F-cluster 7/7 PASS"
date: 2026-05-01
cycle: ω.3.Phase2
status: COMPLETE
parent: "[[program.md]]"
predecessor: "[[Phase1_results.md]]"
successor: "[[Phase3_setup.md]]"
tags:
  - TGP
  - omega3
  - phase2
  - sympy-lock
  - f_a
  - results
  - PASS
---

# ω.3.Phase2 results

**Score: 7/7 PASS** ≥6/7 gate → **Phase 3 ENABLED**.

## Sub-test results

| ID | Test | Result | Detail |
|---|---|---|---|
| O2.1 | f_a sympy LOCK rigorous | **PASS** | f_a = 3125·π²·M_GUT/1273 sympy-rational diff = 0 |
| O2.2 | f_a numerical drift < 1% | **PASS** | UV.2 vs χ.1 drift 0.30% (inherited K-drift) |
| O2.3 | g_aγ direct vs inversion | **PASS** | sympy diff = 0 EXACT, numerical 3·10⁻³⁶ |
| O2.4 | cosmological ALP super-GUT | **PASS** | f_a > 10¹⁶ GeV, anthropic θ_i / kinetic alignment compatible |
| O2.5 | isocurvature H_inf constraint | **PASS** | H_inf < 3.04·10⁷ GeV (θ_i=1) → r < 3.17·10⁻²³ low-scale inflation |
| O2.6 | F-cluster preservation | **PASS** | max drift 0.083% (XS1), F4/F5/F6 untouched < 0.001% |
| O2.7 | 4-channel cascade diff = 0 | **PASS** | sympy 4-channel cascade self-consistent |

## Key derived structural identities

### Identity 1 — f_a sympy-rational LOCK

$$\boxed{\;f_a = \frac{N_A \cdot 2\pi^2 \cdot M_{GUT}}{E_{TGP}} = \frac{(500/57)\cdot 2\pi^2 \cdot M_{GUT}}{536/75} = \frac{3125\,\pi^2\, M_{GUT}}{1273}\;}$$

Sympy diff = 0 EXACT.

**Physical interpretation**:
- 3125 = 5⁵ (numerator inheritance z N_A = 500/57 + 75 z E_TGP denom)
- 1273 = 19 × 67 (denominator inheritance z 57 z N_A + 536 z E_TGP)
- π² geometric inheritance z S³ unit sphere half-volume

### Identity 2 — f_a numerical (UV.2 vs χ.1 cross-check)

| derivation | f_a (GeV) | drift |
|---|---:|---:|
| **UV.2** (K_struct·M_GUT/E_TGP) | **4.8456·10¹⁷** | reference |
| χ.1 (M_TGP_χ.1/E_TGP) | 4.8602·10¹⁷ | 0.30% |

Drift 0.30% inherited z K-drift (UV.2 vs χ.1 anchor).

### Identity 3 — g_aγ structural cross-check

| derivation | g_aγ (GeV⁻¹) |
|---|---:|
| **direct** α_em·E_TGP/(2π·f_a) | **1.7129·10⁻²⁰** |
| **inversion** g_axion(ω.2)/f_a | 1.7129·10⁻²⁰ |
| **sympy diff** | 0 EXACT |

Algebraic identity: α_em·E²_TGP/(2π·M_TGP) sympy-LOCKED.

### Identity 4 — Cosmological ALP super-GUT consistency

f_a = 4.85·10¹⁷ GeV > 10¹⁶ GeV super-GUT threshold:
- **Anthropic θ_i tuning** scenarios compatible (string-axiverse landscape)
- **Kinetic alignment** mechanisms compatible (Co et al. 2020)
- **m_a free parameter** (ALP regime) → Ω_a h² nie locked z f_a

### Identity 5 — Isocurvature constraint (pre-inflation PQ)

$$\Delta_{iso} = \frac{H_{inf}}{2\pi\, f_a\, \theta_i} < 10^{-11}\, (\text{Planck PR4})$$

Dla TGP f_a + θ_i = 1:
$$H_{inf}^{max} = 2\pi\, f_a\, \theta_i\, \Delta_{iso}^{max} = 3.04 \cdot 10^{7}\, \text{GeV}$$

Tensor-to-scalar ratio bound:
$$r^{max} = \frac{16}{\pi}\left(\frac{H_{inf}^{max}}{M_{Pl}}\right)^2 \approx 3.17 \cdot 10^{-23}$$

→ **Low-scale inflation prediction** (r essentially zero, BICEP/LiteBIRD null
expected) IF pre-inflation PQ + θ_i ~ O(1).

**Alternative scenarios**:
- post-inflation PQ breaking → no isocurvature (T_RH > f_a impossible since f_a > M_GUT, but topological defects)
- θ_i tuning to 10⁻⁴–10⁻³ allows H_inf up to ~10¹⁰ GeV

### Identity 6 — F-cluster preservation post-ω.3

| anchor | post-χ.1 | post-UV.2 | post-ω.3 | drift | status |
|---|---|---|---|---:|:-:|
| F4 (α₀) | 4.0447 | 4.0447 | 4.0447 | 0.0009% | ✓ |
| F5 (g̃) | 0.9803 | 0.9803 | 0.9803 | 0.0000% | ✓ |
| F6 (κ) | √(32π) | √(32π) | √(32π) | 0.0001% | ✓ |
| XS1 (√α₀ ≡ κ_TGP) | 0.042% | 0.042% | 0.083% | 0.0831% | ✓ |

→ ω.3 doesn't perturb F-cluster — max drift 0.083% (XS1 reference shift).

### Identity 7 — 4-channel ω.3 cascade

| Channel | Source | Form | Verdict |
|---|---|---|---|
| 1 | UV.1 | g* = 71/100 NGFP | ✓ exact |
| 2 | ξ.1 | N_A = 500/57 | ✓ exact |
| 3 | UV.2 | K_struct = N_A·2π² | ✓ sympy |
| 4 | ω.2 | E_TGP = 536/75 | ✓ sympy |

Sympy diff cascade = 0 EXACT → 4 anchors all flow consistently into f_a.

## Hypothesis status post-Phase 2

**f_a sympy-rational LOCK form:**
$$\boxed{\;f_a = \frac{3125\,\pi^2}{1273}\, M_{GUT} \approx 24.23\, M_{GUT} \approx 4.85 \cdot 10^{17}\, \text{GeV}\;}$$

**Numerical post-Phase 2:**
- f_a = 4.8456·10¹⁷ GeV (drift 0.30% z χ.1 cross-check)
- g_aγ = 1.7129·10⁻²⁰ GeV⁻¹
- H_inf < 3.04·10⁷ GeV (low-scale inflation if pre-inflation PQ + θ_i=1)
- r < 3.17·10⁻²³ tensor-to-scalar ratio bound

## Phase 2 verdict

**SCORE: 7/7 PASS (≥6/7 gate)** → **Phase 3 enabled**.

**Promotion candidates entering Phase 3:**

1. **f_a = 3125·π²·M_GUT/1273 sympy-rational LOCK** — TGP-native ξ.1 + UV.2 + ω.2
2. **g_aγ = α_em·E²_TGP/(2π·M_TGP) sympy** — algebraic identity
3. **Low-scale inflation forecast** — r < 3.17·10⁻²³ if pre-inflation PQ + θ_i=O(1)
4. **F-cluster + cascade preservation** — max drift 0.083%

**Phase 3 plan**: PVLAS/IAXO/CAST/ADMX null forecasts + CMB cosmology +
4-channel convergence + ALP DM band + ω.3 program END verdict.

## Cross-references

- [[program.md]]
- [[Phase1_results.md]]
- [[Phase2_setup.md]]
- [[../op-omega2-axion-coupling-lock/Phase3_results.md]]
- [[../op-uv2-mtgp-absolute-scale/Phase3_results.md]]
