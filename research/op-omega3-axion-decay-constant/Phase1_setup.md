---
title: "ω.3.Phase1 setup — f_a structural form derivation + alt-uniqueness"
date: 2026-05-01
cycle: ω.3.Phase1
status: SETUP
parent: "[[program.md]]"
tags:
  - TGP
  - omega3
  - phase1
  - f_a
  - structural-setup
  - setup
---

# ω.3.Phase1 setup

**Score gate:** ≥4/5 PASS = Phase 2 ENABLED.

## Sub-tests (5)

- **O1.1** **f_a sympy form derivation**:
  Inversion z ω.2 LOCK g_axion = α_em·E_TGP/(2π) + TGP-canonical g_aγ relation:
  $$g_{a\gamma} = \frac{\alpha_{em} E_{TGP}}{2\pi f_a} \quad\Rightarrow\quad f_a = \frac{M_{TGP}}{E_{TGP}}$$
  - Sympy diff(f_a_target − f_a_derived) = 0 EXACT requirement.

- **O1.2** **Alt-f_a candidate uniqueness scan**:
  - (a) f_a = M_TGP/E_TGP (TGP canonical winner expected)
  - (b) f_a = M_TGP (no E_TGP suppression)
  - (c) f_a = M_GUT (gauge-unification scale)
  - (d) f_a = M_Pl (Planck-scale)
  - (e) f_a = M_TGP · E_TGP (multiplied — wrong scaling)
  - **Uniqueness gate**: only (a) satisfies dimensional + ω.2 g_aγ inversion + super-GUT band.

- **O1.3** **Super-GUT band consistency**:
  - TGP f_a ≈ 4.85·10¹⁷ GeV ∈ [10¹⁵, 10¹⁹] GeV super-GUT regime.
  - Classical QCD axion band: f_a ∈ [10⁹, 10¹²] GeV (Pecci-Quinn 1977-79).
  - **Distinct regime gate**: TGP f_a > 10¹⁵ GeV by ≥3 OOM od QCD upper limit.

- **O1.4** **g_aγ structural form post-ω.3**:
  $$g_{a\gamma} = \frac{g_{axion}}{f_a} = \frac{\alpha_{em} E_{TGP}}{2\pi f_a} \approx \frac{8.30 \cdot 10^{-3}}{4.85 \cdot 10^{17}} \approx 1.71 \cdot 10^{-20}\, \text{GeV}^{-1}$$
  - vs PVLAS-IV/OSQAR-II current bound 6.6·10⁻¹¹ GeV⁻¹ → margin 3.9·10⁹×.
  - vs IAXO 2030+ projected 10⁻¹² GeV⁻¹ → margin 5.8·10⁷×.
  - **Bound gate**: TGP g_aγ < PVLAS-IV experimental bound by ≥5 OOM (NULL prediction at all current/2030+ experiments).

- **O1.5** **Type classification ALP vs QCD axion**:
  - ω.2 derives g_axion z **electromagnetic** triangle anomaly E_TGP only (no color N).
  - QCD axion: m_a · f_a = m_π · f_π · √(m_u·m_d)/(m_u+m_d) ≈ 5.7·10⁻³ GeV² (anomaly N).
  - TGP ALP: m_a NIE jest mass-coupled do f_a (no QCD instanton mechanism in framework).
  - **Classification gate**: TGP axion structurally ALP (E-only anomaly), m_a free or future-cycle target (ω.4+); not QCD axion.

## Inputs (LOCKED)

```
g_axion (ω.2 LOCK)  = α_em·E_TGP/(2π) ≈ 8.30·10⁻³
E_TGP (ω.2)         = 536/75 ≈ 7.147
M_TGP (UV.2)        = N_A·2π²·M_GUT ≈ 3.4630·10¹⁸ GeV
M_GUT (input)       = 2.0·10¹⁶ GeV (SM 2-loop)
g* (UV.1)           = 71/100
N_A (ξ.1)           = 500/57
α_em                = 1/137.036
M_Pl (PDG)          = 1.221·10¹⁹ GeV
```

## Strategy

**Forward inversion**: ω.2 fixed g_axion = α_em·E_TGP/(2π); UV.2 fixed M_TGP. The
TGP-canonical axion-photon relation g_aγ = g_axion/f_a then identifies f_a
uniquely jako structural ratio M_TGP/E_TGP. Phase 1 confirms this inversion z
sympy LOCK + alt-uniqueness + dimensional consistency.

## Cross-references

- [[program.md]]
- [[../op-omega2-axion-coupling-lock/Phase3_results.md]]
- [[../op-uv2-mtgp-absolute-scale/Phase3_results.md]]
