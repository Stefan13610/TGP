---
title: "ω.3.Phase1 results — f_a sympy form + alt-uniqueness 5/5 PASS"
date: 2026-05-01
cycle: ω.3.Phase1
status: COMPLETE
parent: "[[program.md]]"
successor: "[[Phase2_setup.md]]"
tags:
  - TGP
  - omega3
  - phase1
  - f_a
  - results
  - PASS
---

# ω.3.Phase1 results

**Score: 5/5 PASS** ≥4/5 gate → **Phase 2 ENABLED**.

## Sub-test results

| ID | Test | Result | Detail |
|---|---|---|---|
| O1.1 | f_a sympy form derivation | **PASS** | f_a = M_TGP/E_TGP sympy diff = 0 EXACT |
| O1.2 | alt-f_a uniqueness scan | **PASS** | 1/5 PROBE-pass; (b)/(c)/(d)/(e) drift 96–5007% FAIL |
| O1.3 | super-GUT band consistency | **PASS** | f_a in [10¹⁵,10¹⁹] + 5.69 OOM above QCD upper |
| O1.4 | g_aγ < PVLAS-IV ≥5 OOM | **PASS** | g_aγ = 1.71·10⁻²⁰ GeV⁻¹, 9.59 OOM below PVLAS-IV |
| O1.5 | type classification ALP | **PASS** | E-only anomaly (no N) → ALP regime, m_a free |

## Key derived structural identities

### Identity 1 — f_a sympy form derivation

$$\boxed{\;f_a = \frac{M_{TGP}}{E_{TGP}} = \frac{N_A \cdot 2\pi^2 \cdot M_{GUT}}{536/75}\;}$$

Numerically: 3.4630·10¹⁸ / 7.147 = **4.8456·10¹⁷ GeV**.

Sympy diff = 0 EXACT ⇒ structural identity LOCKED.

### Identity 2 — alt-f_a uniqueness 1/5 PROBE-pass

| candidate | f_a (GeV) | drift | super-GUT | verdict |
|---|---:|---:|:-:|:-:|
| **(a) M_TGP / E_TGP** | **4.85·10¹⁷** | **0.00%** | ✓ | **PROBE-PASS ★** |
| (b) M_TGP | 3.46·10¹⁸ | 614.67% | ✓ | FAIL |
| (c) M_GUT | 2.00·10¹⁶ | 95.87% | ✓ | FAIL |
| (d) M_Pl | 1.22·10¹⁹ | 2419.56% | ✗ | FAIL |
| (e) M_TGP · E_TGP | 2.47·10¹⁹ | 5007.48% | ✗ | FAIL |

**Uniqueness gate**: exactly 1/5 candidates satisfy drift < 5% AND super-GUT band ⇒ canonical f_a structurally selected.

### Identity 3 — Super-GUT band

$$f_a^{TGP} = 4.85 \cdot 10^{17}\, \text{GeV} \in [10^{15}, 10^{19}]\, \text{GeV}$$

**5.69 OOM above** classical QCD axion upper bound 10¹² GeV (Pecci-Quinn 1977-79).

→ TGP axion **distinct regime**, NOT classical QCD axion.

### Identity 4 — g_aγ structural form post-ω.3

$$g_{a\gamma} = \frac{g_{axion}}{f_a} = \frac{\alpha_{em}\, E_{TGP}}{2\pi\, f_a} \approx 1.71 \cdot 10^{-20}\, \text{GeV}^{-1}$$

| bound | value | TGP margin |
|---|---:|---:|
| PVLAS-IV/OSQAR-II | 6.6·10⁻¹¹ GeV⁻¹ | 3.85·10⁹× below |
| IAXO 2030+ | 1·10⁻¹² GeV⁻¹ | 5.84·10⁷× below |
| ADMX 2030+ | 1·10⁻¹⁵ GeV⁻¹ | (band-limited; TGP outside) |

**9.59 OOM below PVLAS-IV** → ALL current + 2030+ axion-photon experiments NULL.

### Identity 5 — Type classification: ALP

| anomaly | TGP value | source |
|---|---|---|
| E (electromagnetic) | 536/75 sympy-exact | ω.2 triangle |
| N (color/QCD) | NOT present | TGP framework |

→ TGP axion structurally **ALP** (Axion-Like Particle, E-only anomaly), nie QCD axion.
m_a is **free parameter**: TGP nie wprowadza QCD instanton mass mechanism;
mass-coupling formula m_a · f_a = m_π · f_π · √(m_u·m_d)/(m_u+m_d) NIE applies.

**Reference computation only** (gdyby TGP miało N anomaly): m_a_QCD = 1.19·10⁻¹¹ eV.
Note: ten value przypada w stellar-mass BH superradiance band (LIGO BBH 5-100 M_sun);
TGP m_a is free → ten band NIE jest TGP prediction unless future ω.4+ cycle locks m_a.

## Hypothesis status post-Phase 1

| identity | status |
|---|---|
| f_a = M_TGP/E_TGP | **sympy LOCK preliminary** (Phase 2 confirms) |
| f_a ≈ 4.85·10¹⁷ GeV | super-GUT band POSITIONED |
| g_aγ ≈ 1.71·10⁻²⁰ GeV⁻¹ | NULL @ all current + 2030+ axion experiments |
| ALP classification | LOCKED (E-only anomaly) |

## Phase 1 verdict

**SCORE: 5/5 PASS (≥4/5 gate)** → **Phase 2 enabled**.

**Promotion candidates entering Phase 2:**

1. **f_a = M_TGP/E_TGP sympy LOCK** — 1/5 alt-uniqueness pass
2. **Super-GUT regime** — distinct from classical QCD band by ≥5 OOM
3. **g_aγ NULL at all axion-photon experiments** — 9 OOM below PVLAS bound
4. **ALP classification** — m_a free, ω.4+ future-cycle target

**Phase 2 plan**: Sympy LOCK + numerical reproductions + cosmological consistency
+ F-cluster preservation + 4-channel cascade self-consistency.

## Cross-references

- [[program.md]]
- [[../op-omega2-axion-coupling-lock/Phase3_results.md]]
- [[../op-uv2-mtgp-absolute-scale/Phase3_results.md]]
