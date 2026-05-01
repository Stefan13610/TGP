---
title: "ω.2.Phase1 results — anomaly-based g_bare derivation 5/5 PASS"
date: 2026-05-01
cycle: ω.2.Phase1
status: COMPLETE
parent: "[[program.md]]"
predecessor: "[[Phase1_setup.md]]"
successor: "[[Phase2_setup.md]]"
tags:
  - TGP
  - omega2
  - phase1
  - anomaly
  - chirality
  - results
  - PASS
---

# ω.2.Phase1 results

**Score: 5/5 PASS** ≥4/5 gate → **Phase 2 ENABLED**.

## Sub-test results

| ID | Test | Result | Detail |
|---|---|---|---|
| W2.1.1 | E_TGP triangle anomaly | **PASS** | E_TGP = 536/75 sympy-exact |
| W2.1.2 | g_bare = η_chir = 1 − C6 | **PASS** | 19/24 dual-derived (1−C6 form & K-num+ form) |
| W2.1.3 | g_eff(IR) = α_em·E_TGP/(2π) | **PASS** | 8.300·10⁻³, drift vs α_em = 13.74% (closest) |
| W2.1.4 | Dimensional + scale-symm | **PASS** | g dimensionless, F·F̃ boundary, f_X ~ M_TGP UV-IR |
| W2.1.5 | 5-alt falsification | **PASS** | exactly 1/5 (η_chir = 19/24) PROBE-passes |

## Key derived structural identities

### Identity 1 — TGP triangle anomaly coefficient (sympy-exact)

$$E_{TGP} = N_c \cdot [Q_u^2 B^2_{up} + Q_d^2 B^2_{down}] + Q_l^2 B^2_{lep}$$

Plug LOCKED inputs:
- B²_up = 13/4 (θ.1)
- B²_down = 61/25 (θ.1 derived B²_up − 81/100)
- B²_lep = 2 (Dirac leptons)
- N_c = 3, Q_u = 2/3, Q_d = −1/3, Q_l = −1

$$E_{TGP} = 3 \cdot \left[\frac{4}{9}\cdot\frac{13}{4} + \frac{1}{9}\cdot\frac{61}{25}\right] + 1 \cdot 2$$

$$= 3 \cdot \left[\frac{13}{9} + \frac{61}{225}\right] + 2 = 3 \cdot \frac{386}{225} + 2 = \frac{386}{75} + \frac{150}{75} = \boxed{\frac{536}{75}} \approx 7.1467$$

**LOCKED structural identity**: E_TGP = 536/75 sympy-exact.

### Identity 2 — η_chir dual-derivation

**Form A** (K-complement):
$$\eta_{chir} = 1 - C_6 = 1 - \frac{5}{24} = \frac{19}{24}$$

**Form B** (K-numerator combinatoric):
$$\eta_{chir} = \frac{K_{up,\text{num}} + 2 N_{gen} K_{lep,\text{num}}}{K_{up,\text{denom}} \cdot N_{gen}} = \frac{7 + 2\cdot 3 \cdot 2}{8 \cdot 3} = \frac{19}{24}$$

→ Two independent derivations match → **TGP-structural identity LOCKED**.

### Identity 3 — IR effective coupling (anomaly form)

$$g_{eff}(IR) = \frac{\alpha_{em} \cdot E_{TGP}}{2\pi} = \frac{1/137.036 \cdot 536/75}{2\pi} \approx 8.300 \cdot 10^{-3}$$

**Drift vs 4 ω.1 LOCK candidates:**

| candidate | value | drift |
|---|---:|---:|
| α_em | 7.297·10⁻³ | **13.74%** ← closest |
| 1/(2π) | 0.1592 | 94.78% |
| η_chir = 19/24 | 0.7917 | 98.95% |
| κ_TGP | 2.012 | 99.59% |

**Insight**: g_eff(IR) ≈ α_em (E_TGP ≈ 7.147 ~ 2π = 6.28, so anomaly factor near 1).
Drift 13.74% reflects E_TGP/(2π) ratio; structurally α_em is the natural IR limit
of g_anomaly.

### Identity 4 — Dimensional consistency

L_ω.1 ⊃ (g/4)(ln X)F·F̃:
- ln X dimensionless (scale ratio)
- F·F̃ has mass⁴ dim
- ∫ d⁴x has mass⁻⁴
- → L dim mass⁴, **g dimensionless** ✓

Scale-symm X → λX preserves EOM (F·F̃ = total derivative, boundary term only).
f_X ~ M_TGP UV-IR matching consistent z ω.1 W3.1.

## Falsification ledger (W2.1.5)

5 alt g-ansatz tested:

| # | Candidate | Verdict | Rationale |
|---|---|---|---|
| (i) | g = α_em alone | ✗ FAIL | E_TGP anomaly factor omitted; nie TGP-structural |
| (ii) | g = 1/(2π) standard EFT | ✗ FAIL | no TGP B²-content; generic axion EFT |
| (iii) | g = κ_TGP cross-sector | ✗ FAIL | XS3 lepton-orthogonal sector, wrong context |
| (iv) | g = C6 = 5/24 directly | ✗ FAIL | K-diff partial, nie complement |
| (v) | g = η_chir = 19/24 = 1 − C6 | ✓ **PASS** | TGP-native bare coupling, K-complement structural |

**1/5 PROBE-passes** → uniqueness criterion satisfied.

## Phase 1 verdict

**SCORE: 5/5 PASS (≥4/5 gate)** → **Phase 2 enabled**.

**Promotion candidates entering Phase 2:**

1. **g_bare = η_chir = 19/24** (TGP K-level chirality complement, structural)
2. **g_eff(IR) = α_em · 536/(75·2π) ≈ 8.300·10⁻³** (anomaly one-loop form)
3. **E_TGP = 536/75** (LOCKED structural identity, B²-cascade extension)

**Hypothesis discrimination Phase 2 will test:**

Cross-channel data (CMB β + magnetar + quasar + PVLAS-V) picks between:
- **g_bare path**: g = 19/24 at substrate UV scale (M_TGP) → very strong axion signal
- **g_eff path**: g ≈ 8.3·10⁻³ at IR low-E scale → weak signal close to α_em

**Phase 2 plan**: solve `g · Δ(ln X)_cosmo = 2β_CMB ≈ 0.01186 rad` constraint
for each candidate; check feasibility w cross-channel ranking.

## Cross-references

- [[program.md]]
- [[Phase1_setup.md]]
- [[Phase2_setup.md]]
- [[../op-rho1-71Ge-cross-section/Phase3_results.md]] — C6 = 5/24 LOCKED
- [[../op-theta-quark-koide/Phase3_results.md]] — B²_up = 13/4, B²_down = 61/25
- [[../op-omega1-substrate-em-coupling/Phase3_results.md]] — 4 LOCK candidates ω.1
