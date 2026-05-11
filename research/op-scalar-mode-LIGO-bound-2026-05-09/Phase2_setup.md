---
title: "Phase 2 setup — σ-cross cancellation explicit derivation (R5 risk resolution)"
date: 2026-05-09
parent: "[[./README.md]]"
type: phase-setup
phase: 2
status: 🟢 OPEN
predecessor: "[[./Phase1_results.md]] (R5 risk identified, h_S/h_T ≈ 2√π naive)"
---

# Phase 2 setup — σ-cross cancellation explicit derivation

## §0 — Goal

Resolve R5 risk z Phase 1: explicit derivation of scalar polarization
amplitude h_S in emergent-metric framework, z explicit σ-cross-coupling
contribution.

**Decision:** which scenario operates?

| Scenario | Mechanism | Falsifier verdict |
|---|---|---|
| 1 (R5 realized) | h_S/h_T ≈ 2√π ≈ 3.5 | LIGO polarization tests FALSIFY framework |
| 2 (σ-cross suppression) | σ-cross redirects energy to TT mode | h_S << h_T, framework PASSES |

## §1 — Strategy

### §1.1 — Linearized δΦ-EOM with binary source

Setup z [[../op-emergent-metric-from-interaction-2026-05-09/Phase5_sympy.py]]:
```
(∂_t² - c²∇² + m_sp²) δΦ = q · ρ_source / (K_1 · Φ_0)
```

Z m_sp ≪ ω_LIGO (cosmological scale): effectively massless propagation.

For binary BBH (m_1, m_2 at orbit z separation r_12):
```
ρ_source = m_1 δ³(x - x_1(t)) + m_2 δ³(x - x_2(t))
```

### §1.2 — Far-field δΦ amplitude

Standard quadrupole-formula derivation:
```
δΦ_far(r, t) ~ (q / (4π K_1 Φ_0 r)) · (M_total + (1/2c²) d²Q/dt² + ...)
              [Newton]   [quadrupole radiation, retarded]
```

Quadrupole moment of mass distribution:
```
Q^Φ_ij = Σ_k m_k (x_k)_i (x_k)_j   (mass quadrupole)
```

For circular binary at frequency ω: d²Q^Φ/dt² ~ μ a² ω² (μ = reduced mass).

### §1.3 — Decomposition into δg_eff modes

Z emergent-metric ansatz:
```
δg_eff^00 = -a_1·h,    h = δΦ/Φ_0
δg_eff^ij = δ^ij·b_1·h + σ^ij·c_0/(Φ_0²·c²)
```

**Critical:** σ^ij = (∂_iΦ)(∂_jΦ) - (1/3)δ^ij(∇Φ)² jest BILINEAR in δΦ.
Z far-field δΦ ~ (q/r)·d²Q/dt² ~ ω²/r, σ ~ ω⁴/r².

### §1.4 — Mode decomposition

| Mode | Structure | Origin | Frequency scaling |
|---|---|---|---|
| **Scalar (trace)** | tr(δg^ij)/3 = b_1·h | linear δΦ | h ~ ω²/r |
| **Tensor (TT)** | σ^ij projection | bilinear (∂Φ)² | σ ~ ω⁴/r² |

**Key insight:** scalar mode h_S ~ ω²/r (1/distance fall-off, far-field standard).
Tensor (σ-induced) mode h_T^σ ~ ω⁴/r² (1/r² fall-off — NOT radiation, near-field!).

Dla detector at d_L ~ 400 Mpc: σ-induced TT mode falls as 1/r² (NOT proper
radiation). Only scalar mode reaches detector at radiation rate 1/r.

⟹ **σ-coupling DOES NOT generate radiative tensor mode** at this order.
This means **R5 risk REAL** unless higher-order mechanism operates.

## §2 — Refined Phase 2 plan (post-§1 reflection)

Original plan was: σ-cross cancels scalar at LIGO. But §1.4 analysis shows
σ-coupling enters at 1/r² (near-field), NOT 1/r (radiation). This is different
from what I expected.

### §2.1 — Re-examination: what σ-coupling actually provides

σ-coupling C(ψ) modifies g_eff_ij locally (near binary, where σ ≠ 0). This
changes:
- Orbital binding energy (Phase 3 emergent-metric: Δe_2^σ = c_0·κ_σ)
- Quadrupole emission rate (radiation reaction modified)
- BUT does NOT directly source radiative tensor mode at 1/r

### §2.2 — Resolution: σ effect is on RADIATION REACTION, not radiated polarization

If σ-coupling modifies orbital quadrupole evolution, it changes h_T amplitude
INDIRECTLY (via different orbital dynamics). Phase 4 already accounted for
this (β_ppE^new with σ contribution).

But scalar mode h_S amplitude depends on q·M coupling (linear in δΦ), which
is NOT modified by σ-coupling at leading order.

⟹ **R5 risk PERSISTS at linearized level.**

### §2.3 — Genuine suppression mechanisms (post-Phase 1)

| Mechanism | Where it operates | Realistic? |
|---|---|---|
| **(a) Scalar charge cancellation in binary** | Universal scalar charge q_i = q·m_i means dipole = 0 (COM frame), but quadrupole ≠ 0 | NO — quadrupole survives |
| **(b) Mass Yukawa Phi-decay** | m_sp ≈ H_0 ⟹ Yukawa suppression ineffective at LIGO | NO |
| **(c) Vainshtein screening** | Active in strong-field BBH, suppresses scalar by (r/r_V)^α | **YES** — most likely realistic |
| **(d) Gauge-fixing in TGP** | Maybe ψ has constraint structure removing scalar mode | UNKNOWN — requires careful analysis |

### §2.4 — Pivot: focus on Vainshtein analog

Phase 2 will focus on **Vainshtein-style screening estimate** dla TGP
emergent-metric framework. Specifically:

```
Vainshtein radius r_V dla TGP:
  r_V ~ (G·M)^(1/3) · m_sp^(-2/3)        (massive-graviton analog)
  
For BBH source z M_total ~ 30 M_sun, m_sp ~ H_0 ≈ 10^(-33) eV:
  r_V ~ (G·30M_sun)^(1/3) · H_0^(-2/3)
      ~ (Schwarzschild radius)^(1/3) · (Hubble radius)^(2/3)
      ~ 10^(11) m · 10^(28) m^(2/3) ... need careful calc
```

If r_V > observation distance r_obs (Earth-source), Vainshtein suppression
inactive at LIGO band. If r_V < r_obs, scalar mode suppressed.

## §3 — Phase 2 deliverables

- Phase2_sympy.py — explicit calculation
- Phase2_sympy.txt
- Phase2_results.md — verdict on R5 risk

## §4 — Phase 2 gate criteria

| # | Criterion |
|---|---|
| G1 | Linearized scalar mode amplitude h_S derived sympy |
| G2 | TT mode (z σ-coupling) amplitude h_T at far-field derived |
| G3 | h_S/h_T ratio numerical |
| G4 | Vainshtein radius estimate dla TGP |
| G5 | Comparison z LIGO bound (5%) z honest verdict |
| G6 | If R5 realized: identify what TGP framework needs |

## §5 — Probability assessment

| Outcome | Probability |
|---|---|
| Scenario 2 (suppression PASSES): | 25-40% (Vainshtein analog plausible) |
| Scenario 1 (R5 REALIZED): | 40-55% (naive analysis confirmed) |
| Inconclusive (multi-session work needed): | 20-30% |

## §6 — Cross-references

- [[./Phase1_results.md]] — R5 risk identification
- [[../op-emergent-metric-from-interaction-2026-05-09/Phase5_sympy.py]] — linearized Phi-EOM
- [[../op-emergent-metric-from-interaction-2026-05-09/Phase4_results.md]] — N14 deferral
- [[../closure_2026-04-26/sigma_ab_pathB/]] — decoupling regime + σ_ab Path B
