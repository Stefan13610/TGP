---
title: "ω.3 program — Axion decay constant f_a structural derivation post-UV.2"
date: 2026-05-01
cycle: ω.3
status: ACTIVE
type: mini-cycle
parent: "[[../op-uv2-mtgp-absolute-scale/program.md]]"
predecessors:
  - "[[../op-omega2-axion-coupling-lock/Phase3_results.md]]"
  - "[[../op-uv2-mtgp-absolute-scale/Phase3_results.md]]"
tags:
  - TGP
  - omega3
  - axion
  - f_a
  - decay-constant
  - super-GUT
  - ALP
  - active
---

# ω.3 program — f_a axion decay constant structural derivation

## Hypothesis

$$\boxed{\;f_a = \frac{M_{TGP}}{E_{TGP}} = \frac{N_A \cdot 2\pi^2 \cdot M_{GUT}}{536/75}\;}$$

**Numerically:**
$$f_a = \frac{173.15 \cdot 2.0 \cdot 10^{16}}{7.147} \approx 4.847 \cdot 10^{17}\, \text{GeV}$$

**Super-GUT axion** scale, distinct from classical QCD axion band [10⁹, 10¹²] GeV
(Pecci-Quinn 1977-79).

## Inputs (LOCKED)

| anchor | value | source |
|---|---:|---|
| g_axion | α_em·E_TGP/(2π) ≈ 8.30·10⁻³ | ω.2 LOCK |
| E_TGP | 536/75 sympy-exact | ω.2 triangle anomaly |
| M_TGP | N_A·2π²·M_GUT ≈ 3.4630·10¹⁸ GeV | UV.2 K-LOCK |
| M_GUT | 2.0·10¹⁶ GeV (SM 2-loop) | external (gauge unification) |
| g* | 71/100 | UV.1 NGFP |
| N_A | 500/57 | ξ.1 photon-ring |

## Strategy

UV.2 fixes M_TGP DERIVED FULL z M_GUT external anchor; ω.2 fixes E_TGP = 536/75
TGP triangle anomaly. Inversion z **g_aγ = α_em·E_TGP/(2π·f_a)** TGP-canonical
relation gives:

$$f_a = \frac{M_{TGP}}{E_{TGP}}$$

**Type classification**: TGP axion jest **ALP** (Axion-Like Particle), nie QCD axion:
- ω.2 derives g_axion z **electromagnetic** triangle anomaly E_TGP (E·F·F̃)
- **No QCD color anomaly N** appears in TGP framework
- → m_a NIE jest mass-coupled do f_a via QCD instanton (m_a · f_a ≠ m_π · f_π · √(m_u m_d)/(m_u+m_d))
- m_a może być free parameter lub require additional ω.4+ cycle

**Cyrkularność broken**: f_a derived structurally z M_TGP_UV.2 + E_TGP_ω.2 — **no
PVLAS / ADMX / CAST anchor required**. ω.3 KEYSTONE: f_a TGP-native via UV.2 + ω.2.

## Phase structure

### Phase 1 — Structural setup (5 sub-tests, gate ≥4/5)

- **O1.1** f_a sympy form derivation: f_a = M_TGP/E_TGP z dim-less ratio
- **O1.2** Alt-f_a candidate uniqueness scan (5 alternatives FALSIFIED expected)
- **O1.3** Super-GUT band consistency: f_a ∈ [10¹⁵, 10¹⁹] GeV vs classical QCD band [10⁹, 10¹²] GeV — distinct regime
- **O1.4** g_aγ structural form post-ω.3: g_aγ = g_axion/f_a EXACT → numerical 1.71·10⁻²⁰ GeV⁻¹
- **O1.5** Type classification: ALP (E-anomaly only, no QCD N) — m_a decoupled from f_a, free parameter or future-cycle target

### Phase 2 — Sympy LOCK + numerics (7 sub-tests, gate ≥6/7)

- **O2.1** f_a sympy LOCK = (N_A·2π²·M_GUT)/(536/75)
- **O2.2** f_a numerical = 4.847·10¹⁷ GeV (drift target < 1% vs UV.2 M_TGP × 1/E_TGP)
- **O2.3** g_aγ numerical = α_em·E_TGP/(2π·f_a) ≈ 1.71·10⁻²⁰ GeV⁻¹ (way below PVLAS-IV bound 6.6·10⁻¹¹)
- **O2.4** Cosmological consistency: m_a free → ALP misalignment generic; if pre-inflation PQ + θ_i ~ O(1) → Ω_a h² depends on m_a (decoupled)
- **O2.5** Isocurvature constraint: if pre-inflation PQ + ALP θ_i ~ O(1), CMB Δ_iso < 10⁻¹¹ requires H_inf < f_a/(2π·θ_i) ~ 8·10¹⁶ GeV — easily satisfied
- **O2.6** F-cluster post-ω.3 preservation (F4/F5/F6/XS1 untouched)
- **O2.7** 4-channel ω.3 cascade self-consistency (UV.1 g* + ξ.1 N_A + UV.2 K + ω.2 E_TGP wszystkie flow do f_a)

### Phase 3 — Predictions + 4-channel convergence (6 sub-tests, gate ≥5/6)

- **O3.1** PVLAS-V 2030+ Light-Shining-Through-Walls NULL prediction (g_aγ ~ 10⁻²⁰ GeV⁻¹ ≪ PVLAS-V sensitivity ~10⁻¹¹)
- **O3.2** ALPS-II 2030+ haloscope NULL (TGP super-GUT outside ALP band)
- **O3.3** CAST/IAXO 2030+ solar axion NULL (g_aγ < IAXO bound 10⁻¹² by 8 OOM)
- **O3.4** CMB cosmological birefringence cross-check β = (g_aγ/2)·∫(∂(ln X)/∂η) dη — already locked z ω.2 + UV.2 → ω.3 reproduces consistent prediction
- **O3.5** ALP DM cosmology: m_a free → no specific TGP DM band prediction (open ω.4+ target)
- **O3.6** 4-channel ω.3 convergence summary (UV.2 M_TGP + ω.2 E_TGP + ω.1 g_axion + cosmological NULL bounds)

## Score gates

| Phase | Gate | Promotion |
|---|---|---|
| 1 | ≥4/5 PASS | Phase 2 enabled |
| 2 | ≥6/7 PASS | Phase 3 enabled |
| 3 | ≥5/6 PASS | **ω.3 program END** (f_a DERIVED FULL post-UV.2 + ω.2) |

## Expected outputs

1. **f_a = M_TGP/E_TGP sympy LOCK** TGP-native structural derivation (UV.2 + ω.2)
2. **f_a ≈ 4.85·10¹⁷ GeV super-GUT band** distinct od QCD axion band
3. **g_aγ ≈ 1.71·10⁻²⁰ GeV⁻¹** structurally LOCKED (sub-PVLAS by 9 OOM)
4. **All current axion experiments NULL prediction** (PVLAS, ADMX, ALPS, CAST, IAXO)
5. **m_a free parameter ALP regime** — open future cycle target (ω.4+)

## Cross-references

- [[../op-omega2-axion-coupling-lock/Phase3_results.md]] — ω.2 g_axion + E_TGP
- [[../op-uv2-mtgp-absolute-scale/Phase3_results.md]] — UV.2 M_TGP DERIVED FULL
- [[../op-uv-as-ngfp/Phase3_results.md]] — UV.1 g* NGFP
- [[../op-xi-photon-ring/Phase3_results.md]] — ξ.1 N_A
- [[../op-omega1-substrate-em-coupling/Phase3_results.md]] — ω.1 axion-photon coupling
- [[../../INDEX.md]]
- [[../../PREDICTIONS_REGISTRY.md]]
