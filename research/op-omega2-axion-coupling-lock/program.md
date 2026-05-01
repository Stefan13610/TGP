---
title: "ω.2 program — axion-substrate coupling g unique selection (4 candidates → 1 LOCK)"
date: 2026-05-01
cycle: ω.2
status: ACTIVE
parent: "[[../op-omega1-substrate-em-coupling/Phase3_results.md]]"
tags:
  - TGP
  - omega2
  - axion-coupling
  - g-lock
  - chirality-anomaly
  - cross-channel-falsification
---

# ω.2 — axion-substrate coupling g unique selection

**Goal:** Po ω.1 mamy 4 strukturalnie clean LOCK candidates dla
g w `L_ω.1 ⊃ (g/4)(ln X)F·F̃`:

| candidate | wartość | pochodzenie |
|---|---:|---|
| g = κ_TGP | 2.012 | XS1 cross-sector identity √α₀ = κ_TGP |
| g = α_em | 7.297·10⁻³ | EM coupling 1/137.036 |
| g = 1/(2π) | 0.1592 | standard axion EFT normalization |
| g = η_chir = 19/24 | 0.7917 | TGP K-level chirality factor (1 − C6) |

ω.2 łączy **strukturalną derywację anomaly-based** (chirality counting
B²_up, B²_down, B²_lep) z **cross-channel falsification** (CMB β, magnetar
E·B, quasar Δχ, vacuum birefringence) by jednoznacznie wybrać g.

**Score gate:** ≥4/5 + ≥6/7 + ≥5/6 = 18 sub-tests, 4-channel convergence
≥3/4 → ω.2 program END (g LOCKED).

## Hypothesis

Strukturalna anomaly-based prediction: **g_TGP = η_chir = 19/24**
(bare-coupling przy substrate scale, K-level chirality complement
1 − C6 z C6 = K_up − K_lep = 5/24 LOCKED).

Alternative low-energy effective form po wycałkowaniu fermionów:

$$g_{\text{eff}}(IR) = \frac{\alpha_{em} \cdot E_{TGP}}{2\pi}$$

gdzie $E_{TGP} = N_c \cdot [Q_u^2 \cdot B^2_{up} + Q_d^2 \cdot B^2_{down}]
+ Q_l^2 \cdot B^2_{lep}$ jest TGP-native EM anomaly coefficient.

Plug B²_up = 13/4, B²_down = 61/25, B²_lep = 2, Q_u = 2/3, Q_d = -1/3, Q_l = -1:

$$E_{TGP} = 3 \cdot \left[\frac{4}{9}\cdot\frac{13}{4} + \frac{1}{9}\cdot\frac{61}{25}\right] + 2
= 3 \cdot \frac{386}{225} + 2 = \frac{1158}{225} + 2 = \frac{536}{75} \approx 7.147$$

→ g_eff(IR) ≈ α_em · 7.147 / (2π) ≈ 0.00830 (close to α_em, drift ~14%).

**Hypothesis discrimination:** g_bare (UV substrate scale) = 19/24
vs g_eff (IR low-energy effective) ≈ 0.0083. Cross-channel data picks
which one is physically observed.

**Strong structural anchor**: g = η_chir = 19/24 wykorzystuje istniejące
LOCKED quantities (K_up = 7/8, K_lep = 2/3, C6 = 5/24) bez nowych free
params. To czyni g = 19/24 **TGP-prior-favored** structurally.

## Reference frame

| Source | Status quo (post-ω.1) | ω.2 lift |
|---|---|---|
| ω.1 g LOCK | 4 candidates, equal-weighted | unique selection 1/4 |
| CMB β = 0.34±0.09° | LIVE PARTIAL ~3.8σ | g·Δ(ln X) = 0.012 constraint |
| ρ.1 C6 = K_up − K_lep = 5/24 | LOCKED | η_chir = 1 − C6 = 19/24 |
| η.2 α-residual = 9/250 | LOCKED, 2(B²_up−B²_down)/(N_gen²·5) | parallel B²-cascade structural form |
| ι.1 |Q_u|² − |Q_d|² = 1/N_gen | LOCKED | charge-chirality cross-link |
| φ.1 substrate-action | AXIOM `S[X]=½(∂lnX)²` | gives Δ(ln X) cosmological evolution form |
| Phase 2 EFT survival | F5 g̃ = 0.9803 STRUCTURAL | sub-percent EFT corrections |

## Phase plan (5 + 7 + 6 = 18 sub-tests)

### Phase 1 — Anomaly-based g_bare derivation (5 sub-tests)

- **W2.1.1** Triangle anomaly E_TGP from chirality B²: sympy compute
  E_TGP = N_c·[Q_u²·B²_up + Q_d²·B²_down] + Q_l²·B²_lep z B²_up=13/4,
  B²_down=61/25, B²_lep=2 → E_TGP = 536/75 sympy-exact
- **W2.1.2** Bare-coupling structural form: g_bare = η_chir = 1 − C6
  = 1 − 5/24 = 19/24 sympy-exact; cross-check że η_chir = (K_up_num
  + 2·N_gen·K_lep_num)/(K_up_denom·N_gen) = (7+12)/24 = 19/24 alternative
  derivation
- **W2.1.3** IR effective g_eff = α_em·E_TGP/(2π) sympy compute,
  drift vs α_em LOCK candidate; check że g_eff ∈ [α_em, 19/24] band
- **W2.1.4** Dimensional consistency: g dimensionless + scale-symm
  X→λX preserved + f_X ~ M_TGP coupling (UV-IR matching condition)
- **W2.1.5** **Alt-form falsification (5 alts):**
  - (i) g = α_em (no anomaly weighting) — fail E_TGP omitted
  - (ii) g = 1/(2π) (standard axion) — fail no TGP-structural content
  - (iii) g = κ_TGP (cross-sector SC) — fail wrong sector (XS3 ortho)
  - (iv) g = 5/24 = C6 directly — fail wrong sign of K-diff
  - (v) g = η_chir = 19/24 — TEST candidate (PROBE preferred)

**Score gate:** ≥4/5 PASS

### Phase 2 — Sympy LOCK + cross-channel matching (7 sub-tests)

- **W2.2.1** CMB β constraint: g_obs · Δ(ln X)_cosmo = 2β ≈ 0.01186 rad;
  solve dla każdego LOCK candidate co Δ(ln X) wymagane
- **W2.2.2** Δ(ln X) cosmological evaluation z φ.1 substrate-action: solve
  EL eq `□(ln X) = 0` w FRW background z BC X(z=1100), X(z=0); compute
  natural Δ(ln X) ~ ln(z_max/z_now) · suppression_factor
- **W2.2.3** PVLAS-V 2030+ projected sensitivity: g/f_a < 1·10⁻¹¹ GeV⁻¹
  (10× tighter); check că każdy LOCK candidate compatible
- **W2.2.4** Magnetar E·B sensitivity (FAST/SKA 2030+): predict
  pulsar-timing anomaly z E_rot · B / (f_X² · M_TGP²) — ranking dla 4 g
- **W2.2.5** Quasar Δχ ∝ ln(1+z) coefficient lock: SKA 2030+ z>4
  polarimetry; predicted Δχ/ln(1+z) = g/2 dla każdego candidate
- **W2.2.6** Combined-channel χ² ranking: 4 LOCK candidates vs 4 channels
  → unique winner z combined-likelihood
- **W2.2.7** UV-IR matching consistency: g_bare(M_TGP) ↔ g_eff(low E)
  via RG running — czy sympy preserves anomaly relation

**Score gate:** ≥6/7 PASS

### Phase 3 — Predictions + 4-channel convergence (6 sub-tests)

- **W2.3.1** SO 2027+ CMB β POST-CONFIRM target (5σ projected); predict
  central β z g_LOCK + ω.2 Δ(ln X) form
- **W2.3.2** LiteBIRD 2029+ CMB β + B-mode E↔B chirality cross-check
- **W2.3.3** PVLAS-V 2030+ vacuum birefringence forward gate (lab)
- **W2.3.4** FAST/SKA 2030+ magnetar E·B forward gate (astrophysical)
- **W2.3.5** SKA 2030+ quasar Δχ z>4 forward gate (cosmological)
- **W2.3.6** **4-channel ω.2 convergence**:

| # | Channel | Form | Method | Status |
|---|---|---|---|---|
| 1 | Anomaly chirality | g_bare = η_chir = 19/24 | W2.1 sympy | structural |
| 2 | CMB β LIVE PARTIAL | g·Δ(ln X) = 0.01186 | Planck PR4 + ACT 2024 | observational |
| 3 | UV-IR matching | g_eff = α·E/(2π) consistent | RG running | structural |
| 4 | Cross-channel discriminate | unique winner from combined χ² | 4-channel falsification | observational |

**Score gate:** ≥5/6 PASS → **ω.2 program END (g LOCKED)**

## Falsifiability path

ω.2 jest falsyfikowane przez:
1. CMB β central deviation > 5σ od g_LOCK·Δ(ln X)_TGP
2. Multiple LOCK candidates surviving combined χ² (no unique winner) →
   ω.2 NIE rozstrzyga — degeneracja persistent
3. Anomaly E_TGP ≠ 536/75 sympy-exact (B² values revised)
4. PVLAS-V detection inconsistent z g_LOCK
5. Quasar SKA Δχ slope ≠ g_LOCK/2

## Promotions post-ω.2 (jeśli FULL CONVERGENCE)

- **g_ω.1 PARTIALLY DERIVED → DERIVED** (unique selection from 4 candidates)
- **η_chir = 19/24 LOCKED** (analog do C6 = 5/24 LOCKED)
- **B²-cascade extension**: anomaly E_TGP = 536/75 LOCKED structural identity
- **Q-charge × B²-chirality cross-product framework** consistent
- **CMB β LIVE PARTIAL → POST-CONFIRMED** (jeśli SO/LiteBIRD 2027+ confirms)
- Sets stage dla **ω.3 f_a derivation** post-χ.1 M_TGP lock

## Cross-references

- [[../op-omega1-substrate-em-coupling/Phase3_results.md]] — ω.1 4 LOCK candidates
- [[../op-rho1-71Ge-cross-section/Phase3_results.md]] — C6 = 5/24 LOCKED
- [[../op-eta2-denom-derivation/Phase3_results.md]] — η.2 α-residual B²-cascade
- [[../op-iota-charge-pmns-unification/Phase3_results.md]] — ι.1 charge-chirality cross-link
- [[../op-theta-quark-koide/Phase3_results.md]] — θ.1 B²_up=13/4, B²_down=61/25
- [[../op-phi1-substrate-action-variational/Phase3_results.md]] — φ.1 EL eq for Δ(ln X) cosmological
- [[../../INDEX.md]]
- [[../../PREDICTIONS_REGISTRY.md]]
