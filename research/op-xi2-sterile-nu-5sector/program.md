---
title: "ξ.2 program — sterile ν 5-sector extension B²_sterile + reactor/gallium anomaly resolution"
date: 2026-04-30
cycle: ξ.2
status: ACTIVE
parent: "[[../../INDEX.md]]"
predecessor: "[[../op-nu-majorana-phase-mbb/Phase3_results.md]]"
tags:
  - TGP
  - xi2
  - sterile-nu
  - 5-sector
  - reactor-anomaly
  - gallium-anomaly
  - STEREO
  - PROSPECT
  - BEST
  - KATRIN
---

# ξ.2 program — sterile ν 5-sector extension B²_sterile

## Cel

Extend TGP 4-sector chirality framework do 5-sector z B²_sterile;
derive |U_e4|² + Δm²_{41} + sterile mass m_4 z B²-cascade; resolve
short-baseline reactor antineutrino anomaly (RAA) + gallium anomaly
(GALLEX/SAGE/BEST 2022 4σ) z TGP-native prediction; register
STEREO 2023 + PROSPECT-I + MicroBooNE + BEST + KATRIN-TRISTAN
falsifications.

## 4-sector → 5-sector extension

**4-sector baseline (post-ν.1):**

| Sector | B² | Numerator | Status |
|--------|-----|-----------|--------|
| ν (Majorana active) | 1 | — | LOCKED z μ.1/ν.1 (chirality-halving) |
| lep (Dirac charged) | 2 | — | LOCKED z μ.1 (chirality cascade) |
| up | 13/4 | 13 | LOCKED z κ.1 ((ν,up)-pair B²_up_num=13) |
| dn | ? (≈ 78/?) | 78 | LOCKED z η.2/κ.1 (Wolfenstein denom 78) |

**5-sector extension hypothesis:**

| Sector | B²_sterile candidate | Implication |
|--------|---------------------|--------------|
| **A (decoupled)** | B²_sterile = 0 | sterile fully decouples; |U_e4|² → 0; RAA + Gallium NOT sterile |
| **B (half-Majorana)** | B²_sterile = 1/2 | sin²(2θ₁₄) ≈ B²_sterile/B²_lep · λ_C² ≈ 0.013 |
| **C (electroweak-suppressed)** | B²_sterile = 1/4 | sin²(2θ₁₄) ≈ 0.0063 — below STEREO 2023 sensitivity |
| **D (vacuum-only)** | B²_sterile = λ_C² · 1 ≈ 0.0506 | sin²(2θ₁₄) ≈ λ_C⁴ ≈ 0.00256 — null all SBL |

## Hipoteza (A) — sterile sector decouples (B²_sterile = 0)

W TGP NO ordering 4-sector chirality structure z minimalnym algebraic
counting, sterile right-handed neutrino jest **structurally absent** —
nie ma B²-eigenvalue dla sterile sektor (B²_sterile = 0):

$$U_{e4}^{2} = \frac{B^{2}_{\rm sterile}}{B^{2}_{\nu} + B^{2}_{\rm lep}} = \frac{0}{3} = 0$$

→ |U_e4|² = 0 strict ⇒ **no sterile oscillation**.

**Predictions:**
- STEREO 2023: NO sterile signal at Δm² ~ 1 eV² (CONFIRMED, post-prediction)
- PROSPECT-I: NO sterile signal (CONFIRMED, post-prediction)
- BEST 2022 4σ Gallium anomaly: **NOT sterile** — must be from
  ⁷¹Ge cross-section systematics or solar-ν flux modeling
- KATRIN-TRISTAN 2027+: NO m_4 < 1 eV signal
- MicroBooNE 2024+: NO ν_e excess at sterile parameters

## Hipoteza (B) — sterile sector at vacuum suppression (B²_sterile = λ_C²)

Alternative: sterile coupling exists but suppressed by Cabibbo
parameter — analog do CKM mixing-suppression cascade:

$$\sin^{2}(2\theta_{14}) \approx \lambda_{C}^{4} \approx 2.56 \times 10^{-3}$$

→ below STEREO 2023 sensitivity (~0.01 at 95% CL); below PROSPECT-I.

**Predictions:**
- All current SBL experiments: NULL (consistent z 2024-2025 data)
- SBN program (ICARUS+SBND+MicroBooNE) 2030+: potential 2-3σ hint
- KATRIN-TRISTAN 2027+: marginal m_4 < 1 eV signal

## Plan 3-fazowy

### Phase 1 — sterile ν landscape audit + B²_sterile candidates + viability (5 sub-tests)

- X1.1: Reactor antineutrino anomaly (RAA) audit: 2011 Mention et al. 6%
  deficit; 2023 STEREO ruled out at Δm² ~ 1 eV² (95% CL exclusion);
  2023 PROSPECT-I confirmed STEREO; 2024 Daya Bay + RENO updated flux
  reconstruction reduced anomaly to ~2%
- X1.2: Gallium anomaly audit: GALLEX (1995), SAGE (2010), BEST (2022)
  4σ-5σ deficit z ⁵¹Cr/⁵¹Cr⁺ν source; sterile Δm² ~ 1-10 eV², |U_e4|² ~ 0.1
- X1.3: KATRIN sterile bounds: m_4 < 1.6 eV (90% CL) for sin²(2θ_e4) > 0.1
  (KATRIN 2022); TRISTAN 2027+ extends sensitivity to m_4 ~ 0.1 eV
- X1.4: 4 B²_sterile candidates (A 0, B 1/2, C 1/4, D λ_C²) → compute
  predicted |U_e4|² + Δm²_{41} + m_4
- X1.5: Viability gate — wszystkie 4 candidates muszą produkować
  |U_e4|² compatible z STEREO 2023 + PROSPECT-I exclusion (sin²(2θ_14) < 0.01)

### Phase 2 — first-principles B²_sterile + |U_e4|² + Δm²_{41} (7 sub-tests)

- X2.1: Form A (B²_sterile = 0) — sympy-exact null prediction
- X2.2: Form B (B²_sterile = λ_C²) — sympy-exact sin²(2θ_14) = λ_C⁴
- X2.3: Δm²_{41} dla Form A vs B — Form A undefined (no oscillation),
  Form B = m_3 · enhancement factor (RG-stability)
- X2.4: m_4 sterile mass — Form A: 0, Form B: ~ 0.1-1 eV (TGP-native)
- X2.5: 5 alternative forms FALSIFIED:
  (i) sterile-3+1 RAA best fit (Δm² = 1.7 eV², |U_e4|² = 0.02) — fails STEREO
  (ii) gallium 3+1 (Δm² = 4 eV², |U_e4|² = 0.1) — fails STEREO
  (iii) eV-scale ν_R Dirac (m_4 ~ eV) — fails KATRIN
  (iv) light sterile w/ flavor-democratic mixing — fails reactor
  (v) heavy ν_R seesaw (m_4 > 100 GeV) — invisible in SBL/KATRIN, irrelevant
- X2.6: TGP-native Σm_ν 5-sector update: ζ.1 Σm_ν = 59.01 meV unchanged
  if B²_sterile = 0 OR + m_4 sterile if Form B (cosmology constraint)
- X2.7: classification cascade — Form A LOCKED dual: |U_e4|² = 0 LOCKED,
  Δm²_{41} N/A LOCKED, m_4 = 0 LOCKED (3 promotions Form A)
  Form B LOCKED dual: sin²(2θ_14) = λ_C⁴ LOCKED, m_4 ~ 0.1 eV LOCKED

### Phase 3 — predictions + falsification convergence (6 sub-tests)

- X3.1: STEREO 2023 ✓ post-prediction confirmation Form A null signal
- X3.2: PROSPECT-I ✓ post-prediction confirmation Form A null signal
- X3.3: BEST 2022 4σ Gallium anomaly — TGP Form A predicts NOT sterile
  (must be ⁷¹Ge cross-section error); register as **BEST → systematics
  tension** falsifier
- X3.4: KATRIN-TRISTAN 2027+ m_4 sensitivity ~ 0.1 eV
  Form A: predicts NULL m_4 signal
  Form B: predicts marginal m_4 ~ 0.1 eV (1-2σ)
  → discriminator dla Form A vs B
- X3.5: SBN program 2030+ (ICARUS+SBND+MicroBooNE) sin²(2θ) ~ 10⁻³
  sensitivity discriminates Form A (null) vs B (λ_C⁴ ~ 2.56·10⁻³ marginal)
- X3.6: 7-channel ξ.2 falsification convergence — STEREO 2023 + PROSPECT-I
  + BEST 2022 + KATRIN-TRISTAN 2027+ + SBN 2030+ + MicroBooNE 2024-2027 +
  cosmological N_eff (Planck+CMB-S4 2030+)

## Verdict gate

- Phase 1: 5/5 → Phase 2 viable
- Phase 2: 7/7 + classification cascade → Phase 3 viable
- Phase 3: 5/6 PASS minimum → ξ.2 program END; 6/6 → FULL CONVERGENCE

## Cumulative ledger

553 + 5 + 7 + 6 = **571** post-ξ.2 (target).

## Cross-references

- [[../op-nu-majorana-phase-mbb/Phase3_results.md]] — ν.1 closure baseline
- [[../op-mu-pmns-phase-hardening/Phase3_results.md]] — μ.1 PMNS angles
- [[../op-zeta-mass-spectrum/Phase3_results.md]] — ζ.1 NO masses + Σm_ν
- [[../op-iota-charge-pmns-unification/Phase3_results.md]] — ι.1 4-sector chirality
- [[../op-kappa-mixing-numerator/Phase3_results.md]] — κ.1 B²_up=13/4
- [[../../INDEX.md]], [[../../PREDICTIONS_REGISTRY.md]]
