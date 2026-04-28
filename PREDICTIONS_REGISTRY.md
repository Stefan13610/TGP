---
title: "TGP Predictions Registry"
date: 2026-04-28
status: live
type: index
tags:
  - TGP
  - predictions
  - registry
  - falsification
  - flasks
related:
  - "[[research/op-phase3-uv-completion/Phase3_R_final_results.md]]"
  - "[[research/op-phase2-quantum-gravity/Phase2_R_final_results.md]]"
  - "[[research/closure_2026-04-26/KNOWN_ISSUES.md]]"
---

# TGP Predictions Registry

**Single source of truth** for every falsifiable prediction made by the
Theoria Geometriae Pulvis (TGP) program, with numerical anchor,
falsification target, time horizon, and the **Zenodo "prediction flask"**
that pre-registered the prediction with an immutable DOI timestamp.

**Master ledger at registry creation:** 281 cumulative structural
verifications (M9 13 + M10 42 + M11 62 + Phase 1 50 + Phase 2 54 + Phase 3 60).

## How to read this registry

Each row carries:
- **Anchor**: the TGP value with uncertainty / sympy-exact label.
- **Reference**: the experimental/theoretical comparison value.
- **Target**: which experiment or measurement would falsify the entry.
- **Horizon**: the realistic window for that target.
- **Status**: see taxonomy below.
- **Flask**: Zenodo DOI of the deposit that pre-registered the prediction
  (immutable timestamp), or `master-only` if it lives only in the
  workshop tree.
- **Master**: the Phase results file in `TGP/TGP_v1/research/` that
  documents the most recent check.

### Status taxonomy

| Status | Meaning |
|--------|---------|
| **TESTED-PASS** | Already compared to existing data; agreement at the stated drift. |
| **LIVE** | Falsifiable by an experiment that is currently running or scheduled. |
| **LOCKED** | Sympy-exact / algebraic / multi-constraint forced; falsified only by direct contradiction. |
| **STRUCTURAL** | Internal-consistency lock; falsifiable only via downstream theoretical contradiction. |
| **STRUCTURAL-POSTULATE** | Pointer placed in a derivation frame; first-principles derivation pending (long-term). |

### Flask DOIs (timestamps for prediction priority)

| Flask | DOI | Sectors covered |
|-------|-----|-----------------|
| **tgp-core** | [10.5281/zenodo.19670324](https://doi.org/10.5281/zenodo.19670324) | Gravity, GW, BH, dark energy, foundations |
| **tgp-leptons** | [10.5281/zenodo.19706861](https://doi.org/10.5281/zenodo.19706861) | Lepton sector, mass ratios, Koide |
| **tgp-qm** | [10.5281/zenodo.19712596](https://doi.org/10.5281/zenodo.19712596) | QM foundations |
| **tgp-sc** v2 | [10.5281/zenodo.19670557](https://doi.org/10.5281/zenodo.19670557) | Superconductivity closure |

---

## Sector 1 — Gravity, WEP, and equivalence-principle tests

| # | Anchor | Reference | Target / Horizon | Status | Flask | Master |
|---|---|---|---|---|---|---|
| **G1** | **WEP η_TGP = 3.54·10⁻¹⁷** (n=2 forced) | MICROSCOPE-1 bound: η < 1.1·10⁻¹⁵ | MICROSCOPE-2 sensitivity ~10⁻¹⁷, 2027–2028 | **LIVE** | tgp-core | `Phase2_E_results.md` (E.2), `closure_2026-04-26/KNOWN_ISSUES.md` (B.2) |
| **G2** | **n=2 smoothness exponent** unique under (C¹ + WEP) | C⁰ exponent classes ruled out | If MICROSCOPE-2 sees η > 10⁻¹⁶ → n≠2 → TGP scope falsified | **LOCKED** | tgp-core | Phase2.E.2 |
| **G3** | **Margin to MICROSCOPE-1 bound: 3.7·10¹⁶×** | observation ≥ TGP×10¹⁶ | window remains until MICROSCOPE-2 + STEP / next-gen | **TESTED-PASS** | tgp-core | Phase2.E.2 |

---

## Sector 2 — Gravitational waves: speed, polarization, dispersion

| # | Anchor | Reference | Target / Horizon | Status | Flask | Master |
|---|---|---|---|---|---|---|
| **GW1** | **c_T = c_s exactly** (no birefringence) | GW170817: \|c_T−c\|/c < 9·10⁻²² | LIGO O5 / Cosmic Explorer high-frequency phase delay | **TESTED-PASS** | tgp-core | M9.3.5, Phase2.A.5 |
| **GW2** | **3 GW polarization DOF** (h_+, h_×, h_b = h_L) | GR: 2 DOF | LIGO O5 / Einstein Telescope: search for h_b/h_L scalar mode | **LIVE** | tgp-core | M9.3.4, Phase1.F.4 |
| **GW3** | **h_b = h_L degeneracy** (single-Φ axiom) | GR: undefined longitudinal | ET / pulsar-timing scalar polarization tests | **STRUCTURAL** | tgp-core | M9.3.4 |
| **GW4** | **m_σ²/m_s² = 2** sympy-exact integer | OPE-invariant, preserved 4/4 UV | LISA / PTA low-k phase shift ~2.9% signature, 2030–2035 | **LOCKED** | tgp-core | `closure_2026-04-26/` (Path B σ_ab 11/11), Phase3.F |
| **GW5** | **No vector modes** (single-scalar) | GR: 0; bimetric: 2 vector | LIGO O5 polarization fit: detection of h_v ≠ 0 falsifies | **STRUCTURAL** | tgp-core | M9.3.4 |
| **GW6** | **m_s ≪ ω_LIGO** (effective masslessness in LIGO band) | M_eff/ω_LIGO ≈ 7·10⁹ | LIGO O5: c_GW(f) frequency-dependence above precision floor | **STRUCTURAL** | tgp-core | Phase1.F |

---

## Sector 3 — Black-hole shadows / photon rings (OP-M92)

| # | Anchor | Reference | Target / Horizon | Status | Flask | Master |
|---|---|---|---|---|---|---|
| **BH1** | **r_ph^TGP / r_ph^GR = 1.293 ± 0.003%** | GR baseline 1.000 | ngEHT M87* + Sgr A* photon-ring imaging, 2030–2032 | **LIVE** | tgp-core | `op-m92/`, Phase1.B.2 |
| **BH2** | **Δb_crit = +14.56%** vs GR | GR: 0% | ngEHT photon-ring resolution; selection among candidates A/B/D | **LIVE** | tgp-core | `op-m92/` |
| **BH3** | **Universality across multi-BH** (M87*, Sgr A*, future targets) | GR: M-scaling universal | ngEHT 2030+: cross-source consistency at ≥3 sources | **STRUCTURAL** | tgp-core | `op-m92/` |

---

## Sector 4 — Dark energy and cosmological constant

| # | Anchor | Reference | Target / Horizon | Status | Flask | Master |
|---|---|---|---|---|---|---|
| **DE1** | **w = −1.000 exact** (no phantom crossing) | Planck+BAO: w ∈ [−1.03, −0.96] (1σ) | DESI DR2/DR3 CPL fit: phantom w < −1 at 3σ falsifies | **LOCKED** | tgp-core | M10.1.4, Phase1.C |
| **DE2** | **w_a = 0** (no time-evolution) | ΛCDM w_a = 0 | DESI DR3: \|w_a\| > 0.05 (3σ) falsifies | **LOCKED** | tgp-core | M10.1.4 |
| **DE3** | **T-Λ ratio = 1.0203** (TGP/observed vacuum energy) | Weinberg anthropic window [0.5, 2.0] | DESI / Euclid precision: ratio > 2.0 falsifies | **TESTED-PASS** | tgp-core | Phase3.F.4, Phase2.A |
| **DE4** | **Friedmann ratio 0.9808** (3Ω_Λ/(8π) ≈ 0.0817 vs 1/12 = 0.0833) | M9.1″ Path-B vacuum prefactor | Euclid/DESI: deviation from 0.9808 beyond 5% gate falsifies | **LOCKED** | tgp-core | Phase3.F.6 |
| **DE5** | **Λ_E = γ/12** sympy-exact (γ/6 × 1/2) | Path-B vacuum prefactor | First-principles γ identification — long-term | **LOCKED** | tgp-core | Phase3.E.5 (PARTIAL DERIVED) |

---

## Sector 5 — Hubble / S₈ / neutrino tensions

| # | Anchor | Reference | Target / Horizon | Status | Flask | Master |
|---|---|---|---|---|---|---|
| **C1** | **H₀ from T-Λ closure** (Φ_eq ≡ H₀ identification) | Cepheid 73 vs CMB 67 km/s/Mpc | DESI DR3 + local-distance ladder reconciliation | **STRUCTURAL** (B.4 STRENGTHENED) | tgp-core | Phase3.E.4, M10.4 |
| **C2** | **TGP does NOT reduce S₈ tension** (no MOND-like enhancement) | σ₈^obs ≈ 0.811 | DESI DR3 / Euclid weak lensing: TGP ≠ ΛCDM signature in S₈ falsifies | **STRUCTURAL** | tgp-core | `s8_tension/` |
| **C3** | **Σm_ν ∈ structural bracket** (TGP-compat range) | CMB+BAO: Σm_ν < 0.12 eV | DESI + CMB-S4 neutrino constraints | **LIVE** | tgp-core | `neutrino_msw/` |

---

## Sector 6 — Particle / lepton sector

| # | Anchor | Reference | Target / Horizon | Status | Flask | Master |
|---|---|---|---|---|---|---|
| **L1** | **r_21 = (A_μ/A_e)^4 = 206.77** | PDG 206.768 (drift 1·10⁻⁵ = 0.001%) | precision PDG updates | **TESTED-PASS** | tgp-leptons | `tgp-leptons-paper/paper/tgp_leptons.tex` |
| **L2** | **r_31 = (A_τ/A_e)^4 = 3477** | PDG 3477.23 (drift 7·10⁻⁵) | precision m_τ updates | **TESTED-PASS** | tgp-leptons | tgp-leptons |
| **L3** | **k_eff = 4.0008** (k=4 integer selection) | k=4 integer (drift 2·10⁻⁴) | dimension-uniqueness: k≠4 in d=3 falsifies | **LOCKED** | tgp-leptons | Thm. `thm:k4` (tgp-leptons) |
| **L4** | **K_koide = 2/3** sympy-exact | PDG 0.66666 (two indep. extractions of g₀^τ agree to 0.036%) | precision m_τ + Koide test | **LOCKED** | tgp-leptons | Thm. `thm:koide-theorem` |
| **L5** | **N = 3 generations** (metric-singularity barrier at γ_crit = 2.206) | observed 3 light leptons | 4th-generation discovery falsifies | **LOCKED** | tgp-leptons | Thm. `thm:N3` |

---

## Sector 7 — UV completion (4-of-4 structural compatibility)

| # | Anchor | Reference | Target / Horizon | Status | Flask | Master |
|---|---|---|---|---|---|---|
| **UV1** | **Litim invariant g\*·λ\* = 0.1349** | Reuter NGFP: 0.135 (drift 0.07%) | AS NGFP refutation by FRG / lattice — long-term | **STRUCTURAL** | tgp-core | Phase3.A, R.F.1 |
| **UV2** | **String/KKLT compat** (T-Λ ∈ [0.5, 2.0]) | Weinberg anthropic window | string vacuum landscape ruling out window — long-term | **STRUCTURAL** | tgp-core | Phase3.B |
| **UV3** | **LQG γ_Imm ≈ 0.2375** + area gate 61.4 dex | Phase 2.D.5 60.93 dex IR pointer | LQG dynamics anomaly cancellation — long-term | **STRUCTURAL** | tgp-core | Phase3.C |
| **UV4** | **CDT d_s = 4 − 2η = 3.948** (η_Phase1.D = 0.026) | CDT 4.02 ± 0.10 (distance 0.72σ < 3σ gate) | CDT continuum-limit existence proof — long-term | **STRUCTURAL** | tgp-core | Phase3.D |
| **UV5** | **4/4 UV synthesis matrix** (AS / string / LQG / CDT all compat) | independent UV frameworks | empirical UV signature preferring subset of 4 | **STRUCTURAL** | tgp-core | Phase3.F (CAPSTONE) |
| **UV6** | **Spectral dim flow d_s: 4 → 2** (AS + CDT independently) | dimensional reduction signature | non-perturbative RG divergence between AS and CDT falsifies the convergence | **STRUCTURAL** | tgp-core | Phase3.A.6, Phase3.D |
| **UV7** | **Δ_target = 0.114** in heat-kernel a₂ frame | a₂ ⊃ (1/2)V''²; ξ_geom=1, α(α−1)=2 | full a₂-from-first-principles closure — long-term | **STRUCTURAL-POSTULATE** | tgp-core | Phase3.E.6 |

---

## Sector 8 — Foundational locks (parameters that gate everything else)

| # | Anchor | Reference | Target / Horizon | Status | Flask | Master |
|---|---|---|---|---|---|---|
| **F1** | **Single-Φ axiom** (one propagating scalar, no auxiliaries) | TGP foundation | discovery of additional light scalar coupling to gravity | **LOCKED** | tgp-core | TGP_FOUNDATIONS §1 |
| **F2** | **K(φ) = K_geo·φ⁴** (α = 2 selection within class (C1)–(C3)) | Phase 1 Theorem `thm:alpha2` | observable inconsistent with K∝φ⁴ falsifies GL-substrate ansatz | **LOCKED** | tgp-core | Phase1.A.1, Remark `rem:alpha2-closure` |
| **F3** | **β = γ vacuum** (V'(Φ_eq=1) = 0 exact, Φ_eq = 1) | stability + scale-locking | Φ_eq drift / instability falsifies | **LOCKED** | tgp-core | Phase1.F.3 |
| **F4** | **α₀ = 0.114 / 0.168² ≈ 4.04472** (sympy exact rational 1069833/264500) | algebraic structural anchor | direct measurement deviating from rational value | **LOCKED** | tgp-core | Phase2.B.3, Phase1.B.3 |
| **F5** | **g̃ ≈ 0.9803** (Phase 2 EFT survival, drift 0.0306%) | T-Λ entropy scaling | full quantum-gravity entropy refutation | **STRUCTURAL** | tgp-core | Phase2.E.3, Phase3.F |
| **F6** | **κ = √(32πG_N) ≈ 10.0265** (graviton coupling, FP-quantized) | Phase 2.A.1 KEYSTONE | graviton-loop / scattering κ-renormalization deviation | **STRUCTURAL** | tgp-core | Phase2.A.1 |
| **F7** | **14 founding constraints, zero-drift** (preserved across all phases) | Phase 1 founding set | any drift > 0% in audited parameters falsifies | **LOCKED** | tgp-core | every Phase R-final |

---

## Sector 9 — QM foundations and superconductivity (orthogonal flasks)

| # | Anchor | Reference | Target / Horizon | Status | Flask | Master |
|---|---|---|---|---|---|---|
| **QM1** | **Born rule from substrate** (n=2 ψ²-rule emergent, not postulated) | standard QM postulate | alternative quantum formalism with n≠2 normalization | **STRUCTURAL** | tgp-qm | tgp-qm-paper |
| **QM2** | **Decoherence via Φ-substrate interaction** | standard QM has no first-principles collapse | spontaneous-localization / objective-collapse experiments | **STRUCTURAL** | tgp-qm | tgp-qm-paper |
| **SC1** | **L6–L10 refinements** (v2 deposit) — see flask paper for full list | precision SC closure conditions | material-by-material verification | **STRUCTURAL** | tgp-sc v2 | tgp-sc-paper |

---

## Cross-sector falsification roadmap (by horizon)

| Window | Experiment | Key TGP entries at risk |
|--------|-----------|----------------------------|
| **2027–2028** | MICROSCOPE-2 | G1 (η = 3.54·10⁻¹⁷), G2 (n=2) |
| **2027+** | LIGO O5 | GW2 (3 DOF), GW5 (no vector), GW6 (dispersion bound) |
| **2027+** | DESI DR2/DR3 | DE1 (w=−1), DE2 (w_a=0), DE3 (T-Λ), C1 (H₀), C2 (S₈), C3 (Σm_ν) |
| **2028+** | Euclid | DE3, DE4 (Friedmann ratio), C2 |
| **2030–2032** | ngEHT | BH1 (r_ph ratio 1.293), BH2 (Δb_crit +14.56%), BH3 (multi-BH) |
| **~2035** | LISA / pulsar-timing arrays | GW4 (m_σ²/m_s² = 2 → 2.9% low-k phase shift) |
| **long-term** | full QG / UV experiment | UV1–UV7, F5, F6 (UV-route selection) |

---

## Audit policy

- **Every entry has a flask DOI or `master-only` tag.** A claim without a DOI is not yet pre-registered for falsification priority.
- **Every entry has a master cross-reference.** The flask paper is the immutable timestamp; the master `Phase*_results.md` is the live working state.
- **TESTED-PASS rows** carry a current drift figure. **LIVE** rows carry a horizon. **LOCKED** rows carry a sympy/algebraic provenance. **STRUCTURAL** rows carry an explicit "what would falsify".
- **No entry is ever silently retracted.** A falsified entry transitions to **FALSIFIED** with a dated note, and the affected derivations are tracked in `closure_2026-04-26/KNOWN_ISSUES.md`.

## Maintenance

This file is the index. Detailed updates live in:
- `research/op-phase3-uv-completion/Phase3_R_final_results.md` (post-Phase-3 ledger)
- `research/closure_2026-04-26/KNOWN_ISSUES.md` (live audit log)
- Per-flask `POST_PHASE3_*` notes in each flask repo's `research/` folder.

When a new prediction enters: add a row, cite the master file that introduced it, and (if depositing) link the new flask DOI.

When an experimental result lands: update the row's status (LIVE → TESTED-PASS / FALSIFIED) with a dated note, and add a one-line entry in the flask's research-folder changelog.
