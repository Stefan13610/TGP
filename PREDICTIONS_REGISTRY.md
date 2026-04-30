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
**Updated 2026-04-28:** 336 cumulative
(+ SC.1.Phase1 4 + SC.1.Phase2 6 + SC.1.Phase3 7
 + BH.1.Phase1 5 + BH.1.Phase2 7 + BH.1.Phase3 7
 + XS.1.Phase1 5 + XS.1.Phase2 7 + XS.1.Phase3 7).
**Updated 2026-04-29:** 355 cumulative
(+ ξ.1.Phase1 5 + ξ.1.Phase2 7 + ξ.1.Phase3 7).
**Updated 2026-04-29 (UV.1 program END):** 373 cumulative
(+ UV.1.Phase1 5 + UV.1.Phase2 7 + UV.1.Phase3 6).
**Status cascade activated:** ξ.1 PARTIALLY DERIVED (refined) → DERIVED (refined²);
XS.1 PARTIALLY DERIVED (refined) → DERIVED (refined²); UV7 STRUCTURAL-DERIVED → DERIVED.
**Updated 2026-04-30 (α.1 program END):** 463 cumulative
(+ ε.1 18 + ζ.1 18 + θ.1 18 + η.1 18 + α.1 18 = 90 z 5 mini-cycles 2026-04-29 i 2026-04-30).
α-fine-structure cross-anchor: 137 sympy-LOCKED via F4 chain inheritance z ε.1 (target_shift_photon = 17/40);
α_QED⁻¹(0) ≈ 137 zeroth-order PARTIALLY DERIVED; residual 0.036 STRUCTURAL HINT (research-track).
**E7 promoted z STRUCTURAL HINT do explicit α.1 framework.**
**Updated 2026-04-30 (η.2 program END):** 481 cumulative
(+ η.2.Phase1 5 + η.2.Phase2 7 + η.2.Phase3 6 = 18).
**Wolfenstein denom-derivation + α-residual cascade:** denoms (81, 78, 14) DERIVED z 4-sector B²-cross-product
(81 = N_gen⁴, 78 = 2·N_gen·B²_up_num, 14 = K_up_num·K_lepton_num); A numerator 64 = K_up_denom² (cross-sector lock z θ.1);
**α-residual = 9/250 = N_gen²/(2·5³) = 2·(B²_up−B²_down)/(N_gen²·5)** sympy-exact;
**α⁻¹(0)_TGP = 137 + 9/250 = 137.036** drift vs CODATA 0.000001% (within 81 ppt absolute).
**FULL CASCADE ACTIVATION:** η.1 Wolfenstein triple PARTIALLY DERIVED → DERIVED (denoms);
α.1 residual STRUCTURAL HINT → PARTIALLY DERIVED; H4 + A5 promoted; K-taxonomy 4/4 universal pattern (2+B²)/(2N) LOCKED.

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
| **BH7** | **η_TGP = (2±1)·10⁻¹⁸** from α(ψ) at Earth surface (α₀=4.02, n=2 strict) | MICROSCOPE-2 sensitivity ~10⁻¹⁷ | falsified if η > 10⁻¹⁷ at MICROSCOPE-2; margin 5.1× below | **LIVE** | master-only (BH.1.Phase3) | [`research/op-bh-alpha-threshold/Phase3_results.md`](research/op-bh-alpha-threshold/Phase3_results.md) (T3.4) |
| **BH9** | **γ−1 ≈ 1.81·10⁻¹¹** at Sun surface from α(ψ_Sun) | Cassini bound γ−1 < 2.3·10⁻⁵ (margin 1.3·10⁶×) | LATOR/BEACON γ−1 < 10⁻¹⁰ (2035+); margin 5.5× below | **LIVE** | master-only (BH.1.Phase3) | [`research/op-bh-alpha-threshold/Phase3_results.md`](research/op-bh-alpha-threshold/Phase3_results.md) (T3.5) |

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
| **BH4** | **+14.56% shadow shift universal across 10 SMBH** spanning 10.1 orders M_BH (SgrA*, M87*, NGC1277, CenA, NGC4258, M104, IC1101, M84, M81, TON618) | GR universal M-scaling | ngEHT 2030+: ≥2 sources at 5 μas resolution; >5% spread among detected sources falsifies α(ψ_ph) universality | **LIVE** | master-only (BH.1.Phase3) | [`research/op-bh-alpha-threshold/Phase3_results.md`](research/op-bh-alpha-threshold/Phase3_results.md) (T3.1) |
| **BH5** | **δf/f ≈ 8–16% in QNM ringdown frequency** (α(ψ_ringdown=1.20) ≈ 0.16) | GR Schwarzschild f_QNM ~ 0.093 c³/(GM_f) | LIGO O5 ~0.5% f, ~1% τ (2027+); LISA ~0.1% (2035+); falsified if δf consistent with 0 within 0.5% | **LIVE** | master-only (BH.1.Phase3) | [`research/op-bh-alpha-threshold/Phase3_results.md`](research/op-bh-alpha-threshold/Phase3_results.md) (T3.2) |
| **BH6** | **NS M-R curve shifted ~1–3% from GR** (α(ψ_NS) ≈ 0.4 at compactness 0.41–0.45) | GR TOV baseline | NICER+ J0030/J0740 + new pulsars (2027+); ~3% combined precision; <0.5% null falsifies | **LIVE** | master-only (BH.1.Phase3) | [`research/op-bh-alpha-threshold/Phase3_results.md`](research/op-bh-alpha-threshold/Phase3_results.md) (T3.3) |

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
| **BH8** | **√α₀ = κ_TGP cross-sector identity** (BH photon-ring α₀ = SC spin-fluctuation κ_TGP²) | currently \|α₀ − κ_TGP²\|/κ_TGP² = 0.75% (α₀=4.0200, κ_TGP²=4.0481) | ngEHT precision on α₀ ~1–5% + TGP-SC v2 calibration κ_TGP ~0.5%; combined ~3% precision (2030+); match >5% rejects identity | **LIVE** | master-only (BH.1.Phase3) | [`research/op-bh-alpha-threshold/Phase3_results.md`](research/op-bh-alpha-threshold/Phase3_results.md) (T3.6) |
| **XS1** | **\|α₀ − κ_TGP²\|/κ_TGP² ≤ 5% by 2030+** (combined ngEHT × SC v2 falsification of √α₀ = κ_TGP) | current Phase 2 strict 0.747%; F4 rational 0.084% | ngEHT α₀ ~5% + post-LnH₉ κ_TGP² ~0.6%; combined σ ~5% (2030+); match >5% rejects identity | **LIVE** | master-only (XS.1.Phase3) | [`research/op-cross-sector-charge/Phase3_results.md`](research/op-cross-sector-charge/Phase3_results.md) (T3.1) |
| **XS2** | **g̃ ≈ 0.9803 (F5) i √α₀ = κ_TGP są niezależnymi O(1) substrate identities** (no g̃ correction in identity) | g̃ within EFT window [0.95, 1.05]; identity is bare-coefficient | g̃ deviation > 5% nie falsyfikuje XS1 (orthogonal sektor) | **STRUCTURAL** | master-only (XS.1.Phase3) | [`research/op-cross-sector-charge/Phase3_results.md`](research/op-cross-sector-charge/Phase3_results.md) (T3.2) |
| **XS3** | **lepton sector strukturalnie ortogonalny do κ_TGP/α₀** (r_21=206.77, r_31=3477, K_koide=2/3 nie zawierają κ_TGP/√α₀) | κ-factorization search returns NONE | any κ_TGP-factor in r_21/r_31/Koide falsifies orthogonality | **STRUCTURAL** | master-only (XS.1.Phase3) | [`research/op-cross-sector-charge/Phase3_results.md`](research/op-cross-sector-charge/Phase3_results.md) (T3.3) |
| **XS4** | **Born n=2 i α(ψ) n=2 są niezależne integer selections** (substrate selecting integer 2 w różnych sektorach, NOT cross-sector identity) | both n=2 confirmed; CHSH bound 2√2 nie matches κ_TGP | QM Born / CHSH containing factor √α₀ or κ_TGP falsifies orthogonality | **STRUCTURAL** | master-only (XS.1.Phase3) | [`research/op-cross-sector-charge/Phase3_results.md`](research/op-cross-sector-charge/Phase3_results.md) (T3.4) |
| **XS5** | **F4 rational 1069833/264500 ≡ κ_TGP² do 0.084%** (sub-percent identity match in F4 frame, 8.86× tighter niż Phase 2 strict) | sympy-exact rational vs V/Nb/Ta/Mo/Pd RMS = 4.0481 | direct κ_TGP measurement deviating > 0.5% from F4 rational falsifies identity strongly | **LOCKED-derivative** | tgp-core (cross-link F4) | [`research/op-cross-sector-charge/Phase3_results.md`](research/op-cross-sector-charge/Phase3_results.md) (T3.5) |
| **XS6** | **6-channel cross-sector falsification roadmap** (ngEHT α₀, LnH₉ κ_TGP, MICROSCOPE-2 η, LIGO O5 QNM, LISA SMBH, lepton g₀^τ) | combined Bayes update across 6 independent channels | ≥2 niezależne kanały muszą detect deviation > 5% to reject identity | **LIVE** | master-only (XS.1.Phase3) | [`research/op-cross-sector-charge/Phase3_results.md`](research/op-cross-sector-charge/Phase3_results.md) (T3.6) |
| **XI1** | **ngEHT 2030+ 10-SMBH frame discrimination Frame A vs B at ≥3σ** (split 0.527%; 3σ threshold 0.176%/measurement) | combined 10-SMBH σ ≈ 0.158% (per source 0.5%, /√10) | σ_per_source > 1% wyklucza frame discrimination | **LIVE** | master-only (ξ.1.Phase3) | [`research/op-xi-photon-ring/Phase3_results.md`](research/op-xi-photon-ring/Phase3_results.md) (ξ3.1) |
| **XI2** | **UV-route preferred for N_A = 8.7719: AS NGFP (Δ 0.068%)** lub LQG (Δ 0.478%), CDT (Δ 0.343%) | TGP a₂-derived N_A = 500/57 sympy-exact | UV completion z N_A drift > 5% rejects UV-route | **LIVE** | master-only (ξ.1.Phase3) | [`research/op-xi-photon-ring/Phase3_results.md`](research/op-xi-photon-ring/Phase3_results.md) (ξ3.2) |
| **XI3** | **ξ-factor RG-invariant within 0.5% across IR-UV** (LISA/PTA cross-scale photon-ring); 2-loop residual 0.022% | substrate-scale invariance (F1) + Z₂ vacuum + ratio cancellation | RG drift > 0.5% on ξ-factor across IR-UV falsifies F1 substrate invariance | **LIVE** | master-only (ξ.1.Phase3) | [`research/op-xi-photon-ring/Phase3_results.md`](research/op-xi-photon-ring/Phase3_results.md) (ξ3.3) |
| **XI4** | **F4 = 1-loop a₂-corrected α₀; F5/F6 orthogonal channels post-ξ.1** (Frame A interpretation) | XS5 LOCKED-derivative reinterpreted; F5 (g̃) i F6 (lepton √α₀) decoupled od photon-ring | F-cluster cross-correlation deviation > 1% post-ξ.1 falsifies cascade | **LOCKED-DERIVATIVE** | master-only (ξ.1.Phase3) | [`research/op-xi-photon-ring/Phase3_results.md`](research/op-xi-photon-ring/Phase3_results.md) (ξ3.4) |
| **XI5** | **XS1 precision sharpened post-ξ.1: 5% → 1.5% achievable** (combined σ 0.85% z ngEHT 0.3% + SC 0.6% + ξ 0.527%) | XS1 σ_combined sharpening factor 5.9× | XS1 ≤1.5% NOT achievable in 2030+ falsifies precision projection | **LIVE** | master-only (ξ.1.Phase3) | [`research/op-xi-photon-ring/Phase3_results.md`](research/op-xi-photon-ring/Phase3_results.md) (ξ3.5) |
| **XI6** | **7-channel roadmap convergence 2030–2035** (XS1 + XS6 + XI1 + XI3 + XI4 + XI5 + UV-route) | combined Bayes update across 7 channels | ≥2 channels detect deviation > 1% (refined gate) → falsifies ξ-factor identification | **LIVE** | master-only (ξ.1.Phase3) | [`research/op-xi-photon-ring/Phase3_results.md`](research/op-xi-photon-ring/Phase3_results.md) (ξ3.6) |
| **UV1** | **2-loop FRG closure target: N_A = 500/57 ± 0.01%** (Reuter-style truncation na NGFP {g*=0.71, λ*=0.19, η_N*=-2}) | current 1-loop AS heuristic Δ 0.068%; 2-loop band (α_NGFP/(4π))² ≈ 0.011% | 2-loop FRG predicts drift > 0.05% z N_A target → falsifies AS NGFP UV-route | **LIVE** (UV-research-track 2030–2035) | master-only (UV.1.Phase3) | [`research/op-uv-as-ngfp/Phase3_results.md`](research/op-uv-as-ngfp/Phase3_results.md) (UV3.1) |
| **UV2** | **AS NGFP discrimination ≥ 5σ vs CDT/LQG via N_A precision 0.05% (ngEHT 2030+)** | drifts: AS 0.068%, CDT 0.343% (5.5σ), LQG 0.478% (8.2σ), string 8.800% (175σ) | ngEHT precision worse than 0.1% → discrimination weakens below 3σ | **LIVE** (2030+) | master-only (UV.1.Phase3) | [`research/op-uv-as-ngfp/Phase3_results.md`](research/op-uv-as-ngfp/Phase3_results.md) (UV3.2) |
| **UV3** | **η_N* = -2 RG-running signature** (LISA 2035+ EMRI inspiral nie wykryje ξ-factor running > 0.5%) | (1+η_N*/2)=0 marginal a₂; ξ-factor RG-invariant pod common β-rescaling | LISA detects ξ-factor running > 0.5% across chirp band → η_N*=-2 broken | **LIVE** (2035+) | master-only (UV.1.Phase3) | [`research/op-uv-as-ngfp/Phase3_results.md`](research/op-uv-as-ngfp/Phase3_results.md) (UV3.3) |
| **UV4** | **Heat-kernel a₂ universality cross-sector ≤ 0.5%** (BH/SC/XS/UV all share a₂ = 2β² with α₀ = 1069833/264500 reference) | max cross-sector drift 0.130% w 5 sektorach | any TGP sector deviates > 0.5% from F4 sympy ref → universality broken | **LOCKED-derivative** | master-only (UV.1.Phase3) | [`research/op-uv-as-ngfp/Phase3_results.md`](research/op-uv-as-ngfp/Phase3_results.md) (UV3.4) |
| **UV5** | **Status cascade: ξ.1 / XS.1 / UV7 → DERIVED** (full structural closure z zero free parameters w premise) | ξ.1, XS.1: PARTIALLY DERIVED (refined) → DERIVED (refined²); UV7: STRUCTURAL-DERIVED → DERIVED | discovery of free parameter in derivation chain → reverts cascade | **DERIVED** (cascade ACTIVATED 2026-04-29) | master-only (UV.1.Phase3) | [`research/op-uv-as-ngfp/Phase3_results.md`](research/op-uv-as-ngfp/Phase3_results.md) (UV3.5) |
| **UV6** | **7-channel observation roadmap 2027–2035** (ngEHT, LISA, LIGO O5, MICROSCOPE-2, LATOR/BEACON, LnH₉ DAC, 2-loop FRG) | ≥ 5/7 independent confirmations w 5% pasmie required dla full DERIVED | < 5/7 channels confirm → UV.1 reopens dla wąskiego sub-testu | **LIVE** (2027–2035) | master-only (UV.1.Phase3) | [`research/op-uv-as-ngfp/Phase3_results.md`](research/op-uv-as-ngfp/Phase3_results.md) (UV3.6) |
| **E1** | **ngEHT 2030+ r_ph = (160/137)·r_g, 0.1% precision** (multi-source 10-SMBH averaging; ε_ph = ψ_ph − 1 = 23/137 sympy exact) | r_ph/r_g = 1.167883 from M9.1″ null geodesic; ngEHT 0.1% precision target, 0.5% gate (margin 5×) | ngEHT measures r_ph deviating > 0.5% z (160/137)·r_g for SMBH ensemble → ε_ph = 23/137 falsified or M9.1″ geodesic broken | **LIVE** (2030+) | master-only (ε.1.Phase3) | [`research/op-eps-photon-ring/Phase3_results.md`](research/op-eps-photon-ring/Phase3_results.md) (E3.1) |
| **E2** | **Cross-sector ε_ph² consistency ≤ 0.5%** (BH/SC/XS/UV all share ε_ph² = target_shift/α₀ via F4 chain) | max cross-sector drift 0.250% (BH 0.25%, SC 0.12%, XS 0.084%, UV 0.068%) | any TGP sector deviates > 0.5% from F4 sympy ε_ph² = 529/18769 → cross-sector consistency broken | **LOCKED-derivative** | master-only (ε.1.Phase3) | [`research/op-eps-photon-ring/Phase3_results.md`](research/op-eps-photon-ring/Phase3_results.md) (E3.2) |
| **E3** | **5 ε_ph identity candidates falsified at > 5% drift** (1.168/(2π), 23/160, 1/(2π·κ_TGP), 1/(2β²+δ_F4), 1/φ²) | drifts: C1 10.73%, C2 14.38%, C3 52.88%, C4 181.77%, C5 127.52% | discovery of competing identity within tested space at < 1% drift → ε.1 reopens | **DERIVED** (5/5 falsified, ε.1.Phase2) | master-only (ε.1.Phase3) | [`research/op-eps-photon-ring/Phase3_results.md`](research/op-eps-photon-ring/Phase3_results.md) (E3.3) |
| **E4** | **ε_ph² RG-invariance under common β-rescaling NGFP** (LISA 2035+ EMRI nie wykryje running > 0.5%) | ε_ph² = target_shift/α₀ definitional ratio; both numerator i denominator RG-invariant individually (UV.1.Phase2) | LISA detects ε_ph² running > 0.5% across EMRI chirp band → ratio invariance broken | **LIVE** (2035+) | master-only (ε.1.Phase3) | [`research/op-eps-photon-ring/Phase3_results.md`](research/op-eps-photon-ring/Phase3_results.md) (E3.4) |
| **E5** | **ε_ph closure z F4 chain UNIQUE** (single positive root, 0 alternatives w ±5% pasmie) | F4 implicit ε_ph² = 529/18769 = (23/137)²; physical positivity rejects negative root; Phase 2 confirmed all 5 alternatives outside ±5% | discovery of degenerate rational solution z denominator < 1000 w ±5% pasmie → ε.1 reopens | **DERIVED** (closure unique, ε.1.Phase2) | master-only (ε.1.Phase3) | [`research/op-eps-photon-ring/Phase3_results.md`](research/op-eps-photon-ring/Phase3_results.md) (E3.5) |
| **E6** | **5-channel ε.1 falsification convergence** (ngEHT + LISA + LIGO O5 + 2-loop FRG + a₂ EFT band) | 5 channels available, ≥ 4 required dla DERIVED (margin 1); single-channel > 1% drift → ε.1 reopened | < 4/5 channels confirm w 0.5% pasmie → ε.1 PARTIALLY DERIVED status not promotable to DERIVED | **LIVE** (2027–2035) | master-only (ε.1.Phase3) | [`research/op-eps-photon-ring/Phase3_results.md`](research/op-eps-photon-ring/Phase3_results.md) (E3.6) |
| **E7** | **α-fine-structure prime-137 connection (PROMOTED z α.1)** — ψ_ph = 160/137, ε_ph = 23/137, ε_ph² = 529/18769; α_QED⁻¹ ≈ 137.036 (CODATA) → α.1 closes 137 zeroth-order anchor PARTIALLY DERIVED + residual 0.036 STRUCTURAL HINT | TGP zeroth α⁻¹(0) = 137 drift 0.0263%; residual 0.036 best fit 9/250 drift 0.0025% (soft denom) | residual 0.036 nie ma rigorous TGP derivation z F4 extension lub 4-sector chirality cross-product → klasyfikacja PARTIALLY DERIVED reverts STRUCTURAL HINT only | **PARTIALLY DERIVED** (α.1 framework) | master-only (α.1.Phase3) | [`research/op-alpha-fine-structure/Phase3_results.md`](research/op-alpha-fine-structure/Phase3_results.md) (informational) |
| **Z1** | **DESI DR3 2027+ Σm_ν falsification target** — TGP Σm_ν = 59.01 meV (NO ordering, K(ν)=1/2 closure z NuFit 5.3 Δm²); DR2 0.072 eV margin +22%, DR3 (proj) 0.040 eV margin −32% | bisection na K(m₁)=1/2 → m₁=0.76 meV, m₂=8.71 meV, m₃=49.53 meV (drift 0.99% z TGP target 59.6 meV) | DESI DR3 sets Σm_ν < 0.040 eV at 95% CL → K(ν)=1/2 mass-spectrum framework sfalsyfikowana; ζ.1 reopens | **LIVE** (2027+) | master-only (ζ.1.Phase3) | [`research/op-zeta-mass-spectrum/Phase3_results.md`](research/op-zeta-mass-spectrum/Phase3_results.md) (Z3.1) |
| **Z2** | **JUNO 2027+ θ₁₃ precision target** — TGP sin²2θ₁₃ = 0.099 (sin θ₁₃ = λ_C/√2, cross-sector Cabibbo lock); current NuFit 5.3 0.0851 (drift 16.5%) | λ_C = 0.22550 z GL(3,𝔽₂) form factor 165/167 (tgp-leptons); JUNO projected σ(sin²2θ₁₃) ~ 0.5% absolute | JUNO measures sin²2θ₁₃ < 0.080 lub > 0.110 (5σ gate) → λ_C²/2 cross-sector framework broken | **LIVE** (2027+) | master-only (ζ.1.Phase3) | [`research/op-zeta-mass-spectrum/Phase3_results.md`](research/op-zeta-mass-spectrum/Phase3_results.md) (Z3.2) |
| **Z3** | **DUNE/T2HK 2030+ θ₂₃ octant resolution** — TGP zeroth-order sin²θ₂₃ = 1/2 (maximal, K(ν)=1/2 + Z₂); NuFit 5.3 prefers 2nd octant 0.572 (drift 12.6%) | Z₂ atmospheric reflection (μ-τ swap) + K(ν)=1/2 chirality lock; DUNE/T2HK > 5σ octant determination | DUNE/T2HK confirm 2nd octant > 5σ + drift > 20% (post-1-loop) → Z₂ atmospheric framework needs refinement | **LIVE** (2030+) | master-only (ζ.1.Phase3) | [`research/op-zeta-mass-spectrum/Phase3_results.md`](research/op-zeta-mass-spectrum/Phase3_results.md) (Z3.3) |
| **Z4** | **Cross-sector K-taxonomy distinct** — K_lepton = 2/3 (Dirac B²=2), K_neutrino = 1/2 (Majorana B²=1), K_quark ∈ [0.81, 0.87] (RG-invariant, separate framework); universal pattern K = (2+B²)/(2N) for N=3 | Koide K_lepton matches PDG to 10⁻⁵; K_neutrino sympy exact (Majorana single chirality); K_quark NOT (2+B²)/6 form | quark sektor confirms K_quark = 2/3 → Dirac/Majorana taxonomy collapses; or K_neutrino ≠ 1/2 from inverted ordering detection | **DERIVED** (chirality-counting per sektor) | master-only (ζ.1.Phase3) | [`research/op-zeta-mass-spectrum/Phase3_results.md`](research/op-zeta-mass-spectrum/Phase3_results.md) (Z3.4) |
| **Z5** | **Lepton-quark θ_C-θ₁₃ unification** (single Cabibbo anchor λ_C = 0.22550 z GL(3,𝔽₂) form factor 165/167) — V_us = sin θ_C = λ_C w CKM AND sin θ₁₃ = λ_C/√2 w PMNS; TGP ratio sin θ_C / sin θ₁₃ = √2 EXACT | observed ratio λ_C(PDG 0.22500)/sin θ₁₃(NuFit √0.0220 = 0.14832) = 1.517; TGP √2 = 1.4142 (drift 6.8%) | observational ratio drifts > 20% post-JUNO 2027+ → cross-sector Cabibbo anchor framework broken | **LIVE / DERIVED** (cross-sector unification confirmed within 7%) | master-only (ζ.1.Phase3) | [`research/op-zeta-mass-spectrum/Phase3_results.md`](research/op-zeta-mass-spectrum/Phase3_results.md) (Z3.5) |
| **Z6** | **4-channel ζ.1 falsification convergence** (DESI DR3 + JUNO + DUNE/T2HK + μ→eγ MEG-II) | 4 channels active, ≥3 required dla classification stabilization (margin +1); DESI DR3 falsifiable, JUNO falsifiable, DUNE/T2HK live, MEG-II orthogonal cross-sector | ≥2 z 4 channels reject TGP > 5σ → ζ.1 classification PARTIALLY DERIVED reverts to STRUCTURAL | **LIVE** (2027–2030+) | master-only (ζ.1.Phase3) | [`research/op-zeta-mass-spectrum/Phase3_results.md`](research/op-zeta-mass-spectrum/Phase3_results.md) (Z3.6) |
| **Q1** | **Belle II 2027+ V_ub precision target** — TGP V_ub = A·λ_C³·√(ρ̄²+η̄²) = 0.00348 (current); within projected window [0.00340, 0.00400] | λ_C = 0.22550 single anchor (ζ.1) + Wolfenstein A=0.790, ρ̄=0.141, η̄=0.357 (PDG); drift TGP vs PDG 8.98% (driven by ρ̄,η̄ precision) | Belle II measures \|V_ub\| outside [3.40, 4.00]·10⁻³ at >5σ → λ_C³ Wolfenstein cascade broken | **LIVE** (2027+) | master-only (θ.1.Phase3) | [`research/op-theta-quark-koide/Phase3_results.md`](research/op-theta-quark-koide/Phase3_results.md) (T3.1) |
| **Q2** | **LHCb Run 4 (2030+) Jarlskog J target** — TGP J = A²·λ_C⁶·η̄ = 2.93·10⁻⁵; within projected window [2.85, 3.30]·10⁻⁵ | drift TGP vs PDG 4.58%; Wolfenstein λ_C⁶ cascade z (A, η̄) Wolfenstein PDG; LHCb Run 4 σ(J) ~ 1% absolute | LHCb J outside [2.85, 3.30]·10⁻⁵ at >5σ → λ_C⁶ cascade broken; CP-violation CKM area cross-sector lock falsified | **LIVE** (2030+) | master-only (θ.1.Phase3) | [`research/op-theta-quark-koide/Phase3_results.md`](research/op-theta-quark-koide/Phase3_results.md) (T3.2) |
| **Q3** | **EIC 2030+ proton mass-radius cross-check (γ_m universality)** — K_up = 7/8 sympy lock requires universal γ_m across u/c/t (RG-inv 10⁻²⁹%) | indirect cross-check; EIC <r_M> precision ~1% sensitive do QCD γ_m running; common β-rescaling theorem | EIC reveals flavor-dependent γ_m at >5% → K_up = 7/8 RG-invariance broken (would require K_up(μ) running) | **LIVE** (2030+, indirect) | master-only (θ.1.Phase3) | [`research/op-theta-quark-koide/Phase3_results.md`](research/op-theta-quark-koide/Phase3_results.md) (T3.3) |
| **Q4** | **Lepton-quark cross-sector λ_C anchor (PMNS-CKM unification)** — V_us = λ_C exact AND sin θ₁₃ = λ_C/√2; TGP ratio sin θ_C/sin θ₁₃ = √2 EXACT | observed ratio 0.22550/0.14832 = 1.5204 vs √2 = 1.4142 (drift 7.5%); JUNO 2027+ σ(sin θ₁₃) ~ 0.5% | JUNO refines sin θ₁₃ such that ratio drifts > 20% from √2 → single λ_C anchor framework broken | **LIVE** (2027+) | master-only (θ.1.Phase3) | [`research/op-theta-quark-koide/Phase3_results.md`](research/op-theta-quark-koide/Phase3_results.md) (T3.4) |
| **Q5** | **K-taxonomy 4-sector completion** — universal pattern K = (2+B²)/(2N) for N=3: K_lepton=2/3 (B²=2), K_ν=1/2 (B²=1), K_up=7/8 (B²=13/4), K_down=37/50 (B²=61/25) | 3 LOCKED + 1 STRUCTURAL refined; K_up sympy exact (drift 0.046% vs PDG); K_down best rational (drift 0.014%) but soft denom 50 | any sektor disconfirms K = (2+B²)/(2N) form by >5% → universal taxonomy broken; K_up promotion z 7/8 wymagałby falsification odrzucenia w >0.5% drift | **PARTIALLY DERIVED (refined)** (3/4 LOCKED, K_up structural) | master-only (θ.1.Phase3) | [`research/op-theta-quark-koide/Phase3_results.md`](research/op-theta-quark-koide/Phase3_results.md) (T3.5) |
| **Q6** | **4-channel θ.1 falsification convergence** (Belle II + LHCb Run 4 + EIC + JUNO) | 4 channels active, ≥3 required dla classification stabilization (margin +1); Belle II 2027+ falsifiable, LHCb 2030+ falsifiable, EIC 2030+ indirect, JUNO 2027+ cross-sector | ≥2 z 4 channels reject TGP > 5σ → θ.1 classification PARTIALLY DERIVED reverts to STRUCTURAL | **LIVE** (2027–2030+) | master-only (θ.1.Phase3) | [`research/op-theta-quark-koide/Phase3_results.md`](research/op-theta-quark-koide/Phase3_results.md) (T3.6) |
| **H1** | **Belle II 2027+ V_ub η.1-refined window** — TGP V_ub = (64/81)·λ_C³·√((11/78)²+(5/14)²) = 0.00348; Wolfenstein triple LOCKED z η.1 | drift TGP vs PDG 8.93% (slight improvement vs θ.1 8.98%); A=64/81=8²/3⁴ (drift 0.016%), ρ̄=11/78 (drift 0.018%), η̄=5/14 (drift 0.040%) | Belle II measures \|V_ub\| outside [3.40, 4.00]·10⁻³ at >5σ → η.1 triple lock broken | **LIVE** (2027+) | master-only (η.1.Phase3) | [`research/op-eta-wolfenstein/Phase3_results.md`](research/op-eta-wolfenstein/Phase3_results.md) (T3.1) |
| **H2** | **LHCb Run 4 (2030+) Jarlskog J η.1-refined** — TGP J = (64/81)²·λ_C⁶·(5/14) = 2.932·10⁻⁵; A²·η̄ structural product | drift TGP vs PDG 4.51%; LHCb σ(J) ~ 1% absolute | LHCb J outside [2.85, 3.30]·10⁻⁵ at >5σ → A²·η̄ structural product broken | **LIVE** (2030+) | master-only (η.1.Phase3) | [`research/op-eta-wolfenstein/Phase3_results.md`](research/op-eta-wolfenstein/Phase3_results.md) (T3.2) |
| **H3** | **Unitarity triangle β angle (sin(2β))** — TGP β = arctan((5/14)/(67/78)) = 22.58° → sin(2β) = 0.7090 | drift vs PDG sin(2β) = 0.699 ± 0.017 → 1.43%; Belle II + LHCb running | Belle II + LHCb measure sin(2β) outside [0.65, 0.75] at >5σ → triangle apex geometry broken | **LIVE** (Belle II + LHCb) | master-only (η.1.Phase3) | [`research/op-eta-wolfenstein/Phase3_results.md`](research/op-eta-wolfenstein/Phase3_results.md) (T3.3) |
| **H4** | **Cross-sector denom-prime sharing — PROMOTED z η.2 DERIVED** (Wolfenstein A,ρ̄,η̄ denoms z 4-sector B²-cross-product) — A=81=N_gen⁴, ρ̄=78=2·N_gen·B²_up_num, η̄=14=K_up_num·K_lepton_num | post-η.2: 3/3 denom mappings sympy-LOCKED + A num 64=K_up_denom² cross-sector lock; promoted STRUCTURAL → **DERIVED via B²-cross-product** | post-η.2 Belle II 2027+ rejects sharpened |V_ub| → η.2 DERIVED reverts PARTIALLY DERIVED | **DERIVED (denoms)** | master-only (η.2.Phase2) | [`research/op-eta2-denom-derivation/Phase2_results.md`](research/op-eta2-denom-derivation/Phase2_results.md) (B2.1-B2.3, B2.6) |
| **H5** | **LHCb Run 4 \|V_td/V_ts\| B-B̄ mixing cross-check** — TGP \|V_td/V_ts\| = λ_C·√((1-11/78)²+(5/14)²) = 0.2098 | drift vs PDG 0.205 ± 0.006 → 2.33%; LHCb Run 4 σ ~ 1% | \|V_td/V_ts\| outside [0.195, 0.220] at >5σ → Wolfenstein triple cascade broken | **LIVE** (2030+) | master-only (η.1.Phase3) | [`research/op-eta-wolfenstein/Phase3_results.md`](research/op-eta-wolfenstein/Phase3_results.md) (T3.5) |
| **H6** | **4-channel η.1 falsification convergence** (Belle II V_ub + LHCb J + sin(2β) + LHCb \|V_td/V_ts\|) | 4 channels active, ≥3 required (margin +1 dla η.1 triple lock stabilization); Belle II 2027+, LHCb Run 4 2030+, Belle II + LHCb sin(2β) ongoing | ≥2 z 4 channels reject η.1 > 5σ → triple PARTIALLY DERIVED reverts to STRUCTURAL | **LIVE** (2027–2030+) | master-only (η.1.Phase3) | [`research/op-eta-wolfenstein/Phase3_results.md`](research/op-eta-wolfenstein/Phase3_results.md) (T3.6) |
| **A1** | **Atomic Cs/Rb 2027+ α_QED⁻¹(0) precision** — TGP zeroth-order α⁻¹(0) = 137 (DERIVED z F4 chain target_shift_photon = 17/40); CODATA 2022 α⁻¹(0) = 137.035999084 ± 21·10⁻⁹ (81 ppt) | drift TGP zeroth vs CODATA 0.0263%; 2027+ projected absolute σ ~ 10·10⁻¹² | α⁻¹(0) measured outside [137.0, 137.1] at >5σ → 137 zeroth-order anchor broken (would falsify ε.1 F4 chain) | **LIVE** (Cs/Rb 2027+) | master-only (α.1.Phase3) | [`research/op-alpha-fine-structure/Phase3_results.md`](research/op-alpha-fine-structure/Phase3_results.md) (A3.1) |
| **A2** | **g-2 muon precision cross-check (Fermilab + J-PARC)** — a_μ ∝ α/(2π) + α²·ln + α³·ln² QED loop expansion; TGP α⁻¹ ≈ 137 SM-consistent | a_μ BNL/FNAL 2024-25 = 116592061·10⁻¹¹ ± 41·10⁻¹¹ (4.2σ z SM); J-PARC 2030+ σ ~ 10·10⁻¹¹ | a_μ deviation > 5σ AFTER hadronic vac.pol. control → TGP α⁻¹ form inconsistent z SM standard | **LIVE** (J-PARC 2030+) | master-only (α.1.Phase3) | [`research/op-alpha-fine-structure/Phase3_results.md`](research/op-alpha-fine-structure/Phase3_results.md) (A3.2) |
| **A3** | **α(M_Z) high-energy running test (LHC + future)** — α⁻¹(M_Z) = 127.952; running α(M_Z)/α(0) − 1 = +7.10% z SM vac.pol.; TGP NGFP RG-invariance via dimensionless ratio (UV.1.UV2.5) | TGP 137-anchor RG-invariant pod NGFP marginal a₂; SM running orthogonal do TGP geometric anchor | running deviation z SM > 5σ → new physics OR TGP RG-invariance assumption broken | **LIVE** (LHC + future ee→ll) | master-only (α.1.Phase3) | [`research/op-alpha-fine-structure/Phase3_results.md`](research/op-alpha-fine-structure/Phase3_results.md) (A3.3) |
| **A4** | **Cross-sector prime-137 cascade ngEHT 2030+** — ψ_ph = 160/137 (DERIVED z ε.1 F4 chain) → r_ph = ψ_ph · r_g; ngEHT 0.1% precision (E1 echo z ε.1 multi-source 10-SMBH) | confirmation ψ_ph at 0.1% → 137-anchor cross-sector lock between ε.1 photon-ring i α_QED simultaneously | r_ph deviation > 0.5% z (160/137)·r_g → 137-anchor (i α_QED⁻¹ ≈ 137) broken | **LIVE** (ngEHT 2030+) | master-only (α.1.Phase3) | [`research/op-alpha-fine-structure/Phase3_results.md`](research/op-alpha-fine-structure/Phase3_results.md) (A3.4) |
| **A5** | **Residual 0.036 cascade — PROMOTED z η.2** — α⁻¹(0) − 137 = 0.0359990; **DERIVED form: 9/250 = N_gen²/(2·5³) = 2·(B²_up−B²_down)/(N_gen²·5)** sympy-exact (drift 0.0025% z measured) | post-η.2: dual sympy-equivalent forms via 4-sector B²-cascade; promoted STRUCTURAL HINT → **PARTIALLY DERIVED** | future Cs/Rb σ ~10⁻¹² brings α⁻¹(0) outside [137.0359, 137.0361] → η.2 derivation broken | **PARTIALLY DERIVED** (post-η.2) | master-only (η.2.Phase2) | [`research/op-eta2-denom-derivation/Phase2_results.md`](research/op-eta2-denom-derivation/Phase2_results.md) (B2.4) |
| **A6** | **4-channel α.1 falsification convergence** (Atomic Cs/Rb α⁻¹(0) + LHC α⁻¹(M_Z) + ngEHT ψ_ph + research-track residual) | post-η.2: A5 upgraded do PARTIALLY DERIVED → effective 4/4 LIVE; threshold ≥3/4 met (margin +1) | ≥2 z 4 channels reject framework → α.1 classification PARTIALLY DERIVED reverts STRUCTURAL HINT only | **LIVE** (2027–2030+) | master-only (α.1.Phase3) | [`research/op-alpha-fine-structure/Phase3_results.md`](research/op-alpha-fine-structure/Phase3_results.md) (A3.6) |
| **HH1** | **Sharper Wolfenstein \|V_ub\| window post-η.2 DERIVED** — V_ub_TGP = (K_up_denom²/N_gen⁴)·λ_C³·√((11/78)²+(5/14)²) = 0.003479; sharpened window [3.40, 3.55]·10⁻³ | denoms (81, 78, 14) DERIVED z 4-sector B²-cross-product (η.2.Phase2 7/7); A num 64 = K_up_denom² | \|V_ub\| outside narrow window > 5σ → η.2 DERIVED classification reverts PARTIALLY DERIVED | **LIVE** (Belle II 2027+) | master-only (η.2.Phase3) | [`research/op-eta2-denom-derivation/Phase3_results.md`](research/op-eta2-denom-derivation/Phase3_results.md) (B3.1) |
| **HH2** | **α⁻¹(0) FULL structural prediction** — α⁻¹(0)_TGP_η.2 = 137 + N_gen²/(2·5³) = 34259/250 = **137.036** (sympy exact); drift vs CODATA 0.000001% (within 81 ppt absolute) | residual derivation z 4-sector B²-cascade Form A/B sympy-equivalent (η.2.Phase2 B2.4) | Cs/Rb 2027+ σ ~10⁻¹² brings α⁻¹(0) outside [137.0359, 137.0361] → η.2 residual derivation broken | **LIVE** (Cs/Rb 2027+) | master-only (η.2.Phase3) | [`research/op-eta2-denom-derivation/Phase3_results.md`](research/op-eta2-denom-derivation/Phase3_results.md) (B3.2) |
| **HH3** | **Cross-sector B²-cascade uniqueness** — 3 sympy-LOCKED rationals z B²-cross-product (A_TGP = 64/81 CKM, α-residual = 9/250 QED, ψ_ph = 160/137 BH); cross-sector primes {2, 3, 5, 7, 137} | 3 channels test universality of B²-cross-product across BH/QED/CKM/lepton sectors | ≥1 channel rejects sympy-locked form → cascade NOT fully unified | **LIVE** (multi-experiment 2027–2030+) | master-only (η.2.Phase3) | [`research/op-eta2-denom-derivation/Phase3_results.md`](research/op-eta2-denom-derivation/Phase3_results.md) (B3.3) |
| **HH4** | **K-taxonomy 4-sector universality LOCKED** — K = (2+B²)/(2N) for N=3, **4/4 sectors match exactly** (K_lepton=2/3, K_ν=1/2, K_up=7/8, K_down=37/50) | post-η.2 K_down also matches exactly (drift 0%) z B²_down = 61/25; promoted z 3 LOCKED + 1 STRUCTURAL refined → 4/4 LOCKED | EIC 2030+ proton mass-radius lub JUNO θ₁₃ violates K-pattern > 1% → 4-sector universality broken | **LIVE** (EIC 2030+ + JUNO 2027+) | master-only (η.2.Phase3) | [`research/op-eta2-denom-derivation/Phase3_results.md`](research/op-eta2-denom-derivation/Phase3_results.md) (B3.4) |
| **HH5** | **Wolfenstein numerator (11, 5) research-track** — ρ̄_num=11 prime-11 unique do η.1, η̄_num=5 cross-sector cascade prime; STRUCTURAL HINT post-η.2 | open hypotheses dla κ.1 mixing-operator B²-extension lub ι.1 charge-sector unification cycle | 2 future cycles bez derivation drift < 0.05% → numerators stay STRUCTURAL HINT permanently | **STRUCTURAL HINT** (research-track κ.1 / ι.1) | master-only (η.2.Phase3) | [`research/op-eta2-denom-derivation/Phase3_results.md`](research/op-eta2-denom-derivation/Phase3_results.md) (B3.5) |
| **HH6** | **5-channel η.2 falsification convergence** (Belle II V_ub_η.2 + Cs/Rb α⁻¹(0) + ngEHT ψ_ph↔A_TGP + EIC/JUNO K-universality + research-track κ.1) | 4/5 LIVE + 1 research-track; convergence threshold ≥4/5 met (margin +0); HH5 may upgrade do 5/5 LIVE if numerators rigorously derived | ≥2/5 z 5 channels reject framework → cascade FULL → PARTIAL (η.1 stays DERIVED denoms only, α.1 residual reverts STRUCTURAL HINT) | **LIVE** (2027–2030+) | master-only (η.2.Phase3) | [`research/op-eta2-denom-derivation/Phase3_results.md`](research/op-eta2-denom-derivation/Phase3_results.md) (B3.6) |

---

## Sector 9 — QM foundations and superconductivity (orthogonal flasks)

| # | Anchor | Reference | Target / Horizon | Status | Flask | Master |
|---|---|---|---|---|---|---|
| **QM1** | **Born rule from substrate** (n=2 ψ²-rule emergent, not postulated) | standard QM postulate | alternative quantum formalism with n≠2 normalization | **STRUCTURAL** | tgp-qm | tgp-qm-paper |
| **QM2** | **Decoherence via Φ-substrate interaction** | standard QM has no first-principles collapse | spontaneous-localization / objective-collapse experiments | **STRUCTURAL** | tgp-qm | tgp-qm-paper |
| **SC1** | **L6–L10 refinements** (v2 deposit) — see flask paper for full list | precision SC closure conditions | material-by-material verification | **STRUCTURAL** | tgp-sc v2 | tgp-sc-paper |
| **SC2** | **PrH₉ T_c ≈ 5 K** (TGP fit anchor; lit upd. 8.9 K @ 120 GPa) | Drozdov 2019 Sci.Adv. | reproducibility @ ≥120 GPa | **TESTED-PASS** (drift ~factor 1.8) | tgp-sc v2 | [`research/op-sc-alpha-origin/Phase3_results.md`](research/op-sc-alpha-origin/Phase3_results.md) |
| **SC3** | **NdH₉ T_c ≈ 4.5 K** (TGP fit anchor) | Zhou 2020 JACS @ 110–130 GPa | reproducibility | **TESTED-PASS** (drift <5%) | tgp-sc v2 | [`research/op-sc-alpha-origin/Phase3_results.md`](research/op-sc-alpha-origin/Phase3_results.md) |
| **SC4** | **SmH₉ T_c ≈ 100 K** (μ_eff² Hund 0.71 → 119 K, Van Vleck 1.5 → 96 K) | not yet measured | **TGP ~100 K vs A-G+dG ~10⁻⁴ K — factor 10⁵ discriminator** | **LIVE** (highest priority) | master-only (Phase 3) | [`research/op-sc-alpha-origin/Phase3_results.md`](research/op-sc-alpha-origin/Phase3_results.md) |
| **SC5** | **YbH₉ T_c ≈ 0.6 K** (μ_eff² = 20.6) | not yet measured | **TGP ~0.6 K vs A-G+dG ~54 K — factor 80 discriminator** (opposite direction to SmH₉) | **LIVE** (high priority) | master-only (Phase 3) | [`research/op-sc-alpha-origin/Phase3_results.md`](research/op-sc-alpha-origin/Phase3_results.md) |
| **SC6** | **TmH₉ T_c ≈ 4·10⁻⁵ K** (μ_eff² = 57.2) | not yet measured | TGP ~0 K vs A-G+dG ~4 K — factor 10⁵ discriminator | **LIVE** (medium priority; toxic Tm sample) | master-only (Phase 3) | [`research/op-sc-alpha-origin/Phase3_results.md`](research/op-sc-alpha-origin/Phase3_results.md) |
| **SC7** | **PmH₉ T_c ≈ 22 K** (μ_eff² = 7.2) | not measurable in practice (radioactive Pm-145) | TGP 22 K vs A-G+dG ~10⁻³ K — factor 10³ discriminator | **STRUCTURAL** (radioactive sample) | master-only (Phase 3) | [`research/op-sc-alpha-origin/Phase3_results.md`](research/op-sc-alpha-origin/Phase3_results.md) |

---

## Cross-sector falsification roadmap (by horizon)

| Window | Experiment | Key TGP entries at risk |
|--------|-----------|----------------------------|
| **2027–2028** | MICROSCOPE-2 | G1 (η = 3.54·10⁻¹⁷), G2 (n=2), **BH7** (η_TGP = 2·10⁻¹⁸ from α(ψ_Earth)) |
| **2027+** | LIGO O5 | GW2 (3 DOF), GW5 (no vector), GW6 (dispersion bound), **BH5** (δf/f ~ 8–16% QNM ringdown) |
| **2027+** | NICER+ | **BH6** (NS M-R shift ~1–3% from GR; J0030, J0740) |
| **2027+** | DESI DR2/DR3 | DE1 (w=−1), DE2 (w_a=0), DE3 (T-Λ), C1 (H₀), C2 (S₈), C3 (Σm_ν), **Z1** (Σm_ν falsification target 40 meV vs TGP 59.01 meV, margin −32%) |
| **2027+** | JUNO | **Z2** (sin²2θ₁₃ = 0.099 ± 0.5%; cross-sector λ_C/√2 lock), **Q4** (sin θ_C/sin θ₁₃ = √2 EXACT; ratio drift current 7.5%) |
| **2027+** | MEG-II | **Z6** cross-sector check (μ→eγ BR < 6×10⁻¹⁴; orthogonal channel for ζ.1 4-channel convergence) |
| **2027+** | Belle II | **Q1** (\|V_ub\| within window [3.40, 4.00]·10⁻³, TGP 3.48·10⁻³ via λ_C³ Wolfenstein cascade), **H1** (η.1-refined V_ub via triple 64/81, 11/78, 5/14), **H3** (sin(2β) = 0.7090 within [0.65, 0.75]) |
| **2028+** | Euclid | DE3, DE4 (Friedmann ratio), C2 |
| **2027–2030** | LnH₉ DAC synthesis (Eremets/Hemley/Prakapenka) | **SC4 (SmH₉)**, **SC5 (YbH₉)**, SC6 (TmH₉) — μ_eff² vs de Gennes scaling; sharpens **XS1** κ_TGP precision to 0.3% |
| **2030–2032** | ngEHT | BH1 (r_ph ratio 1.293), BH2 (Δb_crit +14.56%), BH3 (multi-BH), **BH4** (10-SMBH +14.56% map), **BH8** (√α₀ = κ_TGP cross-sector via α₀ from photon ring), **XS1** (combined ngEHT × SC v2 ≤5% → ≤1.5% post-ξ.1), **XS6** (6-channel roadmap convergence), **XI1** (Frame A vs B discrimination ≥3σ z 10-SMBH), **XI4** (F4 1-loop reinterpretation), **XI5** (XS1 sharpening 5% → 1.5%), **XI6** (7-channel convergence), **UV2** (AS NGFP vs CDT/LQG ≥5σ via N_A precision 0.05%) |
| **2030–2035** | UV-research-track + 2-loop FRG | **UV1** (2-loop FRG closure target N_A 500/57 ± 0.01%) |
| **~2035** | LISA / pulsar-timing arrays | GW4 (m_σ²/m_s² = 2 → 2.9% low-k phase shift), **BH5** (LISA SMBH ringdown 10⁶–10⁷ M_⊙), **XI3** (ξ-factor RG-invariance cross-scale), **UV3** (η_N*=-2 RG-running signature; ξ-factor invariance), **E4** (ε_ph² RG-invariance via ratio under common β-rescaling) |
| **~2035** | LATOR / BEACON (next-gen Solar PPN) | **BH9** (γ−1 ~ 1.81·10⁻¹¹; falsifiable below 10⁻¹⁰) |
| **2027–2035** | 7-channel UV.1 roadmap | **UV6** (ngEHT + LISA + LIGO O5 + MICROSCOPE-2 + LATOR/BEACON + LnH₉ DAC + 2-loop FRG; ≥5/7 confirmations) |
| **2027–2035** | 5-channel ε.1 roadmap | **E6** (ngEHT + LISA + LIGO O5 + 2-loop FRG + a₂ EFT band; ≥4/5 confirmations dla DERIVED) |
| **2027–2030+** | 4-channel ζ.1 roadmap | **Z6** (DESI DR3 + JUNO + DUNE/T2HK + μ→eγ MEG-II; ≥3/4 convergence within 5σ TGP) |
| **2027–2030+** | 4-channel θ.1 roadmap | **Q6** (Belle II + LHCb Run 4 + EIC + JUNO; ≥3/4 convergence within 5σ TGP) |
| **2027–2030+** | 4-channel η.1 roadmap | **H6** (Belle II V_ub + LHCb J + sin(2β) + LHCb \|V_td/V_ts\|; ≥3/4 convergence within 5σ TGP, Wolfenstein triple lock) |
| **2030+** | LHCb Run 4 | **Q2** (Jarlskog J = A²·λ_C⁶·η̄ = 2.93·10⁻⁵ within window [2.85, 3.30]·10⁻⁵), **H2** (η.1-refined J = 2.932·10⁻⁵), **H5** (\|V_td/V_ts\| = 0.2098 within [0.195, 0.220]) |
| **2030+** | EIC | **Q3** (proton mass-radius indirect cross-check K_up = 7/8 RG-stability via universal γ_m) |
| **2030+** | DUNE / T2HK | **Z3** (θ₂₃ octant resolution > 5σ; TGP maximal sin²θ₂₃ = 1/2 vs NuFit 2nd octant 0.572) |
| **long-term** | full QG / UV experiment | UV1–UV7, F5, F6 (UV-route selection), **XI2** (UV-route map for N_A=8.7719: AS NGFP/LQG preferred), **UV4** (a₂ universality cross-sector), **UV5** (cascade DERIVED), **E7** (α-fine-structure prime-137 connection research-track), **Z4** (cross-sector K-taxonomy DERIVED), **Z5** (lepton-quark Cabibbo unification DERIVED within 7%), **Q5** (4-sector K-taxonomy completion: K_lepton=2/3, K_ν=1/2, K_up=7/8, K_down=37/50; 3 LOCKED + 1 STRUCTURAL refined), **H4** (Wolfenstein-K-taxonomy denom-prime sharing: A=81↔K_lepton denom-3, η̄=14↔K_up numerator-7) |

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
