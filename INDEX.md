---
title: "TGP_v1 — Index hub"
date: 2026-04-28
type: hub
tags:
  - TGP
  - index
  - navigation
related:
  - "[[README.md]]"
  - "[[PREDICTIONS_REGISTRY.md]]"
  - "[[research/op-phase3-uv-completion/Phase3_R_final_results.md]]"
  - "[[research/op-bh-alpha-threshold/Phase3_results.md]]"
  - "[[research/op-sc-alpha-origin/Phase3_results.md]]"
  - "[[research/op-cross-sector-charge/Phase1_results.md]]"
  - "[[research/op-cross-sector-charge/Phase2_results.md]]"
  - "[[research/op-cross-sector-charge/Phase3_results.md]]"
  - "[[research/op-xi-photon-ring/Phase1_results.md]]"
  - "[[research/op-xi-photon-ring/Phase2_results.md]]"
  - "[[research/op-xi-photon-ring/Phase3_results.md]]"
---

# TGP_v1 — Index hub

A one-screen navigation map for the TGP research ecosystem
post-Phase-3 (closure 2026-04-28). For the project overview, see
[`README.md`](README.md). This file is the **navigation entry point**.

## At a glance

- **Master verification ledger:** **355** cumulative closures
  (M9 13 + M10 42 + M11 62 + Phase 1 50 + Phase 2 54 + Phase 3 60 + SC.1.Phase1 4 + SC.1.Phase2 6 + SC.1.Phase3 7 + BH.1.Phase1 5 + BH.1.Phase2 7 + BH.1.Phase3 7 + XS.1.Phase1 5 + XS.1.Phase2 7 + XS.1.Phase3 7 + ξ.1.Phase1 5 + ξ.1.Phase2 7 + ξ.1.Phase3 7;
  closure_2026-04-26's 35 tests are tracked in the phase-ledger table but
  overlap conceptually with Phase 1/2/3 follow-ups so are not double-counted here).
- **Active program:** none — **ξ.1 program END** 2026-04-29 (3 phases, 19 sub-tests, 5/5 + 6/7 + 7/7 PASS); ξ-factor identified jako a₂ EFT 1-loop correction (0.527% F4 vs strict split); 6 new predictions XI1–XI6 registered; XS.1 promoted PARTIALLY DERIVED → PARTIALLY DERIVED (refined); UV7 promoted STRUCTURAL-POSTULATE → STRUCTURAL-DERIVED; AS NGFP najbliższy UV-route match dla N_A = 8.7719 (Δ 0.068%). Promotion do full DERIVED czeka na UV completion w long-term track (op-uv-renormalizability-research/UV1).
- **Recently closed cycles:**
  - **SC.1 program (3 phases) END:** α_PB scaling map registered, 5 new SC predictions (SC2–SC7) in registry.
  - **BH.1 program (3 phases) END:** multi-source ISSUE strukturalnie zamknięty; ψ_th=1 i n=2 promoted to DERIVED (Z₂ + WEP-MICROSCOPE-2); α₀ ≈ 4.02 PARTIALLY DERIVED z cross-sector STRUCTURAL HINT α₀ = κ_TGP² (match 0.75% z TGP-SC); 6 new BH predictions (BH4–BH9) registered with horizons 2027–2035.
  - **XS.1 program (3 phases) END:** cross-sector identity √α₀ = κ_TGP promoted z STRUCTURAL HINT do **PARTIALLY DERIVED** via substrate-action S[Φ, g, J, ψ_e] under TGP single-Φ + Z₂ + K_geo=1 + φ_eq=1; M_BH = M_SC = M_universal ≈ 4.03–4.05; F4 form match 0.084%, Phase 2 strict 0.747%; 6 new XS predictions (XS1–XS6) registered with **6-channel falsification roadmap** 2027–2035 (ngEHT + LnH₉ DAC + MICROSCOPE-2 + LIGO O5 + LISA + lepton precision).
  - **ξ.1 program (3 phases) END:** ξ-factor (0.527% F4 vs strict split) **identified strukturalnie** jako a₂ EFT 1-loop heat-kernel correction; Frame A (F4 0.114) = 1-loop a₂-corrected, Frame B (11/97) = tree-level bare geometric; XS.1 promoted PARTIALLY DERIVED → PARTIALLY DERIVED **(refined)**; UV7 STRUCTURAL-POSTULATE → STRUCTURAL-DERIVED; AS NGFP najbliższy UV-route match dla N_A = 8.7719 (Δ 0.068%); 6 new predictions XI1–XI6 (frame discrimination ngEHT 2030+, UV-route map, RG-invariance LISA/PTA, F4 reinterpretation, XS1 sharpening 5%→1.5%, 7-channel convergence).
- **4 Zenodo flask papers** deposited as immutable timestamped predictions.
- **1 predictions registry** ([`PREDICTIONS_REGISTRY.md`](PREDICTIONS_REGISTRY.md))
  cross-linking every flask to falsification target + horizon.
- **Active queued work:** none — research now in long-term track (op-uv-renormalizability-research/ UV1–UV7 + F5/F6 supplement). Experimental wait on ngEHT 2030+, LIGO O5 2027+, LISA 2035+, MICROSCOPE-2 2030+, NICER+ 2027+, LATOR/BEACON 2035+, SmH₉/YbH₉ DAC synthesis 2027–2030.

## Top-level entry points

| Document | Purpose |
|----------|---------|
| [`README.md`](README.md) | Public project overview, abstract, 40-prediction highlights, repository structure, build instructions |
| [`PREDICTIONS_REGISTRY.md`](PREDICTIONS_REGISTRY.md) | Single source of truth for every falsifiable prediction (sector, anchor, target, horizon, status, flask DOI) |
| [`TGP_FOUNDATIONS.md`](TGP_FOUNDATIONS.md) | Axiomatic content (4 axioms, single-Φ, Z₂, substrate) |
| [`DEPENDENCIES.md`](DEPENDENCIES.md) / [`DEPENDENCIES_REVERSE.md`](DEPENDENCIES_REVERSE.md) | Logical-dependency graph forward/backward across modules |

## The four "prediction-flask" deposits

Each flask is an immutable Zenodo deposit with a DOI timestamp; the
master state may evolve, but the deposit is the timestamped prediction.

| Flask | DOI | Sector | Status | Audit note |
|-------|-----|--------|--------|------------|
| **tgp-core** | [10.5281/zenodo.19670324](https://doi.org/10.5281/zenodo.19670324) | Axioms / gravity / GW / BH / dark energy | clean (POST_PHASE3 addendum, 6 strengthenings) | [tgp-core-paper/research/POST_PHASE3_ADDENDUM_2026-04-28.md](https://github.com/Stefan13610/tgp-core-paper) |
| **tgp-leptons** | [10.5281/zenodo.19706861](https://doi.org/10.5281/zenodo.19706861) | Charged leptons, Koide, Cabibbo | micro-fix (r_21: 206.74 → 206.77, drift 1·10⁻⁵) | [tgp-leptons-paper/research/POST_PHASE3_NOTE_2026-04-28.md](https://github.com/Stefan13610/tgp-leptons-paper) |
| **tgp-qm** | [10.5281/zenodo.19712596](https://doi.org/10.5281/zenodo.19712596) | Emergent QM (Born, CHSH, spin, decoherence) | clean (orthogonal to Phase 1/2/3) | [tgp-qm-paper/research/POST_PHASE3_NOTE_2026-04-28.md](https://github.com/Stefan13610/tgp-qm-paper) |
| **tgp-sc** v2 | [10.5281/zenodo.19670557](https://doi.org/10.5281/zenodo.19670557) | Superconductivity T_c closure (5 families) | clean (v2 already current with L6–L10) | [tgp-sc-paper/research/POST_PHASE3_NOTE_2026-04-28.md](https://github.com/Stefan13610/tgp-sc-paper) |

## Phase ledger (341 cumulative)

| Block | Tests | Status | Master file |
|-------|------:|--------|-------------|
| **M9** classical gravity | 13 | CLOSED | `core/sek*` |
| **M10** FRW cosmology | 42 | CLOSED | `core/sek05*` |
| **M11** quantum closure (9 sub-cycles + R-final) | 62 | CLOSED | `research/op-quantum-closure/` |
| **closure_2026-04-26** (T-FP 12/12 + T-Λ 7/7 + T-α 5/5 + Path B σ_ab 11/11) | 35 | CLOSED | [`research/closure_2026-04-26/`](research/closure_2026-04-26/) |
| **Phase 1** covariant 4D | 50 | CLOSED | [`research/op-phase1-covariant/`](research/op-phase1-covariant/) |
| **Phase 2** quantum gravity / EFT (Donoghue-grade) | 54 | CLOSED | [`research/op-phase2-quantum-gravity/Phase2_R_final_results.md`](research/op-phase2-quantum-gravity/Phase2_R_final_results.md) |
| **Phase 3** UV-completion structural-consistency audit | 60 | CLOSED 2026-04-28 | [`research/op-phase3-uv-completion/Phase3_R_final_results.md`](research/op-phase3-uv-completion/Phase3_R_final_results.md) |
| **SC.1.Phase1** α_PB ↔ α_0 unit-bridge audit (negative) | 4 | CLOSED 2026-04-28 | [`research/op-sc-alpha-origin/Phase1_results.md`](research/op-sc-alpha-origin/Phase1_results.md) |
| **SC.1.Phase2** α_PB Abrikosov–Gorkov first-principles audit (H_AG_PARTIAL; SmH₉/YbH₉ targets) | 6 | CLOSED 2026-04-28 | [`research/op-sc-alpha-origin/Phase2_results.md`](research/op-sc-alpha-origin/Phase2_results.md) |
| **SC.1.Phase3** multi-LnH₉ falsification map (15 lantanowców; TGP RMS 0.42 < AG 1.53) | 7 | CLOSED 2026-04-28 | [`research/op-sc-alpha-origin/Phase3_results.md`](research/op-sc-alpha-origin/Phase3_results.md) |
| **BH.1.Phase1** multi-source dimensional + mass-scaling audit (M² scaling; H₀ rejected; Path E α(ψ) unique) | 5 | CLOSED 2026-04-28 | [`research/op-bh-alpha-threshold/Phase1_results.md`](research/op-bh-alpha-threshold/Phase1_results.md) |
| **BH.1.Phase2** substrate-physics upgrade of (α₀, n, ψ_th); ψ_th=1, n=2 DERIVED; α₀ ≈ κ_TGP² (0.75%) | 7 | CLOSED 2026-04-28 | [`research/op-bh-alpha-threshold/Phase2_results.md`](research/op-bh-alpha-threshold/Phase2_results.md) |
| **BH.1.Phase3** multi-source falsification map of α(ψ); 10-SMBH ngEHT, LIGO/LISA ringdown, NICER NS, MICROSCOPE-2 WEP, Cassini-class PPN, cross-sector √α₀=κ_TGP; 6 new predictions BH4–BH9 — **BH.1 program END** | 7 | CLOSED 2026-04-28 | [`research/op-bh-alpha-threshold/Phase3_results.md`](research/op-bh-alpha-threshold/Phase3_results.md) |
| **XS.1.Phase1** cross-sector identity √α₀ = κ_TGP feasibility audit (dimension, numerical 1σ, multi-anchor, data independence, Bayes factor 600); 5/5 PASS → Phase 2 viable | 5 | CLOSED 2026-04-28 | [`research/op-cross-sector-charge/Phase1_results.md`](research/op-cross-sector-charge/Phase1_results.md) |
| **XS.1.Phase2** substrate-action derivation of √α₀ = κ_TGP; single-Φ + Z₂ unify g_TGP across BH/SC; common-generator M_BH = M_SC = M_universal ≈ 4.03–4.05; RG-invariant ratio under common γ_an; **classification PARTIALLY DERIVED** (ξ unresolved: F4 rational 0.084% vs Phase 2 strict 0.747% rel diff) | 7 | CLOSED 2026-04-28 | [`research/op-cross-sector-charge/Phase2_results.md`](research/op-cross-sector-charge/Phase2_results.md) |
| **XS.1.Phase3** multi-sector falsification map of √α₀ = κ_TGP; XS1 ngEHT × SC v2 (≤5% by 2030+), XS2 g̃ orthogonality, XS3 lepton orthogonality, XS4 QM Born/CHSH orthogonality, XS5 F4 sub-percent (0.084%) match, XS6 6-channel roadmap; 6 new predictions XS1–XS6 — **XS.1 program END** | 7 | CLOSED 2026-04-28 | [`research/op-cross-sector-charge/Phase3_results.md`](research/op-cross-sector-charge/Phase3_results.md) |
| **ξ.1.Phase1** ξ-factor photon-ring frame audit (ξ_geom=1, α(α−1)=2, ψ_ph=1.168, F4 rational provenance, Phase 2 strict 11/97 sympy-exact); 5/5 PASS → all 5 inputs LOCKED w istniejących closurach; F4 vs strict split = 0.527% (z target_shift) | 5 | CLOSED 2026-04-29 | [`research/op-xi-photon-ring/Phase1_results.md`](research/op-xi-photon-ring/Phase1_results.md) |
| **ξ.1.Phase2** heat-kernel a₂ first-principles derivation (Birrell-Davies/Avramidi); a₂_vacuum = (1/2) V''(1)² = 2β² pod F2+F3+β=γ Z₂ vacuum; **Frame A (F4 0.114) = 1-loop a₂-corrected**, **Frame B (strict 11/97) = tree-level bare geometric**; 0.527% split = a₂ EFT 1-loop correction within ~1% EFT band; 6/7 PASS → XS.1 promotion PARTIALLY DERIVED → PARTIALLY DERIVED **(refined)**; full DERIVED czeka na UV-fixing N_A=8.7719 (closest 9, Δ 2.6%) | 7 | CLOSED 2026-04-29 | [`research/op-xi-photon-ring/Phase2_results.md`](research/op-xi-photon-ring/Phase2_results.md) |
| **ξ.1.Phase3** predictions + UV-route map for N_A normalization; ngEHT 2030+ 10-SMBH frame discrimination at ≥3σ (combined 0.158% < 0.176% threshold); AS NGFP najbliższy UV-route match dla N_A=8.7719 (Δ 0.068%); ξ-factor RG-invariant (0.022% drift); F4 1-loop reinterpretation, F5/F6 orthogonal; XS1 sharpening 5%→1.5% (factor 5.9×); 7-channel roadmap convergent 2030–2035; 6 new predictions XI1–XI6; UV7 promoted STRUCTURAL-POSTULATE → STRUCTURAL-DERIVED — **ξ.1 program END** | 7 | CLOSED 2026-04-29 | [`research/op-xi-photon-ring/Phase3_results.md`](research/op-xi-photon-ring/Phase3_results.md) |

Sub-cycle results for the most recent cycle:
[3.A KEYSTONE (AS NGFP)](research/op-phase3-uv-completion/Phase3_A_results.md) ·
[3.B (string / KKLT)](research/op-phase3-uv-completion/Phase3_B_results.md) ·
[3.C (LQG Ashtekar–Lewandowski)](research/op-phase3-uv-completion/Phase3_C_results.md) ·
[3.D (CDT Ambjørn–Loll)](research/op-phase3-uv-completion/Phase3_D_results.md) ·
[3.E (B.4/B.6/Δ_target)](research/op-phase3-uv-completion/Phase3_E_results.md) ·
[3.F CAPSTONE (4-of-4 UV synthesis)](research/op-phase3-uv-completion/Phase3_F_results.md) ·
[3.R-final](research/op-phase3-uv-completion/Phase3_R_final_results.md)

## Long-term research-track (open problems)

`research/op-uv-renormalizability-research/` — seven multi-year open
items, complementary to the structural compatibility closed by Phase 3:

1. Full UV-complete renormalizability proof
2. Selection of one of 4 UV completions
3. String vacuum-landscape selection (~10⁵⁰⁰ vacua)
4. LQG dynamics — Hamiltonian-constraint anomaly
5. CDT continuum-limit existence + Phase C universal class
6. Cosmological-constant problem — first-principles (beyond γ/12)
7. Empirical falsification at Planck-scale energies

## Falsification calendar (active observation windows)

Sourced from [`PREDICTIONS_REGISTRY.md`](PREDICTIONS_REGISTRY.md);
this is just the time-ordered scan of currently-LIVE entries.

| Window | Experiment | Predictions in registry |
|--------|-----------|----------------------------|
| **2027–2028** | MICROSCOPE-2 | G1 (η = 3.54·10⁻¹⁷) — clean shot at n=2; **BH7** (η_TGP = 2·10⁻¹⁸ from α(ψ_Earth), margin 5.1×) |
| **2027+** | LIGO O5 | GW2 (3 polarization DOF), GW5 (no vector), GW6 (dispersion bound), **BH5** (δf/f ~ 8–16% QNM ringdown) |
| **2027+** | NICER+ | **BH6** (NS M-R shift ~1–3% from GR; J0030/J0740) |
| **2027+** | DESI DR2 / DR3 | DE1 (w = −1), DE2 (w_a = 0), DE3 (T-Λ), C1 (H₀), C2 (S₈), C3 (Σm_ν) |
| **2027–2030** | LnH₉ DAC synthesis (Eremets/Hemley) | SC4 (SmH₉ T_c ~100 K), SC5 (YbH₉ T_c ~0.6 K), SC6 (TmH₉) — μ_eff² vs de Gennes; sharpens **XS1** κ_TGP precision to 0.3% |
| **2028+** | Euclid | DE3, DE4 (Friedmann ratio), C2 |
| **2030–2032** | ngEHT | BH1 (r_ph 1.293), BH2 (Δb_crit +14.56%), BH3, **BH4** (10-SMBH +14.56% multi-source map), **BH8** (√α₀ = κ_TGP cross-sector via α₀ from photon ring), **XS1** (combined ngEHT × SC v2 ≤5% → ≤1.5% post-ξ.1), **XS6** (6-channel roadmap convergence), **XI1** (Frame A vs B discrimination at ≥3σ z 10-SMBH), **XI4** (F4 1-loop reinterpretation), **XI5** (XS1 sharpening), **XI6** (7-channel convergence) |
| **~2035** | LISA / PTA | GW4 (m_σ²/m_s² = 2 → 2.9% low-k phase shift), **BH5** (LISA SMBH 10⁶–10⁷ M_⊙ ringdown), **XI3** (ξ-factor RG-invariance cross-scale) |
| **~2035** | LATOR / BEACON | **BH9** (γ−1 ~ 1.81·10⁻¹¹ at Sun surface; falsifiable below 10⁻¹⁰) |
| **research-track** | full QG / UV | UV1–UV7, F5, F6, **XI2** (UV-route map for N_A=8.7719: AS NGFP/LQG preferred) |

## Cross-references between flasks (logical dependencies)

```
tgp-core (DOI 19670324)
  ├── prerequisite for tgp-leptons (DOI 19706861)
  ├── prerequisite for tgp-qm     (DOI 19712596)
  └── prerequisite for tgp-sc v2  (DOI 19670557)
```

All three downstream flasks **inherit** the substrate-field framework
from tgp-core; none of them re-derive the axioms or the metric.
The four flasks are otherwise **independent** of each other.

## Audit trail policy

- **Deposit = immutable.** Once published to Zenodo, the deposit text is
  the timestamped prediction record. It is never silently rewritten.
- **In-repo notes = the audit log.** Each flask's `research/POST_PHASE3_*.md`
  is the live audit trail for that flask.
- **Master = the working state.** Discrepancies between master and any
  flask are logged in the relevant POST_PHASE3 note, with one of three
  outcomes: (a) flask amended in-repo (e.g. tgp-leptons r_21 micro-fix);
  (b) flask supplemented in-repo with addendum (e.g. tgp-core six
  strengthenings); (c) future v2 deposit on Zenodo when material warrants.
- **No silent retraction.** Falsified entries transition to FALSIFIED
  with a dated note and references in
  [`research/closure_2026-04-26/KNOWN_ISSUES.md`](research/closure_2026-04-26/KNOWN_ISSUES.md).

## Keeping this index up to date

When a new closure lands or a new flask deposits:

1. Add row to **Phase ledger** or **The four prediction-flask deposits** above.
2. Add (or upgrade status of) entry in [`PREDICTIONS_REGISTRY.md`](PREDICTIONS_REGISTRY.md).
3. Update grand total in **At a glance** + Phase-ledger table.
4. If a deposit was issued, add `POST_PHASE3_*` note to that flask's `research/`.
