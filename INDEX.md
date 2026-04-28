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
---

# TGP_v1 — Index hub

A one-screen navigation map for the TGP research ecosystem
post-Phase-3 (closure 2026-04-28). For the project overview, see
[`README.md`](README.md). This file is the **navigation entry point**.

## At a glance

- **Master verification ledger:** **291** cumulative closures
  (M9 13 + M10 42 + M11 62 + Phase 1 50 + Phase 2 54 + Phase 3 60 + SC.1.Phase1 4 + SC.1.Phase2 6;
  closure_2026-04-26's 35 tests are tracked in the phase-ledger table but
  overlap conceptually with Phase 1/2/3 follow-ups so are not double-counted here).
- **All major cycles CLOSED** as of 2026-04-28.
- **4 Zenodo flask papers** deposited as immutable timestamped predictions.
- **1 predictions registry** ([`PREDICTIONS_REGISTRY.md`](PREDICTIONS_REGISTRY.md))
  cross-linking every flask to falsification target + horizon.
- **Active queued work:** SC.1.Phase3 (multi-LnH₉ validation — SmH₉/YbH₉ falsification targets).

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

## Phase ledger (291 cumulative)

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
| **2027–2028** | MICROSCOPE-2 | G1 (η = 3.54·10⁻¹⁷) — clean shot at n=2 |
| **2027+** | LIGO O5 | GW2 (3 polarization DOF), GW5 (no vector), GW6 (dispersion bound) |
| **2027+** | DESI DR2 / DR3 | DE1 (w = −1), DE2 (w_a = 0), DE3 (T-Λ), C1 (H₀), C2 (S₈), C3 (Σm_ν) |
| **2028+** | Euclid | DE3, DE4 (Friedmann ratio), C2 |
| **2030–2032** | ngEHT | BH1 (r_ph 1.293), BH2 (Δb_crit +14.56%), BH3 |
| **~2035** | LISA / PTA | GW4 (m_σ²/m_s² = 2 → 2.9% low-k phase shift) |
| **research-track** | full QG / UV | UV1–UV7, F5, F6 |

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
