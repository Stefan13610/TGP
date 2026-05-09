---
title: "Phase 0 balance sheet — op-ppE-mapping"
date: 2026-05-07
parent: "[[README.md]]"
type: phase0-balance
tgp_owner: research/op-ppE-mapping
tags:
  - phase0
  - balance-sheet
  - M03-gate
  - op-ppE-mapping
related:
  - "[[README.md]]"
  - "[[../../meta/CALIBRATION_PROTOCOL.md]]"
  - "[[../../audyt/M03_balance_sheet_missing/README.md]]"
---

# Phase 0 balance sheet — op-ppE-mapping

> **Cel pliku.** ABSOLUTE BINDING gate enforcement zgodnie z M03
> Phase 6 (post-2026-05-06): każdy nowy cykl `research/op-…/` MUSI
> dostarczyć Phase0_balance.md przed registry commit. Ten plik
> dokumentuje: (a) inputy cyklu (z explicit listą plików read +
> sympy/literature artifacts), (b) outputy i ich epistemic class,
> (c) predyktywność ratio.

## §1 — Inputy

### 1.1 Files read (audyt + research baseline)

| Plik | Sekcje | Wkład |
|------|--------|-------|
| [[../op-newton-momentum/M9_1_pp_P1_results.md]] | §2.1 (analytical c_n derivation), §3.2 (PN expansion), §4 (interpretation) | **Primary input** — c_n=2..7, sympy LOCK 5/5 |
| [[../op-newton-momentum/M9_1_pp_setup.md]] | §1 (M9.1'' postulate), §2 (algebraic derivation), §3 (physical interpretation), §2.5 (modified c²(ψ)) | M9.1'' canonical metric f(ψ) = (4−3ψ)/ψ |
| [[../op-newton-momentum/m9_1_pp_p1_higher_pn.py]] | sympy script | Reference computation method |
| [[../../audyt/T01_LIGO3G_falsifier/PPN_TO_PPE_MAPPING.md]] | §1, §2, §3, §5 | ppE framework + conventions + assumptions A1–A6 |
| [[../../audyt/T01_LIGO3G_falsifier/CONVENTION_DECISION.md]] | §1.2, §2 | PHASE convention adoption (b_ppE = −1 dla U³ w g_tt) |
| [[../../audyt/T01_LIGO3G_falsifier/CYCLE_KICKOFF_op-ppE-mapping.md]] | §3 (analytical setup), §5 (validation gates) | Phase 1 setup recipe |

### 1.2 Pre-existing locked anchors

| Anchor | Wartość | Source | Status pre-cykl |
|--------|---------|--------|-----------------|
| α (ODE-kinetic exponent) | 2 | F2 LOCKED w PREDICTIONS_REGISTRY (Phase1.A.1) | LOCKED |
| f(ψ) M9.1'' canonical | (4−3ψ)/ψ | M9_1_pp_setup §2.2 | DERIVED |
| h(ψ) M9.1'' canonical | ψ/(4−3ψ) | f·h = 1 budżet (sek08c) | DERIVED |
| c_n = a_n / a_1^n (n=2..7) | -1, +5/3, -10/3, +22/3, -154/9, +374/9 | M9_1_pp_P1 §2.1 sympy 5/5 | DERIVED |
| α_n^TGP w g_tt expansion | -2, +2, -7/3, +35/12, -91/24, +91/18 | M9_1_pp_P1 §3.2 | DERIVED |
| α_n^GR w g_tt expansion | -2, +2, -3/2, +1, -5/8, +3/8 | Schwarzschild izotropowy std | LOCKED-textbook |
| Δα_n = α_n^TGP − α_n^GR | 0, 0, **−5/6**, +23/12, −19/6, +337/72 | M9_1_pp_P1 §3.2 | DERIVED |

### 1.3 Literature inputs

| Reference | Rola |
|-----------|------|
| Yunes & Pretorius, Phys. Rev. D 80:122003 (2009) | ppE framework definition |
| Cutler & Flanagan, Phys. Rev. D 49:2658 (1994) | SPA inspiral phase formula |
| Blanchet, Living Rev. Relativ. 17:2 (2014) | PN waveform reference |
| Damour, Jaranowski, Schäfer, Phys. Rev. D 89:064058 (2014) | 4PN ADM Hamiltonian |
| Sampson, Yunes, Cornish, Phys. Rev. D 88:064056 (2013) | ppE-PN dictionary |

## §2 — Outputs (predicted contribution)

| Output | Plik | Epistemic class |
|--------|------|------------------|
| Two-body Lagrangian M9.1'' do v⁸ | [[Phase1_results.md]] §1 | DERIVED (sympy assist) |
| E_orb(v) reproduction Newton + 1PN | [[Phase1_results.md]] §2 | DERIVED (sympy LOCK 3/3) |
| E_orb(v) deviation w 2PN (z (5/6) U³) | [[Phase1_results.md]] §2 | DERIVED (sympy LOCK 1/1) |
| dE/dt luminosity z modyfikowaną quadrupole formula | [[Phase1_results.md]] §3 | DERIVED-PRELIMINARY (założenie A2 z PPN_TO_PPE_MAPPING) |
| SPA inversion → β_(N-PN_phase) coefficients | [[Phase1_results.md]] §4 | DERIVED |
| **β_ppE^TGP^(b=−1) liczbowo** | [[Phase1_results.md]] §4 | DERIVED (output głowny) |
| Multi-coefficient pattern {β_2PN, β_3PN, β_4PN} | [[Phase1_results.md]] §5 | DERIVED |
| Literature cross-check vs catalog | [[Phase2_literature_crosscheck.md]] | NUMERICAL (non-locking) |
| Paper-ready output | [[Phase3_paper_ready.md]] | SYNTHESIS |

## §3 — Predyktywność ratio

| Metric | Wartość |
|--------|---------|
| N_locked_inputs | 7 (α=2, f, h, c_n=2..7, α_n^TGP, α_n^GR, Δα_n) |
| N_outputs | 7 (Lagr, E_orb 1PN, E_orb 2PN, dE/dt, β_(2PN-phase), β_(3PN-phase), β_(4PN-phase), multi-coef pattern) |
| Ratio | 7/7 = **1.0** (≈ N=N, wyprowadzenie nie wprowadza nowych free parameters) |

**Komentarz:** ratio 1.0 oznacza że ten cykl jest *czystym propagatorem*
istniejących anchors do nowej domeny (waveform space) bez wprowadzania
nowych założeń liczbowych. Jedyne nowe założenia to A1–A6 z
[[../../audyt/T01_LIGO3G_falsifier/PPN_TO_PPE_MAPPING.md]] §5
(structural, nie numerical).

## §4 — Założenia walidowane w Phase 1

Z [[../../audyt/T01_LIGO3G_falsifier/PPN_TO_PPE_MAPPING.md]] §5:

| ID | Założenie | Walidacja w Phase 1 |
|----|-----------|----------------------|
| A1 | M9.1'' two-body Lagrangian istnieje (no horizon blow-up) | Phase 1.1: linear superposition ψ → 1 + ε_1(r_1) + ε_2(r_2) + ε_int (small interaction term) — VALID dla weak-field |
| A2 | Quadrupole formula struktura zachowana, deviation O((5/6)U³) wchodzi przez metric perturbation | Phase 1.3: dE/dt^TGP/dE/dt^GR − 1 = O((5/6) v⁶) ✓ STRUCTURAL |
| A3 | dE/dt luminosity modyfikuje 2PN+ orbital w sposób spójny z (5/6) U³ | Phase 1.3: ✓ propagacja z E_orb(v⁶) deviation |
| A4 | SPA stacjonarna w M9.1'' (adiabatic inspiral) | Phase 1.4: ✓ żaden mechanizm nie łamie adiabatic; spójne z ω_orb << M_pl |
| A5 | TGP nie wprowadza nowych radiacyjnych DOF | LIVE w PREDICTIONS_REGISTRY GW1, GW2 (3 DOF, c_T = c_s) |
| A6 | Konwencja PN counting "metric N-PN ↔ phase (N−1)-PN" | Phase 1.4: ✓ U³ w g_tt → b_ppE = −1 (2PN-phase) |

## §5 — Phase 0 sign-off

- ✓ Inputy zinwentaryzowane, wszystkie pre-existing locked.
- ✓ Outputy planowane z explicit epistemic class.
- ✓ Predyktywność ratio 1.0 — *propagator*, nie *new physics*.
- ✓ Założenia A1–A6 znane i będą walidowane w Phase 1.
- ✓ Konwencja PN: PHASE primary (z [[../../audyt/T01_LIGO3G_falsifier/CONVENTION_DECISION.md]]).
- ✓ M03 gate enforcement compliant.

**Phase 0 SIGNED 2026-05-07.** → Cykl wchodzi w Phase 1.
