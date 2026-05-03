---
title: "FOLDER_STATUS_INDEX — globalna mapa wszystkich folderów research/"
date: 2026-05-03
type: index
status: GENERATED (Sesja 4) — auto-build z YAML w research/<X>/README.md
parent: "[[meta/PLAN_RESEARCH_WORKFLOW_v1.md]]"
related:
  - "[[meta/research/AGENT_PROTOCOL.md]]"
  - "[[meta/research/AUDIT_RESEARCH_S1.md]]"
  - "[[meta/research/templates/STATUS_BLOCK.yaml]]"
  - "[[meta/research/CORE_CANDIDATES.md]]"
  - "[[meta/research/_classify_s4.py]]"
tags:
  - index
  - folder-status
  - generated
---

# FOLDER_STATUS_INDEX — globalna mapa wszystkich folderów `research/`

> **Auto-generated z `_classify_s4.py` + `_build_status_index.py`** —
> Sesja 4 (2026-05-03). Source of truth: YAML `tgp_status:` w
> `research/<X>/README.md` per folder (po `--apply`).
>
> **Regenerate:** `python meta/research/_build_status_index.py`
>
> **Liczba folderów:** 86 (po wyłączeniu `_sandbox/` i `_archive/`).

## 1. Rozkład statusów

### 1.1 `folder_status`

| Status | Liczba | % |
|--------|-------:|---:|
| `active` | 77 | 89% |
| `needs-bridge` | 8 | 9% |
| `review` | 1 | 1% |

### 1.2 `level`

| Level | Liczba | % |
|-------|-------:|---:|
| `L1` | 31 | 36% |
| `L2` | 18 | 20% |
| `L3` | 10 | 11% |
| `mixed` | 5 | 5% |
| `unknown` | 22 | 25% |

### 1.3 `kind`

| Kind | Liczba | % |
|------|-------:|---:|
| `derivation` | 71 | 82% |
| `phenomenology` | 13 | 15% |
| `closure-aggregator` | 1 | 1% |
| `review` | 1 | 1% |

### 1.4 `core_compatibility`

| Compatibility | Liczba | % |
|---------------|-------:|---:|
| `unknown` | 77 | 89% |
| `current` | 8 | 9% |
| `partial` | 1 | 1% |

## 2. Flagi

- **Polluted-74394a8** (kwarantanna): **4** folderów
- **Banned phrase 'FULL CONVERGENCE'**: **21** folderów (zob. § 5)

## 3. Pełna mapa (sortowana alfabetycznie)

| # | Folder | folder_status | level | kind | core_compat | mtime | flags |
|---|--------|---------------|-------|------|-------------|-------|-------|
| 1 | `atom_from_soliton` | active | unknown | derivation | unknown | 2026-05-03 | — |
| 2 | `atomic_shells_closure` | active | L1 | phenomenology | unknown | 2026-05-03 | — |
| 3 | `brannen_sqrt2` | active | unknown | derivation | unknown | 2026-05-03 | — |
| 4 | `cabibbo_correction` | active | unknown | derivation | unknown | 2026-05-03 | — |
| 5 | `casimir_mof` | active | unknown | phenomenology | unknown | 2026-05-03 | — |
| 6 | `closure_2026-04-26` | active | mixed | closure-aggregator | partial | 2026-05-03 | — |
| 7 | `cohesion_closure` | active | unknown | derivation | unknown | 2026-05-03 | — |
| 8 | `continuum_limit` | active | L1 | derivation | unknown | 2026-05-03 | — |
| 9 | `cosmo_tensions` | active | L1 | phenomenology | unknown | 2026-05-03 | — |
| 10 | `desi_dark_energy` | active | L1 | phenomenology | unknown | 2026-05-03 | — |
| 11 | `em_from_substrate` | active | unknown | derivation | unknown | 2026-05-03 | — |
| 12 | `external_review_2026-04-25` | review | mixed | review | current | 2026-05-03 | — |
| 13 | `galaxy_scaling` | active | L1 | phenomenology | unknown | 2026-05-03 | — |
| 14 | `hubble_tension` | active | L1 | phenomenology | unknown | 2026-05-03 | — |
| 15 | `liquid_viscosity` | active | unknown | phenomenology | unknown | 2026-05-03 | — |
| 16 | `mass_scaling_k4` | active | L1 | derivation | unknown | 2026-05-03 | — |
| 17 | `metric_ansatz` | active | unknown | derivation | unknown | 2026-05-03 | — |
| 18 | `muon_g_minus_2` | active | unknown | phenomenology | unknown | 2026-05-03 | — |
| 19 | `nbody` | active | L1 | phenomenology | current | 2026-05-03 | audit×8 |
| 20 | `neutrino_msw` | needs-bridge | L1 | derivation | unknown | 2026-05-03 | — |
| 21 | `op-alpha-fine-structure` | active | L3 | derivation | unknown | 2026-05-03 | — |
| 22 | `op-bh-alpha-threshold` | active | L2 | derivation | unknown | 2026-05-03 | — |
| 23 | `op-chi1-newton-constant-derivation` | needs-bridge | unknown | derivation | unknown | 2026-05-03 | ⚠POLLUTED ⚠FULLCONV×17 |
| 24 | `op-cosmology-closure` | active | L1 | phenomenology | unknown | 2026-05-03 | — |
| 25 | `op-cross-sector-charge` | active | L3 | derivation | unknown | 2026-05-03 | — |
| 26 | `op-delta1-g-tilde-derivation` | active | L1 | derivation | unknown | 2026-05-03 | — |
| 27 | `op-delta2-Nf-derivation` | active | L1 | derivation | unknown | 2026-05-03 | — |
| 28 | `op-eht` | active | L1 | derivation | unknown | 2026-05-03 | — |
| 29 | `op-eht-A` | active | L1 | derivation | unknown | 2026-05-03 | — |
| 30 | `op-eps-photon-ring` | active | L3 | derivation | unknown | 2026-05-03 | — |
| 31 | `op-eta-wolfenstein` | active | L3 | derivation | unknown | 2026-05-03 | — |
| 32 | `op-eta2-denom-derivation` | active | L3 | derivation | unknown | 2026-05-03 | — |
| 33 | `op-g0-r3-from-canonical-projection` | active | L1 | derivation | current | 2026-05-03 | audit×10 |
| 34 | `op-gamma1-phi-eff-anchor-resolution` | active | L1 | derivation | unknown | 2026-05-03 | — |
| 35 | `op-iota-charge-pmns-unification` | active | L2 | derivation | unknown | 2026-05-03 | ⚠FULLCONV×8 |
| 36 | `op-kappa-mixing-numerator` | active | L2 | derivation | unknown | 2026-05-03 | — |
| 37 | `op-lambda1-e2-amplitude-emergence` | active | L1 | derivation | unknown | 2026-05-03 | ⚠FULLCONV×2 |
| 38 | `op-m92` | active | L1 | derivation | unknown | 2026-05-03 | — |
| 39 | `op-mu-pmns-phase-hardening` | active | L2 | derivation | unknown | 2026-05-03 | ⚠FULLCONV×8 |
| 40 | `op-mu1-minimal-substrate-log-redefinition` | active | L1 | derivation | unknown | 2026-05-03 | — |
| 41 | `op-newton-momentum` | active | L1 | derivation | current | 2026-05-03 | audit×11 |
| 42 | `op-nu-majorana-phase-mbb` | active | L2 | derivation | unknown | 2026-05-03 | ⚠FULLCONV×8 |
| 43 | `op-omega1-substrate-em-coupling` | active | L2 | derivation | unknown | 2026-05-03 | ⚠FULLCONV×6 |
| 44 | `op-omega2-axion-coupling-lock` | needs-bridge | unknown | derivation | unknown | 2026-05-03 | ⚠POLLUTED ⚠FULLCONV×2 |
| 45 | `op-omega3-axion-decay-constant` | needs-bridge | unknown | derivation | unknown | 2026-05-03 | ⚠POLLUTED ⚠FULLCONV×7 audit×1 |
| 46 | `op-omicron1-sigmamnu-cosmo` | active | L2 | phenomenology | unknown | 2026-05-03 | ⚠FULLCONV×8 |
| 47 | `op-phase1-covariant` | active | L1 | derivation | unknown | 2026-05-03 | — |
| 48 | `op-phase2-quantum-gravity` | active | L1 | derivation | unknown | 2026-05-03 | — |
| 49 | `op-phase3-uv-completion` | active | unknown | derivation | unknown | 2026-05-03 | — |
| 50 | `op-phi1-substrate-action-variational` | active | L2 | derivation | unknown | 2026-05-03 | ⚠FULLCONV×6 |
| 51 | `op-pi1-bb0nu-nme-isotope` | active | L2 | derivation | unknown | 2026-05-03 | ⚠FULLCONV×6 |
| 52 | `op-psi1-substrate-light-acceleration` | active | L2 | derivation | unknown | 2026-05-03 | ⚠FULLCONV×4 |
| 53 | `op-quantum-closure` | active | L1 | derivation | unknown | 2026-05-03 | — |
| 54 | `op-rho1-71Ge-cross-section` | active | L2 | derivation | unknown | 2026-05-03 | ⚠FULLCONV×11 |
| 55 | `op-sc-alpha-origin` | active | L3 | derivation | unknown | 2026-05-03 | — |
| 56 | `op-sigma1-substrate-light-dispersion` | active | L2 | derivation | unknown | 2026-05-03 | ⚠FULLCONV×7 |
| 57 | `op-tau1-closure-overlap-coulomb` | active | L2 | derivation | unknown | 2026-05-03 | ⚠FULLCONV×10 |
| 58 | `op-tau2-substrate-time-coupling` | active | L2 | derivation | unknown | 2026-05-03 | ⚠FULLCONV×4 |
| 59 | `op-tau3-substrate-clock-acceleration` | active | L2 | derivation | current | 2026-05-03 | ⚠FULLCONV×3 audit×21 |
| 60 | `op-theta-quark-koide` | active | L3 | derivation | unknown | 2026-05-03 | — |
| 61 | `op-upsilon1-closure-cross-family` | active | L2 | derivation | unknown | 2026-05-03 | ⚠FULLCONV×7 |
| 62 | `op-uv-as-ngfp` | active | L3 | derivation | unknown | 2026-05-03 | — |
| 63 | `op-uv-renormalizability-research` | needs-bridge | L1 | derivation | unknown | 2026-05-03 | — |
| 64 | `op-uv2-mtgp-absolute-scale` | needs-bridge | unknown | derivation | unknown | 2026-05-03 | ⚠POLLUTED ⚠FULLCONV×16 |
| 65 | `op-uv3-phi0-renormalization` | active | L2 | derivation | unknown | 2026-05-03 | ⚠FULLCONV×9 |
| 66 | `op-xi-photon-ring` | active | L3 | derivation | unknown | 2026-05-03 | — |
| 67 | `op-xi2-sterile-nu-5sector` | active | L2 | derivation | unknown | 2026-05-03 | ⚠FULLCONV×9 |
| 68 | `op-zeta-mass-spectrum` | active | L3 | derivation | unknown | 2026-05-03 | — |
| 69 | `op1-op2-op4` | active | mixed | derivation | current | 2026-05-03 | — |
| 70 | `op6` | active | mixed | derivation | current | 2026-05-03 | — |
| 71 | `op7` | active | mixed | derivation | current | 2026-05-03 | — |
| 72 | `particle_sector_closure` | active | L1 | derivation | unknown | 2026-05-03 | — |
| 73 | `qm_born_rule` | needs-bridge | L1 | derivation | unknown | 2026-05-03 | — |
| 74 | `qm_decoherence` | active | unknown | derivation | unknown | 2026-05-03 | — |
| 75 | `qm_entanglement` | active | unknown | derivation | unknown | 2026-05-03 | — |
| 76 | `qm_foundations` | active | unknown | derivation | unknown | 2026-05-03 | — |
| 77 | `qm_measurement` | active | unknown | derivation | unknown | 2026-05-03 | — |
| 78 | `qm_spin` | active | unknown | derivation | unknown | 2026-05-03 | — |
| 79 | `qm_statistics` | active | unknown | derivation | unknown | 2026-05-03 | — |
| 80 | `qm_superposition` | active | unknown | derivation | unknown | 2026-05-03 | — |
| 81 | `rho_normal_state_closure` | active | L1 | derivation | unknown | 2026-05-03 | — |
| 82 | `s8_tension` | active | L1 | phenomenology | unknown | 2026-05-03 | — |
| 83 | `superconductivity_closure` | active | L1 | derivation | unknown | 2026-05-03 | — |
| 84 | `thermal_transport_molecular` | active | unknown | phenomenology | unknown | 2026-05-03 | — |
| 85 | `uv_completion` | needs-bridge | L1 | derivation | unknown | 2026-05-03 | — |
| 86 | `why_n3` | active | L1 | derivation | unknown | 2026-05-03 | — |

## 4. Per `folder_status`

### 4.x. `active` (77)

- `research/atom_from_soliton`
- `research/atomic_shells_closure`
- `research/brannen_sqrt2`
- `research/cabibbo_correction`
- `research/casimir_mof`
- `research/closure_2026-04-26`
- `research/cohesion_closure`
- `research/continuum_limit`
- `research/cosmo_tensions`
- `research/desi_dark_energy`
- `research/em_from_substrate`
- `research/galaxy_scaling`
- `research/hubble_tension`
- `research/liquid_viscosity`
- `research/mass_scaling_k4`
- `research/metric_ansatz`
- `research/muon_g_minus_2`
- `research/nbody`
- `research/op-alpha-fine-structure`
- `research/op-bh-alpha-threshold`
- `research/op-cosmology-closure`
- `research/op-cross-sector-charge`
- `research/op-delta1-g-tilde-derivation`
- `research/op-delta2-Nf-derivation`
- `research/op-eht`
- `research/op-eht-A`
- `research/op-eps-photon-ring`
- `research/op-eta-wolfenstein`
- `research/op-eta2-denom-derivation`
- `research/op-g0-r3-from-canonical-projection`
- `research/op-gamma1-phi-eff-anchor-resolution`
- `research/op-iota-charge-pmns-unification`
- `research/op-kappa-mixing-numerator`
- `research/op-lambda1-e2-amplitude-emergence`
- `research/op-m92`
- `research/op-mu-pmns-phase-hardening`
- `research/op-mu1-minimal-substrate-log-redefinition`
- `research/op-newton-momentum`
- `research/op-nu-majorana-phase-mbb`
- `research/op-omega1-substrate-em-coupling`
- `research/op-omicron1-sigmamnu-cosmo`
- `research/op-phase1-covariant`
- `research/op-phase2-quantum-gravity`
- `research/op-phase3-uv-completion`
- `research/op-phi1-substrate-action-variational`
- `research/op-pi1-bb0nu-nme-isotope`
- `research/op-psi1-substrate-light-acceleration`
- `research/op-quantum-closure`
- `research/op-rho1-71Ge-cross-section`
- `research/op-sc-alpha-origin`
- `research/op-sigma1-substrate-light-dispersion`
- `research/op-tau1-closure-overlap-coulomb`
- `research/op-tau2-substrate-time-coupling`
- `research/op-tau3-substrate-clock-acceleration`
- `research/op-theta-quark-koide`
- `research/op-upsilon1-closure-cross-family`
- `research/op-uv-as-ngfp`
- `research/op-uv3-phi0-renormalization`
- `research/op-xi-photon-ring`
- `research/op-xi2-sterile-nu-5sector`
- `research/op-zeta-mass-spectrum`
- `research/op1-op2-op4`
- `research/op6`
- `research/op7`
- `research/particle_sector_closure`
- `research/qm_decoherence`
- `research/qm_entanglement`
- `research/qm_foundations`
- `research/qm_measurement`
- `research/qm_spin`
- `research/qm_statistics`
- `research/qm_superposition`
- `research/rho_normal_state_closure`
- `research/s8_tension`
- `research/superconductivity_closure`
- `research/thermal_transport_molecular`
- `research/why_n3`

### 4.x. `needs-bridge` (8)

- `research/neutrino_msw`
- `research/op-chi1-newton-constant-derivation`
- `research/op-omega2-axion-coupling-lock`
- `research/op-omega3-axion-decay-constant`
- `research/op-uv-renormalizability-research`
- `research/op-uv2-mtgp-absolute-scale`
- `research/qm_born_rule`
- `research/uv_completion`

### 4.x. `review` (1)

- `research/external_review_2026-04-25`

## 5. Per `level`

### 5.x. `L1` (31)

- `research/atomic_shells_closure`
- `research/continuum_limit`
- `research/cosmo_tensions`
- `research/desi_dark_energy`
- `research/galaxy_scaling`
- `research/hubble_tension`
- `research/mass_scaling_k4`
- `research/nbody`
- `research/neutrino_msw`
- `research/op-cosmology-closure`
- `research/op-delta1-g-tilde-derivation`
- `research/op-delta2-Nf-derivation`
- `research/op-eht`
- `research/op-eht-A`
- `research/op-g0-r3-from-canonical-projection`
- `research/op-gamma1-phi-eff-anchor-resolution`
- `research/op-lambda1-e2-amplitude-emergence`
- `research/op-m92`
- `research/op-mu1-minimal-substrate-log-redefinition`
- `research/op-newton-momentum`
- `research/op-phase1-covariant`
- `research/op-phase2-quantum-gravity`
- `research/op-quantum-closure`
- `research/op-uv-renormalizability-research`
- `research/particle_sector_closure`
- `research/qm_born_rule`
- `research/rho_normal_state_closure`
- `research/s8_tension`
- `research/superconductivity_closure`
- `research/uv_completion`
- `research/why_n3`

### 5.x. `L2` (18)

- `research/op-bh-alpha-threshold`
- `research/op-iota-charge-pmns-unification`
- `research/op-kappa-mixing-numerator`
- `research/op-mu-pmns-phase-hardening`
- `research/op-nu-majorana-phase-mbb`
- `research/op-omega1-substrate-em-coupling`
- `research/op-omicron1-sigmamnu-cosmo`
- `research/op-phi1-substrate-action-variational`
- `research/op-pi1-bb0nu-nme-isotope`
- `research/op-psi1-substrate-light-acceleration`
- `research/op-rho1-71Ge-cross-section`
- `research/op-sigma1-substrate-light-dispersion`
- `research/op-tau1-closure-overlap-coulomb`
- `research/op-tau2-substrate-time-coupling`
- `research/op-tau3-substrate-clock-acceleration`
- `research/op-upsilon1-closure-cross-family`
- `research/op-uv3-phi0-renormalization`
- `research/op-xi2-sterile-nu-5sector`

### 5.x. `L3` (10)

- `research/op-alpha-fine-structure`
- `research/op-cross-sector-charge`
- `research/op-eps-photon-ring`
- `research/op-eta-wolfenstein`
- `research/op-eta2-denom-derivation`
- `research/op-sc-alpha-origin`
- `research/op-theta-quark-koide`
- `research/op-uv-as-ngfp`
- `research/op-xi-photon-ring`
- `research/op-zeta-mass-spectrum`

### 5.x. `mixed` (5)

- `research/closure_2026-04-26`
- `research/external_review_2026-04-25`
- `research/op1-op2-op4`
- `research/op6`
- `research/op7`

### 5.x. `unknown` (22)

- `research/atom_from_soliton`
- `research/brannen_sqrt2`
- `research/cabibbo_correction`
- `research/casimir_mof`
- `research/cohesion_closure`
- `research/em_from_substrate`
- `research/liquid_viscosity`
- `research/metric_ansatz`
- `research/muon_g_minus_2`
- `research/op-chi1-newton-constant-derivation`
- `research/op-omega2-axion-coupling-lock`
- `research/op-omega3-axion-decay-constant`
- `research/op-phase3-uv-completion`
- `research/op-uv2-mtgp-absolute-scale`
- `research/qm_decoherence`
- `research/qm_entanglement`
- `research/qm_foundations`
- `research/qm_measurement`
- `research/qm_spin`
- `research/qm_statistics`
- `research/qm_superposition`
- `research/thermal_transport_molecular`

## 6. Banned phrase 'FULL CONVERGENCE' violators (Q7 finding)

> Per [[meta/research/AGENT_PROTOCOL.md]] §3.0/§3.1: 'FULL CONVERGENCE' jest
> zarezerwowane. Ta lista wymaga **phrase cleanup w future maintenance pass**
> (zalecenie z [[meta/research/SANITY_PASS_S3_5.md]] §4.1).

| Folder | Hits | Files (top 3) | Status |
|--------|-----:|---------------|--------|
| `op-chi1-newton-constant-derivation` | 17 | `CRITIQUE_circular_anchor_2026-05-02.md`, `Phase2_results.md`, `phase3_predictions.py` | POLLUTED |
| `op-iota-charge-pmns-unification` | 8 | `FINDINGS.md`, `phase3_iota_predictions.py`, `phase3_iota_predictions.txt` | active |
| `op-lambda1-e2-amplitude-emergence` | 2 | `README.md` | active |
| `op-mu-pmns-phase-hardening` | 8 | `FINDINGS.md`, `phase3_mu_predictions.py`, `phase3_mu_predictions.txt` | active |
| `op-nu-majorana-phase-mbb` | 8 | `phase3_predictions.py`, `phase3_predictions.txt`, `Phase3_results.md` | active |
| `op-omega1-substrate-em-coupling` | 6 | `phase3_omega1_predictions.py`, `Phase3_results.md` | active |
| `op-omega2-axion-coupling-lock` | 2 | `program.md`, `README.md` | POLLUTED |
| `op-omega3-axion-decay-constant` | 7 | `phase3_predictions.py`, `phase3_predictions.txt`, `Phase3_results.md` | POLLUTED |
| `op-omicron1-sigmamnu-cosmo` | 8 | `Phase3_results.md`, `Phase3_setup.md`, `phase3_sigmamnu_predictions.py` | active |
| `op-phi1-substrate-action-variational` | 6 | `phase3_phi1_predictions.py`, `Phase3_results.md` | active |
| `op-pi1-bb0nu-nme-isotope` | 6 | `phase3_nme_predictions.py`, `phase3_nme_predictions.txt`, `Phase3_results.md` | active |
| `op-psi1-substrate-light-acceleration` | 4 | `phase3_psi1_predictions.py`, `Phase3_results.md` | active |
| `op-rho1-71Ge-cross-section` | 11 | `phase3_71Ge_predictions.py`, `phase3_71Ge_predictions.txt`, `Phase3_results.md` | active |
| `op-sigma1-substrate-light-dispersion` | 7 | `Phase3_results.md`, `Phase3_setup.md`, `phase3_sigma1_predictions.py` | active |
| `op-tau1-closure-overlap-coulomb` | 10 | `phase1_tau1_landscape.py`, `Phase3_results.md`, `phase3_tau1_predictions.py` | active |
| `op-tau2-substrate-time-coupling` | 4 | `Phase3_results.md`, `Phase3_setup.md`, `phase3_tau2_predictions.py` | active |
| `op-tau3-substrate-clock-acceleration` | 3 | `Phase3_results.md`, `phase3_tau3_predictions.py` | active |
| `op-upsilon1-closure-cross-family` | 7 | `Phase3_results.md`, `phase3_upsilon1_predictions.py`, `program.md` | active |
| `op-uv2-mtgp-absolute-scale` | 16 | `CRITIQUE_repackaged_circularity_2026-05-02.md`, `phase3_predictions.py`, `phase3_predictions.txt` | POLLUTED |
| `op-uv3-phi0-renormalization` | 9 | `FINDINGS.md`, `Phase3_results.md`, `README.md` | active |
| `op-xi2-sterile-nu-5sector` | 9 | `phase3_predictions.py`, `phase3_predictions.txt`, `Phase3_results.md` | active |

## 7. Foldery z audit-aware markers (B6/B7/B8/B9/A5 patches)

> Foldery które zawierają patche z audytu 2026-05-01 (B6-CLOSED, B8-CLOSED,
> B9-CLOSED, A5-PATCHED, etc.). To są foldery o **wzmocnionym** zaufaniu —
> ich `core_compatibility: current` jest udokumentowany.

| Folder | Audit hits | folder_status | level |
|--------|-----------:|---------------|-------|
| `op-tau3-substrate-clock-acceleration` | 21 | active | L2 |
| `op-newton-momentum` | 11 | active | L1 |
| `op-g0-r3-from-canonical-projection` | 10 | active | L1 |
| `nbody` | 8 | active | L1 |
| `op-omega3-axion-decay-constant` | 1 | needs-bridge | unknown |

## 8. Anti-overclaim audit (CRITICAL flagi)

- ✅ **0 CRITICAL findings** — wszystkie foldery przechodzą basic anti-overclaim audit.

Twarde reguły potwierdzone:
- 0 folderów z `level: L4` (zakazane bez INTAKE workflow w Sesji 7)
- 0 folderów z niepustym `promoted_to_core` w Sesji 4
- Każdy folder z `level >= L2` ma niepusty `source_of_status`

## 9. Generation metadata

- Generated: 2026-05-03T14:22:33
- Source: `_classify_s4_dryrun.json` (matches applied YAMLs)
- Builder: `_build_status_index.py`
- Folders: 86