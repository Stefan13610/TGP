---
title: "CANDIDATE_BRIDGES — luki analityczne potencjalnie domykalne wynikami z innych folderów"
date: 2026-05-03
type: bridges
status: ACTIVE — Sesja 6 PROPOSED matches; czeka na human review (HUMAN-CONFIRMED gate)
parent: "[[meta/PLAN_RESEARCH_WORKFLOW_v1.md]]"
related:
  - "[[meta/research/RESEARCH_BUS.md]]"
  - "[[meta/research/IMPACT_MATRIX.md]]"
  - "[[meta/research/AGENT_PROTOCOL.md]]"
tags:
  - bridges
  - matching
  - inter-folder
---

# CANDIDATE_BRIDGES — luki domykalne mostami między folderami

> **Sesja 6 auto-generated** (2026-05-03). Heurystyczny match `NEEDS.md` ×
> `FINDINGS.md` przez physics-aware tokenizer (zob.
> `_match_bridges_s6.py`). Każdy wpis to PROPOSED kandydat na bridge —
> czeka na HUMAN-CONFIRMED gate.

## 0. Reguły score (anty-spam)

| Strength | Wymaganie |
|----------|-----------|
| **EXACT** | (a) explicit cross-reference (folder mentioned by name) + ≥1 specific token, lub (b) hotspot ID match + ≥1 specific token |
| **PARTIAL** | explicit cross-ref LUB hotspot LUB phase marker (B6-CLOSED, A5-PATCHED, etc.) |
| **HEURISTIC** | ≥6 specific physics tokens overlap (np. M9, σ_ab, Φ_eq, β_PPN) |
| WEAK | (filtrowane out, nie są publikowane) |

## 1. Statystyka

| Metryka | Liczba |
|---------|-------:|
| PROPOSED total | 4 |
| EXACT | 0 |
| PARTIAL | 2 |
| HEURISTIC | 2 |
| HUMAN-CONFIRMED | 0 (czeka na review) |
| EXECUTED | 0 |

## 2. Top 3 bridges for human review

> Najsilniejsze kandydatury (highest total score). Każdy wymaga manualnej decyzji.

### BR-001: `op-delta1-g-tilde-derivation` ← `op-uv3-phi0-renormalization` (HEURISTIC, score=10)

- **Source NEEDS:** `research/op-delta1-g-tilde-derivation/NEEDS.md`
- **Target FINDINGS:** `research/op-uv3-phi0-renormalization/FINDINGS.md`
- **Strength:** HEURISTIC (xref=0, hotspot=0, phase=0, tokens=10)
- **Common signals:**
  - Specific tokens: ['n_c', 'sek00', 'α_s', 'φ_eff', 'ω_λ']
- **Decision:**
  - [ ] HUMAN-CONFIRMED
  - [ ] EXECUTE (copy F<id> z target FINDINGS do source FINDINGS z notatką pochodzenia)
  - [ ] REJECT (powód: ____)

### BR-002: `op-gamma1-phi-eff-anchor-resolution` ← `op-uv3-phi0-renormalization` (HEURISTIC, score=6)

- **Source NEEDS:** `research/op-gamma1-phi-eff-anchor-resolution/NEEDS.md`
- **Target FINDINGS:** `research/op-uv3-phi0-renormalization/FINDINGS.md`
- **Strength:** HEURISTIC (xref=0, hotspot=0, phase=0, tokens=6)
- **Common signals:**
  - Specific tokens: ['sek00', 'α_s', 'φ_eff']
- **Decision:**
  - [ ] HUMAN-CONFIRMED
  - [ ] EXECUTE (copy F<id> z target FINDINGS do source FINDINGS z notatką pochodzenia)
  - [ ] REJECT (powód: ____)

### BR-003: `atomic_shells_closure` ← `qm_superposition` (PARTIAL, score=5)

- **Source NEEDS:** `research/atomic_shells_closure/NEEDS.md`
- **Target FINDINGS:** `research/qm_superposition/FINDINGS.md`
- **Strength:** PARTIAL (xref=5, hotspot=0, phase=0, tokens=0)
- **Common signals:**
  - NEEDS mentions `qm_superposition` explicitly
- **Decision:**
  - [ ] HUMAN-CONFIRMED
  - [ ] EXECUTE (copy F<id> z target FINDINGS do source FINDINGS z notatką pochodzenia)
  - [ ] REJECT (powód: ____)

## 3. PROPOSED matches (full table)

| ID | Source NEEDS folder | → Target FINDINGS folder | Strength | Score | Top signal |
|----|---------------------|---------------------------|----------|------:|-----------|
| BR-001 | `op-delta1-g-tilde-derivation` | `op-uv3-phi0-renormalization` | HEURISTIC | 10 | Specific tokens: ['n_c', 'sek00', 'α_s', 'φ_eff', 'ω_λ'] |
| BR-002 | `op-gamma1-phi-eff-anchor-resolution` | `op-uv3-phi0-renormalization` | HEURISTIC | 6 | Specific tokens: ['sek00', 'α_s', 'φ_eff'] |
| BR-003 | `atomic_shells_closure` | `qm_superposition` | PARTIAL | 5 | NEEDS mentions `qm_superposition` explicitly |
| BR-004 | `mass_scaling_k4` | `op-lambda1-e2-amplitude-emergence` | PARTIAL | 5 | FINDINGS mentions `mass_scaling_k4` explicitly |

## 4. Per-source folder (consumers)

### `atomic_shells_closure` — 1 kandydat consumers

- BR-003: `qm_superposition` (PARTIAL, score=5)

### `mass_scaling_k4` — 1 kandydat consumers

- BR-004: `op-lambda1-e2-amplitude-emergence` (PARTIAL, score=5)

### `op-delta1-g-tilde-derivation` — 1 kandydat consumers

- BR-001: `op-uv3-phi0-renormalization` (HEURISTIC, score=10)

### `op-gamma1-phi-eff-anchor-resolution` — 1 kandydat consumers

- BR-002: `op-uv3-phi0-renormalization` (HEURISTIC, score=6)

## 5. Orphan NEEDS (45 folderów bez kandydata bridge)

> Foldery z NEEDS, ale żaden FINDINGS w innym folderze nie matchuje powyżej WEAK threshold.
> To są **realnie open** luki — wymagają NOWEJ pracy, nie bridge.

- `research/atom_from_soliton`
- `research/brannen_sqrt2`
- `research/cabibbo_correction`
- `research/closure_2026-04-26`
- `research/cosmo_tensions`
- `research/em_from_substrate`
- `research/galaxy_scaling`
- `research/liquid_viscosity`
- `research/metric_ansatz`
- `research/neutrino_msw`
- `research/op-alpha-fine-structure`
- `research/op-bh-alpha-threshold`
- `research/op-cross-sector-charge`
- `research/op-delta2-Nf-derivation`
- `research/op-eps-photon-ring`
- `research/op-eta-wolfenstein`
- `research/op-g0-r3-from-canonical-projection`
- `research/op-iota-charge-pmns-unification`
- `research/op-kappa-mixing-numerator`
- `research/op-lambda1-e2-amplitude-emergence`
- `research/op-mu-pmns-phase-hardening`
- `research/op-mu1-minimal-substrate-log-redefinition`
- `research/op-omega1-substrate-em-coupling`
- `research/op-phase1-covariant`
- `research/op-phase2-quantum-gravity`
- `research/op-phase3-uv-completion`
- `research/op-psi1-substrate-light-acceleration`
- `research/op-rho1-71Ge-cross-section`
- `research/op-sc-alpha-origin`
- `research/op-sigma1-substrate-light-dispersion`
- `research/op-tau2-substrate-time-coupling`
- `research/op-tau3-substrate-clock-acceleration`
- `research/op-theta-quark-koide`
- `research/op-uv-as-ngfp`
- `research/op-uv-renormalizability-research`
- `research/op-uv3-phi0-renormalization`
- `research/op-xi-photon-ring`
- `research/op1-op2-op4`
- `research/particle_sector_closure`
- `research/qm_foundations`
- `research/qm_measurement`
- `research/qm_superposition`
- `research/superconductivity_closure`
- `research/thermal_transport_molecular`
- `research/why_n3`

## 6. Orphan FINDINGS (67 folderów eksportujących bez konsumenta)

> Foldery z `exports_findings: true`, ale żaden NEEDS nie pyta o ich wyniki.
> Po 90 dniach: degradacja w `RESEARCH_BUS.md` na STALE.

- `research/brannen_sqrt2`
- `research/cabibbo_correction`
- `research/casimir_mof`
- `research/closure_2026-04-26`
- `research/continuum_limit`
- `research/cosmo_tensions`
- `research/external_review_2026-04-25`
- `research/galaxy_scaling`
- `research/hubble_tension`
- `research/liquid_viscosity`
- `research/mass_scaling_k4`
- `research/metric_ansatz`
- `research/muon_g_minus_2`
- `research/neutrino_msw`
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
- `research/op-xi-photon-ring`
- `research/op-xi2-sterile-nu-5sector`
- `research/op-zeta-mass-spectrum`
- `research/op7`
- `research/particle_sector_closure`
- `research/qm_born_rule`
- `research/qm_decoherence`
- `research/qm_entanglement`
- `research/qm_foundations`
- `research/qm_measurement`
- `research/qm_spin`
- `research/qm_statistics`
- `research/s8_tension`
- `research/superconductivity_closure`
- `research/thermal_transport_molecular`
- `research/uv_completion`

## 7. Wyłączone z matching

- 4 polluted-74394a8 folders (kwarantanna; AGENT_PROTOCOL §3.2)
- 33 folderów z empty NEEDS
- 16 folderów z `exports_findings: false` (jako sources nie eksportują dla matchingu)

## 8. Generation metadata

- Generated: 2026-05-03T14:31:44
- Builder: `_match_bridges_s6.py`
- Folders processed: 86
- All-pairs scored: 7310
- Matches above WEAK threshold: 4