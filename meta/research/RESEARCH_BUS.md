---
title: "RESEARCH_BUS — tablica ogłoszeń międzyfolderowych"
date: 2026-05-03
type: bus
status: ACTIVE — Sesja 5 initial broadcasts; Sesja 6 fills consumers via CANDIDATE_BRIDGES
parent: "[[meta/PLAN_RESEARCH_WORKFLOW_v1.md]]"
related:
  - "[[meta/research/AGENT_PROTOCOL.md]]"
  - "[[meta/research/CANDIDATE_BRIDGES.md]]"
  - "[[meta/research/IMPACT_MATRIX.md]]"
  - "[[meta/research/FOLDER_STATUS_INDEX.md]]"
tags:
  - research-bus
  - broadcast
  - inter-folder
---

# RESEARCH_BUS — tablica ogłoszeń międzyfolderowych

> **Cel pliku:** każdy nowy wynik (item w `FINDINGS.md`) generuje wpis
> broadcastowy z polem `consumers:` (foldery, które mogą skonsumować).
> Auto-broadcasted by Sesja 5 (per-folder podsumowanie); Sesja 6 dopisuje
> consumers przez `CANDIDATE_BRIDGES.md` matching.

## 0. Reguły wpisu

1. **Trigger**: dodanie/aktualizacja `FINDINGS.md` w folderze badawczym
   (z `exports_findings: true` i `polluted_74394a8: false`).
2. **Format**: jeden wiersz tabeli na folder source.
3. **`Consumers`**: lista folderów z `tgp_status.depends_on:` matching
   ten folder LUB Sesja 6 bridge matching. Pusta dopuszczalna,
   ale wymaga § 3 entry.
4. **Status**:
   - `BROADCAST` — wpis świeży, nieprzeczytany przez konsumentów
   - `CONSUMED` — co najmniej jeden konsument potwierdził użycie
   - `STALE` — wpis > 90 dni bez akcji

## 1. Statystyka

| Metryka | Liczba |
|---------|-------:|
| Aktywnych BROADCAST | 70 |
| CONSUMED | 0 |
| STALE | 0 |
| **Razem (wszystkie czasy)** | **70** |

Last build: 2026-05-03T14:22:08.

## 2. Aktywne ogłoszenia (BROADCAST)

| Date | Source folder | Headline | Score | Findings link | Consumers | Status |
|------|---------------|----------|-------|---------------|-----------|--------|
| 2026-05-03 | `brannen_sqrt2` | > | — | [[research/brannen_sqrt2/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `cabibbo_correction` | > | — | [[research/cabibbo_correction/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `casimir_mof` | **`PLAN.md`:** | — | [[research/casimir_mof/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `closure_2026-04-26` | > | 11/11 PASS | [[research/closure_2026-04-26/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `continuum_limit` | - `cg_results.txt` — 8/8 PASS — context: …oof)    CG-2 (LPA' kinetic stability):     - Already CLOSED: K_IR/K_UV = 1.000 | 8/8 PASS | [[research/continuum_limit/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `cosmo_tensions` | > | — | [[research/cosmo_tensions/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `external_review_2026-04-25` | **`P1_2_U_bounded_below_results.md`:** | — | [[research/external_review_2026-04-25/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `galaxy_scaling` | > | — | [[research/galaxy_scaling/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `hubble_tension` | > | — | [[research/hubble_tension/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `liquid_viscosity` | **`PLAN.md`:** | — | [[research/liquid_viscosity/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `mass_scaling_k4` | > | 9/9 PASS | [[research/mass_scaling_k4/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `metric_ansatz` | > | 11/11 PASS | [[research/metric_ansatz/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `muon_g_minus_2` | **`PLAN.md`:** | — | [[research/muon_g_minus_2/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `neutrino_msw` | **`PLAN.md`:** | — | [[research/neutrino_msw/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `op-alpha-fine-structure` | \| Plik \| Cycle \| Status \| Verdict \| Score \| Program \| | PASS / — | [[research/op-alpha-fine-structure/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `op-bh-alpha-threshold` | \| Plik \| Cycle \| Status \| Verdict \| Score \| Program \| | H0_REJECTED / — | [[research/op-bh-alpha-threshold/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `op-cosmology-closure` | - `m10_1_de.txt` — 6/6 PASS — context: …lsifiability   [PASS]  M10.1.6 T-Lambda closure consistency    Sub-cycle M10.1:  | 6/6 PASS | [[research/op-cosmology-closure/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `op-cross-sector-charge` | \| Plik \| Cycle \| Status \| Verdict \| Score \| Program \| | PASS / — | [[research/op-cross-sector-charge/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `op-delta1-g-tilde-derivation` | **`PLAN.md`:** | — | [[research/op-delta1-g-tilde-derivation/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `op-delta2-Nf-derivation` | **`PLAN.md`:** | — | [[research/op-delta2-Nf-derivation/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `op-eht` | - `op_eht_t1_pn_truncation.txt` — 4/4 PASS — context: …================================================================= | 4/4 PASS | [[research/op-eht/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `op-eht-A` | - `op_eht_a_T2_1pn_matching.txt` — 3/3 PASS — context: …================================================================ | 3/3 PASS | [[research/op-eht-A/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `op-eps-photon-ring` | \| Plik \| Cycle \| Status \| Verdict \| Score \| Program \| | PASS / — | [[research/op-eps-photon-ring/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `op-eta-wolfenstein` | \| Plik \| Cycle \| Status \| Verdict \| Score \| Program \| | PASS / — | [[research/op-eta-wolfenstein/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `op-eta2-denom-derivation` | \| Plik \| Cycle \| Status \| Verdict \| Score \| Program \| | PASS / — | [[research/op-eta2-denom-derivation/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `op-g0-r3-from-canonical-projection` | \| Plik \| Cycle \| Status \| Verdict \| Score \| Program \| | — / 2 | [[research/op-g0-r3-from-canonical-projection/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `op-gamma1-phi-eff-anchor-resolution` | **`PLAN.md`:** | — | [[research/op-gamma1-phi-eff-anchor-resolution/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `op-iota-charge-pmns-unification` | \| Plik \| Cycle \| Status \| Verdict \| Score \| Program \| | — / — | [[research/op-iota-charge-pmns-unification/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `op-kappa-mixing-numerator` | \| Plik \| Cycle \| Status \| Verdict \| Score \| Program \| | — / — | [[research/op-kappa-mixing-numerator/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `op-lambda1-e2-amplitude-emergence` | \| Plik \| Cycle \| Status \| Verdict \| Score \| Program \| | — / — | [[research/op-lambda1-e2-amplitude-emergence/FINDINGS.md]] | `mass_scaling_k4` | BROADCAST |
| 2026-05-03 | `op-mu-pmns-phase-hardening` | \| Plik \| Cycle \| Status \| Verdict \| Score \| Program \| | — / — | [[research/op-mu-pmns-phase-hardening/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `op-mu1-minimal-substrate-log-redefinition` | - `phase1_psi_soliton_and_mass_formula.txt` — 3/3 PASS — context: …ψ ≡ log g):     PASS   P1.3 (μ/e ratio gate):         | 3/3 PASS | [[research/op-mu1-minimal-substrate-log-redefinition/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `op-newton-momentum` | - `B6_m9x_sqrtg_rerun.txt` — 5/5 PASS — context: …e adds       [psi/(4-3psi)]^{3/2} ~ 1 + 6 eps weighting); structural M | 5/5 PASS | [[research/op-newton-momentum/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `op-nu-majorana-phase-mbb` | \| Plik \| Cycle \| Status \| Verdict \| Score \| Program \| | — / — | [[research/op-nu-majorana-phase-mbb/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `op-omega1-substrate-em-coupling` | \| Plik \| Cycle \| Status \| Verdict \| Score \| Program \| | — / — | [[research/op-omega1-substrate-em-coupling/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `op-omicron1-sigmamnu-cosmo` | \| Plik \| Cycle \| Status \| Verdict \| Score \| Program \| | — / — | [[research/op-omicron1-sigmamnu-cosmo/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `op-phase1-covariant` | \| Plik \| Cycle \| Status \| Verdict \| Score \| Program \| | 12/12 PASS | [[research/op-phase1-covariant/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `op-phase2-quantum-gravity` | \| Plik \| Cycle \| Status \| Verdict \| Score \| Program \| | 16/16 PASS | [[research/op-phase2-quantum-gravity/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `op-phase3-uv-completion` | \| Plik \| Cycle \| Status \| Verdict \| Score \| Program \| | 6/6 PASS | [[research/op-phase3-uv-completion/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `op-phi1-substrate-action-variational` | \| Plik \| Cycle \| Status \| Verdict \| Score \| Program \| | — / — | [[research/op-phi1-substrate-action-variational/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `op-pi1-bb0nu-nme-isotope` | \| Plik \| Cycle \| Status \| Verdict \| Score \| Program \| | — / — | [[research/op-pi1-bb0nu-nme-isotope/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `op-psi1-substrate-light-acceleration` | \| Plik \| Cycle \| Status \| Verdict \| Score \| Program \| | — / — | [[research/op-psi1-substrate-light-acceleration/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `op-quantum-closure` | - `m11_1_audit.txt` — 6/6 PASS — context: …] M11.1.4: PASS     [v] M11.1.5: PASS     [v] M11.1.6: PASS    >> M11.1 CLOSE | 6/6 PASS | [[research/op-quantum-closure/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `op-rho1-71Ge-cross-section` | \| Plik \| Cycle \| Status \| Verdict \| Score \| Program \| | — / — | [[research/op-rho1-71Ge-cross-section/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `op-sc-alpha-origin` | \| Plik \| Cycle \| Status \| Verdict \| Score \| Program \| | 4/4 PASS for H1 (structurally distinct) / — | [[research/op-sc-alpha-origin/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `op-sigma1-substrate-light-dispersion` | \| Plik \| Cycle \| Status \| Verdict \| Score \| Program \| | — / — | [[research/op-sigma1-substrate-light-dispersion/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `op-tau1-closure-overlap-coulomb` | \| Plik \| Cycle \| Status \| Verdict \| Score \| Program \| | — / — | [[research/op-tau1-closure-overlap-coulomb/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `op-tau2-substrate-time-coupling` | \| Plik \| Cycle \| Status \| Verdict \| Score \| Program \| | — / — | [[research/op-tau2-substrate-time-coupling/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `op-tau3-substrate-clock-acceleration` | \| Plik \| Cycle \| Status \| Verdict \| Score \| Program \| | — / — | [[research/op-tau3-substrate-clock-acceleration/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `op-theta-quark-koide` | \| Plik \| Cycle \| Status \| Verdict \| Score \| Program \| | PASS / — | [[research/op-theta-quark-koide/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `op-upsilon1-closure-cross-family` | \| Plik \| Cycle \| Status \| Verdict \| Score \| Program \| | — / — | [[research/op-upsilon1-closure-cross-family/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `op-uv-as-ngfp` | \| Plik \| Cycle \| Status \| Verdict \| Score \| Program \| | PASS / — | [[research/op-uv-as-ngfp/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `op-uv3-phi0-renormalization` | \| Plik \| Cycle \| Status \| Verdict \| Score \| Program \| | PASS / 5 | [[research/op-uv3-phi0-renormalization/FINDINGS.md]] | `op-delta1-g-tilde-derivation`, `op-gamma1-phi-eff-anchor-resolution` | BROADCAST |
| 2026-05-03 | `op-xi-photon-ring` | \| Plik \| Cycle \| Status \| Verdict \| Score \| Program \| | PASS / — | [[research/op-xi-photon-ring/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `op-xi2-sterile-nu-5sector` | \| Plik \| Cycle \| Status \| Verdict \| Score \| Program \| | — / — | [[research/op-xi2-sterile-nu-5sector/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `op-zeta-mass-spectrum` | \| Plik \| Cycle \| Status \| Verdict \| Score \| Program \| | PASS / — | [[research/op-zeta-mass-spectrum/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `op7` | - `eht_photon_ring_m911.txt` — 5/6 PASS — context: …==================================================================== | 5/6 PASS | [[research/op7/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `particle_sector_closure` | **`README.md`:** | — | [[research/particle_sector_closure/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `qm_born_rule` | > | — | [[research/qm_born_rule/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `qm_decoherence` | > | 8/8 PASS | [[research/qm_decoherence/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `qm_entanglement` | > | 3/4 PASS | [[research/qm_entanglement/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `qm_foundations` | > | 22/25 PASS | [[research/qm_foundations/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `qm_measurement` | > | 8/9 PASS | [[research/qm_measurement/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `qm_spin` | > | 7/7 PASS | [[research/qm_spin/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `qm_statistics` | > | 8/8 PASS | [[research/qm_statistics/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `qm_superposition` | > | 7/7 PASS | [[research/qm_superposition/FINDINGS.md]] | `atomic_shells_closure` | BROADCAST |
| 2026-05-03 | `s8_tension` | > | — | [[research/s8_tension/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `superconductivity_closure` | - `ps32_results.txt` — 0/5 PASS — context: …======================================   Szlachetne (Cu,Ag,Au,Pt,Pd)         | 0/5 PASS | [[research/superconductivity_closure/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `thermal_transport_molecular` | **`PLAN.md`:** | — | [[research/thermal_transport_molecular/FINDINGS.md]] | (orphan, no matches) | BROADCAST |
| 2026-05-03 | `uv_completion` | > | — | [[research/uv_completion/FINDINGS.md]] | (orphan, no matches) | BROADCAST |

## 3. Skanowania (consumers: [] — pending Sesja 6 matching)

All 70 broadcasts mają consumers: [] do czasu Sesji 6
(CANDIDATE_BRIDGES matching). To NIE oznacza orfan — oznacza pending.

Po Sesji 6 broadcasts przejdą jeden z trzech stanów:
- (a) Match z `NEEDS.md` innego folderu → wpis w `CANDIDATE_BRIDGES.md` jako PROPOSED
- (b) Brak match (orfan) → tabela § 3 z notatką `# scanned, no consumers found YYYY-MM-DD`
- (c) Pre-existing `depends_on:` w innym folderze cytuje ten folder → CONSUMED

## 4. CONSUMED

(empty — czeka na Sesję 6 matching i acknowledgments)

## 5. STALE (≥ 90 dni bez akcji)

(empty)

## 6. Foldery wyłączone z broadcast

- 4 polluted-74394a8 (kwarantanna; AGENT_PROTOCOL §3.2)
- 12 folderów z `exports_findings: false` (FINDINGS.md zawiera
  `(Brak ekstrahowalnych findings)`):
  - `atom_from_soliton`
  - `atomic_shells_closure`
  - `cohesion_closure`
  - `desi_dark_energy`
  - `em_from_substrate`
  - `nbody`
  - `op-m92`
  - `op-uv-renormalizability-research`
  - `op1-op2-op4`
  - `op6`
  - `rho_normal_state_closure`
  - `why_n3`

Te foldery są kandydatami do (a) manual review FINDINGS w Sesji 6+,
(b) flag w Sesji 8 jako wymagające scaffoldingu.