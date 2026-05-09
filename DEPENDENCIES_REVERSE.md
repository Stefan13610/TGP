# Reverse Dependencies - TGP/TGP_v1

> Generated: 2026-05-09 by `tooling/build_deps_graph.py`
> Sources: \input{} + \ref{} + \cite{} + [[wikilink]]

## Coarse view

### <root>
- depended on by: audyt, meta, papers, research

### audyt
- depended on by: <root>, meta, research

### axioms
- depended on by: <root>, audyt, core, core/_meta_latex, core/formalizm, meta, papers_external, partial_proofs, research

### core
- depended on by: <root>, audyt, axioms, core/_meta_latex, core/formalizm, papers_external, partial_proofs, research

### core/_meta_latex
- depended on by: <root>, audyt, core, research

### core/formalizm
- depended on by: <root>, axioms, core, core/_meta_latex, papers_external, partial_proofs, research

### meta
- depended on by: <root>, audyt, research

### papers
- depended on by: audyt, research

### papers_external
- depended on by: <root>, partial_proofs, research

### partial_proofs
- depended on by: <root>, axioms, core, core/_meta_latex, core/formalizm, research

### research
- depended on by: <root>, audyt, core, meta, papers, partial_proofs

### tooling
- depended on by: research

## Fine view (per subfolder)

### <root>
- depended on by:
  - audyt (12x | wikilink:12)
  - audyt/L01_rho_operational (2x | wikilink:2)
  - audyt/L02_beta_gamma_semantics (6x | wikilink:6)
  - audyt/L07_zero_sum_axiom (1x | wikilink:1)
  - audyt/L08_kink_fermion_closure (3x | wikilink:3)
  - audyt/M01_status_creep (4x | wikilink:4)
  - audyt/S01_metric_four_forms (2x | wikilink:2)
  - audyt/S05_tensor_sector_singleField (6x | wikilink:6)
  - audyt/S07_M911_derivation (1x | wikilink:1)
  - audyt/T01_LIGO3G_falsifier (24x | wikilink:24)
  - meta (4x | wikilink:4)
  - meta/core (1x | wikilink:1)
  - meta/research (1x | wikilink:1)
  - meta/research/_examples (3x | wikilink:3)
  - meta/research/templates (1x | wikilink:1)
  - papers/M911_LIGO3G_paper (1x | wikilink:1)
  - research/_archive (1x | wikilink:1)
  - research/_sandbox (1x | wikilink:1)
  - research/atom_from_soliton (1x | wikilink:1)
  - research/atomic_shells_closure (1x | wikilink:1)
  - research/audyt_cosmology_drift_2026-05-03 (1x | wikilink:1)
  - research/casimir_mof (1x | wikilink:1)
  - research/closure_2026-04-26 (2x | wikilink:2)
  - research/closure_2026-04-26/Lambda_from_Phi0 (2x | wikilink:2)
  - research/closure_2026-04-26/alpha_psi_threshold (1x | wikilink:1)
  - research/closure_2026-04-26/f_psi_principle (2x | wikilink:2)
  - research/closure_2026-04-26/sigma_ab_pathB (3x | wikilink:3)
  - research/cohesion_closure (1x | wikilink:1)
  - research/em_from_substrate (1x | wikilink:1)
  - research/external_review_2026-04-25 (1x | wikilink:1)
  - research/liquid_viscosity (1x | wikilink:1)
  - research/muon_g_minus_2 (1x | wikilink:1)
  - research/neutrino_msw (1x | wikilink:1)
  - research/op-CORE-CLEANUP-B-2026-05-04 (2x | wikilink:2)
  - research/op-FRW-radiation-era-varying-c-2026-05-06 (2x | wikilink:2)
  - research/op-GWTC3-reanalysis (4x | wikilink:4)
  - research/op-L01-rho-stress-energy-bridge-2026-05-04 (1x | wikilink:1)
  - research/op-L03-spectral-stability-2026-05-06 (2x | wikilink:2)
  - research/op-LIGO-3G-deviation (5x | wikilink:5)
  - research/op-M03-balance-sheet-retrofit-2026-05-06 (5x | wikilink:5)
  - research/op-Phi-vacuum-scale-2026-05-09 (2x | wikilink:2)
  - research/op-S07-alternative-f-psi-derivation-2026-05-09 (4x | wikilink:4)
  - research/op-SPIN-SU2-substrate-derivation-2026-05-08 (1x | wikilink:1)
  - research/op-alpha-fine-structure (1x | wikilink:1)
  - research/op-bh-alpha-threshold (3x | wikilink:3)
  - research/op-chi1-newton-constant-derivation (10x | wikilink:10)
  - research/op-cosmology-closure (3x | wikilink:3)
  - research/op-cross-sector-charge (2x | wikilink:2)
  - research/op-eht (1x | wikilink:1)
  - research/op-eht-A (1x | wikilink:1)
  - research/op-emergent-metric-from-interaction-2026-05-09 (3x | wikilink:3)
  - research/op-eps-photon-ring (1x | wikilink:1)
  - research/op-eta-wolfenstein (1x | wikilink:1)
  - research/op-eta2-denom-derivation (2x | wikilink:2)
  - research/op-g0-r3-from-canonical-projection (1x | wikilink:1)
  - research/op-iota-charge-pmns-unification (7x | wikilink:7)
  - research/op-kappa-mixing-numerator (6x | wikilink:6)
  - research/op-m92 (1x | wikilink:1)
  - research/op-mu-pmns-phase-hardening (11x | wikilink:11)
  - research/op-newton-momentum (9x | wikilink:9)
  - research/op-nu-majorana-phase-mbb (7x | wikilink:7)
  - research/op-omega1-substrate-em-coupling (5x | wikilink:5)
  - research/op-omega2-axion-coupling-lock (7x | wikilink:7)
  - research/op-omega3-axion-decay-constant (4x | wikilink:4)
  - research/op-omicron1-sigmamnu-cosmo (6x | wikilink:6)
  - research/op-omicron2-phi-mean-shift-cosmo (1x | wikilink:1)
  - research/op-phase1-covariant (1x | wikilink:1)
  - research/op-phase2-quantum-gravity (1x | wikilink:1)
  - research/op-phase3-uv-completion (1x | wikilink:1)
  - research/op-phi1-substrate-action-variational (5x | wikilink:5)
  - research/op-pi1-bb0nu-nme-isotope (6x | wikilink:6)
  - research/op-ppE-mapping (7x | wikilink:7)
  - research/op-ppE-mapping/TGP/TGP_v1/audyt/T01_LIGO3G_falsifier (3x | wikilink:3)
  - research/op-ppE-mapping/TGP/TGP_v1/papers/M911_LIGO3G_paper (1x | wikilink:1)
  - research/op-ppE-mapping/TGP/TGP_v1/research/op-LIGO-3G-deviation (5x | wikilink:5)
  - research/op-ppE-mapping/TGP/TGP_v1/research/op-ppE-mapping (3x | wikilink:3)
  - research/op-psi1-substrate-light-acceleration (2x | wikilink:2)
  - research/op-quantum-closure (2x | wikilink:2)
  - research/op-rho1-71Ge-cross-section (10x | wikilink:10)
  - research/op-sc-alpha-origin (8x | wikilink:8)
  - research/op-sigma1-substrate-light-dispersion (5x | wikilink:5)
  - research/op-tau1-closure-overlap-coulomb (6x | wikilink:6)
  - research/op-tau2-substrate-time-coupling (3x | wikilink:3)
  - research/op-tau3-substrate-clock-acceleration (4x | wikilink:4)
  - research/op-theta-quark-koide (1x | wikilink:1)
  - research/op-upsilon1-closure-cross-family (6x | wikilink:6)
  - research/op-uv-as-ngfp (2x | wikilink:2)
  - research/op-uv2-mtgp-absolute-scale (6x | wikilink:6)
  - research/op-void-flat-modes-h0-2026-05-06 (2x | wikilink:2)
  - research/op-xi-photon-ring (2x | wikilink:2)
  - research/op-xi2-sterile-nu-5sector (11x | wikilink:11)
  - research/op-zeta-mass-spectrum (1x | wikilink:1)
  - research/op7 (2x | wikilink:2)
  - research/rho_normal_state_closure (1x | wikilink:1)
  - research/thermal_transport_molecular (1x | wikilink:1)
  - research/why_n3 (1x | wikilink:1)

### audyt
- depended on by:
  - <root> (4x | wikilink:4)
  - audyt/D01_drifting_numbers (5x | wikilink:5)
  - audyt/L01_rho_operational (7x | wikilink:7)
  - audyt/L02_beta_gamma_semantics (3x | wikilink:3)
  - audyt/L03_K_phi_stability (5x | wikilink:5)
  - audyt/L04_ODE_dualism_alpha (3x | wikilink:3)
  - audyt/L05_mass_exponent_drift (3x | wikilink:3)
  - audyt/L06_axion_mass_locked (3x | wikilink:3)
  - audyt/L07_zero_sum_axiom (5x | wikilink:5)
  - audyt/L08_kink_fermion_closure (5x | wikilink:5)
  - audyt/M01_status_creep (3x | wikilink:3)
  - audyt/M02_ledger_pollution (3x | wikilink:3)
  - audyt/M03_balance_sheet_missing (11x | wikilink:11)
  - audyt/S01_metric_four_forms (4x | wikilink:4)
  - audyt/S02_volume_element_M9 (3x | wikilink:3)
  - audyt/S03_beta_PPN_convention (3x | wikilink:3)
  - audyt/S04_metric_coupling_axiom (3x | wikilink:3)
  - audyt/S05_tensor_sector_singleField (3x | wikilink:3)
  - audyt/S06_circular_anchors (3x | wikilink:3)
  - audyt/S07_M911_derivation (5x | wikilink:5)
  - audyt/T01_LIGO3G_falsifier (13x | wikilink:13)
  - research/op-CORE-CLEANUP-B-2026-05-04 (4x | wikilink:4)
  - research/op-D01-anchor-lock-2026-05-06 (2x | wikilink:2)
  - research/op-FRW-radiation-era-varying-c-2026-05-06 (9x | wikilink:9)
  - research/op-L03-spectral-stability-2026-05-06 (1x | wikilink:1)
  - research/op-M03-balance-sheet-retrofit-2026-05-06 (2x | wikilink:2)
  - research/op-Phi-decomposition-photon-2026-05-07 (2x | wikilink:2)
  - research/op-S07-alternative-f-psi-derivation-2026-05-09 (1x | wikilink:1)

### audyt/D01_drifting_numbers
- depended on by:
  - <root> (5x | wikilink:5)
  - audyt (73x | wikilink:73)
  - audyt/L01_rho_operational (18x | wikilink:18)
  - audyt/L02_beta_gamma_semantics (5x | wikilink:5)
  - audyt/L03_K_phi_stability (8x | wikilink:8)
  - audyt/L04_ODE_dualism_alpha (12x | wikilink:12)
  - audyt/L05_mass_exponent_drift (1x | wikilink:1)
  - audyt/L06_axion_mass_locked (1x | wikilink:1)
  - audyt/L07_zero_sum_axiom (3x | wikilink:3)
  - audyt/L08_kink_fermion_closure (7x | wikilink:7)
  - audyt/M01_status_creep (1x | wikilink:1)
  - audyt/M02_ledger_pollution (1x | wikilink:1)
  - audyt/M03_balance_sheet_missing (17x | wikilink:17)
  - audyt/S01_metric_four_forms (11x | wikilink:11)
  - audyt/S02_volume_element_M9 (7x | wikilink:7)
  - audyt/S03_beta_PPN_convention (7x | wikilink:7)
  - audyt/S04_metric_coupling_axiom (7x | wikilink:7)
  - audyt/S05_tensor_sector_singleField (9x | wikilink:9)
  - audyt/S06_circular_anchors (7x | wikilink:7)
  - audyt/S07_M911_derivation (7x | wikilink:7)
  - audyt/T01_LIGO3G_falsifier (87x | wikilink:87)
  - meta (3x | wikilink:3)
  - meta/core (1x | wikilink:1)
  - meta/research/_examples (6x | wikilink:6)
  - meta/research/templates (5x | wikilink:5)
  - research/_archive (3x | wikilink:3)
  - research/_sandbox (1x | wikilink:1)
  - research/atom_from_soliton (5x | wikilink:5)
  - research/atomic_shells_closure (8x | wikilink:8)
  - research/audyt_cosmology_drift_2026-05-03 (4x | wikilink:4)
  - research/brannen_sqrt2 (5x | wikilink:5)
  - research/cabibbo_correction (5x | wikilink:5)
  - research/casimir_mof (5x | wikilink:5)
  - research/closure_2026-04-26 (8x | wikilink:8)
  - research/closure_2026-04-26/Lambda_from_Phi0 (2x | wikilink:2)
  - research/cohesion_closure (5x | wikilink:5)
  - research/continuum_limit (5x | wikilink:5)
  - research/cosmo_tensions (6x | wikilink:6)
  - research/desi_dark_energy (5x | wikilink:5)
  - research/em_from_substrate (5x | wikilink:5)
  - research/external_review_2026-04-25 (6x | wikilink:6)
  - research/galaxy_scaling (5x | wikilink:5)
  - research/hubble_tension (6x | wikilink:6)
  - research/liquid_viscosity (5x | wikilink:5)
  - research/mass_scaling_k4 (6x | wikilink:6)
  - research/metric_ansatz (5x | wikilink:5)
  - research/muon_g_minus_2 (5x | wikilink:5)
  - research/nbody (5x | wikilink:5)
  - research/nbody/docs (1x | wikilink:1)
  - research/nbody/paper (1x | wikilink:1)
  - research/neutrino_msw (5x | wikilink:5)
  - research/op-CORE-CLEANUP-B-2026-05-04 (10x | wikilink:10)
  - research/op-D01-anchor-lock-2026-05-06 (17x | wikilink:17)
  - research/op-FRW-radiation-era-varying-c-2026-05-06 (24x | wikilink:24)
  - research/op-GWTC3-reanalysis (12x | wikilink:12)
  - research/op-L01-rho-stress-energy-bridge-2026-05-04 (15x | wikilink:15)
  - research/op-L03-spectral-stability-2026-05-06 (24x | wikilink:24)
  - research/op-L04-ODE-canonicalization-2026-05-04 (18x | wikilink:18)
  - research/op-LIGO-3G-deviation (12x | wikilink:12)
  - research/op-M03-balance-sheet-retrofit-2026-05-06 (186x | wikilink:186)
  - research/op-MAG-Lorentz-A-mu-coupling-2026-05-09 (8x | wikilink:8)
  - research/op-MAG-Phase5-V-reference-clarification-2026-05-09 (8x | wikilink:8)
  - research/op-MAG-anomalous-moment-2026-05-09 (6x | wikilink:6)
  - research/op-MAG-resonance-formalization-2026-05-09 (29x | wikilink:29)
  - research/op-Phase5-MAG-erratum-2026-05-09 (3x | wikilink:3)
  - research/op-Phi-decomposition-photon-2026-05-07 (29x | wikilink:29)
  - research/op-Phi-vacuum-scale-2026-05-09 (26x | wikilink:26)
  - research/op-Phi0-spatial-variation-predictions-2026-05-09 (3x | wikilink:3)
  - research/op-S07-alternative-f-psi-derivation-2026-05-09 (11x | wikilink:11)
  - research/op-SPIN-SU2-substrate-derivation-2026-05-08 (44x | wikilink:44)
  - research/op-V-canonical-consistency-audit-2026-05-09 (9x | wikilink:9)
  - research/op-alpha-fine-structure (5x | wikilink:5)
  - research/op-bh-alpha-threshold (5x | wikilink:5)
  - research/op-chi1-newton-constant-derivation (3x | wikilink:3)
  - research/op-cosmology-closure (8x | wikilink:8)
  - research/op-cross-sector-charge (5x | wikilink:5)
  - research/op-delta1-g-tilde-derivation (9x | wikilink:9)
  - research/op-delta2-Nf-derivation (10x | wikilink:10)
  - research/op-dual-V-structure-clarification-2026-05-09 (8x | wikilink:8)
  - research/op-eht (5x | wikilink:5)
  - research/op-eht-A (5x | wikilink:5)
  - research/op-emergent-metric-from-interaction-2026-05-09 (11x | wikilink:11)
  - research/op-eps-photon-ring (5x | wikilink:5)
  - research/op-eta-wolfenstein (5x | wikilink:5)
  - research/op-eta2-denom-derivation (5x | wikilink:5)
  - research/op-g0-r3-from-canonical-projection (16x | wikilink:16)
  - research/op-gamma1-phi-eff-anchor-resolution (8x | wikilink:8)
  - research/op-iota-charge-pmns-unification (5x | wikilink:5)
  - research/op-kappa-mixing-numerator (5x | wikilink:5)
  - research/op-lambda1-e2-amplitude-emergence (13x | wikilink:13)
  - research/op-m92 (5x | wikilink:5)
  - research/op-mu-pmns-phase-hardening (5x | wikilink:5)
  - research/op-mu1-minimal-substrate-log-redefinition (9x | wikilink:9)
  - research/op-newton-momentum (5x | wikilink:5)
  - research/op-nu-majorana-phase-mbb (5x | wikilink:5)
  - research/op-omega1-substrate-em-coupling (5x | wikilink:5)
  - research/op-omega2-axion-coupling-lock (3x | wikilink:3)
  - research/op-omega3-axion-decay-constant (3x | wikilink:3)
  - research/op-omicron1-sigmamnu-cosmo (5x | wikilink:5)
  - research/op-omicron2-phi-mean-shift-cosmo (9x | wikilink:9)
  - research/op-phase1-covariant (5x | wikilink:5)
  - research/op-phase2-quantum-gravity (5x | wikilink:5)
  - research/op-phase3-uv-completion (7x | wikilink:7)
  - research/op-phi1-substrate-action-variational (5x | wikilink:5)
  - research/op-pi1-bb0nu-nme-isotope (5x | wikilink:5)
  - research/op-ppE-mapping (17x | wikilink:17)
  - research/op-ppE-mapping/TGP/TGP_v1/audyt/T01_LIGO3G_falsifier (7x | wikilink:7)
  - research/op-ppE-mapping/TGP/TGP_v1/research/op-LIGO-3G-deviation (12x | wikilink:12)
  - research/op-ppE-mapping/TGP/TGP_v1/research/op-ppE-mapping (10x | wikilink:10)
  - research/op-psi1-substrate-light-acceleration (5x | wikilink:5)
  - research/op-quantum-closure (6x | wikilink:6)
  - research/op-rho1-71Ge-cross-section (5x | wikilink:5)
  - research/op-sc-alpha-origin (5x | wikilink:5)
  - research/op-sigma1-substrate-light-dispersion (5x | wikilink:5)
  - research/op-tau1-closure-overlap-coulomb (5x | wikilink:5)
  - research/op-tau2-substrate-time-coupling (5x | wikilink:5)
  - research/op-tau3-substrate-clock-acceleration (5x | wikilink:5)
  - research/op-theta-quark-koide (5x | wikilink:5)
  - research/op-upsilon1-closure-cross-family (5x | wikilink:5)
  - research/op-uv-as-ngfp (5x | wikilink:5)
  - research/op-uv-renormalizability-research (8x | wikilink:8)
  - research/op-uv2-mtgp-absolute-scale (3x | wikilink:3)
  - research/op-uv3-phi0-renormalization (9x | wikilink:9)
  - research/op-void-flat-modes-h0-2026-05-06 (17x | wikilink:17)
  - research/op-xi-photon-ring (5x | wikilink:5)
  - research/op-xi2-sterile-nu-5sector (5x | wikilink:5)
  - research/op-zeta-mass-spectrum (5x | wikilink:5)
  - research/op1-op2-op4 (5x | wikilink:5)
  - research/op6 (5x | wikilink:5)
  - research/op7 (5x | wikilink:5)
  - research/particle_sector_closure (10x | wikilink:10)
  - research/qm_born_rule (5x | wikilink:5)
  - research/qm_decoherence (5x | wikilink:5)
  - research/qm_entanglement (5x | wikilink:5)
  - research/qm_foundations (5x | wikilink:5)
  - research/qm_measurement (5x | wikilink:5)
  - research/qm_spin (5x | wikilink:5)
  - research/qm_statistics (5x | wikilink:5)
  - research/qm_superposition (5x | wikilink:5)
  - research/rho_normal_state_closure (5x | wikilink:5)
  - research/s8_tension (5x | wikilink:5)
  - research/superconductivity_closure (5x | wikilink:5)
  - research/thermal_transport_molecular (5x | wikilink:5)
  - research/uv_completion (5x | wikilink:5)
  - research/why_n3 (5x | wikilink:5)

### audyt/L01_rho_operational
- depended on by:
  - audyt (13x | wikilink:13)
  - audyt/S02_volume_element_M9 (3x | wikilink:3)
  - audyt/S03_beta_PPN_convention (3x | wikilink:3)
  - research/op-FRW-radiation-era-varying-c-2026-05-06 (9x | wikilink:9)
  - research/op-L01-rho-stress-energy-bridge-2026-05-04 (2x | wikilink:2)
  - research/op-L04-ODE-canonicalization-2026-05-04 (2x | wikilink:2)
  - research/op-Phi-decomposition-photon-2026-05-07 (4x | wikilink:4)

### audyt/L02_beta_gamma_semantics
- depended on by: --

### audyt/L03_K_phi_stability
- depended on by: --

### audyt/L04_ODE_dualism_alpha
- depended on by: --

### audyt/L05_mass_exponent_drift
- depended on by: --

### audyt/L06_axion_mass_locked
- depended on by: --

### audyt/L07_zero_sum_axiom
- depended on by: --

### audyt/L08_kink_fermion_closure
- depended on by: --

### audyt/M01_status_creep
- depended on by: --

### audyt/M02_ledger_pollution
- depended on by: --

### audyt/M03_balance_sheet_missing
- depended on by:
  - <root> (1x | wikilink:1)
  - audyt (1x | wikilink:1)
  - research/op-M03-balance-sheet-retrofit-2026-05-06 (3x | wikilink:3)

### audyt/S01_metric_four_forms
- depended on by: --

### audyt/S02_volume_element_M9
- depended on by: --

### audyt/S03_beta_PPN_convention
- depended on by: --

### audyt/S04_metric_coupling_axiom
- depended on by: --

### audyt/S05_tensor_sector_singleField
- depended on by: --

### audyt/S06_circular_anchors
- depended on by: --

### audyt/S07_M911_derivation
- depended on by: --

### audyt/T01_LIGO3G_falsifier
- depended on by:
  - <root> (5x | wikilink:5)
  - audyt (2x | wikilink:2)
  - audyt/D01_drifting_numbers (1x | wikilink:1)
  - audyt/L01_rho_operational (1x | wikilink:1)
  - audyt/L03_K_phi_stability (1x | wikilink:1)
  - audyt/L04_ODE_dualism_alpha (1x | wikilink:1)
  - meta/core (1x | wikilink:1)
  - meta/research (70x | wikilink:70)
  - meta/research/templates (1x | wikilink:1)
  - research/atom_from_soliton (1x | wikilink:1)
  - research/atomic_shells_closure (1x | wikilink:1)
  - research/brannen_sqrt2 (1x | wikilink:1)
  - research/cabibbo_correction (1x | wikilink:1)
  - research/casimir_mof (1x | wikilink:1)
  - research/closure_2026-04-26 (1x | wikilink:1)
  - research/cohesion_closure (1x | wikilink:1)
  - research/continuum_limit (1x | wikilink:1)
  - research/cosmo_tensions (1x | wikilink:1)
  - research/desi_dark_energy (1x | wikilink:1)
  - research/em_from_substrate (1x | wikilink:1)
  - research/external_review_2026-04-25 (1x | wikilink:1)
  - research/galaxy_scaling (1x | wikilink:1)
  - research/hubble_tension (1x | wikilink:1)
  - research/liquid_viscosity (1x | wikilink:1)
  - research/mass_scaling_k4 (1x | wikilink:1)
  - research/metric_ansatz (1x | wikilink:1)
  - research/muon_g_minus_2 (1x | wikilink:1)
  - research/nbody (1x | wikilink:1)
  - research/neutrino_msw (1x | wikilink:1)
  - research/op-CORE-CLEANUP-B-2026-05-04 (1x | wikilink:1)
  - research/op-D01-anchor-lock-2026-05-06 (2x | wikilink:2)
  - research/op-FRW-radiation-era-varying-c-2026-05-06 (2x | wikilink:2)
  - research/op-GWTC3-reanalysis (13x | wikilink:13)
  - research/op-L01-rho-stress-energy-bridge-2026-05-04 (2x | wikilink:2)
  - research/op-L03-spectral-stability-2026-05-06 (3x | wikilink:3)
  - research/op-L04-ODE-canonicalization-2026-05-04 (2x | wikilink:2)
  - research/op-LIGO-3G-deviation (14x | wikilink:14)
  - research/op-M03-balance-sheet-retrofit-2026-05-06 (3x | wikilink:3)
  - research/op-Phi-decomposition-photon-2026-05-07 (7x | wikilink:7)
  - research/op-alpha-fine-structure (1x | wikilink:1)
  - research/op-bh-alpha-threshold (1x | wikilink:1)
  - research/op-chi1-newton-constant-derivation (1x | wikilink:1)
  - research/op-cosmology-closure (1x | wikilink:1)
  - research/op-cross-sector-charge (1x | wikilink:1)
  - research/op-delta1-g-tilde-derivation (1x | wikilink:1)
  - research/op-delta2-Nf-derivation (1x | wikilink:1)
  - research/op-eht (1x | wikilink:1)
  - research/op-eht-A (1x | wikilink:1)
  - research/op-eps-photon-ring (1x | wikilink:1)
  - research/op-eta-wolfenstein (1x | wikilink:1)
  - research/op-eta2-denom-derivation (1x | wikilink:1)
  - research/op-g0-r3-from-canonical-projection (1x | wikilink:1)
  - research/op-gamma1-phi-eff-anchor-resolution (1x | wikilink:1)
  - research/op-iota-charge-pmns-unification (1x | wikilink:1)
  - research/op-kappa-mixing-numerator (1x | wikilink:1)
  - research/op-lambda1-e2-amplitude-emergence (1x | wikilink:1)
  - research/op-m92 (1x | wikilink:1)
  - research/op-mu-pmns-phase-hardening (1x | wikilink:1)
  - research/op-mu1-minimal-substrate-log-redefinition (1x | wikilink:1)
  - research/op-newton-momentum (1x | wikilink:1)
  - research/op-nu-majorana-phase-mbb (1x | wikilink:1)
  - research/op-omega1-substrate-em-coupling (1x | wikilink:1)
  - research/op-omega2-axion-coupling-lock (1x | wikilink:1)
  - research/op-omega3-axion-decay-constant (1x | wikilink:1)
  - research/op-omicron1-sigmamnu-cosmo (1x | wikilink:1)
  - research/op-phase1-covariant (1x | wikilink:1)
  - research/op-phase2-quantum-gravity (1x | wikilink:1)
  - research/op-phase3-uv-completion (1x | wikilink:1)
  - research/op-phi1-substrate-action-variational (1x | wikilink:1)
  - research/op-pi1-bb0nu-nme-isotope (1x | wikilink:1)
  - research/op-ppE-mapping (34x | wikilink:34)
  - research/op-ppE-mapping/TGP/TGP_v1/audyt/T01_LIGO3G_falsifier (11x | wikilink:11)
  - research/op-ppE-mapping/TGP/TGP_v1/research/op-LIGO-3G-deviation (14x | wikilink:14)
  - research/op-ppE-mapping/TGP/TGP_v1/research/op-ppE-mapping (13x | wikilink:13)
  - research/op-psi1-substrate-light-acceleration (1x | wikilink:1)
  - research/op-quantum-closure (1x | wikilink:1)
  - research/op-rho1-71Ge-cross-section (1x | wikilink:1)
  - research/op-sc-alpha-origin (1x | wikilink:1)
  - research/op-sigma1-substrate-light-dispersion (1x | wikilink:1)
  - research/op-tau1-closure-overlap-coulomb (1x | wikilink:1)
  - research/op-tau2-substrate-time-coupling (1x | wikilink:1)
  - research/op-tau3-substrate-clock-acceleration (1x | wikilink:1)
  - research/op-theta-quark-koide (1x | wikilink:1)
  - research/op-upsilon1-closure-cross-family (1x | wikilink:1)
  - research/op-uv-as-ngfp (1x | wikilink:1)
  - research/op-uv-renormalizability-research (1x | wikilink:1)
  - research/op-uv2-mtgp-absolute-scale (1x | wikilink:1)
  - research/op-uv3-phi0-renormalization (1x | wikilink:1)
  - research/op-void-flat-modes-h0-2026-05-06 (1x | wikilink:1)
  - research/op-xi-photon-ring (1x | wikilink:1)
  - research/op-xi2-sterile-nu-5sector (1x | wikilink:1)
  - research/op-zeta-mass-spectrum (1x | wikilink:1)
  - research/op1-op2-op4 (1x | wikilink:1)
  - research/op6 (1x | wikilink:1)
  - research/op7 (1x | wikilink:1)
  - research/particle_sector_closure (1x | wikilink:1)
  - research/qm_born_rule (1x | wikilink:1)
  - research/qm_decoherence (1x | wikilink:1)
  - research/qm_entanglement (1x | wikilink:1)
  - research/qm_foundations (1x | wikilink:1)
  - research/qm_measurement (1x | wikilink:1)
  - research/qm_spin (1x | wikilink:1)
  - research/qm_statistics (1x | wikilink:1)
  - research/qm_superposition (1x | wikilink:1)
  - research/rho_normal_state_closure (1x | wikilink:1)
  - research/s8_tension (1x | wikilink:1)
  - research/superconductivity_closure (1x | wikilink:1)
  - research/thermal_transport_molecular (1x | wikilink:1)
  - research/uv_completion (1x | wikilink:1)
  - research/why_n3 (1x | wikilink:1)

### axioms
- depended on by: --

### axioms/notacja
- depended on by:
  - <root> (2x | input:2)
  - audyt (2x | wikilink:2)
  - audyt/L02_beta_gamma_semantics (3x | wikilink:3)
  - audyt/L03_K_phi_stability (1x | wikilink:1)
  - axioms/roznica_N0 (2x | ref:2)
  - core/_meta_latex (2x | ref:2)
  - meta/research (2x | wikilink:2)
  - partial_proofs/bh_ringdown (2x | ref:2)
  - partial_proofs/hierarchia_mas (2x | ref:2)
  - research/op-L03-spectral-stability-2026-05-06 (1x | wikilink:1)

### axioms/roznica_N0
- depended on by:
  - <root> (1x | input:1)
  - core/sek01_ontologia (5x | ref:5)

### axioms/substrat
- depended on by:
  - <root> (1x | input:1)
  - axioms/roznica_N0 (2x | ref:2)
  - core/_meta_latex (1x | ref:1)
  - core/formalizm (13x | ref:13)
  - core/sek00_summary (2x | ref:2)
  - core/sek01_ontologia (3x | ref:3)
  - core/sek06_czarne_dziury (1x | ref:1)
  - core/sek08_formalizm (26x | ref:26)
  - core/sek08b_ghost_resolution (5x | ref:5)
  - core/sek08c_metryka_z_substratu (2x | ref:2)
  - core/sek09_cechowanie (1x | ref:1)
  - core/sek10_N0_wyprowadzenie (3x | ref:3)
  - papers_external/tgp_core_paper (1x | ref:1)
  - partial_proofs/bh_ringdown (2x | ref:2)
  - partial_proofs/hierarchia_mas (1x | ref:1)
  - partial_proofs/superconductivity (1x | ref:1)
  - partial_proofs/wielki_wybuch (8x | ref:8)

### core
- depended on by: --

### core/_meta_latex
- depended on by:
  - <root> (2x | input:2)
  - audyt/M03_balance_sheet_missing (1x | wikilink:1)
  - core/sek00_summary (3x | input:1, ref:2)
  - core/sek08_formalizm (2x | ref:2)
  - research/op-M03-balance-sheet-retrofit-2026-05-06 (1x | wikilink:1)
  - research/op-uv3-phi0-renormalization (1x | wikilink:1)

### core/formalizm
- depended on by:
  - <root> (11x | input:11)
  - axioms/roznica_N0 (1x | ref:1)
  - axioms/substrat (4x | ref:4)
  - core/_meta_latex (9x | ref:9)
  - core/sek00_summary (3x | ref:3)
  - core/sek03_rezimy (1x | ref:1)
  - core/sek04_stale (1x | ref:1)
  - core/sek07_predykcje (14x | ref:14)
  - core/sek07a_wymiar_wzmocniony (2x | ref:2)
  - core/sek08_formalizm (10x | ref:10)
  - core/sek08a_akcja_zunifikowana (1x | ref:1)
  - core/sek08c_metryka_z_substratu (9x | ref:9)
  - core/sek09_cechowanie (1x | ref:1)
  - core/sek10_N0_wyprowadzenie (6x | ref:6)
  - papers_external/tgp_english_summary (1x | ref:1)
  - partial_proofs/bh_ringdown (1x | ref:1)
  - partial_proofs/defect_hierarchy (2x | ref:2)
  - partial_proofs/fermion_from_soliton (2x | wikilink:2)
  - partial_proofs/hierarchia_mas (12x | ref:12)
  - partial_proofs/most_gamma_phi (3x | ref:3)
  - partial_proofs/nuclear_from_soliton (2x | wikilink:2)
  - partial_proofs/particle_sector (7x | ref:7)
  - partial_proofs/superconductivity (3x | ref:3)
  - partial_proofs/wielki_wybuch (2x | ref:2)
  - research/atom_from_soliton (1x | wikilink:1)
  - research/em_from_substrate (1x | wikilink:1)
  - research/nbody/paper (3x | ref:3)
  - research/op-CORE-CLEANUP-B-2026-05-04 (2x | wikilink:2)
  - research/op-M03-balance-sheet-retrofit-2026-05-06 (1x | wikilink:1)
  - research/op-uv3-phi0-renormalization (6x | wikilink:6)
  - research/why_n3 (1x | wikilink:1)

### core/sek00_summary
- depended on by:
  - <root> (1x | input:1)
  - audyt/M03_balance_sheet_missing (1x | wikilink:1)
  - core/_meta_latex (1x | ref:1)
  - core/formalizm (1x | ref:1)
  - core/sek09_cechowanie (1x | ref:1)
  - research/op-M03-balance-sheet-retrofit-2026-05-06 (1x | wikilink:1)
  - research/op-uv3-phi0-renormalization (3x | wikilink:3)

### core/sek01_ontologia
- depended on by:
  - <root> (1x | input:1)
  - audyt/L07_zero_sum_axiom (3x | wikilink:3)
  - axioms/notacja (9x | ref:9)
  - axioms/roznica_N0 (5x | ref:5)
  - axioms/substrat (9x | ref:9)
  - core/formalizm (8x | ref:8)
  - core/sek00_summary (4x | ref:4)
  - core/sek02_pole (2x | ref:2)
  - core/sek04_stale (1x | ref:1)
  - core/sek05_ciemna_energia (1x | ref:1)
  - core/sek06_czarne_dziury (2x | ref:2)
  - core/sek07_predykcje (2x | ref:2)
  - core/sek07a_wymiar_wzmocniony (1x | ref:1)
  - core/sek08_formalizm (52x | ref:52)
  - core/sek09_cechowanie (4x | ref:4)
  - core/sek10_N0_wyprowadzenie (1x | ref:1)
  - papers_external/tgp_core_paper (7x | ref:7)
  - papers_external/tgp_english_summary (4x | ref:4)
  - partial_proofs/particle_sector (2x | ref:2)
  - partial_proofs/wielki_wybuch (2x | ref:2)
  - research/op-L03-spectral-stability-2026-05-06 (1x | wikilink:1)
  - research/op-omicron2-phi-mean-shift-cosmo (1x | wikilink:1)
  - research/op-void-flat-modes-h0-2026-05-06 (1x | wikilink:1)

### core/sek02_pole
- depended on by:
  - <root> (1x | input:1)
  - axioms/notacja (1x | ref:1)
  - core/formalizm (4x | ref:4)
  - core/sek00_summary (2x | ref:2)
  - core/sek01_ontologia (2x | ref:2)
  - core/sek03_rezimy (2x | ref:2)
  - core/sek08_formalizm (5x | ref:5)
  - core/sek08a_akcja_zunifikowana (4x | ref:4)
  - core/sek08b_ghost_resolution (1x | ref:1)
  - core/sek09_cechowanie (4x | ref:4)
  - papers_external/tgp_english_summary (5x | ref:5)
  - partial_proofs/defect_hierarchy (1x | ref:1)
  - partial_proofs/superconductivity (3x | ref:3)

### core/sek03_rezimy
- depended on by:
  - <root> (1x | input:1)
  - axioms/notacja (1x | ref:1)
  - core/_meta_latex (1x | ref:1)
  - core/formalizm (11x | ref:11)
  - core/sek02_pole (2x | ref:2)
  - core/sek07_predykcje (6x | ref:6)
  - core/sek07a_wymiar_wzmocniony (2x | ref:2)
  - core/sek08_formalizm (11x | ref:11)
  - core/sek09_cechowanie (3x | ref:3)
  - papers_external/tgp_english_summary (5x | ref:5)
  - partial_proofs/bounce_topology (1x | ref:1)
  - partial_proofs/hierarchia_mas (6x | ref:6)
  - partial_proofs/particle_sector (1x | ref:1)
  - partial_proofs/quark_sector (2x | ref:2)
  - partial_proofs/trojcialowe_nbody (1x | ref:1)
  - partial_proofs/wielki_wybuch (2x | ref:2)

### core/sek04_stale
- depended on by:
  - <root> (1x | input:1)
  - audyt/L01_rho_operational (2x | wikilink:2)
  - axioms/notacja (5x | ref:5)
  - core/_meta_latex (7x | ref:7)
  - core/formalizm (13x | ref:13)
  - core/sek00_summary (6x | ref:6)
  - core/sek01_ontologia (2x | ref:2)
  - core/sek02_pole (2x | ref:2)
  - core/sek06_czarne_dziury (12x | ref:12)
  - core/sek07_predykcje (12x | ref:12)
  - core/sek07a_wymiar_wzmocniony (1x | ref:1)
  - core/sek08_formalizm (24x | ref:24)
  - core/sek08c_metryka_z_substratu (3x | ref:3)
  - core/sek09_cechowanie (1x | ref:1)
  - papers_external/tgp_english_summary (7x | ref:7)
  - partial_proofs/bh_ringdown (3x | ref:3)
  - partial_proofs/trojcialowe_nbody (2x | ref:2)
  - partial_proofs/wielki_wybuch (2x | ref:2)
  - research/op-FRW-radiation-era-varying-c-2026-05-06 (3x | wikilink:3)
  - research/op-Phi-decomposition-photon-2026-05-07 (4x | wikilink:4)

### core/sek05_ciemna_energia
- depended on by:
  - <root> (1x | input:1)
  - audyt/L07_zero_sum_axiom (3x | wikilink:3)
  - axioms/notacja (1x | ref:1)
  - axioms/substrat (1x | ref:1)
  - core/_meta_latex (1x | ref:1)
  - core/formalizm (3x | ref:3)
  - core/sek00_summary (1x | ref:1)
  - core/sek02_pole (1x | ref:1)
  - core/sek07_predykcje (6x | ref:6)
  - core/sek08_formalizm (4x | ref:4)
  - papers_external/tgp_english_summary (2x | ref:2)
  - partial_proofs/wielki_wybuch (1x | ref:1)
  - research/op-FRW-radiation-era-varying-c-2026-05-06 (3x | wikilink:3)

### core/sek06_czarne_dziury
- depended on by:
  - <root> (1x | input:1)
  - axioms/notacja (2x | ref:2)
  - core/sek00_summary (1x | ref:1)
  - core/sek02_pole (1x | ref:1)
  - core/sek03_rezimy (1x | ref:1)
  - core/sek07_predykcje (3x | ref:3)
  - core/sek08_formalizm (5x | ref:5)
  - partial_proofs/bh_ringdown (1x | ref:1)
  - partial_proofs/wielki_wybuch (3x | ref:3)
  - research/nbody/paper (2x | ref:2)

### core/sek07_predykcje
- depended on by:
  - <root> (1x | input:1)
  - core/_meta_latex (2x | ref:2)
  - core/formalizm (2x | ref:2)
  - core/sek00_summary (2x | ref:2)
  - core/sek05_ciemna_energia (1x | ref:1)
  - core/sek06_czarne_dziury (1x | ref:1)
  - core/sek07a_wymiar_wzmocniony (2x | ref:2)
  - core/sek08_formalizm (7x | ref:7)
  - core/sek09_cechowanie (4x | ref:4)
  - partial_proofs/defect_hierarchy (1x | ref:1)
  - partial_proofs/koide_fp (1x | ref:1)
  - partial_proofs/wielki_wybuch (2x | ref:2)
  - partial_proofs/zero_mode (1x | ref:1)

### core/sek07a_wymiar_wzmocniony
- depended on by:
  - <root> (1x | input:1)
  - core/formalizm (1x | ref:1)

### core/sek08_formalizm
- depended on by:
  - <root> (4x | input:1, ref:3)
  - audyt/L04_ODE_dualism_alpha (2x | wikilink:2)
  - audyt/L07_zero_sum_axiom (2x | wikilink:2)
  - audyt/L08_kink_fermion_closure (2x | wikilink:2)
  - axioms/notacja (35x | ref:35)
  - axioms/roznica_N0 (1x | ref:1)
  - axioms/substrat (15x | ref:15)
  - core/_meta_latex (40x | ref:40)
  - core/formalizm (22x | ref:22)
  - core/sek00_summary (13x | ref:13)
  - core/sek01_ontologia (21x | ref:21)
  - core/sek02_pole (9x | ref:9)
  - core/sek03_rezimy (23x | ref:23)
  - core/sek04_stale (11x | ref:11)
  - core/sek05_ciemna_energia (11x | ref:11)
  - core/sek06_czarne_dziury (7x | ref:7)
  - core/sek07_predykcje (53x | ref:53)
  - core/sek08a_akcja_zunifikowana (20x | ref:20)
  - core/sek08b_ghost_resolution (11x | ref:11)
  - core/sek08c_metryka_z_substratu (7x | ref:7)
  - core/sek09_cechowanie (10x | ref:10)
  - core/sek10_N0_wyprowadzenie (1x | ref:1)
  - papers_external/tgp_core_paper (2x | ref:2)
  - papers_external/tgp_english_summary (18x | ref:18)
  - partial_proofs/alphaK_status (3x | ref:3)
  - partial_proofs/bh_ringdown (10x | ref:10)
  - partial_proofs/defect_hierarchy (1x | ref:1)
  - partial_proofs/hierarchia_mas (3x | ref:3)
  - partial_proofs/most_gamma_phi (1x | ref:1)
  - partial_proofs/quark_sector (3x | ref:3)
  - partial_proofs/trojcialowe_nbody (3x | ref:3)
  - partial_proofs/wielki_wybuch (9x | ref:9)
  - research/nbody/paper (12x | ref:12)
  - research/op-L01-rho-stress-energy-bridge-2026-05-04 (8x | wikilink:8)
  - research/op-L03-spectral-stability-2026-05-06 (7x | wikilink:7)
  - research/op-L04-ODE-canonicalization-2026-05-04 (9x | wikilink:9)
  - research/op-chi1-newton-constant-derivation (1x | wikilink:1)
  - research/op-uv3-phi0-renormalization (1x | wikilink:1)

### core/sek08a_akcja_zunifikowana
- depended on by:
  - <root> (1x | input:1)
  - audyt/L01_rho_operational (4x | wikilink:4)
  - audyt/S07_M911_derivation (3x | wikilink:3)
  - audyt/T01_LIGO3G_falsifier (1x | wikilink:1)
  - axioms/notacja (3x | ref:3)
  - core/_meta_latex (6x | ref:6)
  - core/formalizm (4x | ref:4)
  - core/sek00_summary (3x | ref:3)
  - core/sek06_czarne_dziury (1x | ref:1)
  - core/sek08_formalizm (11x | ref:11)
  - core/sek08b_ghost_resolution (2x | ref:2)
  - core/sek08c_metryka_z_substratu (4x | ref:4)
  - core/sek09_cechowanie (2x | ref:2)
  - partial_proofs/most_gamma_phi (1x | ref:1)
  - research/nbody/paper (1x | ref:1)
  - research/op-CORE-CLEANUP-B-2026-05-04 (3x | wikilink:3)
  - research/op-FRW-radiation-era-varying-c-2026-05-06 (3x | wikilink:3)
  - research/op-L01-rho-stress-energy-bridge-2026-05-04 (4x | wikilink:4)
  - research/op-L03-spectral-stability-2026-05-06 (7x | wikilink:7)
  - research/op-L04-ODE-canonicalization-2026-05-04 (2x | wikilink:2)
  - research/op-MAG-Phase5-V-reference-clarification-2026-05-09 (2x | wikilink:2)
  - research/op-MAG-resonance-formalization-2026-05-09 (1x | wikilink:1)
  - research/op-Phi-decomposition-photon-2026-05-07 (4x | wikilink:4)
  - research/op-Phi-vacuum-scale-2026-05-09 (4x | wikilink:4)
  - research/op-S07-alternative-f-psi-derivation-2026-05-09 (1x | wikilink:1)
  - research/op-SPIN-SU2-substrate-derivation-2026-05-08 (3x | wikilink:3)
  - research/op-V-canonical-consistency-audit-2026-05-09 (3x | wikilink:3)
  - research/op-cosmology-closure (4x | wikilink:4)
  - research/op-dual-V-structure-clarification-2026-05-09 (3x | wikilink:3)
  - research/op-g0-r3-from-canonical-projection (1x | wikilink:1)
  - research/op-omicron2-phi-mean-shift-cosmo (1x | wikilink:1)
  - research/op-quantum-closure (2x | wikilink:2)
  - research/op-void-flat-modes-h0-2026-05-06 (1x | wikilink:1)

### core/sek08b_ghost_resolution
- depended on by:
  - <root> (1x | input:1)
  - audyt (1x | wikilink:1)
  - audyt/L03_K_phi_stability (2x | wikilink:2)
  - core/_meta_latex (3x | ref:3)
  - core/formalizm (3x | ref:3)
  - core/sek00_summary (1x | ref:1)
  - core/sek07_predykcje (1x | ref:1)
  - core/sek09_cechowanie (1x | ref:1)
  - core/sek10_N0_wyprowadzenie (3x | ref:3)
  - partial_proofs/chiralnosc (1x | ref:1)
  - partial_proofs/hierarchia_mas (3x | ref:3)
  - research/op-L03-spectral-stability-2026-05-06 (7x | wikilink:7)
  - research/op-L04-ODE-canonicalization-2026-05-04 (1x | wikilink:1)

### core/sek08c_metryka_z_substratu
- depended on by:
  - <root> (1x | input:1)
  - audyt (2x | wikilink:2)
  - audyt/S01_metric_four_forms (2x | wikilink:2)
  - audyt/S07_M911_derivation (3x | wikilink:3)
  - audyt/T01_LIGO3G_falsifier (4x | wikilink:4)
  - core/_meta_latex (1x | ref:1)
  - core/formalizm (3x | ref:3)
  - core/sek00_summary (2x | ref:2)
  - core/sek02_pole (1x | ref:1)
  - core/sek04_stale (1x | ref:1)
  - core/sek08_formalizm (1x | ref:1)
  - research/op-CORE-CLEANUP-B-2026-05-04 (4x | wikilink:4)
  - research/op-Phi-vacuum-scale-2026-05-09 (1x | wikilink:1)
  - research/op-S07-alternative-f-psi-derivation-2026-05-09 (1x | wikilink:1)
  - research/op-g0-r3-from-canonical-projection (1x | wikilink:1)
  - research/why_n3 (1x | wikilink:1)

### core/sek09_cechowanie
- depended on by:
  - <root> (1x | input:1)
  - core/_meta_latex (12x | ref:12)
  - core/formalizm (51x | ref:51)
  - core/sek00_summary (1x | ref:1)
  - core/sek01_ontologia (1x | ref:1)
  - core/sek07_predykcje (11x | ref:11)
  - core/sek08_formalizm (2x | ref:2)
  - papers_external/paper_lepton_masses (1x | ref:1)
  - partial_proofs/defect_hierarchy (16x | ref:16)
  - partial_proofs/fermion_from_soliton (2x | wikilink:2)
  - partial_proofs/hierarchia_mas (1x | ref:1)
  - research/atom_from_soliton (1x | wikilink:1)

### core/sek10_N0_wyprowadzenie
- depended on by:
  - <root> (1x | input:1)
  - axioms/notacja (1x | ref:1)
  - core/_meta_latex (1x | ref:1)
  - core/formalizm (2x | ref:2)
  - core/sek00_summary (4x | ref:4)
  - core/sek08b_ghost_resolution (2x | ref:2)
  - partial_proofs/most_gamma_phi (1x | ref:1)
  - partial_proofs/nuclear_from_soliton (2x | wikilink:2)

### meta
- depended on by:
  - <root> (8x | wikilink:8)
  - audyt (37x | wikilink:37)
  - audyt/D01_drifting_numbers (5x | wikilink:5)
  - audyt/L01_rho_operational (2x | wikilink:2)
  - audyt/L02_beta_gamma_semantics (5x | wikilink:5)
  - audyt/L03_K_phi_stability (3x | wikilink:3)
  - audyt/L04_ODE_dualism_alpha (6x | wikilink:6)
  - audyt/L05_mass_exponent_drift (4x | wikilink:4)
  - audyt/L06_axion_mass_locked (3x | wikilink:3)
  - audyt/M01_status_creep (7x | wikilink:7)
  - audyt/M02_ledger_pollution (9x | wikilink:9)
  - audyt/M03_balance_sheet_missing (24x | wikilink:24)
  - audyt/S01_metric_four_forms (4x | wikilink:4)
  - audyt/S02_volume_element_M9 (3x | wikilink:3)
  - audyt/S03_beta_PPN_convention (3x | wikilink:3)
  - audyt/S04_metric_coupling_axiom (3x | wikilink:3)
  - audyt/S05_tensor_sector_singleField (3x | wikilink:3)
  - audyt/S06_circular_anchors (12x | wikilink:12)
  - audyt/T01_LIGO3G_falsifier (3x | wikilink:3)
  - meta/core (4x | wikilink:4)
  - meta/core/intake (1x | wikilink:1)
  - meta/research (30x | wikilink:30)
  - meta/research/_examples (10x | wikilink:10)
  - meta/research/templates (3x | wikilink:3)
  - research/_archive (3x | wikilink:3)
  - research/_sandbox (2x | wikilink:2)
  - research/closure_2026-04-26 (5x | wikilink:5)
  - research/op-D01-anchor-lock-2026-05-06 (7x | wikilink:7)
  - research/op-FRW-radiation-era-varying-c-2026-05-06 (7x | wikilink:7)
  - research/op-GWTC3-reanalysis (1x | wikilink:1)
  - research/op-L01-rho-stress-energy-bridge-2026-05-04 (4x | wikilink:4)
  - research/op-L03-spectral-stability-2026-05-06 (1x | wikilink:1)
  - research/op-L04-ODE-canonicalization-2026-05-04 (8x | wikilink:8)
  - research/op-LIGO-3G-deviation (1x | wikilink:1)
  - research/op-M03-balance-sheet-retrofit-2026-05-06 (100x | wikilink:100)
  - research/op-MAG-Lorentz-A-mu-coupling-2026-05-09 (1x | wikilink:1)
  - research/op-MAG-resonance-formalization-2026-05-09 (2x | wikilink:2)
  - research/op-Phi-decomposition-photon-2026-05-07 (8x | wikilink:8)
  - research/op-Phi-vacuum-scale-2026-05-09 (4x | wikilink:4)
  - research/op-Phi0-spatial-variation-predictions-2026-05-09 (1x | wikilink:1)
  - research/op-S07-alternative-f-psi-derivation-2026-05-09 (1x | wikilink:1)
  - research/op-SPIN-MAG-leakage-lean-hypothesis-2026-05-07 (1x | wikilink:1)
  - research/op-SPIN-SU2-substrate-derivation-2026-05-08 (2x | wikilink:2)
  - research/op-V-canonical-consistency-audit-2026-05-09 (2x | wikilink:2)
  - research/op-chi1-newton-constant-derivation (4x | wikilink:4)
  - research/op-cosmology-closure (2x | wikilink:2)
  - research/op-dual-V-structure-clarification-2026-05-09 (2x | wikilink:2)
  - research/op-g0-r3-from-canonical-projection (1x | wikilink:1)
  - research/op-newton-momentum (13x | wikilink:13)
  - research/op-omega2-axion-coupling-lock (7x | wikilink:7)
  - research/op-omega3-axion-decay-constant (8x | wikilink:8)
  - research/op-ppE-mapping (1x | wikilink:1)
  - research/op-ppE-mapping/TGP/TGP_v1/research/op-LIGO-3G-deviation (1x | wikilink:1)
  - research/op-psi1-substrate-light-acceleration (3x | wikilink:3)
  - research/op-tau3-substrate-clock-acceleration (6x | wikilink:6)
  - research/op-tensor-modes-Phi-FUTURE (1x | wikilink:1)
  - research/op-uv2-mtgp-absolute-scale (4x | wikilink:4)

### meta/core
- depended on by:
  - meta/core/intake (23x | wikilink:23)
  - meta/research (7x | wikilink:7)
  - meta/research/_examples (1x | wikilink:1)

### meta/core/intake
- depended on by:
  - meta/core (4x | wikilink:4)
  - meta/research (9x | wikilink:9)

### meta/research
- depended on by:
  - <root> (1x | wikilink:1)
  - audyt (2x | wikilink:2)
  - audyt/M03_balance_sheet_missing (3x | wikilink:3)
  - meta (4x | wikilink:4)
  - meta/core (8x | wikilink:8)
  - meta/research/_examples (7x | wikilink:7)
  - meta/research/templates (6x | wikilink:6)
  - research/_archive (4x | wikilink:4)
  - research/_sandbox (2x | wikilink:2)
  - research/atom_from_soliton (5x | wikilink:5)
  - research/atomic_shells_closure (5x | wikilink:5)
  - research/brannen_sqrt2 (4x | wikilink:4)
  - research/cabibbo_correction (4x | wikilink:4)
  - research/casimir_mof (5x | wikilink:5)
  - research/closure_2026-04-26 (5x | wikilink:5)
  - research/cohesion_closure (5x | wikilink:5)
  - research/continuum_limit (4x | wikilink:4)
  - research/cosmo_tensions (4x | wikilink:4)
  - research/desi_dark_energy (4x | wikilink:4)
  - research/em_from_substrate (5x | wikilink:5)
  - research/external_review_2026-04-25 (5x | wikilink:5)
  - research/galaxy_scaling (4x | wikilink:4)
  - research/hubble_tension (4x | wikilink:4)
  - research/liquid_viscosity (5x | wikilink:5)
  - research/mass_scaling_k4 (4x | wikilink:4)
  - research/metric_ansatz (4x | wikilink:4)
  - research/muon_g_minus_2 (5x | wikilink:5)
  - research/nbody (4x | wikilink:4)
  - research/neutrino_msw (5x | wikilink:5)
  - research/op-alpha-fine-structure (5x | wikilink:5)
  - research/op-bh-alpha-threshold (5x | wikilink:5)
  - research/op-chi1-newton-constant-derivation (1x | wikilink:1)
  - research/op-cosmology-closure (5x | wikilink:5)
  - research/op-cross-sector-charge (5x | wikilink:5)
  - research/op-delta1-g-tilde-derivation (4x | wikilink:4)
  - research/op-delta2-Nf-derivation (4x | wikilink:4)
  - research/op-eht (5x | wikilink:5)
  - research/op-eht-A (5x | wikilink:5)
  - research/op-eps-photon-ring (5x | wikilink:5)
  - research/op-eta-wolfenstein (5x | wikilink:5)
  - research/op-eta2-denom-derivation (5x | wikilink:5)
  - research/op-g0-r3-from-canonical-projection (4x | wikilink:4)
  - research/op-gamma1-phi-eff-anchor-resolution (4x | wikilink:4)
  - research/op-iota-charge-pmns-unification (5x | wikilink:5)
  - research/op-kappa-mixing-numerator (5x | wikilink:5)
  - research/op-lambda1-e2-amplitude-emergence (4x | wikilink:4)
  - research/op-m92 (5x | wikilink:5)
  - research/op-mu-pmns-phase-hardening (5x | wikilink:5)
  - research/op-mu1-minimal-substrate-log-redefinition (4x | wikilink:4)
  - research/op-newton-momentum (5x | wikilink:5)
  - research/op-nu-majorana-phase-mbb (5x | wikilink:5)
  - research/op-omega1-substrate-em-coupling (5x | wikilink:5)
  - research/op-omega2-axion-coupling-lock (1x | wikilink:1)
  - research/op-omega3-axion-decay-constant (1x | wikilink:1)
  - research/op-omicron1-sigmamnu-cosmo (5x | wikilink:5)
  - research/op-phase1-covariant (5x | wikilink:5)
  - research/op-phase2-quantum-gravity (5x | wikilink:5)
  - research/op-phase3-uv-completion (5x | wikilink:5)
  - research/op-phi1-substrate-action-variational (5x | wikilink:5)
  - research/op-pi1-bb0nu-nme-isotope (5x | wikilink:5)
  - research/op-psi1-substrate-light-acceleration (5x | wikilink:5)
  - research/op-quantum-closure (5x | wikilink:5)
  - research/op-rho1-71Ge-cross-section (5x | wikilink:5)
  - research/op-sc-alpha-origin (5x | wikilink:5)
  - research/op-sigma1-substrate-light-dispersion (5x | wikilink:5)
  - research/op-tau1-closure-overlap-coulomb (5x | wikilink:5)
  - research/op-tau2-substrate-time-coupling (5x | wikilink:5)
  - research/op-tau3-substrate-clock-acceleration (5x | wikilink:5)
  - research/op-theta-quark-koide (5x | wikilink:5)
  - research/op-upsilon1-closure-cross-family (5x | wikilink:5)
  - research/op-uv-as-ngfp (5x | wikilink:5)
  - research/op-uv-renormalizability-research (4x | wikilink:4)
  - research/op-uv2-mtgp-absolute-scale (1x | wikilink:1)
  - research/op-uv3-phi0-renormalization (4x | wikilink:4)
  - research/op-void-flat-modes-h0-2026-05-06 (2x | wikilink:2)
  - research/op-xi-photon-ring (5x | wikilink:5)
  - research/op-xi2-sterile-nu-5sector (5x | wikilink:5)
  - research/op-zeta-mass-spectrum (5x | wikilink:5)
  - research/op1-op2-op4 (4x | wikilink:4)
  - research/op6 (4x | wikilink:4)
  - research/op7 (5x | wikilink:5)
  - research/particle_sector_closure (4x | wikilink:4)
  - research/qm_born_rule (4x | wikilink:4)
  - research/qm_decoherence (4x | wikilink:4)
  - research/qm_entanglement (4x | wikilink:4)
  - research/qm_foundations (4x | wikilink:4)
  - research/qm_measurement (4x | wikilink:4)
  - research/qm_spin (4x | wikilink:4)
  - research/qm_statistics (4x | wikilink:4)
  - research/qm_superposition (4x | wikilink:4)
  - research/rho_normal_state_closure (5x | wikilink:5)
  - research/s8_tension (4x | wikilink:4)
  - research/superconductivity_closure (4x | wikilink:4)
  - research/thermal_transport_molecular (5x | wikilink:5)
  - research/uv_completion (4x | wikilink:4)
  - research/why_n3 (4x | wikilink:4)

### meta/research/_examples
- depended on by: --

### meta/research/templates
- depended on by:
  - audyt (1x | wikilink:1)
  - meta/research (6x | wikilink:6)
  - research/_sandbox (1x | wikilink:1)
  - research/atom_from_soliton (1x | wikilink:1)
  - research/atomic_shells_closure (1x | wikilink:1)
  - research/casimir_mof (1x | wikilink:1)
  - research/closure_2026-04-26 (1x | wikilink:1)
  - research/cohesion_closure (1x | wikilink:1)
  - research/em_from_substrate (1x | wikilink:1)
  - research/external_review_2026-04-25 (1x | wikilink:1)
  - research/liquid_viscosity (1x | wikilink:1)
  - research/muon_g_minus_2 (1x | wikilink:1)
  - research/neutrino_msw (1x | wikilink:1)
  - research/op-alpha-fine-structure (1x | wikilink:1)
  - research/op-bh-alpha-threshold (1x | wikilink:1)
  - research/op-chi1-newton-constant-derivation (1x | wikilink:1)
  - research/op-cosmology-closure (1x | wikilink:1)
  - research/op-cross-sector-charge (1x | wikilink:1)
  - research/op-eht (1x | wikilink:1)
  - research/op-eht-A (1x | wikilink:1)
  - research/op-eps-photon-ring (1x | wikilink:1)
  - research/op-eta-wolfenstein (1x | wikilink:1)
  - research/op-eta2-denom-derivation (1x | wikilink:1)
  - research/op-iota-charge-pmns-unification (1x | wikilink:1)
  - research/op-kappa-mixing-numerator (1x | wikilink:1)
  - research/op-m92 (1x | wikilink:1)
  - research/op-mu-pmns-phase-hardening (1x | wikilink:1)
  - research/op-newton-momentum (1x | wikilink:1)
  - research/op-nu-majorana-phase-mbb (1x | wikilink:1)
  - research/op-omega1-substrate-em-coupling (1x | wikilink:1)
  - research/op-omega2-axion-coupling-lock (1x | wikilink:1)
  - research/op-omega3-axion-decay-constant (1x | wikilink:1)
  - research/op-omicron1-sigmamnu-cosmo (1x | wikilink:1)
  - research/op-phase1-covariant (1x | wikilink:1)
  - research/op-phase2-quantum-gravity (1x | wikilink:1)
  - research/op-phase3-uv-completion (1x | wikilink:1)
  - research/op-phi1-substrate-action-variational (1x | wikilink:1)
  - research/op-pi1-bb0nu-nme-isotope (1x | wikilink:1)
  - research/op-psi1-substrate-light-acceleration (1x | wikilink:1)
  - research/op-quantum-closure (1x | wikilink:1)
  - research/op-rho1-71Ge-cross-section (1x | wikilink:1)
  - research/op-sc-alpha-origin (1x | wikilink:1)
  - research/op-sigma1-substrate-light-dispersion (1x | wikilink:1)
  - research/op-tau1-closure-overlap-coulomb (1x | wikilink:1)
  - research/op-tau2-substrate-time-coupling (1x | wikilink:1)
  - research/op-tau3-substrate-clock-acceleration (1x | wikilink:1)
  - research/op-theta-quark-koide (1x | wikilink:1)
  - research/op-upsilon1-closure-cross-family (1x | wikilink:1)
  - research/op-uv-as-ngfp (1x | wikilink:1)
  - research/op-uv2-mtgp-absolute-scale (1x | wikilink:1)
  - research/op-xi-photon-ring (1x | wikilink:1)
  - research/op-xi2-sterile-nu-5sector (1x | wikilink:1)
  - research/op-zeta-mass-spectrum (1x | wikilink:1)
  - research/op7 (1x | wikilink:1)
  - research/rho_normal_state_closure (1x | wikilink:1)
  - research/thermal_transport_molecular (1x | wikilink:1)

### papers/M911_LIGO3G_paper
- depended on by:
  - audyt/T01_LIGO3G_falsifier (5x | wikilink:5)
  - research/op-GWTC3-reanalysis (3x | wikilink:3)
  - research/op-ppE-mapping (3x | wikilink:3)
  - research/op-ppE-mapping/TGP/TGP_v1/audyt/T01_LIGO3G_falsifier (2x | wikilink:2)

### papers_external
- depended on by: --

### papers_external/arxiv_submission
- depended on by:
  - research/nbody/paper (24x | ref:24)

### papers_external/paper_bh_shadow
- depended on by:
  - <root> (2x | ref:2)
  - papers_external/paper_lepton_masses (1x | ref:1)
  - papers_external/tgp_core_paper (6x | ref:6)
  - papers_external/tgp_english_summary (3x | ref:3)
  - papers_external/tgp_sc_paper (3x | ref:3)

### papers_external/paper_lepton_masses
- depended on by:
  - papers_external/tgp_sc_paper (3x | ref:3)

### papers_external/tgp_core_paper
- depended on by:
  - papers_external/tgp_english_summary (2x | ref:2)
  - research/atomic_shells_closure (1x | wikilink:1)
  - research/casimir_mof (1x | wikilink:1)
  - research/closure_2026-04-26/f_psi_principle (1x | wikilink:1)
  - research/closure_2026-04-26/sigma_ab_pathB (1x | wikilink:1)
  - research/external_review_2026-04-25 (2x | wikilink:2)
  - research/liquid_viscosity (1x | wikilink:1)
  - research/muon_g_minus_2 (1x | wikilink:1)
  - research/neutrino_msw (1x | wikilink:1)
  - research/op-eht (5x | wikilink:5)
  - research/op-eht-A (2x | wikilink:2)
  - research/op-m92 (3x | wikilink:3)
  - research/op7 (5x | wikilink:5)
  - research/thermal_transport_molecular (1x | wikilink:1)

### papers_external/tgp_english_summary
- depended on by:
  - papers_external/paper_bh_shadow (1x | ref:1)
  - partial_proofs/hierarchia_mas (1x | ref:1)
  - partial_proofs/zero_mode (1x | ref:1)

### papers_external/tgp_sc_paper
- depended on by:
  - <root> (2x | ref:2)
  - research/superconductivity_closure (1x | wikilink:1)
  - research/thermal_transport_molecular (1x | wikilink:1)

### partial_proofs
- depended on by: --

### partial_proofs/alphaK_status
- depended on by:
  - <root> (1x | input:1)
  - axioms/notacja (1x | ref:1)
  - core/_meta_latex (12x | ref:12)
  - core/sek08_formalizm (5x | ref:5)

### partial_proofs/bh_ringdown
- depended on by:
  - <root> (2x | input:2)
  - axioms/notacja (1x | ref:1)
  - core/_meta_latex (12x | ref:12)
  - core/formalizm (1x | ref:1)
  - core/sek00_summary (4x | ref:4)
  - core/sek06_czarne_dziury (1x | ref:1)
  - core/sek07_predykcje (6x | ref:6)
  - core/sek08_formalizm (2x | ref:2)
  - partial_proofs/hierarchia_mas (22x | ref:22)
  - partial_proofs/koide_fp (1x | ref:1)
  - partial_proofs/zero_mode (1x | ref:1)

### partial_proofs/bounce_topology
- depended on by:
  - <root> (1x | input:1)
  - partial_proofs/superconductivity (2x | ref:2)

### partial_proofs/chiralnosc
- depended on by:
  - <root> (1x | input:1)
  - core/formalizm (3x | ref:3)
  - core/sek09_cechowanie (4x | ref:4)
  - partial_proofs/defect_hierarchy (1x | ref:1)

### partial_proofs/defect_hierarchy
- depended on by:
  - <root> (1x | input:1)
  - core/_meta_latex (4x | ref:4)
  - core/sek09_cechowanie (10x | ref:10)
  - partial_proofs/fermion_from_soliton (1x | wikilink:1)

### partial_proofs/fermion_from_soliton
- depended on by:
  - partial_proofs/nuclear_from_soliton (2x | wikilink:2)
  - research/atom_from_soliton (1x | wikilink:1)
  - research/atomic_shells_closure (1x | wikilink:1)
  - research/casimir_mof (1x | wikilink:1)
  - research/cohesion_closure (1x | wikilink:1)
  - research/em_from_substrate (1x | wikilink:1)
  - research/liquid_viscosity (1x | wikilink:1)
  - research/muon_g_minus_2 (1x | wikilink:1)
  - research/op-CORE-CLEANUP-B-2026-05-04 (1x | wikilink:1)
  - research/op-M03-balance-sheet-retrofit-2026-05-06 (3x | wikilink:3)
  - research/op-gamma1-phi-eff-anchor-resolution (1x | wikilink:1)
  - research/op-mu1-minimal-substrate-log-redefinition (3x | wikilink:3)
  - research/op-omega1-substrate-em-coupling (2x | wikilink:2)
  - research/thermal_transport_molecular (1x | wikilink:1)

### partial_proofs/hierarchia_mas
- depended on by:
  - <root> (6x | input:4, ref:2)
  - axioms/notacja (1x | ref:1)
  - core/_meta_latex (15x | ref:15)
  - core/formalizm (18x | ref:18)
  - core/sek00_summary (1x | ref:1)
  - core/sek04_stale (1x | ref:1)
  - core/sek07_predykcje (3x | ref:3)
  - core/sek08_formalizm (2x | ref:2)
  - core/sek08b_ghost_resolution (11x | ref:11)
  - core/sek09_cechowanie (6x | ref:6)
  - partial_proofs/alphaK_status (4x | ref:4)
  - partial_proofs/bh_ringdown (14x | ref:14)
  - partial_proofs/chiralnosc (1x | ref:1)
  - partial_proofs/koide_fp (8x | ref:8)
  - partial_proofs/most_gamma_phi (5x | ref:5)
  - partial_proofs/nuclear_from_soliton (2x | wikilink:2)
  - partial_proofs/particle_sector (2x | ref:2)
  - partial_proofs/potencjal (1x | ref:1)
  - partial_proofs/quark_sector (5x | ref:5)
  - partial_proofs/superconductivity (2x | ref:2)
  - partial_proofs/zero_mode (8x | ref:8)
  - research/op-delta2-Nf-derivation (1x | wikilink:1)

### partial_proofs/koide_fp
- depended on by:
  - <root> (5x | input:5)
  - core/_meta_latex (4x | ref:4)
  - core/formalizm (4x | ref:4)
  - core/sek09_cechowanie (1x | ref:1)
  - partial_proofs/most_gamma_phi (4x | ref:4)
  - partial_proofs/particle_sector (3x | ref:3)
  - partial_proofs/zero_mode (8x | ref:8)

### partial_proofs/most_gamma_phi
- depended on by:
  - <root> (2x | input:2)
  - core/formalizm (2x | ref:2)
  - partial_proofs/alphaK_status (1x | ref:1)

### partial_proofs/nuclear_from_soliton
- depended on by:
  - partial_proofs/fermion_from_soliton (5x | wikilink:5)

### partial_proofs/particle_sector
- depended on by:
  - <root> (2x | input:2)
  - partial_proofs/hierarchia_mas (1x | ref:1)
  - partial_proofs/koide_fp (2x | ref:2)
  - partial_proofs/superconductivity (1x | ref:1)

### partial_proofs/potencjal
- depended on by:
  - <root> (1x | input:1)
  - axioms/notacja (1x | ref:1)
  - core/_meta_latex (1x | ref:1)
  - core/formalizm (3x | ref:3)
  - core/sek08_formalizm (1x | ref:1)
  - partial_proofs/koide_fp (1x | ref:1)

### partial_proofs/quark_sector
- depended on by:
  - <root> (1x | input:1)
  - core/_meta_latex (12x | ref:12)
  - core/formalizm (3x | ref:3)
  - core/sek09_cechowanie (1x | ref:1)
  - partial_proofs/particle_sector (6x | ref:6)

### partial_proofs/superconductivity
- depended on by:
  - <root> (2x | input:2)

### partial_proofs/trojcialowe_nbody
- depended on by:
  - <root> (2x | input:2)
  - core/sek07_predykcje (1x | ref:1)
  - core/sek08_formalizm (2x | ref:2)
  - research/nbody/paper (1x | ref:1)

### partial_proofs/wielki_wybuch
- depended on by:
  - <root> (1x | input:1)
  - core/_meta_latex (1x | ref:1)
  - core/formalizm (3x | ref:3)
  - core/sek00_summary (1x | ref:1)
  - core/sek08_formalizm (6x | ref:6)

### partial_proofs/zero_mode
- depended on by:
  - <root> (2x | input:2)
  - axioms/roznica_N0 (5x | ref:5)
  - core/_meta_latex (4x | ref:4)
  - core/formalizm (14x | ref:14)
  - partial_proofs/hierarchia_mas (4x | ref:4)
  - partial_proofs/koide_fp (2x | ref:2)
  - partial_proofs/quark_sector (1x | ref:1)
  - partial_proofs/superconductivity (1x | ref:1)

### research
- depended on by:
  - <root> (21x | wikilink:21)
  - research/desi_dark_energy (1x | wikilink:1)
  - research/op-eht (1x | wikilink:1)

### research/_archive
- depended on by: --

### research/_sandbox
- depended on by: --

### research/atom_from_soliton
- depended on by:
  - partial_proofs/fermion_from_soliton (2x | wikilink:2)
  - partial_proofs/nuclear_from_soliton (2x | wikilink:2)

### research/atomic_shells_closure
- depended on by:
  - research/atom_from_soliton (2x | wikilink:2)
  - research/cohesion_closure (2x | wikilink:2)

### research/audyt_cosmology_drift_2026-05-03
- depended on by:
  - research/op-CORE-CLEANUP-B-2026-05-04 (5x | wikilink:5)

### research/brannen_sqrt2
- depended on by:
  - research/particle_sector_closure (1x | wikilink:1)

### research/cabibbo_correction
- depended on by: --

### research/casimir_mof
- depended on by:
  - research/muon_g_minus_2 (3x | wikilink:3)

### research/closure_2026-04-26
- depended on by:
  - <root> (1x | wikilink:1)
  - audyt/S05_tensor_sector_singleField (6x | wikilink:6)
  - meta (5x | wikilink:5)
  - meta/core/intake (2x | wikilink:2)
  - meta/research (1x | wikilink:1)
  - meta/research/_examples (3x | wikilink:3)
  - research/audyt_cosmology_drift_2026-05-03 (1x | wikilink:1)
  - research/closure_2026-04-26/sigma_ab_pathB (1x | wikilink:1)
  - research/desi_dark_energy (1x | wikilink:1)
  - research/external_review_2026-04-25 (1x | wikilink:1)
  - research/op-M03-balance-sheet-retrofit-2026-05-06 (1x | wikilink:1)
  - research/op-cosmology-closure (10x | wikilink:10)
  - research/op-eht (3x | wikilink:3)
  - research/op-eht-A (1x | wikilink:1)
  - research/op-m92 (1x | wikilink:1)
  - research/op-newton-momentum (6x | wikilink:6)
  - research/op-phase1-covariant (5x | wikilink:5)
  - research/op-phase2-quantum-gravity (10x | wikilink:10)
  - research/op-phase3-uv-completion (3x | wikilink:3)
  - research/op-quantum-closure (7x | wikilink:7)
  - research/op-uv-renormalizability-research (2x | wikilink:2)
  - research/op-void-flat-modes-h0-2026-05-06 (1x | wikilink:1)
  - research/op7 (1x | wikilink:1)

### research/closure_2026-04-26/Lambda_from_Phi0
- depended on by:
  - research/op-L01-rho-stress-energy-bridge-2026-05-04 (4x | wikilink:4)
  - research/op-quantum-closure (2x | wikilink:2)

### research/closure_2026-04-26/alpha_psi_threshold
- depended on by:
  - audyt (1x | wikilink:1)
  - audyt/S05_tensor_sector_singleField (3x | wikilink:3)
  - meta/core/intake (8x | wikilink:8)
  - meta/research/_examples (4x | wikilink:4)
  - research/audyt_cosmology_drift_2026-05-03 (2x | wikilink:2)
  - research/closure_2026-04-26 (10x | wikilink:10)
  - research/closure_2026-04-26/Lambda_from_Phi0 (1x | wikilink:1)
  - research/closure_2026-04-26/sigma_ab_pathB (1x | wikilink:1)
  - research/op-M03-balance-sheet-retrofit-2026-05-06 (10x | wikilink:10)
  - research/op-Phi-vacuum-scale-2026-05-09 (2x | wikilink:2)
  - research/op-V-canonical-consistency-audit-2026-05-09 (5x | wikilink:5)
  - research/op-bh-alpha-threshold (1x | wikilink:1)
  - research/op-cosmology-closure (14x | wikilink:14)
  - research/op-delta1-g-tilde-derivation (2x | wikilink:2)
  - research/op-delta2-Nf-derivation (1x | wikilink:1)
  - research/op-dual-V-structure-clarification-2026-05-09 (2x | wikilink:2)
  - research/op-gamma1-phi-eff-anchor-resolution (2x | wikilink:2)
  - research/op-newton-momentum (11x | wikilink:11)
  - research/op-omicron2-phi-mean-shift-cosmo (4x | wikilink:4)
  - research/op-phase1-covariant (9x | wikilink:9)
  - research/op-phase2-quantum-gravity (6x | wikilink:6)
  - research/op-phase3-uv-completion (6x | wikilink:6)
  - research/op-quantum-closure (2x | wikilink:2)
  - research/op-sc-alpha-origin (2x | wikilink:2)
  - research/op-void-flat-modes-h0-2026-05-06 (10x | wikilink:10)

### research/closure_2026-04-26/f_psi_principle
- depended on by: --

### research/closure_2026-04-26/sigma_ab_pathB
- depended on by:
  - research/closure_2026-04-26 (1x | wikilink:1)

### research/cohesion_closure
- depended on by:
  - partial_proofs/nuclear_from_soliton (1x | wikilink:1)
  - research/atom_from_soliton (1x | wikilink:1)

### research/continuum_limit
- depended on by:
  - research/op-quantum-closure (10x | wikilink:10)

### research/cosmo_tensions
- depended on by:
  - research/galaxy_scaling (2x | wikilink:2)
  - research/op-cosmology-closure (16x | wikilink:16)

### research/desi_dark_energy
- depended on by:
  - research/op-cosmology-closure (7x | wikilink:7)

### research/em_from_substrate
- depended on by:
  - research/atom_from_soliton (5x | wikilink:5)
  - research/cohesion_closure (2x | wikilink:2)

### research/external_review_2026-04-25
- depended on by: --

### research/galaxy_scaling
- depended on by:
  - research/desi_dark_energy (3x | wikilink:3)
  - research/op-cosmology-closure (17x | wikilink:17)

### research/hubble_tension
- depended on by: --

### research/liquid_viscosity
- depended on by: --

### research/mass_scaling_k4
- depended on by:
  - audyt/L04_ODE_dualism_alpha (2x | wikilink:2)
  - research/op-L04-ODE-canonicalization-2026-05-04 (12x | wikilink:12)
  - research/op-g0-r3-from-canonical-projection (1x | wikilink:1)
  - research/op-lambda1-e2-amplitude-emergence (8x | wikilink:8)
  - research/why_n3 (1x | wikilink:1)

### research/metric_ansatz
- depended on by: --

### research/muon_g_minus_2
- depended on by:
  - research/op-quantum-closure (1x | wikilink:1)

### research/nbody
- depended on by:
  - partial_proofs/fermion_from_soliton (1x | wikilink:1)
  - partial_proofs/nuclear_from_soliton (8x | wikilink:8)

### research/nbody/docs
- depended on by: --

### research/nbody/examples
- depended on by:
  - research/em_from_substrate (1x | wikilink:1)
  - research/op-cosmology-closure (8x | wikilink:8)

### research/nbody/examples/_outputs
- depended on by: --

### research/nbody/paper
- depended on by:
  - core/sek02_pole (1x | ref:1)
  - core/sek04_stale (2x | input:2)
  - core/sek07_predykcje (2x | input:2)
  - core/sek08_formalizm (2x | input:2)
  - core/sek08c_metryka_z_substratu (2x | input:2)
  - partial_proofs/nuclear_from_soliton (3x | wikilink:3)
  - partial_proofs/trojcialowe_nbody (1x | input:1)

### research/neutrino_msw
- depended on by: --

### research/op-CORE-CLEANUP-B-2026-05-04
- depended on by: --

### research/op-D01-anchor-lock-2026-05-06
- depended on by:
  - <root> (1x | wikilink:1)
  - audyt/D01_drifting_numbers (1x | wikilink:1)
  - meta (1x | wikilink:1)
  - research/op-FRW-radiation-era-varying-c-2026-05-06 (9x | wikilink:9)
  - research/op-GWTC3-reanalysis (4x | wikilink:4)
  - research/op-LIGO-3G-deviation (3x | wikilink:3)
  - research/op-MAG-Lorentz-A-mu-coupling-2026-05-09 (2x | wikilink:2)
  - research/op-MAG-Phase5-V-reference-clarification-2026-05-09 (3x | wikilink:3)
  - research/op-MAG-anomalous-moment-2026-05-09 (1x | wikilink:1)
  - research/op-MAG-resonance-formalization-2026-05-09 (3x | wikilink:3)
  - research/op-Phi-decomposition-photon-2026-05-07 (7x | wikilink:7)
  - research/op-Phi-vacuum-scale-2026-05-09 (6x | wikilink:6)
  - research/op-S07-alternative-f-psi-derivation-2026-05-09 (2x | wikilink:2)
  - research/op-SPIN-SU2-substrate-derivation-2026-05-08 (2x | wikilink:2)
  - research/op-V-canonical-consistency-audit-2026-05-09 (3x | wikilink:3)
  - research/op-dual-V-structure-clarification-2026-05-09 (3x | wikilink:3)
  - research/op-emergent-metric-from-interaction-2026-05-09 (2x | wikilink:2)
  - research/op-ppE-mapping (3x | wikilink:3)
  - research/op-ppE-mapping/TGP/TGP_v1/research/op-LIGO-3G-deviation (3x | wikilink:3)
  - research/op-ppE-mapping/TGP/TGP_v1/research/op-ppE-mapping (1x | wikilink:1)

### research/op-FRW-radiation-era-varying-c-2026-05-06
- depended on by: --

### research/op-GWTC3-reanalysis
- depended on by:
  - <root> (6x | wikilink:6)
  - audyt/T01_LIGO3G_falsifier (9x | wikilink:9)
  - papers/M911_LIGO3G_paper (1x | wikilink:1)
  - research/op-LIGO-3G-deviation (1x | wikilink:1)
  - research/op-Phi-vacuum-scale-2026-05-09 (1x | wikilink:1)
  - research/op-S07-alternative-f-psi-derivation-2026-05-09 (5x | wikilink:5)
  - research/op-emergent-metric-from-interaction-2026-05-09 (3x | wikilink:3)
  - research/op-ppE-mapping (5x | wikilink:5)

### research/op-GWTC3-reanalysis/scripts
- depended on by:
  - research/op-GWTC3-reanalysis (11x | wikilink:11)

### research/op-L01-rho-stress-energy-bridge-2026-05-04
- depended on by:
  - audyt (3x | wikilink:3)
  - audyt/L01_rho_operational (3x | wikilink:3)

### research/op-L03-spectral-stability-2026-05-06
- depended on by:
  - audyt/L03_K_phi_stability (3x | wikilink:3)

### research/op-L04-ODE-canonicalization-2026-05-04
- depended on by:
  - audyt (4x | wikilink:4)
  - audyt/L04_ODE_dualism_alpha (4x | wikilink:4)

### research/op-LIGO-3G-deviation
- depended on by:
  - audyt/T01_LIGO3G_falsifier (7x | wikilink:7)
  - papers/M911_LIGO3G_paper (1x | wikilink:1)
  - research/op-GWTC3-reanalysis (10x | wikilink:10)
  - research/op-ppE-mapping (3x | wikilink:3)
  - research/op-ppE-mapping/TGP/TGP_v1/papers/M911_LIGO3G_paper (1x | wikilink:1)
  - research/op-ppE-mapping/TGP/TGP_v1/research/op-LIGO-3G-deviation (7x | wikilink:7)

### research/op-LIGO-3G-deviation/scripts
- depended on by:
  - research/op-LIGO-3G-deviation (5x | wikilink:5)
  - research/op-ppE-mapping/TGP/TGP_v1/research/op-LIGO-3G-deviation (5x | wikilink:5)

### research/op-M03-balance-sheet-retrofit-2026-05-06
- depended on by:
  - <root> (1x | wikilink:1)
  - audyt (1x | wikilink:1)
  - audyt/M03_balance_sheet_missing (14x | wikilink:14)
  - meta (7x | wikilink:7)
  - meta/research (2x | wikilink:2)
  - research/op-FRW-radiation-era-varying-c-2026-05-06 (4x | wikilink:4)
  - research/op-Phi-decomposition-photon-2026-05-07 (2x | wikilink:2)

### research/op-MAG-Lorentz-A-mu-coupling-2026-05-09
- depended on by:
  - research/op-MAG-anomalous-moment-2026-05-09 (2x | wikilink:2)
  - research/op-MAG-resonance-formalization-2026-05-09 (4x | wikilink:4)
  - research/op-Phi-vacuum-scale-2026-05-09 (3x | wikilink:3)
  - research/op-SPIN-SU2-substrate-derivation-2026-05-08 (1x | wikilink:1)
  - research/op-emergent-metric-from-interaction-2026-05-09 (1x | wikilink:1)

### research/op-MAG-Phase5-V-reference-clarification-2026-05-09
- depended on by:
  - research/op-V-canonical-consistency-audit-2026-05-09 (1x | wikilink:1)
  - research/op-dual-V-structure-clarification-2026-05-09 (5x | wikilink:5)

### research/op-MAG-anomalous-moment-2026-05-09
- depended on by: --

### research/op-MAG-resonance-formalization-2026-05-09
- depended on by:
  - <root> (1x | wikilink:1)
  - research/op-MAG-Lorentz-A-mu-coupling-2026-05-09 (2x | wikilink:2)
  - research/op-MAG-Phase5-V-reference-clarification-2026-05-09 (10x | wikilink:10)
  - research/op-Phase5-MAG-erratum-2026-05-09 (7x | wikilink:7)
  - research/op-Phi-vacuum-scale-2026-05-09 (12x | wikilink:12)
  - research/op-SPIN-SU2-substrate-derivation-2026-05-08 (1x | wikilink:1)
  - research/op-V-canonical-consistency-audit-2026-05-09 (5x | wikilink:5)
  - research/op-dual-V-structure-clarification-2026-05-09 (2x | wikilink:2)

### research/op-Phase5-MAG-erratum-2026-05-09
- depended on by:
  - <root> (1x | wikilink:1)
  - research/op-MAG-resonance-formalization-2026-05-09 (2x | wikilink:2)

### research/op-Phi-decomposition-photon-2026-05-07
- depended on by: --

### research/op-Phi-vacuum-scale-2026-05-09
- depended on by:
  - <root> (2x | wikilink:2)
  - research/op-Phi0-spatial-variation-predictions-2026-05-09 (4x | wikilink:4)
  - research/op-S07-alternative-f-psi-derivation-2026-05-09 (2x | wikilink:2)
  - research/op-V-canonical-consistency-audit-2026-05-09 (7x | wikilink:7)

### research/op-Phi0-spatial-variation-predictions-2026-05-09
- depended on by: --

### research/op-S07-alternative-f-psi-derivation-2026-05-09
- depended on by: --

### research/op-SPIN-MAG-leakage-lean-hypothesis-2026-05-07
- depended on by: --

### research/op-SPIN-SU2-substrate-derivation-2026-05-08
- depended on by:
  - research/op-MAG-anomalous-moment-2026-05-09 (4x | wikilink:4)
  - research/op-MAG-resonance-formalization-2026-05-09 (11x | wikilink:11)
  - research/op-emergent-metric-from-interaction-2026-05-09 (3x | wikilink:3)

### research/op-V-canonical-consistency-audit-2026-05-09
- depended on by:
  - research/op-MAG-Phase5-V-reference-clarification-2026-05-09 (4x | wikilink:4)

### research/op-alpha-fine-structure
- depended on by:
  - <root> (73x | wikilink:73)
  - audyt/S01_metric_four_forms (1x | wikilink:1)
  - audyt/S02_volume_element_M9 (1x | wikilink:1)
  - audyt/S03_beta_PPN_convention (1x | wikilink:1)
  - audyt/T01_LIGO3G_falsifier (12x | wikilink:12)
  - meta (5x | wikilink:5)
  - meta/core/intake (10x | wikilink:10)
  - meta/research (2x | wikilink:2)
  - meta/research/_examples (4x | wikilink:4)
  - meta/research/templates (1x | wikilink:1)
  - research/op-FRW-radiation-era-varying-c-2026-05-06 (5x | wikilink:5)
  - research/op-GWTC3-reanalysis (6x | wikilink:6)
  - research/op-L04-ODE-canonicalization-2026-05-04 (1x | wikilink:1)
  - research/op-LIGO-3G-deviation (13x | wikilink:13)
  - research/op-M03-balance-sheet-retrofit-2026-05-06 (107x | wikilink:107)
  - research/op-MAG-Lorentz-A-mu-coupling-2026-05-09 (1x | wikilink:1)
  - research/op-MAG-resonance-formalization-2026-05-09 (6x | wikilink:6)
  - research/op-Phase5-MAG-erratum-2026-05-09 (4x | wikilink:4)
  - research/op-Phi-decomposition-photon-2026-05-07 (7x | wikilink:7)
  - research/op-Phi-vacuum-scale-2026-05-09 (9x | wikilink:9)
  - research/op-Phi0-spatial-variation-predictions-2026-05-09 (1x | wikilink:1)
  - research/op-S07-alternative-f-psi-derivation-2026-05-09 (2x | wikilink:2)
  - research/op-SPIN-MAG-leakage-lean-hypothesis-2026-05-07 (4x | wikilink:4)
  - research/op-SPIN-SU2-substrate-derivation-2026-05-08 (5x | wikilink:5)
  - research/op-V-canonical-consistency-audit-2026-05-09 (1x | wikilink:1)
  - research/op-bh-alpha-threshold (10x | wikilink:10)
  - research/op-chi1-newton-constant-derivation (52x | wikilink:52)
  - research/op-cross-sector-charge (14x | wikilink:14)
  - research/op-dual-V-structure-clarification-2026-05-09 (8x | wikilink:8)
  - research/op-emergent-metric-from-interaction-2026-05-09 (2x | wikilink:2)
  - research/op-eps-photon-ring (12x | wikilink:12)
  - research/op-eta-wolfenstein (13x | wikilink:13)
  - research/op-eta2-denom-derivation (16x | wikilink:16)
  - research/op-g0-r3-from-canonical-projection (14x | wikilink:14)
  - research/op-iota-charge-pmns-unification (40x | wikilink:40)
  - research/op-kappa-mixing-numerator (35x | wikilink:35)
  - research/op-lambda1-e2-amplitude-emergence (10x | wikilink:10)
  - research/op-mu-pmns-phase-hardening (43x | wikilink:43)
  - research/op-nu-majorana-phase-mbb (40x | wikilink:40)
  - research/op-omega1-substrate-em-coupling (33x | wikilink:33)
  - research/op-omega2-axion-coupling-lock (51x | wikilink:51)
  - research/op-omega3-axion-decay-constant (45x | wikilink:45)
  - research/op-omicron1-sigmamnu-cosmo (38x | wikilink:38)
  - research/op-omicron2-phi-mean-shift-cosmo (1x | wikilink:1)
  - research/op-phi1-substrate-action-variational (34x | wikilink:34)
  - research/op-pi1-bb0nu-nme-isotope (39x | wikilink:39)
  - research/op-ppE-mapping (14x | wikilink:14)
  - research/op-ppE-mapping/TGP/TGP_v1/research/op-LIGO-3G-deviation (13x | wikilink:13)
  - research/op-ppE-mapping/TGP/TGP_v1/research/op-ppE-mapping (3x | wikilink:3)
  - research/op-psi1-substrate-light-acceleration (13x | wikilink:13)
  - research/op-rho1-71Ge-cross-section (42x | wikilink:42)
  - research/op-sc-alpha-origin (13x | wikilink:13)
  - research/op-sigma1-substrate-light-dispersion (26x | wikilink:26)
  - research/op-tau1-closure-overlap-coulomb (35x | wikilink:35)
  - research/op-tau2-substrate-time-coupling (10x | wikilink:10)
  - research/op-tau3-substrate-clock-acceleration (22x | wikilink:22)
  - research/op-tensor-modes-Phi-FUTURE (3x | wikilink:3)
  - research/op-theta-quark-koide (12x | wikilink:12)
  - research/op-upsilon1-closure-cross-family (36x | wikilink:36)
  - research/op-uv-as-ngfp (15x | wikilink:15)
  - research/op-uv2-mtgp-absolute-scale (49x | wikilink:49)
  - research/op-uv3-phi0-renormalization (12x | wikilink:12)
  - research/op-xi-photon-ring (14x | wikilink:14)
  - research/op-xi2-sterile-nu-5sector (39x | wikilink:39)
  - research/op-zeta-mass-spectrum (12x | wikilink:12)

### research/op-bh-alpha-threshold
- depended on by: --

### research/op-chi1-newton-constant-derivation
- depended on by:
  - <root> (3x | wikilink:3)
  - audyt/S06_circular_anchors (2x | wikilink:2)
  - meta (3x | wikilink:3)
  - research/op-M03-balance-sheet-retrofit-2026-05-06 (1x | wikilink:1)
  - research/op-omega2-axion-coupling-lock (1x | wikilink:1)
  - research/op-omega3-axion-decay-constant (1x | wikilink:1)
  - research/op-uv2-mtgp-absolute-scale (1x | wikilink:1)

### research/op-cosmology-closure
- depended on by:
  - research/audyt_cosmology_drift_2026-05-03 (3x | wikilink:3)
  - research/closure_2026-04-26 (7x | wikilink:7)
  - research/cosmo_tensions (1x | wikilink:1)
  - research/desi_dark_energy (3x | wikilink:3)
  - research/hubble_tension (2x | wikilink:2)
  - research/op-FRW-radiation-era-varying-c-2026-05-06 (4x | wikilink:4)
  - research/op-M03-balance-sheet-retrofit-2026-05-06 (9x | wikilink:9)
  - research/op-omicron2-phi-mean-shift-cosmo (7x | wikilink:7)
  - research/op-phase1-covariant (2x | wikilink:2)
  - research/op-phase2-quantum-gravity (2x | wikilink:2)
  - research/op-phase3-uv-completion (1x | wikilink:1)
  - research/op-quantum-closure (9x | wikilink:9)
  - research/op-void-flat-modes-h0-2026-05-06 (9x | wikilink:9)

### research/op-cross-sector-charge
- depended on by: --

### research/op-delta1-g-tilde-derivation
- depended on by:
  - research/op-lambda1-e2-amplitude-emergence (2x | wikilink:2)

### research/op-delta2-Nf-derivation
- depended on by:
  - research/op-lambda1-e2-amplitude-emergence (1x | wikilink:1)

### research/op-dual-V-structure-clarification-2026-05-09
- depended on by:
  - <root> (1x | wikilink:1)

### research/op-eht
- depended on by:
  - research/op-M03-balance-sheet-retrofit-2026-05-06 (10x | wikilink:10)
  - research/op-eht-A (4x | wikilink:4)
  - research/op-m92 (1x | wikilink:1)

### research/op-eht-A
- depended on by:
  - research/op-M03-balance-sheet-retrofit-2026-05-06 (7x | wikilink:7)
  - research/op-m92 (3x | wikilink:3)

### research/op-emergent-metric-from-interaction-2026-05-09
- depended on by:
  - research/op-Phi-decomposition-photon-2026-05-07 (6x | wikilink:6)

### research/op-eps-photon-ring
- depended on by: --

### research/op-eta-wolfenstein
- depended on by: --

### research/op-eta2-denom-derivation
- depended on by: --

### research/op-g0-r3-from-canonical-projection
- depended on by:
  - audyt/S01_metric_four_forms (2x | wikilink:2)
  - research/op-S07-alternative-f-psi-derivation-2026-05-09 (3x | wikilink:3)
  - research/op-V-canonical-consistency-audit-2026-05-09 (4x | wikilink:4)
  - research/op-dual-V-structure-clarification-2026-05-09 (5x | wikilink:5)

### research/op-gamma1-phi-eff-anchor-resolution
- depended on by:
  - research/op-lambda1-e2-amplitude-emergence (2x | wikilink:2)

### research/op-iota-charge-pmns-unification
- depended on by: --

### research/op-kappa-mixing-numerator
- depended on by: --

### research/op-lambda1-e2-amplitude-emergence
- depended on by:
  - research/op-M03-balance-sheet-retrofit-2026-05-06 (2x | wikilink:2)
  - research/op-delta1-g-tilde-derivation (2x | wikilink:2)
  - research/op-gamma1-phi-eff-anchor-resolution (4x | wikilink:4)
  - research/op-mu1-minimal-substrate-log-redefinition (5x | wikilink:5)

### research/op-m92
- depended on by:
  - meta/core/intake (1x | wikilink:1)
  - research/closure_2026-04-26 (1x | wikilink:1)
  - research/closure_2026-04-26/alpha_psi_threshold (8x | wikilink:8)
  - research/op-bh-alpha-threshold (3x | wikilink:3)
  - research/op-newton-momentum (2x | wikilink:2)

### research/op-mu-pmns-phase-hardening
- depended on by: --

### research/op-mu1-minimal-substrate-log-redefinition
- depended on by: --

### research/op-newton-momentum
- depended on by:
  - <root> (2x | wikilink:2)
  - audyt (2x | wikilink:2)
  - audyt/D01_drifting_numbers (2x | wikilink:2)
  - audyt/S01_metric_four_forms (1x | wikilink:1)
  - audyt/S04_metric_coupling_axiom (3x | wikilink:3)
  - audyt/T01_LIGO3G_falsifier (23x | wikilink:23)
  - meta (9x | wikilink:9)
  - meta/core/intake (1x | wikilink:1)
  - papers/M911_LIGO3G_paper (1x | wikilink:1)
  - research/closure_2026-04-26 (2x | wikilink:2)
  - research/closure_2026-04-26/Lambda_from_Phi0 (3x | wikilink:3)
  - research/closure_2026-04-26/f_psi_principle (6x | wikilink:6)
  - research/op-D01-anchor-lock-2026-05-06 (7x | wikilink:7)
  - research/op-FRW-radiation-era-varying-c-2026-05-06 (2x | wikilink:2)
  - research/op-L01-rho-stress-energy-bridge-2026-05-04 (6x | wikilink:6)
  - research/op-M03-balance-sheet-retrofit-2026-05-06 (13x | wikilink:13)
  - research/op-cosmology-closure (11x | wikilink:11)
  - research/op-eht (4x | wikilink:4)
  - research/op-eht-A (2x | wikilink:2)
  - research/op-m92 (1x | wikilink:1)
  - research/op-phase1-covariant (4x | wikilink:4)
  - research/op-phase2-quantum-gravity (3x | wikilink:3)
  - research/op-phase3-uv-completion (3x | wikilink:3)
  - research/op-ppE-mapping (9x | wikilink:9)
  - research/op-ppE-mapping/TGP/TGP_v1/papers/M911_LIGO3G_paper (1x | wikilink:1)
  - research/op-ppE-mapping/TGP/TGP_v1/research/op-ppE-mapping (1x | wikilink:1)
  - research/op-quantum-closure (7x | wikilink:7)
  - research/op7 (2x | wikilink:2)

### research/op-nu-majorana-phase-mbb
- depended on by: --

### research/op-omega1-substrate-em-coupling
- depended on by: --

### research/op-omega2-axion-coupling-lock
- depended on by:
  - <root> (2x | wikilink:2)
  - audyt (1x | wikilink:1)
  - audyt/S06_circular_anchors (3x | wikilink:3)
  - meta (2x | wikilink:2)
  - research/op-omega3-axion-decay-constant (1x | wikilink:1)

### research/op-omega3-axion-decay-constant
- depended on by:
  - <root> (2x | wikilink:2)
  - audyt (1x | wikilink:1)
  - audyt/S06_circular_anchors (3x | wikilink:3)
  - meta (2x | wikilink:2)
  - research/op-omega2-axion-coupling-lock (1x | wikilink:1)

### research/op-omicron1-sigmamnu-cosmo
- depended on by: --

### research/op-omicron2-phi-mean-shift-cosmo
- depended on by:
  - research/audyt_cosmology_drift_2026-05-03 (1x | wikilink:1)
  - research/op-M03-balance-sheet-retrofit-2026-05-06 (5x | wikilink:5)
  - research/op-void-flat-modes-h0-2026-05-06 (4x | wikilink:4)

### research/op-phase1-covariant
- depended on by:
  - research/closure_2026-04-26 (11x | wikilink:11)
  - research/op-M03-balance-sheet-retrofit-2026-05-06 (4x | wikilink:4)
  - research/op-phase2-quantum-gravity (16x | wikilink:16)
  - research/op-phase3-uv-completion (4x | wikilink:4)

### research/op-phase2-quantum-gravity
- depended on by:
  - <root> (1x | wikilink:1)
  - research/closure_2026-04-26 (12x | wikilink:12)
  - research/op-M03-balance-sheet-retrofit-2026-05-06 (4x | wikilink:4)
  - research/op-chi1-newton-constant-derivation (6x | wikilink:6)
  - research/op-phase3-uv-completion (10x | wikilink:10)
  - research/op-uv2-mtgp-absolute-scale (1x | wikilink:1)

### research/op-phase3-uv-completion
- depended on by:
  - <root> (2x | wikilink:2)
  - meta/research/_examples (2x | wikilink:2)
  - research/closure_2026-04-26 (19x | wikilink:19)
  - research/op-M03-balance-sheet-retrofit-2026-05-06 (5x | wikilink:5)
  - research/op-uv-renormalizability-research (6x | wikilink:6)

### research/op-phi1-substrate-action-variational
- depended on by: --

### research/op-pi1-bb0nu-nme-isotope
- depended on by: --

### research/op-ppE-mapping
- depended on by:
  - <root> (6x | wikilink:6)
  - audyt/T01_LIGO3G_falsifier (2x | wikilink:2)
  - papers/M911_LIGO3G_paper (2x | wikilink:2)
  - research/op-GWTC3-reanalysis (6x | wikilink:6)
  - research/op-LIGO-3G-deviation (4x | wikilink:4)
  - research/op-Phi-vacuum-scale-2026-05-09 (1x | wikilink:1)
  - research/op-S07-alternative-f-psi-derivation-2026-05-09 (5x | wikilink:5)
  - research/op-emergent-metric-from-interaction-2026-05-09 (4x | wikilink:4)
  - research/op-ppE-mapping/TGP/TGP_v1/papers/M911_LIGO3G_paper (1x | wikilink:1)
  - research/op-ppE-mapping/TGP/TGP_v1/research/op-LIGO-3G-deviation (3x | wikilink:3)
  - research/op-ppE-mapping/TGP/TGP_v1/research/op-ppE-mapping (4x | wikilink:4)

### research/op-ppE-mapping/TGP/TGP_v1/audyt/T01_LIGO3G_falsifier
- depended on by: --

### research/op-ppE-mapping/TGP/TGP_v1/papers/M911_LIGO3G_paper
- depended on by: --

### research/op-ppE-mapping/TGP/TGP_v1/research/op-LIGO-3G-deviation
- depended on by: --

### research/op-ppE-mapping/TGP/TGP_v1/research/op-ppE-mapping
- depended on by: --

### research/op-ppE-mapping/scripts
- depended on by:
  - <root> (2x | wikilink:2)
  - research/op-GWTC3-reanalysis (1x | wikilink:1)
  - research/op-ppE-mapping (7x | wikilink:7)
  - research/op-ppE-mapping/TGP/TGP_v1/research/op-ppE-mapping (2x | wikilink:2)

### research/op-psi1-substrate-light-acceleration
- depended on by:
  - <root> (3x | wikilink:3)
  - audyt/S01_metric_four_forms (1x | wikilink:1)
  - research/op-L01-rho-stress-energy-bridge-2026-05-04 (1x | wikilink:1)
  - research/op-M03-balance-sheet-retrofit-2026-05-06 (7x | wikilink:7)

### research/op-quantum-closure
- depended on by:
  - research/closure_2026-04-26 (5x | wikilink:5)
  - research/op-M03-balance-sheet-retrofit-2026-05-06 (8x | wikilink:8)
  - research/op-phase1-covariant (14x | wikilink:14)
  - research/op-phase2-quantum-gravity (3x | wikilink:3)
  - research/op-phase3-uv-completion (1x | wikilink:1)

### research/op-rho1-71Ge-cross-section
- depended on by: --

### research/op-sc-alpha-origin
- depended on by: --

### research/op-sigma1-substrate-light-dispersion
- depended on by: --

### research/op-tau1-closure-overlap-coulomb
- depended on by: --

### research/op-tau2-substrate-time-coupling
- depended on by: --

### research/op-tau3-substrate-clock-acceleration
- depended on by:
  - research/op-L01-rho-stress-energy-bridge-2026-05-04 (1x | wikilink:1)

### research/op-tensor-modes-Phi-FUTURE
- depended on by: --

### research/op-theta-quark-koide
- depended on by: --

### research/op-upsilon1-closure-cross-family
- depended on by: --

### research/op-uv-as-ngfp
- depended on by: --

### research/op-uv-renormalizability-research
- depended on by: --

### research/op-uv2-mtgp-absolute-scale
- depended on by:
  - <root> (2x | wikilink:2)
  - audyt/S06_circular_anchors (2x | wikilink:2)
  - meta (3x | wikilink:3)
  - research/op-M03-balance-sheet-retrofit-2026-05-06 (4x | wikilink:4)
  - research/op-omega2-axion-coupling-lock (1x | wikilink:1)
  - research/op-omega3-axion-decay-constant (1x | wikilink:1)
  - research/op-uv3-phi0-renormalization (4x | wikilink:4)

### research/op-uv3-phi0-renormalization
- depended on by: --

### research/op-void-flat-modes-h0-2026-05-06
- depended on by: --

### research/op-xi-photon-ring
- depended on by: --

### research/op-xi2-sterile-nu-5sector
- depended on by: --

### research/op-zeta-mass-spectrum
- depended on by: --

### research/op1-op2-op4
- depended on by:
  - research/external_review_2026-04-25 (2x | wikilink:2)
  - research/op-phase1-covariant (2x | wikilink:2)
  - research/op-quantum-closure (21x | wikilink:21)

### research/op6
- depended on by: --

### research/op7
- depended on by:
  - meta/research/_examples (1x | wikilink:1)
  - research/closure_2026-04-26 (7x | wikilink:7)
  - research/closure_2026-04-26/Lambda_from_Phi0 (3x | wikilink:3)
  - research/closure_2026-04-26/alpha_psi_threshold (2x | wikilink:2)
  - research/closure_2026-04-26/f_psi_principle (3x | wikilink:3)
  - research/closure_2026-04-26/sigma_ab_pathB (8x | wikilink:8)
  - research/external_review_2026-04-25 (1x | wikilink:1)
  - research/op-M03-balance-sheet-retrofit-2026-05-06 (1x | wikilink:1)
  - research/op-eht (10x | wikilink:10)
  - research/op-eht-A (1x | wikilink:1)
  - research/op-m92 (2x | wikilink:2)

### research/particle_sector_closure
- depended on by: --

### research/qm_born_rule
- depended on by: --

### research/qm_decoherence
- depended on by: --

### research/qm_entanglement
- depended on by: --

### research/qm_foundations
- depended on by: --

### research/qm_measurement
- depended on by: --

### research/qm_spin
- depended on by: --

### research/qm_statistics
- depended on by:
  - partial_proofs/fermion_from_soliton (2x | wikilink:2)

### research/qm_superposition
- depended on by: --

### research/rho_normal_state_closure
- depended on by:
  - research/superconductivity_closure (28x | wikilink:28)

### research/s8_tension
- depended on by: --

### research/superconductivity_closure
- depended on by:
  - research/atomic_shells_closure (2x | wikilink:2)
  - research/casimir_mof (1x | wikilink:1)
  - research/liquid_viscosity (1x | wikilink:1)
  - research/muon_g_minus_2 (1x | wikilink:1)

### research/thermal_transport_molecular
- depended on by: --

### research/uv_completion
- depended on by: --

### research/why_n3
- depended on by:
  - audyt (2x | wikilink:2)
  - audyt/L04_ODE_dualism_alpha (4x | wikilink:4)
  - audyt/L05_mass_exponent_drift (4x | wikilink:4)
  - audyt/T01_LIGO3G_falsifier (1x | wikilink:1)
  - research/mass_scaling_k4 (4x | wikilink:4)
  - research/op-L04-ODE-canonicalization-2026-05-04 (49x | wikilink:49)
  - research/op-g0-r3-from-canonical-projection (2x | wikilink:2)
  - research/op-gamma1-phi-eff-anchor-resolution (1x | wikilink:1)
  - research/op-lambda1-e2-amplitude-emergence (3x | wikilink:3)

### tooling
- depended on by: --

### tooling/scripts
- depended on by:
  - research/em_from_substrate (1x | wikilink:1)
  - research/op-quantum-closure (1x | wikilink:1)

### tooling/scripts/gauge
- depended on by:
  - research/em_from_substrate (2x | wikilink:2)

### tooling/scripts/gravity
- depended on by:
  - research/op-eht (1x | wikilink:1)

### tooling/scripts/substrate
- depended on by:
  - research/op7 (1x | wikilink:1)

