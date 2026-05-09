# Dependencies - TGP/TGP_v1 (forward graph)

> Generated: 2026-05-09 by `tooling/build_deps_graph.py`
> Sources: \input{} + \ref{} + \cite{} + [[wikilink]]

## Summary

- Folders analyzed (fine): 213
- Folders analyzed (coarse): 12
- Total dependencies found: 7430
  - `\input`  edges: 70
  - `\ref`    edges: 1469
  - `[[wiki]]` edges: 5891
- `\cite{}` usages counted (bib keys): 100
- Bibliography keys in tgp_main.bib: 25

- Orphans (unresolved): \input=0, \ref=8, wikilink=655

## Coarse view (top-level folders)

### <root>
- depends on: audyt, axioms, core, core/_meta_latex, core/formalizm, meta, papers_external, partial_proofs, research

### audyt
- depends on: <root>, axioms, core, core/_meta_latex, meta, papers, research

### axioms
- depends on: core, core/formalizm, partial_proofs

### core
- depends on: axioms, core/_meta_latex, core/formalizm, partial_proofs, research

### core/_meta_latex
- depends on: axioms, core, core/formalizm, partial_proofs

### core/formalizm
- depends on: axioms, core, partial_proofs

### meta
- depends on: <root>, audyt, axioms, research

### papers
- depends on: <root>, research

### papers_external
- depends on: axioms, core, core/formalizm

### partial_proofs
- depends on: axioms, core, core/formalizm, papers_external, research

### research
- depends on: <root>, audyt, axioms, core, core/_meta_latex, core/formalizm, meta, papers, papers_external, partial_proofs, tooling

### tooling
- depends on: --

## Fine view (per subfolder)

### <root>
- `\input`:
  - axioms/notacja (2x)
  - axioms/roznica_N0 (1x)
  - axioms/substrat (1x)
  - core/_meta_latex (2x)
  - core/formalizm (11x)
  - core/sek00_summary (1x)
  - core/sek01_ontologia (1x)
  - core/sek02_pole (1x)
  - core/sek03_rezimy (1x)
  - core/sek04_stale (1x)
  - core/sek05_ciemna_energia (1x)
  - core/sek06_czarne_dziury (1x)
  - core/sek07_predykcje (1x)
  - core/sek07a_wymiar_wzmocniony (1x)
  - core/sek08_formalizm (1x)
  - core/sek08a_akcja_zunifikowana (1x)
  - core/sek08b_ghost_resolution (1x)
  - core/sek08c_metryka_z_substratu (1x)
  - core/sek09_cechowanie (1x)
  - core/sek10_N0_wyprowadzenie (1x)
  - partial_proofs/alphaK_status (1x)
  - partial_proofs/bh_ringdown (2x)
  - partial_proofs/bounce_topology (1x)
  - partial_proofs/chiralnosc (1x)
  - partial_proofs/defect_hierarchy (1x)
  - partial_proofs/hierarchia_mas (4x)
  - partial_proofs/koide_fp (5x)
  - partial_proofs/most_gamma_phi (2x)
  - partial_proofs/particle_sector (2x)
  - partial_proofs/potencjal (1x)
  - partial_proofs/quark_sector (1x)
  - partial_proofs/superconductivity (2x)
  - partial_proofs/trojcialowe_nbody (2x)
  - partial_proofs/wielki_wybuch (1x)
  - partial_proofs/zero_mode (2x)
- `\ref` to other folders:
  - core/sek08_formalizm (3x)
  - papers_external/paper_bh_shadow (2x)
  - papers_external/tgp_sc_paper (2x)
  - partial_proofs/hierarchia_mas (2x)
- Wikilinks to:
  - audyt (4x)
  - audyt/D01_drifting_numbers (5x)
  - audyt/M03_balance_sheet_missing (1x)
  - audyt/T01_LIGO3G_falsifier (5x)
  - meta (8x)
  - meta/research (1x)
  - research (21x)
  - research/closure_2026-04-26 (1x)
  - research/op-D01-anchor-lock-2026-05-06 (1x)
  - research/op-GWTC3-reanalysis (6x)
  - research/op-M03-balance-sheet-retrofit-2026-05-06 (1x)
  - research/op-MAG-resonance-formalization-2026-05-09 (1x)
  - research/op-Phase5-MAG-erratum-2026-05-09 (1x)
  - research/op-Phi-vacuum-scale-2026-05-09 (2x)
  - research/op-alpha-fine-structure (73x)
  - research/op-chi1-newton-constant-derivation (3x)
  - research/op-dual-V-structure-clarification-2026-05-09 (1x)
  - research/op-newton-momentum (2x)
  - research/op-omega2-axion-coupling-lock (2x)
  - research/op-omega3-axion-decay-constant (2x)
  - research/op-phase2-quantum-gravity (1x)
  - research/op-phase3-uv-completion (2x)
  - research/op-ppE-mapping (6x)
  - research/op-ppE-mapping/scripts (2x)
  - research/op-psi1-substrate-light-acceleration (3x)
  - research/op-uv2-mtgp-absolute-scale (2x)
- Cite count: 6

### audyt
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (12x)
  - audyt/D01_drifting_numbers (73x)
  - audyt/L01_rho_operational (13x)
  - audyt/M03_balance_sheet_missing (1x)
  - audyt/T01_LIGO3G_falsifier (2x)
  - axioms/notacja (2x)
  - core/sek08b_ghost_resolution (1x)
  - core/sek08c_metryka_z_substratu (2x)
  - meta (37x)
  - meta/research (2x)
  - meta/research/templates (1x)
  - research/closure_2026-04-26/alpha_psi_threshold (1x)
  - research/op-L01-rho-stress-energy-bridge-2026-05-04 (3x)
  - research/op-L04-ODE-canonicalization-2026-05-04 (4x)
  - research/op-M03-balance-sheet-retrofit-2026-05-06 (1x)
  - research/op-newton-momentum (2x)
  - research/op-omega2-axion-coupling-lock (1x)
  - research/op-omega3-axion-decay-constant (1x)
  - research/why_n3 (2x)
- Cite count: 0

### audyt/D01_drifting_numbers
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - audyt (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta (5x)
  - research/op-D01-anchor-lock-2026-05-06 (1x)
  - research/op-newton-momentum (2x)
- Cite count: 0

### audyt/L01_rho_operational
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (2x)
  - audyt (7x)
  - audyt/D01_drifting_numbers (18x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - core/sek04_stale (2x)
  - core/sek08a_akcja_zunifikowana (4x)
  - meta (2x)
  - research/op-L01-rho-stress-energy-bridge-2026-05-04 (3x)
- Cite count: 0

### audyt/L02_beta_gamma_semantics
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (6x)
  - audyt (3x)
  - audyt/D01_drifting_numbers (5x)
  - axioms/notacja (3x)
  - meta (5x)
- Cite count: 0

### audyt/L03_K_phi_stability
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - audyt (5x)
  - audyt/D01_drifting_numbers (8x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - axioms/notacja (1x)
  - core/sek08b_ghost_resolution (2x)
  - meta (3x)
  - research/op-L03-spectral-stability-2026-05-06 (3x)
- Cite count: 0

### audyt/L04_ODE_dualism_alpha
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - audyt (3x)
  - audyt/D01_drifting_numbers (12x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - core/sek08_formalizm (2x)
  - meta (6x)
  - research/mass_scaling_k4 (2x)
  - research/op-L04-ODE-canonicalization-2026-05-04 (4x)
  - research/why_n3 (4x)
- Cite count: 0

### audyt/L05_mass_exponent_drift
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - audyt (3x)
  - audyt/D01_drifting_numbers (1x)
  - meta (4x)
  - research/why_n3 (4x)
- Cite count: 0

### audyt/L06_axion_mass_locked
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - audyt (3x)
  - audyt/D01_drifting_numbers (1x)
  - meta (3x)
- Cite count: 0

### audyt/L07_zero_sum_axiom
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (1x)
  - audyt (5x)
  - audyt/D01_drifting_numbers (3x)
  - core/sek01_ontologia (3x)
  - core/sek05_ciemna_energia (3x)
  - core/sek08_formalizm (2x)
- Cite count: 0

### audyt/L08_kink_fermion_closure
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (3x)
  - audyt (5x)
  - audyt/D01_drifting_numbers (7x)
  - core/sek08_formalizm (2x)
- Cite count: 0

### audyt/M01_status_creep
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (4x)
  - audyt (3x)
  - audyt/D01_drifting_numbers (1x)
  - meta (7x)
- Cite count: 0

### audyt/M02_ledger_pollution
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - audyt (3x)
  - audyt/D01_drifting_numbers (1x)
  - meta (9x)
- Cite count: 0

### audyt/M03_balance_sheet_missing
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - audyt (11x)
  - audyt/D01_drifting_numbers (17x)
  - core/_meta_latex (1x)
  - core/sek00_summary (1x)
  - meta (24x)
  - meta/research (3x)
  - research/op-M03-balance-sheet-retrofit-2026-05-06 (14x)
- Cite count: 0

### audyt/S01_metric_four_forms
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (2x)
  - audyt (4x)
  - audyt/D01_drifting_numbers (11x)
  - core/sek08c_metryka_z_substratu (2x)
  - meta (4x)
  - research/op-alpha-fine-structure (1x)
  - research/op-g0-r3-from-canonical-projection (2x)
  - research/op-newton-momentum (1x)
  - research/op-psi1-substrate-light-acceleration (1x)
- Cite count: 0

### audyt/S02_volume_element_M9
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - audyt (3x)
  - audyt/D01_drifting_numbers (7x)
  - audyt/L01_rho_operational (3x)
  - meta (3x)
  - research/op-alpha-fine-structure (1x)
- Cite count: 0

### audyt/S03_beta_PPN_convention
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - audyt (3x)
  - audyt/D01_drifting_numbers (7x)
  - audyt/L01_rho_operational (3x)
  - meta (3x)
  - research/op-alpha-fine-structure (1x)
- Cite count: 0

### audyt/S04_metric_coupling_axiom
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - audyt (3x)
  - audyt/D01_drifting_numbers (7x)
  - meta (3x)
  - research/op-newton-momentum (3x)
- Cite count: 0

### audyt/S05_tensor_sector_singleField
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (6x)
  - audyt (3x)
  - audyt/D01_drifting_numbers (9x)
  - meta (3x)
  - research/closure_2026-04-26 (6x)
  - research/closure_2026-04-26/alpha_psi_threshold (3x)
- Cite count: 0

### audyt/S06_circular_anchors
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - audyt (3x)
  - audyt/D01_drifting_numbers (7x)
  - meta (12x)
  - research/op-chi1-newton-constant-derivation (2x)
  - research/op-omega2-axion-coupling-lock (3x)
  - research/op-omega3-axion-decay-constant (3x)
  - research/op-uv2-mtgp-absolute-scale (2x)
- Cite count: 0

### audyt/S07_M911_derivation
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (1x)
  - audyt (5x)
  - audyt/D01_drifting_numbers (7x)
  - core/sek08a_akcja_zunifikowana (3x)
  - core/sek08c_metryka_z_substratu (3x)
- Cite count: 0

### audyt/T01_LIGO3G_falsifier
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (24x)
  - audyt (13x)
  - audyt/D01_drifting_numbers (87x)
  - core/sek08a_akcja_zunifikowana (1x)
  - core/sek08c_metryka_z_substratu (4x)
  - meta (3x)
  - papers/M911_LIGO3G_paper (5x)
  - research/op-GWTC3-reanalysis (9x)
  - research/op-LIGO-3G-deviation (7x)
  - research/op-alpha-fine-structure (12x)
  - research/op-newton-momentum (23x)
  - research/op-ppE-mapping (2x)
  - research/why_n3 (1x)
- Cite count: 0

### axioms
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - --
- Cite count: 0

### axioms/notacja
- `\input`:
  - --
- `\ref` to other folders:
  - core/sek01_ontologia (9x)
  - core/sek02_pole (1x)
  - core/sek03_rezimy (1x)
  - core/sek04_stale (5x)
  - core/sek05_ciemna_energia (1x)
  - core/sek06_czarne_dziury (2x)
  - core/sek08_formalizm (35x)
  - core/sek08a_akcja_zunifikowana (3x)
  - core/sek10_N0_wyprowadzenie (1x)
  - partial_proofs/alphaK_status (1x)
  - partial_proofs/bh_ringdown (1x)
  - partial_proofs/hierarchia_mas (1x)
  - partial_proofs/potencjal (1x)
- Wikilinks to:
  - --
- Cite count: 0

### axioms/roznica_N0
- `\input`:
  - --
- `\ref` to other folders:
  - axioms/notacja (2x)
  - axioms/substrat (2x)
  - core/formalizm (1x)
  - core/sek01_ontologia (5x)
  - core/sek08_formalizm (1x)
  - partial_proofs/zero_mode (5x)
- Wikilinks to:
  - --
- Cite count: 0

### axioms/substrat
- `\input`:
  - --
- `\ref` to other folders:
  - core/formalizm (4x)
  - core/sek01_ontologia (9x)
  - core/sek05_ciemna_energia (1x)
  - core/sek08_formalizm (15x)
- Wikilinks to:
  - --
- Cite count: 0

### core
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - --
- Cite count: 0

### core/_meta_latex
- `\input`:
  - --
- `\ref` to other folders:
  - axioms/notacja (2x)
  - axioms/substrat (1x)
  - core/formalizm (9x)
  - core/sek00_summary (1x)
  - core/sek03_rezimy (1x)
  - core/sek04_stale (7x)
  - core/sek05_ciemna_energia (1x)
  - core/sek07_predykcje (2x)
  - core/sek08_formalizm (40x)
  - core/sek08a_akcja_zunifikowana (6x)
  - core/sek08b_ghost_resolution (3x)
  - core/sek08c_metryka_z_substratu (1x)
  - core/sek09_cechowanie (12x)
  - core/sek10_N0_wyprowadzenie (1x)
  - partial_proofs/alphaK_status (12x)
  - partial_proofs/bh_ringdown (12x)
  - partial_proofs/defect_hierarchy (4x)
  - partial_proofs/hierarchia_mas (15x)
  - partial_proofs/koide_fp (4x)
  - partial_proofs/potencjal (1x)
  - partial_proofs/quark_sector (12x)
  - partial_proofs/wielki_wybuch (1x)
  - partial_proofs/zero_mode (4x)
- Wikilinks to:
  - --
- Cite count: 0

### core/formalizm
- `\input`:
  - --
- `\ref` to other folders:
  - axioms/substrat (13x)
  - core/sek00_summary (1x)
  - core/sek01_ontologia (8x)
  - core/sek02_pole (4x)
  - core/sek03_rezimy (11x)
  - core/sek04_stale (13x)
  - core/sek05_ciemna_energia (3x)
  - core/sek07_predykcje (2x)
  - core/sek07a_wymiar_wzmocniony (1x)
  - core/sek08_formalizm (22x)
  - core/sek08a_akcja_zunifikowana (4x)
  - core/sek08b_ghost_resolution (3x)
  - core/sek08c_metryka_z_substratu (3x)
  - core/sek09_cechowanie (51x)
  - core/sek10_N0_wyprowadzenie (2x)
  - partial_proofs/bh_ringdown (1x)
  - partial_proofs/chiralnosc (3x)
  - partial_proofs/hierarchia_mas (18x)
  - partial_proofs/koide_fp (4x)
  - partial_proofs/most_gamma_phi (2x)
  - partial_proofs/potencjal (3x)
  - partial_proofs/quark_sector (3x)
  - partial_proofs/wielki_wybuch (3x)
  - partial_proofs/zero_mode (14x)
- Wikilinks to:
  - --
- Cite count: 0

### core/sek00_summary
- `\input`:
  - core/_meta_latex (1x)
- `\ref` to other folders:
  - axioms/substrat (2x)
  - core/_meta_latex (2x)
  - core/formalizm (3x)
  - core/sek01_ontologia (4x)
  - core/sek02_pole (2x)
  - core/sek04_stale (6x)
  - core/sek05_ciemna_energia (1x)
  - core/sek06_czarne_dziury (1x)
  - core/sek07_predykcje (2x)
  - core/sek08_formalizm (13x)
  - core/sek08a_akcja_zunifikowana (3x)
  - core/sek08b_ghost_resolution (1x)
  - core/sek08c_metryka_z_substratu (2x)
  - core/sek09_cechowanie (1x)
  - core/sek10_N0_wyprowadzenie (4x)
  - partial_proofs/bh_ringdown (4x)
  - partial_proofs/hierarchia_mas (1x)
  - partial_proofs/wielki_wybuch (1x)
- Wikilinks to:
  - --
- Cite count: 0

### core/sek01_ontologia
- `\input`:
  - --
- `\ref` to other folders:
  - axioms/roznica_N0 (5x)
  - axioms/substrat (3x)
  - core/sek02_pole (2x)
  - core/sek04_stale (2x)
  - core/sek08_formalizm (21x)
  - core/sek09_cechowanie (1x)
- Wikilinks to:
  - --
- Cite count: 0

### core/sek02_pole
- `\input`:
  - --
- `\ref` to other folders:
  - core/sek01_ontologia (2x)
  - core/sek03_rezimy (2x)
  - core/sek04_stale (2x)
  - core/sek05_ciemna_energia (1x)
  - core/sek06_czarne_dziury (1x)
  - core/sek08_formalizm (9x)
  - core/sek08c_metryka_z_substratu (1x)
  - research/nbody/paper (1x)
- Wikilinks to:
  - --
- Cite count: 0

### core/sek03_rezimy
- `\input`:
  - --
- `\ref` to other folders:
  - core/formalizm (1x)
  - core/sek02_pole (2x)
  - core/sek06_czarne_dziury (1x)
  - core/sek08_formalizm (23x)
- Wikilinks to:
  - --
- Cite count: 2

### core/sek04_stale
- `\input`:
  - research/nbody/paper (2x)
- `\ref` to other folders:
  - core/formalizm (1x)
  - core/sek01_ontologia (1x)
  - core/sek08_formalizm (11x)
  - core/sek08c_metryka_z_substratu (1x)
  - partial_proofs/hierarchia_mas (1x)
- Wikilinks to:
  - --
- Cite count: 1

### core/sek05_ciemna_energia
- `\input`:
  - --
- `\ref` to other folders:
  - core/sek01_ontologia (1x)
  - core/sek07_predykcje (1x)
  - core/sek08_formalizm (11x)
- Wikilinks to:
  - --
- Cite count: 3

### core/sek06_czarne_dziury
- `\input`:
  - --
- `\ref` to other folders:
  - axioms/substrat (1x)
  - core/sek01_ontologia (2x)
  - core/sek04_stale (12x)
  - core/sek07_predykcje (1x)
  - core/sek08_formalizm (7x)
  - core/sek08a_akcja_zunifikowana (1x)
  - partial_proofs/bh_ringdown (1x)
- Wikilinks to:
  - --
- Cite count: 3

### core/sek07_predykcje
- `\input`:
  - research/nbody/paper (2x)
- `\ref` to other folders:
  - core/formalizm (14x)
  - core/sek01_ontologia (2x)
  - core/sek03_rezimy (6x)
  - core/sek04_stale (12x)
  - core/sek05_ciemna_energia (6x)
  - core/sek06_czarne_dziury (3x)
  - core/sek08_formalizm (53x)
  - core/sek08b_ghost_resolution (1x)
  - core/sek09_cechowanie (11x)
  - partial_proofs/bh_ringdown (6x)
  - partial_proofs/hierarchia_mas (3x)
  - partial_proofs/trojcialowe_nbody (1x)
- Wikilinks to:
  - --
- Cite count: 2

### core/sek07a_wymiar_wzmocniony
- `\input`:
  - --
- `\ref` to other folders:
  - core/formalizm (2x)
  - core/sek01_ontologia (1x)
  - core/sek03_rezimy (2x)
  - core/sek04_stale (1x)
  - core/sek07_predykcje (2x)
- Wikilinks to:
  - --
- Cite count: 0

### core/sek08_formalizm
- `\input`:
  - research/nbody/paper (2x)
- `\ref` to other folders:
  - axioms/substrat (26x)
  - core/_meta_latex (2x)
  - core/formalizm (10x)
  - core/sek01_ontologia (52x)
  - core/sek02_pole (5x)
  - core/sek03_rezimy (11x)
  - core/sek04_stale (24x)
  - core/sek05_ciemna_energia (4x)
  - core/sek06_czarne_dziury (5x)
  - core/sek07_predykcje (7x)
  - core/sek08a_akcja_zunifikowana (11x)
  - core/sek08c_metryka_z_substratu (1x)
  - core/sek09_cechowanie (2x)
  - partial_proofs/alphaK_status (5x)
  - partial_proofs/bh_ringdown (2x)
  - partial_proofs/hierarchia_mas (2x)
  - partial_proofs/potencjal (1x)
  - partial_proofs/trojcialowe_nbody (2x)
  - partial_proofs/wielki_wybuch (6x)
- Wikilinks to:
  - --
- Cite count: 0

### core/sek08a_akcja_zunifikowana
- `\input`:
  - --
- `\ref` to other folders:
  - core/formalizm (1x)
  - core/sek02_pole (4x)
  - core/sek08_formalizm (20x)
- Wikilinks to:
  - --
- Cite count: 0

### core/sek08b_ghost_resolution
- `\input`:
  - --
- `\ref` to other folders:
  - axioms/substrat (5x)
  - core/sek02_pole (1x)
  - core/sek08_formalizm (11x)
  - core/sek08a_akcja_zunifikowana (2x)
  - core/sek10_N0_wyprowadzenie (2x)
  - partial_proofs/hierarchia_mas (11x)
- Wikilinks to:
  - --
- Cite count: 0

### core/sek08c_metryka_z_substratu
- `\input`:
  - research/nbody/paper (2x)
- `\ref` to other folders:
  - axioms/substrat (2x)
  - core/formalizm (9x)
  - core/sek04_stale (3x)
  - core/sek08_formalizm (7x)
  - core/sek08a_akcja_zunifikowana (4x)
- Wikilinks to:
  - --
- Cite count: 0

### core/sek09_cechowanie
- `\input`:
  - --
- `\ref` to other folders:
  - axioms/substrat (1x)
  - core/formalizm (1x)
  - core/sek00_summary (1x)
  - core/sek01_ontologia (4x)
  - core/sek02_pole (4x)
  - core/sek03_rezimy (3x)
  - core/sek04_stale (1x)
  - core/sek07_predykcje (4x)
  - core/sek08_formalizm (10x)
  - core/sek08a_akcja_zunifikowana (2x)
  - core/sek08b_ghost_resolution (1x)
  - partial_proofs/chiralnosc (4x)
  - partial_proofs/defect_hierarchy (10x)
  - partial_proofs/hierarchia_mas (6x)
  - partial_proofs/koide_fp (1x)
  - partial_proofs/quark_sector (1x)
- Wikilinks to:
  - --
- Cite count: 0

### core/sek10_N0_wyprowadzenie
- `\input`:
  - --
- `\ref` to other folders:
  - axioms/substrat (3x)
  - core/formalizm (6x)
  - core/sek01_ontologia (1x)
  - core/sek08_formalizm (1x)
  - core/sek08b_ghost_resolution (3x)
- Wikilinks to:
  - --
- Cite count: 0

### meta
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (4x)
  - audyt/D01_drifting_numbers (3x)
  - meta/research (4x)
  - research/closure_2026-04-26 (5x)
  - research/op-D01-anchor-lock-2026-05-06 (1x)
  - research/op-M03-balance-sheet-retrofit-2026-05-06 (7x)
  - research/op-alpha-fine-structure (5x)
  - research/op-chi1-newton-constant-derivation (3x)
  - research/op-newton-momentum (9x)
  - research/op-omega2-axion-coupling-lock (2x)
  - research/op-omega3-axion-decay-constant (2x)
  - research/op-uv2-mtgp-absolute-scale (3x)
- Cite count: 0

### meta/core
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (1x)
  - audyt/D01_drifting_numbers (1x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta (4x)
  - meta/core/intake (4x)
  - meta/research (8x)
- Cite count: 0

### meta/core/intake
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - meta (1x)
  - meta/core (23x)
  - research/closure_2026-04-26 (2x)
  - research/closure_2026-04-26/alpha_psi_threshold (8x)
  - research/op-alpha-fine-structure (10x)
  - research/op-m92 (1x)
  - research/op-newton-momentum (1x)
- Cite count: 0

### meta/research
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (1x)
  - audyt/T01_LIGO3G_falsifier (70x)
  - axioms/notacja (2x)
  - meta (30x)
  - meta/core (7x)
  - meta/core/intake (9x)
  - meta/research/templates (6x)
  - research/closure_2026-04-26 (1x)
  - research/op-M03-balance-sheet-retrofit-2026-05-06 (2x)
  - research/op-alpha-fine-structure (2x)
- Cite count: 0

### meta/research/_examples
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (3x)
  - audyt/D01_drifting_numbers (6x)
  - meta (10x)
  - meta/core (1x)
  - meta/research (7x)
  - research/closure_2026-04-26 (3x)
  - research/closure_2026-04-26/alpha_psi_threshold (4x)
  - research/op-alpha-fine-structure (4x)
  - research/op-phase3-uv-completion (2x)
  - research/op7 (1x)
- Cite count: 0

### meta/research/templates
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (1x)
  - audyt/D01_drifting_numbers (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta (3x)
  - meta/research (6x)
  - research/op-alpha-fine-structure (1x)
- Cite count: 0

### papers/M911_LIGO3G_paper
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (1x)
  - research/op-GWTC3-reanalysis (1x)
  - research/op-LIGO-3G-deviation (1x)
  - research/op-newton-momentum (1x)
  - research/op-ppE-mapping (2x)
- Cite count: 0

### papers_external
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - --
- Cite count: 0

### papers_external/arxiv_submission
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - --
- Cite count: 3

### papers_external/paper_bh_shadow
- `\input`:
  - --
- `\ref` to other folders:
  - papers_external/tgp_english_summary (1x)
- Wikilinks to:
  - --
- Cite count: 24

### papers_external/paper_lepton_masses
- `\input`:
  - --
- `\ref` to other folders:
  - core/sek09_cechowanie (1x)
  - papers_external/paper_bh_shadow (1x)
- Wikilinks to:
  - --
- Cite count: 19

### papers_external/tgp_core_paper
- `\input`:
  - --
- `\ref` to other folders:
  - axioms/substrat (1x)
  - core/sek01_ontologia (7x)
  - core/sek08_formalizm (2x)
  - papers_external/paper_bh_shadow (6x)
- Wikilinks to:
  - --
- Cite count: 1

### papers_external/tgp_english_summary
- `\input`:
  - --
- `\ref` to other folders:
  - core/formalizm (1x)
  - core/sek01_ontologia (4x)
  - core/sek02_pole (5x)
  - core/sek03_rezimy (5x)
  - core/sek04_stale (7x)
  - core/sek05_ciemna_energia (2x)
  - core/sek08_formalizm (18x)
  - papers_external/paper_bh_shadow (3x)
  - papers_external/tgp_core_paper (2x)
- Wikilinks to:
  - --
- Cite count: 26

### papers_external/tgp_sc_paper
- `\input`:
  - --
- `\ref` to other folders:
  - papers_external/paper_bh_shadow (3x)
  - papers_external/paper_lepton_masses (3x)
- Wikilinks to:
  - --
- Cite count: 1

### partial_proofs
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - --
- Cite count: 0

### partial_proofs/alphaK_status
- `\input`:
  - --
- `\ref` to other folders:
  - core/sek08_formalizm (3x)
  - partial_proofs/hierarchia_mas (4x)
  - partial_proofs/most_gamma_phi (1x)
- Wikilinks to:
  - --
- Cite count: 0

### partial_proofs/bh_ringdown
- `\input`:
  - --
- `\ref` to other folders:
  - axioms/notacja (2x)
  - axioms/substrat (2x)
  - core/formalizm (1x)
  - core/sek04_stale (3x)
  - core/sek06_czarne_dziury (1x)
  - core/sek08_formalizm (10x)
  - partial_proofs/hierarchia_mas (14x)
- Wikilinks to:
  - --
- Cite count: 0

### partial_proofs/bounce_topology
- `\input`:
  - --
- `\ref` to other folders:
  - core/sek03_rezimy (1x)
- Wikilinks to:
  - --
- Cite count: 0

### partial_proofs/chiralnosc
- `\input`:
  - --
- `\ref` to other folders:
  - core/sek08b_ghost_resolution (1x)
  - partial_proofs/hierarchia_mas (1x)
- Wikilinks to:
  - --
- Cite count: 0

### partial_proofs/defect_hierarchy
- `\input`:
  - --
- `\ref` to other folders:
  - core/formalizm (2x)
  - core/sek02_pole (1x)
  - core/sek07_predykcje (1x)
  - core/sek08_formalizm (1x)
  - core/sek09_cechowanie (16x)
  - partial_proofs/chiralnosc (1x)
- Wikilinks to:
  - --
- Cite count: 0

### partial_proofs/fermion_from_soliton
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - core/formalizm (2x)
  - core/sek09_cechowanie (2x)
  - partial_proofs/defect_hierarchy (1x)
  - partial_proofs/nuclear_from_soliton (5x)
  - research/atom_from_soliton (2x)
  - research/nbody (1x)
  - research/qm_statistics (2x)
- Cite count: 0

### partial_proofs/hierarchia_mas
- `\input`:
  - --
- `\ref` to other folders:
  - axioms/notacja (2x)
  - axioms/substrat (1x)
  - core/formalizm (12x)
  - core/sek03_rezimy (6x)
  - core/sek08_formalizm (3x)
  - core/sek08b_ghost_resolution (3x)
  - core/sek09_cechowanie (1x)
  - papers_external/tgp_english_summary (1x)
  - partial_proofs/bh_ringdown (22x)
  - partial_proofs/particle_sector (1x)
  - partial_proofs/zero_mode (4x)
- Wikilinks to:
  - --
- Cite count: 2

### partial_proofs/koide_fp
- `\input`:
  - --
- `\ref` to other folders:
  - core/sek07_predykcje (1x)
  - partial_proofs/bh_ringdown (1x)
  - partial_proofs/hierarchia_mas (8x)
  - partial_proofs/particle_sector (2x)
  - partial_proofs/potencjal (1x)
  - partial_proofs/zero_mode (2x)
- Wikilinks to:
  - --
- Cite count: 2

### partial_proofs/most_gamma_phi
- `\input`:
  - --
- `\ref` to other folders:
  - core/formalizm (3x)
  - core/sek08_formalizm (1x)
  - core/sek08a_akcja_zunifikowana (1x)
  - core/sek10_N0_wyprowadzenie (1x)
  - partial_proofs/hierarchia_mas (5x)
  - partial_proofs/koide_fp (4x)
- Wikilinks to:
  - --
- Cite count: 0

### partial_proofs/nuclear_from_soliton
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - core/formalizm (2x)
  - core/sek10_N0_wyprowadzenie (2x)
  - partial_proofs/fermion_from_soliton (2x)
  - partial_proofs/hierarchia_mas (2x)
  - research/atom_from_soliton (2x)
  - research/cohesion_closure (1x)
  - research/nbody (8x)
  - research/nbody/paper (3x)
- Cite count: 0

### partial_proofs/particle_sector
- `\input`:
  - --
- `\ref` to other folders:
  - core/formalizm (7x)
  - core/sek01_ontologia (2x)
  - core/sek03_rezimy (1x)
  - partial_proofs/hierarchia_mas (2x)
  - partial_proofs/koide_fp (3x)
  - partial_proofs/quark_sector (6x)
- Wikilinks to:
  - --
- Cite count: 0

### partial_proofs/potencjal
- `\input`:
  - --
- `\ref` to other folders:
  - partial_proofs/hierarchia_mas (1x)
- Wikilinks to:
  - --
- Cite count: 0

### partial_proofs/quark_sector
- `\input`:
  - --
- `\ref` to other folders:
  - core/sek03_rezimy (2x)
  - core/sek08_formalizm (3x)
  - partial_proofs/hierarchia_mas (5x)
  - partial_proofs/zero_mode (1x)
- Wikilinks to:
  - --
- Cite count: 0

### partial_proofs/superconductivity
- `\input`:
  - --
- `\ref` to other folders:
  - axioms/substrat (1x)
  - core/formalizm (3x)
  - core/sek02_pole (3x)
  - partial_proofs/bounce_topology (2x)
  - partial_proofs/hierarchia_mas (2x)
  - partial_proofs/particle_sector (1x)
  - partial_proofs/zero_mode (1x)
- Wikilinks to:
  - --
- Cite count: 0

### partial_proofs/trojcialowe_nbody
- `\input`:
  - research/nbody/paper (1x)
- `\ref` to other folders:
  - core/sek03_rezimy (1x)
  - core/sek04_stale (2x)
  - core/sek08_formalizm (3x)
- Wikilinks to:
  - --
- Cite count: 2

### partial_proofs/wielki_wybuch
- `\input`:
  - --
- `\ref` to other folders:
  - axioms/substrat (8x)
  - core/formalizm (2x)
  - core/sek01_ontologia (2x)
  - core/sek03_rezimy (2x)
  - core/sek04_stale (2x)
  - core/sek05_ciemna_energia (1x)
  - core/sek06_czarne_dziury (3x)
  - core/sek07_predykcje (2x)
  - core/sek08_formalizm (9x)
- Wikilinks to:
  - --
- Cite count: 0

### partial_proofs/zero_mode
- `\input`:
  - --
- `\ref` to other folders:
  - core/sek07_predykcje (1x)
  - papers_external/tgp_english_summary (1x)
  - partial_proofs/bh_ringdown (1x)
  - partial_proofs/hierarchia_mas (8x)
  - partial_proofs/koide_fp (8x)
- Wikilinks to:
  - --
- Cite count: 0

### research
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - --
- Cite count: 0

### research/_archive
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (1x)
  - audyt/D01_drifting_numbers (3x)
  - meta (3x)
  - meta/research (4x)
- Cite count: 0

### research/_sandbox
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (1x)
  - audyt/D01_drifting_numbers (1x)
  - meta (2x)
  - meta/research (2x)
  - meta/research/templates (1x)
- Cite count: 0

### research/atom_from_soliton
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (1x)
  - audyt/D01_drifting_numbers (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - core/formalizm (1x)
  - core/sek09_cechowanie (1x)
  - meta/research (5x)
  - meta/research/templates (1x)
  - partial_proofs/fermion_from_soliton (1x)
  - research/atomic_shells_closure (2x)
  - research/cohesion_closure (1x)
  - research/em_from_substrate (5x)
- Cite count: 0

### research/atomic_shells_closure
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (1x)
  - audyt/D01_drifting_numbers (8x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (5x)
  - meta/research/templates (1x)
  - papers_external/tgp_core_paper (1x)
  - partial_proofs/fermion_from_soliton (1x)
  - research/superconductivity_closure (2x)
- Cite count: 0

### research/audyt_cosmology_drift_2026-05-03
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (1x)
  - audyt/D01_drifting_numbers (4x)
  - research/closure_2026-04-26 (1x)
  - research/closure_2026-04-26/alpha_psi_threshold (2x)
  - research/op-cosmology-closure (3x)
  - research/op-omicron2-phi-mean-shift-cosmo (1x)
- Cite count: 0

### research/brannen_sqrt2
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - audyt/D01_drifting_numbers (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (4x)
- Cite count: 0

### research/cabibbo_correction
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - audyt/D01_drifting_numbers (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (4x)
- Cite count: 0

### research/casimir_mof
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (1x)
  - audyt/D01_drifting_numbers (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (5x)
  - meta/research/templates (1x)
  - papers_external/tgp_core_paper (1x)
  - partial_proofs/fermion_from_soliton (1x)
  - research/superconductivity_closure (1x)
- Cite count: 0

### research/closure_2026-04-26
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (2x)
  - audyt/D01_drifting_numbers (8x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta (5x)
  - meta/research (5x)
  - meta/research/templates (1x)
  - research/closure_2026-04-26/alpha_psi_threshold (10x)
  - research/closure_2026-04-26/sigma_ab_pathB (1x)
  - research/op-cosmology-closure (7x)
  - research/op-m92 (1x)
  - research/op-newton-momentum (2x)
  - research/op-phase1-covariant (11x)
  - research/op-phase2-quantum-gravity (12x)
  - research/op-phase3-uv-completion (19x)
  - research/op-quantum-closure (5x)
  - research/op7 (7x)
- Cite count: 0

### research/closure_2026-04-26/Lambda_from_Phi0
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (2x)
  - audyt/D01_drifting_numbers (2x)
  - research/closure_2026-04-26/alpha_psi_threshold (1x)
  - research/op-newton-momentum (3x)
  - research/op7 (3x)
- Cite count: 0

### research/closure_2026-04-26/alpha_psi_threshold
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (1x)
  - research/op-m92 (8x)
  - research/op7 (2x)
- Cite count: 0

### research/closure_2026-04-26/f_psi_principle
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (2x)
  - papers_external/tgp_core_paper (1x)
  - research/op-newton-momentum (6x)
  - research/op7 (3x)
- Cite count: 0

### research/closure_2026-04-26/sigma_ab_pathB
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (3x)
  - papers_external/tgp_core_paper (1x)
  - research/closure_2026-04-26 (1x)
  - research/closure_2026-04-26/alpha_psi_threshold (1x)
  - research/op7 (8x)
- Cite count: 0

### research/cohesion_closure
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (1x)
  - audyt/D01_drifting_numbers (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (5x)
  - meta/research/templates (1x)
  - partial_proofs/fermion_from_soliton (1x)
  - research/atomic_shells_closure (2x)
  - research/em_from_substrate (2x)
- Cite count: 0

### research/continuum_limit
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - audyt/D01_drifting_numbers (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (4x)
- Cite count: 0

### research/cosmo_tensions
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - audyt/D01_drifting_numbers (6x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (4x)
  - research/op-cosmology-closure (1x)
- Cite count: 0

### research/desi_dark_energy
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - audyt/D01_drifting_numbers (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (4x)
  - research (1x)
  - research/closure_2026-04-26 (1x)
  - research/galaxy_scaling (3x)
  - research/op-cosmology-closure (3x)
- Cite count: 0

### research/em_from_substrate
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (1x)
  - audyt/D01_drifting_numbers (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - core/formalizm (1x)
  - meta/research (5x)
  - meta/research/templates (1x)
  - partial_proofs/fermion_from_soliton (1x)
  - research/nbody/examples (1x)
  - tooling/scripts (1x)
  - tooling/scripts/gauge (2x)
- Cite count: 0

### research/external_review_2026-04-25
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (1x)
  - audyt/D01_drifting_numbers (6x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (5x)
  - meta/research/templates (1x)
  - papers_external/tgp_core_paper (2x)
  - research/closure_2026-04-26 (1x)
  - research/op1-op2-op4 (2x)
  - research/op7 (1x)
- Cite count: 0

### research/galaxy_scaling
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - audyt/D01_drifting_numbers (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (4x)
  - research/cosmo_tensions (2x)
- Cite count: 0

### research/hubble_tension
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - audyt/D01_drifting_numbers (6x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (4x)
  - research/op-cosmology-closure (2x)
- Cite count: 0

### research/liquid_viscosity
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (1x)
  - audyt/D01_drifting_numbers (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (5x)
  - meta/research/templates (1x)
  - papers_external/tgp_core_paper (1x)
  - partial_proofs/fermion_from_soliton (1x)
  - research/superconductivity_closure (1x)
- Cite count: 0

### research/mass_scaling_k4
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - audyt/D01_drifting_numbers (6x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (4x)
  - research/why_n3 (4x)
- Cite count: 0

### research/metric_ansatz
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - audyt/D01_drifting_numbers (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (4x)
- Cite count: 0

### research/muon_g_minus_2
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (1x)
  - audyt/D01_drifting_numbers (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (5x)
  - meta/research/templates (1x)
  - papers_external/tgp_core_paper (1x)
  - partial_proofs/fermion_from_soliton (1x)
  - research/casimir_mof (3x)
  - research/superconductivity_closure (1x)
- Cite count: 0

### research/nbody
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - audyt/D01_drifting_numbers (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (4x)
- Cite count: 0

### research/nbody/docs
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - audyt/D01_drifting_numbers (1x)
- Cite count: 0

### research/nbody/examples
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - --
- Cite count: 0

### research/nbody/examples/_outputs
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - --
- Cite count: 0

### research/nbody/paper
- `\input`:
  - --
- `\ref` to other folders:
  - core/formalizm (3x)
  - core/sek06_czarne_dziury (2x)
  - core/sek08_formalizm (12x)
  - core/sek08a_akcja_zunifikowana (1x)
  - papers_external/arxiv_submission (24x)
  - partial_proofs/trojcialowe_nbody (1x)
- Wikilinks to:
  - audyt/D01_drifting_numbers (1x)
- Cite count: 3

### research/neutrino_msw
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (1x)
  - audyt/D01_drifting_numbers (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (5x)
  - meta/research/templates (1x)
  - papers_external/tgp_core_paper (1x)
- Cite count: 0

### research/op-CORE-CLEANUP-B-2026-05-04
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (2x)
  - audyt (4x)
  - audyt/D01_drifting_numbers (10x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - core/formalizm (2x)
  - core/sek08a_akcja_zunifikowana (3x)
  - core/sek08c_metryka_z_substratu (4x)
  - partial_proofs/fermion_from_soliton (1x)
  - research/audyt_cosmology_drift_2026-05-03 (5x)
- Cite count: 0

### research/op-D01-anchor-lock-2026-05-06
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - audyt (2x)
  - audyt/D01_drifting_numbers (17x)
  - audyt/T01_LIGO3G_falsifier (2x)
  - meta (7x)
  - research/op-newton-momentum (7x)
- Cite count: 0

### research/op-FRW-radiation-era-varying-c-2026-05-06
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (2x)
  - audyt (9x)
  - audyt/D01_drifting_numbers (24x)
  - audyt/L01_rho_operational (9x)
  - audyt/T01_LIGO3G_falsifier (2x)
  - core/sek04_stale (3x)
  - core/sek05_ciemna_energia (3x)
  - core/sek08a_akcja_zunifikowana (3x)
  - meta (7x)
  - research/op-D01-anchor-lock-2026-05-06 (9x)
  - research/op-M03-balance-sheet-retrofit-2026-05-06 (4x)
  - research/op-alpha-fine-structure (5x)
  - research/op-cosmology-closure (4x)
  - research/op-newton-momentum (2x)
- Cite count: 0

### research/op-GWTC3-reanalysis
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (4x)
  - audyt/D01_drifting_numbers (12x)
  - audyt/T01_LIGO3G_falsifier (13x)
  - meta (1x)
  - papers/M911_LIGO3G_paper (3x)
  - research/op-D01-anchor-lock-2026-05-06 (4x)
  - research/op-GWTC3-reanalysis/scripts (11x)
  - research/op-LIGO-3G-deviation (10x)
  - research/op-alpha-fine-structure (6x)
  - research/op-ppE-mapping (6x)
  - research/op-ppE-mapping/scripts (1x)
- Cite count: 0

### research/op-GWTC3-reanalysis/scripts
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - --
- Cite count: 0

### research/op-L01-rho-stress-energy-bridge-2026-05-04
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (1x)
  - audyt/D01_drifting_numbers (15x)
  - audyt/L01_rho_operational (2x)
  - audyt/T01_LIGO3G_falsifier (2x)
  - core/sek08_formalizm (8x)
  - core/sek08a_akcja_zunifikowana (4x)
  - meta (4x)
  - research/closure_2026-04-26/Lambda_from_Phi0 (4x)
  - research/op-newton-momentum (6x)
  - research/op-psi1-substrate-light-acceleration (1x)
  - research/op-tau3-substrate-clock-acceleration (1x)
- Cite count: 0

### research/op-L03-spectral-stability-2026-05-06
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (2x)
  - audyt (1x)
  - audyt/D01_drifting_numbers (24x)
  - audyt/T01_LIGO3G_falsifier (3x)
  - axioms/notacja (1x)
  - core/sek01_ontologia (1x)
  - core/sek08_formalizm (7x)
  - core/sek08a_akcja_zunifikowana (7x)
  - core/sek08b_ghost_resolution (7x)
  - meta (1x)
- Cite count: 0

### research/op-L04-ODE-canonicalization-2026-05-04
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - audyt/D01_drifting_numbers (18x)
  - audyt/L01_rho_operational (2x)
  - audyt/T01_LIGO3G_falsifier (2x)
  - core/sek08_formalizm (9x)
  - core/sek08a_akcja_zunifikowana (2x)
  - core/sek08b_ghost_resolution (1x)
  - meta (8x)
  - research/mass_scaling_k4 (12x)
  - research/op-alpha-fine-structure (1x)
  - research/why_n3 (49x)
- Cite count: 0

### research/op-LIGO-3G-deviation
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (5x)
  - audyt/D01_drifting_numbers (12x)
  - audyt/T01_LIGO3G_falsifier (14x)
  - meta (1x)
  - research/op-D01-anchor-lock-2026-05-06 (3x)
  - research/op-GWTC3-reanalysis (1x)
  - research/op-LIGO-3G-deviation/scripts (5x)
  - research/op-alpha-fine-structure (13x)
  - research/op-ppE-mapping (4x)
- Cite count: 0

### research/op-LIGO-3G-deviation/scripts
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - --
- Cite count: 0

### research/op-M03-balance-sheet-retrofit-2026-05-06
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (5x)
  - audyt (2x)
  - audyt/D01_drifting_numbers (186x)
  - audyt/M03_balance_sheet_missing (3x)
  - audyt/T01_LIGO3G_falsifier (3x)
  - core/_meta_latex (1x)
  - core/formalizm (1x)
  - core/sek00_summary (1x)
  - meta (100x)
  - partial_proofs/fermion_from_soliton (3x)
  - research/closure_2026-04-26 (1x)
  - research/closure_2026-04-26/alpha_psi_threshold (10x)
  - research/op-alpha-fine-structure (107x)
  - research/op-chi1-newton-constant-derivation (1x)
  - research/op-cosmology-closure (9x)
  - research/op-eht (10x)
  - research/op-eht-A (7x)
  - research/op-lambda1-e2-amplitude-emergence (2x)
  - research/op-newton-momentum (13x)
  - research/op-omicron2-phi-mean-shift-cosmo (5x)
  - research/op-phase1-covariant (4x)
  - research/op-phase2-quantum-gravity (4x)
  - research/op-phase3-uv-completion (5x)
  - research/op-psi1-substrate-light-acceleration (7x)
  - research/op-quantum-closure (8x)
  - research/op-uv2-mtgp-absolute-scale (4x)
  - research/op7 (1x)
- Cite count: 0

### research/op-MAG-Lorentz-A-mu-coupling-2026-05-09
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - audyt/D01_drifting_numbers (8x)
  - meta (1x)
  - research/op-D01-anchor-lock-2026-05-06 (2x)
  - research/op-MAG-resonance-formalization-2026-05-09 (2x)
  - research/op-alpha-fine-structure (1x)
- Cite count: 0

### research/op-MAG-Phase5-V-reference-clarification-2026-05-09
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - audyt/D01_drifting_numbers (8x)
  - core/sek08a_akcja_zunifikowana (2x)
  - research/op-D01-anchor-lock-2026-05-06 (3x)
  - research/op-MAG-resonance-formalization-2026-05-09 (10x)
  - research/op-V-canonical-consistency-audit-2026-05-09 (4x)
- Cite count: 0

### research/op-MAG-anomalous-moment-2026-05-09
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - audyt/D01_drifting_numbers (6x)
  - research/op-D01-anchor-lock-2026-05-06 (1x)
  - research/op-MAG-Lorentz-A-mu-coupling-2026-05-09 (2x)
  - research/op-SPIN-SU2-substrate-derivation-2026-05-08 (4x)
- Cite count: 0

### research/op-MAG-resonance-formalization-2026-05-09
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - audyt/D01_drifting_numbers (29x)
  - core/sek08a_akcja_zunifikowana (1x)
  - meta (2x)
  - research/op-D01-anchor-lock-2026-05-06 (3x)
  - research/op-MAG-Lorentz-A-mu-coupling-2026-05-09 (4x)
  - research/op-Phase5-MAG-erratum-2026-05-09 (2x)
  - research/op-SPIN-SU2-substrate-derivation-2026-05-08 (11x)
  - research/op-alpha-fine-structure (6x)
- Cite count: 0

### research/op-Phase5-MAG-erratum-2026-05-09
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - audyt/D01_drifting_numbers (3x)
  - research/op-MAG-resonance-formalization-2026-05-09 (7x)
  - research/op-alpha-fine-structure (4x)
- Cite count: 0

### research/op-Phi-decomposition-photon-2026-05-07
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - audyt (2x)
  - audyt/D01_drifting_numbers (29x)
  - audyt/L01_rho_operational (4x)
  - audyt/T01_LIGO3G_falsifier (7x)
  - core/sek04_stale (4x)
  - core/sek08a_akcja_zunifikowana (4x)
  - meta (8x)
  - research/op-D01-anchor-lock-2026-05-06 (7x)
  - research/op-M03-balance-sheet-retrofit-2026-05-06 (2x)
  - research/op-alpha-fine-structure (7x)
  - research/op-emergent-metric-from-interaction-2026-05-09 (6x)
- Cite count: 0

### research/op-Phi-vacuum-scale-2026-05-09
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (2x)
  - audyt/D01_drifting_numbers (26x)
  - core/sek08a_akcja_zunifikowana (4x)
  - core/sek08c_metryka_z_substratu (1x)
  - meta (4x)
  - research/closure_2026-04-26/alpha_psi_threshold (2x)
  - research/op-D01-anchor-lock-2026-05-06 (6x)
  - research/op-GWTC3-reanalysis (1x)
  - research/op-MAG-Lorentz-A-mu-coupling-2026-05-09 (3x)
  - research/op-MAG-resonance-formalization-2026-05-09 (12x)
  - research/op-alpha-fine-structure (9x)
  - research/op-ppE-mapping (1x)
- Cite count: 0

### research/op-Phi0-spatial-variation-predictions-2026-05-09
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - audyt/D01_drifting_numbers (3x)
  - meta (1x)
  - research/op-Phi-vacuum-scale-2026-05-09 (4x)
  - research/op-alpha-fine-structure (1x)
- Cite count: 0

### research/op-S07-alternative-f-psi-derivation-2026-05-09
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (4x)
  - audyt (1x)
  - audyt/D01_drifting_numbers (11x)
  - core/sek08a_akcja_zunifikowana (1x)
  - core/sek08c_metryka_z_substratu (1x)
  - meta (1x)
  - research/op-D01-anchor-lock-2026-05-06 (2x)
  - research/op-GWTC3-reanalysis (5x)
  - research/op-Phi-vacuum-scale-2026-05-09 (2x)
  - research/op-alpha-fine-structure (2x)
  - research/op-g0-r3-from-canonical-projection (3x)
  - research/op-ppE-mapping (5x)
- Cite count: 0

### research/op-SPIN-MAG-leakage-lean-hypothesis-2026-05-07
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - meta (1x)
  - research/op-alpha-fine-structure (4x)
- Cite count: 0

### research/op-SPIN-SU2-substrate-derivation-2026-05-08
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (1x)
  - audyt/D01_drifting_numbers (44x)
  - core/sek08a_akcja_zunifikowana (3x)
  - meta (2x)
  - research/op-D01-anchor-lock-2026-05-06 (2x)
  - research/op-MAG-Lorentz-A-mu-coupling-2026-05-09 (1x)
  - research/op-MAG-resonance-formalization-2026-05-09 (1x)
  - research/op-alpha-fine-structure (5x)
- Cite count: 0

### research/op-V-canonical-consistency-audit-2026-05-09
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - audyt/D01_drifting_numbers (9x)
  - core/sek08a_akcja_zunifikowana (3x)
  - meta (2x)
  - research/closure_2026-04-26/alpha_psi_threshold (5x)
  - research/op-D01-anchor-lock-2026-05-06 (3x)
  - research/op-MAG-Phase5-V-reference-clarification-2026-05-09 (1x)
  - research/op-MAG-resonance-formalization-2026-05-09 (5x)
  - research/op-Phi-vacuum-scale-2026-05-09 (7x)
  - research/op-alpha-fine-structure (1x)
  - research/op-g0-r3-from-canonical-projection (4x)
- Cite count: 0

### research/op-alpha-fine-structure
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (1x)
  - audyt/D01_drifting_numbers (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (5x)
  - meta/research/templates (1x)
- Cite count: 0

### research/op-bh-alpha-threshold
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (3x)
  - audyt/D01_drifting_numbers (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (5x)
  - meta/research/templates (1x)
  - research/closure_2026-04-26/alpha_psi_threshold (1x)
  - research/op-alpha-fine-structure (10x)
  - research/op-m92 (3x)
- Cite count: 0

### research/op-chi1-newton-constant-derivation
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (10x)
  - audyt/D01_drifting_numbers (3x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - core/sek08_formalizm (1x)
  - meta (4x)
  - meta/research (1x)
  - meta/research/templates (1x)
  - research/op-alpha-fine-structure (52x)
  - research/op-phase2-quantum-gravity (6x)
- Cite count: 0

### research/op-cosmology-closure
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (3x)
  - audyt/D01_drifting_numbers (8x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - core/sek08a_akcja_zunifikowana (4x)
  - meta (2x)
  - meta/research (5x)
  - meta/research/templates (1x)
  - research/closure_2026-04-26 (10x)
  - research/closure_2026-04-26/alpha_psi_threshold (14x)
  - research/cosmo_tensions (16x)
  - research/desi_dark_energy (7x)
  - research/galaxy_scaling (17x)
  - research/nbody/examples (8x)
  - research/op-newton-momentum (11x)
- Cite count: 0

### research/op-cross-sector-charge
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (2x)
  - audyt/D01_drifting_numbers (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (5x)
  - meta/research/templates (1x)
  - research/op-alpha-fine-structure (14x)
- Cite count: 0

### research/op-delta1-g-tilde-derivation
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - audyt/D01_drifting_numbers (9x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (4x)
  - research/closure_2026-04-26/alpha_psi_threshold (2x)
  - research/op-lambda1-e2-amplitude-emergence (2x)
- Cite count: 0

### research/op-delta2-Nf-derivation
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - audyt/D01_drifting_numbers (10x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (4x)
  - partial_proofs/hierarchia_mas (1x)
  - research/closure_2026-04-26/alpha_psi_threshold (1x)
- Cite count: 0

### research/op-dual-V-structure-clarification-2026-05-09
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - audyt/D01_drifting_numbers (8x)
  - core/sek08a_akcja_zunifikowana (3x)
  - meta (2x)
  - research/closure_2026-04-26/alpha_psi_threshold (2x)
  - research/op-D01-anchor-lock-2026-05-06 (3x)
  - research/op-MAG-Phase5-V-reference-clarification-2026-05-09 (5x)
  - research/op-MAG-resonance-formalization-2026-05-09 (2x)
  - research/op-alpha-fine-structure (8x)
  - research/op-g0-r3-from-canonical-projection (5x)
- Cite count: 0

### research/op-eht
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (1x)
  - audyt/D01_drifting_numbers (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (5x)
  - meta/research/templates (1x)
  - papers_external/tgp_core_paper (5x)
  - research (1x)
  - research/closure_2026-04-26 (3x)
  - research/op-newton-momentum (4x)
  - research/op7 (10x)
  - tooling/scripts/gravity (1x)
- Cite count: 0

### research/op-eht-A
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (1x)
  - audyt/D01_drifting_numbers (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (5x)
  - meta/research/templates (1x)
  - papers_external/tgp_core_paper (2x)
  - research/closure_2026-04-26 (1x)
  - research/op-eht (4x)
  - research/op-newton-momentum (2x)
  - research/op7 (1x)
- Cite count: 0

### research/op-emergent-metric-from-interaction-2026-05-09
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (3x)
  - audyt/D01_drifting_numbers (11x)
  - research/op-D01-anchor-lock-2026-05-06 (2x)
  - research/op-GWTC3-reanalysis (3x)
  - research/op-MAG-Lorentz-A-mu-coupling-2026-05-09 (1x)
  - research/op-SPIN-SU2-substrate-derivation-2026-05-08 (3x)
  - research/op-alpha-fine-structure (2x)
  - research/op-ppE-mapping (4x)
- Cite count: 0

### research/op-eps-photon-ring
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (1x)
  - audyt/D01_drifting_numbers (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (5x)
  - meta/research/templates (1x)
  - research/op-alpha-fine-structure (12x)
- Cite count: 0

### research/op-eta-wolfenstein
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (1x)
  - audyt/D01_drifting_numbers (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (5x)
  - meta/research/templates (1x)
  - research/op-alpha-fine-structure (13x)
- Cite count: 0

### research/op-eta2-denom-derivation
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (2x)
  - audyt/D01_drifting_numbers (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (5x)
  - meta/research/templates (1x)
  - research/op-alpha-fine-structure (16x)
- Cite count: 0

### research/op-g0-r3-from-canonical-projection
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (1x)
  - audyt/D01_drifting_numbers (16x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - core/sek08a_akcja_zunifikowana (1x)
  - core/sek08c_metryka_z_substratu (1x)
  - meta (1x)
  - meta/research (4x)
  - research/mass_scaling_k4 (1x)
  - research/op-alpha-fine-structure (14x)
  - research/why_n3 (2x)
- Cite count: 0

### research/op-gamma1-phi-eff-anchor-resolution
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - audyt/D01_drifting_numbers (8x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (4x)
  - partial_proofs/fermion_from_soliton (1x)
  - research/closure_2026-04-26/alpha_psi_threshold (2x)
  - research/op-lambda1-e2-amplitude-emergence (4x)
  - research/why_n3 (1x)
- Cite count: 0

### research/op-iota-charge-pmns-unification
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (7x)
  - audyt/D01_drifting_numbers (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (5x)
  - meta/research/templates (1x)
  - research/op-alpha-fine-structure (40x)
- Cite count: 0

### research/op-kappa-mixing-numerator
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (6x)
  - audyt/D01_drifting_numbers (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (5x)
  - meta/research/templates (1x)
  - research/op-alpha-fine-structure (35x)
- Cite count: 0

### research/op-lambda1-e2-amplitude-emergence
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - audyt/D01_drifting_numbers (13x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (4x)
  - research/mass_scaling_k4 (8x)
  - research/op-alpha-fine-structure (10x)
  - research/op-delta1-g-tilde-derivation (2x)
  - research/op-delta2-Nf-derivation (1x)
  - research/op-gamma1-phi-eff-anchor-resolution (2x)
  - research/why_n3 (3x)
- Cite count: 0

### research/op-m92
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (1x)
  - audyt/D01_drifting_numbers (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (5x)
  - meta/research/templates (1x)
  - papers_external/tgp_core_paper (3x)
  - research/closure_2026-04-26 (1x)
  - research/op-eht (1x)
  - research/op-eht-A (3x)
  - research/op-newton-momentum (1x)
  - research/op7 (2x)
- Cite count: 0

### research/op-mu-pmns-phase-hardening
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (11x)
  - audyt/D01_drifting_numbers (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (5x)
  - meta/research/templates (1x)
  - research/op-alpha-fine-structure (43x)
- Cite count: 0

### research/op-mu1-minimal-substrate-log-redefinition
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - audyt/D01_drifting_numbers (9x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (4x)
  - partial_proofs/fermion_from_soliton (3x)
  - research/op-lambda1-e2-amplitude-emergence (5x)
- Cite count: 0

### research/op-newton-momentum
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (9x)
  - audyt/D01_drifting_numbers (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta (13x)
  - meta/research (5x)
  - meta/research/templates (1x)
  - research/closure_2026-04-26 (6x)
  - research/closure_2026-04-26/alpha_psi_threshold (11x)
  - research/op-m92 (2x)
- Cite count: 0

### research/op-nu-majorana-phase-mbb
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (7x)
  - audyt/D01_drifting_numbers (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (5x)
  - meta/research/templates (1x)
  - research/op-alpha-fine-structure (40x)
- Cite count: 0

### research/op-omega1-substrate-em-coupling
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (5x)
  - audyt/D01_drifting_numbers (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (5x)
  - meta/research/templates (1x)
  - partial_proofs/fermion_from_soliton (2x)
  - research/op-alpha-fine-structure (33x)
- Cite count: 0

### research/op-omega2-axion-coupling-lock
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (7x)
  - audyt/D01_drifting_numbers (3x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta (7x)
  - meta/research (1x)
  - meta/research/templates (1x)
  - research/op-alpha-fine-structure (51x)
  - research/op-chi1-newton-constant-derivation (1x)
  - research/op-omega3-axion-decay-constant (1x)
  - research/op-uv2-mtgp-absolute-scale (1x)
- Cite count: 0

### research/op-omega3-axion-decay-constant
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (4x)
  - audyt/D01_drifting_numbers (3x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta (8x)
  - meta/research (1x)
  - meta/research/templates (1x)
  - research/op-alpha-fine-structure (45x)
  - research/op-chi1-newton-constant-derivation (1x)
  - research/op-omega2-axion-coupling-lock (1x)
  - research/op-uv2-mtgp-absolute-scale (1x)
- Cite count: 0

### research/op-omicron1-sigmamnu-cosmo
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (6x)
  - audyt/D01_drifting_numbers (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (5x)
  - meta/research/templates (1x)
  - research/op-alpha-fine-structure (38x)
- Cite count: 0

### research/op-omicron2-phi-mean-shift-cosmo
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (1x)
  - audyt/D01_drifting_numbers (9x)
  - core/sek01_ontologia (1x)
  - core/sek08a_akcja_zunifikowana (1x)
  - research/closure_2026-04-26/alpha_psi_threshold (4x)
  - research/op-alpha-fine-structure (1x)
  - research/op-cosmology-closure (7x)
- Cite count: 0

### research/op-phase1-covariant
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (1x)
  - audyt/D01_drifting_numbers (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (5x)
  - meta/research/templates (1x)
  - research/closure_2026-04-26 (5x)
  - research/closure_2026-04-26/alpha_psi_threshold (9x)
  - research/op-cosmology-closure (2x)
  - research/op-newton-momentum (4x)
  - research/op-quantum-closure (14x)
  - research/op1-op2-op4 (2x)
- Cite count: 0

### research/op-phase2-quantum-gravity
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (1x)
  - audyt/D01_drifting_numbers (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (5x)
  - meta/research/templates (1x)
  - research/closure_2026-04-26 (10x)
  - research/closure_2026-04-26/alpha_psi_threshold (6x)
  - research/op-cosmology-closure (2x)
  - research/op-newton-momentum (3x)
  - research/op-phase1-covariant (16x)
  - research/op-quantum-closure (3x)
- Cite count: 0

### research/op-phase3-uv-completion
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (1x)
  - audyt/D01_drifting_numbers (7x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (5x)
  - meta/research/templates (1x)
  - research/closure_2026-04-26 (3x)
  - research/closure_2026-04-26/alpha_psi_threshold (6x)
  - research/op-cosmology-closure (1x)
  - research/op-newton-momentum (3x)
  - research/op-phase1-covariant (4x)
  - research/op-phase2-quantum-gravity (10x)
  - research/op-quantum-closure (1x)
- Cite count: 0

### research/op-phi1-substrate-action-variational
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (5x)
  - audyt/D01_drifting_numbers (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (5x)
  - meta/research/templates (1x)
  - research/op-alpha-fine-structure (34x)
- Cite count: 0

### research/op-pi1-bb0nu-nme-isotope
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (6x)
  - audyt/D01_drifting_numbers (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (5x)
  - meta/research/templates (1x)
  - research/op-alpha-fine-structure (39x)
- Cite count: 0

### research/op-ppE-mapping
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (7x)
  - audyt/D01_drifting_numbers (17x)
  - audyt/T01_LIGO3G_falsifier (34x)
  - meta (1x)
  - papers/M911_LIGO3G_paper (3x)
  - research/op-D01-anchor-lock-2026-05-06 (3x)
  - research/op-GWTC3-reanalysis (5x)
  - research/op-LIGO-3G-deviation (3x)
  - research/op-alpha-fine-structure (14x)
  - research/op-newton-momentum (9x)
  - research/op-ppE-mapping/scripts (7x)
- Cite count: 0

### research/op-ppE-mapping/TGP/TGP_v1/audyt/T01_LIGO3G_falsifier
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (3x)
  - audyt/D01_drifting_numbers (7x)
  - audyt/T01_LIGO3G_falsifier (11x)
  - papers/M911_LIGO3G_paper (2x)
- Cite count: 0

### research/op-ppE-mapping/TGP/TGP_v1/papers/M911_LIGO3G_paper
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (1x)
  - research/op-LIGO-3G-deviation (1x)
  - research/op-newton-momentum (1x)
  - research/op-ppE-mapping (1x)
- Cite count: 0

### research/op-ppE-mapping/TGP/TGP_v1/research/op-LIGO-3G-deviation
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (5x)
  - audyt/D01_drifting_numbers (12x)
  - audyt/T01_LIGO3G_falsifier (14x)
  - meta (1x)
  - research/op-D01-anchor-lock-2026-05-06 (3x)
  - research/op-LIGO-3G-deviation (7x)
  - research/op-LIGO-3G-deviation/scripts (5x)
  - research/op-alpha-fine-structure (13x)
  - research/op-ppE-mapping (3x)
- Cite count: 0

### research/op-ppE-mapping/TGP/TGP_v1/research/op-ppE-mapping
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (3x)
  - audyt/D01_drifting_numbers (10x)
  - audyt/T01_LIGO3G_falsifier (13x)
  - research/op-D01-anchor-lock-2026-05-06 (1x)
  - research/op-alpha-fine-structure (3x)
  - research/op-newton-momentum (1x)
  - research/op-ppE-mapping (4x)
  - research/op-ppE-mapping/scripts (2x)
- Cite count: 0

### research/op-ppE-mapping/scripts
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - --
- Cite count: 0

### research/op-psi1-substrate-light-acceleration
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (2x)
  - audyt/D01_drifting_numbers (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta (3x)
  - meta/research (5x)
  - meta/research/templates (1x)
  - research/op-alpha-fine-structure (13x)
- Cite count: 0

### research/op-quantum-closure
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (2x)
  - audyt/D01_drifting_numbers (6x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - core/sek08a_akcja_zunifikowana (2x)
  - meta/research (5x)
  - meta/research/templates (1x)
  - research/closure_2026-04-26 (7x)
  - research/closure_2026-04-26/Lambda_from_Phi0 (2x)
  - research/closure_2026-04-26/alpha_psi_threshold (2x)
  - research/continuum_limit (10x)
  - research/muon_g_minus_2 (1x)
  - research/op-cosmology-closure (9x)
  - research/op-newton-momentum (7x)
  - research/op1-op2-op4 (21x)
  - tooling/scripts (1x)
- Cite count: 0

### research/op-rho1-71Ge-cross-section
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (10x)
  - audyt/D01_drifting_numbers (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (5x)
  - meta/research/templates (1x)
  - research/op-alpha-fine-structure (42x)
- Cite count: 0

### research/op-sc-alpha-origin
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (8x)
  - audyt/D01_drifting_numbers (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (5x)
  - meta/research/templates (1x)
  - research/closure_2026-04-26/alpha_psi_threshold (2x)
  - research/op-alpha-fine-structure (13x)
- Cite count: 0

### research/op-sigma1-substrate-light-dispersion
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (5x)
  - audyt/D01_drifting_numbers (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (5x)
  - meta/research/templates (1x)
  - research/op-alpha-fine-structure (26x)
- Cite count: 0

### research/op-tau1-closure-overlap-coulomb
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (6x)
  - audyt/D01_drifting_numbers (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (5x)
  - meta/research/templates (1x)
  - research/op-alpha-fine-structure (35x)
- Cite count: 0

### research/op-tau2-substrate-time-coupling
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (3x)
  - audyt/D01_drifting_numbers (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (5x)
  - meta/research/templates (1x)
  - research/op-alpha-fine-structure (10x)
- Cite count: 0

### research/op-tau3-substrate-clock-acceleration
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (4x)
  - audyt/D01_drifting_numbers (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta (6x)
  - meta/research (5x)
  - meta/research/templates (1x)
  - research/op-alpha-fine-structure (22x)
- Cite count: 0

### research/op-tensor-modes-Phi-FUTURE
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - meta (1x)
  - research/op-alpha-fine-structure (3x)
- Cite count: 0

### research/op-theta-quark-koide
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (1x)
  - audyt/D01_drifting_numbers (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (5x)
  - meta/research/templates (1x)
  - research/op-alpha-fine-structure (12x)
- Cite count: 0

### research/op-upsilon1-closure-cross-family
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (6x)
  - audyt/D01_drifting_numbers (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (5x)
  - meta/research/templates (1x)
  - research/op-alpha-fine-structure (36x)
- Cite count: 0

### research/op-uv-as-ngfp
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (2x)
  - audyt/D01_drifting_numbers (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (5x)
  - meta/research/templates (1x)
  - research/op-alpha-fine-structure (15x)
- Cite count: 0

### research/op-uv-renormalizability-research
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - audyt/D01_drifting_numbers (8x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (4x)
  - research/closure_2026-04-26 (2x)
  - research/op-phase3-uv-completion (6x)
- Cite count: 0

### research/op-uv2-mtgp-absolute-scale
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (6x)
  - audyt/D01_drifting_numbers (3x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta (4x)
  - meta/research (1x)
  - meta/research/templates (1x)
  - research/op-alpha-fine-structure (49x)
  - research/op-chi1-newton-constant-derivation (1x)
  - research/op-phase2-quantum-gravity (1x)
- Cite count: 0

### research/op-uv3-phi0-renormalization
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - audyt/D01_drifting_numbers (9x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - core/_meta_latex (1x)
  - core/formalizm (6x)
  - core/sek00_summary (3x)
  - core/sek08_formalizm (1x)
  - meta/research (4x)
  - research/op-alpha-fine-structure (12x)
  - research/op-uv2-mtgp-absolute-scale (4x)
- Cite count: 0

### research/op-void-flat-modes-h0-2026-05-06
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (2x)
  - audyt/D01_drifting_numbers (17x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - core/sek01_ontologia (1x)
  - core/sek08a_akcja_zunifikowana (1x)
  - meta/research (2x)
  - research/closure_2026-04-26 (1x)
  - research/closure_2026-04-26/alpha_psi_threshold (10x)
  - research/op-cosmology-closure (9x)
  - research/op-omicron2-phi-mean-shift-cosmo (4x)
- Cite count: 0

### research/op-xi-photon-ring
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (2x)
  - audyt/D01_drifting_numbers (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (5x)
  - meta/research/templates (1x)
  - research/op-alpha-fine-structure (14x)
- Cite count: 0

### research/op-xi2-sterile-nu-5sector
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (11x)
  - audyt/D01_drifting_numbers (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (5x)
  - meta/research/templates (1x)
  - research/op-alpha-fine-structure (39x)
- Cite count: 0

### research/op-zeta-mass-spectrum
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (1x)
  - audyt/D01_drifting_numbers (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (5x)
  - meta/research/templates (1x)
  - research/op-alpha-fine-structure (12x)
- Cite count: 0

### research/op1-op2-op4
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - audyt/D01_drifting_numbers (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (4x)
- Cite count: 0

### research/op6
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - audyt/D01_drifting_numbers (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (4x)
- Cite count: 0

### research/op7
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (2x)
  - audyt/D01_drifting_numbers (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (5x)
  - meta/research/templates (1x)
  - papers_external/tgp_core_paper (5x)
  - research/closure_2026-04-26 (1x)
  - research/op-newton-momentum (2x)
  - tooling/scripts/substrate (1x)
- Cite count: 0

### research/particle_sector_closure
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - audyt/D01_drifting_numbers (10x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (4x)
  - research/brannen_sqrt2 (1x)
- Cite count: 0

### research/qm_born_rule
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - audyt/D01_drifting_numbers (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (4x)
- Cite count: 0

### research/qm_decoherence
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - audyt/D01_drifting_numbers (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (4x)
- Cite count: 0

### research/qm_entanglement
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - audyt/D01_drifting_numbers (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (4x)
- Cite count: 0

### research/qm_foundations
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - audyt/D01_drifting_numbers (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (4x)
- Cite count: 0

### research/qm_measurement
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - audyt/D01_drifting_numbers (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (4x)
- Cite count: 0

### research/qm_spin
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - audyt/D01_drifting_numbers (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (4x)
- Cite count: 0

### research/qm_statistics
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - audyt/D01_drifting_numbers (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (4x)
- Cite count: 0

### research/qm_superposition
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - audyt/D01_drifting_numbers (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (4x)
- Cite count: 0

### research/rho_normal_state_closure
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (1x)
  - audyt/D01_drifting_numbers (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (5x)
  - meta/research/templates (1x)
- Cite count: 0

### research/s8_tension
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - audyt/D01_drifting_numbers (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (4x)
- Cite count: 0

### research/superconductivity_closure
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - audyt/D01_drifting_numbers (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (4x)
  - papers_external/tgp_sc_paper (1x)
  - research/rho_normal_state_closure (28x)
- Cite count: 0

### research/thermal_transport_molecular
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (1x)
  - audyt/D01_drifting_numbers (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (5x)
  - meta/research/templates (1x)
  - papers_external/tgp_core_paper (1x)
  - papers_external/tgp_sc_paper (1x)
  - partial_proofs/fermion_from_soliton (1x)
- Cite count: 0

### research/uv_completion
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - audyt/D01_drifting_numbers (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - meta/research (4x)
- Cite count: 0

### research/why_n3
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - <root> (1x)
  - audyt/D01_drifting_numbers (5x)
  - audyt/T01_LIGO3G_falsifier (1x)
  - core/formalizm (1x)
  - core/sek08c_metryka_z_substratu (1x)
  - meta/research (4x)
  - research/mass_scaling_k4 (1x)
- Cite count: 0

### tooling
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - --
- Cite count: 0

### tooling/scripts
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - --
- Cite count: 0

### tooling/scripts/gauge
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - --
- Cite count: 0

### tooling/scripts/gravity
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - --
- Cite count: 0

### tooling/scripts/substrate
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - --
- Cite count: 0

## Unresolved references (orphans)

### \ref (8)
- `axioms/substrat/dodatekB_substrat.tex` -> `ax:substrat`
- `axioms/substrat/dodatekB_substrat.tex` -> `app:A-aksjomaty`
- `axioms/substrat/dodatekB_substrat.tex` -> `app:B-mapa-params`
- `core/sek01_ontologia/sek01_ontologia.tex` -> `ax:substrat`
- `core/sek08_formalizm/sek08_formalizm.tex` -> `para:basin-stability`
- `core/sek08_formalizm/sek08_formalizm.tex` -> `ssec:disformal`
- `core/sek08_formalizm/sek08_formalizm.tex` -> `eq:Phi-sigma-action`
- `core/sek08_formalizm/sek08_formalizm.tex` -> `ssec:disformal-spectrum-tests`

### wikilinks (655)
- `audyt/CLEANUP_INVENTORY_2026-05-04.md` -> `[[../research/op-g0-r3-from-canonical-projection]]`
- `audyt/CLEANUP_INVENTORY_2026-05-04.md` -> `[[../research/op-L01-rho-stress-energy-bridge-2026-05-04]]`
- `audyt/CLEANUP_INVENTORY_2026-05-04.md` -> `[[../research/op-L04-ODE-canonicalization-2026-05-04]]`
- `audyt/CLEANUP_INVENTORY_2026-05-04.md` -> `[[../research/op-L01-rho-stress-energy-bridge-2026-05-04]]`
- `audyt/CLEANUP_INVENTORY_2026-05-04.md` -> `[[../research/op-L02-...]]`
- `audyt/CLEANUP_INVENTORY_2026-05-04.md` -> `[[../research/op-L04-ODE-canonicalization-2026-05-04]]`
- `audyt/CLOSURE_SUMMARY_2026-05-06.md` -> `[[audyt]]`
- `audyt/CLOSURE_SUMMARY_2026-05-06.md` -> `[[../research/op-L03-spectral-stability-2026-05-06]]`
- `audyt/CLOSURE_SUMMARY_2026-05-06.md` -> `[[../research/op-D01-anchor-lock-2026-05-06]]`
- `audyt/CLOSURE_SUMMARY_2026-05-06.md` -> `[[../research/op-M03-balance-sheet-retrofit-2026-05-06]]`
- `audyt/D01_drifting_numbers/NEEDS.md` -> `[[../L04_ODE_dualism_alpha]]`
- `audyt/D01_drifting_numbers/NEEDS.md` -> `[[../L05_mass_exponent_drift]]`
- `audyt/D01_drifting_numbers/NEEDS.md` -> `[[../L04_ODE_dualism_alpha]]`
- `audyt/D01_drifting_numbers/README.md` -> `[[../L04_ODE_dualism_alpha]]`
- `audyt/D01_drifting_numbers/README.md` -> `[[../L04_ODE_dualism_alpha]]`
- `audyt/D01_drifting_numbers/README.md` -> `[[../L04_ODE_dualism_alpha]]`
- `audyt/D01_drifting_numbers/README.md` -> `[[../L04_ODE_dualism_alpha]]`
- `audyt/D01_drifting_numbers/README.md` -> `[[../L05_mass_exponent_drift]]`
- `audyt/D01_drifting_numbers/README.md` -> `[[../L04_ODE_dualism_alpha]]`
- `audyt/D01_drifting_numbers/README.md` -> `[[../L05_mass_exponent_drift]]`
- `audyt/EXTERNAL_REVIEW_2026-05-06.md` -> `[[S01_metric_four_forms]]`
- `audyt/EXTERNAL_REVIEW_2026-05-06.md` -> `[[L05_mass_exponent_drift]]`
- `audyt/L01_rho_operational/EXT1_FRW_radiation_era_2026-05-06.md` -> `[[../../research/op-L01-rho-stress-energy-bridge-2026-05-04/]]`
- `audyt/L01_rho_operational/EXT1_FRW_radiation_era_2026-05-06.md` -> `[[../../core/sek08a]]`
- `audyt/L01_rho_operational/NEEDS.md` -> `[[../S04_metric_coupling_axiom]]`
- `audyt/L01_rho_operational/POST_ACTION_UPDATE_2026-05-04.md` -> `[[../../research/op-L01-rho-stress-energy-bridge-2026-05-04]]`
- `audyt/L01_rho_operational/POST_ACTION_UPDATE_2026-05-04.md` -> `[[../../research/op-L01-rho-stress-energy-bridge-2026-05-04]]`
- `audyt/L01_rho_operational/POST_ACTION_UPDATE_2026-05-04.md` -> `[[../../research/op-L01-rho-stress-energy-bridge-2026-05-04/]]`
- `audyt/L01_rho_operational/README.md` -> `[[../S04_metric_coupling_axiom]]`
- `audyt/L01_rho_operational/README.md` -> `[[../S04_metric_coupling_axiom]]`
- `audyt/L01_rho_operational/README.md` -> `[[../S04_metric_coupling_axiom]]`
- `audyt/L03_K_phi_stability/README.md` -> `[[../L04_ODE_dualism_alpha]]`
- `audyt/L04_ODE_dualism_alpha/NEEDS.md` -> `[[../L05_mass_exponent_drift]]`
- `audyt/L04_ODE_dualism_alpha/NEEDS.md` -> `[[../L03_K_phi_stability]]`
- `audyt/L04_ODE_dualism_alpha/NEEDS.md` -> `[[../D01_drifting_numbers]]`
- `audyt/L04_ODE_dualism_alpha/POST_ACTION_UPDATE_2026-05-04.md` -> `[[../../research/op-L04-ODE-canonicalization-2026-05-04]]`
- `audyt/L04_ODE_dualism_alpha/POST_ACTION_UPDATE_2026-05-04.md` -> `[[../../research/op-L04-ODE-canonicalization-2026-05-04]]`
- `audyt/L04_ODE_dualism_alpha/POST_ACTION_UPDATE_2026-05-04.md` -> `[[../../research/op-L04-ODE-canonicalization-2026-05-04/]]`
- `audyt/L04_ODE_dualism_alpha/POST_ACTION_UPDATE_2026-05-04.md` -> `[[../L05_mass_exponent_drift]]`
- `audyt/L04_ODE_dualism_alpha/README.md` -> `[[../L03_K_phi_stability]]`
- `audyt/L04_ODE_dualism_alpha/README.md` -> `[[../L05_mass_exponent_drift]]`
- `audyt/L04_ODE_dualism_alpha/README.md` -> `[[../L05_mass_exponent_drift]]`
- `audyt/L04_ODE_dualism_alpha/README.md` -> `[[../L03_K_phi_stability]]`
- `audyt/L04_ODE_dualism_alpha/README.md` -> `[[../D01_drifting_numbers]]`
- `audyt/L04_ODE_dualism_alpha/README.md` -> `[[../L05_mass_exponent_drift]]`
- `audyt/L04_ODE_dualism_alpha/README.md` -> `[[../L03_K_phi_stability]]`
- `audyt/L04_ODE_dualism_alpha/README.md` -> `[[../D01_drifting_numbers]]`
- `audyt/L05_mass_exponent_drift/README.md` -> `[[../L04_ODE_dualism_alpha]]`
- `audyt/L05_mass_exponent_drift/README.md` -> `[[../L04_ODE_dualism_alpha]]`
- `audyt/L05_mass_exponent_drift/README.md` -> `[[../L04_ODE_dualism_alpha]]`
- ... +605 more

