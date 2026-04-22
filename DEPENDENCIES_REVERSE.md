# Reverse Dependencies - TGP/TGP_v1

> Generated: 2026-04-22 by `tooling/build_deps_graph.py`
> Sources: \input{} + \ref{} + \cite{} + [[wikilink]]

## Coarse view

### <root>
- depended on by: --

### axioms
- depended on by: <root>, core, core/_meta_latex, core/formalizm, partial_proofs, research

### core
- depended on by: <root>, axioms, core/_meta_latex, core/formalizm, papers_external, partial_proofs, research

### core/_meta_latex
- depended on by: <root>, core

### core/formalizm
- depended on by: <root>, axioms, core, core/_meta_latex, papers_external, partial_proofs, research

### meta
- depended on by: --

### papers_external
- depended on by: <root>, partial_proofs, research

### partial_proofs
- depended on by: <root>, axioms, core, core/_meta_latex, core/formalizm, research

### research
- depended on by: core, partial_proofs

### tooling
- depended on by: research

## Fine view (per subfolder)

### <root>
- depended on by: --

### axioms
- depended on by: --

### axioms/notacja
- depended on by:
  - <root> (2x | input:2)
  - axioms/roznica_N0 (2x | ref:2)
  - core/_meta_latex (2x | ref:2)
  - partial_proofs/bh_ringdown (2x | ref:2)
  - partial_proofs/hierarchia_mas (2x | ref:2)
  - research (1x | wikilink:1)
  - research/atomic_shells_closure (3x | wikilink:3)
  - research/particle_sector_closure (5x | wikilink:5)

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
  - core/sek06_czarne_dziury (1x | ref:1)
  - core/sek08_formalizm (24x | ref:24)
  - core/sek08b_ghost_resolution (5x | ref:5)
  - core/sek08c_metryka_z_substratu (2x | ref:2)
  - core/sek09_cechowanie (1x | ref:1)
  - core/sek10_N0_wyprowadzenie (3x | ref:3)
  - partial_proofs/bh_ringdown (2x | ref:2)
  - partial_proofs/hierarchia_mas (1x | ref:1)
  - partial_proofs/superconductivity (1x | ref:1)
  - partial_proofs/wielki_wybuch (8x | ref:8)

### core
- depended on by: --

### core/_meta_latex
- depended on by:
  - <root> (2x | input:2)
  - core/sek00_summary (3x | input:1, ref:2)
  - core/sek08_formalizm (2x | ref:2)

### core/formalizm
- depended on by:
  - <root> (10x | input:10)
  - axioms/roznica_N0 (1x | ref:1)
  - axioms/substrat (4x | ref:4)
  - core/_meta_latex (8x | ref:8)
  - core/sek00_summary (2x | ref:2)
  - core/sek03_rezimy (1x | ref:1)
  - core/sek07_predykcje (14x | ref:14)
  - core/sek07a_wymiar_wzmocniony (2x | ref:2)
  - core/sek08_formalizm (6x | ref:6)
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

### core/sek00_summary
- depended on by:
  - <root> (1x | input:1)

### core/sek01_ontologia
- depended on by:
  - <root> (1x | input:1)
  - axioms/notacja (9x | ref:9)
  - axioms/roznica_N0 (5x | ref:5)
  - axioms/substrat (7x | ref:7)
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
  - axioms/notacja (5x | ref:5)
  - core/_meta_latex (7x | ref:7)
  - core/formalizm (12x | ref:12)
  - core/sek00_summary (6x | ref:6)
  - core/sek01_ontologia (2x | ref:2)
  - core/sek02_pole (2x | ref:2)
  - core/sek06_czarne_dziury (12x | ref:12)
  - core/sek07_predykcje (12x | ref:12)
  - core/sek07a_wymiar_wzmocniony (1x | ref:1)
  - core/sek08_formalizm (24x | ref:24)
  - core/sek08c_metryka_z_substratu (4x | ref:4)
  - core/sek09_cechowanie (1x | ref:1)
  - papers_external/tgp_english_summary (7x | ref:7)
  - partial_proofs/bh_ringdown (3x | ref:3)
  - partial_proofs/trojcialowe_nbody (2x | ref:2)
  - partial_proofs/wielki_wybuch (2x | ref:2)

### core/sek05_ciemna_energia
- depended on by:
  - <root> (1x | input:1)
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
  - research/nbody (2x | ref:2)

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
  - axioms/notacja (33x | ref:33)
  - axioms/roznica_N0 (1x | ref:1)
  - axioms/substrat (15x | ref:15)
  - core/_meta_latex (41x | ref:41)
  - core/formalizm (21x | ref:21)
  - core/sek00_summary (13x | ref:13)
  - core/sek01_ontologia (21x | ref:21)
  - core/sek02_pole (9x | ref:9)
  - core/sek03_rezimy (23x | ref:23)
  - core/sek04_stale (11x | ref:11)
  - core/sek05_ciemna_energia (11x | ref:11)
  - core/sek06_czarne_dziury (7x | ref:7)
  - core/sek07_predykcje (51x | ref:51)
  - core/sek08a_akcja_zunifikowana (24x | ref:24)
  - core/sek08b_ghost_resolution (8x | ref:8)
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
  - research/nbody (12x | ref:12)

### core/sek08a_akcja_zunifikowana
- depended on by:
  - <root> (1x | input:1)
  - axioms/notacja (3x | ref:3)
  - core/_meta_latex (3x | ref:3)
  - core/formalizm (3x | ref:3)
  - core/sek00_summary (3x | ref:3)
  - core/sek06_czarne_dziury (1x | ref:1)
  - core/sek08_formalizm (12x | ref:12)
  - core/sek08c_metryka_z_substratu (4x | ref:4)
  - core/sek09_cechowanie (1x | ref:1)
  - partial_proofs/most_gamma_phi (1x | ref:1)
  - research/nbody (1x | ref:1)

### core/sek08b_ghost_resolution
- depended on by:
  - <root> (1x | input:1)
  - core/_meta_latex (3x | ref:3)
  - core/formalizm (3x | ref:3)
  - core/sek00_summary (1x | ref:1)
  - core/sek07_predykcje (1x | ref:1)
  - core/sek09_cechowanie (1x | ref:1)
  - core/sek10_N0_wyprowadzenie (3x | ref:3)
  - partial_proofs/chiralnosc (1x | ref:1)
  - partial_proofs/hierarchia_mas (3x | ref:3)

### core/sek08c_metryka_z_substratu
- depended on by:
  - <root> (1x | input:1)
  - core/_meta_latex (2x | ref:2)
  - core/formalizm (4x | ref:4)
  - core/sek00_summary (3x | ref:3)
  - core/sek02_pole (1x | ref:1)
  - core/sek04_stale (2x | ref:2)
  - core/sek08_formalizm (3x | ref:3)
  - research/nbody (3x | ref:3)

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
- depended on by: --

### papers_external
- depended on by: --

### papers_external/arxiv_submission
- depended on by:
  - research/nbody (24x | ref:24)

### papers_external/paper_bh_shadow
- depended on by:
  - <root> (2x | ref:2)
  - papers_external/paper_lepton_masses (1x | ref:1)
  - papers_external/tgp_core_paper (6x | ref:6)
  - papers_external/tgp_english_summary (3x | ref:3)
  - papers_external/tgp_sc_paper (3x | ref:3)
  - research (1x | wikilink:1)

### papers_external/paper_lepton_masses
- depended on by:
  - papers_external/tgp_sc_paper (3x | ref:3)

### papers_external/tgp_core_paper
- depended on by:
  - papers_external/tgp_english_summary (2x | ref:2)
  - research (1x | wikilink:1)
  - research/atomic_shells_closure (1x | wikilink:1)
  - research/casimir_mof (1x | wikilink:1)
  - research/liquid_viscosity (1x | wikilink:1)
  - research/muon_g_minus_2 (1x | wikilink:1)
  - research/neutrino_msw (1x | wikilink:1)
  - research/thermal_transport_molecular (1x | wikilink:1)

### papers_external/tgp_english_summary
- depended on by:
  - papers_external/paper_bh_shadow (1x | ref:1)
  - partial_proofs/hierarchia_mas (1x | ref:1)
  - partial_proofs/zero_mode (1x | ref:1)

### papers_external/tgp_sc_paper
- depended on by:
  - <root> (2x | ref:2)
  - research (1x | wikilink:1)
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
  - research (1x | wikilink:1)

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
  - research (5x | wikilink:5)
  - research/atom_from_soliton (1x | wikilink:1)
  - research/atomic_shells_closure (1x | wikilink:1)
  - research/casimir_mof (1x | wikilink:1)
  - research/cohesion_closure (1x | wikilink:1)
  - research/em_from_substrate (1x | wikilink:1)
  - research/liquid_viscosity (1x | wikilink:1)
  - research/muon_g_minus_2 (1x | wikilink:1)
  - research/thermal_transport_molecular (1x | wikilink:1)

### partial_proofs/hierarchia_mas
- depended on by:
  - <root> (6x | input:4, ref:2)
  - axioms/notacja (1x | ref:1)
  - core/_meta_latex (15x | ref:15)
  - core/formalizm (18x | ref:18)
  - core/sek04_stale (1x | ref:1)
  - core/sek07_predykcje (3x | ref:3)
  - core/sek08_formalizm (2x | ref:2)
  - core/sek08b_ghost_resolution (11x | ref:11)
  - core/sek09_cechowanie (5x | ref:5)
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
  - research/nbody (1x | ref:1)

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
  - research/atomic_shells_closure (1x | wikilink:1)
  - research/casimir_mof (1x | wikilink:1)
  - research/desi_dark_energy (4x | wikilink:4)
  - research/galaxy_scaling (1x | wikilink:1)
  - research/liquid_viscosity (1x | wikilink:1)
  - research/muon_g_minus_2 (1x | wikilink:1)
  - research/particle_sector_closure (4x | wikilink:4)

### research/atom_from_soliton
- depended on by:
  - partial_proofs/fermion_from_soliton (2x | wikilink:2)
  - partial_proofs/nuclear_from_soliton (2x | wikilink:2)

### research/atomic_shells_closure
- depended on by:
  - research/atom_from_soliton (2x | wikilink:2)
  - research/cohesion_closure (2x | wikilink:2)

### research/brannen_sqrt2
- depended on by:
  - research/particle_sector_closure (1x | wikilink:1)

### research/cabibbo_correction
- depended on by: --

### research/casimir_mof
- depended on by:
  - research/muon_g_minus_2 (3x | wikilink:3)

### research/cohesion_closure
- depended on by:
  - partial_proofs/nuclear_from_soliton (1x | wikilink:1)
  - research/atom_from_soliton (1x | wikilink:1)

### research/continuum_limit
- depended on by: --

### research/cosmo_tensions
- depended on by:
  - research/galaxy_scaling (2x | wikilink:2)

### research/desi_dark_energy
- depended on by: --

### research/em_from_substrate
- depended on by:
  - research/atom_from_soliton (5x | wikilink:5)
  - research/cohesion_closure (2x | wikilink:2)

### research/em_from_substrate/TGP/TGP_v1/research/em_from_substrate
- depended on by:
  - research/em_from_substrate (1x | wikilink:1)

### research/galaxy_scaling
- depended on by:
  - research (4x | wikilink:4)
  - research/desi_dark_energy (3x | wikilink:3)

### research/hubble_tension
- depended on by: --

### research/liquid_viscosity
- depended on by: --

### research/mass_scaling_k4
- depended on by: --

### research/metric_ansatz
- depended on by: --

### research/muon_g_minus_2
- depended on by: --

### research/nbody
- depended on by:
  - core/sek02_pole (1x | ref:1)
  - core/sek04_stale (2x | input:2)
  - core/sek07_predykcje (2x | input:2)
  - core/sek08_formalizm (2x | input:2)
  - core/sek08c_metryka_z_substratu (2x | input:2)
  - partial_proofs/fermion_from_soliton (1x | wikilink:1)
  - partial_proofs/nuclear_from_soliton (11x | wikilink:11)
  - partial_proofs/trojcialowe_nbody (1x | input:1)

### research/nbody/examples
- depended on by:
  - research/em_from_substrate (1x | wikilink:1)

### research/nbody/examples/_outputs
- depended on by: --

### research/neutrino_msw
- depended on by: --

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
- depended on by: --

### tooling
- depended on by: --

### tooling/scripts
- depended on by:
  - research/em_from_substrate (1x | wikilink:1)

### tooling/scripts/gauge
- depended on by:
  - research/em_from_substrate (2x | wikilink:2)

