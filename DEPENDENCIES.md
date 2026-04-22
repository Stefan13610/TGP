# Dependencies - TGP/TGP_v1 (forward graph)

> Generated: 2026-04-22 by `tooling/build_deps_graph.py`
> Sources: \input{} + \ref{} + \cite{} + [[wikilink]]

## Summary

- Folders analyzed (fine): 91
- Folders analyzed (coarse): 10
- Total dependencies found: 1657
  - `\input`  edges: 69
  - `\ref`    edges: 1436
  - `[[wiki]]` edges: 152
- `\cite{}` usages counted (bib keys): 100
- Bibliography keys in tgp_main.bib: 25

- Orphans (unresolved): \input=0, \ref=0, wikilink=0

## Coarse view (top-level folders)

### <root>
- depends on: axioms, core, core/_meta_latex, core/formalizm, papers_external, partial_proofs

### axioms
- depends on: core, core/formalizm, partial_proofs

### core
- depends on: axioms, core/_meta_latex, core/formalizm, partial_proofs, research

### core/_meta_latex
- depends on: axioms, core, core/formalizm, partial_proofs

### core/formalizm
- depends on: axioms, core, partial_proofs

### meta
- depends on: --

### papers_external
- depends on: core, core/formalizm

### partial_proofs
- depends on: axioms, core, core/formalizm, papers_external, research

### research
- depends on: axioms, core, core/formalizm, papers_external, partial_proofs, tooling

### tooling
- depends on: --

## Fine view (per subfolder)

### <root>
- `\input`:
  - axioms/notacja (2x)
  - axioms/roznica_N0 (1x)
  - axioms/substrat (1x)
  - core/_meta_latex (2x)
  - core/formalizm (10x)
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
  - --
- Cite count: 6

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
  - core/sek08_formalizm (33x)
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
  - core/sek01_ontologia (7x)
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
  - core/formalizm (8x)
  - core/sek03_rezimy (1x)
  - core/sek04_stale (7x)
  - core/sek05_ciemna_energia (1x)
  - core/sek07_predykcje (2x)
  - core/sek08_formalizm (41x)
  - core/sek08a_akcja_zunifikowana (3x)
  - core/sek08b_ghost_resolution (3x)
  - core/sek08c_metryka_z_substratu (2x)
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
  - core/sek01_ontologia (8x)
  - core/sek02_pole (4x)
  - core/sek03_rezimy (11x)
  - core/sek04_stale (12x)
  - core/sek05_ciemna_energia (3x)
  - core/sek07_predykcje (2x)
  - core/sek07a_wymiar_wzmocniony (1x)
  - core/sek08_formalizm (21x)
  - core/sek08a_akcja_zunifikowana (3x)
  - core/sek08b_ghost_resolution (3x)
  - core/sek08c_metryka_z_substratu (4x)
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
  - core/formalizm (2x)
  - core/sek01_ontologia (4x)
  - core/sek02_pole (2x)
  - core/sek04_stale (6x)
  - core/sek05_ciemna_energia (1x)
  - core/sek06_czarne_dziury (1x)
  - core/sek07_predykcje (2x)
  - core/sek08_formalizm (13x)
  - core/sek08a_akcja_zunifikowana (3x)
  - core/sek08b_ghost_resolution (1x)
  - core/sek08c_metryka_z_substratu (3x)
  - core/sek09_cechowanie (1x)
  - core/sek10_N0_wyprowadzenie (4x)
  - partial_proofs/bh_ringdown (4x)
  - partial_proofs/wielki_wybuch (1x)
- Wikilinks to:
  - --
- Cite count: 0

### core/sek01_ontologia
- `\input`:
  - --
- `\ref` to other folders:
  - axioms/roznica_N0 (5x)
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
  - core/sek01_ontologia (1x)
  - core/sek08_formalizm (11x)
  - core/sek08c_metryka_z_substratu (2x)
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
  - core/sek08_formalizm (51x)
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
  - axioms/substrat (24x)
  - core/_meta_latex (2x)
  - core/formalizm (6x)
  - core/sek01_ontologia (52x)
  - core/sek02_pole (5x)
  - core/sek03_rezimy (11x)
  - core/sek04_stale (24x)
  - core/sek05_ciemna_energia (4x)
  - core/sek06_czarne_dziury (5x)
  - core/sek07_predykcje (7x)
  - core/sek08a_akcja_zunifikowana (12x)
  - core/sek08c_metryka_z_substratu (3x)
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
  - core/sek02_pole (4x)
  - core/sek08_formalizm (24x)
- Wikilinks to:
  - --
- Cite count: 0

### core/sek08b_ghost_resolution
- `\input`:
  - --
- `\ref` to other folders:
  - axioms/substrat (5x)
  - core/sek02_pole (1x)
  - core/sek08_formalizm (8x)
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
  - core/sek04_stale (4x)
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
  - core/sek01_ontologia (4x)
  - core/sek02_pole (4x)
  - core/sek03_rezimy (3x)
  - core/sek04_stale (1x)
  - core/sek07_predykcje (4x)
  - core/sek08_formalizm (10x)
  - core/sek08a_akcja_zunifikowana (1x)
  - core/sek08b_ghost_resolution (1x)
  - partial_proofs/chiralnosc (4x)
  - partial_proofs/defect_hierarchy (10x)
  - partial_proofs/hierarchia_mas (5x)
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
  - --
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
  - axioms/notacja (1x)
  - papers_external/paper_bh_shadow (1x)
  - papers_external/tgp_core_paper (1x)
  - papers_external/tgp_sc_paper (1x)
  - partial_proofs/bh_ringdown (1x)
  - partial_proofs/fermion_from_soliton (5x)
  - research/galaxy_scaling (4x)
- Cite count: 0

### research/atom_from_soliton
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - core/formalizm (1x)
  - core/sek09_cechowanie (1x)
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
  - axioms/notacja (3x)
  - papers_external/tgp_core_paper (1x)
  - partial_proofs/fermion_from_soliton (1x)
  - research (1x)
  - research/superconductivity_closure (2x)
- Cite count: 0

### research/brannen_sqrt2
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - --
- Cite count: 0

### research/cabibbo_correction
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - --
- Cite count: 0

### research/casimir_mof
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - papers_external/tgp_core_paper (1x)
  - partial_proofs/fermion_from_soliton (1x)
  - research (1x)
  - research/superconductivity_closure (1x)
- Cite count: 0

### research/cohesion_closure
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
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
  - --
- Cite count: 0

### research/cosmo_tensions
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - --
- Cite count: 0

### research/desi_dark_energy
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - research (4x)
  - research/galaxy_scaling (3x)
- Cite count: 0

### research/em_from_substrate
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - core/formalizm (1x)
  - partial_proofs/fermion_from_soliton (1x)
  - research/em_from_substrate/TGP/TGP_v1/research/em_from_substrate (1x)
  - research/nbody/examples (1x)
  - tooling/scripts (1x)
  - tooling/scripts/gauge (2x)
- Cite count: 0

### research/em_from_substrate/TGP/TGP_v1/research/em_from_substrate
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - --
- Cite count: 0

### research/galaxy_scaling
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - research (1x)
  - research/cosmo_tensions (2x)
- Cite count: 0

### research/hubble_tension
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - --
- Cite count: 0

### research/liquid_viscosity
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - papers_external/tgp_core_paper (1x)
  - partial_proofs/fermion_from_soliton (1x)
  - research (1x)
  - research/superconductivity_closure (1x)
- Cite count: 0

### research/mass_scaling_k4
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - --
- Cite count: 0

### research/metric_ansatz
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - --
- Cite count: 0

### research/muon_g_minus_2
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - papers_external/tgp_core_paper (1x)
  - partial_proofs/fermion_from_soliton (1x)
  - research (1x)
  - research/casimir_mof (3x)
  - research/superconductivity_closure (1x)
- Cite count: 0

### research/nbody
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - --
- Cite count: 0

### research/nbody/docs
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - axioms/notacja (1x)
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
  - core/sek06_czarne_dziury (2x)
  - core/sek08_formalizm (12x)
  - core/sek08a_akcja_zunifikowana (1x)
  - core/sek08c_metryka_z_substratu (3x)
  - papers_external/arxiv_submission (24x)
  - partial_proofs/trojcialowe_nbody (1x)
- Wikilinks to:
  - axioms/notacja (1x)
- Cite count: 3

### research/neutrino_msw
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - papers_external/tgp_core_paper (1x)
- Cite count: 0

### research/particle_sector_closure
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - axioms/notacja (5x)
  - research (4x)
  - research/brannen_sqrt2 (1x)
- Cite count: 0

### research/qm_born_rule
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - --
- Cite count: 0

### research/qm_decoherence
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - --
- Cite count: 0

### research/qm_entanglement
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - --
- Cite count: 0

### research/qm_foundations
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - --
- Cite count: 0

### research/qm_measurement
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - --
- Cite count: 0

### research/qm_spin
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - --
- Cite count: 0

### research/qm_statistics
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - --
- Cite count: 0

### research/qm_superposition
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - --
- Cite count: 0

### research/rho_normal_state_closure
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - --
- Cite count: 0

### research/s8_tension
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - --
- Cite count: 0

### research/superconductivity_closure
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - papers_external/tgp_sc_paper (1x)
  - research/rho_normal_state_closure (28x)
- Cite count: 0

### research/thermal_transport_molecular
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
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
  - --
- Cite count: 0

### research/why_n3
- `\input`:
  - --
- `\ref` to other folders:
  - --
- Wikilinks to:
  - --
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

