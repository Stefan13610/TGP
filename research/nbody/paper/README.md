# research/nbody/paper — fragmenty publikacji N-body

Zrodla LaTeX i artefakty publikacyjne dla monografii N-body TGP.

## Dokumenty glowne (self-contained)

| Plik | Rola |
|------|------|
| `tgp_nbody_results_clean.tex` + `.pdf` | **Glowny draft** (10 str., 8 Results) |
| `tgp_yukawa_exact_reduction.tex` + `.pdf` | Redukcja Yukawa (agreguje ~18 fragmentow) |
| `supplementary_material.tex` | Material uzupelniajacy do publikacji |
| `tgp_nbody.bib` | Bibliografia BibTeX |

## Fragmenty `\input{}`-owane przez **monografie** (`main.tex`)

Pliki ponizej sa wciagane do glownej monografii TGP:

| Plik | Wciagany przez |
|------|----------------|
| `tgp_cosmo_likelihood.tex` | `core/sek07_predykcje/` |
| `tgp_gw_summary_table.tex` | `core/sek07_predykcje/` |
| `tgp_metrology_observables.tex` | `core/sek04_stale/` |
| `tgp_time_glossary.tex` | `core/sek04_stale/` |
| `tgp_perturbations_lss.tex` | `core/sek08_formalizm/` |
| `tgp_lensing_geff.tex` | `core/sek08_formalizm/` |
| `tgp_metric_bridge_table.tex` | `core/sek08c_metryka_z_substratu/` |
| `tgp_ppn_full.tex` | `core/sek08c_metryka_z_substratu/` |
| `tgp_nbody_lagrangian_eom.tex` | `partial_proofs/trojcialowe_nbody/dodatekY_nbody.tex` |

## Fragmenty agregowane przez `tgp_yukawa_exact_reduction.tex`

`\input{}`-owane lokalnie (bez sciezki, rozwiazywane w tym samym folderze):

- `tgp_yukawa_IY_multipole_gegenbauer`, `tgp_nonstationary_equations`,
  `tgp_nbody_predictions`, `tgp_topological_defect`, `tgp_false_vacuum`,
  `tgp_ncrit_exact`, `tgp_scattering`, `tgp_bound_state`, `tgp_efimov_analog`,
  `tgp_optimal_geometry`, `tgp_quantum_window`, `tgp_quantum_window_msp`,
  `tgp_polygon_geometry`, `tgp_wavefunction`, `tgp_physical_interpretation`,
  `tgp_triangle_stability`, `tgp_open_questions`, `tgp_rotation_curve`,
  `tgp_2d_isosceles`

## Samodzielne fragmenty (nie zintegrowane jeszcze)

- `tgp_lyapunov_benettin.tex` — analiza chaosu (Benettin)
- `tgp_yukawa_fourier_feynman_reduction.tex` — alternatywna redukcja
- `tgp_nbody_special_configurations.tex` — specjalne konfiguracje N-body

## Kompilacja

Z folderu `paper/`:

```bash
pdflatex tgp_nbody_results_clean.tex
bibtex   tgp_nbody_results_clean
pdflatex tgp_nbody_results_clean.tex
pdflatex tgp_nbody_results_clean.tex
```

Pliki `.aux/.log/.out/.toc/.bbl/.blg` sa ignorowane przez `.gitignore`.

## Kontekst

Patrz [[../README.md]] dla pelnego opisu pakietu N-body TGP.
