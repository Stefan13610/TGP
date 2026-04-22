# Publication manifest: N-body dynamics of TGP

**Paper:** "N-body dynamics of Topology-Generated Potential: Earnshaw violation, chaos suppression, and the EFT origin of Yukawa coupling"
**Author:** Mateusz Serafin
**Date:** April 2026

## Files included in submission

### LaTeX source
| File | Description |
|------|-------------|
| `tgp_nbody_results_clean.tex` | Main paper (13 pages, 8 Results) |
| `tgp_nbody.bib` | Bibliography (35 entries) |

### Reproducibility package
| File | Description |
|------|-------------|
| `README.md` | Full package documentation |
| `requirements.txt` | Python dependencies (numpy, scipy, matplotlib) |
| `reproduce_all.py` | Master reproduction script (10/10 stages PASS, ~100s) |

### Result-to-script mapping

| Result | Script | Output | Test |
|--------|--------|--------|------|
| 1 (Earnshaw) | `examples/ex201_full_hessian_stability_scan.py` | Normal-mode eigenvalues | 11/11 configs PASS |
| 2 (Chaos) | `examples/ex200_lyapunov_beta_scan_yukawa.py` | λ_max(β) scan | 35-45% suppression |
| 2 (Chaos) | `examples/ex207_lyapunov_multi_ic_beta_scan_v3.py` | Multi-IC scan | IC-dependent |
| 3 (C_eff) | `examples/ex205_path_c_yukawa_from_defect.py` | EFT projection | C_eff formula |
| 4 (Classical≠Yukawa) | `examples/ex206_soliton_interaction_energy.py` | Oscillatory vs Yukawa | Verified |
| 5 (I_Y multipole) | `examples/ex202_multipole_vs_feynman_triple_overlap.py` | Convergence | <0.01% at L=15 |
| 6 (V₃/V₂ map) | `examples/ex209_v3_v2_regime_msp_scan.py` | Regime diagram | Perturbative ~1% |
| 7 (Equilibria) | `examples/ex210_analytical_equilibria_hill.py` | d_rep, d_well, ω², R_H | Closed-form verified |
| 8 (Phase diagram) | `examples/ex211_unequal_mass_phase_diagram.py` | β_crit boundary | 100% agreement |

### Regression tests
| Suite | Scripts | Status |
|-------|---------|--------|
| `verify_all.py` | 59 tests | 59/59 PASS |
| `verify_nbody_eom_quick.py` | EOM + conservation | PASS |
| `verify_nbody_lyapunov_quick.py` | Lyapunov convergence | PASS |

### Core modules (14 files)
- `tgp_field.py`, `pairwise.py`, `three_body_force_exact.py`
- `three_body_terms.py`, `multipole_triple_overlap.py`
- `equilibria.py`, `stability.py`
- `dynamics_v2.py`, `dynamics_backends.py`, `eom_tgp.py`
- `lyapunov.py`, `tgp_nbody_lagrangian_eom.py`
- `soliton_interaction.py`, `__init__.py`

## Reproduction instructions

```bash
# 1. Install dependencies
pip install numpy>=1.26 scipy>=1.12 matplotlib>=3.8

# 2. Navigate to TGP root
cd TGP_v1

# 3. Run master reproduction (10 stages, ~100s)
python -m nbody.reproduce_all --quick

# 4. Run full regression (59 tests, ~2 min)
python -m nbody.examples.verify_all

# Expected output: 10/10 PASS, 59/59 PASS
```

## System requirements
- Python >= 3.10
- OS: tested on Windows 11, should work on Linux/macOS
- RAM: < 500 MB
- Time: ~100s for quick reproduction, ~5 min for full regression
