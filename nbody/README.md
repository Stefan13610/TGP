# TGP N-body package

**Topology-Generated Potential: N-body dynamics, chaos, and analytical structure**

This package implements the complete N-body physics of TGP (Topology-Generated Potential), from the field Lagrangian through pairwise and three-body potentials to equations of motion, Lyapunov chaos analysis, and analytical equilibrium theory.

## Quick start

```bash
# Install dependencies
pip install -r requirements.txt

# Run full regression suite (59 tests, ~2 min)
cd TGP_v1
python -m nbody.examples.verify_all

# Run Lyapunov regression (~3 min with --quick)
python -m nbody.examples.verify_nbody_lyapunov_quick --quick

# Run EOM regression
python -m nbody.examples.verify_nbody_eom_quick --quick
```

All three should exit with code 0 (PASS).

## Requirements

- Python >= 3.10
- numpy >= 1.26
- scipy >= 1.12
- matplotlib >= 3.8 (for plotting scripts only)

See `requirements.txt`.

## Package structure

### Core modules (nbody/)

| Module | Function | Status |
|--------|----------|--------|
| `tgp_field.py` | Yukawa profile, energy density, screening mass m_sp | Stable |
| `pairwise.py` | 2-body potential V_2 (analytical, exact) | Stable |
| `three_body_force_exact.py` | 3-body potential V_3 via Feynman 2D reduction (exact) | Stable |
| `three_body_terms.py` | Triple overlap integral (numerical 3D quadrature) | Stable |
| `multipole_triple_overlap.py` | I_Y multipole Legendre expansion (semi-analytic) | Stable |
| `equilibria.py` | Static equilibria (exact for pairs, approx with 3B) | Stable |
| `stability.py` | Hessian + zero-mode projection + full normal mode analysis | Stable |
| `dynamics_v2.py` | Leapfrog (symplectic) + RK45 integrators | Stable |
| `dynamics_backends.py` | Backend abstraction for integration (Newton/TGP variants) | Stable |
| `lyapunov.py` | Benettin algorithm (RK4/leapfrog+tangent), spectrum QR | Stable |
| `eom_tgp.py` | Vectorial N-body EOM | Stable |
| `nbody_energy.py` | Total energy (V_2 + V_3) | Stable |
| `configurations.py` | IC generators: triangle, tetrahedron, N-gon, Burrau | Stable |
| `soliton_interaction.py` | Classical soliton overlap analysis | Stable |
| `yukawa_from_defect.py` | EFT projection: defect -> C_eff | Stable |

### Example scripts (examples/)

166 active scripts in `examples/ex*.py`, 55 archived in `examples/_archiwum/`.

| Range | Topic | Key results |
|-------|-------|-------------|
| ex55-ex103 | Efimov, rotation curves, V_eff | Historical |
| ex104-ex124 | Soliton, tail, Koide, formalization | Path 9 |
| ex125-ex147 | EOM, regression, 3B forces | Integration |
| ex148-ex194 | Lyapunov, Benettin, matched-H, CSV | Chaos (P1) |
| ex195-ex197 | K_sub full vs LPA, particle predictions | Soliton ODE |
| ex198-ex200 | Beta scan, normalization, V2+V3 scan | Chaos synthesis |
| ex201 | Full Hessian stability scan | Earnshaw (P2) |
| ex202 | Multipole vs Feynman verification | I_Y multipole (P4) |
| ex203-ex204 | Figure-8, Lagrange points, Poincare | Orbits (P3) |
| ex205 | EFT projection, C_eff scaling | Path C (P5) |
| ex206 | Soliton interaction analysis | Soliton (P6) |
| ex207 | Multi-IC beta scan (Burrau, equilateral, random) | IC-dependence |
| ex208 | Full Lyapunov spectrum (k=18=6N) | Hamiltonian check |
| ex209 | V3/V2 regime map vs m_sp*d | V3/V2 (P7) |
| ex210 | Analytical equilibria + Hill criterion | Equilibria (P8) |
| ex211 | Unequal mass, phase diagram, v_esc, breathing | Phase diagram (P9) |

### Regression tests

| Script | Tests | Time |
|--------|-------|------|
| `verify_all.py` | 59 | ~2 min |
| `verify_nbody_lyapunov_quick.py --quick` | ~47 | ~3 min |
| `verify_nbody_eom_quick.py --quick` | ~10 | ~1 min |

## Reproducing publication results

The clean publication document is `tgp_nbody_results_clean.tex` (compile with `pdflatex`, 3 passes).

### Manifest: Result -> Script -> Output

| Result in paper | Script(s) | Key output |
|----------------|-----------|------------|
| **Result 1**: Earnshaw violation | `ex201` | omega^2 > 0 for TGP at d_well (11/11 points) |
| **Result 2**: Chaos suppression 35-45% | `ex200`, `ex207`, `ex208` | `_outputs/ex207_multi_ic_beta_v3.csv` |
| **Result 3**: C_eff from EFT projection | `ex205` | C_eff ~ (1-g0)^1.02 scaling |
| **Result 4**: Classical != Yukawa | `ex206` | delta(d) oscillates, 7 sign changes |
| **Result 5**: I_Y multipole | `ex202` | <1% error at L_max=10 |
| **Result 6**: V3/V2 regime map | `ex209` | `_outputs/ex209_v3v2_vs_msp_d.csv`, `_outputs/ex209_v3v2_vs_beta.csv` |
| **Result 7**: Equilibria + Hill | `ex210` | d_well analytical, R_H(TGP)/R_H(Newton) = 1.0-2.9 |
| **Result 8**: Phase diagram + v_esc + breathing | `ex211` | beta_crit = (3/2)*sqrt(gamma*(C1+C2)); 100% agreement |
| **Supplementary**: Closed orbits | `ex203`, `ex204` | Figure-8 survives; L4 stable |
| **Table 1**: Beta scan (Sec. 3) | `ex200` | 9/10 SUPPRESS with V2+V3 |
| **Table 2**: Multi-IC (Sec. 3) | `ex207` | Burrau 75% SUPPRESS; equilateral 0% |

### Reproducing a single result

```bash
# Example: reproduce Result 2 (chaos suppression)
python -m nbody.examples.ex200_beta_scan_matched_h_v2v3 --quick

# Example: reproduce Result 7 (Hill criterion)
python -m nbody.examples.ex210_analytical_equilibria_hill --quick

# Example: reproduce Result 8 (phase diagram)
python -m nbody.examples.ex211_unequal_mass_phase_diagram --quick
```

All scripts support `--quick` for fast verification (reduced grid/precision).

## Physics summary

TGP generates an effective N-body potential from a single scalar field Lagrangian:

```
V_total = sum_{i<j} V_2(d_ij) + sum_{i<j<k} V_3(d_ij, d_ik, d_jk)
```

where:
- `V_2(d) = -4*pi*C1*C2/d + 8*pi*beta*C1*C2/d^2 - 12*pi*gamma*C1*C2*(C1+C2)/d^3`
- `V_3 = -6*gamma*C1*C2*C3 * I_Y(d12, d13, d23)` (Feynman integral)

Key parameters: `beta` (self-interaction), `gamma` (cubic coupling), `C_i` (source strengths), `m_sp = sqrt(3*gamma - 2*beta)` (screening mass).

### Main results

1. **Earnshaw violation**: TGP admits stable static equilibria (Newton forbids them)
2. **Chaos suppression**: 35-45% reduction in lambda_max for beta > 0.025 (IC-dependent)
3. **EFT closure**: Yukawa coupling C_eff derived from topological defect (not postulated)
4. **Classical != Yukawa**: Soliton interaction has oscillatory tail, Yukawa requires EFT mass gap
5. **Multipole I_Y**: Semi-analytic form with <1% error at L_max=10
6. **V3/V2 regime**: Three-body forces are ~1% in Lyapunov regime, perturbative
7. **Analytical equilibria**: Closed-form d_well, omega^2, Hill sphere
8. **Phase diagram**: beta_crit = (3/2)*sqrt(gamma*(C1+C2)); escape velocity; breathing mode

## LaTeX documents (nbody/)

22 `.tex` files document specific N-body derivations. Key files:
- `tgp_nbody_results_clean.tex` -- publication-ready draft (10 pages, 8 Results)
- `tgp_nbody_lagrangian_eom.tex` -- EOM derivation
- `tgp_lyapunov_benettin.tex` -- Benettin chaos analysis and synthesis

## Analysis documents

- `PLAN_ROZWOJU_NBODY.md` -- Development plan with status of all phases (P1-P9)
- `ANALIZA_NBODY_ROZWOJ.md` -- Full analysis: inventory, consistency, roadmap
- `_archiwum_docs/ODPOWIEDZ_NA_RECENZJE_ZEWNETRZNA.md` -- Response to external review critique (archived)

## Citation

If you use this code, please cite:

```
M. Serafin, "N-body dynamics of Topology-Generated Potential:
Earnshaw violation, chaos suppression, and the EFT origin of Yukawa coupling"
(2026), in preparation.
```
