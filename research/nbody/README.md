# TGP N-body package

**Topology-Generated Potential: N-body dynamics, chaos, and analytical structure**

This package implements the synchronized `nbody` layer of TGP: the bridge from
the classical defect sector to effective Yukawa sources, then to pairwise and
three-body EOM, integrators, and chaos diagnostics.

## Quick start

```bash
# From the repo root:
python -m nbody.examples.verify_nbody_canonical_quick
python -m nbody.examples.verify_nbody_bridge_extended
python -m nbody.examples.verify_nbody_eom_quick
```

Read first:

- `nbody/THEORY_SYNC_NBODY.md`
- `nbody/examples/STATUS_MAP.md`

Recommended interpretation order:

1. `CLASSICAL`: defect / soliton with `sin(r)/r` tail
2. `EFT-DERIVED`: `C_eff`, `m_sp`, bridge from defect to Yukawa source
3. `N-BODY`: effective `V_2 + V_3`, EOM, integrators, Lyapunov diagnostics

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
| `tgp_field.py` | Field profiles, Yukawa helpers, screening mass `m_sp` | Stable |
| `bridge_nbody.py` | Canonical bridge: defect -> `C_eff` -> Yukawa-source inputs | Stable |
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

Do not treat all `exNNN` scripts equally. The canonical map lives in:

- `examples/STATUS_MAP.md`
- `examples/README.md`

| Entry point | Role | Notes |
|-------------|------|-------|
| `verify_nbody_canonical_quick.py` | Short canonical route | FULL/LPA sync -> EFT bridge -> EOM |
| `verify_nbody_bridge_extended.py` | Stable extended route | canonical route + representative Lyapunov checks |
| `verify_nbody_eom_quick.py` | EOM/backend regression | force/potential consistency, COM/P/L checks |
| `verify_nbody_lyapunov_quick.py` | Broad Lyapunov package | larger historical P1 sweep |
| `ex195`-`ex197` | Canonical soliton synchronization | FULL `K_sub = g^2` |
| `ex205` | Canonical EFT bridge | defect -> `C_eff` |
| `ex55`-`ex103` | Historical layer | mostly `LEGACY-TRANSLATIONAL` |
| `ex104+` | Mixed active exploration | see `STATUS_MAP.md` |

### Regression tests

| Script | Tests | Time |
|--------|-------|------|
| `verify_nbody_canonical_quick.py` | canonical chain | short |
| `verify_nbody_bridge_extended.py` | canonical + representative chaos checks | moderate |
| `verify_nbody_eom_quick.py` | EOM/backend regression | short |

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

`nbody` uses the synchronized three-layer interpretation:

- classical defect: oscillatory tail `sin(r)/r`
- EFT bridge: projection to `C_eff` and `m_sp`
- effective N-body layer: point-source `V_2 + V_3`

The effective potential is

```text
V_total = sum_{i<j} V_2(d_ij) + sum_{i<j<k} V_3(d_ij, d_ik, d_jk)
```

where:

- `V_2(d) = -4*pi*C1*C2/d + 8*pi*beta*C1*C2/d^2 - 12*pi*gamma*C1*C2*(C1+C2)/d^3`
- `V_3 = (2*beta - 6*gamma)*C1*C2*C3 * I_Y(d12, d13, d23)`

Key parameters: `beta`, `gamma`, `C_i`, `m_sp = sqrt(3*gamma - 2*beta)`.

Exact Yukawa-overlap scaling identity:

```text
(d12*∂/∂d12 + d13*∂/∂d13 + d23*∂/∂d23 - m_sp*∂/∂m_sp) I_Y = 0
```

This identity is now exposed numerically in `three_body_force_exact.py`.

Useful canonical shape variables:

```text
d_min <= d_mid <= d_max
q1 = d_min / d_max
q2 = d_mid / d_max
t  = m_sp * d_max
I_Y = F(t; q1, q2)
```

For general shape, the large-`t` exponential rate is controlled by

```text
lambda(q1,q2) = min_{alpha in Delta_2} sqrt(Q/Delta)
I_Y(t;q1,q2) ~ exp(-lambda(q1,q2) * t) * [prefactor + corrections]
```

Controlled equilateral large-`t` asymptotic:

```text
I_Y^eq(t) ~ A_eq * t^(-3/2) * exp(-sqrt(3) * t)
A_eq = 4*sqrt(2)*pi^(3/2) / 3^(3/4)
```

First subleading correction:

```text
I_Y^eq(t) ~ A_eq * t^(-3/2) * exp(-sqrt(3) * t)
            * (1 - 5/(8*sqrt(3)*t) + O(t^(-2)))
```

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
