# TGP Numerical Verification Scripts

Numerical verification package for **Teoria Generowanej Przestrzeni (TGP)**
by Mateusz Serafin (March 2026).

## Quick Start

```bash
# Install dependencies
pip install numpy scipy

# Run core verification (~60s)
python run_all.py

# Run quick check (~30s)
python run_all.py --quick

# Run all scripts (~5 min)
python run_all.py --full
```

Requires: Python 3.9+, numpy, scipy.

## Script Categories

### Core Verification (7 scripts)
Key predictions — these are the strongest quantitative tests of TGP.

| Script | Tests | Key Result |
|--------|-------|------------|
| `ex147` | Substrate ODE soliton | g(r) solution, oscillatory tail |
| `ex157` | One-parameter prediction | g0^e=0.869 -> r21, m_mu, m_tau |
| `ex174` | Full predictions summary | 11 observables, all within 1sigma |
| `ex184` | Discrete running | alpha_s(N_f=3,5) from same g0^e |
| `ex186` | Mass-coupling unification | r21 -> g0^e -> alpha_s (0 free params) |
| `ex187` | Cosmology confrontation | H(z), w(z), G(z) vs Planck/DESI/BBN |
| `ex165` | Slow-roll inflation | n_s=0.965, r=0.003 |

### Supporting Verification (16 scripts)
Consistency checks, gauge sector, ERG stability, Koide relations, PPN.

### Exploratory (3 scripts)
Multi-loop running comparison, Phi_0 exact value, alpha_s variants.

## Key Results

### Particle Sector (strongest argument)
- **Lepton masses from 1 parameter**: g0^e = 0.869 (calibrated from r21 = m_mu/m_e)
  - r21 = 206.77 (PDG: 206.768, 0.001%)
  - m_tau = 1776.9 MeV (PDG: 1776.86, 0.002%)
  - Koide K = 2/3 exact
- **Strong coupling constant (0 free parameters)**:
  - alpha_s(M_Z) = 0.1174 (PDG: 0.1179 +/- 0.0009, 0.6 sigma)
  - alpha_s(m_tau) = 0.326 (PDG: 0.330 +/- 0.014, 0.3 sigma)
  - Ratio: (5/3)^2 = 2.778 (PDG: 2.799 +/- 0.121, 0.18 sigma)

### Cosmological Sector (9/9 PASS)
- H(z)/H_LCDM < 1.5% deviation (within DESI BAO limits)
- w_DE = -1 + O(10^-9) (indistinguishable from cosmological constant)
- G(z_BBN)/G_0 within 3% (BBN constraint: <10%)
- n_s = 0.965 (Planck: 0.9649, 0.03 sigma)

### Master Verification (ex194)
45/45 tests PASS across all modules (ex190–ex193). GO status.
- **28 PASS** from 2 free parameters (Phi_eff, g0^e)
- Efficiency: N_pass/N_param = 14.0

## Script Naming Convention

- `ex1XX_name.py` — numbered experiments (chronological order)
- `tgp_name.py` — thematic utility scripts
- `check_refs.py` — LaTeX reference checker

## Correspondence to Manuscript

| Manuscript Section | Verification Scripts |
|-------------------|---------------------|
| Soliton ODE (Sec. 8) | ex147, ex147d, ex150 |
| Lepton masses (App. J2, T) | ex153, ex155, ex157 |
| alpha_s (App. V) | ex175, ex178, ex184, ex186 |
| Koide relations (App. T2-T4) | ex168, ex153 |
| Inflation (Sec. 7) | ex165 |
| Cosmology (Sec. 8, App. H) | ex187 |
| PPN (Sec. 8) | ex167 |
| ERG/stability (Sec. 8b) | ex141, ex142, ex146, ex148 |
| Gauge sector (Sec. 9, App. O,U,V) | ex140, ex145, ex183 |
| Dimensional argument k=4 (App. R) | ex188 |
| Entropy chain Z₃ (App. R2) | ex189 |
| Full consistency chain | ex190 (9/9 PASS) |
| Quark confinement m₀ (App. X) | ex191 |
| Cosmology vs DESI DR1 (Sec. 8a) | ex192 (5/5 PASS) |
| Unified prediction table | ex193 (28 PASS / 1 CAL / 1 PRED) |
| **Master runner (all above)** | **ex194 (45/45 PASS → GO)** |
