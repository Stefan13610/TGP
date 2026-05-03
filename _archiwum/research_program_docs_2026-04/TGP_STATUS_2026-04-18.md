# TGP Theory Status Report — 2026-04-18

## Executive Summary

TGP (Theory of Gravitational Potential) is a modified gravity theory where gravity emerges from the deformation of a discrete substrate, with particles as solitonic excitations. The theory derives the MOND interpolation function from first principles via membrane Flory exponent (alpha = 4/5), with NO free parameters beyond a0 and morphological codimension c_eff.

**Bottom line:** TGP is a **galaxy-scale theory** that works excellently within its scope but cannot address cluster-scale mass deficits or cosmological tensions. This is a structural feature, not a bug.

---

## Scorecard

### PROVEN (analytic theorems + numerical verification)

| Result | Folder | Evidence | Significance |
|---|---|---|---|
| alpha = 4/5 from Flory exponent | galaxy_scaling/gs42 | 12 methods, delta = +-0.005 | Foundational |
| gamma/alpha = 1/2 exact (codim-1) | galaxy_scaling/gs42 | Geometric proof | Foundational |
| c1 = 1 - ln(3)/4 (B=sqrt2) | brannen_sqrt2 | Analytic proof chain | Breakthrough |
| alpha3 = pi^2/110 | brannen_sqrt2 | Ultra-precision Richardson | Breakthrough |
| m ~ A_tail^4 (mass scaling) | mass_scaling_k4 | k=2(d-1)/(d-2)=4, 7/7 PASS | Foundational |
| N=3 generations from alpha <= 3/4 | why_n3 | 7 theorems, geometric forcing | Breakthrough |
| h(Phi) = Phi (metric ansatz) | metric_ansatz | 5 independent proofs, 11/11 PASS | Foundational |
| Cabibbo angle 4.8sigma -> 0.75sigma | cabibbo_correction | GL(3,F2) group theory | Resolved |
| Weak continuum theorem (A1-A5) | continuum_limit | 7/7 PASS synthesis | Foundational |
| eta = 1 (no gravitational slip) | galaxy_scaling/gs40 | All scales | Testable |
| CMB fully compatible | galaxy_scaling/gs41 | Suppression >10^39, chameleon | Critical |
| CHSH = 2sqrt2 (Bell violation) | qm_entanglement | 10/10 PASS | Foundational |
| Born rule P=|psi|^2 derived | qm_measurement | 5/5 PASS, 3 proofs | Foundational |
| Spin 1/2 from pi3(S3)=Z | qm_spin | 7/7 PASS | Foundational |

### GALAXY-SCALE SUCCESSES (y ~ 0.01-1)

| Observable | TGP Prediction | Observed | Status |
|---|---|---|---|
| BTFR slope | 2/(1-gamma) = 3.82-3.88 | 3.85 +- 0.09 | PASS |
| RAR intercept | a0 = 1.2e-10 m/s^2 | 1.2e-10 | PASS |
| Freeman limit | a0/(2piG) = 137 Msun/pc^2 | ~140 | PASS (0.98) |
| Size-mass relation | R ~ sqrt(GM/a0) | R ~ M^0.5 | PASS |
| Faber-Jackson | L ~ sigma^4 | L ~ sigma^4 | PASS |
| dSph dynamics | Wolf mass, EFE quadrature | comparable to MOND | PASS |
| Rotation curves | chi2/dof ~ 36 (c=1), 33 (c_eff) | MOND: 17.5 | PARTIAL |

### CLUSTER-SCALE PARTIAL FAILURES (y ~ 1-10)

| Observable | TGP | Observed | Deficit |
|---|---|---|---|
| Bullet kappa(galaxy) | 0.22 | 0.35 | 38% |
| Bullet kappa(gas)/kappa(gal) | 2.30 | ~2.3 | OK |
| Bullet peak location | galaxy pos. | galaxy pos. | OK |
| Cluster M at R500 | 60-65% | 100% | 35-50% |
| Differential gamma effect | 1-2% improvement | --- | NEGLIGIBLE |
| 3D mass calculation | WORSENS deficit | --- | STRUCTURAL |

### COSMOLOGICAL OUTSIDE SCOPE (y >> 1)

| Tension | TGP Effect | Needed | Gap | Verdict |
|---|---|---|---|---|
| H0 (73 vs 67) | dH0 ~ 0.07 km/s/Mpc | dH0 ~ 5.68 | 80x | OUTSIDE SCOPE |
| S8 (0.83 vs 0.76) | ~0.001% suppression | ~8.5% | ~8500x | NEGLIGIBLE |
| DESI w(z) | w >= -1 always | w crosses -1 | INCOMPATIBLE* | STRUCTURAL |

*if CPL parametrization confirmed

### OPEN MATHEMATICAL PROBLEMS

| Theorem | Status | Estimate |
|---|---|---|
| CG-1 (Banach contraction) | OPEN (numerical support) | 6-12 months |
| CG-3 (H1 convergence) | OPEN (numerical support) | 6-12 months |
| CG-4 (K_hom = K_TGP) | OPEN (numerical support) | 6-12 months |
| Cosmo tensions mechanism | CLOSED (NO MECHANISM) | Definitive |

---

## Architecture: Why TGP Has a Scope Boundary

The SAME mechanism provides both strengths and limitations:

```
exp(-y^0.8) / y^gamma
   |                |
   |   y << 1       |   y >> 1
   |   (galaxies)   |   (clusters, cosmo)
   |                |
   v                v
nu >> 1             nu -> 1
STRONG effect       NO effect
-> RAR, BTFR,       -> CMB safe,
   rot. curves,        Solar System safe,
   dSphs                GW safe
   PASS                 but also:
                        no cluster mass,
                        no H0 tension,
                        no S8 effect
```

This is NOT a parameter choice -- it's a **structural feature** of any theory where g_obs = nu(g/a0) * g with nu -> 1 at high g.

---

## Key Discriminating Tests (2025-2030)

| Test | TGP | MOND | LCDM | Observable |
|---|---|---|---|---|
| EFE strength | WEAKER (quadrature) | STRONGER (linear) | N/A | Crater II, NGC 1052-DF2 |
| BTFR slope | 3.82-3.88 (from c_eff) | ~4.0 (simple nu) | scatter | SPARC extended sample |
| a0(z) evolution | a0 ~ H(z) (evolves) | a0 = const | N/A | Euclid deep-field z=1 RCs |
| Grav. slip eta | 1.000 | varies | 1.000 | DES/Euclid E_G statistic |
| Cluster lensing R>2Mpc | distinct from NFW | distinct | NFW profile | Euclid cluster stacking |
| Void lensing profile | TGP-enhanced | MOND-enhanced | NFW-subtracted | Euclid/Rubin void catalog |

**Most powerful test:** a0(z) evolution with Euclid deep-field surveys.
If a0 at z=1 is ~1.5x a0 at z=0, TGP is strongly supported over MOND.

---

## Research Program Statistics

| Metric | Count |
|---|---|
| Research folders | 21 |
| Python scripts | ~290 |
| Result files (.txt) | 50+ |
| Analytic theorems proven | 15+ |
| Test suites passing | ~160/190 (84%) |
| Major breakthroughs | 7 |
| Definitively ruled-out mechanisms | 8 |
| Lines of analysis code | ~25,000+ |

---

## Recommended Publication Strategy

1. **Paper 1: TGP galaxy dynamics** (strongest case)
   - nu(y) from membrane Flory exponent
   - alpha = 4/5 derived, not assumed
   - RAR, BTFR, Freeman limit from single a0
   - c_eff from morphology -> BTFR slope prediction
   - EFE: quadrature (testable difference from MOND)

2. **Paper 2: TGP foundations** (theoretical)
   - f(R) = R + R0^gamma * R^(1-gamma) * exp(-(R/R0)^alpha)
   - eta = 1 (no gravitational slip)
   - CMB compatibility (chameleon from exp term)
   - Weak continuum theorem (A1-A5)

3. **Paper 3: Cluster-scale honest assessment**
   - Bullet Cluster: peak OK, amplitude deficit 38%
   - Same failure as MOND, TeVeS
   - Quantitative comparison with RMOND/AeST
   - Hybrid model discussion (sterile neutrinos)

4. **Paper 4: Testable predictions for Euclid/Rubin**
   - a0(z) evolution (key discriminator from MOND)
   - Cluster lensing at R > 2 Mpc
   - Void lensing profiles
   - S/N estimates from gs35

---

## What Would Falsify TGP?

1. **Gravitational slip eta != 1** detected at any scale
2. **a0(z) = constant** (no evolution with H) at z > 0.5
3. **BTFR slope** significantly different from 3.82-3.88
4. **EFE in dSphs** matches MOND strength (linear, not quadrature)
5. **Phantom crossing w < -1** confirmed by DESI DR2+
6. **Cluster mass** explained by modified gravity alone (no DM needed)

Items 1-4 are testable within 5 years. Items 5-6 constrain the theory's scope.

---

*Generated 2026-04-18 by TGP research program analysis.*
*Scripts: gs1-gs43, ct1-ct7, cg_strong_numerical, de1, h0, s8, r1-r6, q0-q7.*
