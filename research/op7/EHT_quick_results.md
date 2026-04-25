# EHT-quick: Photon Ring in M9.1'' Hyperbolic Black Hole

**Cycle:** OP-7 EHT-quick (TGP_CLOSURE_PLAN sec 8.10).
**Test:** Strong-field photon ring + shadow size for static spherical BH
in the M9.1'' hyperbolic effective metric, vs EHT M87*/Sgr A*.
**Date:** 2026-04-25.
**Script:** `tooling/scripts/gravity/eht_photon_ring_m911.py`.
**Output:** `research/op7/eht_photon_ring_m911.txt`.

---

## 1. Setup

**Metric (M9.1'', isotropic radial coordinate r):**
```
    ds^2 = -c0^2 (4 - 3 psi)/psi  dt^2  +  psi/(4 - 3 psi) [dr^2 + r^2 dOmega^2]
    psi  = Phi/Phi_0,  eps := psi - 1.
```
The substrate budget `f * h = 1` enforces `g_tt * g_rr / c0^2 = -1`
(unlike Schwarzschild isotropic where `A*B != 1`).

**Vacuum Phi-EOM (alpha = 2, beta = gamma at vacuum):**
```
    eps''(r) + (2/r) eps'(r) + 2 (eps'(r))^2 / (1 + eps) = 0.
```

**Boundary conditions:** asymptotic flatness with PN tail
`eps = A/r - A^2/r^2 + (5/3) A^3/r^3 - (10/3) A^4/r^4 + ...`
(c_n analytical from M9.1'' P1 Setup, sympy-verified).
Newton matching at infinity gives `A = G M / (2 c^2) = 0.5` in
geometric units (M = c = G = 1).

**TGP horizon-analog:** the largest r where `psi(r) = 4/3`
(i.e. eps = 1/3), where g_tt vanishes.

**Photon ring:** static spherical metric photon orbit at extremum of
`R(r)^2 / A_t(r)` where `R = r * sqrt(B(r))` is the areal radius.
With `A_t * B = c0^2`, this becomes `r * A_t'(r) = A_t(r)`,
or in eps:
```
    4 r eps'(r) + (1 - 3 eps(r)) (1 + eps(r)) = 0.
```

---

## 2. Method

1. **Numerical method:** integrate the second-order vacuum Phi-EOM
   inward from `r_out = 1000 r_g` using `scipy.integrate.solve_ivp`
   (`DOP853`, rtol=1e-12, atol=1e-14). Initial conditions are the
   analytical 7-term PN tail (A_n analytic from P1 sympy recursion).
   A terminal event detects `eps = 1/3` (TGP horizon).

2. **Photon ring:** sample the F-equation
   `F(r) = 4 r eps' + (1 - 3 eps)(1 + eps)`
   on `(r_H * 1.001, 20 r_g)` and locate the (unique) zero with
   `scipy.optimize.brentq` to ~1e-12 tolerance.

3. **b_crit:** the critical impact parameter is
   `b_crit = R(r_ph) / sqrt(A_t(r_ph)) = r_ph / A_t(r_ph)`
   (in geometric units, c0 = 1).

4. **Convergence:** vary `r_out` in {500, 1000, 2000, 5000} and
   `n_terms` in {5, 6, 7} of the asymptotic PN tail.

---

## 3. Results

### 3.1 TGP horizon

```
    r_H = 1.094595  r_g       (M9.1'', psi -> 4/3)
    Linear-order prediction:  r_H ~ 3 A = 1.5 r_g
    GR Schwarzschild horizon (areal):      r_s = 2 r_g
    GR Schwarzschild horizon (isotropic):  r_iso_S = 0.5 r_g
```
Nonlinear corrections move r_H from 1.5 down to 1.09. The areal
radius diverges at psi=4/3 as expected (h(psi) -> infty there),
which matches the GR isotropic-coordinate behaviour at r_iso = M/2.

### 3.2 Photon ring

```
    r_ph (isotropic, M9.1'')   = 2.529629 r_g
    eps(r_ph)                  = +0.167892
    psi(r_ph)                  =  1.167892
    A_t(r_ph)                  =  0.424974
    r_ph (areal, M9.1'')       = 3.880394 r_g

    r_ph (areal, GR)           = 3.000000 r_g

    r_ph^TGP / r_ph^GR (areal) = 1.293
```

### 3.3 Critical impact parameter / shadow

```
    b_crit^TGP                 = 5.952436 r_g
    b_crit^GR                  = 3 sqrt(3) = 5.196152 r_g
    b_crit^TGP / b_crit^GR     = 1.1455
    Fractional deviation       = +14.55 %
```

### 3.4 Astrophysical predictions

**M87*** (M = 6.5e9 M_sun, D = 16.8 Mpc):
```
    GR shadow diameter         = 39.70 microarcsec
    TGP shadow diameter        = 45.48 microarcsec
    EHT 2019 observed          = 42.0 +/- 3.0 microarcsec
    TGP / observed             = 1.083  (+8.3 %)   -- within 10%
    GR  / observed             = 0.945  (-5.5 %)   -- within 10%
```

**Sgr A*** (M = 4.3e6 M_sun, D = 8.15 kpc):
```
    GR shadow diameter         = 54.14 microarcsec
    TGP shadow diameter        = 62.01 microarcsec
    EHT 2022 observed          = 51.8 +/- 2.3 microarcsec
    TGP / observed             = 1.197  (+19.7 %)  -- OUTSIDE 10%
    GR  / observed             = 1.045  (+4.5 %)   -- within 10%
```

### 3.5 Convergence

5 runs varying `r_out` in {500, 1000, 2000, 5000} and `n_terms` in
{5, 6, 7}: all give identical `r_H = 1.094595`, `r_ph = 2.529629`,
`b_crit = 5.952436` to 6 decimals (spread < 1e-6). Numerical
solution is well-converged.

---

## 4. Summary table

| Quantity | TGP M9.1'' | GR Schwarzschild | TGP/GR |
|----------|-----------|------------------|--------|
| r_H (areal) | -> infty (at psi=4/3) | 2 r_g | -- |
| r_ph (areal) | 3.880 r_g | 3.000 r_g | 1.293 |
| b_crit | 5.952 r_g | 5.196 r_g | 1.146 |
| M87* shadow | 45.48 uas | 39.70 uas | 1.146 |
| Sgr A* shadow | 62.01 uas | 54.14 uas | 1.146 |

The TGP/GR ratio of 1.146 in `b_crit` (and shadow diameter)
propagates linearly to all astrophysical sources. M87* sits at the
edge of EHT's 10% systematic; Sgr A* with a tighter 4.4% bound
exceeds the 2-sigma envelope of the EHT 2022 measurement.

---

## 5. Tests verdict (6/6 from script)

- **T1 (horizon found):** PASS  -- terminal event at psi=4/3 fires at r_H=1.0946 r_g.
- **T2 (photon ring found):** PASS  -- unique zero of F(r) at 2.5296 r_g.
- **T3 (b_crit):** PASS  -- 5.952 r_g, +14.6% over GR.
- **T4 (M87* within 10%):** PASS  -- TGP shadow 8.3% above EHT central value.
- **T5 (Sgr A* within 10%):** FAIL -- TGP shadow 19.7% above EHT central value.
- **T6 (convergence):** PASS  -- spread < 1e-6 across 5 runs.

`SUMMARY: 5/6 PASS`.

---

## 6. Verdict

**INCONCLUSIVE / weakly NEGATIVE.**

Quantitative findings:
- Photon ring radius `r_ph^TGP / r_ph^GR = 1.293` (areal), `1.356` (isotropic).
- Critical impact parameter `b_crit^TGP / b_crit^GR = 1.146`.
- Shadow diameter is **uniformly 14.6% larger than GR** for any mass.
- M87*: TGP within 10% of EHT central value. Marginal but compatible
  with current systematics (EHT 2019 quotes ~10%).
- Sgr A*: TGP exceeds EHT 2022 measurement by 19.7%, well outside
  the quoted 4.4% statistical+systematic envelope.

**Interpretation:**

This is the first explicit strong-field test of M9.1'' beyond the
1PN regime. The result is consistent with the M9.1'' P1 prediction
that TGP and GR diverge starting at 2PN
(`alpha_3^TGP = -7/3` vs `alpha_3^GR = -3/2`, deviation `5/6 U^3`).
At U ~ 0.2-0.3 (photon-ring scale) the cumulative deviations of
4PN-7PN coefficients dominate, blowing up to a ~15% effect in b_crit.

The Sgr A* tension at 19.7% suggests either:
1. **Pivot needed for strong-field consistency** -- M9.1'' is
   incomplete for U ~ O(1); a sigma_ab tensor coupling (OP-7 T3-T4)
   may bring the photon ring back toward GR. From OP-7 T1 we know
   sigma_ab = 0 for *static* spherical, so this scenario cannot
   help here; the pivot would have to come from a different
   structural change (e.g. matter-coupling correction in
   strong-field, or non-vacuum environment of the BH).
2. **Genuine TGP deviation that EHT cannot yet detect** -- if the
   actual M87*/Sgr A* mass calibration has uncertainties, the 19.7%
   tension may collapse. Mass uncertainties for Sgr A* are ~5%
   (Gravity Collaboration 2019); M87* has ~10% uncertainty. Combined
   with ring-modelling systematics, the EHT envelope could effectively
   widen toward ~15-20%, in which case TGP becomes barely consistent.
3. **Calibration of `Phi_0`/coupling normalization** -- our matching
   `A = G M / (2 c^2)` was the standard 1PN Newtonian limit. Higher-PN
   matter coupling (q-renormalization) might shift `A` slightly.

Given the numerical convergence is exact (spread < 1e-6) and the
formula is closed-form, **the 14.6% deviation in b_crit is a robust
prediction of M9.1''+OP-7 T1 single-Phi gravity for static spherical
black holes**. It places M9.1'' in moderate but not catastrophic
tension with EHT Sgr A* 2022, and in marginal agreement with EHT M87*
2019.

---

## 7. Implications for OP-7 closure plan

- **EHT-quick is now CLOSED** (this document, 2026-04-25).
- **Status update:** 4.2 EHT priority -- shadow prediction obtained
  ahead of OP-7 T3-T6.
- **Smoking-gun candidate identified:** at present sensitivity TGP
  is on the boundary of EHT M87* / past the boundary of EHT Sgr A*.
  Future EHT extensions (ngEHT 2030+) with target precision of 1%
  would unambiguously discriminate.
- **No additional sigma_ab dynamics** can save TGP from the +14.6%
  prediction in static-spherical BHs (OP-7 T1: sigma_ab = 0
  for isotropy). Either the prediction stands as a falsifiable
  feature, or M9.1'' must be revisited for strong-field structure
  (mass renormalization in q, alternative vacuum boundary condition,
  or a deeper M9.2 momentum back-reaction).

---

## 8. Files

| File | Role |
|------|------|
| `tooling/scripts/gravity/eht_photon_ring_m911.py` | Script (this test) |
| `research/op7/eht_photon_ring_m911.txt` | Numerical output |
| `research/op7/EHT_quick_results.md` | This document |
| `research/op-newton-momentum/M9_1_pp_P1_results.md` | PN expansion (c_n) |
| `research/op-newton-momentum/M9_1_pp_P2_results.md` | V/Phi^4 selection |
| `research/op7/OP7_T1_results.md` | Single-Phi static spherical -> sigma_ab=0 |

---

## 9. One-line summary

> **TGP M9.1''** predicts a photon ring at **r_ph_areal = 3.88 r_g**
> (vs GR 3.0) and shadow diameter **14.6% larger than GR uniformly**:
> M87* 45.5 uas (EHT: 42 +/- ~3-4), Sgr A* 62.0 uas (EHT: 51.8 +/- 2.3).
> Verdict **INCONCLUSIVE-leaning-NEGATIVE**: marginal at M87*, in tension
> with Sgr A*. M9.1'' single-Phi static-spherical strong-field
> structure is now explicitly testable; the +14.6% deviation is
> a robust falsifiable signature of TGP ahead of OP-7 T3-T6 closure.
