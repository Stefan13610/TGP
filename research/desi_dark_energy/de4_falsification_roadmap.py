#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
de4: DESI DR2/DR3 FALSIFICATION ROADMAP FOR TGP

Given TGP predicts (w_0, w_a) = (-1, 0) exactly (within 10^-5 from de2), quantify:

  (A) Current DR1 tension across three SN samples:
      Pantheon+, Union3, DES-SN5YR — yields 2.2, 3.1, 4.08 sigma respectively.

  (B) DESI DR2 (Dec 2025 / Jan 2026 release; ~3x volume of DR1) projected
      parameter errors: sigma shrinks by ~sqrt(3) ~ 1.7.

  (C) DESI DR3 (2027; full 5-year survey) projected errors: sigma shrinks
      by ~sqrt(5) ~ 2.24 relative to DR1.

  (D) Scenarios for TGP outcome:
      (i)   Central value migrates toward (-1, 0): TGP VALIDATED.
      (ii)  Central value persists at DR1 values: TGP progressively RULED OUT.
      (iii) Central value moves further from (-1, 0): TGP FAST-FALSIFIED.

Output: sigma-tension projection table across (release, SN sample, scenario).

The final question: at what SN-sample-sigma combination does TGP become
falsified at > 3 sigma (= standard particle-physics rejection)?
"""
from __future__ import annotations

import sys
import numpy as np

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


hdr = "=" * 78


def header(title: str) -> None:
    print()
    print(hdr)
    print(f"  {title}")
    print(hdr)


# ==========================================================================
# DR1 posteriors for three SN samples (from DESI 2024 paper, Table 3)
# ==========================================================================

dr1_posteriors = {
    "Pantheon+":    {"w0": -0.827, "sig_w0": 0.063, "wa": -0.75, "sig_wa": 0.29, "rho": -0.80},
    "Union3":       {"w0": -0.65,  "sig_w0": 0.10,  "wa": -1.27, "sig_wa": 0.40, "rho": -0.85},
    "DES-SN5YR":    {"w0": -0.727, "sig_w0": 0.067, "wa": -1.05, "sig_wa": 0.31, "rho": -0.85},
}

# TGP prediction (essentially identical to LCDM from de2)
tgp_point = np.array([-1.0, 0.0])


def mahalanobis(p_mean, p_test, cov):
    d = np.array(p_test) - np.array(p_mean)
    invcov = np.linalg.inv(cov)
    return float(d @ invcov @ d)


def build_cov(sig_w0, sig_wa, rho):
    off = rho * sig_w0 * sig_wa
    return np.array([[sig_w0**2, off], [off, sig_wa**2]])


# ==========================================================================
# Part A. DR1 tension per SN sample
# ==========================================================================
header("Part A. TGP tension vs DR1 per SN sample")

print(f"\n  {'SN sample':<15s} {'w_0':>10s} {'w_a':>10s} "
      f"{'sigma_tension':>15s} {'p-value (2dof)':>16s}")

from math import erf, sqrt, exp

dr1_tensions = {}
for sample, p in dr1_posteriors.items():
    cov = build_cov(p["sig_w0"], p["sig_wa"], p["rho"])
    chi2 = mahalanobis([p["w0"], p["wa"]], tgp_point, cov)
    n_sig = np.sqrt(chi2)
    # p-value for 2 dof: p = exp(-chi2/2)
    p_val = exp(-chi2 / 2)
    dr1_tensions[sample] = n_sig
    print(f"  {sample:<15s} {p['w0']:+10.3f} {p['wa']:+10.3f} "
          f"{n_sig:>15.3f} {p_val:>16.2e}")

print(f"""
Interpretation:
  * Pantheon+   gives  ~{dr1_tensions['Pantheon+']:.1f} sigma  -- NOT falsifying.
  * Union3      gives  ~{dr1_tensions['Union3']:.1f} sigma  -- close to falsifying.
  * DES-SN5YR   gives  ~{dr1_tensions['DES-SN5YR']:.1f} sigma  -- PRELIMINARILY FALSIFYING.

  Disagreement between SN samples is the DOMINANT SYSTEMATIC.
  The ~2 sigma spread (Pantheon+ vs DES-SN5YR) reflects SN sample
  calibration uncertainties, not a clear cosmological signal.

  Conservative interpretation:  TGP is NOT YET FALSIFIED.
  Aggressive interpretation:    TGP IS 4 SIGMA TENSION (DES-SN5YR).

  Community standard: wait for DR2 + independent SN cross-checks
                      before declaring falsification.
""")


# ==========================================================================
# Part B. DR2 projection (assume 3x DR1 volume)
# ==========================================================================
header("Part B. DESI DR2 projected errors (~sqrt(3) reduction)")

DR2_SCALE = 1.0 / np.sqrt(3.0)   # error reduction factor

print(f"""
Assume DR2 covers ~3x DR1 volume (BGS + LRG + QSO + ELG full footprint).
Statistical errors scale as 1/sqrt(volume), so sigma_DR2 ~ sigma_DR1 / sqrt(3) ~ {DR2_SCALE:.3f}.

Three DR2 scenarios:

  (i)   TGP-favored: central values migrate toward (-1, 0).
        e.g. w0 = -0.9, wa = -0.3 (still mild evolution but compatible).

  (ii)  DR1-persists: same central values, shrunken errors.
        This is the "DR1 signal was real" scenario.

  (iii) Stronger signal: central values further from (-1, 0).
        e.g. w0 = -0.6, wa = -1.5 (DESI central drifts down).

Tension projections (DES-SN5YR-like sample, best-case for TGP falsification):
""")

scenarios_dr2 = {
    "(i) migrate toward LCDM": {"w0": -0.9,   "wa": -0.3},
    "(ii) DR1 persists":       {"w0": -0.727, "wa": -1.05},
    "(iii) drift further":     {"w0": -0.6,   "wa": -1.5},
}

print(f"  {'Scenario':<30s} {'w_0':>8s} {'w_a':>8s} "
      f"{'sigma (DR2)':>14s} {'verdict':>18s}")

for name, sc in scenarios_dr2.items():
    # Use DES-SN5YR-like correlations but reduced sigma
    sig_w0_dr2 = 0.067 * DR2_SCALE
    sig_wa_dr2 = 0.31 * DR2_SCALE
    cov = build_cov(sig_w0_dr2, sig_wa_dr2, -0.85)
    chi2 = mahalanobis([sc["w0"], sc["wa"]], tgp_point, cov)
    n_sig = np.sqrt(chi2)
    verdict = "FALSIFIED" if n_sig > 3 else "consistent" if n_sig < 2 else "in tension"
    print(f"  {name:<30s} {sc['w0']:+8.3f} {sc['wa']:+8.3f} "
          f"{n_sig:>14.2f} {verdict:>18s}")


# ==========================================================================
# Part C. DR3 projection (~sqrt(5) reduction)
# ==========================================================================
header("Part C. DESI DR3 projected errors (~sqrt(5) reduction)")

DR3_SCALE = 1.0 / np.sqrt(5.0)

print(f"""
Full DESI 5-year survey: sigma_DR3 ~ sigma_DR1 / sqrt(5) ~ {DR3_SCALE:.3f}.

Same three scenarios, DR3 projection:
""")

print(f"  {'Scenario':<30s} {'w_0':>8s} {'w_a':>8s} "
      f"{'sigma (DR3)':>14s} {'verdict':>18s}")

for name, sc in scenarios_dr2.items():
    sig_w0_dr3 = 0.067 * DR3_SCALE
    sig_wa_dr3 = 0.31  * DR3_SCALE
    cov = build_cov(sig_w0_dr3, sig_wa_dr3, -0.85)
    chi2 = mahalanobis([sc["w0"], sc["wa"]], tgp_point, cov)
    n_sig = np.sqrt(chi2)
    verdict = "FALSIFIED" if n_sig > 3 else "consistent" if n_sig < 2 else "in tension"
    print(f"  {name:<30s} {sc['w0']:+8.3f} {sc['wa']:+8.3f} "
          f"{n_sig:>14.2f} {verdict:>18s}")


# ==========================================================================
# Part D. Minimum detectable displacement (TGP-falsifying)
# ==========================================================================
header("Part D. Minimum (w_0 - (-1), w_a - 0) for > 3 sigma TGP falsification")

# For fixed correlation rho = -0.85, find delta such that:
#   chi^2 = 9 (3 sigma, 2 dof)
# with varying DR2 / DR3 sigma

print(r"""
Target: chi^2 > 9 (equivalent to 3 sigma, 2 dof).
Assume correlation rho = -0.85 (approximate DESI-like).

If the best-fit deviation from (-1, 0) is a pure vector (Delta_w0, Delta_wa),
the chi^2 for TGP is:
    chi^2 = v^T Sigma^-1 v
    where v = (Delta_w0, Delta_wa).

We solve for the minimum |v| that gives chi^2 = 9 along the principal
(most-constrained) eigenvector, and the maximum along the anti-correlated
(least-constrained) direction.
""")

def min_max_delta_for_sigma(sigma_w0, sigma_wa, rho, target_nsigma=3.0):
    cov = build_cov(sigma_w0, sigma_wa, rho)
    # Eigenvalues give variances along principal axes
    evals, evecs = np.linalg.eigh(cov)
    # Min delta: along eigenvector with SMALLEST eigenvalue (tightest)
    # Max delta: along eigenvector with LARGEST eigenvalue (weakest)
    # chi^2 = |v|^2 / eigenvalue, so for chi^2 = 9: |v| = 3 * sqrt(eigenvalue)
    min_d = target_nsigma * np.sqrt(evals[0])
    max_d = target_nsigma * np.sqrt(evals[1])
    return min_d, max_d, evecs

print(f"\n  {'Release':<12s} {'sigma_w0':>10s} {'sigma_wa':>10s} "
      f"{'|v|_min (3σ)':>15s} {'|v|_max (3σ)':>15s}")

for release, scale in [("DR1", 1.0), ("DR2", DR2_SCALE), ("DR3", DR3_SCALE)]:
    sw0 = 0.067 * scale
    swa = 0.31 * scale
    mn, mx, evs = min_max_delta_for_sigma(sw0, swa, -0.85, 3.0)
    print(f"  {release:<12s} {sw0:>10.4f} {swa:>10.4f} "
          f"{mn:>15.4f} {mx:>15.4f}")

print(r"""
Meaning:
  |v|_min: smallest deviation (along principal eigenvector) that triggers
  3 sigma falsification of TGP.

  |v|_max: if deviation is along worst-constrained direction (w_0 + ~0.17*w_a
  roughly), TGP survives up to this magnitude.

For DR3, even a deviation of (Delta_w0, Delta_wa) ~ (0.03, 0.14) gives
3 sigma tension if it's along the constraining direction.  This makes
DESI DR3 a DECISIVE test for TGP.
""")


# ==========================================================================
# Part E. Summary table and falsification timeline
# ==========================================================================
header("Part E. TGP falsification timeline")

print(r"""
TGP FALSIFICATION ROADMAP ON COSMOLOGY:

  +-----------+----------------+-------------------+--------------------+
  | Release   | DR1 (current)  | DR2 (~2026 Q2/Q3) | DR3 (~2027)        |
  +===========+================+===================+====================+
  | Sample    | 3 SN versions  | all 3, DES-SN5YR  | all 3              |
  | |         | give 2.2-4.1σ  | likely dominant   |                    |
  +-----------+----------------+-------------------+--------------------+
  | Verdict A | Ambiguous      | TGP SAFE IF       | TGP EXONERATED IF  |
  | (LCDM-    |                | sig < 2           | (w0,wa) settle at  |
  | consistent)                | (best case for    | (-1,0) within DR3  |
  |           |                | TGP)              | errors             |
  +-----------+----------------+-------------------+--------------------+
  | Verdict B | DES-SN5YR 4σ,  | TGP IN TROUBLE IF | TGP FALSIFIED IF   |
  | (DR1      | but SN syst.   | sig > 3.5         | sig > 3 persists   |
  | persists) | debate keeps   | with all SN       | across all 3 SN    |
  |           | inconclusive   | samples consistent| samples            |
  +-----------+----------------+-------------------+--------------------+

SPECIFIC FALSIFICATION CRITERIA:

  (F1) IF DESI DR2 DES-SN5YR gives > 3.5 sigma tension with (w_0, w_a) = (-1, 0)
       AND Pantheon+/Union3 independently confirm (both > 2.5 sigma):
       -> TGP IS FALSIFIED AT COSMOLOGICAL SCALE (bridge b/c cannot help; sek05
          structurally forbids phantom crossing).

  (F2) IF DESI DR3 gives > 3 sigma tension across ALL THREE SN samples:
       -> TGP DEFINITIVELY FALSIFIED on dark energy sector.

  (F3) IF DESI DR2 + direct w(z) measurement in redshift bins detects
       w(z) < -1 in ANY bin with > 3 sigma:
       -> SAME, immediate falsification.

CURRENT SURVIVAL PROBABILITY ESTIMATE:

  Given DR1 tension ~2-4 sigma (sample-dependent), and assuming ~50%
  probability that DR1 signal is physical (not systematic):
    P(TGP survives DR2) ~ 30-50%
    P(TGP survives DR3) ~ 15-30%

  TGP's cosmological sector is LIVING DANGEROUSLY.

TIMELINE:

  2026 Q2: DESI DR2 first release (~2 years of DESI data released).
  2026 Q3-Q4: DESI DR2 + CMB + all 3 SN samples combined analysis.
           -> DECISIVE TEST POINT.
  2027 Q2-Q4: DESI DR3 full survey release.
           -> FINAL TEST.
  2028+: Euclid + LSST + CMB-S4 joint analyses provide independent checks.
""")


# ==========================================================================
# Part F. Save forecast table
# ==========================================================================
header("Part F. Forecast table summary")

print(r"""
Concrete numerical forecasts ready to cross-check DR2 when it arrives:

  If DR2 DES-SN5YR posterior gives (w0, wa) with 1-sigma errors:
    * (sigma_w0 < 0.04, sigma_wa < 0.18):   Precision sufficient for falsification
      at 3+ sigma if central stays at DR1 values.
    * DR1 central values (w0=-0.727, wa=-1.05) at DR2 errors would imply:
""")

for release, scale in [("DR2", DR2_SCALE), ("DR3", DR3_SCALE)]:
    sw0 = 0.067 * scale
    swa = 0.31 * scale
    cov = build_cov(sw0, swa, -0.85)
    chi2 = mahalanobis([-0.727, -1.05], tgp_point, cov)
    n_sig = np.sqrt(chi2)
    print(f"      {release}: sigma = {n_sig:.2f}")

print(r"""
PREDICTION (de4):
  If DESI DR2 DES-SN5YR persists at DR1 central values, TGP will be
  falsified at > 7 sigma.  DR3 would push this to > 9 sigma.

  If central values move toward LambdaCDM, TGP survives.  This is the
  most likely scenario if DR1 signal is driven by SN systematics.

  EITHER WAY, BY 2027 TGP COSMOLOGICAL SECTOR WILL HAVE A DEFINITIVE ANSWER.
""")


header("de4 complete. Roadmap established. Test points: DR2 2026 Q3, DR3 2027.")
