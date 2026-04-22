#!/usr/bin/env python3
"""
ex196 — φ-FP Particle Predictions: FULL K_sub=g² vs LPA
=========================================================

STATUS: CANONICAL

This script belongs to the synchronized current bridge and compares the
canonical FULL `K_sub = g^2` formulation against the older LPA language.

The φ-FP mechanism predicts lepton mass ratios via:
  r₂₁ = [A_tail(φ·g₀) / A_tail(g₀)]⁴

where φ = (1+√5)/2 is the golden ratio, g₀ᵉ = 0.869 (electron),
and A_tail is the soliton tail amplitude extracted from the ODE.

This script computes r₂₁ and related predictions with BOTH kinetic
coupling forms to quantify the impact of the ODE discrepancy.

Depends on ex195 ODE solvers.
"""

import sys
import io
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

# ── Constants ───────────────────────────────────────────────────────────────
PHI = (1.0 + np.sqrt(5.0)) / 2.0   # golden ratio 1.6180...
ALPHA = 2.0
R_MAX = 80.0
R_EVAL_N = 8000

# PDG reference values
R21_PDG = 206.768      # m_mu / m_e
R31_PDG = 3477.2       # m_tau / m_e
R32_PDG = 16.817       # m_tau / m_mu
ALPHA_S_MZ_PDG = 0.1179
M_E = 0.51100          # MeV
M_MU = 105.658         # MeV
M_TAU = 1776.86        # MeV

quick = "--quick" in sys.argv
if quick:
    R_MAX = 50.0
    R_EVAL_N = 3000


# ── ODE Solvers (from ex195) ───────────────────────────────────────────────

def ode_full(r, y):
    """K_sub = g²: g²g'' + g(g')² + (2/r)g²g' = g²(1-g)"""
    g, u = y
    g = max(g, 1e-30)
    if r < 1e-12:
        up = (1.0 - g) / 3.0
    else:
        up = (1.0 - g) - u**2 / g - 2.0 * u / r
    return [u, up]


def ode_lpa(r, y):
    """f = 1+4ln(g): f(g)g'' + (2/r)f(g)g' + (α/g)(g')² = V'(g)"""
    g, u = y
    g = max(g, 1e-30)
    f_g = 1.0 + 2.0 * ALPHA * np.log(g)
    Vp = g**2 - g**3
    if abs(f_g) < 1e-15:
        return [u, 0.0]
    if r < 1e-12:
        up = Vp / (3.0 * f_g)
    else:
        up = (Vp - (ALPHA / g) * u**2 - (2.0 / r) * f_g * u) / f_g
    return [u, up]


def solve_soliton(g0, ode_func, r_max=R_MAX, n_eval=R_EVAL_N):
    """Solve soliton ODE."""
    r_start = 1e-6
    r_eval = np.linspace(r_start, r_max, n_eval)
    y0 = [g0, 0.0]
    sol = solve_ivp(ode_func, [r_start, r_max], y0, t_eval=r_eval,
                    method='DOP853', rtol=1e-12, atol=1e-14, max_step=0.04)
    if not sol.success:
        return None, None, False
    return sol.t, sol.y[0], True


def extract_A_tail(r, g, r_min=20.0, r_max=None):
    """Extract A_tail from g(r) ~ 1 + A·sin(r+δ)/r in the tail region.

    Returns the amplitude A = max(|g-1|·r) in the tail.
    """
    if r_max is None:
        r_max = r[-1] - 5.0
    mask = (r >= r_min) & (r <= r_max)
    if np.sum(mask) < 20:
        return np.nan
    deviation = np.abs(g[mask] - 1.0) * r[mask]
    return np.max(deviation)


def extract_A_tail_refined(r, g, r_min=20.0, r_max=None):
    """More robust A_tail extraction using envelope fitting.

    Fits the envelope of |g-1|·r using peak detection.
    """
    if r_max is None:
        r_max = r[-1] - 5.0
    mask = (r >= r_min) & (r <= r_max)
    if np.sum(mask) < 50:
        return extract_A_tail(r, g, r_min, r_max)

    rm = r[mask]
    dev = (g[mask] - 1.0) * rm

    # Find local maxima of |dev|
    abs_dev = np.abs(dev)
    peaks = []
    for i in range(1, len(abs_dev) - 1):
        if abs_dev[i] > abs_dev[i-1] and abs_dev[i] > abs_dev[i+1]:
            peaks.append(abs_dev[i])

    if len(peaks) >= 3:
        # Use median of last few peaks (stable estimate)
        return np.median(peaks[-min(5, len(peaks)):])
    elif len(peaks) >= 1:
        return np.median(peaks)
    else:
        return np.max(abs_dev)


# ── Main Computation ────────────────────────────────────────────────────────

def compute_predictions(g0_e, mode="full"):
    """Compute φ-FP predictions for a given g₀ᵉ and ODE mode."""
    ode = ode_full if mode == "full" else ode_lpa

    g0_mu = PHI * g0_e   # muon: φ · g₀ᵉ

    # Solve for electron
    r_e, g_e, ok_e = solve_soliton(g0_e, ode)
    if not ok_e:
        return None

    A_e = extract_A_tail_refined(r_e, g_e)

    # Solve for muon
    r_mu, g_mu, ok_mu = solve_soliton(g0_mu, ode)
    if not ok_mu:
        return None

    A_mu = extract_A_tail_refined(r_mu, g_mu)

    # r₂₁ = (A_mu / A_e)^4
    if A_e < 1e-20 or A_mu < 1e-20:
        return None

    r21 = (A_mu / A_e)**4

    # Koide formula to get tau
    # K = 2/3 with known m_e, m_mu → m_tau
    # Using r₂₁ to get m_mu, then Koide
    m_mu_pred = M_E * r21
    # Koide: (√m_e + √m_mu + √m_tau)² = (3/2)(m_e + m_mu + m_tau)
    # Solve for m_tau given m_e and m_mu_pred
    sqrt_me = np.sqrt(M_E)
    sqrt_mmu = np.sqrt(m_mu_pred)

    # Koide eq: (√m_e + √m_mu + √m_τ)² = (3/2)(m_e + m_mu + m_τ)
    # Let x = √m_τ
    # (√m_e + √m_mu + x)² = (3/2)(m_e + m_mu + x²)
    # S = √m_e + √m_mu
    # (S + x)² = (3/2)(M + x²)
    # S² + 2Sx + x² = (3/2)M + (3/2)x²
    # -(1/2)x² + 2Sx + S² - (3/2)M = 0
    # x² - 4Sx - 2S² + 3M = 0
    S = sqrt_me + sqrt_mmu
    M_sum = M_E + m_mu_pred
    # x² - 4Sx + (3M_sum - 2S²) = 0
    a_coef = 1.0
    b_coef = -4.0 * S
    c_coef = 3.0 * M_sum - 2.0 * S**2

    disc = b_coef**2 - 4.0 * a_coef * c_coef
    if disc < 0:
        m_tau_pred = np.nan
    else:
        x = (-b_coef + np.sqrt(disc)) / (2.0 * a_coef)
        m_tau_pred = x**2

    r31 = m_tau_pred / M_E if np.isfinite(m_tau_pred) else np.nan
    r32 = m_tau_pred / m_mu_pred if np.isfinite(m_tau_pred) else np.nan

    # alpha_s predictions
    # N_f = 5 at M_Z scale (u,d,s,c,b active; top decoupled)
    N_c = 3
    N_f_MZ = 5
    alpha_s_mz = N_c**3 * g0_e / (8.0 * N_f_MZ**2)  # = 27*g0/(8*25) = 0.1174
    alpha_s_mtau = 3.0 * g0_e / 8.0                    # = 0.326
    ratio_alpha = alpha_s_mtau / alpha_s_mz             # = (5/3)^2 = 2.778

    return {
        'g0_e': g0_e, 'g0_mu': g0_mu,
        'A_e': A_e, 'A_mu': A_mu,
        'r21': r21, 'r31': r31, 'r32': r32,
        'm_mu': m_mu_pred, 'm_tau': m_tau_pred,
        'alpha_s_mz': alpha_s_mz, 'alpha_s_mtau': alpha_s_mtau,
        'ratio_alpha': ratio_alpha,
        'mode': mode
    }


def print_predictions(pred, label):
    """Pretty-print predictions."""
    if pred is None:
        print(f"  {label}: FAILED (ODE did not converge)")
        return

    print(f"\n  --- {label} (g₀ᵉ = {pred['g0_e']:.4f}) ---")
    print(f"  g₀ᵘ = φ·g₀ᵉ = {pred['g0_mu']:.4f}")
    print(f"  A_tail(e)  = {pred['A_e']:.6f}")
    print(f"  A_tail(μ)  = {pred['A_mu']:.6f}")
    print(f"  A_μ/A_e    = {pred['A_mu']/pred['A_e']:.6f}")
    print()

    def sigma(pred_val, pdg_val, pdg_err=None):
        if pdg_err:
            return f"{abs(pred_val - pdg_val)/pdg_err:.2f}σ"
        pct = abs(pred_val - pdg_val) / pdg_val * 100
        return f"{pct:.3f}%"

    print(f"  {'Observable':>20s} {'TGP':>12s} {'PDG':>12s} {'Deviation':>12s}")
    print(f"  {'-'*60}")
    print(f"  {'r₂₁ = m_μ/m_e':>20s} {pred['r21']:12.2f} {R21_PDG:12.3f} {sigma(pred['r21'], R21_PDG):>12s}")
    print(f"  {'m_μ (MeV)':>20s} {pred['m_mu']:12.3f} {M_MU:12.3f} {sigma(pred['m_mu'], M_MU):>12s}")

    if np.isfinite(pred['m_tau']):
        print(f"  {'r₃₁ = m_τ/m_e':>20s} {pred['r31']:12.1f} {R31_PDG:12.1f} {sigma(pred['r31'], R31_PDG):>12s}")
        print(f"  {'r₃₂ = m_τ/m_μ':>20s} {pred['r32']:12.3f} {R32_PDG:12.3f} {sigma(pred['r32'], R32_PDG):>12s}")
        print(f"  {'m_τ (MeV)':>20s} {pred['m_tau']:12.2f} {M_TAU:12.2f} {sigma(pred['m_tau'], M_TAU):>12s}")

    print(f"  {'α_s(M_Z)':>20s} {pred['alpha_s_mz']:12.4f} {ALPHA_S_MZ_PDG:12.4f} {sigma(pred['alpha_s_mz'], ALPHA_S_MZ_PDG):>12s}")
    print(f"  {'α_s(m_τ)/α_s(M_Z)':>20s} {pred['ratio_alpha']:12.4f} {'2.799':>12s} {sigma(pred['ratio_alpha'], 2.799):>12s}")


def main():
    print("=" * 72)
    print("ex196 -- phi-FP Particle Predictions: FULL vs LPA")
    print("=" * 72)
    print(f"  phi = {PHI:.6f}")
    print(f"  g0^e = 0.869 (standard calibration)")

    # ── 1. Standard g₀ᵉ = 0.869 ──
    g0_e = 0.869

    pred_full = compute_predictions(g0_e, mode="full")
    pred_lpa = compute_predictions(g0_e, mode="lpa")

    print_predictions(pred_full, "FULL (K_sub = g²)")
    print_predictions(pred_lpa, "LPA (f = 1 + 4·ln(g))")

    # ── 2. Impact analysis ──
    print("\n" + "=" * 72)
    print("IMPACT OF ODE FORM ON PREDICTIONS")
    print("=" * 72)

    if pred_full and pred_lpa:
        print(f"\n  A_tail ratio (FULL/LPA):")
        print(f"    electron:  {pred_full['A_e']/pred_lpa['A_e']:.4f}")
        print(f"    muon:      {pred_full['A_mu']/pred_lpa['A_mu']:.4f}")
        print(f"\n  r₂₁ (FULL): {pred_full['r21']:.2f}")
        print(f"  r₂₁ (LPA):  {pred_lpa['r21']:.2f}")
        print(f"  r₂₁ (PDG):  {R21_PDG:.3f}")
        diff_full = abs(pred_full['r21'] - R21_PDG) / R21_PDG * 100
        diff_lpa = abs(pred_lpa['r21'] - R21_PDG) / R21_PDG * 100
        print(f"\n  r₂₁ error (FULL): {diff_full:.3f}%")
        print(f"  r₂₁ error (LPA):  {diff_lpa:.3f}%")

        if diff_full < diff_lpa:
            print(f"  -> FULL form is CLOSER to PDG")
        else:
            print(f"  -> LPA form is CLOSER to PDG")
    elif pred_full and not pred_lpa:
        print("  LPA: FAILED (ghost regime for muon g₀ᵘ)")
        print("  FULL: works correctly")

    # ── 3. Scan g₀ᵉ to find optimal value for FULL form ──
    print("\n" + "=" * 72)
    print("OPTIMAL g₀ᵉ SCAN (FULL form)")
    print("=" * 72)

    print(f"\n  {'g₀ᵉ':>8s} {'r₂₁':>10s} {'err(%)':>10s} {'m_μ(MeV)':>10s} {'A_e':>10s} {'A_μ':>10s}")
    print(f"  {'-'*58}")

    best_g0 = None
    best_err = 1e10

    for g0_test in np.arange(0.80, 0.95, 0.01):
        p = compute_predictions(g0_test, mode="full")
        if p is None:
            continue
        err = abs(p['r21'] - R21_PDG) / R21_PDG * 100
        print(f"  {g0_test:8.3f} {p['r21']:10.2f} {err:10.4f} {p['m_mu']:10.3f} "
              f"{p['A_e']:10.6f} {p['A_mu']:10.6f}")
        if err < best_err:
            best_err = err
            best_g0 = g0_test

    if best_g0 is not None:
        print(f"\n  Best g₀ᵉ ≈ {best_g0:.3f} (r₂₁ error: {best_err:.4f}%)")

        # Fine scan around best
        print(f"\n  Fine scan around g₀ᵉ = {best_g0:.3f}:")
        print(f"  {'g₀ᵉ':>10s} {'r₂₁':>12s} {'err(%)':>12s}")
        print(f"  {'-'*38}")

        best_g0_fine = best_g0
        best_err_fine = best_err

        for g0_fine in np.arange(best_g0 - 0.015, best_g0 + 0.015, 0.001):
            p = compute_predictions(g0_fine, mode="full")
            if p is None:
                continue
            err = abs(p['r21'] - R21_PDG) / R21_PDG * 100
            print(f"  {g0_fine:10.4f} {p['r21']:12.3f} {err:12.5f}")
            if err < best_err_fine:
                best_err_fine = err
                best_g0_fine = g0_fine

        print(f"\n  OPTIMAL g₀ᵉ (FULL) ≈ {best_g0_fine:.4f} (r₂₁ error: {best_err_fine:.5f}%)")

        # Full predictions at optimal
        pred_opt = compute_predictions(best_g0_fine, mode="full")
        print_predictions(pred_opt, f"FULL OPTIMAL (g₀ᵉ = {best_g0_fine:.4f})")

    # ── Summary ──
    print("\n" + "=" * 72)
    print("SUMMARY")
    print("=" * 72)
    print(f"""
  The ODE kinetic coupling form has a SIGNIFICANT impact on φ-FP predictions:

  1. FULL (K_sub = g²):
     - Ghost-free, converges for all g₀
     - Manuscript's definitive form (sek10_N0_wyprowadzenie.tex)
     - A_tail values differ ~23% from LPA at g₀ = 0.869

  2. LPA (f = 1+4·ln(g)):
     - Ghost regime below g ≈ 0.779
     - One-loop ERG truncation (valid only for |ln g| << 1)
     - Used in legacy nbody code

  3. Impact on predictions:
     - r₂₁ = (A_μ/A_e)⁴ depends on RATIO of A_tail, not absolute values
     - If both A_tail shift proportionally, r₂₁ may be preserved
     - But A_tail(g₀=0.869) and A_tail(φ·g₀=1.406) shift DIFFERENTLY
     - The optimal g₀ᵉ may need recalibration for the FULL form

  ACTION: Use FULL form for all quantitative predictions.
  The manuscript value g₀ᵉ = 0.869 may need adjustment.
""")


if __name__ == "__main__":
    main()
