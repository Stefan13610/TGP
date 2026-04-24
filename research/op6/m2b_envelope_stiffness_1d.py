#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
OP-6 / M2-b numerical test — envelope stiffness K_eff(Phi) of H_1.
===================================================================

See M2b_numerical_plan.md for the full derivation. This script is a
numerical surrogate for the envelope FRG of H_1 = the canonical core-paper
substrate Hamiltonian. We simulate H_1 in 1D with continuous fields,
measure the connected 2-pt function of the Z_2-even density
Phi(x) = s(x)^2, and extract the effective kinetic stiffness

    K_eff(<Phi>) = xi_Phi^2 / chi_Phi               ... (Ornstein-Zernike)

with

    chi_Phi = sum_x C_Phi(x),          C_Phi(x) = <Phi(0)Phi(x)> - <Phi>^2
    xi_Phi^2 = sum_x x^2 C_Phi(x) / (2 * chi_Phi)    (d = 1 denominator)

Scanning control parameters (beta, m0^2, lambda_0) varies <Phi>. A
power-law fit K_eff = C * <Phi>^p answers the M2-b question:

    p  +1  <=>  K(Phi) ~ Phi  <=>  K(phi) ~ phi^4  <=>  alpha = 2 (TGP target)
    p   0  <=>  K(Phi) = const   (alpha = 1/2 in Phi)
    p  -1  <=>  K(Phi) ~ 1/Phi   (Gaussian sub-limit, alpha = -1/2)

Hamiltonian (per the core paper, H_1; NOT H_3):

    H_1 = sum_i [ (m0^2/2) s_i^2 + (lambda_0/4) s_i^4 ]
          - J sum_<ij> s_i s_j

Part E (added 2026-04-24) also simulates H_3 bilinear (the object
inventoried as "H_3" in M1 §1),

    H_3_bil = sum_i [(m0^2/2) s_i^2 + (lambda_0/4) s_i^4]
              - J sum_<ij> (s_i s_j)^2

and checks the M2c claim that this gives K(Phi) = const (p_Phi = 0),
distinguishing it from the Ginzburg-Landau functional H_GL (which
gives K(phi) = K_geo * phi^4 ⇔ K(Phi) ∝ Phi ⇔ p_Phi = +1).

Dependencies: numpy only. (Does not require scipy / matplotlib.)
"""

import sys, time, warnings
warnings.filterwarnings("ignore")
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import numpy as np


# -----------------------------------------------------------------------------
# Monte Carlo: checkerboard Metropolis for H_1 in 1D with continuous s.
# -----------------------------------------------------------------------------

def mc_sweep_checkerboard(s, beta, m2, lam, J, sigma_prop, rng):
    """One full Metropolis sweep via bipartite (even/odd) updates.

    H_1 bond term: -J Σ_<ij> s_i s_j  (bilinear in s, Z_2-odd).
    """
    N = s.size
    for parity in (0, 1):
        idx = np.arange(parity, N, 2)
        s_old = s[idx]
        s_new = s_old + rng.normal(0.0, sigma_prop, size=idx.size)
        left = s[(idx - 1) % N]
        right = s[(idx + 1) % N]
        dE = ((m2 / 2.0) * (s_new ** 2 - s_old ** 2)
              + (lam / 4.0) * (s_new ** 4 - s_old ** 4)
              - J * (s_new - s_old) * (left + right))
        # Metropolis
        u = rng.random(idx.size)
        accept = (dE <= 0.0) | (u < np.exp(-beta * np.clip(dE, None, 50.0)))
        s[idx[accept]] = s_new[accept]
    return s


def mc_sweep_checkerboard_h3bilinear(s, beta, m2, lam, J, sigma_prop, rng):
    """One sweep for H_3 bilinear (quartic-looking, but Ising-in-Φ bond).

    H_3 bilinear bond: -J Σ_<ij> (s_i s_j)^2 = -J Σ_<ij> Φ_i Φ_j, Φ = s^2.

    Energy change per site (updating s → s'):
      (m^2/2)(s'^2 - s^2) + (λ/4)(s'^4 - s^4) - J (s'^2 - s^2)(left^2 + right^2).
    """
    N = s.size
    for parity in (0, 1):
        idx = np.arange(parity, N, 2)
        s_old = s[idx]
        s_new = s_old + rng.normal(0.0, sigma_prop, size=idx.size)
        left = s[(idx - 1) % N]
        right = s[(idx + 1) % N]
        d_s2 = s_new ** 2 - s_old ** 2
        d_s4 = s_new ** 4 - s_old ** 4
        nn_Phi = left ** 2 + right ** 2
        dE = ((m2 / 2.0) * d_s2
              + (lam / 4.0) * d_s4
              - J * d_s2 * nn_Phi)
        u = rng.random(idx.size)
        accept = (dE <= 0.0) | (u < np.exp(-beta * np.clip(dE, None, 50.0)))
        s[idx[accept]] = s_new[accept]
    return s


def autotune_sigma(s, beta, m2, lam, J, rng,
                   target_acc=0.4, n_tune=400, sigma0=0.5):
    """Adapt Metropolis step size to hit target acceptance ~40%."""
    sigma = sigma0
    for _ in range(n_tune):
        s_before = s.copy()
        mc_sweep_checkerboard(s, beta, m2, lam, J, sigma, rng)
        # Rough acceptance proxy: fraction of sites that moved
        moved = np.mean(np.abs(s - s_before) > 1e-12)
        # Logistic-like update
        if moved > target_acc + 0.05:
            sigma *= 1.1
        elif moved < target_acc - 0.05:
            sigma *= 0.9
    return sigma


# -----------------------------------------------------------------------------
# Observable: Phi correlator via FFT, then chi, xi, K_eff.
# -----------------------------------------------------------------------------

def phi_correlator_stats(Phi_trajectory, n_k_fit=4, n_jk=8):
    """
    Phi_trajectory: array (n_samples, N) of Phi = s^2 per site per sample.
    Extract chi_Phi and K_eff by an Ornstein-Zernike fit in k-space on the
    lowest few momentum modes, with jackknife error estimation:

        1 / S(k_hat)  ≈  1/chi_Phi  +  K_eff * k_hat^2,   k_hat^2 = 2(1 - cos k).

    Returns: mean_Phi, chi_Phi, xi2_Phi, K_eff, K_eff_err, Sk.
    """
    n_samples, N = Phi_trajectory.shape
    Phi_mean_all = Phi_trajectory.mean()

    def _oz_from_samples(traj):
        Phi_shifted = traj - Phi_mean_all
        F = np.fft.rfft(Phi_shifted, axis=1)
        Sk_ = (F * np.conj(F)).real.mean(axis=0) / N
        j_ = np.arange(N // 2 + 1)
        k_ = 2.0 * np.pi * j_ / N
        k_hat2_ = 2.0 * (1.0 - np.cos(k_))
        n_fit = min(n_k_fit + 1, Sk_.size)
        Sk_fit = Sk_[:n_fit]
        k2_fit = k_hat2_[:n_fit]
        if np.any(Sk_fit <= 0):
            return np.nan, np.nan, Sk_
        inv_S = 1.0 / Sk_fit
        A = np.column_stack([np.ones_like(k2_fit), k2_fit])
        coef, *_ = np.linalg.lstsq(A, inv_S, rcond=None)
        inv_chi, K_eff_ = coef
        if inv_chi <= 0 or K_eff_ <= 0:
            return np.nan, np.nan, Sk_
        chi_ = 1.0 / inv_chi
        return chi_, K_eff_, Sk_

    chi_central, K_central, Sk_central = _oz_from_samples(Phi_trajectory)
    if not np.isfinite(K_central):
        return Phi_mean_all, chi_central, np.nan, np.nan, np.nan, Sk_central

    # Jackknife
    block = n_samples // n_jk
    K_jk = []
    for b in range(n_jk):
        idx_keep = np.concatenate([
            np.arange(0, b * block),
            np.arange((b + 1) * block, n_samples),
        ])
        _, Kb, _ = _oz_from_samples(Phi_trajectory[idx_keep])
        if np.isfinite(Kb):
            K_jk.append(Kb)
    K_jk = np.array(K_jk)
    if K_jk.size >= 2:
        mean_jk = K_jk.mean()
        K_err = np.sqrt((K_jk.size - 1) / K_jk.size * np.sum((K_jk - mean_jk) ** 2))
    else:
        K_err = np.nan

    xi2 = K_central * chi_central
    return Phi_mean_all, chi_central, xi2, K_central, K_err, Sk_central


def run_mc_single(beta, m2, lam, J, N, rng_seed,
                  n_therm=2000, n_measure=4000, measure_every=2,
                  verbose=False):
    """Run a single H_1 MC simulation and return diagnostics."""
    rng = np.random.default_rng(rng_seed)
    s = rng.normal(0.0, 0.5, size=N)
    # Auto-tune step size during thermalisation
    sigma_prop = autotune_sigma(s, beta, m2, lam, J, rng,
                                 n_tune=min(500, n_therm // 2))
    for _ in range(n_therm):
        mc_sweep_checkerboard(s, beta, m2, lam, J, sigma_prop, rng)

    samples = []
    for t in range(n_measure):
        mc_sweep_checkerboard(s, beta, m2, lam, J, sigma_prop, rng)
        if t % measure_every == 0:
            samples.append(s * s)  # Phi_i = s_i^2
    Phi_traj = np.array(samples)
    return phi_correlator_stats(Phi_traj) + (sigma_prop,)


def run_mc_single_h3bilinear(beta, m2, lam, J, N, rng_seed,
                             n_therm=2000, n_measure=4000, measure_every=2):
    """Run a single H_3-bilinear MC simulation (bond = -J Σ (s_is_j)^2)."""
    rng = np.random.default_rng(rng_seed)
    s = rng.normal(0.0, 0.5, size=N)
    # Autotune by proxy: reuse the H_1 autotuner as rough heuristic then
    # refine with H_3-bilinear sweeps during the remainder of thermalisation.
    sigma = 0.5
    for _ in range(min(500, n_therm // 2)):
        s_before = s.copy()
        mc_sweep_checkerboard_h3bilinear(s, beta, m2, lam, J, sigma, rng)
        moved = np.mean(np.abs(s - s_before) > 1e-12)
        if moved > 0.45:
            sigma *= 1.1
        elif moved < 0.35:
            sigma *= 0.9
    for _ in range(n_therm):
        mc_sweep_checkerboard_h3bilinear(s, beta, m2, lam, J, sigma, rng)

    samples = []
    for t in range(n_measure):
        mc_sweep_checkerboard_h3bilinear(s, beta, m2, lam, J, sigma, rng)
        if t % measure_every == 0:
            samples.append(s * s)  # Phi_i = s_i^2
    Phi_traj = np.array(samples)
    return phi_correlator_stats(Phi_traj) + (sigma,)


# -----------------------------------------------------------------------------
# Scan and fit.
# -----------------------------------------------------------------------------

def power_law_fit_with_ci(x, y):
    """Least-squares fit y = C * x^p in log-log. Return p, C, stderr(p)."""
    mask = (x > 0) & (y > 0) & np.isfinite(x) & np.isfinite(y)
    lx = np.log(x[mask]); ly = np.log(y[mask])
    n = lx.size
    if n < 3:
        return np.nan, np.nan, np.nan, 0
    # Unweighted linear fit
    xbar = lx.mean(); ybar = ly.mean()
    Sxx = ((lx - xbar) ** 2).sum()
    Sxy = ((lx - xbar) * (ly - ybar)).sum()
    p = Sxy / Sxx
    intercept = ybar - p * xbar
    C = np.exp(intercept)
    # Residuals, standard error
    resid = ly - (p * lx + intercept)
    s2 = (resid ** 2).sum() / max(n - 2, 1)
    stderr_p = np.sqrt(s2 / Sxx)
    return p, C, stderr_p, n


def print_header(title, width=78):
    print()
    print("=" * width)
    print(f"  {title}")
    print("=" * width)


def print_row(cells, widths):
    print("  " + "  ".join(f"{str(c):>{w}s}" for c, w in zip(cells, widths)))


def main():
    print_header("OP-6 / M2-b  —  Envelope stiffness of H_1 (1D MC test)")
    print()
    print("  Hamiltonian: H_1 = Σ (m0^2/2) s^2 + (λ/4) s^4 − J Σ s_i s_j")
    print("  Goal: fit K_eff(<Φ>) = C · <Φ>^p; p=+1 is TGP target (α=2).")
    print()

    t_start = time.time()
    J = 1.0
    N = 1024

    # -----------------------------------------------------------------
    # Part A: Sanity check at J = 0
    # -----------------------------------------------------------------
    print_header("Part A  —  Sanity: J = 0 should give K_eff ~ 0  (no bond)")
    print()
    sanity = run_mc_single(beta=1.0, m2=1.0, lam=1.0, J=0.0, N=N, rng_seed=11,
                            n_therm=1500, n_measure=2000, measure_every=2)
    Phi_mean, chi, xi2, K, K_err, _, sigma_prop = sanity
    print(f"  J=0:   <Phi> = {Phi_mean:.4f}   χ = {chi:.4e}   "
          f"ξ² = {xi2:.4e}   K_eff = {K:.4e}")
    print(f"  (step sigma_prop = {sigma_prop:.3f})")
    print()
    # For J=0, OZ fit should fail (S(k) flat) → K = NaN. That is the right outcome.
    if not np.isfinite(K) or abs(K) < 1e-2:
        print("  ✓ J=0 gives negligible / undefined envelope stiffness, as expected.")
    else:
        print(f"  ⚠ J=0 gives K_eff = {K:.4e}. Statistical artifact; interpret part B / C with caution.")

    # -----------------------------------------------------------------
    # Part B: Gaussian sub-limit (λ = 0)  —  M2-a §2.v.a predicts p = −1
    # -----------------------------------------------------------------
    print_header("Part B  —  Gaussian sub-limit (λ = 0): expect p ≈ −1")
    print()
    print("  Vary (β, m0^2) with λ=0 to scan <Φ>.")
    print()
    print("  Exact 1D Gaussian identities (J=1):")
    print("    cosh(ν) = m²/2,  <Φ>_exact = 1/(2β sinh(ν)),")
    print("    K_exact = (β²/2) tanh(ν),  K*<Φ> = β/(2m²)  (exact, all m²).")
    print()
    print_row(["β", "m0^2", "<Φ>_MC", "<Φ>_ex", "K_MC", "K_err",
               "K_ex", "K_MC/K_ex"],
              [6, 6, 9, 9, 10, 10, 10, 10])
    print("  " + "-" * 86)

    gaussian_scan = []
    # For 1D Gaussian with J=1: ξ_Φ = 1/(2·arccosh(m²/2)). To get ξ_Φ ≥ 1 we
    # need m² ≤ 2·cosh(0.5) = 2.26. We scan m² ∈ [2.02, 2.25].
    beta_list_B = [0.5, 1.0, 2.0, 4.0]
    m2_list_B = [2.02, 2.05, 2.10, 2.20]
    for beta in beta_list_B:
        for m2 in m2_list_B:
            # Analytical
            nu = np.arccosh(m2 / (2.0 * J))
            Phi_ex = 1.0 / (2.0 * beta * J * np.sinh(nu))
            K_ex = (beta ** 2 * J ** 2 / 2.0) * np.tanh(nu)
            out = run_mc_single(beta=beta, m2=m2, lam=0.0, J=J, N=N,
                                rng_seed=int(1000 * beta + 100 * m2),
                                n_therm=3000, n_measure=8000, measure_every=2)
            Phi_mean, chi, xi2, K, K_err, _, _ = out
            gaussian_scan.append((beta, m2, Phi_mean, chi, xi2, K, K_err,
                                  Phi_ex, K_ex))
            ratio = K / K_ex if (np.isfinite(K) and K_ex > 0) else np.nan
            Kerr_str = f"{K_err:.2e}" if np.isfinite(K_err) else "nan"
            ratio_str = f"{ratio:.3f}" if np.isfinite(ratio) else "nan"
            print_row([f"{beta:.2f}", f"{m2:.2f}",
                       f"{Phi_mean:.4f}", f"{Phi_ex:.4f}",
                       f"{K:.4e}", Kerr_str,
                       f"{K_ex:.4e}", ratio_str],
                      [6, 6, 9, 9, 10, 10, 10, 10])
    gs = np.array(gaussian_scan)
    if gs.size > 0:
        Phi_arr = gs[:, 2]; K_arr = gs[:, 5]
        p_G, C_G, err_G, n_G = power_law_fit_with_ci(Phi_arr, K_arr)
        # Also fit the analytical K_ex vs Phi_ex in the same parameter range
        p_G_ex, C_G_ex, err_G_ex, _ = power_law_fit_with_ci(gs[:, 7], gs[:, 8])
        print()
        print(f"  Analytical (Gaussian, same range): K = {C_G_ex:.4f} * <Φ>^{p_G_ex:.3f}")
        print(f"  MC measurement:                    K = {C_G:.4f} * <Φ>^{p_G:.3f} "
              f"(± {err_G:.3f}, n = {n_G})")
        print()
        if abs(p_G - p_G_ex) < max(2 * err_G, 0.3):
            print(f"  ✓ MC agrees with analytical prediction ({p_G_ex:.2f}); "
                  "method validated.")
        else:
            print(f"  ! Discrepancy: MC p = {p_G:.2f} vs analytical {p_G_ex:.2f}.")

    # -----------------------------------------------------------------
    # Part C: Full φ^4 (λ > 0)  —  the actual M2-b test
    # -----------------------------------------------------------------
    print_header("Part C  —  Full H_1 (λ = 1): is K_eff ∝ Φ or not?")
    print()
    print("  Scan (β, m0^2) with λ = 1.0.")
    print("  Interpretation: p = +1 means K(Φ) ~ Φ ⇔ α = 2 (TGP target).")
    print()
    print_row(["β", "m0^2", "λ", "<Φ>", "χ_Φ", "ξ²_Φ", "K_eff", "±K_err"],
              [6, 6, 6, 10, 10, 10, 10, 10])
    print("  " + "-" * 82)

    full_scan = []
    # Choose m² at or below 2J so that at weak λ the theory is near-critical
    # and ξ_Φ is larger. λ=1 stabilises against runaway (m² can be negative).
    beta_list_C = [0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0]
    for beta in beta_list_C:
        for m2, lam in [(-1.0, 1.0), (0.0, 1.0), (1.0, 1.0), (2.0, 1.0), (3.0, 1.0)]:
            out = run_mc_single(beta=beta, m2=m2, lam=lam, J=J, N=N,
                                rng_seed=int(2000 * beta + 7 * (m2 + 2)),
                                n_therm=3000, n_measure=8000, measure_every=2)
            Phi_mean, chi, xi2, K, K_err, _, _ = out
            full_scan.append((beta, m2, lam, Phi_mean, chi, xi2, K, K_err))
            Kerr_str = f"{K_err:.2e}" if np.isfinite(K_err) else "nan"
            print_row([f"{beta:.2f}", f"{m2:.2f}", f"{lam:.2f}",
                       f"{Phi_mean:.4f}", f"{chi:.4f}",
                       f"{xi2:.4f}", f"{K:.4e}", Kerr_str],
                      [6, 6, 6, 10, 10, 10, 10, 10])

    fs = np.array(full_scan)
    Phi_arr = fs[:, 3]; K_arr = fs[:, 6]
    p_F, C_F, err_F, n_F = power_law_fit_with_ci(Phi_arr, K_arr)

    # -----------------------------------------------------------------
    # Part D: Verdict
    # -----------------------------------------------------------------
    print_header("Part D  —  Verdict")
    print()
    print(f"  Full H_1 fit:  K_eff = {C_F:.4f} * <Φ>^{p_F:.3f}  "
          f"(± {err_F:.3f}, n = {n_F})")
    print()
    # Test the key TGP hypothesis: p = +1 (α = 2).
    z_plus1 = (p_F - 1.0) / err_F if err_F > 0 else np.inf
    print(f"  Test of TGP target p = +1:  deviation = {p_F - 1.0:+.2f}, "
          f"  z = {z_plus1:+.1f} σ")
    print()
    print("  Decision table:")
    print("    p ≈ +1   → M2-b PASSES.  K ∝ Φ, α=2 from H_1 supported.")
    print("    p ≈  0   → K ~ const in Φ.  α(Φ)=1/2, α(φ)=1.  Not TGP target.")
    print("    p ≈ −1   → Gaussian-like.  α(Φ)=-1/2.  M1-B FAILS → M1-A.")
    print()
    if p_F > 0.7 and p_F < 1.3:
        print("  ✓ VERDICT: p ≈ +1. Envelope mechanism generates K ∝ Φ.")
        print("    M1-B vindicated. Proceed to M3 (3D FRG / MC confirmation).")
        verdict = "PASS"
    elif abs(p_F) < 0.35:
        print("  ~ VERDICT: p ≈ 0. K_eff roughly constant in Φ.")
        print("    Not the TGP target (+1).  Moderate evidence against M1-B.")
        verdict = "FLAT"
    elif p_F < -0.4:
        print("  ✗ VERDICT: p < 0. K_eff decreases with Φ.")
        print("    Envelope coarse-graining of H_1 does NOT produce α=2.")
        print("    M1-B fails. Recommended pivot: M1-A (change axiom to H_3).")
        verdict = "FAIL"
    else:
        print(f"  ? VERDICT: p = {p_F:.2f} does not match any expected regime.")
        verdict = "UNCLEAR"
    if abs(z_plus1) > 2.0:
        print(f"  → TGP target p = +1 rejected at {abs(z_plus1):.1f} σ.")

    # -----------------------------------------------------------------
    # Part E: H_3 bilinear bond (-J Σ (s_i s_j)^2)  —  M2c check
    # -----------------------------------------------------------------
    # NOTE (interpretation): M2c's algebraic derivation that
    #   −J Σ (φ_iφ_j)² = −J Σ Φ_iΦ_j ⇒ K(Φ) = const (p_Φ = 0)
    # concerns the *bare* bond stiffness of Φ in the continuum. The MC
    # here measures the *dressed* envelope stiffness via the composite
    # Φ = s² correlator, with all fluctuation corrections included. The
    # two agree qualitatively in the weak-coupling regime (where the
    # bare term dominates), but at the strong coupling needed to keep
    # the MC stable (λ > 4J) the OZ-extracted K(<Φ>) may show residual
    # Φ-dependence from non-Gaussian renormalization. The critical
    # distinction from H_GL (which gives p_Φ = +1 cleanly from the bare
    # bond) should still be visible: expect p_Φ ≪ +1 even if p_Φ is
    # not exactly 0. If p_Φ comes out close to +1, that would falsify
    # the M2c analytical claim and would need investigation.
    # -----------------------------------------------------------------
    print_header("Part E  —  H_3 bilinear (-J Σ (s_is_j)^2): expect p_Φ ≈ 0")
    print()
    print("  Hamiltonian: H_3_bil = Σ (m²/2)s² + (λ/4)s⁴ − J Σ (s_is_j)².")
    print("  Analytical prediction (M2c §1 / M1 §4.2 corrected):")
    print("     −J Σ (φ_iφ_j)² = −J Σ Φ_iΦ_j  ⇒  K(Φ) = const (bare)")
    print("  Equivalently: log-log slope p_Φ ≈ 0 in the Φ variable.")
    print("  Distinction from H_GL (p_Φ = +1): H_3 bilinear should give")
    print("  p_Φ ≪ +1 even after fluctuation renormalization.")
    print()
    print_row(["β", "m0^2", "λ", "<Φ>", "χ_Φ", "ξ²_Φ", "K_eff", "±K_err"],
              [6, 6, 6, 10, 10, 10, 10, 10])
    print("  " + "-" * 82)

    # Stability. The bond -J Σ Φ_iΦ_j is ferromagnetic in Φ ≥ 0, so the
    # mean-field energy per site is e(Φ) = (m²/2)Φ + (λ/4 − J·z/2)Φ² (z=2 in
    # 1D). Hence stability requires λ > 4J. We take λ = 6–10 with J = 1.
    h3_scan = []
    beta_list_E = [0.2, 0.3, 0.5, 0.7, 1.0]
    # (m², λ) pairs chosen to keep the system stable and span a range of <Φ>.
    params_E = [(0.5, 6.0), (1.0, 6.0), (2.0, 6.0), (3.0, 8.0), (4.0, 10.0)]
    for beta in beta_list_E:
        for m2E, lamE in params_E:
            out = run_mc_single_h3bilinear(
                beta=beta, m2=m2E, lam=lamE, J=J, N=N,
                rng_seed=int(5000 * beta + 13 * m2E + 31 * lamE),
                n_therm=3000, n_measure=8000, measure_every=2)
            Phi_mean, chi, xi2, K, K_err, _, _ = out
            h3_scan.append((beta, m2E, lamE, Phi_mean, chi, xi2, K, K_err))
            Kerr_str = f"{K_err:.2e}" if np.isfinite(K_err) else "nan"
            print_row([f"{beta:.2f}", f"{m2E:.2f}", f"{lamE:.2f}",
                       f"{Phi_mean:.4f}", f"{chi:.4f}",
                       f"{xi2:.4f}", f"{K:.4e}", Kerr_str],
                      [6, 6, 6, 10, 10, 10, 10, 10])

    hs = np.array(h3_scan)
    Phi_arr_E = hs[:, 3]; K_arr_E = hs[:, 6]
    p_E, C_E, err_E, n_E = power_law_fit_with_ci(Phi_arr_E, K_arr_E)
    print()
    print(f"  H_3 bilinear fit:  K_eff = {C_E:.4f} * <Φ>^{p_E:.3f}  "
          f"(± {err_E:.3f}, n = {n_E})")
    print()
    z_plus1_E = (p_E - 1.0) / err_E if err_E > 0 else np.inf
    z_zero_E = p_E / err_E if err_E > 0 else np.inf
    print(f"  Test of 'H_3 ⇒ α=2' (target p_Φ=+1): z = {z_plus1_E:+.1f} σ")
    print(f"  Test of M2c prediction (target p_Φ= 0):  z = {z_zero_E:+.1f} σ")
    print()
    if abs(p_E) < 0.4:
        print("  ✓ M2c prediction confirmed: bilinear H_3 has K(Φ) ≈ const (p_Φ ≈ 0).")
        print("    H_3 bilinear is NOT the α=2 target; H_GL is.")
        verdict_E = "CONFIRM_M2c"
    elif abs(p_E - 1.0) < 0.4:
        print("  ? Unexpected: p_Φ ≈ +1.  Would indicate bilinear H_3 generates α=2.")
        print("    Contradicts M2c analytical derivation; investigate.")
        verdict_E = "SURPRISE"
    else:
        print(f"  ? Intermediate: p_Φ = {p_E:.2f}, neither 0 nor +1 cleanly.")
        verdict_E = "UNCLEAR"

    elapsed = time.time() - t_start
    print()
    print(f"  Wall time: {elapsed:.1f} s   |   "
          f"H_1 verdict: {verdict}   |   H_3-bil verdict: {verdict_E}")
    print()
    print("  Saving raw scan data to m2b_scan_results.txt ...")
    with open("m2b_scan_results.txt", "w", encoding="utf-8") as f:
        f.write("# OP-6 / M2-b envelope stiffness scan\n")
        f.write("# Part B (Gaussian, λ=0):\n")
        f.write("# beta  m0^2   <Phi>_MC chi       xi2       K_MC      K_err     "
                "<Phi>_ex K_ex\n")
        for row in gaussian_scan:
            vs = [f"{v:10.4e}" if np.isfinite(v) else "       nan" for v in row]
            f.write("  " + "  ".join(vs) + "\n")
        f.write(f"# Gaussian MC fit:  K = {C_G:.4e} * Phi^{p_G:.4f} ± {err_G:.4f}\n")
        f.write(f"# Gaussian analytical fit (same param range): "
                f"K = {C_G_ex:.4e} * Phi^{p_G_ex:.4f}\n")
        f.write("\n# Part C (Full H_1, λ=1):\n")
        f.write("# beta  m0^2   lam    <Phi>    chi       xi2       K_eff     K_err\n")
        for row in full_scan:
            vs = [f"{v:10.4e}" if np.isfinite(v) else "       nan" for v in row]
            f.write("  " + "  ".join(vs) + "\n")
        f.write(f"# Full H_1 fit: K = {C_F:.4e} * Phi^{p_F:.4f} ± {err_F:.4f}\n")
        f.write(f"# H_1 verdict: {verdict}\n")
        f.write("\n# Part E (H_3 bilinear bond -J Σ (s_is_j)^2):\n")
        f.write("# beta  m0^2   lam    <Phi>    chi       xi2       K_eff     K_err\n")
        for row in h3_scan:
            vs = [f"{v:10.4e}" if np.isfinite(v) else "       nan" for v in row]
            f.write("  " + "  ".join(vs) + "\n")
        f.write(f"# H_3 bilinear fit: K = {C_E:.4e} * Phi^{p_E:.4f} ± {err_E:.4f}\n")
        f.write(f"# H_3 bilinear verdict: {verdict_E}\n")
    print("  Done.")


if __name__ == "__main__":
    main()
