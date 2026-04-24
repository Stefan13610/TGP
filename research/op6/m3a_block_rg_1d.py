#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
OP-6 / M3-a — Block-RG bond-form test of H_1 in 1D.
====================================================

See M3a_block_rg_1d_plan.md for the full protocol. Brief:

We simulate H_1 in 1D (same as m2b_envelope_stiffness_1d.py Part C)
but form block-averaged Phi fields

    Phi_{B,I} = (1/B) * sum_{i in block I} s_i^2     (block size B)

and measure the block-level Ornstein-Zernike stiffness K_B(<Phi_B>)
for a ladder of block sizes B in {1, 2, 4, 8, 16}. Then fit

    K_B(<Phi_B>) = C_B * <Phi_B>^{p_B}

and report p_B as a function of B.

Interpretation:
- p_B -> +1 as B grows  =>  block-RG of H_1 generates H_GL
                             (M1-B revival).
- p_B stays near 0/-0.3 across B  =>  K^{(0)} dominates K^{(1)}
                             (M3-c confirmed; M1-B closed).

Dependencies: numpy only.
"""

import sys, time, warnings
warnings.filterwarnings("ignore")
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import numpy as np


# -----------------------------------------------------------------------------
# Reuse MC infrastructure from m2b (duplicated here to keep this standalone).
# -----------------------------------------------------------------------------

def mc_sweep_checkerboard(s, beta, m2, lam, J, sigma_prop, rng):
    """One checkerboard Metropolis sweep for H_1 in 1D."""
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
        u = rng.random(idx.size)
        accept = (dE <= 0.0) | (u < np.exp(-beta * np.clip(dE, None, 50.0)))
        s[idx[accept]] = s_new[accept]
    return s


def autotune_sigma(s, beta, m2, lam, J, rng,
                   target_acc=0.4, n_tune=400, sigma0=0.5):
    sigma = sigma0
    for _ in range(n_tune):
        s_before = s.copy()
        mc_sweep_checkerboard(s, beta, m2, lam, J, sigma, rng)
        moved = np.mean(np.abs(s - s_before) > 1e-12)
        if moved > target_acc + 0.05:
            sigma *= 1.1
        elif moved < target_acc - 0.05:
            sigma *= 0.9
    return sigma


# -----------------------------------------------------------------------------
# Block-averaged OZ fit.
# -----------------------------------------------------------------------------

def block_average(Phi_trajectory, B):
    """
    Phi_trajectory: (n_samples, N) single-site Phi = s^2.
    Returns block-averaged trajectory (n_samples, N // B).
    Φ_B,I = (1/B) * Σ_{i in block I} Φ_i.
    """
    n_s, N = Phi_trajectory.shape
    N_B = N // B
    # Truncate N to multiple of B, then mean over block axis.
    return Phi_trajectory[:, :N_B * B].reshape(n_s, N_B, B).mean(axis=2)


def oz_fit_block(Phi_B_traj, n_k_fit=4, n_jk=8):
    """
    OZ fit on a block-level trajectory. Same protocol as m2b's
    phi_correlator_stats but applied to a pre-blocked field.

    Returns: (mean_Phi_B, chi_B, xi2_B, K_B, K_B_err).
    """
    n_samples, NB = Phi_B_traj.shape
    Phi_mean_all = Phi_B_traj.mean()

    def _oz_from_samples(traj):
        Phi_shifted = traj - Phi_mean_all
        F = np.fft.rfft(Phi_shifted, axis=1)
        Sk_ = (F * np.conj(F)).real.mean(axis=0) / NB
        j_ = np.arange(NB // 2 + 1)
        k_ = 2.0 * np.pi * j_ / NB
        k_hat2_ = 2.0 * (1.0 - np.cos(k_))
        n_fit = min(n_k_fit + 1, Sk_.size)
        Sk_fit = Sk_[:n_fit]
        k2_fit = k_hat2_[:n_fit]
        if np.any(Sk_fit <= 0):
            return np.nan, np.nan
        inv_S = 1.0 / Sk_fit
        A = np.column_stack([np.ones_like(k2_fit), k2_fit])
        coef, *_ = np.linalg.lstsq(A, inv_S, rcond=None)
        inv_chi, K_eff_ = coef
        if inv_chi <= 0 or K_eff_ <= 0:
            return np.nan, np.nan
        return 1.0 / inv_chi, K_eff_

    chi_central, K_central = _oz_from_samples(Phi_B_traj)
    if not np.isfinite(K_central):
        return Phi_mean_all, chi_central, np.nan, np.nan, np.nan

    block = n_samples // n_jk
    K_jk = []
    for b in range(n_jk):
        idx_keep = np.concatenate([
            np.arange(0, b * block),
            np.arange((b + 1) * block, n_samples),
        ])
        _, Kb = _oz_from_samples(Phi_B_traj[idx_keep])
        if np.isfinite(Kb):
            K_jk.append(Kb)
    K_jk = np.array(K_jk)
    if K_jk.size >= 2:
        mean_jk = K_jk.mean()
        K_err = np.sqrt((K_jk.size - 1) / K_jk.size * np.sum((K_jk - mean_jk) ** 2))
    else:
        K_err = np.nan

    xi2 = K_central * chi_central
    return Phi_mean_all, chi_central, xi2, K_central, K_err


def run_mc_and_block_scan(beta, m2, lam, J, N, block_sizes, rng_seed,
                           n_therm=3000, n_measure=8000, measure_every=2):
    """Run a single H_1 simulation, then for each block size B compute
    (<Phi_B>, chi_B, K_B, K_B_err) via OZ on the block-averaged field.
    Returns dict {B: (Phi_mean_B, chi_B, xi2_B, K_B, K_B_err)}.
    """
    rng = np.random.default_rng(rng_seed)
    s = rng.normal(0.0, 0.5, size=N)
    sigma_prop = autotune_sigma(s, beta, m2, lam, J, rng,
                                 n_tune=min(500, n_therm // 2))
    for _ in range(n_therm):
        mc_sweep_checkerboard(s, beta, m2, lam, J, sigma_prop, rng)

    samples = []
    for t in range(n_measure):
        mc_sweep_checkerboard(s, beta, m2, lam, J, sigma_prop, rng)
        if t % measure_every == 0:
            samples.append(s * s)
    Phi_traj = np.array(samples)

    results = {}
    for B in block_sizes:
        Phi_B = block_average(Phi_traj, B)
        results[B] = oz_fit_block(Phi_B)
    return results


# -----------------------------------------------------------------------------
# Power-law fit utility.
# -----------------------------------------------------------------------------

def power_law_fit_with_ci(x, y):
    mask = (x > 0) & (y > 0) & np.isfinite(x) & np.isfinite(y)
    lx = np.log(x[mask]); ly = np.log(y[mask])
    n = lx.size
    if n < 3:
        return np.nan, np.nan, np.nan, n
    xbar = lx.mean(); ybar = ly.mean()
    Sxx = ((lx - xbar) ** 2).sum()
    Sxy = ((lx - xbar) * (ly - ybar)).sum()
    p = Sxy / Sxx
    intercept = ybar - p * xbar
    C = np.exp(intercept)
    resid = ly - (p * lx + intercept)
    s2 = (resid ** 2).sum() / max(n - 2, 1)
    stderr_p = np.sqrt(s2 / Sxx)
    return p, C, stderr_p, n


# -----------------------------------------------------------------------------
# Main.
# -----------------------------------------------------------------------------

def print_header(title, width=78):
    print()
    print("=" * width)
    print(f"  {title}")
    print("=" * width)


def print_row(cells, widths):
    print("  " + "  ".join(f"{str(c):>{w}s}" for c, w in zip(cells, widths)))


def main():
    print_header("OP-6 / M3-a  —  Block-RG bond-form test of H_1 in 1D")
    print()
    print("  Hamiltonian: H_1 = Σ (m²/2) s² + (λ/4) s⁴ − J Σ s_is_j (same as M2-b).")
    print("  Method: run MC, form Φ_{B,I} = (1/B)·Σ_{i in block I} s_i²,")
    print("          measure K_B(<Φ_B>) via OZ at multiple block sizes B.")
    print("  Fit K_B = C_B · <Φ_B>^{p_B} per B, watch p_B(B).")
    print()

    t_start = time.time()
    J = 1.0
    lam = 1.0
    N = 4096
    block_sizes = [1, 2, 4, 8, 16]

    # Same grid as M2-b Part C so the B=1 result can be cross-checked.
    beta_list = [0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0]
    m2_list = [-1.0, 0.0, 1.0, 2.0, 3.0]

    # {B: [(beta, m2, Phi_mean_B, chi_B, xi2_B, K_B, K_B_err), ...]}
    per_block_rows = {B: [] for B in block_sizes}

    n_runs = len(beta_list) * len(m2_list)
    run_idx = 0
    for beta in beta_list:
        for m2 in m2_list:
            run_idx += 1
            t_run = time.time()
            results = run_mc_and_block_scan(
                beta=beta, m2=m2, lam=lam, J=J, N=N,
                block_sizes=block_sizes,
                rng_seed=int(3000 * beta + 17 * (m2 + 2)),
                n_therm=3000, n_measure=8000, measure_every=2,
            )
            dt = time.time() - t_run
            print(f"  [run {run_idx}/{n_runs}]  β={beta:.2f}  m²={m2:+.2f}  "
                  f"({dt:.1f}s)")
            for B in block_sizes:
                Phi_mean_B, chi_B, xi2_B, K_B, K_B_err = results[B]
                per_block_rows[B].append(
                    (beta, m2, Phi_mean_B, chi_B, xi2_B, K_B, K_B_err)
                )

    # ---------------------------------------------------------------------
    # Per-block table + fit
    # ---------------------------------------------------------------------
    print_header("Per-block fits")
    print()
    print_row(["B", "n_valid", "p_B", "±err(p_B)", "C_B", "note"],
              [4, 8, 10, 10, 12, 30])
    print("  " + "-" * 76)

    fit_summary = {}
    for B in block_sizes:
        rows = np.array(per_block_rows[B])
        Phi_arr = rows[:, 2]
        K_arr = rows[:, 5]
        p_B, C_B, err_B, n_B = power_law_fit_with_ci(Phi_arr, K_arr)
        fit_summary[B] = (p_B, C_B, err_B, n_B)
        # Interpretive note
        if np.isfinite(p_B):
            if p_B > 0.7:
                note = "trending toward +1 (M1-B?)"
            elif abs(p_B) < 0.35:
                note = "flat ~0 (K indep. of Phi)"
            elif p_B < -0.4:
                note = "Gaussian-like"
            else:
                note = ""
        else:
            note = "fit failed"
        p_str = f"{p_B:+.3f}" if np.isfinite(p_B) else "nan"
        err_str = f"{err_B:.3f}" if np.isfinite(err_B) else "nan"
        C_str = f"{C_B:.4e}" if np.isfinite(C_B) else "nan"
        print_row([str(B), str(n_B), p_str, err_str, C_str, note],
                  [4, 8, 10, 10, 12, 30])

    # ---------------------------------------------------------------------
    # Flow table
    # ---------------------------------------------------------------------
    print_header("Block-RG flow of p_B")
    print()
    print("  Trend:  p_1 → p_16 traces the running of the exponent under coarse-graining.")
    print("          +1 target = TGP (M1-B). Flat ≈ 0 or negative = M3-c prediction.")
    print()
    ps = [fit_summary[B][0] for B in block_sizes]
    errs = [fit_summary[B][2] for B in block_sizes]
    print("  B:   " + "   ".join(f"{B:>6d}" for B in block_sizes))
    print("  p_B: " + "   ".join(
        f"{p:+6.3f}" if np.isfinite(p) else "   nan" for p in ps))
    print("  err: " + "   ".join(
        f"{e:6.3f}" if np.isfinite(e) else "   nan" for e in errs))
    print()

    # Verdict
    valid_ps = [(B, p, e) for B, p, e in zip(block_sizes, ps, errs)
                if np.isfinite(p) and np.isfinite(e)]
    if len(valid_ps) >= 2:
        B_lo, p_lo, e_lo = valid_ps[0]
        B_hi, p_hi, e_hi = valid_ps[-1]
        drift = p_hi - p_lo
        drift_err = np.sqrt(e_lo ** 2 + e_hi ** 2)
        drift_sigma = drift / drift_err if drift_err > 0 else np.inf

        z_hi_plus1 = (p_hi - 1.0) / e_hi if e_hi > 0 else np.inf
        z_hi_zero = p_hi / e_hi if e_hi > 0 else np.inf

        print(f"  Drift p(B={B_hi}) − p(B={B_lo}) = {drift:+.3f} "
              f"± {drift_err:.3f}  ({drift_sigma:+.1f}σ)")
        print(f"  Largest-B test of p = +1:  z = {z_hi_plus1:+.1f}σ")
        print(f"  Largest-B test of p =  0:  z = {z_hi_zero:+.1f}σ")
        print()

        if p_hi > 0.6 and z_hi_plus1 > -2.0:
            verdict = "M1_B_REVIVED"
            print("  ✓ VERDICT: p_B trending to +1 under block-RG. "
                  "M1-B possibly revived.")
            print("    Recommend: proceed to M3-b (3D cluster MC) to confirm "
                  "in Ising universality class.")
        elif abs(p_hi) < 0.5 and abs(drift_sigma) < 2.0:
            verdict = "M3_C_CONFIRMED"
            print("  ✓ VERDICT: p_B stays near 0/-0.3 across block-RG ladder. "
                  "K^{(0)} dominates.")
            print("    M3-c prediction confirmed in 1D. M1-B analytically + "
                  "numerically closed.")
            print("    Recommend: reopen axiom conversation (M1-A' / axiom "
                  "modification).")
        elif drift_sigma > 2.0 and p_hi > 0.4:
            verdict = "DRIFTING_UP"
            print("  ~ VERDICT: p_B drifting upward with B but not reaching +1.")
            print("    Inconclusive; proceed to M3-b (3D) to see if drift "
                  "continues / completes.")
        else:
            verdict = "UNCLEAR"
            print("  ? VERDICT: non-monotone or unclear flow pattern.")

    # ---------------------------------------------------------------------
    # Save raw data
    # ---------------------------------------------------------------------
    elapsed = time.time() - t_start
    print()
    print(f"  Wall time: {elapsed:.1f} s   |   Verdict: {verdict}")
    print()
    print("  Saving raw data to m3a_block_rg_results.txt ...")
    with open("m3a_block_rg_results.txt", "w", encoding="utf-8") as f:
        f.write("# OP-6 / M3-a block-RG bond-form test of H_1 in 1D\n")
        f.write(f"# N = {N}, lam = {lam}, J = {J}\n")
        f.write(f"# beta_list = {beta_list}\n")
        f.write(f"# m2_list = {m2_list}\n")
        f.write(f"# block_sizes = {block_sizes}\n")
        f.write("\n")
        for B in block_sizes:
            f.write(f"# Block size B = {B}\n")
            f.write("# beta   m2      <Phi_B>   chi_B     xi2_B     K_B       K_err\n")
            for row in per_block_rows[B]:
                vs = [f"{v:10.4e}" if np.isfinite(v) else "       nan" for v in row]
                f.write("  " + "  ".join(vs) + "\n")
            p_B, C_B, err_B, n_B = fit_summary[B]
            f.write(f"# B = {B} fit: K_B = {C_B:.4e} * <Phi_B>^{p_B:.4f} "
                    f"± {err_B:.4f}  (n = {n_B})\n\n")
        f.write(f"# Verdict: {verdict}\n")
    print("  Done.")


if __name__ == "__main__":
    main()
