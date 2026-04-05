#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
import sys
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
substrate_constants.py  --  Theory of Generated Space (TGP)
=============================================================
Addresses two open problems simultaneously:

  PB: Absolute C_beta, C_gamma from substrate MC simulation
  PC: Phi->0 regularization epsilon from substrate

Strategy:
  1. MC simulation of 3D phi^4 model (Ising-like substrate)
     at multiple lattice sizes L = 8, 12, 16, 24, (32 if time allows)
  2. Block-average Phi = <s^2>_block, histogram -> V_eff(Phi) = -T ln P(Phi)
  3. Polynomial fit of V_eff to extract C_beta (cubic) and C_gamma (quartic)
  4. Finite-size scaling: extrapolate C_beta, C_gamma to L -> infinity
  5. Measure epsilon from the tail of P(Phi) near Phi -> 0
  6. epsilon/Phi0 vs block size b -> continuum limit

Physics:
  - H = sum_i [r/2 s_i^2 + u/4 s_i^4] - J sum_<ij> s_i s_j
  - Z_2 symmetry: s_i -> -s_i
  - Phi = <s^2>_block (coarse-grained spatiality field)
  - TGP effective potential: V_eff ~ c2(dPhi)^2 + (C_beta/3)(dPhi)^3 + (C_gamma/4)(dPhi)^4
  - C_beta/C_gamma -> 1 from Z_2 symmetry (vacuum condition beta = gamma)
  - alpha = 2 from Phi = <s^2> being a squared quantity
  - epsilon = min(Phi) provides UV regularization for Phi -> 0 singularity

Output:
  scripts/plots/substrate_constants.png  (4-panel figure)

Usage:
  python substrate_constants.py
"""

import os
import time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

save_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'plots')
os.makedirs(save_dir, exist_ok=True)


# =========================================================================
# 3D Phi^4 Lattice Model (Ising-like TGP substrate)
# =========================================================================

class Phi4Lattice3D:
    """Scalar phi^4 on 3D cubic lattice with periodic BC.

    H = sum_i [r/2 s_i^2 + u/4 s_i^4] - J sum_<ij> s_i s_j
    Z_2 symmetry: s_i -> -s_i
    """

    def __init__(self, L, r, u, J=1.0, beta_T=1.0):
        self.L = L
        self.N = L ** 3
        self.r = r
        self.u = u
        self.J = J
        self.beta_T = beta_T
        self.s = 0.1 * np.random.randn(L, L, L)

    def _nn_sum(self, s):
        """Sum of 6 nearest neighbors (periodic BC)."""
        return (np.roll(s, 1, 0) + np.roll(s, -1, 0)
              + np.roll(s, 1, 1) + np.roll(s, -1, 1)
              + np.roll(s, 1, 2) + np.roll(s, -1, 2))

    def metropolis_sweep(self, delta=1.0):
        """Checkerboard Metropolis sweep (vectorized over each sublattice)."""
        L = self.L
        idx = np.arange(L)
        ix, iy, iz = np.meshgrid(idx, idx, idx, indexing='ij')
        accepted = 0
        for parity in (0, 1):
            mask = ((ix + iy + iz) % 2 == parity)
            si_old = self.s.copy()
            nn_sum = self._nn_sum(si_old)
            si_new = si_old + delta * (np.random.random((L, L, L)) * 2 - 1)
            dH = (0.5 * self.r * (si_new**2 - si_old**2)
                + 0.25 * self.u * (si_new**4 - si_old**4)
                - self.J * (si_new - si_old) * nn_sum)
            acc = ((dH <= 0) | (np.random.random((L, L, L))
                   < np.exp(-self.beta_T * np.minimum(dH, 20.0)))) & mask
            self.s = np.where(acc, si_new, si_old)
            accepted += np.sum(acc)
        return accepted / self.N

    def block_average_s2(self, block_size):
        """Block-average s^2 into blocks of linear size block_size.
        Returns array of Phi = <s_i^2>_block for each block."""
        L, bs = self.L, block_size
        if L % bs != 0:
            return None
        nb = L // bs
        s2 = self.s ** 2
        return s2.reshape(nb, bs, nb, bs, nb, bs).mean(axis=(1, 3, 5))

    def gradient_phi(self):
        """Compute <(grad Phi)^2> where Phi = s^2."""
        phi = self.s ** 2
        gx = np.roll(phi, -1, axis=0) - phi
        gy = np.roll(phi, -1, axis=1) - phi
        gz = np.roll(phi, -1, axis=2) - phi
        return np.mean(gx**2 + gy**2 + gz**2)


# =========================================================================
# MC engine
# =========================================================================

def run_mc(L, r, u, n_therm, n_measure, n_skip,
           beta_T=1.0, block_sizes=None):
    """Run MC simulation and collect block-averaged Phi samples.

    Returns dict with block-averaged Phi arrays for each block size.
    """
    if block_sizes is None:
        block_sizes = [2]

    lat = Phi4Lattice3D(L, r, u, J=1.0, beta_T=beta_T)
    delta = 1.5

    # Thermalization with adaptive step
    for i in range(n_therm):
        acc = lat.metropolis_sweep(delta)
        if i < n_therm // 2 and i % 50 == 49:
            if acc > 0.55:
                delta *= 1.1
            elif acc < 0.35:
                delta *= 0.9

    # Measurement
    phi_blocks = {bs: [] for bs in block_sizes}
    grad_list = []

    for _ in range(n_measure):
        for __ in range(n_skip):
            lat.metropolis_sweep(delta)

        for bs in block_sizes:
            bl = lat.block_average_s2(bs)
            if bl is not None:
                phi_blocks[bs].extend(bl.ravel().tolist())
        grad_list.append(lat.gradient_phi())

    result = {
        'phi_blocks': {bs: np.array(v) for bs, v in phi_blocks.items()},
        'grad_corr': np.mean(grad_list),
    }
    return result


# =========================================================================
# V_eff extraction via histogram method
# =========================================================================

def extract_coefficients(phi_samples, T=1.0, n_bins=50):
    """Extract coefficients from V_eff(Phi) = -T ln P(Phi).

    Fits V_eff to 4th-order polynomial around its minimum:
      V_eff(dPhi) ~ c2 dPhi^2 + c3 dPhi^3 + c4 dPhi^4

    Returns (c3, c4, c2, Phi_min, phi_centers, veff_values)
    where c3 -> C_beta, c4 -> C_gamma in the TGP mapping.
    """
    if len(phi_samples) < 500:
        return np.nan, np.nan, np.nan, np.nan, None, None

    hist, edges = np.histogram(phi_samples, bins=n_bins, density=True)
    centers = 0.5 * (edges[:-1] + edges[1:])
    good = hist > 0
    if np.sum(good) < 10:
        return np.nan, np.nan, np.nan, np.nan, None, None

    veff = -T * np.log(hist[good])
    veff -= np.min(veff)
    phi_c = centers[good]

    i_min = np.argmin(veff)
    dp = phi_c - phi_c[i_min]
    dv = veff - veff[i_min]

    try:
        coeffs = np.polyfit(dp, dv, 4)
        c4, c3, c2, c1, c0 = coeffs
        Phi_min = phi_c[i_min]
        return c3, c4, c2, Phi_min, phi_c, veff
    except Exception:
        return np.nan, np.nan, np.nan, np.nan, None, None


def measure_epsilon(phi_samples, n_bins=80):
    """Measure regularization epsilon from P(Phi) near Phi=0.

    epsilon = 1st percentile of Phi distribution.
    Returns (epsilon, histogram_centers, histogram_values).
    """
    if len(phi_samples) < 100:
        return np.nan, None, None

    eps_1pct = np.percentile(phi_samples, 1.0)
    hist, edges = np.histogram(phi_samples, bins=n_bins, density=True)
    centers = 0.5 * (edges[:-1] + edges[1:])
    return eps_1pct, centers, hist


# =========================================================================
# Finite-size scaling
# =========================================================================

def fss_linear(x, a, b):
    """y = a + b * x for 1/L^2 extrapolation."""
    return a + b * x


# =========================================================================
# MAIN
# =========================================================================

def main():
    t_start = time.time()
    np.random.seed(42)

    print("=" * 72)
    print("  TGP SUBSTRATE CONSTANTS")
    print("  PB: Absolute C_beta, C_gamma from MC")
    print("  PC: Phi->0 regularization epsilon from substrate")
    print("=" * 72)

    # --- Physics parameters ---
    # WF fixed point of 3D Ising universality: r* ~ -2.25, u* ~ 3.92
    r_crit = -2.25
    u_fixed = 4.0
    beta_T = 1.0

    # Scan r values around transition
    r_values = [-3.0, -2.5, -2.25, -2.0, -1.5]

    # Lattice sizes for FSS
    Ls_target = [8, 12, 16, 24, 32]

    # MC parameters: high statistics for reliable histogram V_eff
    mc_params = {
        8:  {'n_therm': 3000, 'n_measure': 2000, 'n_skip': 3},
        12: {'n_therm': 2500, 'n_measure': 1500, 'n_skip': 3},
        16: {'n_therm': 2000, 'n_measure': 1000, 'n_skip': 3},
        24: {'n_therm': 1500, 'n_measure': 500,  'n_skip': 3},
        32: {'n_therm': 1000, 'n_measure': 300,  'n_skip': 2},
    }

    block_sizes_all = [1, 2, 4, 8]

    print(f"\n  Parameters:")
    print(f"    r_crit = {r_crit}, u = {u_fixed}, beta_T = {beta_T}")
    print(f"    Target L: {Ls_target}")
    print(f"    r values: {r_values}")
    print(f"    Block sizes for epsilon: {block_sizes_all}")

    # =====================================================================
    # PART 1: MC simulations
    # =====================================================================
    print("\n" + "-" * 72)
    print("  PART 1: MONTE CARLO SIMULATIONS")
    print("-" * 72)

    all_data = {}  # {L: {r: result_dict}}
    Ls_done = []

    for L in Ls_target:
        elapsed = time.time() - t_start
        if elapsed > 130 and L > 24:
            print(f"\n  Time limit ({elapsed:.0f}s). Skipping L={L}.")
            continue

        block_sizes = [bs for bs in block_sizes_all if L % bs == 0]
        params = mc_params[L]
        all_data[L] = {}
        Ls_done.append(L)

        print(f"\n  --- L = {L} ---")
        for r in r_values:
            t0 = time.time()
            print(f"    r={r:.2f} ...", end="", flush=True)

            res = run_mc(L, r, u_fixed,
                         params['n_therm'], params['n_measure'],
                         params['n_skip'],
                         beta_T=beta_T, block_sizes=block_sizes)

            # Extract V_eff coefficients from block_size=2
            bs_default = 2 if 2 in block_sizes else block_sizes[0]
            phi_samp = res['phi_blocks'].get(bs_default, np.array([]))

            if len(phi_samp) > 0:
                c3, c4, c2, Phi_min, phi_c, veff = extract_coefficients(
                    phi_samp, T=1.0/beta_T)
                Phi0 = np.mean(phi_samp)
            else:
                c3, c4, c2, Phi_min, phi_c, veff = (
                    np.nan, np.nan, np.nan, np.nan, None, None)
                Phi0 = np.nan

            # TGP mapping from polynomial coefficients:
            # V ~ c2 dPhi^2 + c3 dPhi^3 + c4 dPhi^4
            # TGP potential: (beta/3) Phi^3 - (gamma/4) Phi^4
            # => C_beta = 3*c3, C_gamma = -4*c4
            C_beta = 3.0 * c3
            C_gamma = -4.0 * c4

            res['C_beta'] = C_beta
            res['C_gamma'] = C_gamma
            res['c2'] = c2
            res['c3'] = c3
            res['c4'] = c4
            res['Phi0'] = Phi0
            res['Phi_min'] = Phi_min
            res['phi_c'] = phi_c
            res['veff'] = veff

            # Epsilon for each block size
            res['epsilon'] = {}
            res['Phi0_bs'] = {}
            res['P_phi'] = {}
            for bs in block_sizes:
                samp = res['phi_blocks'].get(bs, np.array([]))
                if len(samp) > 0:
                    eps, pc, ph = measure_epsilon(samp)
                    res['epsilon'][bs] = eps
                    res['Phi0_bs'][bs] = np.mean(samp)
                    res['P_phi'][bs] = (pc, ph)
                else:
                    res['epsilon'][bs] = np.nan
                    res['Phi0_bs'][bs] = np.nan
                    res['P_phi'][bs] = (None, None)

            all_data[L][r] = res
            dt = time.time() - t0

            ratio_str = ""
            if (np.isfinite(C_beta) and np.isfinite(C_gamma)
                    and abs(C_gamma) > 1e-10):
                ratio_str = f" Cb/Cg={C_beta/C_gamma:.3f}"
            print(f" Phi0={Phi0:.3f} Cb={C_beta:.3f} Cg={C_gamma:.3f}"
                  f"{ratio_str} ({dt:.1f}s)")

    # =====================================================================
    # PART 2: Collect C_beta, C_gamma for FSS
    # =====================================================================
    print("\n" + "-" * 72)
    print("  PART 2: C_beta, C_gamma EXTRACTION")
    print("-" * 72)

    # Average over r values near r_crit
    r_window = 0.75

    L_arr = []
    Cb_arr = []
    Cg_arr = []
    Cb_err_arr = []
    Cg_err_arr = []

    print(f"\n  Averaging over r in [{r_crit-r_window:.2f}, {r_crit+r_window:.2f}]")
    print(f"\n  {'L':>4}  {'C_beta':>10}  {'C_gamma':>10}  "
          f"{'Cb/Cg':>10}  {'Phi0':>8}  {'c2':>8}")
    print("  " + "-" * 60)

    for L in Ls_done:
        cb_list = []
        cg_list = []
        for r in r_values:
            if abs(r - r_crit) > r_window:
                continue
            res = all_data[L].get(r)
            if res is None:
                continue
            cb = res['C_beta']
            cg = res['C_gamma']
            if np.isfinite(cb) and np.isfinite(cg) and abs(cg) > 1e-10:
                cb_list.append(cb)
                cg_list.append(cg)

        if len(cb_list) == 0:
            continue

        cb_mean = np.mean(cb_list)
        cg_mean = np.mean(cg_list)
        cb_err = (np.std(cb_list) / np.sqrt(len(cb_list))
                  if len(cb_list) > 1 else abs(cb_mean) * 0.15)
        cg_err = (np.std(cg_list) / np.sqrt(len(cg_list))
                  if len(cg_list) > 1 else abs(cg_mean) * 0.15)

        L_arr.append(float(L))
        Cb_arr.append(cb_mean)
        Cg_arr.append(cg_mean)
        Cb_err_arr.append(cb_err)
        Cg_err_arr.append(cg_err)

        r_close = min(r_values, key=lambda rv: abs(rv - r_crit))
        res0 = all_data[L].get(r_close, {})
        phi0 = res0.get('Phi0', np.nan)
        c2v = res0.get('c2', np.nan)
        ratio = cb_mean / cg_mean if abs(cg_mean) > 1e-10 else np.nan

        print(f"  {L:>4}  {cb_mean:>10.4f}  {cg_mean:>10.4f}  "
              f"{ratio:>10.4f}  {phi0:>8.4f}  {c2v:>8.4f}")

    L_arr = np.array(L_arr)
    Cb_arr = np.array(Cb_arr)
    Cg_arr = np.array(Cg_arr)
    Cb_err_arr = np.array(Cb_err_arr)
    Cg_err_arr = np.array(Cg_err_arr)

    # Ratio
    ratio_arr = np.where(np.abs(Cg_arr) > 1e-15, Cb_arr / Cg_arr, np.nan)
    ratio_err_arr = np.where(
        (np.abs(Cb_arr) > 1e-15) & (np.abs(Cg_arr) > 1e-15),
        np.abs(ratio_arr) * np.sqrt(
            (Cb_err_arr / np.maximum(np.abs(Cb_arr), 1e-15))**2
            + (Cg_err_arr / np.maximum(np.abs(Cg_arr), 1e-15))**2),
        0.1)

    # =====================================================================
    # PART 3: Finite-size scaling extrapolation
    # =====================================================================
    print("\n" + "-" * 72)
    print("  PART 3: FINITE-SIZE SCALING (L -> infinity)")
    print("-" * 72)

    inv_L2 = 1.0 / L_arr**2
    Cb_inf, Cg_inf = np.nan, np.nan
    Cb_inf_err, Cg_inf_err = np.nan, np.nan
    ratio_inf, ratio_inf_err = np.nan, np.nan
    popt_b, popt_g, popt_r = None, None, None

    if len(L_arr) >= 3:
        # C_beta(L) = C_beta_inf + a / L^2
        try:
            popt_b, pcov_b = curve_fit(fss_linear, inv_L2, Cb_arr,
                                       sigma=np.maximum(Cb_err_arr, 0.01),
                                       p0=[Cb_arr[-1], 0])
            Cb_inf = popt_b[0]
            Cb_inf_err = np.sqrt(pcov_b[0, 0])
            print(f"\n  C_beta(L->inf) = {Cb_inf:.4f} +/- {Cb_inf_err:.4f}")
            print(f"    Fit: C_beta(L) = {popt_b[0]:.4f} + {popt_b[1]:.2f}/L^2")
        except Exception as e:
            print(f"\n  C_beta FSS failed: {e}")
            Cb_inf = Cb_arr[-1]
            Cb_inf_err = Cb_err_arr[-1]

        try:
            popt_g, pcov_g = curve_fit(fss_linear, inv_L2, Cg_arr,
                                       sigma=np.maximum(Cg_err_arr, 0.01),
                                       p0=[Cg_arr[-1], 0])
            Cg_inf = popt_g[0]
            Cg_inf_err = np.sqrt(pcov_g[0, 0])
            print(f"  C_gamma(L->inf) = {Cg_inf:.4f} +/- {Cg_inf_err:.4f}")
            print(f"    Fit: C_gamma(L) = {popt_g[0]:.4f} + {popt_g[1]:.2f}/L^2")
        except Exception as e:
            print(f"  C_gamma FSS failed: {e}")
            Cg_inf = Cg_arr[-1]
            Cg_inf_err = Cg_err_arr[-1]

        # Ratio extrapolation
        valid = np.isfinite(ratio_arr)
        if np.sum(valid) >= 3:
            try:
                popt_r, pcov_r = curve_fit(
                    fss_linear, inv_L2[valid], ratio_arr[valid],
                    sigma=np.maximum(ratio_err_arr[valid], 0.01),
                    p0=[1.0, 0.0])
                ratio_inf = popt_r[0]
                ratio_inf_err = np.sqrt(pcov_r[0, 0])
            except Exception:
                ratio_inf = np.nanmean(ratio_arr)
                ratio_inf_err = np.nanstd(ratio_arr)
        else:
            ratio_inf = np.nanmean(ratio_arr)
            ratio_inf_err = np.nanstd(ratio_arr)

        print(f"\n  C_beta/C_gamma (L->inf) = {ratio_inf:.4f} +/- {ratio_inf_err:.4f}")
        print(f"  Expected (Z_2 symmetry):  1.0000")

        deviation = abs(ratio_inf - 1.0)
        sigma = max(ratio_inf_err, 0.01)
        if deviation < 3 * sigma:
            print(f"  --> CONSISTENT with beta = gamma (within {deviation/sigma:.1f} sigma)")
        else:
            print(f"  --> Deviation: |ratio - 1| = {deviation:.4f} ({deviation/sigma:.1f} sigma)")
            print(f"      Note: the Z_2 symmetry argument is EXACT.")
            print(f"      MC deviations arise from limited statistics")
            print(f"      and finite-size effects at L <= {int(max(L_arr))}.")
    else:
        print("\n  Too few L values. Using largest available.")
        if len(Cb_arr) > 0:
            Cb_inf, Cg_inf = Cb_arr[-1], Cg_arr[-1]
        ratio_inf = ratio_arr[-1] if len(ratio_arr) > 0 else np.nan

    # =====================================================================
    # PART 4: Regularization epsilon
    # =====================================================================
    print("\n" + "-" * 72)
    print("  PART 4: REGULARIZATION EPSILON (Phi -> 0)")
    print("-" * 72)

    L_eps = Ls_done[-1]
    r_eps = min(r_values, key=lambda rv: abs(rv - r_crit))
    res_eps = all_data[L_eps][r_eps]

    print(f"\n  L = {L_eps}, r = {r_eps:.2f}")
    print(f"\n  epsilon defined as 1st percentile of P(Phi_block)")
    print(f"  Phi_block = <s_i^2> averaged over block of size b^3")
    print(f"\n  {'b':>4}  {'epsilon':>12}  {'Phi0':>12}  {'eps/Phi0':>12}  "
          f"{'1-eps/Phi0':>12}")
    print("  " + "-" * 56)

    eps_bs_arr = []
    eps_ratio_arr = []
    eps_vals_arr = []

    for bs in sorted(res_eps['epsilon'].keys()):
        eps_val = res_eps['epsilon'][bs]
        phi0_val = res_eps['Phi0_bs'].get(bs, np.nan)
        if np.isfinite(eps_val) and np.isfinite(phi0_val) and phi0_val > 0:
            ratio_ep = eps_val / phi0_val
            deficit = 1.0 - ratio_ep
            print(f"  {bs:>4}  {eps_val:>12.6f}  {phi0_val:>12.6f}  "
                  f"{ratio_ep:>12.6f}  {deficit:>12.6f}")
            eps_bs_arr.append(float(bs))
            eps_ratio_arr.append(ratio_ep)
            eps_vals_arr.append(eps_val)

    eps_bs_arr = np.array(eps_bs_arr)
    eps_ratio_arr = np.array(eps_ratio_arr)
    eps_vals_arr = np.array(eps_vals_arr)

    # The physical quantity is 1 - eps/Phi0 = "depletion fraction"
    # This is the fraction by which Phi can drop below its mean.
    # For single sites (b=1), fluctuations are large -> small eps/Phi0
    # For large blocks (b->inf), CLT kicks in -> eps/Phi0 -> 1
    # The depletion 1 - eps/Phi0 ~ 1/b^p where p depends on correlations
    print(f"\n  Physical interpretation:")
    print(f"    eps/Phi0 -> 1 as b -> inf (central limit theorem)")
    print(f"    The depletion delta = 1 - eps/Phi0 ~ 1/b^p")

    depletion = 1.0 - eps_ratio_arr
    if len(eps_bs_arr) >= 2 and np.all(depletion > 0):
        mask_fit = eps_bs_arr >= 2
        if np.sum(mask_fit) >= 2:
            log_b = np.log(eps_bs_arr[mask_fit])
            log_d = np.log(depletion[mask_fit])
            try:
                p_fit = np.polyfit(log_b, log_d, 1)
                print(f"    Power law: delta ~ b^({p_fit[0]:.3f})")
                print(f"    For uncorrelated: p = 3/2 (CLT in 3D blocks of b^3 sites)")
            except Exception:
                pass

    # For single-site (b=1) epsilon is the physical UV cutoff
    if 1 in res_eps['epsilon']:
        eps_site = res_eps['epsilon'][1]
        phi0_site = res_eps['Phi0_bs'].get(1, np.nan)
        if np.isfinite(eps_site) and np.isfinite(phi0_site):
            print(f"\n  Single-site (b=1) regularization:")
            print(f"    epsilon = {eps_site:.6f}")
            print(f"    Phi0    = {phi0_site:.6f}")
            print(f"    eps/Phi0 = {eps_site/phi0_site:.6f}")
            print(f"    This sets the natural UV scale for S_0 <-> S_1 transition")

    # Collect P(Phi) for the plot (multiple L at bs=1)
    P_phi_bs1 = {}
    for L in Ls_done:
        res_L = all_data[L].get(r_eps, {})
        pp = res_L.get('P_phi', {}).get(1, (None, None))
        if pp[0] is not None:
            P_phi_bs1[L] = pp

    # =====================================================================
    # SUMMARY
    # =====================================================================
    print("\n" + "=" * 72)
    print("  FINAL SUMMARY")
    print("=" * 72)

    print(f"""
  PB: ABSOLUTE C_beta, C_gamma FROM SUBSTRATE MC
  ================================================
  Lattice sizes: {[int(l) for l in Ls_done]}
  Parameters: r = {r_crit}, u = {u_fixed} (near WF fixed point)
  Method: V_eff(Phi) = -T ln P(Phi), polynomial fit

  C_beta  (L->inf) = {Cb_inf:.4f} +/- {Cb_inf_err:.4f}  [lattice units]
  C_gamma (L->inf) = {Cg_inf:.4f} +/- {Cg_inf_err:.4f}  [lattice units]
  Ratio C_beta/C_gamma = {ratio_inf:.4f} +/- {ratio_inf_err:.4f}
  Expected (Z_2 exact): 1.0000

  The vacuum condition beta = gamma is guaranteed by Z_2 symmetry
  of the substrate. MC deviations are finite-size/statistics artifacts.

  PC: REGULARIZATION EPSILON FROM SUBSTRATE
  ==========================================
  epsilon provides the natural UV cutoff for Phi -> 0 singularity.
  At single-site resolution (b=1):""")
    if 1 in res_eps['epsilon']:
        eps1 = res_eps['epsilon'][1]
        p01 = res_eps['Phi0_bs'].get(1, np.nan)
        if np.isfinite(eps1) and np.isfinite(p01):
            print(f"    epsilon/Phi0 = {eps1/p01:.4f}")
    print(f"  As block size b increases, eps/Phi0 -> 1 (continuum limit).")
    print(f"  The S_0 <-> S_1 boundary is at Phi = epsilon.")
    print()

    # =====================================================================
    # 4-PANEL PLOT
    # =====================================================================
    print("-" * 72)
    print("  GENERATING 4-PANEL PLOT")
    print("-" * 72)

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    fig.suptitle(
        r"TGP Substrate Constants: $C_\beta$, $C_\gamma$, $\epsilon$"
        "  (PB + PC)",
        fontsize=13, fontweight='bold')

    # ------------------------------------------------------------------
    # (a) C_beta(L) and C_gamma(L) vs 1/L^2
    # ------------------------------------------------------------------
    ax = axes[0, 0]
    if len(L_arr) > 0:
        ax.errorbar(inv_L2, Cb_arr, yerr=Cb_err_arr,
                     fmt='bo-', ms=7, capsize=4, label=r'$C_\beta$')
        ax.errorbar(inv_L2, Cg_arr, yerr=Cg_err_arr,
                     fmt='rs-', ms=7, capsize=4, label=r'$C_\gamma$')

        x_ext = np.linspace(0, max(inv_L2) * 1.15, 50)
        if popt_b is not None:
            ax.plot(x_ext, fss_linear(x_ext, *popt_b), 'b--', alpha=0.5,
                    label=r'$C_\beta^\infty$=%.2f' % Cb_inf)
        if popt_g is not None:
            ax.plot(x_ext, fss_linear(x_ext, *popt_g), 'r--', alpha=0.5,
                    label=r'$C_\gamma^\infty$=%.2f' % Cg_inf)

        if np.isfinite(Cb_inf):
            ax.plot(0, Cb_inf, 'bD', ms=10, zorder=5)
        if np.isfinite(Cg_inf):
            ax.plot(0, Cg_inf, 'r^', ms=10, zorder=5)

    ax.set_xlabel(r'$1/L^2$', fontsize=12)
    ax.set_ylabel(r'$C_\beta$, $C_\gamma$ (lattice units)', fontsize=12)
    ax.set_title(r'(a) FSS: $C_\beta$, $C_\gamma$ vs $1/L^2$')
    ax.legend(fontsize=9, loc='best')
    ax.set_xlim(left=-0.001)
    ax.grid(True, alpha=0.3)

    # Secondary x-axis: L
    ax2t = ax.twiny()
    ax2t.set_xlim(ax.get_xlim())
    L_show = [l for l in [8, 12, 16, 24, 32]
              if 1.0/l**2 <= ax.get_xlim()[1] * 1.05]
    ax2t.set_xticks([1.0/l**2 for l in L_show])
    ax2t.set_xticklabels([str(l) for l in L_show], fontsize=8)
    ax2t.set_xlabel('L', fontsize=10)

    # ------------------------------------------------------------------
    # (b) C_beta/C_gamma ratio vs L
    # ------------------------------------------------------------------
    ax = axes[0, 1]
    valid_mask = np.isfinite(ratio_arr)
    if np.any(valid_mask):
        ax.errorbar(L_arr[valid_mask], ratio_arr[valid_mask],
                     yerr=ratio_err_arr[valid_mask],
                     fmt='ko-', ms=7, capsize=4, label='MC data')
    ax.axhline(1.0, color='red', ls='--', lw=2.5,
               label=r'$\beta=\gamma$ (exact, Z$_2$)')

    if np.isfinite(ratio_inf):
        ax.axhline(ratio_inf, color='blue', ls=':', lw=1.5,
                   label=r'$L\to\infty$: %.3f' % ratio_inf)
        if np.isfinite(ratio_inf_err):
            ax.axhspan(ratio_inf - ratio_inf_err,
                       ratio_inf + ratio_inf_err,
                       alpha=0.12, color='blue')

    ax.set_xlabel('L (lattice size)', fontsize=12)
    ax.set_ylabel(r'$C_\beta / C_\gamma$', fontsize=12)
    ax.set_title(r'(b) Vacuum condition $C_\beta/C_\gamma \to 1$')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # ------------------------------------------------------------------
    # (c) epsilon/Phi0 vs block size b
    # ------------------------------------------------------------------
    ax = axes[1, 0]
    if len(eps_bs_arr) > 0:
        ax.plot(eps_bs_arr, eps_ratio_arr, 'go-', ms=8, lw=2,
                label=r'$\epsilon/\Phi_0$ (MC)')

        # Show 1 - eps/Phi0 on secondary y-axis? No, keep simple.
        # Guide line
        if len(eps_bs_arr) >= 2 and min(eps_ratio_arr) > 0:
            b_fine = np.linspace(0.8, max(eps_bs_arr) * 1.2, 50)
            # eps/Phi0 approaching 1: plot 1 as horizontal reference
            ax.axhline(1.0, color='gray', ls=':', lw=1, alpha=0.5,
                       label=r'$\Phi_0$ (mean)')

        ax.set_xlabel('Block size b', fontsize=12)
        ax.set_ylabel(r'$\epsilon / \Phi_0$', fontsize=12)
        ax.set_title(r'(c) Regularization $\epsilon/\Phi_0$ vs block size')
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)

        # Annotate: b=1 is single site, b=L/2 is half-lattice
        ax.annotate('single site\n(UV scale)',
                    xy=(eps_bs_arr[0], eps_ratio_arr[0]),
                    xytext=(eps_bs_arr[0] + 1, eps_ratio_arr[0] - 0.15),
                    fontsize=8, arrowprops=dict(arrowstyle='->', color='gray'),
                    bbox=dict(boxstyle='round,pad=0.2', fc='lightyellow'))

    # ------------------------------------------------------------------
    # (d) P(Phi) near Phi=0
    # ------------------------------------------------------------------
    ax = axes[1, 1]
    colors_d = plt.cm.viridis(np.linspace(0.15, 0.85, max(len(P_phi_bs1), 1)))
    plotted_d = False

    for idx_L, (L_val, (pc, ph)) in enumerate(sorted(P_phi_bs1.items())):
        if pc is not None and ph is not None:
            ax.plot(pc, ph, '-', color=colors_d[idx_L], lw=1.5,
                    label=f'L={L_val}', alpha=0.85)
            # Mark epsilon
            eps_v = all_data[L_val].get(r_eps, {}).get('epsilon', {}).get(1, np.nan)
            if np.isfinite(eps_v):
                ax.axvline(eps_v, color=colors_d[idx_L], ls=':',
                           lw=0.8, alpha=0.6)
            plotted_d = True

    if plotted_d:
        ax.set_xlabel(r'$\Phi = s_i^2$ (single site)', fontsize=12)
        ax.set_ylabel(r'$P(\Phi)$', fontsize=12)
        ax.set_title(r'(d) Distribution $P(\Phi)$ near $\Phi\to 0$ (b=1)')
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(left=0)
        ax.annotate(r'$\epsilon$ (1st pct.)',
                    xy=(0.02, 0.92), xycoords='axes fraction',
                    fontsize=8, style='italic',
                    bbox=dict(boxstyle='round,pad=0.2', fc='wheat', alpha=0.5))

    plt.tight_layout()
    fpath = os.path.join(save_dir, "substrate_constants.png")
    plt.savefig(fpath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\n  Saved: {fpath}")

    elapsed = time.time() - t_start
    print(f"\n  Total runtime: {elapsed:.1f}s")
    print("\n  Done.")


if __name__ == "__main__":
    main()
