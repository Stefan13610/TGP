#!/usr/bin/env python3
"""
monte_carlo_substrate.py -- Theory of Generated Space (TGP)
============================================================
Monte Carlo simulation of the TGP substrate Hamiltonian on a 3D cubic
lattice.  Validates effective parameters (beta_eff, gamma_eff) by
Metropolis sampling of the static Hamiltonian:

    H = sum_i [m0^2/2 * s_i^2 + lambda0/4 * s_i^4] - J sum_<ij> s_i*s_j

Z2 symmetry s_i -> -s_i breaks spontaneously at T < T_c, giving v != 0.
The spatiality field Phi = <s^2> is block-averaged; its effective potential
V_eff(Phi) is fitted to extract beta_eff, gamma_eff.

Wilson-Fisher fixed point (MK-RG): r* = -2.251, u* = 3.917, |r*|/u* = 0.575

Outputs (to scripts/plots/):
    mc_phase_diagram.png       -- v(r) and chi(r) for different L
    mc_effective_potential.png  -- V_eff(Phi) at several r values
    mc_beta_gamma_ratio.png    -- beta_eff/gamma_eff vs r

Usage:
    python monte_carlo_substrate.py              # default run
    python monte_carlo_substrate.py --production # long production run
"""

import os
import sys
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ===========================================================================
# 3D Phi^4 Lattice Model
# ===========================================================================

class Phi4Lattice3D:
    """Scalar phi^4 on 3D cubic lattice with periodic BC."""

    def __init__(self, L, m0_sq, lam0, J, beta_T):
        self.L, self.N = L, L ** 3
        self.m0_sq, self.lam0, self.J = m0_sq, lam0, J
        self.beta_T = beta_T
        self.s = 0.1 * np.random.randn(L, L, L)

    def _nn_sum(self, s):
        """Sum of nearest neighbors for every site (vectorized)."""
        return (np.roll(s, 1, 0) + np.roll(s, -1, 0)
              + np.roll(s, 1, 1) + np.roll(s, -1, 1)
              + np.roll(s, 1, 2) + np.roll(s, -1, 2))

    def total_energy(self):
        """Total Hamiltonian H."""
        s = self.s
        e_onsite = np.sum(0.5 * self.m0_sq * s**2 + 0.25 * self.lam0 * s**4)
        e_nn = -0.5 * self.J * np.sum(s * self._nn_sum(s))
        return e_onsite + e_nn

    def vectorized_sweep(self, delta=1.0):
        """Checkerboard Metropolis: update even/odd sublattices in parallel."""
        L = self.L
        idx = np.arange(L)
        ix, iy, iz = np.meshgrid(idx, idx, idx, indexing='ij')
        accepted = 0
        for parity in (0, 1):
            mask = ((ix + iy + iz) % 2 == parity)
            si_old = self.s.copy()
            nn_sum = self._nn_sum(si_old)
            si_new = si_old + delta * (np.random.random((L, L, L)) * 2 - 1)
            dH = (0.5 * self.m0_sq * (si_new**2 - si_old**2)
                + 0.25 * self.lam0 * (si_new**4 - si_old**4)
                - self.J * (si_new - si_old) * nn_sum)
            acc = ((dH <= 0) | (np.random.random((L, L, L))
                   < np.exp(-self.beta_T * np.minimum(dH, 20.0)))) & mask
            self.s = np.where(acc, si_new, si_old)
            accepted += np.sum(acc)
        return accepted / self.N

    def measure(self):
        """Return dict of basic observables."""
        s = self.s
        m = np.mean(s)
        return {'mag': m, 'mag_abs': np.abs(m),
                's2': np.mean(s**2), 's4': np.mean(s**4), 'mag2': m**2}

    def block_s2(self, block_size=2):
        """Block-average s^2 into blocks of given linear size."""
        L, bs = self.L, block_size
        if L % bs != 0:
            return None
        nb = L // bs
        return (self.s ** 2).reshape(nb, bs, nb, bs, nb, bs).mean(axis=(1, 3, 5))


# ===========================================================================
# Simulation driver
# ===========================================================================

def run_simulation(L, m0_sq, lam0, J, beta_T,
                   n_therm, n_measure, n_skip):
    """Run MC and return measurements dict."""
    lat = Phi4Lattice3D(L, m0_sq, lam0, J, beta_T)
    delta = 1.5

    # Thermalization with adaptive step-size
    for i in range(n_therm):
        acc = lat.vectorized_sweep(delta)
        if i < n_therm // 2 and i % 50 == 49:
            if acc > 0.55:
                delta *= 1.1
            elif acc < 0.35:
                delta *= 0.9

    # Measurement phase
    obs_list, phi_blocks = [], []
    for _ in range(n_measure):
        for __ in range(n_skip):
            lat.vectorized_sweep(delta)
        obs_list.append(lat.measure())
        bl = lat.block_s2(block_size=2)
        if bl is not None:
            phi_blocks.extend(bl.ravel().tolist())

    # Aggregate observables
    avg = {}
    for k in obs_list[0]:
        vals = np.array([o[k] for o in obs_list])
        avg[k] = np.mean(vals)
        avg[k + '_err'] = np.std(vals) / np.sqrt(len(vals))

    mag_vals = np.array([o['mag'] for o in obs_list])
    avg['v'] = np.mean(np.abs(mag_vals))
    avg['v_err'] = np.std(np.abs(mag_vals)) / np.sqrt(len(mag_vals))
    avg['chi'] = L**3 * (np.mean(mag_vals**2) - np.mean(np.abs(mag_vals))**2)
    avg['phi_blocks'] = np.array(phi_blocks)
    return avg


def extract_veff(phi_samples, T, n_bins=40):
    """Extract V_eff(Phi) = -T * ln P(Phi) from histogram of block-avg s^2."""
    if len(phi_samples) < 100:
        return None, None
    hist, edges = np.histogram(phi_samples, bins=n_bins, density=True)
    centers = 0.5 * (edges[:-1] + edges[1:])
    good = hist > 0
    veff = -T * np.log(hist[good])
    veff -= np.min(veff)
    return centers[good], veff


def fit_veff_coefficients(phi_arr, veff_arr):
    """Fit V_eff to polynomial; extract beta_eff (cubic) and gamma_eff (quartic).

    The TGP effective potential form is:
        V_eff ~ (beta_eff/3)*Phi^3 - (gamma_eff/4)*Phi^4
    So from polyfit coefficients: beta_eff = 3*c3, gamma_eff = -4*c4.
    """
    if phi_arr is None or len(phi_arr) < 6:
        return np.nan, np.nan
    try:
        i_min = np.argmin(veff_arr)
        dp = phi_arr - phi_arr[i_min]
        dv = veff_arr - veff_arr[i_min]
        c4, c3, c2, c1, c0 = np.polyfit(dp, dv, 4)
        return 3.0 * c3, -4.0 * c4
    except Exception:
        return np.nan, np.nan


# ===========================================================================
# Parameter scan
# ===========================================================================

def parameter_scan(Ls, r_values, u, beta_T, n_therm, n_measure, n_skip):
    """Scan over r = m0^2/J at fixed u = lambda0/J and beta_T."""
    J, lam0 = 1.0, u
    results = {}
    total, done = len(Ls) * len(r_values), 0

    for L in Ls:
        res = {'r': [], 'v': [], 'v_err': [], 'chi': [], 's2': [],
               'beta_eff': [], 'gamma_eff': [],
               'phi_centers': [], 'veff_curves': []}
        results[L] = res
        for r in r_values:
            done += 1
            print(f"  [{done}/{total}] L={L}, r={r:.2f}", end="", flush=True)
            avg = run_simulation(L, r * J, lam0, J, beta_T,
                                 n_therm, n_measure, n_skip)
            phi_c, veff = extract_veff(avg['phi_blocks'], 1.0 / beta_T)
            b_eff, g_eff = fit_veff_coefficients(phi_c, veff)

            res['r'].append(r)
            res['v'].append(avg['v'])
            res['v_err'].append(avg['v_err'])
            res['chi'].append(avg['chi'])
            res['s2'].append(avg['s2'])
            res['beta_eff'].append(b_eff)
            res['gamma_eff'].append(g_eff)
            res['phi_centers'].append(phi_c)
            res['veff_curves'].append(veff)

            if g_eff != 0 and np.isfinite(g_eff):
                print(f"  v={avg['v']:.4f}, chi={avg['chi']:.2f}, "
                      f"b/g={b_eff/g_eff:.3f}")
            else:
                print(f"  v={avg['v']:.4f}, chi={avg['chi']:.2f}")
    return results


# ===========================================================================
# Plotting
# ===========================================================================

def plot_phase_diagram(results, plot_dir):
    """Plot v(r) and chi(r) for different L."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    for L in sorted(results):
        d = results[L]
        r, v, v_err = np.array(d['r']), np.array(d['v']), np.array(d['v_err'])
        chi = np.array(d['chi'])
        ax1.errorbar(r, v, yerr=v_err, marker='o', ms=4,
                     label=f'L={L}', capsize=2)
        ax2.plot(r, chi, marker='s', ms=4, label=f'L={L}')
    ax1.set_xlabel('r = m0^2 / J')
    ax1.set_ylabel('v = <|<s>|>')
    ax1.set_title('Order parameter vs bare mass')
    ax1.legend()
    ax1.axhline(0, color='gray', ls='--', lw=0.5)
    ax1.axvline(-2.251, color='red', ls=':', lw=1, label='WF r*')
    ax2.set_xlabel('r = m0^2 / J')
    ax2.set_ylabel('chi (susceptibility)')
    ax2.set_title('Susceptibility vs bare mass')
    ax2.legend()
    plt.tight_layout()
    path = os.path.join(plot_dir, 'mc_phase_diagram.png')
    plt.savefig(path, dpi=150)
    plt.close()
    print(f"  Saved {path}")


def plot_effective_potential(results, plot_dir):
    """Plot V_eff(Phi) at several r values for the largest L."""
    L_max = max(results)
    d = results[L_max]
    fig, ax = plt.subplots(figsize=(8, 6))
    n_r = len(d['r'])
    for idx in np.linspace(0, n_r - 1, min(5, n_r), dtype=int):
        phi_c, veff, r_val = d['phi_centers'][idx], d['veff_curves'][idx], d['r'][idx]
        if phi_c is not None and len(phi_c) > 0:
            ax.plot(phi_c, veff, marker='.', ms=3, label=f'r={r_val:.2f}')
    ax.set_xlabel('Phi = block-avg <s^2>')
    ax.set_ylabel('V_eff(Phi) (shifted)')
    ax.set_title(f'Effective potential V_eff(Phi), L={L_max}')
    ax.legend()
    plt.tight_layout()
    path = os.path.join(plot_dir, 'mc_effective_potential.png')
    plt.savefig(path, dpi=150)
    plt.close()
    print(f"  Saved {path}")


def plot_beta_gamma_ratio(results, plot_dir):
    """Plot beta_eff/gamma_eff vs r."""
    fig, ax = plt.subplots(figsize=(8, 5))
    for L in sorted(results):
        d = results[L]
        r = np.array(d['r'])
        b, g = np.array(d['beta_eff']), np.array(d['gamma_eff'])
        good = np.isfinite(b) & np.isfinite(g) & (np.abs(g) > 1e-10)
        if np.any(good):
            ax.plot(r[good], b[good] / g[good], marker='o', ms=4, label=f'L={L}')
    ax.axhline(1.0, color='red', ls='--', lw=1.5, label='vacuum condition')
    ax.set_xlabel('r = m0^2 / J')
    ax.set_ylabel('beta_eff / gamma_eff')
    ax.set_title('Vacuum condition: beta_eff / gamma_eff vs r')
    ax.legend()
    ax.set_ylim(-5, 5)
    plt.tight_layout()
    path = os.path.join(plot_dir, 'mc_beta_gamma_ratio.png')
    plt.savefig(path, dpi=150)
    plt.close()
    print(f"  Saved {path}")


# ===========================================================================
# Main
# ===========================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Monte Carlo simulation of TGP substrate')
    parser.add_argument('--production', action='store_true',
                        help='Long production run (L=16, 50k measurements)')
    parser.add_argument('--L', type=int, default=None,
                        help='Override lattice size')
    parser.add_argument('--n_therm', type=int, default=None)
    parser.add_argument('--n_measure', type=int, default=None)
    parser.add_argument('--n_skip', type=int, default=None)
    args = parser.parse_args()

    # Parameter sets
    if args.production:
        Ls = [8, 12, 16]
        n_therm, n_measure, n_skip = 20000, 50000, 10
    else:
        Ls = [8, 12, 16]
        n_therm, n_measure, n_skip = 5000, 10000, 5

    if args.L is not None:
        Ls = [args.L]
    if args.n_therm is not None:
        n_therm = args.n_therm
    if args.n_measure is not None:
        n_measure = args.n_measure
    if args.n_skip is not None:
        n_skip = args.n_skip

    # Physics parameters
    J = 1.0
    u = 4.0       # lambda0/J, near Wilson-Fisher u*=3.917
    beta_T = 0.5  # inverse temperature (middle of 0.1--1.0 range)
    r_values = np.linspace(-3.0, 0.0, 13)

    # Output directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    plot_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'plots')
    os.makedirs(plot_dir, exist_ok=True)

    print("=" * 60)
    print("TGP Substrate Monte Carlo Simulation")
    print("=" * 60)
    print(f"Lattice sizes: {Ls}")
    print(f"r = m0^2/J scan: {r_values[0]:.1f} to {r_values[-1]:.1f} "
          f"({len(r_values)} points)")
    print(f"u = lambda0/J = {u:.1f},  beta_T = 1/T = {beta_T:.2f}")
    print(f"Therm: {n_therm}, Meas: {n_measure}, Skip: {n_skip}")
    print("-" * 60)

    results = parameter_scan(Ls, r_values, u, beta_T,
                             n_therm, n_measure, n_skip)

    print("\nGenerating plots...")
    plot_phase_diagram(results, plot_dir)
    plot_effective_potential(results, plot_dir)
    plot_beta_gamma_ratio(results, plot_dir)

    # Summary
    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    L_ref = max(Ls)
    d = results[L_ref]
    r_arr, v_arr = np.array(d['r']), np.array(d['v'])
    chi_arr = np.array(d['chi'])
    b_arr, g_arr = np.array(d['beta_eff']), np.array(d['gamma_eff'])

    i_tc = np.argmax(chi_arr)
    print(f"  Phase transition (max chi) near r = {r_arr[i_tc]:.2f} (L={L_ref})")
    print(f"  v at transition: {v_arr[i_tc]:.4f}")
    print(f"  Wilson-Fisher prediction: r* = -2.251, |r*|/u* = 0.575")

    good = np.isfinite(b_arr) & np.isfinite(g_arr) & (np.abs(g_arr) > 1e-10)
    if np.any(good):
        ratio = b_arr[good] / g_arr[good]
        r_good = r_arr[good]
        i_vac = np.argmin(np.abs(ratio - 1.0))
        print(f"  Vacuum condition (beta/gamma~1) best at r = "
              f"{r_good[i_vac]:.2f}, ratio = {ratio[i_vac]:.3f}")

    print("\nDone.")


if __name__ == '__main__':
    main()
