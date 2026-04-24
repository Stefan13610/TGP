#!/usr/bin/env python3
"""
M2a sanity check
================

Numerically test the tree-level H-S prediction derived analytically in
  M2a_HS_derivation.md
for the composite-field effective potential V_eff(Phi) in the single-site
(on-site only) substrate Hamiltonian

    H_loc(s) = (m0^2 / 2) s^2 + (lambda_0 / 4) s^4 ,

with Phi := s^2.  The analytical tree-level prediction is

    V_eff(Phi) = (m0^2/2) Phi + (lambda_0/4) Phi^2 + (T/2) ln(Phi) + const.

Taylor coefficients at the saddle Phi_0:

    c_2 = (lambda_0 Phi_0^2 - T) / (4 Phi_0^2)
    c_3 = T / (6 Phi_0^3)
    c_4 = -T / (8 Phi_0^4)

Dodatek B identifications:

    beta_eff  = 3 c_3 = T / (2 Phi_0^3)
    gamma_eff = -4 c_4 = T / (2 Phi_0^4)
    beta_eff / gamma_eff = Phi_0

In dimensionless phi = Phi / Phi_0 (TGP natural units):

    beta = beta_eff * Phi_0^3  = T/2
    gamma = gamma_eff * Phi_0^4 = T/2
    beta / gamma = 1            <-- Route 1 + Route 2 unified.

This script:
  (a) computes V_eff(Phi) by direct numerical integration of the single-site
      Boltzmann distribution over s,
  (b) fits a polynomial of Phi around Phi_0, extracts (c_2, c_3, c_4),
  (c) compares against the analytical formulas above,
  (d) verifies beta_eff / gamma_eff = Phi_0.

v2-status (post 2026-04-24 axiom pivot): this file contains no claim about
the full H_Gamma (v2 GL bond); it treats the on-site (m0^2, lambda_0) part
only.  Closing the GL bond contribution to V_eff is M2b work.
"""

import numpy as np
from pathlib import Path


def v_eff_analytical(phi, m0sq, lambda0, T):
    """Tree-level V_eff(Phi) from H-S derivation (eq. 2.5 of M2a)."""
    return 0.5 * m0sq * phi + 0.25 * lambda0 * phi**2 + 0.5 * T * np.log(phi)


def phi_saddle_analytical(m0sq, lambda0, T):
    """Saddle of v_eff: lambda0 Phi_0^2 + m0sq Phi_0 + T = 0, take + branch."""
    disc = m0sq**2 - 4 * lambda0 * T
    if disc < 0:
        raise ValueError(f"No real saddle: m0^4 < 4*lambda0*T (disc={disc})")
    return (-m0sq + np.sqrt(disc)) / (2 * lambda0)


def compute_v_eff_numerical(m0sq, lambda0, T, phi_grid):
    """Compute V_eff(Phi) = -T ln P(Phi) by direct integration of P(s).

    For each Phi in the grid, P(Phi) = P(s=+sqrt(Phi)) / |ds/dPhi|
                                     + P(s=-sqrt(Phi)) / |ds/dPhi|
    = 2 * P(s=sqrt(Phi)) * (1 / (2 sqrt(Phi)))  (by symmetry)
    = P(s=sqrt(Phi)) / sqrt(Phi)

    where P(s) ~ exp(-beta H_loc(s)), beta = 1/T.
    """
    beta = 1.0 / T
    s_pos = np.sqrt(phi_grid)
    h_loc = 0.5 * m0sq * phi_grid + 0.25 * lambda0 * phi_grid**2
    # P(Phi) up to normalisation:
    # P(Phi) propto Phi^{-1/2} * exp(-beta h_loc(Phi))
    log_p = -0.5 * np.log(phi_grid) - beta * h_loc
    v_eff = -T * log_p
    # Drop the overall additive const (will cancel in Taylor coefficients):
    v_eff -= v_eff.min()
    return v_eff


def compute_v_eff_mc(m0sq, lambda0, T, n_samples=2_000_000, n_bins=300):
    """MC sample s from single-site Boltzmann, histogram Phi = s^2,
    return (phi_centers, v_eff_mc = -T log P(Phi))."""
    # Metropolis sampling of s in 1D.
    rng = np.random.default_rng(42)
    beta = 1.0 / T
    s = 1.0
    samples = np.empty(n_samples)
    step = 0.5  # proposal width, tuned by hand for ~50% acceptance
    accepted = 0
    burn = 10_000
    for t in range(-burn, n_samples):
        s_new = s + step * rng.standard_normal()
        h_old = 0.5 * m0sq * s**2 + 0.25 * lambda0 * s**4
        h_new = 0.5 * m0sq * s_new**2 + 0.25 * lambda0 * s_new**4
        if rng.uniform() < np.exp(-beta * (h_new - h_old)):
            s = s_new
            accepted += 1
        if t >= 0:
            samples[t] = s
    acc_rate = accepted / (n_samples + burn)
    # Histogram Phi = s^2:
    phi_samples = samples**2
    # Clip to (0, reasonable max) to avoid tail instability in log:
    phi_max = np.percentile(phi_samples, 99.5)
    phi_min = np.percentile(phi_samples, 0.5)
    phi_edges = np.linspace(phi_min, phi_max, n_bins + 1)
    hist, _ = np.histogram(phi_samples, bins=phi_edges, density=True)
    phi_centers = 0.5 * (phi_edges[1:] + phi_edges[:-1])
    # Avoid log(0) by dropping empty bins:
    mask = hist > 0
    phi_centers = phi_centers[mask]
    p_hist = hist[mask]
    v_eff_mc = -T * np.log(p_hist)
    v_eff_mc -= v_eff_mc.min()
    return phi_centers, v_eff_mc, acc_rate


def fit_taylor_coefficients(phi_centers, v_eff, phi_0, window=0.15):
    """Fit v_eff(Phi) as a polynomial a0 + a1(Phi-Phi_0) + ... + a4(Phi-Phi_0)^4
    in a window of [phi_0*(1-window), phi_0*(1+window)].  Return (a1, a2, a3, a4)."""
    phi_lo = phi_0 * (1 - window)
    phi_hi = phi_0 * (1 + window)
    mask = (phi_centers >= phi_lo) & (phi_centers <= phi_hi)
    if mask.sum() < 6:
        raise RuntimeError(f"Too few points in fit window [{phi_lo}, {phi_hi}] "
                           f"(only {mask.sum()} points; need >=6)")
    x = phi_centers[mask] - phi_0
    y = v_eff[mask]
    # Polynomial fit of degree 4.  numpy.polyfit returns highest degree first.
    coeffs = np.polyfit(x, y, 4)
    # Rearrange to (a4, a3, a2, a1, a0) order then split:
    a4, a3, a2, a1, a0 = coeffs
    return a1, a2, a3, a4


def run_case(m0sq, lambda0, T, do_mc=True, label=""):
    print(f"=== case {label}  (m0^2={m0sq}, lambda0={lambda0}, T={T}) ===")

    # Analytical saddle and coefficients:
    phi0 = phi_saddle_analytical(m0sq, lambda0, T)
    c2_a = (lambda0 * phi0**2 - T) / (4.0 * phi0**2)
    c3_a = T / (6.0 * phi0**3)
    c4_a = -T / (8.0 * phi0**4)
    beta_eff_a = 3.0 * c3_a
    gamma_eff_a = -4.0 * c4_a

    print(f"  Phi_0 (analytical):        {phi0:.6f}")
    print(f"  c_2 (analytical):          {c2_a:+.6e}")
    print(f"  c_3 (analytical):          {c3_a:+.6e}")
    print(f"  c_4 (analytical):          {c4_a:+.6e}")
    print(f"  beta_eff (analytical):     {beta_eff_a:+.6e}")
    print(f"  gamma_eff (analytical):    {gamma_eff_a:+.6e}")
    print(f"  beta_eff/gamma_eff:        {beta_eff_a/gamma_eff_a:.6f}  (expect Phi_0 = {phi0:.6f})")
    print(f"  In natural units (phi=Phi/Phi_0):")
    print(f"    beta  = beta_eff  * Phi_0^3 = {beta_eff_a * phi0**3:+.6f}  (expect T/2 = {T/2:+.6f})")
    print(f"    gamma = gamma_eff * Phi_0^4 = {gamma_eff_a * phi0**4:+.6f}  (expect T/2 = {T/2:+.6f})")
    print(f"    beta/gamma = {(beta_eff_a*phi0**3)/(gamma_eff_a*phi0**4):.6f}  (expect 1.0)")

    # (a) Direct numerical V_eff from analytical marginal:
    phi_grid = np.linspace(0.05 * phi0, 3.0 * phi0, 2000)
    v_eff_direct = compute_v_eff_numerical(m0sq, lambda0, T, phi_grid)
    a1d, a2d, a3d, a4d = fit_taylor_coefficients(phi_grid, v_eff_direct, phi0)
    print(f"  [direct integration]  c_2 = {a2d:+.6e}   c_3 = {a3d:+.6e}   c_4 = {a4d:+.6e}")
    print(f"                        beta_eff = {3*a3d:+.6e}   gamma_eff = {-4*a4d:+.6e}")
    print(f"                        beta_eff/gamma_eff = {(3*a3d)/(-4*a4d):.6f}  "
          f"(rel.err. vs Phi_0: {(3*a3d)/(-4*a4d)/phi0 - 1:+.2e})")

    # (b) MC sampling:
    if do_mc:
        phi_mc, v_eff_mc, acc = compute_v_eff_mc(m0sq, lambda0, T)
        a1m, a2m, a3m, a4m = fit_taylor_coefficients(phi_mc, v_eff_mc, phi0)
        print(f"  [MC, acc={acc:.2f}]   c_2 = {a2m:+.6e}   c_3 = {a3m:+.6e}   c_4 = {a4m:+.6e}")
        print(f"                        beta_eff = {3*a3m:+.6e}   gamma_eff = {-4*a4m:+.6e}")
        print(f"                        beta_eff/gamma_eff = {(3*a3m)/(-4*a4m):.6f}  "
              f"(rel.err. vs Phi_0: {(3*a3m)/(-4*a4m)/phi0 - 1:+.2e})")
    print()


def main():
    # Three parameter cases spanning the (m0^2, lambda0, T) triangle.
    # Use units where J = 1 (dodatek B convention: r = m0^2/J, u = lambda0/J).
    # Constraint for real saddle: m0^4 > 4*lambda0*T  (otherwise thermal melting).
    cases = [
        ("A",  -4.0,  1.0, 1.0),   # ordered phase, moderate T,   m0^4=16 > 4
        ("B",  -5.0,  2.0, 1.0),   # deeper quartic, same T,       m0^4=25 > 8
        ("C",  -2.0,  1.0, 0.1),   # colder, further from critical, m0^4=4 > 0.4
    ]
    print("M2a tree-level H-S sanity check -- testing beta_eff/gamma_eff = Phi_0")
    print("(analytical vs. direct numerical integration vs. MC sampling)")
    print()
    for label, m0sq, lambda0, T in cases:
        run_case(m0sq, lambda0, T, do_mc=True, label=label)


if __name__ == "__main__":
    main()
