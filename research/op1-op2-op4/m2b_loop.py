#!/usr/bin/env python3
"""
M2b one-loop V_eff(Phi)
=======================

Numerical computation of the one-loop effective potential for the v2
substrate Hamiltonian (see M2b_loop_derivation.md):

    H_Gamma = sum_i V_onsite(Phi_i) + J sum_{<ij>} Phi_i Phi_j (Phi_j - Phi_i)^2

with
    V_onsite(Phi) = (m0^2/2) Phi + (lambda0/4) Phi^2 + (T/2) ln Phi

and composite field Phi = s^2.  Background-field one-loop formula
(eq. (2.7) of M2b_loop_derivation.md):

    V_eff(Phi) = V_onsite(Phi) + (T/2) * Integral_BZ [d^3q/(2pi)^3]
                    * ln[ K_Phi(q; Phi) + A(Phi) ]

with
    K_Phi(q; Phi) = 2 J Phi^2 * (3 - cos q_x - cos q_y - cos q_z)      (a = 1)
    A(Phi)        = V_onsite''(Phi) = lambda0/2 - T/(2 Phi^2)

UV cutoff: first Brillouin zone q_mu in [-pi, pi]; by symmetry, integrate
over [0, pi]^3 with factor 8 and divide by (2 pi)^3 = (pi)^3 normalisation.

Physical parameter regime: same three cases as M2a (m0^2, lambda0, T) in
ordered-phase, real-saddle regime (m0^4 > 4 lambda0 T).  J = 1.

Output:
  - Tree-level check (c_2, c_3, c_4 from direct polynomial fit of
    V_onsite(Phi) vs analytical M2a formulas).
  - One-loop shifts delta c_2, delta c_3, delta c_4.
  - Full (tree + 1-loop) Taylor coefficients and beta_eff/gamma_eff.
  - Dodatek-B C_gamma estimator: C_gamma = gamma_eff * v^2 * a^2 / lambda0^2
    at MF v^2 = |m0^2|/lambda0.
"""

import numpy as np
from pathlib import Path


def v_onsite(phi, m0sq, lambda0, T):
    """Onsite potential V_onsite(Phi) including H-S Jacobian."""
    return 0.5*m0sq*phi + 0.25*lambda0*phi**2 + 0.5*T*np.log(phi)


def A_phi(phi, lambda0, T):
    """A(Phi) = V_onsite''(Phi) = lambda0/2 - T/(2 Phi^2)."""
    return 0.5*lambda0 - 0.5*T/phi**2


def phi_saddle(m0sq, lambda0, T):
    """Saddle of V_onsite': lambda0 Phi_0^2 + m0^2 Phi_0 + T = 0."""
    disc = m0sq**2 - 4*lambda0*T
    if disc < 0:
        raise ValueError(f"No real saddle: disc = {disc}")
    return (-m0sq + np.sqrt(disc))/(2*lambda0)


def make_bz_grid(n_q):
    """Gauss-Legendre 3D grid on [0, pi]^3 (8-th of the BZ).
    Returns:
        cos_sum (N,): cos q_x + cos q_y + cos q_z at each grid point.
        weights (N,): associated weights (product of 1D weights).
    The BZ integral is then
        int_BZ d^3q/(2pi)^3 f(q) = (1/pi^3) * sum_i w_i * f(q_i)
    using the symmetry of the integrand (even in each q_mu).
    """
    nodes_1d, w_1d = np.polynomial.legendre.leggauss(n_q)
    # Map [-1, 1] -> [0, pi]:
    q_1d = 0.5*np.pi*(nodes_1d + 1.0)
    dq_1d = 0.5*np.pi*w_1d
    qx, qy, qz = np.meshgrid(q_1d, q_1d, q_1d, indexing='ij')
    wx, wy, wz = np.meshgrid(dq_1d, dq_1d, dq_1d, indexing='ij')
    weights = (wx*wy*wz).flatten()
    cos_sum = (np.cos(qx) + np.cos(qy) + np.cos(qz)).flatten()
    return cos_sum, weights


def v_loop_one(phi, m0sq, lambda0, T, J, cos_sum, weights):
    """One-loop correction (T/2) int_BZ d^3q/(2pi)^3 ln[K_Phi + A]
    evaluated at a single Phi."""
    A = A_phi(phi, lambda0, T)
    K = 2.0 * J * phi**2 * (3.0 - cos_sum)
    arg = K + A
    if np.any(arg <= 0):
        min_arg = arg.min()
        raise ValueError(f"K + A <= 0 somewhere: min={min_arg:+.4e} at phi={phi}, A={A:+.4e}")
    integral = np.sum(weights * np.log(arg)) / np.pi**3
    return 0.5 * T * integral


def v_eff_grid(phi_grid, m0sq, lambda0, T, J, cos_sum, weights):
    """Compute V_tree, V_1loop_correction, V_full on the Phi grid."""
    v_tree = v_onsite(phi_grid, m0sq, lambda0, T)
    v_1loop = np.array([v_loop_one(p, m0sq, lambda0, T, J, cos_sum, weights)
                        for p in phi_grid])
    return v_tree + v_1loop, v_tree, v_1loop


def fit_taylor(phi_grid, v_arr, phi0, window=0.12):
    """Fit a degree-4 polynomial to v(Phi) around Phi_0 and return
    (a1, a2, a3, a4) = coefficients of (Phi - Phi_0)^{1..4}."""
    phi_lo = phi0 * (1 - window)
    phi_hi = phi0 * (1 + window)
    mask = (phi_grid >= phi_lo) & (phi_grid <= phi_hi)
    if mask.sum() < 6:
        raise RuntimeError(f"Too few fit points: {mask.sum()}")
    x = phi_grid[mask] - phi0
    y = v_arr[mask]
    coeffs = np.polyfit(x, y, 4)
    a4, a3, a2, a1, a0 = coeffs
    return a1, a2, a3, a4


def run_case(label, m0sq, lambda0, T, J=1.0, n_q=32, n_phi=400, window=0.12):
    print(f"\n=== case {label}: m0^2={m0sq}, lambda0={lambda0}, T={T}, J={J} ===")
    phi0 = phi_saddle(m0sq, lambda0, T)
    A0 = A_phi(phi0, lambda0, T)
    print(f"  Phi_0 (tree saddle):            {phi0:.6f}")
    print(f"  A(Phi_0) = V''_onsite(Phi_0):   {A0:+.6e}  "
          f"(loop stable: {'YES' if A0 > 0 else 'NO -- tachyonic'})")
    if A0 <= 0:
        print("  Skipping case (A(Phi_0) <= 0).")
        return

    # Tree-level analytical coefficients:
    c2_t = (lambda0*phi0**2 - T) / (4.0 * phi0**2)
    c3_t = T / (6.0 * phi0**3)
    c4_t = -T / (8.0 * phi0**4)
    beta_t = 3*c3_t
    gamma_t = -4*c4_t
    print(f"  TREE analytical:  c_2 = {c2_t:+.6e}, c_3 = {c3_t:+.6e}, c_4 = {c4_t:+.6e}")
    print(f"                    beta_eff = {beta_t:+.6e},  gamma_eff = {gamma_t:+.6e}")

    # Build BZ integration grid and Phi grid:
    cos_sum, weights = make_bz_grid(n_q)
    phi_grid = np.linspace(0.35*phi0, 2.5*phi0, n_phi)
    # Valid region: A(Phi) > 0 => Phi > sqrt(T/lambda0). Otherwise tachyonic.
    phi_thresh = np.sqrt(T/lambda0)
    phi_grid = phi_grid[phi_grid > phi_thresh*1.001]
    if len(phi_grid) < 50:
        print(f"  Warning: small Phi grid ({len(phi_grid)} pts).")
    print(f"  BZ: n_q={n_q}^3 = {n_q**3} nodes;  Phi-grid: {len(phi_grid)} pts "
          f"in [{phi_grid[0]:.4f}, {phi_grid[-1]:.4f}];  A > 0 threshold = {phi_thresh:.4f}")

    v_full, v_tree_arr, v_1loop_arr = v_eff_grid(phi_grid, m0sq, lambda0, T, J, cos_sum, weights)

    # Tree-only fit (sanity check vs analytical):
    a1_t, a2_t, a3_t, a4_t = fit_taylor(phi_grid, v_tree_arr, phi0, window)
    print(f"  Tree FIT vs analytical:")
    print(f"    c_2 fit = {a2_t:+.6e}  (an {c2_t:+.6e}, relerr {(a2_t/c2_t - 1):+.2e})")
    print(f"    c_3 fit = {a3_t:+.6e}  (an {c3_t:+.6e}, relerr {(a3_t/c3_t - 1):+.2e})")
    print(f"    c_4 fit = {a4_t:+.6e}  (an {c4_t:+.6e}, relerr {(a4_t/c4_t - 1):+.2e})")

    # 1-loop correction (fit V_1loop alone to isolate shift):
    a1_L, a2_L, a3_L, a4_L = fit_taylor(phi_grid, v_1loop_arr, phi0, window)
    print(f"  Loop-only FIT (coefficient shifts):")
    print(f"    dc_2 = {a2_L:+.6e}")
    print(f"    dc_3 = {a3_L:+.6e}")
    print(f"    dc_4 = {a4_L:+.6e}")

    # Full (tree + 1-loop) fit:
    a1_f, a2_f, a3_f, a4_f = fit_taylor(phi_grid, v_full, phi0, window)
    beta_1l = 3*a3_f
    gamma_1l = -4*a4_f
    print(f"  1-LOOP TOTAL:")
    print(f"    c_2 = {a2_f:+.6e}  (tree {c2_t:+.6e})")
    print(f"    c_3 = {a3_f:+.6e}  (tree {c3_t:+.6e})")
    print(f"    c_4 = {a4_f:+.6e}  (tree {c4_t:+.6e})")
    print(f"    beta_eff_1loop   = {beta_1l:+.6e}   (tree {beta_t:+.6e},  ratio {beta_1l/beta_t:+.4f})")
    print(f"    gamma_eff_1loop  = {gamma_1l:+.6e}   (tree {gamma_t:+.6e}, ratio {gamma_1l/gamma_t:+.4f})")
    # Dimensional ratio:
    print(f"    beta/gamma (dim, raw ratio)     = {beta_1l/gamma_1l:+.6f}  "
          f"(Phi_0 = {phi0:.4f}; tree: beta/gamma = Phi_0)")
    # Dimensionless ratio:
    print(f"    beta_natural / gamma_natural    = {(beta_1l*phi0**3)/(gamma_1l*phi0**4):+.6f}  "
          f"(tree = 1.0)")

    # Dodatek B matching: gamma_eff = (lambda0^2 / (v^2 a^2)) * C_gamma,
    # with MF v^2 = |m0^2|/lambda0 (a = 1).
    v_sq_mf = abs(m0sq)/lambda0 if m0sq < 0 else np.nan
    if np.isfinite(v_sq_mf):
        C_gamma_tree = gamma_t * v_sq_mf / lambda0**2
        C_gamma_1l = gamma_1l * v_sq_mf / lambda0**2
        C_beta_tree = beta_t * v_sq_mf / lambda0  # (dodatek B: beta = lambda0/a^2 * C_beta)
        C_beta_1l = beta_1l * v_sq_mf / lambda0
        print(f"  Dodatek B matching  (v^2 = |m0^2|/lambda0 = {v_sq_mf:.4f}):")
        print(f"    C_gamma (tree)   = {C_gamma_tree:+.4e}")
        print(f"    C_gamma (1-loop) = {C_gamma_1l:+.4e}   "
              f"(shift: {(C_gamma_1l - C_gamma_tree)/C_gamma_tree:+.3f})")
        print(f"    C_beta  (tree)   = {C_beta_tree:+.4e}")
        print(f"    C_beta  (1-loop) = {C_beta_1l:+.4e}")
        print(f"    C_beta/C_gamma (1-loop) = {C_beta_1l/C_gamma_1l:+.4f}  "
              f"(M3 gives -0.57; M3 vs M2b are different regimes)")


def main():
    print("=" * 72)
    print("M2b one-loop V_eff(Phi) -- numerical correction to M2a tree-level")
    print("v2 bond: K_Phi(q) = 2 J Phi^2 sum (1 - cos q_mu),  Lambda = pi,  a = 1")
    print("=" * 72)
    cases = [
        ("A", -4.0, 1.0, 1.0),
        ("B", -5.0, 2.0, 1.0),
        ("C", -2.0, 1.0, 0.1),
    ]
    for lab, m0sq, l0, T in cases:
        run_case(lab, m0sq, l0, T, J=1.0, n_q=32, n_phi=400)

    # Convergence check at one case (A), varying n_q:
    print("\n" + "=" * 72)
    print("Convergence scan of n_q (BZ 1D Gauss-Legendre nodes) for case A")
    print("=" * 72)
    for nq in (8, 16, 24, 32, 40):
        print(f"\n-- n_q = {nq} ({nq**3} total) --")
        run_case(f"A(n_q={nq})", -4.0, 1.0, 1.0, J=1.0, n_q=nq, n_phi=400)


if __name__ == "__main__":
    main()
