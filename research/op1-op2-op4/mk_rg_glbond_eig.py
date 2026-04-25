#!/usr/bin/env python3
"""
M7 (Track B sketch / P3.2 of external_review_2026-04-25/review_response_plan.md):
J_GL eigenvalue at the M3 fixed point.

Question:  is the v2 GL bond  H_GL = J_GL Sum_<ij> A_ij Phi_i Phi_j (Phi_j-Phi_i)^2
RELEVANT or IRRELEVANT under MK-RG linearised at the 3D Ising WF FP?

If irrelevant (|lambda_GL| < 1): bare J_GL decays under MK; the GL bond
cannot drive closure of OP-2b on its own.  M6's closure-in-principle at
J_GL ~ 5.89 was an artefact of unphysically large bare J_GL.

If relevant (|lambda_GL| > 1): bare J_GL grows under MK; even a small
bare value can drive closure in the IR; M6's finding stands as evidence
that GL is the missing physics.

Method (M7_glbond_eigenvalue.md):

  1. Find M3 FP cb_M3*.
  2. Compute K_new at cb_M3* via one M3 RG step (M3 baseline).
  3. Build a 2D outer (s_1, s_3) grid; for each pair, compute the
     source-shifted moments M_2(h), M_4(h), M_6(h) at h = K_eff*(s_1+s_3)
     under exp(-V(s_2) + h s_2) using the M3 single-site quadrature.
  4. Form
       dF/J_GL_pre = -4 * [(s1^2+s3^2)*M_6(h)
                            - 2(s1^4+s3^4)*M_4(h)
                            + (s1^6+s3^6)*M_2(h)]
     (this is the linear response of the 2D F(s_1, s_3) to J_GL_pre,
     with the factor 4 = b^{d-1} = J_GL_eff/J_GL_pre absorbed.)
  5. Project onto monomials s_1^a s_3^b (a, b >= 0 even, a+b <= 12);
     read off:
       lambda_pre = -[coeff at (2,6)]
       lambda_GL  = lambda_pre / K_new
     plus the consistency check: coeff(2,6) ?= coeff(6,2) (bilateral
     symmetry), and the shape diagnostic c(4,4)/c(2,6) (expected -2
     if pure-GL-shaped).
  6. Numerical-perturbation cross-check: also compute lambda_pre by
     finite difference (F(eps) - F(0)) / eps for small eps, and compare.

Decision criteria (M7_glbond_eigenvalue.md sec. 6):
  |lambda_GL| > 1:  RELEVANT;  GL bond is the missing physics candidate.
  0.5 < |lambda_GL| < 1: IRRELEVANT but slow.
  |lambda_GL| < 0.5: STRONGLY IRRELEVANT; M3 + GL ~ M3 in IR.
"""

import numpy as np
from math import factorial

from mk_rg_bgamma import MigdalKadanoffRGN, B_RESCALE, BD1, D


# =============================================================================
# Core 2D-extraction routines
# =============================================================================

def compute_M3_FP(n_ops=8, n_quad=1200, s_max=10.0,
                  r0=-2.251, u0=3.917, tol=1e-10, max_iter=8000):
    """Find the M3 fixed point at given truncation."""
    rg = MigdalKadanoffRGN(n_ops=n_ops, n_quad=n_quad, s_max=s_max)
    cb0 = np.zeros(n_ops, dtype=float)
    cb0[0] = r0; cb0[1] = u0
    cb_fp, it, conv = rg.find_fp_bar(cb0, tol=tol, max_iter=max_iter)
    return rg, cb_fp, it, conv


def compute_K_new_M3(rg, cb_FP, K_pre=1.0):
    """One M3 RG step at the FP -> returns K_new."""
    _, K_new = rg.rg_step(cb_FP, K_pre)
    return K_new


def make_outer_grid(s_max, n_outer):
    """Outer 1D Gauss-Legendre grid on [-s_max, s_max]."""
    nodes, weights = np.polynomial.legendre.leggauss(n_outer)
    s_outer = s_max * nodes
    ds_outer = s_max * weights
    return s_outer, ds_outer


def compute_F_components_2d(rg, cb_FP, K_pre, s_outer, ds_outer,
                            J_GL_pre=0.0):
    """
    Compute the full F(s_1, s_3) = ln int ds_2 exp[...] on the 2D outer
    grid. Used both at J_GL_pre = 0 (M3 baseline) and at finite J_GL_pre
    (numerical perturbation cross-check).

    Returns F2d of shape (n_outer, n_outer).
    """
    K_eff = BD1 * K_pre
    J_GL_eff = BD1 * J_GL_pre
    n_outer = len(s_outer)
    s2g = rg.s_pts
    s2g_sq = s2g * s2g

    # Single-site V on inner grid.
    log_w0 = np.zeros_like(s2g)
    s2k = np.ones_like(s2g_sq)
    for k, c in enumerate(cb_FP):
        s2k = s2k * s2g_sq
        log_w0 = log_w0 - (c / (2.0 * (k + 1))) * s2k

    # 2D outer grid.
    s1 = s_outer[:, None]               # (n_outer, 1)
    s3 = s_outer[None, :]               # (1, n_outer)
    s1_sq = s1 * s1; s1_4 = s1_sq * s1_sq; s1_6 = s1_4 * s1_sq
    s3_sq = s3 * s3; s3_4 = s3_sq * s3_sq; s3_6 = s3_4 * s3_sq

    h2d = K_eff * (s1 + s3)             # (n_outer, n_outer)
    h_flat = h2d.flatten()              # (n_outer^2,)

    # Source shift + GL bond contributions to log-integrand.
    # log_w_full[i_pair, j_inner] = log_w0[j] + h*s_2[j]
    #     - J_GL_eff * (bond_12 + bond_23)
    log_wh = log_w0[None, :] + h_flat[:, None] * s2g[None, :]

    if J_GL_eff != 0.0:
        # bond_12 + bond_23 at each (s_1, s_3) and inner s_2:
        # = (s1^2 + s3^2) s2^6 - 2(s1^4 + s3^4) s2^4 + (s1^6 + s3^6) s2^2
        s2_p2 = s2g_sq
        s2_p4 = s2_p2 * s2_p2
        s2_p6 = s2_p4 * s2_p2
        sum_sq = (s1_sq + s3_sq).flatten()    # (n_outer^2,)
        sum_4 = (s1_4 + s3_4).flatten()
        sum_6 = (s1_6 + s3_6).flatten()
        bond_term = (sum_sq[:, None] * s2_p6[None, :]
                     - 2.0 * sum_4[:, None] * s2_p4[None, :]
                     + sum_6[:, None] * s2_p2[None, :])
        log_wh = log_wh - J_GL_eff * bond_term

    # Stable log-sum-exp per row.
    max_row = np.max(log_wh, axis=1, keepdims=True)
    wh = rg.ds[None, :] * np.exp(log_wh - max_row)
    Zh = np.sum(wh, axis=1)
    Zh = np.where(Zh > 1e-300, Zh, 1e-300)
    F = np.log(Zh) + max_row[:, 0]      # F(h) at this (s_1, s_3) pair
    F2d = F.reshape(n_outer, n_outer)
    return F2d


def compute_dF_GL_first_order(rg, cb_FP, K_pre, s_outer):
    """
    First-order-in-J_GL_pre linear response of F(s_1, s_3):
        dF[i, j] / J_GL_pre
            = -4 * [(s_1^2+s_3^2) M_6(h)
                    - 2(s_1^4+s_3^4) M_4(h)
                    + (s_1^6+s_3^6) M_2(h)],
        h = K_eff*(s_1+s_3),
        K_eff = b^{d-1} K_pre.

    The factor -4 absorbs J_GL_eff = 4 J_GL_pre.

    Returns (dF, M2, M4, M6) where M_n are reshaped to (n_outer, n_outer).
    """
    K_eff = BD1 * K_pre
    n_outer = len(s_outer)
    s2g = rg.s_pts
    s2g_sq = s2g * s2g

    # Single-site V on inner grid (no J_GL_eff perturbation; we are at
    # J_GL_pre = 0 in this computation, the linear response).
    log_w0 = np.zeros_like(s2g)
    s2k = np.ones_like(s2g_sq)
    for k, c in enumerate(cb_FP):
        s2k = s2k * s2g_sq
        log_w0 = log_w0 - (c / (2.0 * (k + 1))) * s2k

    # 2D outer grid.
    s1 = s_outer[:, None]
    s3 = s_outer[None, :]
    s1_sq = s1 * s1; s1_4 = s1_sq * s1_sq; s1_6 = s1_4 * s1_sq
    s3_sq = s3 * s3; s3_4 = s3_sq * s3_sq; s3_6 = s3_4 * s3_sq
    h2d = K_eff * (s1 + s3)
    h_flat = h2d.flatten()

    # Source-shifted weights.
    log_wh = log_w0[None, :] + h_flat[:, None] * s2g[None, :]
    max_row = np.max(log_wh, axis=1, keepdims=True)
    wh = rg.ds[None, :] * np.exp(log_wh - max_row)
    Zh = np.sum(wh, axis=1)
    Zh = np.where(Zh > 1e-300, Zh, 1e-300)

    # Moments M_2, M_4, M_6 at each (s_1, s_3).
    s2_p2 = s2g_sq
    s2_p4 = s2_p2 * s2_p2
    s2_p6 = s2_p4 * s2_p2
    M2 = (np.sum(wh * s2_p2[None, :], axis=1) / Zh).reshape(n_outer, n_outer)
    M4 = (np.sum(wh * s2_p4[None, :], axis=1) / Zh).reshape(n_outer, n_outer)
    M6 = (np.sum(wh * s2_p6[None, :], axis=1) / Zh).reshape(n_outer, n_outer)

    # dF / J_GL_pre = -4 * <bond_12 + bond_23>.
    bond_avg = ((s1_sq + s3_sq) * M6
                - 2.0 * (s1_4 + s3_4) * M4
                + (s1_6 + s3_6) * M2)
    dF = -4.0 * bond_avg
    return dF, M2, M4, M6


def project_2d_monomials(F2d, s_outer, ds_outer, n_max=12, parity="even"):
    """
    L^2-weighted lstsq projection of F(s_1, s_3) onto monomials
    s_1^a s_3^b with a+b <= n_max.

    parity = "even"  : a, b both even >= 0  (the GL/operator-mixing sector)
    parity = "all"   : a, b >= 0 of either parity
    parity = "even_outer" : a, b both even, AND a + b is even (always true if both even)

    Returns dict {(a, b): coeff} and (basis_pairs, design_matrix, coeffs_array).
    """
    n_outer = len(s_outer)
    if parity == "even":
        pairs = [(a, b)
                 for a in range(0, n_max + 1, 2)
                 for b in range(0, n_max + 1, 2)
                 if a + b <= n_max]
    elif parity == "all":
        pairs = [(a, b)
                 for a in range(0, n_max + 1)
                 for b in range(0, n_max + 1)
                 if a + b <= n_max]
    else:
        raise ValueError(f"Unknown parity={parity}")

    n_pairs = len(pairs)
    W2d = np.outer(np.sqrt(np.abs(ds_outer)),
                   np.sqrt(np.abs(ds_outer)))
    W_flat = W2d.flatten()
    F_flat = F2d.flatten()
    A = np.empty((n_outer * n_outer, n_pairs))
    for idx, (a, b) in enumerate(pairs):
        col = (s_outer ** a)[:, None] * (s_outer ** b)[None, :]
        A[:, idx] = col.flatten()

    A_w = A * W_flat[:, None]
    F_w = F_flat * W_flat
    coeffs, *_ = np.linalg.lstsq(A_w, F_w, rcond=None)
    return {pairs[i]: coeffs[i] for i in range(n_pairs)}, pairs, A, coeffs


def project_onto_GL_direction(F2d, s_outer, ds_outer):
    """
    Direct L^2 projection onto O_GL = s_1^2 s_3^6 - 2 s_1^4 s_3^4 + s_1^6 s_3^2:
        c = <F, O_GL> / <O_GL, O_GL>.

    The L^2 inner product is 2D Gauss-Legendre quadrature on
    [-s_max, s_max]^2.
    """
    s1 = s_outer[:, None]; s3 = s_outer[None, :]
    s1_sq = s1 * s1; s1_4 = s1_sq * s1_sq; s1_6 = s1_4 * s1_sq
    s3_sq = s3 * s3; s3_4 = s3_sq * s3_sq; s3_6 = s3_4 * s3_sq
    O_GL = s1_sq * s3_6 - 2.0 * s1_4 * s3_4 + s1_6 * s3_sq
    W2d = np.outer(ds_outer, ds_outer)        # signed weights are fine (Gauss-Legendre).
    inner_F = np.sum(F2d * O_GL * W2d)
    inner_GL = np.sum(O_GL * O_GL * W2d)
    return inner_F / inner_GL, inner_F, inner_GL


# =============================================================================
# Main eigenvalue computation
# =============================================================================

def compute_J_GL_eigenvalue(n_ops=8, n_quad=1200, s_max=10.0,
                            n_outer=30, eps_pert=1e-3, n_max_proj=12,
                            verbose=True):
    """
    Pipeline of M7 sec.7. Returns a dict of all diagnostics.
    """
    if verbose:
        print("M7 -- J_GL eigenvalue at M3 fixed point")
        print("=" * 78)

    # --- Step 1: M3 FP ---
    if verbose:
        print()
        print("=== (1) M3 fixed point (n_ops={}, n_quad={}, s_max={}) "
              "===".format(n_ops, n_quad, s_max))
    rg, cb_FP, it, conv = compute_M3_FP(n_ops=n_ops, n_quad=n_quad,
                                        s_max=s_max)
    if verbose:
        print(f"  converged: {conv} in {it} iterations")
        print(f"  cb_FP[:4]: r*={cb_FP[0]:+.5f}  u*={cb_FP[1]:+.5f}  "
              f"B*={cb_FP[2]:+.5f}  G*={cb_FP[3]:+.5f}")
    if not conv:
        raise RuntimeError("M3 FP did not converge; eigenvalue undefined.")

    # --- Step 2: K_new at the FP ---
    K_new = compute_K_new_M3(rg, cb_FP, K_pre=1.0)
    if verbose:
        print(f"  K_new at M3 FP = {K_new:+.6f}")

    # --- Step 3-4: outer grid + first-order linear response ---
    if verbose:
        print()
        print(f"=== (2) Linear-response dF/dJ_GL_pre on 2D outer grid "
              f"(n_outer={n_outer}, n_max_proj={n_max_proj}) ===")
    s_outer, ds_outer = make_outer_grid(s_max=s_max, n_outer=n_outer)
    dF, M2, M4, M6 = compute_dF_GL_first_order(rg, cb_FP, K_pre=1.0,
                                               s_outer=s_outer)

    # --- Step 5: project dF onto monomials ---
    coeffs, pairs, A, c_arr = project_2d_monomials(
        dF, s_outer, ds_outer, n_max=n_max_proj, parity="even")

    c_2_6 = coeffs.get((2, 6), 0.0)
    c_6_2 = coeffs.get((6, 2), 0.0)
    c_4_4 = coeffs.get((4, 4), 0.0)
    c_8_4 = coeffs.get((8, 4), 0.0)
    c_4_8 = coeffs.get((4, 8), 0.0)
    c_2_4 = coeffs.get((2, 4), 0.0)
    c_4_2 = coeffs.get((4, 2), 0.0)
    c_6_6 = coeffs.get((6, 6), 0.0)
    c_2_2 = coeffs.get((2, 2), 0.0)
    if verbose:
        print(f"  Projected coefficients of dF/dJ_GL_pre on monomials s_1^a s_3^b:")
        print(f"    (a, b) = (2, 2):  {c_2_2:+.6e}")
        print(f"    (a, b) = (2, 4):  {c_2_4:+.6e}    [(4, 2) = {c_4_2:+.6e}]")
        print(f"    (a, b) = (2, 6):  {c_2_6:+.6e}    [(6, 2) = {c_6_2:+.6e}]")
        print(f"    (a, b) = (4, 4):  {c_4_4:+.6e}")
        if (4, 8) in coeffs:
            print(f"    (a, b) = (4, 8):  {c_4_8:+.6e}    [(8, 4) = {c_8_4:+.6e}]")
            print(f"    (a, b) = (6, 6):  {c_6_6:+.6e}")

    # Bilateral symmetry check.
    sym_err = (abs(c_2_6 - c_6_2) /
               max(abs(c_2_6), abs(c_6_2), 1e-20))
    if verbose:
        print(f"  bilateral symmetry  |c(2,6) - c(6,2)| / max = {sym_err:.3e}")

    # GL shape diagnostic: c(4,4) / c(2,6) -- expected -2 for pure GL.
    if abs(c_2_6) > 1e-20:
        shape_44 = c_4_4 / c_2_6
        if verbose:
            print(f"  shape c(4,4)/c(2,6) = {shape_44:+.4f}   "
                  f"(GL bond expects -2.0000)")
    else:
        shape_44 = float("nan")

    # --- lambda_pre from (2,6) extraction; lambda_GL = lambda_pre / K_new ---
    lambda_pre_26 = -c_2_6
    lambda_GL_26 = lambda_pre_26 / K_new
    lambda_pre_62 = -c_6_2
    lambda_GL_62 = lambda_pre_62 / K_new
    lambda_pre_avg = -0.5 * (c_2_6 + c_6_2)
    lambda_GL_avg = lambda_pre_avg / K_new
    if verbose:
        print()
        print(f"  lambda_pre = -c(2,6) = {lambda_pre_26:+.6e}")
        print(f"  lambda_pre = -c(6,2) = {lambda_pre_62:+.6e}   "
              "(bilateral check)")
        print(f"  lambda_pre (avg)     = {lambda_pre_avg:+.6e}")
        print(f"  lambda_GL = lambda_pre / K_new = "
              f"{lambda_pre_avg:+.6e} / {K_new:+.6f} = {lambda_GL_avg:+.6e}")

    # --- L^2 projection onto O_GL ---
    if verbose:
        print()
        print("=== (3) L^2 projection onto O_GL = s_1^2 s_3^6 - 2 s_1^4 s_3^4 + s_1^6 s_3^2 ===")
    proj_GL, ip_F, ip_GL = project_onto_GL_direction(
        dF, s_outer, ds_outer)
    if verbose:
        print(f"  <dF, O_GL>            = {ip_F:+.6e}")
        print(f"  <O_GL, O_GL>          = {ip_GL:+.6e}")
        print(f"  proj  = <dF, O_GL>/<O_GL,O_GL> = {proj_GL:+.6e}")
    lambda_pre_norm = -proj_GL
    lambda_GL_norm = lambda_pre_norm / K_new
    if verbose:
        print(f"  lambda_pre_norm = -proj  = {lambda_pre_norm:+.6e}")
        print(f"  lambda_GL_norm  = lambda_pre_norm / K_new = "
              f"{lambda_GL_norm:+.6e}")

    # --- Step 6: numerical-perturbation cross-check ---
    if verbose:
        print()
        print(f"=== (4) Numerical-perturbation cross-check (eps={eps_pert}) ===")
    F0 = compute_F_components_2d(rg, cb_FP, K_pre=1.0, s_outer=s_outer,
                                 ds_outer=ds_outer, J_GL_pre=0.0)
    F_eps = compute_F_components_2d(rg, cb_FP, K_pre=1.0, s_outer=s_outer,
                                    ds_outer=ds_outer, J_GL_pre=eps_pert)
    F_neg = compute_F_components_2d(rg, cb_FP, K_pre=1.0, s_outer=s_outer,
                                    ds_outer=ds_outer, J_GL_pre=-eps_pert)
    dF_num_sym = (F_eps - F_neg) / (2.0 * eps_pert)
    coeffs_num, _, _, _ = project_2d_monomials(
        dF_num_sym, s_outer, ds_outer, n_max=n_max_proj, parity="even")
    c_2_6_num = coeffs_num.get((2, 6), 0.0)
    c_6_2_num = coeffs_num.get((6, 2), 0.0)
    c_4_4_num = coeffs_num.get((4, 4), 0.0)
    lambda_pre_num = -0.5 * (c_2_6_num + c_6_2_num)
    lambda_GL_num = lambda_pre_num / K_new
    if verbose:
        print(f"  numerical-perturbation lambda_pre = "
              f"{lambda_pre_num:+.6e}    (analytic: {lambda_pre_avg:+.6e})")
        rel_err = (abs(lambda_pre_num - lambda_pre_avg)
                   / max(abs(lambda_pre_avg), 1e-20))
        print(f"  relative error = {rel_err:.3e}")
        print(f"  numerical-perturbation lambda_GL  = {lambda_GL_num:+.6e}")
        if abs(c_2_6_num) > 1e-20:
            print(f"  shape c(4,4)/c(2,6) [num]: {c_4_4_num/c_2_6_num:+.4f}  "
                  f"(expected -2 if pure GL)")

    # --- Step 7: canonical baseline ---
    lambda_GL_can = BD1 / (K_new ** 4)
    if verbose:
        print()
        print("=== (5) Canonical-scaling baseline ===")
        print(f"  Canonical lambda_GL = b^(d-1) / K_new^4 = "
              f"{BD1} / {K_new:.6f}^4 = {lambda_GL_can:+.6e}")
        print(f"  Numerical lambda_GL = {lambda_GL_avg:+.6e}")
        print(f"  ratio (numerical / canonical) = "
              f"{lambda_GL_avg / lambda_GL_can:+.6f}  "
              f"(if ~ 1, canonical estimate ~ correct)")

    # --- Verdict ---
    abs_lam = abs(lambda_GL_avg)
    if abs_lam > 1.0:
        verdict = "RELEVANT"
        msg = ("J_GL grows under MK; small bare J_GL can flow to "
               "finite IR value. M6 closure-in-principle is supported.")
    elif abs_lam > 0.5:
        verdict = "IRRELEVANT (slow)"
        msg = ("J_GL decays under MK, but slowly. Finite-J_GL "
               "transients could persist across modest RG depth.")
    else:
        verdict = "STRONGLY IRRELEVANT"
        msg = ("J_GL decays fast; M3 + GL ~ M3 in the IR; M6 "
               "closure-in-principle was an artefact of unphysically "
               "large bare J_GL.")
    if verbose:
        print()
        print("=== (6) Verdict ===")
        print(f"  |lambda_GL| = {abs_lam:.5f}  -> {verdict}")
        print(f"  -> {msg}")

    return dict(
        cb_FP=cb_FP, K_new=K_new,
        lambda_pre_26=lambda_pre_26,
        lambda_pre_62=lambda_pre_62,
        lambda_pre_avg=lambda_pre_avg,
        lambda_GL_26=lambda_GL_26,
        lambda_GL_62=lambda_GL_62,
        lambda_GL_avg=lambda_GL_avg,
        lambda_GL_norm=lambda_GL_norm,
        lambda_GL_num=lambda_GL_num,
        lambda_GL_can=lambda_GL_can,
        c_2_6=c_2_6, c_6_2=c_6_2, c_4_4=c_4_4,
        sym_err=sym_err, shape_44=shape_44,
        verdict=verdict,
    )


def main():
    res = compute_J_GL_eigenvalue(n_ops=8, n_quad=1200, s_max=10.0,
                                  n_outer=30, eps_pert=1e-3,
                                  n_max_proj=12, verbose=True)

    # Robustness: also run at finer outer grid + larger n_max_proj.
    print()
    print("=" * 78)
    print("=== ROBUSTNESS scan: vary n_outer in {20, 30, 40, 50}, "
          "n_max_proj in {8, 12, 16} ===")
    print("=" * 78)
    print()
    print(f"  (n_outer, n_max_proj) lambda_GL_avg    lambda_GL_norm    "
          "lambda_GL_num    sym_err     shape(4,4)/(2,6)")
    print(f"  --------------------- --------------   --------------    "
          "--------------   ---------   --------")
    for n_o in (20, 30, 40, 50):
        for n_mp in (8, 12, 16):
            try:
                rr = compute_J_GL_eigenvalue(
                    n_ops=8, n_quad=1200, s_max=10.0,
                    n_outer=n_o, eps_pert=1e-3, n_max_proj=n_mp,
                    verbose=False)
                print(f"  ({n_o:>3d}, {n_mp:>2d})             "
                      f"{rr['lambda_GL_avg']:+.6e}   "
                      f"{rr['lambda_GL_norm']:+.6e}    "
                      f"{rr['lambda_GL_num']:+.6e}   "
                      f"{rr['sym_err']:.2e}   "
                      f"{rr['shape_44']:+.4f}")
            except Exception as e:
                print(f"  ({n_o:>3d}, {n_mp:>2d})             FAILED: {e}")


if __name__ == "__main__":
    main()
