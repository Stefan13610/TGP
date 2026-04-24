#!/usr/bin/env python3
"""
M3: N-operator MK-RG for computing C_beta / C_gamma at the 3D Ising
Wilson-Fisher fixed point.

Extends tooling/scripts/substrate/renormalization_substrate.py
(which tracks r, u, K) to include all even-power operators up to
s^{2 N_ops}.

State: (c_2 = r, c_4 = u, c_6 = B, c_8 = Gamma, c_10 = E, c_12 = L, ..., K).
Single-site potential:
  V(s) = sum_{k=1}^{N_ops} c_{2k} / (2k)  s^{2k}
       = r/2 s^2 + u/4 s^4 + B/6 s^6 + Gamma/8 s^8 + E/10 s^10 + L/12 s^12 + ...

Decimation integral (h = sum of surviving neighbours):
  F(h) = ln int ds exp[-V(s) + K_eff * h * s]
       = F_0 + sum_{k>=1} F_{2k}/(2k)! h^{2k}
with F_{2k} = K_eff^{2k} * kappa_{2k},   K_eff = b^{d-1} K = 4 K.

After bond-move + decimation with TWO surviving neighbours:
  V'(s) = V(s) - 2 [F(s) - F_0].

Matching s^{2k} / (2k) coefficients:
  c_{2k}' = c_{2k} - 2 (2k) F_{2k} / (2k)! = c_{2k} - 2 F_{2k} / (2k-1)!
and K' = F_2.

Explicit for k=1..6:
  r'  = r - 2 F_2              (2·F_2/1!)
  u'  = u - F_4/3              (2·F_4/3!)
  B'  = B - F_6/60             (2·F_6/5!)
  G'  = G - F_8/2520           (2·F_8/7!)
  E'  = E - F_10/181440        (2·F_10/9!)
  L'  = L - F_12/19958400      (2·F_12/11!)

Cumulants from moments m_k = <s^k>_V via the moment-cumulant recursion:
  kappa_n = m_n - sum_{j=1}^{n-1} C(n-1, j-1) kappa_j m_{n-j}.

Validation protocol:
  (1) Gaussian: V = s^2/2 -> kappa_{2k}=0 for k>=2 (cumulant-formula check).
  (2) 2-op sublimit (held fixed): recover r* = -2.251, u* = 3.917.
  (3) N-op full flow for N_ops = 2, 3, 4, 5, 6: check convergence of B*/Gamma*.

Target (dodatek B eq:B-bg-ratio, written in plan's normalisation):
    beta/gamma |_* = (|r*|/u*) * (C_beta / C_gamma) = v*^2 * (B*/Gamma*).
    Vacuum condition beta = gamma  =>  B*/Gamma* = 1/v*^2 ~ 1.74 (at v*^2 = 0.575).
"""

import numpy as np
from math import comb, factorial

D = 3
B_RESCALE = 2
BD1 = B_RESCALE ** (D - 1)   # = 4


class MigdalKadanoffRGN:
    """N-operator MK-RG; state = (c_2, c_4, ..., c_{2 N_ops}, K)."""

    def __init__(self, n_ops, n_quad=800, s_max=8.0):
        if n_ops < 2:
            raise ValueError("Need at least 2 operators (r, u).")
        self.n_ops = n_ops
        self.n_quad = n_quad
        self.s_max = s_max
        nodes, weights = np.polynomial.legendre.leggauss(n_quad)
        self.s_pts = s_max * nodes
        self.ds = s_max * weights
        # Pre-compute factorials 2*F_{2k}/(2k-1)! denominators.
        self._denoms = np.array(
            [factorial(2 * (k + 1) - 1) for k in range(n_ops)],
            dtype=float,
        )

    def _moments(self, couplings):
        """Return moments m_0=1, m_1, ..., m_{2 n_ops} under exp(-V)."""
        s = self.s_pts
        s2 = s * s
        # Build V(s) = sum c_{2k}/(2k) s^{2k}.
        log_w = np.zeros_like(s)
        s2k = np.ones_like(s)
        for k, c in enumerate(couplings):
            s2k = s2k * s2
            log_w = log_w - (c / (2.0 * (k + 1))) * s2k
        log_w = log_w - np.max(log_w)
        w = self.ds * np.exp(log_w)
        Z = np.sum(w)
        n_m = 2 * self.n_ops + 1
        moments = np.empty(n_m, dtype=float)
        moments[0] = 1.0
        s_pow = np.ones_like(s)
        for k in range(1, n_m):
            s_pow = s_pow * s
            moments[k] = np.sum(s_pow * w) / Z
        return moments

    @staticmethod
    def cumulants_from_moments(moments):
        """Moment-cumulant recursion (Wikipedia formula)."""
        N = len(moments) - 1
        k = np.zeros(N + 1, dtype=float)
        for n in range(1, N + 1):
            val = moments[n]
            for j in range(1, n):
                val -= comb(n - 1, j - 1) * k[j] * moments[n - j]
            k[n] = val
        return k

    def rg_step(self, couplings, K):
        couplings = np.asarray(couplings, dtype=float)
        K_eff = BD1 * K
        moments = self._moments(couplings)
        cumulants = self.cumulants_from_moments(moments)
        new_c = np.empty_like(couplings)
        # K-update from F_2.
        F2 = K_eff**2 * cumulants[2]
        K_new = max(F2, 1e-15)
        # Update each coupling: c_{2k}' = c_{2k} - 2 F_{2k}/(2k-1)!
        for k in range(self.n_ops):
            two_k = 2 * (k + 1)
            F = (K_eff ** two_k) * cumulants[two_k]
            new_c[k] = couplings[k] - 2.0 * F / factorial(two_k - 1)
        return new_c, K_new

    def flow_bar(self, couplings_bar_init, n_steps=40, verbose=False):
        """Iterate with K reset to 1 each step (bar-variable flow)."""
        cb = np.asarray(couplings_bar_init, dtype=float).copy()
        traj = [cb.copy()]
        for step in range(n_steps):
            try:
                c_new, K_new = self.rg_step(cb, 1.0)
                cb_new = c_new / K_new
                if not np.all(np.isfinite(cb_new)) or np.max(np.abs(cb_new)) > 1e14:
                    break
                cb = cb_new
                traj.append(cb.copy())
                if verbose:
                    s = "  ".join(f"c{2*(k+1)}={cb[k]:+.4e}" for k in range(self.n_ops))
                    print(f"step {step+1:3d}  {s}")
            except Exception:
                break
        return np.array(traj)

    def find_fp_bar(self, couplings_bar_init, tol=1e-10, max_iter=5000):
        """Fixed-point iteration in bar-variables."""
        cb = np.asarray(couplings_bar_init, dtype=float).copy()
        for it in range(1, max_iter + 1):
            c_new, K_new = self.rg_step(cb, 1.0)
            if K_new < 1e-30:
                return cb, it, False
            cb_new = c_new / K_new
            if not np.all(np.isfinite(cb_new)):
                return cb, it, False
            delta = np.max(np.abs(cb_new - cb))
            cb = cb_new
            if delta < tol:
                return cb, it, True
        return cb, max_iter, False

    def jacobian_bar(self, cb_fp, eps=1e-5):
        """Numerical Jacobian of the bar-RG map at the fixed point."""
        n = self.n_ops
        J = np.zeros((n, n))

        def G(cb):
            c_new, K_new = self.rg_step(cb, 1.0)
            return c_new / K_new

        g0 = G(cb_fp)
        for j in range(n):
            cbp = cb_fp.copy();  cbp[j] += eps
            cbm = cb_fp.copy();  cbm[j] -= eps
            gp = G(cbp)
            gm = G(cbm)
            J[:, j] = (gp - gm) / (2 * eps)
        return J, g0


# =============================================================================
# Convenience for running at a given n_ops
# =============================================================================

def run_one_truncation(n_ops, n_quad=800, s_max=8.0,
                       r0=-2.251, u0=3.917, verbose_flow=False):
    rg = MigdalKadanoffRGN(n_ops=n_ops, n_quad=n_quad, s_max=s_max)
    cb0 = np.zeros(n_ops, dtype=float)
    cb0[0] = r0
    cb0[1] = u0
    cb_fp, it, conv = rg.find_fp_bar(cb0, tol=1e-10, max_iter=8000)
    return rg, cb_fp, it, conv


def validate_gaussian(rg):
    """r=1, all other c=0 -> V=s^2/2. Expect all higher cumulants zero."""
    cb = np.zeros(rg.n_ops);  cb[0] = 1.0
    moments = rg._moments(cb)
    cumulants = rg.cumulants_from_moments(moments)
    return moments, cumulants


def main():
    print("M3 -- N-operator MK-RG for C_beta / C_gamma at 3D Ising WF")
    print("=" * 72)

    # =====================================================================
    # (1) Gaussian validation with largest truncation (6 ops)
    # =====================================================================
    print()
    print("=== (1) Gaussian sanity (V = s^2/2) ===")
    rg_test = MigdalKadanoffRGN(n_ops=6, n_quad=1200, s_max=12.0)
    moms, cums = validate_gaussian(rg_test)
    print("  moments (even only):")
    for k in range(1, 7):
        print(f"    m_{2*k:<2d} = {moms[2*k]:+.8e}   "
              f"(expect {float(factorial(2*k) // (2**k * factorial(k))):.1f})")
    print("  cumulants:")
    for k in range(1, 7):
        exp = 1.0 if k == 1 else 0.0
        print(f"    kappa_{2*k:<2d} = {cums[2*k]:+.3e}   (expect {exp:.1f})")

    # =====================================================================
    # (2) Truncation convergence
    # =====================================================================
    print()
    print("=== (2) Truncation scan ===")
    print("  n_ops  r*         u*         B*          Gamma*       "
          "|r*|/u*   B*/Gamma*    eigenvalues (abs, sorted)")
    print("  -----  ---------  ---------  ----------  ----------   "
          "--------  ----------   ------------------------------")
    last_ratio = None
    last_rb = None
    last_ub = None
    for n_ops in [2, 3, 4, 5, 6, 7, 8]:
        rg, cb_fp, it, conv = run_one_truncation(n_ops, n_quad=1200, s_max=10.0)
        rb = cb_fp[0]
        ub = cb_fp[1]
        Bb = cb_fp[2] if n_ops >= 3 else float("nan")
        Gb = cb_fp[3] if n_ops >= 4 else float("nan")
        v_sq = abs(rb) / ub if ub > 0 else float("nan")
        ratio = Bb / Gb if (n_ops >= 4 and abs(Gb) > 1e-30) else float("nan")

        # Eigenvalues of linearised flow at FP:
        eigs_abs = None
        try:
            J, g0 = rg.jacobian_bar(cb_fp, eps=1e-5)
            eigs = np.linalg.eigvals(J)
            eigs_abs = np.sort(np.abs(eigs))[::-1]
        except Exception:
            eigs_abs = np.array([])

        eigstr = " ".join(f"{e:.4f}" for e in eigs_abs[:min(4, len(eigs_abs))])
        ratio_str = f"{ratio:+.5f}" if np.isfinite(ratio) else "     --    "
        Bb_str = f"{Bb:+.4e}" if np.isfinite(Bb) else "     --    "
        Gb_str = f"{Gb:+.4e}" if np.isfinite(Gb) else "     --    "
        conv_flag = "*" if conv else "!"
        print(f"  {n_ops}{conv_flag}      {rb:+.5f}   {ub:+.5f}   "
              f"{Bb_str}  {Gb_str}   {v_sq:+.5f}   {ratio_str}   {eigstr}")
        last_ratio = ratio
        last_rb = rb
        last_ub = ub

    # =====================================================================
    # (3) Critical exponents from the highest even truncation
    # =====================================================================
    print()
    print("=== (3) Critical exponents at n_ops=8 (even) ===")
    rg6, cb6, it6, conv6 = run_one_truncation(8, n_quad=1200, s_max=10.0)
    J6, g06 = rg6.jacobian_bar(cb6, eps=1e-5)
    eigs = np.linalg.eigvals(J6)
    eigs_sorted = sorted(eigs, key=lambda z: -abs(z))
    ln_b = np.log(B_RESCALE)
    print(f"  Eigenvalues of RG Jacobian (abs, sorted desc):")
    for i, e in enumerate(eigs_sorted):
        mod = abs(e)
        y = np.log(mod) / ln_b if mod > 1e-15 else -np.inf
        relev = "RELEVANT" if mod > 1 else "irrelevant"
        print(f"    lambda_{i} = {e.real:+.5f}{e.imag:+.5f}j   "
              f"|lambda| = {mod:.5f}   y = {y:+.4f}   ({relev})")
    # Thermal eigenvalue (dominant relevant):
    lam_t = abs(eigs_sorted[0])
    if lam_t > 1.0 + 1e-10:
        nu = ln_b / np.log(lam_t)
        y_t = np.log(lam_t) / ln_b
        print(f"  nu_MK  = {nu:.4f}  (3D Ising: 0.630)   mismatch: {abs(nu-0.630)/0.630*100:.1f}%")
        print(f"  y_t    = {y_t:.4f}  (3D Ising: 1.587)")

    # =====================================================================
    # (4) Final ratio at largest truncation
    # =====================================================================
    print()
    print("=== (4) C_beta/C_gamma verdict ===")
    rb = cb6[0]; ub = cb6[1]; Bb = cb6[2]; Gb = cb6[3]
    v_sq = abs(rb) / ub
    ratio = Bb / Gb
    target = 1.0 / v_sq
    print(f"  n_ops=8 converged: {conv6} in {it6} iterations")
    print(f"  r*     = {rb:+.6f}")
    print(f"  u*     = {ub:+.6f}")
    print(f"  B*     = {Bb:+.6e}")
    print(f"  Gamma* = {Gb:+.6e}")
    if cb6.size >= 5:
        print(f"  E*     = {cb6[4]:+.6e}")
    if cb6.size >= 6:
        print(f"  L*     = {cb6[5]:+.6e}")
    print(f"  v*^2 = |r*|/u* = {v_sq:.5f}")
    print(f"  B*/Gamma*      = {ratio:+.5f}   [= C_beta / C_gamma]")
    print(f"  target (beta=gamma): 1/v*^2 = {target:.5f}")
    print(f"  measured/target = {ratio/target:+.5f}")


if __name__ == "__main__":
    main()
