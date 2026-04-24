#!/usr/bin/env python3
"""
M4 (Test A of external_review_2026-04-25/review_response_plan.md):
N-operator MK-RG with the Hubbard-Stratonovich Jacobian explicit
in the s-action.

Extension of mk_rg_bgamma.py. The single-site weight is

    w(s; mu, eps) = (s^2 + eps^2)^{-mu} * exp(-V_poly(s)),

with V_poly(s) = sum_{k>=1} c_{2k}/(2k) s^{2k} as in M3.

Background:
  Phi := s^2 (composite field). The H-S map ds = dPhi / (2 sqrt(Phi))
  introduces a Jacobian Phi^{-1/2} when going from s- to Phi-action.
  Equivalently, the s-action acquires a `(1/2) ln(s^2)` term, i.e.
  mu_HS = 1/2 in the regulated form mu ln(s^2 + eps^2).

Marginality (M4_phi_variable_derivation.md, sec. 3):
  Under bond-move + decimation + bar-rescaling with the bilinear bond
  -K Sum_<ij> s_i s_j, the operator mu ln(s^2+eps^2) is exactly
  marginal (mu' = mu identically). It thus labels a continuous family
  of MK-RG fixed points; we scan mu and report B*(mu)/Gamma*(mu) and
  1/v*^2(mu).

Convergence:
  In the unregulated form (eps=0), the weight (s^2)^{-mu} is integrable
  near s=0 iff 2 mu < 1, i.e. mu < 1/2. The H-S-natural value mu = 1/2
  is the boundary of convergence: a logarithmic divergence at s=0 in
  the symmetric s-form, even though it is finite in the asymmetric
  Phi >= 0 form. We use eps > 0 as a UV regulator and scan
  mu in {0, 0.1, 0.2, 0.3, 0.4}, eps in {0.1, 0.01}.

Decision criteria for OP-2b (M4 sec. 6):
  (1) Positive: B*(mu_*)/Gamma*(mu_*) crosses 1/v*^2(mu_*) for some
      mu_* in (0, 1/2]. The H-S Jacobian is the missing physics.
  (2) Negative: B*/Gamma* stays negative or does not cross 1/v*^2 for
      all mu in [0, 1/2). M3's result is robust, OP-2b open.
  (3) Inconclusive: trend correct but no crossing in the convergent
      regime; extrapolation suggestive but not provable.

Validation:
  mu = 0 (any eps) must reproduce M3 exactly: B*/Gamma* approx -0.57,
  v*^2 approx 0.81 at N_ops = 8.
"""

import numpy as np
from math import comb, factorial

from mk_rg_bgamma import MigdalKadanoffRGN, B_RESCALE, BD1, D


class MigdalKadanoffPhi(MigdalKadanoffRGN):
    """
    MK-RG with H-S Jacobian explicit in s-action.
    Weight: (s^2 + eps^2)^{-mu} * exp(-V_poly(s)).

    mu is exactly marginal under the MK rules (proof in
    M4_phi_variable_derivation.md sec. 3) so we hold it fixed during
    the flow; only the polynomial sector evolves.
    """

    def __init__(self, n_ops, mu=0.0, eps=0.1, n_quad=1200, s_max=10.0):
        super().__init__(n_ops=n_ops, n_quad=n_quad, s_max=s_max)
        self.mu = float(mu)
        self.eps = float(eps)

    def _moments(self, couplings):
        """
        Return moments m_0..m_{2 n_ops} under
        (s^2 + eps^2)^{-mu} exp(-V_poly(s)).

        For mu = 0 this reduces exactly to MigdalKadanoffRGN._moments.
        """
        s = self.s_pts
        s2 = s * s
        # Polynomial part: log_w_poly = -V_poly(s).
        log_w = np.zeros_like(s)
        s2k = np.ones_like(s)
        for k, c in enumerate(couplings):
            s2k = s2k * s2
            log_w = log_w - (c / (2.0 * (k + 1))) * s2k
        # H-S Jacobian (regulated): -mu * ln(s^2 + eps^2).
        if self.mu != 0.0:
            log_w = log_w - self.mu * np.log(s2 + self.eps * self.eps)
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


def find_fp_with_seeds(rg, seeds, tol=1e-10, max_iter=8000):
    """
    Try a list of (r0, u0, ...) seeds and return the first converged FP
    that has u* > 0 (physical WF) and reasonable magnitude (no
    runaway saddle).
    """
    best = None
    for seed in seeds:
        cb0 = np.zeros(rg.n_ops)
        for j, v in enumerate(seed):
            if j < rg.n_ops:
                cb0[j] = v
        cb_fp, it, conv = rg.find_fp_bar(cb0, tol=tol, max_iter=max_iter)
        if (conv and cb_fp[1] > 0 and cb_fp[0] < 0
                and np.all(np.isfinite(cb_fp))
                and np.max(np.abs(cb_fp)) < 1e3):  # exclude runaway saddles
            return cb_fp, it, True
        if best is None:
            best = (cb_fp, it, conv)
    return best


def scan_mu(n_ops=8, mu_list=(0.0, 0.1, 0.2, 0.3, 0.4),
            eps_list=(0.1, 0.01),
            n_quad_for_eps=None,
            verbose=True):
    """
    Scan mu (H-S coefficient) at fixed eps (regulator), using
    continuation seeding: each mu's FP is seeded with the previous
    mu's FP so we track the same branch smoothly.

    For each (mu, eps) find WF FP and report (r*, u*, B*, Gamma*,
    v*^2, B*/Gamma*, 1/v*^2).
    """
    if n_quad_for_eps is None:
        n_quad_for_eps = {0.1: 1200, 0.01: 2400, 0.05: 1600}
    results = []
    # Default fallbacks if continuation fails.
    base_seeds = [
        (-2.45694, 3.01611),  # M3 8-op FP
        (-2.451, 3.000, -5.0, 9.0),
        (-2.251, 3.917),      # M3 2-op WF
    ]
    for eps in eps_list:
        n_quad = n_quad_for_eps.get(eps, 1200)
        prev_seed = None
        for mu in mu_list:
            # Continuation: previous mu's full FP vector (all c_2k) first.
            seeds = []
            if prev_seed is not None:
                seeds.append(prev_seed)
            seeds.extend(base_seeds)
            rg = MigdalKadanoffPhi(n_ops=n_ops, mu=mu, eps=eps,
                                   n_quad=n_quad, s_max=10.0)
            cb_fp, it, conv = find_fp_with_seeds(rg, seeds)
            rb = cb_fp[0]
            ub = cb_fp[1]
            Bb = cb_fp[2] if n_ops >= 3 else float("nan")
            Gb = cb_fp[3] if n_ops >= 4 else float("nan")
            v_sq = abs(rb) / ub if ub > 0 else float("nan")
            ratio = Bb / Gb if (n_ops >= 4 and abs(Gb) > 1e-30) else float("nan")
            target = 1.0 / v_sq if v_sq > 0 else float("nan")
            row = dict(
                mu=mu, eps=eps, n_quad=n_quad, n_ops=n_ops,
                r=rb, u=ub, B=Bb, Gamma=Gb,
                v_sq=v_sq, ratio=ratio, target=target,
                converged=bool(conv), iterations=int(it),
                cb_fp=cb_fp.copy(),
            )
            results.append(row)
            # Update continuation seed only if we got a physical FP.
            if (conv and ub > 0 and rb < 0
                    and np.max(np.abs(cb_fp)) < 1e3):
                prev_seed = tuple(cb_fp)
            if verbose:
                conv_flag = "*" if conv else "!"
                print(f"  mu={mu:.3f}  eps={eps:.3f}  nq={n_quad:4d}  "
                      f"{conv_flag} it={it:4d}  "
                      f"r*={rb:+.5f}  u*={ub:+.5f}  "
                      f"B*={Bb:+.4e}  G*={Gb:+.4e}  "
                      f"v*^2={v_sq:+.5f}  "
                      f"B*/G*={ratio:+.5f}  1/v*^2={target:+.5f}  "
                      f"diff={ratio-target:+.5f}")
    return results


def main():
    print("M4 -- MK-RG with H-S Jacobian (mu ln(s^2+eps^2)) in s-action")
    print("=" * 78)

    # =====================================================================
    # (0) Validation: mu=0 must reproduce M3.
    # =====================================================================
    print()
    print("=== (0) Validation: mu=0 must reproduce M3 ===")
    print("    Expected (M3, n_ops=8):  r*=-2.45694  u*=+3.01611  "
          "B*/Gamma*=-0.5687  v*^2=0.8146")
    rg0 = MigdalKadanoffPhi(n_ops=8, mu=0.0, eps=0.1,
                            n_quad=1200, s_max=10.0)
    cb0 = np.zeros(8); cb0[0] = -2.251; cb0[1] = 3.917
    cb_fp, it, conv = rg0.find_fp_bar(cb0, tol=1e-10, max_iter=8000)
    rb, ub = cb_fp[0], cb_fp[1]
    Bb, Gb = cb_fp[2], cb_fp[3]
    v_sq = abs(rb) / ub
    ratio = Bb / Gb
    print(f"  Got: r*={rb:+.5f}  u*={ub:+.5f}  "
          f"B*={Bb:+.4e}  Gamma*={Gb:+.4e}  "
          f"v*^2={v_sq:.5f}  B*/Gamma*={ratio:+.5f}  "
          f"converged={conv} in {it} steps")
    err_r = abs(rb - (-2.45694))
    err_u = abs(ub - 3.01611)
    err_ratio = abs(ratio - (-0.5687))
    print(f"  |delta r*| = {err_r:.5f}   |delta u*| = {err_u:.5f}   "
          f"|delta ratio| = {err_ratio:.5f}")
    if err_r < 5e-3 and err_u < 5e-3 and err_ratio < 5e-3:
        print("  -> M3 reproduced (mu=0 baseline OK)")
    else:
        print("  -> WARNING: mu=0 does NOT reproduce M3; investigate "
              "before trusting the scan.")

    # =====================================================================
    # (1) mu-scan at n_ops=8 (largest converging M3 truncation).
    # =====================================================================
    print()
    print("=== (1) mu-scan at n_ops=8, eps in {0.1, 0.01} ===")
    print()
    print("  Columns:  r*  u*  B*  Gamma*  v*^2  B*/Gamma*  1/v*^2  "
          "diff = B*/Gamma* - 1/v*^2")
    print()
    # Fine grid with continuation seeding so we track the same FP branch.
    mu_grid = tuple(round(0.025 * i, 3) for i in range(0, 19))  # 0..0.45 step 0.025
    results8 = scan_mu(
        n_ops=8,
        mu_list=mu_grid,
        eps_list=(0.1, 0.05, 0.01),
        verbose=True,
    )

    # =====================================================================
    # (2) Trend table at eps=0.1 (numerically safe regulator).
    # =====================================================================
    print()
    print("=== (2) Trend table eps=0.1, n_ops=8 ===")
    print()
    print("    mu   |   r*       u*        B*        Gamma*    "
          "v*^2     B*/Gamma*  1/v*^2   diff")
    print("    -----+--------------------------------------------"
          "------------------------------------")
    for row in results8:
        if row["eps"] != 0.1:
            continue
        diff = row["ratio"] - row["target"]
        print(f"    {row['mu']:.2f} | "
              f"{row['r']:+.5f}  {row['u']:+.5f}  "
              f"{row['B']:+.4e}  {row['Gamma']:+.4e}  "
              f"{row['v_sq']:+.5f}  {row['ratio']:+.5f}  "
              f"{row['target']:+.5f}  {diff:+.5f}")

    # =====================================================================
    # (3) Lower-truncation cross-check (n_ops=6) to test stability.
    # =====================================================================
    print()
    print("=== (3) n_ops=6 cross-check at eps=0.1 ===")
    results6 = scan_mu(
        n_ops=6,
        mu_list=tuple(round(0.05 * i, 3) for i in range(0, 10)),
        eps_list=(0.1,),
        verbose=True,
    )

    # =====================================================================
    # (4) Verdict.
    # =====================================================================
    print()
    print("=== (4) OP-2b verdict ===")
    print()
    eps_main = 0.1
    rows = [r for r in results8 if r["eps"] == eps_main]
    diffs = [(r["mu"], r["ratio"] - r["target"]) for r in rows]
    sign_changes = []
    for i in range(1, len(diffs)):
        if diffs[i - 1][1] * diffs[i][1] < 0:
            sign_changes.append((diffs[i - 1][0], diffs[i][0]))
    if sign_changes:
        print(f"  POSITIVE: B*/Gamma* - 1/v*^2 changes sign in mu in "
              f"{sign_changes} (eps={eps_main}).")
        print("  -> H-S Jacobian recovers beta=gamma at WF; OP-2b CLOSED.")
    else:
        all_neg = all(d[1] < 0 for d in diffs)
        all_pos = all(d[1] > 0 for d in diffs)
        if all_neg:
            print(f"  NEGATIVE: B*/Gamma* < 1/v*^2 throughout mu in "
                  f"[{diffs[0][0]}, {diffs[-1][0]}] (eps={eps_main}).")
            print("  -> M3 result is robust under H-S Jacobian addition; "
                  "OP-2b OPEN.")
            print("  -> Missing physics must come from another channel "
                  "(Z_Phi, GL bond, NPRG).")
        elif all_pos:
            print(f"  Above target throughout: B*/Gamma* > 1/v*^2 for all "
                  f"scanned mu (eps={eps_main}).")
            print("  -> Unexpected; mu < 0 may be the relevant regime.")
        else:
            print("  Mixed result without sign change; INCONCLUSIVE.")

    # Trend toward mu=1/2:
    print()
    print("  Trend: B*/Gamma* - 1/v*^2 as a function of mu (eps=0.1):")
    for mu, d in diffs:
        bar = "+" * max(0, int(round(d * 10))) + "-" * max(0, int(round(-d * 10)))
        print(f"    mu={mu:.2f}  diff={d:+.5f}   {bar}")


if __name__ == "__main__":
    main()
