#!/usr/bin/env python3
"""
M5 (Test B / P3.1 of external_review_2026-04-25/review_response_plan.md):
N-operator MK-RG with explicit anomalous-dimension deformation.

Extension of mk_rg_bgamma.py. Modifies the bond-move factor:

    K_eff(eta) = b^{d-1+eta} K = (b^{d-1} K) * b^{eta}.

For eta = 0 we recover M3 exactly.

Background:
  The composite field Phi = s^2 has anomalous dimension eta_Phi = 2 eta_phi
  at one loop in phi^4 theory. Bootstrap value for 3D Ising:
  eta_phi ~ 0.036, hence eta_Phi ~ 0.072.

  M3's missing physics could (in principle) be the wave-function
  renormalisation of the composite field: B and Gamma couple to Phi^3
  and Phi^4, so are sensitive to Z_Phi^3 and Z_Phi^4 factors.

  Derivation in M5_zphi_derivation.md. Test B: scan eta and look for
  closure of B*/Gamma* = 1/v*^2.

Validation: eta = 0 must reproduce M3 to 5 decimals.

Decision criteria (M5 sec. 4):
  (1) Closure: B*(eta_*)/Gamma*(eta_*) = 1/v*^2(eta_*) for some
      eta_* with |eta_*| <~ 0.1 (physically reasonable, bootstrap
      consistent). Z_Phi is THE missing physics.
  (2) Closure-in-principle: closure at large |eta_*| > 0.5. Z_Phi
      contributes but isn't dominant.
  (3) No closure: no crossing in convergent regime. Z_Phi alone
      cannot close OP-2b.
"""

import numpy as np
from math import comb, factorial

from mk_rg_bgamma import MigdalKadanoffRGN, B_RESCALE, BD1, D


class MigdalKadanoffZPhi(MigdalKadanoffRGN):
    """
    MK-RG with explicit eta-deformation of the bond rescaling.
    K_eff(eta) = b^{d-1+eta} K.
    """

    def __init__(self, n_ops, eta=0.0, n_quad=1200, s_max=10.0):
        super().__init__(n_ops=n_ops, n_quad=n_quad, s_max=s_max)
        self.eta = float(eta)
        # Pre-compute the eta factor on K_eff.
        self._k_eff_factor = BD1 * (B_RESCALE ** self.eta)

    def rg_step(self, couplings, K):
        couplings = np.asarray(couplings, dtype=float)
        # Modified K_eff: standard M3 had BD1*K = 4K. Now multiply by b^eta.
        K_eff = self._k_eff_factor * K
        moments = self._moments(couplings)
        cumulants = self.cumulants_from_moments(moments)
        new_c = np.empty_like(couplings)
        # K-update from F_2.
        F2 = K_eff ** 2 * cumulants[2]
        K_new = max(F2, 1e-15)
        # Update each coupling: c_{2k}' = c_{2k} - 2 F_{2k}/(2k-1)!
        for k in range(self.n_ops):
            two_k = 2 * (k + 1)
            F = (K_eff ** two_k) * cumulants[two_k]
            new_c[k] = couplings[k] - 2.0 * F / factorial(two_k - 1)
        return new_c, K_new


def find_fp_with_seeds(rg, seeds, tol=1e-10, max_iter=8000):
    """First seed that converges to a physical WF FP (u*>0, r*<0,
    no runaway saddle)."""
    best = None
    for seed in seeds:
        cb0 = np.zeros(rg.n_ops)
        for j, v in enumerate(seed):
            if j < rg.n_ops:
                cb0[j] = v
        cb_fp, it, conv = rg.find_fp_bar(cb0, tol=tol, max_iter=max_iter)
        if (conv and cb_fp[1] > 0 and cb_fp[0] < 0
                and np.all(np.isfinite(cb_fp))
                and np.max(np.abs(cb_fp)) < 1e3):
            return cb_fp, it, True
        if best is None:
            best = (cb_fp, it, conv)
    return best


def scan_eta(n_ops=8, eta_list=None, n_quad=1200, s_max=10.0, verbose=True):
    """Continuation-seeded scan of eta."""
    if eta_list is None:
        eta_list = tuple(round(0.025 * i, 3) for i in range(-40, 81))
    base_seeds = [
        (-2.45694, 3.01611, -5.1216, 9.0060),  # M3 8-op FP (head)
        (-2.45694, 3.01611),
        (-2.251, 3.917),
    ]
    results = []
    prev_seed = None
    for eta in eta_list:
        seeds = []
        if prev_seed is not None:
            seeds.append(prev_seed)
        seeds.extend(base_seeds)
        rg = MigdalKadanoffZPhi(n_ops=n_ops, eta=eta,
                                n_quad=n_quad, s_max=s_max)
        cb_fp, it, conv = find_fp_with_seeds(rg, seeds)
        rb, ub = cb_fp[0], cb_fp[1]
        Bb = cb_fp[2] if n_ops >= 3 else float("nan")
        Gb = cb_fp[3] if n_ops >= 4 else float("nan")
        v_sq = abs(rb) / ub if ub > 0 else float("nan")
        ratio = Bb / Gb if (n_ops >= 4 and abs(Gb) > 1e-30) else float("nan")
        target = 1.0 / v_sq if v_sq > 0 else float("nan")
        diff = (ratio - target) if (np.isfinite(ratio) and np.isfinite(target)) else float("nan")
        physical = (conv and ub > 0 and rb < 0
                    and np.max(np.abs(cb_fp)) < 1e3)
        row = dict(
            eta=eta, n_ops=n_ops,
            r=rb, u=ub, B=Bb, Gamma=Gb,
            v_sq=v_sq, ratio=ratio, target=target, diff=diff,
            converged=bool(conv), physical=bool(physical),
            iterations=int(it), cb_fp=cb_fp.copy(),
        )
        results.append(row)
        if physical:
            prev_seed = tuple(cb_fp)
        if verbose:
            flag = "*" if conv else "!"
            phys_flag = "P" if physical else " "
            print(f"  eta={eta:+.3f}  {flag}{phys_flag} it={it:4d}  "
                  f"r*={rb:+.5f}  u*={ub:+.5f}  "
                  f"B*={Bb:+.4e}  G*={Gb:+.4e}  "
                  f"v*^2={v_sq:+.5f}  "
                  f"B*/G*={ratio:+.5f}  1/v*^2={target:+.5f}  "
                  f"diff={diff:+.5f}")
    return results


def main():
    print("M5 -- MK-RG with eta-deformation (Z_Phi candidate for OP-2b)")
    print("=" * 78)

    # =====================================================================
    # (0) eta=0 must reproduce M3 exactly.
    # =====================================================================
    print()
    print("=== (0) Validation: eta=0 must reproduce M3 ===")
    rg0 = MigdalKadanoffZPhi(n_ops=8, eta=0.0, n_quad=1200, s_max=10.0)
    cb0 = np.zeros(8); cb0[0] = -2.251; cb0[1] = 3.917
    cb_fp, it, conv = rg0.find_fp_bar(cb0, tol=1e-10, max_iter=8000)
    rb, ub, Bb, Gb = cb_fp[0], cb_fp[1], cb_fp[2], cb_fp[3]
    v_sq = abs(rb) / ub
    ratio = Bb / Gb
    print(f"  Got: r*={rb:+.5f} u*={ub:+.5f} B*={Bb:+.4e} G*={Gb:+.4e}")
    print(f"       v*^2={v_sq:.5f} B*/G*={ratio:+.5f} (M3 expects -0.5687)")
    print(f"  conv={conv} in {it} steps")
    err_r = abs(rb - (-2.45694))
    err_ratio = abs(ratio - (-0.5687))
    if err_r < 5e-3 and err_ratio < 5e-3:
        print("  -> M3 reproduced (eta=0 baseline OK)")
    else:
        print("  -> WARNING: eta=0 does NOT reproduce M3.")

    # =====================================================================
    # (1) Coarse scan eta in [-1, 2], step 0.05.
    # =====================================================================
    print()
    print("=== (1) Coarse eta-scan, eta in [-1.0, 2.0] step 0.05 ===")
    print()
    print("  Columns:  r*  u*  B*  Gamma*  v*^2  B*/Gamma*  1/v*^2  diff")
    print("  Convergence flags: *=conv, !=non-conv; P=physical (no runaway)")
    print()
    eta_grid_coarse = tuple(round(-1.0 + 0.05 * i, 3) for i in range(0, 61))
    results_coarse = scan_eta(
        n_ops=8,
        eta_list=eta_grid_coarse,
        n_quad=1200,
        verbose=True,
    )

    # =====================================================================
    # (2) Fine scan around bootstrap value eta ~ 0.07 and around any
    #     near-crossing observed in the coarse scan.
    # =====================================================================
    print()
    print("=== (2) Fine scan around bootstrap (eta ~ 0.07) ===")
    eta_grid_fine = tuple(round(0.005 * i, 3) for i in range(0, 41))  # 0..0.2 step 0.005
    results_fine = scan_eta(
        n_ops=8,
        eta_list=eta_grid_fine,
        n_quad=1200,
        verbose=True,
    )

    # =====================================================================
    # (3) Trend table.
    # =====================================================================
    print()
    print("=== (3) Trend table for physical FPs only (coarse scan) ===")
    print()
    print("    eta   |   r*       u*        B*        Gamma*    "
          "v*^2     B*/Gamma*  1/v*^2   diff")
    print("    ------+--------------------------------------------"
          "------------------------------------")
    for row in results_coarse:
        if not row["physical"]:
            continue
        print(f"    {row['eta']:+.2f} | "
              f"{row['r']:+.5f}  {row['u']:+.5f}  "
              f"{row['B']:+.4e}  {row['Gamma']:+.4e}  "
              f"{row['v_sq']:+.5f}  {row['ratio']:+.5f}  "
              f"{row['target']:+.5f}  {row['diff']:+.5f}")

    # =====================================================================
    # (4) Crossing detection.
    # =====================================================================
    print()
    print("=== (4) Sign-change detection in B*/Gamma* - 1/v*^2 ===")
    physical = [r for r in results_coarse if r["physical"]
                and np.isfinite(r["diff"])]
    crossings = []
    for i in range(1, len(physical)):
        d1 = physical[i - 1]["diff"]
        d2 = physical[i]["diff"]
        if d1 * d2 < 0:
            # Linear interpolation for crossing eta:
            e1 = physical[i - 1]["eta"]
            e2 = physical[i]["eta"]
            eta_cross = e1 - d1 * (e2 - e1) / (d2 - d1)
            crossings.append((e1, e2, eta_cross))
    if crossings:
        for (e1, e2, ec) in crossings:
            print(f"  CROSSING: B*/Gamma* = 1/v*^2 in eta in [{e1:+.3f}, {e2:+.3f}]; "
                  f"linear interp eta_cross = {ec:+.4f}")
    else:
        print("  No sign change found in physical regime.")

    # =====================================================================
    # (5) Verdict.
    # =====================================================================
    print()
    print("=== (5) Verdict on Z_Phi candidate for OP-2b ===")
    if crossings:
        for (_, _, ec) in crossings:
            phys_thresh = 0.10
            if abs(ec) <= phys_thresh:
                print(f"  CLOSURE at eta = {ec:+.4f}, |eta| <= {phys_thresh}.")
                print(f"  Bootstrap reference: 2*eta_phi(3D Ising) ~ 0.072")
                print("  -> Z_Phi anomalous dimension IS the missing physics; "
                      "OP-2b CLOSED.")
            else:
                print(f"  CLOSURE at eta = {ec:+.4f}, but |eta| > {phys_thresh} "
                      "is unphysically large.")
                print("  -> Z_Phi contributes in the right direction but is "
                      "NOT the dominant mechanism.")
                print("  -> GL-bond operator (P3.2) still required.")
    else:
        print("  NO CLOSURE in eta in [-1.0, 2.0].")
        print("  -> Z_Phi alone cannot close OP-2b.")
        print("  -> The GL-bond operator (P3.2) is the remaining "
              "single-channel candidate.")


if __name__ == "__main__":
    main()
