#!/usr/bin/env python3
"""
M6 (Test C / P3.2 of external_review_2026-04-25/review_response_plan.md):
N-operator MK-RG with explicit GL-bond operator perturbation (Track A).

Extension of mk_rg_bgamma.py. Adds a fixed J_GL parameter for the v2
axiom-level GL bond

    H_GL = J_GL * Sum_<ij> A_ij Phi_i Phi_j (Phi_j - Phi_i)^2
         = J_GL * Sum_<ij> [s_i^2 s_j^6 - 2 s_i^4 s_j^4 + s_i^6 s_j^2]

treated as a first-order perturbation of M3's bilinear bond.

Background:
  After M4 (H-S Jacobian, P2/Test A) and M5 (Z_Phi via eta-deformation,
  P3.1/Test B) ruled out two single-channel candidates for closing
  OP-2b, the GL bond is the third (and final) candidate from the
  M3 sec.6 list.

Track A approach (M6_glbond_derivation.md sec. 2-5):
  Project the GL bond onto the on-site V via the s_3 = 0 slice. To
  first order in J_GL, the M3 decimation integral acquires
      DeltaF(s_1) = -J_GL [s_1^2 M_6(K_eff s_1)
                           - 2 s_1^4 M_4(K_eff s_1)
                           + s_1^6 M_2(K_eff s_1)]
  with M_n(h) = <s_2^n>_{exp(-V + h s_2)} (n-th moment under the
  source-shifted single-site weight).
  This DeltaF is then projected onto the even-power basis and added
  via Delta c_{2k} = -2 (2k) a_{2k}, where a_{2k} is the coefficient
  of s^{2k} in DeltaF(s).

Validation: J_GL = 0 must reproduce M3 to 5 decimals.

Decision criteria (M6 sec. 6):
  (1) Closure: B*(J_*)/Gamma*(J_*) = 1/v*^2(J_*) for some J_* in the
      perturbative regime (|J_*| <~ 1). GL bond is the missing physics.
  (2) Closure-in-principle: closure at unphysically large |J_*|.
      Higher-order J_GL^2 terms or NPRG required.
  (3) No closure: no crossing in the convergent regime. All three
      single-channel candidates exhausted.
"""

import numpy as np
from math import comb, factorial

from mk_rg_bgamma import MigdalKadanoffRGN, B_RESCALE, BD1, D


class MigdalKadanoffGL(MigdalKadanoffRGN):
    """
    MK-RG with first-order J_GL perturbation from the GL bond.

    Track A: J_GL fixed (does not flow). The on-site V'-update at
    each MK step receives an additional contribution

        DeltaV'(s) = -2 [DeltaF(s) - DeltaF(0)],
        DeltaF(s) = -J_GL [s^2 M_6(K_eff s) - 2 s^4 M_4(K_eff s)
                                            + s^6 M_2(K_eff s)],

    where M_n(h) is the n-th moment under exp(-V(s_2) + h s_2).
    M_n is computed numerically by reusing the M3 quadrature grid
    for s_2 (inner integral) and a smaller Gauss-Legendre grid for
    s_1 (outer projection grid).

    Note on bond-move: at first order in J_GL the bond-move factor
    b^{d-1} on the GL coupling is implicitly absorbed into the
    user-facing J_GL; closure is reported as a function of J_GL,
    independent of the convention.
    """

    def __init__(self, n_ops, J_GL=0.0, n_quad=1200, s_max=10.0,
                 n_outer=120):
        super().__init__(n_ops=n_ops, n_quad=n_quad, s_max=s_max)
        self.J_GL = float(J_GL)
        self.n_outer = int(n_outer)
        # Outer Gauss-Legendre grid for the s_1 projection.
        nodes_o, weights_o = np.polynomial.legendre.leggauss(self.n_outer)
        self._s_outer = s_max * nodes_o
        self._ds_outer = s_max * weights_o
        # Pre-compute the projection design matrix:
        #   A[i, k] = s_outer[i]^{2(k+1)}  for k = 0..n_ops-1.
        s2_o = self._s_outer ** 2
        A = np.empty((self.n_outer, self.n_ops), dtype=float)
        s2k = np.ones_like(s2_o)
        for k in range(self.n_ops):
            s2k = s2k * s2_o
            A[:, k] = s2k
        self._proj_A = A
        # Use sqrt(|ds_outer|) as L^2 weights for the projection.
        self._proj_W = np.sqrt(np.abs(self._ds_outer))

    def rg_step(self, couplings, K):
        """M3 step plus first-order GL correction in J_GL (Track A)."""
        # Standard M3 update.
        new_c, K_new = super().rg_step(couplings, K)
        if self.J_GL == 0.0:
            return new_c, K_new

        couplings = np.asarray(couplings, dtype=float)
        K_eff = BD1 * K   # 4K; identical to M3.

        # --- M3 single-site log-weight on the inner s_2 grid ---
        s2g = self.s_pts                   # shape (n_quad,)
        s2g_sq = s2g * s2g
        log_w0 = np.zeros_like(s2g)
        s2k_inner = np.ones_like(s2g_sq)
        for k, c in enumerate(couplings):
            s2k_inner = s2k_inner * s2g_sq
            log_w0 = log_w0 - (c / (2.0 * (k + 1))) * s2k_inner

        # --- Outer s_1 grid ---
        s1g = self._s_outer                # shape (n_outer,)
        s1g_sq = s1g * s1g
        s1g_4 = s1g_sq * s1g_sq
        s1g_6 = s1g_4 * s1g_sq

        # --- Source-shifted weight: log_wh[i, j] = log_w0[j] + K_eff s_1[i] s_2[j] ---
        h_outer = K_eff * s1g
        log_wh = log_w0[None, :] + h_outer[:, None] * s2g[None, :]
        # log-sum-exp stabilisation per row.
        max_per_row = np.max(log_wh, axis=1, keepdims=True)
        wh = self.ds[None, :] * np.exp(log_wh - max_per_row)
        Zh = np.sum(wh, axis=1)
        # Guard against accidental zero normalisation (should never trigger).
        Zh = np.where(Zh > 1e-300, Zh, 1e-300)

        # --- Moments M_n(K_eff s_1) for n = 2, 4, 6 ---
        s2_pow2 = s2g_sq
        s2_pow4 = s2g_sq * s2g_sq
        s2_pow6 = s2_pow4 * s2g_sq
        M2 = np.sum(wh * s2_pow2[None, :], axis=1) / Zh
        M4 = np.sum(wh * s2_pow4[None, :], axis=1) / Zh
        M6 = np.sum(wh * s2_pow6[None, :], axis=1) / Zh

        # --- DeltaF(s_1) at first order in J_GL (s_3 = 0 projection) ---
        dF = -self.J_GL * (s1g_sq * M6 - 2.0 * s1g_4 * M4 + s1g_6 * M2)

        # --- Project DeltaF onto even-power basis (Gauss-Legendre L^2) ---
        # dF(s) approx Sum_k a[k] s^{2(k+1)}; weighted lstsq.
        A = self._proj_A
        W = self._proj_W
        AW = A * W[:, None]
        bW = dF * W
        a, *_ = np.linalg.lstsq(AW, bW, rcond=None)

        # --- Apply Delta c_{2(k+1)} = -2 * 2(k+1) * a[k] ---
        for k in range(self.n_ops):
            two_k = 2 * (k + 1)
            new_c[k] = new_c[k] - 2.0 * two_k * a[k]
        return new_c, K_new


def find_fp_with_seeds(rg, seeds, tol=1e-10, max_iter=8000):
    """First seed that converges to a physical WF FP."""
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


def scan_J_GL(n_ops=8, J_list=None, n_quad=1200, s_max=10.0,
              n_outer=120, verbose=True):
    """Continuation-seeded scan of J_GL."""
    if J_list is None:
        J_list = tuple(round(0.025 * i, 4) for i in range(-40, 41))
    base_seeds = [
        (-2.45694, 3.01611, -5.1216, 9.0060),  # M3 8-op FP (head)
        (-2.45694, 3.01611),
        (-2.251, 3.917),
    ]
    results = []
    prev_seed = None
    for J in J_list:
        seeds = []
        if prev_seed is not None:
            seeds.append(prev_seed)
        seeds.extend(base_seeds)
        rg = MigdalKadanoffGL(n_ops=n_ops, J_GL=J,
                              n_quad=n_quad, s_max=s_max,
                              n_outer=n_outer)
        cb_fp, it, conv = find_fp_with_seeds(rg, seeds)
        rb, ub = cb_fp[0], cb_fp[1]
        Bb = cb_fp[2] if n_ops >= 3 else float("nan")
        Gb = cb_fp[3] if n_ops >= 4 else float("nan")
        v_sq = abs(rb) / ub if ub > 0 else float("nan")
        ratio = (Bb / Gb
                 if (n_ops >= 4 and abs(Gb) > 1e-30) else float("nan"))
        target = 1.0 / v_sq if v_sq > 0 else float("nan")
        diff = ((ratio - target)
                if (np.isfinite(ratio) and np.isfinite(target))
                else float("nan"))
        physical = (conv and ub > 0 and rb < 0
                    and np.max(np.abs(cb_fp)) < 1e3)
        row = dict(
            J_GL=J, n_ops=n_ops,
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
            print(f"  J_GL={J:+.4f}  {flag}{phys_flag} it={it:4d}  "
                  f"r*={rb:+.5f}  u*={ub:+.5f}  "
                  f"B*={Bb:+.4e}  G*={Gb:+.4e}  "
                  f"v*^2={v_sq:+.5f}  "
                  f"B*/G*={ratio:+.5f}  1/v*^2={target:+.5f}  "
                  f"diff={diff:+.5f}")
    return results


def main():
    print("M6 -- MK-RG with GL-bond perturbation (Track A; J_GL fixed)")
    print("=" * 78)

    # =====================================================================
    # (0) J_GL = 0 must reproduce M3.
    # =====================================================================
    print()
    print("=== (0) Validation: J_GL=0 must reproduce M3 ===")
    rg0 = MigdalKadanoffGL(n_ops=8, J_GL=0.0,
                           n_quad=1200, s_max=10.0)
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
        print("  -> M3 reproduced (J_GL=0 baseline OK)")
    else:
        print("  -> WARNING: J_GL=0 does NOT reproduce M3.")

    # =====================================================================
    # (0b) Tiny-J_GL sanity: continuity at J_GL = 0.
    # =====================================================================
    print()
    print("=== (0b) Tiny-J_GL continuity check (|J_GL| = 1e-3) ===")
    for J_test in (1e-3, -1e-3):
        rg_t = MigdalKadanoffGL(n_ops=8, J_GL=J_test,
                                n_quad=1200, s_max=10.0)
        cb_t, it_t, conv_t = rg_t.find_fp_bar(cb0, tol=1e-10, max_iter=8000)
        rt, ut = cb_t[0], cb_t[1]
        Bt, Gt = cb_t[2], cb_t[3]
        rt_t = Bt / Gt
        print(f"  J_GL={J_test:+.4f}: r*={rt:+.5f} u*={ut:+.5f} "
              f"B*/G*={rt_t:+.5f} conv={conv_t}")

    # =====================================================================
    # (1) Coarse scan J_GL in [-2, 10], step 0.25.
    #
    # Rationale: the wide-J_GL probe (during P3.2.B development) showed
    # a sign change in (B*/Gamma* - 1/v*^2) at J_GL ~ 5.8, with negative
    # J_GL <~ -2 falling off the convergent branch.
    # =====================================================================
    print()
    print("=== (1) Coarse J_GL-scan, J_GL in [-2.0, 10.0] step 0.25 ===")
    print("  Columns:  r* u* B* Gamma* v*^2 B*/Gamma* 1/v*^2 diff")
    print("  Convergence flags: *=conv, !=non-conv; P=physical (no runaway)")
    print()
    J_grid_coarse = tuple(round(-2.0 + 0.25 * i, 4) for i in range(0, 49))
    results_coarse = scan_J_GL(
        n_ops=8,
        J_list=J_grid_coarse,
        n_quad=1200,
        n_outer=120,
        verbose=True,
    )

    # =====================================================================
    # (2) Fine scan in the *perturbative* regime |J_GL| <= 1, step 0.05.
    #     Track A's first-order-in-J_GL truncation is justified here.
    # =====================================================================
    print()
    print("=== (2) Fine scan in perturbative regime |J_GL| <= 1.0, "
          "step 0.05 ===")
    J_grid_pert = tuple(round(-1.0 + 0.05 * i, 4) for i in range(0, 41))
    results_pert = scan_J_GL(
        n_ops=8,
        J_list=J_grid_pert,
        n_quad=1200,
        n_outer=120,
        verbose=True,
    )

    # =====================================================================
    # (3) Fine scan around the crossing, J_GL in [4.5, 7.0] step 0.05.
    # =====================================================================
    print()
    print("=== (3) Fine scan around crossing, J_GL in [4.5, 7.0] "
          "step 0.05 ===")
    J_grid_cross = tuple(round(4.5 + 0.05 * i, 4) for i in range(0, 51))
    results_cross = scan_J_GL(
        n_ops=8,
        J_list=J_grid_cross,
        n_quad=1200,
        n_outer=120,
        verbose=True,
    )

    # =====================================================================
    # (4) Trend table for physical FPs only (coarse scan).
    # =====================================================================
    print()
    print("=== (4) Trend table for physical FPs only (coarse scan) ===")
    print()
    print("    J_GL   |   r*       u*        B*        Gamma*    "
          "v*^2     B*/Gamma*  1/v*^2   diff")
    print("    -------+--------------------------------------------"
          "------------------------------------")
    for row in results_coarse:
        if not row["physical"]:
            continue
        print(f"    {row['J_GL']:+.3f} | "
              f"{row['r']:+.5f}  {row['u']:+.5f}  "
              f"{row['B']:+.4e}  {row['Gamma']:+.4e}  "
              f"{row['v_sq']:+.5f}  {row['ratio']:+.5f}  "
              f"{row['target']:+.5f}  {row['diff']:+.5f}")

    # =====================================================================
    # (5) Crossing detection on the combined coarse + fine-cross scan.
    # =====================================================================
    print()
    print("=== (5) Sign-change detection in B*/Gamma* - 1/v*^2 ===")
    # Combine coarse and fine-crossing scans, sort by J_GL, dedupe.
    combined_rows = list(results_coarse) + list(results_cross)
    seen = {}
    for r in combined_rows:
        seen[round(r["J_GL"], 4)] = r
    physical = [r for r in sorted(seen.values(), key=lambda r: r["J_GL"])
                if r["physical"] and np.isfinite(r["diff"])]
    crossings = []
    for i in range(1, len(physical)):
        d1 = physical[i - 1]["diff"]
        d2 = physical[i]["diff"]
        if d1 * d2 < 0:
            J1 = physical[i - 1]["J_GL"]
            J2 = physical[i]["J_GL"]
            J_cross = J1 - d1 * (J2 - J1) / (d2 - d1)
            crossings.append((J1, J2, J_cross))
    if crossings:
        for (J1, J2, Jc) in crossings:
            print(f"  CROSSING: B*/Gamma* = 1/v*^2 in J_GL in "
                  f"[{J1:+.3f}, {J2:+.3f}]; "
                  f"linear interp J_cross = {Jc:+.4f}")
    else:
        print("  No sign change found in physical regime.")

    if physical:
        idx_min = int(np.argmin([abs(r["diff"]) for r in physical]))
        rm = physical[idx_min]
        print(f"  Min |diff| = {abs(rm['diff']):.4f}  at J_GL = "
              f"{rm['J_GL']:+.4f}  "
              f"(B*/G*={rm['ratio']:+.5f}, 1/v*^2={rm['target']:+.5f})")

    # =====================================================================
    # (6) Verdict.
    # =====================================================================
    print()
    print("=== (6) Verdict on GL-bond candidate for OP-2b ===")
    # Threshold for "perturbative" Track A regime: |J_GL| <= 1.
    # Track A truncates the J_GL expansion at O(J_GL); the per-step
    # correction is J_GL * <bond_12>.  At the M3 FP <bond_12> ~ O(1-10),
    # so |J_GL| > 1 means the GL correction is no longer a small
    # perturbation of the M3 weight and J_GL^2 terms become important.
    phys_thresh = 1.0
    if crossings:
        any_perturbative = False
        for (_, _, Jc) in crossings:
            if abs(Jc) <= phys_thresh:
                any_perturbative = True
                print(f"  CLOSURE (perturbative) at J_GL = {Jc:+.4f}, "
                      f"|J_GL| <= {phys_thresh}.")
                print("  -> GL-bond operator IS the missing physics; "
                      "OP-2b CLOSED.")
            else:
                print(f"  CLOSURE at J_GL = {Jc:+.4f}, but |J_GL| > "
                      f"{phys_thresh} is outside the perturbative regime.")
                print("  -> GL bond pushes B*/Gamma* in the right "
                      "direction, but at first order in J_GL it is NOT "
                      "the dominant mechanism. CLOSURE-IN-PRINCIPLE only.")
        if not any_perturbative:
            print()
            print("  All crossings are outside |J_GL| <= 1.")
            print("  -> Track A inconclusive; need either (i) full "
                  "two-site Track B with J_GL flowing or (ii) NPRG with "
                  "GL kinetic ansatz to settle whether closure persists "
                  "at the non-perturbative level.")
    else:
        print("  NO CLOSURE in J_GL in [-2.0, 10.0].")
        print("  -> GL bond at first order alone cannot close OP-2b.")
        print("  -> All three single-channel candidates (H-S, Z_Phi, GL) "
              "exhausted at the level of single-operator additions to "
              "single-site MK.")
        print("  -> OP-2b requires NPRG (P3.4) or genuinely new physics.")


if __name__ == "__main__":
    main()
