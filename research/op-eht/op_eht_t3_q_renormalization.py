"""
OP-EHT T3 — Q-renormalization audit at strong field.

Question: czy field-dependent constants G(Φ) = G_0/ψ, c(Φ) = c_0/sqrt(ψ)
prowadzą do strong-field matching correction zmniejszającej +14.56%
b_crit deviation?

Method:
1. T3.1 — Invariance check: G/c² jest invariant pod field-dependent
   constants? Jeśli tak, A = GM/(2c²) niezmienione.
2. T3.2 — Strong-field metric coefficients at photon ring (psi=1.17):
   - g_tt^TGP = -c_0² (4-3·1.17)/1.17 = -c_0² · 0.418
   - Compare to GR Schwarzschild same areal radius
3. T3.3 — Inverse problem: jakie A daje b_crit = b_GR?
4. T3.4 — Q-renormalization plausibility: czy A_required jest realistyczne?
   ADM mass argument vs bare mass.
5. T3.5 — Effective shift of A under various q-renormalization scenarios:
   (a) M_eff = M (no renormalization)
   (b) M_eff = M·sqrt(1+psi)/sqrt(2) (sqrt rescaling)
   (c) M_eff = M·psi (linear in psi)
   (d) M_eff = M_ADM with self-energy

PASS criteria:
- T3.1: G/c² invariance verified analytically
- T3.2: Metric coefficients consistent with EHT-quick + T1
- T3.3: A_required computed
- T3.4: Plausibility verdict
- T3.5: Best q-renorm scenario reduces deviation to <5%? (POSITIVE) or fails?
"""
import numpy as np
import sympy as sp
from scipy.optimize import brentq


def main():
    print("=" * 70)
    print(" OP-EHT T3 — Q-renormalization audit")
    print(" Question: does G(Phi)=G_0/psi reduce +14.56% strong-field deviation?")
    print("=" * 70)

    pass_count = 0
    n_total = 5

    # ================================================================
    # T3.1: G/c^2 invariance
    # ================================================================
    print(f"\n[T3.1] G/c^2 invariance check (analytic):")
    psi = sp.Symbol('psi', positive=True)
    G_TGP = sp.Symbol('G_0') / psi
    c_TGP_sq = sp.Symbol('c_0')**2 / psi
    G_over_c2 = sp.simplify(G_TGP / c_TGP_sq)
    G_over_c2_GR = sp.Symbol('G_0') / sp.Symbol('c_0')**2
    diff = sp.simplify(G_over_c2 - G_over_c2_GR)
    print(f"  G(psi) = G_0/psi")
    print(f"  c(psi)^2 = c_0^2/psi")
    print(f"  G(psi)/c(psi)^2 = {G_over_c2}")
    print(f"  G_0/c_0^2 = {G_over_c2_GR}")
    print(f"  Difference: {diff}")

    t3_1 = (diff == 0)
    pass_count += int(t3_1)
    print(f"  T3.1 G/c^2 invariant: {'PASS' if t3_1 else 'FAIL'}")
    print(f"  ⇒ Standard matching A = GM/(2c^2) = G_0 M / (2 c_0^2) holds in TGP.")
    print(f"  ⇒ Q-renormalization of constants alone CANNOT shift A.")

    # ================================================================
    # T3.2: Strong-field metric coefficients at photon ring
    # ================================================================
    print(f"\n[T3.2] Metric coefficients at photon ring (psi=1.168 from T1):")
    psi_ph = 1.167892  # eps_ph = 0.167892 from T1 / EHT-quick
    g_tt_norm = (4.0 - 3.0 * psi_ph) / psi_ph  # = -g_tt/c_0^2
    print(f"  TGP g_tt/(-c_0^2) at photon ring (psi={psi_ph:.4f}):  {g_tt_norm:.4f}")

    # GR Schwarzschild at areal r = 3M (photon ring): g_tt = -(1 - 2M/r) = -1/3
    g_tt_GR = 1.0/3.0
    print(f"  GR g_tt/(-c^2) at photon ring (r=3M):                 {g_tt_GR:.4f}")
    print(f"  Ratio TGP/GR: {g_tt_norm/g_tt_GR:.4f}")
    print(f"  ⇒ TGP has STRONGER time dilation at photon ring (0.418 vs 0.333)")
    print(f"  ⇒ Source of +14.56% deviation is structural (metric form)")
    t3_2 = abs(g_tt_norm - 0.418) < 0.01
    pass_count += int(t3_2)
    print(f"  T3.2 g_tt^TGP/(-c^2) ~ 0.418: {'PASS' if t3_2 else 'FAIL'}")

    # ================================================================
    # T3.3: Inverse problem - what A gives b_crit = b_GR?
    # ================================================================
    print(f"\n[T3.3] Inverse problem: find A such that b_crit^TGP = b_GR (5.196):")

    # PN tail c_n analytic to n=15 (from T1)
    c_n = {
        1: 1.0, 2: -1.0, 3: 5.0/3.0, 4: -10.0/3.0, 5: 22.0/3.0,
        6: -154.0/9.0, 7: 374.0/9.0, 8: -935.0/9.0, 9: 21505.0/81.0,
        10: -55913.0/81.0, 11: 147407.0/81.0, 12: -1179256.0/243.0,
        13: 3174920.0/243.0, 14: -8617640.0/243.0, 15: 70664648.0/729.0
    }

    def b_crit_at_A(A_val, n_max=15):
        """Compute b_crit for given matching A."""
        def eps(r):
            return sum(c_n[n] * A_val**n / r**n for n in range(1, n_max+1))
        def eps_prime(r):
            return -sum(n * c_n[n] * A_val**n / r**(n+1) for n in range(1, n_max+1))
        def F(r):
            return 4.0 * r * eps_prime(r) + (1.0 - 3.0 * eps(r)) * (1.0 + eps(r))

        # Search for photon ring
        rs = np.linspace(0.1*A_val/0.5, 30.0*A_val/0.5, 5000)  # scaled with A
        Fs = np.array([F(r) for r in rs])
        sign_changes = np.where(np.diff(np.sign(Fs)))[0]
        if len(sign_changes) == 0:
            return None

        idx = sign_changes[0]
        r_ph = brentq(F, rs[idx], rs[idx+1], xtol=1e-10)
        eps_ph = eps(r_ph)
        psi = 1.0 + eps_ph
        denom = 4.0 - 3.0 * psi
        if denom <= 0:
            return None
        A_t = denom / psi
        return r_ph / A_t

    A_baseline = 0.5  # standard 1PN matching
    b_baseline = b_crit_at_A(A_baseline)
    print(f"  A = {A_baseline:.4f} (standard) -> b_crit = {b_baseline:.4f}")

    b_crit_GR = 3.0 * np.sqrt(3.0)
    # Solve b_crit_at_A(A) = b_crit_GR
    # b_crit scales roughly linearly with A in deep PN, so try A_required ~ A_baseline / ratio
    A_grid = np.linspace(0.40, 0.55, 31)
    b_grid = []
    for A in A_grid:
        b = b_crit_at_A(A)
        b_grid.append(b if b is not None else np.nan)
    b_grid = np.array(b_grid)
    print(f"\n  Scan A ∈ [0.40, 0.55]:")
    print(f"    {'A':>6}  {'b_crit':>10}  {'dev vs GR':>10}")
    for A, b in zip(A_grid[::3], b_grid[::3]):
        if not np.isnan(b):
            dev = (b/b_crit_GR - 1.0)*100
            print(f"    {A:>6.3f}  {b:>10.4f}  {dev:>+9.2f}%")

    # Find A that gives b_crit = b_GR via interpolation
    valid = ~np.isnan(b_grid)
    if valid.sum() >= 2:
        from scipy.interpolate import interp1d
        f = interp1d(b_grid[valid], A_grid[valid], kind='linear',
                     fill_value='extrapolate')
        A_required = float(f(b_crit_GR))
        print(f"\n  A_required (b_crit = {b_crit_GR:.4f}): {A_required:.4f}")
        print(f"  A_required / A_baseline: {A_required/A_baseline:.4f}")
        rescaling = A_required / A_baseline
    else:
        A_required = None
        rescaling = None

    t3_3 = A_required is not None
    pass_count += int(t3_3)
    print(f"  T3.3 inverse A computable: {'PASS' if t3_3 else 'FAIL'}")

    # ================================================================
    # T3.4: Plausibility of q-renormalization
    # ================================================================
    print(f"\n[T3.4] Plausibility of q-renormalization rescaling:")
    if A_required is not None:
        delta_M_frac = (rescaling - 1.0)  # M_eff/M_bare - 1
        print(f"  Required M_eff/M_bare = {rescaling:.4f}")
        print(f"  ⇒ effective mass shift {delta_M_frac*100:+.2f}%")
        # ADM mass from BH formation: typically <1% effects, mostly from gravitational binding
        # Self-energy of Schwarzschild BH: ~10% of M (binding energy at horizon)
        # Q-renormalization in TGP: at psi=1.17, the local G/c² shift is 0% (T3.1)
        if abs(delta_M_frac) < 0.02:
            t3_4_class = "TRIVIAL — within ADM bookkeeping"
            t3_4 = True
        elif abs(delta_M_frac) < 0.05:
            t3_4_class = "PLAUSIBLE — within typical self-energy ~5%"
            t3_4 = True
        elif abs(delta_M_frac) < 0.10:
            t3_4_class = "MARGINAL — at edge of self-energy budget"
            t3_4 = False
        else:
            t3_4_class = "IMPLAUSIBLE — exceeds physical self-energy budget"
            t3_4 = False
    else:
        t3_4_class = "INCONCLUSIVE"
        t3_4 = False
    pass_count += int(t3_4)
    print(f"  T3.4 verdict: [{t3_4_class}] {'PASS' if t3_4 else 'FAIL'}")

    # ================================================================
    # T3.5: Effective shifts under q-renorm scenarios
    # ================================================================
    print(f"\n[T3.5] Test q-renormalization scenarios:")
    scenarios = [
        ("(a) No renormalization", 0.5),
        ("(b) sqrt(psi_ph) rescale", 0.5 * np.sqrt(psi_ph)),
        ("(c) 1/psi_ph rescale", 0.5 / psi_ph),
        ("(d) psi_ph linear rescale", 0.5 * psi_ph),
        ("(e) sqrt(g_tt^GR/g_tt^TGP)", 0.5 * np.sqrt(g_tt_GR/g_tt_norm)),
    ]
    print(f"    {'scenario':<45} {'A':>8} {'b_crit':>10} {'dev':>8}")
    print(f"    {'-'*45} {'-'*8} {'-'*10} {'-'*8}")
    best_dev = float('inf')
    best_scenario = None
    for name, A in scenarios:
        b = b_crit_at_A(A)
        if b is not None:
            dev = (b/b_crit_GR - 1.0)*100
            print(f"    {name:<45} {A:>8.4f} {b:>10.4f} {dev:>+7.2f}%")
            if abs(dev) < abs(best_dev):
                best_dev = dev
                best_scenario = name
        else:
            print(f"    {name:<45} {A:>8.4f}      n/a       n/a")

    print(f"\n  Best scenario: {best_scenario} (deviation {best_dev:+.2f}%)")
    if abs(best_dev) < 5.0:
        t3_5_class = "POSITIVE — q-renorm scenario reduces deviation < 5%"
        t3_5 = True
    elif abs(best_dev) < 10.0:
        t3_5_class = "MARGINAL — best deviation in 5-10% range"
        t3_5 = False
    else:
        t3_5_class = "FAIL — no q-renorm scenario reaches < 10%"
        t3_5 = False
    pass_count += int(t3_5)
    print(f"  T3.5 verdict: [{t3_5_class}] {'PASS' if t3_5 else 'FAIL'}")

    # ================================================================
    # Summary
    # ================================================================
    print(f"\n{'=' * 70}")
    print(f" SUMMARY: {pass_count}/{n_total} PASS")
    print(f" OP-EHT T3 verdict: ", end="")
    if pass_count >= 4:
        print("OP-EHT T3 closes POSITIVE — q-renorm absorbs deviation.")
    elif pass_count >= 3:
        print("PARTIAL — invariance + structure verified, but q-renorm cannot save.")
        print("  ⇒ Defers to T4 (environment) or M9.2 pivot.")
    else:
        print("NEGATIVE — q-renormalization path closed.")
        print("  ⇒ T4 (environment) or M9.2 pivot becomes mandatory.")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
