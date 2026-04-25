"""
OP-EHT-A T-A3 — numerical photon ring with Track A self-consistent coupling.

Cel: Sprawdzić czy naive Track A coupling f(psi) = sqrt(|g_tt|/c_0^2)
=  sqrt((4-3 psi)/psi) prowadzi do photon ring deviation <= 5% (POSITIVE)
lub overshoot/breakdown (NEGATIVE).

Method:
1. Self-consistent equation:
   A_eff = A_baseline × f(psi_ph(A_eff))
   gdzie psi_ph wynika z photon ring solution dla given A_eff.

2. Iteracja:
   A_eff^(0) = A_baseline = 0.5
   psi_ph^(0) = solve photon ring equation for A_eff^(0)
   A_eff^(1) = A_baseline × f(psi_ph^(0))
   etc. do convergencji.

3. Compare to scenarios from OP-EHT T3:
   - Scenario (e): A_eff = A × sqrt(g_tt^GR/g_tt^TGP) → deviation +1.46%
   - Track A:     A_eff = A × sqrt(|g_tt^TGP|/c_0^2) → deviation = ?

PASS criteria:
- T-A3.1: self-consistent solution converges
- T-A3.2: deviation |b_crit/b_GR - 1| <= 5%
- T-A3.3: convergence stable (no oscillation)
"""
import numpy as np
from scipy.optimize import brentq

# c_n PN tail to n=15 (from OP-EHT T1)
C_N = {
    1: 1.0, 2: -1.0, 3: 5.0/3.0, 4: -10.0/3.0, 5: 22.0/3.0,
    6: -154.0/9.0, 7: 374.0/9.0, 8: -935.0/9.0, 9: 21505.0/81.0,
    10: -55913.0/81.0, 11: 147407.0/81.0, 12: -1179256.0/243.0,
    13: 3174920.0/243.0, 14: -8617640.0/243.0, 15: 70664648.0/729.0
}

B_CRIT_GR = 3.0 * np.sqrt(3.0)


def photon_ring(A_val, n_max=15):
    """Compute (r_ph, psi_ph, b_crit) for given matching A.

    Returns None if no photon ring exists or PN tail breaks down.
    """
    def eps(r):
        return sum(C_N[n] * A_val**n / r**n for n in range(1, n_max+1))
    def eps_prime(r):
        return -sum(n * C_N[n] * A_val**n / r**(n+1) for n in range(1, n_max+1))
    def F(r):
        return 4.0 * r * eps_prime(r) + (1.0 - 3.0 * eps(r)) * (1.0 + eps(r))

    # Search for photon ring (scaled with A)
    rs = np.linspace(0.1*A_val/0.5, 30.0*A_val/0.5, 5000)
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
    b_crit = r_ph / A_t
    return r_ph, psi, b_crit


def track_A_weight(psi):
    """Track A coupling factor: f(psi) = sqrt((4-3 psi)/psi).

    At psi=1: factor=1 (1PN preserved)
    At psi=1.168: factor=0.652 (strong-field reduction)
    """
    arg = (4.0 - 3.0 * psi) / psi
    if arg <= 0:
        return None
    return np.sqrt(arg)


def main():
    print("=" * 70)
    print(" OP-EHT-A T-A3 — Self-consistent Track A photon ring")
    print(" Question: does naive proper-time coupling reduce deviation < 5%?")
    print("=" * 70)

    pass_count = 0
    n_total = 3

    # --------------------------------------------------------------
    # T-A3.1: Self-consistent iteration
    # --------------------------------------------------------------
    print(f"\n[T-A3.1] Self-consistent iteration:")
    A_baseline = 0.5
    A_eff = A_baseline  # initial guess
    print(f"  A_baseline = {A_baseline:.4f}")
    print(f"  Initial A_eff = {A_eff:.4f}")
    print(f"  {'iter':>4}  {'A_eff':>8}  {'psi_ph':>8}  {'b_crit':>10}  {'dev_GR':>8}  {'weight':>8}")

    history = []
    converged = False
    for k in range(50):
        result = photon_ring(A_eff)
        if result is None:
            print(f"   {k:>3}  {A_eff:>8.4f}  PN tail breakdown — no photon ring")
            break
        r_ph, psi_ph, b_crit = result
        weight = track_A_weight(psi_ph)
        if weight is None:
            print(f"   {k:>3}  {A_eff:>8.4f}  {psi_ph:>8.4f}  weight invalid (psi too large)")
            break
        dev_GR = (b_crit / B_CRIT_GR - 1.0) * 100
        history.append((k, A_eff, psi_ph, b_crit, dev_GR, weight))
        print(f"   {k:>3}  {A_eff:>8.4f}  {psi_ph:>8.4f}  {b_crit:>10.4f}  {dev_GR:>+7.2f}%  {weight:>8.4f}")

        # Update A_eff = A_baseline × weight (self-consistent)
        A_new = A_baseline * weight
        if abs(A_new - A_eff) < 1e-6:
            A_eff = A_new
            converged = True
            break
        A_eff = A_new

    t_a3_1 = converged
    pass_count += int(t_a3_1)
    print(f"  T-A3.1 self-consistent solution converges: {'PASS' if t_a3_1 else 'FAIL'}")

    # --------------------------------------------------------------
    # T-A3.2: Final deviation
    # --------------------------------------------------------------
    if converged:
        final_A_eff = A_eff
        final_result = photon_ring(final_A_eff)
        if final_result is not None:
            r_ph_f, psi_ph_f, b_crit_f = final_result
            dev_final = (b_crit_f / B_CRIT_GR - 1.0) * 100
            print(f"\n[T-A3.2] Final self-consistent state:")
            print(f"  A_eff_final = {final_A_eff:.4f}")
            print(f"  psi_ph_final = {psi_ph_f:.4f}")
            print(f"  b_crit_final = {b_crit_f:.4f} r_g")
            print(f"  deviation vs GR = {dev_final:+.2f}%")
            t_a3_2 = abs(dev_final) <= 5.0
        else:
            print(f"  Final photon ring computation failed")
            dev_final = None
            t_a3_2 = False
    else:
        dev_final = None
        t_a3_2 = False
    pass_count += int(t_a3_2)
    print(f"  T-A3.2 |deviation| <= 5%: {'PASS' if t_a3_2 else 'FAIL'}")

    # --------------------------------------------------------------
    # T-A3.3: Compare with scenario (e) target
    # --------------------------------------------------------------
    print(f"\n[T-A3.3] Compare with OP-EHT T3 scenario (e) target:")
    print(f"  Scenario (e) target:  A_eff = 0.4428, b_crit = 5.272, dev = +1.46%")
    if dev_final is not None:
        print(f"  Track A self-consistent: A_eff = {final_A_eff:.4f}, b_crit = {b_crit_f:.4f}, dev = {dev_final:+.2f}%")
        gap = abs(dev_final - 1.46)
        print(f"  Gap to scenario (e): {gap:.2f}% (in deviation units)")
        if gap < 2.0:
            print(f"  ⇒ Track A APPROXIMATES scenario (e) within 2%")
            t_a3_3 = True
        else:
            print(f"  ⇒ Track A NIE PASUJE scenario (e) (gap > 2%)")
            t_a3_3 = False
    else:
        t_a3_3 = False
    pass_count += int(t_a3_3)
    print(f"  T-A3.3 Track A approximates scenario (e): {'PASS' if t_a3_3 else 'FAIL'}")

    # --------------------------------------------------------------
    # Summary
    # --------------------------------------------------------------
    print(f"\n{'=' * 70}")
    print(f" SUMMARY: T-A3 = {pass_count}/{n_total} PASS")
    print(f" OP-EHT-A T-A3 verdict: ", end="")
    if pass_count == n_total:
        print("POSITIVE — Track A naive coupling matches scenario (e)")
    elif pass_count >= 1 and t_a3_1:
        print("PARTIAL — converges but magnitude OFF from scenario (e)")
        print(" ⇒ Naive proper-time coupling INSUFFICIENT for M9.1'' rescue")
        print(" ⇒ Refined coupling f(psi) = sqrt(g_tt^GR/g_tt^TGP) needed")
        print("    but lacks first-principles derivation in M9.1''.")
    else:
        print("NEGATIVE — naive proper-time coupling fails (no convergence or breakdown)")
        print(" ⇒ M9.2 pivot becomes inevitable")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
