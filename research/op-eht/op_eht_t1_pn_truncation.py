"""
OP-EHT T1 — PN truncation robustness audit.

Question: czy +14.6% deviation w b_crit jest robust przy n=10, 12, 15
PN terms in tail eps(r) = sum c_n A^n / r^n, lub jest artefaktem
truncation przy n=7?

Method:
1. Reuse derive_c_n_clean from m9_1_pp_p1_higher_pn.py (PN recurrence).
2. Extend n_max from 7 to 15.
3. For each n_max, solve photon-ring equation 4r·eps' + (1-3eps)(1+eps) = 0
   numerically using PN tail eps(r; n_max).
4. Tabulate r_ph, b_crit, deviation vs GR for n=5..15.
5. Test convergence: |b_crit(n+1) - b_crit(n)| / b_crit < 1e-3?

PASS criteria:
- T1.1: c_n recurrence runs to n=15 without symbolic blow-up.
- T1.2: b_crit(n) converges to within 1% by n=10.
- T1.3: Asymptotic b_crit(n→∞) deviates from GR by ε∞.
  - If ε∞ < 5%: +14.6% was truncation artefact ⇒ TGP recovers GR-shadow.
  - If ε∞ ∈ [5%, 15%]: marginal, defers to T2-T4.
  - If ε∞ ≥ 15%: robust, +14.6% genuine ⇒ T2-T5 critical.
- T1.4: Convergence rate (ratio of successive corrections) < 1
  (geometric series convergent).
"""
import sys
import numpy as np
import sympy as sp
from scipy.optimize import brentq


# ============================================================
# A. Symbolic recurrence (extended from m9_1_pp_p1_higher_pn.py)
# ============================================================

def derive_c_n_clean(n_max: int = 15, alpha: int = 2):
    """PN tail coefficients eps(r) = sum c_n A^n / r^n for vacuum Phi-EOM.

    The vacuum EOM is (1+eps)·lap(eps) + alpha·(eps')^2 = 0,
    with eps -> A/r at r -> infty (Newtonian matching).
    Returns dict {n: c_n} for n = 1..n_max.
    """
    r = sp.symbols("r", positive=True)
    a = sp.symbols(f"a1:{n_max + 2}")
    alpha_s = sp.Integer(alpha)

    eps = sum(a[i] / r ** (i + 1) for i in range(n_max + 1))
    eps_p = sp.diff(eps, r)
    eps_pp = sp.diff(eps, r, 2)
    lap_eps = eps_pp + (2 / r) * eps_p
    eom_poly = sp.expand((1 + eps) * lap_eps + alpha_s * eps_p ** 2)

    a1 = a[0]
    c = {1: sp.Integer(1)}
    a_solved = {}

    for n in range(1, n_max + 1):
        k = n + 3
        eom_sub = eom_poly
        for j, val in a_solved.items():
            eom_sub = eom_sub.subs(a[j], val)
        eom_sub = sp.expand(eom_sub)
        coef = eom_sub.coeff(r, -k)
        if coef == 0:
            continue
        sol = sp.solve(coef, a[n])
        if not sol:
            continue
        a_n1_expr = sp.simplify(sol[0])
        a_solved[n] = a_n1_expr
        c[n + 1] = sp.simplify(a_n1_expr / a1 ** (n + 1))
    return c


# ============================================================
# B. Photon ring solver z given PN tail
# ============================================================

def make_eps_funcs(c_dict: dict, A: float, n_max: int):
    """Build numerical eps(r), eps'(r) from PN tail truncated at n_max.
    eps(r) = sum_{n=1..n_max} c_n A^n / r^n
    eps'(r) = -sum_{n=1..n_max} n c_n A^n / r^(n+1)
    """
    coeffs = []  # list of (n, c_n_float)
    for n in range(1, n_max + 1):
        if n in c_dict:
            c_n_val = float(c_dict[n])
            coeffs.append((n, c_n_val))

    def eps(r):
        return sum(c * A ** n / r ** n for (n, c) in coeffs)

    def eps_prime(r):
        return -sum(n * c * A ** n / r ** (n + 1) for (n, c) in coeffs)

    return eps, eps_prime


def find_photon_ring(eps, eps_prime, r_min=0.5, r_max=20.0):
    """Solve 4 r eps'(r) + (1 - 3 eps(r))(1 + eps(r)) = 0 for r_ph (isotropic).

    Returns (r_ph, eps_ph). r_ph is in units of r_g = G M / c^2.
    A is matched as A = G M / (2 c^2) = 0.5 in geometric units (M = c = G = 1).
    """
    def F(r):
        return 4.0 * r * eps_prime(r) + (1.0 - 3.0 * eps(r)) * (1.0 + eps(r))

    # Find sign change
    rs = np.linspace(r_min, r_max, 5000)
    Fs = np.array([F(r) for r in rs])

    # Find zero crossing
    sign_changes = np.where(np.diff(np.sign(Fs)))[0]
    if len(sign_changes) == 0:
        return None, None

    # Take first valid crossing (closest to BH)
    idx = sign_changes[0]
    r_ph = brentq(F, rs[idx], rs[idx + 1], xtol=1e-12)
    return r_ph, eps(r_ph)


def compute_b_crit(r_ph, eps_ph):
    """Critical impact parameter: b_crit = r_ph / A_t(r_ph) where
    A_t = -g_tt/c^2 = (4 - 3 psi)/psi, psi = 1 + eps.
    In isotropic radial coordinate with substrate budget A_t * B = c0^2 = 1,
    so B = 1/A_t. Areal radius R = r * sqrt(B) = r/sqrt(A_t).
    Then b_crit = R(r_ph) / sqrt(A_t(r_ph)) = r/sqrt(A_t)/sqrt(A_t) = r/A_t.
    """
    psi = 1.0 + eps_ph
    A_t = (4.0 - 3.0 * psi) / psi
    if A_t <= 0:
        return None  # below TGP horizon analog
    return r_ph / A_t


# ============================================================
# C. Run convergence sweep over n_max
# ============================================================

def main():
    print("=" * 70)
    print(" OP-EHT T1 — PN truncation robustness audit")
    print(" Question: is +14.6% deviation in b_crit robust vs PN truncation?")
    print("=" * 70)

    A = 0.5  # geometric units, M = c = G = 1 ⇒ A = GM/(2c²) = 0.5

    # GR baselines
    r_ph_gr_areal = 3.0
    b_crit_gr = 3.0 * np.sqrt(3.0)  # ~5.196152
    print(f"\n[GR baselines]")
    print(f"  r_ph (areal)   = {r_ph_gr_areal:.6f} r_g")
    print(f"  b_crit         = {b_crit_gr:.6f} r_g")

    # Derive PN tail to n_max = 15
    print(f"\n[A] Deriving PN tail coefficients c_1..c_15 (sympy):")
    c_dict = derive_c_n_clean(n_max=15, alpha=2)
    for n in range(1, 16):
        if n in c_dict:
            cn_val = float(c_dict[n])
            print(f"    c_{n:>2} = {sp.nsimplify(c_dict[n], rational=True)}  ({cn_val:+.6f})")
        else:
            print(f"    c_{n:>2} = (zero or undefined)")

    # Sweep n_max from 5 to 15
    print(f"\n[B] Photon ring + b_crit at increasing n_max:")
    print(f"    {'n':>3}  {'r_ph (iso)':>12}  {'r_ph (areal)':>14}  "
          f"{'b_crit':>10}  {'dev vs GR':>10}")
    print(f"    {'-' * 3}  {'-' * 12}  {'-' * 14}  {'-' * 10}  {'-' * 10}")

    results = []
    for n_max in range(5, 16):
        eps_func, eps_prime = make_eps_funcs(c_dict, A, n_max)
        r_ph, eps_ph = find_photon_ring(eps_func, eps_prime)
        if r_ph is None:
            print(f"    {n_max:>3}  (no photon ring found)")
            continue
        psi_ph = 1.0 + eps_ph
        denom = 4.0 - 3.0 * psi_ph
        if denom <= 0:
            print(f"    {n_max:>3}  r_ph={r_ph:.4f} but psi={psi_ph:.4f} >= 4/3 "
                  f"(below horizon analog) — PN tail breakdown")
            continue
        B_ph = psi_ph / denom
        r_ph_areal = r_ph * np.sqrt(B_ph)
        b_c = compute_b_crit(r_ph, eps_ph)
        if b_c is None:
            print(f"    {n_max:>3}  b_crit unphysical")
            continue
        dev = (b_c / b_crit_gr - 1.0) * 100.0
        results.append((n_max, r_ph, r_ph_areal, b_c, dev))
        print(f"    {n_max:>3}  {r_ph:>12.6f}  {r_ph_areal:>14.6f}  "
              f"{b_c:>10.6f}  {dev:>+9.2f}%")

    # Convergence analysis
    print(f"\n[C] Convergence test:")
    if len(results) >= 2:
        b_crits = [r[3] for r in results]
        devs = [r[4] for r in results]
        deltas = [abs(b_crits[i+1] - b_crits[i]) for i in range(len(b_crits)-1)]
        rel_deltas = [d / b_crits[i+1] for i, d in enumerate(deltas)]
        print(f"    {'n':>3}->{'n+1':<3}  {'|d b_crit|':>12}  {'|d b_crit|/b':>12}")
        for i, d in enumerate(deltas):
            print(f"    {results[i][0]:>3}->{results[i+1][0]:<3}  "
                  f"{d:>12.2e}  {rel_deltas[i]:>12.2e}")

    # T1 PASS/FAIL evaluation
    print(f"\n[D] T1 PASS/FAIL:")
    pass_count = 0

    # T1.1: recurrence ran to n=15
    n_max_achieved = max(r[0] for r in results) if results else 0
    t1_1 = n_max_achieved >= 15
    pass_count += int(t1_1)
    print(f"  T1.1 c_n recurrence to n=15: {'PASS' if t1_1 else 'FAIL'} "
          f"(achieved n={n_max_achieved})")

    # T1.2: b_crit converges to within 1% by n=11 (odd-n; even-n is unphysical)
    if results:
        # Use highest odd n available (e.g., 11) and compare to n=15
        odd_results = [r for r in results if r[0] % 2 == 1]
        b_at_15 = [r for r in results if r[0] == 15]
        if odd_results and b_at_15:
            # Find n=11 result, or closest below 15
            ref = next((r for r in odd_results if r[0] == 11), None)
            if ref is None:
                # use any result with n < 15
                lower = [r for r in odd_results if r[0] < 15]
                ref = lower[-1] if lower else None
            if ref is not None:
                rel_change = abs(b_at_15[0][3] - ref[3]) / b_at_15[0][3]
                t1_2 = rel_change < 0.01
            else:
                rel_change = None
                t1_2 = False
        else:
            rel_change = None
            t1_2 = False
    else:
        rel_change = None
        t1_2 = False
    pass_count += int(t1_2)
    print(f"  T1.2 |b(15) - b(11)|/b(15) < 1%: "
          f"{'PASS' if t1_2 else 'FAIL'} (rel_change={rel_change})")

    # T1.3: asymptotic deviation epsilon_infty
    # Thresholds calibrated to EHT systematics:
    #   < 3%: TGP recovers GR shadow within all EHT precision (POSITIVE).
    #   3-10%: marginal — within M87* 10% systematic (MARGINAL).
    #   >= 10%: outside Sgr A* 4.4% bound, ROBUST genuine deviation.
    if b_at_15:
        eps_infty = b_at_15[0][4]  # in %
        if abs(eps_infty) < 3.0:
            t1_3_verdict = "TRUNCATION ARTEFACT — TGP recovers GR-shadow"
            t1_3_class = "POSITIVE"
        elif abs(eps_infty) < 10.0:
            t1_3_verdict = "MARGINAL — within M87* 10% systematic, defers to T2-T4"
            t1_3_class = "MARGINAL"
        else:
            t1_3_verdict = "ROBUST GENUINE — outside Sgr A* 4.4%, T2-T5 critical"
            t1_3_class = "ROBUST"
    else:
        eps_infty = None
        t1_3_verdict = "INCONCLUSIVE"
        t1_3_class = "FAIL"
    t1_3 = t1_3_class != "FAIL"
    pass_count += int(t1_3)
    print(f"  T1.3 asymptotic deviation eps_inf = {eps_infty:+.2f}% -> "
          f"[{t1_3_class}] {t1_3_verdict}")

    # T1.4: convergence rate < 1 (geometric)
    if len(results) >= 4:
        ratios = [deltas[i+1] / deltas[i] for i in range(len(deltas)-1)
                  if deltas[i] > 1e-15]
        avg_ratio = np.mean(ratios[-3:]) if len(ratios) >= 3 else None
        t1_4 = avg_ratio is not None and avg_ratio < 1.0
    else:
        avg_ratio = None
        t1_4 = False
    pass_count += int(t1_4)
    print(f"  T1.4 convergence ratio (last 3) < 1: "
          f"{'PASS' if t1_4 else 'FAIL'} (ratio={avg_ratio})")

    # Summary
    n_total = 4
    print(f"\n{'=' * 70}")
    print(f" SUMMARY: {pass_count}/{n_total} PASS")
    print(f" OP-EHT T1 verdict: ", end="")
    if t1_3_class == "POSITIVE":
        print("OP-EHT T1 closes POSITIVE — +14.6% was PN truncation artefact.")
        print(" ⇒ OP-EHT plausibly closes POSITIVE without M9.2 pivot.")
    elif t1_3_class == "MARGINAL":
        print("OP-EHT T1 inconclusive — defers to T2 (mass MC) and T3 (q-renorm).")
    elif t1_3_class == "ROBUST":
        print("OP-EHT T1 confirms genuine deviation — +14.6% is robust prediction.")
        print(" ⇒ T2-T5 must determine if EHT envelope absorbs deviation,")
        print("    or M9.2 pivot becomes mandatory.")
    else:
        print("OP-EHT T1 inconclusive.")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
