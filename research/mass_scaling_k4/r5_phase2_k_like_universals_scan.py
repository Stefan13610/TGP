"""
================================================================================
  K-LIKE UNIVERSALS SCAN — bridge §9.2 follow-up
================================================================================

  Pytanie: R5 zweryfikował K = int(1/2)*g^(2*alpha)*(g')^2*r^2 dr ~ A^2
           (slope 1.99989 alpha=1, 2.00004 alpha=2 — verified to 0.006%).
           Czy inne kinetyczne integrale skalują się jak A^k dla całkowitego k?

  Test: I(p, q, s) = int g^p * |g'|^q * r^s dr
  Sprawdz slope_emp = log(I(g0_mu)/I(g0_e)) / log(A_mu/A_e)
  dla siatki (p, q, s) i obu substratów (alpha=1, alpha=2).

  "Universal" def: slope (a) near-integer (b) sam dla alpha=1 i alpha=2.

  Goal:
   1. Identyfikacja KLASY universals (jeśli istnieje)
   2. Pattern: czy slope = f(p, q, s) dla prostego f?
   3. Jeśli K~A^2 jest izolowany → strukturalny lock-in K^2 = A^4
   4. Jeśli klasa szersza → nowy nurt: topological universals w TGP

  Metoda:
   - g0_e = 0.86941, g0_mu = phi*g0_e = 1.40674 (golden ratio scaling)
   - ODE: g'' = (1-g)*g^(2-2*alpha) - (alpha/g)*g'^2 - (d-1)/r * g'
   - Substraty: alpha=1 (R5 substrate), alpha=2 (TGP-canonical)
   - Skan (p, q, s) ∈ {0, 1, 2, 3, 4} × {0, 1, 2, 3, 4} × {0, 1, 2, 3}

  Dimension counting (linear tail g ~ A*sin(r)/r):
    I ~ A^(p+q) * int [sin/r]^p * [stuff/r]^q * r^s dr
    Convergence at infty: s - p - q < -1 (high-r cutoff matters)
    Naive prediction: slope_linear = p + q (jeśli linear tail dominates)

  Empirical (R5 K~A^2): p+q=4 ale slope=2 — NIE linear tail dominate.
  Sugeruje że NONLINEAR core regime jest istotne.
================================================================================
"""

from __future__ import annotations
import math
import numpy as np
from scipy.integrate import solve_ivp, trapezoid
from scipy.optimize import curve_fit

# ================================================================
# Stałe
# ================================================================

PHI = (1.0 + math.sqrt(5.0)) / 2.0   # golden ratio
G0_E = 0.86941                        # electron g0 (R5 substrate calibration)
G0_MU = PHI * G0_E                    # muon via golden ratio scaling


# ================================================================
# ODE solver (z r5_phase2_analytical_bridge)
# ================================================================

def solve_ode(g0, alpha, d=3, r_max=300.0, n_points=30000, g_floor=1e-10):
    """ODE: g'' = (1-g)*g^(2-2alpha) - (alpha/g)*g'^2 - ((d-1)/r)*g'"""
    singular = [False]
    def rhs(r, y):
        g, gp = y
        if g < g_floor:
            singular[0] = True
            g = g_floor
        rhs_val = (1.0 - g) * g**(2 - 2*alpha)
        if r < 1e-12:
            gpp = rhs_val / max(d, 1.0)
        else:
            gpp = rhs_val - (alpha/g) * gp**2 - ((d-1.0)/r) * gp
        return [gp, gpp]
    r_eval = np.linspace(1e-10, r_max, n_points)
    sol = solve_ivp(rhs, (1e-10, r_max), [g0, 0.0],
                    method='RK45', t_eval=r_eval,
                    rtol=1e-10, atol=1e-12, max_step=0.05)
    return sol, singular[0]


def extract_atail(r, g, r_min=80.0, r_max=250.0):
    mask = (r >= r_min) & (r <= r_max)
    r_f = r[mask]
    if len(r_f) < 10:
        return None
    u_f = (g[mask] - 1.0) * r_f
    def model(r, B, C):
        return B * np.cos(r) + C * np.sin(r)
    try:
        popt, _ = curve_fit(model, r_f, u_f, p0=[0.01, 0.01])
        return math.sqrt(popt[0]**2 + popt[1]**2)
    except:
        return None


def integrate_pqs(r, g, gp, p, q, s, r_cutoff=None):
    """I(p, q, s) = int_0^r_cutoff g^p * |g'|^q * r^s dr"""
    if r_cutoff is None:
        r_cutoff = r[-1]
    mask = r <= r_cutoff
    r_use = r[mask]
    g_use = g[mask]
    gp_use = gp[mask]
    # safety: |g| floor żeby uniknąć log(0) lub g^p z g≈0
    g_safe = np.where(np.abs(g_use) > 1e-10, g_use, 1e-10)
    integrand = (g_safe**p) * (np.abs(gp_use)**q) * (r_use**s)
    return trapezoid(integrand, r_use)


def get_solution(g0, alpha, r_max=300.0):
    sol, sing = solve_ode(g0, alpha, r_max=r_max)
    if sing or not sol.success:
        return None
    r = sol.t
    g = sol.y[0]
    gp = np.gradient(g, r)
    A = extract_atail(r, g)
    return {'r': r, 'g': g, 'gp': gp, 'A': A}


# ================================================================
# MAIN
# ================================================================

def main():
    print("=" * 78)
    print("  K-LIKE UNIVERSALS SCAN — bridge §9.2 follow-up")
    print("=" * 78)
    print()
    print("  Pytanie: czy I(p, q, s) = int g^p |g'|^q r^s dr ~ A^k uniwersalnie?")
    print()

    # Solve dla obu substratów i obu lepton g0
    print("  Solving ODE...")
    print(f"    g0_e = {G0_E:.5f}, g0_mu = phi*g0_e = {G0_MU:.5f}")
    print()

    sols = {}
    for alpha in (1, 2):
        sols[(alpha, 'e')] = get_solution(G0_E, alpha)
        sols[(alpha, 'mu')] = get_solution(G0_MU, alpha)
        if sols[(alpha, 'e')] is None or sols[(alpha, 'mu')] is None:
            print(f"  ERROR: ODE singular dla alpha={alpha}")
            return

    for alpha in (1, 2):
        A_e = sols[(alpha, 'e')]['A']
        A_mu = sols[(alpha, 'mu')]['A']
        print(f"  alpha={alpha}: A_e = {A_e:.6f}, A_mu = {A_mu:.6f}")
    print()

    # Skan (p, q, s)
    P_VALS = [0, 1, 2, 3, 4]
    Q_VALS = [0, 1, 2, 3, 4]
    S_VALS = [0, 1, 2, 3]

    print("=" * 78)
    print("  CASE 1: full scan (p, q, s) — alpha=1 substrate")
    print("=" * 78)
    print()
    print(f"  {'p':>3} {'q':>3} {'s':>3} | {'slope_emp':>10} | {'I(e)':>12} | {'I(mu)':>14}")
    print("  " + "-" * 68)

    results_a1 = []
    A_e1 = sols[(1, 'e')]['A']
    A_mu1 = sols[(1, 'mu')]['A']
    log_A_ratio_1 = math.log(A_mu1 / A_e1)

    for p in P_VALS:
        for q in Q_VALS:
            for s in S_VALS:
                try:
                    I_e = integrate_pqs(sols[(1, 'e')]['r'],
                                        sols[(1, 'e')]['g'],
                                        sols[(1, 'e')]['gp'], p, q, s)
                    I_mu = integrate_pqs(sols[(1, 'mu')]['r'],
                                         sols[(1, 'mu')]['g'],
                                         sols[(1, 'mu')]['gp'], p, q, s)
                    if I_e <= 0 or I_mu <= 0:
                        slope = float('nan')
                    else:
                        slope = math.log(I_mu / I_e) / log_A_ratio_1
                    results_a1.append((p, q, s, slope, I_e, I_mu))
                except Exception:
                    results_a1.append((p, q, s, float('nan'), float('nan'), float('nan')))

    # Pokaz tylko slope's bliskie integerowi (within 1%)
    print()
    print("  ONLY near-integer slopes (|slope - round(slope)| < 0.01):")
    print()
    near_integer_a1 = []
    for p, q, s, slope, I_e, I_mu in results_a1:
        if math.isnan(slope):
            continue
        if abs(slope - round(slope)) < 0.01:
            near_integer_a1.append((p, q, s, slope, I_e, I_mu))
            print(f"  {p:>3} {q:>3} {s:>3} | {slope:>10.5f} | {I_e:>12.4e} | {I_mu:>14.4e}")
    print()

    print("=" * 78)
    print("  CASE 2: full scan (p, q, s) — alpha=2 canonical")
    print("=" * 78)
    print()

    A_e2 = sols[(2, 'e')]['A']
    A_mu2 = sols[(2, 'mu')]['A']
    log_A_ratio_2 = math.log(A_mu2 / A_e2)

    results_a2 = []
    for p in P_VALS:
        for q in Q_VALS:
            for s in S_VALS:
                try:
                    I_e = integrate_pqs(sols[(2, 'e')]['r'],
                                        sols[(2, 'e')]['g'],
                                        sols[(2, 'e')]['gp'], p, q, s)
                    I_mu = integrate_pqs(sols[(2, 'mu')]['r'],
                                         sols[(2, 'mu')]['g'],
                                         sols[(2, 'mu')]['gp'], p, q, s)
                    if I_e <= 0 or I_mu <= 0:
                        slope = float('nan')
                    else:
                        slope = math.log(I_mu / I_e) / log_A_ratio_2
                    results_a2.append((p, q, s, slope, I_e, I_mu))
                except Exception:
                    results_a2.append((p, q, s, float('nan'), float('nan'), float('nan')))

    print("  ONLY near-integer slopes (|slope - round(slope)| < 0.01):")
    print()
    near_integer_a2 = []
    for p, q, s, slope, I_e, I_mu in results_a2:
        if math.isnan(slope):
            continue
        if abs(slope - round(slope)) < 0.01:
            near_integer_a2.append((p, q, s, slope, I_e, I_mu))
            print(f"  {p:>3} {q:>3} {s:>3} | {slope:>10.5f} | {I_e:>12.4e} | {I_mu:>14.4e}")
    print()

    # Cross-check: które (p, q, s) są UNIVERSAL (ten sam slope dla obu α)?
    print("=" * 78)
    print("  CASE 3: UNIVERSAL across alpha (slope_alpha=1 ~= slope_alpha=2)")
    print("=" * 78)
    print()
    print(f"  {'p':>3} {'q':>3} {'s':>3} | {'slope_a1':>10} | {'slope_a2':>10} | {'diff':>8}")
    print("  " + "-" * 56)

    universals = []
    for i, (p, q, s, slope_a1, _, _) in enumerate(results_a1):
        slope_a2 = results_a2[i][3]
        if math.isnan(slope_a1) or math.isnan(slope_a2):
            continue
        diff = abs(slope_a1 - slope_a2)
        if diff < 0.05:  # within 5% — universal
            universals.append((p, q, s, slope_a1, slope_a2))
            mark = " ok" if diff < 0.005 else ""
            print(f"  {p:>3} {q:>3} {s:>3} | {slope_a1:>10.5f} | {slope_a2:>10.5f} | {diff:>8.5f}{mark}")
    print()
    print(f"  Total universals (within 5%): {len(universals)}")
    print(f"  Total tight (within 0.5%): {sum(1 for u in universals if abs(u[3]-u[4]) < 0.005)}")
    print()

    # Pattern hunt: czy slope = a*p + b*q + c*s + d?
    print("=" * 78)
    print("  CASE 4: Pattern hunt — czy slope = a*p + b*q + c*s + d?")
    print("=" * 78)
    print()

    # Linear regression na universals: slope ~ a*p + b*q + c*s + d
    if len(universals) >= 4:
        # Use only tight universals (diff < 0.5%)
        tight = [u for u in universals if abs(u[3]-u[4]) < 0.005 and abs(u[3] - round(u[3])) < 0.05]
        if len(tight) >= 4:
            X = np.array([[u[0], u[1], u[2], 1.0] for u in tight])
            y = np.array([(u[3] + u[4])/2 for u in tight])  # average slopes
            try:
                coef, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
                a, b, c, d = coef
                print(f"  Linear fit (tight universals, n={len(tight)}):")
                print(f"    slope ~= {a:.4f}*p + {b:.4f}*q + {c:.4f}*s + {d:.4f}")
                print()
                # Sprawdz residua
                print(f"  {'p':>3} {'q':>3} {'s':>3} | {'slope_emp':>10} | {'slope_fit':>10} | {'residue':>10}")
                print("  " + "-" * 60)
                for u in tight:
                    p, q, s = u[0], u[1], u[2]
                    slope_emp = (u[3] + u[4]) / 2
                    slope_fit = a*p + b*q + c*s + d
                    res = slope_emp - slope_fit
                    print(f"  {p:>3} {q:>3} {s:>3} | {slope_emp:>10.5f} | {slope_fit:>10.5f} | {res:>10.5f}")
                print()
            except np.linalg.LinAlgError:
                print("  Linear fit failed (singular)")

    # Specjalne testy: K-like (p=2*alpha, q=2, s=2) dla obu alpha
    print("=" * 78)
    print("  CASE 5: R5 K-action verification (p, q, s) = (2*alpha, 2, 2)")
    print("=" * 78)
    print()

    for alpha in (1, 2):
        p_K = 2 * alpha
        q_K = 2
        s_K = 2
        idx = next(i for i, r in enumerate(results_a1)
                   if r[0] == p_K and r[1] == q_K and r[2] == s_K)
        if alpha == 1:
            slope = results_a1[idx][3]
        else:
            slope = results_a2[idx][3]
        print(f"  alpha={alpha}: K = int(1/2)*g^{p_K}*(g')^{q_K}*r^{s_K} dr")
        print(f"           slope_emp = {slope:.5f}  (R5 verified ~ A^2)")
    print()

    # Wnioski
    print("=" * 78)
    print("  KONKLUZJA")
    print("=" * 78)
    print()
    print(f"  - Total skanowanych (p, q, s): {len(results_a1)}")
    print(f"  - Near-integer slope alpha=1: {len(near_integer_a1)}")
    print(f"  - Near-integer slope alpha=2: {len(near_integer_a2)}")
    print(f"  - Universal (slope_a1 = slope_a2 within 5%): {len(universals)}")
    print(f"  - Tight universal (within 0.5%): "
          f"{sum(1 for u in universals if abs(u[3]-u[4]) < 0.005)}")
    print()


if __name__ == "__main__":
    main()
