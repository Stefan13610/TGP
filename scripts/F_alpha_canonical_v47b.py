#!/usr/bin/env python3
"""
F_alpha_canonical_v47b.py -- F(alpha) profile in canonical Form A.

In canonical Form A, the ODE is:
  g'' + (d-1)/r * g' + (alpha/g)*g'^2 + g^2(1-g) = 0

For each alpha:
  1. Calibrate g0_e from r_21 = m_mu/m_e = 206.768
     using phi-FP spacing: g0_mu = phi * g0_e
  2. Compute F(alpha) = (A_tail(phi*g0_e) / A_tail(g0_e))^4
  3. F(alpha*) = r_21 at the physical alpha

KEY QUESTION: Does F(alpha) have a minimum? At what alpha?
In Form B (ex92): F_min ~ 201 at alpha_min ~ 2.54.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

PHI = (1 + np.sqrt(5)) / 2
r_21_pdg = 206.7682830  # m_mu / m_e


def make_solver(alpha, d=3, r_max=300):
    def solver(g0):
        def rhs(r, y):
            g, gp = y
            g = max(g, 1e-10)
            source = g**2 * (1.0 - g)
            cross = (alpha / g) * gp**2
            if r < 1e-10:
                return [gp, (source - cross) / float(d)]
            return [gp, source - cross - float(d - 1) * gp / r]
        sol = solve_ivp(rhs, (1e-6, r_max), [g0, 0.0],
                        rtol=1e-12, atol=1e-14, max_step=0.02)
        return sol.t, sol.y[0]
    return solver


def A_tail(solver, g0):
    r, g = solver(g0)
    mask = (r > 50) & (r < 200)
    if np.sum(mask) < 100:
        return 0.0
    rf = r[mask]
    if np.any(np.abs(g[mask] - 1) > 0.5):
        return 0.0
    df = (g[mask] - 1.0) * rf
    M = np.column_stack([np.cos(rf), np.sin(rf)])
    bc = np.linalg.lstsq(M, df, rcond=None)[0]
    return np.sqrt(bc[0]**2 + bc[1]**2)


def compute_F_alpha(alpha):
    """Compute F(alpha) = (A_mu/A_e)^4 for given alpha.

    Method: scan g0_e to find where (A(phi*g0)/A(g0))^4 = any value.
    F(alpha) is the function (A_mu/A_e)^4 evaluated at the calibrated g0_e.

    But we DON'T calibrate to r_21 -- we just compute the ratio for
    a range of g0_e values and report the maximum/minimum/shape.

    Actually, the right approach: F(alpha) depends on which g0_e we choose.
    In the Form B approach, g0_e is determined by B_coeff = 0 (phase condition).
    In canonical Form A, there's no such condition -- g0_e is free.

    The correct analog: calibrate g0_e so that r_21(g0_e) = 206.77.
    If no such g0_e exists, F(alpha) is undefined.
    """
    solver = make_solver(alpha)
    gc = (2*alpha + 4) / (2*alpha + 1)

    # R21 as function of g0_e
    def r21_of_g0e(g0_e):
        if g0_e <= 0.1 or g0_e >= gc/PHI:
            return 0.0
        A_e = A_tail(solver, g0_e)
        g0_mu = PHI * g0_e
        if g0_mu >= gc:
            return 0.0
        A_mu = A_tail(solver, g0_mu)
        if A_e < 1e-15:
            return 0.0
        return (A_mu / A_e) ** 4

    # Scan g0_e range
    g0_max = min(0.98, gc/PHI - 0.01)
    if g0_max <= 0.3:
        return None

    g0_scan = np.linspace(0.3, g0_max, 40)
    r21_vals = []
    for g0 in g0_scan:
        r = r21_of_g0e(g0)
        r21_vals.append(r)

    r21_vals = np.array(r21_vals)

    # Find where r21 = r_21_pdg
    bracket = None
    for i in range(len(r21_vals) - 1):
        if r21_vals[i] > 0 and r21_vals[i+1] > 0:
            if (r21_vals[i] - r_21_pdg) * (r21_vals[i+1] - r_21_pdg) < 0:
                bracket = (g0_scan[i], g0_scan[i+1])
                break

    if bracket is None:
        # Report range of r21
        valid = r21_vals[r21_vals > 0]
        if len(valid) > 0:
            return {'alpha': alpha, 'gc': gc, 'status': 'NO_CALIB',
                    'r21_min': np.min(valid), 'r21_max': np.max(valid)}
        return None

    # Calibrate g0_e
    def r21_res(g0_e):
        return r21_of_g0e(g0_e) - r_21_pdg

    try:
        g0_e = brentq(r21_res, bracket[0], bracket[1], xtol=1e-8)
    except:
        return None

    g0_mu = PHI * g0_e
    A_e = A_tail(solver, g0_e)
    A_mu = A_tail(solver, g0_mu)
    F_val = (A_mu / A_e) ** 4

    # Also check Koide: find g0_tau
    g0_tau = None
    Q_K = None
    g0_lo = g0_mu + 0.005
    g0_hi = gc - 0.001
    if g0_hi > g0_lo:
        def koide_res(g0_3):
            A3 = A_tail(solver, g0_3)
            if A3 < 1e-15:
                return 1e10
            S2 = A_e**2 + A_mu**2 + A3**2
            S4 = A_e**4 + A_mu**4 + A3**4
            return S2**2 / S4 - 1.5

        g0_tau_scan = np.linspace(g0_lo, g0_hi, 80)
        k_resids = [koide_res(g0) for g0 in g0_tau_scan]

        for i in range(len(k_resids) - 1):
            if k_resids[i] != 1e10 and k_resids[i+1] != 1e10:
                if k_resids[i] * k_resids[i+1] < 0:
                    try:
                        g0_tau = brentq(koide_res, g0_tau_scan[i], g0_tau_scan[i+1],
                                        xtol=1e-10)
                        A_tau = A_tail(solver, g0_tau)
                        S2 = A_e**2 + A_mu**2 + A_tau**2
                        S4 = A_e**4 + A_mu**4 + A_tau**4
                        Q_K = S2**2 / S4
                    except:
                        pass
                    break

    r31 = None
    if g0_tau is not None:
        A_tau = A_tail(solver, g0_tau)
        r31 = (A_tau / A_e) ** 4

    return {
        'alpha': alpha, 'gc': gc,
        'g0_e': g0_e, 'g0_mu': g0_mu,
        'A_e': A_e, 'A_mu': A_mu,
        'F': F_val, 'status': 'CALIBRATED',
        'g0_tau': g0_tau, 'Q_K': Q_K, 'r31': r31,
        'margin': gc - (g0_tau if g0_tau else g0_mu),
    }


# ================================================================
print("=" * 70)
print("F(alpha) PROFILE -- CANONICAL FORM A")
print("=" * 70)

# Scan alpha values
alpha_vals = np.arange(0.5, 5.1, 0.25)

print(f"\n  {'alpha':>6s} {'gc':>8s} {'g0_e':>8s} {'g0_mu':>8s} {'F':>10s} "
      f"{'g0_tau':>8s} {'Q_K':>10s} {'r31':>10s} {'margin':>8s}")

results = []
for alpha in alpha_vals:
    res = compute_F_alpha(alpha)
    results.append(res)

    if res is None:
        print(f"  {alpha:6.2f} {'---':>8s}")
    elif res['status'] == 'NO_CALIB':
        print(f"  {alpha:6.2f} {res['gc']:8.4f} {'---':>8s} {'---':>8s} "
              f"{'NO_CALIB':>10s} r21:[{res['r21_min']:.1f},{res['r21_max']:.1f}]")
    elif res['status'] == 'CALIBRATED':
        g0t_str = f"{res['g0_tau']:.4f}" if res['g0_tau'] else "---"
        qk_str = f"{res['Q_K']:.6f}" if res['Q_K'] else "---"
        r31_str = f"{res['r31']:.1f}" if res['r31'] else "---"
        print(f"  {alpha:6.2f} {res['gc']:8.4f} {res['g0_e']:8.4f} "
              f"{res['g0_mu']:8.4f} {res['F']:10.2f} "
              f"{g0t_str:>8s} {qk_str:>10s} {r31_str:>10s} {res['margin']:8.4f}")


# ================================================================
# Analysis
# ================================================================
print("\n\n" + "=" * 70)
print("ANALYSIS")
print("=" * 70)

# F(alpha) should be constant = r_21 if calibration is correct
# So the profile is trivial in the canonical approach!
# The reason: in canonical Form A, we CALIBRATE g0_e to match r_21.
# There's no "phase condition" like B_coeff = 0.

print("""
  KEY INSIGHT: In canonical Form A, F(alpha) = r_21 BY CONSTRUCTION.
  We calibrate g0_e for each alpha so that (A_mu/A_e)^4 = r_21.
  F(alpha) has no profile -- it's a constant.

  The nontrivial information is:
  1. Does g0_e(alpha) exist? (calibration succeeds)
  2. What is g0_tau(alpha) from Koide? (Q_K = 3/2)
  3. What is r31(alpha) = m_tau/m_e?
  4. Does r31 match PDG (3477.23)?

  The Form B "F(alpha) profile" came from a DIFFERENT selection rule
  (B_coeff = 0 fixes g0_e independently of r_21).
  In Form A, there's no such constraint, so g0_e is free.
""")

# Show r31 vs alpha
print("  r31(alpha) from Koide (Form A):")
print(f"  {'alpha':>6s} {'r31':>10s} {'r31/PDG':>10s}")
for res in results:
    if res and res.get('r31'):
        print(f"  {res['alpha']:6.2f} {res['r31']:10.2f} "
              f"{res['r31']/3477.23:10.4f}")


# Show g0_e(alpha) -- how does it change?
print("\n  g0_e(alpha) calibration curve:")
print(f"  {'alpha':>6s} {'g0_e':>10s} {'g0_mu':>10s} {'gc':>10s} {'g0_e/gc':>10s}")
for res in results:
    if res and res.get('g0_e'):
        print(f"  {res['alpha']:6.2f} {res['g0_e']:10.6f} {res['g0_mu']:10.6f} "
              f"{res['gc']:10.6f} {res['g0_e']/res['gc']:10.6f}")


# Summary
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print("""
  In canonical Form A:
  - F(alpha) = r_21 by construction (calibration).
  - The nontrivial quantity is r31(alpha) from Koide Q_K = 3/2.
  - r31 varies with alpha because the ODE dynamics change.
  - The PDG value r31 = 3477 selects a specific alpha.
  - This is the CANONICAL analog of the Form B alpha* determination.
""")
