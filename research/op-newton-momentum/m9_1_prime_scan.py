"""
m9_1_prime_scan.py -- M9.1' numerical scan of (alpha, p) parameter space.

Goal: verify the analytical master formula for TGP statics+PPN:
    c_2(alpha)      = -alpha/2                               (from Phi-EOM)
    beta_PPN(p, c_2) = (p-1)/p + 2 c_2 / p                   (from metric f=psi^p, h=psi^-p)

Equation solved (general alpha; vacuum beta=0; spherical):
    nabla^2 eps + alpha (grad eps)^2 / (1+eps) = -q rho(r)

State: v(r) = r * eps(r).
    dv/dr  = v'
    dv'/dr = r * [- alpha (eps')^2 / (1+eps) - q rho(r)]
where eps = v/r, eps' = v'/r - v/r^2.

We do NOT actually evaluate the metric here; we extract c_2 from the
asymptotic 1/r^2 coefficient of eps(r) and compare to -alpha/2.
beta_PPN is then computed analytically from (p, c_2). The only TGP
input that requires numerics is c_2 vs alpha.

Reference forms (substrate budget f*h = 1, gamma_PPN = 1 automatic):
  p = -1     "boxed" sek08c eq:metric-full-derived (f=1/psi, h=psi)
  p = -1/2   sek08c thm:antipodal-uniqueness (f=psi^-1/2, h=psi^1/2)

GR target: beta_PPN = 1.

Date: 2026-04-25
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import solve_bvp
from scipy.special import erf
from scipy.optimize import curve_fit

from m9_1_static import gaussian_source


def make_rhs_alpha(alpha: float, q: float, rho_func):
    """Vacuum (beta=0) Phi-EOM with general kinetic coefficient alpha."""

    def rhs(r, y):
        r_safe = np.maximum(r, 1.0e-9)
        eps = y[0] / r_safe
        eps_p = y[1] / r_safe - y[0] / r_safe ** 2
        rho = rho_func(r)
        ope = 1.0 + eps
        kin_nl = alpha * (eps_p ** 2) / ope
        dv_dr = y[1]
        dvp_dr = r * (-kin_nl - q * rho)
        return np.vstack([dv_dr, dvp_dr])

    return rhs


def bc_static(ya, yb):
    return np.array([ya[0], yb[1]])


def solve_alpha(alpha, q, rho_func, M_hint=1.0, sigma_hint=1.0,
                r_max=800.0, n_pts=5000):
    r_grid = np.linspace(1.0e-3, r_max, n_pts)
    A_g = q * M_hint / (4.0 * np.pi)
    v_g = A_g * erf(r_grid / (sigma_hint * np.sqrt(2.0)))
    vp_g = (A_g * np.sqrt(2.0 / np.pi) / sigma_hint
            * np.exp(-(r_grid ** 2) / (2.0 * sigma_hint ** 2)))
    y0 = np.vstack([v_g, vp_g])
    rhs = make_rhs_alpha(alpha, q, rho_func)
    sol = solve_bvp(rhs, bc_static, r_grid, y0,
                    tol=1.0e-9, max_nodes=100000, verbose=0)
    if not sol.success:
        raise RuntimeError(f"BVP failed (alpha={alpha}): {sol.message}")
    return sol


def fit_c2(sol, sigma=1.0, r_lo=8.0, r_hi=80.0, n=600):
    r_test = np.linspace(r_lo * sigma, r_hi * sigma, n)
    eps_num = sol.sol(r_test)[0] / r_test

    def ansatz(r, a1, a2, a3, a4, a5):
        return a1/r + a2/r**2 + a3/r**3 + a4/r**4 + a5/r**5

    popt, _ = curve_fit(ansatz, r_test, eps_num, p0=[0.1, 0.0, 0.0, 0.0, 0.0])
    a1, a2 = popt[0], popt[1]
    return a1, a2, a2 / a1**2


def beta_ppn(p, c_2):
    """Master formula: beta_PPN = (p-1)/p + 2 c_2 / p, given f=psi^p and f*h=1."""
    return (p - 1) / p + 2.0 * c_2 / p


def main():
    print("=" * 72)
    print("M9.1'  alpha-scan: verify c_2(alpha) = -alpha/2")
    print("=" * 72)
    print()
    print("Setup: M=q=sigma=1 (so A_0 = 1/(4 pi) ~= 0.0796); R_max=800, n_pts=5000")
    print()
    print(f"{'alpha':>6}  {'a_1 (A_eff)':>13}  {'a_2':>14}  {'c_2 numeric':>13}  "
          f"{'c_2 analytic':>13}  {'rel.err':>10}")
    print("-" * 72)

    sigma, M, q = 1.0, 1.0, 1.0
    rho = gaussian_source(M, sigma)

    rows = []
    for alpha in [0.0, 0.5, 1.0, 2.0, 3.0]:
        try:
            sol = solve_alpha(alpha, q, rho, M_hint=M, sigma_hint=sigma,
                              r_max=800.0, n_pts=5000)
            a1, a2, c2 = fit_c2(sol, sigma=sigma)
            c2_th = -alpha / 2.0
            err = (c2 - c2_th) / abs(c2_th) if c2_th != 0 else abs(c2)
            print(f"{alpha:6.2f}  {a1:13.7f}  {a2:+14.5e}  {c2:+13.5f}  "
                  f"{c2_th:+13.5f}  {err:+10.3e}")
            rows.append((alpha, a1, a2, c2, c2_th))
        except Exception as e:
            print(f"{alpha:6.2f}  -- failed: {e}")

    print()
    print("=" * 72)
    print("Master formula:  beta_PPN(p, alpha) = (p-1)/p - alpha/p")
    print("=" * 72)
    print()
    print(f"{'p':>6}  {'metric form':>26}  "
          f"{'beta_PPN(alpha=2)':>17}  {'with c_2 numeric':>18}")
    print("-" * 72)

    # numeric c_2 at alpha=2
    c2_at_2 = next((r[3] for r in rows if r[0] == 2.0), None)

    for p, label in [(-1.0,  "f=1/psi (boxed sek08c)"),
                     (-0.5,  "f=psi^-1/2 (antipodal)"),
                     (-2.0,  "f=1/psi^2 (DET=-c^2/psi^3)"),
                     (-1.0/3, "f=psi^-1/3")]:
        b_th = beta_ppn(p, -1.0)        # c_2 = -1 (alpha=2 analytic)
        b_num = beta_ppn(p, c2_at_2) if c2_at_2 is not None else float('nan')
        print(f"{p:+6.3f}  {label:>26}  {b_th:+17.4f}  {b_num:+18.4f}")

    print()
    print("Target GR: beta_PPN = 1.0 -- requires c_2 = +1/2, i.e. alpha = -1.")
    print("Excluded by N0-4 (K(0)=0 needs alpha > 0) and by stability.")
    print()
    print("Additionally, alpha = -1 corresponds to K(varphi) = K_0/varphi^2,")
    print("singular at varphi=0 -- contradicts substrate vacuum condition.")


if __name__ == "__main__":
    main()
