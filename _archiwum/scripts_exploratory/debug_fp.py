"""Debug: grid FP via ODE shooting from minimum, check nu."""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.linalg import eigvals

C_D = 1.0 / (6.0 * np.pi**2)
D_DIM = 3


def fp_ode(rho, y):
    """3rd-order ODE for LPA FP: y = [u, u', u'']."""
    u, up, u2v = y
    D = 1.0 + up + 2.0 * rho * u2v
    if abs(D) < 1e-9:
        D = 1e-9
    # Derived from differentiating the FP equation
    # u''' = [u''(rho*D^2 - 3*c3) - 2*u'*D^2] / (2*rho*c3)
    if rho < 1e-8:
        return [up, u2v, 0.0]
    u3 = (u2v * (rho * D**2 - 3.0 * C_D) - 2.0 * up * D**2) / (2.0 * rho * C_D)
    return [up, u2v, u3]


def shoot_from_kappa(kappa, u2_kappa, n_pts=150, rho_max=1.5):
    """
    Integrate FP ODE forward and backward from rho=kappa.
    Returns (rho_grid, u_grid, status_fwd, status_bwd).
    """
    D_kap = 1.0 + 2.0 * kappa * u2_kappa
    u_min = C_D / (3.0 * D_kap)
    y0    = [u_min, 0.0, u2_kappa]

    eps   = 1e-5
    # Forward: kappa -> rho_max
    sol_fwd = solve_ivp(
        fp_ode, [kappa, rho_max], y0,
        method="DOP853", rtol=1e-10, atol=1e-12,
        max_step=0.005, dense_output=True,
    )
    # Backward: kappa -> eps
    sol_bwd = solve_ivp(
        fp_ode, [kappa, eps], y0,
        method="DOP853", rtol=1e-10, atol=1e-12,
        max_step=0.001, dense_output=True,
    )
    # Build grid
    rho_lo = np.linspace(eps,   kappa,   n_pts // 2,  endpoint=False)
    rho_hi = np.linspace(kappa, rho_max, n_pts // 2 + 1)
    rho_all = np.concatenate([rho_lo, rho_hi])

    u_all = np.zeros_like(rho_all)
    for i, r in enumerate(rho_all):
        if r <= kappa:
            y = sol_bwd.sol(r)
        else:
            y = sol_fwd.sol(r)
        u_all[i] = y[0]

    return rho_all, u_all, sol_fwd.status, sol_bwd.status


def lpa_rhs(u, rho, eta=0.0):
    dx  = rho[1] - rho[0]
    du  = np.gradient(u, dx)
    d2u = np.gradient(du, dx)
    den = 1.0 + du + 2.0 * rho * d2u
    den = np.where(np.abs(den) < 1e-10, 1e-10, den)
    rhs = -D_DIM * u + (1.0 + eta / 2.0) * rho * du + C_D * (1.0 - eta / 5.0) / den
    rhs[0]  = rhs[1]
    rhs[-1] = rhs[-2]
    return rhs


def compute_nu(u, rho, eta=0.0):
    eps  = 1e-5
    n    = len(rho)
    M    = np.zeros((n, n))
    rhs0 = lpa_rhs(u, rho, eta)
    step = 4
    for j in range(0, n, step):
        up      = u.copy()
        up[j]  += eps
        M[:, j] = (lpa_rhs(up, rho, eta) - rhs0) / eps
    for j in range(n):
        if j % step != 0:
            j0 = (j // step) * step
            j1 = min(j0 + step, n - 1)
            a  = (j - j0) / max(j1 - j0, 1)
            M[:, j] = (1 - a) * M[:, j0] + a * M[:, j1]
    eigs  = np.sort(eigvals(M).real)[::-1]
    theta = eigs[0]
    return (1.0 / theta if theta > 1e-6 else float("nan")), eigs[:4]


# Scan u2 to find the one with smallest max|RHS| on the grid
print("Scanning u2 for grid-based FP near kappa=0.075...")
print("u2       kap_used  max_RHS   nu        status")
for u2 in [0.5, 0.8, 1.0, 1.2, 1.5, 2.0, 3.0, 5.0]:
    kappa_try = 0.075
    try:
        rho_g, u_g, sf, sb = shoot_from_kappa(kappa_try, u2, n_pts=100, rho_max=1.0)
        max_rhs = float(np.max(np.abs(lpa_rhs(u_g, rho_g))))
        nu, eigs = compute_nu(u_g, rho_g)
        print(f"u2={u2:.2f}  k={kappa_try:.4f}  res={max_rhs:.4e}  nu={nu:.4f}  sf={sf} sb={sb}")
    except Exception as e:
        print(f"u2={u2:.2f}  ERROR {e}")

print()
print("Now scan (kappa, u2) for smallest max|RHS|...")
best_res = 1e10
best_kappa, best_u2 = 0.075, 1.0
for u2 in [0.5, 1.0, 2.0]:
    for kappa in [0.05, 0.07, 0.075, 0.08, 0.09, 0.10]:
        try:
            rho_g, u_g, sf, sb = shoot_from_kappa(kappa, u2, n_pts=80, rho_max=0.8)
            max_rhs = float(np.max(np.abs(lpa_rhs(u_g, rho_g))))
            if max_rhs < best_res:
                best_res = max_rhs
                best_kappa, best_u2 = kappa, u2
            print(f"  k={kappa:.4f} u2={u2:.2f}  res={max_rhs:.4e} sf={sf} sb={sb}")
        except Exception as e:
            print(f"  k={kappa:.4f} u2={u2:.2f}  ERROR {e}")

print(f"Best: kappa={best_kappa}, u2={best_u2}, res={best_res:.4e}")
