"""
ex147c_correct_ode.py
=====================
Verify phi-FP using the CORRECT TGP soliton ODE from the paper:

  f(g)*g'' + (2/r)*g' = V'(g)     [eq:ghost-ode, eq:R-soliton-ode]

NOT the full Euler-Lagrange equation (which has extra f'(g')^2/2 and f(g)*(2/r)g' terms).

Then compare with substrate version:
  K_sub(g)*g'' + (2/r)*g' = V'(g)   with K_sub = g^2
"""
import numpy as np
from scipy.integrate import solve_ivp

PHI = (1 + np.sqrt(5)) / 2

def shoot(g0, ode_func, r_max=80):
    r_span = (1e-10, r_max)
    r_eval = np.linspace(1e-10, r_max, 30000)
    sol = solve_ivp(ode_func, r_span, [g0, 0.0],
                    method='Radau', rtol=1e-10, atol=1e-12,
                    t_eval=r_eval, max_step=0.3)
    if sol.success:
        return sol.t, sol.y[0]
    return None, None

def ode_paper(r, y):
    """TGP soliton ODE as in paper: f(g)g'' + (2/r)g' = V'(g)"""
    g, gp = y
    g = max(g, 1e-30)
    f = 1 + 4 * np.log(g)
    Vp = g**2 * (1 - g)
    if r < 1e-12:
        if abs(f) < 1e-15: f = 1e-15
        return [gp, -Vp / f / 3.0]
    if abs(f) < 1e-15:
        return [gp, 0.0]
    gpp = (Vp - 2 * gp / r) / f
    return [gp, gpp]

def ode_substrate_paper(r, y):
    """Substrate version: g^2 * g'' + (2/r)*g' = V'(g)"""
    g, gp = y
    g = max(g, 1e-30)
    Vp = g**2 * (1 - g)
    if r < 1e-12:
        return [gp, -Vp / g**2 / 3.0]
    gpp = (Vp - 2 * gp / r) / g**2
    return [gp, gpp]

def get_Atail(r, g, r_min=25, r_max=65):
    mask = (r >= r_min) & (r <= r_max)
    rr = r[mask]
    uu = (g[mask] - 1) * rr
    X = np.column_stack([np.sin(rr), np.cos(rr)])
    coeffs = np.linalg.lstsq(X, uu, rcond=None)[0]
    return np.sqrt(coeffs[0]**2 + coeffs[1]**2)

def compute_r21(g0, ode_func):
    r1, g1 = shoot(g0, ode_func)
    if r1 is None: return None
    r2, g2 = shoot(PHI * g0, ode_func)
    if r2 is None: return None
    A1 = get_Atail(r1, g1)
    A2 = get_Atail(r2, g2)
    if A1 < 1e-20: return None
    return (A2 / A1) ** 4

print("=" * 60)
print("ex147c: Paper ODE vs Substrate ODE")
print("  Paper:     f(g)g'' + (2/r)g' = V'(g)")
print("  Substrate: g^2*g'' + (2/r)g' = V'(g)")
print("=" * 60)

# Scan
g0_vals = np.linspace(1.15, 1.45, 25)
print("\n  g0       r21_paper    r21_sub")
print("-" * 45)
for g0 in g0_vals:
    r21_p = compute_r21(g0, ode_paper)
    r21_s = compute_r21(g0, ode_substrate_paper)
    rp = f"{r21_p:.1f}" if r21_p else "FAIL"
    rs = f"{r21_s:.1f}" if r21_s else "FAIL"
    mark = " <-- PDG" if r21_p and abs(r21_p - 206.768) < 20 else ""
    print(f"  {g0:.3f}    {rp:>10s}    {rs:>10s}{mark}")

# Fine scan for paper ODE
print("\n--- Fine scan paper ODE near g0=1.249 ---")
g0_fine = np.linspace(1.24, 1.26, 11)
for g0 in g0_fine:
    r21 = compute_r21(g0, ode_paper)
    if r21:
        mark = " ***" if abs(r21 - 206.768) / 206.768 < 0.01 else ""
        print(f"  g0={g0:.4f}  r21={r21:.2f}  delta={abs(r21-206.768)/206.768*100:.3f}%{mark}")
    else:
        print(f"  g0={g0:.4f}  FAIL")

# Fine scan for substrate ODE
print("\n--- Fine scan substrate ODE near g0=1.249 ---")
for g0 in g0_fine:
    r21 = compute_r21(g0, ode_substrate_paper)
    if r21:
        print(f"  g0={g0:.4f}  r21={r21:.2f}  delta={abs(r21-206.768)/206.768*100:.3f}%")
    else:
        print(f"  g0={g0:.4f}  FAIL")

print("\n[DONE]")
