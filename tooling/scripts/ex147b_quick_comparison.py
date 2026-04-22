"""
ex147b_quick_comparison.py
==========================
Quick comparison: continuum f(g)=1+4ln(g) vs substrate K_sub=g^2
for phi-FP and r_21.

Key question: does the substrate ODE preserve r_21 = 206.768?
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

def ode_continuum(r, y):
    g, gp = y
    g = max(g, 1e-30)
    f = 1 + 4 * np.log(g)
    fp = 4.0 / g
    Vp = g**2 * (1 - g)
    if r < 1e-12:
        return [gp, -Vp / max(abs(f), 1e-15) / 3.0]
    if abs(f) < 1e-15:
        return [gp, 0.0]
    gpp = (Vp - fp * gp**2 / 2 - 2 * gp / r * f) / f
    return [gp, gpp]

def ode_substrate(r, y):
    g, gp = y
    g = max(g, 1e-30)
    Vp = g**2 * (1 - g)
    if r < 1e-12:
        return [gp, -Vp / g**2 / 3.0]
    gpp = (Vp - g * gp**2 - 2 * gp / r * g**2) / g**2
    return [gp, gpp]

def get_Atail(r, g, r_min=25, r_max=65):
    mask = (r >= r_min) & (r <= r_max) & (r > 0.1)
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
print("ex147b: Continuum vs Substrate ODE comparison")
print("=" * 60)

# Scan r_21 vs g0 for both descriptions
g0_vals = np.linspace(1.15, 1.45, 25)

print("\n  g0       r21_cont     r21_sub")
print("-" * 45)
for g0 in g0_vals:
    r21_c = compute_r21(g0, ode_continuum)
    r21_s = compute_r21(g0, ode_substrate)
    rc = f"{r21_c:.1f}" if r21_c else "FAIL"
    rs = f"{r21_s:.1f}" if r21_s else "FAIL"
    print(f"  {g0:.3f}    {rc:>10s}    {rs:>10s}")

# Try to find phi-FP for continuum (g0* ~ 1.249)
print("\n--- Fine scan near g0* = 1.249 (continuum) ---")
g0_fine = np.linspace(1.24, 1.26, 11)
for g0 in g0_fine:
    r21 = compute_r21(g0, ode_continuum)
    if r21:
        print(f"  g0={g0:.4f}  r21={r21:.2f}  delta={abs(r21-206.768)/206.768*100:.2f}%")
    else:
        print(f"  g0={g0:.4f}  FAIL")

print("\n--- Fine scan near g0* = 1.249 (substrate) ---")
for g0 in g0_fine:
    r21 = compute_r21(g0, ode_substrate)
    if r21:
        print(f"  g0={g0:.4f}  r21={r21:.2f}  delta={abs(r21-206.768)/206.768*100:.2f}%")
    else:
        print(f"  g0={g0:.4f}  FAIL")

print("\n[DONE]")
