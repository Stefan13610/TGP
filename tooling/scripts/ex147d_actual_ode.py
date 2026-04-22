"""
ex147d_actual_ode.py
====================
Comparison using the ACTUAL TGP soliton ODE from fermion_mass_spectrum.py:

  g'' + (2/r)g' + (alpha/g)(g')^2 = V'(g)
  alpha = 2, V'(g) = g^2 - g^3

Substrate analogue: replace (alpha/g) with K_sub'/(2*K_sub) where K_sub = g^2:
  K_sub'/(2*K_sub) = 2g/(2g^2) = 1/g
  So substrate ODE: g'' + (2/r)g' + (1/g)(g')^2 = V'(g)

Note: alpha=2 gives (alpha/g) = 2/g. The substrate K_sub=g^2 gives 1/g.
These are DIFFERENT coefficients => different physics.
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

def ode_actual_tgp(r, y):
    """ACTUAL TGP: g'' + (2/r)g' + (2/g)(g')^2 = g^2 - g^3"""
    g, gp = y
    g = max(g, 1e-30)
    Vp = g**2 - g**3
    if r < 1e-12:
        return [gp, Vp / 3.0]  # L'Hopital for (2/r)g' at r=0
    gpp = Vp - 2 * gp / r - (2.0 / g) * gp**2
    return [gp, gpp]

def ode_substrate(r, y):
    """Substrate K_sub=g^2: g'' + (2/r)g' + (1/g)(g')^2 = g^2 - g^3"""
    g, gp = y
    g = max(g, 1e-30)
    Vp = g**2 - g**3
    if r < 1e-12:
        return [gp, Vp / 3.0]
    gpp = Vp - 2 * gp / r - (1.0 / g) * gp**2
    return [gp, gpp]

def ode_substrate_exact(r, y):
    """Substrate exact: full E-L from S=int[g^2/2*(g')^2 + V(g)]r^2 dr
    g'' + (2/r)g' + (1/g)(g')^2 = V'(g)/g^2  [note: V' divided by K_sub!]
    """
    g, gp = y
    g = max(g, 1e-30)
    Vp = g**2 - g**3  # V'(g)
    if r < 1e-12:
        return [gp, (Vp / g**2) / 3.0]
    gpp = Vp / g**2 - 2 * gp / r - (1.0 / g) * gp**2
    return [gp, gpp]

def get_Atail(r, g, r_min=25, r_max=65):
    mask = (r >= r_min) & (r <= r_max)
    if np.sum(mask) < 10:
        return 0
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

print("=" * 65)
print("ex147d: Actual TGP ODE vs Substrate versions")
print("  TGP:   g'' + (2/r)g' + (2/g)(g')^2 = g^2 - g^3")
print("  Sub1:  g'' + (2/r)g' + (1/g)(g')^2 = g^2 - g^3")
print("  Sub2:  g'' + (2/r)g' + (1/g)(g')^2 = (g^2-g^3)/g^2")
print("=" * 65)

g0_vals = np.linspace(1.10, 1.45, 30)
print("\n  g0       r21_TGP      r21_sub1     r21_sub2")
print("-" * 58)
for g0 in g0_vals:
    r21_t = compute_r21(g0, ode_actual_tgp)
    r21_s1 = compute_r21(g0, ode_substrate)
    r21_s2 = compute_r21(g0, ode_substrate_exact)
    rt = f"{r21_t:.1f}" if r21_t else "FAIL"
    rs1 = f"{r21_s1:.1f}" if r21_s1 else "FAIL"
    rs2 = f"{r21_s2:.1f}" if r21_s2 else "FAIL"
    mark = ""
    if r21_t and abs(r21_t - 206.768) < 15:
        mark = " <-- PDG"
    print(f"  {g0:.3f}    {rt:>10s}    {rs1:>10s}    {rs2:>10s}{mark}")

# Fine scan near PDG value for TGP
print("\n--- Fine scan TGP near phi-FP ---")
# First find approximate location
best_g0 = None
best_diff = 1e10
for g0 in np.linspace(1.10, 1.35, 50):
    r21 = compute_r21(g0, ode_actual_tgp)
    if r21 and abs(r21 - 206.768) < best_diff:
        best_diff = abs(r21 - 206.768)
        best_g0 = g0

if best_g0:
    print(f"  Approximate g0* = {best_g0:.3f}")
    g0_fine = np.linspace(best_g0 - 0.02, best_g0 + 0.02, 21)
    for g0 in g0_fine:
        r21 = compute_r21(g0, ode_actual_tgp)
        if r21:
            mark = " ***" if abs(r21 - 206.768) / 206.768 < 0.005 else ""
            print(f"  g0={g0:.5f}  r21={r21:.2f}  delta={abs(r21-206.768)/206.768*100:.3f}%{mark}")

print("\n[DONE]")
