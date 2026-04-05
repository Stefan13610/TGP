#!/usr/bin/env python3
"""Quick phi-FP bisection for g0*."""
import sys, os
os.environ['PYTHONIOENCODING'] = 'utf-8'
if sys.stdout.encoding != 'utf-8':
    sys.stdout.reconfigure(encoding='utf-8')

import numpy as np
from scipy.integrate import solve_ivp

PHI = (1 + np.sqrt(5)) / 2
ALPHA = 2
GG = np.exp(-1 / (2 * ALPHA))

def fk(g):
    return 1 + 2 * ALPHA * np.log(g) if g > 0 else -1e30

def Vp(g):
    return g**2 * (1 - g)

def ss(g0, rm=300, np_=25000):
    fg0 = fk(g0)
    c2 = Vp(g0) / (3 * fg0)
    rs = 0.01
    def rhs(r, y):
        g, p = y
        if g <= 1e-15:
            return [p, 0]
        fg = fk(g)
        if abs(fg) < 1e-10:
            return [p, 0]
        if r < 1e-10:
            return [p, Vp(g) / fg / 3]
        return [p, (Vp(g) - 2/r * p) / fg]
    def ev(r, y):
        return 100 - abs(y[0])
    ev.terminal = True
    s = solve_ivp(rhs, [rs, rm], [g0 + c2*rs**2, 2*c2*rs],
                  method='RK45', rtol=1e-11, atol=1e-13,
                  max_step=0.05, events=[ev], dense_output=True)
    r = np.linspace(rs, min(s.t[-1], rm), np_)
    return r, s.sol(r)[0]

def At(r, g):
    m = (r >= 120) & (r <= 260)
    rf = r[m]
    tl = (g[m] - 1) * r[m]
    if len(rf) < 10:
        return np.nan
    A = np.column_stack([np.cos(rf), np.sin(rf)])
    B, C = np.linalg.lstsq(A, tl, rcond=None)[0]
    return np.sqrt(B**2 + C**2)

def R_func(g0):
    r1, g1 = ss(g0)
    r2, g2 = ss(PHI * g0)
    if r1[-1] < 250 or r2[-1] < 250:
        return np.nan
    Ae = At(r1, g1)
    Am = At(r2, g2)
    if Ae < 1e-15:
        return np.nan
    return (Am / Ae)**4

# Bisekcja: R rosnie z g0
target = 206.768
lo, hi = 0.898, 0.901
for _ in range(30):
    mid = (lo + hi) / 2
    R = R_func(mid)
    if R is None or np.isnan(R):
        break
    if R < target:
        lo = mid
    else:
        hi = mid

g0_star = (lo + hi) / 2
R_star = R_func(g0_star)
r1, g1 = ss(g0_star)
r2, g2 = ss(PHI * g0_star)
Ae = At(r1, g1)
Am = At(r2, g2)

print("=" * 60)
print("phi-FP WYNIK (ODE dodatekJ: f(g)*g'' + 2/r*g' = V'(g))")
print("=" * 60)
print(f"  g_0* (elektron)  = {g0_star:.10f}")
print(f"  g_0^mu (mion)    = {PHI * g0_star:.10f}  [= phi * g_0*]")
print(f"  A_tail(e)        = {Ae:.10f}")
print(f"  A_tail(mu)       = {Am:.10f}")
print(f"  r_21 = (A_mu/A_e)^4 = {R_star:.6f}")
print(f"  r_21 PDG             = {target}")
print(f"  Odchylenie           = {abs(R_star - target) / target * 100:.6f}%")
print(f"  Ghost g*             = {GG:.6f}")
print(f"  g_0* / g*            = {g0_star / GG:.6f}")
print()

# Tau: phi^2 * g0*
g0_tau = PHI**2 * g0_star
print(f"  g_0^tau = phi^2*g0* = {g0_tau:.6f}")
try:
    rt, gt = ss(g0_tau)
    if rt[-1] >= 250:
        A_tau = At(rt, gt)
        r31 = (A_tau / Ae)**4
        print(f"  A_tail(tau)      = {A_tau:.8f}")
        print(f"  r_31 (phi^2)     = {r31:.2f}")
        print(f"  r_31 PDG         = 3477.48")
        print(f"  Odchylenie       = {abs(r31 - 3477.48) / 3477.48 * 100:.2f}%")
    else:
        print(f"  [ODE diverguje przy r={rt[-1]:.0f}]")
except Exception as e:
    print(f"  [Blad: {e}]")

# Tau: 2 * g0*
g0_tau2 = 2 * g0_star
print(f"\n  g_0^tau = 2*g0*   = {g0_tau2:.6f}")
try:
    rt2, gt2 = ss(g0_tau2)
    if rt2[-1] >= 250:
        A_tau2 = At(rt2, gt2)
        r31_h = (A_tau2 / Ae)**4
        print(f"  A_tail(tau)      = {A_tau2:.8f}")
        print(f"  r_31 (harmon.)   = {r31_h:.2f}")
        print(f"  r_31 PDG         = 3477.48")
        print(f"  Odchylenie       = {abs(r31_h - 3477.48) / 3477.48 * 100:.2f}%")
    else:
        print(f"  [ODE diverguje przy r={rt2[-1]:.0f}]")
except Exception as e:
    print(f"  [Blad: {e}]")

# Koide r31
r21 = R_star
a = 1 + np.sqrt(r21)
disc = 6 * a**2 - 3 - 3 * r21
r31_K = (2 * a + np.sqrt(disc))**2
print(f"\n  r_31 (Koide)     = {r31_K:.2f}")
print(f"  r_31 PDG         = 3477.48")
print(f"  Odchylenie       = {abs(r31_K - 3477.48) / 3477.48 * 100:.4f}%")
Q = (1 + np.sqrt(r21) + np.sqrt(r31_K))**2 / (1 + r21 + r31_K)
print(f"  Q (Koide)        = {Q:.10f}")
print()
print("GOTOWE.")
