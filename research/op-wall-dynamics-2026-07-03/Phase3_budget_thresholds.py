#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
op-wall-dynamics (Phase 3) -- budget as common source of both thresholds.

LOCK: Phase0_balance.md sec. 2 Phase 3 (W3a analytic+numeric, W3b
descriptive only, SPECULATIVE, zero claims).

W3a: define E_core[g0] and B_core[g0] (energy and created-space budget
in the profile core) from ODE profiles; check whether g* = e^(-1/4)
(lower, kinetic wall) and g_crit = 8/5 (upper, H7/H8 in
tgp_master_consistency_v47.py) correspond to extrema / loss of
monotonicity of the same budget quantity (sympy + numerics).

Definitions used (fixed here, before the scan):
  core      = ball r < r_core, r_core = first node of g_eq - 1
  B_core    = int_0^{r_core} (g - 1) r^2 dr          (w = 1)
  E_core    = int_0^{r_core} [ (1/2) f(g) g'^2 + W(g) - W(1) ] r^2 dr
  W(g)      = g^3/3 - g^4/4   (W' = g^2(1-g)),  f(g) = 1 + 4 ln g
  touch     = profile reaches wall guard g* + 0.005 (a3d event)

Analytic tool (derived exactly below with sympy): the "mechanical
energy" H(r) = (1/2) f(g) g'^2 - W(g) satisfies
  dH/dr = -(2/r) f(g) g'^2  (monotonically decreasing while f > 0),
hence a NECESSARY condition for wall contact is W(g0) <= W(g*),
i.e. g0 >= g0_nec (frictionless bound).  The numeric scan then finds
the ACTUAL contact onset g0_wall and compares it with 8/5 = 1.6.

W3b (SPECULATIVE, descriptive): table of (g0, bounces, N_localized)
along the scan -- correlation documentation only, no promotion.

Units gamma=beta=K_geo=1.  Localized mode: lam < -1 - 1e-3.
"""
import numpy as np
import sympy as sp
from scipy.linalg import eigh_tridiagonal
from scipy.integrate import solve_ivp

PHI = (1 + 5 ** 0.5) / 2
G0_E = 1.24915
G0_MU = PHI * G0_E
G0_TAU = 3.18912
ALPHA = 2.0
G_GHOST = np.exp(-1 / (2 * ALPHA))
LOC_EDGE = -1.0 - 1e-3

f_log   = lambda g: 1 + 2 * ALPHA * np.log(g)
fp_log  = lambda g: 2 * ALPHA / g
fpp_log = lambda g: -2 * ALPHA / g ** 2
W1_fun  = lambda g: g ** 2 * (1 - g)
W2_fun  = lambda g: 2 * g - 3 * g ** 2
W_fun   = lambda g: g ** 3 / 3 - g ** 4 / 4

# =====================================================================
# W3a -- analytic part (sympy, exact)
# =====================================================================
print("=" * 78)
print("Phase 3: budget thresholds.  W3a analytic (sympy) + numeric scan.")
print("=" * 78)

g_s, gp_s, r_s = sp.symbols('g gp r', positive=True)
alpha_s = sp.Integer(2)
f_s = 1 + 2 * alpha_s * sp.log(g_s)
W_s = g_s ** 3 / 3 - g_s ** 4 / 4
Wp_s = sp.diff(W_s, g_s)

# (a) exact wall location and exact potential values
g_star = sp.exp(sp.Rational(-1, 4))
print("\n[sympy] f(g*) with g* = e^(-1/4):", sp.simplify(f_s.subs(g_s, g_star)))
W_gstar = sp.simplify(W_s.subs(g_s, g_star))
print("[sympy] W(g*) = e^(-3/4)/3 - e^(-1)/4 =", sp.nsimplify(W_gstar),
      "= %.10f" % float(W_gstar))
print("[sympy] W(1)  = 1/12 = %.10f" % float(W_s.subs(g_s, 1)))

# (b) exact monotonicity identity dH/dr = -(2/r) f g'^2 on the EL ODE
#     EL: f g'' + (1/2) f' g'^2 + (2/r) f g' - W' = 0
gpp_from_ode = (Wp_s - sp.Rational(1, 2) * sp.diff(f_s, g_s) * gp_s ** 2
                - 2 * f_s * gp_s / r_s) / f_s
H_s = sp.Rational(1, 2) * f_s * gp_s ** 2 - W_s
dH = (sp.diff(H_s, g_s) * gp_s
      + sp.diff(H_s, gp_s) * gpp_from_ode)      # dH/dr along solutions
identity = sp.simplify(dH + 2 * f_s * gp_s ** 2 / r_s)
print("[sympy] dH/dr + (2/r) f g'^2 along EL solutions =", identity,
      " (exact PASS)" if identity == 0 else " (FAIL)")

# (c) frictionless necessary bound for wall contact: W(g0) <= W(g*)
g0_nec = sp.nsolve(W_s - W_gstar, g_s, 1.2)
print("[sympy] necessary contact bound g0_nec (W(g0)=W(g*), g0>1):"
      " g0 >= %.6f" % float(g0_nec))

# (d) H7 formula (registered, tgp_master_consistency_v47.py):
g0_crit_H7 = sp.Rational(2 * (2 + 2), 2 * (2 + 2) - 3)
print("[H7]    g0_crit = 2(alpha+2)/(2(alpha+2)-d) = %s = %.6f"
      % (g0_crit_H7, float(g0_crit_H7)))

# =====================================================================
# numeric machinery (CP-7 profile, verbatim conventions)
# =====================================================================
def soliton_profile(g0, R=60.0, N=3000):
    def rhs(r_, y):
        g, gp = y
        g = max(g, G_GHOST + 1e-6)
        fk = f_log(g)
        drv = W1_fun(g) - (ALPHA / g) * gp ** 2
        if r_ < 1e-10:
            return [gp, drv / (3 * fk)]
        return [gp, (drv - fk * 2 * gp / r_) / fk]

    def wall(r_, y):
        return y[0] - (G_GHOST + 0.005)
    wall.terminal = True
    wall.direction = -1

    h = R / N
    r_grid = (np.arange(N) + 0.5) * h
    g_arr = np.empty(N); gp_arr = np.empty(N)
    r0, y = 1e-6, [g0, 0.0]
    bounces = 0; idx0 = 0
    for b in range(30):
        sol = solve_ivp(rhs, [r0, R + h], y, method='DOP853', max_step=0.02,
                        rtol=1e-10, atol=1e-13, events=[wall],
                        dense_output=True)
        seg_end = sol.t[-1]
        sel = (r_grid >= r0) & (r_grid <= seg_end)
        if np.any(sel):
            vals = sol.sol(r_grid[sel])
            g_arr[idx0:idx0 + sel.sum()] = vals[0]
            gp_arr[idx0:idx0 + sel.sum()] = vals[1]
            idx0 += sel.sum()
        if sol.t_events[0].size > 0 and b < 29:
            r_b = float(sol.t_events[0][0])
            gp_b = float(sol.y_events[0][0, 1])
            r0 = r_b + 1e-6
            y = [G_GHOST + 0.005 + 1e-5, -gp_b]
            bounces += 1
        else:
            break
    if idx0 < N:
        r_grid = r_grid[:idx0]; g_arr = g_arr[:idx0]; gp_arr = gp_arr[:idx0]
    d2 = np.array([rhs(r_grid[j], [g_arr[j], gp_arr[j]])[1]
                   for j in range(len(r_grid))])
    return r_grid, g_arr, gp_arr, d2, bounces, float(np.min(g_arr))


def touches_wall(g0, R=60.0):
    """True if the integration triggers the wall event at least once."""
    _, _, _, _, nb, gmin = soliton_profile(g0, R=R, N=300)
    return nb >= 1


def core_quantities(r, g, gp):
    s = g - 1.0
    r_core = None
    for j in range(len(r) - 1):
        if s[j] == 0.0 or s[j] * s[j + 1] < 0:
            r_core = r[j] + (r[j + 1] - r[j]) * abs(s[j]) / (abs(s[j]) + abs(s[j + 1]))
            break
    if r_core is None:
        return None, None, None
    m = r <= r_core
    B = float(np.trapezoid((g[m] - 1) * r[m] ** 2, r[m]))
    dens = 0.5 * f_log(g[m]) * gp[m] ** 2 + (W_fun(g[m]) - W_fun(1.0))
    E = float(np.trapezoid(dens * r[m] ** 2, r[m]))
    return r_core, B, E


def build_tridiag(r, g, gp, d2, l=0):
    N = len(r)
    h = r[1] - r[0]
    F_nodes = f_log(g)
    F_mid = f_log(0.5 * (g[:-1] + g[1:]))
    Q = W2_fun(g) - 0.5 * fpp_log(g) * gp ** 2 - fp_log(g) * (d2 + 2 * gp / r)
    r_mid = 0.5 * (r[:-1] + r[1:])
    a = r_mid ** 2 * F_mid
    diag = np.empty(N)
    diag[0] = a[0] / h ** 2
    diag[1:-1] = (a[:-1] + a[1:]) / h ** 2
    diag[-1] = (a[-1] + r[-1] ** 2 * F_nodes[-1]) / h ** 2
    diag += r ** 2 * (Q + F_nodes * l * (l + 1) / r ** 2)
    off = -a / h ** 2
    return diag / r ** 2, off / (r[:-1] * r[1:])


# =====================================================================
# numeric part 1: actual wall-contact onset g0_wall vs 8/5
# =====================================================================
print("\n[numeric] wall-contact onset g0_wall (bisection, guard g*+0.005):")
lo, hi = 1.02, 2.03
# ensure bracket
assert not touches_wall(lo), "lower bracket touches wall"
assert touches_wall(hi), "upper bracket does not touch wall"
while hi - lo > 1e-5:
    mid = 0.5 * (lo + hi)
    if touches_wall(mid):
        hi = mid
    else:
        lo = mid
g0_wall = 0.5 * (lo + hi)
print("  g0_wall = %.5f" % g0_wall)
print("  8/5     = 1.60000   -> |g0_wall - 8/5| = %.5f (%.2f%%)"
      % (abs(g0_wall - 1.6), 100 * abs(g0_wall - 1.6) / 1.6))
print("  g0_nec  = %.5f (frictionless necessary bound)" % float(g0_nec))
print("  reference points: g0_e=%.5f (no contact), g0_mu=%.5f, g0_tau=%.5f"
      % (G0_E, G0_MU, G0_TAU))

# =====================================================================
# numeric part 2: B_core(g0), E_core(g0) scan + monotonicity
# =====================================================================
print("\n[numeric] scan g0 in [1.04, 3.40], step 0.02 (N=3000, R=60):")
print("  %-6s %-4s %-8s %-8s %-12s %-12s" %
      ("g0", "nb", "min_g", "r_core", "B_core", "E_core"))
g0_grid = np.round(np.arange(1.04, 3.4001, 0.02), 4)
scan = []
for g0 in g0_grid:
    r, g, gp, d2, nb, gmin = soliton_profile(g0)
    r_core, B, E = core_quantities(r, g, gp)
    scan.append((g0, nb, gmin, r_core, B, E))
    print("  %-6.2f %-4d %-8.4f %-8s %-12s %-12s"
          % (g0, nb, gmin,
             "%.4f" % r_core if r_core else "None",
             "%.6f" % B if B is not None else "n/a",
             "%.6f" % E if E is not None else "n/a"))

print("\n[numeric] monotonicity / extrema of B_core, E_core (sign of dX/dg0):")
gs = np.array([s[0] for s in scan])
Bs = np.array([s[4] if s[4] is not None else np.nan for s in scan])
Es = np.array([s[5] if s[5] is not None else np.nan for s in scan])
nbs = np.array([s[1] for s in scan])
for name, X in (("B_core", Bs), ("E_core", Es)):
    dX = np.diff(X) / np.diff(gs)
    flips = [0.5 * (gs[i] + gs[i + 1])
             for i in range(1, len(dX))
             if np.isfinite(dX[i - 1]) and np.isfinite(dX[i])
             and dX[i - 1] * dX[i] < 0]
    print("  %s: sign changes of dX/dg0 at g0 ~ %s"
          % (name, ["%.2f" % v for v in flips] if flips else "NONE"))
bounce_up = [0.5 * (gs[i] + gs[i + 1]) for i in range(len(nbs) - 1)
             if nbs[i + 1] > nbs[i]]
print("  bounce-count increments at g0 ~ %s" % ["%.2f" % v for v in bounce_up])

# =====================================================================
# W3b (SPECULATIVE, descriptive only): saddle index vs g0 / bounces
# =====================================================================
print("\n[W3b - SPECULATIVE, correlation documentation ONLY, no claims]")
print("  N_localized (l=0, N=2000, R=60) along the g0 scan:")
print("  %-6s %-4s %-6s  localized modes" % ("g0", "nb", "N_loc"))
for g0 in np.round(np.arange(1.1, 3.4001, 0.1), 4):
    r, g, gp, d2, nb, gmin = soliton_profile(g0, N=2000)
    d, e = build_tridiag(r, g, gp, d2, l=0)
    vals = eigh_tridiagonal(d, e, select='i', select_range=(0, 59),
                            eigvals_only=True)
    loc = [float(v) for v in vals if v < LOC_EDGE]
    print("  %-6.1f %-4d %-6d %s" % (g0, nb, len(loc),
                                     ["%.3f" % v for v in loc]))

print("\n" + "=" * 78)
print("W3 verdict inputs complete (interpretation -> README per LOCK).")
print("=" * 78)
