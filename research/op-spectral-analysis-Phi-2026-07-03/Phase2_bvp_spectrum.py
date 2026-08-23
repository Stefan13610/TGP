#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
op-spectral-analysis-Phi (CP-7 / L03) -- Phase 2: numerical BVP diagonalization.

Operator (Phase 1, verified exactly):
    L[v] = -(1/r^2) d/dr[ r^2 F(u0) dv/dr ] + [ Q(r) + F(u0) l(l+1)/r^2 ] v = lambda v
    Q    = W''(u0) - (1/2) F''(u0) u0'^2 - F'(u0) [u0'' + (2/r) u0']

Discretization: staggered grid r_j=(j+1/2)h, self-adjoint flux form,
natural regularity BC at 0 (r^2 F v' -> 0), Dirichlet at R.
Symmetrized tridiagonal via B^(-1/2) M B^(-1/2), B = diag(r_j^2).

Backgrounds (Phase0 LOCK): B1 vacuum F-A; B2 Yukawa F-A (linear A=0.01,0.1;
nonlinear Newton relaxation for target amplitudes 0.3, 1.0); B3 solitons
F-S (crown log-form alpha=2, wall bounces as in a3d) for g0 in {e, mu, tau};
B4 vacuum F-S (dualism control); B5 solitons F-S' (substrate alpha=1).

Criteria: C2, C3, C4(report), C5, C6a as locked. Units: gamma=beta=K_geo=1.
Tolerance for negative mode: lambda < -1e-6 confirmed on 3 grids.
"""
import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.integrate import solve_ivp

PHI = (1 + 5 ** 0.5) / 2
G0_E = 1.24915
G0_MU = PHI * G0_E
G0_TAU = 3.18912
ALPHA = 2.0
G_GHOST = np.exp(-1 / (2 * ALPHA))          # 0.7788 (F-S log form)
TOL_NEG = -1e-6

# ======================================================================
# formulations: (F, F', F'', W'') as callables
# ======================================================================
FORM = {
    'F-A':  dict(F=lambda u: u**4, Fp=lambda u: 4*u**3, Fpp=lambda u: 12*u**2,
                 W2=lambda u: 7*u**6 - 6*u**5),
    'F-S':  dict(F=lambda u: 1 + 4*np.log(u), Fp=lambda u: 4/u, Fpp=lambda u: -4/u**2,
                 W2=lambda u: 2*u - 3*u**2),
    'F-Sp': dict(F=lambda u: u**2, Fp=lambda u: 2*u, Fpp=lambda u: 2*np.ones_like(u),
                 W2=lambda u: 2*u - 3*u**2),
}

def assemble_and_solve(r, F_nodes, F_mid, Q, l=0, k_eigs=40):
    """Symmetric tridiagonal eigenproblem for L, staggered grid, Dirichlet at R."""
    N = len(r)
    h = r[1] - r[0]
    r_mid = 0.5 * (r[:-1] + r[1:])                     # interior midpoints
    a = r_mid**2 * F_mid                               # flux coeffs (N-1)
    # left flux at r=0 vanishes (r^2 -> 0): natural regularity BC
    diag = np.empty(N)
    diag[0] = a[0] / h**2
    diag[1:-1] = (a[:-1] + a[1:]) / h**2
    diag[-1] = (a[-1] + r[-1]**2 * F_nodes[-1]) / h**2  # Dirichlet ghost at R
    diag += r**2 * (Q + F_nodes * l * (l + 1) / r**2)
    off = -a / h**2
    # symmetrize with B = diag(r^2): C = B^(-1/2) M B^(-1/2)
    d = diag / r**2
    e = off / (r[:-1] * r[1:])
    k = min(k_eigs, N - 1)
    vals = eigh_tridiagonal(d, e, select='i', select_range=(0, k - 1),
                            eigvals_only=True)
    return vals

def report_spectrum(tag, vals, expect=None):
    nneg = int(np.sum(vals < TOL_NEG))
    print("  %-34s lam_min=%+.6f  lam_1=%+.6f  N_neg=%d%s"
          % (tag, vals[0], vals[1] if len(vals) > 1 else np.nan, nneg,
             ("  [expect %s]" % expect) if expect else ""))
    return nneg, vals[0]

# ======================================================================
# B1 vacuum F-A ; B4 vacuum F-S
# ======================================================================
def vacuum_case(form_key, Qval, R=60.0, N=4000, l=0):
    h = R / N
    r = (np.arange(N) + 0.5) * h
    ones = np.ones_like(r)
    F_nodes = FORM[form_key]['F'](ones)
    F_mid = FORM[form_key]['F'](np.ones(N - 1))
    Q = Qval * np.ones_like(r)
    return r, assemble_and_solve(r, F_nodes, F_mid, Q, l=l)

# ======================================================================
# B2 Yukawa F-A: linear profile / nonlinear Newton relaxation
# ======================================================================
def yukawa_linear_background(A, R, N):
    h = R / N
    r = (np.arange(N) + 0.5) * h
    ex = np.exp(-r)
    psi = 1 + A * ex / r
    dpsi = -A * ex * (1 + r) / r**2
    d2psi = A * ex * (r**2 + 2*r + 2) / r**3
    return r, psi, dpsi, d2psi

def newton_background(S, R, N, sig=0.5, maxit=200):
    """Damped Newton for: psi'' + (2/r)psi' + (2/psi)psi'^2 + psi^2 - psi^3 = -rho.
    rho(r) = S exp(-(r/sig)^2)  (q~ folded into S; sign per rem:sign-convention).
    Staggered grid, psi'(0)=0 (mirror), psi(R)=1 (Dirichlet ghost)."""
    h = R / N
    r = (np.arange(N) + 0.5) * h
    rho = S * np.exp(-(r / sig) ** 2)
    psi = 1 + 0.1 * np.exp(-(r / sig) ** 2)            # init guess
    for it in range(maxit):
        psim = np.empty(N + 2)
        psim[1:-1] = psi
        psim[0] = psi[0]                               # mirror (psi'(0)=0)
        psim[-1] = 2.0 - psi[-1]                       # Dirichlet psi(R)=1 ghost
        d1 = (psim[2:] - psim[:-2]) / (2 * h)
        d2 = (psim[2:] - 2 * psim[1:-1] + psim[:-2]) / h**2
        Res = d2 + 2 * d1 / r + 2 * d1**2 / psi + psi**2 - psi**3 + rho
        if np.max(np.abs(Res)) < 1e-12:
            break
        # tridiagonal Jacobian
        Jd = (-2 / h**2) - 2 * d1**2 / psi**2 + 2 * psi - 3 * psi**2
        Ju = (1 / h**2) + (1 / (2 * h)) * (2 / r + 4 * d1 / psi)
        Jl = (1 / h**2) - (1 / (2 * h)) * (2 / r + 4 * d1 / psi)
        Jd0 = Jd.copy()
        Jd0[0] += Jl[0]                                # mirror BC folds lower into diag
        Jd0[-1] -= Ju[-1]                              # ghost: psim[-1]=2-psi -> -Ju
        # Thomas solve J dpsi = -Res
        n = N
        c = np.empty(n); dvec = np.empty(n)
        c[0] = Ju[0] / Jd0[0]; dvec[0] = -Res[0] / Jd0[0]
        for j in range(1, n):
            m = Jd0[j] - Jl[j] * c[j - 1]
            c[j] = Ju[j] / m if j < n - 1 else 0.0
            dvec[j] = (-Res[j] - Jl[j] * dvec[j - 1]) / m
        dpsi_step = np.empty(n)
        dpsi_step[-1] = dvec[-1]
        for j in range(n - 2, -1, -1):
            dpsi_step[j] = dvec[j] - c[j] * dpsi_step[j + 1]
        lam_damp = 1.0
        while lam_damp > 1e-4:
            trial = psi + lam_damp * dpsi_step
            if np.all(trial > 0.05):
                break
            lam_damp *= 0.5
        psi = psi + lam_damp * dpsi_step
    # derivatives of converged background
    psim = np.empty(N + 2)
    psim[1:-1] = psi; psim[0] = psi[0]; psim[-1] = 2.0 - psi[-1]
    d1 = (psim[2:] - psim[:-2]) / (2 * h)
    d2 = (psim[2:] - 2 * psim[1:-1] + psim[:-2]) / h**2
    return r, psi, d1, d2, rho, float(np.max(np.abs(Res)))

def Q_of(form_key, u, du, d2u, r, extra_rho=None):
    f = FORM[form_key]
    Q = f['W2'](u) - 0.5 * f['Fpp'](u) * du**2 - f['Fp'](u) * (d2u + 2 * du / r)
    if extra_rho is not None and form_key == 'F-A':
        # matter part U_rho' = -rho psi^4  =>  U_rho'' = -4 rho psi^3
        Q = Q - 4 * extra_rho * u**3
    return Q

def spectrum_on_background(form_key, r, u, du, d2u, l=0, extra_rho=None):
    f = FORM[form_key]
    F_nodes = f['F'](u)
    u_mid = 0.5 * (u[:-1] + u[1:])
    F_mid = f['F'](u_mid)
    Q = Q_of(form_key, u, du, d2u, r, extra_rho)
    return assemble_and_solve(r, F_nodes, F_mid, Q, l=l)

# ======================================================================
# B3 / B5 solitons
# ======================================================================
def soliton_profile(form_key, g0, R=60.0, N=4000):
    """Integrate radial soliton ODE, return grid arrays + diagnostics."""
    if form_key == 'F-S':
        def rhs(r_, y):
            g, gp = y
            g = max(g, G_GHOST + 1e-6)
            fk = 1 + 2 * ALPHA * np.log(g)
            drv = g**2 * (1 - g) - (ALPHA / g) * gp**2
            if r_ < 1e-10:
                return [gp, drv / (3 * fk)]
            return [gp, (drv - fk * 2 * gp / r_) / fk]
        def wall(r_, y):
            return y[0] - (G_GHOST + 0.005)
        wall.terminal = True
        wall.direction = -1
    else:                                              # F-S' substrate alpha=1
        def rhs(r_, y):
            g, gp = y
            g = max(g, 1e-3)
            drv = (1 - g) - (1 / g) * gp**2
            if r_ < 1e-10:
                return [gp, drv / 3]
            return [gp, drv - 2 * gp / r_]
        def wall(r_, y):
            return y[0] - (np.exp(-1) + 0.005)         # diagnostic only
        wall.terminal = False
        wall.direction = -1

    h = R / N
    r_grid = (np.arange(N) + 0.5) * h
    g_arr = np.empty(N); gp_arr = np.empty(N)
    r0, y = 1e-6, [g0, 0.0]
    bounces = 0
    idx0 = 0
    for b in range(30):
        sol = solve_ivp(rhs, [r0, R + h], y, method='DOP853', max_step=0.02,
                        rtol=1e-10, atol=1e-13, events=[wall], dense_output=True)
        seg_end = sol.t[-1]
        sel = (r_grid >= r0) & (r_grid <= seg_end)
        if np.any(sel):
            vals = sol.sol(r_grid[sel])
            g_arr[idx0:idx0 + sel.sum()] = vals[0]
            gp_arr[idx0:idx0 + sel.sum()] = vals[1]
            idx0 += sel.sum()
        if form_key == 'F-S' and sol.t_events[0].size > 0 and b < 29:
            r_b = float(sol.t_events[0][0])
            gp_b = float(sol.y_events[0][0, 1])
            r0 = r_b + 1e-6
            y = [G_GHOST + 0.005 + 1e-5, -gp_b]        # reflection (a3d convention)
            bounces += 1
        else:
            break
    n_fill = idx0
    if n_fill < N:                                     # trim if wall stopped early
        r_grid = r_grid[:n_fill]; g_arr = g_arr[:n_fill]; gp_arr = gp_arr[:n_fill]
    # g'' from ODE rhs (piecewise-consistent)
    d2 = np.array([rhs(r_grid[j], [g_arr[j], gp_arr[j]])[1] for j in range(len(r_grid))])
    return r_grid, g_arr, gp_arr, d2, bounces, float(np.min(g_arr))

# ======================================================================
# RUN
# ======================================================================
print("=" * 78)
print("Phase 2: BVP spectra.  Units gamma=beta=K_geo=1. TOL_NEG=%.0e" % TOL_NEG)
print("=" * 78)

print("\n[C2] B1 vacuum F-A (expect: all lam >= 0, edge = 1):")
for N in (2000, 4000, 8000):
    r, vals = vacuum_case('F-A', 1.0, N=N)
    nneg, lmin = report_spectrum("N=%d" % N, vals, expect="lam_min ~ 1")
edge_err = abs(vals[0] - 1.0)
print("  C2 verdict: %s (edge err %.2e, N_neg=%d)"
      % ("PASS" if (nneg == 0 and edge_err < 0.02) else "FAIL", edge_err, nneg))

print("\n[C3] B2 Yukawa F-A, linear profiles A=0.01, 0.1 (l=0):")
c3_ok = True
for A in (0.01, 0.1):
    for N in (2000, 4000, 8000):
        r, psi, dpsi, d2psi = yukawa_linear_background(A, 60.0, N)
        vals = spectrum_on_background('F-A', r, psi, dpsi, d2psi)
        nneg, lmin = report_spectrum("A=%.2f N=%d (min psi=%.3f, max psi=%.1f)"
                                     % (A, N, psi.min(), psi.max()), vals)
        c3_ok &= (nneg == 0)

print("\n[C3] B2 Yukawa F-A, NONLINEAR Newton backgrounds (targets ~0.3, ~1.0):")
for S in (0.9, 4.0):
    r, psi, dpsi, d2psi, rho, res = newton_background(S, 60.0, 4000)
    amp = psi.max() - 1
    vals = spectrum_on_background('F-A', r, psi, dpsi, d2psi, extra_rho=rho)
    nneg, lmin = report_spectrum("S=%.1f (amp=%.3f, res=%.1e)" % (S, amp, res), vals)
    c3_ok &= (nneg == 0)
print("  C3 verdict: %s" % ("PASS" if c3_ok else "FAIL"))

print("\n[C5] B4 vacuum F-S (dualism control; expect edge = -1):")
r, vals = vacuum_case('F-S', -1.0, N=4000)
nneg, lmin = report_spectrum("N=4000", vals, expect="lam_min ~ -1")
print("  C5 verdict: %s (measured edge %.4f)"
      % ("CONFIRMED-NEGATIVE (as predicted analytically)" if abs(lmin + 1) < 0.02
         else "UNEXPECTED", lmin))

print("\n[C4/C6a] B3 solitons F-S (crown log-form, alpha=2, wall reflections):")
sol_results = {}
for tag, g0 in (("e", G0_E), ("mu", G0_MU), ("tau", G0_TAU)):
    r, g, gp, d2, nb, gmin = soliton_profile('F-S', g0, R=60.0, N=4000)
    fmin = 1 + 4 * np.log(gmin)
    print("  g0^%-3s=%.5f: bounces=%d, min g=%.4f (g*=%.4f), min f(g)=%.4f, grid pts=%d"
          % (tag, g0, nb, gmin, G_GHOST, fmin, len(r)))
    for l in (0, 1):
        vals = spectrum_on_background('F-S', r, g, gp, d2, l=l)
        nneg, lmin = report_spectrum("    l=%d" % l, vals)
        sol_results[(tag, l)] = (nneg, lmin, nb)

print("\n[C4] B5 solitons F-S' (substrate alpha=1, dictionary-preferred):")
for tag, g0 in (("e", G0_E), ("mu", G0_MU), ("tau", G0_TAU)):
    r, g, gp, d2, nb, gmin = soliton_profile('F-Sp', g0, R=60.0, N=4000)
    print("  g0^%-3s=%.5f: min g=%.4f (ghost_log=%.4f n/a for poly K), grid pts=%d"
          % (tag, g0, gmin, np.exp(-1), len(r)))
    for l in (0, 1):
        vals = spectrum_on_background('F-Sp', r, g, gp, d2, l=l)
        nneg, lmin = report_spectrum("    l=%d" % l, vals)
        sol_results[('sub_' + tag, l)] = (nneg, lmin, nb)

print("\n[C4] grid convergence for soliton-e (F-S, l=0):")
for N in (2000, 4000, 8000):
    r, g, gp, d2, nb, gmin = soliton_profile('F-S', G0_E, R=60.0, N=N)
    vals = spectrum_on_background('F-S', r, g, gp, d2, l=0)
    report_spectrum("N=%d" % N, vals)
for R in (40.0, 80.0):
    r, g, gp, d2, nb, gmin = soliton_profile('F-S', G0_E, R=R, N=int(R / 60 * 4000))
    vals = spectrum_on_background('F-S', r, g, gp, d2, l=0)
    report_spectrum("R=%.0f" % R, vals)

print("\n" + "=" * 78)
print("SUMMARY TABLE (N_neg = modes with lambda < -1e-6, N=4000, R=60)")
for key in sorted(sol_results):
    nneg, lmin, nb = sol_results[key]
    print("  %-14s N_neg=%-3d lam_min=%+.5f  bounces=%d" % (str(key), nneg, lmin, nb))
print("=" * 78)
