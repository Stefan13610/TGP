#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
bs2d.py -- Relativistic baby-Skyrme model in 2+1D (arena for e-p scattering test).

ONTOLOGY / NON-CIRCULARITY RULES (enforced by construction):
  * Field n(x,y,t) in S^2 (unit 3-vector). NOTHING ELSE.
  * NO gauge field, NO pairwise charge law q_i q_j, NO Coulomb propagator inserted.
  * Interaction between solitons arises ONLY from the O(3)+Skyrme field energy overlap.
  * "Charge" is measured OPERATIONALLY as momentum transfer of a soliton pushed
    through another (RAMY sec.4). Sign/selectivity must FALL OUT, not be inserted.

Energy (mostly-minus metric):
  E = INT [ 1/2 |grad n|^2 + kappa^2/2 |d1 n x d2 n|^2 + mu^2 (1 - n3) ] d^2x

Topological charge (pi_2(S^2)=Z):
  Q = 1/(4pi) INT n . (d1 n x d2 n) d^2x

Hamiltonian dynamics on the sphere (relativistic, 2nd order in time, lossless):
  dt n   = pi
  dt pi  = -P_n[ dU/dn ] - (pi.pi) n
  dU/dn  = -lap n - mu^2 e3 + kappa^2 [ dy(dx n x w) - dx(dy n x w) ],  w = dx n x dy n

Energy conservation is MONITORED (not assumed) -> the "lossless" claim is a measurement.
"""
import numpy as np

# ---------- field calculus on a periodic grid (vectorized numpy) ----------

def d_dx(a, dx):
    return (np.roll(a, -1, axis=2) - np.roll(a, 1, axis=2)) / (2.0 * dx)

def d_dy(a, dx):
    return (np.roll(a, -1, axis=1) - np.roll(a, 1, axis=1)) / (2.0 * dx)

def lap(a, dx):
    return (np.roll(a, -1, axis=2) + np.roll(a, 1, axis=2)
            + np.roll(a, -1, axis=1) + np.roll(a, 1, axis=1)
            - 4.0 * a) / (dx * dx)

def cross(a, b):
    """cross product of two (3,Ny,Nx) vector fields along axis 0."""
    return np.stack((
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ), axis=0)

def normalize(n):
    nn = np.sqrt(np.sum(n * n, axis=0))
    return n / nn

# ---------- energy, charge, variational gradient ----------

def energy_density(n, dx, kappa2, mu2):
    dxn = d_dx(n, dx)
    dyn = d_dy(n, dx)
    w = cross(dxn, dyn)
    e_grad = 0.5 * (np.sum(dxn * dxn, axis=0) + np.sum(dyn * dyn, axis=0))
    e_sky = 0.5 * kappa2 * np.sum(w * w, axis=0)
    e_pot = mu2 * (1.0 - n[2])
    return e_grad + e_sky + e_pot, e_grad, e_sky, e_pot

def total_energy(n, dx, kappa2, mu2):
    e, eg, es, ep = energy_density(n, dx, kappa2, mu2)
    cell = dx * dx
    return (e.sum() * cell, eg.sum() * cell, es.sum() * cell, ep.sum() * cell)

def charge_density(n, dx):
    dxn = d_dx(n, dx)
    dyn = d_dy(n, dx)
    w = cross(dxn, dyn)
    q = np.sum(n * w, axis=0) / (4.0 * np.pi)
    return q

def topological_charge(n, dx):
    return charge_density(n, dx).sum() * dx * dx

def dU_dn(n, dx, kappa2, mu2):
    """Variational derivative of the STATIC energy w.r.t. n (unconstrained)."""
    dxn = d_dx(n, dx)
    dyn = d_dy(n, dx)
    w = cross(dxn, dyn)
    grad = -lap(n, dx)
    grad[2] -= mu2  # -mu^2 e3
    sky = kappa2 * (d_dy(cross(dxn, w), dx) - d_dx(cross(dyn, w), dx))
    return grad + sky

def tangential_force(n, dx, kappa2, mu2):
    """-(dU/dn) projected onto tangent space of S^2."""
    g = dU_dn(n, dx, kappa2, mu2)
    gn = np.sum(g * n, axis=0)
    g_tan = g - gn * n          # projection
    return -g_tan                # force = -grad

# ---------- initial configurations ----------

def make_grid(N, L):
    dx = L / N
    xs = (np.arange(N) - N / 2) * dx
    X, Y = np.meshgrid(xs, xs, indexing='xy')  # X varies along axis 2, Y along axis 1
    return X, Y, dx

def hedgehog(X, Y, x0=0.0, y0=0.0, w=3.0, charge=+1, phase=0.0):
    """Baby-skyrmion hedgehog. charge=+1 skyrmion, -1 antiskyrmion.
    f: pi at center -> 0 at infinity (n3: -1 -> +1 vacuum)."""
    xr = X - x0
    yr = Y - y0
    r = np.sqrt(xr * xr + yr * yr) + 1e-12
    theta = np.arctan2(yr, xr)
    f = np.pi * np.exp(-r / w)               # profile, relaxed later
    m = charge
    n = np.stack((
        np.sin(f) * np.cos(m * theta + phase),
        np.sin(f) * np.sin(m * theta + phase),
        np.cos(f),
    ), axis=0)
    return normalize(n)

# ---------- 1D hedgehog profile (robust, avoids 2D parabolic stiffness) ----------

def profile_energy_1d(f, dr, kappa2, mu2):
    """Reduced hedgehog energy E = 2pi INT r dr [ 1/2 f'^2 + 1/2 sin^2 f/r^2
       + 1/2 kappa^2 f'^2 sin^2 f/r^2 + mu^2 (1-cos f) ]. f on nodes r_i=i*dr,
       f[0]=pi, f[-1]=0 pinned. Returns total E."""
    M = f.size - 1
    r = np.arange(M + 1) * dr
    fp = np.zeros_like(f)
    fp[1:M] = (f[2:M + 1] - f[0:M - 1]) / (2 * dr)
    s = np.sin(f); c = np.cos(f)
    i = np.arange(1, M)
    ri = r[i]
    P = ri * dr * (0.5 * fp[i] ** 2 + 0.5 * s[i] ** 2 / ri ** 2
                   + 0.5 * kappa2 * fp[i] ** 2 * s[i] ** 2 / ri ** 2
                   + mu2 * (1 - c[i]))
    return 2 * np.pi * P.sum()

def profile_grad_1d(f, dr, kappa2, mu2):
    """Analytic dE/df_j for interior nodes (boundaries pinned)."""
    M = f.size - 1
    r = np.arange(M + 1) * dr
    fp = np.zeros_like(f)
    fp[1:M] = (f[2:M + 1] - f[0:M - 1]) / (2 * dr)
    s = np.sin(f); c = np.cos(f)
    g = np.zeros_like(f)
    D = np.zeros_like(f)                      # D_i = dP_i/dfp_i on interior nodes
    ii = np.arange(1, M)
    D[ii] = r[ii] * dr * fp[ii] * (1.0 + kappa2 * s[ii] ** 2 / r[ii] ** 2)
    direct = np.zeros_like(f)
    direct[ii] = r[ii] * dr * (s[ii] * c[ii] / r[ii] ** 2
                               + kappa2 * fp[ii] ** 2 * s[ii] * c[ii] / r[ii] ** 2
                               + mu2 * s[ii])
    # chain via fp: node j gets +D[j-1]/(2dr) - D[j+1]/(2dr)
    g[ii] = 2 * np.pi * (direct[ii] + (D[ii - 1] - D[ii + 1]) / (2 * dr))
    g[0] = 0.0; g[M] = 0.0
    return g

def relax_profile_1d(M, Rmax, kappa2, mu2, w0=None):
    """Relax the 1D hedgehog profile with L-BFGS-B (exact analytic gradient).
    Boundaries pinned: f(0)=pi, f(Rmax)=0. Returns r, f(r)."""
    from scipy.optimize import minimize, minimize_scalar
    dr = Rmax / M
    r = np.arange(M + 1) * dr

    def atan_fam(w):
        f = 4.0 * np.arctan(np.exp(-r / w)); f[0] = np.pi; f[-1] = 0.0
        return f
    # pre-optimize the scale on a stable ansatz (avoids discrete collapse)
    sc = minimize_scalar(lambda w: profile_energy_1d(atan_fam(w), dr, kappa2, mu2),
                         bounds=(0.4, 0.4 * Rmax), method='bounded')
    f0 = atan_fam(sc.x)

    def build(x):
        f = np.empty(M + 1)
        f[0] = np.pi; f[-1] = 0.0; f[1:M] = x
        return f

    def fun(x):
        return profile_energy_1d(build(x), dr, kappa2, mu2)

    def jac(x):
        return profile_grad_1d(build(x), dr, kappa2, mu2)[1:M]

    # moderate tolerance: stop at the physical soliton, not the discrete-collapse state
    res = minimize(fun, f0[1:M], jac=jac, method='L-BFGS-B',
                   options=dict(maxiter=5000, ftol=1e-12, gtol=1e-9))
    return r, build(res.x)

def build_from_profile(X, Y, r_prof, f_prof, x0=0.0, y0=0.0, charge=+1, phase=0.0):
    """Embed a relaxed radial profile f(r) as a 2D hedgehog field."""
    xr = X - x0; yr = Y - y0
    r = np.sqrt(xr * xr + yr * yr)
    theta = np.arctan2(yr, xr)
    f = np.interp(r, r_prof, f_prof, right=0.0)
    m = charge
    n = np.stack((np.sin(f) * np.cos(m * theta + phase),
                  np.sin(f) * np.sin(m * theta + phase),
                  np.cos(f)), axis=0)
    return normalize(n)

def bump_Q0(X, Y, x0=0.0, y0=0.0, g0=2.0, w=2.0, chi=0.0):
    """Non-topological Q=0 lump: n tilts from n3=cos(g0) at center to vacuum n3=1,
    with NO azimuthal winding (constant internal phase chi). Candidate 'neutral balloon'
    seed; evolves into an oscillon if the model supports it."""
    r = np.sqrt((X - x0) ** 2 + (Y - y0) ** 2)
    g = g0 * np.exp(-(r / w) ** 2)
    n = np.stack((np.sin(g) * np.cos(chi),
                  np.sin(g) * np.sin(chi),
                  np.cos(g)), axis=0)
    return normalize(n)

# ---------- static relaxation (gradient descent on the sphere) ----------

def relax(n, dx, kappa2, mu2, steps=4000, dt=0.05, mom=0.9, tol=1e-6, verbose=False):
    """Heavy-ball gradient descent on S^2 to find a static soliton."""
    v = np.zeros_like(n)
    for it in range(steps):
        F = tangential_force(n, dx, kappa2, mu2)   # descent direction (toward lower E)
        v = mom * v + dt * F
        n = normalize(n + v)
        if it % 200 == 0 or it == steps - 1:
            fmax = np.sqrt(np.sum(F * F, axis=0)).max()
            if verbose:
                Etot = total_energy(n, dx, kappa2, mu2)[0]
                print(f"  relax it={it:5d}  |F|max={fmax:.3e}  E={Etot:.5f}")
            if fmax < tol:
                break
    return n

# ---------- relativistic dynamics (lossless, 2nd order in time) ----------

def accel(n, pi, dx, kappa2, mu2):
    """dt n = pi ;  dt pi = -P_n[dU/dn] - (pi.pi) n  (geodesic + force on S^2)."""
    g = dU_dn(n, dx, kappa2, mu2)
    gn = np.sum(g * n, axis=0)
    g_tan = g - gn * n
    pi2 = np.sum(pi * pi, axis=0)
    return pi, -g_tan - pi2 * n

def project(n, pi):
    """Renormalize n onto S^2 and remove the pi component along n (keep pi tangent)."""
    n = normalize(n)
    pn = np.sum(pi * n, axis=0)
    pi = pi - pn * n
    return n, pi

def rk4_step(n, pi, dt, dx, kappa2, mu2):
    k1n, k1p = accel(n, pi, dx, kappa2, mu2)
    k2n, k2p = accel(n + 0.5 * dt * k1n, pi + 0.5 * dt * k1p, dx, kappa2, mu2)
    k3n, k3p = accel(n + 0.5 * dt * k2n, pi + 0.5 * dt * k2p, dx, kappa2, mu2)
    k4n, k4p = accel(n + dt * k3n, pi + dt * k3p, dx, kappa2, mu2)
    n = n + (dt / 6.0) * (k1n + 2 * k2n + 2 * k3n + k4n)
    pi = pi + (dt / 6.0) * (k1p + 2 * k2p + 2 * k3p + k4p)
    return project(n, pi)

def kinetic_energy(pi, dx):
    return 0.5 * np.sum(pi * pi) * dx * dx

def total_energy_dyn(n, pi, dx, kappa2, mu2):
    Estat = total_energy(n, dx, kappa2, mu2)[0]
    return Estat + kinetic_energy(pi, dx)

def field_momentum(n, pi, dx):
    """Sigma-model canonical momentum P_i = INT pi . d_i n d^2x (dominant part)."""
    Px = np.sum(pi * d_dx(n, dx)) * dx * dx
    Py = np.sum(pi * d_dy(n, dx)) * dx * dx
    return Px, Py

def com_of_charge(n, X, Y, dx, region=None):
    """Center of mass weighted by |topological charge density|. region=(xmin,xmax) optional."""
    q = np.abs(charge_density(n, dx))
    if region is not None:
        mask = (X >= region[0]) & (X < region[1])
        q = q * mask
    W = q.sum()
    if W < 1e-12:
        return np.nan, np.nan, 0.0
    xc = np.sum(q * X) / W
    yc = np.sum(q * Y) / W
    return xc, yc, W * dx * dx  # last = |Q| in region

def boost_soliton(X, Y, r_prof, f_prof, x0, y0, v, charge=+1, phase=0.0):
    """Build a soliton moving with velocity v in +x (Lorentz-contracted), returns n, pi."""
    gamma = 1.0 / np.sqrt(1.0 - v * v)
    # contract x about x0
    Xc = x0 + (X - x0) * gamma
    n = build_from_profile(Xc, Y, r_prof, f_prof, x0=x0, y0=y0, charge=charge, phase=phase)
    # dt n = -v d_x n  (moving rigidly to the right)
    pi = -v * d_dx(n, dx_of(X))
    n, pi = project(n, pi)
    return n, pi

def dx_of(X):
    return X[0, 1] - X[0, 0]

# ---------- radial profile extraction (for diagnostics) ----------

def radial_profile_f(n, X, Y, x0=0.0, y0=0.0, nbins=60):
    r = np.sqrt((X - x0) ** 2 + (Y - y0) ** 2).ravel()
    f = np.arccos(np.clip(n[2].ravel(), -1.0, 1.0))  # f = arccos(n3)
    rmax = r.max() / np.sqrt(2)
    bins = np.linspace(0, rmax, nbins + 1)
    idx = np.digitize(r, bins) - 1
    rc = 0.5 * (bins[:-1] + bins[1:])
    fbar = np.array([f[idx == k].mean() if np.any(idx == k) else np.nan
                     for k in range(nbins)])
    return rc, fbar
