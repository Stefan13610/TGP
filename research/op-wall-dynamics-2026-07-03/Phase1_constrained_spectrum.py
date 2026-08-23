#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
op-wall-dynamics (Phase 1) -- constrained stability of crown solitons.

LOCK: Phase0_balance.md (2026-07-03), section 2 Phase 1 + section 4.
Question: do the localized negative modes of mu/tau (CP-7: mu=2, tau=3,
l=0) disappear on the subspace of budget-preserving fluctuations?

Constraints (ALL computed, none selected post-hoc; LOCKED):
  K1 (simple budget):    int v * r^2 dr = 0            (w = 1)
  K2 (metric budget):    int v * g_eq^2 * r^2 dr = 0   (w = g_eq^2)
  K3 (kinetic budget):   int v * f(g_eq) * r^2 dr = 0  (w = f(g_eq))
  K4 (core budget):      int_{r<r_core} v * r^2 dr = 0 (w = 1_{r<r_core}),
                         r_core = first node of g_eq - 1
  Pairs K_i ^ K_j: PRE-declared in LOCK (one scalar constraint removes
  at most 1 negative mode by interlacing; mu has 2, tau has 3).

Method (LOCK: projection P L P on the constraint-orthogonal subspace,
"deflacja w dyskretyzacji tridiagonalnej ... uwaga na wage B=r^2"):
  Discretization identical to CP-7 Phase2_bvp_spectrum.py (staggered
  grid, self-adjoint flux form, symmetrization B^(-1/2) M B^(-1/2),
  B = diag(r^2)).  Physical constraint  sum_j v_j w(r_j) r_j^2 h = 0
  becomes, in symmetrized coordinates u_j = r_j v_j:
      c_j = w(r_j) * r_j ,   c^T u = 0.
  Exact spectrum of P C P restricted to span{c_i}^perp via inertia
  (Haynsworth):  N_constr(< t) = n_-(C - tI) + n_-(S(t)) - p,
  S(t) = -Chat^T (C - tI)^(-1) Chat  (banded solve; p = #constraints).
  Constrained eigenvalues located by bisection on the counting
  function; eigenvectors from  u = (C - lam)^(-1) Chat y,  y in ker S.
  Cross-validation: dense eigendecomposition of P C P at N=2000.

Definitions (LOCK section 4): gamma=beta=K_geo=1; TOL_NEG=-1e-6;
localized mode: lam < -1 - 1e-3 (below continuum edge -1);
grids N in {2000,4000,8000}, R=60 (control R in {40,80} borderline).

NOTE (methodological, fixed before run): the tachyonic continuum from
-1 is a property of the F-S vacuum (CP-7 C5, box-count floor(R/pi));
scalar constraints cannot remove it.  W1 counts are therefore the
LOCALIZED modes per LOCK section 4; raw N_neg(< -1e-6) is reported
alongside for completeness.

Reused verbatim from CP-7 (verified by convergence there):
soliton_profile, Q_of, tridiagonal assembly; constants G0_E=1.24915,
G0_MU=phi*G0_E, G0_TAU=3.18912.
"""
import numpy as np
from scipy.linalg import eigh_tridiagonal, eigh, solve_banded, qr
from scipy.integrate import solve_ivp

PHI = (1 + 5 ** 0.5) / 2
G0_E = 1.24915
G0_MU = PHI * G0_E
G0_TAU = 3.18912
ALPHA = 2.0
G_GHOST = np.exp(-1 / (2 * ALPHA))          # 0.7788
TOL_NEG = -1e-6
LOC_EDGE = -1.0 - 1e-3                      # localized-mode threshold (LOCK 4)

# ---------------------------------------------------------------------
# F-S formulation (CP-7)
# ---------------------------------------------------------------------
F_fun   = lambda u: 1 + 4 * np.log(u)
Fp_fun  = lambda u: 4 / u
Fpp_fun = lambda u: -4 / u ** 2
W2_fun  = lambda u: 2 * u - 3 * u ** 2


def soliton_profile(g0, R=60.0, N=4000):
    """CP-7 soliton_profile, F-S branch only (wall reflections as a3d)."""
    def rhs(r_, y):
        g, gp = y
        g = max(g, G_GHOST + 1e-6)
        fk = 1 + 2 * ALPHA * np.log(g)
        drv = g ** 2 * (1 - g) - (ALPHA / g) * gp ** 2
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
    bounces = 0
    idx0 = 0
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


def build_tridiag(r, g, gp, d2, l=0):
    """CP-7 assembly -> symmetric tridiagonal (d, e) acting on u = r*v."""
    N = len(r)
    h = r[1] - r[0]
    F_nodes = F_fun(g)
    g_mid = 0.5 * (g[:-1] + g[1:])
    F_mid = F_fun(g_mid)
    Q = W2_fun(g) - 0.5 * Fpp_fun(g) * gp ** 2 - Fp_fun(g) * (d2 + 2 * gp / r)
    r_mid = 0.5 * (r[:-1] + r[1:])
    a = r_mid ** 2 * F_mid
    diag = np.empty(N)
    diag[0] = a[0] / h ** 2
    diag[1:-1] = (a[:-1] + a[1:]) / h ** 2
    diag[-1] = (a[-1] + r[-1] ** 2 * F_nodes[-1]) / h ** 2
    diag += r ** 2 * (Q + F_nodes * l * (l + 1) / r ** 2)
    off = -a / h ** 2
    d = diag / r ** 2
    e = off / (r[:-1] * r[1:])
    return d, e


# ---------------------------------------------------------------------
# constraints
# ---------------------------------------------------------------------
def first_node_radius(r, g):
    """r_core = first zero of g_eq - 1 (LOCK K4)."""
    s = g - 1.0
    for j in range(len(r) - 1):
        if s[j] == 0.0 or s[j] * s[j + 1] < 0:
            return r[j] + (r[j + 1] - r[j]) * abs(s[j]) / (abs(s[j]) + abs(s[j + 1]))
    return None


def constraint_matrix(keys, r, g, r_core):
    """Columns c_j = w(r_j)*r_j (symmetrized coords), each normalized."""
    cols = []
    for k in keys:
        if k == 'K1':
            w = np.ones_like(r)
        elif k == 'K2':
            w = g ** 2
        elif k == 'K3':
            w = F_fun(g)
        elif k == 'K4':
            w = (r < r_core).astype(float)
        else:
            raise ValueError(k)
        c = w * r
        cols.append(c / np.linalg.norm(c))
    return np.column_stack(cols)


# ---------------------------------------------------------------------
# constrained spectrum via inertia + bisection
# ---------------------------------------------------------------------
def make_counter(d, e, Cmat, eigs_uncon):
    """Return Nc(t) = #{constrained eigenvalues < t} (exact, Haynsworth)."""
    N = len(d)
    p = Cmat.shape[1]
    ab_template = np.zeros((3, N))
    ab_template[0, 1:] = e
    ab_template[2, :-1] = e

    def S_of(t):
        ab = ab_template.copy()
        ab[1, :] = d - t
        X = solve_banded((1, 1), ab, Cmat)
        S = -(Cmat.T @ X)
        return 0.5 * (S + S.T), X

    def Nc(t):
        m = int(np.searchsorted(eigs_uncon, t))
        for shift in (0.0, 1e-9, -1e-9, 1e-8):
            try:
                S, _ = S_of(t + shift)
                sv = np.linalg.eigvalsh(S)
                return m + int(np.sum(sv < 0)) - p
            except Exception:
                continue
        raise RuntimeError('singular solve at t=%r' % t)

    return Nc, S_of


def constrained_modes_below(d, e, Cmat, eigs_uncon, edge, tol=1e-9):
    """All constrained eigenvalues < edge, by bisection on Nc."""
    Nc, S_of = make_counter(d, e, Cmat, eigs_uncon)
    lo0 = float(eigs_uncon[0]) - 0.5
    total = Nc(edge)
    lams = []
    for k in range(1, total + 1):
        lo, hi = lo0, edge
        while hi - lo > tol:
            mid = 0.5 * (lo + hi)
            if Nc(mid) >= k:
                hi = mid
            else:
                lo = mid
        lams.append(0.5 * (lo + hi))
    return total, lams, S_of


def constrained_eigvec(S_of, lam):
    """u = (C-lam)^(-1) Chat y, y = null vector of S(lam); normalized."""
    S, X = S_of(lam + 1e-9)
    w_, V = np.linalg.eigh(S)
    y = V[:, int(np.argmin(np.abs(w_)))]
    u = X @ y
    return u / np.linalg.norm(u)


# ---------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------
def localized(vals, edge=LOC_EDGE):
    return [v for v in vals if v < edge]


def profile_family_direction(g0, r_ref, delta=1e-4, R=60.0, N=4000):
    """dg/dg0 by central difference, in symmetrized coords, normalized."""
    out = []
    for s in (+1, -1):
        r_, g_, _, _, _, _ = soliton_profile(g0 + s * delta, R=R, N=N)
        if len(g_) != len(r_ref):
            return None
        out.append(g_)
    dg = (out[0] - out[1]) / (2 * delta)
    u = dg * r_ref
    n = np.linalg.norm(u)
    return u / n if n > 0 else None


SINGLES = ['K1', 'K2', 'K3', 'K4']
PAIRS = [('K1', 'K2'), ('K1', 'K3'), ('K1', 'K4'),
         ('K2', 'K3'), ('K2', 'K4'), ('K3', 'K4')]
ALL_SETS = [(s,) for s in SINGLES] + PAIRS

# ---------------------------------------------------------------------
# RUN
# ---------------------------------------------------------------------
print("=" * 78)
print("Phase 1: constrained spectra (budget constraints K1-K4 + pairs).")
print("Units gamma=beta=K_geo=1.  TOL_NEG=%.0e  localized: lam < %.4f" %
      (TOL_NEG, LOC_EDGE))
print("=" * 78)

results = {}          # (bg, N, set) -> (nloc, lams, raw_note)
backgrounds = (("e", G0_E), ("mu", G0_MU), ("tau", G0_TAU))

for tag, g0 in backgrounds:
    print("\n" + "-" * 78)
    print("[background %s]  g0 = %.5f" % (tag, g0))
    for N in (2000, 4000, 8000):
        r, g, gp, d2, nb, gmin = soliton_profile(g0, R=60.0, N=N)
        d, e = build_tridiag(r, g, gp, d2, l=0)
        eigs = eigh_tridiagonal(d, e, select='i', select_range=(0, 59),
                                eigvals_only=True)
        loc0 = localized(eigs)
        r_core = first_node_radius(r, g)
        if N == 4000:
            print("  profile: bounces=%d, min g=%.4f, r_core(first node)=%s"
                  % (nb, gmin, ("%.4f" % r_core) if r_core else "None"))
            print("  UNCONSTRAINED check vs CP-7: N_loc=%d  %s"
                  % (len(loc0), ["%.4f" % v for v in loc0]))
        for keys in ALL_SETS:
            if 'K4' in keys and r_core is None:
                results[(tag, N, keys)] = (None, [], 'no r_core')
                continue
            Cmat = constraint_matrix(keys, r, g, r_core)
            if Cmat.shape[1] == 2:
                cosang = float(abs(Cmat[:, 0] @ Cmat[:, 1]))
            else:
                cosang = None
            nloc, lams, S_of = constrained_modes_below(d, e, Cmat, eigs,
                                                       LOC_EDGE)
            nraw = None
            if N == 4000:   # raw count vs -1e-6 needs eigs above edge too
                Nc, _ = make_counter(d, e, Cmat, eigs)
                if eigs[-1] > TOL_NEG:
                    nraw = Nc(TOL_NEG)
            results[(tag, N, keys)] = (nloc, lams, None)
            label = "^".join(keys)
            extra = ""
            if cosang is not None and N == 4000:
                extra = "  cos(c_i,c_j)=%.6f" % cosang
            if nraw is not None:
                extra += "  N_neg_raw=%d" % nraw
            print("    N=%d  %-8s N_loc=%d  %s%s"
                  % (N, label, nloc,
                     ["%.6f" % v for v in lams], extra))

# ---------------------------------------------------------------------
# R-control for borderline cases (LOCK 4): mu/tau, pairs, R in {40,80}
# ---------------------------------------------------------------------
print("\n" + "-" * 78)
print("[R-control] mu/tau, all constraint sets, R in {40, 80} (N ~ R/60*4000)")
for tag, g0 in (("mu", G0_MU), ("tau", G0_TAU)):
    for R in (40.0, 80.0):
        N = int(R / 60.0 * 4000)
        r, g, gp, d2, nb, gmin = soliton_profile(g0, R=R, N=N)
        d, e = build_tridiag(r, g, gp, d2, l=0)
        eigs = eigh_tridiagonal(d, e, select='i', select_range=(0, 59),
                                eigvals_only=True)
        r_core = first_node_radius(r, g)
        for keys in ALL_SETS:
            Cmat = constraint_matrix(keys, r, g, r_core)
            nloc, lams, _ = constrained_modes_below(d, e, Cmat, eigs, LOC_EDGE)
            print("  %-4s R=%.0f  %-8s N_loc=%d  %s"
                  % (tag, R, "^".join(keys), nloc,
                     ["%.6f" % v for v in lams]))

# ---------------------------------------------------------------------
# dense cross-validation at N=2000 (P C P + penalty on span{c})
# ---------------------------------------------------------------------
print("\n" + "-" * 78)
print("[validation] dense P C P eigendecomposition, N=2000:")
for tag, g0, keys in (("mu", G0_MU, ('K1',)), ("mu", G0_MU, ('K1', 'K3')),
                      ("tau", G0_TAU, ('K2',)), ("tau", G0_TAU, ('K2', 'K3'))):
    N = 2000
    r, g, gp, d2, nb, gmin = soliton_profile(g0, R=60.0, N=N)
    d, e = build_tridiag(r, g, gp, d2, l=0)
    r_core = first_node_radius(r, g)
    Cmat = constraint_matrix(keys, r, g, r_core)
    Qc, _ = qr(Cmat, mode='economic')
    M = np.diag(d)
    M[np.arange(N - 1), np.arange(1, N)] = e
    M[np.arange(1, N), np.arange(N - 1)] = e
    P = np.eye(N) - Qc @ Qc.T
    A = P @ M @ P + 1e6 * (Qc @ Qc.T)
    dense_vals = eigh(A, eigvals_only=True, subset_by_index=(0, 9))
    dense_loc = localized(dense_vals)
    eigs = eigh_tridiagonal(d, e, select='i', select_range=(0, 59),
                            eigvals_only=True)
    nloc, lams, _ = constrained_modes_below(d, e, Cmat, eigs, LOC_EDGE)
    dmax = max([abs(a - b) for a, b in zip(sorted(dense_loc), sorted(lams))],
               default=0.0)
    ok = (len(dense_loc) == nloc) and dmax < 1e-6
    print("  %-4s %-8s dense N_loc=%d %s | inertia N_loc=%d %s  max|diff|=%.1e  %s"
          % (tag, "^".join(keys), len(dense_loc),
             ["%.6f" % v for v in dense_loc], nloc,
             ["%.6f" % v for v in lams], dmax, "PASS" if ok else "FAIL"))

# ---------------------------------------------------------------------
# residual-mode identification (LOCK W1: tau residual allowed ONLY if
# overlap with dg/dg0 > 0.9)
# ---------------------------------------------------------------------
print("\n" + "-" * 78)
print("[residual test] overlap of surviving constrained modes with dg/dg0")
print("(central diff, delta in {1e-4, 1e-3}; N=4000, R=60)")
for tag, g0 in (("mu", G0_MU), ("tau", G0_TAU)):
    N = 4000
    r, g, gp, d2, nb, gmin = soliton_profile(g0, R=60.0, N=N)
    d, e = build_tridiag(r, g, gp, d2, l=0)
    eigs = eigh_tridiagonal(d, e, select='i', select_range=(0, 59),
                            eigvals_only=True)
    r_core = first_node_radius(r, g)
    fams = {}
    for delta in (1e-4, 1e-3):
        fams[delta] = profile_family_direction(g0, r, delta=delta, N=N)
    for keys in ALL_SETS:
        nloc, lams, S_of = constrained_modes_below(
            d, e, constraint_matrix(keys, r, g, r_core), eigs, LOC_EDGE)
        if nloc == 0 or nloc > 1:
            continue
        u = constrained_eigvec(S_of, lams[0])
        for delta, fam in fams.items():
            if fam is None:
                print("  %-4s %-8s lam=%.6f  delta=%.0e: dg/dg0 UNAVAILABLE"
                      " (profile length mismatch)" % (tag, "^".join(keys),
                                                      lams[0], delta))
                continue
            ov = float(abs(u @ fam))
            print("  %-4s %-8s lam=%.6f  delta=%.0e: |<v,dg/dg0>| = %.4f  %s"
                  % (tag, "^".join(keys), lams[0], delta, ov,
                     "(>0.9: family direction)" if ov > 0.9 else ""))

# ---------------------------------------------------------------------
# W1 verdict per LOCK
# ---------------------------------------------------------------------
print("\n" + "=" * 78)
print("W1 VERDICT (criteria as LOCKED in Phase0_balance.md sec. 2):")


def conv_ok(tag, keys):
    ns = [results.get((tag, N, keys), (None,))[0] for N in (2000, 4000, 8000)]
    return (None not in ns) and len(set(ns)) == 1, ns


any_single_pass = False
for keys in [('K1',), ('K2',), ('K3',)]:
    okm, nm = conv_ok("mu", keys)
    okt, nt = conv_ok("tau", keys)
    nmu = nm[1]; ntau = nt[1]
    verdict = (nmu == 0 and ntau is not None and ntau <= 1)
    print("  single %-3s : N_loc(mu)=%s %s, N_loc(tau)=%s %s -> %s"
          % (keys[0], nmu, "conv" if okm else "NONCONV",
             ntau, "conv" if okt else "NONCONV",
             "candidate-PASS (tau residual test required)" if verdict
             else "not sufficient"))
    if verdict:
        any_single_pass = True

print("\n  pre-declared extensions (pairs + K4), same criterion:")
for keys in ALL_SETS:
    if keys in [('K1',), ('K2',), ('K3',)]:
        continue
    okm, nm = conv_ok("mu", keys)
    okt, nt = conv_ok("tau", keys)
    nmu = nm[1]; ntau = nt[1]
    verdict = (nmu == 0 and ntau is not None and ntau <= 1)
    print("  %-8s : N_loc(mu)=%s %s, N_loc(tau)=%s %s -> %s"
          % ("^".join(keys), nmu, "conv" if okm else "NONCONV",
             ntau, "conv" if okt else "NONCONV",
             "candidate-PASS (tau residual test required)" if verdict
             else "not sufficient"))

print("\n  (final W1 wording -> README; negative result reported verbatim")
print("   per LOCK if modes survive all K1-K3 and declared combinations)")
print("=" * 78)
