#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
op-wall-dynamics (Phase 2) -- wall as one-sided condition + regularizations.

LOCK: Phase0_balance.md sec. 2 Phase 2 (W2a, W2b), criteria W2.

W2a: spectrum with one-sided condition g_eq + v >= g* on the contact set
     (where g_eq < g* + 0.01), approximated by Dirichlet v = 0 on the
     contact set (exact submatrix decoupling).  Full linear-eigenvalue
     complementarity problem (LCP): NOT well-posed in the strict
     linearization limit -- the obstacle v >= g* - g_eq sits at finite
     distance (g_eq >= g* + 0.005 by construction of the reflected
     profile), so for infinitesimal fluctuation amplitude the cone
     constraint is inactive and the problem reduces to the
     unconstrained one.  This LIMITATION is recorded per LOCK
     ("inaczej odnotowac ograniczenie").  As a finite-amplitude
     supplement we minimize the quadratic form over the obstacle cone
     at fixed amplitude (projected gradient) and report the result.

W2b: soft-wall family f_eps(g) = (1/2)[f(g) + sqrt(f(g)^2 + eps^2)]
     (smooth, > 0), eps in {0.2, 0.1, 0.05, 0.02} (LOCKED):
     re-solve mu/tau profiles (ODE with f_eps) and spectra.
     Locked questions: (i) does lam_min(eps) converge as eps -> 0 or
     diverge; (ii) are r21, r31 (A_tail^4 / Koide chain, a3d) robust
     in eps (allowed drift < 0.1%; otherwise report sensitivity).
     (electron profile is re-solved with the same f_eps because
     r21 = (A_mu/A_e)^4, r31 = (A_tau/A_e)^4 require A_e(eps).)

Units gamma=beta=K_geo=1.  Localized mode: lam < -1 - 1e-3 (LOCK 4).
Baseline hard-wall values from CP-7 / a3d reflection convention.
"""
import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.integrate import solve_ivp

PHI = (1 + 5 ** 0.5) / 2
G0_E = 1.24915
G0_MU = PHI * G0_E
G0_TAU = 3.18912
ALPHA = 2.0
G_GHOST = np.exp(-1 / (2 * ALPHA))
LOC_EDGE = -1.0 - 1e-3
EPS_LIST = (0.2, 0.1, 0.05, 0.02)          # LOCKED
CONTACT_MARGIN = 0.01                      # LOCKED (g_eq < g* + 0.01)

f_log   = lambda g: 1 + 2 * ALPHA * np.log(g)
fp_log  = lambda g: 2 * ALPHA / g
fpp_log = lambda g: -2 * ALPHA / g ** 2
W1_fun  = lambda g: g ** 2 * (1 - g)       # W'
W2_fun  = lambda g: 2 * g - 3 * g ** 2     # W''


def make_feps(eps):
    def s(g):
        fv = f_log(g)
        return np.sqrt(fv ** 2 + eps ** 2)
    F   = lambda g: 0.5 * (f_log(g) + s(g))
    Fp  = lambda g: 0.5 * fp_log(g) * (1 + f_log(g) / s(g))
    Fpp = lambda g: 0.5 * (fpp_log(g) * (1 + f_log(g) / s(g))
                           + fp_log(g) ** 2 * eps ** 2 / s(g) ** 3)
    return F, Fp, Fpp


# ---------------------------------------------------------------------
# profiles
# ---------------------------------------------------------------------
def soliton_profile_hard(g0, R=60.0, N=4000):
    """CP-7 hard-wall profile (reflection g'->-g' at g*+0.005), verbatim."""
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


def soliton_profile_soft(g0, eps, R=60.0, N=4000):
    """EL profile with smooth positive f_eps (no wall event, no reflection).
    EL: g'' = [W'(g) - (1/2) f_eps'(g) g'^2] / f_eps(g) - (2/r) g'."""
    F, Fp, Fpp = make_feps(eps)

    def rhs(r_, y):
        g, gp = y
        g = max(g, 1e-3)
        drv = W1_fun(g) - 0.5 * Fp(g) * gp ** 2
        if r_ < 1e-10:
            return [gp, drv / (3 * F(g))]
        return [gp, drv / F(g) - 2 * gp / r_]

    h = R / N
    r_grid = (np.arange(N) + 0.5) * h
    sol = solve_ivp(rhs, [1e-6, R + h], [g0, 0.0], method='DOP853',
                    max_step=0.02, rtol=1e-10, atol=1e-13,
                    dense_output=True)
    r_end = sol.t[-1]
    sel = r_grid <= r_end
    r_g = r_grid[sel]
    vals = sol.sol(r_g)
    g_arr, gp_arr = vals[0], vals[1]
    d2 = np.array([rhs(r_g[j], [g_arr[j], gp_arr[j]])[1]
                   for j in range(len(r_g))])
    # sub-wall crossings (diagnostic analogue of "bounces")
    below = g_arr < G_GHOST
    crossings = int(np.sum(np.diff(below.astype(int)) == 1))
    return r_g, g_arr, gp_arr, d2, crossings, float(np.min(g_arr)), r_end


# ---------------------------------------------------------------------
# spectra
# ---------------------------------------------------------------------
def build_tridiag_gen(r, g, gp, d2, F, Fp, Fpp, l=0):
    N = len(r)
    h = r[1] - r[0]
    F_nodes = F(g)
    F_mid = F(0.5 * (g[:-1] + g[1:]))
    Q = W2_fun(g) - 0.5 * Fpp(g) * gp ** 2 - Fp(g) * (d2 + 2 * gp / r)
    r_mid = 0.5 * (r[:-1] + r[1:])
    a = r_mid ** 2 * F_mid
    diag = np.empty(N)
    diag[0] = a[0] / h ** 2
    diag[1:-1] = (a[:-1] + a[1:]) / h ** 2
    diag[-1] = (a[-1] + r[-1] ** 2 * F_nodes[-1]) / h ** 2
    diag += r ** 2 * (Q + F_nodes * l * (l + 1) / r ** 2)
    off = -a / h ** 2
    return diag / r ** 2, off / (r[:-1] * r[1:])


def localized_modes(d, e, edge=LOC_EDGE, k=60):
    vals = eigh_tridiagonal(d, e, select='i', select_range=(0, k - 1),
                            eigvals_only=True)
    return [float(v) for v in vals if v < edge], float(vals[0])


# ---------------------------------------------------------------------
# tails and mass ratios (a3d convention: g-1 ~ (B cos r + C sin r)/r,
# fit window r in [20, 35], m ~ A_tail^4)
# ---------------------------------------------------------------------
def fit_tail(r_arr, g_arr, r_L=20.0, r_R=35.0):
    mask = (r_arr >= r_L) & (r_arr <= r_R)
    if np.sum(mask) < 10:
        return float('nan')
    r_f = r_arr[mask]
    y_f = (g_arr[mask] - 1.0) * r_f
    X = np.column_stack([np.cos(r_f), np.sin(r_f)])
    coef, _, _, _ = np.linalg.lstsq(X, y_f, rcond=None)
    return float(np.hypot(coef[0], coef[1]))


# ---------------------------------------------------------------------
# W2a
# ---------------------------------------------------------------------
print("=" * 78)
print("Phase 2: wall models.  W2a one-sided (contact Dirichlet), W2b soft wall.")
print("=" * 78)
print("\n[W2a] Dirichlet v=0 on contact set {r: g_eq < g* + %.2f}:" %
      CONTACT_MARGIN)
print("  (exact: rows/cols of contact nodes removed; spectrum = union over")
print("   decoupled segments; localized modes lam < %.4f)" % LOC_EDGE)

for tag, g0 in (("mu", G0_MU), ("tau", G0_TAU)):
    for N in (2000, 4000, 8000):
        r, g, gp, d2, nb, gmin = soliton_profile_hard(g0, R=60.0, N=N)
        d, e = build_tridiag_gen(r, g, gp, d2, f_log, fp_log, fpp_log)
        contact = g < (G_GHOST + CONTACT_MARGIN)
        ncontact = int(contact.sum())
        # unconstrained reference
        loc0, lmin0 = localized_modes(d, e)
        # segments of non-contact indices
        idx = np.where(~contact)[0]
        segs = np.split(idx, np.where(np.diff(idx) != 1)[0] + 1)
        loc_dir = []
        for s in segs:
            if len(s) < 3:
                continue
            ds = d[s].copy()
            es = e[s[:-1]] if s[-1] < len(d) - 1 else e[s[:-1]]
            es = es[:len(s) - 1]
            lv, _ = localized_modes(ds, es, k=min(60, len(s) - 1))
            loc_dir.extend(lv)
        loc_dir.sort()
        print("  %-4s N=%d: contact pts=%d (segments=%d)  free: N_loc=%d %s"
              % (tag, N, ncontact, len(segs), len(loc0),
                 ["%.4f" % v for v in loc0]))
        print("           Dirichlet-on-contact: N_loc=%d %s"
              % (len(loc_dir), ["%.4f" % v for v in loc_dir]))

print("\n[W2a/LCP] limitation (recorded per LOCK): the obstacle")
print("  v >= g* - g_eq has finite distance (g_eq >= g*+0.005), so in the")
print("  strict linearization limit the complementarity cone is inactive")
print("  and the LCP reduces to the unconstrained spectrum.  Finite-")
print("  amplitude supplement (projected gradient on the obstacle cone):")


def cone_min(tag, g0, amp, N=4000, iters=4000, lr=None):
    """min <v,Lv>/<v,v> over ||v||_B = amp, v >= g* - g_eq (nodewise),
    projected gradient in symmetrized coords u = r v, obstacle
    u >= r*(g* - g_eq)  (negative near contact, very negative elsewhere)."""
    r, g, gp, d2, nb, gmin = soliton_profile_hard(g0, R=60.0, N=N)
    d, e = build_tridiag_gen(r, g, gp, d2, f_log, fp_log, fpp_log)
    h = r[1] - r[0]
    lower = r * (G_GHOST - g) / np.sqrt(h)  # u-coords with ||.||_2 = L2(B)
    # matvec
    def Lu(u):
        out = d * u
        out[:-1] += e * u[1:]
        out[1:] += e * u[:-1]
        return out
    # start from unconstrained ground state
    vals, vecs = eigh_tridiagonal(d, e, select='i', select_range=(0, 0))
    u = vecs[:, 0] * amp
    if np.mean(u) < 0:
        u = -u
    u = np.maximum(u, lower)
    u *= amp / np.linalg.norm(u)
    if lr is None:
        lr = 0.02 / max(1.0, abs(vals[0]))
    q_prev = None
    for it in range(iters):
        gvec = Lu(u)
        q = float(u @ gvec) / float(u @ u)
        gr = gvec - q * u
        u = u - lr * amp ** 2 * gr / max(np.linalg.norm(gr), 1e-30) * 0.05
        u = np.maximum(u, lower)
        u *= amp / np.linalg.norm(u)
        if q_prev is not None and abs(q - q_prev) < 1e-12:
            break
        q_prev = q
    nact = int(np.sum(np.abs(u - lower) < 1e-14))
    return q, vals[0], nact


for tag, g0 in (("mu", G0_MU), ("tau", G0_TAU)):
    for amp in (0.01, 0.05):
        q, l0, nact = cone_min(tag, g0, amp)
        print("  %-4s amp=%.2f: min Rayleigh on cone = %+.6f"
              " (unconstrained lam_min=%+.6f, active nodes=%d)"
              % (tag, amp, q, l0, nact))

# ---------------------------------------------------------------------
# W2b
# ---------------------------------------------------------------------
print("\n" + "=" * 78)
print("[W2b] soft wall f_eps = (f + sqrt(f^2+eps^2))/2, eps in %s"
      % (EPS_LIST,))

# baseline (hard wall, a3d reflection) for drift comparison
base = {}
print("\n  baseline (hard wall, reflection, N=4000, R=60):")
for tag, g0 in (("e", G0_E), ("mu", G0_MU), ("tau", G0_TAU)):
    r, g, gp, d2, nb, gmin = soliton_profile_hard(g0)
    A = fit_tail(r, g)
    d, e = build_tridiag_gen(r, g, gp, d2, f_log, fp_log, fpp_log)
    loc, lmin = localized_modes(d, e)
    base[tag] = dict(A=A, loc=loc, lmin=lmin)
    print("    %-4s bounces=%d  min g=%.4f  A_tail=%.6f  N_loc=%d  lam=%s"
          % (tag, nb, gmin, A, len(loc), ["%.4f" % v for v in loc]))
r21_base = (base["mu"]["A"] / base["e"]["A"]) ** 4
r31_base = (base["tau"]["A"] / base["e"]["A"]) ** 4
print("    r21_base=%.3f  r31_base=%.3f  (PDG 206.77 / 3477.4)"
      % (r21_base, r31_base))

soft = {}
for eps in EPS_LIST:
    F, Fp, Fpp = make_feps(eps)
    print("\n  eps = %.2f:" % eps)
    for tag, g0 in (("e", G0_E), ("mu", G0_MU), ("tau", G0_TAU)):
        r, g, gp, d2, ncr, gmin, r_end = soliton_profile_soft(g0, eps)
        A = fit_tail(r, g)
        row = dict(A=A, gmin=gmin, ncr=ncr, r_end=r_end)
        if len(r) > 100:
            d, e = build_tridiag_gen(r, g, gp, d2, F, Fp, Fpp)
            loc, lmin = localized_modes(d, e)
            row["loc"] = loc; row["lmin"] = lmin
            # grid convergence for mu/tau (LOCK: no mode without 3 grids)
            convs = []
            if tag in ("mu", "tau"):
                for N2 in (2000, 8000):
                    r2, g2, gp2, d22, _, _, _ = soliton_profile_soft(
                        g0, eps, N=N2)
                    d2t, e2t = build_tridiag_gen(r2, g2, gp2, d22, F, Fp, Fpp)
                    loc2, lm2 = localized_modes(d2t, e2t)
                    convs.append((N2, len(loc2), lm2))
            row["convs"] = convs
        soft[(eps, tag)] = row
        msg = ("    %-4s min g=%.4f  sub-wall crossings=%d  r_end=%.1f"
               "  A_tail=%s" % (tag, gmin, ncr, r_end,
                                ("%.6f" % A) if A == A else "n/a"))
        if "loc" in row:
            msg += ("  N_loc=%d  lam_min=%.4f" % (len(row["loc"]),
                                                  row["lmin"]))
        print(msg)
        for (N2, nl2, lm2) in row.get("convs", []):
            print("         conv N=%d: N_loc=%d  lam_min=%.4f" % (N2, nl2, lm2))
        if "loc" in row and tag in ("mu", "tau"):
            print("         localized: %s" % ["%.4f" % v for v in row["loc"]])

print("\n  TABLE lam_min(eps) (l=0, N=4000):")
print("    %-6s %-12s %-12s" % ("eps", "mu", "tau"))
print("    %-6s %-12s %-12s" % ("hard", "%.4f" % base["mu"]["lmin"],
                                "%.4f" % base["tau"]["lmin"]))
for eps in EPS_LIST:
    print("    %-6.2f %-12s %-12s"
          % (eps,
             "%.4f" % soft[(eps, "mu")].get("lmin", float('nan')),
             "%.4f" % soft[(eps, "tau")].get("lmin", float('nan'))))

print("\n  TABLE r21/r31 drift vs hard-wall baseline (criterion: <0.1%):")
print("    %-6s %-12s %-10s %-12s %-10s" % ("eps", "r21", "drift%",
                                            "r31", "drift%"))
for eps in EPS_LIST:
    Ae = soft[(eps, "e")]["A"]; Am = soft[(eps, "mu")]["A"]
    At = soft[(eps, "tau")]["A"]
    if Ae == Ae and Am == Am and At == At and Ae > 0:
        r21 = (Am / Ae) ** 4; r31 = (At / Ae) ** 4
        print("    %-6.2f %-12.3f %-10.3f %-12.1f %-10.3f"
              % (eps, r21, 100 * abs(r21 - r21_base) / r21_base,
                 r31, 100 * abs(r31 - r31_base) / r31_base))
    else:
        print("    %-6.2f tail fit unavailable (profile ended at r=%.1f)"
              % (eps, soft[(eps, "tau")]["r_end"]))

print("\n" + "=" * 78)
print("W2 verdict inputs complete (interpretation -> README per LOCK).")
print("=" * 78)
