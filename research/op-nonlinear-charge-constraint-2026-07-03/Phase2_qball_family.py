#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
op-nonlinear-charge-constraint (Phase 2) -- Q-ball family + VK criterion.

LOCK: Phase0_balance.md sec. 3 Phase 2 (P2a-P2d), criterion V2, sec. 5.
VK sign convention and box renormalization: LOCKED in Phase1 output
(P1c), BEFORE this file was run:
  Q_sol(omega) = omega * int r^2 [ f(phi) phi^2 - f(phi_inf) phi_inf^2 ] dr
  stable branch <=> dQ_sol/domega < 0  (Q > 0, omega > 0 normalization).

From Phase 1 (exact, sympy):
  profile ODE (M1 stationary):  W_eff(g) = W(g) - (omega^2/2) f(g) g^2,
    f g'' = W_eff'(g) - (1/2) f' g'^2 - (2/r) f g',
    W_eff'(g) = W'(g) - omega^2 g (f(g) + alpha)         [f = 1+2 alpha ln g]
  GSS split:  L+ : Q+ = W_eff''(phi) - f''/2 phi'^2 - f'(phi'' + 2 phi'/r)
              L- : Q- = W'/phi + f' phi'^2/(2 phi) - omega^2 (f + f' phi/2)
              L- phi_omega = 0 exactly (phase zero mode).
  omega-vacuum: W_eff'(phi_inf) = 0; phi_inf = 1 - 3 omega^2 + O(omega^4);
  edge(L+) = W_eff''(phi_inf)/f(phi_inf) = -1 - 7 omega^2 + O(omega^4);
  edge(L-) = 0 for every omega (exact);
  phi_inf crosses g* = e^{-1/4} at omega_gh = 0.293488 (f(phi_inf) < 0
  beyond: kinetically ghosted background -- edge undefined there; per
  LOCK V2 note this is reported as 'nierozstrzygalne z powodu tla').

Wall convention (LOCKED, Phase0 sec. 2): baseline hard wall with
reflection at g* + 0.005 (a3d/CP-7); cross-check ONLY f_eps at eps=0.2
(the single N-converged point of #62 W2b) + control eps=0.1.

Grids (LOCKED, sec. 5): omega in [0, 1], step 0.05 (full scan reported);
N in {2000, 4000, 8000}, R = 60 (R in {40, 80} control at V2-candidate
points); localized mode: lam < edge(omega) - 1e-3; DOP853 rtol 1e-10,
atol 1e-13, max_step 0.02; family/charge overlap threshold 0.9; mass
drift at omega -> 0: < 0.1% vs hard-wall baseline r21=206.73,
r31=3479.6 (#62).

Code reuse (verified in #60/#62): soliton profile integrator with
reflections, tridiagonal assembly, Haynsworth inertia + bisection for
constrained spectra, f_eps family, tail fit (generalized here to
oscillation wavenumber k(omega) = sqrt(-edge) around phi_inf(omega),
fit window r in [20, 35] as a3d).
"""
import numpy as np
from scipy.linalg import eigh_tridiagonal, solve_banded
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

PHI_GOLD = (1 + 5 ** 0.5) / 2
G0_E = 1.24915
G0_MU = PHI_GOLD * G0_E
G0_TAU = 3.18912
ALPHA = 2.0
G_GHOST = np.exp(-1 / (2 * ALPHA))          # g* = e^{-1/4}
GUARD = G_GHOST + 0.005
OMEGA_GRID = [round(0.05 * k, 2) for k in range(21)]   # 0.00 .. 1.00
R_BOX = 60.0
N_GRIDS = (2000, 4000, 8000)
LOC_MARGIN = 1e-3
R21_BASE = 206.73
R31_BASE = 3479.6

f_log = lambda g: 1 + 2 * ALPHA * np.log(g)
fp_log = lambda g: 2 * ALPHA / g
fpp_log = lambda g: -2 * ALPHA / g ** 2
W1_fun = lambda g: g ** 2 * (1 - g)         # W'
W2_fun = lambda g: 2 * g - 3 * g ** 2       # W''


def make_feps(eps):
    def s(g):
        fv = f_log(g)
        return np.sqrt(fv ** 2 + eps ** 2)
    F = lambda g: 0.5 * (f_log(g) + s(g))
    Fp = lambda g: 0.5 * fp_log(g) * (1 + f_log(g) / s(g))
    Fpp = lambda g: 0.5 * (fpp_log(g) * (1 + f_log(g) / s(g))
                           + fp_log(g) ** 2 * eps ** 2 / s(g) ** 3)
    return F, Fp, Fpp


# ---------------------------------------------------------------------
# effective potential derivatives (general kinetic F)
# ---------------------------------------------------------------------
def Weff1(g, om, F, Fp):
    """W_eff' = W' - (om^2/2)(F' g^2 + 2 F g)."""
    return W1_fun(g) - 0.5 * om ** 2 * (Fp(g) * g ** 2 + 2 * F(g) * g)


def Weff2(g, om, F, Fp, Fpp):
    """W_eff'' = W'' - (om^2/2)(F'' g^2 + 4 F' g + 2 F)."""
    return W2_fun(g) - 0.5 * om ** 2 * (Fpp(g) * g ** 2 + 4 * Fp(g) * g
                                        + 2 * F(g))


# ---------------------------------------------------------------------
# profiles
# ---------------------------------------------------------------------
def profile_hard(g0, om, R=R_BOX, N=4000):
    """Hard-wall Q-ball profile (reflection at g*+0.005, a3d convention).
    Returns r, g, gp, d2, bounces, gmin, filled(bool)."""
    def rhs(r_, y):
        g, gp = y
        g = max(g, G_GHOST + 1e-6)
        fk = f_log(g)
        drv = Weff1(g, om, f_log, fp_log) - (ALPHA / g) * gp ** 2
        if r_ < 1e-10:
            return [gp, drv / (3 * fk)]
        return [gp, (drv - fk * 2 * gp / r_) / fk]

    def wall(r_, y):
        return y[0] - GUARD
    wall.terminal = True
    wall.direction = -1

    h = R / N
    r_grid = (np.arange(N) + 0.5) * h
    g_arr = np.empty(N)
    gp_arr = np.empty(N)
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
            y = [GUARD + 1e-5, -gp_b]
            bounces += 1
        else:
            break
    filled = (idx0 == N)
    r_grid = r_grid[:idx0]
    g_arr = g_arr[:idx0]
    gp_arr = gp_arr[:idx0]
    d2 = np.array([rhs(r_grid[j], [g_arr[j], gp_arr[j]])[1]
                   for j in range(len(r_grid))])
    return r_grid, g_arr, gp_arr, d2, bounces, float(np.min(g_arr)) \
        if idx0 else np.nan, filled


def profile_soft(g0, om, eps, R=R_BOX, N=4000):
    """Smooth f_eps EL profile (no wall event), with omega^2 term."""
    F, Fp, Fpp = make_feps(eps)

    def rhs(r_, y):
        g, gp = y
        g = max(g, 1e-3)
        drv = Weff1(g, om, F, Fp) - 0.5 * Fp(g) * gp ** 2
        if r_ < 1e-10:
            return [gp, drv / (3 * F(g))]
        return [gp, drv / F(g) - 2 * gp / r_]

    h = R / N
    r_grid = (np.arange(N) + 0.5) * h
    sol = solve_ivp(rhs, [1e-6, R + h], [g0, 0.0], method='DOP853',
                    max_step=0.02, rtol=1e-10, atol=1e-13,
                    dense_output=True)
    r_end = float(sol.t[-1])
    sel = r_grid <= r_end
    r_g = r_grid[sel]
    if len(r_g) == 0:
        return r_g, np.array([]), np.array([]), np.array([]), r_end
    vals = sol.sol(r_g)
    g_arr, gp_arr = vals[0], vals[1]
    d2 = np.array([rhs(r_g[j], [g_arr[j], gp_arr[j]])[1]
                   for j in range(len(r_g))])
    return r_g, g_arr, gp_arr, d2, r_end


# ---------------------------------------------------------------------
# omega-vacuum and continuum edges
# ---------------------------------------------------------------------
def vacuum_of(om, F=f_log, Fp=fp_log, lo=0.4, hi=1.2):
    """phi_inf(omega): root of W_eff' = 0 on the branch through phi=1."""
    if om == 0.0:
        return 1.0
    fun = lambda x: Weff1(x, om, F, Fp)
    try:
        return float(brentq(fun, lo, hi, xtol=1e-14))
    except ValueError:
        return np.nan


def edge_plus(om, F=f_log, Fp=fp_log, Fpp=fpp_log):
    xv = vacuum_of(om, F, Fp)
    if not np.isfinite(xv):
        return np.nan, np.nan, np.nan
    return (Weff2(xv, om, F, Fp, Fpp) / F(xv), xv, F(xv))


# ---------------------------------------------------------------------
# tridiagonal operators L+ / L-
# ---------------------------------------------------------------------
def build_L(r, g, gp, d2, om, F, Fp, Fpp, which):
    N = len(r)
    h = r[1] - r[0]
    F_nodes = F(g)
    F_mid = F(0.5 * (g[:-1] + g[1:]))
    if which == 'plus':
        Qpot = (Weff2(g, om, F, Fp, Fpp) - 0.5 * Fpp(g) * gp ** 2
                - Fp(g) * (d2 + 2 * gp / r))
    else:
        Qpot = (W1_fun(g) / g + Fp(g) * gp ** 2 / (2 * g)
                - om ** 2 * (F(g) + 0.5 * Fp(g) * g))
    r_mid = 0.5 * (r[:-1] + r[1:])
    a = r_mid ** 2 * F_mid
    diag = np.empty(N)
    diag[0] = a[0] / h ** 2
    diag[1:-1] = (a[:-1] + a[1:]) / h ** 2
    diag[-1] = (a[-1] + r[-1] ** 2 * F_nodes[-1]) / h ** 2
    diag += r ** 2 * Qpot
    off = -a / h ** 2
    return diag / r ** 2, off / (r[:-1] * r[1:])


def low_eigs(d, e, k=120):
    k = min(k, len(d) - 1)
    return eigh_tridiagonal(d, e, select='i', select_range=(0, k - 1),
                            eigvals_only=True)


# ---------------------------------------------------------------------
# constrained spectrum via Haynsworth inertia (+ bisection) -- reuse #62
# ---------------------------------------------------------------------
def make_counter(d, e, Cmat, eigs_uncon):
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

    return Nc


def constrained_count_below(d, e, Cmat, eigs_uncon, edge):
    Nc = make_counter(d, e, Cmat, eigs_uncon)
    return Nc(edge)


def norm_col(c):
    n = np.linalg.norm(c)
    return c / n if n > 0 else c


# ---------------------------------------------------------------------
# tails (a3d generalized: oscillation around phi_inf with k(omega))
# ---------------------------------------------------------------------
def fit_tail_k(r_arr, g_arr, phi_inf, kk, r_L=20.0, r_R=35.0):
    mask = (r_arr >= r_L) & (r_arr <= r_R)
    if np.sum(mask) < 10 or not np.isfinite(kk) or kk <= 0:
        return float('nan')
    r_f = r_arr[mask]
    y_f = (g_arr[mask] - phi_inf) * r_f
    X = np.column_stack([np.cos(kk * r_f), np.sin(kk * r_f)])
    coef, _, _, _ = np.linalg.lstsq(X, y_f, rcond=None)
    return float(np.hypot(coef[0], coef[1]))


# ---------------------------------------------------------------------
# charges (box, R=60; raw + vacuum-subtracted per LOCKED convention)
# ---------------------------------------------------------------------
def charges(r, g, gp, om, F=f_log):
    h = r[1] - r[0]
    xv = vacuum_of(om)
    dens_q = F(g) * g ** 2
    dens_q_vac = F(xv) * xv ** 2 if np.isfinite(xv) else np.nan
    Q_raw = om * float(np.sum(r ** 2 * dens_q) * h)
    Q_sol = om * float(np.sum(r ** 2 * (dens_q - dens_q_vac)) * h)
    W_of = lambda x: x ** 3 / 3.0 - x ** 4 / 4.0
    e_dens = 0.5 * om ** 2 * F(g) * g ** 2 + 0.5 * F(g) * gp ** 2 + W_of(g)
    e_vac = 0.5 * om ** 2 * dens_q_vac + W_of(xv)
    E_raw = float(np.sum(r ** 2 * e_dens) * h)
    E_sol = float(np.sum(r ** 2 * (e_dens - e_vac)) * h)
    return Q_raw, Q_sol, E_raw, E_sol, xv


def family_direction(g0, om, N=4000, delta=0.01):
    """d phi_omega / d omega at fixed g0 (central difference), u = r*dphi."""
    outs = []
    for s in (+1, -1):
        omv = om + s * delta
        if omv < 0:
            return None
        r_, g_, _, _, _, _, filled = profile_hard(g0, omv, N=N)
        if not filled:
            return None
        outs.append(g_)
    if len(outs[0]) != len(outs[1]):
        return None
    return (outs[0] - outs[1]) / (2 * delta)


# =====================================================================
# RUN
# =====================================================================
print("=" * 78)
print("Phase 2: Q-ball family (M1, model-extension) + VK.  LOCK: Phase0 sec.3.")
print("Wall: hard (reflection g*+0.005) baseline; cross-check f_eps 0.2/0.1.")
print("omega grid: [0,1] step 0.05 (full).  N in {2000,4000,8000}, R=60.")
print("Localized: lam < edge_+(omega) - 1e-3 [L+], lam < -1e-3 [L-].")
print("=" * 78)

# ---------------- PART A: vacuum / edges ------------------------------
print("\n[A] omega-vacuum and continuum edges (hard-wall kinetic f):")
print("    %-6s %-11s %-11s %-12s %-9s %s"
      % ("omega", "phi_inf", "f(phi_inf)", "edge(L+)", "k_tail", "status"))
edges = {}
for om in OMEGA_GRID:
    ed, xv, fv = edge_plus(om)
    kk = np.sqrt(-ed) if (np.isfinite(ed) and ed < 0 and fv > 0) else np.nan
    ghosted = (fv <= 0) or (xv < GUARD)
    edges[om] = dict(edge=ed, xv=xv, fv=fv, k=kk, ghosted=ghosted)
    print("    %-6.2f %-11.6f %-11.6f %-12.6f %-9s %s"
          % (om, xv, fv, ed, ("%.4f" % kk) if np.isfinite(kk) else "n/a",
             "GHOSTED (f(phi_inf)<0 / vacuum below wall)" if ghosted else ""))
print("    edge(L-) = 0 dla kazdego omega (exact, Phase1).")

# ---------------- PART B: branches (P2a) ------------------------------
print("\n[B] P2a: branches phi_omega (hard wall, N=4000, R=60):")
BACKGROUNDS = (("e", G0_E), ("mu", G0_MU), ("tau", G0_TAU))
prof_cache = {}     # (tag, om, N) -> profile tuple
for tag, g0 in BACKGROUNDS:
    print("  [%s]  g0 = %.5f (FIXED, zero re-fitting)" % (tag, g0))
    print("    %-6s %-8s %-9s %-9s %-12s %s"
          % ("omega", "bounces", "min g", "filled", "A_tail(k)", "uwagi"))
    g_prev = None
    for om in OMEGA_GRID:
        r, g, gp, d2, nb, gmin, filled = profile_hard(g0, om, N=4000)
        prof_cache[(tag, om, 4000)] = (r, g, gp, d2, nb, gmin, filled)
        A = fit_tail_k(r, g, edges[om]['xv'], edges[om]['k']) \
            if filled and not edges[om]['ghosted'] else float('nan')
        note = ""
        if not filled:
            note = "INCOMPLETE (bounce cap / collapse)"
        if om == 0.05 and g_prev is not None and len(g) == len(g_prev):
            note += "  sup|phi_0.05-phi_0|=%.4f" % float(np.max(np.abs(
                g - g_prev)))
        print("    %-6.2f %-8d %-9.4f %-9s %-12s %s"
              % (om, nb, gmin, filled,
                 ("%.6f" % A) if A == A else "n/a", note))
        if om == 0.0:
            g_prev = g

# ---------------- PART C: Q(omega), E(omega), VK (P2b) ----------------
print("\n[C] P2b: Q(omega), E(omega) [pudlo R=60; raw + vacuum-subtracted")
print("    wg konwencji LOCKED w Phase1/P1c], znak dQ_sol/domega:")
qtab = {}
for tag, g0 in BACKGROUNDS:
    rows = []
    for om in OMEGA_GRID:
        r, g, gp, d2, nb, gmin, filled = prof_cache[(tag, om, 4000)]
        if not filled:
            rows.append((om, np.nan, np.nan, np.nan, np.nan))
            continue
        Q_raw, Q_sol, E_raw, E_sol, xv = charges(r, g, gp, om)
        rows.append((om, Q_raw, Q_sol, E_raw, E_sol))
    qtab[tag] = rows
    print("  [%s]" % tag)
    print("    %-6s %-13s %-13s %-13s %-13s %-11s %s"
          % ("omega", "Q_raw", "Q_sol", "E_raw", "E_sol", "dQsol/dom",
             "VK(dQ<0?)"))
    for i, (om, Qr, Qs, Er, Es) in enumerate(rows):
        if i == 0 or i == len(rows) - 1 or not np.isfinite(Qs):
            dq = np.nan
        else:
            Qm = rows[i - 1][2]
            Qp = rows[i + 1][2]
            dq = (Qp - Qm) / 0.10 if np.isfinite(Qm) and np.isfinite(Qp) \
                else np.nan
        vk = ("TAK" if dq < 0 else "nie") if np.isfinite(dq) else "-"
        print("    %-6.2f %-13.5g %-13.5g %-13.5g %-13.5g %-11.4g %s"
              % (om, Qr, Qs, Er, Es, dq, vk))
    # dE = omega dQ diagnostic (not a criterion)
    fin = [(om, Qs, Es) for (om, Qr, Qs, Er, Es) in rows
           if np.isfinite(Qs) and np.isfinite(Es)]
    if len(fin) > 4:
        difs = []
        for i in range(1, len(fin) - 1):
            om = fin[i][0]
            dQ = (fin[i + 1][1] - fin[i - 1][1]) / (fin[i + 1][0]
                                                    - fin[i - 1][0])
            dE = (fin[i + 1][2] - fin[i - 1][2]) / (fin[i + 1][0]
                                                    - fin[i - 1][0])
            if abs(dQ) > 1e-12:
                difs.append(dE / (om * dQ) if om > 0 else np.nan)
        difs = [d for d in difs if np.isfinite(d)]
        if difs:
            print("    diagnostyka dE/domega = omega dQ/domega: mediana "
                  "ratio = %.3f (odchylenia = czlony brzegowe pudla/odbic;"
                  " nie-kryterium)" % float(np.median(difs)))

# ---------------- PART D: spectra L+/L- + deflacja (P2c) --------------
print("\n[D] P2c: spektra L+/L- wokol phi_omega; N_loc z deflacja")
print("    kierunku ladunkowego (dQ/dphi ~ f'phi^2+2f phi) i rodziny")
print("    (dphi/domega); zbieznosc N in {2000,4000,8000}.")
print("    Dla omega z GHOSTED tlem: krawedz nieokreslona -> raport")
print("    surowego lambda_min, N_loc oznaczone 'n/o' (per LOCK V2 note).")

spec = {}
for tag, g0 in BACKGROUNDS:
    print("  [%s]" % tag)
    print("    %-6s %-4s %-28s %-10s %-8s %-8s %-12s %s"
          % ("omega", "N/k", "loc(L+) [lam<edge-1e-3]", "lam_min(L-)",
             "Nc(chg)", "Nc(c+f)", "phase-mode", "uwagi"))
    for om in OMEGA_GRID:
        ed = edges[om]['edge']
        ghosted = edges[om]['ghosted']
        loc_thresh = ed - LOC_MARGIN if not ghosted else np.nan
        fam_u = None
        per_grid = []
        for N in N_GRIDS:
            if (tag, om, N) in prof_cache:
                r, g, gp, d2, nb, gmin, filled = prof_cache[(tag, om, N)]
            else:
                r, g, gp, d2, nb, gmin, filled = profile_hard(g0, om, N=N)
                prof_cache[(tag, om, N)] = (r, g, gp, d2, nb, gmin, filled)
            if not filled:
                per_grid.append((N, None))
                continue
            dP, eP = build_L(r, g, gp, d2, om, f_log, fp_log, fpp_log,
                             'plus')
            dM, eM = build_L(r, g, gp, d2, om, f_log, fp_log, fpp_log,
                             'minus')
            eigsP = low_eigs(dP, eP)
            eigsM = low_eigs(dM, eM, k=8)
            # phase zero mode check: overlap of L- ground pair with r*phi
            valsM, vecsM = eigh_tridiagonal(dM, eM, select='i',
                                            select_range=(0, 2))
            uphi = norm_col(r * g)
            ovl = max(abs(float(vecsM[:, j] @ uphi)) for j in range(3))
            if not ghosted:
                locP = [v for v in eigsP if v < loc_thresh]
                # constrained counts (inertia)
                w_chg = fp_log(g) * g ** 2 + 2 * f_log(g) * g
                C1 = norm_col(w_chg * r)[:, None]
                nc = constrained_count_below(dP, eP, C1, eigsP, loc_thresh)
                ncf = None
                if N == 4000:
                    if fam_u is None:
                        fd = family_direction(g0, om, N=N)
                        fam_u = norm_col(fd * r) if fd is not None else False
                    if fam_u is not False and fam_u is not None:
                        C2 = np.column_stack([C1[:, 0], fam_u])
                        ncf = constrained_count_below(dP, eP, C2, eigsP,
                                                      loc_thresh)
                per_grid.append((N, dict(nloc=len(locP), lams=locP,
                                         lmM=float(eigsM[0]), nc=nc,
                                         ncf=ncf, ovl=ovl)))
            else:
                per_grid.append((N, dict(nloc=None,
                                         lams=[float(eigsP[0])],
                                         lmM=float(eigsM[0]), nc=None,
                                         ncf=None, ovl=ovl)))
        # print per-grid rows
        for N, row in per_grid:
            if row is None:
                print("    %-6.2f %-4d PROFIL NIEKOMPLETNY -- pomijam" %
                      (om, N))
                continue
            if row['nloc'] is None:
                lam_s = "raw lam_min(L+)=%.4f" % row['lams'][0]
                nl = "n/o"
            else:
                lam_s = "%d: %s" % (row['nloc'],
                                    ["%.4f" % v for v in row['lams'][:4]])
                nl = "%d" % row['nloc']
            print("    %-6.2f %-4d %-28s %-10.4f %-8s %-8s %-12.4f %s"
                  % (om, N, lam_s, row['lmM'],
                     str(row['nc']) if row['nc'] is not None else "n/o",
                     str(row['ncf']) if row['ncf'] is not None else "-",
                     row['ovl'],
                     "GHOSTED" if row['nloc'] is None else ""))
        spec[(tag, om)] = per_grid

# ---------------- PART E: mass control (P2d) --------------------------
print("\n[E] P2d: kontrola mas (A_tail^4, okno [20,35], k(omega)):")
print("    baseline #62 (hard wall, omega=0): r21=%.2f r31=%.1f" %
      (R21_BASE, R31_BASE))
print("    %-6s %-12s %-10s %-12s %-10s" % ("omega", "r21", "drift%",
                                            "r31", "drift%"))
mass_rows = []
for om in OMEGA_GRID:
    if edges[om]['ghosted']:
        print("    %-6.2f tlo GHOSTED -- fit ogona niezdefiniowany" % om)
        continue
    As = {}
    ok = True
    for tag, g0 in BACKGROUNDS:
        r, g, gp, d2, nb, gmin, filled = prof_cache[(tag, om, 4000)]
        if not filled:
            ok = False
            break
        As[tag] = fit_tail_k(r, g, edges[om]['xv'], edges[om]['k'])
        if not np.isfinite(As[tag]) or As[tag] <= 0:
            ok = False
            break
    if not ok:
        print("    %-6.2f fit niedostepny (profil/ogon)" % om)
        continue
    r21 = (As['mu'] / As['e']) ** 4
    r31 = (As['tau'] / As['e']) ** 4
    d21 = 100 * abs(r21 - R21_BASE) / R21_BASE
    d31 = 100 * abs(r31 - R31_BASE) / R31_BASE
    mass_rows.append((om, r21, d21, r31, d31))
    print("    %-6.2f %-12.3f %-10.3f %-12.1f %-10.3f"
          % (om, r21, d21, r31, d31))
if mass_rows and mass_rows[0][0] == 0.0:
    gate = (mass_rows[0][2] < 0.1) and (mass_rows[0][4] < 0.1)
    print("    gate omega->0 (<0.1%%): drift r21=%.4f%%, r31=%.4f%% -> %s"
          % (mass_rows[0][2], mass_rows[0][4],
             "PASS" if gate else "FAIL"))

# ---------------- PART F: cross-check soft wall (f_eps) ---------------
print("\n[F] cross-check f_eps (LOCK: WYLACZNIE eps=0.2 + kontrola 0.1):")
for eps in (0.2, 0.1):
    F, Fp, Fpp = make_feps(eps)
    print("  eps = %.1f:" % eps)
    print("    %-4s %-6s %-9s %-8s %-12s %s"
          % ("tag", "omega", "min g", "r_end", "lam_min(L+)", "uwagi"))
    for tag, g0 in BACKGROUNDS:
        for om in (0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30):
            r, g, gp, d2, r_end = profile_soft(g0, om, eps, N=4000)
            if len(r) < 200 or r_end < R_BOX:
                print("    %-4s %-6.2f KOLAPS/urwanie r_end=%.2f" %
                      (tag, om, r_end))
                continue
            dP, eP = build_L(r, g, gp, d2, om, F, Fp, Fpp, 'plus')
            eigsP = low_eigs(dP, eP, k=8)
            # soft-wall vacuum and edge
            xv = vacuum_of(om, F, Fp, lo=0.05, hi=1.2)
            edv = (Weff2(xv, om, F, Fp, Fpp) / F(xv)) \
                if np.isfinite(xv) else np.nan
            print("    %-4s %-6.2f %-9.4f %-8.1f %-12.4f edge=%.4f"
                  % (tag, om, float(np.min(g)), r_end, float(eigsP[0]),
                     edv))

# ---------------- PART G: V2 verdict ----------------------------------
print("\n" + "=" * 78)
print("WERDYKT V2 (koniunkcja (i)-(iii), LOCKED; + klauzula tla):")
print("  KLAUZULA TLA (LOCK): kryterium ma sens TYLKO jesli P1d/P2c")
print("  wykazuja usuniecie/podniesienie kontinuum.  Phase1: edge(L+) =")
print("  -1 - 7 omega^2 + O(omega^4) < -1 dla kazdego omega; dla omega >")
print("  omega_gh=0.2935 tlo zduchowione.  Kontinuum pozostaje tachioniczne")
print("  przy KAZDYM omega.")
for tag, g0 in BACKGROUNDS:
    # (i) existence of omega_k > 0 on a branch continuous from omega=0
    cont = []
    for om in OMEGA_GRID:
        _, _, _, _, _, _, filled = prof_cache[(tag, om, 4000)]
        if not filled:
            break
        cont.append(om)
    # (ii) VK at omega_k, (iii) constrained N_loc = 0 on 3 grids
    cands = []
    rows = qtab[tag]
    for i, om in enumerate(OMEGA_GRID):
        if om == 0.0 or om not in cont:
            continue
        if edges[om]['ghosted']:
            continue
        if i == 0 or i >= len(rows) - 1:
            continue
        Qm, Qp = rows[i - 1][2], rows[i + 1][2]
        if not (np.isfinite(Qm) and np.isfinite(Qp)):
            continue
        vk_ok = (Qp - Qm) < 0
        pg = spec[(tag, om)]
        ncs = [row['nc'] for (N, row) in pg if row is not None
               and row.get('nc') is not None]
        iii_ok = (len(ncs) == 3 and all(n == 0 for n in ncs))
        if vk_ok and iii_ok:
            cands.append(om)
    print("  [%-3s] galaz ciagla (kompletne profile): omega <= %s;"
          % (tag, ("%.2f" % cont[-1]) if cont else "BRAK"))
    print("        kandydaci omega_k (VK i (iii) jednoczesnie): %s"
          % (cands if cands else "BRAK"))
print("\n  (pelne brzmienie werdyktu -> README; jesli kandydatow brak i/lub")
print("   kontinuum tachioniczne dla kazdego omega: wynik per LOCK --")
print("   negatywny/nierozstrzygalny zgloszony wprost)")
print("=" * 78)
