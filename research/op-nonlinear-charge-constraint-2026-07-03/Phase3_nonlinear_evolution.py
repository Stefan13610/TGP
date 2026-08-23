#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
op-nonlinear-charge-constraint (Phase 3) -- nonlinear dynamical test (M0).

LOCK: Phase0_balance.md sec. 3 Phase 3 (P3a-P3b), criterion V3, sec. 5.
Independent of the M1 branch: this is the CANONICAL embedding M0 with the
smooth wall f_eps (W2b, #62), eps = 0.2 (the ONLY N-converged soft-wall
spectral point: mu, lam_min = -1.389) + control eps = 0.1.

tau: NO EL representative in the smooth f_eps family (collapse at every
eps -- #62 W2b); tau is OUT OF SCOPE of Phase 3 per LOCK P3a (recorded).

P3a (LOCKED):
  g(r,0) = g_eq + a * v_deep,  g_t(r,0) = 0,  a in {+-0.01, +-0.03};
  RK4 (allowed by LOCK: 'leapfrog lub RK4'), dt convergence (2 steps);
  energy gate |Delta E|/E < 1e-6 for the whole run; R = 60, CP-7 grid
  (N = 4000, staggered r_j = (j+1/2) h).

P3b (LOCKED):
  a(t) = <v_deep, delta g>_B  (B = r^2 dr, v_deep normalized in B);
  (i) growth rate sigma vs sqrt(|lam_min|) = sqrt(1.389) = 1.1786
      (linear-theory validation);
  (ii) saturation level / recurrences up to t = 200.

V3 (LOCKED): (i) exponential growth consistent with lam_min (+-20%) and
no saturation below ||delta g|| ~ 10% of background -> 'niestabilnosc
potwierdzona nieliniowo w M0-f_eps'; (ii) saturation < 10% + recurrences
-> 'efektywna stabilizacja nieliniowa'.  Either direction = execution PASS.

METHOD NOTE (fixed BEFORE runs, no post-hoc changes):
 1. Semi-discrete Hamiltonian dynamics.  We evolve the EXACT Hamiltonian
    system of the discretized energy
      H = sum_j h r_j^2 [ pi_j^2 / (2 F(g_j)) + W(g_j) ]
        + 1/2 sum_m h r_mid^2 F(g_mid) ((g_{m+1}-g_m)/h)^2,
    g_mid = (g_m + g_{m+1})/2, zero-flux ends (no ghost cells).
    This is the M0 EOM (F = f_eps) in flux form; dH/dt = 0 EXACTLY for
    the ODE system, so the energy gate measures pure integrator error.
 2. Linear-theory rates.  The CP-7/#62 spectral convention solves
    L v = lam v with weight 1 (B = r^2); the dynamical linearization is
    F(g_eq) v_tt = -L v, i.e. growth sigma^2 = -lam_gen of the
    F-weighted problem L v = lam_gen F v.  BOTH are computed and
    reported; the LOCKED comparison (V3 i) is against sqrt(1.389)
    (weight-1, as written in the LOCK), the F-weighted value is
    reported as the exact linear prediction (documentation, not a
    criterion change).
 3. v_deep = weight-1 eigenvector of lam_min (as in #62), normalized
    ||v||_B = 1.  'Background amplitude' for the 10% clause :=
    max_r |g_eq(r) - 1| (the crown amplitude over the F-S vacuum).
 4. dt convergence (LOCK '2 kroki'): the a = +0.01 run is repeated with
    dt/2; remaining runs use the base dt.  Base dt = 0.004, control
    dt = 0.002 (CFL: dt << h = 0.015).
"""
import sys
import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.integrate import solve_ivp

PHI_GOLD = (1 + 5 ** 0.5) / 2
G0_E = 1.24915
G0_MU = PHI_GOLD * G0_E
G0_TAU = 3.18912
ALPHA = 2.0
R_BOX = 60.0
N_GRID = 4000
T_END = 200.0
DT_BASE = 0.004
AMPS = (0.01, -0.01, 0.03, -0.03)
GATE = 1e-6
LAM_LOCK = 1.389                      # #62 W2b: mu, eps=0.2, N-converged

f_log = lambda g: 1 + 2 * ALPHA * np.log(g)
fp_log = lambda g: 2 * ALPHA / g
fpp_log = lambda g: -2 * ALPHA / g ** 2
W_fun = lambda g: g ** 3 / 3.0 - g ** 4 / 4.0
W1_fun = lambda g: g ** 2 * (1 - g)
W2_fun = lambda g: 2 * g - 3 * g ** 2


def make_feps(eps):
    def s(g):
        fv = f_log(g)
        return np.sqrt(fv ** 2 + eps ** 2)
    F = lambda g: 0.5 * (f_log(g) + s(g))
    Fp = lambda g: 0.5 * fp_log(g) * (1 + f_log(g) / s(g))
    Fpp = lambda g: 0.5 * (fpp_log(g) * (1 + f_log(g) / s(g))
                           + fp_log(g) ** 2 * eps ** 2 / s(g) ** 3)
    return F, Fp, Fpp


def profile_soft(g0, eps, R=R_BOX, N=N_GRID):
    """#62 Phase2 soliton_profile_soft (omega = 0), verbatim convention."""
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
    r_end = float(sol.t[-1])
    sel = r_grid <= r_end
    r_g = r_grid[sel]
    vals = sol.sol(r_g)
    g_arr, gp_arr = vals[0], vals[1]
    d2 = np.array([rhs(r_g[j], [g_arr[j], gp_arr[j]])[1]
                   for j in range(len(r_g))])
    return r_g, g_arr, gp_arr, d2, r_end


def build_tridiag_gen(r, g, gp, d2, F, Fp, Fpp):
    """#62 assembly (weight-1 problem on u = r v), verbatim."""
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
    diag += r ** 2 * Q
    off = -a / h ** 2
    return diag / r ** 2, off / (r[:-1] * r[1:])


# ---------------------------------------------------------------------
# semi-discrete Hamiltonian dynamics (exact discrete energy law)
# ---------------------------------------------------------------------
class SoftWallDynamics:
    def __init__(self, r, eps):
        self.r = r
        self.h = r[1] - r[0]
        self.r2 = r ** 2
        self.rm2 = (0.5 * (r[:-1] + r[1:])) ** 2
        self.F, self.Fp, self.Fpp = make_feps(eps)

    def energy(self, g, pi):
        F = self.F
        gm = 0.5 * (g[:-1] + g[1:])
        dg = np.diff(g) / self.h
        Ek = np.sum(self.r2 * (pi ** 2 / (2 * F(g)) + W_fun(g))) * self.h
        Es = 0.5 * np.sum(self.rm2 * F(gm) * dg ** 2) * self.h
        return Ek + Es

    def rhs(self, g, pi):
        F, Fp = self.F, self.Fp
        h, r2, rm2 = self.h, self.r2, self.rm2
        gm = 0.5 * (g[:-1] + g[1:])
        dg = np.diff(g) / h
        Fg = F(g)
        # dg_j/dt = pi_j / F_j
        gdot = pi / Fg
        # dpi_j/dt = -(1/(h r_j^2)) dH/dg_j
        dH = np.zeros_like(g)
        # local terms
        dH += h * r2 * (-pi ** 2 * Fp(g) / (2 * Fg ** 2) + W1_fun(g))
        # spatial terms: sum_m  h rm2_m [ F(gm_m) dg_m^2 / 2 ]
        # d/dg_j: m = j: -F(gm) dg/h + F'(gm)/2 * dg^2 / 2
        #         m = j-1: +F(gm) dg/h + F'(gm)/2 * dg^2 / 2
        t_flux = rm2 * F(gm) * dg          # length N-1
        t_quad = 0.25 * h * rm2 * Fp(gm) * dg ** 2
        dH[:-1] += -t_flux + t_quad
        dH[1:] += t_flux + t_quad
        pidot = -dH / (h * r2)
        return gdot, pidot

    def rk4_step(self, g, pi, dt):
        k1g, k1p = self.rhs(g, pi)
        k2g, k2p = self.rhs(g + 0.5 * dt * k1g, pi + 0.5 * dt * k1p)
        k3g, k3p = self.rhs(g + 0.5 * dt * k2g, pi + 0.5 * dt * k2p)
        k4g, k4p = self.rhs(g + dt * k3g, pi + dt * k3p)
        g2 = g + dt / 6 * (k1g + 2 * k2g + 2 * k3g + k4g)
        pi2 = pi + dt / 6 * (k1p + 2 * k2p + 2 * k3p + k4p)
        return g2, pi2


def run_case(dyn, g_eq, v_deep, amp, dt, r, bg_amp, t_end=T_END,
             sample_dt=0.02):
    """Returns (ts, a(t), ||dg||/bg, gate, E0, done, t_break).
    Breakdown detection (recorded, not interpreted here): the field g
    leaves the model domain (g <= 0: log undefined / non-finite) --
    the run stops at the LAST FINITE sample; t_break reported."""
    h = r[1] - r[0]
    r2h = r ** 2 * h
    g = g_eq + amp * v_deep
    pi = np.zeros_like(g)
    E0 = dyn.energy(g, pi)
    nstep = int(round(t_end / dt))
    every = max(1, int(round(sample_dt / dt)))
    ts, avals, dnorm = [], [], []
    emax = 0.0
    t_break = None
    for k in range(nstep + 1):
        if k % every == 0:
            if (not np.all(np.isfinite(g))) or float(np.min(g)) <= 0.0:
                t_break = k * dt
                break
            dg = g - g_eq
            ts.append(k * dt)
            avals.append(float(np.sum(r2h * v_deep * dg)))
            dnorm.append(float(np.max(np.abs(dg))) / bg_amp)
            E = dyn.energy(g, pi)
            if np.isfinite(E):
                emax = max(emax, abs(E - E0) / abs(E0))
        if k < nstep:
            g, pi = dyn.rk4_step(g, pi, dt)
    done = t_break is None
    return (np.array(ts), np.array(avals), np.array(dnorm), emax, E0,
            done, t_break)


def fit_growth(ts, avals, a0, cap_rel=20.0, floor_rel=3.0):
    """log-linear fit of |a(t)| in the window floor*|a0| < |a| < cap*|a0|,
    monotone segment from the start."""
    aa = np.abs(avals)
    lo, hi = floor_rel * abs(a0), cap_rel * abs(a0)
    idx = []
    for i in range(len(ts)):
        if aa[i] >= hi:
            break
        if aa[i] > lo:
            idx.append(i)
    if len(idx) < 5:
        return np.nan, (np.nan, np.nan)
    i0, i1 = idx[0], idx[-1]
    tt, yy = ts[i0:i1 + 1], np.log(aa[i0:i1 + 1])
    p = np.polyfit(tt, yy, 1)
    return float(p[0]), (float(ts[i0]), float(ts[i1]))


# =====================================================================
# RUN
# =====================================================================
print("=" * 78)
print("Phase 3: nonlinear dynamical test, M0-f_eps.  LOCK: Phase0 sec. 3.")
print("mu only (tau: NO EL representative in f_eps -- collapse, #62 W2b;")
print("OUT OF SCOPE of P3 per LOCK, recorded here explicitly).")
print("RK4 on the exact semi-discrete Hamiltonian system; gate |dE|/E<1e-6.")
print("=" * 78)

EPS_RUN = tuple(float(a) for a in sys.argv[1:]) or (0.2, 0.1)

for eps in EPS_RUN:
    label = "MAIN" if eps == 0.2 else "CONTROL"
    print("\n" + "-" * 78)
    print("[eps = %.1f, %s]" % (eps, label))
    r, g_eq, gp, d2, r_end = profile_soft(G0_MU, eps)
    if len(r) < N_GRID:
        print("  profil mu niekompletny (r_end = %.2f < %.1f) -- eps=%.1f"
              % (r_end, R_BOX, eps))
        print("  przypadek pomijany (zgodnie z #62: brak reprezentanta EL).")
        continue
    h = r[1] - r[0]
    bg_amp = float(np.max(np.abs(g_eq - 1.0)))
    print("  profil: min g = %.4f, max|g_eq - 1| = %.4f (tlo dla progu 10%%)"
          % (float(np.min(g_eq)), bg_amp))

    F, Fp, Fpp = make_feps(eps)
    d, e = build_tridiag_gen(r, g_eq, gp, d2, F, Fp, Fpp)
    vals, vecs = eigh_tridiagonal(d, e, select='i', select_range=(0, 5))
    lam_min = float(vals[0])
    # F-weighted (exact dynamical) eigenvalue: B^{-1/2} A B^{-1/2}
    sF = np.sqrt(F(g_eq))
    dw = d / F(g_eq)
    ew = e / (sF[:-1] * sF[1:])
    vals_w = eigh_tridiagonal(dw, ew, select='i', select_range=(0, 5),
                              eigvals_only=True)
    lam_w = float(vals_w[0])
    sig_lock = np.sqrt(LAM_LOCK)
    sig_lam = np.sqrt(-lam_min) if lam_min < 0 else np.nan
    sig_w = np.sqrt(-lam_w) if lam_w < 0 else np.nan
    print("  spektrum L (waga 1):   lam_min = %.4f  (siatka N=4000; #62: "
          "-1.389 dla eps=0.2)" % lam_min)
    print("  spektrum L (waga F):   lam_min = %.4f" % lam_w)
    print("  przewidywania liniowe: sigma_LOCK = sqrt(1.389) = %.4f;"
          % sig_lock)
    print("    sigma(waga 1) = %.4f;  sigma(waga F, dokladna dynamika) = %.4f"
          % (sig_lam, sig_w))

    # v_deep: weight-1 eigenvector (as #62), u -> v = u/r, ||v||_B = 1
    u = vecs[:, 0]
    v = u / r
    v /= np.sqrt(np.sum(r ** 2 * v ** 2) * h)
    if v[np.argmax(np.abs(v))] < 0:
        v = -v
    print("  v_deep: ||v||_B = 1, max|v| = %.4f przy r = %.2f"
          % (float(np.max(np.abs(v))), float(r[int(np.argmax(np.abs(v)))])))
    print("  UWAGA brzeg: zero-flux przy r=R=60; odbicie wraca do rdzenia")
    print("  ok. t ~ 2R = 120 -- rekurencje po tym czasie moga zawierac")
    print("  skladowa odbita od pudla (raportowane, nie interpretowane).")

    dyn = SoftWallDynamics(r, eps)

    for amp in AMPS:
        dts = (DT_BASE, DT_BASE / 2) if amp == 0.01 else (DT_BASE,)
        for dt in dts:
            ts, avals, dnorm, egate, E0, done, t_break = run_case(
                dyn, g_eq, v, amp, dt, r, bg_amp)
            sig_fit, win = fit_growth(ts, avals, avals[0])
            tag = "a=%+.2f dt=%.4f" % (amp, dt)
            print("\n  [%s]  E0 = %.4f, gate max|dE|/E = %.2e  -> %s%s"
                  % (tag, E0, egate,
                     "PASS" if egate < GATE else "FAIL",
                     "" if done else
                     "  (BREAKDOWN t*=%.2f: g <= 0 -- pole opuszcza"
                     " dziedzine modelu (log); wyjscie w skonczonym"
                     " czasie)" % t_break))
            if np.isfinite(sig_fit):
                dev_lock = 100 * (sig_fit - sig_lock) / sig_lock
                dev_w = 100 * (sig_fit - sig_w) / sig_w if np.isfinite(sig_w) \
                    else np.nan
                print("    wzrost: sigma_fit = %.4f (okno t=[%.1f,%.1f]);"
                      % (sig_fit, win[0], win[1]))
                print("      vs sigma_LOCK: %+.1f%%;  vs sigma(waga F): %+.1f%%"
                      % (dev_lock, dev_w))
            else:
                print("    wzrost: brak okna log-liniowego (|a| nie urosla "
                      "3x nad start)")
            # saturation / recurrence report
            imax = int(np.argmax(np.abs(avals)))
            amax = float(np.abs(avals)[imax])
            tmax = float(ts[imax])
            dmax = float(np.max(dnorm))
            idmax = int(np.argmax(dnorm))
            print("    max|a| = %.4f przy t = %.1f;  max||dg||_inf/tlo = "
                  "%.2f%% przy t = %.1f"
                  % (amax, tmax, 100 * dmax, float(ts[idmax])))
            # recurrence: does |a| drop below 50% of its max afterwards?
            rec = ""
            after = np.abs(avals)[imax:]
            if len(after) > 3 and np.min(after) < 0.5 * amax:
                tr = float(ts[imax + int(np.argmin(after))])
                rec = ("rekurencja: |a| spada do %.1f%% max przy t = %.1f"
                       % (100 * float(np.min(after)) / amax, tr))
            else:
                rec = "brak rekurencji (|a| nie wraca ponizej 50% max)"
            print("    %s" % rec)
            # sampled a(t): step 0.5 for short runs, 10 for long runs
            step = 0.5 if (not done or ts[-1] <= 20.0) else 10.0
            sel = [i for i in range(len(ts))
                   if abs(ts[i] / step - round(ts[i] / step)) < 1e-6]
            line = ", ".join("a(%.1f)=%+.3e" % (ts[i], avals[i])
                             for i in sel)
            print("    a(t) co %.1f: %s" % (step, line))
            line2 = ", ".join("%.1f:%.1f%%" % (ts[i], 100 * dnorm[i])
                              for i in sel)
            print("    ||dg||/tlo co %.1f [%%]: %s" % (step, line2))

print("\n" + "=" * 78)
print("WERDYKT V3 (kryteria LOCKED):")
print("  (i)  wzrost wykladniczy zgodny z lam_min (+-20%) i brak saturacji")
print("       ponizej ||dg|| ~ 10% tla -> 'niestabilnosc potwierdzona")
print("       nieliniowo w M0-f_eps';")
print("  (ii) saturacja < 10% + rekurencje -> 'efektywna stabilizacja")
print("       nieliniowa' (z analiza mechanizmu).")
print("  Interpretacja liczb powyzej -> README (kazdy kierunek = PASS")
print("  wykonania).  tau: poza zakresem P3 (brak EL w f_eps, #62).")
print("=" * 78)
