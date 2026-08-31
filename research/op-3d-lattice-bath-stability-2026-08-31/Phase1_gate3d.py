#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
op-3d-lattice-bath-stability (Phase 1) -- bramka maszynerii 3D.

LOCK: Phase0_balance.md sec. 3 Phase 1 (P1a-P1d); decyzje FROZEN:
Phase_method_decisions.md sec. 0-5.

P1a: dyspersja prozni 3D exact (pudlo L=2pi, N in {32,64}); operator
     kanoniczny (K=g^4 -> K(1)=1) i kontrola stabilna; zalockowane k:
     Gamma, X, M, R. Gate (metryka FROZEN): najnizsza galaz, blad
     |dw2|/max(|w2_ex|,1) <= 1e-3 przy N=32 ORAZ maxerr32/maxerr64
     in [3,5] na sektor. Deskryptywnie: 10 galezi.
P1b: kotwica radialna (model #63 M0-f_eps verbatim): lam_min(w1) =
     -1.3896 +- 1e-3; t* = 3.62 +- 0.05 (a=+-0.01, dt in {0.004,0.002}).
P1c: most radialny->kartezjanski (model #63): pudlo L=30, N in {76,100};
     |lam_min(3D,N76) - lam_rad| <= 5% |lam_rad| ORAZ poprawa przy N=100;
     t*_izo(3D) na N=76 (a=+-0.01, dt in {0.02,0.01}) = 3.62 +- 15%;
     t*_izo(3D) := min po biegach (zasila Phase 4).
P1d: gate energii ewolucji 3D (L=2pi, N=32, a=1e-3, t_end=4, dt in
     {0.004,0.002}); obie dynamiki (kanoniczna K_eps 0.2/0.1 i #63 F_eps
     0.2): |dE|/E <= 1e-6 oraz ||g_dt - g_dt/2||_inf(t_end) <= 1e-6.

KAZDY FAIL bramki => STOP (raport, bez Phase 2-4).
INPUT: g0_mu = 2.02117 (kalibracja mu #63), eps = 0.2 (kontrola 0.1),
beta = gamma = 1.
"""
import time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.linalg import eigh_tridiagonal
from scipy.interpolate import CubicSpline
from scipy.integrate import solve_ivp

BETA = 1.0
GAMMA = 1.0
PHI_GOLD = (1 + 5 ** 0.5) / 2
G0_MU = PHI_GOLD * 1.24915           # = 2.02117... INPUT (#63)
ALPHA = 2.0
EPS = 0.2
GATE_E = 1e-6
G_FLOOR = 0.01                       # korekta 1 poprzednika (bath-two-sectors)
LAM_ANCHOR = -1.3896                 # kotwica #63/P2a
TSTAR_ANCHOR = 3.62
RNG_SEED = 20260831

t0_wall = time.time()


def stamp(msg):
    print("[t=%7.1fs] %s" % (time.time() - t0_wall, msg), flush=True)


# ===================================================================
# model #63 M0-f_eps (verbatim konwencje)
# ===================================================================
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


F_EPS, FP_EPS, FPP_EPS = make_feps(EPS)

# akcja kanoniczna (Phase 2-4): K = g^4, U' = g^6 - g^7; K_eps w ewolucji
K_fun = lambda g: g ** 4
Kp_fun = lambda g: 4 * g ** 3
U_fun = lambda g: BETA * g ** 7 / 7 - GAMMA * g ** 8 / 8
U1_fun = lambda g: BETA * g ** 6 - GAMMA * g ** 7


def make_keps(eps):
    def s(g):
        kv = K_fun(g)
        return np.sqrt(kv ** 2 + eps ** 2)
    Ke = lambda g: 0.5 * (K_fun(g) + s(g))
    Kep = lambda g: 0.5 * Kp_fun(g) * (1 + K_fun(g) / s(g))
    return Ke, Kep


def profile_soft(g0, R, N):
    """#62/#63 soliton_profile_soft (omega=0), verbatim."""
    def rhs(r_, y):
        g, gp = y
        g = max(g, 1e-3)
        drv = W1_fun(g) - 0.5 * FP_EPS(g) * gp ** 2
        if r_ < 1e-10:
            return [gp, drv / (3 * F_EPS(g))]
        return [gp, drv / F_EPS(g) - 2 * gp / r_]
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
    return r_g, g_arr, gp_arr, d2, r_end, len(r_g) == N


def build_tridiag_gen(r, g, gp, d2):
    """#62/#63 assembly (waga 1 na u = r v), verbatim."""
    N = len(r)
    h = r[1] - r[0]
    F_nodes = F_EPS(g)
    F_mid = F_EPS(0.5 * (g[:-1] + g[1:]))
    Q = W2_fun(g) - 0.5 * FPP_EPS(g) * gp ** 2 \
        - FP_EPS(g) * (d2 + 2 * gp / r)
    r_mid = 0.5 * (r[:-1] + r[1:])
    a = r_mid ** 2 * F_mid
    diag = np.empty(N)
    diag[0] = a[0] / h ** 2
    diag[1:-1] = (a[:-1] + a[1:]) / h ** 2
    diag[-1] = (a[-1] + r[-1] ** 2 * F_nodes[-1]) / h ** 2
    diag += r ** 2 * Q
    off = -a / h ** 2
    return diag / r ** 2, off / (r[:-1] * r[1:])


class RadialDynamics:
    """#63 SoftWallDynamics verbatim (F_eps, W)."""
    def __init__(self, r):
        self.r = r
        self.h = r[1] - r[0]
        self.r2 = r ** 2
        self.rm2 = (0.5 * (r[:-1] + r[1:])) ** 2

    def energy(self, g, pi):
        gm = 0.5 * (g[:-1] + g[1:])
        dg = np.diff(g) / self.h
        Ek = np.sum(self.r2 * (pi ** 2 / (2 * F_EPS(g)) + W_fun(g))) * self.h
        Es = 0.5 * np.sum(self.rm2 * F_EPS(gm) * dg ** 2) * self.h
        return Ek + Es

    def rhs(self, g, pi):
        h, r2, rm2 = self.h, self.r2, self.rm2
        gm = 0.5 * (g[:-1] + g[1:])
        dg = np.diff(g) / h
        Fg = F_EPS(g)
        gdot = pi / Fg
        dH = np.zeros_like(g)
        dH += h * r2 * (-pi ** 2 * FP_EPS(g) / (2 * Fg ** 2) + W1_fun(g))
        t_flux = rm2 * F_EPS(gm) * dg
        t_quad = 0.25 * h * rm2 * FP_EPS(gm) * dg ** 2
        dH[:-1] += -t_flux + t_quad
        dH[1:] += t_flux + t_quad
        pidot = -dH / (h * r2)
        return gdot, pidot

    def rk4(self, g, pi, dt):
        k1g, k1p = self.rhs(g, pi)
        k2g, k2p = self.rhs(g + 0.5 * dt * k1g, pi + 0.5 * dt * k1p)
        k3g, k3p = self.rhs(g + 0.5 * dt * k2g, pi + 0.5 * dt * k2p)
        k4g, k4p = self.rhs(g + dt * k3g, pi + dt * k3p)
        return (g + dt / 6 * (k1g + 2 * k2g + 2 * k3g + k4g),
                pi + dt / 6 * (k1p + 2 * k2p + 2 * k3p + k4p))


def evolve_radial(r, g_init, dt, t_end):
    dyn = RadialDynamics(r)
    g = g_init.copy()
    pi = np.zeros_like(g)
    E0 = dyn.energy(g, pi)
    nstep = int(round(t_end / dt))
    every = max(1, int(round(0.02 / dt)))
    emax, t_break = 0.0, None
    for k in range(nstep + 1):
        if (not np.all(np.isfinite(g))) or float(np.min(g)) <= G_FLOOR:
            t_break = k * dt
            break
        if k % every == 0:
            E = dyn.energy(g, pi)
            if np.isfinite(E):
                emax = max(emax, abs(E - E0) / abs(E0))
        if k < nstep:
            g, pi = dyn.rk4(g, pi, dt)
    return t_break, emax


# ===================================================================
# operator 3D (sparse, 7-punktowy stencil; FROZEN method_decisions sec.1-2)
# ===================================================================
def build_op3d(g3, Phi, Q3, h, phases=(1.0, 1.0, 1.0), weight3=None):
    """A = -div(Phi(g) grad) + Q, waga opcjonalna (H = B^-1/2 A B^-1/2).
    phases: e^{i k_a d} in {+1,-1} (zalockowane punkty k -> rzeczywiste)."""
    N = g3.shape[0]
    n = N ** 3
    idx3 = np.arange(n).reshape(N, N, N)
    diag = Q3.astype(float).ravel().copy()
    rows, cols, vals = [], [], []
    for ax in range(3):
        gn = np.roll(g3, -1, axis=ax)
        c = Phi(0.5 * (g3 + gn)) / h ** 2
        diag += (c + np.roll(c, 1, axis=ax)).ravel()
        ph = np.ones_like(g3)
        sl = [slice(None)] * 3
        sl[ax] = N - 1
        ph[tuple(sl)] = phases[ax]
        v = (-c * ph).ravel()
        r_ = idx3.ravel()
        c_ = np.roll(idx3, -1, axis=ax).ravel()
        rows.append(r_)
        cols.append(c_)
        vals.append(v)
        rows.append(c_)
        cols.append(r_)
        vals.append(v)
    rows.append(np.arange(n))
    cols.append(np.arange(n))
    vals.append(diag)
    A = sp.coo_matrix((np.concatenate(vals),
                       (np.concatenate(rows), np.concatenate(cols))),
                      shape=(n, n)).tocsr()
    if weight3 is not None:
        dinv = 1.0 / np.sqrt(weight3.ravel())
        D = sp.diags(dinv)
        A = D @ A @ D
    return A


_RNG = np.random.default_rng(RNG_SEED)
_V0_CACHE = {}


def v0_det(n):
    if n not in _V0_CACHE:
        v = np.random.default_rng(RNG_SEED).standard_normal(n)
        _V0_CACHE[n] = v / np.linalg.norm(v)
    return _V0_CACHE[n]


def lowest_eigs(A, nev=10, vecs=False):
    """eigsh 'SA' (FROZEN: tol=1e-6, ncv=80, maxiter=200000, v0 det.)."""
    n = A.shape[0]
    kw = dict(k=nev, which='SA', tol=1e-6, ncv=min(80, n - 1),
              maxiter=200000, v0=v0_det(n))
    if vecs:
        vals, V = spla.eigsh(A, **kw)
        o = np.argsort(vals)
        return vals[o], V[:, o]
    vals = spla.eigsh(A, return_eigenvectors=False, **kw)
    return np.sort(vals)


def exact_vacuum(kvec, d, m2, nb=10):
    G = 2 * np.pi / d * np.arange(-6, 7)
    KX, KY, KZ = np.meshgrid(kvec[0] + G, kvec[1] + G, kvec[2] + G,
                             indexing='ij')
    w2 = np.sort((KX ** 2 + KY ** 2 + KZ ** 2).ravel())[:nb] + m2
    return w2


def kpoints(d):
    q = np.pi / d
    return [("Gamma", (0.0, 0.0, 0.0), (1, 1, 1)),
            ("X", (q, 0.0, 0.0), (-1, 1, 1)),
            ("M", (q, q, 0.0), (-1, -1, 1)),
            ("R", (q, q, q), (-1, -1, -1))]


# ===================================================================
# dynamika 3D (semi-dyskretny hamiltonian, RK4; FROZEN sec. 5)
# ===================================================================
class Dyn3D:
    def __init__(self, h, Phi, Phip, P, P1):
        self.h = h
        self.Phi, self.Phip, self.P, self.P1 = Phi, Phip, P, P1

    def energy(self, g, pi):
        h = self.h
        E = np.sum(pi ** 2 / (2 * self.Phi(g)) + self.P(g)) * h ** 3
        for ax in range(3):
            gn = np.roll(g, -1, axis=ax)
            dg = (gn - g) / h
            E += 0.5 * np.sum(self.Phi(0.5 * (g + gn)) * dg ** 2) * h ** 3
        return E

    def rhs(self, g, pi):
        h = self.h
        Fg = self.Phi(g)
        gdot = pi / Fg
        dH = -pi ** 2 * self.Phip(g) / (2 * Fg ** 2) + self.P1(g)
        for ax in range(3):
            gn = np.roll(g, -1, axis=ax)
            gm = 0.5 * (g + gn)
            dg = (gn - g) / h
            t_flux = self.Phi(gm) * dg / h
            t_quad = 0.25 * self.Phip(gm) * dg ** 2
            dH += -t_flux + t_quad
            dH += np.roll(t_flux + t_quad, 1, axis=ax)
        return gdot, -dH

    def rk4(self, g, pi, dt):
        k1g, k1p = self.rhs(g, pi)
        k2g, k2p = self.rhs(g + 0.5 * dt * k1g, pi + 0.5 * dt * k1p)
        k3g, k3p = self.rhs(g + 0.5 * dt * k2g, pi + 0.5 * dt * k2p)
        k4g, k4p = self.rhs(g + dt * k3g, pi + dt * k3p)
        return (g + dt / 6 * (k1g + 2 * k2g + 2 * k3g + k4g),
                pi + dt / 6 * (k1p + 2 * k2p + 2 * k3p + k4p))


def evolve3d(dyn, g_init, dt, t_end, detect_floor=True):
    """Zwraca (t_break, egate, g_final)."""
    g = g_init.copy()
    pi = np.zeros_like(g)
    E0 = dyn.energy(g, pi)
    nstep = int(round(t_end / dt))
    every = max(1, int(round(0.02 / dt)))
    emax, t_break = 0.0, None
    for k in range(nstep + 1):
        if detect_floor and ((not np.all(np.isfinite(g)))
                             or float(np.min(g)) <= G_FLOOR):
            t_break = k * dt
            break
        if k % every == 0:
            E = dyn.energy(g, pi)
            if np.isfinite(E):
                emax = max(emax, abs(E - E0) / abs(E0))
        if k < nstep:
            g, pi = dyn.rk4(g, pi, dt)
    return t_break, emax, g


def minimg_r(N, h, L):
    """Odleglosc minimalnego obrazu od centrum L/2 (siatka wezlowa)."""
    x = np.arange(N) * h
    dx = np.mod(x - L / 2 + L / 2, L) - L / 2
    DX, DY, DZ = np.meshgrid(dx, dx, dx, indexing='ij')
    return np.sqrt(DX ** 2 + DY ** 2 + DZ ** 2)


# =====================================================================
# RUN
# =====================================================================
print("=" * 78)
print("Phase 1 -- bramka maszynerii 3D. LOCK sec. 3 Phase 1; decyzje FROZEN")
print("w Phase_method_decisions.md. INPUT: g0_mu=%.5f, eps=0.2 (kontrola")
print("0.1), beta=gamma=1. KAZDY FAIL => STOP bez Phase 2-4." % ())
print("=" * 78)
print("g0_mu = %.5f (INPUT #63)" % G0_MU)

gate_all = {}

# ------------------------------------------------------------------ P1a
stamp("[P1a] dyspersja prozni 3D exact: pudlo L=2pi, N in {32,64}")
L1 = 2 * np.pi
maxerr = {}
for sector, m2 in (("tach", -GAMMA), ("stab", +GAMMA)):
    for N in (32, 64):
        h = L1 / N
        g3 = np.ones((N, N, N))
        if sector == "tach":
            Q3 = K_fun(g3) * (2 * BETA * g3 - 3 * GAMMA * g3 ** 2)
        else:
            Q3 = -K_fun(g3) * (2 * BETA * g3 - 3 * GAMMA * g3 ** 2)
        # (|grad g|^2 = 0 na prozni)
        me_low, table = 0.0, []
        for kname, kv, phs in kpoints(L1):
            A = build_op3d(g3, K_fun, Q3, h, phases=phs, weight3=K_fun(g3))
            vals = lowest_eigs(A, nev=10)
            ex = exact_vacuum(kv, L1, m2, nb=10)
            err = np.abs(vals - ex) / np.maximum(np.abs(ex), 1.0)
            me_low = max(me_low, float(err[0]))
            table.append((kname, float(vals[0]), float(ex[0]),
                          float(err[0]), float(np.max(err))))
        maxerr[(sector, N)] = me_low
        stamp("  sektor %s N=%2d: maxerr(najnizsza galaz po k) = %.3e"
              % (sector, N, me_low))
        for kname, vn, ve, e0, e10 in table:
            print("    k=%-5s: w2_num=%+.6f w2_ex=%+.6f err=%.2e "
                  "(10 galezi: maxerr=%.2e)" % (kname, vn, ve, e0, e10))
p1a_pass = True
for sector in ("tach", "stab"):
    e32, e64 = maxerr[(sector, 32)], maxerr[(sector, 64)]
    ratio = e32 / e64 if e64 > 0 else np.inf
    ok = (e32 <= 1e-3) and (3.0 <= ratio <= 5.0)
    print("  P1a %s: maxerr(N32)=%.3e (gate<=1e-3), ratio 32/64=%.3f "
          "(gate [3,5]) -> %s" % (sector, e32, ratio,
                                  "PASS" if ok else "FAIL"))
    p1a_pass = p1a_pass and ok
gate_all["P1a"] = p1a_pass

# ------------------------------------------------------------------ P1b
stamp("[P1b] kotwica radialna (model #63 verbatim): R=60, N=4000")
r, g_eq, gp, d2, r_end, full = profile_soft(G0_MU, 60.0, 4000)
assert full, "profil mu niekompletny -- STOP"
dgr, egr = build_tridiag_gen(r, g_eq, gp, d2)
vals1, vecs1 = eigh_tridiagonal(dgr, egr, select='i', select_range=(0, 5))
lam_rad = float(vals1[0])
print("  lam_min(w1) = %.5f (kotwica %.4f +- 1e-3)" % (lam_rad, LAM_ANCHOR))
ok_lam = abs(lam_rad - LAM_ANCHOR) <= 1e-3
u = vecs1[:, 0]
v_rad = u / r
hr = r[1] - r[0]
v_rad /= np.sqrt(np.sum(r ** 2 * v_rad ** 2) * hr)
if v_rad[np.argmax(np.abs(v_rad))] < 0:
    v_rad = -v_rad
tstars = []
ok_evo = True
for amp, dt in ((0.01, 0.004), (0.01, 0.002), (-0.01, 0.004)):
    tb, egate = evolve_radial(r, g_eq + amp * v_rad, dt, 12.0)
    okg = egate < GATE_E
    okt = tb is not None and abs(tb - TSTAR_ANCHOR) <= 0.05
    ok_evo = ok_evo and okg and okt
    if tb is not None:
        tstars.append(tb)
    stamp("  a=%+.2f dt=%.4f: gate|dE|/E=%.2e %s; t*=%s (kotwica 3.62"
          "+-0.05) %s" % (amp, dt, egate, "PASS" if okg else "FAIL",
                          ("%.2f" % tb) if tb is not None else "BRAK",
                          "PASS" if okt else "FAIL"))
tstar_rad = min(tstars) if tstars else None
p1b_pass = ok_lam and ok_evo
print("  P1b: lam_min %s; ewolucja %s -> %s (t*_radial = %s)"
      % ("PASS" if ok_lam else "FAIL", "PASS" if ok_evo else "FAIL",
         "PASS" if p1b_pass else "FAIL",
         ("%.2f" % tstar_rad) if tstar_rad else "BRAK"))
gate_all["P1b"] = p1b_pass

# ------------------------------------------------------------------ P1c
stamp("[P1c] most radialny->kartezjanski (model #63): L=30, N in {76,100}")
L3 = 30.0
cs_g = CubicSpline(r, g_eq - 1.0)
cs_gp = CubicSpline(r, gp)
cs_d2 = CubicSpline(r, d2)
d2_0 = float(cs_d2(0.0))


def radial_fields(rr):
    """(g, gp, lap) w odleglosci rr (zero-padding do prozni za r_end)."""
    inside = rr <= r[-1]
    g3 = np.where(inside, 1.0 + cs_g(np.minimum(rr, r[-1])), 1.0)
    gp3 = np.where(inside, cs_gp(np.minimum(rr, r[-1])), 0.0)
    rr_safe = np.maximum(rr, 1e-12)
    lap3 = np.where(inside, cs_d2(np.minimum(rr, r[-1]))
                    + 2 * gp3 / rr_safe, 0.0)
    lap3 = np.where(rr < 1e-9, 3 * d2_0, lap3)
    return g3, gp3, lap3


lam3 = {}
v_mode76 = None
g3_76 = None
for N in (76, 100):
    h = L3 / N
    rr = minimg_r(N, h, L3)
    g3, gp3, lap3 = radial_fields(rr)
    Q3 = W2_fun(g3) - 0.5 * FPP_EPS(g3) * gp3 ** 2 - FP_EPS(g3) * lap3
    A = build_op3d(g3, F_EPS, Q3, h)          # waga 1 (model #63)
    vals, V = lowest_eigs(A, nev=10, vecs=True)
    lam3[N] = float(vals[0])
    stamp("  N=%3d (h=%.4f): lam_min(3D) = %.5f; 10 najnizszych: %s"
          % (N, h, lam3[N], ["%+.4f" % x for x in vals]))
    if N == 76:
        vv = V[:, 0].reshape(N, N, N)
        vv = vv / np.sqrt(np.sum(vv ** 2) * h ** 3)
        if vv.ravel()[np.argmax(np.abs(vv))] < 0:
            vv = -vv
        v_mode76 = vv
        g3_76 = g3
dev76 = abs(lam3[76] - lam_rad)
dev100 = abs(lam3[100] - lam_rad)
ok_l76 = dev76 <= 0.05 * abs(lam_rad)
ok_impr = dev100 < dev76
print("  |lam(N76)-lam_rad| = %.4f (gate <= %.4f): %s; poprawa N100: "
      "%.4f < %.4f: %s" % (dev76, 0.05 * abs(lam_rad),
                           "PASS" if ok_l76 else "FAIL", dev100, dev76,
                           "PASS" if ok_impr else "FAIL"))
# t*_izo(3D): ewolucja #63-3D na N=76
h76 = L3 / 76
dyn63 = Dyn3D(h76, F_EPS, FP_EPS, W_fun, W1_fun)
t3s = []
ok_t3 = True
for amp, dt in ((0.01, 0.02), (0.01, 0.01), (-0.01, 0.02)):
    tb, egate, _ = evolve3d(dyn63, g3_76 + amp * v_mode76, dt, 8.0)
    okg = egate < GATE_E
    okt = tb is not None and abs(tb - TSTAR_ANCHOR) <= 0.15 * TSTAR_ANCHOR
    ok_t3 = ok_t3 and okg and okt
    if tb is not None:
        t3s.append(tb)
    stamp("  evo3D a=%+.2f dt=%.3f: gate=%.2e %s; t*=%s (3.62 +- 15%%:"
          " [3.08,4.16]) %s" % (amp, dt, egate, "PASS" if okg else "FAIL",
                                ("%.2f" % tb) if tb is not None else "BRAK",
                                "PASS" if okt else "FAIL"))
tstar_izo3d = min(t3s) if t3s else None
p1c_pass = ok_l76 and ok_impr and ok_t3
print("  P1c: %s (t*_izo(3D) = %s -> Phase 4)"
      % ("PASS" if p1c_pass else "FAIL",
         ("%.2f" % tstar_izo3d) if tstar_izo3d else "BRAK"))
gate_all["P1c"] = p1c_pass

# ------------------------------------------------------------------ P1d
stamp("[P1d] gate energii ewolucji 3D: L=2pi, N=32, a=1e-3, t_end=4")
N = 32
h = L1 / N
x = np.arange(N) * h
X3, Y3, Z3 = np.meshgrid(x, x, x, indexing='ij')
g_init = 1.0 + 1e-3 * (0.5 + np.cos(X3) + np.cos(Y3) + np.cos(Z3))
p1d_pass = True
variants = []
Ke2, Kep2 = make_keps(0.2)
Ke1, Kep1 = make_keps(0.1)
variants.append(("kanoniczna K_eps=0.2", Dyn3D(h, Ke2, Kep2, U_fun, U1_fun)))
variants.append(("kanoniczna K_eps=0.1", Dyn3D(h, Ke1, Kep1, U_fun, U1_fun)))
variants.append(("#63 F_eps=0.2", Dyn3D(h, F_EPS, FP_EPS, W_fun, W1_fun)))
for name, dyn in variants:
    gfin = {}
    for dt in (0.004, 0.002):
        tb, egate, gf = evolve3d(dyn, g_init, dt, 4.0)
        gfin[dt] = gf
        okg = (egate < GATE_E) and tb is None
        p1d_pass = p1d_pass and okg
        stamp("  [%s] dt=%.3f: gate|dE|/E=%.2e %s%s"
              % (name, dt, egate, "PASS" if okg else "FAIL",
                 "" if tb is None else " BREAKDOWN t=%.2f" % tb))
    dcons = float(np.max(np.abs(gfin[0.004] - gfin[0.002])))
    okc = dcons <= 1e-6
    p1d_pass = p1d_pass and okc
    print("    zbieznosc dt: ||g_dt - g_dt/2||_inf(t=4) = %.2e "
          "(gate<=1e-6) %s" % (dcons, "PASS" if okc else "FAIL"))
gate_all["P1d"] = p1d_pass

# ------------------------------------------------------------------ werdykt
print("\n" + "=" * 78)
print("BRAMKA Phase 1 (LOCK sec. 3; kazdy FAIL => STOP):")
for k in ("P1a", "P1b", "P1c", "P1d"):
    print("  %s: %s" % (k, "PASS" if gate_all[k] else "FAIL"))
allpass = all(gate_all.values())
print("\nP1 GATE: %s" % ("PASS -> Phase 2" if allpass
                         else "FAIL -> STOP (raport bez Phase 2-4)"))
if allpass:
    print("lam_rad = %.5f; t*_radial = %.2f; t*_izo(3D) = %.2f "
          "(zasila Phase 4: t_end = 3 t*_izo = %.2f)"
          % (lam_rad, tstar_rad, tstar_izo3d, 3 * tstar_izo3d))
print("=" * 78)
