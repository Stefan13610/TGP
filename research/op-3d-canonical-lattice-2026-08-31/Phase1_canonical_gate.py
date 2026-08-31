#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
op-3d-canonical-lattice (Phase 1) -- nowa kotwica kanoniczna + most 3D.

LOCK: Phase0_balance.md sec. 3 Phase 1 (P1b-kan, P1c-kan; P1a/P1d
DZIEDZICZONE); decyzje FROZEN: Phase_method_decisions.md sec. 0-2.
Baza kodu: ../op-3d-lattice-bath-stability-2026-08-31/Phase1_gate3d.py
(operator 3D, dynamika, interpolacja radialna->3D verbatim); eigsh tol=0
(korekta 1 poprzednika, przyjeta od startu).

P1b-kan: profil radialny mu kanoniczny (EL: g''+(2/r)g'+(2/g)g'^2 =
    g^2(1-g); INPUT g0_mu = phi*0.90548 = 1.4650974); residuum
    (operacjonalizacja FROZEN): ||g(rtol1e-12)-g(rtol1e-13)||_inf <=
    1e-10; kotwica = lam_min(w1, u-forma: -lap + (4g-5g^2)) przy
    h=0.0125, R=60; gate wewnetrzny |lam(0.025)-lam(0.0125)| <=
    1e-3*|lam|; OBOWIAZKOWO: diagnostyka kieszeni V(r) (min, FWHM);
    t*_ref: ewolucja radialna K_eps=0.2 (kontrola 0.1), a=+-0.01,
    dt in {0.004,0.002}, dt-konsystencja <= 1%; detektor t* (FROZEN):
    min g <= 0.01 LUB max g >= 2 g0 LUB niefinitycznosc.
P1c-kan (bramka glowna): operator wagi-K (forma Phase 3) wokol profilu
    interpolowanego, L=30, N in {76,100}: |lam_min(N100) - kotwica| <=
    5%|kotwica| ORAZ |lam(N76)-kotwica| > |lam(N100)-kotwica|;
    t*_izo(3D, N=100) = t*_ref +- 15% (kazdy bieg; a=+-0.01,
    dt in {0.02,0.01}); t*_izo := min po biegach (zasila Phase 4).
P1a/P1d: DZIEDZICZONE od poprzednika (ten sam model kanoniczny,
    ten sam integrator, to samo srodowisko) -- cytat liczb w output.

FAIL P1c-kan => STOP (raport bez Phase 2-4).
INPUT: g0_e=0.90548, phi (kalibracja mu), eps=0.2 (kontrola 0.1),
beta=gamma=1.
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
G0_E = 0.90548                       # INPUT (r21/phi-FP, bath-two-sectors)
G0_MU = PHI_GOLD * G0_E              # = 1.4650974 (definicja INPUT; MD)
EPS = 0.2
G_FLOOR = 0.01
G_CEIL = 2.0 * G0_MU                 # detektor runaway g->duze (MD sec.1)
RNG_SEED = 20260831
R_RAD = 60.0
H_MASTER = 0.0125
H_LIST = (0.05, 0.025, 0.0125)
L3 = 30.0

t0_wall = time.time()


def stamp(msg):
    print("[t=%7.1fs] %s" % (time.time() - t0_wall, msg), flush=True)


# ===================================================================
# akcja kanoniczna: K = g^4, U = g^7/7 - g^8/8 (beta=gamma=1)
# ===================================================================
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


def profile_canonical(g0, R, h_master, rtol, atol):
    """Profil radialny kanoniczny: g'' = g^2(1-g) - (2/g)g'^2 - (2/r)g'."""
    def rhs(r_, y):
        g, gp = y
        g = max(g, 1e-3)
        drv = g ** 2 * (BETA - GAMMA * g) - (2.0 / g) * gp ** 2
        if r_ < 1e-10:
            return [gp, drv / 3.0]
        return [gp, drv - 2 * gp / r_]
    N = int(round(R / h_master))
    r_grid = (np.arange(N) + 0.5) * h_master
    sol = solve_ivp(rhs, [1e-6, R + h_master], [g0, 0.0], method='DOP853',
                    max_step=0.02, rtol=rtol, atol=atol, dense_output=True)
    assert float(sol.t[-1]) >= r_grid[-1], "profil niekompletny -- STOP"
    vals = sol.sol(r_grid)
    g_arr, gp_arr = vals[0], vals[1]
    d2 = np.array([rhs(r_grid[j], [g_arr[j], gp_arr[j]])[1]
                   for j in range(N)])
    return r_grid, g_arr, gp_arr, d2, sol


def uform_lam_radial(sol, h, R, nev=6):
    """u-forma w1 (FROZEN MD sec.1): -u'' + V u = lam u, V = 4g-5g^2,
    u = r*chi, siatka cell-centered, ghost u(-h/2)=-u(h/2)."""
    N = int(round(R / h))
    r = (np.arange(N) + 0.5) * h
    g = sol.sol(r)[0]
    V = 4 * BETA * g - 5 * GAMMA * g ** 2
    diag = np.full(N, 2.0 / h ** 2) + V
    diag[0] = 3.0 / h ** 2 + V[0]
    off = np.full(N - 1, -1.0 / h ** 2)
    vals, vecs = eigh_tridiagonal(diag, off, select='i',
                                  select_range=(0, nev - 1))
    return vals, vecs, r, g


def wK_lam_radial(sol, h, R, nev=3):
    """Deskryptywnie: radialna forma wagi-K (FV; kontrola rownowaznosci).
    -(1/r^2)(r^2 K v')' + Q v = lam K v, Q = K(2g-3g^2+2g'^2/g^2)."""
    N = int(round(R / h))
    r = (np.arange(N) + 0.5) * h
    vv = sol.sol(r)
    g, gp = vv[0], vv[1]
    K = K_fun(g)
    Q = K * (2 * BETA * g - 3 * GAMMA * g ** 2 + 2 * gp ** 2 / g ** 2)
    r_mid = 0.5 * (r[:-1] + r[1:])
    g_mid = 0.5 * (g[:-1] + g[1:])
    a = r_mid ** 2 * K_fun(g_mid)
    diag = np.empty(N)
    diag[0] = a[0] / h ** 2
    diag[1:-1] = (a[:-1] + a[1:]) / h ** 2
    diag[-1] = (a[-1] + r[-1] ** 2 * K[-1]) / h ** 2
    diag += r ** 2 * Q
    off = -a / h ** 2
    w = r ** 2 * K
    diag_s = diag / w
    off_s = off / np.sqrt(w[:-1] * w[1:])
    vals, _ = eigh_tridiagonal(diag_s, off_s, select='i',
                               select_range=(0, nev - 1))
    return vals


class RadialDynamicsCan:
    """RadialDynamics poprzednika z (F_eps,W)->(K_eps,U) (struktura
    hamiltonowska verbatim; MD sec.1)."""
    def __init__(self, r, Ke, Kep):
        self.r = r
        self.h = r[1] - r[0]
        self.r2 = r ** 2
        self.rm2 = (0.5 * (r[:-1] + r[1:])) ** 2
        self.Ke, self.Kep = Ke, Kep

    def energy(self, g, pi):
        gm = 0.5 * (g[:-1] + g[1:])
        dg = np.diff(g) / self.h
        Ek = np.sum(self.r2 * (pi ** 2 / (2 * self.Ke(g)) + U_fun(g))) * self.h
        Es = 0.5 * np.sum(self.rm2 * self.Ke(gm) * dg ** 2) * self.h
        return Ek + Es

    def rhs(self, g, pi):
        h, r2, rm2 = self.h, self.r2, self.rm2
        gm = 0.5 * (g[:-1] + g[1:])
        dg = np.diff(g) / h
        Kg = self.Ke(g)
        gdot = pi / Kg
        dH = np.zeros_like(g)
        dH += h * r2 * (-pi ** 2 * self.Kep(g) / (2 * Kg ** 2) + U1_fun(g))
        t_flux = rm2 * self.Ke(gm) * dg
        t_quad = 0.25 * h * rm2 * self.Kep(gm) * dg ** 2
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


def evolve_radial(dyn, g_init, dt, t_end):
    """Detektor t* FROZEN (MD sec.1): floor 0.01 / ceil 2 g0 / non-finite."""
    g = g_init.copy()
    pi = np.zeros_like(g)
    E0 = dyn.energy(g, pi)
    nstep = int(round(t_end / dt))
    every = max(1, int(round(0.02 / dt)))
    emax, t_break = 0.0, None
    for k in range(nstep + 1):
        if (not np.all(np.isfinite(g))) or float(np.min(g)) <= G_FLOOR \
                or float(np.max(g)) >= G_CEIL:
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
# operator 3D + dynamika 3D (verbatim Phase1_gate3d.py poprzednika)
# ===================================================================
def build_op3d(g3, Phi, Q3, h, phases=(1.0, 1.0, 1.0), weight3=None):
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
        rows += [r_, c_]
        cols += [c_, r_]
        vals += [v, v]
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


_V0_CACHE = {}


def v0_det(n):
    if n not in _V0_CACHE:
        v = np.random.default_rng(RNG_SEED).standard_normal(n)
        _V0_CACHE[n] = v / np.linalg.norm(v)
    return _V0_CACHE[n]


def lowest_eigs(A, nev=10, vecs=False):
    """eigsh 'SA'; tol=0 (korekta 1 poprzednika -- przyjeta od startu)."""
    n = A.shape[0]
    kw = dict(k=nev, which='SA', tol=0, ncv=min(80, n - 1),
              maxiter=200000, v0=v0_det(n))
    if vecs:
        vals, V = spla.eigsh(A, **kw)
        o = np.argsort(vals)
        return vals[o], V[:, o]
    return np.sort(spla.eigsh(A, return_eigenvectors=False, **kw))


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


def evolve3d(dyn, g_init, dt, t_end):
    g = g_init.copy()
    pi = np.zeros_like(g)
    E0 = dyn.energy(g, pi)
    nstep = int(round(t_end / dt))
    every = max(1, int(round(0.02 / dt)))
    emax, t_break = 0.0, None
    for k in range(nstep + 1):
        if (not np.all(np.isfinite(g))) or float(np.min(g)) <= G_FLOOR \
                or float(np.max(g)) >= G_CEIL:
            t_break = k * dt
            break
        if k % every == 0:
            E = dyn.energy(g, pi)
            if np.isfinite(E):
                emax = max(emax, abs(E - E0) / abs(E0))
        if k < nstep:
            g, pi = dyn.rk4(g, pi, dt)
    return t_break, emax


def minimg_r(N, h, L):
    x = np.arange(N) * h
    dx = np.mod(x - L / 2 + L / 2, L) - L / 2
    DX, DY, DZ = np.meshgrid(dx, dx, dx, indexing='ij')
    return np.sqrt(DX ** 2 + DY ** 2 + DZ ** 2)


# =====================================================================
# RUN
# =====================================================================
print("=" * 78)
print("Phase 1 (op-3d-canonical-lattice) -- kotwica kanoniczna + most 3D.")
print("LOCK sec. 3 Phase 1; decyzje FROZEN: Phase_method_decisions.md.")
print("INPUT: g0_e=0.90548, g0_mu=phi*g0_e=%.7f, eps=0.2 (kontrola 0.1),"
      % G0_MU)
print("d*1=3.0790 [nieuzywane w P1], beta=gamma=1. eigsh tol=0 (korekta 1).")
print("=" * 78)

# --------------------------------------------------- P1a/P1d DZIEDZICZONE
print("""
[DZIEDZICZENIE P1a/P1d -- LOCK sec.1; NIE liczone ponownie]
Zrodlo: ../op-3d-lattice-bath-stability-2026-08-31/Phase1_output.txt
(ten sam model kanoniczny, ten sam integrator RK4/Dyn3D verbatim,
to samo srodowisko: CPython 3.14.2, numpy 2.4.3, scipy 1.17.1).
P1a (dyspersja prozni 3D exact, operator kanoniczny K(1)=1): PASS
  maxerr(najnizsza galaz, N=32) = 6.022e-04 (tach) / 3.441e-04 (stab)
  (gate <= 1e-3); ratio N32/N64 = 3.999 / 3.999 (gate [3,5], rzad 2).
P1d (gate energii ewolucji 3D, L=2pi, N=32, a=1e-3, t_end=4): PASS
  kanoniczna K_eps=0.2: |dE|/E <= 6.82e-15 (dt=0.004), 4.01e-16 (0.002),
    dt-konsystencja 1.23e-13;
  kanoniczna K_eps=0.1: |dE|/E <= 7.42e-15 / 8.02e-16, dt-kons. 1.29e-13
  (gate <= 1e-6 wszedzie).
""")
gate_all = {"P1a(dziedz.)": True, "P1d(dziedz.)": True}

# ------------------------------------------------------------------ P1b-kan
stamp("[P1b-kan] profil radialny mu kanoniczny (g0=%.7f), R=60" % G0_MU)
r_m, g_m, gp_m, d2_m, sol12 = profile_canonical(
    G0_MU, R_RAD, H_MASTER, rtol=1e-12, atol=1e-14)
_, g_m13, _, _, _ = profile_canonical(
    G0_MU, R_RAD, H_MASTER, rtol=1e-13, atol=1e-15)
res_ode = float(np.max(np.abs(g_m - g_m13)))
ok_res = res_ode <= 1e-10
print("  residuum (FROZEN: ||g(1e-12)-g(1e-13)||_inf) = %.2e (gate<=1e-10)"
      " -> %s" % (res_ode, "PASS" if ok_res else "FAIL"))
gate_all["P1b-residuum"] = ok_res
# deskryptywnie: residuum EL ze splajnu g' (informacyjnie, nosi blad splajnu)
gpp_spl = CubicSpline(r_m, gp_m)(r_m, 1)
res_el = np.abs(gpp_spl + 2 * gp_m / r_m + 2 * gp_m ** 2 / g_m
                - g_m ** 2 * (1 - g_m))
print("  residuum EL (deskryptywnie, g'' ze splajnu g'): max=%.2e"
      % float(np.max(res_el)))
print("  profil: g(0+)=%.7f, g in [%.6f, %.6f], g(R)=%.8f"
      % (g_m[0], float(np.min(g_m)), float(np.max(g_m)), g_m[-1]))

# kotwica lam_min(w1) na h in {0.05, 0.025, 0.0125}
stamp("[P1b-kan] kotwica lam_min(w1, u-forma) na h in {0.05,0.025,0.0125}")
lam_h = {}
vecs_h = {}
for h in H_LIST:
    vals, vecs, r_h, g_h = uform_lam_radial(sol12, h, R_RAD)
    lam_h[h] = float(vals[0])
    vecs_h[h] = (vals, vecs, r_h, g_h)
    lamK = wK_lam_radial(sol12, h, R_RAD)
    print("  h=%.4f: lam_min(w1) = %.6f (6 najn.: %s); waga-K (deskr.):"
          " %.6f (|delta|=%.2e)"
          % (h, lam_h[h], ["%+.4f" % x for x in vals],
             float(lamK[0]), abs(float(lamK[0]) - lam_h[h])), flush=True)
LAM_ANCHOR = lam_h[0.0125]
dev_int = abs(lam_h[0.025] - lam_h[0.0125])
ok_anchor = dev_int <= 1e-3 * abs(LAM_ANCHOR)
print("  TABELA ZBIEZNOSCI KOTWICY: h=0.05: %.6f | h=0.025: %.6f |"
      " h=0.0125: %.6f" % (lam_h[0.05], lam_h[0.025], lam_h[0.0125]))
print("  KOTWICA KANONICZNA (h=0.0125): lam_anchor = %.6f" % LAM_ANCHOR)
print("  gate wewnetrzny |lam(0.025)-lam(0.0125)| = %.2e (gate <= %.2e)"
      " -> %s" % (dev_int, 1e-3 * abs(LAM_ANCHOR),
                  "PASS" if ok_anchor else "FAIL"))
gate_all["P1b-kotwica-zbiezn"] = ok_anchor

# diagnostyka kieszeni V(r) (OBOWIAZKOWA -- anty-pulapka po f_eps)
stamp("[P1b-kan] diagnostyka kieszeni Q(r) linearyzacji (pointwise)")
V_m = 4 * BETA * g_m - 5 * GAMMA * g_m ** 2
V_vac = 4 * BETA - 5 * GAMMA
i_min = int(np.argmin(V_m))
V_min = float(V_m[i_min])
depth = V_min - V_vac
half = V_vac + 0.5 * depth
mask = V_m <= half
if mask.any():
    idx = np.where(mask)[0]
    fwhm = float(r_m[idx[-1]] - r_m[idx[0]] + H_MASTER)
    r_lo, r_hi = float(r_m[idx[0]]), float(r_m[idx[-1]])
else:
    fwhm, r_lo, r_hi = 0.0, np.nan, np.nan
print("  V(r) = 4g-5g^2 (u-forma, waga 1): V_vac = %+.3f;"
      " min = %+.6f przy r = %.4f (g=%.6f)"
      % (V_vac, V_min, float(r_m[i_min]), g_m[i_min]))
print("  glebokosc kieszeni D = V_min - V_vac = %+.6f;"
      " FWHM {V <= V_vac + D/2} = %.4f (r in [%.4f, %.4f])"
      % (depth, fwhm, r_lo, r_hi))
print("  wniosek anty-pulapkowy: FWHM/h_mostu(0.3) = %.1f"
      " (>> 1 = kieszen ROZDZIELCZA na siatce mostu; poprzednik f_eps:"
      " powloka ~0.1 przy h=0.4)" % (fwhm / 0.3))
QK_m = 2 * BETA * g_m - 3 * GAMMA * g_m ** 2 + 2 * gp_m ** 2 / g_m ** 2
iK = int(np.argmin(QK_m))
print("  deskryptywnie Q/K (forma wazona): min = %+.6f przy r=%.4f;"
      " Q/K(vac) = -1" % (float(QK_m[iK]), float(r_m[iK])))

# t*_ref: ewolucja radialna kanoniczna K_eps
stamp("[P1b-kan] t*_ref: ewolucja radialna izolowanego mu, K_eps=0.2"
      " (kontrola 0.1)")
vals_a, vecs_a, r_a, g_a = vecs_h[0.0125]
u1 = vecs_a[:, 0]
phi1 = (u1 / r_a) / g_a ** 2
phi1 /= np.sqrt(np.sum(r_a ** 2 * phi1 ** 2) * H_MASTER)
if phi1[np.argmax(np.abs(phi1))] < 0:
    phi1 = -phi1
g_eq = sol12.sol(r_a)[0]
Ke2, Kep2 = make_keps(0.2)
Ke1, Kep1 = make_keps(0.1)
dyn02 = RadialDynamicsCan(r_a, Ke2, Kep2)
dyn01 = RadialDynamicsCan(r_a, Ke1, Kep1)
runs = [("K_eps=0.2", dyn02, +0.01, 0.004),
        ("K_eps=0.2", dyn02, +0.01, 0.002),
        ("K_eps=0.2", dyn02, -0.01, 0.004),
        ("K_eps=0.1 (kontrola)", dyn01, +0.01, 0.004)]
tstars = {}
for name, dyn, amp, dt in runs:
    tb, egate = evolve_radial(dyn, g_eq + amp * phi1, dt, 20.0)
    tstars[(name, amp, dt)] = tb
    stamp("  [%s] a=%+.2f dt=%.4f: t* = %s; |dE|/E = %.2e (deskryptywnie)"
          % (name, amp, dt, ("%.3f" % tb) if tb is not None else "BRAK",
             egate))
t1 = tstars[("K_eps=0.2", +0.01, 0.004)]
t2 = tstars[("K_eps=0.2", +0.01, 0.002)]
t3 = tstars[("K_eps=0.2", -0.01, 0.004)]
ok_dt = (t1 is not None and t2 is not None
         and abs(t1 - t2) <= 0.01 * abs(t2))
print("  dt-konsystencja (LOCK <=1%%): |t*(0.004)-t*(0.002)|/t* = %s -> %s"
      % (("%.4f" % (abs(t1 - t2) / abs(t2))) if ok_dt or
         (t1 is not None and t2 is not None) else "n/d",
         "PASS" if ok_dt else "FAIL"))
ts_ok = [t for t in (t1, t2, t3) if t is not None]
TSTAR_REF = min(ts_ok) if len(ts_ok) == 3 else None
ok_tref = TSTAR_REF is not None
print("  t*_ref := min po biegach K_eps=0.2 = %s; kontrola K_eps=0.1:"
      " t* = %s" % (("%.3f" % TSTAR_REF) if ok_tref else "BRAK",
                    ("%.3f" % tstars[("K_eps=0.1 (kontrola)", 0.01, 0.004)])
                    if tstars[("K_eps=0.1 (kontrola)", 0.01, 0.004)]
                    is not None else "BRAK"))
gate_all["P1b-tref"] = ok_dt and ok_tref
print("  P1b-kan: %s" % ("PASS" if (ok_res and ok_anchor and ok_dt
                                    and ok_tref) else "FAIL"))

# ------------------------------------------------------------------ P1c-kan
stamp("[P1c-kan] most radialny->kartezjanski (operator wagi-K Phase 3):"
      " L=30, N in {76,100}")
cs_g = CubicSpline(r_m, g_m - 1.0)
cs_gp = CubicSpline(r_m, gp_m)
R_END = float(r_m[-1])


def radial_fields3d(rr):
    inside = rr <= R_END
    g3 = np.where(inside, 1.0 + cs_g(np.minimum(rr, R_END)), 1.0)
    gp3 = np.where(inside, cs_gp(np.minimum(rr, R_END)), 0.0)
    return g3, gp3


lam3 = {}
lam3_u = {}
v_mode100 = None
g3_100 = None
for N in (76, 100):
    h = L3 / N
    rr = minimg_r(N, h, L3)
    g3, gp3 = radial_fields3d(rr)
    K3 = K_fun(g3)
    Q3 = K3 * (2 * BETA * g3 - 3 * GAMMA * g3 ** 2 + 2 * gp3 ** 2 / g3 ** 2)
    A = build_op3d(g3, K_fun, Q3, h, weight3=K3)
    vals, V = lowest_eigs(A, nev=10, vecs=True)
    lam3[N] = float(vals[0])
    stamp("  N=%3d (h=%.4f): lam_min(3D, waga-K) = %.6f; 10 najn.: %s"
          % (N, h, lam3[N], ["%+.4f" % x for x in vals]))
    # deskryptywnie: u-forma 3D (waga 1)
    V3 = 4 * BETA * g3 - 5 * GAMMA * g3 ** 2
    Au = build_op3d(g3, lambda x: np.ones_like(x), V3, h)
    vals_u = lowest_eigs(Au, nev=3)
    lam3_u[N] = float(vals_u[0])
    stamp("  N=%3d: u-forma 3D (deskr.): lam_min = %.6f (|delta form| ="
          " %.2e)" % (N, lam3_u[N], abs(lam3_u[N] - lam3[N])))
    if N == 100:
        psi = V[:, 0]
        phi = (psi / np.sqrt(K3.ravel())).reshape(N, N, N)
        phi = phi / np.sqrt(np.sum(phi ** 2) * h ** 3)
        if phi.ravel()[np.argmax(np.abs(phi))] < 0:
            phi = -phi
        v_mode100 = phi
        g3_100 = g3
dev76 = abs(lam3[76] - LAM_ANCHOR)
dev100 = abs(lam3[100] - LAM_ANCHOR)
ok_l100 = dev100 <= 0.05 * abs(LAM_ANCHOR)
ok_impr = dev76 > dev100
print("  |lam(N100)-kotwica| = %.4f (gate <= %.4f): %s;"
      " trend |lam(N76)-kotwica| = %.4f > %.4f: %s"
      % (dev100, 0.05 * abs(LAM_ANCHOR), "PASS" if ok_l100 else "FAIL",
         dev76, dev100, "PASS" if ok_impr else "FAIL"))
gate_all["P1c-lam"] = ok_l100 and ok_impr

# t*_izo(3D): ewolucja kanoniczna K_eps=0.2 na siatce bramki N=100
ok_t3 = False
TSTAR_IZO3D = None
if TSTAR_REF is not None:
    stamp("[P1c-kan] t*_izo(3D): N=100, K_eps=0.2, mod wagi-K, a=+-0.01")
    h100 = L3 / 100
    dyn3 = Dyn3D(h100, Ke2, Kep2, U_fun, U1_fun)
    t_end3 = max(12.0, 2 * TSTAR_REF)
    lo, hi = 0.85 * TSTAR_REF, 1.15 * TSTAR_REF
    t3s = []
    ok_t3 = True
    for amp, dt in ((0.01, 0.02), (0.01, 0.01), (-0.01, 0.02)):
        tb, egate = evolve3d(dyn3, g3_100 + amp * v_mode100, dt, t_end3)
        okt = tb is not None and lo <= tb <= hi
        ok_t3 = ok_t3 and okt
        if tb is not None:
            t3s.append(tb)
        stamp("  evo3D a=%+.2f dt=%.3f: t* = %s (okno t*_ref+-15%%:"
              " [%.3f, %.3f]) %s; |dE|/E = %.2e (deskr.)"
              % (amp, dt, ("%.3f" % tb) if tb is not None else "BRAK",
                 lo, hi, "PASS" if okt else "FAIL", egate))
    TSTAR_IZO3D = min(t3s) if len(t3s) == 3 else None
else:
    print("  t*_izo(3D): POMINIETE (brak t*_ref)")
gate_all["P1c-tizo"] = ok_t3
p1c_pass = gate_all["P1c-lam"] and ok_t3
print("  P1c-kan: %s (t*_izo(3D) = %s -> Phase 4: t_end = 3 t*_izo = %s)"
      % ("PASS" if p1c_pass else "FAIL",
         ("%.3f" % TSTAR_IZO3D) if TSTAR_IZO3D is not None else "BRAK",
         ("%.3f" % (3 * TSTAR_IZO3D)) if TSTAR_IZO3D is not None else "n/d"))

# ------------------------------------------------------------------ werdykt
print("\n" + "=" * 78)
print("BRAMKA Phase 1 (LOCK sec. 3; P1c-kan FAIL => STOP):")
for k in ("P1a(dziedz.)", "P1d(dziedz.)", "P1b-residuum",
          "P1b-kotwica-zbiezn", "P1b-tref", "P1c-lam", "P1c-tizo"):
    print("  %s: %s" % (k, "PASS" if gate_all[k] else "FAIL"))
allpass = all(gate_all.values())
print("\nP1 GATE: %s" % ("PASS -> Phase 2" if allpass
                         else "FAIL -> STOP (raport bez Phase 2-4)"))
if allpass:
    print("lam_anchor(kanoniczna) = %.6f; t*_ref = %.3f; t*_izo(3D) = %.3f"
          " (Phase 4: t_end = %.3f)"
          % (LAM_ANCHOR, TSTAR_REF, TSTAR_IZO3D, 3 * TSTAR_IZO3D))
print("=" * 78)
