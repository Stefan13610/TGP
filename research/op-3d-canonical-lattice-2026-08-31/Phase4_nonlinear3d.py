#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
op-3d-canonical-lattice (Phase 4, warunkowa) -- nieliniowa weryfikacja:
superkomorka 2x2x2, perturbacja wzdluz modu minimalnego (NIE-translacyjnego)
z Phase 3, ewolucja kanoniczna K_eps do t_end = 3 t*_izo(3D).

LOCK: Phase0_balance.md sec. 3 Phase 4; decyzje FROZEN: Phase_method_
decisions.md sec. 5 (= MD-poprzednika sec. 7 z detektorem kanonicznym).
Warunek wejscia: tlo ze ZBIEZNYM Phase 3 (conv=true). Ruling kwantyfikatora
zapisany PRZED tym biegiem: Phase3_verdict_ruling.md.

Superkomorka: N_sc = 2*48 = 96 na wymiar (tlo N=48 powielone x2).
Mod: minimalny NIE-translacyjny z Phase 3 (N=48, argmin k po 4 punktach);
rozszerzenie phi_sc(x+m d) = phi(x) * prod_a (ph_a)^(m_a), ph_a = e^{ik_a d}
in {+1,-1}; normalizacja ||phi_sc||_inf = 1. a = 0.01 * ||tlo||,
||tlo|| := max|g0-1|. Ewolucja: K_eps, RK4, dt=0.01; kontrole: dt/2 przy
(+, eps=0.2), eps=0.1 przy (+); oba znaki przy eps=0.2;
t_end = 3 * t*_izo(3D) = 3 * 4.710 = 14.130 [INPUT z Phase1_output.txt].
UCIECZKA := (g <= 0.01 LUB max g >= 2 max(g0) LUB niefinityczne -- co krok)
LUB max|g - g0| > 1.0*||tlo|| (co probke 0.02). Gate energii <= 1e-6 do
momentu ucieczki. Raport t_esc vs 2 t*_izo i 3 t*_izo.

INPUT: g0_mu=1.4650974 (=phi*0.90548) [posrednio przez tla Phase 2],
d*1=3.0790, eps=0.2 (kontrola 0.1), beta=gamma=1, t*_izo(3D)=4.710.
"""
import json
import time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

BETA = 1.0
GAMMA = 1.0
RNG_SEED = 20260831
TSTAR_IZO3D = 4.710                  # INPUT: Phase1_output.txt (P1c-kan)
T_END = 3 * TSTAR_IZO3D
GATE_E = 1e-6
BASE = ("C:/Users/Mateusz/Documents/ObsydnianMain/TGP/TGP_v1/research/"
        "op-3d-canonical-lattice-2026-08-31/")
NPZ = BASE + "Phase2_backgrounds3d.npz"
JSON_P3 = BASE + "Phase3_results3d.json"

K_fun = lambda g: g ** 4
Kp_fun = lambda g: 4 * g ** 3
U_fun = lambda g: BETA * g ** 7 / 7 - GAMMA * g ** 8 / 8
U1_fun = lambda g: BETA * g ** 6 - GAMMA * g ** 7

t0_wall = time.time()


def stamp(msg):
    print("[t=%7.1fs] %s" % (time.time() - t0_wall, msg), flush=True)


def make_keps(eps):
    def s(g):
        kv = K_fun(g)
        return np.sqrt(kv ** 2 + eps ** 2)
    Ke = lambda g: 0.5 * (K_fun(g) + s(g))
    Kep = lambda g: 0.5 * Kp_fun(g) * (1 + K_fun(g) / s(g))
    return Ke, Kep


def kpoints(d):
    q = np.pi / d
    return {"Gamma": ((0.0, 0.0, 0.0), (1, 1, 1)),
            "X": ((q, 0.0, 0.0), (-1, 1, 1)),
            "M": ((q, q, 0.0), (-1, -1, 1)),
            "R": ((q, q, q), (-1, -1, -1))}


def build_op3d(g3, Phi, Q3, h, phases=(1, 1, 1), weight3=None):
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


_V0 = {}


def v0_det(n):
    if n not in _V0:
        v = np.random.default_rng(RNG_SEED).standard_normal(n)
        _V0[n] = v / np.linalg.norm(v)
    return _V0[n]


def lowest_eigs(A, nev=10, ncv=80, vecs=False):
    n = A.shape[0]
    kw = dict(k=nev, which='SA', tol=0, ncv=min(ncv, n - 1),
              maxiter=200000, v0=v0_det(n))
    try:
        if vecs:
            vals, V = spla.eigsh(A, **kw)
            o = np.argsort(vals)
            return vals[o], V[:, o]
        return np.sort(spla.eigsh(A, return_eigenvectors=False, **kw)), None
    except spla.ArpackNoConvergence:
        return None, None


def grad2_of(g, h):
    s = 0
    for ax in range(3):
        s = s + ((np.roll(g, -1, axis=ax) - np.roll(g, 1, axis=ax))
                 / (2 * h)) ** 2
    return s


def tach_Q(g, h):
    return K_fun(g) * (2 * BETA * g - 3 * GAMMA * g ** 2
                       + 2 * grad2_of(g, h) / g ** 2)


def min_nontrans_mode(g, h, d, kname):
    """Najnizszy NIE-translacyjny mod operatora wagi-K w punkcie kname
    (identyfikacja translacji jak Phase 3: overlap B-wazony >= 0.9)."""
    kv, phs = kpoints(d)[kname]
    K3 = K_fun(g)
    A = build_op3d(g, K_fun, tach_Q(g, h), h, phases=phs, weight3=K3)
    vals, V = lowest_eigs(A, nev=10, vecs=True)
    if vals is None:
        return None, None, None
    sqw = np.sqrt(K3.ravel())
    taus = []
    for ax in range(3):
        t = ((np.roll(g, -1, axis=ax) - np.roll(g, 1, axis=ax))
             / (2 * h)).ravel() * sqw
        taus.append(t)
    T = np.stack(taus, axis=1)
    Qt, _ = np.linalg.qr(T)
    for j in range(len(vals)):
        psi = V[:, j]
        ovl = float(np.linalg.norm(Qt.T @ psi) / np.linalg.norm(psi))
        if ovl < 0.9:
            phi = (psi / sqw).reshape(g.shape)
            return float(vals[j]), phi, phs
    return None, None, None


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


def evolve_sc(dyn, g0_sc, g_init, dt, t_end, amp_bg, gmax0):
    """Ewolucja superkomorki; zwraca (t_esc, mechanizm, egate_do_ucieczki,
    maxdev_final)."""
    g = g_init.copy()
    pi = np.zeros_like(g)
    E0 = dyn.energy(g, pi)
    nstep = int(round(t_end / dt))
    every = max(1, int(round(0.02 / dt)))
    emax = 0.0
    for k in range(nstep + 1):
        if not np.all(np.isfinite(g)):
            return k * dt, "non-finite", emax, np.inf
        gmin, gmax = float(np.min(g)), float(np.max(g))
        if gmin <= 0.01:
            return k * dt, "g<=0.01", emax, float(np.max(np.abs(g - g0_sc)))
        if gmax >= 2.0 * gmax0:
            return k * dt, "g>=2max(g0)", emax, \
                float(np.max(np.abs(g - g0_sc)))
        if k % every == 0:
            dev = float(np.max(np.abs(g - g0_sc)))
            if dev > 1.0 * amp_bg:
                return k * dt, "dev>||tlo||", emax, dev
            E = dyn.energy(g, pi)
            if np.isfinite(E):
                emax = max(emax, abs(E - E0) / abs(E0))
        if k < nstep:
            g, pi = dyn.rk4(g, pi, dt)
    return None, "brak-ucieczki", emax, float(np.max(np.abs(g - g0_sc)))


# ===================================================================== RUN
print("=" * 78)
print("Phase 4 -- nieliniowa weryfikacja: superkomorka 2x2x2, K_eps,")
print("t_end = 3 t*_izo(3D) = %.3f. LOCK sec. 3 Phase 4; MD sec. 5." % T_END)
print("INPUT: t*_izo(3D)=%.3f (Phase1), tla: Phase2_backgrounds3d.npz,"
      % TSTAR_IZO3D)
print("wybor modu: Phase3_results3d.json (tla ze zbieznym Phase 3).")
print("=" * 78)

data = np.load(NPZ)
with open(JSON_P3) as f:
    p3 = json.load(f)

Ke2, Kep2 = make_keps(0.2)
Ke1, Kep1 = make_keps(0.1)
summary = {}
for label, rr in p3["bgs"].items():
    if rr.get("w2min") is None or not rr.get("conv"):
        print("\n[tlo %s]: Phase 3 niezbiezne/INCOMPLETE -> Phase 4 POMINIETA"
              " (warunek LOCKa: zbiezny Phase 3)" % label)
        continue
    dn, tag = label.split("/")
    d = rr["d"]
    argk = rr["argk"]
    g48 = data["%s__%s__N48" % (dn, tag)]
    h = d / 48
    stamp("[tlo %s] d=%.6f, w2_min=%+.6f (argmin %s) -> mod perturbacji"
          % (label, d, rr["w2min"], argk))
    lam_j, phi, phs = min_nontrans_mode(g48, h, d, argk)
    if phi is None:
        print("  eigsh INCOMPLETE przy odtwarzaniu modu -> punkt INCOMPLETE")
        summary[label] = dict(status="INCOMPLETE-mode")
        continue
    print("  mod: lambda = %+.6f (zgodnosc z Phase 3 w2_min: |d|=%.2e)"
          % (lam_j, abs(lam_j - rr["w2min"])))
    # superkomorka 2x2x2 z fazami +-1
    sgn = {1: np.array([1.0, 1.0]), -1: np.array([1.0, -1.0])}
    g_sc = np.tile(g48, (2, 2, 2))
    phi_sc = np.tile(phi, (2, 2, 2))
    N = 48
    for ax, p in enumerate(phs):
        if p == -1:
            sl = [slice(None)] * 3
            sl[ax] = slice(N, 2 * N)
            phi_sc[tuple(sl)] *= -1.0
    phi_sc = phi_sc / float(np.max(np.abs(phi_sc)))
    amp_bg = float(np.max(np.abs(g48 - 1.0)))
    gmax0 = float(np.max(g48))
    a = 0.01 * amp_bg
    print("  ||tlo|| = %.4f, a = %.6f, N_sc = 96, dt = 0.01, t_end = %.3f"
          % (amp_bg, a, T_END))
    runs = [("+a eps=0.2 dt=0.01", +a, Ke2, Kep2, 0.01),
            ("+a eps=0.2 dt=0.005", +a, Ke2, Kep2, 0.005),
            ("-a eps=0.2 dt=0.01", -a, Ke2, Kep2, 0.01),
            ("+a eps=0.1 dt=0.01", +a, Ke1, Kep1, 0.01)]
    res_runs = {}
    for name, aa, Ke, Kep, dt in runs:
        dyn = Dyn3D(h, Ke, Kep, U_fun, U1_fun)
        t1 = time.time()
        t_esc, mech, egate, dev = evolve_sc(dyn, g_sc, g_sc + aa * phi_sc,
                                            dt, T_END, amp_bg, gmax0)
        okE = egate <= GATE_E
        res_runs[name] = dict(t_esc=t_esc, mech=mech, egate=egate)
        stamp("  [%s]: %s (t_esc=%s, mechanizm=%s); vs 2t*=%.2f / 3t*=%.2f;"
              " gate|dE|/E=%.2e %s; maxdev=%.3f [%.0fs]"
              % (name,
                 "UCIECZKA" if t_esc is not None else "BRAK UCIECZKI",
                 ("%.3f" % t_esc) if t_esc is not None else "-", mech,
                 2 * TSTAR_IZO3D, 3 * TSTAR_IZO3D, egate,
                 "PASS" if okE else "FAIL", dev, time.time() - t1))
    summary[label] = dict(status="OK", d=d, w2min=rr["w2min"],
                          runs=res_runs)

print("\n" + "=" * 78)
print("PODSUMOWANIE Phase 4 (kategorie LOCKa sec. 3 Phase 4 -- ocena w"
      " FINAL close wg rulingu):")
for label, s in summary.items():
    if s.get("status") != "OK":
        print("  %s: %s" % (label, s.get("status")))
        continue
    esc02 = [s["runs"][n]["t_esc"] for n in s["runs"] if "eps=0.2" in n]
    both_signs_no_esc = all(t is None for t in
                            (s["runs"]["+a eps=0.2 dt=0.01"]["t_esc"],
                             s["runs"]["-a eps=0.2 dt=0.01"]["t_esc"]))
    esc_fast = all(t is not None and t <= 2 * TSTAR_IZO3D for t in esc02)
    print("  %s (w2_min=%+.5f): brak ucieczki (oba znaki, eps=0.2): %s;"
          " ucieczka <= 2t* (wszystkie eps=0.2): %s"
          % (label, s["w2min"], both_signs_no_esc, esc_fast))
print("=" * 78)
