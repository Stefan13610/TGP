#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
op-bloch-chain-stability (Phase 1) -- bramka maszynerii (osiagalne FAIL-e).

LOCK: Phase0_balance.md sec. 3 Phase 1; decyzje FROZEN:
Phase_method_decisions.md sec. 1-5.

P1a (dyspersja prozni, exact): Bloch na prozni g=1 w tej samej maszynerii
    co Phase 3. Sektor tachionowy (kanoniczny: K=g^4, Q/K = 2g-3g^2+2g'^2/g^2)
    MUSI dac omega^2(k) = k^2 - 1; sektor stabilny (operator poprzednika:
    w=psi^4, potencjal -w(2psi-3psi^2-2psi'^2/psi^2)) MUSI dac k^2 + 1.
    Blad wzgledny := |num-ex| / max(|ex|,1); gate <= 1e-3 na siatce bazowej
    (N=400, komorka d=2pi, 16 k w [0,pi/d], 3 najnizsze galezie)
    ORAZ ratio maxerr(N400)/maxerr(N800) in [3,5] (rzad 2, ~x4).

P1b (gate energii ewolucji): RK4 na semi-dyskretnym ukladzie hamiltonowskim
    akcji kanonicznej z K_eps (eps=0.2 i kontrola 0.1); prozna + mala
    perturbacja; |dE|/E <= 1e-6; zbieznosc dt (x2):
    ||g_dt(t_end)-g_dt/2(t_end)||_inf <= 1e-6.

Gate Phase 1 FAIL => STOP calego cyklu (kod niewazny).
"""
import numpy as np
from scipy.linalg import eigh

BETA = 1.0
GAMMA = 1.0
D_CELL = 2 * np.pi
N_LIST = (400, 800)
NK = 16
NBRANCH = 3
GATE_DISP = 1e-3
RATIO_LO, RATIO_HI = 3.0, 5.0
EPS_LIST = (0.2, 0.1)
DT_BASE = 0.004
T_END = 4.0
GATE_E = 1e-6
GATE_DT = 1e-6

# ---------------------------------------------------------------- operators


def K_fun(g):
    return g ** 4


def Keps(g, eps):
    K = K_fun(g)
    return 0.5 * (K + np.sqrt(K ** 2 + eps ** 2))


def Keps_prime(g, eps):
    K = K_fun(g)
    Kp = 4 * g ** 3
    return 0.5 * Kp * (1 + K / np.sqrt(K ** 2 + eps ** 2))


def U_fun(g):
    return BETA * g ** 7 / 7 - GAMMA * g ** 8 / 8


def U1_fun(g):
    return BETA * g ** 6 - GAMMA * g ** 7


def bloch_lowest(g, h, d, k, sector, nev):
    """Hermitowski problem uogolniony A(k) phi = omega^2 B phi,
    symetryzowany B^-1/2 A B^-1/2.  sector: 'tach' (kanoniczny K=g^4)
    lub 'stab' (operator poprzednika, w=psi^4).  Zwraca nev najnizszych
    omega^2 oraz wektory (kolumny, w zmiennej phi)."""
    N = len(g)
    gp = (np.roll(g, -1) - np.roll(g, 1)) / (2 * h)
    if sector == 'tach':
        w = K_fun(g)
        pot = w * (2 * BETA * g - 3 * GAMMA * g ** 2 + 2 * gp ** 2 / g ** 2)
    elif sector == 'stab':
        w = g ** 4
        pot = -w * (2 * BETA * g - 3 * GAMMA * g ** 2 - 2 * gp ** 2 / g ** 2)
    else:
        raise ValueError(sector)
    wmid = 0.5 * (w + np.roll(w, -1))          # K_{j+1/2}, indeks j
    diag = (wmid + np.roll(wmid, 1)) / h ** 2 + pot
    A = np.zeros((N, N), dtype=complex)
    idx = np.arange(N)
    A[idx, idx] = diag
    off = -wmid / h ** 2                        # sprzeglo (j, j+1)
    A[idx[:-1], idx[:-1] + 1] = off[:-1]
    A[idx[:-1] + 1, idx[:-1]] = off[:-1]
    ph = np.exp(1j * k * d)
    A[N - 1, 0] = off[N - 1] * ph
    A[0, N - 1] = off[N - 1] * np.conj(ph)
    sw = np.sqrt(w)
    H = A / sw[:, None] / sw[None, :]
    vals, vecs = eigh(H, subset_by_index=(0, nev - 1))
    phi = vecs / sw[:, None]
    return vals, phi


def exact_folded(k, d, m2, nb):
    kap = k + 2 * np.pi * np.arange(-6, 7) / d
    om = np.sort(kap ** 2 + m2)
    return om[:nb]


# --------------------------------------------------------------------- P1a
print("=" * 78)
print("Phase 1 -- bramka maszynerii (LOCK sec. 3 Phase 1; decyzje FROZEN")
print("Phase_method_decisions.md).  Komorka d=2pi; metryka bledu:")
print("|num-ex|/max(|ex|,1); gate <= 1e-3 (N=400) + ratio x4 in [3,5].")
print("=" * 78)

p1a_pass = True
kgrid = np.linspace(0.0, np.pi / D_CELL, NK)
for sector, m2 in (('tach', -1.0), ('stab', +1.0)):
    maxerr = {}
    for N in N_LIST:
        h = D_CELL / N
        g = np.ones(N)
        errs = []
        for k in kgrid:
            vals, _ = bloch_lowest(g, h, D_CELL, k, sector, NBRANCH)
            ex = exact_folded(k, D_CELL, m2, NBRANCH)
            errs.append(np.abs(vals - ex) / np.maximum(np.abs(ex), 1.0))
        maxerr[N] = float(np.max(errs))
    ratio = maxerr[N_LIST[0]] / maxerr[N_LIST[1]]
    ok_base = maxerr[N_LIST[0]] <= GATE_DISP
    ok_ratio = RATIO_LO <= ratio <= RATIO_HI
    print("\n[P1a %s] m^2 = %+g: maxerr(N=400) = %.3e (gate <= 1e-3: %s);"
          % (sector, m2, maxerr[N_LIST[0]], "PASS" if ok_base else "FAIL"))
    print("  maxerr(N=800) = %.3e; ratio = %.3f (gate [3,5]: %s)"
          % (maxerr[N_LIST[1]], ratio, "PASS" if ok_ratio else "FAIL"))
    # probki dyspersji (dokumentacja): k = 0 i k = pi/d, galaz najnizsza
    h = D_CELL / N_LIST[0]
    g = np.ones(N_LIST[0])
    for k in (0.0, float(np.pi / D_CELL)):
        vals, _ = bloch_lowest(g, h, D_CELL, k, sector, NBRANCH)
        ex = exact_folded(k, D_CELL, m2, NBRANCH)
        print("  k=%.4f: num=%s vs ex=%s"
              % (k, ["%+.6f" % v for v in vals], ["%+.6f" % v for v in ex]))
    if not (ok_base and ok_ratio):
        p1a_pass = False

# --------------------------------------------------------------------- P1b


class CanonDyn:
    """Semi-dyskretny hamiltonian akcji kanonicznej z K_eps, periodyczny.
    H = sum_j h [pi_j^2/(2 K_eps(g_j)) + U(g_j)]
      + 1/2 sum_j h K_eps(g_mid,j) ((g_{j+1}-g_j)/h)^2,  g_mid=(g_j+g_j+1)/2.
    dH/dt = 0 dokladnie dla ODE -- gate mierzy czysty blad RK4 (jak #63)."""

    def __init__(self, N, h, eps):
        self.N, self.h, self.eps = N, h, eps

    def energy(self, g, pi):
        h, eps = self.h, self.eps
        gm = 0.5 * (g + np.roll(g, -1))
        dg = (np.roll(g, -1) - g) / h
        return float(np.sum(h * (pi ** 2 / (2 * Keps(g, eps)) + U_fun(g))
                            + 0.5 * h * Keps(gm, eps) * dg ** 2))

    def rhs(self, g, pi):
        h, eps = self.h, self.eps
        gm = 0.5 * (g + np.roll(g, -1))
        dg = (np.roll(g, -1) - g) / h
        Kg = Keps(g, eps)
        gdot = pi / Kg
        dH = h * (-pi ** 2 * Keps_prime(g, eps) / (2 * Kg ** 2) + U1_fun(g))
        t_flux = Keps(gm, eps) * dg                     # krawedz j
        t_quad = 0.25 * h * Keps_prime(gm, eps) * dg ** 2
        dH += -t_flux + t_quad                          # udzial krawedzi j
        dH += np.roll(t_flux + t_quad, 1)               # udzial krawedzi j-1
        pidot = -dH / h
        return gdot, pidot

    def rk4(self, g, pi, dt):
        k1g, k1p = self.rhs(g, pi)
        k2g, k2p = self.rhs(g + 0.5 * dt * k1g, pi + 0.5 * dt * k1p)
        k3g, k3p = self.rhs(g + 0.5 * dt * k2g, pi + 0.5 * dt * k2p)
        k4g, k4p = self.rhs(g + dt * k3g, pi + dt * k3p)
        return (g + dt / 6 * (k1g + 2 * k2g + 2 * k3g + k4g),
                pi + dt / 6 * (k1p + 2 * k2p + 2 * k3p + k4p))


def evolve(dyn, g0, dt, t_end):
    g = g0.copy()
    pi = np.zeros_like(g)
    E0 = dyn.energy(g, pi)
    nstep = int(round(t_end / dt))
    emax = 0.0
    for _ in range(nstep):
        g, pi = dyn.rk4(g, pi, dt)
        if not np.all(np.isfinite(g)) or float(np.min(g)) <= 0:
            return g, np.inf, E0, False
        emax = max(emax, abs(dyn.energy(g, pi) - E0) / abs(E0))
    return g, emax, E0, True


print("\n" + "-" * 78)
print("[P1b] ewolucja nieliniowa na prozni + perturbacja a=1e-3;")
print("  g0 = 1 + a(0.5 + cos x), d=2pi, N=400, t_end=4.0, RK4;")
print("  gate |dE|/E <= 1e-6; zbieznosc dt: ||g_dt - g_dt/2||_inf <= 1e-6.")
p1b_pass = True
N = 400
h = D_CELL / N
x = np.arange(N) * h
g0 = 1.0 + 1e-3 * (0.5 + np.cos(x))
for eps in EPS_LIST:
    dyn = CanonDyn(N, h, eps)
    finals = {}
    for dt in (DT_BASE, DT_BASE / 2):
        gf, emax, E0, ok = evolve(dyn, g0, dt, T_END)
        finals[dt] = gf
        gate_ok = ok and emax <= GATE_E
        print("  eps=%.1f dt=%.4f: E0=%.6f, max|dE|/E = %.2e -> %s%s"
              % (eps, dt, E0, emax, "PASS" if gate_ok else "FAIL",
                 "" if ok else " (BREAKDOWN)"))
        if not gate_ok:
            p1b_pass = False
    ddt = float(np.max(np.abs(finals[DT_BASE] - finals[DT_BASE / 2])))
    ok_dt = ddt <= GATE_DT
    print("  eps=%.1f zbieznosc dt: ||g_dt - g_dt/2||_inf = %.2e -> %s"
          % (eps, ddt, "PASS" if ok_dt else "FAIL"))
    print("  eps=%.1f wzrost tachionowy (kontrola sensownosci, nie gate): "
          "max|g-1| t=0: %.2e -> t=4: %.2e"
          % (eps, float(np.max(np.abs(g0 - 1))),
             float(np.max(np.abs(finals[DT_BASE] - 1)))))
    if not ok_dt:
        p1b_pass = False

# ------------------------------------------------------------------ VERDICT
print("\n" + "=" * 78)
print("GATE PHASE 1: P1a %s, P1b %s -> %s"
      % ("PASS" if p1a_pass else "FAIL",
         "PASS" if p1b_pass else "FAIL",
         "PASS (maszyneria wazna)" if (p1a_pass and p1b_pass)
         else "FAIL -> STOP CYKLU (kod niewazny)"))
print("=" * 78)
