#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
op-bath-two-sectors (Phase 2, komparator zadeklarowany w
Phase_method_decisions.md Phase 2 pkt 9, pominiety w pierwszym skrypcie
-- uzupelnienie zadeklarowanego planu, zadna zmiana kryteriow):

IZOLACJA-W-KOMORCE: te same komorki co skan P2b (R_cell in {d/2, d},
h=0.015), obciety profil izolowany BEZ kapieli (c_bath=0); spektrum
(raw/clean) + ewolucja amp=0 do t_end = 3 t*_izo. Porownanie t* z
przebiegami kapieli (Phase2_output.txt) rozdziela: ucieczka od kapieli
vs artefakt obciecia tla w komorce.

Funkcje profile_soft/spectra/evolve: kopiowane verbatim z
Phase2_bath_runaway.py (z korekta 1).
"""
import json
import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.integrate import solve_ivp

CYC = "TGP/TGP_v1/research/op-bath-two-sectors-2026-08-23/"
PHI_GOLD = (1 + 5 ** 0.5) / 2
G0 = {"e": 1.24915, "mu": PHI_GOLD * 1.24915}
ALPHA = 2.0
EPS = 0.2
DT_BASE = 0.004
GATE = 1e-6
G_FLOOR = 0.01
TSTAR = {"e": 7.08, "mu": 3.62}          # z P2a (Phase2_output.txt)

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


def profile_soft(g0, R, N):
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
    sel = r_grid <= float(sol.t[-1])
    r_g = r_grid[sel]
    vals = sol.sol(r_g)
    g_arr, gp_arr = vals[0], vals[1]
    d2 = np.array([rhs(r_g[j], [g_arr[j], gp_arr[j]])[1]
                   for j in range(len(r_g))])
    return r_g, g_arr, gp_arr, d2, len(r_g) == N


def build_tridiag_gen(r, g, gp, d2):
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


def spectrum_F(r, g, gp, d2, nmodes=6):
    d, e = build_tridiag_gen(r, g, gp, d2)
    k = min(nmodes - 1, len(r) - 2)
    sF = np.sqrt(F_EPS(g))
    vw = eigh_tridiagonal(d / F_EPS(g), e / (sF[:-1] * sF[1:]),
                          select='i', select_range=(0, k),
                          eigvals_only=True)
    return np.array(vw)


class Dyn:
    def __init__(self, r):
        self.r = r
        self.h = r[1] - r[0]
        self.r2 = r ** 2
        self.rm2 = (0.5 * (r[:-1] + r[1:])) ** 2

    def energy(self, g, pi):
        gm = 0.5 * (g[:-1] + g[1:])
        dg = np.diff(g) / self.h
        return (np.sum(self.r2 * (pi ** 2 / (2 * F_EPS(g)) + W_fun(g)))
                * self.h
                + 0.5 * np.sum(self.rm2 * F_EPS(gm) * dg ** 2) * self.h)

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
        return gdot, -dH / (h * r2)

    def step(self, g, pi, dt):
        k1g, k1p = self.rhs(g, pi)
        k2g, k2p = self.rhs(g + 0.5 * dt * k1g, pi + 0.5 * dt * k1p)
        k3g, k3p = self.rhs(g + 0.5 * dt * k2g, pi + 0.5 * dt * k2p)
        k4g, k4p = self.rhs(g + dt * k3g, pi + dt * k3p)
        return (g + dt / 6 * (k1g + 2 * k2g + 2 * k3g + k4g),
                pi + dt / 6 * (k1p + 2 * k2p + 2 * k3p + k4p))


def evolve0(r, g_init, dt, t_end):
    dyn = Dyn(r)
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
            g, pi = dyn.step(g, pi, dt)
    return t_break, emax


print("=" * 78)
print("Komparator: IZOLACJA-W-KOMORCE (bez kapieli), amp=0, h=0.015.")
print("Deklaracja: Phase_method_decisions.md Phase 2 pkt 9.")
print("=" * 78)

with open(CYC + "Phase1_results.json") as f:
    P1 = json.load(f)
LAD = P1["ladders"]["MP"]
CONFIGS = [("ee", "e"), ("emu", "mu"), ("mumu", "mu")]
h = 0.015
for cname, central in CONFIGS:
    ds = LAD[cname]["dstars"]
    d_list = [ds[0], ds[1], ds[2], 0.5 * ds[0], 1.5 * ds[0]]
    tstar = TSTAR[central]
    print("\n[%s] central=%s, t*_izo=%.2f, t_end=3 t*_izo=%.2f"
          % (cname, central, tstar, 3 * tstar))
    for d in d_list:
        for R_cell in (0.5 * d, d):
            N = max(24, int(round(R_cell / h)))
            r, g, gp, d2, full = profile_soft(G0[central], R_cell, N)
            if not full:
                print("  d=%.3f R=%.2f: profil niekompletny" % (d, R_cell))
                continue
            vw = spectrum_F(r, g, gp, d2)
            tb, egate = evolve0(r, g, DT_BASE, 3 * tstar)
            print("  d=%6.3f R=%6.2f: w2_min(F,izol-kom)=%+9.4f; amp=0: %s"
                  " gate=%.1e" % (
                      d, R_cell, float(vw[0]),
                      ("BREAKDOWN t*=%.2f (t/t*_izo=%.2f)"
                       % (tb, tb / tstar)) if tb is not None
                      else "bez ucieczki do 3 t*_izo", egate))
print("\n(interpretacja: porownanie z Phase2_output.txt -> "
      "Phase_FINAL_close.md)")
