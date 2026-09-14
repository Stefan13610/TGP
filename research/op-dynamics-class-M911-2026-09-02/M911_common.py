#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
op-dynamics-class-M911 -- wspolny rdzen (model, dyskretny gradient E,
detektory, starty). LOCK: Phase0_balance.md; MD: Phase_method_decisions.md.

Formy FROZEN (dziedziczone z op-metric-pair-M911, cytaty sek08a tamze):
  Ueff = w*V = gamma(psi^4/4 - psi^3/3); Ueff' = gamma psi^2(psi-1)
  Keff = K_geo psi^4 (odczyt B); B = w^2 K = psi^6/(4-3psi)^2 (lapse M9.1'')
  B' = 12 psi^5 (2-psi)/(4-3psi)^3

REJESTR [INPUT]: K_geo=gamma=1; seed=20260904 amp=1e-3; dip(0.45,sigma=1.5);
lat: skala per-siatka do psi_max=1.30; progi 5/6,7/6; pas 4/3-1e-6;
psi_low=1e-6; psi_stiff=0.12.
"""
import os
import numpy as np
from scipy import ndimage

GAMMA = 1.0
KGEO = 1.0
SEED = 20260904                                # INPUT (LOCK sec. 1)
NOISE_AMP = 1e-3
DIP_AMP = 0.45                                 # INPUT (LOCK: psi_min=0.55)
DIP_SIGMA = 1.5
PSI_MAX_START = 1.30
PSI_THR_DN = 5.0 / 6.0
PSI_THR_UP = 7.0 / 6.0
PSI_LIMIT = 4.0 / 3.0
PSI_BAND = PSI_LIMIT - 1e-6
PSI_LOW = 1e-6
PSI_STIFF = 0.12                               # LOCK (pre-rejestracja B)
L_GEN = 4 * np.pi
L_LAT = 2 * np.pi
N_PAIR = (32, 48)                              # LOCK (obie geometrie)
BASE = ("C:/Users/Mateusz/Documents/ObsydnianMain/TGP/TGP_v1/research/"
        "op-dynamics-class-M911-2026-09-02/")
NPZ_BG = ("C:/Users/Mateusz/Documents/ObsydnianMain/TGP/TGP_v1/research/"
          "op-3d-canonical-lattice-2026-08-31/Phase2_backgrounds3d.npz")


def banner(extra=""):
    print("REJESTR [INPUT]: K_geo=gamma=1; seed=%d amp=%g; dip(%.2f,"
          "sigma=%.1f) psi_min=0.55; lat->psi_max=%.2f per-siatka; "
          "M=1 (A); B=w^2K=psi^6/(4-3psi)^2 (B, lapse M9.1''); "
          "thr 5/6,7/6; pas=4/3-1e-6; psi_stiff=%.2f%s"
          % (SEED, NOISE_AMP, DIP_AMP, DIP_SIGMA, PSI_MAX_START,
             PSI_STIFF, extra), flush=True)


# ------------------------------------------------------- model (FROZEN)
def Ueff(p):
    return GAMMA * (p ** 4 / 4.0 - p ** 3 / 3.0)


def Ueffp(p):
    return GAMMA * p * p * (p - 1.0)


def Keff(p):
    return KGEO * p ** 4


def Keffp(p):
    return 4.0 * KGEO * p ** 3


def Bfun(p):
    return p ** 6 / (4.0 - 3.0 * p) ** 2


def Bprime(p):
    return 12.0 * p ** 5 * (2.0 - p) / (4.0 - 3.0 * p) ** 3


# ------------------------------------- geometria 3D + dyskretny gradient E
class Grid3D:
    """Dyskretna energia E_h (struktura strumieniowa dziedziczona) i jej
    DOKLADNY gradient dEdpsi = +delta E/delta psi (per objetosc)."""

    def __init__(self, N, L):
        self.N = N
        self.L = L
        self.h = L / N
        k1 = (2 - 2 * np.cos(2 * np.pi * np.fft.fftfreq(N))) / self.h ** 2
        kr = (2 - 2 * np.cos(2 * np.pi * np.fft.rfftfreq(N))) / self.h ** 2
        self.ksym = (k1[:, None, None] + k1[None, :, None]
                     + kr[None, None, :])          # symbol -nabla^2 >= 0

    def dEdpsi(self, g):
        h = self.h
        dH = Ueffp(g)
        for ax in range(3):
            gn = np.roll(g, -1, axis=ax)
            gm = 0.5 * (g + gn)
            dg = (gn - g) / h
            t_flux = Keff(gm) * dg / h
            t_quad = 0.25 * Keffp(gm) * dg ** 2
            dH += -t_flux + t_quad
            dH += np.roll(t_flux + t_quad, 1, axis=ax)
        return dH

    def energy(self, g):
        E = float(np.sum(Ueff(g))) * self.h ** 3
        for ax in range(3):
            gn = np.roll(g, -1, axis=ax)
            dg = (gn - g) / self.h
            E += 0.5 * float(np.sum(Keff(0.5 * (g + gn))
                                    * dg ** 2)) * self.h ** 3
        return E

    def r2_center(self):
        x = (np.arange(self.N) + 0.5) * self.h - self.L / 2
        return (x[:, None, None] ** 2 + x[None, :, None] ** 2
                + x[None, None, :] ** 2)


# --------------------------------------------- detektory (dziedziczone)
def label_mask(mask):
    lab, n = ndimage.label(mask)
    if n == 0:
        return 0, []
    parent = np.arange(n + 1)

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for ax in range(3):
        a = np.take(lab, 0, axis=ax).ravel()
        b = np.take(lab, -1, axis=ax).ravel()
        for x, y in set(zip(a.tolist(), b.tolist())):
            if x > 0 and y > 0:
                rx, ry = find(x), find(y)
                if rx != ry:
                    parent[max(rx, ry)] = min(rx, ry)
    counts = np.bincount(lab.ravel(), minlength=n + 1)
    agg = {}
    for i in range(1, n + 1):
        rt = find(i)
        agg[rt] = agg.get(rt, 0) + int(counts[i])
    return len(agg), sorted(agg.values())


def detect(psi):
    n_dn, s_dn = label_mask(psi < PSI_THR_DN)
    n_up, s_up = label_mask(psi > PSI_THR_UP)
    return n_dn, s_dn, n_up, s_up


# ------------------------------------------------------------- starty
def noise_field(N):
    """Pasmowy szum dziedziczony (|n_i|<=8), seed=20260904;
    te same wspolczynniki N=32/48; normalizacja z N=48 (MD sec. 4)."""
    rng = np.random.default_rng(SEED)
    C = rng.standard_normal((17, 17, 17, 2))
    Cc = C[..., 0] + 1j * C[..., 1]
    Csym = 0.5 * (Cc + np.conj(Cc[::-1, ::-1, ::-1]))

    def build(Ng):
        F = np.zeros((Ng, Ng, Ng), dtype=complex)
        for i in range(17):
            for j in range(17):
                for k in range(17):
                    F[(i - 8) % Ng, (j - 8) % Ng, (k - 8) % Ng] = \
                        Csym[i, j, k]
        return np.real(np.fft.ifftn(F)) * Ng ** 3

    f48 = build(48)
    scale = NOISE_AMP / float(np.max(np.abs(f48)))
    if N == 48:
        return f48 * scale
    return build(N) * scale


def start_gen(N):
    return 1.0 + noise_field(N)


def start_dip(N):
    g3 = Grid3D(N, L_GEN)
    return 1.0 - DIP_AMP * np.exp(-g3.r2_center() / (2 * DIP_SIGMA ** 2))


def start_lat(N, stamp=print):
    mt0 = os.path.getmtime(NPZ_BG)
    data = np.load(NPZ_BG)
    g = np.array(data["2pi__A1.0__N%d" % N])
    data.close()
    mt1 = os.path.getmtime(NPZ_BG)
    assert mt0 == mt1, "npz mtime zmieniony -- STOP"
    psi_raw = g ** 2
    s = (PSI_MAX_START - 1.0) / (float(np.max(psi_raw)) - 1.0)
    psi0 = 1.0 + s * (psi_raw - 1.0)
    stamp("  [lat N=%d] oryginal g[%.4f,%.4f] psi_raw[%.4f,%.4f] "
          "s=%.7f -> psi0[%.4f,%.4f]; mtime OK"
          % (N, float(np.min(g)), float(np.max(g)),
             float(np.min(psi_raw)), float(np.max(psi_raw)), s,
             float(np.min(psi0)), float(np.max(psi0))))
    return psi0


STARTS = {"gen": (start_gen, L_GEN), "dip": (start_dip, L_GEN),
          "lat": (lambda N: start_lat(N), L_LAT)}
