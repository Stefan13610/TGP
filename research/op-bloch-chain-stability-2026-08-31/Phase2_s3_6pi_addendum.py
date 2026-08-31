#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
op-bloch-chain-stability (Phase 2, ADDENDUM -- dozwolony DODATEK startu).

Powod (incydent udokumentowany, zero zmian kryteriow): dla d=6pi zalockowany
start gleboki zbiegl do tla DWUGARBNEGO (2 minima w komorce = powielona
orbita 3pi), a model M-C LOCKa specyfikuje 'jeden obiekt na komorke'.
LOCK sec. 2: 'wolno dodac starty, NIE usuwac' -- dodajemy start S3-single
(strzelanie na polokres d/2 = 3pi z pierwszej calki), zeby uzyskac tlo
JEDNOGARBNE dla d=6pi. Wynik zalockowanego startu (2-garbny) pozostaje
w Phase2_output.txt i w npz; nic nie usuwamy.

Kryteria identyczne jak Phase 2 (LOCK): ||R||_inf <= 1e-10 (pelna komorka),
||g_N - g_2N||_inf <= 1e-4, ||g-1||_inf >= 0.05.
"""
import numpy as np
from scipy.integrate import solve_ivp

import importlib.util as ilu

BASE = ("C:/Users/Mateusz/Documents/ObsydnianMain/TGP/TGP_v1/research/"
        "op-bloch-chain-stability-2026-08-31/")
spec = ilu.spec_from_file_location("ph2", BASE + "Phase2_relax_chain.py")
# UWAGA: import modulu wykonalby caly skan Phase 2 ponownie; zamiast tego
# duplikujemy tu tylko potrzebne funkcje (identyczne definicje).

BETA = 1.0
GAMMA = 1.0
TOL_RES = 1e-10
TOL_GRID = 1e-4
TOL_AMP = 0.05
D = 6 * np.pi
OUT_NPZ = BASE + "Phase2_backgrounds.npz"


def U_fun(g):
    return BETA * g ** 7 / 7 - GAMMA * g ** 8 / 8


def half_residual(g, h):
    gp = np.empty_like(g)
    gm = np.empty_like(g)
    gp[:-1] = g[1:]
    gp[-1] = g[-2]
    gm[1:] = g[:-1]
    gm[0] = g[1]
    d2 = (gp - 2 * g + gm) / h ** 2
    d1 = (gp - gm) / (2 * h)
    return d2 + 2 * d1 ** 2 / g - g ** 2 * (1 - g), d1


def half_newton(g_init, h):
    from scipy.linalg import solve_banded
    g = g_init.copy()
    M = len(g)
    best = np.inf
    stall = 0
    for it in range(200):
        R, d1 = half_residual(g, h)
        nrm = float(np.max(np.abs(R)))
        if nrm < best - 1e-16:
            best, stall = nrm, 0
        else:
            stall += 1
        if nrm < 5e-12:
            return g, nrm, "OK"
        if stall > 8:
            return g, nrm, "STAGNACJA"
        diag = -2 / h ** 2 - 2 * d1 ** 2 / g ** 2 - 2 * BETA * g \
            + 3 * GAMMA * g ** 2
        up = 1 / h ** 2 + 2 * d1 / (g * h)
        lo = 1 / h ** 2 - 2 * d1 / (g * h)
        up[0] = 2 / h ** 2
        lo[-1] = 2 / h ** 2
        ab = np.zeros((3, M))
        ab[0, 1:] = up[:-1]
        ab[1, :] = diag
        ab[2, :-1] = lo[1:]
        dg = solve_banded((1, 1), ab, -R)
        lam = 1.0
        moved = False
        for _ in range(50):
            trial = g + lam * dg
            if float(np.min(trial)) > 1e-6:
                Rt, _ = half_residual(trial, h)
                if float(np.max(np.abs(Rt))) < nrm:
                    g, moved = trial, True
                    break
            lam *= 0.5
        if not moved:
            return g, nrm, "LINE-SEARCH-FAIL"
    R, _ = half_residual(g, h)
    return g, float(np.max(np.abs(R))), "MAX-ITER"


def mirror_full(g_half, N):
    g = np.empty(N)
    M = N // 2
    g[:M + 1] = g_half
    g[M + 1:] = g_half[M - 1:0:-1]
    return g


def full_residual(g, h):
    gp = np.roll(g, -1)
    gm = np.roll(g, 1)
    d2 = (gp - 2 * g + gm) / h ** 2
    d1 = (gp - gm) / (2 * h)
    return d2 + 2 * d1 ** 2 / g - g ** 2 * (1 - g)


def half_period(gmax):
    """Polokres orbity mechanicznej startujacej z g(0)=gmax, g'(0)=0
    (pierwsze zdarzenie g'=0)."""
    def rhs(x, y):
        g, gp = y
        g = max(g, 1e-9)
        return [gp, g * g * (1 - g) - 2 * gp * gp / g]

    def ev(x, y):
        return y[1]
    ev.terminal = True
    ev.direction = 1.0          # g' rosnie przez zero w minimum
    sol = solve_ivp(rhs, [1e-12, 200.0], [gmax, 0.0], method='DOP853',
                    rtol=1e-12, atol=1e-14, events=ev, dense_output=True)
    if len(sol.t_events[0]) == 0:
        return None, None
    return float(sol.t_events[0][0]), sol


print("=" * 78)
print("Phase 2 ADDENDUM: dodany start S3-single dla d=6pi (jednogarbne tlo).")
print("Zalockowane starty i ich wyniki NIETKNIETE (Phase2_output.txt).")
print("=" * 78)

# bisekcja: half_period(gmax) = d/2 = 3pi
target = D / 2
lo_g, hi_g = 1.02, 8.0 / 7.0 - 1e-9
T_lo, _ = half_period(lo_g)
T_hi, _ = half_period(hi_g)
print("polokres(gmax=%.4f)=%.6f, polokres(gmax=%.8f)=%s (cel %.6f)"
      % (lo_g, T_lo, hi_g, T_hi, target))
for _ in range(80):
    mid = 0.5 * (lo_g + hi_g)
    Tm, _ = half_period(mid)
    if Tm is None:
        hi_g = mid
        continue
    if Tm < target:
        lo_g = mid
    else:
        hi_g = mid
gmax_star = 0.5 * (lo_g + hi_g)
T_star, sol = half_period(gmax_star)
print("bisekcja: gmax* = %.10f, polokres = %.10f (cel %.10f)"
      % (gmax_star, T_star, target))

results = {}
g_half_store = {}
for N in (400, 800):
    h = D / N
    M = N // 2
    x = np.arange(M + 1) * h
    if N == 400:
        # profil ze strzelania (skalowanie x na dokladny polokres)
        gi = sol.sol(np.clip(x * (T_star / target), 0, T_star))[0]
    else:
        from scipy.interpolate import CubicSpline
        x4 = np.arange(201) * (D / 400)
        cs = CubicSpline(x4, g_half_store[400], bc_type=((1, 0.), (1, 0.)))
        gi = cs(x)
    gh, nrm, status = half_newton(gi, h)
    g_half_store[N] = gh
    gf = mirror_full(gh, N)
    res_full = float(np.max(np.abs(full_residual(gf, h))))
    amp = float(np.max(np.abs(gf - 1)))
    nmin = int(np.sum((gf < np.roll(gf, 1)) & (gf < np.roll(gf, -1))))
    print("N=%d: ||R||_inf(pelna)=%.2e (%s), g in [%.6f, %.6f], "
          "||g-1||_inf=%.6f, minima=%d"
          % (N, res_full, status, float(np.min(gf)), float(np.max(gf)),
             amp, nmin))
    results[N] = dict(g_full=gf, res=res_full, amp=amp, nmin=nmin)

dgrid = float(np.max(np.abs(g_half_store[400] - g_half_store[800][::2])))
conv = (results[400]["res"] <= TOL_RES and results[800]["res"] <= TOL_RES
        and dgrid <= TOL_GRID)
noncst = results[800]["amp"] >= TOL_AMP
single = results[800]["nmin"] == 1
klasa = "NIESTALE-ZBIEZNE" if (conv and noncst) else "NIEZBIEZNE/KOLAPS"
print("||g_400-g_800||_inf(wezly wspolne) = %.2e (gate<=1e-4: %s)"
      % (dgrid, "PASS" if dgrid <= TOL_GRID else "FAIL"))
print("KLASA: %s; jednogarbne: %s" % (klasa, "TAK" if single else "NIE"))

if klasa == "NIESTALE-ZBIEZNE" and single:
    data = dict(np.load(OUT_NPZ))
    data["6pi__S3single__N400"] = results[400]["g_full"]
    data["6pi__S3single__N800"] = results[800]["g_full"]
    data["6pi__S3single__d"] = np.array([D])
    np.savez(OUT_NPZ, **data)
    print("dopisano tlo 6pi__S3single do %s (%d tablic)"
          % (OUT_NPZ, len(data)))
else:
    print("tlo jednogarbne dla d=6pi NIEUZYSKANE -> Phase 3 uzyje tylko"
          " tla z zalockowanego startu (2-garbnego); raportowac wprost.")
print("=" * 78)
