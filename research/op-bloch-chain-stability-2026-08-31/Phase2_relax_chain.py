#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
op-bloch-chain-stability (Phase 2) -- istnienie samouzgodnionego tla lancucha.

LOCK: Phase0_balance.md sec. 3 Phase 2; decyzje FROZEN:
Phase_method_decisions.md sec. 3.

M-C: g'' + (2/g)g'^2 = g^2(1-g), komorka [0,d) periodyczna, jeden obiekt.
Relaksacja: Newton tlumiony na POL-komorce [0,d/2] w klasie parzystej
(lustro na obu koncach; zabija dokladne zero translacyjne Jacobianu).
Po zbiegnieciu: rozwiniecie lustrzane i weryfikacja kryterium LOCKa
||residuum EL||_inf <= 1e-10 NA PELNEJ komorce periodycznej.

Punkty (zalockowane): d in {pi, 2pi, 3pi, 4pi, 6pi} x 2 starty
(gauss g=1-A exp(-(x-d/2)^2/(2 sigma^2)), A in {0.7, 0.3}, sigma=d/8).
Siatki N in {400, 800}; zbieznosc ||g_N - g_2N||_inf <= 1e-4 (wezly wspolne).
Kryterium istnienia: NIESTALE (||g-1||_inf >= 0.05) + zbiezne + residuum.
Start pomocniczy S3 (dozwolony dodatek, tylko gdy oba zalockowane starty
zawioda dla danego d): strzelanie z pierwszej calki E=1/2 g^4 g'^2 - U(g).

Brak istnienia dla WSZYSTKICH d => CLOSED-GATE-STOP (pelnoprawny negatyw).
"""
import numpy as np
from scipy.linalg import solve_banded
from scipy.interpolate import CubicSpline
from scipy.integrate import solve_ivp

BETA = 1.0
GAMMA = 1.0
D_LIST = (np.pi, 2 * np.pi, 3 * np.pi, 4 * np.pi, 6 * np.pi)
D_NAMES = ("pi", "2pi", "3pi", "4pi", "6pi")
STARTS = (0.7, 0.3)          # A (deep, shallow); sigma = d/8 (FROZEN)
N_LIST = (400, 800)
TOL_RES = 1e-10              # LOCK
TOL_GRID = 1e-4              # LOCK
TOL_AMP = 0.05               # LOCK
OUT_NPZ = ("C:/Users/Mateusz/Documents/ObsydnianMain/TGP/TGP_v1/research/"
           "op-bloch-chain-stability-2026-08-31/Phase2_backgrounds.npz")


def U_fun(g):
    return BETA * g ** 7 / 7 - GAMMA * g ** 8 / 8


def half_residual(g, h):
    """Residuum M-C na pol-komorce z lustrem na obu koncach."""
    gp = np.empty_like(g)
    gm = np.empty_like(g)
    gp[:-1] = g[1:]
    gp[-1] = g[-2]
    gm[1:] = g[:-1]
    gm[0] = g[1]
    d2 = (gp - 2 * g + gm) / h ** 2
    d1 = (gp - gm) / (2 * h)
    return d2 + 2 * d1 ** 2 / g - g ** 2 * (1 - g), d1


def half_newton(g_init, h, tag=""):
    """Newton tlumiony; zwraca (g, ||R||_inf, status)."""
    g = g_init.copy()
    M = len(g)
    best = np.inf
    stall = 0
    for it in range(200):
        R, d1 = half_residual(g, h)
        nrm = float(np.max(np.abs(R)))
        if nrm < best - 1e-16:
            best = nrm
            stall = 0
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
        up[0] = 2 / h ** 2       # lustro: sprzeglo do g1 podwojone
        lo[-1] = 2 / h ** 2
        ab = np.zeros((3, M))
        ab[0, 1:] = up[:-1]
        ab[1, :] = diag
        ab[2, :-1] = lo[1:]
        try:
            dg = solve_banded((1, 1), ab, -R)
        except Exception:
            return g, nrm, "SOLVER-FAIL"
        lam = 1.0
        moved = False
        for _ in range(50):
            trial = g + lam * dg
            if float(np.min(trial)) > 1e-6:
                Rt, _ = half_residual(trial, h)
                if float(np.max(np.abs(Rt))) < nrm:
                    g = trial
                    moved = True
                    break
            lam *= 0.5
        if not moved:
            return g, nrm, "LINE-SEARCH-FAIL"
    R, _ = half_residual(g, h)
    return g, float(np.max(np.abs(R))), "MAX-ITER"


def mirror_full(g_half, N):
    """Pol-komorka (N/2+1 wezlow) -> pelna komorka (N wezlow)."""
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


def emech(g, h):
    d1 = (np.roll(g, -1) - np.roll(g, 1)) / (2 * h)
    return 0.5 * g ** 4 * d1 ** 2 - U_fun(g)


def count_minima(g):
    lo = (g < np.roll(g, 1)) & (g < np.roll(g, -1))
    return int(np.sum(lo))


def relax_point(d, g_start_half_400):
    """Relaksacja na N=400, transfer na N=800; pelna diagnostyka."""
    out = {}
    g_half = {}
    for N in N_LIST:
        h = d / N
        M = N // 2
        if N == 400:
            gi = g_start_half_400
        else:
            x4 = np.arange(201) * (d / 400)
            cs = CubicSpline(x4, g_half[400],
                             bc_type=((1, 0.0), (1, 0.0)))
            gi = cs(np.arange(M + 1) * h)
        gh, nrm_half, status = half_newton(gi, h)
        g_half[N] = gh
        gf = mirror_full(gh, N)
        res_full = float(np.max(np.abs(full_residual(gf, h))))
        out[N] = dict(g_full=gf, g_half=gh, res=res_full,
                      status=status, gmin=float(np.min(gf)),
                      gmax=float(np.max(gf)),
                      amp=float(np.max(np.abs(gf - 1))))
    dgrid = float(np.max(np.abs(g_half[400] - g_half[800][::2])))
    out["dgrid"] = dgrid
    return out


def s3_start(d, M, h):
    """Start pomocniczy S3: strzelanie gmax -> polokres d/2 (pierwsza calka).
    Zwraca profil na pol-komorce (max w x=0, min w x=d/2) lub None."""
    def shoot(gmax):
        def rhs(x, y):
            g, gp = y
            g = max(g, 1e-9)
            return [gp, g * g * (1 - g) - 2 * gp * gp / g]
        sol = solve_ivp(rhs, [0, d / 2], [gmax, 0.0], method='DOP853',
                        rtol=1e-12, atol=1e-14, dense_output=True)
        return sol

    # bisekcja na g'(d/2) = 0; gmax in (1, 8/7)
    lo_g, hi_g = 1.0 + 1e-6, 8.0 / 7.0 - 1e-9
    f_lo = shoot(lo_g).y[1, -1]
    f_hi = shoot(hi_g).y[1, -1]
    if f_lo * f_hi > 0:
        return None
    for _ in range(80):
        mid = 0.5 * (lo_g + hi_g)
        fm = shoot(mid).y[1, -1]
        if fm == 0:
            break
        if fm * f_lo > 0:
            lo_g, f_lo = mid, fm
        else:
            hi_g = mid
    sol = shoot(0.5 * (lo_g + hi_g))
    x = np.arange(M + 1) * h
    return sol.sol(x)[0]


print("=" * 78)
print("Phase 2 -- istnienie samouzgodnionego tla lancucha (M-C).")
print("LOCK sec. 3 Phase 2; Newton na pol-komorce (method_decisions sec. 3).")
print("Kryteria: ||R||_inf <= 1e-10 (pelna komorka), ||g_N-g_2N||_inf <= 1e-4,")
print("niestalosc ||g-1||_inf >= 0.05.  Starty A in {0.7, 0.3}, sigma=d/8.")
print("=" * 78)

results = {}
store = {}
for dn, d in zip(D_NAMES, D_LIST):
    print("\n[d = %s = %.6f]" % (dn, d))
    exists = False
    sols = {}
    for A in STARTS:
        h4 = d / 400
        x4 = np.arange(201) * h4
        gi = 1.0 - A * np.exp(-0.5 * ((x4 - d / 2) / (d / 8)) ** 2)
        out = relax_point(d, gi)
        o4, o8 = out[400], out[800]
        conv_res = o4["res"] <= TOL_RES and o8["res"] <= TOL_RES
        conv_grid = out["dgrid"] <= TOL_GRID
        noncst = o8["amp"] >= TOL_AMP
        if o8["amp"] < TOL_AMP:
            klasa = "KOLAPS-DO-PROZNI"
        elif o8["gmin"] < 0.01:
            klasa = "UCIECZKA-g->0"
        elif not (conv_res and conv_grid):
            klasa = "NIEZBIEZNE"
        else:
            klasa = "NIESTALE-ZBIEZNE"
            exists = True
        em = emech(o8["g_full"], d / 800)
        print("  start A=%.1f: %s" % (A, klasa))
        print("    N=400: ||R||=%.2e (%s), g in [%.6f, %.6f]"
              % (o4["res"], o4["status"], o4["gmin"], o4["gmax"]))
        print("    N=800: ||R||=%.2e (%s), g in [%.6f, %.6f], "
              "||g-1||_inf=%.6f" % (o8["res"], o8["status"], o8["gmin"],
                                    o8["gmax"], o8["amp"]))
        print("    ||g_400-g_800||_inf(wezly wspolne)=%.2e (gate<=1e-4: %s);"
              " minima=%d; rozrzut E_mech=%.2e"
              % (out["dgrid"], "PASS" if conv_grid else "FAIL",
                 count_minima(o8["g_full"]),
                 float(np.max(em) - np.min(em))))
        sols[A] = (klasa, out)
    # dedup obu startow
    k7, k3 = sols[0.7][0], sols[0.3][0]
    if k7 == k3 == "NIESTALE-ZBIEZNE":
        dd = float(np.max(np.abs(sols[0.7][1][800]["g_full"]
                                 - sols[0.3][1][800]["g_full"])))
        print("  oba starty niestale: ||g(0.7)-g(0.3)||_inf = %.2e -> %s"
              % (dd, "IDENTYCZNE TLO" if dd < 1e-8 else "DWA ROZNE TLA"))
    # S3 fallback (dozwolony dodatek; tylko gdy oba zalockowane zawiodly)
    if not exists:
        print("  oba zalockowane starty bez rozwiazania niestalego ->"
              " start pomocniczy S3 (kwadratura/strzelanie)")
        h4 = d / 400
        gi = s3_start(d, 200, h4)
        if gi is None:
            print("    S3: brak orbity o okresie d (strzelanie bez zmiany"
                  " znaku) -> istnienie NIEPOTWIERDZONE dla tego d")
        else:
            out = relax_point(d, gi)
            o8 = out[800]
            conv = (out[400]["res"] <= TOL_RES and o8["res"] <= TOL_RES
                    and out["dgrid"] <= TOL_GRID)
            noncst = o8["amp"] >= TOL_AMP
            klasa = ("NIESTALE-ZBIEZNE" if (conv and noncst) else
                     ("KOLAPS-DO-PROZNI" if o8["amp"] < TOL_AMP
                      else "NIEZBIEZNE"))
            print("    S3: %s; N=800: ||R||=%.2e, g in [%.6f, %.6f], "
                  "dgrid=%.2e" % (klasa, o8["res"], o8["gmin"], o8["gmax"],
                                  out["dgrid"]))
            if klasa == "NIESTALE-ZBIEZNE":
                exists = True
                sols["S3"] = (klasa, out)
    results[dn] = dict(exists=exists, sols=sols, d=d)
    # zapis tel do npz (tylko niestale-zbiezne, z dedup)
    saved = []
    for key, (klasa, out) in sols.items():
        if klasa != "NIESTALE-ZBIEZNE":
            continue
        dup = False
        for prev in saved:
            if float(np.max(np.abs(out[800]["g_full"] - prev))) < 1e-8:
                dup = True
                break
        if dup:
            continue
        saved.append(out[800]["g_full"])
        tagk = str(key)
        store["%s__%s__N400" % (dn, tagk)] = out[400]["g_full"]
        store["%s__%s__N800" % (dn, tagk)] = out[800]["g_full"]
        store["%s__%s__d" % (dn, tagk)] = np.array([d])

np.savez(OUT_NPZ, **store)
print("\nzapisano tla: %s (%d tablic)" % (OUT_NPZ, len(store)))

print("\n" + "=" * 78)
print("PODSUMOWANIE Phase 2 (kryterium istnienia LOCK, doslownie):")
any_exists = False
for dn in D_NAMES:
    r = results[dn]
    print("  d=%s: istnienie tla niestalego-zbieznego: %s"
          % (dn, "TAK" if r["exists"] else "NIE"))
    if r["exists"]:
        any_exists = True
if any_exists:
    print("\nISTNIENIE POTWIERDZONE dla co najmniej jednego d -> Phase 3.")
else:
    print("\nBRAK istnienia dla WSZYSTKICH d -> CLOSED-GATE-STOP")
    print("(pelnoprawny negatyw istnienia; hipoteza drabiny wymaga 3D).")
print("=" * 78)
