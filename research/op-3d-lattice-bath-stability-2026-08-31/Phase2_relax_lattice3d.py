#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
op-3d-lattice-bath-stability (Phase 2) -- istnienie samouzgodnionego tla
sieci sc 3D (akcja kanoniczna K=g^4, U=g^7/7-g^8/8, beta=gamma=1).

LOCK: Phase0_balance.md sec. 3 Phase 2; decyzje FROZEN:
Phase_method_decisions.md sec. 0, 1, 3, 4.

EL statyka: lap(g) + (2/g)|grad g|^2 = g^2(1-g); komorka [0,d)^3
periodyczna, jeden soliton na komorke (siec sc).
Relaksacja: scipy.optimize.newton_krylov (lgmres, precond spilu(Lap-I)),
f_tol=1e-9; weryfikacja kryterium LOCKa ||R||_inf <= 1e-8.
Punkty (zalockowane): d in {pi, 2pi, 3.0790 [INPUT d*1 mu-mu], 3pi, 4pi}
x 2 starty (superpozycja periodyzowana profili izolowanych mu [INPUT
g0=2.02117]; PRIMARY 1.0x oraz 0.7x amplitudy) x N in {32,48}.
Kryteria: ||R||_inf <= 1e-8; niestalost ||g-1||_inf >= 0.05;
zbieznosc siatkowa ||g_32-g_48||_inf <= 5e-3 (wspolna podsiatka 16^3).
GUARD domenowy (implementacyjny, raportowany gdy aktywny): start obcinany
od dolu do 0.05 (dziedzina g>0 czlonu 2/g); kryteria nietkniete.

Brak istnienia dla WSZYSTKICH d => CLOSED-GATE-STOP (pelnoprawny negatyw).
"""
import time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.optimize import newton_krylov
from scipy.optimize import NoConvergence
from scipy.interpolate import CubicSpline
from scipy.integrate import solve_ivp

BETA = 1.0
GAMMA = 1.0
PHI_GOLD = (1 + 5 ** 0.5) / 2
G0_MU = PHI_GOLD * 1.24915           # INPUT (#63)
ALPHA = 2.0
EPS = 0.2
DSTAR1 = 3.0790                      # INPUT (bath-two-sectors Phase 1)
D_LIST = (np.pi, 2 * np.pi, DSTAR1, 3 * np.pi, 4 * np.pi)
D_NAMES = ("pi", "2pi", "dstar1", "3pi", "4pi")
STARTS = (1.0, 0.7)
N_LIST = (32, 48)
TOL_RES = 1e-8                       # LOCK
TOL_GRID = 5e-3                      # LOCK
TOL_AMP = 0.05                       # LOCK
BASE = ("C:/Users/Mateusz/Documents/ObsydnianMain/TGP/TGP_v1/research/"
        "op-3d-lattice-bath-stability-2026-08-31/")
OUT_NPZ = BASE + "Phase2_backgrounds3d.npz"

t0_wall = time.time()


def stamp(msg):
    print("[t=%7.1fs] %s" % (time.time() - t0_wall, msg), flush=True)


# ------------------------------------------------ profil radialny #63 (start)
f_log = lambda g: 1 + 2 * ALPHA * np.log(g)
fp_log = lambda g: 2 * ALPHA / g
W1_fun = lambda g: g ** 2 * (1 - g)


def make_feps(eps):
    def s(g):
        fv = f_log(g)
        return np.sqrt(fv ** 2 + eps ** 2)
    F = lambda g: 0.5 * (f_log(g) + s(g))
    Fp = lambda g: 0.5 * fp_log(g) * (1 + f_log(g) / s(g))
    return F, Fp


F_EPS, FP_EPS = make_feps(EPS)


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
    return r_g, sol.sol(r_g)[0], len(r_g) == N


stamp("profil radialny mu (INPUT g0=%.5f, eps=0.2)..." % G0_MU)
r_rad, g_rad, full = profile_soft(G0_MU, 60.0, 4000)
assert full, "profil niekompletny"
cs_gm1 = CubicSpline(r_rad, g_rad - 1.0)
R_END = float(r_rad[-1])


def gm1_of_r(rr):
    return np.where(rr <= R_END, cs_gm1(np.minimum(rr, R_END)), 0.0)


def start_field(d, N, scale):
    """Superpozycja periodyzowana 27 obrazow (FROZEN sec. 3)."""
    h = d / N
    x = np.arange(N) * h
    c = d / 2
    g = np.ones((N, N, N))
    for mx in (-1, 0, 1):
        for my in (-1, 0, 1):
            for mz in (-1, 0, 1):
                DX = x - c - mx * d
                DY = x - c - my * d
                DZ = x - c - mz * d
                RR = np.sqrt(DX[:, None, None] ** 2
                             + DY[None, :, None] ** 2
                             + DZ[None, None, :] ** 2)
                g += scale * gm1_of_r(RR)
    return g


# ------------------------------------------------------ residuum i precond
def residual_maker(N, h):
    def residual(gv):
        g = gv.reshape(N, N, N)
        lap = np.zeros_like(g)
        grad2 = np.zeros_like(g)
        for ax in range(3):
            gp_ = np.roll(g, -1, axis=ax)
            gm_ = np.roll(g, 1, axis=ax)
            lap += (gp_ - 2 * g + gm_) / h ** 2
            grad2 += ((gp_ - gm_) / (2 * h)) ** 2
        return (lap + 2 * grad2 / g - g ** 2 * (BETA - GAMMA * g)).ravel()
    return residual


def make_precond(N, h):
    n = N ** 3
    idx3 = np.arange(n).reshape(N, N, N)
    rows, cols, vals = [], [], []
    diag = np.full(n, -6.0 / h ** 2 - 1.0)
    for ax in range(3):
        r_ = idx3.ravel()
        c_ = np.roll(idx3, -1, axis=ax).ravel()
        v = np.full(n, 1.0 / h ** 2)
        rows += [r_, c_]
        cols += [c_, r_]
        vals += [v, v]
    rows.append(np.arange(n))
    cols.append(np.arange(n))
    vals.append(diag)
    A = sp.coo_matrix((np.concatenate(vals),
                       (np.concatenate(rows), np.concatenate(cols))),
                      shape=(n, n)).tocsc()
    ilu = spla.spilu(A, drop_tol=1e-5, fill_factor=15)
    return spla.LinearOperator((n, n), matvec=ilu.solve)


def reflect(g, ax):
    return np.roll(np.flip(g, axis=ax), 1, axis=ax)


def asymmetry(g):
    return max(float(np.max(np.abs(g - reflect(g, ax)))) for ax in range(3))


def fft_prolong(g32, N2):
    """Interpolacja trygonometryczna N1^3 -> N2^3 (pole periodyczne)."""
    N1 = g32.shape[0]
    F1 = np.fft.fftn(g32) / N1 ** 3
    freqs = (np.fft.fftfreq(N1) * N1).astype(int)
    pos = freqs % N2
    F2 = np.zeros((N2, N2, N2), dtype=complex)
    F2[np.ix_(pos, pos, pos)] = F1
    return np.real(np.fft.ifftn(F2) * N2 ** 3)


def solve_point(d, N, g_init):
    """newton_krylov (FROZEN); zwraca (g, ||R||_inf, status)."""
    h = d / N
    res = residual_maker(N, h)
    M = make_precond(N, h)
    try:
        g = newton_krylov(res, g_init.ravel(), method='lgmres',
                          inner_M=M, f_tol=1e-9, maxiter=300,
                          inner_maxiter=200, verbose=False)
        status = "OK"
    except NoConvergence as e:
        g = np.asarray(e.args[0], dtype=float)
        status = "NO-CONV"
    except Exception as e:
        return None, np.inf, "SOLVER-FAIL(%s)" % type(e).__name__
    g = g.reshape(N, N, N)
    if not np.all(np.isfinite(g)):
        return None, np.inf, "NON-FINITE"
    nrm = float(np.max(np.abs(res(g.ravel()))))
    return g, nrm, status


# ===================================================================== RUN
print("=" * 78)
print("Phase 2 -- istnienie tla sieci sc 3D (akcja kanoniczna).")
print("LOCK sec. 3 Phase 2; newton_krylov FROZEN (method_decisions sec. 4).")
print("Kryteria: ||R||_inf <= 1e-8, ||g-1||_inf >= 0.05,")
print("||g_32-g_48||_inf <= 5e-3 (wspolna podsiatka 16^3).")
print("INPUT: g0_mu=%.5f, d*1=%.4f, beta=gamma=1." % (G0_MU, DSTAR1))
print("=" * 78)

results = {}
store = {}
for dn, d in zip(D_NAMES, D_LIST):
    print("\n[d = %s = %.6f]" % (dn, d), flush=True)
    exists = False
    sols = {}
    for A_sc in STARTS:
        tag = "A%.1f" % A_sc
        out = {}
        for N in N_LIST:
            h = d / N
            if N == 32 or out.get(32) is None or out[32]["g"] is None \
                    or out[32]["klasa"] == "KOLAPS-DO-PROZNI":
                gi = start_field(d, N, A_sc)
                gmin0 = float(np.min(gi))
                clipped = gmin0 < 0.05
                if clipped:
                    gi = np.maximum(gi, 0.05)
                    stamp("  %s N=%d: GUARD start: min=%.3f < 0.05 -> clip"
                          % (tag, N, gmin0))
            else:
                gi = fft_prolong(out[32]["g"], N)
                gi = np.maximum(gi, 0.02)
            t1 = time.time()
            g, nrm, status = solve_point(d, N, gi)
            if g is None:
                out[N] = dict(g=None, res=nrm, klasa="NIEZBIEZNE",
                              status=status)
                stamp("  %s N=%d: %s (||R||=%.2e) [%.1fs]"
                      % (tag, N, status, nrm, time.time() - t1))
                continue
            amp = float(np.max(np.abs(g - 1)))
            gmin, gmax = float(np.min(g)), float(np.max(g))
            if amp < TOL_AMP:
                klasa = "KOLAPS-DO-PROZNI"
            elif gmin < 0.01:
                klasa = "UCIECZKA-g->0"
            elif nrm > TOL_RES:
                klasa = "NIEZBIEZNE"
            else:
                klasa = "OK-KANDYDAT"
            out[N] = dict(g=g, res=nrm, klasa=klasa, status=status,
                          amp=amp, gmin=gmin, gmax=gmax)
            stamp("  %s N=%2d: %s ||R||_inf=%.2e (%s), g in [%.4f,%.4f],"
                  " ||g-1||_inf=%.4f, asym=%.1e [%.1fs]"
                  % (tag, N, klasa, nrm, status, gmin, gmax, amp,
                     asymmetry(g), time.time() - t1))
        # zbieznosc siatkowa na wspolnej podsiatce 16^3
        o32, o48 = out.get(32), out.get(48)
        dgrid = None
        if o32 and o48 and o32["g"] is not None and o48["g"] is not None:
            dgrid = float(np.max(np.abs(
                o32["g"][::2, ::2, ::2] - o48["g"][::3, ::3, ::3])))
        conv_grid = dgrid is not None and dgrid <= TOL_GRID
        both_ok = (o32 and o48 and o32["klasa"] == "OK-KANDYDAT"
                   and o48["klasa"] == "OK-KANDYDAT")
        if both_ok and conv_grid:
            klasa_pt = "NIESTALE-ZBIEZNE"
            exists = True
        elif o48 is not None and o48.get("klasa") == "KOLAPS-DO-PROZNI":
            klasa_pt = "KOLAPS-DO-PROZNI"
        elif both_ok:
            klasa_pt = "NIEZBIEZNE-SIATKOWO"
        else:
            klasa_pt = "/".join(o.get("klasa", "?") if o else "?"
                                for o in (o32, o48))
        print("  %s: dgrid(16^3)=%s (gate<=5e-3) -> %s"
              % (tag, ("%.2e" % dgrid) if dgrid is not None else "n/d",
                 klasa_pt), flush=True)
        sols[tag] = (klasa_pt, out, dgrid)
    # dedup istniejacych tel
    saved = []
    for tag, (klasa_pt, out, dgrid) in sols.items():
        if klasa_pt != "NIESTALE-ZBIEZNE":
            continue
        g48 = out[48]["g"]
        dup = None
        for ptag, pg in saved:
            if float(np.max(np.abs(g48 - pg))) < 1e-8:
                dup = ptag
                break
        if dup:
            print("  %s: tlo identyczne z %s (||delta||<1e-8) -> dedup"
                  % (tag, dup))
            continue
        saved.append((tag, g48))
        store["%s__%s__N32" % (dn, tag)] = out[32]["g"]
        store["%s__%s__N48" % (dn, tag)] = out[48]["g"]
        store["%s__%s__d" % (dn, tag)] = np.array([d])
    results[dn] = dict(exists=exists, d=d,
                       klasy={t: s[0] for t, s in sols.items()})

np.savez(OUT_NPZ, **store)
print("\nzapisano tla: %s (%d tablic)" % (OUT_NPZ, len(store)))

print("\n" + "=" * 78)
print("PODSUMOWANIE Phase 2 (kryterium istnienia LOCK, doslownie):")
any_exists = False
for dn in D_NAMES:
    rr = results[dn]
    print("  d=%-6s (%.4f): istnienie tla: %s   [%s]"
          % (dn, rr["d"], "TAK" if rr["exists"] else "NIE",
             "; ".join("%s:%s" % (t, k) for t, k in rr["klasy"].items())))
    any_exists = any_exists or rr["exists"]
if any_exists:
    print("\nISTNIENIE POTWIERDZONE dla co najmniej jednego d -> Phase 3.")
else:
    print("\nBRAK istnienia dla WSZYSTKICH d -> CLOSED-GATE-STOP")
    print("(pelnoprawny negatyw: siec sc 3D nie istnieje w klasie zbadanej;")
    print(" drzewo LOCKa sec. 6 -> NEEDS: bcc/fcc / inna geometria).")
print("=" * 78)
