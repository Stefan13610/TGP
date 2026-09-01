#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
op-metametric-boundary (Phase 2) -- Q2: relaksacja z podloga substratowa
(RACHUNEK CENTRALNY).

LOCK: Phase0_balance.md sec. 3 Phase 2; decyzje FROZEN:
Phase_method_decisions.md sec. 1-6 (mapowanie g_floor=sqrt(Phi_c/Phi_vac),
kara C^2 V_fl=(kappa/3)(g_floor-g)^3 dla g<g_floor, kappa=100;
gradient flow L2 stabilizowany semi-implicit Euler dt=0.01, dt/2=0.005;
detektor: ndimage.label maski g<(1+g_floor)/2 co dt=1, periodyczne
sklejanie, N_seed=liczba obiektow w t=0, NUKLEACJA = N_obj>N_seed
utrzymane >=10 j.cz.).

Uzycie:
  python Phase2_floor_relax.py p2a            # gate maszynerii (STOP przy FAIL)
  python Phase2_floor_relax.py job <id>       # pojedynczy bieg P2b/P2c
  python Phase2_floor_relax.py jobs <id1> ... # sekwencja biegow
  python Phase2_floor_relax.py list           # lista biegow
  python Phase2_floor_relax.py verdict        # zlozenie werdyktu Q2

REJESTR WEJSC [INPUT, flagowane]: progi QB-2 {0.197,0.298,0.331};
seed=20260901, amp=1e-3; g0_mu=phi*0.90548=1.4650974; beta=gamma=1;
kotwica lam_min=-1.646589 (tabela h -> siatki radialne {0.025,0.0125});
tlo 2pi z ../op-3d-canonical-lattice-2026-08-31/Phase2_backgrounds3d.npz
(READ-ONLY; klucze 2pi__A1.0__N{32,48}).
"""
import json
import os
import sys
import time
import numpy as np
from scipy import ndimage
from scipy.integrate import solve_ivp
from scipy.linalg import solve_banded

BETA = 1.0
GAMMA = 1.0
PHI_GOLD = (1 + 5 ** 0.5) / 2
G0_MU = PHI_GOLD * 0.90548           # INPUT
KAPPA = 100.0                        # FROZEN (MD sec. 2)
SEED = 20260901                      # INPUT (LOCK sec. 2)
NOISE_AMP = 1e-3                     # INPUT (LOCK sec. 2)
THRESHOLDS = (0.197, 0.298, 0.331)   # INPUT QB-2
FLOORS = tuple(np.sqrt(t) for t in THRESHOLDS)
R_RAD = 60.0
H_RAD = (0.025, 0.0125)              # FROZEN (MD sec. 3)
L3 = 2 * np.pi
N3 = (32, 48)
DT_MAIN = 0.01
T_MAX = 200.0
STAT_TOL = 1e-8
BASE = ("C:/Users/Mateusz/Documents/ObsydnianMain/TGP/TGP_v1/research/"
        "op-metametric-boundary-2026-09-01/")
NPZ_BG = ("C:/Users/Mateusz/Documents/ObsydnianMain/TGP/TGP_v1/research/"
          "op-3d-canonical-lattice-2026-08-31/Phase2_backgrounds3d.npz")
RESDIR = BASE + "Phase2_results/"

t0_wall = time.time()


def stamp(msg):
    print("[t=%7.1fs] %s" % (time.time() - t0_wall, msg), flush=True)


# ------------------------------------------------------- model kanoniczny
K_fun = lambda g: g ** 4
Kp_fun = lambda g: 4 * g ** 3


def U_parts(sector):
    """sector 'tach': U = b g^7/7 - c g^8/8 (m^2(1) = -gamma);
    sector 'stab': U_stab = -b g^7/7 + c g^8/8 (m^2(1) = +gamma)."""
    s = 1.0 if sector == "tach" else -1.0
    U = lambda g: s * (BETA * g ** 7 / 7 - GAMMA * g ** 8 / 8)
    U1 = lambda g: s * (BETA * g ** 6 - GAMMA * g ** 7)
    return U, U1


def vfl_parts(gf):
    """Kara C^2 (FROZEN MD sec.2): V=(kappa/3)(gf-g)^3 dla g<gf."""
    def V(g):
        d = np.maximum(gf - g, 0.0)
        return (KAPPA / 3.0) * d ** 3
    def V1(g):
        d = np.maximum(gf - g, 0.0)
        return -KAPPA * d ** 2
    return V, V1


def make_keps(eps):
    def s(g):
        kv = K_fun(g)
        return np.sqrt(kv ** 2 + eps ** 2)
    Ke = lambda g: 0.5 * (K_fun(g) + s(g))
    Kep = lambda g: 0.5 * Kp_fun(g) * (1 + K_fun(g) / s(g))
    return Ke, Kep


def profile_canonical(g0, R, rtol=1e-12, atol=1e-14):
    def rhs(r_, y):
        gv, gp = y
        gv = max(gv, 1e-3)
        drv = gv ** 2 * (BETA - GAMMA * gv) - (2.0 / gv) * gp ** 2
        if r_ < 1e-10:
            return [gp, drv / 3.0]
        return [gp, drv - 2 * gp / r_]
    sol = solve_ivp(rhs, [1e-6, R + 0.05], [g0, 0.0], method='DOP853',
                    max_step=0.02, rtol=rtol, atol=atol, dense_output=True)
    assert float(sol.t[-1]) >= R, "profil niekompletny -- STOP"
    return sol


# ------------------------------------------------- szum pasmowy (FROZEN)
def noise_field(N):
    """Pasmowo ograniczony szum |n_i|<=8, rng(20260901); te same
    wspolczynniki na kazdej siatce; skala z max|f| na N=48."""
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


# --------------------------------------------------------------- rhs flow
class FlowRadial:
    def __init__(self, h, U1t, Ut):
        N = int(round(R_RAD / h))
        self.h = h
        self.r = (np.arange(N) + 0.5) * h
        self.r2 = self.r ** 2
        self.rm2 = (0.5 * (self.r[:-1] + self.r[1:])) ** 2
        self.U1t, self.Ut = U1t, Ut
        # macierz semi-implicit: I - dt*A*L_r (trojdiagonalna, K=1)
        self.a_face = self.rm2  # r_{j+1/2}^2

    def rhs(self, g):
        h = self.h
        gm = 0.5 * (g[:-1] + g[1:])
        dg = np.diff(g) / h
        dH = h * self.r2 * self.U1t(g)
        t_flux = self.rm2 * K_fun(gm) * dg
        t_quad = 0.25 * h * self.rm2 * Kp_fun(gm) * dg ** 2
        dH[:-1] += -t_flux + t_quad
        dH[1:] += t_flux + t_quad
        return -dH / (h * self.r2)

    def energy(self, g):
        gm = 0.5 * (g[:-1] + g[1:])
        dg = np.diff(g) / self.h
        return float(np.sum(self.r2 * self.Ut(g)) * self.h
                     + 0.5 * np.sum(self.rm2 * K_fun(gm) * dg ** 2)
                     * self.h)

    def step(self, g, dt):
        r = self.rhs(g)
        A = 1.05 * K_fun(float(np.max(g)))
        N = len(g)
        h2 = self.h ** 2
        lo = np.zeros(N)
        di = np.ones(N)
        up = np.zeros(N)
        c = dt * A / h2
        af = self.a_face
        di[:-1] += c * af / self.r2[:-1]
        di[1:] += c * af / self.r2[1:]
        up[1:] = -c * af / self.r2[:-1]
        lo[:-1] = -c * af / self.r2[1:]
        ab = np.vstack([up, di, lo])
        dgv = solve_banded((1, 1), ab, dt * r)
        return g + dgv, r


class Flow3D:
    def __init__(self, N, U1t, Ut):
        self.N = N
        self.h = L3 / N
        self.U1t, self.Ut = U1t, Ut
        k1 = (2 - 2 * np.cos(2 * np.pi * np.fft.fftfreq(N))) / self.h ** 2
        kr = (2 - 2 * np.cos(2 * np.pi * np.fft.rfftfreq(N))) / self.h ** 2
        self.ksym = (k1[:, None, None] + k1[None, :, None]
                     + kr[None, None, :])

    def rhs(self, g):
        h = self.h
        dH = self.U1t(g)
        for ax in range(3):
            gn = np.roll(g, -1, axis=ax)
            gm = 0.5 * (g + gn)
            dg = (gn - g) / h
            t_flux = K_fun(gm) * dg / h
            t_quad = 0.25 * Kp_fun(gm) * dg ** 2
            dH += -t_flux + t_quad
            dH += np.roll(t_flux + t_quad, 1, axis=ax)
        return -dH

    def energy(self, g):
        E = float(np.sum(self.Ut(g))) * self.h ** 3
        for ax in range(3):
            gn = np.roll(g, -1, axis=ax)
            dg = (gn - g) / self.h
            E += 0.5 * float(np.sum(K_fun(0.5 * (g + gn)) * dg ** 2)) \
                * self.h ** 3
        return E

    def step(self, g, dt):
        r = self.rhs(g)
        A = 1.05 * K_fun(float(np.max(g)))
        rh = np.fft.rfftn(dt * r)
        rh /= (1.0 + dt * A * self.ksym)
        return g + np.fft.irfftn(rh, s=g.shape, axes=(0, 1, 2)), r


# ------------------------------------------------------ detektor (FROZEN)
def label_objects(g, gthr):
    """Zwraca (liczba obiektow, posortowane rozmiary) maski g<gthr;
    3D periodycznie (sklejanie scian union-find), 1D bez periodycznosci."""
    mask = g < gthr
    if mask.ndim == 1:
        lab, n = ndimage.label(mask)
        if n == 0:
            return 0, []
        sizes = np.bincount(lab.ravel())[1:]
        return int(n), sorted(int(x) for x in sizes)
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


# ------------------------------------------------------------ silnik biegu
def run_flow(flow, g0, gf, dt, label, tmax=T_MAX):
    """Gradient flow do stacjonarnosci / nukleacji / t_max / zalamania.
    Zwraca dict wyniku + stan koncowy."""
    np.seterr(over='ignore', invalid='ignore')   # korekta 1: NaN
    gthr = (1.0 + gf) / 2.0
    steps_per_unit = int(round(1.0 / dt))
    nsteps = int(round(tmax / dt))
    g = g0.copy()
    n0, sizes0 = label_objects(g, gthr)
    series = [dict(t=0.0, n=n0, sizes=sizes0,
                   gmin=float(np.min(g)), gmax=float(np.max(g)),
                   gnorm=float(np.max(np.abs(flow.rhs(g)))),
                   E=flow.energy(g))]
    stamp("  [%s] start: N_seed=%d, g in [%.4f,%.4f], E=%.6f, "
          "||gdot||=%.2e" % (label, n0, series[0]["gmin"],
                             series[0]["gmax"], series[0]["E"],
                             series[0]["gnorm"]))
    status, t_end, t_nuc, n_det = "INCOMPLETE", tmax, None, None
    streak = 0
    g_prev = g.copy()
    for k in range(1, nsteps + 1):
        g_prev = g
        try:
            g, _ = flow.step(g, dt)         # korekta 1: try/except
            mx = float(np.max(g)) if np.all(np.isfinite(g)) else np.nan
        except (ValueError, FloatingPointError, OverflowError):
            mx = np.nan
        if not np.isfinite(mx):
            status, t_end = "BREAKDOWN", k * dt
            g = g_prev                       # ostatni finityczny stan
            stamp("  [%s] ZALAMANIE (niefinitycznosc) t=%.3f; ostatni"
                  " finityczny stan: g in [%.4g, %.4g]"
                  % (label, t_end, float(np.min(g)), float(np.max(g))))
            break
        if k % steps_per_unit == 0:
            t = k * dt
            r = flow.rhs(g)
            gnorm = float(np.max(np.abs(r)))
            n, sizes = label_objects(g, gthr)
            series.append(dict(t=t, n=n, sizes=sizes,
                               gmin=float(np.min(g)), gmax=mx,
                               gnorm=gnorm, E=flow.energy(g)))
            if n > n0:
                streak += 1
            else:
                streak = 0
            if streak >= 11:
                t_nuc = t - 10.0
                n_det = series[-11]["n"]
                status, t_end = "NUCLEATION", t
                stamp("  [%s] NUKLEACJA: t0=%.1f, N_obj(t0)=%d > "
                      "N_seed=%d, utrzymana 10 j.cz." % (label, t_nuc,
                                                         n_det, n0))
                break
            if gnorm <= STAT_TOL:
                status, t_end = "STATIONARY", t
                stamp("  [%s] STACJONARNOSC t=%.1f: ||gdot||=%.2e, g in "
                      "[%.4f,%.4f]" % (label, t, gnorm, series[-1]["gmin"],
                                       mx))
                break
    if status == "BREAKDOWN":
        n_b, sizes_b = label_objects(g, gthr)
        series.append(dict(t=t_end, n=n_b, sizes=sizes_b,
                           gmin=float(np.min(g)), gmax=float(np.max(g)),
                           gnorm=float('inf'), E=float('nan')))
    last = series[-1]
    dev_const = 0.5 * (last["gmax"] - last["gmin"])
    res = dict(label=label, status=status, t_end=t_end, t_nuc=t_nuc,
               n_seed=n0, n_det=n_det, gthr=gthr,
               gmin=last["gmin"], gmax=last["gmax"],
               gnorm_end=last["gnorm"], E_end=last["E"],
               dev_const=dev_const, n_end=last["n"],
               sizes_end=last["sizes"], series=series)
    stamp("  [%s] KONIEC: %s (t=%.1f), g in [%.4f,%.4f], ||g-const||=%.4f,"
          " N_obj=%d (seed %d), ||gdot||=%.2e"
          % (label, status, t_end, last["gmin"], last["gmax"], dev_const,
             last["n"], n0, last["gnorm"]))
    return res, g


# ------------------------------------------------------------- starty
_SOL_CACHE = {}


def start_radial_sol(h):
    if "sol" not in _SOL_CACHE:
        _SOL_CACHE["sol"] = profile_canonical(G0_MU, R_RAD)
    sol = _SOL_CACHE["sol"]
    N = int(round(R_RAD / h))
    r = (np.arange(N) + 0.5) * h
    return sol.sol(r)[0].copy()


def start_lat(N):
    mt0 = os.path.getmtime(NPZ_BG)
    data = np.load(NPZ_BG)
    g = np.array(data["2pi__A1.0__N%d" % N])
    data.close()
    mt1 = os.path.getmtime(NPZ_BG)
    assert mt0 == mt1, "npz mtime zmieniony -- STOP"
    return g


def start_noise(N):
    return 1.0 + noise_field(N)


# ------------------------------------------------------------- rejestr biegow
def job_registry():
    jobs = {}
    for sector in ("tach", "stab"):
        for fi, gf in enumerate(FLOORS):
            for h in H_RAD:
                jid = "%s_sol_f%d_h%s" % (sector, fi + 1,
                                          ("025" if h == 0.025 else "0125"))
                jobs[jid] = dict(sector=sector, start="sol", fi=fi, h=h)
            for N in N3:
                jobs["%s_lat_f%d_N%d" % (sector, fi + 1, N)] = \
                    dict(sector=sector, start="lat", fi=fi, N=N)
                jobs["%s_noise_f%d_N%d" % (sector, fi + 1, N)] = \
                    dict(sector=sector, start="noise", fi=fi, N=N)
    return jobs


def run_job(jid, dtfactor=1):
    os.makedirs(RESDIR, exist_ok=True)
    jobs = job_registry()
    base_id = jid[:-4] if jid.endswith("_dt2") else jid
    if jid.endswith("_dt2"):
        dtfactor = 2
    spec = jobs[base_id]
    dt = DT_MAIN / dtfactor
    U, U1 = U_parts(spec["sector"])
    gf = FLOORS[spec["fi"]]
    V, V1 = vfl_parts(gf)
    Ut = lambda g: U(g) + V(g)
    U1t = lambda g: U1(g) + V1(g)
    if spec["start"] == "sol":
        flow = FlowRadial(spec["h"], U1t, Ut)
        g0 = start_radial_sol(spec["h"])
        geom = dict(kind="radial", h=spec["h"], R=R_RAD)
    else:
        flow = Flow3D(spec["N"], U1t, Ut)
        g0 = start_lat(spec["N"]) if spec["start"] == "lat" \
            else start_noise(spec["N"])
        geom = dict(kind="3d", N=spec["N"], L=L3)
    label = jid + (" dt=%.4g" % dt)
    print("REJESTR [INPUT]: prog QB-2=%.3f -> g_floor=%.7f (kappa=%.0f "
          "FROZEN); seed=%d amp=%g; g0_mu=%.7f; beta=gamma=1"
          % (THRESHOLDS[spec["fi"]], gf, KAPPA, SEED, NOISE_AMP, G0_MU),
          flush=True)
    res, gfin = run_flow(flow, g0, gf, dt, label)
    res.update(job=jid, dt=dt, geom=geom, g_floor=gf,
               threshold=THRESHOLDS[spec["fi"]])
    # dekompresja serii do json-serializowalnej
    with open(RESDIR + jid + ".json", "w") as f:
        json.dump(res, f)
    np.savez_compressed(RESDIR + jid + ".npz", g=gfin,
                        meta=np.array([res["g_floor"], dt], dtype=float))
    stamp("zapisano: %s.json / .npz" % (RESDIR + jid))
    return res


# ---------------------------------------------------------------- P2a gate
def run_p2a():
    print("=" * 78)
    print("P2a -- GATE maszynerii podlogi (LOCK sec. 3 P2a; FROZEN MD "
          "sec. 5). FAIL => STOP.")
    print("REJESTR [INPUT]: progi QB-2 %s -> g_floor %s; kappa=%.0f;"
          % (THRESHOLDS, tuple(round(f, 7) for f in FLOORS), KAPPA))
    print("g0_mu=%.7f; beta=gamma=1; seed=%d (nieuzywany w P2a)."
          % (G0_MU, SEED))
    print("=" * 78)
    U, U1 = U_parts("tach")
    allpass = True
    # P2a-i: proznia bez zaburzen zostaje w g=1 (gradient flow z podloga)
    for fi, gf in enumerate(FLOORS):
        V, V1 = vfl_parts(gf)
        Ut = lambda g: U(g) + V(g)
        U1t = lambda g: U1(g) + V1(g)
        fr = FlowRadial(0.0125, U1t, Ut)
        g = np.ones(len(fr.r))
        dev = 0.0
        for k in range(int(round(10.0 / DT_MAIN))):
            g, _ = fr.step(g, DT_MAIN)
            dev = max(dev, float(np.max(np.abs(g - 1.0))))
        ok_r = dev <= 1e-12
        f3 = Flow3D(32, U1t, Ut)
        g3 = np.ones((32, 32, 32))
        dev3 = 0.0
        for k in range(int(round(10.0 / DT_MAIN))):
            g3, _ = f3.step(g3, DT_MAIN)
            dev3 = max(dev3, float(np.max(np.abs(g3 - 1.0))))
        ok_3 = dev3 <= 1e-12
        print("  P2a-i floor=%.7f: radial ||g-1||_inf=%.2e -> %s | "
              "3D N=32 ||g-1||_inf=%.2e -> %s"
              % (gf, dev, "PASS" if ok_r else "FAIL", dev3,
                 "PASS" if ok_3 else "FAIL"), flush=True)
        allpass = allpass and ok_r and ok_3
    # P2a-ii: hamiltonowski gate energii tam, gdzie podloga nieaktywna
    Ke, Kep = make_keps(0.2)
    for fi, gf in enumerate(FLOORS):
        V, V1 = vfl_parts(gf)
        Ut = lambda g: U(g) + V(g)
        U1t = lambda g: U1(g) + V1(g)
        # radialnie
        h = 0.0125
        N = int(round(R_RAD / h))
        r = (np.arange(N) + 0.5) * h
        r2 = r ** 2
        rm2 = (0.5 * (r[:-1] + r[1:])) ** 2

        def energy_r(g, pi):
            gm = 0.5 * (g[:-1] + g[1:])
            dg = np.diff(g) / h
            return (np.sum(r2 * (pi ** 2 / (2 * Ke(g)) + Ut(g))) * h
                    + 0.5 * np.sum(rm2 * Ke(gm) * dg ** 2) * h)

        def rhs_r(g, pi):
            gm = 0.5 * (g[:-1] + g[1:])
            dg = np.diff(g) / h
            Kg = Ke(g)
            dH = h * r2 * (-pi ** 2 * Kep(g) / (2 * Kg ** 2) + U1t(g))
            t_flux = rm2 * Ke(gm) * dg
            t_quad = 0.25 * h * rm2 * Kep(gm) * dg ** 2
            dH[:-1] += -t_flux + t_quad
            dH[1:] += t_flux + t_quad
            return pi / Kg, -dH / (h * r2)

        g = 1.0 + 1e-3 * np.exp(-(r - 10.0) ** 2)
        pi = np.zeros_like(g)
        E0 = energy_r(g, pi)
        emax, gmin_run = 0.0, np.inf
        dt = 0.004
        for k in range(int(round(4.0 / dt))):
            k1g, k1p = rhs_r(g, pi)
            k2g, k2p = rhs_r(g + 0.5 * dt * k1g, pi + 0.5 * dt * k1p)
            k3g, k3p = rhs_r(g + 0.5 * dt * k2g, pi + 0.5 * dt * k2p)
            k4g, k4p = rhs_r(g + dt * k3g, pi + dt * k3p)
            g = g + dt / 6 * (k1g + 2 * k2g + 2 * k3g + k4g)
            pi = pi + dt / 6 * (k1p + 2 * k2p + 2 * k3p + k4p)
            if k % 5 == 0:
                emax = max(emax, abs(energy_r(g, pi) - E0) / abs(E0))
                gmin_run = min(gmin_run, float(np.min(g)))
        ok_er = emax <= 1e-6 and gmin_run > gf
        print("  P2a-ii floor=%.7f radial: |dE|/E=%.2e (gate<=1e-6), "
              "min g=%.4f > g_floor -> %s"
              % (gf, emax, gmin_run, "PASS" if ok_er else "FAIL"),
              flush=True)
        allpass = allpass and ok_er
        # 3D N=32
        N = 32
        h3 = L3 / N
        x = np.arange(N) * h3
        g3 = 1.0 + 1e-3 * np.cos(x)[:, None, None] * np.ones((N, N, N))
        pi3 = np.zeros_like(g3)

        def energy_3(g_, pi_):
            E = np.sum(pi_ ** 2 / (2 * Ke(g_)) + Ut(g_)) * h3 ** 3
            for ax in range(3):
                gn = np.roll(g_, -1, axis=ax)
                dg = (gn - g_) / h3
                E += 0.5 * np.sum(Ke(0.5 * (g_ + gn)) * dg ** 2) * h3 ** 3
            return float(E)

        def rhs_3(g_, pi_):
            Fg = Ke(g_)
            dH = -pi_ ** 2 * Kep(g_) / (2 * Fg ** 2) + U1t(g_)
            for ax in range(3):
                gn = np.roll(g_, -1, axis=ax)
                gm = 0.5 * (g_ + gn)
                dg = (gn - g_) / h3
                t_flux = Ke(gm) * dg / h3
                t_quad = 0.25 * Kep(gm) * dg ** 2
                dH += -t_flux + t_quad
                dH += np.roll(t_flux + t_quad, 1, axis=ax)
            return pi_ / Fg, -dH

        E0 = energy_3(g3, pi3)
        emax3, gmin3 = 0.0, np.inf
        for k in range(int(round(4.0 / dt))):
            k1g, k1p = rhs_3(g3, pi3)
            k2g, k2p = rhs_3(g3 + 0.5 * dt * k1g, pi3 + 0.5 * dt * k1p)
            k3g, k3p = rhs_3(g3 + 0.5 * dt * k2g, pi3 + 0.5 * dt * k2p)
            k4g, k4p = rhs_3(g3 + dt * k3g, pi3 + dt * k3p)
            g3 = g3 + dt / 6 * (k1g + 2 * k2g + 2 * k3g + k4g)
            pi3 = pi3 + dt / 6 * (k1p + 2 * k2p + 2 * k3p + k4p)
            if k % 5 == 0:
                emax3 = max(emax3, abs(energy_3(g3, pi3) - E0) / abs(E0))
                gmin3 = min(gmin3, float(np.min(g3)))
        ok_e3 = emax3 <= 1e-6 and gmin3 > gf
        print("  P2a-ii floor=%.7f 3D N=32: |dE|/E=%.2e (gate<=1e-6), "
              "min g=%.4f > g_floor -> %s"
              % (gf, emax3, gmin3, "PASS" if ok_e3 else "FAIL"),
              flush=True)
        allpass = allpass and ok_e3
    print("\nP2a GATE: %s" % ("PASS -> P2b" if allpass
                              else "FAIL -> STOP (maszyneria niewazna)"))
    return allpass


# ---------------------------------------------------------------- werdykt
def load_res(jid):
    p = RESDIR + jid + ".json"
    if not os.path.exists(p):
        return None
    with open(p) as f:
        return json.load(f)


def common_subgrid_dev(jid_coarse, jid_fine, spec_c, spec_f):
    dc = np.load(RESDIR + jid_coarse + ".npz")["g"]
    df = np.load(RESDIR + jid_fine + ".npz")["g"]
    if spec_c["start"] == "sol":
        hc, hf = spec_c["h"], spec_f["h"]
        rc = (np.arange(len(dc)) + 0.5) * hc
        rf = (np.arange(len(df)) + 0.5) * hf
        fi = np.interp(rc, rf, df)
        return float(np.max(np.abs(dc - fi)))
    gc = dc[::2, ::2, ::2]
    gfi = df[::3, ::3, ::3]
    return float(np.max(np.abs(gc - gfi)))


def verdict():
    jobs = job_registry()
    print("=" * 78)
    print("WERDYKT Q2 (skladanie; litera LOCKa sec. 3 Phase 2)")
    print("=" * 78)
    # tabela biegow
    pairs = []
    for sector in ("tach", "stab"):
        for start in ("sol", "lat", "noise"):
            for fi in range(3):
                if start == "sol":
                    ids = ["%s_sol_f%d_h%s" % (sector, fi + 1, s)
                           for s in ("025", "0125")]
                else:
                    ids = ["%s_%s_f%d_N%d" % (sector, start, fi + 1, N)
                           for N in N3]
                pairs.append((sector, start, fi, ids))
    tach_nuc_convergent = []
    tach_static_candidates = []
    tach_all_homog = True
    tach_has_breakdown_or_inc = False
    stab_alarm = False
    for sector, start, fi, ids in pairs:
        rr = [load_res(j) for j in ids]
        line = "%s %-5s floor=%.7f (prog %.3f): " % (
            sector, start, FLOORS[fi], THRESHOLDS[fi])
        cells = []
        for j, r in zip(ids, rr):
            if r is None:
                cells.append("%s: BRAK" % j)
                continue
            c = "%s: %s t=%.0f g[%.3f,%.3f] dev=%.3f N%d/seed%d" % (
                j.split("_")[-1], r["status"], r["t_end"], r["gmin"],
                r["gmax"], r["dev_const"], r["n_end"], r["n_seed"])
            if r["status"] == "NUCLEATION":
                c += " t0=%.0f Ndet=%d" % (r["t_nuc"], r["n_det"])
            cells.append(c)
        print(line + " | ".join(cells))
        if any(r is None for r in rr):
            continue
        if sector == "stab":
            if any(r["status"] == "NUCLEATION" for r in rr):
                stab_alarm = True
            continue
        # sektor tachionowy: klasyfikacja
        sts = [r["status"] for r in rr]
        if all(s == "NUCLEATION" for s in sts):
            # wymagane dt/2 na obu siatkach
            rr2 = [load_res(j + "_dt2") for j in ids]
            if all(r2 is not None and r2["status"] == "NUCLEATION"
                   for r2 in rr2):
                dets = [r["n_det"] for r in rr] + [r2["n_det"]
                                                  for r2 in rr2]
                okpm1 = max(dets) - min(dets) <= 1
                print("    dt/2: %s; N_det %s -> zgodnosc +-1: %s"
                      % ([r2["status"] for r2 in rr2], dets,
                         "TAK" if okpm1 else "NIE"))
                if okpm1:
                    tach_nuc_convergent.append((start, fi, dets))
                else:
                    tach_has_breakdown_or_inc = True
            else:
                print("    dt/2: %s -> nukleacja NIEZBIEZNA (dt)"
                      % [None if r2 is None else r2["status"]
                         for r2 in rr2])
                tach_has_breakdown_or_inc = True
            tach_all_homog = False
        elif all(s == "STATIONARY" for s in sts):
            devs = [r["dev_const"] for r in rr]
            spec_c = jobs[ids[0]]
            spec_f = jobs[ids[1]]
            dgrid = common_subgrid_dev(ids[0], ids[1], spec_c, spec_f)
            nonconst = min(devs) >= 0.05
            conv = dgrid <= 5e-3
            print("    STACJONARNE: ||g-const||_inf = %s (wymog >=0.05: "
                  "%s); ||g_c-g_f||_inf(podsiatka) = %.2e (wymog <=5e-3: "
                  "%s)" % (["%.4f" % d for d in devs],
                           "TAK" if nonconst else "NIE", dgrid,
                           "TAK" if conv else "NIE"))
            if nonconst and conv:
                tach_static_candidates.append((start, fi))
                tach_all_homog = False
            elif nonconst and not conv:
                tach_has_breakdown_or_inc = True
                tach_all_homog = False
            # jednorodne stacjonarne: zostaje w tach_all_homog
        else:
            if any(s in ("BREAKDOWN", "INCOMPLETE") for s in sts):
                tach_has_breakdown_or_inc = True
            if any(s == "NUCLEATION" for s in sts):
                print("    nukleacja tylko na czesci siatek -> NIEZBIEZNA")
                tach_has_breakdown_or_inc = True
            tach_all_homog = False
    print("\nP2c (sektor stabilny m^2=+gamma): %s"
          % ("FAIL -- FALSZYWY ALARM detektora => STOP" if stab_alarm
             else "PASS -- zero nukleacji"))
    print("\nKLASYFIKACJA Q2 (litera):")
    if stab_alarm:
        print("  P2c FAIL => detektor niewazny => STOP")
        return
    if tach_nuc_convergent:
        print("  Q2-PASS-NUCLEATION: nukleacja ZBIEZNA (obie siatki + dt/2,"
              " liczba obiektow +-1) w: %s"
              % ["%s floor#%d N_det=%s" % (s, fi + 1, d)
                 for s, fi, d in tach_nuc_convergent])
    elif tach_static_candidates:
        print("  Q2-PASS-STATIC: stan stacjonarny NIESTALY zbiezny w: %s"
              % ["%s floor#%d" % (s, fi + 1)
                 for s, fi in tach_static_candidates])
    elif tach_all_homog and not tach_has_breakdown_or_inc:
        print("  Q2-FAIL: wszystkie starty relaksuja do stanu jednorodnego"
              " na wszystkich podlogach")
    else:
        print("  INCONCLUSIVE: zalamanie nie-nukleacyjne / niezbieznosc /"
              " INCOMPLETE (NIE pozytyw -- litera LOCKa)")
    # zbiorczy npz stanow
    arrs = {}
    for jid in list(jobs.keys()):
        for suff in ("", "_dt2"):
            p = RESDIR + jid + suff + ".npz"
            if os.path.exists(p):
                arrs[jid + suff] = np.load(p)["g"]
    np.savez_compressed(BASE + "Phase2_relaxed_states.npz", **arrs)
    print("\nzapisano stany: %sPhase2_relaxed_states.npz (%d tablic)"
          % (BASE, len(arrs)))


# -------------------------------------------------------------------- main
if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "list"
    if mode == "list":
        for k in job_registry():
            print(k)
    elif mode == "p2a":
        ok = run_p2a()
        sys.exit(0 if ok else 1)
    elif mode == "job":
        run_job(sys.argv[2])
    elif mode == "jobs":
        for j in sys.argv[2:]:
            run_job(j)
    elif mode == "verdict":
        verdict()
    else:
        raise SystemExit("nieznany tryb: %s" % mode)
