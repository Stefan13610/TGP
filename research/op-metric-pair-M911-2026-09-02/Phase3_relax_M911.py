#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
op-metric-pair-M911 (Phase 3) -- Q-B RACHUNEK CENTRALNY: relaksacja
w parze metrycznej (w, V_M9.1'', K=K_geo psi^4).

LOCK: Phase0_balance.md sec. 1-2; decyzje FROZEN: Phase_method_decisions.md:
- PRIMARY (odczyt B, rozstrzygniety CYTATEM -- MD sec. 2):
  E[psi] = int [ 1/2 Keff |grad psi|^2 + Ueff ] dx,
  Keff = K_geo psi^4 (w * g^ij(M9.1'') = 1), Ueff = w*V =
  -gamma psi^3 (4-3psi)/12 = gamma(psi^4/4 - psi^3/3) [tozsamosc
  zweryfikowana sympy, Phase1 P1b PASS]; Ueff' = gamma psi^2 (psi-1).
- Gradient flow: dpsi/dt = -[ -div(Keff grad psi) + 1/2 Keff' |grad|^2
  + Ueff' ]; semi-implicit Euler dt=0.01 (kontrola dt/2), t_max=200,
  stacjonarnosc ||rhs||_inf <= 1e-8 co dt=1.
- ZERO podlog/barier (sedno Q-A); pas graniczny psi > 4/3-1e-6 =
  BREAKDOWN-BOUNDARY (klasyfikacja deskryptywna, stop; MD sec. 4);
  min psi < 1e-6 = BREAKDOWN-BOUNDARY-LOWER; niefinitycznosc = BREAKDOWN.
- Detektory FROZEN (MD sec. 5): dolny psi<5/6, gorny psi>7/6;
  ndimage.label + periodyczne sklejanie; N_seed w t=0; nukleacja =
  N>N_seed utrzymane >=10 j.cz. (11 probek).
- Starty (MD sec. 6): (i) geneza psi=1+szum (seed=20260903, amp=1e-3,
  |n_i|<=8), L=4pi, N in {48,64}; (ii) bump gauss psi=1+0.3 exp(-r^2/
  (2*5^2)), R=60, h in {0.025,0.0125}; (iii) siec 2pi z npz READ-ONLY:
  psi_raw=g^2, psi0=1+s(psi_raw-1), s=0.30/(max psi_raw - 1) per siatka.

Uzycie:
  python Phase3_relax_M911.py list
  python Phase3_relax_M911.py job <id> [--resume]
  python Phase3_relax_M911.py jobs <id1> <id2> ... [--resume]
  python Phase3_relax_M911.py dt2needed
  python Phase3_relax_M911.py verdict

REJESTR WEJSC [INPUT, flagowane]: K_geo=gamma=1; seed=20260903 amp=1e-3;
psi_max startow (ii)/(iii) = 1.30; progi detektorow psi {5/6, 7/6};
pas graniczny 4/3-1e-6; sigma_bump=5.0 [INPUT-MD]; skalowanie sieci
per-siatka FROZEN; tlo 2pi z npz READ-ONLY (mtime weryfikowany).
"""
import json
import os
import sys
import time
import numpy as np
from scipy import ndimage
from scipy.linalg import solve_banded

GAMMA = 1.0                                    # INPUT (LOCK sec. 1)
KGEO = 1.0                                     # INPUT (LOCK sec. 1)
SEED = 20260903                                # INPUT (LOCK sec. 1)
NOISE_AMP = 1e-3                               # INPUT (LOCK sec. 1)
PSI_MAX_START = 1.30                           # INPUT (LOCK sec. 1)
PSI_THR_DN = 5.0 / 6.0                         # LOCK (detektor dolny)
PSI_THR_UP = 7.0 / 6.0                         # LOCK (detektor gorny)
PSI_LIMIT = 4.0 / 3.0
PSI_BAND = PSI_LIMIT - 1e-6                    # LOCK (pas graniczny)
PSI_LOW = 1e-6                                 # MD sec. 4 (deskryptywny)
SIGMA_BUMP = 5.0                               # INPUT-MD (MD rejestr)
AMP_BUMP = 0.30
R_RAD = 60.0
H_RAD = (0.025, 0.0125)                        # LOCK
L_LAT = 2 * np.pi
N_LAT = (32, 48)                               # klucze npz
L_GEN = 4 * np.pi
N_GEN = (48, 64)                               # LOCK
DT_MAIN = 0.01
T_MAX = 200.0
STAT_TOL = 1e-8
CHECKPOINT_DT = 10.0
BASE = ("C:/Users/Mateusz/Documents/ObsydnianMain/TGP/TGP_v1/research/"
        "op-metric-pair-M911-2026-09-02/")
NPZ_BG = ("C:/Users/Mateusz/Documents/ObsydnianMain/TGP/TGP_v1/research/"
          "op-3d-canonical-lattice-2026-08-31/Phase2_backgrounds3d.npz")
RESDIR = BASE + "Phase3_results/"

t0_wall = time.time()


def stamp(msg):
    print("[t=%7.1fs] %s" % (time.time() - t0_wall, msg), flush=True)


def registry_banner(extra=""):
    print("REJESTR [INPUT]: K_geo=gamma=1; seed=%d amp=%g; psi_max_start="
          "%.2f; thr_dn=5/6=%.7f thr_up=7/6=%.7f; pas=4/3-1e-6; "
          "sigma_bump=%.1f [INPUT-MD]; formy: w=psi/(4-3psi) "
          "[eq:vol-element-M911], V=-g*psi^2(4-3psi)^2/12 [eq:V-M911], "
          "Keff=psi^4 [odczyt B, MD sec.2]%s"
          % (SEED, NOISE_AMP, PSI_MAX_START, PSI_THR_DN, PSI_THR_UP,
             SIGMA_BUMP, extra), flush=True)


# ------------------------------------------------- model (FROZEN, MD sec. 3)
def Ueff(p):
    """Ueff = w*V = gamma(psi^4/4 - psi^3/3) -- postac wielomianowa
    DOKLADNA (tozsamosc sympy Phase1 P1b PASS)."""
    return GAMMA * (p ** 4 / 4.0 - p ** 3 / 3.0)


def Ueffp(p):
    return GAMMA * p * p * (p - 1.0)


def Keff(p):
    return KGEO * p ** 4


def Keffp(p):
    return 4.0 * KGEO * p ** 3


class Model:
    Keff = staticmethod(Keff)
    Keffp = staticmethod(Keffp)
    Ueff = staticmethod(Ueff)
    Ueffp = staticmethod(Ueffp)


# ------------------------------------------------- silniki (dziedziczone)
class FlowRadial:
    def __init__(self, h):
        N = int(round(R_RAD / h))
        self.h = h
        self.r = (np.arange(N) + 0.5) * h
        self.r2 = self.r ** 2
        self.rm2 = (0.5 * (self.r[:-1] + self.r[1:])) ** 2
        self.v = Model

    def rhs(self, g):
        h = self.h
        gm = 0.5 * (g[:-1] + g[1:])
        dg = np.diff(g) / h
        dH = h * self.r2 * self.v.Ueffp(g)
        t_flux = self.rm2 * self.v.Keff(gm) * dg
        t_quad = 0.25 * h * self.rm2 * self.v.Keffp(gm) * dg ** 2
        dH[:-1] += -t_flux + t_quad
        dH[1:] += t_flux + t_quad
        return -dH / (h * self.r2)

    def energy(self, g):
        gm = 0.5 * (g[:-1] + g[1:])
        dg = np.diff(g) / self.h
        return float(np.sum(self.r2 * self.v.Ueff(g)) * self.h
                     + 0.5 * np.sum(self.rm2 * self.v.Keff(gm) * dg ** 2)
                     * self.h)

    def step(self, g, dt):
        r = self.rhs(g)
        A = 1.05 * float(np.max(self.v.Keff(g)))
        N = len(g)
        h2 = self.h ** 2
        lo = np.zeros(N)
        di = np.ones(N)
        up = np.zeros(N)
        c = dt * A / h2
        af = self.rm2
        di[:-1] += c * af / self.r2[:-1]
        di[1:] += c * af / self.r2[1:]
        up[1:] = -c * af / self.r2[:-1]
        lo[:-1] = -c * af / self.r2[1:]
        ab = np.vstack([up, di, lo])
        dgv = solve_banded((1, 1), ab, dt * r)
        return g + dgv, r


class Flow3D:
    def __init__(self, N, L):
        self.N = N
        self.h = L / N
        self.v = Model
        k1 = (2 - 2 * np.cos(2 * np.pi * np.fft.fftfreq(N))) / self.h ** 2
        kr = (2 - 2 * np.cos(2 * np.pi * np.fft.rfftfreq(N))) / self.h ** 2
        self.ksym = (k1[:, None, None] + k1[None, :, None]
                     + kr[None, None, :])

    def rhs(self, g):
        h = self.h
        dH = self.v.Ueffp(g)
        for ax in range(3):
            gn = np.roll(g, -1, axis=ax)
            gm = 0.5 * (g + gn)
            dg = (gn - g) / h
            t_flux = self.v.Keff(gm) * dg / h
            t_quad = 0.25 * self.v.Keffp(gm) * dg ** 2
            dH += -t_flux + t_quad
            dH += np.roll(t_flux + t_quad, 1, axis=ax)
        return -dH

    def energy(self, g):
        E = float(np.sum(self.v.Ueff(g))) * self.h ** 3
        for ax in range(3):
            gn = np.roll(g, -1, axis=ax)
            dg = (gn - g) / self.h
            E += 0.5 * float(np.sum(self.v.Keff(0.5 * (g + gn))
                                    * dg ** 2)) * self.h ** 3
        return E

    def step(self, g, dt):
        r = self.rhs(g)
        A = 1.05 * float(np.max(self.v.Keff(g)))
        rh = np.fft.rfftn(dt * r)
        rh /= (1.0 + dt * A * self.ksym)
        return g + np.fft.irfftn(rh, s=g.shape, axes=(0, 1, 2)), r


# --------------------------------------------- detektory (FROZEN, MD sec. 5)
def label_mask(mask):
    """(liczba obiektow, rozmiary); 3D periodycznie (union-find przez
    pary scian), 1D bez periodycznosci. Dziedziczone doslownie."""
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


def detect(psi):
    n_dn, s_dn = label_mask(psi < PSI_THR_DN)
    n_up, s_up = label_mask(psi > PSI_THR_UP)
    return n_dn, s_dn, n_up, s_up


# ------------------------------------------------------------ silnik biegu
def run_flow(flow, psi0, dt, label, tmax=T_MAX, ckpt_path=None,
             resume=False):
    """Gradient flow do stacjonarnosci / nukleacji / t_max / pasa
    granicznego / zalamania. Detektory i klasyfikacje FROZEN."""
    np.seterr(over='ignore', invalid='ignore')
    steps_per_unit = int(round(1.0 / dt))
    nsteps = int(round(tmax / dt))
    k_start = 0
    if resume and ckpt_path and os.path.exists(ckpt_path):
        ck = np.load(ckpt_path, allow_pickle=True)
        g = np.array(ck["g"])
        k_start = int(ck["k"])
        series = json.loads(str(ck["series"]))
        n0_dn, n0_up = int(ck["n0_dn"]), int(ck["n0_up"])
        streak_dn, streak_up = int(ck["streak_dn"]), int(ck["streak_up"])
        stamp("  [%s] RESUME z t=%.1f" % (label, k_start * dt))
    else:
        g = psi0.copy()
        n0_dn, s0_dn, n0_up, s0_up = detect(g)
        series = [dict(t=0.0, n_dn=n0_dn, n_up=n0_up,
                       sz_dn=s0_dn[-5:], sz_up=s0_up[-5:],
                       pmin=float(np.min(g)), pmax=float(np.max(g)),
                       pnorm=float(np.max(np.abs(flow.rhs(g)))),
                       E=flow.energy(g))]
        streak_dn = streak_up = 0
        stamp("  [%s] start: N_seed dn=%d up=%d, psi in [%.4f,%.4f], "
              "E=%.6f, ||psidot||=%.2e"
              % (label, n0_dn, n0_up, series[0]["pmin"],
                 series[0]["pmax"], series[0]["E"], series[0]["pnorm"]))
    status, t_end, t_nuc, n_det, nuc_dir = "TMAX", tmax, None, None, None
    g_prev = g.copy()
    for k in range(k_start + 1, nsteps + 1):
        g_prev = g
        try:
            g, _ = flow.step(g, dt)
            finite = bool(np.all(np.isfinite(g)))
            mx = float(np.max(g)) if finite else np.nan
            mn = float(np.min(g)) if finite else np.nan
        except (ValueError, FloatingPointError, OverflowError):
            finite, mx, mn = False, np.nan, np.nan
        if not finite:
            status, t_end = "BREAKDOWN", k * dt
            g = g_prev
            stamp("  [%s] ZALAMANIE (niefinitycznosc) t=%.3f; ostatni "
                  "finityczny: psi in [%.4g,%.4g]"
                  % (label, t_end, float(np.min(g)), float(np.max(g))))
            break
        if mx >= PSI_BAND:
            status, t_end = "BREAKDOWN-BOUNDARY", k * dt
            stamp("  [%s] PAS GRANICZNY t=%.3f: max psi=%.9f >= 4/3-1e-6"
                  " -- klasyfikacja deskryptywna, stop (MD sec. 4)"
                  % (label, t_end, mx))
            break
        if mn <= PSI_LOW:
            status, t_end = "BREAKDOWN-BOUNDARY-LOWER", k * dt
            stamp("  [%s] DOLNY KRANIEC t=%.3f: min psi=%.3g <= 1e-6"
                  % (label, t_end, mn))
            break
        if k % steps_per_unit == 0:
            t = k * dt
            r = flow.rhs(g)
            pnorm = float(np.max(np.abs(r)))
            n_dn, s_dn, n_up, s_up = detect(g)
            series.append(dict(t=t, n_dn=n_dn, n_up=n_up,
                               sz_dn=s_dn[-5:], sz_up=s_up[-5:],
                               pmin=mn, pmax=mx, pnorm=pnorm,
                               E=flow.energy(g)))
            streak_dn = streak_dn + 1 if n_dn > n0_dn else 0
            streak_up = streak_up + 1 if n_up > n0_up else 0
            fired = None
            if streak_dn >= 11:
                fired = "DN"
            elif streak_up >= 11:
                fired = "UP"
            if fired:
                t_nuc = t - 10.0
                n_det = series[-11]["n_dn" if fired == "DN" else "n_up"]
                status, t_end, nuc_dir = "NUCLEATION", t, fired
                stamp("  [%s] NUKLEACJA %s: t0=%.1f, N_obj(t0)=%d > "
                      "N_seed=%d, utrzymana 10 j.cz."
                      % (label, fired, t_nuc, n_det,
                         n0_dn if fired == "DN" else n0_up))
                break
            if pnorm <= STAT_TOL:
                status, t_end = "STATIONARY", t
                stamp("  [%s] STACJONARNOSC t=%.1f: ||psidot||=%.2e, "
                      "psi in [%.6f,%.6f]" % (label, t, pnorm, mn, mx))
                break
            if ckpt_path and (k % int(round(CHECKPOINT_DT / dt)) == 0):
                np.savez_compressed(ckpt_path, g=g, k=k,
                                    series=json.dumps(series),
                                    n0_dn=n0_dn, n0_up=n0_up,
                                    streak_dn=streak_dn,
                                    streak_up=streak_up)
    if status.startswith("BREAKDOWN"):
        n_dn, s_dn, n_up, s_up = detect(g)
        En = flow.energy(g)
        series.append(dict(t=t_end, n_dn=n_dn, n_up=n_up,
                           sz_dn=s_dn[-5:], sz_up=s_up[-5:],
                           pmin=float(np.min(g)), pmax=float(np.max(g)),
                           pnorm=float(np.max(np.abs(flow.rhs(g))))
                           if np.all(np.isfinite(g)) else float('inf'),
                           E=En if np.isfinite(En) else float('nan')))
    last = series[-1]
    dev_const = 0.5 * (last["pmax"] - last["pmin"])
    E_vac = flow.energy(np.ones_like(g))
    res = dict(label=label, status=status, t_end=t_end, t_nuc=t_nuc,
               nuc_dir=nuc_dir, n_seed_dn=n0_dn, n_seed_up=n0_up,
               n_det=n_det, thr_dn=PSI_THR_DN, thr_up=PSI_THR_UP,
               pmin=last["pmin"], pmax=last["pmax"],
               pmax_vs_43=last["pmax"] - PSI_LIMIT,
               pnorm_end=last["pnorm"], E_end=last["E"],
               E_rel=(last["E"] - E_vac) if np.isfinite(last["E"])
               else None, dev_const=dev_const, mean_end=float(np.mean(g)),
               n_dn_end=last["n_dn"], n_up_end=last["n_up"], series=series)
    stamp("  [%s] KONIEC: %s%s (t=%.1f), psi in [%.6f,%.6f] "
          "(pmax-4/3=%+.4f), ||psi-const||=%.6f, mean=%.6f, "
          "N_dn=%d/seed%d N_up=%d/seed%d, ||psidot||=%.2e, E_rel=%s"
          % (label, status, "-" + nuc_dir if nuc_dir else "", t_end,
             last["pmin"], last["pmax"], last["pmax"] - PSI_LIMIT,
             dev_const, res["mean_end"], last["n_dn"], n0_dn,
             last["n_up"], n0_up, last["pnorm"],
             "%.3e" % res["E_rel"] if res["E_rel"] is not None else "nan"))
    return res, g


# ------------------------------------------------------------- starty
def noise_field(N):
    """Szum pasmowy dziedziczony verbatim (MD sec. 6), seed=20260903."""
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

    f64 = build(64)
    scale = NOISE_AMP / float(np.max(np.abs(f64)))
    if N == 64:
        return f64 * scale
    return build(N) * scale


def start_gen(N):
    return 1.0 + noise_field(N)


def start_bump(h):
    N = int(round(R_RAD / h))
    r = (np.arange(N) + 0.5) * h
    return 1.0 + AMP_BUMP * np.exp(-r ** 2 / (2.0 * SIGMA_BUMP ** 2))


def start_lat(N):
    """psi_raw = g^2; psi0 = 1 + s(psi_raw - 1), s = 0.30/(max-1)
    per siatka (procedura FROZEN, MD sec. 6). Raport ksztaltu oryg."""
    mt0 = os.path.getmtime(NPZ_BG)
    data = np.load(NPZ_BG)
    g = np.array(data["2pi__A1.0__N%d" % N])
    data.close()
    mt1 = os.path.getmtime(NPZ_BG)
    assert mt0 == mt1, "npz mtime zmieniony -- STOP"
    psi_raw = g ** 2
    s = (PSI_MAX_START - 1.0) / (float(np.max(psi_raw)) - 1.0)
    psi0 = 1.0 + s * (psi_raw - 1.0)
    stamp("  [lat N=%d] oryginal: g in [%.6f,%.6f], psi_raw in "
          "[%.6f,%.6f]; s=%.7f -> psi0 in [%.6f,%.6f] (max=1.30 FROZEN);"
          " npz mtime OK (READ-ONLY)"
          % (N, float(np.min(g)), float(np.max(g)),
             float(np.min(psi_raw)), float(np.max(psi_raw)), s,
             float(np.min(psi0)), float(np.max(psi0))))
    return psi0


# ------------------------------------------------------------- rejestr biegow
def job_registry():
    jobs = {}
    for N in N_GEN:
        jobs["gen_N%d" % N] = dict(start="gen", N=N)
    for h in H_RAD:
        jobs["bump_h%s" % ("025" if h == 0.025 else "0125")] = \
            dict(start="bump", h=h)
    for N in N_LAT:
        jobs["lat_N%d" % N] = dict(start="lat", N=N)
    return jobs


PAIRS = [("gen", ["gen_N48", "gen_N64"]),
         ("bump", ["bump_h025", "bump_h0125"]),
         ("lat", ["lat_N32", "lat_N48"])]


def build_flow_and_start(spec):
    if spec["start"] == "gen":
        return Flow3D(spec["N"], L_GEN), start_gen(spec["N"]), \
            dict(kind="3d", N=spec["N"], L=float(L_GEN))
    if spec["start"] == "bump":
        return FlowRadial(spec["h"]), start_bump(spec["h"]), \
            dict(kind="radial", h=spec["h"], R=R_RAD)
    return Flow3D(spec["N"], L_LAT), start_lat(spec["N"]), \
        dict(kind="3d", N=spec["N"], L=float(L_LAT))


def run_job(jid, resume=False):
    os.makedirs(RESDIR, exist_ok=True)
    jobs = job_registry()
    base_id = jid[:-4] if jid.endswith("_dt2") else jid
    dt = DT_MAIN / (2 if jid.endswith("_dt2") else 1)
    spec = jobs[base_id]
    flow, psi0, geom = build_flow_and_start(spec)
    registry_banner("; job=%s dt=%g" % (jid, dt))
    label = jid + (" dt=%.4g" % dt)
    ckpt = RESDIR + jid + "_ckpt.npz"
    res, gfin = run_flow(flow, psi0, dt, label, ckpt_path=ckpt,
                         resume=resume)
    res.update(job=jid, dt=dt, geom=geom)
    with open(RESDIR + jid + ".json", "w") as f:
        json.dump(res, f)
    np.savez_compressed(RESDIR + jid + ".npz", psi=gfin,
                        meta=np.array([dt], dtype=float))
    if os.path.exists(ckpt):
        os.remove(ckpt)
    stamp("zapisano: %s.json / .npz" % (RESDIR + jid))
    return res


# ---------------------------------------------------------------- werdykt
def load_res(jid):
    p = RESDIR + jid + ".json"
    if not os.path.exists(p):
        return None
    with open(p) as f:
        return json.load(f)


def common_subgrid_dev(jid_c, jid_f, spec_c, spec_f):
    dc = np.load(RESDIR + jid_c + ".npz")["psi"]
    df = np.load(RESDIR + jid_f + ".npz")["psi"]
    if spec_c["start"] == "bump":
        rc = (np.arange(len(dc)) + 0.5) * spec_c["h"]
        rf = (np.arange(len(df)) + 0.5) * spec_f["h"]
        return float(np.max(np.abs(dc - np.interp(rc, rf, df))))
    sc = spec_c["N"] // 16
    sf = spec_f["N"] // 16
    return float(np.max(np.abs(dc[::sc, ::sc, ::sc]
                               - df[::sf, ::sf, ::sf])))


def dt2_needed():
    """Pary z NUCLEATION lub BREAKDOWN-BOUNDARY* w biegu glownym
    (MD sec. 6: dt/2 'przy zdarzeniach')."""
    need = []
    for start, ids in PAIRS:
        rr = [load_res(j) for j in ids]
        if any(r is not None and (r["status"] == "NUCLEATION"
                                  or r["status"].startswith(
                                      "BREAKDOWN-BOUNDARY"))
               for r in rr):
            need += [j + "_dt2" for j in ids]
    return need


def verdict():
    jobs = job_registry()
    out = []

    def em(s=""):
        print(s, flush=True)
        out.append(s)

    em("=" * 78)
    em("WERDYKT Q-B (skladanie; litera LOCKa sec. 2 Phase 3 + MD sec. 6)")
    em("REJESTR [INPUT]: K_geo=gamma=1; seed=%d amp=%g; psi_max_start=%.2f;"
       % (SEED, NOISE_AMP, PSI_MAX_START))
    em("  thr_dn=5/6=%.7f thr_up=7/6=%.7f; pas=4/3-1e-6; sigma_bump=%.1f"
       " [INPUT-MD]; formy CYTAT sek08a (MD sec.1); odczyt B (MD sec.2)"
       % (PSI_THR_DN, PSI_THR_UP, SIGMA_BUMP))
    em("=" * 78)
    nuc_convergent = []
    static_candidates = []
    boundary_convergent = []
    all_homog = True
    has_inc = False
    lat_note = None
    for start, ids in PAIRS:
        rr = [load_res(j) for j in ids]
        em("%-4s:" % start)
        for j, r in zip(ids, rr):
            if r is None:
                em("    %s: BRAK" % j)
                continue
            c = ("    %s: %s%s t=%.1f psi[%.4f,%.4f] pmax-4/3=%+.4f "
                 "dev=%.4f mean=%.4f Ndn%d/%d Nup%d/%d E_rel=%s"
                 % (j, r["status"],
                    "-" + r["nuc_dir"] if r.get("nuc_dir") else "",
                    r["t_end"], r["pmin"], r["pmax"], r["pmax_vs_43"],
                    r["dev_const"], r["mean_end"], r["n_dn_end"],
                    r["n_seed_dn"], r["n_up_end"], r["n_seed_up"],
                    "%.3e" % r["E_rel"] if r.get("E_rel") is not None
                    else "nan"))
            if r["status"] == "NUCLEATION":
                c += " t0=%.0f Ndet=%d" % (r["t_nuc"], r["n_det"])
            em(c)
        if any(r is None for r in rr):
            has_inc = True
            continue
        sts = [r["status"] for r in rr]
        if all(s == "NUCLEATION" for s in sts):
            all_homog = False
            dirs = set(r["nuc_dir"] for r in rr)
            rr2 = [load_res(j + "_dt2") for j in ids]
            if len(dirs) == 1 and all(
                    r2 is not None and r2["status"] == "NUCLEATION"
                    and r2["nuc_dir"] == rr[0]["nuc_dir"] for r2 in rr2):
                dets = [r["n_det"] for r in rr] + [r2["n_det"]
                                                  for r2 in rr2]
                okpm1 = max(dets) - min(dets) <= 1
                em("    dt/2: NUCLEATION-%s x2; N_det %s -> zgodnosc +-1:"
                   " %s" % (rr[0]["nuc_dir"], dets,
                            "TAK" if okpm1 else "NIE"))
                if okpm1:
                    nuc_convergent.append((start, rr[0]["nuc_dir"], dets))
                else:
                    has_inc = True
            else:
                em("    dt/2/kierunek: %s -> nukleacja NIEZBIEZNA"
                   % [None if r2 is None else (r2["status"],
                                               r2.get("nuc_dir"))
                      for r2 in rr2])
                has_inc = True
        elif all(s == "STATIONARY" for s in sts):
            devs = [r["dev_const"] for r in rr]
            dgrid = common_subgrid_dev(ids[0], ids[1], jobs[ids[0]],
                                       jobs[ids[1]])
            nonconst = min(devs) >= 0.05
            conv = dgrid <= 5e-3
            em("    STACJONARNE: ||psi-const||_inf = %s (>=0.05: %s); "
               "||c-f||_inf(podsiatka) = %.2e (<=5e-3: %s); "
               "sanity pmax-4/3: %s; mean vs psi*=1: %s"
               % (["%.6f" % d for d in devs],
                  "TAK" if nonconst else "NIE", dgrid,
                  "TAK" if conv else "NIE",
                  ["%+.4f" % r["pmax_vs_43"] for r in rr],
                  ["%+.2e" % (r["mean_end"] - 1.0) for r in rr]))
            if nonconst and conv:
                static_candidates.append((start, devs, dgrid))
                all_homog = False
            elif nonconst and not conv:
                has_inc = True
                all_homog = False
        elif all(s == sts[0] and s.startswith("BREAKDOWN-BOUNDARY")
                 for s in sts):
            all_homog = False
            has_inc = True
            rr2 = [load_res(j + "_dt2") for j in ids]
            same_dt2 = all(r2 is not None and r2["status"] == sts[0]
                           for r2 in rr2)
            em("    %s na OBU siatkach (t=%s)%s -> kategoria deskryptywna"
               " 'pole wybiera granice' (NIE pozytyw)"
               % (sts[0], ["%.2f" % r["t_end"] for r in rr],
                  "; dt/2: potwierdzone" if same_dt2 else
                  "; dt/2: %s" % [None if r2 is None else r2["status"]
                                  for r2 in rr2]))
            if same_dt2:
                boundary_convergent.append((start, sts[0]))
        else:
            all_homog = False
            has_inc = True
            em("    statusy mieszane %s -> wklad INCONCLUSIVE" % sts)
        if start == "lat":
            lat_note = (sts if all(r is not None for r in rr) else None)
    em()
    em("DESKRYPTYWNIE (obowiazkowe, LOCK sec.2): los startu sieciowego:")
    for j in ["lat_N32", "lat_N48"]:
        r = load_res(j)
        if r is None:
            em("    %s: BRAK" % j)
            continue
        em("    %s: %s; struktura sieci %s (dev koncowy %.6f, "
           "N_dn/N_up koncowe %d/%d)"
           % (j, r["status"],
              "PRZEZYWA (stan niestaly)" if (r["status"] == "STATIONARY"
                                             and r["dev_const"] >= 0.05)
              else ("NIE przezywa -- relaksacja do jednorodnosci"
                    if r["status"] == "STATIONARY" else "los: " +
                    r["status"]), r["dev_const"], r["n_dn_end"],
              r["n_up_end"]))
    em()
    em("KLASYFIKACJA Q-B (litera):")
    if nuc_convergent:
        em("  Q-B-PASS-NUCLEATION: nukleacja ZBIEZNA (obie siatki + dt/2,"
           " +-1) w: %s" % ["%s kier=%s N_det=%s" % t
                            for t in nuc_convergent])
        qb = "Q-B-PASS-NUCLEATION"
    elif static_candidates:
        em("  Q-B-PASS-STATIC: stan stacjonarny NIESTALY zbiezny w: %s"
           % ["%s dev=%s dgrid=%.2e" % (s, ["%.4f" % d for d in dv], dg)
              for s, dv, dg in static_candidates])
        qb = "Q-B-PASS-STATIC"
    elif all_homog and not has_inc:
        em("  Q-B-FAIL: wszystko relaksuje do jednorodnej prozni")
        qb = "Q-B-FAIL"
    else:
        em("  Q-B-INCONCLUSIVE (NIE pozytyw -- litera LOCKa)%s"
           % ("; kategoria deskryptywna BREAKDOWN-BOUNDARY zbiezna w: %s"
              % boundary_convergent if boundary_convergent else ""))
        qb = "Q-B-INCONCLUSIVE"
    arrs = {}
    for jid in list(jobs.keys()):
        for suff in ("", "_dt2"):
            p = RESDIR + jid + suff + ".npz"
            if os.path.exists(p):
                arrs[jid + suff] = np.load(p)["psi"]
    np.savez_compressed(BASE + "Phase3_relaxed_states.npz", **arrs)
    em()
    em("zapisano stany: %sPhase3_relaxed_states.npz (%d tablic)"
       % (BASE, len(arrs)))
    em("WERDYKT: %s" % qb)
    with open(BASE + "Phase3_output.txt", "w") as f:
        f.write("\n".join(out) + "\n")
    print("zapisano:", BASE + "Phase3_output.txt")


# -------------------------------------------------------------------- main
if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "list"
    if mode == "list":
        for k in job_registry():
            print(k)
    elif mode == "job":
        run_job(sys.argv[2], resume="--resume" in sys.argv)
    elif mode == "jobs":
        for j in sys.argv[2:]:
            if j == "--resume":
                continue
            run_job(j, resume="--resume" in sys.argv)
    elif mode == "dt2needed":
        for j in dt2_needed():
            print(j)
    elif mode == "verdict":
        verdict()
    else:
        raise SystemExit("nieznany tryb: %s" % mode)
