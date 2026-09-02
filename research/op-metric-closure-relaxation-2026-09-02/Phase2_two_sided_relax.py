#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
op-metric-closure-relaxation (Phase 2) -- RACHUNEK CENTRALNY: relaksacja
z OBUSTRONNYM domknieciem (podloga QB-2 + granica metryczna psi=4/3).

LOCK: Phase0_balance.md sec. 2-3; decyzje FROZEN: Phase_method_decisions.md:
- PRIMARY: E[g] = int w(psi)[1/2 K|grad g|^2 + U_b] dx, w(psi)=psi/(4-3psi)
  (eq:vol-element-M911), psi=g^2, K=g^4, U_b = U - U(1) + V_fl,
  U = b g^7/7 - c g^8/8; referencja prozniowa U(1)=1/56 (MD sec. 1);
  wariacja: rhs = -[-div(Keff grad g) + 1/2 Keff' |grad g|^2 + Ueff'],
  Keff = W g^4, Ueff = W U_b, W(g)=g^2/(4-3g^2), W'=8g/(4-3g^2)^2;
  regularizacja bieguna TYLKO powyzej g_c = g_ceil - 1e-6 (plaska
  kontynuacja W, W'=0 -- MD sec. 3).
- C-BAR: w=1, U_b = U - U(1) + V_fl + V_ceil,
  V_ceil = (kappa/3)(g-g_ceil)^3 dla g>g_ceil, kappa=100.
- Podloga dziedziczona: V_fl = (kappa/3)(g_floor-g)^3 dla g<g_floor.
- Detektory: dolny maska g<(1+g_floor)/2 (dziedziczony), gorny maska
  g>(1+g_ceil)/2=1.0773503 (LOCK); ndimage.label, periodyczne sklejanie,
  N_seed w t=0, nukleacja = N>N_seed utrzymane >=10 j.cz.
- Flow: semi-implicit Euler dt=0.01 (dt/2 kontrola), t_max=200,
  stacjonarnosc ||rhs||_inf <= 1e-8.

Uzycie:
  python Phase2_two_sided_relax.py list
  python Phase2_two_sided_relax.py job <id>
  python Phase2_two_sided_relax.py jobs <id1> <id2> ...
  python Phase2_two_sided_relax.py verdict

REJESTR WEJSC [INPUT, flagowane]: g_ceil=sqrt(4/3) [M9.1'' CYTAT:
sqrt(-g_eff)=c0 psi/(4-3psi), eq:vol-element-M911, sek08a]; progi QB-2
{0.197,0.298,0.331} -> g_floor=sqrt(.) (PRIMARY=0.5458938; wrazliwosc
tylko geneza); kappa=100 FROZEN; seed=20260902, amp=1e-3;
g0_mu=phi*0.90548=1.4650974; beta=gamma=1; tlo 2pi z npz READ-ONLY.
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
G0_MU = PHI_GOLD * 0.90548                    # INPUT
KAPPA = 100.0                                 # FROZEN (dziedziczone)
SEED = 20260902                               # INPUT (LOCK sec. 2)
NOISE_AMP = 1e-3                              # INPUT (LOCK sec. 2)
THRESHOLDS = (0.197, 0.298, 0.331)            # INPUT QB-2
FLOORS = tuple(np.sqrt(t) for t in THRESHOLDS)
FI_PRIMARY = 1                                # f2 = 0.5458938 (PRIMARY)
G_CEIL = float(np.sqrt(4.0 / 3.0))            # INPUT M9.1'' (CYTAT w naglowku)
G_C = G_CEIL - 1e-6                           # prog regularizacji (LOCK)
PSI_C = G_C ** 2
W_CLIP = PSI_C / (4.0 - 3.0 * PSI_C)
GTHR_UP = 0.5 * (1.0 + G_CEIL)                # detektor gorny (LOCK)
U_VAC_TACH = BETA / 7.0 - GAMMA / 8.0         # U(1) = 1/56
R_RAD = 60.0
H_RAD = (0.025, 0.0125)                       # FROZEN (dziedziczone)
L_LAT = 2 * np.pi
N_LAT = (32, 48)
L_GEN = 4 * np.pi                             # LOCK: wieksze pudlo
N_GEN = (48, 64)
DT_MAIN = 0.01
T_MAX = 200.0
STAT_TOL = 1e-8
CHECKPOINT_DT = 10.0
BASE = ("C:/Users/Mateusz/Documents/ObsydnianMain/TGP/TGP_v1/research/"
        "op-metric-closure-relaxation-2026-09-02/")
NPZ_BG = ("C:/Users/Mateusz/Documents/ObsydnianMain/TGP/TGP_v1/research/"
          "op-3d-canonical-lattice-2026-08-31/Phase2_backgrounds3d.npz")
RESDIR = BASE + "Phase2_results/"

t0_wall = time.time()


def stamp(msg):
    print("[t=%7.1fs] %s" % (time.time() - t0_wall, msg), flush=True)


def registry_banner(gf, variant, sector):
    print("REJESTR [INPUT]: g_ceil=sqrt(4/3)=%.7f [M9.1'' eq:vol-element-"
          "M911]; g_floor=%.7f; kappa=%.0f FROZEN; seed=%d amp=%g; "
          "g0_mu=%.7f; beta=gamma=1; wariant=%s sektor=%s"
          % (G_CEIL, gf, KAPPA, SEED, NOISE_AMP, G0_MU, variant, sector),
          flush=True)


# ---------------------------------------------------- model: potencjaly
def U_parts(sector):
    """'tach': U = b g^7/7 - c g^8/8 (m^2(1)=-gamma), U(1)=+1/56;
    'stab': U = -b g^7/7 + c g^8/8 (m^2(1)=+gamma), U(1)=-1/56."""
    s = 1.0 if sector == "tach" else -1.0
    Uvac = s * U_VAC_TACH
    U = lambda g: s * (BETA * g ** 7 / 7 - GAMMA * g ** 8 / 8) - Uvac
    U1 = lambda g: s * (BETA * g ** 6 - GAMMA * g ** 7)
    return U, U1                       # U juz z referencja prozniowa


def vfl_parts(gf):
    """Podloga C^2 (dziedziczona): V=(kappa/3)(gf-g)^3 dla g<gf."""
    def V(g):
        d = np.maximum(gf - g, 0.0)
        return (KAPPA / 3.0) * d ** 3
    def V1(g):
        d = np.maximum(gf - g, 0.0)
        return -KAPPA * d ** 2
    return V, V1


def vceil_parts():
    """Bariera C-BAR C^2: V=(kappa/3)(g-g_ceil)^3 dla g>g_ceil."""
    def V(g):
        d = np.maximum(g - G_CEIL, 0.0)
        return (KAPPA / 3.0) * d ** 3
    def V1(g):
        d = np.maximum(g - G_CEIL, 0.0)
        return KAPPA * d ** 2
    return V, V1


def W_parts():
    """W(g)=w(g^2)=g^2/(4-3g^2); regularizacja: dla g^2>=psi_c
    W=W_clip (stala), W'=0 (MD sec. 3; TYLKO powyzej g_ceil-1e-6)."""
    def W(g):
        psi = g * g
        den = 4.0 - 3.0 * np.minimum(psi, PSI_C)
        return np.where(psi < PSI_C, psi / den, W_CLIP)
    def W1(g):
        psi = g * g
        den = 4.0 - 3.0 * np.minimum(psi, PSI_C)
        return np.where(psi < PSI_C, 8.0 * g / den ** 2, 0.0)
    return W, W1


class Variant:
    """Buduje (Keff, Keffp, Ueff, Ueffp) dla wariantu domkniecia.
    variant in {'PRIM','CBAR','NONE'}; NONE = bez domkniec (P1b)."""

    def __init__(self, variant, sector, fi):
        self.name = variant
        U, U1 = U_parts(sector)
        if variant == "NONE":
            Ub, Ub1 = U, U1
        else:
            Vf, Vf1 = vfl_parts(FLOORS[fi])
            if variant == "CBAR":
                Vc, Vc1 = vceil_parts()
                Ub = lambda g: U(g) + Vf(g) + Vc(g)
                Ub1 = lambda g: U1(g) + Vf1(g) + Vc1(g)
            else:
                Ub = lambda g: U(g) + Vf(g)
                Ub1 = lambda g: U1(g) + Vf1(g)
        if variant == "PRIM":
            W, W1 = W_parts()
            self.Keff = lambda g: W(g) * g ** 4
            self.Keffp = lambda g: W1(g) * g ** 4 + W(g) * 4 * g ** 3
            self.Ueff = lambda g: W(g) * Ub(g)
            self.Ueffp = lambda g: W1(g) * Ub(g) + W(g) * Ub1(g)
        else:
            self.Keff = lambda g: g ** 4
            self.Keffp = lambda g: 4 * g ** 3
            self.Ueff = Ub
            self.Ueffp = Ub1


# --------------------------------------------------------- profil solitonu
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
    """Pasmowo ograniczony szum |n_i|<=8 (w L=4pi: k=n/2<=4),
    rng(20260902); te same wspolczynniki na obu siatkach;
    skala z max|f| na N=64 (MD sec. 5)."""
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


# --------------------------------------------------------------- silniki
class FlowRadial:
    def __init__(self, h, var):
        N = int(round(R_RAD / h))
        self.h = h
        self.r = (np.arange(N) + 0.5) * h
        self.r2 = self.r ** 2
        self.rm2 = (0.5 * (self.r[:-1] + self.r[1:])) ** 2
        self.v = var
        self.a_face = self.rm2

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
        af = self.a_face
        di[:-1] += c * af / self.r2[:-1]
        di[1:] += c * af / self.r2[1:]
        up[1:] = -c * af / self.r2[:-1]
        lo[:-1] = -c * af / self.r2[1:]
        ab = np.vstack([up, di, lo])
        dgv = solve_banded((1, 1), ab, dt * r)
        return g + dgv, r


class Flow3D:
    def __init__(self, N, L, var):
        self.N = N
        self.h = L / N
        self.v = var
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


# ------------------------------------------------------ detektory (FROZEN)
def label_mask(mask):
    """(liczba obiektow, rozmiary) maski; 3D periodycznie (union-find
    przez pary scian), 1D bez periodycznosci. Dziedziczone doslownie."""
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


def detect(g, gthr_dn):
    n_dn, s_dn = label_mask(g < gthr_dn)
    n_up, s_up = label_mask(g > GTHR_UP)
    return n_dn, s_dn, n_up, s_up


# ------------------------------------------------------------ silnik biegu
def run_flow(flow, g0, gf, dt, label, tmax=T_MAX, ckpt_path=None,
             resume=False):
    """Gradient flow do stacjonarnosci / nukleacji (dolnej LUB gornej) /
    t_max / zalamania. Detektory FROZEN (MD sec. 4)."""
    np.seterr(over='ignore', invalid='ignore')
    gthr_dn = 0.5 * (1.0 + gf)
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
        g = g0.copy()
        n0_dn, s0_dn, n0_up, s0_up = detect(g, gthr_dn)
        series = [dict(t=0.0, n_dn=n0_dn, n_up=n0_up,
                       sz_dn=s0_dn[-5:], sz_up=s0_up[-5:],
                       gmin=float(np.min(g)), gmax=float(np.max(g)),
                       gnorm=float(np.max(np.abs(flow.rhs(g)))),
                       E=flow.energy(g))]
        streak_dn = streak_up = 0
        stamp("  [%s] start: N_seed dn=%d up=%d, g in [%.4f,%.4f], "
              "E=%.6f, ||gdot||=%.2e"
              % (label, n0_dn, n0_up, series[0]["gmin"], series[0]["gmax"],
                 series[0]["E"], series[0]["gnorm"]))
    status, t_end, t_nuc, n_det, nuc_dir = "INCOMPLETE", tmax, None, None, None
    g_prev = g.copy()
    for k in range(k_start + 1, nsteps + 1):
        g_prev = g
        try:
            g, _ = flow.step(g, dt)
            mx = float(np.max(g)) if np.all(np.isfinite(g)) else np.nan
        except (ValueError, FloatingPointError, OverflowError):
            mx = np.nan
        if not np.isfinite(mx):
            status, t_end = "BREAKDOWN", k * dt
            g = g_prev
            stamp("  [%s] ZALAMANIE (niefinitycznosc) t=%.3f; ostatni "
                  "finityczny stan: g in [%.4g, %.4g]"
                  % (label, t_end, float(np.min(g)), float(np.max(g))))
            break
        if k % steps_per_unit == 0:
            t = k * dt
            r = flow.rhs(g)
            gnorm = float(np.max(np.abs(r)))
            n_dn, s_dn, n_up, s_up = detect(g, gthr_dn)
            series.append(dict(t=t, n_dn=n_dn, n_up=n_up,
                               sz_dn=s_dn[-5:], sz_up=s_up[-5:],
                               gmin=float(np.min(g)), gmax=mx,
                               gnorm=gnorm, E=flow.energy(g)))
            streak_dn = streak_dn + 1 if n_dn > n0_dn else 0
            streak_up = streak_up + 1 if n_up > n0_up else 0
            fired = None
            if streak_dn >= 11:
                fired, nn = "DN", [s["n_dn"] for s in series]
            elif streak_up >= 11:
                fired, nn = "UP", [s["n_up"] for s in series]
            if fired:
                t_nuc = t - 10.0
                n_det = series[-11]["n_dn" if fired == "DN" else "n_up"]
                status, t_end, nuc_dir = "NUCLEATION", t, fired
                stamp("  [%s] NUKLEACJA %s: t0=%.1f, N_obj(t0)=%d > "
                      "N_seed=%d, utrzymana 10 j.cz."
                      % (label, fired, t_nuc, n_det,
                         n0_dn if fired == "DN" else n0_up))
                break
            if gnorm <= STAT_TOL:
                status, t_end = "STATIONARY", t
                stamp("  [%s] STACJONARNOSC t=%.1f: ||gdot||=%.2e, g in "
                      "[%.4f,%.4f]" % (label, t, gnorm,
                                       series[-1]["gmin"], mx))
                break
            if ckpt_path and (k % int(round(CHECKPOINT_DT / dt)) == 0):
                np.savez_compressed(ckpt_path, g=g, k=k,
                                    series=json.dumps(series),
                                    n0_dn=n0_dn, n0_up=n0_up,
                                    streak_dn=streak_dn,
                                    streak_up=streak_up)
    if status == "BREAKDOWN":
        n_dn, s_dn, n_up, s_up = detect(g, gthr_dn)
        series.append(dict(t=t_end, n_dn=n_dn, n_up=n_up,
                           sz_dn=s_dn[-5:], sz_up=s_up[-5:],
                           gmin=float(np.min(g)), gmax=float(np.max(g)),
                           gnorm=float('inf'), E=float('nan')))
    last = series[-1]
    dev_const = 0.5 * (last["gmax"] - last["gmin"])
    res = dict(label=label, status=status, t_end=t_end, t_nuc=t_nuc,
               nuc_dir=nuc_dir, n_seed_dn=n0_dn, n_seed_up=n0_up,
               n_det=n_det, gthr_dn=gthr_dn, gthr_up=GTHR_UP,
               gmin=last["gmin"], gmax=last["gmax"],
               gmax_vs_gceil=last["gmax"] - G_CEIL,
               gnorm_end=last["gnorm"], E_end=last["E"],
               dev_const=dev_const, n_dn_end=last["n_dn"],
               n_up_end=last["n_up"], series=series)
    stamp("  [%s] KONIEC: %s%s (t=%.1f), g in [%.4f,%.4f] "
          "(gmax-gceil=%+.4f), ||g-const||=%.4f, N_dn=%d/seed%d "
          "N_up=%d/seed%d, ||gdot||=%.2e"
          % (label, status, "-" + nuc_dir if nuc_dir else "", t_end,
             last["gmin"], last["gmax"], last["gmax"] - G_CEIL, dev_const,
             last["n_dn"], n0_dn, last["n_up"], n0_up, last["gnorm"]))
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


def start_gen(N):
    return 1.0 + noise_field(N)


# ------------------------------------------------------------- rejestr biegow
def job_registry():
    """Macierz P2 (MD sec. 5): 14 biegow glownych, sektor tachionowy."""
    jobs = {}
    for fi in range(3):
        for N in N_GEN:
            jobs["gen_PRIM_f%d_N%d" % (fi + 1, N)] = dict(
                sector="tach", start="gen", variant="PRIM", fi=fi, N=N)
    for variant in ("PRIM", "CBAR"):
        for h in H_RAD:
            jobs["sol_%s_f2_h%s" % (variant,
                                    "025" if h == 0.025 else "0125")] = \
                dict(sector="tach", start="sol", variant=variant,
                     fi=FI_PRIMARY, h=h)
        for N in N_LAT:
            jobs["lat_%s_f2_N%d" % (variant, N)] = dict(
                sector="tach", start="lat", variant=variant,
                fi=FI_PRIMARY, N=N)
    return jobs


def build_flow_and_start(spec):
    var = Variant(spec["variant"], spec["sector"], spec["fi"])
    if spec["start"] == "sol":
        flow = FlowRadial(spec["h"], var)
        g0 = start_radial_sol(spec["h"])
        geom = dict(kind="radial", h=spec["h"], R=R_RAD)
    elif spec["start"] == "lat":
        flow = Flow3D(spec["N"], L_LAT, var)
        g0 = start_lat(spec["N"])
        geom = dict(kind="3d", N=spec["N"], L=L_LAT)
    else:
        flow = Flow3D(spec["N"], L_GEN, var)
        g0 = start_gen(spec["N"])
        geom = dict(kind="3d", N=spec["N"], L=L_GEN)
    return flow, g0, geom


def run_job(jid, resume=False):
    os.makedirs(RESDIR, exist_ok=True)
    jobs = job_registry()
    base_id = jid[:-4] if jid.endswith("_dt2") else jid
    dt = DT_MAIN / (2 if jid.endswith("_dt2") else 1)
    spec = jobs[base_id]
    flow, g0, geom = build_flow_and_start(spec)
    gf = FLOORS[spec["fi"]]
    registry_banner(gf, spec["variant"], spec["sector"])
    label = jid + (" dt=%.4g" % dt)
    ckpt = RESDIR + jid + "_ckpt.npz"
    res, gfin = run_flow(flow, g0, gf, dt, label, ckpt_path=ckpt,
                         resume=resume)
    res.update(job=jid, dt=dt, geom=geom, g_floor=gf,
               threshold=THRESHOLDS[spec["fi"]], variant=spec["variant"])
    with open(RESDIR + jid + ".json", "w") as f:
        json.dump(res, f)
    np.savez_compressed(RESDIR + jid + ".npz", g=gfin,
                        meta=np.array([gf, dt], dtype=float))
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


def common_subgrid_dev(jid_coarse, jid_fine, spec_c, spec_f):
    dc = np.load(RESDIR + jid_coarse + ".npz")["g"]
    df = np.load(RESDIR + jid_fine + ".npz")["g"]
    if spec_c["start"] == "sol":
        hc, hf = spec_c["h"], spec_f["h"]
        rc = (np.arange(len(dc)) + 0.5) * hc
        rf = (np.arange(len(df)) + 0.5) * hf
        fi = np.interp(rc, rf, df)
        return float(np.max(np.abs(dc - fi)))
    sc = spec_c["N"] // 16
    sf = spec_f["N"] // 16
    gc = dc[::sc, ::sc, ::sc]
    gfi = df[::sf, ::sf, ::sf]
    return float(np.max(np.abs(gc - gfi)))


def verdict():
    jobs = job_registry()
    print("=" * 78)
    print("WERDYKT Q (skladanie; litera LOCKa sec. 3 Phase 2)")
    print("REJESTR [INPUT]: g_ceil=%.7f (M9.1'' eq:vol-element-M911); "
          "g_thr_up=%.7f; g_floor PRIMARY=%.7f; kappa=%.0f; seed=%d"
          % (G_CEIL, GTHR_UP, FLOORS[FI_PRIMARY], KAPPA, SEED))
    print("=" * 78)
    pairs = []
    for fi in range(3):
        pairs.append(("gen", "PRIM", fi,
                      ["gen_PRIM_f%d_N%d" % (fi + 1, N) for N in N_GEN]))
    for variant in ("PRIM", "CBAR"):
        pairs.append(("sol", variant, FI_PRIMARY,
                      ["sol_%s_f2_h%s" % (variant, s)
                       for s in ("025", "0125")]))
        pairs.append(("lat", variant, FI_PRIMARY,
                      ["lat_%s_f2_N%d" % (variant, N) for N in N_LAT]))
    nuc_convergent = []
    static_candidates = []
    all_homog = True
    has_breakdown_or_inc = False
    for start, variant, fi, ids in pairs:
        rr = [load_res(j) for j in ids]
        line = "%-4s %-4s floor=%.7f: " % (start, variant, FLOORS[fi])
        cells = []
        for j, r in zip(ids, rr):
            if r is None:
                cells.append("%s: BRAK" % j)
                continue
            c = ("%s: %s%s t=%.0f g[%.3f,%.3f] gmax-gceil=%+.3f dev=%.3f "
                 "Ndn%d/%d Nup%d/%d"
                 % (j.split("_")[-1], r["status"],
                    "-" + r["nuc_dir"] if r.get("nuc_dir") else "",
                    r["t_end"], r["gmin"], r["gmax"], r["gmax_vs_gceil"],
                    r["dev_const"], r["n_dn_end"], r["n_seed_dn"],
                    r["n_up_end"], r["n_seed_up"]))
            if r["status"] == "NUCLEATION":
                c += " t0=%.0f Ndet=%d" % (r["t_nuc"], r["n_det"])
            cells.append(c)
        print(line)
        for c in cells:
            print("    " + c)
        if any(r is None for r in rr):
            has_breakdown_or_inc = True
            continue
        sts = [r["status"] for r in rr]
        if all(s == "NUCLEATION" for s in sts):
            dirs = set(r["nuc_dir"] for r in rr)
            rr2 = [load_res(j + "_dt2") for j in ids]
            if len(dirs) == 1 and all(
                    r2 is not None and r2["status"] == "NUCLEATION"
                    and r2["nuc_dir"] == rr[0]["nuc_dir"] for r2 in rr2):
                dets = [r["n_det"] for r in rr] + [r2["n_det"]
                                                  for r2 in rr2]
                okpm1 = max(dets) - min(dets) <= 1
                print("    dt/2: NUCLEATION-%s x2; N_det %s -> zgodnosc "
                      "+-1: %s" % (rr[0]["nuc_dir"], dets,
                                   "TAK" if okpm1 else "NIE"))
                if okpm1:
                    nuc_convergent.append((start, variant, fi,
                                           rr[0]["nuc_dir"], dets))
                else:
                    has_breakdown_or_inc = True
            else:
                print("    dt/2/kierunek: %s -> nukleacja NIEZBIEZNA"
                      % [None if r2 is None else
                         (r2["status"], r2.get("nuc_dir")) for r2 in rr2])
                has_breakdown_or_inc = True
            all_homog = False
        elif all(s == "STATIONARY" for s in sts):
            devs = [r["dev_const"] for r in rr]
            dgrid = common_subgrid_dev(ids[0], ids[1], jobs[ids[0]],
                                       jobs[ids[1]])
            nonconst = min(devs) >= 0.05
            conv = dgrid <= 5e-3
            print("    STACJONARNE: ||g-const||_inf = %s (>=0.05: %s); "
                  "||g_c-g_f||_inf(podsiatka) = %.2e (<=5e-3: %s); "
                  "sanity g_max vs g_ceil: %s"
                  % (["%.4f" % d for d in devs],
                     "TAK" if nonconst else "NIE", dgrid,
                     "TAK" if conv else "NIE",
                     ["%+.4f" % r["gmax_vs_gceil"] for r in rr]))
            if nonconst and conv:
                static_candidates.append((start, variant, fi))
                all_homog = False
            elif nonconst and not conv:
                has_breakdown_or_inc = True
                all_homog = False
        else:
            if any(s in ("BREAKDOWN", "INCOMPLETE") for s in sts):
                has_breakdown_or_inc = True
            if any(s == "NUCLEATION" for s in sts):
                print("    nukleacja tylko na czesci siatek -> NIEZBIEZNA")
                has_breakdown_or_inc = True
            all_homog = False
    print("\nKLASYFIKACJA Q (litera; rozjazd PRIMARY<->C-BAR deskryptywnie"
          " powyzej):")
    if nuc_convergent:
        print("  Q-PASS-NUCLEATION: nukleacja ZBIEZNA (obie siatki + dt/2,"
              " +-1) w: %s"
              % ["%s/%s floor#%d kier=%s N_det=%s" % (s, v, fi + 1, d, n)
                 for s, v, fi, d, n in nuc_convergent])
    elif static_candidates:
        print("  Q-PASS-STATIC: stan stacjonarny NIESTALY zbiezny w: %s"
              % ["%s/%s floor#%d" % (s, v, fi + 1)
                 for s, v, fi in static_candidates])
    elif all_homog and not has_breakdown_or_inc:
        print("  Q-FAIL: wszystko relaksuje do stanu jednorodnego")
    else:
        print("  Q-INCONCLUSIVE: zalamanie nie-nukleacyjne / niezbieznosc /"
              " INCOMPLETE (NIE pozytyw -- litera LOCKa)")
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
    elif mode == "job":
        run_job(sys.argv[2], resume="--resume" in sys.argv)
    elif mode == "jobs":
        for j in sys.argv[2:]:
            if j == "--resume":
                continue
            run_job(j, resume="--resume" in sys.argv)
    elif mode == "verdict":
        verdict()
    else:
        raise SystemExit("nieznany tryb: %s" % mode)
