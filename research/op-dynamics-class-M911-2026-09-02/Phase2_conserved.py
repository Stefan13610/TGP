#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
op-dynamics-class-M911 (Phase 2) -- wariant A: dynamika ZACHOWANA
(Cahn-Hilliard) pary (w, V_M9.1'', K=psi^4).

LOCK sec. 1-2; MD sec. 2: dpsi/dt = nabla^2 mu, mu = dE/dpsi (dokladny
gradient dyskretnej E_h, wspolny rdzen M911_common), M=1 FROZEN.
Krok spektralny semi-implicit: dpsi_h = -dt k^2 mu_h/(1+dt A k^4),
A=1.05 max Keff; masa zachowana dokladnie (wiersz k=0). dt=0.01
(dt/2 kontrola), t_max=200, stacjonarnosc ||psidot||<=1e-8 co dt=1.
Klasyfikacje: STATIONARY/TMAX/NUCLEATION/BREAKDOWN-BOUNDARY(-LOWER)/
BREAKDOWN/BREAKDOWN-NUMERIC (wzrost E > 1e-12 miedzy probkami).

Uzycie: list | job <id> [--resume] | jobs <id...> | dt2needed | verdict
"""
import json
import os
import sys
import time
import numpy as np

BASE = ("C:/Users/Mateusz/Documents/ObsydnianMain/TGP/TGP_v1/research/"
        "op-dynamics-class-M911-2026-09-02/")
sys.path.insert(0, BASE)
from M911_common import (Grid3D, Keff, detect, banner, STARTS,  # noqa: E402
                         PSI_BAND, PSI_LOW, PSI_LIMIT, N_PAIR)

DT_MAIN = 0.01
T_MAX = 200.0
STAT_TOL = 1e-8
E_RISE_TOL = 1e-12
DEV_STRUCT = 0.05
SUBGRID_TOL = 5e-3
TAIL = 20.0
RESDIR = BASE + "Phase2_results/"
t0_wall = time.time()


def stamp(msg):
    print("[t=%7.1fs] %s" % (time.time() - t0_wall, msg), flush=True)


class CHFlow:
    def __init__(self, N, L):
        self.g3 = Grid3D(N, L)
        self.k2 = self.g3.ksym
        self.k4 = self.k2 ** 2

    def mu(self, g):
        return self.g3.dEdpsi(g)

    def psidot(self, g):
        muh = np.fft.rfftn(self.mu(g))
        return np.fft.irfftn(-self.k2 * muh, s=g.shape)

    def step(self, g, dt):
        muh = np.fft.rfftn(self.mu(g))
        A = 1.05 * float(np.max(Keff(g)))
        dpsih = -dt * self.k2 * muh / (1.0 + dt * A * self.k4)
        return g + np.fft.irfftn(dpsih, s=g.shape)

    def energy(self, g):
        return self.g3.energy(g)


def run_flow(flow, psi0, dt, label, tmax=T_MAX, ckpt=None, resume=False):
    np.seterr(over='ignore', invalid='ignore')
    spu = int(round(1.0 / dt))
    nsteps = int(round(tmax / dt))
    k0 = 0
    if resume and ckpt and os.path.exists(ckpt):
        ck = np.load(ckpt, allow_pickle=True)
        g = np.array(ck["g"])
        k0 = int(ck["k"])
        series = json.loads(str(ck["series"]))
        n0_dn, n0_up = int(ck["n0_dn"]), int(ck["n0_up"])
        st_dn, st_up = int(ck["st_dn"]), int(ck["st_up"])
        stamp("  [%s] RESUME t=%.1f" % (label, k0 * dt))
    else:
        g = psi0.copy()
        n0_dn, s0d, n0_up, s0u = detect(g)
        series = [dict(t=0.0, n_dn=n0_dn, n_up=n0_up,
                       pmin=float(g.min()), pmax=float(g.max()),
                       dev=0.5 * float(g.max() - g.min()),
                       pnorm=float(np.max(np.abs(flow.psidot(g)))),
                       E=flow.energy(g), mean=float(g.mean()))]
        st_dn = st_up = 0
        stamp("  [%s] start: N_seed dn=%d up=%d, psi[%.4f,%.4f], "
              "mean=%.6f, E=%.6f, ||psidot||=%.2e"
              % (label, n0_dn, n0_up, series[0]["pmin"], series[0]["pmax"],
                 series[0]["mean"], series[0]["E"], series[0]["pnorm"]))
    status, t_end, t_nuc, n_det, ndir = "TMAX", tmax, None, None, None
    E_prev = series[-1]["E"]
    for k in range(k0 + 1, nsteps + 1):
        gp = g
        try:
            g = flow.step(g, dt)
            fin = bool(np.all(np.isfinite(g)))
        except (ValueError, FloatingPointError, OverflowError):
            fin = False
        if not fin:
            status, t_end, g = "BREAKDOWN", k * dt, gp
            stamp("  [%s] ZALAMANIE t=%.3f" % (label, t_end))
            break
        mx, mn = float(g.max()), float(g.min())
        if mx >= PSI_BAND:
            status, t_end = "BREAKDOWN-BOUNDARY", k * dt
            stamp("  [%s] PAS t=%.3f max=%.9f" % (label, t_end, mx))
            break
        if mn <= PSI_LOW:
            status, t_end = "BREAKDOWN-BOUNDARY-LOWER", k * dt
            stamp("  [%s] DOLNY t=%.3f min=%.3g" % (label, t_end, mn))
            break
        if k % spu == 0:
            t = k * dt
            pnorm = float(np.max(np.abs(flow.psidot(g))))
            E = flow.energy(g)
            n_dn, sd, n_up, su = detect(g)
            series.append(dict(t=t, n_dn=n_dn, n_up=n_up, pmin=mn,
                               pmax=mx, dev=0.5 * (mx - mn), pnorm=pnorm,
                               E=E, mean=float(g.mean())))
            if E > E_prev + E_RISE_TOL:
                status, t_end = "BREAKDOWN-NUMERIC", t
                stamp("  [%s] E ROSNIE t=%.1f (dE=%+.3e) -- klasyfikacja"
                      % (label, t, E - E_prev))
                break
            E_prev = E
            st_dn = st_dn + 1 if n_dn > n0_dn else 0
            st_up = st_up + 1 if n_up > n0_up else 0
            fired = "DN" if st_dn >= 11 else ("UP" if st_up >= 11 else None)
            if fired:
                t_nuc = t - 10.0
                n_det = series[-11]["n_dn" if fired == "DN" else "n_up"]
                status, t_end, ndir = "NUCLEATION", t, fired
                stamp("  [%s] NUKLEACJA %s t0=%.1f N=%d>seed"
                      % (label, fired, t_nuc, n_det))
                break
            if pnorm <= STAT_TOL:
                status, t_end = "STATIONARY", t
                stamp("  [%s] STACJONARNOSC t=%.1f ||psidot||=%.2e"
                      % (label, t, pnorm))
                break
            if ckpt and (k % (10 * spu) == 0):
                np.savez_compressed(ckpt, g=g, k=k,
                                    series=json.dumps(series),
                                    n0_dn=n0_dn, n0_up=n0_up,
                                    st_dn=st_dn, st_up=st_up)
    last = series[-1]
    t_dev = next((s["t"] for s in series if s["dev"] < 1e-3), None)
    tail = [s for s in series if s["t"] >= t_end - TAIL]
    occup = (sum(1 for s in tail if s["n_dn"] + s["n_up"] >= 1)
             / max(len(tail), 1))
    res = dict(label=label, status=status, t_end=t_end, t_nuc=t_nuc,
               nuc_dir=ndir, n_det=n_det, n_seed_dn=n0_dn,
               n_seed_up=n0_up, pmin=last["pmin"], pmax=last["pmax"],
               pmax_vs_43=last["pmax"] - PSI_LIMIT, dev=last["dev"],
               mean_end=last["mean"], pnorm_end=last["pnorm"],
               E_end=last["E"], n_dn_end=last["n_dn"],
               n_up_end=last["n_up"], t_dev=t_dev, tail_occup=occup,
               series=series)
    stamp("  [%s] KONIEC: %s%s t=%.1f psi[%.4f,%.4f] dev=%.6f mean=%.6f"
          " Ndn=%d/%d Nup=%d/%d occup(20)=%.2f t_dev=%s"
          % (label, status, "-" + ndir if ndir else "", t_end,
             last["pmin"], last["pmax"], last["dev"], last["mean"],
             last["n_dn"], n0_dn, last["n_up"], n0_up, occup,
             "%.0f" % t_dev if t_dev is not None else ">t_end"))
    return res, g


def job_registry():
    jobs = {}
    for s in ("gen", "dip", "lat"):
        for N in N_PAIR:
            jobs["%s_N%d" % (s, N)] = dict(start=s, N=N)
    return jobs


PAIRS = [(s, ["%s_N%d" % (s, N) for N in N_PAIR])
         for s in ("gen", "dip", "lat")]


def run_job(jid, resume=False):
    os.makedirs(RESDIR, exist_ok=True)
    base_id = jid[:-4] if jid.endswith("_dt2") else jid
    dt = DT_MAIN / (2 if jid.endswith("_dt2") else 1)
    spec = job_registry()[base_id]
    fstart, L = STARTS[spec["start"]]
    psi0 = fstart(spec["N"])
    flow = CHFlow(spec["N"], L)
    banner("; wariant=A(CH) job=%s dt=%g" % (jid, dt))
    res, gfin = run_flow(flow, psi0, dt, jid + " dt=%.4g" % dt,
                         ckpt=RESDIR + jid + "_ckpt.npz", resume=resume)
    res.update(job=jid, dt=dt, N=spec["N"], start=spec["start"])
    with open(RESDIR + jid + ".json", "w") as f:
        json.dump(res, f)
    np.savez_compressed(RESDIR + jid + ".npz", psi=gfin)
    ck = RESDIR + jid + "_ckpt.npz"
    if os.path.exists(ck):
        os.remove(ck)
    stamp("zapisano: %s.json/.npz" % (RESDIR + jid))
    return res


def load_res(jid):
    p = RESDIR + jid + ".json"
    if not os.path.exists(p):
        return None
    with open(p) as f:
        return json.load(f)


def subgrid_dev(jc, jf):
    dc = np.load(RESDIR + jc + ".npz")["psi"]
    df = np.load(RESDIR + jf + ".npz")["psi"]
    sc, sf = dc.shape[0] // 16, df.shape[0] // 16
    return float(np.max(np.abs(dc[::sc, ::sc, ::sc]
                               - df[::sf, ::sf, ::sf])))


def is_event(r):
    return r is not None and (
        r["status"] == "NUCLEATION" or r["status"].startswith("BREAKDOWN")
        or (r["status"] in ("STATIONARY", "TMAX") and r["dev"] >= DEV_STRUCT
            and r["tail_occup"] >= 0.999))


def dt2_needed():
    need = []
    for s, ids in PAIRS:
        if any(is_event(load_res(j)) for j in ids):
            need += [j + "_dt2" for j in ids]
    return need


def verdict():
    out = []

    def em(x=""):
        print(x, flush=True)
        out.append(x)

    em("=" * 78)
    em("WERDYKT Q-CONS (wariant A: Cahn-Hilliard; litera LOCKa sec. 2)")
    em("=" * 78)
    nuc_conv, struct_conv, homog, inc = [], [], [], False
    for s, ids in PAIRS:
        rr = [load_res(j) for j in ids]
        em("%s:" % s)
        for j, r in zip(ids, rr):
            if r is None:
                em("    %s: BRAK" % j)
                continue
            em("    %s: %s%s t=%.1f psi[%.4f,%.4f] dev=%.6f mean=%.6f "
               "Ndn%d/%d Nup%d/%d occup=%.2f t_dev=%s E=%.6f"
               % (j, r["status"],
                  "-" + r["nuc_dir"] if r.get("nuc_dir") else "",
                  r["t_end"], r["pmin"], r["pmax"], r["dev"],
                  r["mean_end"], r["n_dn_end"], r["n_seed_dn"],
                  r["n_up_end"], r["n_seed_up"], r["tail_occup"],
                  r["t_dev"], r["E_end"]))
        if any(r is None for r in rr):
            inc = True
            continue
        sts = [r["status"] for r in rr]
        rr2 = [load_res(j + "_dt2") for j in ids]
        if all(x == "NUCLEATION" for x in sts):
            dirs = set(r["nuc_dir"] for r in rr)
            if len(dirs) == 1 and all(
                    r2 is not None and r2["status"] == "NUCLEATION"
                    and r2["nuc_dir"] == rr[0]["nuc_dir"] for r2 in rr2):
                dets = [r["n_det"] for r in rr] + [r2["n_det"] for r2 in rr2]
                if max(dets) - min(dets) <= 1:
                    nuc_conv.append((s, rr[0]["nuc_dir"], dets))
                    em("    -> NUKLEACJA zbiezna %s N_det=%s"
                       % (rr[0]["nuc_dir"], dets))
                else:
                    inc = True
            else:
                inc = True
                em("    -> nukleacja NIEZBIEZNA (dt/2/kierunek)")
        elif all(x in ("STATIONARY", "TMAX") for x in sts):
            devs = [r["dev"] for r in rr]
            if min(devs) >= DEV_STRUCT:
                dg = subgrid_dev(ids[0], ids[1])
                occ = all(r["tail_occup"] >= 0.999 for r in rr)
                dt2ok = all(r2 is not None and r2["dev"] >= DEV_STRUCT
                            and r2["tail_occup"] >= 0.999 for r2 in rr2)
                em("    -> STRUKTURA-kandydat: dev=%s podsiatka=%.2e "
                   "(<=5e-3:%s) occup:%s dt/2:%s"
                   % (["%.4f" % d for d in devs], dg,
                      "TAK" if dg <= SUBGRID_TOL else "NIE",
                      "TAK" if occ else "NIE",
                      "TAK" if dt2ok else "NIE/BRAK"))
                if dg <= SUBGRID_TOL and occ and dt2ok:
                    struct_conv.append((s, devs, dg))
                else:
                    inc = True
            else:
                homog.append((s, [r["mean_end"] for r in rr],
                              [r["t_dev"] for r in rr]))
                em("    -> JEDNORODNE (dev<0.05): mean=%s t_dev=%s"
                   % (["%.6f" % r["mean_end"] for r in rr],
                      [r["t_dev"] for r in rr]))
        else:
            inc = True
            em("    -> statusy %s: wklad INCONCLUSIVE" % sts)
    em()
    em("DESKRYPTYWNIE (obowiazkowe): t_dev(dev<1e-3) vs poprzednik "
       "(gradient flow: gen 7.0 / bump 16.0 / lat 16.0):")
    for s, ids in PAIRS:
        for j in ids:
            r = load_res(j)
            if r:
                em("    %s: t_dev=%s (status %s, t_end=%.0f)"
                   % (j, r["t_dev"], r["status"], r["t_end"]))
    em()
    em("KLASYFIKACJA Q-CONS (litera):")
    if nuc_conv:
        q = "Q-CONS-PASS-NUCLEATION"
        em("  %s: %s" % (q, nuc_conv))
    elif struct_conv:
        q = "Q-CONS-PASS-STRUCTURE"
        em("  %s: %s" % (q, [(s, d) for s, d, _ in struct_conv]))
    elif len(homog) == len(PAIRS) and not inc:
        q = "Q-CONS-FAIL"
        em("  %s: wszystko dyfunduje do jednorodnej sredniej" % q)
    else:
        q = "Q-CONS-INCONCLUSIVE"
        em("  %s (NIE pozytyw)" % q)
    arrs = {}
    for jid in list(job_registry().keys()):
        for suff in ("", "_dt2"):
            p = RESDIR + jid + suff + ".npz"
            if os.path.exists(p):
                arrs[jid + suff] = np.load(p)["psi"]
    np.savez_compressed(BASE + "Phase2_relaxed_states.npz", **arrs)
    em()
    em("WERDYKT: %s ; stany: Phase2_relaxed_states.npz (%d)"
       % (q, len(arrs)))
    with open(BASE + "Phase2_output.txt", "w") as f:
        f.write("\n".join(out) + "\n")
    print("zapisano:", BASE + "Phase2_output.txt")


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "list"
    if mode == "list":
        for k in job_registry():
            print(k)
    elif mode == "job":
        run_job(sys.argv[2], resume="--resume" in sys.argv)
    elif mode == "jobs":
        for j in sys.argv[2:]:
            if j != "--resume":
                run_job(j, resume="--resume" in sys.argv)
    elif mode == "dt2needed":
        for j in dt2_needed():
            print(j)
    elif mode == "verdict":
        verdict()
    else:
        raise SystemExit("nieznany tryb: %s" % mode)
