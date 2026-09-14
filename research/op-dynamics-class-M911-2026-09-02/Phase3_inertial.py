#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
op-dynamics-class-M911 (Phase 3) -- wariant B: dynamika 2. rzedu
BEZWLADNA z lapse M9.1'' dla pary (w, V_M9.1'', K=psi^4).

LOCK sec. 1-2; MD sec. 1,3: B(psi)=w^2 K=psi^6/(4-3psi)^2 (czlon czasowy
akcji z |g^tt|=psi/(4-3psi), CYTAT metryki w MD poprzednika sec. 1);
EOM: B psidd + 1/2 B' psid^2 = -dE/dpsi (wspolny dyskretny gradient
M911_common); H = int[1/2 B v^2 + 1/2 K|grad|^2 + U] zachowane (eta=0).
Integrator RK4 na (u,v), dt=0.0025 (dt/2=0.00125), t_max=100, start v=0.
Gate energii: eps_H = |H(t)-H(0)|/max(|H_rel(0)|,0.01) <= 0.02
(przekroczenie = flaga INTEGRATOR-DRIFT, blokuje PASS).
Pre-rejestracja: psi_min < psi_stiff=0.12 -> flaga STIFF (sztywnosc
B->0); zalamanie po fladze = INCONCLUSIVE-STIFF (rozstrzyga dt/2).
Okno oceny [t_max-20, t_max]; kryteria na skalarach per siatka
(bez porownan polowych -- dekoherencja fazowa, pre-rejestrowane).

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
from M911_common import (Grid3D, Bfun, Bprime, Ueff, detect, banner,  # noqa: E402
                         STARTS, PSI_BAND, PSI_LOW, PSI_STIFF,
                         PSI_LIMIT, N_PAIR)

DT_MAIN = 0.0025
T_MAX = 100.0
EPS_H_TOL = 0.02
DEV_STRUCT = 0.05
OCCUP_MIN = 0.80
WIN = 20.0
RESDIR = BASE + "Phase3_results/"
t0_wall = time.time()


def stamp(msg):
    print("[t=%7.1fs] %s" % (time.time() - t0_wall, msg), flush=True)


class InertialFlow:
    def __init__(self, N, L):
        self.g3 = Grid3D(N, L)
        self.h3 = self.g3.h ** 3

    def acc(self, u, v):
        return (-self.g3.dEdpsi(u) - 0.5 * Bprime(u) * v * v) / Bfun(u)

    def step(self, u, v, dt):
        k1u = v
        k1v = self.acc(u, v)
        u2, v2 = u + 0.5 * dt * k1u, v + 0.5 * dt * k1v
        k2u = v2
        k2v = self.acc(u2, v2)
        u3, v3 = u + 0.5 * dt * k2u, v + 0.5 * dt * k2v
        k3u = v3
        k3v = self.acc(u3, v3)
        u4, v4 = u + dt * k3u, v + dt * k3v
        k4u = v4
        k4v = self.acc(u4, v4)
        un = u + dt / 6.0 * (k1u + 2 * k2u + 2 * k3u + k4u)
        vn = v + dt / 6.0 * (k1v + 2 * k2v + 2 * k3v + k4v)
        return un, vn

    def H(self, u, v):
        return (self.g3.energy(u)
                + 0.5 * float(np.sum(Bfun(u) * v * v)) * self.h3)

    def Hvac(self, u):
        return float(Ueff(1.0)) * self.g3.L ** 3


def run_flow(flow, psi0, dt, label, tmax=T_MAX, ckpt=None, resume=False):
    np.seterr(over='ignore', invalid='ignore')
    spu = int(round(1.0 / dt))
    nsteps = int(round(tmax / dt))
    k0 = 0
    if resume and ckpt and os.path.exists(ckpt):
        ck = np.load(ckpt, allow_pickle=True)
        u, v = np.array(ck["u"]), np.array(ck["v"])
        k0 = int(ck["k"])
        series = json.loads(str(ck["series"]))
        n0_dn, n0_up = int(ck["n0_dn"]), int(ck["n0_up"])
        st_dn, st_up = int(ck["st_dn"]), int(ck["st_up"])
        H0, Hden = float(ck["H0"]), float(ck["Hden"])
        stiff_t = None if float(ck["stiff_t"]) < 0 else float(ck["stiff_t"])
        stamp("  [%s] RESUME t=%.1f" % (label, k0 * dt))
    else:
        u = psi0.copy()
        v = np.zeros_like(u)
        n0_dn, s0d, n0_up, s0u = detect(u)
        H0 = flow.H(u, v)
        Hden = max(abs(H0 - flow.Hvac(u)), 0.01)
        series = [dict(t=0.0, n_dn=n0_dn, n_up=n0_up,
                       pmin=float(u.min()), pmax=float(u.max()),
                       dev=0.5 * float(u.max() - u.min()),
                       H=H0, epsH=0.0, vmax=0.0)]
        st_dn = st_up = 0
        stiff_t = None
        stamp("  [%s] start: N_seed dn=%d up=%d psi[%.4f,%.4f] H=%.6f"
              " H_rel=%.6f" % (label, n0_dn, n0_up, series[0]["pmin"],
                               series[0]["pmax"], H0, H0 - flow.Hvac(u)))
    status, t_end, t_nuc, n_det, ndir = "TMAX", tmax, None, None, None
    drift = any(s["epsH"] > EPS_H_TOL for s in series)
    for k in range(k0 + 1, nsteps + 1):
        up_, vp_ = u, v
        try:
            u, v = flow.step(u, v, dt)
            fin = bool(np.all(np.isfinite(u)) and np.all(np.isfinite(v)))
        except (ValueError, FloatingPointError, OverflowError):
            fin = False
        if not fin:
            status, t_end, u, v = "BREAKDOWN", k * dt, up_, vp_
            stamp("  [%s] ZALAMANIE t=%.4f (stiff_t=%s)"
                  % (label, t_end, stiff_t))
            break
        mx, mn = float(u.max()), float(u.min())
        if mx >= PSI_BAND:
            status, t_end = "BREAKDOWN-BOUNDARY", k * dt
            stamp("  [%s] PAS t=%.4f max=%.9f" % (label, t_end, mx))
            break
        if mn <= PSI_LOW:
            status, t_end = "BREAKDOWN-BOUNDARY-LOWER", k * dt
            stamp("  [%s] DOLNY t=%.4f min=%.3g" % (label, t_end, mn))
            break
        if mn < PSI_STIFF and stiff_t is None:
            stiff_t = k * dt
            stamp("  [%s] FLAGA STIFF t=%.4f psi_min=%.4f < %.2f"
                  % (label, stiff_t, mn, PSI_STIFF))
        if k % spu == 0:
            t = k * dt
            H = flow.H(u, v)
            epsH = abs(H - H0) / Hden
            n_dn, sd, n_up, su = detect(u)
            series.append(dict(t=t, n_dn=n_dn, n_up=n_up, pmin=mn,
                               pmax=mx, dev=0.5 * (mx - mn), H=H,
                               epsH=epsH, vmax=float(np.max(np.abs(v)))))
            if epsH > EPS_H_TOL and not drift:
                drift = True
                stamp("  [%s] INTEGRATOR-DRIFT t=%.1f epsH=%.3g > %.2f"
                      % (label, t, epsH, EPS_H_TOL))
            st_dn = st_dn + 1 if n_dn > n0_dn else 0
            st_up = st_up + 1 if n_up > n0_up else 0
            fired = "DN" if st_dn >= 11 else ("UP" if st_up >= 11 else None)
            if fired:
                t_nuc = t - 10.0
                n_det = series[-11]["n_dn" if fired == "DN" else "n_up"]
                status, t_end, ndir = "NUCLEATION", t, fired
                stamp("  [%s] NUKLEACJA %s t0=%.1f N=%d"
                      % (label, fired, t_nuc, n_det))
                break
            if ckpt and (k % (10 * spu) == 0):
                np.savez_compressed(
                    ckpt, u=u, v=v, k=k, series=json.dumps(series),
                    n0_dn=n0_dn, n0_up=n0_up, st_dn=st_dn, st_up=st_up,
                    H0=H0, Hden=Hden,
                    stiff_t=-1.0 if stiff_t is None else stiff_t)
    last = series[-1]
    win = [s for s in series if s["t"] >= T_MAX - WIN]
    win_min_dev = min((s["dev"] for s in win), default=float('nan'))
    win_max_dev = max((s["dev"] for s in win), default=float('nan'))
    occup = (sum(1 for s in win if s["n_dn"] + s["n_up"] >= 1)
             / max(len(win), 1))
    eps_max = max(s["epsH"] for s in series)
    res = dict(label=label, status=status, t_end=t_end, t_nuc=t_nuc,
               nuc_dir=ndir, n_det=n_det, n_seed_dn=n0_dn,
               n_seed_up=n0_up, pmin=last["pmin"], pmax=last["pmax"],
               pmax_vs_43=last["pmax"] - PSI_LIMIT, dev=last["dev"],
               win_min_dev=win_min_dev, win_max_dev=win_max_dev,
               win_occup=occup, epsH_max=eps_max, drift=drift,
               stiff_t=stiff_t, n_dn_end=last["n_dn"],
               n_up_end=last["n_up"], series=series)
    stamp("  [%s] KONIEC: %s%s t=%.1f psi[%.4f,%.4f] dev=%.6f "
          "okno[dev min/max]=[%.6f,%.6f] occup=%.2f epsH_max=%.2e "
          "drift=%s stiff_t=%s Ndn=%d/%d Nup=%d/%d"
          % (label, status, "-" + ndir if ndir else "", t_end,
             last["pmin"], last["pmax"], last["dev"], win_min_dev,
             win_max_dev, occup, eps_max, drift, stiff_t,
             last["n_dn"], n0_dn, last["n_up"], n0_up))
    return res, u, v


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
    flow = InertialFlow(spec["N"], L)
    banner("; wariant=B(inercja+lapse) job=%s dt=%g" % (jid, dt))
    res, ufin, vfin = run_flow(flow, psi0, dt, jid + " dt=%.4g" % dt,
                               ckpt=RESDIR + jid + "_ckpt.npz",
                               resume=resume)
    res.update(job=jid, dt=dt, N=spec["N"], start=spec["start"])
    with open(RESDIR + jid + ".json", "w") as f:
        json.dump(res, f)
    np.savez_compressed(RESDIR + jid + ".npz", psi=ufin, v=vfin)
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


def osc_ok(r):
    return (r is not None and r["status"] == "TMAX"
            and r["win_min_dev"] >= DEV_STRUCT
            and r["win_occup"] >= OCCUP_MIN and not r["drift"])


def is_event(r):
    return r is not None and (
        r["status"] != "TMAX" or r.get("stiff_t") is not None or osc_ok(r))


def dt2_needed():
    need = []
    for s, ids in PAIRS:
        rr = [load_res(j) for j in ids]
        if any(is_event(r) for r in rr) or all(osc_ok(r) for r in rr):
            need += [j + "_dt2" for j in ids]
    return need


def verdict():
    out = []

    def em(x=""):
        print(x, flush=True)
        out.append(x)

    em("=" * 78)
    em("WERDYKT Q-INER (wariant B: inercja+lapse M9.1''; litera LOCKa)")
    em("=" * 78)
    nuc_conv, osc_conv, disperse, inc = [], [], [], False
    for s, ids in PAIRS:
        rr = [load_res(j) for j in ids]
        em("%s:" % s)
        for j, r in zip(ids, rr):
            if r is None:
                em("    %s: BRAK" % j)
                continue
            em("    %s: %s%s t=%.1f psi[%.4f,%.4f] okno dev[%.6f,%.6f]"
               " occup=%.2f epsH=%.2e drift=%s stiff=%s Ndn%d/%d Nup%d/%d"
               % (j, r["status"],
                  "-" + r["nuc_dir"] if r.get("nuc_dir") else "",
                  r["t_end"], r["pmin"], r["pmax"], r["win_min_dev"],
                  r["win_max_dev"], r["win_occup"], r["epsH_max"],
                  r["drift"], r["stiff_t"], r["n_dn_end"],
                  r["n_seed_dn"], r["n_up_end"], r["n_seed_up"]))
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
                em("    -> nukleacja NIEZBIEZNA")
        elif all(osc_ok(r) for r in rr):
            dt2ok = all(r2 is not None and osc_ok(r2) for r2 in rr2)
            em("    -> OSCYLON-kandydat obie siatki; dt/2: %s"
               % ("TAK" if dt2ok else "NIE/BRAK"))
            if dt2ok:
                osc_conv.append((s, [r["win_min_dev"] for r in rr]))
            else:
                inc = True
        elif all(x == "TMAX" and r["win_max_dev"] < DEV_STRUCT
                 and not r["drift"] for x, r in zip(sts, rr)):
            disperse.append((s, [r["win_max_dev"] for r in rr]))
            em("    -> DYSPERSJA (okno max dev %s < 0.05)"
               % ["%.2e" % r["win_max_dev"] for r in rr])
        else:
            inc = True
            em("    -> wklad INCONCLUSIVE (statusy %s, drift/stiff: %s)"
               % (sts, [(r["drift"], r["stiff_t"]) for r in rr]))
    em()
    em("DESKRYPTYWNIE (obowiazkowe): psi_min(t) dipa (strefa spinodalna)"
       " i N_obj(t):")
    for j in ["dip_N32", "dip_N48"]:
        r = load_res(j)
        if r:
            ser = r["series"]
            pts = [0, 1, 2, 3, 5, 10, 20, 40, 60, 80, 100]
            em("    %s: " % j + " ".join(
                "t=%g:%.3f" % (ser[i]["t"], ser[i]["pmin"])
                for i in range(len(ser)) if ser[i]["t"] in pts))
            em("      N_dn(t): " + " ".join(
                "%g:%d" % (ser[i]["t"], ser[i]["n_dn"])
                for i in range(len(ser)) if ser[i]["t"] in pts))
    em()
    em("KLASYFIKACJA Q-INER (litera):")
    if nuc_conv:
        q = "Q-INER-PASS-NUCLEATION"
        em("  %s: %s" % (q, nuc_conv))
    elif osc_conv:
        q = "Q-INER-PASS-OSCILLON"
        em("  %s: %s" % (q, osc_conv))
    elif len(disperse) == len(PAIRS) and not inc:
        q = "Q-INER-FAIL-DISPERSE"
        em("  %s: pole dyspersuje we wszystkich parach" % q)
    else:
        q = "Q-INER-INCONCLUSIVE"
        em("  %s (NIE pozytyw; kategorie STIFF/DRIFT/BOUNDARY powyzej)" % q)
    arrs = {}
    for jid in list(job_registry().keys()):
        for suff in ("", "_dt2"):
            p = RESDIR + jid + suff + ".npz"
            if os.path.exists(p):
                arrs[jid + suff] = np.load(p)["psi"]
    np.savez_compressed(BASE + "Phase3_final_states.npz", **arrs)
    em()
    em("WERDYKT: %s ; stany: Phase3_final_states.npz (%d)" % (q, len(arrs)))
    with open(BASE + "Phase3_output.txt", "w") as f:
        f.write("\n".join(out) + "\n")
    print("zapisano:", BASE + "Phase3_output.txt")


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
