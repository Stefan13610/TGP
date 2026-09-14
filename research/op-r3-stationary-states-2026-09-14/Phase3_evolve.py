#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
op-r3-stationary-states (Phase 3) -- Q-E RACHUNEK CENTRALNY:
ewolucje 2. rzedu; detektor oscylonu FROZEN (LOCK sec. 1;
MD sec. 5-6).

Starty (LOCK, deterministyczne, pi0=0):
 (i)  gauss: a in {-0.3,-0.15,+0.15,+0.25} x sigma in {3,6}  [8]
 (ii) quasi-R3: a in {-0.2,+0.2}, sigma_w=15                  [2]
 (iii) kontrola: prozna                                        [1]
Siatki h in {0.05, 0.025}; R=200; dt=0.005; t_max=1000;
sponge ON (MD sec. 4); dt/2 (0.0025, h=0.05) przy zdarzeniach.

Detektor (FROZEN): E_ref=E_core(50); t_hold = max t: E_core>=0.5 E_ref
na calym [50,t]; KANDYDAT: (t_hold-50)>=100 T0 oraz >=50 przejsc
psi(0,t) przez 1 na [50,t_hold]; tau = t_hold-50 (cenzura przy
t_hold=t_max); RADIATED: nie-kandydat i E_core(t_max)<=0.05 E_ref;
pas graniczny => BREAKDOWN-BOUNDARY; niefinitycznosc => BREAKDOWN.

Uzycie:
  python Phase3_evolve.py list
  python Phase3_evolve.py job <id> [--resume]
  python Phase3_evolve.py batch <id1> <id2> ... [--resume]
  python Phase3_evolve.py dt2needed
  python Phase3_evolve.py verdict

REJESTR WEJSC [INPUT]: K_geo=gamma=c0=1; rodziny a,sigma [LOCK];
R=200, sponge smootherstep [160,200] gamma0=1.0 [INPUT-MD];
h {0.05,0.025}; dt=0.005; progi detektora [LOCK]; t_transient=50;
t_max=1000; dt_out=0.1; profil co 50; brak seeda.
"""
import json
import os
import sys
import time
import numpy as np

BASE = ("C:/Users/Mateusz/Documents/ObsydnianMain/TGP/TGP_v1/research/"
        "op-r3-stationary-states-2026-09-14/")
sys.path.insert(0, BASE)
import engine_core as ec

RESDIR = BASE + "Phase3_results/"
DT = 0.005
T0 = 2.0*np.pi
T_MAX = 1000.0
T_TRANS = 50.0
DT_OUT = 0.1
PROF_DT = 50.0
HOLD_MIN = 100.0*T0          # 628.3185...
NCROSS_MIN = 50
RAD_FRac = 0.05
R_BOX = 200.0

t0_wall = time.time()


def stamp(msg):
    print("[t=%7.1fs] %s" % (time.time() - t0_wall, msg), flush=True)


# ------------------------------------------------------- rejestr biegow
def job_registry():
    jobs = {}
    for a in (-0.30, -0.15, 0.15, 0.25):
        for s in (3.0, 6.0):
            jobs["g_a%+0.2f_s%g" % (a, s)] = dict(kind="gauss", a=a,
                                                  sigma=s)
    for a in (-0.20, 0.20):
        jobs["qR3_a%+0.2f" % a] = dict(kind="quasiR3", a=a)
    jobs["vac"] = dict(kind="vac")
    return jobs


STARTS = list(job_registry().keys())
H_TAG = {0.05: "h05", 0.025: "h025"}


def run_ids():
    ids = []
    for s in STARTS:
        for h in (0.05, 0.025):
            ids.append("%s__%s" % (s, H_TAG[h]))
    return ids


def parse_id(rid):
    dt2 = rid.endswith("_dt2")
    if dt2:
        rid = rid[:-4]
    sname, htag = rid.rsplit("__", 1)
    h = 0.05 if htag == "h05" else 0.025
    return sname, h, (DT/2 if dt2 else DT)


def build_start(eng, spec):
    if spec["kind"] == "gauss":
        return ec.start_gauss(eng, spec["a"], spec["sigma"])
    if spec["kind"] == "quasiR3":
        return ec.start_quasiR3(eng, spec["a"], 15.0)
    return ec.start_vacuum(eng)


# ------------------------------------------------------------- bieg
def run_job(rid, resume=False):
    os.makedirs(RESDIR, exist_ok=True)
    sname, h, dt = parse_id(rid)
    spec = job_registry()[sname]
    eng = ec.Engine(h, R_BOX, sponge=True)
    ckpt = RESDIR + rid + "_ckpt.npz"
    nsteps = int(round(T_MAX/dt))
    ns_out = int(round(DT_OUT/dt))
    ns_prof = int(round(PROF_DT/dt))
    ns_ck = int(round(100.0/dt))
    nblock = int(round(PROF_DT/DT_OUT))     # probki na blok obwiedni

    if resume and os.path.exists(ckpt):
        ck = np.load(ckpt)
        g = np.array(ck["g"])
        pi = np.array(ck["pi"])
        k0 = int(ck["k"])
        ts = list(ck["ts"])
        pc = list(ck["pc"])
        Ec = list(ck["Ec"])
        profs = list(ck["profs"])
        prof_ts = list(ck["prof_ts"])
        Asum = np.array(ck["Asum"])
        Acnt = int(ck["Acnt"])
        Ablocks = list(ck["Ablocks"])
        stamp("[%s] RESUME t=%.1f" % (rid, k0*dt))
    else:
        g = build_start(eng, spec)
        pi = np.zeros(eng.N)
        k0 = 0
        ts, pc, Ec = [0.0], [float(g[0])], [eng.energy(g, pi, 80.0)]
        profs, prof_ts = [g.copy()], [0.0]
        Asum = np.zeros(eng.N)
        Acnt = 0
        Ablocks = []
        stamp("[%s] START h=%g dt=%g N=%d E_core(0)=%.6e psi0(0)=%.4f"
              % (rid, h, dt, eng.N, Ec[0], pc[0]))

    status, t_end = "OK", T_MAX
    for k in range(k0 + 1, nsteps + 1):
        try:
            g, pi = eng.step(g, pi, dt)
        except (ec.NonConvergence, FloatingPointError):
            status, t_end = "BREAKDOWN", k*dt
            break
        st = eng.band_status(g)
        if st:
            status, t_end = st, k*dt
            break
        if k % ns_out == 0:
            t = k*dt
            ts.append(t)
            pc.append(float(g[0]))
            Ec.append(eng.energy(g, pi, 80.0))
            Asum += np.abs(g - 1.0)
            Acnt += 1
            if Acnt == nblock:
                Ablocks.append((Asum/Acnt).copy())
                Asum[:] = 0.0
                Acnt = 0
        if k % ns_prof == 0:
            profs.append(g.copy())
            prof_ts.append(k*dt)
        if k % ns_ck == 0:
            np.savez_compressed(ckpt, g=g, pi=pi, k=k, ts=ts, pc=pc,
                                Ec=Ec, profs=profs, prof_ts=prof_ts,
                                Asum=Asum, Acnt=Acnt, Ablocks=Ablocks)
            stamp("[%s] t=%.0f E_core=%.4e psi(0)=%.4f "
                  "psi in [%.4f,%.4f]"
                  % (rid, k*dt, Ec[-1], pc[-1], float(np.min(g)),
                     float(np.max(g))))

    res = analyze(rid, sname, h, dt, status, t_end,
                  np.array(ts), np.array(pc), np.array(Ec))
    np.savez_compressed(RESDIR + rid + ".npz",
                        t=np.array(ts), psi_c=np.array(pc),
                        Ecore=np.array(Ec),
                        profs=np.array(profs),
                        prof_ts=np.array(prof_ts),
                        Ablocks=np.array(Ablocks),
                        g_final=g, pi_final=pi,
                        meta=np.array([h, dt, T_MAX]))
    with open(RESDIR + rid + ".json", "w") as f:
        json.dump(res, f, indent=1)
    if os.path.exists(ckpt):
        os.remove(ckpt)
    stamp("[%s] KONIEC: %s tau=%s kandydat=%s przejscia=%d "
          "E_ref=%.4e E_end=%.4e; zapisano json/npz"
          % (rid, res["status"], res["tau_str"], res["candidate"],
             res["ncross"], res["E_ref"], res["E_end"]))
    return res


# ------------------------------------------------ detektor (FROZEN)
def analyze(rid, sname, h, dt, status, t_end, ts, pc, Ec):
    i50 = int(np.searchsorted(ts, T_TRANS))
    res = dict(rid=rid, start=sname, h=h, dt=dt, status=status,
               t_end=t_end)
    if len(ts) <= i50 or status in ("BREAKDOWN",):
        res.update(E_ref=float("nan"), E_end=float("nan"), tau=0.0,
                   tau_str="n/a", censored=False, candidate=False,
                   radiated=False, ncross=0, t_hold=0.0,
                   omega=None, harmonics=[])
        return res
    E_ref = float(Ec[i50])
    below = Ec[i50:] < 0.5*E_ref
    if not np.any(below):
        t_hold = float(ts[-1])
    else:
        jf = int(np.argmax(below))
        t_hold = float(ts[i50 + jf - 1]) if jf > 0 else T_TRANS
    censored = (status == "OK") and (t_hold >= float(ts[-1]) - 1e-9)
    tau = t_hold - T_TRANS
    mwin = (ts >= T_TRANS) & (ts <= t_hold + 1e-9)
    sgn = np.sign(pc[mwin] - 1.0)
    sgn = sgn[sgn != 0]
    ncross = int(np.sum(sgn[1:] != sgn[:-1]))
    candidate = (tau >= HOLD_MIN) and (ncross >= NCROSS_MIN)
    radiated = (not candidate) and (status == "OK") \
        and (float(Ec[-1]) <= RAD_FRac*E_ref) and E_ref > 0
    omega, harm = None, []
    if candidate:
        w0 = max(300.0, t_hold - 600.0)
        m = (ts >= w0) & (ts <= t_hold + 1e-9)
        x = pc[m] - float(np.mean(pc[m]))
        wnd = np.hanning(len(x))
        sp = np.abs(np.fft.rfft(x*wnd))
        dw = 2.0*np.pi/(len(x)*DT_OUT)
        j = 1 + int(np.argmax(sp[1:-1]))
        lm, l0, lp = np.log(sp[j-1] + 1e-300), np.log(sp[j] + 1e-300), \
            np.log(sp[j+1] + 1e-300)
        dlt = 0.5*(lm - lp)/(lm - 2*l0 + lp)
        omega = float((j + dlt)*dw)
        pk = sp[j]**2
        jharm = []
        for jj in range(2, len(sp) - 2):
            if (sp[jj] > sp[jj-1] and sp[jj] > sp[jj+1]
                    and sp[jj]**2 >= 0.01*pk and abs(jj - j) > 2):
                jharm.append((float(jj*dw), float(sp[jj]**2/pk)))
        harm = sorted(jharm, key=lambda z: -z[1])[:6]
    res.update(E_ref=E_ref, E_end=float(Ec[-1]), tau=tau,
               tau_str=(">=%.0f" % tau) if censored else "%.1f" % tau,
               censored=censored, candidate=bool(candidate),
               radiated=bool(radiated), ncross=ncross, t_hold=t_hold,
               omega=omega, harmonics=harm)
    return res


# --------------------------------------------------------- dt2 / werdykt
def load_res(rid):
    p = RESDIR + rid + ".json"
    if not os.path.exists(p):
        return None
    with open(p) as f:
        return json.load(f)


def dt2_needed():
    need = []
    for s in STARTS:
        if s == "vac":
            continue
        rr = [load_res("%s__%s" % (s, tg)) for tg in ("h05", "h025")]
        ev = any(r is not None and (r["candidate"]
                                    or r["status"].startswith(
                                        "BREAKDOWN"))
                 for r in rr)
        if ev:
            need.append("%s__h05_dt2" % s)
    return need


def tau_agree(ra, rb):
    """Zgodnosc czasow zycia (MD sec. 5): <=10%; cenzura wg reguly."""
    ca, cb = ra["censored"], rb["censored"]
    ta, tb = ra["tau"], rb["tau"]
    if ca and cb:
        return True, "obie cenzurowane"
    if ca or cb:
        tother = tb if ca else ta
        return tother >= 0.9*(T_MAX - T_TRANS), "cenzura+prog 855"
    rel = abs(ta - tb)/max(ta, tb, 1e-30)
    return rel <= 0.10, "rel=%.3f" % rel


def classify_start(s):
    r1 = load_res("%s__h05" % s)
    r2 = load_res("%s__h025" % s)
    if r1 is None or r2 is None:
        return "MISSING", r1, r2, None, ""
    rd = load_res("%s__h05_dt2" % s)
    if r1["status"].startswith("BREAKDOWN") \
            and r2["status"].startswith("BREAKDOWN"):
        conv = (r1["status"] == r2["status"]
                and rd is not None and rd["status"] == r1["status"])
        return (r1["status"] if conv else "INCONCLUSIVE"), r1, r2, rd, \
            ("zbiezne (h,h/2,dt/2)" if conv else "niezbiezne")
    if r1["candidate"] and r2["candidate"]:
        if rd is None:
            return "INCONCLUSIVE", r1, r2, rd, "brak dt/2"
        okg, ig = tau_agree(r1, r2)
        okd, idd = tau_agree(r1, rd)
        okc = rd["candidate"]
        if okg and okd and okc:
            return "OSCILLON", r1, r2, rd, "siatki %s; dt/2 %s" % (ig,
                                                                   idd)
        return "INCONCLUSIVE", r1, r2, rd, \
            "kandydat, ale zbieznosc: siatki(%s)=%s dt2(%s)=%s" \
            % (ig, okg, idd, okd and okc)
    if r1["radiated"] and r2["radiated"]:
        return "RADIATED", r1, r2, rd, "zbieznie (obie siatki)"
    return "INCONCLUSIVE", r1, r2, rd, \
        "statusy: h05 cand=%s rad=%s %s / h025 cand=%s rad=%s %s" \
        % (r1["candidate"], r1["radiated"], r1["status"],
           r2["candidate"], r2["radiated"], r2["status"])


def verdict():
    out = []

    def em(s=""):
        print(s, flush=True)
        out.append(s)

    em("=" * 78)
    em("WERDYKT Q-E (litera LOCKa sec. 2 Phase 3; detektor FROZEN "
       "MD sec. 5)")
    em("REJESTR [INPUT]: K_geo=gamma=c0=1; R=200; h {0.05,0.025}; "
       "dt=0.005 (dt/2=0.0025 przy zdarzeniach); t_max=1000; sponge "
       "gamma0=1.0 smootherstep [160,200]; E_core r<=80; E_ref=E_core"
       "(50); prog 0.5; 100 T0=%.3f; >=50 przejsc; RADIATED prog "
       "0.05 E_ref; pas 4/3-1e-6 / 1e-6; brak seeda" % HOLD_MIN)
    em("=" * 78)
    em("Tabela biegow (start x siatka):")
    em("  %-16s %-5s %-22s %-9s %-6s %-8s %-9s %s"
       % ("start", "h", "status", "tau", "cand", "ncross", "omega",
          "E_end/E_ref"))
    for s in STARTS:
        for tg in ("h05", "h025", "h05_dt2"):
            rid = "%s__%s" % (s, tg)
            r = load_res(rid)
            if r is None:
                continue
            em("  %-16s %-5s %-22s %-9s %-6s %-8d %-9s %s"
               % (s, tg, r["status"], r["tau_str"],
                  "TAK" if r["candidate"] else "nie", r["ncross"],
                  "%.4f" % r["omega"] if r["omega"] else "-",
                  "%.3e" % (r["E_end"]/r["E_ref"])
                  if r["E_ref"] and np.isfinite(r["E_ref"])
                  and r["E_ref"] > 0 else "-"))
    em("")
    em("Klasyfikacja per start (zbieznosc: h,h/2 oraz dt/2 przy "
       "zdarzeniach):")
    classes = {}
    for s in STARTS:
        if s == "vac":
            continue
        cl, r1, r2, rd, info = classify_start(s)
        classes[s] = cl
        em("  %-16s: %-22s [%s]" % (s, cl, info))
    rv = [load_res("vac__%s" % tg) for tg in ("h05", "h025")]
    em("  kontrola vac: alarmy detektora = %s (wymagane: zero)"
       % ["%s cand=%s" % (r["status"], r["candidate"]) if r else "BRAK"
          for r in rv])
    em("")
    em("Deskryptywnie (obowiazkowo, LOCK): los startow quasi-R3 "
       "(ksztalt sin(r)/r):")
    for s in ("qR3_a-0.20", "qR3_a+0.20"):
        cl = classes.get(s, "MISSING")
        r1 = load_res("%s__h05" % s)
        em("  %s: klasa %s; tau(h05)=%s, przejscia=%s, omega=%s"
           % (s, cl, r1["tau_str"] if r1 else "-",
              r1["ncross"] if r1 else "-",
              ("%.4f" % r1["omega"]) if (r1 and r1["omega"]) else "-"))
    em("")
    vals = [classes[s] for s in classes]
    n_osc = sum(1 for v in vals if v == "OSCILLON")
    n_rad = sum(1 for v in vals if v == "RADIATED")
    if n_osc >= 1:
        qe = "Q-E-PASS"
    elif n_rad == len(vals):
        qe = "Q-E-FAIL"
    else:
        qe = "Q-E-INCONCLUSIVE"
    em("KLASY: OSCILLON=%d RADIATED=%d BREAKDOWN*=%d INCONCLUSIVE=%d"
       % (n_osc, n_rad,
          sum(1 for v in vals if v.startswith("BREAKDOWN")),
          sum(1 for v in vals if v == "INCONCLUSIVE")))
    em("WERDYKT: %s (PASS: >=1 OSCILLON zbiezny; FAIL: wszystkie "
       "RADIATED zbieznie; INCONCLUSIVE: reszta)" % qe)
    em("=" * 78)
    with open(RESDIR + "verdict.json", "w") as f:
        json.dump(dict(qe=qe, classes=classes), f, indent=1)
    with open(BASE + "Phase3_output.txt", "w", encoding="ascii") as f:
        f.write("\n".join(out) + "\n")
    print("zapisano:", BASE + "Phase3_output.txt")


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "list"
    if mode == "list":
        for rid in run_ids():
            print(rid)
    elif mode == "job":
        run_job(sys.argv[2], resume="--resume" in sys.argv)
    elif mode == "batch":
        for rid in sys.argv[2:]:
            if rid.startswith("--"):
                continue
            run_job(rid, resume="--resume" in sys.argv)
    elif mode == "dt2needed":
        for rid in dt2_needed():
            print(rid)
    elif mode == "verdict":
        verdict()
    else:
        raise SystemExit("nieznany tryb: %s" % mode)
