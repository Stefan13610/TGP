#!/usr/bin/env python3
# -*- coding: ascii -*-
"""op-matter-induced-creation -- Phase 2: bramka (LOCK sec.3; MD sec.5).
FAIL ktoregokolwiek gate'u => STOP cyklu.

P2a: proznia lam=0, sponge ON, psi=1, 100 T0, obie siatki;
     gate ||psi-1||_inf <= 1e-10 w calym biegu.
P2b: regresje wobec kotwic poprzednika (Phase3_qh1_output.txt):
     (1) lam=0.05 stale, start psi=1, t=300, h=0.05, okno [250,300]:
         psibar(0) = 0.865982 +-1%;
     (2) lam=0.5 stale, start psi=1, h=0.05: COLLAPSE t_end=0.375 +-5%.
P2c: dryf energii przy lam STALYM 0.05 (start gauss a=+0.05 s=3,
     sponge OFF, t=700, obie siatki); estymator SEKULARNY (MD sec.5,
     przejety z correction note 1 poprzednika, zamrozony PRZED biegiem):
     dryf = (100/80)*|<E>_[90T0,100T0] - <E>_[10T0,20T0]|/|<E>_[0,10T0]|
     <= 1e-6 /100T0.
"""
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import engine_core as ec      # noqa: E402
import runlib as rl           # noqa: E402

OUTP = os.path.join(HERE, "Phase2_output.txt")
RES = os.path.join(HERE, "Phase3_results")
T0 = 2.0*np.pi
DT = 0.005


def task_p2a(h):
    t1 = time.time()
    eng = ec.Engine(h, 200.0, sponge=True, lam=0.0)
    g = ec.start_vacuum(eng)
    pi = np.zeros(eng.N)
    nst = int(round(100*T0/DT))
    dev = 0.0
    for n in range(nst):
        g, pi = eng.step(g, pi, DT)
        d = float(np.max(np.abs(g - 1.0)))
        if d > dev:
            dev = d
    return ("p2a", h, dict(dev=dev, ok=bool(dev <= 1e-10),
                           sec=time.time() - t1))


def _vac(h, sponge, lam):
    eng = ec.Engine(h, 200.0, sponge=sponge, lam=lam)
    return ec.start_vacuum(eng), np.zeros(eng.N)


def task_p2b_1():
    t1 = time.time()
    g, pi = _vac(0.05, True, 0.05)
    res = rl.evolve(0.05, DT, 0.05, 0.0, 300.0, g, pi, sponge=True,
                    ramp=False)
    if res["status"] != "OK":
        return ("p2b1", 0.05, dict(status=res["status"],
                                   t_end=res["t_end"], ok=False,
                                   sec=time.time() - t1))
    mp = rl.window_mask(res["prof_t"], 250.0, 300.0)
    psibar = res["prof"][mp].mean(axis=0)
    pb0 = float(psibar[0])
    ok = abs(pb0 - 0.865982) <= 0.01*0.865982
    return ("p2b1", 0.05, dict(status="OK", psibar0=pb0,
                               n_prof=int(np.sum(mp)), ok=bool(ok),
                               sec=time.time() - t1))


def task_p2b_2():
    t1 = time.time()
    g, pi = _vac(0.05, True, 0.5)
    res = rl.evolve(0.05, DT, 0.5, 0.0, 30.0, g, pi, sponge=True,
                    ramp=False, keep_profiles=False)
    te = res["t_end"]
    ok = (res["status"] in ("BREAKDOWN", "BREAKDOWN-BOUNDARY",
                            "BREAKDOWN-BOUNDARY-LOWER")
          and te is not None and abs(te - 0.375) <= 0.05*0.375)
    return ("p2b2", 0.5, dict(status=res["status"], t_end=te,
                              ok=bool(ok), sec=time.time() - t1))


def task_p2c(h):
    t1 = time.time()
    eng = ec.Engine(h, 200.0, sponge=False, lam=0.05)
    g = ec.start_gauss(eng, 0.05, 3.0)
    pi = np.zeros(eng.N)
    nst = int(round(700.0/DT))
    nout = int(round(0.1/DT))
    ts = [0.0]
    Es = [eng.energy(g, pi)]
    bad = None
    for n in range(nst):
        try:
            g, pi = eng.step(g, pi, DT)
        except ec.NonConvergence:
            bad = "NonConvergence@%.3f" % ((n + 1)*DT)
            break
        st = eng.band_status(g)
        if st is not None:
            bad = st + "@%.3f" % ((n + 1)*DT)
            break
        if (n + 1) % nout == 0:
            ts.append((n + 1)*DT)
            Es.append(eng.energy(g, pi))
    if bad is not None:
        return ("p2c", h, dict(bad=bad, ok=False, sec=time.time() - t1))
    ts = np.array(ts)
    Es = np.array(Es)
    E0m = float(np.mean(Es[(ts >= 0.0) & (ts <= 10*T0)]))
    Ee = float(np.mean(Es[(ts >= 10*T0) & (ts <= 20*T0)]))
    El = float(np.mean(Es[(ts >= 90*T0) & (ts <= 100*T0)]))
    drift = (100.0/80.0)*abs(El - Ee)/abs(E0m)
    drift_old = abs(El - E0m)/abs(E0m)
    mf = (ts >= 10*T0) & (ts <= 100*T0)
    A = np.vstack([ts[mf], np.ones(int(np.sum(mf)))]).T
    sol = np.linalg.lstsq(A, Es[mf], rcond=None)[0]
    drift_fit = abs(float(sol[0]))*100*T0/abs(E0m)
    return ("p2c", h, dict(E0=float(Es[0]), Ee=Ee, El=El, E0m=E0m,
                           drift=drift, drift_old=drift_old,
                           drift_fit=drift_fit,
                           mx=float(np.max(np.abs(Es - Es[0]))),
                           ok=bool(drift <= 1e-6), sec=time.time() - t1))


def main():
    rl.ensure_dir(RES)
    out = []

    def w(s=""):
        out.append(str(s))
        with open(OUTP, "w") as fh:
            fh.write("\n".join(out) + "\n")

    BAR = "=" * 78
    t_all = time.time()
    w(BAR)
    w("PHASE 2 -- bramka (FROZEN; MD sec.5). FAIL => STOP")
    w("REJESTR [INPUT]: dt=0.005; R=200; sponge gamma0=1 smootherstep")
    w("  [160,200] (ON: P2a/P2b, OFF: P2c [INPUT-MD]); dt_out=0.1;")
    w("  progi: P2a 1e-10; P2b psibar(0)=0.865982+-1%, t_end=0.375+-5%;")
    w("  P2c dryf sekularny <=1e-6/100T0, okna [10,20]T0 vs [90,100]T0")
    w("  [INPUT-MD: estymator przejety z correction note 1 poprzednika]")
    w(BAR)
    tasks = []
    with ProcessPoolExecutor(max_workers=6) as ex:
        futs = [ex.submit(task_p2a, 0.05), ex.submit(task_p2a, 0.025),
                ex.submit(task_p2b_1), ex.submit(task_p2b_2),
                ex.submit(task_p2c, 0.05), ex.submit(task_p2c, 0.025)]
        for f in futs:
            tasks.append(f.result())
    R = {}
    for kind, key, d in tasks:
        R[(kind, key)] = d

    w("P2a: proznia lam=0, sponge ON, 100 T0 = %.6f" % (100*T0))
    p2a_ok = True
    for h in (0.05, 0.025):
        d = R[("p2a", h)]
        p2a_ok = p2a_ok and d["ok"]
        w("  h=%.3f: max||psi-1||_inf = %.3e (prog 1e-10) %s [%.0f s]"
          % (h, d["dev"], "PASS" if d["ok"] else "FAIL", d["sec"]))
    w("  P2a: " + ("PASS" if p2a_ok else "FAIL"))

    w("P2b-1: regresja lam=0.05 stale, start psi=1, t=300, okno [250,300]")
    d = R[("p2b1", 0.05)]
    w("  status=%s psibar(0)=%s (kotwica 0.865982, +-1%% => [0.857322,"
      " 0.874642]); n_prof=%s"
      % (d["status"], ("%.6f" % d["psibar0"]) if "psibar0" in d else "n/a",
         d.get("n_prof", "n/a")))
    if "psibar0" in d:
        w("  dpsi(0) = %+.6f (kotwica -0.134018; wzgl. odchylenie %.3f%%)"
          % (d["psibar0"] - 1.0,
             100.0*abs(d["psibar0"] - 0.865982)/0.865982))
    w("  P2b-1: %s [%.0f s]" % ("PASS" if d["ok"] else "FAIL", d["sec"]))
    p2b1_ok = d["ok"]

    w("P2b-2: regresja lam=0.5 stale, start psi=1, h=0.05")
    d = R[("p2b2", 0.5)]
    w("  status=%s t_end=%s (kotwica 0.375, +-5%% => [0.35625, 0.39375])"
      % (d["status"], "%.4f" % d["t_end"] if d["t_end"] else "n/a"))
    w("  P2b-2: %s [%.0f s]" % ("PASS" if d["ok"] else "FAIL", d["sec"]))
    p2b2_ok = d["ok"]
    p2b_ok = p2b1_ok and p2b2_ok
    w("  P2b: " + ("PASS" if p2b_ok else "FAIL"))

    w("P2c: energia przy lam STALYM 0.05 (gauss a=+0.05 s=3, sponge OFF,"
      " t=700)")
    p2c_ok = True
    for h in (0.05, 0.025):
        d = R[("p2c", h)]
        p2c_ok = p2c_ok and d["ok"]
        if "bad" in d:
            w("  h=%.3f: PRZERWANE (%s) FAIL" % (h, d["bad"]))
            continue
        w("  h=%.3f: E(0)=%.9e <E>_[10,20]T0=%.9e <E>_[90,100]T0=%.9e"
          % (h, d["E0"], d["Ee"], d["El"]))
        w("          dryf_sekularny=%.3e (prog 1e-6/100T0) %s [%.0f s]"
          % (d["drift"], "PASS" if d["ok"] else "FAIL", d["sec"]))
        w("          deskr.: offset_dC2(stary estym.)=%.3e"
          "  dryf_fitLSQ[10,100]T0=%.3e  max|E-E0|=%.3e"
          % (d["drift_old"], d["drift_fit"], d["mx"]))
    w("  P2c: " + ("PASS" if p2c_ok else "FAIL"))

    w(BAR)
    allok = p2a_ok and p2b_ok and p2c_ok
    w("PHASE 2 PODSUMOWANIE: P2a %s ; P2b %s ; P2c %s => %s [%.0f s]"
      % ("PASS" if p2a_ok else "FAIL", "PASS" if p2b_ok else "FAIL",
         "PASS" if p2c_ok else "FAIL",
         "BRAMKA OTWARTA" if allok else "STOP (FAIL => STOP)",
         time.time() - t_all))
    w(BAR)
    w("PHASE2 DONE")


if __name__ == "__main__":
    main()
    print("PHASE2 DONE")
