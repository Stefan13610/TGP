#!/usr/bin/env python3
# -*- coding: ascii -*-
"""op-collapse-matter-source -- Phase 3 / Q-H2: stabilizacja kolapsu
(LOCK sec. 3 Phase 3 Q-H2; MD sec. 7).

Starty (FROZEN, identyczne definicje jak u poprzednika): qR3 a=-0.20,
qR3 a=+0.20, gauss a=-0.30 s=3, gauss a=+0.15 s=6; pi0=0; zrodlo od
t=0; lam in {0.05, 0.2, 0.5}; t_max=1000; h=0.05; dt=0.005; sponge ON.
Kategorie: COLLAPSE (nadkategoria, podtyp) / RADIATED / STABILIZED /
INCONCLUSIVE-RUN (MD sec. 7). Potwierdzenia (FROZEN): kategoria !=
baseline (COLLAPSE) -> h=0.025 ORAZ dt/2 (0.0025 na h=0.05); kontrola
negatywu: qR3 a=+0.20 przy najnizszym lam dajacym COLLAPSE -> h=0.025.

Uzycie: python Phase3_qh2_stabilize.py {main <lam> | confirm | verdict}
"""
import glob
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, "TGP/TGP_v1/research/op-collapse-matter-source-2026-09-14")
import engine_core as ec

DIR = "TGP/TGP_v1/research/op-collapse-matter-source-2026-09-14"
RES = DIR + "/Phase3_results"
os.makedirs(RES, exist_ok=True)

DT = 0.005
T_MAX = 1000.0
DT_OUT = 0.1
DT_PROF = 50.0
DT_CHK = 100.0
LAMS = [0.05, 0.2, 0.5]
STARTS = ["qR3_a-0.20", "qR3_a+0.20", "g_a-0.30_s3", "g_a+0.15_s6"]
# baseline lam=0 poprzednika (Phase3_output.txt + Phase3_results):
# wszystkie COLLAPSE (BREAKDOWN-BOUNDARY gorny); t_end h05:
BASE_T = {"qR3_a-0.20": 4.74, "qR3_a+0.20": 2.275,
          "g_a-0.30_s3": 6.64, "g_a+0.15_s6": 18.185}


def mk_start(eng, key):
    if key == "qR3_a-0.20":
        return ec.start_quasiR3(eng, -0.20)
    if key == "qR3_a+0.20":
        return ec.start_quasiR3(eng, +0.20)
    if key == "g_a-0.30_s3":
        return ec.start_gauss(eng, -0.30, 3.0)
    if key == "g_a+0.15_s6":
        return ec.start_gauss(eng, +0.15, 6.0)
    raise KeyError(key)


def log(stage, s):
    with open(RES + "/qh2_log_%s.txt" % stage, "a") as fh:
        fh.write(s + "\n")


def run_one(key, lam, h, dt, tag, stage):
    rid = "qh2_%s__lam%g__%s" % (key, lam, tag)
    if os.path.exists(RES + "/" + rid + ".json"):
        with open(RES + "/" + rid + ".json") as fh:
            return json.load(fh)
    eng = ec.Engine(h, 200.0, sponge=True, lam=lam)
    g = mk_start(eng, key)
    pi = np.zeros(eng.N)
    nst = int(round(T_MAX/dt))
    nout = int(round(DT_OUT/dt))
    nprof = int(round(DT_PROF/dt))
    nchk = int(round(DT_CHK/dt))
    m80 = eng.r <= 80.0
    ts = [0.0]
    ecore = [eng.energy(g, pi, 80.0)]
    ecf = [eng.energy(g, pi, 80.0, matter=False)]
    psi0 = [float(g[0])]
    vmax = [0.0]
    pts, profs = [0.0], [g.copy()]
    status, t_end = "OK", None
    t1 = time.time()
    for n in range(nst):
        try:
            g, pi = eng.step(g, pi, dt)
        except ec.NonConvergence:
            status, t_end = "BREAKDOWN", (n + 1)*dt
            break
        st = eng.band_status(g)
        if st is not None:
            status, t_end = st, (n + 1)*dt
            break
        t = (n + 1)*dt
        if (n + 1) % nout == 0:
            ts.append(t)
            psi0.append(float(g[0]))
            ecore.append(eng.energy(g, pi, 80.0))
            ecf.append(eng.energy(g, pi, 80.0, matter=False))
            vmax.append(float(np.max(np.abs((pi/ec.Mfun(g))[m80]))))
        if (n + 1) % nprof == 0:
            pts.append(t)
            profs.append(g.copy())
        if (n + 1) % nchk == 0:
            np.savez(RES + "/checkpoint_" + rid + ".npz", r=eng.r,
                     psi=g, pi=pi, t=t)
    wall = time.time() - t1
    ts = np.array(ts)
    ecore = np.array(ecore)
    ecf = np.array(ecf)
    vmax = np.array(vmax)
    pts = np.array(pts)
    profs = np.array(profs)
    rec = dict(rid=rid, start=key, lam=lam, h=h, dt=dt, tag=tag,
               status=status, t_end=t_end, wall=round(wall, 1))
    # ---- klasyfikacja (MD sec. 7) --------------------------------
    if status != "OK":
        rec["cls"] = "COLLAPSE"
        rec["subtyp"] = status
        rec["tau_str"] = "n/a"
        rec["tau"] = None
        rec["censored"] = False
    else:
        i50 = int(np.argmin(np.abs(ts - 50.0)))
        E_ref = float(ecore[i50])
        E_end = float(ecore[-1])
        rec["E_ref"] = E_ref
        rec["E_end"] = E_end
        rec["E_end_field"] = float(ecf[-1])
        rec["E_ref_field"] = float(ecf[i50])
        if E_ref <= 0.0:
            rec["cls"] = "INCONCLUSIVE-RUN"
            rec["note"] = "E_ref<=0 (progi energetyczne bez sensu)"
            rec["tau"], rec["tau_str"] = None, "n/a"
            rec["censored"] = False
        else:
            m = ts >= 50.0
            sub = ecore[m]
            tt = ts[m]
            below = np.where(sub < 0.5*E_ref)[0]
            if len(below) == 0:
                t_hold = float(tt[-1])
            elif below[0] == 0:
                t_hold = 50.0
            else:
                t_hold = float(tt[below[0] - 1])
            tau = t_hold - 50.0
            cens = t_hold >= T_MAX
            rec["t_hold"] = t_hold
            rec["tau"] = tau
            rec["censored"] = cens
            rec["tau_str"] = (">=950" if cens else "%.1f" % tau)
            if E_end < 0.05*E_ref:
                rec["cls"] = "RADIATED"
            else:
                rec["cls"] = "STABILIZED"
                # podtyp deskryptywny: okno [950,1000]
                mw = pts >= 950.0
                if np.any(mw):
                    psibar = np.mean(profs[mw], axis=0)
                    D = float(np.max(np.abs(psibar[m80] - 1.0)))
                    mv = ts >= 950.0
                    V = float(np.max(vmax[mv]))
                    rec["V"] = V
                    rec["D"] = D
                    rec["subtyp"] = ("osiadly-statycznie"
                                     if V <= 0.01*max(D, 1e-12)
                                     else "oscylujacy")
                    rec["psibar_ctr"] = float(psibar[0])
                    rec["psibar_min"] = float(np.min(psibar))
                    rec["psibar_max"] = float(np.max(psibar))
    np.savez(RES + "/" + rid + ".npz", r=eng.r, ts=ts, psi0=psi0,
             ecore=ecore, ecore_field=ecf, vmax=vmax, prof_ts=pts,
             profs=profs)
    with open(RES + "/" + rid + ".json", "w") as fh:
        json.dump(rec, fh, indent=1, default=str)
    log(stage, "%s: cls=%s%s tau=%s [%.0f s]"
        % (rid, rec["cls"],
           ("" if rec.get("t_end") is None
            else " t_end=%.3f (%s)" % (rec["t_end"],
                                       rec.get("subtyp", "?"))),
           rec.get("tau_str", "-"), wall))
    return rec


def load(rid):
    p = RES + "/" + rid + ".json"
    if os.path.exists(p):
        with open(p) as fh:
            return json.load(fh)
    return None


def tau_agree(a, b):
    """Zgodnosc tau/czasu zdarzenia (MD sec. 7; +-10%, cenzura)."""
    if a["cls"] != b["cls"]:
        return False
    if a["cls"] == "COLLAPSE":
        ta, tb = a["t_end"], b["t_end"]
        return abs(ta - tb)/max(ta, tb) <= 0.10
    ca, cb = a.get("censored", False), b.get("censored", False)
    if ca and cb:
        return True
    if ca != cb:
        tv = b["tau"] if ca else a["tau"]
        return tv is not None and tv >= 855.0
    ta, tb = a.get("tau"), b.get("tau")
    if ta is None or tb is None:
        return True   # np. INCONCLUSIVE-RUN -- kategoria juz zgodna
    if max(ta, tb) == 0.0:
        return True
    return abs(ta - tb)/max(ta, tb) <= 0.10


def stage_main(lam_filter=None):
    for lam in LAMS:
        if lam_filter is not None and lam != lam_filter:
            continue
        for key in STARTS:
            run_one(key, lam, 0.05, DT, "h05", "main")


def stage_confirm(lam_filter=None):
    mains = {}
    for lam in LAMS:
        for key in STARTS:
            rec = load("qh2_%s__lam%g__h05" % (key, lam))
            if rec is None:
                raise RuntimeError("brak main: %s lam=%g" % (key, lam))
            mains[(key, lam)] = rec
    # potwierdzenia kategorii != COLLAPSE: h/2 oraz dt/2
    for (key, lam), rec in sorted(mains.items()):
        if lam_filter is not None and lam != lam_filter:
            continue
        if rec["cls"] != "COLLAPSE":
            run_one(key, lam, 0.025, DT, "h025", "confirm")
            run_one(key, lam, 0.05, DT/2, "h05_dt2", "confirm")
    # kontrola negatywu: qR3_a+0.20, najnizsze lam z COLLAPSE -> h/2
    neg = [lam for lam in sorted(LAMS)
           if mains[("qR3_a+0.20", lam)]["cls"] == "COLLAPSE"]
    if neg and (lam_filter is None or lam_filter == neg[0]):
        run_one("qR3_a+0.20", neg[0], 0.025, DT, "h025", "confirm")
        log("confirm", "kontrola negatywu: qR3_a+0.20 lam=%g" % neg[0])
    elif not neg:
        log("confirm", "kontrola negatywu: BRAK biegu COLLAPSE dla"
            " qR3_a+0.20 (wszystkie lam stabilizuja)")


def stage_verdict():
    OUT = []

    def w(s=""):
        OUT.append(str(s))

    BAR = "=" * 78
    w(BAR)
    w("PHASE 3 / Q-H2 -- stabilizacja kolapsu (LOCK sec. 3-4;"
      " klasyfikatory MD sec. 7)")
    w("REJESTR [INPUT]: starty FROZEN (4 reprezentanci kolapsu"
      " poprzednika, profile")
    w("  identyczne); lam {0.05,0.2,0.5}; t_max=1000; h=0.05"
      " (potw. h=0.025 i dt=0.0025);")
    w("  dt=0.005; sponge ON; E_core r<=80 gestosc PELNA (z U_mat,"
      " prozniowo odjeta);")
    w("  E_ref=E_core(50); progi 0.05/0.5 E_ref; pas 4/3-1e-6 /"
      " 1e-6; tau: t_hold-50;")
    w("  potwierdzenia +-10%; baseline lam=0 = COLLAPSE (wszystkie"
      " 4 starty, sufit)")
    w(BAR)
    w("Tabela biegow (start x lam x konfiguracja):")
    w("  %-13s %-5s %-8s %-13s %-24s %-8s %s"
      % ("start", "lam", "konfig", "klasa", "podtyp/t_end", "tau",
         "E_end/E_ref"))
    allrecs = {}
    for p in sorted(glob.glob(RES + "/qh2_*.json")):
        with open(p) as fh:
            rec = json.load(fh)
        allrecs[rec["rid"]] = rec
    for rid in sorted(allrecs):
        rec = allrecs[rid]
        ev = ("-" if rec.get("t_end") is None
              else "%s@t=%.3f" % (rec.get("subtyp", "?"),
                                  rec["t_end"]))
        if rec["cls"] == "STABILIZED":
            ev = rec.get("subtyp", "?")
        er = ("-" if rec.get("E_ref") in (None, 0)
              or rec.get("E_end") is None
              else "%.3e" % (rec["E_end"]/rec["E_ref"]))
        w("  %-13s %-5g %-8s %-13s %-24s %-8s %s"
          % (rec["start"], rec["lam"], rec["tag"], rec["cls"], ev,
             rec.get("tau_str", "-"), er))
    w("-" * 78)
    w("Zbieznosc i kategoria koncowa par (potwierdzenia MD sec. 7):")
    finals = {}
    neg_pair = None
    neg_lams = [lam for lam in sorted(LAMS)
                if allrecs.get("qh2_qR3_a+0.20__lam%g__h05"
                               % lam, {}).get("cls") == "COLLAPSE"]
    if neg_lams:
        neg_pair = ("qR3_a+0.20", neg_lams[0])
    for lam in LAMS:
        for key in STARTS:
            main = allrecs["qh2_%s__lam%g__h05" % (key, lam)]
            if main["cls"] == "COLLAPSE":
                fin = "COLLAPSE"
                note = "kategoria na h [LOCK]"
                if neg_pair == (key, lam):
                    cf = allrecs.get("qh2_%s__lam%g__h025"
                                     % (key, lam))
                    okc = cf is not None and tau_agree(main, cf)
                    note = ("kontrola negatywu h/2: " +
                            ("zgodna (t_end %.3f vs %.3f)"
                             % (main["t_end"], cf["t_end"]) if okc
                             else "NIEZGODNA"))
                    if not okc:
                        fin = "INCONCLUSIVE-RUN"
            else:
                c1 = allrecs.get("qh2_%s__lam%g__h025" % (key, lam))
                c2 = allrecs.get("qh2_%s__lam%g__h05_dt2"
                                 % (key, lam))
                ok1 = c1 is not None and tau_agree(main, c1)
                ok2 = c2 is not None and tau_agree(main, c2)
                if ok1 and ok2:
                    fin = main["cls"]
                    note = "potwierdzone (h/2 i dt/2)"
                else:
                    fin = "INCONCLUSIVE-RUN"
                    note = ("rozjazd potwierdzen: h025=%s dt2=%s"
                            % (c1["cls"] if c1 else "brak",
                               c2["cls"] if c2 else "brak"))
            finals[(key, lam)] = fin
            base = BASE_T[key]
            w("  %-13s lam=%-5g: %-16s [%s; baseline lam=0:"
              " COLLAPSE t=%.3f]" % (key, lam, fin, note, base))
    w("-" * 78)
    n_stab = sum(1 for v in finals.values() if v == "STABILIZED")
    n_coll = sum(1 for v in finals.values() if v == "COLLAPSE")
    n_rad = sum(1 for v in finals.values() if v == "RADIATED")
    n_inc = sum(1 for v in finals.values()
                if v == "INCONCLUSIVE-RUN")
    w("KLASY KONCOWE: STABILIZED=%d RADIATED=%d COLLAPSE=%d"
      " INCONCLUSIVE-RUN=%d (z 12)" % (n_stab, n_rad, n_coll, n_inc))
    if n_stab >= 1:
        verdict = "Q-H2-PASS"
    elif n_coll == 12:
        verdict = "Q-H2-FAIL"
    else:
        verdict = "Q-H2-INCONCLUSIVE"
    w("WERDYKT Q-H2 (litera LOCKa sec. 4): " + verdict)
    w("  (PASS: >=1 para STABILIZED potwierdzona h/2 i dt/2 przy"
      " baseline COLLAPSE;")
    w("   FAIL: wszystkie 12 par COLLAPSE zbieznie; inaczej"
      " INCONCLUSIVE)")
    w(BAR)
    with open(DIR + "/Phase3_qh2_output.txt", "w") as fh:
        fh.write("\n".join(OUT) + "\n")
    qh1 = {}
    p = RES + "/qh1_summary.json"
    if os.path.exists(p):
        with open(p) as fh:
            qh1 = json.load(fh)
    with open(RES + "/verdict.json", "w") as fh:
        json.dump(dict(QH1=qh1.get("verdict"),
                       QH2=verdict,
                       qh1_final=qh1.get("final"),
                       qh2_final={"%s|%g" % k: v
                                  for k, v in finals.items()},
                       n=dict(STABILIZED=n_stab, RADIATED=n_rad,
                              COLLAPSE=n_coll,
                              INCONCLUSIVE=n_inc)), fh, indent=1)
    print("verdict written:", verdict)


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "main"
    if mode == "main":
        lamf = float(sys.argv[2]) if len(sys.argv) > 2 else None
        stage_main(lamf)
        print("main done", lamf)
    elif mode == "confirm":
        lamf = float(sys.argv[2]) if len(sys.argv) > 2 else None
        stage_confirm(lamf)
        print("confirm done", lamf)
    elif mode == "verdict":
        stage_verdict()
    else:
        raise SystemExit("tryb: main [lam] | confirm | verdict")
