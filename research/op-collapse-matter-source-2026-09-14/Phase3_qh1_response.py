#!/usr/bin/env python3
# -*- coding: ascii -*-
"""op-collapse-matter-source -- Phase 3 / Q-H1: odpowiedz prozni na
statyczne zrodlo (LOCK sec. 3 Phase 3 Q-H1; MD sec. 6).

Biegi: psi=1, pi0=0, zrodlo od t=0; lam in {0.01,0.05,0.2,0.5};
t_max=300; h=0.05 (lam=0.01 i 0.5 dodatkowo h=0.025); dt=0.005;
sponge ON. Klasyfikacja DEFORMATION / THRESHOLD-PULL / COLLAPSE /
INCONCLUSIVE-RUN wg MD sec. 6 (okna FROZEN). Gate liniowy przy
lam=0.01: profil osiadly vs -5 lam (G_Yuk*rho_hat), <=5% na r<=40.
Deskryptywnie: dpsi(0) zmierzone vs liniowe dla wszystkich lam.
"""
import json
import os
import sys
import time

import numpy as np
from scipy.integrate import quad

sys.path.insert(0, "TGP/TGP_v1/research/op-collapse-matter-source-2026-09-14")
import engine_core as ec

DIR = "TGP/TGP_v1/research/op-collapse-matter-source-2026-09-14"
RES = DIR + "/Phase3_results"
os.makedirs(RES, exist_ok=True)

OUT = []


def w(s=""):
    OUT.append(str(s))
    with open(DIR + "/Phase3_qh1_output.txt", "w") as fh:
        fh.write("\n".join(OUT) + "\n")


DT = 0.005
T_MAX = 300.0
DT_OUT = 0.1
DT_PROF = 5.0
DT_CHK = 100.0
LAMS = [0.01, 0.05, 0.2, 0.5]
GRIDS = {0.01: [0.05, 0.025], 0.05: [0.05], 0.2: [0.05],
         0.5: [0.05, 0.025]}
PSI_LO_THR = 5.0/6.0
PSI_HI_THR = 7.0/6.0
WIN = (250.0, 300.0)
RMAX_Q = 60.0


def rho_np(t):
    return np.exp(-t*t/18.0)


def dpsi_lin(rv, lamv):
    """Wzorzec P1-H2 (kwadratura scipy; MD sec. 4)."""
    if rv == 0.0:
        val = quad(lambda t: t*rho_np(t)*np.exp(-t), 0.0, RMAX_Q,
                   limit=200)[0]
        return -5.0*lamv*val

    def integ(t):
        return t*rho_np(t)*(np.exp(-abs(rv - t))
                            - np.exp(-(rv + t)))/2.0

    val = quad(integ, 0.0, RMAX_Q, points=[rv], limit=400)[0]
    return -5.0*lamv*val/rv


def run_one(lam, h):
    rid = "qh1_lam%g__h%s" % (lam, ("05" if h == 0.05 else "025"))
    jp = RES + "/" + rid + ".json"
    if os.path.exists(jp):
        with open(jp) as fh:
            rec = json.load(fh)
        dat = np.load(RES + "/" + rid + ".npz")
        pb = dat["psibar"]
        return (rec, (pb if pb.size else None),
                ec.Engine(h, 200.0, sponge=True, lam=lam))
    eng = ec.Engine(h, 200.0, sponge=True, lam=lam)
    g = ec.start_vacuum(eng)
    pi = np.zeros(eng.N)
    nst = int(round(T_MAX/DT))
    nout = int(round(DT_OUT/DT))
    nprof = int(round(DT_PROF/DT))
    nchk = int(round(DT_CHK/DT))
    m80 = eng.r <= 80.0
    ts, psi0, ecore, vmax = [0.0], [g[0]], [eng.energy(g, pi, 80.0)], \
        [0.0]
    pts, profs = [0.0], [g.copy()]
    status, t_end = "OK", None
    t1 = time.time()
    for n in range(nst):
        try:
            g, pi = eng.step(g, pi, DT)
        except ec.NonConvergence:
            status, t_end = "BREAKDOWN", (n + 1)*DT
            break
        st = eng.band_status(g)
        if st is not None:
            status, t_end = st, (n + 1)*DT
            break
        t = (n + 1)*DT
        if (n + 1) % nout == 0:
            ts.append(t)
            psi0.append(float(g[0]))
            ecore.append(eng.energy(g, pi, 80.0))
            vmax.append(float(np.max(np.abs((pi/ec.Mfun(g))[m80]))))
        if (n + 1) % nprof == 0:
            pts.append(t)
            profs.append(g.copy())
        if (n + 1) % nchk == 0:
            np.savez(RES + "/checkpoint_" + rid + ".npz", r=eng.r,
                     psi=g, pi=pi, t=t)
    wall = time.time() - t1
    ts = np.array(ts)
    psi0 = np.array(psi0)
    ecore = np.array(ecore)
    vmax = np.array(vmax)
    pts = np.array(pts)
    profs = np.array(profs)
    # klasyfikacja (MD sec. 6; priorytety FROZEN)
    det = {}
    if status != "OK":
        cls = "COLLAPSE"
        det["subtyp"] = status
        psibar = None
    else:
        mw = (pts >= WIN[0]) & (pts <= WIN[1])
        pw = profs[mw]
        psibar = np.mean(pw, axis=0)
        mv = (ts >= WIN[0]) & (ts <= WIN[1])
        V = float(np.max(vmax[mv]))
        D = float(np.max(np.abs(psibar[m80] - 1.0)))
        settled = V <= 0.01*max(D, 1e-12)
        pull_lo = bool(np.all(np.min(pw, axis=1) < PSI_LO_THR))
        pull_hi = bool(np.all(np.max(pw, axis=1) > PSI_HI_THR))
        within = bool(np.all((pw >= PSI_LO_THR) & (pw <= PSI_HI_THR)))
        det.update(V=V, D=D, settled=settled, pull_lo=pull_lo,
                   pull_hi=pull_hi, within=within,
                   psibar_min=float(np.min(psibar)),
                   psibar_ctr=float(psibar[0]))
        if settled and (pull_lo or pull_hi):
            cls = "THRESHOLD-PULL"
        elif settled and within:
            cls = "DEFORMATION"
        else:
            cls = "INCONCLUSIVE-RUN"
    np.savez(RES + "/" + rid + ".npz", r=eng.r, ts=ts, psi0=psi0,
             ecore=ecore, vmax=vmax, prof_ts=pts, profs=profs,
             psibar=(psibar if psibar is not None
                     else np.zeros(0)))
    rec = dict(rid=rid, lam=lam, h=h, status=status, cls=cls,
               t_end=t_end, wall=round(wall, 1), **det)
    with open(RES + "/" + rid + ".json", "w") as fh:
        json.dump(rec, fh, indent=1, default=str)
    return rec, psibar, eng


BAR = "=" * 78
w(BAR)
w("PHASE 3 / Q-H1 -- zrodlo na prozni (LOCK sec. 3; klasyfikatory MD"
  " sec. 6)")
w("REJESTR [INPUT]: lam {0.01,0.05,0.2,0.5}; t_max=300; h=0.05"
  " (0.01/0.5 + h=0.025);")
w("  dt=0.005; sponge ON; okno [250,300]; osiadlosc V<=0.01*max(D,"
  "1e-12) r<=80;")
w("  progi 5/6, 7/6 (trwale = wszystkie probki okna); pas 4/3-1e-6 /"
  " 1e-6;")
w("  gate liniowy lam=0.01: <=5% wzgledem max|dpsi_lin| na r<=40")
w(BAR)

I0 = quad(lambda t: t*rho_np(t)*np.exp(-t), 0.0, RMAX_Q, limit=200)[0]
recs = {}
bars = {}
engs = {}
for lam in LAMS:
    for h in GRIDS[lam]:
        rec, psibar, eng = run_one(lam, h)
        recs[(lam, h)] = rec
        bars[(lam, h)] = psibar
        engs[(lam, h)] = eng
        w("bieg lam=%-5g h=%.3f: status=%s cls=%s%s  [%.0f s]"
          % (lam, h, rec["status"], rec["cls"],
             ("" if rec["t_end"] is None else " t_end=%.3f"
              % rec["t_end"]), rec["wall"]))
        if psibar is not None:
            w("    V=%.3e D=%.3e settled=%s pull_lo=%s pull_hi=%s"
              " within=%s" % (rec["V"], rec["D"], rec["settled"],
                              rec["pull_lo"], rec["pull_hi"],
                              rec["within"]))
            w("    min psibar=%.6f  dpsi(0)=%+.6f"
              % (rec["psibar_min"], rec["psibar_ctr"] - 1.0))

w("-" * 78)
w("Gate liniowy (lam=0.01, h=0.05): profil osiadly vs"
  " -5 lam (G_Yuk*rho_hat)")
gate_lin = None
key = (0.01, 0.05)
if bars.get(key) is not None:
    eng = engs[key]
    m40 = eng.r <= 40.0
    rlin = eng.r[m40]
    lin = np.array([dpsi_lin(float(x), 0.01) for x in rlin])
    meas = bars[key][m40] - 1.0
    err = float(np.max(np.abs(meas - lin))/np.max(np.abs(lin)))
    gate_lin = err <= 0.05
    w("  err = max|dpsi_meas - dpsi_lin| / max|dpsi_lin| (r<=40)"
      " = %.4f  (prog 0.05)  %s" % (err, "PASS" if gate_lin
                                    else "FAIL"))
else:
    w("  NIEDOSTEPNY (bieg lam=0.01 h=0.05 bez profilu osiadlego)")

w("-" * 78)
w("Deskryptywnie (obowiazkowo): dpsi(0) zmierzone vs predykcja"
  " liniowa -5*lam*I0")
w("  I0 = %.9f (kwadratura P1-H2)" % I0)
w("  lam     h      klasa            dpsi(0)_meas   dpsi(0)_lin "
  "  meas/lin")
for lam in LAMS:
    for h in GRIDS[lam]:
        rec = recs[(lam, h)]
        linv = -5.0*lam*I0
        if rec["cls"] == "COLLAPSE":
            w("  %-6g  %-5.3f  %-15s  %-13s  %+.6f     -"
              % (lam, h, rec["cls"] + "(" + rec["subtyp"] + ")",
                 "n/a", linv))
        else:
            mv = rec["psibar_ctr"] - 1.0
            w("  %-6g  %-5.3f  %-15s  %+.6f      %+.6f     %.3f"
              % (lam, h, rec["cls"], mv, linv, mv/linv))

w("-" * 78)
w("Zbieznosc klas (MD sec. 6): lam 0.01/0.5 wymagaja h05=h025;"
  " lam 0.05/0.2 klasa z h05")
final = {}
for lam in LAMS:
    if len(GRIDS[lam]) == 2:
        c1 = recs[(lam, 0.05)]
        c2 = recs[(lam, 0.025)]
        agree = c1["cls"] == c2["cls"]
        if agree and c1["cls"] == "COLLAPSE":
            ta, tb = c1["t_end"], c2["t_end"]
            agree = abs(ta - tb)/max(ta, tb) <= 0.10
        final[lam] = c1["cls"] if agree else "INCONCLUSIVE-RUN"
        w("  lam=%-5g: h05=%s h025=%s => %s"
          % (lam, c1["cls"], c2["cls"],
             final[lam] + ("" if agree else " (rozjazd siatek)")))
    else:
        final[lam] = recs[(lam, 0.05)]["cls"]
        w("  lam=%-5g: h05=%s => %s" % (lam, final[lam], final[lam]))

all_def = all(v == "DEFORMATION" for v in final.values())
any_pull = any(v in ("THRESHOLD-PULL", "COLLAPSE")
               for v in final.values())
if all_def and gate_lin:
    verdict = "Q-H1-DEFORMATION"
elif any_pull:
    verdict = "Q-H1-PULL"
else:
    verdict = "Q-H1-INCONCLUSIVE"
w(BAR)
w("WERDYKT Q-H1 (litera LOCKa sec. 4): " + verdict)
w("  (DEFORMATION: wszystkie lam DEFORMATION zbieznie + gate liniowy"
  " <=5%;")
w("   PULL: >=1 lam THRESHOLD-PULL/COLLAPSE zbieznie; inaczej"
  " INCONCLUSIVE)")
w(BAR)
with open(RES + "/qh1_summary.json", "w") as fh:
    json.dump(dict(final={str(k): v for k, v in final.items()},
                   gate_lin=gate_lin, verdict=verdict,
                   runs={r["rid"]: r["cls"]
                         for r in recs.values()}), fh, indent=1)
print("done", verdict)
