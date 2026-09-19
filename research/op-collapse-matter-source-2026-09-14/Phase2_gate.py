#!/usr/bin/env python3
# -*- coding: ascii -*-
"""op-collapse-matter-source -- Phase 2: bramka (LOCK sec. 3 Phase 2;
MD sec. 5). FAIL ktoregokolwiek gate'u => STOP cyklu.

P2a: proznia lam=0, sponge ON (konfiguracja produkcyjna), psi=1,
     100 T0, obie siatki; gate ||psi-1||_inf <= 1e-10 caly bieg.
P2b: regresja lam=0, start qR3 a=-0.20, h=0.05, sponge ON; gate:
     COLLAPSE (BREAKDOWN-BOUNDARY) z czasem zdarzenia 4.74 +/- 2%
     (baseline poprzednika Phase3_results/qR3_a-0.20__h05.json).
P2c: energia ZE zrodlem lam=0.05, start gauss a=+0.05 sigma=3,
     sponge OFF [INPUT-MD], t=700, obie siatki [INPUT-MD];
     dryf SEKULARNY (correction note 1 -- okno [0,10T0] zawiera
     dt^2-owy offset stanu poczatkowego, nie zmiane sekularna):
     dryf = (100/80)*|<E>_[90T0,100T0] - <E>_[10T0,20T0]|
            / |<E>_[0,10T0]|  <= 1e-6 /100T0;
     deskryptywnie: stary estymator (offset dC2) i fit LSQ
     na [10T0,100T0].
"""
import sys
import time

import numpy as np

sys.path.insert(0, "TGP/TGP_v1/research/op-collapse-matter-source-2026-09-14")
import engine_core as ec

OUT = []


def w(s=""):
    OUT.append(str(s))
    with open("TGP/TGP_v1/research/op-collapse-matter-source-2026-09-14/"
              "Phase2_output.txt", "w") as fh:
        fh.write("\n".join(OUT) + "\n")


T0 = 2.0*np.pi
DT = 0.005
BAR = "=" * 78
w(BAR)
w("PHASE 2 -- bramka (FROZEN; MD sec. 5). REJESTR [INPUT]: dt=0.005;"
  " R=200;")
w("  sponge gamma0=1.0 smootherstep [160,200] (ON prod., OFF P2c);"
  " dt_out=0.1;")
w("  okna dryfu [0,10T0]/[90T0,100T0]; progi: 1e-10 / 4.74+-2% / 1e-6")
w(BAR)

t_all0 = time.time()

# ------------------------------------------------------------- P2a
w("P2a: proznia lam=0, sponge ON, 100 T0 = %.6f" % (100*T0))
p2a_ok = True
for h in (0.05, 0.025):
    eng = ec.Engine(h, 200.0, sponge=True, lam=0.0)
    g = ec.start_vacuum(eng)
    pi = np.zeros(eng.N)
    nst = int(round(100*T0/DT))
    dev = 0.0
    t1 = time.time()
    for n in range(nst):
        g, pi = eng.step(g, pi, DT)
        d = float(np.max(np.abs(g - 1.0)))
        if d > dev:
            dev = d
    ok = dev <= 1e-10
    p2a_ok = p2a_ok and ok
    w("  h=%.3f: max||psi-1||_inf = %.3e  (prog 1e-10)  %s  [%.0f s]"
      % (h, dev, "PASS" if ok else "FAIL", time.time() - t1))
w("  P2a: " + ("PASS" if p2a_ok else "FAIL"))

# ------------------------------------------------------------- P2b
w("P2b: regresja lam=0, qR3 a=-0.20, h=0.05, sponge ON")
eng = ec.Engine(0.05, 200.0, sponge=True, lam=0.0)
g = ec.start_quasiR3(eng, -0.20)
pi = np.zeros(eng.N)
nst = int(round(30.0/DT))
t_end = None
status = None
t1 = time.time()
for n in range(nst):
    try:
        g, pi = eng.step(g, pi, DT)
    except ec.NonConvergence:
        status = "BREAKDOWN"
        t_end = (n + 1)*DT
        break
    st = eng.band_status(g)
    if st is not None:
        status = st
        t_end = (n + 1)*DT
        break
p2b_ok = (status == "BREAKDOWN-BOUNDARY" and t_end is not None
          and abs(t_end - 4.74) <= 0.02*4.74)
w("  status=%s t_end=%s (cel 4.74 +-2%%, tj. [4.6452, 4.8348])"
  % (status, "%.4f" % t_end if t_end else "n/a"))
w("  P2b: " + ("PASS" if p2b_ok else "FAIL") + "  [%.0f s]"
  % (time.time() - t1))

# ------------------------------------------------------------- P2c
w("P2c: energia ze zrodlem lam=0.05, gauss a=+0.05 s=3, sponge OFF,"
  " t=700")
p2c_ok = True
for h in (0.05, 0.025):
    eng = ec.Engine(h, 200.0, sponge=False, lam=0.05)
    g = ec.start_gauss(eng, 0.05, 3.0)
    pi = np.zeros(eng.N)
    nst = int(round(700.0/DT))
    nout = int(round(0.1/DT))
    ts = [0.0]
    Es = [eng.energy(g, pi)]
    t1 = time.time()
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
        w("  h=%.3f: PRZERWANE (%s)  FAIL" % (h, bad))
        p2c_ok = False
        continue
    ts = np.array(ts)
    Es = np.array(Es)
    m_0 = (ts >= 0.0) & (ts <= 10*T0)
    m_e = (ts >= 10*T0) & (ts <= 20*T0)
    m_l = (ts >= 90*T0) & (ts <= 100*T0)
    E0m = float(np.mean(Es[m_0]))
    Ee = float(np.mean(Es[m_e]))
    El = float(np.mean(Es[m_l]))
    # estymator SEKULARNY (correction note 1): okna plateau,
    # przeliczenie na 100T0 (odstep srodkow okien = 80T0)
    drift = (100.0/80.0)*abs(El - Ee)/abs(E0m)
    # deskryptywnie: stary estymator (offset dC2) + fit LSQ
    drift_old = abs(El - E0m)/abs(E0m)
    mf = (ts >= 10*T0) & (ts <= 100*T0)
    A = np.vstack([ts[mf], np.ones(int(np.sum(mf)))]).T
    sol = np.linalg.lstsq(A, Es[mf], rcond=None)[0]
    drift_fit = abs(float(sol[0]))*100*T0/abs(E0m)
    mx = float(np.max(np.abs(Es - Es[0])))
    ok = drift <= 1e-6
    p2c_ok = p2c_ok and ok
    w("  h=%.3f: E(0)=%.9e <E>_[10,20]T0=%.9e <E>_[90,100]T0=%.9e"
      % (h, Es[0], Ee, El))
    w("          dryf_sekularny=%.3e (prog 1e-6/100T0)  %s  [%.0f s]"
      % (drift, "PASS" if ok else "FAIL", time.time() - t1))
    w("          deskr.: offset_dC2(stary estym.)=%.3e"
      "  dryf_fitLSQ[10,100]T0=%.3e  max|E-E0|=%.3e"
      % (drift_old, drift_fit, mx))
w("  P2c: " + ("PASS" if p2c_ok else "FAIL"))

w(BAR)
allok = p2a_ok and p2b_ok and p2c_ok
w("PHASE 2 PODSUMOWANIE: P2a %s ; P2b %s ; P2c %s => %s  [%.0f s]"
  % ("PASS" if p2a_ok else "FAIL", "PASS" if p2b_ok else "FAIL",
     "PASS" if p2c_ok else "FAIL",
     "BRAMKA OTWARTA" if allok else "STOP (FAIL => STOP)",
     time.time() - t_all0))
w(BAR)
print("done")
