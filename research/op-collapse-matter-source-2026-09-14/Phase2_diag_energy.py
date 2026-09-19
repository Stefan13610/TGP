#!/usr/bin/env python3
# -*- coding: ascii -*-
"""op-collapse-matter-source -- diagnostyka P2c (deskryptywna; ZERO
zmian progow/okien LOCKa). Pytanie: czy zmierzony 'dryf' okienny
2.28e-6 (>1e-6, FAIL) jest (a) zmiana SEKULARNA energii (=> gate
zasadnie FAIL => STOP), czy (b) aliasingiem OGRANICZONEJ oscylacji
energii wlasciwej metodom symetrycznym (bias estymatora okiennego --
blad pomiaru, nie trajektorii; precedens: correction note 2
poprzednika).

Testy (konfiguracja P2c: lam=0.05, gauss a=+0.05 s=3, sponge OFF,
h=0.05, t=630 >= 100 T0):
  T1: dt in {0.01, 0.005, 0.0025} -- skalowanie amplitudy oscylacji
      E(t) i 'dryfu' okiennego z dt (oczekiwane dt^2 przy (b)).
  T2: fit LSQ liniowy E(t) na [0,100T0] -> dryf sekularny na 100T0
      (estymator odporny na oscylacje ograniczona).
  T3: srednie E po kolejnych oknach 10T0 (monotonicznosc: trend
      sekularny vs wedrowanie oscylacyjne).
"""
import sys
import time

import numpy as np

sys.path.insert(0, "TGP/TGP_v1/research/op-collapse-matter-source-2026-09-14")
import engine_core as ec

DIR = "TGP/TGP_v1/research/op-collapse-matter-source-2026-09-14"
OUT = []


def w(s=""):
    OUT.append(str(s))
    with open(DIR + "/Phase2_diag_output.txt", "w") as fh:
        fh.write("\n".join(OUT) + "\n")


T0 = 2.0*np.pi
T_END = 630.0
w("=" * 78)
w("DIAGNOSTYKA P2c (deskryptywna; progi/okna LOCKa NIETKNIETE)")
w("konfiguracja: lam=0.05, gauss a=+0.05 s=3, sponge OFF, h=0.05,"
  " t=630")
w("=" * 78)

res = {}
for dt in (0.01, 0.005, 0.0025):
    eng = ec.Engine(0.05, 200.0, sponge=False, lam=0.05)
    g = ec.start_gauss(eng, 0.05, 3.0)
    pi = np.zeros(eng.N)
    nst = int(round(T_END/dt))
    nout = max(1, int(round(0.1/dt)))
    ts = [0.0]
    Es = [eng.energy(g, pi)]
    t1 = time.time()
    for n in range(nst):
        g, pi = eng.step(g, pi, dt)
        if (n + 1) % nout == 0:
            ts.append((n + 1)*dt)
            Es.append(eng.energy(g, pi))
    ts = np.array(ts)
    Es = np.array(Es)
    np.save(DIR + "/Phase3_results/diag_p2c_E_dt%g.npy" % dt,
            np.vstack([ts, Es]))
    m_e = (ts >= 0.0) & (ts <= 10*T0)
    m_l = (ts >= 90*T0) & (ts <= 100*T0)
    Ee = float(np.mean(Es[m_e]))
    El = float(np.mean(Es[m_l]))
    drift_win = abs(El - Ee)/abs(Ee)
    # T2: fit LSQ na [0,100T0]
    mf = ts <= 100*T0
    A = np.vstack([ts[mf], np.ones(int(np.sum(mf)))]).T
    sol, _, _, _ = np.linalg.lstsq(A, Es[mf], rcond=None)
    slope = float(sol[0])
    drift_fit = abs(slope)*100*T0/abs(Ee)
    # amplituda oscylacji wokol fitu
    resid = Es[mf] - (A @ sol)
    amp = float(np.max(np.abs(resid)))
    res[dt] = (drift_win, drift_fit, amp, slope)
    w("dt=%-7g: dryf_okienny=%.3e  dryf_fit(100T0)=%.3e"
      % (dt, drift_win, drift_fit))
    w("           amplituda oscylacji wokol fitu=%.3e  slope=%.3e"
      "  [%.0f s]" % (amp, slope, time.time() - t1))
    # T3: kolejne okna 10T0
    means = []
    for k in range(10):
        mk = (ts >= 10*k*T0) & (ts <= 10*(k + 1)*T0)
        means.append(float(np.mean(Es[mk])))
    dm = np.diff(means)
    w("  T3 srednie kolejnych okien 10T0 (odchylki od okna 1,"
      " x1e-6):")
    w("    " + " ".join("%+.2f" % ((m - means[0])/abs(means[0])*1e6)
                        for m in means))
    w("  T3 znaki roznic kolejnych okien: "
      + "".join("+" if d > 0 else "-" for d in dm)
      + "  (monotonia => trend sekularny; zmienne => oscylacja)")

w("-" * 78)
d1 = res[0.01]
d2 = res[0.005]
d3 = res[0.0025]
w("SKALOWANIE (T1): dt 0.01 -> 0.005 -> 0.0025")
w("  dryf_okienny : %.3e -> %.3e -> %.3e  (iloraz %.2f, %.2f;"
  " dt^2 => ~4)" % (d1[0], d2[0], d3[0],
                    d1[0]/max(d2[0], 1e-300),
                    d2[0]/max(d3[0], 1e-300)))
w("  amplituda osc: %.3e -> %.3e -> %.3e  (iloraz %.2f, %.2f)"
  % (d1[2], d2[2], d3[2], d1[2]/max(d2[2], 1e-300),
     d2[2]/max(d3[2], 1e-300)))
w("  dryf_fit     : %.3e -> %.3e -> %.3e" % (d1[1], d2[1], d3[1]))
w("=" * 78)
print("diag done")
