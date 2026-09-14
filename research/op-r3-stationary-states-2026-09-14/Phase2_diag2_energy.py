#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
op-r3-stationary-states -- DIAGNOSTYKA DESKRYPTYWNA 2 (energia P2c).
Cel: zlokalizowac skladnik C1 (~ -7e-6, saturacja t~300, pozornie
niezalezny od dt i h). NIE werdyktotworcza.

T1: skalowanie dt: dt in {0.01, 0.005, 0.00125} (h=0.05), pelne
    100 T0, srednie <E> po oknach 10 T0.
T2: odwracalnosc czasowa: forward 350, pi -> -pi, back 350;
    max|psi_back - psi_0| i E na koncu.
T3: dekompozycja przestrzenna E(t): r<=80 / (80,160] / (160,200]
    oraz <psi>-1 i max|psi-1|, co 10 T0 (dt=0.005).
"""
import sys
import time
import numpy as np

BASE = ("C:/Users/Mateusz/Documents/ObsydnianMain/TGP/TGP_v1/research/"
        "op-r3-stationary-states-2026-09-14/")
sys.path.insert(0, BASE)
import engine_core as ec

T0 = 2.0*np.pi
out = []
t0_wall = time.time()


def em(s=""):
    print(s, flush=True)
    out.append(s)


def stamp(msg):
    print("[t=%7.1fs] %s" % (time.time() - t0_wall, msg), flush=True)


def fresh(h=0.05):
    eng = ec.Engine(h, 200.0, sponge=False)
    g = 1.0 + 1e-3*np.exp(-eng.r**2/(2.0*3.0**2))
    pi = np.zeros(eng.N)
    return eng, g, pi


em("=" * 78)
em("DIAGNOSTYKA 2 -- energia P2c (deskryptywna)")
em("=" * 78)

# ---------------- T1: skalowanie dt --------------------------------
em("")
em("T1: plateau <E>-E0 wg dt (h=0.05, okna 10 T0)")
for dt in (0.01, 0.005, 0.00125):
    eng, g, pi = fresh()
    ts, Es = [0.0], [eng.energy(g, pi)]
    nso = max(1, int(round(0.1/dt)))
    for k in range(1, int(round(100.0*T0/dt)) + 1):
        g, pi = eng.step(g, pi, dt)
        if k % nso == 0:
            ts.append(k*dt)
            Es.append(eng.energy(g, pi))
        if k % 100000 == 0:
            stamp("  T1 dt=%g t=%.0f/628" % (dt, k*dt))
    ts = np.array(ts)
    Es = np.array(Es)
    E0 = Es[0]
    em("  dt=%.5f:" % dt)
    for w in (0, 1, 4, 9):
        m = (ts >= w*10.0*T0) & (ts < (w + 1)*10.0*T0)
        em("    okno [%3d,%3d] T0: (<E>-E0)/E0 = %+.4e"
           % (10*w, 10*(w + 1), (float(np.mean(Es[m])) - E0)/E0))

# ---------------- T2: odwracalnosc ---------------------------------
em("")
em("T2: odwracalnosc (dt=0.005, h=0.05): forward 350, pi->-pi,")
em("    back 350")
eng, g, pi = fresh()
g0 = g.copy()
E0 = eng.energy(g, pi)
for k in range(1, int(round(350.0/0.005)) + 1):
    g, pi = eng.step(g, pi, 0.005)
Emid = eng.energy(g, pi)
pi = -pi
for k in range(1, int(round(350.0/0.005)) + 1):
    g, pi = eng.step(g, pi, 0.005)
Eend = eng.energy(g, pi)
em("  E0=%.9e  E(350)=%.9e  (dE/E=%+.3e)" % (E0, Emid,
                                             (Emid - E0)/E0))
em("  po powrocie: max|psi-psi0| = %.3e ; max|pi| = %.3e ;")
em("  E_end=%.9e (dE/E=%+.3e)"
   % (Eend, (Eend - E0)/E0))
em("  (wiersz wyzej: max|psi-psi0|=%.3e, max|pi|=%.3e)"
   % (float(np.max(np.abs(g - g0))), float(np.max(np.abs(pi)))))

# ---------------- T3: dekompozycja przestrzenna --------------------
em("")
em("T3: dekompozycja E(t) (dt=0.005, h=0.05): core r<=80 /")
em("    mid (80,160] / wall (160,200]; <psi>-1; max|psi-1|")
eng, g, pi = fresh()
E0 = eng.energy(g, pi)


def report(t, g, pi):
    Ec = eng.energy(g, pi, 80.0)
    Em_ = eng.energy(g, pi, 160.0) - Ec
    Et = eng.energy(g, pi)
    Ew = Et - Ec - Em_
    em("  t=%6.1f: E_tot-E0=%+.4e  core=%.6e mid=%.6e wall=%.6e"
       "  <psi>-1=%+.2e max|psi-1|=%.2e"
       % (t, (Et - E0)/E0, Ec, Em_, Ew,
          float(np.mean(g)) - 1.0, float(np.max(np.abs(g - 1.0)))))


report(0.0, g, pi)
nrep = int(round(10.0*T0/0.005))
for k in range(1, int(round(100.0*T0/0.005)) + 1):
    g, pi = eng.step(g, pi, 0.005)
    if k % nrep == 0:
        report(k*0.005, g, pi)

with open(BASE + "Phase2_diag2_output.txt", "w",
          encoding="ascii") as fh:
    fh.write("\n".join(out) + "\n")
print("zapisano:", BASE + "Phase2_diag2_output.txt")
