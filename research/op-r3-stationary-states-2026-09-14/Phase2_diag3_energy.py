#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
op-r3-stationary-states -- DIAGNOSTYKA DESKRYPTYWNA 3 (energia P2c):
lokalizacja skladnika C1 (dt-, h-niezaleznego, odwracalnego).
Hipoteza robocza: mierzone E nie jest scisle hamiltonianem
generujacym zaimplementowany przeplyw (niespojnosc ktoregos czlonu).

T4a: ewolucja bazowa do t=150 (dt=0.005), potem 20 j.cz. z dt in
     {0.02, 0.01, 0.005, 0.00125}: dE na odcinku vs dt.
T4b: bisekcja czlonow z tego samego stanu: warianty M=1 (bez c-term),
     K=1, oba -- ktory wylacza skladnik dt-niezalezny.
T4c: spojnosc gradientu: centralne roznice dE/dpsi_i vs -w_i F_i
     oraz dE/dpi_i vs w_i pi_i/M_i, na stanie t=150.
T4d: rzad amplitudy: pelny bieg a=1e-4 (dt=0.005): C1/E0 vs a.
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


class EngineVar(ec.Engine):
    """Wariant diagnostyczny: opcjonalnie M=1 i/lub K=1 (TYLKO
    diagnostyka; maszyneria produkcyjna nietknieta)."""

    def __init__(self, h, R, m1=False, k1=False):
        super().__init__(h, R, sponge=False)
        self.m1 = m1
        self.k1 = k1

    def Mf(self, p):
        return np.ones_like(p) if self.m1 else ec.Mfun(p)

    def Mpf(self, p):
        return np.zeros_like(p) if self.m1 else ec.Mpfun(p)

    def Kf(self, p):
        return np.ones_like(p) if self.k1 else ec.Kfun(p)

    def Kpf(self, p):
        return np.zeros_like(p) if self.k1 else ec.Kpfun(p)

    def F(self, g):
        h = self.h
        gm = 0.5*(g[:-1] + g[1:])
        dg = np.diff(g)/h
        dH = h*self.r2*ec.Upfun(g)
        t_flux = self.rm2*self.Kf(gm)*dg
        t_quad = 0.25*h*self.rm2*self.Kpf(gm)*dg**2
        dH[:-1] += -t_flux + t_quad
        dH[1:] += t_flux + t_quad
        return -dH/(h*self.r2)

    def energy(self, g, pi, rmax=None):
        gm = 0.5*(g[:-1] + g[1:])
        dg = np.diff(g)/self.h
        dens = self.r2*(pi*pi/(2.0*self.Mf(g))
                        + (ec.Ufun(g) - ec.U_VAC))
        grad = 0.5*self.rm2*self.Kf(gm)*dg**2
        if rmax is None:
            return 4.0*np.pi*self.h*(float(np.sum(dens))
                                     + float(np.sum(grad)))
        mp = self.r <= rmax
        mg = self.rm <= rmax
        return 4.0*np.pi*self.h*(float(np.sum(dens[mp]))
                                 + float(np.sum(grad[mg])))

    def step(self, g, pi, dt):
        F1 = self.F(g)
        M1 = self.Mf(g)
        c1 = self.Mpf(g)/(2.0*M1*M1)
        hdt = 0.5*dt
        ph = pi
        d_prev = np.inf
        for _ in range(ec.FP_MAXIT):
            new = pi + hdt*(F1 + ph*ph*c1)
            d = float(np.max(np.abs(new - ph)))
            ph = new
            if d == 0.0 or d >= d_prev:
                break
            d_prev = d
        invM1 = 1.0/M1
        x = g
        d_prev = np.inf
        for _ in range(ec.FP_MAXIT):
            new = g + hdt*ph*(invM1 + 1.0/self.Mf(x))
            d = float(np.max(np.abs(new - x)))
            x = new
            if d == 0.0 or d >= d_prev:
                break
            d_prev = d
        g2 = x
        M2 = self.Mf(g2)
        c2 = self.Mpf(g2)/(2.0*M2*M2)
        pi2 = ph + hdt*(self.F(g2) + ph*ph*c2)
        return g2, pi2


em("=" * 78)
em("DIAGNOSTYKA 3 -- lokalizacja C1 (deskryptywna)")
em("=" * 78)

# stan bazowy t=150
eng = ec.Engine(0.05, 200.0, sponge=False)
g = 1.0 + 1e-3*np.exp(-eng.r**2/(2.0*3.0**2))
pi = np.zeros(eng.N)
for k in range(1, int(round(150.0/0.005)) + 1):
    g, pi = eng.step(g, pi, 0.005)
gs, ps = g.copy(), pi.copy()
Es = eng.energy(gs, ps)
stamp("stan bazowy t=150 gotowy; E=%.9e" % Es)

em("")
em("T4a: dE na odcinku [150,170] vs dt (pelny model):")
for dt in (0.02, 0.01, 0.005, 0.00125):
    g, pi = gs.copy(), ps.copy()
    for k in range(int(round(20.0/dt))):
        g, pi = eng.step(g, pi, dt)
    dE = eng.energy(g, pi) - Es
    em("  dt=%.5f: (E(170)-E(150))/E = %+.4e" % (dt, dE/Es))

em("")
em("T4b: bisekcja czlonow -- ewolucja [0,150]+[150,170] w wariancie,")
em("  dE na [150,170] vs dt (uwaga: inny model => inna trajektoria;")
em("  diagnozujemy TYLKO obecnosc skladnika dt-niezaleznego):")
for m1, k1, lab in ((True, False, "M=1"), (False, True, "K=1"),
                    (True, True, "M=1,K=1")):
    ev = EngineVar(0.05, 200.0, m1=m1, k1=k1)
    g = 1.0 + 1e-3*np.exp(-ev.r**2/(2.0*3.0**2))
    pi = np.zeros(ev.N)
    for k in range(1, int(round(150.0/0.005)) + 1):
        g, pi = ev.step(g, pi, 0.005)
    gv, pv = g.copy(), pi.copy()
    Ev = ev.energy(gv, pv)
    row = []
    for dt in (0.02, 0.005, 0.00125):
        g, pi = gv.copy(), pv.copy()
        for k in range(int(round(20.0/dt))):
            g, pi = ev.step(g, pi, dt)
        row.append((ev.energy(g, pi) - Ev)/Ev)
    em("  %-8s: dt=0.02: %+.4e  dt=0.005: %+.4e  dt=0.00125: %+.4e"
       % (lab, row[0], row[1], row[2]))
    stamp("T4b %s done" % lab)

em("")
em("T4c: spojnosc gradientu na stanie t=150 (centralne roznice):")
w = 4.0*np.pi*eng.h*eng.r2
F = eng.F(gs)
idx = [0, 1, 2, 100, 500, 1000, 1500, 2000, 2500, 3000, 3998, 3999]
worst = 0.0
for i in idx:
    ep = 1e-7
    gp = gs.copy()
    gp[i] += ep
    gm_ = gs.copy()
    gm_[i] -= ep
    dnum = (eng.energy(gp, ps) - eng.energy(gm_, ps))/(2*ep)
    dana = -w[i]*F[i]
    rel = abs(dnum - dana)/max(abs(dana), 1e-30)
    worst = max(worst, rel)
    em("  i=%4d r=%7.2f  dE/dpsi_num=%+.6e  -w F=%+.6e  rel=%.2e"
       % (i, eng.r[i], dnum, dana, rel))
em("  max rel (psi) = %.2e" % worst)
worst = 0.0
for i in idx:
    ep = 1e-7
    pp = ps.copy()
    pp[i] += ep
    pm = ps.copy()
    pm[i] -= ep
    dnum = (eng.energy(gs, pp) - eng.energy(gs, pm))/(2*ep)
    dana = w[i]*ps[i]/ec.Mfun(gs[i])
    rel = abs(dnum - dana)/max(abs(dana), 1e-30)
    worst = max(worst, rel)
em("  max rel (pi)  = %.2e" % worst)

em("")
em("T4d: rzad amplitudy -- pelny bieg a=1e-4, dt=0.005, h=0.05:")
eng = ec.Engine(0.05, 200.0, sponge=False)
g = 1.0 + 1e-4*np.exp(-eng.r**2/(2.0*3.0**2))
pi = np.zeros(eng.N)
ts, Es_ = [0.0], [eng.energy(g, pi)]
for k in range(1, int(round(100.0*T0/0.005)) + 1):
    g, pi = eng.step(g, pi, 0.005)
    if k % 20 == 0:
        ts.append(k*0.005)
        Es_.append(eng.energy(g, pi))
    if k % 60000 == 0:
        stamp("T4d t=%.0f/628" % (k*0.005))
ts = np.array(ts)
Es_ = np.array(Es_)
E0 = Es_[0]
for wnd in (0, 1, 4, 9):
    m = (ts >= wnd*10.0*T0) & (ts < (wnd + 1)*10.0*T0)
    em("  okno [%3d,%3d] T0: (<E>-E0)/E0 = %+.4e"
       % (10*wnd, 10*(wnd + 1), (float(np.mean(Es_[m])) - E0)/E0))
em("  (a=1e-3 dawalo plateau-C2 ~ -7.0e-6; jesli C1/E0 ~ a =>"
   " teraz ~ -7e-7; jesli ~ a^2 => ~ -7e-8; jesli staly => -7e-6)")

with open(BASE + "Phase2_diag3_output.txt", "w",
          encoding="ascii") as fh:
    fh.write("\n".join(out) + "\n")
print("zapisano:", BASE + "Phase2_diag3_output.txt")
