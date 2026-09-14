#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
op-r3-stationary-states -- DIAGNOSTYKA DESKRYPTYWNA dwoch FAIL-i
Phase 2 (pre-correction). NIE werdyktotworcza; sluzy udokumentowaniu
correction note PRZED poprawka harnessu. Zero zmian progow/definicji
LOCKa; ewolucje deterministyczne identyczne z Phase 2.

D1 (P2b): struktura kwadrantow FFT 2D u(r,t) -- hipoteza: odczyt
    A[+omega, +k] niesie skladowa PRZYCHODZACA; wychodzaca zyje
    w A[+omega, N-k]. Raport pikow w obu kwadrantach i w sumie.
D2 (P2c-E): ksztalt E(t) -- hipoteza: jednorazowy offset O(dt^2)
    (hamiltonian-cien), nie dryf sekularny. Raport srednich E po
    oknach 10 T0 oraz skalowanie offsetu przy dt/2.
"""
import sys
import time
import numpy as np

BASE = ("C:/Users/Mateusz/Documents/ObsydnianMain/TGP/TGP_v1/research/"
        "op-r3-stationary-states-2026-09-14/")
sys.path.insert(0, BASE)
import engine_core as ec

DT = 0.005
T0 = 2.0*np.pi
out = []
t0_wall = time.time()


def em(s=""):
    print(s, flush=True)
    out.append(s)


def stamp(msg):
    print("[t=%7.1fs] %s" % (time.time() - t0_wall, msg), flush=True)


em("=" * 78)
em("DIAGNOSTYKA DESKRYPTYWNA Phase 2 (pre-correction; nie werdykt)")
em("=" * 78)

# ---------------- D1: kwadranty FFT (dyspersja h=0.05) -------------
em("")
em("D1: dyspersja h=0.05 -- piki w kwadrantach (+om,+k), (+om,-k),")
em("    suma; ewolucja identyczna z P2b (deterministyczna)")
eng = ec.Engine(0.05, 200.0, sponge=False)
g = 1.0 + 1e-3*np.exp(-eng.r**2/(2.0*2.0**2))
pi = np.zeros(eng.N)
dt_out = 0.2
nt = int(round(400.0/dt_out)) + 1
U = np.zeros((nt, eng.N), dtype=np.float32)
U[0, :] = (eng.r*(g - 1.0)).astype(np.float32)
ns = int(round(dt_out/DT))
for k in range(1, int(round(400.0/DT)) + 1):
    g, pi = eng.step(g, pi, DT)
    if k % ns == 0:
        U[k//ns, :] = (eng.r*(g - 1.0)).astype(np.float32)
    if k % 20000 == 0:
        stamp("  D1 t=%.0f/400" % (k*DT))
Uu = U.astype(np.float64)
wt = np.hanning(nt)
wr = np.hanning(eng.N)
A = np.abs(np.fft.fft2(Uu*wt[:, None]*wr[None, :]))
dw = 2.0*np.pi/(nt*dt_out)
dk = 2.0*np.pi/(eng.N*0.05)
j0, j1 = int(round(0.5/dw)), int(round(5.0/dw))


def peak(s):
    j = j0 + int(np.argmax(s[j0:j1]))
    lm, l0, lp = np.log(s[j-1]), np.log(s[j]), np.log(s[j+1])
    d = 0.5*(lm - lp)/(lm - 2*l0 + lp)
    return (j + d)*dw, s[j]


for ktar in (0.6, 1.0, 1.4):
    nk = int(round(ktar/dk))
    kbin = nk*dk
    om_th = np.sqrt(kbin**2 + 1.0)
    sP = A[:, nk]                       # kwadrant (+om, +k)
    sM = A[:, (eng.N - nk) % eng.N]     # kwadrant (+om, -k)
    sS = sP + sM
    omP, aP = peak(sP)
    omM, aM = peak(sM)
    omS, aS = peak(sS)
    em("  k_bin=%.4f om_teor=%.6f:" % (kbin, om_th))
    em("    (+om,+k):  om_pik=%.6f  amp=%.3e  |d|/om=%.2e"
       % (omP, aP, abs(omP - om_th)/om_th))
    em("    (+om,-k):  om_pik=%.6f  amp=%.3e  |d|/om=%.2e"
       % (omM, aM, abs(omM - om_th)/om_th))
    em("    suma:      om_pik=%.6f  amp=%.3e  |d|/om=%.2e"
       % (omS, aS, abs(omS - om_th)/om_th))
    em("    stosunek amp (+k)/(-k) w piku sumy: %.3e"
       % (sP[int(round(omS/dw))]/max(sM[int(round(omS/dw))], 1e-300)))

# ---------------- D2: ksztalt E(t) i skalowanie dt -----------------
em("")
em("D2: energia zamknietego pudla h=0.05, puls a=1e-3 sigma=3;")
em("    srednie <E> po oknach 10 T0; dt=0.005 vs dt/2=0.0025")
for dt in (0.005, 0.0025):
    eng = ec.Engine(0.05, 200.0, sponge=False)
    g = 1.0 + 1e-3*np.exp(-eng.r**2/(2.0*3.0**2))
    pi = np.zeros(eng.N)
    ts, Es = [0.0], [eng.energy(g, pi)]
    nso = int(round(0.1/dt))
    for k in range(1, int(round(100.0*T0/dt)) + 1):
        g, pi = eng.step(g, pi, dt)
        if k % nso == 0:
            ts.append(k*dt)
            Es.append(eng.energy(g, pi))
        if k % 40000 == 0:
            stamp("  D2 dt=%g t=%.0f/628" % (dt, k*dt))
    ts = np.array(ts)
    Es = np.array(Es)
    E0 = Es[0]
    em("  dt=%.4f: E(0)=%.9e" % (dt, E0))
    means = []
    for w in range(10):
        m = (ts >= w*10.0*T0) & (ts < (w + 1)*10.0*T0)
        means.append(float(np.mean(Es[m])))
        em("    okno [%3d,%3d] T0: <E>=%.9e  (<E>-E0)/E0=%+.3e"
           % (10*w, 10*(w + 1), means[-1], (means[-1] - E0)/E0))
    d19 = abs(means[9] - means[1])/means[1]
    d09 = abs(means[9] - means[0])/means[0]
    em("    |<E>_[90,100]-<E>_[0,10]|/<E>_[0,10]   = %.3e"
       "  (harness pre-correction)" % d09)
    em("    |<E>_[90,100]-<E>_[10,20]|/<E>_[10,20] = %.3e"
       "  (dryf sekularny post-wirializacja)" % d19)

with open(BASE + "Phase2_diag_output.txt", "w", encoding="ascii") as fh:
    fh.write("\n".join(out) + "\n")
print("zapisano:", BASE + "Phase2_diag_output.txt")
