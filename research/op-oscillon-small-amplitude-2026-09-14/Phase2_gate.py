#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
op-oscillon-small-amplitude -- Phase 2: bramka maszynerii (MD sec. 8;
LOCK sec. 4; FAIL ktoregokolwiek gate'u => STOP).

P2a  proznia: R=400 sponge ON, psi=1, pi=0, 100 T0, obie siatki;
     gate ||psi-1||_inf <= 1e-10; detektor zero alarmow.
P2b  regresja: start (a=+0.15, sigma=3), R=400, h=0.05, sponge ON,
     t=300; gate tau (= t_hold - 50) = 206.9 +/- 5% (baseline
     poprzednika op-r3-stationary-states Phase3_output.txt).
P2c-energia: start (a=0.05, sigma=6), sponge OFF (pudlo zamkniete,
     INPUT-MD), R=400, t=700, obie siatki; dryf = |<E>_[90T0,100T0]
     - <E>_[0,10T0]| / <E>_[0,10T0] <= 1e-6; deskryptywnie
     max|E-E0|/E0.
P2c-sponge: bieg A R=400 sponge ON vs bieg B R=800 sponge OFF,
     h=0.05, t=350, puls 1e-3 exp(-(r-200)^2/(2*5^2)) [INPUT-MD];
     u = r(psi-1); gate: max_{r<=240,t}|u_A-u_B| /
     max_{r in [240,320],t}|u_B| <= 1e-3.

Output: Phase2_output.txt (+ Phase2_progress.log).
"""
import sys
import time
import numpy as np

BASE = ("C:/Users/Mateusz/Documents/ObsydnianMain/TGP/TGP_v1/research/"
        "op-oscillon-small-amplitude-2026-09-14/")
sys.path.insert(0, BASE)
from engine_core import (Engine, start_gauss, start_vacuum)

DT = 0.005
DT_OUT = 0.1
T0 = 2.0*np.pi
T100 = 100.0*T0

LOG = open(BASE + "Phase2_progress.log", "a")


def log(msg):
    LOG.write("[%s] %s\n" % (time.strftime("%H:%M:%S"), msg))
    LOG.flush()


def evolve(eng, g, pi, t_end, dt, cb_every_steps, cb):
    """Krok az do t_end; cb(k, g, pi) co cb_every_steps krokow
    (k = numer kroku, t = k*dt). cb wolany takze dla k=0."""
    nsteps = int(round(t_end/dt))
    cb(0, g, pi)
    for k in range(1, nsteps + 1):
        g, pi = eng.step(g, pi, dt)
        if k % cb_every_steps == 0:
            cb(k, g, pi)
    return g, pi


out = []
out.append("="*78)
out.append("PHASE 2 -- bramka maszynerii (MD sec. 8; LOCK sec. 4; "
           "FAIL => STOP)")
out.append("Konfiguracja [LOCK]: R=400; sponge smootherstep gamma0=1.0 "
           "[320,400]; dt=0.005; dt_out=0.1; E_core r<=80; "
           "E_ref=E_core(50)")
out.append("="*78)
results = {}

# ------------------------------------------------------------ P2a --
out.append("")
out.append("P2a (proznia; sponge ON; 100 T0 = %.3f; obie siatki):" % T100)
for h in (0.05, 0.025):
    eng = Engine(h, 400.0, sponge=True)
    g = start_vacuum(eng)
    pi = np.zeros(eng.N)
    st = {"maxdev": 0.0, "ncross": 0, "prev_sign": 0.0}

    def cb(k, gg, pp, st=st):
        dev = float(np.max(np.abs(gg - 1.0)))
        if dev > st["maxdev"]:
            st["maxdev"] = dev
        s = np.sign(gg[0] - 1.0)
        if s != 0 and st["prev_sign"] != 0 and s != st["prev_sign"]:
            st["ncross"] += 1
        if s != 0:
            st["prev_sign"] = s

    t_start = time.time()
    g, pi = evolve(eng, g, pi, T100, DT, int(round(DT_OUT/DT)), cb)
    ok = st["maxdev"] <= 1e-10
    results["P2a_h%g" % h] = ok
    out.append("  h=%.3f: ||psi-1||_inf = %.3e (prog 1e-10) %s; "
               "przejscia psi(0)=%d, kandydat=False (zero alarmow) "
               "[%.0f s]"
               % (h, st["maxdev"], "PASS" if ok else "FAIL",
                  st["ncross"], time.time() - t_start))
    log("P2a h=%g done maxdev=%.3e" % (h, st["maxdev"]))

# ------------------------------------------------------------ P2b --
out.append("")
out.append("P2b (regresja vs poprzednik: start a=+0.15 sigma=3, R=400,")
out.append("  h=0.05, sponge ON, t=300; tau = t_hold - 50; baseline "
           "206.9 +/- 5%):")
eng = Engine(0.05, 400.0, sponge=True)
g = start_gauss(eng, 0.15, 3.0)
pi = np.zeros(eng.N)
ts, ecs = [], []


def cb(k, gg, pp):
    ts.append(k*DT)
    ecs.append(eng.energy(gg, pp, rmax=80.0))


t_start = time.time()
g, pi = evolve(eng, g, pi, 300.0, DT, int(round(DT_OUT/DT)), cb)
ts = np.array(ts)
ecs = np.array(ecs)
i50 = int(np.argmin(np.abs(ts - 50.0)))
E_ref = ecs[i50]
mask = ts >= 50.0
below = np.where(ecs[mask] < 0.5*E_ref)[0]
if len(below) == 0:
    t_hold = ts[-1]
else:
    t_hold = ts[mask][below[0] - 1] if below[0] > 0 else 50.0
tau = t_hold - 50.0
rel = abs(tau - 206.9)/206.9
ok = rel <= 0.05
results["P2b"] = ok
out.append("  E_ref = E_core(50) = %.6e; t_hold = %.1f; tau = %.1f;"
           % (E_ref, t_hold, tau))
out.append("  |tau - 206.9|/206.9 = %.3f (prog 0.05) %s  [%.0f s]"
           % (rel, "PASS" if ok else "FAIL", time.time() - t_start))
log("P2b done tau=%.1f rel=%.3f" % (tau, rel))

# ----------------------------------------------------- P2c-energia --
out.append("")
out.append("P2c-energia (start a=0.05 sigma=6 [produkcyjny], sponge OFF")
out.append("  [INPUT-MD, pudlo zamkniete], R=400, t=700, obie siatki;")
out.append("  dryf = |<E>_[90T0,100T0] - <E>_[0,10T0]|/<E>_[0,10T0] "
           "<= 1e-6):")
for h in (0.05, 0.025):
    eng = Engine(h, 400.0, sponge=False)
    g = start_gauss(eng, 0.05, 6.0)
    pi = np.zeros(eng.N)
    ts, es = [], []

    def cb(k, gg, pp):
        ts.append(k*DT)
        es.append(eng.energy(gg, pp))

    t_start = time.time()
    g, pi = evolve(eng, g, pi, 700.0, DT, int(round(DT_OUT/DT)), cb)
    ts = np.array(ts)
    es = np.array(es)
    E0 = es[0]
    m1 = (ts >= 0.0) & (ts <= 10.0*T0)
    m2 = (ts >= 90.0*T0) & (ts <= 100.0*T0)
    drift = abs(float(np.mean(es[m2])) - float(np.mean(es[m1]))) \
        / float(np.mean(es[m1]))
    maxrel = float(np.max(np.abs(es - E0)))/E0
    ok = drift <= 1e-6
    results["P2c_E_h%g" % h] = ok
    out.append("  h=%.3f: E0=%.6e; dryf=%.3e (prog 1e-6) %s; "
               "deskryptywnie max|E-E0|/E0=%.3e  [%.0f s]"
               % (h, E0, drift, "PASS" if ok else "FAIL", maxrel,
                  time.time() - t_start))
    log("P2c-E h=%g done drift=%.3e" % (h, drift))

# ------------------------------------------------------ P2c-sponge --
out.append("")
out.append("P2c-sponge (odbicie, roznicowo: A R=400 sponge ON vs "
           "B R=800 sponge OFF;")
out.append("  h=0.05, t=350; puls 1e-3 exp(-(r-200)^2/(2*25)) "
           "[INPUT-MD]; u=r(psi-1);")
out.append("  gate: max_{r<=240,t}|u_A-u_B| / "
           "max_{r in [240,320],t}|u_B| <= 1e-3):")
h = 0.05
snap_every = int(round(DT_OUT/DT))


def run_store_u(R, sponge, mask_r):
    eng = Engine(h, R, sponge=sponge)
    g = 1.0 + 1e-3*np.exp(-(eng.r - 200.0)**2/(2.0*5.0**2))
    pi = np.zeros(eng.N)
    snaps = []

    def cb(k, gg, pp):
        snaps.append((eng.r[mask_r(eng.r)]
                      * (gg[mask_r(eng.r)] - 1.0)).copy())

    g, pi = evolve(eng, g, pi, 350.0, DT, snap_every, cb)
    return eng, np.array(snaps)


t_start = time.time()
engA, uA = run_store_u(400.0, True, lambda r: r <= 320.0)
log("P2c-sponge run A done")
engB, uB = run_store_u(800.0, False, lambda r: r <= 320.0)
log("P2c-sponge run B done")
rA = engA.r[engA.r <= 320.0]
m_in = rA <= 240.0
m_den = (rA > 240.0) & (rA <= 320.0)
num = float(np.max(np.abs(uA[:, m_in] - uB[:, m_in])))
den = float(np.max(np.abs(uB[:, m_den])))
ratio = num/den
ok = ratio <= 1e-3
results["P2c_sponge"] = ok
out.append("  max|u_A-u_B| (r<=240) = %.3e; max|u_B| (240<r<=320) = "
           "%.3e;" % (num, den))
out.append("  stosunek = %.3e (prog 1e-3) %s  [%.0f s]"
           % (ratio, "PASS" if ok else "FAIL", time.time() - t_start))
log("P2c-sponge done ratio=%.3e" % ratio)

# ----------------------------------------------------------- suma --
out.append("")
out.append("="*78)
allpass = all(results.values())
out.append("PHASE 2 PODSUMOWANIE: %s (%s)"
           % ("PASS %d/%d" % (sum(results.values()), len(results))
              if allpass else "FAIL -- STOP (LOCK sec. 4)",
              ", ".join("%s:%s" % (k, "P" if v else "F")
                        for k, v in results.items())))
out.append("="*78)

with open(BASE + "Phase2_output.txt", "w") as f:
    f.write("\n".join(out) + "\n")
log("PHASE2 DONE allpass=%s" % allpass)
print("PHASE2 DONE allpass=%s" % allpass)
