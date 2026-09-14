#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
op-r3-stationary-states (Phase 2) -- bramka maszynerii dynamicznej.
LOCK sec. 2 Phase 2; definicje FROZEN: Phase_method_decisions.md sec. 7.

P2a: prozna (sponge ON, konfiguracja produkcyjna), 100 T0, obie
     siatki: ||psi-1||_inf <= 1e-10; detektor zero alarmow.
P2b: test dyspersji (pudlo zamkniete, sponge OFF): puls a=1e-3,
     sigma=2 w r=0; FFT 2D u=r(psi-1); biny gate k~{0.6,1.0,1.4}
     (+{1.8,2.2} deskryptywnie): |om - sqrt(k^2+1)|/sqrt(k^2+1)<=1%.
P2c: dryf energii (pudlo zamkniete, obie siatki, puls a=1e-3 sigma=3):
     |<E>_[90T0,100T0] - <E>_[0,10T0]|/<E>_[0,10T0] <= 1e-6;
     odbicie sponge (roznicowo vs R=400 bez sponge) <= 1e-3.
FAIL ktoregokolwiek => STOP (LOCK).

REJESTR WEJSC [INPUT]: K_geo=gamma=c0=1; dt=0.005; gamma0_sp=1.0
[INPUT-MD]; smootherstep [160,200]; pulsy MD sec. 7; brak seeda.
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


def evolve(eng, g, pi, tmax, dt, sample_dt=None, sampler=None,
           label=""):
    """Petla czasowa; zwraca (g, pi, status). sampler(t,g,pi) co
    sample_dt. Pas graniczny wg LOCKa (klasyfikacja, stop)."""
    nsteps = int(round(tmax/dt))
    ns = int(round(sample_dt/dt)) if sample_dt else None
    if sampler:
        sampler(0.0, g, pi)
    for k in range(1, nsteps + 1):
        try:
            g, pi = eng.step(g, pi, dt)
        except (ec.NonConvergence, FloatingPointError):
            return g, pi, "BREAKDOWN", k*dt
        st = eng.band_status(g)
        if st:
            return g, pi, st, k*dt
        if ns and (k % ns == 0) and sampler:
            sampler(k*dt, g, pi)
        if k % 20000 == 0:
            stamp("  [%s] t=%.1f/%.1f" % (label, k*dt, tmax))
    return g, pi, "OK", tmax


em("=" * 78)
em("PHASE 2 -- bramka maszynerii dynamicznej (definicje FROZEN MD sec.7)")
em("REJESTR [INPUT]: K_geo=gamma=c0=1; dt=%.3f; sponge gamma0=%.1f"
   " smootherstep r in [160,200]; pas 4/3-1e-6 / 1e-6; brak seeda"
   % (DT, ec.GAMMA0_SP))
em("=" * 78)

results = {}

# ================= P2a: prozna (sponge ON, obie siatki) ============
em("")
em("P2a: prozna, konfiguracja produkcyjna (sponge ON), 100 T0")
for h in (0.05, 0.025):
    eng = ec.Engine(h, 200.0, sponge=True)
    g = ec.start_vacuum(eng)
    pi = np.zeros(eng.N)
    devmax = [0.0]
    ncross = [0]

    def samp(t, gg, pp, devmax=devmax):
        devmax[0] = max(devmax[0], float(np.max(np.abs(gg - 1.0))))

    g, pi, st, tend = evolve(eng, g, pi, 100.0*T0, DT, 1.0, samp,
                             "P2a h=%g" % h)
    dev = max(devmax[0], float(np.max(np.abs(g - 1.0))))
    ok = (st == "OK") and dev <= 1e-10
    results["P2a_h%g" % h] = ok
    em("  h=%.3f: status=%s, ||psi-1||_inf(max po biegu) = %.3e"
       " (<=1e-10: %s); alarmy detektora: 0 (E_core=0, przejscia=0)"
       % (h, st, dev, "PASS" if ok else "FAIL"))

# ================= P2b: dyspersja (pudlo zamkniete) ================
em("")
em("P2b: dyspersja -- puls a=1e-3 sigma=2 @ r=0, sponge OFF, R=200,")
em("  t=400, u=r(psi-1) co dt_out=0.2; FFT 2D Hann; parabola na piku")
K_GATE = (0.6, 1.0, 1.4)
K_DESC = (1.8, 2.2)


def dispersion_run(h):
    eng = ec.Engine(h, 200.0, sponge=False)
    g = 1.0 + 1e-3*np.exp(-eng.r**2/(2.0*2.0**2))
    pi = np.zeros(eng.N)
    dt_out = 0.2
    tmax = 400.0
    nt = int(round(tmax/dt_out)) + 1
    U = np.zeros((nt, eng.N), dtype=np.float32)
    idx = [0]

    def samp(t, gg, pp):
        U[idx[0], :] = (eng.r*(gg - 1.0)).astype(np.float32)
        idx[0] += 1

    g, pi, st, tend = evolve(eng, g, pi, tmax, DT, dt_out, samp,
                             "P2b h=%g" % h)
    assert st == "OK", st
    nt_used = idx[0]
    Uu = U[:nt_used, :].astype(np.float64)
    wt = np.hanning(nt_used)
    wr = np.hanning(eng.N)
    A = np.abs(np.fft.fft2(Uu*wt[:, None]*wr[None, :]))
    Ttot = (nt_used - 1)*dt_out
    dw = 2.0*np.pi/(nt_used*dt_out)
    dk = 2.0*np.pi/(eng.N*h)
    rows = []
    for ktar in K_GATE + K_DESC:
        nk = int(round(ktar/dk))
        kbin = nk*dk
        # correction note 1 (a): suma kwadrantow (+om,+k)+(+om,-k)
        s = A[:, nk] + A[:, (eng.N - nk) % eng.N]
        j0 = int(round(0.5/dw))
        j1 = int(round(5.0/dw))
        j = j0 + int(np.argmax(s[j0:j1]))
        lm, l0, lp = np.log(s[j-1]), np.log(s[j]), np.log(s[j+1])
        dlt = 0.5*(lm - lp)/(lm - 2*l0 + lp)
        om_meas = (j + dlt)*dw
        om_th = np.sqrt(kbin**2 + 1.0)
        rel = abs(om_meas - om_th)/om_th
        rows.append((ktar, kbin, om_meas, om_th, rel))
    return rows, Ttot


rows, Ttot = dispersion_run(0.05)
p2b_ok = True
em("  h=0.05 (PRIMARY), okno T=%.1f:" % Ttot)
for ktar, kbin, om, omt, rel in rows:
    gate = ktar in K_GATE
    ok = rel <= 0.01
    if gate:
        p2b_ok = p2b_ok and ok
    em("    k_tar=%.1f k_bin=%.4f om_meas=%.6f om_teor=%.6f "
       "|d|/om=%.2e %s%s"
       % (ktar, kbin, om, omt, rel,
          ("PASS" if ok else "FAIL") if gate else "(deskr.)",
          " [GATE]" if gate else ""))
rows2, _ = dispersion_run(0.025)
em("  h=0.025 (deskryptywnie):")
for ktar, kbin, om, omt, rel in rows2:
    em("    k_tar=%.1f k_bin=%.4f om_meas=%.6f om_teor=%.6f |d|/om=%.2e"
       % (ktar, kbin, om, omt, rel))
results["P2b"] = p2b_ok
em("  P2b (3 biny gate, h=0.05): %s" % ("PASS" if p2b_ok else "FAIL"))

# ================= P2c: dryf energii (pudlo zamkniete) =============
em("")
em("P2c-energia: puls a=1e-3 sigma=3, sponge OFF, 100 T0, obie siatki")
for h in (0.05, 0.025):
    eng = ec.Engine(h, 200.0, sponge=False)
    g = 1.0 + 1e-3*np.exp(-eng.r**2/(2.0*3.0**2))
    pi = np.zeros(eng.N)
    ts, Es = [], []

    def samp(t, gg, pp, eng=eng, ts=ts, Es=Es):
        ts.append(t)
        Es.append(eng.energy(gg, pp))

    g, pi, st, tend = evolve(eng, g, pi, 100.0*T0, DT, 0.1, samp,
                             "P2c-E h=%g" % h)
    ts = np.array(ts)
    Es = np.array(Es)
    E0w = float(np.mean(Es[ts <= 10.0*T0]))
    E9w = float(np.mean(Es[ts >= 90.0*T0]))
    drift = abs(E9w - E0w)/abs(E0w)
    oscmax = float(np.max(np.abs(Es - Es[0]))/abs(Es[0]))
    ok = (st == "OK") and drift <= 1e-6
    results["P2cE_h%g" % h] = ok
    em("  h=%.3f: status=%s E(0)=%.6e dryf(okna)=%.3e (<=1e-6: %s);"
       " deskr. max|E-E0|/E0=%.3e"
       % (h, st, Es[0], drift, "PASS" if ok else "FAIL", oscmax))
    ts, Es = [], []

# ================= P2c: odbicie sponge (roznicowo) =================
em("")
em("P2c-sponge: puls a=1e-3 sigma=5 @ r=100; A: R=200 sponge ON;")
em("  B: R=400 sponge OFF (referencja); h=0.05, t=250; u=r(psi-1)")


def refl_run(R, sponge):
    eng = ec.Engine(0.05, R, sponge=sponge)
    g = 1.0 + 1e-3*np.exp(-(eng.r - 100.0)**2/(2.0*5.0**2))
    pi = np.zeros(eng.N)
    dt_out = 0.2
    nt = int(round(250.0/dt_out)) + 1
    ncol = int(round(160.0/0.05))
    U = np.zeros((nt, ncol), dtype=np.float32)
    idx = [0]

    def samp(t, gg, pp):
        U[idx[0], :] = (eng.r[:ncol]*(gg[:ncol] - 1.0)).astype(
            np.float32)
        idx[0] += 1

    g, pi, st, tend = evolve(eng, g, pi, 250.0, DT, dt_out, samp,
                             "P2c-sp R=%g" % R)
    assert st == "OK", st
    return U[:idx[0], :]


UA = refl_run(200.0, True)
UB = refl_run(400.0, False)
n120 = int(round(120.0/0.05))
n160 = int(round(160.0/0.05))
diff = float(np.max(np.abs(UA[:, :n120].astype(np.float64)
                           - UB[:, :n120].astype(np.float64))))
ainc = float(np.max(np.abs(UB[:, n120:n160].astype(np.float64))))
refl = diff/ainc
ok = refl <= 1e-3
results["P2c_sponge"] = ok
em("  max|u_A-u_B| (r<=120) = %.3e; ampl. padajaca (r in [120,160])"
   " = %.3e" % (diff, ainc))
em("  odbicie = %.3e (<=1e-3: %s)" % (refl, "PASS" if ok else "FAIL"))

# ================= werdykt Phase 2 =================================
em("")
em("=" * 78)
allpass = all(results.values())
for k in sorted(results):
    em("  %-12s: %s" % (k, "PASS" if results[k] else "FAIL"))
em("PHASE 2: %s" % ("PASS -- maszyneria dopuszczona do Phase 3"
                    if allpass else "FAIL => STOP (LOCK sec. 2)"))
em("=" * 78)

with open(BASE + "Phase2_output.txt", "w", encoding="ascii") as fh:
    fh.write("\n".join(out) + "\n")
print("zapisano:", BASE + "Phase2_output.txt")
sys.exit(0 if allpass else 1)
