#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
op-dynamics-class-M911 (Phase 1) -- bramki maszynerii (LOCK sec. 2 Phase 1).
Testowane SA rzeczywiste steppery silnikow (import z Phase2_conserved /
Phase3_inertial; wspolny gradient dyskretny M911_common).
G1 (A): proznia zostaje; masa zachowana (dip); E niemalejaco nie rosnie.
G2 (B): proznia zostaje; fala liniowa omega vs dyspersja dyskretna
  (formula LOCKa omega=sqrt(K(1)k^2+U''(1))/sqrt(B(1)) z k^2 =
  wartoscia wlasna zaimplementowanego -nabla^2 dla zasianego modu;
  odchylka od kontinuum sqrt(2) raportowana deskryptywnie);
  dryf H <= 1e-3 wzgledem energii fali; B' vs sympy 1e-12.
G3: detektory (dip/bump/proznia) 1+-0.
Dowolny FAIL ==> STOP.
"""
import sys
import numpy as np
import sympy as sp

BASE = ("C:/Users/Mateusz/Documents/ObsydnianMain/TGP/TGP_v1/research/"
        "op-dynamics-class-M911-2026-09-02/")
sys.path.insert(0, BASE)
from M911_common import (Grid3D, Bprime, detect, banner, start_dip,  # noqa: E402
                         L_GEN, L_LAT)
from Phase2_conserved import CHFlow  # noqa: E402
from Phase3_inertial import InertialFlow, DT_MAIN as DT_B  # noqa: E402

lines = []


def em(s=""):
    print(s, flush=True)
    lines.append(s)


em("=" * 78)
em("Phase 1 -- bramki maszynerii (silniki A i B)")
banner("; Phase1 gate")
em("=" * 78)
allpass = True

# ------------------------------------------------------------------ G1
em("")
em("G1 (wariant A, Cahn-Hilliard):")
for L, name in ((L_LAT, "L=2pi"), (L_GEN, "L=4pi")):
    fl = CHFlow(32, L)
    g = np.ones((32, 32, 32))
    drift = 0.0
    for k in range(1000):
        g = fl.step(g, 0.01)
        drift = max(drift, float(np.max(np.abs(g - 1.0))))
    ok = drift <= 1e-10
    allpass &= ok
    em("  proznia N=32 %s t=10: dryf=%.3e -> %s"
       % (name, drift, "OK" if ok else "FAIL"))
fl = CHFlow(32, L_GEN)
g = start_dip(32)
m0 = float(g.mean())
E_prev = fl.energy(g)
mono_ok = True
worst_rise = 0.0
for k in range(1, 1001):
    g = fl.step(g, 0.01)
    if k % 100 == 0:
        E = fl.energy(g)
        rise = E - E_prev
        worst_rise = max(worst_rise, rise)
        if rise > 1e-12:
            mono_ok = False
        E_prev = E
dm = abs(float(g.mean()) - m0) / abs(m0)
ok_m = dm <= 1e-13
allpass &= ok_m and mono_ok
em("  masa (dip N=32, t=10): |d mean|/mean = %.3e (<=1e-13) -> %s"
   % (dm, "OK" if ok_m else "FAIL"))
em("  E niemalejaco (probki co 1): max wzrost = %.3e (<=1e-12) -> %s"
   % (worst_rise, "OK" if mono_ok else "FAIL"))
em("  G1: %s" % ("PASS" if (ok_m and mono_ok) else "FAIL"))

# ------------------------------------------------------------------ G2
em("")
em("G2 (wariant B, inercja+lapse):")
p = sp.symbols('p', positive=True)
Bs = p ** 6 / (4 - 3 * p) ** 2
Bps = sp.diff(Bs, p)
worst = 0.0
for x in (sp.Rational(1, 2), sp.Integer(1), sp.Rational(7, 6),
          sp.Rational(13, 10)):
    ex = float(sp.N(Bps.subs(p, x), 30))
    im = Bprime(float(x))
    worst = max(worst, abs(ex - im) / max(abs(ex), 1.0))
ok_bp = worst <= 1e-12
allpass &= ok_bp
em("  B'(psi)=12psi^5(2-psi)/(4-3psi)^3 vs sympy: max|d|_rel=%.3e -> %s"
   % (worst, "OK" if ok_bp else "FAIL"))

fl = InertialFlow(32, L_GEN)
u = np.ones((32, 32, 32))
v = np.zeros_like(u)
drift = 0.0
for k in range(int(round(10.0 / DT_B))):
    u, v = fl.step(u, v, DT_B)
    drift = max(drift, float(np.max(np.abs(u - 1.0))))
ok_v = drift <= 1e-10
allpass &= ok_v
em("  proznia N=32 L=4pi t=10 (v0=0): dryf=%.3e -> %s"
   % (drift, "OK" if ok_v else "FAIL"))

N = 32
fl = InertialFlow(N, L_GEN)
h = L_GEN / N
kx = 2 * np.pi * 2 / L_GEN                       # mod m=2 -> k=1
x = (np.arange(N) + 0.5) * h
cosx = np.cos(kx * x)[:, None, None] * np.ones((1, N, N))
eps = 1e-6                                       # LOCK
u = 1.0 + eps * cosx
v = np.zeros_like(u)
H0 = fl.H(u, v)
Hrel0 = H0 - fl.Hvac(u)
keff2 = (2 - 2 * np.cos(kx * h)) / h ** 2
om_disc = float(np.sqrt(keff2 * 1.0 + 1.0))      # K(1)=B(1)=U''(1)=1
amps, times, Hmax = [], [], 0.0
nst = int(round(10.0 / DT_B))
for k in range(1, nst + 1):
    u, v = fl.step(u, v, DT_B)
    amps.append(float(np.mean(2 * (u - 1.0) * cosx)))
    times.append(k * DT_B)
    if k % 400 == 0:
        Hmax = max(Hmax, abs(fl.H(u, v) - H0))
amps = np.array(amps)
times = np.array(times)
sgn = np.sign(amps)
idx = np.where(sgn[1:] * sgn[:-1] < 0)[0]
tc = times[idx] - amps[idx] * (times[idx + 1] - times[idx]) / \
    (amps[idx + 1] - amps[idx])
om_meas = float(np.pi / np.mean(np.diff(tc)))
d_disc = abs(om_meas - om_disc) / om_disc
d_cont = abs(om_meas - np.sqrt(2.0)) / np.sqrt(2.0)
epsH = Hmax / max(abs(Hrel0), 1e-300)
ok_w = d_disc <= 1e-3
ok_h = epsH <= 1e-3
allpass &= ok_w and ok_h
em("  fala eps=1e-6 k=1 (m=2, N=32): omega_meas=%.6f vs omega_disc="
   "%.6f (|d|=%.2e <=1e-3) -> %s" % (om_meas, om_disc, d_disc,
                                     "OK" if ok_w else "FAIL"))
em("    deskryptywnie: vs kontinuum sqrt(2)=1.414214: |d|=%.2e "
   "(dyspersja siatkowa -- poza gate'em)" % d_cont)
em("  dryf H: max|H-H0|/H_rel(0) = %.3e (H_rel(0)=%.3e) (<=1e-3) -> %s"
   % (epsH, Hrel0, "OK" if ok_h else "FAIL"))
em("  G2: %s" % ("PASS" if (ok_bp and ok_v and ok_w and ok_h)
                 else "FAIL"))

# ------------------------------------------------------------------ G3
em("")
em("G3 (detektory, N=48 L=4pi):")
g3 = Grid3D(48, L_GEN)
r2 = g3.r2_center()
tests = [("dip", 1.0 - 0.4 * np.exp(-r2 / 2.0), 1, 0),
         ("bump", 1.0 + 0.3 * np.exp(-r2 / 2.0), 0, 1),
         ("proznia", np.ones((48, 48, 48)), 0, 0)]
g3ok = True
for name, f, wd, wu in tests:
    nd, sd, nu, su = detect(f)
    ok = (nd == wd and nu == wu)
    g3ok &= ok
    em("  %-8s psi[%.3f,%.3f]: N_dn=%d/%d N_up=%d/%d -> %s"
       % (name, float(f.min()), float(f.max()), nd, wd, nu, wu,
          "OK" if ok else "FAIL"))
allpass &= g3ok
em("  G3: %s" % ("PASS" if g3ok else "FAIL"))

em("")
em("PODSUMOWANIE Phase 1: %s" % ("PASS" if allpass else "FAIL -- STOP"))
with open(BASE + "Phase1_output.txt", "w") as f:
    f.write("\n".join(lines) + "\n")
print("zapisano:", BASE + "Phase1_output.txt")
if not allpass:
    sys.exit(1)
