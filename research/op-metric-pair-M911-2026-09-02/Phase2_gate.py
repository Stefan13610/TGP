#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
op-metric-pair-M911 (Phase 2) -- bramka maszynerii relaksacyjnej.
LOCK sec. 2 Phase 2; MD sec. 8. Testuje DOKLADNIE silnik Phase 3
(import z Phase3_relax_M911 -- ta sama maszyneria, zero duplikacji).

- P2a: proznia psi=psi*=1 (P1a: psi=1 JEST punktem krytycznym pelnego E,
  U'(1)=0 sympy -- Phase1_output.txt) bez zaburzen, t=10, gate
  max_t ||psi-1||_inf <= 1e-10; geometrie: radialna h=0.0125,
  3D L=2pi N=32, 3D L=4pi N=48.
- P2b: detektory z zasianym obiektem dolnym (dip do 0.6 < 5/6)
  i gornym (bump do 1.3 > 7/6) wykrycie 1+-0; czysta proznia zero
  alarmow. Geometrie: 3D L=4pi N=48 i radialna h=0.025 (6 testow).
- Dowolny FAIL ==> STOP (litera LOCKa).

REJESTR [INPUT]: K_geo=gamma=1; progi detektorow psi {5/6, 7/6};
sigma_t=1.0 (pole testowe P2b, MD sec. 8).
"""
import sys
import numpy as np

BASE = ("C:/Users/Mateusz/Documents/ObsydnianMain/TGP/TGP_v1/research/"
        "op-metric-pair-M911-2026-09-02/")
sys.path.insert(0, BASE)
from Phase3_relax_M911 import (FlowRadial, Flow3D, detect, Ueffp,  # noqa: E402
                               PSI_THR_DN, PSI_THR_UP, DT_MAIN,
                               L_LAT, L_GEN, registry_banner)

SIGMA_T = 1.0                                   # INPUT-MD (MD sec. 8)
OUT = BASE + "Phase2_output.txt"
lines = []


def em(s=""):
    print(s, flush=True)
    lines.append(s)


def r2_box(N, L):
    x = (np.arange(N) + 0.5) * (L / N) - L / 2
    return (x[:, None, None] ** 2 + x[None, :, None] ** 2
            + x[None, None, :] ** 2)


em("=" * 78)
em("Phase 2 -- bramka maszynerii (silnik = import Phase3_relax_M911)")
registry_banner("; Phase2 gate")
em("psi* = 1 (P1a sympy: U'(1)=0 -- proznia punktem krytycznym pelnego"
   " E; kontrola float: Ueffp(1.0)=%r)" % Ueffp(1.0))
em("=" * 78)
em()

# ------------------------------------------------------------------ P2a
em("P2a: proznia psi=1 bez zaburzen, t=10, gate max||psi-1|| <= 1e-10")
geoms = [("radial h=0.0125", FlowRadial(0.0125),
          np.ones(int(round(60.0 / 0.0125)))),
         ("3D L=2pi N=32", Flow3D(32, L_LAT), np.ones((32, 32, 32))),
         ("3D L=4pi N=48", Flow3D(48, L_GEN), np.ones((48, 48, 48)))]
p2a_pass = True
for name, flow, g in geoms:
    drift = 0.0
    for k in range(int(round(10.0 / DT_MAIN))):
        g, _ = flow.step(g, DT_MAIN)
        drift = max(drift, float(np.max(np.abs(g - 1.0))))
    ok = drift <= 1e-10
    p2a_pass = p2a_pass and ok
    em("  %-18s max_t ||psi-1||_inf = %.3e -> %s"
       % (name, drift, "OK" if ok else "FAIL"))
em("  P2a: %s" % ("PASS" if p2a_pass else "FAIL"))
em()

# ------------------------------------------------------------------ P2b
em("P2b: detektory (progi FROZEN: dn psi<5/6=%.7f, up psi>7/6=%.7f);"
   % (PSI_THR_DN, PSI_THR_UP))
em("  pola testowe: dip 1-0.4exp(-r^2/2), bump 1+0.3exp(-r^2/2),"
   " proznia (sigma_t=%.1f)" % SIGMA_T)
N = 48
r2b = r2_box(N, L_GEN)
rr = (np.arange(int(round(60.0 / 0.025))) + 0.5) * 0.025
fields = [
    ("3D dip (dolny)", 1.0 - 0.4 * np.exp(-r2b / (2 * SIGMA_T ** 2)), 1, 0),
    ("3D bump (gorny)", 1.0 + 0.3 * np.exp(-r2b / (2 * SIGMA_T ** 2)), 0, 1),
    ("3D proznia", np.ones((N, N, N)), 0, 0),
    ("rad dip (dolny)", 1.0 - 0.4 * np.exp(-rr ** 2 / (2 * SIGMA_T ** 2)),
     1, 0),
    ("rad bump (gorny)", 1.0 + 0.3 * np.exp(-rr ** 2 / (2 * SIGMA_T ** 2)),
     0, 1),
    ("rad proznia", np.ones_like(rr), 0, 0),
]
p2b_pass = True
for name, f, want_dn, want_up in fields:
    n_dn, s_dn, n_up, s_up = detect(f)
    ok = (n_dn == want_dn and n_up == want_up)
    p2b_pass = p2b_pass and ok
    em("  %-18s psi in [%.3f,%.3f]: N_dn=%d (oczek. %d, rozm. %s), "
       "N_up=%d (oczek. %d, rozm. %s) -> %s"
       % (name, float(np.min(f)), float(np.max(f)), n_dn, want_dn,
          s_dn, n_up, want_up, s_up, "OK" if ok else "FAIL"))
em("  P2b: %s" % ("PASS" if p2b_pass else "FAIL"))
em()

verdict = "PASS" if (p2a_pass and p2b_pass) else "FAIL -- STOP"
em("PODSUMOWANIE Phase 2: %s" % verdict)
with open(OUT, "w") as fo:
    fo.write("\n".join(lines) + "\n")
print("zapisano:", OUT)
if not (p2a_pass and p2b_pass):
    sys.exit(1)
