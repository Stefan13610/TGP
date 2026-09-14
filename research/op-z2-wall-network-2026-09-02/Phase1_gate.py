#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
op-z2-wall-network (Phase 1) -- bramka estymatorow i ciaglosci
(LOCK sec. 2 Phase 1 + Amendment A1). Silnik/estymatory = import
z Phase2_network (ta sama maszyneria).
P1a: paski poprzednika: L_wall = 128.0 DOKLADNIE (2L, sciany osiowe);
     N_dom+=1, N_dom-=1 (sklejanie periodyczne).
P1b: kropla antyfazowa R=8: L_wall(0) = 8R +- 10% (manhattan, A1);
     N_dom-=1; po tau=60 L_wall=0 (ciaglosc: tau_life=52 poprzednika).
P1c: GRF grid-niezalezny: |f128 - interp(f256)| <= 0.05; std zgodne 2%.
P1d: H_Gamma nierosnace (geneza seed=20260906 A=1.0, 2000 krokow).
Dowolny FAIL ==> STOP.
"""
import sys
import numpy as np

BASE = ("C:/Users/Mateusz/Documents/ObsydnianMain/TGP/TGP_v1/research/"
        "op-z2-wall-network-2026-09-02/")
sys.path.insert(0, BASE)
from Phase2_network import (Sub2D, grf, wall_length, label_periodic,  # noqa: E402
                            observables, banner, S_STAR, EPS, L_BOX)

DELTA = 1.07
lines = []


def em(x=""):
    print(x, flush=True)
    lines.append(x)


em("=" * 78)
em("Phase 1 -- bramka estymatorow (silnik = import Phase2_network)")
banner("; Phase1 gate")
em("=" * 78)
allok = True

# ------------------------------------------------------------------ P1a
N, dx = 128, 0.5
x = (np.arange(N) + 0.5) * dx
left = np.tanh((x - L_BOX / 4) / DELTA)
right = np.tanh((3 * L_BOX / 4 - x) / DELTA)
prof = np.where(x < L_BOX / 2, left, right)
s_str = S_STAR * prof[:, None] * np.ones((1, N))
Lw = wall_length(s_str, dx)
ob = observables(s_str, dx)
ok = (Lw == 128.0 and ob["ndp"] == 1 and ob["ndm"] == 1)
allok &= ok
em("P1a paski: L_wall=%.1f (oczek. 128.0) Ndom+=%d Ndom-=%d -> %s"
   % (Lw, ob["ndp"], ob["ndm"], "OK" if ok else "FAIL"))

# ------------------------------------------------------------------ P1b
R = 8.0
xc = (np.arange(N) + 0.5) * dx - L_BOX / 2
r = np.sqrt(xc[:, None] ** 2 + xc[None, :] ** 2)
s_drop = S_STAR * np.tanh((r - R) / DELTA)
Lw0 = wall_length(s_drop, dx)
ob0 = observables(s_drop, dx)
ok0 = abs(Lw0 - 8 * R) / (8 * R) <= 0.10 and ob0["ndm"] == 1
eng = Sub2D(N, dx, 0.02)
s = s_drop.copy()
for k in range(3000):                              # tau=60
    s = eng.step(s)
Lw60 = wall_length(s, dx)
ok60 = (Lw60 == 0.0)
allok &= ok0 and ok60
em("P1b kropla R=8: L_wall(0)=%.1f (oczek. 8R=64 +-10%%) Ndom-=%d; "
   "L_wall(tau=60)=%.1f (oczek. 0) -> %s"
   % (Lw0, ob0["ndm"], Lw60, "OK" if (ok0 and ok60) else "FAIL"))

# ------------------------------------------------------------------ P1c
f128 = grf(20260906, 128)
f256 = grf(20260906, 256)
fi = f256[0::2, 0::2]        # wezly WSPOLPOLOZONE (korekta 1:
# Phase_correction_note_p1c_colocation.md -- ifft2 = rejestr wezlowy)
dmax = float(np.max(np.abs(f128 - fi)))
dstd = abs(float(np.std(f128)) - float(np.std(f256))) \
    / float(np.std(f256))
okc = dmax <= 0.05 and dstd <= 0.02
allok &= okc
em("P1c GRF grid-niezalezny (seed 20260906): max|f128-avg(f256)|=%.4f"
   " (<=0.05); |dstd|=%.4f (<=0.02) -> %s"
   % (dmax, dstd, "OK" if okc else "FAIL"))

# ------------------------------------------------------------------ P1d
s = 1.0 * grf(20260906, 128)
Hp, mono = None, True
for k in range(2000):
    s = eng.step(s)
    if k % 50 == 0:
        Hn = eng.H(s)
        if Hp is not None and Hn > Hp + 1e-10:
            mono = False
        Hp = Hn
allok &= mono
em("P1d H_Gamma nierosnace (geneza A=1.0, 2000 krokow): %s"
   % ("OK" if mono else "FAIL"))

em()
em("PODSUMOWANIE Phase 1: %s" % ("PASS" if allok else "FAIL -- STOP"))
with open(BASE + "Phase1_output.txt", "w") as f:
    f.write("\n".join(lines) + "\n")
print("zapisano:", BASE + "Phase1_output.txt")
if not allok:
    sys.exit(1)
