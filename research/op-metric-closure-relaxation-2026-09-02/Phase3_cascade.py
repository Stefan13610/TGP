#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
op-metric-closure-relaxation (Phase 3, wariant kaskadowy) -- charakterystyka
DESKRYPTYWNA nukleacji przy Q-PASS-NUCLEATION (LOCK sec. 3 Phase 3;
MD sec. 7): BEZ progow -- N_obj(t) OBU detektorow, przyrosty na jednostke
czasu, rozklady rozmiarow w t0, t0+5, t0+10, geometria obiektu ze stanu
koncowego. Phase3_spectrum.py NIE UTWORZONY (warunek Q-PASS-STATIC
niespelniony -- litera LOCKa).

Zrodlo danych: Phase2_results/sol_CBAR_f2_*.json/.npz (para nukleacyjna
zbiezna) + deskryptywnie pozostale serie. ZERO nowych progow/kryteriow.
"""
import json
import numpy as np

BASE = ("C:/Users/Mateusz/Documents/ObsydnianMain/TGP/TGP_v1/research/"
        "op-metric-closure-relaxation-2026-09-02/")
RES = BASE + "Phase2_results/"
G_CEIL = float(np.sqrt(4.0 / 3.0))
JOBS = ["sol_CBAR_f2_h025", "sol_CBAR_f2_h0125",
        "sol_CBAR_f2_h025_dt2", "sol_CBAR_f2_h0125_dt2"]
H_OF = {"sol_CBAR_f2_h025": 0.025, "sol_CBAR_f2_h0125": 0.0125,
        "sol_CBAR_f2_h025_dt2": 0.025, "sol_CBAR_f2_h0125_dt2": 0.0125}

print("=" * 78)
print("PHASE 3 (wariant kaskadowy) -- charakterystyka nukleacji sol x C-BAR")
print("REJESTR [INPUT]: g_ceil=%.7f (M9.1'' eq:vol-element-M911); "
      "g_floor=0.5458938; kappa=100; g0_mu=1.4650974" % G_CEIL)
print("=" * 78)
for jid in JOBS:
    with open(RES + jid + ".json") as f:
        r = json.load(f)
    ser = r["series"]
    t0 = r["t_nuc"]
    print("\n--- %s (dt=%.4g): %s-%s, t0=%.1f, N_det=%d ---"
          % (jid, r["dt"], r["status"], r["nuc_dir"], t0, r["n_det"]))
    print("  N_obj(t) [detektor dolny | gorny], rozmiary dolne "
          "(j. komorek):")
    for s in ser:
        print("    t=%5.1f  N_dn=%d  N_up=%d  sz_dn=%-12s sz_up=%-12s "
              "g in [%.4f, %.4f]  E=%.5f"
              % (s["t"], s["n_dn"], s["n_up"], s["sz_dn"], s["sz_up"],
                 s["gmin"], s["gmax"], s["E"]))
    # przyrosty na jednostke czasu (deskryptywnie)
    ndn = [s["n_dn"] for s in ser]
    dif = [ndn[i + 1] - ndn[i] for i in range(len(ndn) - 1)]
    print("  przyrosty N_dn na j.cz.: %s" % dif)
    for toff in (0.0, 5.0, 10.0):
        tt = t0 + toff
        row = [s for s in ser if abs(s["t"] - tt) < 1e-9]
        if row:
            print("  rozklad rozmiarow t0+%g: dolny %s | gorny %s"
                  % (toff, row[0]["sz_dn"], row[0]["sz_up"]))
    # geometria obiektu ze stanu koncowego (deskryptywnie)
    g = np.load(RES + jid + ".npz")["g"]
    h = H_OF[jid]
    rr = (np.arange(len(g)) + 0.5) * h
    m_dn = g < r["gthr_dn"]
    m_up = g > r["gthr_up"]
    if m_dn.any():
        print("  obiekt dolny (stan koncowy): r in [%.3f, %.3f] "
              "(szer. %.3f), min g=%.4f w r=%.3f"
              % (rr[m_dn].min(), rr[m_dn].max(),
                 rr[m_dn].max() - rr[m_dn].min(), g.min(),
                 rr[int(np.argmin(g))]))
    if m_up.any():
        idx = np.where(m_up)[0]
        segs = np.split(idx, np.where(np.diff(idx) > 1)[0] + 1)
        for sgi, sg in enumerate(segs):
            print("  region gorny #%d: r in [%.3f, %.3f], max g=%.4f "
                  "(g_ceil%+.4f)"
                  % (sgi + 1, rr[sg[0]], rr[sg[-1]], g[sg].max(),
                     g[sg].max() - G_CEIL))
print("\nUWAGA deskryptywna (geometria z danych): obiekt dolny = KULA"
      " rdzeniowa r in [0, ~10] -- rdzen solitonu (start g=1.4651,"
      " powyzej g_ceil) INWERTUJE W DOL i osiada na podlodze QB-2"
      " (min g -> 0.549 ~ g_floor+delta); rownolegle objetosc ZEWNETRZNA"
      " (r ponad ~34) wspina sie do studni barierowej C-BAR"
      " g=1.2541 = g_ceil+0.0994 (kappa=100). Kierunek kaskady wg"
      " detektora: DN (jeden obiekt podprogowy, bez dalszego mnozenia"
      " w oknie detekcji); stan w t=12 wciaz ewoluuje (bieg zakonczony"
      " potwierdzeniem okna detektora -- litera LOCKa).")
print("Charakterystyka BEZ progow -- zadnych nowych kryteriow (litera).")
