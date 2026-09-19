#!/usr/bin/env python3
# -*- coding: ascii -*-
"""op-matter-induced-creation -- PROBE DESKRYPTYWNY (poza protokolem
werdyktowym; uruchamiany PO zapisaniu werdyktu Q-I1).

Cel: charakterystyka ROZJAZDU siatek stwierdzonego w Q-I1
(lam=0.10: SETTLED-SUB na h=0.05, COLLAPSE t=0.95 na h=0.025).
Pytanie deskryptywne (material do NEEDS metodologicznego, NIE do
werdyktu): czy rozjazd dotyczy tylko okolicy lam_crit, i czy jest
efektem siatki czy kroku czasowego?

Biegi (wszystkie faza WLACZONA, t=600, klasyfikacja MD sec.6):
  lam=0.06 h=0.025 ; lam=0.08 h=0.025 ; lam=0.10 h=0.05 dt/2.
WERDYKT Q-I1 POZOSTAJE BEZ ZMIAN (LOCK sec.3 przewiduje potwierdzenie
h=0.025 WYLACZNIE dla najglebszego SETTLED-SUB = lam=0.10; ten
rozjechal sie i litera LOCK sec.4 daje Q-I1-INCONCLUSIVE).
"""
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import runlib as rl                      # noqa: E402
from Phase3_qi1_settle import run_on     # noqa: E402

OUTP = os.path.join(HERE, "Phase3_desc_gridconv_output.txt")
RES = os.path.join(HERE, "Phase3_results")
JOBS = [(0.06, 0.025, 0.005), (0.08, 0.025, 0.005),
        (0.10, 0.05, 0.0025)]


def main():
    rl.ensure_dir(RES)
    out = []

    def w(s=""):
        out.append(str(s))
        with open(OUTP, "w") as fh:
            fh.write("\n".join(out) + "\n")

    BAR = "=" * 78
    t0 = time.time()
    w(BAR)
    w("PROBE DESKRYPTYWNY -- zbieznosc siatkowa w okolicy lam_crit")
    w("(POZA protokolem werdyktowym; werdykt Q-I1 juz zapisany")
    w(" i NIEZMIENIONY: Q-I1-INCONCLUSIVE, litera LOCK sec.4)")
    w("Kontekst: lam=0.10 SETTLED-SUB (h=0.05) vs COLLAPSE t=0.95")
    w("  (h=0.025) -- rozjazd kategorii na potwierdzeniu LOCKowym.")
    w(BAR)
    R = {}
    with ProcessPoolExecutor(max_workers=3) as ex:
        futs = [ex.submit(run_on, lam, h, dt, "desc") for lam, h, dt in JOBS]
        for f in futs:
            d = f.result()
            R[(d["lam"], d["h"], d["dt"])] = d
    w("  lam    h       dt       klasa                    psibar(0)"
      "   t_end   [s]")
    for lam, h, dt in JOBS:
        d = R[(lam, h, dt)]
        if d["cls"] == "COLLAPSE":
            w("  %.3f  %.4f  %.5f  %-24s %-11s %-7.4f %.0f"
              % (lam, h, dt, "COLLAPSE(" + d["subtype"] + ")", "n/a",
                 d["t_end"], d["sec"]))
        else:
            w("  %.3f  %.4f  %.5f  %-24s %-11.6f %-7s %.0f"
              % (lam, h, dt, d["cls"], d["psibar0"], "-", d["sec"]))
    w("")
    w("Odniesienie (h=0.05, dt=0.005, Phase3_qi1_output.txt):")
    w("  lam=0.06 SETTLED-DEF psibar(0)=0.845733")
    w("  lam=0.08 SETTLED-SUB psibar(0)=0.808031")
    w("  lam=0.10 SETTLED-SUB psibar(0)=0.772903")
    w(BAR)
    with open(os.path.join(RES, "desc_gridconv.json"), "w") as fh:
        json.dump({("%.3f_%.4f_%.5f" % k): rl.jsonable(v)
                   for k, v in R.items()}, fh, indent=1)
    w("Czas: %.0f s" % (time.time() - t0))
    w("PROBE-GRIDCONV DONE")


if __name__ == "__main__":
    main()
    print("PROBE-GRIDCONV DONE")
