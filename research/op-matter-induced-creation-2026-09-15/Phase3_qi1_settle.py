#!/usr/bin/env python3
# -*- coding: ascii -*-
"""op-matter-induced-creation -- Phase 3 / Q-I1: faza WLACZONA
(LOCK sec.3; klasyfikatory MD sec.6).

Biegi: start psi=1, pi=0, zrodlo od t=0, lam STALE do t=600
(lam(t) wg rampy LOCKa daje dokladnie lam dla t<=600).
  lista LOCKa : 0.06 0.08 0.10 0.12 0.14 0.16 0.18
  kotwice     : 0.05 0.20
  kontrola Q-I2 (nie wchodzi do werdyktu Q-I1): 0.01
h=0.05, dt=0.005, sponge ON, okno klasyfikacji [500,600].
Potwierdzenie h=0.025: najglebszy SETTLED-SUB (lub najglebszy
SETTLED-DEF, jesli SUB brak).
Bisekcja lam_crit: 6 krokow, protokol IDENTYCZNY (t_max=600).
"""
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import engine_core as ec      # noqa: E402
import runlib as rl           # noqa: E402

OUTP = os.path.join(HERE, "Phase3_qi1_output.txt")
RES = os.path.join(HERE, "Phase3_results")
DT = 0.005
T_ON = 600.0
LIST_LOCK = [0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18]
ANCHORS = [0.05, 0.20]
CONTROL = [0.01]


def tag(lam, h, dt=DT, pre="on"):
    return "%s_lam%.6f_h%.4f_dt%.5f" % (pre, lam, h, dt)


def run_on(lam, h, dt=DT, pre="on", save_state=True, t_max=T_ON):
    """Jeden bieg fazy wlaczonej + klasyfikacja (MD sec.6)."""
    t1 = time.time()
    eng = ec.Engine(h, 200.0, sponge=True, lam=lam)
    g = ec.start_vacuum(eng)
    pi = np.zeros(eng.N)
    tg = tag(lam, h, dt, pre)
    res = rl.evolve(h, dt, lam, 0.0, t_max, g, pi, sponge=True,
                    ramp=True, ckpt_path=os.path.join(RES, tg + "_ckpt.npz"))
    cls = rl.classify_on(res, 500.0, 600.0)
    d = dict(lam=lam, h=h, dt=dt, t_max=t_max, tag=tg,
             status=res["status"], sec=time.time() - t1)
    d.update(rl.jsonable(cls))
    if res["status"] == "OK":
        d["E_core_600"] = float(res["E"][-1])
        d["E_core_field_600"] = float(res["Ef"][-1])
        d["psi0_600"] = float(res["psi0"][-1])
        d["min_psi_600"] = float(res["mn"][-1])
        d["min_psi_run"] = float(np.min(res["mn"]))
        d["max_psi_run"] = float(np.max(res["mx"]))
        if save_state:
            rl.save_state(os.path.join(RES, tg + "_state600.npz"), res,
                          t_max)
        rl.save_series(os.path.join(RES, tg + "_series.npz"), res,
                       extra=dict(psibar=cls["psibar"]))
    else:
        d["min_psi_run"] = float(np.min(res["mn"]))
        d["max_psi_run"] = float(np.max(res["mx"]))
    with open(os.path.join(RES, tg + ".json"), "w") as fh:
        json.dump(d, fh, indent=1)
    return d


def main():
    rl.ensure_dir(RES)
    out = []

    def w(s=""):
        out.append(str(s))
        with open(OUTP, "w") as fh:
            fh.write("\n".join(out) + "\n")

    BAR = "=" * 78
    SEP = "-" * 78
    t_all = time.time()
    w(BAR)
    w("PHASE 3 / Q-I1 -- zrodlo na prozni, faza WLACZONA (LOCK sec.3;"
      " MD sec.6)")
    w("REJESTR [INPUT]: lista lam {0.06,0.08,0.10,0.12,0.14,0.16,0.18}"
      " [LOCK];")
    w("  kotwice {0.05,0.20}; kontrola Q-I2 lam=0.01 (poza werdyktem);")
    w("  h=0.05 (potw. 0.025), dt=0.005, R=200, sponge ON [160,200];")
    w("  start psi=1 pi=0, zrodlo od t=0, t_on=600; okno [500,600];")
    w("  osiadlosc V<=0.01*max(D,1e-12) na r<=80; prog M911 5/6;")
    w("  pas 4/3-1e-6 / 1e-6; bisekcja 6 krokow, t_max=600 [INPUT-MD]")
    w(BAR)

    jobs = sorted(LIST_LOCK + ANCHORS + CONTROL)
    R = {}
    with ProcessPoolExecutor(max_workers=10) as ex:
        futs = {ex.submit(run_on, lam, 0.05): lam for lam in jobs}
        for f in futs:
            d = f.result()
            R[d["lam"]] = d

    w("BIEGI GLOWNE (h=0.05):")
    w("  lam     rola        klasa          psibar(0)   V         D     "
      "    t_end   [s]")
    for lam in jobs:
        d = R[lam]
        rola = ("lista" if lam in LIST_LOCK else
                ("kotwica" if lam in ANCHORS else "kontrola"))
        if d["cls"] == "COLLAPSE":
            w("  %.3f   %-10s  %-14s %-11s %-9s %-9s %.4f  %.0f"
              % (lam, rola, "COLLAPSE(" + d["subtype"] + ")", "n/a",
                 "n/a", "n/a", d["t_end"], d["sec"]))
        else:
            w("  %.3f   %-10s  %-14s %-11.6f %-9.2e %-9.2e %-7s %.0f"
              % (lam, rola, d["cls"], d["psibar0"], d["V"], d["D"],
                 "-", d["sec"]))
    w(SEP)

    # ---- potwierdzenie siatkowe (LOCK sec.3) ---------------------
    subs = [l for l in jobs if R[l]["cls"] == "SETTLED-SUB"]
    defs = [l for l in jobs if R[l]["cls"] == "SETTLED-DEF"]
    if subs:
        lam_conf = min(subs, key=lambda l: R[l]["psibar0"])
        why = "najglebszy SETTLED-SUB"
    elif defs:
        lam_conf = min(defs, key=lambda l: R[l]["psibar0"])
        why = "najglebszy SETTLED-DEF (brak SETTLED-SUB)"
    else:
        lam_conf = None
        why = "brak SETTLED-* -- potwierdzenie siatkowe bezprzedmiotowe"
    w("POTWIERDZENIE SIATKOWE h=0.025 (%s): lam=%s" % (why, lam_conf))

    fut_conf = None
    ex2 = ProcessPoolExecutor(max_workers=2)
    if lam_conf is not None:
        fut_conf = ex2.submit(run_on, lam_conf, 0.025, DT, "onconf")

    # ---- bisekcja lam_crit (6 krokow, FROZEN) --------------------
    cand = sorted(LIST_LOCK + ANCHORS)
    noncol = [l for l in cand if R[l]["cls"] != "COLLAPSE"]
    col = [l for l in cand if R[l]["cls"] == "COLLAPSE"]
    w(SEP)
    if not noncol or not col:
        w("BISEKCJA lam_crit: przedzial startowy nie istnieje"
          " (noncol=%s, col=%s) -- bisekcja bezprzedmiotowa"
          % (noncol, col))
        lam_lo = lam_hi = None
        bis = []
    else:
        lam_lo = max(noncol)
        lam_hi = min(col)
        w("BISEKCJA lam_crit (6 krokow, protokol identyczny t_max=600):")
        w("  przedzial startowy [%.6f, %.6f] (najwyzsze bez COLLAPSE /"
          " najnizsze z COLLAPSE)" % (lam_lo, lam_hi))
        bis = []
        for k in range(1, 7):
            mid = 0.5*(lam_lo + lam_hi)
            d = run_on(mid, 0.05, DT, "bis%d" % k, save_state=False)
            iscol = (d["cls"] == "COLLAPSE")
            bis.append(dict(step=k, lam=mid, cls=d["cls"],
                            subtype=d.get("subtype", ""),
                            t_end=d.get("t_end"),
                            psibar0=d.get("psibar0"), sec=d["sec"]))
            if iscol:
                lam_hi = mid
            else:
                lam_lo = mid
            w("  krok %d: lam=%.7f -> %-24s [lo=%.7f hi=%.7f] (%.0f s)"
              % (k, mid, d["cls"] + ("(" + d["subtype"] + ")" if iscol
                                     else ""), lam_lo, lam_hi, d["sec"]))
        w("  WYNIK deskryptywny: lam_crit = %.7f +- %.7f"
          % (0.5*(lam_lo + lam_hi), 0.5*(lam_hi - lam_lo)))
        w("  (kontekst poprzednika: lam_crit in (0.05, 0.2])")

    d_conf = None
    if fut_conf is not None:
        d_conf = fut_conf.result()
        w(SEP)
        w("POTWIERDZENIE h=0.025 (lam=%.3f): klasa=%s psibar(0)=%s"
          " [%.0f s]"
          % (lam_conf, d_conf["cls"],
             ("%.6f" % d_conf["psibar0"]) if "psibar0" in d_conf
             else "n/a", d_conf["sec"]))
        zg = (d_conf["cls"] == R[lam_conf]["cls"])
        w("  zgodnosc kategorii h=0.05 (%s) vs h=0.025 (%s): %s"
          % (R[lam_conf]["cls"], d_conf["cls"],
             "ZGODNE" if zg else "ROZJAZD"))
    ex2.shutdown()

    # ---- konfrontacja deskryptywna z P1-I2 -----------------------
    w(SEP)
    w("DESKRYPTYWNIE (obowiazkowo): psibar(0) zmierzone vs tabela P1-I2")
    w("  (0D: psi_min(lam) z minimum f=U+lam*psi^2/(4-3psi), rho_hat=1;")
    w("   lam_fold=0.285769734239 -- PRE-REJESTRACJA, bez bramkowania)")
    P1 = {0.01: 0.957824987, 0.05: 0.848184149, 0.06: 0.826611537,
          0.08: 0.786725807, 0.10: 0.749755931, 0.12: 0.714585382,
          0.14: 0.680435158, 0.16: 0.646670971, 0.18: 0.612685373,
          0.20: 0.577786515}
    w("  lam     klasa            psibar(0)_meas  psi_min(0D)  meas-0D")
    for lam in jobs:
        d = R[lam]
        if d["cls"] == "COLLAPSE":
            w("  %.3f   %-16s %-15s %-12.6f %s"
              % (lam, "COLLAPSE", "n/a", P1[lam], "-"))
        else:
            w("  %.3f   %-16s %-15.6f %-12.6f %+.6f"
              % (lam, d["cls"], d["psibar0"], P1[lam],
                 d["psibar0"] - P1[lam]))

    # ---- werdykt Q-I1 (litera LOCK sec.4) ------------------------
    w(SEP)
    subs_list = [l for l in LIST_LOCK if R[l]["cls"] == "SETTLED-SUB"]
    defcol_list = [l for l in LIST_LOCK
                   if R[l]["cls"] in ("SETTLED-DEF", "COLLAPSE")]
    conf_ok = (d_conf is not None and lam_conf is not None
               and d_conf["cls"] == R[lam_conf]["cls"])
    if (subs_list and conf_ok
            and R[lam_conf]["cls"] == "SETTLED-SUB"):
        verdict = "Q-I1-PASS"
        why_v = ("SETTLED-SUB dla lam=%s (lista LOCKa), kategoria"
                 " potwierdzona na h=0.025 dla najglebszego przypadku"
                 " (lam=%.3f)" % (subs_list, lam_conf))
    elif (not subs_list) and len(defcol_list) == len(LIST_LOCK) and conf_ok:
        verdict = "Q-I1-FAIL"
        why_v = ("wszystkie lam z listy: SETTLED-DEF lub COLLAPSE;"
                 " potwierdzenie siatkowe zgodne (lam=%s)" % lam_conf)
    elif (not subs_list) and len(defcol_list) == len(LIST_LOCK) and \
            lam_conf is None:
        verdict = "Q-I1-FAIL"
        why_v = ("wszystkie lam z listy COLLAPSE (brak SETTLED-*;"
                 " potwierdzenie siatkowe bezprzedmiotowe)")
    else:
        verdict = "Q-I1-INCONCLUSIVE"
        why_v = ("brak zbieznego SETTLED-SUB i brak kompletu"
                 " SETTLED-DEF/COLLAPSE na liscie")
    w(BAR)
    w("WERDYKT Q-I1 (litera LOCK sec.4): %s" % verdict)
    w("  uzasadnienie: %s" % why_v)
    w("  SETTLED-SUB na liscie: %s" % (subs_list or "brak"))
    w("  SETTLED-DEF na liscie: %s"
      % ([l for l in LIST_LOCK if R[l]["cls"] == "SETTLED-DEF"] or "brak"))
    w("  COLLAPSE na liscie   : %s"
      % ([l for l in LIST_LOCK if R[l]["cls"] == "COLLAPSE"] or "brak"))
    w(BAR)

    summ = dict(verdict_qi1=verdict, why=why_v,
                lam_conf=lam_conf, conf_ok=bool(conf_ok),
                conf_cls=(d_conf["cls"] if d_conf else None),
                lam_crit=(0.5*(lam_lo + lam_hi) if lam_lo is not None
                          else None),
                lam_crit_halfwidth=(0.5*(lam_hi - lam_lo)
                                    if lam_lo is not None else None),
                bisection=bis,
                runs={("%.3f" % l): rl.jsonable(R[l]) for l in jobs},
                sec=time.time() - t_all)
    with open(os.path.join(RES, "qi1_summary.json"), "w") as fh:
        json.dump(summ, fh, indent=1)
    w("Czas calkowity Q-I1: %.0f s" % (time.time() - t_all))
    w("PHASE3-QI1 DONE")


if __name__ == "__main__":
    main()
    print("PHASE3-QI1 DONE")
