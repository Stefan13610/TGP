#!/usr/bin/env python3
# -*- coding: ascii -*-
"""op-matter-induced-creation -- Phase 3 / Q-I2: WYGASZANIE zrodla
(LOCK sec.3; klasyfikatory MD sec.7).

Dla KAZDEGO biegu SETTLED-* z Q-I1 (lista + kotwice): kontynuacja
ze stanu t=600, rampa lam(t)=lam*S((700-t)/100) (lam=0 od t=700),
dalej ewolucja swobodna do t_max=1700. Klasyfikacja od t=700:
PERSISTENT-OBJECT / RETURN-TO-VACUUM / COLLAPSE / INCONCLUSIVE-RUN.
Potwierdzenia: kazdy PERSISTENT-OBJECT -> h=0.025 ORAZ dt/2 (biegi
od t=0); obowiazkowo jeden RETURN-TO-VACUUM (najnizsze lam SETTLED)
-> h=0.025. Kontrola czystosci wygaszania: lam=0.01 (RETURN
oczekiwany, max_{r<=40}|psi-1| < 1e-3 na koncu).
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

OUTP = os.path.join(HERE, "Phase3_qi2_output.txt")
RES = os.path.join(HERE, "Phase3_results")
DT = 0.005
T_ON = 600.0
T_MAX = 1700.0
T0 = 2.0*np.pi
LIST_LOCK = [0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18]
ANCHORS = [0.05, 0.20]
CONTROL = 0.01


def stag(lam, h, dt, pre):
    return "%s_lam%.6f_h%.4f_dt%.5f" % (pre, lam, h, dt)


def run_off_from_ckpt(lam):
    """Kontynuacja z checkpointu t=600 (h=0.05, dt=0.005)."""
    t1 = time.time()
    src = os.path.join(RES, "on_lam%.6f_h%.4f_dt%.5f_state600.npz"
                       % (lam, 0.05, DT))
    z = np.load(src)
    g = z["g"].copy()
    pi = z["pi"].copy()
    tg = stag(lam, 0.05, DT, "off")
    res = rl.evolve(0.05, DT, lam, T_ON, T_MAX, g, pi, sponge=True,
                    ramp=True,
                    ckpt_path=os.path.join(RES, tg + "_ckpt.npz"))
    cls = rl.classify_off(res)
    d = dict(lam=lam, h=0.05, dt=DT, tag=tg, mode="from_ckpt600",
             status=res["status"], sec=time.time() - t1)
    d.update(rl.jsonable(cls))
    d["E_core_600"] = float(res["E"][0])
    if res["status"] == "OK":
        d["min_psi_off"] = float(np.min(res["mn"][res["t"] >= 700.0]))
        d["max_psi_off"] = float(np.max(res["mx"][res["t"] >= 700.0]))
        d["psi0_1700"] = float(res["psi0"][-1])
    rl.save_series(os.path.join(RES, tg + "_series.npz"), res)
    with open(os.path.join(RES, tg + ".json"), "w") as fh:
        json.dump(d, fh, indent=1)
    return d


def run_full(lam, h, dt, pre):
    """Bieg potwierdzajacy od t=0 do t=1700 (rampa wg LOCKa)."""
    t1 = time.time()
    eng = ec.Engine(h, 200.0, sponge=True, lam=lam)
    g = ec.start_vacuum(eng)
    pi = np.zeros(eng.N)
    tg = stag(lam, h, dt, pre)
    res = rl.evolve(h, dt, lam, 0.0, T_MAX, g, pi, sponge=True,
                    ramp=True,
                    ckpt_path=os.path.join(RES, tg + "_ckpt.npz"))
    cls_on = rl.classify_on(res, 500.0, 600.0)
    cls = rl.classify_off(res)
    d = dict(lam=lam, h=h, dt=dt, tag=tg, mode="full_from_0",
             status=res["status"], sec=time.time() - t1,
             cls_on=cls_on["cls"], psibar0_on=cls_on.get("psibar0"))
    d.update(rl.jsonable(cls))
    if res["status"] == "OK":
        d["psi0_1700"] = float(res["psi0"][-1])
    rl.save_series(os.path.join(RES, tg + "_series.npz"), res)
    with open(os.path.join(RES, tg + ".json"), "w") as fh:
        json.dump(d, fh, indent=1)
    return d


def fmt_off(d):
    if d["cls"] == "COLLAPSE":
        return "COLLAPSE(%s) t_end=%.4f" % (d["subtype"], d["t_end"])
    return d["cls"]


def main():
    rl.ensure_dir(RES)
    with open(os.path.join(RES, "qi1_summary.json")) as fh:
        Q1 = json.load(fh)
    out = []

    def w(s=""):
        out.append(str(s))
        with open(OUTP, "w") as fh:
            fh.write("\n".join(out) + "\n")

    BAR = "=" * 78
    SEP = "-" * 78
    t_all = time.time()
    w(BAR)
    w("PHASE 3 / Q-I2 -- WYGASZANIE zrodla (LOCK sec.3; MD sec.7)")
    w("REJESTR [INPUT]: rampa lam(t)=lam*S((600+100-t)/100), Delta=100,")
    w("  t_off=600 => lam=0 od t=700 [LOCK]; t_max=1700; h=0.05,")
    w("  dt=0.005, sponge ON; E_ref_off=E_core(700), E_core r<=80;")
    w("  okno trwalosci W_P=[700, 700+100T0]=[700, %.4f];" % (700+100*T0))
    w("  okno koncowe W_F=[1600,1700] [INPUT-MD]; progi: 0.5*E_ref_off,")
    w("  max_{r<=40}|psi-1|>=0.02 (PERSISTENT-OBJECT); E(1700)<0.05*")
    w("  E_ref_off lub max|psi-1|<1e-3 w W_F (RETURN-TO-VACUUM);")
    w("  potwierdzenia: PERSISTENT -> h=0.025 + dt/2; jeden RETURN")
    w("  (najnizsze lam SETTLED) -> h=0.025; kontrola lam=0.01")
    w(BAR)
    w("PREDYKCJA PRE-REJESTROWANA (LOCK sec.0, NIENARUSZALNA):")
    w("  \"w galezi zdrowej BEZ zrodla nie ma statycznych solitonow")
    w("   (Q-B-FAIL) ani oscylonow (Q-E/Q-G) ==> oczekiwany wynik Q-I2:")
    w("   RETURN-TO-VACUUM lub COLLAPSE\"")
    w(BAR)

    runs = Q1["runs"]
    settled = sorted([float(k) for k, v in runs.items()
                      if v["cls"].startswith("SETTLED")
                      and float(k) in LIST_LOCK + ANCHORS])
    w("Biegi SETTLED-* z Q-I1 (lista+kotwice), do wygaszenia: %s"
      % (settled or "BRAK"))
    ctrl_settled = (runs.get("%.3f" % CONTROL, {}).get("cls", "")
                    .startswith("SETTLED"))
    w("Kontrola czystosci wygaszania lam=0.01: klasa Q-I1 = %s"
      % runs.get("%.3f" % CONTROL, {}).get("cls", "BRAK"))
    w(SEP)

    jobs = list(settled) + ([CONTROL] if ctrl_settled else [])
    R = {}
    with ProcessPoolExecutor(max_workers=10) as ex:
        futs = [ex.submit(run_off_from_ckpt, lam) for lam in jobs]
        for f in futs:
            d = f.result()
            R[d["lam"]] = d

    w("BIEGI WYGASZANIA (h=0.05, dt=0.005, kontynuacja z t=600):")
    w("  lam     rola      klasa po rampie      E(700)      E(1700)"
      "     E(1700)/E(700)  max|psi-1|(r<=40,t=1700)  tau     [s]")
    for lam in jobs:
        d = R[lam]
        rola = ("lista" if lam in LIST_LOCK else
                ("kotwica" if lam in ANCHORS else "kontrola"))
        if d["cls"] == "COLLAPSE":
            w("  %.3f   %-8s  %-20s %-11s %-11s %-15s %-24s %-7s %.0f"
              % (lam, rola, fmt_off(d), "n/a", "n/a", "n/a", "n/a",
                 "n/a", d["sec"]))
        else:
            w("  %.3f   %-8s  %-20s %-+11.4e %-+11.4e %-15.4e %-24.4e"
              " %-7.1f %.0f"
              % (lam, rola, d["cls"], d["Eref700"], d["E_end"],
                 d["E_ratio"], d["D_end"], d["tau"], d["sec"]))
    w(SEP)
    w("BILANS ENERGII (deskryptywnie; podczas rampy energia NIE jest")
    w("zachowana -- praca zrodla; LOCK: nie bramkowac):")
    w("  lam     E_core(600)   E_core(700)   E_core(1700)  dE_rampa")
    for lam in jobs:
        d = R[lam]
        if d["cls"] == "COLLAPSE":
            w("  %.3f   %+.5e  n/a           n/a           n/a"
              % (lam, d["E_core_600"]))
        else:
            w("  %.3f   %+.5e  %+.5e  %+.5e  %+.5e"
              % (lam, d["E_core_600"], d["Eref700"], d["E_end"],
                 d["Eref700"] - d["E_core_600"]))
    w(SEP)

    # ---- kontrola czystosci wygaszania (LOCK sec.3) --------------
    if CONTROL in R:
        d = R[CONTROL]
        pure_ok = (d["cls"] == "RETURN-TO-VACUUM"
                   and d.get("D_end", 9e9) < 1e-3)
        w("KONTROLA CZYSTOSCI WYGASZANIA (lam=0.01, identyczna rampa):")
        w("  klasa = %s (oczekiwane RETURN-TO-VACUUM)" % d["cls"])
        w("  max_{r<=40}|psi-1| (t=1700) = %.3e  (prog 1e-3)  %s"
          % (d.get("D_end", float("nan")),
             "PASS" if d.get("D_end", 9e9) < 1e-3 else "FAIL"))
        w("  GATE maszynerii wygaszania: %s"
          % ("PASS" if pure_ok else "FAIL"))
    else:
        pure_ok = False
        w("KONTROLA CZYSTOSCI WYGASZANIA: bieg lam=0.01 nie byl"
          " SETTLED-* w Q-I1 -- kontrola bezprzedmiotowa (FAIL formalny)")
    w(SEP)

    # ---- potwierdzenia (LOCK sec.3, FROZEN) ----------------------
    pers = [l for l in settled if R[l]["cls"] == "PERSISTENT-OBJECT"]
    rets = [l for l in settled if R[l]["cls"] == "RETURN-TO-VACUUM"]
    conf_jobs = []
    for l in pers:
        conf_jobs.append((l, 0.025, DT, "conf_h025"))
        conf_jobs.append((l, 0.05, DT/2.0, "conf_dt2"))
    lam_ret = min(rets) if rets else None
    if lam_ret is not None:
        conf_jobs.append((lam_ret, 0.025, DT, "conf_h025"))
    w("POTWIERDZENIA (LOCK sec.3):")
    w("  PERSISTENT-OBJECT do potwierdzenia (h=0.025 + dt/2): %s"
      % (pers or "brak"))
    w("  obowiazkowa kontrola negatywu RETURN-TO-VACUUM (najnizsze lam"
      " SETTLED): %s" % (lam_ret if lam_ret is not None else
                         "brak RETURN-TO-VACUUM"))
    C = {}
    if conf_jobs:
        with ProcessPoolExecutor(max_workers=max(1, len(conf_jobs))) as ex:
            futs = [ex.submit(run_full, l, h, dt, pre)
                    for (l, h, dt, pre) in conf_jobs]
            for f in futs:
                d = f.result()
                C[(d["lam"], d["h"], d["dt"])] = d
        w("  bieg potwierdzajacy        klasa_on   klasa_off"
          "            tau      E(1700)/E(700)  [s]")
        for (l, h, dt, pre) in conf_jobs:
            d = C[(l, h, dt)]
            w("  lam=%.3f h=%.4f dt=%.5f  %-10s %-20s %-8s %-15s %.0f"
              % (l, h, dt, d.get("cls_on", "n/a"), fmt_off(d),
                 ("%.1f" % d["tau"]) if "tau" in d else "n/a",
                 ("%.4e" % d["E_ratio"]) if "E_ratio" in d else "n/a",
                 d["sec"]))
    else:
        w("  (brak biegow potwierdzajacych)")

    def conf_agree(lam_main, h, dt):
        d0 = R[lam_main]
        d1 = C.get((lam_main, h, dt))
        if d1 is None:
            return False, "brak biegu"
        if d1["cls"] != d0["cls"]:
            return False, "kategoria %s vs %s" % (d0["cls"], d1["cls"])
        if d0["cls"] == "PERSISTENT-OBJECT":
            t0, t1_ = d0["tau"], d1["tau"]
            cens0 = d0.get("tau_censored"),
            if d0.get("tau_censored") and d1.get("tau_censored"):
                return True, "obie cenzurowane (tau>=1000)"
            if d0.get("tau_censored"):
                return (t1_ >= 900.0), "tau_conf=%.1f (prog >=900)" % t1_
            ok = abs(t1_ - t0) <= 0.10*t0
            return ok, "tau %.1f vs %.1f (+-10%%)" % (t0, t1_)
        if d0["cls"] == "COLLAPSE":
            ok = abs(d1["t_end"] - d0["t_end"]) <= 0.10*d0["t_end"]
            return ok, "t_end %.4f vs %.4f" % (d0["t_end"], d1["t_end"])
        return True, "kategoria zgodna"

    w(SEP)
    pers_conf = []
    for l in pers:
        a1, m1 = conf_agree(l, 0.025, DT)
        a2, m2 = conf_agree(l, 0.05, DT/2.0)
        w("  lam=%.3f PERSISTENT: h=0.025 %s (%s); dt/2 %s (%s)"
          % (l, "ZGODNE" if a1 else "ROZJAZD", m1,
             "ZGODNE" if a2 else "ROZJAZD", m2))
        if a1 and a2:
            pers_conf.append(l)
    ret_conf = None
    if lam_ret is not None:
        a1, m1 = conf_agree(lam_ret, 0.025, DT)
        ret_conf = a1
        w("  lam=%.3f RETURN (kontrola negatywu): h=0.025 %s (%s)"
          % (lam_ret, "ZGODNE" if a1 else "ROZJAZD", m1))

    # ---- werdykt Q-I2 (litera LOCK sec.4) ------------------------
    cls_all = [R[l]["cls"] for l in settled]
    all_neg = all(c in ("RETURN-TO-VACUUM", "COLLAPSE") for c in cls_all) \
        and len(cls_all) > 0
    if pers_conf:
        verdict = "Q-I2-PASS"
        why = ("PERSISTENT-OBJECT potwierdzony (h/2 i dt/2) dla lam=%s"
               % pers_conf)
    elif all_neg and (ret_conf is True or lam_ret is None):
        verdict = "Q-I2-FAIL"
        why = ("wszystkie biegi SETTLED-* po wygaszeniu:"
               " RETURN-TO-VACUUM lub COLLAPSE; kontrola negatywu"
               " (lam=%s, h=0.025) %s"
               % (lam_ret, "zgodna" if ret_conf else "bezprzedmiotowa"))
    elif not settled:
        verdict = "Q-I2-INCONCLUSIVE"
        why = "brak biegow SETTLED-* do wygaszenia"
    else:
        verdict = "Q-I2-INCONCLUSIVE"
        why = ("brak potwierdzonego PERSISTENT-OBJECT i brak kompletu"
               " RETURN/COLLAPSE zbieznie (klasy: %s)" % cls_all)
    w(BAR)
    w("WERDYKT Q-I2 (litera LOCK sec.4): %s" % verdict)
    w("  uzasadnienie: %s" % why)
    w("  PERSISTENT-OBJECT: %s ; RETURN-TO-VACUUM: %s ; COLLAPSE: %s ;"
      " INCONCLUSIVE-RUN: %s"
      % ([l for l in settled if R[l]["cls"] == "PERSISTENT-OBJECT"] or "brak",
         [l for l in settled if R[l]["cls"] == "RETURN-TO-VACUUM"] or "brak",
         [l for l in settled if R[l]["cls"] == "COLLAPSE"] or "brak",
         [l for l in settled if R[l]["cls"] == "INCONCLUSIVE-RUN"] or "brak"))
    w(SEP)
    w("KONFRONTACJA Z PREDYKCJA PRE-REJESTROWANA (bez reinterpretacji):")
    if verdict == "Q-I2-FAIL":
        w("  predykcja \"RETURN-TO-VACUUM lub COLLAPSE\": TRAFIONA")
    elif verdict == "Q-I2-PASS":
        w("  predykcja \"RETURN-TO-VACUUM lub COLLAPSE\": OBALONA")
        w("  (pozytyw WBREW predykcji -- litera rozstrzyga, LOCK sec.4)")
    else:
        w("  predykcja \"RETURN-TO-VACUUM lub COLLAPSE\": NIEROZSTRZYGNIETA")
        w("  (INCONCLUSIVE != pozytyw)")
    w(BAR)

    summ = dict(verdict_qi2=verdict, why=why,
                settled=settled, classes={("%.3f" % l): R[l]["cls"]
                                          for l in jobs},
                purity_ok=bool(pure_ok), lam_ret=lam_ret,
                ret_conf=ret_conf, pers=pers, pers_conf=pers_conf,
                runs={("%.3f" % l): rl.jsonable(R[l]) for l in jobs},
                confirmations={("%.3f_%.4f_%.5f" % k): rl.jsonable(v)
                               for k, v in C.items()},
                sec=time.time() - t_all)
    with open(os.path.join(RES, "qi2_summary.json"), "w") as fh:
        json.dump(summ, fh, indent=1)
    vj = dict(verdict_qi1=Q1["verdict_qi1"], verdict_qi2=verdict,
              lam_crit=Q1.get("lam_crit"),
              lam_crit_halfwidth=Q1.get("lam_crit_halfwidth"),
              qi1_runs={k: v["cls"] for k, v in Q1["runs"].items()},
              qi2_classes={("%.3f" % l): R[l]["cls"] for l in jobs},
              purity_ok=bool(pure_ok))
    with open(os.path.join(RES, "verdict.json"), "w") as fh:
        json.dump(vj, fh, indent=1)
    w("Czas calkowity Q-I2: %.0f s" % (time.time() - t_all))
    w("PHASE3-QI2 DONE")


if __name__ == "__main__":
    main()
    print("PHASE3-QI2 DONE")
