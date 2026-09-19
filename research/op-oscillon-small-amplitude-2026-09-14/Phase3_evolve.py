#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
op-oscillon-small-amplitude -- Phase 3: Q-G (LOCK sec. 3, 5;
MD sec. 3-5). Driver biegow produkcyjnych z checkpointami npz
(co 500 j.cz.) i wznawianiem.

Uzycie:
  python Phase3_evolve.py run <name> <grid> <t_target>
      name: g_a{a}_s{s} (np. g_a0.05_s6) lub vac
      grid: h05 (h=0.05, dt=0.005) | h025 (h=0.025, dt=0.005)
            | h05_dt2 (h=0.05, dt=0.0025)
  python Phase3_evolve.py triage      -> tabela triage Etapu A (FROZEN)
  python Phase3_evolve.py report      -> klasyfikacja + werdykt Q-G
                                         (Phase3_output.txt, verdict.json)

Protokol FROZEN (MD sec. 3-4): R=400; sponge [320,400] gamma0=1.0;
dt_out=0.1; profil co 100 j.cz.; checkpoint co 500 j.cz.; E_core r<=80;
E_ref=E_core(50); Etap A t=2000 -> triage E_core(2000)>=0.2*E_ref
i zero zdarzen brzegowych -> Etap B t_max=10000.
"""
import json
import os
import sys
import time
import numpy as np

BASE = ("C:/Users/Mateusz/Documents/ObsydnianMain/TGP/TGP_v1/research/"
        "op-oscillon-small-amplitude-2026-09-14/")
sys.path.insert(0, BASE)
from engine_core import (Engine, start_gauss, start_vacuum,
                         NonConvergence)

RES = BASE + "Phase3_results/"
R_BOX = 400.0
DT_OUT = 0.1
T0 = 2.0*np.pi
HOLD_MIN = 100.0*T0          # 628.319 [LOCK]
T_STAGE_A = 2000.0
T_MAX = 10000.0
A_LIST = (0.02, 0.05, 0.08, 0.10)
S_LIST = (3.0, 6.0, 10.0)
STARTS = ["g_a%.2f_s%g" % (a, s) for a in A_LIST for s in S_LIST]
GRIDS = {"h05": (0.05, 0.005), "h025": (0.025, 0.005),
         "h05_dt2": (0.05, 0.0025)}


def rdir(name, grid):
    return RES + "%s__%s/" % (name, grid)


def make_start(eng, name):
    if name == "vac":
        return start_vacuum(eng)
    part = name.split("_")
    a = float(part[1][1:])
    s = float(part[2][1:])
    return start_gauss(eng, a, s)


def load_meta(name, grid):
    p = rdir(name, grid) + "meta.json"
    if not os.path.exists(p):
        return None
    with open(p) as f:
        return json.load(f)


def save_meta(name, grid, meta):
    with open(rdir(name, grid) + "meta.json", "w") as f:
        json.dump(meta, f, indent=1)


def log(msg):
    with open(BASE + "Phase3_progress.log", "a") as f:
        f.write("[%s] %s\n" % (time.strftime("%H:%M:%S"), msg))


# ------------------------------------------------------------- run --
def do_run(name, grid, t_target):
    h, dt = GRIDS[grid]
    d = rdir(name, grid)
    os.makedirs(d, exist_ok=True)
    n_out = int(round(DT_OUT/dt))          # probka co dt_out
    n_prof = int(round(100.0/dt))          # profil co 100 j.cz.
    n_ckpt = int(round(500.0/dt))          # checkpoint co 500 j.cz.
    eng = Engine(h, R_BOX, sponge=True)
    meta = load_meta(name, grid)
    if meta is None or not os.path.exists(d + "state.npz"):
        g = make_start(eng, name)
        pi = np.zeros(eng.N)
        k = 0
        psi0 = [float(g[0])]
        ecore = [eng.energy(g, pi, rmax=80.0)]
        prof_t = [0.0]
        profs = [g.copy()]
        meta = {"name": name, "grid": grid, "h": h, "dt": dt,
                "status": "RUNNING", "event": None, "t_event": None,
                "t_reached": 0.0}
    else:
        if meta["status"] != "RUNNING" or meta["t_reached"] >= t_target:
            print("SKIP %s %s status=%s t=%.1f"
                  % (name, grid, meta["status"], meta["t_reached"]))
            return
        st = np.load(d + "state.npz")
        g = st["g"]
        pi = st["pi"]
        k = int(st["k"])
        psi0 = list(st["psi0"])
        ecore = list(st["ecore"])
        prof_t = list(st["prof_t"])
        profs = [p for p in st["profs"]]

    def save_state():
        np.savez_compressed(
            d + "state.npz", g=g, pi=pi, k=k,
            psi0=np.array(psi0), ecore=np.array(ecore),
            prof_t=np.array(prof_t), profs=np.array(profs))
        meta["t_reached"] = k*dt
        save_meta(name, grid, meta)

    k_target = int(round(t_target/dt))
    t_wall = time.time()
    log("START %s %s k=%d -> %d" % (name, grid, k, k_target))
    while k < k_target:
        try:
            g, pi = eng.step(g, pi, dt)
        except NonConvergence:
            meta["status"] = "EVENT"
            meta["event"] = "NONFINITE"
            meta["t_event"] = (k + 1)*dt
            save_state()
            log("EVENT %s %s NONCONV t=%.3f" % (name, grid, (k+1)*dt))
            return
        k += 1
        bs = eng.band_status(g)
        if bs is not None:
            meta["status"] = "EVENT"
            meta["event"] = {"BREAKDOWN": "NONFINITE",
                             "BREAKDOWN-BOUNDARY": "BOUNDARY-UPPER",
                             "BREAKDOWN-BOUNDARY-LOWER":
                                 "BOUNDARY-LOWER"}[bs]
            meta["t_event"] = k*dt
            save_state()
            log("EVENT %s %s %s t=%.3f" % (name, grid, bs, k*dt))
            return
        if k % n_out == 0:
            psi0.append(float(g[0]))
            ecore.append(eng.energy(g, pi, rmax=80.0))
        if k % n_prof == 0:
            prof_t.append(k*dt)
            profs.append(g.copy())
        if k % n_ckpt == 0:
            save_state()
            dev = float(np.max(np.abs(g - 1.0)))
            cmax = float(np.max((4.0 - 3.0*g)/g))
            log("CKPT %s %s t=%.0f Ecore=%.4e max|psi-1|=%.3e "
                "maxc=%.3f wall=%.0fs"
                % (name, grid, k*dt, ecore[-1], dev, cmax,
                   time.time() - t_wall))
    if k*dt >= T_MAX - 1e-9:
        meta["status"] = "DONE"
    save_state()
    log("END %s %s t=%.0f wall=%.0fs" % (name, grid, k*dt,
                                         time.time() - t_wall))
    print("RUN DONE %s %s t=%.0f" % (name, grid, k*dt))


# ------------------------------------------------- analiza pomocna --
def series(name, grid):
    st = np.load(rdir(name, grid) + "state.npz")
    psi0 = st["psi0"]
    ecore = st["ecore"]
    t = np.arange(len(psi0))*DT_OUT
    return t, psi0, ecore


def hold_windows(t, ecore):
    """Maksymalne przedzialy [50,t_end] z E_core >= 0.5 E_ref.
    Zwraca (E_ref, lista (t1,t2), tau)."""
    i50 = int(round(50.0/DT_OUT))
    if i50 >= len(t):
        return None, [], 0.0
    E_ref = float(ecore[i50])
    m = ecore[i50:] >= 0.5*E_ref
    tt = t[i50:]
    wins = []
    start = None
    for i in range(len(m)):
        if m[i] and start is None:
            start = i
        if (not m[i]) and start is not None:
            wins.append((tt[start], tt[i-1]))
            start = None
    if start is not None:
        wins.append((tt[start], tt[-1]))
    # tau (konwencja poprzednika, MD sec. 4): koniec okna od t=50
    tau = 0.0
    if wins and abs(wins[0][0] - 50.0) < 1e-9:
        tau = wins[0][1] - 50.0
    return E_ref, wins, tau


def fft_peak(t, psi0, t1, t2):
    """FFT psi(0,t)-srednia, okno Hanna, segment [t1,t2];
    pik paraboliczny; zwraca (omega, domega, harmoniki, moc)."""
    m = (t >= t1 - 1e-9) & (t <= t2 + 1e-9)
    x = psi0[m] - np.mean(psi0[m])
    n = len(x)
    if n < 8:
        return None
    w = np.hanning(n)
    P = np.abs(np.fft.rfft(x*w))**2
    om = 2.0*np.pi*np.fft.rfftfreq(n, d=DT_OUT)
    dom_bin = om[1] - om[0]
    sel = om >= 0.1          # odcina DC/najnizsze biny
    if not np.any(sel):
        return None
    i0 = np.where(sel)[0][0]
    ip = i0 + int(np.argmax(P[i0:]))
    if not (P[ip] > 0):
        return None   # sygnal zerowy: brak piku (correction note 1a)
    # interpolacja paraboliczna (log-moc) [INPUT-MD]
    if 0 < ip < len(P) - 1 and P[ip-1] > 0 and P[ip+1] > 0:
        la, lb, lc = np.log(P[ip-1]), np.log(P[ip]), np.log(P[ip+1])
        dd = la - 2*lb + lc
        delta = 0.5*(la - lc)/dd if dd != 0 else 0.0
    else:
        delta = 0.0
    om_pk = om[ip] + delta*dom_bin
    # harmoniki: maksima lokalne o mocy >= 1% dominujacego
    harm = []
    for i in range(i0 + 1, len(P) - 1):
        if P[i] > P[i-1] and P[i] > P[i+1] and P[i] >= 0.01*P[ip] \
                and i != ip:
            harm.append(round(float(om[i]), 4))
    return {"omega": float(om_pk), "domega": float(dom_bin),
            "harmonics": harm[:8], "power": float(P[ip])}


def crossings(t, psi0, t1, t2):
    m = (t >= t1 - 1e-9) & (t <= t2 + 1e-9)
    s = np.sign(psi0[m] - 1.0)
    s = s[s != 0]
    return int(np.sum(s[1:] != s[:-1]))


def analyze_run(name, grid):
    """Pelna analiza pojedynczego biegu -> dict."""
    meta = load_meta(name, grid)
    if meta is None:
        return None
    t, psi0, ecore = series(name, grid)
    t_end = meta["t_reached"]
    out = {"name": name, "grid": grid, "status": meta["status"],
           "event": meta["event"], "t_event": meta["t_event"],
           "t_end": t_end}
    if meta["status"] == "EVENT":
        return out
    E_ref, wins, tau = hold_windows(t, ecore)
    out["E_ref"] = E_ref
    out["tau"] = tau
    out["censored"] = bool(wins and abs(wins[0][0] - 50.0) < 1e-9
                           and abs(wins[0][1] - t_end) < DT_OUT
                           and t_end >= T_MAX - 1e-9)
    out["E_end_over_ref"] = (float(ecore[-1])/E_ref
                             if E_ref and E_ref > 0 else None)
    # najdluzsze okno podtrzymania (moze zaczynac sie > 50)
    best = max(wins, key=lambda w: w[1] - w[0]) if wins else None
    cond1 = bool(best and (best[1] - best[0]) >= HOLD_MIN)
    out["hold_window"] = list(best) if best else None
    out["cond1"] = cond1
    out["ncross"] = crossings(t, psi0, best[0], best[1]) if best else 0
    out["cond2"] = out["ncross"] >= 50
    # segment FFT: T = koniec okna podtrzymania (jesli jest) albo t_end
    T = best[1] if best else t_end
    t1 = max(50.0, T - 2000.0)
    out["fft_segment"] = [t1, T]
    fp = fft_peak(t, psi0, t1, T) if (T - t1) >= 1000.0 else None
    out["fft"] = fp
    out["cond3"] = bool(fp and fp["omega"] <= 0.99)
    # omega_desc: deskryptywnie na CALYM oknie podtrzymania (>=200
    # j.cz.); NIE wchodzi do detektora (correction note 1b)
    out["fft_desc"] = (fft_peak(t, psi0, best[0], best[1])
                       if best and (best[1] - best[0]) >= 200.0
                       else None)
    out["candidate"] = cond1 and out["cond2"] and out["cond3"]
    # quasi-stacjonarnosc (OSCILLON-WEAK; MD sec. 4)
    qs = None
    if t_end >= 1000.0 + 50.0 and E_ref and E_ref > 0:
        m1 = (t >= t_end - 1000.0) & (t <= t_end - 900.0)
        m2 = (t >= t_end - 100.0)
        e1 = float(np.mean(ecore[m1]))
        e2 = float(np.mean(ecore[m2]))
        qs = abs(e2 - e1)/abs(e2) if e2 != 0 else None
    out["quasi_stat"] = qs
    out["weak_flag"] = bool(
        (not out["candidate"]) and out["E_end_over_ref"] is not None
        and out["E_end_over_ref"] < 0.5 and qs is not None
        and qs <= 1e-3 and out["cond3"])
    out["radiated_flag"] = bool(
        out["E_end_over_ref"] is not None
        and out["E_end_over_ref"] < 0.05 and not out["weak_flag"])
    return out


def triage_pass(name):
    """Triage FROZEN: E_core(2000) >= 0.2 E_ref ORAZ zero zdarzen."""
    meta = load_meta(name, "h05")
    if meta is None:
        return None
    if meta["status"] == "EVENT":
        return False
    t, psi0, ecore = series(name, "h05")
    i50 = int(round(50.0/DT_OUT))
    i2000 = int(round(2000.0/DT_OUT))
    if i2000 >= len(ecore):
        return None
    E_ref = float(ecore[i50])
    if name == "vac":
        return True   # E_ref=0: 0 >= 0.2*0, zero zdarzen [litera]
    return bool(float(ecore[i2000]) >= 0.2*E_ref)


def cmd_triage():
    print("TRIAGE Etapu A (FROZEN: E_core(2000)>=0.2*E_ref i zero "
          "zdarzen brzegowych):")
    for name in STARTS + ["vac"]:
        tp = triage_pass(name)
        meta = load_meta(name, "h05")
        if meta is None:
            print("  %-14s BRAK BIEGU" % name)
            continue
        t, psi0, ecore = series(name, "h05")
        i50 = int(round(50.0/DT_OUT))
        E_ref = float(ecore[i50]) if len(ecore) > i50 else float("nan")
        frac = (float(ecore[-1])/E_ref if E_ref > 0 else float("nan"))
        print("  %-14s zywy=%s E_ref=%.4e E_core(end)/E_ref=%.4f "
              "status=%s" % (name, tp, E_ref, frac, meta["status"]))


# ---------------------------------------------------------- report --
def _cat_and_conv(name, runs):
    """Klasyfikacja per start wg MD sec. 4 (litera LOCKa sec. 3)."""
    base = runs.get("h05")
    if base is None:
        return "BRAK", "brak biegu bazowego", {}
    det = {}
    if base["status"] == "EVENT":
        # nadkategoria COLLAPSE: zdarzenie kolapsowe na bazie i KAZDYM
        # wykonanym biegu kontrolnym, czasy parami <= 1 j.cz.
        ctrls = [runs[gk] for gk in ("h05_dt2", "h025") if gk in runs]
        if not ctrls:
            return "COLLAPSE?", "brak kontroli dt/2", det
        times = [base["t_event"]] + [c.get("t_event") for c in ctrls]
        evs = [base["event"]] + [c.get("event") for c in ctrls]
        if any(tv is None for tv in times):
            return "INCONCLUSIVE-RUN", \
                "kontrola bez zdarzenia kolapsowego", det
        ok = all(abs(times[0] - tv) <= 1.0 for tv in times[1:])
        sub = "/".join(sorted(set(e for e in evs if e)))
        det["t_events"] = times
        det["subtypes"] = evs
        if ok:
            return "COLLAPSE", "podtyp %s; t=%s" % (
                sub, ",".join("%.3f" % tv for tv in times)), det
        return "INCONCLUSIVE-RUN", \
            "zdarzenia kolapsowe niezbieznie w czasie (>1 j.cz.)", det
    if base["candidate"]:
        # wymagane pelne potwierdzenie h025 + dt/2
        need = [gk for gk in ("h025", "h05_dt2") if gk not in runs]
        if need:
            return "OSCILLON?", "kandydat; brak potwierdzen %s" % need, det
        ok = True
        why = []
        for gk in ("h025", "h05_dt2"):
            c = runs[gk]
            if c["status"] == "EVENT" or not c.get("candidate"):
                ok = False
                why.append("%s: nie-kandydat" % gk)
                continue
            # tau +/-10% z regula cenzurowania
            tb, tc = base["tau"], c["tau"]
            cb, cc = base["censored"], c["censored"]
            if cb and cc:
                pass
            elif cb or cc:
                other = tc if cb else tb
                if other < 0.9*(T_MAX - 50.0):
                    ok = False
                    why.append("%s: tau niezgodne (cenzura)" % gk)
            elif abs(tb - tc)/max(tb, tc) > 0.10:
                ok = False
                why.append("%s: tau %.0f vs %.0f" % (gk, tb, tc))
            ob = base["fft"]["omega"]
            oc = c["fft"]["omega"] if c.get("fft") else None
            if oc is None or abs(ob - oc)/ob > 0.02:
                ok = False
                why.append("%s: omega niezgodna" % gk)
        if ok:
            return "OSCILLON", "potwierdzony h/2 i dt/2", det
        return "INCONCLUSIVE-RUN", "; ".join(why), det
    if base["weak_flag"]:
        need = [gk for gk in ("h025", "h05_dt2") if gk not in runs]
        if need:
            return "OSCILLON-WEAK?", "brak potwierdzen %s" % need, det
        ok = all(runs[gk]["status"] != "EVENT"
                 and runs[gk].get("weak_flag") for gk in
                 ("h025", "h05_dt2"))
        if ok:
            ob = base["fft"]["omega"]
            for gk in ("h025", "h05_dt2"):
                oc = (runs[gk]["fft"] or {}).get("omega")
                if oc is None or abs(ob - oc)/ob > 0.02:
                    ok = False
        if ok:
            return "OSCILLON-WEAK", "potwierdzony h/2 i dt/2", det
        return "INCONCLUSIVE-RUN", "WEAK niepotwierdzony", det
    if base["radiated_flag"]:
        # zbieznosc kategorii/tau sprawdzana tam, gdzie sa kontrole
        if "h025" in runs:
            c = runs["h025"]
            if c["status"] == "EVENT" or not c.get("radiated_flag"):
                return "INCONCLUSIVE-RUN", \
                    "h025 nie potwierdza RADIATED", det
            tb, tc = base["tau"], c["tau"]
            if max(tb, tc) > 0 and abs(tb - tc)/max(tb, tc) > 0.10:
                return "INCONCLUSIVE-RUN", \
                    "tau RADIATED niezbiezne: %.1f vs %.1f" % (tb, tc), det
            return "RADIATED", "zbieznie h,h/2 (tau %.1f/%.1f)" % (
                tb, tc), det
        return "RADIATED", "tau=%.1f (bez kontroli h/2 wg regul " \
            "potwierdzen)" % base["tau"], det
    return "INCONCLUSIVE-RUN", "poza literami kategorii " \
        "(E_end/E_ref=%.3f, cand=False, weak=False)" % \
        (base["E_end_over_ref"] or -1), det


def cmd_report():
    all_an = {}
    for name in STARTS + ["vac"]:
        all_an[name] = {}
        for grid in GRIDS:
            if os.path.exists(rdir(name, grid) + "meta.json"):
                a = analyze_run(name, grid)
                if a:
                    all_an[name][grid] = a
                    with open(rdir(name, grid) + "analysis.json",
                              "w") as f:
                        json.dump(a, f, indent=1)
    lines = []
    lines.append("="*78)
    lines.append("PHASE 3 -- Q-G (litera LOCKa sec. 5; detektor/triage "
                 "FROZEN MD sec. 3-4)")
    lines.append("REJESTR [INPUT]: K_geo=gamma=c0=1; R=400; h {0.05,"
                 "0.025}; dt=0.005 (dt/2=0.0025); t_max=10000; staging "
                 "A t=2000/triage 0.2*E_ref; sponge gamma0=1.0 "
                 "smootherstep [320,400]; E_core r<=80; E_ref="
                 "E_core(50); detektor 0.5*E_ref/100T0/50 przejsc/"
                 "omega_peak<=0.99; potwierdzenia tau 10%, omega 2%; "
                 "RADIATED 0.05; WEAK 1e-3/1000 j.cz.; COLLAPSE pas "
                 "4/3-1e-6, 1e-6, okno 1 j.cz.; brak seeda")
    lines.append("="*78)
    lines.append("Tabela biegow (start x siatka):")
    lines.append("  %-14s %-8s %-9s %-22s %-9s %-6s %-18s %s"
                 % ("start", "siatka", "status", "zdarzenie(t)",
                    "tau", "cross", "omega_peak+/-dom", "E_end/E_ref"))
    for name in STARTS + ["vac"]:
        for grid in ("h05", "h025", "h05_dt2"):
            a = all_an[name].get(grid)
            if a is None:
                continue
            if a["status"] == "EVENT":
                ev = "%s(%.3f)" % (a["event"], a["t_event"])
                lines.append("  %-14s %-8s %-9s %-22s %-9s %-6s %-18s %s"
                             % (name, grid, "EVENT", ev, "-", "-", "-",
                                "-"))
            else:
                fo = a.get("fft")
                om = ("%.4f+/-%.4f" % (fo["omega"], fo["domega"])
                      if fo else "-")
                tau_s = (">=%.0f" % a["tau"] if a["censored"]
                         else "%.1f" % a["tau"])
                ee = ("%.3e" % a["E_end_over_ref"]
                      if a["E_end_over_ref"] is not None else "-")
                lines.append("  %-14s %-8s %-9s %-22s %-9s %-6d %-18s %s"
                             % (name, grid, "OK", "-", tau_s,
                                a["ncross"], om, ee))
    lines.append("")
    lines.append("Triage Etapu A (FROZEN):")
    for name in STARTS:
        tp = triage_pass(name)
        lines.append("  %-14s -> %s" % (name,
                     "ETAP B" if tp else "martwy w Etapie A"))
    lines.append("")
    lines.append("Klasyfikacja per start:")
    cats = {}
    for name in STARTS:
        cat, why, det = _cat_and_conv(name, all_an[name])
        cats[name] = cat
        lines.append("  %-14s : %-16s [%s]" % (name, cat, why))
    a = all_an["vac"].get("h05")
    vac_ok = a is not None and a["status"] != "EVENT" \
        and not a["candidate"]
    lines.append("  kontrola vac : %s (kandydat=False wymagane)"
                 % ("OK" if vac_ok else "ALARM"))
    lines.append("")
    lines.append("Konfrontacja omega_peak vs mapa LP (MIEKKA, "
                 "deskryptywna, BEZ progu -- MD sec. 6):")
    lp = {0.02: 0.997683, 0.05: 0.985521, 0.08: 0.962933,
          0.10: 0.942083}
    lines.append("  (* = omega_desc: deskryptywnie, pik FFT na calym "
                 "oknie podtrzymania < wymogu 1000 j.cz. detektora --")
    lines.append("   correction note 1b; NIE wchodzi do detektora)")
    for name in STARTS:
        aa = float(name.split("_")[1][1:])
        b = all_an[name].get("h05")
        fo = b.get("fft") if b and b["status"] != "EVENT" else None
        star = ""
        if fo is None and b and b["status"] != "EVENT":
            fo = b.get("fft_desc")
            star = "*"
        om = ("%.4f+/-%.4f%s" % (fo["omega"], fo["domega"], star)
              if fo else "-")
        lines.append("  %-14s omega_peak=%-18s omega_LP(a=%.2f)=%.4f"
                     % (name, om, aa, lp[aa]))
    # ---- werdykt Q-G wg litery LOCKa sec. 5 ----
    n_osc = sum(1 for c in cats.values() if c == "OSCILLON")
    n_weak = sum(1 for c in cats.values()
                 if c.startswith("OSCILLON-WEAK"))
    all_rad = all(c == "RADIATED" for c in cats.values())
    if n_osc >= 1:
        verdict = "Q-G-PASS"
    elif all_rad and n_weak == 0:
        verdict = "Q-G-FAIL"
    else:
        verdict = "Q-G-INCONCLUSIVE"
    lines.append("")
    kl = {}
    for c in cats.values():
        kl[c] = kl.get(c, 0) + 1
    lines.append("KLASY: " + " ".join("%s=%d" % kv
                                      for kv in sorted(kl.items())))
    lines.append("WERDYKT: %s (PASS: >=1 OSCILLON potwierdzony h/2 i "
                 "dt/2; FAIL: WSZYSTKIE 12 RADIATED zbieznie ORAZ zero "
                 "OSCILLON-WEAK; INCONCLUSIVE: reszta)" % verdict)
    lines.append("="*78)
    with open(BASE + "Phase3_output.txt", "w") as f:
        f.write("\n".join(lines) + "\n")
    with open(RES + "verdict.json", "w") as f:
        json.dump({"verdict": verdict, "categories": cats,
                   "n_oscillon": n_osc, "n_weak": n_weak}, f, indent=1)
    print("\n".join(lines))


if __name__ == "__main__":
    os.makedirs(RES, exist_ok=True)
    if sys.argv[1] == "run":
        do_run(sys.argv[2], sys.argv[3], float(sys.argv[4]))
    elif sys.argv[1] == "triage":
        cmd_triage()
    elif sys.argv[1] == "report":
        cmd_report()
    else:
        raise SystemExit("unknown cmd")
