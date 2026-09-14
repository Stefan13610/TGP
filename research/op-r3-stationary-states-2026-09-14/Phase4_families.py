#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
op-r3-stationary-states (Phase 4) -- Q-F: struktura dyskretna stanow
dlugozyciowych (TYLKO przy Q-E-PASS). LOCK sec. 2 Phase 4;
MD sec. 9 (FROZEN): obwiednia A(r)=<|psi-1|>_t po oknie stabilnym
(bloki 50 j.cz. wewnatrz okna FFT), r<=60, wygladzenie srednia
ruchoma 0.45 j.dl. (9 pkt h=0.05 / 17 pkt h=0.025); wezel = minimum
lokalne A z prominencja A(min) < 0.3*min(sasiednie maksima).
Werdykt Q-F wg litery LOCKa:
  PASS: >=2 rodziny wezlowe ORAZ monotonia omega(n) ORAZ brak rodzin
        n>=3 o czasie zycia >= progu (100 T0);
  PARTIAL: tylko n=0; FAIL: brak struktury dyskretnej.
"""
import json
import sys
import numpy as np

BASE = ("C:/Users/Mateusz/Documents/ObsydnianMain/TGP/TGP_v1/research/"
        "op-r3-stationary-states-2026-09-14/")
RESDIR = BASE + "Phase3_results/"
T0 = 2.0*np.pi
out = []


def em(s=""):
    print(s, flush=True)
    out.append(s)


def load_json(rid):
    with open(RESDIR + rid + ".json") as f:
        return json.load(f)


def envelope_nodes(rid, h):
    """A(r) z blokow Ablocks w oknie stabilnym; liczba wezlow wg MD
    sec. 9. Zwraca (n_nodes, r_nodes, A, r)."""
    d = np.load(RESDIR + rid + ".npz")
    res = load_json(rid)
    Ab = d["Ablocks"]
    t_hold = res["t_hold"]
    w0 = max(300.0, t_hold - 600.0)
    blocks = [j for j in range(Ab.shape[0])
              if 50.0*j >= w0 - 1e-9 and 50.0*(j + 1) <= t_hold + 1e-9]
    if not blocks:
        blocks = [Ab.shape[0] - 1]
    A = np.mean(Ab[blocks, :], axis=0)
    N = len(A)
    r = (np.arange(N) + 0.5)*h
    npts = 9 if h == 0.05 else 17
    ker = np.ones(npts)/npts
    As = np.convolve(A, ker, mode="same")
    m = r <= 60.0
    As = As[m]
    rr = r[m]
    # minima lokalne wewnetrzne
    nodes = []
    for i in range(1, len(As) - 1):
        if As[i] < As[i-1] and As[i] <= As[i+1]:
            # sasiednie maksima
            jl = i
            while jl > 0 and As[jl-1] >= As[jl]:
                jl -= 1
            while jl > 0 and not (As[jl] > As[jl-1]
                                  and As[jl] >= As[jl+1]):
                jl -= 1
            jr = i
            while jr < len(As) - 1 and As[jr+1] >= As[jr]:
                jr += 1
            while jr < len(As) - 1 and not (As[jr] > As[jr+1]
                                            and As[jr] >= As[jr-1]):
                jr += 1
            amaxl, amaxr = As[jl], As[jr]
            if As[i] < 0.3*min(amaxl, amaxr):
                nodes.append(float(rr[i]))
    # scal wezly blizsze niz 0.9 (artefakty plaskich minimow)
    merged = []
    for rn in nodes:
        if merged and rn - merged[-1] < 0.9:
            continue
        merged.append(rn)
    return len(merged), merged, As, rr


def main():
    with open(RESDIR + "verdict.json") as f:
        v = json.load(f)
    if v["qe"] != "Q-E-PASS":
        em("Q-E != PASS (%s) -- Phase 4 nie uruchamiana wg LOCKa."
           % v["qe"])
        with open(BASE + "Phase4_output.txt", "w",
                  encoding="ascii") as fh:
            fh.write("\n".join(out) + "\n")
        return
    osc = [s for s, c in v["classes"].items() if c == "OSCILLON"]
    em("=" * 78)
    em("PHASE 4 -- Q-F: rodziny (n, omega, tau); oscylony: %s" % osc)
    em("MD sec. 9 FROZEN: obwiednia r<=60, wygladzenie 0.45,")
    em("  prominencja 0.3; siatka PRIMARY strukturalnie: h=0.025;")
    em("  h=0.05 kontrolnie")
    em("=" * 78)
    fam = {}
    rows = []
    for s in osc:
        r025 = load_json("%s__h025" % s)
        r05 = load_json("%s__h05" % s)
        n025, rn025, _, _ = envelope_nodes("%s__h025" % s, 0.025)
        n05, rn05, _, _ = envelope_nodes("%s__h05" % s, 0.05)
        om = r025["omega"]
        tau = r025["tau_str"]
        em("  %-16s: n(h025)=%d %s | n(h05)=%d %s | omega=%.5f"
           " (h05: %.5f) | tau=%s | harmoniki(h025): %s"
           % (s, n025, ["%.1f" % x for x in rn025],
              n05, ["%.1f" % x for x in rn05],
              om, r05["omega"], tau,
              [("%.3f" % w, "%.3f" % p)
               for w, p in r025["harmonics"][:3]]))
        rows.append((s, n025, n05, om, r025["tau"],
                     r025["censored"]))
        fam.setdefault(n025, []).append((s, om, r025["tau"],
                                         r025["censored"]))
    em("")
    em("Rodziny wg n (h=0.025):")
    for n in sorted(fam):
        oms = [o for _, o, _, _ in fam[n]]
        em("  n=%d: %d stan(y); omega: min=%.5f max=%.5f srednia=%.5f"
           % (n, len(fam[n]), min(oms), max(oms),
              sum(oms)/len(oms)))
    em("")
    em("Klasa czestosci (litera Q-F(b)): %s"
       % ("; ".join("%s: omega=%.5f %s prog kontinuum m=1"
                    % (s, o, "<" if o < 1.0 else ">")
                    for s, _, _, o, _, _ in [(a, b, c, d, e, f)
                                             for a, b, c, d, e, f
                                             in rows])))
    ns = sorted(fam)
    grid_consistent = all(n25 == n05 for _, n25, n05, _, _, _ in rows)
    em("Zgodnosc n miedzy siatkami: %s"
       % ("TAK" if grid_consistent else
          "NIE -- raport per przypadek wyzej"))
    # monotonia omega(n) po srednich rodzin
    means = [sum(o for _, o, _, _ in fam[n])/len(fam[n]) for n in ns]
    mono = all(means[i] < means[i+1] for i in range(len(means) - 1)) \
        or all(means[i] > means[i+1] for i in range(len(means) - 1))
    n3_long = any(n >= 3 and any((c or t >= 100.0*T0)
                                 for _, _, t, c in fam[n])
                  for n in ns)
    if len(ns) >= 2 and mono and not n3_long:
        qf = "Q-F-PASS"
    elif ns == [0]:
        qf = "Q-F-PARTIAL"
    else:
        qf = "Q-F-FAIL"
    em("")
    em("Kryteria litery: rodziny=%d (>=2: %s); monotonia omega(n): %s;"
       % (len(ns), len(ns) >= 2, mono))
    em("  rodziny n>=3 dlugozyciowe: %s (wymagane: brak)" % n3_long)
    em("WERDYKT: %s" % qf)
    em("=" * 78)
    with open(RESDIR + "verdict_qf.json", "w") as f:
        json.dump(dict(qf=qf, families={str(k): [(a, b, c, d)
                                                 for a, b, c, d
                                                 in fam[k]]
                                        for k in fam}), f, indent=1)
    with open(BASE + "Phase4_output.txt", "w", encoding="ascii") as fh:
        fh.write("\n".join(out) + "\n")
    print("zapisano:", BASE + "Phase4_output.txt")


if __name__ == "__main__":
    main()
