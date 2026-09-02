#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
op-metric-closure-relaxation (Phase 1) -- bramka maszynerii (LOCK sec. 3
Phase 1; specyfikacja FROZEN: Phase_method_decisions.md sec. 6).

P1a: proznia g=1 bez zaburzen zostaje w g=1 (PRIMARY i C-BAR),
     gate max_t ||g-1||_inf <= 1e-10, t=10; geometrie: radial h=0.0125,
     3D L=2pi N=32, 3D L=4pi N=48.
P1b: reprodukcja BREAKDOWN g->+inf BEZ domkniec (w=1, bez kar), start
     solitonowy, h in {0.025,0.0125}; PASS <=> BREAKDOWN z t w [2.7,3.2]
     na obu siatkach (poprzednik: 2.75/3.13) -- osiagalny FAIL ciaglosci.
P1c: sektor stabilny m^2=+gamma: relaksacja do prozni, ZERO alarmow OBU
     detektorow; macierz po korekcie (Phase_correction_note_p1c_matrix.md
     -- zasada domenowa MD sec. 6; tlo 2pi ma g_max=1.4734>g_ceil, wiec
     lat x PRIMARY poza domena kontroli, jak sol x PRIMARY):
     lat N=32 x CBAR, gen N=48 x {PRIM,CBAR}, sol h=0.025 x CBAR.
     Pierwotny bieg (z lat x PRIM): Phase1_output_pre_correction.txt.
Dowolny FAIL => STOP.
"""
import sys
import numpy as np

sys.path.insert(0, "C:/Users/Mateusz/Documents/ObsydnianMain/TGP/TGP_v1/"
                   "research/op-metric-closure-relaxation-2026-09-02")
import Phase2_two_sided_relax as eng


def p1a():
    print("=" * 78)
    print("P1a -- stacjonarnosc prozni (oba warianty; gate <= 1e-10, t=10)")
    print("=" * 78)
    eng.registry_banner(eng.FLOORS[eng.FI_PRIMARY], "PRIM+CBAR", "tach")
    allpass = True
    for vname in ("PRIM", "CBAR"):
        var = eng.Variant(vname, "tach", eng.FI_PRIMARY)
        geoms = []
        fr = eng.FlowRadial(0.0125, var)
        geoms.append(("radial h=0.0125", fr, np.ones(len(fr.r))))
        f32 = eng.Flow3D(32, eng.L_LAT, var)
        geoms.append(("3D L=2pi N=32", f32, np.ones((32, 32, 32))))
        f48 = eng.Flow3D(48, eng.L_GEN, var)
        geoms.append(("3D L=4pi N=48", f48, np.ones((48, 48, 48))))
        for gname, flow, g in geoms:
            dev = 0.0
            for k in range(int(round(10.0 / eng.DT_MAIN))):
                g, _ = flow.step(g, eng.DT_MAIN)
                dev = max(dev, float(np.max(np.abs(g - 1.0))))
            ok = dev <= 1e-10
            allpass = allpass and ok
            print("  P1a %s %-16s: max||g-1||_inf = %.3e -> %s"
                  % (vname, gname, dev, "PASS" if ok else "FAIL"),
                  flush=True)
    print("P1a: %s" % ("PASS" if allpass else "FAIL -> STOP"))
    return allpass


def p1b():
    print("=" * 78)
    print("P1b -- reprodukcja BREAKDOWN bez domkniec (w=1, bez kar); "
          "oczekiwane t w [2.7,3.2]")
    print("=" * 78)
    var = eng.Variant("NONE", "tach", eng.FI_PRIMARY)
    allpass = True
    for h in eng.H_RAD:
        flow = eng.FlowRadial(h, var)
        g0 = eng.start_radial_sol(h)
        res, _ = eng.run_flow(flow, g0, eng.FLOORS[eng.FI_PRIMARY],
                              eng.DT_MAIN, "p1b_sol_h%g" % h, tmax=20.0)
        ok = res["status"] == "BREAKDOWN" and 2.7 <= res["t_end"] <= 3.2
        allpass = allpass and ok
        print("  P1b h=%g: %s t=%.3f (wymog BREAKDOWN, t in [2.7,3.2]) "
              "-> %s" % (h, res["status"], res["t_end"],
                         "PASS" if ok else "FAIL"), flush=True)
    print("P1b: %s" % ("PASS" if allpass else "FAIL -> STOP"))
    return allpass


def p1c():
    print("=" * 78)
    print("P1c -- sektor stabilny m^2=+gamma: relaksacja do prozni, ZERO "
          "alarmow OBU detektorow (macierz FROZEN MD sec. 6)")
    print("=" * 78)
    # macierz po korekcie (Phase_correction_note_p1c_matrix.md):
    # lat x PRIMARY wylaczony z KONTROLI (start czesciowo poza domena
    # metryki PRIMARY: g_max=1.4734 > g_ceil; zasada domenowa MD sec. 6)
    runs = [("lat", "CBAR", dict(N=32)),
            ("gen", "PRIM", dict(N=48)), ("gen", "CBAR", dict(N=48)),
            ("sol", "CBAR", dict(h=0.025))]
    allpass = True
    for start, vname, kw in runs:
        var = eng.Variant(vname, "stab", eng.FI_PRIMARY)
        if start == "sol":
            flow = eng.FlowRadial(kw["h"], var)
            g0 = eng.start_radial_sol(kw["h"])
        elif start == "lat":
            flow = eng.Flow3D(kw["N"], eng.L_LAT, var)
            g0 = eng.start_lat(kw["N"])
        else:
            flow = eng.Flow3D(kw["N"], eng.L_GEN, var)
            g0 = eng.start_gen(kw["N"])
        lbl = "p1c_%s_%s" % (start, vname)
        res, gfin = eng.run_flow(flow, g0, eng.FLOORS[eng.FI_PRIMARY],
                                 eng.DT_MAIN, lbl)
        devvac = float(np.max(np.abs(gfin - 1.0)))
        ok = (res["status"] == "STATIONARY" and devvac <= 1e-3
              and res["nuc_dir"] is None)
        allpass = allpass and ok
        print("  P1c %s %s: %s t=%.1f, ||g-1||_inf=%.2e, alarmy: %s "
              "-> %s" % (start, vname, res["status"], res["t_end"],
                         devvac, res["nuc_dir"] or "ZERO",
                         "PASS" if ok else "FAIL"), flush=True)
    print("P1c: %s" % ("PASS" if allpass else "FAIL -> STOP"))
    return allpass


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    ok = True
    if which in ("all", "p1a"):
        ok = p1a() and ok
        if not ok:
            print("STOP (P1a FAIL)")
            sys.exit(1)
    if which in ("all", "p1b"):
        ok = p1b() and ok
        if not ok:
            print("STOP (P1b FAIL)")
            sys.exit(1)
    if which in ("all", "p1c"):
        ok = p1c() and ok
        if not ok:
            print("STOP (P1c FAIL)")
            sys.exit(1)
    print("\nPHASE 1 GATE: %s" % ("PASS -> Phase 2" if ok else "FAIL"))
    sys.exit(0 if ok else 1)
