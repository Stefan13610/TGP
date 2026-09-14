#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
op-action-audit-spectrum-insert (Phase 3) -- Q-D2: DeltaE_insert(A;R,h)
na WSPOLNYM tle prozni psi=1 (korekta metodologiczna P0.3 audytu).

LOCK: Phase0_balance.md sec. 1-2 Phase 3; MD (FROZEN) sec. 4:
- E[psi] = 4pi int_0^R r^2 [ 1/2 Keff (psi')^2 + Ueff ] dr (trapez),
  Keff = psi^4, Ueff = gamma(psi^4/4 - psi^3/3) [PRIMARY, odczyt B];
- siatka WEZLOWA r_i = i h (wezel r=0 istnieje); Neumann w 0 i R przez
  strukture strumieniowa (zerowy strumien przez brzegi);
- gradient flow psi_t = -dE/dpsi (masa wezla m_i = c_i h r_i^2,
  m_0 = h^3/24 -- dokladna calka r^2 na [0,h/2]; parametryzacja flow,
  nie zmienia E ani stanow stacjonarnych);
- semi-implicit Euler (I - dt A_t L) dpsi = dt rhs, A_t=1.05 max Keff,
  dt=0.01, t_max=200, stacjonarnosc ||psidot||_inf <= 1e-8 co dt=1;
- WIEZ pinu (jedyny; LOCK): po KAZDYM kroku projekcja psi[0] <- A;
  residuum pinu |psi0_pre-projekcja - A| raportowane;
- ZERO podlog/barier; pas psi > 4/3-1e-6 => BREAKDOWN-BOUNDARY;
  min psi < 1e-6 => BREAKDOWN-BOUNDARY-LOWER; niefinitycznosc =>
  BREAKDOWN (klasyfikacje deskryptywne, stop);
- DeltaE_insert = E[psi_A^relax] - E[psi==1] na IDENTYCZNEJ siatce;
  implementacja dobrze uwarunkowana: sumowanie POINTWISE roznicy
  Ueff(psi)-Ueff(1) (algebraicznie ta sama wielkosc; kontrola krzyzowa
  surowa roznica E raportowana);
- P3a: proznia bez wiezu t=10 dryf <= 1e-10 (4 pudla); pin A=1 =>
  DeltaE = 0 +- 1e-10 (4 pudla); FAIL => STOP;
- P3b: A in {0.50, 0.70, 5/6, 7/6, 1.25, 1.30} x R in {60,120} x
  h in {0.025, 0.0125} = 24 relaksacje; start psi0 = 1 +
  (A-1) exp(-r^2/(2 sigma^2)), sigma = 5.0 [INPUT-MD];
- zbieznosc h: |DeltaE(hc)-DeltaE(hf)| / max(|DeltaE(hf)|,1e-6)
  <= 5e-3; R: |DeltaE(120)-DeltaE(60)| raportowane.

Uzycie:
  python Phase3_insert_cost.py gates
  python Phase3_insert_cost.py matrix
  python Phase3_insert_cost.py verdict
  python Phase3_insert_cost.py all

REJESTR WEJSC [INPUT, flagowane]: K_geo=gamma=c0=1; rodzina A
{0.50,0.70,5/6,7/6,1.25,1.30}; R {60,120}; h {0.025,0.0125}; dt=0.01;
t_max=200; prog stacjonarnosci 1e-8; sigma=5.0 [INPUT-MD]; pas
4/3-1e-6; progi P3a 1e-10; zbieznosc h 5e-3; brak seeda
(deterministycznie). Formy CYTAT sek08a/sek08c (MD sec. 1); odczyt B
DZIEDZICZONY (MD sec. 2).
"""
import json
import os
import sys
import time
import numpy as np
from scipy.linalg import solve_banded

GAMMA = 1.0                                # INPUT (LOCK sec. 1)
KGEO = 1.0                                 # INPUT (LOCK sec. 1)
A_FAMILY = [0.50, 0.70, 5.0 / 6.0, 7.0 / 6.0, 1.25, 1.30]  # LOCK
A_LABELS = ["A050", "A070", "A56", "A76", "A125", "A130"]
R_LIST = [60.0, 120.0]                     # LOCK
H_LIST = [0.025, 0.0125]                   # LOCK
DT = 0.01                                  # LOCK
T_MAX = 200.0                              # LOCK
STAT_TOL = 1e-8                            # LOCK
SIGMA = 5.0                                # INPUT-MD (dziedziczone)
PSI_LIMIT = 4.0 / 3.0
PSI_BAND = PSI_LIMIT - 1e-6                # LOCK (pas graniczny)
PSI_LOW = 1e-6                             # MD sec. 4 (deskryptywny)
GATE_TOL = 1e-10                           # LOCK (P3a)
CONV_H = 5e-3                              # LOCK (P3b)
BASE = ("C:/Users/Mateusz/Documents/ObsydnianMain/TGP/TGP_v1/research/"
        "op-action-audit-spectrum-insert-2026-09-13/")
RESDIR = BASE + "Phase3_results/"

t0_wall = time.time()


def stamp(msg):
    print("[t=%7.1fs] %s" % (time.time() - t0_wall, msg), flush=True)


def registry_banner(extra=""):
    print("REJESTR [INPUT]: K_geo=gamma=c0=1; A=%s; R=%s; h=%s; dt=%g; "
          "t_max=%g; stat=1e-8; sigma=%.1f [INPUT-MD]; pas=4/3-1e-6; "
          "P3a=1e-10; conv_h=5e-3; formy: w=psi/(4-3psi) "
          "[eq:vol-element-M911], V=-g*psi^2(4-3psi)^2/12 [eq:V-M911], "
          "Keff=psi^4 [odczyt B DZIEDZICZONY], Ueff=w*V=psi^4/4-psi^3/3 "
          "[tozsamosc P1b]%s"
          % (A_FAMILY, R_LIST, H_LIST, DT, T_MAX, SIGMA, extra),
          flush=True)


# ------------------------------------------------- model (FROZEN, MD sec. 4)
def Ueff(p):
    return GAMMA * (p ** 4 / 4.0 - p ** 3 / 3.0)


def Ueffp(p):
    return GAMMA * p * p * (p - 1.0)


def Keff(p):
    return KGEO * p ** 4


def Keffp(p):
    return 4.0 * KGEO * p ** 3


UVAC = Ueff(1.0)                            # = -1/12 (proznia)


# --------------------------------------------- silnik radialny (wezlowy)
class FlowRadialPin:
    """Siatka wezlowa r_i = i h, i=0..N; trapez + 4pi; Neumann przez
    zerowy strumien na brzegach; masa m_0 = h^3/24 (MD sec. 4)."""

    def __init__(self, R, h):
        N = int(round(R / h))
        assert abs(N * h - R) < 1e-12
        self.N = N
        self.h = h
        self.R = R
        self.r = np.arange(N + 1) * h                  # wezly, r_0 = 0
        self.rm2 = ((self.r[:-1] + self.r[1:]) * 0.5) ** 2   # midpointy^2
        c = np.ones(N + 1)
        c[0] = 0.5
        c[N] = 0.5
        self.ctrap = c                                  # wagi trapezu
        self.m = c * h * self.r ** 2                    # masy wezlow
        self.m[0] = h ** 3 / 24.0                       # MD sec. 4
        self.a = self.rm2 / h                           # przewodnosci L

    def energy_rel(self, g):
        """DeltaE wzgledem prozni psi==1 (pointwise, dobrze uwarunkowane):
        4pi [ sum_h r_mid^2 1/2 Keff(mid) (dpsi/h)^2
              + sum trapez r^2 (Ueff(psi)-Ueff(1)) ]."""
        h = self.h
        gm = 0.5 * (g[:-1] + g[1:])
        dg = np.diff(g) / h
        grad = 0.5 * np.sum(self.rm2 * Keff(gm) * dg ** 2) * h
        pot = np.sum(self.ctrap * self.r ** 2 * (Ueff(g) - UVAC)) * h
        return 4.0 * np.pi * (grad + pot)

    def energy_raw(self, g):
        h = self.h
        gm = 0.5 * (g[:-1] + g[1:])
        dg = np.diff(g) / h
        grad = 0.5 * np.sum(self.rm2 * Keff(gm) * dg ** 2) * h
        pot = np.sum(self.ctrap * self.r ** 2 * Ueff(g)) * h
        return 4.0 * np.pi * (grad + pot)

    def rhs(self, g):
        """psidot = -dE_h/dpsi_i / (4pi m_i); struktura strumieniowa
        dziedziczona (t_flux, t_quad); 4pi kasuje sie w parametryzacji
        flow (wchodzi wspolnie do dE i normy L2 z miara 4pi r^2 dr)."""
        h = self.h
        gm = 0.5 * (g[:-1] + g[1:])
        dg = np.diff(g) / h
        dH = self.ctrap * h * self.r ** 2 * Ueffp(g)
        t_flux = self.rm2 * Keff(gm) * dg
        t_quad = 0.25 * h * self.rm2 * Keffp(gm) * dg ** 2
        dH[:-1] += -t_flux + t_quad
        dH[1:] += t_flux + t_quad
        return -dH / self.m

    def step(self, g, dt):
        r = self.rhs(g)
        A = 1.05 * float(np.max(Keff(g)))
        N1 = self.N + 1
        lo = np.zeros(N1)
        di = np.ones(N1)
        up = np.zeros(N1)
        cda = dt * A * self.a
        di[:-1] += cda / self.m[:-1]
        di[1:] += cda / self.m[1:]
        up[1:] = -cda / self.m[:-1]
        lo[:-1] = -cda / self.m[1:]
        ab = np.vstack([up, di, lo])
        dgv = solve_banded((1, 1), ab, dt * r)
        return g + dgv, r


def start_profile(flow, A):
    """Start deterministyczny LOCKa: gauss sigma=5.0 wokol 0,
    psi(inf)=1, amplituda dopasowana do A; psi0(0) = A dokladnie."""
    return 1.0 + (A - 1.0) * np.exp(-flow.r ** 2 / (2.0 * SIGMA ** 2))


def run_relax(flow, psi0, pinA=None, tmax=T_MAX, label=""):
    """Gradient flow do stacjonarnosci / pasa / zalamania / t_max.
    pinA=None: bez wiezu. pinA=A: projekcja psi[0] <- A po KAZDYM
    kroku (MD sec. 4); residuum pinu raportowane."""
    np.seterr(over='ignore', invalid='ignore')
    g = psi0.copy()
    if pinA is not None:
        g[0] = pinA
    steps_per_unit = int(round(1.0 / DT))
    nsteps = int(round(tmax / DT))
    status, t_end = "TMAX", tmax
    pin_res_max = 0.0
    pin_res_last = 0.0
    pnorm = float(np.max(np.abs(flow.rhs(g)[1:] if pinA is not None
                                else flow.rhs(g))))
    drift_max = 0.0
    E_series = [(0.0, flow.energy_rel(g))]
    for k in range(1, nsteps + 1):
        try:
            g_new, _ = flow.step(g, DT)
            finite = bool(np.all(np.isfinite(g_new)))
        except (ValueError, FloatingPointError, OverflowError):
            finite = False
            g_new = g
        if not finite:
            status, t_end = "BREAKDOWN", k * DT
            break
        if pinA is not None:
            res = abs(float(g_new[0]) - pinA)
            pin_res_max = max(pin_res_max, res)
            pin_res_last = res
            g_new[0] = pinA                      # PROJEKCJA PINU (MD s.4)
        g = g_new
        mx = float(np.max(g))
        mn = float(np.min(g))
        if mx >= PSI_BAND:
            status, t_end = "BREAKDOWN-BOUNDARY", k * DT
            break
        if mn <= PSI_LOW:
            status, t_end = "BREAKDOWN-BOUNDARY-LOWER", k * DT
            break
        drift_max = max(drift_max, float(np.max(np.abs(g - 1.0))))
        if k % steps_per_unit == 0:
            t = k * DT
            r = flow.rhs(g)
            if pinA is not None:
                r = r[1:]                        # wezel pinowany: psidot=0
            pnorm = float(np.max(np.abs(r)))
            E_series.append((t, flow.energy_rel(g)))
            if pnorm <= STAT_TOL:
                status, t_end = "STATIONARY", t
                break
    dE = flow.energy_rel(g)
    res = dict(label=label, status=status, t_end=t_end, pnorm_end=pnorm,
               dE_rel=dE, E_raw=flow.energy_raw(g),
               E_raw_vac=flow.energy_raw(np.ones_like(g)),
               pin_res_max=pin_res_max, pin_res_last=pin_res_last,
               pmin=float(np.min(g)), pmax=float(np.max(g)),
               pmax_vs_43=float(np.max(g)) - PSI_LIMIT,
               drift_max=drift_max, E_series=E_series[-5:])
    return res, g


def profile_stats(flow, g, A):
    """Rozciaglosc: najwieksze r z |psi-1| >= 0.5|A-1|; + psi w probkach."""
    dev = np.abs(g - 1.0)
    thr = 0.5 * abs(A - 1.0)
    idx = np.where(dev >= thr)[0]
    r_half = float(flow.r[idx[-1]]) if len(idx) else 0.0
    n_half = int(len(idx))
    samples = {}
    for rr in (0.0, flow.h, 2 * flow.h, 0.5, 1.0, 2.0, 5.0, 10.0):
        i = int(round(rr / flow.h))
        if i <= flow.N:
            samples["r=%g" % rr] = float(g[i])
    return dict(r_half=r_half, n_nodes_half=n_half, samples=samples)


# ----------------------------------------------------------------- P3a
def gates():
    os.makedirs(RESDIR, exist_ok=True)
    registry_banner("; tryb=gates(P3a)")
    out = {"vac": {}, "pin1": {}}
    ok_all = True
    for R in R_LIST:
        for h in H_LIST:
            key = "R%g_h%g" % (R, h)
            flow = FlowRadialPin(R, h)
            # (i) proznia bez wiezu, t=10, dryf <= 1e-10
            res, g = run_relax(flow, np.ones(flow.N + 1), pinA=None,
                               tmax=10.0, label="vac_" + key)
            drift = res["drift_max"]
            ok = drift <= GATE_TOL and not res["status"].startswith(
                "BREAKDOWN")
            ok_all = ok_all and ok
            out["vac"][key] = dict(drift=drift, status=res["status"],
                                   ok=bool(ok))
            stamp("P3a(i) %s: dryf=%.3e status=%s => %s"
                  % (key, drift, res["status"], "PASS" if ok else "FAIL"))
            # (ii) pin A=1: DeltaE = 0 +- 1e-10
            res2, g2 = run_relax(flow, start_profile(flow, 1.0), pinA=1.0,
                                 tmax=T_MAX, label="pin1_" + key)
            dE = res2["dE_rel"]
            ok2 = abs(dE) <= GATE_TOL and not res2["status"].startswith(
                "BREAKDOWN")
            ok_all = ok_all and ok2
            out["pin1"][key] = dict(dE=dE, status=res2["status"],
                                    t_end=res2["t_end"],
                                    pin_res_max=res2["pin_res_max"],
                                    ok=bool(ok2))
            stamp("P3a(ii) %s: DeltaE=%.3e status=%s t=%.1f pin_res=%.2e"
                  " => %s" % (key, dE, res2["status"], res2["t_end"],
                              res2["pin_res_max"],
                              "PASS" if ok2 else "FAIL"))
    out["P3a"] = "PASS" if ok_all else "FAIL"
    with open(RESDIR + "gates_P3a.json", "w") as f:
        json.dump(out, f, indent=1)
    stamp("P3a: %s (zapisano %sgates_P3a.json)"
          % (out["P3a"], RESDIR))
    if not ok_all:
        raise SystemExit("P3a FAIL => STOP (litera LOCKa)")


# ----------------------------------------------------------------- P3b
def matrix():
    os.makedirs(RESDIR, exist_ok=True)
    registry_banner("; tryb=matrix(P3b)")
    for A, lab in zip(A_FAMILY, A_LABELS):
        for R in R_LIST:
            for h in H_LIST:
                jid = "%s_R%g_h%g" % (lab, R, h)
                fp = RESDIR + jid + ".json"
                if os.path.exists(fp):
                    stamp("%s: JUZ policzony (pomijam)" % jid)
                    continue
                flow = FlowRadialPin(R, h)
                res, g = run_relax(flow, start_profile(flow, A), pinA=A,
                                   tmax=T_MAX, label=jid)
                res.update(A=A, R=R, h=h,
                           profile=profile_stats(flow, g, A))
                with open(fp, "w") as f:
                    json.dump(res, f, indent=1)
                np.savez_compressed(RESDIR + jid + ".npz", psi=g,
                                    r=flow.r)
                stamp("%s: %s t=%.1f DeltaE=%+.8e psi[%.6f,%.6f] "
                      "pin_res_max=%.2e r_half=%.3f"
                      % (jid, res["status"], res["t_end"], res["dE_rel"],
                         res["pmin"], res["pmax"], res["pin_res_max"],
                         res["profile"]["r_half"]))


# --------------------------------------------------------------- werdykt
def verdict():
    out = []

    def em(s=""):
        print(s, flush=True)
        out.append(s)

    em("=" * 78)
    em("WERDYKT Q-D2 (litera LOCKa sec. 2 Phase 3)")
    em("REJESTR [INPUT]: K_geo=gamma=c0=1; A=%s;" % A_FAMILY)
    em("  R=%s h=%s dt=%g t_max=%g stat=1e-8 sigma=%.1f [INPUT-MD];"
       % (R_LIST, H_LIST, DT, T_MAX, SIGMA))
    em("  pas=4/3-1e-6; P3a tol=1e-10; conv_h<=5e-3; formy CYTAT MD s.1;")
    em("  odczyt B DZIEDZICZONY (MD s.2); wiez: pin psi(0)=A (projekcja)")
    em("=" * 78)
    with open(RESDIR + "gates_P3a.json") as f:
        g3a = json.load(f)
    em("P3a: %s (dryf prozni i pin A=1 na 4 pudlach; szczegoly json)"
       % g3a["P3a"])
    em()
    em("P3b: tabela DeltaE_insert(A;R,h) [E(psi_A^relax) - E(psi==1),")
    em("     wspolne tlo/pudlo/siatka/brzeg]:")
    hdr = ("%-6s %-6s | %-14s %-14s | %-9s %-6s | %-14s | %s"
           % ("A", "R", "dE(h=0.025)", "dE(h=0.0125)", "conv_h",
              "<=5e-3", "|dE_R120-R60|", "status/prof"))
    em(hdr)
    em("-" * len(hdr))
    per_A = {}
    all_ok = True
    any_boundary = False
    for A, lab in zip(A_FAMILY, A_LABELS):
        per_A[lab] = {}
        row_R = {}
        for R in R_LIST:
            rr = {}
            for h in H_LIST:
                jid = "%s_R%g_h%g" % (lab, R, h)
                with open(RESDIR + jid + ".json") as f:
                    rr[h] = json.load(f)
            dEc, dEf = rr[H_LIST[0]]["dE_rel"], rr[H_LIST[1]]["dE_rel"]
            conv = abs(dEc - dEf) / max(abs(dEf), 1e-6)
            statuses = [rr[h]["status"] for h in H_LIST]
            stat_ok = all(s == "STATIONARY" for s in statuses)
            if any(s.startswith("BREAKDOWN-BOUNDARY") for s in statuses):
                any_boundary = True
            row_R[R] = dict(dEc=dEc, dEf=dEf, conv=conv,
                            statuses=statuses, stat_ok=stat_ok,
                            prof=rr[H_LIST[1]].get("profile", {}),
                            pmin=rr[H_LIST[1]]["pmin"],
                            pmax=rr[H_LIST[1]]["pmax"])
        dR = abs(row_R[120.0]["dEf"] - row_R[60.0]["dEf"])
        for R in R_LIST:
            d = row_R[R]
            okh = d["conv"] <= CONV_H and d["stat_ok"]
            all_ok = all_ok and okh
            pr = d["prof"]
            em("%-6.4f %-6g | %+.7e %+.7e | %.3e %-6s | %-14s | %s "
               "r_half=%.2f psi[%.4f,%.4f]"
               % (A, R, d["dEc"], d["dEf"], d["conv"],
                  "TAK" if d["conv"] <= CONV_H else "NIE",
                  ("%.3e" % dR) if R == 120.0 else "",
                  ",".join(s[:10] for s in d["statuses"]),
                  pr.get("r_half", -1), d["pmin"], d["pmax"]))
        per_A[lab] = dict(A=A, row_R={str(k): v for k, v in row_R.items()},
                          dR=dR)
    em()
    em("Werdykt per A (litera: zbieznosc h + znak):")
    signs_conv = []
    for A, lab in zip(A_FAMILY, A_LABELS):
        pa = per_A[lab]
        convs = [pa["row_R"][str(R)]["conv"] for R in R_LIST]
        stats_ok = all(pa["row_R"][str(R)]["stat_ok"] for R in R_LIST)
        conv_ok = all(c <= CONV_H for c in convs) and stats_ok
        dEs = [pa["row_R"][str(R)]["dEf"] for R in R_LIST]
        sgn = "DODATNI" if min(dEs) > 0 else (
            "NIEDODATNI(<=0)" if max(dEs) <= 0 else "MIESZANY")
        signs_conv.append((A, conv_ok, min(dEs), max(dEs), pa["dR"]))
        em("  A=%-6.4f: conv_h=%s (max %.2e); dE(R60,hf)=%+.6e "
           "dE(R120,hf)=%+.6e |dR|=%.2e => %s%s"
           % (A, "TAK" if conv_ok else "NIE", max(convs), dEs[0], dEs[1],
              pa["dR"], sgn, "" if conv_ok else " [NIEZBIEZNY h]"))
    em()
    conv_all = all(c for _, c, _, _, _ in signs_conv)
    pos_all = all(lo > 0 for _, c, lo, _, _ in signs_conv)
    neg_conv = [(A, lo) for A, c, lo, hi, _ in signs_conv
                if c and hi <= 0]
    em("KLASYFIKACJA Q-D2 (litera):")
    if conv_all and pos_all and not any_boundary:
        em("  Q-D2-COST: DeltaE_insert(A) > 0 dla wszystkich A != 1,"
           " zbieznie (h i R)")
        q = "Q-D2-COST"
    elif neg_conv:
        em("  Q-D2-CHANNEL: istnieje A z DeltaE_insert <= 0 zbieznie: %s"
           % neg_conv)
        q = "Q-D2-CHANNEL"
    else:
        em("  Q-D2-INCONCLUSIVE (brak zbieznosci h lub zaleznosc od R"
           " lub BREAKDOWN-BOUNDARY) -- NIE pozytyw; definicja wymaga"
           " poprawki (LOCK sec. 2), raport bez werdyktu znaku")
        q = "Q-D2-INCONCLUSIVE"
    if any_boundary:
        em("  BREAKDOWN-BOUNDARY: osobna kategoria deskryptywna"
           " (wystapil)")
    em()
    em("DESKRYPTYWNIE (obowiazkowe, LOCK sec. 2):")
    em("  monotonia DeltaE_insert(A) (R=120, h=0.0125): %s"
       % ["A=%.4f: %+.6e" % (A, per_A[lab]["row_R"]["120.0"]["dEf"])
          for A, lab in zip(A_FAMILY, A_LABELS)])
    em("  Odniesienie definicyjne do Q1-POS poprzednika"
       " (op-metametric-boundary, NIE reprodukcja): tam DeltaE_create =")
    em("  {-0.179 (soliton mu | proznia), +16156.6 (soliton | stan"
       " pusty)} poro-")
    em("  wnywalo ten sam profil wzgledem ROZNYCH tel (mieszane znaki,")
    em("  skladnik objetosciowy U(1)V(R) w drugim); tu DeltaE_insert"
       " jest roznica")
    em("  energii DWOCH stanow na TYM SAMYM tle/pudle/siatce/brzegu --")
    em("  czlon objetosciowy kasuje sie z definicji (kontrola |dR|"
       " wyzej);")
    em("  sektor rowniez inny: tam kanoniczny tachionowy (U'=g^6(1-g)),")
    em("  tu wlasciwa para (w, V_M9.1'') z Q-A-PASS.")
    em()
    em("WERDYKT: %s" % q)
    with open(BASE + "Phase3_output.txt", "w") as f:
        f.write("\n".join(out) + "\n")
    with open(RESDIR + "verdict.json", "w") as f:
        json.dump(dict(verdict=q, per_A=per_A), f, indent=1)
    print("zapisano:", BASE + "Phase3_output.txt")


# ------------------------------------------------------------------ main
if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "all"
    if mode in ("gates", "all"):
        gates()
    if mode in ("matrix", "all"):
        matrix()
    if mode in ("verdict", "all"):
        verdict()
