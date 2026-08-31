#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
op-bath-two-sectors (Phase 1) -- zmierzone (A, delta) ogonow + drabina d*.

LOCK: Phase0_balance.md sec. 3 Phase 1 (P1a-P1c); decyzje implementacyjne:
Phase_method_decisions.md (FROZEN przed uruchomieniem).

P1a: fit (g-1)*r ~ B cos r + C sin r, okna [50,150] i [120,260] (PRIMARY
     dla Phase 2: [120,260], R=450); delta = atan2(C,B); R-kontrola
     R=300 vs 450 (gate |d delta| <= 3.6 deg = 1% z 360); R2 >= 0.999.
     Rozjazd okien > 3 deg => flaga WINDOW-SENSITIVE.
P1b: E_int(d) = -Ai*Aj*exp(-(ki+kj)d/2)*cos(d - di - dj)/d, d in (2,60);
     pierwsze minimum d*, odstepy drabiny; predykcja LOCKED: 2pi +- 5%.
P1c: kontrola negatywna Yukawa E = -Ai*Aj*exp(-d)/d: wymagane 0 minimow;
     minimum => STOP i debug PRZED interpretacja.

Modele (LOCK sec. 2):
  M-P (PRIMARY): g'' + (2/r)g' + (2/g)g'^2 = g^2(1-g)   [EL Noty kanon.]
  M-L (CROSS):   F(g)g'' + (2/r)g' = g^2(1-g), F = 1+2*a_eff*ln g,
                 a_eff = 2/(1+eta_K*(g-1)^2), eta_K = 181/15  [INPUT]
Gatunki: g0_e=0.90548 [INPUT r21/phi-FP], g0_mu=phi*g0_e,
  M-P tau: g0=1.5696 [INPUT Q_K=3/2] + wrazliwosc g0*(1+-0.02).

Kod solverow: reuse N1_provenance_check.py (solve_m2 / solve_p131).
Wyniki maszynowe: Phase1_results.json (jedyne zrodlo (A,delta,kappa,d*)
dla Phase 2 -- LOCK sec. 4.3).
"""
import json
import numpy as np
from scipy.integrate import solve_ivp

PHI = (1 + np.sqrt(5)) / 2
ETA_K = 181.0 / 15.0                    # INPUT
G0_E = 0.90548                          # INPUT (kalibracja r21/phi-FP)
G0_MU = PHI * G0_E
G0_TAU = 1.5696                         # INPUT (Q_K=3/2)
WIN1 = (50.0, 150.0)
WIN2 = (120.0, 260.0)                   # PRIMARY
R_MAIN = 450.0
R_CTRL = 300.0
GATE_R2 = 0.999
GATE_RSTAB_DEG = 3.6                    # 1% z 360 deg
FLAG_WINDOW_DEG = 3.0
D_LO, D_HI, D_N = 2.0, 60.0, 580001


def solve_mp(g0, rmax):
    """M-P: g'' + (2/r)g' = g^2(1-g) - (2/g)g'^2 (reuse solve_m2, N1)."""
    def rhs(r, y):
        g, p = y
        if g <= 1e-12:
            return [p, 0]
        s = g * g * (1 - g) - (2.0 / g) * p * p
        if r < 1e-10:
            return [p, s / 3]
        return [p, s - 2 / r * p]

    def ev(r, y):
        return 100 - abs(y[0])
    ev.terminal = True
    s = solve_ivp(rhs, [0.01, rmax], [g0, 0.0], method='RK45',
                  rtol=1e-11, atol=1e-13, max_step=0.05,
                  events=[ev], dense_output=True)
    r = np.linspace(0.01, min(s.t[-1], rmax), int(rmax * 70))
    return r, s.sol(r)[0]


def solve_ml(g0, rmax):
    """M-L: F(g)g'' + (2/r)g' = g^2(1-g) (reuse solve_p131, N1)."""
    def fk(g):
        a = 2.0 / (1 + ETA_K * (g - 1) ** 2)
        return 1 + 2 * a * np.log(g) if g > 0 else -1e30

    def Vp(g):
        return g ** 2 * (1 - g)
    fg0 = fk(g0)
    if abs(fg0) < 1e-15:
        return None, None
    c2 = Vp(g0) / (3 * fg0)
    rs = 0.01

    def rhs(r, y):
        g, p = y
        if g <= 1e-15:
            return [p, 0]
        fg = fk(g)
        if abs(fg) < 1e-10:
            return [p, 0]
        if r < 1e-10:
            return [p, Vp(g) / fg / 3]
        return [p, (Vp(g) - 2 / r * p) / fg]

    def ev(r, y):
        return 100 - abs(y[0])
    ev.terminal = True
    s = solve_ivp(rhs, [rs, rmax], [g0 + c2 * rs ** 2, 2 * c2 * rs],
                  method='RK45', rtol=1e-11, atol=1e-13,
                  max_step=0.05, events=[ev], dense_output=True)
    r = np.linspace(rs, min(s.t[-1], rmax), int(rmax * 70))
    return r, s.sol(r)[0]


def fit_window(r, g, lo, hi):
    """(g-1)*r ~ B cos r + C sin r; zwraca (A, delta_rad, R2)."""
    m = (r >= lo) & (r <= hi)
    rf, tl = r[m], (g[m] - 1) * r[m]
    if len(rf) < 10:
        return np.nan, np.nan, np.nan
    M = np.column_stack([np.cos(rf), np.sin(rf)])
    coeff, _, _, _ = np.linalg.lstsq(M, tl, rcond=None)
    B, C = coeff
    resid = tl - M @ coeff
    ss_res = float(np.sum(resid ** 2))
    ss_tot = float(np.sum((tl - tl.mean()) ** 2))
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
    return float(np.hypot(B, C)), float(np.arctan2(C, B)), r2


def wrap_deg(x):
    return (x + 180.0) % 360.0 - 180.0


def find_minima(fun, lo=D_LO, hi=D_HI, n=D_N):
    d = np.linspace(lo, hi, n)
    E = fun(d)
    h = d[1] - d[0]
    mins = []
    for i in range(1, n - 1):
        if E[i] < E[i - 1] and E[i] < E[i + 1]:
            vpp = (E[i + 1] - 2 * E[i] + E[i - 1]) / h ** 2
            if vpp > 0:
                mins.append((float(d[i]), float(E[i]), float(vpp)))
    return mins


def measure(model_tag, solver, tag, g0):
    """Pelny pomiar P1a dla jednego gatunku; None gdy kolaps/brak ogona."""
    r4, g4 = solver(g0, R_MAIN)
    if r4 is None or r4[-1] < WIN2[1]:
        rend = 0.0 if r4 is None else float(r4[-1])
        print("  %s %-3s g0=%.5f  KOLAPS/NIE DOBIEGL (r_end=%.1f)"
              % (model_tag, tag, g0, rend))
        return None
    A1, d1, r2_1 = fit_window(r4, g4, *WIN1)
    A2, d2, r2_2 = fit_window(r4, g4, *WIN2)
    r3, g3 = solver(g0, R_CTRL)
    if r3 is None or r3[-1] < WIN2[1]:
        print("  %s %-3s R-kontrola nie dobiegla" % (model_tag, tag))
        return None
    A2c, d2c, _ = fit_window(r3, g3, *WIN2)
    ddeg_R = abs(wrap_deg(np.degrees(d2 - d2c)))
    ddeg_W = abs(wrap_deg(np.degrees(d1 - d2)))
    dA_R = abs(A2 - A2c) / A2 if A2 > 0 else np.nan
    kappa = np.log(A1 / A2) / 90.0 if (A1 > 0 and A2 > 0) else np.nan
    kappa_used = kappa if abs(kappa) > 1e-3 else 0.0
    ws = ddeg_W > FLAG_WINDOW_DEG
    ok_r2 = (r2_1 >= GATE_R2) and (r2_2 >= GATE_R2)
    ok_rs = ddeg_R <= GATE_RSTAB_DEG
    print("  %s %-3s g0=%.5f:" % (model_tag, tag, g0))
    print("    [50,150]:  A=%.6f  delta=%+8.2f deg  R2=%.7f"
          % (A1, np.degrees(d1), r2_1))
    print("    [120,260]: A=%.6f  delta=%+8.2f deg  R2=%.7f  (PRIMARY)"
          % (A2, np.degrees(d2), r2_2))
    print("    R-kontrola (300 vs 450): |d delta|=%.4f deg (gate<=3.6) "
          "dA/A=%.2e -> %s" % (ddeg_R, dA_R, "PASS" if ok_rs else "FAIL"))
    print("    rozjazd okien: %.2f deg %s;  kappa_eff=%+.2e -> uzyte "
          "kappa=%.4f" % (ddeg_W,
                          "-> FLAGA WINDOW-SENSITIVE" if ws else "(<=3, ok)",
                          kappa, kappa_used))
    print("    gate R2>=0.999: %s" % ("PASS" if ok_r2 else "FAIL"))
    return dict(g0=g0, A=A2, delta_rad=d2, delta_deg=np.degrees(d2),
                A_w1=A1, delta_w1_deg=np.degrees(d1), r2_w1=r2_1,
                r2_w2=r2_2, ddelta_R_deg=ddeg_R, dA_R=dA_R,
                window_shift_deg=ddeg_W, window_sensitive=bool(ws),
                kappa_eff=float(kappa), kappa_used=float(kappa_used),
                gate_r2=bool(ok_r2), gate_rstab=bool(ok_rs))


print("=" * 78)
print("Phase 1: zmierzone (A, delta) + drabina d*.  LOCK: Phase0 sec. 3.")
print("INPUT: g0_e=0.90548 (r21/phi-FP), g0_mu=phi*g0_e=%.5f," % G0_MU)
print("       g0_tau=1.5696 (Q_K=3/2), eta_K=181/15 (M-L).")
print("Okno PRIMARY dla Phase 2: [120,260], R=450 (method_decisions p.2).")
print("=" * 78)

results = {"MP": {}, "ML": {}}

print("\n--- P1a, M-P (potegowy, PRIMARY): g''+(2/r)g'+(2/g)g'^2=g^2(1-g) ---")
for tag, g0 in [("e", G0_E), ("mu", G0_MU), ("tau", G0_TAU)]:
    res = measure("M-P", solve_mp, tag, g0)
    if res:
        results["MP"][tag] = res

print("\n--- M-P tau: kontrola wrazliwosci +-2% (method_decisions P1 p.8) ---")
for tag, g0 in [("tau-2%", 0.98 * G0_TAU), ("tau+2%", 1.02 * G0_TAU)]:
    measure("M-P", solve_mp, tag, g0)   # deskryptywnie, nie wchodzi do JSON

print("\n--- P1a, M-L (logarytmiczny, CROSS-CHECK): F g''+(2/r)g'=g^2(1-g) ---")
for tag, g0 in [("e", G0_E), ("mu", G0_MU)]:
    res = measure("M-L", solve_ml, tag, g0)
    if res:
        results["ML"][tag] = res

# ---------------------------------------------------------------- P1b + P1c
PAIRS = [("ee", "e", "e"), ("emu", "e", "mu"), ("mumu", "mu", "mu")]
TWO_PI = 2 * np.pi
ladders = {"MP": {}, "ML": {}}
p1c_clean = True

for mod in ("MP", "ML"):
    print("\n--- P1b/P1c, model %s: E_int(d) = -Ai Aj e^{-(ki+kj)d/2} "
          "cos(d-di-dj)/d ---" % mod)
    for pname, si, sj in PAIRS:
        if si not in results[mod] or sj not in results[mod]:
            print("  %s %-4s: brak gatunku -> pomijam" % (mod, pname))
            continue
        Ri, Rj = results[mod][si], results[mod][sj]
        Ai, Aj = Ri["A"], Rj["A"]
        di, dj = Ri["delta_rad"], Rj["delta_rad"]
        kk = 0.5 * (Ri["kappa_used"] + Rj["kappa_used"])

        def E_osc(d, Ai=Ai, Aj=Aj, di=di, dj=dj, kk=kk):
            return -Ai * Aj * np.exp(-2 * kk * d / 2.0) \
                * np.cos(d - di - dj) / d

        def E_yuk(d, Ai=Ai, Aj=Aj):
            return -Ai * Aj * np.exp(-d) / d

        mins = find_minima(E_osc)
        nyuk = len(find_minima(E_yuk))
        dstars = [m[0] for m in mins[:6]]
        gaps = np.diff(dstars) if len(dstars) >= 2 else np.array([])
        dev = (np.abs(gaps - TWO_PI) / TWO_PI * 100
               if len(gaps) else np.array([np.nan]))
        ok2pi = bool(len(gaps) and np.nanmax(dev) <= 5.0)
        print("  %s %-4s: %d minimow; d*1..3 = %s" % (
            mod, pname, len(mins),
            ["%.4f" % x for x in dstars[:3]]))
        if len(gaps):
            print("    odstepy = %s; max odch. od 2pi = %.2f%% "
                  "(predykcja LOCKED 2pi+-5%%: %s)"
                  % (["%.4f" % x for x in gaps], float(np.max(dev)),
                     "PASS" if ok2pi else "FAIL"))
        print("    P1c Yukawa (te same A): %d minimow -> %s"
              % (nyuk, "CZYSTA" if nyuk == 0 else
                 "SKAZONA -> STOP I DEBUG"))
        if nyuk != 0:
            p1c_clean = False
        ladders[mod][pname] = dict(
            dstars=dstars, gaps=[float(x) for x in gaps],
            max_dev_2pi_pct=(float(np.max(dev)) if len(gaps) else None),
            pred_2pi_pass=ok2pi, n_minima=len(mins), yukawa_minima=nyuk,
            kappa_pair=float(kk))

# ------------------------------------------------------------------- GATE
all_meas = [v for m in results.values() for v in m.values()]
gate_fits = all(v["gate_r2"] for v in all_meas)
gate_rst = all(v["gate_rstab"] for v in all_meas)
gate = gate_fits and gate_rst and p1c_clean
print("\n" + "=" * 78)
print("GATE PHASE 1 (LOCK: fity R2>=0.999 i R-stabilne; kontrola P1c "
      "czysta):")
print("  R2 wszystkie >= 0.999: %s" % ("PASS" if gate_fits else "FAIL"))
print("  R-stabilnosc |d delta|<=3.6 deg: %s"
      % ("PASS" if gate_rst else "FAIL"))
print("  P1c (Yukawa, wszystkie pary/modele): %s"
      % ("CZYSTA" if p1c_clean else "SKAZONA"))
print("  => GATE: %s" % ("PASS" if gate else "FAIL -> STOP przed Phase 2"))
print("=" * 78)

out = dict(lock="Phase0_balance.md", primary_window=[120, 260],
           primary_model_for_phase2="MP",
           inputs=dict(g0_e=G0_E, g0_mu=G0_MU, g0_tau=G0_TAU,
                       eta_K=ETA_K, note="Q_K=3/2, r21/phi-FP = INPUT"),
           species=results, ladders=ladders,
           gate=dict(fits=gate_fits, rstab=gate_rst, p1c=p1c_clean,
                     overall=gate))
with open("TGP/TGP_v1/research/op-bath-two-sectors-2026-08-23/"
          "Phase1_results.json", "w") as f:
    json.dump(out, f, indent=1)
print("zapisano Phase1_results.json")
