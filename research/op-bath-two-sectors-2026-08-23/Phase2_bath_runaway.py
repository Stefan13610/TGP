#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
op-bath-two-sectors (Phase 2, Q1) -- runaway w kapieli o skonczonej
gestosci, na fazach ZMIERZONYCH (Phase1_results.json; LOCK sec. 4.3).

LOCK: Phase0_balance.md sec. 3 Phase 2; decyzje: Phase_method_decisions.md
(FROZEN). Kod dynamiki/spektrum: reuse Phase3_nonlinear_evolution.py (#63,
V3): M0-f_eps, eps=0.2, RK4 na dokladnym ukladzie hamiltonowskim,
gate |dE|/E <= 1e-6.

P2a (gate): reprodukcja runaway izolacji (#63: mu, t*=3.62) na R=60,
  N=4000; a=+-0.01, dt=0.004 i 0.002. Gate FAIL => STOP Q1.
P2b: konfiguracje {ee, emu, mumu}; d in {d1, d2, d3, 0.5 d1, 1.5 d1}
  z drabiny M-P (PRIMARY); kapiel: dg_bath(r) = N_c A_s e^{-k_s d}
  cos(d-delta_s)/d * sin(r)/r (usredniona po sferze sasiadow, N_c=12);
  komorki R_cell in {d/2, d} x siatki h in {0.015, 0.0075} (spektra);
  ewolucja (h=0.015): amp in {+0.01, -0.01, 0 (kontrola sloshingu)}
  wzdluz modu runaway komorki, t_end = 3 t*_izo.
P2c: sprektrum prozni g=1 w TEJ SAMEJ komorce; mod klasyfikowany jako
  CELL-MODE gdy |lam - lam_vac| < 0.05 max(1,|lam_vac|) ORAZ
  lokalizacja(r < min(2, 0.8 R_cell)) < 0.5. omega^2_min(clean) =
  najnizsza wartosc wlasna F-wazona nie-CELL (kryterium werdyktu);
  raw raportowane obok.

Zbieznosc "do 2 cyfr": rozrzut omega^2_min(clean) po 4 kombinacjach
{siatka x komorka} <= 0.01 * max(|srednia|, 0.1).

INPUT: g0_e=1.24915 (galaz f_eps/a3d), phi-FP, eps=0.2, N_c=12;
(A, delta, kappa) WYLACZNIE z Phase 1 (json).
"""
import json
import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.integrate import solve_ivp

CYC = "TGP/TGP_v1/research/op-bath-two-sectors-2026-08-23/"
PHI_GOLD = (1 + 5 ** 0.5) / 2
G0_E = 1.24915                       # INPUT (galaz f_eps/a3d)
G0_MU = PHI_GOLD * G0_E
ALPHA = 2.0
EPS = 0.2
DT_BASE = 0.004
GATE = 1e-6
N_C = 12                             # INPUT (koordynacja WS/FCC)
H_LIST = (0.015, 0.0075)
LAM_LOCK = -1.389

f_log = lambda g: 1 + 2 * ALPHA * np.log(g)
fp_log = lambda g: 2 * ALPHA / g
fpp_log = lambda g: -2 * ALPHA / g ** 2
W_fun = lambda g: g ** 3 / 3.0 - g ** 4 / 4.0
W1_fun = lambda g: g ** 2 * (1 - g)
W2_fun = lambda g: 2 * g - 3 * g ** 2


def make_feps(eps):
    def s(g):
        fv = f_log(g)
        return np.sqrt(fv ** 2 + eps ** 2)
    F = lambda g: 0.5 * (f_log(g) + s(g))
    Fp = lambda g: 0.5 * fp_log(g) * (1 + f_log(g) / s(g))
    Fpp = lambda g: 0.5 * (fpp_log(g) * (1 + f_log(g) / s(g))
                           + fp_log(g) ** 2 * eps ** 2 / s(g) ** 3)
    return F, Fp, Fpp


F_EPS, FP_EPS, FPP_EPS = make_feps(EPS)


def profile_soft(g0, R, N):
    """#62/#63 soliton_profile_soft (omega=0), verbatim convention."""
    def rhs(r_, y):
        g, gp = y
        g = max(g, 1e-3)
        drv = W1_fun(g) - 0.5 * FP_EPS(g) * gp ** 2
        if r_ < 1e-10:
            return [gp, drv / (3 * F_EPS(g))]
        return [gp, drv / F_EPS(g) - 2 * gp / r_]
    h = R / N
    r_grid = (np.arange(N) + 0.5) * h
    sol = solve_ivp(rhs, [1e-6, R + h], [g0, 0.0], method='DOP853',
                    max_step=0.02, rtol=1e-10, atol=1e-13,
                    dense_output=True)
    r_end = float(sol.t[-1])
    sel = r_grid <= r_end
    r_g = r_grid[sel]
    vals = sol.sol(r_g)
    g_arr, gp_arr = vals[0], vals[1]
    d2 = np.array([rhs(r_g[j], [g_arr[j], gp_arr[j]])[1]
                   for j in range(len(r_g))])
    return r_g, g_arr, gp_arr, d2, r_end, len(r_g) == N


def build_tridiag_gen(r, g, gp, d2):
    """#62/#63 assembly (weight-1 problem on u = r v), verbatim."""
    N = len(r)
    h = r[1] - r[0]
    F_nodes = F_EPS(g)
    F_mid = F_EPS(0.5 * (g[:-1] + g[1:]))
    Q = W2_fun(g) - 0.5 * FPP_EPS(g) * gp ** 2 \
        - FP_EPS(g) * (d2 + 2 * gp / r)
    r_mid = 0.5 * (r[:-1] + r[1:])
    a = r_mid ** 2 * F_mid
    diag = np.empty(N)
    diag[0] = a[0] / h ** 2
    diag[1:-1] = (a[:-1] + a[1:]) / h ** 2
    diag[-1] = (a[-1] + r[-1] ** 2 * F_nodes[-1]) / h ** 2
    diag += r ** 2 * Q
    off = -a / h ** 2
    return diag / r ** 2, off / (r[:-1] * r[1:])


def spectra(r, g, gp, d2, nmodes=10):
    """Zwraca (lam1, vec1_u, lamF, vecF_u): waga-1 i F-wazone (u = r v)."""
    d, e = build_tridiag_gen(r, g, gp, d2)
    k = min(nmodes - 1, len(r) - 2)
    v1, u1 = eigh_tridiagonal(d, e, select='i', select_range=(0, k))
    sF = np.sqrt(F_EPS(g))
    dw = d / F_EPS(g)
    ew = e / (sF[:-1] * sF[1:])
    vw, uw = eigh_tridiagonal(dw, ew, select='i', select_range=(0, k))
    uw_u = uw / sF[:, None]          # powrot do u z w~ = B^(1/2) u
    return v1, u1, vw, uw_u


def localization(r, u, r_loc):
    w = u ** 2
    tot = float(np.sum(w))
    if tot <= 0:
        return 0.0
    return float(np.sum(w[r < r_loc])) / tot


def classify_clean(lams, us, lams_vac, r, r_loc):
    """Zwraca (lam_clean_min, idx, flags); CELL-MODE wg method_decisions."""
    flags = []
    for i, lam in enumerate(lams):
        near_vac = bool(np.any(np.abs(lams_vac - lam)
                               < 0.05 * np.maximum(1.0, np.abs(lams_vac))))
        loc = localization(r, us[:, i], r_loc)
        cell = near_vac and (loc < 0.5)
        flags.append((float(lam), loc, cell))
    clean = [(lam, i) for i, (lam, loc, cell) in enumerate(flags)
             if not cell]
    if not clean:
        return None, None, flags
    lam_min, idx = min(clean)
    return float(lam_min), idx, flags


class SoftWallDynamics:
    def __init__(self, r):
        self.r = r
        self.h = r[1] - r[0]
        self.r2 = r ** 2
        self.rm2 = (0.5 * (r[:-1] + r[1:])) ** 2

    def energy(self, g, pi):
        gm = 0.5 * (g[:-1] + g[1:])
        dg = np.diff(g) / self.h
        Ek = np.sum(self.r2 * (pi ** 2 / (2 * F_EPS(g)) + W_fun(g))) * self.h
        Es = 0.5 * np.sum(self.rm2 * F_EPS(gm) * dg ** 2) * self.h
        return Ek + Es

    def rhs(self, g, pi):
        h, r2, rm2 = self.h, self.r2, self.rm2
        gm = 0.5 * (g[:-1] + g[1:])
        dg = np.diff(g) / h
        Fg = F_EPS(g)
        gdot = pi / Fg
        dH = np.zeros_like(g)
        dH += h * r2 * (-pi ** 2 * FP_EPS(g) / (2 * Fg ** 2) + W1_fun(g))
        t_flux = rm2 * F_EPS(gm) * dg
        t_quad = 0.25 * h * rm2 * FP_EPS(gm) * dg ** 2
        dH[:-1] += -t_flux + t_quad
        dH[1:] += t_flux + t_quad
        pidot = -dH / (h * r2)
        return gdot, pidot

    def rk4_step(self, g, pi, dt):
        k1g, k1p = self.rhs(g, pi)
        k2g, k2p = self.rhs(g + 0.5 * dt * k1g, pi + 0.5 * dt * k1p)
        k3g, k3p = self.rhs(g + 0.5 * dt * k2g, pi + 0.5 * dt * k2p)
        k4g, k4p = self.rhs(g + dt * k3g, pi + dt * k3p)
        return (g + dt / 6 * (k1g + 2 * k2g + 2 * k3g + k4g),
                pi + dt / 6 * (k1p + 2 * k2p + 2 * k3p + k4p))


G_FLOOR = 0.01   # korekta 1 (Phase2_correction_note.md): podloga modelu


def evolve(r, g_init, dt, t_end, g_ref, bg_amp):
    """Zwraca (t_break, gate, max||dg||_inf/bg, t_at_max).
    Korekta 1: wyjscie z dziedziny (min g <= G_FLOOR lub non-finite)
    sprawdzane W KAZDYM kroku RK4; energia tylko przed-breakdown."""
    dyn = SoftWallDynamics(r)
    g = g_init.copy()
    pi = np.zeros_like(g)
    E0 = dyn.energy(g, pi)
    nstep = int(round(t_end / dt))
    every = max(1, int(round(0.02 / dt)))
    emax, dmax, t_dmax, t_break = 0.0, 0.0, 0.0, None
    for k in range(nstep + 1):
        if (not np.all(np.isfinite(g))) or float(np.min(g)) <= G_FLOOR:
            t_break = k * dt
            break
        if k % every == 0:
            dn = float(np.max(np.abs(g - g_ref))) / bg_amp
            if dn > dmax:
                dmax, t_dmax = dn, k * dt
            E = dyn.energy(g, pi)
            if np.isfinite(E):
                emax = max(emax, abs(E - E0) / abs(E0))
        if k < nstep:
            g, pi = dyn.rk4_step(g, pi, dt)
    return t_break, emax, dmax, t_dmax


def bath_terms(r, c):
    """c*sin(r)/r i pochodne analityczne."""
    dg = c * np.sin(r) / r
    dgp = c * (r * np.cos(r) - np.sin(r)) / r ** 2
    dg2 = c * (-np.sin(r) / r - 2 * np.cos(r) / r ** 2
               + 2 * np.sin(r) / r ** 3)
    return dg, dgp, dg2


# ===================================================================== RUN
print("=" * 78)
print("Phase 2 (Q1): runaway w kapieli.  LOCK: Phase0 sec. 3 Phase 2.")
print("Model dynamiczny: M0-f_eps eps=0.2 (#63 V3).  INPUT: g0_e=1.24915")
print("(galaz f_eps/a3d), phi-FP, N_c=12.  (A,delta,kappa,d*): WYLACZNIE")
print("Phase1_results.json, model M-P PRIMARY, okno [120,260].")
print("=" * 78)

with open(CYC + "Phase1_results.json") as f:
    P1 = json.load(f)
assert P1["gate"]["p1c"] and P1["gate"]["rstab"], "Phase1 gate: patrz ruling"
SP = P1["species"]["MP"]
LAD = P1["ladders"]["MP"]
print("wczytano Phase1: A_e=%.6f d_e=%+.2f deg; A_mu=%.6f d_mu=%+.2f deg"
      % (SP["e"]["A"], SP["e"]["delta_deg"],
         SP["mu"]["A"], SP["mu"]["delta_deg"]))

# ------------------------------------------------------------ P2a baseline
print("\n" + "-" * 78)
print("[P2a] baseline izolacji, R=60, N=4000 (siatka #63)")
BASE = {}
p2a_fail = False
for sp, g0 in [("mu", G0_MU), ("e", G0_E)]:
    r, g, gp, d2, r_end, full = profile_soft(g0, 60.0, 4000)
    if not full:
        print("  %s: profil NIEKOMPLETNY (r_end=%.2f) -- brak reprezentanta"
              % (sp, r_end))
        BASE[sp] = None
        continue
    v1, u1, vw, uw = spectra(r, g, gp, d2)
    lam1, lamw = float(v1[0]), float(vw[0])
    print("  %s: g0=%.5f  lam_min(w1)=%.4f (#63: %.3f dla mu)  "
          "lam_min(F)=%.4f" % (sp, g0, lam1, LAM_LOCK, lamw))
    u = u1[:, 0]
    v = u / r
    h = r[1] - r[0]
    v /= np.sqrt(np.sum(r ** 2 * v ** 2) * h)
    if v[np.argmax(np.abs(v))] < 0:
        v = -v
    bg_amp = float(np.max(np.abs(g - 1.0)))
    tstars, gates = [], []
    for amp, dt in [(0.01, DT_BASE), (0.01, DT_BASE / 2),
                    (-0.01, DT_BASE)]:
        tb, egate, dmax, tdm = evolve(r, g + amp * v, dt, 12.0, g, bg_amp)
        ok = egate < GATE
        gates.append(ok)
        print("    a=%+.2f dt=%.4f: gate=%.2e %s; %s" % (
            amp, dt, egate, "PASS" if ok else "FAIL",
            ("BREAKDOWN t*=%.2f" % tb) if tb is not None else
            ("bez ucieczki do t=12 (max||dg||/tlo=%.1f%% @t=%.1f)"
             % (100 * dmax, tdm))))
        if tb is not None:
            tstars.append(tb)
    if sp == "mu" and not all(gates):
        # zalockowany gate P2a = reprodukcja #63 V3 (mu); e = dodatek
        p2a_fail = True
    elif sp == "e" and not all(gates):
        print("    UWAGA: gate komparatora e FAIL (dodatek, poza gate P2a)")
    BASE[sp] = dict(r=r, g=g, gp=gp, d2=d2, lam1=lam1, lamw=lamw,
                    tstar=(min(tstars) if tstars else None),
                    bg_amp=bg_amp, g0=g0)
    print("    t*_izo(%s) = %s" % (sp, ("%.2f" % min(tstars))
                                   if tstars else "BRAK (nie ucieka)"))

if p2a_fail:
    print("\nP2a GATE FAIL -> STOP Q1 (kod niewazny; LOCK drzewo sec. 6).")
    raise SystemExit(0)
print("  P2a: gate energii PASS wszystkie przebiegi.")

# ------------------------------------------------------------ P2b + P2c
CONFIGS = [("ee", "e", "e"), ("emu", "mu", "e"), ("mumu", "mu", "mu")]
results = {}
for cname, central, nb in CONFIGS:
    print("\n" + "-" * 78)
    base = BASE.get(central)
    lad = LAD["emu" if cname == "emu" else cname]
    if base is None:
        print("[%s] central=%s NIEREPREZENTOWALNY w f_eps -> konfiguracja "
              "raportowana jako niereprezentowalna (deskryptywnie)" %
              (cname, central))
        results[cname] = "NOT-REPRESENTABLE"
        continue
    ds = lad["dstars"]
    d_list = [ds[0], ds[1], ds[2], 0.5 * ds[0], 1.5 * ds[0]]
    A_s = SP[nb]["A"]
    del_s = SP[nb]["delta_rad"]
    kap_s = SP[nb]["kappa_used"]
    tstar = base["tstar"]
    t_end = 3 * tstar if tstar else 12.0
    print("[%s] central=%s (g0=%.5f, t*_izo=%s), sasiedzi=%s "
          "(A=%.4f, delta=%+.1f deg, kappa=%.3f), N_c=%d" %
          (cname, central, base["g0"],
           ("%.2f" % tstar) if tstar else "brak", nb, A_s,
           np.degrees(del_s), kap_s, N_C))
    print("  d punkty (LOCKED): %s" % ["%.3f" % d for d in d_list])
    conf_pts = []
    for d in d_list:
        c = N_C * A_s * np.exp(-kap_s * d) * np.cos(d - del_s) / d
        print("  -- d=%.3f: c_bath=%+.4f (=N_c A cos(d-delta)/d)" % (d, c))
        lam_grid = {}
        pt = dict(d=d, c=c, combos={}, evo={})
        for R_cell in (0.5 * d, d):
            for h in H_LIST:
                N = max(24, int(round(R_cell / h)))
                r, g, gp, d2, r_end, full = profile_soft(
                    base["g0"], R_cell, N)
                if not full:
                    print("     R=%.2f h=%.4f: profil niekompletny -> skip"
                          % (R_cell, h))
                    continue
                bg, bgp, bg2 = bath_terms(r, c)
                gB, gBp, gB2 = g + bg, gp + bgp, d2 + bg2
                if float(np.min(gB)) <= 0.0:
                    print("     R=%.2f h=%.4f: min(g_bath)=%.3f <= 0 -> "
                          "DOMAIN-VIOLATION (poza dziedzina log)" %
                          (R_cell, h, float(np.min(gB))))
                    pt["combos"][(R_cell, h)] = "DOMAIN"
                    continue
                v1, u1, vw, uw = spectra(r, gB, gBp, gB2)
                # P2c: proznia w tej samej komorce
                ones = np.ones_like(r)
                zer = np.zeros_like(r)
                v1v, _, vwv, _ = spectra(r, ones, zer, zer)
                r_loc = min(2.0, 0.8 * R_cell)
                lamc, idx, flags = classify_clean(
                    vw, uw, np.array(vwv), r, r_loc)
                nneg_vac = int(np.sum(np.array(vwv) < 0))
                lam_grid[(R_cell, h)] = lamc
                pt["combos"][(R_cell, h)] = dict(
                    lam_raw=float(vw[0]), lam_clean=lamc,
                    lam1_raw=float(v1[0]),
                    vac_min=float(vwv[0]), nneg_vac=nneg_vac,
                    nneg_bath=int(np.sum(np.array(vw) < 0)))
                print("     R=%6.2f h=%.4f N=%5d: w2_raw=%+8.4f  "
                      "w2_clean=%s  vac_min=%+8.4f (N_neg_vac=%d, "
                      "floor(R/pi)=%d)" %
                      (R_cell, h, N, float(vw[0]),
                       ("%+8.4f" % lamc) if lamc is not None else "  BRAK ",
                       float(vwv[0]), nneg_vac, int(R_cell / np.pi)))
        # zbieznosc omega^2_clean po kombach
        vals = [v for v in lam_grid.values() if v is not None]
        conv = None
        if len(vals) == 4:
            spread = max(vals) - min(vals)
            mean = float(np.mean(vals))
            conv = spread <= 0.01 * max(abs(mean), 0.1)
            print("     zbieznosc w2_clean: srednia=%+.4f rozrzut=%.4f -> %s"
                  % (mean, spread, "ZBIEZNE" if conv else "NIEZBIEZNE"))
            pt["w2_mean"] = mean
        pt["converged"] = conv
        # ewolucja (h=0.015, oba rozmiary komorki)
        for R_cell in (0.5 * d, d):
            h = 0.015
            N = max(24, int(round(R_cell / h)))
            r, g, gp, d2, r_end, full = profile_soft(base["g0"], R_cell, N)
            if not full:
                continue
            bg, bgp, bg2 = bath_terms(r, c)
            gB = g + bg
            if float(np.min(gB)) <= 0.0:
                pt["evo"][R_cell] = "DOMAIN"
                continue
            v1, u1, vw, uw = spectra(r, gB, gp + bgp, d2 + bg2)
            ones = np.ones_like(r)
            zer = np.zeros_like(r)
            v1v, _, _, _ = spectra(r, ones, zer, zer)
            r_loc = min(2.0, 0.8 * R_cell)
            lamc1, idx1, _ = classify_clean(
                v1, u1, np.array(v1v), r, r_loc)
            iu = idx1 if idx1 is not None else 0
            u = u1[:, iu]
            v = u / r
            hh = r[1] - r[0]
            v /= np.sqrt(np.sum(r ** 2 * v ** 2) * hh)
            if v[np.argmax(np.abs(v))] < 0:
                v = -v
            bg_amp = float(np.max(np.abs(gB - 1.0)))
            outs = []
            for amp in (0.01, -0.01, 0.0):
                tb, egate, dmax, tdm = evolve(r, gB + amp * v, DT_BASE,
                                              t_end, gB, bg_amp)
                outs.append((amp, tb, egate, dmax))
                print("     evo R=%6.2f a=%+5.2f: gate=%.1e %s; %s" % (
                    R_cell, amp, egate, "PASS" if egate < GATE else "FAIL",
                    ("BREAKDOWN t*=%.2f (t/t*_izo=%.2f)" %
                     (tb, tb / tstar if tstar else np.nan))
                    if tb is not None else
                    "bez ucieczki do t=%.1f (3 t*_izo); max||dg||/tlo=%.0f%%"
                    % (t_end, 100 * dmax)))
            pt["evo"][R_cell] = outs
        conf_pts.append(pt)
    results[cname] = conf_pts

# ------------------------------------------------------------------ WERDYKT
print("\n" + "=" * 78)
print("PODSUMOWANIE Q1 (kryteria LOCK sec. 3, doslownie):")
any_pass_pt = False
all_fail_pt = True
n_pts = 0
for cname, pts in results.items():
    if not isinstance(pts, list):
        print("  [%s] %s" % (cname, pts))
        all_fail_pt = False
        continue
    for pt in pts:
        n_pts += 1
        evo_runs = [o for v in pt["evo"].values() if isinstance(v, list)
                    for o in v if o[0] != 0.0]
        no_runaway = (len(evo_runs) > 0
                      and all(o[1] is None for o in evo_runs))
        esc2t = [o for o in evo_runs if o[1] is not None
                 and BASE.get("mu") and o[1] <= 2 * (BASE["mu"]["tstar"]
                                                     or 1e9)]
        w2 = pt.get("w2_mean")
        conv = pt.get("converged")
        pass_pt = bool(conv and w2 is not None and w2 > 0 and no_runaway)
        fail_pt = bool(w2 is not None and w2 < 0 and len(esc2t) > 0
                       and conv)
        if pass_pt:
            any_pass_pt = True
        if not fail_pt:
            all_fail_pt = False
        print("  [%s d=%.3f] w2_clean=%s conv=%s; ewolucja +-: %s => %s"
              % (cname, pt["d"],
                 ("%+.4f" % w2) if w2 is not None else "brak",
                 conv,
                 "bez ucieczki" if no_runaway else
                 "ucieczka/domain", "PASS-pt" if pass_pt else
                 ("FAIL-pt" if fail_pt else "nierozstrzygajacy")))
if any_pass_pt:
    verdict = "Q1-PASS"
elif all_fail_pt and n_pts > 0:
    verdict = "Q1-FAIL"
else:
    verdict = "Q1-INCONCLUSIVE (albo mieszane, albo niezbiezne/domain)"
print("\nWERDYKT: %s" % verdict)
print("(interpretacja liczb i kombinacji -> Phase_FINAL_close.md)")
print("=" * 78)
