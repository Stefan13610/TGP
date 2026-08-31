#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
op-bloch-chain-stability (Phase 4, WARUNKOWA -- uruchomiona, bo Phase 3
ZBIEZNA dla wszystkich policzonych tel) -- test nieliniowy.

LOCK: Phase0_balance.md sec. 3 Phase 4; decyzje FROZEN:
Phase_method_decisions.md sec. 1, 6.

Superkomorka 4d; perturbacja +- wzdluz modu minimalnego z Phase 3
(argmin k=0 wszedzie -> wspolmierny; mod NIE-translacyjny, regula FROZEN),
amplituda 0.01*||tlo||, ||tlo|| := max|g0-1|. Ewolucja: RK4 na
semi-dyskretnym hamiltonianie akcji kanonicznej z K_eps
(eps=0.2; kontrola eps=0.1 przy znaku +; kontrola dt/2 przy (+, 0.2));
t_end = 3*t*_ref = 10.86 (t*_ref = 3.62 [INPUT #63]).
UCIECZKA := (g<=0 lub niefinityczne) LUB max|g-g0| > 1.0*||tlo||.
Raport t_escape vs 2*t*_ref = 7.24 i 3*t*_ref = 10.86; gate energii
raportowany; wzrost sigma_fit vs sqrt(|omega^2_min|) (dokumentacja).
"""
import json
import numpy as np
from scipy.linalg import eigh

BETA = 1.0
GAMMA = 1.0
BASE = ("C:/Users/Mateusz/Documents/ObsydnianMain/TGP/TGP_v1/research/"
        "op-bloch-chain-stability-2026-08-31/")
NPZ = BASE + "Phase2_backgrounds.npz"
RES3 = BASE + "Phase3_results.json"
T_REF = 3.62
T_END = 3 * T_REF
DT_BASE = 0.004
AMP_REL = 0.01
M_SUPER = 4
GATE_E = 1e-6


def Keps(g, eps):
    K = g ** 4
    return 0.5 * (K + np.sqrt(K ** 2 + eps ** 2))


def Keps_prime(g, eps):
    K = g ** 4
    Kp = 4 * g ** 3
    return 0.5 * Kp * (1 + K / np.sqrt(K ** 2 + eps ** 2))


def U1_fun(g):
    return BETA * g ** 6 - GAMMA * g ** 7


def U_fun(g):
    return BETA * g ** 7 / 7 - GAMMA * g ** 8 / 8


def mode_k0(g, h, d):
    """Najnizszy NIE-translacyjny mod k=0 (regula FROZEN sec. 2/6)."""
    N = len(g)
    gp = (np.roll(g, -1) - np.roll(g, 1)) / (2 * h)
    w = g ** 4
    pot = w * (2 * BETA * g - 3 * GAMMA * g ** 2 + 2 * gp ** 2 / g ** 2)
    wmid = 0.5 * (w + np.roll(w, -1))
    A = np.zeros((N, N))
    idx = np.arange(N)
    A[idx, idx] = (wmid + np.roll(wmid, 1)) / h ** 2 + pot
    off = -wmid / h ** 2
    A[idx, (idx + 1) % N] += off
    A[(idx + 1) % N, idx] += off
    sw = np.sqrt(w)
    H = A / sw[:, None] / sw[None, :]
    vals, vecs = eigh(H, subset_by_index=(0, 3))
    for j in range(4):
        phi = vecs[:, j] / sw
        num = abs(float(np.sum(w * phi * gp) * h))
        den = np.sqrt(float(np.sum(w * phi ** 2) * h)
                      * float(np.sum(w * gp ** 2) * h))
        if den == 0 or num / den < 0.9:
            return float(vals[j]), phi, j
    raise RuntimeError("brak modu nie-translacyjnego w 4 najnizszych")


class Dyn:
    def __init__(self, h, eps):
        self.h, self.eps = h, eps

    def energy(self, g, pi):
        h, eps = self.h, self.eps
        gm = 0.5 * (g + np.roll(g, -1))
        dg = (np.roll(g, -1) - g) / h
        return float(np.sum(h * (pi ** 2 / (2 * Keps(g, eps)) + U_fun(g))
                            + 0.5 * h * Keps(gm, eps) * dg ** 2))

    def rhs(self, g, pi):
        h, eps = self.h, self.eps
        gm = 0.5 * (g + np.roll(g, -1))
        dg = (np.roll(g, -1) - g) / h
        Kg = Keps(g, eps)
        gdot = pi / Kg
        dH = h * (-pi ** 2 * Keps_prime(g, eps) / (2 * Kg ** 2) + U1_fun(g))
        t_flux = Keps(gm, eps) * dg
        t_quad = 0.25 * h * Keps_prime(gm, eps) * dg ** 2
        dH += -t_flux + t_quad
        dH += np.roll(t_flux + t_quad, 1)
        return gdot, -dH / h

    def rk4(self, g, pi, dt):
        k1g, k1p = self.rhs(g, pi)
        k2g, k2p = self.rhs(g + 0.5 * dt * k1g, pi + 0.5 * dt * k1p)
        k3g, k3p = self.rhs(g + 0.5 * dt * k2g, pi + 0.5 * dt * k2p)
        k4g, k4p = self.rhs(g + dt * k3g, pi + dt * k3p)
        return (g + dt / 6 * (k1g + 2 * k2g + 2 * k3g + k4g),
                pi + dt / 6 * (k1p + 2 * k2p + 2 * k3p + k4p))


def run(dyn, g0_sc, mode_sc, amp, dt, bg_amp):
    g = g0_sc + amp * mode_sc
    pi = np.zeros_like(g)
    E0 = dyn.energy(g, pi)
    nstep = int(round(T_END / dt))
    every = max(1, int(round(0.02 / dt)))
    emax = 0.0
    t_esc = None
    ts, prj = [], []
    wB = g0_sc ** 4
    nrm = float(np.sum(wB * mode_sc ** 2) * dyn.h)
    for kstep in range(nstep + 1):
        if kstep % every == 0:
            t = kstep * dt
            if (not np.all(np.isfinite(g))) or float(np.min(g)) <= 0.0:
                t_esc = t
                reason = "BREAKDOWN g<=0/nonfinite"
                break
            dgv = g - g0_sc
            ratio = float(np.max(np.abs(dgv))) / bg_amp
            ts.append(t)
            prj.append(float(np.sum(wB * mode_sc * dgv) * dyn.h) / nrm)
            E = dyn.energy(g, pi)
            if np.isfinite(E) and E0 != 0:
                emax = max(emax, abs(E - E0) / abs(E0))
            if ratio > 1.0:
                t_esc = t
                reason = "||dg||>100%% tla (ratio=%.2f)" % ratio
                break
        if kstep < nstep:
            g, pi = dyn.rk4(g, pi, dt)
    else:
        reason = "bez ucieczki do t_end"
    return np.array(ts), np.array(prj), emax, t_esc, reason


def fit_sigma(ts, prj, a0):
    aa = np.abs(prj)
    lo, hi = 3 * abs(a0), 20 * abs(a0)
    idx = [i for i in range(len(ts)) if lo < aa[i] < hi]
    if len(idx) < 5:
        return np.nan
    i0, i1 = idx[0], idx[-1]
    p = np.polyfit(ts[i0:i1 + 1], np.log(aa[i0:i1 + 1]), 1)
    return float(p[0])


print("=" * 78)
print("Phase 4 -- test nieliniowy (warunek: Phase 3 ZBIEZNA -> spelniony).")
print("Superkomorka 4d, perturbacja +-0.01||tlo|| wzdluz modu minimalnego;")
print("K_eps: eps=0.2 (+- ), kontrola eps=0.1 (+), kontrola dt/2 (+, 0.2);")
print("t*_ref=3.62; 2t*=7.24; t_end=3t*=10.86. UCIECZKA: g<=0 lub")
print("max|g-g0| > 100%% tla.")
print("=" * 78)

data = np.load(NPZ)
res3 = json.load(open(RES3))
verdict_rows = []
for label, r3 in res3.items():
    dn, tag = label.split("/")
    d = r3["d"]
    g0 = data["%s__%s__N400" % (dn, tag)]
    h = d / 400
    lam0, phi, jidx = mode_k0(g0, h, d)
    print("\n--- tlo %s (d=%.4f): mod minimalny k=0, lambda=%.6f "
          "(Phase3: %.6f), indeks %d ---"
          % (label, d, lam0, r3["w2min"], jidx))
    phi = np.real(phi)
    phi /= float(np.max(np.abs(phi)))
    bg_amp = float(np.max(np.abs(g0 - 1)))
    amp = AMP_REL * bg_amp
    g0_sc = np.tile(g0, M_SUPER)
    mode_sc = np.tile(phi, M_SUPER)
    sig_lin = np.sqrt(-r3["w2min"]) if r3["w2min"] < 0 else np.nan
    print("  ||tlo|| = %.6f; amplituda = %.6f; sigma_lin = "
          "sqrt(|omega^2_min|) = %.4f" % (bg_amp, amp, sig_lin))
    runs = [(+1, 0.2, DT_BASE), (-1, 0.2, DT_BASE),
            (+1, 0.1, DT_BASE), (+1, 0.2, DT_BASE / 2)]
    esc_all = []
    for sgn, eps, dt in runs:
        dyn = Dyn(h, eps)
        ts, prj, emax, t_esc, reason = run(dyn, g0_sc, mode_sc,
                                           sgn * amp, dt, bg_amp)
        sig = fit_sigma(ts, prj, prj[0] if len(prj) else np.nan)
        tagr = "znak %+d eps=%.1f dt=%.4f" % (sgn, eps, dt)
        print("  [%s] %s; t_esc=%s; gate |dE|/E=%.2e (%s); sigma_fit=%s"
              % (tagr, reason,
                 "%.2f" % t_esc if t_esc is not None else "--",
                 emax, "PASS" if emax <= GATE_E else "FAIL",
                 "%.4f (vs lin %.4f, odchyl %+.1f%%)"
                 % (sig, sig_lin, 100 * (sig - sig_lin) / sig_lin)
                 if np.isfinite(sig) else "brak okna"))
        if sgn == +1 and eps == 0.2 and dt == DT_BASE:
            base_esc = t_esc
        esc_all.append((tagr, t_esc))
    verdict_rows.append((label, r3["w2min"], esc_all))

print("\n" + "=" * 78)
print("PODSUMOWANIE Phase 4 (kryteria LOCK sec. 3 Phase 4, doslownie):")
all_neg = True
all_escape_2t = True
any_no_escape_3t = False
for label, w2, esc_all in verdict_rows:
    if not (w2 < 0):
        all_neg = False
    main = [t for tag, t in esc_all if "znak" in tag]
    esc_main = [t for tag, t in esc_all
                if tag.startswith(("znak +1 eps=0.2 dt=0.0040",
                                   "znak -1 eps=0.2 dt=0.0040"))]
    row_esc = all(t is not None and t <= 2 * T_REF for t in esc_main)
    row_noesc = all(t is None for t in esc_main)
    if not row_esc:
        all_escape_2t = False
    if row_noesc:
        any_no_escape_3t = True
    print("  %s: omega^2_min=%+.6f; ucieczki (znaki +-, eps=0.2): %s"
          % (label, w2,
             ", ".join("%s" % ("%.2f" % t if t is not None else "brak")
                       for t in esc_main)))
print("\n  wszystkie tla omega^2_min<0 zbieznie: %s" % all_neg)
print("  ucieczka w t <= 2t*_ref=7.24 we wszystkich glownych biegach: %s"
      % all_escape_2t)
print("  jakiekolwiek tlo bez ucieczki do 3t*=10.86 przy obu znakach: %s"
      % any_no_escape_3t)
print("  (zlozenie z kryteriami Q -> Phase_FINAL_close.md; odczyt")
print("   kwantyfikatora d wg Phase3_verdict_ruling.md)")
print("=" * 78)
