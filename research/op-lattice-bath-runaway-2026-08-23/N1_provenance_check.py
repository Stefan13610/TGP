#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
N1_provenance_check.py - weryfikacja pochodzenia faz z dodatekH lin. 1126-1129
POST-CLOSE (cykl CLOSED-GATE-STOP; sledztwo NEEDS N1, nie kontynuacja faz).

HIPOTEZA N1: fazy delta_e=-81.4, delta_mu=+38.6, Delta=120.01 pochodza ze
skryptu _archiwum/scripts_exploratory/advanced/p131_eta_refinement.py, ktory
uzywa INNEGO ukladu niz dokumentuje dodatekH/status_map O-L5:
  (P131)  F(g)*g'' + (2/r)*g' = V'(g),  F(g)=1+2*alpha_eff(g)*ln(g),
          alpha_eff(g)=2/(1+eta_K*(g-1)^2),  V'(g)=g^2*(1-g)
zamiast audytowanego (M2): g'' + (2/r)g' + (alpha/g)g'^2 = g^2*(1-g).

TEST: dokladna re-implementacja P131 (rownanie, okno [120,260], konwencja
atan2(C,B) na (g-1)*r ~ B cos r + C sin r, rmax=300, te same g0) i porownanie
z dodatekH: delta_e=-81.4, delta_mu=+38.6, Delta=120.01, A_mu(1.4550)=0.3861.
eta_K = 181/15 (analityczne p139-p140; p131 bisekowal [12,18] do r31=3477).
"""
import io, sys
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from scipy.integrate import solve_ivp

PHI = (1 + np.sqrt(5)) / 2

def solve_p131(g0, eta_K, rmax=300):
    def fk(g):
        a = 2.0 / (1 + eta_K * (g - 1) ** 2)
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
    r = np.linspace(rs, min(s.t[-1], rmax), 20000)
    return r, s.sol(r)[0]

def extract_tail_p131(r, g):
    m = (r >= 120) & (r <= 260)
    rf, tl = r[m], (g[m] - 1) * r[m]
    if len(rf) < 10:
        return np.nan, np.nan, np.nan, np.nan
    A = np.column_stack([np.cos(rf), np.sin(rf)])
    coeff, _, _, _ = np.linalg.lstsq(A, tl, rcond=None)
    B, C = coeff
    return B, C, np.sqrt(B ** 2 + C ** 2), np.arctan2(C, B)

def report(tag, g0, eta):
    r, g = solve_p131(g0, eta)
    if r is None or r[-1] < 260:
        print("  %s: g0=%.5f  -> NIE DOBIEGL (r_end=%.1f)" % (tag, g0, 0 if r is None else r[-1]))
        return None
    B, C, A, d = extract_tail_p131(r, g)
    print("  %s: g0=%.5f  A=%.6f  delta=%+8.2f deg" % (tag, g0, A, np.degrees(d)))
    return d, A

ETA = 181.0 / 15.0
print("=" * 72)
print("[N1] Re-implementacja ukladu P131 (F=1+2*a_eff*ln g; tarcie poza F)")
print("     eta_K = 181/15 = %.6f" % ETA)
print("=" * 72)

for g0e_tag, g0e in [("baseline p131", 0.89926559), ("dodatekH", 0.90548)]:
    print("\n--- g0_e = %s (%s); g0_mu = phi*g0_e = %.5f ---" % (g0e, g0e_tag, PHI * g0e))
    res_e = report("e ", g0e, ETA)
    res_mu = report("mu", PHI * g0e, ETA)
    if res_e and res_mu:
        d_em = np.degrees(res_mu[0] - res_e[0]) % 360
        print("  Delta(e->mu) = %.2f deg   [dodatekH: 120.01]" % d_em)

print("\n--- mapa fazowa p127-p128: A_mu przy g0=1.4550 [dodatekH: 0.3861] ---")
report("mu", 1.4550, ETA)

print("\n--- czulosc na eta_K ---")
for eta in [12.0, ETA, 12.5, 13.0, 14.0]:
    print("  eta_K = %.4f:" % eta)
    res_a = report("   e ", 0.89926559, eta)
    res_b = report("   mu", PHI * 0.89926559, eta)
    if res_a and res_b:
        print("     Delta(e->mu) = %.2f deg" % (np.degrees(res_b[0] - res_a[0]) % 360))

print("\n--- kontrola krzyzowa: uklad M2 (audytowany) na tych samych g0 ---")
def solve_m2(g0, alpha=2.0, rmax=300):
    def rhs(r, y):
        g, p = y
        if g <= 1e-12:
            return [p, 0]
        s = g * g * (1 - g) - (alpha / g) * p * p
        if r < 1e-10:
            return [p, s / 3]
        return [p, s - 2 / r * p]
    def ev(r, y):
        return 100 - abs(y[0])
    ev.terminal = True
    s = solve_ivp(rhs, [0.01, rmax], [g0, 0.0], method='RK45',
                  rtol=1e-11, atol=1e-13, max_step=0.05,
                  events=[ev], dense_output=True)
    r = np.linspace(0.01, min(s.t[-1], rmax), 20000)
    return r, s.sol(r)[0]

for tag, g0 in [("e ", 0.90548), ("mu", PHI * 0.90548)]:
    r, g = solve_m2(g0)
    if r[-1] >= 260:
        B, C, A, d = extract_tail_p131(r, g)
        print("  M2 %s: g0=%.5f  A=%.6f  delta=%+8.2f deg" % (tag, g0, A, np.degrees(d)))
    else:
        print("  M2 %s: g0=%.5f  nie dobiegl (r_end=%.1f)" % (tag, g0, r[-1]))

print("\n[WNIOSEK] patrz ANALIZA_N1_pochodzenie-faz_2026-08-23.md")
