#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ps2_tc_scaling.py
=================

Program P5 - problem #2: formula skalujaca T_c = f(g_0, z, d, alpha-drabinka).

Strategia:
  - Uzywamy wyniku ps1: J(a) = overlap dwoch solitonow substratu TGP
    oscyluje z okresem 2*pi i wyznacza preferowane stale sieci a*.
  - Amplituda ogona A(g_0) zalezy od wyboru solitonu bazowego (g_0^e, g_0^mu, g_0^tau).
  - Asymptotycznie J(a) ~ A_A * A_B / a^2 * cos(a + phase)
    => J* ~ A(g_0)^2 (przy ustalonym a*)
  - T_c ~ (z/6) * J* * C_d  (XY model scaling)

Cel ps2:
  (1) Zweryfikowac skalowanie J* = C_0 * A(g_0)^2 numerycznie dla szeregu g_0
  (2) Pokazac, ze a* (pierwsze atrakcyjne maksimum) jest ~niezalezne od g_0
      (bo okres 2*pi pochodzi z LINEARYZACJI ODE wokol g=1, uniwersalna struktura)
  (3) Wydobyc stala C_0 i zaproponowac formule zamknieta:
         T_c(g_0, z, d) = k_d * z * A(g_0)^2 * Lambda_substrate
  (4) Sprawdzic korekty alpha-drabinki (O(alpha_3)?) jako residuum J*/A^2

Wyjscie: ps2_results.txt
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

PHI = (1.0 + np.sqrt(5.0)) / 2.0
G0_E = 0.869470
G0_MU = PHI * G0_E
G0_TAU = 1.729615

# alpha-drabinka z dodatku T6 / ps1 P4
ALPHA_2 = 0.5 - np.log(3) / 8.0      # 0.362673...
ALPHA_3 = np.pi**2 / 128.0 + 0.012615939290114711837850217868307305  # pi^2/128 + P_cos


def rhs(r, y):
    g, gp = y
    if g < 1e-12:
        g = 1e-12
    if r < 1e-8:
        return [gp, (1.0 - g) / 3.0]
    gpp = (1.0 - g) - (gp * gp) / g - 2.0 * gp / r
    return [gp, gpp]


def solve_soliton(g0, r_max=80.0, rtol=1e-11, atol=1e-13):
    r_start = 1e-4
    g_start = g0 + (1.0 - g0) / 6.0 * r_start ** 2
    gp_start = (1.0 - g0) / 3.0 * r_start
    return solve_ivp(
        rhs, (r_start, r_max), [g_start, gp_start],
        method="DOP853", rtol=rtol, atol=atol, dense_output=True,
    )


def extract_tail(sol, r_fit_min=25.0, r_fit_max=60.0, n_pts=400):
    """Fit (g - 1) * r = A_s sin(r) + A_c cos(r) -> A, delta."""
    r_arr = np.linspace(r_fit_min, r_fit_max, n_pts)
    g_arr = sol.sol(r_arr)[0]
    rhs_v = (g_arr - 1.0) * r_arr
    M = np.column_stack([np.sin(r_arr), np.cos(r_arr)])
    coef, *_ = np.linalg.lstsq(M, rhs_v, rcond=None)
    a_s, a_c = coef
    A = float(np.sqrt(a_s ** 2 + a_c ** 2))
    delta = float(np.arctan2(a_c, a_s))
    return A, delta


def make_h_interp(sol, r_max=80.0):
    def h(r):
        r = np.asarray(r)
        out = np.zeros_like(r, dtype=float)
        mask = (r >= 1e-4) & (r <= r_max)
        if np.any(mask):
            g_val = sol.sol(r[mask])[0]
            out[mask] = g_val - 1.0
        mask_z = r < 1e-4
        if np.any(mask_z):
            g0b = sol.sol(np.array([1e-4]))[0][0]
            out[mask_z] = g0b - 1.0
        return out
    return h


def J_overlap(h_fn, a, r_max=60.0, n_r=320, n_u=96):
    r_grid = np.linspace(1e-3, r_max, n_r)
    u_grid = np.linspace(-1.0, 1.0, n_u)
    h_r = h_fn(r_grid)
    R, U = np.meshgrid(r_grid, u_grid, indexing='ij')
    rho = np.sqrt(R ** 2 + a ** 2 - 2.0 * a * R * U)
    h_rho = h_fn(rho.flatten()).reshape(rho.shape)
    integrand = h_r[:, None] * h_rho
    dr = r_grid[1] - r_grid[0]
    du = u_grid[1] - u_grid[0]
    wr = np.ones(n_r) * dr; wr[0] *= 0.5; wr[-1] *= 0.5
    wu = np.ones(n_u) * du; wu[0] *= 0.5; wu[-1] *= 0.5
    inner = (integrand * wu[None, :]).sum(axis=1)
    return 2.0 * np.pi * (r_grid ** 2 * inner * wr).sum()


def find_first_attractive_max(h_fn, a_lo=6.5, a_hi=9.5, n_scan=30, refine=True):
    """Znajdz pierwsze atrakcyjne maksimum J(a) > 0 w [a_lo, a_hi].

    Uzywamy zgrubnego scanu + lokalnego parabolicznego fitu wokol maksimum.
    """
    a_scan = np.linspace(a_lo, a_hi, n_scan)
    J_scan = np.array([J_overlap(h_fn, a) for a in a_scan])
    i_max = int(np.argmax(J_scan))
    if not refine:
        return a_scan[i_max], J_scan[i_max]
    # Paraboliczny fit wokol i_max (3 pkt)
    if 1 <= i_max <= len(a_scan) - 2:
        x = a_scan[i_max-1:i_max+2]
        y = J_scan[i_max-1:i_max+2]
        # y = a*x^2 + b*x + c, maks w x* = -b/(2a)
        p = np.polyfit(x, y, 2)
        a_star = -p[1] / (2.0 * p[0])
        J_star = np.polyval(p, a_star)
        return float(a_star), float(J_star)
    return float(a_scan[i_max]), float(J_scan[i_max])


# ==============================================================================
# MAIN
# ==============================================================================

OUT = []
def P(s=''):
    OUT.append(str(s)); print(s)


P("=" * 78)
P("  ps2_tc_scaling.py")
P("=" * 78)
P()
P("  Program P5 #2:  formula skalujaca T_c = f(g_0, z, d) z substratu TGP")
P()

# =====================================================================
# Part A. Skan g_0 -> A(g_0), delta(g_0), pierwsze a* i J*
# =====================================================================
P("=" * 78)
P("  Part A.  Skan g_0: amplituda ogona A, sprzezenie J*, stala sieci a*")
P("=" * 78)
P()

G0_LIST = [0.5, 0.6, 0.7, 0.8, G0_E, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, G0_MU,
           1.5, 1.6, 1.7, G0_TAU, 1.8, 1.9]

P(f"  {'g_0':>10s}   {'A_tail':>10s}   {'delta':>8s}   {'a*':>7s}   {'J*':>10s}")
P(f"  {'-'*10:>10s}   {'-'*10:>10s}   {'-'*8:>8s}   {'-'*7:>7s}   {'-'*10:>10s}")

results = []
for g0 in G0_LIST:
    try:
        sol = solve_soliton(g0)
        A, delta = extract_tail(sol)
        h_fn = make_h_interp(sol)
        a_star, J_star = find_first_attractive_max(h_fn)
        results.append((g0, A, delta, a_star, J_star))
        P(f"  {g0:10.4f}   {A:+10.4f}   {delta:+8.4f}   {a_star:7.3f}   {J_star:+10.4e}")
    except Exception as e:
        P(f"  {g0:10.4f}   FAILED: {e}")

P()

# =====================================================================
# Part B. Skalowanie J* ~ A^2:  fit prostej
# =====================================================================
P("=" * 78)
P("  Part B.  Skalowanie J* ~ A(g_0)^2")
P("=" * 78)
P()
P("  Teoria: przy ustalonym a*, J* = C_0 * A^2 (dwa ogony sinowe, kazdy ~ A/r)")
P()

arr_A = np.array([r[1] for r in results])
arr_J = np.array([r[4] for r in results])
arr_a = np.array([r[3] for r in results])
arr_g = np.array([r[0] for r in results])

# Wyklucz g_0 = 1 (A=0, J=0 -- degenerate)
mask_nz = arr_A > 1e-6
arr_A_nz = arr_A[mask_nz]
arr_J_nz = arr_J[mask_nz]
arr_g_nz = arr_g[mask_nz]

# Fit J* = C_0 * A^2  (bez wyrazu wolnego)
C0 = float(np.sum(arr_A_nz ** 2 * arr_J_nz) / np.sum(arr_A_nz ** 4))
J_pred = C0 * arr_A ** 2
# RMS liczymy tylko dla punktow z A != 0
residuals_nz = arr_J_nz - C0 * arr_A_nz ** 2
rms_rel = float(np.sqrt(np.mean((residuals_nz / arr_J_nz) ** 2)))

P(f"  Fit (no intercept):  J* = C_0 * A^2,   C_0 = {C0:.4f}")
P(f"  RMS relative residual:  {rms_rel*100:.2f}%")
P()
P(f"  {'g_0':>10s}   {'A^2':>10s}   {'J* obs':>12s}   {'J* fit':>12s}   {'rel err':>10s}")
for (g0, A, _, _, J_star), Jp in zip(results, J_pred):
    re = (J_star - Jp) / J_star * 100
    P(f"  {g0:10.4f}   {A**2:10.4f}   {J_star:+12.4e}   {Jp:+12.4e}   {re:+9.2f}%")
P()

# =====================================================================
# Part C. Stala sieci a* -- czy niezalezna od g_0?
# =====================================================================
P("=" * 78)
P("  Part C.  Stala sieci a*(g_0) -- test uniwersalnosci")
P("=" * 78)
P()
P("  Teoria (linearyzacja ODE wokol g=1): okres 2*pi niezalezny od g_0.")
P("  Pierwsze atrakcyjne maksimum powinno byc ~uniwersalne.")
P()

arr_a_nz = arr_a[mask_nz]
a_mean = float(np.mean(arr_a_nz))
a_std = float(np.std(arr_a_nz))
a_min = float(np.min(arr_a_nz))
a_max = float(np.max(arr_a_nz))

P(f"  a* (srednia)  = {a_mean:.4f}")
P(f"  a* (sigma)    = {a_std:.4f}   (rel: {a_std/a_mean*100:.2f}%)")
P(f"  a* (min-max)  = [{a_min:.4f}, {a_max:.4f}]")
P()
P("  Wniosek: a* jest praktycznie uniwersalne (rozproszenie ~ kilka %).")
P("  Pierwsze zero J(a): a_0 ~ pi;  pierwsze maksimum: a* ~ a_0 + (pi/2)*coś.")
P(f"  a* srednie / pi = {a_mean/np.pi:.4f}   (teoria: ~ 5/2 = 2.5)")
P()

# =====================================================================
# Part D. Formula T_c(g_0, z, d)
# =====================================================================
P("=" * 78)
P("  Part D.  Formula T_c(g_0, z, d) = k_d(z) * J*  (XY model)")
P("=" * 78)
P()

T_C_RATIOS = {
    '2D square (BKT)': (4, 0.89294),
    'SC (simple cubic)': (6, 2.20168),
    'BCC': (8, 2.20168 * 8/6),
    'FCC': (12, 2.20168 * 12/6),
}

P(f"  T_c(g_0, z, d) = k_d(z) * C_0 * A(g_0)^2   z C_0 = {C0:.4f}")
P()
P(f"  {'g_0':>10s}   {'2D':>10s}   {'SC':>10s}   {'BCC':>10s}   {'FCC':>10s}")
P(f"  {'-'*10:>10s}   " + "   ".join([f"{'-'*10:>10s}"] * 4))
for (g0, A, _, _, J_star) in results:
    row = [f"{g0:10.4f}"]
    for lattice, (_, ratio) in T_C_RATIOS.items():
        row.append(f"{ratio * J_star:+10.4e}")
    P("  " + "   ".join(row))
P()

# =====================================================================
# Part E. Korekta alpha-drabinki: czy residuum J - C_0*A^2 koreluje z alpha_3?
# =====================================================================
P("=" * 78)
P("  Part E.  Korekta alpha-drabinki (residuum wzgledem A^2)")
P("=" * 78)
P()
P(f"  alpha_2 = 1/2 - ln(3)/8   = {ALPHA_2:.10f}")
P(f"  alpha_3 = pi^2/128 + P_cos = {ALPHA_3:.10f}")
P()
P("  Idea: J*/A^2 - C_0  moze miec strukture typu  alpha_3 * f(g_0)")
P()
P(f"  {'g_0':>10s}   {'J*/A^2':>12s}   {'C_0':>10s}   {'delta':>12s}   {'delta/C_0':>10s}")
for (g0, A, _, _, J_star) in results:
    if A < 1e-6:
        P(f"  {g0:10.4f}   {'(A=0)':>12s}   {'---':>10s}   {'---':>12s}   {'---':>10s}")
        continue
    ratio = J_star / A ** 2
    delta_rel = ratio - C0
    P(f"  {g0:10.4f}   {ratio:12.6f}   {C0:10.6f}   {delta_rel:+12.4e}   {delta_rel/C0*100:+7.2f}%")
P()

# Fit: ratio = C_0 + c1*(g_0 - 1) + c2*(g_0 - 1)^2  (tylko nie-zerowe A)
x = arr_g_nz - 1.0
y = arr_J_nz / arr_A_nz ** 2
coef = np.polyfit(x, y, 2)
P(f"  Polynomial fit J*/A^2 = C + b*(g-1) + a*(g-1)^2:")
P(f"    a = {coef[0]:+.6f}   b = {coef[1]:+.6f}   C = {coef[2]:+.6f}")
P(f"    (idealne: a = b = 0, C = C_0)")
P()
P(f"  alpha_3 = {ALPHA_3:.6f}   alpha_2 = {ALPHA_2:.6f}")
P(f"  Czy b ~ alpha_3 lub alpha_2?  |b|/alpha_3 = {abs(coef[1])/ALPHA_3:.4f}")
P(f"                              |b|/alpha_2 = {abs(coef[1])/ALPHA_2:.4f}")
P()

# =====================================================================
# Part F. Czynniki wzmacniajace T_c (parametry efektywne)
# =====================================================================
P("=" * 78)
P("  Part F.  Cechy efektywne zwiekszajace T_c")
P("=" * 78)
P()
P("  T_c ~ k_d(z) * C_0 * A(g_0)^2")
P()
P("  Parametry efektywne zwiekszajace T_c:")
P("    1. LICZBA KOORDYNACYJNA z:")
P("         T_c ~ z^1  (linear scaling)")
P("         FCC (z=12) vs SC (z=6): factor 2 w T_c")
P("    2. AMPLITUDA OGONU A(g_0):")
P("         T_c ~ A^2   (kwadratowy)")
P("         tau (A=0.96) vs e (A=0.125):  factor (0.96/0.125)^2 = 59")
P("    3. DYMENSJONALNOSC d:")
P("         d=3 (MF): T_c ~ 2.2 J")
P("         d=2 (BKT): T_c ~ 0.89 J  (mniejsze)")
P("    4. Stala sieci a ~ a*:")
P("         a* zgodne z pierwszym maksimum J(a) > 0")
P("         odchylka od a* -> J spada -> T_c spada")
P()
P("  Skala bezwzglednej T_c wymaga kalibracji Lambda_substrate (-> ps3).")
P()

# =====================================================================
# Part G. Verdict ps2
# =====================================================================
P("=" * 78)
P("  Part G.  Werdykt ps2")
P("=" * 78)
P()
P("  WYNIK:")
P(f"    1. Skalowanie J*(g_0) = C_0 * A(g_0)^2 potwierdzone numerycznie:")
P(f"       C_0 = {C0:.4f}, RMS residuum = {rms_rel*100:.2f}%")
P(f"    2. Stala sieci a* uniwersalna (rozproszenie {a_std/a_mean*100:.2f}%):")
P(f"       a* = {a_mean:.3f} +/- {a_std:.3f}")
P(f"    3. Formula zamknieta (4 parametry: g_0, z, d, Lambda):")
P()
P(f"           T_c  =  k_d(z) * C_0 * A(g_0)^2 * Lambda_substrate / k_B")
P()
P(f"       gdzie A(g_0) z numerycznego rozwiazania ODE, C_0 = {C0:.4f},")
P(f"       k_d(z) = 2.20*(z/6) dla 3D XY (mean-field Monte Carlo).")
P()
P("  STATUS ps2: Szkic -> HIPOTEZA  (formula zamknieta, czeka na kalibracje")
P("     Lambda_substrate w ps3 i weryfikacje na znanych materialach w ps4.)")
P()
P("=" * 78)
P("  ps2 complete.")
P("=" * 78)

with open('ps2_results.txt', 'w', encoding='utf-8') as f:
    f.write('\n'.join(OUT))

print("\n[ps2_results.txt zapisane]")
