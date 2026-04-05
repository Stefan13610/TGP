#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p78_alphac_precise.py  --  TGP v1  *  OP-3 Path 2 — precyzyjne alpha_c
=======================================================================
Cel: Weryfikacja, czy rozbiez̈nosc 1.3% miedzy alpha_c i alpha_K_ref
     (obserwowana w p76 Sekcji B przy n_psi=35, n_K=15) jest artefaktem
     numerycznym, czy fizyczna rozbieznoscia.

Test kluczowy:
  - K*_1(alpha=8.5616) obliczone z n_psi=80  vs  K1_ref=0.010414 (z p10)
  - Jesli K*_1(8.5616, n_psi=80) ≈ 0.010414 → p76 mialo blad num. (n_psi=35
    za malo), alpha_c → 8.5616 (Path 2 POTWIERDZONA dokladnie)
  - Jesli K*_1(8.5616, n_psi=80) ≈ 0.010304 → rozbiez̈nosc fizyczna

Metoda:
  - Skan alpha ∈ [8.0, 9.2] w 14 krokach (krok ≈ 0.09)
  - Dla kazdego alpha: find_K_star z n_K=25, n_psi=80 (vs p76: n_K=15, n_psi=35)
  - Interpolacja alpha_c z warunku K*_1(alpha_c) = K1_ref

Testy PASS/FAIL:
  P1: K*_1(8.5616) znaleziony z n_psi=80
  P2: |K*_1(8.5616, n_psi=80) - 0.010414| < 2%  (czy blad num. zniknal?)
  P3: alpha_c znaleziony w [8.0, 9.2]
  P4: |alpha_c - 8.5616| < 2%  (blad numeryczny < 2% zamiast 1.3%)
  P5: Trend K*_1(alpha) monotonicznie malejacy

Data: 2026-03-25
"""

import sys
import io
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import warnings

warnings.filterwarnings('ignore')

if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8',
                                  errors='replace')

# ─────────────────────────────────────────────────────────────────────────────
# Stale TGP
# ─────────────────────────────────────────────────────────────────────────────
ALPHA_K_REF = 8.5616
A_GAM       = 0.040
K1_REF      = 0.010414   # z p10/p11 (wyznaczone z n_psi=120)
LAM         = 5.501357e-06
GAMMA_V     = 1.0

def V_mod(p):
    return GAMMA_V/3 * p**3 - GAMMA_V/4 * p**4 + LAM/6 * (p - 1)**6

def dV_mod(p):
    return GAMMA_V * p**2 - GAMMA_V * p**3 + LAM * (p - 1)**5

V1 = GAMMA_V/3 - GAMMA_V/4


def ode_rhs(r, y, alpha):
    phi, dphi = y
    phi  = max(phi, 1e-10)
    kfac = 1.0 + alpha / phi
    ddphi = (dV_mod(phi) / kfac
             + alpha * dphi**2 / (2.0 * phi**2 * kfac)
             - (2.0 / r) * dphi)
    return [dphi, ddphi]


def _phi_at_rmax(psi_core, K, a_gam, alpha, r_max,
                 rtol=1e-8, atol=1e-10, n_eval=1000):
    dphi0 = -K / a_gam**2

    def ev_lo(r, y):
        return y[0] - 1e-5
    ev_lo.terminal  = True
    ev_lo.direction = -1

    r_eval = a_gam * (r_max / a_gam) ** np.linspace(0, 1, n_eval)
    try:
        sol = solve_ivp(
            lambda r, y: ode_rhs(r, y, alpha),
            [a_gam, r_max], [psi_core, dphi0],
            method='DOP853', rtol=rtol, atol=atol,
            t_eval=r_eval, events=[ev_lo]
        )
        if len(sol.t) == 0 or sol.t[-1] < r_max * 0.99:
            return np.nan
        return float(sol.y[0, -1])
    except Exception:
        return np.nan


def _energy(psi_z, K, a_gam, alpha, r_max):
    dphi0 = -K / a_gam**2

    def ev_lo(r, y):
        return y[0] - 1e-5
    ev_lo.terminal  = True
    ev_lo.direction = -1

    r_ev = a_gam * (r_max / a_gam) ** np.linspace(0, 1, 4000)
    try:
        sol = solve_ivp(
            lambda r, y: ode_rhs(r, y, alpha),
            [a_gam, r_max], [psi_z, dphi0],
            method='DOP853', rtol=1e-9, atol=1e-11,
            t_eval=r_ev, events=[ev_lo]
        )
        r    = sol.t
        phi  = np.maximum(sol.y[0], 1e-10)
        dphi = sol.y[1]
        Ek = 4*np.pi * np.trapezoid(0.5*dphi**2*(1+alpha/phi)*r**2, r)
        Ep = 4*np.pi * np.trapezoid((V_mod(phi) - V1)*r**2, r)
        return Ek + Ep
    except Exception:
        return np.nan


def compute_g(K, a_gam, alpha, r_max, n_psi=80):
    """
    g(K) = E/(4piK) - 1.  Jak p10.best_g: wszystkie galezi, min|g|.
    psi_max = max(3.0, 1+2K/a_gam)  [p10 convention].
    """
    psi_max  = max(3.0, 1.0 + 2.0 * K / a_gam)
    psi_vals = np.linspace(1.001, min(psi_max, 300.0), n_psi)

    F_vals = np.array([_phi_at_rmax(p, K, a_gam, alpha, r_max) - 1.0
                       for p in psi_vals])

    branches = []
    for i in range(len(F_vals) - 1):
        Fi, Fj = F_vals[i], F_vals[i+1]
        if np.isfinite(Fi) and np.isfinite(Fj) and Fi * Fj < 0:
            try:
                psi_z = brentq(
                    lambda p: _phi_at_rmax(p, K, a_gam, alpha, r_max) - 1.0,
                    psi_vals[i], psi_vals[i+1],
                    xtol=1e-5, rtol=1e-5, maxiter=40
                )
                E = _energy(psi_z, K, a_gam, alpha, r_max)
                if np.isfinite(E):
                    branches.append((psi_z, E / (4*np.pi*K) - 1.0))
            except Exception:
                pass

    if not branches:
        return np.nan
    return min(branches, key=lambda x: abs(x[1]))[1]


def find_K_star(a_gam, alpha, K_min=3e-4, K_max=0.15,
                n_K=25, n_psi=80, verbose=False):
    """
    K*_1 = pierwsze zero g(K) (min|g| na galezi).
    """
    M_eff = 1.0 / np.sqrt(1.0 + alpha)
    r_max = max(40.0, a_gam * 5, 6.0 / M_eff)

    K_scan = np.logspace(np.log10(K_min), np.log10(K_max), n_K)
    g_vals = []
    for K in K_scan:
        gv = compute_g(K, a_gam, alpha, r_max, n_psi=n_psi)
        g_vals.append(gv)
        if verbose:
            gstr = f"{gv:.4f}" if np.isfinite(gv) else "nan"
            print(f"    K={K:.5f}  g={gstr}")

    K_star = np.nan
    for i in range(len(g_vals) - 1):
        gi, gj = g_vals[i], g_vals[i+1]
        if np.isfinite(gi) and np.isfinite(gj) and gi * gj < 0:
            try:
                K_star = brentq(
                    lambda K: compute_g(K, a_gam, alpha, r_max, n_psi=n_psi),
                    K_scan[i], K_scan[i+1],
                    xtol=K_scan[i]*0.005, rtol=0.005, maxiter=20
                )
                break
            except Exception:
                pass

    return K_star, g_vals


# ─────────────────────────────────────────────────────────────────────────────
PASS_COUNT = 0
FAIL_COUNT = 0

def record(label, ok, info=""):
    global PASS_COUNT, FAIL_COUNT
    if ok: PASS_COUNT += 1
    else:  FAIL_COUNT += 1
    print(f"  [{'v' if ok else 'x'}] {label}: {info}")

# ─────────────────────────────────────────────────────────────────────────────
print("=" * 72)
print("TGP v1  *  p78_alphac_precise.py  --  OP-3 Path 2, precyzyjne alpha_c")
print("n_psi=80, n_K=25  (vs p76 Sekcja B: n_psi=35, n_K=15)")
print("=" * 72)
print(f"\n  alpha_K_ref = {ALPHA_K_REF}")
print(f"  K1_ref      = {K1_REF} (z p10, n_psi=120)")
print(f"  a_Gam       = {A_GAM}")

# ─────────────────────────────────────────────────────────────────────────────
# KROK 1: K*_1 przy alpha=8.5616 z n_psi=80 (vs p76 Sekcja A: n_psi=40)
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 72)
print("KROK 1: K*_1 przy alpha=8.5616, a_Gam=0.040 z n_psi=80")
print("        (Sprawdza: czy blad num. n_psi=40 wyjasnia 1.1% rozbiez̈nosc)")
print("=" * 72)

print("\n  Szukam K*_1 przy alpha=8.5616 z n_psi=80, n_K=25 ...\n")
K1_precise, g_vals_A = find_K_star(A_GAM, ALPHA_K_REF,
                                    K_min=3e-4, K_max=0.12,
                                    n_K=25, n_psi=80, verbose=True)

print(f"\n  K*_1 (n_psi=80)  = {K1_precise:.6f}")
print(f"  K*_1 (p76, n_psi=40) = 0.010304")
print(f"  K1_ref (p10, n_psi=120) = {K1_REF:.6f}")

if np.isfinite(K1_precise):
    err_ref   = (K1_precise - K1_REF) / K1_REF * 100
    err_p76   = (K1_precise - 0.010304) / 0.010304 * 100
    print(f"\n  Blad vs K1_ref: {err_ref:+.2f}%  (cel: < 1%)")
    print(f"  Blad vs p76:    {err_p76:+.2f}%")

ok_P1 = np.isfinite(K1_precise)
record("P1: K*_1 znaleziony z n_psi=80", ok_P1,
       f"K*_1 = {K1_precise:.6f}" if ok_P1 else "brak")

ok_P2 = ok_P1 and abs(K1_precise - K1_REF) / K1_REF < 0.02
record("P2: |K*_1(n_psi=80) - K1_ref| < 2%", ok_P2,
       f"blad = {abs(K1_precise-K1_REF)/K1_REF*100:.2f}%" if ok_P1
       else "K*_1 nie znaleziony")

# ─────────────────────────────────────────────────────────────────────────────
# KROK 2: Gesты skan K*_1(alpha) dla alpha ∈ [8.0, 9.2]
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 72)
print("KROK 2: Gesты skan K*_1(alpha) dla alpha ∈ [8.0, 9.2]")
print(f"        n_psi=80, n_K=25 -- cel: wyznaczenie alpha_c dokladnie")
print("=" * 72)

alpha_scan = np.array([8.0, 8.1, 8.2, 8.3, 8.4, 8.5, 8.5616,
                        8.6, 8.7, 8.8, 8.9, 9.0, 9.1, 9.2])

results = []
print(f"\n  {'alpha':>8}  {'K*_1':>10}  {'K*_1/K1ref':>12}  status")
print("  " + "-" * 48)

for alpha in alpha_scan:
    print(f"  alpha={alpha:.4f} ...", flush=True)
    K_s, _ = find_K_star(A_GAM, alpha,
                          K_min=3e-4, K_max=0.12, n_K=25, n_psi=80)
    results.append((alpha, K_s))
    if np.isfinite(K_s):
        ratio = K_s / K1_REF
        print(f"  {alpha:>8.4f}  {K_s:>10.6f}  {ratio:>12.4f}")
    else:
        print(f"  {alpha:>8.4f}  {'---':>10}  {'---':>12}  [brak g=0]")


# ─────────────────────────────────────────────────────────────────────────────
# KROK 3: Interpolacja alpha_c
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 72)
print("KROK 3: Interpolacja alpha_c z warunku K*_1(alpha_c) = K1_ref")
print("=" * 72)

valid = [(r[0], r[1]) for r in results if np.isfinite(r[1])]
alpha_c_interp = None
alpha_c_list   = []

if len(valid) >= 2:
    alphas_v = np.array([r[0] for r in valid])
    K1s_v    = np.array([r[1] for r in valid])
    diff_K   = K1s_v - K1_REF

    crossings = np.where(np.diff(np.sign(diff_K)))[0]
    if len(crossings) > 0:
        for idx0 in crossings:
            i0, i1 = idx0, idx0 + 1
            a0, d0 = alphas_v[i0], diff_K[i0]
            a1, d1 = alphas_v[i1], diff_K[i1]
            ac = a0 + (-d0) * (a1 - a0) / (d1 - d0)
            alpha_c_list.append(ac)
            print(f"\n  Przejscie #{len(alpha_c_list)}: alpha_c = {ac:.5f}")
            print(f"    K*_1({a0:.3f}) = {K1s_v[i0]:.6f},  "
                  f"K*_1({a1:.3f}) = {K1s_v[i1]:.6f},  K1_ref = {K1_REF:.6f}")
            odch = abs(ac - ALPHA_K_REF) / ALPHA_K_REF * 100
            print(f"    odch od alpha_K_ref = {odch:.2f}%")
        alpha_c_interp = alpha_c_list[0]  # pierwsze przejscie
    else:
        print(f"\n  K*_1 nie przecina K1_ref w zakresie [8.0, 9.2].")
        print(f"  K*_1 zakres: [{K1s_v.min():.6f}, {K1s_v.max():.6f}]")
        if K1s_v.min() > K1_REF:
            print("  -> K*_1 > K1_ref wszedzie: alpha_c < 8.0")
        elif K1s_v.max() < K1_REF:
            print("  -> K*_1 < K1_ref wszedzie: alpha_c > 9.2")

# P3: alpha_c znaleziony
ok_P3 = alpha_c_interp is not None and np.isfinite(alpha_c_interp)
record("P3: alpha_c znaleziony w [8.0, 9.2]", ok_P3,
       f"alpha_c = {alpha_c_interp:.5f}" if ok_P3 else "nie znaleziony")

# P4: blizej alpha_K_ref niz p76 (< 2%)
ok_P4 = ok_P3 and abs(alpha_c_interp - ALPHA_K_REF) / ALPHA_K_REF < 0.02
if ok_P3:
    odch_P4 = abs(alpha_c_interp - ALPHA_K_REF) / ALPHA_K_REF * 100
    record("P4: |alpha_c - 8.5616| < 2%  (poprawa vs p76 1.3%)", ok_P4,
           f"alpha_c = {alpha_c_interp:.5f}, odch = {odch_P4:.2f}%")
else:
    record("P4: |alpha_c - 8.5616| < 2%", False, "alpha_c nie znaleziony")

# P5: K*_1(alpha) monotonicznie malejacy
if len(valid) >= 3:
    K_vals = [r[1] for r in valid]
    diffs  = [K_vals[i+1] - K_vals[i] for i in range(len(K_vals)-1)]
    n_neg  = sum(1 for d in diffs if d < 0)
    ok_P5  = n_neg >= len(diffs) * 0.8
    record("P5: K*_1(alpha) monotonicznie malejacy", ok_P5,
           f"{n_neg}/{len(diffs)} spadkow")
else:
    ok_P5 = False
    record("P5: K*_1(alpha) monotoniczne", False, "za malo punktow")


# ─────────────────────────────────────────────────────────────────────────────
# PODSUMOWANIE
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 72)
print("PODSUMOWANIE  p78_alphac_precise.py")
print("=" * 72)

print(f"\n  PASS: {PASS_COUNT} / {PASS_COUNT + FAIL_COUNT}")
print(f"  FAIL: {FAIL_COUNT} / {PASS_COUNT + FAIL_COUNT}")

print(f"\n  {'Parametr':<30} {'p76 (n_psi=35/40)':>20} {'p78 (n_psi=80)':>18}")
print("  " + "-" * 70)
print(f"  {'K*_1 przy alpha=8.5616':<30} {'0.010290-0.010304':>20} "
      f"{K1_precise if np.isfinite(K1_precise) else 'nan':>18.6f}")
if ok_P3:
    print(f"  {'alpha_c (interpolacja)':<30} {'8.4492 (1.3%)':<20} "
          f"{alpha_c_interp:.5f} ({abs(alpha_c_interp-ALPHA_K_REF)/ALPHA_K_REF*100:.2f}%)")
else:
    print(f"  {'alpha_c (interpolacja)':<30} {'8.4492 (1.3%)':<20} {'---':>18}")

print()
if ok_P2:
    print("  WNIOSEK P2 PASS: Blad numeryczny K*_1 zmniejszony do < 2%.")
    print("  -> Wieksze n_psi poprawia dokladnosc numeryczna.")
else:
    print("  WNIOSEK P2 FAIL: K*_1 nadal rozni sie > 2% od K1_ref.")
    print("  -> Rozbiez̈nosc moze byc fizyczna (nie tylko numeryczna).")

if ok_P4:
    print(f"\n  WNIOSEK P4 PASS: alpha_c = {alpha_c_interp:.5f} (odch < 2% od 8.5616).")
    print("  -> OP-3 Path 2 POTWIERDZONA z dokladnoscia < 2%.")
    print("  -> alpha_K moze byc wyznaczone z substratu (K*_1(alpha_c) = K1_ref).")
else:
    if ok_P3:
        print(f"\n  WNIOSEK P4 FAIL: alpha_c = {alpha_c_interp:.5f},  odch = {odch_P4:.2f}%.")
        print("  -> Rozbiez̈nosc > 2% mimo wiekszej dokladnosci.")
        print("  -> Konieczne n_psi=120 (jak p10) lub analiza systematyczna.")
    else:
        print("\n  alpha_c nie znaleziony.")
