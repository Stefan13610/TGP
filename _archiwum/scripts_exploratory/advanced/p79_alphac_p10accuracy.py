#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p79_alphac_p10accuracy.py  --  TGP v1  *  OP-3 Path 2 — alpha_c z dokladnoscia p10
====================================================================================
Cel: Wyznaczenie alpha_c z warunkiem K*_1(alpha_c) = K1_ref,
     uzywajac DOKLADNIE tych samych parametrow numerycznych co p10:
       n_eval=2000 (phi_at_rmax)
       rtol=1e-9, atol=1e-11 (ODE DOP853)
       n_psi=120 (skan psi_core)
       r_max = max(40.0, a_gam*5, 6.0/M_eff)

Diagnoza p78 FAIL:
  p78 uzywalo n_eval=1000 i rtol=1e-8 (vs p10: 2000 i 1e-9)
  -> systematyczny bias +1%: K*_1(p78)=0.010523 vs K1_ref=0.010414
  -> K*_1 > K1_ref wszedzie w [8.0, 9.2] => alpha_c nie znaleziony

Hipoteza p79:
  Z n_eval=2000 i rtol=1e-9 (jak p10), K*_1(alpha=8.5616) ~ K1_ref=0.010414
  => alpha_c znajdzie sie w okolicach 8.5616 (Path 2 potwierdzona)

Testy PASS/FAIL:
  P1: K*_1(alpha=8.5616) znaleziony z dokladnoscia p10
  P2: |K*_1(8.5616, p10-acc) - K1_ref| < 0.5%  (cel: lepiej niz p78's 1%)
  P3: alpha_c znaleziony w [7.5, 9.5]
  P4: |alpha_c - 8.5616| < 1%  (cel: < 1% od alpha_K)
  P5: K*_1(alpha) zachowuje sie monotoniczne (przynajmniej 70% krokow malejacych)

Notatka o nowej fizyce (odkrytej w p78):
  K*_1(alpha) ma MINIMUM przy alpha ~ 8.65 (wg p78)
  Dla alpha >= 8.9 soliton znika (wg p78 z n_eval=1000)
  p79 sprawdzi, czy to fizyczne czy numeryczne

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
K1_REF      = 0.010414   # z p10/p11 (wyznaczone z n_psi=120, n_eval=2000, rtol=1e-9)
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


# ─────────────────────────────────────────────────────────────────────────────
# Kluczowa zmiana: n_eval=2000, rtol=1e-9, atol=1e-11  (jak p10)
# ─────────────────────────────────────────────────────────────────────────────

def _phi_at_rmax(psi_core, K, a_gam, alpha, r_max,
                 rtol=1e-9, atol=1e-11, n_eval=2000):
    """Identyczne parametry numeryczne jak p10.phi_at_rmax."""
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
    """Identyczne parametry jak p10: rtol=1e-9, atol=1e-11, 4000 punktow."""
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


def compute_g(K, a_gam, alpha, r_max, n_psi=120):
    """
    g(K) = E/(4piK) - 1.  Portowane 1:1 z p10.best_g:
      - psi_max = max(3.0, 1+2K/a_gam)
      - n_psi=120 (domyslnie)
      - wszystkie galezi, min|g|
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
                n_K=25, n_psi=120, verbose=False):
    """
    K*_1 = pierwsze zero g(K).
    r_max: identyczna formula jak p10.best_g.
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
                    xtol=K_scan[i]*0.003, rtol=0.003, maxiter=25
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
print("TGP v1  *  p79_alphac_p10accuracy.py  --  OP-3 Path 2, dokladnosc p10")
print("n_psi=120, n_eval=2000, rtol=1e-9  (identyczne z p10)")
print("=" * 72)
print(f"\n  alpha_K_ref = {ALPHA_K_REF}")
print(f"  K1_ref      = {K1_REF} (z p10, n_eval=2000, rtol=1e-9)")
print(f"  a_Gam       = {A_GAM}")
print(f"\n  Diagnoza p78: n_eval=1000/rtol=1e-8 dawalo K*_1=0.010523 (+1.04%)")
print(f"  Cel p79: K*_1(alpha=8.5616) ~ 0.010414 (< 0.5% bledu)")

# ─────────────────────────────────────────────────────────────────────────────
# KROK 1: K*_1 przy alpha=8.5616 z dokladnoscia p10
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 72)
print("KROK 1: K*_1 przy alpha=8.5616 z parametrami p10 (n_eval=2000, rtol=1e-9)")
print("=" * 72)

print("\n  Szukam K*_1 przy alpha=8.5616, n_psi=120, n_K=25 ...\n")
K1_precise, g_vals_A = find_K_star(A_GAM, ALPHA_K_REF,
                                    K_min=3e-4, K_max=0.15,
                                    n_K=25, n_psi=120, verbose=True)

print(f"\n  K*_1 (p79, p10-acc)   = {K1_precise:.6f}" if np.isfinite(K1_precise)
      else "\n  K*_1 (p79, p10-acc)   = nie znaleziony")
print(f"  K*_1 (p78, n_eval=1000)= 0.010523  (+1.04%)")
print(f"  K*_1 (p76, n_eval=800) = 0.010304  (-1.06%)")
print(f"  K1_ref (p10)           = {K1_REF:.6f}")

if np.isfinite(K1_precise):
    err_ref = (K1_precise - K1_REF) / K1_REF * 100
    print(f"\n  Blad vs K1_ref: {err_ref:+.2f}%  (cel: < 0.5%)")

ok_P1 = np.isfinite(K1_precise)
record("P1: K*_1 znaleziony z dokladnoscia p10", ok_P1,
       f"K*_1 = {K1_precise:.6f}" if ok_P1 else "brak")

ok_P2 = ok_P1 and abs(K1_precise - K1_REF) / K1_REF < 0.005
record("P2: |K*_1(p10-acc) - K1_ref| < 0.5%", ok_P2,
       f"blad = {abs(K1_precise-K1_REF)/K1_REF*100:.2f}%" if ok_P1
       else "K*_1 nie znaleziony")

# ─────────────────────────────────────────────────────────────────────────────
# KROK 2: Skan K*_1(alpha) — gesty, alpha in [7.5, 9.5]
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 72)
print("KROK 2: Skan K*_1(alpha) dla alpha w [7.5, 9.5]")
print(f"        n_psi=120, n_eval=2000, rtol=1e-9  (jak p10)")
print("=" * 72)

# Gestsza siatka, szerszy zakres niz p78
alpha_scan = np.array([7.5, 7.8, 8.0, 8.2, 8.4, 8.5, 8.5616,
                        8.6, 8.7, 8.8, 8.9, 9.0, 9.2, 9.5])

results = []
print(f"\n  {'alpha':>8}  {'K*_1':>10}  {'K*_1/K1ref':>12}  status")
print("  " + "-" * 50)

for alpha in alpha_scan:
    print(f"  alpha={alpha:.4f} ...", flush=True)
    K_s, _ = find_K_star(A_GAM, alpha,
                          K_min=3e-4, K_max=0.15, n_K=25, n_psi=120)
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

if len(valid) >= 2:
    alphas_v = np.array([r[0] for r in valid])
    K1s_v    = np.array([r[1] for r in valid])
    diff_K   = K1s_v - K1_REF

    crossings = np.where(np.diff(np.sign(diff_K)))[0]
    if len(crossings) > 0:
        for idx, idx0 in enumerate(crossings):
            i0, i1 = idx0, idx0 + 1
            a0, d0 = alphas_v[i0], diff_K[i0]
            a1, d1 = alphas_v[i1], diff_K[i1]
            ac = a0 + (-d0) * (a1 - a0) / (d1 - d0)
            print(f"\n  Przejscie #{idx+1}: alpha_c = {ac:.5f}")
            print(f"    K*_1({a0:.3f}) = {K1s_v[i0]:.6f},  "
                  f"K*_1({a1:.3f}) = {K1s_v[i1]:.6f},  K1_ref = {K1_REF:.6f}")
            odch = abs(ac - ALPHA_K_REF) / ALPHA_K_REF * 100
            print(f"    odch od alpha_K_ref = {odch:.2f}%")
            if alpha_c_interp is None:
                alpha_c_interp = ac
    else:
        K_min_v, K_max_v = K1s_v[np.isfinite(K1s_v)].min(), K1s_v[np.isfinite(K1s_v)].max()
        print(f"\n  K*_1 nie przecina K1_ref w zakresie [7.5, 9.5].")
        print(f"  K*_1 zakres (finite): [{K_min_v:.6f}, {K_max_v:.6f}]")
        if K_min_v > K1_REF:
            print(f"  -> K*_1 > K1_ref wszedzie (min K*_1={K_min_v:.6f} > {K1_REF})")
            print(f"     Sugestia: alpha_c < 7.5 lub K1_ref jest niedoszacowane")
        elif K_max_v < K1_REF:
            print(f"  -> K*_1 < K1_ref wszedzie: alpha_c > 9.5")

# Testy P3, P4
ok_P3 = alpha_c_interp is not None and np.isfinite(alpha_c_interp)
record("P3: alpha_c znaleziony w [7.5, 9.5]", ok_P3,
       f"alpha_c = {alpha_c_interp:.5f}" if ok_P3 else "nie znaleziony")

ok_P4 = False
if ok_P3:
    odch_P4 = abs(alpha_c_interp - ALPHA_K_REF) / ALPHA_K_REF * 100
    ok_P4 = odch_P4 < 1.0
    record("P4: |alpha_c - 8.5616| < 1%", ok_P4,
           f"alpha_c = {alpha_c_interp:.5f}, odch = {odch_P4:.2f}%")
else:
    record("P4: |alpha_c - 8.5616| < 1%", False, "alpha_c nie znaleziony")

# Test P5: monotycznosc (z akceptacja minimum)
if len(valid) >= 3:
    K_vals_v = [r[1] for r in valid]
    diffs  = [K_vals_v[i+1] - K_vals_v[i] for i in range(len(K_vals_v)-1)]
    n_neg  = sum(1 for d in diffs if d < 0)
    ok_P5  = n_neg >= len(diffs) * 0.7
    record("P5: K*_1(alpha) ma >= 70% krokow malejacych", ok_P5,
           f"{n_neg}/{len(diffs)} spadkow")
else:
    ok_P5 = False
    record("P5: K*_1(alpha) monotonicznie", False, "za malo punktow")

# ─────────────────────────────────────────────────────────────────────────────
# PODSUMOWANIE
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 72)
print("PODSUMOWANIE  p79_alphac_p10accuracy.py")
print("=" * 72)

print(f"\n  PASS: {PASS_COUNT} / {PASS_COUNT + FAIL_COUNT}")
print(f"  FAIL: {FAIL_COUNT} / {PASS_COUNT + FAIL_COUNT}")

print(f"\n  {'Parametr':<30} {'p76':>12} {'p78':>12} {'p79':>12}")
print("  " + "-" * 68)
print(f"  {'n_eval (phi_at_rmax)':<30} {'800':>12} {'1000':>12} {'2000':>12}")
print(f"  {'rtol':<30} {'1e-8':>12} {'1e-8':>12} {'1e-9':>12}")
print(f"  {'n_psi':<30} {'40':>12} {'80':>12} {'120':>12}")
k1_str = f"{K1_precise:.6f}" if np.isfinite(K1_precise) else "---"
print(f"  {'K*_1 (alpha=8.5616)':<30} {'0.010304':>12} {'0.010523':>12} {k1_str:>12}")
k1_err = f"{(K1_precise-K1_REF)/K1_REF*100:+.2f}%" if np.isfinite(K1_precise) else "---"
print(f"  {'  blad vs K1_ref':<30} {'-1.06%':>12} {'+1.04%':>12} {k1_err:>12}")
if ok_P3:
    ac_str = f"{alpha_c_interp:.5f}"
    ac_err = f"{abs(alpha_c_interp-ALPHA_K_REF)/ALPHA_K_REF*100:.2f}%"
else:
    ac_str = "---"
    ac_err = "---"
print(f"  {'alpha_c':<30} {'8.449 (1.3%)':>12} {'--- (>1%)':>12} {ac_str:>12}")
print(f"  {'  blad vs alpha_K':<30} {'1.3%':>12} {'---':>12} {ac_err:>12}")

print()
if ok_P2:
    print("  WNIOSEK P2 PASS: K*_1(p10-acc) ~ K1_ref do < 0.5%.")
    print("  -> p76 i p78 mialy blad numeryczny (n_eval za male).")
    print("  -> p79 z dokladnoscia p10 jest samospojne.")
else:
    print("  WNIOSEK P2 FAIL: K*_1 rozni sie > 0.5% od K1_ref.")
    if np.isfinite(K1_precise):
        err = (K1_precise - K1_REF) / K1_REF * 100
        if abs(err) < 2.0:
            print(f"  -> Rozbiez̈nosc {err:+.2f}% — poprawa vs p78 (+1.04%).")
        else:
            print(f"  -> Rozbiez̈nosc {err:+.2f}% nadal duza.")
    print("  -> Sprawdz r_max, V1, LAM w porownaniu z p10.")

if ok_P4:
    print(f"\n  WNIOSEK P4 PASS: alpha_c = {alpha_c_interp:.5f},  odch = {odch_P4:.2f}%.")
    print("  -> OP-3 Path 2 POTWIERDZONA: alpha_K wyznaczalne z substratu do < 1%.")
elif ok_P3:
    print(f"\n  WNIOSEK P4 FAIL: alpha_c = {alpha_c_interp:.5f},  odch > 1%.")
    print("  -> Rozbiez̈nosc numeryczna nadal istotna.")
else:
    print("\n  WNIOSEK P3 FAIL: alpha_c nie znaleziony w [7.5, 9.5].")
    print("  -> Sprawdz K*_1(alpha) — czy K1_ref lezzy w zakresie K*_1(alpha)?")

# Nowe odkrycie z p78: minimum K*_1(alpha)
print("\n" + "-" * 72)
print("  NOWE ODKRYCIE (z p78 + p79):")
if len(valid) >= 3:
    k_arr = np.array([r[1] for r in valid])
    a_arr = np.array([r[0] for r in valid])
    idx_min = np.argmin(k_arr)
    print(f"  K*_1(alpha) ma MINIMUM:")
    print(f"    alpha_min = {a_arr[idx_min]:.4f},  K*_1_min = {k_arr[idx_min]:.6f}")
    disappeared = [r[0] for r in results if not np.isfinite(r[1])]
    if disappeared:
        print(f"  Soliton znika dla alpha >= {min(disappeared):.1f}")
        print(f"  -> Istnieje alpha_max ~ {min(disappeared):.1f} powyzej ktorego soliton nie istnieje")
