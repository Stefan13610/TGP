#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p84_branch1_scan.py  --  TGP v1  *  Minimum E* podgalezi B1 (alpha < 8.595)
============================================================================
Cel: W p83 stwierdzono ze podgalaz B1 (psi_core~1.24, E*/K*~12.571)
     wystepuje dla alpha <= 8.595 i tam E*(alpha) rosnie w gore
     (minimum jest NIZEJ, tj. alpha < 8.50).
     Pytanie: gdzie jest minimum E* na podgalezi B1?
     Jesli argmin B1 ~ alpha_K = 8.5616 -> OP-3 Sciezka 1 POTWIERDZONA.

Kontekst (p83):
  - B1: alpha in [8.50, 8.59], E*/K* in [12.5707, 12.5713], rosnace
  - B2: alpha in [8.60, 8.71], E*/K* in [12.5614, 12.5615], minimum ~8.65
  - Przeskak B1->B2 przy alpha~8.595
  - B1 ciagnie sie w dol (alpha < 8.50); minimum B1 nieznane

Metoda:
  - Skan alpha in [5.0, 8.60] krok 0.05 (wstepny)
  - Dla kazdego alpha: znajdz K* i WSZYSTKIE galezie psi_core
  - Wybierz galaz B1: ta z psi_core ~ 1.24 i E*/K* ~ 12.57
  - Jesl wiele galezi, selekcja przez: E*/K* in [12.45, 12.60] (B-range)
  - Znajdz minimum E* na tej galezi

Testy:
  P1: B1 znaleziona dla >= 5 wartosci alpha < 8.50
  P2: E*(alpha) obliczone dla B1 w calym zakresie
  P3: Minimum B1: |argmin - alpha_K| < 0.10
  P4: Minimum B1: |argmin - alpha_K| < 0.02 (mocne - OP-3 Sciezka 1)
  P5: E* ma wyrazne minimum (nie monotoniczna) na B1

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
ALPHA_K_REF = 8.5616
A_GAM       = 0.040
K1_REF      = 0.010414
LAM         = 5.501357e-06
GAMMA_V     = 1.0

def V_mod(p):
    return GAMMA_V/3*p**3 - GAMMA_V/4*p**4 + LAM/6*(p-1)**6
def dV_mod(p):
    return GAMMA_V*p**2 - GAMMA_V*p**3 + LAM*(p-1)**5
V1 = GAMMA_V/3 - GAMMA_V/4


def ode_rhs(r, y, alpha):
    phi, dphi = y; phi = max(phi, 1e-10)
    kfac = 1.0 + alpha/phi
    return [dphi, dV_mod(phi)/kfac + alpha*dphi**2/(2*phi**2*kfac) - 2/r*dphi]


def _phi_at_rmax(psi, K, a_gam, alpha, r_max=40.0, n_eval=2000):
    dphi0 = -K/a_gam**2
    def ev(r, y): return y[0] - 1e-5
    ev.terminal = True; ev.direction = -1
    r_ev = a_gam*(r_max/a_gam)**np.linspace(0, 1, n_eval)
    try:
        sol = solve_ivp(lambda r,y: ode_rhs(r,y,alpha), [a_gam, r_max],
                        [psi, dphi0], method='DOP853', rtol=1e-9, atol=1e-11,
                        t_eval=r_ev, events=[ev])
        return float(sol.y[0,-1]) if sol.t[-1] >= r_max*0.99 else np.nan
    except: return np.nan


def energy_at_psi(psi_z, K, a_gam, alpha, r_max=40.0):
    dphi0 = -K/a_gam**2
    r_ev = a_gam*(r_max/a_gam)**np.linspace(0, 1, 4000)
    try:
        sol = solve_ivp(lambda r,y: ode_rhs(r,y,alpha), [a_gam, r_max],
                        [psi_z, dphi0], method='DOP853', rtol=1e-9, atol=1e-11, t_eval=r_ev)
        r = sol.t; phi = np.maximum(sol.y[0], 1e-10); dphi = sol.y[1]
        Ek = 4*np.pi*np.trapezoid(0.5*dphi**2*(1+alpha/phi)*r**2, r)
        Ep = 4*np.pi*np.trapezoid((V_mod(phi)-V1)*r**2, r)
        return Ek + Ep
    except: return np.nan


def find_all_psi_zeros(K, a_gam, alpha, n_psi=150, r_max=40.0):
    psi_max = max(3.0, 1.0 + 2.0*K/a_gam)
    psis = np.linspace(1.001, min(psi_max, 400), n_psi)
    Fv = [_phi_at_rmax(p, K, a_gam, alpha, r_max) - 1 for p in psis]
    roots = []
    for i in range(len(Fv)-1):
        if np.isfinite(Fv[i]) and np.isfinite(Fv[i+1]) and Fv[i]*Fv[i+1] < 0:
            try:
                pz = brentq(lambda p: _phi_at_rmax(p, K, a_gam, alpha, r_max)-1,
                            psis[i], psis[i+1], xtol=1e-7, rtol=1e-7, maxiter=60)
                roots.append(pz)
            except: pass
    return roots


def find_Kstar_branches(alpha, a_gam=A_GAM, K_min=4e-3, K_max=0.12,
                         n_K=25, n_psi=120, r_max=40.0):
    """Znajdz K* i WSZYSTKIE galezie psi_core z g~0."""
    K_scan = np.logspace(np.log10(K_min), np.log10(K_max), n_K)
    g_vals = []
    for K in K_scan:
        roots = find_all_psi_zeros(K, a_gam, alpha, n_psi=n_psi, r_max=r_max)
        if roots:
            # Znajdz galaz o najblizszej g=0 wsrod wszystkich
            best_g = np.nan
            for pz in roots:
                E = energy_at_psi(pz, K, a_gam, alpha, r_max)
                if np.isfinite(E):
                    g = E/(4*np.pi*K) - 1
                    if np.isnan(best_g) or abs(g) < abs(best_g):
                        best_g = g
            g_vals.append(best_g if np.isfinite(best_g) else np.nan)
        else:
            g_vals.append(np.nan)

    K_star = np.nan
    for i in range(len(g_vals)-1):
        gi, gj = g_vals[i], g_vals[i+1]
        if np.isfinite(gi) and np.isfinite(gj) and gi*gj < 0:
            try:
                def g_interp(K):
                    roots = find_all_psi_zeros(K, a_gam, alpha, n_psi=n_psi)
                    if not roots: return 0.0
                    best_g = np.nan
                    for pz in roots:
                        E = energy_at_psi(pz, K, a_gam, alpha, r_max)
                        if np.isfinite(E):
                            g = E/(4*np.pi*K)-1
                            if np.isnan(best_g) or abs(g) < abs(best_g):
                                best_g = g
                    return float(best_g) if np.isfinite(best_g) else 0.0
                K_star = brentq(g_interp, K_scan[i], K_scan[i+1],
                                xtol=K_scan[i]*0.001, rtol=0.001, maxiter=25)
                break
            except: pass

    if not np.isfinite(K_star):
        return []

    # Przy K_star: wszystkie galezie
    roots_star = find_all_psi_zeros(K_star, a_gam, alpha, n_psi=n_psi, r_max=r_max)
    branches = []
    for pz in roots_star:
        E = energy_at_psi(pz, K_star, a_gam, alpha, r_max)
        if np.isfinite(E):
            g = E/(4*np.pi*K_star) - 1
            EoK = E/K_star
            branches.append({'K': K_star, 'psi': pz, 'E': E, 'g': g, 'EoK': EoK})
    return branches


# Zakres B1: E/K in [12.40, 12.65] (fizyczny soliton)
B1_EOK_MIN = 12.40
B1_EOK_MAX = 12.65

def select_B1(branches, prev_psi=None):
    """Wybierz galaz B1: E/K in [12.40,12.65], najblizej prev_psi jezeli znana."""
    b1_cands = [b for b in branches
                if B1_EOK_MIN <= b['EoK'] <= B1_EOK_MAX]
    if not b1_cands:
        return None
    if prev_psi is None:
        return min(b1_cands, key=lambda b: abs(b['g']))
    return min(b1_cands, key=lambda b: abs(b['psi'] - prev_psi))


# ─────────────────────────────────────────────────────────────────────────────
PASS_COUNT = 0; FAIL_COUNT = 0

def record(label, ok, info=""):
    global PASS_COUNT, FAIL_COUNT
    if ok: PASS_COUNT += 1
    else:  FAIL_COUNT += 1
    print(f"  [{'v' if ok else 'x'}] {label}: {info}")

# ─────────────────────────────────────────────────────────────────────────────
print("=" * 72)
print("TGP v1  *  p84_branch1_scan.py  --  Minimum E* galezi B1 (alpha<8.595)")
print("Zakres: alpha in [5.0, 8.60] krok 0.05;  selekcja E*/K* in [12.40,12.65]")
print("=" * 72)
print(f"\n  alpha_K_ref = {ALPHA_K_REF}")
print(f"  a_Gam = {A_GAM},  K in [4e-3, 0.12], n_K=25, n_psi=120\n")

# ─────────────────────────────────────────────────────────────────────────────
# SEKCJA A: Wstepny skan krok 0.05
# ─────────────────────────────────────────────────────────────────────────────
print("=" * 72)
print("SEKCJA A: Skan alpha in [5.0, 8.60] krok 0.05 — galaz B1")
print("=" * 72)

alpha_coarse = np.arange(5.00, 8.61, 0.05)

b1_coarse = {}
prev_psi = None

print(f"\n  {'alpha':>7}  {'K*':>10}  {'psi_B1':>8}  {'E_B1':>14}  {'E1/K':>9}  {'g':>8}")
print("  " + "-" * 62)

for alpha in alpha_coarse:
    brs = find_Kstar_branches(alpha)
    b1 = select_B1(brs, prev_psi)
    b1_coarse[alpha] = b1
    if b1 is not None:
        prev_psi = b1['psi']
        print(f"  {alpha:>7.4f}  {b1['K']:>10.6f}  {b1['psi']:>8.4f}  "
              f"{b1['E']:>14.8f}  {b1['EoK']:>9.5f}  {b1['g']:>8.5f}")
    else:
        print(f"  {alpha:>7.4f}  {'---':>10}  {'---':>8}  {'---':>14}  {'---':>9}  {'---':>8}")

# ─────────────────────────────────────────────────────────────────────────────
# SEKCJA B: Analiza minimum
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 72)
print("SEKCJA B: Analiza minimum E*(alpha) galezi B1")
print("=" * 72)

valid = [(a, b['E'], b['K'], b['psi']) for a, b in b1_coarse.items()
         if b is not None and np.isfinite(b['E'])]
n_B1 = len(valid)
print(f"\n  Liczba punktow B1: {n_B1}")

alpha_min_B1_dir = np.nan
alpha_min_B1_fit = np.nan
odch_B1 = np.nan

if n_B1 >= 5:
    a_arr = np.array([x[0] for x in valid])
    E_arr = np.array([x[1] for x in valid])

    idx_min = np.argmin(E_arr)
    alpha_min_B1_dir = a_arr[idx_min]
    print(f"  Minimum bezposrednie: alpha = {alpha_min_B1_dir:.4f},  "
          f"E* = {E_arr[idx_min]:.8f}")

    # Fit paraboliczny
    i0 = max(1, idx_min-2); i1 = min(n_B1-2, idx_min+2)
    a_fit = a_arr[i0:i1+1]; E_fit = E_arr[i0:i1+1]
    if len(a_fit) >= 3:
        coeffs = np.polyfit(a_fit, E_fit, 2)
        if coeffs[0] > 0:
            alpha_min_B1_fit = -coeffs[1]/(2*coeffs[0])
            E_min_fit = np.polyval(coeffs, alpha_min_B1_fit)
            print(f"  Fit paraboliczny:     alpha = {alpha_min_B1_fit:.5f},  "
                  f"E* = {E_min_fit:.8f}")

    best = alpha_min_B1_fit if np.isfinite(alpha_min_B1_fit) else alpha_min_B1_dir
    odch_B1 = abs(best - ALPHA_K_REF)
    print(f"\n  Najlepszy argmin B1: alpha = {best:.5f}")
    print(f"  Odchylenie od alpha_K={ALPHA_K_REF}: {odch_B1:.5f}")

    # Czy E* ma minimum (nie monotoniczna)?
    dE = np.diff(E_arr)
    # Minimum istnieje jesli sa zmiany znaku
    sign_changes = np.sum(dE[:-1]*dE[1:] < 0)
    monotone_left = np.all(dE[:idx_min] <= 0) if idx_min > 0 else True
    has_min = (idx_min > 0 and idx_min < n_B1-1)
    print(f"  Monotoniczne do minimum: {monotone_left}")
    print(f"  Minimum wewnatrz zakresu (nie na krawedzi): {has_min}")

    # Sprawdz czy E* bedzie dalej malen (jesli na lewej krawedzi)
    if idx_min == 0:
        print("  UWAGA: minimum na lewej krawedzi alpha=5.0 -- moze byc poza zakresem!")

# ─────────────────────────────────────────────────────────────────────────────
# SEKCJA C: Jesli minimum w srodku — precyzyjny skan ±0.5 wokol argmin
# ─────────────────────────────────────────────────────────────────────────────
alpha_min_B1_precise = np.nan

if np.isfinite(alpha_min_B1_dir) and 5.05 <= alpha_min_B1_dir <= 8.55:
    print("\n" + "=" * 72)
    print(f"SEKCJA C: Precyzyjny skan wokol minimum (alpha={alpha_min_B1_dir:.2f} +/-0.5, krok 0.02)")
    print("=" * 72)

    a_lo_p = max(5.0, alpha_min_B1_dir - 0.5)
    a_hi_p = min(8.60, alpha_min_B1_dir + 0.5)
    alpha_fine = np.arange(a_lo_p, a_hi_p + 0.01, 0.02)

    b1_fine = {}
    prev_psi2 = None
    print(f"\n  {'alpha':>7}  {'psi_B1':>8}  {'E_B1':>14}  {'E1/K':>9}")
    print("  " + "-" * 44)
    for alpha in alpha_fine:
        brs = find_Kstar_branches(alpha)
        b1 = select_B1(brs, prev_psi2)
        b1_fine[alpha] = b1
        if b1 is not None:
            prev_psi2 = b1['psi']
            print(f"  {alpha:>7.4f}  {b1['psi']:>8.4f}  {b1['E']:>14.8f}  {b1['EoK']:>9.5f}")
        else:
            print(f"  {alpha:>7.4f}  {'---':>8}  {'---':>14}  {'---':>9}")

    valid_f = [(a, b['E']) for a, b in b1_fine.items()
               if b is not None and np.isfinite(b['E'])]
    if len(valid_f) >= 5:
        a_f = np.array([x[0] for x in valid_f])
        E_f = np.array([x[1] for x in valid_f])
        idx_f = np.argmin(E_f)
        alpha_min_B1_precise = a_f[idx_f]
        print(f"\n  Minimum precyzyjne (bezposrednie): alpha = {alpha_min_B1_precise:.4f}")
        i0 = max(1,idx_f-2); i1 = min(len(valid_f)-2,idx_f+2)
        a_pfit = a_f[i0:i1+1]; E_pfit = E_f[i0:i1+1]
        if len(a_pfit) >= 3:
            cf = np.polyfit(a_pfit, E_pfit, 2)
            if cf[0] > 0:
                alpha_min_B1_precise = -cf[1]/(2*cf[0])
                print(f"  Minimum precyzyjne (fit para.):   alpha = {alpha_min_B1_precise:.5f}")
        odch_B1 = abs(alpha_min_B1_precise - ALPHA_K_REF)
        print(f"  Odchylenie od alpha_K: {odch_B1:.5f}")

# ─────────────────────────────────────────────────────────────────────────────
# TESTY
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 72)
print("TESTY")
print("=" * 72)

ok_P1 = n_B1 >= 5
record("P1: B1 znaleziona dla >= 5 wartosci alpha", ok_P1, f"{n_B1} punktow")

ok_P2 = n_B1 >= 10
record("P2: E*(alpha) obliczone dla >= 10 wartosci", ok_P2, f"{n_B1} punktow")

best_argmin = (alpha_min_B1_precise if np.isfinite(alpha_min_B1_precise)
               else alpha_min_B1_fit if np.isfinite(alpha_min_B1_fit)
               else alpha_min_B1_dir)

ok_P3 = np.isfinite(odch_B1) and odch_B1 < 0.10
record("P3: |argmin B1 - alpha_K| < 0.10", ok_P3,
       f"odch={odch_B1:.5f}" if np.isfinite(odch_B1) else "brak minimum")

ok_P4 = np.isfinite(odch_B1) and odch_B1 < 0.02
record("P4: |argmin B1 - alpha_K| < 0.02 (mocne - OP-3 Sciezka 1)", ok_P4,
       f"odch={odch_B1:.5f}" if np.isfinite(odch_B1) else "brak minimum")

# P5: czy E* ma minimum wewnatrz zakresu (nie monotoniczna)
if valid and len(valid) >= 3:
    a_arr = np.array([x[0] for x in valid])
    E_arr = np.array([x[1] for x in valid])
    idx_min = np.argmin(E_arr)
    ok_P5 = (idx_min > 0 and idx_min < len(valid)-1)
    record("P5: E*(B1) ma minimum wewnatrz zakresu (nie monot.)", ok_P5,
           f"argmin przy alpha={a_arr[idx_min]:.4f} (idx={idx_min}/{len(valid)-1})")
else:
    ok_P5 = False
    record("P5: E*(B1) ma minimum wewnatrz zakresu", ok_P5, "za malo danych")

# ─────────────────────────────────────────────────────────────────────────────
# PODSUMOWANIE
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 72)
print("PODSUMOWANIE  p84_branch1_scan.py")
print("=" * 72)
print(f"\n  PASS: {PASS_COUNT} / {PASS_COUNT + FAIL_COUNT}")
print(f"  FAIL: {FAIL_COUNT} / {PASS_COUNT + FAIL_COUNT}")

if np.isfinite(best_argmin):
    print(f"""
  WYNIKI:
    argmin E*(galaz B1) = {best_argmin:.5f}
    alpha_K_ref         = {ALPHA_K_REF}
    Odchylenie          = {odch_B1:.5f}
""")

if ok_P4:
    print(f"""  WNIOSEK P4 PASS:
  Galaz B1 daje argmin E*(alpha) = {best_argmin:.5f}
  (odchylenie {odch_B1:.5f} = {odch_B1/ALPHA_K_REF*100:.3f}% od alpha_K)
  -> OP-3 Sciezka 1 POTWIERDZONA na galezi B1!
""")
elif ok_P3:
    print(f"""  WNIOSEK P3 PASS:
  Galaz B1: argmin E* = {best_argmin:.5f}, odch = {odch_B1:.5f}
  Blisko alpha_K (< 0.10), ale P4 nie spelniony.
""")
elif ok_P5 and not ok_P3:
    print(f"""  WNIOSEK: Galaz B1 ma minimum wewnatrz zakresu, ale offset = {odch_B1:.5f} > 0.10.
  OP-3 Sciezka 1 nie potwierdzona.
""")
elif not ok_P5:
    print("""  WNIOSEK: Galaz B1 monotoniczna w skanowanym zakresie.
  Minimum poza zakresem [5.0, 8.60] lub galaz nie istnieje.
""")
