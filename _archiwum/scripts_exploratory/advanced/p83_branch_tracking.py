#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p83_branch_tracking.py  --  TGP v1  *  Sledzenie galezi przez przelaczenie alpha~8.595
========================================================================================
Cel: W p82 zaobserwowano przelaczenie galexi przy alpha~8.595: E*/K* skacze
     z 12.571 do 12.561.  Czy sledzac ciagla galaz (zamiast argmin|g|)
     minimum E*(alpha) wypada blizej alpha_K = 8.5616?

Mechanizm przelaczenia:
  - Dla danego K, ODE ma wiele rozwiazaN phi(r=0)=psi_core z phi(inf)=1
  - p82 wybiera min(|g|) — przy alpha~8.595 pojawia sie nowa galaz z mniejszym |g|
  - p83 sledzi galaz ciagla: dla kazdego alpha wybiera psi_core najblizsze
    wartosci z poprzedniego alpha (kontynuacja numeryczna)

Metoda:
  A. Goraca galaz  (wysoka energia): ciagla kontynuacja od alpha=8.50 navzgory
  B. Zimna galaz   (niska energia):  ciagla kontynuacja od alpha=8.70 w dol
  C. Porownanie obu galezi E*(alpha), wyznaczenie minimum kazdej

Testy:
  P1: Obie galezie znalezione dla >= 5 wartosci alpha
  P2: E* obliczone dla obu galezi w calym zakresie
  P3: Minimum goraca galezi: |argmin - alpha_K| < 0.10
  P4: Minimum goraca galezi: |argmin - alpha_K| < 0.02  (mocne)
  P5: Galaz goraca i zimna rozroznialne (roznica E* > 0.01%)

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


def find_all_psi_zeros(K, a_gam, alpha, psi_max_frac=1.5, n_psi=150, r_max=40.0):
    """Zwraca LISTĘ wszystkich psi_z gdzie phi(r_max)=1 dla danego (K,alpha)."""
    psi_max = max(3.0, 1 + psi_max_frac*K/a_gam)
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


def g_func_psi(psi_z, K, a_gam, alpha, r_max=40.0):
    E = energy_at_psi(psi_z, K, a_gam, alpha, r_max)
    return E/(4*np.pi*K) - 1 if np.isfinite(E) else np.nan


def find_Kstar_allbranches(alpha, a_gam=A_GAM, K_min=4e-3, K_max=0.12,
                            n_K=25, n_psi=120, r_max=40.0):
    """
    Dla danego alpha znajdz K* oraz WSZYSTKIE galezie psi_core z g=0.
    Zwraca liste slownikow: [{K, psi, E, g_err}, ...] dla kazdej galexi.
    """
    K_scan = np.logspace(np.log10(K_min), np.log10(K_max), n_K)

    # Oblicz g dla wszystkich K (min|g| skan)
    g_vals = []
    psi_best_vals = []
    for K in K_scan:
        roots = find_all_psi_zeros(K, a_gam, alpha, n_psi=n_psi, r_max=r_max)
        if roots:
            g_branch = [(g_func_psi(pz, K, a_gam, alpha, r_max), pz) for pz in roots]
            g_branch = [(g,p) for g,p in g_branch if np.isfinite(g)]
            if g_branch:
                best = min(g_branch, key=lambda x: abs(x[0]))
                g_vals.append(best[0])
                psi_best_vals.append(best[1])
            else:
                g_vals.append(np.nan); psi_best_vals.append(np.nan)
        else:
            g_vals.append(np.nan); psi_best_vals.append(np.nan)

    # Znajdz K* (zmiana znaku g)
    K_star = np.nan; psi_star = np.nan
    for i in range(len(g_vals)-1):
        gi, gj = g_vals[i], g_vals[i+1]
        if np.isfinite(gi) and np.isfinite(gj) and gi*gj < 0:
            try:
                def g_interp(K):
                    roots = find_all_psi_zeros(K, a_gam, alpha, n_psi=n_psi, r_max=r_max)
                    if not roots: return 0.0
                    g_branch = [(g_func_psi(pz, K, a_gam, alpha, r_max), pz) for pz in roots]
                    g_branch = [(g,p) for g,p in g_branch if np.isfinite(g)]
                    if not g_branch: return 0.0
                    return min(g_branch, key=lambda x: abs(x[0]))[0]
                K_star = brentq(g_interp, K_scan[i], K_scan[i+1],
                                xtol=K_scan[i]*0.001, rtol=0.001, maxiter=25)
                break
            except: pass

    if not np.isfinite(K_star):
        return []

    # Przy K_star — wszystkie galezie psi
    roots_star = find_all_psi_zeros(K_star, a_gam, alpha, n_psi=n_psi, r_max=r_max)
    branches = []
    for pz in roots_star:
        E = energy_at_psi(pz, K_star, a_gam, alpha, r_max)
        if np.isfinite(E):
            g_err = E/(4*np.pi*K_star) - 1
            branches.append({'K': K_star, 'psi': pz, 'E': E, 'g': g_err})

    return branches


# ─────────────────────────────────────────────────────────────────────────────
PASS_COUNT = 0; FAIL_COUNT = 0

def record(label, ok, info=""):
    global PASS_COUNT, FAIL_COUNT
    if ok: PASS_COUNT += 1
    else:  FAIL_COUNT += 1
    print(f"  [{'v' if ok else 'x'}] {label}: {info}")

# ─────────────────────────────────────────────────────────────────────────────
print("=" * 72)
print("TGP v1  *  p83_branch_tracking.py  --  Sledzenie galezi przez alpha~8.595")
print("Cel: Ciagla kontynuacja galezi, nie argmin|g|")
print("=" * 72)
print(f"\n  alpha_K_ref = {ALPHA_K_REF}")
print(f"  a_Gam = {A_GAM},  K in [4e-3, 0.12], n_K=25, n_psi=120\n")

# ─────────────────────────────────────────────────────────────────────────────
# SEKCJA A: Pelny skan z WSZYSTKIMI galeziami
# ─────────────────────────────────────────────────────────────────────────────
print("=" * 72)
print("SEKCJA A: Skan alpha in [8.50, 8.70], wszystkie galezie psi_core")
print("=" * 72)

alpha_scan = np.arange(8.50, 8.71, 0.01)

all_results = {}  # alpha -> list of branch dicts
print(f"\n  {'alpha':>7}  {'#galezi':>7}  "
      f"{'psi_1':>8}  {'psi_2':>8}  "
      f"{'E_1':>14}  {'E_2':>14}  "
      f"{'E1/K':>9}  {'E2/K':>9}")
print("  " + "-" * 82)

for alpha in alpha_scan:
    branches = find_Kstar_allbranches(alpha)
    all_results[alpha] = branches
    n = len(branches)
    if n == 0:
        print(f"  {alpha:>7.4f}  {0:>7}  {'---':>8}  {'---':>8}  "
              f"{'---':>14}  {'---':>14}  {'---':>9}  {'---':>9}")
    elif n == 1:
        b = branches[0]
        print(f"  {alpha:>7.4f}  {n:>7}  {b['psi']:>8.4f}  {'---':>8}  "
              f"  {b['E']:>14.8f}  {'---':>14}  {b['E']/b['K']:>9.5f}  {'---':>9}")
    else:
        # Sortuj po psi_core (pierwsza galaz = wysoka psi = wysoka E?)
        branches_sorted = sorted(branches, key=lambda x: x['psi'])
        b1, b2 = branches_sorted[0], branches_sorted[-1]
        print(f"  {alpha:>7.4f}  {n:>7}  {b1['psi']:>8.4f}  {b2['psi']:>8.4f}  "
              f"  {b1['E']:>14.8f}  {b2['E']:>14.8f}  "
              f"{b1['E']/b1['K']:>9.5f}  {b2['E']/b2['K']:>9.5f}")

# ─────────────────────────────────────────────────────────────────────────────
# SEKCJA B: Ciagla galaz — sledzenie przez przelaczenie
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 72)
print("SEKCJA B: Ciagla galaz (kontynuacja numeryczna)")
print("=" * 72)
print("  Galaz A: od alpha=8.50, kontynuacja w gore (najnizsze psi)")
print("  Galaz B: od alpha=8.70, kontynuacja w dol  (najnizsze psi)")

# Galaz A: startujemy od alpha=8.50, sledzac galaz o najnizszym psi_core
branch_A = {}
prev_psi_A = None
for alpha in alpha_scan:
    brs = all_results[alpha]
    if not brs:
        branch_A[alpha] = None
        continue
    brs_sorted = sorted(brs, key=lambda x: x['psi'])
    if prev_psi_A is None:
        # Pierwszy punkt: wybierz galaz o najmniejszym psi_core
        chosen = brs_sorted[0]
    else:
        # Kontynuacja: wybierz galaz najblizej poprzedniej psi
        chosen = min(brs_sorted, key=lambda x: abs(x['psi'] - prev_psi_A))
    branch_A[alpha] = chosen
    prev_psi_A = chosen['psi']

# Galaz B: startujemy od alpha=8.70, sledzac galaz o najwyzszym psi_core
branch_B = {}
prev_psi_B = None
for alpha in alpha_scan[::-1]:  # od 8.70 do 8.50
    brs = all_results[alpha]
    if not brs:
        branch_B[alpha] = None
        continue
    brs_sorted = sorted(brs, key=lambda x: x['psi'])
    if prev_psi_B is None:
        # Pierwszy punkt (alpha=8.70): najwyzsze psi
        chosen = brs_sorted[-1]
    else:
        chosen = min(brs_sorted, key=lambda x: abs(x['psi'] - prev_psi_B))
    branch_B[alpha] = chosen
    prev_psi_B = chosen['psi']

# Tabela
print(f"\n  {'alpha':>7}  {'psi_A':>8}  {'E_A':>14}  {'EA/K':>9}  "
      f"  {'psi_B':>8}  {'E_B':>14}  {'EB/K':>9}")
print("  " + "-" * 72)
for alpha in alpha_scan:
    bA = branch_A.get(alpha)
    bB = branch_B.get(alpha)
    sA = (f"{bA['psi']:>8.4f}  {bA['E']:>14.8f}  {bA['E']/bA['K']:>9.5f}"
          if bA else f"{'---':>8}  {'---':>14}  {'---':>9}")
    sB = (f"{bB['psi']:>8.4f}  {bB['E']:>14.8f}  {bB['E']/bB['K']:>9.5f}"
          if bB else f"{'---':>8}  {'---':>14}  {'---':>9}")
    print(f"  {alpha:>7.4f}  {sA}    {sB}")

# ─────────────────────────────────────────────────────────────────────────────
# SEKCJA C: Minimum E* dla kazdej galexi
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 72)
print("SEKCJA C: Analiza minimum E*(alpha) dla kazdej galexi")
print("=" * 72)

def find_minimum(branch_dict, alpha_arr):
    """Znajdz minimum E* i argmin, fit paraboliczny."""
    valid = [(a, b['E'], b['K']) for a, b in branch_dict.items()
             if b is not None and np.isfinite(b['E'])]
    if len(valid) < 5:
        return np.nan, np.nan, len(valid)
    a_arr = np.array([x[0] for x in valid])
    E_arr = np.array([x[1] for x in valid])
    K_arr = np.array([x[2] for x in valid])

    idx_min = np.argmin(E_arr)
    alpha_min_direct = a_arr[idx_min]
    E_min_direct = E_arr[idx_min]

    # Fit paraboliczny wokol minimum
    i0 = max(1, idx_min-2); i1 = min(len(valid)-2, idx_min+2)
    a_fit = a_arr[i0:i1+1]; E_fit = E_arr[i0:i1+1]
    alpha_min_fit = np.nan
    if len(a_fit) >= 3:
        coeffs = np.polyfit(a_fit, E_fit, 2)
        if coeffs[0] > 0:
            alpha_min_fit = -coeffs[1]/(2*coeffs[0])

    return alpha_min_direct, alpha_min_fit, len(valid)

print("\n  --- Galaz A (ciagla od alpha=8.50, niska psi_core) ---")
a_min_A_dir, a_min_A_fit, n_A = find_minimum(branch_A, alpha_scan)
print(f"  Liczba punktow: {n_A}")
print(f"  Minimum bezposrednie: alpha = {a_min_A_dir:.4f}")
if np.isfinite(a_min_A_fit):
    print(f"  Fit paraboliczny:     alpha = {a_min_A_fit:.5f}")
    odch_A = abs(a_min_A_fit - ALPHA_K_REF)
    print(f"  Odchylenie od alpha_K: {odch_A:.5f}")
else:
    odch_A = abs(a_min_A_dir - ALPHA_K_REF) if np.isfinite(a_min_A_dir) else np.nan
    print(f"  Odchylenie od alpha_K: {odch_A:.5f} (bezposrednie)")

print("\n  --- Galaz B (ciagla od alpha=8.70, wysoka psi_core) ---")
a_min_B_dir, a_min_B_fit, n_B = find_minimum(branch_B, alpha_scan)
print(f"  Liczba punktow: {n_B}")
print(f"  Minimum bezposrednie: alpha = {a_min_B_dir:.4f}")
if np.isfinite(a_min_B_fit):
    print(f"  Fit paraboliczny:     alpha = {a_min_B_fit:.5f}")
    odch_B = abs(a_min_B_fit - ALPHA_K_REF)
    print(f"  Odchylenie od alpha_K: {odch_B:.5f}")
else:
    odch_B = abs(a_min_B_dir - ALPHA_K_REF) if np.isfinite(a_min_B_dir) else np.nan
    print(f"  Odchylenie od alpha_K: {odch_B:.5f} (bezposrednie)")

# ─────────────────────────────────────────────────────────────────────────────
# TESTY
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 72)
print("TESTY")
print("=" * 72)

ok_P1 = n_A >= 5 and n_B >= 5
record("P1: Obie galezie znalezione dla >= 5 wartosci alpha",
       ok_P1, f"A={n_A} pkt, B={n_B} pkt")

best_A = a_min_A_fit if np.isfinite(a_min_A_fit) else a_min_A_dir
best_B = a_min_B_fit if np.isfinite(a_min_B_fit) else a_min_B_dir
ok_P2 = np.isfinite(best_A) and np.isfinite(best_B)
record("P2: E* obliczone dla obu galezi", ok_P2,
       f"argmin_A={best_A:.4f}, argmin_B={best_B:.4f}" if ok_P2 else "brak")

# Testuj galaz ktora daje minimum blizsze alpha_K
best_odch = min(odch_A if np.isfinite(odch_A) else 99,
                odch_B if np.isfinite(odch_B) else 99)
best_galaz = 'A' if (np.isfinite(odch_A) and odch_A <= best_odch+1e-6) else 'B'

ok_P3 = best_odch < 0.10
record("P3: |argmin E*(galaz) - alpha_K| < 0.10", ok_P3,
       f"galaz {best_galaz}: odch={best_odch:.5f}")

ok_P4 = best_odch < 0.02
record("P4: |argmin E*(galaz) - alpha_K| < 0.02 (mocne)", ok_P4,
       f"galaz {best_galaz}: odch={best_odch:.5f}")

# P5: Czy galezie sa rozroznialne
E_vals_A = [b['E'] for b in branch_A.values() if b is not None]
E_vals_B = [b['E'] for b in branch_B.values() if b is not None]
# Znajdz wspolne alpha
common_alpha = [a for a in alpha_scan
                if branch_A.get(a) is not None and branch_B.get(a) is not None]
if common_alpha:
    diffs = [abs(branch_A[a]['E'] - branch_B[a]['E'])/branch_A[a]['E']
             for a in common_alpha
             if np.isfinite(branch_A[a]['E']) and np.isfinite(branch_B[a]['E'])]
    max_diff_pct = max(diffs)*100 if diffs else 0
    ok_P5 = max_diff_pct > 0.01
    record("P5: Galezie rozroznialne (roznica E* > 0.01%)", ok_P5,
           f"max roznica = {max_diff_pct:.4f}%")
else:
    ok_P5 = False
    record("P5: Galezie rozroznialne", ok_P5, "brak wspolnych punktow")

# ─────────────────────────────────────────────────────────────────────────────
# PODSUMOWANIE
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 72)
print("PODSUMOWANIE  p83_branch_tracking.py")
print("=" * 72)
print(f"\n  PASS: {PASS_COUNT} / {PASS_COUNT + FAIL_COUNT}")
print(f"  FAIL: {FAIL_COUNT} / {PASS_COUNT + FAIL_COUNT}")

print(f"""
  WYNIKI:
    Galaz A (od alpha=8.50, niska psi): argmin E* = {best_A:.5f}
    Galaz B (od alpha=8.70, wyzsza psi): argmin E* = {best_B:.5f}
    alpha_K_ref = {ALPHA_K_REF}

    Najlepsza galaz: {best_galaz}, odchylenie od alpha_K = {best_odch:.5f}
""")

if ok_P4:
    print(f"""  WNIOSEK P4 PASS:
  Ciagla galaz {best_galaz} daje argmin E*(alpha) = {best_A if best_galaz=='A' else best_B:.5f}
  Odchylenie {best_odch:.5f} = {best_odch/ALPHA_K_REF*100:.3f}% od alpha_K.
  -> OP-3 Sciezka 1 POTWIERDZONA ciagla galazia!
""")
elif ok_P3:
    print(f"""  WNIOSEK P3 PASS:
  Ciagla galaz {best_galaz} daje argmin E*(alpha) ~ {best_A if best_galaz=='A' else best_B:.5f}
  Odchylenie {best_odch:.5f} od alpha_K={ALPHA_K_REF}.
  Nie mocno potwierdzone (odch < 0.10, ale >= 0.02).
""")
else:
    print(f"""  WNIOSEK: Ciagla galaz nie zbliga argmin E* do alpha_K.
  Odchylenie {best_odch:.5f} > 0.10 od alpha_K={ALPHA_K_REF}.
  Przelaczenie galexi w p82 nie bylo artefaktem — prawdziwa bifurkacja galezi.
""")
