"""
EX64: PRECYZYJNE alpha* I WARUNEK M_raw
========================================
STATUS: LEGACY-TRANSLATIONAL

This script belongs to the older `alpha*` / `M_raw` selection program. It may
still be useful diagnostically, but it should not be treated as the canonical
source for current `nbody` assumptions.

Cele:
  1. Precyzyjne wyznaczenie alpha* (gdzie czyste phi daje r21=207)
     Hipoteza: alpha* = sqrt(6) = 2.4495?

  2. Weryfikacja M_raw ∝ A_tail^4 dla pary (g0*, phi*g0*).
     Jesli M_raw(phi*g0*)/M_raw(g0*) = 207, to M_raw jest lepsza miara masy.

  3. Jaki jest wykladnik potegowy lokalny gamma_eff(g0)?
     d log A_tail / d log g0 — czy g0* ma specjalna wartosc gamma_eff?

  4. Porownanie wszystkich kryteriow kwantowania:
     H1, phi-FP, B-S, M_raw-FP — ktory daje r21?

  5. Czy istnieje warunek "izoenergii": E_e + E_mu = const?
     (E = M_raw lub A_tail^2)

Testy:
  T1: |alpha* - sqrt(6)| < 0.01  (hipoteza alpha* = sqrt(6))
  T2: M_raw(phi*g0*)/M_raw(g0*) in [150, 270]  (M_raw daje ~207?)
  T3: gamma_eff(g0*) in [2.0, 10.0]  (sensowny wykladnik)
  T4: M_raw ∝ A_tail^4 dla n=0 sektora (R^2 > 0.95)
  T5: |E_e + E_mu - cokolwiek_specjalnego| < 10%  (izoenergii?)
"""

import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, curve_fit
import warnings
warnings.filterwarnings('ignore')

PHI      = (1.0 + np.sqrt(5.0)) / 2.0
R21_EXP  = 206.768
R_MAX    = 100.0
G0_STAR  = 1.24771102   # phi-FP z ex62
WIN_LIST  = [16, 22, 28, 36, 46, 58, 72]
WIN_WIDTH = 14.0

# ============================================================
# ODE helpers (parametryzowane alpha)
# ============================================================

def make_ode(alpha, g_off=0.005):
    g_ghost  = np.exp(-1.0 / (2.0 * alpha))
    g_bounce = g_ghost + g_off
    def rhs(r, y):
        g, gp = y
        g  = max(g, g_bounce + 1e-7)
        fg = 1.0 + 2.0 * alpha * np.log(g)
        if abs(fg) < 1e-10: return [gp, 0.0]
        dr = g**2 * (1.0 - g)
        cr = (alpha / g) * gp**2
        if r < 1e-10: return [gp, (dr - cr) / (3.0 * fg)]
        return [gp, (dr - cr - fg * 2.0 * gp / r) / fg]
    def ev(r, y): return y[0] - g_bounce
    ev.terminal = True; ev.direction = -1
    return rhs, ev, g_ghost, g_bounce

def integrate(g0, alpha=2.0, g_off=0.005, maxb=8):
    rhs, ev, _, g_bounce = make_ode(alpha, g_off)
    y0 = [g0, 0.0]; r0 = 0.0; ra, ga, gpa = [], [], []
    for _ in range(maxb + 1):
        sol = solve_ivp(rhs, (r0, R_MAX), y0, events=ev,
                        dense_output=True, rtol=1e-9, atol=1e-11, max_step=0.05)
        ra.append(sol.t); ga.append(sol.y[0]); gpa.append(sol.y[1])
        if sol.status == 1 and len(sol.t_events[0]) > 0:
            rh = sol.t_events[0][-1]; st = sol.sol(rh)
            y0 = [st[0], -st[1]]; r0 = rh
        else: break
    r = np.concatenate(ra); g = np.concatenate(ga); gp = np.concatenate(gpa)
    idx = np.argsort(r)
    return r[idx], g[idx], gp[idx]

def fit_win(r, g, rL, rR):
    mask = (r >= rL) & (r <= rR)
    if np.sum(mask) < 20: return 0., 0., 0.
    rf = r[mask]; df = (g[mask] - 1.0) * rf
    M  = np.column_stack([np.cos(rf), np.sin(rf)])
    bc = np.linalg.lstsq(M, df, rcond=None)[0]
    B, C = bc
    return float(np.sqrt(B**2+C**2)), float(B), float(C)

def A_inf(g0, alpha=2.0, g_off=0.005):
    r, g, _ = integrate(g0, alpha, g_off)
    Av, rv = [], []
    for rL in WIN_LIST:
        if rL + WIN_WIDTH > r[-1]: break
        A, _, _ = fit_win(r, g, rL, rL+WIN_WIDTH)
        if A > 0.002: Av.append(A); rv.append(rL)
    if len(Av) < 3: return Av[-1] if Av else 0.
    try:
        p, _ = curve_fit(lambda x,ai,a: ai*(1+a/x), rv, Av, p0=[Av[-1],0.], maxfev=2000)
        return float(p[0])
    except: return float(Av[-1])

def M_raw(g0, alpha=2.0, g_off=0.005):
    """M_raw = 4pi * int r^2 * (g')^2/2 dr  (kinetyczna)"""
    r, g, gp = integrate(g0, alpha, g_off)
    integrand = r**2 * gp**2 / 2.0
    return float(np.trapezoid(integrand, r) * 4 * np.pi)

def ratio_phi_A(g0_e, alpha=2.0):
    Ae = A_inf(g0_e, alpha); Am = A_inf(PHI*g0_e, alpha)
    return (Am/Ae)**4 if Ae > 1e-6 else np.nan

def ratio_phi_M(g0_e, alpha=2.0):
    Me = M_raw(g0_e, alpha); Mm = M_raw(PHI*g0_e, alpha)
    return Mm/Me if Me > 1e-10 else np.nan

def find_z0(alpha=2.0, lo=1.05, hi=1.55):
    def B_of(g0):
        r, g, _ = integrate(g0, alpha)
        _, B, _ = fit_win(r, g, WIN_LIST[2], WIN_LIST[2]+WIN_WIDTH)
        return B
    gs = np.linspace(lo, hi, 50)
    Bs = [B_of(g) for g in gs]
    for i in range(len(Bs)-1):
        if Bs[i]*Bs[i+1] < 0:
            try: return brentq(B_of, gs[i], gs[i+1], xtol=1e-6)
            except: pass
    return None

# ============================================================
print("="*68)
print("EX64: PRECYZYJNE alpha* I WARUNEK M_raw")
print("="*68)
print(f"  phi = {PHI:.7f},  g0* = {G0_STAR:.7f}")
print(f"  sqrt(6) = {np.sqrt(6):.7f}")
print()

# ============================================================
# CZESC 1: Precyzyjne alpha* (fine scan + brentq)
# ============================================================
print("-"*68)
print("--- Czesc 1: Precyzyjne alpha* ---")
print()

# Fine scan wokol 2.43-2.47
alpha_fine = np.linspace(2.35, 2.55, 20)
print(f"  {'alpha':>7}  {'z0':>9}  {'R_phi':>10}  {'delta%':>9}")
print("  " + "-"*45)
alpha_R_fine = []
for alpha in alpha_fine:
    z0 = find_z0(alpha)
    if z0 is None:
        print(f"  {alpha:7.4f}  z0=None")
        alpha_R_fine.append((alpha, None))
        continue
    R = ratio_phi_A(z0, alpha)
    d = (R - R21_EXP)/R21_EXP*100 if not np.isnan(R) else np.nan
    print(f"  {alpha:7.4f}  {z0:9.5f}  {R:10.3f}  {d:9.3f}%")
    alpha_R_fine.append((alpha, R))

# Find alpha* by brentq
av = [a for a,r in alpha_R_fine if r is not None and not np.isnan(r)]
rv = [r for a,r in alpha_R_fine if r is not None and not np.isnan(r)]
alpha_star = None
for i in range(len(rv)-1):
    if (rv[i]-R21_EXP)*(rv[i+1]-R21_EXP) < 0:
        try:
            alpha_star = brentq(
                lambda a: ratio_phi_A(find_z0(a), a) - R21_EXP if find_z0(a) else 1e6,
                av[i], av[i+1], xtol=1e-4)
        except Exception as e:
            print(f"  brentq: {e}")
            f = (R21_EXP-rv[i])/(rv[i+1]-rv[i])
            alpha_star = av[i] + f*(av[i+1]-av[i])
        break

if alpha_star:
    z0_star = find_z0(alpha_star)
    R_check = ratio_phi_A(z0_star, alpha_star) if z0_star else np.nan
    print()
    print(f"  alpha* = {alpha_star:.6f}")
    print(f"  sqrt(6) = {np.sqrt(6):.6f}  diff = {alpha_star - np.sqrt(6):.6f}")
    print(f"  phi^(3/2) = {PHI**1.5:.6f}")
    print(f"  1 + 1/phi + 1/phi^2 = {1+1/PHI+1/PHI**2:.6f}")
    print(f"  sqrt(2*3) = sqrt(6) = {np.sqrt(6):.6f}")
    print(f"  z0(alpha*) = {z0_star:.6f}")
    print(f"  R_phi(z0, alpha*) = {R_check:.4f}  (vs {R21_EXP})")
print()

# ============================================================
# CZESC 2: M_raw vs A_tail^4 — weryfikacja proporcjonalnosci
# ============================================================
print("-"*68)
print("--- Czesc 2: M_raw vs A_tail^4 (sektor n=0) ---")
print()

g0_scan2 = np.linspace(1.10, 1.55, 20)
M_scan2  = []
A_scan2  = []
print(f"  {'g0':>8}  {'M_raw':>12}  {'A_inf':>10}  {'A^4':>12}  {'M/A^4':>10}")
print("  " + "-"*58)
for g0 in g0_scan2:
    Me = M_raw(g0)
    Ae = A_inf(g0)
    A4 = Ae**4
    ratio_MA = Me/A4 if A4 > 1e-10 else np.nan
    M_scan2.append(Me); A_scan2.append(Ae)
    print(f"  {g0:8.4f}  {Me:12.5f}  {Ae:10.6f}  {A4:12.7f}  {ratio_MA:10.4f}")

# Power law fit: M_raw = c * A_tail^k
logA = np.log(np.array(A_scan2))
logM = np.log(np.array(M_scan2))
try:
    p = np.polyfit(logA, logM, 1)
    k_fit = p[0]; logc_fit = p[1]
    M_pred = np.exp(logc_fit) * np.array(A_scan2)**k_fit
    ss_res = np.sum((logM - (logc_fit + k_fit*logA))**2)
    ss_tot = np.sum((logM - logM.mean())**2)
    R2 = 1 - ss_res/ss_tot
    print()
    print(f"  Power law: M_raw = c * A_tail^k")
    print(f"    k = {k_fit:.4f}  (oczekiwane: 4.0)")
    print(f"    c = {np.exp(logc_fit):.6f}")
    print(f"    R^2 = {R2:.6f}")
except Exception as e:
    k_fit = np.nan; R2 = np.nan
    print(f"  Fit error: {e}")
print()

# ============================================================
# CZESC 3: M_raw(phi*g0*) / M_raw(g0*) — warunek M_raw
# ============================================================
print("-"*68)
print("--- Czesc 3: M_raw condition przy g0* ---")
print()

Me_star = M_raw(G0_STAR)
Mm_star = M_raw(PHI * G0_STAR)
Ae_star = A_inf(G0_STAR)
Am_star = A_inf(PHI * G0_STAR)

print(f"  g0_e = g0*       = {G0_STAR:.7f}")
print(f"  g0_mu = phi*g0*  = {PHI*G0_STAR:.7f}")
print(f"  M_raw(g0_e)      = {Me_star:.6f}")
print(f"  M_raw(g0_mu)     = {Mm_star:.6f}")
print(f"  M_raw ratio      = {Mm_star/Me_star:.4f}  (vs r21={R21_EXP})")
print(f"  A_inf(g0_e)      = {Ae_star:.6f}")
print(f"  A_inf(g0_mu)     = {Am_star:.6f}")
print(f"  (A_mu/A_e)^4     = {(Am_star/Ae_star)**4:.4f}")
print()
print(f"  M_raw = c*A^4 daje: M_mu/M_e = (A_mu/A_e)^4 = {(Am_star/Ae_star)**4:.4f}")
print(f"  Bezposrednie M_mu/M_e          = {Mm_star/Me_star:.4f}")
print(f"  Roznica: {abs(Mm_star/Me_star - (Am_star/Ae_star)**4)/(Am_star/Ae_star)**4*100:.3f}%")
print()

# M_raw FP: znajdz g0** takie ze M_raw(phi*g0**)/M_raw(g0**) = r21
print("  Szukam M_raw FP: M_raw(phi*g0)/M_raw(g0) = 207 ...")

# Quick scan
g0_Msc = np.linspace(1.10, 1.60, 20)
M_ratios = [ratio_phi_M(g0) for g0 in g0_Msc]
print(f"  {'g0':>8}  {'M_ratio':>10}")
Mfp = None
for i, (g0, mr) in enumerate(zip(g0_Msc, M_ratios)):
    if not np.isnan(mr):
        print(f"  {g0:8.4f}  {mr:10.3f}")
    if i > 0 and not np.isnan(M_ratios[i-1]) and not np.isnan(mr):
        if (M_ratios[i-1]-R21_EXP)*(mr-R21_EXP) < 0:
            try:
                g0_Mfp = brentq(lambda x: ratio_phi_M(x)-R21_EXP, g0_Msc[i-1], g0_Msc[i], xtol=1e-5)
                Mfp = g0_Mfp
                Mr_check = ratio_phi_M(g0_Mfp)
                print(f"  M_raw FP: g0 = {g0_Mfp:.7f}  M_ratio = {Mr_check:.4f}")
                print(f"  Porownanie: g0* (A-FP) = {G0_STAR:.7f}, g0 (M-FP) = {g0_Mfp:.7f}")
                print(f"  Roznica A-FP vs M-FP: {G0_STAR - g0_Mfp:.7f}")
            except Exception as e:
                print(f"  brentq: {e}")
print()

# ============================================================
# CZESC 4: gamma_eff(g0) = d log A / d log g0
# ============================================================
print("-"*68)
print("--- Czesc 4: Lokalny wykladnik gamma_eff(g0) = d(lnA)/d(ln g0) ---")
print()

eps = 0.005
g0_gamma = np.linspace(1.12, 1.60, 25)
print(f"  {'g0':>8}  {'A_inf':>10}  {'gamma_eff':>11}  {'sector':>8}")
print("  " + "-"*45)
for g0 in g0_gamma:
    Ap  = A_inf(g0 + eps)
    Am2 = A_inf(g0 - eps)
    A0  = A_inf(g0)
    if Ap > 1e-5 and Am2 > 1e-5:
        gamma_eff = (np.log(Ap) - np.log(Am2)) / (2 * eps) * g0
        nb_label = "n=0" if g0 < 1.61 else "n=1"
        mark = " <== g0*" if abs(g0 - G0_STAR) < 0.015 else ""
        print(f"  {g0:8.4f}  {A0:10.6f}  {gamma_eff:11.4f}  {nb_label:>8}{mark}")

print()

# ============================================================
# CZESC 5: Izoenergii — E_e + E_mu = const?
# ============================================================
print("-"*68)
print("--- Czesc 5: Czy E_e + E_mu ma specjalna wartosc? ---")
print()

# Definiujemy E = A_tail^2 (amplituda kwadratowa) lub M_raw
# Szukamy g0 gdzie A_e^2 + A_mu^2 = min lub const
g0_iso = np.linspace(1.10, 1.50, 25)
E_sum  = []
print(f"  {'g0_e':>8}  {'A_e':>8}  {'A_mu':>8}  {'A_e^2+A_mu^2':>14}  {'A_e*A_mu':>10}")
for g0e in g0_iso:
    Ae = A_inf(g0e); Am = A_inf(PHI*g0e)
    E_sum.append(Ae**2 + Am**2)
    print(f"  {g0e:8.4f}  {Ae:8.5f}  {Am:8.5f}  {Ae**2+Am**2:14.6f}  {Ae*Am:10.6f}")

# Czy E_sum ma minimum?
E_arr = np.array(E_sum)
g0_arr = np.array(g0_iso)
idx_min = np.argmin(E_arr)
print(f"\n  Minimum A_e^2+A_mu^2 przy g0 = {g0_arr[idx_min]:.5f}  (vs g0* = {G0_STAR:.5f})")

# Gradient = 0?
if 1 < idx_min < len(E_arr)-1:
    dE = np.diff(E_arr) / np.diff(g0_arr)
    # interpolate
    for i in range(len(dE)-1):
        if dE[i]*dE[i+1] < 0:
            g0_min = g0_arr[i+1]
            print(f"  Przejscie znaku dE/dg0 przy g0 ~ {g0_min:.5f}")
            break
print()

# ============================================================
# TESTY
# ============================================================
print("="*68)
print("TESTY")
print("="*68)

# T1: |alpha* - sqrt(6)| < 0.01
if alpha_star:
    T1 = abs(alpha_star - np.sqrt(6)) < 0.01
    print(f"  T1: |alpha* - sqrt(6)| < 0.01: {'PASS' if T1 else 'FAIL'}  "
          f"(alpha*={alpha_star:.5f}, sqrt(6)={np.sqrt(6):.5f}, diff={alpha_star-np.sqrt(6):.5f})")
else:
    T1 = False
    print(f"  T1: |alpha* - sqrt(6)| < 0.01: FAIL  (brak alpha*)")

# T2: M_raw(phi*g0*)/M_raw(g0*) in [150, 270]
T2 = 150 <= Mm_star/Me_star <= 270
print(f"  T2: M_ratio in [150,270]:           {'PASS' if T2 else 'FAIL'}  (M_ratio = {Mm_star/Me_star:.2f})")

# T3: gamma_eff(g0*) in [2.0, 10.0]
Ap  = A_inf(G0_STAR + eps); Am2 = A_inf(G0_STAR - eps)
gamma_star = (np.log(Ap) - np.log(Am2)) / (2*eps) * G0_STAR if Ap > 1e-5 and Am2 > 1e-5 else np.nan
T3 = 2.0 <= gamma_star <= 10.0 if not np.isnan(gamma_star) else False
print(f"  T3: gamma_eff(g0*) in [2,10]:       {'PASS' if T3 else 'FAIL'}  (gamma_eff = {gamma_star:.4f})")

# T4: M_raw ∝ A_tail^4 (R^2 > 0.95)
T4 = not np.isnan(R2) and R2 > 0.95
print(f"  T4: M_raw prop A^k (R^2 > 0.95):   {'PASS' if T4 else 'FAIL'}  "
      f"(k={k_fit:.3f}, R^2={R2:.5f})")

# T5: A-FP i M-FP sa bliskie (< 0.01 roznica)
if Mfp:
    T5 = abs(G0_STAR - Mfp) < 0.01
    print(f"  T5: |A-FP - M-FP| < 0.01:          {'PASS' if T5 else 'FAIL'}  "
          f"(A-FP={G0_STAR:.6f}, M-FP={Mfp:.6f}, diff={abs(G0_STAR-Mfp):.6f})")
else:
    T5 = False
    print(f"  T5: |A-FP - M-FP| < 0.01:          FAIL  (brak M-FP)")

n_pass = sum([T1, T2, T3, T4, T5])
print(f"\nWYNIK: {n_pass}/5 testow przeszlo")

print()
print("="*68)
print("WNIOSEK (EX64)")
print("="*68)
print(f"  1. alpha* = {alpha_star:.5f}" if alpha_star else "  1. alpha* = N/A")
print(f"     sqrt(6) = {np.sqrt(6):.5f}  diff = {alpha_star - np.sqrt(6):.5f}" if alpha_star else "")
print(f"  2. M_raw(phi*g0*)/M_raw(g0*) = {Mm_star/Me_star:.3f}  (vs r21={R21_EXP})")
print(f"     (A_mu/A_e)^4 = {(Am_star/Ae_star)**4:.3f}")
print(f"     Roznica: {abs(Mm_star/Me_star - (Am_star/Ae_star)**4)/(Am_star/Ae_star)**4*100:.3f}%")
if not np.isnan(k_fit):
    print(f"  3. M_raw ∝ A_tail^{k_fit:.3f}  (R^2={R2:.4f})")
print(f"  4. gamma_eff(g0*) = {gamma_star:.4f}  (lokalny wykladnik)")
if Mfp:
    print(f"  5. M-FP: g0 = {Mfp:.7f}  (vs A-FP: {G0_STAR:.7f}, diff={G0_STAR-Mfp:.7f})")
print("="*68)
