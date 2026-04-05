"""
EX68: MASA RENORMALIZOWANA M_ren = M_total(R) - pi*A_inf^2*R

Hipoteza: Ogon solitonu daje dywergencje liniowa M_total ~ pi*A^2*R.
Subtrakcja: M_ren = M_total(R) - pi*A_inf^2 * R -> stala ~ c*A^4.

Jesli M_ren ~ c*A^4, to M_ren(phi*g0*) / M_ren(g0*) = (A_mu/A_e)^4 = r21.

Tests:
  T1: M_total(R)/pi/A^2/R -> 1  dla R->inf (liniowe zdominowanie)
  T2: M_ren(R) jest zbiezne (wzgledny zakres < 50%)
  T3: M_ren ~ A^k, k in [3.5, 4.5] w sektorze n=0
  T4: |M_ren ratio - 207| < 20%
  T5: M_ren/A^4 CV < 40%
"""

import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import curve_fit, brentq
from scipy.stats import linregress
from numpy.linalg import lstsq
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# STALE TGP (z ex62b)
# ============================================================
ALPHA     = 2.0
G_OFF     = 0.005
G_GHOST   = np.exp(-1.0 / (2.0 * ALPHA))    # ~0.77880
G_BOUNCE  = G_GHOST + G_OFF                  # ~0.78380
PHI       = (1.0 + np.sqrt(5.0)) / 2.0      # 1.6180340
G0_STAR   = 1.24771102
R21_EXP   = 206.768
R_MAX_SOL = 200.0   # dla integracji ODE
WIN_LIST  = [16, 22, 28, 36, 46, 58, 72]
WIN_WIDTH = 14.0

# ============================================================
# ODE z elastycznym odbiciem od g* (z ex62b)
# ============================================================
def make_ode(alpha=ALPHA, g_off=G_OFF):
    g_bounce = np.exp(-1.0/(2.0*alpha)) + g_off
    def rhs(r, y):
        g, gp = y
        g  = max(g, g_bounce + 1e-7)
        fg = 1.0 + 2.0*alpha*np.log(g)
        if abs(fg) < 1e-10: return [gp, 0.0]
        dr = g**2*(1.0-g)
        cr = (alpha/g)*gp**2
        if r < 1e-10: return [gp, (dr-cr)/(3.0*fg)]
        return [gp, (dr-cr-fg*2.0*gp/r)/fg]
    def ev(r, y): return y[0] - g_bounce
    ev.terminal = True; ev.direction = -1
    return rhs, ev

def integrate(g0, alpha=ALPHA, g_off=G_OFF, maxb=8, r_end=R_MAX_SOL):
    """Rozwiaz ODE z elastycznymi odbiicami od granicy g*."""
    rhs, ev = make_ode(alpha, g_off)
    y0 = [g0, 0.0]; r0 = 1e-10; ra, ga, gpa = [], [], []
    for _ in range(maxb+1):
        sol = solve_ivp(rhs, (r0, r_end), y0, events=ev,
                        dense_output=True, rtol=1e-9, atol=1e-11, max_step=0.05)
        ra.append(sol.t); ga.append(sol.y[0]); gpa.append(sol.y[1])
        if sol.status == 1 and len(sol.t_events[0]) > 0:
            rh = sol.t_events[0][-1]
            st = sol.sol(rh)
            y0 = [st[0], -st[1]]   # odbicie elastyczne: gp -> -gp
            r0 = rh
        else:
            break
    r   = np.concatenate(ra)
    g   = np.concatenate(ga)
    gp  = np.concatenate(gpa)
    idx = np.argsort(r)
    return r[idx], g[idx], gp[idx]

# ============================================================
# EKSTRAKCJA A_inf (okienkowa, z ex62b)
# ============================================================
def fit_win(r, g, rL, rR):
    mask = (r >= rL) & (r <= rR)
    if np.sum(mask) < 20: return 0., 0., 0.
    rf = r[mask]; df = (g[mask]-1.0)*rf
    M  = np.column_stack([np.cos(rf), np.sin(rf)])
    bc = np.linalg.lstsq(M, df, rcond=None)[0]
    B, C = bc
    return float(np.sqrt(B**2+C**2)), float(B), float(C)

def A_inf_val(g0, alpha=ALPHA, g_off=G_OFF):
    r, g, _ = integrate(g0, alpha, g_off)
    Av, rv = [], []
    for rL in WIN_LIST:
        if rL + WIN_WIDTH > r[-1]: break
        A, _, _ = fit_win(r, g, rL, rL+WIN_WIDTH)
        if A > 0.002: Av.append(A); rv.append(float(rL))
    if not Av: return 0.
    if len(Av) < 3: return Av[-1]
    try:
        p, _ = curve_fit(lambda x, ai, a: ai*(1+a/x), rv, Av, p0=[Av[-1], 0.], maxfev=2000)
        return float(p[0])
    except:
        return float(Av[-1])

# ============================================================
# OBLICZANIE M_total(R) ze ZSZYWANEGO rozwiazania
# ============================================================
def f_tgp(g):
    return 1.0 + 2.0*ALPHA*np.log(np.clip(g, 1e-15, None))

def V_pot(g):
    return g**3/3.0 - g**4/4.0 - 1.0/12.0

def M_total_R(g0, R_max, n_pts=6000):
    """Calkowita energia do promienia R_max (ze zszywanym rozwiazaniem)."""
    r_arr, g_arr, gp_arr = integrate(g0, r_end=R_max+5)
    # Interpolacja na rownomierna siatke
    r_eval = np.linspace(1e-3, R_max, n_pts)
    # Odcinamy poza zasiegiem
    r_eval = r_eval[r_eval <= r_arr[-1] - 0.01]
    if len(r_eval) < 100:
        return np.nan
    # Interpolacja g i gp
    g_eval  = np.interp(r_eval, r_arr, g_arr)
    gp_eval = np.interp(r_eval, r_arr, gp_arr)
    fg = f_tgp(g_eval)
    dens = fg * gp_eval**2 / 2.0 + V_pot(g_eval)
    return float(np.trapezoid(r_eval**2 * dens, r_eval) * 4.0 * np.pi)

# ============================================================
# WSTEPNE DANE
# ============================================================
print("=" * 68)
print("EX68: MASA RENORMALIZOWANA M_ren = M_total(R) - pi*A_inf^2*R")
print("=" * 68)

A_e  = A_inf_val(G0_STAR)
A_mu = A_inf_val(PHI * G0_STAR)

print(f"\n  g0_e   = {G0_STAR:.7f},   A_e  = {A_e:.7f}")
print(f"  g0_mu  = {PHI*G0_STAR:.7f},   A_mu = {A_mu:.7f}")
print(f"  (A_mu/A_e)^4 = {(A_mu/A_e)**4:.4f}  (r21_exp = {R21_EXP})")

# ============================================================
# CZESC 1: M_total(R) vs R dla elektronu
# ============================================================
print("\n" + "-" * 68)
print("--- Czesc 1: M_total(R) vs R dla elektronu (g0*) ---")
print()

R_list = [20, 40, 60, 80, 100, 130, 160, 200]
M_e_list  = []
M_mu_list = []

print(f"  {'R_max':>7}  {'M_tot_e':>11}  {'M/pi/A2/R':>10}  {'M_ren_e':>11}")
print(f"  {'-'*7}  {'-'*11}  {'-'*10}  {'-'*11}")
M_ren_e_list = []
for R in R_list:
    M = M_total_R(G0_STAR, R)
    M_e_list.append(M)
    coeff = M / (np.pi * A_e**2 * R) if (A_e > 0 and R > 0) else np.nan
    Mren  = M - np.pi * A_e**2 * R
    M_ren_e_list.append(Mren)
    print(f"  {R:7d}  {M:11.5f}  {coeff:10.5f}  {Mren:11.5f}")

# ============================================================
# CZESC 2: M_total(R) vs R dla miona
# ============================================================
print("\n" + "-" * 68)
print("--- Czesc 2: M_total(R) vs R dla miona (phi*g0*) ---")
print()

print(f"  {'R_max':>7}  {'M_tot_mu':>12}  {'M/pi/A2/R':>10}  {'M_ren_mu':>12}")
print(f"  {'-'*7}  {'-'*12}  {'-'*10}  {'-'*12}")
M_ren_mu_list = []
for R in R_list:
    M = M_total_R(PHI * G0_STAR, R)
    M_mu_list.append(M)
    coeff = M / (np.pi * A_mu**2 * R) if (A_mu > 0 and R > 0) else np.nan
    Mren  = M - np.pi * A_mu**2 * R
    M_ren_mu_list.append(Mren)
    print(f"  {R:7d}  {M:12.5f}  {coeff:10.5f}  {Mren:12.5f}")

# ============================================================
# CZESC 3: Ekstrapolacja M_ren(R->inf) z fitem oscylujacym
# ============================================================
print("\n" + "-" * 68)
print("--- Czesc 3: Ekstrapolacja M_ren(R->inf) ---")
print()

def fit_Mren(R_arr, Mren_arr):
    """Fit M_ren(R) = M_inf + b/R + c*cos(2R)/R + d*sin(2R)/R"""
    R = np.asarray(R_arr, dtype=float)
    A_mat = np.column_stack([
        np.ones_like(R),
        1.0 / R,
        np.cos(2.0*R) / R,
        np.sin(2.0*R) / R
    ])
    coef, _, _, _ = lstsq(A_mat, Mren_arr, rcond=None)
    return float(coef[0]), coef

# Uzywamy 5 ostatnich punktow (wieksze R)
R_big   = np.array(R_list[-5:], dtype=float)
Mren_e  = np.array(M_ren_e_list[-5:])
Mren_mu = np.array(M_ren_mu_list[-5:])

M_inf_e,  coef_e  = fit_Mren(R_big, Mren_e)
M_inf_mu, coef_mu = fit_Mren(R_big, Mren_mu)

print(f"  Elektron g0* = {G0_STAR:.7f}:")
print(f"    M_ren (R=20..200): {[f'{x:.4f}' for x in M_ren_e_list]}")
print(f"    M_ren_inf (fit):   {M_inf_e:.6f}")

print(f"\n  Mion phi*g0* = {PHI*G0_STAR:.7f}:")
print(f"    M_ren (R=20..200): {[f'{x:.4f}' for x in M_ren_mu_list]}")
print(f"    M_ren_inf (fit):   {M_inf_mu:.6f}")

ratio_Mren = M_inf_mu / M_inf_e if abs(M_inf_e) > 1e-8 else np.nan
print(f"\n  Stosunek M_ren_inf(mu/e) = {ratio_Mren:.4f}  (r21 = {R21_EXP})")
print(f"  A^4 ratio               = {(A_mu/A_e)**4:.4f}")

# ============================================================
# CZESC 4: M_ren vs A_inf dla calego sektora n=0
# ============================================================
print("\n" + "-" * 68)
print("--- Czesc 4: M_ren/A^4 w sektorze n=0 (R=100) ---")
print()

g0_scan = np.linspace(1.10, 1.58, 16)
R_fixed = 100

Mren_scan = []
A_scan    = []

print(f"  {'g0':>8}  {'A_inf':>10}  {'M_tot':>11}  {'M_ren':>11}  {'Mren/A^4':>11}")
print(f"  {'-'*8}  {'-'*10}  {'-'*11}  {'-'*11}  {'-'*11}")

for g0 in g0_scan:
    try:
        Av = A_inf_val(g0)
        Mt = M_total_R(g0, R_fixed)
        if np.isnan(Mt): continue
        Mr = Mt - np.pi * Av**2 * R_fixed
        Mren_scan.append(Mr)
        A_scan.append(Av)
        c4 = Mr / Av**4 if Av > 1e-6 else np.nan
        marker = " <g0*" if abs(g0 - G0_STAR) < 0.025 else ""
        print(f"  {g0:8.4f}  {Av:10.6f}  {Mt:11.5f}  {Mr:11.5f}  {c4:11.3f}{marker}")
    except Exception as ex:
        print(f"  {g0:8.4f}  BLAD: {ex}")

A_scan    = np.array(A_scan)
Mren_scan = np.array(Mren_scan)

# Power-law fit: log(|Mren|) ~ k*log(A)
valid_pos = np.isfinite(Mren_scan) & (Mren_scan > 0) & (A_scan > 1e-6)
if valid_pos.sum() >= 4:
    sl, ic, r2, _, _ = linregress(np.log(A_scan[valid_pos]), np.log(Mren_scan[valid_pos]))
    k_fit, R2_fit = sl, r2**2
    print(f"\n  Fit M_ren~A^k (M_ren>0): k = {k_fit:.3f}, R^2 = {R2_fit:.5f}")
else:
    valid_all = np.isfinite(Mren_scan) & (A_scan > 1e-6)
    if valid_all.sum() >= 4:
        sl2, _, r2_2, _, _ = linregress(np.log(A_scan[valid_all]), np.log(np.abs(Mren_scan[valid_all])))
        k_fit, R2_fit = sl2, r2_2**2
        print(f"\n  Fit |M_ren|~A^k (wszystkie): k = {k_fit:.3f}, R^2 = {R2_fit:.5f}")
    else:
        k_fit, R2_fit = np.nan, np.nan

# CV
with np.errstate(invalid='ignore', divide='ignore'):
    c4_arr = Mren_scan / A_scan**4
valid_cv = np.isfinite(c4_arr)
if valid_cv.sum() >= 4:
    mean_c4 = np.mean(c4_arr[valid_cv])
    std_c4  = np.std(c4_arr[valid_cv])
    cv_val  = abs(std_c4 / mean_c4) * 100 if mean_c4 != 0 else np.nan
    print(f"  M_ren/A^4: mean = {mean_c4:.2f}, std = {std_c4:.2f}, CV = {cv_val:.1f}%")
else:
    cv_val = np.nan

# ============================================================
# CZESC 5: Zbieznosc M_ren (dodatkowy test)
# ============================================================
print("\n" + "-" * 68)
print("--- Czesc 5: Zbieznosc M_ren(R) dla g0* ---")
print()

R_test = [60, 80, 100, 130, 160, 200]
Mren_cv_e  = []
Mren_cv_mu = []
print(f"  {'R_max':>7}  {'M_ren_e':>11}  {'M_ren_mu':>12}  {'ratio':>8}")
print(f"  {'-'*7}  {'-'*11}  {'-'*12}  {'-'*8}")
for R in R_test:
    Me  = M_total_R(G0_STAR,       R) - np.pi*A_e**2 *R
    Mmu = M_total_R(PHI*G0_STAR,   R) - np.pi*A_mu**2*R
    Mren_cv_e.append(Me)
    Mren_cv_mu.append(Mmu)
    ratio_r = Mmu/Me if abs(Me) > 1e-8 else np.nan
    ratio_str = f"{ratio_r:.3f}" if not np.isnan(ratio_r) else "N/A"
    print(f"  {R:7d}  {Me:11.5f}  {Mmu:12.5f}  {ratio_str:>8}")

Mren_cv_e  = np.array(Mren_cv_e)
Mren_cv_mu = np.array(Mren_cv_mu)
rng_e   = np.ptp(Mren_cv_e)
mean_e  = abs(np.mean(Mren_cv_e))
rel_rng = rng_e / mean_e * 100 if mean_e > 1e-10 else np.nan
print(f"\n  Zakres M_ren_e: {rng_e:.5f}  (wzgledny: {rel_rng:.1f}%)")

M_inf_e3,  _ = fit_Mren(np.array(R_test[-4:], dtype=float), Mren_cv_e[-4:])
M_inf_mu3, _ = fit_Mren(np.array(R_test[-4:], dtype=float), Mren_cv_mu[-4:])
ratio_inf3 = M_inf_mu3 / M_inf_e3 if abs(M_inf_e3) > 1e-8 else np.nan
ratio_inf3_str = f"{ratio_inf3:.4f}" if not np.isnan(ratio_inf3) else "N/A"
print(f"  Ekstrapolacja: M_ren_e_inf = {M_inf_e3:.5f}, M_ren_mu_inf = {M_inf_mu3:.5f}")
print(f"  Stosunek M_ren_inf(mu/e) = {ratio_inf3_str}")

# ============================================================
# TESTY
# ============================================================
print("\n" + "=" * 68)
print("TESTY")
print("=" * 68)

# T1: koefficient M_tot/pi/A^2/R -> 1 dla duzych R
coeff_e_big = []
for R, M in zip(R_list[-3:], M_e_list[-3:]):
    if A_e > 0 and R > 0 and not np.isnan(M):
        coeff_e_big.append(M / (np.pi * A_e**2 * R))
T1 = len(coeff_e_big) >= 2 and all(abs(c-1.0) < 0.30 for c in coeff_e_big)
coeff_str = ', '.join(f'{c:.3f}' for c in coeff_e_big)
print(f"  T1: M_tot ~ pi*A^2*R (|c-1|<0.30):  {'PASS' if T1 else 'FAIL'}  "
      f"(c = {coeff_str})")

# T2: M_ren_e zbiega (wzgledny zakres < 50%)
T2 = not np.isnan(rel_rng) and rel_rng < 50.0
print(f"  T2: M_ren_e zbiega (zakres<50%):    {'PASS' if T2 else 'FAIL'}  "
      f"(zakres/mean = {rel_rng:.1f}%)")

# T3: k in [3.5, 4.5]
T3 = not np.isnan(k_fit) and 3.5 <= k_fit <= 4.5
k_str  = f"{k_fit:.3f}"  if not np.isnan(k_fit)  else "N/A"
r2_str = f"{R2_fit:.4f}" if not np.isnan(R2_fit) else "N/A"
print(f"  T3: M_ren ~ A^k, k in [3.5,4.5]:    {'PASS' if T3 else 'FAIL'}  "
      f"(k={k_str}, R^2={r2_str})")

# T4: |stosunek - 207| < 20%
T4 = not np.isnan(ratio_inf3) and abs(ratio_inf3 - R21_EXP) / R21_EXP < 0.20
print(f"  T4: |M_ren_ratio - 207| < 20%:      {'PASS' if T4 else 'FAIL'}  "
      f"(ratio={ratio_inf3_str})")

# T5: CV < 40%
T5 = not np.isnan(cv_val) and cv_val < 40.0
cv_str = f"{cv_val:.1f}%" if not np.isnan(cv_val) else "N/A"
print(f"  T5: M_ren/A^4 CV < 40%:              {'PASS' if T5 else 'FAIL'}  "
      f"(CV={cv_str})")

n_pass = sum([T1, T2, T3, T4, T5])
print(f"\nWYNIK: {n_pass}/5 testow przeszlo")

print("\n" + "=" * 68)
print("WNIOSEK (EX68)")
print("=" * 68)
print(f"  A_e   = {A_e:.7f}, A_mu = {A_mu:.7f}")
print(f"  A^4 ratio = {(A_mu/A_e)**4:.4f}  (r21 = {R21_EXP})")
print(f"  M_ren_inf_e   = {M_inf_e3:.5f}")
print(f"  M_ren_inf_mu  = {M_inf_mu3:.5f}")
print(f"  M_ren stosunek = {ratio_inf3_str}")
print("=" * 68)
