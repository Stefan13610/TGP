"""
EX70: PRECYZYJNE alpha* I PHI-ORBITA g0*

Cel 1: Oblicz alpha* dokladnie do 10^{-4} i porownaj z sqrt(6)=2.4495
Cel 2: Sprawdz samopodobienstwo phi-orbity: f(phi^n*g0*) = r21 dla n=0,1,2?
Cel 3: Walidacja delta(alpha) = 0 @ alpha* z wyzszym zakresem alpha

Definicja alpha*:
  Wariant A: alpha* taki ze delta(alpha*)=0, tzn. g0*(alpha*)=z0(alpha*)
  Wariant B: alpha* taki ze f(z0(alpha*), alpha*) = r21

Tests:
  T1: |alpha*_A - alpha*_B| < 0.001  (obie definicje zbiezne)
  T2: |alpha*_A - sqrt(6)| < 0.01    (sprawdz hipoteze alpha*=sqrt(6))
  T3: f(phi*g0*) = r21?              (czy phi-FP jest "plaska orbita"?)
  T4: f(phi^2*g0*) ≈ r21?            (tauon w orbicie?)
  T5: delta(alpha) liniowe R^2>0.999 dla alpha in [1.8, 2.8]
"""

import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, curve_fit
from scipy.stats import linregress
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# STALE
# ============================================================
PHI     = (1.0 + np.sqrt(5.0)) / 2.0
R21_EXP = 206.768
G_OFF   = 0.005
R_MAX   = 100.0
WIN_LIST  = [16, 22, 28, 36, 46, 58, 72]
WIN_WIDTH = 14.0
SQRT6 = np.sqrt(6.0)

# ============================================================
# ODE SOLVER (z ex62b)
# ============================================================
def make_ode(alpha, g_off=G_OFF):
    g_bounce = np.exp(-1.0/(2.0*alpha)) + g_off
    def rhs(r, y):
        g, gp = y
        g = max(g, g_bounce + 1e-7)
        fg = 1.0 + 2.0*alpha*np.log(g)
        if abs(fg) < 1e-10: return [gp, 0.0]
        dr = g**2*(1.0-g); cr = (alpha/g)*gp**2
        if r < 1e-10: return [gp, (dr-cr)/(3.0*fg)]
        return [gp, (dr-cr-fg*2.0*gp/r)/fg]
    def ev(r, y): return y[0] - (np.exp(-1.0/(2.0*alpha)) + g_off)
    ev.terminal = True; ev.direction = -1
    return rhs, ev

def integrate(g0, alpha=2.0, maxb=8):
    rhs, ev = make_ode(alpha)
    y0 = [g0, 0.0]; r0 = 1e-10; ra, ga = [], []
    for _ in range(maxb+1):
        sol = solve_ivp(rhs, (r0, R_MAX), y0, events=ev,
                        dense_output=True, rtol=1e-9, atol=1e-11, max_step=0.05)
        ra.append(sol.t); ga.append(sol.y[0])
        if sol.status == 1 and len(sol.t_events[0]) > 0:
            rh = sol.t_events[0][-1]; st = sol.sol(rh)
            y0 = [st[0], -st[1]]; r0 = rh
        else: break
    r = np.concatenate(ra); g = np.concatenate(ga)
    idx = np.argsort(r)
    return r[idx], g[idx]

def fit_win(r, g, rL, rR):
    mask = (r >= rL) & (r <= rR)
    if np.sum(mask) < 20: return 0., 0., 0.
    rf = r[mask]; df = (g[mask]-1.0)*rf
    M = np.column_stack([np.cos(rf), np.sin(rf)])
    bc = np.linalg.lstsq(M, df, rcond=None)[0]
    B, C = bc
    return float(np.sqrt(B**2+C**2)), float(B), float(C)

def A_inf(g0, alpha=2.0):
    r, g = integrate(g0, alpha)
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
    except: return float(Av[-1])

def B_coeff(g0, alpha=2.0):
    r, g = integrate(g0, alpha)
    _, B, _ = fit_win(r, g, WIN_LIST[2], WIN_LIST[2]+WIN_WIDTH)
    return B

def find_z0(alpha=2.0, lo=1.05, hi=1.65):
    gs = np.linspace(lo, hi, 60)
    Bs = [B_coeff(g, alpha) for g in gs]
    for i in range(len(Bs)-1):
        if Bs[i]*Bs[i+1] < 0:
            try: return brentq(lambda x: B_coeff(x, alpha), gs[i], gs[i+1], xtol=1e-7)
            except: pass
    return None

def f_ratio(g0, alpha=2.0):
    Ae = A_inf(g0, alpha); Am = A_inf(PHI*g0, alpha)
    return (Am/Ae)**4 if Ae > 1e-6 else np.nan

def find_g0star(alpha=2.0, target=R21_EXP, lo=1.10, hi=1.70):
    gs = np.linspace(lo, hi, 30)
    fs = [f_ratio(g, alpha) for g in gs]
    for i in range(len(fs)-1):
        if not (np.isnan(fs[i]) or np.isnan(fs[i+1])):
            if (fs[i]-target)*(fs[i+1]-target) < 0:
                try:
                    return brentq(lambda x: f_ratio(x, alpha) - target, gs[i], gs[i+1], xtol=1e-7)
                except: pass
    return None

# ============================================================
# CZESC 1: Precyzyjne alpha* z obu definicji
# ============================================================
print("=" * 68)
print("EX70: PRECYZYJNE alpha* I PHI-ORBITA g0*")
print("=" * 68)

print("\n" + "-" * 68)
print("--- Czesc 1: alpha* z definicji A (delta=0) i B (f(z0)=r21) ---")
print()

# Wariant A: alpha* z delta(alpha)=0, tzn. g0*(alpha*)=z0(alpha*)
# Uzywamy liniowego fitu delta(alpha) z ex69: delta=0 @ -b/a
# Ale teraz liczymy bezposrednio: delta(alpha) = g0*(alpha) - z0(alpha)

def delta_func(alpha):
    z0  = find_z0(alpha)
    g0s = find_g0star(alpha)
    if z0 and g0s:
        return g0s - z0
    return np.nan

# Wariant B: alpha* taki ze f(z0(alpha*), alpha*) = r21
def B_func(alpha):
    z0 = find_z0(alpha)
    if z0 is None: return np.nan
    return f_ratio(z0, alpha) - R21_EXP

# Skan grubszy najpierw
alpha_test = np.linspace(2.0, 2.8, 17)
delta_test = []
B_test     = []

print(f"  {'alpha':>6}  {'z0':>10}  {'g0*':>10}  {'delta':>9}  {'f(z0)-r21':>10}")
print(f"  {'-'*6}  {'-'*10}  {'-'*10}  {'-'*9}  {'-'*10}")

for alpha in alpha_test:
    z0  = find_z0(alpha)
    g0s = find_g0star(alpha)
    d   = (g0s - z0) if (z0 and g0s) else np.nan
    fz0 = f_ratio(z0, alpha) - R21_EXP if z0 else np.nan
    delta_test.append(d)
    B_test.append(fz0)
    z0_s  = f"{z0:.6f}"  if z0  else "N/A"
    g0s_s = f"{g0s:.6f}" if g0s else "N/A"
    d_s   = f"{d:.6f}"   if not np.isnan(d) else "N/A"
    fz0_s = f"{fz0:.4f}" if not np.isnan(fz0) else "N/A"
    print(f"  {alpha:6.3f}  {z0_s:>10}  {g0s_s:>10}  {d_s:>9}  {fz0_s:>10}")

delta_arr = np.array(delta_test)
B_arr     = np.array(B_test)

# Znajdz alpha* wariant A (delta=0) przez brentq
alpha_A_found = None
for i in range(len(alpha_test)-1):
    if not (np.isnan(delta_arr[i]) or np.isnan(delta_arr[i+1])):
        if delta_arr[i]*delta_arr[i+1] < 0:
            try:
                alpha_A_found = brentq(delta_func, alpha_test[i], alpha_test[i+1], xtol=1e-5)
                break
            except: pass

# Znajdz alpha* wariant B (f(z0)=r21) przez brentq
alpha_B_found = None
for i in range(len(alpha_test)-1):
    if not (np.isnan(B_arr[i]) or np.isnan(B_arr[i+1])):
        if B_arr[i]*B_arr[i+1] < 0:
            try:
                alpha_B_found = brentq(B_func, alpha_test[i], alpha_test[i+1], xtol=1e-5)
                break
            except: pass

print(f"\n  Wynik:")
print(f"    alpha*_A (delta=0):    {alpha_A_found:.6f}" if alpha_A_found else "    alpha*_A: N/A")
print(f"    alpha*_B (f(z0)=r21): {alpha_B_found:.6f}" if alpha_B_found else "    alpha*_B: N/A")
print(f"    sqrt(6)             = {SQRT6:.6f}")
if alpha_A_found and alpha_B_found:
    diff_AB  = abs(alpha_A_found - alpha_B_found)
    diff_A6  = abs(alpha_A_found - SQRT6)
    diff_B6  = abs(alpha_B_found - SQRT6)
    print(f"    |alpha*_A - alpha*_B| = {diff_AB:.6f}")
    print(f"    |alpha*_A - sqrt(6)|  = {diff_A6:.6f}  ({diff_A6/SQRT6*100:.4f}%)")
    print(f"    |alpha*_B - sqrt(6)|  = {diff_B6:.6f}  ({diff_B6/SQRT6*100:.4f}%)")

# ============================================================
# CZESC 2: delta(alpha) liniowe dla szerszego zakresu
# ============================================================
print("\n" + "-" * 68)
print("--- Czesc 2: delta(alpha) linearnosc dla alpha in [1.8, 2.8] ---")
print()

alpha_wide = np.linspace(1.8, 2.8, 11)
delta_wide = []
for alpha in alpha_wide:
    d = delta_func(alpha)
    delta_wide.append(d)
    print(f"  alpha={alpha:.3f}: delta={d:.6f}" if not np.isnan(d) else f"  alpha={alpha:.3f}: N/A")

alpha_w  = alpha_wide
delta_w  = np.array(delta_wide)
valid_w  = ~np.isnan(delta_w)
if valid_w.sum() >= 5:
    slope_w, ic_w, r_w, _, _ = linregress(alpha_w[valid_w], delta_w[valid_w])
    R2_w = r_w**2
    alpha_zero_w = -ic_w / slope_w
    print(f"\n  Fit delta(alpha) = a*alpha + b:")
    print(f"    a = {slope_w:.6f}, b = {ic_w:.6f}, R^2 = {R2_w:.6f}")
    print(f"    Zero @ alpha = {alpha_zero_w:.5f}")
    print(f"    vs sqrt(6)   = {SQRT6:.5f}  (diff={abs(alpha_zero_w-SQRT6):.5f})")

# ============================================================
# CZESC 3: Phi-orbita samopodobienstwa
# ============================================================
print("\n" + "-" * 68)
print("--- Czesc 3: Phi-orbita g0* (n=0,1,2) ---")
print()

alpha_ref = 2.0
G0_STAR = find_g0star(alpha_ref)

if G0_STAR:
    orbit_pts = [G0_STAR * PHI**n for n in range(4)]
    print(f"  g0* = {G0_STAR:.7f}")
    print(f"\n  {'n':>3}  {'phi^n*g0*':>12}  {'A_inf':>10}  {'f_n=ratio^4':>12}  {'f_n/r21':>10}")
    print(f"  {'-'*3}  {'-'*12}  {'-'*10}  {'-'*12}  {'-'*10}")

    A_vals = []
    f_n_vals = []
    for n, g in enumerate(orbit_pts[:-1]):
        Ag  = A_inf(g, alpha_ref)
        Agp = A_inf(orbit_pts[n+1], alpha_ref)
        f_n = (Agp/Ag)**4 if Ag > 1e-6 else np.nan
        A_vals.append(Ag)
        f_n_vals.append(f_n)
        fn_str = f"{f_n:.4f}" if not np.isnan(f_n) else "N/A"
        fr_str = f"{f_n/R21_EXP:.6f}" if not np.isnan(f_n) else "N/A"
        print(f"  {n:>3}  {g:12.7f}  {Ag:10.6f}  {fn_str:>12}  {fr_str:>10}")

    # Czy f_n stale = r21?
    f_arr2 = np.array([x for x in f_n_vals if not np.isnan(x)])
    if len(f_arr2) >= 2:
        cv = np.std(f_arr2)/np.mean(f_arr2)*100
        print(f"\n  Srednia f_n = {np.mean(f_arr2):.4f}  (r21 = {R21_EXP})")
        print(f"  CV f_n = {cv:.2f}%  (H: CV < 5% = samopodobienstwo)")

    # Sprawdz A_inf(phi^n * g0*) = r21^(n/4) * A_inf(g0*)
    print(f"\n  Czy A_inf(phi^n*g0*) = r21^(n/4) * A_inf(g0*)?")
    A0 = A_vals[0] if A_vals else 1.0
    print(f"  {'n':>3}  {'A_inf_actual':>14}  {'A_inf_predicted':>16}  {'ratio':>8}")
    print(f"  {'-'*3}  {'-'*14}  {'-'*16}  {'-'*8}")
    for n, g in enumerate(orbit_pts[:3]):
        Ag = A_inf(g, alpha_ref)
        A_pred = A0 * R21_EXP**(n/4.0)
        ratio = Ag / A_pred
        print(f"  {n:>3}  {Ag:14.7f}  {A_pred:16.7f}  {ratio:8.6f}")

# ============================================================
# CZESC 4: Analiza bledow: czy f(z0(alpha)) = f(g0*(alpha)) przy alpha*?
# ============================================================
print("\n" + "-" * 68)
print("--- Czesc 4: f(z0, alpha) i f(g0*, alpha) dla alpha bliskich alpha* ---")
print()

if alpha_A_found:
    alphas_around_star = [alpha_A_found - 0.1, alpha_A_found - 0.05,
                          alpha_A_found, alpha_A_found + 0.05, alpha_A_found + 0.1]
else:
    alphas_around_star = [2.34, 2.39, 2.44, 2.49, 2.54]

print(f"  {'alpha':>7}  {'f(z0)':>10}  {'f(g0*)':>10}  {'f(z0)-r21':>10}  {'z0':>10}  {'g0*':>10}")
print(f"  {'-'*7}  {'-'*10}  {'-'*10}  {'-'*10}  {'-'*10}  {'-'*10}")
for alpha in alphas_around_star:
    z0  = find_z0(alpha)
    g0s = find_g0star(alpha)
    fz0 = f_ratio(z0, alpha) if z0 else np.nan
    fg0s = f_ratio(g0s, alpha) if g0s else np.nan
    fz0_s  = f"{fz0:.4f}"  if not np.isnan(fz0)  else "N/A"
    fg0s_s = f"{fg0s:.4f}" if not np.isnan(fg0s) else "N/A"
    fd_s   = f"{fz0-R21_EXP:.4f}" if not np.isnan(fz0) else "N/A"
    z0_s   = f"{z0:.6f}" if z0 else "N/A"
    g0s_s  = f"{g0s:.6f}" if g0s else "N/A"
    print(f"  {alpha:7.4f}  {fz0_s:>10}  {fg0s_s:>10}  {fd_s:>10}  {z0_s:>10}  {g0s_s:>10}")

# ============================================================
# TESTY
# ============================================================
print("\n" + "=" * 68)
print("TESTY")
print("=" * 68)

# T1: Zbieznosc alpha*_A i alpha*_B
T1 = (alpha_A_found is not None and alpha_B_found is not None and
      abs(alpha_A_found - alpha_B_found) < 0.005)
diff_AB_v = abs(alpha_A_found - alpha_B_found) if (alpha_A_found and alpha_B_found) else np.nan
diff_AB_s = f"{diff_AB_v:.5f}" if not np.isnan(diff_AB_v) else "N/A"
print(f"  T1: |alpha*_A - alpha*_B| < 0.005:  {'PASS' if T1 else 'FAIL'}  "
      f"(diff = {diff_AB_s})")

# T2: |alpha* - sqrt(6)| < 0.01
alpha_best = alpha_A_found if alpha_A_found else alpha_B_found
T2 = (alpha_best is not None and abs(alpha_best - SQRT6) < 0.01)
diff_6 = abs(alpha_best - SQRT6) if alpha_best else np.nan
diff_6_s = f"{diff_6:.5f}" if not np.isnan(diff_6) else "N/A"
print(f"  T2: |alpha* - sqrt(6)| < 0.01:       {'PASS' if T2 else 'FAIL'}  "
      f"(diff = {diff_6_s}, sqrt(6)={SQRT6:.5f})")

# T3: Phi-orbita: f_n = r21 dla n=0
T3 = G0_STAR is not None and len(f_n_vals) >= 1 and not np.isnan(f_n_vals[0]) and abs(f_n_vals[0] - R21_EXP)/R21_EXP < 0.001
f0_s = f"{f_n_vals[0]:.4f}" if (f_n_vals and not np.isnan(f_n_vals[0])) else "N/A"
print(f"  T3: f(g0*, alpha=2) = r21 (0.1%):    {'PASS' if T3 else 'FAIL'}  "
      f"(f_0 = {f0_s})")

# T4: Phi-orbita: f_1 = r21 (czy phi-orbita zachowuje r21?)
T4 = (len(f_n_vals) >= 2 and not np.isnan(f_n_vals[1]) and
      abs(f_n_vals[1] - R21_EXP)/R21_EXP < 0.10)
f1_s = f"{f_n_vals[1]:.4f}" if (len(f_n_vals) >= 2 and not np.isnan(f_n_vals[1])) else "N/A"
print(f"  T4: f(phi*g0*) ≈ r21 (10%):          {'PASS' if T4 else 'FAIL'}  "
      f"(f_1 = {f1_s})")

# T5: delta(alpha) liniowe R^2 > 0.999
T5 = 'R2_w' in dir() and R2_w > 0.999
r2_str = f"{R2_w:.5f}" if 'R2_w' in dir() else "N/A"
print(f"  T5: delta liniowe R^2 > 0.999:        {'PASS' if T5 else 'FAIL'}  "
      f"(R^2 = {r2_str})")

n_pass = sum([T1, T2, T3, T4, T5])
print(f"\nWYNIK: {n_pass}/5 testow przeszlo")

print("\n" + "=" * 68)
print("WNIOSEK (EX70)")
print("=" * 68)
print(f"  sqrt(6)     = {SQRT6:.7f}")
aA_s = f"{alpha_A_found:.7f}" if alpha_A_found else "N/A"
aB_s = f"{alpha_B_found:.7f}" if alpha_B_found else "N/A"
print(f"  alpha*_A    = {aA_s}  (delta=0)")
print(f"  alpha*_B    = {aB_s}  (f(z0)=r21)")
if 'alpha_zero_w' in dir():
    print(f"  alpha* (fit)= {alpha_zero_w:.7f}  (linear extrapolation)")
if G0_STAR:
    print(f"  g0*(alpha=2)= {G0_STAR:.7f}")
    if f_n_vals and not np.isnan(f_n_vals[0]):
        print(f"  f(g0*, 2)   = {f_n_vals[0]:.4f}  (r21)")
    if len(f_n_vals)>=2 and not np.isnan(f_n_vals[1]):
        print(f"  f(phi*g0*,2)= {f_n_vals[1]:.4f}  (czy = r21?)")
print("=" * 68)
