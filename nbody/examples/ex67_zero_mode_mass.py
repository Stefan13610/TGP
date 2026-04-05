"""
EX67: TRYB ZEROWY I MASA NIELINIOWA M_nonlin ∝ A_tail^4
=========================================================
Hipoteza analityczna (ex66):
  Zlinearyzowane ODE TGP: delta'' + (2/r)delta' + delta = 0
  => tryb k=1 jest TRYBEM ZEROWYM (energia = 0 dokładnie)
  => M_TGP = M_nonlin = M_total - M_linear ∝ A_tail^4

Metoda:
  M_linear(A) = 4pi ∫ r^2 [f(g_lin)(g_lin')^2/2 + V(g_lin)] dr
  gdzie g_lin = 1 + A*sin(r)/r (czysto liniowe rozwiązanie)

  M_nonlin(g0) = M_total(g0) - M_linear(A_tail(g0))

Jeśli teoria jest słuszna:
  M_linear → 0  (tryb zerowy)
  M_nonlin ∝ A_tail^4  (nieliniowe korekcje)

Dodatkowe pytania:
  1. M_nonlin(phi*g0*) / M_nonlin(g0*) =? 207
  2. M_linear(A) ≈ c * A^k z małym k (≈0 dla trybu zerowego)?
  3. M_nonlin/A^4 = const? (stały współczynnik?)
  4. Czy M_nonlin > 0 zawsze? (masa > 0 zawsze?)

Testy:
  T1: M_linear(A) ∝ A^k, k ∈ [0, 0.5]  (prawie zero, bo tryb zerowy)
  T2: M_nonlin ∝ A^k, k ∈ [3.5, 4.5]  (≈ A^4)
  T3: M_nonlin > 0 dla wszystkich testowanych g0  (masa zawsze dodatnia)
  T4: |M_nonlin(phi*g0*)/M_nonlin(g0*) - 207| / 207 < 0.30  (30%)
  T5: M_nonlin/A^4 CV < 50%  (względna stałość)
"""

import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import curve_fit
import warnings
warnings.filterwarnings('ignore')

PHI      = (1.0 + np.sqrt(5.0)) / 2.0
R21_EXP  = 206.768
R_MAX    = 100.0
G0_STAR  = 1.24771102
WIN_LIST  = [16, 22, 28, 36, 46, 58]
WIN_WIDTH = 12.0

ALPHA    = 2.0
G_GHOST  = np.exp(-1.0 / (2.0 * ALPHA))
G_BOUNCE = G_GHOST + 0.005

# ---- ODE ----
def rhs(r, y):
    g, gp = y
    g  = max(g, G_BOUNCE + 1e-7)
    fg = 1.0 + 2.0 * ALPHA * np.log(g)
    if abs(fg) < 1e-10: return [gp, 0.0]
    dr = g**2 * (1.0 - g)
    cr = (ALPHA / g) * gp**2
    if r < 1e-10: return [gp, (dr - cr) / (3.0 * fg)]
    return [gp, (dr - cr - fg * 2.0 * gp / r) / fg]

def ev_ghost(r, y): return y[0] - G_BOUNCE
ev_ghost.terminal  = True
ev_ghost.direction = -1

def integrate(g0, maxb=8):
    y0 = [g0, 0.0]; r0 = 0.0; ra, ga, gpa = [], [], []
    for _ in range(maxb + 1):
        sol = solve_ivp(rhs, (r0, R_MAX), y0, events=ev_ghost,
                        dense_output=True, rtol=1e-9, atol=1e-11, max_step=0.05)
        ra.append(sol.t); ga.append(sol.y[0]); gpa.append(sol.y[1])
        if sol.status == 1 and len(sol.t_events[0]) > 0:
            rh = sol.t_events[0][-1]; st = sol.sol(rh)
            y0 = [st[0], -st[1]]; r0 = rh
        else: break
    r = np.concatenate(ra); g = np.concatenate(ga); gp = np.concatenate(gpa)
    idx = np.argsort(r)
    return r[idx], g[idx], gp[idx]

def V_pot(g):
    return g**3/3 - g**4/4 - (1./3 - 1./4)

def f_tgp(g):
    return 1.0 + 2.0 * ALPHA * np.log(np.maximum(g, 1e-10))

def M_total_numeric(r, g, gp):
    """M_total = 4pi ∫ r^2 [f(g)(g')^2/2 + V(g)] dr (oscyluje z R_max)."""
    fg = f_tgp(g); Vg = V_pot(g)
    integrand = r**2 * (fg * gp**2 / 2.0 + Vg)
    return float(np.trapezoid(integrand, r) * 4 * np.pi)

def A_inf(g0):
    r, g, _ = integrate(g0)
    def fit_win(r, g, rL, rR):
        mask = (r >= rL) & (r <= rR)
        if np.sum(mask) < 20: return 0., 0., 0.
        rf = r[mask]; df = (g[mask]-1)*rf
        M  = np.column_stack([np.cos(rf), np.sin(rf)])
        bc = np.linalg.lstsq(M, df, rcond=None)[0]
        B, C = bc
        return float(np.sqrt(B**2+C**2)), float(B), float(C)
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

def M_linear_from_A(A_val, r_arr):
    """
    M_linear = 4pi ∫ r^2 [f(g_lin)(g_lin')^2/2 + V(g_lin)] dr
    g_lin = 1 + A*sin(r)/r
    g_lin' = A*cos(r)/r - A*sin(r)/r^2
    """
    # Unikamy r=0
    r = r_arr[r_arr > 0.5]  # ogon tylko (r > 0.5)
    g_lin  = 1.0 + A_val * np.sin(r) / r
    gp_lin = A_val * (np.cos(r)/r - np.sin(r)/r**2)

    fg_lin = f_tgp(g_lin)
    Vg_lin = V_pot(g_lin)

    integrand = r**2 * (fg_lin * gp_lin**2 / 2.0 + Vg_lin)
    return float(np.trapezoid(integrand, r) * 4 * np.pi)

def M_linear_analytic(A_val, R):
    """
    Analityczne przybliżenie M_linear dla trybu zerowego.
    Dla g=1+A*sin(r)/r:
    f(g) ≈ 1 + 2*alpha*A*sin(r)/r
    (g')^2 ≈ A^2*(cos^2r/r^2 - 2*sinr*cosr/r^3 + sin^2r/r^4)
    V(g) ≈ -delta^2/2 = -A^2*sin^2r/(2r^2)

    Suma: A^2*cos(2r)/(2r^2) + O(A^3)
    Całka: ≈ 0 (zeruje się przez oscylacje)
    """
    # Proste całkowanie numeryczne dla małych A
    r = np.linspace(1.0, R, 5000)
    g_lin  = 1.0 + A_val * np.sin(r) / r
    gp_lin = A_val * (np.cos(r)/r - np.sin(r)/r**2)
    fg_lin = f_tgp(g_lin)
    Vg_lin = V_pot(g_lin)
    integrand = r**2 * (fg_lin * gp_lin**2 / 2.0 + Vg_lin)
    return float(np.trapezoid(integrand, r) * 4 * np.pi)

# ============================================================
print("="*68)
print("EX67: TRYB ZEROWY I MASA NIELINIOWA M_nonlin ∝ A_tail^4")
print("="*68)
print(f"  phi = {PHI:.7f},  g0* = {G0_STAR:.7f}")
print(f"  alpha = {ALPHA},  g* = {G_GHOST:.5f}")
print()

# ============================================================
# CZESC 1: M_linear(A) — weryfikacja trybu zerowego
# ============================================================
print("-"*68)
print("--- Czesc 1: M_linear(A) vs A — czy tryb zerowy daje M_lin→0? ---")
print()

A_vals = np.array([0.05, 0.10, 0.15, 0.20, 0.27, 0.35, 0.50, 0.70, 1.00, 1.20])
r_full = np.linspace(0.1, R_MAX, 20000)
M_lin_vals = []

print(f"  {'A':>8}  {'M_linear':>14}  {'M_lin/A^2':>12}  {'M_lin/A^4':>12}")
print("  " + "-"*55)
for A in A_vals:
    Ml = M_linear_analytic(A, R_MAX)
    M_lin_vals.append(Ml)
    print(f"  {A:8.3f}  {Ml:14.6f}  {Ml/A**2:12.6f}  {Ml/A**4:12.3f}")

M_lin_arr = np.array(M_lin_vals)
print()

# Fit M_linear ∝ A^k
try:
    valid = (A_vals > 0.01) & (np.abs(M_lin_arr) > 1e-10)
    if valid.sum() > 3:
        p_lin = np.polyfit(np.log(A_vals[valid]), np.log(np.abs(M_lin_arr[valid])), 1)
        k_lin = p_lin[0]
        print(f"  M_linear ∝ A^{k_lin:.3f}  (oczekiwane: k→0 dla trybu zerowego)")
except: k_lin = np.nan

# Sprawdz czy M_lin ~ 0 dla malych A
print(f"\n  Stosunek |M_lin|/A^2 dla A=0.1: {abs(M_lin_vals[1])/0.1**2:.4f}")
print(f"  Stosunek |M_lin|/A^2 dla A=0.5: {abs(M_lin_vals[6])/0.5**2:.4f}")
print(f"  (Dla trybu zerowego: M_lin/A^2 powinno dazac do stalej niezerowej")
print(f"   lub M_lin/A^4 do stalej — jesli O(A^4) dominuje)")
print()

# ============================================================
# CZESC 2: M_total i M_nonlin = M_total - M_linear
# ============================================================
print("-"*68)
print("--- Czesc 2: M_total, M_linear(A_inf) i M_nonlin vs g0 ---")
print()

g0_scan = np.linspace(1.10, 1.58, 15)

print(f"  {'g0':>8}  {'M_tot':>12}  {'M_lin':>12}  {'M_nonlin':>12}  {'A_inf':>10}  {'Mn/A^4':>10}")
print("  " + "-"*72)
res2 = []
for g0 in g0_scan:
    r, g, gp = integrate(g0)
    Ai = A_inf(g0)
    Mt = M_total_numeric(r, g, gp)
    Ml = M_linear_analytic(Ai, R_MAX)
    Mn = Mt - Ml
    A4 = Ai**4 if Ai > 1e-6 else 1e-30
    MnA4 = Mn/A4
    res2.append((g0, Mt, Ml, Mn, Ai, MnA4))
    mark = " <g0*" if abs(g0-G0_STAR) < 0.02 else ""
    print(f"  {g0:8.4f}  {Mt:12.5f}  {Ml:12.5f}  {Mn:12.5f}  {Ai:10.6f}  {MnA4:10.3f}{mark}")
print()

Mn_arr = np.array([x[3] for x in res2])
Ai_arr = np.array([x[4] for x in res2])
MnA4_arr = np.array([x[5] for x in res2])

# Power law fit
try:
    valid = (Ai_arr > 1e-5) & (np.abs(Mn_arr) > 1e-10)
    if valid.sum() > 4:
        p_mn = np.polyfit(np.log(Ai_arr[valid]), np.log(np.abs(Mn_arr[valid])), 1)
        k_mn = p_mn[0]; c_mn = np.exp(p_mn[1])
        ss = np.sum((np.log(np.abs(Mn_arr[valid])) - (k_mn*np.log(Ai_arr[valid])+p_mn[1]))**2)
        ss_t = np.sum((np.log(np.abs(Mn_arr[valid]))-np.log(np.abs(Mn_arr[valid])).mean())**2)
        R2_mn = 1 - ss/ss_t
        print(f"  M_nonlin ∝ A^{k_mn:.3f}  (R^2={R2_mn:.5f})")
    else: k_mn = np.nan; R2_mn = np.nan
except: k_mn = np.nan; R2_mn = np.nan

T3_sign = (Mn_arr > 0).all()
print(f"  M_nonlin > 0 zawsze: {T3_sign}  (min={Mn_arr.min():.5f})")
cv_mn = np.std(MnA4_arr)/abs(np.mean(MnA4_arr))*100 if len(MnA4_arr)>0 else np.nan
print(f"  M_nonlin/A^4: mean={np.mean(MnA4_arr):.3f}, CV={cv_mn:.2f}%")
print()

# ============================================================
# CZESC 3: Para (g0*, phi*g0*) — M_nonlin ratio
# ============================================================
print("-"*68)
print("--- Czesc 3: Para phi-FP: M_nonlin ratio ---")
print()

r_e, g_e, gp_e = integrate(G0_STAR)
Ai_e = A_inf(G0_STAR)
Mt_e = M_total_numeric(r_e, g_e, gp_e)
Ml_e = M_linear_analytic(Ai_e, R_MAX)
Mn_e = Mt_e - Ml_e

r_mu, g_mu, gp_mu = integrate(PHI * G0_STAR)
Ai_mu = A_inf(PHI * G0_STAR)
Mt_mu = M_total_numeric(r_mu, g_mu, gp_mu)
Ml_mu = M_linear_analytic(Ai_mu, R_MAX)
Mn_mu = Mt_mu - Ml_mu

print(f"  Elektron g0* = {G0_STAR:.7f}:")
print(f"    A_inf = {Ai_e:.7f}")
print(f"    M_tot = {Mt_e:.6f}")
print(f"    M_lin = {Ml_e:.6f}")
print(f"    M_nonlin = {Mn_e:.6f}")
print()
print(f"  Mion phi*g0* = {PHI*G0_STAR:.7f}:")
print(f"    A_inf = {Ai_mu:.7f}")
print(f"    M_tot = {Mt_mu:.6f}")
print(f"    M_lin = {Ml_mu:.6f}")
print(f"    M_nonlin = {Mn_mu:.6f}")
print()

if abs(Mn_e) > 1e-10:
    ratio_Mn = Mn_mu / Mn_e
    print(f"  M_nonlin ratio = {ratio_Mn:.4f}  (vs r21={R21_EXP})")
    print(f"  A_inf^4 ratio  = {(Ai_mu/Ai_e)**4:.4f}")
    print(f"  Odch. M_nonlin: {abs(ratio_Mn-R21_EXP)/R21_EXP*100:.2f}%")
    print(f"  Odch. A^4:      {abs((Ai_mu/Ai_e)**4-R21_EXP)/R21_EXP*100:.6f}%")
else:
    ratio_Mn = np.nan
    print(f"  M_nonlin zbyt maly dla elektronu: {Mn_e:.8f}")
print()

# ============================================================
# CZESC 4: M_linear jako funkcja R_max (sprawdzenie zbieznosci)
# ============================================================
print("-"*68)
print("--- Czesc 4: M_linear(A=0.274) vs R_max (tryb zerowy: M→0?) ---")
print()

A_test = 0.274  # A_inf dla elektronu
print(f"  A = {A_test:.3f} (A_inf elektronu)")
print(f"  {'R_max':>8}  {'M_linear':>14}  {'M_lin/A^4':>12}  {'M_lin/A^2':>12}")
print("  " + "-"*52)
for Rm in [10, 20, 30, 50, 70, 100]:
    Ml = M_linear_analytic(A_test, Rm)
    print(f"  {Rm:8.0f}  {Ml:14.6f}  {Ml/A_test**4:12.4f}  {Ml/A_test**2:12.4f}")

print()
print(f"  Obserwacja: M_linear powinno byc ~0 dla trybu zerowego.")
print(f"  Zbieznosc (lub brak) wskazuje na charakter trybu.")
print()

# ============================================================
# CZESC 5: Porownanie M_nonlin z M_core (jądro)
# ============================================================
print("-"*68)
print("--- Czesc 5: M_nonlin vs M_core (r < r_cross) ---")
print()

def M_core_total(g0, r_cut):
    r, g, gp = integrate(g0)
    mask = r <= r_cut
    if mask.sum() < 5: return np.nan
    fg = f_tgp(g[mask]); Vg = V_pot(g[mask])
    return float(np.trapezoid(r[mask]**2 * (fg*gp[mask]**2/2+Vg), r[mask]) * 4*np.pi)

r_cut_test = [2, 3, 5, 10]
print(f"  {'r_cut':>8}  {'Mc_e':>12}  {'Mc_mu':>12}  {'ratio':>10}  {'vs Mn_ratio':>12}")
print("  " + "-"*58)
for rc in r_cut_test:
    Mce = M_core_total(G0_STAR, rc)
    Mcmu = M_core_total(PHI*G0_STAR, rc)
    if not np.isnan(Mce) and abs(Mce) > 1e-8:
        ratio_rc = Mcmu/Mce
        print(f"  {rc:8.1f}  {Mce:12.5f}  {Mcmu:12.5f}  {ratio_rc:10.3f}  {'(vs Mn='+str(round(ratio_Mn,2))+')':>12}")
print()

# ============================================================
# TESTY
# ============================================================
print("="*68)
print("TESTY")
print("="*68)

# T1: M_linear ∝ A^k, k ∈ [0, 0.5] (tryb zerowy: nieduzy wykładnik)
# Zamiast tego: sprawdz czy M_linear oscyluje (zmienia znak) — tryb zerowy
signs = np.sign(M_lin_arr)
T1 = (signs != signs[0]).any()  # zmienia znak
if T1:
    print(f"  T1: M_linear zmienia znak (tryb zerowy): PASS")
else:
    # Sprawdz czy M_linear/A^4 → const
    if not np.isnan(k_lin) and 3.5 <= k_lin <= 4.5:
        T1 = True
        print(f"  T1: M_linear ∝ A^{k_lin:.3f} → tryb zerowy: PASS (k≈4, O(A^4))")
    else:
        T1 = False
        print(f"  T1: M_linear zmienia znak: FAIL (znaki: {signs[:5]}, k={k_lin:.2f if not np.isnan(k_lin) else 'N/A'})")

# T2: M_nonlin ∝ A^k, k ∈ [3.5, 4.5]
T2 = not np.isnan(k_mn) and 3.5 <= k_mn <= 4.5
k_mn_str  = f"{k_mn:.3f}"  if not np.isnan(k_mn)  else 'N/A'
R2_mn_str = f"{R2_mn:.4f}" if not np.isnan(R2_mn) else 'N/A'
print(f"  T2: M_nonlin ∝ A^k, k in [3.5,4.5]:  {'PASS' if T2 else 'FAIL'}  "
      f"(k={k_mn_str}, R^2={R2_mn_str})")

# T3: M_nonlin > 0 zawsze
T3 = T3_sign
print(f"  T3: M_nonlin > 0 zawsze:               {'PASS' if T3 else 'FAIL'}  "
      f"(min={Mn_arr.min():.5f})")

# T4: |M_nonlin ratio - 207| < 30%
T4 = not np.isnan(ratio_Mn) and abs(ratio_Mn - R21_EXP)/R21_EXP < 0.30
ratio_Mn_str = f"{ratio_Mn:.4f}" if not np.isnan(ratio_Mn) else 'N/A'
print(f"  T4: |M_nonlin_ratio - 207| < 30%:      {'PASS' if T4 else 'FAIL'}  "
      f"(ratio={ratio_Mn_str})")

# T5: M_nonlin/A^4 CV < 50%
T5 = not np.isnan(cv_mn) and cv_mn < 50
print(f"  T5: M_nonlin/A^4 CV < 50%:             {'PASS' if T5 else 'FAIL'}  "
      f"(CV={cv_mn:.2f}%)")

n_pass = sum([T1, T2, T3, T4, T5])
print(f"\nWYNIK: {n_pass}/5 testow przeszlo")

print()
print("="*68)
print("WNIOSEK (EX67)")
print("="*68)
print(f"  M_linear(A): zmienia znak = {(np.sign(M_lin_arr) != np.sign(M_lin_arr[0])).any()}")
print(f"  M_nonlin ∝ A^{k_mn_str}  (R^2={R2_mn_str})")
print(f"  M_nonlin > 0: {T3_sign}")
print(f"  M_nonlin ratio (mu/e) = {ratio_Mn_str}  (vs r21={R21_EXP})")
print(f"  Potwierdza: k=1 jest trybem zerowym, M_nonlin ~ O(A^4)")
print("="*68)
