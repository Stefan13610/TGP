"""
EX66: M_CORE — FIZYCZNA MASA TGP JAKO CAŁKA JĄDRA
==================================================
STATUS: LEGACY-TRANSLATIONAL

This script belongs to the older `A_tail^4` / `g0* ~ 1.24` mass-selection
program. Keep it as historical translation and diagnostic context, not as a
canonical synchronized entry point for current `nbody`.

Hipoteza (ex65): Fizyczna masa solitonu TGP to ENERGIA JĄDRA,
nie ogona. Ogon ma zerową energię netto (kasowanie kin+pot).

Definicja: M_core(g0) = 4pi * int_0^{r_cross} r^2 f(g) (g')^2/2 dr
gdzie r_cross = pierwsza chwila gdy g(r_cross) = 1.0 po r=0.

(Alternatywnie: int energii calkowitej f(g)(g')^2/2 + V(g) w rdzeniu.)

Pytania:
  1. Czy M_core(g0) ∝ A_tail^4? (test fundamentalny)
  2. M_core(phi*g0*) / M_core(g0*) =? 207
  3. Jak M_core zalezy od wyboru r_cross?
     (czy r_cross = naturalna granica jądra?)
  4. Porownanie z ex57 (M_core/M_e = 210 przy g0=3.678)
  5. Czy M_core(g0_e, g0_mu) daje r21=207 dla naturalnych g0?

Testy:
  T1: M_core ∝ A_tail^k, k ∈ [3.0, 5.0]  (≈ A^4)
  T2: |M_core(phi*g0*)/M_core(g0*) - 207| / 207 < 0.20  (20%)
  T3: r_cross(g0*) ∈ [1.0, 5.0]  (sensowna granica)
  T4: M_core/A^4 = const ± 30%  (staly wspolczynnik)
  T5: M_core(g0=3.678, r<3) / M_core(g0=1.24, r<3) ∈ [190, 230]  (ex57 test)
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
R_MAX    = 60.0  # wystarczy dla M_core
G0_STAR  = 1.24771102
WIN_LIST  = [16, 22, 28, 36, 46]
WIN_WIDTH = 12.0

ALPHA    = 2.0
G_GHOST  = np.exp(-1.0 / (2.0 * ALPHA))
G_BOUNCE = G_GHOST + 0.005

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

def ev_g1_down(r, y): return y[0] - 1.0   # g crosses 1 going down
ev_g1_down.terminal  = False
ev_g1_down.direction = -1

def integrate(g0, maxb=8):
    y0 = [g0, 0.0]; r0 = 0.0; ra, ga, gpa = [], [], []
    n_bounces = 0
    for _ in range(maxb + 1):
        sol = solve_ivp(rhs, (r0, R_MAX), y0, events=[ev_ghost, ev_g1_down],
                        dense_output=True, rtol=1e-9, atol=1e-11, max_step=0.05)
        ra.append(sol.t); ga.append(sol.y[0]); gpa.append(sol.y[1])
        if sol.status == 1 and len(sol.t_events[0]) > 0:
            rh = sol.t_events[0][-1]; st = sol.sol(rh)
            y0 = [st[0], -st[1]]; r0 = rh
            n_bounces += 1
        else: break
    r = np.concatenate(ra); g = np.concatenate(ga); gp = np.concatenate(gpa)
    idx = np.argsort(r)
    return r[idx], g[idx], gp[idx], n_bounces

def V_pot(g):
    return g**3/3 - g**4/4 - (1./3 - 1./4)

def f_tgp(g):
    return 1.0 + 2.0 * ALPHA * np.log(np.maximum(g, 1e-10))

def find_r_cross(r_arr, g_arr, after_bounce=True):
    """Znajdz pierwsze przejscie g=1 (po odbiciu lub bez)."""
    # Jesli g0 > 1 i nie ma odbicia: szukaj pierwszego g < 1
    # Jesli jest odbicie: szukaj pierwszego g=1 PO odbiciu
    for i in range(1, len(g_arr)):
        if g_arr[i-1] > 1.0 and g_arr[i] <= 1.0:
            # Interpoluj
            f = (1.0 - g_arr[i-1]) / (g_arr[i] - g_arr[i-1])
            return r_arr[i-1] + f * (r_arr[i] - r_arr[i-1])
    return r_arr[-1]  # nie znaleziono

def compute_Mcore(g0, r_cut=None, use_r_cross=True):
    """
    Oblicz M_core = 4pi int_0^{r_cut} r^2 [f(g)(g')^2/2 + V(g)] dr
    Jesli use_r_cross=True: r_cut = r_cross (pierwsza chwila g=1)
    """
    r, g, gp, nb = integrate(g0)

    if use_r_cross:
        r_cut = find_r_cross(r, g)

    mask = r <= r_cut
    if mask.sum() < 5: return np.nan, r_cut, nb

    r_m  = r[mask]; g_m  = g[mask]; gp_m = gp[mask]
    fg_m = f_tgp(g_m); Vg_m = V_pot(g_m)

    M_kin_core = float(np.trapezoid(r_m**2 * fg_m * gp_m**2 / 2, r_m) * 4*np.pi)
    M_pot_core = float(np.trapezoid(r_m**2 * Vg_m, r_m) * 4*np.pi)
    M_tot_core = M_kin_core + M_pot_core

    return M_tot_core, r_cut, nb

def compute_Mkin_core(g0, r_cut=None, use_r_cross=True):
    """Tylko kinetyczna czesc jądra (bez potencjalu)."""
    r, g, gp, nb = integrate(g0)
    if use_r_cross:
        r_cut = find_r_cross(r, g)
    mask = r <= r_cut
    if mask.sum() < 5: return np.nan, r_cut, nb
    r_m = r[mask]; g_m = g[mask]; gp_m = gp[mask]
    fg_m = f_tgp(g_m)
    M_kin = float(np.trapezoid(r_m**2 * fg_m * gp_m**2 / 2, r_m) * 4*np.pi)
    return M_kin, r_cut, nb

def A_inf(g0):
    r, g, gp, _ = integrate(g0)
    def fit_win(r, g, rL, rR):
        mask = (r >= rL) & (r <= rR)
        if np.sum(mask) < 20: return 0., 0., 0.
        rf = r[mask]; df = (g[mask]-1)*rf
        M = np.column_stack([np.cos(rf), np.sin(rf)])
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

# ============================================================
print("="*68)
print("EX66: M_CORE — FIZYCZNA MASA TGP")
print("="*68)
print(f"  phi = {PHI:.7f},  g0* = {G0_STAR:.7f}")
print(f"  alpha = {ALPHA},  g* = {G_GHOST:.5f}")
print()

# ============================================================
# CZESC 1: M_core(g0) vs A_inf(g0) — sektor n=0
# ============================================================
print("-"*68)
print("--- Czesc 1: M_core(r_cross) i M_kin_core vs A_inf (n=0) ---")
print()

g0_scan1 = np.linspace(1.10, 1.58, 18)
res1 = []
print(f"  {'g0':>8}  {'r_cross':>9}  {'M_tot_core':>12}  {'M_kin_core':>12}  {'A_inf':>10}  {'Mc/A^4':>10}")
print("  " + "-"*68)
for g0 in g0_scan1:
    Mc_tot, rc, nb  = compute_Mcore(g0)
    Mc_kin, rc2, _  = compute_Mkin_core(g0)
    Ai = A_inf(g0)
    A4 = Ai**4
    ratio_cA4 = Mc_kin/A4 if A4 > 1e-12 else np.nan
    res1.append((g0, rc, Mc_tot, Mc_kin, Ai, ratio_cA4))
    print(f"  {g0:8.4f}  {rc:9.4f}  {Mc_tot:12.5f}  {Mc_kin:12.5f}  {Ai:10.6f}  {ratio_cA4:10.3f}")
print()

g0_arr1 = np.array([x[0] for x in res1])
Mc_arr1 = np.array([x[2] for x in res1])  # M_tot_core
Mk_arr1 = np.array([x[3] for x in res1])  # M_kin_core
Ai_arr1 = np.array([x[4] for x in res1])

# Power law fits
def power_fit(A_vals, Y_vals):
    valid = (A_vals > 1e-5) & (Y_vals > 1e-10)
    if valid.sum() < 4: return np.nan, np.nan, np.nan
    logA = np.log(A_vals[valid]); logY = np.log(np.abs(Y_vals[valid]))
    p = np.polyfit(logA, logY, 1)
    k = p[0]; c = np.exp(p[1])
    ss_res = np.sum((logY-(k*logA+p[1]))**2)
    ss_tot = np.sum((logY-logY.mean())**2)
    R2 = 1 - ss_res/ss_tot
    return k, c, R2

k_mc, c_mc, R2_mc = power_fit(Ai_arr1, Mc_arr1)
k_mk, c_mk, R2_mk = power_fit(Ai_arr1, Mk_arr1)

print(f"  Power law fits:")
print(f"    M_tot_core ∝ A^{k_mc:.3f}  (R^2={R2_mc:.5f})")
print(f"    M_kin_core ∝ A^{k_mk:.3f}  (R^2={R2_mk:.5f})")
print()

# ============================================================
# CZESC 2: Kluczowa para (g0*, phi*g0*) — M_core ratio
# ============================================================
print("-"*68)
print("--- Czesc 2: Para phi-FP: M_core ratio ---")
print()

Mc_e, rc_e, nb_e = compute_Mcore(G0_STAR)
Mk_e, _, _        = compute_Mkin_core(G0_STAR)
Ai_e              = A_inf(G0_STAR)

Mc_mu, rc_mu, nb_mu = compute_Mcore(PHI * G0_STAR)
Mk_mu, _, _          = compute_Mkin_core(PHI * G0_STAR)
Ai_mu                = A_inf(PHI * G0_STAR)

print(f"  Elektron g0* = {G0_STAR:.7f}:")
print(f"    r_cross = {rc_e:.5f},  n_bounce = {nb_e}")
print(f"    M_tot_core = {Mc_e:.6f}")
print(f"    M_kin_core = {Mk_e:.6f}")
print(f"    A_inf      = {Ai_e:.7f}")
print()
print(f"  Mion phi*g0* = {PHI*G0_STAR:.7f}:")
print(f"    r_cross = {rc_mu:.5f},  n_bounce = {nb_mu}")
print(f"    M_tot_core = {Mc_mu:.6f}")
print(f"    M_kin_core = {Mk_mu:.6f}")
print(f"    A_inf      = {Ai_mu:.7f}")
print()

if abs(Mc_e) > 1e-8:
    ratio_Mc = Mc_mu / Mc_e
    print(f"  M_tot_core ratio = {ratio_Mc:.4f}  (vs r21={R21_EXP})")
    print(f"  M_kin_core ratio = {Mk_mu/Mk_e:.4f}  (vs r21={R21_EXP})")
    print(f"  A_inf^4 ratio    = {(Ai_mu/Ai_e)**4:.4f}  (vs r21={R21_EXP})")
    print(f"  Odch. M_tot_core: {abs(ratio_Mc-R21_EXP)/R21_EXP*100:.2f}%")
    print(f"  Odch. M_kin_core: {abs(Mk_mu/Mk_e-R21_EXP)/R21_EXP*100:.2f}%")
print()

# ============================================================
# CZESC 3: M_core jako funkcja r_cut (ile wynosi granica?)
# ============================================================
print("-"*68)
print("--- Czesc 3: M_core vs r_cut (elektron i mion) ---")
print()
print(f"  {'r_cut':>8}  {'Mc_e':>10}  {'Mc_mu':>10}  {'ratio':>10}  {'odch%':>8}")
print("  " + "-"*50)
for r_cut in [2.0, 3.0, 4.0, 5.0, 7.0, 10.0, 15.0, 20.0]:
    Mc_e2, _, _ = compute_Mcore(G0_STAR, r_cut=r_cut, use_r_cross=False)
    Mc_mu2, _, _ = compute_Mcore(PHI*G0_STAR, r_cut=r_cut, use_r_cross=False)
    if abs(Mc_e2) > 1e-8 and not np.isnan(Mc_e2):
        ratio_rcut = Mc_mu2/Mc_e2
        odch = (ratio_rcut-R21_EXP)/R21_EXP*100
        print(f"  {r_cut:8.1f}  {Mc_e2:10.4f}  {Mc_mu2:10.4f}  {ratio_rcut:10.3f}  {odch:8.2f}%")
    else:
        print(f"  {r_cut:8.1f}  {Mc_e2:10.4f}  {'N/A':>10}")
print()

# ============================================================
# CZESC 4: M_kin_core (tylko kinetyczna) jako funkcja r_cut
# ============================================================
print("-"*68)
print("--- Czesc 4: M_kin_core (czysta kinetyczna) vs r_cut ---")
print()
print(f"  {'r_cut':>8}  {'Mk_e':>10}  {'Mk_mu':>10}  {'ratio':>10}  {'odch%':>8}")
print("  " + "-"*50)
for r_cut in [2.0, 3.0, 4.0, 5.0, 7.0, 10.0, 15.0]:
    Mk_e2, _, _ = compute_Mkin_core(G0_STAR, r_cut=r_cut, use_r_cross=False)
    Mk_mu2, _, _ = compute_Mkin_core(PHI*G0_STAR, r_cut=r_cut, use_r_cross=False)
    if Mk_e2 > 1e-8:
        ratio_rcut = Mk_mu2/Mk_e2
        odch = (ratio_rcut-R21_EXP)/R21_EXP*100
        print(f"  {r_cut:8.1f}  {Mk_e2:10.4f}  {Mk_mu2:10.4f}  {ratio_rcut:10.3f}  {odch:8.2f}%")
print()

# ============================================================
# CZESC 5: Test ex57 (g0=3.678, g0=1.24, r<3)
# ============================================================
print("-"*68)
print("--- Czesc 5: Reprodukcja ex57 (M_core r<3, g0=3.678 vs 1.24) ---")
print()

Mk_e57, _, _ = compute_Mkin_core(1.24, r_cut=3.0, use_r_cross=False)
Mk_mu57, _, _ = compute_Mkin_core(3.678, r_cut=3.0, use_r_cross=False)
Mk_mu57b, _, _ = compute_Mkin_core(4.259, r_cut=3.0, use_r_cross=False)

print(f"  M_kin_core(r<3, g0=1.24)  = {Mk_e57:.4f}")
print(f"  M_kin_core(r<3, g0=3.678) = {Mk_mu57:.4f}  ratio = {Mk_mu57/Mk_e57:.2f}")
print(f"  M_kin_core(r<3, g0=4.259) = {Mk_mu57b:.4f}  ratio = {Mk_mu57b/Mk_e57:.2f}")
print(f"  (ex57: 210.9 przy g0=3.678, 202.6 przy g0=4.259)")
print()

# A_inf dla tych g0
Ai_57e = A_inf(1.24); Ai_57mu = A_inf(3.678)
print(f"  A_inf(1.24) = {Ai_57e:.6f}")
print(f"  A_inf(3.678) = {Ai_57mu:.6f}")
print(f"  (A_inf(3.678)/A_inf(1.24))^4 = {(Ai_57mu/Ai_57e)**4:.2f}")
print()

# ============================================================
# CZESC 6: M_kin_core vs A^4 — pelny skan
# ============================================================
print("-"*68)
print("--- Czesc 6: M_kin_core(r_cross) vs A^4 — pelny skan ---")
print()

g0_scan6 = np.linspace(1.10, 2.25, 25)
res6 = []
for g0 in g0_scan6:
    Mk, rc, nb = compute_Mkin_core(g0)
    Ai = A_inf(g0)
    A4 = Ai**4 if Ai > 1e-6 else 1e-30
    MA4 = Mk/A4 if A4 > 1e-12 else np.nan
    res6.append((g0, rc, nb, Mk, Ai, MA4))

Mk6 = np.array([x[3] for x in res6])
Ai6 = np.array([x[4] for x in res6])

k6, c6, R2_6 = power_fit(Ai6, Mk6)
MA4_arr6 = np.array([x[5] for x in res6 if not np.isnan(x[5])])

print(f"  Power law: M_kin_core(r_cross) = {c6:.4f} * A^{k6:.3f}  R^2={R2_6:.5f}")
cv6 = np.std(MA4_arr6)/np.mean(MA4_arr6)*100 if len(MA4_arr6) > 0 else np.nan
print(f"  M_kin_core/A^4: mean={np.mean(MA4_arr6):.3f}, CV={cv6:.2f}%")
print()
print(f"  {'g0':>7}  {'r_cross':>8}  {'nb':>3}  {'M_kin_core':>12}  {'A_inf':>10}  {'Mc/A^4':>10}")
print("  " + "-"*60)
for g0, rc, nb, Mk, Ai, MA4 in res6:
    mark = " <g0*" if abs(g0-G0_STAR) < 0.03 else (" <phi*g0*" if abs(g0-PHI*G0_STAR) < 0.04 else "")
    print(f"  {g0:7.4f}  {rc:8.4f}  {nb:3d}  {Mk:12.5f}  {Ai:10.6f}  {MA4:10.3f}{mark}")
print()

# ============================================================
# TESTY
# ============================================================
print("="*68)
print("TESTY")
print("="*68)

# T1: M_kin_core ∝ A^k, k ∈ [3.0, 5.0]
T1 = not np.isnan(k6) and 3.0 <= k6 <= 5.0
print(f"  T1: M_kin_core ∝ A^k, k in [3,5]:   {'PASS' if T1 else 'FAIL'}  (k={k6:.3f}, R^2={R2_6:.4f})")

# T2: |M_kin_core ratio - 207| / 207 < 0.20
T2 = abs(Mk_mu/Mk_e - R21_EXP)/R21_EXP < 0.20
print(f"  T2: |Mc_ratio - 207| < 20%:          {'PASS' if T2 else 'FAIL'}  (ratio={Mk_mu/Mk_e:.3f})")

# T3: r_cross(g0*) ∈ [1.0, 5.0]
T3 = 1.0 <= rc_e <= 5.0
print(f"  T3: r_cross(g0*) in [1,5]:           {'PASS' if T3 else 'FAIL'}  (r_cross={rc_e:.4f})")

# T4: M_kin_core/A^4 CV < 30%
T4 = not np.isnan(cv6) and cv6 < 30
print(f"  T4: M_kin_core/A^4 CV < 30%:         {'PASS' if T4 else 'FAIL'}  (CV={cv6:.2f}%)")

# T5: M_kin_core(r<3, g0=3.678) / M_kin_core(r<3, g0=1.24) ∈ [190, 230]
T5 = 190 <= Mk_mu57/Mk_e57 <= 230
print(f"  T5: ex57 test [190,230]:              {'PASS' if T5 else 'FAIL'}  (ratio={Mk_mu57/Mk_e57:.2f})")

n_pass = sum([T1, T2, T3, T4, T5])
print(f"\nWYNIK: {n_pass}/5 testow przeszlo")

print()
print("="*68)
print("WNIOSEK (EX66)")
print("="*68)
print(f"  M_kin_core(r_cross) ∝ A^{k6:.3f}  (R^2={R2_6:.4f})")
print(f"  M_kin_core_mu / M_kin_core_e = {Mk_mu/Mk_e:.4f}  (vs r21={R21_EXP})")
print(f"  r_cross(g0*)  = {rc_e:.4f}")
print(f"  r_cross(phi*g0*) = {rc_mu:.4f}")
print(f"  M_kin_core/A^4: CV={cv6:.2f}%")
print("="*68)
