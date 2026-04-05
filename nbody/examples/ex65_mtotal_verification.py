"""
EX65: WERYFIKACJA M_total = M_kin + M_pot ∝ A_tail^4
======================================================
Hipoteza analityczna (ex64):
  - M_raw (kinetyczna) ∝ A_tail^2 (ogon dominuje, rozbiezny)
  - M_pot (potencjalna) ∝ A_tail^2 z opozytnym znakiem w ogonie
  - Suma M_total = M_kin + M_pot: ogon kasuje sie (∝ cos(2r)/r^2 → 0)
  - Pozostaje rdzen: M_total ∝ A_tail^4

Definicje:
  M_kin = 4pi ∫ r^2 f(g) (g')^2/2 dr
  M_pot = 4pi ∫ r^2 V(g) dr   gdzie V(g) = g^3/3 - g^4/4 - 1/12
  M_total = M_kin + M_pot

Testujemy:
  1. M_kin ∝ A_tail^2  (potwierdzenie ex64)
  2. M_pot ∝ A_tail^2  z UJEMNA wartosia (ogon: V(1+d) ≈ -d^2/2)
  3. M_total ∝ A_tail^4  (po supresji ogona)
  4. M_total(phi*g0*) / M_total(g0*) ≈ 207
  5. M_total / A_tail^4 = const (staly wspolczynnik proporcjonalnosci)

Testy:
  T1: M_kin ∝ A^k1 gdzie k1 ∈ [1.8, 2.2]  (≈ A^2)
  T2: M_pot ∝ A^k2 gdzie k2 ∈ [1.8, 2.2] i M_pot < 0  (ujemna!)
  T3: M_total ∝ A^k3 gdzie k3 ∈ [3.5, 4.5]  (≈ A^4)
  T4: |M_total(phi*g0*)/M_total(g0*) - 207| / 207 < 0.10  (10%)
  T5: M_total/A^4 ma stabilny wspolczynnik: zmiennosc < 30%
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
G0_STAR  = 1.24771102
WIN_LIST  = [16, 22, 28, 36, 46, 58, 72]
WIN_WIDTH = 14.0

ALPHA    = 2.0
G_GHOST  = np.exp(-1.0 / (2.0 * ALPHA))   # ~0.7788
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
    """V(g) = int_1^g s^2(1-s)ds = g^3/3 - g^4/4 - (1/3 - 1/4)"""
    return g**3/3 - g**4/4 - (1.0/3 - 1.0/4)

def f_tgp(g):
    """f(g) = 1 + 2*alpha*ln(g)"""
    return 1.0 + 2.0 * ALPHA * np.log(np.maximum(g, 1e-10))

def compute_energies(g0):
    """Zwraca (M_kin, M_pot, M_total, A_inf)."""
    r, g, gp = integrate(g0)

    # M_kin = 4pi * int r^2 * f(g) * gp^2 / 2 dr
    fg     = f_tgp(g)
    M_kin  = float(np.trapezoid(r**2 * fg * gp**2 / 2.0, r) * 4 * np.pi)

    # M_pot = 4pi * int r^2 * V(g) dr
    Vg     = V_pot(g)
    M_pot  = float(np.trapezoid(r**2 * Vg, r) * 4 * np.pi)

    M_tot  = M_kin + M_pot

    # A_inf
    def fit_win(r, g, rL, rR):
        mask = (r >= rL) & (r <= rR)
        if np.sum(mask) < 20: return 0., 0., 0.
        rf = r[mask]; df = (g[mask] - 1.0) * rf
        M  = np.column_stack([np.cos(rf), np.sin(rf)])
        bc = np.linalg.lstsq(M, df, rcond=None)[0]
        B, C = bc
        return float(np.sqrt(B**2+C**2)), float(B), float(C)

    Av, rv = [], []
    for rL in WIN_LIST:
        if rL + WIN_WIDTH > r[-1]: break
        A, _, _ = fit_win(r, g, rL, rL+WIN_WIDTH)
        if A > 0.002: Av.append(A); rv.append(rL)
    if len(Av) < 3:
        Ainf = Av[-1] if Av else 0.
    else:
        try:
            p, _ = curve_fit(lambda x,ai,a: ai*(1+a/x), rv, Av,
                             p0=[Av[-1],0.], maxfev=2000)
            Ainf = float(p[0])
        except: Ainf = float(Av[-1])

    return M_kin, M_pot, M_tot, Ainf

# ============================================================
print("="*68)
print("EX65: WERYFIKACJA M_total = M_kin + M_pot ∝ A_tail^4")
print("="*68)
print(f"  phi = {PHI:.7f},  g0* = {G0_STAR:.7f}")
print(f"  alpha = {ALPHA},  g* = {G_GHOST:.5f},  R_MAX = {R_MAX}")
print()

# ============================================================
# CZESC 1: Skan M_kin, M_pot, M_total, A_inf w sektorze n=0
# ============================================================
print("-"*68)
print("--- Czesc 1: Energetyka vs A_inf (sektor n=0, g0 ∈ [1.10, 1.60]) ---")
print()

g0_scan1 = np.linspace(1.10, 1.60, 20)
res1 = []
print(f"  {'g0':>8}  {'M_kin':>12}  {'M_pot':>12}  {'M_tot':>12}  {'A_inf':>10}  {'M_tot/A^4':>12}")
print("  " + "-"*75)
for g0 in g0_scan1:
    Mk, Mp, Mt, Ai = compute_energies(g0)
    A4 = Ai**4
    ratio_MA4 = Mt/A4 if A4 > 1e-12 else np.nan
    res1.append((g0, Mk, Mp, Mt, Ai, ratio_MA4))
    print(f"  {g0:8.4f}  {Mk:12.4f}  {Mp:12.4f}  {Mt:12.4f}  {Ai:10.6f}  {ratio_MA4:12.4f}")
print()

# Wyodrebnij tablice
g0_arr  = np.array([x[0] for x in res1])
Mk_arr  = np.array([x[1] for x in res1])
Mp_arr  = np.array([x[2] for x in res1])
Mt_arr  = np.array([x[3] for x in res1])
Ai_arr  = np.array([x[4] for x in res1])
MA4_arr = np.array([x[5] for x in res1])

# Power law fits
def power_fit(A_vals, Y_vals):
    """Fit Y = c * A^k."""
    valid = (A_vals > 1e-5) & (Y_vals > 1e-10)
    if valid.sum() < 4: return np.nan, np.nan, np.nan
    logA = np.log(A_vals[valid]); logY = np.log(np.abs(Y_vals[valid]))
    p = np.polyfit(logA, logY, 1)
    k = p[0]; c = np.exp(p[1])
    pred = logY - (k*logA + p[1])
    ss_res = np.sum(pred**2); ss_tot = np.sum((logY-logY.mean())**2)
    R2 = 1 - ss_res/ss_tot
    return k, c, R2

k_kin, c_kin, R2_kin = power_fit(Ai_arr, Mk_arr)
k_pot, c_pot, R2_pot = power_fit(Ai_arr, np.abs(Mp_arr))
k_tot, c_tot, R2_tot = power_fit(Ai_arr, Mt_arr)

print(f"  Power law fits (Y = c * A_inf^k):")
print(f"    M_kin:  k={k_kin:.4f}  c={c_kin:.4f}  R^2={R2_kin:.6f}")
print(f"    |M_pot|: k={k_pot:.4f}  c={c_pot:.4f}  R^2={R2_pot:.6f}")
print(f"    M_tot:  k={k_tot:.4f}  c={c_tot:.4f}  R^2={R2_tot:.6f}")
print()

print(f"  M_tot/A^4 range: [{np.nanmin(MA4_arr):.3f}, {np.nanmax(MA4_arr):.3f}]")
if not np.isnan(MA4_arr).all():
    valid_MA4 = MA4_arr[~np.isnan(MA4_arr)]
    cv = np.std(valid_MA4)/np.mean(valid_MA4) * 100
    print(f"  Wspolczynnik zmiennosci M_tot/A^4: {cv:.2f}%")
print()

# ============================================================
# CZESC 2: Analiza ogona — supresja kinetyczna + potencjalna
# ============================================================
print("-"*68)
print("--- Czesc 2: Analiza energetyki ogona (r > r_cross) ---")
print()

def energy_tail_contribution(g0, r_cut=10.0):
    """Wklady M_kin, M_pot od ogona (r > r_cut)."""
    r, g, gp = integrate(g0)
    mask_full = r > 0.5
    mask_tail = r > r_cut

    fg = f_tgp(g)
    Vg = V_pot(g)

    Mk_full = float(np.trapezoid(r[mask_full]**2 * fg[mask_full] * gp[mask_full]**2 / 2, r[mask_full]) * 4*np.pi)
    Mp_full = float(np.trapezoid(r[mask_full]**2 * Vg[mask_full], r[mask_full]) * 4*np.pi)

    Mk_tail = float(np.trapezoid(r[mask_tail]**2 * fg[mask_tail] * gp[mask_tail]**2 / 2, r[mask_tail]) * 4*np.pi)
    Mp_tail = float(np.trapezoid(r[mask_tail]**2 * Vg[mask_tail], r[mask_tail]) * 4*np.pi)

    Mk_core = Mk_full - Mk_tail
    Mp_core = Mp_full - Mp_tail

    return (Mk_full, Mp_full, Mk_tail, Mp_tail, Mk_core, Mp_core)

print(f"  Analiza g0* = {G0_STAR:.5f}:")
Mk_f, Mp_f, Mk_t, Mp_t, Mk_c, Mp_c = energy_tail_contribution(G0_STAR, r_cut=10.0)
print(f"  M_kin(full)  = {Mk_f:.6f}")
print(f"  M_pot(full)  = {Mp_f:.6f}")
print(f"  M_tot(full)  = {Mk_f+Mp_f:.6f}")
print(f"  M_kin(tail r>10) = {Mk_t:.6f}  ({Mk_t/Mk_f*100:.1f}% M_kin)")
print(f"  M_pot(tail r>10) = {Mp_t:.6f}  ({Mp_t/Mp_f*100:.1f}% M_pot)")
print(f"  Mtail(kin+pot)   = {Mk_t+Mp_t:.6f}  (supresja: {(Mk_t+Mp_t)/(Mk_f+Mp_f)*100:.1f}% M_tot)")
print(f"  M_kin(core r<10) = {Mk_c:.6f}")
print(f"  M_pot(core r<10) = {Mp_c:.6f}")
print(f"  M_tot(core r<10) = {Mk_c+Mp_c:.6f}")
print()

# Test supresji jako funkcji r_cut
print("  Supresja ogona jako funkcja r_cut:")
print(f"  {'r_cut':>8}  {'Mk_tail':>12}  {'Mp_tail':>12}  {'Mtail':>12}  {'frac%':>8}")
for rc in [5.0, 10.0, 20.0, 40.0, 60.0]:
    _, _, Mkt, Mpt, _, _ = energy_tail_contribution(G0_STAR, r_cut=rc)
    Mtail = Mkt + Mpt
    _, _, _, _, Mk_c2, Mp_c2 = energy_tail_contribution(G0_STAR, r_cut=rc)
    Mtot2 = Mk_c2 + Mp_c2 + Mtail
    frac  = abs(Mtail)/(abs(Mk_f+Mp_f)) * 100 if abs(Mk_f+Mp_f) > 1e-10 else np.nan
    print(f"  {rc:8.1f}  {Mkt:12.6f}  {Mpt:12.6f}  {Mtail:12.6f}  {frac:8.2f}%")
print()

# ============================================================
# CZESC 3: Kluczowa para (g0*, phi*g0*) — stosunek energii
# ============================================================
print("-"*68)
print("--- Czesc 3: Para phi-FP: M_total(phi*g0*) / M_total(g0*) ---")
print()

Mk_e, Mp_e, Mt_e, Ai_e = compute_energies(G0_STAR)
Mk_mu, Mp_mu, Mt_mu, Ai_mu = compute_energies(PHI * G0_STAR)

print(f"  Elektron g0* = {G0_STAR:.7f}:")
print(f"    M_kin = {Mk_e:.6f},  M_pot = {Mp_e:.6f},  M_tot = {Mt_e:.6f}")
print(f"    A_inf = {Ai_e:.7f},  A_inf^4 = {Ai_e**4:.6f}")
print()
print(f"  Mion phi*g0* = {PHI*G0_STAR:.7f}:")
print(f"    M_kin = {Mk_mu:.6f},  M_pot = {Mp_mu:.6f},  M_tot = {Mt_mu:.6f}")
print(f"    A_inf = {Ai_mu:.7f},  A_inf^4 = {Ai_mu**4:.6f}")
print()

ratio_M_tot = Mt_mu / Mt_e if abs(Mt_e) > 1e-10 else np.nan
ratio_A4    = (Ai_mu / Ai_e)**4 if Ai_e > 1e-6 else np.nan
print(f"  M_total ratio = {ratio_M_tot:.4f}  (vs r21={R21_EXP})")
print(f"  A_inf^4 ratio = {ratio_A4:.4f}  (vs r21={R21_EXP})")
print(f"  Odch. M_tot:  {abs(ratio_M_tot-R21_EXP)/R21_EXP*100:.3f}%")
print(f"  Odch. A^4:    {abs(ratio_A4-R21_EXP)/R21_EXP*100:.6f}%")
print()
print(f"  M_total/A^4 (elektron): {Mt_e/Ai_e**4:.5f}")
print(f"  M_total/A^4 (mion):     {Mt_mu/Ai_mu**4:.5f}")
print(f"  Stosunek wspolczynnikow: {(Mt_mu/Ai_mu**4)/(Mt_e/Ai_e**4):.5f}  (oczekiwane: 1.0000)")
print()

# ============================================================
# CZESC 4: Pelen skan M_total/A^4 dla obu sektorow
# ============================================================
print("-"*68)
print("--- Czesc 4: M_total/A^4 dla obu sektorow ---")
print()

g0_scan4_n0 = np.linspace(1.10, 1.55, 10)
g0_scan4_n1 = np.linspace(1.65, 2.25, 10)
g0_all = list(g0_scan4_n0) + list(g0_scan4_n1)

print(f"  {'g0':>8}  {'M_tot':>12}  {'A_inf':>10}  {'M/A^4':>12}  {'sector':>8}")
print("  " + "-"*55)
all_MA4 = []
all_Ai  = []
all_Mt  = []
for g0 in g0_all:
    Mk, Mp, Mt, Ai = compute_energies(g0)
    A4  = Ai**4 if Ai > 1e-6 else 1e-30
    MA4 = Mt/A4
    nb  = "n=0" if g0 <= 1.61 else "n=1"
    mark = " <g0*" if abs(g0-G0_STAR) < 0.02 else ""
    all_MA4.append(MA4); all_Ai.append(Ai); all_Mt.append(Mt)
    print(f"  {g0:8.4f}  {Mt:12.4f}  {Ai:10.6f}  {MA4:12.4f}  {nb:>8}{mark}")

# Overall power law
k_all, c_all, R2_all = power_fit(np.array(all_Ai), np.array(all_Mt))
print()
print(f"  Power law (obie sektory): M_tot = {c_all:.4f} * A_inf^{k_all:.4f}  R^2={R2_all:.5f}")
print()

# ============================================================
# TESTY
# ============================================================
print("="*68)
print("TESTY")
print("="*68)

# T1: M_kin ∝ A^k1, k1 ∈ [1.8, 2.2]
T1 = not np.isnan(k_kin) and 1.8 <= k_kin <= 2.2
print(f"  T1: M_kin ∝ A^k, k in [1.8,2.2]:  {'PASS' if T1 else 'FAIL'}  (k={k_kin:.4f}, R^2={R2_kin:.4f})")

# T2: M_pot negative and ∝ A^k2, k2 ∈ [1.8, 2.2]
T2a = not np.isnan(k_pot) and 1.8 <= k_pot <= 2.2
T2b = Mp_arr.mean() < 0  # M_pot ujemna
T2  = T2a and T2b
print(f"  T2: M_pot < 0 and ∝ A^k, k in [1.8,2.2]: {'PASS' if T2 else 'FAIL'}  "
      f"(k={k_pot:.4f}, mean_Mpot={Mp_arr.mean():.4f})")

# T3: M_total ∝ A^k3, k3 ∈ [3.5, 4.5]
T3 = not np.isnan(k_tot) and 3.5 <= k_tot <= 4.5
print(f"  T3: M_tot ∝ A^k, k in [3.5,4.5]:  {'PASS' if T3 else 'FAIL'}  (k={k_tot:.4f}, R^2={R2_tot:.4f})")

# T4: |M_tot ratio - 207| < 10%
T4 = not np.isnan(ratio_M_tot) and abs(ratio_M_tot - R21_EXP)/R21_EXP < 0.10
print(f"  T4: |M_tot_ratio - 207| < 10%:     {'PASS' if T4 else 'FAIL'}  "
      f"(ratio = {ratio_M_tot:.4f})")

# T5: M_tot/A^4 stabilne (CV < 30%)
cv_all = np.std(np.array(all_MA4)[~np.isnan(all_MA4)]) / np.mean(np.array(all_MA4)[~np.isnan(all_MA4)]) * 100
T5 = cv_all < 30
print(f"  T5: M_tot/A^4 CV < 30%:            {'PASS' if T5 else 'FAIL'}  (CV={cv_all:.2f}%)")

n_pass = sum([T1, T2, T3, T4, T5])
print(f"\nWYNIK: {n_pass}/5 testow przeszlo")

print()
print("="*68)
print("WNIOSEK (EX65)")
print("="*68)
print(f"  M_kin ∝ A^{k_kin:.3f}  (oczekiwane: A^2)")
print(f"  |M_pot| ∝ A^{k_pot:.3f}  (oczekiwane: A^2, ujemna)")
print(f"  M_tot ∝ A^{k_tot:.3f}  (oczekiwane: A^4)")
print(f"  M_tot(phi*g0*)/M_tot(g0*) = {ratio_M_tot:.4f}  (vs r21=206.768)")
print(f"  M_tot/A^4 (elektron) = {Mt_e/Ai_e**4:.4f}")
print(f"  M_tot/A^4 (mion)     = {Mt_mu/Ai_mu**4:.4f}")
print(f"  Ogon supresja (r>10, g0*): M_tail/M_tot = {(Mk_t+Mp_t)/(Mk_f+Mp_f)*100:.2f}%")
print("="*68)
