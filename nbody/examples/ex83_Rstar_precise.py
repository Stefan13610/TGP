"""
EX83 v3: PRECYZYJNE alpha*(R_CODATA) przez brentq w alpha
STATUS: LEGACY-TRANSLATIONAL

This file remains in the older `alpha*` / `R21` matching chain using pre-sync
selection language. It should be read as legacy transition material, not as a
canonical current reference.

Problem z v1/v2: fz0(alpha) ma systematyczne przesunięcie ~0.250 w R
(błąd amplitudy A_inf). Pochodne sa prawidlowe, wartosci bezwzgledne - nie.

Rozwiazanie: dla R = r21_CODATA znajdz alpha*(R) przez brentq w alpha,
uzywajac funkcji delta(alpha, R) = g0star(R, alpha) - z0(alpha).
Uzywamy watskich przedziałow [A_i - 5e-3, A_i + 5e-3] jako nawias.

Wynik: precyzyjne S(r21_CODATA) = alpha*_1(CODATA) + alpha*_2(CODATA)
Porownanie z S_formula = 2pi - 11/10.
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, curve_fit
from multiprocessing import Pool, cpu_count
import warnings
warnings.filterwarnings('ignore')

# ===== STALE =====
PHI        = (1.0 + np.sqrt(5.0)) / 2.0
G_OFF      = 0.005
R_MAX      = 100.0
WIN_LIST   = [16, 22, 28, 36, 46, 58, 72]
WIN_WIDTH  = 14.0
N_CORES    = min(16, cpu_count())

A1_KNOWN   = 2.44143051
A2_KNOWN   = 2.74175411
R21_CODE   = 206.768
R21_CODATA = 206.7682830
R21_UNC    = 0.0000046
S_R21_EX80 = 5.18318462
S_FORMULA  = 2*np.pi - 1.1
DSDR_EX82  = 0.007967

# ===== ODE =====

def _integrate(g0, ak):
    gb = np.exp(-1.0/(2.0*ak)) + G_OFF
    def rhs(r, y):
        g, gp = y; g = max(g, gb + 1e-7)
        fg = 1.0 + 2.0*ak*np.log(g)
        if abs(fg) < 1e-10: return [gp, 0.0]
        dr = g**2*(1.0 - g); cr = (ak/g)*gp**2
        if r < 1e-10: return [gp, (dr - cr)/(3.0*fg)]
        return [gp, (dr - cr - fg*2.0*gp/r)/fg]
    def ev(r, y): return y[0] - gb
    ev.terminal = True; ev.direction = -1
    y0 = [g0, 0.0]; r0 = 1e-10; ra, ga = [], []
    for _ in range(9):
        sol = solve_ivp(rhs, (r0, R_MAX), y0, events=ev,
                        dense_output=True, rtol=1e-10, atol=1e-12,
                        max_step=0.04)
        ra.append(sol.t); ga.append(sol.y[0])
        if sol.status == 1 and len(sol.t_events[0]) > 0:
            rh = sol.t_events[0][-1]; st = sol.sol(rh)
            y0 = [st[0], -st[1]]; r0 = rh
        else: break
    r = np.concatenate(ra); g = np.concatenate(ga)
    idx = np.argsort(r)
    return r[idx], g[idx]

def _fit_win(r, g, rL, rR):
    mask = (r >= rL) & (r <= rR)
    if mask.sum() < 20: return 0., 0., 0.
    rf = r[mask]; df = (g[mask] - 1.0)*rf
    M = np.column_stack([np.cos(rf), np.sin(rf)])
    bc = np.linalg.lstsq(M, df, rcond=None)[0]
    return float(np.sqrt(bc[0]**2 + bc[1]**2)), float(bc[0]), float(bc[1])

def _A_inf(g0, ak):
    r, g = _integrate(g0, ak)
    Av, rv = [], []
    for rL in WIN_LIST:
        if rL + WIN_WIDTH > r[-1]: break
        A, _, _ = _fit_win(r, g, rL, rL + WIN_WIDTH)
        if A > 0.002: Av.append(A); rv.append(float(rL))
    if not Av: return 0.
    if len(Av) < 3: return Av[-1]
    try:
        p, _ = curve_fit(lambda x, ai, a: ai*(1 + a/x), rv, Av,
                         p0=[Av[-1], 0.], maxfev=2000)
        return float(p[0])
    except: return float(Av[-1])

def _B_coeff(g0, ak):
    r, g = _integrate(g0, ak)
    _, B, _ = _fit_win(r, g, WIN_LIST[2], WIN_LIST[2] + WIN_WIDTH)
    return B

def _find_z0(ak, lo=1.05, hi=2.0):
    gs = np.linspace(lo, hi, 40)
    Bs = [_B_coeff(g, ak) for g in gs]
    for i in range(len(Bs)-1):
        if Bs[i]*Bs[i+1] < 0:
            try:
                return brentq(lambda x: _B_coeff(x, ak), gs[i], gs[i+1],
                              xtol=1e-9)
            except: pass
    return None

def _find_g0star(ak, R_target, lo=1.05, hi=2.0):
    """Znajdz g0 takie ze (A(phi*g0)/A(g0))^4 = R_target."""
    def ratio(g0):
        Ae = _A_inf(g0, ak)
        Am = _A_inf(PHI * g0, ak)
        return (Am/Ae)**4 - R_target if Ae > 1e-6 else np.nan
    gs = np.linspace(lo, hi, 25)
    rs = []
    for g in gs:
        try: rs.append(ratio(g))
        except: rs.append(np.nan)
    for i in range(len(rs)-1):
        if not (np.isnan(rs[i]) or np.isnan(rs[i+1])):
            if rs[i]*rs[i+1] < 0:
                try:
                    return brentq(ratio, gs[i], gs[i+1], xtol=1e-9)
                except: pass
    return None

def delta_alpha(ak, R_target):
    """delta = g0star(R, alpha) - z0(alpha). Zero przy alpha*(R)."""
    z0  = _find_z0(ak)
    g0s = _find_g0star(ak, R_target)
    if z0 is None or g0s is None: return np.nan
    return g0s - z0

def find_alpha_star_brentq(R_target, a_lo, a_hi):
    """
    Znajdz alpha*(R) przez brentq w [a_lo, a_hi].
    Najpierw sprawdz ze delta zmienia znak na krawędziach.
    """
    d_lo = delta_alpha(a_lo, R_target)
    d_hi = delta_alpha(a_hi, R_target)
    if np.isnan(d_lo) or np.isnan(d_hi): return None
    if d_lo * d_hi > 0:
        a_mid = (a_lo + a_hi) / 2
        d_mid = delta_alpha(a_mid, R_target)
        if np.isnan(d_mid): return None
        if d_lo * d_mid < 0: a_hi = a_mid
        elif d_mid * d_hi < 0: a_lo = a_mid
        else: return None
    try:
        return brentq(lambda a: delta_alpha(a, R_target), a_lo, a_hi, xtol=1e-9)
    except: return None

def worker_alpha1(R_target):
    """Wyznacz alpha*_1(R) przez brentq. Na poziomie modulu dla multiprocessing."""
    bracket = 5e-3
    return R_target, 1, find_alpha_star_brentq(R_target, A1_KNOWN - bracket, A1_KNOWN + bracket)

def worker_alpha2(R_target):
    """Wyznacz alpha*_2(R) przez brentq. Na poziomie modulu dla multiprocessing."""
    bracket = 5e-3
    return R_target, 2, find_alpha_star_brentq(R_target, A2_KNOWN - bracket, A2_KNOWN + bracket)

def worker(args):
    """Dispatcher na poziomie modulu (wymagane przez Windows multiprocessing)."""
    R_target, branch = args
    if branch == 1: return worker_alpha1(R_target)
    else:           return worker_alpha2(R_target)


if __name__ == '__main__':
    print("="*68)
    print("EX83 v3: PRECYZYJNE alpha*(R_CODATA) PRZEZ BRENTQ W ALPHA")
    print("="*68)
    print(f"  S_formula  = {S_FORMULA:.10f}")
    print(f"  S(R21_CODE)= {S_R21_EX80:.10f}  (ex80 brentq)")
    print(f"  r21_CODATA = {R21_CODATA:.7f} +/- {R21_UNC:.7f}")
    print(f"  dS/dR_ex82 = {DSDR_EX82:.6f}")
    print(f"  N_CORES    = {N_CORES}")
    print()

    # R_SCAN: R21_CODE (kontrola) + kilka punktow blisko CODATA
    R_SCAN = np.array([
        R21_CODE,
        R21_CODATA - 0.010,
        R21_CODATA,
        R21_CODATA + 0.010,
    ])

    # Dla kazdego R wyznacz alpha*_1 i alpha*_2 rownolegle
    tasks = [(R, 1) for R in R_SCAN] + [(R, 2) for R in R_SCAN]

    print(f"  Obliczam alpha*(R) dla {len(R_SCAN)} wartosci R (2 galezie)...")
    print(f"  Laczna liczba brentq-alpha: {len(tasks)} (kazdy ~80 ODE ewal)")
    print()

    with Pool(N_CORES) as pool:
        results = pool.map(worker, tasks)

    # Zbierz wyniki
    alpha1_dict = {}; alpha2_dict = {}
    for R, branch, a_star in results:
        if branch == 1: alpha1_dict[R] = a_star
        else:           alpha2_dict[R] = a_star

    print(f"  {'R':>12}  {'alpha*_1':>13}  {'alpha*_2':>13}  {'S(R)':>14}  {'S-formula':>12}")
    print("-"*72)

    R_vals, S_vals = [], []
    for R in sorted(R_SCAN):
        a1 = alpha1_dict.get(R)
        a2 = alpha2_dict.get(R)
        if a1 is not None and a2 is not None:
            S = a1 + a2
            diff = S - S_FORMULA
            marker = " <-CODATA" if abs(R - R21_CODATA) < 1e-6 else (
                     " <-CODE"   if abs(R - R21_CODE)   < 1e-6 else "")
            print(f"  {R:12.7f}  {a1:13.9f}  {a2:13.9f}  {S:14.10f}  {diff:+.4e}{marker}")
            R_vals.append(R); S_vals.append(S)
        else:
            print(f"  {R:12.7f}  {'FAIL':>13}  {'FAIL':>13}  {'nan':>14}  {'nan':>12}")

    R_vals = np.array(R_vals)
    S_vals = np.array(S_vals)

    print()
    print("="*68)
    print("WYNIKI I POROWNANIE Z EX80/EX82")
    print("="*68)

    # S przy R21_CODE (kontrola ex80)
    S_code = None
    for R, S in zip(R_vals, S_vals):
        if abs(R - R21_CODE) < 1e-6: S_code = S; break

    # S przy R21_CODATA
    S_codata = None
    for R, S in zip(R_vals, S_vals):
        if abs(R - R21_CODATA) < 1e-6: S_codata = S; break

    if S_code is not None:
        print(f"  S(R21_CODE)  = {S_code:.10f}")
        print(f"  S(R21_CODE) ex80 = {S_R21_EX80:.10f}")
        print(f"  Roznica      = {S_code - S_R21_EX80:+.4e}")

    if S_codata is not None:
        print(f"\n  S(R21_CODATA) = {S_codata:.10f}")
        print(f"  S_formula     = {S_FORMULA:.10f}")
        delta_S_codata = S_codata - S_FORMULA
        ppm_codata = delta_S_codata / S_FORMULA * 1e6
        print(f"  Roznica       = {delta_S_codata:+.4e}  ({ppm_codata:+.3f} ppm)")

    # R* z interpolacji i liniowej aproksymacji
    print()
    if len(R_vals) >= 3:
        R0 = R21_CODATA
        dR_arr = R_vals - R0
        dS_arr = S_vals - S_FORMULA

        c_lin = np.polyfit(dR_arr, dS_arr, 1)
        dSdR  = c_lin[0]
        S_at_R0 = np.polyval(c_lin, 0.) + S_FORMULA
        Rstar   = R0 - np.polyval(c_lin, 0.) / c_lin[0]

        print(f"  Liniowe: dS/dR = {dSdR:.6f}  (ex82: {DSDR_EX82:.6f})")
        print(f"           S(CODATA) aproks. = {S_at_R0:.10f}")
        print(f"           R*_linear = {Rstar:.9f}")
        print(f"           R*-CODATA = {(Rstar-R21_CODATA)*1e6:+.3f} ppm")

        # Bezposrednia interpolacja
        Rstar_interp = None
        for i in range(len(dS_arr)-1):
            if dS_arr[i]*dS_arr[i+1] < 0:
                Rstar_interp = R_vals[i] + (-dS_arr[i])*(R_vals[i+1]-R_vals[i])/(dS_arr[i+1]-dS_arr[i])
                break
        if Rstar_interp is not None:
            print(f"\n  Interpolacja: R*_interp = {Rstar_interp:.9f}")
            print(f"                R*-CODATA = {(Rstar_interp-R21_CODATA)*1e6:+.3f} ppm")
        else:
            print(f"\n  Interpolacja: brak (S nie przecina S_formula w siatce)")
            # Extrapolacja liniowa z R21_CODE jesli znamy S(CODE)
            if S_code is not None:
                delta_S_code = S_code - S_FORMULA
                Rstar_from_code = R21_CODE - delta_S_code / dSdR
                print(f"  Ekstrapolacja z R21_CODE: R* = {Rstar_from_code:.9f}")
                print(f"  R*-CODATA = {(Rstar_from_code-R21_CODATA)*1e6:+.3f} ppm")

    print()
    print("="*68)
    print("TESTY")
    print("="*68)

    # T1: S(R21_CODE) ~ S_R21_EX80 (1e-5 tolerancja)
    T1 = S_code is not None and abs(S_code - S_R21_EX80) < 1e-5
    print(f"  T1: S(R_CODE) ~ {S_R21_EX80} (1e-5 tol): {'PASS' if T1 else 'FAIL'}"
          + (f"  diff={S_code-S_R21_EX80:.2e}" if S_code else ""))

    # T2: dS/dR blisko ex82 (10% tol)
    T2 = 'dSdR' in dir() and abs(dSdR - DSDR_EX82) < 0.1*DSDR_EX82
    print(f"  T2: dS/dR ~ {DSDR_EX82} (10% tol): {'PASS' if T2 else 'FAIL'}"
          + (f"  dS/dR={dSdR:.6f}" if 'dSdR' in dir() else ""))

    # T3: |R*-CODATA| < 10 ppm
    Rstar_best = Rstar_interp if Rstar_interp is not None else (
                 Rstar if 'Rstar' in dir() else np.nan)
    T3 = not np.isnan(Rstar_best) and abs(Rstar_best-R21_CODATA) < 10e-6*R21_CODATA
    print(f"  T3: |R*-CODATA| < 10 ppm: {'PASS' if T3 else 'FAIL'}"
          + (f"  {abs(Rstar_best-R21_CODATA)/R21_CODATA*1e6:.3f} ppm" if not np.isnan(Rstar_best) else ""))

    # T4: S monotonicznie
    T4 = len(S_vals) > 1 and all(S_vals[i] <= S_vals[i+1] for i in range(len(S_vals)-1))
    print(f"  T4: S(R) monotonicznie: {'PASS' if T4 else 'FAIL'}")

    # T5: S(CODATA) blisko formula (5 ppm)
    T5 = S_codata is not None and abs(S_codata - S_FORMULA) / S_FORMULA * 1e6 < 5.0
    print(f"  T5: |S(CODATA)-formula| < 5 ppm: {'PASS' if T5 else 'FAIL'}"
          + (f"  {(S_codata-S_FORMULA)/S_FORMULA*1e6:+.3f} ppm" if S_codata else ""))

    n_pass = sum([T1, T2, T3, T4, T5])
    print(f"\nWYNIK: {n_pass}/5")

    print()
    print("="*68)
    print("WNIOSEK (EX83 v3)")
    print("="*68)
    if S_codata is not None:
        print(f"  S(r21_CODATA)  = {S_codata:.10f}")
        print(f"  2pi - 11/10    = {S_FORMULA:.10f}")
        print(f"  Roznica        = {S_codata-S_FORMULA:+.4e}  ({(S_codata-S_FORMULA)/S_FORMULA*1e6:+.3f} ppm)")
    if not np.isnan(Rstar_best):
        print(f"  R* (TGP)       = {Rstar_best:.9f}")
        print(f"  r21_CODATA     = {R21_CODATA:.7f} +/- {R21_UNC:.7f}")
        print(f"  R* - CODATA    = {(Rstar_best-R21_CODATA)*1e6:+.3f} ppm")
    print("="*68)
