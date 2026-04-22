"""
EX81: S(R) = alpha*_1(R) + alpha*_2(R) jako funkcja targetu R
STATUS: LEGACY-TRANSLATIONAL

This file continues the older `alpha*` / sum-rule selection program in legacy
variables. Treat it as historical transition material, not as a canonical
modern entry point.

Pytanie: czy wzor S = 2pi - 11/10 jest specyficzny dla R = r21 = 206.768,
czy jest ogolna wlasnoscia geometrii f(z0, alpha)?

Plan:
  1. Oblicz f(z0, alpha) dla gestej siatki alpha ~ [2.30, 2.90] (szeroko)
  2. Dla kazdego R: znajdz zera f(z0,alpha)=R przez interpolacje (bezkoszlowo)
  3. Oblicz S(R) i D(R) dla R w {180, 190, 198, 200, 203, 206.768, 210, 215, 220, 230}
  4. Sprawdz: czy S(R) jest stale = 2pi-11/10, czy zmienia sie?

Metoda:
  - Gesta siatka alpha: 50 punktow w [2.25, 2.90]
  - Dla kazdego alpha: f_val = f(z0(alpha), alpha)
  - Potem interpolacja f_val vs alpha, szukanie zer dla kazdego R
  - Bez zadnych dodatkowych ODE po wyznaczeniu siatki f(z0,alpha)
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, curve_fit
from scipy.interpolate import interp1d
from multiprocessing import Pool, cpu_count
import warnings
warnings.filterwarnings('ignore')

PHI       = (1.0 + np.sqrt(5.0)) / 2.0
R21_EXP   = 206.768
G_OFF     = 0.005
R_MAX     = 100.0
WIN_LIST  = [16, 22, 28, 36, 46, 58, 72]
WIN_WIDTH = 14.0
N_CORES   = min(16, cpu_count())

def _integrate(g0, ak, maxb=8):
    gb = np.exp(-1.0/(2.0*ak)) + G_OFF
    def rhs(r, y):
        g, gp = y; g = max(g, gb+1e-7)
        fg = 1.0 + 2.0*ak*np.log(g)
        if abs(fg) < 1e-10: return [gp, 0.0]
        dr = g**2*(1.0-g); cr = (ak/g)*gp**2
        if r < 1e-10: return [gp, (dr-cr)/(3.0*fg)]
        return [gp, (dr-cr-fg*2.0*gp/r)/fg]
    def ev(r, y): return y[0]-gb
    ev.terminal = True; ev.direction = -1
    y0 = [g0, 0.0]; r0 = 1e-10; ra, ga = [], []
    for _ in range(maxb+1):
        sol = solve_ivp(rhs, (r0, R_MAX), y0, events=ev,
                        dense_output=True, rtol=1e-10, atol=1e-12, max_step=0.04)
        ra.append(sol.t); ga.append(sol.y[0])
        if sol.status == 1 and len(sol.t_events[0]) > 0:
            rh = sol.t_events[0][-1]; st = sol.sol(rh)
            y0 = [st[0], -st[1]]; r0 = rh
        else: break
    r = np.concatenate(ra); g = np.concatenate(ga)
    return r[np.argsort(r)], g[np.argsort(r)]

def _fit_win(r, g, rL, rR):
    mask = (r >= rL) & (r <= rR)
    if mask.sum() < 20: return 0., 0., 0.
    rf = r[mask]; df = (g[mask]-1.0)*rf
    M = np.column_stack([np.cos(rf), np.sin(rf)])
    bc = np.linalg.lstsq(M, df, rcond=None)[0]
    return float(np.sqrt(bc[0]**2+bc[1]**2)), float(bc[0]), float(bc[1])

def _A_inf(g0, ak):
    r, g = _integrate(g0, ak)
    Av, rv = [], []
    for rL in WIN_LIST:
        if rL + WIN_WIDTH > r[-1]: break
        A, _, _ = _fit_win(r, g, rL, rL+WIN_WIDTH)
        if A > 0.002: Av.append(A); rv.append(float(rL))
    if not Av: return 0.
    if len(Av) < 3: return Av[-1]
    try:
        p, _ = curve_fit(lambda x, ai, a: ai*(1+a/x), rv, Av, p0=[Av[-1], 0.], maxfev=2000)
        return float(p[0])
    except: return float(Av[-1])

def _B_coeff(g0, ak):
    r, g = _integrate(g0, ak)
    _, B, _ = _fit_win(r, g, WIN_LIST[2], WIN_LIST[2]+WIN_WIDTH)
    return B

def _find_z0(ak, lo=1.05, hi=2.0):
    gs = np.linspace(lo, hi, 40)
    Bs = [_B_coeff(g, ak) for g in gs]
    for i in range(len(Bs)-1):
        if Bs[i]*Bs[i+1] < 0:
            try: return brentq(lambda x: _B_coeff(x, ak), gs[i], gs[i+1], xtol=1e-8)
            except: pass
    return None

def _f_ratio(g0, ak):
    Ae = _A_inf(g0, ak); Am = _A_inf(PHI*g0, ak)
    return (Am/Ae)**4 if Ae > 1e-6 else np.nan

def fz0(ak):
    """f(z0(alpha), alpha) dla danego alpha."""
    z0 = _find_z0(ak)
    if z0 is None: return np.nan
    return _f_ratio(z0, ak)

def worker(ak):
    return ak, fz0(ak)

def find_zeros_by_interp(a_arr, f_arr, R_target):
    """
    Znajdz zera f(z0,alpha) = R przez interpolacje.
    Zwraca liste zer (wartosci alpha).
    """
    valid = [(a, f) for a, f in zip(a_arr, f_arr) if not np.isnan(f)]
    if len(valid) < 3: return []
    av = np.array([x[0] for x in valid])
    fv = np.array([x[1] for x in valid])

    zeros = []
    for i in range(len(av)-1):
        if (fv[i]-R_target)*(fv[i+1]-R_target) < 0:
            # Interpolacja liniowa
            a_zero = av[i] + (R_target-fv[i])*(av[i+1]-av[i])/(fv[i+1]-fv[i])
            zeros.append(float(a_zero))
    return zeros

if __name__ == '__main__':
    print("=" * 70)
    print("EX81: S(R) = alpha*_1(R) + alpha*_2(R) jako funkcja R")
    print("=" * 70)
    print(f"  Uzywam {N_CORES} rdzeni")
    print()

    # Gesta siatka alpha: 50 punktow w [2.30, 2.90]
    alpha_arr = np.linspace(2.30, 2.90, 50)
    print(f"  Obliczam f(z0, alpha) dla {len(alpha_arr)} punktow alpha...")
    print(f"  alpha in [{alpha_arr[0]:.3f}, {alpha_arr[-1]:.3f}]")
    print()

    with Pool(N_CORES) as pool:
        results = pool.map(worker, alpha_arr)

    avals = np.array([r[0] for r in results])
    fvals = np.array([r[1] for r in results])

    # Wypisz siatke
    print(f"  {'alpha':>7}  {'f(z0)':>10}  {'f-r21':>9}")
    print(f"  {'-'*7}  {'-'*10}  {'-'*9}")
    for a, f in zip(avals, fvals):
        if not np.isnan(f):
            diff_str = f"{f-R21_EXP:+.3f}"
            mark = " <<<" if abs(f-R21_EXP) < 1 else ""
            print(f"  {a:7.4f}  {f:10.4f}  {diff_str:>9}{mark}")
    print()

    # Minimum f(z0)
    valid_mask = ~np.isnan(fvals)
    if valid_mask.sum() > 3:
        av = avals[valid_mask]
        fv = fvals[valid_mask]
        idx_min = np.argmin(fv)
        f_min_approx = fv[idx_min]
        a_min_approx = av[idx_min]
        print(f"  Minimum f(z0) ~ {f_min_approx:.4f} przy alpha ~ {a_min_approx:.4f}")
        print()

    # --- Skan R ---
    R_targets = [180, 190, f_min_approx + 0.5, f_min_approx + 1.0,
                 f_min_approx + 3.0, f_min_approx + 5.0,
                 200, 203, 206.768, 210, 215, 220, 230, 250, 280]
    R_targets = sorted(set([round(r, 4) for r in R_targets]))

    print("=" * 70)
    print("S(R) = alpha*_1(R) + alpha*_2(R) dla roznych R")
    print("=" * 70)
    print(f"  {'R':>8}  {'a1(R)':>10}  {'a2(R)':>10}  {'S(R)':>12}  "
          f"{'S-2pi+11/10':>12}  {'D(R)':>9}")
    print(f"  {'-'*8}  {'-'*10}  {'-'*10}  {'-'*12}  {'-'*12}  {'-'*9}")

    ref_S = 2*np.pi - 1.1
    results_scan = []

    for R in R_targets:
        zeros = find_zeros_by_interp(avals, fvals, R)
        if len(zeros) == 2:
            a1r, a2r = zeros[0], zeros[1]
            Sr = a1r + a2r
            Dr = a2r - a1r
            diff_S = Sr - ref_S
            results_scan.append((R, a1r, a2r, Sr, Dr, diff_S))
            mark = " <<< r21" if abs(R-R21_EXP) < 0.5 else ""
            print(f"  {R:8.3f}  {a1r:10.6f}  {a2r:10.6f}  {Sr:12.8f}  "
                  f"{diff_S:+12.8f}  {Dr:9.6f}{mark}")
        elif len(zeros) == 1:
            print(f"  {R:8.3f}  tylko 1 zero: {zeros[0]:.6f}")
        else:
            print(f"  {R:8.3f}  brak zer")

    print()

    if results_scan:
        Rs = np.array([x[0] for x in results_scan])
        Ss = np.array([x[3] for x in results_scan])
        Ds = np.array([x[4] for x in results_scan])
        dS = np.array([x[5] for x in results_scan])

        # Pochodna dS/dR w okolicach r21
        # Znajdz punkty blisko r21
        idx_r21 = np.argmin(np.abs(Rs - R21_EXP))
        print(f"  Przy R=r21: S = {Ss[idx_r21]:.8f}")
        print(f"  Odchylenie S(r21) - (2pi-11/10) = {dS[idx_r21]:+.4e}")
        print()

        # Sprawdz czy S jest stala
        print(f"  Zakres S(R) dla R w [{Rs[0]:.1f}, {Rs[-1]:.1f}]:")
        print(f"    S_min = {np.min(Ss):.8f} przy R={Rs[np.argmin(Ss)]:.3f}")
        print(f"    S_max = {np.max(Ss):.8f} przy R={Rs[np.argmax(Ss)]:.3f}")
        print(f"    Zmiana S na caly zakres: {np.max(Ss)-np.min(Ss):.6f}")
        print()

        # Gradient dS/dR
        if len(Rs) >= 3:
            dS_dR = np.gradient(Ss, Rs)
            print(f"  dS/dR (gradientowe aproks.):")
            for i, (Rr, ds) in enumerate(zip(Rs, dS_dR)):
                if not np.isnan(ds):
                    print(f"    R={Rr:.3f}: dS/dR = {ds:+.6f}")

        # Czy S jest liniowe w R?
        print()
        print(f"  Interpolacja S(R):")
        if len(Rs) >= 4:
            from numpy.polynomial import polynomial as P_poly
            coeffs = np.polyfit(Rs, Ss, 2)
            print(f"    S(R) ~ {coeffs[0]:.4e}*R^2 + {coeffs[1]:.4e}*R + {coeffs[2]:.4f}")
            # Ewaluuj w r21
            S_interp = np.polyval(coeffs, R21_EXP)
            print(f"    S(r21) z dopasowania: {S_interp:.8f}")
            print(f"    2pi-11/10:             {ref_S:.8f}")

        # Sprawdz czy D(R) ma wzor
        print()
        print(f"  D(R) = a2(R) - a1(R):")
        for R, a1r, a2r, Sr, Dr, _ in results_scan:
            if abs(R - R21_EXP) < 20:
                print(f"    R={R:.3f}: D={Dr:.8f}")

    # TESTY
    print()
    print("=" * 70)
    print("TESTY")
    print("=" * 70)

    T1 = any(abs(x[0]-R21_EXP)<0.01 for x in results_scan)
    print(f"  T1: f(z0) osiaga r21 (2 zera): {'PASS' if T1 else 'FAIL'}")

    T2 = False
    if results_scan:
        # Czy S silnie zmienia sie z R?
        Ss_all = np.array([x[3] for x in results_scan])
        T2 = np.max(Ss_all) - np.min(Ss_all) > 0.01
        print(f"  T2: S(R) zmienia sie z R (>0.01): {'PASS' if T2 else 'FAIL'}  "
              f"(zakres={np.max(Ss_all)-np.min(Ss_all):.4f})")

    # T3: S(r21) jest blisko 2pi-11/10 (< 1 ppm)
    T3 = False
    for R, a1r, a2r, Sr, Dr, dS_val in results_scan:
        if abs(R-R21_EXP) < 0.5:
            T3 = abs(dS_val)/Sr < 1e-5
            print(f"  T3: S(r21) ~ 2pi-11/10 (<10ppm): {'PASS' if T3 else 'FAIL'}  "
                  f"(diff={dS_val:.3e})")
            break

    T4 = len([x for x in results_scan if x[0] < f_min_approx]) == 0
    print(f"  T4: Brak zer dla R < f_min: {'PASS' if T4 else 'FAIL'}")

    if results_scan:
        Ss_all = np.array([x[3] for x in results_scan])
        T5 = abs(Ss_all[0] - Ss_all[-1]) > 0.05
        print(f"  T5: S nie jest stala (zmienia sie >0.05): {'PASS' if T5 else 'FAIL'}  "
              f"(delta_S={abs(Ss_all[0]-Ss_all[-1]):.4f})")
    else:
        T5 = False
        print(f"  T5: T5 - brak danych")

    n_pass = sum([T1, T2, T3, T4, T5])
    print(f"\nWYNIK: {n_pass}/5")

    print()
    print("=" * 70)
    print("WNIOSEK (EX81)")
    print("=" * 70)
    if results_scan:
        Ss_all = np.array([x[3] for x in results_scan])
        Rs_all = np.array([x[0] for x in results_scan])
        print(f"  S(R) NIE jest stala: zmienia sie od {np.min(Ss_all):.6f} do {np.max(Ss_all):.6f}")
        print(f"  => Wzor S=2pi-11/10 jest specyficzny dla R=r21")
        print()
        r21_res = [(R,Sr,dS_v) for R,_,_,Sr,_,dS_v in results_scan if abs(R-R21_EXP)<1]
        if r21_res:
            R, Sr, dS_v = r21_res[0]
            print(f"  S(r21={R21_EXP}) = {Sr:.8f}")
            print(f"  2pi-11/10        = {ref_S:.8f}")
            print(f"  Roznica          = {dS_v:+.3e}")
        print()

        # Szukaj czy S(R) = jakas funkcja R, ktora dla R=r21 daje 2pi-11/10
        print(f"  Szukam wzoru S(R) zawierajacego R=r21...")
        # Prosta liniowa zaleznosc w okolicach r21
        r21_nearby = [(R,Sr) for R,_,_,Sr,_,_ in results_scan if 195 <= R <= 220]
        if len(r21_nearby) >= 3:
            Rs_near = np.array([x[0] for x in r21_nearby])
            Ss_near = np.array([x[1] for x in r21_nearby])
            slope, intercept = np.polyfit(Rs_near, Ss_near, 1)
            print(f"    Liniowe dopasowanie S(R) ~ {slope:.6e}*R + {intercept:.6f}")
            print(f"    dS/dR ~ {slope:.4e} (przy R=r21)")
            S_at_r21_fit = slope*R21_EXP + intercept
            print(f"    S(r21) z dopasowania = {S_at_r21_fit:.8f}")
            print(f"    2pi-11/10            = {ref_S:.8f}")
    print("=" * 70)
