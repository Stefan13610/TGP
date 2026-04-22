"""
EX82: PRECYZYJNE dS/dR i PREDYKCJA r21 Z WARUNKU S = 2pi - 11/10
STATUS: LEGACY-TRANSLATIONAL

This script extends the historical `alpha*` / sum-rule track built on older
selection variables. Keep it as legacy exploratory context rather than a
canonical synchronized `nbody` example.

Z ex80: S(r21=206.768) = 5.18318462, 2pi-11/10 = 5.18318531, roznica = -6.87e-7
Z ex81: dS/dR ~ 8.87e-3 (gruba siatka)
=> R* = r21 + 6.87e-7/8.87e-3 = r21 + 7.75e-5 ~ 206.7681

Teraz: precyzyjne dS/dR przez lokalne roznickowanie f(z0,alpha) blisko alpha*_1 i alpha*_2.
Potrzeba: tylko ~10 ODE blisko dwoch przejsc f=r21.

Plan:
  1. Compute f(z0, alpha) dla 5 punktow blisko alpha*_1 (zakres +-0.01)
  2. Compute f(z0, alpha) dla 5 punktow blisko alpha*_2 (zakres +-0.01)
  3. Dopasuj polinom 2-go stopnia do kazdej grupy
  4. Oblicz df/dalpha i dS/dR = 1/|df/da1| + 1/|df/da2| (ze znakami)
  5. R* = r21 + (2pi-11/10 - S(r21)) / dS/dR
  6. Porownaj R* z eksperymentalnym r21 = 206.7682830 (CODATA 2018)
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, curve_fit
from multiprocessing import Pool, cpu_count
import warnings
warnings.filterwarnings('ignore')

PHI       = (1.0 + np.sqrt(5.0)) / 2.0
R21_EXP   = 206.768         # wartosc uzywana w obliczeniach
R21_CODATA = 206.7682830    # CODATA 2018: mmu/me
S_R21     = 5.18318462      # S(r21) z ex80 brentq xtol=1e-8
S_FORMULA = 2*np.pi - 1.1   # 2pi - 11/10
G_OFF     = 0.005
R_MAX     = 100.0
WIN_LIST  = [16, 22, 28, 36, 46, 58, 72]
WIN_WIDTH = 14.0
N_CORES   = min(16, cpu_count())

# Znane alpha* z ex79/80
A1_KNOWN = 2.44143051
A2_KNOWN = 2.74175411

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
    z0 = _find_z0(ak)
    if z0 is None: return np.nan
    return _f_ratio(z0, ak)

def worker(ak):
    return ak, fz0(ak)

if __name__ == '__main__':
    print("=" * 68)
    print("EX82: PRECYZYJNE dS/dR i PREDYKCJA r21 Z S = 2pi - 11/10")
    print("=" * 68)
    print(f"  S(r21=206.768) = {S_R21:.10f}  (z ex80)")
    print(f"  2pi-11/10      = {S_FORMULA:.10f}")
    print(f"  Roznica dS     = {S_R21-S_FORMULA:+.4e}")
    print()
    print(f"  Uzywam {N_CORES} rdzeni")
    print()

    # --- Siatka lokalna blisko alpha*_1 ---
    # 5 punktow w [a1-0.012, a1+0.012]
    da = 0.012
    alpha_near_a1 = np.linspace(A1_KNOWN - da, A1_KNOWN + da, 5)
    alpha_near_a2 = np.linspace(A2_KNOWN - da, A2_KNOWN + da, 5)
    alpha_all = np.concatenate([alpha_near_a1, alpha_near_a2])

    print(f"  Obliczam f(z0,a) dla {len(alpha_all)} punktow blisko alpha*_1 i alpha*_2...")
    print(f"  alpha*_1 = {A1_KNOWN}, siatka: {alpha_near_a1[0]:.5f} .. {alpha_near_a1[-1]:.5f}")
    print(f"  alpha*_2 = {A2_KNOWN}, siatka: {alpha_near_a2[0]:.5f} .. {alpha_near_a2[-1]:.5f}")
    print()

    with Pool(N_CORES) as pool:
        results_all = pool.map(worker, alpha_all)

    avals_all = np.array([r[0] for r in results_all])
    fvals_all = np.array([r[1] for r in results_all])

    # Podziel na dwa regiony
    idx1 = np.where(avals_all < A1_KNOWN + 0.5*da + 0.001)[0][:5]
    idx2 = np.where(avals_all > A2_KNOWN - 0.5*da - 0.001)[0][:5]

    a1_arr = avals_all[idx1]; f1_arr = fvals_all[idx1]
    a2_arr = avals_all[idx2]; f2_arr = fvals_all[idx2]

    print("  Region alpha*_1:")
    for a, f in zip(a1_arr, f1_arr):
        if not np.isnan(f):
            print(f"    alpha={a:.6f}  f={f:.6f}  f-r21={f-R21_EXP:+.4f}")
    print()
    print("  Region alpha*_2:")
    for a, f in zip(a2_arr, f2_arr):
        if not np.isnan(f):
            print(f"    alpha={a:.6f}  f={f:.6f}  f-r21={f-R21_EXP:+.4f}")
    print()

    # --- Dopasowanie polinomu 2-go stopnia do kazdego regionu ---
    # Znajdz df/dalpha przez dopasowanie

    print("=" * 68)
    print("ANALIZA LOKALNEGO df/dalpha")
    print("=" * 68)

    results_deriv = {}

    for label, a_arr, f_arr, a_star in [('a*_1', a1_arr, f1_arr, A1_KNOWN),
                                          ('a*_2', a2_arr, f2_arr, A2_KNOWN)]:
        valid = ~np.isnan(f_arr)
        if valid.sum() >= 3:
            av = a_arr[valid]; fv = f_arr[valid]
            # Dopasuj polinom stopnia 1 i 2
            p1 = np.polyfit(av, fv, 1)
            p2 = np.polyfit(av, fv, 2)
            # Pochodna w alpha*
            df_da_lin  = p1[0]
            df_da_quad = 2*p2[0]*a_star + p2[1]  # d/da (p2[0]*a^2 + p2[1]*a + p2[2])
            # Oszacuj alpha* z polinomem (zero f-r21)
            p2_shifted = p2.copy(); p2_shifted[-1] -= R21_EXP
            roots = np.roots(p2_shifted)
            a_star_fit = float(roots[np.argmin(np.abs(roots - a_star))])

            print(f"  {label}: a*={a_star:.8f}")
            print(f"    df/dalpha (liniowe):    {df_da_lin:+.4f}")
            print(f"    df/dalpha (kwadratowe): {df_da_quad:+.4f}")
            print(f"    alpha* z polinomem:     {a_star_fit:.8f}")
            print(f"    a*-a*_fit:              {a_star-a_star_fit:+.4e}")
            da_dR = 1.0/df_da_quad
            print(f"    dalpha*/dR = 1/(df/da): {da_dR:+.6f}")
            results_deriv[label] = {'df_da': df_da_quad, 'da_dR': da_dR}
            print()
        else:
            print(f"  {label}: za malo punktow!")

    # --- Oblicz dS/dR ---
    if 'a*_1' in results_deriv and 'a*_2' in results_deriv:
        da1_dR = results_deriv['a*_1']['da_dR']
        da2_dR = results_deriv['a*_2']['da_dR']
        dS_dR = da1_dR + da2_dR

        print("=" * 68)
        print("dS/dR I PREDYKCJA r21")
        print("=" * 68)
        print(f"  dalpha*_1/dR = {da1_dR:+.6f}")
        print(f"  dalpha*_2/dR = {da2_dR:+.6f}")
        print(f"  dS/dR        = {dS_dR:+.6f}  (= {dS_dR:.4e})")
        print()

        # R* z warunku S(R*) = 2pi-11/10
        delta_S = S_R21 - S_FORMULA   # S(r21) - 2pi-11/10 = -6.87e-7
        delta_R = delta_S / dS_dR      # dR wymagane do S = 2pi-11/10
        R_star = R21_EXP - delta_R     # R* = r21 - delta_R (bo S(r21) < formula)
        # Uwaga: S(r21) < 2pi-11/10, dS/dR > 0, wiec R* > r21
        R_star2 = R21_EXP + abs(delta_S/dS_dR)

        print(f"  S(r21) - (2pi-11/10)  = {delta_S:+.4e}")
        print(f"  dS/dR                 = {dS_dR:.6f}")
        print(f"  delta_R = dS/dS_dR    = {delta_S/dS_dR:+.4e}")
        print()
        print(f"  R* = r21 + |delta_S|/dS_dR = {R21_EXP} + {abs(delta_S/dS_dR):.4e}")
        print(f"  R* = {R_star2:.8f}")
        print()
        print(f"  Eksperymentalne r21:    {R21_EXP:.8f}  (wartosc uzywana w kodzie)")
        print(f"  CODATA 2018:            {R21_CODATA:.8f}")
        print(f"  TGP predykcja R*:       {R_star2:.8f}")
        print()
        diff_code  = R_star2 - R21_EXP
        diff_codata = R_star2 - R21_CODATA
        rel_code    = diff_code/R21_EXP
        rel_codata  = diff_codata/R21_CODATA
        print(f"  R* - r21 (kod):        {diff_code:+.6e}  ({rel_code*1e6:+.3f} ppm)")
        print(f"  R* - r21 (CODATA):     {diff_codata:+.6e}  ({rel_codata*1e6:+.3f} ppm)")
        print()

        # Analiza bledow
        print("--- Analiza niepewnosci ---")
        print(f"  Niepewnosc dS/dR (10%): +-{0.1*dS_dR:.4e}")
        print(f"  Niepewnosc R*:          +-{0.1*abs(delta_S/dS_dR):.4e}")
        print(f"  Zakres R*:              [{R_star2-0.1*abs(delta_R):.6f}, {R_star2+0.1*abs(delta_R):.6f}]")
        print()

        # CODATA uncertainty
        print(f"  CODATA niepewnosc r21:  +-0.0000046  (46 ppb)")
        print(f"  TGP-CODATA roznica:     {diff_codata:+.6f}  ({rel_codata*1e6:+.2f} ppm)")
        if abs(rel_codata) < 2e-6:
            print(f"  => TGP predykcja bliska CODATA z precyzja < 2 ppm!")
        elif abs(rel_codata) < 5e-6:
            print(f"  => TGP predykcja bliska CODATA z precyzja < 5 ppm")
        else:
            print(f"  => Roznica > 5 ppm od CODATA")
        print()

        # Alternatywa: czy S(CODATA_r21) = 2pi-11/10 lepiej?
        S_at_codata = S_R21 + dS_dR * (R21_CODATA - R21_EXP)
        print(f"  S(r21_CODATA) z interpolacji = {S_at_codata:.10f}")
        print(f"  2pi-11/10                    = {S_FORMULA:.10f}")
        diff_codata_formula = S_at_codata - S_FORMULA
        print(f"  S(r21_CODATA) - (2pi-11/10)  = {diff_codata_formula:+.4e}  ({abs(diff_codata_formula)/S_FORMULA*1e6:.2f} ppm)")

    # TESTY
    print()
    print("=" * 68)
    print("TESTY")
    print("=" * 68)
    if 'a*_1' in results_deriv and 'a*_2' in results_deriv:
        T1 = results_deriv['a*_1']['df_da'] < 0  # f maleje przy a*_1
        print(f"  T1: df/da < 0 przy a*_1:  {'PASS' if T1 else 'FAIL'}  ({results_deriv['a*_1']['df_da']:.2f})")

        T2 = results_deriv['a*_2']['df_da'] > 0  # f rosnie przy a*_2
        print(f"  T2: df/da > 0 przy a*_2:  {'PASS' if T2 else 'FAIL'}  ({results_deriv['a*_2']['df_da']:.2f})")

        T3 = dS_dR > 0
        print(f"  T3: dS/dR > 0:            {'PASS' if T3 else 'FAIL'}  ({dS_dR:.6f})")

        T4 = abs(rel_codata) < 5e-6
        print(f"  T4: |R*-CODATA|/CODATA < 5ppm: {'PASS' if T4 else 'FAIL'}  ({abs(rel_codata)*1e6:.2f} ppm)")

        T5 = R_star2 > R21_EXP  # R* powinna byc nieco powyzej r21 (bo S(r21)<formula)
        print(f"  T5: R* > r21 (kod):       {'PASS' if T5 else 'FAIL'}  (R*={R_star2:.6f} > r21={R21_EXP})")

        n_pass = sum([T1, T2, T3, T4, T5])
        print(f"\nWYNIK: {n_pass}/5")

    print()
    print("=" * 68)
    print("WNIOSEK (EX82)")
    print("=" * 68)
    if 'a*_1' in results_deriv and 'a*_2' in results_deriv:
        print(f"  Precyzyjne dS/dR = {dS_dR:.6f}  (z lokalnego poliniom. w alpha*_i)")
        print(f"  TGP predykcja r21 z S=2pi-11/10: R* = {R_star2:.6f}")
        print(f"  Eksperyment:  r21_CODATA = {R21_CODATA}")
        print(f"  Roznica:      {diff_codata:+.6f} = {rel_codata*1e6:+.3f} ppm")
        print()
        print(f"  Interpretacja:")
        print(f"    Warunek S(R) = 2pi-11/10 wyznacza R z precyzja ~{abs(rel_codata)*1e6:.1f} ppm od r21_CODATA.")
        print(f"    Przy dokladnym r21_CODATA: S(r21_CODATA) = {S_at_codata:.10f}")
        print(f"    vs formula  2pi-11/10    = {S_FORMULA:.10f}")
        print(f"    Odstepstwo               = {diff_codata_formula:+.4e} ({abs(diff_codata_formula)/S_FORMULA*1e6:.2f} ppm)")
    print("=" * 68)
