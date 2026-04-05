"""
EX79: MINIMUM f(z0, alpha) — geometria doliny i symetria alpha*

Z ex78: f(z0, alpha) ma minimum ~198 przy alpha~2.57.
Pytanie:
  1. Dokladna wartosc alpha_min i f_min
  2. Czy alpha_min = (alpha*_1 + alpha*_2)/2 = 2.5916?
  3. Czy symetria: alpha*_1 = alpha_min - D, alpha*_2 = alpha_min + D?
  4. Jaka jest wartosc f_min? Czy ma sens fizyczny?
  5. Zwiazek: f_min ~ r21 * cos^2(theta) lub inny?
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar, curve_fit
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

def _f_ratio(g0, ak):
    Ae = _A_inf(g0, ak); Am = _A_inf(PHI*g0, ak)
    return (Am/Ae)**4 if Ae > 1e-6 else np.nan

def _find_z0(ak, lo=1.05, hi=2.0):
    gs = np.linspace(lo, hi, 40)
    Bs = [_B_coeff(g, ak) for g in gs]
    for i in range(len(Bs)-1):
        if Bs[i]*Bs[i+1] < 0:
            try: return brentq(lambda x: _B_coeff(x, ak), gs[i], gs[i+1], xtol=1e-8)
            except: pass
    return None

def fz0(ak):
    z0 = _find_z0(ak)
    if z0 is None: return np.nan
    return _f_ratio(z0, ak)

def worker(ak):
    return ak, fz0(ak)

if __name__ == '__main__':
    print("=" * 64)
    print("EX79: MINIMUM f(z0,alpha) — geometria doliny")
    print("=" * 64)
    print(f"  Uzywam {N_CORES} rdzeni, beta=1")
    print()

    # Gesta siatka w okolicach minimum (alpha ~ 2.55-2.65)
    alpha_fine = np.linspace(2.42, 2.80, 40)
    print(f"  Obliczam {len(alpha_fine)} punktow w [2.42, 2.80]...")
    with Pool(N_CORES) as pool:
        results = pool.map(worker, alpha_fine)

    avals = np.array([r[0] for r in results])
    fvals = np.array([r[1] for r in results])

    # Znajdz minimum
    valid = [(a, f) for a, f in zip(avals, fvals) if not np.isnan(f)]
    if valid:
        a_arr_v = np.array([x[0] for x in valid])
        f_arr_v = np.array([x[1] for x in valid])
        idx_min = np.argmin(f_arr_v)
        a_min_coarse = a_arr_v[idx_min]
        f_min_coarse = f_arr_v[idx_min]
        print(f"\n  Minimum (siatka): f_min = {f_min_coarse:.6f} przy alpha = {a_min_coarse:.5f}")

    print()
    print(f"  {'alpha':>7}  {'f(z0)':>10}  {'f-r21':>9}")
    print(f"  {'-'*7}  {'-'*10}  {'-'*9}")
    for a, f in zip(avals, fvals):
        if not np.isnan(f):
            mark = " <<MIN" if abs(a-a_min_coarse) < 0.02 else ""
            print(f"  {a:7.4f}  {f:10.5f}  {f-R21_EXP:+9.4f}{mark}")

    # Precyzyjne minimum przez brentq na df/dalpha = 0
    # Uzywamy interpolacji
    print()
    print("--- Precyzyjne minimum (brentq na pochodna) ---")
    # Centralne roznice dla pochodnej
    if len(a_arr_v) >= 5:
        # Szukaj gdzie pochodna zmienia znak
        for i in range(1, len(a_arr_v)-1):
            df_lo = f_arr_v[i] - f_arr_v[i-1]
            df_hi = f_arr_v[i+1] - f_arr_v[i]
            if df_lo < 0 and df_hi > 0:
                # Minimum miedzy i-1 a i+1
                try:
                    def neg_fz0(a): return fz0(a)
                    from scipy.optimize import minimize_scalar
                    res = minimize_scalar(neg_fz0,
                                         bounds=(a_arr_v[i-1], a_arr_v[i+1]),
                                         method='bounded',
                                         options={'xatol': 1e-5})
                    a_min = res.x
                    f_min = res.fun
                    print(f"  alpha_min = {a_min:.8f}")
                    print(f"  f_min     = {f_min:.8f}")
                    print(f"  r21       = {R21_EXP:.8f}")
                    print(f"  f_min/r21 = {f_min/R21_EXP:.8f}")
                    print()

                    # Porownaj z centrum alpha*
                    a1s = 2.44143051; a2s = 2.74175411
                    center = (a1s + a2s) / 2
                    print(f"  Centrum [alpha*_1, alpha*_2] = {center:.8f}")
                    print(f"  alpha_min                    = {a_min:.8f}")
                    print(f"  |alpha_min - centrum|        = {abs(a_min-center):.8f}")
                    print(f"  Symetria L = a_min - a*_1 = {a_min-a1s:.6f}")
                    print(f"  Symetria R = a*_2 - a_min = {a2s-a_min:.6f}")
                    print(f"  L/R        = {(a_min-a1s)/(a2s-a_min):.6f}")
                    print()

                    # Kandydaci dla f_min
                    print(f"  Kandydaci dla f_min = {f_min:.6f}:")
                    cands_f = [
                        ('r21*(phi-1)',    R21_EXP*(PHI-1)),
                        ('r21/phi',        R21_EXP/PHI),
                        ('r21-1/phi^2*r21', R21_EXP*(1-1/PHI**2)),
                        ('r21*(2-phi)',    R21_EXP*(2-PHI)),
                        ('r21*e/pi^2',    R21_EXP*np.e/np.pi**2),
                        ('r21-10',         R21_EXP-10),
                        ('r21*(1-1/e)',    R21_EXP*(1-1/np.e)),
                        ('r21*3/pi',      R21_EXP*3/np.pi),
                        ('200',            200.0),
                        ('r21*ln(r21)/pi^2', R21_EXP*np.log(R21_EXP)/np.pi**2),
                        ('sqrt(r21)*10',  np.sqrt(R21_EXP)*10),
                    ]
                    for name, val in sorted(cands_f, key=lambda c: abs(c[1]-f_min)):
                        print(f"    {name:<30} = {val:.5f}  diff={abs(val-f_min):.4f}")

                    # Kandydaci dla alpha_min
                    print(f"\n  Kandydaci dla alpha_min = {a_min:.8f}:")
                    cands_a = [
                        ('(a1+a2)/2',        center),
                        ('pi-phi/2',         np.pi - PHI/2),
                        ('e-1/e+1',          np.e - 1/np.e + 1),
                        ('(2pi-11/10)/2',    (2*np.pi-1.1)/2),
                        ('phi+1',            PHI+1),
                        ('sqrt(phi^3)',       np.sqrt(PHI**3)),
                        ('2+phi/phi^2',      2 + PHI/PHI**2),
                        ('pi-phi/pi',        np.pi - PHI/np.pi),
                        ('e/phi+1',          np.e/PHI + 1),
                        ('5/phi-phi/2',      5/PHI - PHI/2),
                    ]
                    for name, val in sorted(cands_a, key=lambda c: abs(c[1]-a_min)):
                        print(f"    {name:<30} = {val:.8f}  diff={abs(val-a_min):.2e}")

                except Exception as e:
                    print(f"  Blad: {e}")
                break

    # Testy
    print()
    print("=" * 64)
    print("TESTY")
    print("=" * 64)

    a1s = 2.44143051; a2s = 2.74175411
    center = (a1s + a2s) / 2

    T1 = not np.isnan(f_min_coarse) and f_min_coarse < R21_EXP
    print(f"  T1: f_min < r21:            {'PASS' if T1 else 'FAIL'}  (f_min={f_min_coarse:.4f})")

    T2 = abs(a_min_coarse - center) < 0.05
    print(f"  T2: alpha_min ~ centrum:    {'PASS' if T2 else 'FAIL'}  (|a_min-cen|={abs(a_min_coarse-center):.5f})")

    # Symetria: L == R?
    L = a_min_coarse - a1s; R = a2s - a_min_coarse
    T3 = abs(L - R) / ((L+R)/2) < 0.1   # < 10% asymetria
    print(f"  T3: Symetria L~R (<10%):   {'PASS' if T3 else 'FAIL'}  (L={L:.4f}, R={R:.4f}, L/R={L/R:.4f})")

    T4 = not np.isnan(f_min_coarse) and abs(f_min_coarse/R21_EXP - 1) > 0.01
    print(f"  T4: f_min != r21 (>1%):    {'PASS' if T4 else 'FAIL'}  (ratio={f_min_coarse/R21_EXP:.5f})")

    T5 = not np.isnan(f_min_coarse) and f_min_coarse > 100
    print(f"  T5: f_min > 100:           {'PASS' if T5 else 'FAIL'}  (f_min={f_min_coarse:.4f})")

    n_pass = sum([T1, T2, T3, T4, T5])
    print(f"\nWYNIK: {n_pass}/5")

    print()
    print("=" * 64)
    print("WNIOSEK (EX79)")
    print("=" * 64)
    print(f"  f_min        ~ {f_min_coarse:.4f}")
    print(f"  alpha_min    ~ {a_min_coarse:.5f}")
    print(f"  centrum      = {center:.5f}")
    print(f"  |a_min-cen|  = {abs(a_min_coarse-center):.5f}")
    if T2:
        print(f"  >> alpha_min ZBLIZONE do centrum [alpha*_1, alpha*_2]!")
        print(f"     Dolina f(z0) jest symetryczna wokol centrum.")
    if T3:
        print(f"  >> Symetria L~R potwierdzona (asymetria < 10%)")
    print("=" * 64)
