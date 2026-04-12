"""
EX77: CZY OKNO BIFURKACJI JEST SYMETRYCZNE WOKOL beta=1?
STATUS: LEGACY-TRANSLATIONAL

This file remains in the older bifurcation/symmetry exploration around
`delta(alpha)` and legacy selection variables. Keep it as historical
exploration, not as a canonical current `nbody` path.

Z ex76:
  - beta_c_upper ~ 1.0003 (min delta[2.5,2.72] = 0)
  - beta_c_lower ~ 0.9997 (max delta[2.4,2.82] = 0)

Pytanie: czy beta_c_lower + beta_c_upper = 2, tzn. srednia = 1 dokladnie?
Jezeli tak => beta=1 jest centrum okna bifurkacji (nawet jesli nie krawedzia).
Jezeli nie => beta=1 jest przypadkowe.

Plan:
  1. Precyzyjna dolna bifurkacja: gdzie max(delta w [2.40,2.82]) = 0
  2. Precyzyjna gorna bifurkacja: gdzie min(delta w [2.50,2.72]) = 0
  3. Oblicz centrum i asymetrie
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
R21_EXP   = 206.768
G_OFF     = 0.005
R_MAX     = 100.0
WIN_LIST  = [16, 22, 28, 36, 46, 58, 72]
WIN_WIDTH = 14.0
N_CORES   = min(16, cpu_count())

def _integrate(g0, ak, bp, maxb=8):
    gb = np.exp(-1.0/(2.0*ak)) + G_OFF
    def rhs(r, y):
        g, gp = y; g = max(g, gb+1e-7)
        fg = 1.0 + 2.0*ak*np.log(g)
        if abs(fg) < 1e-10: return [gp, 0.0]
        dr = bp*g**2*(1.0-g); cr = (ak/g)*gp**2
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

def _A_inf(g0, ak, bp):
    r, g = _integrate(g0, ak, bp)
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

def _B_coeff(g0, ak, bp):
    r, g = _integrate(g0, ak, bp)
    _, B, _ = _fit_win(r, g, WIN_LIST[2], WIN_LIST[2]+WIN_WIDTH)
    return B

def _f_ratio(g0, ak, bp):
    Ae = _A_inf(g0, ak, bp); Am = _A_inf(PHI*g0, ak, bp)
    return (Am/Ae)**4 if Ae > 1e-6 else np.nan

def _find_z0(ak, bp, lo=1.05, hi=2.0):
    gs = np.linspace(lo, hi, 40)
    Bs = [_B_coeff(g, ak, bp) for g in gs]
    for i in range(len(Bs)-1):
        if Bs[i]*Bs[i+1] < 0:
            try: return brentq(lambda x: _B_coeff(x, ak, bp), gs[i], gs[i+1], xtol=1e-8)
            except: pass
    return None

def _find_g0star(ak, bp, lo=1.05, hi=2.0):
    gs = np.linspace(lo, hi, 25)
    fs = [_f_ratio(g, ak, bp) for g in gs]
    for i in range(len(fs)-1):
        if not (np.isnan(fs[i]) or np.isnan(fs[i+1])):
            if (fs[i]-R21_EXP)*(fs[i+1]-R21_EXP) < 0:
                try: return brentq(lambda x: _f_ratio(x, ak, bp)-R21_EXP, gs[i], gs[i+1], xtol=1e-8)
                except: pass
    return None

def delta(ak, bp):
    z0 = _find_z0(ak, bp); g0s = _find_g0star(ak, bp)
    if z0 and g0s: return g0s - z0
    return np.nan

def max_delta_outer(bp):
    """Max delta w szerszym oknie [2.40, 2.82] — dla dolnej bifurkacji."""
    a_arr = np.linspace(2.40, 2.82, 14)
    ds = [delta(a, bp) for a in a_arr]
    valid = [d for d in ds if not np.isnan(d)]
    return max(valid) if valid else -np.inf

def min_delta_inner(bp):
    """Min delta w waskim oknie [2.52, 2.70] — dla gornej bifurkacji."""
    a_arr = np.linspace(2.52, 2.70, 10)
    ds = [delta(a, bp) for a in a_arr]
    valid = [d for d in ds if not np.isnan(d)]
    return min(valid) if valid else np.inf

def worker_outer(bp):
    return bp, max_delta_outer(bp)

def worker_inner(bp):
    return bp, min_delta_inner(bp)

if __name__ == '__main__':
    print("=" * 64)
    print("EX77: SYMETRIA OKNA BIFURKACJI — czy centrum = beta=1?")
    print("=" * 64)
    print(f"  Uzywam {N_CORES} rdzeni, rtol=1e-10")
    print()

    # --- Dolna bifurkacja: max_delta_outer = 0 ---
    print("--- Dolna bifurkacja (max_delta[2.40,2.82] = 0) ---")
    beta_lo = np.linspace(0.9960, 1.0005, 20)
    with Pool(N_CORES) as pool:
        res_lo = pool.map(worker_outer, beta_lo)

    print(f"  {'beta':>7}  {'max_d_outer':>12}")
    for bp, md in res_lo:
        mark = " +" if md > 0 else " -"
        print(f"  {bp:7.4f}  {md:12.7f}{mark}")

    # Znajdz beta_c_lower przez brentq
    vals_lo = dict(res_lo)
    sign_changes = [(beta_lo[i], beta_lo[i+1])
                    for i in range(len(beta_lo)-1)
                    if vals_lo[beta_lo[i]] * vals_lo[beta_lo[i+1]] < 0]
    bc_lower = None
    if sign_changes:
        a, b = sign_changes[0]
        bc_lower = brentq(max_delta_outer, a, b, xtol=1e-5)
        print(f"\n  beta_c_lower = {bc_lower:.7f}")
        print(f"  |beta_c_lower - 1| = {abs(bc_lower-1):.7f}")

    print()

    # --- Gorna bifurkacja: min_delta_inner = 0 ---
    print("--- Gorna bifurkacja (min_delta[2.52,2.70] = 0) ---")
    beta_hi = np.linspace(0.9995, 1.0015, 20)
    with Pool(N_CORES) as pool:
        res_hi = pool.map(worker_inner, beta_hi)

    print(f"  {'beta':>7}  {'min_d_inner':>12}")
    for bp, md in res_hi:
        mark = " +" if md > 0 else " -"
        print(f"  {bp:7.4f}  {md:12.7f}{mark}")

    vals_hi = dict(res_hi)
    sign_changes_hi = [(beta_hi[i], beta_hi[i+1])
                       for i in range(len(beta_hi)-1)
                       if vals_hi[beta_hi[i]] * vals_hi[beta_hi[i+1]] < 0]
    bc_upper = None
    if sign_changes_hi:
        a, b = sign_changes_hi[0]
        bc_upper = brentq(min_delta_inner, a, b, xtol=1e-5)
        print(f"\n  beta_c_upper = {bc_upper:.7f}")
        print(f"  |beta_c_upper - 1| = {abs(bc_upper-1):.7f}")

    # --- Analiza symetrii ---
    print()
    print("=" * 64)
    print("ANALIZA SYMETRII")
    print("=" * 64)
    if bc_lower and bc_upper:
        center = (bc_lower + bc_upper) / 2
        width  = bc_upper - bc_lower
        asym   = (1.0 - bc_lower) / width  # pozycja beta=1 w oknie (0=dolna, 1=gorna)
        print(f"  beta_c_lower = {bc_lower:.7f}")
        print(f"  beta_c_upper = {bc_upper:.7f}")
        print(f"  centrum      = {center:.7f}")
        print(f"  szerokosc    = {width:.7f}")
        print(f"  |centrum - 1| = {abs(center-1):.7f}")
        print(f"  pozycja beta=1 w oknie: {asym:.4f} (0=dolna, 1=gorna)")
        print()
        if abs(center - 1.0) < 0.0005:
            print("  >> CENTRUM OKNA ≈ beta=1 (dokladnosc < 5e-4)!")
            print("  >> beta=1 jest centrum bifurkacji, nie krawedzia.")
        else:
            print(f"  >> Centrum = {center:.6f} != 1.0 (diff = {abs(center-1):.6f})")

        # Kandydaci dla bc_lower i bc_upper
        print()
        print(f"  Kandydaci dla beta_c_lower = {bc_lower:.7f}:")
        cands_lo = [
            ('1 - 3e-4',    1 - 3e-4),
            ('1 - 1/3000',  1 - 1/3000),
            ('1 - width/2', 1 - width/2),
            ('1 - 1/pi^3',  1 - 1/np.pi**3),
            ('phi^(-10)',   PHI**(-10)),
        ]
        for name, val in sorted(cands_lo, key=lambda c: abs(c[1]-bc_lower)):
            print(f"    {name:<20} = {val:.7f}  diff={abs(val-bc_lower):.2e}")

        print(f"\n  Kandydaci dla beta_c_upper = {bc_upper:.7f}:")
        cands_hi = [
            ('1 + 3e-4',    1 + 3e-4),
            ('1 + 1/3000',  1 + 1/3000),
            ('1 + width/2', 1 + width/2),
            ('1 + 1/pi^3',  1 + 1/np.pi**3),
        ]
        for name, val in sorted(cands_hi, key=lambda c: abs(c[1]-bc_upper)):
            print(f"    {name:<20} = {val:.7f}  diff={abs(val-bc_upper):.2e}")

    print()
    print("=" * 64)
    print("WNIOSEK (EX77)")
    print("=" * 64)
    if bc_lower and bc_upper:
        if abs((bc_lower+bc_upper)/2 - 1.0) < 0.001:
            print(f"  Centrum okna bifurkacji = {(bc_lower+bc_upper)/2:.6f} ≈ 1.0")
            print(f"  => beta=1 jest CENTRUM okna (symetria!)")
            print(f"  => Wyroznienie beta=1 zachowane, ale jest centrum, nie krawedzia.")
        else:
            print(f"  Centrum = {(bc_lower+bc_upper)/2:.6f}, beta=1 nie jest centrum.")
    print("=" * 64)
