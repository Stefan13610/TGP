"""
EX74b: PRECYZYJNA BIFURKACJA — jak waska jest okolica beta=1 z dwoma zerami?
STATUS: LEGACY-TRANSLATIONAL

This script continues the older bifurcation/selection chain around
`delta(alpha)` and `alpha*`. It is useful historically, but it is not part of
the canonical synchronized `nbody` path.

Z ex74: w siatce co 0.02 tylko beta=1.00 ma 2 zera.
Teraz skanujemy [0.97, 1.03] co 0.002 aby znalezc dokladna szerokosc okna.
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
        g, gp = y
        g = max(g, gb + 1e-7)
        fg = 1.0 + 2.0*ak*np.log(g)
        if abs(fg) < 1e-10: return [gp, 0.0]
        dr = bp*g**2*(1.0-g); cr = (ak/g)*gp**2
        if r < 1e-10: return [gp, (dr-cr)/(3.0*fg)]
        return [gp, (dr-cr-fg*2.0*gp/r)/fg]
    def ev(r, y): return y[0] - gb
    ev.terminal = True; ev.direction = -1
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
        p, _ = curve_fit(lambda x, ai, a: ai*(1+a/x), rv, Av,
                         p0=[Av[-1], 0.], maxfev=2000)
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
            try: return brentq(lambda x: _B_coeff(x, ak, bp),
                               gs[i], gs[i+1], xtol=1e-7)
            except: pass
    return None

def _find_g0star(ak, bp, lo=1.05, hi=2.0):
    gs = np.linspace(lo, hi, 25)
    fs = [_f_ratio(g, ak, bp) for g in gs]
    for i in range(len(fs)-1):
        if not (np.isnan(fs[i]) or np.isnan(fs[i+1])):
            if (fs[i]-R21_EXP)*(fs[i+1]-R21_EXP) < 0:
                try: return brentq(lambda x: _f_ratio(x, ak, bp)-R21_EXP,
                                   gs[i], gs[i+1], xtol=1e-7)
                except: pass
    return None

def delta(ak, bp):
    z0 = _find_z0(ak, bp); g0s = _find_g0star(ak, bp)
    if z0 and g0s: return g0s - z0
    return np.nan

def scan_beta_fine(bp):
    """Skanuj alpha w [2.2, 3.4] (okolica znanych zer), zwroc liczbe zer i max delta."""
    a_arr = np.linspace(2.2, 3.4, 24)
    ds = [delta(a, bp) for a in a_arr]
    crossings = 0
    max_d = -np.inf
    zeros_found = []
    for i in range(len(a_arr)-1):
        d1, d2 = ds[i], ds[i+1]
        if not np.isnan(d1): max_d = max(max_d, d1)
        if not (np.isnan(d1) or np.isnan(d2)) and d1*d2 < 0:
            crossings += 1
            try:
                z = brentq(lambda a: delta(a, bp), a_arr[i], a_arr[i+1], xtol=1e-5)
                zeros_found.append(z)
            except: pass
    return bp, crossings, max_d, zeros_found

if __name__ == '__main__':
    print("=" * 62)
    print("EX74b: PRECYZYJNA BIFURKACJA — szerokosc okna beta z 2 zerami")
    print("=" * 62)
    print(f"  Uzywam {N_CORES} rdzeni")
    print()

    beta_arr = np.array([
        0.970, 0.974, 0.978, 0.982, 0.986, 0.988, 0.990, 0.992,
        0.994, 0.996, 0.998, 1.000, 1.002, 1.004, 1.006, 1.008,
        1.010, 1.012, 1.014, 1.016, 1.018, 1.020, 1.024, 1.028
    ])

    print(f"  Skanuje {len(beta_arr)} wartosci beta w [0.97, 1.03]...")
    with Pool(N_CORES) as pool:
        results = pool.map(scan_beta_fine, beta_arr)

    print()
    print(f"  {'beta':>6}  {'zera':>5}  {'max_d':>9}  {'alpha*_1':>10}  {'alpha*_2':>10}")
    print(f"  {'-'*6}  {'-'*5}  {'-'*9}  {'-'*10}  {'-'*10}")

    two_betas   = []
    two_sums    = []
    two_zeros   = []

    for bp, n, md, zf in results:
        a1s = f"{zf[0]:.6f}" if len(zf) >= 1 else "N/A"
        a2s = f"{zf[1]:.6f}" if len(zf) >= 2 else "N/A"
        mark = " <<" if n == 2 else ""
        print(f"  {bp:6.4f}  {n:5d}  {md:9.5f}  {a1s:>10}  {a2s:>10}{mark}")
        if n == 2 and len(zf) >= 2:
            two_betas.append(bp)
            two_sums.append(zf[0]+zf[1])
            two_zeros.append((zf[0], zf[1]))

    print()
    print(f"  Betas z 2 zerami: {two_betas}")
    print()

    if two_betas:
        print(f"  Okno beta z 2 zerami: [{min(two_betas):.4f}, {max(two_betas):.4f}]")
        print(f"  Szerokosc okna: {max(two_betas)-min(two_betas):.4f}")
        print()
        print("  Sumy alpha*_1 + alpha*_2 dla kazdego beta w oknie:")
        S_target = 2*np.pi - 1.1
        for bp, s, (z1, z2) in zip(two_betas, two_sums, two_zeros):
            print(f"    beta={bp:.4f}: suma={s:.7f}  2pi-11/10={S_target:.7f}  diff={s-S_target:+.4e}")
        print()

        # Test czy centrum okna to beta=1
        beta_center = (min(two_betas) + max(two_betas)) / 2
        print(f"  Centrum okna: beta_center = {beta_center:.4f}")
        print(f"  |beta_center - 1| = {abs(beta_center-1):.4f}")

    print()
    print("=" * 62)
    print("WNIOSEK (EX74b)")
    print("=" * 62)
    if two_betas:
        width = max(two_betas) - min(two_betas)
        if width < 0.05:
            print(f"  beta=1 lezy w WYJATKOWO WASKIM oknie {width:.4f}")
            print(f"  szerokosci z dwoma zerami delta(alpha).")
            if abs((min(two_betas)+max(two_betas))/2 - 1.0) < 0.01:
                print(f"  Centrum okna ~ beta=1.0 z dokladnoscia {abs((min(two_betas)+max(two_betas))/2-1):.4f}")
                print(f"  => beta=1 jest WYROZNIONYM PUNKTEM BIFURKACJI!")
    else:
        print("  Brak przypadkow z 2 zerami w przeskanowanym zakresie.")
    print("=" * 62)
