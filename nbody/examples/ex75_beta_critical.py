"""
EX75: PRECYZYJNA WARTOSC beta_c — czy beta=1 jest dokladnym punktem bifurkacji?

Z ex74b: tylko beta=1.000 (rozdzielczosc 0.002) daje 2 zera delta(alpha).
Teraz:
  1. Siatka Deltabeta=0.0002 wokol beta=1 => czy okno szerokosci 0?
  2. Szukaj dokladnego beta_c (dolna granica) przez brentq na max_delta(beta)=0
  3. Sprawdz czy beta_c = 1 dokladnie czy np. 1 +/- epsilon

Kluczowy test: max_delta(beta) jest funkcja beta. Okolo beta_c ta funkcja
zmienia znak (z ujemnej na dodatnia). brentq na tej funkcji daje beta_c.
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

def max_delta_in_range(bp, a_lo=2.3, a_hi=2.9, n=16):
    """Maksimum delta(alpha, bp) w przedziale alpha okolo zer."""
    a_arr = np.linspace(a_lo, a_hi, n)
    ds = [delta(a, bp) for a in a_arr]
    valid = [d for d in ds if not np.isnan(d)]
    return max(valid) if valid else -np.inf

def scan_worker(bp):
    n_cross = 0
    a_arr = np.linspace(2.2, 3.2, 20)
    ds = [delta(a, bp) for a in a_arr]
    md = max((d for d in ds if not np.isnan(d)), default=-np.inf)
    for i in range(len(a_arr)-1):
        d1, d2 = ds[i], ds[i+1]
        if not (np.isnan(d1) or np.isnan(d2)) and d1*d2 < 0:
            n_cross += 1
    return bp, n_cross, md

if __name__ == '__main__':
    print("=" * 64)
    print("EX75: PRECYZYJNE beta_c — czy beta=1 jest dokladnym punktem?")
    print("=" * 64)
    print(f"  Uzywam {N_CORES} rdzeni")
    print()

    # Siatka Deltabeta = 0.0002 w [0.994, 1.008]
    beta_arr = np.arange(0.9940, 1.0085, 0.0002)
    print(f"  Sekcja 1: Siatka Deltabeta=0.0002, {len(beta_arr)} punktow...")

    with Pool(N_CORES) as pool:
        results = pool.map(scan_worker, beta_arr)

    print()
    print(f"  {'beta':>7}  {'zera':>5}  {'max_d':>9}  {'uwaga'}")
    print(f"  {'-'*7}  {'-'*5}  {'-'*9}  {'-'*15}")

    betas_2 = []
    for bp, n, md in results:
        mark = " <<< DWIE!" if n == 2 else ""
        print(f"  {bp:7.4f}  {n:5d}  {md:9.5f}{mark}")
        if n == 2:
            betas_2.append(bp)

    print()
    print(f"  Betas z 2 zerami (Deltabeta=0.0002): {betas_2}")

    # Sekcja 2: max_delta(beta) jako funkcja beta — szukaj zer (bifurkacji)
    print()
    print("-" * 64)
    print("  Sekcja 2: Szukaj beta_c gdzie max_delta(beta) = 0")
    print()

    # max_delta w oknie alpha=[2.5, 2.8] (miedzy dwoma zerami dla beta=1)
    def max_d_mid(bp):
        return max_delta_in_range(bp, a_lo=2.50, a_hi=2.72, n=10)

    # Sprawdz znak na kilku punktach
    test_betas = [0.995, 0.998, 1.000, 1.002, 1.005]
    print("  max_delta w oknie [2.5, 2.72]:")
    for b in test_betas:
        md = max_d_mid(b)
        print(f"    beta={b:.4f}: max_delta={md:.6f}")

    # Sprobuj znalezc beta_c dolna i gorna
    print()
    print("  Szukam beta_c (dolna) przez brentq...")
    try:
        # Dolna: gdzie max_delta przechodzi przez 0 z dolu
        # Sprawdz przedzialy
        md_low  = max_d_mid(0.996)
        md_high = max_d_mid(1.000)
        print(f"    max_d(0.996)={md_low:.5f}, max_d(1.000)={md_high:.5f}")
        if md_low * md_high < 0:
            bc_lo = brentq(max_d_mid, 0.996, 1.000, xtol=1e-4)
            print(f"    beta_c_lower = {bc_lo:.6f}")
            print(f"    |beta_c_lower - 1| = {abs(bc_lo-1):.6f}")
        else:
            print("    Brak zmiany znaku w [0.996, 1.000] — niemozliwy brentq")
    except Exception as e:
        print(f"    Blad: {e}")

    print()
    print("  Szukam beta_c (gorna) przez brentq...")
    try:
        md_lo = max_d_mid(1.000)
        md_hi = max_d_mid(1.003)
        print(f"    max_d(1.000)={md_lo:.5f}, max_d(1.003)={md_hi:.5f}")
        if md_lo * md_hi < 0:
            bc_up = brentq(max_d_mid, 1.000, 1.003, xtol=1e-4)
            print(f"    beta_c_upper = {bc_up:.6f}")
            print(f"    |beta_c_upper - 1| = {abs(bc_up-1):.6f}")
        else:
            print("    Brak zmiany znaku w [1.000, 1.003]")
    except Exception as e:
        print(f"    Blad: {e}")

    # Sekcja 3: max_delta w oknie ZAWIERAJACYM oba zera (dla beta=1)
    print()
    print("-" * 64)
    print("  Sekcja 3: max_delta w szerszym oknie [2.3, 3.0] vs beta")
    def max_d_wide(bp):
        return max_delta_in_range(bp, a_lo=2.3, a_hi=3.0, n=14)

    beta_fine2 = np.arange(0.997, 1.004, 0.0005)
    print(f"  {'beta':>7}  {'max_d_wide':>12}  {'uwaga'}")
    print(f"  {'-'*7}  {'-'*12}  {'-'*10}")
    for b in beta_fine2:
        md = max_d_wide(b)
        mark = " +" if md > 0 else " -"
        print(f"  {b:7.4f}  {md:12.6f}{mark}")

    print()
    print("=" * 64)
    print("WNIOSEK (EX75)")
    print("=" * 64)
    if betas_2:
        print(f"  Betas z 2 zerami (Deltabeta=0.0002): {betas_2}")
        if len(betas_2) == 1 and abs(betas_2[0] - 1.0) < 1e-9:
            print(f"  => beta=1 jest JEDYNYM punktem z 2 zerami")
            print(f"     (rozdzielczosc 0.0002, okno alpha=[2.2, 3.2])")
            print(f"  => beta_c = 1 dokladnie (lub z dokladnoscia < 0.0002)")
        else:
            print(f"  Okno beta_c: [{min(betas_2):.4f}, {max(betas_2):.4f}]")
            print(f"  Szerokosc: {max(betas_2)-min(betas_2):.4f}")
    else:
        print("  Brak przypadkow z 2 zerami w Deltabeta=0.0002.")
        print("  beta=1 nie jest punktem 2-zerowym przy tej rozdzielczosci?")
    print("=" * 64)
