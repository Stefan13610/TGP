"""
EX74: BIFURKACJA — przy jakim beta dwa zera delta(alpha) zlewaja sie w jedno?

Cel: znalezc krytyczne beta_c (w poblizu beta=1) gdzie delta(alpha) ma
dokladnie jeden dublet zero (d/dalpha = 0 i delta = 0 jednoczesnie).
Dla beta < beta_c: 0 zer, dla beta = beta_c: 1 dublet, dla beta w (beta_c1, beta_c2): 2 zera.
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
        dr = bp*g**2*(1.0-g)
        cr = (ak/g)*gp**2
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
        else:
            break
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
    Ae = _A_inf(g0, ak, bp)
    Am = _A_inf(PHI*g0, ak, bp)
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
    z0  = _find_z0(ak, bp)
    g0s = _find_g0star(ak, bp)
    if z0 and g0s: return g0s - z0
    return np.nan

def scan_beta(bp):
    """Dla ustalonego bp: skanuj alpha w [2.0, 3.8], zwroc liczbe zer i max delta."""
    a_arr = np.linspace(2.0, 3.8, 30)
    ds = [delta(a, bp) for a in a_arr]
    crossings = 0
    max_d = -np.inf
    for i in range(len(a_arr)-1):
        d1, d2 = ds[i], ds[i+1]
        if not np.isnan(d1):
            max_d = max(max_d, d1)
        if not (np.isnan(d1) or np.isnan(d2)) and d1*d2 < 0:
            crossings += 1
    return bp, crossings, max_d, ds, list(a_arr)

if __name__ == '__main__':
    print("=" * 60)
    print("EX74: BIFURKACJA beta — gdzie dwa zera delta(alpha) sie zlewaja?")
    print("=" * 60)
    print(f"  Uzywam {N_CORES} rdzeni")
    print()

    beta_arr = np.array([
        0.80, 0.84, 0.86, 0.88, 0.90, 0.92, 0.94, 0.96, 0.98,
        1.00, 1.02, 1.04, 1.06, 1.08, 1.10, 1.15, 1.20
    ])

    print(f"  Skanuje {len(beta_arr)} wartosci beta...")
    with Pool(N_CORES) as pool:
        results = pool.map(scan_beta, beta_arr)

    print()
    print(f"  {'beta':>6}  {'zera':>5}  {'max_delta':>10}  {'uwaga'}")
    print(f"  {'-'*6}  {'-'*5}  {'-'*10}  {'-'*20}")

    prev_n = None
    for bp, n, md, ds, a_arr in results:
        if n == 2:
            uwaga = "<<< DWIE ZERA"
        elif prev_n == 2 and n == 1:
            uwaga = "<-- bifurkacja gorna?"
        elif prev_n == 1 and n == 2:
            uwaga = "<-- bifurkacja dolna?"
        elif prev_n == 0 and n == 1:
            uwaga = "<-- nowe zero pojawia sie"
        elif prev_n == 1 and n == 0:
            uwaga = "<-- zero znika"
        else:
            uwaga = ""
        print(f"  {bp:6.3f}  {n:5d}  {md:10.5f}  {uwaga}")
        prev_n = n

    # Szukaj bifurkacji: gdzie n zmienia sie z 1 na 2 (od dolu)
    bif_lower = None
    bif_upper = None
    for i in range(len(results)-1):
        n1 = results[i][1]
        n2 = results[i+1][1]
        if n1 < 2 and n2 == 2:
            bif_lower = (results[i][0], results[i+1][0])
        if n1 == 2 and n2 < 2:
            bif_upper = (results[i][0], results[i+1][0])

    print()
    if bif_lower:
        print(f"  Bifurkacja dolna (beta gdzie 2 zera pojawiaja sie): "
              f"[{bif_lower[0]:.3f}, {bif_lower[1]:.3f}]")
    if bif_upper:
        print(f"  Bifurkacja gorna (beta gdzie 2 zera znikaja):       "
              f"[{bif_upper[0]:.3f}, {bif_upper[1]:.3f}]")

    # Szczegol dla beta=1.0
    b1 = [r for r in results if abs(r[0]-1.0) < 1e-9][0]
    print()
    print(f"  Dla beta=1.0: {b1[1]} zer, max_delta = {b1[2]:.6f}")

    # Kandydaci na beta_c
    print()
    print("  Kandydaci na beta_c (dolna bifurkacja):")
    if bif_lower:
        bc_mid = (bif_lower[0] + bif_lower[1]) / 2
        print(f"    ~{bc_mid:.3f}")
        cands = [
            ('1-1/10',    0.9),
            ('9/10',      0.9),
            ('sqrt(4/5)', np.sqrt(4/5)),
            ('phi-1',     PHI-1),
            ('1/phi',     1/PHI),
            ('e/3',       np.e/3),
            ('2/pi-0.17', 2/np.pi - 0.17),
        ]
        for name, val in sorted(cands, key=lambda c: abs(c[1]-bc_mid)):
            print(f"    {name:<20} = {val:.5f}  diff={abs(val-bc_mid):.4f}")

    print()
    print("=" * 60)
    print("WNIOSEK (EX74)")
    print("=" * 60)
    two_betas = [r[0] for r in results if r[1] == 2]
    print(f"  Zakres beta z 2 zerami: {min(two_betas):.2f} do {max(two_betas):.2f}" if two_betas else "  Brak")
    if bif_lower and bif_upper:
        print(f"  Bifurkacja dolna: beta_c1 in [{bif_lower[0]:.3f}, {bif_lower[1]:.3f}]")
        print(f"  Bifurkacja gorna: beta_c2 in [{bif_upper[0]:.3f}, {bif_upper[1]:.3f}]")
        width = bif_upper[0] - bif_lower[1]
        print(f"  Szerokosc obszaru 2-zerowego: ~{width:.3f}")
