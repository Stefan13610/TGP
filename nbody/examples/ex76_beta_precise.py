"""
EX76: CZY beta_c = 1 SCISLE? Weryfikacja przy zwiekszonym rtol/atol

Z ex75: beta_c_upper = 1.000050 (gdzie min delta[2.5,2.72] = 0).
Pytanie: czy to prawdziwy wynik czy artefakt numeryczny?

Plan:
  1. Oblicz delta(alpha, beta) dla alpha w [alpha*_1, alpha*_2] przy beta=1.000 i beta=1.0001
     z trojkrotnie wyzszej dokladnoscia (rtol=1e-11, atol=1e-13).
  2. Znajdz min delta w przedziale i porownaj z wynikiem ex75.
  3. Jezeli min_delta(beta=1) jest WYRAZNIE UJEMNE przy wysokiej precyzji:
     => beta_c > 1, bifurkacja nie jest dokladnie w beta=1.
  4. Jezeli min_delta(beta=1) dazy do 0 przy zwiekszaniu precyzji:
     => beta_c = 1 moze byc scisle.
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, curve_fit, minimize_scalar
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

def _integrate(g0, ak, bp, rtol=1e-10, atol=1e-12, maxb=8):
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
                        dense_output=True, rtol=rtol, atol=atol, max_step=0.04)
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

def _A_inf(g0, ak, bp, rtol=1e-10, atol=1e-12):
    r, g = _integrate(g0, ak, bp, rtol, atol)
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

def _B_coeff(g0, ak, bp, rtol=1e-10, atol=1e-12):
    r, g = _integrate(g0, ak, bp, rtol, atol)
    _, B, _ = _fit_win(r, g, WIN_LIST[2], WIN_LIST[2]+WIN_WIDTH)
    return B

def _f_ratio(g0, ak, bp, rtol=1e-10, atol=1e-12):
    Ae = _A_inf(g0, ak, bp, rtol, atol)
    Am = _A_inf(PHI*g0, ak, bp, rtol, atol)
    return (Am/Ae)**4 if Ae > 1e-6 else np.nan

def _find_z0(ak, bp, rtol=1e-10, atol=1e-12, lo=1.05, hi=2.0):
    gs = np.linspace(lo, hi, 40)
    Bs = [_B_coeff(g, ak, bp, rtol, atol) for g in gs]
    for i in range(len(Bs)-1):
        if Bs[i]*Bs[i+1] < 0:
            try:
                return brentq(lambda x: _B_coeff(x, ak, bp, rtol, atol),
                              gs[i], gs[i+1], xtol=1e-8)
            except: pass
    return None

def _find_g0star(ak, bp, rtol=1e-10, atol=1e-12, lo=1.05, hi=2.0):
    gs = np.linspace(lo, hi, 25)
    fs = [_f_ratio(g, ak, bp, rtol, atol) for g in gs]
    for i in range(len(fs)-1):
        if not (np.isnan(fs[i]) or np.isnan(fs[i+1])):
            if (fs[i]-R21_EXP)*(fs[i+1]-R21_EXP) < 0:
                try:
                    return brentq(lambda x: _f_ratio(x, ak, bp, rtol, atol)-R21_EXP,
                                  gs[i], gs[i+1], xtol=1e-8)
                except: pass
    return None

def delta_precise(ak, bp, rtol=1e-10, atol=1e-12):
    z0  = _find_z0(ak, bp, rtol, atol)
    g0s = _find_g0star(ak, bp, rtol, atol)
    if z0 and g0s: return g0s - z0
    return np.nan

# ============================================================
# Funkcja robocza: oblicz delta w gestej siatce alpha dla danego beta
# ============================================================
def compute_delta_profile(args):
    bp, rtol, atol = args
    a_arr = np.linspace(2.40, 2.82, 22)
    ds = []
    for a in a_arr:
        d = delta_precise(a, bp, rtol, atol)
        ds.append(d)
    return bp, list(a_arr), ds

if __name__ == '__main__':
    print("=" * 66)
    print("EX76: CZY beta_c = 1 SCISLE? Wysokoprecyzyjna weryfikacja")
    print("=" * 66)
    print(f"  Uzywam {N_CORES} rdzeni")
    print()

    # Sekcja 1: profil delta(alpha) dla beta=1.000 vs beta=1.0001
    # przy dwoch poziomach precyzji: standardowej i wysokiej
    print("--- Sekcja 1: Profil delta(alpha) dla beta=1.000 ---")
    print()

    betas_test = [0.9995, 0.9998, 1.0000, 1.0002, 1.0005]
    cases = [(bp, 1e-10, 1e-12) for bp in betas_test]

    with Pool(N_CORES) as pool:
        profiles = pool.map(compute_delta_profile, cases)

    for bp, a_arr, ds in profiles:
        valid = [(a, d) for a, d in zip(a_arr, ds) if not np.isnan(d)]
        if valid:
            min_d = min(d for _, d in valid)
            max_d = max(d for _, d in valid)
            # Crossings
            crossings = 0
            for i in range(len(a_arr)-1):
                d1, d2 = ds[i], ds[i+1]
                if not (np.isnan(d1) or np.isnan(d2)) and d1*d2 < 0:
                    crossings += 1
            print(f"  beta={bp:.4f}: min={min_d:+.6f}  max={max_d:+.6f}  crossings={crossings}")
        else:
            print(f"  beta={bp:.4f}: brak danych")

    # Sekcja 2: min_delta w [2.5, 2.72] vs precyzja dla beta=1.000
    print()
    print("--- Sekcja 2: min delta[2.5,2.72] dla beta=1.000 vs precyzja ---")
    print()

    def min_delta_midrange(bp, rtol, atol):
        a_arr = np.linspace(2.52, 2.70, 10)
        ds = [delta_precise(a, bp, rtol, atol) for a in a_arr]
        valid = [d for d in ds if not np.isnan(d)]
        return min(valid) if valid else np.nan

    configs = [
        (1e-9,  1e-11, "standard (ex71)"),
        (1e-10, 1e-12, "wysoka"),
        (1e-11, 1e-13, "bardzo wysoka"),
    ]

    print(f"  {'rtol':>8}  {'atol':>8}  {'min_delta[2.5,2.72]':>20}  opis")
    print(f"  {'-'*8}  {'-'*8}  {'-'*20}  {'-'*20}")
    for rtol, atol, desc in configs:
        md = min_delta_midrange(1.000, rtol, atol)
        print(f"  {rtol:8.0e}  {atol:8.0e}  {md:+20.8f}  {desc}")

    # Sekcja 3: Precyzyjne brentq na beta_c przy wysokiej precyzji
    print()
    print("--- Sekcja 3: Brentq beta_c przy rtol=1e-10 ---")
    print()

    def min_delta_mid_hp(bp):
        return min_delta_midrange(bp, 1e-10, 1e-12)

    # Najpierw sprawdz znak na kilku punktach
    test_betas2 = [0.998, 0.999, 1.000, 1.001, 1.002]
    print("  Sprawdzam znaki:")
    for b in test_betas2:
        md = min_delta_mid_hp(b)
        print(f"    beta={b:.4f}: min_delta={md:+.6f}")

    print()
    # Znajdz dolna i gorna bifurkacje
    try:
        # Dolna: delta idzie z ujemnej na dodatnia
        # Szukaj w [0.999, 1.001]
        md_lo = min_delta_mid_hp(0.999)
        md_hi = min_delta_mid_hp(1.001)
        if md_lo * md_hi < 0:
            bc = brentq(min_delta_mid_hp, 0.999, 1.001, xtol=1e-5)
            print(f"  beta_c (brentq, rtol=1e-10) = {bc:.7f}")
            print(f"  |beta_c - 1| = {abs(bc-1):.7f}")
            # Testuj rozne hipotezy
            print()
            print(f"  Kandydaci dla beta_c:")
            cands = [
                ('1',              1.0),
                ('1+1e-4',         1.0001),
                ('1+1/(2*pi^2)',   1.0 + 1/(2*np.pi**2)),
                ('1+1/10000',      1.0001),
                ('phi-phi+1',      1.0),
            ]
            for name, val in sorted(cands, key=lambda c: abs(c[1]-bc)):
                print(f"    {name:<20} = {val:.7f}  diff={abs(val-bc):.2e}")
        else:
            print(f"  Brak przejscia znaku w [0.999, 1.001]: md_lo={md_lo:+.5f}, md_hi={md_hi:+.5f}")
    except Exception as e:
        print(f"  Blad: {e}")

    # Sekcja 4: Sprawdz czy delta(alpha*_mid, beta=1) jest dokladnie 0
    print()
    print("--- Sekcja 4: delta w punkcie alfa_mid = (a*1+a*2)/2 ---")
    a_mid = (2.44143051 + 2.74175411) / 2  # = 2.59159
    print(f"  alpha_mid = (a*1 + a*2)/2 = {a_mid:.8f}")
    print()
    for rtol, atol, desc in configs:
        d = delta_precise(a_mid, 1.000, rtol, atol)
        print(f"  rtol={rtol:.0e}: delta(alpha_mid, beta=1) = {d:+.8f}  [{desc}]")

    print()
    print("=" * 66)
    print("WNIOSEK (EX76)")
    print("=" * 66)
    print("  Jesli min_delta(beta=1) dazy do 0 przy zwiekszaniu precyzji:")
    print("    => beta_c = 1 SCISLE (dokl. numeryczna)")
    print("  Jesli min_delta(beta=1) < 0 nawet przy rtol=1e-11:")
    print("    => beta_c > 1, a beta=1 lezy tuz ponizej bifurkacji")
    print("=" * 66)
