"""
EX78: CZY alpha_TGP=2 JEST CENTRUM OKNA W PRZESTRZENI alpha?
STATUS: LEGACY-TRANSLATIONAL

This script continues the older `alpha*` / window-center exploration in legacy
selection variables. It should be read as transitional historical material, not
as a canonical synchronized example.

Analogia do ex77: beta=1 jest centrum okna bifurkacji w przestrzeni beta.
Pytanie: czy alpha_TGP=2 jest centrum jakiegos okna w przestrzeni alpha?

Kontekst: ODE TGP ma f(g)=1+2*alpha*ln(g). Przy alpha=2 mamy dwa zera delta.
Ale "alpha" w tym kontekscie jest PARAMETREM KINETYCZNYM, a nie "scanning variable".

Nowe pytanie: jezeli zmieniamy PARAMETR KINETYCZNY (alpha_kin_fixed) przy FIXED beta=1,
to przy jakich wartosciach alpha_kin_fixed istnieja dwa zera delta(alpha_scan)?

Ale to jest tautologia — w ex71/ex72 alpha_kin_fixed = alpha_scan.

Inny pomysl: Czy alpha_TGP=2 jest wyjatkowy w przestrzeni samej ODE?
Sprawdzmy: dla ODE z f(g)=1+2*alpha_K*ln(g) ale zmiennym alpha_K (nie alpha_scan),
czy istnieje analogiczne okno alpha_K z dwoma zerami?

UWAGA: W ex71-ex77 "alpha" byl jednoczesnie parametrem kinetycznym I zmienna skanowania.
Separacja mozliwa przez: ustalamy alpha_K (kinetyczny, definiuje f(g)) ale skanujemy
cos innego... Ale nie ma nic innego do skanowania w prostej TGP.

Zamiast tego: zbadamy CZY alpha_TGP=2 jest centrum okna WARTOSCI f(z0,alpha) = r21.
Tzn. przy jakich alpha (kinetycznych) istnieje paraboloidalna struktura f(z0,alpha)
i czy minimum/maksimum f(z0,alpha) jest w okolicach alpha=2?

Plan alternatywny (realizowalny):
  Dla stalego beta=1, skanuj alpha_kin w [1.5, 3.5]:
  1. Oblicz f(z0(alpha), alpha) dla kazdego alpha
  2. Sprawdz gdzie f=r21 (znamy: alpha*_1=2.44, alpha*_2=2.74)
  3. Sprawdz czy alpha_TGP=2 jest centrum [alpha*_1, alpha*_2] lub innego okna
  4. Oblicz f(z0, alpha=2) i sprawdz jego distance od r21
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

def _integrate(g0, ak, maxb=8):
    bp = 1.0  # beta=1 fixed
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

def compute_fz0(ak):
    """Oblicz f(z0(alpha), alpha) dla danego alpha."""
    z0 = _find_z0(ak)
    if z0 is None: return ak, None, np.nan
    fval = _f_ratio(z0, ak)
    return ak, z0, fval

if __name__ == '__main__':
    print("=" * 66)
    print("EX78: f(z0, alpha) vs alpha — gdzie centrum i symetria?")
    print("=" * 66)
    print(f"  Uzywam {N_CORES} rdzeni, beta=1 fixed")
    print()

    # Skan alpha w [1.6, 3.2]
    alpha_arr = np.array([
        1.6, 1.8, 2.0, 2.1, 2.2, 2.3, 2.35, 2.40, 2.42, 2.44,
        2.441, 2.443, 2.445, 2.50, 2.55, 2.60, 2.65, 2.70,
        2.741, 2.742, 2.743, 2.75, 2.80, 2.90, 3.0, 3.2
    ])

    print(f"  Obliczam {len(alpha_arr)} punktow alpha...")
    with Pool(N_CORES) as pool:
        results = pool.map(compute_fz0, alpha_arr)

    print()
    print(f"  {'alpha':>7}  {'z0':>10}  {'f(z0)':>10}  {'f-r21':>9}  uwaga")
    print(f"  {'-'*7}  {'-'*10}  {'-'*10}  {'-'*9}  {'-'*15}")

    fvals = []
    avals = []
    for ak, z0, fv in results:
        z0s  = f"{z0:.6f}" if z0 else "N/A"
        fvs  = f"{fv:.4f}" if not np.isnan(fv) else "N/A"
        dfs  = f"{fv-R21_EXP:+.4f}" if not np.isnan(fv) else "N/A"
        mark = ""
        if not np.isnan(fv):
            if abs(fv - R21_EXP) < 1.0: mark = " << r21!"
            elif abs(fv - R21_EXP) < 10: mark = " < bliski"
            fvals.append(fv); avals.append(ak)
        print(f"  {ak:7.4f}  {z0s:>10}  {fvs:>10}  {dfs:>9}  {mark}")

    # Analiza
    print()
    print("-" * 66)
    print("--- Analiza f(z0, alpha) ---")
    print()

    avals = np.array(avals); fvals = np.array(fvals)

    # Gdzie f = r21?
    crosses = []
    for i in range(len(avals)-1):
        if not (np.isnan(fvals[i]) or np.isnan(fvals[i+1])):
            if (fvals[i]-R21_EXP)*(fvals[i+1]-R21_EXP) < 0:
                crosses.append((avals[i], avals[i+1], fvals[i], fvals[i+1]))
    print(f"  Przejscia f(z0)=r21:")
    for a1, a2, f1, f2 in crosses:
        print(f"    alpha in [{a1:.4f}, {a2:.4f}]  (f: {f1:.3f} -> {f2:.3f})")

    # f przy alpha=2 (TGP)
    r2 = [r for r in results if abs(r[0]-2.0) < 1e-9]
    if r2:
        ak2, z02, fv2 = r2[0]
        print(f"\n  f(z0, alpha=2)   = {fv2:.6f}  (r21={R21_EXP})")
        print(f"  f/r21 - 1        = {fv2/R21_EXP - 1:+.6f}")
        print(f"  |f - r21|        = {abs(fv2-R21_EXP):.4f}")

    # Centrum [alpha*_1, alpha*_2]
    a1s = 2.44143051; a2s = 2.74175411
    center_alpha = (a1s + a2s) / 2
    print(f"\n  Centrum [alpha*_1, alpha*_2] = {center_alpha:.8f}")
    print(f"  |centrum - alpha_TGP=2| = {abs(center_alpha-2):.8f}")
    print(f"  Relacja: centrum/2 = {center_alpha/2:.6f}")

    # Czy alpha=2 jest specjalne wzgledem f(z0)?
    # f(z0, alpha) jest monotoniczna czy ma extremum?
    if len(fvals) >= 3:
        # Sprawdz gradient
        grads = np.diff(fvals) / np.diff(avals)
        print(f"\n  df/dalpha w roznych punktach:")
        for i in range(0, min(len(grads), 10), 2):
            a_mid = (avals[i]+avals[i+1])/2
            print(f"    alpha~{a_mid:.3f}: df/da = {grads[i]:+.3f}")

    # TESTY
    print()
    print("=" * 66)
    print("TESTY")
    print("=" * 66)

    T1 = len(crosses) == 2
    print(f"  T1: Dokladnie 2 przejscia f=r21: {'PASS' if T1 else 'FAIL'}  ({len(crosses)} przejsc)")

    T2 = False
    if r2:
        T2 = abs(fv2 - R21_EXP) > 10
        print(f"  T2: f(z0,alpha=2) != r21 (diff>10): {'PASS' if T2 else 'FAIL'}  (f={fv2:.4f})")

    T3 = abs(center_alpha - 2.0) > 0.2
    print(f"  T3: Centrum [a*1,a*2] != alpha_TGP=2: {'PASS' if T3 else 'FAIL'}  "
          f"(centrum={center_alpha:.5f}, diff={abs(center_alpha-2):.5f})")

    # T4: alpha=2 lezy POZA oknem [alpha*_1, alpha*_2]
    T4 = 2.0 < a1s
    print(f"  T4: alpha_TGP=2 < alpha*_1 (poza oknem): {'PASS' if T4 else 'FAIL'}")

    # T5: f(z0) jest monotonicznie rosnaca
    T5 = all(grads[i] > 0 for i in range(len(grads)) if not np.isnan(grads[i])) if len(fvals)>=3 else False
    print(f"  T5: f(z0,alpha) monotonicznie rosnaca: {'PASS' if T5 else 'FAIL'}")

    n_pass = sum([T1, T2, T3, T4, T5])
    print(f"\nWYNIK: {n_pass}/5 testow przeszlo")

    print()
    print("=" * 66)
    print("WNIOSEK (EX78)")
    print("=" * 66)
    print(f"  Centrum [alpha*_1, alpha*_2] = {center_alpha:.6f}")
    print(f"  alpha_TGP = 2.000000")
    print(f"  Roznica   = {abs(center_alpha-2):.6f}")
    print(f"  f(z0, alpha=2) = {fv2:.4f}  (r21={R21_EXP})")
    if abs(center_alpha - 2.0) > 0.3:
        print(f"  => alpha_TGP=2 NIE jest centrum [alpha*_1, alpha*_2].")
        print(f"     Ale alpha*_1 + alpha*_2 = 2pi-11/10 z 1ppm precyzja.")
    print("=" * 66)
