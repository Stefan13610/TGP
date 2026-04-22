"""
EX71: MAPA delta(alpha) — DWA ZERA (ROWNOLEGLE, max 16 rdzeni)
STATUS: LEGACY-TRANSLATIONAL

This script belongs to the historical `delta(alpha)`, `alpha*`, and old
selection-language chain. Read it as legacy/exploratory transition material,
not as a canonical current `nbody` source.
"""
import sys, io, os
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, curve_fit
from multiprocessing import Pool, cpu_count
import warnings
warnings.filterwarnings('ignore')

PHI     = (1.0 + np.sqrt(5.0)) / 2.0
R21_EXP = 206.768
G_OFF   = 0.005
R_MAX   = 100.0
WIN_LIST  = [16, 22, 28, 36, 46, 58, 72]
WIN_WIDTH = 14.0
N_CORES   = min(16, cpu_count())

# ============================================================
# SOLVER
# ============================================================
def _integrate(g0, alpha, maxb=8):
    gb = np.exp(-1.0/(2.0*alpha)) + G_OFF
    def rhs(r, y):
        g, gp = y
        g = max(g, gb + 1e-7)
        fg = 1.0 + 2.0*alpha*np.log(g)
        if abs(fg)<1e-10: return [gp, 0.0]
        dr=g**2*(1.0-g); cr=(alpha/g)*gp**2
        if r<1e-10: return [gp,(dr-cr)/(3.0*fg)]
        return [gp,(dr-cr-fg*2.0*gp/r)/fg]
    def ev(r,y): return y[0]-gb
    ev.terminal=True; ev.direction=-1
    y0=[g0,0.0]; r0=1e-10; ra,ga=[],[]
    for _ in range(maxb+1):
        sol=solve_ivp(rhs,(r0,R_MAX),y0,events=ev,
                      dense_output=True,rtol=1e-9,atol=1e-11,max_step=0.05)
        ra.append(sol.t); ga.append(sol.y[0])
        if sol.status==1 and len(sol.t_events[0])>0:
            rh=sol.t_events[0][-1]; st=sol.sol(rh)
            y0=[st[0],-st[1]]; r0=rh
        else: break
    r=np.concatenate(ra); g=np.concatenate(ga)
    idx=np.argsort(r); return r[idx],g[idx]

def _fit_win(r,g,rL,rR):
    mask=(r>=rL)&(r<=rR)
    if mask.sum()<20: return 0.,0.,0.
    rf=r[mask]; df=(g[mask]-1.0)*rf
    M=np.column_stack([np.cos(rf),np.sin(rf)])
    bc=np.linalg.lstsq(M,df,rcond=None)[0]
    return float(np.sqrt(bc[0]**2+bc[1]**2)),float(bc[0]),float(bc[1])

def _A_inf(g0, alpha):
    r,g=_integrate(g0,alpha)
    Av,rv=[],[]
    for rL in WIN_LIST:
        if rL+WIN_WIDTH>r[-1]: break
        A,_,_=_fit_win(r,g,rL,rL+WIN_WIDTH)
        if A>0.002: Av.append(A); rv.append(float(rL))
    if not Av: return 0.
    if len(Av)<3: return Av[-1]
    try:
        p,_=curve_fit(lambda x,ai,a:ai*(1+a/x),rv,Av,p0=[Av[-1],0.],maxfev=2000)
        return float(p[0])
    except: return float(Av[-1])

def _B_coeff(g0, alpha):
    r,g=_integrate(g0,alpha)
    _,B,_=_fit_win(r,g,WIN_LIST[2],WIN_LIST[2]+WIN_WIDTH)
    return B

def _find_z0(alpha, lo=1.05, hi=1.70):
    gs=np.linspace(lo,hi,40)
    Bs=[_B_coeff(g,alpha) for g in gs]
    for i in range(len(Bs)-1):
        if Bs[i]*Bs[i+1]<0:
            try: return brentq(lambda x:_B_coeff(x,alpha),gs[i],gs[i+1],xtol=1e-6)
            except: pass
    return None

def _f_ratio(g0, alpha):
    Ae=_A_inf(g0,alpha); Am=_A_inf(PHI*g0,alpha)
    return (Am/Ae)**4 if Ae>1e-6 else np.nan

def _find_g0star(alpha, target=R21_EXP, lo=1.05, hi=1.80):
    gs=np.linspace(lo,hi,25)
    fs=[_f_ratio(g,alpha) for g in gs]
    for i in range(len(fs)-1):
        if not(np.isnan(fs[i]) or np.isnan(fs[i+1])):
            if (fs[i]-target)*(fs[i+1]-target)<0:
                try: return brentq(lambda x:_f_ratio(x,alpha)-target,gs[i],gs[i+1],xtol=1e-6)
                except: pass
    return None

# ============================================================
# FUNKCJA ROBOCZA dla kazdego alpha (wywolywana rownoleglem)
# ============================================================
def compute_alpha_point(alpha):
    """Oblicz z0, g0*, delta, f(z0) dla jednego alpha."""
    z0  = _find_z0(alpha)
    g0s = _find_g0star(alpha)
    d   = (g0s - z0) if (z0 and g0s) else np.nan
    fz0 = _f_ratio(z0, alpha) if z0 else np.nan
    return alpha, z0, g0s, d, fz0

# ============================================================
# GLOWNY PROGRAM
# ============================================================
if __name__ == '__main__':
    print("="*68)
    print("EX71: MAPA delta(alpha) — DWA ZERA")
    print("="*68)
    print(f"  Uzywam {N_CORES} rdzeni (dostepne: {cpu_count()})")
    print()

    # Punkty do obliczen
    alpha_arr = np.array([
        1.6, 1.8, 2.0, 2.1, 2.2, 2.3, 2.35, 2.40, 2.42, 2.44,
        2.45, 2.47, 2.50, 2.55, 2.60, 2.65, 2.70, 2.72, 2.73,
        2.74, 2.75, 2.77, 2.80, 2.90, 3.0
    ])

    print(f"  Obliczam {len(alpha_arr)} punktow alpha...")
    with Pool(processes=N_CORES) as pool:
        results = pool.map(compute_alpha_point, alpha_arr)

    # Sortuj wyniki
    results.sort(key=lambda x: x[0])

    # Wydruk
    print()
    print(f"  {'alpha':>6}  {'z0':>10}  {'g0*':>10}  {'delta':>9}  {'f(z0)':>9}  {'f-r21':>8}")
    print(f"  {'-'*6}  {'-'*10}  {'-'*10}  {'-'*9}  {'-'*9}  {'-'*8}")

    alpha_v  = []
    delta_v  = []
    fz0_v    = []
    z0_v     = []
    g0s_v    = []

    for alpha, z0, g0s, d, fz0 in results:
        alpha_v.append(alpha); delta_v.append(d)
        fz0_v.append(fz0); z0_v.append(z0); g0s_v.append(g0s)
        z0s  = f"{z0:.6f}"  if z0         else "N/A"
        g0ss = f"{g0s:.6f}" if g0s        else "N/A"
        ds   = f"{d:.6f}"   if not np.isnan(d)   else "N/A"
        fzs  = f"{fz0:.4f}" if not np.isnan(fz0) else "N/A"
        fds  = f"{fz0-R21_EXP:+.4f}" if not np.isnan(fz0) else "N/A"
        marker = ""
        if not np.isnan(d) and abs(d)<0.002: marker = " <ZERO"
        print(f"  {alpha:6.3f}  {z0s:>10}  {g0ss:>10}  {ds:>9}  {fzs:>9}  {fds:>8}{marker}")

    alpha_v = np.array(alpha_v)
    delta_v = np.array(delta_v)
    fz0_v   = np.array(fz0_v)

    # ============================================================
    # Szukaj zer
    # ============================================================
    print("\n"+"-"*68)
    print("--- Zera delta(alpha) ---")
    print()

    zeros_found = []
    for i in range(len(alpha_v)-1):
        d1, d2 = delta_v[i], delta_v[i+1]
        if not(np.isnan(d1) or np.isnan(d2)) and d1*d2<0:
            try:
                # Lokalna funkcja delta (nie mozna uzyc closure w pool,
                # ale tutaj robimy brentq sekwencyjnie na malym przedziale)
                def _delta_local(a, _a1=alpha_v[i], _a2=alpha_v[i+1]):
                    _, _, _, d, _ = compute_alpha_point(a)
                    return d
                z = brentq(_delta_local, alpha_v[i], alpha_v[i+1], xtol=1e-4)
                zeros_found.append(z)
                _, z0_z, g0s_z, d_z, fz0_z = compute_alpha_point(z)
                print(f"  Zero #{len(zeros_found)}: alpha = {z:.5f}")
                print(f"    z0(alpha)   = {z0_z:.6f}" if z0_z else "    z0 = N/A")
                print(f"    g0*(alpha)  = {g0s_z:.6f}" if g0s_z else "    g0* = N/A")
                print(f"    delta       = {d_z:.8f}" if not np.isnan(d_z) else "    delta = N/A")
                print(f"    f(z0)       = {fz0_z:.4f}  (r21={R21_EXP})" if not np.isnan(fz0_z) else "    f(z0) = N/A")
                print()
            except Exception as e:
                print(f"  Blad [{alpha_v[i]:.3f},{alpha_v[i+1]:.3f}]: {e}")

    # ============================================================
    # Analiza zer
    # ============================================================
    if len(zeros_found) >= 2:
        z1, z2 = zeros_found[0], zeros_found[1]
        print("-"*68)
        print("--- Analiza zer ---")
        print()
        print(f"  zero_1           = {z1:.5f}")
        print(f"  zero_2           = {z2:.5f}")
        print(f"  zero_2/zero_1    = {z2/z1:.5f}  (phi={PHI:.5f})")
        print(f"  phi*zero_1       = {PHI*z1:.5f}  (vs zero_2={z2:.5f})")
        print(f"  zero_1+zero_2    = {z1+z2:.5f}")
        print(f"  zero_1*zero_2    = {z1*z2:.5f}")
        print(f"  sqrt(z1*z2)      = {np.sqrt(z1*z2):.5f}")
        print(f"  (z1+z2)/2        = {(z1+z2)/2:.5f}")
        print(f"  |z2 - sqrt(6)|   = {abs(z2-np.sqrt(6)):.5f}  (sqrt6={np.sqrt(6):.5f})")
        print(f"  |z2 - phi^2|     = {abs(z2-PHI**2):.5f}  (phi^2={PHI**2:.5f})")
        print(f"  |z2/z1 - phi|    = {abs(z2/z1-PHI):.5f}")

    # ============================================================
    # TESTY
    # ============================================================
    print("\n"+"="*68)
    print("TESTY")
    print("="*68)

    T1 = len(zeros_found) >= 2
    print(f"  T1: Dwa zera delta:  {'PASS' if T1 else 'FAIL'}  ({len(zeros_found)} zer)")

    T2 = False
    if zeros_found:
        _, z0_z1, _, _, fz0_z1 = compute_alpha_point(zeros_found[0])
        T2 = not np.isnan(fz0_z1) and abs(fz0_z1-R21_EXP)/R21_EXP < 0.001
        print(f"  T2: f(z0(zero_1))=r21 (0.1%):  {'PASS' if T2 else 'FAIL'}  "
              f"(f={fz0_z1:.4f})" if not np.isnan(fz0_z1) else "  T2: FAIL")

    T3 = False
    if len(zeros_found) >= 2:
        _, z0_z2, _, _, fz0_z2 = compute_alpha_point(zeros_found[1])
        T3 = not np.isnan(fz0_z2) and abs(fz0_z2-R21_EXP)/R21_EXP > 0.01
        print(f"  T3: f(z0(zero_2))!=r21:         {'PASS' if T3 else 'FAIL'}  "
              f"(f={fz0_z2:.4f})" if not np.isnan(fz0_z2) else "  T3: FAIL")

    T4 = False
    if len(zeros_found) >= 2:
        ratio_z = zeros_found[1]/zeros_found[0]
        diff_phi = abs(ratio_z - PHI)
        T4 = diff_phi < 0.05
        print(f"  T4: zero_2/zero_1 ≈ phi:        {'PASS' if T4 else 'FAIL'}  "
              f"(ratio={ratio_z:.5f}, phi={PHI:.5f}, diff={diff_phi:.5f})")

    T5 = False
    valid_fz0 = fz0_v[~np.isnan(fz0_v)]
    cross = sum(1 for i in range(len(valid_fz0)-1)
                if (valid_fz0[i]-R21_EXP)*(valid_fz0[i+1]-R21_EXP)<0)
    T5 = cross >= 2
    print(f"  T5: f(z0)=r21 dwukrotnie:       {'PASS' if T5 else 'FAIL'}  "
          f"(przekroczenia={cross})")

    n_pass = sum([T1, T2, T3, T4, T5])
    print(f"\nWYNIK: {n_pass}/5 testow przeszlo")

    print("\n"+"="*68)
    print("WNIOSEK (EX71)")
    print("="*68)
    if len(zeros_found) >= 2:
        z1, z2 = zeros_found
        print(f"  zero_1 = {z1:.5f}")
        print(f"  zero_2 = {z2:.5f}")
        print(f"  zero_2/zero_1 = {z2/z1:.5f}  vs phi={PHI:.5f}")
    elif zeros_found:
        print(f"  zero_1 = {zeros_found[0]:.5f}")
    print("="*68)
