"""
EX72: ILOCZYN alpha*_1 * alpha*_2 I RELACJE ALGEBRAICZNE
STATUS: LEGACY-TRANSLATIONAL

This file continues the historical `alpha*` algebra/selection program built on
older variables and matching rules. Keep as legacy translation and exploratory
material rather than a canonical modern entry point.

Z ex71:
  alpha*_1 = 2.44142  (pierwsze zero delta)
  alpha*_2 = 2.74175  (drugie zero delta)
  Oba daja f(z0(alpha*)) = r21 = 206.768

Pytania:
  1. Czy alpha*_1 * alpha*_2 = dokladnie jakies stala?
     Propozycje: 2*pi+0.41?, pi+e+1?, 3*e-1?, phi^4?
  2. Czy (alpha*_1 + alpha*_2)/2 = phi^2/cos(something)?
  3. Czy sa kolejne zera delta dla alpha in [3, 5]?
  4. Czy z0(alpha*_1) * z0(alpha*_2) = phi^2?
  5. Jaka jest fizyczna roznica miedzy alpha*_1 a alpha*_2?
     (rozne sektory? rozne wartosci z0?)

Tests:
  T1: |alpha*_1 * alpha*_2 - 2*pi - C| < 0.01 dla prostego C
  T2: z0(alpha*_1) < z0(alpha*_2) (monotonicznie rosnie)
  T3: Czy sa kolejne zera? (f(z0, alpha)=r21 dla alpha in [3,5]?)
  T4: |z0*_1 * z0*_2 - phi^2| < 0.01 (iloczyn z0 = phi^2 ?)
  T5: alpha*_2 = alpha*_1 + delta_explicit (znana stala?)
"""
import sys, io
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

# Znane wartosci z ex71
ALPHA1 = 2.44142
ALPHA2 = 2.74175
Z0_A1  = 1.258516
Z0_A2  = 1.263349

def _integrate(g0, alpha, maxb=8):
    gb = np.exp(-1.0/(2.0*alpha)) + G_OFF
    def rhs(r,y):
        g,gp=y; g=max(g,gb+1e-7)
        fg=1.0+2.0*alpha*np.log(g)
        if abs(fg)<1e-10: return [gp,0.0]
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

def _A_inf(g0,alpha):
    r,g=_integrate(g0,alpha); Av,rv=[],[]
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

def _B_coeff(g0,alpha):
    r,g=_integrate(g0,alpha)
    _,B,_=_fit_win(r,g,WIN_LIST[2],WIN_LIST[2]+WIN_WIDTH)
    return B

def _find_z0(alpha,lo=1.05,hi=1.75):
    gs=np.linspace(lo,hi,40)
    Bs=[_B_coeff(g,alpha) for g in gs]
    for i in range(len(Bs)-1):
        if Bs[i]*Bs[i+1]<0:
            try: return brentq(lambda x:_B_coeff(x,alpha),gs[i],gs[i+1],xtol=1e-7)
            except: pass
    return None

def _f_ratio(g0,alpha):
    Ae=_A_inf(g0,alpha); Am=_A_inf(PHI*g0,alpha)
    return (Am/Ae)**4 if Ae>1e-6 else np.nan

def compute_fz0(alpha):
    """Oblicz f(z0(alpha), alpha) dla danego alpha."""
    z0 = _find_z0(alpha)
    if z0 is None: return alpha, np.nan, np.nan
    fz0 = _f_ratio(z0, alpha)
    return alpha, z0, fz0

if __name__ == '__main__':
    print("="*68)
    print("EX72: ILOCZYN alpha*_1 * alpha*_2 I KOLEJNE ZERA")
    print("="*68)
    print(f"  Uzywam {N_CORES} rdzeni")
    print()

    # ============================================================
    # CZESC 1: Analiza algebraiczna znanych wartosci
    # ============================================================
    print("-"*68)
    print("--- Czesc 1: Relacje algebraiczne alpha*_1 i alpha*_2 ---")
    print()

    a1, a2 = ALPHA1, ALPHA2
    z1, z2 = Z0_A1, Z0_A2

    print(f"  alpha*_1 = {a1:.5f}")
    print(f"  alpha*_2 = {a2:.5f}")
    print(f"  z0(a1)   = {z1:.6f}")
    print(f"  z0(a2)   = {z2:.6f}")
    print()

    product = a1 * a2
    summa   = a1 + a2
    diff    = a2 - a1
    geom    = np.sqrt(a1 * a2)

    print(f"  Suma:            {summa:.6f}")
    print(f"  Roznica:         {diff:.6f}")
    print(f"  Iloczyn:         {product:.6f}")
    print(f"  Srednia geom.:   {geom:.6f}")
    print(f"  Srednia arytm.:  {summa/2:.6f}")
    print()

    # Testy roznych kandydatow na iloczyn
    print("  Kandydaci na iloczyn a1*a2:")
    cands = {
        '2*pi':       2*np.pi,
        '2*pi+0.41':  2*np.pi + 0.41,
        'pi+e':       np.pi + np.e,
        '3*e-1':      3*np.e - 1,
        'phi^4':      PHI**4,
        'phi^3+phi':  PHI**3 + PHI,
        'e^2-0.7':    np.e**2 - 0.7,
        'e^2-2/3':    np.e**2 - 2/3,
        '6+2/3':      6 + 2/3,
        '6.69':       6.69,
        'pi+e-e/pi':  np.pi + np.e - np.e/np.pi,
        '2*pi+1/e':   2*np.pi + 1/np.e,
        'ln(1000)':   np.log(1000),
        '4*phi-1':    4*PHI - 1,
        '8/phi^2+0.5':8/PHI**2 + 0.5,
    }
    best_name, best_diff, best_val = None, 1e9, np.nan
    for name, val in sorted(cands.items(), key=lambda x: abs(x[1]-product)):
        d = abs(val - product)
        marker = " <---" if d < 0.01 else ""
        print(f"    {name:20s} = {val:.6f}   diff = {d:.6f}{marker}")
        if d < best_diff: best_diff, best_name, best_val = d, name, val

    print(f"\n  Najblizszy kandydat: {best_name} = {best_val:.6f}  (diff={best_diff:.6f})")

    # Kandydaci na sume
    print("\n  Kandydaci na sume a1+a2:")
    cands_s = {
        '5+1/5':    5.2,
        'pi+2':     np.pi + 2,
        '2*e-1/3':  2*np.e - 1/3,
        '4*phi/sqrt(3)': 4*PHI/np.sqrt(3),
        '5+phi-1':  4 + PHI,
        '3*phi':    3*PHI,
        '2*phi+2':  2*PHI + 2,
        'e+pi/2':   np.e + np.pi/2,
        'phi^3+1':  PHI**3 + 1,
        '5+1/phi':  5 + 1/PHI,
        '2*pi-1.1': 2*np.pi - 1.1,
        'sqrt(27)': np.sqrt(27),
    }
    for name, val in sorted(cands_s.items(), key=lambda x: abs(x[1]-summa)):
        d = abs(val - summa)
        marker = " <---" if d < 0.01 else ""
        print(f"    {name:20s} = {val:.6f}   diff = {d:.6f}{marker}")

    # Kandydaci na roznice
    print(f"\n  Roznica a2-a1 = {diff:.6f}")
    print(f"    1/(2*phi^2)   = {1/(2*PHI**2):.6f}  diff={abs(diff-1/(2*PHI**2)):.6f}")
    print(f"    1/pi          = {1/np.pi:.6f}  diff={abs(diff-1/np.pi):.6f}")
    print(f"    phi/5         = {PHI/5:.6f}  diff={abs(diff-PHI/5):.6f}")
    print(f"    1/(2*e)       = {1/(2*np.e):.6f}  diff={abs(diff-1/(2*np.e)):.6f}")
    print(f"    3/10          = 0.300000  diff={abs(diff-0.3):.6f}")

    # ============================================================
    # CZESC 2: z0 wartosci przy obu alpha*
    # ============================================================
    print("\n"+"-"*68)
    print("--- Czesc 2: Relacje z0(alpha*_1) i z0(alpha*_2) ---")
    print()

    print(f"  z0(a1) = {z1:.6f}")
    print(f"  z0(a2) = {z2:.6f}")
    print(f"  z0(a2) - z0(a1) = {z2-z1:.6f}")
    print(f"  z0(a2) / z0(a1) = {z2/z1:.6f}")
    print(f"  z0(a1) * z0(a2) = {z1*z2:.6f}")
    print(f"  phi^2           = {PHI**2:.6f}  diff={abs(z1*z2-PHI**2):.6f}")
    print(f"  (z1+z2)/2       = {(z1+z2)/2:.6f}")
    print(f"  sqrt(z1*z2)     = {np.sqrt(z1*z2):.6f}")
    print(f"  (phi^2+1)/2     = {(PHI**2+1)/2:.6f}  diff={abs((z1+z2)/2-(PHI**2+1)/2):.6f}")
    print(f"  z0(a2)/z0(a1)-1 = {z2/z1-1:.6f}")
    print(f"  1/(a1*a2)^0.5   = {1/np.sqrt(a1*a2):.6f}")
    print(f"  G_ghost(a=2)    = {np.exp(-1/(2*2.0)):.6f}  (g* at alpha=2)")

    # ============================================================
    # CZESC 3: Szukaj kolejnych zer f(z0,alpha)=r21 dla alpha in [3,5]
    # ============================================================
    print("\n"+"-"*68)
    print("--- Czesc 3: Kolejne zera f(z0,alpha)=r21 dla alpha in [3,5] ---")
    print()

    alpha_ext = np.array([3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.3, 4.6, 5.0])
    print(f"  Obliczam {len(alpha_ext)} punktow alpha (rownolegle)...")

    with Pool(processes=N_CORES) as pool:
        results_ext = pool.map(compute_fz0, alpha_ext)

    print(f"\n  {'alpha':>6}  {'z0':>10}  {'f(z0)':>10}  {'f-r21':>9}")
    print(f"  {'-'*6}  {'-'*10}  {'-'*10}  {'-'*9}")
    fz0_ext = []
    for alpha, z0, fz0 in sorted(results_ext):
        z0_s  = f"{z0:.6f}" if z0 else "N/A"
        fz0_s = f"{fz0:.4f}" if not np.isnan(fz0) else "N/A"
        fds   = f"{fz0-R21_EXP:+.4f}" if not np.isnan(fz0) else "N/A"
        marker = " <r21?" if (not np.isnan(fz0) and abs(fz0-R21_EXP)<10) else ""
        print(f"  {alpha:6.2f}  {z0_s:>10}  {fz0_s:>10}  {fds:>9}{marker}")
        fz0_ext.append(fz0)

    # Czy sa przejscia przez r21?
    crossings_ext = 0
    for i in range(len(fz0_ext)-1):
        f1, f2 = fz0_ext[i], fz0_ext[i+1]
        if not(np.isnan(f1) or np.isnan(f2)):
            if (f1-R21_EXP)*(f2-R21_EXP)<0:
                crossings_ext += 1
                print(f"  -> Crossing r21 miedzy alpha={alpha_ext[i]:.2f} i {alpha_ext[i+1]:.2f}")
    if crossings_ext == 0:
        print(f"  Brak przekroczen r21 w [3, 5]")

    # ============================================================
    # CZESC 4: Precyzyjniejsze alpha*_1 i alpha*_2 (z ex71)
    # ============================================================
    print("\n"+"-"*68)
    print("--- Czesc 4: Precyzyjna weryfikacja alpha*_1 i alpha*_2 ---")
    print()

    # Oblicz f(z0) z wysoka precyzja przy alpha*_1 i alpha*_2
    for alpha_name, alpha_val in [("alpha*_1", ALPHA1), ("alpha*_2", ALPHA2)]:
        alpha_c, z0_c, fz0_c = compute_fz0(alpha_val)
        print(f"  {alpha_name} = {alpha_val}:")
        print(f"    z0   = {z0_c:.7f}" if z0_c else "    z0 = N/A")
        print(f"    f(z0) = {fz0_c:.6f}  (r21={R21_EXP})" if not np.isnan(fz0_c) else "    f = N/A")
        print(f"    diff  = {fz0_c - R21_EXP:.6f}" if not np.isnan(fz0_c) else "")
        print()

    # ============================================================
    # TESTY
    # ============================================================
    print("="*68)
    print("TESTY")
    print("="*68)

    # T1: iloczyn bliski prostej stalej
    T1 = best_diff < 0.01
    print(f"  T1: a1*a2 = prosta stala (<0.01): {'PASS' if T1 else 'FAIL'}  "
          f"(best: {best_name}={best_val:.5f}, diff={best_diff:.5f})")

    # T2: z0(a2) > z0(a1)
    T2 = z2 > z1
    print(f"  T2: z0(a2) > z0(a1):              {'PASS' if T2 else 'FAIL'}  "
          f"({z2:.5f} > {z1:.5f})")

    # T3: brak kolejnych zer w [3,5]
    T3 = crossings_ext == 0
    print(f"  T3: Brak kolejnych zer w [3,5]:   {'PASS' if T3 else 'FAIL'}  "
          f"(crossings={crossings_ext})")

    # T4: z0(a1)*z0(a2) bliskie phi^2
    T4 = abs(z1*z2 - PHI**2) < 0.02
    print(f"  T4: z0*z0 ≈ phi^2:                {'PASS' if T4 else 'FAIL'}  "
          f"(z1*z2={z1*z2:.5f}, phi^2={PHI**2:.5f}, diff={abs(z1*z2-PHI**2):.5f})")

    # T5: a2 = a1 + prosta stala
    best_diff_diff = min(abs(diff - v) for v in [0.3, 1/np.pi, PHI/5, 1/(2*np.e)])
    T5 = best_diff_diff < 0.005
    print(f"  T5: a2-a1 = prosta stala (<0.005): {'PASS' if T5 else 'FAIL'}  "
          f"(diff={diff:.6f})")

    n_pass = sum([T1, T2, T3, T4, T5])
    print(f"\nWYNIK: {n_pass}/5 testow przeszlo")

    print("\n"+"="*68)
    print("WNIOSEK (EX72)")
    print("="*68)
    print(f"  alpha*_1 * alpha*_2 = {product:.6f}")
    print(f"  alpha*_1 + alpha*_2 = {summa:.6f}")
    print(f"  alpha*_2 - alpha*_1 = {diff:.6f}")
    print(f"  z0(a1)*z0(a2)       = {z1*z2:.6f}  (phi^2={PHI**2:.6f})")
    print(f"  Najblizszy dla iloczyn: {best_name}={best_val:.5f} (diff={best_diff:.5f})")
    if crossings_ext > 0:
        print(f"  Kolejne zera f(z0)=r21: TAK (alpha in [3,5])")
    else:
        print(f"  Kolejne zera f(z0)=r21: NIE w [3,5]")
    print("="*68)
