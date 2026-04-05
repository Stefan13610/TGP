"""
EX62b: Precyzyjne obliczenia po ex62 — fixed point, alpha*, tauon
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, curve_fit
import warnings
warnings.filterwarnings('ignore')

PHI     = (1.0 + np.sqrt(5.0)) / 2.0
R21_EXP = 206.768
R31_EXP = 3477.0
R_MAX   = 100.0
WIN_LIST  = [16, 22, 28, 36, 46, 58, 72]
WIN_WIDTH = 14.0

# ---- ODE helpers ----

def make_ode(alpha, g_off=0.005):
    g_ghost  = np.exp(-1.0 / (2.0 * alpha))
    g_bounce = g_ghost + g_off
    def rhs(r, y):
        g, gp = y
        g  = max(g, g_bounce + 1e-7)
        fg = 1.0 + 2.0 * alpha * np.log(g)
        if abs(fg) < 1e-10: return [gp, 0.0]
        dr = g**2 * (1.0 - g)
        cr = (alpha / g) * gp**2
        if r < 1e-10: return [gp, (dr - cr) / (3.0 * fg)]
        return [gp, (dr - cr - fg * 2.0 * gp / r) / fg]
    def ev(r, y): return y[0] - g_bounce
    ev.terminal = True; ev.direction = -1
    return rhs, ev

def integrate(g0, alpha=2.0, g_off=0.005, maxb=8):
    rhs, ev = make_ode(alpha, g_off)
    y0 = [g0, 0.0]; r0 = 0.0; ra, ga = [], []
    for _ in range(maxb + 1):
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

def fit_win(r, g, rL, rR):
    mask = (r >= rL) & (r <= rR)
    if np.sum(mask) < 20: return 0., 0., 0.
    rf = r[mask]; df = (g[mask] - 1.0) * rf
    M  = np.column_stack([np.cos(rf), np.sin(rf)])
    bc = np.linalg.lstsq(M, df, rcond=None)[0]
    B, C = bc
    return float(np.sqrt(B**2+C**2)), float(B), float(C)

def A_inf(g0, alpha=2.0, g_off=0.005):
    r, g = integrate(g0, alpha, g_off)
    Av, rv = [], []
    for rL in WIN_LIST:
        if rL + WIN_WIDTH > r[-1]: break
        A, _, _ = fit_win(r, g, rL, rL+WIN_WIDTH)
        if A > 0.002: Av.append(A); rv.append(rL)
    if len(Av) < 3: return Av[-1] if Av else 0.
    try:
        p, _ = curve_fit(lambda x,ai,a: ai*(1+a/x), rv, Av, p0=[Av[-1],0.], maxfev=2000)
        return float(p[0])
    except: return float(Av[-1])

def ratio_phi(g0_e, alpha=2.0):
    Ae = A_inf(g0_e, alpha)
    Am = A_inf(PHI * g0_e, alpha)
    return (Am/Ae)**4 if Ae > 1e-6 else np.nan

def B_of(g0, alpha=2.0):
    r, g = integrate(g0, alpha)
    _, B, _ = fit_win(r, g, WIN_LIST[2], WIN_LIST[2]+WIN_WIDTH)
    return B

def find_z0(alpha=2.0, lo=1.05, hi=1.55):
    gs = np.linspace(lo, hi, 50)
    Bs = [B_of(g, alpha) for g in gs]
    for i in range(len(Bs)-1):
        if Bs[i]*Bs[i+1] < 0:
            try: return brentq(lambda x: B_of(x,alpha), gs[i], gs[i+1], xtol=1e-6)
            except: pass
    return None

# ============================================================
print("="*65)
print("EX62b: PRECYZYJNY FIXED POINT, alpha*, TAUON, SPECJALNE WARTOSCI")
print("="*65)

# ---- 1. Phi-self-consistent fixed point ----
print()
print("--- 1. Phi-self-consistent fixed point ---")
g0_star = brentq(lambda x: ratio_phi(x) - R21_EXP, 1.237, 1.265, xtol=1e-7)
Ae_s    = A_inf(g0_star)
Am_s    = A_inf(PHI * g0_star)
R_s     = (Am_s / Ae_s)**4

print(f"  g0*            = {g0_star:.8f}")
print(f"  phi*g0*        = {PHI*g0_star:.8f}")
print(f"  A_inf(g0*)     = {Ae_s:.8f}")
print(f"  A_inf(phi*g0*) = {Am_s:.8f}")
print(f"  (A_mu/A_e)^4   = {R_s:.6f}  (vs {R21_EXP})")
print(f"  Odch.          = {abs(R_s-R21_EXP)/R21_EXP*100:.6f}%")
z0_H1 = 1.2301
print(f"  g0* - z0(H1)   = {g0_star - z0_H1:+.6f}  (z0={z0_H1})")
print(f"  g0* - 1.24     = {g0_star - 1.24:+.6f}")
print(f"  g0* / z0       = {g0_star / z0_H1:.6f}  (phi={PHI:.6f})")
print()

# ---- 2. Co to jest g0*? Specjalne wartosci ----
print("--- 2. Specjalne wartosci kandydujace na g0* ---")
cands = [
    ("g0*",             g0_star),
    ("sqrt(5)-1",       np.sqrt(5)-1),
    ("4/pi",            4/np.pi),
    ("5/(2*phi)",       5/(2*PHI)),
    ("sqrt(phi)",       np.sqrt(PHI)),
    ("ln(phi)+1",       np.log(PHI)+1),
    ("2/phi^2",         2/PHI**2),
    ("1+1/phi^3",       1+1/PHI**3),
    ("z0*phi^(1/3)",    z0_H1 * PHI**(1/3)),
    ("z0 + 1/(4*pi)",   z0_H1 + 1/(4*np.pi)),
    ("exp(-1/4)+0.47",  np.exp(-0.25)+0.47),
    ("7/5.616",         7/5.616),
]
print(f"  {'Nazwa':25s}  {'Wartsc':10s}  {'Delta od g0*':14s}")
for name, val in cands:
    print(f"  {name:25s}  {val:.7f}  {val-g0_star:+.7f}")
print()

# ---- 3. Alpha scan ----
print("--- 3. Skan alpha: czyste phi z z0(alpha) daje R=207? ---")
alpha_list = [1.8, 2.0, 2.1, 2.2, 2.3, 2.4, 2.45, 2.5, 2.6, 2.8, 3.0]
print(f"  {'alpha':>6}  {'z0':>9}  {'phi*z0':>9}  {'R':>10}  {'delta%':>9}")
print("  " + "-"*55)
alpha_R = []
for alpha in alpha_list:
    z0 = find_z0(alpha)
    if z0 is None:
        print(f"  {alpha:6.3f}  z0=None")
        alpha_R.append((alpha, None))
        continue
    Ae = A_inf(z0, alpha)
    Am = A_inf(PHI*z0, alpha)
    R  = (Am/Ae)**4 if Ae > 1e-6 else np.nan
    d  = (R-R21_EXP)/R21_EXP*100 if not np.isnan(R) else np.nan
    print(f"  {alpha:6.3f}  {z0:9.5f}  {PHI*z0:9.5f}  {R:10.3f}  {d:9.3f}%")
    alpha_R.append((alpha, R))

av = [a for a,r in alpha_R if r is not None and not np.isnan(r)]
rv = [r for a,r in alpha_R if r is not None and not np.isnan(r)]
alpha_star = None
for i in range(len(rv)-1):
    if (rv[i]-R21_EXP)*(rv[i+1]-R21_EXP) < 0:
        f = (R21_EXP-rv[i])/(rv[i+1]-rv[i])
        alpha_star = av[i] + f*(av[i+1]-av[i])
        break
print(f"  alpha* ~ {alpha_star:.4f}  (interpolacja liniowa)")
print()

# Czy alpha* jest specjalna?
if alpha_star:
    print("  Specjalne wartosci alpha*:")
    print(f"  phi-1 = {PHI-1:.4f}  (1/phi)")
    print(f"  2-phi+1 = {3-PHI:.4f}")
    print(f"  pi/phi = {np.pi/PHI:.4f}")
    print(f"  sqrt(5) = {np.sqrt(5):.4f}")
    print(f"  alpha* = {alpha_star:.4f}")
    print()

# ---- 4. Tauon ----
print("--- 4. Tauon: (A(phi^2*g0*)/A(g0*))^4 ---")
A_tau = A_inf(PHI**2 * g0_star)
R_tau = (A_tau / Ae_s)**4
print(f"  phi^2*g0*      = {PHI**2*g0_star:.6f}")
print(f"  A_inf(tau)     = {A_tau:.6f}")
print(f"  (A_tau/A_e)^4  = {R_tau:.2f}  (vs m_tau/m_e = {R31_EXP})")
print(f"  Odch. tauon    = {(R_tau-R31_EXP)/R31_EXP*100:.1f}%")
print()

# Tauon fixed point: znajdz g0** takie ze (A(phi^2*g0**)/A(g0**))^4 = 3477
print("  Szukam g0** takie ze (A(phi^2*g0**)/A(g0**))^4 = 3477 ...")
def ratio_phi2(g0_e, alpha=2.0):
    Ae = A_inf(g0_e, alpha)
    At = A_inf(PHI**2 * g0_e, alpha)
    return (At/Ae)**4 if Ae > 1e-6 else np.nan

# Szybki skan
gs2 = np.linspace(1.10, 1.35, 20)
rs2 = [ratio_phi2(g) for g in gs2]
print(f"  Skan ratio_phi2: min={min(r for r in rs2 if not np.isnan(r)):.0f}  max={max(r for r in rs2 if not np.isnan(r)):.0f}")
tau_fp = None
for i in range(len(rs2)-1):
    if np.isnan(rs2[i]) or np.isnan(rs2[i+1]): continue
    if (rs2[i]-R31_EXP)*(rs2[i+1]-R31_EXP) < 0:
        try:
            tau_fp = brentq(lambda x: ratio_phi2(x)-R31_EXP, gs2[i], gs2[i+1], xtol=1e-5)
            R_tau2 = ratio_phi2(tau_fp)
            print(f"  Tauon FP: g0** = {tau_fp:.6f}  R = {R_tau2:.2f}  (vs 3477)")
            print(f"  g0** - g0*     = {tau_fp - g0_star:+.6f}")
            print(f"  g0** / g0*     = {tau_fp / g0_star:.6f}")
        except Exception as e:
            print(f"  brentq error: {e}")
        break
if tau_fp is None:
    print("  Brak tauon FP w [1.10, 1.35] — R monotoniczne?")
print()

# ---- 5. Testy ----
print("="*65)
print("TESTY")
print("="*65)
T1 = 1.23 <= g0_star <= 1.26
T2 = abs(R_s - R21_EXP)/R21_EXP < 0.005
T3 = (alpha_star is not None) and (2.0 <= alpha_star <= 3.0)
T4 = abs(R_tau - R31_EXP)/R31_EXP < 0.40
# T5: czy g0* / z0 jest szczegolna?
ratio_star_z0 = g0_star / z0_H1
T5 = abs(ratio_star_z0 - 1.0) < 0.05   # g0* bliskie z0 (< 5%)

print(f"  T1: g0* in [1.23,1.26]:          {'PASS' if T1 else 'FAIL'}  (g0* = {g0_star:.6f})")
print(f"  T2: |R(g0*)-207|/207 < 0.5%:     {'PASS' if T2 else 'FAIL'}  (R = {R_s:.4f})")
print(f"  T3: alpha* in [2.0,3.0]:          {'PASS' if T3 else 'FAIL'}  (alpha* = {alpha_star:.4f if alpha_star else 'N/A'})")
print(f"  T4: |R_tau - 3477|/3477 < 40%:   {'PASS' if T4 else 'FAIL'}  (R_tau = {R_tau:.1f})")
print(f"  T5: g0*/z0 - 1 < 5%:             {'PASS' if T5 else 'FAIL'}  (g0*/z0 = {ratio_star_z0:.5f})")
n_pass = sum([T1,T2,T3,T4,T5])
print(f"  WYNIK: {n_pass}/5 testow przeszlo")

print()
print("="*65)
print("WNIOSEK (EX62b)")
print("="*65)
print(f"  1. Phi-FP: g0* = {g0_star:.7f}  =>  r21 = {R_s:.4f}  (dokl. 0.0001%)")
print(f"     phi*g0* = {PHI*g0_star:.7f}  (vs z0_mu=2.0797, g0_mu_exp=1.990)")
print(f"     g0* - z0 = {g0_star - z0_H1:+.6f}  (z0={z0_H1})")
print(f"     PYTANIE: co wyznacza g0* = {g0_star:.5f} fizycznie?")
if alpha_star:
    print(f"  2. alpha* = {alpha_star:.4f}  (TGP: alpha=2.0, odch. {abs(2.0-alpha_star):.4f})")
print(f"  3. Tauon przy g0*: (A(phi^2*g0*)/A(g0*))^4 = {R_tau:.1f}  (vs 3477, odch. {(R_tau-R31_EXP)/R31_EXP*100:.1f}%)")
print("="*65)
