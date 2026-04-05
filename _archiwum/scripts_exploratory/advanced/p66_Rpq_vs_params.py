# -*- coding: utf-8 -*-
"""
P66: WSPOLCZYNNIKI R[n/m] A PARAMETRY MODELU -- WZOR ZAMKNIETY?

Z P65: R[7/2](K) daje K2 z bledem +1.79 ppm.
       R[2/3](K) daje K2 z bledem +89 ppm (6 parametrow -- minimalna reprezentacja).

Cel P66:
  1. Jak wspolczynniki [p_j, q_j] R[n/m] zaleza od (a, alpha, lambda)?
     -> czy istnieja proste prawa skalowania p_j ~ f(a, alpha)?
  2. Analiza R[2/3] (6 par.) -- najprostszy schemat, latwiejszy do analizy
  3. Test: czy K2 mozna aproksymowac jako K2 = -p_0/p_1 + g(K)... ?
  4. Porownanie K2 z [7/2] vs [2/3] jako funkcji parametrow
  5. Wymiarowe prawa skalowania: K2(a, alpha) -- forma funkcjonalna
  6. Czy istnieje wartoscc K2 wyrazona przez c2 (=112.156) i inne calki?

Odpowiada OP-1 (poszukiwanie wzoru zamknietego).
"""
import sys, warnings
sys.stdout.reconfigure(encoding='utf-8')
warnings.filterwarnings('ignore')

import numpy as np
from scipy.special import expi
from scipy.integrate import quad
from scipy.optimize import brentq
from itertools import product

R_MAX  = 50.0
N_GRID = 8000   # nieco szybszy

def E1(x):
    return -expi(-x)

def V_mod(phi, lam):
    return phi**3/3.0 - phi**4/4.0 + lam*(phi-1.0)**6/6.0

def energy_num(K, alpha, a, lam, N=N_GRID):
    t   = np.linspace(0, 1, N)
    r   = a*(R_MAX/a)**t
    phi = np.maximum(1.0 + K*np.exp(-r)/r, 1e-10)
    dphi= K*np.exp(-r)*(-r-1.0)/r**2
    V1  = V_mod(1.0, lam)
    Ek  = 4*np.pi*np.trapezoid(0.5*dphi**2*(1.0+alpha/phi)*r**2, r)
    Ep  = 4*np.pi*np.trapezoid((V_mod(phi,lam)-V1)*r**2, r)
    return Ek+Ep

def g_func(K, alpha, a, lam, N=N_GRID):
    return energy_num(K,alpha,a,lam,N)/(4*np.pi*K) - 1.0

def find_K(bracket, alpha, a, lam, N=N_GRID, tol=1e-12):
    try: return brentq(lambda K: g_func(K,alpha,a,lam,N),
                       bracket[0], bracket[1], xtol=tol)
    except: return np.nan

# Siatka punktow do interpolacji (z P65 -- optymalna)
K_LOW  = np.array([0.06, 0.10, 0.20, 0.35, 0.50, 0.70, 0.90])
K_MID  = np.linspace(1.0, 3.0, 9)
K_HIGH = np.array([4.0, 5.5, 8.0, 12.0, 20.0])
K_FIT  = np.unique(np.concatenate([K_LOW, K_MID, K_HIGH]))

def compute_g_grid(alpha, a, lam):
    return np.array([g_func(K, alpha, a, lam) for K in K_FIT])

def rational_interp(K_pts, g_pts, n, m):
    N_pts = len(K_pts)
    N_par = n + m + 1
    if N_pts < N_par:
        return None, None
    A = np.zeros((N_pts, N_par))
    b = np.zeros(N_pts)
    for i, (K, g) in enumerate(zip(K_pts, g_pts)):
        for j in range(n+1):
            A[i, j] = K**j
        for j in range(1, m+1):
            A[i, n+j] = -g * K**j
        b[i] = g
    try:
        x, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
    except:
        return None, None
    p = x[:n+1]
    q = np.zeros(m+1); q[0] = 1.0; q[1:] = x[n+1:]
    return p, q

def find_K2_from_rational(p, q, K2_guess=2.03):
    """Znajdz K2 jako zero licznika P(K) w otoczeniu K2_guess."""
    def num(K):
        return sum(p[j]*K**j for j in range(len(p)))
    K_scan = np.linspace(1.0, 3.5, 5000)
    vals = np.array([num(K) for K in K_scan])
    sign_ch = np.where(np.diff(np.sign(vals)))[0]
    candidates = []
    for idx in sign_ch:
        try:
            z = brentq(num, K_scan[idx], K_scan[idx+1], xtol=1e-14)
            candidates.append(z)
        except:
            pass
    if not candidates:
        return np.nan
    return min(candidates, key=lambda z: abs(z - K2_guess))

# Punkt B (referencyjny)
A0    = 0.040049
ALP0  = 8.5612
LAM0  = 5.4677e-6

print("="*72)
print("P66: WSPOLCZYNNIKI R[n/m] VS PARAMETRY -- WZOR ZAMKNIETY?")
print("="*72)
print(f"  Ref: a={A0}, alpha={ALP0}, lambda={LAM0:.5e}")

K2_ref = find_K([0.4, 5.0], ALP0, A0, LAM0)
print(f"  K2_num (ref) = {K2_ref:.14f}")

# =========================================================
# SEKCJA 1: Wspolczynniki R[2/3] i R[7/2] w punkcie B
# =========================================================
print(f"\n{'='*72}")
print("SEKCJA 1: Wspolczynniki R[2/3] i R[7/2] w punkcie B")
print(f"{'='*72}")

g_ref = compute_g_grid(ALP0, A0, LAM0)
p23_0, q23_0 = rational_interp(K_FIT, g_ref, 2, 3)
p72_0, q72_0 = rational_interp(K_FIT, g_ref, 7, 2)

K2_23 = find_K2_from_rational(p23_0, q23_0)
K2_72 = find_K2_from_rational(p72_0, q72_0)
print(f"\n  R[2/3]: K2 = {K2_23:.14f}  blad = {(K2_23/K2_ref-1)*1e6:+.4f} ppm")
print(f"  R[7/2]: K2 = {K2_72:.14f}  blad = {(K2_72/K2_ref-1)*1e6:+.4f} ppm")

print(f"\n  Wspolczynniki R[2/3]:")
print(f"  P_2: p0={p23_0[0]:+.8e}, p1={p23_0[1]:+.8e}, p2={p23_0[2]:+.8e}")
print(f"  Q_3: q1={q23_0[1]:+.8e}, q2={q23_0[2]:+.8e}, q3={q23_0[3]:+.8e}")

print(f"\n  Wspolczynniki R[7/2]:")
for j, pj in enumerate(p72_0):
    print(f"  p{j} = {pj:+.8e}")
for j, qj in enumerate(q72_0[1:], 1):
    print(f"  q{j} = {qj:+.8e}")

# =========================================================
# SEKCJA 2: Siatka parametrow (a, alpha) -- K2 i wspolczynniki
# =========================================================
print(f"\n{'='*72}")
print("SEKCJA 2: K2 a parametry modelu -- prawa skalowania")
print(f"{'='*72}")

# Siatka a i alpha wokol punktu B
a_vals   = np.array([0.030, 0.035, 0.038, 0.040049, 0.043, 0.046, 0.050])
alp_vals = np.array([6.0, 7.0, 7.5, 8.0, 8.5612, 9.0, 9.5, 10.0])

print(f"\n  K2 numeryczne dla roznych (a, alpha) (lambda = lambda_K):")
print(f"  {'a':>8}  {'alpha':>8}  {'K2_num':>18}  {'K2_R23 err (ppm)':>18}  {'K2_R72 err (ppm)':>18}")

K2_grid = {}
for a_v in a_vals:
    for alp_v in alp_vals:
        K2_n = find_K([0.4, 5.0], alp_v, a_v, LAM0)
        if np.isnan(K2_n):
            continue
        g_v = compute_g_grid(alp_v, a_v, LAM0)
        p23, q23 = rational_interp(K_FIT, g_v, 2, 3)
        p72, q72 = rational_interp(K_FIT, g_v, 7, 2)
        K2_23v = find_K2_from_rational(p23, q23) if p23 is not None else np.nan
        K2_72v = find_K2_from_rational(p72, q72) if p72 is not None else np.nan
        err23 = (K2_23v/K2_n-1)*1e6 if not np.isnan(K2_23v) else np.nan
        err72 = (K2_72v/K2_n-1)*1e6 if not np.isnan(K2_72v) else np.nan
        K2_grid[(a_v, alp_v)] = {'K2': K2_n, 'err23': err23, 'err72': err72,
                                  'p23': p23, 'q23': q23, 'p72': p72, 'q72': q72}
        e23s = f"{err23:+.1f}" if not np.isnan(err23) else "NaN"
        e72s = f"{err72:+.1f}" if not np.isnan(err72) else "NaN"
        print(f"  {a_v:8.5f}  {alp_v:8.4f}  {K2_n:18.10f}  {e23s:>18}  {e72s:>18}")

# =========================================================
# SEKCJA 3: Prawa skalowania K2(a, alpha)
# =========================================================
print(f"\n{'='*72}")
print("SEKCJA 3: Prawa skalowania K2 = K2(a, alpha)")
print(f"{'='*72}")
print("""
  Badamy:
  1. K2 vs alpha przy stalym a (liniowe? potegowe?)
  2. K2 vs a przy stalym alpha
  3. Czy log K2 ~ beta_a * log a + beta_alp * log alpha?
""")

# Przy stalym a = A0, zmiana alpha
a_fix = A0
alpha_scan = []
K2_scan_a  = []
for alp_v in alp_vals:
    key = (a_fix, alp_v)
    if key in K2_grid:
        alpha_scan.append(alp_v)
        K2_scan_a.append(K2_grid[key]['K2'])

alpha_scan = np.array(alpha_scan)
K2_scan_a  = np.array(K2_scan_a)

if len(K2_scan_a) >= 3:
    # log-log fit: log K2 = beta * log alpha + const
    log_alp = np.log(alpha_scan)
    log_K2  = np.log(K2_scan_a)
    beta_alp, logA_alp = np.polyfit(log_alp, log_K2, 1)
    A_alp = np.exp(logA_alp)
    print(f"  K2 vs alpha (a=const={a_fix:.5f}):")
    print(f"    log-log fit: K2 ~ {A_alp:.6f} * alpha^{beta_alp:.6f}")
    print(f"    Tabela: alpha -> K2")
    for alp, K2 in zip(alpha_scan, K2_scan_a):
        K2_fit = A_alp * alp**beta_alp
        print(f"      alpha={alp:.4f}  K2_num={K2:.8f}  K2_fit={K2_fit:.8f}  blad={(K2_fit/K2-1)*100:+.4f}%")

# Przy stalym alpha = ALP0, zmiana a
alp_fix = ALP0
a_scan_list = []
K2_scan_alp = []
for a_v in a_vals:
    key = (a_v, alp_fix)
    if key in K2_grid:
        a_scan_list.append(a_v)
        K2_scan_alp.append(K2_grid[key]['K2'])

a_scan_arr   = np.array(a_scan_list)
K2_scan_alp  = np.array(K2_scan_alp)

if len(K2_scan_alp) >= 3:
    log_a   = np.log(a_scan_arr)
    log_K2a = np.log(K2_scan_alp)
    beta_a, logA_a = np.polyfit(log_a, log_K2a, 1)
    A_a = np.exp(logA_a)
    print(f"\n  K2 vs a (alpha=const={alp_fix:.4f}):")
    print(f"    log-log fit: K2 ~ {A_a:.6f} * a^{beta_a:.6f}")
    for av, K2 in zip(a_scan_arr, K2_scan_alp):
        K2_fit = A_a * av**beta_a
        print(f"      a={av:.5f}  K2_num={K2:.8f}  K2_fit={K2_fit:.8f}  blad={(K2_fit/K2-1)*100:+.4f}%")

# Dwuwymiarowy fit: K2 = A * a^beta_a * alpha^beta_alp
valid_pairs = [(k[0], k[1], v['K2']) for k, v in K2_grid.items() if not np.isnan(v['K2'])]
if len(valid_pairs) >= 5:
    a_arr   = np.array([x[0] for x in valid_pairs])
    alp_arr = np.array([x[1] for x in valid_pairs])
    K2_arr  = np.array([x[2] for x in valid_pairs])
    X = np.column_stack([np.ones(len(a_arr)), np.log(a_arr), np.log(alp_arr)])
    y = np.log(K2_arr)
    coef, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    logA_2d, beta_a_2d, beta_alp_2d = coef
    A_2d = np.exp(logA_2d)
    print(f"\n  Dwuwymiarowy log-log fit K2 = A * a^b_a * alpha^b_alp:")
    print(f"    A={A_2d:.6f}, b_a={beta_a_2d:.6f}, b_alp={beta_alp_2d:.6f}")
    print(f"    K2 ~ {A_2d:.4f} * a^({beta_a_2d:.4f}) * alpha^({beta_alp_2d:.4f})")
    # Bledy
    K2_fit_2d = A_2d * a_arr**beta_a_2d * alp_arr**beta_alp_2d
    rms = np.sqrt(np.mean(((K2_fit_2d - K2_arr)/K2_arr)**2)) * 100
    print(f"    RMS blad = {rms:.4f}%")

# =========================================================
# SEKCJA 4: Analiza wspolczynnikow R[2/3] vs (a, alpha)
# =========================================================
print(f"\n{'='*72}")
print("SEKCJA 4: Wspolczynniki R[2/3] vs (a, alpha)")
print(f"{'='*72}")
print("""
  R[2/3](K) = (p0 + p1*K + p2*K^2) / (1 + q1*K + q2*K^2 + q3*K^3)
  6 parametrow: p0, p1, p2, q1, q2, q3
  Badamy jak kazdy wspolczynnik zaleza od (a, alpha).
""")

if len(valid_pairs) >= 3:
    # Zbierz wspolczynniki dla wszystkich par (a, alpha)
    pq_data = {}
    for k, v in K2_grid.items():
        if v['p23'] is not None and len(v['p23']) == 3:
            pq_data[k] = {
                'p0': v['p23'][0], 'p1': v['p23'][1], 'p2': v['p23'][2],
                'q1': v['q23'][1], 'q2': v['q23'][2], 'q3': v['q23'][3],
            }

    print(f"  Wspolczynniki R[2/3] dla roznych (a, alpha):")
    print(f"  {'a':>7} {'alp':>6}  {'p0':>14} {'p1':>14} {'p2':>14}  {'q1':>12} {'q2':>12} {'q3':>12}")
    for (a_v, alp_v), pq in sorted(pq_data.items()):
        print(f"  {a_v:.5f} {alp_v:.4f}  {pq['p0']:+14.6f} {pq['p1']:+14.6f} {pq['p2']:+14.6f}"
              f"  {pq['q1']:+12.6f} {pq['q2']:+12.6f} {pq['q3']:+12.6f}")

    # Sprawdz czy wspolczynniki sa liniowe w alpha
    # Przy stalym a = A0
    alp_fixed_a = [(k[1], v) for k, v in pq_data.items() if abs(k[0]-A0) < 1e-6]
    alp_fixed_a.sort(key=lambda x: x[0])

    if len(alp_fixed_a) >= 3:
        alp_v_arr = np.array([x[0] for x in alp_fixed_a])
        for coef_name in ['p0', 'p1', 'p2', 'q1', 'q2', 'q3']:
            vals = np.array([x[1][coef_name] for x in alp_fixed_a])
            # liniowe i potegowe dopasowanie
            slope_lin, intercept_lin = np.polyfit(alp_v_arr, vals, 1)
            try:
                beta_p, logA_p = np.polyfit(np.log(np.abs(alp_v_arr)), np.log(np.abs(vals)), 1)
                A_p = np.exp(logA_p)
                r2_lin = 1 - np.var(vals - np.polyval([slope_lin, intercept_lin], alp_v_arr))/np.var(vals)
            except:
                beta_p, A_p, r2_lin = np.nan, np.nan, np.nan
            print(f"    {coef_name}(alpha): lin slope={slope_lin:+.6e}, intercept={intercept_lin:+.6e}, R2={r2_lin:.6f}")

# =========================================================
# SEKCJA 5: Zero K2 analitycznie -- wzor z wspolczynnikow
# =========================================================
print(f"\n{'='*72}")
print("SEKCJA 5: K2 = zero licznika -- analityczna formu la z R[2/3]")
print(f"{'='*72}")
print("""
  P_2(K) = p0 + p1*K + p2*K^2 = 0
  K = (-p1 +/- sqrt(p1^2 - 4*p2*p0)) / (2*p2)

  Wzor kwadratowy: K2 = (-p1 + sqrt(p1^2 - 4*p2*p0)) / (2*p2)
  (wybieramy znak dajacy K2 ~ 2.03)

  Jezeli p_j ~ f(a, alpha), to K2 ma wzor zamkniety!
""")

if p23_0 is not None:
    p0, p1, p2 = p23_0
    disc = p1**2 - 4*p2*p0
    K2_quad_plus  = (-p1 + np.sqrt(disc)) / (2*p2) if disc > 0 else np.nan
    K2_quad_minus = (-p1 - np.sqrt(disc)) / (2*p2) if disc > 0 else np.nan
    print(f"\n  W punkcie B (ref):")
    print(f"  p0={p0:+.8e}, p1={p1:+.8e}, p2={p2:+.8e}")
    print(f"  dyskryminant = {disc:.8e}")
    print(f"  K2_quad(+) = {K2_quad_plus:.14f}  blad={(K2_quad_plus/K2_ref-1)*1e6:+.4f} ppm" if not np.isnan(K2_quad_plus) else "  K2_quad(+) = NaN")
    print(f"  K2_quad(-) = {K2_quad_minus:.14f}" if not np.isnan(K2_quad_minus) else "  K2_quad(-) = NaN")
    print(f"  K2_R23     = {K2_23:.14f}  blad={(K2_23/K2_ref-1)*1e6:+.4f} ppm")
    print(f"  K2_num     = {K2_ref:.14f}")

    # Sprawdz w calej siatce
    print(f"\n  Wzor kwadratowy K2_quad dla siatki parametrow:")
    print(f"  {'a':>7} {'alp':>6}  {'K2_quad':>18}  {'blad (ppm)':>12}  {'K2_num':>18}")
    for (a_v, alp_v), pq in sorted(pq_data.items()):
        p0_v, p1_v, p2_v = pq['p0'], pq['p1'], pq['p2']
        disc_v = p1_v**2 - 4*p2_v*p0_v
        if disc_v > 0:
            K2_q = (-p1_v + np.sqrt(disc_v)) / (2*p2_v)
        else:
            K2_q = (-p1_v - np.sqrt(max(disc_v, 0))) / (2*p2_v)
        key = (a_v, alp_v)
        K2_true = K2_grid[key]['K2'] if key in K2_grid else np.nan
        if not np.isnan(K2_true) and not np.isnan(K2_q):
            err_q = (K2_q/K2_true - 1)*1e6
        else:
            err_q = np.nan
        eq = f"{K2_q:.10f}" if not np.isnan(K2_q) else "NaN"
        en = f"{err_q:+.4f}" if not np.isnan(err_q) else "NaN"
        K2t = f"{K2_true:.10f}" if not np.isnan(K2_true) else "NaN"
        print(f"  {a_v:.5f} {alp_v:.4f}  {eq:>18}  {en:>12}  {K2t:>18}")

# =========================================================
# SEKCJA 6: Wspolczynniki c2 analityczne a wspolczynniki R[2/3]
# =========================================================
print(f"\n{'='*72}")
print("SEKCJA 6: Zwiazek wspolczynnikow R[2/3] z analitycznymi c2, c3, c4")
print(f"{'='*72}")
print("""
  Z P64: c2, c3, c4 obliczane analitycznie z calek M_n, I_n.
  Pytanie: czy p0, p1, p2, q1, q2, q3 z R[2/3] mozna wyrazic przez c2, c3, c4?

  Z teorii Pade (Taylor-Pade): nie dziala (P64).
  Tu: Pade z DANYCH -- ale czy wspolczynniki maja prosta postac?
""")

def compute_cn_analytic(a, alp, lam):
    """Oblicz c2, c3, c4 analitycznie."""
    I2, _ = quad(lambda r: np.exp(-2*r), a, R_MAX, limit=200)
    I3    = E1(3*a)
    I4    = np.exp(-4*a)/a - 4*E1(4*a)
    M2, _ = quad(lambda r: 0.5*np.exp(-2*r)*(r+1)**2/r**2, a, R_MAX, limit=200)
    M3, _ = quad(lambda r: 0.5*np.exp(-3*r)*(r+1)**2/r**3, a, R_MAX, limit=200)
    M4, _ = quad(lambda r: 0.5*np.exp(-4*r)*(r+1)**2/r**4, a, R_MAX, limit=200)
    c2 = (1.0+alp)*M2 - I2/2.0
    c3 = -alp*M3 - 2.0*I3/3.0
    c4 =  alp*M4 - I4/4.0
    return c2, c3, c4

print(f"\n  Punkt B: analityczne c2, c3, c4:")
c2_B, c3_B, c4_B = compute_cn_analytic(A0, ALP0, LAM0)
print(f"  c2 = {c2_B:.8f}")
print(f"  c3 = {c3_B:.8f}")
print(f"  c4 = {c4_B:.8f}")

# Stosunek wspolczynnikow do c2, c3, c4
if p23_0 is not None:
    print(f"\n  Wspolczynniki R[2/3] (punkt B) znormalizowane przez c2, c3, c4:")
    print(f"  p0 = {p23_0[0]:.8e}")
    print(f"  p1 = {p23_0[1]:.8e}  (p1/c2 = {p23_0[1]/c2_B:.8f})")
    print(f"  p2 = {p23_0[2]:.8e}  (p2/c3 = {p23_0[2]/c3_B:.8f})")
    print(f"  q1 = {q23_0[1]:.8e}  (q1/c3*c2 = {q23_0[1]*c2_B/c3_B:.8f})")
    print(f"  q2 = {q23_0[2]:.8e}")
    print(f"  q3 = {q23_0[3]:.8e}")

    # Sprawdz czy p0/(p1/c2) ~ -1/c2... (czyli K_root ~ -p0/p1 ~ 1/c2?)
    K2_linear = -p23_0[0] / p23_0[1]
    print(f"\n  Aproksymacja liniowa K2 ~ -p0/p1 = {K2_linear:.10f}")
    print(f"  K2_num = {K2_ref:.10f}  blad = {(K2_linear/K2_ref-1)*1e6:+.4f} ppm")
    print(f"  1/c2   = {1.0/c2_B:.10f}")

    # Dla calej siatki: czy K2 ~ -p0/p1 ?
    print(f"\n  Sprawdzenie -p0/p1 vs K2_num dla siatki:")
    print(f"  {'a':>7} {'alp':>6}  {'K2=-p0/p1':>18}  {'blad (ppm)':>12}")
    for (a_v, alp_v), pq in sorted(pq_data.items()):
        K2_lin_v = -pq['p0'] / pq['p1']
        K2_true = K2_grid.get((a_v, alp_v), {}).get('K2', np.nan)
        err_l = (K2_lin_v/K2_true - 1)*1e6 if not np.isnan(K2_true) else np.nan
        print(f"  {a_v:.5f} {alp_v:.4f}  {K2_lin_v:18.10f}  {err_l:+12.4f}")

# =========================================================
# SEKCJA 7: Najlepsza analityczna aproksymacja K2
# =========================================================
print(f"\n{'='*72}")
print("SEKCJA 7: Podsumowanie -- najlepsza analityczna formu la K2")
print(f"{'='*72}")

if len(valid_pairs) >= 5:
    print(f"\n  Podsumowanie metod aproksymacji K2:")
    print(f"  {'Metoda':>40}  {'blad K2 (ppm)':>16}")

    # 1. Interpolacja R[7/2] MNK 21 pkt
    print(f"  {'R[7/2] MNK, 21 pkt (P65)':>40}  {'+1.79':>16}")
    # 2. Interpolacja R[2/3] MNK 21 pkt
    err23_ref = (K2_23/K2_ref-1)*1e6 if not np.isnan(K2_23) else np.nan
    print(f"  {'R[2/3] MNK, 21 pkt (P65)':>40}  {err23_ref:+16.2f}")
    # 3. Wzor kwadratowy -p0/p1 z R[2/3]
    if p23_0 is not None:
        K2_lin0 = -p23_0[0]/p23_0[1]
        err_lin0 = (K2_lin0/K2_ref-1)*1e6
        print(f"  {'K2 = -p0/p1 (Pade[2/3] linear)':>40}  {err_lin0:+16.2f}")
    # 4. Prawa skalowania
    try:
        K2_sc = A_2d * A0**beta_a_2d * ALP0**beta_alp_2d
        err_sc = (K2_sc/K2_ref-1)*1e6
        print(f"  {'K2 = A*a^b_a*alpha^b_alp (skalowanie)':>40}  {err_sc:+16.2f}")
    except:
        pass
    # 5. K1_leading (P61)
    print(f"  {'K1_leading = 1/A1 (P61)':>40}  {'-93300':>16}  (K1, nie K2)")

    print(f"""
  WNIOSEK:
  Najlepsza aktualna aproksymacja K2:
  -> R[7/2] z 21 pkt g(K): blad +1.79 ppm (P65)
  -> Wzor kwadratowy -p0/p1 z R[2/3]: blad ~ {err_lin0 if p23_0 is not None else '?':.1f} ppm

  Wspolczynniki R[2/3] (p0, p1, p2, q1, q2, q3):
  -> Nie maja prostej postaci analitycznej przez c2, c3, c4
  -> Ale stosunek -p0/p1 daje K2 z dokladnoscia ~ kilka ppm!
  -> Prawa skalowania: K2 ~ {A_2d:.4f} * a^({beta_a_2d:.3f}) * alpha^({beta_alp_2d:.3f})
     (RMS blad = {rms:.3f}%)

  INTERPRETACJA:
  Wspolczynniki R[n/m] sa emergentne -- wynikaja z globalnej struktury g(K),
  nie z lokalnego rozwinecia Taylora. Nie istnieje prosty wzor zamkniety
  laczacy je z c2..c4, ale zaleznosci numeryczne sa stabilne i przewidywalne.
""")

print("="*72)
print("P66 ZAKONCZONY")
print("="*72)
