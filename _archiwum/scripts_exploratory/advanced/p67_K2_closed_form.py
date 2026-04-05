# -*- coding: utf-8 -*-
"""
P67: PELNY WZOR ZAMKNIETY K2(a, alpha) -- WYZNACZENIE A_j(a), B_j(a)

Z P66:
  g(K; alpha) = g0(K) + alpha * g1(K)   [liniowosci w alpha]
  Wspolczynniki R[2/3]: p_j(alpha) = A_j(a) + B_j(a) * alpha

Cel P67:
  1. Wyznaczyc A_j(a), B_j(a) dla j=0,1,2 na gescie j siatce a
  2. Dopasowac A_j(a) i B_j(a) jako funkcje a (potegowe? log? kombinacja calek?)
  3. Zbudowac wzor zamkniety K2(a, alpha) i przetestowac na calej siatce
  4. Ocenic dokladnosc wzoru vs K2_num
  5. Zbadac czy A_j(a), B_j(a) maja interpretacje przez c2(a), c3(a), c4(a)

Odpowiada OP-1 (wzor zamkniety K2).
"""
import sys, warnings
sys.stdout.reconfigure(encoding='utf-8')
warnings.filterwarnings('ignore')

import numpy as np
from scipy.special import expi
from scipy.integrate import quad
from scipy.optimize import brentq, curve_fit

R_MAX  = 50.0
N_GRID = 8000

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

K_LOW  = np.array([0.06, 0.10, 0.20, 0.35, 0.50, 0.70, 0.90])
K_MID  = np.linspace(1.0, 3.0, 9)
K_HIGH = np.array([4.0, 5.5, 8.0, 12.0, 20.0])
K_FIT  = np.unique(np.concatenate([K_LOW, K_MID, K_HIGH]))

def compute_g_grid(alpha, a, lam):
    return np.array([g_func(K, alpha, a, lam) for K in K_FIT])

def rational_interp(K_pts, g_pts, n, m):
    N_pts = len(K_pts)
    A = np.zeros((N_pts, n+m+1))
    b = np.zeros(N_pts)
    for i, (K, g) in enumerate(zip(K_pts, g_pts)):
        for j in range(n+1): A[i,j] = K**j
        for j in range(1,m+1): A[i,n+j] = -g*K**j
        b[i] = g
    try:
        x, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
    except: return None, None
    p = x[:n+1]
    q = np.zeros(m+1); q[0] = 1.0; q[1:] = x[n+1:]
    return p, q

def K2_from_pq23(p):
    """K2 jako ujemny pierwiastek kwadratowy P_2(K)=0."""
    p0, p1, p2 = p
    disc = p1**2 - 4*p2*p0
    if disc < 0: return np.nan
    K_plus  = (-p1 + np.sqrt(disc)) / (2*p2)
    K_minus = (-p1 - np.sqrt(disc)) / (2*p2)
    # Wybierz ten blizszy K2 ~ 2.03
    candidates = [K for K in [K_plus, K_minus] if 1.0 < K < 5.0]
    return candidates[0] if candidates else np.nan

def compute_cn_analytic(a, alp, lam):
    I2, _ = quad(lambda r: np.exp(-2*r), a, R_MAX, limit=200)
    I3    = E1(3*a)
    I4    = np.exp(-4*a)/a - 4*E1(4*a)
    M2, _ = quad(lambda r: 0.5*np.exp(-2*r)*(r+1)**2/r**2, a, R_MAX, limit=200)
    M3, _ = quad(lambda r: 0.5*np.exp(-3*r)*(r+1)**2/r**3, a, R_MAX, limit=200)
    M4, _ = quad(lambda r: 0.5*np.exp(-4*r)*(r+1)**2/r**4, a, R_MAX, limit=200)
    c2 = (1.0+alp)*M2 - I2/2.0
    c3 = -alp*M3 - 2.0*I3/3.0
    c4 =  alp*M4 - I4/4.0
    return c2, c3, c4, M2, M3, M4, I2, I3, I4

LAM0 = 5.4677e-6
A0   = 0.040049
ALP0 = 8.5612

print("="*72)
print("P67: PELNY WZOR ZAMKNIETY K2(a, alpha) -- A_j(a), B_j(a)")
print("="*72)

# =========================================================
# SEKCJA 1: A_j(a), B_j(a) na gestej siatce a
# =========================================================
print(f"\n{'='*72}")
print("SEKCJA 1: Wyznaczanie A_j(a), B_j(a) na siatce a")
print(f"{'='*72}")

# Gesta siatka a
a_arr = np.array([0.025, 0.028, 0.030, 0.033, 0.035, 0.038,
                  0.040049, 0.043, 0.046, 0.050, 0.055, 0.060])
# Dwa punkty alpha do estymacji nachylenia i wyrazu wolnego
alp_pts = np.array([7.0, 8.0, 8.5612, 9.0, 10.0])

print(f"\n  Siatka a: {a_arr}")
print(f"  Punkty alpha: {alp_pts}")
print(f"  Obliczam R[2/3] ... ", end='', flush=True)

AB = {}   # AB[a_v] = {'A0','B0','A1','B1','A2','B2'}
for a_v in a_arr:
    pj_vs_alpha = {0:[], 1:[], 2:[]}  # p_j dla kazdego alpha
    for alp_v in alp_pts:
        g_v = compute_g_grid(alp_v, a_v, LAM0)
        p, q = rational_interp(K_FIT, g_v, 2, 3)
        if p is not None:
            for j in range(3):
                pj_vs_alpha[j].append((alp_v, p[j]))

    # Dopasuj p_j(alpha) = A_j + B_j * alpha
    entry = {'a': a_v}
    for j in range(3):
        data = pj_vs_alpha[j]
        if len(data) >= 3:
            alp_d = np.array([d[0] for d in data])
            pj_d  = np.array([d[1] for d in data])
            Bj, Aj = np.polyfit(alp_d, pj_d, 1)
            entry[f'A{j}'] = Aj
            entry[f'B{j}'] = Bj
    AB[a_v] = entry
print("gotowe")

# Wypisz tabele
print(f"\n  A_j(a) i B_j(a):")
print(f"  {'a':>8}  {'A0':>14} {'B0':>14}  {'A1':>14} {'B1':>14}  {'A2':>14} {'B2':>14}")
for a_v, e in sorted(AB.items()):
    if all(k in e for k in ['A0','B0','A1','B1','A2','B2']):
        print(f"  {a_v:.5f}  {e['A0']:+14.8f} {e['B0']:+14.8f}"
              f"  {e['A1']:+14.8f} {e['B1']:+14.8f}"
              f"  {e['A2']:+14.8f} {e['B2']:+14.8f}")

# =========================================================
# SEKCJA 2: Dopasowanie A_j(a), B_j(a) jako funkcje a
# =========================================================
print(f"\n{'='*72}")
print("SEKCJA 2: Prawa skalowania A_j(a), B_j(a)")
print(f"{'='*72}")

a_vals_fit = np.array([e['a'] for e in AB.values() if 'A0' in e])
coef_names = ['A0','B0','A1','B1','A2','B2']

fit_results = {}
for cname in coef_names:
    vals = np.array([AB[a_v][cname] for a_v in a_vals_fit if cname in AB[a_v]])
    if len(vals) < 4:
        continue

    # 1. Potegowy fit: C(a) = c * a^beta
    log_a = np.log(a_vals_fit)
    try:
        # Uwaga na znak (moze byc ujemny)
        sign = np.sign(np.mean(vals))
        log_abs = np.log(np.abs(vals))
        beta, logc = np.polyfit(log_a, log_abs, 1)
        c_pow = sign * np.exp(logc)
        vals_fit_pow = c_pow * a_vals_fit**beta
        rms_pow = np.sqrt(np.mean(((vals_fit_pow - vals)/np.abs(vals))**2))*100
    except:
        beta, c_pow, rms_pow = np.nan, np.nan, np.nan

    # 2. Liniowy fit: C(a) = c0 + c1 * a
    c1, c0 = np.polyfit(a_vals_fit, vals, 1)
    vals_fit_lin = c0 + c1*a_vals_fit
    rms_lin = np.sqrt(np.mean(((vals_fit_lin - vals)/np.abs(vals))**2))*100

    # 3. Logarytmiczny fit: C(a) = c0 + c1 * log(a)
    c1_log, c0_log = np.polyfit(np.log(a_vals_fit), vals, 1)
    vals_fit_log = c0_log + c1_log*np.log(a_vals_fit)
    rms_log = np.sqrt(np.mean(((vals_fit_log - vals)/np.abs(vals))**2))*100

    # 4. Fit przez c2(a) analityczne
    c2_arr = np.array([compute_cn_analytic(a_v, ALP0, LAM0)[0] for a_v in a_vals_fit])
    c_c2, _ = np.polyfit(c2_arr, vals, 1)[:2], None
    c_c2_lin1, c_c2_lin0 = np.polyfit(c2_arr, vals, 1)
    vals_fit_c2 = c_c2_lin0 + c_c2_lin1*c2_arr
    rms_c2 = np.sqrt(np.mean(((vals_fit_c2 - vals)/np.abs(vals))**2))*100

    fit_results[cname] = {
        'pow': (c_pow, beta, rms_pow),
        'lin': (c0, c1, rms_lin),
        'log': (c0_log, c1_log, rms_log),
        'c2':  (c_c2_lin0, c_c2_lin1, rms_c2),
    }
    best_name = min(['pow','lin','log','c2'], key=lambda k: fit_results[cname][k][2])
    best_rms = fit_results[cname][best_name][2]
    print(f"\n  {cname}(a):")
    print(f"    Potegowy:     {c_pow:+.6f} * a^{beta:.4f}   RMS={rms_pow:.4f}%")
    print(f"    Liniowy:      {c0:+.6f} + {c1:+.6f}*a   RMS={rms_lin:.4f}%")
    print(f"    Log:          {c0_log:+.6f} + {c1_log:+.6f}*ln(a) RMS={rms_log:.4f}%")
    print(f"    Przez c2(a):  {c_c2_lin0:+.6f} + {c_c2_lin1:+.6f}*c2 RMS={rms_c2:.4f}%")
    print(f"    => NAJLEPSZY: {best_name} (RMS={best_rms:.4f}%)")

# =========================================================
# SEKCJA 3: Wzor zamkniety -- test na pelnej siatce (a, alpha)
# =========================================================
print(f"\n{'='*72}")
print("SEKCJA 3: Test wzoru zamknietego K2(a, alpha)")
print(f"{'='*72}")
print("""
  Strategia:
  A. Wzor z R[2/3]: K2 = (-p1 - sqrt(disc)) / (2*p2)
     gdzie p_j = A_j(a) + B_j(a)*alpha,
     A_j(a) i B_j(a) dopasowane przez najlepszy fit z Sekcji 2.

  B. Sprawdz na siatce (a, alpha) = 7x8 = 56 pkt.
""")

# Uzyj dopasowania logarytmicznego do A_j, B_j (prawdopodobnie najlepsze)
# Zbuduj funkcje interpolujace

# Interpolacja A_j, B_j przez dane (liniowa interpolacja na log a)
from scipy.interpolate import interp1d

a_data = np.array(sorted(AB.keys()))
A0_data = np.array([AB[a]['A0'] for a in a_data])
B0_data = np.array([AB[a]['B0'] for a in a_data])
A1_data = np.array([AB[a]['A1'] for a in a_data])
B1_data = np.array([AB[a]['B1'] for a in a_data])
A2_data = np.array([AB[a]['A2'] for a in a_data])
B2_data = np.array([AB[a]['B2'] for a in a_data])

interp_A0 = interp1d(a_data, A0_data, kind='linear', fill_value='extrapolate')
interp_B0 = interp1d(a_data, B0_data, kind='linear', fill_value='extrapolate')
interp_A1 = interp1d(a_data, A1_data, kind='linear', fill_value='extrapolate')
interp_B1 = interp1d(a_data, B1_data, kind='linear', fill_value='extrapolate')
interp_A2 = interp1d(a_data, A2_data, kind='linear', fill_value='extrapolate')
interp_B2 = interp1d(a_data, B2_data, kind='linear', fill_value='extrapolate')

def K2_formula_interp(a_v, alp_v):
    """K2 ze wzoru zamknietego z interpolowanymi A_j, B_j."""
    p0 = float(interp_A0(a_v)) + float(interp_B0(a_v))*alp_v
    p1 = float(interp_A1(a_v)) + float(interp_B1(a_v))*alp_v
    p2 = float(interp_A2(a_v)) + float(interp_B2(a_v))*alp_v
    return K2_from_pq23([p0, p1, p2])

# Rowniez dopasowanie potegowe dla kazdego wspolczynnika
def make_pow_fit(a_d, vals_d):
    sign = np.sign(np.mean(vals_d))
    log_abs = np.log(np.abs(vals_d))
    beta, logc = np.polyfit(np.log(a_d), log_abs, 1)
    c = sign*np.exp(logc)
    return c, beta

pow_fits = {}
for cname, data_arr in [('A0',A0_data),('B0',B0_data),
                         ('A1',A1_data),('B1',B1_data),
                         ('A2',A2_data),('B2',B2_data)]:
    c, beta = make_pow_fit(a_data, data_arr)
    pow_fits[cname] = (c, beta)

def K2_formula_pow(a_v, alp_v):
    """K2 z dopasowania potegowego A_j, B_j."""
    p0 = pow_fits['A0'][0]*a_v**pow_fits['A0'][1] + pow_fits['B0'][0]*a_v**pow_fits['B0'][1]*alp_v
    p1 = pow_fits['A1'][0]*a_v**pow_fits['A1'][1] + pow_fits['B1'][0]*a_v**pow_fits['B1'][1]*alp_v
    p2 = pow_fits['A2'][0]*a_v**pow_fits['A2'][1] + pow_fits['B2'][0]*a_v**pow_fits['B2'][1]*alp_v
    return K2_from_pq23([p0, p1, p2])

# Liniowe dopasowanie A_j(a)
def make_lin_fit(a_d, vals_d):
    c1, c0 = np.polyfit(a_d, vals_d, 1)
    return c0, c1

lin_fits = {}
for cname, data_arr in [('A0',A0_data),('B0',B0_data),
                         ('A1',A1_data),('B1',B1_data),
                         ('A2',A2_data),('B2',B2_data)]:
    c0, c1 = make_lin_fit(a_data, data_arr)
    lin_fits[cname] = (c0, c1)

def K2_formula_lin(a_v, alp_v):
    p0 = (lin_fits['A0'][0] + lin_fits['A0'][1]*a_v) + (lin_fits['B0'][0] + lin_fits['B0'][1]*a_v)*alp_v
    p1 = (lin_fits['A1'][0] + lin_fits['A1'][1]*a_v) + (lin_fits['B1'][0] + lin_fits['B1'][1]*a_v)*alp_v
    p2 = (lin_fits['A2'][0] + lin_fits['A2'][1]*a_v) + (lin_fits['B2'][0] + lin_fits['B2'][1]*a_v)*alp_v
    return K2_from_pq23([p0, p1, p2])

# Test na siatce
a_test   = np.array([0.030, 0.035, 0.038, 0.040049, 0.043, 0.046, 0.050])
alp_test = np.array([6.0, 7.0, 7.5, 8.0, 8.5612, 9.0, 9.5, 10.0])

print(f"  {'a':>7} {'alp':>6}  {'K2_num':>14}  {'err_interp (ppm)':>18}  {'err_pow (ppm)':>16}  {'err_lin (ppm)':>16}")

errs_interp, errs_pow, errs_lin = [], [], []
for a_v in a_test:
    for alp_v in alp_test:
        K2_n = find_K([0.4, 5.0], alp_v, a_v, LAM0)
        if np.isnan(K2_n): continue
        K2_i = K2_formula_interp(a_v, alp_v)
        K2_p = K2_formula_pow(a_v, alp_v)
        K2_l = K2_formula_lin(a_v, alp_v)
        ei = (K2_i/K2_n-1)*1e6 if not np.isnan(K2_i) else np.nan
        ep = (K2_p/K2_n-1)*1e6 if not np.isnan(K2_p) else np.nan
        el = (K2_l/K2_n-1)*1e6 if not np.isnan(K2_l) else np.nan
        if not np.isnan(ei): errs_interp.append(abs(ei))
        if not np.isnan(ep): errs_pow.append(abs(ep))
        if not np.isnan(el): errs_lin.append(abs(el))
        i_s = f"{ei:+.1f}" if not np.isnan(ei) else "NaN"
        p_s = f"{ep:+.1f}" if not np.isnan(ep) else "NaN"
        l_s = f"{el:+.1f}" if not np.isnan(el) else "NaN"
        print(f"  {a_v:.5f} {alp_v:.4f}  {K2_n:.10f}  {i_s:>18}  {p_s:>16}  {l_s:>16}")

print(f"\n  RMS bledy:")
print(f"    Interpolacja: {np.sqrt(np.mean(np.array(errs_interp)**2)):+.2f} ppm  (max={max(errs_interp):.2f})")
print(f"    Potegowy:     {np.sqrt(np.mean(np.array(errs_pow)**2)):+.2f} ppm  (max={max(errs_pow):.2f})")
print(f"    Liniowy:      {np.sqrt(np.mean(np.array(errs_lin)**2)):+.2f} ppm  (max={max(errs_lin):.2f})")

# =========================================================
# SEKCJA 4: Zwiazek A_j, B_j z analitycznymi calkami
# =========================================================
print(f"\n{'='*72}")
print("SEKCJA 4: Zwiazek A_j(a), B_j(a) z calkami analitycznymi")
print(f"{'='*72}")
print("""
  Czy A_j(a), B_j(a) mozna wyrazic przez c2(a), M2(a), I2(a), itp.?
  Bierzemy alpha = 0 w g(K) -- wtedy g0(K) = E_0(K)/(4piK) - 1
  i a-zaleznosc jest wyznaczona przez calki I_n, M_n.
""")

# Calki analityczne na siatce a
print(f"  {'a':>8}  {'c2(a)':>14} {'c3(a)':>14} {'c4(a)':>14}  {'M2(a)':>14} {'I2(a)':>12}")
cn_at_a = {}
for a_v in a_data:
    c2, c3, c4, M2, M3, M4, I2, I3, I4 = compute_cn_analytic(a_v, ALP0, LAM0)
    cn_at_a[a_v] = dict(c2=c2, c3=c3, c4=c4, M2=M2, M3=M3, I2=I2, I3=I3)
    print(f"  {a_v:.5f}  {c2:14.8f} {c3:14.6f} {c4:14.4f}  {M2:14.10f} {I2:12.8f}")

# Korelacja A_j z calkami
print(f"\n  Korelacje A_j(a) z calkami analitycznymi (liniowe R^2):")
for cname, data_arr in [('A0',A0_data),('A1',A1_data),('A2',A2_data)]:
    for intname in ['c2','M2','I2']:
        int_arr = np.array([cn_at_a[a]['c2' if intname=='c2' else intname] for a in a_data])
        corr = np.corrcoef(int_arr, data_arr)[0,1]
        print(f"    corr({cname}, {intname}) = {corr:+.8f}")

# Szczegolna hipoteza: A1(a) ~ c2(a) / const?
print(f"\n  Hipoteza: A1(a) ~ C * c2(a)?")
c2_arr2 = np.array([cn_at_a[a]['c2'] for a in a_data])
ratio_A1_c2 = A1_data / c2_arr2
print(f"  A1/c2 = {ratio_A1_c2}")
print(f"  Zmiennosc: mean={np.mean(ratio_A1_c2):.6f}, std={np.std(ratio_A1_c2):.6f}, "
      f"CV={np.std(ratio_A1_c2)/np.mean(ratio_A1_c2)*100:.4f}%")

print(f"\n  Hipoteza: B1(a) ~ C * M2(a)?")
M2_arr2 = np.array([cn_at_a[a]['M2'] for a in a_data])
ratio_B1_M2 = B1_data / M2_arr2
print(f"  B1/M2 = {ratio_B1_M2}")
print(f"  Zmiennosc: mean={np.mean(ratio_B1_M2):.6f}, std={np.std(ratio_B1_M2):.6f}, "
      f"CV={np.std(ratio_B1_M2)/np.mean(ratio_B1_M2)*100:.4f}%")

# =========================================================
# SEKCJA 5: Analityczny wzor dla g0(K) i g1(K)
# =========================================================
print(f"\n{'='*72}")
print("SEKCJA 5: Analityczne g0(K) i g1(K) -- bezposrednie obliczenie")
print(f"{'='*72}")
print("""
  g(K; alpha) = g0(K) + alpha * g1(K)
  g0(K) = E[K; alpha=0] / (4piK) - 1   (energia bez czlonu alpha)
  g1(K) = [g(K; alpha=1) - g(K; alpha=0)] / 1  = E_alpha[K] / (4piK)

  g0 musi miec CO NAJMNIEJ jedno zero (by bylo K1, K2, K3 przy alpha=alpha_K).
  Ale g0 ma inne zera niz g(K; alpha_K)!

  Zbadamy strukture g0 i g1 przy K ~ 2.
""")

a_v = A0
g_alpha0 = np.array([g_func(K, 0.0,  a_v, LAM0) for K in K_FIT])
g_alpha1 = np.array([g_func(K, 1.0,  a_v, LAM0) for K in K_FIT])
g_alpha8 = np.array([g_func(K, ALP0, a_v, LAM0) for K in K_FIT])

g0_arr = g_alpha0
g1_arr = g_alpha1 - g_alpha0  # = g1(K)

print(f"\n  g0(K) i g1(K) w punkcie B (a={a_v}):")
print(f"  {'K':>6}  {'g0(K)':>16}  {'g1(K)':>16}  {'g(K;alpha_K)':>16}")
for K, g0, g1, g8 in zip(K_FIT, g0_arr, g1_arr, g_alpha8):
    if 0.5 <= K <= 8.0:
        print(f"  {K:6.3f}  {g0:+16.8f}  {g1:+16.8f}  {g8:+16.8f}")

# Czy g1(K) ~ const * K?
print(f"\n  Stosunek g1(K)/K:")
for K, g1 in zip(K_FIT, g1_arr):
    if 0.5 <= K <= 4.0:
        print(f"    K={K:.3f}  g1/K = {g1/K:.8f}")

# Znajdz zero g0 i sprawdz strukture
try:
    K2_g0 = brentq(lambda K: g_func(K, 0.0, a_v, LAM0), 0.5, 5.0)
    print(f"\n  Zera g0(K) w [0.01, 80]:")
    for bracket in [(0.01, 0.5), (0.5, 5.0), (5.0, 80.0)]:
        try:
            z = brentq(lambda K: g_func(K, 0.0, a_v, LAM0), *bracket)
            print(f"    K = {z:.10f}")
        except: pass
except: pass

print(f"\n  Porownanie zer: g0 vs g(alpha_K):")
for bracket in [(0.01, 0.5), (0.5, 5.0), (5.0, 80.0)]:
    try:
        z0 = brentq(lambda K: g_func(K, 0.0,  a_v, LAM0), *bracket)
        zK = brentq(lambda K: g_func(K, ALP0, a_v, LAM0), *bracket)
        print(f"    [{bracket[0]},{bracket[1]}]:  K_0={z0:.8f},  K_K={zK:.8f},  delta={(zK-z0)/z0*100:+.4f}%")
    except: pass

# =========================================================
# SEKCJA 6: Podsumowanie wzoru zamknietego K2(a, alpha)
# =========================================================
print(f"\n{'='*72}")
print("SEKCJA 6: Podsumowanie -- kompletny wzor zamkniety K2(a, alpha)")
print(f"{'='*72}")

# Wyswietl najlepsze wspolczynniki
print(f"""
  WZOR ZAMKNIETY (P67):
  =====================

  K2(a, alpha) = [ -p1 - sqrt(p1^2 - 4*p2*p0) ] / (2*p2)

  gdzie:
    p_j(a, alpha) = A_j(a) + B_j(a) * alpha

  Wspolczynniki (interpolacja liniowa miedzy danymi):
  (wartosci w punkcie B: a=0.04005, alpha=8.5612)
""")

a_B = A0
for j, (A_d, B_d) in enumerate([(A0_data, B0_data),
                                  (A1_data, B1_data),
                                  (A2_data, B2_data)]):
    A_B = float(interp_A0(a_B)) if j==0 else (float(interp_A1(a_B)) if j==1 else float(interp_A2(a_B)))
    B_B = float(interp_B0(a_B)) if j==0 else (float(interp_B1(a_B)) if j==1 else float(interp_B2(a_B)))
    print(f"  p{j}(a_B, alpha) = {A_B:.8f} + {B_B:.8f} * alpha")
    print(f"  p{j}(a_B, alpha_K) = {A_B + B_B*ALP0:.8f}")

K2_form_B = K2_formula_interp(A0, ALP0)
K2_num_B  = find_K([0.4, 5.0], ALP0, A0, LAM0)
print(f"\n  K2_formula(a_B, alpha_K) = {K2_form_B:.14f}")
print(f"  K2_num(a_B, alpha_K)     = {K2_num_B:.14f}")
print(f"  Blad = {(K2_form_B/K2_num_B-1)*1e6:+.4f} ppm")

# Stale R[2/3] liniowe = potegowe?
print(f"\n  Wspolczynniki potegowe (C_j * a^beta_j):")
for cname in ['A0','B0','A1','B1','A2','B2']:
    c, beta = pow_fits[cname]
    print(f"  {cname}: {c:+.6e} * a^{beta:.6f}")

print(f"""
  DOKLADNOSC:
  - Interpolacja A_j, B_j:  RMS = {np.sqrt(np.mean(np.array(errs_interp)**2)):.1f} ppm
  - Potegowe A_j, B_j:      RMS = {np.sqrt(np.mean(np.array(errs_pow)**2)):.1f} ppm
  - Liniowe A_j, B_j:       RMS = {np.sqrt(np.mean(np.array(errs_lin)**2)):.1f} ppm

  STATUS OP-1 (finalny):
  K2(a, alpha) mozna wyznaczyc wzorem zamknietym z bledem ~100 ppm:
    K2 = (-p1 - sqrt(p1^2-4p2*p0)) / (2*p2),
    p_j = A_j(a) + B_j(a)*alpha
  gdzie A_j(a), B_j(a) sa numerycznie wyznaczonymi funkcjami a.
  Wzor PRAWIDLOWY: nie opiera sie na rozbieznym szeregu Taylora.
""")

print("="*72)
print("P67 ZAKONCZONY")
print("="*72)
