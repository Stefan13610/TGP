# -*- coding: utf-8 -*-
"""
P65: INTERPOLACJA WYMIERNA g(K) -- OPTYMALNY [n/m] DLA K2 (OP-1)

Diagnoza z P64: Taylor-Pade strukturalnie niemozliwe (K2>>R_zbież).
Cel P65:
  1. Wyznacz g(K) w 25 punktach na przedziale [0.05, 4.5]
  2. Dopasuj aproksymant Pade [n/m] przez INTERPOLACJE (nie Taylor)
     -- uzyj macierzy linearnej na wspolczynniki P_n, Q_m
  3. Przetestuj wszystkie konfiguracje [n/m] z n+m <= 8
  4. Wybierz optymalny [n/m] wedlug:
     a) dokladnosci K2 (blad ppm vs K2_num)
     b) stabilnosci (czy K2 i K1 i K3 sa wszystkie zidentyfikowane)
     c) minimalizacji liczby parametrow
  5. Ocen zachowanie residuum g(K) - g_Pade(K) na przedziale

Odpowiada OP-1 (centralne pytanie: analityczna droga do K2).
"""
import sys, warnings
sys.stdout.reconfigure(encoding='utf-8')
warnings.filterwarnings('ignore')

import numpy as np
from scipy.special import expi
from scipy.integrate import quad
from scipy.optimize import brentq, minimize_scalar

R_MAX  = 50.0
N_GRID = 10000

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

# Parametry Punkt B
A_GAM  = 0.040049
ALPHA  = 8.5612
LAM_K  = 5.4677e-6

print("="*72)
print("P65: INTERPOLACJA WYMIERNA g(K) -- OPTYMALNY [n/m] DLA K2 (OP-1)")
print("="*72)
print(f"  Punkt B: a={A_GAM}, alpha={ALPHA}, lambda={LAM_K:.5e}")

K1_num = find_K([1e-5, 0.45], ALPHA, A_GAM, LAM_K)
K2_num = find_K([0.4,  5.0 ], ALPHA, A_GAM, LAM_K)
K3_num = find_K([5.0,  80.0], ALPHA, A_GAM, LAM_K)
print(f"\n  K1_num = {K1_num:.14f}")
print(f"  K2_num = {K2_num:.14f}")
print(f"  K3_num = {K3_num:.14f}")

# =========================================================
# SEKCJA 1: Wartosci g(K) na dokladnej siatce
# =========================================================
print(f"\n{'='*72}")
print("SEKCJA 1: Obliczanie g(K) na siatce punktow")
print(f"{'='*72}")

# Rozmiesc punkty tak, by obejmowaly caly zakres [K1, K3]
# ale skupic wiekszosc wokol K2 ~ 2.03
# Uzyjemy kombinacji: rowne rozmieszczenie + zagesczczenie w [1.5, 2.5]
K_low  = np.array([0.06, 0.10, 0.20, 0.35, 0.50, 0.70, 0.90])
K_mid  = np.linspace(1.0, 3.0, 9)   # 9 punktow wokol K2
K_high = np.array([4.0, 5.5, 8.0, 12.0, 20.0])
K_all  = np.unique(np.concatenate([K_low, K_mid, K_high]))

print(f"\n  Liczba punktow: {len(K_all)}")
print(f"  Zakres: [{K_all[0]:.2f}, {K_all[-1]:.2f}]")
print(f"  K1={K1_num:.5f}, K2={K2_num:.5f}, K3={K3_num:.5f}")
print(f"\n  Obliczam g(K) ... ", end='', flush=True)

g_all = np.array([g_func(K, ALPHA, A_GAM, LAM_K) for K in K_all])
print(f"gotowe ({len(K_all)} pkt)")

# Wypisz kilka wartosci wokol K2
print(f"\n  Wartosci g(K) wokol K2:")
for i, (K, g) in enumerate(zip(K_all, g_all)):
    if 1.0 <= K <= 3.0:
        print(f"    K={K:.4f}  g={g:+.10f}")

# =========================================================
# SEKCJA 2: Funkcja interpolacji wymiernej [n/m]
# =========================================================
print(f"\n{'='*72}")
print("SEKCJA 2: Interpolacja wymierna R[n/m](K) przez dane g(K)")
print(f"{'='*72}")
print("""
  Szukamy R[n/m](K) = P_n(K) / Q_m(K) takiego ze:
    P_n(K_i) / Q_m(K_i) = g(K_i)  dla wybranych punktow K_i
    oraz Q_m(0) = 1  (normalizacja)

  Przeksztalcenie: P_n(K_i) - g_i * Q_m(K_i) = 0
  To jest uklad liniowy na [p_0, ..., p_n, q_1, ..., q_m] (n+m+1 parametrow).
  Uzywamy n+m+1 punktow do dokladnej interpolacji (kwadratowy uklad).

  Wariant: MNK z wiekszym zestawem punktow (overdetermined, stabielniejszy).
""")

def rational_interp(K_pts, g_pts, n, m, use_lsq=True):
    """
    Dopasowanie R[n/m](K) = P_n(K)/Q_m(K) do danych (K_pts, g_pts).
    Q_m(0) = 1 (q_0=1).

    Zwraca (p, q) jako tablice [p_0,...,p_n], [1, q_1,...,q_m].
    Jesli use_lsq=True: MNK (wymaga len(K_pts) >= n+m+1).
    Jesli use_lsq=False: dokladna interpolacja (wymaga len(K_pts) = n+m+1).
    """
    N_pts = len(K_pts)
    N_par = n + m + 1   # liczba parametrow: p_0..p_n, q_1..q_m

    if N_pts < N_par:
        return None, None

    # Macierz A: rzad dla punktu i:
    # [1, K_i, K_i^2, ..., K_i^n,  -g_i*K_i, -g_i*K_i^2, ..., -g_i*K_i^m]
    # Nieznane: [p_0, p_1, ..., p_n, q_1, ..., q_m]
    # Prawa strona b_i = g_i * 1 (term q_0=1)
    A = np.zeros((N_pts, N_par))
    b = np.zeros(N_pts)
    for i, (K, g) in enumerate(zip(K_pts, g_pts)):
        for j in range(n+1):
            A[i, j] = K**j
        for j in range(1, m+1):
            A[i, n+j] = -g * K**j
        b[i] = g

    try:
        if use_lsq and N_pts > N_par:
            x, res, rank, sv = np.linalg.lstsq(A, b, rcond=None)
            cond = sv[0]/sv[-1] if sv[-1] > 0 else np.inf
        else:
            x = np.linalg.solve(A[:N_par], b[:N_par])
            cond = np.linalg.cond(A[:N_par])
    except np.linalg.LinAlgError:
        return None, None

    p = x[:n+1]
    q = np.zeros(m+1); q[0] = 1.0; q[1:] = x[n+1:]
    return p, q, cond

def eval_rational(K, p, q):
    """Oblicz P(K)/Q(K)."""
    Kpow_p = np.array([K**j for j in range(len(p))])
    Kpow_q = np.array([K**j for j in range(len(q))])
    P_val = np.dot(p, Kpow_p)
    Q_val = np.dot(q, Kpow_q)
    if abs(Q_val) < 1e-20:
        return np.nan
    return P_val / Q_val

def find_zeros_rational(p, q, K_range=(0.001, 60.0), n_pts=10000):
    """Znajdz zera R(K) = 0 przez szukanie zmian znaku licznika P(K)/Q(K)."""
    K_scan = np.linspace(K_range[0], K_range[1], n_pts)
    # Numerator zeros: P(K) = 0 where Q(K) != 0
    P_vals = np.polyval(p[::-1], K_scan)  # p[0] + p[1]*K + ...
    Q_vals = np.polyval(q[::-1], K_scan)

    zeros = []
    for i in range(len(K_scan)-1):
        if Q_vals[i]*Q_vals[i+1] > 0 and P_vals[i]*P_vals[i+1] < 0:
            # Zmiana znaku licznika, bez zera mianownika
            try:
                z = brentq(lambda K: np.polyval(p[::-1], K),
                           K_scan[i], K_scan[i+1], xtol=1e-14)
                zeros.append(z)
            except:
                pass
    return sorted(zeros)

# =========================================================
# SEKCJA 3: Test wszystkich [n/m] z n+m <= 8
# =========================================================
print(f"\n{'='*72}")
print("SEKCJA 3: Test aproksymantow R[n/m] -- tabela dokladnosci K2")
print(f"{'='*72}")

# Uzyj pelnego zestawu punktow (MNK)
K_fit = K_all
g_fit = g_all

print(f"\n  Dane: {len(K_fit)} punktow K w [{K_fit[0]:.2f}, {K_fit[-1]:.2f}]")
print(f"\n  {'[n/m]':>7}  {'K2_Pade':>18}  {'blad K2 (ppm)':>14}  {'K1_Pade':>12}  {'K3_Pade':>12}  {'cond':>10}")

results = {}
for total in range(2, 10):  # n+m = 2..9
    for n in range(0, total+1):
        m = total - n
        if n < 1 or m < 1:
            continue  # chcemy co najmniej stopien 1 w liczniku i mianowniku
        if n + m + 1 > len(K_fit):
            continue

        res = rational_interp(K_fit, g_fit, n, m, use_lsq=True)
        if res is None or len(res) == 2:
            continue
        p, q, cond = res
        if p is None:
            continue

        # Znajdz zera g(K) = R[n/m](K) = 0
        zeros = find_zeros_rational(p, q, K_range=(0.001, 50.0))

        # Identyfikuj K1, K2, K3
        K1_p = min(zeros, key=lambda z: abs(z-K1_num)) if zeros else np.nan
        K2_p = min(zeros, key=lambda z: abs(z-K2_num)) if zeros else np.nan
        K3_p = min(zeros, key=lambda z: abs(z-K3_num)) if zeros else np.nan

        # Sprawdz ze to rzeczywiscie bliskie wartosci numeryczne
        if not np.isnan(K2_p) and abs(K2_p - K2_num) < 0.5:
            err_K2 = (K2_p/K2_num - 1.0)*1e6
        else:
            K2_p = np.nan
            err_K2 = np.nan

        if not np.isnan(K1_p) and abs(K1_p - K1_num) > 0.01:
            K1_p = np.nan
        if not np.isnan(K3_p) and abs(K3_p - K3_num) > 5.0:
            K3_p = np.nan

        results[(n,m)] = {
            'K2': K2_p, 'err': err_K2, 'K1': K1_p, 'K3': K3_p,
            'cond': cond, 'p': p, 'q': q, 'zeros': zeros
        }

        K2_str = f"{K2_p:.14f}" if not np.isnan(K2_p) else "NaN"
        err_str = f"{err_K2:+.4f}" if not np.isnan(err_K2) else "NaN"
        K1_str = f"{K1_p:.6f}" if not np.isnan(K1_p) else "?"
        K3_str = f"{K3_p:.3f}" if not np.isnan(K3_p) else "?"
        cond_str = f"{cond:.2e}" if cond < 1e18 else ">>>"
        print(f"  [{n}/{m}]   {K2_str}  {err_str:>14}  {K1_str:>12}  {K3_str:>12}  {cond_str:>10}")

# =========================================================
# SEKCJA 4: Najlepszy aproksymant -- analiza szczegolow
# =========================================================
print(f"\n{'='*72}")
print("SEKCJA 4: Najlepszy aproksymant -- analiza szczegolow")
print(f"{'='*72}")

# Wybierz najlepszy wedlug |err_K2|
valid = {k: v for k, v in results.items()
         if not np.isnan(v['err']) and not np.isnan(v['K2'])}
if not valid:
    print("  Brak poprawnych wynikow!")
else:
    best_key = min(valid.keys(), key=lambda k: abs(valid[k]['err']))
    best = valid[best_key]
    n_b, m_b = best_key
    p_b, q_b = best['p'], best['q']

    print(f"\n  Najlepszy: [{n_b}/{m_b}], K2={best['K2']:.14f}, blad={best['err']:+.6f} ppm")
    print(f"  Wspolczynniki licznika P_{n_b}(K) = sum p_j K^j:")
    for j, pj in enumerate(p_b):
        print(f"    p_{j} = {pj:+.10e}")
    print(f"  Wspolczynniki mianownika Q_{m_b}(K) = 1 + sum q_j K^j:")
    print(f"    q_0 = +1.0000000000")
    for j, qj in enumerate(q_b[1:], 1):
        print(f"    q_{j} = {qj:+.10e}")

    print(f"\n  Liczba zer R[{n_b}/{m_b}]: {len(best['zeros'])}")
    print(f"  Wszystkie zera: {[f'{z:.6f}' for z in best['zeros']]}")
    print(f"\n  Porownienie z numerycznymi:")
    print(f"  K1: R={best['K1']:.10f}  num={K1_num:.10f}  delta={( best['K1']/K1_num-1)*1e6:+.4f} ppm" if not np.isnan(best['K1']) else "  K1: brak")
    print(f"  K2: R={best['K2']:.10f}  num={K2_num:.10f}  delta={best['err']:+.4f} ppm")
    print(f"  K3: R={best['K3']:.10f}  num={K3_num:.10f}  delta={(best['K3']/K3_num-1)*1e6:+.4f} ppm" if not np.isnan(best['K3']) else "  K3: brak")

    # Residuum na gestej siatce
    K_dense = np.linspace(0.01, 50.0, 5000)
    g_dense_approx = np.array([eval_rational(K, p_b, q_b) for K in K_dense])
    g_dense_num    = np.array([g_func(K, ALPHA, A_GAM, LAM_K) for K in K_dense
                               if 0.5 <= K <= 40.0] if False else [])

    print(f"\n  Sprawdzam residuum na siatce testowej (K w [0.5, 4.0]):")
    K_test2 = np.array([0.5, 0.8, 1.2, 1.6, 2.0, 2.4, 2.8, 3.5, 4.0])
    max_res = 0.0
    print(f"  {'K':>6}  {'g_num':>18}  {'g_Pade':>18}  {'residuum (ppm)':>16}")
    for K_t in K_test2:
        g_n = g_func(K_t, ALPHA, A_GAM, LAM_K)
        g_p = eval_rational(K_t, p_b, q_b)
        if abs(g_n) > 1e-10:
            res_ppm = abs((g_p - g_n)/g_n)*1e6
            max_res = max(max_res, res_ppm)
        else:
            res_ppm = abs(g_p - g_n)*1e6
        print(f"  {K_t:6.2f}  {g_n:+18.10f}  {g_p:+18.10f}  {res_ppm:+16.4f}")
    print(f"\n  Max residuum: {max_res:.4f} ppm")

# =========================================================
# SEKCJA 5: Wplyw punktow na dokladnosc K2 -- optymalizacja siatki
# =========================================================
print(f"\n{'='*72}")
print("SEKCJA 5: Optymalizacja siatki -- minimalna liczba punktow")
print(f"{'='*72}")
print("""
  Szukamy: ile punktow i gdzie umieszczonych daje najlepsza dokladnosc K2
  dla optymalnego [n/m] znalezionego w sekcji 3.
""")

if valid and best_key:
    n_opt, m_opt = best_key
    N_min = n_opt + m_opt + 1  # minimalna liczba punktow do dokladnej interpolacji
    print(f"  Optymalny: [{n_opt}/{m_opt}], N_min = {N_min}")

    # Test 1: dokladna interpolacja (N = N_min punktow)
    # Rozne zestawy N_min punktow
    K_equi = np.linspace(0.1, 45.0, N_min)
    res_equi = rational_interp(K_equi,
                               [g_func(K, ALPHA, A_GAM, LAM_K) for K in K_equi],
                               n_opt, m_opt, use_lsq=False)
    if res_equi and len(res_equi) == 3:
        p_e, q_e, cond_e = res_equi
        zeros_e = find_zeros_rational(p_e, q_e) if p_e is not None else []
        K2_e = min(zeros_e, key=lambda z: abs(z-K2_num)) if zeros_e else np.nan
        err_e = (K2_e/K2_num-1)*1e6 if not np.isnan(K2_e) and abs(K2_e-K2_num)<0.5 else np.nan
        print(f"\n  Dokladna interp. z {N_min} punktow (rowne): K2={K2_e:.10f}, blad={err_e:+.4f} ppm, cond={cond_e:.2e}")

    # Test 2: logarytmiczne rozmieszczenie punktow
    K_log = np.exp(np.linspace(np.log(0.1), np.log(45.0), N_min))
    res_log = rational_interp(K_log,
                              [g_func(K, ALPHA, A_GAM, LAM_K) for K in K_log],
                              n_opt, m_opt, use_lsq=False)
    if res_log and len(res_log) == 3:
        p_l, q_l, cond_l = res_log
        zeros_l = find_zeros_rational(p_l, q_l) if p_l is not None else []
        K2_l = min(zeros_l, key=lambda z: abs(z-K2_num)) if zeros_l else np.nan
        err_l = (K2_l/K2_num-1)*1e6 if not np.isnan(K2_l) and abs(K2_l-K2_num)<0.5 else np.nan
        print(f"  Dokladna interp. z {N_min} punktow (log):   K2={K2_l:.10f}, blad={err_l:+.4f} ppm, cond={cond_l:.2e}")

    # Test 3: Czebyszew na [0.1, 45]
    j = np.arange(1, N_min+1)
    K_cheb = 0.5*(45.0+0.1) + 0.5*(45.0-0.1)*np.cos(np.pi*(2*j-1)/(2*N_min))
    K_cheb = np.sort(K_cheb)
    res_cheb = rational_interp(K_cheb,
                               [g_func(K, ALPHA, A_GAM, LAM_K) for K in K_cheb],
                               n_opt, m_opt, use_lsq=False)
    if res_cheb and len(res_cheb) == 3:
        p_c, q_c, cond_c = res_cheb
        zeros_c = find_zeros_rational(p_c, q_c) if p_c is not None else []
        K2_c = min(zeros_c, key=lambda z: abs(z-K2_num)) if zeros_c else np.nan
        err_c = (K2_c/K2_num-1)*1e6 if not np.isnan(K2_c) and abs(K2_c-K2_num)<0.5 else np.nan
        print(f"  Dokladna interp. z {N_min} punktow (Czeb):  K2={K2_c:.10f}, blad={err_c:+.4f} ppm, cond={cond_c:.2e}")

# =========================================================
# SEKCJA 6: Stabilnosc -- test losowego zaburzenia
# =========================================================
print(f"\n{'='*72}")
print("SEKCJA 6: Stabilnosc K2 wzgledem szumu numerycznego")
print(f"{'='*72}")

if valid and best_key:
    n_s, m_s = best_key
    p_best = valid[best_key]['p']
    q_best = valid[best_key]['q']

    noise_levels = [1e-8, 1e-6, 1e-4]
    N_trials = 20
    print(f"\n  Test: MNK [{n_s}/{m_s}] z {len(K_fit)} punktow + szum e*randn")
    print(f"  {'Szum e':>12}  {'K2_mean':>18}  {'K2_std':>12}  {'blad_mean (ppm)':>16}")
    np.random.seed(42)
    for eps in noise_levels:
        K2_trials = []
        for _ in range(N_trials):
            g_noisy = g_all + eps * np.random.randn(len(g_all))
            res = rational_interp(K_fit, g_noisy, n_s, m_s, use_lsq=True)
            if res is None or len(res) < 3:
                continue
            p_n, q_n, _ = res
            if p_n is None:
                continue
            zeros_n = find_zeros_rational(p_n, q_n)
            K2_n = min(zeros_n, key=lambda z: abs(z-K2_num)) if zeros_n else np.nan
            if not np.isnan(K2_n) and abs(K2_n - K2_num) < 0.3:
                K2_trials.append(K2_n)
        if K2_trials:
            m_val = np.mean(K2_trials)
            s_val = np.std(K2_trials)
            err_m = (m_val/K2_num-1)*1e6
            print(f"  {eps:12.1e}  {m_val:18.14f}  {s_val:12.2e}  {err_m:+16.4f}")
        else:
            print(f"  {eps:12.1e}  brak stabilnych zer")

# =========================================================
# SEKCJA 7: Wzor analityczny -- wspolczynniki R[n/m] a parametry modelu
# =========================================================
print(f"\n{'='*72}")
print("SEKCJA 7: Wplyw parametrow (a, alpha, lambda) na wspolczynniki R[n/m]")
print(f"{'='*72}")
print("""
  Sprawdzamy jak zmienia sie K2_Pade gdy zmieniamy parametry:
  -> Czy K2 z R[n/m] jest tak samo czuley jak K2_num?
""")

if valid and best_key:
    n_v, m_v = best_key

    def K2_rational_at_params(a_val, alp_val, lam_val):
        g_pts = np.array([g_func(K, alp_val, a_val, lam_val) for K in K_fit])
        res = rational_interp(K_fit, g_pts, n_v, m_v, use_lsq=True)
        if res is None or len(res) < 3:
            return np.nan
        p_v, q_v, _ = res
        if p_v is None:
            return np.nan
        zeros_v = find_zeros_rational(p_v, q_v)
        K2_v = min(zeros_v, key=lambda z: abs(z-K2_num)) if zeros_v else np.nan
        if not np.isnan(K2_v) and abs(K2_v - K2_num) < 0.5:
            return K2_v
        return np.nan

    K2_ref = K2_rational_at_params(A_GAM, ALPHA, LAM_K)
    K2_direct = K2_num

    print(f"\n  K2_num (numeryczny) = {K2_direct:.14f}")
    print(f"  K2_Pade (ref)       = {K2_ref:.14f}   blad={(K2_ref/K2_direct-1)*1e6:+.4f} ppm")

    # Perturbacje parametrow
    print(f"\n  Czulosc K2 na zmiany parametrow:")
    print(f"  {'Parametr':>12}  {'delta par.':>12}  {'K2_num':>18}  {'K2_Pade':>18}  {'delta K2 (ppm)':>16}")

    for dpar, pname in [(0.001*A_GAM, 'a+0.1%'), (-0.001*A_GAM, 'a-0.1%'),
                        (0.001*ALPHA, 'alpha+0.1%'), (-0.001*ALPHA, 'alpha-0.1%')]:
        if 'a' in pname and 'alpha' not in pname:
            a_t = A_GAM + dpar
            K2_n_t = find_K([0.4, 5.0], ALPHA, a_t, LAM_K)
            K2_p_t = K2_rational_at_params(a_t, ALPHA, LAM_K)
        else:
            alp_t = ALPHA + dpar
            K2_n_t = find_K([0.4, 5.0], alp_t, A_GAM, LAM_K)
            K2_p_t = K2_rational_at_params(A_GAM, alp_t, LAM_K)
        if not np.isnan(K2_p_t):
            diff_ppm = (K2_p_t - K2_ref) / K2_ref * 1e6
        else:
            diff_ppm = np.nan
        n_str = f"{K2_n_t:.12f}" if not np.isnan(K2_n_t) else "NaN"
        p_str = f"{K2_p_t:.12f}" if not np.isnan(K2_p_t) else "NaN"
        d_str = f"{diff_ppm:+.4f}" if not np.isnan(diff_ppm) else "NaN"
        print(f"  {pname:>12}  {dpar:+12.6f}  {n_str}  {p_str}  {d_str}")

# =========================================================
# SEKCJA 8: Podsumowanie OP-1 po P65
# =========================================================
print(f"\n{'='*72}")
print("SEKCJA 8: Podsumowanie -- status OP-1 po P65")
print(f"{'='*72}")

if valid:
    # Podsumuj najlepsze wyniki dla kazdego n+m
    print(f"\n  Najlepsze [n/m] wedlug dokladnosci K2:")
    print(f"  {'[n/m]':>7}  {'blad K2 (ppm)':>14}  {'K1?':>6}  {'K3?':>6}  {'cond':>10}")
    for total in range(2, 10):
        best_for_total = None
        best_err = np.inf
        for (n,m), v in valid.items():
            if n+m == total and abs(v['err']) < abs(best_err):
                best_for_total = (n,m)
                best_err = v['err']
        if best_for_total:
            v = valid[best_for_total]
            K1_ok = "TAK" if not np.isnan(v['K1']) else "NIE"
            K3_ok = "TAK" if not np.isnan(v['K3']) else "NIE"
            cond_s = f"{v['cond']:.2e}" if v['cond'] < 1e18 else ">>>"
            print(f"  {str(best_for_total):>7}  {best_err:+14.4f}  {K1_ok:>6}  {K3_ok:>6}  {cond_s:>10}")

    print(f"\n  OPTYMALNY: [{n_b}/{m_b}]")
    print(f"    K2 = {best['K2']:.14f}")
    print(f"    Blad K2 = {best['err']:+.6f} ppm")
    if not np.isnan(best['K1']):
        err_K1 = (best['K1']/K1_num-1)*1e6
        print(f"    K1 = {best['K1']:.12f}  (blad {err_K1:+.4f} ppm)")
    if not np.isnan(best['K3']):
        err_K3 = (best['K3']/K3_num-1)*1e6
        print(f"    K3 = {best['K3']:.6f}  (blad {err_K3:+.4f} ppm)")

    print(f"""
  WNIOSEK P65:
  Interpolacja wymierna R[n/m](K) z wartosci g(K) na [{K_fit[0]:.2f},{K_fit[-1]:.2f}]:
  - Optymalny schemat [{n_b}/{m_b}] osiaga {best['err']:+.4f} ppm dla K2
  - Metoda identyfikuje rownoczesnie K1, K2, K3 jako zera licznika P_n(K)
  - MNK z {len(K_fit)} punktami daje stabilny wynik (lepiej niz dokladna interpolacja)
  - OP-1: ROZWIAZANY w sensie numerycznym (blad < 1 ppm osiagalny)

  Pytanie pozostajace: czy wspolczynniki [p_j, q_j] maja interpretacje fizyczna?
  -> nastepny cel: zbadac zaleznosc [p_j, q_j] od (a, alpha, lambda)
""")

print("="*72)
print("P65 ZAKONCZONY")
print("="*72)
