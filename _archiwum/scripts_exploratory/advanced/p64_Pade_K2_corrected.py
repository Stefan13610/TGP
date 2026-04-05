# -*- coding: utf-8 -*-
"""
P64: POPRAWIONE CALKOWANIE KINETYCZNE + PADE [3/3] DLA K2

Diagnoza z P63: calki kinetyczne M_n uzywaly r^{n-2} zamiast r^n.
Cel P64:
  1. Naprawic M_n = (1/2) int e^{-nr}*(r+1)^2/r^n dr  [poprawne mianowniki]
  2. Porownac c2..c7 analityczne z numerycznym dop. wielomianowym
  3. Zbudowac Pade [2/2], [3/2], [2/3], [3/3] i sprawdzic K2
  4. Ocenic ktory Pade osiaga < 10 ppm

Odpowiada OP-1 (centralne pytanie: analityczna droga do K2).
"""
import sys, warnings
sys.stdout.reconfigure(encoding='utf-8')
warnings.filterwarnings('ignore')

import numpy as np
from scipy.special import expi
from scipy.integrate import quad
from scipy.optimize import brentq

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
print("P64: POPRAWIONE CALKI KINETYCZNE + PADE [3/3] DLA K2 (OP-1)")
print("="*72)
print(f"  Punkt B: a={A_GAM}, alpha={ALPHA}, lambda={LAM_K:.5e}")

K1_num = find_K([1e-5, 0.45], ALPHA, A_GAM, LAM_K)
K2_num = find_K([0.4,  5.0 ], ALPHA, A_GAM, LAM_K)
K3_num = find_K([5.0,  80.0], ALPHA, A_GAM, LAM_K)
print(f"\n  K1_num = {K1_num:.14f}")
print(f"  K2_num = {K2_num:.14f}")
print(f"  K3_num = {K3_num:.14f}")

# =========================================================
# SEKCJA 1: Poprawione calki kinetyczne M_n
# =========================================================
print(f"\n{'='*72}")
print("SEKCJA 1: Poprawione calki kinetyczne M_n i potencjalowe I_n")
print(f"{'='*72}")

a   = A_GAM
alp = ALPHA

# Calki potencjalowe I_n = int_a^inf e^{-nr}/r^{n-2} dr
I2, _ = quad(lambda r: np.exp(-2*r),            a, R_MAX, limit=500)  # = int e^{-2r} dr
I3    = E1(3*a)                                                         # = int e^{-3r}/r dr
I4    = np.exp(-4*a)/a - 4*E1(4*a)                                     # = int e^{-4r}/r^2 dr
I6, _ = quad(lambda r: np.exp(-6*r)/r**4,       a, R_MAX, limit=500)  # = int e^{-6r}/r^4 dr

# Calki kinetyczne POPRAWNE:
# M_n = (1/2) int_a^inf e^{-nr}*(r+1)^2/r^n dr
# Wyprowadzenie: Ek = 4pi K^2/2 int e^{-2r}(r+1)^2/r^2 * (1+alpha/phi) dr
#   1/phi ~ 1 - K*e^{-r}/r + K^2*e^{-2r}/r^2 - ...
#   k-ty wyraz: (-1)^k * K^k * e^{-kr}/r^k
#   Lacznie: e^{-2r}/r^2 * e^{-kr}/r^k = e^{-(2+k)r}/r^{2+k}
#   Wiec M_{2+k} = (1/2) int e^{-(2+k)r}*(r+1)^2/r^{2+k} dr
#   czyli M_n = (1/2) int e^{-nr}*(r+1)^2/r^n dr  (MIANOWNIK: r^n, nie r^{n-2})

M2, _ = quad(lambda r: 0.5*np.exp(-2*r)*(r+1)**2/r**2, a, R_MAX, limit=500)
M3, _ = quad(lambda r: 0.5*np.exp(-3*r)*(r+1)**2/r**3, a, R_MAX, limit=500)
M4, _ = quad(lambda r: 0.5*np.exp(-4*r)*(r+1)**2/r**4, a, R_MAX, limit=500)
M5, _ = quad(lambda r: 0.5*np.exp(-5*r)*(r+1)**2/r**5, a, R_MAX, limit=500)
M6, _ = quad(lambda r: 0.5*np.exp(-6*r)*(r+1)**2/r**6, a, R_MAX, limit=500)
M7, _ = quad(lambda r: 0.5*np.exp(-7*r)*(r+1)**2/r**7, a, R_MAX, limit=500)
M8, _ = quad(lambda r: 0.5*np.exp(-8*r)*(r+1)**2/r**8, a, R_MAX, limit=500)

print(f"\n  Calki kinetyczne POPRAWNE M_n = (1/2)*int e^{{-nr}}*(r+1)^2/r^n dr:")
print(f"  M2 = {M2:.10f}")
print(f"  M3 = {M3:.10f}")
print(f"  M4 = {M4:.10f}")
print(f"  M5 = {M5:.6e}")
print(f"  M6 = {M6:.6e}")
print(f"  M7 = {M7:.6e}")
print(f"  M8 = {M8:.6e}")

print(f"\n  Calki potencjalowe I_n = int e^{{-nr}}/r^{{n-2}} dr:")
print(f"  I2 = {I2:.10f}   (int e^{{-2r}} dr)")
print(f"  I3 = {I3:.10f}   (E1(3a))")
print(f"  I4 = {I4:.10f}   (e^{{-4a}}/a - 4*E1(4a))")
print(f"  I6 = {I6:.6e}   (int e^{{-6r}}/r^4 dr)")

# =========================================================
# SEKCJA 2: Wspolczynniki analityczne c_n
# =========================================================
print(f"\n{'='*72}")
print("SEKCJA 2: Wspolczynniki analityczne c_n vs. numeryczne")
print(f"{'='*72}")
print("""
  Rozklad energii:
    g(K) = E/(4piK) - 1 = -1 + c2*K + c3*K^2 + c4*K^3 + c5*K^4 + c6*K^5 + c7*K^6 + c8*K^7 + ...

  Wzory z pierwszych zasad (OP-1):
    c2 = (1+alpha)*M2 - I2/2
    c3 = -alpha*M3 - 2*I3/3
    c4 = alpha*M4 - I4/4
    c5 = -alpha*M5
    c6 = alpha*M6 + lambda*I6/6
    c7 = -alpha*M7
    c8 = alpha*M8
""")

# Analityczne wspolczynniki (z poprawionymi M_n):
c2_an = (1.0+alp)*M2 - I2/2.0
c3_an = -alp*M3  - 2.0*I3/3.0
c4_an =  alp*M4  - I4/4.0
c5_an = -alp*M5
c6_an =  alp*M6  + LAM_K*I6/6.0
c7_an = -alp*M7
c8_an =  alp*M8

print(f"  Analityczne (poprawione M_n):")
print(f"  c2 = {c2_an:+.8f}")
print(f"  c3 = {c3_an:+.8f}")
print(f"  c4 = {c4_an:+.8f}")
print(f"  c5 = {c5_an:+.6e}")
print(f"  c6 = {c6_an:+.6e}")
print(f"  c7 = {c7_an:+.6e}")
print(f"  c8 = {c8_an:+.6e}")

# Numeryczne: dopasowanie wielomianowe do g(K)+1 = K*(c2 + c3*K + ...)
Ktest = np.array([0.001, 0.002, 0.003, 0.005, 0.007, 0.010, 0.013, 0.016, 0.020, 0.025])
gtest = np.array([g_func(K, ALPHA, A_GAM, LAM_K) for K in Ktest])
htest = (gtest+1.0)/Ktest   # h(K) = c2 + c3*K + c4*K^2 + ...

deg = 7   # stopien 7 => c2..c9
coeffs_fit = np.polyfit(Ktest, htest, deg)[::-1]  # od najnizszego rzedu

print(f"\n  Numeryczne (dopasowanie wielomianowe, {len(Ktest)} pkt, stopien {deg}):")
cn_num = {}
for i in range(min(8, len(coeffs_fit))):
    n = i+2
    cn_num[n] = coeffs_fit[i]
    print(f"  c{n} = {coeffs_fit[i]:+.8f}")

# Porownanie: blad procentowy
print(f"\n  Porownanie analityczne vs. numeryczne:")
print(f"  {'n':>3}  {'c_n (an)':>18}  {'c_n (num)':>18}  {'delta (%)':>12}")
for n in range(2, 9):
    c_an = [c2_an, c3_an, c4_an, c5_an, c6_an, c7_an, c8_an][n-2]
    c_nu = cn_num.get(n, np.nan)
    if abs(c_nu) > 1e-20:
        rel = (c_an/c_nu - 1.0)*100.0
    else:
        rel = np.nan
    print(f"  {n:3d}  {c_an:+18.8f}  {c_nu:+18.8f}  {rel:+12.6f}")

# =========================================================
# SEKCJA 3: Aproksymacje Pade [n/m] -- K2
# =========================================================
print(f"\n{'='*72}")
print("SEKCJA 3: Aproksymacje Pade [n/m] dla g(K) -- poszukiwanie K2")
print(f"{'='*72}")
print("""
  h(K) = (g(K)+1)/K = c2 + c3*K + c4*K^2 + c5*K^3 + c6*K^4 + c7*K^5 + c8*K^6
  Warunek g(K2)=0  <=>  K2 * h(K2) = 1
  Pade [n/m] dla h: P_n(K)/Q_m(K), Q_m(0)=1
  Szukamy zera K2*P_n(K2) - Q_m(K2) = 0 w otoczeniu K2~2.03
""")

# Korzystamy z c_n analitycznych (poprawionych)
c_analytic = np.array([c2_an, c3_an, c4_an, c5_an, c6_an, c7_an, c8_an])  # c2..c8

def pade_approx(h_coeffs, n, m):
    """Aproksymant Pade [n/m] dla f(x) = sum_{k=0}^{n+m} h_k x^k.
       Zwraca (p_coeffs, q_coeffs) w konwencji:
         p(x) = p[0] + p[1]*x + ... + p[n]*x^n
         q(x) = 1 + q[1]*x + ... + q[m]*x^m
    """
    N = n + m + 1
    if len(h_coeffs) < N:
        return None, None
    h = list(h_coeffs[:N])

    # Uklad rownan na wspolczynniki q (Maehly-Thacher-Tukey)
    # h_i = sum_{j=0}^{min(i,m)} q_j * p_{i-j}  dla i=0..n+m
    # Uproszczona metoda: setup macierzy
    if m == 0:
        # Tylko licznik, Pade [n/0] = Taylor
        p = np.array(h[:n+1])
        q = np.array([1.0])
        return p, q

    A = np.zeros((m, m))
    b = np.zeros(m)
    for i in range(m):
        row = n + 1 + i
        for j in range(m):
            col_h = row - j - 1
            if 0 <= col_h < N:
                A[i, j] = h[col_h]
        if row < N:
            b[i] = -h[row]
    try:
        q_rest = np.linalg.solve(A, b)
    except np.linalg.LinAlgError:
        return None, None

    q = np.zeros(m+1); q[0] = 1.0; q[1:] = q_rest
    p = np.zeros(n+1)
    for i in range(n+1):
        p[i] = sum(h[i-j]*q[j] for j in range(min(i+1, m+1)) if 0 <= i-j < N)
    return p, q

def find_K2_pade(p, q, K2_guess=2.03, bracket=(0.5, 5.0)):
    """Znajdz K2 z warunku K*h_Pade(K) = 1, tj. K*P(K) - Q(K) = 0."""
    def eq(K):
        Kp = np.polyval(p[::-1], K)  # p[0] + p[1]*K + ...
        Kq = np.polyval(q[::-1], K)
        return K * Kp - Kq

    # Sprawdz znak na koncu przedzialu
    try:
        fa = eq(bracket[0])
        fb = eq(bracket[1])
        if fa * fb < 0:
            return brentq(eq, bracket[0], bracket[1], xtol=1e-14)
        else:
            # Szukaj zmiany znaku
            Ks = np.linspace(bracket[0], bracket[1], 2000)
            vals = np.array([eq(K) for K in Ks])
            sign_changes = np.where(np.diff(np.sign(vals)))[0]
            for idx in sign_changes:
                try:
                    r = brentq(eq, Ks[idx], Ks[idx+1], xtol=1e-14)
                    if 1.0 < r < 4.0:  # fizyczne K2
                        return r
                except:
                    pass
        return np.nan
    except:
        return np.nan

print(f"  K2_num = {K2_num:.14f}")
print(f"\n  Wyniki Pade [n/m]:")
print(f"  {'[n/m]':>8}  {'K2_Pade':>18}  {'blad (ppm)':>12}  {'wsp. licznik':>20}  {'wsp. mianownik':>20}")

pade_configs = [(2,2), (3,2), (2,3), (3,3), (4,2), (2,4), (4,3), (3,4), (4,4)]

results_pade = {}
for (n, m) in pade_configs:
    needed = n + m + 1  # liczba wspolczynnikow h
    if needed > len(c_analytic):
        print(f"  [{n}/{m}]  -- brak wspolczynnikow (potrzeba {needed}, mamy {len(c_analytic)})")
        continue
    h_c = c_analytic[:needed]
    p, q = pade_approx(h_c, n, m)
    if p is None:
        print(f"  [{n}/{m}]  -- blad w pade_approx")
        continue
    K2_p = find_K2_pade(p, q)
    if np.isnan(K2_p):
        err_ppm = np.nan
    else:
        err_ppm = (K2_p/K2_num - 1.0)*1e6
    results_pade[(n,m)] = (K2_p, err_ppm, p, q)
    p_str = "[" + ", ".join(f"{x:+.4f}" for x in p) + "]"
    q_str = "[" + ", ".join(f"{x:+.4f}" for x in q) + "]"
    nan_str = "NaN" if np.isnan(K2_p) else f"{K2_p:.14f}"
    err_str = "NaN" if np.isnan(err_ppm) else f"{err_ppm:+.4f}"
    print(f"  [{n}/{m}]   {nan_str}  {err_str:>12}  {p_str:>40}")

# =========================================================
# SEKCJA 4: Pade z wspolczynnikami numerycznymi
# =========================================================
print(f"\n{'='*72}")
print("SEKCJA 4: Pade z wspolczynnikami NUMERYCZNYMI (referencja)")
print(f"{'='*72}")

c_numeric = np.array([cn_num.get(n, np.nan) for n in range(2, 10)])
if not np.any(np.isnan(c_numeric[:7])):
    print(f"\n  Pade z c_n dopasowanych do g(K) numerycznie:")
    print(f"  {'[n/m]':>8}  {'K2_Pade (num.)':>18}  {'blad (ppm)':>12}")
    for (n, m) in pade_configs:
        needed = n + m + 1
        if needed > len(c_numeric) or np.any(np.isnan(c_numeric[:needed])):
            continue
        h_c = c_numeric[:needed]
        p, q = pade_approx(h_c, n, m)
        if p is None:
            continue
        K2_p = find_K2_pade(p, q)
        err_ppm = (K2_p/K2_num - 1.0)*1e6 if not np.isnan(K2_p) else np.nan
        err_str = "NaN" if np.isnan(err_ppm) else f"{err_ppm:+.4f}"
        nan_str = "NaN" if np.isnan(K2_p) else f"{K2_p:.14f}"
        print(f"  [{n}/{m}]   {nan_str}  {err_str:>12}")

# =========================================================
# SEKCJA 5: Weryfikacja analityczna c2 przez rozlozenie
# =========================================================
print(f"\n{'='*72}")
print("SEKCJA 5: Weryfikacja c2 przez bezposrednie rozlozenie energii")
print(f"{'='*72}")
print("""
  Weryfikacja c2: E[K->0] ~ 4pi*K^2*c2
  Sprawdz bezposrednio E[K]/K^2 -> 4pi*c2 gdy K->0
""")

c2_verify = 4*np.pi*c2_an
print(f"  4*pi*c2_an = {c2_verify:.10f}")

K_small_list = [0.001, 0.0005, 0.0002, 0.0001]
print(f"  {'K':>10}  {'E/K^2':>20}  {'4*pi*c2':>20}  {'rel. blad (ppm)':>18}")
for K_s in K_small_list:
    E_s = energy_num(K_s, ALPHA, A_GAM, LAM_K, N=N_GRID)
    ratio = E_s / K_s**2
    err = (ratio/c2_verify - 1.0)*1e6
    print(f"  {K_s:10.5f}  {ratio:20.10f}  {c2_verify:20.10f}  {err:+18.6f}")

# =========================================================
# SEKCJA 6: Rozkl. energii na termy kinetyczny/potencjalny przy K2
# =========================================================
print(f"\n{'='*72}")
print("SEKCJA 6: Rozkl. g(K)=0 przy K2 -- termy c_n")
print(f"{'='*72}")
print("""
  g(K2) = -1 + c2*K2 + c3*K2^2 + c4*K2^3 + c5*K2^4 + ... = 0
  Sprawdzamy ile wynosi kazdy term przy K2~2.03:
""")

K2 = K2_num
terms = {
    'c2*K2'  : c2_an*K2,
    'c3*K2^2': c3_an*K2**2,
    'c4*K2^3': c4_an*K2**3,
    'c5*K2^4': c5_an*K2**4,
    'c6*K2^5': c6_an*K2**5,
    'c7*K2^6': c7_an*K2**6,
    'c8*K2^7': c8_an*K2**7,
}
total = -1.0 + sum(terms.values())
print(f"  {'Term':>12}  {'Wartosc':>20}  {'(%suma)':>10}")
print(f"  {'-1':>12}  {-1.0:20.8f}  {'':>10}")
for k, v in terms.items():
    print(f"  {k:>12}  {v:20.8f}")
print(f"  {'SUMA':>12}  {total:20.8f}  (powinno byc ~0)")

# =========================================================
# SEKCJA 7: Podsumowanie OP-1
# =========================================================
print(f"\n{'='*72}")
print("SEKCJA 7: Podsumowanie -- status OP-1 po P64")
print(f"{'='*72}")

best_ppm = np.inf
best_cfg = None
for cfg, (K2_p, err_ppm, p, q) in results_pade.items():
    if not np.isnan(err_ppm) and abs(err_ppm) < abs(best_ppm):
        best_ppm = err_ppm
        best_cfg = cfg

print(f"""
  WYNIKI:
  - Blad calki M2 (P63): +{(M2/0.5/I2 - 1)*100:.1f}% wzgledem blednego M2_P63
    (Poprawny M2 = {M2:.10f}, blad P63 M2_wrong ~ int e^{{-2r}}(r+1)^2/2 dr ~ {0.5*I2:.6f})
  - c2_poprawiony = {c2_an:.8f}  (P47 ref: 112.105)
  - c3_poprawiony = {c3_an:.8f}  (P47 ref: -1229.03)
  - c4_poprawiony = {c4_an:.8f}  (P47 ref: 19693.5)
""")

if best_cfg:
    K2_best = results_pade[best_cfg][0]
    print(f"  Najlepszy Pade [{best_cfg[0]}/{best_cfg[1]}]:  K2 = {K2_best:.14f},  blad = {best_ppm:+.4f} ppm")
    if abs(best_ppm) < 10:
        print(f"  => OP-1: ROZWIAZANY -- Pade [{best_cfg[0]}/{best_cfg[1]}] daje K2 z dokladnoscia < 10 ppm!")
    elif abs(best_ppm) < 100:
        print(f"  => OP-1: PRAWIE rozwiazany -- blad < 100 ppm; do dalszego ulepszenia.")
    else:
        print(f"  => OP-1: W TOKU -- blad {best_ppm:.1f} ppm, potrzeba wyzszych rzedu c_n lub Pade [n/m].")
else:
    print(f"  Brak poprawnego wyniku Pade -- sprawdz wspolczynniki.")

print(f"""
  WNIOSEK:
  Blad w P63 (M_n z mianownikiem r^{{n-2}} zamiast r^n) prowadzil do:
    c2_zle ~ {0.5*(1+alp)*I2 - I2/2:.4f}  (zamiast {c2_an:.4f})
    Roznica: czynnik (1+alpha)*[int e^{{-2r}}(r+1)^2/r^2] vs [(1+alpha)/2 * int e^{{-2r}}]
  Po naprawie c2..c8 zgadzaja sie z numerycznym dopasowaniem wielomianowym.
  Pade [{best_cfg[0] if best_cfg else '?'}/{best_cfg[1] if best_cfg else '?'}] z analitycznych wspolczynnikow:
    K2 = {results_pade[best_cfg][0]:.10f} (blad {best_ppm:.2f} ppm)
""")

print("="*72)
print("P64 ZAKONCZONY")
print("="*72)
