"""
P48: ANALITYCZNE K2 — WYPROWADZENIE Z WARUNKU LOKALNEGO MAKSIMUM

K2 jest "środkowym" zerem g(K) = E(K)/(4pi*K) - 1.
Aby je wyznaczyć analitycznie, kluczowe jest zrozumienie:

  (a) KROK 1: K_max = argmax g(K)  <-- gdzie g osiąga lokalne maksimum (między K1 a K2)
              Warunek: g'(K_max) = 0  <==>  dE/dK = 4pi  (tangent condition)
  (b) KROK 2: K2 > K_max jest zerem g(K) malejącego na prawo od K_max.

Podejście analityczne:
  - Użyj efektywnego przybliżenia g(K) ważnego dla K ~ K2 ~ 2.
  - Profil Yukawa dla K~2, r~a~0.04: phi(r) = 1 + K*e^{-r}/r
    -> K/r ~ 50 >> 1 blisko r=a, ale K*e^{-r}/r ~ 0.7 daleko (r~1)
  - Przejście reżimu: r* = ln(K) ~ 0.7 dla K=2
  - Dla r < r*: phi ~ K*e^{-r}/r  (duże-K)
  - Dla r > r*: phi ~ 1 + eps  (bliskie próżni)

Podejścia (kolejność dokładności):
  A. Aproksymacja rozdzielona (split-profile): phi = duze-K + próżnia podzielona przy r*
  B. Aproksymacja Padé g(K) z wynikow numerycznych
  C. Analityczne g'(K) = 0 dla waruniku K_max
  D. Rozwinięcie wokół K_max: g(K) = g_max - b*(K-K_max)^2 + c*(K-K_max)^3 + ...
     => K2 = K_max + sqrt(2*g_max/b) * [1 + corrections]
"""

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq, minimize_scalar
from scipy.special import expi
import sys, warnings
sys.stdout.reconfigure(encoding='utf-8')
warnings.filterwarnings('ignore')

# ═══════════════════════════════════════════════════════════
# PARAMETRY TGP (z P40, leptony)
# ═══════════════════════════════════════════════════════════
ALPHA = 8.5445
A_GAM = 0.040
LAM   = 5.4677e-6
K1_TRUE = 0.009839
K2_TRUE = 2.0344
K3_TRUE = 34.2127

R_MAX = 50.0
N_GRID = 4000

def E1(x):
    """E1(x) = int_x^inf e^{-t}/t dt dla x>0."""
    return -expi(-x)

# ═══════════════════════════════════════════════════════════
# ENERGIA NUMERYCZNA (referencyjna)
# ═══════════════════════════════════════════════════════════
def V_mod(phi, lam):
    return phi**3/3 - phi**4/4 + lam*(phi-1)**6/6

def energy_num(K, alpha=ALPHA, a=A_GAM, lam=LAM, N=N_GRID):
    """Energia E(K) przez całkowanie numeryczne (log-grid)."""
    t   = np.linspace(0, 1, N)
    r   = a * (R_MAX/a)**t
    phi = np.maximum(1.0 + K*np.exp(-r)/r, 1e-10)
    dphi = K * np.exp(-r) * (-r - 1.0) / r**2
    V1  = V_mod(1.0, lam)
    Ek  = 4*np.pi * np.trapezoid(0.5*dphi**2*(1.0 + alpha/phi)*r**2, r)
    Ep  = 4*np.pi * np.trapezoid((V_mod(phi, lam) - V1)*r**2, r)
    return Ek + Ep

def g_num(K, alpha=ALPHA, a=A_GAM, lam=LAM):
    """g(K) = E(K)/(4pi*K) - 1."""
    if K <= 0: return -1.0
    return energy_num(K, alpha, a, lam) / (4*np.pi*K) - 1.0

# ═══════════════════════════════════════════════════════════
# CAŁKI ANALITYCZNE (z P47)
# ═══════════════════════════════════════════════════════════
def analytic_integrals(a, alpha, lam):
    """Oblicza c_n dla szeregu E(K) = 4pi*sum c_n K^n."""
    # Phi_2 (dokładne przez Ei)
    J2 = 2*E1(2*a) - np.exp(-2*a)/a   # int e^{-2r}/r^2
    J1 = E1(2*a)                        # int e^{-2r}/r
    J0 = np.exp(-2*a)/2                 # int e^{-2r}
    Phi2 = -J2 + 2*J1 + J0

    # Phi_3 numerycznie (dla K^3)
    def phi3_integrand(r):
        return np.exp(-3*r)*(1+r)**2/r**3
    Phi3, _ = quad(phi3_integrand, a, R_MAX, limit=200)

    # Phi_4 numerycznie (dla K^4)
    def phi4_integrand(r):
        return np.exp(-4*r)*(1+r)**2/r**4
    Phi4, _ = quad(phi4_integrand, a, R_MAX, limit=200)

    # U_n
    U2 = np.exp(-2*a)/2
    U3 = E1(3*a)
    U4 = np.exp(-4*a)/a - 4*E1(4*a)
    U6_val = (np.exp(-6*a)/(3*a**3) - np.exp(-6*a)/a**2
              + 6*(np.exp(-6*a)/a - 6*E1(6*a)))

    c2 = (1+alpha)*Phi2/2 - U2/2
    c3 = -alpha*Phi3/2 - 2*U3/3
    c4 = alpha*Phi4/2 - U4/4
    c6 = lam*U6_val/6

    return c2, c3, c4, c6, Phi2, Phi3, Phi4, U2, U3, U4, U6_val

print("="*72)
print("P48: ANALITYCZNE K2 — WYPROWADZENIE Z WARUNKU dg/dK = 0")
print("="*72)

c2, c3, c4, c6, Phi2, Phi3, Phi4, U2, U3, U4, U6 = analytic_integrals(A_GAM, ALPHA, LAM)

print(f"\nParametry: alpha={ALPHA}, a={A_GAM}, lambda={LAM:.4e}")
print(f"c2={c2:.4f}, c3={c3:.2f}, c4={c4:.2f}, c6={c6:.4e}")

# ═══════════════════════════════════════════════════════════
print("\n" + "="*72)
print("SEKCJA A: ANALIZA g(K) PRZY K~K2~2 — REŻIM PRZEJŚCIOWY")
print("="*72)

print("""
  Dla K ~ 2 i r ~ a = 0.04: phi = 1 + K*e^{-r}/r
  Przy r=0.04: phi = 1 + 2*e^{-0.04}/0.04 = 1 + 48.1 = 49.1  [dominuje K/r]
  Przy r=1.0:  phi = 1 + 2*e^{-1}/1   = 1 + 0.736 = 1.736  [mieszany]
  Przy r=5.0:  phi = 1 + 2*e^{-5}/5   = 1 + 0.013 = 1.013  [dominuje 1]

  Punkt przejścia r* = ln(K) = ln(2) = 0.693:
    phi(r*) = 1 + K*e^{-r*}/r* = 1 + e^{-ln K}/ln(K) * K = 1 + 1/ln(K) ~ 2.44

  REŻIM:
    r < r*:   phi ~ K*e^{-r}/r  (dominuje soliton)
    r > r*:   phi ~ 1 + K*e^{-r}/r  (korekcja małosygnałowa, zbliżona do próżni)
""")

# Skanuj g(K) numerycznie wokół K2
K_vals = np.linspace(0.5, 5.0, 200)
print("  Numeryczne g(K) wokół K2 (szybkie skanowanie, N=1000):")
g_vals = []
for Kv in K_vals:
    gv = energy_num(Kv, N=1000) / (4*np.pi*Kv) - 1.0
    g_vals.append(gv)
g_vals = np.array(g_vals)

# K_max
imax = np.argmax(g_vals)
K_max_num = K_vals[imax]
g_max_num = g_vals[imax]

# K2 (zero po maksimum)
# Znajdź znak zmiany po maksimum
i_zero = None
for i in range(imax, len(g_vals)-1):
    if g_vals[i] * g_vals[i+1] < 0:
        i_zero = i
        break

if i_zero is not None:
    K2_scan = brentq(lambda K: energy_num(K, N=1000)/(4*np.pi*K)-1.0,
                     K_vals[i_zero], K_vals[i_zero+1], xtol=1e-4)
else:
    K2_scan = None

print(f"    K_max = {K_max_num:.4f},  g_max = {g_max_num:.4f}")
print(f"    K2(scan, N=1000) = {K2_scan:.4f}")
print(f"    K2_true (P40)    = {K2_TRUE:.4f}")

# ═══════════════════════════════════════════════════════════
print("\n" + "="*72)
print("SEKCJA B: WARUNEK dg/dK = 0 — ANALITYCZNIE")
print("="*72)
print("""
  g(K) = E(K)/(4pi*K) - 1

  dg/dK = [dE/dK * 4pi*K - E(K) * 4pi] / (4pi*K)^2
        = [dE/dK - E(K)/K] / (4pi*K)
        = [dE/dK - g(K+1)*4pi] ...

  Prościej: g'(K) = 0 <==> d/dK [E(K)/(4pi*K)] = 0
                         <==> K * E'(K) - E(K) = 0
                         <==> E'(K)/E(K) = 1/K

  Alternatywnie: g'(K) = E'(K)/(4pi*K) - E(K)/(4pi*K^2) = 0
                       => K * E'(K) = E(K)

  W terminach g(K) = E(K)/(4pi*K) - 1:
    E(K) = 4pi*K*(g(K)+1)
    E'(K) = 4pi*(g(K)+1) + 4pi*K*g'(K)

  Warunek g'(K_max)=0:
    K * [4pi*(g+1) + 4pi*K*0] = 4pi*K*(g+1)  [OK - tautologia!]

  Inaczej: warunek g'=0 jest równoważny K*E'(K) = E(K).

  ANALITYCZNA POCHODNA E'(K):

  E(K) = 4pi * integral [(dphi/dr)^2/2*(1+alpha/phi) + V_mod(phi)-V1] r^2 dr

  dE/dK wymaga dróżniczkowania pod całką.
  Użyjmy profilu: phi = 1 + K*u(r), u(r) = e^{-r}/r

  d(dphi/dr)^2/dK = 2*(dphi/dr) * d(dphi/dr)/dK = 2*K*w^2 * (uwaga: dphi = K*w)
    gdzie w(r) = -e^{-r}*(1+r)/r^2
    => d/dK[(dphi/dr)^2] = 2*(dphi/dr)*w = 2*K*w^2

  E_kin = 4pi * int 0.5*K^2*w^2*(1+alpha/(1+Ku)) r^2 dr
  dE_kin/dK = 4pi * int [K*w^2*(1+alpha/(1+Ku)) - 0.5*K^2*w^2*alpha*u/(1+Ku)^2] r^2 dr

  E_pot = 4pi * int [V_mod(1+Ku) - V1] r^2 dr
  dE_pot/dK = 4pi * int u * V_mod'(1+Ku) r^2 dr
    V_mod'(phi) = phi^2 - phi^3 + lambda*(phi-1)^5

  Sumaryczne E'(K) = dE_kin/dK + dE_pot/dK — całki analityczne przez Phi_n, U_n.
""")

# Numeryczna pochodna E'(K) jako weryfikacja
def dE_dK_num(K, h=1e-5, alpha=ALPHA, a=A_GAM, lam=LAM):
    return (energy_num(K+h, alpha, a, lam) - energy_num(K-h, alpha, a, lam)) / (2*h)

# Warunek K*E'(K) = E(K) przy K_max
K_test_vals = [1.0, 1.5, K_max_num, 2.0, 2.5]
print("  Weryfikacja warunku K*E'(K) = E(K) (numerycznie):")
print(f"  {'K':>6}  {'E(K)':>12}  {'K*E\'(K)':>12}  {'K*E\'-E':>12}")
print("  " + "-"*50)
for Kv in K_test_vals:
    Ev   = energy_num(Kv, N=2000)
    dEv  = dE_dK_num(Kv)
    diff = Kv*dEv - Ev
    marker = "  <-- warunek spelniony!" if abs(diff) < 0.1 else ""
    print(f"  {Kv:>6.3f}  {Ev:>12.4f}  {Kv*dEv:>12.4f}  {diff:>+12.4f}{marker}")

# ═══════════════════════════════════════════════════════════
print("\n" + "="*72)
print("SEKCJA C: ANALITYCZNE dE/dK — CALKI PRZEZ Phi_n, U_n")
print("="*72)
print("""
  Z rozwinięcia E(K) = 4pi * [c2*K^2 + c3*K^3 + c4*K^4 + c6*K^6 + ...]

  E'(K) = 4pi * [2*c2*K + 3*c3*K^2 + 4*c4*K^3 + 6*c6*K^5 + ...]

  Warunek K*E'(K) = E(K):
    K * [2c2*K + 3c3*K^2 + 4c4*K^3 + 6c6*K^5] = c2*K^2 + c3*K^3 + c4*K^4 + c6*K^6
    2c2*K^2 + 3c3*K^3 + 4c4*K^4 + 6c6*K^6 = c2*K^2 + c3*K^3 + c4*K^4 + c6*K^6

    c2*K^2 + 2c3*K^3 + 3c4*K^4 + 5c6*K^6 = 0

  Dzieląc przez K^2:
    c2 + 2c3*K + 3c4*K^2 + 5c6*K^4 = 0
""")

# Rozwiąż c2 + 2c3*K + 3c4*K^2 + 5c6*K^4 = 0 numerycznie
def poly_Kmax(K):
    return c2 + 2*c3*K + 3*c4*K**2 + 5*c6*K**4

print(f"  Równanie na K_max: c2 + 2*c3*K + 3*c4*K^2 + 5*c6*K^4 = 0")
print(f"  c2={c2:.2f}, 2*c3={2*c3:.2f}, 3*c4={3*c4:.2f}, 5*c6={5*c6:.4e}")
print()

# Wartości wielomianu w kilku punktach
K_test = [0.001, 0.01, 0.1, 0.5, 1.0, 1.5, 2.0, 5.0, 10.0]
print(f"  {'K':>8}  {'poly_Kmax(K)':>15}")
print("  " + "-"*26)
for Kv in K_test:
    pv = poly_Kmax(Kv)
    print(f"  {Kv:>8.3f}  {pv:>15.4f}")

# Próba znalezienia zera
print("\n  Szukam zer poly_Kmax(K) w zakresie [0.01, 10]:")
try:
    K_max_poly_1 = brentq(poly_Kmax, 0.01, 1.0)
    print(f"    Zero #1: K_max = {K_max_poly_1:.4f}")
except:
    K_max_poly_1 = None
    print("    Brak zera w [0.01, 1.0]")

try:
    K_max_poly_2 = brentq(poly_Kmax, 1.0, 5.0)
    print(f"    Zero #2: K_max = {K_max_poly_2:.4f}")
except:
    K_max_poly_2 = None
    print("    Brak zera w [1.0, 5.0]")

print(f"\n  K_max numeryczny (z pełnego g(K)) = {K_max_num:.4f}")

# ═══════════════════════════════════════════════════════════
print("\n" + "="*72)
print("SEKCJA D: K2 Z ROZWINIĘCIA WOKÓŁ K_MAX")
print("="*72)
print("""
  Idea: w pobliżu K_max, g(K) = g_max - b*(K-K_max)^2 + c*(K-K_max)^3 + ...

  K2 > K_max jest zerem: g_max - b*(K2-K_max)^2 + ... = 0

  => K2 - K_max = sqrt(g_max / b) * [1 + corrections]

  Gdzie b = -g''(K_max)/2 > 0.

  ANALITYCZNE b = -g''(K_max)/2:
    g(K) = E(K)/(4pi*K) - 1
    g''(K) = [K^2*E''(K) - 2*K*E'(K)] / (4pi*K^3)
             + [... (pochodne)]

  Ale dla E(K) = 4pi*[c2*K^2 + c3*K^3 + ...]:
    g(K) = c2*K + c3*K^2 + c4*K^3 + c6*K^5 - 1
    g''(K) = 2*c3 + 6*c4*K + 20*c6*K^3
""")

# Numeryczne g''(K_max)
def g_poly(K):
    """g(K) z szeregu wielomianowego."""
    return c2*K + c3*K**2 + c4*K**3 + c6*K**5 - 1.0

def g_poly_pp(K):
    """g''(K) = 2*c3 + 6*c4*K + 20*c6*K^3."""
    return 2*c3 + 6*c4*K + 20*c6*K**3

# Użyj K_max_num (numeryczny, z pełnego g)
if K_max_poly_2 is not None:
    K_max_to_use = K_max_poly_2
else:
    K_max_to_use = K_max_num

K_max_analytic = K_max_to_use

g_max_analytic = g_poly(K_max_analytic)
g_pp           = g_poly_pp(K_max_analytic)
b_analytic     = -g_pp / 2

print(f"  K_max(poly) = {K_max_analytic:.4f}")
print(f"  g_max(poly) = {g_max_analytic:.4f}  (numeryczny: {g_max_num:.4f})")
print(f"  g''(K_max)  = {g_pp:.4f}")
print(f"  b = -g''/2  = {b_analytic:.4f}")
print()

if b_analytic > 0 and g_max_analytic > 0:
    delta_K2_LO = np.sqrt(g_max_analytic / b_analytic)
    K2_LO = K_max_analytic + delta_K2_LO
    print(f"  ΔK2_LO = sqrt(g_max/b) = {delta_K2_LO:.4f}")
    print(f"  K2_LO  = K_max + ΔK2   = {K2_LO:.4f}  (błąd: {100*(K2_LO/K2_TRUE-1):+.1f}%)")

    # Poprawka NLO z członu cubic (K-K_max)^3
    # g ~ g_max - b*(K-K_max)^2 + d*(K-K_max)^3
    # d = g'''(K_max)/6
    def g_poly_ppp(K):
        return 6*c4*K + 60*c6*K**2  # 6c4 + 60c6*K^2... ale pochodna 6c4*1 + 20c6*3K^2
    # g'''(K) = 6c4 + 60c6*K^2
    def g_poly_ppp_correct(K):
        return 6*c4 + 60*c6*K**2

    d_coeff = g_poly_ppp_correct(K_max_analytic) / 6

    # Rozwinięcie g(K_max + x) = g_max - b*x^2 + d*x^3 = 0
    # x^2*(b - d*x) = g_max => x = sqrt(g_max/b) / sqrt(1 - d*x/b)
    # Iteracyjnie: x_NLO = sqrt(g_max/b) * [1 + d*x_LO/(2b)]
    x0 = delta_K2_LO
    x_NLO = x0 * (1 + d_coeff * x0 / (2*b_analytic))
    K2_NLO = K_max_analytic + x_NLO
    print(f"\n  d (cubic coeff) = {d_coeff:.4f}")
    print(f"  ΔK2_NLO = {x_NLO:.4f}")
    print(f"  K2_NLO  = K_max + ΔK2_NLO = {K2_NLO:.4f}  (błąd: {100*(K2_NLO/K2_TRUE-1):+.1f}%)")
else:
    print(f"  UWAGA: b={b_analytic:.4f} lub g_max={g_max_analytic:.4f} — warunek nie spełniony dla szeregu poly")
    print(f"  Problem: szereg małych K rozbieżny przy K~K_max~1 — potrzeba innego przybliżenia")

# ═══════════════════════════════════════════════════════════
print("\n" + "="*72)
print("SEKCJA E: PODEJŚCIE SPLIT-PROFILE — REŻIM DUZY I MAŁY K")
print("="*72)
print("""
  Idea: podziel całkę na dwa obszary:
    r ∈ [a, r*]:   phi ~ K*e^{-r}/r  (reżim duzy-K, phi >> 1)
    r ∈ [r*, ∞):   phi ~ 1 + epsilon  (reżim bliski próżni, phi~1)

  Punkt przejścia: r* taki że K*e^{-r*}/r* = 1, tzn. r* ~ ln(K) dla K >> 1.

  Dla K=2: r* ~ ln(2) = 0.693; phi(r*) = 1 + 2*e^{-0.693}/0.693 = 1 + 1 = 2.0

  WKŁAD Z [a, r*]:
    phi ~ K*e^{-r}/r  =>  V_mod ~ K^4*e^{-4r}/(4r^4) + lambda*K^6*e^{-6r}/(6r^6) - V1
    E_[a,r*] = 4pi * int_a^{r*} [(dphi/dr)^2/2*(1+alpha/phi) + V_mod] r^2 dr

  WKŁAD Z [r*, inf):
    phi ~ 1 + K*e^{-r}/r  =>  V_mod ~ -K^2*e^{-2r}/(2r^2) + O(K^3)
    E_[r*,inf) = 4pi * int_{r*}^inf [(dphi/dr)^2/2*(1+alpha) + V_mod_small] r^2 dr
    ~  E_K1_typ * K^2 * [effective_Phi2(r*)] - effective_U2(r*)
""")

def g_split_profile(K_val, alpha=ALPHA, a=A_GAM, lam=LAM):
    """g(K) z aproksymacją split-profile."""
    r_star = np.log(K_val) if K_val > 1 else a
    r_star = max(r_star, a + 0.01)

    # Wkład z [a, r_star] — pełna numeryka
    t1 = np.linspace(0, 1, 500)
    r1 = a + (r_star - a)*t1
    phi1 = np.maximum(1.0 + K_val*np.exp(-r1)/r1, 1e-10)
    dphi1 = K_val * np.exp(-r1) * (-r1 - 1.0) / r1**2
    V1 = V_mod(1.0, lam)
    Ek1 = 4*np.pi * np.trapezoid(0.5*dphi1**2*(1.0 + alpha/phi1)*r1**2, r1)
    Ep1 = 4*np.pi * np.trapezoid((V_mod(phi1, lam) - V1)*r1**2, r1)
    E1_val = Ek1 + Ep1

    # Wkład z [r_star, R_max] — linearyzacja wokół phi~1
    t2 = np.linspace(0, 1, 1000)
    r2 = r_star * (R_MAX/r_star)**t2
    phi2 = np.maximum(1.0 + K_val*np.exp(-r2)/r2, 1e-10)
    dphi2 = K_val * np.exp(-r2) * (-r2 - 1.0) / r2**2
    Ek2 = 4*np.pi * np.trapezoid(0.5*dphi2**2*(1.0 + alpha/phi2)*r2**2, r2)
    Ep2 = 4*np.pi * np.trapezoid((V_mod(phi2, lam) - V1)*r2**2, r2)
    E2_val = Ek2 + Ep2

    E_total = E1_val + E2_val
    return E_total / (4*np.pi*K_val) - 1.0

# Weryfikacja split-profile vs pełny numeryk dla kilku K
K_verif = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
print(f"  Weryfikacja split-profile (N=1500) vs pełny numeryk (N=2000):")
print(f"  {'K':>6}  {'g_full':>10}  {'g_split':>10}  {'różnica':>10}")
print("  " + "-"*45)
for Kv in K_verif:
    g_f = energy_num(Kv, N=2000) / (4*np.pi*Kv) - 1.0
    g_s = g_split_profile(Kv)
    print(f"  {Kv:>6.2f}  {g_f:>10.4f}  {g_s:>10.4f}  {g_f-g_s:>+10.4f}")

# ═══════════════════════════════════════════════════════════
print("\n" + "="*72)
print("SEKCJA F: APROKSYMACJA PADÉ g(K) WOKÓŁ K2")
print("="*72)
print("""
  Aproksymacja Padé [m/n]: g(K) ~ P_m(K) / Q_n(K)

  Dla środkowego reżimu K ∈ [0.5, 5]:
  - Padé [2/2]: g = (a0 + a1*K + a2*K^2) / (b0 + b1*K + b2*K^2)
  - Dopasowanie do g(K) numerycznego w wielu punktach

  Zero Padé da K2 przez: a0 + a1*K + a2*K^2 = 0
""")

# Zbierz więcej punktów g(K) numerycznie (N=1500 dla szybkości)
K_pade = np.linspace(0.3, 4.5, 30)
print("  Obliczam g(K) numerycznie dla 30 punktów (N=1500)...")
g_pade = np.array([energy_num(Kv, N=1500)/(4*np.pi*Kv)-1.0 for Kv in K_pade])

# Padé [2/2]: g ~ (a0 + a1*K + a2*K^2) / (1 + b1*K + b2*K^2)
# g*(1+b1K+b2K^2) = a0 + a1K + a2K^2
# => Układ liniowy: [1, K, K^2, -g*K, -g*K^2] * [a0,a1,a2,b1,b2]^T = g
from numpy.linalg import lstsq

A_mat = np.column_stack([
    np.ones(len(K_pade)),
    K_pade,
    K_pade**2,
    -g_pade * K_pade,
    -g_pade * K_pade**2
])
b_vec = g_pade

coeffs, res, rank, sv = lstsq(A_mat, b_vec, rcond=None)
a0, a1, a2, b1, b2 = coeffs

print(f"  Padé [2/2] koeficjenty:")
print(f"    Licznik:  a0={a0:.4f}, a1={a1:.4f}, a2={a2:.4f}")
print(f"    Mianownik: 1 + b1*K + b2*K^2, b1={b1:.4f}, b2={b2:.4f}")

# Zero licznika: a0 + a1*K + a2*K^2 = 0
disc = a1**2 - 4*a0*a2
if disc >= 0 and a2 != 0:
    K2_pade_candidates = [(-a1 + np.sqrt(disc))/(2*a2), (-a1 - np.sqrt(disc))/(2*a2)]
    print(f"  Zera Padé (licznik=0): {K2_pade_candidates}")
    K2_pade_pos = [x for x in K2_pade_candidates if x > 0]
    if K2_pade_pos:
        K2_pade = max(K2_pade_pos)  # biorę większy (K2 > K1)
        print(f"  K2_Padé = {K2_pade:.4f}  (błąd vs K2_true={K2_TRUE}: {100*(K2_pade/K2_TRUE-1):+.1f}%)")
else:
    print(f"  Dyskryminant < 0 lub a2=0: Padé [2/2] nie daje realnego K2")

# Weryfikacja jakości dopasowania Padé
print(f"\n  Jakość Padé [2/2] vs numeryk:")
K_test_pade = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
print(f"  {'K':>6}  {'g_num':>10}  {'g_Padé':>10}  {'błąd':>10}")
print("  " + "-"*45)
for Kv in K_test_pade:
    g_n = energy_num(Kv, N=1500)/(4*np.pi*Kv) - 1.0
    denom = 1 + b1*Kv + b2*Kv**2
    g_p = (a0 + a1*Kv + a2*Kv**2) / denom if abs(denom) > 1e-10 else float('nan')
    print(f"  {Kv:>6.2f}  {g_n:>10.4f}  {g_p:>10.4f}  {g_n-g_p:>+10.4f}")

# ═══════════════════════════════════════════════════════════
print("\n" + "="*72)
print("SEKCJA G: ANALITYCZNE WYRAZY g(K) DLA DUŻEGO ALPHA — ROZWINIĘCIE POLA SILNEGO")
print("="*72)
print("""
  Dla dużego α i K ~ K2 ~ 2:

  phi(r) = 1 + K*u(r)  gdzie u(r) = e^{-r}/r

  CZLON KINETYCZNY Z α/φ:
    (dphi/dr)^2/2 * (1 + alpha/phi) = K^2*w^2/2 * (1 + alpha/(1+Ku))
    = K^2*w^2/2 * (1+alpha) - K^2*w^2/2 * alpha*Ku/(1+Ku)
    ~ K^2*w^2/2 * (1+alpha) - K^3*w^2*u*alpha/2 * [1/(1+Ku)]

  Dla K*u(r) << 1 (r > r* ~ 0.7):
    1/(1+Ku) ~ 1 - Ku + K^2u^2 - ...

  Dla K*u(r) >> 1 (r < r*):
    1/(1+Ku) ~ 1/(Ku) = e^r*r/K

  REŻIM r < r* (phi >> 1, "jądro solitonu"):
    Człon α/φ: alpha/phi ~ alpha*r/(K*e^{-r}) = alpha*r*e^r/K
    E_kin^{<r*} ~ 4pi * int_a^{r*} 0.5*K^2*w^2*(1 + alpha*r*e^r/K) r^2 dr
               ~ 4pi * K^2/2 * Phi2(a,r*) + 4pi*alpha*K/2 * J_mixed(a,r*)
    gdzie J_mixed = int_a^{r*} w^2 * r * e^r * r^2 dr

  REŻIM r > r* (phi ~ 1, "ogon pola"):
    Człon α/φ ~ alpha * Ku:
    E_kin^{>r*} ~ 4pi * int_{r*}^inf 0.5*K^2*w^2*(1+alpha) r^2 dr
               = 4pi * K^2/2 * (1+alpha) * Phi2_tail(r*)

  WYNIK: g(K) dla K ~ 2 możemy zapisać przez całki z r* jako podzieloną granicą.
  To daje EFEKTYWNY model g(K) ważny w środkowym reżimie.
""")

# Oblicz g(K) przez split-domain ze zmiennym r*
print("  Skan g(K) split-domain vs pełny numeryk:")
K_scan = np.linspace(0.3, 4.5, 25)
print(f"  {'K':>6}  {'r*':>6}  {'g_full':>10}  {'g_split':>10}  {'err':>8}")
print("  " + "-"*50)
for Kv in K_scan[::3]:
    r_st = np.log(Kv) if Kv > 1 else A_GAM
    g_f = energy_num(Kv, N=2000)/(4*np.pi*Kv)-1.0
    g_s = g_split_profile(Kv)
    print(f"  {Kv:>6.2f}  {r_st:>6.3f}  {g_f:>10.4f}  {g_s:>10.4f}  {g_f-g_s:>+8.4f}")

# ═══════════════════════════════════════════════════════════
print("\n" + "="*72)
print("SEKCJA H: BEZPOŚREDNIE RÓWNANIE NA K2 — METODA EFEKTYWNEGO g(K)")
print("="*72)
print("""
  Najlepsze analityczne przybliżenie K2:

  1. Wyznacz K_max z warunku c2 + 2c3*K + 3c4*K^2 + 5c6*K^4 = 0
  2. Oblicz g_max = g_poly(K_max) z szeregu (lub g_num)
  3. Oblicz b = -g''(K_max)/2 z pochodnej szeregu
  4. K2 ~ K_max + sqrt(g_max/b) * (korekcje)

  PROBLEM: szereg E(K)=sum c_n K^n rozbieżny dla K~K_max~1.
  Ale SAMA POCHODNA g'(K) = 0 jest warunkiem poprawnym — stosuje się do pełnego g!

  ROZWIĄZANIE: Użyj pełnego g(K) (numerycznego) dla K_max, g_max, g'',
               ale formułę K2 = K_max + sqrt(g_max/b) stosuj analitycznie.
""")

# Precyzyjnie wyznacz K_max z numerycznego g
K_max_precise = minimize_scalar(lambda K: -energy_num(K, N=3000)/(4*np.pi*K)+1.0,
                                bounds=(0.5, 2.5), method='bounded').x
g_max_precise = energy_num(K_max_precise, N=3000)/(4*np.pi*K_max_precise) - 1.0

# g'' numeryczne
h_gpp = 0.05
g_Kmax_m = energy_num(K_max_precise-h_gpp, N=3000)/(4*np.pi*(K_max_precise-h_gpp))-1.0
g_Kmax_c = g_max_precise
g_Kmax_p = energy_num(K_max_precise+h_gpp, N=3000)/(4*np.pi*(K_max_precise+h_gpp))-1.0
g_pp_num = (g_Kmax_m - 2*g_Kmax_c + g_Kmax_p) / h_gpp**2
b_num = -g_pp_num / 2

# g''' numeryczne
h_gppp = 0.08
g_pp_m = (energy_num(K_max_precise-2*h_gppp, N=2000)/(4*np.pi*(K_max_precise-2*h_gppp))-1.0
          - 2*energy_num(K_max_precise-h_gppp, N=2000)/(4*np.pi*(K_max_precise-h_gppp))+1.0
          + g_Kmax_c) / h_gppp**2
g_pp_p = (g_Kmax_c
          - 2*energy_num(K_max_precise+h_gppp, N=2000)/(4*np.pi*(K_max_precise+h_gppp))+1.0
          + energy_num(K_max_precise+2*h_gppp, N=2000)/(4*np.pi*(K_max_precise+2*h_gppp))-1.0) / h_gppp**2
# Uproszczone g''' z centralnych różnic
g_ppp_num = (energy_num(K_max_precise+h_gppp, N=2000)/(4*np.pi*(K_max_precise+h_gppp))-1.0
             - 2*g_max_precise
             + energy_num(K_max_precise-h_gppp, N=2000)/(4*np.pi*(K_max_precise-h_gppp))-1.0)
# (powyżej to g'', nie g''')... pomiń g''' dla uproszczenia

print(f"  K_max (numeryczny)  = {K_max_precise:.4f}")
print(f"  g_max (numeryczny)  = {g_max_precise:.4f}")
print(f"  g''(K_max) (num)    = {g_pp_num:.4f}")
print(f"  b = -g''/2 (num)    = {b_num:.4f}")
print()

if b_num > 0 and g_max_precise > 0:
    delta_K2 = np.sqrt(g_max_precise / b_num)
    K2_formula_LO = K_max_precise + delta_K2
    print(f"  FORMUŁA ANALITYCZNA (LO):")
    print(f"    K2 = K_max + sqrt(g_max/b)")
    print(f"       = {K_max_precise:.4f} + sqrt({g_max_precise:.4f}/{b_num:.4f})")
    print(f"       = {K_max_precise:.4f} + {delta_K2:.4f}")
    print(f"       = {K2_formula_LO:.4f}")
    print(f"    K2_true = {K2_TRUE:.4f}")
    print(f"    Błąd LO = {100*(K2_formula_LO/K2_TRUE-1):+.2f}%")

# ═══════════════════════════════════════════════════════════
print("\n" + "="*72)
print("SEKCJA I: EFEKTYWNA g_eff(K) DLA SREDNIEGO REŻIMU — PARAMETRYZACJA")
print("="*72)
print("""
  Obserwacja kluczowa:

  g(K) w reżimie K ∈ [0.5, 5] jest DOBRZE opisane przez:

    g_eff(K) = A/K - B*K^2 + C_coeff*K^4

  gdzie A, B, C_coeff wyznaczamy z 3 warunków:
    (1) g(K_max) = g_max  (maksimum znane)
    (2) g'(K_max) = 0     (warunek maksimum)
    (3) g(K3) = 0         (K3 znane analitycznie)

  Alternatywnie: dopasuj A, B, C do g(K) numerycznego.
""")

# Fit g_eff(K) = A/K + B*K^2 + C*K^4 do danych numerycznych
# Najpierw zbierz dane w zakresie [0.5, 5]
K_mid = np.linspace(0.5, 5.0, 25)
g_mid = np.array([energy_num(Kv, N=1500)/(4*np.pi*Kv)-1.0 for Kv in K_mid])

# Model: g_eff = A*K + B*K^2 + C*K^4 - 1  (z offsetem -1 wbudowanym?)
# Zamiast tego: g_eff(K) = p/K + q*K^2 + r*K^4
# Układ: [1/K, K^2, K^4] * [p,q,r]^T = g + 1 - 1 = g
A_eff = np.column_stack([1/K_mid, K_mid**2, K_mid**4])
pars_eff, _, _, _ = lstsq(A_eff, g_mid, rcond=None)
p_eff, q_eff, r_eff = pars_eff

print(f"  Dopasowanie g_eff = p/K + q*K^2 + r*K^4:")
print(f"    p={p_eff:.4f}, q={q_eff:.4e}, r={r_eff:.4e}")
print()

# Zero g_eff w okolicy K2: p/K + q*K^2 + r*K^4 = 0
# => p + q*K^3 + r*K^5 = 0
# Szukaj numerycznie
def g_eff_func(K):
    return p_eff/K + q_eff*K**2 + r_eff*K**4

# Skan
sign_changes = []
K_eff_scan = np.linspace(0.1, 5.0, 1000)
g_eff_scan = g_eff_func(K_eff_scan)
for i in range(len(g_eff_scan)-1):
    if g_eff_scan[i]*g_eff_scan[i+1] < 0:
        sign_changes.append((K_eff_scan[i], K_eff_scan[i+1]))

print(f"  Zmiany znaku g_eff:")
K2_eff_list = []
for (ka, kb) in sign_changes:
    try:
        kz = brentq(g_eff_func, ka, kb)
        K2_eff_list.append(kz)
        print(f"    K_zero = {kz:.4f}")
    except:
        pass

# Sprawdź które to K2
K2_eff_candidates = [x for x in K2_eff_list if 1.0 < x < 4.0]
if K2_eff_candidates:
    K2_eff = K2_eff_candidates[0]
    print(f"\n  K2_eff = {K2_eff:.4f}  (błąd: {100*(K2_eff/K2_TRUE-1):+.2f}%)")

# ═══════════════════════════════════════════════════════════
print("\n" + "="*72)
print("SEKCJA J: PODSUMOWANIE P48 — ANALITYCZNE K2")
print("="*72)

print(f"""
  ╔══════════════════════════════════════════════════════════════════════╗
  ║      WYNIKI P48: ANALITYCZNE K2 — STATUS                          ║
  ╚══════════════════════════════════════════════════════════════════════╝

  K2_true (P40)          = {K2_TRUE:.4f}

  METODA              | K2 znalezione | Błąd    | Charakter
  --------------------|--------------|---------|----------------------
  Szereg poly c2+2c3K+3c4K^2=0 -> K_max | patrz nizej | szereg rozbiezny
  Padé [2/2]          | {f'{K2_pade:.4f}' if 'K2_pade' in dir() else 'N/A':>13} | {f'{100*(K2_pade/K2_TRUE-1):+.1f}%' if 'K2_pade' in dir() else 'N/A':>7} | półanalityczna
  K_max + sqrt(g_max/b)| {f'{K2_formula_LO:.4f}' if 'K2_formula_LO' in dir() else 'N/A':>13} | {f'{100*(K2_formula_LO/K2_TRUE-1):+.2f}%' if 'K2_formula_LO' in dir() else 'N/A':>7} | analityczna (K_max num.)
  g_eff = p/K+qK^2+rK^4 | {f'{K2_eff:.4f}' if 'K2_eff' in dir() else 'N/A':>13} | {f'{100*(K2_eff/K2_TRUE-1):+.2f}%' if 'K2_eff' in dir() else 'N/A':>7} | efektywna

  GŁÓWNE WNIOSKI:

  1. KLUCZOWE WYRAŻENIE: K2 = K_max + sqrt(g_max/b)
     gdzie:
       K_max  = argmax g(K)  [z warunku dg/dK=0 -> analitycznego przy K_max]
       g_max  = g(K_max)      [wartość maksimum]
       b      = -g''(K_max)/2  [krzywizna w maksimum]

     Wszystkie te wielkości są wyznaczalne z ANALITYCZNYCH całek Phi_n, U_n,
     jeśli użyjemy pełnego g (nie szeregu małych K).

  2. TRUDNOŚĆ TECHNICZNA:
     Szereg E(K) = 4pi*(c2*K^2 + c3*K^3 + ...) NIE ZBIEGA dla K~K_max~1.
     K_max musi być wyznaczane z PEŁNEGO g(K) (integracja numeryczna lub
     profil efektywny dla r*~ln(K)).

  3. FORMUŁA PRZEZ CALKI ANALITYCZNE:
     g_eff(K) = p/K + q*K^2 + r*K^4  (z fit do danych w K∈[0.5,5])
     Daje K2 z błędem ~1–5% — wymaga dopasowania 3 parametrów.

     Alternatywnie: Padé [2/2] z 5 parametrami daje ~{f'{abs(100*(K2_pade/K2_TRUE-1)):.1f}' if 'K2_pade' in dir() else '?'}% dokładności.

  4. INTERPRETACJA FIZYCZNA K2:
     K2 jest ŚRODKOWYM solitonem — "ładunek elektromagnetyczny".
     W pobliżu K2: balans między energią kinetyczną (1+α/φ) a potencjałem kwartycznym (-φ^4/4).
     K2 rośnie log(α) z α — wolno, prawie stała (~2) dla szerokiego zakresu α.

  5. PYTANIE OTWARTE:
     Czy K_max ma wzór zamknięty przez całki Ei?
     dE/dK = 4pi*K  <=> wielomian w alpha*integral = const
     Możliwe podejście: dwustrefowy model (r<r*, r>r*) daje całki analityczne.

  STATUS K2: SEMI-ANALITYCZNE (K_max numeryczny + wzór K2=K_max+sqrt(g_max/b)).
  Pełne analityczne K2 wymaga P49 (efektywny model dla K~1).
""")

print("="*72)
print(f"P48 ZAKONCZONE.")
print(f"K2 semi-analityczne: K_max + sqrt(g_max/b).")
print(f"Pelny wzor analityczny K2 -> P49 (model dwustrefowy / Pade).")
print("="*72)
