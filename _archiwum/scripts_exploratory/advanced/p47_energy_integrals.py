#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')
"""
p47_energy_integrals.py
========================
Analiza P47: Analityczne calki energii TGP — wyprowadzenie K1, K3 z pierwszych zasad

Cel: Wyprowadzic analitycznie K1, K2, K3 z funktionalu energii TGP,
     bez zadnych zalozonych wartosci numerycznych.

Punkt wyjscia (Lagrangian TGP):
  E[phi] = 4pi * int_a^inf [ 0.5*(d phi/dr)^2 * (1 + alpha/phi)
                             + V_mod(phi) - V_mod(1) ] * r^2 dr

  phi(r) = 1 + K * e^{-r}/r   (profil Yukawa — ansatz ruchliwosci solitonu)
  V_mod(phi) = phi^3/3 - phi^4/4 + lambda*(phi-1)^6/6

  Warunek samokonzystencji (P28):  g(K) = E(K)/(4*pi*K) - 1 = 0

Podejscie analityczne (trzy rezimy):
  A. Rozwiniciecze E(K) w szereg potegowy w K — calki I_n(a) wyrazalne przez Ei
  B. Rezimy K:
     K << 1  (K1):  g ~ -1 + K*A1(alpha,a) + K^2*A2(alpha,a)   => K1 analitycznie
     K >> 1  (K3):  g ~ K^5*C5(a,lambda) - K^3*C3(a) + ...     => K3 = C*a/sqrt(lambda)
     K ~ K2:        g = 0 z warunku lokalnego ekstremum
  C. Calki Phi_n i U_n — wyrazenie przez Ei(x) = -E1(-x)
  D. Porownanienumeryczne: formula analityczna vs numeryk

Kluczowe calki:
  Phi_n(a) = int_a^inf e^{-n*r}*(1+r)^2/r^n dr   [z czlonu kinetycznego]
  U_n(a)   = int_a^inf e^{-n*r}/r^{n-2} dr        [z czlonu potencjalnego]

Rozwinicieenergii:
  E(K) = sum_n c_n(alpha, a, lambda) * K^n

  c2 = 2*pi*(1+alpha)*Phi_2 - 4*pi*U_2/2       = 4*pi * [(1+alpha)*Phi_2/2 - U_2/2]
  c3 = -2*pi*alpha*Phi_3 + 4*pi*(-2*U_3/3)     = 4*pi * [-alpha*Phi_3/2 - 2*U_3/3]
  c4 = 2*pi*alpha*Phi_4 + 4*pi*(-U_4/4)        = 4*pi * [alpha*Phi_4/2 - U_4/4]
  c5 = ...                                       ...
  c6 = 4*pi * lambda * U_6/6                     [czlon seksyczny]
"""

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import expi   # Ei(x) = -E1(-x) dla x>0

# =========================================================================
# PARAMETRY TGP
# =========================================================================
R_MAX = 60.0
a     = 0.040      # a_Gamma
alpha_lep = 8.5445
lam_K     = 5.4677e-6
EULER_MAS = 0.5772156649  # stala Eulera-Mascheroniego

# =========================================================================
# DEFINICJE CALEK ANALITYCZNYCH
# =========================================================================

def Phi_n_numeric(n, a_val, R=50.0):
    """
    Phi_n = int_a^R e^{-n*r}*(1+r)^2/r^n dr
    Calka numeryczna (referencyjna).
    """
    def f(r):
        return np.exp(-n*r)*(1+r)**2 / r**n
    result, _ = quad(f, a_val, R, limit=200)
    return result

def U_n_numeric(n, a_val, R=50.0):
    """
    U_n = int_a^R e^{-n*r}/r^{n-2} dr
    Calka numeryczna (referencyjna).
    """
    def f(r):
        return np.exp(-n*r) / r**(n-2)
    result, _ = quad(f, a_val, R, limit=200)
    return result

def Ei_safe(x):
    """Ei(x) = integral representation; dla x<0: Ei(x) = -E1(-x)."""
    if x > 0:
        return expi(x)
    else:
        return np.nan

def E1(x):
    """E1(x) = int_x^inf e^{-t}/t dt dla x>0."""
    return -expi(-x)

# =========================================================================
# SEKCJA A: Calki Phi_n i U_n — formy analityczne
# =========================================================================
def section_A():
    print("\n" + "=" * 72)
    print("SEKCJA A: CALKI KLUCZOWE Phi_n I U_n — WYRAZENIA ANALITYCZNE")
    print("=" * 72)
    print("""
  Profil Yukawa: phi(r) = 1 + K*exp(-r)/r
  dPhi/dr = -K*exp(-r)*(1+r)/r^2

  CALKI KINETYCZNE: Phi_n(a) = int_a^inf exp(-nr)*(1+r)^2/r^n dr

  Rozwiniecieexponencjalne dla malego a:
    exp(-nr) = exp(-n*a) * exp(-n*(r-a)) ~ 1 - nr dla r~a

  Rozwinice (1+r)^2 = 1 + 2r + r^2:

    Phi_n(a) = int_a^inf exp(-nr)/r^n dr + 2*int_a^inf exp(-nr)/r^{n-1} dr
               + int_a^inf exp(-nr)/r^{n-2} dr

  Oznaczamy: J_k(a,n) = int_a^inf e^{-nr}/r^k dr

  J_k(a,n) = n^{k-1} * E_k(na)     [uogolniona calka wykladnicza]
  gdzie E_k(z) = int_1^inf e^{-zt}/t^k dt

  Dla n=2:
    J_2(a,2) = int_a^inf e^{-2r}/r^2 dr = 2*E1(2a) - e^{-2a}/a
                                         ~ 1/a - 2*ln(2a) - 2*gamma + O(a)
    J_1(a,2) = int_a^inf e^{-2r}/r dr   = E1(2a)
                                         ~ -ln(2a) - gamma + O(a)
    J_0(a,2) = int_a^inf e^{-2r} dr     = e^{-2a}/2

  => Phi_2(a) = J_2 + 2*J_1 + J_0
              ~ (1/a - 2*ln(2a) - 2*gamma) + 2*(-ln(2a) - gamma) + 1/2
              = 1/a - 4*ln(2a) - 4*gamma + 1/2
    """)

    print(f"  Weryfikacja Phi_2 dla a={a}:")
    Phi2_analytic = 1/a - 4*np.log(2*a) - 4*EULER_MAS + 0.5
    Phi2_numeric  = Phi_n_numeric(2, a)
    # Bardziej dokladna formula:
    Phi2_exact = E1(2*a) - np.exp(-2*a)/a + 2*E1(2*a) + np.exp(-2*a)/2
    # Phi2 = J2 + 2*J1 + J0
    J2 = 2*E1(2*a) - np.exp(-2*a)/a   # int_a^inf e^{-2r}/r^2 dr = d/da[-E1(2a)] - ...
    # Actually: int_a^inf e^{-2r}/r^2 dr  using IBP: let u=e^{-2r}, dv=dr/r^2
    # = [-e^{-2r}/r]_a^inf + int_a^inf (-2e^{-2r}/r)(-1) dr
    # = e^{-2a}/a - 2 * E1(2a)
    # Wait, E1(x) = int_x^inf e^{-t}/t dt, so E1(2a) = int_2a^inf e^{-t}/t dt
    # int_a^inf e^{-2r}/r dr = int_2a^inf e^{-s}/(s/2) ds/2 = int_2a^inf e^{-s}/s ds = E1(2a)
    J2_v2 = np.exp(-2*a)/a - 2*E1(2*a)
    J1_v2 = E1(2*a)
    J0_v2 = np.exp(-2*a)/2
    Phi2_exact_v2 = J2_v2 + 2*J1_v2 + J0_v2

    print(f"    Phi_2 numerycznie    = {Phi2_numeric:.6f}")
    print(f"    Phi_2 przyblizenie   = {Phi2_analytic:.6f}  (blad: {(Phi2_analytic-Phi2_numeric)/Phi2_numeric*100:+.2f}%)")
    print(f"    Phi_2 dokladna (Ei)  = {Phi2_exact_v2:.6f}  (blad: {(Phi2_exact_v2-Phi2_numeric)/Phi2_numeric*100:+.2f}%)")

    print(f"""
  CALKI POTENCJALNE: U_n(a) = int_a^inf e^{{-nr}}/r^{{n-2}} dr

  U_2(a) = int_a^inf e^{{-2r}} dr         = e^{{-2a}}/2
  U_3(a) = int_a^inf e^{{-3r}}/r dr       = E1(3a)
  U_4(a) = int_a^inf e^{{-4r}}/r^2 dr    = e^{{-4a}}/a - 4*E1(4a)
  U_6(a) = int_a^inf e^{{-6r}}/r^4 dr    = [formula przez Ei i potegi]
    """)

    print(f"  Weryfikacja dla a={a}:")
    for n, name in [(2, 'U_2'), (3, 'U_3'), (4, 'U_4'), (6, 'U_6')]:
        U_num = U_n_numeric(n, a)
        if n == 2:
            U_ana = np.exp(-2*a)/2
        elif n == 3:
            U_ana = E1(3*a)
        elif n == 4:
            U_ana = np.exp(-4*a)/a - 4*E1(4*a)
        elif n == 6:
            # U_6 = int_a^inf e^{-6r}/r^4 dr
            # IBP: = [-e^{-6r}/(3r^3)]_a^inf + int_a^inf (-6e^{-6r}/(3r^3))(-1/(3-1))...
            # Bardziej systematycznie: int e^{-6r}/r^4 dr
            # Iterujac IBP 3 razy:
            # int_a^inf e^{-6r}/r^4 dr = e^{-6a}/(3a^3) - (6/(3)) * int_a^inf e^{-6r}/r^3 dr
            #                          = e^{-6a}/(3a^3) - 2*[e^{-6a}/(2a^2) - 3*int_a^inf e^{-6r}/r^2 dr]
            #                          = e^{-6a}/(3a^3) - e^{-6a}/a^2 + 6*[e^{-6a}/a - 6*E1(6a)]
            U_ana = (np.exp(-6*a)/(3*a**3)
                     - np.exp(-6*a)/a**2
                     + 6*(np.exp(-6*a)/a - 6*E1(6*a)))
        else:
            U_ana = np.nan
        err = (U_ana - U_num)/U_num*100 if abs(U_num) > 1e-15 else 0
        print(f"    {name}(a={a}) numerycznie={U_num:.6e},  analitycznie={U_ana:.6e},  blad={err:+.4f}%")

    return J2_v2, J1_v2, J0_v2


# =========================================================================
# SEKCJA B: Rozwiniciecze E(K) — wspolczynniki c_n
# =========================================================================
def section_B():
    print("\n" + "=" * 72)
    print("SEKCJA B: ROZWINICIECZE E(K) = sum c_n * K^n — ANALITYCZNE c_n")
    print("=" * 72)
    print("""
  Substitucja phi(r) = 1 + K*u(r), u(r) = e^{-r}/r:

  Potencjal:  V_mod(phi) - V_mod(1) = -K^2*u^2/2 - 2K^3*u^3/3 - K^4*u^4/4
                                       + lambda*K^6*u^6/6   + O(K^8)
  Kinetyczny: 0.5*(dphi/dr)^2 * (1+alpha/phi)
    = 0.5*K^2*w^2 * (1+alpha/(1+K*u))
    gdzie w = dphi/dr / K = -e^{-r}*(1+r)/r^2

  Rozwiniciecze (1+alpha/(1+Ku)) w K:
    1 + alpha/(1+Ku) = (1+alpha) - alpha*K*u + alpha*K^2*u^2 - alpha*K^3*u^3 + ...

  =>
    0.5*K^2*w^2*(1+alpha) * [ czlon K^2 ]
    - 0.5*alpha*K^3*w^2*u  [ czlon K^3 ]
    + 0.5*alpha*K^4*w^2*u^2[ czlon K^4 ]

  Po calkowaniu 4*pi*int[...] r^2 dr:

  E(K) = 4*pi*[ c2*K^2 + c3*K^3 + c4*K^4 + c5*K^5 + c6*K^6 + ...]

  gdzie:
    c2 = (1+alpha)*Phi_2/2 - U_2/2
    c3 = -alpha*Phi_3/2 - 2*U_3/3
    c4 = alpha*Phi_4/2 - U_4/4
    c6 = lambda*U_6/6

  (czlon c5 pochodzi z czlonu lambda*(phi-1)^5... ale V_mod nie ma phi^5, wiec c5=0)
    """)

    # Oblicz wspolczynniki numerycznie dla a=0.040, alpha=8.5445, lambda=5.4677e-6
    print(f"  Parametry: a={a}, alpha={alpha_lep}, lambda={lam_K:.4e}")
    print()

    Phi2 = Phi_n_numeric(2, a)
    Phi3 = Phi_n_numeric(3, a)
    Phi4 = Phi_n_numeric(4, a)

    U2 = U_n_numeric(2, a)
    U3 = U_n_numeric(3, a)
    U4 = U_n_numeric(4, a)
    U6 = U_n_numeric(6, a)

    c2 = (1+alpha_lep)*Phi2/2 - U2/2
    c3 = -alpha_lep*Phi3/2 - 2*U3/3
    c4 = alpha_lep*Phi4/2 - U4/4
    c6 = lam_K*U6/6

    print(f"  Calki (numeryczne, referencyjna weryfikacja):")
    print(f"    Phi_2 = {Phi2:.6f}")
    print(f"    Phi_3 = {Phi3:.6f}")
    print(f"    Phi_4 = {Phi4:.6f}")
    print(f"    U_2   = {U2:.6f}")
    print(f"    U_3   = {U3:.6f}")
    print(f"    U_4   = {U4:.6f}")
    print(f"    U_6   = {U6:.6e}")
    print()
    print(f"  Wspolczynniki energii E(K) = 4*pi * sum c_n * K^n:")
    print(f"    c2 = (1+alpha)*Phi2/2 - U2/2     = {c2:.6f}")
    print(f"    c3 = -alpha*Phi3/2 - 2*U3/3      = {c3:.6f}")
    print(f"    c4 = alpha*Phi4/2 - U4/4          = {c4:.6f}")
    print(f"    c6 = lambda*U6/6                   = {c6:.6e}")

    # Weryfikacja: E(K) = 4*pi*(c2*K^2 + c3*K^3 + c4*K^4 + c6*K^6)
    # vs. numeryczna calkowanie dla roznych K
    print(f"\n  Weryfikacja E(K)/4pi: szereg vs numeryk:")
    print(f"  {'K':>8}  {'E_szereg/4pi':>14}  {'E_numer/4pi':>14}  {'blad':>8}")
    print(f"  {'-'*52}")

    def E_numeric(K):
        N = 2000
        t = np.linspace(0, 1, N)
        r = a * (R_MAX/a)**t
        phi = np.maximum(1 + K*np.exp(-r)/r, 1e-10)
        dphi = K*np.exp(-r)*(-r-1)/r**2
        V1 = 1/3 - 1/4  # V_mod(1)
        Ek = np.trapezoid(0.5*dphi**2*(1+alpha_lep/phi)*r**2, r)
        Ep = np.trapezoid(((phi**3/3 - phi**4/4 + lam_K*(phi-1)**6/6) - V1)*r**2, r)
        return Ek + Ep

    for K_test in [0.001, 0.005, 0.01, 0.05, 0.1, 0.5]:
        E_ser = c2*K_test**2 + c3*K_test**3 + c4*K_test**4 + c6*K_test**6
        E_num = E_numeric(K_test)
        err = (E_ser - E_num)/abs(E_num)*100 if abs(E_num) > 1e-15 else 0
        print(f"  {K_test:>8.4f}  {E_ser:>14.6e}  {E_num:>14.6e}  {err:>+8.2f}%")

    print(f"""
  UWAGA: Dla malych K szereg K^2 + K^3 + K^4 szybko sie zbiega.
  Dla duzych K (K~K3~34) szereg jest rozbiezny — potrzebne inne podejscie.
    """)

    return c2, c3, c4, c6, Phi2, Phi3, Phi4, U2, U3, U4, U6


# =========================================================================
# SEKCJA C: Analityczne K1 — male K
# =========================================================================
def section_C(c2, c3, c4):
    print("\n" + "=" * 72)
    print("SEKCJA C: ANALITYCZNE K1 — APROKSYMACJA MALYCH K")
    print("=" * 72)
    print(f"""
  g(K) = E(K)/(4*pi*K) - 1 = c2*K + c3*K^2 + c4*K^3 + ... - 1

  Rownanie g(K1) = 0:
    c2*K1 + c3*K1^2 + c4*K1^3 + ... = 1

  Rzad LO (c2*K1 = 1):
    K1^(LO) = 1/c2

  Rzad NLO (c2*K + c3*K^2 = 1):
    K1^(NLO) = (-c2 + sqrt(c2^2 + 4*c3)) / (2*c3)
               [z rownania kwadratowego, z wlasciwym znakiem]

  Rzad NNLO (c2*K + c3*K^2 + c4*K^3 = 1):
    Rozwiazujemy numerycznie.
    """)

    K1_LO  = 1.0 / c2
    K1_NLO = (-c2 + np.sqrt(c2**2 + 4*c3)) / (2*c3)
    # NNLO: numerycznie
    def g_trunc(K):
        return c2*K + c3*K**2 + c4*K**3 - 1.0
    try:
        K1_NNLO = brentq(g_trunc, 1e-6, K1_LO*3, xtol=1e-10)
    except Exception:
        K1_NNLO = np.nan

    K1_true = 0.009839  # z P40

    print(f"  c2 = {c2:.6f}")
    print(f"  c3 = {c3:.6f}")
    print(f"  c4 = {c4:.6f}")
    print()
    print(f"  K1^(LO)   = 1/c2                = {K1_LO:.6f}  (blad: {(K1_LO-K1_true)/K1_true*100:+.2f}%)")
    print(f"  K1^(NLO)  = z rownania kwadr.    = {K1_NLO:.6f}  (blad: {(K1_NLO-K1_true)/K1_true*100:+.2f}%)")
    print(f"  K1^(NNLO) = z rownania sz. 3     = {K1_NNLO:.6f}  (blad: {(K1_NNLO-K1_true)/K1_true*100:+.2f}%)")
    print(f"  K1^(true) = {K1_true:.6f}  (numeryk P40)")

    # Analityczna formula explicite dla K1^NLO:
    print(f"""
  FORMULA ANALITYCZNA dla K1^(NLO):

    K1 = (-c2 + sqrt(c2^2 + 4*c3)) / (2*c3)

  gdzie:
    c2 = (1+alpha)*Phi_2/2 - U_2/2
    c3 = -alpha*Phi_3/2 - 2*U_3/3

  i Phi_n, U_n sa calkami wykladniczymi Phi_n(a), U_n(a).

  Substituujac formuly analityczne (Phi_2 ~ 1/a, Phi_3 ~ 1/(2a^2)):
    c2 ~ (1+alpha)/(2a)
    K1^(LO) ~ 2*a / (1+alpha)

  Dla alpha=8.5445, a=0.040:
    K1^(LO) ~ 2*0.04/(1+8.5445) = 0.08/9.5445 = {2*a/(1+alpha_lep):.6f}
    K1^(true) = {K1_true:.6f}
    """)

    # Asymptotics dla alpha duze:
    print(f"  Asymptotyka dla duzego alpha:")
    print(f"    K1 ~ 2*a/alpha  (dla alpha >> 1)")
    print(f"    K1 ~ 2*{a:.3f}/{alpha_lep:.4f} = {2*a/alpha_lep:.6f}")
    print(f"    K1_true = {K1_true:.6f}  (roznica: {abs(2*a/alpha_lep - K1_true)/K1_true*100:.1f}%)")

    return K1_LO, K1_NLO, K1_NNLO


# =========================================================================
# SEKCJA D: Analityczne K3 — duze K, korekcja do C=sqrt(4.5)
# =========================================================================
def section_D(U4, U6):
    print("\n" + "=" * 72)
    print("SEKCJA D: ANALITYCZNE K3 — DUZE K, SKAD C=2.000?")
    print("=" * 72)
    print(f"""
  Dla duzego K: phi(r) = K*e^{{-r}}/r >> 1 przy r ~ a.

  Czlony dominujace w g(K) dla K >> 1:
    Kinetyczny (czlon alpha/phi): ~ alpha * K * J_alpha(a) / (4*pi*K) = alpha*J_a/4pi
    Potencjal kw.  (phi^4/4):    ~ -K^3 * I4(a) / (4*pi*K) = -K^2*I4/(4pi)
    Potencjal seks. (phi^6*lam/6): ~ lambda*K^5 * I6(a) / (4*pi*K) = lam*K^4*I6/(4pi)

  Warunek g(K3) = 0:
    -K3^2 * I4/4 + lambda*K3^4 * I6/6 = 1 + ...

  Dla LO (pomijamy 1 i czlon kinetyczny):
    K3^2 = 3*I4 / (2*lambda*I6)

  CALKI DOKLADNE (przez Ei):
    I4(a) = int_a^inf e^{{-4r}}/r^2 dr = e^{{-4a}}/a - 4*E1(4a)
    I6(a) = int_a^inf e^{{-6r}}/r^4 dr  [formula wyzej]

  LO asymptotyka dla a->0:
    I4 ~ 1/a,  I6 ~ 1/(3*a^3)
    K3^2 ~ 3*(1/a) / (2*lambda * 1/(3a^3)) = 9*a^2 / (2*lambda)
    K3^(LO-0) = sqrt(4.5) * a/sqrt(lambda) = {np.sqrt(4.5):.4f} * a/sqrt(lambda)

  DOKLADNE (Ei):
    K3^2 = 3*I4_exact / (2*lambda*I6_exact)
    """)

    I4_exact = np.exp(-4*a)/a - 4*E1(4*a)
    I6_exact = (np.exp(-6*a)/(3*a**3)
                - np.exp(-6*a)/a**2
                + 6*(np.exp(-6*a)/a - 6*E1(6*a)))

    I4_naive = 1.0/a
    I6_naive = 1.0/(3*a**3)

    K3_LO0_sq = 9*a**2 / (2*lam_K)
    K3_exact_sq = 3*I4_exact / (2*lam_K*I6_exact)

    K3_LO0   = np.sqrt(K3_LO0_sq)
    K3_exact = np.sqrt(K3_exact_sq)
    K3_true  = 2.000 * a / np.sqrt(lam_K)   # C=2.000, P31

    print(f"  Dla a={a}, lambda={lam_K:.4e}:")
    print(f"    I4(a) naive = 1/a           = {I4_naive:.4f}")
    print(f"    I4(a) exact (Ei)            = {I4_exact:.4f}  (stosunek: {I4_exact/I4_naive:.4f})")
    print(f"    I6(a) naive = 1/(3a^3)      = {I6_naive:.4e}")
    print(f"    I6(a) exact (Ei)            = {I6_exact:.4e}  (stosunek: {I6_exact/I6_naive:.4f})")
    print()
    print(f"    K3^(LO-0,  C=sqrt(4.5))     = {K3_LO0:.5f}  (C_LO0  = {np.sqrt(4.5):.4f})")
    print(f"    K3^(exact-I, C z dokl.I4,I6) = {K3_exact:.5f}  (C_exI  = {K3_exact*np.sqrt(lam_K)/a:.4f})")
    print(f"    K3^(true, C=2.000)            = {K3_true:.5f}  (C_true = 2.0000)")
    print()
    print(f"    Blad K3^(LO-0)  od prawdy: {(K3_LO0 - K3_true)/K3_true*100:+.2f}%")
    print(f"    Blad K3^(exact) od prawdy: {(K3_exact - K3_true)/K3_true*100:+.2f}%")

    # Wyznaczyc C analitycznie ze stanu K3^2 = 3*I4/(2lam*I6)
    C_analytic = np.sqrt(3*I4_exact/(2*I6_exact)) / a
    print(f"""
  FORMULA ANALITYCZNA:
    C = sqrt( 3*I4(a) / (2*I6(a)) ) / a
      = sqrt( 3*{I4_exact:.4f} / (2*{I6_exact:.4e}) ) / {a}
      = {C_analytic:.4f}

  Porownanie:
    C_LO (a->0)  = sqrt(4.5) = {np.sqrt(4.5):.4f}  (bez eksponenntow)
    C_exact(I)   = {C_analytic:.4f}  (z dokladnymi I4, I6)
    C_numeryczna = 2.000           (z P31)
    """)

    # Skad pochodzi roznica miedzy C_analytic a 2.000?
    # Wyznaczyc co jeszcze wplywa na K3 (czlon staly w g=0)
    print(f"  ANALIZA KOREKTY: dlaczego C_exact={C_analytic:.4f} != 2.000?")
    print(f"""
  Rownanie g(K3) = 0 ma DODATKOWE czlony:
    g(K3) = -1 + E_kin(K3)/4piK3 - K3^2*I4/4 + lam*K3^4*I6/6 = 0

  1. Czlon STALY: -1
     Wplyw na K3: jesli wlaczamy -1, K3 rosnie nieco (musi kompensowac -1)
     Estymacja: deltaK3^2 ~ 2/(K3^2*(lam*I6/6 - I4/4*dln/dK)) ... skomplikowane

  2. Czlon KINETYCZNY: alpha * J_alpha / (2*K3)
     Dla K3=34, alpha=8.5, J_alpha ~ 1/a = 25:
     ~ 8.5*25/(2*34) ~ 3.1  (nie zaniedbywalny!)

  3. Czlon KINETYCZNY K^2: (1+alpha)*Phi2/(2*K3)
     ~ 9.5*27/(2*34) ~ 3.8

  Wiec pomin WSZYSTKICH korekt daje blad ~{abs(K3_exact - K3_true)/K3_true*100:.0f}%.
  Z uwzglednieniem czlonu stalego i kinetycznego mozna zblizycdo C=2.000.
    """)

    # Pelne rownanie na K3: NLO z czlonem stalym i kinetycznym
    # g(K3) = -1 + A_kin(K3) - K3^2*I4/4 + lam*K3^4*I6/6 = 0
    # Gdzie A_kin(K3) ~ alpha * J_alpha / (2*K3) [dominujacy czlon kinetyczny dla duzego K]

    # Calka J_alpha
    def J_alpha_integrand(r):
        return np.exp(-r)*(1+r)**2/r

    J_alpha, _ = quad(J_alpha_integrand, a, 50.0)
    A_alpha = alpha_lep * J_alpha / 2    # korekta z czlonu alpha/phi

    print(f"  J_alpha(a={a}) = int_a^inf e^{{-r}}*(1+r)^2/r dr = {J_alpha:.4f}")
    print(f"  A_alpha = alpha*J_alpha/2 = {A_alpha:.4f}")
    print()
    print(f"  Pelne LO rownanie z korektami dla K3:")
    print(f"    -1 + A_alpha/K3 - K3^2*I4/4 + lam*K3^4*I6/6 = 0")
    print()

    # Rozwiaz numerycznie to rownanie
    def g_approx_large_K(K3):
        return -1 + A_alpha/K3 - K3**2*I4_exact/4 + lam_K*K3**4*I6_exact/6

    try:
        K3_NLO = brentq(g_approx_large_K, 20.0, 60.0, xtol=1e-8)
        C_NLO = K3_NLO * np.sqrt(lam_K) / a
        print(f"  K3^(NLO z I4_exact, I6_exact, A_alpha) = {K3_NLO:.5f}  (C_NLO = {C_NLO:.4f})")
        print(f"  Blad od prawdy K3=34.213: {(K3_NLO - K3_true)/K3_true*100:+.2f}%")
    except Exception as e:
        print(f"  Blad NLO: {e}")
        K3_NLO = np.nan
        C_NLO = np.nan

    return K3_exact, K3_NLO, C_analytic


# =========================================================================
# SEKCJA E: Struktura trzech zer — topologia g(K)
# =========================================================================
def section_E(c2, c3, c4, c6):
    print("\n" + "=" * 72)
    print("SEKCJA E: TOPOLOGIA g(K) — DLACZEGO DOKLADNIE TRZY ZERA?")
    print("=" * 72)
    print(f"""
  g(K) = E(K)/(4*pi*K) - 1
       = c2*K + c3*K^2 + c4*K^3 + c6*K^5 - 1 + [wyzsze rzedy]

  Wspolczynniki dla a={a}, alpha={alpha_lep}:
    c2 = {c2:.4f}  > 0  (zawsze, bo Phi2 >> U2/2)
    c3 = {c3:.4f}  < 0  (sila odcinkajaca: alpha*Phi3 dominuje)
    c4 = {c4:.4f}  > 0  (poprawki NLO)
    c6 = {c6:.4e} > 0  (czlon stabilizujacy lambda*U6/6)

  ANALIZA ZNAKOW:
    g(K->0+) = -1  < 0
    g(K~K1) = 0  (pierwsze przejscie przez 0: K1 ~ 1/c2 ~ 2a/(1+alpha))
    g(K > K1): rosnie az do lokalnego maksimum
    g(K_max) > 0  (jesli max > 0)  [warunek istnienia K2]
    g(K2) = 0  (drugie przejscie: suma c3*K^2 + ... < 0)
    g(K > K2): maleje (czlon c3*K^2 dominuje nad c2*K dla K~K2~2)
    g(K ~ K3): c6*K^5 dominuje, rosnie -> trzecie przejscie przez 0
    """)

    # Narysuj g(K) numerycznie i zaznacz zera
    def E_numeric_full(K):
        N = 2000
        t = np.linspace(0, 1, N)
        r = a * (R_MAX/a)**t
        phi = np.maximum(1 + K*np.exp(-r)/r, 1e-10)
        dphi = K*np.exp(-r)*(-r-1)/r**2
        V1 = 1/3 - 1/4
        Ek = np.trapezoid(0.5*dphi**2*(1+alpha_lep/phi)*r**2, r)
        Ep = np.trapezoid(((phi**3/3 - phi**4/4 + lam_K*(phi-1)**6/6) - V1)*r**2, r)
        return Ek + Ep

    def g_numeric(K):
        if K <= 0:
            return np.nan
        return E_numeric_full(K) / (4*np.pi*K) - 1

    # Szybki skan
    K_vals = np.concatenate([
        np.logspace(-4, -2, 20),
        np.linspace(0.01, 0.05, 20),
        np.linspace(0.05, 0.3, 20),
        np.logspace(-0.5, 0, 10),
        np.logspace(0, 1, 20),
        np.logspace(1, np.log10(50), 15)
    ])
    K_vals = np.unique(np.sort(K_vals))

    g_vals = np.array([g_numeric(K) for K in K_vals])

    # Znajdz zera
    sign_changes = np.where(np.diff(np.sign(g_vals)))[0]
    K_zeros = []
    for idx in sign_changes:
        if np.isfinite(g_vals[idx]) and np.isfinite(g_vals[idx+1]):
            try:
                K_z = brentq(g_numeric, K_vals[idx], K_vals[idx+1], xtol=1e-8)
                K_zeros.append(K_z)
            except Exception:
                pass

    print(f"  Zera g(K) numerycznie: {len(K_zeros)}")
    for i, Kz in enumerate(K_zeros):
        print(f"    K{i+1} = {Kz:.6f}")

    # Klucz: warunek istnienia K2 (g_max > 0)
    # Szukamy maksimum g(K) miedzy K1 i K3
    if len(K_zeros) >= 2:
        K1_num = K_zeros[0]
        K3_num = K_zeros[-1]
        # Maksimum g miedzy K1 i K2
        g_between = [(K, g_numeric(K)) for K in np.logspace(np.log10(K1_num), 0, 30)]
        K_max_pair = max(g_between, key=lambda x: x[1] if np.isfinite(x[1]) else -np.inf)
        g_max = K_max_pair[1]
        K_at_gmax = K_max_pair[0]

        print(f"\n  Maximum g(K) miedzy K1 i K2:")
        print(f"    K_at_gmax ~ {K_at_gmax:.4f},  g_max ~ {g_max:.4f}")
        print(f"""
  TWIERDZENIE ISTNIENIA (P47):
    Niech g(K) = c2*K + c3*K^2 + ... - 1 dla potencjalu V_mod.
    Wtedy g(K) ma DOKLADNIE TRZY ZERA jesli i tylko jesli:

    (i)  c2 > 0  (zawsze prawda dla alpha > 0, a > 0)   ← K1 istnieje
    (ii) g ma lokalne maksimum g_max > 0 przy K_max      ← K2 istnieje
    (iii) g ma lokalne minimum g_min < 0 przy K_min > K_max  ← K3 > K2 istnieje

    Warunek (iii) z pojawieniem sie czlonu c6*K^5 > 0 (lambda > 0):
      g(K->inf) -> +inf  =>  K3 zawsze istnieje gdy lambda > 0.

    Warunek (ii) zalezy od ZNAKU lokalnego maksimum.
    Z P38/P39: warunek bifurkacji alpha > alpha_c(a) gwarantuje g_max > 0.

  FIZYCZNA INTERPRETACJA:
    K1: rownowaga kinetyczna (Ek ~ 4pi*K) — "najslabszy" soliton
    K2: rownowaga EM (czlon alpha/phi nasica Ek) — "sredni" soliton
    K3: rownowaga seksyczna (V_sex ~ lambda*K^6 ~ 4pi*K) — "cieski" soliton
        """)

    # Zbadaj warunek alpha_c
    print(f"\n  Warunek na trzy zera: alpha > alpha_c(a)")
    print(f"  Szukamy alpha_c dla a={a}:")
    # Chodzi o to ze przy alpha=alpha_c, g_max = 0 (K2 = K1 bifurkacja)
    def g_max_value(alpha_test):
        """Wartosc maksimum g(K) dla danego alpha."""
        def g_test(K):
            N = 1000
            t = np.linspace(0, 1, N)
            r = a * (R_MAX/a)**t
            phi = np.maximum(1 + K*np.exp(-r)/r, 1e-10)
            dphi = K*np.exp(-r)*(-r-1)/r**2
            V1 = 1/3 - 1/4
            Ek = np.trapezoid(0.5*dphi**2*(1+alpha_test/phi)*r**2, r)
            Ep = np.trapezoid(((phi**3/3 - phi**4/4 + lam_K*(phi-1)**6/6) - V1)*r**2, r)
            return (Ek+Ep)/(4*np.pi*K) - 1

        # Maksimum g miedzy K=0.01 i K=5
        K_scan = np.linspace(0.01, 3.0, 80)
        g_scan = np.array([g_test(K) for K in K_scan])
        return np.max(g_scan[np.isfinite(g_scan)])

    # Gruba bifurkacja (szybko)
    alpha_c_approx = None
    for alpha_test in [0.1, 0.5, 1.0, 1.5, 2.0, 3.0]:
        gm = g_max_value(alpha_test)
        if gm > 0:
            alpha_c_approx = alpha_test
            break

    if alpha_c_approx is not None:
        print(f"  Pierwsze alpha z g_max>0: alpha ~ {alpha_c_approx}")
        print(f"  (Dokladna bifurkacja z P39: alpha_c ~ 1.1–1.2 dla a=0.040)")

    return K_zeros


# =========================================================================
# SEKCJA F: Analityczne r21(alpha) — skad r21=207?
# =========================================================================
def section_F(K1_NLO):
    print("\n" + "=" * 72)
    print("SEKCJA F: ANALITYCZNE r21(alpha) — SKAD r21=207 DLA LEPTONOW?")
    print("=" * 72)
    print(f"""
  r21 = K2/K1

  K1 jest analityczne (Sekcja C): K1 ~ 2a/(1+alpha) dla duzego alpha.
  K2 ~ const(alpha, a) ~ niezalezne od alpha dla duzych alpha (P31/P37).

  Asymptotics dla duzego alpha:
    K1 ~ 2a/alpha  =>  r21 = K2/K1 ~ K2*alpha/(2a)

  Dla K2 ~ 2.034 (P40), a=0.040:
    r21 ~ 2.034 * alpha / (2*0.040) = 25.43 * alpha

  Dla r21 = 207: alpha ~ 207/25.43 = {207/25.43:.2f}
  Rzeczywiste alpha_Koide = 8.5445

  Skad roznica? Asymptotyka K1~2a/alpha jest LO; K1^(NLO) z sekcji C:
    K1^(NLO) = {K1_NLO:.6f}
  =>  K2/K1^(NLO) ~ {2.0344/K1_NLO:.1f}   (powinno byc 207)
    """)

    # Oblicz r21 dla roznych alpha z formuły analitycznej
    print(f"  Skalowanie r21(alpha) z LO i NLO formul:")
    print(f"  {'alpha':>8}  {'K1_LO':>10}  {'r21_LO':>8}  {'K1_NLO':>10}  {'r21_NLO':>10}  {'K1_num':>10}  {'r21_num':>8}")
    print(f"  {'-'*72}")

    # c2, c3 zaleza od alpha — przelicz dla roznych alpha
    for alpha_t in [0.22, 1.0, 2.0, 5.0, 8.5445, 10.0, 20.343]:
        # Calki dla danego alpha
        Phi2 = Phi_n_numeric(2, a)  # niezalezne od alpha
        Phi3 = Phi_n_numeric(3, a)
        U2   = U_n_numeric(2, a)
        U3   = U_n_numeric(3, a)

        c2_t = (1+alpha_t)*Phi2/2 - U2/2
        c3_t = -alpha_t*Phi3/2 - 2*U3/3

        K1_LO_t  = 1.0/c2_t if c2_t > 0 else np.nan
        disc = c2_t**2 + 4*c3_t
        K1_NLO_t = (-c2_t + np.sqrt(max(disc,0))) / (2*c3_t) if disc > 0 and c3_t != 0 else np.nan

        # K2 ~ const (tu uzywamy wartosci z P40 jako ref)
        K2_ref = 2.0344  # dla leptonow

        r21_LO_t  = K2_ref/K1_LO_t  if not np.isnan(K1_LO_t) else np.nan
        r21_NLO_t = K2_ref/K1_NLO_t if not np.isnan(K1_NLO_t) else np.nan

        # Numeryczne K1 dla danego alpha (szybkie szukanie)
        def g_t(K):
            N = 1500
            t = np.linspace(0, 1, N)
            r = a * (R_MAX/a)**t
            phi = np.maximum(1 + K*np.exp(-r)/r, 1e-10)
            dphi = K*np.exp(-r)*(-r-1)/r**2
            V1 = 1/3 - 1/4
            Ek = np.trapezoid(0.5*dphi**2*(1+alpha_t/phi)*r**2, r)
            Ep = np.trapezoid(((phi**3/3 - phi**4/4 + lam_K*(phi-1)**6/6) - V1)*r**2, r)
            return (Ek+Ep)/(4*np.pi*K) - 1

        try:
            K1_num_t = brentq(g_t, 1e-5, 0.5, xtol=1e-8)
            r21_num_t = K2_ref / K1_num_t
        except Exception:
            K1_num_t = np.nan
            r21_num_t = np.nan

        print(f"  {alpha_t:>8.4f}  {K1_LO_t:>10.6f}  {r21_LO_t:>8.1f}  "
              f"{K1_NLO_t:>10.6f}  {r21_NLO_t:>10.1f}  "
              f"{K1_num_t:>10.6f}  {r21_num_t:>8.1f}")

    print(f"""
  WNIOSEK SEKCJI F:
    K1^(NLO) daje r21 ~ K2/K1^(NLO)  — blad {abs(2.0344/K1_NLO - 207)/207*100:.0f}% od 207.
    K1^(NNLO) z 3-wyrazowym szeregiem daje lepsze przyblizenie.

    Asymptotics: r21 ~ K2 * (1+alpha) / (2 * Phi2) gdzie Phi2 ~ 1/a + log terms.
    Dla r21=207 nie ma specjalnej zasady — r21 wyznaczone przez alpha_Koide=8.5445,
    ktore samo pochodzi z dopasowania do mas leptonow (P28/P29).

    PYTANIE OTWARTE: Czy istnieje zasada wyzszego rzedu wiazaca alpha_Koide=8.5445?
    (Teoria grup? Topologia? Kwantowanie topologiczne?)
    """)


# =========================================================================
# SEKCJA G: Podsumowanie i wyniki analityczne
# =========================================================================
def section_G(K1_LO, K1_NLO, K1_NNLO, K3_exact, K3_NLO, C_analytic, K_zeros):
    print("\n" + "=" * 72)
    print("SEKCJA G: PODSUMOWANIE P47 — WYPROWADZENIA Z PIERWSZYCH ZASAD")
    print("=" * 72)

    K1_true = 0.009839
    K3_true = 2.000 * a / np.sqrt(lam_K)
    K2_true = 2.0344

    print(f"""
  ╔══════════════════════════════════════════════════════════════════════╗
  ║      WYNIKI P47: ANALITYCZNE CALKI ENERGII TGP                     ║
  ╚══════════════════════════════════════════════════════════════════════╝

  PUNKT WYJSCIA: E[K] = 4pi*int [(dphi/dr)^2/2*(1+alpha/phi) + V_mod-V_mod(1)] r^2 dr

  ROZWINICIECZE E(K) = 4pi * [c2*K^2 + c3*K^3 + c4*K^4 + c6*K^6 + ...]

  1) ANALITYCZNE K1 (male-K):
     K1^(LO)   = 1/c2 = 2a/[(1+alpha)*Phi2 - U2]  ~  2a/(1+alpha)
     K1^(NLO)  = z rownania kwadratowego c2*K + c3*K^2 = 1
               = [-c2 + sqrt(c2^2 + 4*c3)] / (2*c3)
     Wyniki dla alpha={alpha_lep}, a={a}:
       K1^(LO)   = {K1_LO:.6f}  (blad: {(K1_LO-K1_true)/K1_true*100:+.1f}%)
       K1^(NLO)  = {K1_NLO:.6f}  (blad: {(K1_NLO-K1_true)/K1_true*100:+.1f}%)
       K1^(NNLO) = {K1_NNLO:.6f}  (blad: {(K1_NNLO-K1_true)/K1_true*100:+.1f}%)
       K1_true   = {K1_true:.6f}

  2) ANALITYCZNE K3 (duze-K, przez I4(a), I6(a)):
     K3^2 = 3*I4(a) / (2*lambda*I6(a))
     I4(a) = e^{{-4a}}/a - 4*E1(4a)
     I6(a) = e^{{-6a}}/(3a^3) - e^{{-6a}}/a^2 + 6*(e^{{-6a}}/a - 6*E1(6a))
     C = sqrt(3*I4/(2*I6)) / a
     Wyniki:
       K3^(LO, a->0, C=sqrt(4.5)) = {np.sqrt(4.5)*a/np.sqrt(lam_K):.5f}  (C={np.sqrt(4.5):.4f})
       K3^(exact-I, z Ei)          = {K3_exact:.5f}  (C={C_analytic:.4f})
       K3^(NLO, +czl.kin)          = {K3_NLO:.5f}  (C={K3_NLO*np.sqrt(lam_K)/a:.4f})
       K3_true (P31)               = {K3_true:.5f}  (C=2.0000)

  3) TOPOLOGIA g(K) — TRZY ZERA:
     Warunki konieczne i wystarczajace:
       (i)  c2 > 0  <==> alpha > 0  [zawsze spelone]
       (ii) g_max > 0 <==> alpha > alpha_c(a)  [bifurkacja P38/P39]
       (iii) lambda > 0  =>  g(K->inf) = +inf  [K3 zawsze istnieje]
     Zera numeryczne: K1={K_zeros[0]:.5f}, K2={K_zeros[1]:.4f}, K3={K_zeros[2]:.4f}

  4) ROZWINICIECZE r21(alpha):
     K1 ~ 2a/(1+alpha) * [1 + corrections]  =>  r21 = K2/K1 ~ K2*(1+alpha)/(2a)
     Ale K2 ~ const(alpha) ≈ 2.034 (slab. zal. od alpha)
     => r21 ~ 2.034*(1+alpha)/(2*0.040) = 25.43*(1+alpha)

  STAN WIEDZY (P47):
    K1^(NLO):  {(K1_NLO-K1_true)/K1_true*100:+.1f}% od prawdy  [wzor zamkniety]
    K3^(exact):  {(K3_exact-K3_true)/K3_true*100:+.1f}% od prawdy  [wzor Ei]
    K3^(NLO):  {(K3_NLO-K3_true)/K3_true*100:+.1f}% od prawdy  [+czl.kin]

  BRAKI ANALITYCZNE (na P48):
    - K2 nie ma wzoru analitycznego (pochodzi z lokalnego maksimum dE/dK)
    - Korrekta K1 do <1% wymaga resummacji szeregu lub lepszego profilu
    - Wartosci parametrow alpha_Koide=8.5445 nie sa wyprowadzone z pierwszych zasad
    """)


# =========================================================================
# MAIN
# =========================================================================
if __name__ == '__main__':
    print("=" * 72)
    print("P47: ANALITYCZNE CALKI ENERGII TGP — WYPROWADZENIE Z PIERWSZYCH ZASAD")
    print("     E[K] = 4pi*int [(dphi/dr)^2*(1+alpha/phi)/2 + V_mod-V1] r^2 dr")
    print("=" * 72)

    section_A()
    c2, c3, c4, c6, Phi2, Phi3, Phi4, U2, U3, U4, U6 = section_B()
    K1_LO, K1_NLO, K1_NNLO = section_C(c2, c3, c4)
    K3_exact, K3_NLO, C_analytic = section_D(U4, U6)
    K_zeros = section_E(c2, c3, c4, c6)
    section_F(K1_NLO)
    section_G(K1_LO, K1_NLO, K1_NNLO, K3_exact, K3_NLO, C_analytic, K_zeros)

    print("\n" + "=" * 72)
    print("P47 ZAKONCZONE.")
    print("Calki E(K) analitycznie przez Ei; K1^NLO z wzoru kwadr.; K3^NLO z I4,I6,J_alpha.")
    print("=" * 72)
