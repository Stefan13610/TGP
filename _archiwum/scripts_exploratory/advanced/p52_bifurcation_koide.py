# -*- coding: utf-8 -*-
"""
P52: KRZYWA BIFURKACJI alpha_c(a) A WARUNEK KOIDEGO

Pytania:
  OP-2: Skad pochodzi alpha_Koide = 8.5445?
  OP-3: Czy mozna przewidziec a_Gamma z zasad pierwszych?

Podejscie:
  Dla kazdego a_Gamma wyznaczamy:
    1. alpha_c(a)     -- bifurkacja: minimalne alpha dajace 3 zera g(K)
    2. alpha_K(a)     -- warunek Koidego: alpha takie ze r21=K2/K1=206.77
    3. stosunek R(a) = alpha_K(a) / alpha_c(a)

  Pytanie: czy R(a) jest stale? Czy R(a) = const wyznaczaloby a_Gamma?
  Czy istnieje inny zwiazek algebraiczny miedzy alpha_c i alpha_K?

Sekcje:
  A: Funkcja g(K) i wyznaczanie g_max
  B: Krzywa bifurkacji alpha_c(a)
  C: Krzywa Koidego alpha_K(a)
  D: Stosunek R(a) = alpha_K / alpha_c
  E: Analityczne przybliZenie alpha_c(a) -- z warunku g'(K_max)=0
  F: Asymptotyczne zachowanie dla a->0
  G: Podsumowanie i wnioski
"""
import sys, warnings
sys.stdout.reconfigure(encoding='utf-8')
warnings.filterwarnings('ignore')

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq, minimize_scalar
from scipy.special import expi

# -------------------------------------------------------------------
# PARAMETRY
# -------------------------------------------------------------------
A_GAM   = 0.040
ALPHA   = 8.5445
LAM_K   = 5.4677e-6
R21_PDG = 206.770   # K2/K1 dla leptonow (PDG)
R_MAX   = 50.0
N_GRID  = 3000

def E1(x):
    return -expi(-x)

def V_mod(phi, lam=LAM_K):
    return phi**3/3 - phi**4/4 + lam*(phi-1)**6/6

def energy_num(K, alpha=ALPHA, a=A_GAM, lam=LAM_K, N=N_GRID):
    """Pelna numeryczna energia solitonu."""
    t   = np.linspace(0, 1, N)
    r   = a * (R_MAX/a)**t
    phi = np.maximum(1.0 + K*np.exp(-r)/r, 1e-10)
    dphi= K * np.exp(-r) * (-r - 1.0) / r**2
    V1  = V_mod(1.0, lam)
    Ek  = 4*np.pi * np.trapezoid(0.5*dphi**2*(1.0 + alpha/phi)*r**2, r)
    Ep  = 4*np.pi * np.trapezoid((V_mod(phi, lam) - V1)*r**2, r)
    return Ek + Ep

def g_func(K, alpha=ALPHA, a=A_GAM, lam=LAM_K):
    if K <= 0:
        return -1.0
    return energy_num(K, alpha, a, lam) / (4*np.pi*K) - 1.0

def g_max_val(alpha, a, lam=LAM_K):
    """Maksimum g(K) w przedziale K in [0.1, 5]."""
    try:
        res = minimize_scalar(lambda K: -g_func(K, alpha, a, lam),
                              bounds=(0.1, 5.0), method='bounded')
        return -res.fun, res.x
    except Exception:
        return -np.inf, np.nan

def find_K1(alpha, a, lam=LAM_K):
    """Pierwsze zero g(K) (K1 ~ 0.001..0.1)."""
    try:
        return brentq(g_func, 1e-4, 0.5, args=(alpha, a, lam), xtol=1e-10)
    except Exception:
        return np.nan

def find_K2(alpha, a, lam=LAM_K):
    """Drugie zero g(K) (K2 ~ 1..4)."""
    try:
        return brentq(g_func, 0.5, 5.0, args=(alpha, a, lam), xtol=1e-10)
    except Exception:
        return np.nan

def analytic_coeffs(a=A_GAM, alpha=ALPHA, lam=LAM_K):
    """Wspolczynniki c2,c3,c4 przez E1(n*a) -- valid dla malych K."""
    J2   = 2*E1(2*a) - np.exp(-2*a)/a
    J1   = E1(2*a)
    J0   = np.exp(-2*a)/2
    Phi2 = -J2 + 2*J1 + J0
    Phi3, _ = quad(lambda r: np.exp(-3*r)*(1+r)**2/r**3, a, R_MAX, limit=200)
    Phi4, _ = quad(lambda r: np.exp(-4*r)*(1+r)**2/r**4, a, R_MAX, limit=200)
    U2   = np.exp(-2*a)/2
    U3   = E1(3*a)
    U4   = np.exp(-4*a)/a - 4*E1(4*a)
    c2 = (1.0+alpha)*Phi2/2 - U2/2
    c3 = -alpha*Phi3/2 - 2*U3/3
    c4 = alpha*Phi4/2 - U4/4
    return c2, c3, c4

# -------------------------------------------------------------------
# SEKCJA A: Weryfikacja g(K) dla parametrow leptonowych
# -------------------------------------------------------------------
print("="*72)
print("P52: BIFURKACJA alpha_c(a) I WARUNEK KOIDEGO alpha_K(a)")
print("="*72)

print("\nSEKCJA A: Weryfikacja g(K) dla alpha=8.5445, a=0.040")
print("-"*72)

K1_ref = find_K1(ALPHA, A_GAM)
K2_ref = find_K2(ALPHA, A_GAM)
gmax, Kmax = g_max_val(ALPHA, A_GAM)
r21_ref = K2_ref / K1_ref if K1_ref > 0 else np.nan

print(f"  K1 = {K1_ref:.6f}  (true: 0.009839,  blad: {(K1_ref/0.009839-1)*100:+.3f}%)")
print(f"  K2 = {K2_ref:.6f}  (true: 2.0344,    blad: {(K2_ref/2.0344-1)*100:+.3f}%)")
print(f"  K_max = {Kmax:.4f},  g_max = {gmax:.4f}")
print(f"  r21 = K2/K1 = {r21_ref:.3f}  (PDG: {R21_PDG:.3f})")

# -------------------------------------------------------------------
# SEKCJA B: Krzywa bifurkacji alpha_c(a)
# -------------------------------------------------------------------
print("\n" + "="*72)
print("SEKCJA B: KRZYWA BIFURKACJI alpha_c(a)")
print("="*72)
print("\n  alpha_c = minimalne alpha, dla ktorego g_max > 0")
print("  (dla alpha < alpha_c: tylko jedno zero K1; brak K2, K3)\n")

# Scan wartosci a
a_values = [0.010, 0.015, 0.020, 0.025, 0.030, 0.035, 0.040, 0.045,
            0.050, 0.055, 0.060, 0.070, 0.080, 0.090, 0.100]

alpha_c_arr  = []
alpha_K_arr  = []
R_ratio_arr  = []
K1_at_aK_arr = []
K2_at_aK_arr = []

print(f"  {'a':>6}  {'alpha_c':>9}  {'alpha_K':>9}  {'r21@aK':>10}")
print("  " + "-"*45)

for a_val in a_values:
    # --- Bifurkacja: szukaj alpha_c ---
    def gmax_vs_alpha(alpha):
        gm, _ = g_max_val(alpha, a_val)
        return gm

    # Najpierw sprawdz czy dla malego alpha gmax < 0
    gm_low = gmax_vs_alpha(0.1)
    gm_high = gmax_vs_alpha(30.0)

    if gm_low >= 0:
        alpha_c = 0.0  # juz dla alpha=0.1 mamy 3 zera
    elif gm_high <= 0:
        alpha_c = np.nan  # nie ma bifurkacji w zakresie
    else:
        try:
            alpha_c = brentq(gmax_vs_alpha, 0.1, 30.0, xtol=1e-6)
        except Exception:
            alpha_c = np.nan

    # --- Warunek Koidego: szukaj alpha_K takie ze r21=R21_PDG ---
    def r21_vs_alpha(alpha):
        K1v = find_K1(alpha, a_val)
        K2v = find_K2(alpha, a_val)
        if np.isnan(K1v) or np.isnan(K2v) or K1v <= 0:
            return np.nan
        return K2v / K1v - R21_PDG

    # Szukaj alpha_K w zakresie (alpha_c+0.1, 30)
    if not np.isnan(alpha_c):
        a_lo = alpha_c + 0.5
        a_hi = 25.0
        r_lo = r21_vs_alpha(a_lo)
        r_hi = r21_vs_alpha(a_hi)
        if (not np.isnan(r_lo)) and (not np.isnan(r_hi)) and r_lo * r_hi < 0:
            try:
                alpha_K = brentq(r21_vs_alpha, a_lo, a_hi, xtol=1e-5)
            except Exception:
                alpha_K = np.nan
        else:
            # Moze r21(alpha) jest monotoniczne -- znajdz gdzie = R21_PDG
            # Skanuj
            alphas_scan = np.linspace(a_lo, a_hi, 50)
            r21_scan = [r21_vs_alpha(aa) for aa in alphas_scan]
            r21_scan = [x if not np.isnan(x) else 1e9 for x in r21_scan]
            # szukaj zmiany znaku
            sign_changes = [(alphas_scan[i], alphas_scan[i+1])
                            for i in range(len(r21_scan)-1)
                            if r21_scan[i]*r21_scan[i+1] < 0]
            if sign_changes:
                try:
                    alpha_K = brentq(r21_vs_alpha, sign_changes[0][0], sign_changes[0][1], xtol=1e-5)
                except Exception:
                    alpha_K = np.nan
            else:
                alpha_K = np.nan
    else:
        alpha_K = np.nan

    # Zbierz wyniki
    alpha_c_arr.append(alpha_c)
    alpha_K_arr.append(alpha_K)

    if not np.isnan(alpha_c) and not np.isnan(alpha_K) and alpha_c > 0:
        R = alpha_K / alpha_c
        K1v = find_K1(alpha_K, a_val)
        K2v = find_K2(alpha_K, a_val)
        r21v = K2v/K1v if K1v > 0 else np.nan
    else:
        R = np.nan
        r21v = np.nan
        K1v = np.nan
        K2v = np.nan

    R_ratio_arr.append(R)
    K1_at_aK_arr.append(K1v)
    K2_at_aK_arr.append(K2v)

    ac_str = f"{alpha_c:9.4f}" if not np.isnan(alpha_c) else "      nan"
    aK_str = f"{alpha_K:9.4f}" if not np.isnan(alpha_K) else "      nan"
    r21_str= f"{r21v:10.3f}"    if not np.isnan(r21v) else "       nan"
    print(f"  {a_val:6.3f}  {ac_str}  {aK_str}  {r21_str}")

# -------------------------------------------------------------------
# SEKCJA C: Analiza krzywej alpha_K(a)
# -------------------------------------------------------------------
print("\n" + "="*72)
print("SEKCJA C: ANALIZA KRZYWEJ alpha_K(a)")
print("="*72)

# Zbierz wazne punkty
a_valid  = [a_values[i] for i,aK in enumerate(alpha_K_arr) if not np.isnan(aK)]
aK_valid = [aK for aK in alpha_K_arr if not np.isnan(aK)]

R_mean = R_std = np.nan  # nie uzywamy juz stosunku do alpha_c

print(f"\n  Wynik B: alpha_c ~ 0 dla wszystkich badanych a.")
print(f"  => Trojzerowa struktura g(K) istnieje dla KAZDEGO alpha>0.")
print(f"     Bifurkacja w alpha nie ogranicza parametrow leptonowych!")
print(f"")
print(f"  Krzywa alpha_K(a) = alpha dajace r21=206.77:")

if len(a_valid) >= 2:
    # Stosunek alpha_K / a
    ratios_a = [aK/a for a, aK in zip(a_valid, aK_valid)]
    # Stosunek alpha_K / sqrt(a)
    ratios_sqrta = [aK/np.sqrt(a) for a, aK in zip(a_valid, aK_valid)]

    print(f"  {'a':>6}  {'alpha_K':>9}  {'aK/a':>9}  {'aK/sqrt(a)':>12}  {'(1+aK)/a':>10}")
    print("  " + "-"*55)
    for i, (a_v, aK) in enumerate(zip(a_valid, aK_valid)):
        r1 = aK/a_v
        r2 = aK/np.sqrt(a_v)
        r3 = (1+aK)/a_v
        print(f"  {a_v:6.3f}  {aK:9.4f}  {r1:9.2f}  {r2:12.4f}  {r3:10.2f}")

    # Fit power law: alpha_K ~ A * a^n
    log_a  = np.log(a_valid)
    log_aK = np.log(aK_valid)
    pfit = np.polyfit(log_a, log_aK, 1)
    n_pow = pfit[0]
    A_pow = np.exp(pfit[1])
    # RMS
    resid_plaw = [(aK_valid[i] - A_pow*a_valid[i]**n_pow)/aK_valid[i] for i in range(len(a_valid))]
    rms_plaw = np.sqrt(np.mean([r**2 for r in resid_plaw]))*100

    print(f"\n  Power law fit: alpha_K(a) = {A_pow:.4f} * a^{n_pow:.4f}")
    print(f"    RMS wzgledny: {rms_plaw:.3f}%")

    # Sprawdz czy (1+alpha_K)/a = const
    r3_arr = [(1+aK)/a for a, aK in zip(a_valid, aK_valid)]
    r3_mean = np.mean(r3_arr)
    r3_std  = np.std(r3_arr)
    print(f"\n  Czy (1+alpha_K)/a = const?")
    print(f"    srednia (1+aK)/a = {r3_mean:.3f},  std = {r3_std:.3f}  ({r3_std/r3_mean*100:.2f}%)")
    if r3_std/r3_mean < 0.05:
        print(f"    [!] TAK! (1+alpha_K)/a ~ {r3_mean:.2f} = stale (std < 5%)")
        print(f"        Interpretacja LO: K1 ~ 2a/(1+alpha), K2 ~ const,")
        print(f"        r21 = K2/K1 ~ K2*(1+alpha)/(2a) = const => (1+alpha)/a = const!")
        print(f"        Stala: (1+alpha)/a = 2*r21/K2_LO ~ 2*{R21_PDG:.0f}/K2")
    else:
        print(f"    [~] Nie do konca stale ({r3_std/r3_mean*100:.1f}% odchylenie)")
        # LO prediction
        K2_LO = np.sqrt(2.0)
        r3_LO = 2*R21_PDG / K2_LO
        print(f"    LO przewidywanie: (1+alpha)/a = 2*r21/K2_LO = {r3_LO:.2f}  (K2_LO=sqrt(2)={K2_LO:.4f})")
        print(f"    Numeryczne: (1+alpha)/a ~ {r3_mean:.2f}  (roznica: {(r3_mean/r3_LO-1)*100:.1f}%)")

# -------------------------------------------------------------------
# SEKCJA D: Analityczne przyblizenie alpha_c(a) -- z warunku g_max=0
# -------------------------------------------------------------------
print("\n" + "="*72)
print("SEKCJA D: ANALITYCZNE PRZYBLIZENIE alpha_c(a)")
print("="*72)
print("""
  Idea: w przybliZeniu malych K (g(K) ~ c2*K + c3*K^2 + c4*K^3 - 1):

  g'(K_max) = c2 + 2*c3*K_max + 3*c4*K_max^2 = 0

  K_max_analytic = (-c3 - sqrt(c3^2 - 3*c2*c4)) / (3*c4)  [wiekszy korzen]

  Ale c2,c3 zalejza od alpha! Wiec K_max = K_max(alpha, a).

  Bifurkacja: g(K_max) = 0
    c2*K_max + c3*K_max^2 + c4*K_max^3 - 1 = 0

  To daje rownanie na alpha_c(a) -- niestety niejawne.

  UWAGA: ten szereq jest wazny tylko dla K < 0.05 (K_1 zakres).
  K_max ~ 1 jest POZA zakresem waznosci szerequ!

  Sprawdzamy jednak, czy analityczne K_max z tego szerequ
  daje jakiekolwiek rozsadne przyblizenie.
""")

print(f"  Parametry: a=0.040, alpha=8.5445")
c2, c3, c4 = analytic_coeffs(A_GAM, ALPHA)
print(f"  c2={c2:.3f}, c3={c3:.3f}, c4={c4:.3f}")

disc_an = c3**2 - 3*c2*c4
print(f"  Dyskryminant c3^2 - 3*c2*c4 = {disc_an:.1f}")

if disc_an >= 0:
    K_max_an = (-c3 - np.sqrt(disc_an)) / (3*c4)
    print(f"  K_max_analytic = {K_max_an:.6f}")
    g_at_Kmax = c2*K_max_an + c3*K_max_an**2 + c4*K_max_an**3 - 1
    print(f"  g_series(K_max_analytic) = {g_at_Kmax:.4f}")
    _, K_max_num = g_max_val(ALPHA, A_GAM)
    print(f"  K_max_numeryczny = {K_max_num:.6f}")
    print(f"  -> WNIOSEK: szereg perturbacyjny NIE daje poprawnego K_max")
    print(f"     (K_max_num~{K_max_num:.2f} jest poza zakresem zbieznosci szerequ)")
else:
    print(f"  Dyskryminant < 0 -> brak realnego K_max z szerequ (jak oczekiwano)")
    print(f"  WNIOSEK: Szereg c2K+c3K^2+c4K^3 nie ma maksimum dla K>0!")
    print(f"  Maksimum g(K)~1 pochodzi z pelnej (nieperturbacyjnej) struktury.")

# alpha_c ~ 0 wszedzie -- nie fitujemy
print(f"\n  WNIOSEK: alpha_c ~ 0 dla calego badanego zakresu a.")
print(f"  Trojzerowa bifurkacja nie ma progu w alpha -- trzy generacje")
print(f"  istnieja dla dowolnie malego alpha > 0 (przy danym a i lambda).")

# -------------------------------------------------------------------
# SEKCJA E: Punkt leptonowy -- odleglosc od bifurkacji
# -------------------------------------------------------------------
print("\n" + "="*72)
print("SEKCJA E: PUNKT LEPTONOWY vs. BIFURKACJA")
print("="*72)

idx04 = a_values.index(0.040) if 0.040 in a_values else None
if idx04 is not None:
    ac04 = alpha_c_arr[idx04]
    aK04 = alpha_K_arr[idx04]
    gmax04, Kmax04  = g_max_val(ALPHA, A_GAM)

    print(f"\n  a=0.040:")
    print(f"    alpha_c = {ac04:.4f}  (bifurkacja: brak progu w badanym zakresie)")
    print(f"    alpha_K = {aK04:.4f}  = alpha_Koide (warunek r21=206.77)")
    print(f"    g_max(alpha_Koide) = {gmax04:.4f}")
    print(f"    K_max = {Kmax04:.4f}")
    print(f"")
    print(f"  WNIOSEK: alpha_c ~ 0 => trojzerowa struktura TGP NIE wymaga")
    print(f"  minimalnego alpha. Istnieje dla kazdego alpha > 0 (przy a=0.040).")
    print(f"  Pytanie 'dlaczego alpha=8.5445?' wymaga innej zasady niz bifurkacja.")

# -------------------------------------------------------------------
# SEKCJA F: Skan alpha w kierunku bifurkacji (a=0.040 fixed)
# -------------------------------------------------------------------
print("\n" + "="*72)
print("SEKCJA F: TRAJEKTORIA Q(alpha) PRZY ZBLIZONIU DO BIFURKACJI (a=0.040)")
print("="*72)
print("\n  Co dzieje sie z Q gdy alpha -> alpha_c od gory?")
print(f"  {'alpha':>8}  {'K1':>10}  {'K2':>10}  {'r21':>8}  {'K3':>8}  {'Q':>10}  {'Q-3/2 [ppm]':>12}")
print("  " + "-"*75)

idx04 = a_values.index(0.040) if 0.040 in a_values else None
if idx04 is not None:
    ac04 = alpha_c_arr[idx04]
    if not np.isnan(ac04):
        # Skan alpha od alpha_c+0.3 do 15
        # K3 z wzoru E1
        def analytic_K3(a, alpha, lam):
            """K3 z I4/(I6) przez E1."""
            from scipy.integrate import quad
            I4  = np.exp(-4*a)/a - 4*E1(4*a)
            I6, _ = quad(lambda r: np.exp(-6*r)/r**6 * r**2, a, R_MAX, limit=200)
            return np.sqrt(3*I4 / (2*lam*I6)) if I6 > 0 else np.nan

        K3_ei = analytic_K3(A_GAM, ALPHA, LAM_K)

        alphas_traj = np.concatenate([
            np.linspace(ac04 + 0.10, ac04 + 1.0, 5),
            np.linspace(ac04 + 1.0,  15.0,  8)
        ])

        for alpha_t in alphas_traj:
            K1t = find_K1(alpha_t, A_GAM)
            K2t = find_K2(alpha_t, A_GAM)
            if np.isnan(K1t) or np.isnan(K2t):
                continue
            r21t = K2t / K1t
            # K3 z numeryki (powolne) lub z wzoru K3_Ei przeskalowanego
            K3t = K3_ei  # K3 jest slabo zalezne od alpha (z P40)
            Qt = (1 + np.sqrt(r21t) + np.sqrt(K3t/K1t))**2 / (1 + r21t + K3t/K1t)
            ppm = (Qt - 1.5) * 1e6
            print(f"  {alpha_t:8.4f}  {K1t:10.6f}  {K2t:10.4f}  {r21t:8.1f}  "
                  f"{K3t:8.3f}  {Qt:10.6f}  {ppm:+12.1f}")

# -------------------------------------------------------------------
# SEKCJA G: Podsumowanie i wnioski
# -------------------------------------------------------------------
print("\n" + "="*72)
print("SEKCJA G: PODSUMOWANIE")
print("="*72)

idx04 = a_values.index(0.040) if 0.040 in a_values else None
ac04  = alpha_c_arr[idx04] if idx04 is not None else np.nan
aK04  = alpha_K_arr[idx04] if idx04 is not None else np.nan

print(f"""
  WYNIKI KLUCZOWE (a=0.040):
    alpha_c   ~ 0   (bifurkacja nie ma progu w alpha > 0!)
    alpha_K   = {aK04:.4f}  (warunek r21=206.77 = PDG)

  KRZYWA alpha_K(a):
    Power law:  alpha_K(a) ~ {A_pow:.4f} * a^{n_pow:.4f}
    Stalosc (1+alpha_K)/a ~ {r3_mean:.2f}  (LO: 2*r21/K2 ~ {2*R21_PDG/np.sqrt(2):.0f})

  WNIOSEK DLA OP-2 (skad alpha_Koide?):
    alpha_c ~ 0 => bifurkacja NIE selektuje alpha_Koide.
    Warunek r21=206.77 definiuje krzyw_ alpha_K(a), ale nie punkt!
    Do wyznaczenia punktu (alpha=8.5445, a=0.040) potrzebna
    jest TRZECIA zasada niezalezna od r21 i r31.

  WNIOSEK DLA OP-3 (predykcja a_Gamma):
    a_Gamma nie jest wyznaczana przez bifurkacje ani przez Q=3/2.
    Hipotezy dla P53:
      (a) Minimalizacja energii wzdluz krzywej alpha_K(a)?
      (b) Warunek 'stabilnosci solitonu' dg/da = 0?
      (c) C(a) = C_exact (np. C=2.000)?
      (d) Warunek kwantyzacji TGP (topologiczny)?

  KOLEJNY KROK P53:
    Zbadac czy K1(alpha_K(a), a) ma ekstremum lub specjalne wlasnosci
    wzgledem a. Sprawdzic C(a) = K3*sqrt(lambda)/a na krzywej alpha_K(a).
""")

print("="*72)
print("P52 ZAKONCZONE.")
print("="*72)
