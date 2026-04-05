"""
p30_exact_lepton_fit.py
=======================
CEL: KLUCZOWY TEST -- czy TGP z dokladnymi masami leptonowymi daje Q=3/2?

KONTEKST (P29):
  - alpha_0 = 8.5526 -> Q=3/2 przy (r21=206.746, r31=3477.10)
  - Masy leptonowe:                  r21=206.768, r31=3477.65
  - Roznica: 0.011%, 0.016%

PYTANIE KLUCZOWE (P30):
  Czy istnieje (alpha, a_gam, lambda) takie ze JEDNOCZESNIE:
    (1) r21 = 206.768  (dokladna wartosc leptonowa)
    (2) r31 = 3477.65  (dokladna wartosc leptonowa)
    (3) Q  = 3/2 dokladnie

  Jesli TAK: TGP dokladnie reprodukuje masy I Q=3/2 to konsekwencja!
  Jesli NIE: Q=3/2 jest bledne o delta dla dokladnych mas.

PLAN:
  Czesc A: Dla a_gam=0.040 -- znajdz (alpha, lambda) spelniajace (1)+(2)
            Czyli: dopasuj TGP do OBU mas jednoczesnie
  Czesc B: Oblicz Q przy dokladnym dopasowaniu
  Czesc C: Skan a_gam -- czy Q=3/2 przy dokladnych masach dla jakiegos a_gam?
  Czesc D: Porownanie -- punkt TGP a punkt dokladny a krzywa Koide
  Czesc E: Wnioski analityczne

UWAGA: Dotychczas optymalizowalismy:
  - ustalano r31=3477.1 (TGP wartossc, nie dokladna leptonowa 3477.65)
  - szukano alpha minimalizujacego |r21 - 207|
  - Teraz: optymalizuj do (r21, r31) = (206.768, 3477.65)
"""

import numpy as np
from scipy.optimize import brentq, fsolve
import warnings
warnings.filterwarnings('ignore')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

R_MAX  = 60.0
GAMMA  = 1.0

# Dokladne wartosci leptonowe
R21_LEPTON = 206.768    # m_mu^2 / m_e^2
R31_LEPTON = 3477.65    # m_tau^2 / m_e^2

print("P30: Dokladne dopasowanie TGP do mas leptonowych -- test Q=3/2")
print("=" * 70)
print()
print(f"  Cel: znalezc (alpha, a_gam, lambda) dajace r21={R21_LEPTON}, r31={R31_LEPTON}")
print(f"  Pytanie: czy wynikowe Q = 3/2 dokladnie?")
print()

# ============================================================
# NARZEDZIA (jak p27-p29)
# ============================================================
def V_mod(phi, lam):
    return GAMMA/3*phi**3 - GAMMA/4*phi**4 + lam/6*(phi-1)**6

def energy_log(K, alpha, a_gam, lam, N=2000):
    t    = np.linspace(0, 1, N)
    r    = a_gam * (R_MAX/a_gam)**t
    phi  = np.maximum(1.0 + K*np.exp(-r)/r, 1e-10)
    dphi = K*np.exp(-r)*(-r - 1.0)/r**2
    V1   = V_mod(1.0, lam)
    Ek   = 4*np.pi * np.trapezoid(0.5*dphi**2 * (1+alpha/phi) * r**2, r)
    Ep   = 4*np.pi * np.trapezoid((V_mod(phi, lam) - V1) * r**2, r)
    return Ek + Ep

def g_func(K, alpha, a_gam, lam):
    return energy_log(K, alpha, a_gam, lam) / (4*np.pi*K) - 1.0

def find_zero_K(alpha, a_gam, lam, K_lo, K_hi):
    try:
        g_lo = g_func(K_lo, alpha, a_gam, lam)
        g_hi = g_func(K_hi, alpha, a_gam, lam)
        if not (np.isfinite(g_lo) and np.isfinite(g_hi)):
            return np.nan
        if g_lo * g_hi >= 0:
            return np.nan
        return brentq(lambda K: g_func(K, alpha, a_gam, lam), K_lo, K_hi, xtol=1e-9)
    except Exception:
        return np.nan

def find_all_zeros(alpha, a_gam, lam):
    intervals = [(0.001, 0.050), (0.050, 0.500), (0.500, 5.0),
                 (5.0, 50.0), (50.0, 200.0)]
    zeros = []
    for K_lo, K_hi in intervals:
        K = find_zero_K(alpha, a_gam, lam, K_lo, K_hi)
        if not np.isnan(K):
            if not zeros or abs(K - zeros[-1])/zeros[-1] > 0.01:
                zeros.append(K)
    return sorted(zeros)

def koide_Q(k1, k2, k3):
    s = np.sqrt(k1) + np.sqrt(k2) + np.sqrt(k3)
    return s**2 / (k1 + k2 + k3)

def Q_from_ratios(r21, r31):
    s = 1.0 + np.sqrt(r21) + np.sqrt(r31)
    return s**2 / (1.0 + r21 + r31)

def full_solve(alpha, a_gam, lam):
    zeros = find_all_zeros(alpha, a_gam, lam)
    if len(zeros) < 3:
        return np.nan, np.nan, np.nan
    K1, K2, K3 = zeros[0], zeros[1], zeros[2]
    return koide_Q(K1, K2, K3), K2/K1, K3/K1

# ============================================================
# CZESC A: Dokladne dopasowanie (alpha, lambda) -> (r21, r31) = leptony
# ============================================================
print("CZESC A: Dopasowanie (alpha, lambda) do dokladnych mas leptonowych")
print("-" * 70)
print()
print("  Szukam (alpha, lambda) dla a_gam=0.040 takich ze:")
print(f"    r21(alpha, lambda) = {R21_LEPTON}  ORAZ")
print(f"    r31(alpha, lambda) = {R31_LEPTON}")
print()
print("  Strategia: dla kazdego alpha, znajdz lambda s.t. r31=3477.65,")
print("  nastepnie sprawdz r21. Szukaj alpha gdzie r21=206.768.")
print()

AGAM_FIX = 0.040

def find_lambda_for_r31_exact(alpha, agam, r31_target=R31_LEPTON,
                               lam_lo=3e-6, lam_hi=8e-6):
    """Znajdz lambda takie ze r31(alpha, agam, lambda) = r31_target."""
    def f(lam):
        _, _, r31 = full_solve(alpha, agam, lam)
        if np.isnan(r31):
            return np.nan
        return r31 - r31_target
    try:
        f_lo = f(lam_lo)
        f_hi = f(lam_hi)
        if np.isnan(f_lo) or np.isnan(f_hi):
            return np.nan
        if f_lo * f_hi >= 0:
            return np.nan
        return brentq(f, lam_lo, lam_hi, xtol=lam_lo*1e-5, maxiter=60)
    except Exception:
        return np.nan

# Skan alpha -- szukaj gdzie r21=206.768 przy r31=3477.65
alpha_scan = np.linspace(8.40, 8.70, 25)
r21_vals   = []
Q_vals     = []
lam_vals   = []
alpha_ok   = []

print(f"  {'alpha':>8}  {'lambda':>12}  {'r21':>12}  {'r31':>12}  "
      f"{'Q':>12}  {'Q-3/2':>12}")
print("  " + "-"*74)

lam_prev = 5.5e-6
for alpha_s in alpha_scan:
    lam_s = find_lambda_for_r31_exact(alpha_s, AGAM_FIX,
                                       r31_target=R31_LEPTON,
                                       lam_lo=lam_prev*0.7,
                                       lam_hi=lam_prev*1.5)
    if np.isnan(lam_s):
        lam_s = find_lambda_for_r31_exact(alpha_s, AGAM_FIX,
                                           r31_target=R31_LEPTON,
                                           lam_lo=3e-6, lam_hi=9e-6)
    if np.isnan(lam_s):
        continue
    lam_prev = lam_s
    Q_s, r21_s, r31_s = full_solve(alpha_s, AGAM_FIX, lam_s)
    if np.isnan(Q_s):
        continue
    r21_vals.append(r21_s)
    Q_vals.append(Q_s)
    lam_vals.append(lam_s)
    alpha_ok.append(alpha_s)
    mark = " <<" if abs(r21_s - R21_LEPTON) < 0.1 else ""
    print(f"  {alpha_s:>8.4f}  {lam_s:>12.5e}  {r21_s:>12.5f}  "
          f"{r31_s:>12.4f}  {Q_s:>12.8f}  {Q_s-1.5:>+12.8f}{mark}")

print()

# Znajdz alpha gdzie r21 = R21_LEPTON (interpolacja lub brentq)
r21_arr  = np.array(r21_vals)
Q_arr    = np.array(Q_vals)
lam_arr  = np.array(lam_vals)
alp_arr  = np.array(alpha_ok)

alpha_exact = np.nan
lam_exact   = np.nan
Q_exact     = np.nan
r21_exact   = np.nan
r31_exact   = np.nan

# Szukaj przejscia r21 = R21_LEPTON
for i in range(len(r21_arr)-1):
    if (r21_arr[i] - R21_LEPTON) * (r21_arr[i+1] - R21_LEPTON) < 0:
        try:
            def r21_minus_target(alpha):
                lam_t = find_lambda_for_r31_exact(alpha, AGAM_FIX,
                                                   r31_target=R31_LEPTON,
                                                   lam_lo=3e-6, lam_hi=9e-6)
                if np.isnan(lam_t):
                    return np.nan
                _, r21_t, _ = full_solve(alpha, AGAM_FIX, lam_t)
                if np.isnan(r21_t):
                    return np.nan
                return r21_t - R21_LEPTON

            alpha_exact = brentq(r21_minus_target,
                                  alp_arr[i], alp_arr[i+1],
                                  xtol=1e-6, maxiter=50)
            lam_exact = find_lambda_for_r31_exact(alpha_exact, AGAM_FIX,
                                                   r31_target=R31_LEPTON,
                                                   lam_lo=3e-6, lam_hi=9e-6)
            Q_exact, r21_exact, r31_exact = full_solve(alpha_exact, AGAM_FIX, lam_exact)
        except Exception as e:
            print(f"  BLAD brentq: {e}")
        break

print()
if not np.isnan(alpha_exact):
    print(f"  *** DOKLADNE DOPASOWANIE ZNALEZIONE ***")
    print(f"  alpha_exact  = {alpha_exact:.8f}")
    print(f"  lambda_exact = {lam_exact:.8e}")
    print(f"  r21_exact    = {r21_exact:.8f}  (target: {R21_LEPTON})")
    print(f"  r31_exact    = {r31_exact:.4f}   (target: {R31_LEPTON})")
    print()
    print(f"  *** Q PRZY DOKLADNYCH MASACH: Q = {Q_exact:.10f} ***")
    print(f"  *** Q - 3/2 = {Q_exact-1.5:+.10f}              ***")
    print()
    if abs(Q_exact - 1.5) < 1e-4:
        print("  WYNIK: Q jest bardzo bliskie 3/2 przy dokladnych masach!")
        print(f"         Odchylenie: {(Q_exact-1.5)*1e4:.4f} x 10^-4")
    elif abs(Q_exact - 1.5) < 1e-3:
        print("  WYNIK: Q jest bliskie ale nie rowne 3/2.")
    else:
        print(f"  WYNIK: Q = {Q_exact:.6f}, roznica od 3/2 = {Q_exact-1.5:+.6f}")
else:
    print("  UWAGA: Nie znaleziono dokladnego dopasowania w skanie alpha.")
    print("  Sprawdz zakres alpha lub zmien a_gam.")

print()

# ============================================================
# CZESC B: Weryfikacja -- bezposrednie sprawdzenie Q z r21=206.768, r31=3477.65
# ============================================================
print("CZESC B: Weryfikacja -- Q z dokladnych wartosci leptonowych")
print("-" * 70)
print()
print(f"  Q_formula(r21={R21_LEPTON}, r31={R31_LEPTON}) = ", end="")
Q_lept_exact = Q_from_ratios(R21_LEPTON, R31_LEPTON)
print(f"{Q_lept_exact:.12f}")
print(f"  Q - 3/2 = {Q_lept_exact-1.5:+.2e}")
print()

# Krzywa Koide dla r21=206.768
def r31_koide_exact(r21):
    a = 1.0 + np.sqrt(r21)
    disc = 24*a**2 - 12 - 12*r21
    if disc < 0:
        return np.nan
    x = (4*a + np.sqrt(disc)) / 2.0
    return x**2

r31_K_lept = r31_koide_exact(R21_LEPTON)
print(f"  Krzywa Koide dla r21={R21_LEPTON}: r31_K = {r31_K_lept:.6f}")
print(f"  r31_lepton = {R31_LEPTON}")
print(f"  Roznica:    {R31_LEPTON - r31_K_lept:+.6f}")
print()
if abs(R31_LEPTON - r31_K_lept) < 0.01:
    print("  => Masy leptonowe leza NA krzywej Koide (jak oczekiwano)")
else:
    print(f"  => Masy leptonowe leza POZA krzywq Koide o {R31_LEPTON - r31_K_lept:.4f}")

print()
print(f"  Porownanie Q przy roznych punktach:")
points = [
    ("Leptony (exp)",     R21_LEPTON, R31_LEPTON),
    ("TGP glowny (P28)",  206.999,    3477.134),
    ("alpha_0 (P29)",     206.746,    3477.100),
]
if not np.isnan(r21_exact):
    points.append(("TGP dokladne (P30)", r21_exact, r31_exact))

print(f"  {'Punkt':>22}  {'r21':>10}  {'r31':>10}  {'Q':>14}  {'Q-3/2':>14}")
print("  " + "-"*76)
for name, r21_p, r31_p in points:
    Q_p = Q_from_ratios(r21_p, r31_p)
    print(f"  {name:>22}  {r21_p:>10.4f}  {r31_p:>10.4f}  "
          f"{Q_p:>14.10f}  {Q_p-1.5:>+14.10f}")

print()

# ============================================================
# CZESC C: Skan a_gam przy dokladnych masach leptonowych
# ============================================================
print("CZESC C: Skan a_gam -- Q(a_gam) przy dokladnych masach")
print("-" * 70)
print()

agam_scan_c = np.array([0.025, 0.030, 0.035, 0.038, 0.040, 0.042, 0.045, 0.050])
Q_agam_c    = []
alpha_agam_c = []
lam_agam_c   = []

print(f"  {'a_gam':>8}  {'alpha*':>10}  {'lambda*':>12}  {'r21':>10}  {'r31':>10}  "
      f"{'Q':>12}  {'Q-3/2':>12}")
print("  " + "-"*80)

for agam_c in agam_scan_c:
    # Dla kazdego a_gam, znajdz (alpha, lambda) takie ze r21=lept, r31=lept
    alpha_scan_c2 = np.linspace(5.5, 9.5, 30)
    best_alpha = np.nan
    best_lam   = np.nan
    best_r21   = np.nan
    best_r31   = np.nan
    best_Q     = np.nan

    r21_scan_c2 = []
    alp_c2 = []
    lam_c2 = []

    lam_g = 5.5e-6 * (agam_c / 0.040)**1.5  # skalowanie z a_gam
    for alpha_cc in alpha_scan_c2:
        lam_t = find_lambda_for_r31_exact(alpha_cc, agam_c,
                                           r31_target=R31_LEPTON,
                                           lam_lo=lam_g*0.3,
                                           lam_hi=lam_g*3.0)
        if np.isnan(lam_t):
            continue
        lam_g = lam_t
        _, r21_t, _ = full_solve(alpha_cc, agam_c, lam_t)
        if np.isnan(r21_t):
            continue
        r21_scan_c2.append(r21_t)
        alp_c2.append(alpha_cc)
        lam_c2.append(lam_t)

    r21_c2 = np.array(r21_scan_c2)
    alp_c2 = np.array(alp_c2)
    lam_c2 = np.array(lam_c2)

    # Szukaj przejscia r21 = R21_LEPTON
    if len(r21_c2) < 2:
        print(f"  {agam_c:>8.4f}  {'brak danych':>50}")
        continue

    alpha_found = np.nan
    lam_found   = np.nan
    for ii in range(len(r21_c2)-1):
        if (r21_c2[ii]-R21_LEPTON) * (r21_c2[ii+1]-R21_LEPTON) < 0:
            try:
                lam_hint = (lam_c2[ii] + lam_c2[ii+1]) / 2

                def r21_minus_tgt(a):
                    lt = find_lambda_for_r31_exact(a, agam_c,
                                                   r31_target=R31_LEPTON,
                                                   lam_lo=lam_hint*0.5,
                                                   lam_hi=lam_hint*2.0)
                    if np.isnan(lt):
                        return np.nan
                    _, r21v, _ = full_solve(a, agam_c, lt)
                    return (r21v if not np.isnan(r21v) else 1e10) - R21_LEPTON

                alpha_found = brentq(r21_minus_tgt, alp_c2[ii], alp_c2[ii+1],
                                     xtol=1e-5, maxiter=40)
                lam_found = find_lambda_for_r31_exact(alpha_found, agam_c,
                                                       r31_target=R31_LEPTON,
                                                       lam_lo=lam_hint*0.5,
                                                       lam_hi=lam_hint*2.0)
            except Exception:
                pass
            break

    if not np.isnan(alpha_found) and not np.isnan(lam_found):
        Q_f, r21_f, r31_f = full_solve(alpha_found, agam_c, lam_found)
        Q_agam_c.append(Q_f)
        alpha_agam_c.append(alpha_found)
        lam_agam_c.append(lam_found)
        mark = " <" if abs(Q_f-1.5) < 0.0005 else ""
        print(f"  {agam_c:>8.4f}  {alpha_found:>10.5f}  {lam_found:>12.5e}  "
              f"{r21_f:>10.4f}  {r31_f:>10.4f}  "
              f"{Q_f:>12.8f}  {Q_f-1.5:>+12.8f}{mark}")
    else:
        print(f"  {agam_c:>8.4f}  {'nie znaleziono':>55}")

print()

# Szukaj a_gam gdzie Q=3/2 przy dokladnych masach
Q_ag_c = np.array(Q_agam_c)
if len(Q_ag_c) >= 2:
    print("  Przeciecia Q=3/2 wzgledem a_gam (przy dokladnych masach):")
    found_crossing = False
    for ii in range(len(Q_ag_c)-1):
        if (Q_ag_c[ii]-1.5) * (Q_ag_c[ii+1]-1.5) < 0:
            found_crossing = True
            ag_lo = agam_scan_c[ii]
            ag_hi = agam_scan_c[ii+1]
            print(f"    Przejscie miedzy a_gam={ag_lo:.4f} a a_gam={ag_hi:.4f}")

            # Doprecyzuj bisection dla a_gam
            try:
                def Q_minus_3_2_agam(ag):
                    # Znajdz alpha dla tego ag i dokladnych mas
                    alp_scan_q = np.linspace(5.5, 9.5, 15)
                    r21_sq = []
                    alp_sq = []
                    lam_sq = []
                    lam_hq = 5.5e-6 * (ag/0.040)**1.5
                    for a_q in alp_scan_q:
                        lt_q = find_lambda_for_r31_exact(a_q, ag,
                                                          r31_target=R31_LEPTON,
                                                          lam_lo=lam_hq*0.3,
                                                          lam_hi=lam_hq*3.0)
                        if np.isnan(lt_q):
                            continue
                        lam_hq = lt_q
                        _, r21_q, _ = full_solve(a_q, ag, lt_q)
                        if np.isnan(r21_q):
                            continue
                        r21_sq.append(r21_q)
                        alp_sq.append(a_q)
                        lam_sq.append(lt_q)
                    if len(r21_sq) < 2:
                        return np.nan
                    r21_sq = np.array(r21_sq)
                    for jj in range(len(r21_sq)-1):
                        if (r21_sq[jj]-R21_LEPTON)*(r21_sq[jj+1]-R21_LEPTON) < 0:
                            try:
                                lh2 = (lam_sq[jj]+lam_sq[jj+1])/2
                                def r21_q2(a):
                                    lt2 = find_lambda_for_r31_exact(a, ag,
                                              r31_target=R31_LEPTON,
                                              lam_lo=lh2*0.5, lam_hi=lh2*2.0)
                                    if np.isnan(lt2):
                                        return np.nan
                                    _, rv, _ = full_solve(a, ag, lt2)
                                    return (rv if not np.isnan(rv) else 1e10) - R21_LEPTON
                                a_fin = brentq(r21_q2, alp_sq[jj], alp_sq[jj+1],
                                               xtol=1e-5, maxiter=30)
                                l_fin = find_lambda_for_r31_exact(a_fin, ag,
                                            r31_target=R31_LEPTON,
                                            lam_lo=lh2*0.5, lam_hi=lh2*2.0)
                                if np.isnan(l_fin):
                                    return np.nan
                                Q_fin, _, _ = full_solve(a_fin, ag, l_fin)
                                return Q_fin - 1.5 if not np.isnan(Q_fin) else np.nan
                            except Exception:
                                pass
                    return np.nan

                ag_star = brentq(Q_minus_3_2_agam, ag_lo, ag_hi,
                                 xtol=1e-4, maxiter=20)
                print(f"    a_gam_Koide = {ag_star:.6f}")
            except Exception as ee:
                print(f"    Blad precyzji: {ee}")
    if not found_crossing:
        print("    Brak przeciecia Q=3/2 w zbadanym zakresie a_gam")

print()

# ============================================================
# CZESC D: Zestawienie koncowe
# ============================================================
print("CZESC D: Zestawienie koncowe -- co wiemy")
print("-" * 70)
print()

print("  WSZYSTKIE KLUCZOWE PUNKTY:")
print()
print(f"  {'Punkt':>25}  {'alpha':>8}  {'a_gam':>6}  {'r21':>10}  {'r31':>10}  "
      f"{'Q':>12}  {'Q-3/2':>12}")
print("  " + "-"*86)

# Leptony
Q_l = Q_from_ratios(R21_LEPTON, R31_LEPTON)
print(f"  {'Leptony (exp)':>25}  {'---':>8}  {'---':>6}  "
      f"{R21_LEPTON:>10.4f}  {R31_LEPTON:>10.4f}  "
      f"{Q_l:>12.10f}  {Q_l-1.5:>+12.2e}")

# TGP glowny z P28
print(f"  {'TGP glowny (P28)':>25}  {8.5616:>8.4f}  {0.040:>6.4f}  "
      f"{206.999:>10.4f}  {3477.134:>10.4f}  "
      f"{1.50025:>12.8f}  {+0.00025:>+12.2e}")

# alpha_0 z P29
print(f"  {'alpha_0=3/2 (P29)':>25}  {8.5526:>8.4f}  {0.040:>6.4f}  "
      f"{206.746:>10.4f}  {3477.10:>10.4f}  "
      f"{1.50000:>12.8f}  {0.0:>+12.2e}")

# Dokladne dopasowanie (P30)
if not np.isnan(alpha_exact):
    print(f"  {'Dokladne r21+r31 (P30)':>25}  {alpha_exact:>8.5f}  {0.040:>6.4f}  "
          f"{r21_exact:>10.5f}  {r31_exact:>10.4f}  "
          f"{Q_exact:>12.10f}  {Q_exact-1.5:>+12.2e}")

# Punkt Koide dla r21=206.768
r31_K_lepton = r31_koide_exact(R21_LEPTON)
print(f"  {'Koide(r21=206.768)':>25}  {'---':>8}  {'---':>6}  "
      f"{R21_LEPTON:>10.4f}  {r31_K_lepton:>10.4f}  "
      f"{1.50000000:>12.8f}  {0.0:>+12.2e}")

print()

# ============================================================
# CZESC E: Wnioski
# ============================================================
print("CZESC E: Wnioski i interpretacja fizyczna")
print("-" * 70)
print()

print(f"  1. Dokladne masy leptonowe (r21={R21_LEPTON}, r31={R31_LEPTON}):")
print(f"     Q_formula = {Q_from_ratios(R21_LEPTON, R31_LEPTON):.12f}")
print(f"     Q - 3/2   = {Q_from_ratios(R21_LEPTON, R31_LEPTON)-1.5:+.2e}")
print()

if not np.isnan(Q_exact):
    print(f"  2. TGP dokladnie dopasowane do (r21={R21_LEPTON}, r31={R31_LEPTON}):")
    print(f"     wymaga alpha={alpha_exact:.6f}, lambda={lam_exact:.4e}")
    print(f"     daje Q = {Q_exact:.12f}")
    print(f"     Q - 3/2 = {Q_exact-1.5:+.2e}")
    print()

    if abs(Q_exact - 1.5) < 1e-5:
        print("  *** WYNIK KLUCZOWY: TGP dopasowane do DOKLADNYCH mas leptonowych ***")
        print("  *** daje Q = 3/2 z dokladnoscia < 1e-5 !!                       ***")
        print()
        print("  Implikacja: Q = 3/2 JEST KONSEKWENCJA tego, ze (r21, r31)")
        print("  odpowiadaja masom leptonowym w TGP.")
        print("  Nie jest to dodatkowe zalozenie -- wynika z dopasowania mas!")
    elif abs(Q_exact - 1.5) < 1e-3:
        print("  WYNIK: TGP przy dokladnych masach daje Q BLISKIE 3/2.")
        print(f"  Odchylenie: {abs(Q_exact-1.5):.2e}  (mniejsze niz przy r31=3477.1: {0.00025:.2e})")
    else:
        print(f"  WYNIK: Q = {Q_exact:.6f} przy dokladnych masach.")
        print("  Q nie jest 3/2 nawet przy dokladnym dopasowaniu.")

print()
print("  STATUS PROBLEMU O-K1 po P30:")
print("  OK-1 (ograniczenie Q>=3/2): ODPOWIEDZIANE -- NIE (P28)")
print("  OK-2 (predykcja Q=3/2):     -> patrz wynik P30 powyzej")
print("  OK-3 (szczelina stala):      ODPOWIEDZIANE -- NIE (P28)")
print("  OK-4 (punkt specjalny):      CZESCIOWO -- P29,P30")
print()

# ============================================================
# CZESC F: Wykres
# ============================================================
print("CZESC F: Wykres")
print("-" * 70)

try:
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle("P30: Q przy dokladnych masach leptonowych", fontsize=13)

    # Panel 1: Q(alpha) przy r31=3477.65 (Czesc A)
    ax1 = axes[0]
    if len(alp_arr) > 1:
        ax1.plot(alp_arr, Q_arr, 'b-o', markersize=5, label='Q(alpha) r31=3477.65')
        ax1.axhline(1.5, color='r', linestyle='--', label='Q=3/2')
        ax1.axhline(Q_from_ratios(R21_LEPTON, R31_LEPTON), color='purple',
                   linestyle=':', label=f'Q_lepton={Q_from_ratios(R21_LEPTON, R31_LEPTON):.6f}')
        ax1.axvline(8.5616, color='g', linestyle=':', alpha=0.7, label='alpha_opt')
        if not np.isnan(alpha_exact):
            ax1.axvline(alpha_exact, color='orange', linestyle='-.',
                       label=f'alpha_exact={alpha_exact:.4f}')
    ax1.set_xlabel('alpha')
    ax1.set_ylabel('Q')
    ax1.set_title('Q(alpha) przy r31=3477.65 (dokladne masy)')
    ax1.legend(fontsize=7)
    ax1.grid(True, alpha=0.3)

    # Panel 2: Porownanie punktow w (r21, r31)
    ax2 = axes[1]
    # Krzywa Koide
    r21_line = np.linspace(200, 215, 300)
    r31_K_line = np.array([r31_koide_exact(r) for r in r21_line])
    mask_K = ~np.isnan(r31_K_line)
    ax2.plot(r21_line[mask_K], r31_K_line[mask_K], 'r-', lw=2, label='Krzywa Koide Q=3/2')

    # Trajektoria TGP (Czesc A: r31=3477.65, r21 zmienny)
    if len(alp_arr) > 1:
        # Zbierz r21 z skanu
        r21_trj = r21_arr
        r31_trj = np.full_like(r21_trj, R31_LEPTON)
        ax2.plot(r21_trj, r31_trj, 'b--', lw=1, alpha=0.5, label='Trajektoria (r31=3477.65)')

    # Punkty kluczowe
    ax2.scatter([R21_LEPTON], [R31_LEPTON], s=200, c='red', marker='*',
               zorder=10, label=f'Leptony ({R21_LEPTON}, {R31_LEPTON})')
    ax2.scatter([206.999], [3477.134], s=100, c='green', marker='s',
               zorder=10, label='TGP glowny')
    ax2.scatter([206.746], [3477.10], s=100, c='orange', marker='^',
               zorder=10, label='alpha_0 P29')
    if not np.isnan(r21_exact):
        ax2.scatter([r21_exact], [r31_exact], s=150, c='purple', marker='D',
                   zorder=10, label=f'TGP dokladne P30')

    ax2.set_xlim(205, 210)
    ax2.set_ylim(3470, 3490)
    ax2.set_xlabel('r21')
    ax2.set_ylabel('r31')
    ax2.set_title('Punkty TGP vs krzywa Koide (pow. r21=206-210)')
    ax2.legend(fontsize=7)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('p30_exact_lepton.png', dpi=120, bbox_inches='tight')
    print("  Wykres zapisany: p30_exact_lepton.png")
except Exception as e:
    print(f"  UWAGA: Blad wykresu: {e}")

print()

# ============================================================
# PODSUMOWANIE P30
# ============================================================
print()
print("=" * 70)
print("PODSUMOWANIE P30: Czy Q=3/2 przy dokladnych masach leptonowych?")
print("=" * 70)
print()
Q_lept_v = Q_from_ratios(R21_LEPTON, R31_LEPTON)
print(f"  Q(r21={R21_LEPTON}, r31={R31_LEPTON}) = {Q_lept_v:.12f}")
print(f"  Q - 3/2 = {Q_lept_v-1.5:+.4e}")
print()
if not np.isnan(Q_exact):
    print(f"  TGP przy dokladnych masach:")
    print(f"    alpha = {alpha_exact:.6f},  lambda = {lam_exact:.5e}")
    print(f"    Q     = {Q_exact:.12f}")
    print(f"    Q-3/2 = {Q_exact-1.5:+.4e}")
print()
print("P30 zakonczone.")
