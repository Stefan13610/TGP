"""
p29_koide_zero.py
=================
CEL: Zlokalizowanie "punktu Koide" w przestrzeni parametrow TGP.

KLUCZOWE ODKRYCIE Z P28:
  alpha=8.50  -> Q(lambda*) = 1.49906  (Q < 3/2)
  alpha=8.5616 -> Q(lambda*) = 1.50025  (Q > 3/2)

  => Miedzy alpha=8.50 a alpha=8.5616 istnieje PRZEJSCIE przez Q=3/2!
  => Istnieje alpha_0 (punkt zerowy) gdzie Q(lambda*(alpha_0)) = 3/2 dokladnie.

PLAN P29:
  Czesc A: Znalezienie alpha_0 (bisection + brentq) -- Q(alpha) = 3/2
  Czesc B: Skan Q(alpha) przy stalym a_gam=0.040 -- pelna krzywa Q(alpha)
  Czesc C: Skan Q(a_gam) przy stalym alpha=8.5616 -- Q(a_gam)
  Czesc D: Mapa 2D Q(alpha, a_gam) -- krzywa zerowa Q=3/2 w przestrzeni par.
  Czesc E: Co oznacza punkt alpha_0?
           - Jakie r21, r31 daje TGP przy (alpha_0, lambda*(alpha_0))?
           - Czy (r21, r31) lezy na krzywej Koide?
           - Jak daleko jest od mas leptonowych?

HIPOTEZA:
  Punkt (alpha_0, a_gam_0) gdzie Q=3/2 dokladnie odpowiada dokladnym masom
  leptonowym (r21=206.768, r31=3477.65). Jesli tak: TGP PRZEWIDUJE Q=3/2
  dla dokladnych mas (ale nie dla TGP-optymalnego r21=207.0).
"""

import numpy as np
from scipy.optimize import brentq
import warnings
warnings.filterwarnings('ignore')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

R_MAX  = 60.0
GAMMA  = 1.0

print("P29: Punkt zerowy Q=3/2 w przestrzeni (alpha, a_gam)")
print("=" * 70)
print()

# ============================================================
# NARZEDZIA PODSTAWOWE (jak p27, p28)
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
    E = energy_log(K, alpha, a_gam, lam)
    return E / (4*np.pi*K) - 1.0

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
    """Zwraca (Q, r21, r31, K1, K2, K3) lub (nan,...) jesli nie ma 3 zer."""
    zeros = find_all_zeros(alpha, a_gam, lam)
    if len(zeros) < 3:
        return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
    K1, K2, K3 = zeros[0], zeros[1], zeros[2]
    Q  = koide_Q(K1, K2, K3)
    return Q, K2/K1, K3/K1, K1, K2, K3

def find_lambda_for_r31(alpha, a_gam, r31_target, lam_lo, lam_hi):
    """
    Dla danych (alpha, a_gam) znajdz lambda takie ze r31(lambda) = r31_target.
    r31 rosnie z lambda (wieksza lambda -> mniejsze K*3 -> ale K*1 rowniez maleje
    -> r31 rosnie lub maleje zalezne od parametrow).
    Sprawdz najpierw monotonocznosc.
    """
    def f(lam):
        Q, r21, r31, K1, K2, K3 = full_solve(alpha, a_gam, lam)
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

# Wartosc referencyjna leptonow
R21_LEPTON = 206.768   # m_mu / m_e
R31_LEPTON = 3477.65   # m_tau / m_e

def r31_koide_from_r21(r21):
    """r31 na krzywej Q=3/2 dla zadanego r21."""
    a = 1.0 + np.sqrt(r21)
    disc = 24*a**2 - 12 - 12*r21
    if disc < 0:
        return np.nan
    x = (4*a + np.sqrt(disc)) / 2.0
    return x**2

# ============================================================
# CZESC A: PUNKT ZEROWY alpha_0 -- Q(alpha, lambda*(alpha)) = 3/2
# ============================================================
print("CZESC A: Znalezienie alpha_0 -- punkt gdzie Q(lambda*) = 3/2")
print("-" * 70)
print()
print("  Z P28 wiemy:")
print("    alpha=8.50,   lambda*=5.45e-6:  Q=1.49906 (Q < 3/2)")
print("    alpha=8.5616, lambda*=5.501e-6: Q=1.50025 (Q > 3/2)")
print()
print("  => Zero Q(alpha)-3/2 istnieje miedzy 8.50 a 8.5616")
print()

AGAM_FIX = 0.040

def Q_at_optimal(alpha, agam, r21_target=207.0, r31_target=3477.1,
                 lam_guess=5.5e-6):
    """
    Znajdz lambda takie ze r31 = r31_target (przy danym alpha, agam),
    zwroc Q przy tym lambda.
    Strategia: skanuj lambda wokol lam_guess.
    """
    # Szukaj lambda gdzie r31 = r31_target
    lam = find_lambda_for_r31(alpha, agam, r31_target,
                               lam_guess*0.5, lam_guess*3.0)
    if np.isnan(lam):
        # Probuj szerszy zakres
        lam = find_lambda_for_r31(alpha, agam, r31_target,
                                   lam_guess*0.1, lam_guess*10.0)
    if np.isnan(lam):
        return np.nan, np.nan, np.nan, np.nan
    Q, r21, r31, K1, _, _ = full_solve(alpha, agam, lam)
    return Q, r21, r31, lam

# Weryfikacja punktow znanych z P28
print("  Weryfikacja punktow z P28:")
for alpha_test, lam_guess_t in [(8.50, 5.45e-6), (8.5616, 5.501e-6)]:
    Q_t, r21_t, r31_t, lam_t = Q_at_optimal(alpha_test, AGAM_FIX, lam_guess=lam_guess_t)
    print(f"    alpha={alpha_test}: lambda={lam_t:.4e}, r21={r21_t:.3f}, "
          f"r31={r31_t:.3f}, Q={Q_t:.8f}, Q-3/2={Q_t-1.5:+.8f}")

print()

# Szukanie alpha_0 metoda brentq
def Q_minus_3_2_alpha(alpha, agam):
    Q_v, _, _, _ = Q_at_optimal(alpha, agam)
    if np.isnan(Q_v):
        return np.nan
    return Q_v - 1.5

print("  Szukanie alpha_0 (brentq, a_gam=0.040):")
alpha_lo, alpha_hi = 8.50, 8.5616

try:
    alpha_0 = brentq(
        lambda a: Q_minus_3_2_alpha(a, AGAM_FIX),
        alpha_lo, alpha_hi,
        xtol=1e-6, maxiter=40
    )
    Q_0, r21_0, r31_0, lam_0 = Q_at_optimal(alpha_0, AGAM_FIX)
    r31_K_0 = r31_koide_from_r21(r21_0)

    print(f"  ZNALEZIONO alpha_0 = {alpha_0:.8f}")
    print(f"  lambda*(alpha_0)   = {lam_0:.6e}")
    print(f"  Q(alpha_0, lambda*) = {Q_0:.10f}")
    print(f"  r21 = {r21_0:.6f}  (leptony: {R21_LEPTON})")
    print(f"  r31 = {r31_0:.4f}  (leptony: {R31_LEPTON})")
    print(f"  r31_Koide(r21_0) = {r31_K_0:.4f}")
    print(f"  Czy r31 = r31_Koide?  delta = {r31_0 - r31_K_0:+.6f}")
    print()
    print(f"  INTERPRETACJA:")
    if abs(r21_0 - R21_LEPTON) < 0.5:
        print(f"  => r21_0 = {r21_0:.3f} BLISKIE wartosci leptonowej ({R21_LEPTON})!")
    if abs(r31_0 - R31_LEPTON) < 2.0:
        print(f"  => r31_0 = {r31_0:.2f} BLISKIE wartosci leptonowej ({R31_LEPTON})!")
    print(f"  => Odleglosc od mas leptonowych: dr21={r21_0-R21_LEPTON:+.4f}, "
          f"dr31={r31_0-R31_LEPTON:+.4f}")
    alpha_0_found = alpha_0
    lam_0_found   = lam_0
    r21_0_found   = r21_0
    r31_0_found   = r31_0
except Exception as e:
    print(f"  BLAD: {e}")
    alpha_0_found = np.nan

print()

# ============================================================
# CZESC B: SKAN Q(alpha) przy stalym a_gam = 0.040
# ============================================================
print("CZESC B: Skan Q(alpha) -- pelna krzywa (a_gam=0.040)")
print("-" * 70)
print()

alpha_scan = np.linspace(8.2, 9.0, 20)
Q_alpha    = []
r21_alpha  = []
r31_alpha  = []
lam_alpha  = []
alpha_valid = []

print(f"  {'alpha':>8}  {'lambda*':>12}  {'r21':>10}  {'r31':>12}  "
      f"{'Q':>10}  {'Q-3/2':>12}")
print("  " + "-"*70)

lam_prev = 5.5e-6
for alpha_s in alpha_scan:
    Q_s, r21_s, r31_s, lam_s = Q_at_optimal(alpha_s, AGAM_FIX,
                                              lam_guess=lam_prev)
    if not np.isnan(Q_s):
        lam_prev = lam_s  # tracking dla lepszego guess
        Q_alpha.append(Q_s)
        r21_alpha.append(r21_s)
        r31_alpha.append(r31_s)
        lam_alpha.append(lam_s)
        alpha_valid.append(alpha_s)
        marker = " <<< Q=3/2" if abs(Q_s - 1.5) < 0.001 else ""
        print(f"  {alpha_s:>8.4f}  {lam_s:>12.4e}  {r21_s:>10.4f}  {r31_s:>12.3f}  "
              f"{Q_s:>10.6f}  {Q_s-1.5:>+12.6f}{marker}")

print()
if Q_alpha:
    Q_arr    = np.array(Q_alpha)
    r21_arr  = np.array(r21_alpha)
    r31_arr  = np.array(r31_alpha)
    lam_arr  = np.array(lam_alpha)
    alp_arr  = np.array(alpha_valid)

    # Szukaj gdzie Q zmienia znak (przejscia przez 3/2)
    crossings = []
    for i in range(len(Q_arr)-1):
        if (Q_arr[i] - 1.5) * (Q_arr[i+1] - 1.5) < 0:
            crossings.append((alp_arr[i], alp_arr[i+1]))
    print(f"  Liczba przeciec Q=3/2 w skanowanym zakresie: {len(crossings)}")
    for a_lo_c, a_hi_c in crossings:
        print(f"    Przejscie miedzy alpha={a_lo_c:.4f} a alpha={a_hi_c:.4f}")
    print()
    print(f"  Q_min na siatce: {np.min(Q_arr):.8f} (alpha={alp_arr[np.argmin(Q_arr)]:.4f})")
    print(f"  Q_max na siatce: {np.max(Q_arr):.8f} (alpha={alp_arr[np.argmax(Q_arr)]:.4f})")

print()

# ============================================================
# CZESC C: SKAN Q(a_gam) przy stalym alpha = 8.5616
# ============================================================
print("CZESC C: Skan Q(a_gam) -- (alpha=8.5616)")
print("-" * 70)
print()

ALPHA_FIX = 8.5616
agam_scan = np.linspace(0.030, 0.055, 15)
Q_agam    = []
agam_valid = []

print(f"  {'a_gam':>8}  {'lambda*':>12}  {'r21':>10}  {'r31':>12}  "
      f"{'Q':>10}  {'Q-3/2':>12}")
print("  " + "-"*70)

lam_prev = 5.5e-6
for agam_s in agam_scan:
    Q_s, r21_s, r31_s, lam_s = Q_at_optimal(ALPHA_FIX, agam_s, lam_guess=lam_prev)
    if not np.isnan(Q_s):
        lam_prev = lam_s
        Q_agam.append(Q_s)
        agam_valid.append(agam_s)
        marker = " <<< Q=3/2" if abs(Q_s - 1.5) < 0.001 else ""
        print(f"  {agam_s:>8.4f}  {lam_s:>12.4e}  {r21_s:>10.4f}  {r31_s:>12.3f}  "
              f"{Q_s:>10.6f}  {Q_s-1.5:>+12.6f}{marker}")

print()
if len(Q_agam) > 2:
    Q_ag_arr = np.array(Q_agam)
    crossings_ag = []
    for i in range(len(Q_ag_arr)-1):
        if (Q_ag_arr[i] - 1.5) * (Q_ag_arr[i+1] - 1.5) < 0:
            crossings_ag.append((agam_valid[i], agam_valid[i+1]))
    print(f"  Przeciecia Q=3/2 wzgl. a_gam: {len(crossings_ag)}")
    for c in crossings_ag:
        print(f"    Przejscie miedzy a_gam={c[0]:.4f} a a_gam={c[1]:.4f}")

print()

# Szukaj a_gam_0
print("  Szukanie a_gam_0 (jesli istnieje):")
if len(Q_agam) >= 2:
    Q_ag_arr = np.array(Q_agam)
    for i in range(len(Q_ag_arr)-1):
        if (Q_ag_arr[i] - 1.5) * (Q_ag_arr[i+1] - 1.5) < 0:
            try:
                agam_0 = brentq(
                    lambda ag: Q_at_optimal(ALPHA_FIX, ag)[0] - 1.5,
                    agam_valid[i], agam_valid[i+1],
                    xtol=1e-6, maxiter=40
                )
                Q_ag0, r21_ag0, r31_ag0, lam_ag0 = Q_at_optimal(ALPHA_FIX, agam_0)
                print(f"  ZNALEZIONO a_gam_0 = {agam_0:.8f}")
                print(f"  lambda*(a_gam_0) = {lam_ag0:.6e}")
                print(f"  Q(a_gam_0) = {Q_ag0:.10f}")
                print(f"  r21 = {r21_ag0:.6f},  r31 = {r31_ag0:.4f}")
                print(f"  Leptony: r21={R21_LEPTON}, r31={R31_LEPTON}")
                print(f"  delta r21 = {r21_ag0-R21_LEPTON:+.4f}")
                print(f"  delta r31 = {r31_ag0-R31_LEPTON:+.4f}")
            except Exception as e:
                print(f"  Blad brentq: {e}")
            break
    else:
        print("  Brak przeciecia Q=3/2 w przeskanowanym zakresie a_gam.")

print()

# ============================================================
# CZESC D: MAPA 2D Q(alpha, a_gam) -- krzywa zerowa
# ============================================================
print("CZESC D: Mapa 2D |Q-3/2| w przestrzeni (alpha, a_gam)")
print("-" * 70)
print()

# Mniejsza siatka 2D -- ekonomiczna ale pouczajaca
alpha_2d = np.linspace(8.35, 8.75, 9)
agam_2d  = np.array([0.035, 0.038, 0.040, 0.042, 0.045])

print(f"  Siatka: {len(alpha_2d)} x {len(agam_2d)} = {len(alpha_2d)*len(agam_2d)} punktow")
print()

Q_map_2d = np.full((len(agam_2d), len(alpha_2d)), np.nan)

print(f"  {'':>8}", end="")
for ag in agam_2d:
    print(f"  a_gam={ag:.3f}", end="")
print()
print("  " + "-" * 70)

for j, alpha_j in enumerate(alpha_2d):
    print(f"  alpha={alpha_j:.3f}", end="")
    lam_g = 5.5e-6
    for i, agam_i in enumerate(agam_2d):
        Q_v, _, _, lam_v = Q_at_optimal(alpha_j, agam_i, lam_guess=lam_g)
        Q_map_2d[i, j] = Q_v
        if not np.isnan(lam_v):
            lam_g = lam_v
        dQ = Q_v - 1.5 if not np.isnan(Q_v) else np.nan
        marker = "*" if (not np.isnan(dQ) and abs(dQ) < 0.001) else " "
        val = f"{dQ:+.5f}" if not np.isnan(dQ) else "  ----  "
        print(f"  {marker}{val}", end="")
    print()

print()
print("  (* oznacza |Q-3/2| < 0.001)")
print()

# Znajdz minimum |Q-3/2| na siatce 2D
min_dQ = np.inf
min_alpha = np.nan
min_agam  = np.nan
for j, alpha_j in enumerate(alpha_2d):
    for i, agam_i in enumerate(agam_2d):
        Q_v = Q_map_2d[i, j]
        if not np.isnan(Q_v) and abs(Q_v - 1.5) < min_dQ:
            min_dQ = abs(Q_v - 1.5)
            min_alpha = alpha_j
            min_agam  = agam_i

print(f"  Minimum |Q-3/2| na siatce: {min_dQ:.8f}")
print(f"  Przy: alpha={min_alpha:.4f}, a_gam={min_agam:.4f}")
print(f"  Punkt optymalny (dopasowany do mas): alpha=8.5616, a_gam=0.040")
print()

# ============================================================
# CZESC E: INTERPRETACJA FIZYCZNA
# ============================================================
print("CZESC E: Interpretacja fizyczna")
print("-" * 70)
print()
print("  ZESTAWIENIE PUNKTOW:")
print(f"  {'Punkt':>20}  {'alpha':>8}  {'a_gam':>8}  {'r21':>10}  {'r31':>12}  {'Q':>10}")
print("  " + "-"*72)

# Leptony (eksperymentalne)
Q_lept = Q_from_ratios(R21_LEPTON, R31_LEPTON)
print(f"  {'Leptony (exp)':>20}  {'---':>8}  {'---':>8}  "
      f"{R21_LEPTON:>10.4f}  {R31_LEPTON:>12.4f}  {Q_lept:>10.8f}")

# TGP optymalny (P28)
print(f"  {'TGP glowny':>20}  {8.5616:>8.4f}  {0.040:>8.4f}  "
      f"{207.0:>10.4f}  {3477.10:>12.4f}  {1.50025:>10.8f}")

# Punkt alpha_0
if not np.isnan(alpha_0_found):
    print(f"  {'TGP alpha_0 (Q=3/2)':>20}  {alpha_0_found:>8.4f}  {0.040:>8.4f}  "
          f"{r21_0_found:>10.4f}  {r31_0_found:>12.4f}  {1.50000000:>10.8f}")

# Krzywa Koide dla r21=207
print(f"  {'Koide(r21=207)':>20}  {'---':>8}  {'---':>8}  "
      f"{207.0:>10.4f}  {r31_koide_from_r21(207):>12.4f}  {1.50000000:>10.8f}")

print()
print("  KLUCZOWE WNIOSKI:")
print()
print("  1. Punkt alpha_0 (Q=3/2 dokladnie):")
if not np.isnan(alpha_0_found):
    dr21 = r21_0_found - R21_LEPTON
    dr31 = r31_0_found - R31_LEPTON
    print(f"     alpha_0 = {alpha_0_found:.6f}  (vs optimal 8.5616, roznica {alpha_0_found-8.5616:+.6f})")
    print(f"     r21_0 = {r21_0_found:.4f} (leptony: {R21_LEPTON}, delta={dr21:+.4f})")
    print(f"     r31_0 = {r31_0_found:.2f} (leptony: {R31_LEPTON}, delta={dr31:+.4f})")
    if abs(dr21) < 1.0 and abs(dr31) < 5.0:
        print()
        print("  *** Q=3/2 jest osiagane przy (r21, r31) BLISKICH masom leptonowym! ***")
        print("  *** Oznacza to: jesli dopasujemy TGP do DOKLADNYCH mas ->  Q=3/2 ***")
    elif abs(dr21) < 2.0 and abs(dr31) < 20.0:
        print()
        print("  -> Q=3/2 wymaga nieznacznej korekty parametrow wzgledem optimum.")
    else:
        print()
        print("  -> Q=3/2 daleko od punktu optymalnego -- przypadkowe przejscie.")
print()
print("  2. Monotonizacja Q(alpha):")
if Q_alpha:
    dQ_darr = np.diff(np.array(Q_alpha)) / np.diff(np.array(alpha_valid))
    print(f"     dQ/dalpha (skan) ~ {np.mean(dQ_darr):.5f} (>0: Q rosnie z alpha)")
    if np.all(dQ_darr > 0):
        print("     Q(alpha) MONOTONICZNIE ROSNIE z alpha (dla a_gam=0.040)")
    else:
        print("     Q(alpha) nie jest monotoniczne w tym zakresie")
print()
print("  3. Podsumowanie O-K1:")
print("     O-K1a: ODPOWIEDZIANE (P28): Q_min = 1.215 != 3/2")
print("     O-K1b: CZESCIOWO: r31/r21 -> 7+4sqrt(3) + NLO")
print("     O-K1c: ODPOWIEDZIANE (P28): delta nie jest stale")
print("     O-K1d (P29): Czy punkt optymalny min(chi^2) = punkt min|Q-3/2|?")

if not np.isnan(alpha_0_found):
    print(f"            alpha_0={alpha_0_found:.5f} vs alpha_opt=8.5616")
    frac = (alpha_0_found - 8.5616) / (8.62 - 8.5616) * 100
    print(f"            alpha_0 lezy na {alpha_0_found-8.5616:+.5f} od alpha_opt")
    if abs(r21_0_found - R21_LEPTON) < 0.5 and abs(r31_0_found - R31_LEPTON) < 5.0:
        print()
        print("  HIPOTEZA POTWIERDZONA: TGP daje Q=3/2 gdy parametry")
        print("  dokladnie reprodukuja masy leptonowe (nie TGP-optymalne).")
        print("  Innymi slowy: Q=3/2 jest KONSEKWENCJA DOKLADNOSCI DOPASOWANIA,")
        print("  nie dodatkowym zalozeniem TGP.")

print()

# ============================================================
# CZESC F: WYKRES
# ============================================================
print("CZESC F: Rysowanie wykresow")
print("-" * 70)

try:
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    fig.suptitle("P29: Punkt zerowy Q=3/2 w przestrzeni parametrow TGP", fontsize=13)

    # Panel 1: Q(alpha) dla a_gam=0.040
    ax1 = axes[0]
    if alpha_valid:
        ax1.plot(alpha_valid, Q_alpha, 'b-o', markersize=5)
        ax1.axhline(1.5, color='r', linestyle='--', label='Q=3/2')
        ax1.axvline(8.5616, color='g', linestyle=':', label='alpha_opt=8.5616', alpha=0.7)
        if not np.isnan(alpha_0_found):
            ax1.axvline(alpha_0_found, color='orange', linestyle='-.',
                       label=f'alpha_0={alpha_0_found:.4f}')
        ax1.set_xlabel('alpha')
        ax1.set_ylabel('Q')
        ax1.set_title('Q(alpha) dla a_gam=0.040, r31=3477')
        ax1.legend(fontsize=8)
        ax1.grid(True, alpha=0.3)

    # Panel 2: Q(a_gam) dla alpha=8.5616
    ax2 = axes[1]
    if agam_valid and Q_agam:
        ax2.plot(agam_valid, Q_agam, 'b-o', markersize=5)
        ax2.axhline(1.5, color='r', linestyle='--', label='Q=3/2')
        ax2.axvline(0.040, color='g', linestyle=':', label='a_gam_opt=0.040')
        ax2.set_xlabel('a_gam')
        ax2.set_ylabel('Q')
        ax2.set_title('Q(a_gam) dla alpha=8.5616, r31=3477')
        ax2.legend(fontsize=8)
        ax2.grid(True, alpha=0.3)

    # Panel 3: Mapa |Q-3/2| 2D
    ax3 = axes[2]
    valid_mask = ~np.isnan(Q_map_2d)
    if np.any(valid_mask):
        dQ_map = np.abs(Q_map_2d - 1.5)
        im = ax3.contourf(alpha_2d, agam_2d, dQ_map,
                          levels=20, cmap='RdYlGn_r')
        plt.colorbar(im, ax=ax3, label='|Q - 3/2|')
        # Koide zero line
        try:
            cs = ax3.contour(alpha_2d, agam_2d, Q_map_2d - 1.5,
                            levels=[0], colors='blue', linewidths=2)
            ax3.clabel(cs, fmt='Q=3/2', fontsize=9)
        except Exception:
            pass
        # Punkt optymalny
        ax3.scatter([8.5616], [0.040], color='white', marker='*', s=200,
                   zorder=10, label='punkt opt.')
        ax3.set_xlabel('alpha')
        ax3.set_ylabel('a_gam')
        ax3.set_title('Mapa |Q-3/2|(alpha, a_gam)')
        ax3.legend(fontsize=8)

    plt.tight_layout()
    plt.savefig('p29_koide_zero.png', dpi=120, bbox_inches='tight')
    print("  Wykres zapisany: p29_koide_zero.png")
except Exception as e:
    print(f"  UWAGA: Blad wykresu: {e}")

print()

# ============================================================
# PODSUMOWANIE P29
# ============================================================
print()
print("=" * 70)
print("PODSUMOWANIE P29: Punkt zerowy Q=3/2")
print("=" * 70)
print()

if not np.isnan(alpha_0_found):
    print(f"  alpha_0 = {alpha_0_found:.6f}  (gdzie Q(lambda*) = 3/2 dokladnie)")
    print(f"  lambda_0 = {lam_0_found:.6e}")
    print(f"  (r21_0, r31_0) = ({r21_0_found:.4f}, {r31_0_found:.4f})")
    print(f"  Leptony:         ({R21_LEPTON}, {R31_LEPTON})")
    print()
    delta_r21 = r21_0_found - R21_LEPTON
    delta_r31 = r31_0_found - R31_LEPTON
    print(f"  Odleglosc punktu Q=3/2 od mas leptonowych:")
    print(f"    delta_r21 = {delta_r21:+.4f}  ({delta_r21/R21_LEPTON*100:+.4f}%)")
    print(f"    delta_r31 = {delta_r31:+.4f}  ({delta_r31/R31_LEPTON*100:+.4f}%)")
    print()
    if abs(delta_r21) < 1.0 and abs(delta_r31) < 5.0:
        print("  STATUS: Bliskosc krzywej Koide wynika z bliskosci")
        print("  dokladnych mas leptonowych. Punkt TGP optymalny lezy na")
        print("  'skraju' krzywej Koide -- przy drobnej korekcie alpha")
        print("  TGP mogloby odtworzyc Q=3/2 dokladnie WRAZ z r21,r31 leptonu.")
    elif abs(delta_r21) < 3.0 and abs(delta_r31) < 20.0:
        print("  STATUS: Punkt Q=3/2 jest bliski masom leptonowym, ale")
        print("  nie identyczny. Dalsze badania potrzebne.")
    else:
        print("  STATUS: Punkt Q=3/2 jest daleko od mas leptonowych.")
        print("  Bliskosc krzywej Koide przy alpha_opt=8.5616 jest przypadkowa.")

print()
print("  NASTEPNE KROKI (P30):")
print("  - Dokladne dopasowanie TGP do (r21, r31) = (206.768, 3477.65)")
print("  - Czy dokladne dopasowanie wymaga (r21,r31) na krzywej Koide?")
print("  - Analiza: jak zmienia sie Q gdy parametry sa dokladnie dopasowane")
print("    do mas leptonowych zamiast do uproszczonego r31=3477?")
print()
print("P29 zakonczone.")
