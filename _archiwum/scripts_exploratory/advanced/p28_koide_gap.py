"""
p28_koide_gap.py
================
CEL: Systematyczna analiza "szczeliny Koide" w TGP.
     Dla kazdego punktu optymalnej rodziny parametrow:
       - znajdz lambda_Koide (Q=3/2 dokladnie)
       - porownaj z lambda* (punkt optymalny dopasowany do mas leptonow)
       - oblicz delta_Q = Q(lambda*) - 3/2
       - oblicz delta_r31 = r31_Koide - r31(lambda*)
     Analityczna formula na r31_Koide(r21).
     Skan Q(lambda) w okolicach lambda -> 0 (Q_min gdy K*3 ledwo istnieje).
     Mapa 2D: Q(r21, lambda) z krzywq Koide.

KONTEKST (Problem otwarty O-K1):
  P27 wykazal:
    (1) Q NIE jest zawsze >= 3/2; moze byc tak niskie jak ~1.22
    (2) Q(lambda) monotonicznie rosnie z lambda
    (3) lambda_Koide ~ 5.1e-6 (gdzie Q=3/2 dokladnie)
    (4) r31_Koide = 3477.65 vs TGP r31 = 3477.10 (roznica ~0.55)

  P28 badania:
    (A) Dokladna wartosc lambda_Koide i jej stosunek do lambda*
    (B) Analityczna formula r31_Koide(r21) -- weryfikacja
    (C) Zachowanie Q gdy lambda -> 0 (Q_min asymptotyczne)
    (D) Czy stosunek lambda_Koide/lambda* jest universalny?

ANALITYCZNA FORMULA (wyprowadzona w P27):
  Krzywa Koide Q=3/2: 2(1+sqrt(r21)+sqrt(r31))^2 = 3(1+r21+r31)
  Rozwiazanie dla r31:
    a = 1 + sqrt(r21)
    r31_Koide = (2a + sqrt(3 + 12*sqrt(r21) + 3*r21))^2 / 4
              = (2a + sqrt(6*a^2 - 3*(1-r21+r21) - 3*r21))^2 / 4
  Asymptotyki dla duzego r21:
    r31_Koide / r21 -> (2+sqrt(3))^2 = 7+4*sqrt(3) ~ 13.928   (wiodacy)
    NLO: + 4*(5+3*sqrt(3))/sqrt(r21)                 ~ +2.84/sqrt(r21) * r21
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

print("P28: Analiza szczeliny Koide -- Problem otwarty O-K1")
print("=" * 70)
print()

# ============================================================
# PODSTAWOWE NARZEDZIA (identyczne jak p27)
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

def find_zero(alpha, a_gam, lam, K_lo, K_hi):
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
    intervals = [
        (0.001, 0.050),
        (0.050, 0.500),
        (0.500, 5.0),
        (5.0,   50.0),
        (50.0, 200.0),
    ]
    zeros = []
    for K_lo, K_hi in intervals:
        K = find_zero(alpha, a_gam, lam, K_lo, K_hi)
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

# ============================================================
# CZESC A: ANALITYCZNA FORMULA r31_Koide(r21)
# ============================================================
print("CZESC A: Analityczna formula r31_Koide(r21)")
print("-" * 70)
print()
print("  Krzywa Koide: 2*(1+sqrt(r21)+sqrt(r31))^2 = 3*(1+r21+r31)")
print("  Podstawienie x = sqrt(r31), a = 1+sqrt(r21):")
print("    -x^2 + 4*a*x + (2*a^2 - 3 - 3*r21) = 0")
print("    x = 2*a + sqrt(2*a^2 + 3*r21 - 3)")
print()

def r31_koide_analytic(r21):
    """Dokladna r31 na krzywej Koide Q=3/2."""
    a = 1.0 + np.sqrt(r21)
    # Dyskryminant: 16*a^2 + 4*(2*a^2 - 3 - 3*r21) = 24*a^2 - 12 - 12*r21
    disc = 24*a**2 - 12 - 12*r21
    if disc < 0:
        return np.nan
    x = (4*a + np.sqrt(disc)) / 2.0   # sqrt(r31)
    return x**2

def r31_koide_asymptotic_LO(r21):
    """Wiodacy czlon asymptotyczny: r31/r21 -> (2+sqrt(3))^2"""
    return (2.0 + np.sqrt(3.0))**2 * r21

def r31_koide_asymptotic_NLO(r21):
    """NLO: r31/r21 -> (2+sqrt(3))^2 + 4*(5+3*sqrt(3))/sqrt(r21)"""
    LO  = (2.0 + np.sqrt(3.0))**2
    NLO = 4.0*(5.0 + 3.0*np.sqrt(3.0)) / np.sqrt(r21)
    return (LO + NLO) * r21

print(f"  (2+sqrt(3))^2 = {(2.0+np.sqrt(3.0))**2:.6f}  (wiodacy wspolczynnik)")
print(f"  NLO prefactor = {4.0*(5.0+3.0*np.sqrt(3.0)):.6f}")
print()

# Weryfikacja dla kilku r21
print(f"  {'r21':>8}  {'r31_exact':>12}  {'r31_LO':>12}  {'r31_NLO':>12}  {'Q(exact)':>10}")
print("  " + "-"*60)
for r21 in [9, 25, 50, 100, 207, 500, 1000]:
    r31_e = r31_koide_analytic(r21)
    r31_lo = r31_koide_asymptotic_LO(r21)
    r31_nlo = r31_koide_asymptotic_NLO(r21)
    Q_e = Q_from_ratios(r21, r31_e) if not np.isnan(r31_e) else np.nan
    print(f"  {r21:>8.0f}  {r31_e:>12.4f}  {r31_lo:>12.4f}  {r31_nlo:>12.4f}  {Q_e:>10.8f}")

print()
print("  OBSERWACJA: NLO jest dobrym przyblizeniem dla r21 >= 50.")
print(f"  Dla r21=207: r31_Koide = {r31_koide_analytic(207):.4f}")
print(f"               r31_NLO   = {r31_koide_asymptotic_NLO(207):.4f}")
print()

# ============================================================
# CZESC B: OPTYMALNA RODZINA PARAMETROW TGP
# ============================================================
print("CZESC B: Optymalna rodzina parametrow -- lambda_Koide vs lambda*")
print("-" * 70)
print()

# Optymalna rodzina parametrow (z poprzednich skryptow p21b, p22-p25)
# Cztery punkty z okolic minimum chi^2
optimal_family = [
    # (alpha, a_gam, lambda*, r21_TGP, r31_TGP, nazwa)
    (8.5616, 0.040, 5.501e-6, 207.0, 3477.1,  "glowny"),
    (8.50,   0.040, 5.45e-6,  205.0, 3460.0,  "alphaL"),
    (8.62,   0.040, 5.55e-6,  209.0, 3494.0,  "alphaH"),
    (8.5616, 0.042, 5.80e-6,  208.0, 3490.0,  "agamH"),
]

print("  Dla kazdego punktu rodziny szukam lambda_Koide:")
print("  (lambda takie ze Q(K1,K2,K3) = 3/2 dokladnie)")
print()

def Q_at_lambda(alpha, a_gam, lam):
    """Oblicz Q dla podanych parametrow. Zwraca (Q, r21, r31) lub (nan,nan,nan)."""
    zeros = find_all_zeros(alpha, a_gam, lam)
    if len(zeros) < 3:
        return np.nan, np.nan, np.nan
    K1, K2, K3 = zeros[0], zeros[1], zeros[2]
    Q  = koide_Q(K1, K2, K3)
    r21 = K2/K1
    r31 = K3/K1
    return Q, r21, r31

def Q_minus_3_2(lam, alpha, a_gam):
    Q, _, _ = Q_at_lambda(alpha, a_gam, lam)
    if np.isnan(Q):
        return np.nan
    return Q - 1.5

results = []

for alpha, agam, lam_star, r21_star, r31_star, name in optimal_family:
    print(f"  [{name}] alpha={alpha}, a_gam={agam}, lambda*={lam_star:.3e}")

    # Sprawdz Q przy lambda*
    Q_star, r21_s, r31_s = Q_at_lambda(alpha, agam, lam_star)
    print(f"    Q(lambda*) = {Q_star:.8f},  delta_Q = {Q_star-1.5:+.8f}")
    print(f"    r21={r21_s:.3f}, r31={r31_s:.3f}")

    # Szukaj lambda_Koide metodq brentq
    # Wiemy ze Q rosnie z lambda (P27), wiec:
    # - jezeli Q(lambda*) < 3/2: lambda_Koide > lambda*
    # - jezeli Q(lambda*) > 3/2: lambda_Koide < lambda*
    lam_koide = np.nan
    Q_lo_val = Q_minus_3_2(lam_star * 0.5, alpha, agam)
    Q_hi_val = Q_minus_3_2(lam_star * 2.0, alpha, agam)

    if np.isfinite(Q_lo_val) and np.isfinite(Q_hi_val) and Q_lo_val * Q_hi_val < 0:
        try:
            lam_koide = brentq(
                Q_minus_3_2, lam_star*0.5, lam_star*2.0,
                args=(alpha, agam), xtol=lam_star*1e-6, maxiter=50
            )
        except Exception as e:
            print(f"    UWAGA: brentq nieudany: {e}")
    else:
        # Probuj szerszy zakres
        for factor_lo, factor_hi in [(0.1, 5.0), (0.01, 10.0)]:
            Q_lo_val = Q_minus_3_2(lam_star * factor_lo, alpha, agam)
            Q_hi_val = Q_minus_3_2(lam_star * factor_hi, alpha, agam)
            if np.isfinite(Q_lo_val) and np.isfinite(Q_hi_val) and Q_lo_val * Q_hi_val < 0:
                try:
                    lam_koide = brentq(
                        Q_minus_3_2, lam_star*factor_lo, lam_star*factor_hi,
                        args=(alpha, agam), xtol=lam_star*1e-6, maxiter=50
                    )
                    break
                except Exception:
                    pass

    if not np.isnan(lam_koide):
        Q_k, r21_k, r31_k = Q_at_lambda(alpha, agam, lam_koide)
        r31_analytic = r31_koide_analytic(r21_k) if not np.isnan(r21_k) else np.nan
        ratio = lam_koide / lam_star
        delta_r31 = r31_k - r31_s
        print(f"    lambda_Koide = {lam_koide:.6e}")
        print(f"    lambda_Koide / lambda* = {ratio:.6f}")
        print(f"    Q(lambda_Koide) = {Q_k:.10f}  (powinno byc 1.5)")
        print(f"    r21_Koide = {r21_k:.4f},  r31_Koide = {r31_k:.4f}")
        print(f"    r31_analytic(r21_Koide) = {r31_analytic:.4f}")
        print(f"    delta_r31 = r31_Koide - r31* = {delta_r31:+.4f}")
        results.append({
            'name': name, 'alpha': alpha, 'agam': agam,
            'lam_star': lam_star, 'lam_koide': lam_koide,
            'ratio': ratio,
            'Q_star': Q_star, 'delta_Q': Q_star - 1.5,
            'r21_star': r21_s, 'r31_star': r31_s,
            'r21_koide': r21_k, 'r31_koide': r31_k,
            'delta_r31': delta_r31,
        })
    else:
        print(f"    OSTRZEZENIE: Nie znaleziono lambda_Koide!")
        results.append({
            'name': name, 'alpha': alpha, 'agam': agam,
            'lam_star': lam_star, 'lam_koide': np.nan,
            'ratio': np.nan, 'Q_star': Q_star,
            'delta_Q': Q_star - 1.5 if not np.isnan(Q_star) else np.nan,
            'r21_star': r21_s, 'r31_star': r31_s,
            'r21_koide': np.nan, 'r31_koide': np.nan, 'delta_r31': np.nan,
        })
    print()

# Podsumowanie rodziny
print()
print("  PODSUMOWANIE RODZINY:")
print(f"  {'nazwa':>8}  {'lam_K/lam*':>12}  {'delta_Q':>12}  {'delta_r31':>12}")
print("  " + "-"*50)
for r in results:
    print(f"  {r['name']:>8}  {r['ratio']:>12.6f}  {r['delta_Q']:>+12.6f}  {r['delta_r31']:>+12.4f}")

print()

# ============================================================
# CZESC C: ZACHOWANIE Q gdy lambda -> 0 (Q_min)
# ============================================================
print("CZESC C: Zachowanie Q gdy lambda -> 0")
print("-" * 70)
print()
print("  Pytanie: co sie dzieje z Q gdy lambda -> 0?")
print("  Czy istnieje Q_min > 1 gdy K*3 ledwo egzystuje?")
print()

alpha_c = 8.5616
agam_c  = 0.040

# Skan ku malym lambda
lam_small = np.logspace(-8, -5.5, 30)
Q_small   = []
r21_small = []
r31_small = []
lam_valid_small = []

print("  Skan lambda od 1e-8 do 3e-6:")
for lam in lam_small:
    zeros = find_all_zeros(alpha_c, agam_c, lam)
    if len(zeros) >= 3:
        K1, K2, K3 = zeros[0], zeros[1], zeros[2]
        Q_val  = koide_Q(K1, K2, K3)
        r21_v  = K2/K1
        r31_v  = K3/K1
        Q_small.append(Q_val)
        r21_small.append(r21_v)
        r31_small.append(r31_v)
        lam_valid_small.append(lam)

if Q_small:
    Q_arr  = np.array(Q_small)
    r21_arr = np.array(r21_small)
    r31_arr = np.array(r31_small)
    lam_arr = np.array(lam_valid_small)

    print()
    print(f"  {'lambda':>12}  {'Q':>10}  {'r21':>10}  {'r31':>12}  {'Q-3/2':>12}")
    print("  " + "-"*60)
    for i in range(0, len(lam_arr), max(1, len(lam_arr)//10)):
        print(f"  {lam_arr[i]:>12.3e}  {Q_arr[i]:>10.6f}  "
              f"{r21_arr[i]:>10.3f}  {r31_arr[i]:>12.3f}  "
              f"{Q_arr[i]-1.5:>+12.6f}")

    Q_min_idx = np.argmin(Q_arr)
    print()
    print(f"  Q_min = {Q_arr[Q_min_idx]:.8f}  przy lambda = {lam_arr[Q_min_idx]:.3e}")
    print(f"  r21 przy Q_min = {r21_arr[Q_min_idx]:.4f}")
    print(f"  r31 przy Q_min = {r31_arr[Q_min_idx]:.4f}")
    print()

    # Szacunek lambda progu (gdzie K*3 znika)
    print("  Szacunek lambda_prog (gdzie K*3 przestaje istniec):")
    lam_test_lo = lam_arr[0]
    # Sprobuj jeszcze mniejsze lambda
    lam_prog_est = lam_arr[0]
    for lam_t in np.logspace(-9, -8, 10):
        zeros_t = find_all_zeros(alpha_c, agam_c, lam_t)
        if len(zeros_t) >= 3:
            lam_prog_est = lam_t
        else:
            break
    print(f"  lambda_prog ~ {lam_prog_est:.3e} (K*3 ostatni raz wykryty)")
    print()
else:
    print("  UWAGA: Brak 3 zer dla lambda < 3e-6 (K*3 nie istnieje)")
    print()

# ============================================================
# CZESC D: MAPA 2D Q(r21, lambda) z krzywq Koide
# ============================================================
print("CZESC D: Mapa struktury Q w przestrzeni (r21, lambda)")
print("-" * 70)
print()
print("  Pytanie: Czy Q=3/2 tworzy izoline w przestrzeni parametrow?")
print("  Buduje mape Q(r21, lambda) dla alpha=8.5616, a_gam=0.040")
print()

alpha_d = 8.5616
agam_d  = 0.040

lam_map = np.logspace(-6.5, -5.3, 20)
Q_map   = []
r21_map = []
r31_map = []
lam_map_valid = []

for lam in lam_map:
    zeros = find_all_zeros(alpha_d, agam_d, lam)
    if len(zeros) >= 3:
        K1, K2, K3 = zeros[0], zeros[1], zeros[2]
        Q_map.append(koide_Q(K1, K2, K3))
        r21_map.append(K2/K1)
        r31_map.append(K3/K1)
        lam_map_valid.append(lam)

if Q_map:
    Q_m   = np.array(Q_map)
    r21_m = np.array(r21_map)
    r31_m = np.array(r31_map)
    lam_m = np.array(lam_map_valid)

    print(f"  {'lambda':>12}  {'r21':>10}  {'r31':>12}  {'Q':>10}  {'Q-3/2':>12}")
    print("  " + "-"*62)
    for i in range(len(lam_m)):
        marker = " <-- Q=3/2 przekroczone" if abs(Q_m[i]-1.5) < 0.01 else ""
        print(f"  {lam_m[i]:>12.3e}  {r21_m[i]:>10.3f}  {r31_m[i]:>12.3f}  "
              f"{Q_m[i]:>10.6f}  {Q_m[i]-1.5:>+12.6f}{marker}")

    # r31_Koide na tej siatce r21
    r31_koide_line = np.array([r31_koide_analytic(r) for r in r21_m])

    print()
    print("  Porownanie r31_TGP vs r31_Koide (na tej samej siatce r21):")
    print(f"  {'lambda':>12}  {'r21':>10}  {'r31_TGP':>12}  {'r31_Koide':>12}  {'delta':>10}")
    print("  " + "-"*60)
    for i in range(len(lam_m)):
        delta = r31_koide_line[i] - r31_m[i]
        print(f"  {lam_m[i]:>12.3e}  {r21_m[i]:>10.3f}  {r31_m[i]:>12.3f}  "
              f"{r31_koide_line[i]:>12.3f}  {delta:>+10.3f}")

    print()
    print("  KLUCZOWE PYTANIE O-K1c:")
    print("  Czy delta = r31_Koide - r31_TGP jest stale (niezaleznie od lambda)?")
    deltas = r31_koide_line - r31_m
    if len(deltas) > 3:
        print(f"  delta_min = {np.min(deltas):.4f}")
        print(f"  delta_max = {np.max(deltas):.4f}")
        print(f"  delta_mean = {np.mean(deltas):.4f}")
        print(f"  delta_std  = {np.std(deltas):.4f}")
        if np.std(deltas) < 5:
            print("  -> delta jest PRAWIE STALE! (~niezmiennicze wzgledem lambda)")
        else:
            print("  -> delta zmienia sie z lambda (nie jest stale)")
    print()

# ============================================================
# CZESC E: POROWNANIE STOSUNEK r31/r21 wzdluz trajektorii TGP
# ============================================================
print("CZESC E: Ewolucja stosunku r31/r21 z lambda")
print("-" * 70)
print()
print("  Asymptotyczny stosunek Koide: (2+sqrt(3))^2 = 7+4*sqrt(3) =",
      f"{7+4*np.sqrt(3):.6f}")
print("  Ale dla r21~207 korekta NLO daje ~16.77")
print()

if Q_map:
    ratio_TGP   = r31_m / r21_m
    ratio_Koide = r31_koide_line / r21_m

    print(f"  {'lambda':>12}  {'r31/r21 TGP':>14}  {'r31/r21 Koide':>15}  {'roznica':>10}")
    print("  " + "-"*56)
    for i in range(len(lam_m)):
        diff = ratio_Koide[i] - ratio_TGP[i]
        print(f"  {lam_m[i]:>12.3e}  {ratio_TGP[i]:>14.6f}  {ratio_Koide[i]:>15.6f}  {diff:>+10.6f}")

    print()
    print("  Sprawdzam asymptotyki NLO:")
    for r21_test in [100, 207, 500, 1000, 5000]:
        nlo = r31_koide_asymptotic_NLO(r21_test) / r21_test
        lo  = (2+np.sqrt(3))**2
        print(f"    r21={r21_test:5.0f}: LO={lo:.5f}, NLO wspolcz={nlo:.5f}, "
              f"delta_NLO = +{4*(5+3*np.sqrt(3))/np.sqrt(r21_test):.5f}")

print()

# ============================================================
# CZESC F: WYKRES (opcjonalny)
# ============================================================
print("CZESC F: Rysowanie wykresu Q(lambda)")
print("-" * 70)

try:
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle("P28: Analiza szczeliny Koide -- TGP vs Q=3/2", fontsize=13)

    # Panel 1: Q(lambda) z zaznaczeniem Q=3/2
    ax1 = axes[0]
    if lam_map_valid and Q_map:
        ax1.semilogx(lam_map_valid, Q_map, 'b-o', markersize=4)
        ax1.axhline(1.5, color='r', linestyle='--', label='Q=3/2 (Koide)')
        # Zaznacz lambda*
        ax1.axvline(5.501e-6, color='g', linestyle=':', label='lambda*=5.501e-6')
        ax1.set_xlabel('lambda')
        ax1.set_ylabel('Q')
        ax1.set_title('Q(lambda) dla alpha=8.5616, a_gam=0.040')
        ax1.legend(fontsize=8)
        ax1.grid(True, alpha=0.3)

    # Panel 2: krzywa Koide w (r21, r31)
    ax2 = axes[1]
    r21_line = np.linspace(4, 1200, 500)
    r31_line = np.array([r31_koide_analytic(r) for r in r21_line])
    mask = ~np.isnan(r31_line)
    ax2.plot(r21_line[mask], r31_line[mask], 'r-', label='Krzywa Koide Q=3/2')
    if lam_map_valid:
        sc = ax2.scatter(r21_m, r31_m, c=np.log10(lam_m), cmap='viridis',
                        s=30, zorder=5, label='Trajektoria TGP')
        plt.colorbar(sc, ax=ax2, label='log10(lambda)')
    # Punkt leptonowy
    ax2.scatter([206.768], [3477.65], color='red', marker='*', s=200,
                zorder=10, label='Leptony (exp.)')
    ax2.set_xlabel('r21 = K*2/K*1')
    ax2.set_ylabel('r31 = K*3/K*1')
    ax2.set_title('Trajektoria TGP vs krzywa Koide')
    ax2.legend(fontsize=8)
    ax2.grid(True, alpha=0.3)

    # Panel 3: delta = r31_Koide - r31_TGP vs lambda
    ax3 = axes[2]
    if lam_map_valid and len(r31_koide_line) > 0:
        deltas_plot = r31_koide_line - r31_m
        ax3.semilogx(lam_m, deltas_plot, 'b-o', markersize=5)
        ax3.axhline(0, color='r', linestyle='--')
        ax3.set_xlabel('lambda')
        ax3.set_ylabel('r31_Koide - r31_TGP')
        ax3.set_title('Szczelina Koide (delta_r31)')
        ax3.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('p28_koide_gap.png', dpi=120, bbox_inches='tight')
    print("  Wykres zapisany: p28_koide_gap.png")
except Exception as e:
    print(f"  UWAGA: Blad wykresu: {e}")

print()

# ============================================================
# PODSUMOWANIE P28
# ============================================================
print()
print("=" * 70)
print("PODSUMOWANIE P28: Szczelina Koide")
print("=" * 70)
print()
print("1. FORMULA ANALITYCZNA r31_Koide(r21):")
print("     r31_K = (2a + sqrt(6a^2 - 3 - 3*r21))^2/4,  a=1+sqrt(r21)")
print(f"   Dla r21=207: r31_K = {r31_koide_analytic(207):.4f}")
print(f"   Asymptotyka LO: r31/r21 -> {(2+np.sqrt(3))**2:.6f} = 7+4*sqrt(3)")
print()
print("2. SZCZELINA KOIDE:")
if results:
    r = results[0]  # glowny punkt
    print(f"   Punkt glowny: lambda* = {r['lam_star']:.3e}")
    print(f"   Q(lambda*)  = {r['Q_star']:.8f}")
    print(f"   delta_Q     = Q - 3/2 = {r['delta_Q']:+.8f}")
    if not np.isnan(r.get('lam_koide', np.nan)):
        print(f"   lambda_Koide = {r['lam_koide']:.6e}")
        print(f"   lambda_K/lambda* = {r['ratio']:.6f}")
        print(f"   delta_r31    = r31_K - r31* = {r['delta_r31']:+.4f}")
print()
print("3. STATUS PROBLEMU O-K1:")
print("   O-K1a: Q NIE jest zawsze >= 3/2 (P27). Q min ~ 1.22.")
print("   O-K1b: r31/r21 -> 7+4*sqrt(3)=13.928 + NLO (P27, P28).")
print("   O-K1c: Czy delta_r31 jest STALE (niezmiennicze wzgl. lambda)?")
print("          -> zbadano w Czesci D powyzej.")
print()
print("4. NASTEPNE KROKI (P29):")
print("   - Analiza analityczna: skad pochodzi szczelina delta_r31 ~ 0.55?")
print("   - Czy istnieje symetria V_mod wymuszajqca r31_K/r21 = 7+4*sqrt(3)?")
print("   - Czyzba energia bifurkacji przy lambda_Koide jest specjalna?")
print("   - Zwiazek z fizycznq interpretacjq: generacja leptonow = K*1, K*2, K*3")
print()
print("P28 zakonczone.")
