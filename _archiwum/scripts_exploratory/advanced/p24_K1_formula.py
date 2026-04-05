"""
p24_K1_formula.py
==================
CEL: Analityczna formula na K*1(alpha, a_gam) z teorii perturbacji.

HIPOTEZA (z obserwacji P23):
  K*1 ≈ C_K * a_gam / (1 + alpha),   C_K ≈ 2.35

UZASADNIENIE ANALITYCZNE:
  Warunek samospójności: g(K) = E(K)/(4*pi*K) - 1 = 0

  Dla malego K (K << 1), profil Yukawa: phi = 1 + K*exp(-r)/r

  Energia kinetyczna (dominant term):
    E_k = 2*pi*K^2*(1+alpha) * I_k(a_gam)
    gdzie I_k(a) = integral_a^inf e^{-2r}*(r+1)^2/r^2 dr

  Energia potencjalna:
    E_p ≈ -2*pi*K^2 * I_p(a_gam)
    gdzie I_p(a) = integral_a^inf e^{-2r} dr = e^{-2a}/2

  Warunek E = 4*pi*K:
    K*1 ≈ 2 / [(1+alpha)*I_k(a) - I_p(a)]

KROKI:
  1. Oblicz numerycznie K*1 dla siatki (alpha, a_gam)
  2. Sprawdz skalowanie K*1 ∝ a_gam/(1+alpha)
  3. Oblicz C_K numerycznie i analitycznie
  4. Porownaj z wynikami P21b/P23
"""

import numpy as np
from scipy.optimize import brentq
from scipy.integrate import quad
import warnings
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')

# ============================================================
# STALE
# ============================================================
R_MAX = 60.0
GAMMA = 1.0

print("P24: Analityczna formula K*1(alpha, a_gam)")
print("=" * 65)
print()

# ============================================================
# DEFINICJE
# ============================================================
def V_mod(phi, lam):
    return GAMMA/3*phi**3 - GAMMA/4*phi**4 + lam/6*(phi-1)**6

def energy_log(K, alpha, a_gam, lam, N=1500):
    """Energia Yukawa na siatce log (identyczna z p21b)."""
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

def find_K1(alpha, a_gam, lam, K_lo=0.003, K_hi=0.025):
    """Znajdz K*1 metoda brentq."""
    try:
        g_lo = g_func(K_lo, alpha, a_gam, lam)
        g_hi = g_func(K_hi, alpha, a_gam, lam)
        if not (np.isfinite(g_lo) and np.isfinite(g_hi)):
            return np.nan
        if g_lo * g_hi > 0:
            # Szersza okolica
            for K_lo2, K_hi2 in [(0.001, 0.030), (0.001, 0.050)]:
                g1 = g_func(K_lo2, alpha, a_gam, lam)
                g2 = g_func(K_hi2, alpha, a_gam, lam)
                if np.isfinite(g1) and np.isfinite(g2) and g1*g2 < 0:
                    return brentq(lambda K: g_func(K, alpha, a_gam, lam),
                                  K_lo2, K_hi2, xtol=1e-8)
            return np.nan
        return brentq(lambda K: g_func(K, alpha, a_gam, lam),
                      K_lo, K_hi, xtol=1e-8)
    except Exception:
        return np.nan

# ============================================================
# ANALITYCZNE PRZYBLIŻENIE PIERWSZEGO RZĘDU
# ============================================================
def I_k_analytic(a):
    """
    I_k(a) = integral_a^inf e^{-2r}*(r+1)^2/r^2 dr
    = integral e^{-2r} dr + 2*integral e^{-2r}/r dr + integral e^{-2r}/r^2 dr
    """
    # Term 1: e^{-2a}/2
    I1 = np.exp(-2*a)/2.0
    # Term 2: 2 * E1(2a)  gdzie E1 = integral_x^inf e^{-t}/t dt
    # Przybliżenie E1(x) = exp(-x)*(-1/x + 1/x^2 - 2/x^3 ...) + ln(x) + gamma
    # Dla malych a: E1(2a) ≈ -ln(2a) - gamma + 2a - (2a)^2/4 + ...
    from scipy.special import exp1
    I2 = 2.0 * exp1(2*a)
    # Term 3: integral_a^inf e^{-2r}/r^2 dr ≈ e^{-2a}/a - 2*E1(2a)
    # (przez cześci: integral e^{-2r}/r^2 dr = -e^{-2r}/r + 2*integral e^{-2r}/r dr)
    I3 = np.exp(-2*a)/a - 2.0*exp1(2*a)
    return I1 + I2 + I3

def I_p_analytic(a):
    """I_p(a) = integral_a^inf e^{-2r} dr = e^{-2a}/2."""
    return np.exp(-2*a)/2.0

def K1_perturbative(alpha, a_gam):
    """
    Przybliżenie K*1 z perturbacji pierwszego rzędu:
    K*1 ≈ 2 / [(1+alpha)*I_k(a) - I_p(a)]
    """
    Ik = I_k_analytic(a_gam)
    Ip = I_p_analytic(a_gam)
    denom = (1+alpha)*Ik - Ip
    if denom <= 0:
        return np.nan
    return 2.0 / denom

# ============================================================
# KROK 1: Weryfikacja hipotezy K*1 ∝ a_gam/(1+alpha)
# ============================================================
print("KROK 1: Weryfikacja hipotezy K*1 ≈ C_K * a_gam/(1+alpha)")
print("-" * 65)

# Punkty z P23 (poprawna lambda*)
family_data = [
    (5.9148, 0.025, 2.883e-6),
    (6.8675, 0.030, 3.743e-6),
    (7.7449, 0.035, 4.618e-6),
    (8.5616, 0.040, 5.501e-6),
]

print(f"  {'alpha':>8} {'a_gam':>7}  {'K1_num':>10}  {'K1_pert':>10}  {'ratio':>8}  {'C_K':>8}")
print("  " + "-"*60)

C_K_vals = []
for alpha, agam, lam in family_data:
    K1_num  = find_K1(alpha, agam, lam)
    K1_pert = K1_perturbative(alpha, agam)
    ratio   = K1_num / K1_pert if (not np.isnan(K1_num) and K1_pert > 0) else np.nan
    # C_K z definicji K*1 = C_K * a_gam/(1+alpha)
    C_K = K1_num * (1+alpha) / agam if not np.isnan(K1_num) else np.nan
    C_K_vals.append(C_K)
    print(f"  {alpha:>8.4f} {agam:>7.4f}  {K1_num:>10.7f}  {K1_pert:>10.7f}  {ratio:>8.5f}  {C_K:>8.5f}")

C_K_vals = [c for c in C_K_vals if not np.isnan(c)]
C_K_mean = np.mean(C_K_vals)
C_K_std  = np.std(C_K_vals)
print(f"\n  C_K = {C_K_mean:.5f} ± {C_K_std:.5f}  (std/mean = {C_K_std/C_K_mean*100:.3f}%)")
print()

# ============================================================
# KROK 2: Analityczna wartosc C_K
# ============================================================
print("KROK 2: Analityczna wartosc C_K")
print("-" * 65)

# C_K = K*1 * (1+alpha) / a_gam = 2*(1+alpha)/a_gam / [(1+alpha)*I_k - I_p]
# = 2 / [a_gam*(I_k - I_p/(1+alpha))]
# Dla dominujacego czlonu (I_k >> I_p/(1+alpha) dla malego a_gam):
# C_K ≈ 2 / [a_gam * I_k(a_gam)]

# Oblicz C_K_analytic dla a_gam w zakresie
a_test = [0.025, 0.030, 0.035, 0.040]
print(f"  {'a_gam':>8}  {'I_k':>12}  {'I_p':>12}  {'K1_pert':>12}  {'C_K_pert':>12}")
print("  " + "-"*65)
for a, (alpha, agam, lam) in zip(a_test, family_data):
    Ik = I_k_analytic(a)
    Ip = I_p_analytic(a)
    K1p = K1_perturbative(alpha, a)
    CK_pert = K1p * (1+alpha) / a
    print(f"  {a:>8.4f}  {Ik:>12.4f}  {Ip:>12.6f}  {K1p:>12.7f}  {CK_pert:>12.5f}")
print()
print(f"  C_K numeryczny: {C_K_mean:.5f}")
print(f"  C_K z perturbacji (srednia): {np.mean([K1_perturbative(a, ag)*(1+a)/ag for a, ag, _ in family_data]):.5f}")
print()

# Stosunek C_K_num / C_K_pert (korekcja wyzszego rzedu)
ratio_mean = np.mean([c / (K1_perturbative(a, ag)*(1+a)/ag)
                       for c, (a, ag, _) in zip(C_K_vals, family_data)])
print(f"  Korekcja wyzszego rzedu: C_K_num/C_K_pert = {ratio_mean:.5f}")
print()

# ============================================================
# KROK 3: Skan K*1(alpha, a_gam) - szeroka siatka
# ============================================================
print("KROK 3: Skan K*1 na siatce (alpha, a_gam) -- weryfikacja skalowania")
print("-" * 65)

# Lam z interpolacji rodziny (lub z reg. empirycznej)
from scipy.interpolate import interp1d
FAM_ALPHA = np.array([5.9148, 6.8675, 7.7449, 8.5616])
FAM_LAM   = np.array([2.883e-6, 3.743e-6, 4.618e-6, 5.501e-6])
FAM_AGAM  = np.array([0.025, 0.030, 0.035, 0.040])
f_lam  = interp1d(FAM_ALPHA, FAM_LAM,  kind='linear', fill_value='extrapolate')
f_agam = interp1d(FAM_ALPHA, FAM_AGAM, kind='linear', fill_value='extrapolate')

alpha_vals = np.array([5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 15.0])
print(f"  {'alpha':>8} {'a_gam':>8} {'K1_num':>12} {'K1_pred':>12} {'C_K':>8} {'err%':>8}")
print("  " + "-"*65)
C_K_grid = []
for alpha in alpha_vals:
    agam = float(f_agam(alpha))
    lam  = float(f_lam(alpha))
    K1_num = find_K1(alpha, agam, lam)
    C_K_val = K1_num*(1+alpha)/agam if not np.isnan(K1_num) else np.nan
    K1_pred = C_K_mean * agam / (1+alpha)
    err = (K1_num - K1_pred)/K1_pred*100 if not np.isnan(K1_num) else np.nan
    C_K_grid.append(C_K_val)
    print(f"  {alpha:>8.1f} {agam:>8.4f} {K1_num:>12.7f} {K1_pred:>12.7f} {C_K_val:>8.5f} {err:>8.3f}%")
print()

# Sprawdz czy C_K jest stale takze poza rodzina
# Zmieniamy tylko alpha przy stalym a_gam=0.040
print("  Skan alpha przy stalym a_gam=0.040, lam=5.501e-6:")
print(f"  {'alpha':>8} {'K1_num':>12}  {'C_K':>8}  {'K1_pred':>12} {'err%':>8}")
print("  " + "-"*60)
alpha_fix_scan = [4.0, 5.0, 6.0, 7.0, 8.0, 8.5616, 9.0, 10.0, 12.0]
for alpha in alpha_fix_scan:
    agam, lam = 0.040, 5.501e-6
    K1_num = find_K1(alpha, agam, lam)
    C_K_val = K1_num*(1+alpha)/agam if not np.isnan(K1_num) else np.nan
    K1_pred = C_K_mean * agam / (1+alpha)
    err = (K1_num - K1_pred)/K1_pred*100 if not np.isnan(K1_num) else np.nan
    print(f"  {alpha:>8.1f} {K1_num:>12.7f}  {C_K_val:>8.5f}  {K1_pred:>12.7f} {err:>8.3f}%")
print()

# Skan a_gam przy stalym alpha=8.5616
print("  Skan a_gam przy stalym alpha=8.5616, lam=5.501e-6:")
print(f"  {'a_gam':>8} {'K1_num':>12}  {'C_K':>8}  {'K1_pred':>12} {'err%':>8}")
print("  " + "-"*60)
agam_scan = [0.010, 0.020, 0.030, 0.040, 0.050, 0.060, 0.080, 0.100]
for agam in agam_scan:
    alpha, lam = 8.5616, 5.501e-6
    K1_num = find_K1(alpha, agam, lam, K_lo=0.001, K_hi=0.050)
    C_K_val = K1_num*(1+alpha)/agam if not np.isnan(K1_num) else np.nan
    K1_pred = C_K_mean * agam / (1+alpha)
    err = (K1_num - K1_pred)/K1_pred*100 if not np.isnan(K1_num) else np.nan
    print(f"  {agam:>8.4f} {K1_num:>12.7f}  {C_K_val:>8.5f}  {K1_pred:>12.7f} {err:>8.3f}%")
print()

# ============================================================
# KROK 4: Analityczne wyjaśnienie C_K
# ============================================================
print("KROK 4: Analityczne wyznaczenie C_K")
print("-" * 65)

# K*1 * (1+alpha) / a_gam = C_K
# Warunek samospójności: K = 2 / [(1+alpha)*I_k(a) - I_p(a)]
# K*(1+alpha)/a = 2*(1+alpha) / {a * [(1+alpha)*I_k(a) - I_p(a)]}
# = 2 / {a * [I_k(a) - I_p(a)/(1+alpha)]}
# ≈ 2 / [a * I_k(a)]  dla (1+alpha) >> 1

print("  C_K = K*1*(1+alpha)/a_gam ≈ 2 / [a_gam * I_k(a_gam)]")
print()
for agam in [0.025, 0.030, 0.035, 0.040]:
    Ik = I_k_analytic(agam)
    CK_approx = 2.0 / (agam * Ik)
    print(f"  a_gam={agam:.3f}: I_k={Ik:.4f}, C_K_analytic = 2/(a*I_k) = {CK_approx:.5f}")

# Asymptotyka dla malego a:
# I_k(a) ≈ e^{-2a}/a + (2E1(2a) - 2E1(2a)) + e^{-2a}/2
# Dominujacy czlon: I_k(a) ≈ 1/a  dla a << 1
# Wiec: C_K ≈ 2 / [a * 1/a] = 2
# Korekcja: I_k(a) = 1/a - 2*ln(a) - gamma - 1/2 + ...
# C_K ≈ 2 / [1 - a*(2*ln(a) + gamma + 1/2) + ...]

print()
print("  Asymptotyka dla malego a_gam:")
print("    I_k(a) ≈ 1/a - 2*ln(a) - 0.577 + O(a)")
for agam in [0.025, 0.030, 0.035, 0.040]:
    gamma_E = 0.5772
    Ik_asym = 1/agam - 2*np.log(agam) - gamma_E
    CK_asym = 2.0 / (agam * Ik_asym)
    Ik_exact = I_k_analytic(agam)
    print(f"  a_gam={agam:.3f}: I_k(asymp)={Ik_asym:.4f}, I_k(exact)={Ik_exact:.4f}, C_K_asym={CK_asym:.5f}")
print()

# ============================================================
# KROK 5: Podsumowanie -- formula K*1
# ============================================================
print("KROK 5: Ostateczna formula K*1")
print("=" * 65)
print()
print("  FORMULA PERTURBACYJNA:")
print()
print("    K*1(alpha, a_gam) ≈ C_K * a_gam / (1 + alpha)")
print()
print(f"    C_K = {C_K_mean:.5f} ± {C_K_std:.5f}  (numerycznie)")
print()
print("  WYPROWADZENIE:")
print("    Warunek samospójności: E(K*1)/(4*pi*K*1) = 1")
print("    Dla malego K i profilu Yukawa phi = 1 + K*exp(-r)/r:")
print("      E_k ≈ 2*pi*(1+alpha)*K^2 * I_k(a_gam)")
print("      I_k(a) = e^{-2a}/a + ... ≈ 1/a  dla a << 1")
print("    => E_k ≈ 2*pi*(1+alpha)*K^2/a_gam")
print("    => g(K) = K*(1+alpha)/(2*a_gam) - 1 = 0")
print("    => K*1_LO = 2*a_gam/(1+alpha)  [L.O. = leading order]")
print()
print(f"    K*1_LO(a=0.040, alpha=8.5616) = {2*0.040/9.5616:.7f}")
print(f"    K*1_num(a=0.040, alpha=8.5616) = 0.0098200  (z p21b)")
print(f"    Korekcja: {0.009820 / (2*0.040/9.5616):.5f}  = C_K/2")
print()
print(f"    Pelna formula: K*1 ≈ {C_K_mean:.4f} * a_gam/(1+alpha)")
print()

# Sprawdzenie dla czterech punktow rodziny
print("  Weryfikacja na rodzinie optymalnej:")
print(f"  {'alpha':>8} {'a_gam':>7}  {'K1_num':>10}  {'K1_form':>10}  {'err%':>8}")
for alpha, agam, lam in family_data:
    K1_num  = find_K1(alpha, agam, lam)
    K1_form = C_K_mean * agam / (1+alpha)
    err = (K1_num-K1_form)/K1_form*100
    print(f"  {alpha:>8.4f} {agam:>7.4f}  {K1_num:>10.7f}  {K1_form:>10.7f}  {err:>8.4f}%")
print()

# ============================================================
# WYKRES
# ============================================================
print("Tworze wykres...")
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
fig.suptitle(f'P24: Analityczna formula K*1 ≈ {C_K_mean:.4f} * a_gam/(1+alpha)\n'
             'Weryfikacja skalowania z teorii perturbacji', fontsize=11, fontweight='bold')

# Panel 1: K*1 vs a_gam/(1+alpha) -- skalowanie
ax = axes[0]
ratio_vals_x = [agam/(1+alpha) for alpha, agam, _ in family_data]
K1_num_vals  = [find_K1(alpha, agam, lam) for alpha, agam, lam in family_data]
ax.scatter(ratio_vals_x, K1_num_vals, color='red', s=120, zorder=5, label='TGP (rodzina)')
x_lin = np.linspace(0, max(ratio_vals_x)*1.3, 50)
ax.plot(x_lin, C_K_mean*x_lin, 'b-', lw=2, label=f'K*1 = {C_K_mean:.4f} * a_gam/(1+alpha)')
ax.set_xlabel('a_gam/(1+alpha)'); ax.set_ylabel('K*1')
ax.set_title('Skalowanie K*1', fontsize=9)
ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

# Panel 2: C_K(alpha) przy stalym a_gam
ax = axes[1]
alpha_arr = np.array([4,5,6,7,8,9,10,12])
CK_arr = []
for al in alpha_arr:
    k1 = find_K1(al, 0.040, 5.501e-6)
    CK_arr.append(k1*(1+al)/0.040 if not np.isnan(k1) else np.nan)
CK_arr = np.array(CK_arr)
mask = np.isfinite(CK_arr)
ax.plot(alpha_arr[mask], CK_arr[mask], 'ro-', lw=2, ms=8, label='C_K(alpha, a=0.04)')
ax.axhline(C_K_mean, color='blue', lw=2, linestyle='--', label=f'C_K mean={C_K_mean:.4f}')
ax.set_xlabel('alpha'); ax.set_ylabel('C_K = K*1*(1+alpha)/a_gam')
ax.set_title('Stalosciosc C_K przy zmiennym alpha', fontsize=9)
ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

# Panel 3: C_K(a_gam) przy stalym alpha
ax = axes[2]
agam_arr = np.array([0.010, 0.020, 0.030, 0.040, 0.050, 0.060, 0.080])
CK_agam = []
for ag in agam_arr:
    k1 = find_K1(8.5616, ag, 5.501e-6, K_lo=0.001, K_hi=0.060)
    CK_agam.append(k1*(1+8.5616)/ag if not np.isnan(k1) else np.nan)
CK_agam = np.array(CK_agam)
mask2 = np.isfinite(CK_agam)
ax.semilogx(agam_arr[mask2], CK_agam[mask2], 'gs-', lw=2, ms=8, label='C_K(a_gam, alpha=8.56)')
ax.axhline(C_K_mean, color='blue', lw=2, linestyle='--', label=f'C_K mean={C_K_mean:.4f}')
ax.set_xlabel('a_gam'); ax.set_ylabel('C_K = K*1*(1+alpha)/a_gam')
ax.set_title('Stalosciosc C_K przy zmiennym a_gam', fontsize=9)
ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

plt.tight_layout()
out = __file__.replace('.py', '.png')
plt.savefig(out, dpi=120, bbox_inches='tight')
plt.close()
print(f"Zapisano: {out}")

print()
print("=" * 65)
print("PODSUMOWANIE P24:")
print("=" * 65)
print()
print(f"  Wynik: K*1 ≈ C_K * a_gam / (1+alpha)")
print(f"  C_K = {C_K_mean:.5f} ± {C_K_std:.5f}  [{C_K_std/C_K_mean*100:.2f}% rozrzut]")
print()
print(f"  Analitycznie:")
print(f"    K*1_LO = 2*a_gam/(1+alpha)  [c_K=2, LO]")
print(f"    Korekcja: C_K/2 = {C_K_mean/2:.5f}")
print(f"    C_K ≈ 2 / [a_gam * I_k(a_gam)]  (ze wzoru perturbacyjnego)")
print(f"    Dla a_gam -> 0: C_K -> 2 (czysta teoria)")
print()
print("GOTOWE: P24 zakonczone.")
