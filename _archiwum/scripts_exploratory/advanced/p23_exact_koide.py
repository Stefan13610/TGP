"""
p23_exact_koide.py
==================
CEL: Znalezc dokladny punkt TGP na krzywej rodziny optymalnej
     gdzie Q = 3/2 DOKLADNIE (formula Koide spelniona perfekcyjnie).

MOTYWACJA (z P22):
  Rodzina optymalna (r21~207, r31~3477):
    alpha=6.8675, a_gam=0.030 -> Q=1.499892  (ponizej 3/2)
    alpha=7.7449, a_gam=0.035 -> Q=1.500107  (powyzej 3/2)
  => Istnieje punkt (alpha*, a_gam*) gdzie Q=1.5 DOKLADNIE.
  Szacunek: alpha*~7.3, a_gam*~0.0325

PYTANIE: W tym punkcie:
  - Czy r21 = 206.768 (dokladna wartosc leptonu)?
  - Czy theta = theta_lepton = 132.7328 deg?

METODA:
  1. Interpoluj lambda*(alpha) i a_gam(alpha) z rodziny
  2. Dla alpha w [6.8, 7.8]: oblicz K*(alpha) -> Q(alpha)
  3. Wyznacz alpha* gdzie Q=1.5 (brentq)
  4. Oblicz r21, r31, theta dla alpha*
  5. Porownaj z leptonami

DEFINICJA MODELU (identyczna z p21b):
  g(K, alpha, a_gam, lam) = E_Yukawa(K) / K - 4*pi = 0
  E_Yukawa z siatka logarytmiczna (jak v2_bisekcja_alpha)
"""

import numpy as np
from scipy.optimize import brentq, minimize_scalar
from scipy.interpolate import interp1d
import warnings
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')

# ============================================================
# STALE
# ============================================================
R_MAX    = 60.0
GAMMA    = 1.0
N_ENERGY = 1500   # jak p21b

# Znane punkty rodziny optymalnej (z p21b)
FAM_ALPHA = np.array([5.9148, 6.8675, 7.7449, 8.5616])
FAM_AGAM  = np.array([0.025,  0.030,  0.035,  0.040 ])
FAM_LAM   = np.array([2.883e-6, 3.743e-6, 4.618e-6, 5.501e-6])
FAM_K1    = np.array([0.009820]*4)
FAM_K2    = np.array([2.03303,  2.03289, 2.03279, 2.03272])
FAM_K3    = np.array([34.245,   34.202,  34.168,  34.145 ])

# Interpolanty 1D (liniowe)
f_agam_of_alpha = interp1d(FAM_ALPHA, FAM_AGAM,  kind='linear', fill_value='extrapolate')
f_lam_of_alpha  = interp1d(FAM_ALPHA, FAM_LAM,   kind='linear', fill_value='extrapolate')

# Funkcje Koide
def koide_Q(k1, k2, k3):
    s = np.sqrt(k1) + np.sqrt(k2) + np.sqrt(k3)
    return s**2 / (k1 + k2 + k3)

def koide_params(k1, k2, k3):
    s1, s2, s3 = np.sqrt(k1), np.sqrt(k2), np.sqrt(k3)
    A = (s1 + s2 + s3) / 3.0
    cos_th = (s1/A - 1.0) / np.sqrt(2)
    sin_th2 = (-0.5*cos_th - (s2/A - 1.0)/np.sqrt(2)) / (np.sqrt(3)/2)
    sin_th3 = ((s3/A - 1.0)/np.sqrt(2) + 0.5*cos_th) / (np.sqrt(3)/2)
    sin_th = (sin_th2 + sin_th3) / 2.0
    theta = np.arctan2(sin_th, cos_th)
    return A, theta

print("P23: Poszukiwanie punktu TGP gdzie Q = 3/2 DOKLADNIE")
print("=" * 65)
print()
print("Rodzina optymalna (interpolacja parametrow):")
print(f"  a_gam(alpha) ~ {f_agam_of_alpha(7.0):.4f} dla alpha=7.0")
print(f"  lam(alpha)   ~ {f_lam_of_alpha(7.0):.4e} dla alpha=7.0")
print()

# ============================================================
# FUNKCJE ENERGII (siatka log, jak p21b)
# ============================================================
def V_mod(phi, lam):
    return GAMMA/3*phi**3 - GAMMA/4*phi**4 + lam/6*(phi-1)**6

def energy_log(K, alpha, a_gam, lam, N=N_ENERGY):
    t    = np.linspace(0, 1, N)
    r    = a_gam * (R_MAX/a_gam)**t
    phi  = np.maximum(1.0 + K*np.exp(-r)/r, 1e-10)
    dphi = K*np.exp(-r)*(-r - 1.0)/r**2
    V1   = V_mod(1.0, lam)
    Ek   = 4*np.pi * np.trapezoid(0.5*dphi**2 * (1+alpha/phi) * r**2, r)
    Ep   = 4*np.pi * np.trapezoid((V_mod(phi, lam) - V1) * r**2, r)
    return Ek + Ep

def g_func(K, alpha, a_gam, lam, N=N_ENERGY):
    return energy_log(K, alpha, a_gam, lam, N) / (4*np.pi*K) - 1.0

def find_K_star(K_lo, K_hi, alpha, a_gam, lam, label=''):
    """Znajdz zero g(K) w [K_lo, K_hi] metoda brentq."""
    try:
        g_lo = g_func(K_lo, alpha, a_gam, lam)
        g_hi = g_func(K_hi, alpha, a_gam, lam)
        if not (np.isfinite(g_lo) and np.isfinite(g_hi)):
            return np.nan
        if g_lo * g_hi > 0:
            return np.nan
        Kstar = brentq(lambda K: g_func(K, alpha, a_gam, lam),
                       K_lo, K_hi, xtol=1e-7, rtol=1e-7, maxiter=50)
        return Kstar
    except Exception:
        return np.nan

def compute_Kstars(alpha, a_gam, lam, verbose=False):
    """Oblicz K*1, K*2, K*3 dla zadanych parametrow."""
    # K*1 ~ 0.0098 (waska okolica)
    K1 = find_K_star(0.006, 0.014, alpha, a_gam, lam)
    # K*2 ~ 2.03 (z wiekszym zakresem)
    K2 = find_K_star(1.0, 3.0, alpha, a_gam, lam)
    # K*3 ~ 34 (szeroki zakres)
    K3 = find_K_star(20.0, 60.0, alpha, a_gam, lam)
    if verbose:
        print(f"    K1={K1:.6f}, K2={K2:.6f}, K3={K3:.5f}")
    return K1, K2, K3

# ============================================================
# WERYFIKACJA PUNKTOW RODZINY
# ============================================================
print("KROK 0: Weryfikacja znanych punktow rodziny")
print("-" * 65)
for i, (alpha, agam, lam) in enumerate(zip(FAM_ALPHA, FAM_AGAM, FAM_LAM)):
    K1, K2, K3 = compute_Kstars(alpha, agam, lam, verbose=False)
    if np.any(np.isnan([K1,K2,K3])):
        print(f"  alpha={alpha:.4f}: NaN!")
        continue
    Q = koide_Q(K1, K2, K3)
    _, theta = koide_params(K1, K2, K3)
    print(f"  alpha={alpha:.4f}, a_gam={agam:.4f}: K1={K1:.6f}, K2={K2:.5f}, K3={K3:.4f}")
    print(f"    r21={K2/K1:.4f}, r31={K3/K1:.2f}, Q={Q:.6f}, theta={np.degrees(theta):.4f} deg")
print()

# ============================================================
# SKAN ALPHA: szukaj Q(alpha) = 1.5
# ============================================================
print("KROK 1: Skan alpha w [6.7, 7.9] -> Q(alpha), szukamy Q=1.5")
print("-" * 65)

alpha_scan = np.linspace(6.7, 7.9, 13)  # gesty skan
Q_scan = []
r21_scan = []
theta_scan_deg = []
K_data = []

for alpha in alpha_scan:
    agam = float(f_agam_of_alpha(alpha))
    lam  = float(f_lam_of_alpha(alpha))
    K1, K2, K3 = compute_Kstars(alpha, agam, lam)
    if np.any(np.isnan([K1,K2,K3])):
        Q_scan.append(np.nan); r21_scan.append(np.nan)
        theta_scan_deg.append(np.nan); K_data.append((np.nan,)*3)
        continue
    Q = koide_Q(K1, K2, K3)
    r21 = K2/K1
    _, theta = koide_params(K1, K2, K3)
    Q_scan.append(Q)
    r21_scan.append(r21)
    theta_scan_deg.append(np.degrees(theta))
    K_data.append((K1, K2, K3))
    tag = '  <Q=1.5>' if abs(Q-1.5) < 0.0002 else ''
    print(f"  alpha={alpha:.4f}: Q={Q:.6f}, r21={r21:.4f}, theta={np.degrees(theta):.4f} deg{tag}")

Q_scan = np.array(Q_scan)
r21_scan = np.array(r21_scan)
theta_scan_deg = np.array(theta_scan_deg)
print()

# ============================================================
# ZNAJDZ ALPHA* GDZIE Q = 1.5
# ============================================================
print("KROK 2: Wyznaczenie alpha* (brentq) gdzie Q=1.5")
print("-" * 65)

# Znajdz przedzial z zeropsem
alpha_lo_Q = None
alpha_hi_Q = None
for i in range(len(alpha_scan)-1):
    q1, q2 = Q_scan[i], Q_scan[i+1]
    if np.isfinite(q1) and np.isfinite(q2):
        if (q1-1.5)*(q2-1.5) < 0:
            alpha_lo_Q = alpha_scan[i]
            alpha_hi_Q = alpha_scan[i+1]
            break

alpha_star = None
K1_star = K2_star = K3_star = np.nan

if alpha_lo_Q is not None:
    def Q_minus_15(alpha):
        agam = float(f_agam_of_alpha(alpha))
        lam  = float(f_lam_of_alpha(alpha))
        K1, K2, K3 = compute_Kstars(alpha, agam, lam)
        if np.any(np.isnan([K1,K2,K3])): return np.nan
        return koide_Q(K1, K2, K3) - 1.5

    try:
        alpha_star = brentq(Q_minus_15, alpha_lo_Q, alpha_hi_Q, xtol=1e-4, maxiter=30)
        agam_star = float(f_agam_of_alpha(alpha_star))
        lam_star  = float(f_lam_of_alpha(alpha_star))
        K1_star, K2_star, K3_star = compute_Kstars(alpha_star, agam_star, lam_star, verbose=True)
        Q_check = koide_Q(K1_star, K2_star, K3_star)
        _, theta_star = koide_params(K1_star, K2_star, K3_star)
        r21_star = K2_star / K1_star
        r31_star = K3_star / K1_star
        print(f"  alpha* = {alpha_star:.6f}")
        print(f"  a_gam* = {agam_star:.6f}")
        print(f"  lam*   = {lam_star:.4e}")
        print(f"  K*1 = {K1_star:.8f}")
        print(f"  K*2 = {K2_star:.8f}")
        print(f"  K*3 = {K3_star:.6f}")
        print(f"  Q   = {Q_check:.8f}  (cel: 1.500000)")
        print(f"  r21 = {r21_star:.6f}  (lepton: 206.768311)")
        print(f"  r31 = {r31_star:.4f}   (lepton: 3477.228)")
        print(f"  theta = {np.degrees(theta_star):.6f} deg  (lepton: 132.7328 deg)")
    except Exception as e:
        print(f"  BLAD brentq: {e}")
else:
    print("  Brak przejscia Q=1.5 w zadanym zakresie alpha!")
print()

# ============================================================
# KROK 3: Tabela Q=1.5 vs Leptony
# ============================================================
ME, MMU, MTAU = 0.510999, 105.6584, 1776.86
Q_lep = koide_Q(ME, MMU, MTAU)
r21_lep = MMU/ME
r31_lep = MTAU/ME
_, theta_lep = koide_params(ME, MMU, MTAU)

if alpha_star is not None and not np.isnan(K1_star):
    Q_check2 = koide_Q(K1_star, K2_star, K3_star)
    _, theta_star = koide_params(K1_star, K2_star, K3_star)
    r21_star = K2_star/K1_star
    r31_star = K3_star/K1_star

    print("KROK 3: Tabela -- TGP(Q=1.5) vs Leptony")
    print("-" * 65)
    print(f"  {'Wielkosc':>25}  {'TGP(Q=1.5)':>14}  {'Leptony':>14}  {'Diff%':>10}")
    print("  " + "-"*67)
    print(f"  {'Q':>25}  {Q_check2:>14.6f}  {Q_lep:>14.6f}  {(Q_check2-Q_lep)/1.5*100:>+10.4f}%")
    print(f"  {'r21':>25}  {r21_star:>14.6f}  {r21_lep:>14.6f}  {(r21_star-r21_lep)/r21_lep*100:>+10.4f}%")
    print(f"  {'r31':>25}  {r31_star:>14.4f}  {r31_lep:>14.4f}  {(r31_star-r31_lep)/r31_lep*100:>+10.4f}%")
    print(f"  {'theta [deg]':>25}  {np.degrees(theta_star):>14.6f}  {np.degrees(theta_lep):>14.6f}  {np.degrees(theta_star-theta_lep):>+10.4f} deg")
    print()

# ============================================================
# WYKRES
# ============================================================
print("Tworze wykres...")

fig, axes = plt.subplots(1, 3, figsize=(15, 5))
fig.suptitle('P23: Szukanie punktu TGP z Q=3/2 (Koide dokladny)\n'
             'Interpolacja rodziny optymalnej (alpha, a_gam, lam)',
             fontsize=11, fontweight='bold')

mask = np.isfinite(Q_scan)

ax = axes[0]
ax.plot(alpha_scan[mask], Q_scan[mask], 'bo-', lw=2, ms=7, label='Q(alpha) TGP')
ax.axhline(1.5, color='red', lw=2, linestyle='--', label='Q=3/2 (Koide)')
ax.axhline(Q_lep, color='green', lw=1.5, linestyle=':', label=f'Q_lep={Q_lep:.6f}')
if alpha_star is not None:
    ax.axvline(alpha_star, color='orange', lw=2, linestyle='-.',
               label=f'alpha*={alpha_star:.3f}')
ax.set_xlabel('alpha'); ax.set_ylabel('Q')
ax.set_title('Q(alpha) wzdluz rodziny optymalnej', fontsize=9)
ax.legend(fontsize=7); ax.grid(True, alpha=0.3)

ax = axes[1]
ax.plot(alpha_scan[mask], r21_scan[mask], 'ro-', lw=2, ms=7, label='r21(alpha) TGP')
ax.axhline(r21_lep, color='green', lw=2, linestyle='--', label=f'r21_lep={r21_lep:.3f}')
ax.axhline(207.0, color='orange', lw=1.5, linestyle=':', label='r21=207')
if alpha_star is not None:
    ax.axvline(alpha_star, color='orange', lw=2, linestyle='-.', label=f'alpha*={alpha_star:.3f}')
    ax.scatter([alpha_star], [r21_star], color='orange', s=120, zorder=5)
ax.set_xlabel('alpha'); ax.set_ylabel('r21 = K2/K1')
ax.set_title('r21(alpha) wzdluz rodziny optymalnej', fontsize=9)
ax.legend(fontsize=7); ax.grid(True, alpha=0.3)

ax = axes[2]
ax.plot(alpha_scan[mask], theta_scan_deg[mask], 'gs-', lw=2, ms=7, label='theta_TGP(alpha)')
ax.axhline(np.degrees(theta_lep), color='red', lw=2, linestyle='--',
           label=f'theta_lep={np.degrees(theta_lep):.4f} deg')
if alpha_star is not None:
    ax.axvline(alpha_star, color='orange', lw=2, linestyle='-.',
               label=f'alpha*={alpha_star:.3f}')
    ax.scatter([alpha_star], [np.degrees(theta_star)], color='orange', s=120, zorder=5)
ax.set_xlabel('alpha'); ax.set_ylabel('theta [deg]')
ax.set_title('theta_Koide(alpha) wzdluz rodziny', fontsize=9)
ax.legend(fontsize=7); ax.grid(True, alpha=0.3)

plt.tight_layout()
out = __file__.replace('.py', '.png')
plt.savefig(out, dpi=120, bbox_inches='tight')
plt.close()
print(f"Zapisano: {out}")

# ============================================================
# PODSUMOWANIE
# ============================================================
print()
print("=" * 65)
print("PODSUMOWANIE P23:")
print("=" * 65)
if alpha_star is not None and not np.isnan(K1_star):
    Q_final = koide_Q(K1_star, K2_star, K3_star)
    _, theta_final = koide_params(K1_star, K2_star, K3_star)
    print(f"  Punkt TGP z Q=3/2:")
    print(f"    alpha* = {alpha_star:.4f},  a_gam* = {agam_star:.4f},  lam* = {lam_star:.4e}")
    print(f"    Q      = {Q_final:.8f}  (odchylenie od 3/2: {abs(Q_final-1.5)/1.5*100:.6f}%)")
    print(f"    r21    = {K2_star/K1_star:.6f}  (lepton: 206.768,  roznica: {(K2_star/K1_star-r21_lep)/r21_lep*100:+.4f}%)")
    print(f"    r31    = {K3_star/K1_star:.4f}   (lepton: 3477.23, roznica: {(K3_star/K1_star-r31_lep)/r31_lep*100:+.4f}%)")
    print(f"    theta  = {np.degrees(theta_final):.6f} deg  (lepton: {np.degrees(theta_lep):.4f} deg)")
    print()
    print(f"  Wniosek:")
    if abs(K2_star/K1_star - r21_lep)/r21_lep < 0.002:
        print(f"    ** r21(TGP, Q=1.5) ≈ r21_lepton z dokladnoscia <0.2% **")
    else:
        print(f"    Roznica r21 = {abs(K2_star/K1_star-r21_lep)/r21_lep*100:.4f}% -- TGP i lepton roznia sie")
    if abs(np.degrees(theta_final-theta_lep)) < 0.01:
        print(f"    ** theta(TGP, Q=1.5) ≈ theta_lepton z dokladnoscia {abs(np.degrees(theta_final-theta_lep))*60:.2f} arcmin **")
else:
    print("  Brak punktu Q=1.5 w zbadanym zakresie.")
print()
print("GOTOWE: P23 zakonczone.")
