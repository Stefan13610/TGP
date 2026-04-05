"""
p22_koide.py  (v2 -- poprawiony wzor na A)
==========================================
CEL: Weryfikacja formuly Koide dla trojki K*1, K*2, K*3 z TGP.

FORMULA KOIDE (oryginal Koide 1982):
  sqrt(m_n) = A * (1 + sqrt(2)*cos(theta + 2*pi*(n-1)/3))
  => Q = (sum sqrt(m_n))^2 / sum(m_n) = 3/2  identycznie

POPRAWNY WZOR NA A:
  A = (sqrt(k1) + sqrt(k2) + sqrt(k3)) / 3
  (NIE: sqrt((k1+k2+k3)/3) -- to dawaloby A*sqrt(2), bledny theta)

KROKI:
  1. Oblicz Q dla K* z p21b (log-grid)
  2. Wyznacz A i theta z poprawnej parametryzacji
  3. Porownaj theta_TGP z theta_lepton (me, mmu, mtau)
  4. Wyznacz theta_Koide takie ze r21=206.768 (dokladna wartosc leptonu)
  5. Rownanie masy TGP: K_n^* = A^2*(1+sqrt(2)*cos(theta+2pi*(n-1)/3))^2
"""

import numpy as np
import warnings
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import brentq
warnings.filterwarnings('ignore')

# ============================================================
# STALE
# ============================================================
ALPHA  = 8.5616
LAM    = 5.501357e-06
A_GAM  = 0.040

# Referencyjne K* z p21b (log-grid, <0.003% od referencji)
K1_REF = 0.009820
K2_REF = 2.032724
K3_REF = 34.145108

# Masy leptonow [MeV]
ME   = 0.510999
MMU  = 105.6584
MTAU = 1776.86

print("P22: Weryfikacja Koide dla TGP K*_n  (v2 -- poprawiony wzor na A)")
print("=" * 65)
print()

# ============================================================
# DEFINICJE
# ============================================================
def koide_Q(k1, k2, k3):
    """Q = (sqrt(k1)+sqrt(k2)+sqrt(k3))^2 / (k1+k2+k3)."""
    s = np.sqrt(k1) + np.sqrt(k2) + np.sqrt(k3)
    return s**2 / (k1 + k2 + k3)

def koide_params(k1, k2, k3):
    """
    Wyznacz A i theta z parametryzacji Koide:
      sqrt(k_n) = A * (1 + sqrt(2)*cos(theta + 2*pi*(n-1)/3))
    Poprawny wzor: A = (sqrt(k1)+sqrt(k2)+sqrt(k3)) / 3
    """
    s1, s2, s3 = np.sqrt(k1), np.sqrt(k2), np.sqrt(k3)
    A = (s1 + s2 + s3) / 3.0
    # x_n = sqrt(k_n)/A = 1 + sqrt(2)*cos(theta + 2pi(n-1)/3)
    x1, x2, x3 = s1/A, s2/A, s3/A
    # cos(theta) z n=1
    cos_th = (x1 - 1.0) / np.sqrt(2)
    # sin(theta) z n=2: cos(theta+2pi/3) = -1/2*cos - sqrt(3)/2*sin
    # => sin = [-1/2*cos - (x2-1)/sqrt(2)] / (sqrt(3)/2)
    sin_th_2 = (-0.5*cos_th - (x2 - 1.0)/np.sqrt(2)) / (np.sqrt(3)/2)
    # sin(theta) z n=3: cos(theta+4pi/3) = -1/2*cos + sqrt(3)/2*sin
    # => sin = [(x3-1)/sqrt(2) + 1/2*cos] / (sqrt(3)/2)
    sin_th_3 = ((x3 - 1.0)/np.sqrt(2) + 0.5*cos_th) / (np.sqrt(3)/2)
    sin_th = (sin_th_2 + sin_th_3) / 2.0
    theta = np.arctan2(sin_th, cos_th)
    return A, theta

def koide_reconstruct(A, theta):
    """Rekonstruuj K* z A i theta."""
    k_vals = [A**2 * (1 + np.sqrt(2)*np.cos(theta + 2*np.pi*n/3))**2 for n in range(3)]
    return sorted(k_vals)

def r21_from_theta(theta):
    """r21 = k2/k1 = (f2/f1)^2  dla parametryzacji Koide."""
    f1 = 1 + np.sqrt(2)*np.cos(theta)
    f2 = 1 + np.sqrt(2)*np.cos(theta + 2*np.pi/3)
    if f1 == 0 or f2 == 0: return np.nan
    return (f2/f1)**2

def r31_from_theta(theta):
    """r31 = k3/k1 = (f3/f1)^2."""
    f1 = 1 + np.sqrt(2)*np.cos(theta)
    f3 = 1 + np.sqrt(2)*np.cos(theta + 4*np.pi/3)
    if f1 == 0 or f3 == 0: return np.nan
    return (f3/f1)**2

# ============================================================
# KROK 1: Weryfikacja Q = 3/2
# ============================================================
print("KROK 1: Weryfikacja formuly Koide (Q = 3/2)")
print("-" * 65)

Q_ref = koide_Q(K1_REF, K2_REF, K3_REF)
Q_lep = koide_Q(ME, MMU, MTAU)
r21_TGP = K2_REF / K1_REF
r31_TGP = K3_REF / K1_REF
r21_lep = MMU / ME
r31_lep = MTAU / ME

print(f"  TGP:     K1={K1_REF:.6f}, K2={K2_REF:.6f}, K3={K3_REF:.6f}")
print(f"           Q = {Q_ref:.6f}  (cel 1.5,  odchylenie: {abs(Q_ref-1.5)/1.5*100:.4f}%)")
print(f"           r21 = {r21_TGP:.3f},  r31 = {r31_TGP:.2f}")
print()
print(f"  Leptony: me={ME:.3f}, mmu={MMU:.4f}, mtau={MTAU:.2f} MeV")
print(f"           Q = {Q_lep:.6f}  (cel 1.5,  odchylenie: {abs(Q_lep-1.5)/1.5*100:.4f}%)")
print(f"           r21 = {r21_lep:.3f},  r31 = {r31_lep:.2f}")
print()

# Rodzina optymalna
family = [
    (5.9148, 0.025, 0.009820, 2.03303, 34.245),
    (6.8675, 0.030, 0.009820, 2.03289, 34.202),
    (7.7449, 0.035, 0.009820, 2.03279, 34.168),
    (8.5616, 0.040, 0.009820, 2.03272, 34.145),
]
print("  Rodzina optymalna (r21=207, r31=3477):")
print(f"  {'alpha':>8} {'a_gam':>7}  {'Q':>10}  {'dev(Q)%':>9}")
for alpha, agam, k1, k2, k3 in family:
    Q = koide_Q(k1, k2, k3)
    print(f"  {alpha:>8.4f} {agam:>7.4f}  {Q:>10.6f}  {abs(Q-1.5)/1.5*100:>9.4f}%")
print()

# ============================================================
# KROK 2: Parametryzacja Koide -- wyznaczenie A i theta
# ============================================================
print("KROK 2: Parametryzacja Koide -- A i theta")
print("-" * 65)
print("  Poprawna definicja: sqrt(k_n) = A*(1+sqrt(2)*cos(theta+2pi*(n-1)/3))")
print("  Poprawny wzor:      A = (sqrt(k1)+sqrt(k2)+sqrt(k3)) / 3")
print()

A_ref, theta_ref = koide_params(K1_REF, K2_REF, K3_REF)
A_lep, theta_lep = koide_params(ME, MMU, MTAU)

print(f"  TGP:     A = {A_ref:.6f},  theta = {theta_ref:.6f} rad = {np.degrees(theta_ref):.4f} deg")
print(f"  Leptony: A = {A_lep:.4f} MeV^0.5, theta = {theta_lep:.6f} rad = {np.degrees(theta_lep):.4f} deg")
print()

# Weryfikacja rekonstrukcji
K_rec = koide_reconstruct(A_ref, theta_ref)
print(f"  Rekonstrukcja TGP z (A, theta):")
for i, (K_orig, K_r, label) in enumerate(zip([K1_REF,K2_REF,K3_REF], K_rec, ['K1','K2','K3'])):
    err = abs(K_r-K_orig)/K_orig*100
    print(f"    K*{i+1} = {K_r:.6f}  (oryginal: {K_orig:.6f}, blad: {err:.4f}%)")
print()

# Rodzina optymalna -- theta
print(f"  Theta Koide dla rodziny optymalnej:")
print(f"  {'alpha':>8} {'a_gam':>7}  {'A':>10}  {'theta[rad]':>11}  {'theta[deg]':>11}")
for alpha, agam, k1, k2, k3 in family:
    A, th = koide_params(k1, k2, k3)
    print(f"  {alpha:>8.4f} {agam:>7.4f}  {A:>10.6f}  {th:>11.6f}  {np.degrees(th):>11.4f}")
print()

# Roznica theta
dtheta_deg = np.degrees(theta_ref) - np.degrees(theta_lep)
dtheta_rad = theta_ref - theta_lep
print(f"  theta_TGP    = {np.degrees(theta_ref):.6f} deg")
print(f"  theta_lepton = {np.degrees(theta_lep):.6f} deg")
print(f"  Delta_theta  = {dtheta_deg:+.6f} deg = {dtheta_rad*1000:.4f} mrad")
print()

# ============================================================
# KROK 3: Theta dla r21=206.768 (dokladna wartosc leptonu)
# ============================================================
print("KROK 3: Theta dla r21 = 206.768 (dokladna wartosc leptonu)")
print("-" * 65)

r21_target = r21_lep  # 206.768

theta_scan = np.linspace(1.5, 2.5, 4000)
r21_scan = np.array([r21_from_theta(t) for t in theta_scan])

theta_exact = None
for i in range(len(theta_scan)-1):
    r1, r2 = r21_scan[i], r21_scan[i+1]
    if np.isfinite(r1) and np.isfinite(r2) and np.isfinite(r1*r2):
        if (r1 - r21_target) * (r2 - r21_target) < 0:
            try:
                theta_exact = brentq(lambda t: r21_from_theta(t) - r21_target,
                                     theta_scan[i], theta_scan[i+1], xtol=1e-12)
            except Exception:
                pass
            break

if theta_exact is not None:
    r21_chk = r21_from_theta(theta_exact)
    r31_chk = r31_from_theta(theta_exact)
    # A z K1_fixed
    f1_ex = 1 + np.sqrt(2)*np.cos(theta_exact)
    A_exact = np.sqrt(K1_REF) / abs(f1_ex)
    K_exact = koide_reconstruct(A_exact, theta_exact)
    print(f"  theta_exact = {theta_exact:.8f} rad = {np.degrees(theta_exact):.6f} deg")
    print(f"  Sprawdzenie: r21 = {r21_chk:.6f} (cel: {r21_target:.6f})")
    print(f"               r31 = {r31_chk:.4f}  (leptony: {r31_lep:.4f})")
    print(f"  K* dla theta_exact (z K1 = {K1_REF}):")
    for i, (K_e, label) in enumerate(zip(K_exact, ['K1','K2','K3'])):
        print(f"    K*{i+1} = {K_e:.6f}  (TGP: {[K1_REF,K2_REF,K3_REF][i]:.6f})")
    print()
    print(f"  Roznica: theta_TGP - theta_exact = {np.degrees(theta_ref-theta_exact):+.6f} deg")
    print(f"                                   = {(theta_ref-theta_exact)*1000:+.4f} mrad")
else:
    print("  BRAK przeciecia!")
print()

# ============================================================
# KROK 4: Tabela zbiorcza
# ============================================================
print("KROK 4: Tabela zbiorcza -- TGP vs Leptony")
print("-" * 65)
print()
print(f"  {'Wielkosc':>28}  {'TGP':>14}  {'Leptony':>14}  {'Roznica':>14}")
print("  " + "-" * 72)
print(f"  {'Q (Koide = 3/2)':>28}  {Q_ref:>14.6f}  {Q_lep:>14.6f}  {(Q_ref-Q_lep)*100/1.5:>+12.4f}%")
print(f"  {'r21 = K*2/K*1':>28}  {r21_TGP:>14.3f}  {r21_lep:>14.3f}  {(r21_TGP-r21_lep)/r21_lep*100:>+12.4f}%")
print(f"  {'r31 = K*3/K*1':>28}  {r31_TGP:>14.2f}  {r31_lep:>14.2f}  {(r31_TGP-r31_lep)/r31_lep*100:>+12.4f}%")
print(f"  {'theta [deg]':>28}  {np.degrees(theta_ref):>14.4f}  {np.degrees(theta_lep):>14.4f}  {np.degrees(theta_ref-theta_lep):>+12.4f} deg")
print(f"  {'theta [mrad]':>28}  {theta_ref*1000:>14.4f}  {theta_lep*1000:>14.4f}  {(theta_ref-theta_lep)*1000:>+12.4f} mrad")
print()

# ============================================================
# KROK 5: Kandydat na rowranie masy TGP
# ============================================================
print("KROK 5: Kandydat na 'rownanie masy' TGP")
print("=" * 65)
print()
print("  ROWNANIE MASY TGP (Koide-Radiative Form):")
print()
print("    K_n^* = A^2 * (1 + sqrt(2) * cos(theta + 2*pi*(n-1)/3))^2,   n=1,2,3")
print()
print(f"  Parametry TGP (alpha={ALPHA}, a_gam={A_GAM}, lam={LAM:.4e}):")
print(f"    A      = {A_ref:.6f}")
print(f"    theta  = {theta_ref:.6f} rad = {np.degrees(theta_ref):.4f} deg")
print()
print(f"  Parametry leptonowe (Koide dokladny, Q={Q_lep:.6f}):")
print(f"    A      = {A_lep:.6f} MeV^0.5")
print(f"    theta  = {theta_lep:.6f} rad = {np.degrees(theta_lep):.4f} deg")
print()
print(f"  ODCHYLENIE theta TGP od leptonu:")
print(f"    Delta_theta = {(theta_ref-theta_lep)*1000:.3f} mrad = {np.degrees(theta_ref-theta_lep)*60:.2f} arcmin")
print()
print(f"  WNIOSKI:")
print(f"    1. Q_TGP = {Q_ref:.6f} ≈ 3/2  [{abs(Q_ref-1.5)/1.5*100:.4f}% odchylenia]")
print(f"    2. theta_TGP = {np.degrees(theta_ref):.4f} deg  vs  theta_lepton = {np.degrees(theta_lep):.4f} deg")
print(f"       Roznica: {np.degrees(theta_ref-theta_lep):.4f} deg = {(theta_ref-theta_lep)*1000:.3f} mrad")
print(f"    3. r21: TGP={r21_TGP:.3f} vs lepton={r21_lep:.3f}  (odchylenie {abs(r21_TGP-r21_lep)/r21_lep*100:.3f}%)")
print(f"    4. Rownanie masy: K_n = A^2*(1+sqrt(2)*cos(theta+2pi*n/3))^2")
print(f"    5. TGP wyznacza theta ≈ theta_lepton z dokladnoscia ~{abs(np.degrees(theta_ref-theta_lep))*60:.1f}'")
print()

# ============================================================
# WYKRES
# ============================================================
print("Tworze wykres...")
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle(
    f'P22: Formula Koide dla TGP   (alpha={ALPHA}, a_gam={A_GAM})\n'
    f'K_n^* = A^2*(1+sqrt(2)*cos(theta+2pi(n-1)/3))^2'
    f'   theta_TGP={np.degrees(theta_ref):.3f} deg,  theta_lep={np.degrees(theta_lep):.3f} deg',
    fontsize=10, fontweight='bold')

# --- Panel 1: Krzywa Koide r31(r21) ---
ax = axes[0, 0]
th_arr = np.linspace(0.5, 2.4, 5000)
r21_arr = np.array([r21_from_theta(t) for t in th_arr])
r31_arr = np.array([r31_from_theta(t) for t in th_arr])
mask = (np.isfinite(r21_arr) & np.isfinite(r31_arr) &
        (r21_arr > 0) & (r31_arr > 0) & (r21_arr < 1000) & (r31_arr < 20000))
ax.plot(r21_arr[mask], r31_arr[mask], 'b-', lw=1.5, label='Krzywa Koide Q=3/2', alpha=0.6)
ax.scatter([r21_TGP], [r31_TGP], color='red', s=150, zorder=6,
           label=f'TGP: ({r21_TGP:.1f}, {r31_TGP:.0f})')
ax.scatter([r21_lep], [r31_lep], color='green', s=150, marker='*', zorder=6,
           label=f'Lepton: ({r21_lep:.1f}, {r31_lep:.0f})')
ax.set_xlabel('r21 = K2/K1'); ax.set_ylabel('r31 = K3/K1')
ax.set_title('Krzywa Koide: TGP vs Leptony', fontsize=9)
ax.legend(fontsize=8); ax.grid(True, alpha=0.3)
ax.set_xlim(100, 300); ax.set_ylim(1000, 6000)

# --- Panel 2: theta dla rodziny TGP ---
ax = axes[0, 1]
alphas = [f[0] for f in family]
thetas = [np.degrees(koide_params(f[2],f[3],f[4])[1]) for f in family]
ax.plot(alphas, thetas, 'ro-', lw=2, ms=8, label='theta_TGP(alpha)')
ax.axhline(np.degrees(theta_lep), color='green', lw=2, linestyle='--', label=f'theta_lepton={np.degrees(theta_lep):.4f} deg')
ax.axhline(np.degrees(theta_ref), color='red', lw=1.5, linestyle=':', alpha=0.5)
ax.set_xlabel('alpha'); ax.set_ylabel('theta [deg]')
ax.set_title('Kat Koide theta vs alpha (rodzina optymalna)', fontsize=9)
ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

# --- Panel 3: r21(theta) z zaznaczeniem TGP i leptonu ---
ax = axes[1, 0]
th2 = np.linspace(2.2, 2.4, 2000)
r21_2 = np.array([r21_from_theta(t) for t in th2])
mask2 = np.isfinite(r21_2) & (r21_2 > 0) & (r21_2 < 500)
ax.plot(np.degrees(th2[mask2]), r21_2[mask2], 'b-', lw=2)
ax.axhline(r21_TGP, color='red', lw=1.5, linestyle='--', label=f'r21_TGP={r21_TGP:.3f}')
ax.axhline(r21_lep, color='green', lw=1.5, linestyle=':', label=f'r21_lep={r21_lep:.3f}')
ax.axvline(np.degrees(theta_ref), color='red', lw=1.5, linestyle='--', alpha=0.5)
ax.axvline(np.degrees(theta_lep), color='green', lw=1.5, linestyle=':', alpha=0.5)
ax.set_xlabel('theta [deg]'); ax.set_ylabel('r21 = K2/K1')
ax.set_title('r21 vs theta (Koide)', fontsize=9)
ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

# --- Panel 4: Diagram Koide (wektory) ---
ax = axes[1, 1]
th_circle = np.linspace(0, 2*np.pi, 200)
ax.plot(np.cos(th_circle), np.sin(th_circle), 'k-', lw=0.8, alpha=0.3)
labels_tgp = ['K*1', 'K*2', 'K*3']
labels_lep = ['e', 'mu', 'tau']
for n, lab in enumerate(labels_tgp):
    ang = theta_ref + 2*np.pi*n/3
    r_n = (1 + np.sqrt(2)*np.cos(ang))
    ax.annotate('', xy=(r_n*np.cos(ang), r_n*np.sin(ang)), xytext=(0, 0),
                arrowprops=dict(arrowstyle='->', color='red', lw=2))
    ax.text(r_n*np.cos(ang)*1.08, r_n*np.sin(ang)*1.08, f'TGP {lab}', fontsize=7, color='red')
for n, lab in enumerate(labels_lep):
    ang = theta_lep + 2*np.pi*n/3
    r_n = (1 + np.sqrt(2)*np.cos(ang))
    ax.annotate('', xy=(r_n*np.cos(ang), r_n*np.sin(ang)), xytext=(0, 0),
                arrowprops=dict(arrowstyle='->', color='darkgreen', lw=2, linestyle='--'))
    ax.text(r_n*np.cos(ang)*1.08, r_n*np.sin(ang)*1.08, f'LEP {lab}', fontsize=7, color='darkgreen')
ax.set_aspect('equal'); ax.grid(True, alpha=0.25)
ax.set_xlim(-2.5, 2.5); ax.set_ylim(-2.5, 2.5)
ax.set_title(f'Diagram Koide: czerwone=TGP, zielone=Leptony\n'
             f'Delta_theta={np.degrees(theta_ref-theta_lep)*60:.2f} arcmin', fontsize=8)

plt.tight_layout()
out = __file__.replace('.py', '.png')
plt.savefig(out, dpi=120, bbox_inches='tight')
plt.close()
print(f"Zapisano: {out}")

# ============================================================
# PODSUMOWANIE KONCOWE
# ============================================================
print()
print("=" * 65)
print("PODSUMOWANIE P22 (v2):")
print("=" * 65)
print()
print(f"  1. TGP K* spelniaja Koide: Q = {Q_ref:.6f}  [{abs(Q_ref-1.5)/1.5*100:.4f}% od 3/2]")
print(f"     Leptony:                Q = {Q_lep:.6f}  [{abs(Q_lep-1.5)/1.5*100:.4f}% od 3/2]")
print()
print(f"  2. Kat Koide theta:")
print(f"     theta_TGP    = {np.degrees(theta_ref):.6f} deg = {theta_ref:.8f} rad")
print(f"     theta_lepton = {np.degrees(theta_lep):.6f} deg = {theta_lep:.8f} rad")
print(f"     Roznica:       {np.degrees(theta_ref-theta_lep):.6f} deg = {(theta_ref-theta_lep)*1000:.3f} mrad")
print()
print(f"  3. ROWNANIE MASY TGP:")
print(f"     K_n^* = A^2 * (1 + sqrt(2)*cos(theta + 2*pi*(n-1)/3))^2")
print(f"     TGP:  A = {A_ref:.6f},  theta = {np.degrees(theta_ref):.4f} deg")
print(f"     Lept: A = {A_lep:.6f} MeV^0.5, theta = {np.degrees(theta_lep):.4f} deg")
print()
print(f"  4. r21: TGP = {r21_TGP:.3f},  lepton = {r21_lep:.3f},  odchylenie = {abs(r21_TGP-r21_lep)/r21_lep*100:.4f}%")
print(f"     r31: TGP = {r31_TGP:.2f}, lepton = {r31_lep:.2f}, odchylenie = {abs(r31_TGP-r31_lep)/r31_lep*100:.4f}%")
print()
print(f"  5. KLUCZOWY WYNIK:")
print(f"     theta_TGP = theta_lepton z dokladnoscia {abs(np.degrees(theta_ref-theta_lep))*60:.2f} arcmin")
print(f"     Rownanie masy K_n^* = A^2*(1+sqrt(2)*cos(theta+2pi*n/3))^2")
print(f"     jest wspolne dla TGP i leptow ze ZBIEZNYM theta!")
print()
print("GOTOWE: P22 (v2) zakonczone.")
