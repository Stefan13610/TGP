"""
p25_theta_mass_eq.py
====================
CEL: Analityczne wyprowadzenie kata Koide theta z parametrow TGP.

KLUCZOWA OBSERWACJA (z P22-P24):
  W parametryzacji Koide: sqrt(K_n) = A*(1+sqrt(2)*cos(theta+2pi(n-1)/3))
  Gdzie: A = (sqrt(K1)+sqrt(K2)+sqrt(K3))/3

  Mamy DOKLADNIE:
    f1 = sqrt(K1)/A = 3/(1 + sqrt(r21) + sqrt(r31))
    cos(theta) = (f1 - 1)/sqrt(2)

  WIEC: theta jest WYZNACZONY PRZEZ r21 I r31:

    cos(theta) = [3/(1 + sqrt(r21) + sqrt(r31)) - 1] / sqrt(2)

PYTANIA:
  1. Jaka jest wrazliwosc theta na r21 vs r31?
  2. Ile theta odbiega od leptonu ze wzgledu na Dr21=0.23 vs Dr31=0.13?
  3. Co "wyznacza" r21 i r31 w TGP?

PODSUMOWANIE ROWNANIA MASY TGP:
  K_n* = K1* * r_n1
  K1* = C_K * a_gam / (1+alpha),   C_K = 2.351
  r21 ≈ 207.000, r31 ≈ 3477.1  <- wyznaczone przez V_mod (potencjal)
  theta = arccos([3/(1+sqrt(r21)+sqrt(r31))-1]/sqrt(2))
        ≈ 132.74 deg (TGP) vs 132.73 deg (leptony)
"""

import numpy as np
import warnings
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')

# ============================================================
# STALE REFERENCYJNE
# ============================================================
# TGP (alpha=8.5616, a_gam=0.040, lam=5.501e-6)
K1_TGP = 0.009820
K2_TGP = 2.032724
K3_TGP = 34.145108
r21_TGP = K2_TGP / K1_TGP
r31_TGP = K3_TGP / K1_TGP

# Leptony
ME, MMU, MTAU = 0.510999, 105.6584, 1776.86
r21_lep = MMU / ME
r31_lep = MTAU / ME

print("P25: Rownanie masy TGP -- analityczna formula theta(r21, r31)")
print("=" * 65)
print()

# ============================================================
# DEFINICJE
# ============================================================
def theta_from_ratios(r21, r31):
    """
    Dokladna formula:
      f1 = 3 / (1 + sqrt(r21) + sqrt(r31))
      cos(theta) = (f1 - 1) / sqrt(2)
    """
    f1 = 3.0 / (1.0 + np.sqrt(r21) + np.sqrt(r31))
    cos_th = (f1 - 1.0) / np.sqrt(2)
    return np.arccos(np.clip(cos_th, -1, 1))

def koide_Q_from_ratios(r21, r31):
    """Q = (1 + sqrt(r21) + sqrt(r31))^2 / (1 + r21 + r31)."""
    return (1 + np.sqrt(r21) + np.sqrt(r31))**2 / (1 + r21 + r31)

def koide_params_exact(k1, k2, k3):
    """Koide params z dokladnego wzoru (A = sum_sqrt/3)."""
    s1, s2, s3 = np.sqrt(k1), np.sqrt(k2), np.sqrt(k3)
    A = (s1 + s2 + s3) / 3.0
    f1 = s1 / A
    cos_th = (f1 - 1.0) / np.sqrt(2)
    return A, np.arccos(np.clip(cos_th, -1, 1))

# ============================================================
# KROK 1: Weryfikacja formuly theta(r21, r31)
# ============================================================
print("KROK 1: Weryfikacja formuly theta = f(r21, r31)")
print("-" * 65)
print("  Formula: cos(theta) = [3/(1+sqrt(r21)+sqrt(r31)) - 1] / sqrt(2)")
print()

theta_TGP_form = theta_from_ratios(r21_TGP, r31_TGP)
theta_TGP_exact = koide_params_exact(K1_TGP, K2_TGP, K3_TGP)[1]
theta_lep_form  = theta_from_ratios(r21_lep, r31_lep)
theta_lep_exact = koide_params_exact(ME, MMU, MTAU)[1]

print(f"  TGP:    theta_formula = {np.degrees(theta_TGP_form):.8f} deg")
print(f"          theta_exact   = {np.degrees(theta_TGP_exact):.8f} deg")
print(f"          Roznica       = {abs(np.degrees(theta_TGP_form-theta_TGP_exact)):.2e} deg  (IDENTYCZNE)")
print()
print(f"  Lepton: theta_formula = {np.degrees(theta_lep_form):.8f} deg")
print(f"          theta_exact   = {np.degrees(theta_lep_exact):.8f} deg")
print(f"          Roznica       = {abs(np.degrees(theta_lep_form-theta_lep_exact)):.2e} deg  (IDENTYCZNE)")
print()

# ============================================================
# KROK 2: Wrazliwosc theta na r21 i r31
# ============================================================
print("KROK 2: Wrazliwosc theta na r21 i r31")
print("-" * 65)

# dtheta/dr21 numerycznie
dr = 1e-4
dth_dr21 = (theta_from_ratios(r21_TGP+dr*r21_TGP, r31_TGP) -
            theta_from_ratios(r21_TGP-dr*r21_TGP, r31_TGP)) / (2*dr*r21_TGP)
dth_dr31 = (theta_from_ratios(r21_TGP, r31_TGP+dr*r31_TGP) -
            theta_from_ratios(r21_TGP, r31_TGP-dr*r31_TGP)) / (2*dr*r31_TGP)

print(f"  dtheta/dr21 = {np.degrees(dth_dr21)*1000:.4f} mdeg/unit_r21")
print(f"  dtheta/dr31 = {np.degrees(dth_dr31)*1000:.6f} mdeg/unit_r31")
print()

# Ile odchylenia theta pochodzi z Delta_r21 vs Delta_r31?
dr21_actual = r21_TGP - r21_lep  # 207.000 - 206.768 = +0.232
dr31_actual = r31_TGP - r31_lep  # 3477.10 - 3477.23 = -0.13

dtheta_from_r21 = dth_dr21 * dr21_actual
dtheta_from_r31 = dth_dr31 * dr31_actual
dtheta_total = theta_TGP_exact - theta_lep_exact

print(f"  TGP vs lepton:")
print(f"    Delta_r21 = {dr21_actual:+.4f}  -> Delta_theta = {np.degrees(dtheta_from_r21)*1000:+.3f} mdeg")
print(f"    Delta_r31 = {dr31_actual:+.4f}  -> Delta_theta = {np.degrees(dtheta_from_r31)*1000:+.3f} mdeg")
print(f"    Suma lin.  = {np.degrees(dtheta_from_r21+dtheta_from_r31)*1000:+.3f} mdeg")
print(f"    Dokladna   = {np.degrees(dtheta_total)*1000:+.3f} mdeg")
print()
frac_r21 = abs(dtheta_from_r21)/abs(dtheta_total) * 100
frac_r31 = abs(dtheta_from_r31)/abs(dtheta_total) * 100
print(f"    r21 odpowiada za {frac_r21:.1f}% rozniacy theta TGP-lepton")
print(f"    r31 odpowiada za {frac_r31:.1f}% rozniacy theta TGP-lepton")
print()

# ============================================================
# KROK 3: Kompletne rownanie masy TGP
# ============================================================
print("KROK 3: Kompletne rownanie masy TGP")
print("=" * 65)
print()

C_K = 2.351  # z P24
alpha = 8.5616
a_gam = 0.040
lam   = 5.501e-6

K1_form = C_K * a_gam / (1 + alpha)
r21_form = r21_TGP  # ≈ 207 (wyznaczone przez V_mod)
r31_form = r31_TGP  # ≈ 3477 (wyznaczone przez lam*)
K2_form = K1_form * r21_form
K3_form = K1_form * r31_form
theta_form = theta_from_ratios(r21_form, r31_form)

print("  ROWNANIE MASY TGP:")
print()
print("  POZIOM 1 -- Skala masowa:")
print(f"    K*1 = C_K * a_gam / (1+alpha)")
print(f"         = {C_K:.4f} * {a_gam:.3f} / (1+{alpha:.4f})")
print(f"         = {K1_form:.7f}  [wynik P24]")
print()
print("  POZIOM 2 -- Hierarchia mas:")
print(f"    K*2 = r21 * K*1,   r21 = {r21_form:.3f}  [wyznaczone przez V_mod]")
print(f"    K*3 = r31 * K*1,   r31 = {r31_form:.2f} [wyznaczone przez lam*]")
print(f"    K*2 = {K2_form:.6f},  K*3 = {K3_form:.4f}")
print()
print("  POZIOM 3 -- Formula Koide:")
print(f"    theta = arccos([3/(1+sqrt(r21)+sqrt(r31))-1]/sqrt(2))")
print(f"          = arccos([3/{1+np.sqrt(r21_form)+np.sqrt(r31_form):.4f}-1]/sqrt(2))")
print(f"          = {np.degrees(theta_form):.6f} deg")
print()
print("  POZIOM 4 -- Masa jako Koide form:")
print(f"    K_n = A^2 * (1 + sqrt(2)*cos(theta + 2*pi*(n-1)/3))^2")
print(f"    A = sqrt(K1) * (1+sqrt(r21)+sqrt(r31))/3")
print(f"      = {np.sqrt(K1_form) * (1+np.sqrt(r21_form)+np.sqrt(r31_form))/3:.6f}")
print()

# Porownanie z leptonami
print("  POROWNANIE TGP vs LEPTONY:")
print(f"  {'Wielkosc':>20}  {'TGP':>14}  {'Leptony':>14}  {'Roznica':>14}")
print("  " + "-"*66)
print(f"  {'r21':>20}  {r21_TGP:>14.4f}  {r21_lep:>14.4f}  {r21_TGP-r21_lep:>+14.4f} ({(r21_TGP-r21_lep)/r21_lep*100:+.4f}%)")
print(f"  {'r31':>20}  {r31_TGP:>14.2f}  {r31_lep:>14.2f}  {r31_TGP-r31_lep:>+14.2f} ({(r31_TGP-r31_lep)/r31_lep*100:+.4f}%)")
print(f"  {'theta [deg]':>20}  {np.degrees(theta_TGP_exact):>14.6f}  {np.degrees(theta_lep_exact):>14.6f}  {np.degrees(theta_TGP_exact-theta_lep_exact):>+14.6f}")
print(f"  {'Q':>20}  {koide_Q_from_ratios(r21_TGP, r31_TGP):>14.6f}  {koide_Q_from_ratios(r21_lep,r31_lep):>14.6f}  {koide_Q_from_ratios(r21_TGP,r31_TGP)-koide_Q_from_ratios(r21_lep,r31_lep):>+14.6f}")
print()

# ============================================================
# KROK 4: Co wyznacza r21 i r31 w TGP?
# ============================================================
print("KROK 4: Fizyczne wyznaczniki r21 i r31 w TGP")
print("-" * 65)
print()
print("  r21 = K*2/K*1 ≈ 207:")
print("    K*2 wyznaczone przez: drugie zero g2(K)=0 w modelu Yukawa")
print("    K*1 wyznaczone przez: pierwsze zero g1(K)=0 w modelu Yukawa")
print("    Oba zera wynikaja z ksztaltu V_mod(phi) (M_eff = 1/sqrt(1+alpha))")
print("    BRAK prostej analitycznej formuly dla K*2 (nieliniowe ODE)")
print()
print("  r31 = K*3/K*1 ≈ 3477:")
print("    K*3 wyznaczone przez: trzecie zero g3(K)=0, zdominowane przez lam*")
print("    K*3 ≈ a_gam * (e/sqrt(2)) * exp(a_gam) / sqrt(lam*)   [P21b]")
print(f"    Dla a_gam=0.04, lam=5.501e-6: K*3 ≈ {0.04*(np.e/np.sqrt(2))*np.exp(0.04)/np.sqrt(5.501e-6):.4f}")
print(f"    K*3 (numeryczne) = {K3_TGP:.4f}")
print()

# Sprawdz czy r31 wyznacza lam* (albo lam* wyznacza r31)
# K*3 = C3 * 1/sqrt(lam*)
# r31 = K*3/K*1 = C3/sqrt(lam*) / (C_K * a_gam/(1+alpha))
# = C3*(1+alpha) / (C_K * a_gam * sqrt(lam*))
C3 = K3_TGP * np.sqrt(lam)
print(f"  K*3 * sqrt(lam*) = {C3:.5f}  [stala z P21b: 0.08003]")
print(f"  r31 = K*3/K*1 = {C3:.5f} / sqrt(lam*) / (C_K * a_gam/(1+alpha))")
print(f"       = {C3/(np.sqrt(lam) * C_K * a_gam/(1+alpha)):.2f}")
print()

# r31 jako funkcja alpha, a_gam, lam
print("  Analityczna formula dla r31:")
print(f"    r31 = K*3/K*1 = [C3/sqrt(lam*)] / [C_K * a_gam/(1+alpha)]")
print(f"        = C3 * (1+alpha) / (C_K * a_gam * sqrt(lam*))")
C3_val = 0.08003  # z P21b
r31_form_an = C3_val * (1+alpha) / (C_K * a_gam * np.sqrt(lam))
print(f"        = {C3_val:.5f} * {1+alpha:.4f} / ({C_K:.4f} * {a_gam:.3f} * sqrt({lam:.4e}))")
print(f"        = {r31_form_an:.2f}  [cf. numeryczne: {r31_TGP:.2f}]")
print()

# r21 analitycznie: trudniejsze, wymaga K*2
# Ale mozemy zauwazyc ze r21 ≈ 207 jest praktycznie stale w rodzinie
print("  r21 ≈ 207: stala charakterystyczna V_mod dla gamma=1, lambda<<1")
print("    (nie ma prostej zamknietej formuly; wymaga rozwiazania g2(K)=0)")
print()

# ============================================================
# KROK 5: Zmiana theta przy modyfikacji parametrow
# ============================================================
print("KROK 5: Jak zmienic theta przez modyfikacje parametrow TGP?")
print("-" * 65)

# theta jest okreslone przez r21 i r31
# dtheta = dtheta/dr21 * dr21 + dtheta/dr31 * dr31

print(f"  Aby przejsc theta: {np.degrees(theta_TGP_exact):.4f} deg -> {np.degrees(theta_lep_exact):.4f} deg")
print(f"  Delta_theta = {np.degrees(theta_lep_exact-theta_TGP_exact)*1000:.3f} mdeg")
print()
print(f"  Wymagana zmiana r21 (przy stalym r31): Dr21 = {-np.degrees(theta_TGP_exact-theta_lep_exact)/np.degrees(dth_dr21):.4f}")
target_r21 = r21_TGP + (-np.degrees(theta_TGP_exact-theta_lep_exact)/np.degrees(dth_dr21))
print(f"    r21* = {target_r21:.4f}  (lepton: {r21_lep:.4f})")
print()
print(f"  Wymagana zmiana r31 (przy stalym r21): Dr31 = {-np.degrees(theta_TGP_exact-theta_lep_exact)/np.degrees(dth_dr31):.1f}")
target_r31 = r31_TGP + (-np.degrees(theta_TGP_exact-theta_lep_exact)/np.degrees(dth_dr31))
print(f"    r31* = {target_r31:.2f}  (lepton: {r31_lep:.2f})")
print()

# ============================================================
# WYKRES
# ============================================================
print("Tworze wykres...")

fig, axes = plt.subplots(1, 3, figsize=(15, 5))
fig.suptitle(
    'P25: Rownanie masy TGP -- theta = arccos([3/(1+sqrt(r21)+sqrt(r31))-1]/sqrt(2))\n'
    'theta okreslone przez r21 i r31 (TGP vs Leptony)',
    fontsize=10, fontweight='bold')

# Panel 1: theta(r21) przy stalym r31
ax = axes[0]
r21_arr = np.linspace(150, 280, 300)
theta_arr = np.degrees([theta_from_ratios(r, r31_TGP) for r in r21_arr])
ax.plot(r21_arr, theta_arr, 'b-', lw=2, label=f'theta(r21), r31={r31_TGP:.0f}')
ax.axvline(r21_TGP, color='red', lw=2, linestyle='--', label=f'TGP r21={r21_TGP:.2f}')
ax.axvline(r21_lep, color='green', lw=2, linestyle=':', label=f'Lepton r21={r21_lep:.3f}')
ax.axhline(np.degrees(theta_TGP_exact), color='red', lw=1, linestyle='--', alpha=0.5)
ax.axhline(np.degrees(theta_lep_exact), color='green', lw=1, linestyle=':', alpha=0.5)
ax.set_xlabel('r21 = K*2/K*1'); ax.set_ylabel('theta [deg]')
ax.set_title('theta(r21)', fontsize=9)
ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

# Panel 2: theta(r31) przy stalym r21
ax = axes[1]
r31_arr = np.linspace(2000, 6000, 300)
theta_arr2 = np.degrees([theta_from_ratios(r21_TGP, r) for r in r31_arr])
ax.plot(r31_arr, theta_arr2, 'b-', lw=2, label=f'theta(r31), r21={r21_TGP:.1f}')
ax.axvline(r31_TGP, color='red', lw=2, linestyle='--', label=f'TGP r31={r31_TGP:.0f}')
ax.axvline(r31_lep, color='green', lw=2, linestyle=':', label=f'Lepton r31={r31_lep:.0f}')
ax.set_xlabel('r31 = K*3/K*1'); ax.set_ylabel('theta [deg]')
ax.set_title('theta(r31)  (mniejsza wrazliwosc)', fontsize=9)
ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

# Panel 3: Diagram Koide (schemat rownania masy)
ax = axes[2]
# Narysuj schemat hierarchii
levels = {'K*1': K1_TGP, 'K*2': K2_TGP, 'K*3': K3_TGP}
levels_lep = {'me': ME, 'mmu': MMU, 'mtau': MTAU}
ax.barh(['K*3\n(tau)', 'K*2\n(mu)', 'K*1\n(e)'],
        [K3_TGP/K3_TGP, K2_TGP/K3_TGP, K1_TGP/K3_TGP],
        color=['#d62728', '#1f77b4', '#2ca02c'], alpha=0.7)
ax.set_xscale('log')
ax.set_xlabel('K*_n / K*3')
ax.set_title(f'Hierarchia mas TGP\nr21={r21_TGP:.1f}, r31={r31_TGP:.0f}\ntheta={np.degrees(theta_TGP_exact):.3f} deg', fontsize=9)
ax.grid(True, alpha=0.3)
ax.text(1.1, 2, f'r21={r21_TGP:.1f}', fontsize=9, color='red')
ax.text(1.1, 1, f'√r21={np.sqrt(r21_TGP):.2f}', fontsize=9, color='blue')
ax.text(1.1, 0, f'r31={r31_TGP:.0f}', fontsize=9, color='green')

plt.tight_layout()
out = __file__.replace('.py', '.png')
plt.savefig(out, dpi=120, bbox_inches='tight')
plt.close()
print(f"Zapisano: {out}")

# ============================================================
# PODSUMOWANIE FINALNE -- ROWNANIE MASY TGP
# ============================================================
print()
print("=" * 65)
print("ROWNANIE MASY TGP -- KOMPLETNE PODSUMOWANIE")
print("=" * 65)
print()
print("  1. SKALA MASOWA:")
print(f"     K*1 = 2.351 * a_gam / (1+alpha)   [P24]")
print()
print("  2. HIERARCHIA:")
print(f"     K*n = K*1 * r_n1  (n=1,2,3)")
print(f"     r21 ≈ 207  [Yukawa + V_mod]")
print(f"     r31 ≈ 0.08003 * (1+alpha) / (C_K * a_gam * sqrt(lam*))  [P21b]")
print()
print("  3. KAT KOIDE (DOKLADNIE):")
print(f"     cos(theta) = [3/(1+sqrt(r21)+sqrt(r31)) - 1] / sqrt(2)")
print(f"     theta_TGP   = {np.degrees(theta_TGP_exact):.4f} deg")
print(f"     theta_lepton = {np.degrees(theta_lep_exact):.4f} deg")
print(f"     Roznica: {np.degrees(theta_TGP_exact-theta_lep_exact)*1000:.2f} mdeg  ({(r21_TGP-r21_lep)/r21_lep*100:.3f}% z r21)")
print()
print("  4. WZOR KOIDE (zbiorczy):")
print(f"     K_n* = A^2 * (1+sqrt(2)*cos(theta + 2pi*(n-1)/3))^2")
print(f"     A    = sqrt(K*1) * (1+sqrt(r21)+sqrt(r31))/3 = {np.sqrt(K1_TGP)*(1+np.sqrt(r21_TGP)+np.sqrt(r31_TGP))/3:.6f}")
print(f"     theta = {np.degrees(theta_TGP_exact):.4f} deg  (od parametrow TGP)")
print()
print("  5. SENS FIZYCZNY:")
print(f"     - Skala K*1 ∝ a_gam/(1+alpha): rozmiar solitonu / nieliniowosci kinetycznej")
print(f"     - r21 ≈ 207: stala V_mod przy gamma=1 (nie zalezna od alpha!)")
print(f"     - r31 = f(alpha, a_gam, lam*): trzecia generacja wyznaczona przez czlon lam*")
print(f"     - theta ≈ 132.74 deg: kąt Koide wyznaczony przez (r21, r31) dokladnie")
print()
print("GOTOWE: P25 zakonczone.")
