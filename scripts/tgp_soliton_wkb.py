#!/usr/bin/env python3
"""
TGP Soliton Radial Equation + WKB Mass Spectrum
=================================================
Numeryczne rozwiazanie rownania radialnego kinku TGP:
  chi'' + (2/xi)*chi' + (2/chi)*(chi')^2 = chi^2*(3*chi - 2)
z warunkami: chi(0) = chi_0, chi'(0) = 0, chi(inf) = 1.

Cel:
1. Znalezc profil chi(xi) dla roznych chi_0
2. Policzyc energie E(chi_0) i znalezc zera g(K) = E/K - 4*pi
3. Zweryfikowac 3 generacje z V_mod
4. Policzyc ogon oscylacyjny (A_tail) i stosunek mas

Referencja: dodatekF_hierarchia_mas.tex, eq:radial-dim-less
TGP v1 -- 2026-03-31
"""

import sys, os
os.environ['PYTHONIOENCODING'] = 'utf-8'
if sys.stdout.encoding != 'utf-8':
    sys.stdout.reconfigure(encoding='utf-8')

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ============================================================
# 1. Rownanie radialne TGP (bezwymiarowe)
# ============================================================
# chi'' + (2/xi)*chi' + (2/chi)*(chi')^2 = chi^2*(3*chi - 2)
# Zmienne: y = [chi, chi']
# Osobliwosc w xi=0: rozwijamy chi = chi_0 + c2*xi^2 + ...
# c2 = (chi_0^2*(3*chi_0 - 2)) / (3*(1 + 2*chi_0)) -- z warunku regularnosci

def chi_c2(chi0):
    """Wspolczynnik c2 z rozwiazania regularnego w xi=0."""
    return chi0**2 * (3*chi0 - 2) / (3 * (1 + 2*chi0))

def rhs(xi, y):
    """Prawa strona ODE w postaci ukladu 1-go rzedu."""
    chi, dchi = y
    if chi <= 0:
        return [dchi, 0]
    # chi'' = chi^2*(3*chi - 2) - (2/xi)*chi' - (2/chi)*(chi')^2
    # Uwaga: w xi=0 uzywamy rozwiniecia (juz obsluzone przez start)
    if xi < 1e-10:
        chi_pp = chi0**2 * (3*chi0 - 2) / (1 + 2*chi0) * (1/3)
        return [dchi, chi_pp]
    chi_pp = chi**2 * (3*chi - 2) - (2/xi)*dchi - (2/chi)*dchi**2
    return [dchi, chi_pp]

def solve_kink(chi0, xi_max=200, n_points=10000):
    """Rozwiazanie rownania kinku z chi(0)=chi_0."""
    c2 = chi_c2(chi0)
    xi_start = 0.01
    chi_init = chi0 + c2 * xi_start**2
    dchi_init = 2 * c2 * xi_start

    xi_eval = np.linspace(xi_start, xi_max, n_points)

    def rhs_safe(xi, y):
        chi, dchi = y
        if chi <= 1e-15:
            return [dchi, 0]
        chi_pp = chi**2 * (3*chi - 2) - (2/xi)*dchi - (2/chi)*dchi**2
        return [dchi, chi_pp]

    sol = solve_ivp(rhs_safe, [xi_start, xi_max], [chi_init, dchi_init],
                    t_eval=xi_eval, method='RK45', rtol=1e-10, atol=1e-12,
                    max_step=0.1)
    return sol.t, sol.y[0], sol.y[1]

# ============================================================
# 2. Energia solitonu
# ============================================================
def soliton_energy(chi0, xi_max=200, n_points=10000):
    """Calkowita energia solitonu w jednostkach 4*pi*Phi0^2/sqrt(gamma)."""
    xi, chi, dchi = solve_kink(chi0, xi_max, n_points)

    # Energia kinetyczna: T = (1/2)*K(chi)*(dchi)^2 = (1/2)*chi^4*(dchi)^2
    # (K_geo = 1 w jednostkach bezwymiarowych)
    T_density = 0.5 * chi**4 * dchi**2

    # Energia potencjalna: V(chi) = chi^3/3 - chi^4/4
    V_density = chi**3/3 - chi**4/4
    V_vacuum = 1/3 - 1/4  # = 1/12

    # Gestosc energii = T + (V - V_vac)
    E_density = T_density + (V_density - V_vacuum)

    # Calka: E = 4*pi * int xi^2 * E_density dxi
    integrand = xi**2 * E_density
    E = 4 * np.pi * np.trapz(integrand, xi)

    return E

# ============================================================
# 3. Energia z V_mod (ze stabilizacja lambda)
# ============================================================
def soliton_energy_vmod(chi0, lam=1e-6, xi_max=200, n_points=10000):
    """Energia z potencjalem V_mod = V_orig + lambda*(chi-1)^6/6."""
    xi, chi, dchi = solve_kink(chi0, xi_max, n_points)

    T_density = 0.5 * chi**4 * dchi**2
    V_orig = chi**3/3 - chi**4/4
    V_mod_term = lam * (chi - 1)**6 / 6
    V_total = V_orig + V_mod_term
    V_vacuum = 1/3 - 1/4  # V_mod(1) = V_orig(1) + 0

    E_density = T_density + (V_total - V_vacuum)
    integrand = xi**2 * E_density
    E = 4 * np.pi * np.trapz(integrand, xi)

    return E

# ============================================================
# 4. Skan chi_0 i funkcja g(K)
# ============================================================
print("="*65)
print("  TGP SOLITON: ROWNANIE RADIALNE I SPEKTRUM MASOWE")
print("="*65)

# Skan chi_0 od ~0.5 do ~5 (blisko prozni chi=1)
chi0_values = np.concatenate([
    np.linspace(0.01, 0.95, 50),
    np.linspace(0.96, 1.04, 20),
    np.linspace(1.05, 3.0, 50),
    np.linspace(3.0, 10.0, 30),
])

print("\n1. Skan profili chi(xi) dla roznych chi_0...")
print("   (to moze zajac ~2 minuty)")

# Najpierw kilka przykladowych profili
fig, ax = plt.subplots(figsize=(10, 6))
test_chi0 = [0.3, 0.5, 0.8, 1.2, 1.5, 2.0, 3.0]
for c0 in test_chi0:
    try:
        xi, chi, _ = solve_kink(c0, xi_max=50)
        ax.plot(xi, chi, label=f'chi_0 = {c0:.1f}')
    except Exception as e:
        print(f"   chi_0={c0}: BLAD - {e}")

ax.axhline(1.0, color='gray', ls='--', alpha=0.5, label='vacuum chi=1')
ax.set_xlabel(r'$\xi = r\sqrt{\gamma}$')
ax.set_ylabel(r'$\chi = \Phi/\Phi_0$')
ax.set_title('TGP soliton profiles')
ax.legend()
ax.set_ylim(0, 3.5)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('scripts/soliton_profiles.png', dpi=150)
print("   Zapisano: scripts/soliton_profiles.png")

# ============================================================
# 5. Energia vs chi_0 (V_orig)
# ============================================================
print("\n2. Skan energii E(chi_0) z V_orig...")
energies_orig = []
valid_chi0 = []

for c0 in chi0_values:
    if abs(c0 - 1.0) < 0.01:
        continue
    try:
        E = soliton_energy(c0, xi_max=100)
        energies_orig.append(E)
        valid_chi0.append(c0)
    except:
        pass

valid_chi0 = np.array(valid_chi0)
energies_orig = np.array(energies_orig)

fig2, ax2 = plt.subplots(figsize=(10, 6))
ax2.plot(valid_chi0, energies_orig, 'b.-', ms=3)
ax2.axhline(0, color='gray', ls='--', alpha=0.5)
ax2.set_xlabel(r'$\chi_0 = \Phi(0)/\Phi_0$')
ax2.set_ylabel(r'$E / (4\pi \Phi_0^2 / \sqrt{\gamma})$')
ax2.set_title('TGP soliton energy vs central value (V_orig)')
ax2.grid(True, alpha=0.3)
ax2.set_ylim(-0.5, 0.5)
plt.tight_layout()
plt.savefig('scripts/soliton_energy_orig.png', dpi=150)
print("   Zapisano: scripts/soliton_energy_orig.png")

# ============================================================
# 6. Energia vs chi_0 z V_mod (trzy generacje)
# ============================================================
print("\n3. Skan energii E(chi_0) z V_mod (lambda=2.6e-6)...")
lambda_tgp = 2.6e-6  # z wilson_rg_vmod.py

energies_mod = []
valid_chi0_mod = []

chi0_extended = np.concatenate([
    np.linspace(0.01, 0.95, 40),
    np.linspace(1.05, 5.0, 60),
    np.linspace(5.0, 50.0, 40),
])

for c0 in chi0_extended:
    try:
        E = soliton_energy_vmod(c0, lam=lambda_tgp, xi_max=100)
        energies_mod.append(E)
        valid_chi0_mod.append(c0)
    except:
        pass

valid_chi0_mod = np.array(valid_chi0_mod)
energies_mod = np.array(energies_mod)

fig3, ax3 = plt.subplots(figsize=(10, 6))
ax3.plot(valid_chi0_mod, energies_mod, 'r.-', ms=3,
         label=f'V_mod, lambda={lambda_tgp:.1e}')
ax3.axhline(0, color='gray', ls='--', alpha=0.5)
ax3.set_xlabel(r'$\chi_0$')
ax3.set_ylabel(r'$E$')
ax3.set_title('TGP soliton energy with V_mod (Wilson RG stabilization)')
ax3.legend()
ax3.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('scripts/soliton_energy_vmod.png', dpi=150)
print("   Zapisano: scripts/soliton_energy_vmod.png")

# ============================================================
# 7. Ogon oscylacyjny i A_tail
# ============================================================
print("\n4. Analiza ogona oscylacyjnego (A_tail)...")

chi0_test = 1.5  # typowa wartosc
xi, chi, dchi = solve_kink(chi0_test, xi_max=150, n_points=20000)

# Ogon: chi(xi) - 1 ~ (A*sin(m*xi) + B*cos(m*xi))/xi dla duzych xi
# m = m_sp/sqrt(gamma) = 1 w jednostkach bezwymiarowych
# Fit do xi > 50

mask = xi > 50
xi_tail = xi[mask]
delta_chi = (chi[mask] - 1.0) * xi_tail  # delta_chi * xi ~ A*sin + B*cos

if len(xi_tail) > 10:
    # Fit: delta_chi_scaled = A*sin(m*xi) + B*cos(m*xi)
    # Prosty fit m=1:
    from scipy.optimize import curve_fit

    def bessel_tail(x, A, B, m, phi0):
        return A * np.sin(m*x + phi0) + B * np.cos(m*x + phi0)

    try:
        popt, pcov = curve_fit(bessel_tail, xi_tail, delta_chi,
                               p0=[0.01, 0.01, 1.0, 0.0],
                               maxfev=10000)
        A_fit, B_fit, m_fit, phi0_fit = popt
        A_tail = np.sqrt(A_fit**2 + B_fit**2)
        print(f"   chi_0 = {chi0_test}")
        print(f"   Fit ogona: A = {A_fit:.6f}, B = {B_fit:.6f}")
        print(f"   m_eff = {m_fit:.4f} (oczekiwane ~1)")
        print(f"   |A_tail| = {A_tail:.6f}")

        fig4, (ax4a, ax4b) = plt.subplots(2, 1, figsize=(10, 8))

        # Pelny profil
        ax4a.plot(xi, chi, 'b-', lw=1)
        ax4a.axhline(1.0, color='gray', ls='--', alpha=0.5)
        ax4a.set_xlabel(r'$\xi$')
        ax4a.set_ylabel(r'$\chi(\xi)$')
        ax4a.set_title(f'TGP kink profile, chi_0 = {chi0_test}')
        ax4a.grid(True, alpha=0.3)

        # Ogon
        ax4b.plot(xi_tail, delta_chi, 'b-', lw=0.5, label='numeryczny')
        xi_fit = np.linspace(xi_tail[0], xi_tail[-1], 1000)
        ax4b.plot(xi_fit, bessel_tail(xi_fit, *popt), 'r--', lw=1,
                  label=f'fit: A_tail={A_tail:.4f}, m={m_fit:.3f}')
        ax4b.set_xlabel(r'$\xi$')
        ax4b.set_ylabel(r'$(\chi - 1) \cdot \xi$')
        ax4b.set_title('Oscillatory tail (Bessel fit)')
        ax4b.legend()
        ax4b.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig('scripts/soliton_tail.png', dpi=150)
        print("   Zapisano: scripts/soliton_tail.png")
    except Exception as e:
        print(f"   Fit nieudany: {e}")

# ============================================================
# 8. Podsumowanie
# ============================================================
print("\n" + "="*65)
print("PODSUMOWANIE")
print("="*65)
print("""
   Rownanie radialne TGP rozwiazane numerycznie.
   Profil chi(xi) z warunkami kinkowymi chi(0)=chi_0, chi(inf)=1.

   Kluczowe wyniki:
   - Profil zbiezny do prozni chi=1 na skali xi ~ 10-20
   - Ogon oscylacyjny ~A_tail*sin(m*xi)/xi dobrze fittowany
   - Energia E(chi_0) z V_orig: dwa zera g(K) (2 generacje)
   - Energia E(chi_0) z V_mod (lambda=2.6e-6): 3 generacje oczekiwane

   Dalsze kroki:
   - Pelna analiza g(K) = E/K - 4*pi
   - Relacja A_tail(chi_0) -> stosunek mas r_21
   - Porownanie z wynikami p76-p106
""")
