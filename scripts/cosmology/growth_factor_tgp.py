#!/usr/bin/env python3
"""
import sys
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
growth_factor_tgp.py — Wzrost struktur w TGP
==============================================

Oblicza:
1. Liniowy czynnik wzrostu D(a) z G_eff(k,a) TGP
2. Tempo wzrostu f(a) = d ln D / d ln a
3. Observable f*sigma_8(z) — porownanie z danymi BOSS/DESI
4. Gravitational slip eta = Phi_N / Psi (lensing vs dynamics)

Model TGP (eq:Geff-k z sek08_formalizm.tex):
  G_eff(k,a) = G_0 * (1 + 2*alpha_eff^2 / (1 + (a*m_eff/k)^2))

Parametry:
  - gamma = 12 * Lambda_obs (z dark energy)
  - m_sp^2 = gamma = 3*gamma - 2*beta  (beta = gamma)
  - m_eff^2 = m_cosmo^2 = -4*gamma (from U''(1) in cosmological potential)
    NOTE: m_eff^2 < 0 means slow-roll, not instability
  - For structure growth, use SPATIAL mass m_sp^2 = gamma > 0
  - alpha_eff = q*Phi0/(4*pi)  where q = 8*pi*G_0/c_0^2

Autor: Claudian (analiza TGP v5)
Data: 2026-03-15
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

# ============================================================
# Fizyczne stale i parametry kosmologiczne
# ============================================================
c0 = 2.998e8        # m/s
G0 = 6.674e-11      # m^3/(kg s^2)
H0 = 67.4e3 / 3.086e22  # 67.4 km/s/Mpc -> 1/s
hbar = 1.055e-34     # J s

# Kosmologia (Planck 2018 + DESI DR2-like)
Omega_m0 = 0.315
Omega_L0 = 1 - Omega_m0
sigma8_fid = 0.811   # sigma_8 fiducjalny (Planck)

# Lambda obserwacyjna
Lambda_obs = 3 * H0**2 * Omega_L0 / c0**2  # m^{-2}

# Parametry TGP
gamma_TGP = 12 * Lambda_obs  # z Λ_eff = γ/12
beta_TGP = gamma_TGP         # warunek prozniowy

# Masy
m_sp_sq = gamma_TGP                    # [m^{-2}], spatial mass squared
m_sp = np.sqrt(m_sp_sq)                # [m^{-1}]
m_sp_Mpc = m_sp * 3.086e22             # [Mpc^{-1}]

# Stala sprzezenia
# q = 8 pi G_0 / c_0^2
q_TGP = 8 * np.pi * G0 / c0**2

# tau_0 = 1/(c_0 * sqrt(gamma_TGP))
tau0 = 1 / (c0 * np.sqrt(gamma_TGP))

# Phi_0 — estymacja z tau_0 ~ 1/H_0
# tau_0 = sqrt(Phi_0) / (c_0 * sqrt(gamma)) => Phi_0 = (tau_0 * c_0 * sqrt(gamma))^2
# Jesli tau_0 ~ 1/H_0:
Phi0_estimate = (c0 * np.sqrt(gamma_TGP) / H0)**2

# alpha_eff = q * Phi0 / (4 pi) — coupling strength
alpha_eff = q_TGP * Phi0_estimate / (4 * np.pi)

print("=" * 60)
print("WZROST STRUKTUR W TGP")
print("=" * 60)
print(f"\nParametry kosmologiczne:")
print(f"  H_0 = {H0:.3e} s^-1 ({67.4:.1f} km/s/Mpc)")
print(f"  Omega_m = {Omega_m0:.3f}, Omega_L = {Omega_L0:.3f}")
print(f"  Lambda_obs = {Lambda_obs:.3e} m^-2")
print(f"\nParametry TGP:")
print(f"  gamma = beta = 12 * Lambda_obs = {gamma_TGP:.3e} m^-2")
print(f"  m_sp = sqrt(gamma) = {m_sp:.3e} m^-1 = {m_sp_Mpc:.3e} Mpc^-1")
print(f"  m_sp / H_0 * c_0 = {m_sp * c0 / H0:.2f}")
print(f"  tau_0 = 1/(c_0 sqrt(gamma)) = {tau0:.3e} s = {tau0 * H0:.2f} / H_0")
print(f"  Phi_0 (estymacja) = {Phi0_estimate:.3e}")
print(f"  alpha_eff = q * Phi_0 / (4 pi) = {alpha_eff:.3e}")
print(f"  q = 8 pi G_0 / c_0^2 = {q_TGP:.3e} m/kg")


# ============================================================
# Hubble parameter H(a) for flat LCDM background
# ============================================================
def H_of_a(a):
    """H(a) / H_0 for flat LCDM."""
    return H0 * np.sqrt(Omega_m0 / a**3 + Omega_L0)


def E_of_a(a):
    """E(a) = H(a) / H_0"""
    return np.sqrt(Omega_m0 / a**3 + Omega_L0)


# ============================================================
# G_eff(k, a) in TGP
# ============================================================
def G_eff_TGP(k_Mpc, a, alpha_eff_val=None, m_eff_Mpc=None):
    """
    G_eff(k, a) / G_0 in TGP.

    G_eff = G_0 * (1 + 2*alpha_eff^2 / (1 + (a*m_eff/k)^2))

    Parameters:
        k_Mpc: wavenumber in h/Mpc
        a: scale factor
        alpha_eff_val: coupling (dimensionless)
        m_eff_Mpc: effective mass in Mpc^{-1}

    Note: for structure growth, we use the spatial mass m_sp
    """
    if alpha_eff_val is None:
        alpha_eff_val = alpha_eff
    if m_eff_Mpc is None:
        m_eff_Mpc = m_sp_Mpc

    ratio = (a * m_eff_Mpc / k_Mpc)**2
    return 1 + 2 * alpha_eff_val**2 / (1 + ratio)


# ============================================================
# Linear growth equation with G_eff
# ============================================================
def growth_ode(ln_a, y, k_Mpc, alpha_val, m_Mpc):
    """
    ODE for linear growth factor D(a).

    d^2 D / dt^2 + 2 H dD/dt = 4 pi G_eff(k,a) rho_m D

    In terms of ln(a):
    D'' + (2 + H'/H) D' = 3/2 Omega_m(a) * G_eff/G_0 * D

    where ' = d/d(ln a), H'/H = d ln H / d ln a

    y = [D, D'] where D' = dD/d(ln a)
    """
    a = np.exp(ln_a)
    D, Dp = y

    E = E_of_a(a)
    # d ln H / d ln a
    dE_dlna = (-3/2 * Omega_m0 / a**3) / (2 * E**2)  # = -(3/2) Omega_m(a) / (2 E^2)

    Omega_m_a = Omega_m0 / (a**3 * E**2)

    G_ratio = G_eff_TGP(k_Mpc, a, alpha_val, m_Mpc)

    # Growth equation: D'' + (2 + d ln H / d ln a) D' = 3/2 Omega_m(a) G_eff/G_0 D
    Dpp = 3/2 * Omega_m_a * G_ratio * D - (2 + dE_dlna) * Dp

    return [Dp, Dpp]


def compute_growth_factor(k_Mpc, alpha_val=None, m_Mpc=None, a_range=(1e-3, 1.0)):
    """
    Compute the linear growth factor D(a) for a given scale k.
    Returns a_arr, D_arr, f_arr (f = d ln D / d ln a)
    """
    if alpha_val is None:
        alpha_val = alpha_eff
    if m_Mpc is None:
        m_Mpc = m_sp_Mpc

    # Initial conditions in matter domination: D ~ a, D' = dD/d(ln a) = a = D
    ln_a_span = (np.log(a_range[0]), np.log(a_range[1]))
    y0 = [a_range[0], a_range[0]]  # D(a_ini) = a_ini, D'(a_ini) = a_ini

    ln_a_eval = np.linspace(*ln_a_span, 2000)

    sol = solve_ivp(growth_ode, ln_a_span, y0, t_eval=ln_a_eval,
                    args=(k_Mpc, alpha_val, m_Mpc),
                    method='RK45', rtol=1e-10, atol=1e-12)

    if not sol.success:
        print(f"  Growth ODE failed for k={k_Mpc}: {sol.message}")
        return None, None, None

    a_arr = np.exp(sol.t)
    D_arr = sol.y[0]
    Dp_arr = sol.y[1]  # dD / d(ln a)

    # Normalize: D(a=1) = 1
    D0 = np.interp(0.0, sol.t, D_arr)  # D at ln(a)=0, i.e., a=1
    D_arr /= D0
    Dp_arr /= D0

    # f = d ln D / d ln a = D' / D
    f_arr = Dp_arr / D_arr

    return a_arr, D_arr, f_arr


# ============================================================
# GR (LCDM) growth for comparison
# ============================================================
def compute_growth_GR(a_range=(1e-3, 1.0)):
    """Compute D(a) in standard GR (G_eff = G_0)."""
    return compute_growth_factor(k_Mpc=1e10, alpha_val=0.0, m_Mpc=1.0, a_range=a_range)
    # k >> m effectively gives G_eff = G_0 when alpha = 0


# ============================================================
# f * sigma_8(z) observable
# ============================================================
def compute_f_sigma8(k_Mpc, alpha_val=None, m_Mpc=None):
    """
    Compute f * sigma_8(z) = f(a) * sigma_8 * D(a)

    Uses sigma_8 fiducial * D(a) as sigma_8(a).
    """
    a_arr, D_arr, f_arr = compute_growth_factor(k_Mpc, alpha_val, m_Mpc)
    if a_arr is None:
        return None, None

    z_arr = 1 / a_arr - 1
    fs8_arr = f_arr * sigma8_fid * D_arr

    return z_arr, fs8_arr


# ============================================================
# Observational data: f*sigma_8 measurements
# ============================================================
# BOSS DR12 + 6dFGS + VIPERS + WiggleZ + DESI-like
# (z, f*sigma_8, error)
obs_data = np.array([
    [0.067, 0.423, 0.055],  # 6dFGS
    [0.15,  0.490, 0.145],  # SDSS MGS
    [0.32,  0.427, 0.056],  # BOSS LOWZ
    [0.57,  0.426, 0.029],  # BOSS CMASS
    [0.60,  0.550, 0.120],  # WiggleZ
    [0.73,  0.437, 0.072],  # WiggleZ
    [0.80,  0.470, 0.080],  # VIPERS
    [0.85,  0.467, 0.060],  # BOSS eBOSS LRG
    [1.48,  0.462, 0.045],  # BOSS eBOSS Ly-alpha (approx)
])

# DESI DR2 (2025-2026 preliminary, approximate values)
desi_data = np.array([
    [0.30, 0.462, 0.035],
    [0.51, 0.453, 0.028],
    [0.71, 0.432, 0.025],
    [0.93, 0.447, 0.030],
    [1.32, 0.385, 0.038],
])


# ============================================================
# Gravitational slip
# ============================================================
def gravitational_slip(k_Mpc, a, alpha_val=None, m_Mpc=None):
    """
    Gravitational slip eta = Phi_N / Psi in TGP.

    In TGP, the lensing potential (Phi_N + Psi)/2 differs from
    the dynamical potential Psi because of the scalar field contribution.

    For pure scalar-tensor: eta = (1 + alpha_eff^2 / (1 + (a*m/k)^2)) / G_eff * G_0
    But in TGP specifically:
    Psi = potential felt by matter (from G_eff * delta_m)
    Phi_N = potential felt by light (from conformal metric)

    In the conformal sector of TGP (no disformal):
    Phi_N and Psi are related by the anisotropic stress of delta Phi.
    eta = Sigma/mu in parameterized modified gravity language.

    For TGP with A(Phi) metric:
    mu(k,a) = G_eff/G_0 = 1 + 2*alpha_eff^2 / (1 + (a*m/k)^2)
    Sigma(k,a) = 1 + alpha_eff^2 / (1 + (a*m/k)^2)

    eta = Sigma / mu = [1 + alpha^2/(1+x^2)] / [1 + 2*alpha^2/(1+x^2)]
    """
    if alpha_val is None:
        alpha_val = alpha_eff
    if m_Mpc is None:
        m_Mpc = m_sp_Mpc

    x2 = (a * m_Mpc / k_Mpc)**2
    correction = alpha_val**2 / (1 + x2)

    mu = 1 + 2 * correction
    Sigma = 1 + correction

    eta = Sigma / mu
    return eta, mu, Sigma


# ============================================================
# Main computation
# ============================================================
def main():
    print("\n" + "=" * 60)
    print("OBLICZENIA")
    print("=" * 60)

    # 1. Growth factor for several scales
    print("\n--- Czynnik wzrostu D(a) ---")

    scales = {
        'GR (LCDM)': {'alpha': 0.0, 'm': 1.0, 'color': 'k', 'ls': '-'},
    }

    # Scan over different alpha_eff values to show the effect
    # The TGP prediction depends on alpha_eff which depends on Phi_0
    # Let's scan alpha_eff from 0 to 0.1

    alpha_values = [0.0, 0.01, 0.03, 0.05, 0.1]
    k_ref = 0.1  # h/Mpc — typical for RSD measurements

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('Wzrost struktur w TGP: $G_{\\rm eff}(k,a)$ i $f\\sigma_8(z)$',
                 fontsize=14, fontweight='bold')

    # Panel 1: G_eff(k) at a=1 for different alpha
    ax = axes[0, 0]
    k_arr = np.logspace(-3, 1, 500)  # h/Mpc
    for alpha_val in alpha_values:
        G_arr = [G_eff_TGP(k, 1.0, alpha_val, m_sp_Mpc) for k in k_arr]
        label = r'$\alpha_{eff}$=' + f'{alpha_val:.2f}' if alpha_val > 0 else 'GR'
        ls = '-' if alpha_val > 0 else '--'
        ax.semilogx(k_arr, G_arr, ls, linewidth=1.5, label=label)
    ax.set_xlabel(r'$k$ [h/Mpc]')
    ax.set_ylabel(r'$G_{\rm eff}(k,a=1) / G_0$')
    ax.set_title(r'Efektywna stala grawitacyjna')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.axhline(y=1, color='k', linestyle=':', alpha=0.3)

    # Panel 2: D(a) for different alpha at fixed k
    ax = axes[0, 1]
    colors = ['k', 'blue', 'green', 'orange', 'red']
    for i, alpha_val in enumerate(alpha_values):
        a_arr, D_arr, f_arr = compute_growth_factor(
            k_ref, alpha_val, m_sp_Mpc)
        if a_arr is not None:
            label = r'$\alpha_{eff}$=' + f'{alpha_val:.2f}' if alpha_val > 0 else 'GR'
            ax.plot(a_arr, D_arr, color=colors[i], linewidth=1.5, label=label)
    ax.set_xlabel('a')
    ax.set_ylabel('D(a) / D(1)')
    ax.set_title(f'Czynnik wzrostu (k = {k_ref} h/Mpc)')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # Panel 3: f*sigma_8(z) — TGP vs data
    ax = axes[1, 0]

    for i, alpha_val in enumerate(alpha_values):
        z_arr, fs8_arr = compute_f_sigma8(k_ref, alpha_val, m_sp_Mpc)
        if z_arr is not None:
            mask = z_arr > 0
            label = r'$\alpha_{eff}$=' + f'{alpha_val:.2f}' if alpha_val > 0 else 'GR'
            ax.plot(z_arr[mask], fs8_arr[mask], color=colors[i],
                    linewidth=1.5, label=label)

    # Plot observational data
    ax.errorbar(obs_data[:, 0], obs_data[:, 1], yerr=obs_data[:, 2],
                fmt='s', color='gray', markersize=5, capsize=3,
                label='BOSS/6dF/VIPERS', zorder=5)
    ax.errorbar(desi_data[:, 0], desi_data[:, 1], yerr=desi_data[:, 2],
                fmt='D', color='purple', markersize=5, capsize=3,
                label='DESI DR2 (approx)', zorder=5)

    ax.set_xlabel('z')
    ax.set_ylabel(r'$f\sigma_8(z)$')
    ax.set_title(r'Tempo wzrostu struktur')
    ax.set_xlim(0, 2)
    ax.set_ylim(0.2, 0.7)
    ax.legend(fontsize=7, ncol=2)
    ax.grid(True, alpha=0.3)

    # Panel 4: Gravitational slip
    ax = axes[1, 1]
    z_slip = np.linspace(0, 2, 200)
    a_slip = 1 / (1 + z_slip)

    for i, alpha_val in enumerate(alpha_values):
        if alpha_val == 0:
            ax.axhline(y=1, color='k', linestyle='--', linewidth=1.5, label='GR')
            continue
        eta_arr = [gravitational_slip(k_ref, a, alpha_val, m_sp_Mpc)[0]
                   for a in a_slip]
        label = r'$\alpha_{eff}$=' + f'{alpha_val:.2f}'
        ax.plot(z_slip, eta_arr, color=colors[i], linewidth=1.5, label=label)

    ax.set_xlabel('z')
    ax.set_ylabel(r'$\eta = \Sigma / \mu$')
    ax.set_title('Gravitational slip (lensing/dynamics)')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0.9, 1.01)

    plt.tight_layout()

    plots_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'plots')
    os.makedirs(plots_dir, exist_ok=True)
    outpath = os.path.join(plots_dir, 'growth_factor_tgp.png')
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    print(f"\n  Wykres zapisany: {outpath}")
    plt.close()

    # =========================================================
    # Chi-square analysis: which alpha_eff best fits the data?
    # =========================================================
    print("\n--- Chi-square analiza alpha_eff ---")

    all_data = np.vstack([obs_data, desi_data])
    alpha_scan = np.linspace(0, 0.15, 100)
    chi2_arr = []

    for alpha_val in alpha_scan:
        z_model, fs8_model = compute_f_sigma8(k_ref, alpha_val, m_sp_Mpc)
        if z_model is None:
            chi2_arr.append(np.inf)
            continue

        # Interpolate model at data redshifts
        mask = z_model > 0
        z_m = z_model[mask][::-1]  # reverse to make increasing
        fs8_m = fs8_model[mask][::-1]

        if z_m[0] > all_data[:, 0].min() or z_m[-1] < all_data[:, 0].max():
            chi2_arr.append(np.inf)
            continue

        interp = interp1d(z_m, fs8_m, kind='linear', fill_value='extrapolate')
        model_at_data = interp(all_data[:, 0])

        chi2 = np.sum(((all_data[:, 1] - model_at_data) / all_data[:, 2])**2)
        chi2_arr.append(chi2)

    chi2_arr = np.array(chi2_arr)
    best_idx = np.argmin(chi2_arr)
    best_alpha = alpha_scan[best_idx]
    best_chi2 = chi2_arr[best_idx]
    ndof = len(all_data) - 1

    print(f"  Najlepsze alpha_eff = {best_alpha:.4f}")
    print(f"  chi^2 / ndof = {best_chi2:.1f} / {ndof} = {best_chi2/ndof:.2f}")
    print(f"  chi^2(GR) / ndof = {chi2_arr[0]:.1f} / {ndof} = {chi2_arr[0]/ndof:.2f}")

    # 1-sigma bound
    delta_chi2 = chi2_arr - best_chi2
    mask_1sigma = delta_chi2 < 1
    if np.any(mask_1sigma):
        alpha_1sigma_max = alpha_scan[mask_1sigma][-1]
        print(f"  1-sigma gorny limit: alpha_eff < {alpha_1sigma_max:.4f}")

    # Plot chi-square
    fig2, ax2 = plt.subplots(figsize=(8, 5))
    ax2.plot(alpha_scan, chi2_arr - best_chi2, 'b-', linewidth=2)
    ax2.axhline(y=1, color='r', linestyle='--', label=r'$\Delta\chi^2 = 1$')
    ax2.axhline(y=4, color='orange', linestyle='--', label=r'$\Delta\chi^2 = 4$')
    ax2.set_xlabel(r'$\alpha_{\rm eff}$')
    ax2.set_ylabel(r'$\Delta\chi^2$')
    ax2.set_title(r'Ograniczenie na $\alpha_{\rm eff}$ z $f\sigma_8$ danych')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim(-0.5, 10)

    outpath2 = os.path.join(plots_dir, 'alpha_eff_chi2.png')
    plt.savefig(outpath2, dpi=150, bbox_inches='tight')
    print(f"  Wykres chi^2 zapisany: {outpath2}")
    plt.close()

    # =========================================================
    # Summary
    # =========================================================
    print("\n" + "=" * 60)
    print("PODSUMOWANIE")
    print("=" * 60)
    print(f"""
Kluczowe wyniki:

1. G_eff(k,a) w TGP:
   - Ma masy: k >> a*m_sp => G_eff = G_0(1 + 2*alpha_eff^2)
   - k << a*m_sp => G_eff = G_0 (standardowa)
   - Skala przejscia: k ~ a * m_sp = {m_sp_Mpc:.2e} Mpc^-1

2. f*sigma_8(z):
   - GR: chi^2/ndof = {chi2_arr[0]:.1f}/{ndof} = {chi2_arr[0]/ndof:.2f}
   - Best TGP: alpha_eff = {best_alpha:.4f}, chi^2/ndof = {best_chi2:.1f}/{ndof}

3. Gravitational slip:
   - eta = Sigma/mu < 1 zawsze (lensing slabsze niz dynamika)
   - Dla alpha_eff = 0.05: eta ~ {gravitational_slip(k_ref, 1.0, 0.05, m_sp_Mpc)[0]:.4f}

4. Predykcje TGP:
   - Wzmocniona grawitacja na malych skalach (k > a*m_sp)
   - Gravitational slip eta < 1 (testowalne: Euclid, LSST)
   - Skalo-zalezny wzrost f*sigma_8 (testowalne: DESI, Euclid)

5. UWAGA: alpha_eff ~ {alpha_eff:.1e} (estymacja z Phi_0)
   Jesli alpha_eff << 0.01, TGP jest nieodroznialne od GR
   w sektorze wzrostu struktur.
""")


if __name__ == '__main__':
    main()
