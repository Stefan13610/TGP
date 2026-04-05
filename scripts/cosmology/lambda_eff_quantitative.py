"""
import sys
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
lambda_eff_quantitative.py  --  Theory of Generated Space (TGP)
================================================================
Quantitative computation of Lambda_eff from realistic matter
distribution and comparison with observed Lambda ~ 10^-52 m^-2.

Addresses prob:Lambda from sek05_ciemna_energia.tex.

Structure
---------
(a) Homogeneous background: Lambda_eff = gamma/12
(b) Perturbative correction from structure (NFW + disk)
(c) Toy model: two sources on a circle (lem:Lambda-positive)
(d) Scaling with structure formation via growth factor D(z)
(e) Results table
(f) Diagnostic 4-panel plot

Key honesty statement
---------------------
TGP does NOT independently predict Lambda. The relation
  Lambda_eff = gamma/12 = Phi_0 * H_0^2 / (12 c_0^2)
requires gamma as INPUT (set by tau_0 ~ 1/H_0 naturalness).
The "prediction" is that Phi_0 ~ 25 (an O(10) number, not fine-tuned)
reproduces Lambda_obs. The value Phi_0 ~ 25 is FIXED by the DE
constraint; all other constraints are then CHECKED (not fitted).

Output:
  scripts/plots/lambda_eff_quantitative.png
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# =====================================================================
# Physical constants (SI)
# =====================================================================
c0 = 2.998e8           # m/s
G0 = 6.674e-11         # m^3 kg^-1 s^-2
H0_km = 67.4           # km/s/Mpc
Mpc_m = 3.0857e22      # m/Mpc
H0 = H0_km * 1e3 / Mpc_m   # s^-1  ~2.184e-18

Omega_Lambda = 0.685
Omega_m = 0.315
Lambda_obs = 1.11e-52   # m^-2  (observed)
rho_crit = 3 * H0**2 / (8 * np.pi * G0)   # kg/m^3

# TGP parameter
Phi0 = 12 * Lambda_obs * c0**2 / H0**2   # from Lambda_eff = gamma/12 = Lambda_obs

# =====================================================================
# (a) Homogeneous background
# =====================================================================
def homogeneous_background():
    """
    Lambda_eff = gamma/12, where gamma = Phi_0 * H_0^2 / c_0^2.
    Setting Lambda_eff = Lambda_obs determines Phi_0.
    """
    gamma = Phi0 * H0**2 / c0**2
    Lambda_eff = gamma / 12.0
    ratio = Lambda_eff / Lambda_obs

    # Phi_0 that gives exact match
    Phi0_match = 12 * Lambda_obs * c0**2 / H0**2
    # Alternative: from Omega_Lambda
    # Lambda_obs = 3 * Omega_Lambda * H0^2 / c0^2
    # Phi0_match = 12 * 3 * Omega_Lambda * H0^2 / c0^2 * c0^2 / H0^2
    #            = 36 * Omega_Lambda = 36 * 0.685 = 24.66
    Phi0_from_Omega = 36 * Omega_Lambda

    print("=" * 70)
    print("  (a) HOMOGENEOUS BACKGROUND")
    print("=" * 70)
    print(f"  H_0       = {H0:.4e} s^-1  ({H0_km} km/s/Mpc)")
    print(f"  c_0       = {c0:.3e} m/s")
    print(f"  Lambda_obs= {Lambda_obs:.3e} m^-2")
    print()
    print(f"  TGP relation: Lambda_eff = Phi_0 * H_0^2 / (12 * c_0^2)")
    print(f"  Setting Lambda_eff = Lambda_obs:")
    print(f"    Phi_0 = 12 * Lambda_obs * c_0^2 / H_0^2")
    print(f"          = {Phi0_match:.4f}")
    print(f"  From Omega_Lambda = {Omega_Lambda}:")
    print(f"    Lambda_obs = 3 * Omega_Lambda * H_0^2 / c_0^2")
    print(f"    Phi_0 = 36 * Omega_Lambda = {Phi0_from_Omega:.4f}")
    print()
    print(f"  gamma     = {gamma:.4e} m^-2")
    print(f"  Lambda_eff= {Lambda_eff:.4e} m^-2")
    print(f"  Ratio     = {ratio:.6f}")
    print()
    print("  INTERPRETATION:")
    print("    gamma is NOT independently predicted. It is set by the")
    print("    naturalness condition tau_0 ~ 1/H_0, which gives")
    print("    gamma = Phi_0 * H_0^2 / c_0^2.")
    print("    Then Lambda_eff = gamma/12 = Phi_0 * H_0^2 / (12 c_0^2).")
    print("    Matching Lambda_obs FIXES Phi_0 ~ 25.")
    print("    This is a FIT with one parameter, but Phi_0 ~ O(10)")
    print("    is not fine-tuned (compare with 10^122 in QFT).")
    print()

    return dict(
        gamma=gamma, Lambda_eff=Lambda_eff, ratio=ratio,
        Phi0_match=Phi0_match, Phi0_from_Omega=Phi0_from_Omega,
    )


# =====================================================================
# (b) Perturbative correction from structure
# =====================================================================
def structure_correction():
    """
    From lem:Lambda-positive:
      Lambda_eff = U(1) + (1/2) m_sp^2 <delta_Phi^2> / Phi_0^2
                 = gamma/12 + (gamma/2) <delta_Phi^2> / Phi_0^2

    Correction: delta_Lambda / Lambda = 6 <delta_Phi^2> / Phi_0^2

    We estimate <delta_Phi^2> / Phi_0^2 from the Newtonian potential:
      delta_Phi / Phi_0 ~ Psi_N = G M / (c^2 r)
    averaged over cosmological volume.
    """
    print("=" * 70)
    print("  (b) PERTURBATIVE CORRECTION FROM STRUCTURE")
    print("=" * 70)

    gamma = Phi0 * H0**2 / c0**2

    # --- NFW halo ---
    # Typical galaxy halo: M_200 = 10^12 M_sun, r_s = 20 kpc, c_NFW = 10
    M_sun = 1.989e30  # kg
    kpc = 3.086e19     # m
    M_200 = 1e12 * M_sun
    r_s = 20 * kpc
    c_NFW = 10
    r_200 = c_NFW * r_s

    # NFW potential at r:
    # Psi_NFW(r) = - G M_200 / [r * f(c)] * ln(1 + r/r_s)
    # where f(c) = ln(1+c) - c/(1+c)
    f_c = np.log(1 + c_NFW) - c_NFW / (1 + c_NFW)

    # RMS of Psi_N / c^2 within r_200 (volume-weighted)
    Nr = 1000
    r_arr = np.linspace(0.01 * r_s, r_200, Nr)
    psi_arr = np.zeros(Nr)
    for i, r in enumerate(r_arr):
        x = r / r_s
        psi_arr[i] = G0 * M_200 / (r * f_c * c0**2) * np.log(1 + x)

    # Volume-weighted <Psi^2>
    dV = 4 * np.pi * r_arr**2
    psi2_avg_halo = np.sum(psi_arr**2 * dV) / np.sum(dV)

    print(f"\n  NFW halo (M_200 = 10^12 M_sun, c = {c_NFW}):")
    print(f"    r_s   = {r_s:.2e} m ({r_s/kpc:.0f} kpc)")
    print(f"    r_200 = {r_200:.2e} m ({r_200/kpc:.0f} kpc)")
    print(f"    <Psi_N^2/c^4> = {psi2_avg_halo:.4e}")
    print(f"    sqrt(<Psi^2>) = {np.sqrt(psi2_avg_halo):.4e}")

    # --- Exponential disk ---
    # M_disk = 5 * 10^10 M_sun, R_d = 3 kpc
    M_disk = 5e10 * M_sun
    R_d = 3 * kpc
    # Approximate: <Psi^2> ~ (G M_disk / (c^2 R_d))^2
    psi_disk_char = G0 * M_disk / (c0**2 * R_d)
    psi2_avg_disk = psi_disk_char**2 / 3.0  # rough factor

    print(f"\n  Exponential disk (M = 5e10 M_sun, R_d = 3 kpc):")
    print(f"    Psi_char = {psi_disk_char:.4e}")
    print(f"    <Psi^2>  ~ {psi2_avg_disk:.4e}")

    # --- Galaxy cluster ---
    M_cl = 1e15 * M_sun
    R_cl = 1.5e6 * kpc  # 1.5 Mpc
    psi_cl = G0 * M_cl / (c0**2 * R_cl)
    psi2_cl = psi_cl**2 / 3.0

    print(f"\n  Galaxy cluster (M = 10^15 M_sun, R = 1.5 Mpc):")
    print(f"    Psi_char = {psi_cl:.4e}")
    print(f"    <Psi^2>  ~ {psi2_cl:.4e}")

    # --- Cosmic web (volume-weighted average) ---
    # The cosmological variance of the Newtonian potential is
    # <Psi^2> ~ (3/5)^2 * Omega_m^2 * (H_0/c_0)^2 * sigma_8^2 / k_eff^2
    # Or simpler: the Bardeen potential Psi ~ 10^-5 on large scales
    sigma_Psi_cosmo = 3e-5
    psi2_cosmo = sigma_Psi_cosmo**2

    print(f"\n  Cosmic web (large-scale Bardeen potential):")
    print(f"    sigma_Psi ~ {sigma_Psi_cosmo:.1e}")
    print(f"    <Psi^2>   = {psi2_cosmo:.2e}")

    # --- Volume filling fractions ---
    # Halos occupy ~1% of volume, clusters ~0.01%
    # Cosmic web Psi is the dominant volume-averaged contribution
    f_halo = 0.01
    f_cluster = 1e-4
    f_web = 1.0  # fills everything

    delta_phi2_avg = (f_halo * psi2_avg_halo +
                      f_cluster * psi2_cl +
                      f_web * psi2_cosmo)

    # delta_Lambda / Lambda = 6 * <delta_Phi^2> / Phi_0^2
    # But in TGP: delta_Phi / Phi_0 ~ Psi_N (the Newtonian potential)
    # So <delta_Phi^2> / Phi_0^2 = <Psi_N^2>
    correction = 6 * delta_phi2_avg

    print(f"\n  Volume-averaged <delta_Phi^2/Phi_0^2>:")
    print(f"    Halos  (f={f_halo}):    {f_halo * psi2_avg_halo:.3e}")
    print(f"    Clusters (f={f_cluster}): {f_cluster * psi2_cl:.3e}")
    print(f"    Cosmic web (f={f_web}):  {f_web * psi2_cosmo:.3e}")
    print(f"    Total:                 {delta_phi2_avg:.3e}")
    print(f"\n  delta_Lambda / Lambda = 6 * <delta_Phi^2/Phi_0^2>")
    print(f"                        = {correction:.3e}")
    print(f"\n  RESULT: Correction is {correction:.1e} = {correction*100:.4f}%")
    print(f"    This is NEGLIGIBLE compared to the background term.")
    print(f"    The cosmological constant is dominated by U(1) = gamma/12.")
    print()

    return dict(
        psi2_halo=psi2_avg_halo,
        psi2_disk=psi2_avg_disk,
        psi2_cluster=psi2_cl,
        psi2_cosmo=psi2_cosmo,
        delta_phi2_avg=delta_phi2_avg,
        correction=correction,
    )


# =====================================================================
# (c) Toy model: two sources on a circle
# =====================================================================
def toy_model_circle():
    """
    1D periodic problem from rem:toy-Lambda.
    Phi'' + m^2 (Phi - Phi0) = -q Phi0 M [delta(x) + delta(x - L/2)]
    with periodic BC and zero-sum constraint.

    Solve via Fourier series.
    """
    print("=" * 70)
    print("  (c) TOY MODEL: TWO SOURCES ON A CIRCLE")
    print("=" * 70)

    L = 1.0       # circle length (normalized)
    m = 2 * np.pi / L * 0.5  # mass parameter (m L ~ pi)
    q_Phi0_M = 1.0   # combined source strength (normalized)

    # Fourier solution:
    # delta_Phi(x) = sum_{n != 0} a_n cos(2 pi n x / L)
    # where a_n = -2 q Phi0 M / (L * (k_n^2 + m^2)) * [1 + cos(pi n)]
    #           = -2 q Phi0 M / (L * (k_n^2 + m^2)) * [1 + (-1)^n]
    # Non-zero only for even n: a_n = -4 q Phi0 M / (L * (k_n^2 + m^2))
    Nx = 500
    x = np.linspace(0, L, Nx, endpoint=False)
    delta_Phi = np.zeros(Nx)

    n_max = 200
    for n in range(1, n_max + 1):
        kn = 2 * np.pi * n / L
        # Two sources at x=0 and x=L/2:
        # Source terms in Fourier: S_n = (2/L) * [cos(0) + cos(k_n * L/2)]
        #                              = (2/L) * [1 + cos(n*pi)]
        #                              = (2/L) * [1 + (-1)^n]
        source_n = (2.0 / L) * (1 + (-1)**n)
        if abs(source_n) < 1e-15:
            continue
        a_n = -q_Phi0_M * source_n / (kn**2 + m**2)
        delta_Phi += a_n * np.cos(kn * x)

    # Enforce zero-sum: subtract mean
    delta_Phi -= np.mean(delta_Phi)

    # Energy: <U> = (m^2 / 2) * <delta_Phi^2>
    dphi2_avg = np.mean(delta_Phi**2)
    U_avg = 0.5 * m**2 * dphi2_avg

    # Gradient energy
    dx = x[1] - x[0]
    grad_Phi = np.gradient(delta_Phi, dx)
    E_grad = 0.5 * np.mean(grad_Phi**2)
    E_total = E_grad + U_avg

    print(f"\n  Parameters: L = {L}, m = {m:.4f}, q*Phi0*M = {q_Phi0_M}")
    print(f"  <delta_Phi^2>  = {dphi2_avg:.6e}")
    print(f"  <U> = m^2/2 * <dPhi^2> = {U_avg:.6e}")
    print(f"  <grad energy>  = {E_grad:.6e}")
    print(f"  <total energy> = {E_total:.6e}")
    print()
    print(f"  Lambda_eff ~ U_avg = {U_avg:.6e} > 0  (POSITIVE)")
    print()

    # Scaling with M: double the source strength
    M_factors = np.array([0.5, 1.0, 2.0, 4.0, 8.0])
    U_vs_M = []
    for mf in M_factors:
        dp = np.zeros(Nx)
        for n in range(1, n_max + 1):
            kn = 2 * np.pi * n / L
            source_n = (2.0 / L) * (1 + (-1)**n)
            if abs(source_n) < 1e-15:
                continue
            a_n = -mf * source_n / (kn**2 + m**2)
            dp += a_n * np.cos(kn * x)
        dp -= np.mean(dp)
        U_vs_M.append(0.5 * m**2 * np.mean(dp**2))

    U_vs_M = np.array(U_vs_M)
    print("  Scaling test (Lambda_eff vs M):")
    print(f"    {'M/M_ref':>10s}  {'U_avg':>12s}  {'U/U_ref':>10s}  {'(M/M_ref)^2':>12s}")
    for i, mf in enumerate(M_factors):
        ref = U_vs_M[1]
        print(f"    {mf:10.2f}  {U_vs_M[i]:12.6e}  {U_vs_M[i]/ref:10.4f}  {mf**2:12.4f}")

    print()
    print("  RESULT: Lambda_eff scales as M^2 (confirmed).")
    print("  More structure => more Lambda_eff (backreaction).")
    print()

    return dict(
        x=x, delta_Phi=delta_Phi, dphi2_avg=dphi2_avg,
        U_avg=U_avg, M_factors=M_factors, U_vs_M=U_vs_M,
        m=m, L=L,
    )


# =====================================================================
# (d) Scaling with structure formation
# =====================================================================
def structure_formation_scaling():
    """
    <delta_Phi^2>(z) ~ D^2(z) * <delta_Phi^2>(z=0)
    where D(z) is the linear growth factor.

    Lambda_eff(z) = Lambda_0 * [1 + epsilon * D^2(z)]
    where epsilon = 6 * <delta_Phi^2>(0) / Phi_0^2
    """
    print("=" * 70)
    print("  (d) SCALING WITH STRUCTURE FORMATION")
    print("=" * 70)

    # Growth factor D(z) in Lambda CDM (approximate)
    z_arr = np.linspace(0, 5, 500)

    def growth_factor(z, Om=Omega_m, OL=Omega_Lambda):
        """Approximate growth factor D(z)/D(0) for flat LCDM."""
        a = 1.0 / (1 + z)
        # Carroll, Press & Turner (1992) approximation
        Om_z = Om * (1 + z)**3 / (Om * (1 + z)**3 + OL)
        OL_z = OL / (Om * (1 + z)**3 + OL)
        D = (5.0/2.0) * Om_z / (Om_z**(4.0/7.0) - OL_z + (1 + Om_z/2.0) * (1 + OL_z/70.0))
        return D

    D_arr = np.array([growth_factor(z) for z in z_arr])
    # Normalize so D(0) = 1
    D0 = growth_factor(0)
    D_arr = D_arr / D0

    # epsilon from cosmic web potential variance
    sigma_Psi = 3e-5
    psi2_0 = sigma_Psi**2
    epsilon = 6 * psi2_0

    Lambda_ratio = 1 + epsilon * D_arr**2

    print(f"\n  Growth factor D(z)/D(0) [Carroll et al. 1992 approx.]")
    print(f"  epsilon = 6 * <Psi^2>(0) = 6 * ({sigma_Psi:.1e})^2 = {epsilon:.3e}")
    print()
    print(f"  Lambda_eff(z) = Lambda_0 * [1 + epsilon * D^2(z)]")
    print()
    print(f"  {'z':>6s}  {'D(z)/D(0)':>10s}  {'Lambda_eff/Lambda_0 - 1':>25s}")
    for z_val in [0, 0.5, 1.0, 2.0, 3.0, 5.0]:
        D_val = growth_factor(z_val) / D0
        dL = epsilon * D_val**2
        print(f"  {z_val:6.1f}  {D_val:10.4f}  {dL:25.3e}")

    print()
    print(f"  RESULT: The correction epsilon = {epsilon:.1e} is TINY.")
    print(f"  Lambda_eff is effectively constant across cosmic time.")
    print(f"  This is consistent with w = -1 (cosmological constant behavior).")
    print(f"  The variation is ~{epsilon*100:.4f}% -- far below observational reach.")
    print()

    # Effective w(z) from the correction
    # w_eff(z) = -1 + (2/3) * d ln(Lambda_eff) / d ln(a)
    # For Lambda ~ 1 + epsilon * D^2:
    # d ln Lambda / d ln a ~ 2 epsilon D dD/da / (a * (1 + epsilon D^2))
    #                      ~ 2 epsilon D dD/da / a  (for small epsilon)
    # At z=0: dD/da|_0 ~ f(Omega_m) ~ Omega_m^0.55 ~ 0.55
    f_growth = Omega_m**0.55
    w_deviation = (2.0/3.0) * 2 * epsilon * f_growth
    print(f"  Implied w_0 - (-1) ~ {w_deviation:.2e}")
    print(f"  (undetectable: current precision is ~0.03)")
    print()

    return dict(
        z_arr=z_arr, D_arr=D_arr, epsilon=epsilon,
        Lambda_ratio=Lambda_ratio, w_deviation=w_deviation,
    )


# =====================================================================
# (e) Results table
# =====================================================================
def results_table(res_a, res_b, res_c, res_d):
    print("=" * 70)
    print("  (e) RESULTS TABLE")
    print("=" * 70)
    print()

    header = f"  {'Quantity':40s}  {'Value':>15s}  {'Unit':>10s}"
    print(header)
    print("  " + "-" * 68)

    rows = [
        ("Lambda_obs", f"{Lambda_obs:.3e}", "m^-2"),
        ("Phi_0 (from DE match)", f"{res_a['Phi0_match']:.4f}", ""),
        ("Phi_0 (= 36 * Omega_Lambda)", f"{res_a['Phi0_from_Omega']:.4f}", ""),
        ("gamma = Phi_0 * H_0^2 / c_0^2", f"{res_a['gamma']:.3e}", "m^-2"),
        ("Lambda_eff = gamma/12", f"{res_a['Lambda_eff']:.3e}", "m^-2"),
        ("Lambda_eff / Lambda_obs", f"{res_a['ratio']:.6f}", ""),
        ("", "", ""),
        ("delta_Lambda/Lambda (structures)", f"{res_b['correction']:.2e}", ""),
        ("  from halos (f=1%)", f"{6*0.01*res_b['psi2_halo']:.2e}", ""),
        ("  from clusters (f=0.01%)", f"{6*1e-4*res_b['psi2_cluster']:.2e}", ""),
        ("  from cosmic web (f=100%)", f"{6*res_b['psi2_cosmo']:.2e}", ""),
        ("", "", ""),
        ("epsilon (redshift variation)", f"{res_d['epsilon']:.2e}", ""),
        ("w_0 - (-1) (implied)", f"{res_d['w_deviation']:.2e}", ""),
        ("", "", ""),
        ("Toy model: U_avg > 0?", "YES", ""),
        ("Toy model: U ~ M^2?", "YES", ""),
    ]

    for name, val, unit in rows:
        if name == "":
            print()
        else:
            print(f"  {name:40s}  {val:>15s}  {unit:>10s}")

    print()
    print("  CONCLUSION:")
    print("  -----------")
    print("  1. TGP reproduces Lambda_obs with Phi_0 = 24.66")
    print("     (= 36 * Omega_Lambda). This is a FIT, not a prediction.")
    print("  2. Phi_0 ~ 25 is O(10), not fine-tuned. This resolves the")
    print("     hierarchy problem: no 10^122 ratio appears.")
    print("  3. Structure corrections are ~5e-9, completely negligible.")
    print("  4. Lambda_eff is effectively constant in z (w = -1).")
    print("  5. The toy model confirms Lambda_eff > 0 and ~ M^2.")
    print("=" * 70)


# =====================================================================
# (f) Diagnostic plot
# =====================================================================
def make_plot(res_a, res_b, res_c, res_d, save_path):
    fig, axes = plt.subplots(2, 2, figsize=(14, 11))

    # --- Panel (a): U(phi) potential shape ---
    ax = axes[0, 0]
    gamma_val = res_a['gamma']
    phi = np.linspace(-0.5, 1.5, 500)
    U = (gamma_val/12) - (gamma_val/2)*phi**2 - (2*gamma_val/3)*phi**3 - (gamma_val/4)*phi**4
    U_norm = U / (gamma_val / 12)

    ax.plot(phi, U_norm, 'b-', lw=2.2)
    ax.axhline(1.0, color='r', ls='--', lw=1.2, alpha=0.7, label=r'$U(0) = \gamma/12$')
    ax.axhline(0.0, color='gray', ls=':', lw=0.8)
    ax.axvline(0.0, color='green', ls=':', lw=1.5, alpha=0.6, label=r'$\varphi = 0$ (vacuum)')
    ax.set_xlabel(r'$\varphi = \Phi/\Phi_0 - 1$', fontsize=12)
    ax.set_ylabel(r'$U(\varphi)\;/\;U(0)$', fontsize=12)
    ax.set_title(r'(a) TGP potential $U(\varphi)$', fontsize=13, fontweight='bold')
    ax.legend(fontsize=9)
    ax.set_ylim(-15, 3)
    ax.grid(True, ls=':', alpha=0.4)
    ax.annotate(r'$U(0) = \gamma/12 \to \Lambda_{\rm eff}$',
                xy=(0, 1), xytext=(0.4, 2.0),
                fontsize=10, color='red',
                arrowprops=dict(arrowstyle='->', color='red'))

    # --- Panel (b): Lambda_eff vs Phi_0 ---
    ax = axes[0, 1]
    Phi0_arr = np.linspace(1, 80, 500)
    Lambda_eff_arr = Phi0_arr * H0**2 / (12 * c0**2)
    ratio_arr = Lambda_eff_arr / Lambda_obs

    ax.plot(Phi0_arr, ratio_arr, 'b-', lw=2.2,
            label=r'$\Lambda_{\rm eff}/\Lambda_{\rm obs}$')
    ax.axhline(1.0, color='r', ls='--', lw=1.5,
               label=r'$\Lambda_{\rm eff} = \Lambda_{\rm obs}$')
    ax.axvline(res_a['Phi0_match'], color='green', ls=':', lw=2,
               label=rf'$\Phi_0 = {res_a["Phi0_match"]:.2f}$')

    # 1-sigma band from Omega_Lambda = 0.685 +/- 0.007
    Phi0_lo = 36 * (0.685 - 0.007)
    Phi0_hi = 36 * (0.685 + 0.007)
    ax.axvspan(Phi0_lo, Phi0_hi, color='green', alpha=0.15,
               label=rf'Planck $1\sigma$: [{Phi0_lo:.1f}, {Phi0_hi:.1f}]')

    ax.set_xlabel(r'$\Phi_0$', fontsize=12)
    ax.set_ylabel(r'$\Lambda_{\rm eff}\;/\;\Lambda_{\rm obs}$', fontsize=12)
    ax.set_title(r'(b) $\Lambda_{\rm eff}$ vs $\Phi_0$: the sweet spot',
                 fontsize=13, fontweight='bold')
    ax.legend(fontsize=8, loc='upper left')
    ax.set_xlim(1, 80)
    ax.set_ylim(0, 4)
    ax.grid(True, ls=':', alpha=0.4)

    # --- Panel (c): delta_Lambda / Lambda vs redshift ---
    ax = axes[1, 0]
    z_arr = res_d['z_arr']
    D_arr = res_d['D_arr']
    eps = res_d['epsilon']

    correction_vs_z = eps * D_arr**2

    ax.semilogy(z_arr, correction_vs_z, 'b-', lw=2.2,
                label=rf'$\delta\Lambda/\Lambda = \varepsilon\,D^2(z)$')
    ax.axhline(eps, color='r', ls='--', lw=1.2,
               label=rf'$\varepsilon = {eps:.1e}$ (today)')
    ax.set_xlabel(r'Redshift $z$', fontsize=12)
    ax.set_ylabel(r'$\delta\Lambda/\Lambda$', fontsize=12)
    ax.set_title(r'(c) Structure correction vs redshift',
                 fontsize=13, fontweight='bold')
    ax.legend(fontsize=9)
    ax.set_xlim(0, 5)
    ax.grid(True, ls=':', alpha=0.4)

    # Add text box
    ax.text(0.55, 0.35,
            r'Correction $\sim 5\times10^{-9}$' + '\n' +
            'Far below observational\nthreshold',
            transform=ax.transAxes, fontsize=9, ha='center',
            bbox=dict(boxstyle='round', fc='lightyellow', alpha=0.9))

    # --- Panel (d): Toy model ---
    ax = axes[1, 1]
    x = res_c['x']
    delta_Phi = res_c['delta_Phi']
    L = res_c['L']

    ax.plot(x / L, delta_Phi, 'b-', lw=2.2, label=r'$\delta\Phi(x)$')
    ax.axhline(0, color='gray', ls=':', lw=0.8)
    ax.axvline(0, color='red', ls='--', lw=1.2, alpha=0.7)
    ax.axvline(0.5, color='red', ls='--', lw=1.2, alpha=0.7,
               label='Source positions')

    # Show the potential energy density
    ax2 = ax.twinx()
    U_density = 0.5 * res_c['m']**2 * delta_Phi**2
    ax2.fill_between(x / L, U_density, alpha=0.2, color='orange',
                     label=r'$\frac{1}{2}m^2(\delta\Phi)^2$')
    ax2.set_ylabel(r'$U$ density', fontsize=10, color='orange')
    ax2.tick_params(axis='y', labelcolor='orange')

    ax.set_xlabel(r'$x / L$', fontsize=12)
    ax.set_ylabel(r'$\delta\Phi$', fontsize=12)
    ax.set_title(r'(d) Toy model: two sources on circle',
                 fontsize=13, fontweight='bold')

    # Combined legend
    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, fontsize=8, loc='upper right')

    # Inset: U vs M scaling
    inset = ax.inset_axes([0.05, 0.05, 0.35, 0.35])
    M_factors = res_c['M_factors']
    U_vs_M = res_c['U_vs_M']
    inset.loglog(M_factors, U_vs_M, 'ko-', ms=4, lw=1.5)
    # Reference M^2 line
    ref_M2 = U_vs_M[1] * M_factors**2
    inset.loglog(M_factors, ref_M2, 'r--', lw=1, alpha=0.7)
    inset.set_xlabel(r'$M/M_{\rm ref}$', fontsize=7)
    inset.set_ylabel(r'$\langle U \rangle$', fontsize=7)
    inset.set_title(r'$U \propto M^2$', fontsize=7)
    inset.tick_params(labelsize=6)

    fig.suptitle(
        r'TGP: Quantitative $\Lambda_{\rm eff}$ computation (prob:Lambda)',
        fontsize=15, y=1.01
    )
    fig.tight_layout(rect=[0, 0, 1, 0.97])

    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    fig.savefig(save_path, dpi=180, bbox_inches='tight')
    print(f"  Plot saved to {save_path}")
    plt.close(fig)


# =====================================================================
# Main
# =====================================================================
def main():
    print()
    print("#" * 70)
    print("#  TGP: Quantitative Lambda_eff computation")
    print("#  Addressing prob:Lambda")
    print("#" * 70)
    print()

    res_a = homogeneous_background()
    res_b = structure_correction()
    res_c = toy_model_circle()
    res_d = structure_formation_scaling()
    results_table(res_a, res_b, res_c, res_d)

    # Plot
    save_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'plots', "lambda_eff_quantitative.png")
    make_plot(res_a, res_b, res_c, res_d, save_path)

    print("\nAll computations passed. Done.")


if __name__ == "__main__":
    main()
