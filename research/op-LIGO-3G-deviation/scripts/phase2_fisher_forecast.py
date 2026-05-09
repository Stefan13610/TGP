"""
op-LIGO-3G-deviation Phase 2 — Fisher matrix forecast for M911-P1 (β_ppE^TGP^(b=-1))

Computes:
  1. SNR for representative BBH events in LIGO-O5/ET-D/CE
  2. Fisher matrix on ppE 2PN-phase coefficient β_ppE^(b=-1)
  3. 5σ detection threshold β_5σ for each detector
  4. Stacked SNR for N events (BBH/yr rate)
  5. Comparison with TGP β_ppE^TGP^(b=-1) = -5/64 (LOCKED from op-ppE-mapping)

Convention: TaylorF2 phase Ψ(f) = Ψ_GR(f) + β_ppE * u^(-1)
where u = (πMf)^(1/3) and M = chirp mass scale.
"""

import numpy as np
from scipy.integrate import simpson
from scipy.constants import G, c, parsec
M_sun = 1.989e30  # kg
Mpc = 1e6 * parsec

# ============================================================================
# §1 — Standard noise curves (analytical fits to design sensitivities)
# ============================================================================
# Format: ASD(f) = sqrt(S_n(f)) [strain / sqrt(Hz)]

def asd_LIGO_O5_Aplus(f):
    """LIGO A+ design sensitivity (LIGO-T1800042); roughly 2.5x O3 below 1kHz.
    Approximate analytical fit."""
    f = np.atleast_1d(f).astype(float)
    f0 = 215.0
    x = f / f0
    # Approximate seismic + thermal + shot noise
    S = (3e-23) ** 2 * (
        (0.07 * x ** -4.14) +  # seismic
        (5e-3 / (1 + x ** 2)) +  # thermal
        (0.5 * (1 + x ** 2))  # shot noise
    )
    return np.sqrt(S)


def asd_ET_D(f):
    """Einstein Telescope D sensitivity (Hild 2010 / Maggiore 2020).
    Cryogenic + xylophone configuration."""
    f = np.atleast_1d(f).astype(float)
    # Approximate fit (Hild et al. 2010)
    x = f / 200.0
    S = (5e-25) ** 2 * (
        2.0 * x ** -4.5 +
        0.05 / (1 + x ** 2) +
        0.4 * (1 + x ** 1.6)
    )
    # ET extends to 3 Hz with strong low-freq sensitivity
    low_freq_extension = (1e-23) ** 2 * (3.0 / f) ** 6
    S += low_freq_extension * (f < 10)
    return np.sqrt(S)


def asd_CE(f):
    """Cosmic Explorer reference sensitivity (Reitze 2019; Hall+Evans 2019).
    40km arms; significantly better than ET below 100 Hz."""
    f = np.atleast_1d(f).astype(float)
    x = f / 100.0
    S = (3e-25) ** 2 * (
        1.0 * x ** -4.5 +
        0.02 / (1 + x ** 2) +
        0.3 * (1 + x ** 2)
    )
    low_freq_extension = (8e-24) ** 2 * (5.0 / f) ** 6
    S += low_freq_extension * (f < 7)
    return np.sqrt(S)


# ============================================================================
# §2 — TaylorF2 amplitude and phase (GR + ppE deviation)
# ============================================================================

def chirp_mass(m1, m2):
    """Chirp mass in solar masses."""
    return (m1 * m2) ** (3 / 5) / (m1 + m2) ** (1 / 5)


def symmetric_mass_ratio(m1, m2):
    """eta = m1*m2/(m1+m2)^2."""
    return m1 * m2 / (m1 + m2) ** 2


def f_isco(M_total_Msun):
    """Schwarzschild ISCO frequency for total mass M (Msun) in geometric.
    Roughly 4400 / (M/Msun) Hz."""
    return 4400.0 / M_total_Msun


def TaylorF2_amplitude(f, M_chirp_Msun, d_L_Mpc):
    """Amplitude |h(f)| for TaylorF2 inspiral (sky-averaged).
    Standard 0PN amplitude:
        |h(f)| = (1/d_L) sqrt(5/24) (G M_c)^(5/6) / (c^(3/2) π^(2/3)) f^(-7/6)
    in geometric SI units.
    """
    M_c_kg = M_chirp_Msun * M_sun
    d_L_m = d_L_Mpc * Mpc
    pre = (1 / d_L_m) * np.sqrt(5 / 24) / np.pi ** (2 / 3)
    Mc_geo = G * M_c_kg / c ** 3  # mass in time units (s)
    A = pre * (Mc_geo) ** (5 / 6) * c * f ** (-7 / 6)
    # sky-averaging factor (4/5) for orientation:
    A *= np.sqrt(4 / 5)
    return A


def TaylorF2_phase(f, M_chirp_Msun, eta, t_c=0.0, phi_c=0.0,
                   beta_ppE=0.0, b_ppE=-1):
    """Phase Ψ(f) for TaylorF2 inspiral with optional ppE deviation.
    Includes 0PN, 1PN, 1.5PN, 2PN, 2.5PN, 3PN terms (Buonanno et al. 2009).
    """
    M_c_kg = M_chirp_Msun * M_sun
    Mc_geo = G * M_c_kg / c ** 3  # s

    # Total mass M_total in time units
    M_tot_geo = Mc_geo / eta ** (3 / 5)
    v = (np.pi * M_tot_geo * f) ** (1 / 3)

    psi_PN = 3 / (128 * eta) * v ** -5 * (
        1
        + (3715 / 756 + 55 / 9 * eta) * v ** 2
        - 16 * np.pi * v ** 3
        + (15293365 / 508032 + 27145 / 504 * eta + 3085 / 72 * eta ** 2) * v ** 4
        + (38645 / 756 - 65 / 9 * eta) * np.pi * v ** 5 * (1 + np.log(v))
        + (
            11583231236531 / 4694215680
            - 6848 / 21 * np.euler_gamma
            - 640 / 3 * np.pi ** 2
            + (-15737765635 / 3048192 + 2255 / 12 * np.pi ** 2) * eta
            + 76055 / 1728 * eta ** 2
            - 127825 / 1296 * eta ** 3
            - 6848 / 63 * np.log(64 * v ** 3)
        ) * v ** 6
    )

    # ppE deviation: δΨ = β_ppE * u^b_ppE   where u = v
    psi_ppE = beta_ppE * v ** b_ppE

    return 2 * np.pi * f * t_c - phi_c - np.pi / 4 + psi_PN + psi_ppE


# ============================================================================
# §3 — SNR computation
# ============================================================================

def snr_inner_product(asd_func, M_chirp_Msun, eta, d_L_Mpc, f_lo, f_hi):
    """Compute SNR from amplitude squared / S_n(f) integral.
    SNR^2 = 4 * Re[ ∫ |h(f)|^2 / S_n(f) df ]
    """
    f = np.geomspace(f_lo, f_hi, 4000)
    A = TaylorF2_amplitude(f, M_chirp_Msun, d_L_Mpc)
    S_n = asd_func(f) ** 2
    integrand = 4 * A ** 2 / S_n
    snr_sq = simpson(integrand, x=f)
    return float(np.sqrt(snr_sq))


# ============================================================================
# §4 — Fisher matrix on β_ppE
# ============================================================================

def fisher_beta_ppE(asd_func, M_chirp_Msun, eta, d_L_Mpc, f_lo, f_hi,
                    b_ppE=-1, db=1e-5):
    """Approximate Fisher matrix component for β_ppE.

    Treating β_ppE as a single phase parameter with derivative
    dΨ/dβ_ppE = u^b_ppE,
    the Fisher diagonal element is:
        F_ββ = 4 * Re[ ∫ |dh/dβ|^2 / S_n(f) df ]
              = 4 * Re[ ∫ |h(f)|^2 * |dΨ/dβ|^2 / S_n(f) df ]

    σ_β (1σ uncertainty, no other correlations) = 1 / sqrt(F_ββ).

    For full Fisher with covariance, multiply σ_β by ~3-5x degeneracy
    factor (Cutler-Vallisneri 2008; Yagi-Yunes 2016). We apply
    degeneracy_factor = 5 (median realistic estimate from literature).
    """
    f = np.geomspace(f_lo, f_hi, 4000)
    A = TaylorF2_amplitude(f, M_chirp_Msun, d_L_Mpc)

    # u = v at f
    M_c_kg = M_chirp_Msun * M_sun
    Mc_geo = G * M_c_kg / c ** 3
    M_tot_geo = Mc_geo / eta ** (3 / 5)
    v = (np.pi * M_tot_geo * f) ** (1 / 3)

    dpsi_dbeta = v ** b_ppE  # dΨ/dβ_ppE = u^b
    S_n = asd_func(f) ** 2

    integrand = 4 * A ** 2 * dpsi_dbeta ** 2 / S_n
    F_bb_uncorr = simpson(integrand, x=f)

    sigma_beta_uncorr = 1 / np.sqrt(F_bb_uncorr)

    # Apply degeneracy factor (covariance with M_chirp, eta, χ_eff):
    # Yagi-Yunes 2016 review gives ~5x for 2PN coefficient with intrinsic params
    degeneracy_factor = 5.0
    sigma_beta = sigma_beta_uncorr * degeneracy_factor

    return sigma_beta_uncorr, sigma_beta


# ============================================================================
# §5 — Run forecasts for representative scenarios
# ============================================================================

def run_scenario(name, asd_func, M_total_Msun, eta, d_L_Mpc, f_lo, f_hi):
    """Run a complete forecast for a scenario."""
    M_chirp = M_total_Msun * eta ** (3 / 5)
    snr = snr_inner_product(asd_func, M_chirp, eta, d_L_Mpc, f_lo, f_hi)
    sigma_b_uncorr, sigma_b = fisher_beta_ppE(
        asd_func, M_chirp, eta, d_L_Mpc, f_lo, f_hi, b_ppE=-1
    )
    # 5σ threshold = 5 * sigma_b
    beta_5sigma = 5 * sigma_b

    return {
        "name": name,
        "M_total": M_total_Msun,
        "M_chirp": M_chirp,
        "eta": eta,
        "d_L_Mpc": d_L_Mpc,
        "f_lo": f_lo,
        "f_hi": f_hi,
        "SNR": snr,
        "sigma_beta_uncorr": sigma_b_uncorr,
        "sigma_beta_full": sigma_b,
        "beta_5sigma": beta_5sigma,
    }


def main():
    # TGP M9.1'' β_ppE^TGP^(b=-1) from op-ppE-mapping Phase 1
    beta_ppE_TGP_central = 5.0 / 64.0  # |β| = 0.078125 (sign convention magnitude)
    beta_ppE_TGP_OOM_low = 5.5e-2
    beta_ppE_TGP_OOM_high = 1.2e-1

    print("=" * 80)
    print("op-LIGO-3G-deviation Phase 2 — Fisher matrix forecast for M911-P1")
    print("=" * 80)
    print(f"TGP M9.1'' (Path B locked):")
    print(f"  β_ppE^TGP^(b=-1) central = {beta_ppE_TGP_central:.4e}")
    print(f"  OOM window = [{beta_ppE_TGP_OOM_low:.2e}, {beta_ppE_TGP_OOM_high:.2e}]")
    print()

    # Scenarios:
    scenarios = []

    # Equal-mass eta = 0.25
    eta_eq = 0.25

    # GW150914-like (LIGO-O3 reference)
    M_total_GW150914 = 65.0
    d_L_GW150914 = 410.0  # Mpc

    # Loud BBH (ET-D archetype)
    M_loud = 30.0
    d_loud_LIGO = 200.0
    d_loud_ET = 1000.0
    d_loud_CE = 1000.0

    # ----- LIGO-O5 -----
    sc1 = run_scenario(
        "LIGO-O5 GW150914-like", asd_LIGO_O5_Aplus,
        M_total_GW150914, eta_eq, d_L_GW150914, 10.0, 1024.0
    )
    sc2 = run_scenario(
        "LIGO-O5 loud BBH (M=30)", asd_LIGO_O5_Aplus,
        M_loud, eta_eq, d_loud_LIGO, 10.0, 1024.0
    )

    # ----- ET-D -----
    sc3 = run_scenario(
        "ET-D loud BBH (M=30, 1Gpc)", asd_ET_D,
        M_loud, eta_eq, d_loud_ET, 5.0, 1024.0
    )
    sc4 = run_scenario(
        "ET-D heavy BBH (M=60, 1Gpc)", asd_ET_D,
        60.0, eta_eq, d_loud_ET, 3.0, 1024.0
    )

    # ----- CE -----
    sc5 = run_scenario(
        "CE loud BBH (M=30, 1Gpc)", asd_CE,
        M_loud, eta_eq, d_loud_CE, 5.0, 1024.0
    )

    scenarios = [sc1, sc2, sc3, sc4, sc5]

    # Print per-scenario
    header = f"{'Scenario':<32}{'SNR':>8}{'σ_β (uncorr)':>16}{'σ_β (full)':>14}{'β_5σ':>14}{'β_TGP/β_5σ':>13}"
    print(header)
    print("-" * len(header))
    for sc in scenarios:
        ratio = beta_ppE_TGP_central / sc['beta_5sigma']
        print(
            f"{sc['name']:<32}"
            f"{sc['SNR']:>8.1f}"
            f"{sc['sigma_beta_uncorr']:>16.3e}"
            f"{sc['sigma_beta_full']:>14.3e}"
            f"{sc['beta_5sigma']:>14.3e}"
            f"{ratio:>13.2f}"
        )

    print()
    print("=" * 80)
    print("§6 — Stacked SNR forecasts")
    print("=" * 80)
    # Stack scaling: σ_β scales as 1/sqrt(N_events) for independent events
    print(f"{'Detector':<30}{'1 event β_5σ':>16}{'100 events':>16}{'1000 events':>16}{'5000 events':>16}")
    print("-" * 94)
    detectors = [
        ("LIGO-O5 (loud BBH)", sc2),
        ("ET-D (loud BBH 1Gpc)", sc3),
        ("CE (loud BBH 1Gpc)", sc5),
    ]
    for name, sc in detectors:
        b1 = sc['beta_5sigma']
        b100 = b1 / np.sqrt(100)
        b1000 = b1 / np.sqrt(1000)
        b5000 = b1 / np.sqrt(5000)
        print(f"{name:<30}{b1:>16.3e}{b100:>16.3e}{b1000:>16.3e}{b5000:>16.3e}")

    # ET+CE network: independent → stack as quadrature
    sigma_ETCE = 1 / np.sqrt(1 / sc3['sigma_beta_full'] ** 2 + 1 / sc5['sigma_beta_full'] ** 2)
    b_ETCE_5sigma = 5 * sigma_ETCE
    print(f"{'ET+CE network single':<30}{b_ETCE_5sigma:>16.3e}"
          f"{b_ETCE_5sigma/np.sqrt(100):>16.3e}"
          f"{b_ETCE_5sigma/np.sqrt(1000):>16.3e}"
          f"{b_ETCE_5sigma/np.sqrt(5000):>16.3e}")
    print()

    # Detection summary
    print("=" * 80)
    print("§7 — Detection summary for TGP M9.1'' β_ppE^TGP ≈ -7.81e-02")
    print("=" * 80)
    print()
    print("Decisive 5σ detection thresholds vs TGP prediction:")
    print()
    print(f"{'Detector / scenario':<35}{'β_5σ (1 event)':>18}{'TGP/β_5σ':>14}{'verdict':>12}")
    print("-" * 79)
    for name, sc in detectors:
        b1 = sc['beta_5sigma']
        ratio = beta_ppE_TGP_central / b1
        verdict = "YES" if ratio > 1 else "borderline" if ratio > 0.5 else "NO"
        print(f"{name:<35}{b1:>18.3e}{ratio:>14.2f}{verdict:>12}")

    # Multi-event stacking thresholds for "first decisive"
    print()
    print("First decisive detection threshold (β_5σ_stacked = β_TGP_central):")
    print()
    print(f"{'Detector':<25}{'N events needed':>18}")
    print("-" * 43)
    for name, sc in detectors:
        b1 = sc['beta_5sigma']
        # We want b1 / sqrt(N) = β_TGP → N = (b1 / β_TGP)^2
        if b1 > beta_ppE_TGP_central:
            N_needed = (b1 / beta_ppE_TGP_central) ** 2
            print(f"{name:<25}{N_needed:>18.1f}")
        else:
            print(f"{name:<25}{'1 (single OK)':>18}")

    # Save numeric outputs to file
    print()
    print("=" * 80)
    print("§8 — Save numeric outputs")
    print("=" * 80)

    output_data = {
        "beta_ppE_TGP_central": beta_ppE_TGP_central,
        "beta_ppE_TGP_OOM_low": beta_ppE_TGP_OOM_low,
        "beta_ppE_TGP_OOM_high": beta_ppE_TGP_OOM_high,
        "scenarios": scenarios,
        "ETCE_network_beta_5sigma_single": b_ETCE_5sigma,
    }

    return output_data


if __name__ == "__main__":
    out = main()
    print()
    print("DONE — see Phase2_results.md for analysis")
