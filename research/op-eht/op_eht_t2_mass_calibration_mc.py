"""
OP-EHT T2 — Mass calibration Monte Carlo.

Question: czy effective EHT envelope (combining mass, distance, ring-modeling
uncertainties) absorbuje +14.6% TGP deviation?

Method:
1. For each BH (M87*, Sgr A*):
   - Sample (M, D) from published priors
   - Sample ring modeling systematic from published EHT systematics
   - Predict shadow diameter under (a) GR Schwarzschild, (b) TGP M9.1''
   - Compute deviation (TGP_predicted - EHT_observed) / sigma_total
2. Run N=10000 samples per source.
3. Report 95% CL of (TGP - EHT) deviation in units of sigma_total.

Priors:
- M87* mass: Gebhardt et al 2011 (stellar dynamics) M = (6.6 +/- 0.4)e9 M_sun
            Walsh et al 2013 (gas dynamics) M = (3.5 +/- 0.7)e9 M_sun
            EHT 2019 inferred M = (6.5 +/- 0.7)e9 M_sun (using GR shadow!)
            We use M = (6.5 +/- 0.7)e9 M_sun (most direct)
- M87* distance: D = 16.8 +/- 0.8 Mpc (TRGB; Bird et al 2010)
- M87* ring diameter (EHT 2019): theta = 42.0 +/- 3.0 µas (10% systematic)

- Sgr A* mass: Gravity 2019: M = (4.297 +/- 0.013)e6 M_sun (very tight)
              Stellar orbits: M = (4.0 +/- 0.6)e6 M_sun
              Combined: ~(4.3 +/- 0.1)e6 M_sun
- Sgr A* distance: GRAVITY 2019: D = 8.277 +/- 0.033 kpc
- Sgr A* ring diameter (EHT 2022): theta = 51.8 +/- 2.3 µas (4.4% precision)

Prediction:
- Shadow diameter (uas) = 2 * b_crit * (G M / c^2) / D * (180/pi) * (3600 * 1e6)
                       = 2 * b_crit * r_g_uas
- TGP/GR ratio: 5.952/5.196 = 1.1455 (uniform across sources)

PASS criteria:
- T2.1: M87* posterior consistent with TGP at 95% CL.
- T2.2: Sgr A* posterior consistent with TGP at 95% CL.
- T2.3: Combined evidence (M87* + Sgr A*): TGP within 95% CL.
"""
import numpy as np

# Universal TGP prediction (from EHT-quick / T1 robust)
B_CRIT_GR = 3.0 * np.sqrt(3.0)       # GR Schwarzschild
B_CRIT_TGP = 5.952474                  # TGP M9.1'' from T1 n=15
RATIO_TGP_GR = B_CRIT_TGP / B_CRIT_GR  # ~1.1456

# Constants
G = 6.6743e-11        # m^3/kg/s^2
c = 2.998e8           # m/s
M_sun = 1.989e30      # kg
Mpc = 3.0857e22       # m
kpc = 3.0857e19       # m
uas_per_rad = (180.0/np.pi) * 3600.0 * 1e6  # microarcsec per radian


def shadow_uas(M_kg, D_m, b_crit_units):
    """Shadow diameter in microarcseconds.
    theta = 2 * b_crit * r_g / D, where r_g = G M / c^2.
    """
    r_g = G * M_kg / c**2
    theta_rad = 2.0 * b_crit_units * r_g / D_m
    return theta_rad * uas_per_rad


def main():
    print("=" * 70)
    print(" OP-EHT T2 — Mass calibration Monte Carlo")
    print(" Question: does effective EHT envelope absorb +14.6% TGP deviation?")
    print("=" * 70)

    print(f"\n[Constants]")
    print(f"  b_crit_GR    = {B_CRIT_GR:.6f} r_g")
    print(f"  b_crit_TGP   = {B_CRIT_TGP:.6f} r_g")
    print(f"  TGP/GR       = {RATIO_TGP_GR:.4f} ({(RATIO_TGP_GR-1)*100:+.2f}%)")

    rng = np.random.default_rng(seed=20260425)
    N = 10000

    pass_count = 0
    n_total = 3

    # --------------------------------------------------------------
    # M87*
    # --------------------------------------------------------------
    print(f"\n[T2.1] M87* posterior (N={N}):")
    M87_M_mean = 6.5e9 * M_sun
    M87_M_sigma = 0.7e9 * M_sun
    M87_D_mean = 16.8 * Mpc
    M87_D_sigma = 0.8 * Mpc
    M87_theta_obs = 42.0  # µas
    M87_theta_sigma = 3.0  # µas

    M87_M_samples = rng.normal(M87_M_mean, M87_M_sigma, N)
    M87_D_samples = rng.normal(M87_D_mean, M87_D_sigma, N)

    # GR predictions for sample
    M87_theta_GR = np.array([
        shadow_uas(M, D, B_CRIT_GR)
        for M, D in zip(M87_M_samples, M87_D_samples)
    ])
    # TGP predictions
    M87_theta_TGP = np.array([
        shadow_uas(M, D, B_CRIT_TGP)
        for M, D in zip(M87_M_samples, M87_D_samples)
    ])

    # Add ring modeling systematic to observation (theta_obs sample around 42 ± 3)
    M87_theta_obs_samples = rng.normal(M87_theta_obs, M87_theta_sigma, N)

    # Residuals: TGP_pred - obs (in sigma_obs units)
    M87_residuals_TGP = (M87_theta_TGP - M87_theta_obs_samples) / M87_theta_sigma
    M87_residuals_GR = (M87_theta_GR - M87_theta_obs_samples) / M87_theta_sigma

    print(f"  GR mean predicted theta:  {np.mean(M87_theta_GR):.2f} +/- {np.std(M87_theta_GR):.2f} uas")
    print(f"  TGP mean predicted theta: {np.mean(M87_theta_TGP):.2f} +/- {np.std(M87_theta_TGP):.2f} uas")
    print(f"  EHT 2019 observed:        {M87_theta_obs:.2f} +/- {M87_theta_sigma:.2f} uas")

    # 95% CL: |residual| < 1.96
    M87_TGP_in_95 = np.mean(np.abs(M87_residuals_TGP) < 1.96)
    M87_GR_in_95 = np.mean(np.abs(M87_residuals_GR) < 1.96)
    M87_TGP_pvalue = M87_TGP_in_95
    print(f"  Fraction TGP within 95% CL of EHT obs: {M87_TGP_in_95*100:.1f}%")
    print(f"  Fraction GR  within 95% CL of EHT obs: {M87_GR_in_95*100:.1f}%")

    # PASS if TGP at least 50% within 95% CL (i.e., consistent with data)
    t2_1 = M87_TGP_in_95 > 0.5
    pass_count += int(t2_1)
    print(f"  T2.1 M87* TGP consistent at 95% CL: {'PASS' if t2_1 else 'FAIL'}")

    # --------------------------------------------------------------
    # Sgr A*
    # --------------------------------------------------------------
    print(f"\n[T2.2] Sgr A* posterior (N={N}):")
    SgrA_M_mean = 4.297e6 * M_sun
    SgrA_M_sigma = 0.013e6 * M_sun  # GRAVITY 2019 tight prior
    SgrA_D_mean = 8.277 * kpc
    SgrA_D_sigma = 0.033 * kpc
    SgrA_theta_obs = 51.8  # µas
    SgrA_theta_sigma = 2.3  # µas (4.4% precision)

    SgrA_M_samples = rng.normal(SgrA_M_mean, SgrA_M_sigma, N)
    SgrA_D_samples = rng.normal(SgrA_D_mean, SgrA_D_sigma, N)

    SgrA_theta_GR = np.array([
        shadow_uas(M, D, B_CRIT_GR)
        for M, D in zip(SgrA_M_samples, SgrA_D_samples)
    ])
    SgrA_theta_TGP = np.array([
        shadow_uas(M, D, B_CRIT_TGP)
        for M, D in zip(SgrA_M_samples, SgrA_D_samples)
    ])

    SgrA_theta_obs_samples = rng.normal(SgrA_theta_obs, SgrA_theta_sigma, N)
    SgrA_residuals_TGP = (SgrA_theta_TGP - SgrA_theta_obs_samples) / SgrA_theta_sigma
    SgrA_residuals_GR = (SgrA_theta_GR - SgrA_theta_obs_samples) / SgrA_theta_sigma

    print(f"  GR mean predicted theta:  {np.mean(SgrA_theta_GR):.2f} +/- {np.std(SgrA_theta_GR):.2f} uas")
    print(f"  TGP mean predicted theta: {np.mean(SgrA_theta_TGP):.2f} +/- {np.std(SgrA_theta_TGP):.2f} uas")
    print(f"  EHT 2022 observed:        {SgrA_theta_obs:.2f} +/- {SgrA_theta_sigma:.2f} uas")

    SgrA_TGP_in_95 = np.mean(np.abs(SgrA_residuals_TGP) < 1.96)
    SgrA_GR_in_95 = np.mean(np.abs(SgrA_residuals_GR) < 1.96)
    print(f"  Fraction TGP within 95% CL of EHT obs: {SgrA_TGP_in_95*100:.1f}%")
    print(f"  Fraction GR  within 95% CL of EHT obs: {SgrA_GR_in_95*100:.1f}%")

    # Mean residual (positive = TGP overshoots)
    print(f"  Mean TGP residual (sigma units): {np.mean(SgrA_residuals_TGP):+.2f}")
    print(f"  Mean GR  residual (sigma units): {np.mean(SgrA_residuals_GR):+.2f}")

    t2_2 = SgrA_TGP_in_95 > 0.5
    pass_count += int(t2_2)
    print(f"  T2.2 Sgr A* TGP consistent at 95% CL: {'PASS' if t2_2 else 'FAIL'}")

    # --------------------------------------------------------------
    # Combined evidence
    # --------------------------------------------------------------
    print(f"\n[T2.3] Combined M87* + Sgr A* (joint chi^2):")
    chi2_TGP = M87_residuals_TGP**2 + SgrA_residuals_TGP**2
    chi2_GR = M87_residuals_GR**2 + SgrA_residuals_GR**2
    # 95% for chi2 with 2 dof: 5.991
    chi2_95_dof2 = 5.991
    TGP_in_joint_95 = np.mean(chi2_TGP < chi2_95_dof2)
    GR_in_joint_95 = np.mean(chi2_GR < chi2_95_dof2)
    print(f"  chi2 95% CL (2 dof) = {chi2_95_dof2:.3f}")
    print(f"  Mean chi2_TGP: {np.mean(chi2_TGP):.2f}")
    print(f"  Mean chi2_GR:  {np.mean(chi2_GR):.2f}")
    print(f"  Fraction TGP in joint 95% CL: {TGP_in_joint_95*100:.1f}%")
    print(f"  Fraction GR  in joint 95% CL: {GR_in_joint_95*100:.1f}%")

    t2_3 = TGP_in_joint_95 > 0.5
    pass_count += int(t2_3)
    print(f"  T2.3 Combined TGP at joint 95% CL: {'PASS' if t2_3 else 'FAIL'}")

    # --------------------------------------------------------------
    # Summary
    # --------------------------------------------------------------
    print(f"\n{'=' * 70}")
    print(f" SUMMARY: {pass_count}/{n_total} PASS")
    print(f" OP-EHT T2 verdict: ", end="")

    if pass_count == n_total:
        print("OP-EHT T2 closes POSITIVE — EHT envelope absorbs TGP deviation.")
    elif t2_1 and not t2_2:
        print("M87* OK, Sgr A* TENSION — TGP outside Sgr A* 95% CL.")
        print(" ⇒ Sgr A* is the discriminator. Defers to T3-T4.")
    elif not t2_1 and not t2_2:
        print("BOTH SOURCES TENSION — TGP outside 95% CL for M87* and Sgr A*.")
        print(" ⇒ T3-T4 must reduce deviation, otherwise OP-EHT NEGATIVE.")
    else:
        print(f"PARTIAL ({pass_count}/{n_total}) — see T3-T4 for resolution path.")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
