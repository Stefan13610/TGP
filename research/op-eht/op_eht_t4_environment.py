"""
OP-EHT T4 — Non-vacuum environment correction.

Question: czy accretion disk + plasma scattering moga przesunac apparent
shadow diameter na tyle ze +14.56% TGP deviation jest absorbowana przez
astrophysical environment?

Method:
1. Literature-based estimate of delta_theta / theta from:
   a) emission asymmetry (Doppler beaming, photon trajectory bending)
   b) synchrotron self-absorption (SSA)
   c) plasma scattering (interstellar + sublim. radius)
   d) emission ring vs photon ring offset (geometric thickness)

2. M87* (high accretion, M_dot ~ 1e-4 M_sun/yr):
   - SANE/MAD GRMHD simulations: Psaltis 2020 finds delta_theta/theta ~ 2-4%
   - Younsi 2023: emission ring radius differs from photon ring by ~5%
   - Plasma scattering: negligible (M87 is jetted, low ISM column)

3. Sgr A* (low accretion, M_dot ~ 1e-8 M_sun/yr):
   - Plasma scattering through galactic center: substantial (theta_scatt ~ 22 µas
     at 230 GHz; deconvolved by EHT 2022)
   - GRMHD: emission asymmetry ~ 1-2%
   - Disk/jet morphology uncertainty: ~3% (variability dominated)

Conservative envelope:
- M87*: delta_theta_env / theta = +/- 5% (combined emission+geometry)
- Sgr A*: delta_theta_env / theta = +/- 4% (after scattering deconv.)

PASS criteria:
- T4.1: M87* environment can absorb residual +8.3% deviation? Need >= 8%.
- T4.2: Sgr A* environment can absorb residual +19.7% deviation? Need >= 19.7%.
- T4.3: Combined OP-EHT NEGATIVE (Sgr A* tension persists) or POSITIVE?

CRITICAL: T4 tests whether environment is sufficient PHYSICAL mechanism;
if not (likely Sgr A*), then deferral to ngEHT (T5) becomes mandatory.
"""
import numpy as np

# Universal TGP prediction
B_CRIT_GR = 3.0 * np.sqrt(3.0)
B_CRIT_TGP = 5.952474
RATIO_TGP_GR = B_CRIT_TGP / B_CRIT_GR  # ~1.1456


def main():
    print("=" * 70)
    print(" OP-EHT T4 — Non-vacuum environment correction")
    print(" Question: does accretion disk + plasma absorb +14.56% deviation?")
    print("=" * 70)

    print(f"\n[Constants]")
    print(f"  TGP/GR ratio in b_crit:        {RATIO_TGP_GR:.4f}")
    print(f"  Universal +deviation:          +{(RATIO_TGP_GR-1)*100:.2f}%")

    pass_count = 0
    n_total = 3

    # --------------------------------------------------------------
    # T4.1 M87* environment budget
    # --------------------------------------------------------------
    print(f"\n[T4.1] M87* environment budget:")
    print(f"  Source: Psaltis 2020, Younsi 2023, EHT 2019 paper IV/V")
    print(f"  Mass accretion rate: M_dot ~ 1e-4 M_sun/yr (high)")
    print(f"  GRMHD emission ring vs photon ring: delta = +/- 5%")
    print(f"  SSA correction at 230 GHz: ~ +/- 1% (transparent disk)")
    print(f"  Plasma scattering: negligible (M87 ISM column low)")
    print(f"  Net: delta_theta / theta = +/- 5% (1-sigma)")
    M87_env_sigma_pct = 5.0  # percent
    M87_residual_after_T2 = 8.3  # +8.3% (M87* TGP already 1.5σ tension)
    print(f"  M87* residual after T2 (best M, D priors): +{M87_residual_after_T2:.1f}%")

    # Can environment absorb residual? Yes if env >= residual at 1.96σ
    # i.e. 5% * 1.96 = 9.8% reach. That's > 8.3%.
    M87_env_reach_2sigma = M87_env_sigma_pct * 1.96
    t4_1 = M87_env_reach_2sigma >= M87_residual_after_T2
    pass_count += int(t4_1)
    print(f"  Environment 95%-CL reach (2-sigma): +/- {M87_env_reach_2sigma:.1f}%")
    print(f"  T4.1 M87* env can absorb residual: {'PASS' if t4_1 else 'FAIL'}")

    # --------------------------------------------------------------
    # T4.2 Sgr A* environment budget
    # --------------------------------------------------------------
    print(f"\n[T4.2] Sgr A* environment budget:")
    print(f"  Source: EHT 2022 paper II/III/V, Bower 2014 (scattering)")
    print(f"  Mass accretion rate: M_dot ~ 1e-8 M_sun/yr (very low)")
    print(f"  Galactic-center plasma scattering: theta_scatt ~ 22 uas at 230 GHz")
    print(f"    (deconvolved by EHT 2022; residual systematic ~ 1-2%)")
    print(f"  GRMHD emission ring vs photon ring: delta = +/- 2-3%")
    print(f"  Variability + flaring: +/- 2% on ring diameter")
    print(f"  Net: delta_theta / theta = +/- 4% (1-sigma; conservative)")
    SgrA_env_sigma_pct = 4.0
    SgrA_residual_full = 19.7  # full deviation TGP vs EHT 2022
    print(f"  Sgr A* full deviation (TGP vs EHT obs): +{SgrA_residual_full:.1f}%")

    SgrA_env_reach_2sigma = SgrA_env_sigma_pct * 1.96
    t4_2 = SgrA_env_reach_2sigma >= SgrA_residual_full
    pass_count += int(t4_2)
    print(f"  Environment 95%-CL reach (2-sigma): +/- {SgrA_env_reach_2sigma:.1f}%")
    print(f"  T4.2 Sgr A* env can absorb deviation: {'PASS' if t4_2 else 'FAIL'}")
    if not t4_2:
        gap = SgrA_residual_full - SgrA_env_reach_2sigma
        print(f"  Gap remaining after env: +{gap:.1f}% (~{gap/SgrA_env_sigma_pct:.1f} sigma)")

    # --------------------------------------------------------------
    # T4.3 Combined verdict
    # --------------------------------------------------------------
    print(f"\n[T4.3] Combined OP-EHT environment verdict:")
    print(f"  M87* residual absorbed by env: {'YES' if t4_1 else 'NO'}")
    print(f"  Sgr A* residual absorbed by env: {'YES' if t4_2 else 'NO'}")

    # Combined PASS only if both M87* and Sgr A* can be reconciled
    t4_3 = t4_1 and t4_2
    pass_count += int(t4_3)
    print(f"  T4.3 Both sources reconciled by env: {'PASS' if t4_3 else 'FAIL'}")

    if not t4_3:
        if t4_1 and not t4_2:
            print(f"  ⇒ M87* OK, Sgr A* CANNOT be saved by environment alone.")
            print(f"  ⇒ Deferral to T5 (ngEHT 2030+) MANDATORY.")
            print(f"  ⇒ Or: Q-renorm scenario (e) from T3 must be derived.")

    # --------------------------------------------------------------
    # Summary
    # --------------------------------------------------------------
    print(f"\n{'=' * 70}")
    print(f" SUMMARY: {pass_count}/{n_total} PASS")
    print(f" OP-EHT T4 verdict: ", end="")
    if pass_count == n_total:
        print("OP-EHT T4 closes POSITIVE — env absorbs deviation.")
    elif t4_1 and not t4_2:
        print("M87* env OK, Sgr A* env INSUFFICIENT.")
        print(" ⇒ Sgr A* tension persists — T5 ngEHT needed for falsification.")
    elif not t4_1 and not t4_2:
        print("BOTH sources env INSUFFICIENT.")
        print(" ⇒ OP-EHT NEGATIVE on env path — T3 q-renorm or T5 falsification.")
    else:
        print(f"PARTIAL ({pass_count}/{n_total}) — see T5 for resolution path.")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
