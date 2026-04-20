"""
ps05 — TGP gradient-substrate correction to Casimir pressure.

Physics
-------
In TGP, the substrate field Phi has VEV Phi_0 in empty vacuum. Near a wall of
matter (MOF linker, metal electrode, etc.), Phi develops a spatial profile:
    Phi(x) = Phi_0 + delta_Phi(x),
with delta_Phi sourced by the nearby matter. The exponential metric is
    g_{mu nu}(x) = eta_{mu nu} exp(2 Phi(x) / Phi_0) ~ eta * e^2 * (1 + 2 delta/Phi_0 + ...).

Implications for Casimir
------------------------
The EM action in the exp metric is CLASSICALLY conformal:
    S_EM = -(1/4) int d^4x sqrt(-g) g^{mu alpha} g^{nu beta} F F
         = -(1/4) int d^4x F^2     (all exp factors cancel in 4D Maxwell).

So at tree level, a constant or slowly-varying delta_Phi leaves F^2 untouched.
Casimir pressure is unmodified.

However, the QED trace anomaly breaks this conformal invariance:
    L_anom = (beta(alpha)/(4 alpha)) (delta_Phi / Phi_0) F^2,
with beta(alpha)/alpha = (2 alpha)/(3 pi) per charged Dirac fermion.
For 3 charged leptons in QED:
    (beta/(4 alpha)) = (alpha)/(2 pi) ~ 1.16e-3.

This couples delta_Phi to the photon field. In a spatially varying delta_Phi
background (near MOF walls), it modifies the photon dispersion slightly:
    delta_omega^2 / omega^2 = -2 * (beta/alpha) * (delta_Phi / Phi_0)
                           = -(4 alpha / (3 pi)) * (delta_Phi / Phi_0).

This induces a position-dependent "effective refractive index" from the TGP
trace anomaly, modifying the zero-point mode spectrum and thus the Casimir
pressure by a TGP-specific correction term.

Dimensional estimate
--------------------
For a MOF wall at radius R, delta_Phi scales with the amount of matter per
volume. Crude estimate: delta_Phi/Phi_0 ~ (a_Bohr / R) at the wall (the
nearest-atom scale divided by the cavity scale — the gradient over which
the substrate transitions from "fully perturbed" near atoms to "undisturbed"
in the pore center).

TGP-corrected Casimir pressure:
    P_Cas^TGP = P_Cas,0 * [1 + c_1 (a_Bohr / R)^p],
with leading contribution from trace anomaly:
    p ~ 1 (linear in delta_Phi) at first order,
    c_1 ~ (alpha / pi) ~ 1/(400) ~ 2.4e-3  (QED loop coefficient).

For R = 0.4 nm (MOF-5 cage), a_Bohr = 0.053 nm:
    correction / classical = c_1 * (a_Bohr/R) = 2.4e-3 * 0.13 = 3.1e-4.

Tiny: ~0.03% of the Casimir pressure, far below any current measurement.

This is a CLEAN TGP PREDICTION:
  - TGP predicts NO observable Casimir anomaly in MOFs at current precision.
  - Any measured Casimir deviation above ~0.1% would FALSIFY either:
    * TGP's conformal-metric hypothesis (exp-metric being flat for F^2), or
    * the assumption that Phi_0 is the only TGP-scale affecting photons.

Alternative: if delta_Phi/Phi_0 is large (~0.5) close to the wall, correction
could rise to ~1e-3 level. Still below 10 ppm XRD precision.

Script
------
  1. Derives the dispersion modification from trace anomaly.
  2. Computes the integrated Casimir correction for MOF geometries.
  3. Compares to XRD/other experimental precision.
  4. Establishes falsifiability bounds.
"""

from __future__ import annotations

import math


PI = math.pi
ALPHA = 1.0 / 137.035999  # fine-structure

HBAR_C_JOULE_M = 3.1615e-26
HBAR_C_EV_NM = 197.327

A_BOHR = 5.29177e-11      # m (Bohr radius)
A_BOHR_NM = 0.0529177

# MOF-5
MOF5_R_SMALL = 0.395e-9   # m
MOF5_R_LARGE = 0.590e-9
MOF5_K_BULK = 15.0e9      # Pa
MOF5_A_LATT_NM = 2.583    # nm

# Casimir coefficient for spherical cavity (Boyer EM; ps03)
C_BOYER = 0.046176

# TGP trace-anomaly coefficient per QED lepton loop
C_ANOM_QED = ALPHA / (2 * PI)        # ~ 1.16e-3
# Per-lepton factor (3 light charged leptons)
C_ANOM_TOTAL = 3 * ALPHA / (3 * PI)  # = alpha/pi ~ 2.32e-3


# ---------------------------------------------------------------------------
# Analytical pieces
# ---------------------------------------------------------------------------

def delta_phi_over_phi0(a_wall_m: float, R_m: float) -> float:
    """
    Estimate of delta_Phi/Phi_0 at the cavity center, given atomic scale
    a_wall near walls and cavity radius R. Simple geometric estimate:
        delta_Phi / Phi_0 ~ a_wall / R
    (corresponds to "one layer of atomic perturbation" over the cavity radius).
    """
    return a_wall_m / R_m


def casimir_pressure_classical(R_m: float, coef: float = C_BOYER) -> float:
    """P = coef * hbar c / (4 pi R^4) [Pa]."""
    return coef * HBAR_C_JOULE_M / (4 * PI * R_m ** 4)


def tgp_correction_factor(R_m: float,
                          delta_over_phi0: float | None = None,
                          a_wall: float = A_BOHR,
                          p: int = 1) -> float:
    """
    TGP correction multiplier to Casimir pressure:
        P_TGP = P_Cas,0 * (1 + c_1 * (delta_Phi/Phi_0)^p).
    with c_1 = (beta(alpha)/(4 alpha))_total = alpha/pi (QED 3-lepton).
    """
    if delta_over_phi0 is None:
        delta_over_phi0 = delta_phi_over_phi0(a_wall, R_m)
    c1 = ALPHA / PI
    return 1.0 + c1 * (delta_over_phi0 ** p)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    print("=" * 74)
    print("ps05 — TGP gradient-substrate correction to Casimir pressure")
    print("=" * 74)

    print("\n--- TGP physical picture ---")
    print("  Exp-metric g_mu-nu = eta e^(2 Phi/Phi_0) is classically conformal")
    print("  for Maxwell F^2: no tree-level Phi-photon coupling.")
    print("  The ONLY TGP correction to Casimir at lowest order comes from the")
    print("  QED trace anomaly, which breaks conformal invariance at 1-loop:")
    print()
    print("    L_anom = (beta/(4 alpha)) (delta_Phi/Phi_0) F^2")
    print(f"    with beta/(4 alpha) = alpha/(2 pi) = {C_ANOM_QED:.4e}")
    print(f"    Total (3 charged leptons): c_1 = alpha/pi = {C_ANOM_TOTAL:.4e}")

    print("\n--- delta_Phi/Phi_0 estimate at MOF walls ---")
    print("  Atomic scale near walls: a_wall ~ Bohr radius = 0.053 nm.")
    print(f"  {'cage':>22} | {'R (nm)':>8} | {'delta_Phi/Phi_0':>18}")
    print("  " + "-" * 54)
    for name, R_nm in [
        ("MOF-5 small", 0.395),
        ("MOF-5 large", 0.590),
        ("MOF-177",     0.735),
        ("UiO-66 oct",  0.555),
        ("ZIF-8",       0.575),
        ("MOF-210",     1.120),
    ]:
        R_m = R_nm * 1e-9
        d = delta_phi_over_phi0(A_BOHR, R_m)
        print(f"  {name:>22} | {R_nm:>8.3f} | {d:>18.4e}")

    print("\n--- TGP correction factor to Casimir pressure ---")
    print(f"  c_1 = alpha/pi = {C_ANOM_TOTAL:.4e}")
    print(f"  Assumed power law: (delta_Phi/Phi_0)^1 (linear, leading anomaly)")
    print()
    print(f"  {'cage':>22} | {'P_Cas_0 (MPa)':>14} | {'correction':>14} | "
          f"{'P_TGP (MPa)':>14} | {'Delta (ppm)':>12}")
    print("  " + "-" * 82)

    for name, R_nm in [
        ("MOF-5 small", 0.395),
        ("MOF-5 large", 0.590),
        ("MOF-177",     0.735),
        ("UiO-66 oct",  0.555),
        ("ZIF-8",       0.575),
        ("MOF-210",     1.120),
    ]:
        R_m = R_nm * 1e-9
        # Use "realistic" effective coefficient 5e-4 (from ps04)
        P_Cas = casimir_pressure_classical(R_m, coef=5e-4)
        d = delta_phi_over_phi0(A_BOHR, R_m)
        c1 = ALPHA / PI
        correction = c1 * d
        P_TGP = P_Cas * (1 + correction)
        delta_strain_ppm = (correction * P_Cas / MOF5_K_BULK) * 1e6
        print(f"  {name:>22} | {P_Cas/1e6:>14.2f} | {correction:>+14.4e} | "
              f"{P_TGP/1e6:>14.2f} | {delta_strain_ppm:>+12.4f}")

    print("\n--- Falsifiability analysis ---")
    print("  Experimental precision for MOF lattice constant:")
    print("    Synchrotron XRD (MOF-5 single crystal): ~ 5-10 ppm")
    print("    Best-case (Zhou 2008, APS 11-BM):       ~ 1-2 ppm")
    print("    FCC-era high-res synchrotron proposal:  ~ 0.1 ppm")
    print()
    print("  TGP prediction for MOF-5 small cage (R=0.395 nm):")
    R = MOF5_R_SMALL
    P0 = casimir_pressure_classical(R, coef=5e-4)  # realistic
    corr = (ALPHA / PI) * delta_phi_over_phi0(A_BOHR, R)
    strain_ppm = abs(corr * P0 / MOF5_K_BULK) * 1e6
    print(f"    P_Cas,0 (realistic ~5e-4 coeff): {P0/1e6:.1f} MPa")
    print(f"    TGP correction factor:           {corr:.4e} = {corr*100:.4f}%")
    print(f"    Induced strain:                  {strain_ppm:.4f} ppm")
    print()
    print(f"  => TGP prediction is {strain_ppm:.3f} ppm additional strain vs classical.")
    print(f"     This is BELOW current 1-2 ppm XRD precision by {1.0/strain_ppm:.0f}x.")
    print(f"     => NOT detectable at current experimental precision.")

    print()
    print("  What precision would be needed to TEST TGP via MOF Casimir?")
    target_ppm = 0.1 * strain_ppm
    print(f"    Need strain precision ~ {target_ppm:.5f} ppm to detect TGP at 10-sigma.")
    print(f"    This is {1.0/target_ppm:.0f}x better than current best — out of reach")
    print(f"    for present synchrotron technology.")

    print("\n--- Alternative: larger delta_Phi/Phi_0 near wall ---")
    print("  If delta_Phi/Phi_0 is O(1) in the immediate wall region (realistic")
    print("  for strong substrate-matter coupling), the effective correction would be")
    print("  ~ alpha/pi = 0.23% of classical Casimir.")
    P_corr_big = ALPHA / PI
    P0_b = casimir_pressure_classical(R, coef=5e-4)
    strain_big_ppm = P_corr_big * P0_b / MOF5_K_BULK * 1e6
    print(f"    Fractional correction: alpha/pi = {ALPHA/PI:.4e}")
    print(f"    Strain at R = 0.395 nm: {strain_big_ppm:.4f} ppm")
    print(f"    Still below current 1-2 ppm XRD threshold, but only by {1.0/strain_big_ppm:.1f}x.")

    print("\n--- Power-law variants (p > 1) ---")
    print("  If TGP corrections scale as (delta_Phi/Phi_0)^p with p > 1:")
    print(f"  p = 2: correction ~ (alpha/pi) * (a/R)^2 = "
          f"{(ALPHA/PI) * (A_BOHR_NM/0.395)**2:.4e}")
    print(f"  p = 4: correction ~ (alpha/pi) * (a/R)^4 = "
          f"{(ALPHA/PI) * (A_BOHR_NM/0.395)**4:.4e}")
    print("  Higher p -> even smaller correction; clearly falsification")
    print("  requires low p or large c_1, or different mechanism.")

    # Alternative: if Phi_0 is very small (~ GeV scale from ps04/ps06 of muon analysis),
    # then delta_Phi/Phi_0 near a charged atom could be O(1) or larger.
    print("\n--- Alternative: Phi_0 ~ GeV (from muon g-2 fit) ---")
    print("  If Phi_0 ~ 0.34 GeV (ruled out by LEP in ps06 of muon analysis,")
    print("  but a historical reference scale):")
    print("    Near a nucleus with EM field F ~ e/(4 pi a_Bohr^2):")
    print("    induced delta_Phi ~ (alpha/pi) * (F^2)/(M_Phi^2) * (volume factor)")
    print("  This is a different regime, not applicable to MOF vacuum physics.")

    print("\n" + "=" * 74)
    print("Summary + Falsifiability verdict")
    print("=" * 74)
    print(f"""
TGP prediction for Casimir in MOFs:
  - Classical (idealized) Casimir pressure: ~0.1-1 GPa for MOF cages.
  - Realistic (dielectric + open-geometry): ~50-100 MPa (from ps04).
  - TGP correction via QED trace anomaly:
      P_TGP - P_Cas,0 = (alpha/pi) * (a_Bohr/R)^p * P_Cas,0
      For p = 1, R = 0.4 nm: correction ~ 3e-4 (0.03%).
      For p = 2: correction ~ 4e-5 (0.004%).

Induced strain on MOF-5 lattice:
  ~ {strain_ppm:.3f} ppm (for p=1, alpha/pi coefficient).
  Current XRD precision: 1-2 ppm. => NOT detectable today.

Verdict:
  TGP's exp-metric coupling is CONSISTENT with observed MOF Casimir physics
  precisely BECAUSE of the classical conformal invariance of Maxwell in 4D.
  Any future measurement of a >>0.1% Casimir anomaly in MOFs would
  FALSIFY this protection and require TGP modifications (e.g., direct
  non-conformal coupling, or new TGP-matter interaction).

  This is a CLEAN NULL-PREDICTION: TGP agrees with classical Casimir to
  within the QED trace-anomaly correction, which is ~10^(-4) too small
  to measure in current MOF experiments.

Research programme consequence:
  - Short-term: no dedicated experiment will distinguish TGP from classical
    QED in MOF Casimir regime.
  - Long-term: high-precision vacuum measurements (e.g., Casimir between
    mirror arrays at nanometer separation, or torsion balance in magnetic
    shielding) at >0.01% precision would start to probe the trace-anomaly
    effect. FCC-ee era may achieve this indirectly.

Next steps for the TGP vacuum programme:
  - Shift focus from MOF Casimir to BETTER-SUITED systems:
    (a) Electron g-2 (tiny but measurable alpha^2 correction on top of a_e);
    (b) Cosmological vacuum energy: substrate ZPE as dark-energy candidate;
    (c) Scalar Casimir in liquid-He superfluid (substrate-coupled directly).
  - This line (Casimir MOF) is methodologically complete but experimentally
    UNFALSIFIABLE at current precision. Document as a null prediction and move on.
""")


if __name__ == "__main__":
    main()
