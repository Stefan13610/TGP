"""
ps04 — TGP exponential-metric derivation of the lepton-Phi coupling.

Central question
----------------
ps01/ps03 used the ANSATZ Delta a_l = (alpha/2 pi) xi (m_l/M_TGP)^2,
which is a FLAVOR-UNIVERSAL coupling scaling. Does the minimal TGP
exponential metric g_{mu nu} = eta_{mu nu} * exp(2 Phi/Phi_0) actually
produce this m^2 scaling, or something else?

Method
------
Expand the QED+matter Lagrangian in a conformally-rescaled metric, field-
redefine to canonical kinetic terms, and read off:
  (1) the photon-Phi coupling (c_vert, relevant for ps03 ansatz)
  (2) the fermion-Phi coupling (relevant for direct vertex correction)

Known 1-loop formulae (Jegerlehner, Leveille):
  - Universal VECTOR exchange (massive photon-like Z'): a_l ~ (1/12 pi^2)(g/M)^2 m_l^2
  - Flavor-universal SCALAR (y = g): a_l = (y^2/8 pi^2)(m_l/M)^2 [log - 7/6 + O(m_l/M)^2]
  - Mass-weighted SCALAR (y_l = m_l/v): a_l ~ (m_l^2/v^2)/(8 pi^2) * (m_l/M)^2 [log]
    = (1/8 pi^2) * m_l^4/(v^2 M^2) * log — i.e. m_l^4 scaling!

=> TGP exponential-metric minimal coupling gives a DILATON-LIKE coupling
   y_l = m_l / Phi_0, which leads to m_l^4 scaling, NOT m_l^2.

Implications
------------
The ps01/ps03 m^2 ansatz is NOT the prediction of minimal TGP. Three
possibilities:

  (A) TGP really predicts m_l^4 scaling. Then:
      - Delta a_e / Delta a_mu = (m_e/m_mu)^4 = 5.5e-10
      - Delta a_tau / Delta a_mu = (m_tau/m_mu)^4 = 80,000
      - For Delta a_mu = 1e-9, Delta a_tau = 8e-5 (close to DELPHI bound 1.3e-2,
        well above FCC-ee future sensitivity 1e-6).
      Strongly falsifiable and DIFFERENT from ps01/ps03 predictions.

  (B) TGP has an additional non-minimal UNIVERSAL coupling (beyond minimal
      exp metric), giving flavor-universal y. Then m_l^2 scaling is recovered
      but we need a second TGP parameter.

  (C) The photon propagator dressing picture of ps03 is a DIFFERENT mechanism
      from the scalar Yukawa of minimal metric coupling. If the ps03
      epsilon_TGP(k^2) is produced by a non-minimal Phi*F^2 operator (not by
      minimal coupling), it is INDEPENDENT from the fermion dilaton Yukawa.

Script
------
Evaluates the scalar-exchange contribution to a_mu in both the mass-weighted
(dilaton) and flavor-universal scenarios, and compares to ps01/ps03 results.
"""

from __future__ import annotations

import math

import sympy as sp


# Constants
ALPHA = 1 / 137.035999
PI = math.pi
M_E_GEV = 0.0005109989
M_MU_GEV = 0.1056583745
M_TAU_GEV = 1.77686
V_HIGGS_GEV = 246.22  # SM Higgs VEV

# Target: observed muon anomaly
DELTA_A_MU_OBS = 249e-11  # WP20 data-driven (sensitivity region)
DELTA_A_MU_LATTICE = 99e-11  # BMW24 lattice


def scalar_loop_function(x: float) -> float:
    """
    Scalar-exchange contribution to a_l, Leveille 1978 / Jegerlehner:

        F_S(x) = integral_0^1 dz (1-z)^2 (1+z) / [(1-z)^2 + z x^2]

    where x = m_l/m_S. For m_S >> m_l (x<<1):
        F_S(x) ~ -7/6 - ln(x^2) + O(x^2)
    """
    import scipy.integrate as si

    def integrand(z):
        return (1 - z) ** 2 * (1 + z) / ((1 - z) ** 2 + z * x * x)

    val, _ = si.quad(integrand, 0, 1)
    return val


def delta_a_scalar(y_coupling: float, m_lepton: float, m_scalar: float) -> float:
    """
    a_l from scalar exchange (Jegerlehner eq. 4.24):
        a_l = y^2 / (8 pi^2) * (m_l / m_S)^2 * F_S(m_l/m_S)
    """
    x = m_lepton / m_scalar
    F = scalar_loop_function(x)
    return y_coupling**2 / (8 * PI * PI) * x * x * F


def scenario_A_dilaton(Phi_0_GeV: float, M_Phi_GeV: float):
    """
    Dilaton-like coupling: y_l = m_l / Phi_0, giving m_l^4 in a_l.
    Phi_0 is the substrate VEV, M_Phi is the substrate quantum mass.
    """
    out = {}
    for name, m_l in (("e", M_E_GEV), ("mu", M_MU_GEV), ("tau", M_TAU_GEV)):
        y = m_l / Phi_0_GeV
        out[name] = delta_a_scalar(y, m_l, M_Phi_GeV)
    return out


def scenario_B_universal(g_universal: float, M_Phi_GeV: float):
    """
    Flavor-universal scalar: y_l = g (same for all leptons), giving m_l^2 scaling.
    """
    out = {}
    for name, m_l in (("e", M_E_GEV), ("mu", M_MU_GEV), ("tau", M_TAU_GEV)):
        out[name] = delta_a_scalar(g_universal, m_l, M_Phi_GeV)
    return out


def fit_Phi0_dilaton(M_Phi_GeV: float, target_a_mu: float):
    """Solve a_mu_dilaton(Phi_0) = target for Phi_0."""
    from scipy.optimize import brentq

    def residual(log_Phi0):
        Phi0 = math.exp(log_Phi0)
        a_mu = scenario_A_dilaton(Phi0, M_Phi_GeV)["mu"]
        return a_mu - target_a_mu

    # a_mu ~ 1/Phi_0^2, so Phi_0 too small => a_mu too big (positive residual)
    # Phi_0 too big => a_mu too small (negative residual)
    log_Phi0 = brentq(residual, math.log(1e-6), math.log(1e8))
    return math.exp(log_Phi0)


def fit_g_universal(M_Phi_GeV: float, target_a_mu: float):
    """Solve a_mu_universal(g) = target for g."""
    from scipy.optimize import brentq

    def residual(log_g):
        g = math.exp(log_g)
        a_mu = scenario_B_universal(g, M_Phi_GeV)["mu"]
        return a_mu - target_a_mu

    log_g = brentq(residual, math.log(1e-6), math.log(100))
    return math.exp(log_g)


def main() -> None:
    print("=" * 72)
    print("ps04 — TGP exponential metric: m^2 or m^4 scaling for lepton g-2?")
    print("=" * 72)

    print("\n--- Symbolic expansion of exp metric ---")
    phi = sp.symbols("phi", real=True)
    sqrt_g = sp.exp(4 * phi)
    g_inv_squared = sp.exp(-4 * phi)
    L_maxwell = sqrt_g * g_inv_squared
    print(f"Maxwell Lagrangian prefactor: sqrt(-g) * (g^{{-1}})^2 = {sp.simplify(L_maxwell)}")
    print(f"=> Classical photon-Phi coupling vanishes (conformal invariance).")
    print(f"=> Tree-level c_vert = 0 in photon sector.")

    print(f"\nFermion mass term sqrt(-g) m psibar psi after field redefinition:")
    mass_prefactor = sp.series(sp.exp(phi), phi, 0, 3).removeO()
    print(f"  m * (1 + phi + phi^2/2 + ...) psibar psi")
    print(f"=> Phi couples to leptons with Yukawa y_l = m_l / Phi_0.")
    print(f"=> This is a DILATON-LIKE mass-weighted coupling (scenario A).")

    print("\n--- Scenario A: dilaton-like (m^4 scaling) ---")
    # Fit Phi_0 at a few M_Phi:
    for M_Phi in (10.0, 91.2, 500.0):
        Phi_0_fit = fit_Phi0_dilaton(M_Phi, DELTA_A_MU_OBS)
        a_values = scenario_A_dilaton(Phi_0_fit, M_Phi)
        ratio_tau_mu = a_values["tau"] / a_values["mu"]
        ratio_e_mu = a_values["e"] / a_values["mu"]
        print(f"  M_Phi = {M_Phi:6.1f} GeV => Phi_0 = {Phi_0_fit:8.3f} GeV")
        print(
            f"    a_e = {a_values['e']:.3e}, a_mu = {a_values['mu']:.3e}, a_tau = {a_values['tau']:.3e}"
        )
        print(
            f"    ratios: a_tau/a_mu = {ratio_tau_mu:.3e} (m^4 pred = {(M_TAU_GEV/M_MU_GEV)**4:.3e})"
        )
        print(f"            a_e/a_mu   = {ratio_e_mu:.3e}  (m^4 pred = {(M_E_GEV/M_MU_GEV)**4:.3e})")

    print("\n--- Scenario B: flavor-universal (m^2 scaling) ---")
    for M_Phi in (10.0, 91.2, 500.0):
        g_fit = fit_g_universal(M_Phi, DELTA_A_MU_OBS)
        a_values = scenario_B_universal(g_fit, M_Phi)
        ratio_tau_mu = a_values["tau"] / a_values["mu"]
        ratio_e_mu = a_values["e"] / a_values["mu"]
        print(f"  M_Phi = {M_Phi:6.1f} GeV => g_univ = {g_fit:.3e}")
        print(
            f"    a_e = {a_values['e']:.3e}, a_mu = {a_values['mu']:.3e}, a_tau = {a_values['tau']:.3e}"
        )
        print(
            f"    ratios: a_tau/a_mu = {ratio_tau_mu:.3e} (m^2 pred = {(M_TAU_GEV/M_MU_GEV)**2:.3e})"
        )
        print(f"            a_e/a_mu   = {ratio_e_mu:.3e}  (m^2 pred = {(M_E_GEV/M_MU_GEV)**2:.3e})")

    print("\n--- Dilaton-specific prediction for electron and tau ---")
    # Use Phi_0 fitted at M_Phi = 91.2 GeV:
    M_Phi = 91.2
    Phi_0_fit = fit_Phi0_dilaton(M_Phi, DELTA_A_MU_OBS)
    a_vals = scenario_A_dilaton(Phi_0_fit, M_Phi)
    print(f"  Fitting Phi_0 from Delta a_mu = {DELTA_A_MU_OBS:.2e} (WP20):")
    print(f"    Phi_0 = {Phi_0_fit:.2f} GeV")
    print(f"    Delta a_e   (dilaton) = {a_vals['e']:.3e}")
    print(f"                   [cf. experimental precision 1.3e-13, not excluded]")
    print(f"    Delta a_tau (dilaton) = {a_vals['tau']:.3e}")
    print(f"                   [cf. DELPHI 2004 bound |a_tau| < 1.3e-2]")

    # Check consistency with electron g-2 experiment
    ae_bound = 5e-13
    if abs(a_vals["e"]) > ae_bound:
        print(f"\n    WARNING: a_e = {a_vals['e']:.3e} exceeds current electron bound {ae_bound:.2e}")
        print(f"    Dilaton scenario A may be EXCLUDED by electron g-2!")
    else:
        print(f"\n    a_e = {a_vals['e']:.3e} < experimental {ae_bound:.2e}: dilaton survives e-bound.")

    print("\n--- Comparison with Phi_0 = v_Higgs (SM Higgs VEV) ---")
    Phi_0_Higgs = V_HIGGS_GEV
    # Solve for required M_Phi given Phi_0 = v:
    from scipy.optimize import brentq

    def residual(log_M):
        M = math.exp(log_M)
        return scenario_A_dilaton(Phi_0_Higgs, M)["mu"] - DELTA_A_MU_OBS

    try:
        log_M = brentq(residual, math.log(1e-3), math.log(1e5))
        M_Phi_Higgs = math.exp(log_M)
        print(f"  If Phi_0 = v_Higgs = {V_HIGGS_GEV} GeV, then M_Phi = {M_Phi_Higgs:.4f} GeV")
        print(f"  (extracting the substrate scale from the Higgs-VEV identification)")
    except Exception as e:
        print(f"  Brentq failed for Phi_0 = v_Higgs case: {e}")

    print("\n--- Key takeaway ---")
    print("If TGP substrate couples ONLY through the exp metric, then:")
    print("  * Prediction is m^4 scaling (scenario A), NOT m^2.")
    print("  * ps01/ps03 m^2 ansatz requires an ADDITIONAL non-minimal coupling.")
    print("  * Scenario A a_tau ~ 8e-5 >> a_tau^{m^2} ~ 7e-7")
    print("  * FCC-ee (sensitivity ~1e-6) can distinguish scenario A from m^2 ansatz.")
    print("")
    print("This refines the muon g-2 program to a 3-way test:")
    print("  1. Flat (BSM LFU, a_tau/a_mu = 1): excluded by existing bounds if confirmed")
    print("  2. m^2 (SM-loop-like): ps01/ps03 ansatz, needs non-minimal TGP coupling")
    print("  3. m^4 (TGP exp metric minimal): natural from exp(2 Phi/Phi_0) alone")


if __name__ == "__main__":
    main()
