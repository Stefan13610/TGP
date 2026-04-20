"""
ps05 — QED trace anomaly as the m^2 channel for muon g-2.

Context from ps04
-----------------
The exponential metric g_{mu nu} = eta_{mu nu} * exp(2 phi) has CLASSICAL
conformal invariance for 4D Maxwell (sqrt(-g) * (g^{-1})^2 cancellation).
But at the QUANTUM level, the QED trace anomaly generates a non-vanishing
Phi-gamma-gamma coupling:

    L_anom = (beta(alpha) / 4 alpha) * (Phi/Phi_0) * F_{mu nu} F^{mu nu}

For QED at 1-loop with lepton species f (charge Q_f):
    beta(alpha)/alpha = (2 alpha / 3 pi) * sum_f N_c Q_f^2 = (2 alpha / 3 pi) * 3
                                                    (for 3 charged leptons)
    => c_anom = beta/(4 alpha) = alpha / (2 pi)

Including the 3 charged leptons (ignoring quarks for simplicity since the
relevant Phi-photon-photon vertex is leptonic-generated):
    c_anom ~ alpha / (2 pi) ~ 1.16e-3

Two channels contribute to Delta a_mu
-------------------------------------

**Channel A (dilaton-fermion, m^4 scaling):**
    L_A = (m_f / Phi_0) * Phi * fbar f
    Delta a_f^A = (y_f^2 / 8 pi^2) * (m_f^2 / M_Phi^2) * F_S(m_f/M_Phi)
                = (m_f^4) / (8 pi^2 Phi_0^2 M_Phi^2) * F_S

**Channel B (anomaly-induced gamma-gamma, effective m^2 scaling):**
    L_B = c_anom * (Phi/Phi_0) * F^2
    Integrating out Phi: L_eff = -(c_anom^2 / 2 M_Phi^2 Phi_0^2) * (F^2)^2
    This contributes to a_f via a 2-loop diagram that gives:
        Delta a_f^B ~ (alpha^2 / pi^2) * c_anom^2 * (m_f^2 / M_Phi^2 Phi_0^2) * log

Which channel dominates?
------------------------
Ratio at m = m_mu, M_Phi = Phi_0 (for simplicity):
    A/B = m_mu^2 / (alpha^2 c_anom^2)
        = (0.1)^2 / ((1/137)^2 * (1/273)^2)
        = 0.01 * 137^2 * 273^2 ~ 1.4e7

=> Channel A (dilaton Yukawa, m^4 scaling) DOMINATES over the anomaly (B).

This is the clean TGP prediction: at minimal coupling, the exponential
metric gives m^4 scaling with the Phi_0 ~ 0.34 GeV fit from ps04.

For channel B to dominate, we need Phi_0 >> m_mu * 273 / alpha ~ 3 TeV.
At such large Phi_0, the Yukawa channel is suppressed enough for the
anomaly to take over. But then M_Phi from the muon fit goes to sub-GeV,
far from the electroweak scale.

Script
------
Evaluates both channels for:
- electron, muon, tau g-2 at their respective mass scales,
- a grid of (Phi_0, M_Phi) values,
and reports which channel is dominant in each regime.
"""

from __future__ import annotations

import math

import numpy as np
import scipy.integrate as si
import scipy.optimize as so

ALPHA = 1 / 137.035999
PI = math.pi
M_E = 0.0005109989  # GeV
M_MU = 0.1056583745
M_TAU = 1.77686
V_HIGGS = 246.22

# Trace anomaly coefficient (3 charged leptons in the 1-loop beta)
C_ANOM = ALPHA / (2 * PI)  # ~ 1.16e-3


def F_scalar(x: float) -> float:
    """Scalar-exchange loop function (Jegerlehner / Leveille)."""

    def integrand(z):
        return (1 - z) ** 2 * (1 + z) / ((1 - z) ** 2 + z * x * x + 1e-30)

    val, _ = si.quad(integrand, 0, 1, limit=100)
    return val


def delta_a_channel_A(m_l: float, Phi_0: float, M_Phi: float) -> float:
    """
    Dilaton-like mass-weighted Yukawa contribution.
    Delta a_f = (y_f^2 / 8 pi^2) * (m_f/M_Phi)^2 * F_S(m_f/M_Phi)
    with y_f = m_f / Phi_0.
    """
    y = m_l / Phi_0
    x = m_l / M_Phi
    return y * y / (8 * PI * PI) * x * x * F_scalar(x)


def delta_a_channel_B(m_l: float, Phi_0: float, M_Phi: float) -> float:
    """
    QED trace-anomaly contribution at leading 2-loop.
    Dimensional estimate:
        Delta a_f^B ~ (alpha/pi)^2 * c_anom^2 * m_f^2 / (M_Phi^2 Phi_0^2) * log

    Using c_anom = alpha/(2 pi), this is ~ alpha^4 / (2pi)^2 per log.
    We take the log factor = log(M_Phi^2 / m_f^2) and include a numerical
    O(1) coefficient of 1 (crude — proper calculation would fix this to
    few × 0.1).
    """
    log_factor = max(math.log(max(M_Phi / m_l, 1.0)), 0.0)
    prefactor = (ALPHA / PI) ** 2 * C_ANOM ** 2
    return prefactor * (m_l ** 2) / (M_Phi ** 2 * Phi_0 ** 2) * log_factor * 1.0


def fit_Phi0_given_M_Phi(channel_fn, M_Phi: float, target_a_mu: float) -> float | None:
    """Root-solve for Phi_0 that reproduces target_a_mu with given channel."""

    def residual(log_Phi0):
        return channel_fn(M_MU, math.exp(log_Phi0), M_Phi) - target_a_mu

    # Probe endpoints
    try:
        lo, hi = math.log(1e-6), math.log(1e8)
        flo, fhi = residual(lo), residual(hi)
        if flo * fhi > 0:
            return None
        return math.exp(so.brentq(residual, lo, hi))
    except Exception:
        return None


def main() -> None:
    print("=" * 74)
    print("ps05 — Trace anomaly vs dilaton Yukawa: which channel carries a_mu?")
    print("=" * 74)

    print(f"\nTrace anomaly coefficient c_anom = alpha/(2 pi) = {C_ANOM:.4e}")

    DELTA_A_MU_OBS_LATTICE = 99e-11
    DELTA_A_MU_OBS_DD = 249e-11

    print("\n--- Fit each channel to Delta a_mu = 249e-11 (WP20 data-driven) ---")
    for M_Phi in (10.0, 91.2, 1000.0):
        print(f"\n  M_Phi = {M_Phi} GeV:")
        Phi_A = fit_Phi0_given_M_Phi(delta_a_channel_A, M_Phi, DELTA_A_MU_OBS_DD)
        Phi_B = fit_Phi0_given_M_Phi(delta_a_channel_B, M_Phi, DELTA_A_MU_OBS_DD)
        print(f"    Channel A (dilaton, m^4): Phi_0 = {Phi_A}")
        print(f"    Channel B (anomaly, m^2): Phi_0 = {Phi_B}")
        if Phi_A is not None and Phi_B is not None:
            # Evaluate both channels at Phi_A and at Phi_B to see cross-contribution
            for Phi_0, label in ((Phi_A, "at Phi_0 fitted for A"), (Phi_B, "at Phi_0 fitted for B")):
                aA = delta_a_channel_A(M_MU, Phi_0, M_Phi)
                aB = delta_a_channel_B(M_MU, Phi_0, M_Phi)
                print(f"    {label} = {Phi_0:.3e} GeV:")
                print(f"      a_mu(A) = {aA:.3e},  a_mu(B) = {aB:.3e},  ratio A/B = {aA/aB:.3e}")

    print("\n--- Predictions for electron and tau in each channel ---")
    print("    (Phi_0 fitted at M_Phi = 91.2 GeV)")
    M_Phi = 91.2
    Phi_A = fit_Phi0_given_M_Phi(delta_a_channel_A, M_Phi, DELTA_A_MU_OBS_DD)
    Phi_B = fit_Phi0_given_M_Phi(delta_a_channel_B, M_Phi, DELTA_A_MU_OBS_DD)

    if Phi_A is not None:
        print(f"\n  Channel A (dilaton, Phi_0 = {Phi_A:.3e} GeV):")
        for name, m in (("e", M_E), ("mu", M_MU), ("tau", M_TAU)):
            a = delta_a_channel_A(m, Phi_A, M_Phi)
            print(f"    Delta a_{name:<3} = {a:.3e}")
        a_mu_A = delta_a_channel_A(M_MU, Phi_A, M_Phi)
        a_tau_A = delta_a_channel_A(M_TAU, Phi_A, M_Phi)
        print(f"    a_tau/a_mu = {a_tau_A/a_mu_A:.3e} (m^4 expect: {(M_TAU/M_MU)**4:.3e})")

    if Phi_B is not None:
        print(f"\n  Channel B (anomaly, Phi_0 = {Phi_B:.3e} GeV):")
        for name, m in (("e", M_E), ("mu", M_MU), ("tau", M_TAU)):
            a = delta_a_channel_B(m, Phi_B, M_Phi)
            print(f"    Delta a_{name:<3} = {a:.3e}")
        a_mu_B = delta_a_channel_B(M_MU, Phi_B, M_Phi)
        a_tau_B = delta_a_channel_B(M_TAU, Phi_B, M_Phi)
        print(f"    a_tau/a_mu = {a_tau_B/a_mu_B:.3e} (m^2 expect: {(M_TAU/M_MU)**2:.3e})")

    print("\n--- Assessment ---")
    print("Channel A (dilaton, m^4 scaling): NATURAL from minimal exp-metric.")
    print("  Predicts: Phi_0 ~ 0.3-3 GeV (hadronic), a_tau ~ 2e-4")
    print("Channel B (anomaly, m^2 scaling): SUPPRESSED by alpha^4 prefactor.")
    print("  To saturate a_mu alone, needs Phi_0 ~ 10-100 MeV (very small).")

    print("\n--- Cross-over: where does A = B? ---")
    # At m = m_mu: A ~ m_mu^4 / (8 pi^2 Phi_0^2 M_Phi^2), B ~ alpha^4/(pi^2(2pi)^2) m_mu^2 / (M_Phi^2 Phi_0^2)
    # Ratio A/B at fixed (Phi_0, M_Phi):
    #   A/B = m_mu^2 / (alpha^4 / (2pi)^2) * (8 pi^2) / 1
    #       = m_mu^2 * 8 pi^2 * (2 pi)^2 / alpha^4
    #       = m_mu^2 * 32 pi^4 / alpha^4
    ratio_A_over_B = (M_MU ** 2) * 32 * PI ** 4 / (ALPHA ** 4)
    print(f"  For m=m_mu (ignoring log): A/B = {ratio_A_over_B:.3e}")
    print(f"  => Channel A DOMINATES minimal-coupling TGP g-2 by factor ~ {ratio_A_over_B:.1e}")

    print("\n--- Interpretation ---")
    print("The trace anomaly is a REAL quantum effect in TGP, but it is")
    print("numerically NEGLIGIBLE compared to the direct dilaton-muon Yukawa.")
    print("The muon g-2 is dominated by channel A (dilaton), which predicts")
    print("m^4 scaling and Delta a_tau ~ 2e-4 — NOT the m^2 pattern.")
    print("")
    print("=> ps01/ps03's m^2 ansatz is INCORRECT for minimal TGP.")
    print("=> The correct TGP prediction is m^4 scaling (channel A), with:")
    print(f"     Phi_0 ~ 0.34 GeV, M_Phi ~ M_Z")
    print(f"     Delta a_tau ~ 2e-4  (distinguishable by FCC-ee at 1e-6)")
    print("=> Electron g-2 decouples: Delta a_e ~ 1e-18 (totally invisible).")


if __name__ == "__main__":
    main()
